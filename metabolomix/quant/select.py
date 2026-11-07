"""付加体の選別と除外（Task 13。定量の最後の工程）。

人のキュレーションを当ててから、規則プラグインの `select`（内部標準 → 名前・コメント・MS/MS の
除外 → 付加体の採用表）で各特徴の残す・除外を決め、残る行だけの派生データセットを作る。

- 内部標準の行は、正規化の記録が持つ内部標準の特徴 ID（`quant_state["internal_standard_ids"]`）で
  決める。正規化していなければ role は "unassigned"（内部標準を除かない）。
- 統合後のデータセットは、極性ごとのキュレーション（`quant_state["curation_by_polarity"]`）を
  `feature_metadata[fid]["origin"]` の極性・元の特徴 ID で当て、mode も origin の極性。統合前の
  データセットは `ds.curation_flags` と、そのデータセット自身の極性。
- `assign` は name / ontology を差し替えてから規則へ。`wrong` / `redundant` は、規則が残すと
  した行だけを理由 "curation_excluded" の除外にする（msdial-workflow の `_curated` と同じ。
  内部標準など規則がもう除外した行の理由は書き換えない）。
- 全判断（と、キュレーションの audit）は記録の横の `<record_id>.audit.json` に書く。戻り値は
  理由ごとの件数と例 5 件だけ。
"""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import missing_state
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.quant import derive as qderive
from metabolomix.quant.inputs import ReferenceInput, adduct_key, adduct_table, quant_error
from metabolomix.quant.normalize import analysis_folder
from metabolomix.quant.plugin import call, load_rules, rules_settings
from metabolomix.quant.records import AUDIT_SUFFIX, write_audit, write_record
from metabolomix.quant.table import apply_curation as curate
from metabolomix.quant.table import feature_rows, resolve_quant_polarity

__all__ = ["AUDIT_SCHEMA", "is_merged", "select", "select_tool"]

AUDIT_SCHEMA = "lipidmix.quant-select-audit.v1"
_MAX_EXAMPLES = 5
_STATE_KEYS = {"schema", "normalized", "method", "unit", "record_id", "steps", "skipped"}
_HOW_ADDUCTS = ("lipidmix.local.toml の [quant] adducts に付加体の採用表（CSV。列: Class, Adduct, "
                "Ion mode, IsSelected）を設定する。付加体の選別を飛ばすなら "
                "skip_adduct_filter=True を渡す（内部標準・名前・MS/MS の除外は行う）")


def is_merged(ds) -> bool:
    """特徴がすべて統合の出自（`origin`）を持つか。"""
    metadata = getattr(ds, "feature_metadata", None) or {}
    ids = list(getattr(ds, "feature_ids", None) or [])
    return bool(ids) and all(isinstance(metadata.get(f, {}).get("origin"), dict) for f in ids)


def _curation_groups(ds, rows, polarity, apply_curation: bool):
    """(行ごとの mode, 差し替え後の行, 除外する特徴 ID の集合, audit, 極性ごとの適用の報告)。

    報告は統合後で apply_curation=True のときだけ（それ以外は None）。判断の無い極性は
    {"applied": False, "reason": "no_flags"}。"""
    state = getattr(ds, "quant_state", None) or {}
    modes: dict[str, str] = {}
    groups: dict[str, list] = {}
    origin_ids: dict[str, str] = {}
    merged = is_merged(ds)
    for row in rows:
        if merged:
            origin = ds.feature_metadata[row.feature_id]["origin"]
            key = str(origin["polarity"]).strip().lower()
            origin_ids[row.feature_id] = str(origin["feature_id"])
        else:
            key = polarity
            origin_ids[row.feature_id] = row.feature_id
        modes[row.feature_id] = key
        groups.setdefault(key, []).append(row)
    replaced: dict[str, object] = {}
    excluded: set[str] = set()
    audit: list[dict] = []
    if not apply_curation:
        return modes, {r.feature_id: r for r in rows}, excluded, audit, None
    report: dict | None = {} if merged else None
    for key, members in groups.items():
        flags = (state.get("curation_by_polarity") or {}).get(key) if merged else ds.curation_flags
        if not flags:
            if report is not None:
                report[key] = {"applied": False, "reason": "no_flags"}
            continue
        if report is not None:
            report[key] = {"applied": True}
        back = {origin_ids[r.feature_id]: r.feature_id for r in members}
        out, sub_audit = curate([replace(r, feature_id=origin_ids[r.feature_id]) for r in members],
                                flags)
        kept = {back[r.feature_id]: r for r in out}
        for row in members:
            if row.feature_id in kept:
                replaced[row.feature_id] = replace(kept[row.feature_id],
                                                   feature_id=row.feature_id)
            else:
                excluded.add(row.feature_id)
        audit.extend({**a, "feature_id": back[a["feature_id"]]} for a in sub_audit)
    final = {r.feature_id: replaced.get(r.feature_id, r) for r in rows}
    return modes, final, excluded, audit, report


def _decide(plugin, rows, modes, adduct, ds):
    state = getattr(ds, "quant_state", None) or {}
    normalized = bool(state.get("normalized"))
    is_ids = {str(i) for i in state.get("internal_standard_ids") or []} if normalized else set()
    default_role = "target" if normalized else "unassigned"
    feats = [{"feature_id": r.feature_id, "name": r.name, "comment": r.comment,
              "msms_assigned": r.msms_assigned, "ontology": r.ontology, "adduct": r.adduct,
              "mode": modes[r.feature_id],
              "role": "internal_standard" if r.feature_id in is_ids else default_role}
             for r in rows]
    result = call(plugin, "select", features=feats,
                  adduct_text=adduct.text if adduct is not None else None)
    by_id = {str(d.get("feature_id")): d for d in result or []}
    if len(by_id) != len(feats) or set(by_id) != {f["feature_id"] for f in feats}:
        raise DomainError("QUANT_RULES_BAD_RESULT",
                          "規則プラグインの select の戻りが、渡した特徴と 1 対 1 になっていません。")
    return [dict(by_id[f["feature_id"]]) for f in feats]


def select(ds, *, adduct: ReferenceInput | None, apply_curation: bool = True,
           polarity: str | None = None, reuse: dict | None = None):
    """選別して記録と audit を書き、残る行だけの派生データセットを返す。`(new, info)`。

    `adduct` が None なら付加体は見ない（規則の filter_disabled。内部標準・名前・MS/MS の除外は
    行う）。`polarity` は統合前のデータセットの mode（省くとデータセットから決め、決まらなければ
    QUANT_POLARITY_REQUIRED）。統合後は各特徴の出自の極性を使い、引数は見ない。
    `reuse` は既存の記録（復元用）: 記録と audit を書かず、その `record_id` で派生を作る。
    """
    plugin = load_rules()
    merged = is_merged(ds)
    if not merged:
        polarity = resolve_quant_polarity(ds, None, polarity)
    rows = feature_rows(ds)
    state = getattr(ds, "quant_state", None)
    folder = Path(state["record_folder"]) if state and state.get("record_folder") \
        else analysis_folder(ds)

    modes, final, excluded, curation_audit, by_polarity = _curation_groups(ds, rows, polarity, apply_curation)
    ordered = [final[r.feature_id] for r in rows]
    decisions = _decide(plugin, ordered, modes, adduct, ds)
    for d in decisions:
        if d["action"] == "retained" and d["feature_id"] in excluded:
            d["action"], d["reason"] = "excluded", "curation_excluded"
    keep = [i for i, d in enumerate(decisions) if d["action"] == "retained"]
    reasons = Counter(d["reason"] for d in decisions)
    curation = dict(Counter(a["action"] for a in curation_audit)) if apply_curation else None

    adduct_info = {"sha256": adduct.sha256, "label": adduct.label} if adduct is not None else None
    inputs = dict(ds.source_files)
    if adduct is not None:
        inputs[adduct_key(adduct.label)] = adduct.sha256
    parents = [state["record_id"]] if state and state.get("record_id") else []
    counts = {"input": len(rows), "retained": len(keep), "excluded": len(rows) - len(keep)}
    if reuse is not None:
        record_id = reuse["record_id"]
        audit_name = f"{record_id}{AUDIT_SUFFIX}"
    else:
        record = write_record(
            folder, "select", inputs=inputs, parent_ids=parents,
            settings={"adducts": adduct_info, "apply_curation": bool(apply_curation),
                      "polarity": None if merged else polarity, "rules": rules_settings(plugin)},
            summary={"n_features": counts, "reasons": dict(reasons), "curation": curation})
        record_id = record["record_id"]
        audit_name = write_audit(folder, record_id, {
            "schema": AUDIT_SCHEMA, "record_id": record_id, "decisions": decisions,
            "curation": curation_audit}).name

    kept_ids = {rows[i].feature_id for i in keep}
    extras = {k: v for k, v in (state or {}).items() if k not in _STATE_KEYS}
    if "internal_standard_ids" in extras:
        extras["internal_standard_ids"] = [i for i in extras["internal_standard_ids"]
                                           if str(i) in kept_ids]
    skipped = list((state or {}).get("skipped") or [])
    if adduct is None:
        skipped.append("adduct_filter")
    new_state = qderive.quant_state(
        normalized=bool((state or {}).get("normalized")), method=(state or {}).get("method"),
        unit=(state or {}).get("unit"), record_id=record_id,
        steps=list((state or {}).get("steps") or []) + ["select"], skipped=skipped, **extras)
    new = qderive.derive(
        ds, matrix=np.asarray(ds.feature_matrix, dtype=float)[keep],
        detected=np.asarray(ds.detected_mask, dtype=bool)[keep],
        feature_ids=[rows[i].feature_id for i in keep], state=new_state)
    for i in keep:  # キュレーションの assign は残す行の同定に反映する
        fid = rows[i].feature_id
        if final[fid] is not rows[i]:
            new.feature_metadata[fid]["name"] = final[fid].name
            new.feature_metadata[fid]["ontology"] = final[fid].ontology
    examples = [{"feature_id": d["feature_id"], "name": final[d["feature_id"]].name,
                 "reason": d["reason"]}
                for d in decisions if d["action"] == "excluded"][:_MAX_EXAMPLES]
    info = {"record_id": record_id, "audit_file": audit_name,
            "n_retained": counts["retained"], "n_excluded": counts["excluded"],
            "reasons": dict(reasons), "curation": curation, "curation_by_polarity": by_polarity,
            "adducts": adduct_info,
            "skipped": skipped, "examples": examples}
    return new, info


def select_tool(holder, *, parent_handle_id: str | None, skip_adduct_filter: bool = False,
                apply_curation: bool = True, polarity: str | None = None) -> str:
    """dataset_select_features の実体（採用表・極性の確認・登録・封筒）。"""
    ds = getattr(holder, "dataset", None)
    if ds is None:
        return missing_state(
            "dataset", ["dataset_load", "dataset_handle_open"],
            "データセットがありません。dataset_load を実行するか、dataset_handle_open で開いた"
            "ハンドルを渡してください。")
    try:
        adduct = None if skip_adduct_filter else adduct_table()
        if adduct is None and not skip_adduct_filter:
            return json_payload({
                "status": "needs_input", "reason": "adduct_table_required",
                "message": "付加体の採用表が設定されていません。", "how_to_answer": _HOW_ADDUCTS})
        if not is_merged(ds):
            polarity = resolve_quant_polarity(ds, holder, polarity)
        new, info = select(ds, adduct=adduct, apply_curation=apply_curation, polarity=polarity)
        handle_id = qderive.register(holder, new, step="select",
                                     parent_handle_ids=[parent_handle_id] if parent_handle_id else [])
    except DomainError as exc:
        return quant_error(exc.code, exc.message, exc.details or None)
    return json_payload(round_floats({"status": "selected", "dataset_handle": handle_id, **info}))
