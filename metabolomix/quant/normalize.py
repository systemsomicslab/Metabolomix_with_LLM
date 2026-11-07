"""内部標準の確認（preview）と正規化（apply / skip）。Task 11。

規則はプラグイン（`plugin.call`）の `normalize` に任せ、ここは `DatasetState` ↔ 規則の入出力の
変換・記録・派生データセットの組み立てだけをする。値は gap-fill も計算し（利用者の判断 2）、
`detected_mask` で「内部標準も対象も検出された」セルだけを True にする。

記録の置き場所は、データセット自身の解析フォルダ（`reports/quant/`）。`mcp_core._resolve_report_dir`
や DATA_DIR ではない（B-3 の規約: 統合前の記録は各極性の解析フォルダの `reports/`）。

規則の `feature_id` は int なので `int()` で渡し、戻すときに str にする。統合後の ID
（`positive:...`）は正規化の対象にしない（正規化は統合の前）。
"""
from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path

import numpy as np

from metabolomix.core.atomic_io import DomainError, atomic_write_bytes
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.quant import derive as qderive
from metabolomix.quant.inputs import ReferenceInput, quant_error, reference_key
from metabolomix.quant.plugin import call, load_rules, rules_settings
from metabolomix.quant.records import record_dir, write_record
from metabolomix.quant.table import feature_rows, require_msdial_source

__all__ = ["analysis_folder", "apply", "apply_tool", "preview", "preview_tool", "skip"]

_UNRESOLVED = ("standard_ambiguous", "standard_missing")
#: 値を残す status（それ以外のセルは NaN）。
_KEEP = ("ok", "gap_filled")
_ANY_OTHERS = "any others"  # 参照表の「それ以外のクラス」の行（assignments の鍵）


_PRIMARY_SUFFIXES = ((".arf2",), (".arf",))


def analysis_folder(ds) -> Path:
    """記録を置く解析フォルダ＝データセット自身の解析フォルダ（`.arf2`、無ければ `.arf`）の
    置き場所。サイドカー（`.mddata` `_tags.xml` `.dcl` `.aef`）は見ない。
    主ファイルが無い・別のフォルダに散るときは QUANT_RECORD_LOCATION_REQUIRED。
    """
    paths = [Path(p) for p in (getattr(ds, "source_files", None) or {})]
    for suffixes in _PRIMARY_SUFFIXES:
        found = {p.parent for p in paths if p.suffix.lower() in suffixes}
        if len(found) == 1:
            return found.pop()
        if found:
            break
    raise DomainError("QUANT_RECORD_LOCATION_REQUIRED",
                      "データセットの解析フォルダを 1 つに決められないため、記録を置けません。"
                      "MS-DIAL の解析フォルダ（.arf2 / .arf）から読んだデータセットで"
                      "実行してください。")


def _check_not_normalized(ds) -> None:
    """正規化を掛けてよいか: 解析フォルダ由来・未正規化・未選別。"""
    require_msdial_source(ds)
    state = getattr(ds, "quant_state", None)
    if state and "select" in (state.get("steps") or []):
        raise DomainError("QUANT_ALREADY_SELECTED",
                          "このデータセットは選別済みです（行を落としたあとでは内部標準の行が無い"
                          "ことがある）。正規化は選別の前に、元のハンドルから行ってください。",
                          {"record_id": state.get("record_id")})
    if state and state.get("normalized"):
        raise DomainError("QUANT_ALREADY_NORMALIZED",
                          "このデータセットは正規化済みです（二重に割らない）。元のハンドルから"
                          "やり直してください。", {"record_id": state.get("record_id")})


def _bindings(bindings) -> dict[str, int]:
    out = {}
    for key, value in (bindings or {}).items():
        try:
            out[str(key)] = int(value)
        except (TypeError, ValueError):
            raise DomainError("QUANT_BAD_REQUEST",
                              "bindings は {内部標準の名前: 特徴 ID（整数）} です（規則の対応指定は内部標準の名前が鍵）。",
                              {"field": "bindings", "standard_name": str(key)}) from None
    return out


def _features(ds, rows) -> list[dict]:
    feats = []
    for row in rows:
        try:
            fid = int(row.feature_id)
        except ValueError:
            raise DomainError(
                "QUANT_FEATURE_ID_NOT_INTEGER",
                f"特徴 ID {row.feature_id!r} は整数ではありません。正規化は統合の前に行います"
                "（統合後の ID は正規化できません）。",
                {"feature_id": row.feature_id}) from None
        feats.append({"feature_id": fid, "name": row.name, "lipid_class": row.ontology,
                      "mz": row.mz, "rt": row.rt, "adduct": row.adduct,
                      "heights": list(row.heights),
                      "detected": [True] * len(row.heights)})  # gap-fill も計算する
    return feats


def _reference(plugin, ref: ReferenceInput) -> tuple[str, str, str]:
    """(規則に渡す表の本文, その sha256, 規則に渡す method)。濃度の形なら係数表を作る。"""
    if ref.concentrations is None:
        return ref.text, ref.sha256, ref.method
    text = call(plugin, "coefficient_reference_text", reference_text=ref.text,
                concentrations=ref.concentrations, solvent_volume=ref.solvent_volume,
                cells_millions=ref.cells_millions, unit=ref.unit)
    return text, hashlib.sha256(text.encode("utf-8")).hexdigest(), "cell_coefficient"


def _compute(ds, ref: ReferenceInput, bindings):
    plugin = load_rules()
    _check_not_normalized(ds)
    rows = feature_rows(ds)
    feats = _features(ds, rows)
    text, sha, method = _reference(plugin, ref)
    result = call(plugin, "normalize", features=feats, n_samples=len(ds.sample_names),
                  reference_text=text, method=method, bindings=_bindings(bindings))
    return plugin, rows, text, sha, result


def _by_id(rows):
    return {r.feature_id: r for r in rows}


def _preview_dict(rows, ref: ReferenceInput, sha: str, polarity: str, result: dict) -> dict:
    by_id = _by_id(rows)
    any_others = sorted(result.get("any_others_classes") or ())
    classes, grouped = [], {}
    present = {c for c in result.get("target_class_by_feature") or () if c is not None}
    for cls, a in sorted(result["assignments"].items()):
        if cls not in present:  # 参照表の全行が返るが、データセットに現れたクラスだけを扱う
            continue
        if a["status"] == "ok":
            fid = str(a["feature_id"])
            std = by_id.get(fid)
            classes.append({
                "target_class": cls, "standard_name": a["standard_name"], "feature_id": fid,
                "adduct": a.get("standard_adduct"), "mz": std.mz if std else None,
                "rt": std.rt if std else None, "coefficient": a["coefficient"],
                "via_any_others": cls.strip().casefold() == _ANY_OTHERS,
                **({"dataset_classes": any_others}
                   if cls.strip().casefold() == _ANY_OTHERS else {})})
        elif a["status"] in _UNRESOLVED:
            entry = grouped.setdefault((a["standard_name"], a["status"]), {
                "standard_name": a["standard_name"], "target_classes": [], "status": a["status"],
                "candidates": [{**c, "feature_id": str(c["feature_id"])}
                               for c in a.get("candidates") or []]})
            entry["target_classes"].append(cls)
    unresolved = list(grouped.values())
    missing: Counter = Counter()
    for role, row_statuses in zip(result["roles"], result["statuses"]):
        if role == "internal_standard":
            continue
        missing.update(s for s in row_statuses if s not in _KEEP)
    return round_floats({
        "status": "needs_input" if unresolved else "ready",
        "polarity": polarity, "unit": result.get("unit"),
        "reference": {"sha256": sha, "label": ref.label, "method": ref.method},
        "classes": classes, "unresolved": unresolved, "missing_values": dict(missing)})


def preview(ds, ref: ReferenceInput, bindings: dict | None, polarity: str) -> dict:
    """どの内部標準・係数でどのクラスを割るかを、値を作らずに示す。"""
    _plugin, rows, _text, sha, result = _compute(ds, ref, bindings)
    return _preview_dict(rows, ref, sha, polarity, result)


def _parent_ids(ds) -> list[str]:
    state = getattr(ds, "quant_state", None) or {}
    return [state["record_id"]] if state.get("record_id") else []


def _used_classes(preview_dict: dict) -> list[dict]:
    return [{k: c[k] for k in ("target_class", "standard_name", "feature_id", "coefficient",
                               "via_any_others")} for c in preview_dict["classes"]]


def apply(ds, ref: ReferenceInput, bindings: dict | None, polarity: str, *,
          folder: Path | None = None, reuse: dict | None = None):
    """正規化して記録を書き、派生データセットを返す。`(new, info)`。

    未解決の内部標準が残るときは何も書かず `(None, preview)`。`reuse` は既存の記録（復元用）:
    記録を書かず、その `record_id` で派生を作る。
    """
    plugin, rows, text, sha, result = _compute(ds, ref, bindings)
    pv = _preview_dict(rows, ref, sha, polarity, result)
    if pv["status"] == "needs_input":
        return None, pv
    n_f, n_s = len(rows), len(ds.sample_names)
    assignments = result["assignments"]
    index = {r.feature_id: i for i, r in enumerate(rows)}
    matrix = np.full((n_f, n_s), np.nan)
    mask = np.zeros((n_f, n_s), dtype=bool)
    counts: Counter = Counter()
    for i, row in enumerate(rows):
        if result["roles"][i] == "internal_standard":
            continue
        key = result["target_class_by_feature"][i]
        std = assignments.get(key) or {} if key is not None else {}
        j = index.get(str(std.get("feature_id")))
        std_detected = rows[j].detected if j is not None else (False,) * n_s
        for s in range(n_s):
            status, value = result["statuses"][i][s], result["values"][i][s]
            counts[status] += 1
            if status in _KEEP and value is not None:
                matrix[i, s] = float(value)
                mask[i, s] = row.detected[s] and std_detected[s]
    is_ids = [r.feature_id for r, role in zip(rows, result["roles"]) if role == "internal_standard"]

    reference = {"sha256": sha, "label": ref.label, "method": ref.method, "unit": result.get("unit")}
    if reuse is not None:
        record_id = reuse["record_id"]
    else:
        folder = Path(folder) if folder is not None else analysis_folder(ds)
        qderive.check_capacity()  # 登録できない記録を残さない
        inputs = dict(ds.source_files)
        if ref.path is not None:
            inputs[reference_key(ref.label)] = ref.sha256
        record = write_record(
            folder, "normalize", inputs=inputs, parent_ids=_parent_ids(ds),
            settings={"reference": reference, "rules": rules_settings(plugin),
                      "polarity": polarity, "bindings": _bindings(bindings),
                      "gap_fill": "computed"},
            summary={"status_counts": dict(counts), "unit": result.get("unit"),
                     "n_features": n_f, "n_samples": n_s, "internal_standard_ids": is_ids,
                     "classes": _used_classes(pv)})
        record_id = record["record_id"]
        if ref.concentrations is not None:
            atomic_write_bytes(record_dir(folder) / f"{record_id}.reference.csv",
                               text.encode("utf-8"))
    state = qderive.quant_state(
        normalized=True, method="internal_standard", unit=result.get("unit"),
        record_id=record_id, steps=["normalize"], skipped=[], internal_standard_ids=is_ids)
    new = qderive.derive(ds, matrix=matrix, detected=mask, state=state)
    info = {"record_id": record_id, "status_counts": dict(counts), "unit": result.get("unit"),
            "polarity": polarity, "reference": {k: reference[k] for k in ("sha256", "label", "method")},
            "classes": pv["classes"], "internal_standard_ids": is_ids}
    return new, info


def skip(ds, polarity: str, *, folder: Path | None = None, reuse: dict | None = None):
    """値を変えず、「正規化を飛ばした」記録と派生データセットを作る。規則は読まない。

    `reuse` は既存の記録（復元用）: 記録を書かず、その `record_id` で派生を作る。
    """
    _check_not_normalized(ds)
    if reuse is not None:
        record = reuse
    else:
        folder = Path(folder) if folder is not None else analysis_folder(ds)
        qderive.check_capacity()
        record = write_record(
            folder, "normalize", inputs=dict(ds.source_files), parent_ids=_parent_ids(ds),
            settings={"reference": None, "rules": None, "polarity": polarity},
            summary={"skipped": True}, skipped=True)
    state = qderive.quant_state(normalized=False, method=None, unit=None,
                                record_id=record["record_id"], steps=[], skipped=["normalize"])
    new = qderive.derive(ds, matrix=ds.feature_matrix, detected=ds.detected_mask, state=state)
    return new, {"record_id": record["record_id"], "polarity": polarity, "skipped": True}


def _error(exc: DomainError) -> str:
    return quant_error(exc.code, exc.message, exc.details or None)


def preview_tool(ds, ref: ReferenceInput, bindings, polarity: str) -> str:
    try:
        return json_payload(preview(ds, ref, bindings, polarity))
    except DomainError as exc:
        return _error(exc)


def apply_tool(ds, ref: ReferenceInput | None, bindings, polarity: str, *, holder,
               parent_handle_id: str | None, skip_step: bool = False) -> str:
    """apply / skip の封筒つきの実体（ハンドルの登録まで）。"""
    if ref is None and not skip_step:
        return quant_error("QUANT_BAD_REQUEST",
                           "reference_path は必須です（省略できるのは skip=True のときだけ）。",
                           {"field": "reference_path"})
    try:
        if skip_step:
            new, info = skip(ds, polarity)
        else:
            new, info = apply(ds, ref, bindings, polarity)
            if new is None:
                return json_payload(info)  # needs_input: preview と同じ封筒
        handle_id = qderive.register(holder, new, step="normalize", polarity=polarity,
                                     parent_handle_ids=[parent_handle_id] if parent_handle_id else [])
    except DomainError as exc:
        return _error(exc)
    return json_payload(round_floats({"status": "skipped" if skip_step else "applied",
                                      "dataset_handle": handle_id, **info}))
