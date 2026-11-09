"""定量の鎖の復元（再開の部品）。

解析フォルダの `reports/quant/` の記録から、古くない最新の鎖（正規化 → 統合前の選別）を、記録の
設定から計算し直して派生データセットとして開く。記録も audit も書かない（記録は既にある）。

- 参照表は設定にも記録にも場所を残さない。記録の `settings.reference.sha256` と一致する表を、
  保存した `<record_id>.reference.csv`（濃度から作った係数表）と、呼び出しで渡された
  `reference_paths` の中から sha256 で探す。見つからなければ `reference_missing` でその工程から先を戻さない。
- 規則プラグインが未設定なら `rules_unavailable`、`digest` が記録と違えば `rules_digest_changed`
  （古い扱い）で、その工程から先を戻さない。
- 統合の記録と、統合後の選別は v1 では戻さない（両極性のデータセットが要る）。
  `merge_restore_needs_both_polarities` として skipped に挙げる。
- 選別の記録の `settings.curation_digest` が、いまのキュレーション（`flags.jsonl`）の digest と
  違えば `curation_changed`（古い扱い）。一致するときだけ、いまのフラグで計算し直す。
- 途中の工程を戻せなければ、鎖のそれより後ろの記録は `prior_step_not_restored` として
  skipped に挙げる（黙って消さない）。
- `load_dataset` からの呼び出しは B-2 の入口と衝突するので、統合時に足す。
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.quant import derive as qderive
from metabolomix.quant import normalize as qnormalize
from metabolomix.quant import select as qselect
from metabolomix.quant.inputs import ReferenceInput, adduct_table, adduct_key
from metabolomix.quant.plugin import load_rules
from metabolomix.quant.records import is_stale, read_records, record_dir

__all__ = ["restore_chain"]

MERGE_REASON = "merge_restore_needs_both_polarities"
CURATION_CHANGED = "curation_changed"
AFTER_BREAK = "prior_step_not_restored"
_REFERENCE_PREFIX = "reference:"


def _entry(record: dict, reason: str | None = None, **extra) -> dict:
    out = {"record_id": record["record_id"], "step": record["step"]}
    if reason is not None:
        out["reason"] = reason
    out.update(extra)
    return out


def _find_reference(folder: Path, record: dict, reference_paths) -> ReferenceInput | None:
    """記録の参照表（sha256 が一致するもの）を保存した表と渡された表から探す。"""
    wanted = ((record.get("settings") or {}).get("reference") or {})
    sha = wanted.get("sha256")
    if not sha:
        return None
    candidates = [record_dir(folder) / f"{record['record_id']}.reference.csv",
                  *(Path(p) for p in reference_paths)]
    for path in candidates:
        try:
            data = path.read_bytes()
        except OSError:
            continue
        if hashlib.sha256(data).hexdigest() == sha:
            return ReferenceInput(
                text=data.decode("utf-8-sig"), sha256=sha, label=str(wanted.get("label") or path.name),
                path=None, method=wanted.get("method") or "cell_coefficient", concentrations=None,
                solvent_volume=None, cells_millions=None, unit=None)
    return None


def _current_adduct_inputs() -> dict[str, str]:
    try:
        table = adduct_table()
    except DomainError:
        return {}
    return {adduct_key(table.label): table.sha256} if table is not None else {}


def _curation_reason(record: dict, base_ds) -> str | None:
    """選別の記録のキュレーションが、いまの判断と違えば理由（同じ・記録に digest が無ければ None）。"""
    if record.get("step") != "select":
        return None
    recorded = (record.get("settings") or {}).get("curation_digest")
    if not isinstance(recorded, dict):
        return None
    try:
        current = qselect.current_curation_digests(base_ds, recorded, merged=False)
    except DomainError:
        return "curation_unreadable"
    return CURATION_CHANGED if current != recorded else None


def _classify(records: list[dict], base_ds, out: dict) -> list[dict]:
    """古くない・統合に絡まない記録を作った順に返す（古いものは stale、統合は skipped へ）。"""
    current = dict(base_ds.source_files or {})
    if any(k.startswith("adducts:") for r in records for k in r.get("inputs", {})):
        current.update(_current_adduct_inputs())
    merge_ids: set[str] = set()
    known: set[str] = set()
    stale_ids: set[str] = set()
    live: list[dict] = []
    for record in records:
        rid = record["record_id"]
        if record.get("step") == "merge" or any(p in merge_ids for p in record.get("parent_ids", [])):
            merge_ids.add(rid)
            out["skipped"].append(_entry(record, MERGE_REASON))
            continue
        # 参照表は sha256 で別に解く（濃度形では記録の入力の指紋が元の表のもので、保存した係数表とは違う）
        view = {**record, "inputs": {k: v for k, v in record.get("inputs", {}).items()
                                     if not k.startswith(_REFERENCE_PREFIX)}}
        reason = is_stale(view, current, known, stale_ids=stale_ids) or _curation_reason(record, base_ds)
        if reason is not None:
            stale_ids.add(rid)
            out["stale"].append(_entry(record, reason))
            continue
        known.add(rid)
        live.append(record)
    return live


def _chain_to(live: list[dict]) -> list[dict]:
    """最新の記録から親をたどった鎖（根から順）。"""
    if not live:
        return []
    by_id = {r["record_id"]: r for r in live}
    chain, current = [], live[-1]
    while current is not None:
        chain.append(current)
        parents = current.get("parent_ids") or []
        current = by_id.get(parents[0]) if parents else None
    return chain[::-1]


def _restore_step(record: dict, ds, folder: Path, reference_paths):
    """1 工程を計算し直す。`(new, None)` か `(None, (out のキー, 理由))`。"""
    settings = record.get("settings") or {}
    step = record["step"]
    needs_rules = step == "select" or not record.get("skipped")
    if needs_rules:
        try:
            plugin = load_rules()
        except DomainError:
            return None, ("skipped", "rules_unavailable")
        recorded = ((settings.get("rules") or {}).get("digest"))
        if recorded != plugin.digest:
            return None, ("stale", "rules_digest_changed")
    try:
        if step == "normalize":
            polarity = settings.get("polarity")
            if record.get("skipped"):
                new, _info = qnormalize.skip(ds, polarity, reuse=record)
                return new, None
            ref = _find_reference(folder, record, reference_paths)
            if ref is None:
                return None, ("skipped", "reference_missing")
            new, _info = qnormalize.apply(ds, ref, settings.get("bindings"), polarity, reuse=record)
            if new is None:
                return None, ("skipped", "needs_input")
            return new, None
        if step == "select":
            if qselect.is_merged(ds):
                return None, ("skipped", MERGE_REASON)
            adduct = adduct_table() if settings.get("adducts") else None
            if settings.get("adducts") and adduct is None:
                return None, ("skipped", "adduct_table_unavailable")
            new, _info = qselect.select(
                ds, adduct=adduct, apply_curation=bool(settings.get("apply_curation", True)),
                polarity=settings.get("polarity"), reuse=record)
            return new, None
    except DomainError as exc:
        if exc.code == qselect.CURATION_CHANGED:  # 分類のあとで判断が変わった
            return None, ("stale", CURATION_CHANGED)
        return None, ("skipped", exc.code)
    return None, ("skipped", f"unsupported_step:{step}")


def restore_chain(folder: Path, base_ds, base_holder, *, fresh: bool = False,
                  reference_paths=()) -> dict:
    """古くない最新の鎖を計算し直して派生ハンドルを作る。

    `{"restored": [{record_id, step, dataset_handle}], "stale": [...], "skipped": [...]}`。
    `fresh=True` は記録を見ない。途中で戻せなければ、その工程から先は戻さず理由を書く。
    """
    out: dict = {"restored": [], "stale": [], "skipped": []}
    if fresh:
        return out
    folder = Path(folder)
    records, warnings = read_records(folder)
    if warnings:
        out["warnings"] = warnings
    chain = _chain_to(_classify(records, base_ds, out))
    ds, parent = base_ds, getattr(base_holder, "handle_id", None)
    for position, record in enumerate(chain):
        new, problem = _restore_step(record, ds, folder, reference_paths)
        if problem is not None:
            out[problem[0]].append(_entry(record, problem[1]))
            out["skipped"].extend(_entry(later, AFTER_BREAK) for later in chain[position + 1:])
            break
        handle_id = qderive.register(base_holder, new, step=record["step"],
                                     parent_handle_ids=[parent] if parent else [])
        out["restored"].append(_entry(record, dataset_handle=handle_id))
        ds, parent = new, handle_id
    return out
