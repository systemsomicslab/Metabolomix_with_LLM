"""統計エンジンの行列の登録と前処理の記録（`dataset_preprocess` / `dataset_compare` / `dataset_pca` /
`pathway_project` が共有する、ツールではない部品。`tools/_holder.py` と同じ位置づけ）。

行列は `ds.stat_matrices` へ登録し、前処理の記録を `reports/stats/` に 1 件追記する
（書けないときも登録は成功し、`record.written=false` で理由を返す）。
"""
from __future__ import annotations

from metabolomix.analysis import stats_records, stats_results
from metabolomix.core import mcp_core
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.serialization import json_safe, round_floats, truncate_lists

__all__ = ["store_matrix", "register_default", "preprocess_summary"]

_STEP_LIST_LIMIT = 20
#: 要約に載せる文字列の長さの上限。超えたら上限の長さで切り、末尾を「…」にする（黙って落とさない）。
_STR_LIMIT = 80


def _compact(value: dict) -> dict:
    """工程の dict を小さな要約にする（入れ子も）。スカラーはそのまま、一覧は件数（`n_<key>`）、
    値がすべて一覧の dict（`excluded_by_role` など）は `{キー: 件数}`、入れ子の dict（`filters` など）は
    再帰する。特徴 ID の一覧（`dropped_ids`）は件数が `dropped` にあるので載せない。文字列は
    `_STR_LIMIT` 字までで、長いものは `_STR_LIMIT` 字に切って末尾を「…」にする。"""
    if value and all(isinstance(v, (list, tuple, set, frozenset)) for v in value.values()):
        return {k: len(v) for k, v in value.items()}
    out: dict = {}
    for key, item in value.items():
        if key == "dropped_ids":
            continue
        if isinstance(item, bool) or item is None or isinstance(item, (int, float)):
            out[key] = item
        elif isinstance(item, str):
            out[key] = item if len(item) <= _STR_LIMIT else item[:_STR_LIMIT - 1] + "…"
        elif isinstance(item, (list, tuple, set, frozenset)):
            out[f"n_{key}"] = len(item)
        elif isinstance(item, dict):
            if len(item) > 8 and all(isinstance(v, (int, float)) for v in item.values()):
                out[f"n_{key}"] = len(item)
                out[f"sum_{key}"] = sum(item.values())
            else:
                out[key] = _compact(item)
    return out


def _step_summary(step: dict) -> dict:
    """工程 1 つの要約: 名前・状態・理由を先に、あとに `_compact` した件数と足切りごとの内訳。"""
    head = {"step": step.get("step"), "status": step.get("status")}
    if step.get("reason") is not None:
        head["reason"] = step["reason"]
    rest = {k: v for k, v in step.items() if k not in ("step", "status", "reason")}
    return {**head, **_compact(rest)}


def _group_counts(prepared) -> dict:
    counts: dict[str, int] = {}
    for unit in prepared.units:
        label = unit.group if unit.group is not None else "(none)"
        counts[label] = counts.get(label, 0) + 1
    return dict(sorted(counts.items()))


def _warning_summary(warnings: list[dict]) -> list[dict]:
    """注意の要約（コードと、短い欄だけ。試料名の長い一覧は切る）。"""
    return truncate_lists(json_safe(warnings), _STEP_LIST_LIMIT, keys=("samples", "sample_ids"))


def _write_record(holder, ds, prepared, config, inputs, summary) -> dict:
    """記録を書く。書き先が決まらない・書けないときは `{"written": False, "reason": ...}`。"""
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        parents = []
        record_id = (getattr(ds, "quant_state", None) or {}).get("record_id")
        if record_id:
            parents.append(record_id)
        written = stats_records.write_stats_record(
            reports_dir, "preprocess",
            inputs={"dataset_hash": prepared.dataset_hash, "matrix_id": prepared.matrix_id},
            parent_ids=parents,
            settings={"config": config, "group_overrides": inputs["group_overrides"],
                      "exclusions": {"samples": inputs["exclude_samples"],
                                     "features": inputs["exclude_features"]}},
            summary=summary)
    except DomainError as exc:
        return {"written": False, "reason": f"{exc.code}: {exc.message}"}
    except OSError as exc:
        return {"written": False, "reason": f"OSError: {exc}"}
    return {"written": True, "record_id": written["record_id"], "path": written["path"]}


def preprocess_summary(study, prepared) -> dict:
    """前処理の戻り値の要約（工程・群ごとの個体数・注意。記録の summary にもなる）。"""
    quant = prepared.quant
    groups = _group_counts(prepared)
    warnings = list(prepared.warnings)
    if len(groups) == 1:
        warnings.append({"code": "single_group", "group": next(iter(groups)),
                         "message": "すべての個体が 1 つの群に入っている（比較できない）。"
                                    "group_levels / group_factors の綴りか試料の表の群を確かめてください。"})
    summary = {
        "matrix_id": prepared.matrix_id,
        "n_units": len(prepared.units),
        "n_qc": len(prepared.qc_ids),
        "n_features_before": len(study.feature_ids),
        "n_features": len(prepared.feature_ids),
        "subject_source": prepared.subject_source,
        "groups": groups,
        "quant": ({"normalized": quant.normalized, "method": quant.method, "unit": quant.unit}
                  if quant is not None else None),
        "steps": [_step_summary(s) for s in prepared.steps],
        "warnings": _warning_summary(warnings),
    }
    return truncate_lists(json_safe(round_floats(summary)), _STEP_LIST_LIMIT)


def store_matrix(holder, ds, study, prepared, config, inputs):
    """行列を登録し、前処理の記録を 1 件書く。`(登録, 要約, record 欄)` を返す。

    記録が書けたときは登録に `record_id` を持たせる（後続の比較の記録がここを親にする）。
    別のデータから作った行列は DomainError(STUDY_DATA_INVALID)。
    """
    try:
        entry = stats_results.register_matrix(
            ds, prepared, config, inputs["group_overrides"],
            exclusions={"samples": inputs["exclude_samples"],
                        "features": inputs["exclude_features"]})
    except ValueError as exc:
        raise DomainError("STUDY_DATA_INVALID", str(exc)) from exc
    summary = preprocess_summary(study, prepared)
    record = _write_record(holder, ds, prepared, config, inputs, summary)
    if record.get("written"):
        entry["record_id"] = record["record_id"]
    return entry, summary, record


def register_default(holder, ds, study, prepared, config, inputs):
    """既定の前処理の行列を登録する（書き込みを伴う部分）。同じ行列が登録済みなら再利用する。
    `(登録, 前処理の記録 or None)` を返す（記録を書いたのは新しく登録したときだけ）。"""
    entry = stats_results.get_matrix(ds, prepared.matrix_id)
    if entry is not None:
        return entry, None
    entry, _summary, record = store_matrix(holder, ds, study, prepared, config, inputs)
    return entry, record
