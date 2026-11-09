"""`dataset_analysis_run` の v2（analysis-result.v2）の要約と行の表（MCP にもセッションにも依存しない）。

エンジンの結果の全量は `output_path` の JSON に書き、戻り値は要約だけにする（戻り値を肥大させない）。
行を返すのは、どの解析も `MAX_INLINE_ROWS` 以下のときだけ。
"""
from __future__ import annotations

from metabolomix.core.serialization import json_safe, round_floats, truncate_lists

#: output_path なしで行を返す上限（解析ごと）。v1 と同じ。
MAX_INLINE_ROWS = 50
#: 注意・一覧を切る件数（戻り値用）。
LIST_LIMIT = 20

ROW_COLUMNS = ("contrast_id", "feature_id", "status", "reason", "log2fc", "ci_low", "ci_high",
               "p", "q", "verdict")
_VERDICTS = ("changed", "unchanged", "inconclusive")


def tables_of(entry: dict) -> dict[str, list[dict]]:
    """解析 1 件の行の表。compare は 1 つ（`""`）、composition は clr と絶対量、pca は無し。"""
    kind, res = entry["kind"], entry["result"]
    if kind == "compare":
        return {"": res["rows"]}
    if kind == "composition":
        return {name: sub["rows"] for name in ("clr", "absolute")
                if (sub := res.get(name)) is not None}
    return {}


def n_rows(entry: dict) -> int:
    return sum(len(rows) for rows in tables_of(entry).values())


def _counts(rows: list[dict], q_threshold: float) -> dict:
    verdicts = dict.fromkeys(_VERDICTS, 0)
    for r in rows:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    return {"n_rows": len(rows), "n_tested": sum(1 for r in rows if r["status"] == "tested"),
            "verdicts": verdicts,
            "n_q_below": sum(1 for r in rows if r["q"] is not None and r["q"] <= q_threshold)}


def _cut_warnings(warnings: list[dict]) -> list[dict]:
    kept = truncate_lists(json_safe(warnings), LIST_LIMIT,
                          keys=("samples", "sample_ids", "classes", "feature_ids"))
    return kept


def analysis_summary(entry: dict, q_threshold: float) -> dict:
    kind, res = entry["kind"], entry["result"]
    base = {"analysis_id": entry["analysis_id"], "kind": kind}
    if kind == "pca":
        return truncate_lists(json_safe(round_floats({
            **base, "scaling": res["scaling"], "imputation": res["imputation"],
            "n_fit": res["n_fit"], "n_units": len(res["unit_ids"]),
            "n_qc_projected": len(res["qc_ids"]), "n_features": len(res["kept_features"]),
            "explained_variance_ratio": res["explained_variance_ratio"]})), LIST_LIMIT)
    tables = tables_of(entry)
    all_rows = [r for rows in tables.values() for r in rows]
    out = {**base, **_counts(all_rows, q_threshold)}
    if kind == "composition":
        out["tables"] = {name: _counts(rows, q_threshold) for name, rows in tables.items()}
        out["class_label_kind"] = res["class_label_kind"]
        out["n_classes"] = len(res["classes"])
        out["warnings"] = _cut_warnings(res["warnings"])
    else:
        out["test"] = res["test"]
        out["warnings"] = _cut_warnings(res["warnings"])
    return out


def inline_rows(entry: dict, rows_table) -> str:
    """行の TSV（列名 1 回）。composition は先頭に `table` 列（clr / absolute）。

    `rows_table(rows, columns)` は呼び出し側（v1 と共通の表の整形）が渡す。"""
    tables = tables_of(entry)
    composition = entry["kind"] == "composition"
    rows = []
    for name, items in tables.items():
        rows += [{**r, "table": name} for r in items] if composition else items
    columns = (("table",) + ROW_COLUMNS) if composition else ROW_COLUMNS
    return rows_table(rows, columns) if rows else ""


def run_summary(result: dict, *, extra_warnings: list[dict]) -> dict:
    """analysis-result.v2 の要約（行は含めない）。"""
    q = {a["analysis_id"]: a.get("verdict", {}).get("q") for a in result["spec"]["analyses"]}
    analyses = [analysis_summary(e, q.get(e["analysis_id"]) if q.get(e["analysis_id"]) is not None else 0.05) for e in result["analyses"]]
    warnings = _cut_warnings([*result["preprocess"]["warnings"], *extra_warnings])
    return {"schema": result["schema"], "rule_hash": result["rule_hash"],
            "matrix_id": result["matrix_id"], "input_fingerprint": result["input_fingerprint"],
            "primary": result["primary"], "warnings": warnings, "analyses": analyses}
