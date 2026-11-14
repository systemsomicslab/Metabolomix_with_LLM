"""`dataset_analysis_run` の v2（analysis-result.v2）の要約と行の表（MCP にもセッションにも依存しない）。

エンジンの結果の全量は `output_path` の JSON に書き、戻り値は要約だけにする（戻り値を肥大させない）。
`output_path` が無ければ解析ごとに行を TSV で返し、`MAX_INLINE_ROWS` を超える解析は p の小さい順の
先頭 `MAX_INLINE_ROWS` 行だけを載せる（断らない。`first_by_p`）。
"""
from __future__ import annotations

import math

from metabolomix.core.serialization import json_safe, round_floats, truncate_lists

#: output_path なしで行を返す上限（解析ごと）。v1 と同じ。
MAX_INLINE_ROWS = 50
#: 注意・一覧を切る件数（戻り値用）。
LIST_LIMIT = 20

ROW_COLUMNS = ("contrast_id", "feature_id", "status", "reason", "log2fc", "ci_low", "ci_high",
               "p", "q", "verdict")
_VERDICTS = ("changed", "unchanged", "inconclusive")


#: plsda の要約に載せる VIP の上位の数（全量は output_path の JSON）。
VIP_TOP = 10


def tables_of(entry: dict) -> dict[str, list[dict]]:
    """解析 1 件の行の表。compare は 1 つ（`""`）、composition は clr と絶対量、pca・plsda は無し。"""
    kind, res = entry["kind"], entry["result"]
    if kind == "compare":
        return {"": res["rows"]}
    if kind == "composition":
        return {name: sub["rows"] for name in ("clr", "absolute")
                if (sub := res.get(name)) is not None}
    return {}


def n_rows(entry: dict) -> int:
    return sum(len(rows) for rows in tables_of(entry).values())


def first_by_p(rows: list[dict], limit: int) -> list[dict]:
    """行が `limit` を超えれば、p の小さい順（p の無い・非有限の行は後ろ、同じなら元の順）の先頭 `limit` 行。

    `limit` 以下なら元の順のまま全部。`dataset_analysis_run` の v1 と v2 が、`output_path` の無い戻り値を
    切るのに使う（どの行を載せたかを `rows_order` で返す）。
    """
    if len(rows) <= limit:
        return list(rows)

    def key(index: int):
        p = rows[index].get("p")
        usable = isinstance(p, (int, float)) and not isinstance(p, bool) and math.isfinite(p)
        return (not usable, p if usable else 0.0, index)

    return [rows[i] for i in sorted(range(len(rows)), key=key)[:limit]]


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


def confounding_summary(results) -> dict | None:
    """比較の結果（`stats-compare.v1`）の交絡 `{batch_status, by_contrast}`。最初に交絡を持つ結果のもの。

    組成は内側の比較（`absolute` / `clr`）を渡す。どちらも同じ個体で組むので同じ値で、`absolute` を先に読む
    （CLR が計算できないときも絶対量はある）。交絡を持つ結果が無ければ None。既定の `allow_confounded=true` は
    完全交絡でも比べるので、交絡は要約に印として残す（止める代わり）。"""
    for sub in results:
        conf = sub.get("confounding")
        if conf:
            return {"batch_status": conf.get("batch_status"), "by_contrast": conf.get("by_contrast")}
    return None


def composition_warnings(own: list[dict], inner: dict[str, dict]) -> list[dict]:
    """組成の注意: 組成そのもの（`own`）と内側の比較（`absolute` / `clr`）の注意をコードで重ねずにまとめる。

    内側の注意（完全交絡・注入順の偏り・片群一定ほか）は、どの表から来たかを `tables` に付ける。同じコードが
    両方の表にあれば先に現れたもの（`absolute`）を残し、`tables` に両方を並べる。"""
    out = [dict(w) for w in own]
    seen = {w.get("code"): w for w in out}
    for name, sub in inner.items():
        for w in sub.get("warnings") or []:
            code = w.get("code")
            if code in seen:
                tables = seen[code].get("tables")
                if tables is not None and name not in tables:
                    tables.append(name)
                continue
            item = {**w, "tables": [name]}
            out.append(item)
            seen[code] = item
    return out


def factorial_summary(res: dict, q_threshold: float) -> dict | None:
    """要因の F（`stats-compare.v1` の `factorial`）の効果ごとの件数。無ければ None。"""
    factorial = res.get("factorial")
    if not factorial:
        return None
    out = {}
    for name, eff in factorial["effects"].items():
        rows = [r for r in factorial["rows"] if r["effect"] == name and r["f_p"] is not None]
        out[name] = {"df1": eff["df1"], "n_tested": len(rows),
                     "n_q_below": sum(1 for r in rows if r["f_q"] is not None and r["f_q"] <= q_threshold)}
    return out


def _plsda_summary(base: dict, res: dict) -> dict:
    return truncate_lists(json_safe(round_floats({
        **base, "classes": res["classes"], "n_units": res["n_units"], "n_features": res["n_features"],
        "scaling": res["scaling"], "imputation": res["imputation"],
        "n_components": res["n_components"], "r2x_cum": res["r2x_cum"], "r2y_cum": res["r2y_cum"],
        "q2_cum": res["q2_cum"], "q2": res["q2"], "cv": res["cv"],
        "permutation": res["permutation"], "vip_top": res["vip"][:VIP_TOP],
        "warnings": res["warnings"]})), LIST_LIMIT, keys=("unit_ids",))


def analysis_summary(entry: dict, q_threshold: float) -> dict:
    kind, res = entry["kind"], entry["result"]
    base = {"analysis_id": entry["analysis_id"], "kind": kind}
    if kind == "plsda":
        return _plsda_summary(base, res)
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
        inner = {name: sub for name in ("absolute", "clr") if (sub := res.get(name)) is not None}
        confounding = confounding_summary(inner.values())
        if confounding is not None:
            out["confounding"] = confounding
        out["warnings"] = _cut_warnings(composition_warnings(res["warnings"], inner))
    else:
        out["test"] = res["test"]
        out["p_adjust_method"] = res.get("p_adjust_method", "bh")
        factorial = factorial_summary(res, q_threshold)
        if factorial is not None:
            out["factorial"] = factorial
        confounding = confounding_summary([res])
        if confounding is not None:
            out["confounding"] = confounding
        out["warnings"] = _cut_warnings(res["warnings"])
    return out


def inline_rows(entry: dict, rows_table, limit: int | None = None) -> str:
    """行の TSV（列名 1 回）。composition は先頭に `table` 列（clr / absolute）。

    `rows_table(rows, columns)` は呼び出し側（v1 と共通の表の整形）が渡す。`limit` を渡し、行がそれを
    超えれば `first_by_p` の先頭 `limit` 行だけ（composition は 2 つの表を合わせて選ぶ）。"""
    tables = tables_of(entry)
    composition = entry["kind"] == "composition"
    rows = []
    for name, items in tables.items():
        rows += [{**r, "table": name} for r in items] if composition else items
    if limit is not None:
        rows = first_by_p(rows, limit)
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
