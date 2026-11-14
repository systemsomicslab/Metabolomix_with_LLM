"""`dataset_compare(output_path=...)` が書く全行の表（全特徴 × 全対比の TSV と、F 検定の TSV）。

差次的エクスポート（`export_contract.py` の 15 列。別リポジトリ massbank-context との契約）とは別物で、
あちらの列・版は変えない。こちらは比較の行（`stats-compare.v1` の `rows`）をそのまま全列で、同定の列を
添えて書く:

- `<output_path>`: 列は `COLUMNS`（対比 × 特徴。行の順はエンジンの行の順）。
- `<stem>.ftests<suffix>`: 多群の omnibus F と要因の F（主効果・交互作用。`stats-compare.v1` の
  `omnibus` / `factorial`）があるときだけ。列は `FTEST_COLUMNS`（`kind` = `omnibus` / `factorial`、
  `effect` は要因の効果の名前（`group` / `group:B` など）で、omnibus では空）。

値の書き方: 数は有効 10 桁（`.10g`）、真偽は `true` / `false`、無い値（None）は空欄、一覧は `,` で結ぶ。
無限大は `inf` と書く（F の `df2` は事前分布の自由度 `df_prior` が ∞ のとき ∞＝χ² の極限で、空欄の「計算していない」と
区別する）。`confounded` は `true`（その対比が完全交絡）/ `false` / 空欄（バッチ不明で評価できない。交絡なしではない）。
同定の列は `feature_metadata`（`.arf2` / mzTab-M SME）を先に、無ければ `feature_annotations`（SML。
MS1 注釈）から取る。`id_level` は `msms` / `ms1` / `none`（`overview.id_level` と同じ 3 分類）、`msi_level`
は道具が付けた MSI の水準（解析フォルダのデータセットだけが持つ。mzTab-M では空欄）。
MCP・セッションには依存しない（`ds` と比較の記録を受け取って書くだけ）。
"""
from __future__ import annotations

import math
from pathlib import Path

from metabolomix.analysis.overview import id_level
from metabolomix.core.atomic_io import atomic_write_bytes

__all__ = ["COLUMNS", "FTEST_COLUMNS", "IDENTITY_COLUMNS", "ftests_path", "identity_of",
           "write_compare_tables"]

IDENTITY_COLUMNS = ("name", "ontology", "inchikey", "id_level", "annotation_source",
                    "msms_assigned", "msi_level")
_STAT_COLUMNS = ("reference", "test", "status", "reason", "n_reference", "n_test",
                 "detected_reference", "detected_test", "quantity_judgeable", "constant_in_group",
                 "confounded", "mean_reference", "mean_test", "log2fc", "se", "ci_low", "ci_high",
                 "ci_unavailable_reason", "statistic", "df", "p", "q", "presence_p", "presence_q",
                 "verdict", "verdict_reason")
COLUMNS = ("contrast_id", "feature_id", *IDENTITY_COLUMNS, *_STAT_COLUMNS)
FTEST_COLUMNS = ("kind", "effect", "feature_id", "name", "id_level", "f", "df1", "df2", "f_p", "f_q")


def identity_of(ds, fid: str) -> dict:
    """特徴 1 つの同定の列（`IDENTITY_COLUMNS` のキー）。"""
    meta = (getattr(ds, "feature_metadata", None) or {}).get(fid) or {}
    ann = (getattr(ds, "feature_annotations", None) or {}).get(fid) or {}
    ident = meta.get("identity") or {}
    return {"name": meta.get("name") or ann.get("name") or "",
            "ontology": meta.get("ontology") or ann.get("ontology") or "",
            "inchikey": meta.get("inchikey") or ann.get("inchikey") or "",
            "id_level": id_level(ds, fid),
            "annotation_source": meta.get("annotation_source") or "",
            "msms_assigned": meta.get("msms_assigned"),
            "msi_level": ident.get("msi_level")}


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return format(value, ".10g") if math.isfinite(value) else ""
    if isinstance(value, (list, tuple)):
        return ",".join(_cell(v) for v in value)
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _tsv(columns, rows) -> bytes:
    lines = ["\t".join(columns)]
    lines.extend("\t".join(_cell(row.get(c)) for c in columns) for row in rows)
    return ("\n".join(lines) + "\n").encode("utf-8")


def ftests_path(path: Path) -> Path:
    """F 検定の表の置き場所（`<stem>.ftests<suffix>`。拡張子が無ければ `.tsv`）。"""
    path = Path(path)
    return path.with_name(path.stem + ".ftests" + (path.suffix or ".tsv"))


def _ftest_rows(engine: dict, ident) -> list[dict]:
    rows = []
    for o in engine.get("omnibus") or []:
        rows.append({"kind": "omnibus", "effect": None, **o, **ident(o["feature_id"])})
    factorial = engine.get("factorial") or {}
    for o in factorial.get("rows") or []:
        rows.append({"kind": "factorial", **o, **ident(o["feature_id"])})
    return rows


def write_compare_tables(path, ds, record: dict) -> dict:
    """比較の記録（`register_compare` の戻り）の全行を TSV に書く。F 検定があれば隣にもう 1 つ。

    戻り値は `{"path", "n_rows", "ftests_path", "n_ftest_rows"}`（F 検定が無ければ `ftests_path` は None）。
    書けなければ OSError をそのまま上げる（呼び出し側が戻り値に理由を載せる）。
    """
    engine = record["engine"]
    cache: dict[str, dict] = {}

    def ident(fid):
        if fid not in cache:
            cache[fid] = identity_of(ds, fid)
        return cache[fid]

    path = Path(path)
    rows = [{**r, **ident(r["feature_id"])} for r in engine["rows"]]
    atomic_write_bytes(path, _tsv(COLUMNS, rows))
    out = {"path": str(path), "n_rows": len(rows), "ftests_path": None, "n_ftest_rows": 0}
    ftests = _ftest_rows(engine, ident)
    if ftests:
        fpath = ftests_path(path)
        atomic_write_bytes(fpath, _tsv(FTEST_COLUMNS, ftests))
        out.update({"ftests_path": str(fpath), "n_ftest_rows": len(ftests)})
    return out
