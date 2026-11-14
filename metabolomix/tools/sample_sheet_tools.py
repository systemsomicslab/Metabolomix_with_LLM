"""生データのフォルダに実験情報シートを書く MCP ツール: pipeline_sample_sheet。

ISA の試料の表（`s_*.txt` と `a_*.txt`）か mwTab（`SUBJECT_SAMPLE_FACTORS`）から、因子の列付きの
sample-manifest.v2 を raw のフォルダに書く（組み立ては `metabolomix/rawmeta/sample_table.py`）。pipeline の
受付（`pipeline_plan` / `pipeline_run`）は、要求で別のシートを指さなければ raw のフォルダの
`sample-manifest.tsv` を読む。書く前に、pipeline の受付と同じ厳密な読み手（`sample_manifest.parse_manifest`）で
読み直して確かめる。
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core.atomic_io import DomainError, atomic_write_bytes
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload, truncate_lists

__all__ = ["pipeline_sample_sheet"]

_LIST_LIMIT = 20
_LIST_KEYS = ("unreadable", "qc_samples", "blank_samples", "unassigned_samples", "files_not_in_assay")


def _envelope(exc: DomainError) -> str:
    details = truncate_lists(exc.details or {}, _LIST_LIMIT,
                             keys=tuple((exc.details or {}).keys())) if exc.details else None
    return mztab_error(exc.code, exc.message, details)


def _bad(message: str, **details) -> DomainError:
    return DomainError("SAMPLE_TABLE_INVALID", message, details or None)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True),
          structured_output=False)
def pipeline_sample_sheet(
    raw_dir: str,
    sample_table: str,
    factor_columns: list[str],
    assay_table: str | None = None,
    sample_table_format: str = "auto",
    qc_selector: list[str] | None = None,
    blank_selector: list[str] | None = None,
    keep_extension: str | None = None,
    output_name: str = "sample-manifest.tsv",
    overwrite: bool = False,
) -> str:
    """ISA / mwTab の試料の表から、因子の列付きの実験情報シート（sample-manifest.v2）を raw のフォルダに書く。

    群は `factor_columns`（試料の表の因子の列。`Factor Value[名]`・`Characteristics[名]` の名前でも引ける）の
    値の組（値を `_` でつなぐ）。値が空・表に行の無い試料の群は `unassigned`。因子の列はシートの後ろに足す
    （`dataset_compare(factors=...)` がその列で計画を組める）。注入順は raw の取得時刻（メタデータだけを読む）の順で、
    1 つでも読めなければ空にする（`injection_order_source: "none"` と読めなかったファイルを返す）。バッチは空。
    1 行は raw のフォルダの主ファイル 1 つ（pipeline の受付と同じ一覧。形式が混ざるときは `keep_extension`）。

    raw_dir: 生データのフォルダ（pipeline の `raw_dir`）。sample_table: ISA の `s_*.txt` か mwTab の
        ファイル。assay_table: ISA の `a_*.txt`（`Sample Name` と `Raw Spectral Data File` でファイルと試料を
        結ぶ）。省けば mwTab の `RAW_FILE_NAME`、それも無ければファイルの幹を試料名とみなす。
    sample_table_format: "auto"（既定。1 行目が `#METABOLOMICS WORKBENCH` なら mwTab）/ "isa" / "mwtab"。
    qc_selector / blank_selector: [列名, 正規表現]。当たった試料を role=qc（群 QC）/ role=blank（群 blank）に
        する。列が `Sample Name` なら表に行の無い試料（ファイルの幹）にも当てる。
    output_name: 書くファイル名（既定 `sample-manifest.tsv`。pipeline が既定で読む名前）。既にあれば
        `SAMPLE_SHEET_EXISTS`（`overwrite=True` で書き直す）。
    """
    from metabolomix.analysis import sample_manifest
    from metabolomix.pipeline.inputs import list_primary_names
    from metabolomix.rawmeta.acquisition import gaps_over
    from metabolomix.rawmeta.sample_table import build_sample_sheet

    try:
        root = Path(raw_dir).expanduser()
        if not root.is_dir():
            raise _bad(f"raw のフォルダが無い: {root}", raw_dir=str(root))
        if (not isinstance(output_name, str) or not output_name.strip()
                or Path(output_name).name != output_name or output_name in (".", "..")):
            raise _bad("output_name はフォルダを含まないファイル名にする", output_name=output_name)
        if not isinstance(overwrite, bool):
            raise _bad("overwrite は真偽値", overwrite=repr(overwrite))
        target = root / output_name
        if target.exists() and not overwrite:
            raise DomainError("SAMPLE_SHEET_EXISTS",
                              f"{output_name} が既にある。書き直すなら overwrite=True を渡す",
                              {"path": str(target)})
        names = list_primary_names(root, {"keep_extension": keep_extension})
        built = build_sample_sheet(
            root, names, Path(sample_table).expanduser(), factor_columns,
            assay_table=None if assay_table is None else Path(assay_table).expanduser(),
            sample_table_format=sample_table_format, qc_selector=qc_selector,
            blank_selector=blank_selector)
        # pipeline の受付と同じ読み手で確かめてから書く（読めないシートを置かない）
        with tempfile.TemporaryDirectory() as tmp:
            probe = Path(tmp) / output_name
            probe.write_text(built["text"], encoding="utf-8")
            sample_manifest.parse_manifest(probe, source_root=root, expected_sources=names)
        atomic_write_bytes(target, built["text"].encode("utf-8"))
    except DomainError as exc:
        return _envelope(exc)
    except OSError as exc:
        return mztab_error("SAMPLE_TABLE_INVALID", f"読み書きできない: {exc}", None)

    notes = ["batch は空（分からない）。分かっていれば batch 列を埋める。"]
    if built["injection_order_source"] != "raw_acquisition_time":
        notes.append("取得時刻の読めない raw があるので injection_order は空にした（ファイル名の順は注入順とは限らない）。")
    if built["unassigned_samples"]:
        notes.append("因子の値が空か試料の表に行の無い試料は群 unassigned にした。")
    summary = {
        "status": "success",
        "path": str(target),
        "n_rows": len(built["rows"]),
        "columns": built["columns"],
        "groups": built["groups"],
        "factor_columns": built["factor_columns"],
        "sample_table_format": built["sample_table_format"],
        "injection_order_source": built["injection_order_source"],
        "gaps_over_40min": gaps_over(built["acquisition"]),
        "unreadable": built["unreadable"],
        "qc_samples": built["qc_samples"],
        "blank_samples": built["blank_samples"],
        "unassigned_samples": built["unassigned_samples"],
        "files_not_in_assay": built["files_not_in_assay"],
        "notes": notes,
        "next": (f"pipeline_plan / pipeline_run は raw_dir の {output_name} を読む（既定名でなければ要求の "
                 "sample_manifest で指す）。読み込み済みのデータセットには dataset_set_sample_metadata で当てる"),
    }
    return json_payload(truncate_lists(summary, _LIST_LIMIT, keys=_LIST_KEYS))
