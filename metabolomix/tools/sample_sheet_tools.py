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
_LIST_KEYS = ("unreadable", "qc_samples", "blank_samples", "unassigned_samples", "files_not_in_assay",
              "ambiguous_files", "conflicting_samples")


def _envelope(exc: DomainError) -> str:
    details = truncate_lists(exc.details or {}, _LIST_LIMIT,
                             keys=tuple((exc.details or {}).keys())) if exc.details else None
    return mztab_error(exc.code, exc.message, details)


def _capped_reasons(reasons: dict[str, int]) -> tuple[dict[str, int], int]:
    """読めなかった理由 → 件数を、件数の多い順（同数は先に出た順）に先頭 `_LIST_LIMIT` 種類だけ返す。

    理由はファイルによらない文面にしてある（`rawmeta.acquisition`）が、OS やベンダーの例外の文面が
    ファイルごとに違うこともあるので、一覧と同じ上限で切る。2 つ目は切った理由の種類の数。
    """
    ranked = sorted(reasons.items(), key=lambda item: -item[1])
    return dict(ranked[:_LIST_LIMIT]), max(len(ranked) - _LIST_LIMIT, 0)


def _bad(message: str, **details) -> DomainError:
    return DomainError("SAMPLE_TABLE_INVALID", message, details or None)


def _notes(built: dict) -> list[str]:
    """シートを組んだ結果のうち、利用者に確かめてもらうこと。"""
    notes = ["batch は空（分からない）。分かっていれば batch 列を埋める。"]
    source, n_unreadable = built["injection_order_source"], len(built["unreadable"])
    if source == "none":
        notes.append("どの raw も取得時刻を読めないので injection_order はすべて空にした"
                     "（ファイル名の順は注入順とは限らない）。unreadable_reasons に理由。")
    elif n_unreadable:
        note = (f"取得時刻の読めない raw {n_unreadable} 件（unreadable）は injection_order を空にし、読めた raw は"
                "取得時刻の順（読めたものの中の順位）にした。空欄は不明のまま扱われる（MS-DIAL の番号では埋めない）。")
        unread = set(built["unreadable"])
        # 前処理は include=false の試料を先に外す（`_select`）ので、欠けた順が補正と判定を止めるのは統計に入る試料だけ
        missing = [r for r in built["rows"] if r["source_file"] in unread and r["role"] in ("sample", "qc")]
        if any(r["include"] == "true" for r in missing):
            note += ("生体試料か QC の順が欠けるので、ドリフト補正と注入順の偏りの判定は、その順を埋めるまで"
                     "行われない（分かれば injection_order を埋める）。")
        elif missing:
            note += ("順の欠けた生体試料は統計から外した試料（include=false）だけなので、ドリフト補正と注入順の偏りの"
                     "判定は止まらない（その試料を include=true にするなら injection_order も埋める）。")
        notes.append(note)
    n_unassigned = len(built["unassigned_samples"])
    if n_unassigned:
        notes.append(f"群の決まらない生体試料 {n_unassigned} 件（unassigned_samples）は群 unassigned・include=false に"
                     "した（統計に入れない）。因子の値を埋めて include を true にすれば使える。")
    if built["files_not_in_assay"]:
        notes.append(f"試料の表のファイル名（と幹）に当たらない raw が {len(built['files_not_in_assay'])} 件ある"
                     "（files_not_in_assay）。")
    if built["ambiguous_files"]:
        notes.append("別々の試料に当たる raw がある（ambiguous_files）。どちらの群にも入れていない。")
    if built["conflicting_samples"]:
        notes.append("試料の表に同じ試料の行が複数あり因子の値が食い違う（conflicting_samples）。どちらの群にも入れていない。")
    for given, written in built["renamed_factor_columns"].items():
        notes.append(f"因子の列 {given!r} はシートの既知の列と同じ名前なので {written!r} として書いた"
                     f"（dataset_compare(factors=...) にはこの名前を渡す）。")
    return notes


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
    値の組（値を `_` でつなぐ）。群の決まらない生体試料（因子の値が空・表に行が無い・1 つのファイルが別々の試料に
    当たる・同じ試料の行どうしで値が食い違う）は群 `unassigned`・`include=false`（統計の母集団から外す。値を埋めて
    `include` を true にすれば使える）で、`unassigned_samples`・`files_not_in_assay`・`ambiguous_files`・
    `conflicting_samples` と `notes` に出す。生体試料のどれも群が決まらなければ書かずに `SAMPLE_TABLE_INVALID`
    （ファイルと試料の対応か因子の列が合っていない）。因子の列はシートの後ろに渡した名前で足す
    （`dataset_compare(factors=...)` がその列で計画を組める）。シートの既知の列と**完全に同じ**名前（`group`・`batch`
    など）だけは `factor_<名>` に改めて書き、`renamed_factor_columns` と `notes` に出す（`Group` のように大文字小文字
    だけ違う名前はそのまま）。注入順は raw の取得時刻（メタデータだけを読む。`.wiff2` は SCIEX OS が並べて書く
    同じ幹の `.wiff` から読む）の順で、シートの schema の行の直後に
    `# injection_order_source = raw_acquisition_time` を書く（読み込んだ試料の表・統計の注入順の出所が
    `raw_acquisition_time` になる）。読めない raw があっても読めた raw の順は残し、読めない raw の注入順だけを空に
    する（`injection_order_source: "raw_acquisition_time_partial"`。空欄は MS-DIAL の番号で埋まらず不明のまま）。
    どれも読めなければ `none`（宣言も書かない）。読めなかったファイルは `unreadable`、理由は `unreadable_reasons`
    （理由 → 件数。件数の多い順に先頭 20 種類で、切った種類の数は `unreadable_reasons_n_more`。一覧も先頭 20 件で、
    切った件数は `<一覧>_n_more`）。順の欠けた生体試料か QC が統計に入る（`include=true`）ときだけ、`notes` がドリフト補正と注入順の偏りの判定は
    行われないと書く（`include=false` の試料は前処理が先に外すので止めない）。
    バッチは空。試料の表は UTF-8（BOM 可）。
    1 行は raw のフォルダの主ファイル 1 つ（pipeline の受付と同じ一覧。形式が混ざるときは `keep_extension`）。

    raw_dir: 生データのフォルダ（pipeline の `raw_dir`）。sample_table: ISA の `s_*.txt` か mwTab の
        ファイル。assay_table: ISA の `a_*.txt`（`Sample Name` と `Raw Spectral Data File`・
        `Derived Spectral Data File` でファイルと試料を結ぶ）。省けば mwTab の `RAW_FILE_NAME`、それも無ければ
        ファイルの幹を試料名とみなす。表のファイル名と合わなければ幹（拡張子を除いた名前）で引く。
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

    notes = _notes(built)
    reasons, n_more_reasons = _capped_reasons(built["unreadable_reasons"])
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
        "unreadable_reasons": reasons,
        **({"unreadable_reasons_n_more": n_more_reasons} if n_more_reasons else {}),
        "unassigned_samples": built["unassigned_samples"],
        "files_not_in_assay": built["files_not_in_assay"],
        "ambiguous_files": built["ambiguous_files"],
        "conflicting_samples": built["conflicting_samples"],
        "renamed_factor_columns": built["renamed_factor_columns"],
        "notes": notes,
        "next": (f"pipeline_plan / pipeline_run は raw_dir の {output_name} を読む（既定名でなければ要求の "
                 "sample_manifest で指す）。読み込み済みのデータセットには dataset_set_sample_metadata で当てる"),
    }
    return json_payload(truncate_lists(summary, _LIST_LIMIT, keys=_LIST_KEYS))
