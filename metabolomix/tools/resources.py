"""MCP リソース（@mcp.resource ×2。静的 1・テンプレート 1）: output-format 参照
（共通核とトピック別）。同じ本文を読むツール `docs_read` も、リソースを読めないクライアント向けに置く。

このモジュールを import すると副作用でリソースが mcp に登録される。server は
`from metabolomix.tools import resources as tools_resources` するだけでよい。tools_* / server は import しない。
"""
from mcp.types import ToolAnnotations

from metabolomix.core import mcp_core
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core import session_state
from metabolomix.core.mcp_core import mcp, OUTPUT_FORMAT_DOC


@mcp.resource(
    "lipidmix://docs/output-format",
    name="output_format",
    title="MS-DIAL parser output format and ontology",
    description=(
        "Shared core of the output-format reference (row granularity, lipid-name "
        "grammar, mandatory caveats) plus the index of per-topic sections: dataset, "
        "dcl, eic, identity, mztab, library, curation, pathway, reanalysis, quant. Read before "
        "interpreting dataset results."
    ),
    mime_type="text/markdown",
)
def output_format_reference() -> str:
    """Return the shared ontology core; per-parser topics live in sibling resources."""
    try:
        text = OUTPUT_FORMAT_DOC.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Output format reference was not found: {OUTPUT_FORMAT_DOC}"
        ) from exc
    # リソースが読まれた＝意味論が既に届いているので、以後パーサー出力へ
    # ダイジェストを前置しない（条件付きガードの「未読」条件を解除）。
    session_state.session.output_format_seen = True
    return text


@mcp.resource(
    "lipidmix://docs/output-format/{topic}",
    name="output_format_section",
    title="Per-parser output format section",
    description=(
        "Field-by-field reference for one tool family: dataset (a dataset loaded from "
        "a MS-DIAL analysis folder), dcl (.dcl MS/MS), eic (EIC figures of dataset_plot), "
        "identity, mztab (a dataset loaded from mzTab-M), library (reference-library "
        "MS/MS matching), curation (alignment curation: curation_review's verdicts "
        "and reason codes), pathway (MSLipidMapper pathway projection: node table and "
        "adjacent opposite pairs), reanalysis (dataset handles and the spec-driven "
        "reanalysis tools), quant (internal-standard normalization, polarity "
        "merge, feature selection), or core for the shared ontology. Fetch the topic matching "
        "the output you are about to interpret instead of the whole document."
    ),
    mime_type="text/markdown",
)
def output_format_section(topic: str) -> str:
    """トピック別の出力定義を返す（未知トピックは有効一覧付きで ValueError）。"""
    path = mcp_core.output_format_section_path(topic)
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Output format section was not found: {path}") from exc
    if topic == "core":
        session_state.session.output_format_seen = True
    return text


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def docs_read(topic: str | None = None) -> str:
    """output-format の文書（出力の意味の定義）を Markdown で返す。読み取り専用。

    MCP のリソース `lipidmix://docs/output-format`（topic なし＝共通核）と
    `lipidmix://docs/output-format/{topic}` のツール版。リソースを読めないクライアント
    （Workflow の役など）は、こちらで同じ本文を読む。リソースと同じ文字列を返す。

    topic: 省略すると共通核（行の粒度・脂質名の文法・必須の注意）。有効な topic は
    `core`・`dataset`・`dcl`・`eic`・`identity`・`mztab`・`library`・`curation`・`pathway`・
    `reanalysis`・`quant`。解釈する出力に合う topic を読む。知らない topic は
    `DATASET_BAD_REQUEST` の誤りで、有効な topic の一覧を返す。
    """
    if topic is None:
        return output_format_reference()
    try:
        mcp_core.output_format_section_path(topic)
    except ValueError as exc:
        return mztab_error("DATASET_BAD_REQUEST", str(exc),
                           {"field": "topic", "valid": list(mcp_core.OUTPUT_FORMAT_SECTIONS)})
    return output_format_section(topic)
