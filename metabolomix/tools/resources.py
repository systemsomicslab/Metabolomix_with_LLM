"""MCP リソース（@mcp.resource ×2。静的 1・テンプレート 1）: output-format 参照
（共通核とトピック別）。

このモジュールを import すると副作用でリソースが mcp に登録される。server は
`from metabolomix.tools import resources as tools_resources` するだけでよい。tools_* / server は import しない。
"""
from metabolomix.core import mcp_core
from metabolomix.core import session_state
from metabolomix.core.mcp_core import mcp, OUTPUT_FORMAT_DOC


@mcp.resource(
    "lipidmix://docs/output-format",
    name="output_format",
    title="MS-DIAL parser output format and ontology",
    description=(
        "Shared core of the output-format reference (row granularity, lipid-name "
        "grammar, mandatory caveats) plus the index of per-topic sections: dataset, "
        "dcl, eic, identity, mztab, library, curation, pathway, reanalysis. Read before "
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
    session_state.session.sections_seen.add("core")
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
        "reanalysis tools), or core for the shared ontology. Fetch the topic matching "
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
    session_state.session.sections_seen.add(topic)
    if topic == "core":
        session_state.session.output_format_seen = True
    return text
