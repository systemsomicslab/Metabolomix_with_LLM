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


def _split_sections(text: str) -> tuple[str, list[tuple[str, str]]]:
    """Markdown を、最初の `## ` 見出しより前（題と前置き）と `## ` 見出しごとの節に分ける。

    節は次の `## ` 見出しの直前まで（`###` の小節を含む）。コードブロックの中の `## ` は見出しとみなさない。
    """
    head: list[str] = []
    sections: list[tuple[str, list[str]]] = []
    in_code = False
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith("```"):
            in_code = not in_code
        if not in_code and line.startswith("## "):
            sections.append((line[3:].strip(), [line]))
        elif sections:
            sections[-1][1].append(line)
        else:
            head.append(line)
    return "".join(head), [(title, "".join(body)) for title, body in sections]


def _outline(topic: str, text: str) -> str:
    _, sections = _split_sections(text)
    rows = "".join(f"{title}\t{len(body)}\n" for title, body in sections)
    return (f"# {topic} の節（`docs_read(topic=..., section=<見出し>)` で 1 節だけ読める）\n\n"
            f"heading\tchars\n{rows}")


def _pick_section(topic: str, text: str, section: str):
    """見出しの完全一致 → 大文字小文字を無視した部分一致の順で 1 節を選ぶ。0 個か複数なら誤りの封筒。"""
    _, sections = _split_sections(text)
    titles = [title for title, _ in sections]
    exact = [body for title, body in sections if title == section]
    if len(exact) == 1:
        return exact[0]
    needle = section.strip().lower()
    hits = [(title, body) for title, body in sections if needle and needle in title.lower()]
    if len(hits) == 1:
        return hits[0][1]
    reason = "一致する見出しがありません" if not hits else "複数の見出しに一致します"
    return mztab_error(
        "DATASET_BAD_REQUEST", f"topic {topic!r} の section {section!r}: {reason}。",
        {"field": "section", "headings": [t for t, _ in hits] if hits else titles})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def docs_read(topic: str | None = None, section: str | None = None, outline: bool = False) -> str:
    """output-format の文書（出力の意味の定義）を Markdown で返す。読み取り専用。

    MCP のリソース `lipidmix://docs/output-format`（topic なし＝共通核）と
    `lipidmix://docs/output-format/{topic}` のツール版。リソースを読めないクライアント
    （Workflow の役など）は、こちらで同じ本文を読む。`section` と `outline` を省けば、
    リソースと同じ文字列を返す。

    topic: 省略すると共通核（行の粒度・脂質名の文法・必須の注意）。有効な topic は
    `core`・`dataset`・`dcl`・`eic`・`identity`・`mztab`・`library`・`curation`・`pathway`・
    `reanalysis`・`quant`。解釈する出力に合う topic を読む。知らない topic は
    `DATASET_BAD_REQUEST` の誤りで、有効な topic の一覧を返す。
    section: `## ` 見出しを 1 つ指して、その節（小節を含む）だけを返す。見出しの完全一致、無ければ
        大文字小文字を無視した部分一致で選ぶ。0 個か複数に当たれば `DATASET_BAD_REQUEST` で、
        候補の見出しを返す。大きな topic（dataset・reanalysis は約 5 万字）を分けて読むときに使う。
    outline: True なら本文の代わりに `## ` 見出しと字数の TSV を返す（section より先に効く）。
    """
    if topic is None and section is None and not outline:
        return output_format_reference()
    name = "core" if topic is None else topic
    try:
        mcp_core.output_format_section_path(name)
    except ValueError as exc:
        return mztab_error("DATASET_BAD_REQUEST", str(exc),
                           {"field": "topic", "valid": list(mcp_core.OUTPUT_FORMAT_SECTIONS)})
    if section is None and not outline:
        return output_format_section(name)
    text = mcp_core.output_format_section_path(name).read_text(encoding="utf-8")
    if outline:
        return _outline(name, text)
    return _pick_section(name, text, section)
