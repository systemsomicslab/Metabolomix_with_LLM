"""MCP リソース（@mcp.resource ×7。静的 4・テンプレート 3）: output-format 参照
（共通核とトピック別）、knowledge/playbook の index・expand、knowledge inbox。

このモジュールを import すると副作用でリソースが mcp に登録される。server は
`from lipidmix.tools import resources as tools_resources` するだけでよい。tools_* / server は import しない。
"""
from lipidmix.corpus import knowledge_store
from lipidmix.core import mcp_core
from lipidmix.core import session_state
from lipidmix.core.mcp_core import mcp, OUTPUT_FORMAT_DOC, PLAYBOOK_DIR


@mcp.resource(
    "lipidmix://docs/output-format",
    name="output_format",
    title="MS-DIAL parser output format and ontology",
    description=(
        "Shared core of the output-format reference (row granularity, lipid-name "
        "grammar, mandatory caveats) plus the index of per-topic sections: arf, "
        "arf2, pai2, dcl, eic, identity, mztab, library, curation. Read before "
        "interpreting parser or dataset results."
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
        "Field-by-field reference for one parser/tool family: arf, arf2, pai2, dcl, "
        "eic, identity, mztab (dataset_* / mzTab-M path), library (reference-library "
        "MS/MS matching), curation (alignment curation: curation_review's verdicts "
        "and reason codes), or core for the shared ontology. Fetch the topic matching "
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


# --- 知識・ワークフロー蓄積層（knowledge / playbook） ---
# 索引は frontmatter から動的生成（実INDEXファイルは持たない）。展開は [[link]]
# グラフを構造予算内（max 1 hop / 5本体 / 約15kトークン）で束ねて返す。
# 関連性の判断（どのノートを採用するか）は LLM 側に委ねる。
@mcp.resource(
    "lipidmix://knowledge/index",
    name="knowledge_index",
    title="Knowledge note index (literature-derived)",
    description=(
        "One line per knowledge note (description + claim_strength). Consult this "
        "before interpreting; expand only relevant slugs."
    ),
    mime_type="text/markdown",
)
def knowledge_index() -> str:
    """論文由来の宣言的知識ノートの1行索引を返す。"""
    return knowledge_store.build_index(mcp_core.KNOWLEDGE_DIR, "knowledge")


@mcp.resource(
    "lipidmix://playbook/index",
    name="playbook_index",
    title="Playbook index (reusable analysis workflows)",
    description=(
        "One line per playbook note (when_to_use). Consult this before proposing a "
        "workflow; expand only relevant slugs."
    ),
    mime_type="text/markdown",
)
def playbook_index() -> str:
    """再利用可能な解析手順ノートの1行索引を返す。"""
    return knowledge_store.build_index(PLAYBOOK_DIR, "playbook")


@mcp.resource(
    "lipidmix://knowledge/expand/{slug}",
    name="knowledge_expand",
    title="Expand a knowledge note with its 1-hop neighbors",
    description=(
        "Returns the note body plus directly-linked neighbors within a structural "
        "budget (1 hop, 5 bodies, ~15k tokens). Overflow is demoted to index lines."
    ),
    mime_type="text/markdown",
)
def knowledge_expand(slug: str) -> str:
    """knowledge ノートを1ホップ展開して予算内で返す。"""
    return knowledge_store.expand(slug, [mcp_core.KNOWLEDGE_DIR, PLAYBOOK_DIR])


@mcp.resource(
    "lipidmix://playbook/expand/{slug}",
    name="playbook_expand",
    title="Expand a playbook note with its 1-hop neighbors",
    description=(
        "Returns the playbook body plus directly-linked neighbors within a "
        "structural budget (1 hop, 5 bodies, ~15k tokens). Overflow is demoted."
    ),
    mime_type="text/markdown",
)
def playbook_expand(slug: str) -> str:
    """playbook ノートを1ホップ展開して予算内で返す。"""
    return knowledge_store.expand(slug, [PLAYBOOK_DIR, mcp_core.KNOWLEDGE_DIR])


@mcp.resource(
    "lipidmix://knowledge/inbox",
    name="knowledge_inbox",
    title="Pending literature notes awaiting review",
    description=(
        "Quarantined (speculative) notes from gap-driven discovery, grouped by "
        "analysis_id/Qi. Promote with ingest_promote or discard with ingest_reject."
    ),
    mime_type="text/markdown",
)
def knowledge_inbox() -> str:
    """_inbox の保留中ノートを found_for/query/score 付きで一覧する。"""
    return knowledge_store.build_inbox_index(mcp_core.KNOWLEDGE_DIR)
