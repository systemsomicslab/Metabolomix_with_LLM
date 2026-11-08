"""`MCP_INSTRUCTIONS` と意味論の要約（`SEMANTICS_CAVEAT`）が、いまの公開面だけを語ること。

B-2b で旧 ARF / ARF2 / PAI2 / EIC のツールと `dataset_handle_*` などが消えた。LLM に毎回
届く指示が消えたツールの名前を挙げると、存在しないツールを呼ばせる。ここで縛る。
"""
import asyncio
import re

import server  # noqa: F401  - ツールを登録する副作用
from metabolomix.core import mcp_core, session_state

#: 名前の形でも、ツールではない識別子（引数名・フィールド名・設定のキーなど）。
NON_TOOL_IDENTIFIERS = {
    "assay_kind", "dataset_handle", "all_files", "peak_top", "detected_mask",
    "max_intensity", "rule_hash", "analysis_id", "caveat_mode", "has_msms",
    "feature_matrix", "analysis_dir", "reports_dir", "mztab_path", "job_path",
    "pipeline_path", "pipeline_id", "dataset_root", "file_path", "excluded_samples",
    "needs_input",  # `pipeline_status` の状態名
    "add_subquestions",  # `record_objective` の引数
    "not_found",  # `dcl_find_msms` の status の値
    "missing_state",  # エラー封筒の code
    "matrix_id",  # `dataset_preprocess` が返す行列の ID（`dataset_compare` / `dataset_pca` の引数）
    "retry_tools",  # `needs_input` の封筒の欄
}
WILDCARD_RE = re.compile(r"^(dataset|library|curation|dcl|pathway|console|pipeline)_\*$")
#: バッククォートの中身が「名前」か「名前(...)」か「名前=...」のどれかで、名前に `_` を含むもの。
SPAN_RE = re.compile(r"`([a-z][a-z0-9]*(?:_[a-z0-9]+)+)(?:[(=][^`]*)?`")
WILDCARD_SPAN_RE = re.compile(r"`([a-z]+_\*)`")

REMOVED = (
    "dataset_load", "dataset_handle_open", "dataset_handle_list", "dataset_handle_close",
    "arf_parser", "arf_preprocess", "arf_differential", "arf_exclude", "arf2_parser",
    "arf2_annotate_identities", "pai2_parser", "pai2_inspect_peak", "dcl_parser",
    "eic_parser", "eic_search", "sample_search", "plot_pca_loadings",
    "verify_peak_annotation", "knowledge_coverage", "paper_search", "log_search",
    "update_objective", "list_reports",
    "dataset_differential", "dataset_statistic", "dataset_build_matrix",
)


def _registered_tools() -> set[str]:
    return {tool.name for tool in asyncio.run(mcp_core.mcp.list_tools())}


def _texts():
    return {"MCP_INSTRUCTIONS": mcp_core.MCP_INSTRUCTIONS,
            "SEMANTICS_CAVEAT": session_state.SEMANTICS_CAVEAT}


def test_instructions_name_only_registered_tools():
    registered = _registered_tools()
    for label, text in _texts().items():
        for name in SPAN_RE.findall(text):
            if name in NON_TOOL_IDENTIFIERS:
                continue
            assert name in registered, f"{label} が登録されていないツールを挙げている: {name}"
        for wildcard in WILDCARD_SPAN_RE.findall(text):
            assert WILDCARD_RE.match(wildcard), f"{label} の想定外のワイルドカード: {wildcard}"


def test_instructions_never_mention_removed_tools():
    for label, text in _texts().items():
        for name in REMOVED:
            assert not re.search(rf"\b{name}\b", text), f"{label} に消したツールの名前: {name}"


def test_entry_point_is_list_data_files_then_load_dataset():
    text = mcp_core.MCP_INSTRUCTIONS
    assert "ENTRY POINT" in text
    section = text[text.index("ENTRY POINT"):]
    assert section.index("list_data_files") < section.index("load_dataset")
    # MS-DIAL を回した形跡があれば下流へ、無ければ許可を取ってから Console（変えない規則）。
    assert "ASK THE USER FOR PERMISSION" in section


def test_assay_kind_is_settled_after_load_dataset_with_record_objective():
    text = mcp_core.MCP_INSTRUCTIONS
    section = text[text.index("ASSAY KIND"):]
    assert "record_objective(assay_kind=" in section
    assert "AFTER" in section and "load_dataset" in section


def test_records_need_a_loaded_dataset():
    text = mcp_core.MCP_INSTRUCTIONS
    section = text[text.index("RECORDS"):text.index("MS/MS EVIDENCE")]
    for name in ("record_objective", "write_report", "save_figure", "[paths]", "reports_dir", "loaded dataset"):
        assert name in section


def test_removed_sections_are_gone():
    text = mcp_core.MCP_INSTRUCTIONS
    for heading in ("GATEWAY", "CONFLICTS", "LITERATURE DISCOVERY"):
        assert heading not in text
    for word in ("knowledge/", "playbook"):
        assert word not in text


def test_parallel_reanalysis_names_the_surviving_handle_tools():
    text = mcp_core.MCP_INSTRUCTIONS
    section = text[text.index("PARALLEL REANALYSIS"):]
    for name in ("load_dataset", "dataset_list", "dataset_close", "dataset_handle"):
        assert name in section


def test_msms_evidence_is_written_in_terms_of_dcl_and_surviving_tools():
    text = mcp_core.MCP_INSTRUCTIONS
    section = text[text.index("MS/MS EVIDENCE"):]
    for name in ("dcl_*", "library_msms_qualify", "dataset_claim_evidence"):
        assert name in section
    assert ".dcl" in section


def test_semantics_caveat_describes_the_dataset_state_model():
    caveat = session_state.SEMANTICS_CAVEAT
    assert "特徴×試料" in caveat          # 行の粒度
    assert "detected_mask" in caveat      # gap-fill の印
    assert ".dcl" in caveat               # MS/MS の実スペクトルの所在
    assert "peak_top" in caveat           # EIC の頂点は座標
    # 旧パーサの粒度（ARF 行 / ARF2 行 / PAI2）は語らない。
    for stale in ("ARF行", "ARF2行", "PAI2", "IsGapFilled"):
        assert stale not in caveat
    assert "lipidmix://docs/output-format" in caveat
