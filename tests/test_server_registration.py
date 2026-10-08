"""特性化テスト（安全網）: server の MCP 登録面を固定する。

server.py の分割リファクタ中、ツール/リソースの登録漏れや名前変化を即検出するための
スナップショット。挙動は変えない前提なので、ここが赤くなったら「外形が壊れた」合図。
`import server` が例外を出さないこと自体も回帰対象（循環 import 等の早期検出）。
"""
import asyncio

import server


# MCP 登録ツールの正準スナップショット（sorted）。
# EIC ツールは eicaef_* → eic_* に改称、eicaef_top_peak_tops は強度基準の
# eic_rank_by_max_intensity に置換。PAI2 の PCA 依存 2 ツール
# （pai2_get_top_metabolites / pai2_update_analysis_filter）は撤去。
# arf_re_pca は arf_parser（min_intensity/annotation_keyword を吸収）へ統合し撤去。
# pai2_inspect_metabolite_details は peak 語彙へ統一し pai2_inspect_peak に改称。
# Task18: 生データフォルダ起点の pipeline_plan/run/status/resume/cancel を追加(56→61)。
# v2メタボロミクス Task13: dataset_statistic を追加(61→62)。
# 対話経路S1: pipeline を経ずに解析行列を作る dataset_build_matrix を追加(62→63)。
# Task10: 参照ライブラリの MS/MS 照合 3 本を追加(63→66)。
# 配布更新: 通知だけでなく適用まで行う server_update を追加(66→67)。
# アラインメントのキュレーション: curation_review/submit/flags/view_data を追加(67→71)。
# キュレーションの候補付け: curation_suggest を追加(71→72)。
# 群別強度プロット: arf_plot_group_intensity を追加(72→73)。
# 群別強度プロットの保存: save_group_intensity_figure を追加(73→74)。
# 図の保存の一本化: save_figure に統合し旧 4 ツールを撤去(74→71)。
# 分子種ごとの図: arf_plot_species を追加(71→72)。
EXPECTED_TOOLS = sorted([
    "library_load",
    "library_match_feature",
    "library_plot_mirror",
    "curation_review",
    "curation_submit",
    "curation_flags",
    "curation_view_data",
    "curation_suggest",
    "arf2_annotate_identities",
    "arf2_parser",
    "arf_differential",
    "arf_exclude",
    "arf_export_differential",
    "arf_list_classes",
    "arf_list_sample_roles",
    "arf_list_tags",
    "arf_parser",
    "arf_pca_preprocessed",
    "arf_plot_volcano",
    "arf_plot_group_intensity",
    "arf_plot_species",
    "arf_preprocess",
    "eic_parser",
    "eic_plot_chromatograms",
    "eic_plot_compounds",
    "eic_rank_by_max_intensity",
    "eic_search_by_mz_range",
    "eic_search_by_rt_range",
    "ingest_promote",
    "ingest_reject",
    "ingest_review_queue",
    "ingest_stage",
    "knowledge_coverage",
    "list_data_files",
    "list_reports",
    "load_dataset",
    "log_search",
    "pai2_inspect_peak",
    "pipeline_cancel",
    "pipeline_plan",
    "pipeline_resume",
    "pipeline_run",
    "pipeline_status",
    "console_cleanup",
    "console_method_candidates",
    "console_method_template",
    "console_plan",
    "console_prepare_input",
    "console_run",
    "console_status",
    "dataset_load",
    "dataset_status",
    "job_list",
    "dataset_differential",
    "dataset_export_differential",
    "dataset_pca",
    "dataset_preprocess",
    "dataset_set_sample_metadata",
    "dataset_build_matrix",
    "dataset_statistic",
    "dcl_find_msms",
    "dcl_parser",
    "pai2_parser",
    "paper_search",
    "read_report",
    "record_objective",
    "sample_search",
    "server_update",
    "save_figure",
    "update_objective",
    "verify_peak_annotation",
    "write_report",
])

EXPECTED_RESOURCES = sorted([
    "lipidmix://docs/output-format",
    "lipidmix://knowledge/index",
    "lipidmix://playbook/index",
    "lipidmix://knowledge/inbox",
    "ui://ms-data-parser/curation-viewer",
])

EXPECTED_TEMPLATES = sorted([
    "lipidmix://docs/output-format/{topic}",
    "lipidmix://knowledge/expand/{slug}",
    "lipidmix://playbook/expand/{slug}",
])


def test_tool_count_is_stable():
    # ツール数の正準は EXPECTED_TOOLS。ここに数値リテラルを置くと、
    # ツール追加時に EXPECTED_TOOLS だけ更新されて数値が取り残される
    # （f981cf3 で実際に起きた）。
    tools = asyncio.run(server.mcp.list_tools())
    assert len(tools) == len(EXPECTED_TOOLS)


def test_tool_names_snapshot():
    tools = asyncio.run(server.mcp.list_tools())
    assert sorted(t.name for t in tools) == EXPECTED_TOOLS


def test_static_resources_snapshot():
    resources = asyncio.run(server.mcp.list_resources())
    assert sorted(str(r.uri) for r in resources) == EXPECTED_RESOURCES


def test_resource_templates_snapshot():
    templates = asyncio.run(server.mcp.list_resource_templates())
    assert sorted(t.uriTemplate for t in templates) == EXPECTED_TEMPLATES


def test_import_server_does_not_raise():
    # import server 自体が副作用で例外を出さないことの回帰（循環 import の早期検出）。
    import importlib

    importlib.reload(server)
