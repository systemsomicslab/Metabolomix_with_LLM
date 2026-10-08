"""docs/workflow/ の呼び出し連鎖が実装と一致していることを検証する。

文書には行番号を書かない代わりに、参照された関数が本当に存在することを
AST で確かめる。実装をリファクタして文書が取り残されたらここが落ちる。
"""
import ast
import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_DIR = REPO_ROOT / "docs" / "workflow"

# 呼び出し連鎖行の書式:  `1. └─ metabolomix/arf/reader.py  run_pca()`
# 分岐（排他的な if/else 経路）は `├─` / `│  └─` で表し、経路の説明を角括弧で
# 添えてよい: `4. ├─ [specs 省略時] metabolomix/tools/samples.py  _apply_role_filter()`
CHAIN_RE = re.compile(
    r"^\s*\d+\.\s*(?:[└├│─]\s*)*(?:\[[^\]]*\]\s*)?"
    r"(?P<path>[\w/]+\.py)\s+(?P<func>[\w.]+)\(\)"
)
# 節見出し: `## arf_parser`
HEADING_RE = re.compile(r"^##\s+(?P<tool>[a-z][a-z0-9_]*)\s*$", re.MULTILINE)

# 設計書 §4 の対象範囲表と一対一で対応する。増減させる場合は設計書も直すこと。
IN_SCOPE: dict[str, tuple[str, ...]] = {
    "dataset.md": ("list_data_files", "load_dataset", "sample_search"),
    "arf.md": (
        "arf_parser", "arf_list_classes", "arf_list_tags", "arf_list_sample_roles",
        "arf_exclude", "arf_preprocess", "arf_pca_preprocessed", "arf_differential",
        "arf_export_differential", "arf_plot_volcano", "arf_plot_group_intensity", "arf_plot_species",
        "arf_pca_species",
    ),
    "arf2.md": ("arf2_parser", "arf2_annotate_identities"),
    "pai2.md": ("pai2_parser", "pai2_inspect_peak", "verify_peak_annotation"),
    "dcl.md": ("dcl_parser", "dcl_find_msms"),
    "eic.md": (
        "eic_parser", "eic_search_by_mz_range", "eic_search_by_rt_range",
        "eic_rank_by_max_intensity", "eic_plot_chromatograms", "eic_plot_compounds",
    ),
    "plots.md": ("save_figure",),
    "mztab.md": ("dataset_load", "dataset_status"),
    "dataset_analysis.md": (
        "dataset_preprocess", "dataset_pca", "dataset_differential",
        "dataset_export_differential", "dataset_set_sample_metadata",
    ),
    "pipeline.md": (
        "pipeline_plan", "pipeline_run", "pipeline_status",
        "pipeline_resume", "pipeline_cancel",
    ),
    "metabolomics.md": ("dataset_statistic", "dataset_build_matrix"),
    "library.md": ("library_load", "library_match_feature", "library_plot_mirror"),
    "curation.md": (
        "curation_review", "curation_suggest", "curation_submit", "curation_flags",
        "curation_view_data",
    ),
}

# 今回の範囲外。文書に混入したら落とす（線引きを固定するため）。
OUT_OF_SCOPE: tuple[str, ...] = (
    "record_objective", "update_objective", "log_search", "knowledge_coverage",
    "paper_search", "ingest_stage", "ingest_review_queue", "ingest_promote",
    "ingest_reject", "write_report", "read_report", "list_reports",
    # Console 実行層（Phase 2 追加・ワークフロー文書は Phase 2 完了後に追加）
    "console_plan", "console_prepare_input", "console_method_template",
    "console_method_candidates",
    "console_run", "console_status", "console_cleanup", "job_list",
    # サーバ自身の保守（解析の流れには現れない）
    "server_update",
)


def _defined_names(py_path: Path) -> set[str]:
    """モジュール内で参照可能な名前を集める。

    トップレベルの def / async def / class に加え、`ClassName.method` 形式と、
    `list_data_files = mcp.tool(...)(path_resolvers.list_data_files)` のような
    モジュールレベル代入も拾う。
    """
    tree = ast.parse(py_path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.ClassDef):
            names.add(node.name)
            for sub in node.body:
                if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    names.add(f"{node.name}.{sub.name}")
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    names.add(target.id)
    return names


def _iter_chain_refs():
    """全文書の呼び出し連鎖行を (文書名, 行番号, パス, 関数名) で返す。"""
    for md in sorted(WORKFLOW_DIR.rglob("*.md")):
        for lineno, line in enumerate(md.read_text(encoding="utf-8").splitlines(), 1):
            m = CHAIN_RE.match(line)
            if m:
                yield md.name, lineno, m.group("path"), m.group("func")


class TestWorkflowDocsExist(unittest.TestCase):
    def test_all_expected_documents_exist(self):
        self.assertTrue(WORKFLOW_DIR.is_dir(), f"{WORKFLOW_DIR} がない")
        expected = {"index.md", *IN_SCOPE}
        actual = {p.name for p in WORKFLOW_DIR.rglob("*.md")}
        self.assertEqual(expected, actual)


class TestChainReferencesResolve(unittest.TestCase):
    def test_every_reference_points_at_an_existing_file(self):
        for doc, lineno, path, _func in _iter_chain_refs():
            with self.subTest(doc=doc, lineno=lineno, path=path):
                self.assertTrue(
                    (REPO_ROOT / path).is_file(),
                    f"{doc}:{lineno} が存在しないファイルを参照: {path}",
                )

    def test_every_reference_points_at_a_defined_name(self):
        cache: dict[str, set[str]] = {}
        for doc, lineno, path, func in _iter_chain_refs():
            with self.subTest(doc=doc, lineno=lineno, ref=f"{path}:{func}"):
                target = REPO_ROOT / path
                if not target.is_file():
                    self.skipTest("ファイル不在は別テストで報告済み")
                if path not in cache:
                    cache[path] = _defined_names(target)
                self.assertIn(
                    func,
                    cache[path],
                    f"{doc}:{lineno} の {func}() が {path} に定義されていない",
                )

    def test_documents_contain_at_least_one_chain(self):
        counts: dict[str, int] = {name: 0 for name in IN_SCOPE}
        for doc, _lineno, _path, _func in _iter_chain_refs():
            if doc in counts:
                counts[doc] += 1
        for doc, n in counts.items():
            if not (WORKFLOW_DIR / doc).is_file():
                continue  # 不在は test_all_expected_documents_exist が報告する
            with self.subTest(doc=doc):
                self.assertGreater(n, 0, f"{doc} に呼び出し連鎖が 1 行もない")


class TestScopeBoundary(unittest.TestCase):
    def test_every_in_scope_tool_has_a_section(self):
        for doc, tools in IN_SCOPE.items():
            path = WORKFLOW_DIR / doc
            if not path.is_file():
                continue  # 不在は test_all_expected_documents_exist が報告する
            headings = set(HEADING_RE.findall(path.read_text(encoding="utf-8")))
            for tool in tools:
                with self.subTest(doc=doc, tool=tool):
                    self.assertIn(tool, headings, f"{doc} に `## {tool}` 節がない")

    def test_out_of_scope_tools_have_no_section(self):
        for md in sorted(WORKFLOW_DIR.rglob("*.md")):
            headings = set(HEADING_RE.findall(md.read_text(encoding="utf-8")))
            for tool in OUT_OF_SCOPE:
                with self.subTest(doc=md.name, tool=tool):
                    self.assertNotIn(
                        tool, headings,
                        f"{md.name} に範囲外ツール `## {tool}` の節がある",
                    )

    def test_scope_totals_match_registered_tool_count(self):
        """対象範囲の分割が、実際に登録されているツール数と一致することを縛る。

        46 / 20 という分割は設計上の線引きだが、その合計は「全ツールを漏れなく
        分類した」という主張でもある。リテラルどうしの比較では自分自身を検証して
        しまうので、実際の登録数を引いて突き合わせる。
        """
        import asyncio

        import server

        registered = {t.name for t in asyncio.run(server.mcp.list_tools())}
        classified = {t for tools in IN_SCOPE.values() for t in tools} | set(OUT_OF_SCOPE)
        self.assertEqual(
            classified,
            registered,
            "対象範囲の分類と登録済みツールが一致しない",
        )
        self.assertEqual(sum(len(v) for v in IN_SCOPE.values()), 52)
