"""record_objective ツール（作成と更新を 1 本で担う）の検証。

純ロジック（frontmatter・小問の採番）は tests/test_notes.py。ここは server ツールの振る舞い。
"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from metabolomix.tools import notes
import server
from metabolomix.core import mcp_core
from metabolomix.core import session_state


class ObjectiveServerToolTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self.dataset_dir = base / "dataset"
        self.dataset_dir.mkdir()
        self.reports = self.dataset_dir / "reports"
        self._env = mock.patch.dict("os.environ")
        self._env.start()
        session_state.session = session_state.AnalysisSession()
        # 記録の置き場所は直近のデータセットの解析フォルダ配下 reports/。
        from metabolomix.mztab.dataset_state import DatasetState
        ds = DatasetState()
        ds.analysis_dir = str(self.dataset_dir)
        session_state.session.dataset = ds

    def tearDown(self):
        self._env.stop()
        self._tmp.cleanup()

    def _create(self, analysis_id, **over):
        args = dict(dataset="NEG", polarity="NEG", groups=["a", "b"], comparison="c",
                    sub_questions=["最初の問い"])
        args.update(over)
        return server.record_objective(analysis_id, **args)

    def test_create_writes_the_objective_and_does_not_point_at_removed_tools(self):
        out = self._create("exp-2", biological_context="RAW macrophage",
                           sub_questions=["糖代謝のグルコース取り込み", "エーテル脂質は減るか"])
        meta, subqs, _body = notes.parse_objective(self.reports / "exp-2.objective.md")
        self.assertEqual(meta["biological_context"], "RAW macrophage")
        self.assertEqual([label for label, _ in subqs], ["Q1", "Q2"])
        self.assertIn("小問2件", out)
        self.assertNotIn("knowledge_coverage", out)

    def test_create_without_required_fields_names_them(self):
        out = server.record_objective("exp-new", confirmed_objective="確定")
        self.assertIn("見つかりません", out)
        for name in ("dataset", "polarity", "groups", "comparison", "sub_questions"):
            self.assertIn(name, out)
        self.assertFalse((self.reports / "exp-new.objective.md").exists())

    def test_update_changes_only_the_passed_fields(self):
        self._create("exp-3", biological_context="RAW macrophage")
        out = server.record_objective(
            "exp-3", confirmed_objective="確定", status="done", add_subquestions=["創発の問い"])
        self.assertIn("更新", out)
        meta, subqs, _ = notes.parse_objective(self.reports / "exp-3.objective.md")
        self.assertEqual(meta["confirmed_objective"], "確定")
        self.assertEqual(meta["status"], "done")
        self.assertIs(meta["confirmed"], True)
        self.assertEqual(meta["biological_context"], "RAW macrophage")  # 渡さなかった欄は温存
        self.assertEqual(meta["dataset"], "NEG")
        self.assertEqual([label for label, _ in subqs], ["Q1", "Q2"])
        self.assertIn("創発の問い", dict(subqs)["Q2"])

    def test_update_preserves_multiline_biological_context(self):
        context = "RAW264.7 細胞\nLPS 100 ng/mL, 24 h\n対照は無処理"
        self._create("exp-ml", biological_context=context)
        server.record_objective("exp-ml", status="done", confirmed_objective="確定\n二行目")
        meta, _subqs, _ = notes.parse_objective(self.reports / "exp-ml.objective.md")
        self.assertEqual(meta["biological_context"], context)
        self.assertEqual(meta["confirmed_objective"], "確定\n二行目")
        self.assertEqual(meta["status"], "done")

    def test_update_does_not_rewrite_existing_subquestions(self):
        self._create("exp-3b")
        out = server.record_objective("exp-3b", sub_questions=["別の問い"])
        _meta, subqs, _ = notes.parse_objective(self.reports / "exp-3b.objective.md")
        self.assertEqual([text for _label, text in subqs], ["最初の問い"])
        self.assertIn("add_subquestions", out)

    def test_record_objective_persists_and_applies_assay_kind(self):
        self._create("exp-4", assay_kind="metabolite")
        # セッションの解釈規則が切り替わる（脂質名文法を当てない）
        self.assertEqual(session_state.session.assay_kind, "metabolite")
        meta, _subqs, _body = notes.parse_objective(self.reports / "exp-4.objective.md")
        self.assertEqual(meta["assay_kind"], "metabolite")

    def test_record_objective_defaults_to_unknown_kind(self):
        self._create("exp-5")
        self.assertEqual(session_state.session.assay_kind, "unknown")
        meta, _subqs, _body = notes.parse_objective(self.reports / "exp-5.objective.md")
        self.assertEqual(meta["assay_kind"], "unknown")

    def test_record_objective_rejects_bad_assay_kind(self):
        out = self._create("exp-6", assay_kind="proteomics")
        self.assertIn("assay_kind", out)
        self.assertFalse((self.reports / "exp-6.objective.md").exists())

    def test_update_switches_assay_kind_only_when_passed(self):
        self._create("exp-7", assay_kind="metabolite")
        server.record_objective("exp-7", status="active")
        self.assertEqual(session_state.session.assay_kind, "metabolite")  # 渡さなければ触らない
        server.record_objective("exp-7", assay_kind="lipid")
        self.assertEqual(session_state.session.assay_kind, "lipid")
        meta, _subqs, _body = notes.parse_objective(self.reports / "exp-7.objective.md")
        self.assertEqual(meta["assay_kind"], "lipid")

    def test_objective_lives_with_the_data_not_the_project(self):
        """objective は解析フォルダ配下 reports/ に書き、リポジトリ側に analyses/ を作らない。"""
        self._create("exp-8")
        self.assertTrue((self.reports / "exp-8.objective.md").is_file())
        self.assertFalse((mcp_core.BASE_DIR / "analyses" / "exp-8.md").exists())
        self.assertFalse((mcp_core.BASE_DIR / "analyses" / "exp-8.objective.md").exists())

    def test_report_with_same_id_does_not_clobber_objective(self):
        """同じ analysis_id のレポートと objective は同じ reports/ に並んで共存する。"""
        self._create("exp-9")
        server.write_report("exp-9", "NEG", "## 結論\nx")
        self.assertTrue((self.reports / "exp-9.md").is_file())
        meta, subqs, _body = notes.parse_objective(self.reports / "exp-9.objective.md")
        self.assertEqual(meta["type"], "objective")
        self.assertEqual([label for label, _ in subqs], ["Q1"])
        listing = server.read_report()
        self.assertEqual(listing.count("exp-9"), 1)


if __name__ == "__main__":
    unittest.main()
