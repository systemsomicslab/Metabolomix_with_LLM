"""objective ライフサイクル（knowledge_store 関数 ＋ server ツール）の検証。"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from metabolomix.corpus import knowledge_store as ks
import server
from metabolomix.core import mcp_core
from metabolomix.corpus import paper_ingest


class ObjectiveStoreTests(unittest.TestCase):
    def _write(self, directory: Path) -> Path:
        return ks.write_objective(
            directory,
            "exp-1",
            {
                "dataset": "NEG / x",
                "polarity": "NEG",
                "groups": ["control", "treatment"],
                "comparison": "群間差",
                "biological_context": "",
                "inferred_objective": "推測",
                "confirmed_objective": "",
                "expected_biology": [],
            },
            ["どのクラスに差が出るか", "エーテル脂質は減るか"],
        )

    def test_write_and_parse_labeled_subquestions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            self.assertTrue(path.is_file())
            meta, subqs, _body = ks.parse_objective(path)
            self.assertEqual(meta["analysis_id"], "exp-1")
            self.assertEqual([label for label, _ in subqs], ["Q1", "Q2"])
            self.assertIn("どのクラス", subqs[0][1])
            self.assertEqual(str(meta["confirmed"]), "False")

    def test_update_meta_preserves_body(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            ks.update_objective_meta(path, {"confirmed_objective": "確定した目的", "biological_context": "RAW macrophage"})
            meta, subqs, _ = ks.parse_objective(path)
            self.assertEqual(meta["confirmed_objective"], "確定した目的")
            self.assertEqual(meta["biological_context"], "RAW macrophage")
            self.assertEqual(str(meta["confirmed"]), "True")
            self.assertEqual(len(subqs), 2)  # 本文の小問は温存

    def test_add_subquestions_numbers_after_existing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            ks.add_subquestions(path, ["創発した問い"])
            _meta, subqs, _ = ks.parse_objective(path)
            self.assertEqual([label for label, _ in subqs], ["Q1", "Q2", "Q3"])
            self.assertIn("創発した問い", dict(subqs)["Q3"])

    def test_search_log_roundtrip_with_commas(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp))
            ks.append_search_log(path, "Q2", "2026-06-14", "plasmalogen, oxidation, macrophage", 8, 0)
            searched = ks.searched_labels(path)
            self.assertIn("Q2", searched)
            self.assertIn("hits=8", searched["Q2"])
            # 小問パースは探索ログ行に汚染されない
            _meta, subqs, _ = ks.parse_objective(path)
            self.assertEqual([label for label, _ in subqs], ["Q1", "Q2"])


class ObjectiveServerToolTests(unittest.TestCase):
    def setUp(self):
        self._orig_data = mcp_core.DATA_DIR
        self._orig_knowledge = mcp_core.KNOWLEDGE_DIR
        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        mcp_core.DATA_DIR = base / "dataset"
        mcp_core.KNOWLEDGE_DIR = base / "knowledge"
        mcp_core.DATA_DIR.mkdir()
        mcp_core.KNOWLEDGE_DIR.mkdir()
        self.reports = mcp_core.DATA_DIR / "reports"
        self._env = mock.patch.dict("os.environ")
        self._env.start()
        os.environ.pop("LIPIDMIX_REPORTS_DIR", None)

    def tearDown(self):
        self._env.stop()
        mcp_core.DATA_DIR = self._orig_data
        mcp_core.KNOWLEDGE_DIR = self._orig_knowledge
        self._tmp.cleanup()

    def test_record_then_coverage_then_log_churn(self):
        server.record_objective(
            "exp-2", "NEG / x", "NEG", ["control", "treatment"], "群間差",
            ["糖代謝のグルコース取り込み", "エーテル脂質は減るか"],
            biological_context="RAW macrophage",
        )
        cov = server.knowledge_coverage("exp-2")
        self.assertIn("Q1", cov)
        self.assertIn("Q2", cov)
        self.assertIn("GAP", cov)  # knowledge 空なので GAP

        # 探索を記録すると coverage に already searched 注記が出る
        server.log_search("exp-2", "Q1", "glucose uptake", 5, 0)
        cov2 = server.knowledge_coverage("exp-2")
        self.assertIn("already searched", cov2)

    def test_update_adds_subquestion(self):
        server.record_objective("exp-3", "NEG", "NEG", ["a", "b"], "c", ["最初の問い"])
        server.update_objective("exp-3", add_subquestions=["創発の問い"], confirmed_objective="確定")
        cov = server.knowledge_coverage("exp-3")
        self.assertIn("Q2", cov)
        self.assertIn("創発の問い", cov)


    def test_record_objective_persists_and_applies_assay_kind(self):
        from metabolomix.core import session_state
        session_state.session = session_state.AnalysisSession()
        server.record_objective(
            "exp-4", "HILIC / POS", "POS", ["control", "treatment"], "群間差",
            ["どの経路が動くか"], assay_kind="metabolite",
        )
        # セッションの解釈規則が切り替わる（脂質名文法を当てない）
        self.assertEqual(session_state.session.assay_kind, "metabolite")
        meta, _subqs, _body = ks.parse_objective(self.reports / "exp-4.objective.md")
        self.assertEqual(meta["assay_kind"], "metabolite")

    def test_record_objective_defaults_to_unknown_kind(self):
        from metabolomix.core import session_state
        session_state.session = session_state.AnalysisSession()
        server.record_objective("exp-5", "NEG", "NEG", ["a", "b"], "c", ["問い"])
        self.assertEqual(session_state.session.assay_kind, "unknown")
        meta, _subqs, _body = ks.parse_objective(self.reports / "exp-5.objective.md")
        self.assertEqual(meta["assay_kind"], "unknown")

    def test_record_objective_rejects_bad_assay_kind(self):
        out = server.record_objective(
            "exp-6", "NEG", "NEG", ["a", "b"], "c", ["問い"], assay_kind="proteomics")
        self.assertIn("assay_kind", out)
        self.assertFalse((self.reports / "exp-6.objective.md").exists())

    def test_update_objective_switches_assay_kind(self):
        from metabolomix.core import session_state
        session_state.session = session_state.AnalysisSession()
        server.record_objective("exp-7", "NEG", "NEG", ["a", "b"], "c", ["問い"])
        server.update_objective("exp-7", assay_kind="lipid")
        self.assertEqual(session_state.session.assay_kind, "lipid")
        meta, _subqs, _body = ks.parse_objective(self.reports / "exp-7.objective.md")
        self.assertEqual(meta["assay_kind"], "lipid")

    def test_objective_lives_with_the_data_not_the_project(self):
        """objective は解析フォルダ配下 reports/ に書き、リポジトリ側に analyses/ を作らない。"""
        server.record_objective("exp-8", "NEG", "NEG", ["a", "b"], "c", ["問い"])
        self.assertTrue((self.reports / "exp-8.objective.md").is_file())
        self.assertFalse((mcp_core.BASE_DIR / "analyses" / "exp-8.md").exists())
        self.assertFalse((mcp_core.BASE_DIR / "analyses" / "exp-8.objective.md").exists())

    def test_report_with_same_id_does_not_clobber_objective(self):
        """同じ analysis_id のレポートと objective は同じ reports/ に並んで共存する。"""
        server.record_objective("exp-9", "NEG", "NEG", ["a", "b"], "c", ["問い"])
        server.write_report("exp-9", "NEG", "## 結論\nx")
        self.assertTrue((self.reports / "exp-9.md").is_file())
        meta, subqs, _body = ks.parse_objective(self.reports / "exp-9.objective.md")
        self.assertEqual(meta["type"], "objective")
        self.assertEqual([label for label, _ in subqs], ["Q1"])
        self.assertIn("Q1", server.knowledge_coverage("exp-9"))
        listing = server.list_reports()
        self.assertEqual(listing.count("exp-9"), 1)

    def test_missing_objective_message(self):
        self.assertIn("見つかりません", server.knowledge_coverage("nope"))
        self.assertIn("見つかりません", server.log_search("nope", "Q1", "q", 0, 0))


class PaperSearchSlimTests(unittest.TestCase):
    def test_paper_search_caps_abstract_and_keeps_titles(self):
        fake = [
            {"title": "T1", "abstract": "L" * 1200, "journal": "J",
             "year": 2024, "doi": "d1", "pmid": "1"},
            {"title": "T2", "abstract": "short abstract", "journal": "J",
             "year": 2024, "doi": "d2", "pmid": "2"},
        ]
        with mock.patch.object(paper_ingest, "search_europepmc", lambda q, n: fake), \
             mock.patch.object(paper_ingest, "check_retraction", lambda c: c), \
             mock.patch.object(paper_ingest, "deduplicate", lambda c, e: c):
            out = server.paper_search("q", max_results=5)
        # 全候補の title は残る（breadth 保持）
        self.assertIn("## T1", out)
        self.assertIn("## T2", out)
        # 長い抄録は上限＋截断マーカー、短い抄録はそのまま
        self.assertIn("…（截断）", out)
        self.assertIn("short abstract", out)
        for line in out.splitlines():
            if line.startswith("- abstract:"):
                self.assertLessEqual(
                    len(line), len("- abstract: ") + 500 + len("…（截断）"))


if __name__ == "__main__":
    unittest.main()
