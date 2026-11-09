import json
import os
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # import server が pyplot を読む前にヘッドレス指定

import server
from metabolomix.core import session_state
from metabolomix.core import mcp_core
from metabolomix.tools import notes
from tests.config_fixtures import config_env


def _open_analysis_dir(directory):
    """解析フォルダが `directory` のデータセットを直近にする（記録の置き場所の基準）。"""
    from metabolomix.mztab.dataset_state import DatasetState

    session_state.session = session_state.AnalysisSession()
    ds = DatasetState()
    ds.analysis_dir = str(directory)
    session_state.session.dataset = ds


class WriteLocationHelpers(unittest.TestCase):
    def test_dir_is_writable_true_for_new_dir(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "reports"
            self.assertTrue(server._dir_is_writable(target))
            self.assertTrue(target.is_dir())

    def test_dir_is_writable_false_when_parent_is_file(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            blocker = Path(d) / "blocker"
            blocker.write_text("x", encoding="utf-8")
            self.assertFalse(server._dir_is_writable(blocker / "reports"))

    def test_first_writable_dir_skips_unwritable(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            blocker = Path(d) / "blocker"
            blocker.write_text("x", encoding="utf-8")
            bad = blocker / "reports"
            good = Path(d) / "ok"
            self.assertEqual(server._first_writable_dir([bad, good]), good)

    def test_first_writable_dir_raises_when_none(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            blocker = Path(d) / "blocker"
            blocker.write_text("x", encoding="utf-8")
            with self.assertRaises(OSError):
                server._first_writable_dir([blocker / "a", blocker / "b"])

    def test_build_report_meta_shape(self):
        meta = server._build_report_meta("a-1", "DS", "draft")
        self.assertEqual(meta["type"], "report")
        self.assertEqual(meta["analysis_id"], "a-1")
        self.assertEqual(meta["dataset"], "DS")
        self.assertEqual(meta["status"], "draft")
        self.assertNotIn("knowledge_refs", meta)
        self.assertRegex(meta["date"], r"^\d{4}-\d{2}-\d{2}$")


class ReportToolTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        _open_analysis_dir(self.tmp)
        self.enterContext(config_env(self.tmp, {"paths.reports_dir": self.tmp / "reports_fallback"}))

    def tearDown(self):
        self._tmp.cleanup()

    def test_write_report_creates_file_with_frontmatter_and_body(self):
        msg = server.write_report("a-1", "NEG / DS", "## 目的\nグループ比較")
        path = self.tmp / "reports_fallback" / "a-1.md"
        self.assertTrue(path.is_file())
        self.assertIn("reports", msg)
        meta, body = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["type"], "report")
        self.assertEqual(meta["analysis_id"], "a-1")
        self.assertEqual(meta["status"], "draft")
        self.assertIn("## 目的", body)

    def test_write_report_overwrites_on_second_call(self):
        server.write_report("a-1", "DS", "## 目的\n古い本文")
        server.write_report("a-1", "DS", "## 目的\n新しい本文", status="final")
        path = self.tmp / "reports_fallback" / "a-1.md"
        meta, body = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["status"], "final")
        self.assertIn("新しい本文", body)
        self.assertNotIn("古い本文", body)

    def test_read_report_round_trips(self):
        server.write_report("a-1", "DS", "## 結論\nXがYより高い")
        text = server.read_report("a-1")
        self.assertIn("## 結論", text)
        self.assertIn("XがYより高い", text)

    def test_read_report_missing_returns_guidance(self):
        text = server.read_report("does-not-exist")
        self.assertIn("見つかりません", text)

    def test_read_report_without_id_shows_written_report(self):
        server.write_report("a-1", "DS", "## 目的\nx", status="final")
        listing = server.read_report()
        self.assertIn("a-1", listing)
        self.assertIn("status=final", listing)

    def test_read_report_without_id_empty(self):
        listing = server.read_report()
        self.assertIn("まだありません", listing)


class PcaPlotHelperTests(unittest.TestCase):
    def test_pca_scatter_arrays_extracts_points(self):
        plot = {
            "title": "T",
            "x_label": "PC1 (50.00%)",
            "y_label": "PC2 (30.00%)",
            "points": [
                {"x": 1.0, "y": 2.0, "label": "s1"},
                {"x": -1.0, "y": 0.5, "label": "s2"},
            ],
        }
        xs, ys, labels, x_label, y_label, title = server._pca_scatter_arrays(plot)
        self.assertEqual(xs, [1.0, -1.0])
        self.assertEqual(ys, [2.0, 0.5])
        self.assertEqual(labels, ["s1", "s2"])
        self.assertEqual(x_label, "PC1 (50.00%)")
        self.assertEqual(title, "T")

class ReportEdgeCaseTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        _open_analysis_dir(self.tmp)
        self.enterContext(config_env(self.tmp, {"paths.reports_dir": self.tmp / "reports_fallback"}))

    def tearDown(self):
        self._tmp.cleanup()

    def test_write_report_refuses_slug_collision_with_different_id(self):
        # Both analysis_ids collapse to the same slug.
        slug = notes.make_slug("Group A vs B")
        self.assertEqual(slug, notes.make_slug("group-a-vs-b"))

        server.write_report("Group A vs B", "DS", "## 目的\nfirst")
        msg = server.write_report("group-a-vs-b", "DS", "## 目的\nsecond")

        self.assertIn("衝突", msg)
        # The original report must be intact (not overwritten).
        path = self.tmp / "reports_fallback" / f"{slug}.md"
        meta, body = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["analysis_id"], "Group A vs B")
        self.assertIn("first", body)
        self.assertNotIn("second", body)

    def test_write_report_same_id_still_overwrites(self):
        # Regression: the collision guard must NOT block same-id overwrites.
        server.write_report("a-1", "DS", "## 目的\n古い")
        msg = server.write_report("a-1", "DS", "## 目的\n新しい", status="final")
        self.assertIn("保存", msg)
        path = self.tmp / "reports_fallback" / "a-1.md"
        meta, body = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["status"], "final")
        self.assertIn("新しい", body)
        self.assertNotIn("古い", body)

    def _write_report_file(self, directory, slug, analysis_id, marker, status, mtime):
        import os as _os
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{slug}.md"
        path.write_text(
            f"---\ntype: report\nanalysis_id: {analysis_id}\ndate: 2026-06-24\n"
            f"status: {status}\n---\n\n## 結論\n{marker}\n",
            encoding="utf-8",
        )
        _os.utime(path, (mtime, mtime))
        return path

    def test_read_report_prefers_newest_across_candidates(self):
        slug = notes.make_slug("a-1")
        self._write_report_file(self.tmp / "reports", slug, "a-1", "OLD primary", "draft", 1000)
        self._write_report_file(self.tmp / "reports_fallback", slug, "a-1", "NEW fallback", "final", 2000)
        text = server.read_report("a-1")
        self.assertIn("NEW fallback", text)
        self.assertNotIn("OLD primary", text)

    def test_read_report_without_id_dedupes_keeping_newest(self):
        slug = notes.make_slug("a-1")
        self._write_report_file(self.tmp / "reports", slug, "a-1", "old", "draft", 1000)
        self._write_report_file(self.tmp / "reports_fallback", slug, "a-1", "new", "final", 2000)
        listing = server.read_report()
        self.assertEqual(listing.count("a-1"), 1)
        self.assertIn("status=final", listing)
        self.assertNotIn("status=draft", listing)


if __name__ == "__main__":
    unittest.main()


class ReportDirStaysWithDataTests(unittest.TestCase):
    """解析記録は解析フォルダ側に置く契約。リポジトリ側へ退避する候補を持たない。"""

    def test_candidates_without_override_are_the_analysis_dir_only(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp, config_env(Path(tmp), {}):
            _open_analysis_dir(Path(tmp) / "analysis")
            candidates = mcp_core._report_dir_candidates(session_state.session)  # [paths] reports_dir が無い
        self.assertEqual(candidates, [Path(tmp) / "analysis" / "reports"])
        self.assertNotIn(mcp_core.BASE_DIR / "reports", candidates)

    def test_no_dataset_and_no_setting_is_missing_state(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp, config_env(Path(tmp), {}):
            session_state.session = session_state.AnalysisSession()
            out = json.loads(server.write_report("a-1", "DS", "本文"))
        self.assertEqual(out["error"]["code"], "missing_state")
        self.assertEqual(out["error"]["state"], "dataset")
        self.assertEqual(out["error"]["required_tools"], ["load_dataset"])
