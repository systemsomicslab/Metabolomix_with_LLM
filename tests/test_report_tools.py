import json
import os
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # import server が pyplot を読む前にヘッドレス指定

import server
from metabolomix.core import session_state
from metabolomix.core import mcp_core
from metabolomix.corpus import knowledge_store


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
        meta = server._build_report_meta("a-1", "DS", "draft", ["s1", "s2"])
        self.assertEqual(meta["type"], "report")
        self.assertEqual(meta["analysis_id"], "a-1")
        self.assertEqual(meta["dataset"], "DS")
        self.assertEqual(meta["status"], "draft")
        self.assertEqual(meta["knowledge_refs"], ["s1", "s2"])
        self.assertRegex(meta["date"], r"^\d{4}-\d{2}-\d{2}$")

    def test_build_report_meta_defaults_refs_to_empty_list(self):
        meta = server._build_report_meta("a-1", "DS", "draft", None)
        self.assertEqual(meta["knowledge_refs"], [])


class ReportToolTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._saved_data_dir = mcp_core.DATA_DIR
        mcp_core.DATA_DIR = self.tmp
        self._saved_env = os.environ.get("LIPIDMIX_REPORTS_DIR")
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "reports_fallback")

    def tearDown(self):
        mcp_core.DATA_DIR = self._saved_data_dir
        if self._saved_env is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_env
        self._tmp.cleanup()

    def test_write_report_creates_file_with_frontmatter_and_body(self):
        msg = server.write_report("a-1", "NEG / DS", "## 目的\nグループ比較")
        path = self.tmp / "reports_fallback" / "a-1.md"
        self.assertTrue(path.is_file())
        self.assertIn("reports", msg)
        meta, body = knowledge_store.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["type"], "report")
        self.assertEqual(meta["analysis_id"], "a-1")
        self.assertEqual(meta["status"], "draft")
        self.assertIn("## 目的", body)

    def test_write_report_overwrites_on_second_call(self):
        server.write_report("a-1", "DS", "## 目的\n古い本文")
        server.write_report("a-1", "DS", "## 目的\n新しい本文", status="final")
        path = self.tmp / "reports_fallback" / "a-1.md"
        meta, body = knowledge_store.parse_frontmatter(path.read_text(encoding="utf-8"))
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

    def test_list_reports_shows_written_report(self):
        server.write_report("a-1", "DS", "## 目的\nx", status="final")
        listing = server.list_reports()
        self.assertIn("a-1", listing)
        self.assertIn("status=final", listing)

    def test_list_reports_empty(self):
        listing = server.list_reports()
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

    def test_remember_arf_pca_plot_builds_session_state(self):
        saved = session_state.session.arf.last_pca_plot
        try:
            server._remember_arf_pca_plot(
                {"components": [[1.0, 2.0], [3.0, 4.0]],
                 "explained_variance_ratio": [0.5, 0.3]},
                ["s1", "s2"],
                "PCA Score Plot (x.arf)",
            )
            plot = session_state.session.arf.last_pca_plot
            self.assertEqual(plot["title"], "PCA Score Plot (x.arf)")
            self.assertEqual(plot["points"][0], {"x": 1.0, "y": 2.0, "label": "s1"})
            self.assertIn("50.00%", plot["x_label"])
        finally:
            session_state.session.arf.last_pca_plot = saved

    def test_format_pca_plot_block_includes_points_and_summary(self):
        import json as _json
        from metabolomix.core import tool_helpers
        block = tool_helpers._format_pca_plot_block(
            {"components": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]],
             "explained_variance_ratio": [0.24, 0.15]},
            ["s1", "s2", "s3", "s4"],
            title="PCA Score Plot (x.arf)",
            intro="\n#### 📊 PCA スコア要約\n",
            groups={"s1": "A", "s2": "A", "s3": "B", "s4": "B"},
        )
        # 座標点列が JSON フェンスで同梱される
        self.assertIn("```json", block)
        self.assertIn('"sample": "s1"', block)
        self.assertIn('"pc1"', block)
        self.assertIn('"pc2"', block)
        self.assertIn('"group": "A"', block)
        # 結論・群別サンプル数・図示 note も残る
        self.assertIn("24.00%", block)
        self.assertIn("15.00%", block)
        self.assertIn("A=2", block)
        self.assertIn("B=2", block)
        self.assertIn("save_figure", block)
        # JSON は valid で points が 4 件、座標が一致
        payload = _json.loads(block.split("```json")[1].split("```")[0].strip())
        self.assertEqual(len(payload["points"]), 4)
        self.assertEqual(payload["points"][0]["pc1"], 1.0)
        self.assertEqual(payload["points"][0]["pc2"], 2.0)

    def test_format_pca_plot_block_no_groups_includes_points_without_group(self):
        from metabolomix.core import tool_helpers
        block = tool_helpers._format_pca_plot_block(
            {"components": [[1.0, 2.0], [3.0, 4.0]],
             "explained_variance_ratio": [0.5, 0.3]},
            ["s1", "s2"],
            title="T",
            intro="\n#### PCA\n",
        )
        self.assertIn("サンプル数: 2", block)
        self.assertIn('"pc1"', block)
        self.assertNotIn('"group"', block)


class SavePcaFigureTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._saved_data_dir = mcp_core.DATA_DIR
        mcp_core.DATA_DIR = self.tmp
        self._saved_env = os.environ.get("LIPIDMIX_REPORTS_DIR")
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "reports_fallback")
        self._saved_plot = session_state.session.arf.last_pca_plot

    def tearDown(self):
        mcp_core.DATA_DIR = self._saved_data_dir
        session_state.session.arf.last_pca_plot = self._saved_plot
        if self._saved_env is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_env
        self._tmp.cleanup()

    def test_save_pca_figure_writes_png_and_returns_relpath(self):
        session_state.session.arf.last_pca_plot = {
            "title": "T", "x_label": "PC1", "y_label": "PC2",
            "points": [
                {"x": 1.0, "y": 2.0, "label": "s1"},
                {"x": -1.0, "y": 0.5, "label": "s2"},
            ],
        }
        msg = server.save_figure("pca", "a-1")
        png = self.tmp / "reports_fallback" / "figures" / "a-1_pca.png"
        self.assertTrue(png.is_file())
        self.assertIn("figures/a-1_pca.png", msg)

    def test_save_pca_figure_guidance_when_no_plot(self):
        session_state.session.arf.last_pca_plot = None
        msg = server.save_figure("pca", "a-1")
        self.assertIn("PCA", msg)
        self.assertFalse((self.tmp / "reports_fallback" / "figures").exists())


class ReportEdgeCaseTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._saved_data_dir = mcp_core.DATA_DIR
        mcp_core.DATA_DIR = self.tmp
        self._saved_env = os.environ.get("LIPIDMIX_REPORTS_DIR")
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "reports_fallback")

    def tearDown(self):
        mcp_core.DATA_DIR = self._saved_data_dir
        if self._saved_env is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_env
        self._tmp.cleanup()

    def test_write_report_refuses_slug_collision_with_different_id(self):
        # Both analysis_ids collapse to the same slug.
        slug = knowledge_store.make_slug("Group A vs B")
        self.assertEqual(slug, knowledge_store.make_slug("group-a-vs-b"))

        server.write_report("Group A vs B", "DS", "## 目的\nfirst")
        msg = server.write_report("group-a-vs-b", "DS", "## 目的\nsecond")

        self.assertIn("衝突", msg)
        # The original report must be intact (not overwritten).
        path = self.tmp / "reports_fallback" / f"{slug}.md"
        meta, body = knowledge_store.parse_frontmatter(path.read_text(encoding="utf-8"))
        self.assertEqual(meta["analysis_id"], "Group A vs B")
        self.assertIn("first", body)
        self.assertNotIn("second", body)

    def test_write_report_same_id_still_overwrites(self):
        # Regression: the collision guard must NOT block same-id overwrites.
        server.write_report("a-1", "DS", "## 目的\n古い")
        msg = server.write_report("a-1", "DS", "## 目的\n新しい", status="final")
        self.assertIn("保存", msg)
        path = self.tmp / "reports_fallback" / "a-1.md"
        meta, body = knowledge_store.parse_frontmatter(path.read_text(encoding="utf-8"))
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
        slug = knowledge_store.make_slug("a-1")
        self._write_report_file(self.tmp / "reports", slug, "a-1", "OLD primary", "draft", 1000)
        self._write_report_file(self.tmp / "reports_fallback", slug, "a-1", "NEW fallback", "final", 2000)
        text = server.read_report("a-1")
        self.assertIn("NEW fallback", text)
        self.assertNotIn("OLD primary", text)

    def test_list_reports_dedupes_keeping_newest(self):
        slug = knowledge_store.make_slug("a-1")
        self._write_report_file(self.tmp / "reports", slug, "a-1", "old", "draft", 1000)
        self._write_report_file(self.tmp / "reports_fallback", slug, "a-1", "new", "final", 2000)
        listing = server.list_reports()
        self.assertEqual(listing.count("a-1"), 1)
        self.assertIn("status=final", listing)
        self.assertNotIn("status=draft", listing)


if __name__ == "__main__":
    unittest.main()


class DatasetFigureFallbackTests(unittest.TestCase):
    """mzTab-M（DatasetState）経路からも図を保存できること。

    session.arf 側だけを見ていたため、ds.last_pca / ds.last_differential に
    全量があるのに mzTab-M 経路からはレポート用の図が作れなかった。
    """

    def setUp(self):
        import tempfile
        from metabolomix.mztab.dataset_state import DatasetState
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._saved_data_dir = mcp_core.DATA_DIR
        mcp_core.DATA_DIR = self.tmp
        self._saved_env = os.environ.get("LIPIDMIX_REPORTS_DIR")
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "reports_fallback")
        self._saved_plot = session_state.session.arf.last_pca_plot
        self._saved_diff = session_state.session.arf.last_differential
        self._saved_ds = session_state.session.dataset
        session_state.session.arf.last_pca_plot = None
        session_state.session.arf.last_differential = None
        session_state.session.dataset = DatasetState()

    def tearDown(self):
        mcp_core.DATA_DIR = self._saved_data_dir
        session_state.session.arf.last_pca_plot = self._saved_plot
        session_state.session.arf.last_differential = self._saved_diff
        session_state.session.dataset = self._saved_ds
        if self._saved_env is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_env
        self._tmp.cleanup()

    def test_save_pca_figure_falls_back_to_dataset_state(self):
        session_state.session.dataset.last_pca = {
            "explained_variance_ratio": [0.31, 0.22],
            "scores": [
                {"name": "s1", "role": "sample", "PC1": 1.0, "PC2": 2.0},
                {"name": "s2", "role": "sample", "PC1": -1.0, "PC2": 0.5},
            ],
            "n_samples": 2, "n_features": 10, "log_transform": False,
        }
        msg = server.save_figure("pca", "ds-1")
        png = self.tmp / "reports_fallback" / "figures" / "ds-1_pca.png"
        self.assertTrue(png.is_file())
        self.assertIn("figures/ds-1_pca.png", msg)
        # どちらの経路で描いたかを明示する（ARF と混同させない）
        self.assertIn("mztab", msg)

    def test_save_volcano_figure_falls_back_to_dataset_state(self):
        session_state.session.dataset.last_differential = {
            "kind": "two_group", "a": "ctrl", "b": "treat", "n_a": 3, "n_b": 3,
            "q_threshold": 0.05, "log2fc_threshold": 1.0, "log_transform": True,
            # 実形状: differential.volcano_data() の戻り値そのまま（点の list）
            "volcano": [
                {"feature": "f1", "log2fc": 1.5, "neg_log10_p": 2.0, "sig": "up"},
                {"feature": "f2", "log2fc": -1.8, "neg_log10_p": 2.2, "sig": "down"},
                {"feature": "f3", "log2fc": 0.1, "neg_log10_p": 0.2, "sig": "ns"},
            ],
        }
        msg = server.save_figure("volcano", "ds-2")
        png = self.tmp / "reports_fallback" / "figures" / "ds-2_volcano.png"
        self.assertTrue(png.is_file())
        self.assertIn("figures/ds-2_volcano.png", msg)
        self.assertIn("mztab", msg)

    def _both_paths_have_pca(self):
        session_state.session.arf.last_pca_plot = {
            "title": "ARF", "x_label": "PC1", "y_label": "PC2",
            "points": [{"x": 0.0, "y": 0.0, "label": "arf-sample"}],
        }
        session_state.session.dataset.last_pca = {
            "explained_variance_ratio": [0.9, 0.05],
            "scores": [{"name": "ds-sample", "role": "sample", "PC1": 5.0, "PC2": 5.0}],
            "n_samples": 1, "n_features": 3, "log_transform": False,
        }

    def test_two_valid_results_stop_instead_of_picking_one(self):
        """両方載っているときに黙って ARF を採ると、前のデータの図が保存される。

        どちらを描くかは図の数字そのものを変える。優先順位を決め打ちすると、
        それが暗黙の既定になって「なぜこの図なのか」を誰も説明できなくなる。
        """
        self._both_paths_have_pca()
        parsed = json.loads(server.save_figure("pca", "both-1"))
        self.assertEqual(parsed["error"]["code"], "AMBIGUOUS_RESULT_SOURCE")
        self.assertEqual(
            sorted(c["source"] for c in parsed["error"]["details"]["candidates"]),
            ["arf", "mztab"])
        self.assertFalse((self.tmp / "reports_fallback" / "figures" / "both-1_pca.png").exists())

    def test_an_explicit_source_resolves_the_ambiguity(self):
        self._both_paths_have_pca()
        msg = server.save_figure("pca", "both-2", source="mztab")
        self.assertIn("source=mztab", msg)
        self.assertTrue((self.tmp / "reports_fallback" / "figures" / "both-2_pca.png").is_file())

    def test_guidance_mentions_dataset_tools_when_nothing_available(self):
        session_state.session.dataset = None
        msg = server.save_figure("pca", "none-1")
        self.assertIn("dataset_pca", msg)
