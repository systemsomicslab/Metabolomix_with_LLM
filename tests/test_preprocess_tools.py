# tests/test_preprocess_tools.py
import json
import unittest
import numpy as np
import server
from metabolomix.core import session_state


class TestPreprocessTools(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_list_sample_roles_requires_data(self):
        out = json.loads(server.arf_list_sample_roles())
        self.assertEqual(out["error"]["code"], "missing_state")

    def test_preprocess_sets_session_matrix(self):
        # minimal fake ARF matrix path: inject via monkeypatch of build helper
        session_state.session.arf.filtered_features = [{"AlignedPeakProperties": []}]
        session_state.session.arf.class_index = None

        def fake_matrix(*a, **k):
            m = np.array([[10.0, 1.0], [20.0, 2.0], [15.0, 1.5]])
            return m, ["s1", "s2", "q1"], ["Spot_0_height", "Spot_1_height"]

        from metabolomix.core import tool_helpers
        # _pp_build_matrix の正準定義元は tool_helpers。arf_preprocess はそこを module 修飾で
        # 参照するため、差し替え・復元も tool_helpers 側で行う（他テストへの漏れを防ぐ）。
        _orig_build = tool_helpers._pp_build_matrix
        tool_helpers._pp_build_matrix = fake_matrix
        self.addCleanup(setattr, tool_helpers, "_pp_build_matrix", _orig_build)
        out = json.loads(server.arf_preprocess(normalize="median", impute="half_min"))
        self.assertEqual(out["status"], "success")
        self.assertIsNotNone(session_state.session.arf.feature_matrix)
        self.assertEqual(session_state.session.arf.preprocessing_recipe["normalize"], "median")


class TestPcaPreprocessed(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_requires_preprocessed_matrix(self):
        out = server.arf_pca_preprocessed()
        self.assertIn("前処理", out)

    def test_runs_pca_on_preprocessed_matrix(self):
        import numpy as np
        session_state.session.arf.feature_matrix = np.array(
            [[1.0, 2.0, 3.0], [2.0, 1.0, 0.0], [3.0, 3.0, 3.0], [0.0, 1.0, 2.0]])
        session_state.session.arf.pp_sample_names = ["a", "b", "c", "d"]
        session_state.session.arf.pp_feature_names = ["Spot_0_height", "Spot_1_height", "Spot_2_height"]
        session_state.session.arf.class_index = None
        session_state.session.arf.features = []  # get_pca_loading_features tolerates empty spots
        session_state.session.arf.preprocessing_recipe = {"normalize": "median"}
        self.assertTrue(server._pp_has_preprocessed())
        out = server.arf_pca_preprocessed()
        self.assertIn("PCA", out)
        self.assertIn("前処理レシピ", out)


if __name__ == "__main__":
    unittest.main()


class TestBlankExclusionFromAnalysisMatrix(unittest.TestCase):
    """ブランクは背景除去に使ったあと解析行列から外す。QC は PCA 用に残す。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def _prime(self):
        session_state.session.arf.filtered_features = [{"AlignedPeakProperties": []}]
        session_state.session.arf.class_index = None

        def fake_matrix(*a, **k):
            m = np.array([
                [10.0, 1.0],    # s1   sample
                [12.0, 1.2],    # s2   sample
                [11.0, 1.1],    # QC_1 qc
                [0.05, 0.01],   # blank_1 blank
            ])
            return m, ["s1", "s2", "QC_1", "blank_1"], ["Spot_0_height", "Spot_1_height"]

        from metabolomix.core import tool_helpers
        _orig = tool_helpers._pp_build_matrix
        tool_helpers._pp_build_matrix = fake_matrix
        self.addCleanup(setattr, tool_helpers, "_pp_build_matrix", _orig)

    def test_blank_row_is_removed_but_qc_row_is_kept(self):
        self._prime()
        out = json.loads(server.arf_preprocess(normalize="none", impute="none"))
        self.assertEqual(out["status"], "success")
        self.assertEqual(session_state.session.arf.pp_sample_names, ["s1", "s2", "QC_1"])
        self.assertEqual(session_state.session.arf.feature_matrix.shape[0], 3)

    def test_removal_is_reported_and_caveated(self):
        self._prime()
        out = json.loads(server.arf_preprocess(normalize="none", impute="none"))
        self.assertEqual(out["excluded_from_matrix"], {"blank": ["blank_1"]})
        self.assertTrue(
            any("ブランク" in c for c in out["caveats"]),
            f"ブランク除外の caveat が無い: {out['caveats']}",
        )

    def test_blank_still_usable_for_background_filtering(self):
        """行を外すのは blank_filter を通したあと（背景除去の参照は失わない）。"""
        self._prime()
        out = json.loads(
            server.arf_preprocess(normalize="none", impute="none", blank_min_fold=3.0)
        )
        self.assertIn("blank_filter", out["recipe_applied"])
        self.assertNotIn(
            "ブランクまたは生体試料が無いため",
            json.dumps(out["steps"]["blank_filter"], ensure_ascii=False),
        )
