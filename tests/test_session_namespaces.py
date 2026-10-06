"""パーサ間で解析状態が混線しないことの回帰テスト。

背景: PAI2（単一サンプル）用ツールと ARF（多サンプル）用ツールが AnalysisSession の
同じ `features` / `filtered_features` / `current_file_path` を共有していたため、
`pai2_parser` を1回呼ぶだけで進行中の ARF 解析（前処理行列・差次的結果・手動除外）が
無言で消えていた。実際の WebUI セッションでは

    arf_preprocess(成功) → arf_differential(成功) → pai2_parser → arf_pca_preprocessed
    → 「前処理後の行列がありません」

という順で再現し、しかも LLM には破棄の事実が伝わらないため「手動除外が過度」という
誤った原因が報告された。ここではその順序をそのまま固定する。
"""
import json
import os
import tempfile
import unittest

import numpy as np

from metabolomix.pai2 import reader as pai2_reader
import server
from metabolomix.core import session_state
from metabolomix.pai2 import tools as tools_pai2


FAKE_PAI2_PEAKS = [
    {
        "id": 1,
        "name": "PC 34:1",
        "m/z": 760.5851,
        "time": {"rt": 5.12},
        "peak_height": 12000.0,
        "ion_mode": "Negative",
    },
    {
        "id": 2,
        "name": "Unknown",
        "m/z": 281.2486,
        "time": {"rt": 6.40},
        "peak_height": 3100.0,
        "ion_mode": "Negative",
    },
]

# loadings 整形（get_pca_loading_features）が spot.Name/MassCenter/RT を参照するため
# メタデータも付す。MasterAlignmentID は features のリスト位置に一致させる。
ARF_SPOTS = [
    {"MasterAlignmentID": 0, "Name": "PE 36:2", "MassCenter": 742.5386, "RT": 7.81,
     "AlignedPeakProperties": []},
    {"MasterAlignmentID": 1, "Name": "PS 38:4", "MassCenter": 810.5291, "RT": 6.44,
     "AlignedPeakProperties": []},
    {"MasterAlignmentID": 2, "Name": "Unknown", "MassCenter": 885.5499, "RT": 5.02,
     "AlignedPeakProperties": []},
]


class Pai2ArfIsolationTestCase(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        handle, path = tempfile.mkstemp(suffix=".pai2")
        os.write(handle, b"fake-pai2-bytes")
        os.close(handle)
        self.pai2_path = path
        self.addCleanup(os.unlink, path)

        # tools_pai2 は resolve_pai2_file_path を名前で束縛済みなので、そちらを差し替える。
        original_resolver = tools_pai2.resolve_pai2_file_path
        tools_pai2.resolve_pai2_file_path = lambda file_path=None: file_path
        self.addCleanup(
            setattr, tools_pai2, "resolve_pai2_file_path", original_resolver
        )

        # deserialize は関数内で from pai2_reader import deserialize されるため module 側を差し替える。
        original_deserialize = pai2_reader.deserialize
        pai2_reader.deserialize = lambda _stream: [dict(p) for p in FAKE_PAI2_PEAKS]
        self.addCleanup(setattr, pai2_reader, "deserialize", original_deserialize)

    def prime_completed_arf_analysis(self):
        """arf_parser → arf_preprocess まで終わった状態を再現する。"""
        s = session_state.session.arf
        s.features = [dict(spot) for spot in ARF_SPOTS]
        s.filtered_features = s.features
        s.current_file_path = "C:/data/AlignmentResult_2026_07_09.arf"
        s.class_index = None
        s.feature_matrix = np.array(
            [[1.0, 2.0, 3.0], [2.0, 1.0, 0.0], [3.0, 3.0, 3.0], [0.0, 1.0, 2.0]]
        )
        s.pp_sample_names = ["9w_a", "9w_b", "24M_a", "24M_b"]
        s.pp_feature_names = ["Spot_0_height", "Spot_1_height", "Spot_2_height"]
        s.preprocessing_recipe = {"normalize": "none", "impute": "half_min"}
        s.last_differential = {"kind": "two_group", "group_a": "9w", "group_b": "24M"}
        s.excluded_samples = {"190628_QC_Neg_05"}

    def run_pai2_parser(self):
        out = server.pai2_parser(file_path=self.pai2_path)
        self.assertNotIn("[ERROR]", out, f"pai2_parser 自体が失敗している: {out}")
        return out


class TestPai2ParserPreservesArfState(Pai2ArfIsolationTestCase):
    def test_preprocessed_matrix_survives_pai2_parser(self):
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        self.assertIsNotNone(
            session_state.session.arf.feature_matrix,
            "pai2_parser が ARF の前処理後行列を破棄した",
        )
        self.assertEqual(
            session_state.session.arf.pp_sample_names, ["9w_a", "9w_b", "24M_a", "24M_b"]
        )

    def test_differential_result_survives_pai2_parser(self):
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        self.assertIsNotNone(
            session_state.session.arf.last_differential,
            "pai2_parser が差次的解析の結果を破棄した（volcano 再描画が不能になる）",
        )

    def test_manual_exclusions_survive_pai2_parser(self):
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        self.assertEqual(session_state.session.arf.excluded_samples, {"190628_QC_Neg_05"})

    def test_arf_pca_preprocessed_still_runs_after_pai2_parser(self):
        """WebUI セッション 52ecf212 の #79 → #85 → #86 をそのまま再現する。"""
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        out = server.arf_pca_preprocessed()
        self.assertNotIn("前処理後の行列がありません", out)
        self.assertIn("PCA", out)

    def test_arf_list_classes_still_runs_after_pai2_parser(self):
        """current_file_path が .pai2 に差し替わると ARF 系ツールが総崩れする。"""
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        out = server.arf_list_classes()
        self.assertNotIn("先に arf_parser を実行して", out)

    def test_arf_list_sample_roles_still_runs_after_pai2_parser(self):
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        out = server.arf_list_sample_roles()
        self.assertNotIn(
            "missing_state", out,
            "pai2_parser 後に ARF のサンプル一覧が引けなくなっている",
        )
        self.assertIn("# サンプル役割一覧", out)


class TestPai2ToolsUseTheirOwnState(Pai2ArfIsolationTestCase):
    def test_pai2_inspect_peak_reads_pai2_peaks_not_arf_spots(self):
        self.prime_completed_arf_analysis()
        self.run_pai2_parser()
        out = json.loads(server.pai2_inspect_peak(peak_name="PC 34:1"))
        self.assertNotEqual(out.get("status"), "not_found")

    def test_pai2_inspect_peak_without_pai2_parser_is_an_error(self):
        """ARF を読んだだけで PAI2 ツールが「データあり」と誤認してはいけない。"""
        self.prime_completed_arf_analysis()
        out = json.loads(server.pai2_inspect_peak(peak_name="PC 34:1"))
        self.assertEqual(out["error"]["code"], "missing_state")
        self.assertIn("pai2_parser", out["error"].get("message", ""))


class TestArf2ParserPreservesArfState(unittest.TestCase):
    """arf2_parser も features / current_file_path を上書きしていた。"""

    def setUp(self):
        session_state.session = server.AnalysisSession()

    def test_arf2_parser_does_not_overwrite_arf_features(self):
        from metabolomix.arf2 import reader as arf2_reader
        from metabolomix.arf2 import tools as tools_arf2

        handle, path = tempfile.mkstemp(suffix=".arf2")
        os.write(handle, b"fake-arf2-bytes")
        os.close(handle)
        self.addCleanup(os.unlink, path)

        original_resolver = tools_arf2.resolve_arf2_file_path
        tools_arf2.resolve_arf2_file_path = lambda file_path=None: file_path
        self.addCleanup(
            setattr, tools_arf2, "resolve_arf2_file_path", original_resolver
        )
        original_deserialize = arf2_reader.deserialize
        arf2_reader.deserialize = lambda _stream: [
            {"MasterAlignmentID": 9, "Name": "PE 36:2"}
        ]
        self.addCleanup(setattr, arf2_reader, "deserialize", original_deserialize)

        s = session_state.session.arf
        s.features = [dict(spot) for spot in ARF_SPOTS]
        s.filtered_features = s.features
        s.current_file_path = "C:/data/AlignmentResult_2026_07_09.arf"
        s.feature_matrix = np.array([[1.0, 2.0], [3.0, 4.0]])

        server.arf2_parser(file_path=path)

        self.assertIsNotNone(
            session_state.session.arf.feature_matrix,
            "arf2_parser が ARF の前処理後行列を破棄した",
        )
        self.assertEqual(
            len(session_state.session.arf.features),
            len(ARF_SPOTS),
            "arf2_parser が ARF の features を上書きした",
        )


if __name__ == "__main__":
    unittest.main()
