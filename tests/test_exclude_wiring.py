import json
import tempfile
import unittest
from pathlib import Path

import server
from metabolomix.core import session_state


class _FakeParserArfState:
    """arf_parser を実 build_pca_matrix で回すための最小 ARF スロット。

    load_data はファイルを読まず fixture をそのまま返す（discover_* を回避）。
    arf_parser が触る属性のみ用意する。
    """

    def __init__(self, spots):
        self.features = spots
        self.filtered_features = spots
        self.current_file_path = "dummy.arf"
        self.current_tag_directory = None
        self.tag_index = {}
        self.class_index = None
        self.excluded_samples = set()
        self.excluded_spots = set()
        self.last_pca_plot = None

    def load_data(self, file_path, tag_directory=None):
        self.current_file_path = file_path
        return self.features


class _FakeParserSession:
    def __init__(self, spots):
        self.arf = _FakeParserArfState(spots)

    def maybe_prepend_caveat(self, text, topic=None):
        return text  # 意味論 caveat / トピック誘導は本テストの対象外（恒等パススルー）


def _spot(master_id, entries):
    # loadings 整形（get_pca_loading_features）が spot.Name/MassCenter/RT を参照するため、
    # Spot メタデータも付す。MasterAlignmentID は features のリスト位置に一致させる。
    return {
        "MasterAlignmentID": master_id,
        "Name": f"Lipid_{master_id}",
        "MassCenter": 700.0 + master_id,
        "RT": 5.0 + master_id,
        "AlignedPeakProperties": list(entries),
    }


def _fixture():
    # 4 サンプル x 3 スポット。build_pca_matrix が実際に消費できる生 list 行。
    # MasterAlignmentID は 0..2（実データ同様 features のリスト位置＝group_index に一致）。
    samples = ["sA", "sB", "sC", "sD"]
    spots = []
    for mid in (0, 1, 2):
        rows = []
        for i, name in enumerate(samples):
            row = [i] * 40  # data[18]=height を確保する長さ
            row[1] = name           # file_name（先頭付近の文字列）
            row[2] = 100 + mid       # master_peak_id（>=0 → 非ギャップフィル）
            row[18] = float(10 * (mid + 1) + i)  # height（サンプル間で分散を持たせる）
            rows.append(row)
        spots.append(_spot(mid, rows))
    return spots


class TestPreprocessExcludeWiring(unittest.TestCase):
    def setUp(self):
        session_state.session = server.AnalysisSession()
        session_state.session.arf.filtered_features = _fixture()
        session_state.session.arf.class_index = None

    def test_preprocess_without_exclusion_uses_all_samples(self):
        out = json.loads(server.arf_preprocess())
        # matrix_shape = [n_samples, n_features]
        self.assertEqual(out["matrix_shape"][0], 4)

    def test_preprocess_honors_excluded_sample(self):
        session_state.session.arf.excluded_samples.add("sB")
        out = json.loads(server.arf_preprocess())
        self.assertEqual(out["matrix_shape"][0], 3)
        self.assertTrue(any("手動除外" in c for c in out.get("caveats", [])))

    def test_preprocess_discloses_zero_exclusions_when_none_applied(self):
        # #6: excluded_samples/excluded_spots が空のとき、開示自体が消えては
        # いけない。arf_parser/load_dataset の再実行で reset_analysis() が
        # 除外集合を黙ってゼロ化しても、「除外はゼロ件」という事実そのものは
        # ここで積極的に述べられる必要がある。
        out = json.loads(server.arf_preprocess())
        self.assertTrue(
            any("現在なし" in c for c in out.get("caveats", [])),
            f"ゼロ件の手動除外が開示されていない: {out.get('caveats')}",
        )

    def test_preprocess_honors_excluded_spot(self):
        session_state.session.arf.excluded_spots.add(1)
        out = json.loads(server.arf_preprocess())
        # 3 スポット → 2 スポット（列数が減る。各スポット1プロパティ height）
        self.assertEqual(out["matrix_shape"][1], 2)

    def test_preprocess_all_samples_excluded_flags_empty_matrix(self):
        for name in ("sA", "sB", "sC", "sD"):
            session_state.session.arf.excluded_samples.add(name)
        out = json.loads(server.arf_preprocess())
        self.assertEqual(out["matrix_shape"][0], 0)
        self.assertTrue(any("行列が空" in c for c in out.get("caveats", [])))


class TestParserExcludeWiring(unittest.TestCase):
    """arf_parser も手動除外(arf_exclude)を PCA 前に反映すること（#1 で統一）。"""

    def setUp(self):
        session_state.session = _FakeParserSession(_fixture())

    def test_parser_honors_excluded_sample(self):
        session_state.session.arf.excluded_samples.add("sB")
        with tempfile.TemporaryDirectory() as tmp:
            arf_path = Path(tmp) / "dummy.arf"
            arf_path.touch()
            text = server.arf_parser(str(arf_path))
        # 行列形状 (サンプル数 x 特徴量数) に 3 サンプルが反映される
        self.assertIn("(3,", text)
        self.assertIn("手動除外", text)


if __name__ == "__main__":
    unittest.main()
