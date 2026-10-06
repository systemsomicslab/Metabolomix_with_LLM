"""ARF 読み込みのホットパスに対する計算量の契約。

`.arf` は spot × 注入の直積で行が増えるため、1 行あたりの定数倍と「同じ仕事を
何度やるか」がそのまま所要時間になる（実測: 393MB / 16,596 spot の `.arf` で
`arf_parser` が 96.6 秒。うち大半が本ファイルで固定する 3 経路）。

ここで固定するのは速度そのものではなく、速度を決めている**構造**:

1. `normalize_sample_name` は純関数なのでメモ化されていること
   （spot × 注入ぶん呼ばれ、異なり値は注入数しかない）。
2. 行数を知るためだけに pandas の DataFrame を組まないこと
   （`count_peak_property_rows` が `extract_peak_properties` と同じ数を返す）。
3. `exclusions._entry_file_name` が file_name 1 個のために
   `_convert_to_alignment_feature` の dict を丸ごと組まないこと。

いずれも「同じ答えを返し続けること」を対で固定してある。速い経路が違う答えを
返したらテストが落ちる。
"""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from metabolomix.arf import exclusions
from metabolomix.arf import reader as arf_reader
from metabolomix.msdial import tags as msdial_tags


def _row(file_id, name, height):
    """`AlignmentChromPeakFeature` の生 row（len>25・data[18] が数値）。

    tests/test_exclusions.py の `_row` と同じ形。これでないと
    `_convert_to_alignment_feature` が実経路に入らない。
    """
    row = [0] * 26
    row[0] = file_id
    row[1] = name
    row[2] = file_id
    row[18] = float(height)
    return row


def _spot(master_id, rows):
    return {"MasterAlignmentID": master_id,
            "AlignedPeakProperties": [list(r) for r in rows]}


def _fixture():
    return [
        _spot(1, [_row(0, "sA.wiff", 10), _row(1, "sB.wiff", 20)]),
        _spot(2, [_row(0, "sA.wiff", 11), _row(1, "sB.wiff", 21)]),
    ]


class NormalizeSampleNameMemoTests(unittest.TestCase):
    """spot × 注入ぶん呼ばれる純関数なので、同じ入力を毎回計算し直さない。"""

    def test_repeated_call_with_same_argument_is_served_from_cache(self):
        msdial_tags.normalize_sample_name.cache_clear()
        msdial_tags.normalize_sample_name("20240207_A1_POS_1.wiff")
        before = msdial_tags.normalize_sample_name.cache_info().hits
        msdial_tags.normalize_sample_name("20240207_A1_POS_1.wiff")
        after = msdial_tags.normalize_sample_name.cache_info().hits
        self.assertEqual(after, before + 1)

    def test_cache_does_not_confuse_the_two_timestamp_policies(self):
        """キーワード引数はキーの一部。混ざると片方の答えが他方に漏れる。"""
        msdial_tags.normalize_sample_name.cache_clear()
        name = "sample_A_202605151012_tags.xml"
        self.assertEqual(msdial_tags.normalize_sample_name(name), "sample_a")
        self.assertEqual(
            msdial_tags.normalize_sample_name(name, strip_processing_timestamp=False),
            "sample_a_202605151012",
        )
        self.assertEqual(msdial_tags.normalize_sample_name(name), "sample_a")

    def test_cached_results_match_the_documented_normalizations(self):
        msdial_tags.normalize_sample_name.cache_clear()
        cases = {
            "sample_A": "sample_a",
            "sample_A.wiff": "sample_a",
            "C:/data/QC01.d": "qc01",
            "sample_A_202605151012_tags.xml": "sample_a",
        }
        for _ in range(2):  # 2 周目はキャッシュ経由。同じ答えでなければならない
            for value, expected in cases.items():
                self.assertEqual(msdial_tags.normalize_sample_name(value), expected)

    def test_none_and_path_arguments_still_work(self):
        from pathlib import Path
        msdial_tags.normalize_sample_name.cache_clear()
        self.assertEqual(msdial_tags.normalize_sample_name(None), "")
        self.assertEqual(msdial_tags.normalize_sample_name(Path("x/sample_A.wiff")),
                         "sample_a")


class CountPeakPropertyRowsTests(unittest.TestCase):
    """`arf_parser` は行数しか使わないのに DataFrame を組んでいた。"""

    def test_count_matches_extract_peak_properties_length(self):
        spots = _fixture()
        self.assertEqual(arf_reader.count_peak_property_rows(spots),
                         len(arf_reader.extract_peak_properties(spots)))

    def test_count_skips_the_same_rows_extract_skips(self):
        """短すぎる row・list でない row・AlignedPeakProperties 無しの spot。"""
        spots = [
            _spot(1, [_row(0, "sA", 10)]),
            {"MasterAlignmentID": 2},                       # キー自体が無い
            {"MasterAlignmentID": 3, "AlignedPeakProperties": None},
            {"MasterAlignmentID": 4, "AlignedPeakProperties": [
                [1, "short"],                               # len < 3
                "not-a-list",
                _row(1, "sB", 20),
            ]},
        ]
        self.assertEqual(arf_reader.count_peak_property_rows(spots),
                         len(arf_reader.extract_peak_properties(spots)))
        self.assertEqual(arf_reader.count_peak_property_rows(spots), 2)

    def test_empty_input_counts_zero(self):
        self.assertEqual(arf_reader.count_peak_property_rows([]), 0)


def _row_battery():
    """`peak_field` / `file_name_of` の同値性を試す行の一式。

    実経路に入る行と、入らない行（短い・list でない・数値でない data[18]）と、
    欠けた行（data[37] が無い・文字列が 1 つも無い）を混ぜる。
    """
    full = [0] * 40
    full[0], full[1], full[2], full[3] = 7, "sA.wiff", -2, 11   # -2 = ギャップフィル
    full[9] = 900
    full[10] = {1: 35.0}
    full[15] = [[1, [4.25]], [3, [700.5]]]
    full[18], full[20], full[21], full[22] = 1234.5, 2000.0, 1800.0, 700.5
    full[37] = [1.5, 22.0]

    ri_only = list(full)
    ri_only[15] = 0                      # data[15] が list でない → data[16] を見る
    ri_only[16] = [[1, [9.75]]]

    no_shape = [0] * 30                  # len <= 37 → PeakShape が取れない
    no_shape[1], no_shape[2], no_shape[18] = b"sB.raw", 5, 99.0

    nameless = [0] * 26                  # 先頭付近に文字列が無い
    nameless[18] = 1.0

    return [full, ri_only, no_shape, nameless,
            [0] * 26,                    # 数値 data[18] だが名前なし（別インスタンス）
            [1, "short"],                # len < 3
            [0] * 30,                     # data[18] が 0（int なので実経路に入る）
            "not-a-list",
            None]


class PeakFieldTests(unittest.TestCase):
    """1 フィールドのために 16 キーの dict を組まない。答えは dict と同じ。"""

    _KEYS = ("peak_id", "file_id", "file_name", "master_peak_id", "height", "area",
             "area_above_baseline", "m_z", "rt", "signal_to_noise", "estimated_noise",
             "ms2_raw_id", "is_msms", "is_gap_filled", "is_msms_matched", "is_matched")

    def test_every_key_matches_the_full_feature_dict(self):
        for row in _row_battery():
            feature = (arf_reader.alignment_feature_row(row)
                       if isinstance(row, list) else {})
            for key in self._KEYS:
                with self.subTest(row=row, key=key):
                    self.assertEqual(arf_reader.peak_field(row, key), feature.get(key))

    def test_unknown_key_is_none_like_the_dict(self):
        row = _row_battery()[0]
        self.assertIsNone(arf_reader.peak_field(row, "no_such_field"))

    def test_does_not_build_the_full_feature_dict(self):
        row = _row_battery()[0]
        with patch.object(arf_reader, "_convert_to_alignment_feature",
                          side_effect=AssertionError(
                              "peak_field が feature dict を組んでいる")):
            self.assertEqual(arf_reader.peak_field(row, "height"), 1234.5)
            self.assertEqual(arf_reader.peak_field(row, "file_name"), "sA.wiff")
            self.assertIs(arf_reader.peak_field(row, "is_gap_filled"), True)


class BuildPcaMatrixAvoidsFeatureDictTests(unittest.TestCase):
    """行列組み立ては file_name / is_gap_filled / props しか読まない。"""

    def test_matrix_is_unchanged_without_the_feature_dict(self):
        spots = _pca_fixture()
        expected = arf_reader.build_pca_matrix(spots, use_properties=["height"])
        with patch.object(arf_reader, "_convert_to_alignment_feature",
                          side_effect=AssertionError(
                              "build_pca_matrix が feature dict を組んでいる")):
            got = arf_reader.build_pca_matrix(spots, use_properties=["height"])
        self.assertEqual(got[1], expected[1])
        self.assertEqual(got[2], expected[2])
        self.assertTrue((got[0] == expected[0]).all())

    def test_detection_rate_filter_still_counts_gap_filled_rows(self):
        """`min_detection_rate` は is_gap_filled を数える。速い経路でも同じ。"""
        spots = _pca_fixture()
        for row in spots[0]["AlignedPeakProperties"]:
            row[2] = -2  # spot 0 の全注入をギャップフィルにする
        expected = arf_reader.build_pca_matrix(
            spots, use_properties=["height"], min_detection_rate=0.5)
        with patch.object(arf_reader, "_convert_to_alignment_feature",
                          side_effect=AssertionError("feature dict を組んでいる")):
            got = arf_reader.build_pca_matrix(
                spots, use_properties=["height"], min_detection_rate=0.5)
        self.assertEqual(got[2], expected[2])
        self.assertNotIn("Spot_0_height", got[2])  # 検出率 0 で落ちている

    def test_non_default_property_takes_the_same_path(self):
        spots = _pca_fixture()
        expected = arf_reader.build_pca_matrix(spots, use_properties=["m_z"])
        with patch.object(arf_reader, "_convert_to_alignment_feature",
                          side_effect=AssertionError("feature dict を組んでいる")):
            got = arf_reader.build_pca_matrix(spots, use_properties=["m_z"])
        self.assertEqual(got[1], expected[1])
        self.assertEqual(got[2], expected[2])


def _pca_fixture():
    """`build_pca_matrix` が実際に消費できる 4 サンプル × 3 スポット。

    形は tests/test_exclude_wiring.py の fixture と同じ（生 list 行・data[18] が
    height・サンプル間で分散を持つ）。
    """
    spots = []
    for mid in (0, 1, 2):
        rows = []
        for i, name in enumerate(["sA", "sB", "sC", "sD"]):
            row = [i] * 40
            row[1] = name
            row[2] = 100 + mid
            row[18] = float(10 * (mid + 1) + i)
            rows.append(row)
        spots.append({"MasterAlignmentID": mid, "Name": f"Lipid_{mid}",
                      "MassCenter": 700.0 + mid, "RT": 5.0 + mid,
                      "AlignedPeakProperties": rows})
    return spots


class _FakeArfState:
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


class _FakeSession:
    def __init__(self, spots):
        self.arf = _FakeArfState(spots)

    def maybe_prepend_caveat(self, text, topic=None):
        return text


class ArfParserAvoidsThrowawayDataFrameTests(unittest.TestCase):
    """`arf_parser` は行数しか使わないので DataFrame を組んではいけない。"""

    def _run(self):
        import server  # ツール登録の副作用を踏む（他テストと同じ入口）
        from metabolomix.core import session_state as state_mod
        original = state_mod.session
        state_mod.session = _FakeSession(_pca_fixture())
        try:
            with tempfile.TemporaryDirectory() as tmp:
                arf_path = Path(tmp) / "dummy.arf"
                arf_path.touch()
                return server.arf_parser(str(arf_path))
        finally:
            state_mod.session = original

    def test_reports_the_same_record_count_without_calling_extract(self):
        expected = len(arf_reader.extract_peak_properties(_pca_fixture()))
        with patch.object(arf_reader, "extract_peak_properties",
                          side_effect=AssertionError(
                              "arf_parser が行数のために DataFrame を組んでいる")):
            out = self._run()
        self.assertIn(f"**抽出された総ピークレコード数**: {expected}", out)
        self.assertIn("**平均サンプル数/スポット**: 4.00", out)


class EntryFileNameTests(unittest.TestCase):
    """除外判定のキー導出は `_convert_to_alignment_feature` と同じ答えのまま、
    dict を組まずに済ませる。"""

    def test_agrees_with_alignment_feature_row(self):
        cases = [
            _row(0, "sA.wiff", 10),
            _row(3, b"sB.wiff", 20),                 # msgpack の bytes
            [0] * 26,                                # 文字列が無い（file_name None）
            [1, "sC"],                               # 短すぎる（feature が {}）
            "not-a-list",
        ]
        for row in cases:
            with self.subTest(row=row):
                self.assertEqual(
                    exclusions._entry_file_name(row),
                    arf_reader.alignment_feature_row(row).get("file_name")
                    if isinstance(row, list) else None,
                )

    def test_does_not_build_the_full_feature_dict(self):
        row = _row(0, "sA.wiff", 10)
        with patch.object(arf_reader, "_convert_to_alignment_feature",
                          side_effect=AssertionError(
                              "_entry_file_name が feature dict を組んでいる")):
            self.assertEqual(exclusions._entry_file_name(row), "sA.wiff")

    def test_roster_stays_correct_without_the_feature_dict(self):
        names, ids = exclusions.roster(_fixture())
        self.assertEqual(names, {"sA.wiff", "sB.wiff"})
        self.assertEqual(ids, {1, 2})


if __name__ == "__main__":
    unittest.main()
