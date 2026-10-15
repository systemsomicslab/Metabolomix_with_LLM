"""アラインメント単位の兄弟ファイルとスポット別 MS/MS（`metabolomix.msdial.alignment_spot`）。"""
import tempfile
import unittest
from pathlib import Path

from metabolomix.dcl.reader import deserialize_dcl, read_dcl_entry
from metabolomix.msdial import alignment_spot
from tests.dcl_fixture import build_dcl_bytes


def _dcl(path: Path) -> Path:
    path.write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 3.0, "spectrum": [(184.07, 999.0), (86.1, 50.0)]},
        {"precursor_mz": 885.5546, "rt": 3.35, "spectrum": []},
        {"precursor_mz": 500.0, "rt": 4.0, "spectrum": [(100.0, 10.0)]},
    ]))
    return path


class SiblingTests(unittest.TestCase):
    def test_every_alignment_file_points_at_its_own_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("AlignmentResult_1.arf2", "AlignmentResult_1.dcl",
                         "AlignmentResult_1.EIC.aef", "AlignmentResult_2.dcl"):
                (root / name).write_bytes(b"x")
            for start in ("AlignmentResult_1.arf2", "AlignmentResult_1.EIC.aef",
                          "AlignmentResult_1.dcl"):
                self.assertEqual(
                    alignment_spot.alignment_sibling(root / start, ".dcl"),
                    root / "AlignmentResult_1.dcl",
                )
            self.assertEqual(
                alignment_spot.alignment_sibling(root / "AlignmentResult_1.EIC.aef", ".arf2"),
                root / "AlignmentResult_1.arf2",
            )
            self.assertIsNone(
                alignment_spot.alignment_sibling(root / "AlignmentResult_2.dcl", ".arf2")
            )
            self.assertIsNone(alignment_spot.alignment_sibling(root / "other.txt", ".dcl"))


class DclEntryTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dcl = _dcl(Path(self._tmp.name) / "AlignmentResult_1.dcl")

    def test_single_entry_read_matches_the_full_parse(self):
        full = deserialize_dcl(str(self.dcl), include_spectrum=True)
        for index in range(3):
            self.assertEqual(read_dcl_entry(str(self.dcl), index), full[index])
        with self.assertRaises(IndexError):
            read_dcl_entry(str(self.dcl), 3)

    def test_msms_lookup_distinguishes_every_outcome(self):
        ok = alignment_spot.alignment_msms(self.dcl, 0, expected_mz=760.585)
        self.assertEqual(ok["status"], "ok")
        self.assertEqual(ok["n_peaks"], 2)
        self.assertEqual(ok["spectrum"][0], [184.07, 999.0])

        self.assertEqual(
            alignment_spot.alignment_msms(self.dcl, 1, expected_mz=885.5546)["status"],
            "no_msms",
        )
        mismatch = alignment_spot.alignment_msms(self.dcl, 2, expected_mz=600.0)
        self.assertEqual(mismatch["status"], "precursor_mismatch")
        self.assertNotIn("spectrum", mismatch)
        self.assertEqual(
            alignment_spot.alignment_msms(self.dcl, 7, expected_mz=1.0)["status"],
            "out_of_range",
        )

    def test_tolerance_follows_the_measured_spread_between_dcl_and_arf2(self):
        # 実データ 6 バッチ・約 2.1 万スポットで、代表試料の precursor と ARF2 MassCenter の
        # 差は最大 0.027 Da、RT 差は最大 0.106 分。0.02 Da ずれた正しい対応は通す。
        ok = alignment_spot.alignment_msms(self.dcl, 0, expected_mz=760.6051, expected_rt=3.1)
        self.assertEqual(ok["status"], "ok")
        # m/z が合っていても RT が離れていれば別スポットとして拒む（隣接スポットは m/z が近い）
        far = alignment_spot.alignment_msms(self.dcl, 0, expected_mz=760.5851, expected_rt=3.5)
        self.assertEqual(far["status"], "precursor_mismatch")
        self.assertTrue(alignment_spot.precursor_agrees(760.5851, 760.63, 3.0, 3.19))
        self.assertFalse(alignment_spot.precursor_agrees(760.5851, 760.64, 3.0, 3.0))

    def test_peak_counts_for_listing(self):
        counts = alignment_spot.msms_peak_counts(self.dcl)
        self.assertEqual(
            [(c["precursor_mz"], c["rt"], c["n_peaks"]) for c in counts],
            [(760.5851, 3.0, 2), (885.5546, 3.35, 0), (500.0, 4.0, 1)],
        )


if __name__ == "__main__":
    unittest.main()
