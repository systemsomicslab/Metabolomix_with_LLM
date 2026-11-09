"""m/z を欠いた `.EIC.aef`（MS-DIAL 2023-11-20 以前の出力）を ARF2 で補う。

上流 `DataObjConverter` は 2023-11-20（MsdialWorkbench 20dfa60b4）まで EIC の
スポット中心 m/z を `ChromXsTop.Mz` の平均で書いていた。LC-MS ではこの値が未設定の
-1 なので、その頃の `.EIC.aef` は全スポットで m/z=-1 になる。RT は正しい。
"""
import json
import os
import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from metabolomix.eic.mz_backfill import backfill_missing_mz, sibling_arf2_path
from metabolomix.eic.reader import summarize_eic_data
from tests.config_fixtures import config_env


def _spot(spot_id, rt, mz):
    return {"spot_id": spot_id, "rt": rt, "mz": mz, "num_samples": 0, "samples": []}


def _record(spot_id, rt, mz):
    return {"AlignmentID": spot_id, "RT": rt, "MassCenter": mz}


class BackfillTests(unittest.TestCase):
    def test_fills_missing_mz_from_arf2_when_rt_agrees(self):
        spots = [_spot(0, 3.0, -1.0), _spot(1, 5.0, -1.0)]
        result = backfill_missing_mz(spots, [_record(0, 3.0, 760.5), _record(1, 5.0, 885.5)])
        self.assertEqual([s["mz"] for s in spots], [760.5, 885.5])
        self.assertEqual(result, {"missing": 2, "filled": 2})

    def test_leaves_spot_missing_when_rt_disagrees_or_record_absent(self):
        spots = [_spot(0, 3.0, -1.0), _spot(1, 5.0, -1.0)]
        result = backfill_missing_mz(spots, [_record(0, 3.5, 760.5)])
        self.assertEqual([s["mz"] for s in spots], [-1.0, -1.0])
        self.assertEqual(result, {"missing": 2, "filled": 0})

    def test_valid_mz_is_never_overwritten(self):
        spots = [_spot(0, 3.0, 760.4)]
        result = backfill_missing_mz(spots, [_record(0, 3.0, 999.0)])
        self.assertEqual(spots[0]["mz"], 760.4)
        self.assertEqual(result, {"missing": 0, "filled": 0})

    def test_sibling_arf2_is_the_same_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            eic = Path(tmp) / "AlignmentResult_1.EIC.aef"
            self.assertIsNone(sibling_arf2_path(eic))
            (Path(tmp) / "AlignmentResult_1.arf2").write_bytes(b"x")
            (Path(tmp) / "AlignmentResult_2.arf2").write_bytes(b"x")
            self.assertEqual(
                sibling_arf2_path(eic), Path(tmp) / "AlignmentResult_1.arf2"
            )

    def test_summary_mz_range_ignores_missing_values(self):
        spots = [_spot(0, 3.0, -1.0), _spot(1, 5.0, 885.5), _spot(2, 6.0, 700.0)]
        self.assertEqual(summarize_eic_data(spots)["mz_range"], (700.0, 885.5))
        self.assertIsNone(summarize_eic_data([_spot(0, 3.0, -1.0)])["mz_range"])


if __name__ == "__main__":
    unittest.main()
