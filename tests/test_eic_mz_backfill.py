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


class MissingMzToolTests(unittest.TestCase):
    """全スポット m/z=-1 の EIC を、隣の ARF2 で補って各ツールが使えること。"""

    def setUp(self):
        from metabolomix.core import mcp_core, session_state
        from metabolomix.eic import tools as tools_eic
        from tests.test_eic_plot import _write_css1

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.path = self.tmp / "alignment.EIC.aef"
        self.arf2 = self.tmp / "alignment.arf2"
        self.arf2.write_bytes(b"placeholder")
        _write_css1(self.path, [
            {"rt": 13.5, "mz": -1.0, "samples": [{
                "file_id": 7, "peak_top": 13.5, "peak_left": 13.4, "peak_right": 13.6,
                "points": [(13.4, 2.0), (13.5, 40.0), (13.6, 3.0)],
            }]},
            {"rt": 22.1, "mz": -1.0, "samples": [{
                "file_id": 7, "peak_top": 22.1, "peak_left": 22.0, "peak_right": 22.2,
                "points": [(22.0, 5.0), (22.1, 90.0), (22.2, 4.0)],
            }]},
        ])
        self.records = [
            {"AlignmentID": 0, "Name": "PC 34:1", "Ontology": "PC",
             "AdductType": "[M+H]+", "RT": 13.5, "MassCenter": 760.5851,
             "HeightAverage": 10.0},
            {"AlignmentID": 1, "Name": "Cer 18:1;O2/24:0", "Ontology": "Cer",
             "AdductType": "[M+H]+", "RT": 22.1, "MassCenter": 650.6,
             "HeightAverage": 20.0},
        ]

        self._saved_data_dir = mcp_core.DATA_DIR
        self._saved_reports = os.environ.get("LIPIDMIX_REPORTS_DIR")
        self._saved_eic = session_state.session.eic
        mcp_core.DATA_DIR = self.tmp
        os.environ["LIPIDMIX_REPORTS_DIR"] = str(self.tmp / "fallback")
        session_state.session.eic = session_state.EicState()

        self._saved_loader = tools_eic.load_arf2_records
        tools_eic.load_arf2_records = lambda _path: self.records

    def tearDown(self):
        from metabolomix.core import mcp_core, session_state
        from metabolomix.eic import tools as tools_eic

        tools_eic.load_arf2_records = self._saved_loader
        mcp_core.DATA_DIR = self._saved_data_dir
        session_state.session.eic = self._saved_eic
        if self._saved_reports is None:
            os.environ.pop("LIPIDMIX_REPORTS_DIR", None)
        else:
            os.environ["LIPIDMIX_REPORTS_DIR"] = self._saved_reports
        self._tmp.cleanup()

    def test_plot_compounds_draws_with_arf2_mz(self):
        from metabolomix.eic.tools import eic_plot_compounds

        payload = json.loads(eic_plot_compounds(
            file_id=7, ontologies=["PC", "Cer"], file_path=str(self.path),
            arf2_path=str(self.arf2), output="payload",
        ))
        self.assertEqual(payload["selection"]["plotted"], 2)
        self.assertEqual(payload["series"][0]["mz"], 760.5851)
        self.assertTrue(any("ARF2" in note for note in payload["caveats"]))

    def test_mz_search_finds_backfilled_spots_and_says_so(self):
        from metabolomix.eic.tools import eic_search_by_mz_range

        text = eic_search_by_mz_range(file_path=str(self.path), min_mz=700, max_mz=800)
        self.assertIn("一致件数: 1", text)
        self.assertIn("mz=760.5851", text)
        self.assertIn("ARF2", text)
        # 2 回目（キャッシュ済み）も注記を落とさない
        again = eic_search_by_mz_range(file_path=str(self.path), min_mz=700, max_mz=800)
        self.assertIn("ARF2", again)

    def test_parser_summary_reports_real_mz_range(self):
        from metabolomix.eic.tools import eic_parser

        text = eic_parser(file_path=str(self.path))
        self.assertIn("650.6", text)
        self.assertNotIn("-1.0", text)
        self.assertIn("ARF2", text)

    def test_single_spot_plot_title_uses_arf2_mz(self):
        from metabolomix.eic.tools import eic_plot_chromatograms

        payload = json.loads(eic_plot_chromatograms(spot_id=0, file_path=str(self.path)))
        self.assertIn("m/z 760.5851", payload["title"])
        self.assertEqual(payload["spot"]["mz"], 760.5851)
        self.assertTrue(any("ARF2" in note for note in payload["caveats"]))

    def test_without_sibling_arf2_the_gap_is_reported_not_silent(self):
        from metabolomix.eic.tools import eic_search_by_mz_range

        self.arf2.unlink()
        text = eic_search_by_mz_range(file_path=str(self.path), min_mz=700, max_mz=800)
        self.assertIn("一致件数: 0", text)
        self.assertIn("m/z", text)
        self.assertIn("2 件", text)


if __name__ == "__main__":
    unittest.main()
