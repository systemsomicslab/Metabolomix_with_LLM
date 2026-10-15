"""`eic_search`: spot ID・m/z・RT の AND 検索と、ARF2 同定・MS/MS ピーク数の付与。"""
import os
import tempfile
import unittest
from pathlib import Path

from tests.dcl_fixture import build_dcl_bytes
from tests.test_eic_plot import _write_css1


def _trace(rt, intensity):
    return {"file_id": 0, "peak_top": rt, "peak_left": rt - 0.1, "peak_right": rt + 0.1,
            "points": [(rt - 0.1, 1.0), (rt, intensity), (rt + 0.1, 1.0)]}


def _rows(text: str) -> list[dict]:
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("spot_id\t"))
    header = lines[start].split("\t")
    rows = []
    for line in lines[start + 1:]:
        if "\t" not in line:
            break
        rows.append(dict(zip(header, line.split("\t"))))
    return rows


class EicSearchTests(unittest.TestCase):
    def setUp(self):
        from metabolomix.core import mcp_core, session_state
        from metabolomix.eic import tools as tools_eic

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.path = self.tmp / "AlignmentResult_1.EIC.aef"
        (self.tmp / "AlignmentResult_1.arf2").write_bytes(b"placeholder")
        _write_css1(self.path, [
            {"rt": 3.0, "mz": 760.5851, "samples": [_trace(3.0, 500.0)]},
            {"rt": 3.35, "mz": 885.5546, "samples": [_trace(3.35, 9000.0)]},
            {"rt": 6.5, "mz": 500.0, "samples": [_trace(6.5, 40.0)]},
        ])
        self.dcl = self.tmp / "AlignmentResult_1.dcl"
        self.dcl.write_bytes(build_dcl_bytes([
            {"precursor_mz": 760.5851, "rt": 3.0, "spectrum": [(184.07, 999.0), (86.1, 5.0)]},
            {"precursor_mz": 885.5546, "rt": 3.35, "spectrum": []},
            {"precursor_mz": 612.0, "rt": 6.5, "spectrum": [(100.0, 10.0)]},
        ]))
        self.records = [
            {"AlignmentID": 0, "Name": "PC 34:1", "Ontology": "PC",
             "AdductType": "[M+H]+", "RT": 3.0, "MassCenter": 760.5851},
            {"AlignmentID": 1, "Name": "Unknown", "Ontology": "",
             "AdductType": "[M-H]-", "RT": 3.35, "MassCenter": 885.5546},
            {"AlignmentID": 2, "Name": "", "Ontology": "",
             "AdductType": "", "RT": 6.5, "MassCenter": 500.0},
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

    def search(self, **kwargs):
        from metabolomix.eic.tools import eic_search

        return eic_search(file_path=str(self.path), **kwargs)

    def test_no_condition_lists_every_spot_in_id_order_as_tsv(self):
        text = self.search()
        self.assertIn("一致件数: 3", text)
        rows = _rows(text)
        self.assertEqual([r["spot_id"] for r in rows], ["0", "1", "2"])
        self.assertEqual(
            list(rows[0]),
            ["spot_id", "rt", "mz", "max_intensity", "num_samples",
             "name", "ontology", "adduct", "n_msms"],
        )

    def test_spot_ids_pick_exact_spots_and_report_unknown_ids(self):
        text = self.search(spot_ids=[2, 0, 99])
        self.assertEqual([r["spot_id"] for r in _rows(text)], ["0", "2"])
        self.assertIn("99", text)
        self.assertIn("範囲外", text)

    def test_conditions_are_combined_with_and(self):
        rows = _rows(self.search(min_mz=700, max_mz=900, min_rt=3.2, max_rt=4.0))
        self.assertEqual([r["spot_id"] for r in rows], ["1"])
        rows = _rows(self.search(spot_ids=[0, 1], min_mz=800))
        self.assertEqual([r["spot_id"] for r in rows], ["1"])

    def test_sort_by_intensity_replaces_the_old_ranking(self):
        text = self.search(sort="max_intensity", max_results=2)
        self.assertEqual([r["spot_id"] for r in _rows(text)], ["1", "0"])
        self.assertIn("表示上限", text)

    def test_rows_carry_arf2_identity_and_msms_peak_count(self):
        rows = {r["spot_id"]: r for r in _rows(self.search())}
        self.assertEqual(
            (rows["0"]["name"], rows["0"]["ontology"], rows["0"]["adduct"]),
            ("PC 34:1", "PC", "[M+H]+"),
        )
        self.assertEqual(rows["0"]["n_msms"], "2")
        # MS/MS 未取得と precursor 不一致は 0 ではなく「-」
        self.assertEqual(rows["1"]["n_msms"], "-")
        self.assertEqual(rows["2"]["n_msms"], "-")

    def test_without_alignment_dcl_msms_is_unknown_and_said_so(self):
        self.dcl.unlink()
        text = self.search()
        self.assertEqual({r["n_msms"] for r in _rows(text)}, {"-"})
        self.assertIn(".dcl", text)

    def test_invalid_arguments_are_rejected(self):
        self.assertIn("sort", self.search(sort="rt"))
        self.assertIn("min_mz", self.search(min_mz=900, max_mz=700))
        self.assertIn("min_rt", self.search(min_rt=5, max_rt=1))

    def test_nothing_matched_says_so(self):
        text = self.search(min_mz=2000)
        self.assertIn("一致件数: 0", text)
        self.assertEqual(_rows(text + "\nspot_id\tx"), [])


class OldToolsAreGoneTests(unittest.TestCase):
    def test_the_three_search_tools_were_folded_into_eic_search(self):
        from metabolomix.eic import tools as tools_eic

        self.assertIn("eic_search", tools_eic.__all__)
        for name in ("eic_search_by_mz_range", "eic_search_by_rt_range",
                     "eic_rank_by_max_intensity"):
            self.assertNotIn(name, tools_eic.__all__)
            self.assertFalse(hasattr(tools_eic, name), name)


if __name__ == "__main__":
    unittest.main()
