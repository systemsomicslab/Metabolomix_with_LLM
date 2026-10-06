"""dcl_reader（MSDecResult / MS-MS パーサ）の検証。

.dcl は MCP 面から未公開のまま長らくテストが無かった層。MS/MS は同定確度の
最強証拠なので、まずパーサ自身の往復（バイト列→dict）を固定する。
"""
import tempfile
import unittest
from pathlib import Path

from metabolomix.dcl import reader as dcl_reader

from tests.dcl_fixture import build_dcl_bytes


def _write(tmp: str, results, name: str = "sample.dcl") -> str:
    path = Path(tmp) / name
    path.write_bytes(build_dcl_bytes(results))
    return str(path)


SPECTRUM = [(255.2330, 900.0), (281.2486, 1500.0), (283.2643, 400.0)]


class DeserializeDclTests(unittest.TestCase):
    def test_roundtrips_scan_and_spectrum_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [{
                "precursor_mz": 885.5499, "rt": 12.34, "ion_mode": 1,
                "scan_id": 7, "raw_spec_id": 42, "sn": 25.0, "spectrum": SPECTRUM,
            }])
            results = dcl_reader.deserialize_dcl(path)
        self.assertEqual(len(results), 1)
        row = results[0]
        self.assertEqual(row["dcl_index"], 0)
        self.assertEqual(row["scan_id"], 7)
        self.assertEqual(row["raw_spec_id"], 42)
        self.assertAlmostEqual(row["precursor_mz"], 885.5499, places=4)
        self.assertAlmostEqual(row["rt"], 12.34, places=2)
        self.assertEqual(row["n_msms_peaks"], 3)
        self.assertAlmostEqual(row["signal_to_noise"], 25.0, places=3)
        self.assertEqual(len(row["msms_spectrum"]), 3)
        self.assertAlmostEqual(row["msms_spectrum"][1][0], 281.2486, places=4)

    def test_rejects_non_dcl_magic(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.dcl"
            path.write_bytes(b"XX" + b"\x00" * 40)
            with self.assertRaises(ValueError):
                dcl_reader.deserialize_dcl(str(path))

    def test_include_spectrum_false_skips_peaks(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [{"precursor_mz": 885.5, "spectrum": SPECTRUM}])
            results = dcl_reader.deserialize_dcl(path, include_spectrum=False)
        self.assertEqual(results[0]["n_msms_peaks"], 3)  # 元本数は保つ
        self.assertEqual(results[0]["msms_spectrum"], [])

    def test_top_n_peaks_keeps_strongest_in_mz_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [{"precursor_mz": 885.5, "spectrum": SPECTRUM}])
            results = dcl_reader.deserialize_dcl(path, top_n_peaks=2)
        spectrum = results[0]["msms_spectrum"]
        self.assertEqual(len(spectrum), 2)
        self.assertEqual([round(m, 4) for m, _ in spectrum], [255.2330, 281.2486])
        self.assertEqual(results[0]["n_msms_peaks"], 3)  # 元本数は間引きで変わらない

    def test_multiple_results_are_indexed_in_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [
                {"precursor_mz": 100.0, "spectrum": []},
                {"precursor_mz": 200.0, "spectrum": SPECTRUM},
                {"precursor_mz": 300.0, "spectrum": SPECTRUM[:1]},
            ])
            results = dcl_reader.deserialize_dcl(path)
        self.assertEqual([r["dcl_index"] for r in results], [0, 1, 2])
        self.assertEqual([r["n_msms_peaks"] for r in results], [0, 3, 1])


class SummarizeDclTests(unittest.TestCase):
    def test_counts_only_results_carrying_msms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, [
                {"precursor_mz": 100.0, "rt": 1.0, "spectrum": []},
                {"precursor_mz": 200.0, "rt": 2.0, "spectrum": SPECTRUM},
                {"precursor_mz": 300.0, "rt": 3.0, "spectrum": SPECTRUM[:1]},
            ])
            summary = dcl_reader.summarize_dcl(dcl_reader.deserialize_dcl(path))
        self.assertEqual(summary["total_results"], 3)
        self.assertEqual(summary["with_msms"], 2)
        self.assertAlmostEqual(summary["msms_rate_pct"], 66.7, places=1)
        self.assertEqual(summary["msms_peak_count_max"], 3)

    def test_empty_input_reports_error(self):
        self.assertIn("error", dcl_reader.summarize_dcl([]))


class AttachMsmsTests(unittest.TestCase):
    def test_attaches_by_index_when_precursor_matches(self):
        features = [{"m/z": 100.0}, {"m/z": 200.0}]
        results = [
            {"precursor_mz": 100.0, "msms_spectrum": [], "n_msms_peaks": 0},
            {"precursor_mz": 200.0, "msms_spectrum": SPECTRUM, "n_msms_peaks": 3},
        ]
        attached = dcl_reader.attach_msms_to_features(features, results)
        self.assertEqual(attached, 1)
        self.assertEqual(features[1]["n_msms_peaks"], 3)

    def test_skips_when_precursor_disagrees_with_feature_mz(self):
        """索引対応が崩れている疑いがあるときは黙って付けない。"""
        features = [{"m/z": 100.0}]
        results = [{"precursor_mz": 999.0, "msms_spectrum": SPECTRUM, "n_msms_peaks": 3}]
        attached = dcl_reader.attach_msms_to_features(features, results)
        self.assertEqual(attached, 0)
        self.assertNotIn("msms_spectrum", features[0])


if __name__ == "__main__":
    unittest.main()
