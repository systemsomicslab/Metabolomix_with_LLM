"""`.dcl` の解決（最新バッチの選択）の検証。

`dcl_parser` と、MS/MS を PAI2 の同定検証へ載せる経路（`pai2_parser` / `verify_peak_annotation`）は
B-2b で消えた。`.dcl` を読む `dcl_find_msms` / `dcl_spot_spectrum` / `dcl_diagnostic_ions` の
検証は `test_library_detail.py` `test_reanalysis_tools.py` が持つ。
"""
import tempfile
import unittest
from pathlib import Path

from metabolomix.core import mcp_core

from tests.dcl_fixture import build_dcl_bytes

SPECTRUM = [(255.2330, 900.0), (281.2486, 1500.0), (283.2643, 400.0)]


def _write_dcl(directory: Path, name: str, results) -> Path:
    path = directory / name
    path.write_bytes(build_dcl_bytes(results))
    return path


class DclResolverTests(unittest.TestCase):
    def test_resolver_picks_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            _write_dcl(directory, "AlignmentResult_2024_01_01_00_00_00.dcl",
                       [{"precursor_mz": 1.0, "spectrum": []}])
            newer = _write_dcl(directory, "AlignmentResult_2026_05_15_10_13_35.dcl",
                               [{"precursor_mz": 2.0, "spectrum": SPECTRUM}])
            original = mcp_core.DATA_DIR
            mcp_core.DATA_DIR = directory
            self.addCleanup(setattr, mcp_core, "DATA_DIR", original)
            from metabolomix.core.path_resolvers import resolve_dcl_file_path
            self.assertEqual(Path(resolve_dcl_file_path()), newer)


if __name__ == "__main__":
    unittest.main()
