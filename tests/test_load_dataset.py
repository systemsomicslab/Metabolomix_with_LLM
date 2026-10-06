"""データフォルダ入口ツール load_dataset と arf 自動選択の検証。

実バイナリ(.arf2/.arf)が無くても通る範囲を確認する:
- PeakProperties.arf を DriftSpots.arf より優先すること
- load_dataset が未登録/欠損を握りつぶさず警告で返すこと
- load_dataset がツールとして登録されていること
"""

import asyncio
import os
import unittest
from pathlib import Path
import tempfile

import server
from metabolomix.core import mcp_core


class ResolveArfPreferenceTests(unittest.TestCase):
    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def _touch(self, directory: Path, name: str) -> None:
        (directory / name).write_bytes(b"")  # 中身は不要（パス解決のみ検証）

    def test_prefers_peakproperties_over_driftspots(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "Sample_DriftSpots.arf")
            self._touch(directory, "Sample_PeakProperties.arf")
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_arf_file_path()
            self.assertIsNotNone(resolved)
            self.assertTrue(resolved.lower().endswith("peakproperties.arf"))

    def test_falls_back_when_no_peakproperties(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "Sample_DriftSpots.arf")
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_arf_file_path()
            self.assertTrue(resolved.lower().endswith("driftspots.arf"))

    def test_explicit_path_is_respected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "Sample_DriftSpots.arf")
            explicit = str(directory / "Sample_DriftSpots.arf")
            self.assertEqual(server.resolve_arf_file_path(explicit), explicit)


class LoadDatasetTests(unittest.TestCase):
    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def test_registered_as_tool(self):
        tools = asyncio.run(server.mcp.list_tools())
        names = {t.name for t in tools}
        self.assertIn("load_dataset", names)

    def test_missing_directory_returns_error(self):
        out = server.load_dataset(directory=str(Path(tempfile.gettempdir()) / "no_such_dir_xyz"))
        self.assertIn("存在しません", out)

    def test_empty_directory_warns_for_both_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = server.load_dataset(directory=tmp)
            self.assertIn(".arf2 ファイルが見つかりませんでした", out)
            self.assertIn(".arf", out)
            # 既定探索先が指定フォルダに更新される
            self.assertEqual(str(mcp_core.DATA_DIR), tmp)


class PickLatestDuplicateTests(unittest.TestCase):
    """同種ファイルが重複（旧版/新版）する場合に最新版を選ぶこと。"""

    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def _touch(self, directory: Path, name: str, mtime: float | None = None) -> Path:
        path = directory / name
        path.write_bytes(b"")
        if mtime is not None:
            os.utime(path, (mtime, mtime))
        return path

    def test_arf2_picks_latest_embedded_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35.arf2")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00.arf2")
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_arf2_file_path()
            self.assertTrue(resolved.endswith("2026_06_01_09_00_00.arf2"))

    def test_arf_embedded_timestamp_beats_mtime(self):
        # 旧タイムスタンプ版に新しい mtime を与えても、埋め込みタイムスタンプが優先される
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(
                directory,
                "AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf",
                mtime=10_000_000_000,  # わざと新しい mtime
            )
            self._touch(
                directory,
                "AlignmentResult_2026_06_01_09_00_00_PeakProperties.arf",
                mtime=1_000_000_000,  # わざと古い mtime
            )
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_arf_file_path()
            self.assertTrue(resolved.endswith("2026_06_01_09_00_00_PeakProperties.arf"))

    def test_falls_back_to_mtime_without_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "older.EIC.aef", mtime=1_000_000_000)
            self._touch(directory, "newer.EIC.aef", mtime=2_000_000_000)
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_eicaef_file_path()
            self.assertTrue(resolved.endswith("newer.EIC.aef"))

    def test_eicaef_resolver_uses_canonical_uppercase_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "newer.eic.aef", mtime=2_000_000_000)
            expected = self._touch(
                directory, "older.EIC.aef", mtime=1_000_000_000
            )
            mcp_core.DATA_DIR = directory

            resolved = server.resolve_eicaef_file_path()

            self.assertEqual(resolved, str(expected.absolute()))


class LoadDatasetBatchAnnounceTests(unittest.TestCase):
    """複数バッチ混在時に load_dataset が選択結果を明示すること。"""

    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def _touch(self, directory: Path, name: str) -> None:
        (directory / name).write_bytes(b"")

    def test_announces_selected_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35.arf2")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00.arf2")
            out = server.load_dataset(directory=str(directory))
            # 検出した最新バッチのタイムスタンプを明示する
            self.assertIn("2026_06_01_09_00_00", out)
            # 旧バッチをスキップした旨が分かる（バッチ数に言及）
            self.assertIn("バッチ", out)

    def test_batch_note_matches_analyzed_arf_not_persample_files(self):
        # A newer per-sample file (.pai2 with a compact 12-digit timestamp) must
        # NOT be announced as the selected batch. The note must name the latest
        # ALIGNMENT (.arf/.arf2) batch that is actually resolved and analyzed.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2024_06_12_18_21_56_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2024_06_12_18_21_56.arf2")
            self._touch(directory, "AlignmentResult_2024_06_13_09_21_45_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2024_06_13_09_21_45.arf2")
            # newer per-sample processing produced only .pai2 (compact ts), no .arf
            self._touch(directory, "sample_x_NEG_202606171407.pai2")
            # project/metadata files carry an AlignmentResult-style timestamp too,
            # but are not analyzable matrices and must be ignored.
            self._touch(directory, "Dataset_2026_06_17_00_34_26.mddata")
            self._touch(directory, "2026_06_17_00_33_16.mdproject")
            note = server._describe_batch_selection(directory)
            self.assertIsNotNone(note)
            self.assertIn("2024_06_13_09_21_45", note)
            self.assertNotIn("202606171407", note)
            self.assertNotIn("2026_06_17", note)

    def test_no_batch_note_for_single_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00.arf2")
            out = server.load_dataset(directory=str(directory))
            self.assertNotIn("複数バッチ", out)


class SelectLatestBatchTests(unittest.TestCase):
    """複数バッチ（処理タイムスタンプ）混在時に最新バッチだけを採用すること。"""

    def test_keeps_only_latest_batch(self):
        paths = [
            "/d/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf",
            "/d/AlignmentResult_2026_05_15_10_13_35.arf2",
            "/d/AlignmentResult_2026_06_01_09_00_00_PeakProperties.arf",
            "/d/AlignmentResult_2026_06_01_09_00_00.arf2",
        ]
        selected = server._select_latest_batch(paths)
        self.assertTrue(all("2026_06_01_09_00_00" in p for p in selected))
        self.assertEqual(len(selected), 2)

    def test_untimestamped_files_are_preserved(self):
        paths = [
            "/d/AlignmentResult_2026_05_15_10_13_35.arf2",
            "/d/AlignmentResult_2026_06_01_09_00_00.arf2",
            "/d/notes.aef",  # タイムスタンプ無し
        ]
        selected = server._select_latest_batch(paths)
        self.assertIn("/d/notes.aef", selected)
        self.assertIn("/d/AlignmentResult_2026_06_01_09_00_00.arf2", selected)
        self.assertNotIn("/d/AlignmentResult_2026_05_15_10_13_35.arf2", selected)

    def test_all_untimestamped_returns_all(self):
        paths = ["/d/a.aef", "/d/b.aef"]
        self.assertEqual(set(server._select_latest_batch(paths)), set(paths))


class ResolvePai2LatestBatchTests(unittest.TestCase):
    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def _touch(self, directory: Path, name: str) -> None:
        (directory / name).write_bytes(b"")

    def test_resolve_pai2_picks_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35.pai2")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00.pai2")
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_pai2_file_path()
            self.assertIsNotNone(resolved)
            self.assertTrue(resolved.endswith("2026_06_01_09_00_00.pai2"))

    def test_resolve_pai2_explicit_path_is_respected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35.pai2")
            explicit = str(directory / "AlignmentResult_2026_05_15_10_13_35.pai2")
            self.assertEqual(server.resolve_pai2_file_path(explicit), explicit)


class ResolveArfCrossBatchTests(unittest.TestCase):
    """旧バッチの PeakProperties.arf があっても最新バッチの方を選ぶこと。"""

    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir

    def _touch(self, directory: Path, name: str) -> None:
        (directory / name).write_bytes(b"")

    def test_arf_picks_peakproperties_from_latest_batch(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2026_05_15_10_13_35_DriftSpots.arf")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00_PeakProperties.arf")
            self._touch(directory, "AlignmentResult_2026_06_01_09_00_00_DriftSpots.arf")
            mcp_core.DATA_DIR = directory
            resolved = server.resolve_arf_file_path()
            self.assertTrue(resolved.endswith("2026_06_01_09_00_00_PeakProperties.arf"))


class ListDataFilesTests(unittest.TestCase):
    """入口ツールの一覧が、解析できないファイルで文脈を埋めないことを縛る。

    MS-DIAL の出力フォルダには測定生データ（.wiff / .wiff.scan / .wiff2 /
    .timeseries.data / .txt）が同居し、実データでは 495 ファイル中 7 割以上が
    それだった。既定でそこまで返すと入口ツールの戻り値だけで 46,977 字になる。
    """

    def setUp(self):
        self._orig_data_dir = mcp_core.DATA_DIR
        self._tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self._tmp.name)
        for name in ("AlignmentResult_2026_01_01_00_00_00_PeakProperties.arf",
                     "AlignmentResult_2026_01_01_00_00_00.arf2",
                     "s1.pai2", "s1.dcl", "study.mddata",
                     "s1.wiff", "s1.wiff.scan", "s1.wiff2",
                     "s1.timeseries.data", "s1.txt", "s1_tags.xml"):
            (self.directory / name).write_bytes(b"")
        mcp_core.DATA_DIR = self.directory

    def tearDown(self):
        mcp_core.DATA_DIR = self._orig_data_dir
        self._tmp.cleanup()

    def test_raw_acquisition_files_are_excluded_by_default(self):
        out = server.list_data_files()
        for parseable in (".arf", ".arf2", ".pai2", ".dcl", ".mddata"):
            self.assertIn(parseable, out)
        for raw in ("s1.wiff", "s1.timeseries.data", "s1.txt"):
            self.assertNotIn(raw, out)

    def test_all_files_opts_back_in(self):
        out = server.list_data_files(all_files=True)
        self.assertIn("s1.wiff", out)

    def test_extension_filter_still_narrows(self):
        out = server.list_data_files(extension=".pai2")
        self.assertIn("s1.pai2", out)
        self.assertNotIn("s1.dcl", out)

    def test_vendor_directories_appear_only_when_opted_in(self):
        """Agilent/Bruker の `.d` は**フォルダ**が 1 検体。

        フォルダだからと一律に除外すると、`.d` だけが入った生データフォルダが
        「空」に見え、入口で「生データが無い」と誤読される。
        """
        (self.directory / "s2.d").mkdir()
        (self.directory / "s2.d" / "AcqData").mkdir()

        self.assertNotIn("s2.d", server.list_data_files())
        self.assertIn("s2.d", server.list_data_files(all_files=True))
        self.assertIn("s2.d", server.list_data_files(extension=".d"))

    def test_non_vendor_directories_never_appear(self):
        """`.d` / `.raw` 以外はフォルダでは計測データにならない（上流の規則）。"""
        (self.directory / "backup.mzml").mkdir()
        self.assertNotIn("backup.mzml", server.list_data_files(all_files=True))

    def test_waters_raw_directory_appears_when_opted_in(self):
        (self.directory / "s3.raw").mkdir()
        self.assertIn("s3.raw", server.list_data_files(all_files=True))

    def test_missing_directory_is_reported_not_returned_as_a_path(self):
        """純関数側は空リストを返し、文面はツール側だけが持つ。

        以前は純関数がエラー文面を1要素のリストで返しており、リゾルバ5箇所が
        「メッセージをパスとして掴まない」よう防御していた。
        """
        from metabolomix.core import path_resolvers

        missing = str(self.directory / "no-such-dir")
        self.assertEqual(path_resolvers.list_data_files(directory=missing), [])
        self.assertIn("存在しません", server.list_data_files(directory=missing))


if __name__ == "__main__":
    unittest.main()
