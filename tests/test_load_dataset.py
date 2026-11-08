"""入口ツール `load_dataset`（B-2b Task 2）と、ファイル解決・一覧の検証。

`load_dataset` は新しいハンドルを**いつも直近にする**。戻り値は読んだ中身の要約・群と役割の
件数・次の一手で、`.arf2` の概観も前処理前の PCA も返さない。解析の記録の置き場所の基準
`DatasetState.analysis_dir` はここで決まる。fixture はすべてテスト自身が tmp に作る。
"""

import asyncio
import os
import unittest
from pathlib import Path
import tempfile

import server


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


class ListDataFilesTests(unittest.TestCase):
    """入口ツールの一覧が、解析できないファイルで文脈を埋めないことを縛る。

    MS-DIAL の出力フォルダには測定生データ（.wiff / .wiff.scan / .wiff2 /
    .timeseries.data / .txt）が同居し、実データでは 495 ファイル中 7 割以上が
    それだった。既定でそこまで返すと入口ツールの戻り値だけで 46,977 字になる。
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self._tmp.name)
        for name in ("AlignmentResult_2026_01_01_00_00_00_PeakProperties.arf",
                     "AlignmentResult_2026_01_01_00_00_00.arf2",
                     "s1.pai2", "s1.dcl", "study.mddata",
                     "s1.wiff", "s1.wiff.scan", "s1.wiff2",
                     "s1.timeseries.data", "s1.txt", "s1_tags.xml"):
            (self.directory / name).write_bytes(b"")

    def tearDown(self):
        self._tmp.cleanup()

    def _list(self, **kwargs):
        return server.list_data_files(directory=str(self.directory), **kwargs)

    def test_raw_acquisition_files_are_excluded_by_default(self):
        out = self._list()
        for parseable in (".arf", ".arf2", ".pai2", ".dcl", ".mddata"):
            self.assertIn(parseable, out)
        for raw in ("s1.wiff", "s1.timeseries.data", "s1.txt"):
            self.assertNotIn(raw, out)

    def test_all_files_opts_back_in(self):
        out = self._list(all_files=True)
        self.assertIn("s1.wiff", out)

    def test_extension_filter_still_narrows(self):
        out = self._list(extension=".pai2")
        self.assertIn("s1.pai2", out)
        self.assertNotIn("s1.dcl", out)

    def test_vendor_directories_appear_only_when_opted_in(self):
        """Agilent/Bruker の `.d` は**フォルダ**が 1 検体。

        フォルダだからと一律に除外すると、`.d` だけが入った生データフォルダが
        「空」に見え、入口で「生データが無い」と誤読される。
        """
        (self.directory / "s2.d").mkdir()
        (self.directory / "s2.d" / "AcqData").mkdir()

        self.assertNotIn("s2.d", self._list())
        self.assertIn("s2.d", self._list(all_files=True))
        self.assertIn("s2.d", self._list(extension=".d"))

    def test_non_vendor_directories_never_appear(self):
        """`.d` / `.raw` 以外はフォルダでは計測データにならない（上流の規則）。"""
        (self.directory / "backup.mzml").mkdir()
        self.assertNotIn("backup.mzml", self._list(all_files=True))

    def test_waters_raw_directory_appears_when_opted_in(self):
        (self.directory / "s3.raw").mkdir()
        self.assertIn("s3.raw", self._list(all_files=True))

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


# ---------- load_dataset ----------

import json

import pytest

from metabolomix.core import handle_dispatch, handles, session_state
from metabolomix.tools import handle_tools
from tests.dbs_fixture import record, write_dbs
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.reanalysis_fixtures import write_mztab
from tests.config_fixtures import use_config

_SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
            for i in range(4)]
_SPOTS = [{"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0]},
          {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]}]


@pytest.fixture
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


def _folder(path):
    write_msdial_folder(path, spots=_SPOTS, samples=_SAMPLES)
    return path


def _load(**kwargs) -> dict:
    out = json.loads(handle_tools.load_dataset(**kwargs))
    assert out["status"] == "success", out
    return out


def test_load_dataset_makes_current_handle(tmp_path, fresh_state):
    first = _load(directory=str(_folder(tmp_path / "a")))
    assert session_state.session.current_handle_id == first["dataset_handle"]
    second = _load(directory=str(_folder(tmp_path / "b")))
    # 何度開いても、開いたハンドルがいつも直近になる（前のは閉じずに残る）。
    assert session_state.session.current_handle_id == second["dataset_handle"]
    assert {h.handle_id for h in handles.registry.list()} == {
        first["dataset_handle"], second["dataset_handle"]}
    handle = handles.registry.peek(second["dataset_handle"])
    assert handle.dataset is session_state.session.dataset


def test_load_dataset_opens_remembered_library(tmp_path, fresh_state, monkeypatch):
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    folder = _folder(tmp_path / "gui")
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record(mz=500.0, ion_mode=0)], ["k"])])
    out = _load(directory=str(folder))
    assert out["library"]["file"] == "Proj_Loaded.msp2.dbs" and out["library"]["record_count"] == 1
    assert str(folder) not in json.dumps(out, ensure_ascii=False)
    assert handles.registry.get(out["dataset_handle"]).library.store is not None


def test_load_dataset_returns_group_and_role_counts(tmp_path, fresh_state):
    out = _load(directory=str(_folder(tmp_path / "a")))
    counts = out["sample_counts"]
    assert counts["groups"] == {"A": 2, "B": 2}
    assert sum(counts["roles"].values()) == 4
    assert counts["n_ungrouped"] == 0
    assert len(out["next"]) == 3
    # 共有の部品（`dataset_status` の要約がこれを使う）が数えたものと同じ。
    from metabolomix.tools import mztab_tools
    assert mztab_tools.sample_group_role_counts(session_state.session.dataset) == counts


def test_load_dataset_has_no_initial_pca(tmp_path, fresh_state):
    """読み込みの後の概観（`.arf2` の要約・前処理前の PCA）は返さず、PCA の状態も作らない。"""
    out = _load(directory=str(_folder(tmp_path / "a")))
    text = json.dumps(out, ensure_ascii=False)
    assert "PCA" not in text.replace("dataset_pca", "")
    assert "概観" not in text
    ds = session_state.session.dataset
    assert ds.last_pca is None and ds.pp_matrix is None


def test_load_dataset_is_always_threaded(tmp_path, fresh_state):
    assert handle_dispatch.ALWAYS_THREADED == frozenset({"load_dataset", "dataset_merge_polarities"})
    assert server.mcp._tool_manager._tools["load_dataset"].is_async


def test_old_entry_names_are_gone():
    names = set(server.mcp._tool_manager._tools)
    assert {"load_dataset", "dataset_list", "dataset_close"} <= names
    assert not names & {"dataset_load", "dataset_handle_open", "dataset_handle_list",
                        "dataset_handle_close"}


def test_analysis_dir_per_adapter(tmp_path, fresh_state):
    # 解析フォルダ → 選んだ .arf のフォルダ（ファイルを渡しても同じ）
    folder = _folder(tmp_path / "gui")
    _load(directory=str(folder))
    assert session_state.session.dataset.analysis_dir == str(folder.resolve())
    arf = next(folder.glob("*_PeakProperties.arf"))
    _load(directory=str(arf))
    assert Path(session_state.session.dataset.analysis_dir).resolve() == folder.resolve()
    # mzTab-M → mzTab-M のフォルダ
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab")
    _load(mztab_path=str(mztab))
    assert Path(session_state.session.dataset.analysis_dir).resolve() == mztab.parent.resolve()


def test_analysis_dir_is_not_part_of_the_fingerprint(tmp_path, fresh_state):
    from metabolomix.analysis.result_state import dataset_fingerprint

    _load(directory=str(_folder(tmp_path / "a")))
    ds = session_state.session.dataset
    before = dataset_fingerprint(ds)
    ds.analysis_dir = str(tmp_path / "elsewhere")
    assert dataset_fingerprint(ds) == before


def test_analysis_dir_for_job_and_pipeline_is_the_chosen_arf_folder(tmp_path, fresh_state):
    """job / pipeline の解析フォルダ経路は、mzTab-M との照合に通った .arf のフォルダ。
    `source="mztab"` で読めば mzTab-M のフォルダ。"""
    from tests.test_dataset_opening import _job, _pipeline

    job = _job(tmp_path)
    _load(job_path=str(job))
    ds = session_state.session.dataset
    assert Path(ds.analysis_dir).resolve() == (job.parent / "msdial").resolve()
    _load(job_path=str(job), source="mztab")
    assert Path(session_state.session.dataset.analysis_dir).resolve() == (
        job.parent / "msdial").resolve()

    run = _pipeline(tmp_path, job)
    _load(pipeline_path=str(run))
    assert Path(session_state.session.dataset.analysis_dir).resolve() == (
        job.parent / "msdial").resolve()


def test_library_loaded_before_load_dataset_returns_only_the_file_name(tmp_path, fresh_state,
                                                                      monkeypatch):
    """`library_load`（ハンドルなし）が先に読んだライブラリがハンドルへ移っても、置き場所は戻り値に出ない。"""
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    lab = tmp_path / "lab_secret_dir"
    lab.mkdir()
    msp = lab / "pos.msp"
    msp.write_text("NAME: Choline\nPRECURSORMZ: 104.107\nIONMODE: Positive\nNum Peaks: 2\n"
                   "60.0808 100\n104.107 50\n", encoding="ascii")
    assert json.loads(library_tools.library_load(file_path=str(msp)))["status"] == "success"
    raw = handle_tools.load_dataset(mztab_path=str(write_mztab(tmp_path / "run" / "A-1.mzTab")))
    out = json.loads(raw)
    assert out["library"]["file"] == "pos.msp"
    assert "lab_secret_dir" not in raw  # JSON の中ではバックスラッシュが重ねて書かれるので、区切りに依らない名前で見る
    assert str(lab) not in raw and str(lab).replace(chr(92), chr(92) * 2) not in raw
    assert "lab_secret_dir" not in handle_tools.dataset_list()


def test_bare_library_load_refuses_when_the_recorded_file_changed(tmp_path, fresh_state,
                                                                  monkeypatch):
    """`library="auto"` で開いたあとに記憶のファイルが変わっても、引数なしの library_load は
    不一致を理由に拒む（開いた store はそのまま残る）。"""
    from metabolomix.library import tools as library_tools

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    folder = _folder(tmp_path / "gui")
    dbs = write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record(mz=500.0, ion_mode=0)], ["k"])])
    out = _load(directory=str(folder))
    assert out["library"]["file"] == "Proj_Loaded.msp2.dbs"
    write_dbs(dbs, [("lib", [record(mz=501.0, ion_mode=0)], ["k"])])
    res = json.loads(library_tools.library_load())
    assert res["status"] == "error" and "不一致" in res["message"]


class ListDataFilesEntryTraceTests(unittest.TestCase):
    """ENTRY POINT が「MS-DIAL を回した形跡」と呼ぶ mzTab-M と analysis-job.json も既定の一覧に出る。"""

    def test_default_listing_shows_mztab_and_job_record(self):
        from metabolomix.core import path_resolvers

        with tempfile.TemporaryDirectory() as tmp:
            for name in ("Height_AlignmentResult_1.mzTab", "analysis-job.json", "s1.wiff"):
                (Path(tmp) / name).write_bytes(b"")
            listed = [Path(p).name for p in path_resolvers.list_data_files(directory=tmp)]
        self.assertIn("Height_AlignmentResult_1.mzTab", listed)
        self.assertIn("analysis-job.json", listed)
        self.assertNotIn("s1.wiff", listed)
