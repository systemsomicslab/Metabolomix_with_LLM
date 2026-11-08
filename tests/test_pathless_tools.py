"""パスを省いたツールは直近のデータセットのファイルを使い、`data/` を探さない（B-2b Task 9）。

fixture は tmp に作った解析フォルダで、`load_dataset(directory=...)` で開く。
"""
import inspect
import json
import shutil

import pytest

from metabolomix.core import handles, session_state
from metabolomix.dcl import tools as dcl_tools
from metabolomix.library import tools as library_tools
from metabolomix.tools import curation_tools, handle_tools, reanalysis_tools, reports
from metabolomix.tools import dataset as dataset_tools
from tests.config_fixtures import use_config
from tests.curation_fixtures import write_alignment_set
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _folder(root, name="neg"):
    paths = write_alignment_set(root / name, n_files=4)
    write_mddata(paths["arf2"].with_suffix(".mddata"), [
        analysis_file_row(i, f"s{i}", "g1" if i < 2 else "g2", file_type=0, order=None, batch=1)
        for i in range(4)])
    return paths


def _open(paths, *, library="auto"):
    body = json.loads(handle_tools.load_dataset(directory=str(paths["arf2"].parent),
                                                library=library))
    assert body["status"] == "success", body
    return handles.registry.get(body["dataset_handle"])


def _with_library(handle, paths):
    from metabolomix.library import store as library_store

    handle.library.store = library_store.open_store(paths["msp"])


def test_pathless_tools_use_current_dataset_files(tmp_path):
    paths = _folder(tmp_path)
    handle = _open(paths)
    _with_library(handle, paths)

    spectrum = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=0))
    assert spectrum["status"] == "ok" and spectrum["n_peaks"] == 2

    found = json.loads(dcl_tools.dcl_find_msms(precursor_mz=760.5851))
    assert found["status"] == "success"
    assert found["file"] == paths["arf2"].with_suffix(".dcl").name

    matched = json.loads(library_tools.library_match_feature(spot_id=0))
    assert matched["status"] == "success", matched
    assert "PC 34:1" in matched["candidates_table"]


def test_pathless_tools_do_not_fall_back_after_close(tmp_path, monkeypatch):
    paths = _folder(tmp_path)
    handle = _open(paths)
    # 閉じた後に落ちうる先: 昔の探索先だった data/ に同名のファイルを置く。
    data = tmp_path / "data"
    data.mkdir()
    for suffix in (".arf2", ".dcl"):
        shutil.copy(paths["arf2"].with_suffix(suffix), data / f"AlignmentResult_x{suffix}")
    monkeypatch.setenv("LIPIDMIX_DATA_DIR", str(data))
    monkeypatch.chdir(tmp_path)
    # もう 1 本、別のハンドルを開いておく（それも直近の代わりにしない）。
    other = _folder(tmp_path, "other")
    other_handle = _open(other)
    handles.registry.close(handle.handle_id)
    session_state.session.current_handle_id = handle.handle_id  # 閉じた直近を指したまま

    def code(raw):
        return json.loads(raw)["error"]

    for raw in (dcl_tools.dcl_spot_spectrum(spot_id=0),
                dcl_tools.dcl_find_msms(precursor_mz=760.5851),
                dcl_tools.dcl_diagnostic_ions(spot_ids=[0], ions=[{"mz": 184.0733}]),
                curation_tools.curation_review(),
                reports.write_report("a-1", "DS", "本文")):
        error = code(raw)
        assert error["code"] == "missing_state" and error["state"] == "dataset", raw
        assert error["required_tools"] == ["load_dataset"]
    # 照合の 3 本（と対向プロット）は、閉じたハンドルのライブラリ（記憶）や照合結果を読むので、データセット・
    # ライブラリ・照合結果のどれかが無い missing_state で止まる。別のハンドルのライブラリや `data/` へは落ちない。
    for raw in (library_tools.library_match_feature(precursor_mz=760.5851),
                reanalysis_tools.library_msms_qualify(rules={}, spot_id=0),
                library_tools.library_plot_mirror()):
        error = code(raw)
        assert error["code"] == "missing_state" and error["state"] in ("dataset", "library_store", "library_match"), raw
        assert set(error["required_tools"]) <= {"load_dataset", "library_load", "library_match_feature"}, raw
    assert not (data / "reports").exists() and not (tmp_path / "reports").exists()
    assert not (other["arf2"].parent / "reports").exists()  # 別のハンドルの解析フォルダにも書いていない
    assert other_handle.dataset.analysis_dir == str(other["arf2"].parent)


def test_list_data_files_requires_directory(tmp_path):
    param = inspect.signature(dataset_tools.list_data_files).parameters["directory"]
    assert param.default is inspect.Parameter.empty
    from metabolomix.core import path_resolvers

    param = inspect.signature(path_resolvers.list_data_files).parameters["directory"]
    assert param.default is inspect.Parameter.empty
    (tmp_path / "a.arf2").write_bytes(b"x")
    assert str(tmp_path / "a.arf2") in dataset_tools.list_data_files(directory=str(tmp_path))


def test_reports_go_to_analysis_dir(tmp_path, monkeypatch):
    paths = _folder(tmp_path)
    handle = _open(paths)
    reports_dir = paths["arf2"].parent / "reports"
    saved = reports.write_report("a-1", "DS", "本文")
    assert "レポートを保存" in saved
    assert (reports_dir / "a-1.md").is_file()
    assert "本文" in reports.read_report("a-1")
    assert handle.dataset.analysis_dir == str(paths["arf2"].parent)

    # `[paths] reports_dir` があればそちらが先。
    elsewhere = tmp_path / "elsewhere"
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache",
                                       "paths.reports_dir": elsewhere})
    reports.write_report("a-2", "DS", "別の本文")
    assert (elsewhere / "a-2.md").is_file()
    assert not (reports_dir / "a-2.md").exists()


def test_library_load_without_memory_is_refused(tmp_path, monkeypatch):
    paths = _folder(tmp_path)
    msp = tmp_path / "pos.msp"
    msp.write_text(paths["msp"].read_text(encoding="utf-8"), encoding="utf-8")
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache",
                                       "library.msp_positive": msp})
    handle = _open(paths)
    assert handle.dataset.library_record["status"] == "unknown"
    body = json.loads(library_tools.library_load())
    assert body["status"] == "error" and body["code"] == "LIBRARY_NOT_REMEMBERED"
    assert "file_path" in body["message"]
    assert str(tmp_path) not in json.dumps(body, ensure_ascii=False)  # 名前だけで、置き場所は出さない
    assert handle.library.store is None  # 設定の極性の MSP を使っていない
    # 明示すれば開く。
    loaded = json.loads(library_tools.library_load(file_path=str(paths["msp"])))
    assert loaded["status"] == "success"


def test_library_load_has_no_ion_mode():
    params = inspect.signature(library_tools.library_load).parameters
    assert "ion_mode" not in params
    from metabolomix.core import path_resolvers

    assert "ion_mode" not in inspect.signature(path_resolvers.resolve_library_path).parameters
