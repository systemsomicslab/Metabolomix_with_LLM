"""キュレーションは DatasetState（ハンドル）の上で動く（B-2b Task 4）。

`.arf2` はデータセットのものを使い、`data/` は探さない。fixture は tmp に作る。
"""
import inspect
import json
import threading

import pytest

from metabolomix.core import handles, session_state
from metabolomix.library import store as library_store
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import curation_tools, handle_tools
from tests.config_fixtures import use_config
from tests.curation_fixtures import write_alignment_set, write_suggest_set
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata
from tests.engine_chain import engine_compare

CURATION_TOOLS = ("curation_review", "curation_suggest", "curation_submit", "curation_flags",
                  "curation_view_data")


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
    """4 試料（g1 / g2 が 2 つずつ）の解析フォルダ。"""
    paths = write_alignment_set(root / name, n_files=4)
    write_mddata(paths["arf2"].with_suffix(".mddata"), [
        analysis_file_row(i, f"s{i}", "g1" if i < 2 else "g2", file_type=0, order=None, batch=1)
        for i in range(4)])
    return paths


def _open(paths, *, library=True):
    body = json.loads(handle_tools.load_dataset(directory=str(paths["arf2"].parent)))
    assert body["status"] == "success", body
    handle = handles.registry.get(body["dataset_handle"])
    if library:
        handle.library.store = library_store.open_store(paths["msp"])
    return handle


def test_curation_tools_take_a_dataset_handle_and_no_file_path():
    for name in CURATION_TOOLS:
        params = inspect.signature(getattr(curation_tools, name)).parameters
        assert "dataset_handle" in params, name
        assert "file_path" not in params, name


def test_curation_review_reads_dataset_arf2(tmp_path):
    first = _folder(tmp_path, "first")
    second = _folder(tmp_path, "second")
    h1 = _open(first)
    h2 = _open(second)                                  # 直近は h2
    body = json.loads(curation_tools.curation_review(dataset_handle=h1.handle_id))
    assert body["counts"]["ok"] == 1, body
    assert (first["arf2"].parent / "curation").is_dir()
    assert not (second["arf2"].parent / "curation").exists()    # 渡したハンドルの .arf2 だけ
    # ハンドルを省けば直近（second）。ビューアのデータは同じ review_id から引ける
    body2 = json.loads(curation_tools.curation_review())
    assert (second["arf2"].parent / "curation").is_dir()
    paged = json.loads(curation_tools.curation_view_data(body2["review_id"], 0))
    assert paged["review"]["review_id"] == body2["review_id"]
    assert h1.handle_id != h2.handle_id


def test_curation_does_not_search_data_dir_without_a_dataset():
    for call in (curation_tools.curation_review, curation_tools.curation_suggest,
                 curation_tools.curation_flags):
        err = json.loads(call())["error"]
        assert err["code"] == "missing_state" and err["state"] == "dataset"
        assert err["required_tools"] == ["load_dataset"]
    err = json.loads(curation_tools.curation_submit(
        review_id="cr-20260101-000000-abcd", flags=[{"spot_id": 0, "flag": "wrong"}]))["error"]
    assert err["state"] == "dataset"


def test_curation_after_close_does_not_fall_back(tmp_path):
    handle = _open(_folder(tmp_path))
    handles.registry.close(handle.handle_id)
    assert json.loads(curation_tools.curation_flags())["error"]["state"] == "dataset"
    body = json.loads(curation_tools.curation_flags(dataset_handle=handle.handle_id))
    assert "DATASET_HANDLE_NOT_FOUND" in json.dumps(body)


def test_curation_without_arf2_is_explained():
    session_state.session.dataset = DatasetState()        # mzTab-M だけ（.arf2 を持たない）データセット
    for call in (curation_tools.curation_review, curation_tools.curation_suggest,
                 curation_tools.curation_flags):
        body = json.loads(call())
        assert body["status"] == "error" and "arf2" in body["message"]
        assert "load_dataset" in body["message"]


def test_curation_uses_dataset_library(tmp_path):
    with_library = _open(_folder(tmp_path, "a"))
    without = _open(_folder(tmp_path, "b"), library=False)         # 直近。ライブラリが無い
    assert session_state.session.library.store is None
    ok = json.loads(curation_tools.curation_review(dataset_handle=with_library.handle_id))
    assert ok["counts"]["ok"] == 1, ok
    missing = json.loads(curation_tools.curation_review(dataset_handle=without.handle_id))
    assert missing["error"]["required_tools"] == ["library_load"]
    assert json.loads(curation_tools.curation_review())["error"]["required_tools"] == ["library_load"]


def test_submit_refreshes_flags_on_open_handles(tmp_path):
    from metabolomix.tools import dataset_analysis_tools as dat

    paths = _folder(tmp_path)
    first = _open(paths)
    second = _open(paths)
    other = _open(write_suggest_set(tmp_path / "elsewhere"))        # 別の .arf2（sha256 が違う）
    assert first.dataset.curation_flags["wrong"] == set()
    review_id = json.loads(curation_tools.curation_review(dataset_handle=first.handle_id))["review_id"]
    body = json.loads(curation_tools.curation_submit(
        review_id=review_id, dataset_handle=first.handle_id,
        flags=[{"spot_id": 1, "flag": "wrong", "note": ""}]))
    assert body["status"] == "ok", body
    assert sorted(body["open_handles"]["refreshed"]) == sorted([first.handle_id, second.handle_id])
    assert first.dataset.curation_flags["wrong"] == {1}
    assert second.dataset.curation_flags["wrong"] == {1}
    assert other.dataset.curation_flags["wrong"] == set()          # 別の .arf2 は触らない
    # 書き出しは出したばかりの判断を使う（同じ .arf2 を持つ別のハンドルでも）
    for handle in (first, second):
        engine_compare("g1", "g2", dataset_handle=handle.handle_id)
        out = json.loads(dat.dataset_export_differential(
            str(tmp_path / f"{handle.handle_id}.tsv"), apply_curation=True,
            dataset_handle=handle.handle_id))
        assert out["curation"]["state"] == "applied", out
        assert out["curation"]["wrong_excluded"] == 1


def test_submit_reports_a_handle_it_could_not_refresh(tmp_path, monkeypatch):
    paths = _folder(tmp_path)
    first = _open(paths)
    second = _open(paths)
    monkeypatch.setattr(curation_tools, "REFRESH_LOCK_TIMEOUT", 0.2)
    review_id = json.loads(curation_tools.curation_review(dataset_handle=first.handle_id))["review_id"]
    held, release = threading.Event(), threading.Event()

    def hold():
        with second.lock:
            held.set()
            release.wait(5)

    thread = threading.Thread(target=hold)
    thread.start()
    held.wait(5)
    try:
        body = json.loads(curation_tools.curation_submit(
            review_id=review_id, dataset_handle=first.handle_id,
            flags=[{"spot_id": 1, "flag": "wrong"}]))
    finally:
        release.set()
        thread.join()
    assert body["status"] == "ok"
    assert body["open_handles"]["refreshed"] == [first.handle_id]
    assert body["open_handles"]["not_refreshed"] == [second.handle_id]
