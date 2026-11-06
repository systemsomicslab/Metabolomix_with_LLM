"""データセットのハンドル（spec 2026-10-22 §4 の 1）: 開閉・上限・ライブラリ・スレッドへの退避。"""
from __future__ import annotations

import asyncio
import json
import threading
import time

import pytest

import server
from metabolomix.core import handle_dispatch, handles, session_state
from metabolomix.tools import handle_tools
from tests.reanalysis_fixtures import NEG_FEATURES, write_mztab


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def pos_mztab(tmp_path):
    return write_mztab(tmp_path / "pos" / "AlignResult-1.mzTab")


@pytest.fixture
def neg_mztab(tmp_path):
    return write_mztab(tmp_path / "neg" / "AlignResult-2.mzTab", features=NEG_FEATURES,
                       sml_adduct="[M-H]1-")


def _open(path, **kwargs) -> dict:
    kwargs.setdefault("library", "none")
    return json.loads(handle_tools.dataset_handle_open(mztab_path=str(path), **kwargs))


def _write_msp(path, *, ion_mode):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"NAME: fake {ion_mode}\nPRECURSORMZ: 104.107\nIONMODE: {ion_mode}\n"
                    "Num Peaks: 1\n60.08 100\n", encoding="ascii")
    return path


def test_open_two_handles_do_not_touch_global_session(pos_mztab, neg_mztab):
    a = _open(pos_mztab)
    b = _open(neg_mztab)
    assert a["status"] == b["status"] == "success"
    assert a["dataset_handle"] != b["dataset_handle"]
    assert a["polarity"] == "positive" and b["polarity"] == "negative"
    assert a["n_features"] == 3 and a["n_samples"] == 8
    assert session_state.session.dataset is None


def test_same_input_reuses_base_but_not_state(pos_mztab):
    a = _open(pos_mztab)
    b = _open(pos_mztab)
    assert a["reused_base"] is False and b["reused_base"] is True
    ha = handles.registry.get(a["dataset_handle"])
    hb = handles.registry.get(b["dataset_handle"])
    assert ha.dataset is not hb.dataset
    assert ha.dataset.feature_matrix is not hb.dataset.feature_matrix


def test_changed_input_is_reloaded(pos_mztab):
    _open(pos_mztab)
    text = pos_mztab.read_text(encoding="utf-8")
    time.sleep(0.01)
    pos_mztab.write_text(text.replace("500.2", "500.3"), encoding="utf-8")
    assert _open(pos_mztab)["reused_base"] is False


def test_handle_limit(monkeypatch, pos_mztab):
    monkeypatch.setenv(handles.MAX_HANDLES_ENV, "1")
    _open(pos_mztab)
    err = _open(pos_mztab)
    assert err["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert len(err["error"]["details"]["open_handles"]) == 1


def test_close_then_list_and_get(pos_mztab):
    handle_id = _open(pos_mztab)["dataset_handle"]
    listed = json.loads(handle_tools.dataset_handle_list())
    assert listed["n_open"] == 1 and handle_id in listed["handles"]
    closed = json.loads(handle_tools.dataset_handle_close(handle_id))
    assert closed["status"] == "success"
    assert json.loads(handle_tools.dataset_handle_list())["n_open"] == 0
    again = json.loads(handle_tools.dataset_handle_close(handle_id))
    assert again["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"


def test_open_requires_exactly_one_path(pos_mztab):
    out = json.loads(handle_tools.dataset_handle_open())
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=str(pos_mztab),
                                                      job_path=str(pos_mztab)))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"


def test_library_auto_on_mztab_opens_nothing_and_ignores_the_polarity_setting(
        monkeypatch, tmp_path, pos_mztab):
    # mzTab-M は同定に使ったライブラリを持てない（unknown）。極性の設定があっても auto は引かない。
    pos = _write_msp(tmp_path / "lab" / "pos.msp", ion_mode="Positive")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(pos))
    out = _open(pos_mztab, library="auto")
    assert out["status"] == "success" and out["polarity"] == "positive"
    assert out["library"] is None and out["library_record"]["status"] == "unknown"
    assert "library_note" in out
    assert handles.registry.get(out["dataset_handle"]).library.store is None


def test_library_auto_uses_only_the_remembered_library(monkeypatch, tmp_path):
    # 解析フォルダの `*_Loaded.msp2.dbs`（MS-DIAL が同定に使った実体）を開く。極性の設定が
    # 別のファイルを指していても、そちらは引かない。
    from tests.dbs_fixture import record, write_dbs
    from tests.msdial_folder_fixtures import write_msdial_folder

    folder = tmp_path / "gui"
    write_msdial_folder(folder, spots=[{"name": "GABA"}], samples=[
        {"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(3)])
    write_dbs(folder / "Proj_Loaded.msp2.dbs", [("lib", [record()], ["k"])])
    lib = tmp_path / "lab"
    setting_file = _write_msp(lib / "setting_pos.msp", ion_mode="Positive")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(setting_file))
    monkeypatch.setenv("MSDIAL_MSP_NEG", str(setting_file))
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    out = json.loads(handle_tools.dataset_handle_open(directory=str(folder), library="auto"))
    assert out["library"]["file"] == "Proj_Loaded.msp2.dbs"
    assert len(out["library"]["sha256"]) == 64 and out["library"]["record_count"] == 1
    dumped = json.dumps(out)
    assert str(folder) not in dumped and folder.as_posix() not in dumped
    assert "setting_pos" not in dumped
    store = handles.registry.get(out["dataset_handle"]).library.store
    with pytest.raises(Exception):
        store._conn.execute("INSERT INTO meta(key, value) VALUES ('x', '1')")


def test_library_not_loaded_when_polarity_unknown(monkeypatch, tmp_path):
    path = write_mztab(tmp_path / "x.mzTab", sml_adduct="null")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(_write_msp(tmp_path / "pos.msp", ion_mode="Positive")))
    out = _open(path, library="auto")
    assert out["polarity"] is None and out["library"] is None
    assert "極性" in out["library_note"]


# ---------- スレッドへの退避（mcp.call_tool 経由） ----------

def _call(name, args):
    return server.mcp.call_tool(name, args)


def test_install_wraps_handle_tools_only():
    tool = server.mcp._tool_manager._tools["dataset_handle_open"]
    assert tool.is_async
    assert not server.mcp._tool_manager._tools["load_dataset"].is_async
    assert handle_dispatch.install(server.mcp) == []          # 2 回目は何もしない


def test_handle_calls_run_off_loop_and_in_parallel(monkeypatch, pos_mztab, neg_mztab):
    """2 つのハンドルの呼び出しが別スレッドで同時に進む（直列なら a が b を待って詰まる）。"""
    a = _open(pos_mztab)["dataset_handle"]
    b = _open(neg_mztab)["dataset_handle"]
    entered_a, release_a = threading.Event(), threading.Event()
    threads = {}
    original = handles.registry.close

    def slow_close(handle_id):
        threads[handle_id] = threading.current_thread()
        if handle_id == a:
            entered_a.set()
            assert release_a.wait(5), "a と b が直列に回っている"
        else:
            assert entered_a.wait(5), "a と b が直列に回っている"
            release_a.set()
        return original(handle_id)

    monkeypatch.setattr(handles.registry, "close", slow_close)
    loop_thread = []

    async def main():
        loop_thread.append(threading.current_thread())
        return await asyncio.gather(
            _call("dataset_handle_close", {"dataset_handle": a}),
            _call("dataset_handle_close", {"dataset_handle": b}))

    asyncio.run(main())
    assert threads[a] is not loop_thread[0] and threads[b] is not loop_thread[0]
    assert threads[a] is not threads[b]
