"""対話のセッションとハンドルの一本化（B-1 Task 9 / spec 2026-11-06 §5）。

読み込んだデータセットはすべてハンドルの登録簿（`handles.registry`）に入り、対話のセッションは
直近のハンドル（`current_handle_id`）を指すだけになる。`dataset_handle` を省いた呼び出しも
直近のハンドルのロックの下で動き、閉じた直近のハンドルから別のハンドルへ黙って落ちない。

登録簿はテストごとに conftest の `_fresh_handle_registry` が新しくする。
fixture はすべてこのテストが tmp に作る。
"""
from __future__ import annotations

import asyncio
import json
import threading
import time

import pytest

import server
from metabolomix.core import handles, session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import handle_tools, mztab_tools
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.reanalysis_fixtures import NEG_FEATURES, write_mztab


@pytest.fixture(autouse=True)
def fresh_session():
    session_state.session = session_state.AnalysisSession()
    yield
    for handle in handles.registry.list():   # ライブラリの接続を閉じる（登録簿は conftest が戻す）
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


def _folder(path):
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(4)]
    spots = [{"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0]},
             {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]}]
    write_msdial_folder(path, spots=spots, samples=samples)
    return path


def _open(**kwargs) -> str:
    out = json.loads(handle_tools.dataset_handle_open(library="none", **kwargs))
    assert out["status"] == "success", out
    return out["dataset_handle"]


def _current():
    return session_state.session.current_handle_id


def _list_rows() -> list[dict]:
    table = json.loads(handle_tools.dataset_handle_list())["handles"]
    header, *lines = table.splitlines()
    columns = header.split("\t")
    return [dict(zip(columns, line.split("\t"))) for line in lines]


def _text(result) -> str:
    """`mcp.call_tool` の戻り（content の列）から本文を取り出す。"""
    return "".join(block.text for block in result)


def test_dataset_load_registers_handle_and_sets_current(tmp_path):
    out = mztab_tools.dataset_load(directory=str(_folder(tmp_path / "a")))
    hid = _current()
    assert hid and hid in out
    assert handles.registry.peek(hid).dataset is session_state.session.dataset


def test_closing_current_handle_does_not_fall_back(tmp_path):  # Review Focus 3
    mztab_tools.dataset_load(directory=str(_folder(tmp_path / "a")))
    first = _current()
    other = _open(directory=str(_folder(tmp_path / "b")))
    assert json.loads(handle_tools.dataset_handle_close(first))["status"] == "success"
    assert session_state.session.dataset is None
    assert '"missing_state"' in mztab_tools.dataset_status()
    # mcp 経由（包みを通る）でも、開いている別のハンドルへは落ちない。
    via_mcp = _text(asyncio.run(server.mcp.call_tool("dataset_status", {})))
    assert '"missing_state"' in via_mcp
    assert handles.registry.peek(other) is not None
    assert _current() is None


def test_library_follows_current_handle(tmp_path):
    sess = session_state.session
    own = sess.library
    assert isinstance(own, session_state.LibraryState)
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    handle = handles.registry.peek(_current())
    assert sess.library is handle.library and sess.library is not own
    replacement = session_state.LibraryState()
    sess.library = replacement          # 直近のハンドルがあれば、そのハンドルへ
    assert handle.library is replacement
    sess.dataset = None                 # 直近のハンドルが無ければ、セッションの LibraryState
    assert sess.library is own
    standalone = session_state.LibraryState()
    sess.library = standalone
    assert sess.library is standalone and handle.library is replacement


def test_handle_limit_refuses_instead_of_evicting(monkeypatch, tmp_path):
    monkeypatch.setenv(handles.MAX_HANDLES_ENV, "2")
    paths = [write_mztab(tmp_path / str(i) / "AlignResult-1.mzTab") for i in range(3)]
    mztab_tools.dataset_load(mztab_path=str(paths[0]))
    first = _current()
    mztab_tools.dataset_load(mztab_path=str(paths[1]))
    second = _current()
    err = json.loads(mztab_tools.dataset_load(mztab_path=str(paths[2])))["error"]
    assert err["code"] == "DATASET_HANDLE_LIMIT"
    assert {h.handle_id for h in handles.registry.list()} == {first, second}
    # 失敗した読み込みは直近のデータセットを動かさない。
    assert _current() == second
    assert session_state.session.dataset is handles.registry.peek(second).dataset


def test_interactive_call_waits_for_agent_call_on_same_handle(tmp_path):
    """省略の呼び出しと、同じハンドルを明示した呼び出しが重ならない（同じロックを取る）。"""
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    hid = _current()
    handle = handles.registry.get(hid)
    real_lock = handle.lock
    spans = []

    class SpyLock:
        def __enter__(self):
            real_lock.acquire()
            start = time.perf_counter()
            time.sleep(0.2)
            spans.append(start)
            return self

        def __exit__(self, *exc):
            spans.append(time.perf_counter())
            real_lock.release()

    handle.lock = SpyLock()

    async def main():
        return await asyncio.gather(
            server.mcp.call_tool("dataset_status", {"dataset_handle": hid}),
            server.mcp.call_tool("dataset_status", {}))

    agent, interactive = asyncio.run(main())
    assert len(spans) == 4, "省略の呼び出しが直近のハンドルのロックを取っていない"
    start1, end1, start2, end2 = spans
    assert end1 <= start2
    assert '"n_features":3' in _text(agent) and '"n_features":3' in _text(interactive)


def test_interactive_call_with_current_handle_runs_off_event_loop(monkeypatch, tmp_path):
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    seen = []
    original = mztab_tools.json_payload

    def spy(*args, **kwargs):
        seen.append(threading.current_thread())
        return original(*args, **kwargs)

    monkeypatch.setattr(mztab_tools, "json_payload", spy)

    async def main():
        seen.append(threading.current_thread())
        await server.mcp.call_tool("dataset_status", {})

    asyncio.run(main())
    assert len(seen) == 2 and seen[0] is not seen[1]


def test_interactive_call_keeps_the_handle_it_locked(tmp_path):
    """呼び出しの途中で直近のデータセットが変わっても、ロックを取ったハンドルのデータを読む
    （ロックを取っていない別のハンドルの状態を触らない）。"""
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    handle = handles.registry.get(_current())
    other = DatasetState()
    other.source_format = "other"
    real_lock = handle.lock
    entered, release = threading.Event(), threading.Event()

    class GateLock:
        def __enter__(self):
            real_lock.acquire()
            entered.set()
            assert release.wait(5)
            return self

        def __exit__(self, *exc):
            real_lock.release()

    handle.lock = GateLock()

    async def main():
        task = asyncio.create_task(server.mcp.call_tool("dataset_status", {}))
        assert await asyncio.to_thread(entered.wait, 5)
        session_state.session.dataset = other   # dataset_load と同じく、イベントループ上で切り替える
        release.set()
        return await task

    out = json.loads(_text(asyncio.run(main())))
    assert out["source_format"] == "mztab" and out["n_features"] == 3
    assert session_state.session.dataset is other


def test_library_loaded_by_interactive_call_is_usable_from_another_thread(tmp_path, monkeypatch):
    """省略の `library_load` はワーカースレッドで開く。その store をイベントループ側の
    ツール（`curation_*` / `verify_peak_annotation`）や別のワーカーが読めること
    （sqlite3 の接続は既定では開いたスレッドでしか使えない）。"""
    from metabolomix.library import tools as library_tools

    msp = tmp_path / "lab" / "pos.msp"
    msp.parent.mkdir()
    msp.write_text("NAME: Choline\nPRECURSORMZ: 104.107\nIONMODE: Positive\nNum Peaks: 2\n"
                   "60.0808 100\n104.107 50\n", encoding="ascii")
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    loaded = json.loads(_text(asyncio.run(
        server.mcp.call_tool("library_load", {"file_path": str(msp)}))))
    assert loaded["status"] == "success", loaded
    store = session_state.session.library.store
    assert store is handles.registry.peek(_current()).library.store
    # このスレッド（開いたワーカーではない）から接続を使う。
    records, _ = store.search(name="Choline")
    assert [row["name"] for row in records] == ["Choline"]


def test_handle_list_marks_current(tmp_path):
    mztab_tools.dataset_load(mztab_path=str(write_mztab(tmp_path / "a" / "AlignResult-1.mzTab")))
    current = _current()
    other = _open(mztab_path=str(write_mztab(tmp_path / "b" / "AlignResult-2.mzTab",
                                             features=NEG_FEATURES, sml_adduct="[M-H]1-")))
    rows = {row["dataset_handle"]: row for row in _list_rows()}
    assert rows[current]["current"] == "True" and rows[other]["current"] == "False"
    assert rows[current]["source_kind"] == "session"


def test_reassigning_the_current_dataset_does_not_open_a_handle():  # R13
    sess = session_state.session
    ds = DatasetState()
    sess.dataset = ds
    hid = sess.current_handle_id
    sess.dataset = ds
    sess.dataset = sess.dataset
    assert sess.current_handle_id == hid
    assert [h.handle_id for h in handles.registry.list()] == [hid]
