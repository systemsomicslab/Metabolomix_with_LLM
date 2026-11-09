"""データセットの読み込みをイベントループから外す（B-2a Task 2）。

- `load_dataset` はワーカースレッドで走り、読み込みの間もイベントループが進む。
- 同じ入力を同時に開いても読むのは 1 回（別の入力どうしは並ぶ）。
- セッションのライブラリを新しいハンドルへ移す箇所と `library_load` の書き込みは排他。
"""
from __future__ import annotations

import asyncio
import json
import threading
import time

import pytest

import server
from metabolomix.core import handles, session_state
from metabolomix.dataset import opening
from metabolomix.library import tools as library_tools
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


def test_dataset_load_does_not_block_loop(monkeypatch, pos_mztab):
    original = opening._open

    def slow_open(*args, **kwargs):
        time.sleep(0.5)
        return original(*args, **kwargs)

    monkeypatch.setattr(opening, "_open", slow_open)
    ticks = 0

    async def main():
        nonlocal ticks
        done = asyncio.Event()

        async def ticker():
            nonlocal ticks
            while not done.is_set():
                await asyncio.sleep(0.05)
                ticks += 1

        task = asyncio.create_task(ticker())
        try:
            return await server.mcp.call_tool("load_dataset", {"mztab_path": str(pos_mztab)})
        finally:
            done.set()
            await task

    asyncio.run(main())
    assert ticks >= 5, f"読み込みの間にイベントループが進んでいない（{ticks} 回）"
    assert session_state.session.dataset is not None


def test_load_dataset_runs_in_a_worker_thread():
    """`load_dataset` は常にワーカースレッドで走る（`handle_dispatch.ALWAYS_THREADED`）。"""
    assert server.mcp._tool_manager._tools["load_dataset"].is_async


def test_concurrent_open_same_input_reads_once(monkeypatch, pos_mztab):
    original = opening._open
    calls = []
    first_in = threading.Event()

    def counting_open(*args, **kwargs):
        calls.append(args[1])
        first_in.set()
        time.sleep(0.3)  # 2 つ目のスレッドがキーのロックで待つ時間
        return original(*args, **kwargs)

    monkeypatch.setattr(opening, "_open", counting_open)
    results = {}

    def run(name):
        results[name] = opening.open_dataset(mztab_path=str(pos_mztab))

    t1 = threading.Thread(target=run, args=("a",))
    t1.start()
    assert first_in.wait(10)
    t2 = threading.Thread(target=run, args=("b",))
    t2.start()
    t1.join(30)
    t2.join(30)
    assert len(calls) == 1
    assert results["a"].reused_base is False
    assert results["b"].reused_base is True
    assert results["a"].dataset is not results["b"].dataset  # 複製


def test_concurrent_open_different_inputs_run_in_parallel(monkeypatch, pos_mztab, neg_mztab):
    original = opening._open
    # 両方が入ってからでないとどちらも出られない。直列なら Barrier が時間切れで壊れる。
    barrier = threading.Barrier(2, timeout=10)
    errors = []

    def meeting_open(*args, **kwargs):
        try:
            barrier.wait()
        except threading.BrokenBarrierError as exc:
            errors.append(exc)
            raise
        return original(*args, **kwargs)

    monkeypatch.setattr(opening, "_open", meeting_open)
    threads = [threading.Thread(target=lambda p=p: opening.open_dataset(mztab_path=str(p)))
               for p in (pos_mztab, neg_mztab)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(30)
    assert not errors, "別の入力の読み込みが直列になっている"


def test_inflight_lock_is_per_key():
    assert opening._inflight_lock("a") is opening._inflight_lock("a")
    assert opening._inflight_lock("a") is not opening._inflight_lock("b")


# ---------- セッションのライブラリの移動と library_load ----------

class _FakeStore:
    def __init__(self, name, on_close=None):
        self.name = name
        self.closed = False
        self._on_close = on_close

    def close(self):
        if self._on_close:
            self._on_close()
        self.closed = True

    def summary(self):
        return {"source_sha256": "0" * 64, "record_count": 1, "ion_modes": ["positive"],
                "records_without_ion_mode": 0, "non_utf8_lines": 0}

    def compound_class_counts(self, n):
        return {}


@pytest.fixture
def fake_library_open(monkeypatch):
    new_store = _FakeStore("new")
    monkeypatch.setattr(library_tools, "resolve_library_path", lambda *a, **k: "fake.msp")
    monkeypatch.setattr(library_tools, "open_store", lambda *a, **k: new_store)
    return new_store


def _dataset(mztab):
    return opening.open_dataset(mztab_path=str(mztab)).dataset


def _wait_not_done(thread, seconds=0.4):
    thread.join(seconds)
    return not thread.is_alive()


def test_library_move_and_library_load_do_not_race(monkeypatch, pos_mztab, fake_library_open):
    """setter が移す決定をしてから入れ替えるまでの間に library_load が来ても、
    ライブラリはちょうど 1 つのハンドルに属する（setter が先の順）。"""
    ds = _dataset(pos_mztab)
    session = session_state.session
    in_add, release = threading.Event(), threading.Event()
    real_add = handles.registry.add

    def paused_add(handle):
        in_add.set()
        assert release.wait(10)
        return real_add(handle)

    monkeypatch.setattr(handles.registry, "add", paused_add)
    setter = threading.Thread(target=lambda: setattr(session, "dataset", ds))
    setter.start()
    assert in_add.wait(10)
    loader = threading.Thread(target=lambda: library_tools.library_load(file_path="fake.msp"))
    loader.start()
    finished_early = _wait_not_done(loader)
    release.set()
    setter.join(30)
    loader.join(30)
    assert not finished_early, "library_load が移動の途中で割り込んだ"
    assert session.library.store is fake_library_open     # 見えるライブラリに書かれた
    assert session._library.store is None                 # セッション側に孤児が残らない


def test_library_move_waits_for_library_load_in_progress(pos_mztab, fake_library_open):
    """library_load が書いている最中は、移動の決定が待つ（library_load が先の順）。"""
    ds = _dataset(pos_mztab)
    session = session_state.session
    in_close, release = threading.Event(), threading.Event()

    def slow_close():
        in_close.set()
        assert release.wait(10)

    old = _FakeStore("old", on_close=slow_close)
    session._library.store = old
    loader = threading.Thread(target=lambda: library_tools.library_load(file_path="fake.msp"))
    loader.start()
    assert in_close.wait(10)
    setter = threading.Thread(target=lambda: setattr(session, "dataset", ds))
    setter.start()
    finished_early = _wait_not_done(setter)
    release.set()
    loader.join(30)
    setter.join(30)
    assert not finished_early, "移動が library_load の書き込みの途中で進んだ"
    assert old.closed
    assert session.library.store is fake_library_open
    assert session._library.store is None
