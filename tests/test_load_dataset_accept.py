"""`load_dataset` の受付番号（`ld_…`）と `dataset_load_status`（枠組みの見直し 第 1 部 Task A1）。

Claude のデスクトップアプリの中継は MCP の呼び出しを約 80 秒で打ち切る。`load_dataset` は
`wait_seconds` の中で終われば今までの戻り値をそのまま返し、終わらなければ受付番号を返す。
クライアントは `dataset_load_status(load_id)` で結果を引く。

読み込みを止めるテストは `opening.open_dataset` を `threading.Event` で止める（`gate`）。
`gate` は後片付けで必ず Event を放し、`fresh_state` は残った読み込みが終わるのを待ってから
ハンドルを閉じる（落ちたテストが全体を止めない・次のテストへ読み込みを持ち越さない）。
待つところにはすべて期限を付ける。
"""
from __future__ import annotations

import asyncio
import json
import math
import threading
from types import SimpleNamespace

import pytest

import server
from metabolomix.core import handle_dispatch, handles, pending_loads, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.dataset import opening
from metabolomix.tools import handle_tools
from tests.config_fixtures import use_config
from tests.msdial_folder_fixtures import write_msdial_folder

#: テストの中で待つ上限（秒）。これを超えたら壊れているとみなす。
_TIMEOUT = 30.0

_SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
            for i in range(4)]
_SPOTS = [{"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0]},
          {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]}]


@pytest.fixture
def fresh_state(monkeypatch):
    session_state.session = session_state.AnalysisSession()
    loads = pending_loads.PendingLoads()
    monkeypatch.setattr(pending_loads, "registry", loads)
    yield loads
    # 残った読み込みが終わるまで待つ（終わる前に片付けると、次のテストの登録簿へ入りうる）。
    for load_id in list(loads._entries):
        loads.wait(load_id, _TIMEOUT)
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def gate(monkeypatch):
    """`opening.open_dataset` を Event で止める。`release` を放すまで読み込みは進まない。"""
    release = threading.Event()
    real = opening.open_dataset
    state = SimpleNamespace(calls=0, error=None)

    def blocked(**kwargs):
        state.calls += 1
        assert release.wait(_TIMEOUT), "テストが Event を放さなかった"
        if state.error is not None:
            raise state.error
        return real(**kwargs)

    monkeypatch.setattr(opening, "open_dataset", blocked)
    try:
        yield SimpleNamespace(release=release, state=state)
    finally:
        release.set()


def _folder(path):
    write_msdial_folder(path, spots=_SPOTS, samples=_SAMPLES)
    return path


def _accept(**kwargs) -> dict:
    """期限の外へ出る読み込み（`gate` で止めた）を受け付けさせ、受付番号の戻り値を返す。"""
    kwargs.setdefault("overview", False)
    kwargs.setdefault("wait_seconds", 0.05)
    out = json.loads(handle_tools.load_dataset(**kwargs))
    assert out["status"] == "loading", out
    return out


def _status_when_done(load_id: str) -> str:
    """読み込みが終わるのを（期限つきで）待ってから `dataset_load_status` を引く。"""
    assert pending_loads.registry.wait(load_id, _TIMEOUT) is not None, "読み込みが終わらない"
    return handle_tools.dataset_load_status(load_id)


# ---------- 計画の 6 本 ----------


def test_fast_load_returns_payload_unchanged(tmp_path, fresh_state):
    folder = _folder(tmp_path / "a")
    out = json.loads(handle_tools.load_dataset(directory=str(folder), overview=False,
                                               wait_seconds=20))
    assert out["status"] == "success"
    assert "load_id" not in out
    assert handles.registry.peek(out["dataset_handle"]) is not None
    assert session_state.session.current_handle_id == out["dataset_handle"]
    # 戻り値の形は本体（`_load_now`）をそのまま呼んだときと同じ。
    direct = json.loads(handle_tools._load_now(directory=str(folder), overview=False))
    assert set(out) - {"caveat"} == set(direct) - {"caveat"}
    # 期限の中で終わった受付は覚えておかない（受付番号を誰にも渡していない）。
    assert fresh_state._entries == {}


def test_slow_load_returns_load_id_then_payload(tmp_path, fresh_state, gate):
    out = _accept(directory=str(_folder(tmp_path / "a")))
    load_id = out["load_id"]
    assert load_id.startswith("ld_")
    assert out["next"] == f'dataset_load_status(load_id="{load_id}")'
    assert set(out) == {"status", "load_id", "next"}
    assert handles.registry.list() == []
    assert session_state.session.current_handle_id is None

    loading = json.loads(handle_tools.dataset_load_status(load_id))
    assert loading["status"] == "loading" and loading["load_id"] == load_id
    assert set(loading) == {"status", "load_id", "elapsed_s"}
    assert loading["elapsed_s"] >= 0
    assert handles.registry.list() == []  # 引いても何も登録しない

    gate.release.set()
    done = json.loads(_status_when_done(load_id))
    assert done["status"] == "success"
    assert session_state.session.current_handle_id == done["dataset_handle"]
    assert [h.handle_id for h in handles.registry.list()] == [done["dataset_handle"]]


def test_status_is_repeatable(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    first = _status_when_done(load_id)
    assert json.loads(first)["status"] == "success"
    assert handle_tools.dataset_load_status(load_id) == first

    # 引いても直近は動かない: 別のハンドルを開いた後に引き直しても、直近はその別のハンドルのまま。
    other = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "b")),
                                                 overview=False))
    assert other["status"] == "success"
    assert handle_tools.dataset_load_status(load_id) == first
    assert session_state.session.current_handle_id == other["dataset_handle"]
    assert len(handles.registry.list()) == 2


def test_failure_after_accept_returns_envelope(tmp_path, fresh_state, gate):
    folder = str(_folder(tmp_path / "a"))
    gate.state.error = DomainError("ALIGNMENT_NOT_FOUND", "合う .arf がありません",
                                   {"candidates": ["x.arf"]})
    load_id = _accept(directory=folder)["load_id"]
    gate.release.set()
    raw = _status_when_done(load_id)
    out = json.loads(raw)
    assert "status" not in out
    assert out["error"]["code"] == "ALIGNMENT_NOT_FOUND"
    assert out["error"]["details"] == {"candidates": ["x.arf"]}
    # 同期の道（本体をそのまま呼ぶ）と同じ封筒の文字列。
    assert raw == mztab_error("ALIGNMENT_NOT_FOUND", "合う .arf がありません",
                              {"candidates": ["x.arf"]})
    assert raw == handle_tools._load_now(directory=folder, overview=False)
    assert handles.registry.list() == []
    assert session_state.session.current_handle_id is None


def test_unknown_load_id(fresh_state):
    out = json.loads(handle_tools.dataset_load_status("ld_nope"))
    assert out["error"]["code"] == "DATASET_LOAD_NOT_FOUND"
    assert out["error"]["details"] == {"load_id": "ld_nope"}


def test_capacity_checked_before_accept(tmp_path, fresh_state, monkeypatch):
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    first = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "a")),
                                                 overview=False))
    assert first["status"] == "success"

    def must_not_accept(fn):
        raise AssertionError("上限なのに受け付けた")

    monkeypatch.setattr(fresh_state, "submit", must_not_accept)
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "b")),
                                               overview=False, wait_seconds=0.05))
    assert "load_id" not in out
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert session_state.session.current_handle_id == first["dataset_handle"]


# ---------- 引数・失敗・同時の受付 ----------


@pytest.mark.parametrize("value", [-1, -0.001, math.nan, math.inf, -math.inf, True, "20", None])
def test_invalid_wait_seconds_is_bad_request(tmp_path, fresh_state, value):
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "a")),
                                               overview=False, wait_seconds=value))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert out["error"]["details"] == {"field": "wait_seconds"}
    assert handles.registry.list() == []
    assert fresh_state._entries == {}


def test_wait_zero_accepts_at_once(tmp_path, fresh_state, gate):
    out = _accept(directory=str(_folder(tmp_path / "a")), wait_seconds=0)
    gate.release.set()
    assert json.loads(_status_when_done(out["load_id"]))["status"] == "success"


def test_argument_errors_are_returned_before_accept(tmp_path, fresh_state, gate):
    """引数の誤りは受け付けずにすぐ返す（`wait_seconds=0` でも受付番号にしない）。"""
    folder = str(_folder(tmp_path / "a"))
    out = json.loads(handle_tools.load_dataset(directory=folder, overview="no", wait_seconds=0))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert out["error"]["details"] == {"field": "overview"}
    out = json.loads(handle_tools.load_dataset(directory=folder, library="lab", wait_seconds=0))
    assert out["error"]["details"] == {"field": "library"}
    assert gate.state.calls == 0
    assert fresh_state._entries == {}


def test_unexpected_exception_after_accept_returns_envelope(tmp_path, fresh_state, gate):
    """読み込みが DomainError 以外で落ちても、受付は誤りの封筒で終わる（読み込み中のまま残らない）。
    例外の文言（置き場所を含みうる）は戻り値に出さず、種類だけを書く。"""
    gate.state.error = RuntimeError(r"boom C:\lab_secret_dir\pos.msp")
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    raw = _status_when_done(load_id)
    out = json.loads(raw)
    assert out["error"]["code"] == "DATASET_LOAD_FAILED"
    assert out["error"]["details"] == {"exception": "RuntimeError"}
    assert "lab_secret_dir" not in raw
    assert handles.registry.list() == []
    assert handle_tools.dataset_load_status(load_id) == raw


def test_capacity_race_after_accept_returns_envelope(tmp_path, fresh_state, gate, monkeypatch):
    """上限 1 で 2 本を同時に受け付けると、先に終わった方だけがハンドルになり、もう片方は
    `DATASET_HANDLE_LIMIT` の封筒で終わる（ハンドルは登録しない・直近も動かさない）。"""
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    a = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    b = _accept(directory=str(_folder(tmp_path / "b")))["load_id"]
    gate.release.set()
    outs = [json.loads(_status_when_done(load_id)) for load_id in (a, b)]
    successes = [o for o in outs if o.get("status") == "success"]
    errors = [o for o in outs if "error" in o]
    assert len(successes) == 1 and len(errors) == 1
    assert errors[0]["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert [h.handle_id for h in handles.registry.list()] == [successes[0]["dataset_handle"]]
    assert session_state.session.current_handle_id == successes[0]["dataset_handle"]


# ---------- 受付の控え（PendingLoads） ----------


class _Clock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


def _blocked(release: threading.Event, value: str):
    def fn():
        assert release.wait(_TIMEOUT), "テストが Event を放さなかった"
        return value
    return fn


def test_submit_returns_ld_ids_and_wait_returns_the_payload():
    loads = pending_loads.PendingLoads()
    load_id = loads.submit(lambda: "x")
    assert load_id.startswith("ld_")
    assert loads.wait(load_id, _TIMEOUT) == "x"
    # 終わった後は 0 秒でも返す（「既に終わっていれば返す」）。
    assert loads.wait(load_id, 0) == "x"
    assert loads.status(load_id) == ("done", "x")


def test_wait_returns_none_while_loading():
    loads = pending_loads.PendingLoads()
    release = threading.Event()
    try:
        load_id = loads.submit(_blocked(release, "x"))
        assert loads.wait(load_id, 0) is None
        assert loads.wait(load_id, 0.05) is None
        assert loads.status(load_id) == ("loading", None)
    finally:
        release.set()
    assert loads.wait(load_id, _TIMEOUT) == "x"


def test_finished_entries_expire_after_an_hour():
    clock = _Clock()
    loads = pending_loads.PendingLoads(clock=clock)
    load_id = loads.submit(lambda: "x")
    assert loads.wait(load_id, _TIMEOUT) == "x"
    clock.now = pending_loads.TTL_SECONDS - 0.5
    assert loads.status(load_id) == ("done", "x")
    clock.now = pending_loads.TTL_SECONDS
    with pytest.raises(DomainError) as info:
        loads.status(load_id)
    assert info.value.code == "DATASET_LOAD_NOT_FOUND"


def test_loading_entries_never_expire():
    clock = _Clock()
    loads = pending_loads.PendingLoads(clock=clock)
    release = threading.Event()
    try:
        load_id = loads.submit(_blocked(release, "x"))
        clock.now = 10 * pending_loads.TTL_SECONDS
        assert loads.status(load_id) == ("loading", None)
    finally:
        release.set()
    assert loads.wait(load_id, _TIMEOUT) == "x"
    # 1 時間は終わったときから数える。
    assert loads.status(load_id) == ("done", "x")


def test_cap_evicts_oldest_finished_first_and_never_a_loading_one():
    clock = _Clock()
    loads = pending_loads.PendingLoads(clock=clock, max_entries=3)
    releases = [threading.Event() for _ in range(4)]
    release_a, release_c, release_d, release_e = releases
    try:
        # a は先に受け付けたが後で終わる。b は後に受け付けて先に終わる。
        a = loads.submit(_blocked(release_a, "a"))
        clock.now = 1.0
        b = loads.submit(lambda: "b")
        assert loads.wait(b, _TIMEOUT) == "b"
        clock.now = 2.0
        release_a.set()
        assert loads.wait(a, _TIMEOUT) == "a"
        c = loads.submit(_blocked(release_c, "c"))  # 3 件目（まだ空きがある）
        # 4 件目: 終わったものから、終わったのが古い b を捨てる（受け付けた順の a ではない）。
        d = loads.submit(_blocked(release_d, "d"))
        with pytest.raises(DomainError):
            loads.status(b)
        assert loads.status(a) == ("done", "a")
        # 5 件目: 終わったものは a だけ。a を捨てる。
        e = loads.submit(_blocked(release_e, "e"))
        with pytest.raises(DomainError):
            loads.status(a)
        # c・d・e がすべて読み込み中なら、受け付けずに断る（読み込み中のものは捨てない）。
        with pytest.raises(DomainError) as info:
            loads.submit(lambda: "f")
        assert info.value.code == "DATASET_LOAD_LIMIT"
        assert {loads.status(x)[0] for x in (c, d, e)} == {"loading"}
    finally:
        for event in releases:
            event.set()
    for x, value in ((c, "c"), (d, "d"), (e, "e")):
        assert loads.wait(x, _TIMEOUT) == value
    # 1 つ終われば受け付ける。
    assert loads.wait(loads.submit(lambda: "g"), _TIMEOUT) == "g"


def test_worker_exception_finishes_the_entry_with_an_envelope():
    loads = pending_loads.PendingLoads()

    def boom():
        raise ValueError("C:/lab_secret_dir/x")

    load_id = loads.submit(boom)
    raw = loads.wait(load_id, _TIMEOUT)
    assert json.loads(raw)["error"]["code"] == "DATASET_LOAD_FAILED"
    assert "lab_secret_dir" not in raw

    def domain():
        raise DomainError("CONFIG_INVALID", "設定ファイルが読めません", {"key": "x"})

    load_id = loads.submit(domain)
    assert loads.wait(load_id, _TIMEOUT) == mztab_error("CONFIG_INVALID", "設定ファイルが読めません",
                                                        {"key": "x"})

    # 文字列でない戻り値も封筒にする（`wait` の None＝終わっていない、と取り違えない）。
    load_id = loads.submit(lambda: None)
    raw = loads.wait(load_id, _TIMEOUT)
    assert json.loads(raw)["error"]["details"] == {"exception": "TypeError"}


def test_discard_forgets_the_entry():
    loads = pending_loads.PendingLoads()
    load_id = loads.submit(lambda: "x")
    assert loads.wait(load_id, _TIMEOUT) == "x"
    loads.discard(load_id)
    loads.discard(load_id)  # 2 回目も誤りにしない
    with pytest.raises(DomainError):
        loads.status(load_id)


# ---------- 公開面 ----------


def test_dataset_load_status_is_registered_and_threaded():
    assert "dataset_load_status" in handle_dispatch.ALWAYS_THREADED
    tool = server.mcp._tool_manager._tools["dataset_load_status"]
    assert tool.is_async
    assert tool.annotations.readOnlyHint is True


def test_accept_through_the_mcp_layer(tmp_path, fresh_state, gate):
    """MCP の呼び出し（引数の検証とスレッドへの退避）を通しても、受付番号から結果を引ける。"""
    folder = str(_folder(tmp_path / "a"))

    def call(name, args):
        content = asyncio.run(server.mcp.call_tool(name, args))
        blocks = content[0] if isinstance(content, tuple) else content
        return json.loads(blocks[0].text)

    out = call("load_dataset", {"directory": folder, "overview": False, "wait_seconds": 0.05})
    assert out["status"] == "loading"
    gate.release.set()
    assert pending_loads.registry.wait(out["load_id"], _TIMEOUT) is not None
    done = call("dataset_load_status", {"load_id": out["load_id"]})
    assert done["status"] == "success"
    assert session_state.session.current_handle_id == done["dataset_handle"]
