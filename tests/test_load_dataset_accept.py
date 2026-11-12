"""`load_dataset` の受付番号（`ld_…`）と `dataset_load_status`（枠組みの見直し 第 1 部 Task A1）。

Claude のデスクトップアプリの中継は MCP の呼び出しを約 80 秒で打ち切る。`load_dataset` は
`wait_seconds`（0〜60 秒）の中で終われば今までの戻り値をそのまま返し、終わらなければ受付番号を
返す。クライアントは `dataset_load_status(load_id)` で結果を引く。受付番号を返した後に別の
ハンドルが直近になっていたら、終わった読み込みは直近にならない（`current: false`）。上限は
読み込み中の受付も 1 本と数える。

読み込みを止めるテストは `opening.open_dataset` を `threading.Event` で止める（`gate`。`only` に
入れたフォルダだけを止められる）。`gate` は後片付けで必ず Event を放し、`fresh_state` は残った
読み込みが終わるのを待ってから（終わらなければ落とす）ハンドルを閉じる——落ちたテストが全体を
止めない・次のテストへ読み込みを持ち越さない。待つところにはすべて期限を付ける。
"""
from __future__ import annotations

import asyncio
import json
import math
import threading
import time
from pathlib import Path
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
    for entry in list(loads._entries.values()):
        assert entry.wait(_TIMEOUT) is not None, f"読み込みが終わらない: {entry.load_id}"
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def gate(monkeypatch):
    """`opening.open_dataset` を Event で止める。`release` を放すまで読み込みは進まない。

    `state.only` にフォルダ（文字列）の集合を入れると、そのフォルダの読み込みだけを止める。
    """
    release = threading.Event()
    entered = threading.Event()
    real = opening.open_dataset
    state = SimpleNamespace(calls=0, error=None, only=None)

    def blocked(**kwargs):
        if state.only is not None and kwargs.get("directory") not in state.only:
            return real(**kwargs)
        state.calls += 1
        entered.set()
        assert release.wait(_TIMEOUT), "テストが Event を放さなかった"
        if state.error is not None:
            raise state.error
        return real(**kwargs)

    monkeypatch.setattr(opening, "open_dataset", blocked)
    try:
        yield SimpleNamespace(release=release, entered=entered, state=state)
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


def _load(**kwargs) -> dict:
    kwargs.setdefault("overview", False)
    out = json.loads(handle_tools.load_dataset(**kwargs))
    assert out["status"] == "success", out
    return out


def _status_when_done(load_id: str) -> str:
    """読み込みが終わるのを（期限つきで）待ってから `dataset_load_status` を引く。"""
    assert pending_loads.registry.wait(load_id, _TIMEOUT) is not None, "読み込みが終わらない"
    return handle_tools.dataset_load_status(load_id)


class _FakeStore:
    """セッション自身のライブラリに置く偽の store（移されたかを見るだけ）。"""

    def close(self):
        pass

    def summary(self):
        return {"source_sha256": "0" * 64, "record_count": 1}


# ---------- 計画の 6 本 ----------


def test_fast_load_returns_payload_unchanged(tmp_path, fresh_state):
    folder = _folder(tmp_path / "a")
    out = json.loads(handle_tools.load_dataset(directory=str(folder), overview=False,
                                               wait_seconds=20))
    assert out["status"] == "success"
    assert "load_id" not in out and "current" not in out and "note" not in out
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
    # ほかに直近になったハンドルが無いので、受け付けた読み込みが直近になる（`current` の鍵は無い）。
    assert "current" not in done and "note" not in done
    assert session_state.session.current_handle_id == done["dataset_handle"]
    assert [h.handle_id for h in handles.registry.list()] == [done["dataset_handle"]]


def test_status_is_repeatable(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    first = _status_when_done(load_id)
    assert json.loads(first)["status"] == "success"
    assert handle_tools.dataset_load_status(load_id) == first

    # 引いても直近は動かない: 別のハンドルを開いた後に引き直しても、直近はその別のハンドルのまま。
    other = _load(directory=str(_folder(tmp_path / "b")))
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
    # 待つ秒を渡しても同じ封筒（待たずに返す）。
    waited = json.loads(handle_tools.dataset_load_status("ld_nope", wait_seconds=5))
    assert waited["error"]["code"] == "DATASET_LOAD_NOT_FOUND"


# ---------- dataset_load_status の wait_seconds（Ruling 29） ----------


def test_status_wait_returns_payload_when_load_finishes_within_wait(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    # 待っている間に読み込みを終わらせる: 1 回の呼び出しで中身が返る。
    timer = threading.Timer(0.2, gate.release.set)
    timer.start()
    try:
        raw = handle_tools.dataset_load_status(load_id, wait_seconds=_TIMEOUT)
    finally:
        timer.cancel()
    out = json.loads(raw)
    assert out["status"] == "success", out
    assert session_state.session.current_handle_id == out["dataset_handle"]
    # 待たずに引き直しても同じ文字列（受付の結果は 1 つ）。
    assert handle_tools.dataset_load_status(load_id) == raw


def test_status_wait_returns_loading_after_about_the_wait(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    start = time.monotonic()
    out = json.loads(handle_tools.dataset_load_status(load_id, wait_seconds=0.3))
    took = time.monotonic() - start
    assert out["status"] == "loading" and out["load_id"] == load_id
    assert set(out) == {"status", "load_id", "elapsed_s"}
    assert took >= 0.25, took  # 待った
    assert took < _TIMEOUT
    assert handles.registry.list() == []
    gate.release.set()
    assert json.loads(_status_when_done(load_id))["status"] == "success"


def test_status_wait_zero_does_not_wait(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    start = time.monotonic()
    out = json.loads(handle_tools.dataset_load_status(load_id, wait_seconds=0))
    assert out["status"] == "loading"
    assert time.monotonic() - start < 0.25
    gate.release.set()
    assert json.loads(_status_when_done(load_id))["status"] == "success"


@pytest.mark.parametrize("bad", [-1, -0.001, 60.0001, 61, math.nan, math.inf, -math.inf,
                                 True, False, "20", None])
def test_status_wait_rejects_bad_values(tmp_path, fresh_state, gate, bad):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    out = json.loads(handle_tools.dataset_load_status(load_id, wait_seconds=bad))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST", out
    assert out["error"]["details"] == {"field": "wait_seconds"}
    gate.release.set()
    assert json.loads(_status_when_done(load_id))["status"] == "success"


def test_status_wait_accepts_the_limit(tmp_path, fresh_state, gate):
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    assert pending_loads.registry.wait(load_id, _TIMEOUT) is not None
    out = json.loads(handle_tools.dataset_load_status(load_id, wait_seconds=60))
    assert out["status"] == "success"


def test_capacity_checked_before_accept(tmp_path, fresh_state, monkeypatch):
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    first = _load(directory=str(_folder(tmp_path / "a")))

    def must_not_accept(fn):
        raise AssertionError("上限なのに受け付けた")

    monkeypatch.setattr(fresh_state, "submit", must_not_accept)
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "b")),
                                               overview=False, wait_seconds=0.05))
    assert "load_id" not in out
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert "loading" not in out["error"]["details"]  # 読み込み中の受付は無い
    assert session_state.session.current_handle_id == first["dataset_handle"]


# ---------- 直近の切り替え（受け付けた後に別のハンドルが直近になったら） ----------


def test_slow_load_does_not_steal_current_from_a_newer_load(tmp_path, fresh_state, gate):
    folder_a, folder_b = str(_folder(tmp_path / "a")), str(_folder(tmp_path / "b"))
    gate.state.only = {folder_a}
    load_id = _accept(directory=folder_a)["load_id"]
    b = _load(directory=folder_b)  # 止めていないので期限の中で終わり、直近になる
    assert session_state.session.current_handle_id == b["dataset_handle"]

    gate.release.set()
    raw = _status_when_done(load_id)
    a = json.loads(raw)
    assert a["status"] == "success"
    assert a["current"] is False and a["note"]
    assert session_state.session.current_handle_id == b["dataset_handle"]  # B のまま
    assert {h.handle_id for h in handles.registry.list()} == {a["dataset_handle"],
                                                                b["dataset_handle"]}
    # 足したのは `current` と `note` の 2 つだけ（ほかの欄は期限の中で終わったときと同じ）。
    assert set(a) - {"caveat", "current", "note"} == set(b) - {"caveat"}
    assert handle_tools.dataset_load_status(load_id) == raw


def test_not_current_load_does_not_move_the_session_library(tmp_path, fresh_state, gate):
    """直近にしない読み込みは、セッション自身のライブラリを移さない（直近が空いていても）。"""
    folder_a, folder_b = str(_folder(tmp_path / "a")), str(_folder(tmp_path / "b"))
    gate.state.only = {folder_a}
    load_id = _accept(directory=folder_a)["load_id"]
    b = _load(directory=folder_b)
    handles.registry.close(b["dataset_handle"])  # 直近が空く（世代は戻らない）
    store = _FakeStore()
    session_state.session._library.store = store  # `library_load`（ハンドルなし）が読んだ形

    gate.release.set()
    a = json.loads(_status_when_done(load_id))
    assert a["status"] == "success" and a["current"] is False
    assert session_state.session._library.store is store
    assert handles.registry.peek(a["dataset_handle"]).library.store is None
    assert session_state.session.current_handle_id != a["dataset_handle"]
    assert session_state.session.dataset is None  # 閉じた B の後、直近は空のまま


def test_fast_load_still_becomes_current_even_if_another_did_meanwhile(tmp_path, fresh_state,
                                                                      gate):
    """期限の中で終わった読み込みは今までどおり（終わった方が直近になる）。"""
    folder_a, folder_b = str(_folder(tmp_path / "a")), str(_folder(tmp_path / "b"))
    gate.state.only = {folder_a}
    results = []
    thread = threading.Thread(target=lambda: results.append(json.loads(
        handle_tools.load_dataset(directory=folder_a, overview=False, wait_seconds=30))))
    thread.start()
    try:
        assert gate.entered.wait(_TIMEOUT)
        b = _load(directory=folder_b)
        assert session_state.session.current_handle_id == b["dataset_handle"]
    finally:
        gate.release.set()
        thread.join(_TIMEOUT)
    assert not thread.is_alive()
    a = results[0]
    assert a["status"] == "success" and "current" not in a and "load_id" not in a
    assert session_state.session.current_handle_id == a["dataset_handle"]


# ---------- 上限（読み込み中の受付も数える） ----------


def test_in_flight_load_fills_the_last_slot(tmp_path, fresh_state, gate, monkeypatch):
    """上限 1 で 1 本を受け付けると、読み込み中のそれが最後の枠を取り、次は受け付けずに断る。"""
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    a = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "b")),
                                               overview=False, wait_seconds=0.05))
    assert "load_id" not in out
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert out["error"]["details"]["loading"] == 1
    assert out["error"]["details"]["open_handles"] == []
    assert gate.state.calls == 1  # 2 本目は読み込みを始めていない

    gate.release.set()
    done = json.loads(_status_when_done(a))
    assert done["status"] == "success"
    assert session_state.session.current_handle_id == done["dataset_handle"]
    # 終わった後は開いたハンドルが枠を使う（予約は登録で消えている。数え直さない）。
    assert handles.registry._reserved == 0
    again = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "c")),
                                                 overview=False))
    assert again["error"]["code"] == "DATASET_HANDLE_LIMIT" and "loading" not in again["error"]["details"]


def test_failed_load_returns_its_slot(tmp_path, fresh_state, gate, monkeypatch):
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    gate.state.error = DomainError("ALIGNMENT_NOT_FOUND", "合う .arf がありません")
    a = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    assert json.loads(_status_when_done(a))["error"]["code"] == "ALIGNMENT_NOT_FOUND"
    assert handles.registry._reserved == 0
    gate.state.error = None
    _load(directory=str(_folder(tmp_path / "b")))  # 枠が返っているので開ける


def test_adopt_recheck_is_a_backstop(tmp_path, fresh_state, gate, monkeypatch):
    """受付の後に（予約を通らない道で）枠が埋まったら、登録の検査が `DATASET_HANDLE_LIMIT` の封筒で
    止め、ハンドルを登録せず直近も動かさない。"""
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 1})
    a = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    assert gate.entered.wait(_TIMEOUT)  # 本体の前置きの検査は通り過ぎている
    other = handles.DatasetHandle(object(), polarity=None, polarity_source=None, source={})
    handles.registry.add(other)
    gate.release.set()
    out = json.loads(_status_when_done(a))
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT"
    assert [h.handle_id for h in handles.registry.list()] == [other.handle_id]
    assert session_state.session.current_handle_id is None
    assert handles.registry._reserved == 0


# ---------- 引数・失敗・受付簿 ----------


@pytest.mark.parametrize("value", [-1, -0.001, 60.0001, 61, math.nan, math.inf, -math.inf,
                                   True, "20", None])
def test_invalid_wait_seconds_is_bad_request(tmp_path, fresh_state, value):
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "a")),
                                               overview=False, wait_seconds=value))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert out["error"]["details"] == {"field": "wait_seconds"}
    assert handles.registry.list() == []
    assert handles.registry._reserved == 0
    assert fresh_state._entries == {}


def test_wait_seconds_60_is_accepted(tmp_path, fresh_state):
    assert _load(directory=str(_folder(tmp_path / "a")), wait_seconds=60)["status"] == "success"


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


def test_domain_error_with_unserializable_details_after_accept(tmp_path, fresh_state, gate):
    """`details` に JSON にならない値（`Path`）を持つ DomainError でも、受付は誤りの封筒で終わる。"""
    gate.state.error = DomainError("ALIGNMENT_NOT_FOUND", "合う .arf がありません",
                                   {"path": Path(tmp_path)})
    load_id = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    gate.release.set()
    out = json.loads(_status_when_done(load_id))
    assert "error" in out and out.get("status") != "loading"
    assert handles.registry.list() == []
    assert handles.registry._reserved == 0


def test_result_survives_eviction_between_submit_and_wait(tmp_path, fresh_state, monkeypatch):
    """終わった受付が、`load_dataset` が待つ前に受付簿から追い出されても、結果は失われない。"""
    loads = pending_loads.PendingLoads(max_entries=1)
    monkeypatch.setattr(pending_loads, "registry", loads)
    real_submit = loads.submit

    def submit_then_evict(fn):
        pending = real_submit(fn)
        assert pending.wait(_TIMEOUT) is not None
        real_submit(lambda: "other")  # 終わった最初の受付を追い出す
        with pytest.raises(DomainError) as info:
            loads.status(pending.load_id)
        assert info.value.code == "DATASET_LOAD_NOT_FOUND"
        return pending

    monkeypatch.setattr(loads, "submit", submit_then_evict)
    out = _load(directory=str(_folder(tmp_path / "a")))
    assert session_state.session.current_handle_id == out["dataset_handle"]


def test_capacity_with_concurrent_accepts(tmp_path, fresh_state, gate, monkeypatch):
    """上限 2 で 2 本を受け付けると 3 本目は断り、2 本とも終われば両方がハンドルになる。"""
    use_config(monkeypatch, tmp_path, {"server.max_dataset_handles": 2})
    a = _accept(directory=str(_folder(tmp_path / "a")))["load_id"]
    b = _accept(directory=str(_folder(tmp_path / "b")))["load_id"]
    out = json.loads(handle_tools.load_dataset(directory=str(_folder(tmp_path / "c")),
                                               overview=False, wait_seconds=0.05))
    assert out["error"]["code"] == "DATASET_HANDLE_LIMIT" and out["error"]["details"]["loading"] == 2
    gate.release.set()
    outs = [json.loads(_status_when_done(x)) for x in (a, b)]
    assert all(o["status"] == "success" for o in outs)
    assert len(handles.registry.list()) == 2


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


def _not_found(loads, load_id):
    with pytest.raises(DomainError) as info:
        loads.status(load_id)
    assert info.value.code == "DATASET_LOAD_NOT_FOUND"


def test_submit_returns_the_entry_and_wait_returns_the_payload():
    loads = pending_loads.PendingLoads()
    pending = loads.submit(lambda: "x")
    assert pending.load_id.startswith("ld_")
    assert pending.wait(_TIMEOUT) == "x"
    assert loads.wait(pending.load_id, _TIMEOUT) == "x"
    # 終わった後は 0 秒でも返す（「既に終わっていれば返す」）。
    assert loads.wait(pending.load_id, 0) == "x"
    assert loads.status(pending.load_id) == ("done", "x")


def test_wait_returns_none_while_loading():
    loads = pending_loads.PendingLoads()
    release = threading.Event()
    try:
        load_id = loads.submit(_blocked(release, "x")).load_id
        assert loads.wait(load_id, 0) is None
        assert loads.wait(load_id, 0.05) is None
        assert loads.status(load_id) == ("loading", None)
    finally:
        release.set()
    assert loads.wait(load_id, _TIMEOUT) == "x"


def test_finished_entries_expire_after_an_hour():
    clock = _Clock()
    loads = pending_loads.PendingLoads(clock=clock)
    load_id = loads.submit(lambda: "x").load_id
    assert loads.wait(load_id, _TIMEOUT) == "x"
    clock.now = pending_loads.TTL_SECONDS - 0.5
    assert loads.status(load_id) == ("done", "x")
    clock.now = pending_loads.TTL_SECONDS
    _not_found(loads, load_id)


def test_loading_entries_never_expire():
    clock = _Clock()
    loads = pending_loads.PendingLoads(clock=clock)
    release = threading.Event()
    try:
        load_id = loads.submit(_blocked(release, "x")).load_id
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
        a = loads.submit(_blocked(release_a, "a")).load_id
        clock.now = 1.0
        b = loads.submit(lambda: "b").load_id
        assert loads.wait(b, _TIMEOUT) == "b"
        clock.now = 2.0
        release_a.set()
        assert loads.wait(a, _TIMEOUT) == "a"
        c = loads.submit(_blocked(release_c, "c")).load_id  # 3 件目（まだ空きがある）
        # 4 件目: 終わったものから、終わったのが古い b を捨てる（受け付けた順の a ではない）。
        d = loads.submit(_blocked(release_d, "d")).load_id
        _not_found(loads, b)
        assert loads.status(a) == ("done", "a")
        # 5 件目: 終わったものは a だけ。a を捨てる。
        e = loads.submit(_blocked(release_e, "e")).load_id
        _not_found(loads, a)
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
    assert loads.submit(lambda: "g").wait(_TIMEOUT) == "g"


def test_worker_exception_finishes_the_entry_with_an_envelope():
    loads = pending_loads.PendingLoads()

    def boom():
        raise ValueError("C:/lab_secret_dir/x")

    raw = loads.submit(boom).wait(_TIMEOUT)
    assert json.loads(raw)["error"]["code"] == "DATASET_LOAD_FAILED"
    assert "lab_secret_dir" not in raw

    def domain():
        raise DomainError("CONFIG_INVALID", "設定ファイルが読めません", {"key": "x"})

    assert loads.submit(domain).wait(_TIMEOUT) == mztab_error(
        "CONFIG_INVALID", "設定ファイルが読めません", {"key": "x"})

    # 文字列でない戻り値も封筒にする（`wait` の None＝終わっていない、と取り違えない）。
    raw = loads.submit(lambda: None).wait(_TIMEOUT)
    assert json.loads(raw)["error"]["details"] == {"exception": "TypeError"}


@pytest.mark.parametrize("error", [
    DomainError("ALIGNMENT_NOT_FOUND", "合う .arf がありません", {"path": Path("x")}),  # JSON にならない
    DomainError("", ""),  # mztab_error が ValueError を出す
])
def test_unconvertible_domain_error_still_finishes_and_can_be_evicted(error):
    """封筒づくり（`mztab_error`）が落ちても受付は終わり、終わったものとして追い出せる。"""
    loads = pending_loads.PendingLoads(max_entries=1)

    def raise_it():
        raise error

    pending = loads.submit(raise_it)
    raw = pending.wait(_TIMEOUT)
    assert raw is not None and "error" in json.loads(raw)
    assert loads.status(pending.load_id) == ("done", raw)
    # 終わっているので、次の受付が追い出せる（読み込み中のまま残っていれば DATASET_LOAD_LIMIT）。
    assert loads.submit(lambda: "next").wait(_TIMEOUT) == "next"
    _not_found(loads, pending.load_id)


def test_details_are_dropped_but_code_kept_when_not_serializable():
    loads = pending_loads.PendingLoads()

    def raise_it():
        raise DomainError("ALIGNMENT_NOT_FOUND", "合う .arf がありません", {"path": Path("x")})

    out = json.loads(loads.submit(raise_it).wait(_TIMEOUT))
    assert out["error"]["code"] == "ALIGNMENT_NOT_FOUND"
    assert out["error"]["details"] == {"details_dropped": True}


def test_discard_forgets_the_entry():
    loads = pending_loads.PendingLoads()
    pending = loads.submit(lambda: "x")
    assert pending.wait(_TIMEOUT) == "x"
    loads.discard(pending.load_id)
    loads.discard(pending.load_id)  # 2 回目も誤りにしない
    _not_found(loads, pending.load_id)
    assert pending.wait(0) == "x"  # 受付そのものは結果を持ち続ける


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


def test_load_status_keeps_library_names_choice(tmp_path, fresh_state, gate, monkeypatch):
    """遅い読み込み（受付番号）でも、`library_names=False` の選択が最終の戻り値に引き継がれる。"""
    from tests.dbs_fixture import record, write_dbs

    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    name = "Zebra_SecretLab_Loaded.msp2.dbs"
    folder = _folder(tmp_path / "a")
    write_dbs(folder / name, [("lib", [record(mz=500.0, ion_mode=0)], ["k"])])
    load_id = _accept(directory=str(folder), library_names=False)["load_id"]
    gate.release.set()
    raw = _status_when_done(load_id)
    done = json.loads(raw)
    assert done["status"] == "success" and done["library"]["sha256"]
    assert name not in raw and "Zebra" not in raw
