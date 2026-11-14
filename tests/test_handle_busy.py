"""ハンドルのロックの待ちに上限を置く（レビュー指摘 R13）。

長い仕事（`handle_dispatch.LONG_RUNNING`）は仕事が終わるまでハンドルのロックを持つ。同じハンドルへのほかの
呼び出しが際限なく待つと、デスクトップアプリの中継が打ち切った呼び出しが、仕事の後でロックを取って黙って効く
（`anyio.to_thread.run_sync` はスレッドを見捨てない）。ロックは最長 `HANDLE_LOCK_WAIT` 秒だけ待ち、空かなければ
何もせずに `HANDLE_BUSY`（ロックを持つ仕事の `job_id` と道具）を返す。仕事のスレッド自身は待ち続ける。

包みの振る舞いは、同じ名前の小さなツールを登録した FastMCP に `install` を当てて確かめる（本物の道具の中身に
依らない）。待つところにはすべて期限を付ける。
"""
from __future__ import annotations

import asyncio
import json
import threading
import time

import pytest
from mcp.server.fastmcp import FastMCP

import server
from metabolomix.core import handle_dispatch, handles, mcp_errors, pending_loads, session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import handle_tools

_TIMEOUT = 30.0
_WAIT = 0.3   # テストのロックの待ちの上限（秒）


@pytest.fixture
def fresh(monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    monkeypatch.setattr(handles, "registry", handles.HandleRegistry())
    book = pending_loads.PendingJobs()
    monkeypatch.setattr(pending_loads, "jobs", book)
    monkeypatch.setattr(handle_dispatch, "HANDLE_LOCK_WAIT", _WAIT)
    yield book
    for entry in list(book._entries.values()):
        assert entry.wait(_TIMEOUT) is not None, f"仕事が終わらない: {entry.load_id}"


@pytest.fixture
def app():
    """長い仕事（`dataset_homolog_check` の名前を借りる。`release` を放すまで終わらない）と、状態を書き換える
    小さなツール（`dataset_exclude` の名前を借りる。呼ばれた回数を数える）。"""
    mcp = FastMCP("busy-test")
    state = {"release": threading.Event(), "started": threading.Event(), "mutations": []}

    @mcp.tool(structured_output=False)
    def dataset_homolog_check(spec: str = "x", dataset_handle: str | None = None,
                              wait_seconds: float = 20.0) -> str:
        state["started"].set()
        assert state["release"].wait(_TIMEOUT), "テストが Event を放さなかった"
        return json.dumps({"status": "success", "spec": spec})

    @mcp.tool(structured_output=False)
    def dataset_exclude(samples: list[str], dataset_handle: str | None = None) -> str:
        state["mutations"].append((dataset_handle, list(samples)))
        return json.dumps({"status": "success", "excluded": samples})

    assert handle_dispatch.install(mcp) == ["dataset_exclude", "dataset_homolog_check"]
    try:
        yield mcp, state
    finally:
        state["release"].set()


def _run(mcp, name, args) -> dict:
    content = asyncio.run(mcp.call_tool(name, args))
    blocks = content[0] if isinstance(content, tuple) else content
    return json.loads(blocks[0].text)


def _run_with_deadline(mcp, name, args, seconds: float = 5.0) -> dict:
    """別のスレッドで呼び、`seconds` 秒で返らなければ落とす（際限なく待つ実装でテストが止まらないように）。"""
    box = {}
    thread = threading.Thread(target=lambda: box.setdefault("out", _run(mcp, name, args)), daemon=True)
    thread.start()
    thread.join(seconds)
    assert not thread.is_alive(), f"{name} がハンドルのロックを {seconds} 秒たっても待ち続けている"
    return box["out"]


def _open_handle() -> str:
    session_state.session.dataset = DatasetState()
    return session_state.session.current_handle_id


def test_call_on_a_handle_held_by_a_job_returns_busy_and_does_nothing(fresh, app):
    mcp, state = app
    hid = _open_handle()
    running = _run(mcp, "dataset_homolog_check", {"dataset_handle": hid, "wait_seconds": 0.05})
    assert running["status"] == "running" and running["dataset_handle"] == hid
    job_id = running["job_id"]
    assert state["started"].wait(_TIMEOUT)
    assert handles.registry.peek(hid).running_job["job_id"] == job_id

    for args in ({"samples": ["S3"], "dataset_handle": hid}, {"samples": ["S3"]}):  # 渡す・省く
        started = time.monotonic()
        out = _run_with_deadline(mcp, "dataset_exclude", args)
        error = out["error"]
        assert error["code"] == "HANDLE_BUSY"
        assert error["details"]["dataset_handle"] == hid
        assert error["details"]["job_id"] == job_id
        assert error["details"]["tool"] == "dataset_homolog_check"
        assert error["details"]["next"] == f'job_status(job_id="{job_id}", wait_seconds=30)'
        assert error["required_tools"] == ["job_status"]
        assert "何もしていません" in error["message"]
        assert time.monotonic() - started >= _WAIT * 0.9   # 上限まで待ってから断る
    assert state["mutations"] == []

    state["release"].set()
    done = json.loads(handle_tools.job_status(job_id, wait_seconds=_TIMEOUT))
    assert done["status"] == "success" and done["job_id"] == job_id
    assert handles.registry.peek(hid).running_job is None      # 仕事が終わればロックと印を放す
    time.sleep(_WAIT)
    assert state["mutations"] == []    # 待ちきれなかった呼び出しは、後から黙って効かない
    out = _run(mcp, "dataset_exclude", {"samples": ["S3"], "dataset_handle": hid})
    assert out["status"] == "success" and state["mutations"] == [(hid, ["S3"])]


def test_dataset_close_on_a_busy_handle_closes_nothing(fresh, app):
    mcp, state = app
    hid = _open_handle()
    job_id = _run(mcp, "dataset_homolog_check", {"wait_seconds": 0.05})["job_id"]
    assert state["started"].wait(_TIMEOUT)

    out = _run_with_deadline(server.mcp, "dataset_close", {"dataset_handle": hid})
    assert out["error"]["code"] == "HANDLE_BUSY" and out["error"]["details"]["job_id"] == job_id
    assert handles.registry.peek(hid) is not None and session_state.session.current_handle_id == hid

    state["release"].set()
    assert json.loads(handle_tools.job_status(job_id, wait_seconds=_TIMEOUT))["status"] == "success"
    time.sleep(_WAIT)
    assert handles.registry.peek(hid) is not None       # 断った閉じ方が後から効いていない
    closed = _run(server.mcp, "dataset_close", {"dataset_handle": hid})
    assert closed["status"] == "success" and handles.registry.peek(hid) is None


def test_lock_held_by_another_call_is_busy_without_a_job(fresh, app):
    mcp, state = app
    hid = _open_handle()
    handle = handles.registry.peek(hid)
    held, release = threading.Event(), threading.Event()

    def hold():
        with handle.lock:
            held.set()
            assert release.wait(_TIMEOUT)

    holder = threading.Thread(target=hold, daemon=True)
    holder.start()
    try:
        assert held.wait(_TIMEOUT)
        out = _run_with_deadline(mcp, "dataset_exclude", {"samples": ["S1"], "dataset_handle": hid})
    finally:
        release.set()
        holder.join(_TIMEOUT)
    error = out["error"]
    assert error["code"] == "HANDLE_BUSY" and error["details"]["dataset_handle"] == hid
    assert "job_id" not in error["details"] and "required_tools" not in error
    assert state["mutations"] == []


def test_a_second_job_waits_for_the_lock_instead_of_failing(fresh, app):
    """仕事のスレッド自身は上限を使わずに待つ（受付番号を返してあるので、呼び出し側は job_status で引ける）。"""
    mcp, state = app
    hid = _open_handle()
    first = _run(mcp, "dataset_homolog_check", {"spec": "a", "wait_seconds": 0.05})["job_id"]
    assert state["started"].wait(_TIMEOUT)
    second = _run(mcp, "dataset_homolog_check", {"spec": "b", "wait_seconds": 0.05})
    assert second["status"] == "running"
    assert second["queued_behind"] == {"job_id": first, "tool": "dataset_homolog_check"}
    time.sleep(_WAIT * 2)    # 上限を過ぎても、2 本目の仕事は HANDLE_BUSY で終わらない
    running = json.loads(handle_tools.job_status(second["job_id"]))
    assert running["status"] == "running" and running["queued_behind"]["job_id"] == first

    state["release"].set()
    out_first = json.loads(handle_tools.job_status(first, wait_seconds=_TIMEOUT))
    out_second = json.loads(handle_tools.job_status(second["job_id"], wait_seconds=_TIMEOUT))
    assert out_first["spec"] == "a" and out_second["spec"] == "b"
    assert out_second["status"] == "success" and out_second["dataset_handle"] == hid


def test_handle_busy_envelope_shape():
    with_job = json.loads(mcp_errors.handle_busy(
        "dh_x", 25.04, {"job_id": "job_1", "tool": "pathway_project", "running_s": 40.0}))["error"]
    assert with_job["code"] == "HANDLE_BUSY"
    assert with_job["details"] == {"dataset_handle": "dh_x", "waited_s": 25.0, "job_id": "job_1",
                                   "tool": "pathway_project", "running_s": 40.0,
                                   "next": 'job_status(job_id="job_1", wait_seconds=30)'}
    assert with_job["required_tools"] == ["job_status"]
    plain = json.loads(mcp_errors.handle_busy("dh_x", 25))["error"]
    assert plain["details"] == {"dataset_handle": "dh_x", "waited_s": 25.0}
    assert "required_tools" not in plain
    assert "HANDLE_BUSY" in mcp_errors.MZTAB_ERROR_CODES


def test_lock_wait_stays_under_the_relay_cut():
    # 中継は 36〜60 秒ほどで呼び出しを打ち切る。その前に断る。
    assert 0 < handle_dispatch.HANDLE_LOCK_WAIT < 36
