"""長い道具の受付番号（`job_…`）と `job_status`、`wait_seconds` の丸め（棚卸し B10-T27-T43・C6-M09-M10）。

`handle_dispatch.LONG_RUNNING` のツール（`pathway_project` / `dataset_homolog_check`）は、`wait_seconds` の中で
終わらなければ仕事を続けたまま `{"status": "running", "job_id"}` を返す。包みの振る舞いは、同じ名前の小さな
ツールを登録した FastMCP に `install` を当てて確かめる（本物の道具の中身に依らない）。待つところにはすべて
期限を付ける。
"""
from __future__ import annotations

import asyncio
import json
import math
import threading

import pytest
from mcp.server.fastmcp import FastMCP

import server
from metabolomix.core import handle_dispatch, handles, pending_loads, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import handle_tools

_TIMEOUT = 30.0


@pytest.fixture
def fresh(monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    book = pending_loads.PendingJobs()
    monkeypatch.setattr(pending_loads, "jobs", book)
    yield book
    for entry in list(book._entries.values()):
        assert entry.wait(_TIMEOUT) is not None, f"仕事が終わらない: {entry.load_id}"


@pytest.fixture
def mini():
    """`dataset_homolog_check` の名前を借りた小さなツール（`release` を放すまで終わらない）。"""
    app = FastMCP("job-test")
    state = {"release": threading.Event(), "started": threading.Event(), "error": None,
             "seen_handle": None, "lock_held": None}

    @app.tool(structured_output=False)
    def dataset_homolog_check(spec: str = "x", dataset_handle: str | None = None,
                              wait_seconds: float = 20.0) -> str:
        state["started"].set()
        handle = handles.registry.peek(session_state.session.current_handle_id)
        if handle is not None:
            # 受付のスレッドがハンドルのロックを持っている（別のスレッドからは取れない）。
            probe = {}
            thread = threading.Thread(
                target=lambda: probe.setdefault("got", handle.lock.acquire(timeout=0.01)))
            thread.start()
            thread.join(_TIMEOUT)
            if probe.get("got"):
                handle.lock.release()
            state["lock_held"] = not probe.get("got")
        assert state["release"].wait(_TIMEOUT), "テストが Event を放さなかった"
        if state["error"] is not None:
            raise state["error"]
        return json.dumps({"status": "success", "spec": spec, "wait_seconds": wait_seconds})

    assert handle_dispatch.install(app) == ["dataset_homolog_check"]
    try:
        yield app, state
    finally:
        state["release"].set()


def _call(app, name, args) -> dict:
    content = asyncio.run(app.call_tool(name, args))
    blocks = content[0] if isinstance(content, tuple) else content
    return json.loads(blocks[0].text)


# ---------- wait_seconds の丸め ----------


@pytest.mark.parametrize("value, used, noted", [
    (0, 0.0, False), (5, 5.0, False), (30, 30.0, False), (2.5, 2.5, False),
    (-1, 0.0, True), (-math.inf, 0.0, True), (31, 30.0, True), (60, 30.0, True),
    (math.inf, 30.0, True)])
def test_clamp_wait_clamps_numbers(value, used, noted):
    got, note = pending_loads.clamp_wait(value)
    assert got == used
    assert (note is not None) is noted
    if noted:
        assert f"{used:g} 秒に丸めた" in note


@pytest.mark.parametrize("value", [math.nan, True, False, "20", None, [1]])
def test_clamp_wait_rejects_non_numbers(value):
    with pytest.raises(DomainError) as info:
        pending_loads.clamp_wait(value)
    assert info.value.code == "DATASET_BAD_REQUEST"
    assert info.value.details == {"field": "wait_seconds"}


def test_with_wait_note_adds_to_json_objects_only():
    assert json.loads(pending_loads.with_wait_note('{"a": 1}', "丸めた")) == {
        "a": 1, "wait_seconds_note": "丸めた"}
    assert pending_loads.with_wait_note('{"a": 1}', None) == '{"a": 1}'
    assert pending_loads.with_wait_note("[1, 2]", "丸めた") == "[1, 2]"
    assert pending_loads.with_wait_note("not json", "丸めた") == "not json"
    blocks = [object()]
    assert pending_loads.with_wait_note(blocks, "丸めた") is blocks


# ---------- 包み（_wrap_job） ----------


def test_fast_job_returns_the_result_and_is_not_remembered(fresh, mini):
    app, state = mini
    state["release"].set()
    out = _call(app, "dataset_homolog_check", {"spec": "s", "wait_seconds": 5})
    assert out == {"status": "success", "spec": "s", "wait_seconds": 5.0}
    assert fresh._entries == {}


def test_slow_job_returns_job_id_then_the_result(fresh, mini):
    app, state = mini
    out = _call(app, "dataset_homolog_check", {"spec": "s", "wait_seconds": 0.05})
    assert out["status"] == "running" and out["tool"] == "dataset_homolog_check"
    job_id = out["job_id"]
    assert job_id.startswith("job_")
    assert out["next"] == f'job_status(job_id="{job_id}", wait_seconds=30)'

    running = json.loads(handle_tools.job_status(job_id))
    assert running["status"] == "running" and running["job_id"] == job_id
    assert running["tool"] == "dataset_homolog_check" and running["elapsed_s"] >= 0

    state["release"].set()
    done = json.loads(handle_tools.job_status(job_id, wait_seconds=_TIMEOUT))
    assert done == {"status": "success", "spec": "s", "wait_seconds": 0.05}
    # 番号を渡せば何度でも同じ結果。省いた問い合わせからは外れる（結果を渡した）。
    assert json.loads(handle_tools.job_status(job_id)) == done
    missing = json.loads(handle_tools.job_status())
    assert missing["error"]["code"] == "missing_state" and missing["error"]["state"] == "job"


def test_job_wait_is_clamped_with_a_note(fresh, mini):
    app, state = mini
    state["release"].set()
    out = _call(app, "dataset_homolog_check", {"wait_seconds": 45})
    assert out["status"] == "success" and out["wait_seconds"] == 30.0
    assert "30 秒に丸めた" in out["wait_seconds_note"]


def test_job_bad_wait_is_bad_request_and_does_not_start(fresh, mini):
    app, state = mini
    out = _call(app, "dataset_homolog_check", {"wait_seconds": math.nan})   # 型の検査は MCP 層が先にする
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert not state["started"].is_set()
    assert fresh._entries == {}


def test_job_failure_is_an_envelope(fresh, mini):
    app, state = mini
    state["error"] = RuntimeError("secret path C:/lab/x")
    out = _call(app, "dataset_homolog_check", {"wait_seconds": 0.05})
    state["release"].set()
    done = json.loads(handle_tools.job_status(out["job_id"], wait_seconds=_TIMEOUT))
    assert done["error"]["code"] == "JOB_FAILED"
    assert "secret" not in json.dumps(done, ensure_ascii=False)   # 例外の文言は出さない


def test_job_holds_the_current_handle_lock_until_it_finishes(fresh, mini):
    app, state = mini
    ds = DatasetState()
    session_state.session.dataset = ds
    handle = handles.registry.peek(session_state.session.current_handle_id)
    out = _call(app, "dataset_homolog_check", {"wait_seconds": 0.05})
    assert out["status"] == "running"
    assert state["started"].wait(_TIMEOUT)
    assert not handle.lock.acquire(timeout=0.05)   # 走っている間はロックを取れない
    state["release"].set()
    assert json.loads(handle_tools.job_status(out["job_id"], wait_seconds=_TIMEOUT))[
        "status"] == "success"
    assert state["lock_held"] is True
    assert handle.lock.acquire(timeout=_TIMEOUT)    # 終われば放す
    handle.lock.release()


# ---------- job_status の受付番号の省き方 ----------


def test_job_status_without_jobs_is_missing_state(fresh):
    out = json.loads(handle_tools.job_status())
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["state"] == "job"
    assert out["error"]["required_tools"] == ["pathway_project", "dataset_homolog_check"]


def test_job_status_unknown_id_is_not_found(fresh):
    out = json.loads(handle_tools.job_status("job_000000000000"))
    assert out["error"]["code"] == "JOB_NOT_FOUND"


def test_job_status_picks_the_only_open_job(fresh, mini):
    app, state = mini
    job_id = _call(app, "dataset_homolog_check", {"wait_seconds": 0.05})["job_id"]
    out = json.loads(handle_tools.job_status())
    assert out["status"] == "running" and out["job_id"] == job_id
    state["release"].set()
    assert json.loads(handle_tools.job_status(wait_seconds=_TIMEOUT))["status"] == "success"


def test_job_status_asks_when_several_jobs_are_open(fresh):
    gate = threading.Event()
    ids = [fresh.submit(lambda: (gate.wait(_TIMEOUT), "{}")[1], tool="pathway_project",
                        label=f"dh_{i}").job_id for i in range(2)]
    try:
        out = json.loads(handle_tools.job_status())
        body = out["needs_input"]
        assert body["argument"] == "job_id" and body["retry_tools"] == ["job_status"]
        assert [o["value"] for o in body["options"]] == ids
        assert all("pathway_project" in o["label"] for o in body["options"])
    finally:
        gate.set()


@pytest.mark.parametrize("value", [math.nan, True, "1"])
def test_job_status_bad_wait_is_bad_request(fresh, mini, value):
    app, state = mini
    job_id = _call(app, "dataset_homolog_check", {"wait_seconds": 0.05})["job_id"]
    out = json.loads(handle_tools.job_status(job_id, wait_seconds=value))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    state["release"].set()


# ---------- 公開面 ----------


def test_long_running_tools_are_wrapped_in_the_server():
    assert handle_dispatch.LONG_RUNNING == frozenset({"pathway_project", "dataset_homolog_check"})
    tools = server.mcp._tool_manager._tools
    for name in handle_dispatch.LONG_RUNNING:
        assert tools[name].is_async, name
        assert "wait_seconds" in tools[name].parameters["properties"], name
    assert tools["job_status"].is_async
    assert tools["job_status"].annotations.readOnlyHint is True
