"""データセットへの呼び出しをワーカースレッドへ退避する（spec 2026-10-22 §4 の 1・
spec 2026-11-06 §5）。

FastMCP は同期関数のツールをイベントループのスレッドでそのまま呼ぶ。数分かかる読み込みが
1 本走ると、ほかの呼び出しはすべて待たされる。そこで、引数に `dataset_handle` を持つツールの
呼び出しをスレッドで回し、同じハンドルへの呼び出しはそのハンドルのロックで直列化する。

- `dataset_handle` を渡した呼び出し: そのハンドルのロックの下でワーカースレッドに回す。
- 渡さない呼び出し: 対話のセッションの直近のハンドル（`session.current_handle_id`）があれば、
  **そのロックの下で**ワーカースレッドに回す（同じデータセットへの対話と再解析の呼び出しが
  重ならない）。呼び出しの間は `session_state.pinned_handle` でそのハンドルに固定するので、
  途中で `load_dataset` が直近を切り替えても、ロックを取っていない別のハンドルの状態を
  触らない。ロックを待つ間にそのハンドルが閉じられたら、ツールを呼ばずに `missing_state`
  を返す。直近のハンドルが無ければ今までどおりイベントループのスレッドで同期に呼ぶ。
- `SESSION_READERS`（`dataset_handle` を持たないが `session.dataset` / `session.library` を
  読むツール）も、渡さない呼び出しと同じ道を通す。イベントループでロックを取らずに読むと、
  ワーカーで動く同じデータセットへの呼び出しが行列を書き換えたり、`library_load` が
  store を閉じたりする最中に読んでしまう。

- `LONG_RUNNING`（長くかかりうる道具。`wait_seconds` を持つ）は、`load_dataset` の受付と同じ形で
  回す: 呼び出しを受付（`core.pending_loads.jobs`）のスレッドで始め（ハンドルのロックと固定はその
  スレッドが取り、仕事が終わるまで持つ）、`wait_seconds`（既定 20 秒。0〜30 秒の外は丸めて
  `wait_seconds_note` で知らせる）の中で終われば戻り値をそのまま返し、終わらなければ
  `{"status": "running", "job_id": "job_…"}` を返して仕事は続ける。結果は `job_status(job_id)` で引く。

`install(mcp)` は全ツールの登録が済んだ後（server.py の末尾）で 1 回呼ぶ。包むのは
引数に `dataset_handle` を持つツールと `SESSION_READERS` と `ALWAYS_THREADED`（グローバルな
セッションに触れない重いツール）と `LONG_RUNNING` だけ。包んだ関数は `tool.fn` を差し替えるだけで、module の
関数はそのまま——ツールの中から別のツールの関数を直に呼んでも、同じスレッドでそのまま走る
（もう一度スレッドへ回したりロックを取り直したりしない）。
"""
from __future__ import annotations

import functools
import inspect

import anyio

from metabolomix.core import handles, pending_loads, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload

__all__ = ["ALWAYS_THREADED", "LONG_RUNNING", "SESSION_READERS", "install"]

#: ハンドルの有無に依らずスレッドで回すツール。読み込みは数分かかりうるので、イベントループ
#: から外す。`load_dataset` は読み込みを受付（`core.pending_loads`）のスレッドで始め、
#: ワーカースレッドで最長 `wait_seconds` 待つ。セッションへ書くのは、読み込みが成功して
#: 終わるときの直近の切り替え（受付のスレッドの `session.adopt_handle`。`_current_lock` の下）
#: だけ。`dataset_load_status` は受付簿のロックを取って読み、最長 `wait_seconds`（30 秒まで）
#: 読み込みの終わりを待つので、イベントループを止めないようにこちらへ載せる。
ALWAYS_THREADED = frozenset({
    "load_dataset",
    "dataset_load_status",
    "job_status",  # 受付簿のロックを取って読み、最長 `wait_seconds`（30 秒まで）待つ
    "dataset_merge_polarities",  # 2 つのハンドルのロックを自分で取る（最大 LOCK_TIMEOUT 待つ）
})

#: 長くかかりうるので、`wait_seconds` の中で終わらなければ受付番号（`job_…`）を返して続けるツール
#: （`load_dataset` の受付と同じ形。どれも引数に `wait_seconds` を持つ）。
LONG_RUNNING = frozenset({
    "pathway_project",        # MSLipidMapper（研究室サーバ → ローカルの R。最悪で 600 秒を超える）
    "dataset_homolog_check",  # auto は全特徴の名前を GOSLIN で読む（分のオーダー）
})

#: `dataset_handle` を持たないが、直近のハンドルの `session.dataset` / `session.library` を
#: （自身か、ヘルパー経由で）読むツール。直近のハンドルがあれば、そのロックと固定の下で
#: ワーカースレッドに回す。自身のソースで読むものは
#: `tests/test_session_handle_unification.py` の検査が捕まえる（ヘルパー経由は手で載せる）。
SESSION_READERS: frozenset[str] = frozenset()  # 今は無い（全部が `dataset_handle` を持つ）


def _call_locked(fn, kwargs: dict):
    handle = handles.registry.peek(kwargs.get("dataset_handle"))
    if handle is None:
        # 未知のハンドルはロックを取らずに呼ぶ（ツール側が封筒で答える）。
        return fn(**kwargs)
    with handle.lock:
        return fn(**kwargs)


def _call_pinned(fn, kwargs: dict, handle):
    """`dataset_handle` を省いた呼び出しを、直近だったハンドルのロックと固定の下で呼ぶ。

    ロックを待つ間にそのハンドルが閉じられていれば、ツールを呼ばずに `missing_state` を返す
    （閉じたハンドルの代わりに、直近になった別のハンドルやセッション自身のライブラリを
    黙って使わない。何もしていないので、呼び直せば今の直近で動く）。
    """
    with handle.lock:
        if handles.registry.peek(handle.handle_id) is not handle:
            return missing_state(
                "dataset", ["load_dataset"],
                f"この呼び出しが動く前に、直近のデータセットのハンドル {handle.handle_id} が"
                "閉じられました。何もしていません（別のハンドルへは落としていません）。"
                "load_dataset で読み直すか、呼び直してください。")
        with session_state.pinned_handle(handle.handle_id):
            return fn(**kwargs)


def _current_handle():
    """対話のセッションの直近のハンドル（無ければ None）。登録簿は module 修飾で読む。"""
    handle_id = getattr(session_state.session, "current_handle_id", None)
    return handles.registry.peek(handle_id) if handle_id is not None else None


def _wrap(fn, *, always: bool):
    @functools.wraps(fn)
    async def runner(**kwargs):
        if always or kwargs.get("dataset_handle") is not None:
            return await anyio.to_thread.run_sync(functools.partial(_call_locked, fn, kwargs))
        handle = _current_handle()
        if handle is None:
            return fn(**kwargs)
        return await anyio.to_thread.run_sync(functools.partial(_call_pinned, fn, kwargs, handle))

    return runner


def _job_call(fn, kwargs: dict):
    """受付のスレッドで回す呼び出し（ハンドルのロックと固定は、ここで取って仕事が終わるまで持つ）。

    `dataset_handle` を渡した呼び出しはそのハンドルのロックの下、省いた呼び出しは受け付けたときの
    直近のハンドルのロックと固定の下（`_call_pinned`。待つ間に閉じられたら `missing_state`）、
    直近が無ければそのまま呼ぶ。
    """
    if kwargs.get("dataset_handle") is not None:
        return functools.partial(_call_locked, fn, kwargs)
    handle = _current_handle()
    if handle is None:
        return functools.partial(fn, **kwargs)
    return functools.partial(_call_pinned, fn, kwargs, handle)


def _wrap_job(fn, name: str):
    """`LONG_RUNNING` のツールを、受付番号（`job_…`）を返せる形で包む。"""
    default_wait = inspect.signature(fn).parameters["wait_seconds"].default

    @functools.wraps(fn)
    async def runner(**kwargs):
        try:
            wait, note = pending_loads.clamp_wait(kwargs.get("wait_seconds", default_wait))
        except DomainError as exc:
            return mztab_error(exc.code, exc.message, exc.details or None)
        kwargs = {**kwargs, "wait_seconds": wait}
        handle_id = kwargs.get("dataset_handle") or getattr(
            session_state.session, "current_handle_id", None)
        try:
            entry = pending_loads.jobs.submit(_job_call(fn, kwargs), tool=name,
                                              label=handle_id)
        except DomainError as exc:  # 走っている受付が上限（JOB_LIMIT）
            return mztab_error(exc.code, exc.message, exc.details or None)
        payload = await anyio.to_thread.run_sync(entry.wait, wait)
        if payload is not None:
            pending_loads.jobs.discard(entry.job_id)  # 受付番号を誰にも渡していないので覚えない
            return pending_loads.with_wait_note(payload, note)
        return json_payload({
            "status": "running", "job_id": entry.job_id, "tool": name,
            **({"wait_seconds_note": note} if note else {}),
            "note": ("仕事は続けている（ハンドルのロックは仕事が終わるまでこの仕事が持つ）。"
                     "結果は job_status で引く。"),
            "next": f'job_status(job_id="{entry.job_id}", wait_seconds=30)'})

    return runner


def install(mcp) -> list[str]:
    """`dataset_handle` を持つ同期ツールと `SESSION_READERS` と `LONG_RUNNING` を包む。包んだツール名を返す
    （2 回目は空）。"""
    wrapped = []
    for name, tool in mcp._tool_manager._tools.items():
        if tool.is_async:  # 元から async のツールと、包み済みのツール
            continue
        if name in LONG_RUNNING and "wait_seconds" in inspect.signature(tool.fn).parameters:
            tool.fn = _wrap_job(tool.fn, name)
            tool.is_async = True
            wrapped.append(name)
            continue
        always = name in ALWAYS_THREADED
        if (not always and name not in SESSION_READERS
                and "dataset_handle" not in inspect.signature(tool.fn).parameters):
            continue
        tool.fn = _wrap(tool.fn, always=always)
        tool.is_async = True
        wrapped.append(name)
    return sorted(wrapped)
