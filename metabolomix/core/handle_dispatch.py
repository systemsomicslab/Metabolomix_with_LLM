"""ハンドルつきの呼び出しをワーカースレッドへ退避する（spec 2026-10-22 §4 の 1）。

FastMCP は同期関数のツールをイベントループのスレッドでそのまま呼ぶ。数分かかる読み込みが
1 本走ると、ほかの呼び出しはすべて待たされる。並行の再解析のエージェントは別々の
ハンドルで動くので、**`dataset_handle` を渡した呼び出しだけ**をスレッドで回し、
同じハンドルへの呼び出しはそのハンドルのロックで直列化する。

`dataset_handle` を渡さない呼び出しは今までどおりイベントループのスレッドで同期に呼ぶ
（対話のセッションは 1 つの可変状態なので、並行にすると壊れる。挙動を変えない）。

`install(mcp)` は全ツールの登録が済んだ後（server.py の末尾）で 1 回呼ぶ。包むのは
引数に `dataset_handle` を持つツールと `ALWAYS_THREADED`（グローバルなセッションに
触れない重いツール）だけ。
"""
from __future__ import annotations

import functools
import inspect

import anyio

from metabolomix.core import handles

__all__ = ["ALWAYS_THREADED", "install"]

#: ハンドルの有無に依らずスレッドで回すツール（対話のセッションを読み書きしないもの）。
ALWAYS_THREADED = frozenset({"dataset_handle_open"})


def _call_locked(fn, kwargs: dict):
    handle = handles.registry.peek(kwargs.get("dataset_handle"))
    if handle is None:
        # 未知のハンドルはロックを取らずに呼ぶ（ツール側が封筒で答える）。
        return fn(**kwargs)
    with handle.lock:
        return fn(**kwargs)


def _wrap(fn, *, always: bool):
    @functools.wraps(fn)
    async def runner(**kwargs):
        if not always and kwargs.get("dataset_handle") is None:
            return fn(**kwargs)
        return await anyio.to_thread.run_sync(functools.partial(_call_locked, fn, kwargs))

    return runner


def install(mcp) -> list[str]:
    """`dataset_handle` を持つ同期ツールを包む。包んだツール名を返す（2 回目は空）。"""
    wrapped = []
    for name, tool in mcp._tool_manager._tools.items():
        if tool.is_async:  # 元から async のツールと、包み済みのツール
            continue
        always = name in ALWAYS_THREADED
        if not always and "dataset_handle" not in inspect.signature(tool.fn).parameters:
            continue
        tool.fn = _wrap(tool.fn, always=always)
        tool.is_async = True
        wrapped.append(name)
    return sorted(wrapped)
