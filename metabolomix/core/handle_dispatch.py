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

`install(mcp)` は全ツールの登録が済んだ後（server.py の末尾）で 1 回呼ぶ。包むのは
引数に `dataset_handle` を持つツールと `SESSION_READERS` と `ALWAYS_THREADED`（グローバルな
セッションに触れない重いツール）だけ。包んだ関数は `tool.fn` を差し替えるだけで、module の
関数はそのまま——ツールの中から別のツールの関数を直に呼んでも、同じスレッドでそのまま走る
（もう一度スレッドへ回したりロックを取り直したりしない）。
"""
from __future__ import annotations

import functools
import inspect

import anyio

from metabolomix.core import handles, session_state
from metabolomix.core.mcp_errors import missing_state

__all__ = ["ALWAYS_THREADED", "SESSION_READERS", "install"]

#: ハンドルの有無に依らずスレッドで回すツール。読み込みは数分かかりうるので、イベントループ
#: から外す。`load_dataset` がセッションへ書くのは成功時の直近の切り替え
#: （`session.adopt_handle`。`_current_lock` の下）だけ。
ALWAYS_THREADED = frozenset({"load_dataset"})

#: `dataset_handle` を持たないが、直近のハンドルの `session.dataset` / `session.library` を
#: （自身か、ヘルパー経由で）読むツール。直近のハンドルがあれば、そのロックと固定の下で
#: ワーカースレッドに回す。自身のソースで読むものは
#: `tests/test_session_handle_unification.py` の検査が捕まえる（ヘルパー経由は手で載せる）。
SESSION_READERS = frozenset({
    "curation_review",          # session.library.store
    "curation_suggest",         # session.library.store / source_path
    "pathway_project",          # _build_dataset → session.dataset
    "verify_peak_annotation",   # peak_verification.msms_evidence → session.library.store
    "save_figure",              # reports._pca_candidates / _differential_candidates → session.dataset
    "plot_pca_loadings",        # pca_loadings_tools._candidates → session.dataset
})


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


def install(mcp) -> list[str]:
    """`dataset_handle` を持つ同期ツールと `SESSION_READERS` を包む。包んだツール名を返す
    （2 回目は空）。"""
    wrapped = []
    for name, tool in mcp._tool_manager._tools.items():
        if tool.is_async:  # 元から async のツールと、包み済みのツール
            continue
        always = name in ALWAYS_THREADED
        if (not always and name not in SESSION_READERS
                and "dataset_handle" not in inspect.signature(tool.fn).parameters):
            continue
        tool.fn = _wrap(tool.fn, always=always)
        tool.is_async = True
        wrapped.append(name)
    return sorted(wrapped)
