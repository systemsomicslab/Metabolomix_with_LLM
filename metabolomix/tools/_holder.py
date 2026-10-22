"""ツールが状態を読む先（対話のセッションかハンドル）を決める共通ヘルパ。

`dataset_handle=None` なら `session_state.session`、そうでなければそのハンドル
（`core.handles`）。開いていないハンドルは `DATASET_HANDLE_NOT_FOUND` の封筒を返し、
対話のセッションへ黙って落とさない。
"""
from __future__ import annotations

from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error

__all__ = ["resolve_holder"]


def resolve_holder(dataset_handle: str | None):
    """`.dataset` / `.library` を持つもの、または封筒（str）。"""
    try:
        return handles.holder_for(dataset_handle)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
