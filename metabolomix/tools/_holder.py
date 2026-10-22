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


def alignment_paths(holder, dataset_handle: str | None, arf2_path: str | None,
                    dcl_file: str | None, *, need: tuple = ("arf2", "dcl")):
    """スポットの照合に使う `.arf2` / `.dcl` と、足りないときの封筒（dict。足りていれば None）。

    ハンドルなら、明示が無ければハンドルの `.arf2` / `.dcl` を使う。**ハンドルに無いときは
    対話のデータフォルダへ落ちない**——落ちると、別のデータセットの同じ spot_id の
    スペクトルを黙って読む（`_spot_query` は `.arf2` と `.dcl` の一致しか確かめない）。
    対話のセッション（`dataset_handle=None`）では渡されたまま返す（従来どおりデータフォルダを探す）。
    """
    if dataset_handle is None:
        return arf2_path, dcl_file, None
    arf2_path = arf2_path or holder.arf2_path
    dcl_file = dcl_file or holder.dcl_path
    lacking = [name for name, value in (("arf2", arf2_path), ("dcl", dcl_file))
               if name in need and not value]
    if not lacking:
        return arf2_path, dcl_file, None
    return arf2_path, dcl_file, {
        "status": "no_alignment_files",
        "message": ("このハンドルにはスポットの照合に要る "
                    + " / ".join("." + name for name in lacking)
                    + " が見つかりません（dataset_handle_open の arf2_path で指定するか、"
                    "arf2_path / dcl_file を渡してください）。対話のデータフォルダは探しません。"),
    }
