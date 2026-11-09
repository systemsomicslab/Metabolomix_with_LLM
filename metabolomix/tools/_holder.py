"""ツールが状態を読む先（対話のセッションかハンドル）を決める共通ヘルパ。

`dataset_handle=None` なら `session_state.session`、そうでなければそのハンドル
（`core.handles`）。開いていないハンドルは `DATASET_HANDLE_NOT_FOUND` の封筒を返し、
対話のセッションへ黙って落とさない。
"""
from __future__ import annotations

import json

from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import missing_state, mztab_error

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

    明示（`arf2_path` / `dcl_file`）があればそれを使い、無いものはデータセット（ハンドル、
    省略時は直近のハンドル）の `.arf2` / `.dcl` を使う。**データフォルダは探さない**
    （B-2b。探すと別のデータセットの同じ spot_id のスペクトルを黙って読む——`_spot_query` は
    `.arf2` と `.dcl` の一致しか確かめない）。要るものが明示で足りていればデータセットが
    無くてもよい。データセットが無い（直近が無い・閉じた）なら `missing_state("dataset")`、
    データセットにそのファイルが無ければ `no_alignment_files`。
    """
    from metabolomix.dataset.artifacts import artifact_path

    wanted = {"arf2": arf2_path, "dcl": dcl_file}
    if all(wanted[name] for name in need):
        return arf2_path, dcl_file, None
    ds = getattr(holder, "dataset", None)
    if ds is None:
        return arf2_path, dcl_file, json.loads(missing_state(
            "dataset", ["load_dataset"],
            "先に load_dataset を実行してください（スポットの照合に使う .arf2 / .dcl は、"
            "直近のデータセットのものを使います。データフォルダは探しません）。"))
    if not arf2_path:
        found = getattr(holder, "arf2_path", None) or artifact_path(ds, "arf2")
        arf2_path = str(found) if found else None
    if not dcl_file:
        found = getattr(holder, "dcl_path", None) or artifact_path(ds, "dcl")
        dcl_file = str(found) if found else None
    lacking = [name for name, value in (("arf2", arf2_path), ("dcl", dcl_file))
               if name in need and not value]
    if not lacking:
        return arf2_path, dcl_file, None
    return arf2_path, dcl_file, {
        "status": "no_alignment_files",
        "message": ("このデータセットにはスポットの照合に要る "
                    + " / ".join("." + name for name in lacking)
                    + " が見つかりません（load_dataset の arf2_path で指定するか、"
                    "arf2_path / dcl_file を渡してください）。データフォルダは探しません。"),
    }
