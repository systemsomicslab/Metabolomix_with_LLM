"""MS-DIAL の解析フォルダを tmp に作るテスト用 fixture（後続の Task が足していく）。"""
from __future__ import annotations

from pathlib import Path

import lz4.block
import msgpack


def pack_msdial(value) -> bytes:
    unpacked = msgpack.packb(value, use_bin_type=True)
    size_header = msgpack.packb(len(unpacked))
    compressed = lz4.block.compress(unpacked, store_size=False)
    extension = msgpack.ExtType(99, size_header + compressed)
    return msgpack.packb(extension, use_bin_type=True)


def analysis_file_row(
    file_id: int,
    name: str,
    class_id: str,
    *,
    file_type: int = 0,
    order: int | None = None,
    batch: int = 1,
) -> list:
    """AnalysisFileBean 1 行（Key0〜11。docs/schema/AnalysisFileBean.md）。"""
    return [
        f"C:/data/{name}.wiff",  # Key0
        name,  # Key1
        file_type,  # Key2
        class_id,  # Key3
        file_id + 1 if order is None else order,  # Key4
        file_id,  # Key5
        True,  # Key6
        "", [], "", "",  # Key7〜10（未使用）
        batch,  # Key11
    ]


def write_mddata(path, rows: list[list]) -> Path:
    path = Path(path)
    path.write_bytes(pack_msdial([rows, [], [], [], [], None, None, None, None]))
    return path
