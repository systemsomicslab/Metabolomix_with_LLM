"""キュレーション用の合成入力（.arf2 / .dcl / .EIC.aef / .arf / .msp）。

実データはリポジトリに含めないので、各形式のレイアウトをここで組み立てる。
Key 番号の正準は docs/schema/*.md（MsScanMatchResult.md / AlignmentSpotProperty.md）。
"""
from __future__ import annotations

from pathlib import Path

import lz4.block
import msgpack


def at_keys(size: int, values: dict) -> list:
    row = [None] * size
    for key, value in values.items():
        row[key] = value
    return row


def chromxs(rt: float, mz: float) -> list:
    """ChromXs の入れ子表現（1=RT, 3=m/z）。"""
    return [[1, [rt]], [3, [mz]]]


def match_result(overrides=None, **kwargs) -> list:
    """MsScanMatchResult の生配列（39 要素）。既定は「参照一致した MS/MS あり」。

    overrides: dict with numeric keys, or use **kwargs for compatibility.
    """
    if overrides is None:
        overrides = kwargs
    values = {
        0: "PC 34:1", 1: "KEY-PC341", 2: 3.2,
        3: 0.81, 4: 0.64, 5: 0.9, 6: 12.0, 7: 0.7, 8: 0.0,
        9: 0.99, 10: 0.0, 11: 0.0, 12: -1.0, 13: 0.98,
        14: 5, 15: True, 16: True, 17: True, 18: False,
        19: True, 20: True, 21: False, 22: False, 23: False,
        24: -1, 26: 4, 27: "Msp20260101000000_lib_1", 28: -1,
        29: 0.0, 30: False, 31: 1, 32: 0.0, 33: True, 34: False,
        35: False, 36: 0.0, 37: float("nan"), 38: -1.0,
    }
    values.update(overrides)
    return at_keys(39, values)


def arf2_spot_raw(*, spot_id: int, name: str = "PC 34:1", mz: float = 760.5851,
                  rt: float = 12.0, ontology: str = "PC", adduct: str = "[M+H]+",
                  formula: str = "C42H82NO8P", ion_mode: int = 0,
                  representative_file_id: int = 0,
                  matches: list | None = None) -> list:
    """AlignmentSpotProperty の生配列（Key 0..59）。matches は match_result() の並び。"""
    values = {
        0: spot_id, 1: spot_id, 3: representative_file_id,
        4: chromxs(rt, mz), 5: mz, 11: ion_mode, 12: name,
        13: [formula, 0.0], 14: ontology, 15: "", 16: "",
        31: 10000.0, 32: 100.0, 33: 20000.0, 34: 0.2,
        35: 30.0, 36: 50.0, 37: 10.0, 43: mz - 0.001, 44: mz + 0.001,
        49: 1.0, 51: 1.0, 54: [0.0, 1, adduct],
        56: [list(matches or []), {}, []],
        59: -1,
    }
    return at_keys(60, values)


def write_arf2(path: Path, spots: list[list]) -> Path:
    """MS-DIAL の外側コンテナ msgpack([header, msgpack(非圧縮長) + LZ4]) で書く。"""
    inner = b"".join(msgpack.packb(spot, use_bin_type=True) for spot in spots)
    payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
        inner, store_size=False)
    path = Path(path)
    path.write_bytes(msgpack.packb(["hdr", payload], use_bin_type=True))
    return path
