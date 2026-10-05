"""`.dbs` / `.lbm2` をテスト自身が組むためのヘルパ（実データに依存しない規約）。"""
from __future__ import annotations

import struct
import zipfile

import lz4.block
import msgpack


def pack_chunk(records: list[list]) -> bytes:
    """LargeListMessagePack のチャンク 1 個を組む。

    要素数に応じて array ヘッダが短形式になる点を再現する（本番の罠）。
    要素は常にオフセット 5 から始まる。
    """
    body = b"".join(msgpack.packb(r, use_bin_type=True) for r in records)
    n = len(records)
    if n < 16:
        header = bytes([0x90 | n])
    elif n < 65536:
        header = b"\xdc" + struct.pack(">H", n)
    else:
        header = b"\xdd" + struct.pack(">I", n)
    raw = header + b"\x00" * (5 - len(header)) + body
    comp = lz4.block.compress(raw, store_size=False)
    return (b"\xc9" + struct.pack(">I", len(comp) + 5) + b"\x63"
            + b"\xd2" + struct.pack(">i", len(raw)) + comp)


def record(name="GABA", mz=104.0706, ion_mode=1, peaks=((87.04, 999.0), (69.03, 500.0)),
           inchikey="BTCSSZJGUNDROE-UHFFFAOYSA-N"):
    """MoleculeMsReference 1 件（Key 0..28）。"""
    r = [None] * 29
    r[0] = 0
    r[1] = mz
    r[2] = [[1, [1.23, 0, 0]]]
    r[3] = ion_mode
    r[4] = [[m, i] + [None] * 4 + [0, 0] + [None] * 3 + [0, False] for m, i in peaks]
    r[5] = name
    r[6] = None
    r[7] = ""
    r[8] = "NCCCC(=O)O"
    r[9] = inchikey
    r[10] = [1.00782503207, 1, "[M+H]+", 1, 1, True, 0.0, 0.0, False, False]
    r[14] = "AminoAcid"
    return r


#: 注釈器の Parameter（Key 0..19）。値は read_storage_meta の検索パラメータ読み取り用。
PARAMETER = [0.0, 2000.0, 2.0, 100.0, 20.0, 0.01, 0.025, 0.0, 0.0, 0.0225,
             0.0225, 0.09, 0.0, 0.8, 1.0, True, True, False, False, 0.1]


def _wrap_storage(storage: dict) -> bytes:
    packed = msgpack.packb(storage, use_bin_type=True)
    comp = lz4.block.compress(packed, store_size=False)
    return (b"\xc9" + struct.pack(">I", len(comp) + 5) + b"\x63"
            + b"\xd2" + struct.pack(">i", len(packed)) + comp)


def write_dbs(path, databases, *, parameter=PARAMETER):
    """`databases` = [(ライブラリ名, レコード列, [注釈器の Key, ...]), ...] から `.dbs` を組む。

    MS-DIAL の `DataBaseStorage` と同じく、`Storage` の `MetabolomicsDataBases` に
    ライブラリごとの `DataBase`（名前・元パス）と注釈器の組（`Pairs`）を並べ、
    実体は `MetabolomicsDB/<名前>/DataBase` に置く。
    """
    storage = {"MetabolomicsDataBases": [
        {"DataBase": [name, 4, 2, f"C:/x/{name}.lbm2"],
         "Pairs": [[0, {"SerializableAnnotatorKey":
                        [3, {"Parameter": parameter, "SourceType": 4, "Key": key,
                             "Priority": 1}]}] for key in keys]}
        for name, _records, keys in databases]}
    with zipfile.ZipFile(path, "w") as z:
        for name, records, _keys in databases:
            z.writestr(f"MetabolomicsDB/{name}/DataBase", pack_chunk(records))
        z.writestr("Storage", _wrap_storage(storage))
    return path
