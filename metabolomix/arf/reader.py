import argparse
import io
import os
from enum import Enum
import sys
import re
from typing import Iterator

import lz4.block
import msgpack
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 非インタラクティブバックエンドを使用
import matplotlib.pyplot as plt
import pprint

from metabolomix.msdial.classes import get_sample_class_id
from metabolomix.msdial.tags import get_sample_peak_tag_info


class IonMode(Enum):
    Positive = 0
    Negative = 1
    Both = 2
    Unknown = -1

    @classmethod
    def from_int(cls, val):
        try:
            return cls(val).name
        except Exception:
            return cls.Unknown.name


def _convert_to_times(data: list) -> dict:
    result = {}
    for t in data[:4]:
        if t[0] == 1 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["rt"] = t[1][0]
        elif t[0] == 2 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["ri"] = t[1][0]
        elif t[0] == 3 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["m/z"] = t[1][0]
        elif t[0] == 4 and isinstance(t[1], (list, tuple)) and len(t[1]) > 0:
            result["dt"] = t[1][0]
    return result


def _decode(value):
    return value.decode('utf-8', errors='ignore') if isinstance(value, bytes) else value


def _convert_arf_peak_group(
    group: list,
    group_index: int | None = None,
    source_block_index: int | None = None,
    source_local_index: int | None = None,
) -> dict | None:
    if not _is_arf_peak_group(group):
        return None

    representative = group[0]
    times = _convert_to_times(representative[15]) if len(representative) > 15 and isinstance(representative[15], list) else {}

    return {
        "MasterAlignmentID": group_index,
        "AlignmentID": group_index,
        "RT": times.get("rt", None),
        "MassCenter": representative[22] if len(representative) > 22 else None,
        "IonMode": IonMode.from_int(representative[23]) if len(representative) > 23 else IonMode.Unknown.name,
        "Name": _decode(representative[24]) if len(representative) > 24 else None,
        "HeightAverage": representative[18] if len(representative) > 18 else None,
        "AlignedPeakProperties": group,
        "SourceBlockIndex": source_block_index,
        "SourceLocalIndex": source_local_index,
    }

def _msgpack_stream(data: bytes):
    return msgpack.Unpacker(
        io.BytesIO(data),
        raw=False,
        # POSモードの .arf はヘッダ相当フィールドに不正なUTF-8バイトを含む
        # ことがある（NEGは含まない）。strict だと1バイトで全体が
        # UnicodeDecodeError で停止するため、既存の _decode(errors='ignore')
        # 規約に合わせて不正バイトを無視して読み進める。
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )


def _decompress_lz4_payload(compressed_data: bytes) -> bytes:
    stream = io.BytesIO(compressed_data)
    unpacker = msgpack.Unpacker(
        stream,
        raw=False,
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )
    size = next(unpacker)
    read_size = unpacker.tell()
    # スライスの写しを作らない（圧縮ペイロードの写しが 1 つ増える。展開後の大きさに比べ小さいが、
    # 1 ブロックが数百 MB の実ファイルでは無視できない）。
    return lz4.block.decompress(
        memoryview(compressed_data)[read_size:],
        uncompressed_size=size,
    )


def _decode_lz4_msgpack_payload(compressed_data: bytes) -> list:
    return list(_msgpack_stream(_decompress_lz4_payload(compressed_data)))


def _payload_from_top_level_object(obj):
    if isinstance(obj, msgpack.ExtType):
        return obj.data

    if isinstance(obj, (list, tuple)) and len(obj) == 2:
        _header, payload = obj
        if isinstance(payload, msgpack.ExtType):
            return payload.data
        if isinstance(payload, (bytes, bytearray)):
            return bytes(payload)

    return None


def _iter_lz4_msgpack_blocks(data: bytes):
    for block_index, obj in enumerate(_msgpack_stream(data)):
        payload = _payload_from_top_level_object(obj)
        if payload is None:
            print(
                f"[WARNING] ARF top-level object {block_index} is not a supported LZ4 MessagePack block; skipped",
                file=sys.stderr,
            )
            continue
        yield block_index, _decode_lz4_msgpack_payload(payload)


def _is_arf_peak_group(item) -> bool:
    """item が 1スポット分の「グループ（= サンプル行の並び）」か判定する。

    グループは空でないリストで、要素がすべてピークの行（`_is_alignment_peak_row`）。
    行数では見分けない（1〜2 試料のアラインメントは行が 1〜2 本）。ヘッダの入れ子リスト
    （ChromXs の形 `[[1.0, 0, 0], ...]` など）は行の形（25 要素超・Key18 が数）を満たさず、
    群に数えない。
    """
    return (
        isinstance(item, list)
        and len(item) > 0
        and all(_is_alignment_peak_row(row) for row in item)
    )


def _is_arf_feature_container_object(item) -> bool:
    """item が特徴量コンテナ（グループを内包する）か判定する。

    ヘッダ長は ion mode で異なる（NEG は [id, id] の2要素、POS は先頭に
    プロジェクト文字列を含む99要素などになる）ため、位置を仮定せず
    『グループ要素を1つ以上含むか』で判定する。単独グループ自身は
    グループ要素を持たないためコンテナ扱いされない。
    """
    return isinstance(item, list) and any(_is_arf_peak_group(e) for e in item)


def _top_level_stream(source):
    """トップレベルのオブジェクト（LZ4 ブロック）を順に解く Unpacker。

    source はファイル（バイナリで開いたもの）か bytes。ファイルなら全体を 1 つの bytes に
    読み込まず、ブロック 1 つ分ずつ読む（圧縮されたブロックの写しは 1 つずつしか持たない）。
    """
    if isinstance(source, (bytes, bytearray, memoryview)):
        return _msgpack_stream(bytes(source))
    return msgpack.Unpacker(
        source,
        raw=False,
        unicode_errors="ignore",  # _msgpack_stream と同じ規約
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
        read_size=1024 * 1024,
    )


def _iter_lz4_decompressed_blocks(source):
    """トップレベルのブロックを展開済みバイト列のまま返す（msgpack は解かない）。

    実ファイルは全体が 1 ブロックのことがある（185 MB → 展開 966 MB。issue #4）。
    展開後を `list(...)` で一度に Python オブジェクト化すると全行が同時に生きるので、
    解くのは呼び出し側がグループ単位で行う。LZ4 のブロックは丸ごとでしか伸長できないので、
    展開後のバイト列はブロック 1 つ分を持つ（それより細かくはできない）。
    source はファイルか bytes（`_top_level_stream`）。
    """
    for block_index, obj in enumerate(_top_level_stream(source)):
        payload = _payload_from_top_level_object(obj)
        obj = None  # 圧縮ペイロードの写しを展開後まで抱えない
        if payload is None:
            print(
                f"[WARNING] ARF top-level object {block_index} is not a supported LZ4 MessagePack block; skipped",
                file=sys.stderr,
            )
            continue
        decompressed = _decompress_lz4_payload(payload)
        payload = None
        yield block_index, decompressed


def _starts_array(first_byte: int) -> bool:
    """msgpack の型バイトが配列（fixarray / array16 / array32）か。"""
    return 0x90 <= first_byte <= 0x9F or first_byte in (0xDC, 0xDD)


def _iter_peak_groups_in_stream(decompressed: bytes):
    """展開済みストリームからグループ（1 スポット分の行の並び）を 1 個ずつ返す。

    トップレベルの各項目を `read_array_header` で開き、要素を 1 個ずつ解く。
    判定は旧来の「項目を丸ごと解いてから調べる」実装と同じ:
    - グループ要素を 1 つ以上含む項目はコンテナ。グループ要素だけを返す
      （先頭のヘッダのスカラー群は読み捨てる）。
    - グループ要素を含まず、要素が全部行（list）なら、項目自体が単独グループ
      （コンテナ外に置かれた末尾スポット等）。
    - それ以外（グループを持たないメタデータ・配列でない項目）は読み飛ばす。
    """
    unpacker = msgpack.Unpacker(
        io.BytesIO(decompressed),
        raw=False,
        # POSモードの .arf はヘッダ相当フィールドに不正なUTF-8バイトを含む
        # ことがある（_msgpack_stream と同じ規約）。
        unicode_errors="ignore",
        strict_map_key=False,
        max_buffer_size=1024 * 1024 * 1024,
    )
    end = len(decompressed)
    try:
        while unpacker.tell() < end:
            if not _starts_array(decompressed[unpacker.tell()]):
                unpacker.skip()
                continue
            length = unpacker.read_array_header()
            found_group = False
            others = []  # 単独グループの判定に要る。コンテナと分かったら溜めない
            for _ in range(length):
                element = unpacker.unpack()
                if _is_arf_peak_group(element):
                    found_group = True
                    yield element
                elif not found_group:
                    others.append(element)
            if not found_group and _is_arf_peak_group(others):
                yield others
    except msgpack.OutOfData:
        # 末尾が切れている。そこまでに読めたグループは返し済み。
        return


def _iter_arf_peak_groups(source):
    for block_index, decompressed in _iter_lz4_decompressed_blocks(source):
        local_index = 0
        for group in _iter_peak_groups_in_stream(decompressed):
            yield block_index, local_index, group
            local_index += 1


def deserialize_lz4_packed_msgpack(data: bytes) -> list:
    """LZ4圧縮されたMsgPackデータを解凍し、内部のリストを返します。"""
    all_spots = []
    for _block_index, items in _iter_lz4_msgpack_blocks(data):
        all_spots.extend(items)
    return all_spots


def extract_arf_data(
    data,
    group_index: int | None = None,
    source_block_index: int | None = None,
    source_local_index: int | None = None,
) -> dict | None:
    if _is_arf_peak_group(data):
        return _convert_arf_peak_group(
            data,
            group_index=group_index,
            source_block_index=source_block_index,
            source_local_index=source_local_index,
        )

    # .arf ファイルのみを扱うため、他の形式はサポートしない
    return None


class ArfMemoryError(MemoryError):
    """`.arf` を読み切る前に RAM が尽きたことと、そこまでの規模を運ぶ。

    `str(MemoryError())` は空文字なので、素の `MemoryError` のままだと呼び出し側は
    理由なしのエラーしか返せない。400 検体規模で最初に当たる失敗がこれなので、
    どこまで読めたか（スポット数・行数）を付けて投げ直す。
    """

    def __init__(self, *, spots: int, rows: int):
        super().__init__(f"{spots} スポット / {rows} 行を読んだところで RAM が尽きた")
        self.spots = spots
        self.rows = rows


#: `deserialize` が行（`AlignedPeakProperties` の要素）にそのまま残す Key。
#: 0〜9 は `_file_name_from_row` が先頭の文字列を探して走査し、10 は is_msms、
#: 18/20/21/22 は height / area / area_above_baseline / m_z。15（rt）と
#: 37（PeakShape）は読まれる部分だけに縮めて残す（`_slim_row`）。
#: **行から新しい Key を読むときは、ここに足すこと**——足さないと None が返る。
#: 同値性は `tests/test_arf_slim_rows.py` が公開境界ごとに固定している。
_RETAINED_ROW_KEYS = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 18, 20, 21, 22)


def _slim_row(row, shared_names: dict):
    """行を、どの読み手にも同じ答えを返す範囲で細くする（issue #4）。

    生の行は 1 行 6.5 KB あり、大半は誰も読まない入れ子フィールド（Key14 / 16 /
    25 / 36 / 47 ほか）だった。長さは変えない（`len(data) > 25` 等の判定を保つ）。
    """
    if not isinstance(row, list):
        return row
    slim = [None] * len(row)
    for key in _RETAINED_ROW_KEYS:
        if key < len(row):
            slim[key] = row[key]
    name = slim[1] if len(slim) > 1 else None
    if isinstance(name, str):
        # 異なり値は注入の数しか無いので、spot をまたいで同じ str を共有する。
        slim[1] = shared_names.setdefault(name, name)
    rt = _rt_of(row)
    if rt is not None:
        # `_rt_of` は Key15 が無ければ Key16 を見る。縮めた形は Key15 に置く。
        slim[15] = [[1, [rt]]]
    if len(row) > 37 and isinstance(row[37], list):
        slim[37] = row[37][:2]  # [EstimatedNoise, SignalToNoise] しか読まれない
    return slim


def iter_groups(file_like_object) -> Iterator[dict]:
    """`.arf` のスポットを 1 つずつ返す（`deserialize` の 1 要素と同じ形）。

    `{"MasterAlignmentID", "AlignmentID", "RT", "MassCenter", "IonMode", "Name",
    "HeightAverage", "AlignedPeakProperties": [細くした行, ...], "SourceBlockIndex",
    "SourceLocalIndex"}`。ファイルを丸ごと 1 つの bytes に読まず、群は msgpack の
    Unpacker が 1 つずつ取り出す。全体のリストは作らない（呼び出し側が要るものだけ残す）。

    RAM が尽きたら、そこまでに返したスポット数・行数を `ArfMemoryError` に載せて投げ直す。
    """
    spots = rows = 0
    shared_names: dict[str, str] = {}
    try:
        for group_index, (block_index, local_index, group) in enumerate(
            _iter_arf_peak_groups(file_like_object),
            start=0,
        ):
            formatted = extract_arf_data(
                group,
                group_index=group_index,
                source_block_index=block_index,
                source_local_index=local_index,
            )
            if formatted is None:
                continue
            formatted["MasterAlignmentID"] = group_index
            formatted["AlignmentID"] = group_index
            # spot 単位の値（RT / Name 等）は生の代表行から取り終えている。
            formatted["AlignedPeakProperties"] = [
                _slim_row(row, shared_names) for row in group
            ]
            spots += 1
            rows += len(group)
            group = None  # 返したあとまで生の行を抱えない
            yield formatted
    except MemoryError as exc:
        # 規模を運ぶ例外へ詰め替える。ここで確保するのは小さなオブジェクト 1 つで、
        # 失敗した確保（数 GB 単位）とは桁が違う。
        raise ArfMemoryError(spots=spots, rows=rows) from exc


def deserialize(file_like_object) -> list[dict]:
    """`.arf` の全スポットのリスト（`iter_groups` を最後まで回したもの）。"""
    return list(iter_groups(file_like_object))


#: msgpack の ext の型バイト → データ長の読み方（fixext は長さ固定、ext8/16/32 は長さのバイト数）。
_FIXEXT_LENGTH = {0xD4: 1, 0xD5: 2, 0xD6: 4, 0xD7: 8, 0xD8: 16}
_EXT_LENGTH_BYTES = {0xC7: 1, 0xC8: 2, 0xC9: 4}


def _leading_size(payload: bytes) -> int:
    """LZ4 ブロックの先頭の整数（伸長後の大きさ）。整数でなければ 0（そのブロックは読まれない）。"""
    try:
        size = next(msgpack.Unpacker(io.BytesIO(payload), raw=False, strict_map_key=False))
    except (StopIteration, msgpack.OutOfData, ValueError):
        return 0
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        return 0
    return size


def _skip(file_like_object, count: int) -> int:
    """count バイト読み飛ばし、実際に飛ばせたバイト数を返す（末尾が切れていれば少ない）。"""
    if count <= 0:
        return 0
    seekable = getattr(file_like_object, "seekable", None)
    if seekable is not None and seekable():
        here = file_like_object.tell()
        end = file_like_object.seek(0, io.SEEK_END)
        target = min(here + count, end)
        file_like_object.seek(target)
        return target - here
    skipped = 0
    while skipped < count:
        chunk = file_like_object.read(min(count - skipped, 1024 * 1024))
        if not chunk:
            break
        skipped += len(chunk)
    return skipped


def uncompressed_size(file_like_object) -> int:
    """`.arf` の LZ4 ブロックの伸長後の大きさの和（バイト）。伸長しない。

    各ブロックの先頭の整数（MessagePack-CSharp の Lz4Block が書く伸長後の長さ）を足す。
    ext のブロックはヘッダと先頭の整数だけを読み、残りは読み飛ばす。ext でない
    トップレベル（`[ヘッダ, 本体]` の配列）は msgpack に解かせる（`_payload_from_top_level_object`
    と同じ見分け方）。末尾が切れたブロックは、reader も読まないので数えない。
    """
    total = 0
    while True:
        first = file_like_object.read(1)
        if not first:
            return total
        tag = first[0]
        if tag in _FIXEXT_LENGTH:
            length = _FIXEXT_LENGTH[tag]
        elif tag in _EXT_LENGTH_BYTES:
            width = _EXT_LENGTH_BYTES[tag]
            raw_length = file_like_object.read(width)
            if len(raw_length) < width:
                return total
            length = int.from_bytes(raw_length, "big")
        else:
            file_like_object.seek(-1, io.SEEK_CUR)
            for obj in _top_level_stream(file_like_object):
                payload = _payload_from_top_level_object(obj)
                if payload is not None:
                    total += _leading_size(payload)
            return total
        if not file_like_object.read(1):  # ext の型コード
            return total
        head = file_like_object.read(min(length, 9))  # msgpack の整数は最大 9 バイト
        if _skip(file_like_object, length - len(head)) + len(head) < length:
            return total
        total += _leading_size(head)


def _is_alignment_peak_row(data) -> bool:
    """`_convert_to_alignment_feature` が実経路に入る行かどうか。

    file_name だけ欲しい呼び出し（`exclusions._entry_file_name`）と、dict を組む
    本体とで判定が食い違わないよう、ここ 1 か所に置く。
    """
    return (isinstance(data, list)
            and len(data) > 25 and isinstance(data[18], (int, float)))


def _file_name_from_row(data: list) -> str | None:
    """行が実経路のものだと分かっている前提で、先頭付近の文字列を拾う。

    呼び出し側が `_is_alignment_peak_row` を通したあとに使う（ここで再判定すると、
    spot × 注入ぶんの isinstance/len が二重に走る）。
    """
    for item in data[:10]:
        if isinstance(item, (bytes, str)):
            decoded = _decode(item)
            if len(decoded) > 0:
                return decoded
    return None


#: `_convert_to_alignment_feature` がそのまま添字を引くだけのフィールド。
#: 実経路の行は len(data) > 25 が保証されているので、範囲判定は要らない。
_DIRECT_FIELD_INDEX = {
    "file_id": 0,            # 真の FileID は Key0
    "master_peak_id": 2,
    "peak_id": 3,
    "ms2_raw_id": 9,         # Key9
    "height": 18,
    "area": 20,
    "area_above_baseline": 21,  # Key21
    "m_z": 22,
}

#: 行に依らず定数のフィールド（同定側があとから書き換える前提の初期値）。
_CONSTANT_FIELDS = {"is_msms_matched": False, "is_matched": False}


def peak_field(data: list, name: str):
    """`_convert_to_alignment_feature(data).get(name)` と同じ値を dict 無しで返す。

    `.arf` を走査する呼び出し（`mztab.evidence.scan_arf`）は数個の欄しか読まないのに、
    1 行ごとに 16 キーの dict を組むと spot × 注入ぶんの無駄になる（393MB の `.arf` で
    294 万回。`rt` のための `_convert_to_times` まで毎行走る）。

    **同値性は `tests/test_arf_hot_paths.py` が全キー × 行の種類で固定している。**
    フィールドを増やすときは両方に足すこと。
    """
    if not _is_alignment_peak_row(data):
        return None
    index = _DIRECT_FIELD_INDEX.get(name)
    if index is not None:
        return data[index]
    if name == "file_name":
        return _file_name_from_row(data)
    if name == "is_gap_filled":
        # MasterPeakID(Key2)が負(-2)のサンプルはギャップフィル（未検出→補間値）
        return isinstance(data[2], int) and data[2] < 0
    if name == "is_msms":
        # MS2: Key10 MS2RawSpectrumID2CE が非空ならMS/MS取得済み
        return bool(data[10]) if isinstance(data[10], dict) else False
    if name == "rt":
        return _rt_of(data)
    if name in ("estimated_noise", "signal_to_noise"):
        # Key37 PeakShape = [EstimatedNoise, SignalToNoise, ...]（msgpack配列）
        peak_shape = data[37] if len(data) > 37 and isinstance(data[37], list) else []
        position = 0 if name == "estimated_noise" else 1
        return peak_shape[position] if len(peak_shape) > position else None
    return _CONSTANT_FIELDS.get(name)


def _rt_of(data: list):
    rt_value = None
    if len(data) > 15 and isinstance(data[15], list):
        rt_value = _convert_to_times(data[15]).get("rt")
    elif len(data) > 16 and isinstance(data[16], list):
        rt_value = _convert_to_times(data[16]).get("rt")
    return rt_value


def _convert_to_alignment_feature(data: list) -> dict:
    if not isinstance(data, list) or len(data) < 3:
        return {}

    # .arf の AlignmentChromPeakFeature row
    if _is_alignment_peak_row(data):
        rt_value = _rt_of(data)

        file_name = _file_name_from_row(data)

        # Key37 PeakShape = [EstimatedNoise, SignalToNoise, ...]（msgpack配列）
        peak_shape = data[37] if len(data) > 37 and isinstance(data[37], list) else []
        estimated_noise = peak_shape[0] if len(peak_shape) > 0 else None
        signal_to_noise = peak_shape[1] if len(peak_shape) > 1 else None

        # MasterPeakID(Key2)が負(-2)のサンプルはギャップフィル（未検出→補間値）
        master_peak_id = data[2] if len(data) > 2 else None
        is_gap_filled = isinstance(master_peak_id, int) and master_peak_id < 0

        # MS2: Key10 MS2RawSpectrumID2CE が非空ならMS/MS取得済み
        ms2_ce = data[10] if len(data) > 10 else None
        is_msms = bool(ms2_ce) if isinstance(ms2_ce, dict) else False

        return {
            "peak_id": data[3] if len(data) > 3 else None,
            "file_id": data[0] if len(data) > 0 else None,  # 修正: 真のFileIDはKey0
            "file_name": file_name, # 追加
            "master_peak_id": master_peak_id,
            "height": data[18] if len(data) > 18 else None,
            "area": data[20] if len(data) > 20 else None,
            "area_above_baseline": data[21] if len(data) > 21 else None,  # 追加(Key21)
            "m_z": data[22] if len(data) > 22 else None,
            "rt": rt_value,
            "signal_to_noise": signal_to_noise,    # 追加(Key37[1])
            "estimated_noise": estimated_noise,    # 追加(Key37[0])
            "ms2_raw_id": data[9] if len(data) > 9 else None,  # 追加(Key9)
            "is_msms": is_msms,                    # 追加: MS/MS取得有無
            "is_gap_filled": is_gap_filled,        # 追加: 補間値フラグ
            "is_msms_matched": False,
            "is_matched": False,
        }

    # .arf ファイルのみを扱うため、他の形式はサポートしない
    return {}


def alignment_feature_row(data: list) -> dict:
    """`AlignedPeakProperties` の1行（注入1件分）を辞書へ落とす公開境界。

    `.arf` の msgpack Key 番号を知っているのはこのモジュールだけ、という約束を
    層をまたいで保つための入口（`docs/schema/AlignmentChromPeakFeature.md` が
    Key 番号の正準）。`metabolomix/mztab/evidence.py` は private 名を掴まず、ここを呼ぶ。

    読めない行（別形式・短すぎる配列）では空 dict を返す——例外にしないのは、
    `.arf` に非ピーク行が混ざるのが常態だから。呼び出し側は空を「この行からは
    証拠を取れない」として扱う。
    """
    return _convert_to_alignment_feature(data)


def extract_peak_properties(deserialized_list: list[dict]) -> pd.DataFrame:
    rows = []
    for spot in deserialized_list:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            continue
        for sample_index, sample in enumerate(aligned):
            if not isinstance(sample, list) or len(sample) < 3:
                continue
            
            feature = _convert_to_alignment_feature(sample)
            
            # 【追加】抽出したファイル名を優先し、無い場合はSample_0などにフォールバック
            file_name = feature.get("file_name")
            sample_label = file_name if file_name else f"Sample_{sample_index}"
            sample_tag_info = get_sample_peak_tag_info(
                spot, feature.get("file_id"), file_name,
            )
            sample_class_id = get_sample_class_id(
                spot, feature.get("file_id"), file_name,
            )

            row = {
                "MasterAlignmentID": spot.get("MasterAlignmentID"),
                "AlignmentID": spot.get("AlignmentID"),
                "SpotRT": spot.get("RT"),
                "SpotMassCenter": spot.get("MassCenter"),
                "IonMode": spot.get("IonMode"),
                "CompoundName": spot.get("Name"),
                "AlignmentTags": list(spot.get("Tags") or []),
                "SampleIndex": sample_index,
                "FileName": sample_label, # CSVにもファイル名を追加
                "ClassID": sample_class_id,
                "PeakID": feature.get("peak_id"),
                "FileID": feature.get("file_id"),
                "MasterPeakID": feature.get("master_peak_id"),
                "PeakHeight": feature.get("height"),
                "PeakArea": feature.get("area"),
                "PeakAreaAboveBaseline": feature.get("area_above_baseline"),
                "PeakMZ": feature.get("m_z"),
                "PeakRT": feature.get("rt"),
                "SignalToNoise": feature.get("signal_to_noise"),
                "IsMsms": feature.get("is_msms"),
                "IsGapFilled": feature.get("is_gap_filled"),
                "PeakTags": list(sample_tag_info.get("Tags") or []),
            }
            rows.append(row)
    
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def count_peak_property_rows(deserialized_list: list[dict]) -> int:
    """`extract_peak_properties` が返す行数だけを、DataFrame を組まずに数える。

    旧 ARF 要約ツールは「抽出された総ピークレコード数」と「平均サンプル数/スポット」
    しか使っていないのに、そのために spot × 注入ぶんの dict を作って pandas の
    DataFrame へ積んでいた（393MB の `.arf` で 881 万行）。行の採否条件は
    `extract_peak_properties` と同一に保つ——ずれたら報告値が嘘になる。
    """
    total = 0
    for spot in deserialized_list:
        aligned = spot.get("AlignedPeakProperties")
        if not isinstance(aligned, list):
            continue
        for sample in aligned:
            if not isinstance(sample, list) or len(sample) < 3:
                continue
            total += 1
    return total


def plot_peak_height_distribution(peak_df: pd.DataFrame, filter_keyword: str = None, group_regex: str = None, output_file: str = "peak_height_dist.png"):
    """
    指定したキーワード（脂質クラス等）でフィルタリングし、
    横軸にサンプル（またはグループ）、縦軸にPeakHeightをとったジッタープロットを作成します。
    """
    if peak_df.empty:
        print("[WARNING] データフレームが空のためプロットできません。", file=sys.stderr)
        return

    # 1. フィルタリング処理
    if filter_keyword:
        # Warning回避のため .copy() を付与
        df_filtered = peak_df[peak_df['CompoundName'].fillna("").astype(str).str.contains(filter_keyword, case=False)].copy()
        title_text = f"Class/Name: {filter_keyword}"
    else:
        df_filtered = peak_df.copy()
        title_text = "All Annotated Lipids"

    if df_filtered.empty:
        print(f"[WARNING] キーワード '{filter_keyword}' に一致するデータがありません。", file=sys.stderr)
        return

    # 2. レプリケートのグループ化
    if group_regex:
        # 正規表現にマッチした部分（_1 など）を消去してグループ名を作る
        df_filtered['GroupLabel'] = df_filtered['FileName'].apply(lambda x: re.sub(group_regex, '', str(x)))
        title_text += " (Grouped by Replicates)"
    else:
        df_filtered['GroupLabel'] = df_filtered['FileName']

    # 元の順序を保ったままユニークなグループ名を取得
    groups = df_filtered['GroupLabel'].unique()
    
    # 3. 描画設定
    plt.figure(figsize=(max(6, len(groups) * 1.2), 6))
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    
    for i, group in enumerate(groups):
        # グループに属する全レプリケートのデータを抽出
        group_data = df_filtered[df_filtered['GroupLabel'] == group]
        
        y_vals = pd.to_numeric(group_data['PeakHeight'], errors='coerce').dropna().values
        if len(y_vals) == 0:
            continue
            
        # データ点を少し散らす（ジッター）
        x_vals = np.random.normal(i, 0.08, size=len(y_vals))
        
        # プロット（レプリケートがまとまって表示される）
        plt.scatter(x_vals, y_vals, color=colors[i % 10], alpha=0.7, edgecolors='black', linewidths=0.5, s=25)
        
        # 平均値のバーを描画
        mean_val = np.mean(y_vals)
        plt.hlines(mean_val, i - 0.3, i + 0.3, colors='black', linewidth=2.5)

    # 軸とタイトルの設定
    plt.xticks(range(len(groups)), groups, rotation=45, ha="right")
    plt.ylabel("Peak Height")
    plt.title(title_text)
    
    # スタイル調整（上・右の枠線を消す）
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"[INFO] 脂質分布プロットを保存しました: {output_file} (データ件数: {len(df_filtered)}, グループ数: {len(groups)})", file=sys.stderr)


def summarize_arf_data(deserialized_list):
    """
    ARFデータの統計的要約を生成
    """
    if not deserialized_list:
        
        return "データがありません。"

    df = pd.DataFrame(deserialized_list)
    
    summary = {
        "total_peaks": len(df),
        "rt_range": (df['RT'].min(), df['RT'].max()) if 'RT' in df.columns else None,
        "mass_range": (df['MassCenter'].min(), df['MassCenter'].max()) if 'MassCenter' in df.columns else None,
        "height_average_mean": df['HeightAverage'].mean() if 'HeightAverage' in df.columns else None,
        "height_average_max": df['HeightAverage'].max() if 'HeightAverage' in df.columns else None,
        "ion_modes": df['IonMode'].value_counts().to_dict() if 'IonMode' in df.columns else {},
        "named_compounds": df['Name'].notna().sum() if 'Name' in df.columns else 0
    }
    
    return summary

def main():
    parser = argparse.ArgumentParser(description="MS-DIAL .arf ファイルを解析します。")
    parser.add_argument("--file", "-f", required=True, help="解析する .arf ファイルのパス（必須。データフォルダは探さない）")
    parser.add_argument("--export", "-e", help="CSV に書き出すパス")
    parser.add_argument("--group-replicates", action="store_true", help="ファイル名からレプリケート番号（_1, _2など）を自動で除外してグループ化する")
    parser.add_argument("--group-regex", type=str, help="特殊なファイル名の場合に、削除したい文字列の正規表現を直接指定")
    parser.add_argument("--plot-distribution", action="store_true", help="脂質クラスごとのPeakHeight分布プロットを作成する")
    parser.add_argument("--filter-name", type=str, help="分布プロットで抽出する脂質名やクラスのキーワード (例: 'EtherPE_P', 'TG')")
    parser.add_argument("--output-dist-plot", type=str, default=None, help="分布プロットの PNG 保存先。未指定なら PNG を生成しない")
    args = parser.parse_args()

    file_path = args.file
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f".arf ファイルが見つかりません: {file_path}")

    print(f"[INFO] 解析対象ファイル: {file_path}")
    with open(file_path, "rb") as f:
        deserialized = deserialize(io.BytesIO(f.read()))

    print(f"[INFO] 読み込んだピーク数（Spots）: {len(deserialized)}")

    peak_df = extract_peak_properties(deserialized)
    print(f"[INFO] 抽出した Peak Property レコード数: {len(peak_df)}")
    
    if len(peak_df) > 0:
        avg_samples = len(peak_df) / len(deserialized)
        print(f"[INFO] 平均サンプル数/スポット: {avg_samples:.2f} (抽出が成功していればサンプル数と一致します)")

    # ===== デシリアライズ結果（抽出した per-sample ピークプロパティ）をCLI表示 =====
    if len(peak_df) > 0:
        print("\n===== デシリアライズ結果: 抽出ピークプロパティ (先頭20行) =====")
        with pd.option_context("display.max_columns", None, "display.width", 220):
            print(peak_df.head(20).to_string(index=False))
        print(f"... 全 {len(peak_df):,} レコード ({len(deserialized):,} スポット x 平均{avg_samples:.0f}サンプル)")

    if args.export:
        peak_df.to_csv(args.export, index=False)
        print(f"[INFO] CSV に出力しました: {args.export}")


    if args.plot_distribution:
        if len(peak_df) > 0:
            # 正規表現パターンの決定
            regex_pattern = None
            if args.group_regex:
                regex_pattern = args.group_regex
            elif args.group_replicates:
                # デフォルトの優秀な正規表現: 「_数字」の後に「_」か「行末」が続く場合のみ削除
                # 例: "..._0h_1_NEG" -> "..._0h_NEG" に自動変換される。(_0hは削除されない)
                regex_pattern = r'_\d+(?=_|$)'
                
            # PNG は明示要求時のみ。CLI が黙って画像ファイルを作らないようにする。
            if args.output_dist_plot:
                plot_peak_height_distribution(
                    peak_df=peak_df,
                    filter_keyword=args.filter_name,
                    group_regex=regex_pattern,  # 追加
                    output_file=args.output_dist_plot
                )
            else:
                print("[INFO] --output-dist-plot 未指定のため分布プロット PNG は生成しません")
        else:
            print("[WARNING] 抽出されたPeak Propertyレコードが存在しないため、分布プロットをスキップします。", file=sys.stderr)
    

    summary = summarize_arf_data(deserialized)
    pprint.pprint(summary)


if __name__ == "__main__":
    main()
