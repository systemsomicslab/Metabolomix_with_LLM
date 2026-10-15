"""`.arf` を読むときに、どの消費者も読まない行フィールドを保持しないことを固定する（issue #4）。

実測（185 MB / 35,803 spot × 41 注入 = 147 万行）では、生の行を丸ごと抱えると 1 行
6.5 KB・ピーク 12.2 GB になり、その大半は誰も読まない入れ子フィールド（Key14 / 16 /
25 / 36 / 47 など）だった。3.5 万 spot × 226 注入では RAM 31 GB を超えて読み終わらない。

固定するのは 2 つ:

1. `deserialize` が返す行は、公開境界（`alignment_feature_row` / `peak_field` /
   `file_name_of`）と、行を添字で読む msdial 側の識別ヘルパにとって、生の行と
   **同じ答え**を返すこと。答えが変わらない範囲でだけ捨てている。
2. 誰も読まない重いフィールドが実際に捨てられていること、そしてブロック全体を
   1 個の list として展開しないこと（行を細くしても、展開の時点で全行が同時に
   Python オブジェクトになればピークは下がらない）。
"""
import io

import lz4.block
import msgpack
import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.msdial import classes as msdial_classes
from metabolomix.msdial import tags as msdial_tags

ROW_LENGTH = 53
#: 実データで重かった、どの消費者も読まない添字。
UNREAD_HEAVY_INDICES = (14, 16, 24, 25, 36, 47)


def _raw_row(file_id: int, name: str, height: float, *, rt_in_key16: bool = False) -> list:
    row = [None] * ROW_LENGTH
    row[0] = file_id
    row[1] = name
    row[2] = -2 if file_id % 2 else file_id + 100     # 奇数は gap-fill
    row[3] = file_id + 1000
    row[4:9] = [7, 8, 9, 10, 11]
    row[9] = 42
    row[10] = {35: 3} if file_id % 2 == 0 else {}
    times = [[1, [3.25 + file_id]], [2, [-1.0]], [3, [701.5]], [4, [-1.0]]]
    if rt_in_key16:
        row[15] = None
        row[16] = times
    else:
        row[15] = times
        row[16] = [[1, [9.9]], [3, [701.5]]]
    row[18] = height
    row[20] = height * 2
    row[21] = height * 1.5
    row[22] = 701.5
    row[23] = 0
    row[24] = "PC 34:1"
    row[37] = [12.0, 34.0, 1.0, 2.0, [5.0, 6.0]]
    for index in (14, 25, 36, 47):
        row[index] = [[float(i), float(i) * 2.0] for i in range(40)]
    return row


def _raw_group(seed: int, *, rt_in_key16: bool = False) -> list:
    return [_raw_row(i, f"sample_{i}.wiff", 1000.0 + seed + i, rt_in_key16=rt_in_key16)
            for i in range(4)]


def _ext_block(inner_objects: list) -> bytes:
    stream = b"".join(msgpack.packb(o, use_bin_type=True) for o in inner_objects)
    payload = msgpack.packb(len(stream), use_bin_type=True) + lz4.block.compress(
        stream, store_size=False,
    )
    return msgpack.packb(msgpack.ExtType(99, payload), use_bin_type=True)


def _load(groups: list) -> list[dict]:
    return arf_reader.deserialize(io.BytesIO(_ext_block([[0, 0, *groups]])))


_PEAK_FIELD_NAMES = (
    *arf_reader._DIRECT_FIELD_INDEX,
    "file_name", "is_gap_filled", "is_msms", "rt",
    "estimated_noise", "signal_to_noise", "is_msms_matched", "is_matched",
)


@pytest.mark.parametrize("rt_in_key16", [False, True])
def test_kept_rows_answer_every_reader_like_the_raw_rows(rt_in_key16):
    raw_groups = [_raw_group(1, rt_in_key16=rt_in_key16), _raw_group(2, rt_in_key16=rt_in_key16)]

    spots = _load(raw_groups)

    assert len(spots) == 2
    for spot, raw_group in zip(spots, raw_groups):
        kept_rows = spot["AlignedPeakProperties"]
        assert len(kept_rows) == len(raw_group)
        for kept, raw in zip(kept_rows, raw_group):
            assert len(kept) == len(raw)
            assert arf_reader.alignment_feature_row(kept) == arf_reader.alignment_feature_row(raw)
            for name in _PEAK_FIELD_NAMES:
                assert arf_reader.peak_field(kept, name) == arf_reader.peak_field(raw, name), name
            assert arf_reader.file_name_of(kept) == arf_reader.file_name_of(raw)
            assert msdial_tags._sample_row_identity(kept) == msdial_tags._sample_row_identity(raw)
            assert msdial_classes._arf_sample_identity(kept) == msdial_classes._arf_sample_identity(raw)


def test_spot_level_fields_still_come_from_the_raw_representative_row():
    spots = _load([_raw_group(1)])

    spot = spots[0]
    assert spot["RT"] == 3.25
    assert spot["MassCenter"] == 701.5
    assert spot["IonMode"] == "Positive"
    assert spot["Name"] == "PC 34:1"
    assert spot["HeightAverage"] == 1001.0


def test_unread_heavy_fields_are_not_retained():
    spots = _load([_raw_group(1)])

    for row in spots[0]["AlignedPeakProperties"]:
        for index in UNREAD_HEAVY_INDICES:
            assert row[index] is None, index
        # PeakShape(Key37) は先頭 2 要素（EstimatedNoise, SignalToNoise）しか読まれない。
        assert row[37] == [12.0, 34.0]


def test_sample_name_strings_are_shared_across_spots():
    """file_name は注入の数しか異なり値が無い。spot ごとに別の str を抱えない。"""
    spots = _load([_raw_group(1), _raw_group(2)])

    first = spots[0]["AlignedPeakProperties"][0][1]
    second = spots[1]["AlignedPeakProperties"][0][1]
    assert first == second
    assert first is second


def test_deserialize_does_not_materialize_a_whole_block_as_one_list(monkeypatch):
    """展開済みブロック全体を `list(...)` にする経路を通らないこと。"""
    def refuse(_payload):
        raise AssertionError("deserialize must stream groups, not list a whole block")

    monkeypatch.setattr(arf_reader, "_decode_lz4_msgpack_payload", refuse)

    spots = _load([_raw_group(1), _raw_group(2)])

    assert [spot["MasterAlignmentID"] for spot in spots] == [0, 1]
