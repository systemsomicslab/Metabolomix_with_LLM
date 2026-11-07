""".arf を 1 回のストリーミングで読み、解析フォルダの経路でセルの dict を作らない（B-2a Task 4）。

fixture はすべてこのテストが tmp に作る。
"""
from __future__ import annotations

import inspect
import math
from pathlib import Path

import lz4.block
import msgpack
import numpy as np
import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.mztab import evidence
from tests.curation_fixtures import arf_row, write_arf
from tests.msdial_folder_fixtures import write_msdial_folder


def _groups(n_spots: int = 5, n_samples: int = 4, base: float = 300.0) -> list[list[list]]:
    groups = []
    for k in range(n_spots):
        groups.append([
            arf_row(file_id=i, mz=base + k + 0.001 * i, rt=1.0 + k,
                    height=100.0 * (i + 1) + k, gap_filled=(i == 1 and k % 2 == 0))
            for i in range(n_samples)])
    return groups


def _block(groups: list[list[list]]) -> bytes:
    inner = msgpack.packb([0, 0, *groups], use_bin_type=True)
    payload = msgpack.packb(len(inner), use_bin_type=True) + lz4.block.compress(
        inner, store_size=False)
    return msgpack.packb(msgpack.ExtType(99, payload), use_bin_type=True)


def _write_two_blocks(path: Path) -> Path:
    """トップレベルのブロックが 2 つある `.arf`（スポットの通し番号がブロックをまたぐ）。"""
    path.write_bytes(_block(_groups(3)) + _block(_groups(4, base=600.0)))
    return path


def _reference(data: bytes) -> list[dict]:
    """変更前の `deserialize` の組み立て（全体を bytes で持って群を変換する）。"""
    out = []
    shared: dict[str, str] = {}
    for index, (block, local, group) in enumerate(arf_reader._iter_arf_peak_groups(data)):
        formatted = arf_reader.extract_arf_data(
            group, group_index=index, source_block_index=block, source_local_index=local)
        if formatted is not None:
            formatted["AlignedPeakProperties"] = [arf_reader._slim_row(r, shared) for r in group]
            out.append(formatted)
    return out


@pytest.mark.parametrize("make", [
    lambda p: write_arf(p, _groups()),
    _write_two_blocks,
])
def test_iter_groups_matches_deserialize(tmp_path, make):
    path = make(tmp_path / "a.arf")
    with open(path, "rb") as handle:
        streamed = list(arf_reader.iter_groups(handle))
    with open(path, "rb") as handle:
        whole = arf_reader.deserialize(handle)
    assert streamed == whole
    assert streamed == _reference(path.read_bytes())
    assert [g["MasterAlignmentID"] for g in streamed] == list(range(len(streamed)))


def test_iter_groups_is_lazy(tmp_path, monkeypatch):
    path = write_arf(tmp_path / "a.arf", _groups(n_spots=6))
    calls: list[int] = []
    real = arf_reader.extract_arf_data

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(arf_reader, "extract_arf_data", counting)
    with open(path, "rb") as handle:
        groups = arf_reader.iter_groups(handle)
        assert inspect.isgenerator(groups)
        first = next(groups)
        assert first["MasterAlignmentID"] == 0
        # 1 つ目を取った時点では、最後の群（6 つ目）をまだ変換していない。
        assert len(calls) == 1
        assert len(list(groups)) == 5
    assert len(calls) == 6


def test_uncompressed_size_sums_block_headers_without_decompressing(tmp_path, monkeypatch):
    path = _write_two_blocks(tmp_path / "a.arf")
    expected = sum(len(msgpack.packb([0, 0, *g], use_bin_type=True))
                   for g in (_groups(3), _groups(4, base=600.0)))

    def no_decompress(*args, **kwargs):
        raise AssertionError("伸長した")

    monkeypatch.setattr(lz4.block, "decompress", no_decompress)
    with open(path, "rb") as handle:
        assert arf_reader.uncompressed_size(handle) == expected


def _spec() -> dict:
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(4)]
    spots = [
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0],
         "gap_filled": [False, True, False, False], "msms": [True, False, False, False]},
        {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]},
        {"name": "Unknown", "mz": 501.1, "heights": [5.0, 6.0, 7.0, 8.0],
         "gap_filled": [True, True, True, True]},
    ]
    return {"spots": spots, "samples": samples}


def test_folder_load_builds_no_cell_dicts(tmp_path, monkeypatch):
    write_msdial_folder(tmp_path, **_spec())

    def boom(*args, **kwargs):
        raise AssertionError("セルの dict を作った")

    monkeypatch.setattr(arf_reader, "alignment_feature_row", boom)
    monkeypatch.setattr(arf_reader, "_convert_to_alignment_feature", boom)
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_matrix.shape == (3, 4)
    assert ds.detected_mask.sum() == 7


# ---------- scan_arf と、変更前のセルの dict の組み立ての一致 ----------

def _number(value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return math.nan
    return float(value)


def _old_expectation(path: Path) -> dict:
    """変更前の `msdial_adapter._cells` と `group_representative_mz` の組み立て。"""
    with open(path, "rb") as handle:
        groups = arf_reader.deserialize(handle)
    cells_by_spot, spot_mz = [], []
    for group in groups:
        rows = group.get("AlignedPeakProperties") or []
        cells: dict[int, dict] = {}
        features = []
        for raw in rows:
            feature = arf_reader.alignment_feature_row(raw)
            if feature:
                features.append(feature)
            if feature and isinstance(feature.get("file_id"), int):
                cells[feature["file_id"]] = feature
        cells_by_spot.append(cells)
        spot_mz.append(evidence.spot_representative_mz(features, group.get("MassCenter")))
    file_ids = sorted({fid for cells in cells_by_spot for fid in cells})
    names: dict[int, str | None] = {}
    for cells in cells_by_spot:
        for fid, cell in cells.items():
            names.setdefault(fid, cell.get("file_name"))
    shape = (len(groups), len(file_ids))
    height, area = np.full(shape, np.nan), np.full(shape, np.nan)
    mz = np.full(shape, np.nan)
    detected = np.zeros(shape, dtype=bool)
    for i, cells in enumerate(cells_by_spot):
        for j, fid in enumerate(file_ids):
            cell = cells.get(fid)
            if cell is None:
                continue
            height[i, j] = _number(cell.get("height"))
            area[i, j] = _number(cell.get("area"))
            if not cell.get("is_gap_filled"):
                detected[i, j] = True
                mz[i, j] = _number(cell.get("m_z"))
    msms = np.array([any(c.get("is_msms") for c in cells.values()) for cells in cells_by_spot],
                    dtype=bool)
    return {"file_ids": file_ids, "file_names": names, "height": height, "area": area,
            "mz": mz, "detected": detected, "msms": msms, "spot_mz": spot_mz}


def _odd_groups() -> list[list[list]]:
    """行の順が FileID 順でない・同じ FileID が 2 行・ピークでない行・名前の無い行・数でない値。"""
    def row(fid, mz, height, gap=False):
        return arf_row(file_id=fid, mz=mz, rt=1.0, height=height, gap_filled=gap)

    shuffled = [row(2, 400.02, 30.0), row(0, 400.00, 10.0), row(1, 400.01, 50.0, gap=True),
                row(3, 400.03, 20.0)]
    duplicated = [row(0, 410.0, 10.0), row(1, 410.1, 20.0), row(1, 410.2, 5.0),
                  row(3, 410.3, 40.0)]
    not_a_peak = [0] * 30
    not_a_peak[18] = "not a number"  # Key18 が数でない行はピークの行でない（読まない）
    nameless = row(5, 420.5, 60.0)  # FileID 5 はこのスポットにしか無い（他のスポットでは欠ける列）
    nameless[1] = None
    odd_values = row(6, 420.0, 70.0)  # 自分だけの FileID（非ピークの行や他の行に上書きされない）
    odd_values[20] = None
    odd_values[22] = True
    with_extras = [odd_values, not_a_peak, row(1, 420.1, 15.0), nameless]
    all_gap = [row(i, 430.0 + i, 5.0, gap=True) for i in range(4)]
    with_msms = [row(i, 440.0, 1.0 + i) for i in range(4)]
    with_msms[2][10] = {1: 30.0}
    return [shuffled, duplicated, with_extras, all_gap, with_msms]


def _shifted() -> dict:
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(3)]
    spots = [{"mz": 300.0 + 10 * k, "heights": [100.0, 110.0 + k, 90.0],
              "gap_filled": [k % 3 == 0, False, False]} for k in range(8)]
    return {"spots": spots, "samples": samples}


@pytest.mark.parametrize("make", [
    lambda d: write_msdial_folder(d, **_spec())["arf"],
    lambda d: write_msdial_folder(d, **_shifted())["arf"],
    lambda d: write_arf(d / "odd.arf", _odd_groups()),
    lambda d: _write_two_blocks(d / "two.arf"),
], ids=["spec", "shifted", "odd", "two_blocks"])
def test_scan_equals_old_cells(tmp_path, make):
    path = make(tmp_path)
    old = _old_expectation(path)
    scan = evidence.scan_arf(path, measures=("height", "area"))
    assert scan.path == Path(path)
    assert scan.file_ids == old["file_ids"]
    assert scan.file_names == old["file_names"]
    assert scan.spot_mz == old["spot_mz"]
    np.testing.assert_array_equal(scan.values["height"], old["height"])
    np.testing.assert_array_equal(scan.values["area"], old["area"])
    np.testing.assert_array_equal(scan.detected, old["detected"])
    np.testing.assert_array_equal(scan.per_sample_mz, old["mz"])
    np.testing.assert_array_equal(scan.msms_assigned, old["msms"])
