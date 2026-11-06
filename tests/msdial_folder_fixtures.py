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


def write_msdial_folder(
    folder,
    *,
    spots: list[dict],
    samples: list[dict],
    with_mddata: bool = True,
    stem: str = "AlignmentResult_2026_02_01_00_00_00",
) -> dict:
    """`.arf` / `.arf2`（と `.mddata`）を 1 組、フォルダへ書く。

    spots[i]: name, mz, rt, ontology, adduct, heights, gap_filled, msms, matches, ion_mode（すべて任意）。
    samples[j]: file_id, name, class_id（必須）と file_type / order / batch、
    `.arf` に書く名前 arf_name（既定は name。改名後の `.mddata` との突き合わせ用）。
    """
    from tests.curation_fixtures import arf2_spot_raw, arf_row, write_arf, write_arf2

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    groups, raws = [], []
    for spot_id, spot in enumerate(spots):
        mz = spot.get("mz", 500.0)
        rt = spot.get("rt", 12.0)
        heights = spot.get("heights") or [100.0] * len(samples)
        gap = spot.get("gap_filled") or [False] * len(samples)
        msms = spot.get("msms") or [False] * len(samples)
        rows = []
        for sample, height, gap_filled, has_msms in zip(samples, heights, gap, msms):
            row = arf_row(file_id=sample["file_id"], mz=mz, rt=rt, height=height,
                          gap_filled=gap_filled)
            row[1] = sample.get("arf_name", sample["name"])
            if has_msms:
                row[10] = {1: 30.0}
            rows.append(row)
        groups.append(rows)
        raws.append(arf2_spot_raw(
            spot_id=spot_id, name=spot.get("name", "Unknown"), mz=mz, rt=rt,
            ontology=spot.get("ontology", ""), adduct=spot.get("adduct", "[M+H]+"),
            ion_mode=spot.get("ion_mode", 0), matches=spot.get("matches")))
    paths = {
        "arf": write_arf(folder / f"{stem}_PeakProperties.arf", groups),
        "arf2": write_arf2(folder / f"{stem}.arf2", raws),
    }
    if with_mddata:
        paths["mddata"] = write_mddata(folder / f"{stem}.mddata", [
            analysis_file_row(s["file_id"], s["name"], s["class_id"],
                              file_type=s.get("file_type", 0), order=s.get("order"),
                              batch=s.get("batch", 1))
            for s in samples])
    return paths
