"""metabolomix.rawmeta: raw の取得時刻（メタデータだけ）・取得の順・群との重なり、ISA / mwTab からのシート。

fixture はテスト自身が tmp に作る（合成の mzML・Thermo のヘッダ・Bruker の tdf・最小の OLE 複合ファイル）。
"""
from __future__ import annotations

import sqlite3
import struct
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from metabolomix.analysis import sample_manifest
from metabolomix.core.atomic_io import DomainError
from metabolomix.rawmeta import acquisition as A
from metabolomix.rawmeta import sample_table as T

T0 = datetime(2026, 3, 1, 9, 0, tzinfo=timezone.utc)


# ---------- 合成の raw ----------

def write_mzml(path: Path, start: datetime | None) -> Path:
    stamp = f' startTimeStamp="{start.isoformat()}"' if start is not None else ""
    path.write_text('<?xml version="1.0" encoding="utf-8"?>\n'
                    '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">'
                    f'<run id="r"{stamp}><spectrumList count="0"/></run></mzML>', encoding="utf-8")
    return path


def _filetime(t: datetime) -> int:
    return int((t - A.EPOCH) / timedelta(microseconds=1)) * 10


def write_thermo_raw(path: Path, start: datetime) -> Path:
    header = bytearray(64)
    struct.pack_into("<H", header, 0, 0xA101)
    header[2:20] = "Finnigan\x00".encode("utf-16-le")
    struct.pack_into("<Q", header, 40, _filetime(start))
    path.write_bytes(bytes(header))
    return path


def write_bruker_d(path: Path, stamp: str, *, tdf: bool = True) -> Path:
    path.mkdir()
    db, table = ("analysis.tdf", "GlobalMetadata") if tdf else ("analysis.sqlite", "Properties")
    con = sqlite3.connect(path / db)
    con.execute(f"CREATE TABLE {table} (Key TEXT, Value TEXT)")
    con.execute(f"INSERT INTO {table} VALUES ('AcquisitionDateTime', ?)", (stamp,))
    con.commit()
    con.close()
    return path


def _dir_entry(name: str, etype: int, ctime: datetime | None, mtime: datetime | None) -> bytes:
    e = bytearray(128)
    raw = (name + "\x00").encode("utf-16-le")
    e[:len(raw)] = raw
    struct.pack_into("<H", e, 64, len(raw))
    e[66] = etype
    struct.pack_into("<III", e, 68, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
    struct.pack_into("<QQ", e, 100, _filetime(ctime) if ctime else 0, _filetime(mtime) if mtime else 0)
    struct.pack_into("<I", e, 116, 0xFFFFFFFE)
    return bytes(e)


def write_wiff(path: Path, start: datetime, end: datetime) -> Path:
    """最小の OLE 複合ファイル（セクタ 512 バイト・FAT 1・ディレクトリ 1）。`SampleSubtree` の作成時刻が開始。"""
    header = bytearray(512)
    header[:8] = A.OLE_MAGIC
    struct.pack_into("<HHH", header, 24, 0x3E, 3, 0xFFFE)
    struct.pack_into("<HH", header, 30, 9, 6)
    struct.pack_into("<I", header, 44, 1)              # FAT のセクタ数
    struct.pack_into("<I", header, 48, 1)              # ディレクトリの最初のセクタ
    struct.pack_into("<I", header, 56, 4096)
    struct.pack_into("<I", header, 60, 0xFFFFFFFE)     # mini FAT 無し
    struct.pack_into("<I", header, 68, 0xFFFFFFFE)     # DIFAT の続き無し
    difat = [0] + [0xFFFFFFFF] * 108
    struct.pack_into("<109I", header, 76, *difat)
    fat = [0xFFFFFFFD, 0xFFFFFFFE] + [0xFFFFFFFF] * 126
    directory = (_dir_entry("Root Entry", 5, None, end) + _dir_entry("SampleSubtree", 1, start, start)
                 + bytes(256))
    path.write_bytes(bytes(header) + struct.pack("<128I", *fat) + directory)
    return path


# ---------- 取得時刻 ----------

def test_read_mzml_start_and_timezone(tmp_path):
    got = A.read_acquisition_times(write_mzml(tmp_path / "a.mzML", T0))
    assert got == {"start_utc": T0.isoformat(), "end_utc": None, "source": "mzml", "error": None}
    local = T0.astimezone(timezone(timedelta(hours=9)))
    assert A.read_acquisition_times(write_mzml(tmp_path / "b.mzML", local))["start_utc"] == T0.isoformat()
    naive = tmp_path / "c.mzML"
    naive.write_text('<mzML><run id="r" startTimeStamp="2026-03-01T09:00:00"/></mzML>', encoding="utf-8")
    assert A.read_acquisition_times(naive)["start_utc"] == T0.isoformat()      # 時刻帯が無ければ UTC
    missing = A.read_acquisition_times(write_mzml(tmp_path / "d.mzML", None))
    assert missing["source"] == "unreadable" and missing["start_utc"] is None and missing["error"]


def test_read_thermo_raw_header(tmp_path):
    got = A.read_acquisition_times(write_thermo_raw(tmp_path / "x.raw", T0))
    assert got["source"] == "thermo_raw" and got["start_utc"] == T0.isoformat()
    bad = tmp_path / "y.raw"
    bad.write_bytes(b"\x00" * 64)
    assert A.read_acquisition_times(bad)["source"] == "unreadable"


def test_read_bruker_tdf_and_baf(tmp_path):
    got = A.read_acquisition_times(write_bruker_d(tmp_path / "s1.d", "2026-03-01T18:00:00+09:00"))
    assert got["source"] == "bruker_tdf" and got["start_utc"] == T0.isoformat()
    baf = A.read_acquisition_times(write_bruker_d(tmp_path / "s2.d", "2026-03-01T09:00:00", tdf=False))
    assert baf["source"] == "bruker_baf" and baf["start_utc"] == T0.isoformat()
    (tmp_path / "s3.d").mkdir()
    assert A.read_acquisition_times(tmp_path / "s3.d")["source"] == "unreadable"


def test_read_wiff_ole_sample_subtree(tmp_path):
    end = T0 + timedelta(minutes=25)
    got = A.read_acquisition_times(write_wiff(tmp_path / "s.wiff", T0, end))
    assert got == {"start_utc": T0.isoformat(), "end_utc": end.isoformat(), "source": "wiff_ole", "error": None}
    not_ole = tmp_path / "t.wiff"
    not_ole.write_bytes(b"not an OLE file" * 40)
    assert A.read_acquisition_times(not_ole)["source"] == "unreadable"


def test_other_formats_are_unreadable_without_raising(tmp_path):
    (tmp_path / "waters.raw").mkdir()                     # Waters の .raw はフォルダ
    (tmp_path / "x.wiff2").write_bytes(b"\x00")
    for name in ("waters.raw", "x.wiff2"):
        got = A.read_acquisition_times(tmp_path / name)
        assert got["source"] == "unreadable" and got["error"]


def test_acquisition_order_ranks_gaps_and_ties(tmp_path):
    files = [write_mzml(tmp_path / "c.mzML", T0),
             write_mzml(tmp_path / "a.mzML", T0 + timedelta(minutes=50)),
             write_mzml(tmp_path / "b.mzML", T0),               # c と同じ時刻 → ファイル名で b が先
             write_mzml(tmp_path / "d.mzML", None)]
    rows = A.acquisition_order(files)
    assert [r["file"] for r in rows] == ["c.mzML", "a.mzML", "b.mzML", "d.mzML"]   # 並びは入力のまま
    assert [r["rank"] for r in rows] == [2, 3, 1, None]
    assert A.gaps_over(rows) == 1 and A.gaps_over(rows, minutes=60) == 0


def _rows(times: dict[str, datetime | None]) -> list[dict]:
    rows = [{"file": f, "start_utc": t.isoformat() if t else None, "rank": None} for f, t in times.items()]
    for rank, r in enumerate(sorted((r for r in rows if r["start_utc"]),
                                    key=lambda r: (r["start_utc"], r["file"])), start=1):
        r["rank"] = rank
    return rows


def test_group_overlap_complete_independent_partial_unknown():
    day2 = T0 + timedelta(days=1)
    complete = _rows({"a1": T0, "a2": T0 + timedelta(minutes=10), "b1": T0 + timedelta(minutes=20),
                      "b2": T0 + timedelta(minutes=30), "qc": T0 + timedelta(minutes=5)})
    groups = {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    out = A.group_overlap(complete, groups)
    assert out["overlap"] == "complete" and out["comparisons"] == [{"groups": ["A", "B"], "overlap": "complete"}]
    assert out["group_rank_ranges"] == {"A": [1, 3], "B": [4, 5]}          # QC も順位に数える
    mixed = _rows({"a1": T0, "b1": T0 + timedelta(minutes=10), "a2": day2, "b2": day2 + timedelta(minutes=5)})
    assert A.group_overlap(mixed, groups)["overlap"] == "independent"
    partial = _rows({"a1": T0, "b1": T0 + timedelta(minutes=10), "a2": T0 + timedelta(minutes=20),
                     "b2": day2})
    assert A.group_overlap(partial, groups)["overlap"] == "partial"
    unknown = _rows({"a1": T0, "a2": None, "b1": T0 + timedelta(minutes=10), "b2": day2})
    out = A.group_overlap(unknown, groups)
    assert out["overlap"] == "unknown" and out["unreadable_grouped"] == ["a2"]
    three = _rows({"a1": T0, "b1": T0 + timedelta(minutes=10), "c1": T0 + timedelta(minutes=20),
                   "a2": day2, "b2": day2 + timedelta(minutes=5), "c2": day2 + timedelta(minutes=50)})
    out = A.group_overlap(three, {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "c1": "C", "c2": "C"})
    assert {tuple(c["groups"]): c["overlap"] for c in out["comparisons"]} == {
        ("A", "B"): "independent", ("A", "C"): "independent", ("B", "C"): "independent"}
    assert A.group_overlap([], {})["overlap"] == "unknown"


# ---------- ISA / mwTab からのシート ----------

def test_manifest_columns_match_sample_manifest_v2():
    assert T.MANIFEST_COLUMNS == sample_manifest.FIELDS_V2


def _isa(tmp_path: Path) -> tuple[Path, Path, list[str]]:
    raw = tmp_path / "raw"
    raw.mkdir()
    times = {"r03.mzML": 0, "r01.mzML": 10, "r02.mzML": 20, "qc1.mzML": 5, "r04.mzML": 30}
    for name, minutes in times.items():
        write_mzml(raw / name, T0 + timedelta(minutes=minutes))
    s_table = tmp_path / "s_study.txt"
    s_table.write_text("Source Name\tSample Name\tFactor Value[Genotype]\tFactor Value[Diet]\n"
                       "m1\tS1\tWT\tchow\nm2\tS2\tKO\tchow\nm3\tS3\tWT\thigh fat\nm4\tS4\t\tchow\n",
                       encoding="utf-8")
    a_table = tmp_path / "a_assay.txt"
    a_table.write_text("Sample Name\tRaw Spectral Data File\n"
                       "S1\traw/r01.mzML\nS2\tr02.mzML\nS3\tr03.mzML\nS4\tr04.mzML\n", encoding="utf-8")
    return raw, s_table, a_table


def test_build_sheet_from_isa_groups_factors_and_acquisition_order(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    names = sorted(p.name for p in raw.iterdir())
    built = T.build_sample_sheet(raw, names, s_table, ["Genotype", "Diet"], assay_table=a_table,
                                 qc_selector=["Sample Name", "^qc"])
    assert built["injection_order_source"] == "raw_acquisition_time"
    assert built["columns"] == [*T.MANIFEST_COLUMNS, "Genotype", "Diet"]
    by_file = {r["source_file"]: r for r in built["rows"]}
    assert by_file["r01.mzML"]["group"] == "WT_chow" and by_file["r01.mzML"]["Genotype"] == "WT"
    assert by_file["r03.mzML"]["group"] == "WT_high_fat"
    assert by_file["r04.mzML"]["group"] == T.UNASSIGNED          # 因子の値が空
    assert by_file["qc1.mzML"]["role"] == "qc" and by_file["qc1.mzML"]["group"] == "QC"
    # 注入順は取得時刻の順（r03 → qc1 → r01 → r02 → r04）
    assert [r["source_file"] for r in built["rows"]] == ["r03.mzML", "qc1.mzML", "r01.mzML", "r02.mzML",
                                                         "r04.mzML"]
    assert [r["injection_order"] for r in built["rows"]] == ["1", "2", "3", "4", "5"]
    assert all(r["batch"] == "" for r in built["rows"])
    assert built["files_not_in_assay"] == ["qc1.mzML"] and built["unassigned_samples"] == ["S4"]
    assert built["factor_columns"] == ["Factor Value[Genotype]", "Factor Value[Diet]"]
    assert built["text"].startswith("# schema = sample-manifest.v2\n")


def test_build_sheet_leaves_order_empty_when_a_raw_is_unreadable(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    (raw / "r02.mzML").write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")   # 開始時刻が無い
    names = sorted(p.name for p in raw.iterdir())
    built = T.build_sample_sheet(raw, names, s_table, ["Genotype"], assay_table=a_table)
    assert built["injection_order_source"] == "none" and built["unreadable"] == ["r02.mzML"]
    assert all(r["injection_order"] == "" for r in built["rows"])
    assert [r["source_file"] for r in built["rows"]] == sorted(names)


def test_build_sheet_from_mwtab_raw_file_names(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, name in enumerate(["b.mzML", "a.mzML", "c.mzML"]):
        write_mzml(raw / name, T0 + timedelta(minutes=10 * k))
    mwtab = tmp_path / "ST000001.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000001\n"
                     "SUBJECT_SAMPLE_FACTORS\tP1\tS_a\tTreatment:ctrl | Time:0h\tRAW_FILE_NAME=a.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\tP2\tS_b\tTreatment:drug | Time:0h\tRAW_FILE_NAME=C:\\raw\\b.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\tP3\tS_c\tTreatment:drug | Time:6h\tRAW_FILE_NAME=c\n",
                     encoding="utf-8")
    built = T.build_sample_sheet(raw, ["a.mzML", "b.mzML", "c.mzML"], mwtab, ["Treatment", "Time"])
    assert built["sample_table_format"] == "mwtab"
    groups = {r["source_file"]: (r["biological_sample_id"], r["group"]) for r in built["rows"]}
    assert groups == {"a.mzML": ("S_a", "ctrl_0h"), "b.mzML": ("S_b", "drug_0h"), "c.mzML": ("S_c", "drug_6h")}
    assert [r["source_file"] for r in built["rows"]] == ["b.mzML", "a.mzML", "c.mzML"]


def test_build_sheet_rejects_bad_factor_columns_and_formats(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    names = sorted(p.name for p in raw.iterdir())
    for bad in ([], ["group"], ["Genotype", "genotype"], ["Missing"], ["Geno\ttype"]):
        with pytest.raises(DomainError) as e:
            T.build_sample_sheet(raw, names, s_table, bad, assay_table=a_table)
        assert e.value.code == "SAMPLE_TABLE_INVALID"
    with pytest.raises(DomainError) as e:
        T.build_sample_sheet(raw, names, s_table, ["Genotype"], sample_table_format="csv")
    assert e.value.code == "SAMPLE_TABLE_INVALID"
    with pytest.raises(DomainError):
        T.build_sample_sheet(raw, names, s_table, ["Genotype"], qc_selector=["Sample Name", "("])


def test_factor_group_names():
    assert T.factor_group(["WT", "high fat"]) == "WT_high_fat"
    assert T.factor_group(["WT", ""]) == T.UNASSIGNED and T.factor_group([]) == T.UNASSIGNED
