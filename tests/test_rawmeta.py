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


#: MassHunter が書く `AcqData/Contents.xml` の形（BOM・名前空間の宣言・7 桁の小数秒・時刻帯つき）。
AGILENT_CONTENTS = (
    '﻿<?xml version="1.0" encoding="UTF-8"?>\n'
    '<Contents xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
    '  <Version>4</Version>\n'
    '  <AcquiredTime>{stamp}</AcquiredTime>\n'
    '  <AcqStatus>2</AcqStatus>\n'
    '  <InstrumentName>LC-QTOF</InstrumentName>\n'
    '</Contents>\n')


def write_agilent_d(path: Path, stamp: str) -> Path:
    (path / "AcqData").mkdir(parents=True)
    (path / "AcqData" / "Contents.xml").write_text(AGILENT_CONTENTS.format(stamp=stamp), encoding="utf-8")
    (path / "AcqData" / "MSScan.bin").write_bytes(b"\x00" * 16)
    return path


def test_read_agilent_d_acquired_time(tmp_path):
    got = A.read_acquisition_times(write_agilent_d(tmp_path / "a.d", "2024-02-07T17:18:28.2110274+09:00"))
    start = datetime(2024, 2, 7, 8, 18, 28, 211027, tzinfo=timezone.utc)   # 17:18 JST = 08:18 UTC（7 桁目は切る）
    assert got == {"start_utc": start.isoformat(), "end_utc": None, "source": "agilent_acqdata", "error": None}
    # 独立に計算した値と突き合わせる（fromisoformat は 6 桁の小数秒まで読める）
    assert got["start_utc"] == datetime.fromisoformat("2024-02-07T17:18:28.211027+09:00").astimezone(
        timezone.utc).isoformat()
    # Z・時が 1 桁の形と、時刻帯が「:」なしの形
    z = A.read_acquisition_times(write_agilent_d(tmp_path / "b.d", "2005-11-20T5:53:39Z"))
    assert z["start_utc"] == datetime(2005, 11, 20, 5, 53, 39, tzinfo=timezone.utc).isoformat()
    west = A.read_acquisition_times(write_agilent_d(tmp_path / "c.d", "2024-02-07T17:18:28-0530"))
    assert west["start_utc"] == datetime(2024, 2, 7, 22, 48, 28, tzinfo=timezone.utc).isoformat()
    naive = A.read_acquisition_times(write_agilent_d(tmp_path / "e.d", "2024-02-07T08:18:28"))
    assert naive["start_utc"] == datetime(2024, 2, 7, 8, 18, 28, tzinfo=timezone.utc).isoformat()  # UTC とみなす


def test_bruker_wins_over_agilent_and_missing_layouts_name_both_vendors(tmp_path):
    both = write_bruker_d(tmp_path / "s.d", "2026-03-01T18:00:00+09:00")
    write_agilent_d(both, "2020-01-01T00:00:00Z")
    assert A.read_acquisition_times(both)["source"] == "bruker_tdf"
    (tmp_path / "empty.d").mkdir()
    got = A.read_acquisition_times(tmp_path / "empty.d")
    assert got["source"] == "unreadable"
    assert "Bruker" in got["error"] and "Agilent" in got["error"] and "Contents.xml" in got["error"]
    broken = write_agilent_d(tmp_path / "broken.d", "not a time")
    assert A.read_acquisition_times(broken)["source"] == "unreadable"
    bad_xml = tmp_path / "bad.d"
    (bad_xml / "AcqData").mkdir(parents=True)
    (bad_xml / "AcqData" / "Contents.xml").write_text("<Contents><AcquiredTime>", encoding="utf-8")
    assert A.read_acquisition_times(bad_xml)["source"] == "unreadable"


def test_wiff2_is_unreadable_with_a_format_specific_reason(tmp_path):
    """SCIEX OS の .wiff2 は stdlib で読める取得時刻の欄が無い。理由は形式に即して書き、変換の道を示す。"""
    (tmp_path / "s1.wiff2").write_bytes(b"\x00" * 64)
    got = A.read_acquisition_times(tmp_path / "s1.wiff2")
    assert got["source"] == "unreadable" and got["start_utc"] is None
    assert ".wiff2" in got["error"] and "mzML" in got["error"]


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
                   "a2": day2, "b2": day2 + timedelta(minutes=5), "c2": day2 + timedelta(minutes=30)})
    out = A.group_overlap(three, {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "c1": "C", "c2": "C"})
    assert {tuple(c["groups"]): c["overlap"] for c in out["comparisons"]} == {
        ("A", "B"): "independent", ("A", "C"): "independent", ("B", "C"): "independent"}
    assert A.group_overlap([], {})["overlap"] == "unknown"


JST = timezone(timedelta(hours=9))


def test_sessions_split_at_gaps_not_at_utc_midnight():
    """UTC+9 の朝 8:50–9:50 に交互に測った 1 本の取得は、UTC の日付の境目（9:00 JST）で割れず independent。"""
    morning = {"b1": datetime(2026, 11, 10, 8, 50, tzinfo=JST), "a1": datetime(2026, 11, 10, 9, 10, tzinfo=JST),
               "b2": datetime(2026, 11, 10, 9, 30, tzinfo=JST), "a2": datetime(2026, 11, 10, 9, 50, tzinfo=JST)}
    rows = _rows({f: t.astimezone(timezone.utc) for f, t in morning.items()})
    groups = {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    assert len({r["start_utc"][:10] for r in rows}) == 2          # UTC の日付では 2 日に割れる取得
    assert A.gaps_over(rows) == 0 and set(A.acquisition_sessions(rows).values()) == {1}
    out = A.group_overlap(rows, groups)
    assert out["overlap"] == "independent"
    assert [s["groups"] for s in out["by_session"]] == [{"A": 2, "B": 2}]


def test_sessions_split_at_long_gaps_within_one_utc_day():
    """逆に、UTC では同じ日でも 40 分を超える中断で区切れたまとまりに片方の群しかいなければ partial。"""
    local = {"a1": datetime(2026, 11, 10, 9, 10, tzinfo=JST), "b1": datetime(2026, 11, 10, 9, 30, tzinfo=JST),
             "a2": datetime(2026, 11, 10, 14, 0, tzinfo=JST), "a3": datetime(2026, 11, 10, 14, 20, tzinfo=JST),
             "b2": datetime(2026, 11, 11, 8, 0, tzinfo=JST)}
    rows = _rows({f: t.astimezone(timezone.utc) for f, t in local.items()})
    assert len({r["start_utc"][:10] for r in rows}) == 1          # UTC ではすべて 11-10
    groups = {"a1": "A", "a2": "A", "a3": "A", "b1": "B", "b2": "B"}
    sessions = A.acquisition_sessions(rows)
    assert sessions == {"a1": 1, "b1": 1, "a2": 2, "a3": 2, "b2": 3}
    out = A.group_overlap(rows, groups)
    assert out["overlap"] == "partial"
    assert [s["groups"] for s in out["by_session"]] == [{"A": 1, "B": 1}, {"A": 2}, {"B": 1}]


def test_qc_injections_keep_a_session_continuous():
    """群の無い QC も間隔の判定に使う（群の注入だけを並べると 40 分を超えて見える取得を割らない）。"""
    times = {"a1": T0, "qc1": T0 + timedelta(minutes=30), "b1": T0 + timedelta(minutes=60),
             "qc2": T0 + timedelta(minutes=90), "a2": T0 + timedelta(minutes=120), "b2": T0 + timedelta(minutes=125)}
    rows = _rows(times)
    out = A.group_overlap(rows, {"a1": "A", "a2": "A", "b1": "B", "b2": "B"})
    assert out["overlap"] == "independent" and len(out["by_session"]) == 1


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
    # 群の決まらない生体試料は統計の母集団から外す（`unassigned` を 1 つの群にしない）
    assert by_file["r04.mzML"]["include"] == "false"
    assert {r["include"] for f, r in by_file.items() if f != "r04.mzML"} == {"true"}
    assert by_file["qc1.mzML"]["role"] == "qc" and by_file["qc1.mzML"]["group"] == "QC"
    # 注入順は取得時刻の順（r03 → qc1 → r01 → r02 → r04）
    assert [r["source_file"] for r in built["rows"]] == ["r03.mzML", "qc1.mzML", "r01.mzML", "r02.mzML",
                                                         "r04.mzML"]
    assert [r["injection_order"] for r in built["rows"]] == ["1", "2", "3", "4", "5"]
    assert all(r["batch"] == "" for r in built["rows"])
    assert built["files_not_in_assay"] == ["qc1.mzML"] and built["unassigned_samples"] == ["S4"]
    assert built["factor_columns"] == ["Factor Value[Genotype]", "Factor Value[Diet]"]
    # schema の行の直後に、注入順を取得時刻から作った宣言（読み手が出所 raw_acquisition_time にする）
    assert built["text"].startswith("# schema = sample-manifest.v2\n" + T.ORDER_SOURCE_LINE + "\n")


def test_one_unreadable_raw_keeps_the_order_of_the_rest(tmp_path):
    """読めない raw が 1 つあっても、読めた raw の取得時刻の順は消さない（読めない raw の注入順だけ空）。"""
    raw, s_table, a_table = _isa(tmp_path)
    (raw / "r02.mzML").write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")   # 開始時刻が無い
    names = sorted(p.name for p in raw.iterdir())
    built = T.build_sample_sheet(raw, names, s_table, ["Genotype"], assay_table=a_table)
    assert built["injection_order_source"] == T.ORDER_FROM_ACQUISITION_PARTIAL
    assert built["unreadable"] == ["r02.mzML"] and built["unreadable_reasons"] == {"開始時刻が無い": 1}
    # 読めた 4 つは取得時刻（r03 0 分・qc1 5 分・r01 10 分・r04 30 分）の順位、読めない r02 は空で後ろ
    order = {r["source_file"]: r["injection_order"] for r in built["rows"]}
    assert order == {"r03.mzML": "1", "qc1.mzML": "2", "r01.mzML": "3", "r04.mzML": "4", "r02.mzML": ""}
    assert [r["source_file"] for r in built["rows"]][-1] == "r02.mzML"
    assert built["text"].startswith("# schema = sample-manifest.v2\n" + T.ORDER_SOURCE_LINE + "\n")
    # 読み手は空の r02 に「取得時刻では分からない」印を付け、読めた行は出所 raw_acquisition_time のまま
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(built["text"], encoding="utf-8")
    rows = {r["source_file"].replace("\\", "/").rsplit("/", 1)[-1]: r
            for r in sample_manifest.parse_manifest(sheet, source_root=raw, expected_sources=names)}
    assert rows["r02.mzML"]["injection_order"] is None and sample_manifest.sheet_order_unknown(rows["r02.mzML"])
    assert rows["r03.mzML"]["injection_order"] == 1 and not sample_manifest.sheet_order_unknown(rows["r03.mzML"])
    assert rows["r03.mzML"]["provenance"]["injection_order"]["derived_from"] == "raw_acquisition_time"


def test_no_readable_raw_leaves_every_order_empty_without_the_declaration(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    for name in ("r01.mzML", "r02.mzML", "r03.mzML", "r04.mzML", "qc1.mzML"):
        (raw / name).write_text("<mzML><run id='r'/></mzML>", encoding="utf-8")
    names = sorted(p.name for p in raw.iterdir())
    built = T.build_sample_sheet(raw, names, s_table, ["Genotype"], assay_table=a_table)
    assert built["injection_order_source"] == "none" and built["unreadable"] == names
    assert all(r["injection_order"] == "" for r in built["rows"])
    assert [r["source_file"] for r in built["rows"]] == sorted(names)
    assert T.ORDER_SOURCE_LINE not in built["text"]              # 注入順が無ければ宣言も書かない


def test_build_sheet_reads_agilent_d_folders(tmp_path):
    """Agilent の .d（AcqData/Contents.xml）から注入順を作る（Bruker 専用の読み手で全体の順を失わない）。"""
    raw = tmp_path / "raw"
    raw.mkdir()
    jst = timezone(timedelta(hours=9))
    stamps = {"b.d": datetime(2026, 3, 1, 18, 0, tzinfo=jst), "a.d": datetime(2026, 3, 1, 18, 20, tzinfo=jst),
              "c.d": datetime(2026, 3, 1, 17, 40, tzinfo=jst)}
    for name, t in stamps.items():
        write_agilent_d(raw / name, t.isoformat())
    mwtab = tmp_path / "ST000002.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000002\n"
                     + "".join(f"SUBJECT_SAMPLE_FACTORS\tP{k}\tS_{n[0]}\tTreatment:{t}\tRAW_FILE_NAME={n}\n"
                               for k, (n, t) in enumerate([("a.d", "ctrl"), ("b.d", "drug"), ("c.d", "ctrl")])),
                     encoding="utf-8")
    built = T.build_sample_sheet(raw, ["a.d", "b.d", "c.d"], mwtab, ["Treatment"])
    assert built["injection_order_source"] == "raw_acquisition_time" and built["unreadable"] == []
    assert {r["source_file"]: r["injection_order"] for r in built["rows"]} == {"c.d": "1", "b.d": "2", "a.d": "3"}
    assert {r["source"] for r in built["acquisition"]} == {"agilent_acqdata"}


def test_order_source_line_is_a_declaration_the_manifest_reader_knows():
    """rawmeta は sample_manifest を import しないので、宣言の行の形と値を突き合わせる。"""
    key, value = T.ORDER_SOURCE_LINE.removeprefix("# ").split(" = ")
    assert key == sample_manifest.ORDER_SOURCE_DECLARATION
    assert value in sample_manifest.ORDER_SOURCE_DECLARED_VALUES


def test_tables_that_are_not_utf8_are_sample_table_invalid(tmp_path):
    """UTF-8 でない試料の表・assay の表・mwTab は素の UnicodeDecodeError でなく SAMPLE_TABLE_INVALID。"""
    raw, s_table, a_table = _isa(tmp_path)
    names = sorted(p.name for p in raw.iterdir())
    latin = tmp_path / "s_latin.txt"
    latin.write_bytes("Sample Name\tFactor Value[Dose]\nS1\t5 \u00b5g\n".encode("latin-1"))
    for table, assay, factors in ((latin, None, ["Dose"]), (s_table, latin, ["Genotype"])):
        with pytest.raises(DomainError) as e:
            T.build_sample_sheet(raw, names, table, factors, assay_table=assay)
        assert e.value.code == "SAMPLE_TABLE_INVALID" and "UTF-8" in e.value.message
        assert e.value.details["path"] == str(latin)
    mwtab = tmp_path / "ST_latin.txt"
    mwtab.write_bytes("#METABOLOMICS WORKBENCH\nSUBJECT_SAMPLE_FACTORS\tP1\tS1\tDose:5 \u00b5g\t\n"
                      .encode("latin-1"))
    for fmt in ("auto", "mwtab"):
        with pytest.raises(DomainError) as e:
            T.build_sample_sheet(raw, names, mwtab, ["Dose"], sample_table_format=fmt)
        assert e.value.code == "SAMPLE_TABLE_INVALID"


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


# ---------- 群の決まらない試料・列名の衝突・ファイルと試料の対応 ----------

def test_unassigned_samples_are_left_out_of_the_statistics(tmp_path):
    """群 unassigned の生体試料は include=false で、統計エンジンの検出率の足切り・群の数に入らない。"""
    from types import SimpleNamespace

    import numpy as np

    from metabolomix.analysis.stats import data as stats_data
    from metabolomix.analysis.stats import preprocess as P

    raw = tmp_path / "raw"
    raw.mkdir()
    names = ["c1.mzML", "c2.mzML", "d1.mzML", "d2.mzML", "x1.mzML", "x2.mzML"]
    for k, name in enumerate(names):
        write_mzml(raw / name, T0 + timedelta(minutes=5 * k))
    s_table = tmp_path / "s.txt"
    s_table.write_text("Sample Name\tFactor Value[Treatment]\n"
                       "c1\tCtrl\nc2\tCtrl\nd1\tCase\nd2\tCase\nx1\t\n", encoding="utf-8")   # x1 は値が空、x2 は行が無い
    built = T.build_sample_sheet(raw, names, s_table, ["Treatment"])
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(built["text"], encoding="utf-8")
    rows = sample_manifest.parse_manifest(sheet, source_root=raw, expected_sources=names)
    by_file = {Path(r["source_file"]).name: r for r in rows}
    assert {f: (r["group"], r["include"]) for f, r in by_file.items() if f.startswith("x")} == {
        "x1.mzML": (T.UNASSIGNED, False), "x2.mzML": (T.UNASSIGNED, False)}
    assert built["unassigned_samples"] == ["x1", "x2"]
    # 特徴 f は unassigned の 2 試料にだけ検出される。群として数えれば検出率の足切り（群ごと 50%）を通ってしまう
    detected = np.array([[True, False], [True, False], [True, False], [True, False], [False, True], [False, True]])
    ordered = [by_file[n] for n in names]
    ds = SimpleNamespace(feature_matrix=np.where(detected, 100.0, 1.0).T.tolist(), detected_mask=detected.T.tolist(),
                         feature_ids=["g", "f"], sample_names=[r["sample_id"] for r in ordered],
                         sample_metadata_rows=ordered)
    prep = P.preprocess(stats_data.from_dataset(ds), {"filters": {"blank_ratio": None, "qc_rsd": None,
                                                                  "d_ratio": None},
                                                      "normalize": {"method": "none"}, "drift": {"mode": "off"}})
    select = next(s for s in prep.steps if s["step"] == "select")
    assert sorted(select["excluded_by_include"]) == ["x1", "x2"]
    assert {u.group for u in prep.units} == {"Ctrl", "Case"}
    assert prep.feature_ids == ["g"]                              # f は Ctrl・Case のどちらにも検出されない


def test_factor_names_that_only_differ_in_case_from_sheet_columns_are_kept(tmp_path):
    """mwTab の因子 `Group`・`Batch` は列 group・batch と大文字小文字だけ違う: そのままの名前で書く。"""
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, name in enumerate(["a.mzML", "b.mzML", "c.mzML", "d.mzML"]):
        write_mzml(raw / name, T0 + timedelta(minutes=10 * k))
    mwtab = tmp_path / "ST000003.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000003\n"
                     "SUBJECT_SAMPLE_FACTORS\tP1\tS_a\tGroup:ctrl | Time:0h | Batch:1\tRAW_FILE_NAME=a.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\tP2\tS_b\tGroup:drug | Time:0h | Batch:1\tRAW_FILE_NAME=b.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\tP3\tS_c\tGroup:ctrl | Time:6h | Batch:2\tRAW_FILE_NAME=c.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\tP4\tS_d\tGroup:drug | Time:6h | Batch:2\tRAW_FILE_NAME=d.mzML\n",
                     encoding="utf-8")
    names = ["a.mzML", "b.mzML", "c.mzML", "d.mzML"]
    built = T.build_sample_sheet(raw, names, mwtab, ["Group", "Time", "Batch"])
    assert built["columns"][-3:] == ["Group", "Time", "Batch"] and built["renamed_factor_columns"] == {}
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(built["text"], encoding="utf-8")
    rows = sample_manifest.parse_manifest(sheet, source_root=raw, expected_sources=names)
    d = next(r for r in rows if r["source_file"].endswith("d.mzML"))
    # 列 group は値の組、因子の列 Group はその因子だけ（2 因子の計画が引ける）。列 batch は空のまま
    assert (d["group"], d["Group"], d["Time"], d["batch"], d["Batch"]) == ("drug_6h_2", "drug", "6h", None, "2")


def test_factor_names_equal_to_sheet_columns_are_renamed_not_refused(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, name in enumerate(["a.mzML", "b.mzML"]):
        write_mzml(raw / name, T0 + timedelta(minutes=10 * k))
    s_table = tmp_path / "s.txt"
    s_table.write_text("Sample Name\tgroup\tbatch\tfactor_batch\n"
                       "a\tctrl\t1\tx\nb\tdrug\t2\ty\n", encoding="utf-8")
    built = T.build_sample_sheet(raw, ["a.mzML", "b.mzML"], s_table, ["group", "batch", "factor_batch"])
    assert built["renamed_factor_columns"] == {"group": "factor_group", "batch": "factor_batch_2"}
    assert built["columns"][-3:] == ["factor_group", "factor_batch_2", "factor_batch"]
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(built["text"], encoding="utf-8")
    rows = sample_manifest.parse_manifest(sheet, source_root=raw, expected_sources=["a.mzML", "b.mzML"])
    b = next(r for r in rows if r["source_file"].endswith("b.mzML"))
    assert (b["group"], b["factor_group"], b["batch"], b["factor_batch_2"], b["factor_batch"]) == (
        "drug_2_y", "drug", None, "2", "y")


def _assay(tmp_path: Path, header: str, lines: list[str]) -> Path:
    path = tmp_path / "a_assay.txt"
    path.write_text(header + "\n" + "".join(line + "\n" for line in lines), encoding="utf-8")
    return path


def _two_mice(tmp_path: Path, files=("QE_001.mzML", "QE_002.mzML")) -> tuple[Path, Path]:
    raw = tmp_path / "raw"
    raw.mkdir()
    for k, name in enumerate(files):
        write_mzml(raw / name, T0 + timedelta(minutes=10 * k))
    s_table = tmp_path / "s.txt"
    s_table.write_text("Sample Name\tFactor Value[Genotype]\nMouse_A\tWT\nMouse_B\tKO\n", encoding="utf-8")
    return raw, s_table


def test_isa_assay_maps_derived_spectral_data_files(tmp_path):
    """MetaboLights は変換した mzML を Derived Spectral Data File に載せる。"""
    raw, s_table = _two_mice(tmp_path)
    a_table = _assay(tmp_path, "Sample Name\tRaw Spectral Data File\tDerived Spectral Data File",
                     ["Mouse_A\tQE_001.raw\tFILES/QE_001.mzML", "Mouse_B\tQE_002.raw\tFILES/QE_002.mzML"])
    built = T.build_sample_sheet(raw, ["QE_001.mzML", "QE_002.mzML"], s_table, ["Genotype"], assay_table=a_table)
    assert {r["source_file"]: r["group"] for r in built["rows"]} == {"QE_001.mzML": "WT", "QE_002.mzML": "KO"}
    assert built["files_not_in_assay"] == [] and built["unassigned_samples"] == []


def test_isa_assay_falls_back_to_the_stem_like_mwtab(tmp_path):
    """assay の表が .raw の名前だけで、渡すのが変換した .mzML でも、幹が 1 つの試料に当たれば結ぶ。"""
    raw, s_table = _two_mice(tmp_path)
    a_table = _assay(tmp_path, "Sample Name\tRaw Spectral Data File", ["Mouse_A\tQE_001.raw", "Mouse_B\tQE_002.raw"])
    built = T.build_sample_sheet(raw, ["QE_001.mzML", "QE_002.mzML"], s_table, ["Genotype"], assay_table=a_table)
    assert {r["source_file"]: r["biological_sample_id"] for r in built["rows"]} == {
        "QE_001.mzML": "Mouse_A", "QE_002.mzML": "Mouse_B"}


def test_a_sheet_where_no_file_matches_is_not_written(tmp_path):
    raw, s_table = _two_mice(tmp_path)
    a_table = _assay(tmp_path, "Sample Name\tRaw Spectral Data File", ["Mouse_A\tX_001.raw", "Mouse_B\tX_002.raw"])
    with pytest.raises(DomainError) as e:
        T.build_sample_sheet(raw, ["QE_001.mzML", "QE_002.mzML"], s_table, ["Genotype"], assay_table=a_table)
    assert e.value.code == "SAMPLE_TABLE_INVALID" and "どれも試料の表の行に当たらない" in e.value.message
    assert e.value.details["files_without_table_row"] == ["QE_001.mzML", "QE_002.mzML"]
    # 当たっても因子の値がすべて空なら、同じく書かない
    empty = tmp_path / "s_empty.txt"
    empty.write_text("Sample Name\tFactor Value[Genotype]\nQE_001\t\nQE_002\t\n", encoding="utf-8")
    with pytest.raises(DomainError) as e:
        T.build_sample_sheet(raw, ["QE_001.mzML", "QE_002.mzML"], empty, ["Genotype"])
    assert e.value.code == "SAMPLE_TABLE_INVALID" and "値が空" in e.value.message


def test_conflicting_duplicates_are_unassigned_in_isa_and_mwtab_alike(tmp_path):
    """同じ試料が因子の違う 2 行にある・1 つのファイルが別々の試料に当たる: ISA も mwTab も、どちらの群にも
    入れず unassigned・include=false にして返す（黙って最初の行を採らない）。値の同じ重複は通す。"""
    raw = tmp_path / "raw"
    raw.mkdir()
    files = ["P07.mzML", "P08.mzML", "P09.mzML", "X.mzML"]
    for k, name in enumerate(files):
        write_mzml(raw / name, T0 + timedelta(minutes=10 * k))
    s_table = tmp_path / "s.txt"
    s_table.write_text("Source Name\tSample Name\tFactor Value[Treatment]\n"
                       "src1\tP07\tCtrl\nsrc2\tP07\tDrug\nsrc3\tP08\tCtrl\nsrc4\tP08\tCtrl\nsrc5\tP09\tDrug\n"
                       "src6\tA\tCtrl\nsrc7\tB\tDrug\n", encoding="utf-8")
    a_table = _assay(tmp_path, "Sample Name\tRaw Spectral Data File",
                     ["P07\tP07.mzML", "P08\tP08.mzML", "P09\tP09.mzML", "A\tX.mzML", "B\tX.mzML"])
    isa = T.build_sample_sheet(raw, files, s_table, ["Treatment"], assay_table=a_table)
    mwtab = tmp_path / "ST000004.txt"
    mwtab.write_text("#METABOLOMICS WORKBENCH STUDY_ID:ST000004\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tP07\tTreatment:Ctrl\tRAW_FILE_NAME=P07.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tP07\tTreatment:Drug\tRAW_FILE_NAME=P07.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tP08\tTreatment:Ctrl\tRAW_FILE_NAME=P08.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tP08\tTreatment:Ctrl\tRAW_FILE_NAME=P08.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tP09\tTreatment:Drug\tRAW_FILE_NAME=P09.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tA\tTreatment:Ctrl\tRAW_FILE_NAME=X.mzML\n"
                     "SUBJECT_SAMPLE_FACTORS\t-\tB\tTreatment:Drug\tRAW_FILE_NAME=X.mzML\n", encoding="utf-8")
    mw = T.build_sample_sheet(raw, files, mwtab, ["Treatment"])
    for built in (isa, mw):
        got = {r["source_file"]: (r["group"], r["include"]) for r in built["rows"]}
        assert got == {"P07.mzML": (T.UNASSIGNED, "false"), "P08.mzML": ("Ctrl", "true"),
                       "P09.mzML": ("Drug", "true"), "X.mzML": (T.UNASSIGNED, "false")}
        assert built["ambiguous_files"] == [{"file": "X.mzML", "samples": ["A", "B"]}]
        assert [(c["sample"], c["values"]) for c in built["conflicting_samples"]] == [
            ("P07", [["Ctrl"], ["Drug"]])]


def test_build_sheet_rejects_bad_factor_columns_and_formats(tmp_path):
    raw, s_table, a_table = _isa(tmp_path)
    names = sorted(p.name for p in raw.iterdir())
    for bad in ([], ["Genotype", "genotype"], ["Missing"], ["Geno\ttype"]):
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
