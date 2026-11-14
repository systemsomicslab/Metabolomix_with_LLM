"""raw の取得時刻（メタデータだけ）と、取得の順・間隔・群との重なり（stdlib だけ）。

スペクトルは解さない。読むのは取得の開始時刻の 1 欄だけ:

- `.wiff`（SCIEX）は OLE 複合ファイル。Analyst が取得時に書く `SampleSubtree` の作成時刻が開始、ルートの更新時刻が終了
  （ファイルシステムの更新時刻は取得ではなく複製・ダウンロードの時刻なので使わない）。
- Thermo の `.raw`（ファイル）は先頭のファイルヘッダ（0xA101・UTF-16 の `Finnigan`）のオフセット 40 の FILETIME
  （RawFileReader の `FileHeader.CreationDate` と一致することを reanalysis-study で確かめた値）。
- Bruker の `.d`（フォルダ）は `analysis.tdf`（表 GlobalMetadata）か `analysis.sqlite`（表 Properties）の鍵
  `AcquisitionDateTime`（読み取り専用・immutable で開く）。
- mzML は `<run startTimeStamp>`。

ほかの形式（Waters の `.raw` フォルダ・`.wiff2`・`.abf` など）は止まらずに `unreadable` とする。
reanalysis-study の `reanalysis/acquisition.py` から移した（判定の数字は道具が出す）。
"""
from __future__ import annotations

import sqlite3
import struct
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timedelta, timezone
from itertools import combinations
from pathlib import Path

__all__ = ["GAP_MINUTES", "read_acquisition_times", "acquisition_order", "gaps_over",
           "group_overlap"]

EPOCH = datetime(1601, 1, 1, tzinfo=timezone.utc)
OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
NO_SECTOR = (0xFFFFFFFE, 0xFFFFFFFF)
#: これを超える取得の間隔を数える（分）。長い中断はドリフトとバッチの切れ目の目安。
GAP_MINUTES = 40
_UNREADABLE = "unreadable"


def _filetime(v: int) -> datetime | None:
    if v == 0:
        return None
    try:
        t = EPOCH + timedelta(microseconds=v // 10)
    except OverflowError:
        return None
    # 取得の時刻として妥当な範囲だけを採る（壊れた欄を時刻として読まない）
    return t if 1990 <= t.year <= 2100 else None


def _ole_entries(path: Path) -> list[tuple[str, int, datetime | None, datetime | None]]:
    """OLE 複合ファイルのディレクトリエントリ（名前・種類・作成時刻・更新時刻）。"""
    with path.open("rb") as fh:
        hdr = fh.read(512)
        if len(hdr) < 512 or hdr[:8] != OLE_MAGIC:
            raise ValueError("OLE 複合ファイルではない")
        sector_shift, mini_shift = struct.unpack_from("<HH", hdr, 30)
        if sector_shift not in (9, 12) or mini_shift != 6:  # 壊れたヘッダで 0 割りや巨大な読み込みをしない
            raise ValueError("OLE のヘッダが壊れている（セクタサイズが不正）")
        ssz = 1 << sector_shift
        n_fat = struct.unpack_from("<I", hdr, 44)[0]
        first_dir = struct.unpack_from("<I", hdr, 48)[0]
        difat = list(struct.unpack_from("<109I", hdr, 76))
        nxt = struct.unpack_from("<I", hdr, 68)[0]
        guard = 0
        while nxt not in NO_SECTOR and len(difat) < n_fat + 109 and guard < 100000:
            fh.seek((nxt + 1) * ssz)
            blk = struct.unpack(f"<{ssz // 4}I", fh.read(ssz))
            difat.extend(blk[:-1])
            nxt = blk[-1]
            guard += 1
        fat_sectors = [s for s in difat if s not in NO_SECTOR][:n_fat]

        def fat_entry(i: int) -> int:
            per = ssz // 4
            fh.seek((fat_sectors[i // per] + 1) * ssz + (i % per) * 4)
            return struct.unpack("<I", fh.read(4))[0]

        entries = []
        s, guard = first_dir, 0
        while s not in NO_SECTOR and guard < 100000:
            fh.seek((s + 1) * ssz)
            blk = fh.read(ssz)
            for off in range(0, len(blk) - 127, 128):
                e = blk[off:off + 128]
                nlen = struct.unpack_from("<H", e, 64)[0]
                etype = e[66]
                if etype == 0:
                    continue
                name = e[: max(nlen - 2, 0)].decode("utf-16-le", errors="replace")
                ct, mt = struct.unpack_from("<QQ", e, 100)
                entries.append((name, etype, _filetime(ct), _filetime(mt)))
            s = fat_entry(s)
            guard += 1
        return entries


def _read_wiff(path: Path) -> tuple[datetime | None, datetime | None]:
    entries = _ole_entries(path)
    subtree = [e for e in entries if e[0] == "SampleSubtree" and e[1] == 1]
    root = [e for e in entries if e[1] == 5]
    return (subtree[0][2] if subtree else None), (root[0][3] if root else None)


def _read_mzml(path: Path) -> tuple[datetime | None, datetime | None]:
    for _, el in ET.iterparse(path, events=("start",)):
        if el.tag.rsplit("}", 1)[-1] == "run":
            stamp = el.get("startTimeStamp")
            if not stamp:
                return None, None
            t = datetime.fromisoformat(stamp)
            if t.tzinfo is None:  # 時刻帯が無ければ UTC とみなす
                t = t.replace(tzinfo=timezone.utc)
            return t.astimezone(timezone.utc), None
    return None, None


def _read_thermo_raw(path: Path) -> tuple[datetime | None, datetime | None]:
    with path.open("rb") as fh:
        hdr = fh.read(48)
    if (len(hdr) < 48 or struct.unpack_from("<H", hdr, 0)[0] != 0xA101
            or hdr[2:20].decode("utf-16-le", errors="replace").rstrip("\x00") != "Finnigan"):
        raise ValueError("Thermo .raw のヘッダではない")
    return _filetime(struct.unpack_from("<Q", hdr, 40)[0]), None


def _read_bruker_d(path: Path) -> tuple[str, datetime]:
    """.d フォルダの AcquisitionDateTime。(source, 開始時刻)。.tdf を先に見る。"""
    for db, table, source in (("analysis.tdf", "GlobalMetadata", "bruker_tdf"),
                              ("analysis.sqlite", "Properties", "bruker_baf")):
        f = path / db
        if not f.is_file():
            continue
        con = sqlite3.connect(f.as_uri() + "?mode=ro&immutable=1", uri=True)
        try:
            row = con.execute(f"SELECT Value FROM {table} WHERE Key = 'AcquisitionDateTime'").fetchone()
        finally:
            con.close()
        if not row or not row[0]:
            raise ValueError(f"{db} に AcquisitionDateTime が無い")
        t = datetime.fromisoformat(str(row[0]))
        if t.tzinfo is None:  # 時刻帯が無ければ UTC とみなす
            t = t.replace(tzinfo=timezone.utc)
        return source, t.astimezone(timezone.utc)
    raise ValueError(".d に analysis.tdf も analysis.sqlite も無い")


def read_acquisition_times(path: Path) -> dict:
    """1 つの raw の取得時刻。`{"start_utc", "end_utc", "source", "error"}`。

    `source` は `wiff_ole` / `mzml` / `thermo_raw` / `bruker_tdf` / `bruker_baf` / `unreadable`。
    読めないときは例外にせず `unreadable` と理由（`error`）を返す。時刻は UTC の ISO 8601。
    """
    path = Path(path)
    unreadable = {"start_utc": None, "end_utc": None, "source": _UNREADABLE}
    suffix = path.suffix.lower()
    if suffix == ".d" and path.is_dir():
        try:
            source, start = _read_bruker_d(path)
        except (OSError, ValueError, sqlite3.Error) as exc:
            return {**unreadable, "error": str(exc) or type(exc).__name__}
        return {"start_utc": start.isoformat(), "end_utc": None, "source": source, "error": None}
    if suffix == ".raw" and path.is_file():
        try:
            start, _ = _read_thermo_raw(path)
        except (OSError, ValueError, struct.error) as exc:
            return {**unreadable, "error": str(exc) or type(exc).__name__}
        if start is None:
            return {**unreadable, "error": "開始時刻が無い"}
        return {"start_utc": start.isoformat(), "end_utc": None, "source": "thermo_raw", "error": None}
    if suffix not in (".wiff", ".mzml") or not path.is_file():
        what = "フォルダ" if path.is_dir() else "形式"
        return {**unreadable, "error": f"取得時刻を読めない{what}: {path.suffix or '(拡張子なし)'}"}
    try:
        start, end = (_read_wiff if suffix == ".wiff" else _read_mzml)(path)
    except (OSError, ValueError, struct.error, ET.ParseError, IndexError) as exc:
        return {**unreadable, "error": str(exc) or type(exc).__name__}
    source = "wiff_ole" if suffix == ".wiff" else "mzml"
    if start is None:
        return {**unreadable, "error": "開始時刻が無い"}
    return {"start_utc": start.isoformat(), "end_utc": end.isoformat() if end else None,
            "source": source, "error": None}


def acquisition_order(paths) -> list[dict]:
    """raw ごとの取得時刻と、そこからの順（1 始まり。読めないものは `rank` が None）。

    並びは入力の順のまま。順は読めた時刻の昇順で、同じ時刻はファイル名で決める（ファイル名の順は
    副キーにだけ使う）。行は `{"file", "path", "start_utc", "end_utc", "source", "error", "rank"}`。
    """
    rows = []
    for p in paths:
        p = Path(p)
        rows.append({"file": p.name, "path": str(p), **read_acquisition_times(p), "rank": None})
    readable = sorted((r for r in rows if r["start_utc"]), key=lambda r: (r["start_utc"], r["file"]))
    for rank, row in enumerate(readable, start=1):
        row["rank"] = rank
    return rows


def gaps_over(rows: list[dict], minutes: float = GAP_MINUTES) -> int:
    """読めた時刻を順に並べたとき、隣どうしの間隔が `minutes` を超える箇所の数。"""
    starts = sorted(datetime.fromisoformat(r["start_utc"]) for r in rows if r.get("start_utc"))
    return sum(1 for a, b in zip(starts, starts[1:]) if (b - a) > timedelta(minutes=minutes))


def _judge_pair(ranks: dict[str, list[int]], dates: dict[str, Counter], a: str, b: str) -> str:
    """complete: 2 群の順位の範囲が重ならない。independent: 重なり、かつ取得日（UTC）ごとに両群が居る。
    それ以外は partial。"""
    ra, rb = ranks[a], ranks[b]
    if max(ra) < min(rb) or max(rb) < min(ra):
        return "complete"
    days = [d for d, c in dates.items() if c[a] or c[b]]
    return "independent" if all(dates[d][a] and dates[d][b] for d in days) else "partial"


def group_overlap(rows: list[dict], group_of: dict[str, str]) -> dict:
    """群（`group_of`: ファイル名 → 群）と取得の順の重なりを、群の組ごとに判定する。

    組ごとの `overlap`: `complete`（2 群の取得順の範囲が重ならない＝取得順と群が完全に交絡）/
    `independent`（重なり、かつどちらかの群がいる取得日のすべてに両群がいる）/ `partial`（それ以外）/
    `unknown`（群のある試料に取得時刻の読めないものがある）。群の無いファイル（QC など）は判定に
    使わない（順位は全ファイルの中の順のまま）。全体の `overlap` は組の判定が 1 種類ならそれ、
    混ざれば `partial`、組が無ければ `unknown`。
    """
    grouped = [r for r in rows if group_of.get(r["file"])]
    unreadable = sorted(r["file"] for r in grouped if r["rank"] is None)
    ranks: dict[str, list[int]] = {}
    dates: dict[str, Counter] = {}
    for r in grouped:
        if r["rank"] is None:
            continue
        g = group_of[r["file"]]
        ranks.setdefault(g, []).append(r["rank"])
        dates.setdefault(r["start_utc"][:10], Counter())[g] += 1
    groups = sorted({group_of[r["file"]] for r in grouped})
    pairs = list(combinations(groups, 2))
    if unreadable:
        comparisons = [{"groups": [a, b], "overlap": "unknown"} for a, b in pairs]
    else:
        comparisons = [{"groups": [a, b], "overlap": _judge_pair(ranks, dates, a, b)}
                       for a, b in pairs]
    states = {c["overlap"] for c in comparisons}
    overall = (next(iter(states)) if len(states) == 1 else "partial") if states else "unknown"
    return {"overlap": overall, "comparisons": comparisons,
            "group_rank_ranges": {g: [min(v), max(v)] for g, v in sorted(ranks.items())},
            "by_date": {d: dict(sorted(c.items())) for d, c in sorted(dates.items())},
            "unreadable_grouped": unreadable}
