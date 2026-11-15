"""raw の取得時刻（メタデータだけ）と、取得の順・間隔・群との重なり（stdlib だけ）。

スペクトルは解さない。読むのは取得の開始時刻の 1 欄だけ:

- `.wiff`（SCIEX）は OLE 複合ファイル。Analyst が取得時に書く `SampleSubtree` の作成時刻が開始、ルートの更新時刻が終了
  （ファイルシステムの更新時刻は取得ではなく複製・ダウンロードの時刻なので使わない）。
- Thermo の `.raw`（ファイル）は先頭のファイルヘッダ（0xA101・UTF-16 の `Finnigan`）のオフセット 40 の FILETIME
  （RawFileReader の `FileHeader.CreationDate` と一致することを reanalysis-study で確かめた値）。
- `.d`（フォルダ）は Bruker と Agilent を見分ける。Bruker は `analysis.tdf`（表 GlobalMetadata）か
  `analysis.sqlite`（表 Properties）の鍵 `AcquisitionDateTime`（読み取り専用・immutable で開く）。どちらも無ければ
  Agilent（MassHunter）の `AcqData/Contents.xml` の `AcquiredTime`（時刻帯つきの ISO 8601。小数秒は 7 桁まで）。
  Agilent の終了時刻は読まない（`end_utc` は None）。
- mzML は `<run startTimeStamp>`。
- SCIEX OS の `.wiff2` は中身の形式が公開されておらず、stdlib で読める取得時刻の欄が無い。SCIEX OS は同じ取得に
  同じ幹の `.wiff`（OLE 複合ファイル。`.wiff.scan`・`.timeseries.data` と並ぶ）も書くので、その `.wiff` の時刻を読む
  （`source` は `wiff_ole_sibling`）。同じ幹の `.wiff` が無ければ `unreadable`（mzML に変換すれば `startTimeStamp` を読む）。

ほかの形式（Waters の `.raw` フォルダ・`.abf` など）は止まらずに `unreadable` とし、
形式ごとに何が無かったかを `error` に書く。
reanalysis-study の `reanalysis/acquisition.py` から移した（判定の数字は道具が出す）。

群と取得順の重なり（`group_overlap`）は、取得の時刻を間隔で区切った取得のまとまり（`acquisition_sessions`）
ごとに見る。暦の日付（UTC でも現地でも）で区切ると、夜をまたぐ 1 本の連続した取得や、UTC+9 の朝 9 時をまたぐ
取得が 2 つに割れる。区切る間隔（`session_gap_minutes`）は 40 分（`GAP_MINUTES`）と、隣どうしの取得の間隔の
中央値（LC の 1 周期の目安）の `SESSION_GAP_FACTOR` 倍の大きいほう——1 周期が 40 分を超える長いグラジエントの
メソッドで、連続した取得の 1 本 1 本をまとまりに割らない。開始が `SAME_ACQUISITION_SECONDS` 秒より近い取得
（ベンダーの raw と変換した mzML、`.wiff2` と同じ幹の `.wiff` のように同じ取得を指すファイル）は中央値では 1 つに数える。
"""
from __future__ import annotations

import re
import sqlite3
import statistics
import struct
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timedelta, timezone
from itertools import combinations
from pathlib import Path

__all__ = ["GAP_MINUTES", "SAME_ACQUISITION_SECONDS", "SESSION_GAP_FACTOR", "read_acquisition_times",
           "acquisition_order", "gaps_over", "session_gap_minutes", "acquisition_sessions", "group_overlap"]

EPOCH = datetime(1601, 1, 1, tzinfo=timezone.utc)
OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
NO_SECTOR = (0xFFFFFFFE, 0xFFFFFFFF)
#: これを超える取得の間隔を数える（分）。長い中断はドリフトとバッチの切れ目の目安。取得のまとまりを
#: 区切る間隔の下限でもある（`session_gap_minutes`）。
GAP_MINUTES = 40
#: 取得のまとまりを区切る間隔は、隣どうしの取得の間隔の中央値（LC の 1 周期の目安）のこの倍より長いもの。
#: 1 本・2 本の注入が抜けた（失敗・別のフォルダの QC）だけの間隔では区切らない。
SESSION_GAP_FACTOR = 3
#: 開始がこれ（秒）より近い 2 つの取得は同じ取得を指すファイルとみなし、1 周期の中央値の間隔に数えない。
#: ベンダーの raw（Thermo の FILETIME・Agilent の 7 桁の小数秒）と、それを変換した mzML の `startTimeStamp`
#: （整数秒）は 1 秒足らず違うので、完全に同じ時刻だけを 1 つにすると、N 組の双子が 2N-1 個の間隔の N 個を
#: 1 秒未満にして中央値を潰し、区切りを `GAP_MINUTES` に戻す。LC-MS の 1 周期は 1 分を下回らない。
#: `GAP_MINUTES / SESSION_GAP_FACTOR`（13.3 分）以下なら、捨てた間隔はもともと下限より短い区切りしか作らない。
SAME_ACQUISITION_SECONDS = 60
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


#: Agilent の `AcquiredTime` の形。`Contents.xml` の注記どおり時が 1 桁のことがあり、小数秒は 7 桁まで出る
#: （`2024-02-07T17:18:28.2110274+09:00`）。
_AGILENT_TIME = re.compile(
    r"(\d{4})-(\d{2})-(\d{2})T(\d{1,2}):(\d{2}):(\d{2})(?:\.(\d+))?(Z|[+-]\d{2}:?\d{2})?")


def _agilent_time(text: str) -> datetime:
    """`AcquiredTime` の値を UTC の時刻にする（時刻帯が無ければ UTC とみなす。小数秒は 6 桁で切る）。"""
    m = _AGILENT_TIME.fullmatch(text.strip())
    if not m:
        raise ValueError(f"Agilent の AcquiredTime の形が読めない: {text.strip()[:40]!r}")
    year, month, day, hour, minute, second, frac, zone = m.groups()
    tz = timezone.utc
    if zone and zone != "Z":
        digits = zone[1:].replace(":", "")
        offset = timedelta(hours=int(digits[:2]), minutes=int(digits[2:]))
        tz = timezone(-offset if zone[0] == "-" else offset)
    t = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second),
                 int((frac or "0")[:6].ljust(6, "0")), tzinfo=tz)
    return t.astimezone(timezone.utc)


def _read_agilent_d(contents: Path) -> datetime:
    """Agilent（MassHunter）の `AcqData/Contents.xml` の `AcquiredTime`（名前空間は問わない）。"""
    root = ET.parse(contents).getroot()
    stamp = next(((el.text or "").strip() for el in root.iter()
                  if el.tag.rsplit("}", 1)[-1] == "AcquiredTime"), "")
    if not stamp:
        raise ValueError("Agilent の AcqData/Contents.xml に AcquiredTime が無い")
    return _agilent_time(stamp)


def _read_d_folder(path: Path) -> tuple[str, datetime]:
    """.d フォルダの取得の開始時刻。(source, 開始)。Bruker（.tdf を先に）→ Agilent の順に見る。"""
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
            raise ValueError(f"Bruker の {db} に AcquisitionDateTime が無い")
        t = datetime.fromisoformat(str(row[0]))
        if t.tzinfo is None:  # 時刻帯が無ければ UTC とみなす
            t = t.replace(tzinfo=timezone.utc)
        return source, t.astimezone(timezone.utc)
    contents = path / "AcqData" / "Contents.xml"
    if contents.is_file():
        return "agilent_acqdata", _read_agilent_d(contents)
    raise ValueError(".d に Bruker の analysis.tdf / analysis.sqlite も Agilent の AcqData/Contents.xml も無い")


def _reason(exc: Exception) -> str:
    """読めなかった理由の文面。OSError の文字列はパスを含むので、その理由（`strerror`）だけにする
    （どのファイルかは行の `file` が持つ。理由ごとに件数を数える `unreadable_reasons` を 1 ファイル 1 種類にしない）。"""
    text = exc.strerror if isinstance(exc, OSError) else str(exc)
    return text or type(exc).__name__


def _wiff_sibling(path: Path) -> Path | None:
    """`.wiff2` と同じ幹の `.wiff`（SCIEX OS が同じ取得に並べて書く）。無ければ None。"""
    for suffix in (".wiff", ".WIFF"):
        sibling = path.with_suffix(suffix)
        if sibling.is_file():
            return sibling
    return None


def _read_wiff2(path: Path, unreadable: dict) -> dict:
    """SCIEX OS の `.wiff2` の取得時刻を、同じ幹の `.wiff`（OLE 複合ファイル）から読む。

    `.wiff2` 自身は中身の形式が公開されておらず、stdlib で読める取得時刻の欄が無い。SCIEX OS は同じ取得に
    `<幹>.wiff`・`<幹>.wiff.scan`・`<幹>.timeseries.data` も書くので、その `.wiff` の `SampleSubtree` を読む
    （`source` は `wiff_ole_sibling`）。

    理由（`error`）にファイル名を入れない: どのファイルかは行（`file`）が持ち、試料の表の `unreadable_reasons` は
    理由ごとに件数を数えるので、名前入りの理由は SCIEX OS のフォルダ（.wiff2 だけ）で 1 ファイル 1 種類になる。
    """
    sibling = _wiff_sibling(path)
    if sibling is None:
        return {**unreadable, "error": ("SCIEX OS の .wiff2 は stdlib で取得時刻を読めない形式で、同じ幹の .wiff"
                                        "（SCIEX OS が並べて書く）も無い（mzML に変換すれば startTimeStamp を読む）")}
    try:
        start, end = _read_wiff(sibling)
    except (OSError, ValueError, struct.error, IndexError) as exc:
        return {**unreadable, "error": f"同じ幹の .wiff を読めない: {_reason(exc)}"}
    if start is None:
        return {**unreadable, "error": "同じ幹の .wiff に開始時刻が無い"}
    return {"start_utc": start.isoformat(), "end_utc": end.isoformat() if end else None,
            "source": "wiff_ole_sibling", "error": None}


def read_acquisition_times(path: Path) -> dict:
    """1 つの raw の取得時刻。`{"start_utc", "end_utc", "source", "error"}`。

    `source` は `wiff_ole` / `wiff_ole_sibling`（`.wiff2` の時刻を同じ幹の `.wiff` から読んだ）/ `mzml` /
    `thermo_raw` / `bruker_tdf` / `bruker_baf` / `agilent_acqdata` / `unreadable`。読めないときは例外にせず
    `unreadable` と理由（`error`。形式ごとに何が無かったか。ファイル名・パスは入れない）を返す。時刻は UTC の ISO 8601。
    """
    path = Path(path)
    unreadable = {"start_utc": None, "end_utc": None, "source": _UNREADABLE}
    suffix = path.suffix.lower()
    if suffix == ".d" and path.is_dir():
        try:
            source, start = _read_d_folder(path)
        except (OSError, ValueError, sqlite3.Error, ET.ParseError) as exc:
            return {**unreadable, "error": _reason(exc)}
        return {"start_utc": start.isoformat(), "end_utc": None, "source": source, "error": None}
    if suffix == ".wiff2" and path.is_file():
        return _read_wiff2(path, unreadable)
    if suffix == ".raw" and path.is_file():
        try:
            start, _ = _read_thermo_raw(path)
        except (OSError, ValueError, struct.error) as exc:
            return {**unreadable, "error": _reason(exc)}
        if start is None:
            return {**unreadable, "error": "開始時刻が無い"}
        return {"start_utc": start.isoformat(), "end_utc": None, "source": "thermo_raw", "error": None}
    if suffix not in (".wiff", ".mzml") or not path.is_file():
        what = "フォルダ" if path.is_dir() else "形式"
        return {**unreadable, "error": f"取得時刻を読めない{what}: {path.suffix or '(拡張子なし)'}"}
    try:
        start, end = (_read_wiff if suffix == ".wiff" else _read_mzml)(path)
    except (OSError, ValueError, struct.error, ET.ParseError, IndexError) as exc:
        return {**unreadable, "error": _reason(exc)}
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


def session_gap_minutes(rows: list[dict], floor: float = GAP_MINUTES,
                        factor: float = SESSION_GAP_FACTOR) -> float:
    """取得のまとまりを区切る間隔（分）: `floor` と、隣どうしの取得の間隔の中央値の `factor` 倍の大きいほう。

    開始から次の開始までの間隔は LC の 1 周期（グラジエント・平衡化・オートサンプラ）より短くならないので、
    1 周期が `floor` を超える長いグラジエントのメソッドでは、固定の `floor` だと連続した取得の 1 本 1 本が
    まとまりに割れる。中央値は `statistics.median_low`（間隔が偶数個なら短いほう。中断が間隔の半分近くを
    占める少ない注入でも、中断の長さを 1 周期と読まない）。開始が `SAME_ACQUISITION_SECONDS` 秒より近い取得
    （ベンダーの raw と変換した mzML、`.wiff2` と同じ幹の `.wiff` のように同じ取得を指すファイル）は 1 つに数える
    （その間の 1 秒足らずの間隔で中央値を縮めない）。数える間隔が無ければ `floor`。
    """
    times = sorted(datetime.fromisoformat(r["start_utc"]) for r in rows if r.get("start_utc"))
    same = timedelta(seconds=SAME_ACQUISITION_SECONDS)
    intervals = [(b - a) / timedelta(minutes=1) for a, b in zip(times, times[1:]) if (b - a) >= same]
    if not intervals:
        return float(floor)
    return max(float(floor), factor * statistics.median_low(intervals))


def acquisition_sessions(rows: list[dict], minutes: float | None = None) -> dict[str, int]:
    """取得のまとまりの番号（1 始まり）。ファイル名 → 番号。

    読めた時刻を順に並べ、隣どうしの間隔が `minutes`（省けば `session_gap_minutes`: 40 分と、間隔の
    中央値の 3 倍の大きいほう）を超えるところで区切る。群の無いファイル（QC など）も並べるので、QC が
    間を埋める連続した取得は 1 つのまとまりになる。暦の日付で区切らないので、夜をまたぐ取得や UTC の
    日付の境目（UTC+9 なら朝 9 時）で割れない。時刻の読めないファイルは入らない。
    """
    if minutes is None:
        minutes = session_gap_minutes(rows)
    timed = sorted(((datetime.fromisoformat(r["start_utc"]), r["file"]) for r in rows if r.get("start_utc")))
    out: dict[str, int] = {}
    session, previous = 0, None
    for t, name in timed:
        if previous is None or (t - previous) > timedelta(minutes=minutes):
            session += 1
        out[name] = session
        previous = t
    return out


def _judge_pair(ranks: dict[str, list[int]], sessions: dict[int, Counter], a: str, b: str) -> str:
    """complete: 2 群の順位の範囲が重ならない。independent: 重なり、かつどちらかの群がいる取得のまとまり
    （`acquisition_sessions`）のすべてに両群がいる。それ以外は partial。"""
    ra, rb = ranks[a], ranks[b]
    if max(ra) < min(rb) or max(rb) < min(ra):
        return "complete"
    used = [k for k, c in sessions.items() if c[a] or c[b]]
    return "independent" if all(sessions[k][a] and sessions[k][b] for k in used) else "partial"


def group_overlap(rows: list[dict], group_of: dict[str, str]) -> dict:
    """群（`group_of`: ファイル名 → 群）と取得の順の重なりを、群の組ごとに判定する。

    組ごとの `overlap`: `complete`（2 群の取得順の範囲が重ならない＝取得順と群が完全に交絡）/
    `independent`（重なり、かつどちらかの群がいる取得のまとまりのすべてに両群がいる。まとまりは
    `session_gap_minutes`（40 分と、隣どうしの間隔の中央値の 3 倍の大きいほう）を超える間隔で区切る:
    `acquisition_sessions`）/ `partial`（それ以外）/ `unknown`（群のある試料に取得時刻の読めないものがある）。
    群の無いファイル（QC など）は判定に使わない（順位は全ファイルの中の順のまま。まとまりの区切りには使う）。
    全体の `overlap` は組の判定が 1 種類ならそれ、混ざれば `partial`、組が無ければ `unknown`。`by_session` は
    まとまりごとの最初と最後の取得時刻（UTC）と群ごとの件数、`session_gap_minutes` は区切った間隔（分）。
    """
    grouped = [r for r in rows if group_of.get(r["file"])]
    unreadable = sorted(r["file"] for r in grouped if r["rank"] is None)
    gap = session_gap_minutes(rows)
    session_of = acquisition_sessions(rows, gap)
    ranks: dict[str, list[int]] = {}
    sessions: dict[int, Counter] = {}
    for r in grouped:
        if r["rank"] is None:
            continue
        g = group_of[r["file"]]
        ranks.setdefault(g, []).append(r["rank"])
        sessions.setdefault(session_of[r["file"]], Counter())[g] += 1
    groups = sorted({group_of[r["file"]] for r in grouped})
    pairs = list(combinations(groups, 2))
    if unreadable:
        comparisons = [{"groups": [a, b], "overlap": "unknown"} for a, b in pairs]
    else:
        comparisons = [{"groups": [a, b], "overlap": _judge_pair(ranks, sessions, a, b)}
                       for a, b in pairs]
    states = {c["overlap"] for c in comparisons}
    overall = (next(iter(states)) if len(states) == 1 else "partial") if states else "unknown"
    bounds: dict[int, list[str]] = {}
    for r in rows:
        k = session_of.get(r["file"])
        if k is not None:
            bounds.setdefault(k, []).append(r["start_utc"])
    by_session = [{"session": k, "first_utc": min(bounds[k], key=datetime.fromisoformat),
                   "last_utc": max(bounds[k], key=datetime.fromisoformat),
                   "groups": dict(sorted(c.items()))}
                  for k, c in sorted(sessions.items())]
    return {"overlap": overall, "comparisons": comparisons,
            "group_rank_ranges": {g: [min(v), max(v)] for g, v in sorted(ranks.items())},
            "by_session": by_session, "session_gap_minutes": round(gap, 1), "unreadable_grouped": unreadable}
