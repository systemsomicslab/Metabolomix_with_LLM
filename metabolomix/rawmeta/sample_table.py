"""ISA の試料の表・mwTab から、因子の列付きの実験情報シート（sample-manifest.v2）を組む（stdlib だけ）。

- 群は試料の表の**指定した因子の列の値の組**（値を `_` でつなぎ、空白と記号は `_`）。論文の比較の選び方は持ち込まない。
  因子の値が空の試料と、表に行の無い試料の群は `unassigned`（黙って別の群へ落とさない）。
- 列名は完全一致 → 大文字小文字を無視した一致 → `Factor Value[<名>]` → `Characteristics[<名>]` の順に引く。
- 試料の表は ISA の TSV（`s_*.txt`）か Metabolomics Workbench の mwTab（`SUBJECT_SAMPLE_FACTORS`）。既定の
  `auto` は 1 行目が `#METABOLOMICS WORKBENCH` なら mwTab。
- ファイルと試料の対応は ISA の assay の表（`a_*.txt`）の `Sample Name` と `Raw Spectral Data File`、mwTab なら
  `RAW_FILE_NAME`。どちらも無ければファイルの幹を試料名とみなす。
- 注入順は raw の取得時刻の順（`acquisition.read_acquisition_times`）で、schema の行の直後に宣言の行
  `# injection_order_source = raw_acquisition_time` を書く。1 つでも読めなければ注入順を空にし、宣言も書かない
  （ファイル名の順は注入順とは限らない）。バッチは分からないので空のまま。

reanalysis-study の `reanalysis/inputs.py`（blind-inputs）から移した。ハードリンクで極性ごとの入力フォルダを
作る部分は移していない（pipeline の受付が raw を run の中へ配置する）。
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.rawmeta.acquisition import acquisition_order

__all__ = ["MANIFEST_COLUMNS", "SAMPLE_TABLE_FORMATS", "UNASSIGNED", "build_sample_sheet",
           "factor_group", "read_mwtab_samples", "read_table", "resolve_column"]

_CODE = "SAMPLE_TABLE_INVALID"
SAMPLE_TABLE_FORMATS = ("auto", "isa", "mwtab")
MWTAB_MAGIC = "#METABOLOMICS WORKBENCH"
UNASSIGNED = "unassigned"
QC_GROUP = "QC"
BLANK_GROUP = "blank"
#: sample-manifest.v2 の列（`sample_manifest.FIELDS_V2` と同じ順。テストが突き合わせる）。
MANIFEST_COLUMNS = ("sample_id", "source_file", "role", "group", "batch", "injection_order",
                    "qc_pool", "include", "biological_sample_id")
#: 注入順を取得時刻から作ったシートの宣言の行（`sample_manifest.ORDER_SOURCE_DECLARATION`。schema の行の直後）。
#: 読み手（`parse_manifest`）は注入順の出所に `derived_from: raw_acquisition_time` を足し、試料の表・統計は
#: 出所を `raw_acquisition_time` と示す。この層は sample_manifest を import しない（テストが突き合わせる）。
ORDER_SOURCE_LINE = "# injection_order_source = raw_acquisition_time"
#: 因子の列の名前に使えない名前（v2 の列・任意の subject_id・行の帳簿の欄）。
_RESERVED = frozenset([*MANIFEST_COLUMNS, "subject_id", "provenance", "conflicts"])


def _err(message: str, **details) -> DomainError:
    return DomainError(_CODE, message, details or None)


def _norm(value: object) -> str:
    return str(value).strip()


def _unreadable(path: Path, exc: Exception) -> DomainError:
    """文字コード（UTF-8 でない）・TSV の形（csv.Error）で読めない表の SAMPLE_TABLE_INVALID。"""
    if isinstance(exc, UnicodeDecodeError):
        return _err(f"{path.name} は UTF-8 として読めない（{exc.reason}。UTF-8 で保存し直す）",
                    path=str(path), encoding="utf-8")
    return _err(f"{path.name} を TSV として読めない（{exc}）", path=str(path))


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """ISA の TSV を（見出し, 行）にする。行は見出しの列数に足りなければ空文字で埋める。

    UTF-8（BOM 可）でない・TSV として読めない表は SAMPLE_TABLE_INVALID。"""
    path = Path(path)
    if not path.is_file():
        raise _err(f"試料の表が無い: {path.name}", path=str(path))
    try:
        with path.open(encoding="utf-8-sig", newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd, None)
            if header is None:
                raise _err(f"試料の表が空: {path.name}", path=str(path))
            rows = [r + [""] * (len(header) - len(r)) for r in rd if any(c.strip() for c in r)]
    except (UnicodeDecodeError, csv.Error) as exc:
        raise _unreadable(path, exc) from exc
    return [_norm(h) for h in header], rows


def _is_mwtab(path: Path) -> bool:
    try:
        with Path(path).open(encoding="utf-8-sig", newline="") as fh:
            return fh.readline().startswith(MWTAB_MAGIC)
    except UnicodeDecodeError as exc:
        raise _unreadable(Path(path), exc) from exc


def read_mwtab_samples(path: Path) -> tuple[list[str], list[list[str]], dict[str, str]]:
    """mwTab の `SUBJECT_SAMPLE_FACTORS` の行を試料の表にする。

    行はタブ区切りで（行名・SUBJECT・SAMPLE・FACTORS・追加の欄）。見出しは `Sample Name`・`Subject`・因子の名前
    （初出順）・追加の欄の鍵（初出順）。因子は `|` で分け最初の `:` で名前と値に、追加の欄は `;` で分け最初の `=` で
    鍵と値に分ける。戻り値の 3 つ目は `RAW_FILE_NAME` の対応（ファイル名と、曖昧でなければその幹 → Sample Name）。
    """
    path = Path(path)
    parsed = []
    try:
        text = path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as exc:
        raise _unreadable(path, exc) from exc
    for line in text.splitlines():
        if line.startswith("#") or not line.startswith("SUBJECT_SAMPLE_FACTORS"):
            continue
        cells = line.split("\t")
        cells += [""] * (5 - len(cells))
        subject, sample, factors_text, extra_text = (c.strip() for c in cells[1:5])
        if not sample:
            continue
        factors = {}
        for part in factors_text.split("|"):
            if ":" in part:
                k, v = part.split(":", 1)
                if k.strip():
                    factors[k.strip()] = v.strip()
        extra = {}
        for part in extra_text.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                if k.strip():
                    extra[k.strip()] = v.strip()
        parsed.append((subject, sample, factors, extra))
    if not parsed:
        raise _err(f"{path.name} に SUBJECT_SAMPLE_FACTORS の行が 1 つも無い（mwTab の試料の表として読めない）",
                   path=str(path))
    factor_names = list(dict.fromkeys(k for _s, _n, f, _e in parsed for k in f))
    extra_names = list(dict.fromkeys(k for _s, _n, _f, e in parsed for k in e))
    header = ["Sample Name", "Subject", *factor_names, *extra_names]
    rows = [[sample, subject, *(f.get(k, "") for k in factor_names), *(e.get(k, "") for k in extra_names)]
            for subject, sample, f, e in parsed]
    raw_map: dict[str, str] = {}
    stems: dict[str, set[str]] = {}
    for _subject, sample, _f, e in parsed:
        fname = re.split(r"[\\/]", e.get("RAW_FILE_NAME", ""))[-1].strip()
        if not fname:
            continue
        if raw_map.setdefault(fname, sample) != sample:
            raise _err(f"{path.name} の RAW_FILE_NAME {fname!r} が複数の試料（{raw_map[fname]!r}・{sample!r}）にある",
                       raw_file=fname)
        stems.setdefault(Path(fname).stem, set()).add(sample)
    for stem, samples in stems.items():   # 幹の鍵は名前の鍵と重ならず、1 つの試料だけのものに限る
        if stem not in raw_map and len(samples) == 1:
            raw_map[stem] = next(iter(samples))
    return header, rows, raw_map


def resolve_column(header: list[str], name: str) -> int:
    """列名を見出しの位置にする（完全一致 → 大文字小文字を無視 → `Factor Value[名]` → `Characteristics[名]`）。"""
    want = _norm(name)
    if want in header:
        return header.index(want)
    lowered = [h.lower() for h in header]
    for candidate in (want.lower(), f"factor value[{want.lower()}]", f"characteristics[{want.lower()}]"):
        if candidate in lowered:
            return lowered.index(candidate)
    raise _err(f"試料の表に列 {name!r} が無い（Factor Value[{name}]・Characteristics[{name}] も無い）",
               column=name, columns=header[:50])


def _slug(text: str) -> str:
    return re.sub(r"[^\w]+", "_", text).strip("_")


def factor_group(values: list[str]) -> str:
    """因子の値の組から群の名前を作る（値を `_` でつなぎ、空白と記号は `_`）。空の値が 1 つでもあれば `unassigned`。"""
    vals = [_norm(v) for v in values]
    if not vals or not all(vals):
        return UNASSIGNED
    return _slug("_".join(vals)) or UNASSIGNED


def _breaks_tsv(s: str) -> bool:
    """シートの TSV で崩れる文字列か（タブ・行に割る文字・先頭の引用符）。"""
    return "\t" in s or len(("|" + s + "|").splitlines()) > 1 or s.startswith('"')


def _check_factor_names(factor_columns: list[str]) -> list[str]:
    if not factor_columns:
        raise _err("factor_columns に列名を 1 つ以上渡す")
    names = [_norm(c) for c in factor_columns]
    reserved = {r.casefold() for r in _RESERVED}
    for n in names:
        if not n:
            raise _err("factor_columns に空の列名がある")
        if _breaks_tsv(n):
            raise _err(f"因子の列名 {n!r} にタブ・改行・先頭の引用符がある（シートの列が崩れる）", column=n)
        if n.casefold() in reserved:
            raise _err(f"因子の列名 {n!r} はシートの既知の列名と重なる（試料の表の別の列を選ぶ）",
                       column=n, reserved=sorted(_RESERVED))
    folded = [n.casefold() for n in names]
    dup = sorted({n for n, f in zip(names, folded) if folded.count(f) > 1})
    if dup:
        raise _err(f"factor_columns の列名が重複している: {dup}", duplicated=dup)
    return names


def _load_sample_table(path: Path, fmt: str) -> tuple[str, list[str], list[list[str]], dict[str, str]]:
    if fmt not in SAMPLE_TABLE_FORMATS:
        raise _err(f"sample_table_format は {list(SAMPLE_TABLE_FORMATS)} のどれか: {fmt!r}", value=fmt)
    path = Path(path)
    if not path.is_file():
        raise _err(f"試料の表が無い: {path.name}", path=str(path))
    resolved = ("mwtab" if _is_mwtab(path) else "isa") if fmt == "auto" else fmt
    if resolved == "mwtab":
        header, rows, raw_map = read_mwtab_samples(path)
    else:
        header, rows = read_table(path)
        raw_map = {}
    return resolved, header, rows, raw_map


def _assay_map(path: Path) -> dict[str, str]:
    """ISA の assay の表の `Raw Spectral Data File`（フォルダを除いた名前）→ `Sample Name`。"""
    header, rows = read_table(path)
    if "Sample Name" not in header:
        raise _err(f"{Path(path).name} に列 Sample Name が無い", path=str(path))
    i_name = header.index("Sample Name")
    raw_cols = [i for i, h in enumerate(header) if h == "Raw Spectral Data File"]
    if not raw_cols:
        raise _err(f"{Path(path).name} に列 Raw Spectral Data File が無い", path=str(path))
    mapping: dict[str, str] = {}
    for r in rows:
        for i in raw_cols:
            fname = re.split(r"[\\/]", r[i].strip())[-1]
            if fname and r[i_name].strip():
                mapping.setdefault(fname, r[i_name].strip())
    return mapping


def _compile_selector(selector, label: str, header: list[str]):
    """（列名, 正規表現）を読む。戻り値は（列名, 正規表現, Sample Name の列か, 列の位置）。"""
    if (not isinstance(selector, (list, tuple)) or len(selector) != 2
            or not all(isinstance(x, str) and x.strip() for x in selector)):
        raise _err(f"{label} は [列名, 正規表現] で、どちらも空にしない: {selector!r}", field=label)
    col, pattern = _norm(selector[0]), selector[1]
    try:
        rx = re.compile(pattern)
    except re.error as exc:
        raise _err(f"{label} の正規表現 {pattern!r} が読めない（{exc}）", field=label) from exc
    on_name = col == "Sample Name"
    return col, rx, on_name, None if on_name else resolve_column(header, col)


def _hit(compiled, sample: str, trow: list[str] | None) -> bool:
    if compiled is None:
        return False
    _col, rx, on_name, idx = compiled
    value = sample if on_name else (trow[idx] if trow is not None else None)
    return value is not None and rx.search(value.strip()) is not None


def _sample_id(name: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", name).strip("_") or "sample"


def build_sample_sheet(raw_dir: Path, raw_names: list[str], sample_table: Path, factor_columns: list[str],
                       *, assay_table: Path | None = None, sample_table_format: str = "auto",
                       qc_selector=None, blank_selector=None) -> dict:
    """raw のフォルダの主ファイル（`raw_names`）ごとに 1 行の sample-manifest.v2 を組む（書かない）。

    戻り値: `{"text"（シートの全文）, "rows"（行の辞書）, "columns", "injection_order_source"
    （`raw_acquisition_time` / `none`）, "acquisition"（ファイルごとの取得時刻と順）, "unreadable",
    "groups"（群 → 件数）, "qc_samples", "blank_samples", "unassigned_samples", "files_not_in_assay",
    "sample_table_format", "factor_columns"（試料の表の列名）}`。誤りは `SAMPLE_TABLE_INVALID`。
    """
    raw_dir = Path(raw_dir)
    factor_names = _check_factor_names(list(factor_columns or []))
    fmt, header, table, raw_map = _load_sample_table(Path(sample_table), sample_table_format)
    i_sample = resolve_column(header, "Sample Name")
    by_name: dict[str, list[str]] = {}
    for r in table:
        by_name.setdefault(r[i_sample].strip(), r)
    factor_idx = [resolve_column(header, c) for c in factor_names]
    if len(set(factor_idx)) != len(factor_idx):
        raise _err("factor_columns に試料の表の同じ列に当たるものが複数ある",
                   columns=[header[i] for i in factor_idx])
    qc = _compile_selector(qc_selector, "qc_selector", header) if qc_selector is not None else None
    blank = (_compile_selector(blank_selector, "blank_selector", header)
             if blank_selector is not None else None)
    assay = _assay_map(Path(assay_table)) if assay_table is not None else raw_map

    entries, files_not_in_assay = [], []
    origins: dict[str, object] = {}
    for name in raw_names:
        stem = Path(name).stem
        sample = assay.get(name)
        if sample is None and fmt == "mwtab" and assay_table is None:
            sample = assay.get(stem)
        if sample is None:
            if assay:
                files_not_in_assay.append(name)
            sample = stem
        trow = by_name.get(sample)
        is_qc, is_blank = _hit(qc, sample, trow), _hit(blank, sample, trow)
        if is_qc and is_blank:
            raise _err(f"試料 {sample!r} は qc_selector にも blank_selector にも当たる（選び方を確かめる）",
                       sample=sample)
        empty = [""] * len(factor_idx)
        if is_qc:
            role, group, origin, factors = "qc", QC_GROUP, "qc_selector に当たった試料", empty
        elif is_blank:
            role, group, origin, factors = "blank", BLANK_GROUP, "blank_selector に当たった試料", empty
        else:
            role = "sample"
            values = tuple(_norm(trow[i]) for i in factor_idx) if trow is not None else ()
            factors = list(values) if values else empty
            for fname, v in zip(factor_names, factors):
                if _breaks_tsv(v):
                    raise _err(f"試料 {sample!r} の因子の列 {fname!r} の値にタブ・改行・先頭の引用符がある",
                               sample=sample, column=fname)
            group = factor_group(list(values)) if trow is not None else UNASSIGNED
            origin = values if values and all(values) else "因子の値が空か、表に行の無い試料"
        first = origins.setdefault(group, origin)
        if first != origin:
            raise _err(f"群の名前 {group!r} が別々の値の組から作られる（{first} と {origin}。群が混ざる）",
                       group=group)
        if _breaks_tsv(sample):
            raise _err(f"試料名 {sample!r} にタブ・改行・先頭の引用符がある", sample=sample)
        entries.append({"file": name, "sample": sample, "group": group, "role": role, "factors": factors})

    acq = acquisition_order([raw_dir / e["file"] for e in entries])
    by_file = {row["file"]: row for row in acq}
    unreadable = sorted(r["file"] for r in acq if r["rank"] is None)
    if unreadable:
        ordered, source = sorted(entries, key=lambda e: e["file"]), "none"
    else:
        ordered = sorted(entries, key=lambda e: by_file[e["file"]]["rank"])
        source = "raw_acquisition_time"

    rows, used = [], set()
    for rank, e in enumerate(ordered, start=1):
        sid = _sample_id(e["sample"])
        if sid in used:          # 同じ試料の複数ファイル（技術反復など）はファイルの幹で見分ける
            sid = _sample_id(Path(e["file"]).stem)
        if sid in used:
            raise _err(f"試料 ID {sid!r} が重複する（試料名とファイルの幹のどちらでも見分けられない）",
                       sample_id=sid)
        used.add(sid)
        row = {"sample_id": sid, "source_file": e["file"], "role": e["role"], "group": e["group"],
               "batch": "", "injection_order": str(rank) if source == "raw_acquisition_time" else "",
               "qc_pool": QC_GROUP if e["role"] == "qc" else "", "include": "true",
               "biological_sample_id": e["sample"]}
        for fname, value in zip(factor_names, e["factors"]):
            row[fname] = value
        rows.append(row)
    columns = [*MANIFEST_COLUMNS, *factor_names]
    lines = ["# schema = sample-manifest.v2"]
    if source == "raw_acquisition_time":
        # 注入順が raw の取得時刻の順であることをシートに残す（読み手が試料の表の出所に出す）
        lines.append(ORDER_SOURCE_LINE)
    lines.append("\t".join(columns))
    lines += ["\t".join(row[c] for c in columns) for row in rows]
    groups: dict[str, int] = {}
    for row in rows:
        groups[row["group"]] = groups.get(row["group"], 0) + 1
    return {
        "text": "\n".join(lines) + "\n", "rows": rows, "columns": columns,
        "injection_order_source": source,
        "acquisition": [{"file": r["file"], "start_utc": r["start_utc"], "rank": r["rank"],
                         "source": r["source"], "error": r.get("error")} for r in acq],
        "unreadable": unreadable, "groups": dict(sorted(groups.items())),
        "qc_samples": [e["sample"] for e in entries if e["role"] == "qc"],
        "blank_samples": [e["sample"] for e in entries if e["role"] == "blank"],
        "unassigned_samples": [e["sample"] for e in entries
                               if e["role"] == "sample" and e["group"] == UNASSIGNED],
        "files_not_in_assay": files_not_in_assay, "sample_table_format": fmt,
        "factor_columns": [header[i] for i in factor_idx],
    }
