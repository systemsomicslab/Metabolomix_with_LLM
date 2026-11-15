"""ISA の試料の表・mwTab から、因子の列付きの実験情報シート（sample-manifest.v2）を組む（stdlib だけ）。

- 群は試料の表の**指定した因子の列の値の組**（値を `_` でつなぎ、空白と記号は `_`）。論文の比較の選び方は持ち込まない。
  因子の値が空の試料・表に行の無い試料・対応が決まらない試料（下記）の群は `unassigned` で、`include=false`
  （統計の母集団から外す。黙って別の群へ落とさず、`unassigned` という群にもしない）。
- 列名は完全一致 → 大文字小文字を無視した一致 → `Factor Value[<名>]` → `Characteristics[<名>]` の順に引く。
  因子の列はシートの後ろに、渡した名前で足す。シートの既知の列（`group`・`batch` など）と**完全に同じ**名前だけは
  `factor_<名>` に改める（`Group` のように大文字小文字だけが違う名前はそのまま。読み手も統計も列名を区別する）。
- 試料の表は ISA の TSV（`s_*.txt`）か Metabolomics Workbench の mwTab（`SUBJECT_SAMPLE_FACTORS`）。既定の
  `auto` は 1 行目が `#METABOLOMICS WORKBENCH` なら mwTab。
- ファイルと試料の対応は ISA の assay の表（`a_*.txt`）の `Sample Name` と `Raw Spectral Data File`・
  `Derived Spectral Data File`（変換した mzML はこちらに載る）、mwTab なら `RAW_FILE_NAME`。名前が合わなければ
  幹（拡張子を除いた名前）で引く（`.raw` と書かれた試料を `.mzML` で渡した場合など）。どちらも無ければファイルの幹を
  試料名とみなす。1 つのファイル（幹）が別々の試料に当たる・同じ Sample Name の行どうしで因子の値が食い違う
  ときは、どちらの群にも入れず `unassigned` にして返す（ISA も mwTab も同じ扱い）。
- 生体試料（role=sample）のどれも群が決まらないシート（表の行に 1 つも当たらない・当たった行の因子の値が
  すべて空）は組まずに `SAMPLE_TABLE_INVALID`（全試料を統計から外すシートを既定の名前で置かない）。
- 注入順は raw の取得時刻の順（`acquisition.read_acquisition_times`）で、schema の行の直後に宣言の行
  `# injection_order_source = raw_acquisition_time` を書く。読めない raw があっても読めた raw の順は残し
  （読めたものの中の順位）、読めない raw の注入順だけを空にする（位置が分からないものを推して埋めない。ファイル名の
  順は注入順とは限らない）。読み手（`sample_manifest.parse_manifest`）はその空欄を「不明」として印を付け、
  MS-DIAL の番号（`.mddata`・mzTab-M）で埋めない（番号の付け方が違う値を混ぜない）。どれも読めなければ注入順は
  すべて空で、宣言も書かない。バッチは分からないので空のまま。

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
#: `build_sample_sheet` の `injection_order_source`: 全 raw の取得時刻が読めた / 一部だけ読めた（読めない raw の
#: 注入順は空）/ どれも読めない（`none`。宣言の行も書かない）。
ORDER_FROM_ACQUISITION = "raw_acquisition_time"
ORDER_FROM_ACQUISITION_PARTIAL = "raw_acquisition_time_partial"
#: 因子の列にそのままの名前では書けない名前（v2 の列・任意の subject_id・行の帳簿の欄）。完全一致だけを見る
#: （`sample_manifest._check_factor_columns` と同じ）。当たれば `factor_<名>` に改めて書く。
_RESERVED = frozenset([*MANIFEST_COLUMNS, "subject_id", "provenance", "conflicts"])
#: ISA の assay の表でファイルと試料を結ぶ列（MetaboLights は変換した mzML を Derived に載せる）。
FILE_COLUMNS = ("Raw Spectral Data File", "Derived Spectral Data File")


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


def read_mwtab_samples(path: Path) -> tuple[list[str], list[list[str]], list[tuple[str, str]]]:
    """mwTab の `SUBJECT_SAMPLE_FACTORS` の行を試料の表にする。

    行はタブ区切りで（行名・SUBJECT・SAMPLE・FACTORS・追加の欄）。見出しは `Sample Name`・`Subject`・因子の名前
    （初出順）・追加の欄の鍵（初出順）。因子は `|` で分け最初の `:` で名前と値に、追加の欄は `;` で分け最初の `=` で
    鍵と値に分ける。戻り値の 3 つ目は `RAW_FILE_NAME`（フォルダを除いた名前）と Sample Name の組の並び
    （対応は `_file_index` が作る。同じファイルが複数の試料にあっても止めない）。
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
    pairs = [(_file_name(e.get("RAW_FILE_NAME", "")), sample) for _subject, sample, _f, e in parsed]
    return header, rows, [(f, sample) for f, sample in pairs if f]


def _file_name(text: str) -> str:
    """表に書かれたファイルのパスから、フォルダを除いた名前（`\\` と `/` のどちらの区切りも外す）。"""
    return re.split(r"[\\/]", text.strip())[-1].strip()


def _file_index(pairs: list[tuple[str, str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    """（ファイル名, 試料名）の組から、ファイル名 → 試料名の対応と、別々の試料に当たる名前。

    幹（拡張子を除いた名前）の鍵も足す（名前の鍵と重ならないものだけ）。1 つの名前・幹が別々の試料に当たれば
    対応に入れず 2 つ目に返す（どちらの群にも黙って入れない。ISA と mwTab で同じ扱い）。同じ試料に何度
    当たっても曖昧ではない（`.raw` と変換した `.mzML` を同じ試料に並べる assay の表）。
    """
    by_name: dict[str, dict[str, None]] = {}
    for fname, sample in pairs:
        if fname and sample:
            by_name.setdefault(fname, {})[sample] = None
    by_stem: dict[str, dict[str, None]] = {}
    for fname, samples in by_name.items():
        stem = Path(fname).stem
        if stem != fname and stem not in by_name:
            by_stem.setdefault(stem, {}).update(samples)
    mapping: dict[str, str] = {}
    ambiguous: dict[str, list[str]] = {}
    for key, samples in (*by_name.items(), *by_stem.items()):
        if len(samples) == 1:
            mapping[key] = next(iter(samples))
        else:
            ambiguous[key] = list(samples)
    return mapping, ambiguous


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
    for n in names:
        if not n:
            raise _err("factor_columns に空の列名がある")
        if _breaks_tsv(n):
            raise _err(f"因子の列名 {n!r} にタブ・改行・先頭の引用符がある（シートの列が崩れる）", column=n)
    folded = [n.casefold() for n in names]
    dup = sorted({n for n, f in zip(names, folded) if folded.count(f) > 1})
    if dup:
        raise _err(f"factor_columns の列名が重複している: {dup}", duplicated=dup)
    return names


def _sheet_factor_names(names: list[str]) -> tuple[list[str], dict[str, str]]:
    """シートに書く因子の列名と、改めた名前（渡した名前 → シートの列名）。

    シートの既知の列（`_RESERVED`）と完全に同じ名前だけを `factor_<名>`（それも塞がっていれば `_2` から番号）に
    改める。大文字小文字だけが違う名前（mwTab の因子 `Group` と列 `group`）はそのまま書く——読み手
    （`sample_manifest._check_factor_columns`）も統計の計画（列名で引く）も完全一致で区別する。
    """
    taken = set(_RESERVED) | set(names)
    out: list[str] = []
    renamed: dict[str, str] = {}
    for n in names:
        if n not in _RESERVED:
            out.append(n)
            continue
        new, k = f"factor_{n}", 2
        while new in taken:
            new, k = f"factor_{n}_{k}", k + 1
        taken.add(new)
        renamed[n] = new
        out.append(new)
    return out, renamed


def _load_sample_table(path: Path, fmt: str) -> tuple[str, list[str], list[list[str]], list[tuple[str, str]]]:
    if fmt not in SAMPLE_TABLE_FORMATS:
        raise _err(f"sample_table_format は {list(SAMPLE_TABLE_FORMATS)} のどれか: {fmt!r}", value=fmt)
    path = Path(path)
    if not path.is_file():
        raise _err(f"試料の表が無い: {path.name}", path=str(path))
    resolved = ("mwtab" if _is_mwtab(path) else "isa") if fmt == "auto" else fmt
    if resolved == "mwtab":
        header, rows, pairs = read_mwtab_samples(path)
    else:
        header, rows = read_table(path)
        pairs = []
    return resolved, header, rows, pairs


def _assay_map(path: Path) -> list[tuple[str, str]]:
    """ISA の assay の表の（ファイル名, `Sample Name`）の組。ファイル名は `FILE_COLUMNS` の列（同じ見出しが何列
    あってもすべて）の値からフォルダを除いた名前。対応にするのは `_file_index`。"""
    header, rows = read_table(path)
    if "Sample Name" not in header:
        raise _err(f"{Path(path).name} に列 Sample Name が無い", path=str(path))
    i_name = header.index("Sample Name")
    file_cols = [i for i, h in enumerate(header) if h in FILE_COLUMNS]
    if not file_cols:
        raise _err(f"{Path(path).name} に列 Raw Spectral Data File も Derived Spectral Data File も無い",
                   path=str(path), file_columns=list(FILE_COLUMNS))
    return [(_file_name(r[i]), r[i_name].strip()) for r in rows for i in file_cols]


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


def _table_rows(table: list[list[str]], i_sample: int, compared: list[int],
                header: list[str]) -> tuple[dict[str, list[str]], list[dict]]:
    """Sample Name → 試料の表の行と、同じ Sample Name の行どうしで `compared` の列の値が食い違う試料。

    値が同じ重複（1 つの試料が複数の Source Name から来る ISA の表）は最初の行を使う。食い違えば群を決められない
    ので、2 つ目の戻り値（`{"sample", "columns", "values"}`）に載せ、呼び出し側はその試料を表に行の無い試料と
    同じに扱う（ISA も mwTab も同じ）。
    """
    by_name: dict[str, list[str]] = {}
    values: dict[str, dict[tuple[str, ...], None]] = {}
    for r in table:
        key = r[i_sample].strip()
        by_name.setdefault(key, r)
        values.setdefault(key, {})[tuple(_norm(r[i]) for i in compared)] = None
    columns = [header[i] for i in compared]
    conflicts = [{"sample": key, "columns": columns, "values": [list(v) for v in seen]}
                 for key, seen in values.items() if len(seen) > 1]
    return by_name, conflicts


def _unmatched_error(entries: list[dict], fmt: str, assay_table, by_name: dict) -> DomainError:
    """生体試料のどれも群が決まらないときの SAMPLE_TABLE_INVALID（何が合わなかったかを details に）。"""
    bio = [e for e in entries if e["role"] == "sample"]
    no_row = [e["file"] for e in bio if not e["matched"]]
    if len(no_row) == len(bio):
        where = ("assay の表の Raw / Derived Spectral Data File の名前（か幹）と Sample Name" if assay_table is not None
                 else "RAW_FILE_NAME の名前（か幹）と Sample Name" if fmt == "mwtab"
                 else "ファイルの幹と Sample Name（ISA なら assay_table を渡すとファイル名で結べる）")
        message = (f"raw の {len(bio)} 件のどれも試料の表の行に当たらない（{where} を確かめる）。"
                   "全試料を群 unassigned・include=false にするシートは書かない")
    else:
        message = ("試料の表の行に当たった試料のどれも、選んだ因子の列の値が空か、同じ試料の行どうしで食い違う"
                   "（factor_columns を確かめる）。全試料を群 unassigned・include=false にするシートは書かない")
    return _err(message, n_samples=len(bio), files_without_table_row=no_row[:20],
                n_files_without_table_row=len(no_row), table_sample_names=list(by_name)[:10])


def build_sample_sheet(raw_dir: Path, raw_names: list[str], sample_table: Path, factor_columns: list[str],
                       *, assay_table: Path | None = None, sample_table_format: str = "auto",
                       qc_selector=None, blank_selector=None) -> dict:
    """raw のフォルダの主ファイル（`raw_names`）ごとに 1 行の sample-manifest.v2 を組む（書かない）。

    戻り値: `{"text"（シートの全文）, "rows"（行の辞書）, "columns", "injection_order_source"
    （`raw_acquisition_time` / `raw_acquisition_time_partial`（読めない raw の注入順だけ空）/ `none`）,
    "acquisition"（ファイルごとの取得時刻と順）, "unreadable"（取得時刻の読めない raw。注入順が空）,
    "unreadable_reasons"（読めなかった理由 → 件数）, "groups"（群 → 件数）, "qc_samples", "blank_samples",
    "unassigned_samples"（群 unassigned・include=false の生体試料）, "files_not_in_assay", "ambiguous_files"
    （別々の試料に当たるファイル: `{"file", "samples"}`）, "conflicting_samples"（行どうしで因子の値が食い違う試料:
    `{"sample", "columns", "values"}`）, "sample_table_format", "factor_columns"（試料の表の列名）,
    "renamed_factor_columns"（シートの既知の列と同じ名前を改めたもの: 渡した名前 → シートの列名）}`。
    誤りは `SAMPLE_TABLE_INVALID`（生体試料のどれも群が決まらないときも）。
    """
    raw_dir = Path(raw_dir)
    factor_names = _check_factor_names(list(factor_columns or []))
    sheet_names, renamed = _sheet_factor_names(factor_names)
    fmt, header, table, pairs = _load_sample_table(Path(sample_table), sample_table_format)
    i_sample = resolve_column(header, "Sample Name")
    factor_idx = [resolve_column(header, c) for c in factor_names]
    if len(set(factor_idx)) != len(factor_idx):
        raise _err("factor_columns に試料の表の同じ列に当たるものが複数ある",
                   columns=[header[i] for i in factor_idx])
    qc = _compile_selector(qc_selector, "qc_selector", header) if qc_selector is not None else None
    blank = (_compile_selector(blank_selector, "blank_selector", header)
             if blank_selector is not None else None)
    # 役割（qc / blank）を決める列も、因子と同じく重複した行どうしで食い違えば決められない
    compared = list(dict.fromkeys([*factor_idx, *(c[3] for c in (qc, blank) if c is not None and c[3] is not None)]))
    by_name, conflicting = _table_rows(table, i_sample, compared, header)
    conflicted = {c["sample"] for c in conflicting}
    if assay_table is not None:
        pairs = _assay_map(Path(assay_table))
    mapping, ambiguous = _file_index(pairs)

    entries, files_not_in_assay, ambiguous_files = [], [], []
    origins: dict[str, object] = {}
    for name in raw_names:
        stem = Path(name).stem
        key = next((k for k in (name, stem) if k in mapping or k in ambiguous), None)
        if key in ambiguous:
            ambiguous_files.append({"file": name, "samples": ambiguous[key]})
            sample, trow = stem, None      # どの試料の行かが決まらない（どちらの群にも入れない）
        else:
            sample = mapping.get(key) if key is not None else None
            if sample is None:
                if mapping or ambiguous:
                    files_not_in_assay.append(name)
                sample = stem
            trow = None if sample in conflicted else by_name.get(sample)
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
            origin = values if values and all(values) else "因子の値が空か、表に行の無い・対応の決まらない試料"
        first = origins.setdefault(group, origin)
        if first != origin:
            raise _err(f"群の名前 {group!r} が別々の値の組から作られる（{first} と {origin}。群が混ざる）",
                       group=group)
        if _breaks_tsv(sample):
            raise _err(f"試料名 {sample!r} にタブ・改行・先頭の引用符がある", sample=sample)
        entries.append({"file": name, "sample": sample, "group": group, "role": role, "factors": factors,
                        "matched": trow is not None})

    used_samples = {e["sample"] for e in entries}
    conflicting = [c for c in conflicting if c["sample"] in used_samples]   # raw の当たる試料だけを返す
    bio = [e for e in entries if e["role"] == "sample"]
    if bio and all(e["group"] == UNASSIGNED for e in bio):
        raise _unmatched_error(entries, fmt, assay_table, by_name)

    acq = acquisition_order([raw_dir / e["file"] for e in entries])
    rank_of = {row["file"]: row["rank"] for row in acq}
    unreadable = sorted(r["file"] for r in acq if r["rank"] is None)
    reasons: dict[str, int] = {}
    for r in acq:
        if r["rank"] is None:
            reason = r.get("error") or "理由不明"
            reasons[reason] = reasons.get(reason, 0) + 1
    # 読めた raw は取得時刻の順（読めたものの中の順位）、読めない raw は注入順を空にして後ろへ（位置は分からない）
    ordered = sorted(entries, key=lambda e: (rank_of[e["file"]] is None, rank_of[e["file"]] or 0, e["file"]))
    n_ranked = len(entries) - len(unreadable)
    source = ("none" if n_ranked == 0 else ORDER_FROM_ACQUISITION if not unreadable
              else ORDER_FROM_ACQUISITION_PARTIAL)

    rows, used = [], set()
    for e in ordered:
        rank = rank_of[e["file"]]
        sid = _sample_id(e["sample"])
        if sid in used:          # 同じ試料の複数ファイル（技術反復など）はファイルの幹で見分ける
            sid = _sample_id(Path(e["file"]).stem)
        if sid in used:
            raise _err(f"試料 ID {sid!r} が重複する（試料名とファイルの幹のどちらでも見分けられない）",
                       sample_id=sid)
        used.add(sid)
        # 群の決まらない生体試料は統計の母集団から外す（`unassigned` を 1 つの群として数えさせない）
        excluded = e["role"] == "sample" and e["group"] == UNASSIGNED
        row = {"sample_id": sid, "source_file": e["file"], "role": e["role"], "group": e["group"],
               "batch": "", "injection_order": "" if rank is None else str(rank),
               "qc_pool": QC_GROUP if e["role"] == "qc" else "", "include": "false" if excluded else "true",
               "biological_sample_id": e["sample"]}
        for column, value in zip(sheet_names, e["factors"]):
            row[column] = value
        rows.append(row)
    columns = [*MANIFEST_COLUMNS, *sheet_names]
    lines = ["# schema = sample-manifest.v2"]
    if n_ranked:
        # 注入順が raw の取得時刻の順であることをシートに残す（読み手が試料の表の出所に出す。一部が空の
        # シートでは、読み手は空欄を「取得時刻が読めず不明」として扱い、MS-DIAL の番号で埋めない）
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
        "unreadable": unreadable, "unreadable_reasons": reasons, "groups": dict(sorted(groups.items())),
        "qc_samples": [e["sample"] for e in entries if e["role"] == "qc"],
        "blank_samples": [e["sample"] for e in entries if e["role"] == "blank"],
        "unassigned_samples": [e["sample"] for e in bio if e["group"] == UNASSIGNED],
        "files_not_in_assay": files_not_in_assay, "ambiguous_files": ambiguous_files,
        "conflicting_samples": conflicting, "sample_table_format": fmt,
        "factor_columns": [header[i] for i in factor_idx], "renamed_factor_columns": renamed,
    }
