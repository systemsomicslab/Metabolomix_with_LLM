"""`.msp`（テキスト）参照ライブラリを読む。

別名表は上流 `MspFileParcer.cs`（`MsdialWorkbench` リポジトリ、
`src/Common/CommonStandard/Parser/MspFileParcer.cs` の `SetMspField` 内
`switch (fieldName.ToLower())`）の分岐を写している。`case "DB#":` は
switch のキーが小文字化された後に比較されるため到達不能な死んだ枝なので
写していない。

ファイルは 1 行ずつ読む（研究室の参照ライブラリは 0.5〜1.2 GB あり、全体を
メモリへ載せると初回の store 構築が落ちる）。文字コードは行ごとに
UTF-8 → cp932 → latin-1 の順で試す。

`.dbs`（`metabolomix/library/dbs.py`）と同じ正規化レコード（`record.make_record`）
を返す。`spectrum` は `[[mz, intensity], ...]`、`ion_mode` は
`"positive"` / `"negative"` の文字列、値が無いフィールドは `None`。

deps: record のみ。mcp_core / session_state / tools_* / dbs を import しない。
"""
from __future__ import annotations

import codecs
from pathlib import Path
from typing import Iterator

from metabolomix.library.record import make_record

#: 小文字化したフィールド名 → 正規化レコードのキー。
#: `None` は「上流では意味を持つが正規化レコードには対応する枠が無いので読み捨てる」
#: （`comment` / `comments` が該当）。ここに無いフィールド名も同様に読み捨てる。
_SPECTRUM_MARKER = "__spectrum__"

_FIELD_ALIASES: dict[str, str | None] = {
    "precursormz": "precursor_mz",
    "precursor_mz": "precursor_mz",
    "precursor_m/z": "precursor_mz",
    "precursortype": "adduct",
    "precursor_type": "adduct",
    "ionmode": "ion_mode",
    "ion_mode": "ion_mode",
    "formula": "formula",
    "inchikey": "inchikey",
    "inchi_key": "inchikey",
    "inchi key": "inchikey",
    "smiles": "smiles",
    "ontology": "ontology",
    "compoundclass": "compound_class",
    "retentiontime": "rt",
    "retention_time": "rt",
    "rt": "rt",
    "num peaks": _SPECTRUM_MARKER,
    "numpeaks": _SPECTRUM_MARKER,
    "num_peaks": _SPECTRUM_MARKER,
    "comment": None,
    "comments": None,
}

#: `float` へ変換するべき正規化レコードのキー。
_FLOAT_FIELDS = {"precursor_mz", "rt"}


def _parse_field_line(line: str) -> tuple[str | None, str | None]:
    """`フィールド名: 値` を分ける。`:` が無ければ `(None, None)`。"""
    if ":" not in line:
        return None, None
    name, value = line.split(":", 1)
    return name.strip(), value.strip()


def _parse_ion_mode(value: str) -> str | None:
    normalized = value.strip().lower()
    return normalized if normalized in ("positive", "negative") else None


#: 行ごとに試す文字コード。latin-1 は全バイト列を受けるので最後の砦になる
#: （化合物名が化けても、1 行のせいでライブラリ全体が読めなくなるよりよい）。
_FALLBACK_ENCODINGS = ("cp932", "latin-1")
_UTF8_BOM = codecs.BOM_UTF8


def _decode_line(raw: bytes, stats: dict) -> str:
    """1 行を UTF-8 → cp932 → latin-1 の順で解く。UTF-8 以外で解いた行は数える。"""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        pass
    stats["non_utf8_lines"] += 1
    for encoding in _FALLBACK_ENCODINGS:
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise AssertionError("latin-1 は全バイト列を解けるので到達しない")


def _iter_lines(path: str | Path, stats: dict) -> Iterator[str]:
    """ファイルを 1 行ずつ解いて返す。全体をメモリへ載せない。

    研究室の参照ライブラリは 0.5〜1.2 GB ある。以前の
    `read_text().splitlines()` は文字列本体に加えて数千万の行オブジェクトを
    抱え込み、初回の store 構築で MemoryError になりえた。
    """
    with open(path, "rb") as fh:
        first = True
        for raw in fh:
            if first:
                first = False
                if raw.startswith(_UTF8_BOM):
                    raw = raw[len(_UTF8_BOM):]
            yield _decode_line(raw, stats)


def _parse_peak(line: str) -> list[float] | None:
    """`m/z intensity` の 1 行。ピークとして解釈できなければ `None`。"""
    parts = line.split()
    if len(parts) < 2:
        return None
    try:
        return [float(parts[0]), float(parts[1])]
    except ValueError:
        return None


def _finish_record(fields: dict, spectrum: list[list[float]], record_index: int) -> dict:
    spectrum.sort(key=lambda peak: peak[0])
    return make_record(record_index=record_index, library_id=None, spectrum=spectrum, **fields)


def iter_records(path: str | Path, *, stats: dict | None = None) -> Iterator[dict]:
    """`.msp` を 1 行ずつ読んで正規化レコードを順に返す。

    レコードの区切りは `NAME:` で始まる行（大小無視）と空行。`#` で始まる行は
    読み飛ばす。最初の `NAME:` より前の行は捨てる。`Num Peaks:` の宣言本数は
    信用しない（実際がそれより少ないファイルがある）——ピークとして解釈できない
    行に当たった時点でピーク列を閉じ、その行はフィールド行として読み直す。

    `stats` を渡すと `non_utf8_lines`（UTF-8 以外で解いた行数）を書き込む。
    store がこれを meta に残し、`library_load` の note に出す。
    """
    if stats is None:
        stats = {}
    stats["non_utf8_lines"] = 0
    fields: dict[str, object] | None = None
    spectrum: list[list[float]] = []
    in_peaks = False
    record_index = 0

    for line in _iter_lines(path, stats):
        stripped = line.strip()
        if not stripped:
            if fields is not None:
                yield _finish_record(fields, spectrum, record_index)
                record_index += 1
                fields = None
            in_peaks = False
            continue
        if stripped.startswith("#"):
            in_peaks = False
            continue
        field_name, field_value = _parse_field_line(stripped)
        if field_name is not None and field_name.lower() == "name":
            if fields is not None:
                yield _finish_record(fields, spectrum, record_index)
                record_index += 1
            fields = {"name": field_value}
            spectrum = []
            in_peaks = False
            continue
        if fields is None:
            continue
        if in_peaks:
            peak = _parse_peak(stripped)
            if peak is not None:
                spectrum.append(peak)
                continue
            in_peaks = False  # ピーク列の終わり。この行はフィールド行として読み直す。
        if field_name is None:
            continue
        canonical = _FIELD_ALIASES.get(field_name.lower())
        if canonical == _SPECTRUM_MARKER:
            in_peaks = True
            continue
        if canonical is None:
            continue
        if canonical in _FLOAT_FIELDS:
            try:
                fields[canonical] = float(field_value)
            except ValueError:
                fields[canonical] = None
        elif canonical == "ion_mode":
            fields[canonical] = _parse_ion_mode(field_value)
        else:
            fields[canonical] = field_value

    if fields is not None:
        yield _finish_record(fields, spectrum, record_index)
