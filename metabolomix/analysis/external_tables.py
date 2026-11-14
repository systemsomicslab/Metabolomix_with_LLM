"""論文側の値の表（MetaboLights の MAF・Metabolomics Workbench の mwTab・補足の xlsx）を
特徴 × 試料の行列として読む（analysis-spec.v1 の外部の入力。reanalysis-study の
`paper_values.py` / `xlsx.py` を移したもの）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。標準ライブラリと numpy だけを使う
（xlsx は zip の中の XML を読む。式は計算せず、キャッシュされた値を読む）。

- **MAF**（`m_*.tsv`）: 見出しの `smallmolecule_abundance_std_error_sub` より右の列が試料。特徴の ID は
  `row_<n>:<metabolite_identification>`（n は見出しの次の行を 1 とするデータの行の番号。空の行も数える）。
- **mwTab**（テキスト形式）: `MS_METABOLITE_DATA_START` 〜 `MS_METABOLITE_DATA_END` の区画。1 行目が
  `Samples`、`Factors` の行は読まない。特徴の ID は 1 列目の名前。
- **xlsx**: 1 枚のシートの横長の表。見出しの行は `feature_column` と `first_sample_column` の両方に値がある
  最初の行で、それより上の表題・注記は読まない。試料は見出しの行の `first_sample_column` から右の空でない
  セル。試料のセルに数が 1 つも無く、欠測の印（`ND`・`<LOQ`・`n.d.`・`-` など）でない文字がある行（群の
  見出しなど）は特徴にしない。飛ばした行の数と名前の例を返す（`n_rows_skipped_as_text`）。欠測の印だけの
  行（どの試料でも検出されなかった物質）は特徴として読み、そのセルを読めないセルに数える。

数として読めないセル（空・`NA`・`ND`・`<LOQ` など）と非有限の値は欠測（NaN）にし、空でない読めない
セルの数と例を返す（黙って捨てない）。試料・特徴の重複は `ANALYSIS_INPUT_INVALID`。
"""
from __future__ import annotations

import csv
import math
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = ["MAF_FIXED_LAST", "SOURCES", "ExternalTable", "column_index", "read_maf", "read_mwtab",
           "read_table", "read_xlsx_matrix", "read_xlsx_rows"]

#: 読める形式（analysis-spec.v1 の `input.source` の値）。
SOURCES = ("maf", "mwtab", "xlsx")
#: MAF の固定の列の最後（これより右が試料の列）。
MAF_FIXED_LAST = "smallmolecule_abundance_std_error_sub"
_MISSING = frozenset({"", "na", "nan", "null", "none"})
#: 論文の表の欠測の印（読めないセルとして数える）。xlsx の行が群の見出しかを決めるとき、これらは文字と
#: 見なさない。`<`・`>`・`≤`・`≥` で始まるもの（`<LOQ`・`<0.1`）と `#` で始まる Excel のエラー値も印。
_MISSING_MARKERS = frozenset({
    "nd", "n.d.", "n.d", "nd.", "n/d", "nq", "n.q.", "n.q", "bdl", "bql", "blq", "lod", "loq",
    "tr", "tr.", "trace", "-", "--", "---", "–", "—", "n/a", "not detected", "below lod",
    "below loq"})
_INVALID = "ANALYSIS_INPUT_INVALID"
_REL_OFFICE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
#: 読めないセルの例として返す種類の数。
_EXAMPLES = 5


@dataclass
class ExternalTable:
    """特徴 × 試料の値（NaN は欠測）と、読めなかったセルの数・例。"""
    features: list[str]
    samples: list[str]
    values: np.ndarray                  # shape (n_features, n_samples)
    n_unreadable: int = 0
    unreadable_examples: list[str] = field(default_factory=list)
    #: 文字の行として飛ばした行の特徴の名前（xlsx だけ。ほかの形式は None）。
    skipped_rows: list[str] | None = None

    def info(self) -> dict:
        out = {"n_features": len(self.features), "n_samples": len(self.samples),
               "n_unreadable_cells": self.n_unreadable,
               "unreadable_examples": list(self.unreadable_examples)}
        if self.skipped_rows is not None:
            out["n_rows_skipped_as_text"] = len(self.skipped_rows)
            out["rows_skipped_as_text"] = self.skipped_rows[:_EXAMPLES]
        return out


class _Cells:
    """セルの文字列を数にし、空でない読めないセルを数える。"""

    def __init__(self):
        self.n_unreadable = 0
        self.examples: list[str] = []

    def number(self, text) -> float:
        cell = (text or "").strip()
        if cell.lower() in _MISSING:
            return math.nan
        try:
            value = float(cell)
        except ValueError:
            value = math.nan
        if math.isfinite(value):
            return value
        self.n_unreadable += 1
        if cell not in self.examples and len(self.examples) < _EXAMPLES:
            self.examples.append(cell)
        return math.nan


def _is_missing_marker(text: str) -> bool:
    """空・`NA` などの欠測の綴りか、論文の表の欠測の印（`_MISSING_MARKERS`・`<`/`>` 始まり・Excel のエラー値）。"""
    cell = text.strip().lower()
    return cell in _MISSING or cell in _MISSING_MARKERS or cell[:1] in ("<", ">", "≤", "≥", "#")


def _fail(label: str, message: str, **details) -> DomainError:
    return DomainError(_INVALID, f"{label}: {message}", {"field": label, **details})


def _require_file(path_text: str, label: str) -> Path:
    path = Path(path_text)
    if not path.is_file():
        raise DomainError("ANALYSIS_INPUT_NOT_FOUND", f"{label} のファイルが見つからない: {path}",
                          {"field": label, "path": str(path)})
    return path


def _unique(names: list[str], what: str, label: str) -> None:
    duplicated = sorted({n for n in names if names.count(n) > 1})
    if duplicated:
        raise _fail(label, f"{what}が重複している: {duplicated[:10]}", duplicated=duplicated[:50])


def _table(features, samples, rows, cells: _Cells, label: str) -> ExternalTable:
    _unique(samples, "試料の名前", label)
    _unique(features, "特徴の ID", label)
    if not samples or not features:
        raise _fail(label, f"試料（{len(samples)}）か特徴（{len(features)}）が 1 つも無い")
    values = np.array(rows, dtype=float).reshape(len(features), len(samples))
    return ExternalTable(list(features), list(samples), values, cells.n_unreadable, cells.examples)


def _tsv_rows(path: Path) -> list[list[str]]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


# ---------- MAF ----------

def read_maf(path_text: str, label: str = "input") -> ExternalTable:
    """MetaboLights の MAF（`m_*.tsv`）。試料は `smallmolecule_abundance_std_error_sub` より右の列。"""
    path = _require_file(path_text, label)
    rows = _tsv_rows(path)
    head = [c.strip() for c in rows[0]] if rows else []
    if MAF_FIXED_LAST not in head:
        raise _fail(label, f"MAF の見出しに {MAF_FIXED_LAST} が無い（MAF ではない）")
    start = head.index(MAF_FIXED_LAST) + 1
    columns = [i for i in range(start, len(head)) if head[i]]
    samples = [head[i] for i in columns]
    name_col = head.index("metabolite_identification") if "metabolite_identification" in head else None
    cells, features, values = _Cells(), [], []
    for n, row in enumerate(rows[1:], 1):
        if not any(c.strip() for c in row):
            continue
        name = row[name_col].strip() if name_col is not None and name_col < len(row) else ""
        features.append(f"row_{n}:{name}")
        values.extend(cells.number(row[i] if i < len(row) else "") for i in columns)
    return _table(features, samples, values, cells, label)


# ---------- mwTab ----------

def read_mwtab(path_text: str, label: str = "input") -> ExternalTable:
    """Metabolomics Workbench の mwTab（テキスト形式）の `MS_METABOLITE_DATA` の区画。"""
    path = _require_file(path_text, label)
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    stripped = [line.strip() for line in lines]
    try:
        start = stripped.index("MS_METABOLITE_DATA_START")
        end = stripped.index("MS_METABOLITE_DATA_END", start)
    except ValueError as exc:
        raise _fail(label, "MS_METABOLITE_DATA_START 〜 _END の区画が無い（mwTab ではない）") from exc
    block = [line.split("\t") for line in lines[start + 1:end]]
    if not block or block[0][0].strip() != "Samples":
        raise _fail(label, "MS_METABOLITE_DATA の区画の 1 行目が Samples でない")
    samples = [s.strip() for s in block[0][1:]]
    while samples and not samples[-1]:
        samples.pop()
    cells, features, values = _Cells(), [], []
    for row in block[1:]:
        name = row[0].strip() if row else ""
        if not name or name == "Factors":
            continue
        features.append(name)
        values.extend(cells.number(row[1 + k] if 1 + k < len(row) else "")
                      for k in range(len(samples)))
    return _table(features, samples, values, cells, label)


# ---------- xlsx ----------

def column_index(letters: str, label: str = "input") -> int:
    """列の文字（A・Z・AA）を 0 始まりの番号にする。"""
    text = str(letters or "").strip().upper()
    if not text or not re.fullmatch(r"[A-Z]+", text):
        raise _fail(label, f"列の文字として読めない: {letters!r}")
    n = 0
    for ch in text:
        n = n * 26 + ord(ch) - 64
    return n - 1


def _text(element: ET.Element) -> str:
    """`<si>` や `<is>` の文字列。`<t>` 直下と `<r><t>` を連結し、ふりがな（`<rPh>`）は読まない。"""
    parts = []
    for child in element:
        tag = child.tag.rpartition("}")[2]
        if tag == "t":
            parts.append(child.text or "")
        elif tag == "r":
            parts.extend(t.text or "" for t in child if t.tag.rpartition("}")[2] == "t")
    return "".join(parts)


def _read_xml(archive: zipfile.ZipFile, name: str, label: str) -> ET.Element:
    try:
        return ET.fromstring(archive.read(name))
    except KeyError as exc:
        raise _fail(label, f"xlsx の中に {name} が無い") from exc
    except ET.ParseError as exc:
        raise _fail(label, f"xlsx の {name} が XML として読めない（{exc}）") from exc


def _sheets(archive: zipfile.ZipFile, label: str) -> list[tuple[str, str]]:
    """[(シート名, zip の中のパス)]。workbook の並び順。"""
    workbook = _read_xml(archive, "xl/workbook.xml", label)
    rels = _read_xml(archive, "xl/_rels/workbook.xml.rels", label)
    target = {r.get("Id"): r.get("Target", "") for r in rels}
    out = []
    for sheet in workbook.iterfind("{*}sheets/{*}sheet"):
        part = target.get(sheet.get(f"{{{_REL_OFFICE}}}id"))
        if part is None:
            continue
        out.append((sheet.get("name", ""), part.lstrip("/") if part.startswith("/") else f"xl/{part}"))
    return out


def _cell_value(cell: ET.Element, shared: list[str], label: str) -> str:
    kind = cell.get("t", "n")
    if kind == "inlineStr":
        inline = cell.find("{*}is")
        return _text(inline) if inline is not None else ""
    value = cell.find("{*}v")
    raw = (value.text or "") if value is not None else ""
    if kind == "s":
        try:
            return shared[int(raw)]
        except (ValueError, IndexError) as exc:
            raise _fail(label, f"セル {cell.get('r')} の共有文字列の番号 {raw!r} が読めない") from exc
    if kind == "b":
        return "TRUE" if raw.strip() == "1" else "FALSE"
    return raw  # n（数）・str（式の文字列の結果）・e（エラー値）・d（日付の文字列）はそのまま


def read_xlsx_rows(path_text: str, sheet: str | None = None,
                   label: str = "input") -> list[list[str]]:
    """シート（名前を省くと先頭）の行を文字列の行の並びにする。空のセルは ""、空の行は []。

    セルの参照の飛びは空の文字列で、行の飛びは空の行で埋める（行の番号と並びの位置がずれない）。
    """
    path = _require_file(path_text, label)
    try:
        archive = zipfile.ZipFile(path)
    except (zipfile.BadZipFile, OSError) as exc:
        raise _fail(label, f"xlsx として開けない（{exc}）") from exc
    with archive:
        sheets = _sheets(archive, label)
        if not sheets:
            raise _fail(label, "xlsx にシートが無い")
        if sheet is None:
            part = sheets[0][1]
        else:
            part = next((p for name, p in sheets if name == sheet), None)
            if part is None:
                raise _fail(label, f"シート {sheet!r} が無い（あるシート: {[n for n, _ in sheets]}）",
                            sheets=[n for n, _ in sheets])
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = [_text(si) for si in
                      _read_xml(archive, "xl/sharedStrings.xml", label).iterfind("{*}si")]
        root = _read_xml(archive, part, label)

    rows: list[list[str]] = []
    next_row = 1
    for row in root.iterfind("{*}sheetData/{*}row"):
        number = int(row.get("r") or next_row)
        rows.extend([] for _ in range(number - next_row))
        cells: list[str] = []
        col = 0
        for cell in row.iterfind("{*}c"):
            ref = re.match(r"([A-Za-z]+)", cell.get("r") or "")
            if ref:
                col = column_index(ref.group(1), label)
            cells.extend([""] * (col + 1 - len(cells)))
            cells[col] = _cell_value(cell, shared, label)
            col += 1
        rows.append(cells)
        next_row = number + 1
    return rows


def read_xlsx_matrix(path_text: str, *, sheet: str | None = None, feature_column: str = "A",
                     first_sample_column: str = "B", label: str = "input") -> ExternalTable:
    """xlsx の 1 枚のシート（特徴 × 試料の横長の表）を行列にする。列は Excel の文字（A・B・AA）で指す。"""
    fcol = column_index(feature_column, label)
    scol = column_index(first_sample_column, label)
    if scol <= fcol:
        raise _fail(label, "first_sample_column は feature_column より右にする（"
                           f"{first_sample_column} と {feature_column}）")
    rows = read_xlsx_rows(path_text, sheet, label)

    def cell(row, i):
        return row[i].strip() if i < len(row) else ""

    head_at = next((n for n, r in enumerate(rows) if cell(r, fcol) and cell(r, scol)), None)
    if head_at is None:
        raise _fail(label, f"列 {feature_column} と {first_sample_column} の両方に値がある見出しの行が無い"
                           f"（シート {sheet or '先頭'}）")
    head = rows[head_at]
    columns = [i for i in range(scol, len(head)) if cell(head, i)]
    samples = [cell(head, i) for i in columns]
    cells, features, values, skipped = _Cells(), [], [], []
    for row in rows[head_at + 1:]:
        name = cell(row, fcol)
        if not name:
            continue
        texts = [cell(row, i) for i in columns]
        probe = _Cells()
        numbers = [probe.number(t) for t in texts]
        if all(math.isnan(v) for v in numbers) and not all(_is_missing_marker(t) for t in texts):
            skipped.append(name)  # 群の見出しなどの文字の行（欠測の印だけの行は特徴として読む）
            continue
        features.append(name)
        values.extend(cells.number(t) for t in texts)
    table = _table(features, samples, values, cells, label)
    table.skipped_rows = skipped
    return table


def read_table(spec_input: dict, label: str = "input") -> ExternalTable:
    """正規化した入力（`source` が maf / mwtab / xlsx）を読む。"""
    source = spec_input["source"]
    if source == "maf":
        return read_maf(spec_input["path"], label)
    if source == "mwtab":
        return read_mwtab(spec_input["path"], label)
    if source == "xlsx":
        return read_xlsx_matrix(spec_input["path"], sheet=spec_input["sheet"],
                                feature_column=spec_input["feature_column"],
                                first_sample_column=spec_input["first_sample_column"],
                                label=label)
    raise ValueError(f"unknown source: {source}")
