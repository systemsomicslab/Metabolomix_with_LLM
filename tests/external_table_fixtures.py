"""論文側の値の表（MAF・mwTab・xlsx）の合成ファイルを書く（`external_tables` のテスト用）。

実在の論文の表は使わない。xlsx は zip に最小の XML（workbook・rels・シート・共有文字列）を入れて作る。
"""
from __future__ import annotations

import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

#: MAF の固定の列（試料はこれより右）。
MAF_FIXED = ("database_identifier", "chemical_formula", "smiles", "inchi",
             "metabolite_identification", "mass_to_charge", "retention_time",
             "smallmolecule_abundance_sub", "smallmolecule_abundance_std_dev_sub",
             "smallmolecule_abundance_std_error_sub")


def write_maf(path: Path, samples, rows) -> Path:
    """rows は (名前, [試料ごとの値の文字列]) の並び。None は空の行（行の番号は進む）。"""
    path = Path(path)
    lines = ["\t".join((*MAF_FIXED, *samples))]
    for row in rows:
        if row is None:
            lines.append("\t" * (len(MAF_FIXED) + len(samples) - 1))
            continue
        name, values = row
        fixed = ["" for _ in MAF_FIXED]
        fixed[MAF_FIXED.index("metabolite_identification")] = name
        lines.append("\t".join((*fixed, *map(str, values))))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_mwtab(path: Path, samples, rows, *, factors=True) -> Path:
    """rows は (名前, [試料ごとの値の文字列]) の並び。"""
    path = Path(path)
    lines = ["#METABOLOMICS WORKBENCH STUDY_ID:ST000000 ANALYSIS_ID:AN000000",
             "MS_METABOLITE_DATA:UNITS\tpeak area",
             "MS_METABOLITE_DATA_START",
             "\t".join(("Samples", *samples))]
    if factors:
        lines.append("\t".join(("Factors", *(f"Group:{s[:2]}" for s in samples))))
    lines += ["\t".join((name, *map(str, values))) for name, values in rows]
    lines += ["MS_METABOLITE_DATA_END", "#END"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _letters(index: int) -> str:
    out = ""
    index += 1
    while index:
        index, rest = divmod(index - 1, 26)
        out = chr(65 + rest) + out
    return out


def write_xlsx(path: Path, sheets: dict[str, list[list]]) -> Path:
    """sheets は シート名 → 行（セルの並び。None は空のセル・[] は空の行）。

    文字列は共有文字列（最初の 1 つだけインラインの文字列）、数は数のセルで書く。空のセルと
    空の行は書かない（参照の飛びとして読む側が埋める）。
    """
    path = Path(path)
    shared: list[str] = []
    inline_used = False
    parts: dict[str, str] = {}
    sheet_entries, rel_entries = [], []
    for k, (name, rows) in enumerate(sheets.items(), start=1):
        xml_rows = []
        for r, row in enumerate(rows, start=1):
            cells = []
            for c, value in enumerate(row):
                if value is None or value == "":
                    continue
                ref = f"{_letters(c)}{r}"
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    cells.append(f'<c r="{ref}"><v>{value}</v></c>')
                elif not inline_used:
                    inline_used = True
                    cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>')
                else:
                    if value not in shared:
                        shared.append(value)
                    cells.append(f'<c r="{ref}" t="s"><v>{shared.index(value)}</v></c>')
            if cells:
                xml_rows.append(f'<row r="{r}">{"".join(cells)}</row>')
        parts[f"xl/worksheets/sheet{k}.xml"] = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>')
        sheet_entries.append(f'<sheet name="{escape(name)}" sheetId="{k}" r:id="rId{k}"/>')
        rel_entries.append(
            f'<Relationship Id="rId{k}" Type="http://schemas.openxmlformats.org/officeDocument/'
            f'2006/relationships/worksheet" Target="worksheets/sheet{k}.xml"/>')
    parts["xl/workbook.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets>{"".join(sheet_entries)}</sheets></workbook>')
    parts["xl/_rels/workbook.xml.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'{"".join(rel_entries)}</Relationships>')
    parts["xl/sharedStrings.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        + "".join(f"<si><t>{escape(s)}</t></si>" for s in shared) + "</sst>")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, text in parts.items():
            archive.writestr(name, text)
    return path
