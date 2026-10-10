"""研究室 MSP の NAME → INCHIKEY 表（spec 2026-10-10 msp-inchikey-design §3.1）。

研究室ライブラリには触れない。MSP は tmp_path に数エントリだけ書く。
"""
from __future__ import annotations

import os

import pytest

from metabolomix.mztab import msp_inchikey
from metabolomix.mztab.msp_inchikey import (
    MSP_DB_PREFIX, MspInchikeyTable, load_msp_inchikey_table, lookup)

CIT = "RHGKLRLOHDJJDR-BYPYZUCNSA-N"   # L-citrulline
ORN = "AHLPHDHHMVZTML-BYPYZUCNSA-N"   # L-ornithine
SUC = "KDYFGRWQOYBRFD-UHFFFAOYSA-N"   # succinic acid
OTHER = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"


def _entry(name, inchikey=None):
    lines = [f"NAME: {name}", "PRECURSORMZ: 176.1030", "PRECURSORTYPE: [M+H]+"]
    if inchikey is not None:
        lines.append(f"INCHIKEY: {inchikey}")
    lines += ["Num Peaks: 1", "70.065\t1000", ""]
    return "\n".join(lines) + "\n"


def _write(tmp_path, text, name="lib.msp", newline="\n"):
    p = tmp_path / name
    p.write_bytes(text.replace("\n", newline).encode("utf-8"))
    return p


@pytest.fixture(autouse=True)
def _fresh_memo(monkeypatch):
    monkeypatch.setattr(msp_inchikey, "_memo", {})


def test_names_map_to_their_inchikey(tmp_path):
    p = _write(tmp_path, _entry("Citrulline", CIT) + _entry("Ornithine", ORN))
    table = load_msp_inchikey_table(p)
    assert table.by_name == {"Citrulline": CIT, "Ornithine": ORN}
    assert table.conflicting == frozenset()
    assert table.file_name == "lib.msp"


def test_empty_and_malformed_inchikeys_are_skipped(tmp_path):
    p = _write(tmp_path, _entry("Empty", "") + _entry("Short", "ABC-DEF")
               + _entry("Lower", CIT.lower()) + _entry("Citrulline", CIT))
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT}


def test_same_name_same_key_is_kept_once(tmp_path):
    p = _write(tmp_path, _entry("Citrulline", CIT) + _entry("Citrulline", CIT))
    table = load_msp_inchikey_table(p)
    assert table.by_name == {"Citrulline": CIT}
    assert table.conflicting == frozenset()


def test_same_name_different_keys_is_conflicting(tmp_path):
    p = _write(tmp_path, _entry("Leucine", CIT) + _entry("Leucine", OTHER)
               + _entry("Leucine", CIT))
    table = load_msp_inchikey_table(p)
    assert "Leucine" not in table.by_name
    assert table.conflicting == frozenset({"Leucine"})


def test_field_names_are_case_insensitive(tmp_path):
    p = _write(tmp_path, f"Name: Citrulline\nInChIKey: {CIT}\nNum Peaks: 0\n\n")
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT}


def test_crlf_and_bom_are_tolerated(tmp_path):
    text = "﻿" + _entry("Citrulline", CIT) + _entry("Ornithine", ORN)
    p = _write(tmp_path, text, newline="\r\n")
    assert load_msp_inchikey_table(p).by_name == {"Citrulline": CIT, "Ornithine": ORN}


def test_entry_without_inchikey_does_not_borrow_the_next_key(tmp_path):
    # 2 件目は NAME を持たない壊れたエントリ。1 件目の名前を引き継がせない。
    text = _entry("NoKey") + f"PRECURSORMZ: 1\nINCHIKEY: {CIT}\n\n"
    assert load_msp_inchikey_table(_write(tmp_path, text)).by_name == {}


def test_lookup_resolves_only_mspdb_identifiers(tmp_path):
    table = MspInchikeyTable("lib.msp", {"Citrulline": CIT}, frozenset({"Leucine"}))
    assert lookup(table, "MspDB:Citrulline") == (CIT, "resolved")
    assert lookup(table, "MspDB:Leucine") == (None, "conflicting")
    assert lookup(table, "MspDB:Unknown") == (None, "unmatched")
    assert lookup(table, "TextDB:Citrulline") == (None, "not_msp")
    assert lookup(table, "Citrulline") == (None, "not_msp")
    assert lookup(table, None) == (None, "not_msp")
    assert lookup(table, "") == (None, "not_msp")
    assert MSP_DB_PREFIX == "MspDB:"


def test_lookup_keeps_colons_in_the_name(tmp_path):
    table = MspInchikeyTable("lib.msp", {"LPC 18:0": CIT, "Citrulline": ORN}, frozenset())
    assert lookup(table, "MspDB:LPC 18:0") == (CIT, "resolved")
    assert lookup(table, "  MspDB:Citrulline  ") == (ORN, "resolved")


def test_table_is_memoized_until_the_file_changes(tmp_path, monkeypatch):
    p = _write(tmp_path, _entry("Citrulline", CIT))
    calls = []
    real_scan = msp_inchikey._scan
    monkeypatch.setattr(msp_inchikey, "_scan",
                        lambda path: calls.append(path) or real_scan(path))
    first = load_msp_inchikey_table(p)
    assert load_msp_inchikey_table(p) is first
    assert len(calls) == 1
    p.write_text(_entry("Ornithine", ORN), encoding="utf-8")
    st = p.stat()
    os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    assert load_msp_inchikey_table(p).by_name == {"Ornithine": ORN}
    assert len(calls) == 2


def test_tables_for_two_files_are_memoized_independently(tmp_path, monkeypatch):
    pos = _write(tmp_path, _entry("Citrulline", CIT), name="pos.msp")
    neg = _write(tmp_path, _entry("Succinate", SUC), name="neg.msp")
    calls = []
    real_scan = msp_inchikey._scan
    monkeypatch.setattr(msp_inchikey, "_scan",
                        lambda path: calls.append(path) or real_scan(path))
    for _ in range(2):
        assert load_msp_inchikey_table(pos).by_name == {"Citrulline": CIT}
        assert load_msp_inchikey_table(neg).by_name == {"Succinate": SUC}
    assert len(calls) == 2


def test_missing_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        load_msp_inchikey_table(tmp_path / "absent.msp")


def test_module_imports_no_mcp_or_library_store():
    import ast
    from pathlib import Path
    tree = ast.parse(Path(msp_inchikey.__file__).read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not any(m.startswith(("metabolomix.core.mcp_core", "metabolomix.core.session_state",
                                 "metabolomix.library", "metabolomix.tools"))
                   for m in imported), imported
