"""論文側の値の表（MAF・mwTab・xlsx）の読み取り（`external_tables`）。fixture は合成のファイル。"""
from __future__ import annotations

import math
import zipfile

import numpy as np
import pytest

from metabolomix.analysis import external_tables as X
from metabolomix.core.atomic_io import DomainError
from tests.external_table_fixtures import write_maf, write_mwtab, write_xlsx

SAMPLES = ["S1", "S2", "S3"]


def _nan_equal(actual, expected):
    a = np.asarray(actual, dtype=float)
    b = np.asarray(expected, dtype=float)
    return a.shape == b.shape and bool(np.all((a == b) | (np.isnan(a) & np.isnan(b))))


# ---------- MAF ----------

def test_maf_samples_are_right_of_the_fixed_columns_and_ids_keep_the_row_number(tmp_path):
    path = write_maf(tmp_path / "m_study.tsv", SAMPLES,
                     [("Glucose", ["1.5", "2", "3"]), None, ("Lactate", ["ND", "", "<LOQ"]),
                      ("Lactate", ["4", "5", "6"])])
    table = X.read_maf(str(path))
    assert table.samples == SAMPLES
    # 空の行は飛ばすが番号は進む（同じ名前でも行の番号で別の特徴）
    assert table.features == ["row_1:Glucose", "row_3:Lactate", "row_4:Lactate"]
    assert _nan_equal(table.values, [[1.5, 2, 3], [math.nan] * 3, [4, 5, 6]])
    # 空のセルは欠測、空でない読めないセル（ND・<LOQ）は欠測にして数える
    assert table.n_unreadable == 2 and table.unreadable_examples == ["ND", "<LOQ"]
    assert table.info() == {"n_features": 3, "n_samples": 3, "n_unreadable_cells": 2,
                            "unreadable_examples": ["ND", "<LOQ"]}


def test_maf_without_the_fixed_column_is_invalid(tmp_path):
    path = tmp_path / "m_bad.tsv"
    path.write_text("metabolite_identification\tS1\nGlucose\t1\n", encoding="utf-8")
    with pytest.raises(DomainError) as info:
        X.read_maf(str(path))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    assert X.MAF_FIXED_LAST in info.value.message


def test_duplicated_sample_names_are_invalid(tmp_path):
    path = write_maf(tmp_path / "m_dup.tsv", ["S1", "S1"], [("Glucose", ["1", "2"])])
    with pytest.raises(DomainError) as info:
        X.read_maf(str(path))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    assert info.value.details["duplicated"] == ["S1"]


def test_missing_file_is_not_found(tmp_path):
    with pytest.raises(DomainError) as info:
        X.read_maf(str(tmp_path / "nope.tsv"), "correlate.external")
    assert info.value.code == "ANALYSIS_INPUT_NOT_FOUND"
    assert info.value.details["field"] == "correlate.external"


# ---------- mwTab ----------

def test_mwtab_reads_the_data_block_and_skips_factors(tmp_path):
    path = write_mwtab(tmp_path / "ST000000.txt", SAMPLES,
                       [("Glucose", ["1", "2", "3"]), ("Lactate", ["4", "NA", "inf"])])
    table = X.read_mwtab(str(path))
    assert table.samples == SAMPLES and table.features == ["Glucose", "Lactate"]
    assert _nan_equal(table.values, [[1, 2, 3], [4, math.nan, math.nan]])
    # NA は欠測（数えない）、非有限の inf は読めない値として数える
    assert table.n_unreadable == 1 and table.unreadable_examples == ["inf"]


def test_mwtab_duplicated_feature_is_invalid(tmp_path):
    path = write_mwtab(tmp_path / "dup.txt", SAMPLES,
                       [("Glucose", ["1", "2", "3"]), ("Glucose", ["4", "5", "6"])])
    with pytest.raises(DomainError) as info:
        X.read_mwtab(str(path))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    assert info.value.details["duplicated"] == ["Glucose"]


def test_mwtab_without_the_block_is_invalid(tmp_path):
    path = tmp_path / "plain.txt"
    path.write_text("Samples\tS1\nGlucose\t1\n", encoding="utf-8")
    with pytest.raises(DomainError) as info:
        X.read_mwtab(str(path))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"


# ---------- xlsx ----------

def _paper_sheet():
    return [["Table S1. Plasma metabolites"],          # 表題（見出しより上は読まない）
            [],
            ["Metabolite", "S1", "S2", "S3"],          # 見出しの行
            ["Amino acids", "group", None, None],      # 文字だけの行（群の見出し）は特徴にしない
            ["Alanine", 1.5, 2, "ND"],
            ["Glycine", 4, None, 6]]


def test_xlsx_matrix_finds_the_header_and_skips_title_and_group_rows(tmp_path):
    path = write_xlsx(tmp_path / "supp.xlsx", {"Table S1": _paper_sheet()})
    table = X.read_xlsx_matrix(str(path))
    assert table.samples == SAMPLES
    assert table.features == ["Alanine", "Glycine"]
    assert _nan_equal(table.values, [[1.5, 2, math.nan], [4, math.nan, 6]])
    assert table.n_unreadable == 1 and table.unreadable_examples == ["ND"]
    # 飛ばした文字の行は数と名前の例を返す（黙って捨てない）
    assert table.skipped_rows == ["Amino acids"]
    assert table.info() == {"n_features": 2, "n_samples": 3, "n_unreadable_cells": 1,
                            "unreadable_examples": ["ND"], "n_rows_skipped_as_text": 1,
                            "rows_skipped_as_text": ["Amino acids"]}


def test_xlsx_rows_of_missing_markers_are_features_not_headings(tmp_path):
    sheet = [["Metabolite", "S1", "S2", "S3"],
             ["Glucose", 1, 2, 3],
             ["Lactate", "ND", "n.d.", "<LOQ"],          # どの試料でも検出されなかった物質
             ["Pyruvate", "-", "#DIV/0!", "≤0.1"],     # 印・Excel のエラー値・≤ 始まり
             ["Citrate", 4, "ND", 6],
             ["Treated", "group", "ND", None]]          # 印でない文字がある行は群の見出し
    table = X.read_xlsx_matrix(str(write_xlsx(tmp_path / "supp.xlsx", {"T": sheet})))
    assert table.features == ["Glucose", "Lactate", "Pyruvate", "Citrate"]
    assert _nan_equal(table.values[1:3], [[math.nan] * 3] * 2)
    assert table.n_unreadable == 7                       # 3 + 3 + 1（飛ばした行のセルは数えない）
    assert table.unreadable_examples == ["ND", "n.d.", "<LOQ", "-", "#DIV/0!"]
    assert table.info()["n_rows_skipped_as_text"] == 1
    assert table.info()["rows_skipped_as_text"] == ["Treated"]
    # 例は最大 5 つ、数は全件
    many = [["Metabolite", "S1"], *[[f"Group {i}", "label"] for i in range(7)], ["Glucose", 1]]
    info = X.read_xlsx_matrix(str(write_xlsx(tmp_path / "many.xlsx", {"T": many}))).info()
    assert info["n_rows_skipped_as_text"] == 7 and len(info["rows_skipped_as_text"]) == 5


def test_maf_and_mwtab_have_no_skipped_rows_fields(tmp_path):
    mwtab = write_mwtab(tmp_path / "mw.txt", SAMPLES, [("Lactate", ["ND", "ND", "ND"])])
    table = X.read_mwtab(str(mwtab))
    assert table.skipped_rows is None and "n_rows_skipped_as_text" not in table.info()
    assert table.features == ["Lactate"] and table.n_unreadable == 3


def test_xlsx_rows_fill_gaps_and_read_inline_and_shared_strings(tmp_path):
    path = write_xlsx(tmp_path / "supp.xlsx", {"Table S1": _paper_sheet()})
    rows = X.read_xlsx_rows(str(path))
    assert rows[0] == ["Table S1. Plasma metabolites"]   # インラインの文字列
    assert rows[1] == []                                 # 書かれていない行は空の行で埋める
    assert rows[2] == ["Metabolite", "S1", "S2", "S3"]   # 共有文字列
    assert rows[5] == ["Glycine", "4", "", "6"]          # 参照の飛びは空の文字列


def test_xlsx_named_sheet_and_columns(tmp_path):
    other = [["note"], ["", "ID", "x", "A", "B"], ["", "M1", "skip", 1, 2], ["", "M2", "skip", 3, 4]]
    path = write_xlsx(tmp_path / "two.xlsx", {"Table S1": _paper_sheet(), "Table S2": other})
    table = X.read_xlsx_matrix(str(path), sheet="Table S2", feature_column="b",
                               first_sample_column="D")
    assert table.samples == ["A", "B"] and table.features == ["M1", "M2"]
    assert _nan_equal(table.values, [[1, 2], [3, 4]])
    with pytest.raises(DomainError) as info:
        X.read_xlsx_matrix(str(path), sheet="Table S9")
    assert info.value.details["sheets"] == ["Table S1", "Table S2"]


def test_xlsx_errors(tmp_path):
    path = write_xlsx(tmp_path / "supp.xlsx", {"Table S1": _paper_sheet()})
    with pytest.raises(DomainError) as info:
        X.read_xlsx_matrix(str(path), feature_column="C", first_sample_column="B")
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    with pytest.raises(DomainError) as info:
        X.read_xlsx_matrix(str(path), feature_column="F", first_sample_column="G")
    assert "見出しの行が無い" in info.value.message
    broken = tmp_path / "broken.xlsx"
    broken.write_bytes(b"not a zip")
    with pytest.raises(DomainError) as info:
        X.read_xlsx_rows(str(broken))
    assert info.value.code == "ANALYSIS_INPUT_INVALID"
    empty = tmp_path / "empty.xlsx"
    with zipfile.ZipFile(empty, "w") as archive:
        archive.writestr("xl/other.xml", "<x/>")
    with pytest.raises(DomainError) as info:
        X.read_xlsx_rows(str(empty))
    assert "workbook.xml" in info.value.message


def test_column_index():
    assert [X.column_index(c) for c in ("A", "b", "Z", "AA", "AZ", "BA")] == [0, 1, 25, 26, 51, 52]
    with pytest.raises(DomainError):
        X.column_index("A1")


def test_read_table_dispatches_by_source(tmp_path):
    maf = write_maf(tmp_path / "m.tsv", SAMPLES, [("Glucose", ["1", "2", "3"])])
    xlsx = write_xlsx(tmp_path / "s.xlsx", {"T": _paper_sheet()})
    assert X.read_table({"source": "maf", "path": str(maf)}).features == ["row_1:Glucose"]
    table = X.read_table({"source": "xlsx", "path": str(xlsx), "sheet": None,
                          "feature_column": "A", "first_sample_column": "B"})
    assert table.features == ["Alanine", "Glycine"]
