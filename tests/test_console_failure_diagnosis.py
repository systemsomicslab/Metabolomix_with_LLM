"""開いた形式の切れたファイルの検査（`console.raw_integrity`）と、Console が完了しなかった
ときの見立て（`console.failure_diagnosis`）。B7-T63・B8-T66。

計測ファイルはすべてテストが作る合成の XML（実 raw は使わない）。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metabolomix.console import failure_diagnosis, raw_integrity

_PAD = "<spectrum/>\n" * 4000  # 下限（32 KiB）を超える本体


def _mzml(path: Path, *, closed: bool = True, indexed: bool = True, pad: str = _PAD,
          encoding: str = "utf-8") -> Path:
    head = '<?xml version="1.0"?>\n' + ("<indexedmzML>\n" if indexed else "") + "<mzML>\n"
    tail = "</mzML>\n" + ("</indexedmzML>\n" if indexed else "")
    text = head + pad + (tail if closed else "")
    path.write_bytes(text.encode(encoding))
    return path


# ---------- raw_integrity ----------

@pytest.mark.parametrize("indexed", [True, False])
def test_closed_mzml_passes(tmp_path, indexed):
    assert raw_integrity.open_format_problem(_mzml(tmp_path / "a.mzML", indexed=indexed)) is None


def test_trailing_whitespace_and_utf16_are_accepted(tmp_path):
    path = _mzml(tmp_path / "a.mzML")
    path.write_bytes(path.read_bytes() + b"\r\n  \n")
    assert raw_integrity.open_format_problem(path) is None
    assert raw_integrity.open_format_problem(
        _mzml(tmp_path / "b.mzML", encoding="utf-16-le")) is None


def test_truncated_mzml_is_named(tmp_path):
    problem = raw_integrity.open_format_problem(_mzml(tmp_path / "cut.mzML", closed=False))
    assert problem["name"] == "cut.mzML"
    assert problem["reason"] == "truncated"


def test_header_only_mzml_is_too_small(tmp_path):
    """閉じタグはあるが中身の無い（4 KB の ROSMAP のような）ファイル。"""
    problem = raw_integrity.open_format_problem(_mzml(tmp_path / "tiny.mzML", pad=""))
    assert problem["reason"] == "too_small"
    assert problem["size"] < raw_integrity.MIN_OPEN_FORMAT_BYTES


def test_small_imzml_is_not_held_to_the_size_floor(tmp_path):
    """imzML はメタデータだけで、スペクトルは隣の .ibd にある（画素の少ない正しい imzML は小さい）。"""
    small = _mzml(tmp_path / "x.imzML", pad="", indexed=False)
    assert small.stat().st_size < raw_integrity.MIN_OPEN_FORMAT_BYTES
    assert raw_integrity.open_format_problem(small) is None
    cut = _mzml(tmp_path / "y.imzML", pad="", indexed=False, closed=False)
    assert raw_integrity.open_format_problem(cut)["reason"] == "truncated"


def test_mzxml_needs_its_own_closing_tag(tmp_path):
    good = tmp_path / "a.mzXML"
    good.write_text("<mzXML>\n" + _PAD + "</mzXML>\n", encoding="utf-8")
    bad = tmp_path / "b.mzXML"
    bad.write_text("<mzXML>\n" + _PAD + "</mzML>\n", encoding="utf-8")
    assert raw_integrity.open_format_problem(good) is None
    assert raw_integrity.open_format_problem(bad)["reason"] == "truncated"


def test_vendor_formats_and_folders_are_not_checked(tmp_path):
    wiff = tmp_path / "a.wiff"
    wiff.write_bytes(b"x")
    folder = tmp_path / "b.mzML"
    folder.mkdir()
    assert raw_integrity.broken_open_format_files([wiff, folder]) == []


# ---------- failure_diagnosis ----------

def test_exit_code_1450_is_a_resource_shortage(tmp_path):
    result = failure_diagnosis.diagnose_failure(
        tmp_path / "missing.log", exit_code=1450, raw_paths=[])
    assert result["cause"] == "system_resources"
    assert "1450" in result["message"]


def test_out_of_memory_in_the_log_is_a_resource_shortage(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_text("Loading library\nUnhandled Exception: System.OutOfMemoryException\n",
                   encoding="utf-8")
    result = failure_diagnosis.diagnose_failure(log, exit_code=-532462766, raw_paths=[])
    assert result["cause"] == "system_resources"


def test_broken_staged_input_is_reported_first(tmp_path):
    good = _mzml(tmp_path / "S1.mzML")
    cut = _mzml(tmp_path / "S2.mzML", closed=False)
    result = failure_diagnosis.diagnose_failure(
        tmp_path / "msdial.log", exit_code=1, raw_paths=[good, cut])
    assert result["cause"] == "broken_input"
    assert result["suspect_files"] == [{"name": "S2.mzML", "reason": "truncated"}]


def test_file_named_in_an_error_line_is_suspected(tmp_path):
    raws = [_mzml(tmp_path / f"S{i}.mzML") for i in range(3)]
    log = tmp_path / "msdial.log"
    log.write_text(
        "Start\n"
        f"System.Xml.XmlException: Unexpected end of file in {tmp_path / 'S1.mzML'} line 9\n",
        encoding="utf-8")
    result = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=raws)
    assert result["cause"] == "input_suspected"
    assert result["suspect_files"] == [{"name": "S1.mzML", "reason": "named_in_error"}]
    # 返すログの行は絶対パスをファイル名にしてある（置き場所を戻り値に出さない）。
    [line] = result["log_lines"]
    assert "S1.mzML" in line and str(tmp_path) not in line


def test_last_file_named_before_the_error_is_a_weaker_suspect(tmp_path):
    raws = [_mzml(tmp_path / f"S{i}.mzML") for i in range(3)]
    log = tmp_path / "msdial.log"
    log.write_text("Reading S0.mzML\nReading S2.mzML\nUnhandled exception: index out of range\n",
                   encoding="utf-8")
    result = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=raws)
    assert result["suspect_files"] == [{"name": "S2.mzML", "reason": "last_named_before_error"}]


def test_nothing_found_is_an_empty_diagnosis(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_text("Start\nDone\n", encoding="utf-8")
    assert failure_diagnosis.diagnose_failure(log, exit_code=2, raw_paths=[]) == {}


def test_log_lines_hide_library_locations(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_text(r"Error: cannot open C:\lab\secret\library\lab.msp" + "\n", encoding="utf-8")
    result = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=[])
    [line] = result["log_lines"]
    assert "lab.msp" in line
    assert "secret" not in line


def test_a_name_does_not_match_inside_a_longer_name(tmp_path):
    raws = [_mzml(tmp_path / "S1.mzML"), _mzml(tmp_path / "S11.mzML"),
            _mzml(tmp_path / "AS1.mzML")]
    log = tmp_path / "msdial.log"
    log.write_text("Error while reading S11.mzML\n", encoding="utf-8")
    result = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=raws)
    assert result["suspect_files"] == [{"name": "S11.mzML", "reason": "named_in_error"}]


def test_paths_with_spaces_and_forward_slashes_are_reduced(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_text(r"Exception: C:\Program Files\lab\data\lib.lbm2 and C:/lab/private/x.msp"
                   + "\n", encoding="utf-8")
    [line] = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=[])["log_lines"]
    # 2 つのパスが 1 行にあっても、どちらもファイル名だけが残る（間の文は残す）。
    assert line == "Exception: lib.lbm2 and x.msp"


@pytest.mark.parametrize("path, fragments", [
    (r"C:\lab\Lab Library\x.msp", ["lab", "Library"]),
    (r"C:\Users\me\OneDrive - Secret Univ\Lab Library\x.msp",
     ["Users", "OneDrive", "Secret", "Univ", "Library"]),
    ('"D:\\研究室 データ\\ライブラリ\\x.msp"', ["研究室", "データ", "ライブラリ"]),
    (r"\\nas01\share name\lib\x.msp", ["nas01", "share", "name", "lib"]),
    ("file:///C:/Users/me/My Docs/x.msp", ["file", "Users", "My", "Docs"]),
    ('{"path": "C:\\\\Users\\\\me\\\\Lab Library\\\\x.msp"}', ["Users", "Lab", "Library"]),
])
def test_folder_names_with_spaces_leave_no_fragment(tmp_path, path, fragments):
    """フォルダの名前に空白があっても、置き場所の断片を残さずファイル名だけにする。"""
    log = tmp_path / "msdial.log"
    log.write_text(f"Error: cannot open {path} now\n", encoding="utf-8")
    [line] = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=[])["log_lines"]
    assert "x.msp" in line
    assert line.startswith("Error: cannot open ") and line.endswith(" now")
    for fragment in fragments:
        assert fragment not in line


def test_a_slash_in_the_text_after_a_backslash_path_is_kept(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_text(r"Error reading C:\lab\x.mzML failed (and/or) retry" + "\n", encoding="utf-8")
    [line] = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=[])["log_lines"]
    assert line == "Error reading x.mzML failed (and/or) retry"


def test_cp932_log_is_read(tmp_path):
    log = tmp_path / "msdial.log"
    log.write_bytes("エラー: システム リソースが不足しています。\r\n".encode("cp932"))
    result = failure_diagnosis.diagnose_failure(log, exit_code=1, raw_paths=[])
    assert result["cause"] == "system_resources"
