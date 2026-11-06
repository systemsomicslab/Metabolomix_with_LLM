"""`.mddata` が読めない・解けないときと、`.mddata` に無い FileID の扱い、`.mddata` の後出し
（B-1 最終レビュー R38・Task 4・指摘 6）。fixture は tmp に作る。
"""
from __future__ import annotations

import pytest

from metabolomix.core import handles, session_state
from metabolomix.dataset import msdial_adapter
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.dataset.opening import open_dataset
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata, write_msdial_folder


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _samples(n=3):
    return [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(n)]


def _folder(path, **kwargs):
    return write_msdial_folder(path, spots=[{"name": "PC 34:1"}], samples=_samples(), **kwargs)


# ---------- R38: 読めない .mddata は止めずに注記して、無いときと同じに続ける ----------

@pytest.mark.parametrize("payload", [b"not msgpack at all \xff\xfe", b"", b"\x93\x01\x02"])
def test_a_corrupt_mddata_is_a_load_note_not_an_exception(tmp_path, payload):
    paths = _folder(tmp_path)
    paths["mddata"].write_bytes(payload)
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert any(".mddata" in note and "読めなかった" in note for note in ds.load_notes)
    assert ds.msdial_samples is None
    # 宣言ではなく推定: 群は付かず、役割の出所も mddata ではない。
    rows = {r["sample_id"]: r for r in ds.sample_metadata_rows}
    assert rows["s0"]["group"] is None
    assert rows["s0"]["provenance"]["role"]["source"] != "msdial_mddata"
    assert rows["s0"]["provenance"]["role"]["confidence"] == "inferred"


def test_a_mddata_that_cannot_be_resolved_is_a_load_note(tmp_path, monkeypatch):
    """`resolve_mddata_path` が読み取りの失敗で止まっても、選択は成り立つ（注記を残す）。"""
    _folder(tmp_path)

    def broken(_path, _explicit=None):
        raise ValueError("Invalid MS-DIAL project: broken")

    monkeypatch.setattr("metabolomix.dataset.alignment_files.resolve_mddata_path", broken)
    selection = select_alignment(tmp_path)
    assert selection.mddata is None and "broken" in selection.mddata_error
    ds = load_msdial_folder(selection)
    assert any("broken" in note for note in ds.load_notes)


def test_a_programming_error_in_the_mddata_reader_is_not_swallowed(tmp_path, monkeypatch):
    _folder(tmp_path)

    def boom(_path):
        raise TypeError("a bug, not a bad file")

    monkeypatch.setattr(msdial_adapter, "parse_analysis_file_classes", boom)
    with pytest.raises(TypeError):
        load_msdial_folder(select_alignment(tmp_path))


# ---------- Task 4: .mddata に無い FileID は宣言ではなく推定 ----------

def test_a_file_id_missing_from_mddata_is_inferred_not_declared(tmp_path):
    paths = _folder(tmp_path)
    write_mddata(paths["mddata"], [analysis_file_row(0, "s0", "A"), analysis_file_row(1, "s1", "A")])
    ds = load_msdial_folder(select_alignment(tmp_path))
    rows = {r["sample_id"]: r for r in ds.sample_metadata_rows}
    assert rows["s0"]["provenance"]["role"]["source"] == "msdial_mddata"
    assert rows["s2"]["role"] == "sample"
    assert rows["s2"]["provenance"]["role"] == {
        "value": "sample", "source": "default", "confidence": "inferred"}
    assert any("FileID 2" in note for note in ds.load_notes)


# ---------- 指摘 6: .mddata の後出しで基のキャッシュが無効になる ----------

def test_the_base_cache_reloads_when_a_mddata_appears_later(tmp_path):
    folder = tmp_path / "gui"
    _folder(folder, with_mddata=False)
    first = open_dataset(directory=str(folder))
    assert first.dataset.msdial_samples is None
    assert open_dataset(directory=str(folder)).reused_base is True
    write_mddata(folder / "AlignmentResult_2026_02_01_00_00_00.mddata",
                 [analysis_file_row(i, f"s{i}", "A" if i < 2 else "B") for i in range(3)])
    reopened = open_dataset(directory=str(folder))
    assert reopened.reused_base is False
    assert reopened.dataset.msdial_samples is not None
