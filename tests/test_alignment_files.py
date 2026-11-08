import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment


def _touch(folder, *names):
    for name in names:
        (folder / name).write_bytes(b"x")


def test_prefers_peak_properties_and_latest_batch(tmp_path):
    for stem in ("AlignmentResult_2026_01_01_00_00_00", "AlignmentResult_2026_02_01_00_00_00"):
        for suffix in ("_PeakProperties.arf", "_DriftSpots.arf", ".arf2", ".dcl"):
            (tmp_path / f"{stem}{suffix}").write_bytes(b"x")
    sel = select_alignment(tmp_path)
    assert sel.arf.name == "AlignmentResult_2026_02_01_00_00_00_PeakProperties.arf"
    assert sel.arf2.name == "AlignmentResult_2026_02_01_00_00_00.arf2"
    assert sel.dcl.name == "AlignmentResult_2026_02_01_00_00_00.dcl"
    assert "AlignmentResult_2026_01_01_00_00_00_PeakProperties.arf" in sel.not_chosen
    assert "AlignmentResult_2026_02_01_00_00_00_DriftSpots.arf" in sel.not_chosen


def test_arf_without_arf2_is_incomplete(tmp_path):
    _touch(tmp_path, "AlignmentResult_2026_02_01_00_00_00_PeakProperties.arf")
    with pytest.raises(DomainError) as exc:
        select_alignment(tmp_path)
    assert exc.value.code == "ALIGNMENT_INCOMPLETE"


def test_no_arf_is_not_found(tmp_path):
    with pytest.raises(DomainError) as exc:
        select_alignment(tmp_path)
    assert exc.value.code == "ALIGNMENT_NOT_FOUND"


def test_console_layout_alignresult_names(tmp_path):
    for suffix in (".arf", ".arf2", ".dcl"):
        (tmp_path / f"AlignResult-2026_03_01_10_00_00{suffix}").write_bytes(b"x")
    sel = select_alignment(tmp_path / "AlignResult-2026_03_01_10_00_00.arf2")
    assert sel.arf.name == "AlignResult-2026_03_01_10_00_00.arf"


def test_driftsopts_spelling_is_a_drift_file(tmp_path):
    stem = "AlignmentResult_2026_02_01_00_00_00"
    _touch(tmp_path, f"{stem}_DriftSopts.arf", f"{stem}_PeakProperties.arf", f"{stem}.arf2")
    sel = select_alignment(tmp_path)
    assert sel.arf.name == f"{stem}_PeakProperties.arf"
    assert f"{stem}_DriftSopts.arf" in sel.not_chosen
    # DriftSopts だけでも同じ幹の .arf2 に結び付く
    (tmp_path / f"{stem}_PeakProperties.arf").unlink()
    assert select_alignment(tmp_path).arf.name == f"{stem}_DriftSopts.arf"


def test_optional_files_none_and_sidecars(tmp_path):
    stem = "AlignmentResult_2026_02_01_00_00_00"
    _touch(tmp_path, f"{stem}.arf", f"{stem}.arf2")
    sel = select_alignment(tmp_path / f"{stem}.arf")
    assert (sel.dcl, sel.eic, sel.tags, sel.mddata) == (None, None, None, None)
    assert sel.loaded_libraries == ()
    _touch(tmp_path, f"{stem}.EIC.aef", f"{stem}_tags.xml", "x_Loaded.lbm2.dbs")
    sel = select_alignment(tmp_path)
    assert sel.eic.name == f"{stem}.EIC.aef"
    assert sel.tags.name == f"{stem}_tags.xml"
    assert [p.name for p in sel.loaded_libraries] == ["x_Loaded.lbm2.dbs"]


def test_subfolders_are_not_searched(tmp_path):
    sub = tmp_path / "other"
    sub.mkdir()
    _touch(sub, "AlignmentResult_2026_02_01_00_00_00.arf", "AlignmentResult_2026_02_01_00_00_00.arf2")
    with pytest.raises(DomainError) as exc:
        select_alignment(tmp_path)
    assert exc.value.code == "ALIGNMENT_NOT_FOUND"


def test_timestamped_batch_beats_untimestamped(tmp_path):
    """時刻の無い `.arf`（PeakProperties）が、時刻を持つ別バッチ（DriftSopts）に勝たない。"""
    stem = "AlignmentResult_2026_02_01_00_00_00"
    _touch(tmp_path, "legacy_PeakProperties.arf", "legacy.arf2",
           f"{stem}_DriftSopts.arf", f"{stem}.arf2")
    # 時刻の無い側のほうが mtime は新しい。
    import os
    old = (tmp_path / f"{stem}_DriftSopts.arf").stat().st_mtime - 1000
    os.utime(tmp_path / f"{stem}_DriftSopts.arf", (old, old))
    assert select_alignment(tmp_path).arf.name == f"{stem}_DriftSopts.arf"
    from metabolomix.core.path_resolvers import _pick_latest
    assert _pick_latest([str(tmp_path / "legacy_PeakProperties.arf"),
                         str(tmp_path / f"{stem}_DriftSopts.arf")]).endswith("DriftSopts.arf")


def test_arf2_input_prefers_peakproperties(tmp_path):
    stem = "AlignmentResult_2026_02_01_00_00_00"
    _touch(tmp_path, f"{stem}_DriftSopts.arf", f"{stem}_PeakProperties.arf", f"{stem}.arf2")
    assert select_alignment(tmp_path / f"{stem}.arf2").arf.name == f"{stem}_PeakProperties.arf"
