"""1〜2 試料のアラインメントを読む（群はピークの行の形で見分ける。行数では見分けない）。"""
import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.mztab import evidence
from tests.curation_fixtures import arf2_spot_raw, arf_row, write_arf, write_arf2
from tests.msdial_folder_fixtures import write_msdial_folder


def _read(path):
    with open(path, "rb") as handle:
        return arf_reader.deserialize(handle)


@pytest.mark.parametrize("n_samples", [1, 2])
def test_one_and_two_sample_alignments_read(tmp_path, n_samples):
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(n_samples)]
    paths = write_msdial_folder(tmp_path, spots=[{"name": "a"}, {"name": "b", "mz": 600.0}],
                                samples=samples)
    groups = _read(paths["arf"])
    assert len(groups) == 2
    assert all(len(g["AlignedPeakProperties"]) == n_samples for g in groups)
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_matrix.shape == (2, n_samples)
    assert len(ds.sample_names) == n_samples


def test_header_lists_are_not_groups(tmp_path):
    header = [[1.0, 0, 0], [2.0, 0, 0], [3.0, 0, 0]]          # ChromXs の形の入れ子リスト
    real = [[arf_row(file_id=0, mz=500.0, rt=12.0, height=100.0, gap_filled=False)],
            [arf_row(file_id=0, mz=600.0, rt=12.0, height=100.0, gap_filled=False)]]
    path = write_arf(tmp_path / "x.arf", [header, *real])
    groups = _read(path)
    assert [g["MassCenter"] for g in groups] == [500.0, 600.0]
    with open(path, "rb") as handle:
        assert groups == list(arf_reader.iter_groups(handle))
    assert evidence.scan_arf(path, measures=()).spot_mz == [500.0, 600.0]


def test_spot_count_mismatch_wording(tmp_path):
    paths = write_msdial_folder(tmp_path, spots=[{"name": "a"}, {"name": "b", "mz": 600.0}],
                                samples=[{"file_id": i, "name": f"s{i}", "class_id": "A"}
                                         for i in range(3)])
    write_arf2(paths["arf2"], [arf2_spot_raw(spot_id=0)])
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ALIGNMENT_MISMATCH"
    assert exc.value.message == (
        ".arf のスポット数（2）と .arf2 のスポット数（1）が違う。"
        "同じアラインメントの出力でないか、.arf を読み切れていない。")
    assert "別のバッチ" not in exc.value.message


def test_feature_count_mismatch_note_does_not_assert_other_batch():
    import numpy as np
    result = evidence._join([500.0, 600.0], ["s0"], np.ones((2, 1), dtype=bool),
                            np.full((2, 1), np.nan), feature_mz=[500.0], sample_names=["s0"],
                            mz_tolerance=evidence.MZ_TOLERANCE)
    assert result["reason"] == "feature_count_mismatch"
    assert "読み切れていない" in result["note"]
    assert "別のバッチ" not in result["note"]


def test_scan_evidence_missing_row_for_file_id_is_not_detected(tmp_path):
    """あるスポットに FileID 1 の行が無ければ、そのセルは未検出（検出のままにしない）。"""
    rows0 = [arf_row(file_id=i, mz=500.0, rt=12.0, height=100.0, gap_filled=False)
             for i in (0, 1)]
    rows1 = [arf_row(file_id=0, mz=600.0, rt=12.0, height=100.0, gap_filled=False)]
    path = write_arf(tmp_path / "x.arf", [rows0, rows1])
    scan = evidence.scan_arf(path, measures=())
    result = evidence.scan_evidence(scan, feature_mz=[500.0, 600.0], sample_names=["sample_0", "sample_1"])
    assert result["status"] == "ok"
    assert result["detected_mask"].tolist() == [[True, True], [True, False]]
