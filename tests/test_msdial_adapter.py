import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.mztab.polarity import dataset_polarity
from tests.curation_fixtures import arf2_spot_raw, write_arf2
from tests.msdial_folder_fixtures import write_msdial_folder


def _samples(n=3, **extra):
    return [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B", **extra}
            for i in range(n)]


def test_values_detection_and_ids(tmp_path):
    write_msdial_folder(tmp_path, spots=[
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0],
         "gap_filled": [False, True, False]},
        {"name": "Unknown", "mz": 500.1, "heights": [10.0, 20.0, 30.0],
         "gap_filled": [False, False, False]}], samples=_samples())
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.source_format == "msdial_folder" and ds.feature_ids == ["0", "1"]
    assert ds.sample_names == ["s0", "s1", "s2"]
    assert ds.feature_matrix[0].tolist() == [100.0, 0.0, 50.0]
    assert ds.detected_mask[0].tolist() == [True, False, True]
    assert "s1" not in ds.per_sample_mz["0"]
    assert ds.quantification_measure == "peak_height"


def test_sample_roles_batches_and_groups_from_mddata(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], samples=[
        {"file_id": 0, "name": "QC1", "class_id": "QC", "file_type": 2, "order": 1, "batch": 1},
        {"file_id": 1, "name": "x", "class_id": "A", "file_type": 0, "order": 2, "batch": 2},
        {"file_id": 2, "name": "std", "class_id": "S", "file_type": 1, "order": 3, "batch": 2},
        {"file_id": 3, "name": "blk", "class_id": "B", "file_type": 3, "order": 4, "batch": 2}])
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.roles == {"QC1": "qc", "x": "sample", "std": "standard", "blk": "blank"}
    rows = {r["sample_id"]: r for r in ds.sample_metadata_rows}
    assert rows["x"]["group"] == "A" and rows["x"]["batch"] == "2" and rows["x"]["injection_order"] == 2
    assert rows["x"]["provenance"]["group"]["source"] == "msdial_mddata"


def test_without_mddata_falls_back_to_name_roles(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], with_mddata=False, samples=[
        {"file_id": 0, "name": "QC_1", "class_id": "x"},
        {"file_id": 1, "name": "blank_1", "class_id": "x"},
        {"file_id": 2, "name": "plain", "class_id": "x"}])
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert any("mddata" in note for note in ds.load_notes)
    assert ds.msdial_samples is None
    assert ds.roles == {"QC_1": "qc", "blank_1": "blank", "plain": "sample"}


def test_mismatched_arf_and_arf2_refused(tmp_path):          # Review Focus 1
    paths = write_msdial_folder(tmp_path, spots=[{"name": "a"}, {"name": "b", "mz": 600.0}],
                                samples=_samples())
    write_arf2(paths["arf2"], [arf2_spot_raw(spot_id=0)])
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ALIGNMENT_MISMATCH"


def test_renamed_sample_joined_by_file_id(tmp_path):         # Review Focus 2
    write_msdial_folder(tmp_path, spots=[{"name": "a"}], samples=[
        {"file_id": 0, "name": "s0", "class_id": "A"},
        {"file_id": 1, "name": "s1", "class_id": "B", "arf_name": "old_b"},
        {"file_id": 2, "name": "s2", "class_id": "B"}])   # .arf の群は 3 行未満だと読まれない
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.sample_names[1] == "s1"
    assert ds.msdial_samples[1]["arf_file_name"] == "old_b"
    assert {r["sample_id"]: r["group"] for r in ds.sample_metadata_rows}["s1"] == "B"


def test_area_measure_and_bad_measure(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "a", "heights": [100.0, 200.0, 300.0]}],
                        samples=_samples())
    selection = select_alignment(tmp_path)
    ds = load_msdial_folder(selection, measure="peak_area_above_zero")
    assert ds.feature_matrix[0].tolist() == [200.0, 400.0, 600.0]   # arf_row は Key20 = height * 2
    assert ds.quantification_measure == "peak_area_above_zero"
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(selection, measure="foo")
    assert exc.value.code == "DATASET_BAD_REQUEST"


def test_polarity_from_arf2_ion_mode(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "a", "ion_mode": 1},
                                         {"name": "b", "mz": 600.0, "ion_mode": 1}],
                        samples=_samples())
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert dataset_polarity(ds) == ("negative", "arf2_ion_mode", [])


def test_msms_assigned_and_internal_standard_sentinel(tmp_path):
    write_msdial_folder(tmp_path, spots=[
        {"name": "a", "msms": [False, True, False]}, {"name": "b", "mz": 600.0}],
        samples=_samples())
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_metadata["0"]["msms_assigned"] is True
    assert ds.feature_metadata["1"]["msms_assigned"] is False
    assert ds.feature_metadata["0"]["internal_standard_id"] is None
