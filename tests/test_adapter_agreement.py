from metabolomix.dataset.agreement import compare_states
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.mztab.loading import load_dataset_state
from tests.curation_fixtures import match_result
from tests.msdial_folder_fixtures import write_matching_mztab, write_msdial_folder


def _spec():
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(4)]
    spots = [
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0, 70.0],
         "gap_filled": [False, True, False, False], "msms": [True, False, False, False],
         "matches": [match_result({0: "PC 34:1"})]},
        {"name": "Glucose|Hexose", "mz": 179.0561, "heights": [10.0, 20.0, 30.0, 40.0]},
        {"name": "Unknown", "mz": 500.1, "heights": [1.0, 2.0, 3.0, 4.0]},
    ]
    return {"spots": spots, "samples": samples}


def _states(tmp_path, spec, mztab_spec=None):
    write_msdial_folder(tmp_path, **spec)
    mz = write_matching_mztab(tmp_path / "out.mzTab", mztab_spec or spec)
    return load_msdial_folder(select_alignment(tmp_path)), load_dataset_state(mztab_path=mz)


def _variant(spec, spot_index, **changes):
    spots = [dict(s) for s in spec["spots"]]
    spots[spot_index].update(changes)
    return {"spots": spots, "samples": spec["samples"]}


def test_same_alignment_agrees(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    assert list(folder.feature_candidates) == ["0"] and mztab.feature_candidates["0"]
    assert set(folder.feature_annotations) == {"0", "1"} == set(mztab.feature_annotations)
    report = compare_states(folder, mztab)
    assert report["ok"] is True, report
    assert report["checks"]["values"]["n_compared"] == 12
    assert report["checks"]["detected"]["n_compared"] == 12
    assert report["known_differences"]


def test_value_difference_is_reported(tmp_path):
    spec = _spec()
    other = _variant(spec, 1, heights=[10.0, 22.0, 30.0, 40.0])   # 1 セルだけ 10% 違う
    report = compare_states(*_states(tmp_path, spec, other))
    assert report["ok"] is False and report["checks"]["values"]["n_over"] == 1
    assert report["checks"]["evidence"]["n_diff"] == 0


def test_evidence_level_difference_is_reported(tmp_path):
    spec = _spec()
    other = _variant(spec, 0, msms=[False] * 4)     # mzTab 側は SME を書かない
    report = compare_states(*_states(tmp_path, spec, other))
    assert report["ok"] is False and report["checks"]["evidence"]["n_diff"] == 1
    assert report["checks"]["evidence"]["examples"][0] == {
        "feature_id": "0", "folder": "sme", "mztab": "sml"}
    assert report["checks"]["values"]["n_over"] == 0
