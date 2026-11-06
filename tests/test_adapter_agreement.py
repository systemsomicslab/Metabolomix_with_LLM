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
    ] + [
        # MS/MS と名前付きの代表を持つが、上流の条件 (b)(c)(e)(f) のどれかで SME にならないスポット。
        {"name": "PC 36:1", "mz": 788.6, "msms": [True] * 4,
         "matches": [match_result({0: "PC 36:1", 26: 64})]},                 # (b) 手動修正
        {"name": "PC 38:4", "mz": 810.6, "msms": [True] * 4,
         "matches": [match_result({0: "PC 38:4", 26: 16})]},                 # (c) Text DB
        {"name": "PC 40:6", "mz": 834.6, "msms": [True] * 4, "blank_filtered": True,
         "matches": [match_result({0: "PC 40:6"})]},                          # (e) ブランク除外
        {"name": "PC 32:0 no MS2", "mz": 734.6, "msms": [True] * 4,
         "matches": [match_result({0: "PC 32:0"})]},                          # (f)
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
    # SME は 0 だけ。(b)(c)(e)(f) のスポット 3〜6 は MS/MS も代表の名前もあるが SME にならない。
    assert list(folder.feature_candidates) == ["0"] and mztab.feature_candidates["0"]
    assert [f for f, c in mztab.feature_candidates.items() if c] == ["0"]
    assert set(folder.feature_annotations) == {"0", "1", "3", "4", "5", "6"}         == set(mztab.feature_annotations)
    report = compare_states(folder, mztab)
    assert report["ok"] is True, report
    assert report["checks"]["values"]["n_compared"] == 28
    assert report["checks"]["detected"]["n_compared"] == 28
    assert report["known_differences"]


def test_exact_half_heights_agree_after_half_to_even_rounding(tmp_path):
    spec = _spec()
    exact = _variant(spec, 1, heights=[2.5, 3.5, 0.5, 1.5])   # 偶数丸めで 2 / 4 / 0 / 2
    folder, mztab = _states(tmp_path, exact)
    assert mztab.feature_matrix[1].tolist() == [2.0, 4.0, 0.0, 2.0]
    assert folder.feature_matrix[1].tolist() == [2.5, 3.5, 0.5, 1.5]     # アダプタは丸めない
    report = compare_states(folder, mztab)
    assert report["ok"] is True and report["checks"]["values"]["n_over"] == 0
    assert report["checks"]["values"]["folder_rounded"] == "half_to_even"
    assert any("偶数丸め" in text for text in report["known_differences"])


def test_difference_beyond_rounding_is_still_counted(tmp_path):
    spec = _spec()
    folder_spec = _variant(spec, 1, heights=[2.5, 3.5, 30.0, 40.0])
    other = _variant(spec, 1, heights=[3.0, 3.5, 30.0, 40.0])   # 3 は 2.5 の丸め（2）と違う
    report = compare_states(*_states(tmp_path, folder_spec, other))
    assert report["ok"] is False and report["checks"]["values"]["n_over"] == 1


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


def test_feature_id_on_one_side_only_is_reported(tmp_path):
    spec = _spec()
    extra = {"spots": spec["spots"] + [{"name": "Unknown", "mz": 900.0}], "samples": spec["samples"]}
    report = compare_states(*_states(tmp_path / "a", spec, extra))
    assert report["ok"] is False
    assert report["checks"]["feature_ids"]["only_mztab"] == ["7"]
    assert report["checks"]["feature_ids"]["only_folder"] == []
    folder, mztab = _states(tmp_path / "b", spec)
    report = compare_states(mztab, folder)      # 引数を逆にして only_folder 側
    assert report["ok"] is True
    mztab_less = _variant(spec, 0)
    mztab_less["spots"] = mztab_less["spots"][:-1]
    folder, mztab = _states(tmp_path / "c", spec, mztab_less)
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["feature_ids"]["only_folder"] == ["6"]


def test_sample_order_difference_is_reported(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    mztab.sample_names = mztab.sample_names[::-1]
    mztab.feature_matrix = mztab.feature_matrix[:, ::-1]
    mztab.detected_mask = mztab.detected_mask[:, ::-1]
    report = compare_states(folder, mztab)
    assert report["checks"]["samples"] == {
        "same_set": True, "same_order": False, "same_length": True, "duplicates": []}
    assert report["checks"]["values"]["n_over"] == 0     # 名前で突き合わせるので値は一致


def test_detected_difference_is_reported(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    folder.detected_mask[0, 0] = not folder.detected_mask[0, 0]
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["detected"]["n_diff"] == 1
    assert report["checks"]["detected"]["skipped"] is False


def test_detected_is_skipped_without_mztab_mask(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    mztab.detected_mask = None
    report = compare_states(folder, mztab)
    assert report["checks"]["detected"] == {"n_compared": 0, "n_diff": 0, "skipped": True}
    assert report["ok"] is True


def test_name_difference_is_reported(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    folder.feature_annotations["1"] = {**folder.feature_annotations["1"], "name": "Other"}
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["names"]["n_diff"] == 1
    assert report["checks"]["names"]["examples"][0]["folder"] == "Other"
    assert report["checks"]["evidence"]["n_diff"] == 0


def test_sample_list_length_or_duplicate_makes_not_ok(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    mztab.sample_names = mztab.sample_names[:3]
    mztab.feature_matrix = mztab.feature_matrix[:, :3]
    mztab.detected_mask = mztab.detected_mask[:, :3]
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["samples"]["same_length"] is False
    folder, mztab = _states(tmp_path / "d", _spec())
    folder.sample_names = folder.sample_names + [folder.sample_names[0]]
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["samples"]["same_set"] is True
    assert report["checks"]["samples"]["duplicates"] == ["s0"]


def test_one_sided_nan_is_counted_and_worst_rel_stays_finite(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    folder.feature_matrix[1, 1] = float("nan")
    report = compare_states(folder, mztab)["checks"]["values"]
    assert report["n_over"] == 1 and report["n_nan_one_side"] == 1 and report["worst_rel"] == 0.0


def test_no_common_features_is_not_ok(tmp_path):
    folder, mztab = _states(tmp_path, _spec())
    mztab.feature_ids = [f"x{f}" for f in mztab.feature_ids]
    report = compare_states(folder, mztab)
    assert report["ok"] is False and report["checks"]["feature_ids"]["n_common"] == 0
    assert report["checks"]["values"]["worst_rel"] is None
