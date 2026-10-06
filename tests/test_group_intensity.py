"""群別強度プロット（metabolomix/plots/group_intensity.py ほか）の純関数テスト。"""
from __future__ import annotations

import pytest

from metabolomix.msdial import analysis_params, sample_factors
from metabolomix.plots import group_intensity as gi


def test_param_file_min_peak_height_is_read(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\nMinimum peak height: 1000\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] == 1000.0


def test_param_file_without_min_peak_height_gives_none(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] is None


CATALOG = [
    {"MasterAlignmentID": 0, "Name": "PG 34:1|PG 16:0_18:1", "Ontology": "PG"},
    {"MasterAlignmentID": 1, "Name": "low score: PG 35:1|PG 16:0_19:1", "Ontology": "PG"},
    {"MasterAlignmentID": 2, "Name": "PC 33:1(d7)|PC 15:0_18:1(d7)", "Ontology": "PC"},
    {"MasterAlignmentID": 3, "Name": "PC 34:1", "Ontology": "PC"},
    {"MasterAlignmentID": 4, "Name": "LPG 34:1", "Ontology": "LPG"},
    {"MasterAlignmentID": 5, "Name": "PG 34:1", "Ontology": "PG"},          # 同じ名前の別スポット
    {"MasterAlignmentID": 6, "Name": "PE O-34:1", "Ontology": "EtherPE"},
    {"MasterAlignmentID": 7, "Name": "Unknown", "Ontology": ""},
]


def _ids(item):
    return [s["spot_id"] for s in item["spots"]]


def test_class_item_matches_ontology_case_insensitively():
    items, _ = gi.resolve_items(["pg"], CATALOG)
    assert _ids(items[0]) == [0, 1, 5]
    assert items[0]["parts"] == [{"part": "pg", "kind": "class", "n_spots": 3}]


def test_name_item_matches_any_candidate_exactly_and_sums_duplicates():
    items, _ = gi.resolve_items(["PG 34:1", "pg 16:0_19:1"], CATALOG)
    assert _ids(items[0]) == [0, 5]               # LPG 34:1 には当たらない（部分一致しない）
    assert _ids(items[1]) == [1]                  # 接頭辞 low score: を除いた候補名
    assert items[0]["parts"][0]["kind"] == "name"


def test_plus_joins_parts_and_unknown_part_is_none():
    items, _ = gi.resolve_items(["PE+EtherPE", "SQDG"], CATALOG)
    assert _ids(items[0]) == [6]
    assert [p["kind"] for p in items[0]["parts"]] == ["none", "class"]
    assert items[1]["spots"] == [] and items[1]["parts"][0]["kind"] == "none"


def test_internal_standard_is_dropped_from_class_but_kept_by_name():
    items, excluded = gi.resolve_items(["PC", "PC 33:1(d7)"], CATALOG)
    assert _ids(items[0]) == [3]
    assert _ids(items[1]) == [2]
    assert excluded["internal_standard"] == ["#2 PC 33:1(d7)|PC 15:0_18:1(d7)"]


def test_curation_wrong_redundant_and_auto_likely_wrong_are_excluded():
    curation = {"wrong": {0}, "redundant": {5}, "assign": {}, "auto_likely_wrong": {1}}
    items, excluded = gi.resolve_items(["PG"], CATALOG, curation=curation)
    assert items[0]["spots"] == []
    assert excluded["curation"] == ["#0 PG 34:1|PG 16:0_18:1", "#5 PG 34:1"]
    assert excluded["auto_likely_wrong"] == ["#1 low score: PG 35:1|PG 16:0_19:1"]


def test_assign_moves_spot_to_the_recorded_class():
    curation = {"assign": {3: {"name": "PE 34:1", "ontology": "PE"}}}
    items, _ = gi.resolve_items(["PE", "PC"], CATALOG, curation=curation)
    assert items[0]["spots"] == [{"spot_id": 3, "name": "PE 34:1", "ontology": "PE"}]
    assert _ids(items[1]) == []                   # 内部標準 #2 は除かれ、#3 は PE へ移った


def test_standard_only_is_dropped_from_class_only():
    items, excluded = gi.resolve_items(["PG", "PG 34:1"], CATALOG, standard_only=frozenset({0}))
    assert _ids(items[0]) == [1, 5]
    assert _ids(items[1]) == [0, 5]
    assert excluded["standard_only"] == ["#0 PG 34:1|PG 16:0_18:1"]


def test_item_count_limits():
    with pytest.raises(ValueError):
        gi.resolve_items([], CATALOG)
    with pytest.raises(ValueError):
        gi.resolve_items(["PG"] * 31, CATALOG)
    with pytest.raises(ValueError):
        gi.resolve_items(["+"], CATALOG)


def test_standard_only_spots_compares_max_heights():
    rows = {0: [{"file_name": "std_1", "height": 5000.0}, {"file_name": "ctrl_1", "height": 100.0}],
            1: [{"file_name": "std_1", "height": 5000.0}, {"file_name": "ctrl_1", "height": 900.0}]}
    assert gi.standard_only_spots(rows, {"std_1"}, {"ctrl_1"}) == frozenset({0})


def _facets(names):
    return sample_factors.build_sample_facets(names)


def test_groups_keep_order_and_skip_blank_unless_named():
    facets = _facets(["x_blank_1", "x_ctrl_1", "x_ctrl_2", "x_ko_1", "x_ko_2"])
    groups, caveats = gi.resolve_groups(["ko", "ctrl"], facets)
    assert [g["label"] for g in groups] == ["ko", "ctrl"]
    assert groups[1]["samples"] == ["x_ctrl_1", "x_ctrl_2"]
    groups, _ = gi.resolve_groups(["blank", "ctrl"], facets)
    assert groups[0]["samples"] == ["x_blank_1"]
    assert caveats == []


def test_overlapping_groups_are_drawn_in_both_with_a_caveat():
    facets = _facets(["x_ko_9w_1", "x_ko_24m_1", "x_ctrl_9w_1"])
    groups, caveats = gi.resolve_groups(["ko", "ko_9w"], facets)
    assert groups[1]["samples"] == ["x_ko_9w_1"]
    assert any("x_ko_9w_1" in c for c in caveats)


def test_unmatched_or_empty_group_raises():
    facets = _facets(["x_ctrl_1"])
    with pytest.raises(ValueError):
        gi.resolve_groups(["ko"], facets)
    with pytest.raises(ValueError):
        gi.resolve_groups([], facets)
    with pytest.raises(ValueError):
        gi.resolve_groups(["_"], facets)


def test_low_reliability_specs_accept_names_and_tokens():
    facets = _facets(["x_ctrl_1", "x_ko_1", "x_ko_2"])
    assert gi.resolve_sample_specs(["x_ko_1"], facets) == {"x_ko_1"}
    assert gi.resolve_sample_specs(["ko"], facets) == {"x_ko_1", "x_ko_2"}


def _rows():
    return {
        0: [{"file_name": "c1", "height": 100.0, "is_gap_filled": False},
            {"file_name": "c2", "height": 1000.0, "is_gap_filled": False},
            {"file_name": "k1", "height": 10.0, "is_gap_filled": True}],
        1: [{"file_name": "c1", "height": 900.0, "is_gap_filled": True},
            {"file_name": "c2", "height": 0.0, "is_gap_filled": True},
            {"file_name": "k1", "height": 0.0, "is_gap_filled": True}],
    }


def _payload(**kw):
    items = [{"item": "PG", "parts": [{"part": "PG", "kind": "class", "n_spots": 2}],
              "spots": [{"spot_id": 0, "name": "PG 34:1", "ontology": "PG"},
                        {"spot_id": 1, "name": "no MS2: PG 35:1", "ontology": "PG"}]},
             {"item": "PE", "parts": [{"part": "PE", "kind": "none", "n_spots": 0}], "spots": []}]
    groups = [{"label": "ctrl", "samples": ["c1", "c2"]}, {"label": "ko", "samples": ["k1"]}]
    return gi.build_group_intensity_payload(items, groups, _rows(), msms={0: True, 1: False}, **kw)


def test_values_are_sums_with_gap_fill_fraction():
    p = _payload()
    ctrl = p["items"][0]["groups"][0]["samples"]
    assert ctrl[0] == {"sample": "c1", "value": 1000.0, "gap_filled_fraction": 0.9, "low_reliability": False}
    assert ctrl[1]["value"] == 1000.0 and ctrl[1]["gap_filled_fraction"] == 0.0
    assert p["plot_schema"] == "lipidmix.group_intensity.v1"


def test_log10_mean_sd_and_low_reliability_exclusion():
    p = _payload()
    g = p["items"][0]["groups"][0]
    assert g["log10_mean"] == 3.0 and g["log10_sd"] == 0.0 and g["n_in_stats"] == 2
    p = _payload(low_reliability=frozenset({"c2"}))
    g = p["items"][0]["groups"][0]
    assert g["n_in_stats"] == 1 and g["log10_sd"] is None
    assert g["samples"][1]["low_reliability"] is True


def test_zero_values_are_kept_and_not_in_stats():
    p = _payload()
    k = p["items"][0]["groups"][1]
    assert k["samples"][0]["value"] == 10.0
    assert k["n_in_stats"] == 1
    pe = p["items"][1]
    assert pe["detected"] is False and pe["n_spots"] == 0
    assert pe["groups"][0]["samples"][0]["value"] == 0.0 and pe["groups"][0]["log10_mean"] is None


def test_msms_count_and_metadata_are_carried():
    p = _payload(detection_limit=1000.0, detection_limit_source="argument",
                 excluded={"internal_standard": ["#2 x"]}, caveats=["c"])
    assert p["items"][0]["n_spots_msms"] == 1
    assert p["detection_limit"] == 1000.0 and p["detection_limit_source"] == "argument"
    assert p["excluded"]["internal_standard"] == ["#2 x"] and p["caveats"] == ["c"]
    assert p["groups"] == [{"label": "ctrl", "samples": ["c1", "c2"]}, {"label": "ko", "samples": ["k1"]}]
