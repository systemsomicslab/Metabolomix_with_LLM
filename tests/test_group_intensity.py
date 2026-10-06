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


def test_manual_exclusion_applies_to_class_and_name_parts():
    items, excluded = gi.resolve_items(["PG", "PG 34:1"], CATALOG, manual={0})
    assert 0 not in _ids(items[0]) and 0 not in _ids(items[1])
    assert _ids(items[1]) == [5]
    assert excluded["manual"] == ["#0 PG 34:1|PG 16:0_18:1"]


def test_render_extends_axis_to_show_a_detection_limit_above_all_values():
    import math
    import matplotlib.pyplot as plt
    fig = gi.render_group_intensity_plot(_payload(detection_limit=1e9))
    try:
        axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
        assert all(ax.get_ylim()[1] >= math.log10(1e9) for ax in axes)
    finally:
        plt.close(fig)


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
    assert len(caveats) == 1
    assert "ko と ko_9w" in caveats[0] and "1 試料" in caveats[0] and "x_ko_9w_1" in caveats[0]


def test_overlap_caveat_is_one_per_group_pair_with_at_most_three_names():
    names = [f"x_ko_9w_{i}" for i in range(1, 6)] + ["x_ctrl_1"]
    groups, caveats = gi.resolve_groups(["ko", "ko_9w"], _facets(names))
    assert len(caveats) == 1
    assert "5 試料" in caveats[0] and "ほか" in caveats[0]
    assert caveats[0].count("x_ko_9w_") == 3


def test_duplicate_group_specs_get_a_caveat():
    groups, caveats = gi.resolve_groups(["ko", "ko"], _facets(["x_ko_1", "x_ctrl_1"]))
    assert len(groups) == 2 and len(caveats) == 1 and "重複" in caveats[0]


def test_role_filtered_group_error_hints_how_to_include_blank():
    facets = _facets(["x_blank_1", "x_blank_2", "x_ctrl_1"])
    with pytest.raises(ValueError) as hinted:
        gi.resolve_groups(["2"], facets)       # 当たるのはブランクだけ（役割で落ちる）
    message = str(hinted.value)
    assert "after role filtering" in message and "blank / qc を書くと" in message
    with pytest.raises(ValueError) as plain:
        gi.resolve_groups(["nonexistent"], facets)
    assert "blank / qc" not in str(plain.value)


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


def _texts(ax):
    return [t.get_text() for t in ax.texts]


def test_render_shares_the_y_axis_and_marks_nd_and_ms1_only():
    import matplotlib.pyplot as plt
    p = _payload(detection_limit=1000.0)
    fig = gi.render_group_intensity_plot(p)
    try:
        axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
        assert len(axes) == 2
        assert axes[0].get_ylim() == axes[1].get_ylim()
        assert any("no annotated species" in t for t in _texts(axes[1]))
        assert any(l.get_linestyle() == "--" for l in axes[0].get_lines())   # 検出下限の破線
    finally:
        plt.close(fig)


def test_render_distinguishes_all_zero_from_no_species():
    import matplotlib.pyplot as plt
    rows = {9: [{"file_name": "c1", "height": 0.0, "is_gap_filled": True}]}
    items = [{"item": "PS", "parts": [{"part": "PS", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 9, "name": "PS 34:1", "ontology": "PS"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "ctrl", "samples": ["c1"]}], rows,
                                         msms={9: False})
    fig = gi.render_group_intensity_plot(p)
    try:
        assert any("intensity 0 in all samples" in t for t in _texts(fig.axes[0]))
    finally:
        plt.close(fig)


def test_render_marks_ms1_only_panels():
    import matplotlib.pyplot as plt
    rows = {1: [{"file_name": "c1", "height": 500.0, "is_gap_filled": False}]}
    items = [{"item": "PG 35:1", "parts": [{"part": "PG 35:1", "kind": "name", "n_spots": 1}],
              "spots": [{"spot_id": 1, "name": "no MS2: PG 35:1", "ontology": "PG"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "ctrl", "samples": ["c1"]}], rows,
                                         msms={1: False})
    fig = gi.render_group_intensity_plot(p)
    try:
        assert any("MS1-only" in t for t in _texts(fig.axes[0]))
    finally:
        plt.close(fig)


def test_render_empty_items_raises():
    p = {"items": [], "groups": [], "plot_schema": "lipidmix.group_intensity.v1",
         "value": "test", "detection_limit": None, "detection_limit_source": None,
         "groups": [], "low_reliability_samples": [], "excluded": {}, "caveats": []}
    with pytest.raises(ValueError, match="描く項目がありません"):
        gi.render_group_intensity_plot(p)


def test_render_diamond_vs_circle_markers():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = {
        0: [{"file_name": "c1", "height": 100.0, "is_gap_filled": False},  # c1: gap_fraction = 0 <= 0.5 → circle
            {"file_name": "c2", "height": 100.0, "is_gap_filled": True}],    # c2: gap_fraction = 1.0 > 0.5 → diamond
        1: [{"file_name": "c1", "height": 50.0, "is_gap_filled": False},
            {"file_name": "c2", "height": 50.0, "is_gap_filled": False}],
    }
    items = [{"item": "item", "parts": [{"part": "item", "kind": "class", "n_spots": 2}],
              "spots": [{"spot_id": 0, "name": "s0", "ontology": "x"},
                        {"spot_id": 1, "name": "s1", "ontology": "x"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "g", "samples": ["c1", "c2"]}], rows,
                                         msms={0: True, 1: True})
    fig = gi.render_group_intensity_plot(p)
    try:
        ax = fig.axes[0]
        collections = ax.collections
        assert len(collections) > 0, "No scatter plots drawn"
        # Find circle (26 vertices) and diamond (5 vertices) by path vertex count
        has_circle = any(len(coll.get_paths()) > 0 and
                        any(len(p.vertices) == 26 for p in coll.get_paths())
                        for coll in collections)
        has_diamond = any(len(coll.get_paths()) > 0 and
                         any(len(p.vertices) == 5 for p in coll.get_paths())
                         for coll in collections)
        assert has_circle and has_diamond, "Both circle and diamond markers should be present"
    finally:
        plt.close(fig)


def test_render_low_reliability_sample_white_fill():
    import matplotlib.pyplot as plt
    import numpy as np
    # Two samples: one normal (c1) and one low-reliability (c2)
    rows = {0: [{"file_name": "c1", "height": 100.0, "is_gap_filled": False},
                {"file_name": "c2", "height": 100.0, "is_gap_filled": False}]}
    items = [{"item": "item", "parts": [{"part": "item", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 0, "name": "s0", "ontology": "x"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "g", "samples": ["c1", "c2"]}], rows,
                                         msms={0: True}, low_reliability=frozenset({"c2"}))
    fig = gi.render_group_intensity_plot(p)
    try:
        ax = fig.axes[0]
        collections = ax.collections
        found_white = False
        found_non_white = False
        for coll in collections:
            facecolors = coll.get_facecolors()
            if len(facecolors) > 0:
                for fc in facecolors:
                    # Check for white fill
                    if np.allclose(fc[:3], [1, 1, 1], atol=0.01):
                        found_white = True
                    # Check for non-white fill (should be the group color)
                    elif not np.allclose(fc[:3], [1, 1, 1], atol=0.01) and not np.allclose(fc[:3], [0, 0, 0], atol=0.01):
                        # Not white, not black (black is edge color)
                        found_non_white = True
        assert found_white, "Low-reliability sample should have white facecolor"
        assert found_non_white, "Normal sample should have non-white facecolor"
    finally:
        plt.close(fig)


def test_render_no_dashed_line_without_detection_limit():
    import matplotlib.pyplot as plt
    rows = {0: [{"file_name": "c1", "height": 1000.0, "is_gap_filled": False}]}
    items = [{"item": "item", "parts": [{"part": "item", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 0, "name": "s0", "ontology": "x"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "g", "samples": ["c1"]}], rows, msms={0: True})
    assert p["detection_limit"] is None
    fig = gi.render_group_intensity_plot(p)
    try:
        ax = fig.axes[0]
        for line in ax.get_lines():
            assert line.get_linestyle() != "--", "Dashed line should not be present when detection_limit is None"
    finally:
        plt.close(fig)


def test_render_group_mean_line_drawn():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Use 2 samples so log10_sd is not None (otherwise errorbar is not drawn)
    rows = {0: [{"file_name": "c1", "height": 100.0, "is_gap_filled": False},
                {"file_name": "c2", "height": 1000.0, "is_gap_filled": False}]}
    items = [{"item": "item", "parts": [{"part": "item", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 0, "name": "s0", "ontology": "x"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "g", "samples": ["c1", "c2"]}], rows, msms={0: True})
    # Verify log10_sd is not None (i.e., we have 2+ samples for statistics)
    assert p["items"][0]["groups"][0]["log10_sd"] is not None, "Test fixture must have 2+ samples for SD calculation"
    log10_mean = p["items"][0]["groups"][0]["log10_mean"]
    log10_sd = p["items"][0]["groups"][0]["log10_sd"]

    fig = gi.render_group_intensity_plot(p)
    try:
        ax = fig.axes[0]
        # Look for LineCollection with horizontal segment at y == log10_mean (from hlines)
        mean_line_found = False
        for coll in ax.collections:
            if hasattr(coll, 'get_segments'):
                for seg in coll.get_segments():
                    # A horizontal line has same y-coords at both endpoints
                    if len(seg) == 2 and pytest.approx(seg[0][1], abs=0.01) == pytest.approx(seg[1][1], abs=0.01):
                        if pytest.approx(seg[0][1], abs=0.01) == pytest.approx(log10_mean, abs=0.01):
                            mean_line_found = True
                            break
        assert mean_line_found, "Mean line (horizontal segment at y=log10_mean) should be drawn"

        # Look for ErrorbarContainer or vertical lines for SD bars
        sd_bar_found = False
        # Check containers for ErrorbarContainer
        if hasattr(ax, 'containers'):
            for container in ax.containers:
                if hasattr(container, 'lines') and len(container.lines) > 0:
                    sd_bar_found = True
                    break
        # Also check for vertical segments in LineCollections
        if not sd_bar_found:
            for coll in ax.collections:
                if hasattr(coll, 'get_segments'):
                    for seg in coll.get_segments():
                        # A vertical line has same x-coords at both endpoints
                        if len(seg) == 2 and pytest.approx(seg[0][0], abs=0.01) == pytest.approx(seg[1][0], abs=0.01):
                            y_min, y_max = min(seg[0][1], seg[1][1]), max(seg[0][1], seg[1][1])
                            # SD bar should span mean ± SD
                            if (pytest.approx(y_min, abs=0.01) == pytest.approx(log10_mean - log10_sd, abs=0.01) and
                                pytest.approx(y_max, abs=0.01) == pytest.approx(log10_mean + log10_sd, abs=0.01)):
                                sd_bar_found = True
                                break
        assert sd_bar_found, "SD bar should be drawn when log10_sd is not None"
    finally:
        plt.close(fig)


def test_render_shared_yaxis_with_different_magnitudes():
    import matplotlib.pyplot as plt
    rows = {
        0: [{"file_name": "c1", "height": 10.0, "is_gap_filled": False}],         # small: log10(10) = 1
        1: [{"file_name": "c1", "height": 100000.0, "is_gap_filled": False}],     # large: log10(100000) = 5
    }
    items = [{"item": "small", "parts": [{"part": "small", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 0, "name": "s0", "ontology": "x"}]},
             {"item": "large", "parts": [{"part": "large", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 1, "name": "s1", "ontology": "x"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "g", "samples": ["c1"]}], rows, msms={0: True, 1: True})
    fig = gi.render_group_intensity_plot(p)
    try:
        axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
        assert len(axes) == 2
        ylim0 = axes[0].get_ylim()
        ylim1 = axes[1].get_ylim()
        # Both should have the same limits (shared y-axis)
        assert ylim0 == ylim1, f"Y-axes should be shared but got {ylim0} vs {ylim1}"
        # Upper limit should be >= log10(100000) = 5
        assert ylim0[1] >= 5.0, f"Upper y-limit {ylim0[1]} should be >= 5.0 to accommodate both items"
    finally:
        plt.close(fig)
