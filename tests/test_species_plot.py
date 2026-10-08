"""plots/species.py: 分子種ごとの図の payload と描画（純関数）。"""
from __future__ import annotations

import math

import pytest

from metabolomix.plots import species


def _spot(sid, label, adduct="[M-H]-"):
    return {"spot_id": sid, "name": label, "label": label, "ontology": label.split()[0], "item": label.split()[0],
            "adduct": adduct, "mz": 700.0, "rt": 3.0, "msms": True}


def _rows(values, gap=()):
    return [{"file_name": n, "height": h, "is_gap_filled": n in gap} for n, h in values.items()]


GROUPS = [{"label": "C", "samples": ["c1", "c2"]}, {"label": "P", "samples": ["p1", "p2"]}]
ROWS = {
    0: _rows({"c1": 30, "c2": 20, "p1": 10, "p2": 0}, gap=("p2",)),
    1: _rows({"c1": 70, "c2": 80, "p1": 90, "p2": 0}),
}


def test_share_is_height_over_basis_total_and_zero_total_is_listed():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1"), _spot(1, "PG 16:0_19:1")], GROUPS, ROWS,
                                      value="share")
    c = p["spots"][0]["groups"][0]
    assert [s["share"] for s in c["samples"]] == [30.0, 20.0]
    assert c["mean"] == 25.0 and c["sd"] == pytest.approx(7.0711, abs=1e-4)
    assert p["zero_denominator_samples"] == ["p2"]
    p2 = p["spots"][0]["groups"][1]["samples"][1]
    assert p2["share"] is None and p2["gap_filled"] is True
    assert p["share_basis"]["spot_ids"] == [0, 1]


def test_basis_spots_change_the_denominator():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="share",
                                      basis_spots=[_spot(1, "PG 16:0_19:1")], basis_items=["PG 16:0_19:1"])
    assert p["spots"][0]["groups"][0]["samples"][0]["share"] == pytest.approx(30 / 70 * 100, abs=1e-3)
    assert p["share_basis"] == {"items": ["PG 16:0_19:1"], "spot_ids": [1]}


def test_height_stats_are_log10_and_skip_zero_and_low_reliability():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="height",
                                      low_reliability=frozenset({"c2"}))
    c, pg = p["spots"][0]["groups"]
    assert c["n_in_stats"] == 1 and c["mean"] == pytest.approx(math.log10(30), abs=1e-4) and c["sd"] is None
    assert pg["n_in_stats"] == 1        # p2 は 0 なので統計に入らない
    assert p["share_basis"] is None and p["stat"].startswith("mean ± SD of log10")


def test_group_with_only_low_reliability_samples_has_no_mean():
    p = species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="share",
                                      low_reliability=frozenset({"p1", "p2"}))
    assert p["spots"][0]["groups"][1]["mean"] is None
    species.render_species_plot(p)      # 落ちない


def test_unknown_value_is_rejected():
    with pytest.raises(ValueError):
        species.build_species_payload([_spot(0, "PG 16:0_18:1")], GROUPS, ROWS, value="area")


@pytest.mark.parametrize("value", ["share", "height"])
def test_render_draws_one_panel_per_spot_with_adduct_in_title(value):
    import matplotlib.pyplot as plt
    spots = [_spot(0, "DG 16:0_18:1", "[M+NH4]+"), _spot(1, "DG 16:0_18:1", "[M+Na]+")]
    p = species.build_species_payload(spots, GROUPS, ROWS, value=value)
    fig = species.render_species_plot(p, ncols=2)
    titles = [ax.get_title(loc="left") for ax in fig.axes if ax.get_visible() and ax.get_title(loc="left")]
    assert any("[M+NH4]+" in t for t in titles) and any("[M+Na]+" in t for t in titles)
    plt.close(fig)


def test_render_rejects_empty_payload():
    p = species.build_species_payload([], GROUPS, ROWS, value="share")
    with pytest.raises(ValueError):
        species.render_species_plot(p)
