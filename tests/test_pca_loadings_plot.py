"""plots/pca_loadings.py: ローディングの選び方と描画（純関数）。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from metabolomix.plots import pca_loadings as pl


def _feat(i, c1, c2, ontology="PG"):
    return {"feature_id": str(i), "label": f"PG {i}:0", "ontology": ontology, "adduct": "[M-H]-",
            "mz": 700.0, "rt": 3.0, "coefficient": [c1, c2], "r": [c1 * 2, c2 * 2]}


FEATS = [_feat(i, c1, c2) for i, (c1, c2) in enumerate(
    [(0.5, 0.1), (0.4, -0.2), (0.1, 0.3), (-0.1, 0.05), (-0.3, -0.4), (-0.45, 0.2)])]


def _payload(**kw):
    args = dict(source="species", result_id="r1", explained_variance_ratio=[0.6, 0.2], pcs=[1, 2],
                top_n=2, value="r", r_available=True)
    args.update(kw)
    return pl.build_loadings_payload(FEATS, **args)


def test_top_n_takes_positive_and_negative_per_pc():
    p = _payload()
    pc1 = [row["feature_id"] for row in p["panels"][0]["rows"]]
    assert pc1 == ["0", "1", "4", "5"] or set(pc1) == {"0", "1", "4", "5"}
    assert p["panels"][0]["explained_pct"] == 60.0 and p["aligned"] is False
    assert p["panels"][0]["rows"][0]["value"] == p["panels"][0]["rows"][0]["r"]


def test_top_n_larger_than_half_dedupes():
    p = _payload(top_n=5)
    ids = [row["feature_id"] for row in p["panels"][0]["rows"]]
    assert len(ids) == len(set(ids)) == 6


def test_all_features_are_aligned_by_the_first_pc():
    p = _payload(top_n=None)
    first = [row["feature_id"] for row in p["panels"][0]["rows"]]
    second = [row["feature_id"] for row in p["panels"][1]["rows"]]
    assert p["aligned"] is True and first == second
    values = [row["value"] for row in p["panels"][0]["rows"]]
    assert values == sorted(values, reverse=True)


def test_all_features_over_the_limit_is_rejected(monkeypatch):
    monkeypatch.setattr(pl, "MAX_ALL_FEATURES", 5)
    with pytest.raises(ValueError):
        _payload(top_n=None)


def test_r_requested_but_unavailable_falls_back_to_coefficient():
    p = _payload(r_available=False)
    assert p["value"] == "coefficient" and p["value_requested"] == "r"
    assert any("coefficient" in c for c in p["caveats"])


@pytest.mark.parametrize("pcs", [[3], [0], [1, 1], [1, 2, 3, 4]])
def test_invalid_pcs_are_rejected(pcs):
    with pytest.raises(ValueError):
        _payload(pcs=pcs)


def test_render_one_panel_per_pc_with_class_legend():
    fig = pl.render_loadings_plot(_payload(top_n=None))
    assert len([ax for ax in fig.axes if ax.get_visible()]) == 2
    plt.close(fig)


def test_class_palette_has_20_distinct_colours_and_keeps_the_first_five():
    from metabolomix.plots import pca_loadings as pl
    assert len(set(pl._CLASS_PALETTE)) == len(pl._CLASS_PALETTE) == 20
    assert pl._CLASS_PALETTE[:5] == ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
