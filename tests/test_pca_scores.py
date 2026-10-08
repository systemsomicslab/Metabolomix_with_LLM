"""plots/pca_scores.py: PCA スコア図（群の色分け・投影点の白抜き・枠外の矢印）。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from metabolomix.plots.pca_scores import render_pca_scores


def _plot(points):
    return {"title": "t", "x_label": "PC1 (50.0%)", "y_label": "PC2 (20.0%)", "points": points}


def test_groups_get_a_legend_and_projected_points_are_open():
    fig = render_pca_scores(_plot([
        {"x": 1, "y": 1, "label": "c1", "group": "C"}, {"x": 1.2, "y": 0.8, "label": "c2", "group": "C"},
        {"x": -1, "y": -1, "label": "p1", "group": "P"}, {"x": -1.1, "y": -0.9, "label": "p2", "group": "P"},
        {"x": -0.5, "y": 0.2, "label": "p3", "group": "P", "fitted": False},
    ]))
    texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert "C" in texts and "P" in texts and any("projected" in t for t in texts)
    plt.close(fig)


def test_projected_point_outside_the_fitted_range_is_drawn_as_an_arrow():
    fig = render_pca_scores(_plot([
        {"x": 1, "y": 1, "label": "a", "group": "C"}, {"x": -1, "y": -1, "label": "b", "group": "C"},
        {"x": 0.5, "y": -0.5, "label": "c", "group": "P"}, {"x": 40, "y": -60, "label": "far", "group": "P", "fitted": False},
    ]))
    notes = [t.get_text() for t in fig.axes[0].texts]
    assert any("far" in n and "off-scale" in n for n in notes)
    plt.close(fig)


def test_points_without_groups_still_render():
    fig = render_pca_scores(_plot([{"x": 1, "y": 2, "label": "a"}, {"x": -1, "y": 0, "label": "b"}]), title="X")
    assert fig.axes[0].get_title() == "X"
    plt.close(fig)


def test_empty_points_are_rejected():
    with pytest.raises(ValueError):
        render_pca_scores(_plot([]))
