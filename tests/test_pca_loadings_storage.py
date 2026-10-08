"""既存の PCA（ARF / mzTab）がローディング図に要る情報をセッションに残すこと。"""
from __future__ import annotations

import numpy as np

from metabolomix.core import session_state
from metabolomix.core.tool_helpers import _remember_arf_pca_plot


def test_remember_arf_pca_plot_keeps_loadings_and_feature_names():
    session_state.session = session_state.AnalysisSession()
    pca = {"components": [[1.0, 0.5], [-1.0, -0.5], [0.2, 0.1]], "explained_variance_ratio": [0.7, 0.2],
           "singular_values": [3.0, 1.0], "loadings": [[0.6, 0.8], [0.8, -0.6]]}
    _remember_arf_pca_plot(pca, ["a", "b", "c"], title="t", feature_names=["Spot_1_height", "Spot_2_height"])
    plot = session_state.session.arf.last_pca_plot
    assert plot["loadings"] == [[0.6, 0.8], [0.8, -0.6]] and plot["singular_values"] == [3.0, 1.0]
    assert plot["feature_names"] == ["Spot_1_height", "Spot_2_height"]
    assert plot["n_fit"] == 3 and plot["scaling"] == "autoscale"
    assert plot["explained_variance_ratio"] == [0.7, 0.2]
    assert len(plot["points"]) == 3                        # 既存のスコアはそのまま


def test_remember_arf_pca_plot_without_feature_names_stores_none():
    session_state.session = session_state.AnalysisSession()
    pca = {"components": [[1.0, 0.5], [-1.0, -0.5]], "explained_variance_ratio": [0.7, 0.3],
           "singular_values": [3.0, 1.0], "loadings": [[1.0, 0.0], [0.0, 1.0]]}
    _remember_arf_pca_plot(pca, ["a", "b"], title="t")
    assert session_state.session.arf.last_pca_plot["feature_names"] is None


def test_run_dataset_pca_returns_singular_values():
    from metabolomix.analysis.dataset_analysis import run_dataset_pca

    class DS:
        pp_matrix = np.array([[1.0, 2.0, 3.0], [2.0, 1.0, 0.5], [3.0, 5.0, 1.0], [0.5, 0.2, 2.0]])
        pp_sample_names = ["s1", "s2", "s3", "s4"]
        pp_feature_names = ["f1", "f2", "f3"]
        roles = {}
        sample_meta = {}

    out = run_dataset_pca(DS(), n_components=2)
    assert len(out["singular_values"]) == 2 and len(out["loadings"]) == 2
