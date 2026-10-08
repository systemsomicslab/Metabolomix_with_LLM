"""plot_pca_loadings（ツール層）: 3 経路の PCA のローディング図。"""
from __future__ import annotations

import json
import math

import pytest

import server
from metabolomix.core import session_state


def _arf_plot(n_feat=4, with_loadings=True):
    plot = {"title": "PCA", "x_label": "PC1", "y_label": "PC2",
            "points": [{"x": 1.0, "y": 0.0, "label": "a"}, {"x": -1.0, "y": 0.5, "label": "b"},
                       {"x": 0.0, "y": -0.5, "label": "c"}, {"x": 0.3, "y": 0.2, "label": "d"}],
            "provenance": {"result_id": "arf-1", "dataset_id": "d"}}
    if with_loadings:
        plot.update({"loadings": [[0.5, -0.5, 0.5, -0.5][:n_feat], [0.1, 0.2, -0.3, 0.9][:n_feat]],
                     "singular_values": [4.0, 2.0], "explained_variance_ratio": [0.7, 0.2],
                     "feature_names": [f"Spot_{i}_height" for i in range(n_feat)], "n_fit": 4,
                     "scaling": "autoscale"})
    return plot


@pytest.fixture(autouse=True)
def _fresh():
    session_state.session = server.AnalysisSession()
    session_state.session.arf.features = [{"MasterAlignmentID": i, "Name": f"PG {30 + i}:0",
                                           "MassCenter": 700.0 + i, "RT": 3.0} for i in range(4)]
    yield


def _call(**kw):
    from metabolomix.tools.pca_loadings_tools import plot_pca_loadings
    return json.loads(plot_pca_loadings(output="payload", **kw))


def test_no_pca_is_missing_state():
    err = _call()["error"]
    assert err["code"] == "missing_state" and "arf_pca_species" in err["required_tools"]


def test_old_arf_result_without_loadings_is_not_a_candidate():
    session_state.session.arf.last_pca_plot = _arf_plot(with_loadings=False)
    assert _call()["error"]["code"] == "missing_state"


def test_arf_source_uses_reader_selection_and_computes_r():
    session_state.session.arf.last_pca_plot = _arf_plot()
    p = _call(top_n=1)
    assert p["source"] == "arf" and p["value"] == "r" and p["result_id"] == "arf-1"
    rows = {row["feature_id"]: row for row in p["panels"][0]["rows"]}
    assert set(rows) == {"0", "1"} or len(rows) == 2
    some = next(iter(rows.values()))
    assert some["r"] == pytest.approx(some["coefficient"] * 4.0 / math.sqrt(4), abs=1e-3)
    assert some["label"].startswith("PG")


def test_species_source_uses_stored_rows():
    session_state.session.arf.last_species_pca = {
        "points": [{"x": 0, "y": 0, "label": "a"}], "explained_variance_ratio": [0.8, 0.1],
        "loadings_rows": [{"feature_id": "7", "label": "DGDG 16:0_18:1", "ontology": "DGDG", "adduct": "[M+CH3COO]-",
                           "mz": 977.6, "rt": 4.1, "coefficient": [0.7, 0.1], "r": [0.99, 0.05]},
                          {"feature_id": "8", "label": "PG 16:0_19:1", "ontology": "PG", "adduct": "[M-H]-",
                           "mz": 761.5, "rt": 3.7, "coefficient": [-0.7, 0.1], "r": [-0.98, 0.04]}],
        "provenance": {"result_id": "sp-1", "dataset_id": "d"}}
    p = _call(top_n=None)
    assert p["source"] == "species" and p["aligned"] is True
    assert [row["feature_id"] for row in p["panels"][0]["rows"]] == ["7", "8"]


def test_two_sources_need_a_choice_and_bad_pcs_is_an_error():
    session_state.session.arf.last_pca_plot = _arf_plot()
    session_state.session.arf.last_species_pca = {
        "points": [{"x": 0, "y": 0, "label": "a"}], "explained_variance_ratio": [0.8, 0.1],
        "loadings_rows": [{"feature_id": "7", "label": "x", "ontology": None, "adduct": "", "mz": None, "rt": None,
                           "coefficient": [0.7, 0.1], "r": [0.9, 0.1]},
                          {"feature_id": "8", "label": "y", "ontology": None, "adduct": "", "mz": None, "rt": None,
                           "coefficient": [-0.7, 0.1], "r": [-0.9, 0.1]}],
        "provenance": {"result_id": "sp-1", "dataset_id": "d"}}
    assert _call()["error"]["code"] == "AMBIGUOUS_RESULT_SOURCE"
    out = _call(source="arf", pcs=[4])
    assert out["status"] == "error" and "pcs" in out["message"]


def test_image_and_save_figure(tmp_path, monkeypatch):
    from mcp.server.fastmcp import Image
    from metabolomix.tools.pca_loadings_tools import plot_pca_loadings
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    session_state.session.arf.last_pca_plot = _arf_plot()
    out = plot_pca_loadings()
    assert isinstance(out[1], Image) and "source=arf" in out[0]
    msg = server.save_figure("pca_loadings", "L1")
    # make_slug は analysis_id を小文字にする
    png = tmp_path / "reports" / "figures" / "l1_pca_loadings.png"
    assert png.is_file() and png.with_suffix(".svg").is_file() and str(png) in msg


def test_failed_call_clears_previous_figure():
    session_state.session.arf.last_pca_plot = _arf_plot()
    _call()
    assert session_state.session.last_loadings_plot is not None
    _call(value="loading")
    assert session_state.session.last_loadings_plot is None


def test_mztab_source_draws_dataset_pca():
    from tests.test_dataset_analysis_tools import _load_ds
    from metabolomix.tools.dataset_analysis_tools import dataset_pca, dataset_preprocess
    _load_ds()
    dataset_preprocess()
    dataset_pca(n_components=2)
    p = _call(source="mztab", top_n=3)
    assert p["source"] == "mztab" and p["value"] == "r"
    assert p["panels"][0]["rows"][0]["label"].startswith("Compound")
