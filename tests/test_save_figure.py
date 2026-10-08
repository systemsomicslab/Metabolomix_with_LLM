"""save_figure: 図の保存の一本化。旧 4 ツールと同じ出力になることを縛る。"""
from __future__ import annotations

import json

import pytest

import server
from metabolomix.core import session_state

_PCA_PLOT = {"title": "PCA", "x_label": "PC1 (60.00%)", "y_label": "PC2 (20.00%)",
             "points": [{"x": 1.0, "y": 2.0, "label": "a"}, {"x": -1.0, "y": 0.5, "label": "b"},
                        {"x": 0.2, "y": -1.0, "label": "c"}]}


@pytest.fixture(autouse=True)
def _fresh(tmp_path, monkeypatch):
    session_state.session = server.AnalysisSession()
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    yield


def test_unknown_kind_is_an_error():
    out = json.loads(server.save_figure("heatmap", "x"))
    assert out["status"] == "error" and "kind" in out["message"]


@pytest.mark.parametrize("kind", ["eic", "group_intensity", "species", "pca_loadings"])
def test_source_is_rejected_for_kinds_without_sources(kind):
    out = json.loads(server.save_figure(kind, "x", source="arf"))
    assert out["status"] == "error" and "source" in out["message"]
    out = json.loads(server.save_figure(kind, "x", result_id="r1"))
    assert out["status"] == "error"


@pytest.mark.parametrize("kind,state,tools", [
    ("pca", "pca_result", ["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "load_dataset", "dataset_pca"]),
    ("volcano", "differential_result", ["arf_differential", "dataset_differential"]),
    ("eic", "eic_plot", ["eic_plot_chromatograms", "eic_plot_compounds"]),
    ("group_intensity", "group_intensity_plot", ["arf_plot_group_intensity"]),
    ("species", "species_plot", ["arf_plot_species"]),
    ("pca_loadings", "pca_loadings_plot", ["plot_pca_loadings"]),
])
def test_missing_state_per_kind(kind, state, tools):
    err = json.loads(server.save_figure(kind, "x"))["error"]
    assert err["code"] == "missing_state" and err["state"] == state and err["required_tools"] == tools


def test_pca_writes_png_and_reports_absolute_path(tmp_path):
    session_state.session.arf.last_pca_plot = dict(_PCA_PLOT)
    msg = server.save_figure("pca", "a-1")
    png = tmp_path / "reports" / "figures" / "a-1_pca.png"
    assert png.is_file() and str(png) in msg and "figures/a-1_pca.png" in msg


def test_group_intensity_writes_png_and_svg(tmp_path):
    payload = {"plot_schema": "lipidmix.group_intensity.v1", "detection_limit": None,
               "groups": [{"label": "g", "samples": ["s1"]}],
               "items": [{"item": "PG", "n_spots": 1, "n_spots_msms": None, "detected": True,
                          "groups": [{"label": "g", "samples": [
                              {"sample": "s1", "value": 100.0, "gap_filled_fraction": 0.0,
                               "low_reliability": False}],
                              "log10_mean": 2.0, "log10_sd": None, "n_in_stats": 1}]}]}
    session_state.session.arf.last_group_intensity = {"payload": payload, "title": None, "ncols": None}
    msg = server.save_figure("group_intensity", "EV membrane")
    pngs = list((tmp_path / "reports" / "figures").glob("*_group_intensity.png"))
    assert len(pngs) == 1 and pngs[0].with_suffix(".svg").is_file()
    assert "group_intensity.png" in msg
