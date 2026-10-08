"""arf_plot_species / arf_pca_species（ツール層）。fixture はテストが作る。"""
from __future__ import annotations

import json

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ctrl_3", "ko_1", "ko_2", "ko_3"]


def _group(heights, gap=()):
    rows = []
    for fid, (name, h) in enumerate(zip(SAMPLES, heights)):
        row = arf_row(file_id=fid, mz=700.0, rt=3.0, height=h, gap_filled=name in gap)
        row[1] = f"20261006_{name}"
        rows.append(row)
    return rows


@pytest.fixture()
def loaded(tmp_path):
    arf = write_arf(tmp_path / f"{STEM}_PeakProperties.arf", [
        _group([10, 5000, 6000, 5500, 2000, 2500, 2200]),     # 0 PG 16:0_18:1（MS/MS）
        _group([10, 3000, 3500, 3200, 4000, 4200, 4100]),     # 1 PG 16:0_19:1（MS/MS なし）
        _group([0, 2000, 2100, 1900, 4000, 4300, 300]),       # 2 DGDG 16:0_18:1（MS/MS）
        _group([0, 900, 950, 920, 800, 700, 760]),            # 3 DG 16:0_18:1 [M+NH4]+（MS/MS）
        _group([0, 300, 320, 310, 280, 260, 270]),            # 4 DG 16:0_18:1 [M+Na]+（MS/MS）
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG", adduct="[M-H]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG", adduct="[M-H]-"),
        arf2_spot_raw(spot_id=2, name="DGDG 34:1|DGDG 16:0_18:1", ontology="DGDG", adduct="[M+CH3COO]-",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=3, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+NH4]+",
                      matches=[match_result()]),
        arf2_spot_raw(spot_id=4, name="DG 34:1|DG 16:0_18:1", ontology="DG", adduct="[M+Na]+",
                      matches=[match_result()]),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _species(**kw):
    from metabolomix.arf.species_tools import arf_plot_species
    kw.setdefault("groups", ["ctrl", "ko"])
    return json.loads(arf_plot_species(output="payload", **kw))


def test_missing_arf_returns_missing_state():
    from metabolomix.arf.species_tools import arf_plot_species
    session_state.session = session_state.AnalysisSession()
    out = json.loads(arf_plot_species(items=["PG"], groups=["ctrl"], output="payload"))
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["required_tools"] == ["arf_parser", "load_dataset"]


def test_share_payload_expands_classes_and_uses_items_as_default_basis(loaded):
    p = _species(items=["PG", "DGDG"])
    assert p["plot_schema"] == "lipidmix.species_intensity.v1" and p["value"] == "share"
    assert [s["spot_id"] for s in p["spots"]] == [0, 1, 2]
    assert p["share_basis"]["spot_ids"] == [0, 1, 2]
    c1 = p["spots"][0]["groups"][0]["samples"][0]
    assert c1["share"] == pytest.approx(5000 / (5000 + 3000 + 2000) * 100, abs=1e-3)


def test_share_basis_and_require_msms(loaded):
    p = _species(items=["PG"], share_basis=["PG", "DGDG"], require_msms=True)
    assert [s["spot_id"] for s in p["spots"]] == [0]
    assert p["share_basis"]["spot_ids"] == [0, 2]        # 分母にも MS/MS の条件が効く
    assert p["excluded"]["no_msms"] == ["#1 PG 35:1|PG 16:0_19:1"]


def test_basis_without_the_plotted_species_adds_caveat(loaded):
    p = _species(items=["PG"], share_basis=["DGDG"])
    assert any("分母に含まれない" in c for c in p["caveats"])


def test_same_name_different_adduct_are_separate_panels(loaded):
    p = _species(items=["DG"], value="height")
    assert [(s["spot_id"], s["adduct"]) for s in p["spots"]] == [(3, "[M+NH4]+"), (4, "[M+Na]+")]


def test_too_many_panels_is_an_error_and_clears_slot(loaded, monkeypatch):
    from metabolomix.plots import species
    monkeypatch.setattr(species, "MAX_PANELS", 2)
    _species(items=["PG"])
    assert session_state.session.arf.last_species_plot is not None
    out = _species(items=["PG", "DGDG"])
    assert out["status"] == "error" and "40" not in out["message"] and "2" in out["message"]
    assert session_state.session.arf.last_species_plot is None


@pytest.mark.parametrize("kw", [{"value": "area"}, {"ncols": 0}, {"groups": ["nonexistent"]}])
def test_invalid_input_returns_error(loaded, kw):
    out = _species(items=["PG"], **kw)
    assert out["status"] == "error" and out["message"]


def test_image_output_caption_and_session(loaded):
    from mcp.server.fastmcp import Image
    from metabolomix.arf.species_tools import arf_plot_species
    out = arf_plot_species(items=["PG"], groups=["ctrl", "ko"], low_reliability_samples=["ko_3"])
    assert isinstance(out, list) and isinstance(out[1], Image)
    assert "2 パネル" in out[0] and "ctrl n=3" in out[0]
    assert session_state.session.arf.last_species_plot["payload"]["spots"][0]["spot_id"] == 0


def test_save_figure_species_writes_png_and_svg(loaded, tmp_path, monkeypatch):
    import server
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    _species(items=["PG"])
    msg = server.save_figure("species", "EV-PG")
    png = tmp_path / "reports" / "figures" / "ev-pg_species.png"      # make_slug は小文字にする
    assert png.is_file() and png.with_suffix(".svg").is_file() and str(png) in msg


def _pca(**kw):
    from metabolomix.arf.species_tools import arf_pca_species
    kw.setdefault("groups", ["ctrl", "ko"])
    return json.loads(arf_pca_species(output="payload", **kw))


def test_pca_species_projects_low_reliability_and_orients(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], low_reliability_samples=["ko_3"], orient_by="ko")
    pts = {pt["label"]: pt for pt in p["points"]}
    assert pts["20261006_ko_3"]["fitted"] is False and pts["20261006_ctrl_1"]["fitted"] is True
    ko_x = [pts[f"20261006_ko_{i}"]["x"] for i in (1, 2)]
    assert sum(ko_x) / 2 > 0                                    # orient_by の群が正
    last = session_state.session.arf.last_species_pca
    assert last["n_fit"] == 5 and last["scaling"] == "autoscale"
    assert {row["feature_id"] for row in last["loadings_rows"]} == {"0", "1", "2", "3", "4"}
    assert all(-1.0001 <= v <= 1.0001 for row in last["loadings_rows"] for v in row["r"])


def test_pca_species_total_normalization_with_zero_total_names_the_sample(loaded):
    session_state.session.arf.excluded_spots = set()
    out = _pca(items=["DG"], groups=["blank", "ctrl"], normalize="total")
    assert out["status"] == "error" and "blank_1" in out["message"]


def test_pca_species_sample_in_two_groups_is_counted_once(loaded):
    p = _pca(items=["PG", "DGDG", "DG"], groups=["ctrl", "ctrl_1", "ko"])
    labels = [pt["label"] for pt in p["points"]]
    assert labels.count("20261006_ctrl_1") == 1


@pytest.mark.parametrize("kw", [{"normalize": "median"}, {"orient_by": "nope"}, {"items": ["PG 16:0_18:1"]},
                                {"groups": ["ctrl"], "low_reliability_samples": ["ctrl_1"]}])
def test_pca_species_invalid_input_returns_error_and_clears_slot(loaded, kw):
    _pca(items=["PG", "DGDG", "DG"])
    assert session_state.session.arf.last_species_pca is not None
    out = _pca(**{"items": ["PG", "DGDG", "DG"], **kw})
    assert out["status"] == "error" and out["message"]
    assert session_state.session.arf.last_species_pca is None


def test_pca_species_image_and_save_figure_needs_source_when_ambiguous(loaded, tmp_path, monkeypatch):
    from mcp.server.fastmcp import Image
    import server
    from metabolomix.arf.species_tools import arf_pca_species
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    out = arf_pca_species(items=["PG", "DGDG", "DG"], groups=["ctrl", "ko"])
    assert isinstance(out[1], Image) and "PC1" in out[0]
    session_state.session.arf.last_pca_plot = {"title": "x", "x_label": "PC1", "y_label": "PC2",
                                               "points": [{"x": 0, "y": 0, "label": "a"}]}
    ambiguous = json.loads(server.save_figure("pca", "s1"))
    assert ambiguous["error"]["code"] == "AMBIGUOUS_RESULT_SOURCE"
    msg = server.save_figure("pca", "s1", source="species")
    assert (tmp_path / "reports" / "figures" / "s1_pca.png").is_file() and "source=species" in msg
