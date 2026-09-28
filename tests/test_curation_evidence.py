import pytest

from lipidmix.arf2.reader import load_catalog
from lipidmix.curation import evidence, judge
from lipidmix.library import store as library_store
from tests.curation_fixtures import write_alignment_set


@pytest.fixture()
def dataset(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_alignment_set(tmp_path / "neg")
    s = library_store.open_store(paths["msp"])
    yield paths, s
    s.close()


def test_sibling_files_share_the_alignment_stem(dataset):
    paths, _ = dataset
    files = evidence.sibling_files(paths["arf2"])
    assert files["dcl"].name == "AlignmentResult_x.dcl"
    assert files["eic"].name == "AlignmentResult_x.EIC.aef"
    assert files["arf"].name == "AlignmentResult_x_PeakProperties.arf"


def test_select_spots_skips_unknown_and_filters(dataset):
    paths, _ = dataset
    catalog = load_catalog(paths["arf2"])
    assert [s["MasterAlignmentID"] for s in evidence.select_spots(catalog, ontology=None, name_contains=None)] == [0, 1]
    assert [s["MasterAlignmentID"] for s in evidence.select_spots(catalog, ontology=None, name_contains="36:2")] == [1]
    assert evidence.select_spots(catalog, ontology=["PE"], name_contains=None) == []


def test_collect_builds_the_evidence_contract(dataset):
    paths, s = dataset
    catalog = load_catalog(paths["arf2"])
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    th = judge.resolve_thresholds(None)
    evs, stats = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025, th=th)
    by_id = {e["spot_id"]: e for e in evs}
    first = by_id[0]
    assert first["reference"]["name"] == "PC 34:1"
    assert first["ppm_basis"] == "reference"
    assert abs(first["ppm"]) < 1.0
    assert first["drt"] == pytest.approx(0.0, abs=1e-6)
    assert first["mirror"]["measured"] and first["mirror"]["reference"]
    assert first["rescore"]["weighted_dot_product"] > 0.5
    assert [smp["detected"] for smp in first["eic"]["samples"]] == [True, True, False]
    assert first["eic"]["samples"][0]["representative"] is True
    assert first["eic_shape"]["band"] == "PASS"
    second = by_id[1]
    assert second["name_prefix"] == "low score"
    assert second["drt"] == pytest.approx(1.5, abs=1e-3)
    assert stats["n_reference_resolved"] == 2


def test_eic_points_are_trimmed_around_the_peak(dataset):
    paths, s = dataset
    catalog = load_catalog(paths["arf2"])
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    evs, _ = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                              th=judge.resolve_thresholds(None))
    sample = evs[0]["eic"]["samples"][0]
    xs = [p[0] for p in sample["points"]]
    width = sample["right"] - sample["left"]
    assert min(xs) >= sample["left"] - 1.5 * width - 1e-6
    assert max(xs) <= sample["right"] + 1.5 * width + 1e-6


def test_missing_reference_is_reported_not_raised(dataset, tmp_path):
    paths, _ = dataset
    other = tmp_path / "other.msp"
    other.write_text("NAME: X\nPRECURSORMZ: 100.0\nIONMODE: Positive\nNum Peaks: 0\n", encoding="utf-8")
    s = library_store.open_store(other)
    try:
        catalog = load_catalog(paths["arf2"])
        spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
        evs, stats = evidence.collect(paths["arf2"], spots, store=s, ms2_tol=0.025,
                                      th=judge.resolve_thresholds(None))
        assert all(e["reference"] is None for e in evs)
        assert evs[0]["ppm_basis"] == "formula"
        assert stats["n_reference_resolved"] == 0
    finally:
        s.close()
