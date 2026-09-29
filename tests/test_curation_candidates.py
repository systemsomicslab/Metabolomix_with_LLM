import pytest

from lipidmix.curation import candidates, judge
from lipidmix.library import store as library_store

TH = judge.resolve_thresholds(None)
SCORING = {"mz_tol": 0.01, "ms2_tol": 0.025, "mass_begin": 0.0, "mass_end": 2000.0,
           "relative_amp_cutoff": 0.0, "absolute_amp_cutoff": 0.0, "use_rt": False}

MSP = """NAME: PC 16:0_18:1
PRECURSORMZ: 804.5760
PRECURSORTYPE: [M+HCOO]-
IONMODE: Negative
INCHIKEY: KEY-PC341
FORMULA: C42H82NO8P
Num Peaks: 2
255.23 999
281.25 800

NAME: PE 18:0_18:2
PRECURSORMZ: 804.5750
PRECURSORTYPE: [M+H]+
INCHIKEY: KEY-PE-POS
Num Peaks: 1
184.07 999

NAME: PS 16:0_18:1
PRECURSORMZ: 804.5700
PRECURSORTYPE: [M-H]-
IONMODE: Negative
INCHIKEY: KEY-PS
Num Peaks: 1
100.00 999
"""


@pytest.fixture()
def store(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "lib.msp"
    path.write_text(MSP, encoding="utf-8")
    s = library_store.open_store(path)
    yield s
    s.close()


def _ev(measured, *, mz=804.5762, ion_mode="Negative"):
    return {"spot_id": 7, "rep_mz": mz, "mz": mz, "rep_rt": 12.0, "rt": 12.0, "ion_mode": ion_mode,
            "_measured": measured}


def test_research_scores_all_in_window_and_filters_polarity(store):
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0], [281.25, 700.0]]), msdial_matches=[], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    names = [c["name"] for c in out["candidates"]]
    # PE [M+H]+ は IONMODE 欄が無い（極性不明）ので store の絞り込みを抜けて検索に掛かり、
    # アダクトの極性が負イオンのスポットと矛盾する（ハード）ので消える。窓は ±0.01 で 3 件とも掛かる
    assert names == ["PC 16:0_18:1", "PS 16:0_18:1"]
    assert out["n_hard_removed"] == 1
    top = out["candidates"][0]
    assert top["candidate_id"] == "L1" and top["source"] == "research"
    assert top["sum_name"] == "PC 34:1" and top["scores"]["matched_peaks_count"] == 2
    assert "no_matched_peaks" in out["candidates"][1]["soft"]
    assert "trend_unknown" in top["info"]


def test_msdial_candidate_merges_with_the_same_research_record(store):
    rec = store.candidates(804.5762, mz_tol=0.01, ion_mode="Negative")
    pc = next(r for r in rec if r["name"] == "PC 16:0_18:1")
    match = {"name": "PC 16:0_18:1", "inchikey": "KEY-PC341", "library_id": pc["record_index"],
             "annotator_id": "lib_1", "total_score": 3.1, "has_msms": True}
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0]]), msdial_matches=[match], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    merged = [c for c in out["candidates"] if c["name"] == "PC 16:0_18:1"]
    assert len(merged) == 1 and merged[0]["source"] == "msdial+research"
    assert merged[0]["msdial_total_score"] == 3.1


def test_current_representative_is_excluded_when_requested(store):
    rec = store.candidates(804.5762, mz_tol=0.01, ion_mode="Negative")
    pc = next(r for r in rec if r["name"] == "PC 16:0_18:1")
    representative = {"name": "PC 16:0_18:1", "library_id": pc["record_index"], "annotator_id": "lib_1",
                      "inchikey": "KEY-PC341"}
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0]]), msdial_matches=[], representative=representative, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=True, top_n=5)
    assert [c["name"] for c in out["candidates"]] == ["PS 16:0_18:1"]


def test_unresolved_msdial_candidate_is_kept_with_its_own_identity(store):
    match = {"name": "PG 34:1", "inchikey": "KEY-PG", "library_id": 999, "annotator_id": "lib_1",
             "total_score": 2.0, "has_msms": True}
    out = candidates.build_library_candidates(
        _ev([]), msdial_matches=[match], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    pg = next(c for c in out["candidates"] if c["name"] == "PG 34:1")
    assert pg["source"] == "msdial" and "reference_unresolved" in pg["info"]
    assert "msms_absent" in pg["soft"] and pg["inchikey"] == "KEY-PG"


def test_no_research_without_msms(store):
    out = candidates.build_library_candidates(
        _ev([]), msdial_matches=[], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    assert out["candidates"] == []


def test_trend_outlier_is_soft_and_ranks_below_clean():
    clean = {"soft": [], "scores": {"total_score": 1.0}}
    outlier = {"soft": ["trend_outlier"], "scores": {"total_score": 3.0}}
    assert candidates.rank([outlier, clean]) == [clean, outlier]


def test_dmz_out_is_hard():
    reasons = candidates.constraint_reasons(
        {"adduct": "[M-H]-", "ontology": "PS", "precursor_mz": 804.555}, rep_mz=804.5762,
        ion_mode="Negative", measured=[[1.0, 1.0]], scores={"matched_peaks_count": 1}, trend_entry=None,
        th=TH)
    assert "dmz_out" in reasons["hard"]
