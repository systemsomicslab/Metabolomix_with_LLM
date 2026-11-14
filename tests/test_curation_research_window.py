"""`curation_suggest` の再検索（`research`）の前駆体の窓は、`library_match_feature` と同じ MS-DIAL の
`FixMassTolerance`（m/z 500 を超えたら同じ ppm に広げる。事実確認 A3）。"""
from __future__ import annotations

import pytest

from metabolomix.analysis.spectral_match import fix_mass_tolerance
from metabolomix.curation import candidates, judge
from metabolomix.library import store as library_store
from tests.config_fixtures import use_config

SCORING = {"mz_tol": 0.01, "ms2_tol": 0.025, "mass_begin": 0.0, "mass_end": 2000.0,
           "relative_amp_cutoff": 0.0, "absolute_amp_cutoff": 0.0, "use_rt": False}

# 合成の記録。m/z 900 で 0.015 Da 離れた参照（0.01 Da の窓の外、20 ppm ＝ 0.018 Da の窓の中）
MSP = """NAME: PC 40:6
PRECURSORMZ: 900.0150
PRECURSORTYPE: [M+HCOO]-
IONMODE: Negative
INCHIKEY: KEY-WIDE
Num Peaks: 1
255.23 999
"""


@pytest.fixture()
def store(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    path = tmp_path / "lib.msp"
    path.write_text(MSP, encoding="utf-8")
    s = library_store.open_store(path)
    yield s
    s.close()


def _ev():
    return {"spot_id": 7, "rep_mz": 900.0, "mz": 900.0, "rep_rt": 12.0, "rt": 12.0,
            "ion_mode": "Negative", "_measured": [[255.23, 900.0]]}


def test_research_window_matches_the_library_match_window(store):
    assert store.candidates(900.0, mz_tol=0.01, ion_mode="Negative") == []
    assert fix_mass_tolerance(0.01, 900.0) == pytest.approx(0.018)
    # 既定のしきい値では窓に入るが |Δm/z| 15 mDa ≥ dmz_fail_mda（10 mDa）でハードに落ちる
    out = candidates.build_library_candidates(
        _ev(), msdial_matches=[], representative=None, store=store, scoring=SCORING,
        trends={"classes": {}}, th=judge.resolve_thresholds(None), exclude_current=False, top_n=5)
    assert out["candidates"] == [] and out["n_hard_removed"] == 1
    # しきい値を緩めれば、0.01 Da の窓の外の候補も再検索に出る
    loose = judge.resolve_thresholds({"dmz_fail_mda": 20.0})
    out = candidates.build_library_candidates(
        _ev(), msdial_matches=[], representative=None, store=store, scoring=SCORING,
        trends={"classes": {}}, th=loose, exclude_current=False, top_n=5)
    assert [c["name"] for c in out["candidates"]] == ["PC 40:6"]
    assert out["candidates"][0]["source"] == "research"
