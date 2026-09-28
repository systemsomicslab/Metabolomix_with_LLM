import pytest

from lipidmix.curation import judge

TH = judge.resolve_thresholds(None)


def match(**kw):
    base = {"has_msms": True, "is_reference_matched": True, "is_annotation_suggested": False,
            "is_precursor_mz_match": True, "is_spectrum_match": True,
            "is_manually_modified": False, "squared_weighted_dot_product": 0.81}
    base.update(kw)
    return base


def ev(**kw):
    base = {"name_prefix": None, "match": match(), "reference": {"rt": 12.0},
            "ppm": 1.5, "adduct_band": "PASS", "drt": 0.1, "rescore": {"weighted_dot_product": 0.9},
            "eic_shape": {"band": "PASS", "apex_rt_sd": 0.01}, "notes": [],
            "adduct": "[M+H]+", "reference_adduct": "[M+H]+"}
    base.update(kw)
    return base


def test_clean_spot_is_ok():
    result = judge.judge_spot(ev(), {"outlier": False, "reliable": True}, TH)
    assert result["verdict"] == "ok"
    assert result["reasons"] == []


def test_no_ms2_is_unknown_not_fail():
    result = judge.judge_spot(ev(name_prefix="no MS2", match=match(has_msms=False,
                              is_reference_matched=False, is_annotation_suggested=True)), None, TH)
    assert result["checks"]["msms"]["band"] == "UNKNOWN"
    assert result["verdict"] == "ok"


def test_low_score_is_suspect():
    result = judge.judge_spot(ev(name_prefix="low score", match=match(is_reference_matched=False,
                              is_annotation_suggested=True)), None, TH)
    assert result["checks"]["msms"]["band"] == "FAIL"
    assert result["verdict"] == "suspect"
    assert "low_score" in result["reasons"]


def test_large_ppm_is_likely_wrong():
    result = judge.judge_spot(ev(ppm=25.0), None, TH)
    assert result["verdict"] == "likely_wrong"
    assert result["reasons"][0] == "ppm_out"


def test_polarity_mismatch_is_likely_wrong():
    assert judge.judge_spot(ev(adduct_band="FAIL"), None, TH)["verdict"] == "likely_wrong"


def test_precursor_unmatched_flag_from_msdial_is_strong():
    result = judge.judge_spot(ev(match=match(is_precursor_mz_match=False)), None, TH)
    assert result["verdict"] == "likely_wrong"


def test_trend_alone_never_raises_the_verdict():
    result = judge.judge_spot(ev(), {"outlier": True, "reliable": True, "z": 5.0}, TH)
    assert result["checks"]["trend"]["band"] == "BORDERLINE"
    assert result["verdict"] == "ok"


def test_trend_reinforces_a_single_borderline():
    result = judge.judge_spot(ev(ppm=7.0), {"outlier": True, "reliable": True, "z": 5.0}, TH)
    assert result["verdict"] == "suspect"


def test_unreliable_trend_is_information_only():
    result = judge.judge_spot(ev(ppm=7.0), {"outlier": True, "reliable": False, "z": 5.0}, TH)
    assert result["checks"]["trend"]["band"] == "PASS"
    assert "trend_outlier_unreliable" in result["info"]
    assert result["verdict"] == "ok"


def test_missing_reference_makes_rt_unknown():
    result = judge.judge_spot(ev(reference=None, drt=None, rescore=None), None, TH)
    assert result["checks"]["rt"]["band"] == "UNKNOWN"
    assert "reference_not_found" in result["info"]


def test_poor_eic_is_suspect():
    result = judge.judge_spot(ev(eic_shape={"band": "FAIL", "apex_rt_sd": 0.01}), None, TH)
    assert result["verdict"] == "suspect"
    assert "eic_poor" in result["reasons"]


def test_rescore_discrepancy_is_information():
    result = judge.judge_spot(ev(rescore={"weighted_dot_product": 0.2}), None, TH)
    assert "rescore_discrepancy" in result["info"]


def test_adduct_differing_from_reference_is_information():
    result = judge.judge_spot(ev(reference_adduct="[M+Na]+"), None, TH)
    assert "adduct_differs_from_reference" in result["info"]
    assert result["verdict"] == "ok"


def test_unknown_threshold_key_is_rejected():
    with pytest.raises(ValueError):
        judge.resolve_thresholds({"ppm_pas": 3})
