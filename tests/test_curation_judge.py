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


def test_large_ppm_alone_is_suspect_not_likely_wrong():
    # ppm > 10 は adduct 非依存(median 約 -0.8 ppm、全 adduct 共通)で、実データで
    # likely_wrong が ppm_out 単独からしか出ていなかった。弱い理由に格下げ(2026-09-29 ユーザー決定)。
    result = judge.judge_spot(ev(ppm=25.0), None, TH)
    assert result["verdict"] == "suspect"
    assert result["reasons"][0] == "ppm_out"
    assert result["checks"]["mz"]["band"] == "FAIL"


def test_large_ppm_with_precursor_unmatched_is_likely_wrong():
    result = judge.judge_spot(ev(ppm=25.0, match=match(is_precursor_mz_match=False)), None, TH)
    assert result["verdict"] == "likely_wrong"
    assert "precursor_unmatched" in result["reasons"]
    assert "ppm_out" in result["reasons"]


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


# --- MS-DIAL の脂質規則フラグ（KB facts/msdial-lcms-lipidomics-match-flags-semantics、
# ユーザー承認 2026-09-29）。規則フラグは Lipidomics 採点器でしか立たないので、
# lipid_rules=True（そのレビューで規則が走った証拠がある）ときだけ効く。---

def rejected(**kw):
    return match(is_reference_matched=False, is_annotation_suggested=True,
                 is_lipid_class_match=False, is_lipid_chains_match=False,
                 is_other_lipid_match=False, **kw)


def test_class_rule_rejected_is_likely_wrong():
    result = judge.judge_spot(ev(name_prefix="low score", match=rejected()), None, TH,
                              lipid_rules=True)
    assert result["verdict"] == "likely_wrong"
    assert result["reasons"][0] == "class_rule_rejected"
    assert "low_score" in result["reasons"]


def test_class_rule_rejected_is_ignored_without_lipid_rules():
    # メタボロミクス採点器では規則フラグが全部 False になる——low score が全件 likely_wrong に化けない。
    result = judge.judge_spot(ev(name_prefix="low score", match=rejected()), None, TH)
    assert result["verdict"] == "suspect"
    assert "class_rule_rejected" not in result["reasons"]


def test_class_rule_rejected_needs_msms():
    # MS/MS 無しなら規則フラグは全部 False（規則を評価した結果ではない）。
    result = judge.judge_spot(ev(name_prefix="no MS2", match=rejected(has_msms=False)), None, TH,
                              lipid_rules=True)
    assert "class_rule_rejected" not in result["reasons"]
    assert result["verdict"] == "ok"


def test_rules_not_run_is_info_only():
    m = match(is_lipid_class_match=False, is_lipid_chains_match=False, is_other_lipid_match=True)
    result = judge.judge_spot(ev(match=m), None, TH, lipid_rules=True)
    assert result["verdict"] == "ok"
    assert "class_rules_not_run" in result["info"]


def test_chain_level_name_without_chain_support_is_info():
    m = match(is_lipid_class_match=True, is_lipid_chains_match=False, is_other_lipid_match=False)
    result = judge.judge_spot(ev(name="PC 34:1|PC 16:0_18:1", match=m), None, TH, lipid_rules=True)
    assert result["verdict"] == "ok"
    assert "chains_unsupported" in result["info"]


@pytest.mark.parametrize("name,chains", [("PC 34:1", False), ("LPC 16:0", False),
                                         ("PC 16:0_18:1", True), ("Cer 18:1;O2/16:0", True)])
def test_chains_unsupported_only_for_multi_chain_names_without_support(name, chains):
    m = match(is_lipid_class_match=True, is_lipid_chains_match=chains, is_other_lipid_match=False)
    result = judge.judge_spot(ev(name=name, match=m), None, TH, lipid_rules=True)
    assert "chains_unsupported" not in result["info"]


def test_chain_level_name_is_detected():
    assert judge.is_chain_level_name("PC 16:0_18:1")
    assert judge.is_chain_level_name("low score: PC 34:1|PC 16:0_18:1")
    assert judge.is_chain_level_name("Cer 18:1;O2/16:0")
    assert not judge.is_chain_level_name("PC 34:1")
    assert not judge.is_chain_level_name("LPC 16:0")
    assert not judge.is_chain_level_name(None)


def test_lipid_rule_infos_need_lipid_rules():
    m = match(is_lipid_class_match=False, is_lipid_chains_match=False, is_other_lipid_match=True)
    result = judge.judge_spot(ev(name="PC 16:0_18:1", match=m), None, TH)
    assert "class_rules_not_run" not in result["info"]
    assert "chains_unsupported" not in result["info"]


def test_lipid_rules_active_detects_any_rule_flag():
    assert judge.lipid_rules_active([{"match": match(is_lipid_class_match=True)}, {"match": None}])
    assert judge.lipid_rules_active([{"match": match(is_other_lipid_match=True)}])
    assert not judge.lipid_rules_active([{"match": rejected()}, {"match": None}])
    assert not judge.lipid_rules_active([])
