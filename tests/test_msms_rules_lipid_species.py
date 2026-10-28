"""msms-rules.v1 の種の目標（`target.lipid_species`）と、証拠の束の種の主張（純ロジック）。

InChIKey の規則・主張の hash は、この口を足す前の値で固定する（既存の挙動を変えない）。
"""
from __future__ import annotations

import pytest

from metabolomix.analysis import claim_evidence as C
from metabolomix.analysis import msms_rules as M
from metabolomix.core.atomic_io import DomainError

PC_FORMULA = "C42H76NO8P"   # PC 34:4 と PE 37:4 は同じ分子式
PRECURSOR = 754.5381
MEASURED = [[184.0733, 999.0], [478.33, 80.0], [754.5381, 300.0]]


def _record(name, inchikey=None, formula=PC_FORMULA):
    return {"name": name, "inchikey": inchikey, "formula": formula}


def _result(score, matched_mz):
    alignment = [{"mz": mz, "measured": 1.0, "reference": 1.0, "matched": True}
                 for mz in matched_mz]
    return {"total_score": score, "weighted_dot_product": score / 3, "alignment": alignment}


def _scored():
    return [
        (_record("PE 37:4"), _result(2.6, [754.5381])),
        (_record("PC 16:1_18:3"), _result(2.4, [184.0733, 478.33, 754.5381])),
        (_record("PC O-34:4"), _result(1.0, [184.0733])),
    ]


def _species_rules(**extra):
    return M.normalize_rules({"target": {"lipid_species": "PC 18:3_16:1"}, **extra})


def test_species_target_is_normalized_to_species_level():
    rules = _species_rules()
    assert rules["target"] == {"lipid_species": "PC 34:4", "formula": None}


@pytest.mark.parametrize("target", [
    {"lipid_species": "PC 34:4", "inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N"},
    {"lipid_species": "PC"},
    {"lipid_species": "Glucose"},
    {},
])
def test_invalid_species_targets(target):
    with pytest.raises(DomainError) as exc:
        M.normalize_rules({"target": target})
    assert exc.value.code == "MSMS_RULES_INVALID"


def test_inchikey_rules_are_unchanged():
    rules = M.normalize_rules({"target": {"inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N",
                                          "formula": "C42H82NO8P"},
                               "lipid_class": "PC", "same_formula_margin": {"min_delta": 0.1}})
    assert rules["target"] == {"inchikey_block": "JLPULHDHAOZNQI", "formula": "C42H82NO8P"}
    assert M.rules_hash(rules) == "932ae02d1c890c8bfcb87eec9b1345bb8cd99a8a864b45f1ef0f08116d132e69"


def test_species_target_rank_is_found_by_the_record_name():
    out = M.qualify(_scored(), MEASURED, PRECURSOR, _species_rules(top_rank=1), ms2_tol=0.025)
    assert out["target_rank"] == 2
    assert out["checks"][0] == {"check": "target_in_top_rank", "status": "fail", "value": 2,
                                "detail": {"top_rank": 1}}
    out = M.qualify(_scored(), MEASURED, PRECURSOR,
                    _species_rules(min_fragments_excl_precursor=2), ms2_tol=0.025)
    assert out["verdict"] == "qualified" and out["target_rank"] == 2


def test_ether_record_is_not_the_diacyl_species():
    rules = M.normalize_rules({"target": {"lipid_species": "PC O-34:4"}})
    assert M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)["target_rank"] == 3


def test_same_formula_margin_compares_against_other_species():
    rules = _species_rules(same_formula_margin={"min_delta": 0.0})
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    check = next(c for c in out["checks"] if c["check"] == "same_formula_margin")
    # 目標の分子式は記録から（PC 16:1_18:3 の C42H76NO8P）。PE 37:4 が同じ分子式で上にいる
    assert check["status"] == "fail" and check["value"] == pytest.approx(-0.2)
    assert check["detail"]["best_alternative"] == "PE 37:4"
    # 別の分子種の同じ種（PC 14:0_20:4）は「別物」ではない
    same_species = [(_record("PC 14:0_20:4"), _result(3.0, [754.5381]))] + _scored()[1:2]
    out = M.qualify(same_species, MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    check = next(c for c in out["checks"] if c["check"] == "same_formula_margin")
    assert check["status"] == "pass" and check["detail"]["reason"] == "no_same_formula_alternative"


# ---------- 証拠の束の spec ----------

GABA_KEY = "BTCSSZJGUNDROE-UHFFFAOYSA-N"


def test_inchikey_claim_rule_hash_is_unchanged():
    spec = {"claim_id": "C1", "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0},
            "scan": {"ppm": [10, 20], "rt_tol_min": [0.3]},
            "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}},
            "relations": {"rt_window": 0.05}}
    # merge_split=false（v7 と同じ規則）は v7 の hash のまま。既定（true）は merge_split が規則に入る
    normalized = C.normalize_claim_spec({**spec, "merge_split": False}, dataset_polarity="positive")
    assert C.rule_hash(normalized) == \
        "fb597799f7388814f93093b472ee773769e40aabaaf6ad0bb0fabb56338f33db"
    merged = C.normalize_claim_spec(spec, dataset_polarity="positive")
    assert merged["merge_split"] is True
    assert C.rule_hash(merged) != C.rule_hash(normalized)


def test_species_claim_passes_the_species_to_the_msms_rules():
    spec = {"claim_id": "C1",
            "target": {"lipid_species": "PC 18:3_16:1", "formula": PC_FORMULA, "rt_min": 8.4},
            "msms": {"rules": {"lipid_class": "PC"}}}
    normalized = C.normalize_claim_spec(spec, dataset_polarity="positive")
    assert normalized["target"]["lipid_species"] == "PC 34:4"
    assert normalized["msms"]["rules"]["target"] == {"lipid_species": "PC 34:4",
                                                     "formula": PC_FORMULA}
    same = C.normalize_claim_spec(
        {**spec, "target": {**spec["target"], "lipid_species": "PC 34:4"}},
        dataset_polarity="positive")
    assert C.rule_hash(normalized) == C.rule_hash(same)
