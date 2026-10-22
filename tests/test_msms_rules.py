"""MS/MS の合否・診断イオン・中性損失（spec 2026-10-22 §4 の 6）。候補とスペクトルは手で作る。"""
from __future__ import annotations

import pytest

from metabolomix.analysis import msms_rules as M
from metabolomix.core.atomic_io import DomainError

PC = "JLPULHDHAOZNQI-ZTIMHPMXSA-N"
ISO = "WTJKGGKOPKCXLL-VYOBOKEXSA-N"
PRECURSOR = 760.5851
MEASURED = [[184.0733, 999.0], [496.34, 120.0], [760.5851, 400.0]]


def _record(name, inchikey, formula="C42H82NO8P"):
    return {"name": name, "inchikey": inchikey, "formula": formula}


def _result(score, matched_mz):
    alignment = [{"mz": mz, "measured": 1.0, "reference": 1.0, "matched": True}
                 for mz in matched_mz]
    return {"total_score": score, "weighted_dot_product": score / 3, "alignment": alignment}


def _scored():
    return [
        (_record("PC 34:1", PC), _result(2.5, [184.0733, 496.34, 760.5851])),
        (_record("PC 16:0_18:1 isomer", ISO), _result(2.0, [760.5851])),
    ]


def _rules(**extra):
    return M.normalize_rules({"target": {"inchikey": PC}, **extra})


def test_defaults_and_unknown_key():
    rules = _rules()
    assert rules["schema"] == "msms-rules.v1"
    assert rules["target"] == {"inchikey_block": "JLPULHDHAOZNQI", "formula": None}
    assert rules["top_rank"] == 5 and rules["diagnostic_ions"] == []
    with pytest.raises(DomainError):
        M.normalize_rules({"target": {"inchikey": PC}, "magic": 1})
    with pytest.raises(DomainError):
        M.normalize_rules({"target": {"inchikey": "not-a-key"}})


def test_qualified_when_all_configured_checks_pass():
    rules = _rules(min_fragments_excl_precursor=2, same_formula_margin={"min_delta": 0.1},
                   lipid_class="PC",
                   diagnostic_ions=[{"mz": 184.0733, "label": "phosphocholine"}],
                   neutral_losses=[{"delta": 264.245, "label": "-264"}])
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    assert out["verdict"] == "qualified" and out["target_rank"] == 1
    checks = {c["check"]: c for c in out["checks"]}
    assert checks["min_fragments"]["value"] == 2
    assert checks["same_formula_margin"]["value"] == pytest.approx(0.5)
    assert checks["lipid_class"]["value"] == "PC"


def test_precursor_only_match_fails_min_fragments():
    rules = M.normalize_rules({"target": {"inchikey": ISO}, "min_fragments_excl_precursor": 1})
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    checks = {c["check"]: c for c in out["checks"]}
    assert checks["min_fragments"]["status"] == "fail" and checks["min_fragments"]["value"] == 0
    assert out["verdict"] == "not_qualified"


def test_same_formula_alternative_above_target_fails_margin():
    rules = M.normalize_rules({"target": {"inchikey": ISO},
                               "same_formula_margin": {"min_delta": 0.0}})
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    check = next(c for c in out["checks"] if c["check"] == "same_formula_margin")
    assert check["status"] == "fail" and check["value"] == pytest.approx(-0.5)
    assert check["detail"]["best_alternative"] == "PC 34:1"


def test_target_outside_top_rank():
    rules = M.normalize_rules({"target": {"inchikey": ISO}, "top_rank": 1})
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    assert out["checks"][0] == {"check": "target_in_top_rank", "status": "fail", "value": 2,
                                "detail": {"top_rank": 1}}


def test_unconfigured_checks_are_not_reported():
    out = M.qualify(_scored(), MEASURED, PRECURSOR, _rules(), ms2_tol=0.025)
    assert [c["check"] for c in out["checks"]] == ["target_in_top_rank"]


def test_target_absent_makes_dependent_checks_not_evaluable():
    rules = M.normalize_rules({"target": {"inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
                               "min_fragments_excl_precursor": 1})
    out = M.qualify(_scored(), MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    assert out["verdict"] == "not_qualified"            # target_in_top_rank が fail
    assert out["checks"][1]["status"] == "not_evaluable"


def test_find_ions_relative_intensity_and_neutral_loss():
    ions = M.normalize_ions([{"mz": 184.07, "label": "PC", "min_relative_intensity": 0.5},
                             {"mz": 496.34, "label": "LPC", "min_relative_intensity": 0.5}],
                            kind="ion", label="ions")
    found = M.find_ions(MEASURED, ions, precursor_mz=PRECURSOR)
    assert [f["found"] for f in found] == [True, False]          # 496 は 0.12 しかない
    assert found[0]["relative_intensity"] == pytest.approx(1.0)
    loss = M.normalize_ions([{"delta": 576.5118, "label": "x"}], kind="loss", label="nl")
    hit = M.find_ions(MEASURED, loss, precursor_mz=PRECURSOR)[0]
    assert hit["kind"] == "neutral_loss" and hit["found"] and hit["target_mz"] == pytest.approx(184.0733)


def test_rules_hash_ignores_key_order():
    a = _rules(top_rank=3, lipid_class="PC")
    b = M.normalize_rules({"lipid_class": "PC", "top_rank": 3, "target": {"inchikey": PC}})
    assert M.rules_hash(a) == M.rules_hash(b)
