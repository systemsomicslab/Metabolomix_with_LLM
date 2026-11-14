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


#: 照合の alignment の `measured` は測定の基準ピークを 100 とした強度（`match_spectrum` と同じ）
RELATIVE = {mz: intensity / 999.0 * 100.0 for mz, intensity in MEASURED}


def _result(score, matched_mz):
    alignment = [{"mz": mz, "measured": RELATIVE.get(mz, 50.0), "reference": 1.0, "matched": True}
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
    # 既定でも断片の一致を 1 本求め、基準ピークの 1% 未満の窓は数えない（2026-11-13 の事実確認 A1）
    assert rules["min_fragments_excl_precursor"] == 1
    assert rules["min_fragment_relative_intensity"] == 0.01
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
    # 既定の規則では順位と断片の本数だけ（質量差・クラス・イオンは設定したときだけ）
    out = M.qualify(_scored(), MEASURED, PRECURSOR, _rules(), ms2_tol=0.025)
    assert [c["check"] for c in out["checks"]] == ["target_in_top_rank", "min_fragments"]
    unset = M.qualify(_scored(), MEASURED, PRECURSOR, _rules(min_fragments_excl_precursor=None),
                      ms2_tol=0.025)
    assert [c["check"] for c in unset["checks"]] == ["target_in_top_rank"]


def test_default_rules_fail_a_precursor_only_match():
    # 既定（min_fragments 1）で、前駆体だけが合った候補は合格にしない
    out = M.qualify(_scored(), MEASURED, PRECURSOR,
                    M.normalize_rules({"target": {"inchikey": ISO}}), ms2_tol=0.025)
    check = next(c for c in out["checks"] if c["check"] == "min_fragments")
    assert check["status"] == "fail" and check["value"] == 0
    assert check["detail"]["min_relative_intensity"] == 0.01
    assert out["verdict"] == "not_qualified"


def test_noise_level_fragment_windows_are_not_counted():
    # 基準ピークの 0.5% の窓が「合った」だけでは断片の一致に数えない
    noisy = [(_record("PC 34:1", PC),
              {"total_score": 2.0, "weighted_dot_product": 0.6,
               "alignment": [{"mz": 300.0, "measured": 0.5, "reference": 80.0, "matched": True},
                             {"mz": PRECURSOR, "measured": 40.0, "reference": 100.0,
                              "matched": True}]})]
    out = M.qualify(noisy, MEASURED, PRECURSOR, _rules(), ms2_tol=0.025)
    check = next(c for c in out["checks"] if c["check"] == "min_fragments")
    assert check["status"] == "fail" and check["value"] == 0
    relaxed = M.qualify(noisy, MEASURED, PRECURSOR, _rules(min_fragment_relative_intensity=0),
                        ms2_tol=0.025)
    assert next(c for c in relaxed["checks"] if c["check"] == "min_fragments")["value"] == 1


@pytest.mark.parametrize("bad", [-0.1, 1.5, True, "0.01", float("nan")])
def test_min_fragment_relative_intensity_must_be_a_ratio(bad):
    with pytest.raises(DomainError) as info:
        _rules(min_fragment_relative_intensity=bad)
    assert info.value.details["field"] == "min_fragment_relative_intensity"


def test_top1_is_target_is_reported_next_to_target_rank():
    first = M.qualify(_scored(), MEASURED, PRECURSOR, _rules(), ms2_tol=0.025)
    assert first["target_rank"] == 1 and first["top1_is_target"] is True
    second = M.qualify(_scored(), MEASURED, PRECURSOR,
                       M.normalize_rules({"target": {"inchikey": ISO}}), ms2_tol=0.025)
    # 既定の top_rank 5 では 2 位でも順位の検査は pass する。1 位でないことは別の欄で読む
    assert second["checks"][0]["status"] == "pass"
    assert second["target_rank"] == 2 and second["top1_is_target"] is False
    absent = M.qualify(_scored(), MEASURED, PRECURSOR,
                       M.normalize_rules({"target": {"inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"}}),
                       ms2_tol=0.025)
    assert absent["target_rank"] is None and absent["top1_is_target"] is False


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


def test_neutral_loss_takes_min_relative_intensity():
    # 中性損失にも強度の下限を書ける（比は生成イオン＝前駆体 − delta の測定強度。E7）
    loss = M.normalize_ions([{"delta": 264.2451, "min_relative_intensity": 0.5}],
                            kind="loss", label="nl")
    assert loss[0]["min_relative_intensity"] == 0.5
    hit = M.find_ions(MEASURED, loss, precursor_mz=PRECURSOR)[0]
    assert hit["relative_intensity"] == pytest.approx(120.0 / 999.0) and hit["found"] is False
    # 省いた中性損失は従来の形のまま（rule_hash を変えない）
    assert "min_relative_intensity" not in M.normalize_ions([{"delta": 1.0}], kind="loss",
                                                            label="nl")[0]
    rules = _rules(neutral_losses=[{"delta": 264.2451, "min_relative_intensity": 0.05}])
    assert rules["neutral_losses"][0]["min_relative_intensity"] == 0.05


def test_ion_errors_name_the_offending_key():
    with pytest.raises(DomainError) as info:
        _rules(neutral_losses=[{"delta": 141.0, "intensity": 0.1}])
    assert "intensity" in info.value.message and "neutral_losses[0]" in info.value.message
    assert info.value.details == {"field": "neutral_losses[0]", "unknown": ["intensity"]}
    with pytest.raises(DomainError) as info:
        _rules(diagnostic_ions=[{"label": "PC"}])
    assert info.value.details["field"] == "diagnostic_ions[0].mz"
    # dcl_diagnostic_ions は中性損失を neutral_loss のキーで受ける（誤りの文もその名で）
    with pytest.raises(DomainError) as info:
        M.normalize_ion({"neutral_loss": 141.0, "minimum": 1}, kind="loss", label="ions[0]",
                        loss_key="neutral_loss")
    assert "neutral_loss" in info.value.message and "minimum" in info.value.message
    ok = M.normalize_ion({"neutral_loss": 141.0}, kind="loss", label="ions[0]",
                         loss_key="neutral_loss")
    assert ok["delta"] == 141.0


def test_rules_hash_ignores_key_order():
    a = _rules(top_rank=3, lipid_class="PC")
    b = M.normalize_rules({"lipid_class": "PC", "top_rank": 3, "target": {"inchikey": PC}})
    assert M.rules_hash(a) == M.rules_hash(b)
