"""msms-rules.v1 の種の目標（`target.lipid_species`）と、証拠の束の種の主張（純ロジック）。

InChIKey の規則・主張の hash は固定値で縛る。2026-11-14 に既定（`min_fragments_excl_precursor` 1・
`min_fragment_relative_intensity` 0.01）を変えたので、そのとき値を引き直した（意図した変化。事実確認 A1）。
陽イオンの既定の付加体の第 2 層を足したとき（A12）にも、付加体の段が hash に入る主張の値を引き直した。
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
    alignment = [{"mz": mz, "measured": 50.0, "reference": 1.0, "matched": True}
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
    {"lipid_species": "PC 15:0_18:1(d7)"},          # 内部標準（標識付きの記録は種に読まないので合否を出せない。レビュー V3）
    {"lipid_species": "PC 33:1(d7)"},
    {},
])
def test_invalid_species_targets(target):
    with pytest.raises(DomainError) as exc:
        M.normalize_rules({"target": target})
    assert exc.value.code == "MSMS_RULES_INVALID"


def test_inchikey_rules_hash_is_pinned():
    rules = M.normalize_rules({"target": {"inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N",
                                          "formula": "C42H82NO8P"},
                               "lipid_class": "PC", "same_formula_margin": {"min_delta": 0.1}})
    assert rules["target"] == {"inchikey_block": "JLPULHDHAOZNQI", "formula": "C42H82NO8P"}
    assert M.rules_hash(rules) == "46345c9781e7f5ea28cd95da198ff17f0424c3090fb3135bad987d6124700408"


def test_species_target_rank_is_found_by_the_record_name():
    out = M.qualify(_scored(), MEASURED, PRECURSOR, _species_rules(top_rank=1), ms2_tol=0.025)
    assert out["target_rank"] == 2
    assert out["checks"][0] == {"check": "target_in_top_rank", "status": "fail", "value": 2,
                                "detail": {"top_rank": 1, "species_route": "name"}}
    out = M.qualify(_scored(), MEASURED, PRECURSOR,
                    _species_rules(min_fragments_excl_precursor=2), ms2_tol=0.025)
    assert out["verdict"] == "qualified" and out["target_rank"] == 2


@pytest.mark.parametrize("record, species, route", [
    ({"name": "PC 16:0_18:1"}, "PC 34:1", "name"),
    ({"name": "low score: PC 34:1"}, "PC 34:1", "name"),
    # 注記・英語のクラス名・古い書式は直してから読む
    ({"name": "PC 34:1; ExampleDB ID-0000"}, "PC 34:1", "name_cleaned"),
    ({"name": "Phosphatidylcholine 16:0-18:1"}, "PC 34:1", "name_cleaned"),
    ({"name": "Cer[NS] d18:1/16:0"}, "Cer 34:1;O2", "name_cleaned"),
    # 酸化の修飾（;O ;O2 ;O3）は空白の注記が続いても残す（注記だけを切る。R10）
    ({"name": "PE 36:2;O (oxPE)"}, "PE 36:2;O", "name_cleaned"),
    ({"name": "FA 18:2;O (9-HODE)"}, "FA 18:2;O", "name_cleaned"),
    ({"name": "TG 52:2;O3 (x)"}, "TG 52:2;O3", "name_cleaned"),
    ({"name": "PE 36:2;O (oxPE) [M+H]+"}, "PE 36:2;O", "name_cleaned"),
    ({"name": "Cer 18:1;O2/16:0 (x)"}, "Cer 34:1;O2", "name_cleaned"),
    # `;O` 以外の修飾も注記の前で残す（レビュー V1）。和: C 18+24=42・二重結合 1・O 2+1=3、
    # エポキシの環は GOSLIN が二重結合 1 つに数える（20:4 の 4 + 環 1 = 和 36:5、O 1）
    ({"name": "Cer 18:1;O2/24:0;2OH (x)"}, "Cer 42:1;O3", "name_cleaned"),
    ({"name": "PC 16:0_20:4;Ep (x)"}, "PC 36:5;O", "name_cleaned"),
    ({"name": "PC 34:2;2OH (x)"}, "PC 34:2;O", "name_cleaned"),
    # GOSLIN が読めない修飾（`;(2OH)` `;(OH)`）は残して読まない（修飾を落とした種に読み違えない）
    ({"name": "Cer 18:1;O2/16:0;(2OH) (Cer AS)"}, None, None),
    ({"name": "PC 16:0_18:2;(OH) (x)"}, None, None),
    # 名前で読めなければ ClassyFire の Ontology と分子式から（LPS / LPG も）
    ({"name": "1-oleoyl-glycero-3-phosphoserine", "ontology": "1-acyl-glycerol-3-phosphoserines",
      "formula": "C24H46NO9P"}, "LPS 18:1", "ontology_formula"),
    ({"name": "lysophosphatidylglycerol", "ontology": "Lysophosphatidylglycerols",
      "formula": "C22H45O9P"}, "LPG 16:0", "ontology_formula"),
    ({"name": "LPS 18:1"}, "LPS 18:1", "name"),
    # 同位体の標識（内部標準）・脂質でない記録は種に読まない
    ({"name": "PC 15:0-18:1(d7)"}, None, None),
    ({"name": "Glucose", "ontology": "Hexoses", "formula": "C6H12O6"}, None, None),
    ({"name": ""}, None, None),
])
def test_record_lipid_species_routes(record, species, route):
    assert M.record_lipid_species(record) == (species, route)


OXPE_FORMULA = "C41H78NO9P"   # PE 36:2;O（PE 36:2 の C41H78NO8P に O が 1 つ）


def test_oxidized_species_target_is_found_through_a_record_with_a_note():
    """`PE 36:2;O (oxPE)` の記録は種 PE 36:2;O（修飾を落とした PE 36:2 ではない。R10）。

    目標の候補になり、同じ種の別の分子種（`PE 18:0_18:2;O`）は `same_formula_margin` の別物に数えない。"""
    scored = [({"name": "PE 36:2;O (oxPE)", "formula": OXPE_FORMULA}, _result(2.6, [754.5381])),
              ({"name": "PE 18:0_18:2;O", "formula": OXPE_FORMULA}, _result(2.4, [754.5381]))]
    rules = M.normalize_rules({"target": {"lipid_species": "PE 36:2;O"},
                               "same_formula_margin": {"min_delta": 0.1}})
    out = M.qualify(scored, MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    assert out["target_rank"] == 1 and out["top1_is_target"] is True
    assert out["checks"][0]["detail"]["species_route"] == "name_cleaned"
    margin = next(c for c in out["checks"] if c["check"] == "same_formula_margin")
    assert margin["status"] == "pass" and margin["detail"]["reason"] == "no_same_formula_alternative"
    # 修飾の無い PE 36:2 は別の種（同じ記録を目標に数えない）
    plain = M.qualify(scored, MEASURED, PRECURSOR,
                      M.normalize_rules({"target": {"lipid_species": "PE 36:2"}}), ms2_tol=0.025)
    assert plain["target_rank"] is None


def test_msdial_agreement_reads_the_oxidized_record_as_the_same_compound():
    from metabolomix.library.matching import msdial_agreement

    scored = [({"name": "PE 36:2;O (oxPE)"}, {}), ({"name": "PE 36:2"}, {})]
    assert msdial_agreement(scored, {"name": "PE 36:2;O"}) == {
        "msdial_name": "PE 36:2;O", "msdial_rank": 1, "top_differs_from_msdial": False}
    assert msdial_agreement(scored, {"name": "PE 36:2"})["msdial_rank"] == 2


CER_NS_FORMULA = "C34H67NO3"   # Cer 18:1;O2/16:0（Cer 34:1;O2）
CER_AS_FORMULA = "C34H67NO4"   # Cer 18:1;O2/16:0;(2OH)（脂肪酸に水酸基。O が 1 つ多い）


def test_a_hydroxy_fatty_acyl_record_with_a_note_is_not_the_non_hydroxy_species():
    """`Cer 18:1;O2/16:0;(2OH) (Cer AS)` の記録を、注記を切るときに `;(2OH)` まで落として
    Cer NS の種 Cer 34:1;O2 に読まない（レビュー V1）。目標の候補にも MS-DIAL の注釈と同じ物質にも数えない。"""
    from metabolomix.library.matching import msdial_agreement

    scored = [({"name": "Cer 18:1;O2/16:0;(2OH) (Cer AS)", "formula": CER_AS_FORMULA}, _result(2.6, [754.5381])),
              ({"name": "Cer 18:1;O2/16:0", "formula": CER_NS_FORMULA}, _result(2.4, [754.5381]))]
    out = M.qualify(scored, MEASURED, PRECURSOR,
                    M.normalize_rules({"target": {"lipid_species": "Cer 34:1;O2"}}), ms2_tol=0.025)
    assert out["target_rank"] == 2 and out["top1_is_target"] is False
    assert msdial_agreement(scored, {"name": "Cer 34:1;O2"}) == {
        "msdial_name": "Cer 34:1;O2", "msdial_rank": 2, "top_differs_from_msdial": True}


@pytest.mark.parametrize("name", ["Cer 18:1;O2/24:0;2OH", "PC 16:0_20:4;Ep", "PC 34:2;2OH", "FA 18:2;9OH,13OH",
                                  "FA 18:1;oxo", "TG 16:0_18:1_18:2;O(FA 16:0)", "Cer 18:1;O2/32:0;O(FA 18:2)"])
def test_a_note_does_not_change_the_record_species(name):
    bare = M.record_lipid_species({"name": name})
    assert bare[1] == "name" and bare[0] is not None
    assert M.record_lipid_species({"name": name + " (x)"}) == (bare[0], "name_cleaned")


def test_species_target_found_through_the_ontology_route():
    scored = [(_record("Glucose", formula="C6H12O6"), _result(2.0, [754.5381])),
              ({"name": "1-oleoyl-glycero-3-phosphoserine",
                "ontology": "1-acyl-glycerol-3-phosphoserines", "formula": "C24H46NO9P"},
               _result(1.5, [184.0733, 754.5381]))]
    rules = M.normalize_rules({"target": {"lipid_species": "LPS 18:1"}})
    out = M.qualify(scored, MEASURED, PRECURSOR, rules, ms2_tol=0.025)
    assert out["target_rank"] == 2 and out["top1_is_target"] is False
    assert out["checks"][0]["detail"]["species_route"] == "ontology_formula"


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


def test_inchikey_claim_rule_hash_is_pinned():
    spec = {"claim_id": "C1", "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0},
            "scan": {"ppm": [10, 20], "rt_tol_min": [0.3]},
            "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}},
            "relations": {"rt_window": 0.05}}
    # merge_split=false は v7 と同じ規則。`msms` があると MS/MS の規則の既定（2026-11-14 に変えた）が
    # hash に入るので、v7 の値（fb5977…）から変わった。既定（true）は merge_split が規則に入る。
    # 陽イオンの既定の付加体の第 2 層に [M+CH3OH+H]+ と [M+H-2H2O]+ を足したとき（A12）にも段の並びが
    # 変わったので値を書き直した（旧い段では 9be780…・c17e3e… が再現することを確かめた）
    normalized = C.normalize_claim_spec({**spec, "merge_split": False}, dataset_polarity="positive")
    assert C.rule_hash(normalized) == \
        "518eae70c026ef5bf94fd511bfa7127d281d5cc8798400b954f63cc206b28e7e"
    # `msms` が無ければ MS/MS の既定は hash に入らない（付加体の段の変更だけを受ける）
    no_msms = {k: v for k, v in spec.items() if k != "msms"}
    assert C.rule_hash(C.normalize_claim_spec({**no_msms, "merge_split": False},
                                              dataset_polarity="positive")) == \
        "0c2297849bae893d85d724a03a45b55a062792a0492d05365d1618938e63d7fb"
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


def test_isotope_labelled_target_names_the_field():
    with pytest.raises(DomainError) as exc:
        M.normalize_rules({"target": {"lipid_species": "PC 15:0_18:1(d7)"}})
    assert exc.value.details == {"field": "target.lipid_species"}
    assert "標識" in exc.value.message
    assert M.isotope_labelled("PC 15:0_18:1(d7)") and not M.isotope_labelled("SM d18:1/16:0")
