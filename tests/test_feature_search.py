# tests/test_feature_search.py
"""feature_search の土台（目標の正規化・理論 m/z・規則 hash）のテスト。"""
import pytest

from metabolomix.analysis.feature_search import (
    DEFAULT_ADDUCT_TIERS,
    MAX_TARGETS,
    normalize_targets,
    rule_hash,
    theoretical_mz,
)
from metabolomix.core.atomic_io import DomainError

GABA_KEY = "BTCSSZJGUNDROE-UHFFFAOYSA-N"


def _gaba(**extra):
    return {"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", **extra}


def test_theoretical_mz_of_gaba_protonated():
    assert theoretical_mz("C4H9NO2", "[M+H]+") == pytest.approx(104.0706, abs=1e-4)


def test_theoretical_mz_unreadable_returns_none():
    assert theoretical_mz("", "[M+H]+") is None
    assert theoretical_mz("Xx9", "[M+H]+") is None
    # parse_formula は読めない部分を黙って捨てるので、全体が読めることを先に確かめる
    for junk in ("C4H9NO2!!", "Ca(OH)2", "C4H9NO2-3", "C4H9NO2 junk", "c4h9no2"):
        assert theoretical_mz(junk, "[M+H]+") is None, junk
    assert theoretical_mz("C4H9NO2", "[bogus]") is None


def test_default_tiers_follow_dataset_polarity():
    [t] = normalize_targets([_gaba()], dataset_polarity="negative")
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["negative"]
    assert t["inchikey_block"] == "BTCSSZJGUNDROE" and t["ppm"] == 10 and t["rt_tol_min"] == 0.3
    assert t["polarity"] == "negative" and t["applicable"] is True
    # 理論 m/z は全層の全付加体について計算される
    flat = [a for tier in DEFAULT_ADDUCT_TIERS["negative"] for a in tier]
    assert set(t["theoretical_mz"]) == set(flat)
    assert t["theoretical_mz"]["[M-H]-"] == pytest.approx(102.0561, abs=1e-4)


def test_default_tiers_values_are_the_spec_ones():
    assert DEFAULT_ADDUCT_TIERS["positive"] == [
        ["[M+H]+", "[M]+"],
        ["[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M+H-NH3]+", "[2M+H]+", "[M+2H]2+"],
    ]
    assert DEFAULT_ADDUCT_TIERS["negative"] == [
        ["[M-H]-"],
        ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[2M-H]-", "[M+Na-2H]-"],
    ]
    assert MAX_TARGETS == 50


def test_normalized_tiers_are_copies_of_defaults():
    [t] = normalize_targets([_gaba()], dataset_polarity="positive")
    t["adduct_tiers"][0].append("[M+K]+")
    assert DEFAULT_ADDUCT_TIERS["positive"][0] == ["[M+H]+", "[M]+"]


def test_expected_mz_form_fixes_a_single_adduct():
    [t] = normalize_targets(
        [{"target_id": "c1", "inchikey": GABA_KEY.lower(), "expected_mz": 104.0706,
          "adduct": "[M+H]+", "rt_min": 3.2, "name": "GABA"}],
        dataset_polarity="positive")
    assert t["adduct_tiers"] == [["[M+H]+"]]
    assert t["theoretical_mz"] == {"[M+H]+": 104.0706}
    assert t["formula"] is None and t["expected_mz"] == 104.0706
    assert t["inchikey_block"] == "BTCSSZJGUNDROE"
    assert t["rt_min"] == 3.2 and t["name"] == "GABA"


def test_explicit_adduct_tiers_are_kept_and_computed():
    [t] = normalize_targets([_gaba(adduct_tiers=[["[M+H]+"], ["[M+Na]+"]])], dataset_polarity="positive")
    assert t["adduct_tiers"] == [["[M+H]+"], ["[M+Na]+"]]
    assert set(t["theoretical_mz"]) == {"[M+H]+", "[M+Na]+"}


def test_rule_hash_is_the_same_whether_defaults_are_explicit():
    omitted = normalize_targets([_gaba()], dataset_polarity="positive")
    explicit = normalize_targets(
        [_gaba(ppm=10, rt_tol_min=0.3, adduct_tiers=DEFAULT_ADDUCT_TIERS["positive"],
               polarity="positive")],
        dataset_polarity="positive")
    assert rule_hash(omitted) == rule_hash(explicit)


def test_rule_hash_ignores_int_float_and_list_tuple_spelling():
    omitted = normalize_targets([_gaba()], dataset_polarity="positive")
    spelled = normalize_targets(
        [_gaba(ppm=10.0, rt_tol_min=0.3,
               adduct_tiers=tuple(tuple(t) for t in DEFAULT_ADDUCT_TIERS["positive"]))],
        dataset_polarity="positive")
    assert rule_hash(omitted) == rule_hash(spelled)
    assert isinstance(spelled[0]["ppm"], float) and isinstance(spelled[0]["adduct_tiers"][0], list)
    a = normalize_targets([_gaba(rt_min=3)], dataset_polarity="positive")
    b = normalize_targets([_gaba(rt_min=3.0)], dataset_polarity="positive")
    assert rule_hash(a) == rule_hash(b)
    mz = {"target_id": "c1", "inchikey": GABA_KEY, "adduct": "[M+H]+"}
    assert (rule_hash(normalize_targets([{**mz, "expected_mz": 104}], dataset_polarity="positive"))
            == rule_hash(normalize_targets([{**mz, "expected_mz": 104.0}], dataset_polarity="positive")))


def test_rule_hash_changes_with_the_rule():
    base = rule_hash(normalize_targets([_gaba()], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba(ppm=5)], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba(rt_min=3.0)], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba()], dataset_polarity="negative"))


def test_target_of_other_polarity_is_not_applicable():
    [t] = normalize_targets([_gaba(polarity="negative")], dataset_polarity="positive")
    assert t["applicable"] is False
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["negative"]
    [t2] = normalize_targets([_gaba(polarity="positive")], dataset_polarity="positive")
    assert t2["applicable"] is True


@pytest.mark.parametrize("bad, field", [
    ({"target_id": "c1", "formula": "C4H9NO2"}, "inchikey"),
    ({"target_id": "c1", "inchikey": GABA_KEY}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07}, "adduct"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "adduct": "[M+H]+"}, "expected_mz"),
    ({"target_id": "c1", "inchikey": "short", "formula": "C4H9NO2"}, "inchikey"),
    ({"target_id": "", "inchikey": GABA_KEY, "formula": "C4H9NO2"}, "target_id"),
    ({"inchikey": GABA_KEY, "formula": "C4H9NO2"}, "target_id"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "bogus": 1}, "bogus"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "Xx9"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2!!"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "Ca(OH)2"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2-3"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "ppm": 0}, "ppm"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "ppm": True}, "ppm"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_tol_min": -1}, "rt_tol_min"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": "3"}, "rt_min"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "polarity": "pos"}, "polarity"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[nope]"}, "adduct"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "expected_mz": 104.07,
      "adduct": "[M+H]+"}, "expected_mz"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[M+H]+",
      "adduct_tiers": [["[M+H]+"]]}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
      "adduct_tiers": [["[M+Zz]+"]]}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "adduct_tiers": []}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
      "adduct_tiers": "[M+H]+"}, "adduct_tiers"),
])
def test_invalid_targets_name_the_field(bad, field):
    with pytest.raises(DomainError) as exc:
        normalize_targets([bad], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == field
    assert "target_id" in exc.value.details


def test_adduct_tiers_contradicting_polarity_is_invalid():
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(adduct_tiers=[["[M-H]-"]])], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == "adduct_tiers"
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(polarity="negative", adduct_tiers=[["[M+H]+"]])],
                          dataset_polarity="negative")
    assert exc.value.details["field"] == "adduct_tiers"


def test_duplicate_target_id_is_invalid():
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(), _gaba()], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["target_id"] == "c1" and exc.value.details["field"] == "target_id"


def test_more_than_fifty_targets_is_invalid():
    many = [{"target_id": f"c{i}", "inchikey": GABA_KEY, "formula": "C4H9NO2"} for i in range(51)]
    with pytest.raises(DomainError) as exc:
        normalize_targets(many, dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == "targets"
    assert len(normalize_targets(many[:50], dataset_polarity="positive")) == 50


def test_unknown_polarity_without_explicit_adducts_is_refused():
    explicit = {"target_id": "c2", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[M+H]+"}
    tiers = _gaba(target_id="c3", adduct_tiers=[["[M+H]+"]])
    bare1 = _gaba(target_id="c4")
    bare2 = _gaba(target_id="c5")
    with pytest.raises(DomainError) as exc:
        normalize_targets([explicit, tiers, bare1, bare2], dataset_polarity=None)
    assert exc.value.code == "POLARITY_UNKNOWN"
    assert exc.value.details["target_ids"] == ["c4", "c5"]
    # 明示の付加体だけなら極性不明でも通る
    out = normalize_targets([explicit, tiers], dataset_polarity=None)
    assert [t["applicable"] for t in out] == [True, True]
    # 目標自身が極性を持てば、データセットの極性が不明でも既定の層が決まる
    [t] = normalize_targets([_gaba(polarity="positive")], dataset_polarity=None)
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["positive"] and t["applicable"] is True
