# tests/test_feature_search_lipid_species.py
"""脂質の種の目標（`lipid_species`。InChIKey が一意に決まらない主張）の照合。

InChIKey の目標の正規化・rule_hash・照合の出力は、この口を足す前の値で固定する
（既存の挙動を変えない。ユーザー判断 2026-10-22）。
"""
import pytest

from metabolomix.analysis.feature_search import (
    classify_identity,
    match_target,
    normalize_targets,
    rule_hash,
    theoretical_mz,
)
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.msdial.lipid_identity import lipid_species_of
from tests.test_feature_search import GABA_KEY, MH, _candidate, _dataset, _row, _reasons

PC_FORMULA = "C42H76NO8P"          # PC 34:4
PC_MH = theoretical_mz(PC_FORMULA, "[M+H]+")
PC_KEY = "ABCDEFGHIJKLMN-UHFFFAOYSA-N"


def _species(**extra):
    target = {"target_id": "C1", "lipid_species": "PC 18:3_16:1", "formula": PC_FORMULA,
              "rt_min": 8.4, **extra}
    [normalized] = normalize_targets([target], dataset_polarity="positive")
    return normalized


def _pc(name, **overrides):
    fields = {"chemical_name": name, "inchikey": None, "inchikey_source": "none",
              "exp_mass_to_charge": PC_MH, "theoretical_mass_to_charge": PC_MH, **overrides}
    return _candidate(**fields)


# ---------- 種のレベルの読み ----------

@pytest.mark.parametrize("name, species", [
    ("PC 34:4", "PC 34:4"),
    ("PC 18:3_16:1", "PC 34:4"),
    ("PC 16:1/18:3", "PC 34:4"),
    ("PC 34:4|PC 18:3_16:1", "PC 34:4"),
    ("w/o MS2: PC 18:3_16:1", "PC 34:4"),
    ("PC O-34:1", "PC O-34:1"),
    ("PC P-34:1", "PC O-34:2"),        # GOSLIN の規約: P- はビニルエーテルの二重結合を数えて O- に
    ("PE 16:0_18:1;O", "PE 34:1;O"),
    ("Cer 18:1;O2/16:0", "Cer 34:1;O2"),
])
def test_lipid_species_of_reduces_to_species_level(name, species):
    assert lipid_species_of(name) == species


@pytest.mark.parametrize("name", [None, "", "PC", "Glucose", "Unknown", 34])
def test_lipid_species_of_rejects_coarser_or_unreadable_names(name):
    assert lipid_species_of(name) is None


# ---------- 目標の正規化 ----------

def test_species_target_is_reduced_and_keeps_the_input():
    target = _species()
    assert target["lipid_species"] == "PC 34:4"
    assert target["lipid_species_input"] == "PC 18:3_16:1"
    assert target["inchikey_block"] is None
    assert target["theoretical_mz"]["[M+H]+"] == pytest.approx(PC_MH)


def test_species_rule_hash_ignores_the_molecular_spelling_but_not_the_species():
    assert rule_hash([_species()]) == rule_hash([_species(lipid_species="PC 34:4")])
    assert rule_hash([_species()]) != rule_hash([_species(lipid_species="PC 34:3")])
    assert rule_hash([_species()]) != rule_hash([_species(lipid_species="PC O-34:4")])


@pytest.mark.parametrize("extra, field", [
    ({"inchikey": PC_KEY}, "lipid_species"),          # 排他
    ({"lipid_species": "PC"}, "lipid_species"),        # クラスだけ
    ({"lipid_species": "Glucose"}, "lipid_species"),
    ({"lipid_species": 34}, "lipid_species"),
    ({"lipid_species": ""}, "lipid_species"),
])
def test_invalid_species_targets_name_the_field(extra, field):
    raw = {"target_id": "C1", "lipid_species": "PC 34:4", "formula": PC_FORMULA, **extra}
    with pytest.raises(DomainError) as exc:
        normalize_targets([raw], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details == {"target_id": "C1", "field": field}


def test_target_without_any_identifier_still_names_inchikey():
    with pytest.raises(DomainError) as exc:
        normalize_targets([{"target_id": "c1", "formula": "C4H9NO2"}], dataset_polarity="positive")
    assert exc.value.details["field"] == "inchikey"


# ---------- InChIKey の目標は変わらない（変更前の値で固定） ----------

def test_inchikey_targets_are_byte_identical_to_before():
    [t] = normalize_targets([{"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
                              "rt_min": 1.2}], dataset_polarity="positive")
    assert sorted(t) == ["adduct_tiers", "applicable", "expected_mz", "formula", "inchikey_block",
                         "name", "polarity", "ppm", "rt_min", "rt_tol_min", "target_id",
                         "theoretical_mz"]
    assert rule_hash([t]) == "1890f39fd2d7408cba18962e9c6028180bad10b1f11fd6cdc2d216e4d4324d19"
    t2 = normalize_targets([{"target_id": "c2", "inchikey": GABA_KEY, "expected_mz": 104.0706,
                             "adduct": "[M+H]+"}], dataset_polarity=None)
    assert rule_hash(t2) == "33131cfe7e1b022b54412d0d6d10b0fe8a7bfdfcf947f65a4ef3cfd30082b7f1"


def test_inchikey_matching_output_is_unchanged():
    ds = _dataset({
        "1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]},
        "2": {"mz": MH, "rt": 1.25, "candidates": [_candidate(inchikey="QWERTYUIOPASDF-UHFFFAOYSA-N")]},
        "3": {"mz": MH, "rt": 1.3},
    })
    [t] = normalize_targets([{"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
                              "rt_min": 1.2}], dataset_polarity="positive")
    out = match_target(ds, t, max_results=20)
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"
    assert canonical_hash(out) == PINNED_MATCH_HASH


PINNED_MATCH_HASH = "5ff0be55dc08b1bfc345b3155bd94de05e1b7395739da1715a3fb1c95265a4ba"


# ---------- 種の目標の照合 ----------

def test_molecular_species_candidate_matches_the_species_target():
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [_pc("PC 16:1_18:3")]}})
    out = match_target(ds, _species(), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "sme_match" and row["qualified"] is True
    assert row["annotated_name"] == "PC 16:1_18:3"
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"
    assert out["selected_evidence"] == "sme_match"


def test_species_target_ignores_the_candidates_inchikey():
    # 種の目標は名前で照らす。候補の InChIKey（どの分子種の構造か）は見ない
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4,
                         "candidates": [_pc("PC 34:4", inchikey=PC_KEY, inchikey_source="msp")]}})
    out = match_target(ds, _species(), max_results=20)
    assert _row(out, "1")["evidence"] == "sme_match"


def test_other_class_with_the_same_mass_is_contradicted():
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [_pc("PE 37:4")]}})
    out = match_target(ds, _species(), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "contradicted" and row["qualified"] is False
    assert _reasons(row) == ["identifier_mismatch"]
    assert out["status"] == "needs_input" and out["reason"] == "no_qualified_candidate"


@pytest.mark.parametrize("name", ["PC O-34:4", "PC 34:3", "PC 34:4;O", "PC 36:4"])
def test_ether_modified_or_other_species_is_not_the_target(name):
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [_pc(name)]}})
    assert _row(match_target(ds, _species(), max_results=20), "1")["evidence"] == "contradicted"


def test_ether_target_matches_ether_and_plasmalogen_spellings():
    target = _species(lipid_species="PC O-34:2")
    assert target["lipid_species"] == "PC O-34:2"
    for name in ("PC O-16:0_18:2", "PC P-34:1"):
        ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [_pc(name)]}})
        assert _row(match_target(ds, target, max_results=20), "1")["evidence"] == "sme_match", name
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [_pc("PC 34:2")]}})
    assert _row(match_target(ds, target, max_results=20), "1")["evidence"] == "contradicted"


def test_two_features_of_the_species_are_needs_input():
    # 種のレベルでは区別できない 2 つの分子種（RT が近い）を規則は選ばない
    ds = _dataset({
        "1": {"mz": PC_MH, "rt": 8.35, "candidates": [_pc("PC 16:1_18:3")]},
        "2": {"mz": PC_MH, "rt": 8.50, "candidates": [_pc("PC 14:0_20:4")]},
    })
    out = match_target(ds, _species(), max_results=20)
    assert out["status"] == "needs_input" and out["reason"] == "multiple_qualified_candidates"
    assert out["selected_evidence"] == "sme_match" and out["selected_feature_id"] is None
    assert out["counts"]["by_evidence"]["sme_match"] == 2


def test_rank2_species_candidate_is_found():
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4,
                         "candidates": [_pc("PE 37:4"), _pc("PC 34:4", rank=2, sme_id="2")]}})
    row = _row(match_target(ds, _species(), max_results=20), "1")
    assert row["evidence"] == "sme_match" and row["sme_rank"] == 2


def test_unreadable_candidate_name_is_unresolved_and_falls_to_sml():
    unreadable = _pc("Unknown")
    sml = {"sml_id": "1", "ambiguous": False, "name": "PC 34:4", "inchikey": None}
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4, "candidates": [unreadable], "annotation": sml}})
    row = _row(match_target(ds, _species(), max_results=20), "1")
    assert row["evidence"] == "sml_match"
    assert "sme_identity_unresolved" in _reasons(row)
    # SML が別の種なら矛盾ではなく unannotated（SML は MS1 照合だけ）
    other = dict(sml, name="PE 37:4")
    assert classify_identity([], other, lipid_species="PC 34:4", adduct="[M+H]+") \
        == ("unannotated", [], None)


def test_ambiguous_sml_cannot_be_read_for_species():
    ambiguous = {"ambiguous": True, "sml_ids": ["1", "2"], "name": None, "inchikeys": []}
    assert classify_identity([], ambiguous, lipid_species="PC 34:4", adduct="[M+H]+") \
        == ("unannotated", [], None)


def test_species_identity_still_checks_the_adduct():
    ds = _dataset({"1": {"mz": PC_MH, "rt": 8.4,
                         "candidates": [_pc("PC 34:4", adduct="[M+Na]1+")]}})
    row = _row(match_target(ds, _species(), max_results=20), "1")
    assert row["evidence"] == "contradicted" and _reasons(row) == ["adduct_mismatch"]


def test_classify_identity_needs_exactly_one_identifier():
    with pytest.raises(ValueError):
        classify_identity([], None, adduct="[M+H]+")
    with pytest.raises(ValueError):
        classify_identity([], None, inchikey_block="BTCSSZJGUNDROE", lipid_species="PC 34:4",
                          adduct="[M+H]+")
