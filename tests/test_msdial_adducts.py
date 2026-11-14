import pytest

from metabolomix.msdial.adducts import (
    adduct_key, ion_counts, mz_from_neutral, neutral_from_mz, parse_adduct, same_adduct,
)
from metabolomix.msdial.peak_verification import ADDUCT_SHIFTS


@pytest.mark.parametrize("name", sorted(ADDUCT_SHIFTS))
def test_parse_agrees_with_the_existing_shift_table(name):
    sign, shift, charge, n_mol = ADDUCT_SHIFTS[name]
    adduct = parse_adduct(name)
    assert adduct is not None
    assert adduct.polarity == sign and adduct.charge == charge and adduct.n_mol == n_mol
    assert adduct.shift == pytest.approx(shift, abs=2e-5)


@pytest.mark.parametrize("name, n_mol, charge, polarity", [
    ("[M+C2H3N+Na-2H]-", 1, 1, "-"), ("[M+CH3COONa-H]-", 1, 1, "-"), ("[2M+FA-H]-", 2, 1, "-"),
    ("[2M+Hac-H]-", 2, 1, "-"), ("[3M-H]-", 3, 1, "-"), ("[M-3H]3-", 1, 3, "-"),
    ("[M+TFA-H]-", 1, 1, "-"), ("[M-C6H10O5-H]-", 1, 1, "-"), ("[M+2H]2+", 1, 2, "+"),
    ("[M+H-H2O]+", 1, 1, "+"), ("[M+K]+", 1, 1, "+"), ("[M+2NH4]2+", 1, 2, "+"),
])
def test_parse_msdial_searched_adducts(name, n_mol, charge, polarity):
    adduct = parse_adduct(name)
    assert (adduct.n_mol, adduct.charge, adduct.polarity) == (n_mol, charge, polarity)


def test_formic_acid_alias_equals_formate():
    # [M+FA-H]- と [M+HCOO]- は同じイオン
    assert parse_adduct("[M+FA-H]-").shift == pytest.approx(parse_adduct("[M+HCOO]-").shift, abs=1e-6)


@pytest.mark.parametrize("name", [None, "", "Unknown", "[M+Xx]+", "M+H", "[M+H]"])
def test_unparseable_is_none(name):
    assert parse_adduct(name) is None


def test_round_trip_between_mz_and_neutral():
    adduct = parse_adduct("[2M+HCOO]-")
    assert neutral_from_mz(mz_from_neutral(759.5778, adduct), adduct) == pytest.approx(759.5778)


# Extra parametrized test cases for real MS-DIAL param files (Task 1 findings)
@pytest.mark.parametrize("name, expected", [
    ("[M+IsoProp+Na+H]+", None),  # Unknown token -> None
    ("[2M+3H2O+2H]+", {"n_mol": 2, "charge": 1, "polarity": "+"}),  # Unusual but valid
])
def test_unusual_adducts(name, expected):
    adduct = parse_adduct(name)
    if expected is None:
        assert adduct is None
    else:
        assert adduct is not None
        for key, value in expected.items():
            assert getattr(adduct, key) == value


# ---------- 同じイオンかを組成で比べる（same_adduct）とイオンの組成（ion_counts） ----------


@pytest.mark.parametrize("a, b", [
    ("[M+FA-H]-", "[M+HCOO]-"), ("[M+FA-H]-", "[M+HCOO]1-"), ("[M+FA-H]-", "[M+CHO2]1-"),
    ("[M+CH3COO]-", "[M+Hac-H]-"), ("[M+CH3COO]-", "[M+C2H3O2]1-"),
    ("[M-H2O-H]-", "[M-H-H2O]-"), ("[M+H]+", "[M+H]1+"), ("[M+H]+", " [M + H]+ "),
    ("[2M+FA-H]-", "[2M+HCOO]1-"), ("[M-2H]2-", "[M-2H]2-"),
])
def test_same_ion_written_differently(a, b):
    assert same_adduct(a, b) and same_adduct(b, a)
    assert adduct_key(a) == adduct_key(b) is not None


@pytest.mark.parametrize("a, b", [
    ("[M+H]+", "[M]+"), ("[M]+", "[M]-"), ("[M-H]-", "[2M-H]-"), ("[M-H]-", "[M+Na-2H]-"),
    ("[M+H]+", "[M+2H]2+"), ("[M+HCOO]-", "[M+CH3COO]-"), ("[M+H]+", None), (None, None), ("", ""),
])
def test_different_ions_are_not_the_same(a, b):
    assert not same_adduct(a, b)


def test_unparseable_adducts_fall_back_to_the_spelling():
    # 略記を解けない表記（IsoProp）は綴りで比べる。価数 1 の `1` と空白・大文字小文字だけを無視する。
    assert adduct_key("[M+IsoProp+H]+") is None
    assert same_adduct("[M+IsoProp+H]1+", "[M+IsoProp+H]+")
    assert same_adduct("[m+isoprop+h]+", "[M+IsoProp+H]+")
    assert not same_adduct("[M+IsoProp+H]+", "[M+IsoProp+Na]+")
    assert not same_adduct("[M+IsoProp+H]2+", "[M+IsoProp+H]+")


def test_ion_counts_multiply_by_n_mol_and_add_the_adduct_atoms():
    neutral = {"C": 20, "H": 40, "O": 2}
    assert ion_counts(neutral, "[2M+H]+") == {"C": 40, "H": 81, "O": 4}
    assert ion_counts(neutral, "[M+Cl]-") == {"C": 20, "H": 40, "O": 2, "Cl": 1}
    assert ion_counts(neutral, "[M+HCOO]1-") == {"C": 21, "H": 41, "O": 4}
    assert ion_counts({"C": 1, "H": 4}, "[M-H2O+H]+") is None      # O が負になる
    assert ion_counts(neutral, "[M+Xx]+") is None
    assert ion_counts(neutral, None) is None
