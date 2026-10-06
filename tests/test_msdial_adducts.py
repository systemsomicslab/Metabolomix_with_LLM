import pytest

from metabolomix.msdial.adducts import mz_from_neutral, neutral_from_mz, parse_adduct
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
