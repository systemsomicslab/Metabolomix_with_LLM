import pytest

from metabolomix.curation import isotope
from metabolomix.msdial.adducts import adduct_composition


def test_adduct_composition_counts_the_added_and_removed_atoms():
    assert adduct_composition("[M+H]+") == (1, 1, {"H": 1})
    assert adduct_composition("[M+NH4]+") == (1, 1, {"N": 1, "H": 4})
    assert adduct_composition("[M+FA-H]-") == (1, 1, {"C": 1, "H": 1, "O": 2})
    assert adduct_composition("[M-H2O+H]+") == (1, 1, {"H": -1, "O": -1})
    assert adduct_composition("[2M-H]-") == (2, 1, {"H": -1})
    assert adduct_composition("[M-2H]2-") == (1, 2, {"H": -2})
    assert adduct_composition("[M+Xx]+") is None
    assert adduct_composition("") is None


def test_charged_formula_gives_the_same_envelope_as_the_uncharged_one():
    charged = isotope.theoretical_envelope("C5H14NO+", "[M]+")
    assert charged is not None
    assert charged == isotope.theoretical_envelope("C5H14NO", "[M]+")


def test_pc_34_1_protonated_has_the_textbook_isotope_ratios():
    # C42H83NO8P+: M+1 ≈ 47 %（炭素 42 個 × 1.07 % が主）、M+2 ≈ 12 %。
    env = isotope.theoretical_envelope("C42H82NO8P", "[M+H]+")
    assert env["basis"] == "formula+adduct"
    m0, m1, m2 = env["relative"]
    assert m0 == 100.0
    assert m1 == pytest.approx(47.1, abs=0.5)
    assert 11.0 < m2 < 14.0


def test_chlorine_adduct_raises_m_plus_2():
    with_cl = isotope.theoretical_envelope("C42H82NO8P", "[M+Cl]-")["relative"]
    without = isotope.theoretical_envelope("C42H82NO8P", "[M-H]-")["relative"]
    assert with_cl[2] > without[2] + 25          # 37Cl ≈ 32 % of 35Cl


def test_unknown_adduct_falls_back_to_the_neutral_formula():
    env = isotope.theoretical_envelope("C42H82NO8P", "[M+Xx]+")
    assert env["basis"] == "formula"
    assert env["relative"][1] == pytest.approx(47.0, abs=0.5)


def test_dimer_doubles_the_atom_counts():
    mono = isotope.theoretical_envelope("C16H32O2", "[M-H]-")["relative"]
    dimer = isotope.theoretical_envelope("C16H32O2", "[2M-H]-")["relative"]
    assert dimer[1] == pytest.approx(2 * mono[1], rel=0.03)


@pytest.mark.parametrize("formula", [None, "", "Unknown", "C10Xq3", "C5H11NO2Se"])
def test_unreadable_formula_gives_none(formula):
    # Se は同位体の表に無い（最も多い 80Se の下に軽い同位体がある）ので計算しない
    assert isotope.theoretical_envelope(formula, "[M+H]+") is None


def test_monoisotopic_iodine_does_not_block_the_envelope():
    # チロキシン C15H11I4NO4 [M-H]-（イオン C15H10I4NO4）。ヨウ素は 127I だけで同位体のピークを作らない。
    # M+1 = 13C 15·q_C + 2H 10·q_H + 15N q_N + 17O 4·q_O（q = a / (1 − a)）
    hand = (15 * 0.0107 / 0.9893 + 10 * 0.000115 / 0.999885 + 0.00364 / 0.99636 + 4 * 0.00038 / 0.99757)
    env = isotope.theoretical_envelope("C15H11I4NO4", "[M-H]-")
    assert env["relative"][1] == pytest.approx(100 * hand, abs=0.01)


def test_every_heavy_isotope_has_its_mass_shift():
    for element, dist in isotope.ABUNDANCES.items():
        heavy = {k for k in dist if k > 0}
        assert heavy == set(isotope.ISOTOPE_SHIFTS.get(element, {})), element
        for k, shift in isotope.ISOTOPE_SHIFTS.get(element, {}).items():
            assert abs(shift - k) < 0.01, (element, k)       # 質量差は整数質量のずれに近い


def test_measured_envelope_keeps_mz_relative_and_absolute_abundance():
    raw = [[100.0, 760.5851, 0.0, "", 10000.4], [46.123456, 761.5898, 1.0047, "", 4612.0],
           [12.5, 762.5945, 2.0093, "", 1250.0]]
    assert isotope.measured_envelope(raw) == [[760.5851, 100.0, 10000], [761.5898, 46.12, 4612],
                                              [762.5945, 12.5, 1250]]
    # 絶対強度（Key 4）が無い・読めなければ null（相対強度は出す）
    assert isotope.measured_envelope([[100.0, 760.5851, 0.0, ""]]) == [[760.5851, 100.0, None]]
    assert isotope.measured_envelope(None) == []
    assert isotope.measured_envelope([["bad"]]) == []
