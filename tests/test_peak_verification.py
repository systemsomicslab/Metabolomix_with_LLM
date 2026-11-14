"""peak_verification 純ロジックの検証。"""

import textwrap
import unittest

from metabolomix.msdial import peak_verification as pv
from tests.config_fixtures import use_config


class FormulaMassTests(unittest.TestCase):
    def test_parse_formula_multi_element_and_multi_digit(self):
        self.assertEqual(
            pv.parse_formula("C42H82NO8P"),
            {"C": 42, "H": 82, "N": 1, "O": 8, "P": 1},
        )

    def test_parse_formula_single_atom_defaults_to_one(self):
        self.assertEqual(pv.parse_formula("CH4"), {"C": 1, "H": 4})

    def test_parse_formula_empty_raises(self):
        with self.assertRaises(ValueError):
            pv.parse_formula("")

    def test_monoisotopic_mass_pc_34_1(self):
        mass = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        self.assertAlmostEqual(mass, 759.5778, places=2)

    def test_monoisotopic_mass_unknown_element_raises(self):
        with self.assertRaises(KeyError):
            pv.monoisotopic_mass({"Xx": 1})


class AdductMassErrorTests(unittest.TestCase):
    def test_adduct_mz_protonated(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        self.assertAlmostEqual(pv.adduct_mz(neutral, "[M+H]+"), 760.5851, places=2)

    def test_adduct_mz_ammonium(self):
        self.assertAlmostEqual(pv.adduct_mz(100.0, "[M+NH4]+"), 118.03383, places=4)

    def test_adduct_mz_deprotonated(self):
        self.assertAlmostEqual(pv.adduct_mz(100.0, "[M-H]-"), 98.99272, places=4)

    def test_adduct_mz_unknown_returns_none(self):
        self.assertIsNone(pv.adduct_mz(100.0, "[M+ZZ]+"))

    def test_mass_error_pass(self):
        result = pv.mass_error_ppm(760.5851, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "PASS")
        self.assertAlmostEqual(result["ppm"], 0.0, delta=5.0)

    def test_mass_error_borderline(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        theo = pv.adduct_mz(neutral, "[M+H]+")
        observed = theo * (1 + 8e-6)  # +8 ppm
        result = pv.mass_error_ppm(observed, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "BORDERLINE")

    def test_mass_error_fail(self):
        neutral = pv.monoisotopic_mass(pv.parse_formula("C42H82NO8P"))
        theo = pv.adduct_mz(neutral, "[M+H]+")
        observed = theo * (1 + 20e-6)  # +20 ppm
        result = pv.mass_error_ppm(observed, "C42H82NO8P", "[M+H]+")
        self.assertEqual(result["band"], "FAIL")

    def test_mass_error_unknown_formula(self):
        result = pv.mass_error_ppm(760.0, "Unknown", "[M+H]+")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["ppm"])

    def test_mass_error_unknown_adduct(self):
        result = pv.mass_error_ppm(760.0, "C42H82NO8P", "Unknown")
        self.assertEqual(result["band"], "UNKNOWN")


class AdductConsistencyTests(unittest.TestCase):
    def test_polarity_match_positive(self):
        result = pv.adduct_consistency("[M+H]+", "Positive", "PC")
        self.assertTrue(result["polarity_ok"])
        self.assertEqual(result["band"], "PASS")
        self.assertTrue(result["class_typical"])

    def test_polarity_mismatch(self):
        result = pv.adduct_consistency("[M+H]+", "Negative", "PC")
        self.assertFalse(result["polarity_ok"])
        self.assertEqual(result["band"], "FAIL")

    def test_class_atypical_adduct_is_advisory_not_fail(self):
        result = pv.adduct_consistency("[M+Na]+", "Positive", "PC")
        self.assertTrue(result["polarity_ok"])
        self.assertEqual(result["band"], "PASS")
        self.assertFalse(result["class_typical"])

    def test_unknown_class_gives_none_typical(self):
        result = pv.adduct_consistency("[M+H]+", "Positive", "ZZZ")
        self.assertIsNone(result["class_typical"])

    def test_unknown_adduct_is_unknown_band(self):
        result = pv.adduct_consistency("Unknown", "Positive", "PC")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["polarity_ok"])

    def test_unrecognized_ion_mode_is_unknown_not_fail(self):
        result = pv.adduct_consistency("[M+H]+", "Both", "PC")
        self.assertEqual(result["band"], "UNKNOWN")
        self.assertIsNone(result["polarity_ok"])
        self.assertTrue(result["class_typical"])  # class advisory still computed


class ClassTokenTests(unittest.TestCase):
    def test_extract_class_token_prefers_ontology(self):
        self.assertEqual(pv.extract_class_token("PC 34:1", "PC"), "pc")

    def test_extract_class_token_falls_back_to_name(self):
        self.assertEqual(pv.extract_class_token("TG 52:2", ""), "tg")

    def test_extract_class_token_empty(self):
        self.assertIsNone(pv.extract_class_token("", ""))


class TestExpandedTables(unittest.TestCase):
    def test_dimer_adduct_mz(self):
        # [2M-H]- : (2*neutral - proton) / 1
        neutral = 100.0
        mz = pv.adduct_mz(neutral, "[2M-H]-")
        self.assertAlmostEqual(mz, 2 * 100.0 - pv.PROTON_MASS, places=4)

    def test_doubly_charged_divides(self):
        neutral = 800.0
        mz = pv.adduct_mz(neutral, "[M-2H]2-")
        self.assertAlmostEqual(mz, (800.0 - 2 * pv.PROTON_MASS) / 2, places=4)

    def test_new_elements_present(self):
        for el in ("D", "F", "Br", "C13"):
            self.assertIn(el, pv.ELEMENT_MASSES)

    def test_singly_charged_unchanged(self):
        # 既存の1価挙動が数値的に不変であることを保証。
        self.assertAlmostEqual(pv.adduct_mz(500.0, "[M+H]+"), 500.0 + pv.PROTON_MASS, places=6)
        self.assertAlmostEqual(pv.adduct_mz(500.0, "[M+HCOO]-"), 500.0 + 44.99820286, places=6)


class CompositionAdductTests(unittest.TestCase):
    """表（ADDUCT_SHIFTS）に無い MS-DIAL の綴りも組成で読む（最終検証 F28）。

    理論 m/z は元素質量から手で組み立てた値と比べる（実装の parse_adduct を使わない独立の計算）。
    """

    H, O, K, E = (pv.ELEMENT_MASSES["H"], pv.ELEMENT_MASSES["O"], pv.ELEMENT_MASSES["K"],
                  pv.ELECTRON_MASS)

    def _neutral(self, formula):
        return pv.monoisotopic_mass(pv.parse_formula(formula))

    def test_msdial_spellings_get_a_theoretical_mz(self):
        cer = self._neutral("C34H67NO3")
        water = 2 * self.H + self.O
        cases = {
            "[M+H-H2O]+": cer + self.H - water - self.E,
            "[M+K]+": cer + self.K - self.E,
            "[M-H2O-H]-": cer - water - self.H + self.E,
            "[M+H]1+": cer + self.H - self.E,
            "[2M+H]+": 2 * cer + self.H - self.E,
            "[M+2H]2+": (cer + 2 * self.H - 2 * self.E) / 2,
        }
        for adduct, theoretical in cases.items():
            observed = theoretical * (1 + 8e-6)                      # +8 ppm
            result = pv.mass_error_ppm(observed, "C34H67NO3", adduct)
            self.assertEqual(result["band"], "BORDERLINE", adduct)
            self.assertAlmostEqual(result["ppm"], 8.0, places=2)
            self.assertAlmostEqual(result["theoretical_mz"], round(theoretical, 4), places=4)

    def test_msdial_water_loss_spelling_matches_the_table_spelling(self):
        observed = 520.5088
        a = pv.mass_error_ppm(observed, "C34H67NO3", "[M+H-H2O]+")
        b = pv.mass_error_ppm(observed, "C34H67NO3", "[M-H2O+H]+")
        self.assertEqual(a, b)
        self.assertEqual(a["band"], "PASS")

    def test_polarity_is_read_from_the_composition(self):
        for adduct, mode, band in (("[M+H-H2O]+", "Negative", "FAIL"), ("[M+K]+", "Negative", "FAIL"),
                                   ("[M-H2O-H]-", "Positive", "FAIL"), ("[M+CHO2]-", "Negative", "PASS"),
                                   ("[M+H]1+", "Positive", "PASS"), ("[M+H-H2O]+", "Positive", "PASS")):
            self.assertEqual(pv.adduct_consistency(adduct, mode)["band"], band, (adduct, mode))
        self.assertEqual(pv.adduct_consistency("[M+ZZ]+", "Positive")["band"], "UNKNOWN")

    def test_class_typical_adducts_are_compared_by_composition(self):
        for adduct in ("[M+HCOO]-", "[M+FA-H]-", "[M+CHO2]-", "[M+HCOO]1-", "[M+Hac-H]-"):
            self.assertIs(pv.adduct_consistency(adduct, "Negative", "PC")["class_typical"], True, adduct)
        self.assertIs(pv.adduct_consistency("[M-H]-", "Negative", "PC")["class_typical"], False)
