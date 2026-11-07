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
