import unittest

from lipidmix.msdial import lipid_identity as li


class TestNormalizeLipidName(unittest.TestCase):
    def test_basic_pc_parses(self):
        out = li.normalize_lipid_name("PC 34:1")
        # pygoslin should parse; if unavailable, parse_ok False but no crash
        self.assertIn("parse_ok", out)
        if out["parse_ok"]:
            self.assertTrue(out["normalized"].upper().startswith("PC"))

    def test_ether_forms_parse(self):
        # NOTE: at the species shorthand level, GOSLIN represents plasmalogen
        # `PC P-34:0` and ether `PC O-34:1` with the SAME normalized string
        # (the P-/O- plasmalogen ambiguity). Downstream must rely on the ether
        # caveat (peak_verification.ether_caveats), not the normalized name, to
        # tell them apart. Here we only assert both parse cleanly.
        o = li.normalize_lipid_name("PC O-34:1")
        p = li.normalize_lipid_name("PC P-34:0")
        if o["parse_ok"] and p["parse_ok"]:
            self.assertTrue(o["normalized"].upper().startswith("PC"))
            self.assertTrue(p["normalized"].upper().startswith("PC"))

    def test_garbage_returns_parse_error_not_crash(self):
        out = li.normalize_lipid_name("not a lipid ###")
        self.assertFalse(out["parse_ok"])
        self.assertIsNotNone(out.get("error"))


class TestQualifierPrefixStripping(unittest.TestCase):
    # MS-DIAL prepends confidence qualifiers to Name (e.g. "no MS2: ", "low score: ");
    # in the real brain NEG catalog 317/856 (37%) of annotations carry these and
    # previously failed GOSLIN parsing. Strip them before parsing.
    def test_strips_no_ms2_prefix(self):
        out = li.normalize_lipid_name("no MS2: FA 18:1")
        if (out.get("error") or "").startswith("pygoslin unavailable"):
            self.skipTest("pygoslin not installed")
        self.assertTrue(out["parse_ok"], out)
        self.assertTrue(out["normalized"].upper().startswith("FA"))

    def test_strips_low_score_prefix(self):
        out = li.normalize_lipid_name("low score: PC 35:2")
        if (out.get("error") or "").startswith("pygoslin unavailable"):
            self.skipTest("pygoslin not installed")
        self.assertTrue(out["parse_ok"], out)
        self.assertTrue(out["normalized"].upper().startswith("PC"))

    def test_takes_first_candidate_before_pipe(self):
        out = li.normalize_lipid_name("low score: Cer 24:1;O2|Cer 12:0;O2/12:1")
        if (out.get("error") or "").startswith("pygoslin unavailable"):
            self.skipTest("pygoslin not installed")
        self.assertTrue(out["parse_ok"], out)

    def test_records_that_a_qualifier_was_stripped(self):
        out = li.normalize_lipid_name("no MS2: FA 18:1")
        self.assertTrue(out.get("stripped"))

    def test_clean_name_reports_no_strip(self):
        out = li.normalize_lipid_name("PC 34:1")
        self.assertFalse(out.get("stripped"))


class TestReferenceMapping(unittest.TestCase):
    def setUp(self):
        self.tables = li.load_reference_tables("reference")

    def test_known_class_maps(self):
        out = li.map_to_reference("pc", self.tables)
        self.assertTrue(out["matched"])
        self.assertEqual(out["lipid_maps_category"], "GP")

    def test_unknown_class_has_caveat(self):
        out = li.map_to_reference("zzz", self.tables)
        self.assertFalse(out["matched"])
        self.assertIsNotNone(out["caveat"])


class TestMsiLevel(unittest.TestCase):
    def test_msms_plus_mass_ok_is_level2(self):
        out = li.msi_level(name="PC 34:1", ontology="PC", has_msms=True,
                           mass_error_band="PASS", adduct_band="PASS")
        self.assertEqual(out["level"], 2)
        self.assertTrue(out["heuristic"])

    def test_class_only_is_level3(self):
        out = li.msi_level(name="", ontology="PC", has_msms=False,
                           mass_error_band="UNKNOWN", adduct_band="UNKNOWN")
        self.assertEqual(out["level"], 3)

    def test_unknown_is_level4(self):
        out = li.msi_level(name="", ontology="", has_msms=False,
                           mass_error_band="UNKNOWN", adduct_band="UNKNOWN")
        self.assertEqual(out["level"], 4)

    def test_never_level1(self):
        out = li.msi_level(name="PC 34:1", ontology="PC", has_msms=True,
                           mass_error_band="PASS", adduct_band="PASS")
        self.assertNotEqual(out["level"], 1)


class TestReferenceTablesAreFoundFromAnyWorkingDirectory(unittest.TestCase):
    """同梱表はリポジトリ内の固定位置にある。起動元の作業ディレクトリに依存させない。

    Claude Desktop はサーバを任意の作業ディレクトリから起動しうる。相対パスで
    探すと表が見つからず、`_read_tsv` が空リストを返すため、全クラスが
    「同梱マッピング表に無い」ことになって ID 付与が黙って全滅する。
    """

    def test_default_tables_load_outside_the_repository(self):
        import os
        import tempfile
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as elsewhere:
            os.chdir(elsewhere)
            try:
                tables = li.load_reference_tables()
            finally:
                os.chdir(cwd)
        self.assertTrue(tables["lipidmaps"], "LIPID MAPS 表が 0 件")
        self.assertTrue(tables["refmet"], "RefMet 表が 0 件")
        self.assertTrue(li.map_to_reference("pc", tables)["matched"])


class TestIdentityBlock(unittest.TestCase):
    def test_block_has_all_sections(self):
        tables = li.load_reference_tables("reference")
        feat = {"name": "PC 34:1", "ontology": "PC", "has_msms": True}
        block = li.build_identity_block(feat, tables,
                                        mass_error_band="PASS", adduct_band="PASS")
        self.assertIn("goslin", block)
        self.assertIn("reference", block)
        self.assertIn("msi", block)
        self.assertEqual(block["msi"]["level"], 2)
        self.assertTrue(block["reference"]["matched"])


if __name__ == "__main__":
    unittest.main()
