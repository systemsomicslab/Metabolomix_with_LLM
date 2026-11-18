import unittest

from metabolomix.msdial import lipid_identity as li


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
        # (the P-/O- plasmalogen ambiguity). Downstream must not rely on the
        # normalized name alone to tell them apart. Here we only assert both parse cleanly.
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


# ---------- GOSLIN のパーサを使い回す（B-1 実データ照合） ----------

def test_goslin_parser_is_built_once_and_output_is_unchanged(monkeypatch):
    import pygoslin.parser.Parser as goslin

    from metabolomix.msdial import lipid_identity as li

    names = ["PC 34:1", "PE 16:0_18:1", "no MS2: FA 5:0", "Cer 24:1;O2|Cer 12:0;O2/12:1"]
    real = goslin.LipidParser
    # 名前ごとに新しいパーサを作った場合の答え（基準）。共有のパーサを名前のたびに捨てて、
    # 前の名前の状態が次の名前へ漏れない基準を作る（共有のパーサで取ると漏れを検出できない）。
    expected_norm, expected_species = [], []
    for name in names:
        monkeypatch.setattr(li, "_SHARED_PARSER", None)
        li.lipid_species_of.cache_clear()
        expected_norm.append(li.normalize_lipid_name(name))
        expected_species.append(li.lipid_species_of(name))

    built = []

    class CountingParser(real):
        def __init__(self, *args, **kwargs):
            built.append(1)
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(goslin, "LipidParser", CountingParser)
    monkeypatch.setattr(li, "_SHARED_PARSER", None)
    li.lipid_species_of.cache_clear()
    for _ in range(3):
        assert [li.normalize_lipid_name(n) for n in names] == expected_norm
        li.lipid_species_of.cache_clear()
        assert [li.lipid_species_of(n) for n in names] == expected_species
    assert len(built) == 1
    li.lipid_species_of.cache_clear()



# ---------- 脂質名の文法で読めるか（脂質アッセイのクラスの集計から外す判定） ----------

def test_lipid_name_parses_accepts_msdial_lipid_names_and_rejects_natural_products():
    lipids = ["PC 34:1", "TG 16:0_18:1_18:2", "Cer 18:1;O2/16:0", "SM d34:1", "CAR 18:1", "ST 27:1;O",
              "w/o MS2: PC 34:1", "low score: TG 52:2", "PC 34:1|PC 16:0_18:1", "PE 36:2;O (oxPE)",
              "LPC 16:0; ExampleDB ID-0000", "PC 34:1 [M+H]+", "Cer[NS] d30:1", "Phosphatidylcholine 16:0-18:1",
              "PC 15:0-18:1(d7)", "Oleic acid", "Cholesterol"]
    non_lipids = ["Ginsenoside Rb1", "Soyasaponin I", "Enniatin B", "Valinomycin", "Beauvericin",
                  "w/o MS2: Saikosaponin A", "Triton X-100", "Xanthine"]
    assert [n for n in lipids if li.lipid_name_parses(n) is not True] == []
    assert [n for n in non_lipids if li.lipid_name_parses(n) is not False] == []


def test_lipid_annotation_accepts_msdial_lipid_classes_even_when_goslin_cannot_read_the_name():
    # MS-DIAL の略記はクラストークン＝Ontology。GOSLIN が読めない形も脂質の注釈
    for name, ontology in [("Cer 12:0;O3/8:0(2OH)", "Cer_AP"), ("AHexCer (O-14:1)12:1;O2/12:1", "AHexCer"),
                           ("SM 24:1;O2(FA 12:0)", "ASM"), ("NAGly 17:0", "NAGly"), ("SE 28:4/16:1", "DEGSE")]:
        assert li.lipid_name_parses(name) is False, name
        assert li.is_lipid_annotation(name, ontology) is True, name
    # MS-DIAL の脂質クラスでなければ名前で決める（ClassyFire 系のクラス・Unknown）
    assert li.is_lipid_annotation("Oleic acid", "Long-chain fatty acids") is True
    assert li.is_lipid_annotation("Saikosaponin A", "Triterpene saponins") is False
    # 組成式が無ければ ClassyFire のクラスからは救えない（名前の IUPAC 名は GOSLIN で読めない）
    assert li.is_lipid_annotation("1,2-Dioleoyl-sn-glycerol", "1,2-diacylglycerols") is False
    assert li.is_lipid_annotation("no MS2: Unknown (carbon number 10); PlaSMA ID-657", "Unknown") is False
    assert li.is_lipid_annotation("Enniatin B", None) is False


def test_lab_msp_lipids_with_unreadable_names_are_rescued_by_classyfire_class_and_formula():
    # 規則 3: ClassyFire の脂質クラスと組成式から作った和の組成の略記を GOSLIN が読めば脂質。略記のクラスへ移す
    cases = [("no MS2: 1,2-Dioleoyl-sn-glycero-3-phosphocholine", "Phosphatidylcholines", "C44H84NO8P", "PC"),
             ("1,2-Dioleoyl-sn-glycerol", "1,2-diacylglycerols", "C39H72O5", "DG"),
             ("low score: MassSpecGymID0000001", "Triacylglycerols", "C57H104O6", "TG"),
             ("MassSpecGymID0000002", "1-alkyl,2-acylglycero-3-phosphocholines", "C42H84NO7P", "EtherPC"),
             ("L-Propionylcarnitine", "Acyl carnitines", "C10H19NO4", "CAR")]
    for name, ontology, formula, lbm in cases:
        call = li.lipid_annotation_class(name, ontology, formula)
        assert call == (True, lbm, li.ROUTE_CLASS_FORMULA, True), name
        assert li.is_lipid_annotation(name, ontology, formula) is True
    # 組成式がクラスの形に合わない（クラスの付け違い）・MS-DIAL で複数のクラスに分かれる・表に無いクラスは救わない
    for name, ontology, formula in [("7.alpha.,24(S)-Dihydroxy-4-cholesten-3-one", "Phosphatidylcholines", "C27H44O3"),
                                    ("N-Octanoylsphingosine", "Ceramides", "C26H51NO3"),
                                    ("Saikosaponin A", "Triterpene saponins", "C42H68O13"),
                                    ("1,2-Dioleoyl-sn-glycerol", "1,2-diacylglycerols", None)]:
        assert li.lipid_annotation_class(name, ontology, formula) == (False, None, li.ROUTE_UNMATCHED, False), name


def test_lipid_annotation_class_moves_classyfire_lipids_to_the_msdial_class():
    # MS-DIAL の脂質クラスはそのまま。読める名前は頭部（鎖のある名前だけ）か、クラスと組成式から MS-DIAL のクラスへ移す
    assert li.lipid_annotation_class("PC 34:1", "PC") == (True, "PC", li.ROUTE_LBM_CLASS, False)
    assert li.lipid_annotation_class("no MS2: PC 18:0/18:1(9Z)", "Phosphatidylcholines", "C44H86NO8P") == (
        True, "PC", li.ROUTE_NAME, False)
    assert li.lipid_annotation_class("PC(16:0/18:1(9Z))", "Phosphatidylcholines", "C42H82NO8P") == (
        True, "PC", li.ROUTE_CLASS_FORMULA, False)
    assert li.lipid_annotation_class("Pentadecanoic acid", "Long-chain fatty acids", "C15H30O2") == (
        True, "FA", li.ROUTE_CLASS_FORMULA, False)
    # 読めるがクラスが決まらない注釈は元のクラスのまま
    assert li.lipid_annotation_class("16-Hydroxyhexadecanoic acid", "Long-chain fatty acids", "C16H32O3") == (
        True, "Long-chain fatty acids", li.ROUTE_NAME_GRAMMAR, False)
    assert li.lipid_annotation_class("Oleic acid", "Fatty acids and conjugates", "null") == (
        True, "Fatty acids and conjugates", li.ROUTE_NAME_GRAMMAR, False)


def test_splash_internal_standard_names_with_dash_separators_are_lipids():
    # 内部標準の鎖の区切りの `-` は `/` に直して読む（同位体の標識は外す）
    for name in ("SM d18:1-18:1(d9)", "Cer d18:1-15:0(d7)", "TG 15:0-18:1(d7)-15:0", "LPC 18:1(d7)"):
        assert li.lipid_name_parses(name) is True, name
    assert li.lipid_annotation_class("SM d18:1-18:1(d9)", "Sphingomyelins") == (True, "SM", li.ROUTE_NAME, False)
    assert li.lipid_annotation_class("Cer d18:1-15:0(d7)", "Long-chain ceramides") == (
        True, "Cer_NS", li.ROUTE_NAME, False)
    # エーテルの `O-` は鎖の区切りではない
    assert li._CHAIN_DASH_RE.sub(r"\1/", "PE O-16:0/18:1") == "PE O-16:0/18:1"


def test_may_move_lipid_class_is_a_cheap_superset_of_the_moves():
    assert li.may_move_lipid_class("MassSpecGymID0000001", "Triacylglycerols") is True
    assert li.may_move_lipid_class("PC 34:1", "null") is True
    assert li.may_move_lipid_class("PC 34:1", "PC") is False           # MS-DIAL の脂質クラスは移さない
    assert li.may_move_lipid_class("Oleic acid", "Fatty acids and conjugates") is False
    assert li.lipid_annotation_class("PC 34:1", "null") == (True, "PC", li.ROUTE_NAME, False)


def test_annotation_name_presence():
    assert all(li.has_annotation_name(n) for n in ("PC 34:1", "Ginsenoside Rb1", "no MS2: FA 5:0"))
    assert not any(li.has_annotation_name(n) for n in (None, "", "  ", "Unknown", "unknown", "null",
                                                        "low score: Unknown", 3))
    # 名前の無いものは読めない（呼び出し側が先に分ける）
    assert li.lipid_name_parses("Unknown") is False and li.lipid_name_parses(None) is False


def test_lipid_name_parses_is_none_without_pygoslin(monkeypatch):
    def missing(_text):
        raise ImportError("no pygoslin")

    monkeypatch.setattr(li, "_parse_lipid", missing)
    li.lipid_name_parses.cache_clear()
    li._lipid_call.cache_clear()
    try:
        assert li.lipid_name_parses("Ginsenoside Rb1") is None
        # 判定できない（規則 1 の MS-DIAL の脂質クラスだけは GOSLIN なしで脂質）
        assert li.lipid_annotation_class("MassSpecGymID0000001", "Triacylglycerols", "C57H104O6").is_lipid is None
        assert li.lipid_annotation_class("PC 34:1", "PC").is_lipid is True
    finally:
        li.lipid_name_parses.cache_clear()
        li._lipid_call.cache_clear()


if __name__ == "__main__":
    unittest.main()
