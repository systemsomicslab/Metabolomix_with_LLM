from metabolomix.pathway import msp_lipid as ml
from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES


def test_parse_formula():
    assert ml.parse_formula("C57H104O6") == {"C": 57, "H": 104, "O": 6}
    assert ml.parse_formula("C39H79N2O6P") == {"C": 39, "H": 79, "N": 2, "O": 6, "P": 1}
    for bad in (None, "", "C39H67[2H]7NO10P", "C6H12O6+", "abc"):
        assert ml.parse_formula(bad) is None


def test_from_class_formula_examples():
    assert ml.from_class_formula("Triacylglycerols", "C57H104O6") == ("TG 54:3", "TG")
    assert ml.from_class_formula("1,2-diacylglycerols", "C37H70O5") == ("DG 34:1", "DG")
    assert ml.from_class_formula("Phosphatidylcholines", "C36H72NO8P") == ("PC 28:0", "PC")
    assert ml.from_class_formula("Very long-chain fatty acids", "C22H32O2") == ("FA 22:6", "FA")
    assert ml.from_class_formula("1-(1Z-alkenyl),2-acylglycerophosphoethanolamines", "C37H72NO7P") == ("PE P-32:1", "EtherPE")
    assert ml.from_class_formula("1-alkyl,2-acylglycero-3-phosphocholines", "C42H84NO7P") == ("PC O-34:1", "EtherPC")
    assert ml.from_class_formula("Sphingomyelins", "C39H77N2O6P") == ("SM 34:2;O2", "SM")
    assert ml.from_class_formula("  Triacylglycerols ", "C57H104O6") == ("TG 54:3", "TG")   # 前後の空白


def test_from_class_formula_reasons():
    assert ml.from_class_formula("Benzoic acid esters", "C16H22O4") == "unmapped_ontology"
    assert ml.from_class_formula("Ceramides", "C34H67NO3") == "ambiguous_ontology"
    assert ml.from_class_formula("Triacylglycerols", None) == "no_formula"
    assert ml.from_class_formula("Lineolic acids and derivatives", "C37H68O5") == "formula_mismatch"   # DG は FA の形でない
    assert ml.from_class_formula("Phosphatidylserines", "C39H67[2H]7NO10P") == "formula_mismatch"
    assert ml.from_class_formula("Triacylglycerols", "C57H103O6") == "formula_mismatch"               # 二重結合が整数でない


def test_tables_are_consistent():
    for key, form in ml.CLASS_FORMS.items():
        assert form.lbm_class in MSDIAL_LIPID_CLASSES, key
        assert ml.from_class_formula_key(key, form.ref_formula) == (
            f"{form.prefix}{form.ref_n}:{form.ref_db}{form.suffix}", form.lbm_class)
    assert set(ml.CLASSYFIRE_TO_CLASS.values()) <= set(ml.CLASS_FORMS)
    assert not set(ml.CLASSYFIRE_TO_CLASS) & ml.AMBIGUOUS_CLASSYFIRE


def test_independent_known_species():
    # 基準の種とは別の、既知の種（基準の式を基準自身に通す一貫性の検査では捕まらない誤りを見る）
    cases = [
        ("Phosphatidylcholines", "C44H84NO8P", "PC 36:2", "PC"),                       # DOPC
        ("Phosphatidylethanolamines", "C41H78NO8P", "PE 36:2", "PE"),                  # DOPE
        ("Phosphatidylserines", "C42H80NO10P", "PS 36:1", "PS"),                       # SOPS
        ("Phosphatidylglycerols", "C40H77O10P", "PG 34:1", "PG"),                      # POPG
        ("Phosphatidylglycerols", "C42H79O10P", "PG 36:2", "PG"),                      # DOPG
        ("1-phosphatidyl-1D-myo-inositols", "C47H83O13P", "PI 38:4", "PI"),            # SAPI
        ("1,2-diacylglycerol-3-phosphates", "C39H73O8P", "PA 36:2", "PA"),             # DOPA
        ("1-acyl-sn-glycero-3-phosphocholines", "C26H52NO7P", "LPC 18:1", "LPC"),
        ("1-acyl-sn-glycero-3-phosphoethanolamines", "C23H46NO7P", "LPE 18:1", "LPE"),
        ("1-(1Z-alkenyl),2-acylglycerophosphoethanolamines", "C43H78NO7P", "PE P-38:4", "EtherPE"),
        ("1-alkyl,2-acylglycerophosphoethanolamines", "C43H80NO7P", "PE O-38:4", "EtherPE"),
        ("1-(1Z-alkenyl),2-acyl-glycerophosphocholines", "C44H80NO7P", "PC P-36:4", "EtherPC"),
        ("1-alkyl,2-acylglycero-3-phosphocholines", "C44H82NO7P", "PC O-36:4", "EtherPC"),
        ("Triacylglycerols", "C55H102O6", "TG 52:2", "TG"),                            # POO
        ("1,2-diacylglycerols", "C39H72O5", "DG 36:2", "DG"),                          # DOG
        ("1-monoacylglycerols", "C19H38O4", "MG 16:0", "MG"),
        ("Long-chain fatty acids", "C16H32O2", "FA 16:0", "FA"),                       # パルミチン酸
        ("Long-chain fatty acids", "C20H32O2", "FA 20:4", "FA"),                       # アラキドン酸
        ("Acyl carnitines", "C25H47NO4", "CAR 18:1", "CAR"),                           # オレオイルカルニチン
        ("N-acylethanolamines", "C20H39NO2", "NAE 18:1", "NAE"),                       # オレオイルエタノールアミド
        ("N-acylethanolamines", "C22H37NO2", "NAE 20:4", "NAE"),                       # アナンダミド
        ("Sphingomyelins", "C41H83N2O6P", "SM 36:1;O2", "SM"),
    ]
    for ontology, formula, name, lbm_class in cases:
        assert ml.from_class_formula(ontology, formula) == (name, lbm_class), (ontology, formula)
