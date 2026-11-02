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
