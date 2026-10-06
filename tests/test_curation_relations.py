import pytest

from metabolomix.curation import relations
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

PC342 = monoisotopic_mass(parse_formula("C42H80NO8P"))
Y_MZ = mz_from_neutral(PC342, parse_adduct("[M+HCOO]-"))
Y = {"spot_id": 0, "name": "PC 34:2", "mz": Y_MZ, "rt": 10.00, "adduct": "[M+HCOO]-",
     "formula": "C42H80NO8P", "heights": [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0]}
ADDUCTS = ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Cl]-"]
KW = {"adducts": ADDUCTS, "rt_window": 0.1, "mz_tol": 0.010, "min_r": 0.8}


def _x(mz, *, rt=10.01, heights=None):
    return {"spot_id": 9, "name": None, "mz": mz, "rt": rt, "adduct": None, "formula": None,
            "heights": heights if heights is not None else [h * 0.1 for h in Y["heights"]]}


def test_expected_isotope_ratio_is_binomial():
    assert relations.expected_isotope_ratio(43, 1) == pytest.approx(43 * 0.0107 / 0.9893, rel=1e-6)
    assert relations.expected_isotope_ratio(43, 2) == pytest.approx(903 * (0.0107 / 0.9893) ** 2, rel=1e-6)


def test_m_plus_2_with_a_consistent_ratio_is_a_strong_explanation():
    found = relations.find_relations(_x(Y_MZ + 2.006710), [Y], links={}, **KW)
    iso = next(r for r in found if r["relation"] == "isotope_M+2")
    assert iso["isotope_consistent"] is True and iso["strong"] is True and iso["of"] == 0
    assert iso["r"] == pytest.approx(1.0, abs=1e-2) and iso["candidate_id"] == "R1"   # log1p なので厳密に 1 ではない


def test_m_plus_2_that_is_too_intense_is_informational_only():
    found = relations.find_relations(_x(Y_MZ + 2.006710, heights=list(Y["heights"])), [Y], links={}, **KW)
    iso = next(r for r in found if r["relation"] == "isotope_M+2")
    assert iso["isotope_consistent"] is False and iso["informational"] is True and iso["strong"] is False


def test_adduct_pair_is_explained():
    x_mz = mz_from_neutral(PC342, parse_adduct("[M+CH3COO]-"))
    found = relations.find_relations(_x(x_mz), [Y], links={}, **KW)
    assert [r["relation"] for r in found] == ["adduct:[M+CH3COO]-/[M+HCOO]-"]


def test_insource_demethylation_from_formate():
    found = relations.find_relations(_x(Y_MZ - 60.021129), [Y], links={}, **KW)
    assert [r["relation"] for r in found] == ["insource:-HCOOCH3"]


def test_found_in_upper_msms_link_alone_makes_a_strong_candidate():
    found = relations.find_relations(_x(500.0), [Y], links={9: [{"spot_id": 0, "kind": "found_in_upper_msms"}]}, **KW)
    assert [(r["relation"], r["strong"]) for r in found] == [("found_in_upper_msms", True)]


def test_found_in_upper_msms_link_needs_the_partner_to_be_heavier():
    # MS-DIAL のリンクは無向の対。断片は前駆体より軽いので、X が Y より重いときは関係にしない。
    found = relations.find_relations(_x(900.0), [Y], links={9: [{"spot_id": 0, "kind": "found_in_upper_msms"}]}, **KW)
    assert found == []


def test_not_coeluting_is_not_related():
    assert relations.find_relations(_x(Y_MZ + 2.006710, rt=10.5), [Y], links={}, **KW) == []


def test_unsupported_relation_is_soft():
    found = relations.find_relations(_x(Y_MZ - 60.021129, heights=[5, 1, 4, 2, 6, 3]), [Y], links={}, **KW)
    assert found[0]["soft"] == ["no_support"] and found[0]["strong"] is False


def test_correlation_needs_five_samples_with_both_values():
    assert relations.profile_correlation([1, 2, 3, 4], [1, 2, 3, 4]) is None
    assert relations.profile_correlation([1, 2, None, 4, 5, 6], [2, 4, 6, 8, 10, 12]) == pytest.approx(1.0, abs=2e-2)
