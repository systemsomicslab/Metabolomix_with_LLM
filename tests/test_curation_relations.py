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


# ---------- 付加体・源内断片の強い説明（A4-T15） ----------

def test_correlated_adduct_pair_is_strong_without_links():
    """Δm/z・同時溶出に加えて試料間の相関 r ≥ min_r がそろえば、MS-DIAL のリンクが無くても強い説明。"""
    x_mz = mz_from_neutral(PC342, parse_adduct("[M+CH3COO]-"))
    [rel] = relations.find_relations(_x(x_mz), [Y], links={}, **KW)
    assert rel["r"] >= KW["min_r"] and rel["strong"] is True and rel["soft"] == []


def test_adduct_pair_without_correlation_is_not_strong():
    x_mz = mz_from_neutral(PC342, parse_adduct("[M+CH3COO]-"))
    [rel] = relations.find_relations(_x(x_mz, heights=[5, 1, 4, 2, 6, 3]), [Y], links={}, **KW)
    assert rel["strong"] is False and rel["soft"] == ["no_support"]
    # 相関のリンク（correl_similar）は裏付けにはなるが、付加体の組を強い説明にするのは相関 r だけ
    linked = relations.find_relations(_x(x_mz, heights=[5, 1, 4, 2, 6, 3]), [Y],
                                      links={9: [{"spot_id": 0, "kind": "correl_similar"}]}, **KW)
    assert linked[0]["soft"] == [] and linked[0]["strong"] is False


def test_correlated_insource_fragment_is_strong():
    [rel] = relations.find_relations(_x(Y_MZ - 60.021129), [Y], links={}, **KW)
    assert rel["relation"] == "insource:-HCOOCH3" and rel["strong"] is True


# ---------- 37Cl・81Br の M+2 と、13C2 との取り違え（A5-T16） ----------

CL = {"spot_id": 1, "name": "Cl drug", "mz": 300.0, "rt": 10.0, "adduct": "[M-H]-",
      "formula": "C10H8ClNO3", "heights": [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0]}


def _m2(y, delta, ratio, **kw):
    x = {"spot_id": 9, "name": None, "mz": y["mz"] + delta, "rt": 10.01, "adduct": None, "formula": None,
         "heights": [h * ratio for h in y["heights"]]}
    return relations.find_relations(x, [y], links={}, **{**KW, **kw})


def test_expected_ratio_for_one_chlorine_and_one_bromine():
    assert relations.expected_isotope_ratio(1, 1, relations.ISOTOPE_ABUNDANCE["Cl"]) == pytest.approx(
        0.2424 / 0.7576, rel=1e-6)
    assert relations.expected_isotope_ratio(2, 1, relations.ISOTOPE_ABUNDANCE["Br"]) == pytest.approx(
        2 * 0.4931 / 0.5069, rel=1e-6)


def test_chlorine_m_plus_2_is_explained_when_the_window_is_narrow():
    found = _m2(CL, 1.997050, 0.32, mz_tol=0.005)
    [rel] = found
    assert rel["relation"] == "isotope_M+2_37Cl" and rel["ambiguous"] is False
    assert rel["expected_ratio"] == pytest.approx(0.32, abs=1e-3)
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_overlapping_m_plus_2_windows_are_ambiguous_not_a_separate_compound():
    """既定の 10 mDa では 13C2 と 37Cl の窓が重なる。一方に決めず、13C2 の期待比に合わないことだけで
    「別物質」と言わない（37Cl の期待比には合う）。"""
    found = {r["relation"]: r for r in _m2(CL, 1.997050, 0.32)}
    assert set(found) == {"isotope_M+2", "isotope_M+2_37Cl"}      # 組成式に Br が無いので 81Br は出さない
    assert all(r["ambiguous"] and not r["strong"] for r in found.values())
    assert found["isotope_M+2"]["isotope_consistent"] is False      # 炭素 10 個の 13C2 には強すぎる
    assert found["isotope_M+2_37Cl"]["isotope_consistent"] is True
    assert not any(r["informational"] for r in found.values())


def test_overlapping_m_plus_2_too_intense_for_every_explanation_is_informational():
    # 強度比 1.0 は 13C2（炭素 10）にも塩素 1 個（約 0.32、1.5 倍で 0.48）にも合わない → 別物質
    found = _m2(CL, 1.99684, 1.0)
    assert {r["relation"] for r in found} == {"isotope_M+2", "isotope_M+2_37Cl"}
    assert all(r["informational"] and r["isotope_consistent"] is False for r in found)


def test_without_a_formula_all_m_plus_2_explanations_are_listed_and_ambiguous():
    unknown = {**CL, "formula": None}
    found = {r["relation"]: r for r in _m2(unknown, 2.006710, 0.05)}
    assert set(found) == {"isotope_M+2", "isotope_M+2_37Cl", "isotope_M+2_81Br"}
    assert all(r["ambiguous"] and r["isotope_consistent"] is None and not r["informational"]
               for r in found.values())


def test_bromine_m_plus_2():
    br = {**CL, "formula": "C6H5Br"}
    [rel] = _m2(br, 1.997952, 0.97, mz_tol=0.0005)
    assert rel["relation"] == "isotope_M+2_81Br" and rel["isotope_consistent"] is True
    assert rel["strong"] is True


def test_m_plus_1_is_never_ambiguous():
    [rel] = relations.find_relations(_x(Y_MZ + 1.003355), [Y], links={}, **KW)
    assert rel["relation"] == "isotope_M+1" and rel["ambiguous"] is False


# ---------- 解析フォルダのデータセット（SME 行が無い）でも同位体の強い説明が立つ（A4-T15） ----------

def _folder_dataset():
    """`.arf2` から開いたときと同じ形: SME 行（sme_rows）は無く、組成式と付加体は SML 注釈が持つ。

    特徴 1 は PC 34:2 [M+HCOO]-（注釈あり）、2 はその M+1（強度 0.4 倍）、3 は [M+CH3COO]-（強度 0.2 倍）。
    """
    import numpy as np

    from metabolomix.mztab.dataset_state import DatasetState

    ds = DatasetState()
    ds.feature_ids = ["1", "2", "3"]
    base = np.array([1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0])
    ds.feature_matrix = np.vstack([base, base * 0.4, base * 0.2])
    ds.sample_names = [f"S{i}" for i in range(6)]
    acetate = mz_from_neutral(PC342, parse_adduct("[M+CH3COO]-"))
    ds.feature_metadata = {
        "1": {"mz": Y_MZ, "rt": 10.0, "formula": "C42H80NO8P"},
        "2": {"mz": Y_MZ + 1.003355, "rt": 10.0, "formula": ""},
        "3": {"mz": acetate, "rt": 10.01, "formula": ""},
    }
    ds.feature_candidates = {}
    ds.feature_annotations = {"1": {"sml_id": "arf2:1", "ambiguous": False, "name": "PC 34:2",
                                    "chemical_formula": "C42H80NO8P", "adduct": "[M+HCOO]-"}}
    return ds


def test_isotope_and_adduct_are_strong_on_a_folder_dataset():
    from metabolomix.analysis import claim_evidence

    out = claim_evidence.feature_relations(_folder_dataset(), "1", rt_window=0.05, mz_tol=0.01,
                                           min_r=0.7, adducts=ADDUCTS)
    ions = {e["ion"]: e for e in out["ions_of_this"]}
    assert ions["2"]["relation"] == "isotope_M+1"
    assert ions["2"]["expected_ratio"] is not None                 # 組成式は SML 注釈から
    assert ions["2"]["isotope_consistent"] is True and ions["2"]["strong"] is True
    assert ions["3"]["relation"] == "adduct:[M+CH3COO]-/[M+HCOO]-" and ions["3"]["strong"] is True
    # 未注釈の特徴から見ると、注釈（SML、MS1 だけの同定）の付加体を相手の付加体として使う
    as_ion = claim_evidence.feature_relations(
        _folder_dataset(), "3", rt_window=0.05, mz_tol=0.01, min_r=0.7, adducts=ADDUCTS)["as_ion_of"]
    assert [(e["of"], e["relation"]) for e in as_ion] == [("1", "adduct:[M+CH3COO]-/[M+HCOO]-")]
    assert as_ion[0]["of_name"] == "PC 34:2" and as_ion[0]["strong"] is True
