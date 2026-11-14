import math
from itertools import product

import numpy as np
import pytest
from scipy.stats import multinomial

from metabolomix.analysis.feature_search import DEFAULT_ADDUCT_TIERS
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


#: 整数質量のずれ → 存在比（IUPAC の代表値）。実装の表（curation.isotope.ABUNDANCES）とは別にここで書き、
#: 原子を 1 個ずつ畳み込む（実装の二乗法を使わない）独立の計算で期待比を確かめる。
IUPAC = {"H": [0.999885, 0.000115], "C": [0.9893, 0.0107], "N": [0.99636, 0.00364],
         "O": [0.99757, 0.00038, 0.00205], "S": [0.9499, 0.0075, 0.0425], "P": [1.0],
         "Cl": [0.7576, 0.0, 0.2424], "Br": [0.5069, 0.0, 0.4931]}


def _brute_ratios(ion: dict) -> tuple[float, float]:
    """イオンの組成（手で書く）から (M+1 / M, M+2 / M)。36S（M+4）は M..M+2 の比に効かない。"""
    dist = np.array([1.0])
    for element, n in ion.items():
        for _ in range(n):
            dist = np.convolve(dist, IUPAC[element])[:3]
    return dist[1] / dist[0], dist[2] / dist[0]


#: 同位体の原子質量（NIST の Atomic Weights and Isotopic Compositions）と存在比。実装の表（質量差の表）とは
#: 別にここで書き、多項分布（scipy）で元素ごとの同位体の組み合わせを数え上げる独立の計算で、窓に入る同位体の
#: 期待比を確かめる。
NIST = {
    "H": [(1.00782503223, 0.999885), (2.01410177812, 0.000115)],
    "C": [(12.0, 0.9893), (13.00335483507, 0.0107)],
    "N": [(14.00307400443, 0.99636), (15.00010889888, 0.00364)],
    "O": [(15.99491461957, 0.99757), (16.99913175650, 0.00038), (17.99915961286, 0.00205)],
    "S": [(31.9720711744, 0.9499), (32.9714589098, 0.0075), (33.967867004, 0.0425), (35.96708071, 0.0001)],
    "Cl": [(34.968852682, 0.7576), (36.965902602, 0.2424)],
    "P": [(30.97376199842, 1.0)],
    "K": [(38.9637064864, 0.932581), (39.963998166, 0.000117), (40.9618252579, 0.067302)],
    "I": [(126.904473, 1.0)],
}

#: 窓（mz_tol）に入る範囲で、Y からのずれ dmz の X に入る同位体だけの和と比べる（`_fine(... mz_tol=WIDE)` は M+k の和）
WIDE = 1.0


def _fine(ion: dict, k: int, dmz: float, mz_tol: float, charge: int = 1) -> float:
    """イオンの組成の M+k の同位体の組み合わせのうち、Y からの m/z のずれ（質量差 / 価数）が dmz から mz_tol 以内に
    来るものの、M に対する比の和。"""
    per_element = []
    for element, n in ion.items():
        isotopes = NIST[element]
        light = isotopes[0][0]
        options = []
        for heavy_counts in product(range(min(n, 2) + 1), repeat=len(isotopes) - 1):
            nominal = sum(c * round(isotopes[i + 1][0] - light) for i, c in enumerate(heavy_counts))
            if sum(heavy_counts) > n or nominal > 2:
                continue
            p = multinomial.pmf([n - sum(heavy_counts), *heavy_counts], n, [a for _, a in isotopes])
            shift = sum(c * (isotopes[i + 1][0] - light) for i, c in enumerate(heavy_counts))
            options.append((nominal, shift, p))
        per_element.append(options)
    mono, total = None, 0.0
    for combo in product(*per_element):
        nominal = sum(o[0] for o in combo)
        p = math.prod(o[2] for o in combo)
        shift = sum(o[1] for o in combo)
        if nominal == 0:
            mono = p
        elif nominal == k and abs(dmz - shift / charge) <= mz_tol:
            total += p
    return float(total / mono)


def test_expected_ratios_of_pure_carbon_are_binomial():
    q = 0.0107 / 0.9893
    assert relations.expected_isotope_ratio("C43", None, 1, 1.003355, WIDE) == pytest.approx(43 * q, rel=1e-9)
    assert relations.expected_isotope_ratio("C43", None, 2, 2.006710, WIDE) == pytest.approx(
        math.comb(43, 2) * q ** 2, rel=1e-9)


@pytest.mark.parametrize("formula, adduct, ion", [
    ("C26H45NO7S", "[M-H]-", {"C": 26, "H": 44, "N": 1, "O": 7, "S": 1}),      # タウロコール酸
    ("C10H17N3O6S", "[M-H]-", {"C": 10, "H": 16, "N": 3, "O": 6, "S": 1}),     # グルタチオン
    ("C6H12O6", "[M-H]-", {"C": 6, "H": 11, "O": 6}),                           # ヘキソース
    ("C20H40O2", "[2M+H]+", {"C": 40, "H": 81, "O": 4}),                        # 二量体は 2 倍
    ("C20H40O2", "[M+Cl]-", {"C": 20, "H": 40, "O": 2, "Cl": 1}),              # 付加体の Cl
    ("C42H80NO8P", "[M+HCOO]-", {"C": 43, "H": 81, "N": 1, "O": 10, "P": 1}),  # ギ酸の原子
    ("C6H5Br", "[M-H]-", {"C": 6, "H": 4, "Br": 1}),
    ("C4H9NO2", "[M+H]1+", {"C": 4, "H": 10, "N": 1, "O": 2}),                 # mzTab-M の表記
])
def test_expected_ratios_are_the_sums_of_the_ion_envelope(formula, adduct, ion):
    """窓が M+k の全体を覆えば、期待比はイオンの同位体分布の M+1 の和・M+2 の和（原子を 1 個ずつ畳み込む計算）。"""
    assert relations.ion_element_counts(formula, adduct) == ion
    charge = parse_adduct(adduct).charge
    m1 = relations.expected_isotope_ratio(formula, adduct, 1, 1.0 / charge, WIDE)
    m2 = relations.expected_isotope_ratio(formula, adduct, 2, 2.0 / charge, WIDE)
    b1, b2 = _brute_ratios(ion)
    assert m1 == pytest.approx(b1, rel=1e-9) and m2 == pytest.approx(b2, rel=1e-9)


def test_taurocholate_m_plus_2_is_about_13c2_plus_34s_plus_18o():
    # 手計算の一次の項: 13C2 C(26,2)·q_C² + 34S 1·q_S + 18O 7·q_O（交差項 13C·15N などが残りの約 5 %）
    q_c, q_s, q_o = 0.0107 / 0.9893, 0.0425 / 0.9499, 0.00205 / 0.99757
    first_order = math.comb(26, 2) * q_c ** 2 + q_s + 7 * q_o
    m2 = relations.expected_isotope_ratio("C26H45NO7S", "[M-H]-", 2, 2.0013, WIDE)
    assert first_order == pytest.approx(0.0971, abs=1e-4)
    assert first_order < m2 < first_order * 1.07


@pytest.mark.parametrize("formula", [None, "", "Unknown", "None", "C10Xq3", "C42H82NO8P (PC)"])
def test_unreadable_formula_has_no_counts_and_no_expected_ratio(formula):
    assert relations.expected_isotope_ratio(formula, "[M-H]-", 1, 1.003355, WIDE) is None
    assert relations.ion_element_counts(formula, "[M-H]-") is None


def test_formula_with_an_element_without_isotope_data_keeps_its_counts():
    """セレノメチオニンの Se は同位体の表に無い（最も多い 80Se の下に 74Se〜78Se がある）。期待比は求めないが、
    組成式は読めるので元素の数は返す（Cl・Br・S の説明の門に使う）。"""
    assert relations.ion_element_counts("C5H11NO2Se", "[M+H]+") == {"C": 5, "H": 12, "N": 1, "O": 2, "Se": 1}
    assert relations.expected_isotope_ratio("C5H11NO2Se", "[M+H]+", 1, 1.003355, WIDE) is None



def test_formula_with_a_trailing_charge_is_read():
    # 第四級アンモニウムの組成式はライブラリによって電荷の記号を付ける（コリン `C5H14NO+`）
    assert relations.ion_element_counts("C5H14NO+", "[M]+") == {"C": 5, "H": 14, "N": 1, "O": 1}


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


def test_chlorine_m_plus_2_is_explained_when_the_window_is_narrow():
    found = _m2(CL, 1.997050, 0.32, mz_tol=0.005)
    [rel] = found
    assert rel["relation"] == "isotope_M+2_37Cl" and rel["ambiguous"] is False
    # 期待比は X の窓（±5 mDa）に入る M+2: 37Cl 0.2424/0.7576 と 13C+15N。13C2（9.66 mDa）・18O（7.2 mDa）は窓の外
    expected = _fine({"C": 10, "H": 7, "Cl": 1, "N": 1, "O": 3}, 2, 1.997050, 0.005)
    assert expected == pytest.approx(0.2424 / 0.7576 + 10 * 0.0107 / 0.9893 * 0.00364 / 0.99636, rel=1e-3)
    assert rel["expected_ratio"] == pytest.approx(expected, abs=1e-4)
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_overlapping_m_plus_2_windows_are_ambiguous_not_a_separate_compound():
    """既定の 10 mDa では 13C2 と 37Cl の窓が重なる。一方に決めず（強い説明にしない）、どの行も X の窓に入る
    M+2 の和（37Cl も 13C2 も入る）を期待比にするので、13C2 だけの期待比に合わないことで「別物質」と言わない。"""
    found = {r["relation"]: r for r in _m2(CL, 1.997050, 0.32)}
    # 組成式に S・Br が無いので 34S・81Br は出さない
    assert set(found) == {"isotope_M+2", "isotope_M+2_37Cl"}
    assert all(r["ambiguous"] and not r["strong"] for r in found.values())
    assert found["isotope_M+2"]["expected_ratio"] == found["isotope_M+2_37Cl"]["expected_ratio"]
    assert all(r["isotope_consistent"] is True for r in found.values())
    assert not any(r["informational"] for r in found.values())


def test_overlapping_m_plus_2_too_intense_for_every_explanation_is_informational():
    # 強度比 1.0 は 13C2（炭素 10）にも塩素 1 個（約 0.32、1.5 倍で 0.48）にも合わない → 別物質
    found = _m2(CL, 1.99684, 1.0)
    assert {r["relation"] for r in found} == {"isotope_M+2", "isotope_M+2_37Cl"}
    assert all(r["informational"] and r["isotope_consistent"] is False for r in found)


def test_without_a_formula_all_m_plus_2_explanations_are_listed_and_ambiguous():
    unknown = {**CL, "formula": None}
    found = {r["relation"]: r for r in _m2(unknown, 2.006710, 0.05)}
    # 34S（+1.995796）は 13C2 から 10.9 mDa で窓の外
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


# ---------- イオンの組成で期待比を求める（二量体・付加体の原子。レビュー R8） ----------

HEIGHTS = [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0]


def _ion_y(formula, adduct, mz=None):
    """注釈付きの Y（組成式と付加体から m/z を計算する。元素質量の表に無い元素を含むなら `mz` で渡す）。"""
    if mz is None:
        mz = mz_from_neutral(monoisotopic_mass(parse_formula(formula)), parse_adduct(adduct))
    return {"spot_id": 1, "name": "Y", "mz": mz, "rt": 10.0, "adduct": adduct, "formula": formula,
            "heights": list(HEIGHTS)}


def test_dimer_m_plus_1_uses_twice_the_carbons():
    """[2M+H]+ の M+1 は炭素 40 個ぶん。中性の組成式（炭素 20）で見ると強すぎて別物質に見えていた。"""
    found = _m2(_ion_y("C20H40O2", "[2M+H]+"), 1.003355, 0.437, adducts=[])
    [rel] = [r for r in found if r["relation"] == "isotope_M+1"]
    assert rel["expected_ratio"] == pytest.approx(_fine({"C": 40, "H": 81, "O": 4}, 1, 1.003355, 0.010), abs=1e-4)
    assert rel["expected_ratio"] == pytest.approx(_brute_ratios({"C": 40, "H": 81, "O": 4})[0], abs=1e-4)
    assert 0.437 > 1.5 * _brute_ratios({"C": 20, "H": 41, "O": 2})[0]      # 中性の炭素数なら別物質
    assert rel["isotope_consistent"] is True and rel["informational"] is False and rel["strong"] is True


def test_chloride_adduct_m_plus_2_is_the_37cl_isotope():
    """[M+Cl]- のイオンは Cl を持つので 37Cl の説明が出る（中性の組成式に Cl が無くても）。"""
    found = {r["relation"]: r for r in _m2(_ion_y("C20H40O2", "[M+Cl]-"), 1.997050, 0.34, adducts=[])}
    assert set(found) == {"isotope_M+2", "isotope_M+2_37Cl"}             # 13C2 は 9.66 mDa で窓の内
    expected = _fine({"C": 20, "H": 40, "O": 2, "Cl": 1}, 2, 1.997050, 0.010)
    for rel in found.values():
        assert rel["expected_ratio"] == pytest.approx(expected, abs=1e-4)
        assert rel["isotope_consistent"] is True and rel["informational"] is False
        assert rel["ambiguous"] is True and rel["strong"] is False
    # 同じ分子の [M-H]- は Cl を持たない: 37Cl の説明は出ず、強度比 0.34 は M+2 には強すぎる（別物質）
    [rel] = _m2(_ion_y("C20H40O2", "[M-H]-"), 1.997050, 0.34, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["informational"] is True


# ---------- 34S・18O の M+2（レビュー R9） ----------

def test_unresolved_m_plus_2_of_a_sulfur_compound_is_not_a_separate_compound():
    """タウロコール酸 [M-H]- の分離できない M+2（13C2・34S・18O の束。重心 +2.0013、強度比 0.097）。
    13C2 だけの期待比（約 0.038）では 1.5 倍を超えて別物質に見えていた。"""
    found = {r["relation"]: r for r in _m2(_ion_y("C26H45NO7S", "[M-H]-"), 2.0013, 0.097, adducts=[])}
    assert set(found) == {"isotope_M+2", "isotope_M+2_34S"}              # Cl・Br は出さない
    expected = _fine({"C": 26, "H": 44, "N": 1, "O": 7, "S": 1}, 2, 2.0013, 0.010)
    assert expected == pytest.approx(_brute_ratios({"C": 26, "H": 44, "N": 1, "O": 7, "S": 1})[1], rel=0.01)
    assert 0.097 > 1.5 * math.comb(26, 2) * (0.0107 / 0.9893) ** 2       # 13C2 だけなら別物質
    for rel in found.values():
        assert rel["expected_ratio"] == pytest.approx(expected, abs=1e-4)
        assert rel["isotope_consistent"] is True and rel["informational"] is False
        assert rel["ambiguous"] is True and rel["strong"] is False


def test_resolved_34s_peak_is_its_own_explanation():
    """分離された 34S（+1.995796）は 13C2 から 10.9 mDa で窓の外。34S の説明が 1 つだけ立つ。"""
    [rel] = _m2(_ion_y("C26H45NO7S", "[M-H]-"), 1.995796, 0.045, adducts=[])
    assert rel["relation"] == "isotope_M+2_34S" and rel["ambiguous"] is False
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_sulfur_isotope_without_a_formula_is_not_called_a_halogen():
    unknown = {**_ion_y("C26H45NO7S", "[M-H]-"), "formula": None}
    found = {r["relation"]: r for r in _m2(unknown, 1.995796, 0.045, adducts=[])}
    assert set(found) == {"isotope_M+2_34S", "isotope_M+2_37Cl", "isotope_M+2_81Br"}
    assert all(r["ambiguous"] and r["isotope_consistent"] is None for r in found.values())


def test_oxygen_rich_m_plus_2_counts_18o():
    """ヘキソース [M-H]- の M+2 は 18O（6 個）が主。13C2 だけの期待比では別物質に見えていた。"""
    [rel] = _m2(_ion_y("C6H12O6", "[M-H]-"), 2.0046, 0.0145, adducts=[])
    assert rel["relation"] == "isotope_M+2"
    assert 0.0145 > 1.5 * math.comb(6, 2) * (0.0107 / 0.9893) ** 2
    assert rel["expected_ratio"] == pytest.approx(_fine({"C": 6, "H": 11, "O": 6}, 2, 2.0046, 0.010), abs=1e-4)
    assert rel["isotope_consistent"] is True and rel["informational"] is False and rel["strong"] is True


def test_resolved_18o_peak_is_its_own_explanation_when_the_window_is_narrow():
    """窓が 13C2 と 18O の差（2.46 mDa）より狭ければ、分離された 18O（+2.004245）は 13C2 の窓の外。
    18O の説明が 1 つだけ立ち、期待比はその窓に入る 18O（6 個）だけ。"""
    y = _ion_y("C6H12O6", "[M-H]-")
    [rel] = _m2(y, 2.004246, 0.0145, adducts=[], mz_tol=0.002)
    assert rel["relation"] == "isotope_M+2_18O" and rel["ambiguous"] is False
    expected = _fine({"C": 6, "H": 11, "O": 6}, 2, 2.004246, 0.002)
    assert expected == pytest.approx(6 * 0.00205 / 0.99757, rel=1e-3)
    assert rel["expected_ratio"] == pytest.approx(expected, abs=1e-4)
    assert rel["isotope_consistent"] is True and rel["strong"] is True
    # 既定の 10 mDa では 13C2 と分けられないので、13C2 の行（isotope_M+2）だけ（18O の行と並べて曖昧にしない）
    [rel] = _m2(y, 2.004246, 0.0145, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False and rel["strong"] is True


def test_18o_needs_oxygen_and_follows_the_charge():
    # O を持たないイオン: 狭い窓で 18O の位置にある X は説明できない
    assert _m2(_ion_y("C10H16N2", "[M+H]+"), 2.004246, 0.01, adducts=[], mz_tol=0.002) == []
    # 2 価では 13C2 と 18O の差は 1.23 mDa。2 mDa の窓では分けられず 13C2 の行になる
    y = _ion_y("C42H80NO10P", "[M-2H]2-")
    [rel] = _m2(y, 2.004246 / 2, 0.05, adducts=[], mz_tol=0.002)
    assert rel["relation"] == "isotope_M+2"


def test_lipid_without_sulfur_has_no_34s_explanation():
    found = relations.find_relations(_x(Y_MZ + 1.995796 + 0.006), [Y], links={}, **KW)
    assert [r["relation"] for r in found] == ["isotope_M+2"]


# ---------- Y 自身の付加体は組成で見分ける（レビュー R11） ----------

POSITIVE = [a for tier in DEFAULT_ADDUCT_TIERS["positive"] for a in tier]
NEGATIVE = [a for tier in DEFAULT_ADDUCT_TIERS["negative"] for a in tier]


def _isobar(y, delta=0.001):
    """Y と同時に溶出し、強度が相関する、Y と同じ m/z（+delta）の別の特徴。"""
    return {"spot_id": 9, "name": None, "mz": y["mz"] + delta, "rt": 10.02, "adduct": None,
            "formula": None, "heights": [h * 0.4 for h in y["heights"]]}


@pytest.mark.parametrize("adduct, adducts", [
    ("[M+H]1+", POSITIVE), ("[M+H]+", POSITIVE),
    ("[M+HCOO]1-", NEGATIVE), ("[M+HCOO]-", NEGATIVE), ("[M+CHO2]1-", NEGATIVE),
    ("[M+CH3COO]1-", NEGATIVE), ("[M+Hac-H]-", NEGATIVE),
])
def test_coeluting_isobar_is_not_an_adduct_of_its_own_ion(adduct, adducts):
    y = _ion_y("C42H82NO8P", adduct)
    assert relations.find_relations(_isobar(y), [y], links={}, **{**KW, "adducts": adducts}) == []


def test_adduct_with_the_same_mz_as_y_is_skipped_but_a_real_partner_is_kept():
    y = _ion_y("C42H82NO8P", "[M+H]1+")
    adducts = ["[2M+2H]2+", "[M+Na]+"]                  # [2M+2H]2+ は [M+H]+ と同じ m/z
    assert relations.find_relations(_isobar(y), [y], links={}, **{**KW, "adducts": adducts}) == []
    sodium = {**_isobar(y), "mz": mz_from_neutral(monoisotopic_mass(parse_formula("C42H82NO8P")),
                                                  parse_adduct("[M+Na]+"))}
    [rel] = relations.find_relations(sodium, [y], links={}, **{**KW, "adducts": adducts})
    assert rel["relation"] == "adduct:[M+Na]+/[M+H]1+" and rel["strong"] is True


# ---------- 分離された 1 つの M+2 はその位置の窓に入る同位体だけと比べる（レビュー V1） ----------

def test_resolved_13c2_of_a_chloride_ion_is_judged_against_its_own_position():
    """[M+Cl]- の 13C2 の位置（37Cl から 9.66 mDa）だけに当たる X。5 mDa の窓では 37Cl と分離されているので、
    期待比は 13C2 とその近く（18O・13C+2H など）の和で、37Cl（0.32）を足した M+2 の全体ではない。
    その位置の同位体の 5 倍の X を「同位体で説明できる」（strong）と言っていた。"""
    ion = {"C": 20, "H": 40, "O": 2, "Cl": 1}
    y = _ion_y("C20H40O2", "[M+Cl]-")
    own = _fine(ion, 2, 2.006710, 0.005)
    assert own == pytest.approx(0.0275, abs=5e-4)                         # 13C2 0.0222 + 18O 0.0041 + 交差項
    assert 5 * own < 1.5 * _brute_ratios(ion)[1]                          # M+2 の全体（0.347）なら「同位体」になる
    [rel] = _m2(y, 2.006710, 5 * own, adducts=[], mz_tol=0.005)
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False
    assert rel["expected_ratio"] == pytest.approx(own, abs=1e-4)
    assert rel["isotope_consistent"] is False and rel["informational"] is True and rel["strong"] is False
    [rel] = _m2(y, 2.006710, own, adducts=[], mz_tol=0.005)
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_resolved_34s_of_a_sulfatide_is_judged_against_its_own_position():
    """スルファチド [M-H]- の 34S の位置（13C2 から 10.91 mDa で窓の外）だけに当たる X。期待比は 34S と、窓に入る
    18O（8.45 mDa）・13C+15N・13C+33S などの和で、13C2（炭素 42 個で約 0.10）は含まない。"""
    ion = {"C": 42, "H": 80, "N": 1, "O": 11, "S": 1}
    y = _ion_y("C42H81NO11S", "[M-H]-")
    own = _fine(ion, 2, 1.995796, 0.010)
    assert 0.0425 / 0.9499 < own < 0.08                                   # 34S 0.0447 + 18O 0.0226 + 交差項
    [rel] = _m2(y, 1.995796, 5 * 0.0425 / 0.9499, adducts=[])
    assert rel["relation"] == "isotope_M+2_34S" and rel["ambiguous"] is False
    assert rel["expected_ratio"] == pytest.approx(own, abs=1e-4)
    assert rel["informational"] is True and rel["strong"] is False
    [rel] = _m2(y, 1.995796, 0.0425 / 0.9499, adducts=[])
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_resolved_13c2_of_glutathione_does_not_borrow_the_34s():
    """グルタチオン [M-H]- の 13C2 の位置に当たる X（34S は 10.91 mDa で窓の外）。期待比は約 0.02 で、
    34S（0.045）を足した M+2 の全体（0.065）ではない。"""
    ion = {"C": 10, "H": 16, "N": 3, "O": 6, "S": 1}
    own = _fine(ion, 2, 2.006710, 0.010)
    assert own == pytest.approx(0.020, abs=2e-3)
    [rel] = _m2(_ion_y("C10H17N3O6S", "[M-H]-"), 2.006710, 0.09, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["expected_ratio"] == pytest.approx(own, abs=1e-4)
    assert rel["informational"] is True and rel["strong"] is False


def test_isotopologue_terms_sum_to_the_envelope_and_carry_their_masses():
    from metabolomix.curation import isotope

    ion = {"C": 26, "H": 44, "N": 1, "O": 7, "S": 1, "K": 1, "Cl": 1}
    terms = isotope.isotopologues(ion)
    envelope = isotope.envelope_ratios(ion)
    for k in (1, 2):
        assert sum(r for kk, _, r in terms if kk == k) == pytest.approx(envelope[k], rel=1e-12)
        assert sum(r for kk, _, r in terms if kk == k) == pytest.approx(_fine(ion, k, float(k), WIDE), rel=1e-9)
    # 質量差つきの項: 窓を狭めた和も数え上げと一致する（13C2 の位置 ±3 mDa、34S の位置 ±1 mDa）
    for k, dmz, tol in ((2, 2.006710, 0.003), (2, 1.995796, 0.001), (1, 0.997035, 0.002)):
        window = sum(r for kk, shift, r in terms if kk == k and abs(dmz - shift) <= tol)
        assert window == pytest.approx(_fine(ion, k, dmz, tol), rel=1e-9)


# ---------- 同位体の表に無い元素を含んでも組成式で門を通す（レビュー V2） ----------

def test_iodine_compound_keeps_its_gate_and_expected_ratios():
    """チロキシン C15H11I4NO4（ヨウ素は単一同位体）。組成式が分かっているので Cl・Br の説明は出さず、M+1・M+2 の
    期待比も求まる（ヨウ素を表に無い元素として組成式ごと捨てると、M+2 は 37Cl・81Br と並んで曖昧になっていた）。"""
    ion = {"C": 15, "H": 10, "I": 4, "N": 1, "O": 4}
    y = _ion_y("C15H11I4NO4", "[M-H]-", mz=775.6801)
    [rel] = _m2(y, 1.003355, 0.17, adducts=[])
    assert rel["relation"] == "isotope_M+1"
    assert rel["expected_ratio"] == pytest.approx(_fine(ion, 1, 1.003355, 0.010), abs=1e-4)
    assert rel["isotope_consistent"] is True and rel["strong"] is True
    [rel] = _m2(y, 2.006710, 0.012, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False
    assert rel["expected_ratio"] == pytest.approx(_fine(ion, 2, 2.006710, 0.010), abs=1e-4)
    assert rel["isotope_consistent"] is True and rel["strong"] is True


def test_element_without_isotope_data_still_gates_the_halogen_explanations():
    """Se の同位体は表に無いので期待比は空だが、組成式に Cl・Br が無いことは分かる: 37Cl の位置の X に
    37Cl・81Br の説明を並べない（13C2 の行だけ、期待比は空で強い説明にも別物質にもならない）。"""
    y = _ion_y("C5H11NO2Se", "[M+H]+", mz=198.0028)
    [rel] = _m2(y, 1.997050, 0.3, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False
    assert rel["expected_ratio"] is None and rel["isotope_consistent"] is None
    assert rel["informational"] is False and rel["strong"] is False


# ---------- 41K の M+2（K を持つイオン。最終検証 F30） ----------

def test_resolved_41k_of_a_potassium_adduct_is_its_own_explanation():
    """ヘキソース [M+K]+ の M+2 は 41K（K 1 個で 0.067302 / 0.932581 ≈ 0.0722）が主。窓が 13C2 との差（8.59 mDa）より
    狭ければ 41K の行（`isotope_M+2_41K`）が 1 つだけ立ち、期待比はその窓に入る同位体だけ（以前は行が無く、説明が付かなかった）。"""
    y = _ion_y("C6H12O6", "[M+K]+")
    ion = {"C": 6, "H": 12, "O": 6, "K": 1}
    assert _fine(ion, 2, 1.998119, 0.005) == pytest.approx(0.067302 / 0.932581, rel=0.01)
    for tol in (0.003, 0.005, 0.008):
        [rel] = _m2(y, 1.998119, 0.072, adducts=[], mz_tol=tol)
        assert rel["relation"] == "isotope_M+2_41K" and rel["ambiguous"] is False
        assert rel["dmz_mda"] == pytest.approx(0.0, abs=0.01)
        assert rel["expected_ratio"] == pytest.approx(_fine(ion, 2, 1.998119, tol), abs=1e-4)
        assert rel["isotope_consistent"] is True and rel["strong"] is True
    # 既定の 10 mDa では 13C2 と分けられないので、13C2 の行（isotope_M+2）が 41K も含む（2 行で曖昧にしない）
    [rel] = _m2(y, 1.998119, 0.072, adducts=[])
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False and rel["strong"] is True
    assert rel["dmz_mda"] == pytest.approx(-8.59, abs=0.01)
    assert rel["expected_ratio"] == pytest.approx(_fine(ion, 2, 1.998119, 0.010), abs=1e-4)


def test_41k_needs_potassium_and_follows_the_charge():
    # K を持たないイオン（[M+H]+）: 狭い窓で 41K の位置にある X は説明できない
    assert _m2(_ion_y("C6H12O6", "[M+H]+"), 1.998119, 0.07, adducts=[], mz_tol=0.005) == []
    # 2 価では 13C2 と 41K の差は 4.30 mDa。5 mDa の窓では分けられず 13C2 の行になる
    y = _ion_y("C42H80NO8P", "[M+2K]2+")
    [rel] = _m2(y, 1.998119 / 2, 0.2, adducts=[], mz_tol=0.005)
    assert rel["relation"] == "isotope_M+2"
    [rel] = _m2(y, 1.998119 / 2, 0.2, adducts=[], mz_tol=0.003)
    assert rel["relation"] == "isotope_M+2_41K"


# ---------- 分離できない M+2 の束と窓の幅（文書の条件を縛る。最終検証 F29） ----------

def test_unresolved_bundle_is_compared_with_the_whole_bundle_only_when_the_window_covers_it():
    """タウロコール酸 [M-H]- の束（重心 +2.0015）が束全体の和と比べられるのは、重心から 13C2（5.21 mDa）と
    34S（5.70 mDa）の両方が窓に入るときだけ。その間の窓では近い 13C2 の 1 行だけに当たり、期待比が 34S を欠くので
    Y 自身の束が別物質（informational）に見える。近い方より狭ければ M+2 の行は出ない（curation.md の期待比の段落）。"""
    ion = {"C": 26, "H": 44, "N": 1, "O": 7, "S": 1}
    y = _ion_y("C26H45NO7S", "[M-H]-")
    total = _fine(ion, 2, 2.0, WIDE)
    assert total == pytest.approx(_brute_ratios(ion)[1], rel=0.01)
    assert _m2(y, 2.0015, total, adducts=[], mz_tol=0.005) == []
    [rel] = _m2(y, 2.0015, total, adducts=[], mz_tol=0.0055)
    assert rel["relation"] == "isotope_M+2" and rel["ambiguous"] is False
    assert rel["expected_ratio"] == pytest.approx(_fine(ion, 2, 2.0015, 0.0055), abs=1e-4)
    assert total > 1.5 * rel["expected_ratio"]
    assert rel["informational"] is True and rel["strong"] is False
    for tol in (0.006, 0.010):
        found = {r["relation"]: r for r in _m2(y, 2.0015, total, adducts=[], mz_tol=tol)}
        assert set(found) == {"isotope_M+2", "isotope_M+2_34S"}
        assert all(r["ambiguous"] and r["isotope_consistent"] is True and not r["informational"]
                   for r in found.values())
    [wide] = [r for r in _m2(y, 2.0015, total, adducts=[]) if r["relation"] == "isotope_M+2"]
    assert wide["expected_ratio"] == pytest.approx(total, abs=1e-4)
