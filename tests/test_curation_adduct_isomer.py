"""別アダクトの取り違え（`adduct_isomer_of` / `adduct_isomer_minor_of`）の純関数。

数値は 20260930_EV の実例（MS-DIAL Console 出力）の MassCenter / RT / Formula / HeightAverage を
写したもの。実データは読まない（fixture はテスト自身が作る規約）。
"""
import pytest

from metabolomix.curation import adduct_isomer, judge
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

TH = judge.resolve_thresholds(None)


def rep(*, msms=True, matched=True, chains=False):
    return {"has_msms": msms, "is_reference_matched": matched, "is_lipid_chains_match": chains}


def spot(spot_id, name, adduct, mz, rt, formula, ontology, height, match=None):
    """`.arf2` のカタログ行と `load_spot_annotations` の 1 件を、プールの 1 要素にする。"""
    catalog = {"MasterAlignmentID": spot_id, "Name": name, "AdductType": adduct, "MassCenter": mz,
               "RT": rt, "Formula": formula, "Ontology": ontology, "HeightAverage": height}
    return adduct_isomer.pool_entry(catalog, {"representative": match})


CHAIN_CONFIRMED = rep(chains=True)

# neg #173 / #189
PI_173 = spot(173, "PI 41:2", "[M-H]-", 931.6355, 4.3637, "C50H93O13P", "PI", 5702.2, rep())
DGDG_189 = spot(189, "DGDG 35:1|DGDG 16:0_19:1", "[M+CH3COO]-", 991.65778, 4.3706,
                "C50H92O15", "DGDG", 58607.5, CHAIN_CONFIRMED)
# pos #660 / #657
PI_660 = spot(660, "no MS2: PI 41:2", "[M+Na]+", 955.63183, 4.323, "C50H93O13P", "PI", 1472.1,
              rep(msms=False, matched=False))
DGDG_657 = spot(657, "DGDG 35:1|DGDG 16:0_19:1", "[M+NH4]+", 950.67661, 4.345,
                "C50H92O15", "DGDG", 114658.3, CHAIN_CONFIRMED)
# pos #470 / #485
SM_470 = spot(470, "no MS2: SM 35:7;O3", "[M+H]+", 721.50146, 3.132, "C40H69N2O7P", "SM", 1655.3,
              rep(msms=False, matched=False))
PG_485 = spot(485, "PG 32:1|PG 16:0_16:1", "[M+NH4]+", 738.52852, 3.132, "C38H73O10P", "PG",
              5070.8, CHAIN_CONFIRMED)


def _fa18_mz(adduct):
    return round(mz_from_neutral(monoisotopic_mass(parse_formula("C18H36O2")), parse_adduct(adduct)), 5)


# neg #97 / #42: FA 18:0 の二量体 [2M-H]- を別物質として注釈したもの。両方とも参照一致（同じ段階）。
RIKEN_97 = spot(97, "RIKEN N-VS1 ID-2804", "[M-H]-", _fa18_mz("[2M-H]-"), 2.566, "C36H70O4",
                "Fatty esters", 1713.0, rep())
FA_42 = spot(42, "FA 18:0", "[M-H]-", _fa18_mz("[M-H]-"), 2.565, "C18H36O2", "FA", 71201.0, rep())


@pytest.mark.parametrize("target,partner,as_adduct", [
    (PI_173, DGDG_189, "[M-H]-"),
    (PI_660, DGDG_657, "[M+Na]+"),
    (SM_470, PG_485, "[M+H]+"),
])
def test_real_cases_are_explained_by_the_stronger_partner(target, partner, as_adduct):
    hit = adduct_isomer.find_adduct_isomer(target, [target, partner], TH)
    assert hit["of"] == partner["spot_id"]
    assert hit["as_adduct"] == as_adduct
    assert hit["severity"] == "strong"
    assert abs(hit["ppm"]) <= TH["adduct_isomer_ppm"]


def test_the_stronger_partner_is_not_flagged_in_return():
    assert adduct_isomer.find_adduct_isomer(DGDG_189, [PI_173, DGDG_189], TH) is None


def test_same_tier_partner_needs_the_height_ratio():
    hit = adduct_isomer.find_adduct_isomer(RIKEN_97, [RIKEN_97, FA_42], TH)
    assert (hit["of"], hit["as_adduct"], hit["severity"]) == (42, "[2M-H]-", "minor")
    assert hit["height_ratio"] == pytest.approx(71201.0 / 1713.0, rel=1e-3)


def test_same_tier_partner_of_similar_height_does_not_fire():
    similar = {**FA_42, "height": RIKEN_97["height"] * 1.04}
    assert adduct_isomer.find_adduct_isomer(RIKEN_97, [RIKEN_97, similar], TH) is None


def test_rt_outside_the_window_does_not_fire():
    far = {**DGDG_189, "rt": PI_173["rt"] + TH["adduct_isomer_drt"] + 0.01}
    assert adduct_isomer.find_adduct_isomer(PI_173, [PI_173, far], TH) is None


def test_mz_outside_the_ppm_window_does_not_fire():
    shifted = {**PI_173, "mz": PI_173["mz"] * (1 + 8e-6)}   # 実例は -1.2 ppm → +6.8 ppm
    assert adduct_isomer.find_adduct_isomer(shifted, [shifted, DGDG_189], TH) is None


def test_same_formula_and_class_is_the_same_species_not_a_misannotation():
    # DGDG 35:1 の [M-H]- を正しく注釈したスポット（組成式もクラスも同じ）は対象外。
    same = {**PI_173, "name": "DGDG 35:1", "formula": "C50H92O15", "ontology": "DGDG"}
    assert adduct_isomer.find_adduct_isomer(same, [same, DGDG_189], TH) is None


def test_partner_of_opposite_polarity_is_ignored():
    positive = spot(189, "DGDG 35:1|DGDG 16:0_19:1", "[M+NH4]+", 950.67661, 4.3706,
                    "C50H92O15", "DGDG", 58607.5, CHAIN_CONFIRMED)
    assert adduct_isomer.find_adduct_isomer(PI_173, [PI_173, positive], TH) is None


def test_partner_under_its_own_adduct_is_not_a_different_adduct():
    # X が Y と同じアダクトで同じ m/z なら重複スポットで、この理由の対象ではない。
    twin = {**DGDG_189, "spot_id": 190, "adduct": "[M-H]-", "mz": PI_173["mz"]}
    assert adduct_isomer.find_adduct_isomer(PI_173, [PI_173, twin], TH) is None


def test_strongest_partner_wins_when_several_explain_the_spot():
    weaker = {**DGDG_657, "spot_id": 300, "tier": 2, "height": 10 ** 6}
    hit = adduct_isomer.find_adduct_isomer(PI_660, [PI_660, weaker, DGDG_657], TH)
    assert hit["of"] == 657


def test_thresholds_can_be_overridden():
    th = judge.resolve_thresholds({"adduct_isomer_drt": 0.001})
    assert adduct_isomer.find_adduct_isomer(PI_173, [PI_173, DGDG_189], th) is None


@pytest.mark.parametrize("name,match,tier", [
    ("PG 32:1|PG 16:0_16:1", rep(chains=True), 3),
    ("PG 32:1", rep(chains=True), 2),                    # 鎖レベルの名前でなければ 3 にしない
    ("PI 41:2", rep(), 2),
    ("low score: PI 41:2", rep(matched=False), 1),
    ("no MS2: PI 41:2", rep(msms=False, matched=False), 0),
    ("no MS2: PI 41:2", rep(msms=True, matched=False), 0),  # 名前接頭辞が MS2 無しを言う
    ("PI 41:2", None, 0),
])
def test_evidence_tier(name, match, tier):
    assert adduct_isomer.evidence_tier(name, match) == tier


def test_pool_skips_unannotated_spots_and_unparseable_adducts():
    rows = [({"MasterAlignmentID": 1, "Name": "Unknown", "AdductType": "[M-H]-", "MassCenter": 1.0,
              "RT": 1.0, "Formula": "", "Ontology": "", "HeightAverage": 1.0}, {}),
            ({"MasterAlignmentID": 2, "Name": "PC 34:1", "AdductType": "bogus", "MassCenter": 1.0,
              "RT": 1.0, "Formula": "", "Ontology": "", "HeightAverage": 1.0}, {})]
    assert [adduct_isomer.pool_entry(c, a) for c, a in rows] == [None, None]
