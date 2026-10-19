# tests/test_feature_search.py
"""feature_search のテスト: 土台（目標の正規化・理論 m/z・規則 hash）と目標照合の判定（spec §5）。"""
import numpy as np
import pytest

from metabolomix.analysis.feature_search import (
    DEFAULT_ADDUCT_TIERS,
    MAX_TARGETS,
    _feature_mz_array,
    classify_identity,
    feature_stats,
    match_target,
    normalize_targets,
    result_summaries,
    rule_hash,
    search_window,
    select_results,
    theoretical_mz,
)
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState

GABA_KEY = "BTCSSZJGUNDROE-UHFFFAOYSA-N"


def _gaba(**extra):
    return {"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", **extra}


def test_theoretical_mz_of_gaba_protonated():
    assert theoretical_mz("C4H9NO2", "[M+H]+") == pytest.approx(104.0706, abs=1e-4)


def test_theoretical_mz_unreadable_returns_none():
    assert theoretical_mz("", "[M+H]+") is None
    assert theoretical_mz("Xx9", "[M+H]+") is None
    # parse_formula は読めない部分を黙って捨てるので、全体が読めることを先に確かめる
    for junk in ("C4H9NO2!!", "Ca(OH)2", "C4H9NO2-3", "C4H9NO2 junk", "c4h9no2"):
        assert theoretical_mz(junk, "[M+H]+") is None, junk
    assert theoretical_mz("C4H9NO2", "[bogus]") is None


def test_default_tiers_follow_dataset_polarity():
    [t] = normalize_targets([_gaba()], dataset_polarity="negative")
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["negative"]
    assert t["inchikey_block"] == "BTCSSZJGUNDROE" and t["ppm"] == 10 and t["rt_tol_min"] == 0.3
    assert t["polarity"] == "negative" and t["applicable"] is True
    # 理論 m/z は全層の全付加体について計算される
    flat = [a for tier in DEFAULT_ADDUCT_TIERS["negative"] for a in tier]
    assert set(t["theoretical_mz"]) == set(flat)
    assert t["theoretical_mz"]["[M-H]-"] == pytest.approx(102.0561, abs=1e-4)


def test_default_tiers_values_are_the_spec_ones():
    assert DEFAULT_ADDUCT_TIERS["positive"] == [
        ["[M+H]+", "[M]+"],
        ["[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M+H-NH3]+", "[2M+H]+", "[M+2H]2+"],
    ]
    assert DEFAULT_ADDUCT_TIERS["negative"] == [
        ["[M-H]-"],
        ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[2M-H]-", "[M+Na-2H]-"],
    ]
    assert MAX_TARGETS == 50


def test_normalized_tiers_are_copies_of_defaults():
    [t] = normalize_targets([_gaba()], dataset_polarity="positive")
    t["adduct_tiers"][0].append("[M+K]+")
    assert DEFAULT_ADDUCT_TIERS["positive"][0] == ["[M+H]+", "[M]+"]


def test_expected_mz_form_fixes_a_single_adduct():
    [t] = normalize_targets(
        [{"target_id": "c1", "inchikey": GABA_KEY.lower(), "expected_mz": 104.0706,
          "adduct": "[M+H]+", "rt_min": 3.2, "name": "GABA"}],
        dataset_polarity="positive")
    assert t["adduct_tiers"] == [["[M+H]+"]]
    assert t["theoretical_mz"] == {"[M+H]+": 104.0706}
    assert t["formula"] is None and t["expected_mz"] == 104.0706
    assert t["inchikey_block"] == "BTCSSZJGUNDROE"
    assert t["rt_min"] == 3.2 and t["name"] == "GABA"


def test_explicit_adduct_tiers_are_kept_and_computed():
    [t] = normalize_targets([_gaba(adduct_tiers=[["[M+H]+"], ["[M+Na]+"]])], dataset_polarity="positive")
    assert t["adduct_tiers"] == [["[M+H]+"], ["[M+Na]+"]]
    assert set(t["theoretical_mz"]) == {"[M+H]+", "[M+Na]+"}


def test_rule_hash_is_the_same_whether_defaults_are_explicit():
    omitted = normalize_targets([_gaba()], dataset_polarity="positive")
    explicit = normalize_targets(
        [_gaba(ppm=10, rt_tol_min=0.3, adduct_tiers=DEFAULT_ADDUCT_TIERS["positive"],
               polarity="positive")],
        dataset_polarity="positive")
    assert rule_hash(omitted) == rule_hash(explicit)


def test_rule_hash_ignores_int_float_and_list_tuple_spelling():
    omitted = normalize_targets([_gaba()], dataset_polarity="positive")
    spelled = normalize_targets(
        [_gaba(ppm=10.0, rt_tol_min=0.3,
               adduct_tiers=tuple(tuple(t) for t in DEFAULT_ADDUCT_TIERS["positive"]))],
        dataset_polarity="positive")
    assert rule_hash(omitted) == rule_hash(spelled)
    assert isinstance(spelled[0]["ppm"], float) and isinstance(spelled[0]["adduct_tiers"][0], list)
    a = normalize_targets([_gaba(rt_min=3)], dataset_polarity="positive")
    b = normalize_targets([_gaba(rt_min=3.0)], dataset_polarity="positive")
    assert rule_hash(a) == rule_hash(b)
    mz = {"target_id": "c1", "inchikey": GABA_KEY, "adduct": "[M+H]+"}
    assert (rule_hash(normalize_targets([{**mz, "expected_mz": 104}], dataset_polarity="positive"))
            == rule_hash(normalize_targets([{**mz, "expected_mz": 104.0}], dataset_polarity="positive")))


def test_rule_hash_changes_with_the_rule():
    base = rule_hash(normalize_targets([_gaba()], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba(ppm=5)], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba(rt_min=3.0)], dataset_polarity="positive"))
    assert base != rule_hash(normalize_targets([_gaba()], dataset_polarity="negative"))


def test_target_of_other_polarity_is_not_applicable():
    [t] = normalize_targets([_gaba(polarity="negative")], dataset_polarity="positive")
    assert t["applicable"] is False
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["negative"]
    [t2] = normalize_targets([_gaba(polarity="positive")], dataset_polarity="positive")
    assert t2["applicable"] is True


@pytest.mark.parametrize("bad, field", [
    ({"target_id": "c1", "formula": "C4H9NO2"}, "inchikey"),
    ({"target_id": "c1", "inchikey": GABA_KEY}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07}, "adduct"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "adduct": "[M+H]+"}, "expected_mz"),
    ({"target_id": "c1", "inchikey": "short", "formula": "C4H9NO2"}, "inchikey"),
    ({"target_id": "", "inchikey": GABA_KEY, "formula": "C4H9NO2"}, "target_id"),
    ({"inchikey": GABA_KEY, "formula": "C4H9NO2"}, "target_id"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "bogus": 1}, "bogus"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "Xx9"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2!!"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "Ca(OH)2"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2-3"}, "formula"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "ppm": 0}, "ppm"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "ppm": True}, "ppm"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_tol_min": -1}, "rt_tol_min"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": "3"}, "rt_min"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "polarity": "pos"}, "polarity"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[nope]"}, "adduct"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "expected_mz": 104.07,
      "adduct": "[M+H]+"}, "expected_mz"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[M+H]+",
      "adduct_tiers": [["[M+H]+"]]}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
      "adduct_tiers": [["[M+Zz]+"]]}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "adduct_tiers": []}, "adduct_tiers"),
    ({"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2",
      "adduct_tiers": "[M+H]+"}, "adduct_tiers"),
])
def test_invalid_targets_name_the_field(bad, field):
    with pytest.raises(DomainError) as exc:
        normalize_targets([bad], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == field
    assert "target_id" in exc.value.details


def test_adduct_tiers_contradicting_polarity_is_invalid():
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(adduct_tiers=[["[M-H]-"]])], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == "adduct_tiers"
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(polarity="negative", adduct_tiers=[["[M+H]+"]])],
                          dataset_polarity="negative")
    assert exc.value.details["field"] == "adduct_tiers"


def test_duplicate_target_id_is_invalid():
    with pytest.raises(DomainError) as exc:
        normalize_targets([_gaba(), _gaba()], dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["target_id"] == "c1" and exc.value.details["field"] == "target_id"


def test_more_than_fifty_targets_is_invalid():
    many = [{"target_id": f"c{i}", "inchikey": GABA_KEY, "formula": "C4H9NO2"} for i in range(51)]
    with pytest.raises(DomainError) as exc:
        normalize_targets(many, dataset_polarity="positive")
    assert exc.value.code == "FEATURE_TARGET_INVALID"
    assert exc.value.details["field"] == "targets"
    assert len(normalize_targets(many[:50], dataset_polarity="positive")) == 50


def test_unknown_polarity_without_explicit_adducts_is_refused():
    explicit = {"target_id": "c2", "inchikey": GABA_KEY, "expected_mz": 104.07, "adduct": "[M+H]+"}
    tiers = _gaba(target_id="c3", adduct_tiers=[["[M+H]+"]])
    bare1 = _gaba(target_id="c4")
    bare2 = _gaba(target_id="c5")
    with pytest.raises(DomainError) as exc:
        normalize_targets([explicit, tiers, bare1, bare2], dataset_polarity=None)
    assert exc.value.code == "POLARITY_UNKNOWN"
    assert exc.value.details["target_ids"] == ["c4", "c5"]
    # 明示の付加体だけなら極性不明でも通る
    out = normalize_targets([explicit, tiers], dataset_polarity=None)
    assert [t["applicable"] for t in out] == [True, True]
    # 目標自身が極性を持てば、データセットの極性が不明でも既定の層が決まる
    [t] = normalize_targets([_gaba(polarity="positive")], dataset_polarity=None)
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["positive"] and t["applicable"] is True


# ---------- 目標照合の判定（spec §5） ----------
#
# ds は tests/test_feature_bindings.py の _dataset / _candidate に倣った最小の合成。
# SME 候補には Task 4 の「候補自身の」inchikey を入れる。

OTHER_KEY = "QWERTYUIOPASDF-UHFFFAOYSA-N"
THIRD_KEY = "ZXCVBNMLKJHGFD-UHFFFAOYSA-N"
MH = theoretical_mz("C4H9NO2", "[M+H]+")
MNA = theoretical_mz("C4H9NO2", "[M+Na]+")


def _shift(mz, ppm):
    return mz * (1 + ppm * 1e-6)


def _candidate(**overrides) -> dict:
    candidate = {
        "sme_id": "1", "rank": 1, "charge": 1,
        "exp_mass_to_charge": MH, "theoretical_mass_to_charge": MH,
        "confidence_value": 0.92, "confidence_measure": "MS-DIAL total score",
        "chemical_name": "GABA", "database_identifier": "HMDB:HMDB0000112",
        "inchi": None, "smiles": None, "adduct": "[M+H]1+",
        "identification_method": "matched to library", "spectra_ref": None,
        "inchikey": GABA_KEY, "inchikey_source": "msp",
    }
    candidate.update(overrides)
    return candidate


def _dataset(features: dict, *, assays=("assay[1]", "assay[2]")) -> DatasetState:
    """`{fid: {"mz":, "rt":, "inchikey":, "candidates": [...], "annotation": {...}}}` から最小の ds。"""
    ds = DatasetState()
    ds.feature_ids = list(features)
    ds.sample_assay_ids = list(assays)
    ds.sample_names = [f"S{i + 1}" for i in range(len(assays))]
    for fid, spec in features.items():
        ds.feature_metadata[fid] = {
            "mz": spec.get("mz"), "rt": spec.get("rt"),
            "name": spec.get("name"), "inchikey": spec.get("inchikey"),
        }
        ds.feature_candidates[fid] = list(spec.get("candidates", []))
        if spec.get("annotation") is not None:
            ds.feature_annotations[fid] = spec["annotation"]
    return ds


def _target(**extra) -> dict:
    [target] = normalize_targets([_gaba(**extra)], dataset_polarity="positive")
    return target


def _rows(out, feature_id):
    return [r for r in out["candidates"] if r["feature_id"] == feature_id]


def _row(out, feature_id, adduct="[M+H]+"):
    [row] = [r for r in _rows(out, feature_id) if r["adduct"] == adduct]
    return row


def _reasons(row):
    return row["reasons"].split(";") if row["reasons"] else []


EMPTY_COUNTS = {"mass": 0, "rt": 0, "qualified": 0,
                "by_evidence": {"sme_match": 0, "sml_match": 0, "unannotated": 0}}


def test_main_ion_tier_wins_over_other_adducts():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2}, "2": {"mz": MNA, "rt": 1.2}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["status"] == "resolved" and out["reason"] is None
    assert out["selected_feature_id"] == "1"
    assert out["selected_adduct"] == "[M+H]+"
    assert out["selected_evidence"] == "unannotated"
    assert out["adduct_tier_used"] == 1
    na = _row(out, "2", "[M+Na]+")
    assert na["qualified"] is True and na["adduct_tier"] == 2
    assert "lower_adduct_tier" in _reasons(na)
    assert "lower_adduct_tier" not in _reasons(_row(out, "1"))
    assert out["counts"] == {"mass": 2, "rt": 2, "qualified": 2,
                             "by_evidence": {"sme_match": 0, "sml_match": 0, "unannotated": 2}}


def test_sme_match_beats_sml_and_unannotated():
    ds = _dataset({
        "1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]},
        "2": {"mz": _shift(MH, 2), "rt": 1.2,
              "annotation": {"sml_id": "5", "ambiguous": False, "name": "GABA",
                             "inchikey": GABA_KEY}},
        "3": {"mz": _shift(MH, -3), "rt": 1.2},
    })
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["status"] == "resolved"
    assert out["selected_feature_id"] == "1"
    assert out["selected_evidence"] == "sme_match"
    row1 = _row(out, "1")
    assert row1["evidence"] == "sme_match" and row1["sme_rank"] == 1
    assert row1["annotated_name"] == "GABA" and row1["annotated_inchikey"] == GABA_KEY
    row2, row3 = _row(out, "2"), _row(out, "3")
    assert row2["evidence"] == "sml_match" and row3["evidence"] == "unannotated"
    assert row2["annotated_name"] == "GABA" and row2["annotated_inchikey"] == GABA_KEY
    for row in (row2, row3):
        assert row["qualified"] is True
        assert "lower_evidence_tier" in _reasons(row)
    assert out["counts"]["by_evidence"] == {"sme_match": 1, "sml_match": 1, "unannotated": 1}


def test_two_qualified_in_best_tier_is_needs_input_even_if_one_is_far_more_intense():
    ds = _dataset({
        "1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]},
        "2": {"mz": _shift(MH, 8), "rt": 1.45, "candidates": [_candidate(sme_id="2")]},
    })
    ds.feature_matrix = np.array([[1000.0, 1000.0], [1.0, 1.0]])
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["status"] == "needs_input"
    assert out["reason"] == "multiple_qualified_candidates"
    assert out["selected_feature_id"] is None and out["selected_adduct"] is None
    assert out["counts"]["qualified"] == 2
    # 強度を逆にしても結論は変わらない
    ds.feature_matrix = ds.feature_matrix[::-1].copy()
    again = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert again["status"] == "needs_input"
    assert again["reason"] == "multiple_qualified_candidates"


def test_no_candidate_is_needs_input_with_counts():
    ds = _dataset({"1": {"mz": 300.1234, "rt": 1.2, "candidates": [_candidate()]}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["status"] == "needs_input"
    assert out["reason"] == "no_qualified_candidate"
    assert out["selected_feature_id"] is None and out["adduct_tier_used"] is None
    assert out["counts"] == EMPTY_COUNTS
    assert out["candidates"] == [] and out["truncated"] == 0
    assert out["theoretical_mz"]["[M+H]+"] == pytest.approx(MH)


def test_rt_outside_tolerance_is_listed_but_not_qualified():
    ds = _dataset({"1": {"mz": MH, "rt": 2.0, "candidates": [_candidate()]},
                   "2": {"mz": MH, "rt": None}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    for fid in ("1", "2"):
        row = _row(out, fid)
        assert row["qualified"] is False
        assert "rt_outside_tolerance" in _reasons(row)
    assert _row(out, "1")["delta_rt"] == pytest.approx(0.8)
    assert _row(out, "2")["delta_rt"] is None
    assert out["counts"]["mass"] == 2 and out["counts"]["rt"] == 0
    assert out["reason"] == "no_qualified_candidate"


def test_feature_without_mz_is_skipped():
    ds = _dataset({"1": {"mz": None, "rt": 1.2}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["candidates"] == [] and out["counts"] == EMPTY_COUNTS


def test_sme_of_another_compound_is_dropped():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2,
                         "candidates": [_candidate(chemical_name="Other",
                                                   inchikey=OTHER_KEY)]}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "contradicted"
    assert "identifier_mismatch" in _reasons(row)
    assert row["qualified"] is False
    assert row["annotated_name"] == "Other" and row["annotated_inchikey"] == OTHER_KEY
    assert out["status"] == "needs_input" and out["reason"] == "no_qualified_candidate"
    assert out["counts"]["rt"] == 1 and out["counts"]["qualified"] == 0


def test_identity_adduct_and_charge_must_match():
    block = GABA_KEY[:14]
    # 付加体の書き方の違い（[M+H]1+ と [M+H]+）は正規化して同じとみなす
    assert classify_identity([_candidate(adduct="[M+H]1+")], None,
                             inchikey_block=block, adduct="[M+H]+") == ("sme_match", [], 1)
    assert classify_identity([_candidate(adduct="[M+Na]1+")], None,
                             inchikey_block=block, adduct="[M+H]+") \
        == ("contradicted", ["adduct_mismatch"], None)
    assert classify_identity([_candidate(charge=2)], None,
                             inchikey_block=block, adduct="[M+H]+") \
        == ("contradicted", ["charge_mismatch"], None)
    # 電荷は符号の書き方に依存しない（陰イオンで -1 と書かれていても 1 価として一致）
    assert classify_identity([_candidate(adduct="[M-H]1-", charge=-1)], None,
                             inchikey_block=block, adduct="[M-H]-")[0] == "sme_match"
    # 電荷が読めない候補は、付加体が一致していれば電荷で落とさない（欠測は反証ではない）。
    # 付加体が違えば電荷が読めなくても付加体の不一致で落ちる。
    assert classify_identity([_candidate(charge=None)], None,
                             inchikey_block=block, adduct="[M+H]+") == ("sme_match", [], 1)
    assert classify_identity([_candidate(charge=None, adduct="[M+Na]1+")], None,
                             inchikey_block=block, adduct="[M+H]+") \
        == ("contradicted", ["adduct_mismatch"], None)
    # 読めた電荷が違えば（None でなければ）従来どおり不一致
    assert classify_identity([_candidate(charge=0)], None,
                             inchikey_block=block, adduct="[M+H]+")[0] == "contradicted"


def test_sml_of_another_isomer_is_not_a_contradiction():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2,
                         "annotation": {"sml_id": "9", "ambiguous": False,
                                        "name": "Other isomer", "inchikey": OTHER_KEY}}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "unannotated"
    assert row["qualified"] is True
    assert row["annotated_name"] is None and row["annotated_inchikey"] is None
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"


def test_unannotated_needs_rt():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2}})
    out = match_target(ds, _target(), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "unannotated"
    assert row["qualified"] is False
    assert "mass_only_not_identified" in _reasons(row)
    assert row["delta_rt"] is None
    assert out["status"] == "needs_input" and out["reason"] == "no_qualified_candidate"
    assert out["counts"]["mass"] == 1 and out["counts"]["rt"] == 1
    # 同定があれば RT 無しでも qualified
    ds2 = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]}})
    assert match_target(ds2, _target(), max_results=20)["status"] == "resolved"


def test_rank2_sme_match_is_found_and_rank1_inchikey_not_borrowed():
    block = GABA_KEY[:14]
    rank1_other = _candidate(rank=1, chemical_name="Other", inchikey=OTHER_KEY)
    rank2_gaba = _candidate(sme_id="2", rank=2, chemical_name="GABA")
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "inchikey": OTHER_KEY,
                         "candidates": [rank1_other, rank2_gaba]}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "sme_match" and row["sme_rank"] == 2
    assert row["annotated_name"] == "GABA" and row["annotated_inchikey"] == GABA_KEY
    assert out["status"] == "resolved"

    # rank1 が GABA でも付加体が違い、rank2 は別物質 → rank1 の InChIKey を rank2 に貸さない
    rank1_gaba_na = _candidate(rank=1, adduct="[M+Na]1+")
    rank2_other = _candidate(sme_id="2", rank=2, chemical_name="Other", inchikey=OTHER_KEY)
    assert classify_identity([rank1_gaba_na, rank2_other], None,
                             inchikey_block=block, adduct="[M+H]+") \
        == ("contradicted", ["adduct_mismatch", "identifier_mismatch"], None)
    # 特徴の rank1 InChIKey（feature_metadata）が目標と一致していても、候補自身が別物質なら矛盾
    ds2 = _dataset({"1": {"mz": MH, "rt": 1.2, "inchikey": GABA_KEY,
                          "candidates": [_candidate(chemical_name="Other",
                                                    inchikey=OTHER_KEY)]}})
    out2 = match_target(ds2, _target(rt_min=1.2), max_results=20)
    assert _row(out2, "1")["evidence"] == "contradicted"


def test_unresolved_sme_identity_is_not_a_contradiction():
    unresolved = _candidate(chemical_name="Something", database_identifier="MspDB:Something",
                            inchikey=None, inchikey_source="none")
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [unresolved]}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "unannotated"
    assert "sme_identity_unresolved" in _reasons(row)
    assert row["qualified"] is True and row["sme_rank"] is None
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"
    # 識別子の無い候補しか無い特徴は SML へ進む
    sml = {"sml_id": "1", "ambiguous": False, "name": "GABA", "inchikey": GABA_KEY}
    assert classify_identity([unresolved], sml, inchikey_block=GABA_KEY[:14],
                             adduct="[M+H]+") \
        == ("sml_match", ["sme_identity_unresolved"], None)
    # 識別子の無い候補と、識別子を持ち一致する候補が並ぶなら一致が勝つ（理由は付けない）
    assert classify_identity([unresolved, _candidate(rank=2)], None,
                             inchikey_block=GABA_KEY[:14], adduct="[M+H]+") \
        == ("sme_match", [], 2)


def test_ambiguous_sml_with_a_matching_member_is_sml_match():
    ambiguous = {"ambiguous": True, "sml_ids": ["1", "2"], "name": None,
                 "inchikeys": [OTHER_KEY, GABA_KEY]}
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "annotation": ambiguous}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "sml_match"
    assert "sml_ambiguous" in _reasons(row)
    assert row["annotated_inchikey"] == f"{OTHER_KEY};{GABA_KEY}"
    assert out["status"] == "resolved" and out["selected_evidence"] == "sml_match"
    # 曖昧な注釈のどれも一致しなければ unannotated（矛盾ではない）
    none_match = dict(ambiguous, inchikeys=[OTHER_KEY, THIRD_KEY])
    assert classify_identity([], none_match, inchikey_block=GABA_KEY[:14],
                             adduct="[M+H]+") == ("unannotated", [], None)


def test_same_feature_under_two_adducts_counts_once_per_tier():
    # 広い ppm で、1 つの特徴が第 1 層の [M+H]+ と [M]+ の両方の窓に入る。
    # 特徴は [M]+ の理論値のほうがずっと近いが、選ぶ付加体は層に並べた順（[M+H]+）。
    m_plus = theoretical_mz("C4H9NO2", "[M]+")
    target = _target(rt_min=1.2, ppm=9000)
    ds = _dataset({"1": {"mz": 103.2, "rt": 1.2}})
    out = match_target(ds, target, max_results=20)
    assert abs(103.2 - m_plus) < abs(103.2 - MH)
    assert sorted(r["adduct"] for r in _rows(out, "1")) == ["[M+H]+", "[M]+"]
    assert all(r["qualified"] for r in _rows(out, "1"))
    assert out["status"] == "resolved"
    assert out["selected_feature_id"] == "1"
    assert out["selected_adduct"] == "[M+H]+"
    assert out["counts"]["mass"] == 1 and out["counts"]["qualified"] == 1
    assert out["counts"]["by_evidence"]["unannotated"] == 1
    # 2 つ目の特徴が同じ層に入れば、行は 4 本でも特徴は 2 つ → needs_input
    ds2 = _dataset({"1": {"mz": 103.2, "rt": 1.2}, "2": {"mz": 103.5, "rt": 1.2}})
    out2 = match_target(ds2, target, max_results=20)
    assert len(out2["candidates"]) == 4
    assert out2["status"] == "needs_input"
    assert out2["reason"] == "multiple_qualified_candidates"
    assert out2["adduct_tier_used"] == 1


def test_candidates_are_truncated_with_count():
    ds = _dataset({str(i): {"mz": _shift(MH, ppm), "rt": 1.2}
                   for i, ppm in enumerate([-9, 1, 5, -3])})
    out = match_target(ds, _target(), max_results=2)
    assert len(out["candidates"]) == 2
    assert out["truncated"] == 2
    assert out["counts"]["mass"] == 4
    # 打ち切りは並べたあと（|ppm| の小さい 2 件が残る）
    assert [r["feature_id"] for r in out["candidates"]] == ["1", "3"]


def test_candidate_rows_have_keys_and_order():
    ds = _dataset({
        "far": {"mz": _shift(MH, -9), "rt": 1.2},
        "near": {"mz": _shift(MH, 1), "rt": 1.2, "candidates": [_candidate()]},
        "na": {"mz": _shift(MNA, 0.5), "rt": 1.2},
        "late": {"mz": _shift(MH, 0.1), "rt": 3.0},
    })
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert list(out["candidates"][0]) == [
        "target_id", "feature_id", "adduct", "adduct_tier", "theoretical_mz", "mz", "ppm",
        "rt", "delta_rt", "evidence", "sme_rank", "annotated_name", "annotated_inchikey",
        "qualified", "reasons"]
    # qualified が先、層の昇順、|ppm| の昇順
    assert [r["feature_id"] for r in out["candidates"]] == ["near", "far", "na", "late"]
    near = _row(out, "near")
    assert near["target_id"] == "c1" and near["reasons"] == ""
    assert near["ppm"] == pytest.approx(1.0, abs=1e-6)
    assert _row(out, "far")["ppm"] == pytest.approx(-9.0, abs=1e-6)


def test_target_of_other_polarity_is_not_applicable_in_matching():
    [target] = normalize_targets([_gaba(polarity="negative")], dataset_polarity="positive")
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]}})
    out = match_target(ds, target, max_results=20)
    assert out["status"] == "not_applicable" and out["reason"] == "polarity_mismatch"
    assert out["candidates"] == [] and out["truncated"] == 0
    assert out["selected_feature_id"] is None
    assert out["counts"] == EMPTY_COUNTS


# ---------- Task 7 の持ち越し（match_target の補足） ----------


def test_needs_input_tie_records_the_evidence_level_where_it_happened():
    # 同じ層・同じ水準（sme_match）に 2 特徴 → 選ばないが、決着が付かなかった水準は残す
    ds = _dataset({
        "1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]},
        "2": {"mz": _shift(MH, 3), "rt": 1.2, "candidates": [_candidate(sme_id="2")]},
        "3": {"mz": _shift(MH, -3), "rt": 1.2},          # 下位の unannotated は水準に影響しない
    })
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    assert out["status"] == "needs_input" and out["reason"] == "multiple_qualified_candidates"
    assert out["selected_feature_id"] is None and out["selected_adduct"] is None
    assert out["selected_evidence"] == "sme_match"
    assert out["adduct_tier_used"] == 1
    # 候補が 1 つも無いときは決着した水準が無いので None のまま
    empty = match_target(_dataset({"1": {"mz": 300.0, "rt": 1.2}}), _target(rt_min=1.2),
                         max_results=20)
    assert empty["selected_evidence"] is None
    # unannotated 同士の同点は unannotated
    tied = match_target(
        _dataset({"1": {"mz": MH, "rt": 1.2}, "2": {"mz": _shift(MH, 2), "rt": 1.2}}),
        _target(rt_min=1.2), max_results=20)
    assert tied["status"] == "needs_input" and tied["selected_evidence"] == "unannotated"


def test_precomputed_mz_array_is_used_when_given():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]}})
    target = _target(rt_min=1.2)
    same = match_target(ds, target, max_results=20, mz_array=_feature_mz_array(ds))
    assert same == match_target(ds, target, max_results=20)
    # 渡した配列がそのまま使われる（ds の m/z ではなく）ことを、窓の外へずらして確かめる
    shifted = match_target(ds, target, max_results=20, mz_array=np.array([MH + 50.0]))
    assert shifted["candidates"] == [] and shifted["status"] == "needs_input"


def test_expected_mz_target_matches_only_its_fixed_adduct():
    # 理論値を expected_mz で与える目標は付加体が 1 つに固定される（MNA の特徴は拾わない）
    [target] = normalize_targets(
        [{"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": MH, "adduct": "[M+H]+",
          "rt_min": 1.2}],
        dataset_polarity="positive")
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [_candidate()]},
                   "2": {"mz": MNA, "rt": 1.2}})
    out = match_target(ds, target, max_results=20)
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"
    assert out["selected_adduct"] == "[M+H]+" and out["selected_evidence"] == "sme_match"
    assert out["theoretical_mz"] == {"[M+H]+": MH}
    assert [r["feature_id"] for r in out["candidates"]] == ["1"]
    assert out["candidates"][0]["sme_rank"] == 1


def test_negative_polarity_target_uses_the_negative_default_tiers():
    mh_neg = theoretical_mz("C4H9NO2", "[M-H]-")
    [target] = normalize_targets([_gaba(rt_min=1.2)], dataset_polarity="negative")
    cand = _candidate(adduct="[M-H]1-", charge=-1, exp_mass_to_charge=mh_neg)
    ds = _dataset({"1": {"mz": mh_neg, "rt": 1.2, "candidates": [cand]},
                   "2": {"mz": theoretical_mz("C4H9NO2", "[M+Cl]-"), "rt": 1.2}})
    out = match_target(ds, target, max_results=20)
    assert out["status"] == "resolved" and out["selected_feature_id"] == "1"
    assert out["selected_adduct"] == "[M-H]-" and out["adduct_tier_used"] == 1
    assert out["selected_evidence"] == "sme_match"
    cl = _row(out, "2", "[M+Cl]-")
    assert cl["adduct_tier"] == 2 and "lower_adduct_tier" in _reasons(cl)
    # 陽イオンの付加体は候補に出ない
    negative = {a for tier in DEFAULT_ADDUCT_TIERS["negative"] for a in tier}
    assert all(r["adduct"] in negative for r in out["candidates"])


def test_sme_rank_is_none_when_the_matched_candidate_has_no_rank():
    ds = _dataset({"1": {"mz": MH, "rt": 1.2, "candidates": [_candidate(rank=None)]}})
    out = match_target(ds, _target(rt_min=1.2), max_results=20)
    row = _row(out, "1")
    assert row["evidence"] == "sme_match" and row["sme_rank"] is None
    assert out["status"] == "resolved"


# ---------- 窓検索と統計の結合（spec §4.2・§6） ----------

def _stats_dataset():
    """4 特徴 × 4 試料。A=S1,S2 / B=S3,S4。特徴 4 は前処理で落ちた（結果の行列に無い）。"""
    ds = _dataset({
        "1": {"mz": 100.0, "rt": 1.0, "candidates": [_candidate(chemical_name="Alpha")]},
        "2": {"mz": 200.0, "rt": 2.0,
              "annotation": {"sml_id": "3", "ambiguous": False, "name": "Beta",
                             "inchikey": OTHER_KEY, "adduct": "[M+Na]+"}},
        "3": {"mz": 300.0, "rt": 3.0},
        "4": {"mz": 400.0, "rt": 4.0},
    }, assays=("a1", "a2", "a3", "a4"))
    ds.dataset_id = "ds1"
    ds.detected_mask = np.array([
        [True, True, True, False],
        [True, False, False, False],
        [False, False, True, True],
        [True, True, True, True],
    ])
    return ds


def _two_group(ds, rows, *, result_id="diff_1", **extra):
    result = {
        "kind": "two_group", "a": "ctrl", "b": "treat",
        "samples_a": ["S1", "S2"], "samples_b": ["S3", "S4"], "n_a": 2, "n_b": 2,
        "log2fc_sign": "log2(b/a)", "caveats": ["小n"],
        "results": rows,
        "provenance": {"result_id": result_id, "kind": "differential", "dataset_id": "ds1",
                       "parent_ids": []},
    }
    result.update(extra)
    ds.results[result_id] = result
    return result


def _row_of(feature, **extra):
    row = {"feature": feature, "mean_a": 1.0, "mean_b": 2.0, "log2fc": 1.0,
           "t": 2.5, "p": 0.01, "q": 0.03}
    row.update(extra)
    return row


def _window(ds, **kw):
    args = dict(feature_ids=None, min_mz=None, max_mz=None, min_rt=None, max_rt=None,
                max_results=50)
    args.update(kw)
    return search_window(ds, **args)


def test_window_and_filters_are_anded_and_inclusive():
    ds = _stats_dataset()
    rows, truncated = _window(ds, min_mz=200.0, max_mz=300.0)
    assert [r["feature_id"] for r in rows] == ["2", "3"] and truncated == 0   # 境界を含む
    rows, _ = _window(ds, min_mz=100.0, max_mz=400.0, min_rt=2.0, max_rt=3.0)
    assert [r["feature_id"] for r in rows] == ["2", "3"]
    rows, _ = _window(ds, min_mz=100.0, max_mz=400.0, min_rt=2.0, max_rt=3.0,
                      feature_ids=["3", "4"])
    assert [r["feature_id"] for r in rows] == ["3"]                          # 条件は AND
    rows, _ = _window(ds, feature_ids=["4", "1"])
    assert [r["feature_id"] for r in rows] == ["1", "4"]                     # m/z 昇順
    rows, _ = _window(ds, min_rt=3.5)                                        # 片側だけでも条件
    assert [r["feature_id"] for r in rows] == ["4"]
    assert list(rows[0]) == ["feature_id", "mz", "rt", "evidence", "name", "inchikey",
                             "adduct", "n_sme"]


def test_window_evidence_name_inchikey_adduct_and_truncation():
    ds = _stats_dataset()
    ds.feature_candidates["1"] = [
        _candidate(rank=2, chemical_name="Second", inchikey=OTHER_KEY, adduct="[M+Na]1+"),
        _candidate(rank=1, chemical_name="Alpha", adduct="[M+H]1+"),
    ]
    ds.feature_annotations["3"] = {"ambiguous": True, "sml_ids": ["1", "2"], "name": None,
                                   "inchikeys": [OTHER_KEY, THIRD_KEY]}
    rows, truncated = _window(ds, min_mz=0)
    by = {r["feature_id"]: r for r in rows}
    assert by["1"] == {"feature_id": "1", "mz": 100.0, "rt": 1.0, "evidence": "sme",
                       "name": "Alpha", "inchikey": GABA_KEY, "adduct": "[M+H]1+", "n_sme": 2}
    assert by["2"]["evidence"] == "sml" and by["2"]["name"] == "Beta"
    assert by["2"]["inchikey"] == OTHER_KEY and by["2"]["adduct"] == "[M+Na]+"
    assert by["2"]["n_sme"] == 0
    assert by["3"]["evidence"] == "sml" and by["3"]["name"] is None
    assert by["3"]["inchikey"] == f"{OTHER_KEY};{THIRD_KEY}"
    assert by["4"] == {"feature_id": "4", "mz": 400.0, "rt": 4.0, "evidence": "none",
                       "name": None, "inchikey": None, "adduct": None, "n_sme": 0}
    rows, truncated = _window(ds, min_mz=0, max_results=3)
    assert [r["feature_id"] for r in rows] == ["1", "2", "3"] and truncated == 1


def test_empty_query_is_refused():
    ds = _stats_dataset()
    with pytest.raises(DomainError) as exc:
        _window(ds)
    assert exc.value.code == "FEATURE_SEARCH_EMPTY_QUERY"
    for kw in ({"min_mz": 500.0, "max_mz": 400.0}, {"min_rt": 3.0, "max_rt": 2.0}):
        with pytest.raises(DomainError) as exc:
            _window(ds, **kw)
        assert exc.value.code == "FEATURE_SEARCH_BAD_RANGE"


def test_unannotated_feature_gets_its_statistics():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1"), _row_of("3", p=0.2, q=0.4, log2fc=-0.5)])
    [row] = feature_stats(ds, ["3"], [result])
    assert row["feature_id"] == "3" and row["result_id"] == "diff_1"
    assert row["status"] == "tested"
    assert (row["log2fc"], row["p"], row["q"]) == (-0.5, 0.2, 0.4)
    assert (row["mean_a"], row["mean_b"]) == (1.0, 2.0)
    assert list(row) == ["feature_id", "result_id", "status", "log2fc", "p", "q", "mean_a",
                         "mean_b", "detected_a", "detected_b"]
    # 行は特徴が外側、結果が内側の順
    other = _two_group(ds, [_row_of("1")], result_id="diff_2")
    rows = feature_stats(ds, ["1", "3"], [result, other])
    assert [(r["feature_id"], r["result_id"], r["status"]) for r in rows] == [
        ("1", "diff_1", "tested"), ("1", "diff_2", "tested"),
        ("3", "diff_1", "tested"), ("3", "diff_2", "not_in_matrix")]


def test_feature_dropped_by_preprocess_is_not_in_matrix():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1")])
    [row] = feature_stats(ds, ["4"], [result])
    assert row["status"] == "not_in_matrix"
    for key in ("log2fc", "p", "q", "mean_a", "mean_b"):
        assert row[key] is None
    # 検出例数はマスクから取れるので、なぜ落ちたかの手がかりとして残す
    assert (row["detected_a"], row["detected_b"]) == (2, 2)


def test_null_p_from_persisted_result_is_not_testable():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1", p=None, q=None),
                             _row_of("2", p=float("nan"), q=None)])
    rows = feature_stats(ds, ["1", "2"], [result])
    assert [r["status"] for r in rows] == ["not_testable", "not_testable"]
    for r in rows:
        for key in ("log2fc", "p", "q", "mean_a", "mean_b"):
            assert r[key] is None


def test_q_is_taken_from_the_result_not_recomputed():
    ds = _stats_dataset()
    # BH を掛け直せば別の値になる p でも、結果の q をそのまま返す
    result = _two_group(ds, [_row_of("1", p=0.04, q=0.987654), _row_of("2", p=0.01, q=0.5)])
    rows = feature_stats(ds, ["1", "2"], [result])
    assert [r["q"] for r in rows] == [0.987654, 0.5]
    # q が無い（None）が p は有限なら tested のまま q は None
    result = _two_group(ds, [_row_of("1", p=0.04, q=None)], result_id="diff_3")
    [row] = feature_stats(ds, ["1"], [result])
    assert row["status"] == "tested" and row["p"] == 0.04 and row["q"] is None


def test_bh_population_counts():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1"), _row_of("2", p=None, q=None),
                             _row_of("3", p=float("nan"), q=float("nan"))])
    [summary] = result_summaries(ds, [result])
    assert summary == {
        "result_id": "diff_1", "group_a": "ctrl", "group_b": "treat", "n_a": 2, "n_b": 2,
        "log2fc_sign": "log2(b/a)",
        "bh_population": {"n_features_total": 4, "n_in_matrix": 3, "n_tested": 1},
        "caveats": ["小n"],
    }
    summary["caveats"].append("x")
    assert result["caveats"] == ["小n"]                         # 元の結果を汚さない


def test_default_selection_skips_stale_and_non_two_group():
    ds = _stats_dataset()
    good = _two_group(ds, [_row_of("1")], result_id="diff_ok")
    _two_group(ds, [_row_of("1")], result_id="diff_old")
    ds.results["diff_old"]["provenance"]["dataset_id"] = "other"      # 古い
    ds.results["stat_1"] = {"kind": "welch", "results": [], "provenance": {
        "result_id": "stat_1", "kind": "statistic", "dataset_id": "ds1", "parent_ids": []}}
    selected = select_results(ds, None)
    assert selected["results"] == [good]
    assert selected["stale"] == ["diff_old"] and selected["missing"] == []
    assert select_results(_stats_dataset(), None) == {"results": [], "stale": [], "missing": []}


def test_named_stale_result_is_reported():
    ds = _stats_dataset()
    good = _two_group(ds, [_row_of("1")], result_id="diff_ok")
    _two_group(ds, [_row_of("1")], result_id="diff_old")
    ds.results["diff_old"]["provenance"]["dataset_id"] = "other"
    ds.results["stat_1"] = {"kind": "welch", "results": []}
    selected = select_results(ds, ["diff_ok", "diff_old", "nope", "stat_1"])
    assert selected["results"] == [good]
    assert selected["stale"] == ["diff_old"]
    # 存在しない ID と、2 群比較でない結果（v2 の統計）は missing（2 群比較としては引けない）
    assert selected["missing"] == ["nope", "stat_1"]


def test_detection_counts_per_group():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1"), _row_of("2"), _row_of("3")])
    rows = {r["feature_id"]: r for r in feature_stats(ds, ["1", "2", "3"], [result])}
    # マスクの行は ds.feature_ids の位置、列は samples_a/b を sample_names の位置に直す
    assert (rows["1"]["detected_a"], rows["1"]["detected_b"]) == (2, 1)
    assert (rows["2"]["detected_a"], rows["2"]["detected_b"]) == (1, 0)
    assert (rows["3"]["detected_a"], rows["3"]["detected_b"]) == (0, 2)
    # 結果の試料の順が sample_names と違っても名前で引く
    swapped = _two_group(ds, [_row_of("1")], result_id="diff_sw",
                         samples_a=["S4", "S3"], samples_b=["S2", "S1"])
    [row] = feature_stats(ds, ["1"], [swapped])
    assert (row["detected_a"], row["detected_b"]) == (1, 2)


def test_detection_counts_skipped_when_mask_shape_mismatches():
    ds = _stats_dataset()
    result = _two_group(ds, [_row_of("1")])
    ds.detected_mask = ds.detected_mask[:, :3]                 # 列数が試料数と合わない
    [row] = feature_stats(ds, ["1"], [result])
    assert row["status"] == "tested"
    assert row["detected_a"] is None and row["detected_b"] is None
    ds.detected_mask = None
    [row] = feature_stats(ds, ["1"], [result])
    assert row["detected_a"] is None and row["detected_b"] is None
