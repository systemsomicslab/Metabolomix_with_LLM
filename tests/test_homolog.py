"""同族列の RT の単調性（spec 2026-10-22 §4 の 6）。"""
from __future__ import annotations

import pytest

from metabolomix.analysis import homolog as H
from metabolomix.core.atomic_io import DomainError

RULE = H.normalize_spec({"members": [{"feature_id": "a", "carbons": 1, "double_bonds": 0},
                                     {"feature_id": "b", "carbons": 2, "double_bonds": 0}]})["rule"]


def _members(*rows):
    return [{"feature_id": f, "carbons": c, "double_bonds": d} for f, c, d in rows]


def test_defaults_and_validation():
    spec = H.normalize_spec({"auto": {"lipid_class": "PC"}})
    assert spec["rule"] == {"carbon_direction": "increasing", "double_bond_direction": "decreasing",
                            "max_violations": 0}
    assert spec["auto"] == {"lipid_class": "PC", "evidence": ["sme"], "adduct": None}
    with pytest.raises(DomainError):
        H.normalize_spec({"members": _members(("a", 34, 1), ("b", 36, 1)),
                          "auto": {"lipid_class": "PC"}})
    with pytest.raises(DomainError):
        H.normalize_spec({"auto": {"lipid_class": "PC"}, "rule": {"wrong": 1}})


def test_monotone_series_passes():
    members = _members(("a", 34, 1), ("b", 36, 1), ("c", 38, 1), ("d", 36, 2))
    out = H.check(members, {"a": 7.0, "b": 8.0, "c": 9.0, "d": 7.5}, RULE)
    assert out["verdict"] == "pass" and out["n_violations"] == 0


def test_carbon_and_double_bond_violations_counted():
    members = _members(("a", 34, 1), ("b", 36, 1), ("c", 36, 2))
    # 36:1 が 34:1 より早い（炭素の向きに反する）、36:2 が 36:1 より遅い（二重結合の向きに反する）
    out = H.check(members, {"a": 8.0, "b": 7.0, "c": 7.5}, RULE)
    assert out["verdict"] == "fail" and out["n_violations"] == 2
    axes = sorted(v["axis"] for v in out["violations"])
    assert axes == ["carbons", "double_bonds"]
    status = {m["feature_id"]: m["status"] for m in out["members"]}
    assert status == {"a": "violation", "b": "violation", "c": "violation"}


def test_max_violations_allows_outliers():
    rule = {**RULE, "max_violations": 1}
    out = H.check(_members(("a", 34, 1), ("b", 36, 1)), {"a": 8.0, "b": 7.0}, rule)
    assert out["verdict"] == "pass" and out["n_violations"] == 1


def test_duplicates_use_median_and_are_listed():
    members = _members(("a", 34, 1), ("a2", 34, 1), ("b", 36, 1))
    out = H.check(members, {"a": 7.0, "a2": 7.4, "b": 7.3}, RULE)
    assert out["duplicates"] == {"34:1": ["a", "a2"]}
    assert out["n_violations"] == 0                      # 中央値 7.2 < 7.3


def test_species_from_names():
    assert H.species_of("PC 34:1") == ("PC", 34, 1)
    assert H.species_of("PC 16:0_18:1") == ("PC", 34, 1)
    assert H.species_of("PC O-34:1") == ("PC O-", 34, 1)
    assert H.species_of("Glucose") is None


def test_members_from_names_filters_class_and_adduct():
    features = [
        {"feature_id": "1", "names": {"sme": "PC 34:1"}, "adduct": "[M+H]+"},
        {"feature_id": "2", "names": {"sme": "PE 34:1"}, "adduct": "[M+H]+"},
        {"feature_id": "3", "names": {"sme": "PC 36:2"}, "adduct": "[M+Na]+"},
        {"feature_id": "4", "names": {"sml": "PC 38:4"}, "adduct": "[M+H]+"},
    ]
    auto = {"lipid_class": "PC", "evidence": ["sme"], "adduct": "[M+H]+"}
    assert [m["feature_id"] for m in H.members_from_names(features, auto)] == ["1"]
    auto = {"lipid_class": "PC", "evidence": ["sme", "sml"], "adduct": None}
    assert [m["feature_id"] for m in H.members_from_names(features, auto)] == ["1", "3", "4"]


# ---------- 付加体は組成で比べ、空の列は pass にしない（最終検証 F26） ----------

@pytest.mark.parametrize("feature_adduct, wanted", [
    ("[M+H]1+", "[M+H]+"), ("[M+H]+", "[M+H]1+"), ("[M+HCOO]1-", "[M+FA-H]-"), ("[M+CHO2]-", "[M+HCOO]-"),
])
def test_members_from_names_compares_adducts_by_composition(feature_adduct, wanted):
    features = [{"feature_id": "1", "names": {"sme": "PC 34:1"}, "adduct": feature_adduct},
                {"feature_id": "2", "names": {"sme": "PC 36:1"}, "adduct": "[M+Na]1+"}]
    auto = {"lipid_class": "PC", "evidence": ["sme"], "adduct": wanted}
    assert [m["feature_id"] for m in H.members_from_names(features, auto)] == ["1"]


def test_nothing_to_compare_is_not_evaluable_not_pass():
    assert H.check([], {}, RULE)["verdict"] == "not_evaluable"
    one = H.check(_members(("a", 34, 1)), {"a": 7.0}, RULE)
    assert one["verdict"] == "not_evaluable" and one["n_pairs"] == 0
    # 隣どうしにならない 2 種（炭素数も二重結合数も違う）・RT の無い種も比べられない
    apart = H.check(_members(("a", 34, 1), ("b", 36, 2)), {"a": 8.0, "b": 7.0}, RULE)
    assert apart["verdict"] == "not_evaluable" and apart["n_pairs"] == 0
    no_rt = H.check(_members(("a", 34, 1), ("b", 36, 1)), {"a": 8.0}, RULE)
    assert no_rt["verdict"] == "not_evaluable"
    # 比べた対があれば従来どおり（34:1–36:1・36:1–38:1・36:1–36:2 の 3 対）
    out = H.check(_members(("a", 34, 1), ("b", 36, 1), ("c", 38, 1), ("d", 36, 2)),
                  {"a": 7.0, "b": 8.0, "c": 9.0, "d": 7.5}, RULE)
    assert out["verdict"] == "pass" and out["n_pairs"] == 3
