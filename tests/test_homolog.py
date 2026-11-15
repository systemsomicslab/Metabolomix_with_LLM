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
    assert H.species_of("PC 34:1") == ("PC", 34, 1, "")
    assert H.species_of("PC 16:0_18:1") == ("PC", 34, 1, "")
    assert H.species_of("PC O-34:1") == ("PC O-", 34, 1, "")
    assert H.species_of("Glucose") is None


@pytest.mark.parametrize("name, species", [
    ("PC 34:1; ExampleDB ID-0001", ("PC", 34, 1, "")),    # MSP の記録名（注記付き）がそのまま SME の名前になる
    ("PE 36:2;O (oxPE)", ("PE", 36, 2, "O")),             # 酸化 PE の修飾は落とさない（PE 36:2 の重複ではない）
    ("Phosphatidylcholine 16:0-18:1", ("PC", 34, 1, "")),
    ("SM d18:1/16:0", ("SM", 34, 1, "O2")),               # スフィンゴ脂質の種のレベルは ;O2 を持つ
    ("Cer 18:0;O3/16:0", ("Cer", 34, 0, "O3")),           # Cer NP（;O3）は Cer NDS（;O2）と別の列
    ("PE P-38:4", ("PE O-", 38, 5, "")),                  # GOSLIN の種のレベル: P- は O- に直る
    ("PE P-18:0_20:4", ("PE O-", 38, 5, "")),             # 分子種でも種のレベルと同じ（以前は PE P- 38:4 で別の列）
    ("PE O-36:2;O", ("PE O-", 36, 2, "O")),
    ("PC 15:0_18:1(d7)", None),                           # 同位体の標識のある内部標準は読まない
])
def test_species_are_read_like_the_scan_and_the_msms_side(name, species):
    """同族列の名前は、目標照合の種の目標・MS/MS の合否と同じ読み手（`msms_rules.name_lipid_species`）で読む
    （レビュー V2）。種の文字列を独立に（GOSLIN の種のレベルの鎖の和として）確かめる。"""
    from metabolomix.analysis.msms_rules import name_lipid_species

    assert H.species_of(name) == species
    read, _route = name_lipid_species(name)
    if species is None:
        assert read is None
    else:
        head = species[0].replace(" O-", "")
        ether = "O-" if species[0].endswith(" O-") else ""
        core, _, modification = read.partition(";")
        assert core == f"{head} {ether}{species[1]}:{species[2]}"
        assert modification == species[3]


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


# ---------- 修飾（;O…）の違う種は別の列（最終検証 WB） ----------

def _auto(lipid_class, adduct=None):
    return {"lipid_class": lipid_class, "evidence": ["sme"], "adduct": adduct}


def _features(names):
    return [{"feature_id": fid, "names": {"sme": name}, "adduct": "[M+H]+"} for fid, name in names]


def _plain_series_violations(names_rt):
    """独立の計算: 修飾（`;`）もエーテル（`O-` / `P-`）も無い `PE a:b[_c:d]` の名前だけを鎖の和の種に直し、
    同じ種の RT の中央値で、二重結合の隣どうし（同じ炭素数）が逆相の向き（二重結合が増えると早まる）に
    反する対の数を数える。"""
    from statistics import median as med

    by_species = {}
    for name, rt in names_rt:
        tail = name.split(" ", 1)[1]
        if ";" in tail or tail[:2] in ("O-", "P-"):
            continue
        chains = [tuple(map(int, c.split(":"))) for c in tail.split(" ")[0].split("_")]
        key = (sum(c for c, _ in chains), sum(d for _, d in chains))
        by_species.setdefault(key, []).append(rt)
    rts = {k: med(v) for k, v in by_species.items()}
    return sum(1 for (c, d), rt in rts.items() if (c, d + 1) in rts and rts[(c, d + 1)] > rt)


def test_oxidized_species_do_not_hide_a_violation_of_the_parent_series():
    """`PE 36:2;O (oxPE)` は `PE 36:2` と同じ（炭素数, 二重結合数）ではない。混ぜると RT の中央値が下がり、
    36:1（11.0）→ 36:2 の本当の違反が消えていた（最終検証 WB）。エーテル（`PE O-36:2`）と
    プラスマローゲン（`PE P-36:2` は種のレベルで `PE O-36:3`）は無修飾の PE の列に入らない。"""
    names_rt = [("PE 18:0_18:2", 12.0), ("PE 36:2;O (oxPE)", 9.0), ("PE 18:0_18:1", 11.0),
                ("PE O-36:2", 13.0), ("PE P-36:2", 12.5), ("PE 36:2", 12.2)]
    ids = [f"F{i}" for i in range(len(names_rt))]
    features = _features(zip(ids, (n for n, _ in names_rt)))
    rt_of = dict(zip(ids, (rt for _, rt in names_rt)))
    members = H.members_from_names(features, _auto("PE", "[M+H]+"))
    assert {m["feature_id"]: m["modification"] for m in members} == {"F0": "", "F1": "O", "F2": "", "F5": ""}
    out = H.check(members, rt_of, RULE)
    # 無修飾の列: 36:1 = 11.0、36:2 = 中央値(12.0, 12.2) = 12.1 → 二重結合が増えて遅れる 1 対が違反
    assert _plain_series_violations(names_rt) == 1
    assert out["verdict"] == "fail" and out["n_violations"] == 1 and out["n_pairs"] == 1
    [violation] = out["violations"]
    assert (violation["low"], violation["high"]) == ("36:1", "36:2")
    assert violation["rt_high"] == pytest.approx(12.1)
    assert violation["features_high"] == ["F0", "F5"]          # 酸化 PE（F1）は重複に入らない
    assert out["duplicates"] == {"36:2": ["F0", "F5"]}
    status = {m["feature_id"]: (m["modification"], m["status"]) for m in out["members"]}
    assert status["F1"] == ("O", "ok")
    # エーテルの列: O-36:2（13.0）と P-36:2 = O-36:3（12.5）は無修飾の PE と混ざらず、互いに 36:2 → 36:3
    ether = H.members_from_names(features, _auto("PE O-", "[M+H]+"))
    assert {(m["feature_id"], m["carbons"], m["double_bonds"]) for m in ether} == {("F3", 36, 2), ("F4", 36, 3)}
    assert H.check(ether, rt_of, RULE)["verdict"] == "pass"


def test_an_early_oxidized_species_does_not_fail_a_clean_series():
    """PE 34:1 → 36:1 → 38:1 は炭素数とともに遅れる（違反なし）。早く溶出する PE 36:1;O を 36:1 に数えると
    中央値が 8.0 に下がり、きれいな列が fail になっていた（最終検証 WB）。"""
    features = _features([("a", "PE 34:1"), ("b", "PE 36:1"), ("c", "PE 38:1"), ("d", "PE 36:1;O")])
    rt_of = {"a": 10.0, "b": 11.0, "c": 12.0, "d": 5.0}
    out = H.check(H.members_from_names(features, _auto("PE")), rt_of, RULE)
    assert out["verdict"] == "pass" and out["n_violations"] == 0 and out["duplicates"] == {}
    assert out["n_pairs"] == 2                       # 34:1–36:1・36:1–38:1。36:1;O は 1 種だけの列


def test_ceramides_with_different_hydroxylation_are_separate_series():
    """Cer NP（`Cer 18:0;O3/16:0` = Cer 34:0;O3）は Cer NDS（Cer 34:0;O2）の重複ではない。混ぜると 34:0 の
    中央値が 9.7 に下がり、NS 34:1（10.0）との対が違反になっていた（最終検証 WB）。"""
    features = _features([("NS34", "Cer 18:1;O2/16:0"), ("NS36", "Cer 18:1;O2/18:0"),
                          ("NDS34", "Cer 18:0;O2/16:0"), ("NP34", "Cer 18:0;O3/16:0")])
    rt_of = {"NS34": 10.0, "NS36": 11.0, "NDS34": 10.4, "NP34": 9.0}
    out = H.check(H.members_from_names(features, _auto("Cer")), rt_of, RULE)
    assert out["verdict"] == "pass" and out["duplicates"] == {}
    assert out["n_pairs"] == 2                       # ;O2 の 34:1–36:1 と 34:0–34:1
    # 違反の種は修飾を添えて書く（;O2 の列の中の 34:0 → 34:1）
    late = H.check(H.members_from_names(features, _auto("Cer")), {**rt_of, "NDS34": 9.5}, RULE)
    [violation] = late["violations"]
    assert (violation["low"], violation["high"]) == ("34:0;O2", "34:1;O2")
    # lipid_class に修飾を書けばその修飾の種だけを集める
    assert [m["feature_id"] for m in H.members_from_names(features, _auto("Cer;O3"))] == ["NP34"]
    assert [m["feature_id"] for m in H.members_from_names(features, _auto("Cer;O2"))] == [
        "NS34", "NS36", "NDS34"]


def test_hand_written_members_keep_the_old_series():
    # members（手で並べる）には修飾が無い: 鍵も表記も前と同じ
    out = H.check(_members(("a", 34, 1), ("a2", 34, 1), ("b", 36, 1)), {"a": 7.0, "a2": 7.4, "b": 7.3}, RULE)
    assert out["duplicates"] == {"34:1": ["a", "a2"]}
    assert all(m["modification"] == "" for m in out["members"])
