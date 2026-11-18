"""クラスの組成 `composition`（統計エンジン spec §11）。PreparedMatrix はテスト内で直接組み立てる。"""
from __future__ import annotations

import numpy as np
from pytest import approx

from metabolomix.analysis.stats.composition import composition
from metabolomix.analysis.stats.data import QuantState
from metabolomix.analysis.stats.preprocess import PreparedMatrix, UnitInfo

BVA = [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"}]
FEATS = ["pc1", "pc2", "pc3", "tg1", "tg2", "tg3", "cer1", "cer2"]
CLASSES = {"pc1": "PC", "pc2": "PC", "pc3": "PC", "tg1": "TG", "tg2": "TG", "tg3": "TG",
           "cer1": "Cer", "cer2": "Cer"}


def _prep(log2, n_a, *, classes=None, quant=None, feature_ids=None, names=None, formulas=None) -> PreparedMatrix:
    log2 = np.asarray(log2, dtype=float)
    n, g = log2.shape
    feature_ids = feature_ids or FEATS[:g]
    groups = ["A"] * n_a + ["B"] * (n - n_a)
    units = [UnitInfo(unit_id=f"u{i:02d}", group=groups[i], batch="b1", subject_id=None,
                      injection_order=None, n_injections=1, sample_ids=[f"u{i:02d}_inj0"])
             for i in range(n)]
    return PreparedMatrix(
        matrix_id="mx", log2=log2, units=units, feature_ids=list(feature_ids), detected_units=None,
        qc_log2=np.empty((0, g)), qc_ids=[], steps=[], config={}, quant=quant, dataset_hash=None,
        subject_source="injection", warnings=[],
        feature_classes=dict(CLASSES if classes is None else classes), annotated_features=frozenset(),
        feature_names=dict(names or {}), feature_formulas=dict(formulas or {}))


def _data(n_a=5, n_b=5, *, factor=None, sd=0.05, seed=0):
    """全特徴 ~ log2 20。factor: {特徴 index: B 群への加算（log2）}、None なら加算なし。"""
    rng = np.random.default_rng(seed)
    x = rng.normal(20.0, sd, (n_a + n_b, len(FEATS)))
    for j, d in (factor or {}).items():
        x[n_a:, j] += d
    return x


def _verdicts(res):
    return {r["feature_id"]: r["verdict"] for r in res["rows"]}


def test_clr_rows_sum_to_zero_per_unit():
    prep = _prep(_data(), 5)
    res = composition(prep, assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["schema"] == "stats-composition.v1"
    assert res["classes"] == {"Cer": ["cer1", "cer2"], "PC": ["pc1", "pc2", "pc3"],
                              "TG": ["tg1", "tg2", "tg3"]}
    clr = res["clr"]
    assert clr["matrix_id"] == "mx:clr" and res["absolute"]["matrix_id"] == "mx:class_total"
    # 検算用に値を直接組み直す: 個体ごとの CLR 合計は 0
    x = prep.log2
    tot = np.stack([np.log2(np.exp2(x[:, [FEATS.index(f) for f in fs]]).sum(axis=1))
                    for fs in res["classes"].values()], axis=1)
    c = tot - tot.mean(axis=1, keepdims=True)
    assert c.sum(axis=1) == approx(np.zeros(10), abs=1e-9)
    # compare に渡った CLR と一致（A 群・B 群の平均差が log2fc）
    want = c[5:].mean(axis=0) - c[:5].mean(axis=0)
    got = {r["feature_id"]: r["log2fc"] for r in clr["rows"]}
    assert [got[k] for k in res["classes"]] == approx(list(want), abs=1e-9)


def test_class_with_zero_total_is_excluded_and_reported():
    x = _data()
    x[[1, 7], 6] = np.nan          # Cer は 2 個体で全メンバー欠測
    x[[1, 7], 7] = np.nan
    res = composition(_prep(x, 5), assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["zero_handling"] == {
        "rule": "exclude_class_with_nonpositive_total", "excluded_classes": {"Cer": 2},
        "member_missing_rule": "missing_member_counts_as_zero",
        "missing_member_cells": {"Cer": {"A": 2, "B": 2}, "PC": {"A": 0, "B": 0},
                                 "TG": {"A": 0, "B": 0}}}
    assert {r["feature_id"] for r in res["clr"]["rows"]} == {"PC", "TG"}
    # 絶対量には残り、欠測個体は n に数えない
    ab = {r["feature_id"]: r for r in res["absolute"]["rows"]}
    assert set(ab) == {"Cer", "PC", "TG"}
    assert ab["Cer"]["n_reference"] == 4 and ab["Cer"]["n_test"] == 4


def test_partial_missing_members_count_as_zero_in_total():
    x = _data()
    x[0, 6] = np.nan               # Cer の片方だけ欠測: 合計は残りの 1 つ
    res = composition(_prep(x, 5), assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["zero_handling"]["excluded_classes"] == {}
    assert {r["feature_id"] for r in res["clr"]["rows"]} == {"Cer", "PC", "TG"}
    # 欠測のメンバーを 0 と数えたことを開示する（spec §11）
    zh = res["zero_handling"]
    assert zh["member_missing_rule"] == "missing_member_counts_as_zero"
    assert zh["missing_member_cells"]["Cer"] == {"A": 1, "B": 0}
    assert zh["missing_member_cells"]["PC"] == {"A": 0, "B": 0}
    warn = [w for w in res["warnings"] if w["code"] == "class_total_includes_missing_members"]
    assert len(warn) == 1 and warn[0]["classes"] == ["Cer"]


def test_no_missing_members_no_disclosure_warning():
    res = composition(_prep(_data(), 5), assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["zero_handling"]["missing_member_cells"]["Cer"] == {"A": 0, "B": 0}
    assert not any(w["code"] == "class_total_includes_missing_members" for w in res["warnings"])


def test_no_feature_classes_warns_and_returns_nothing():
    res = composition(_prep(_data(), 5, classes={}), assay_kind="lipid", factors=["group"],
                      comparisons=BVA)
    assert res["clr"] is None and res["absolute"] is None
    assert res["classes"] == {} and res["unclassified_features"] == len(FEATS)
    assert any(w["code"] == "no_feature_classes" for w in res["warnings"])


def test_label_kind_follows_assay_kind():
    prep = _prep(_data(), 5)
    kinds = {k: composition(prep, assay_kind=k, factors=["group"], comparisons=BVA)["class_label_kind"]
             for k in ("lipid", "metabolite", "unknown")}
    assert kinds == {"lipid": "lipid_class", "metabolite": "msdial_ontology",
                     "unknown": "msdial_ontology"}


def test_assay_kind_changes_only_the_label_when_all_names_are_lipids():
    # クラスはいつも ontology（feature_classes）で切る。MS-DIAL の脂質クラスの脂質名（と名前の無い特徴）だけなら、
    # 種別が変えるのはラベルと non_lipid・lipid_class_relabel の入れ物（lipid では 0 件）だけ（reanalysis.md §2.1）
    import json

    names = {"pc1": "PC 34:1", "pc2": "w/o MS2: PC 16:0_18:1", "tg1": "TG 52:2"}
    for prep in (_prep(_data(), 5), _prep(_data(), 5, names=names)):
        runs = {k: composition(prep, assay_kind=k, factors=["group"], comparisons=BVA)
                for k in ("lipid", "metabolite", "unknown")}
        blocks = ("class_label_kind", "non_lipid", "lipid_class_relabel")
        bodies = {k: json.dumps({key: v for key, v in r.items() if key not in blocks},
                                sort_keys=True, default=str)
                  for k, r in runs.items()}
        assert len(set(bodies.values())) == 1
        assert runs["lipid"]["classes"] and runs["lipid"]["classes"] == runs["metabolite"]["classes"]
        assert runs["lipid"]["non_lipid"]["n_features"] == 0
        assert runs["lipid"]["lipid_class_relabel"]["n_features"] == 0
        for k in ("metabolite", "unknown"):
            assert runs[k]["non_lipid"] is None and runs[k]["lipid_class_relabel"] is None


# 研究室 MSP の天然物・医薬品の名前が付いた特徴（脂質名の文法に乗らない）を足した行列。脂質の側には、MS-DIAL の
# 脂質クラスで GOSLIN が読めない略記（Cer_AP の `(2OH)`）と、ClassyFire 系のクラスで GOSLIN が読む慣用名を置く
MIXED_FEATS = FEATS + ["sap1", "sap2", "dep1", "noname", "cerap", "ole"]
MIXED_CLASSES = {**CLASSES, "sap1": "Triterpene saponins", "sap2": "Triterpene saponins",
                 "dep1": "Cyclic depsipeptides", "noname": "PC", "cerap": "Cer_AP",
                 "ole": "Long-chain fatty acids"}
MIXED_NAMES = {"pc1": "PC 34:1", "pc2": "PC 16:0_18:1", "pc3": "low score: PC 36:2", "tg1": "TG 52:2",
               "tg2": "TG 16:0_18:1_18:2", "tg3": "TG 54:3", "cer1": "Cer 18:1;O2/16:0", "cer2": "Cer d18:1/24:0",
               "sap1": "Ginsenoside Rb1", "sap2": "w/o MS2: Soyasaponin I", "dep1": "Enniatin B",
               "noname": "Unknown", "cerap": "Cer 12:0;O3/8:0(2OH)", "ole": "Oleic acid"}
NON_LIPID = ("sap1", "sap2", "dep1")
# 脂質アッセイでは `Cer`（MS-DIAL の脂質クラスでない）の 2 つは名前の頭部から `Cer_NS` へ移る
LIPID_CLASSES = {"Cer_NS", "Cer_AP", "Long-chain fatty acids", "PC", "TG"}
ONTOLOGY_CLASSES = {"Cer", "Cer_AP", "Long-chain fatty acids", "PC", "TG"}


def _mixed(seed=0, factor=None):
    rng = np.random.default_rng(seed)
    x = rng.normal(20.0, 0.05, (10, len(MIXED_FEATS)))
    for j, d in (factor or {}).items():
        x[5:, j] += d
    return _prep(x, 5, classes=MIXED_CLASSES, feature_ids=MIXED_FEATS, names=MIXED_NAMES)


def test_lipid_assay_puts_unparsable_names_in_non_lipid_and_keeps_them_out_of_classes():
    # サポニンの量が B 群で 8 倍: 脂質クラスの数字には効かない
    prep = _mixed(factor={MIXED_FEATS.index("sap1"): 3.0, MIXED_FEATS.index("sap2"): 3.0})
    res = composition(prep, assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["class_label_kind"] == "lipid_class"
    # 名前の無い特徴（Unknown）は判定せず、クラス（PC）に数える
    assert res["classes"] == {"Cer_NS": ["cer1", "cer2"], "Cer_AP": ["cerap"], "Long-chain fatty acids": ["ole"],
                              "PC": ["pc1", "pc2", "pc3", "noname"], "TG": ["tg1", "tg2", "tg3"]}
    assert res["lipid_class_relabel"] == {
        "rule": "msdial_lipid_class_from_name_or_class_formula", "n_features": 2, "n_rescued": 0,
        "moves": [{"from": "Cer", "to": "Cer_NS", "route": "name", "n_features": 2, "n_rescued": 0,
                   "example": {"feature_id": "cer1", "name": "Cer 18:1;O2/16:0"}}]}
    nl = res["non_lipid"]
    assert nl["label"] == "non-lipid" and nl["rule"] == "lipid_annotation_unmatched"
    assert nl["n_features"] == 3 and nl["feature_ids"] == ["sap1", "sap2", "dep1"]
    assert nl["by_ontology"] == {"Triterpene saponins": 2, "Cyclic depsipeptides": 1}
    assert list(nl["by_ontology"]) == ["Triterpene saponins", "Cyclic depsipeptides"]   # 多い順
    assert nl["examples"] == [
        {"feature_id": "sap1", "name": "Ginsenoside Rb1", "ontology": "Triterpene saponins"},
        {"feature_id": "dep1", "name": "Enniatin B", "ontology": "Cyclic depsipeptides"}]
    # 非脂質はクラスの合計・CLR・欠測の数えのどれにも出ない
    for table in ("absolute", "clr"):
        assert {r["feature_id"] for r in res[table]["rows"]} == LIPID_CLASSES
    assert set(res["zero_handling"]["missing_member_cells"]) == LIPID_CLASSES
    assert res["unclassified_features"] == 0
    # 非脂質を除いた行列で計算したものと同じ数字（脂質クラスの数字は非脂質の値に左右されない）
    lipid_only = [f for f in MIXED_FEATS if f not in NON_LIPID]
    cols = [MIXED_FEATS.index(f) for f in lipid_only]
    ref = composition(_prep(prep.log2[:, cols], 5, classes={f: MIXED_CLASSES[f] for f in lipid_only},
                            feature_ids=lipid_only, names={f: MIXED_NAMES[f] for f in lipid_only}),
                      assay_kind="lipid", factors=["group"], comparisons=BVA)
    for table in ("absolute", "clr"):
        got = {r["feature_id"]: r["log2fc"] for r in res[table]["rows"]}
        want = {r["feature_id"]: r["log2fc"] for r in ref[table]["rows"]}
        assert got == approx(want, abs=1e-12)


def test_metabolite_assay_keeps_every_ontology_as_a_class():
    res = composition(_mixed(), assay_kind="metabolite", factors=["group"], comparisons=BVA)
    assert res["non_lipid"] is None and res["lipid_class_relabel"] is None
    assert set(res["classes"]) == ONTOLOGY_CLASSES | {"Triterpene saponins", "Cyclic depsipeptides"}
    assert {r["feature_id"] for r in res["absolute"]["rows"]} == set(res["classes"])


def test_only_non_lipid_classes_gives_no_composition_and_says_why():
    prep = _mixed()
    keep = list(NON_LIPID)
    cols = [MIXED_FEATS.index(f) for f in keep]
    only = _prep(prep.log2[:, cols], 5, classes={f: MIXED_CLASSES[f] for f in keep}, feature_ids=keep,
                 names={f: MIXED_NAMES[f] for f in keep})
    res = composition(only, assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["classes"] == {} and res["absolute"] is None and res["clr"] is None
    assert res["non_lipid"]["n_features"] == 3
    warn = next(w for w in res["warnings"] if w["code"] == "no_feature_classes")
    assert "non_lipid" in warn["message"] and "3" in warn["message"]


def test_lipid_grammar_unavailable_keeps_classes_and_warns(monkeypatch):
    from metabolomix.analysis.stats import composition as comp_mod

    from metabolomix.msdial.lipid_identity import LipidCall

    monkeypatch.setattr(comp_mod, "lipid_annotation_class",
                        lambda name, ontology, formula=None: LipidCall(None, None, "grammar_unavailable", False))
    res = comp_mod.composition(_mixed(), assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert res["non_lipid"] is None and res["lipid_class_relabel"] is None
    assert "Triterpene saponins" in res["classes"] and "Cer" in res["classes"]   # 外さず、移さない
    assert any(w["code"] == "lipid_grammar_unavailable" for w in res["warnings"])


# 研究室 MSP の脂質で名前が IUPAC 名・ID のもの（規則 3 で ClassyFire のクラスと組成式から救う）と、内部標準の名前
RESCUE_FEATS = ["pc1", "pc2", "tg1", "pcx", "dgx", "badpc", "sm", "ole"]
RESCUE_CLASSES = {"pc1": "PC", "pc2": "PC", "tg1": "TG", "pcx": "Phosphatidylcholines", "dgx": "1,2-diacylglycerols",
                  "badpc": "Phosphatidylcholines", "sm": "Sphingomyelins", "ole": "Long-chain fatty acids"}
RESCUE_NAMES = {"pc1": "PC 34:1", "pc2": "PC 16:0_18:1", "tg1": "TG 52:2",
                "pcx": "no MS2: 1,2-Dioleoyl-sn-glycero-3-phosphocholine", "dgx": "1,2-Dioleoyl-sn-glycerol",
                "badpc": "7.alpha.,24(S)-Dihydroxy-4-cholesten-3-one", "sm": "SM d18:1-18:1(d9)",
                "ole": "Oleic acid"}
RESCUE_FORMULAS = {"pcx": "C44H84NO8P", "dgx": "C39H72O5", "badpc": "C27H44O3", "ole": "C18H34O2"}


def _rescue(seed=0, factor=None):
    rng = np.random.default_rng(seed)
    x = rng.normal(20.0, 0.05, (10, len(RESCUE_FEATS)))
    for j, d in (factor or {}).items():
        x[5:, j] += d
    return _prep(x, 5, classes=RESCUE_CLASSES, feature_ids=RESCUE_FEATS, names=RESCUE_NAMES,
                 formulas=RESCUE_FORMULAS)


def test_lipid_assay_rescues_lab_msp_lipids_by_class_and_formula_into_the_msdial_class():
    prep = _rescue(factor={RESCUE_FEATS.index("pcx"): 2.0})
    res = composition(prep, assay_kind="lipid", factors=["group"], comparisons=BVA)
    # IUPAC 名の PC は LBM の PC と同じクラスに入り、DG は DG に。内部標準の SM は名前の頭部から SM に。
    # 慣用名の脂肪酸は名前で読め、クラスと組成式から FA に。組成式がクラスの形に合わない行は非脂質
    assert res["classes"] == {"DG": ["dgx"], "FA": ["ole"], "PC": ["pc1", "pc2", "pcx"], "SM": ["sm"], "TG": ["tg1"]}
    assert res["non_lipid"]["feature_ids"] == ["badpc"]
    assert res["non_lipid"]["by_ontology"] == {"Phosphatidylcholines": 1}
    relabel = res["lipid_class_relabel"]
    assert relabel["n_features"] == 4 and relabel["n_rescued"] == 2
    assert [(m["from"], m["to"], m["route"], m["n_features"], m["n_rescued"]) for m in relabel["moves"]] == [
        ("1,2-diacylglycerols", "DG", "class_formula", 1, 1),
        ("Long-chain fatty acids", "FA", "class_formula", 1, 0),
        ("Phosphatidylcholines", "PC", "class_formula", 1, 1),
        ("Sphingomyelins", "SM", "name", 1, 0)]
    # PC の合計は pc1・pc2・pcx の和: クラスを手で移した行列（代謝物のアッセイ＝移さない）と同じ数字
    moved = {**RESCUE_CLASSES, "pcx": "PC", "dgx": "DG", "sm": "SM", "ole": "FA"}
    keep = [f for f in RESCUE_FEATS if f != "badpc"]
    cols = [RESCUE_FEATS.index(f) for f in keep]
    ref = composition(_prep(prep.log2[:, cols], 5, classes={f: moved[f] for f in keep}, feature_ids=keep),
                      assay_kind="metabolite", factors=["group"], comparisons=BVA)
    for table in ("absolute", "clr"):
        got = {r["feature_id"]: r["log2fc"] for r in res[table]["rows"]}
        want = {r["feature_id"]: r["log2fc"] for r in ref[table]["rows"]}
        assert got == approx(want, abs=1e-12)
    # 注釈そのものは変えない
    assert prep.feature_classes["pcx"] == "Phosphatidylcholines"


def test_metabolite_assay_does_not_move_classyfire_classes():
    res = composition(_rescue(), assay_kind="metabolite", factors=["group"], comparisons=BVA)
    assert res["lipid_class_relabel"] is None and res["non_lipid"] is None
    assert res["classes"]["Phosphatidylcholines"] == ["pcx", "badpc"]
    assert res["classes"]["PC"] == ["pc1", "pc2"]


def test_composition_shift_detected_by_clr_not_by_absolute():
    # 全特徴が B 群で厳密に 2 倍（log2 で +1）: 総量だけが動く
    prep = _prep(_data(factor={j: 1.0 for j in range(8)}), 5)
    res = composition(prep, assay_kind="lipid", factors=["group"], comparisons=BVA)
    assert set(_verdicts(res["absolute"]).values()) == {"changed"}
    assert "changed" not in _verdicts(res["clr"]).values()
    for r in res["absolute"]["rows"]:
        assert r["log2fc"] == approx(1.0, abs=0.2)


def test_class_share_change_detected_by_clr():
    # TG のメンバーだけが 4 倍: TG の取り分が増え、他は相対的に減る
    prep = _prep(_data(factor={3: 2.0, 4: 2.0, 5: 2.0}), 5)
    res = composition(prep, assay_kind="lipid", factors=["group"], comparisons=BVA)
    v = _verdicts(res["clr"])
    assert v["TG"] == "changed"
    fc = {r["feature_id"]: r["log2fc"] for r in res["clr"]["rows"]}
    assert fc["TG"] > 0 > fc["PC"]


def test_absolute_analysis_reports_unit_from_quant_state():
    q = QuantState(normalized=True, method="internal_standard", unit="nmol/mL")
    res = composition(_prep(_data(), 5, quant=q), assay_kind="lipid", factors=["group"],
                      comparisons=BVA)
    assert res["unit"] == "nmol/mL"
    assert composition(_prep(_data(), 5), assay_kind="lipid", factors=["group"],
                       comparisons=BVA)["unit"] is None


def test_unclassified_features_counted_and_ignored():
    classes = {k: v for k, v in CLASSES.items() if k != "cer2"}
    res = composition(_prep(_data(), 5, classes=classes), assay_kind="lipid", factors=["group"],
                      comparisons=BVA)
    assert res["unclassified_features"] == 1
    assert res["classes"]["Cer"] == ["cer1"]


def test_clr_not_computable_with_fewer_than_two_classes():
    classes = {k: "PC" for k in FEATS}
    res = composition(_prep(_data(), 5, classes=classes), assay_kind="lipid", factors=["group"],
                      comparisons=BVA)
    assert res["clr"] is None
    assert any(w["code"] == "clr_not_computable" for w in res["warnings"])
    assert res["absolute"] is not None


def test_bh_population_is_forced_to_matrix():
    res = composition(_prep(_data(), 5), assay_kind="lipid", factors=["group"], comparisons=BVA,
                      bh_population="annotated")
    assert res["clr"]["bh_population"]["name"] == "matrix"
    assert res["absolute"]["bh_population"]["name"] == "matrix"
    assert any(w["code"] == "bh_population_overridden" for w in res["warnings"])
