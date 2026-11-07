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


def _prep(log2, n_a, *, classes=None, quant=None, feature_ids=None) -> PreparedMatrix:
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
        feature_classes=dict(CLASSES if classes is None else classes), annotated_features=frozenset())


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
    assert res["zero_handling"] == {"rule": "exclude_class_with_nonpositive_total",
                                    "excluded_classes": {"Cer": 2}}
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


def test_label_kind_follows_assay_kind():
    prep = _prep(_data(), 5)
    kinds = {k: composition(prep, assay_kind=k, factors=["group"], comparisons=BVA)["class_label_kind"]
             for k in ("lipid", "metabolite", "unknown")}
    assert kinds == {"lipid": "lipid_class", "metabolite": "msdial_ontology",
                     "unknown": "msdial_ontology"}


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
