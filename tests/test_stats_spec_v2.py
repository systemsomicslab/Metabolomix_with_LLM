"""analysis-spec.v2 の正規化と実行（spec §8・§12）。"""
import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from metabolomix.analysis.stats import spec_v2
from metabolomix.analysis.stats.spec_v2 import (ANALYSIS_SCHEMA_V2, RESULT_SCHEMA_V2,
                                                normalize_spec_v2, rule_hash_v2, run_spec_v2)
from metabolomix.core.atomic_io import DomainError

CMP = {"id": "T_vs_C", "factor": "group", "reference": "C", "test": "T"}


def _spec(**over):
    s = {"schema": "analysis-spec.v2", "comparisons": [dict(CMP)],
         "analyses": [{"analysis_id": "a1", "kind": "compare"}]}
    s.update(over)
    return s


def _code(fn, *a, **k):
    with pytest.raises(DomainError) as e:
        fn(*a, **k)
    return e.value


def test_defaults_filled_and_unknown_keys_rejected():
    n = normalize_spec_v2(_spec())
    assert n["schema"] == ANALYSIS_SCHEMA_V2 and n["input"] == {
        "source": "dataset", "exclude_samples": [], "exclude_features": [], "group_overrides": {}}
    assert n["design"] == {"factors": ["group"], "covariates": [], "block": None}
    assert n["preprocess"]["drift"]["mode"] == "auto"
    a = n["analyses"][0]
    assert a["test"] == "moderated_t" and a["effect_size"] == "log2_mean_difference"
    assert a["bh_population"] == "matrix" and a["min_detected_fraction"] == 0.5
    assert a["verdict"] == {"q": 0.05, "equivalence_fc": 1.2} and a["allow_confounded"] is True
    assert "p_adjust" not in a                       # 既定（bh）は書かない（rule_hash を前と同じに保つ）
    assert n["primary"] == "a1"
    assert normalize_spec_v2(copy.deepcopy(n)) == n
    assert rule_hash_v2(n) == rule_hash_v2(normalize_spec_v2(_spec()))
    for bad in (_spec(extra=1), _spec(design={"factors": ["group"], "x": 1}),
                _spec(analyses=[{"analysis_id": "a1", "kind": "compare", "foo": 1}]),
                _spec(analyses=[{"analysis_id": "a1", "kind": "pca", "test": "welch"}]),
                _spec(comparisons=[{**CMP, "zzz": 1}])):
        assert _code(normalize_spec_v2, bad).code == "ANALYSIS_SPEC_INVALID"
    e = _code(normalize_spec_v2, _spec(extra=1))
    assert "extra" in json.dumps(e.details)


def test_v1_schema_rejected_here():
    assert _code(normalize_spec_v2, _spec(schema="analysis-spec.v1")).code == "ANALYSIS_SPEC_INVALID"
    s = _spec()
    del s["schema"]
    assert _code(normalize_spec_v2, s).code == "ANALYSIS_SPEC_INVALID"


def test_primary_required_when_several_analyses():
    two = [{"analysis_id": "a1", "kind": "compare"}, {"analysis_id": "a2", "kind": "composition"}]
    assert _code(normalize_spec_v2, _spec(analyses=two)).code == "ANALYSIS_SPEC_INVALID"
    assert normalize_spec_v2(_spec(analyses=two, primary="a2"))["primary"] == "a2"
    assert _code(normalize_spec_v2, _spec(analyses=two, primary="nope")).code == "ANALYSIS_SPEC_INVALID"
    # pca を足しても primary の必須性は compare/composition の数で決まる
    one_pca = [{"analysis_id": "a1", "kind": "compare"}, {"analysis_id": "p", "kind": "pca"}]
    n = normalize_spec_v2(_spec(analyses=one_pca))
    assert n["primary"] == "a1" and n["analyses"][1]["n_components"] == 5
    assert _code(normalize_spec_v2, _spec(analyses=one_pca, primary="p")).code == "ANALYSIS_SPEC_INVALID"


def test_reported_requires_feature_ids():
    a = {"analysis_id": "a1", "kind": "compare", "bh_population": "reported"}
    assert _code(normalize_spec_v2, _spec(analyses=[a])).code == "ANALYSIS_SPEC_INVALID"
    ok = normalize_spec_v2(_spec(analyses=[{**a, "reported_feature_ids": ["f1"]}]))
    assert ok["analyses"][0]["reported_feature_ids"] == ["f1"]


@pytest.mark.parametrize("bh", [None, "matrix", "annotated"])
def test_reported_ids_with_other_population_are_dropped_with_a_note(bh):
    """reported 以外の母集団で書いた reported_feature_ids は効かない: 誤りにせず捨てて notes（dataset_compare と同じ）。
    捨てた ID は正規化した spec にも rule_hash にも入らない。"""
    a = {"analysis_id": "a1", "kind": "compare", **({"bh_population": bh} if bh else {})}
    notes: list[str] = []
    n = normalize_spec_v2(_spec(analyses=[{**a, "reported_feature_ids": ["f1"]}]), notes=notes)
    assert "reported_feature_ids" not in n["analyses"][0]
    assert len(notes) == 1 and "analyses[0].reported_feature_ids" in notes[0]
    assert rule_hash_v2(n) == rule_hash_v2(normalize_spec_v2(_spec(analyses=[a])))
    normalize_spec_v2(_spec(analyses=[a]), notes=(clean := []))
    assert clean == []


# main（p_adjust を足す前）の normalize_spec_v2 → rule_hash_v2 が返した値。p_adjust を足しても、既定（bh）の spec は
# 同じ hash でなければならない（同じ規則を前後で別の規則と読ませない）。値は git の main の実装で計算した。
_MAIN_RULE_HASH = {
    True: "bfa7ec48f6bae8287c6d07bac936040b141697b9aa60561758cff6979df46fff",
    False: "6cc266eb227fc93caab9dab98650ae22dc3753150eac22f9063c81598a9eb909",
}


@pytest.mark.parametrize("confounded", [True, False])
def test_default_p_adjust_keeps_the_rule_hash_from_before_p_adjust(confounded):
    base = {"analysis_id": "a1", "kind": "compare", "allow_confounded": confounded}
    omitted = normalize_spec_v2(_spec(analyses=[base]))
    explicit = normalize_spec_v2(_spec(analyses=[{**base, "p_adjust": "bh"}]))
    assert "p_adjust" not in omitted["analyses"][0] and "p_adjust" not in explicit["analyses"][0]
    assert rule_hash_v2(omitted) == rule_hash_v2(explicit) == _MAIN_RULE_HASH[confounded]
    holm = normalize_spec_v2(_spec(analyses=[{**base, "p_adjust": "holm"}]))
    assert holm["analyses"][0]["p_adjust"] == "holm"         # 既定でない方法は書き、hash も変わる
    assert rule_hash_v2(holm) != _MAIN_RULE_HASH[confounded]


def test_omitted_allow_confounded_hashes_as_the_new_default_true():
    """allow_confounded の既定を false → true に変えたのは規則の変更: 省いた spec は明示の true と同じ hash
    （main で省いた spec の hash＝明示の false の hash とは違う）。明示した true / false は main と同じ。"""
    omitted = rule_hash_v2(normalize_spec_v2(_spec()))
    assert omitted == _MAIN_RULE_HASH[True] != _MAIN_RULE_HASH[False]


def test_ids_and_factors_validated():
    dup = [{"analysis_id": "a", "kind": "pca"}, {"analysis_id": "a", "kind": "compare"}]
    assert _code(normalize_spec_v2, _spec(analyses=dup)).code == "ANALYSIS_SPEC_INVALID"
    assert _code(normalize_spec_v2, _spec(comparisons=[CMP, CMP])).code == "ANALYSIS_SPEC_INVALID"
    assert _code(normalize_spec_v2, _spec(comparisons=[{**CMP, "id": "a b"}])).code == "ANALYSIS_SPEC_INVALID"
    assert _code(normalize_spec_v2, _spec(comparisons=[{**CMP, "factor": "sex"}])).code == "ANALYSIS_SPEC_INVALID"
    assert _code(normalize_spec_v2, _spec(comparisons=[])).code == "ANALYSIS_SPEC_INVALID"


def _ds(n_per=6, qc_per_batch=8, with_subject=True):
    rng = np.random.default_rng(3)
    names, rows, order = [], [], 0
    for b in ("B1", "B2"):
        for i in range(qc_per_batch):
            order += 1
            names.append(f"QC_{b}_{i}")
            rows.append({"role": "qc", "group": "QC", "batch": b, "injection_order": order, "include": True})
        for g in ("C", "T"):
            for i in range(n_per):
                order += 1
                names.append(f"{g}_{b}_{i}")
                r = {"role": "sample", "group": g, "batch": b, "injection_order": order, "include": True}
                if with_subject:
                    r["subject_id"] = f"{g}{b}{i}"
                rows.append(r)
    nf = 12
    fm = np.exp2(rng.normal(20, 1, (nf, len(names))))
    for j, nm in enumerate(names):
        if nm.startswith("T_"):
            fm[:3, j] *= 4
    return SimpleNamespace(
        feature_matrix=fm, detected_mask=np.ones_like(fm, dtype=bool),
        sample_names=names, sample_ids=names, sample_metadata_rows=rows,
        feature_ids=[f"f{k}" for k in range(nf)],
        feature_metadata={f"f{k}": {"ontology": "PC" if k < 6 else "TG"} for k in range(nf)},
        feature_candidates={"f0": [{"name": "x"}], "f1": [{"name": "y"}]})


_FULL = {"preprocess": {"drift": {"mode": "off"}, "filters": {"blank_ratio": None, "d_ratio": None, "qc_rsd": None}},
         "analyses": [{"analysis_id": "cmp", "kind": "compare"},
                      {"analysis_id": "comp", "kind": "composition"},
                      {"analysis_id": "pc", "kind": "pca"}],
         "primary": "cmp"}


def test_end_to_end_on_fake_dataset_and_all_analyses_reported():
    n = normalize_spec_v2(_spec(**_FULL))
    r = run_spec_v2(n, _ds(), assay_kind="lipid")
    assert r["schema"] == RESULT_SCHEMA_V2 and r["spec"] == n and r["rule_hash"] == rule_hash_v2(n)
    assert set(r["input_fingerprint"]) == {"dataset_hash", "sample_table_hash"}
    assert len(r["input_fingerprint"]["sample_table_hash"]) == 64
    assert r["matrix_id"] and set(r["preprocess"]) == {"steps", "warnings", "subject_source"}
    assert r["primary"] == "cmp"
    assert [(a["analysis_id"], a["kind"]) for a in r["analyses"]] == [
        ("cmp", "compare"), ("comp", "composition"), ("pc", "pca")]
    assert all("result" in a for a in r["analyses"])
    assert r["analyses"][0]["result"]["schema"] == "stats-compare.v1"
    assert r["analyses"][1]["result"]["schema"] == "stats-composition.v1"
    assert [a["analysis_id"] for a in r["analyses"]] == [a["analysis_id"] for a in n["analyses"]]
    cmp_r = r["analyses"][0]["result"]
    comp_r = r["analyses"][1]["result"]
    defs = [(cmp_r, r["matrix_id"]), (comp_r["clr"], r["matrix_id"] + ":clr"),
            (comp_r["absolute"], r["matrix_id"] + ":class_total")]
    for res, mid in defs:
        assert res["effect_size_definition"]["name"] == "log2_mean_difference"
        assert res["test"] == "moderated_t"
        assert res["bh_population"]["name"] == "matrix"
        assert res["ci_level"] == 0.95
        assert res["matrix_id"] == mid
    pr = r["analyses"][2]["result"]
    assert len(pr["scores"]["units"]) == len(pr["unit_ids"]) > 0
    assert pr["explained_variance_ratio"]


def test_same_spec_same_dataset_same_result():
    n = normalize_spec_v2(_spec(**_FULL))
    ds = _ds()
    a = json.dumps(run_spec_v2(n, ds), sort_keys=True, default=str)
    b = json.dumps(run_spec_v2(n, ds), sort_keys=True, default=str)
    assert a == b


def test_is_normalized_dataset_with_qc_runs_without_drift_and_offers_on():
    """C3: 換算済み・QC が挟まった注入順でも drift=auto は尋ねず（needs_input を返さない）、補正なしで回して
    当てなかった理由と drift_mode="on" での呼び直しを前処理の段と warnings に残す。"""
    ds = _ds()
    qc_slots = {1, 3, 6, 9, 12, 15, 18, 20}           # 各バッチ 20 注入のうち QC の位置
    for b in range(2):
        rows = ds.sample_metadata_rows[20 * b:20 * (b + 1)]
        qc = iter(sorted(qc_slots))
        others = iter(o for o in range(1, 21) if o not in qc_slots)
        for row in rows:
            row["injection_order"] = 20 * b + (next(qc) if row["role"] == "qc" else next(others))
    ds.quant_state = {"normalized": True, "method": "is", "unit": "pmol"}
    n = normalize_spec_v2(_spec(preprocess={"filters": {"blank_ratio": None, "d_ratio": None, "qc_rsd": None}}))
    assert n["preprocess"]["drift"]["mode"] == "auto"
    r = run_spec_v2(n, ds)
    assert "status" not in r and r["rule_hash"] == rule_hash_v2(n) and r["spec"] == n
    drift = next(st for st in r["preprocess"]["steps"] if st["step"] == "drift")
    assert drift["status"] == "skipped" and drift["reason"] == "quant_normalized"
    assert drift["rerun_with"] == {"preprocess.drift.mode": "on"}


def test_rule_hash_ignores_key_order_but_tracks_values():
    a = _spec(**_FULL)
    b = json.loads(json.dumps(a))
    b = {k: b[k] for k in reversed(list(b))}
    b["analyses"] = [{k: x[k] for k in reversed(list(x))} for x in b["analyses"]]
    b["comparisons"] = [{k: x[k] for k in reversed(list(x))} for x in b["comparisons"]]
    b["preprocess"] = {k: b["preprocess"][k] for k in reversed(list(b["preprocess"]))}
    assert rule_hash_v2(normalize_spec_v2(a)) == rule_hash_v2(normalize_spec_v2(b))
    c = copy.deepcopy(a)
    c["analyses"][0]["verdict"] = {"q": 0.1}
    assert rule_hash_v2(normalize_spec_v2(a)) != rule_hash_v2(normalize_spec_v2(c))


def test_pca_knn_impute_and_interaction_comparison_normalize():
    n = normalize_spec_v2(_spec(analyses=[{"analysis_id": "p", "kind": "pca", "impute": "knn"}]))
    assert n["analyses"][0]["impute"] == "knn"
    inter = {"id": "ix", "kind": "interaction", "factors": ["group", "sex"],
             "levels": {"group": ["C", "T"], "sex": ["F", "M"]}}
    n = normalize_spec_v2(_spec(design={"factors": ["group", "sex"]}, comparisons=[inter]))
    assert n["comparisons"][0] == inter


def test_preprocess_error_details_are_flat():
    e = _code(normalize_spec_v2, _spec(preprocess={"drift": {"mode": "bogus"}}))
    assert e.code == "ANALYSIS_SPEC_INVALID" and "details" not in e.details
    assert e.details["field"].startswith("preprocess.")


# --- 補正の方法・PLS-DA（仕分けの一覧 2026-11-13 の B2-T21） ----------------------------------

def test_p_adjust_normalized_validated_and_reaches_engine():
    n = normalize_spec_v2(_spec(analyses=[{"analysis_id": "a1", "kind": "compare", "p_adjust": "holm"}]))
    assert n["analyses"][0]["p_adjust"] == "holm"
    assert rule_hash_v2(n) != rule_hash_v2(normalize_spec_v2(_spec()))
    e = _code(normalize_spec_v2, _spec(analyses=[{"analysis_id": "a1", "kind": "compare",
                                                  "p_adjust": "fdr"}]))
    assert e.code == "ANALYSIS_SPEC_INVALID" and e.details["field"] == "analyses[0].p_adjust"
    r = run_spec_v2(normalize_spec_v2(_spec(preprocess=_FULL["preprocess"], analyses=[
        {"analysis_id": "a1", "kind": "compare", "p_adjust": "bonferroni"}])), _ds())
    res = r["analyses"][0]["result"]
    assert res["p_adjust_method"] == "bonferroni"
    assert res["bh_population"]["q_interpretation"] == "fwer"


def test_allow_confounded_default_true_and_explicit_false_still_stops():
    n = normalize_spec_v2(_spec())
    assert n["analyses"][0]["allow_confounded"] is True
    off = normalize_spec_v2(_spec(analyses=[{"analysis_id": "a1", "kind": "compare",
                                             "allow_confounded": False}]))
    assert off["analyses"][0]["allow_confounded"] is False


def test_plsda_normalized_with_defaults_and_validated():
    n = normalize_spec_v2(_spec(analyses=[{"analysis_id": "pl", "kind": "plsda"}]))
    assert n["analyses"][0] == {"analysis_id": "pl", "kind": "plsda", "scaling": "autoscale",
                                "impute": "half_min", "n_components": 2, "cv_folds": 7,
                                "n_permutations": 200, "seed": 0}
    assert n["primary"] is None                       # plsda は primary にならない
    for bad in ({"n_components": 0}, {"cv_folds": 1}, {"n_permutations": 0}, {"seed": -1},
                {"scaling": "unit"}, {"impute": "zero"}, {"classes": ["C"]}, {"classes": ["C", "C"]},
                {"test": "welch"}):
        e = _code(normalize_spec_v2, _spec(analyses=[{"analysis_id": "pl", "kind": "plsda", **bad}]))
        assert e.code == "ANALYSIS_SPEC_INVALID", bad
    e = _code(normalize_spec_v2, _spec(analyses=[{"analysis_id": "pl", "kind": "plsda"}],
                                       primary="pl"))
    assert e.details["field"] == "primary"


def test_plsda_permutations_above_the_cap_are_lowered_with_a_note():
    """上限（10000）を超える n_permutations は断らず、上限に下げて notes に書く。spec と hash は回した回数。"""
    notes: list[str] = []
    n = normalize_spec_v2(_spec(analyses=[{"analysis_id": "pl", "kind": "plsda", "n_permutations": 50000}]),
                          notes=notes)
    assert n["analyses"][0]["n_permutations"] == 10000
    assert len(notes) == 1 and "n_permutations=50000" in notes[0] and "10000" in notes[0]
    capped = normalize_spec_v2(_spec(analyses=[{"analysis_id": "pl", "kind": "plsda", "n_permutations": 10000}]),
                               notes=(none := []))
    assert none == [] and rule_hash_v2(n) == rule_hash_v2(capped)


def test_inert_preprocess_keys_are_dropped_with_notes_and_keep_the_hash():
    """normalize.method=none の係数の特徴・drift.mode=off の loess_degree は効かない: 捨てて notes、hash は書かない spec と同じ。"""
    plain = {"drift": {"mode": "off"}, "normalize": {"method": "none"}}
    inert = {"drift": {"mode": "off", "loess_degree": 2},
             "normalize": {"method": "none", "features": ["f1"], "exclude_features": ["f2"]}}
    notes: list[str] = []
    n = normalize_spec_v2(_spec(preprocess=inert), notes=notes)
    assert n["preprocess"] == normalize_spec_v2(_spec(preprocess=plain))["preprocess"]
    assert rule_hash_v2(n) == rule_hash_v2(normalize_spec_v2(_spec(preprocess=plain)))
    assert len(notes) == 2 and all(x.startswith("preprocess: ") for x in notes)
    assert any("loess_degree" in x for x in notes) and any("normalization_features" in x for x in notes)


def test_plsda_runs_in_spec_and_returns_q2_permutation_and_vip():
    spec = _spec(preprocess=_FULL["preprocess"],
                 analyses=[{"analysis_id": "cmp", "kind": "compare"},
                           {"analysis_id": "pl", "kind": "plsda", "n_permutations": 20,
                            "classes": ["C", "T"]}])
    r = run_spec_v2(normalize_spec_v2(spec), _ds())
    assert [a["kind"] for a in r["analyses"]] == ["compare", "plsda"] and r["primary"] == "cmp"
    pl = r["analyses"][1]["result"]
    assert pl["schema"] == "stats-plsda.v1" and pl["n_components"] == 2
    assert len(pl["q2_cum"]) == 2 and pl["q2"] == pl["q2_cum"][-1]
    assert pl["permutation"]["n"] == 20 and 0 < pl["permutation"]["p"] <= 1
    assert {v["feature_id"] for v in pl["vip"]} <= {f"f{k}" for k in range(12)}
    assert pl["vip"][0]["feature_id"] in {"f0", "f1", "f2"}           # 仕込んだ差の特徴が上位
    again = run_spec_v2(normalize_spec_v2(spec), _ds())
    assert again["analyses"][1]["result"] == pl                         # seed で決まる
