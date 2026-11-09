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
    assert a["verdict"] == {"q": 0.05, "equivalence_fc": 1.2} and a["allow_confounded"] is False
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
    bad = {"analysis_id": "a1", "kind": "compare", "reported_feature_ids": ["f1"]}
    assert _code(normalize_spec_v2, _spec(analyses=[bad])).code == "ANALYSIS_SPEC_INVALID"


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


def test_needs_input_passthrough_for_is_normalized_dataset_with_qc():
    ds = _ds()
    ds.quant_state = {"normalized": True, "method": "is", "unit": "pmol"}
    n = normalize_spec_v2(_spec())   # drift mode auto
    r = run_spec_v2(n, ds)
    assert r["status"] == "needs_input"
    assert r["rule_hash"] == rule_hash_v2(n) and r["spec"] == n


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
