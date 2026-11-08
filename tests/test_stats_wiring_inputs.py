"""配線向けの受け口（除外・群の上書き・注入の列・注入ごとの値・平均・割れのまとめ）。"""
from __future__ import annotations

import copy
import json
import math

import numpy as np
import pytest

from metabolomix.analysis.stats import data as D
from metabolomix.analysis.stats import preprocess as P
from metabolomix.analysis.stats.compare import compare
from metabolomix.analysis.stats.pca import pca
from metabolomix.analysis.stats.preprocess import PreparedMatrix, UnitInfo
from metabolomix.analysis.stats.spec_v2 import (RESULT_SCHEMA_V2, normalize_spec_v2, prepare_spec_v2,
                                                run_spec_v2)
from metabolomix.core.atomic_io import DomainError
from tests.test_stats_preprocess import OFF_FILTERS, _study
from tests.test_stats_spec_v2 import _FULL, _ds, _spec

OFF = {"drift": {"mode": "off"}, "filters": OFF_FILTERS, "normalize": {"method": "none"}}
BVA = [{"id": "T_vs_C", "factor": "group", "reference": "C", "test": "T"}]


def _code(fn, *a, **k):
    with pytest.raises(DomainError) as e:
        fn(*a, **k)
    return e.value


def _pre(study, **cfg):
    return P.preprocess(study, {**OFF, **cfg})


def _dump(x) -> str:
    return json.dumps(x, default=str, sort_keys=True)


# --- from_dataset: 除外・上書き ---------------------------------------------

def test_exclusions_drop_samples_and_features_and_change_matrix_id():
    ds = _ds()
    base = D.from_dataset(ds)
    st = D.from_dataset(ds, exclude_samples=["C_B1_0", "T_B2_1"], exclude_features=["f3", "f5"])
    assert [s.include for s in st.samples].count(False) == 2
    assert not next(s for s in st.samples if s.sample_id == "C_B1_0").include
    assert st.values.shape == (len(base.sample_ids), len(base.feature_ids) - 2)
    assert st.detected.shape == st.values.shape
    assert st.feature_ids == [f for f in base.feature_ids if f not in ("f3", "f5")]
    assert all(s.include for s in base.samples) and len(ds.feature_ids) == 12  # ds は変わらない
    prep = _pre(st)
    assert not any("C_B1_0" in u.sample_ids for u in prep.units)
    assert _pre(base).matrix_id != prep.matrix_id
    # 名前でも ID でも引ける（ID != 表示名のとき）
    ds2 = copy.copy(ds)
    ds2.sample_ids = [f"id{i}" for i in range(len(ds.sample_names))]
    pos = ds.sample_names.index("C_B1_0")
    a = D.from_dataset(ds2, exclude_samples=["C_B1_0"])
    b = D.from_dataset(ds2, exclude_samples=[f"id{pos}"])
    assert [s.include for s in a.samples] == [s.include for s in b.samples]
    assert a.samples[pos].include is False
    for kw in ({"exclude_samples": ["nope"]}, {"exclude_features": ["nope"]},
               {"group_overrides": {"nope": "X"}}):
        e = _code(D.from_dataset, ds, **kw)
        assert e.code == "STUDY_DATA_INVALID" and "nope" in json.dumps(e.details)


def test_group_overrides_replace_sample_table_group():
    ds = _ds()
    st = D.from_dataset(ds, group_overrides={"C_B1_0": "T"})
    assert next(s for s in st.samples if s.sample_id == "C_B1_0").group == "T"
    assert next(s for s in st.samples if s.sample_id == "C_B1_1").group == "C"
    assert D.from_dataset(ds).samples[ds.sample_names.index("C_B1_0")].group == "C"
    assert _code(D.from_dataset, ds, group_overrides={"C_B1_0": " "}).code == "STUDY_DATA_INVALID"


def test_matrix_id_changes_when_only_feature_list_changes():
    ds = _ds()
    a = D.from_dataset(ds)
    b = D.from_dataset(ds, exclude_features=["f11"])
    assert a.dataset_hash == b.dataset_hash
    assert _pre(a).matrix_id != _pre(b).matrix_id
    assert _pre(a).matrix_id == _pre(D.from_dataset(ds)).matrix_id


# --- 注入の列・割れた単位の注意 ---------------------------------------------

def _rep_study(attr_key="run_day"):
    rows = []
    for s, g in (("S0", "C"), ("S1", "C"), ("S2", "T"), ("S3", "T")):
        for r in range(2):
            rows.append((f"{s}_{r}", "sample", g, None, None, s, True, {attr_key: f"d{r}"}))
    vals = 1000.0 * (1 + 0.1 * np.random.default_rng(1).standard_normal((len(rows), 3)))
    return _study(vals, rows)


def test_injection_columns_keep_replicates_in_one_unit():
    st = _rep_study()
    plain = _pre(st)
    assert len(plain.units) == 8
    prep = _pre(st, injection_columns=["run_day"])
    assert len(prep.units) == 4 and all(u.n_injections == 2 for u in prep.units)
    assert not [w for w in prep.warnings if w["code"] == "units_split_by_attribute"]
    assert prep.config["injection_columns"] == ["run_day"]
    assert plain.matrix_id != prep.matrix_id
    for bad in ([1], [""], "run_day", None):
        e = _code(P.normalize_preprocess_config, {"injection_columns": bad})
        assert e.code == "PREPROCESS_CONFIG_INVALID"
    assert P.PREPROCESS_DEFAULTS["injection_columns"] == []


def test_attribute_splitting_a_subject_warns_units_split_by_attribute():
    prep = _pre(_rep_study())
    w = [w for w in prep.warnings if w["code"] == "units_split_by_attribute"]
    assert len(w) == 1 and w[0]["attributes"] == ["run_day"] and w[0]["n_subjects"] == 4
    assert "injection_columns" in w[0]["message"]
    assert w[0]["suggestion"] == {"preprocess.injection_columns": ["run_day"]}
    # collapse を切ったときは「分かれた」のではない（単位は注入のまま）ので警告しない
    off = _pre(_rep_study(), collapse_replicates=False)
    assert not [x for x in off.warnings if x["code"] == "units_split_by_attribute"]
    # 長い一覧は 20 で切る
    rows = []
    for k in range(25):
        for r in range(2):
            rows.append((f"S{k}_{r}", "sample", "C" if k % 2 else "T", None, None, f"S{k}", True,
                         {"run_day": f"d{r}"}))
    vals = np.full((len(rows), 2), 100.0) + np.arange(len(rows))[:, None]
    bw = next(w for w in _pre(_study(vals, rows)).warnings if w["code"] == "units_split_by_attribute")
    assert bw["n_subjects"] == 25 and len(bw["subjects"]) == 20 and "ほか 5 件" in bw["message"]


# --- 注入ごとの値 -------------------------------------------------------------

def test_prepared_sample_values_are_pre_log2_and_match_log2():
    rows = [(f"a{i}", "sample", "C", None, None, None) for i in range(3)]
    rows += [(f"b{i}", "sample", "T", None, None, None) for i in range(3)]
    rows += [(f"q{i}", "qc", None, None, None, None) for i in range(3)]
    rows += [("blk", "blank", None, None, None, None)]
    rng = np.random.default_rng(2)
    vals = 1000.0 * (1 + 0.1 * rng.standard_normal((len(rows), 4)))
    vals[0, 1] = 0.0
    vals[1, 2] = np.nan
    vals[2, 3] = -5.0
    det = np.ones(vals.shape, bool)
    det[3, 0] = False
    prep = P.preprocess(_study(vals, rows, detected=det), {**OFF, "normalize": {"method": "tic"}})
    assert prep.sample_ids == [f"a{i}" for i in range(3)] + [f"b{i}" for i in range(3)]
    assert prep.sample_values.shape == (6, len(prep.feature_ids))
    with np.errstate(all="ignore"):
        expect = np.where(prep.sample_values > 0, np.log2(prep.sample_values), np.nan)
    np.testing.assert_allclose(prep.log2, expect, equal_nan=True)  # 個体 = 注入（subject 無し）
    assert math.isnan(prep.sample_values[1, 2]) and prep.sample_values[2, 3] < 0
    assert prep.sample_values[0, 1] == 0.0  # 0 以下は欠測にせずそのまま
    assert prep.sample_detected.shape == prep.sample_values.shape
    assert not prep.sample_detected[3, 0] and prep.sample_detected[0, 0]
    assert P.preprocess(_study(vals, rows), OFF).sample_detected is None
    assert not np.allclose(prep.sample_values[:, 0], vals[:6, 0])  # tic の係数が掛かった後の値


# --- 比較の行 -------------------------------------------------------------------

def test_compare_rows_carry_linear_means_for_simple_contrasts():
    prep = _pre(D.from_dataset(_ds()))
    res = compare(prep, factors=["group"], comparisons=BVA)
    j = prep.feature_ids.index("f0")
    ref = [i for i, u in enumerate(prep.units) if u.group == "C"]
    tst = [i for i, u in enumerate(prep.units) if u.group == "T"]
    row = next(r for r in res["rows"] if r["feature_id"] == "f0")
    assert row["mean_reference"] == pytest.approx(float(np.mean(2.0 ** prep.log2[ref, j])))
    assert row["mean_test"] == pytest.approx(float(np.mean(2.0 ** prep.log2[tst, j])))
    log2 = prep.log2.copy()
    log2[ref[0], j] = np.nan
    log2[:, 1] = np.nan
    p2 = copy.copy(prep)
    p2.log2 = log2
    r2 = compare(p2, factors=["group"], comparisons=BVA)
    row = next(r for r in r2["rows"] if r["feature_id"] == "f0")
    assert row["mean_reference"] == pytest.approx(float(np.mean(2.0 ** log2[ref[1:], j])))
    none_row = next(r for r in r2["rows"] if r["feature_id"] == prep.feature_ids[1])
    assert none_row["mean_reference"] is None and none_row["mean_test"] is None


def _matrix(log2, units) -> PreparedMatrix:
    g = log2.shape[1]
    return PreparedMatrix(matrix_id="mx", log2=log2, units=units,
                          feature_ids=[f"f{j}" for j in range(g)], detected_units=None,
                          qc_log2=np.empty((0, g)), qc_ids=[], steps=[], config={}, quant=None,
                          dataset_hash=None, subject_source="injection", warnings=[],
                          feature_classes={}, annotated_features=frozenset())


def _unit(uid, group, **attrs) -> UnitInfo:
    return UnitInfo(unit_id=uid, group=group, batch=None, subject_id=None, injection_order=None,
                    n_injections=1, sample_ids=[uid], attributes=attrs)


def test_interaction_rows_have_no_means():
    rng = np.random.default_rng(5)
    units = [_unit(f"{g}{s}{i}", g, sex=s) for g in ("C", "T") for s in ("M", "F") for i in range(4)]
    prep = _matrix(rng.normal(10, 1, (len(units), 3)), units)
    comps = [{"id": "ix", "kind": "interaction", "factors": ["group", "sex"],
              "levels": {"group": ["C", "T"], "sex": ["M", "F"]}},
             {"id": "s", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "M"}}]
    res = compare(prep, factors=["group", "sex"], comparisons=comps)
    ix = [r for r in res["rows"] if r["contrast_id"] == "ix"]
    assert ix and all(r["mean_reference"] is None and r["mean_test"] is None for r in ix)
    s = [r for r in res["rows"] if r["contrast_id"] == "s"]
    assert s and all(r["mean_reference"] is not None for r in s)


def _two(a, b, extra=None) -> PreparedMatrix:
    units = [_unit(f"u{i}", "C" if i < len(a) else "T") for i in range(len(a) + len(b))]
    cols = [np.concatenate([a, b])] + ([extra] if extra is not None else [])
    return _matrix(np.column_stack(cols), units)


def test_welch_near_constant_groups_are_not_tested_degenerate_variance():
    rng = np.random.default_rng(7)
    a = np.full(4, 10.0)
    b = np.full(4, 12.0) + 1e-15 * rng.standard_normal(4)
    ok = rng.normal(10, 1, 8)
    for test in ("welch", "student"):
        res = compare(_two(a, b, ok), factors=["group"], comparisons=BVA, test=test)
        r0, r1 = res["rows"]
        assert r0["status"] == "not_tested" and r0["reason"] == "degenerate_variance"
        assert r0["verdict"] == "inconclusive" and r0["p"] is None
        assert r1["status"] == "tested"
    res = compare(_two(np.full(4, 10.0), rng.normal(12, 1, 4)), factors=["group"],
                  comparisons=BVA, test="welch")
    assert res["rows"][0]["status"] == "tested"  # 片方の群だけが定数なら検定する


def test_pca_knn_failure_becomes_pca_imputation_failed():
    rng = np.random.default_rng(1)
    log2 = rng.normal(10, 1, (6, 10))
    log2[0, :9] = np.nan  # 1 個体が 90% 欠測: knn の colmax=0.8 を超える
    prep = _matrix(log2, [_unit(f"u{i}", "C" if i < 3 else "T") for i in range(6)])
    e = _code(pca, prep, impute="knn")
    assert e.code == "PCA_IMPUTATION_FAILED" and "half_min" in e.message
    assert e.details["impute"] == "knn" and e.details["cause"]["code"] == "RECIPE_PRECONDITION"
    assert pca(prep, impute="half_min", n_components=2)["imputation"]["used"] is True


# --- spec v2 --------------------------------------------------------------------

def _norm(**over):
    return normalize_spec_v2(_spec(**{**_FULL, **over}))


def test_spec_v2_input_defaults_and_exclusions_reproduce():
    n = _norm()
    assert n["input"] == {"source": "dataset", "exclude_samples": [], "exclude_features": [],
                          "group_overrides": {}}
    assert normalize_spec_v2(copy.deepcopy(n)) == n
    ds = _ds()
    r1 = run_spec_v2(n, ds, assay_kind="lipid")
    assert _dump(r1) == _dump(run_spec_v2(copy.deepcopy(n), ds, assay_kind="lipid"))
    ex = _norm(input={"exclude_samples": ["C_B1_0", "C_B1_0"], "exclude_features": ["f11"],
                      "group_overrides": {"T_B1_0": "C"}})
    assert ex["input"]["exclude_samples"] == ["C_B1_0"]
    rx = run_spec_v2(ex, ds, assay_kind="lipid")
    assert rx["matrix_id"] != r1["matrix_id"] and rx["rule_hash"] != r1["rule_hash"]
    assert _dump(rx) == _dump(run_spec_v2(copy.deepcopy(ex), ds, assay_kind="lipid"))
    rx_cmp = next(a for a in rx["analyses"] if a["analysis_id"] == "cmp")["result"]
    r1_cmp = next(a for a in r1["analyses"] if a["analysis_id"] == "cmp")["result"]
    assert "f11" not in {r["feature_id"] for r in rx_cmp["rows"]}
    # C から 1 本除外し、T の 1 本を C へ移した: C は同数、T は 1 本減る（除外と上書きが数字に効く）
    assert rx_cmp["group_n"]["C"] == r1_cmp["group_n"]["C"]
    assert rx_cmp["group_n"]["T"] == r1_cmp["group_n"]["T"] - 1
    assert _code(normalize_spec_v2, _spec(input={"exclude_samples": "x"})).code == "ANALYSIS_SPEC_INVALID"
    assert _code(normalize_spec_v2, _spec(input={"foo": 1})).code == "ANALYSIS_SPEC_INVALID"
    assert _code(run_spec_v2, _norm(input={"exclude_samples": ["zzz"]}), ds).code == "STUDY_DATA_INVALID"


def test_run_with_prepared_equals_run_without():
    n = _norm()
    ds = _ds()
    prepared = prepare_spec_v2(n, ds)
    assert isinstance(prepared, PreparedMatrix)
    a = run_spec_v2(n, ds, assay_kind="lipid")
    b = run_spec_v2(n, ds, assay_kind="lipid", prepared=prepared)
    assert a["schema"] == RESULT_SCHEMA_V2
    assert _dump(a) == _dump(b)


def test_run_spec_v2_accepts_study_data_source_and_does_not_reapply_input():
    n = _norm()
    ds = _ds()
    st = D.from_dataset(ds)
    assert _dump(run_spec_v2(n, ds, assay_kind="lipid")) == _dump(run_spec_v2(n, st, assay_kind="lipid"))
    # StudyData には spec の input を再適用しない（作った側が反映済みという前提）
    ex = _norm(input={"exclude_samples": ["C_B1_0"]})
    pre = prepare_spec_v2(ex, st)
    assert any("C_B1_0" in u.sample_ids for u in pre.units)
    assert not any("C_B1_0" in u.sample_ids for u in prepare_spec_v2(ex, ds).units)


# --- 割れのまとめ ---------------------------------------------------------------

def test_merge_feature_columns_takes_max_and_or_detected():
    rows = [(f"a{i}", "sample", "C", None, None, None) for i in range(4)]
    vals = np.array([[10.0, 20.0, 5.0, 1.0], [30.0, np.nan, 6.0, 2.0],
                     [0.0, -3.0, 7.0, 3.0], [8.0, 4.0, 8.0, 4.0]])
    det = np.array([[1, 0, 1, 0], [0, 0, 1, 1], [0, 0, 1, 0], [1, 1, 0, 0]], bool)
    st = _study(vals, rows, detected=det)
    st.feature_classes = {"f1": "PC", "f3": "TG"}
    st.annotated_features = frozenset({"f1"})
    out = D.merge_feature_columns(st, [["f0", "f1"], ["f2", "missing"]])
    assert out.feature_ids == ["f0", "f2", "f3"]
    np.testing.assert_array_equal(out.values[:, 0], [20.0, 30.0, np.nan, 8.0])
    np.testing.assert_array_equal(out.values[:, 1:], vals[:, 2:])
    np.testing.assert_array_equal(out.detected[:, 0], [True, False, False, True])
    assert out.feature_classes == {"f0": "PC", "f3": "TG"} and out.annotated_features == {"f0"}
    assert st.feature_ids == ["f0", "f1", "f2", "f3"] and st.values.shape == (4, 4)  # 元は不変
    assert _code(D.merge_feature_columns, st, [["f0", "f1"], ["f1", "f2"]]).code == "STUDY_DATA_INVALID"
    assert D.merge_feature_columns(st, []).feature_ids == st.feature_ids


# --- 修正ラウンド 1 -------------------------------------------------------------

def test_merge_and_exclude_give_different_matrix_ids():
    ds = _ds()
    merged = D.merge_feature_columns(D.from_dataset(ds), [["f0", "f1"]])
    excluded = D.from_dataset(ds, exclude_features=["f1"])
    assert merged.feature_ids == excluded.feature_ids
    assert _pre(merged).matrix_id != _pre(excluded).matrix_id
    again = D.merge_feature_columns(D.from_dataset(ds), [["f0", "f1"]])
    assert _pre(merged).matrix_id == _pre(again).matrix_id


def test_injection_columns_dropped_from_unit_attributes_ids_and_warning():
    rows = []
    for s, g in (("S0", "C"), ("S1", "C"), ("S2", "T"), ("S3", "T")):
        for r in range(2):
            rows.append((f"{s}_{r}", "sample", g, None, None, s, True,
                         {"run_day": f"d{r}", "plate": f"p{r}"}))
    vals = 1000.0 * (1 + 0.1 * np.random.default_rng(1).standard_normal((len(rows), 3)))
    prep = _pre(_study(vals, rows), injection_columns=["run_day"])
    assert len(prep.units) == 8  # plate がまだ割る
    assert all("run_day" not in u.attributes and "plate" in u.attributes for u in prep.units)
    assert "S0|C|plate=p0" in {u.unit_id for u in prep.units}
    assert not any("run_day" in u.unit_id for u in prep.units)
    w = next(w for w in prep.warnings if w["code"] == "units_split_by_attribute")
    assert w["attributes"] == ["plate"]
    # 注入の列を因子に使う設計は、最初の注入の値を黙って使わず欠けた値として外れる
    e = _code(compare, prep, factors=["group"], comparisons=BVA,
              covariates=[{"name": "run_day", "type": "categorical"}])
    assert e.code == "DESIGN_INVALID"  # 全単位で値が欠ける: 最初の注入の値を黙って使わない


def test_split_warning_reports_count_of_the_rest():
    rows = []
    for k in range(25):
        for r in range(2):
            rows.append((f"S{k}_{r}", "sample", "C" if k % 2 else "T", None, None, f"S{k}", True,
                         {"run_day": f"d{r}"}))
    vals = np.full((len(rows), 2), 100.0) + np.arange(len(rows))[:, None]
    w = next(w for w in _pre(_study(vals, rows)).warnings if w["code"] == "units_split_by_attribute")
    assert w["subjects_n_more"] == 5 and len(w["subjects"]) == 20


def test_pca_value_error_from_imputation_becomes_imputation_failed(monkeypatch):
    from metabolomix.analysis.stats import pca as pca_mod
    rng = np.random.default_rng(1)
    log2 = rng.normal(10, 1, (6, 5))
    log2[0, 0] = np.nan
    prep = _matrix(log2, [_unit(f"u{i}", "C" if i < 3 else "T") for i in range(6)])

    def boom(*a, **k):
        raise ValueError("boom")
    monkeypatch.setattr(pca_mod, "impute_for_pca", boom)
    e = _code(pca, prep, impute="half_min")
    assert e.code == "PCA_IMPUTATION_FAILED" and "half_min" in e.message


def test_merge_dedupes_ids_within_group():
    rows = [(f"a{i}", "sample", "C", None, None, None) for i in range(3)]
    st = _study(np.arange(1.0, 10.0).reshape(3, 3), rows)
    out = D.merge_feature_columns(st, [["f0", "f0"]])
    assert out.feature_ids == st.feature_ids
    np.testing.assert_array_equal(out.values, st.values)


def test_equivalent_inputs_canonicalize_to_same_hash():
    assert (P.normalize_preprocess_config({"injection_columns": [" b", "a", "b"]})["injection_columns"]
            == ["a", "b"])
    a = _norm(input={"group_overrides": {" C_B1_0 ": " T "}})
    b = _norm(input={"group_overrides": {"C_B1_0": "T"}})
    assert a == b and a["input"]["group_overrides"] == {"C_B1_0": "T"}
