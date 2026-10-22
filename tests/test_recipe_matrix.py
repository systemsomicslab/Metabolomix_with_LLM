"""行列の recipe v2（spec 2026-10-22 §4 の 4）。値はすべて手で計算できる小さな行列で確かめる。"""
from __future__ import annotations

import math

import numpy as np
import pytest

from metabolomix.analysis import recipe_matrix as R
from metabolomix.core.atomic_io import DomainError

NAN = math.nan


def _build(values, recipe, *, roles=None, detected=None, order=None, samples=None):
    values = np.asarray(values, dtype=float)
    samples = samples or [f"S{i}" for i in range(values.shape[0])]
    features = [f"F{j}" for j in range(values.shape[1])]
    roles = roles or {s: "sample" for s in samples}
    return R.build_matrix(values, samples, features, roles, R.normalize_recipe(recipe),
                          detected=detected, injection_order=order)


# ---------- 正規化 ----------

def test_defaults_are_filled_and_echoed():
    out = R.normalize_recipe({"steps": ["impute"]})
    assert out == {"schema": "matrix-recipe.v2", "steps": ["impute"],
                   "impute": {"method": "half_min", "missing": "zero", "k": 10, "scale": "log2",
                              "rowmax": 0.5, "colmax": 0.8}}


def test_empty_steps_is_no_preprocessing():
    assert R.normalize_recipe({}) == {"schema": "matrix-recipe.v2", "steps": []}


@pytest.mark.parametrize("recipe", [
    {"steps": ["impute"], "bogus": 1},                       # 未知のキー
    {"steps": ["impute"], "normalize": {"method": "pqn"}},    # steps に無い段の設定
    {"steps": ["impute", "impute"]},                          # 重複
    {"steps": ["sort"]},                                      # 未知の段
    {"steps": ["impute"], "impute": {"method": "mice"}},      # 未知の方法
    {"steps": ["impute"], "impute": {"k": 0, "method": "knn"}},
    {"steps": ["detection_filter"], "detection_filter": {"min_rate": {"plasma": 0.5}}},
])
def test_invalid_recipe_rejected(recipe):
    with pytest.raises(DomainError) as info:
        R.normalize_recipe(recipe)
    assert info.value.code == "MATRIX_RECIPE_INVALID"


# ---------- 順序 ----------

def test_steps_run_in_given_order():
    values = [[1, 2], [3, 0], [5, 6]]           # S1 の F1 が欠測（0）
    a = _build(values, {"steps": ["impute", "normalize"], "impute": {"scale": "none"},
                        "normalize": {"method": "tic", "reference": "sample"}})
    b = _build(values, {"steps": ["normalize", "impute"], "impute": {"scale": "none"},
                        "normalize": {"method": "tic", "reference": "sample"}})
    assert not np.allclose(a["values"], b["values"])
    assert [s["step"] for s in a["steps"]] == ["impute", "normalize"]


def test_inputs_not_mutated():
    values = np.array([[1.0, 0.0], [3.0, 4.0]])
    before = values.copy()
    R.build_matrix(values, ["A", "B"], ["F0", "F1"], {"A": "sample", "B": "sample"},
                   R.normalize_recipe({"steps": ["impute"]}))
    assert np.array_equal(values, before)


# ---------- 検出率 ----------

def test_detection_filter_by_role_positive_vs_detected():
    samples = ["S1", "S2", "S3", "Q1", "Q2"]
    roles = {"S1": "sample", "S2": "sample", "S3": "sample", "Q1": "qc", "Q2": "qc"}
    values = np.array([[5, 1], [5, 1], [0, 1], [5, 1], [5, 1]], dtype=float)
    detected = np.array([[1, 0], [1, 0], [0, 0], [1, 1], [1, 1]], dtype=bool)  # F1 は gap-fill
    recipe = {"steps": ["detection_filter"],
              "detection_filter": {"basis": "positive", "min_rate": {"qc": 0.5, "sample": 0.5}}}
    pos = _build(values, recipe, roles=roles, samples=samples, detected=detected)
    assert pos["feature_ids"] == ["F0", "F1"]
    recipe["detection_filter"]["basis"] = "detected"
    det = _build(values, recipe, roles=roles, samples=samples, detected=detected)
    assert det["feature_ids"] == ["F0"]                     # F1 は試料で実検出 0/3
    assert det["steps"][0]["removed"] == ["F1"]


def test_detected_basis_requires_detection_state():
    with pytest.raises(DomainError) as info:
        _build([[1.0]], {"steps": ["detection_filter"],
                         "detection_filter": {"basis": "detected", "min_rate": {"sample": 0.5}}})
    assert info.value.code == "DETECTION_UNAVAILABLE"


# ---------- 補完 ----------

def test_half_min_uses_positive_minimum_and_keeps_all_missing_feature():
    out = _build([[4, 0], [0, 0], [8, 0]], {"steps": ["impute"], "impute": {"scale": "none"}})
    assert out["values"][1, 0] == 2.0
    assert math.isnan(out["values"][0, 1])                  # 正の値が無い特徴は欠測のまま


def _knn(features_by_samples, **impute):
    values = np.asarray(features_by_samples, dtype=float).T   # 試料 × 特徴へ
    recipe = {"steps": ["impute"],
              "impute": {"method": "knn", "missing": "nan", "scale": "none", **impute}}
    return _build(values, recipe)["values"].T


KNN_DATA = [
    [1, 2, 3, NAN],          # F0: 欠測 1 つ
    [1, 2, 3, 4],            # F1: 距離 0
    [2, 3, 4, 5],            # F2: 距離 1
    [10, 20, 30, 40],        # F3: 距離 378
    [1.1, 2.1, 2.9, 3.9],    # F4: 距離 0.01
]


def test_knn_matches_hand_computed_impute_knn():
    out = _knn(KNN_DATA, k=2)
    assert out[0, 3] == pytest.approx((4 + 3.9) / 2)       # 近傍は F1・F4
    assert np.allclose(out[1:], np.asarray(KNN_DATA[1:], dtype=float))


def test_knn_skips_neighbours_missing_in_that_column():
    data = [row[:] for row in KNN_DATA]
    data[4][3] = NAN                                        # F4 も列 3 が欠測
    out = _knn(data, k=2)
    # F0 の近傍 F1(0)・F4(0.01) のうち F4 は列 3 を持たないので F1 だけの平均
    assert out[0, 3] == pytest.approx(4.0)


def test_knn_rowmax_uses_column_mean():
    data = [[1, NAN, NAN, NAN]] + [row for row in KNN_DATA[1:]]
    out = _knn(data, k=2, rowmax=0.5)
    col_means = np.nanmean(np.asarray(KNN_DATA[1:], dtype=float), axis=0)
    assert np.allclose(out[0, 1:], col_means[1:])


def test_knn_colmax_rejects():
    data = [[1, NAN], [2, NAN], [3, 4]]
    with pytest.raises(DomainError) as info:
        _knn(data, k=1, colmax=0.5)
    assert info.value.code == "RECIPE_PRECONDITION"


def test_knn_log2_scale_back_transforms():
    data = [[2, 4, 8, NAN], [2, 4, 8, 16], [1000, 1, 1, 1]]
    out = _knn(data, k=1, scale="log2")
    assert out[0, 3] == pytest.approx(16.0)


# ---------- 正規化 ----------

def test_pqn_reference_qc():
    samples = ["S1", "S2", "Q1", "Q2"]
    roles = {"S1": "sample", "S2": "sample", "Q1": "qc", "Q2": "qc"}
    values = np.array([[2, 4, 6], [1, 1, 1], [1, 2, 3], [1, 2, 3]], dtype=float)
    out = _build(values, {"steps": ["normalize"], "normalize": {"method": "pqn", "reference": "qc"}},
                 roles=roles, samples=samples)
    assert np.allclose(out["values"][0], [1, 2, 3])          # 商の中央値 2 で割る
    assert out["steps"][0]["factors"]["S1"] == pytest.approx(2.0)


# ---------- ドリフト補正 ----------

def test_loess_matches_statsmodels_direct_call():
    from statsmodels.nonparametric.smoothers_lowess import lowess

    n = 16
    samples = [f"I{i}" for i in range(n)]
    roles = {s: ("qc" if i % 2 == 0 else "sample") for i, s in enumerate(samples)}
    order = {s: float(i + 1) for i, s in enumerate(samples)}
    column = np.array([100 + 3 * i + (5 if i % 3 == 0 else 0) for i in range(n)], dtype=float)
    out = _build(column[:, None], {"steps": ["drift_correct"],
                                   "drift_correct": {"span": 0.75, "iterations": 0, "min_qc": 4}},
                 roles=roles, samples=samples, order=order)
    qc = [i for i in range(n) if i % 2 == 0]
    fit = lowess(column[qc], np.array(qc, dtype=float) + 1, frac=0.75, it=0, delta=0.0,
                 xvals=np.arange(1, n + 1, dtype=float))
    assert np.all(np.isfinite(fit))
    expected = column / fit * np.median(column[qc])
    assert np.allclose(out["values"][:, 0], expected)
    assert out["steps"][0]["uncorrected_features"] == 0


def test_loess_fit_that_is_not_finite_leaves_feature_uncorrected():
    """statsmodels は QC が少ないと当てはめが NaN になる点を返す。そのときは補正しない。"""
    samples = [f"I{i}" for i in range(8)]
    roles = {s: ("qc" if i % 2 == 0 else "sample") for i, s in enumerate(samples)}
    order = {s: float(i + 1) for i, s in enumerate(samples)}
    column = np.array([100, 50, 110, 55, 125, 60, 130, 70], dtype=float)
    out = _build(column[:, None], {"steps": ["drift_correct"],
                                   "drift_correct": {"span": 0.75, "iterations": 0, "min_qc": 4}},
                 roles=roles, samples=samples, order=order)
    assert np.array_equal(out["values"][:, 0], column)
    assert out["steps"][0]["uncorrected_features"] == 1


def test_drift_requires_order_for_every_sample():
    samples = ["A", "B", "C", "D", "E"]
    roles = {s: "qc" for s in samples}
    with pytest.raises(DomainError) as info:
        _build(np.ones((5, 1)), {"steps": ["drift_correct"]}, roles=roles, samples=samples,
               order={"A": 1, "B": 2, "C": 3, "D": 4})
    assert info.value.code == "INJECTION_ORDER_MISSING"
    assert info.value.details["samples"] == ["E"]


def test_drift_leaves_feature_with_too_few_qc_uncorrected():
    samples = ["Q1", "Q2", "S1"]
    roles = {"Q1": "qc", "Q2": "qc", "S1": "sample"}
    out = _build([[1.0], [2.0], [3.0]], {"steps": ["drift_correct"]}, roles=roles,
                 samples=samples, order={"Q1": 1, "Q2": 2, "S1": 3})
    assert np.allclose(out["values"][:, 0], [1, 2, 3])
    assert out["steps"][0]["uncorrected_features"] == 1


# ---------- QC の RSD ----------

def test_qc_rsd_filter_drops_high_and_nonfinite_rsd():
    samples = ["Q1", "Q2", "Q3", "S1"]
    roles = {"Q1": "qc", "Q2": "qc", "Q3": "qc", "S1": "sample"}
    values = np.array([[10, 1, NAN], [10, 10, NAN], [10, 100, NAN], [5, 5, 5]], dtype=float)
    out = _build(values, {"steps": ["qc_rsd_filter"], "qc_rsd_filter": {"max_rsd": 0.3}},
                 roles=roles, samples=samples)
    assert out["feature_ids"] == ["F0"]
    step = out["steps"][0]
    assert step["removed_high_rsd"] == ["F1"] and step["removed_no_rsd"] == ["F2"]


def test_filtered_features_do_not_feed_later_steps():
    """篩で落とした特徴は直後に外れ、後の PQN の商に混ざらない。"""
    samples = ["Q1", "Q2", "S1"]
    roles = {"Q1": "qc", "Q2": "qc", "S1": "sample"}
    values = np.array([[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 100, 100]], dtype=float)
    detected = np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 0, 0]], dtype=bool)
    out = _build(values, {"steps": ["detection_filter", "normalize"],
                          "detection_filter": {"basis": "detected", "min_rate": {"sample": 1.0}},
                          "normalize": {"method": "pqn", "reference": "qc"}},
                 roles=roles, samples=samples, detected=detected)
    assert out["feature_ids"] == ["F0", "F1"]
    assert out["steps"][1]["factors"]["S1"] == pytest.approx(2.0)
