"""再解析の spec が呼ぶ検定・FC・BH（spec 2026-10-22 §4 の 3）。"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from metabolomix.analysis import spec_statistics as S

REF = [10.0, 12.0, 11.0, 13.0]
TEST = [20.0, 22.0, 19.0, 25.0]


def test_student_on_log2_matches_scipy():
    out = S.two_group(REF, TEST, kind="student", transform="log2")
    expected = stats.ttest_ind(np.log2(TEST), np.log2(REF), equal_var=True)
    assert out["status"] == "tested"
    assert out["p"] == pytest.approx(expected.pvalue)
    assert out["statistic"] == pytest.approx(expected.statistic)
    assert (out["n_ref"], out["n_test"]) == (4, 4)


def test_welch_differs_from_student_with_unequal_variance():
    a, b = [1.0, 1.1, 0.9, 1.0], [2.0, 8.0, 0.5, 5.0]
    student = S.two_group(a, b, kind="student", transform="none")["p"]
    welch = S.two_group(a, b, kind="welch", transform="none")["p"]
    assert welch == pytest.approx(stats.ttest_ind(b, a, equal_var=False).pvalue)
    assert student != pytest.approx(welch)


def test_mannwhitney_two_sided():
    out = S.two_group(REF, TEST, kind="mannwhitney", transform="log2")
    assert out["p"] == pytest.approx(stats.mannwhitneyu(TEST, REF, alternative="two-sided").pvalue)


def test_log2_drops_non_positive_values_and_needs_two():
    out = S.two_group([0.0, -1.0, 5.0], TEST, kind="student", transform="log2")
    assert out["status"] == "not_testable" and out["reason"] == "insufficient_finite_values"
    assert out["n_ref"] == 1


def test_zero_variance_is_not_testable():
    out = S.two_group([5.0, 5.0], [5.0, 5.0], kind="welch", transform="none")
    assert out["status"] == "not_testable" and out["reason"] == "zero_variance"


def test_ols_group_p_matches_closed_form():
    y = np.array(REF + TEST)
    group = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)
    age = np.array([30, 40, 35, 50, 45, 33, 60, 41], dtype=float)
    out = S.ols_group(y, group, age[:, None], transform="log2")
    X = np.column_stack([np.ones(8), group, age])
    ly = np.log2(y)
    beta = np.linalg.solve(X.T @ X, X.T @ ly)
    resid = ly - X @ beta
    se = math.sqrt(resid @ resid / 5 * np.linalg.inv(X.T @ X)[1, 1])
    assert out["beta_log2"] == pytest.approx(beta[1])
    assert out["p"] == pytest.approx(2 * stats.t.sf(abs(beta[1] / se), 5))


def test_ols_without_covariates_equals_student():
    y = np.array(REF + TEST)
    group = np.array([0] * 4 + [1] * 4, dtype=float)
    ols = S.ols_group(y, group, np.zeros((8, 0)), transform="log2")
    student = S.two_group(REF, TEST, kind="student", transform="log2")
    assert ols["p"] == pytest.approx(student["p"])


def test_ols_rank_deficient_is_not_testable():
    y = np.array(REF + TEST)
    group = np.array([0] * 4 + [1] * 4, dtype=float)
    out = S.ols_group(y, group, group[:, None], transform="log2")   # 共変量 = 群
    assert out["status"] == "not_testable" and out["reason"] == "rank_deficient"


def test_fold_change_definitions():
    assert S.fold_change([1, 3], [4, 4], "arithmetic_mean_ratio") == pytest.approx(2.0)
    assert S.fold_change([1, 4], [4, 4], "geometric_mean_ratio") == pytest.approx(2.0)
    assert S.fold_change([0, 0], [1, 1], "arithmetic_mean_ratio") is None


def test_bh_skips_missing_p():
    q = S.bh([0.01, None, 0.04, math.nan, 0.03])
    assert q[1] is None and q[3] is None
    assert q[0] == pytest.approx(0.03) and q[2] == pytest.approx(0.04)
    assert q[4] == pytest.approx(0.04)
