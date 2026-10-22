"""再解析の spec が呼ぶ検定・効果量・補正（spec 2026-10-22 §4 の 3）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。v1 の
`differential.two_group_test`（log2(x+1) の Welch）と v2 の `statistics_v2`（Welch・ANOVA）は
変えない——ここは論文の方法を写すための別の入口で、定義は呼び出しごとに名前で選ぶ:

- `student`: 等分散の t（`scipy.stats.ttest_ind(equal_var=True)`）、両側。
- `welch`: 不等分散の t、両側。
- `mannwhitney`: `scipy.stats.mannwhitneyu(alternative="two-sided")`。
- `ols`: `transform(y) ~ 1 + group(test=1) + 共変量` の最小二乗。群の係数の t の両側 p。

`transform="log2"` は有限の正値だけに当て、ほかは欠測にする（clip も pseudocount も
入れない。statistics_v2 と同じ規則）。
"""
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy import stats

__all__ = ["FC_DEFINITIONS", "KINDS", "bh", "fold_change", "ols_group", "transform_values",
           "two_group"]

KINDS = ("student", "welch", "mannwhitney", "ols")
FC_DEFINITIONS = ("arithmetic_mean_ratio", "geometric_mean_ratio")
TRANSFORMS = ("log2", "none")
_MIN_N = 2


def transform_values(values, transform: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if transform == "none":
        return np.where(np.isfinite(array), array, np.nan)
    out = np.full(array.shape, np.nan)
    usable = np.isfinite(array) & (array > 0)
    out[usable] = np.log2(array[usable])
    return out


def _result(status="not_testable", reason=None, statistic=None, p=None, **extra) -> dict:
    return {"status": status, "reason": reason, "statistic": statistic, "p": p, **extra}


def two_group(ref_raw, test_raw, *, kind: str, transform: str) -> dict:
    """2 群の検定（student / welch / mannwhitney）。`n_ref`/`n_test` は使った値の数。"""
    a = transform_values(ref_raw, transform)
    b = transform_values(test_raw, transform)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    counts = {"n_ref": int(a.size), "n_test": int(b.size)}
    if a.size < _MIN_N or b.size < _MIN_N:
        return _result(reason="insufficient_finite_values", **counts)
    if kind in ("student", "welch") and np.var(a) == 0 and np.var(b) == 0:
        return _result(reason="zero_variance", **counts)
    # ほぼ同じ値の群で scipy が出す桁落ちの警告は、特徴数ぶん並ぶだけなので抑える
    # （p はそのまま返す。分散 0 の群は上で not_testable にしてある）。
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if kind == "mannwhitney":
            statistic, p = stats.mannwhitneyu(b, a, alternative="two-sided")
        else:
            statistic, p = stats.ttest_ind(b, a, equal_var=(kind == "student"))
    if not math.isfinite(p):
        return _result(reason="test_statistic_not_finite", **counts)
    return _result(status="tested", statistic=float(statistic), p=float(p), **counts)


def ols_group(y_raw, group, covariates, *, transform: str) -> dict:
    """`transform(y) ~ 1 + group + covariates` の群の係数。

    `group` は 0（基準）/1（比較）、`covariates` は（n × p）の数値行列（符号化済み）。
    y が欠測の行は除く。自由度が 0 以下・計画行列が階数不足なら not_testable。
    """
    y = transform_values(y_raw, transform)
    group = np.asarray(group, dtype=float)
    cov = np.asarray(covariates, dtype=float).reshape(len(group), -1)
    design = np.column_stack([np.ones(len(group)), group, cov])
    ok = np.isfinite(y) & np.all(np.isfinite(design), axis=1)
    y, design = y[ok], design[ok]
    n, p = design.shape
    counts = {"n_ref": int((design[:, 1] == 0).sum()), "n_test": int((design[:, 1] == 1).sum()),
              "beta_log2": None}
    if counts["n_ref"] < _MIN_N or counts["n_test"] < _MIN_N:
        return _result(reason="insufficient_finite_values", **counts)
    dof = n - p
    if dof <= 0:
        return _result(reason="no_residual_degrees_of_freedom", **counts)
    if np.linalg.matrix_rank(design) < p:
        return _result(reason="rank_deficient", **counts)
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    s2 = float(residual @ residual) / dof
    se = math.sqrt(s2 * np.linalg.inv(design.T @ design)[1, 1])
    if se == 0:
        return _result(reason="zero_variance", **counts)
    t = float(beta[1]) / se
    counts["beta_log2"] = float(beta[1])
    return _result(status="tested", statistic=t, p=float(2 * stats.t.sf(abs(t), dof)), **counts)


def fold_change(ref_raw, test_raw, definition: str) -> float | None:
    """test / ref の倍率。平均が 0 以下・値が無ければ None。"""
    a = np.asarray(ref_raw, dtype=float)
    b = np.asarray(test_raw, dtype=float)
    if definition == "arithmetic_mean_ratio":
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if not a.size or not b.size or a.mean() <= 0 or b.mean() <= 0:
            return None
        return float(b.mean() / a.mean())
    la, lb = transform_values(a, "log2"), transform_values(b, "log2")
    la, lb = la[np.isfinite(la)], lb[np.isfinite(lb)]
    if not la.size or not lb.size:
        return None
    return float(2.0 ** (lb.mean() - la.mean()))


def bh(pvalues: list) -> list:
    """Benjamini–Hochberg の q。None・非有限は母集団に入れず None を返す。"""
    index = [i for i, p in enumerate(pvalues) if p is not None and math.isfinite(p)]
    out: list = [None] * len(pvalues)
    if not index:
        return out
    p = np.array([pvalues[i] for i in index], dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    q = np.empty(len(p))
    q[order] = adjusted
    for i, value in zip(index, q):
        out[i] = float(value)
    return out
