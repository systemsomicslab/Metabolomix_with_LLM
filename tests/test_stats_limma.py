"""limma 型の moderated t / F と区間を R の limma と突き合わせる（C §6・§7）。

参照値 `stats_reference/limma_expected.json` は `make_limma_expected.R` が R の limma で一度だけ作った
もの（版は JSON の `_meta`）。このテストは JSON と入力の TSV だけを読み、R を要らない。
計画行列と対比の行列は JSON に書かれたもの（R の model.matrix の出力）をそのまま使う。
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest
from pytest import approx
from scipy.special import polygamma

from metabolomix.analysis.stats.limma import (
    contrasts_fit,
    ebayes,
    fit_f_dist,
    fit_f_dist_unequal_df1,
    lm_fit,
    squeeze_var,
    trigamma_inverse,
)

REF = Path(__file__).parent / "stats_reference"
EXPECTED = json.loads((REF / "limma_expected.json").read_text(encoding="utf-8"))


def _read_input() -> tuple[list[str], list[str], np.ndarray]:
    with (REF / "limma_input.tsv").open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    samples = rows[0][1:]
    features = [r[0] for r in rows[1:]]
    values = np.array([[np.nan if v == "NA" else float(v) for v in r[1:]] for r in rows[1:]])
    return features, samples, values


FEATURES, SAMPLES, VALUES = _read_input()


def _num(x) -> float:
    return math.inf if x == "Inf" else float(x)


def _case(design_id: str):
    exp = EXPECTED[design_id]
    fi = [FEATURES.index(f) for f in exp["features"]]
    si = [SAMPLES.index(s) for s in exp["samples"]]
    y = VALUES[np.ix_(fi, si)]
    x = np.array(exp["design"], dtype=float)
    c = np.array(exp["contrasts"], dtype=float)
    return exp, y, x, c


DESIGNS = ["cell_means", "cell_means_complete", "batch_covariate", "paired_block"]
TOL = dict(rel=1e-6, abs=1e-10)


@pytest.mark.parametrize("design_id", DESIGNS)
def test_moderated_t_matches_limma(design_id):
    exp, y, x, c = _case(design_id)
    fit = lm_fit(y, x, exp["column_names"])
    assert fit.df_residual == approx(np.array(exp["df_residual"], dtype=float))
    assert fit.sigma == approx(np.array(exp["sigma"]), **TOL)

    cfit = contrasts_fit(fit, c, exp["contrast_names"])
    assert cfit.stdev_unscaled == approx(np.array(exp["stdev_unscaled"]), **TOL)

    mod = ebayes(cfit)
    assert mod.contrast_names == exp["contrast_names"]
    assert mod.s2_prior == approx(_num(exp["s2_prior"]), rel=1e-6)
    assert mod.df_prior == approx(_num(exp["df_prior"]), rel=1e-6)
    assert mod.s2_post == approx(np.array(exp["s2_post"]), **TOL)
    assert mod.df_total == approx(np.array(exp["df_total"]), **TOL)
    for k, name in enumerate(exp["contrast_names"]):
        ref = exp["contrast"][name]
        assert mod.coefficients[:, k] == approx(np.array(ref["logFC"]), **TOL)
        assert mod.t[:, k] == approx(np.array(ref["t"]), **TOL)
        assert mod.p[:, k] == approx(np.array(ref["P.Value"]), **TOL)
        assert mod.ci_low[:, k] == approx(np.array(ref["CI.L"]), **TOL)
        assert mod.ci_high[:, k] == approx(np.array(ref["CI.R"]), **TOL)
    assert mod.f == approx(np.array(exp["F"]), rel=1e-6)
    assert mod.f_p == approx(np.array(exp["F_p_value"]), rel=1e-6, abs=1e-300)


def test_reference_exercises_both_prior_estimators():
    """欠測で自由度が揃わない設計と揃う設計の両方が参照に入っている（推定の経路が 2 つあるため）。"""
    assert len(set(EXPECTED["cell_means"]["df_residual"])) > 1
    assert len(set(EXPECTED["cell_means_complete"]["df_residual"])) == 1


def test_fit_f_dist_degenerate_inputs():
    assert fit_f_dist(np.array([0.3]), np.array([4.0])) == (0.3, 0.0)
    s2p, d0 = fit_f_dist(np.zeros(5), np.full(5, 4.0))
    assert math.isfinite(s2p) and s2p > 0
    s2p, d0 = fit_f_dist(np.full(5, 0.2), np.full(5, 4.0))
    assert s2p == approx(0.2) and d0 == math.inf

    # 全特徴が同じ分散（df_prior=inf）でも p は有限
    rng = np.random.default_rng(7)
    noise = rng.normal(0, 1, (1, 6))
    y = np.vstack([noise + i for i in range(5)])       # 残差は全特徴で同一
    x = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]])
    mod = ebayes(lm_fit(y, x, ["Intercept", "g"]))
    assert mod.df_prior == math.inf
    assert np.all(np.isfinite(mod.p))
    assert np.all(np.isfinite(mod.f_p))
    assert np.all(mod.df_total == approx(4.0 * 5))     # Σ df_residual で頭打ち


def test_squeeze_var_few_features_has_zero_prior_df():
    """limma squeezeVar と同じく、特徴が 3 未満なら事前分布で縮めない（df_prior=0）。"""
    s2_prior, df_prior, s2_post = squeeze_var(np.array([0.1, 0.4]), np.array([3.0, 3.0]))
    assert df_prior == 0.0
    assert s2_post == approx([0.1, 0.4])
    assert math.isfinite(s2_prior)


def test_fit_f_dist_unequal_df1_needs_two_informative():
    s2p, d0 = fit_f_dist_unequal_df1(np.array([0.0, 0.2, 0.0]), np.array([3.0, 4.0, 5.0]))
    assert math.isnan(s2p) and math.isnan(d0)


def test_trigamma_inverse_roundtrip():
    x = np.array([0.05, 0.5, 2.0, 50.0])
    assert polygamma(1, trigamma_inverse(x)) == approx(x, rel=1e-8)
