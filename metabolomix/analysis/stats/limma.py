"""limma 型の moderated t / F と区間（spec §6・§7）。

R の limma 3.68.5 の `lmFit`（`lm.series`）・`contrasts.fit`・`eBayes`（`.ebayes`・`squeezeVar`・
`fitFDist`・`fitFDistUnequalDF1`・`trigammaInverse`・`classifyTestsF`）と `topTable(confint=)` を
移植したもの。`trend=FALSE`・`robust=FALSE` だけを扱う。R の参照値との突き合わせは
`tests/test_stats_limma.py`（参照値の作り方は `tests/stats_reference/make_limma_expected.R`）。

行は特徴、列は試料（`y` は G×n）。欠測は NaN（非有限）で、特徴ごとにその試料を除いて当てはめる
（limma `lm.series` の欠測の扱い）。

limma と意図して違えている点（どれも limma なら止まる／意味の無い値を返す縁の場合だけ）:
- 計画行列が列フルランクでなければ `ValueError`（limma は推定できない係数を落として続ける）。
- 事前分布の自由度が推定できない（情報のある分散が 2 つ未満）とき、limma の `squeezeVar` は止まるが、
  ここでは `df_prior=0`（縮めない）として続ける。
- 特徴が 3 未満のとき limma の `var.prior` は特徴ごとの分散そのもの（ベクトル）になる。ここでは
  `s2_prior` を有限な分散の中央値で返す（`df_prior=0` なので事後分散には効かない）。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import special, stats

__all__ = [
    "LinearFit",
    "Moderated",
    "lm_fit",
    "contrasts_fit",
    "trigamma_inverse",
    "fit_f_dist",
    "fit_f_dist_unequal_df1",
    "squeeze_var",
    "ebayes",
]

# R の `optimize` の既定の許容（.Machine$double.eps^0.25）
_OPTIMIZE_TOL = float(np.finfo(float).eps) ** 0.25
# LINPACK dqrdc2（R の lm.fit）が列を推定不能とみなす相対の閾値
_QR_TOL = 1e-7


@dataclass
class LinearFit:
    """特徴ごとの線形モデルの当てはめ（limma の MArrayLM の必要な部分）。

    `coefficients`・`stdev_unscaled` は G×p（推定できない係数は NaN）、`sigma`・`df_residual` は (G,)。
    `cov_coefficients` は欠測なしの計画行列から求めた p×p の (XᵀX)⁻¹（対比を当てたあとは CᵀVC）。
    """

    coefficients: np.ndarray
    stdev_unscaled: np.ndarray
    sigma: np.ndarray
    df_residual: np.ndarray
    cov_coefficients: np.ndarray
    column_names: list[str]


@dataclass
class Moderated:
    """`ebayes` の結果。G×k の行列は対比（または係数）ごとの列。

    `f`・`f_p` は limma `classifyTestsF(fstat.only=TRUE)` の F と上側確率。F の分母の自由度
    `f_df2` は limma と同じく `df_prior + df_residual`（`df_total` と違い Σdf_residual で頭打ちにしない）。
    """

    t: np.ndarray
    p: np.ndarray
    coefficients: np.ndarray
    ci_low: np.ndarray
    ci_high: np.ndarray
    df_total: np.ndarray
    s2_post: np.ndarray
    s2_prior: float
    df_prior: float
    f: np.ndarray
    f_p: np.ndarray
    contrast_names: list[str]
    f_df1: int = 0
    f_df2: np.ndarray | None = None


# ---------------------------------------------------------------------------
# 当てはめ（lmFit / lm.series）


def _is_full_rank(x: np.ndarray) -> bool:
    """limma `is.fullrank`: XᵀX の固有値の比で判定する。"""
    e = np.linalg.eigvalsh(x.T @ x)[::-1]
    return bool(e[0] > 0 and abs(e[-1] / e[0]) > 1e-13)


def _estimable_columns(x: np.ndarray) -> list[int]:
    """R の lm.fit（LINPACK dqrdc2）と同じ規則で推定できる列を元の順に選ぶ。

    先に選んだ列で説明したあとの残りのノルムが、元のノルムの `_QR_TOL` 倍未満の列は落とす。
    """
    selected: list[int] = []
    basis = np.zeros((x.shape[0], 0))
    for j in range(x.shape[1]):
        col = x[:, j]
        norm0 = float(np.linalg.norm(col))
        if norm0 == 0.0:
            continue
        resid = col - basis @ (basis.T @ col)
        norm_r = float(np.linalg.norm(resid))
        if norm_r < _QR_TOL * norm0:
            continue
        basis = np.column_stack([basis, resid / norm_r])
        selected.append(j)
    return selected


def lm_fit(y: np.ndarray, design: np.ndarray, column_names: list[str]) -> LinearFit:
    """特徴ごとに最小二乗で当てはめる（limma `lmFit` → `lm.series`、重みなし）。

    特徴ごとに非有限の試料を除いて当てはめ、除いた後の計画で推定できない係数は NaN にする。
    残差の自由度が 0 の特徴の `sigma` は NaN。`cov_coefficients` は欠測なしの計画行列から。
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(design, dtype=float)
    if y.ndim != 2 or x.ndim != 2:
        raise ValueError("y は G×n、design は n×p の 2 次元配列であること")
    n_genes, n_arrays = y.shape
    if x.shape[0] != n_arrays:
        raise ValueError("design の行数が y の列（試料）数と一致しない")
    p = x.shape[1]
    if len(column_names) != p:
        raise ValueError("column_names の数が design の列数と一致しない")
    if not np.all(np.isfinite(x)):
        raise ValueError("design に非有限の値がある")
    if not _is_full_rank(x):
        raise ValueError("design が列フルランクでない（推定できない係数がある）")

    _, r_full = np.linalg.qr(x)
    r_inv = np.linalg.solve(r_full, np.eye(p))
    cov = r_inv @ r_inv.T

    beta = np.full((n_genes, p), np.nan)
    stdev = np.full((n_genes, p), np.nan)
    sigma = np.full(n_genes, np.nan)
    df_res = np.zeros(n_genes)

    finite = np.isfinite(y)
    patterns: dict[bytes, list[int]] = {}
    for i in range(n_genes):
        patterns.setdefault(finite[i].tobytes(), []).append(i)
    for key, rows in patterns.items():
        obs = np.frombuffer(key, dtype=bool)
        if not obs.any():
            continue
        xo = x[obs]
        cols = _estimable_columns(xo)
        idx = np.asarray(rows)
        if not cols:
            df_res[idx] = obs.sum()
            continue
        xs = xo[:, cols]
        q, r = np.linalg.qr(xs)
        yo = y[np.ix_(idx, np.flatnonzero(obs))].T          # nobs × m
        effects = q.T @ yo
        b = np.linalg.solve(r, effects)                     # rank × m
        ri = np.linalg.solve(r, np.eye(len(cols)))
        su = np.sqrt(np.sum(ri * ri, axis=1))               # diag((XᵀX)⁻¹) の平方根
        df = int(obs.sum()) - len(cols)
        beta[np.ix_(idx, cols)] = b.T
        stdev[np.ix_(idx, cols)] = su
        df_res[idx] = df
        if df > 0:
            resid = yo - q @ effects
            sigma[idx] = np.sqrt(np.sum(resid * resid, axis=0) / df)
    return LinearFit(
        coefficients=beta,
        stdev_unscaled=stdev,
        sigma=sigma,
        df_residual=df_res,
        cov_coefficients=cov,
        column_names=list(column_names),
    )


# ---------------------------------------------------------------------------
# 対比（contrasts.fit）


def _cov2cor(v: np.ndarray) -> np.ndarray:
    d = np.sqrt(np.diag(v))
    return v / np.outer(d, d)


def contrasts_fit(fit: LinearFit, contrasts: np.ndarray, names: list[str]) -> LinearFit:
    """係数に対比を当てる（limma `contrasts.fit`）。

    特徴ごとの `stdev_unscaled` は、limma と同じく欠測なしの計画から作った係数の相関
    （`cov_coefficients` の相関行列の Cholesky 因子 R）を使って
    `√colSums((R·(stdev_unscaled[i]·C))²)` で求める（係数が直交なら `√(stdev_unscaled²·C²)`）。
    特徴ごとの欠測を除いた計画から計算し直すことはしない。
    """
    c = np.asarray(contrasts, dtype=float)
    if c.ndim == 1:
        c = c[:, None]
    coef = np.array(fit.coefficients, dtype=float)
    stdev = np.array(fit.stdev_unscaled, dtype=float)
    cov = np.array(fit.cov_coefficients, dtype=float)
    ncoef = coef.shape[1]
    if c.shape[0] != ncoef:
        raise ValueError("対比の行列の行数が係数の数と一致しない")
    if not np.all(np.isfinite(c)):
        raise ValueError("対比の行列に非有限の値がある")
    if len(names) != c.shape[1]:
        raise ValueError("names の数が対比の数と一致しない")
    cor = _cov2cor(cov)

    keep = np.abs(c).sum(axis=1) != 0
    if not keep.all():
        c = c[keep]
        coef = coef[:, keep]
        stdev = stdev[:, keep]
        cov = cov[np.ix_(keep, keep)]
        cor = cor[np.ix_(keep, keep)]

    na_coef = np.isnan(coef)
    has_na = bool(na_coef.any())
    if has_na:
        coef[na_coef] = 0.0
        stdev[na_coef] = 1e30

    new_coef = coef @ c
    if cor.size < 2:
        orthog = True
    else:
        orthog = bool(np.all(np.abs(cor[np.tril_indices_from(cor, k=-1)]) < 1e-14))
    new_cov = c.T @ cov @ c
    if orthog:
        new_stdev = np.sqrt((stdev ** 2) @ (c ** 2))
    else:
        r_upper = np.linalg.cholesky(cor).T                  # R の chol（上三角、RᵀR = cor）
        ruc = np.einsum("ab,gb,bk->gak", r_upper, stdev, c)  # 特徴ごとに R %*% (su_i * C)
        new_stdev = np.sqrt(np.sum(ruc ** 2, axis=1))
    if has_na:
        big = new_stdev > 1e20
        new_coef[big] = np.nan
        new_stdev[big] = np.nan
    return LinearFit(
        coefficients=new_coef,
        stdev_unscaled=new_stdev,
        sigma=np.array(fit.sigma, dtype=float),
        df_residual=np.array(fit.df_residual, dtype=float),
        cov_coefficients=new_cov,
        column_names=list(names),
    )


# ---------------------------------------------------------------------------
# 事前分布の推定（squeezeVar / fitFDist / fitFDistUnequalDF1）


def _logmdigamma(x: np.ndarray) -> np.ndarray:
    """log(x) − ψ(x)（statmod `logmdigamma` の移植。x ≤ 0 は NaN）。"""
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x) & (x > 0)
    z = x[ok]
    acc = np.zeros_like(z)
    small = z < 5
    # 小さい z は 5 ずらして漸近展開へ（log(z/(z+5)) + f(z+5) + Σ 1/(z+k)）
    zz = z.copy()
    acc[small] = (np.log(z[small] / (z[small] + 5)) + 1 / z[small] + 1 / (z[small] + 1)
                  + 1 / (z[small] + 2) + 1 / (z[small] + 3) + 1 / (z[small] + 4))
    zz[small] = z[small] + 5
    w = 1 / zz ** 2
    tail = w * (-1 / 12 + w * (1 / 120 + w * (-1 / 252 + w * (1 / 240 + w * (
        -1 / 132 + w * (691 / 32760 + w * (-1 / 12 + (3617 * w) / 8160)))))))
    out[ok] = acc + 1 / (2 * zz) - tail
    return out


def trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """ψ'(y) = x を解く y（limma `trigammaInverse`。Newton 法）。"""
    x = np.asarray(x, dtype=float)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    y = np.full(x.shape, np.nan)
    ok = ~np.isnan(x)
    big = ok & (x > 1e7)
    y[big] = 1 / np.sqrt(x[big])
    tiny = ok & (x >= 0) & (x < 1e-6)
    with np.errstate(divide="ignore"):
        y[tiny] = 1 / x[tiny]
    mid = ok & (x >= 1e-6) & (x <= 1e7)
    if np.any(mid):
        xm = x[mid]
        ym = 0.5 + 1 / xm
        for _ in range(51):
            tri = special.polygamma(1, ym)
            dif = tri * (1 - tri / xm) / special.polygamma(2, ym)
            ym = ym + dif
            if np.max(-dif / ym) < 1e-8:
                break
        y[mid] = ym
    return y[0] if scalar else y


def fit_f_dist(s2: np.ndarray, df1: np.ndarray) -> tuple[float, float]:
    """`(s2_prior, df_prior)` を積率で推定する（limma `fitFDist`、`covariate=NULL`）。

    `x = max(x, 1e-5·median)`、`e = log x − ψ(df1/2) + log(df1/2)`、
    `evar = Σ(e−ē)²/(n−1) − mean ψ'(df1/2)`。`evar > 0` なら `df_prior = 2·ψ'⁻¹(evar)`、
    `s2_prior = exp(ē − log(df_prior/2) + ψ(df_prior/2))`。さもなくば `df_prior = inf`、`s2_prior = mean x`。
    n=1（使える値が 1 つ）は `(その値, 0.0)`。中央値が 0 なら中央値を 1 とみなす（limma と同じ）。
    """
    x = np.atleast_1d(np.asarray(s2, dtype=float))
    d = np.atleast_1d(np.asarray(df1, dtype=float))
    n = x.size
    if n == 0:
        return (math.nan, math.nan)
    if n == 1:
        return (float(x[0]), 0.0)
    if d.size == 1:
        if not (math.isfinite(d[0]) and d[0] > 1e-15):
            return (math.nan, math.nan)
        ok = np.ones(n, dtype=bool)
    else:
        if d.size != n:
            raise ValueError("s2 と df1 の長さが違う")
        ok = np.isfinite(d) & (d > 1e-15)
    ok = ok & np.isfinite(x) & (x > -1e-15)
    nok = int(ok.sum())
    if nok == 0:
        return (math.nan, math.nan)
    if nok == 1:
        return (float(x[ok][0]), 0.0)
    x = x[ok]
    if d.size > 1:
        d = d[ok]
    x = np.maximum(x, 0.0)
    m = float(np.median(x))
    if m == 0.0:
        m = 1.0
    x = np.maximum(x, 1e-5 * m)
    e = np.log(x) + _logmdigamma(d / 2)
    emean = float(np.mean(e))
    evar = float(np.sum((e - emean) ** 2) / (nok - 1))
    evar -= float(np.mean(special.polygamma(1, d / 2)))
    if evar > 0:
        df2 = 2 * float(trigamma_inverse(np.array([evar]))[0])
        s20 = math.exp(emean - float(_logmdigamma(np.array([df2 / 2]))[0]))
    else:
        df2 = math.inf
        s20 = float(np.mean(x))
    return (s20, df2)


def _brent_fmin(f, ax: float, bx: float, tol: float) -> float:
    """R の `optimize`（C の Brent_fmin）の移植。黄金分割と放物線補間で最小点を返す。"""
    c = (3.0 - math.sqrt(5.0)) * 0.5
    eps = math.sqrt(float(np.finfo(float).eps))
    a, b = ax, bx
    v = a + c * (b - a)
    w = x = v
    d = 0.0
    e = 0.0
    fx = f(x)
    fv = fw = fx
    tol3 = tol / 3.0
    while True:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2.0
        if abs(x - xm) <= t2 - (b - a) * 0.5:
            break
        p = q = r = 0.0
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            else:
                q = -q
            r = e
            e = d
        if abs(p) >= abs(q * 0.5 * r) or p <= q * (a - x) or p >= q * (b - x):
            e = (b - x) if x < xm else (a - x)
            d = c * e
        else:
            d = p / q
            u = x + d
            if u - a < t2 or b - u < t2:
                d = tol1
                if x >= xm:
                    d = -d
        if abs(d) >= tol1:
            u = x + d
        elif d > 0.0:
            u = x + tol1
        else:
            u = x - tol1
        fu = f(u)
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x


def fit_f_dist_unequal_df1(s2: np.ndarray, df1: np.ndarray) -> tuple[float, float]:
    """自由度が揃わないときの `(s2_prior, df_prior)`（limma `fitFDistUnequalDF1`、robust=FALSE）。

    重み 1/ψ'(df1/2) の加重平均で位置を決め、事前の自由度は最尤（R の `optimize` を
    `par ∈ [1/2, 0.9998]`、`d2 = par/(1−par)` の上で）で求める。情報のある分散（正の値）が
    2 つ未満なら `(nan, nan)`。
    """
    x = np.array(np.atleast_1d(s2), dtype=float)
    d1all = np.array(np.atleast_1d(df1), dtype=float)
    n = x.size
    if d1all.size == 1:
        d1all = np.full(n, d1all[0])
    if d1all.size != n:
        raise ValueError("s2 と df1 の長さが違う")
    if np.any(np.isnan(d1all)):
        raise ValueError("df1 に NaN がある")
    weights: np.ndarray | None = None
    na = np.isnan(x)
    if na.any():
        weights = (~na).astype(float)
        x[na] = 0.0
    if d1all.min() < 0.01:
        low = d1all < 0.01
        weights = (~low).astype(float) if weights is None else np.where(low, 0.0, weights)
        d1all[low] = 1.0
    info = x > 0
    if weights is not None:
        info[weights == 0] = False
    if int(info.sum()) < 2:
        return (math.nan, math.nan)
    m = float(np.median(x[info]))
    xpos = np.maximum(x, 1e-12 * m)
    d1 = d1all / 2
    e = np.log(xpos) + _logmdigamma(d1)
    w = 1 / special.polygamma(1, d1)
    if weights is not None:
        w = w * weights
    emean = float(np.sum(w * e) / np.sum(w))
    d1x = d1 * xpos
    pw = np.ones(n) if weights is None else weights

    def minus_twice_loglik(par: float) -> float:
        d2 = par / (1 - par)
        d2s20 = d2 * math.exp(emean - float(_logmdigamma(np.array([d2]))[0]))
        ll = (-(d1 + d2) * np.log1p(d1x / d2s20) - d1 * math.log(d2s20)
              + special.gammaln(d1 + d2) - math.lgamma(d2))
        return float(-2 * np.sum(pw * ll))

    par = _brent_fmin(minus_twice_loglik, 0.5, 0.9998, _OPTIMIZE_TOL)
    d2 = par / (1 - par)
    s20 = math.exp(emean - float(_logmdigamma(np.array([d2]))[0]))
    return (s20, 2 * d2)


def squeeze_var(var: np.ndarray, df: np.ndarray) -> tuple[float, float, np.ndarray]:
    """分散を事前分布へ縮める（limma `squeezeVar`、covariate=NULL・robust=FALSE）。

    戻り値は `(s2_prior, df_prior, s2_post)`。特徴が 3 未満なら縮めない（`df_prior=0`）。
    正の自由度がすべて等しければ `fit_f_dist`、揃わなければ `fit_f_dist_unequal_df1`（limma の
    `legacy` の自動判定と同じ）。自由度 0 の特徴の分散は 0 として扱う。
    """
    v = np.array(np.atleast_1d(var), dtype=float)
    d = np.array(np.atleast_1d(df), dtype=float)
    n = v.size
    if n == 0:
        raise ValueError("var が空")
    if n < 3:
        finite = v[np.isfinite(v)]
        prior = float(np.median(finite)) if finite.size else math.nan
        return (prior, 0.0, v.copy())
    if d.size > 1:
        v[d == 0] = 0.0
    dfp = d[d > 0]
    legacy = dfp.size > 0 and float(dfp.min()) == float(dfp.max())
    if legacy:
        s2_prior, df_prior = fit_f_dist(v, d)
    else:
        s2_prior, df_prior = fit_f_dist_unequal_df1(v, d)
    if math.isnan(df_prior):
        # limma はここで止まる。分散の情報が足りないので縮めない（df_prior=0）。
        finite = v[np.isfinite(v) & (v > 0)]
        prior = float(np.median(finite)) if finite.size else 0.0
        return (prior, 0.0, v.copy())
    if math.isinf(df_prior):
        post = np.full(n, s2_prior)
    else:
        post = (d * v + df_prior * s2_prior) / (d + df_prior)
    return (float(s2_prior), float(df_prior), post)


# ---------------------------------------------------------------------------
# 経験ベイズ（eBayes）


def _f_sf(f: np.ndarray, df1: int, df2: np.ndarray) -> np.ndarray:
    """F(df1, df2) の上側確率。df2 = inf は χ²(df1)/df1 の極限（R の `pf` と同じ）。"""
    f = np.asarray(f, dtype=float)
    df2 = np.broadcast_to(np.asarray(df2, dtype=float), f.shape)
    out = np.full(f.shape, np.nan)
    inf = np.isinf(df2)
    fin = ~inf
    out[fin] = stats.f.sf(f[fin], df1, df2[fin])
    out[inf] = stats.chi2.sf(f[inf] * df1, df1)
    return out


def _classify_tests_f(t: np.ndarray, cov: np.ndarray) -> tuple[np.ndarray, int]:
    """limma `classifyTestsF(fstat.only=TRUE)`: 対比の相関で直交化した F と分子の自由度。"""
    ntests = t.shape[1]
    if ntests == 1:
        return t[:, 0] ** 2, 1
    cov = np.array(cov, dtype=float)
    diag = np.diag(cov).copy()
    diag[diag == 0] = 1.0
    np.fill_diagonal(cov, diag)
    cor = _cov2cor(cov)
    values, vectors = np.linalg.eigh(cor)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    r = int(np.sum(values / values[0] > 1e-8))
    q = vectors[:, :r] / np.sqrt(values[:r]) / math.sqrt(r)
    return np.sum((t @ q) ** 2, axis=1), r


def ebayes(fit: LinearFit, *, ci_level: float = 0.95) -> Moderated:
    """moderated t / F と区間（limma `eBayes(trend=FALSE, robust=FALSE)` と `topTable(confint=)`）。

    `df_total = min(df_residual + df_prior, Σdf_residual)`、`p = 2·P(T_{df_total} > |t|)`。
    区間は `coef ± t_{df_total,(1+level)/2}·√s2_post·stdev_unscaled`。
    F の分母の自由度は limma と同じく `df_prior + df_residual`（頭打ちにしない。`df_prior = inf` なら χ² の極限）。
    """
    if not 0 < ci_level < 1:
        raise ValueError("ci_level は 0 と 1 の間であること")
    coef = np.asarray(fit.coefficients, dtype=float)
    stdev = np.asarray(fit.stdev_unscaled, dtype=float)
    sigma = np.asarray(fit.sigma, dtype=float)
    df_res = np.asarray(fit.df_residual, dtype=float)
    if df_res.size == 0 or float(np.max(df_res)) == 0.0:
        raise ValueError("残差の自由度がどの特徴にも無い")
    if not np.any(np.isfinite(sigma)):
        raise ValueError("有限の残差標準偏差が無い")

    s2_prior, df_prior, s2_post = squeeze_var(sigma ** 2, df_res)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = coef / stdev / np.sqrt(s2_post)[:, None]
    df_total = np.minimum(df_res + df_prior, float(np.nansum(df_res)))
    p = 2 * stats.t.sf(np.abs(t), df_total[:, None])
    margin = (np.sqrt(s2_post)[:, None] * stdev
              * stats.t.ppf((1 + ci_level) / 2, df_total)[:, None])
    f, df1 = _classify_tests_f(t, fit.cov_coefficients)
    f_df2 = df_prior + df_res
    f_p = _f_sf(f, df1, f_df2)
    return Moderated(
        t=t,
        p=p,
        coefficients=coef.copy(),
        ci_low=coef - margin,
        ci_high=coef + margin,
        df_total=df_total,
        s2_post=s2_post,
        s2_prior=float(s2_prior),
        df_prior=float(df_prior),
        f=f,
        f_p=f_p,
        contrast_names=list(fit.column_names),
        f_df1=df1,
        f_df2=f_df2,
    )
