"""再解析の spec が呼ぶ検定・効果量・補正（spec 2026-10-22 §4 の 3）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。統計エンジン
（`analysis/stats/`）とは別の、論文の方法を写すための入口で、定義は呼び出しごとに名前で選ぶ:

- `student`: 等分散の t（`scipy.stats.ttest_ind(equal_var=True)`）、両側。
- `welch`: 不等分散の t、両側。
- `mannwhitney`: `scipy.stats.mannwhitneyu(alternative="two-sided")`。
- `ols`: `transform(y) ~ 1 + group(test=1) + 共変量` の最小二乗。群の係数の t の両側 p。

`transform="log2"` は有限の正値だけに当て、ほかは欠測にする（clip も pseudocount も
入れない）。

ほかに、外部の値との Spearman の ρ（`spearman`。spec の `correlate`）と、群の値が
すべて同じか（`constant_groups` / `constant_columns`。埋め値の印。統計エンジンの `compare` と共有する）。
"""
from __future__ import annotations

import functools
import itertools
import math
import warnings

import numpy as np
from scipy import stats

__all__ = ["FC_DEFINITIONS", "KINDS", "bh", "constant_columns", "constant_groups",
           "detection_driven",
           "detection_fisher_p", "fold_change", "ols_group", "spearman", "transform_values",
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


def detection_fisher_p(detected_ref: int, n_ref: int, detected_test: int, n_test: int) -> float:
    """検出の有無（検出 / 未検出）の 2×2 の Fisher の正確検定、両側 p。"""
    return float(stats.fisher_exact([[detected_test, n_test - detected_test],
                                     [detected_ref, n_ref - detected_ref]])[1])


def detection_driven(direction: str | None, detected_ref: int | None, n_ref: int,
                     detected_test: int | None, n_test: int) -> bool | None:
    """群の差が、値の大小でなく検出の有無から来ているか（機械的な当てはめで、判定は変えない）。

    高い側の群の検出が群の例数の半分に満たない（`<`）か、低い側の群の検出が 0 なら True。
    向きが無いか検出の例数が無ければ None。
    """
    if direction is None or detected_ref is None or detected_test is None:
        return None
    if direction == "test_higher":
        high, n_high, low = detected_test, n_test, detected_ref
    else:
        high, n_high, low = detected_ref, n_ref, detected_test
    return bool(high < n_high / 2 or low == 0)


#: Spearman の ρ を出す最少の組の数（これ未満は ρ を null にする）。
SPEARMAN_MIN_N = 5
#: 組がこの数以下なら、同順位の有無によらず全順列（9! = 362,880 通り）を数えた正確な並べ替えの p。
#: R の `cor.test(method="spearman")` が正確な分布を数える境界と同じ。
SPEARMAN_EXACT_MAX_N = 9
#: 組がこの数以下なら、同順位が無ければ AS 89 の Edgeworth 展開（R の `cor.test` と同じ）、あれば
#: モンテカルロの並べ替え（`SPEARMAN_PERMUTATIONS` 回）。超えれば t 近似（大きい n では十分に正確）。
SPEARMAN_EDGEWORTH_MAX_N = 1290
#: 同順位のある組の並べ替えの回数（両側 p の最小は 2 / (回数 + 1) ≈ 1e-4）と、その乱数の種（同じ入力なら同じ p）。
SPEARMAN_PERMUTATIONS = 20000
_SPEARMAN_SEED = 0
#: 並べ替えの帰無分布を一度に作る行の数（n が大きくても一時の行列を数十 MB に抑える）。
_PERMUTATION_CHUNK = 2000
#: Spearman の p の出し方（`spearman` の `p_method`）。
SPEARMAN_P_METHODS = ("exact_permutation", "edgeworth", "monte_carlo_permutation", "t_approx")
#: AS 89（Best & Roberts 1975）の Edgeworth 展開の係数（R の `prho` と同じ）。
_AS89 = (0.2274, 0.2531, 0.1745, 0.0758, 0.1033, 0.3932, 0.0879, 0.0151, 0.0072, 0.0831, 0.0131,
         4.6e-4)


@functools.lru_cache(maxsize=None)
def _permutations(n: int) -> np.ndarray:
    """0..n−1 の全順列（n! × n、int8）。n ≤ `SPEARMAN_EXACT_MAX_N` だけで呼ぶ。"""
    count = math.factorial(n)
    flat = np.fromiter(itertools.chain.from_iterable(itertools.permutations(range(n))),
                       dtype=np.int8, count=count * n)
    return flat.reshape(count, n)


@functools.lru_cache(maxsize=None)
def _s_null_counts(n: int) -> np.ndarray:
    """同順位の無い n 組の `S = Σ d²`（順位の差の 2 乗和）の帰無分布: S の値ごとの順列の数。"""
    perms = _permutations(n).astype(np.int32)
    return np.bincount(((perms - np.arange(n)) ** 2).sum(axis=1))


def _two_sided(upper: float, lower: float) -> float:
    return float(min(1.0, 2.0 * min(upper, lower)))


def _tail_counts(null: np.ndarray, observed: float) -> tuple[int, int]:
    """帰無分布のうち観測以上・観測以下の数（丸めの誤差で同じ値を取りこぼさない幅つき）。"""
    tol = 1e-9 * max(1.0, abs(observed), float(np.abs(null).max()) if null.size else 0.0)
    return int((null >= observed - tol).sum()), int((null <= observed + tol).sum())


def _exact_p(ra: np.ndarray, rb: np.ndarray, ties: bool) -> float:
    """全順列を数えた正確な両側 p（`min(1, 2·min(上側, 下側))`。R・scipy の並べ替えの両側と同じ）。

    同順位が無ければ S の帰無分布（n ごとに 1 回だけ数える）から、あれば平均の順位のまま `rb` を
    並べ替えて `Σ (ra − ā)·rb` を数える（ρ は周辺を変えない並べ替えでこの量に比例する）。
    """
    n = int(ra.size)
    if not ties:
        counts = _s_null_counts(n)
        total = float(counts.sum())
        s = int(round(float(((ra - rb) ** 2).sum())))
        # S が小さいほど正の相関。上側（S が大きい側）は P[S′ ≥ s]、下側は P[S′ ≤ s]
        return _two_sided(counts[s:].sum() / total, counts[:s + 1].sum() / total)
    centered = ra - ra.mean()
    null = rb[_permutations(n)] @ centered
    upper, lower = _tail_counts(null, float(rb @ centered))
    return _two_sided(upper / null.size, lower / null.size)


def _monte_carlo_p(ra: np.ndarray, rb: np.ndarray) -> float:
    """モンテカルロの並べ替えの両側 p（`rb` を `SPEARMAN_PERMUTATIONS` 回並べ替える。種は固定）。

    片側ごとに `(1 + 観測以上の数) / (1 + 回数)`（観測の並びも帰無分布の 1 つに数える。p が 0 にならない）を
    とり、`min(1, 2·min(上側, 下側))`。同順位は平均の順位のまま並べ替える（周辺を変えない）。
    """
    centered = ra - ra.mean()
    observed = float(rb @ centered)
    rng = np.random.default_rng(_SPEARMAN_SEED)
    upper = lower = 0
    for start in range(0, SPEARMAN_PERMUTATIONS, _PERMUTATION_CHUNK):
        size = min(_PERMUTATION_CHUNK, SPEARMAN_PERMUTATIONS - start)
        null = rng.permuted(np.tile(rb, (size, 1)), axis=1) @ centered
        up, low = _tail_counts(null, observed)
        upper, lower = upper + up, lower + low
    total = SPEARMAN_PERMUTATIONS + 1
    return _two_sided((1 + upper) / total, (1 + lower) / total)


def _as89_upper(n: int, s: int) -> float:
    """P[S′ ≥ s]（同順位の無い n 組。AS 89 の Edgeworth 展開。R の `prho` と同じ式）。"""
    if s <= 0:
        return 1.0
    if s > n * (n * n - 1) / 3:
        return 0.0
    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = _AS89
    b = 1.0 / n
    x = (6.0 * (s - 1) * b / (n * n - 1) - 1.0) * math.sqrt(1.0 / b - 1.0)
    y = x * x
    u = x * b * (c1 + b * (c2 + c3 * b) + y * (-c4 + b * (c5 + c6 * b) - y * b * (
        c7 + c8 * b - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))))
    return float(min(1.0, max(0.0, u / math.exp(y / 2.0) + stats.norm.sf(x))))


def _edgeworth_p(ra: np.ndarray, rb: np.ndarray) -> float:
    """同順位の無い組の両側 p（R の `cor.test(method="spearman")` と同じ: 観測の側の裾の 2 倍）。"""
    n = int(ra.size)
    s = int(round(float(((ra - rb) ** 2).sum())))
    if s > (n ** 3 - n) / 6:
        one = _as89_upper(n, s)
    else:
        one = 1.0 - _as89_upper(n, s + 2)        # P[S′ ≤ s]（S は偶数の値だけをとる）
    return float(min(1.0, 2.0 * one))


def _p_floor(ra: np.ndarray, rb: np.ndarray) -> float:
    """並べ替えの両側 p がとりうる最小の下限（最も極端な並びの確率の 2 倍）。

    同じ ρ を与える並びは、少なくとも一方の同順位の中の入れ替えの数だけある（同順位が無ければ 1 つ）。
    近似（Edgeworth・t）の p がこれを下回らないように抑える（|ρ| = 1 で p = 0 にしない）。
    """
    n = int(ra.size)
    log_same = max(sum(math.lgamma(int(c) + 1) for c in np.unique(r, return_counts=True)[1])
                   for r in (ra, rb))
    return float(min(1.0, 2.0 * math.exp(log_same - math.lgamma(n + 1))))


def spearman(x, y) -> dict:
    """両方が有限の組だけで Spearman の ρ（同順位は平均の順位）と両側 p。

    t 近似（自由度 n−2）は小さい n で p を小さく出しすぎる（n=5 の |ρ|=1 で 0、正確には 2/5! = 0.0167）。
    p の出し方（`p_method`）は組の数と同順位で決まる:

    - 組が `SPEARMAN_EXACT_MAX_N`（9）以下: 全順列を数えた正確な並べ替えの p（`exact_permutation`）。同順位が
      無ければ R の `cor.test(method="spearman")` と同じ値、あれば平均の順位のまま並べ替える条件付きの正確な p
      （R は同順位があると t 近似に落ちる）。
    - 10〜`SPEARMAN_EDGEWORTH_MAX_N`（1290）組で同順位が無い: AS 89 の Edgeworth 展開（`edgeworth`。R と同じ）。
    - 10〜1290 組で同順位がある: モンテカルロの並べ替え（`monte_carlo_permutation`。`SPEARMAN_PERMUTATIONS` 回、
      種は固定で同じ入力なら同じ p。最小は約 1e-4）。同順位（未検出の 0 など）があると t 近似は p を最大で
      半分ほどに小さく出す。
    - 1290 組を超える: t 近似（`t_approx`）。

    近似の 2 つ（Edgeworth・t）は、並べ替えでとりうる最小の p（`_p_floor`）で下から抑える（|ρ| = 1 でも 0 に
    しない）。組が `SPEARMAN_MIN_N` 未満なら ρ と p は null（`reason: too_few_pairs`）、どちらかの値がすべて
    同じなら null（`constant`）。どちらも `p_method` は null。
    """
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = int(a.size)
    if n < SPEARMAN_MIN_N:
        return {"n": n, "rho": None, "p": None, "reason": "too_few_pairs", "p_method": None}
    ra, rb = stats.rankdata(a), stats.rankdata(b)
    if np.all(ra == ra[0]) or np.all(rb == rb[0]):
        return {"n": n, "rho": None, "p": None, "reason": "constant", "p_method": None}
    rho = float(np.corrcoef(ra, rb)[0, 1])
    rho = max(-1.0, min(1.0, rho))
    ties = len(np.unique(ra)) < n or len(np.unique(rb)) < n
    if n <= SPEARMAN_EXACT_MAX_N:
        p, method = _exact_p(ra, rb, ties), "exact_permutation"
    elif ties and n <= SPEARMAN_EDGEWORTH_MAX_N:
        p, method = _monte_carlo_p(ra, rb), "monte_carlo_permutation"
    else:
        if n <= SPEARMAN_EDGEWORTH_MAX_N:
            p, method = _edgeworth_p(ra, rb), "edgeworth"
        else:
            t = rho * math.sqrt((n - 2) / (1.0 - rho * rho)) if abs(rho) < 1.0 else math.inf
            p, method = float(2 * stats.t.sf(abs(t), n - 2)), "t_approx"
        # 近似は裾で 0 に潰れうる（|ρ| = 1 なら t 近似は 0、Edgeworth も 0 で頭打ち）。並べ替えでとりうる最小で抑える
        p = max(p, _p_floor(ra, rb))
    return {"n": n, "rho": rho, "p": p, "reason": None, "p_method": method}


#: 「すべて同じ」の許容（相対）。値の大きさ × 1e3·eps（統計エンジンの退化の分散と同じ基準）。
CONSTANT_RTOL = 1e3 * float(np.finfo(float).eps)


def constant_columns(block, *, positive_only: bool = False) -> np.ndarray:
    """列ごとに、観測の値が 2 つ以上あり、すべて同じ（相対の丸めの誤差 `CONSTANT_RTOL` の内）か。

    観測の値は有限の値（`positive_only` なら有限の正値。0 以下を欠測とする統計エンジンの読み方）。
    `block` は（値 × 列）。1 次元なら 1 列として扱う。
    """
    x = np.asarray(block, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    if x.shape[0] == 0:
        return np.zeros(x.shape[1], dtype=bool)
    with np.errstate(invalid="ignore"):
        observed = np.isfinite(x) & (x > 0) if positive_only else np.isfinite(x)
    n = observed.sum(axis=0)
    hi = np.where(observed, x, -np.inf).max(axis=0)
    lo = np.where(observed, x, np.inf).min(axis=0)
    scale = np.where(observed, np.abs(x), 0.0).max(axis=0)
    with np.errstate(invalid="ignore"):
        return (n >= 2) & (hi - lo <= CONSTANT_RTOL * scale)


def constant_groups(values_by_group: dict[str, np.ndarray], *,
                    positive_only: bool = False) -> list[str] | None:
    """観測の値が 2 つ以上あり、すべて同じ（`constant_columns`）群の名前（埋め値の印）。無ければ None。

    渡すのは**読んだ値**（試料・注入ごと。正規化・ドリフト補正・反復のまとめの前）。片群一定の印
    （`constant_in_group`）の意味は v1（analysis_spec）と統計エンジン（stats/compare）で同じ。
    """
    out = [label for label, raw in values_by_group.items()
           if bool(constant_columns(np.asarray(raw, dtype=float).ravel(),
                                    positive_only=positive_only)[0])]
    return out or None
