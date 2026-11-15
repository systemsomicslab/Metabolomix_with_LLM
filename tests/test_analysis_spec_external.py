"""analysis-spec.v1 の外部の表（MAF・mwTab・xlsx）・`correlate`・和の寄与・群の中で一定の印・
効かない設定の notes・aggregate と features.feature_ids の併用の誤り。fixture は合成のファイル。"""
from __future__ import annotations

import functools
import itertools
import math

import numpy as np
import pytest
from scipy import stats

from metabolomix.analysis import analysis_spec as A
from metabolomix.analysis import spec_statistics as S
from metabolomix.core.atomic_io import DomainError
from tests.external_table_fixtures import write_maf, write_mwtab, write_xlsx

SAMPLES = ["HC1", "HC2", "HC3", "RA1", "RA2", "RA3", "QC1", "QC2"]
VALUES = {
    "F1": [10, 12, 11, 20, 22, 21, 15, 16],
    "F2": [5, 6, 5, 5, 6, 5, 5, 5],
    "F3": [0, 0, 0, 1, 0, 0, 0, 0],          # HC の 3 試料がすべて 0（群の中で一定）
    "F4": [100, 110, 90, 50, 55, 45, 75, 76],
}
#: 外部の表の値（試料は SAMPLES の順）。
GLUCOSE = [1, 3, 2, 5, 8, 6, 4, 7]
LACTATE = [9, 8, 7, 6, 5, 4, 3, 2]


def _write(path, header, rows):
    path.write_text("\n".join(["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]) + "\n",
                    encoding="utf-8")
    return path


@pytest.fixture
def files(tmp_path):
    matrix = _write(tmp_path / "matrix.tsv", ["feature_id", *SAMPLES],
                    [[f, *v] for f, v in VALUES.items()])
    groups = _write(tmp_path / "groups.tsv", ["sample", "group", "role"],
                    [[s, s[:2] if not s.startswith("QC") else "", "qc" if s.startswith("QC") else "sample"]
                     for s in SAMPLES])
    mwtab = write_mwtab(tmp_path / "paper.txt", SAMPLES,
                        [("Glucose", GLUCOSE), ("Lactate", LACTATE)])
    return {"matrix": matrix, "groups": groups, "mwtab": mwtab, "dir": tmp_path}


def _spec(files, analyses=None, **extra):
    spec = {"input": {"source": "tsv", "matrix_path": str(files["matrix"])},
            "groups": {"source": "tsv", "path": str(files["groups"]), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": analyses or [{"analysis_id": "main", "test": {"kind": "student"}}]}
    spec.update(extra)
    return spec


def _run(spec):
    return A.run_spec(A.normalize_spec(spec))


def _rows(result, analysis=0):
    return {r["feature_id"]: r for r in result["analyses"][analysis]["rows"]}


# ---------- Spearman と群の中で一定の印（spec_statistics） ----------

def _exact_spearman_p(x, y):
    """両側の正確な並べ替えの p（独立の計算: y の平均の順位の全 n! の並びの ρ を数える）。"""
    ra, rb = stats.rankdata(x), stats.rankdata(y)
    perms = np.array(list(itertools.permutations(range(len(ra)))))
    shuffled = rb[perms]
    centered = shuffled - shuffled.mean(axis=1, keepdims=True)
    rhos = (centered @ (ra - ra.mean())) / (np.linalg.norm(centered, axis=1) * np.linalg.norm(ra - ra.mean()))
    observed = float(np.corrcoef(ra, rb)[0, 1])
    upper, lower = (rhos >= observed - 1e-9).mean(), (rhos <= observed + 1e-9).mean()
    return min(1.0, 2.0 * min(upper, lower))


#: R 4.6.0 の `cor.test(x, y, method = "spearman")$p.value`（同順位なし）。9 組以下は R も正確な分布。
R_SPEARMAN_EXACT = [
    ([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], 0.01666666667),
    ([1, 2, 3, 5, 4], [1, 2, 3, 4, 5], 0.08333333333),
    ([2, 1, 3, 4, 6, 5], [1, 2, 3, 4, 5, 6], 0.03333333333),
    ([2, 1, 3, 4, 5, 7, 6], list(range(1, 8)), 0.006746031746),
    ([2, 1, 3, 4, 5, 6, 8, 7], list(range(1, 9)), 0.001140873016),
    ([3, 1, 2, 5, 4, 7, 6, 9, 8], list(range(1, 10)), 0.002028218695),
]
#: 同じく 15 組以上で p ≥ 1e-3（R も AS 89 の Edgeworth 展開。この範囲では正確な分布との差が数 % 以内）。
R_SPEARMAN_EDGEWORTH = [
    ([5, 1, 9, 8, 4, 7, 3, 10, 2, 6, 12, 11, 15, 13, 14], list(range(1, 16)), 0.004315957134),
    ([11, 14, 15, 12, 13, 2, 7, 8, 9, 10, 5, 4, 6, 1, 3], list(range(1, 16)), 0.00108003331),
    ([7, 2, 13, 1, 4, 11, 6, 5, 8, 3, 12, 16, 15, 9, 10, 14], list(range(1, 17)), 0.02957777048),
    ([10, 1, 2, 4, 5, 19, 8, 6, 12, 3, 7, 9, 18, 17, 11, 14, 15, 16, 20, 13], list(range(1, 21)),
     0.002001916004),
    ([23, 8, 10, 22, 21, 16, 20, 9, 18, 25, 24, 7, 14, 3, 5, 13, 19, 12, 15, 11, 1, 17, 4, 6, 2],
     list(range(1, 26)), 0.01005877021),
    ([19, 23, 15, 10, 25, 30, 11, 3, 2, 1, 8, 5, 4, 22, 14, 31, 29, 21, 6, 12, 17, 9, 38, 27, 13, 28, 35,
      26, 36, 33, 20, 24, 32, 7, 40, 16, 18, 37, 39, 34], list(range(1, 41)), 0.001198919954),
]


@pytest.mark.parametrize("x, y, r_p", R_SPEARMAN_EXACT)
def test_spearman_p_up_to_nine_pairs_matches_r_exact(x, y, r_p):
    out = S.spearman(x, y)
    assert out["rho"] == pytest.approx(stats.spearmanr(x, y).statistic)
    assert out["p_method"] == "exact_permutation"
    assert out["p"] == pytest.approx(r_p, rel=1e-8)
    assert out["p"] == pytest.approx(_exact_spearman_p(x, y), rel=1e-12)


@pytest.mark.parametrize("x, y, r_p", R_SPEARMAN_EDGEWORTH)
def test_spearman_p_from_fifteen_pairs_matches_r_edgeworth(x, y, r_p):
    """15 組以上で同順位が無く p ≥ 1e-3 なら R の cor.test と同じ AS 89 の Edgeworth 展開（負の ρ も）。"""
    out = S.spearman(x, y)
    assert out["rho"] == pytest.approx(stats.spearmanr(x, y).statistic)
    assert out["p_method"] == "edgeworth"
    assert out["p"] == pytest.approx(r_p, rel=1e-8)


def _count_s_at_most(n, s_max):
    """`S = Σ (i − π(i))² ≤ s_max` の順列 π の数（独立の計算: 深さ優先で数え、途中の和が超えた枝を切る）。"""
    used = [False] * n

    def walk(position, partial):
        if position == n:
            return 1
        total = 0
        for value in range(n):
            step = partial + (position - value) ** 2
            if not used[value] and step <= s_max:
                used[value] = True
                total += walk(position + 1, step)
                used[value] = False
        return total

    return walk(0, 0)


def _tail_p_by_count(x, y):
    """同順位の無い組の正確な両側 p（`_count_s_at_most` で観測の側の裾を数える。S の分布は対称）。"""
    n = len(x)
    s = int(round(float(((stats.rankdata(x) - stats.rankdata(y)) ** 2).sum())))
    s_top = n * (n * n - 1) // 3
    return min(1.0, 2 * _count_s_at_most(n, min(s, s_top - s)) / math.factorial(n))


@functools.lru_cache(maxsize=None)
def _s_distribution_by_enumeration(n):
    """S の帰無分布（独立の計算: 先頭の順位ごとに残りの (n−1)! 通りを itertools で並べて数える）。"""
    counts = np.zeros(n * (n * n - 1) // 3 + 1, dtype=np.int64)
    for first in range(n):
        rest = [v for v in range(n) if v != first]
        perms = np.fromiter(itertools.chain.from_iterable(itertools.permutations(rest)), dtype=np.int8,
                            count=math.factorial(n - 1) * (n - 1)).reshape(-1, n - 1).astype(np.int32)
        counts += np.bincount(first ** 2 + ((perms - np.arange(1, n)) ** 2).sum(axis=1),
                              minlength=counts.size)
    return counts


def test_spearman_null_distribution_without_ties_is_exact():
    """S の正確な分布（部分集合の動的計画法）は、10 組で全 10! 通りの数え上げと一致し、11〜14 組で
    合計 n!・Var(ρ) = 1/(n−1)・E[ρ⁴] = 3(25n³ − 38n² − 35n + 72) / (25n(n+1)(n−1)³)（Kendall の閉じた形）。"""
    by_enumeration = _s_distribution_by_enumeration(10)
    lower_half = S._s_null_lower_half(10)
    assert np.array_equal(lower_half, by_enumeration[:lower_half.size])
    for n in range(10, S.SPEARMAN_EXACT_NO_TIES_MAX_N + 1):
        lower_half = S._s_null_lower_half(n)
        full = np.concatenate([lower_half, lower_half[-2::-1]])        # n(n²−1)/6 を中心に対称
        assert int(full.sum()) == math.factorial(n)
        weights = full / full.sum()
        rho = 1 - 6 * np.arange(full.size) / (n ** 3 - n)
        assert float(weights @ rho ** 2) == pytest.approx(1 / (n - 1), rel=1e-12)
        fourth = 3 * (25 * n ** 3 - 38 * n ** 2 - 35 * n + 72) / (25 * n * (n + 1) * (n - 1) ** 3)
        assert float(weights @ rho ** 4) == pytest.approx(fourth, rel=1e-12)


@pytest.mark.parametrize("x, y", [
    # 10 組の ρ = 0.95: Edgeworth（R の cor.test も）は裾で潰れ、正確な p 1.08e-4 の 1/200（R は 0）
    (list(range(1, 11)), [2, 1, 4, 3, 6, 5, 7, 8, 10, 9]),
    ([3, 1, 2, 5, 4, 7, 6, 10, 9, 8], list(range(1, 11))),
    ([3, 1, 2, 5, 4, 7, 6, 10, 9, 8], list(range(10, 0, -1))),
    ([5, 1, 2, 8, 4, 7, 3, 10, 9, 6], list(range(1, 11))),
    (list(range(1, 13)), [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11]),
    (list(range(1, 15)), [13, 14, 12, 11, 9, 10, 8, 7, 5, 6, 4, 3, 1, 2]),   # 負の ρ
])
def test_spearman_without_ties_up_to_fourteen_pairs_is_exact(x, y):
    out = S.spearman(x, y)
    assert out["rho"] == pytest.approx(stats.spearmanr(x, y).statistic)
    assert out["p_method"] == "exact_permutation"
    if len(x) == 10:
        n = 10
        counts = _s_distribution_by_enumeration(n)
        s = int(round(float(((stats.rankdata(x) - stats.rankdata(y)) ** 2).sum())))
        s_top = n * (n * n - 1) // 3
        expected = min(1.0, 2 * counts[:min(s, s_top - s) + 1].sum() / math.factorial(n))
    else:
        expected = _tail_p_by_count(x, y)
    assert out["p"] == pytest.approx(expected, rel=1e-12)


def test_spearman_strong_correlation_is_not_the_collapsed_edgeworth_tail():
    """10 組の |ρ| ≥ 0.94 は旧来どれも 5.51e-7（Edgeworth が 0 に潰れ 2/10! で抑えた値）。正確な p を数える。"""
    x, y = list(range(1, 11)), [2, 1, 4, 3, 6, 5, 7, 8, 10, 9]
    out = S.spearman(x, y)
    exact = 2 * _count_s_at_most(10, 8) / math.factorial(10)     # S = 8 以下の並びは 196 通り
    assert out["p"] == pytest.approx(exact, rel=1e-12) and out["p"] == pytest.approx(1.0802e-4, rel=1e-4)
    assert S._edgeworth_p(stats.rankdata(x), stats.rankdata(y)) < exact / 100
    assert stats.spearmanr(x, y).pvalue < exact / 4                 # t 近似も小さく出しすぎる


def test_spearman_approximation_tail_falls_back_to_the_permutation():
    """近似（Edgeworth・t）の両側 p が 1e-3 未満なら並べ替えの p（裾で何桁も小さく出す・0 に潰れるため）。"""
    # 15 組・S = 106: 正確な p 4.34e-4（S の正確な分布）に対し Edgeworth は 3.86e-4。並べ替えの推定は正確な p に近い
    y15 = list(range(15))
    y15[0], y15[7], y15[8], y15[10] = y15[7], y15[0], y15[10], y15[8]
    out = S.spearman(list(range(15)), y15)
    lower_half = S._s_null_lower_half(15)
    exact = 2 * int(lower_half[:107].sum()) / math.factorial(15)
    assert out["p_method"] == "monte_carlo_permutation"
    assert S._edgeworth_p(np.arange(15.0), np.array(y15, dtype=float)) < S.SPEARMAN_APPROX_MIN_P
    assert out["p"] == pytest.approx(exact, rel=0.25)
    # 20 組・S = 192: Edgeworth は 8.2e-19、t 近似は 1.4e-6。2×10⁷ 回の無作為な並べ替えでは約 4e-6。
    # 並べ替えの p はとりうる最小 2/(回数+1) を下回らない（本当の p はそれ以下＝上限として読む）
    y20 = list(range(20))
    for a, b in ((0, 4), (1, 5), (2, 6), (3, 7), (8, 12), (9, 13)):
        y20[a], y20[b] = y20[b], y20[a]
    strong = S.spearman(list(range(20)), y20)
    assert strong["p_method"] == "monte_carlo_permutation"
    assert strong["p"] >= 2 / (S.SPEARMAN_PERMUTATIONS + 1) > 4e-6
    assert strong["p"] > stats.spearmanr(list(range(20)), y20).pvalue
    assert S._edgeworth_p(np.arange(20.0), np.array(y20, dtype=float)) < 1e-15
    # 20 組・R の cor.test は 1.98e-5（2×10⁷ 回の並べ替えでは約 4.0e-5＝R が半分に小さい）
    r_row = S.spearman([10, 1, 9, 8, 4, 7, 3, 5, 2, 6, 12, 11, 15, 13, 14, 20, 17, 16, 19, 18],
                       list(range(1, 21)))
    assert r_row["p_method"] == "monte_carlo_permutation" and r_row["p"] > 1.98e-5
    # 200 組の ρ = −0.78: Edgeworth は 0、2/200! も 0 に潰れる。p は 0 にしない
    rng = np.random.default_rng(1)
    x200 = np.arange(200.0)
    y200 = -(x200 + rng.normal(0, 50, 200))
    wide = S.spearman(x200, y200)
    assert wide["rho"] < -0.7 and wide["p_method"] == "monte_carlo_permutation"
    assert wide["p"] == pytest.approx(2 / (S.SPEARMAN_PERMUTATIONS + 1))


def test_spearman_beyond_the_edgeworth_range_uses_t_or_the_permutation(monkeypatch):
    """1290 組を超えれば t 近似。ただし両側 p が 1e-3 未満（|ρ| = 1 の 0 など）なら並べ替え。"""
    monkeypatch.setattr(S, "SPEARMAN_PERMUTATIONS", 2000)          # 1300 組の並べ替えを軽くする
    rng = np.random.default_rng(3)
    x = np.arange(1300.0)
    y = x + rng.normal(0, 6000, 1300)
    moderate = S.spearman(x, y)
    assert moderate["p_method"] == "t_approx"
    assert moderate["p"] == pytest.approx(stats.spearmanr(x, y).pvalue, rel=1e-9)
    assert moderate["p"] >= S.SPEARMAN_APPROX_MIN_P
    perfect = S.spearman(x, 2 * x)
    assert perfect["p_method"] == "monte_carlo_permutation" and perfect["p"] == pytest.approx(2 / 2001)


def test_spearman_small_n_is_not_the_t_approximation():
    """t 近似は小さい n で p を小さく出しすぎる: n=5 の |ρ|=1 は 0 でなく 2/5!、ρ=0.9 は 0.05 を超える。"""
    perfect = S.spearman([1, 2, 3, 4, 5], [10, 20, 30, 40, 50])
    assert perfect["rho"] == pytest.approx(1.0) and perfect["p"] == pytest.approx(2 / math.factorial(5))
    reverse = S.spearman([1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1])
    assert reverse["rho"] == pytest.approx(-1.0) and reverse["p"] == pytest.approx(2 / math.factorial(6))
    swap = S.spearman([1, 2, 3, 5, 4], [1, 2, 3, 4, 5])
    assert swap["p"] == pytest.approx(1 / 12) and swap["p"] > 0.05
    assert stats.spearmanr([1, 2, 3, 5, 4], [1, 2, 3, 4, 5]).pvalue < 0.05   # t 近似は「有意」に見せる
    # 14 組までの同順位なしの |ρ| = 1 も正確な 2/n!（0 にしない）
    big = S.spearman(list(range(12)), list(range(12)))
    assert big["p_method"] == "exact_permutation"
    assert big["p"] == pytest.approx(2 / math.factorial(12), rel=1e-12) and big["p"] > 0


def test_spearman_small_n_with_ties_is_the_exact_conditional_permutation():
    x = [1.0, 2.0, 2.0, 4.0, 5.0, 7.0, math.nan]
    y = [2.0, 1.0, 4.0, 3.0, 6.0, 5.0, 1.0]
    out = S.spearman(x, y)
    assert out["n"] == 6 and out["reason"] is None and out["p_method"] == "exact_permutation"
    assert out["rho"] == pytest.approx(stats.spearmanr(x[:6], y[:6]).statistic)
    assert out["p"] == pytest.approx(_exact_spearman_p(x[:6], y[:6]), rel=1e-12)
    zeros = S.spearman([1, 2, 3, 4, 5, 6, 7, 8], [0, 0, 0, 1, 3, 2, 5, 4])   # 未検出の 0 の同順位
    assert zeros["p"] == pytest.approx(_exact_spearman_p([1, 2, 3, 4, 5, 6, 7, 8],
                                                          [0, 0, 0, 1, 3, 2, 5, 4]), rel=1e-12)


def _exact_p_with_zero_ties(x, y):
    """y の 0（同順位）の置き場所と残りの並びをすべて数えた正確な両側 p（区別できる並びは等確率）。"""
    ra, rb = stats.rankdata(x), stats.rankdata(y)
    n, low = len(ra), rb.min()
    others = [v for v in rb if v != low]
    centered = ra - ra.mean()
    null = []
    for zeros in itertools.combinations(range(n), n - len(others)):
        rest = [i for i in range(n) if i not in zeros]
        for perm in itertools.permutations(others):
            arrangement = np.full(n, low)
            arrangement[rest] = perm
            null.append(float(arrangement @ centered))
    null = np.array(null)
    observed = float(rb @ centered)
    return min(1.0, 2.0 * min((null >= observed - 1e-9).mean(), (null <= observed + 1e-9).mean()))


def test_spearman_ties_from_ten_pairs_use_a_seeded_permutation():
    """10 組以上で同順位（未検出の 0 など）があれば t 近似でなくモンテカルロの並べ替え（種は固定）。"""
    x = list(range(1, 11))
    y = [0, 0, 0, 2, 0, 0, 5, 1, 4, 3]
    out = S.spearman(x, y)
    exact = _exact_p_with_zero_ties(x, y)                   # 10!/5! = 30,240 通りを数える
    assert out["p_method"] == "monte_carlo_permutation"
    assert out["p"] == pytest.approx(exact, rel=0.15)        # 100,000 回の並べ替えの誤差の内
    assert stats.spearmanr(x, y).pvalue < 0.8 * exact        # t 近似は小さく出しすぎる
    assert S.spearman(x, y) == out                           # 同じ入力なら同じ p
    perfect = S.spearman(x, [0, 0, 0, 1, 2, 3, 4, 5, 6, 7])
    assert perfect["p"] >= 2 / (S.SPEARMAN_PERMUTATIONS + 1)


def test_spearman_too_few_pairs_and_constant_are_null_with_reason():
    assert S.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == {
        "n": 4, "rho": None, "p": None, "reason": "too_few_pairs", "p_method": None}
    out = S.spearman([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])
    assert out["reason"] == "constant" and out["p_method"] is None


def test_constant_groups():
    assert S.constant_groups({"A": [0, 0, 0], "B": [1, 0, 0]}) == ["A"]
    assert S.constant_groups({"A": [5, math.nan, 5], "B": [math.nan, 2, math.nan]}) == ["A"]
    assert S.constant_groups({"A": [1, 2], "B": [3, 4]}) is None
    # 0 以下を欠測とする読み方（log2 の検定・統計エンジン）では 0 の並びは観測の値でない
    assert S.constant_groups({"A": [0, 0, 0], "B": [4, 0, 4]}, positive_only=True) == ["B"]
    # 丸めの誤差の内は同じ値（相対 1e3·eps）、それを超えれば違う値
    assert S.constant_groups({"A": [1000.0, 1000.0 * (1 + 1e-15)]}) == ["A"]
    assert S.constant_groups({"A": [1000.0, 1000.0 * (1 + 1e-9)]}) is None
    flags = S.constant_columns(np.array([[7.0, 0.0, 1.0], [7.0, 0.0, 2.0], [np.nan, 0.0, 1.0]]),
                               positive_only=True)
    assert flags.tolist() == [True, False, False]


def test_constant_in_group_marks_rows_without_changing_the_test(files):
    result = _run(_spec(files))
    rows = _rows(result)
    # log2 の検定では 0 は欠測（検定に入らない）。HC の 0 だけの F3 は観測の値が無いので印を付けない
    assert rows["F3"]["constant_in_group"] is None and rows["F3"]["status"] == "not_testable"
    assert rows["F1"]["constant_in_group"] is None
    assert list(rows["F1"])[-1] == "constant_in_group"        # 行の最後の欄
    # 変換なしの検定では 0 が値のまま入るので、HC の 0 の並びは片群一定
    linear = _rows(_run(_spec(files, [{"analysis_id": "lin", "test": {"kind": "student",
                                                                    "transform": "none"}}])))
    assert linear["F3"]["constant_in_group"] == ["HC"] and linear["F1"]["constant_in_group"] is None


FILL_VALUES = {
    "F1": [10, 12, 11, 20, 22, 21, 15, 16],
    "F2": [5, 6, 5, 5, 6, 5, 5, 5],
    "FILL": [7, 7, 7, 30, 35, 28, 20, 21],     # HC を同じ値（7）で埋めた
    "ZERO": [0, 0, 0, 30, 35, 28, 20, 21],     # HC が未検出（0）。補完の half_min が同じ値で埋める
    "F4": [100, 110, 90, 50, 55, 45, 75, 76],
}


def test_constant_in_group_reads_values_before_normalization(files):
    """埋め値は読んだ値（recipe の前）と補完の直後の値で探す。正規化（PQN）の試料ごとの係数は埋め値を散らす。"""
    from metabolomix.analysis import recipe_matrix

    matrix = _write(files["dir"] / "fill.tsv", ["feature_id", *SAMPLES],
                    [[f, *v] for f, v in FILL_VALUES.items()])
    recipe = {"steps": ["impute", "normalize"], "impute": {"method": "half_min"},
              "normalize": {"method": "pqn"}}
    spec = _spec(files, [{"analysis_id": "pqn", "test": {"kind": "student"}, "matrix": recipe}],
                 input={"source": "tsv", "matrix_path": str(matrix)})
    rows = _rows(_run(spec))
    assert rows["FILL"]["constant_in_group"] == ["HC"]
    assert rows["ZERO"]["constant_in_group"] == ["HC"]
    assert rows["F1"]["constant_in_group"] is None and rows["F4"]["constant_in_group"] is None
    # 検定に入れる値（recipe の後）では HC の埋め値がもう同じでない（ここで探すと見落とす）
    values = np.array(list(FILL_VALUES.values()), dtype=float).T
    roles = {s: ("qc" if s.startswith("QC") else "sample") for s in SAMPLES}
    built = recipe_matrix.build_matrix(values, SAMPLES, list(FILL_VALUES), roles,
                                       recipe_matrix.normalize_recipe(recipe))
    after = built["values"][:3, [built["feature_ids"].index(f) for f in ("FILL", "ZERO")]]
    assert np.ptp(after, axis=0).min() > 0.01
    assert S.constant_groups({"HC": after[:, 0]}) is None and S.constant_groups({"HC": after[:, 1]}) is None


# ---------- aggregate（和）の寄与・E9 ----------

def test_sum_aggregate_returns_member_contributions(files):
    spec = _spec(files, aggregate={"kind": "sum", "label": "total", "feature_ids": ["F1", "F2", "F4"]})
    info = _run(spec)["analyses"][0]["aggregate"]
    assert info["n_samples_for_shares"] == 6                 # 比べる 2 群の試料だけ（QC は入れない）
    assert [c["feature_id"] for c in info["contributions"]] == ["F1", "F2", "F4"]
    block = np.array([VALUES[f][:6] for f in ("F1", "F2", "F4")], dtype=float)
    shares = block / block.sum(axis=0)
    by_id = {c["feature_id"]: c for c in info["contributions"]}
    for k, fid in enumerate(("F1", "F2", "F4")):
        assert by_id[fid]["share_median"] == pytest.approx(np.median(shares[k]))
        assert by_id[fid]["share_median_ref"] == pytest.approx(np.median(shares[k][:3]))
        assert by_id[fid]["share_median_test"] == pytest.approx(np.median(shares[k][3:]))
        assert by_id[fid]["share_max"] == pytest.approx(shares[k].max())
    assert "_block" not in info
    # 試料ごとの割合の和は 1（寄与は和を分けたもの）
    assert np.allclose(shares.sum(axis=0), 1.0)


def test_ratio_aggregate_has_no_contributions(files):
    spec = _spec(files, aggregate={"kind": "ratio", "label": "F1/F4", "numerator": "F1",
                                   "denominator": "F4"})
    info = _run(spec)["analyses"][0]["aggregate"]
    assert "contributions" not in info and "_block" not in info


def test_aggregate_with_feature_ids_is_an_error(files):
    spec = _spec(files, aggregate={"kind": "sum", "label": "total", "feature_ids": ["F1", "F2"]},
                 features={"feature_ids": ["F1"]})
    with pytest.raises(DomainError) as info:
        A.normalize_spec(spec)
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == "features"
    assert "aggregate" in info.value.message


# ---------- 効かない設定は notes ----------

def test_covariates_on_a_non_ols_test_are_dropped_with_a_note(files):
    covariates = _write(files["dir"] / "cov.tsv", ["sample", "age"],
                        [[s, 30 + i] for i, s in enumerate(SAMPLES)])
    cov = {"path": str(covariates), "columns": [{"name": "age", "type": "numeric"}]}
    with_cov = _spec(files, [{"analysis_id": "main", "test": {"kind": "welch", "covariates": ["age"]}}],
                     covariates=cov)
    notes: list[str] = []
    normalized = A.normalize_spec(with_cov, notes=notes)
    assert normalized["analyses"][0]["test"]["covariates"] == []
    assert "analyses[0].test.covariates" in notes[0] and "ols" in notes[0]
    # 共変量で調整する解析が無いので、covariates の表そのものも効かない（spec には残し、読まない）
    assert len(notes) == 2 and notes[1].startswith("covariates は")
    plain = A.normalize_spec(_spec(files, [{"analysis_id": "main", "test": {"kind": "welch"}}],
                                   covariates=cov))
    assert A.rule_hash(normalized) == A.rule_hash(plain)
    # 知らない共変量はこれまでどおり誤り（効かない設定でも名前は確かめる、ではなく ols でだけ確かめる）
    ols = _spec(files, [{"analysis_id": "main", "test": {"kind": "ols", "covariates": ["sex"]}}],
                covariates=cov)
    with pytest.raises(DomainError):
        A.normalize_spec(ols)


def test_recipe_settings_not_in_steps_are_notes_not_errors(files):
    notes: list[str] = []
    spec = _spec(files, [{"analysis_id": "main", "test": {"kind": "student"},
                          "matrix": {"steps": [], "normalize": {"method": "pqn"}}}])
    normalized = A.normalize_spec(spec, notes=notes)
    assert normalized["analyses"][0]["matrix"] == {"schema": "matrix-recipe.v2", "steps": []}
    assert notes and notes[0].startswith("analyses[0].matrix:")


def test_correction_none_with_bh_population_and_q_threshold_is_a_note(files):
    """C7: 補正しない解析の bh_population と thresholds.q は効かない。spec には残し（rule_hash は前と同じ）、
    notes に書く。q は null なので、q の閾値は p・fc を満たす行でも満たせない（meets_thresholds=false）。"""
    analysis = {"analysis_id": "main", "test": {"kind": "student"}, "correction": "none",
                "bh_population": "reported", "thresholds": {"p": 0.05, "q": 0.05}}
    notes: list[str] = []
    normalized = A.normalize_spec(_spec(files, [analysis], features={"feature_ids": ["F1", "F4"]}),
                                  notes=notes)
    a = normalized["analyses"][0]
    assert (a["correction"], a["bh_population"], a["thresholds"]["q"]) == ("none", "reported", 0.05)
    assert len(notes) == 2
    assert notes[0].startswith("analyses[0].bh_population") and "correction=none" in notes[0]
    assert notes[1].startswith("analyses[0].thresholds.q") and "meets_thresholds" in notes[1]
    row = _rows(A.run_spec(normalized))["F1"]
    p = stats.ttest_ind(np.log2(VALUES["F1"][3:6]), np.log2(VALUES["F1"][:3])).pvalue   # 独立の計算
    assert row["p"] == pytest.approx(p) and p < 0.05 and row["q"] is None
    assert row["meets_thresholds"] is False                                 # 書いたとおり満たせない
    # 既定の bh_population・q の閾値の無い補正なしの解析は何も言わない
    quiet: list[str] = []
    A.normalize_spec(_spec(files, [{"analysis_id": "main", "test": {"kind": "student"}, "correction": "none",
                                    "thresholds": {"p": 0.05}}]), notes=quiet)
    assert quiet == []


def test_unused_injection_order_and_covariates_are_notes_and_not_read(files):
    """C7: drift_correct の無い spec の injection_order と、共変量で調整する ols の無い spec の covariates は効かない。
    表は読まず（欠けた試料・無いファイルで断らない）、notes に書く。spec には残す（rule_hash は前と同じ）。"""
    partial = _write(files["dir"] / "order.tsv", ["sample", "order"],
                     [[s, i] for i, s in enumerate(SAMPLES[:-1], start=1)])        # QC2 が欠ける
    order = {"source": "tsv", "path": str(partial)}
    cov = {"path": str(files["dir"] / "no_such_covariates.tsv"),
           "columns": [{"name": "age", "type": "numeric"}]}
    spec = _spec(files, injection_order=order, covariates=cov)
    notes: list[str] = []
    normalized = A.normalize_spec(spec, notes=notes)
    assert normalized["injection_order"]["path"] == str(partial) and normalized["covariates"]["path"]
    assert [n.split(" ")[0] for n in notes] == ["injection_order", "covariates"]
    result = A.run_spec(normalized)
    assert set(result["input_fingerprint"]["files"]) == {"groups", "matrix"}
    # 効かない表を足しても数字は同じ
    plain = A.run_spec(A.normalize_spec(_spec(files)))
    assert result["analyses"][0]["rows"] == plain["analyses"][0]["rows"]
    # drift_correct のある解析があれば、その表は従来どおり読み・確かめる（欠けた試料は SAMPLE_MISMATCH）
    drift = {"analysis_id": "main", "test": {"kind": "student"},
             "matrix": {"steps": ["drift_correct"]}}
    used = A.normalize_spec(_spec(files, [drift], injection_order=order), notes=(used_notes := []))
    assert not any(n.startswith("injection_order") for n in used_notes)
    with pytest.raises(DomainError) as exc:
        A.run_spec(used)
    assert exc.value.code == "SAMPLE_MISMATCH"


def test_alignment_qc_on_external_input_is_a_note(files):
    notes: list[str] = []
    normalized = A.normalize_spec(_spec(files, features={"all": True, "alignment_qc": "x.json"}),
                                  notes=notes)
    assert "alignment_qc" not in normalized["features"]
    assert any("外部の行列" in n for n in notes)


def test_unused_dataset_only_settings_on_external_input_are_notes_not_errors(files):
    """最終検証 WE: 外部の行列で効かない `correction=none` の `bh_population=annotated` と、drift_correct の無い spec の
    `injection_order.source=metadata` は、自分の notes が「効かない」と書くとおり断らずに進む。正規化した spec には
    書いたまま残り（`rule_hash` は前に通った spec と同じ規則）、数字は書かないときと同じ。使うとき（correction=bh・
    drift_correct のあるとき）だけ従来どおり ANALYSIS_SPEC_INVALID。"""
    analysis = {"analysis_id": "main", "test": {"kind": "student"}, "correction": "none",
                "bh_population": "annotated"}
    notes: list[str] = []
    normalized = A.normalize_spec(_spec(files, [analysis], injection_order={"source": "metadata"}), notes=notes)
    assert normalized["analyses"][0]["bh_population"] == "annotated"
    assert normalized["injection_order"] == {"source": "metadata"}
    assert [n.split(" ")[0] for n in notes] == ["analyses[0].bh_population", "injection_order"]
    result = A.run_spec(normalized)
    rows = _rows(result)
    p = stats.ttest_ind(np.log2(VALUES["F1"][3:6]), np.log2(VALUES["F1"][:3])).pvalue   # 独立の計算
    assert rows["F1"]["p"] == pytest.approx(p) and rows["F1"]["q"] is None
    assert result["analyses"][0]["bh_population_n"] is None
    assert set(result["input_fingerprint"]["files"]) == {"groups", "matrix"}       # 注入順は読まない
    plain = A.run_spec(A.normalize_spec(_spec(files, [{**analysis, "bh_population": "matrix"}])))
    assert result["analyses"][0]["rows"] == plain["analyses"][0]["rows"]
    # 注入順だけ（既定の correction=bh の解析）でも同じ: 断らず notes、数字は書かないときと同じ
    order_only = A.normalize_spec(_spec(files, injection_order={"source": "metadata"}),
                                  notes=(order_notes := []))
    assert order_only["injection_order"] == {"source": "metadata"}
    assert [n.split(" ")[0] for n in order_notes] == ["injection_order"]
    assert A.run_spec(order_only)["analyses"][0]["rows"] == _run(_spec(files))["analyses"][0]["rows"]
    # 使うときは従来どおり誤り（外部の行列に SME も試料の表も無い）
    for spec, field in (
            (_spec(files, [{**analysis, "correction": "bh"}]), "analyses"),
            (_spec(files, [{"analysis_id": "main", "test": {"kind": "student"},
                            "matrix": {"steps": ["drift_correct"]}}],
                   injection_order={"source": "metadata"}), "injection_order")):
        with pytest.raises(DomainError) as info:
            A.normalize_spec(spec)
        assert info.value.code == "ANALYSIS_SPEC_INVALID" and info.value.details["field"] == field


def test_unused_dataset_only_settings_on_external_input_reach_the_tool(files):
    """dataset_analysis_run（v1）も同じ spec を断らず、notes を返す（最終検証 WE）。"""
    import json

    from metabolomix.tools import reanalysis_tools as R

    spec = _spec(files, [{"analysis_id": "main", "test": {"kind": "welch"}, "correction": "none",
                          "bh_population": "annotated"}], injection_order={"source": "metadata"})
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["status"] == "success", out
    assert any(n.startswith("analyses[0].bh_population") for n in out["notes"])
    assert any(n.startswith("injection_order") for n in out["notes"])


# ---------- 外部の入力（MAF・mwTab・xlsx） ----------

def test_maf_input_reads_like_the_tsv_and_reports_unreadable_cells(files):
    maf_rows = [(f, [*v[:2], "ND" if f == "F2" else v[2], *v[3:]]) for f, v in VALUES.items()]
    maf = write_maf(files["dir"] / "m_paper.tsv", SAMPLES, maf_rows)
    result = _run(_spec(files, input={"source": "maf", "path": str(maf)}))
    assert result["input_read"] == {"source": "maf", "n_features": 4, "n_samples": 8,
                                    "n_unreadable_cells": 1, "unreadable_examples": ["ND"]}
    rows = _rows(result)
    assert sorted(rows) == ["row_1:F1", "row_2:F2", "row_3:F3", "row_4:F4"]
    tsv_rows = _rows(_run(_spec(files)))
    assert rows["row_1:F1"]["p"] == pytest.approx(tsv_rows["F1"]["p"])
    assert rows["row_2:F2"]["n_ref"] == 2                    # ND の 1 セルは欠測
    assert len(result["input_fingerprint"]["files"]["matrix"]) == 64
    # TSV の入力には input_read を作らない
    assert "input_read" not in _run(_spec(files))


def test_mwtab_and_xlsx_inputs(files):
    mwtab = write_mwtab(files["dir"] / "mw.txt", SAMPLES, list(VALUES.items()))
    via_mwtab = _rows(_run(_spec(files, input={"source": "mwtab", "path": str(mwtab)})))
    sheet = [["Supplementary Table"], ["Metabolite", *SAMPLES],
             *[[f, *v] for f, v in VALUES.items()]]
    xlsx = write_xlsx(files["dir"] / "supp.xlsx", {"S1": sheet})
    result = _run(_spec(files, input={"source": "xlsx", "path": str(xlsx)}))
    via_xlsx = _rows(result)
    assert result["spec"]["input"] == {"source": "xlsx", "path": str(xlsx), "sheet": None,
                                       "feature_column": "A", "first_sample_column": "B"}
    tsv = _rows(_run(_spec(files)))
    for fid in VALUES:
        assert via_mwtab[fid]["p"] == pytest.approx(tsv[fid]["p"], nan_ok=True)
        assert via_xlsx[fid]["p"] == pytest.approx(tsv[fid]["p"], nan_ok=True)


def test_external_input_rule_hash_ignores_the_path(files, tmp_path):
    a = write_mwtab(files["dir"] / "a.txt", SAMPLES, list(VALUES.items()))
    b = write_mwtab(tmp_path / "b.txt", SAMPLES, list(VALUES.items()))
    ha = A.rule_hash(A.normalize_spec(_spec(files, input={"source": "mwtab", "path": str(a)})))
    hb = A.rule_hash(A.normalize_spec(_spec(files, input={"source": "mwtab", "path": str(b)})))
    assert ha == hb


@pytest.mark.parametrize("raw, field", [
    ({"source": "maf"}, "input"),
    ({"source": "xlsx", "path": "x.xlsx", "feature_column": "C", "first_sample_column": "B"}, "input"),
    ({"source": "xlsx", "path": "x.xlsx", "feature_column": "1"}, "input.feature_column"),
    ({"source": "csv", "path": "x.csv"}, "input.source"),
])
def test_external_input_spec_errors(files, raw, field):
    with pytest.raises(DomainError) as info:
        A.normalize_spec(_spec(files, input=raw))
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == field


def test_external_input_needs_groups_tsv(files):
    spec = _spec(files, input={"source": "maf", "path": "m.tsv"})
    del spec["groups"]
    with pytest.raises(DomainError) as info:
        A.normalize_spec(spec)
    assert "input.source=maf" in info.value.message


# ---------- correlate ----------

def _corr(files, **extra):
    return {"external": {"source": "mwtab", "path": str(files["mwtab"])},
            "pairs": [{"feature_id": "F1", "external_id": "Glucose"},
                      {"feature_id": "F4", "external_id": "Lactate"}], **extra}


def test_correlate_input_stage_is_spearman_against_the_external_values(files):
    result = _run(_spec(files, correlate=_corr(files)))
    block = result["correlate"]
    assert block["stage"] == "input" and block["n_samples_matched"] == 8
    assert block["unmatched_samples"] == [] and block["unmatched_external_samples"] == []
    assert block["external"] == {"source": "mwtab", "n_features": 2, "n_samples": 8,
                                 "n_unreadable_cells": 0, "unreadable_examples": []}
    first, second = block["rows"]
    expected = stats.spearmanr(VALUES["F1"], GLUCOSE)
    assert first["analysis_id"] is None and first["merged_into"] is None
    assert (first["feature_id"], first["external_id"], first["n"]) == ("F1", "Glucose", 8)
    assert first["rho"] == pytest.approx(expected.statistic)
    # 8 組・同順位なし: 全順列を数えた正確な p（t 近似の scipy.stats.spearmanr ではない）
    assert first["p_method"] == "exact_permutation"
    assert first["p"] == pytest.approx(_exact_spearman_p(VALUES["F1"], GLUCOSE), rel=1e-12)
    assert second["rho"] == pytest.approx(stats.spearmanr(VALUES["F4"], LACTATE).statistic)
    assert len(result["input_fingerprint"]["files"]["correlate.external"]) == 64
    # 報告の数字は変えない・rule_hash は correlate の有無で分かれる
    plain = _run(_spec(files))
    assert result["analyses"] == plain["analyses"]
    assert result["rule_hash"] != plain["rule_hash"]
    assert "correlate" not in plain


def test_correlate_sample_map_and_unmatched(files):
    external = [f"P{i}" for i in range(1, 10)]
    mwtab = write_mwtab(files["dir"] / "renamed.txt", external,
                        [("Glucose", [*GLUCOSE, 99]), ("Lactate", [*LACTATE, 99])])
    sample_map = _write(files["dir"] / "map.tsv", ["sample", "external_sample"],
                        [[s, f"P{i}"] for i, s in enumerate(SAMPLES[:7], start=1)])
    spec = _spec(files, correlate={"external": {"source": "mwtab", "path": str(mwtab)},
                                   "pairs": [{"feature_id": "F1", "external_id": "Glucose"}],
                                   "sample_map": {"path": str(sample_map)}})
    block = _run(spec)["correlate"]
    assert block["n_samples_matched"] == 7
    assert block["unmatched_samples"] == ["QC2"]            # 対応表に無い試料は使わない
    assert block["unmatched_external_samples"] == ["P8", "P9"]
    [row] = block["rows"]
    assert row["n"] == 7
    assert row["rho"] == pytest.approx(stats.spearmanr(VALUES["F1"][:7], GLUCOSE[:7]).statistic)


def test_correlate_matrix_stage_reads_each_analysis_matrix(files):
    analyses = [{"analysis_id": "raw", "test": {"kind": "student"}},
                {"analysis_id": "filtered", "test": {"kind": "student"},
                 "matrix": {"steps": ["detection_filter"],
                            "detection_filter": {"basis": "positive", "min_rate": {"sample": 0.5}}}}]
    corr = _corr(files, stage="matrix",
                 pairs=[{"feature_id": "F1", "external_id": "Glucose"},
                        {"feature_id": "F3", "external_id": "Lactate"}])
    block = _run(_spec(files, analyses, correlate=corr))["correlate"]
    assert block["stage"] == "matrix"
    by_key = {(r["analysis_id"], r["feature_id"]): r for r in block["rows"]}
    assert set(by_key) == {("raw", "F1"), ("raw", "F3"), ("filtered", "F1"), ("filtered", "F3")}
    assert by_key[("raw", "F1")]["rho"] == pytest.approx(
        stats.spearmanr(VALUES["F1"], GLUCOSE).statistic)
    # recipe が落とした特徴（F3）は not_in_matrix
    assert by_key[("filtered", "F3")] == {
        "analysis_id": "filtered", "feature_id": "F3", "external_id": "Lactate",
        "merged_into": None, "n": 0, "rho": None, "p": None, "reason": "not_in_matrix", "p_method": None}


def test_correlate_errors(files):
    unknown_ext = _corr(files, pairs=[{"feature_id": "F1", "external_id": "Citrate"}])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=unknown_ext))
    assert info.value.code == "FEATURE_NOT_FOUND" and info.value.details["feature_ids"] == ["Citrate"]
    unknown = _corr(files, pairs=[{"feature_id": "F9", "external_id": "Glucose"}])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=unknown))
    assert info.value.code == "FEATURE_NOT_FOUND" and info.value.details["feature_ids"] == ["F9"]
    other = write_mwtab(files["dir"] / "other.txt", ["X1", "X2"], [("Glucose", [1, 2])])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate={"external": {"source": "mwtab", "path": str(other)},
                                     "pairs": [{"feature_id": "F1", "external_id": "Glucose"}]}))
    assert info.value.code == "SAMPLE_MISMATCH"


def test_correlate_sample_map_many_to_one_is_refused(files):
    # 外部の同じ試料に 2 つの試料を当てると、外部の 1 つの値を 2 回数えて n と p が膨らむ
    sample_map = _write(files["dir"] / "dup.tsv", ["sample", "external_sample"],
                        [["HC1", "HC1"], ["HC2", "HC1"], *[[s, s] for s in SAMPLES[2:]]])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=_corr(files, sample_map={"path": str(sample_map)})))
    assert info.value.code == "SAMPLE_MISMATCH"
    assert info.value.details["field"] == "correlate.sample_map"
    assert info.value.details["duplicated_external"] == ["HC1"]
    assert info.value.details["samples"] == {"HC1": ["HC1", "HC2"]}
    # この側の試料の重複も従来どおり SAMPLE_MISMATCH
    twice = _write(files["dir"] / "twice.tsv", ["sample", "external_sample"],
                   [["HC1", "HC1"], ["HC1", "HC2"]])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=_corr(files, sample_map={"path": str(twice)})))
    assert info.value.code == "SAMPLE_MISMATCH" and info.value.details["duplicated"] == ["HC1"]


def test_correlate_xlsx_external_keeps_undetected_rows_and_names_skipped_rows(files):
    sheet = [["Metabolite", *SAMPLES],
             ["Organic acids", *(["group"] * len(SAMPLES))],     # 群の見出しの行（飛ばす）
             ["Glucose", *GLUCOSE],
             ["Lactate", *(["ND"] * (len(SAMPLES) - 1)), "<LOQ"]]   # どの試料でも検出されなかった物質
    xlsx = write_xlsx(files["dir"] / "supp.xlsx", {"T": sheet})
    external = {"source": "xlsx", "path": str(xlsx)}
    block = _run(_spec(files, correlate={
        "external": external,
        "pairs": [{"feature_id": "F1", "external_id": "Glucose"},
                  {"feature_id": "F4", "external_id": "Lactate"}]}))["correlate"]
    assert block["external"] == {"source": "xlsx", "n_features": 2, "n_samples": 8,
                                 "n_unreadable_cells": 8, "unreadable_examples": ["ND", "<LOQ"],
                                 "n_rows_skipped_as_text": 1,
                                 "rows_skipped_as_text": ["Organic acids"]}
    lactate = block["rows"][1]
    assert lactate["external_id"] == "Lactate" and lactate["reason"] == "too_few_pairs"
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate={"external": external,
                                     "pairs": [{"feature_id": "F1", "external_id": "Organic acids"}]}))
    assert info.value.code == "FEATURE_NOT_FOUND"
    assert info.value.details["skipped_as_text"] == ["Organic acids"]
    assert "群の見出し" in info.value.message


def test_xlsx_input_feature_ids_name_skipped_rows(files):
    sheet = [["Metabolite", *SAMPLES], ["Lipids", *(["x"] * len(SAMPLES))],
             *[[f, *v] for f, v in VALUES.items()]]
    xlsx = write_xlsx(files["dir"] / "supp.xlsx", {"T": sheet})
    spec = _spec(files, input={"source": "xlsx", "path": str(xlsx)},
                 features={"feature_ids": ["F1", "Lipids"]})
    with pytest.raises(DomainError) as info:
        _run(spec)
    assert info.value.code == "FEATURE_NOT_FOUND"
    assert info.value.details["feature_ids"] == ["Lipids"]
    assert info.value.details["skipped_as_text"] == ["Lipids"]
    result = _run(_spec(files, input={"source": "xlsx", "path": str(xlsx)}))
    assert result["input_read"]["n_rows_skipped_as_text"] == 1


def test_correlate_too_few_matched_samples_is_null(files):
    sample_map = _write(files["dir"] / "few.tsv", ["sample", "external_sample"],
                        [[s, s] for s in SAMPLES[:4]])
    block = _run(_spec(files, correlate=_corr(files, sample_map={"path": str(sample_map)})))["correlate"]
    assert block["n_samples_matched"] == 4
    assert all(r["rho"] is None and r["reason"] == "too_few_pairs" for r in block["rows"])


@pytest.mark.parametrize("raw, field", [
    ({"external": {"source": "mwtab", "path": "x"}, "pairs": []}, "correlate.pairs"),
    ({"external": {"source": "mwtab", "path": "x"}, "pairs": [{"feature_id": "F1"}]},
     "correlate.pairs[0]"),
    ({"external": {"source": "mwtab", "path": "x"},
      "pairs": [{"feature_id": "F1", "external_id": "G"}], "stage": "after"}, "correlate.stage"),
    ({"external": {"source": "dataset"}, "pairs": [{"feature_id": "F1", "external_id": "G"}]},
     "correlate.external.source"),
])
def test_correlate_spec_errors(files, raw, field):
    with pytest.raises(DomainError) as info:
        A.normalize_spec(_spec(files, correlate=raw))
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == field
