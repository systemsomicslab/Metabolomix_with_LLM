"""3群以上の比較: 一元配置ANOVA と Tukey HSD（spec §10）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

**ANOVAとTukeyは同じ変換後行列で、同時に計算する。** ANOVAが有意なfeatureだけに
Tukeyを絞ると、Tukeyの補正が「ANOVAで選ばれた」という条件付きの分布になり、
報告した補正pが実際の誤り率と合わなくなる。ここでは全featureで全群対を一度に
評価する。

**Tukeyの補正はfeature内の群対に閉じている。** feature横断のFDRではない。両者を
同じ「補正済みp」として並べると、読み手はfeature数で補正されたと誤読する。
呼び出し側（`statistics_v2`）はこの区別を出力に明記する。

**SciPyのAPIが無ければ止める。** 別の検定へ黙って落ちない——「Tukeyのつもりが
Bonferroniだった」は結果を見ても分からない。
"""
from __future__ import annotations

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = ["test_feature"]

_SPEC_INVALID = "STATISTIC_SPECIFICATION_INVALID"

#: 各群に必要な有限値の数。1点の群は群内分散を持たない。
_MIN_PER_GROUP = 2


def _require_scipy():
    """ANOVA/Tukey に必要な SciPy API を取りに行く（無ければ明示的に失敗）。"""
    try:
        from scipy import stats
    except ImportError as exc:                    # pragma: no cover - 環境依存
        raise DomainError(
            "STATISTIC_DEPENDENCY_MISSING",
            "多群比較には SciPy が必要です（別の検定へは切り替えません）。",
            {"missing": "scipy"}) from exc
    missing = [name for name in ("f_oneway", "tukey_hsd")
               if not hasattr(stats, name)]
    if missing:                                   # pragma: no cover - 環境依存
        raise DomainError(
            "STATISTIC_DEPENDENCY_MISSING",
            f"SciPy に必要な関数がありません: {', '.join(missing)}"
            "（別の検定へは切り替えません）。",
            {"missing": missing})
    return stats


def _not_testable(reason: str, labels: list[str], groups) -> dict:
    return {
        "status": "not_testable", "reason": reason,
        "f_statistic": None, "p_value": None,
        "df_between": None, "df_within": None,
        "group_n": {label: int(np.isfinite(g).sum())
                    for label, g in zip(labels, groups)},
        "tukey": [], "alpha": None,
    }


def test_feature(groups: list[np.ndarray], alpha: float,
                 labels: list[str] | None = None) -> dict:
    """1 feature ぶんの ANOVA と Tukey HSD を返す。

    `groups` は群ごとの値（統計変換後）。`labels` は同じ並びの群名で、
    Tukey の各対は `test_group - reference_group` の向きで報告する
    ——SciPy の `statistic[i][j]` は `mean_i - mean_j` なので、そのまま出すと
    「処置群が下がった」と読める符号が付く。

    検定できない feature は例外にせず `status="not_testable"` と理由を返す
    （1 feature の欠測で検定全体を落とさない）。NaN の p 値は `None` にして、
    0 や「有意」へ化けさせない。
    """
    labels = list(labels or [f"group{i + 1}" for i in range(len(groups))])
    if len(labels) != len(groups):
        raise DomainError(
            _SPEC_INVALID,
            f"群名の数が群の数と一致しません（labels={len(labels)}, "
            f"groups={len(groups)}）。",
            {"labels": len(labels), "groups": len(groups)})
    if len(groups) < 3:
        raise DomainError(
            _SPEC_INVALID,
            f"ANOVA/Tukey には3群以上が必要です（現在 {len(groups)} 群）。"
            "n不足の群を落として別のANOVAにはしません。",
            {"n_groups": len(groups)})
    if not (0.0 < float(alpha) < 1.0):
        raise DomainError(
            _SPEC_INVALID, f"alpha は 0 と 1 の間で指定してください: {alpha!r}",
            {"alpha": alpha})

    finite = [np.asarray(g, dtype=float)[np.isfinite(np.asarray(g, dtype=float))]
              for g in groups]
    if any(g.size < _MIN_PER_GROUP for g in finite):
        return _not_testable("insufficient_finite_values", labels, groups)
    if all(float(np.var(g, ddof=1)) == 0.0 for g in finite):
        # 全群で群内分散0。F は 0/0 で NaN になる。
        return _not_testable("zero_within_group_variance", labels, groups)

    stats = _require_scipy()
    anova = stats.f_oneway(*finite)
    if not np.isfinite(anova.statistic) or not np.isfinite(anova.pvalue):
        return _not_testable("test_statistic_not_finite", labels, groups)

    n_total = sum(g.size for g in finite)
    df_between = len(finite) - 1
    df_within = n_total - len(finite)

    result = stats.tukey_hsd(*finite)
    interval = result.confidence_interval(1.0 - float(alpha))
    tukey = []
    for i in range(len(finite)):
        for j in range(i + 1, len(finite)):
            # SciPy の [a][b] は mean_a - mean_b。reference=i, test=j の差
            # （test - reference）が欲しいので [j][i] を読む。
            tukey.append({
                "reference_group": labels[i],
                "test_group": labels[j],
                "mean_difference": float(result.statistic[j][i]),
                "ci_low": float(interval.low[j][i]),
                "ci_high": float(interval.high[j][i]),
                "p_adjusted": float(result.pvalue[j][i]),
            })

    return {
        "status": "tested", "reason": None,
        "f_statistic": float(anova.statistic),
        "p_value": float(anova.pvalue),
        "df_between": int(df_between),
        "df_within": int(df_within),
        "group_n": {label: int(g.size) for label, g in zip(labels, finite)},
        "tukey": tukey,
        "alpha": float(alpha),
    }


#: pytest は `test_` で始まる名前をテスト関数として収集する。これは plan が
#: 定めた公開APIの名前なので変えられないが、import しただけで「引数 groups の
#: fixture が無い」というエラーになるのを防ぐ。
test_feature.__test__ = False
