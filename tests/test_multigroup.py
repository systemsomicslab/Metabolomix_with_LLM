"""3群以上の比較（spec §10）: 一元配置ANOVA と Tukey HSD。

ANOVAは「どこかに差がある」しか言わない。どの対かはTukeyが答えるが、両者は
**同じ変換後行列**で計算し、ANOVAが有意なfeatureだけにTukeyを絞らない
（絞ると、Tukeyの補正がANOVAの選択に条件付けられた別物になる）。

期待値は独立に作る。F と df は手計算、Tukey の信頼区間は**公表された
studentized range の臨界値**から手計算した値と照合する——計算に使うのと同じ
ライブラリで期待値を作ると、実装が間違っていても一致してしまう。
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from metabolomix.analysis.multigroup import test_feature as run_feature_test
from metabolomix.core.atomic_io import DomainError

#: 固定小例 [1,2,3] / [2,3,4] / [4,5,6] の独立な期待値。
#:
#: ANOVA（手計算）: 群平均 2, 3, 5、全体平均 10/3。
#:   SSB = 3×((2−10/3)² + (3−10/3)² + (5−10/3)²) = 14、df=2 → MSB = 7
#:   SSW = 2 + 2 + 2 = 6、df=6 → MSW = 1     ⇒ F = 7/1 = 7
#:
#: Tukey（公表の臨界値から手計算）: q_{0.05}(k=3, df=6) = 4.339
#:   （studentized range の標準的な臨界値表。scipy では生成していない）
#:   半幅 = q × √(MSW/n) = 4.339 × √(1/3) = 2.5052
#:   差は test − reference。(g1→g3) の差 = 5 − 2 = 3 ⇒ CI = [0.4948, 5.5052]
_EXPECTED = {
    "f_statistic": 7.0,
    "df_between": 2,
    "df_within": 6,
    "tukey_half_width": 2.5052,
    "pairs": {
        ("g1", "g2"): {"difference": 1.0, "significant": False},
        ("g1", "g3"): {"difference": 3.0, "significant": True},
        ("g2", "g3"): {"difference": 2.0, "significant": False},
    },
}


def _groups():
    return [np.array([1., 2., 3.]), np.array([2., 3., 4.]), np.array([4., 5., 6.])]


# ---------- ANOVA ----------

def test_anova_matches_the_hand_computed_example():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    assert out["status"] == "tested"
    assert out["f_statistic"] == pytest.approx(_EXPECTED["f_statistic"])
    assert out["df_between"] == _EXPECTED["df_between"]
    assert out["df_within"] == _EXPECTED["df_within"]
    assert 0.0 < out["p_value"] < 1.0


def test_group_sizes_are_reported():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    assert out["group_n"] == {"g1": 3, "g2": 3, "g3": 3}


# ---------- Tukey ----------

def test_tukey_differences_are_test_minus_reference():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    pairs = {(p["reference_group"], p["test_group"]): p for p in out["tukey"]}
    for key, expected in _EXPECTED["pairs"].items():
        assert pairs[key]["mean_difference"] == pytest.approx(expected["difference"])


def test_tukey_confidence_intervals_match_the_published_critical_value():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    pairs = {(p["reference_group"], p["test_group"]): p for p in out["tukey"]}
    pair = pairs[("g1", "g3")]
    half = _EXPECTED["tukey_half_width"]
    assert pair["ci_low"] == pytest.approx(3.0 - half, abs=1e-3)
    assert pair["ci_high"] == pytest.approx(3.0 + half, abs=1e-3)


def test_only_the_separated_pair_is_significant():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    pairs = {(p["reference_group"], p["test_group"]): p for p in out["tukey"]}
    for key, expected in _EXPECTED["pairs"].items():
        assert (pairs[key]["p_adjusted"] < 0.05) is expected["significant"]
        # CI が0をまたぐかどうかと、補正pの判定は一致していなければならない。
        excludes_zero = pairs[key]["ci_low"] > 0 or pairs[key]["ci_high"] < 0
        assert excludes_zero is expected["significant"]


def test_every_pair_is_evaluated_in_one_call():
    """ANOVA有意なfeatureだけにTukeyを絞らない。全群対を一度に見る。"""
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    assert len(out["tukey"]) == 3
    assert out["alpha"] == 0.05


def test_alpha_widens_the_interval():
    narrow = run_feature_test(_groups(), 0.01, labels=["g1", "g2", "g3"])
    wide = run_feature_test(_groups(), 0.10, labels=["g1", "g2", "g3"])
    pick = lambda out: next(p for p in out["tukey"]                 # noqa: E731
                            if (p["reference_group"], p["test_group"]) == ("g1", "g3"))
    assert pick(narrow)["ci_low"] < pick(wide)["ci_low"]


# ---------- 検定できない場合 ----------

def test_fewer_than_three_groups_is_rejected():
    with pytest.raises(DomainError) as caught:
        run_feature_test([np.array([1., 2.]), np.array([3., 4.])], 0.05,
                     labels=["a", "b"])
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_a_group_with_one_finite_value_is_not_testable():
    groups = [np.array([1., np.nan]), np.array([2., 3.]), np.array([4., 5.])]
    out = run_feature_test(groups, 0.05, labels=["g1", "g2", "g3"])
    assert out["status"] == "not_testable"
    assert out["reason"] == "insufficient_finite_values"
    assert out["p_value"] is None
    assert out["tukey"] == []


def test_a_constant_feature_is_not_testable():
    groups = [np.array([5., 5.]), np.array([5., 5.]), np.array([5., 5.])]
    out = run_feature_test(groups, 0.05, labels=["g1", "g2", "g3"])
    assert out["status"] == "not_testable"
    assert out["reason"] == "zero_within_group_variance"
    assert out["p_value"] is None


def test_a_nan_p_value_is_never_reported_as_significant():
    groups = [np.array([5., 5.]), np.array([5., 5.]), np.array([5., 5.])]
    out = run_feature_test(groups, 0.05, labels=["g1", "g2", "g3"])
    assert out["p_value"] is None
    assert not any(p["p_adjusted"] == 0 for p in out["tukey"])


def test_label_count_must_match_the_groups():
    with pytest.raises(DomainError) as caught:
        run_feature_test(_groups(), 0.05, labels=["g1", "g2"])
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_alpha_outside_the_unit_interval_is_rejected():
    with pytest.raises(DomainError) as caught:
        run_feature_test(_groups(), 1.5, labels=["g1", "g2", "g3"])
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_f_statistic_is_finite_when_reported():
    out = run_feature_test(_groups(), 0.05, labels=["g1", "g2", "g3"])
    assert math.isfinite(out["f_statistic"])
