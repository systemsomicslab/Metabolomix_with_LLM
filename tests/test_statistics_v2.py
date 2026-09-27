"""v2統計: 変換・効果量・検定（spec §10）。

v1 の two_group_test は `log2(x + pseudo_count)` で変換し、効果量も log 空間の
平均差から作る。それは v1 の契約として**そのまま変えない**が、v2 は別の定義を
持つ:

- `log2` は有限の正値だけに適用する。0以下を1へ clip したり pseudocount を足して
  「変換できた」ことにしない——clip は 0 と 1 を同じ値に潰し、pseudocount は
  比の大きさを定数の大きさにすり替える。変換できない値は欠損にして理由を残す。
- 効果量は**統計変換前**の test/reference の算術平均比の log2
  （`effect_size_definition = log2_arithmetic_mean_ratio`）。log 空間の平均差は
  幾何平均比であって、算術平均比ではない。
- BH の母集団は、同じ statistic_id で**検定できた** feature 全体。検定不能を
  母集団へ入れると q が不当に緩み、除いた集合を書かないと再現できない。
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from lipidmix.analysis.statistics_v2 import (
    arithmetic_log2fc,
    run_statistic,
    transform_values,
)
from lipidmix.core.atomic_io import DomainError


# ---------- brief記載のRED ----------

def test_v2_transform_and_effect_size():
    assert transform_values(np.array([.25, .5]), "log2").tolist() == [-2., -1.]
    assert math.isclose(arithmetic_log2fc(np.array([1., 9.]), np.array([4., 4.])),
                        math.log2(4 / 5))


# ---------- transform_values ----------

def test_log2_does_not_clip_non_positive_values():
    out = transform_values(np.array([0., -1., 4.]), "log2")
    assert np.isnan(out[0]) and np.isnan(out[1])
    assert out[2] == 2.0


def test_log2_adds_no_pseudocount():
    out = transform_values(np.array([1.]), "log2")
    assert out[0] == 0.0          # log2(1+1)=1 になっていたら pseudocount がある


def test_none_transform_is_identity():
    values = np.array([0., -1., 4.])
    out = transform_values(values, "none")
    assert np.array_equal(out, values, equal_nan=True)
    assert out is not values      # 元配列を書き換えない


def test_unknown_transform_is_rejected():
    with pytest.raises(DomainError) as caught:
        transform_values(np.array([1.]), "ln")
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


# ---------- arithmetic_log2fc ----------

def test_effect_size_is_the_arithmetic_not_the_geometric_ratio():
    reference = np.array([1., 9.])      # 算術平均 5、幾何平均 3
    test = np.array([4., 4.])
    assert math.isclose(arithmetic_log2fc(reference, test), math.log2(4 / 5))
    assert not math.isclose(arithmetic_log2fc(reference, test), math.log2(4 / 3))


def test_non_positive_mean_gives_na_not_a_clipped_number():
    assert math.isnan(arithmetic_log2fc(np.array([-1., 1.]), np.array([4., 4.])))
    assert math.isnan(arithmetic_log2fc(np.array([1., 1.]), np.array([0., 0.])))


def test_empty_group_gives_na():
    assert math.isnan(arithmetic_log2fc(np.array([np.nan]), np.array([4.])))


# ---------- run_statistic の土台 ----------

def _rows(groups: list[str], *, roles=None, bio=None, include=None) -> list[dict]:
    rows = []
    for index, group in enumerate(groups):
        rows.append({
            "sample_id": f"s{index}", "source_file": f"S{index}.wiff",
            "role": (roles or {}).get(index, "sample"), "group": group,
            "batch": "B1", "injection_order": index + 1, "qc_pool": None,
            "include": (include or {}).get(index, True),
            "biological_sample_id": (bio or {}).get(index, f"bio{index}"),
        })
    return rows


def _matrix(values, *, feature_ids=None, eligibility=None, units=None) -> dict:
    values = np.asarray(values, dtype=float)
    n_assays, n_features = values.shape
    feature_ids = feature_ids or [str(i + 1) for i in range(n_features)]
    return {
        "schema": "analysis-matrix.v1", "matrix_id": "mat_test",
        "parent_id": None, "stage": "finalized",
        "recipe_id": "default",
        "assay_ids": [f"assay[{i + 1}]" for i in range(n_assays)],
        "feature_ids": feature_ids,
        "values": values,
        "units": units or ["peak_height"] * n_features,
        "eligibility_mask": (np.ones(n_features, dtype=bool) if eligibility is None
                             else np.asarray(eligibility, dtype=bool)),
        "detected_mask": None,
        "imputed_mask": np.zeros(values.shape, dtype=bool),
        "locked_mask": np.zeros(values.shape, dtype=bool),
        "missing_reasons": {}, "support_feature_ids": [],
        "correction_history": [], "caveats": [],
    }


def _welch_spec(**overrides) -> dict:
    spec = {"statistic_id": "w1", "kind": "welch", "matrix_recipe_id": "default",
            "transform": "none", "feature_scope": {"mode": "all_eligible"},
            "reference_group": "control", "test_group": "treated",
            "q_threshold": 0.05, "log2fc_threshold": 1.0}
    spec.update(overrides)
    return spec


def _two_group_matrix():
    """control 3 / treated 3。feature1 は明確に上昇、feature2 は定数。"""
    values = [[10., 5.], [11., 5.], [9., 5.],
              [40., 5.], [42., 5.], [38., 5.]]
    metadata = _rows(["control"] * 3 + ["treated"] * 3)
    return _matrix(values), metadata


# ---------- welch ----------

def test_welch_reports_p_q_and_effect_size():
    matrix, metadata = _two_group_matrix()
    out = run_statistic(matrix, _welch_spec(), metadata)

    assert out["status"] == "completed"
    assert out["effect_size_definition"] == "log2_arithmetic_mean_ratio"
    by_id = {f["feature_id"]: f for f in out["features"]}
    assert by_id["1"]["p_value"] < 0.05
    assert by_id["1"]["q_value"] is not None
    assert by_id["1"]["log2fc"] == pytest.approx(math.log2(40 / 10), abs=1e-6)


def test_welch_t_statistic_has_the_same_sign_as_log2fc():
    # t は「test − reference」。log2FC（test / reference）や Tukey の
    # mean_difference（test − reference）と同じ向きでないと、同じ TSV の同じ行で
    # 逆符号の 2 列が並ぶ（spec A14「群方向が一致する」）。
    matrix, metadata = _two_group_matrix()
    out = run_statistic(matrix, _welch_spec(), metadata)

    feature = {f["feature_id"]: f for f in out["features"]}["1"]
    assert feature["log2fc"] > 0
    assert feature["t_statistic"] > 0


def test_effect_size_uses_untransformed_values_even_with_log2():
    matrix, metadata = _two_group_matrix()
    plain = run_statistic(matrix, _welch_spec(), metadata)
    logged = run_statistic(matrix, _welch_spec(transform="log2"), metadata)

    plain_fc = {f["feature_id"]: f["log2fc"] for f in plain["features"]}
    logged_fc = {f["feature_id"]: f["log2fc"] for f in logged["features"]}
    assert logged_fc["1"] == pytest.approx(plain_fc["1"])
    assert logged["transform"] == "log2"
    # 変換は検定側にだけ効く（p値は変わってよい）。
    logged_p = {f["feature_id"]: f["p_value"] for f in logged["features"]}
    plain_p = {f["feature_id"]: f["p_value"] for f in plain["features"]}
    assert logged_p["1"] != plain_p["1"]


def test_constant_feature_is_not_tested_and_says_why():
    matrix, metadata = _two_group_matrix()
    out = run_statistic(matrix, _welch_spec(), metadata)
    feature = next(f for f in out["features"] if f["feature_id"] == "2")
    assert feature["status"] == "not_testable"
    assert feature["p_value"] is None
    assert feature["q_value"] is None
    assert feature["reason"] == "zero_variance"


def test_bh_population_is_only_the_testable_features():
    """検定不能featureをBHの母集団に入れない（入れるとqが不当に緩む）。"""
    values = [[10., 5.], [11., 5.], [9., 5.],
              [40., 5.], [42., 5.], [38., 5.]]
    matrix = _matrix(values)
    metadata = _rows(["control"] * 3 + ["treated"] * 3)
    out = run_statistic(matrix, _welch_spec(), metadata)

    tested = [f for f in out["features"] if f["status"] == "tested"]
    assert out["bh_population"] == len(tested) == 1
    # 母集団1件なので q == p。母集団に定数featureを混ぜていれば q = 2p になる。
    assert tested[0]["q_value"] == pytest.approx(tested[0]["p_value"])


def test_ineligible_features_are_excluded_from_the_test():
    matrix, metadata = _two_group_matrix()
    matrix["eligibility_mask"] = np.array([False, True])
    out = run_statistic(matrix, _welch_spec(), metadata)
    assert [f["feature_id"] for f in out["features"]] == ["2"]


def test_feature_scope_targets_restricts_the_population():
    matrix, metadata = _two_group_matrix()
    spec = _welch_spec(feature_scope={"mode": "targets", "target_ids": ["gaba"]})
    out = run_statistic(matrix, spec, metadata,
                        target_features={"gaba": "1"})
    assert [f["feature_id"] for f in out["features"]] == ["1"]


def test_all_features_untestable_makes_the_statistic_not_evaluable():
    values = [[5.], [5.], [5.], [5.], [5.], [5.]]
    matrix = _matrix(values)
    metadata = _rows(["control"] * 3 + ["treated"] * 3)
    out = run_statistic(matrix, _welch_spec(), metadata)
    assert out["status"] == "not_evaluable"
    assert out["reason"] == "no_testable_feature"


# ---------- 母集団の検査 ----------

def test_a_group_with_one_sample_stops_the_statistic():
    values = [[10.], [11.], [40.]]
    matrix = _matrix(values)
    metadata = _rows(["control", "control", "treated"])
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _welch_spec(), metadata)
    assert caught.value.code == "STATISTIC_GROUP_TOO_SMALL"


def test_repeated_injections_of_one_donor_stop_the_statistic():
    matrix, metadata = _two_group_matrix()
    metadata[0]["biological_sample_id"] = "bio1"      # s0 と s1 が同一個体
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _welch_spec(), metadata)
    assert caught.value.code == "REPEATED_MEASURES_UNSUPPORTED"


def test_qc_and_standard_injections_never_enter_a_group():
    values = [[10.], [11.], [9.], [40.], [42.], [38.], [999.], [999.]]
    matrix = _matrix(values)
    metadata = _rows(["control"] * 3 + ["treated"] * 3 + ["control", "treated"],
                     roles={6: "qc", 7: "standard"})
    out = run_statistic(matrix, _welch_spec(), metadata)
    assert out["groups"] == {"control": 3, "treated": 3}


def test_excluded_injections_never_enter_a_group():
    values = [[10.], [11.], [9.], [40.], [42.], [38.], [999.]]
    matrix = _matrix(values)
    metadata = _rows(["control"] * 3 + ["treated"] * 3 + ["treated"],
                     include={6: False})
    out = run_statistic(matrix, _welch_spec(), metadata)
    assert out["groups"]["treated"] == 3


def test_same_group_twice_is_rejected():
    matrix, metadata = _two_group_matrix()
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _welch_spec(test_group="control"), metadata)
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


# ---------- anova_tukey ----------

def _three_group_matrix():
    values = [[1.], [2.], [3.], [2.], [3.], [4.], [4.], [5.], [6.]]
    metadata = _rows(["g1"] * 3 + ["g2"] * 3 + ["g3"] * 3)
    return _matrix(values), metadata


def _anova_spec(**overrides) -> dict:
    spec = {"statistic_id": "a1", "kind": "anova_tukey",
            "matrix_recipe_id": "default", "transform": "none",
            "feature_scope": {"mode": "all_eligible"},
            "groups": ["g1", "g2", "g3"], "alpha": 0.05}
    spec.update(overrides)
    return spec


def test_anova_result_carries_f_df_p_q_and_group_n():
    matrix, metadata = _three_group_matrix()
    out = run_statistic(matrix, _anova_spec(), metadata)
    feature = out["features"][0]
    assert feature["f_statistic"] == pytest.approx(7.0)
    assert (feature["df_between"], feature["df_within"]) == (2, 6)
    assert feature["q_value"] is not None
    assert feature["group_n"] == {"g1": 3, "g2": 3, "g3": 3}


def test_tukey_is_run_for_every_feature_not_just_significant_anova():
    values = [[1., 5.], [2., 5.1], [3., 4.9], [2., 5.], [3., 5.1], [4., 4.9],
              [4., 5.], [5., 5.1], [6., 4.9]]
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 3 + ["g2"] * 3 + ["g3"] * 3)
    out = run_statistic(matrix, _anova_spec(), metadata)
    for feature in out["features"]:
        assert len(feature["tukey"]) == 3


def test_tukey_is_not_reported_as_a_cross_feature_fdr():
    matrix, metadata = _three_group_matrix()
    out = run_statistic(matrix, _anova_spec(), metadata)
    assert out["tukey_correction_scope"] == "within_feature_group_pairs"


def test_fewer_than_three_groups_is_rejected_for_anova():
    matrix, metadata = _three_group_matrix()
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _anova_spec(groups=["g1", "g2"]), metadata)
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


# ---------- pca ----------

def _pca_spec(**overrides) -> dict:
    spec = {"statistic_id": "p1", "kind": "pca", "matrix_recipe_id": "default",
            "transform": "none", "feature_scope": {"mode": "all_eligible"},
            "scaling": "autoscale", "n_components": 5}
    spec.update(overrides)
    return spec


def test_pca_limits_components_to_the_matrix_rank():
    values = np.array([[1., 2., 3.], [2., 4., 5.], [3., 5., 9.], [7., 1., 2.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 2 + ["g2"] * 2)
    out = run_statistic(matrix, _pca_spec(), metadata)
    assert out["status"] == "completed"
    assert out["n_components"] <= min(values.shape)
    assert len(out["explained_variance_ratio"]) == out["n_components"]


def test_pca_excludes_non_finite_features_and_counts_them():
    values = np.array([[1., 2., np.nan], [2., 4., 1.], [3., 5., 2.], [7., 1., 3.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 2 + ["g2"] * 2)
    out = run_statistic(matrix, _pca_spec(), metadata)
    assert out["n_features_excluded"] == 1
    assert out["n_features_used"] == 2


def test_pca_does_not_learn_from_qc_or_standards():
    values = np.array([[1., 2.], [2., 4.], [3., 5.], [7., 1.], [999., 999.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 2 + ["g2"] * 2 + ["g2"], roles={4: "qc"})
    out = run_statistic(matrix, _pca_spec(), metadata)
    assert out["n_samples"] == 4


def test_pca_with_too_few_samples_is_not_evaluable():
    values = np.array([[1., 2.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"])
    out = run_statistic(matrix, _pca_spec(), metadata)
    assert out["status"] == "not_evaluable"
    assert out["reason"] == "insufficient_samples_or_features"


@pytest.mark.parametrize("scaling", ["none", "autoscale"])
def test_pca_components_are_limited_by_centered_rank(scaling):
    out = run_statistic(_matrix(np.array([[1., 2.], [2., 3.]])),
                        _pca_spec(scaling=scaling), _rows(["g1", "g2"]))
    assert out["status"] == "completed"
    assert out["n_components"] == 1
    assert out["explained_variance_ratio"] == [1.0]
    assert np.asarray(out["scores"]).shape == (2, 1)


# ---------- 共通の記録 ----------

def test_every_result_records_the_matrix_it_used():
    matrix, metadata = _two_group_matrix()
    out = run_statistic(matrix, _welch_spec(), metadata)
    assert out["matrix_id"] == "mat_test"
    assert out["matrix_recipe_id"] == "default"
    assert out["statistic_id"] == "w1"


def test_unknown_kind_is_rejected():
    matrix, metadata = _two_group_matrix()
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _welch_spec(kind="t_test"), metadata)
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_pca_scaling_none_is_not_silently_autoscaled():
    """profileが scaling="none" と書いたなら、分散で割らない。"""
    values = np.array([[1., 100.], [2., 400.], [3., 500.], [7., 120.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 2 + ["g2"] * 2)
    auto = run_statistic(matrix, _pca_spec(scaling="autoscale"), metadata)
    plain = run_statistic(matrix, _pca_spec(scaling="none"), metadata)
    assert auto["scaling"] == "autoscale" and plain["scaling"] == "none"
    assert auto["explained_variance_ratio"] != plain["explained_variance_ratio"]


def test_unknown_pca_scaling_is_rejected():
    values = np.array([[1., 2.], [2., 4.], [3., 5.], [7., 1.]])
    matrix = _matrix(values)
    metadata = _rows(["g1"] * 2 + ["g2"] * 2)
    with pytest.raises(DomainError) as caught:
        run_statistic(matrix, _pca_spec(scaling="pareto"), metadata)
    assert caught.value.code == "STATISTIC_SPECIFICATION_INVALID"

