"""固定母集団のQCと、解析filterの分離（spec §9）。

QCと解析filterは**別の軸**。同じ mask を兼用すると、閾値で feature を落とすほど
バッチが合格に近づく——落とした分だけ fail が分母から消えるからで、これは
「品質が上がった」ではなく「測らなかった」。だからここで縛るのは:

1. 評価集合（母集団）は filter 前に固定し、以後の工程でも**同じ集合**を使う。
   処理後QCで集合が縮まない。
2. 除外は解析行列の eligibility にだけ反映し、QC結果を再集計しない。
3. 評価不能（U）は分母から消さない。`P/N≥t` で pass、`(P+U)/N<t` で fail、
   それ以外は not_evaluable。
4. 前提が無い metric は not_evaluable。標準品しか無いバッチに pooled QC の
   合格を付けない。補完後の行列で QC-RSD を合格にしない。
"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.analysis.assay_qc import aggregate_counts, evaluate_qc


# ---------- brief記載のRED ----------

def test_exclusion_does_not_turn_failed_batch_into_pass():
    assert aggregate_counts(60, 40, 0, .8) == "fail"
    assert aggregate_counts(60, 10, 30, .8) == "not_evaluable"
    assert aggregate_counts(0, 0, 0, .8) == "not_evaluable"


def test_pass_fraction_counts_unevaluable_in_the_denominator():
    assert aggregate_counts(80, 0, 20, .8) == "pass"
    assert aggregate_counts(79, 0, 21, .8) == "not_evaluable"
    assert aggregate_counts(79, 21, 0, .8) == "fail"


def test_blank_fold_cannot_borrow_blanks_from_another_batch():
    matrix = _matrix([[100.], [100.], [1.], [100.]])
    metadata = _rows([{"batch": "A"}, {"batch": "B"},
                      {"role": "blank", "batch": "A"},
                      {"role": "blank", "batch": "B"}])
    policy = [_policy("blank_fold", threshold={"operator": ">=", "value": 1.5})]
    out = evaluate_qc(matrix, _evidence(), metadata, policy, None)
    assert out["batch_status"] == "fail"
    assert sorted(e["value"] for e in out["metrics"]["blank_fold"]["elements"]) == [1., 100.]


@pytest.mark.parametrize("n_second_qc", [0, 2])
def test_qc_insufficient_batch_remains_in_evaluation(n_second_qc):
    metadata = _rows([{"role": "qc", "qc_pool": "pool", "batch": "A"}] * 3
                     + [{"role": "qc", "qc_pool": "pool", "batch": "B"}] * n_second_qc
                     + [{"batch": "B"}])
    matrix = _matrix([[100.]] * len(metadata))
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    assert out["batch_status"] == "not_evaluable"
    assert out["metrics"]["pooled_qc_rsd"]["counts"]["not_evaluable"] == 1


# ---------- fixture ----------

def _rows(specs: list[dict]) -> list[dict]:
    rows = []
    for index, spec in enumerate(specs):
        row = {"sample_id": f"s{index}", "source_file": f"S{index}.wiff",
               "role": "sample", "group": "A", "batch": "B1",
               "injection_order": index + 1, "qc_pool": None, "include": True,
               "biological_sample_id": f"bio{index}"}
        row.update(spec)
        rows.append(row)
    return rows


def _matrix(values, *, feature_ids=None, assay_ids=None, imputed=None,
            locked=None, detected=None, eligibility=None, units=None) -> dict:
    values = np.asarray(values, dtype=float)
    n_assays, n_features = values.shape
    feature_ids = feature_ids or [str(i + 1) for i in range(n_features)]
    assay_ids = assay_ids or [f"assay[{i + 1}]" for i in range(n_assays)]
    return {
        "schema": "analysis-matrix.v1", "matrix_id": "mat_test",
        "stage": "preprocessed", "parent_id": None,
        "assay_ids": assay_ids, "feature_ids": feature_ids,
        "values": values,
        "units": units or ["peak_height"] * n_features,
        "eligibility_mask": (np.ones(n_features, dtype=bool)
                             if eligibility is None
                             else np.asarray(eligibility, dtype=bool)),
        "detected_mask": None if detected is None else np.asarray(detected, dtype=bool),
        "imputed_mask": (np.zeros(values.shape, dtype=bool) if imputed is None
                         else np.asarray(imputed, dtype=bool)),
        "locked_mask": (np.zeros(values.shape, dtype=bool) if locked is None
                        else np.asarray(locked, dtype=bool)),
        "missing_reasons": {}, "support_feature_ids": [],
        "correction_history": [], "caveats": [],
    }


def _policy(metric, **overrides) -> dict:
    entry = {"metric": metric, "scope": "all_features", "target_ids": None,
             "required": True, "threshold": {"operator": "<=", "value": 30.0},
             "evidence_requirement": False, "minimum_pass_fraction": None}
    entry.update(overrides)
    return entry


def _evidence(rows=(), *, availability=True) -> dict:
    return {"schema": "assay-feature-evidence.v1", "availability": availability,
            "reasons": [], "rows": list(rows)}


# ---------- pooled_qc_rsd ----------

def _qc_matrix(qc_values) -> tuple[dict, list[dict]]:
    """QC3本 + 生物試料2本。QC値だけを引数で変える。"""
    values = [[v] for v in qc_values] + [[100.], [110.]]
    metadata = _rows([{"role": "qc", "qc_pool": "pool1"}] * len(qc_values)
                     + [{"role": "sample"}, {"role": "sample"}])
    return _matrix(values), metadata


def test_pooled_qc_rsd_passes_when_qc_is_tight():
    matrix, metadata = _qc_matrix([100., 101., 99.])
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    assert out["metrics"]["pooled_qc_rsd"]["status"] == "pass"


def test_pooled_qc_rsd_fails_when_qc_scatters():
    matrix, metadata = _qc_matrix([1., 10., 100.])
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    assert out["metrics"]["pooled_qc_rsd"]["status"] == "fail"
    assert out["batch_status"] == "fail"


def test_two_qc_injections_are_not_enough():
    matrix, metadata = _qc_matrix([100., 101.])
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    metric = out["metrics"]["pooled_qc_rsd"]
    assert metric["status"] == "not_evaluable"
    assert out["batch_status"] == "not_evaluable"


def test_standards_are_not_counted_as_pooled_qc():
    """標準品だけのバッチにpooled QCの合格を付けない。"""
    values = [[100.], [101.], [99.], [100.]]
    metadata = _rows([{"role": "standard"}] * 3 + [{"role": "sample"}])
    out = evaluate_qc(_matrix(values), _evidence(), metadata,
                      [_policy("pooled_qc_rsd")], None)
    assert out["metrics"]["pooled_qc_rsd"]["status"] == "not_evaluable"
    assert out["metrics"]["pooled_qc_rsd"]["reason"] == "insufficient_qc_injections"


def test_qc_without_a_declared_pool_is_not_evaluable():
    """pool不明はpooled QCの必須判定を満たさない（spec §7）。"""
    values = [[100.], [101.], [99.], [100.]]
    metadata = _rows([{"role": "qc"}] * 3 + [{"role": "sample"}])
    out = evaluate_qc(_matrix(values), _evidence(), metadata,
                      [_policy("pooled_qc_rsd")], None)
    assert out["metrics"]["pooled_qc_rsd"]["status"] == "not_evaluable"


def test_rsd_on_an_imputed_matrix_is_refused():
    """補完後の行列でQC-RSDを合格にしない（補完は分散を縮める）。"""
    matrix, metadata = _qc_matrix([1., 10., 100.])
    matrix["imputed_mask"] = np.ones(matrix["values"].shape, dtype=bool)
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    metric = out["metrics"]["pooled_qc_rsd"]
    assert metric["status"] == "not_evaluable"
    assert metric["reason"] == "imputed_matrix_not_allowed"


# ---------- blank_fold ----------

def _blank_matrix(sample_values, blank_values) -> tuple[dict, list[dict]]:
    values = [[v] for v in sample_values] + [[v] for v in blank_values]
    metadata = _rows([{"role": "sample"}] * len(sample_values)
                     + [{"role": "blank"}] * len(blank_values))
    return _matrix(values), metadata


def test_blank_fold_compares_medians():
    matrix, metadata = _blank_matrix([100., 120.], [10.])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("blank_fold",
                               threshold={"operator": ">=", "value": 3.0})], None)
    metric = out["metrics"]["blank_fold"]
    assert metric["status"] == "pass"
    assert metric["elements"][0]["value"] == pytest.approx(11.0)


def test_zero_blank_with_signal_is_positive_infinity_not_a_json_number():
    matrix, metadata = _blank_matrix([100., 120.], [0.])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("blank_fold",
                               threshold={"operator": ">=", "value": 3.0})], None)
    element = out["metrics"]["blank_fold"]["elements"][0]
    assert element["value"] is None
    assert element["special_value"] == "positive_infinity"
    assert element["status"] == "pass"


def test_zero_over_zero_is_not_evaluable():
    matrix, metadata = _blank_matrix([0., 0.], [0.])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("blank_fold",
                               threshold={"operator": ">=", "value": 3.0})], None)
    element = out["metrics"]["blank_fold"]["elements"][0]
    assert element["status"] == "not_evaluable"
    assert element["special_value"] is None


def test_blank_fold_without_blanks_is_not_evaluable():
    matrix, metadata = _blank_matrix([100., 120.], [])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("blank_fold",
                               threshold={"operator": ">=", "value": 3.0})], None)
    assert out["metrics"]["blank_fold"]["status"] == "not_evaluable"


# ---------- detection_rate ----------

def test_detection_rate_uses_the_detection_mask():
    detected = [[True], [True], [False], [False]]
    matrix = _matrix([[1.], [1.], [1.], [1.]], detected=detected)
    metadata = _rows([{"role": "sample"}] * 4)
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("detection_rate",
                               threshold={"operator": ">=", "value": 0.6})], None)
    metric = out["metrics"]["detection_rate"]
    assert metric["elements"][0]["value"] == pytest.approx(0.5)
    assert metric["status"] == "fail"


def test_detection_rate_without_a_mask_is_not_evaluable():
    """maskが無いのを「全部検出」と読まない。"""
    matrix = _matrix([[1.], [1.]])
    metadata = _rows([{"role": "sample"}] * 2)
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("detection_rate",
                               threshold={"operator": ">=", "value": 0.6})], None)
    metric = out["metrics"]["detection_rate"]
    assert metric["status"] == "not_evaluable"
    assert metric["reason"] == "detection_unknown"


def test_detection_rate_denominator_keeps_every_included_injection():
    """検出0の注入を分母から外さない（外すと検出率が常に1になる）。"""
    detected = [[True], [False], [False], [False]]
    matrix = _matrix([[1.], [np.nan], [np.nan], [np.nan]], detected=detected)
    metadata = _rows([{"role": "sample"}] * 4)
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("detection_rate",
                               threshold={"operator": ">=", "value": 0.6})], None)
    assert out["metrics"]["detection_rate"]["elements"][0]["value"] \
        == pytest.approx(0.25)


def test_excluded_injections_leave_the_denominator():
    detected = [[True], [True], [False]]
    matrix = _matrix([[1.], [1.], [1.]], detected=detected)
    metadata = _rows([{"role": "sample"}, {"role": "sample"},
                      {"role": "sample", "include": False}])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("detection_rate",
                               threshold={"operator": ">=", "value": 0.6})], None)
    assert out["metrics"]["detection_rate"]["elements"][0]["value"] == 1.0


# ---------- 標準品のRT/質量誤差 ----------

def _standard_population(metric) -> dict:
    return {metric: {
        "feature_ids": ["1"], "assay_ids": ["assay[1]", "assay[2]"],
        "targets": {"gaba": {"feature_id": "1", "expected_rt_min": 1.2,
                             "expected_mz": 104.0706}},
        "hash": "h" * 64}}


def _standard_rows(rt_values, mz_values=None) -> list[dict]:
    mz_values = mz_values or [104.0706] * len(rt_values)
    return [{"feature_id": "1", "assay_id": f"assay[{i + 1}]",
             "observed_rt_min": rt, "observed_mz": mz,
             "detection_status": "detected"}
            for i, (rt, mz) in enumerate(zip(rt_values, mz_values))]


def test_standard_rt_error_is_per_injection():
    matrix = _matrix([[1.], [1.]])
    metadata = _rows([{"role": "standard"}, {"role": "standard"}])
    policy = [_policy("standard_rt_error", scope="targets", target_ids=["gaba"],
                      evidence_requirement=True,
                      threshold={"operator": "<=", "value": 0.1})]
    out = evaluate_qc(matrix, _evidence(_standard_rows([1.22, 1.5])), metadata,
                      policy, _standard_population("standard_rt_error"))
    metric = out["metrics"]["standard_rt_error"]
    assert [round(e["value"], 3) for e in metric["elements"]] == [0.02, 0.3]
    assert metric["counts"] == {"pass": 1, "fail": 1, "not_evaluable": 0, "total": 2}
    assert metric["status"] == "fail"


def test_standard_mass_error_is_signed_but_judged_on_absolute_value():
    matrix = _matrix([[1.], [1.]])
    metadata = _rows([{"role": "standard"}, {"role": "standard"}])
    policy = [_policy("standard_mass_error_ppm", scope="targets",
                      target_ids=["gaba"], evidence_requirement=True,
                      threshold={"operator": "<=", "value": 10.0})]
    rows = _standard_rows([1.2, 1.2], [104.0701, 104.0706])
    out = evaluate_qc(matrix, _evidence(rows), metadata, policy,
                      _standard_population("standard_mass_error_ppm"))
    metric = out["metrics"]["standard_mass_error_ppm"]
    assert metric["elements"][0]["value"] < 0          # 符号は残す
    assert metric["status"] == "pass"                  # 判定は絶対値


def test_standard_metric_without_evidence_is_not_evaluable():
    matrix = _matrix([[1.], [1.]])
    metadata = _rows([{"role": "standard"}, {"role": "standard"}])
    policy = [_policy("standard_rt_error", scope="targets", target_ids=["gaba"],
                      evidence_requirement=True,
                      threshold={"operator": "<=", "value": 0.1})]
    out = evaluate_qc(matrix, _evidence(availability=False), metadata, policy,
                      _standard_population("standard_rt_error"))
    metric = out["metrics"]["standard_rt_error"]
    assert metric["status"] == "not_evaluable"
    assert metric["reason"] == "assay_evidence_unavailable"


def test_standard_metric_without_expected_values_is_not_evaluable():
    matrix = _matrix([[1.], [1.]])
    metadata = _rows([{"role": "standard"}, {"role": "standard"}])
    policy = [_policy("standard_rt_error", scope="targets", target_ids=["gaba"],
                      evidence_requirement=True,
                      threshold={"operator": "<=", "value": 0.1})]
    out = evaluate_qc(matrix, _evidence(_standard_rows([1.2])), metadata,
                      policy, None)
    assert out["metrics"]["standard_rt_error"]["status"] == "not_evaluable"


# ---------- internal_standard_valid_fraction ----------

def test_internal_standard_valid_fraction_counts_locked_denominators():
    locked = [[False], [True], [False], [False]]
    matrix = _matrix([[1.], [np.nan], [1.], [1.]], locked=locked,
                     units=["internal_standard_ratio"])
    metadata = _rows([{"role": "sample"}] * 4)
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("internal_standard_valid_fraction",
                               threshold={"operator": ">=", "value": 0.8})], None)
    element = out["metrics"]["internal_standard_valid_fraction"]["elements"][0]
    assert element["value"] == pytest.approx(0.75)
    assert element["status"] == "fail"


# ---------- 母集団の固定 ----------

def test_population_is_fixed_before_filters_and_reused():
    detected = [[True, True], [True, False], [True, False]]
    matrix = _matrix([[1., 1.], [1., 1.], [1., 1.]], detected=detected)
    metadata = _rows([{"role": "sample"}] * 3)
    policy = [_policy("detection_rate",
                      threshold={"operator": ">=", "value": 0.6})]

    first = evaluate_qc(matrix, _evidence(), metadata, policy, None)
    assert first["population"]["detection_rate"]["feature_ids"] == ["1", "2"]

    # feature 2 を解析対象から外しても、QC の母集団は縮まない。
    narrowed = dict(matrix)
    narrowed["eligibility_mask"] = np.array([True, False])
    second = evaluate_qc(narrowed, _evidence(), metadata, policy,
                         first["population"])
    assert second["population"]["detection_rate"]["feature_ids"] == ["1", "2"]
    assert second["metrics"]["detection_rate"]["counts"] \
        == first["metrics"]["detection_rate"]["counts"]


def test_filtering_out_failures_cannot_turn_the_batch_into_a_pass():
    """除外はeligibilityにだけ反映し、QC結果を再集計しない。"""
    detected = [[True, True], [True, False], [True, False]]
    matrix = _matrix([[1., 1.], [1., 1.], [1., 1.]], detected=detected)
    metadata = _rows([{"role": "sample"}] * 3)
    policy = [_policy("detection_rate", minimum_pass_fraction=1.0,
                      threshold={"operator": ">=", "value": 0.6})]

    before = evaluate_qc(matrix, _evidence(), metadata, policy, None)
    assert before["batch_status"] == "fail"

    narrowed = dict(matrix)
    narrowed["eligibility_mask"] = np.array([True, False])
    after = evaluate_qc(narrowed, _evidence(), metadata, policy,
                        before["population"])
    assert after["batch_status"] == "fail"


def test_eligibility_suggestion_lists_failed_features_without_applying_them():
    detected = [[True, True], [True, False], [True, False]]
    matrix = _matrix([[1., 1.], [1., 1.], [1., 1.]], detected=detected)
    metadata = _rows([{"role": "sample"}] * 3)
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("detection_rate",
                               threshold={"operator": ">=", "value": 0.6})], None)
    assert out["eligibility_suggestion"]["exclude_feature_ids"] == ["2"]
    assert matrix["eligibility_mask"].tolist() == [True, True]


# ---------- バッチ判定 ----------

def test_optional_metric_failure_is_a_warning_not_a_batch_failure():
    matrix, metadata = _qc_matrix([1., 10., 100.])
    out = evaluate_qc(matrix, _evidence(), metadata,
                      [_policy("pooled_qc_rsd", required=False)], None)
    assert out["metrics"]["pooled_qc_rsd"]["status"] == "fail"
    assert out["batch_status"] == "not_evaluable"      # 必須項目がゼロ
    assert any("pooled_qc_rsd" in w for w in out["warnings"])


def test_batch_passes_only_when_every_required_metric_passes():
    matrix, metadata = _qc_matrix([100., 101., 99.])
    out = evaluate_qc(matrix, _evidence(), metadata, [_policy("pooled_qc_rsd")], None)
    assert out["batch_status"] == "pass"


def test_unknown_metric_is_rejected():
    matrix, metadata = _qc_matrix([100., 101., 99.])
    from metabolomix.core.atomic_io import DomainError
    with pytest.raises(DomainError) as caught:
        evaluate_qc(matrix, _evidence(), metadata, [_policy("made_up_metric")], None)
    assert caught.value.code == "QC_POLICY_INVALID"


def test_metadata_must_line_up_with_the_matrix_assays():
    matrix, metadata = _qc_matrix([100., 101., 99.])
    from metabolomix.core.atomic_io import DomainError
    with pytest.raises(DomainError) as caught:
        evaluate_qc(matrix, _evidence(), metadata[:-1], [_policy("pooled_qc_rsd")],
                    None)
    assert caught.value.code == "QC_POLICY_INVALID"


# ---------- v2 は自動policyを通さない ----------

def test_explicit_policy_adds_no_automatic_decisions():
    """profileが明示したstepを、証拠不足を理由に黙って落とさない。

    落とす判断（前提不足）は `matrix_state` が QC_PREREQUISITE_MISSING で
    **止める**役目で、conservative-v1 のように skip して続行しない。
    """
    from metabolomix.analysis.preprocess_policy import (
        EXPLICIT_POLICY_VERSION,
        explicit_policy,
    )

    recipe = {"base": "peak_height", "normalize": "median", "drift_correct": True,
              "filter": {"min_detection_rate": 0.5}, "impute": "half_min"}
    plan = explicit_policy(recipe)
    assert plan["policy_version"] == EXPLICIT_POLICY_VERSION
    assert plan["skipped_steps"] == []
    assert plan["assumptions"] == []
    assert plan["resolved_recipe"]["normalize"] == "median"
    assert plan["resolved_recipe"]["drift_correct"] is True
    assert plan["requested_recipe"] == recipe

