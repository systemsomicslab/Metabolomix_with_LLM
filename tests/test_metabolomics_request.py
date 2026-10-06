"""pipeline-request.v2 の解決・更新契約（spec §6, §6.2）を検証する。

request_v2.resolve/validate_statistics/merge_updates はすべて
metabolomix.pipeline.request_v2 にある。既存v1 (metabolomix.pipeline.request) は
schemaディスパッチだけを追加で持ち、v1自身の挙動（既定値・null・更新）は
一切変えない——tests/test_pipeline_request.py の既存assertは変更しない。

profileはlcms-profile.v1のうちrequest_v2が実際に参照するキー（matrix_recipes・
feature_targets・analysis_recipe.statistics）だけを持つ最小dictで足りる
（validate_profileを通す義務はない——request_v2はprofileを検証済みとして受け取る側）。
"""
from __future__ import annotations

import copy
import json

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import request as request_v1
from metabolomix.pipeline import request_v2
from metabolomix.pipeline.request_v2 import merge_updates, resolve, validate_statistics


def _profile(matrix_recipes=None, feature_targets=None, statistics=None):
    return {
        "matrix_recipes": matrix_recipes if matrix_recipes is not None else {
            "default": {"base": "peak_height", "normalize": "none",
                        "drift_correct": False, "filter": None, "impute": "none"},
        },
        "feature_targets": feature_targets if feature_targets is not None else {},
        "analysis_recipe": {
            "statistics": statistics if statistics is not None else [],
            "internal_standards": [],
        },
    }


def _pca(statistic_id="pca", matrix_recipe_id="default", **extra):
    base = {
        "statistic_id": statistic_id, "kind": "pca", "matrix_recipe_id": matrix_recipe_id,
        "transform": "none", "feature_scope": {"mode": "all_eligible"},
    }
    base.update(extra)
    return base


def _welch(statistic_id="welch1", reference_group="control", test_group="treated", **extra):
    base = {
        "statistic_id": statistic_id, "kind": "welch", "matrix_recipe_id": "default",
        "transform": "log2", "feature_scope": {"mode": "all_eligible"},
        "reference_group": reference_group, "test_group": test_group,
    }
    base.update(extra)
    return base


# ---------- brief記載のRED（pca rejects groups. 手を加えず、原文を保つ） ----------

def test_pca_rejects_groups():
    s = {"statistic_id": "p", "kind": "pca", "matrix_recipe_id": "default",
         "transform": "none", "feature_scope": {"mode": "all_eligible"},
         "groups": ["a", "b"]}
    with pytest.raises(DomainError):
        validate_statistics([s], {"matrix_recipes": {"default": {}}})


# ---------- unknown ----------

def test_unknown_top_level_key_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": "profile.json", "bogus_key": 1}, _profile())


def test_unknown_key_in_statistic_item_rejected():
    profile = _profile()
    stat = _pca()
    stat["bogus"] = "x"
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics([stat], profile)


def test_unknown_statistic_kind_rejected():
    profile = _profile()
    stat = _pca()
    stat["kind"] = "boxplot"
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics([stat], profile)


# ---------- null ----------

def test_explicit_null_target_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": "profile.json", "target": None}, _profile())


def test_explicit_null_statistics_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": "profile.json", "statistics": None}, _profile())


def test_explicit_null_profile_file_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": None}, _profile())


def test_sample_manifest_explicit_null_is_allowed_like_v1():
    resolved = resolve(
        {"profile_file": "profile.json", "sample_manifest": None,
         "statistics": [_pca()]},
        _profile(),
    )
    assert resolved["sample_manifest"] is None
    assert resolved["value_sources"]["sample_manifest"] == "explicit"


# ---------- 旧comparisons ----------

def test_old_comparisons_field_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json",
             "comparisons": [{"comparison_id": "a_vs_b",
                               "reference_group": "a", "test_group": "b"}]},
            _profile(),
        )


# ---------- method_file直接指定 ----------

def test_method_file_direct_specification_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": "profile.json", "method_file": "method.txt"}, _profile())


def test_lbm_file_direct_specification_rejected():
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve({"profile_file": "profile.json", "lbm_file": "lib.lbm"}, _profile())


# ---------- 重複statistic_id ----------

def test_duplicate_statistic_id_rejected():
    profile = _profile()
    stats = [_pca(statistic_id="dup"), _pca(statistic_id="dup")]
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics(stats, profile)


# ---------- target不整合 ----------

def test_exploratory_target_with_test_statistic_rejected():
    profile = _profile()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "target": "exploratory",
             "statistics": [_welch()]},
            profile,
        )


def test_differential_target_without_test_statistic_rejected():
    profile = _profile()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "target": "differential",
             "statistics": [_pca()]},
            profile,
        )


def test_target_auto_resolves_to_differential_with_welch():
    profile = _profile()
    resolved = resolve(
        {"profile_file": "profile.json", "statistics": [_welch()]}, profile,
    )
    assert resolved["target"] == "auto"
    assert resolved["effective_target"] == "differential"


def test_target_auto_resolves_to_exploratory_with_only_pca():
    profile = _profile()
    resolved = resolve(
        {"profile_file": "profile.json", "statistics": [_pca()]}, profile,
    )
    assert resolved["effective_target"] == "exploratory"


# ---------- schema省略はv1 ----------

def test_schema_omission_defaults_to_v1(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    resolved = request_v1.resolve_request(root, {})
    assert resolved["schema"] == "pipeline-request.v1"


def test_schema_omission_in_explicit_dict_still_dispatches_v1_even_with_profile(tmp_path):
    """profileを渡しても、schemaを明示しない限りv2を推測しない(A01)。"""
    root = tmp_path / "source"
    root.mkdir()
    resolved = request_v1.resolve_request(root, {}, profile=_profile())
    assert resolved["schema"] == "pipeline-request.v1"


# ---------- routine範囲外override ----------

def _profile_with_two_recipes():
    return _profile(matrix_recipes={
        "default": {"base": "peak_height", "normalize": "none",
                    "drift_correct": False, "filter": None, "impute": "none"},
        "alt": {"base": "peak_height", "normalize": "tic",
                "drift_correct": True, "filter": None, "impute": "half_min"},
    })


def test_routine_preprocess_override_omitted_routine_overrides_fails_closed():
    """routine_overridesを渡さない(既定None)なら、routineでは何も上書き
    できない——fail-closed。構造的にどれだけ妥当な値でも通らない。"""
    profile = _profile_with_two_recipes()
    with pytest.raises(DomainError, match="PROFILE_SCOPE_MISMATCH"):
        resolve(
            {"profile_file": "profile.json",
             "preprocess": {"default": {"normalize": "tic"}},
             "statistics": [_pca()]},
            profile,
        )


def test_routine_preprocess_override_outside_allowance_rejected():
    profile = _profile_with_two_recipes()
    routine_overrides = {
        "preprocess": {"default": {"normalize": {"mode": "values", "values": ["tic"]}}},
    }
    with pytest.raises(DomainError, match="PROFILE_SCOPE_MISMATCH"):
        resolve(
            {"profile_file": "profile.json",
             "preprocess": {"default": {"normalize": "median"}},
             "statistics": [_pca()]},
            profile,
            routine_overrides=routine_overrides,
        )


def test_routine_preprocess_override_allowed_by_any_mode():
    profile = _profile_with_two_recipes()
    routine_overrides = {"preprocess": {"default": {"normalize": {"mode": "any"}}}}
    resolved = resolve(
        {"profile_file": "profile.json",
         "preprocess": {"default": {"normalize": "median"}},
         "statistics": [_pca()]},
        profile,
        routine_overrides=routine_overrides,
    )
    assert resolved["preprocess"]["default"]["normalize"] == "median"


def test_routine_preprocess_override_allowed_by_values_mode():
    profile = _profile_with_two_recipes()
    routine_overrides = {
        "preprocess": {"default": {"normalize": {"mode": "values", "values": ["tic"]}}},
    }
    resolved = resolve(
        {"profile_file": "profile.json",
         "preprocess": {"default": {"normalize": "tic"}},
         "statistics": [_pca()]},
        profile,
        routine_overrides=routine_overrides,
    )
    assert resolved["preprocess"]["default"]["normalize"] == "tic"


def test_routine_preprocess_override_field_not_named_in_allowance_rejected():
    """recipe_idは許容されていても、そのfield自体が宣言に無ければ拒否
    （fieldごとの粒度）。"""
    profile = _profile_with_two_recipes()
    routine_overrides = {
        "preprocess": {"default": {"impute": {"mode": "any"}}},  # normalizeは宣言されていない
    }
    with pytest.raises(DomainError, match="PROFILE_SCOPE_MISMATCH"):
        resolve(
            {"profile_file": "profile.json",
             "preprocess": {"default": {"normalize": "tic"}},
             "statistics": [_pca()]},
            profile,
            routine_overrides=routine_overrides,
        )


def test_validation_purpose_ignores_routine_overrides_entirely():
    """execution_purpose='validation'ではroutine_overridesを一切参照しない
    ——まだ証明書が検証していない値を試すのがvalidationの目的だから。
    routine_overrides=None（何も許容していない）でも構造的に妥当な上書きは
    通る。"""
    profile = _profile_with_two_recipes()
    resolved = resolve(
        {"profile_file": "profile.json", "execution_purpose": "validation",
         "preprocess": {"default": {"normalize": "median"}},
         "statistics": [_pca()]},
        profile,
        routine_overrides=None,
    )
    assert resolved["preprocess"]["default"]["normalize"] == "median"


def test_resolve_request_dispatch_passes_routine_overrides_through(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    profile = _profile_with_two_recipes()
    routine_overrides = {"preprocess": {"default": {"normalize": {"mode": "any"}}}}
    resolved = request_v1.resolve_request(
        root,
        {"schema": "pipeline-request.v2", "profile_file": "profile.json",
         "preprocess": {"default": {"normalize": "median"}}, "statistics": [_pca()]},
        profile=profile, routine_overrides=routine_overrides,
    )
    assert resolved["preprocess"]["default"]["normalize"] == "median"


def test_merge_updates_dispatch_passes_routine_overrides_through():
    profile = _profile_with_two_recipes()
    routine_overrides = {"preprocess": {"default": {"normalize": {"mode": "any"}}}}
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    updated = request_v1.merge_updates(
        resolved, {"preprocess": {"default": {"normalize": "median"}}},
        profile=profile, routine_overrides=routine_overrides,
    )
    assert updated["preprocess"]["default"]["normalize"] == "median"


def test_preprocess_override_unknown_recipe_rejected():
    profile = _profile()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json",
             "preprocess": {"nonexistent": {"normalize": "none"}},
             "statistics": [_pca()]},
            profile,
        )


def test_preprocess_override_double_normalization_rejected():
    profile = _profile(matrix_recipes={
        "default": {"base": "internal_standard_ratio", "normalize": "none",
                    "drift_correct": False, "filter": None, "impute": "none"},
    })
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json",
             "preprocess": {"default": {"normalize": "tic"}},
             "statistics": [_pca()]},
            profile,
        )


# ---------- 既定値 ----------

def test_statistics_omitted_falls_back_to_profile_default():
    profile = _profile(statistics=[_welch(statistic_id="welch1")])
    resolved = resolve({"profile_file": "profile.json"}, profile)
    assert resolved["value_sources"]["statistics"] == "profile_default"
    assert [s["statistic_id"] for s in resolved["statistics"]] == ["welch1"]
    assert resolved["effective_target"] == "differential"


def test_statistics_omitted_and_profile_has_none_uses_global_pca_default():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json"}, profile)
    assert resolved["value_sources"]["statistics"] == "v2_default"
    assert resolved["statistics"] == [{
        "statistic_id": "pca", "kind": "pca", "matrix_recipe_id": "default",
        "transform": "none", "feature_scope": {"mode": "all_eligible"},
        "scaling": "autoscale", "n_components": 2,
    }]
    assert resolved["effective_target"] == "exploratory"


def test_pca_defaults_scaling_and_n_components():
    profile = _profile()
    stat = {"statistic_id": "p", "kind": "pca", "matrix_recipe_id": "default",
            "transform": "none", "feature_scope": {"mode": "all_eligible"}}
    normalized = validate_statistics([stat], profile)
    assert normalized[0]["scaling"] == "autoscale"
    assert normalized[0]["n_components"] == 2


def test_welch_defaults_thresholds():
    profile = _profile()
    stat = {"statistic_id": "w", "kind": "welch", "matrix_recipe_id": "default",
            "transform": "log2", "feature_scope": {"mode": "all_eligible"},
            "reference_group": "a", "test_group": "b"}
    normalized = validate_statistics([stat], profile)
    assert normalized[0]["q_threshold"] == 0.05
    assert normalized[0]["log2fc_threshold"] == 1.0


def test_anova_requires_three_or_more_groups():
    profile = _profile()
    stat = {"statistic_id": "an", "kind": "anova_tukey", "matrix_recipe_id": "default",
            "transform": "none", "feature_scope": {"mode": "all_eligible"},
            "groups": ["a", "b"]}
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics([stat], profile)


def test_feature_scope_targets_must_reference_existing_target_id():
    profile = _profile(feature_targets={"gaba": {}})
    stat = _pca()
    stat["feature_scope"] = {"mode": "targets", "target_ids": ["not_gaba"]}
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics([stat], profile)


def test_matrix_recipe_id_must_exist_in_profile():
    profile = _profile()
    stat = _pca(matrix_recipe_id="nonexistent")
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        validate_statistics([stat], profile)


# ---------- standard_assays ----------

def test_standard_assays_unknown_target_id_rejected():
    profile = _profile(feature_targets={"gaba": {}})
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "standard_assays": {"not_gaba": ["sample-1"]}},
            profile,
        )


def test_standard_assays_valid_structure_accepted():
    profile = _profile(feature_targets={"gaba": {}})
    resolved = resolve(
        {"profile_file": "profile.json", "statistics": [_pca()],
         "standard_assays": {"gaba": ["sample-1", "sample-2"]}},
        profile,
    )
    assert resolved["standard_assays"] == {"gaba": ["sample-1", "sample-2"]}


def test_standard_assays_omitted_defaults_to_empty():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    assert resolved["standard_assays"] == {}


# ---------- 更新不可キー ----------

def test_update_profile_file_requires_new_pipeline():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(resolved, {"profile_file": "other.json"}, profile)


def test_update_execution_purpose_requires_new_pipeline():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    with pytest.raises(DomainError, match="NEW_PIPELINE_REQUIRED"):
        merge_updates(resolved, {"execution_purpose": "validation"}, profile)


def test_update_statistics_is_allowed_and_revalidated():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    updated = merge_updates(resolved, {"statistics": [_welch()]}, profile)
    assert [s["statistic_id"] for s in updated["statistics"]] == ["welch1"]
    assert updated["effective_target"] == "differential"
    assert updated["value_sources"]["statistics"] == "explicit_update"


def test_update_preserves_value_sources_for_untouched_fields():
    profile = _profile()
    resolved = resolve(
        {"profile_file": "profile.json", "statistics": [_pca()],
         "sample_manifest": "manifest.tsv"},
        profile,
    )
    assert resolved["value_sources"]["sample_manifest"] == "explicit"
    updated = merge_updates(resolved, {"target": "differential", "statistics": [_welch()]}, profile)
    # targetを更新してもsample_manifestの出所は変わらない
    assert updated["value_sources"]["sample_manifest"] == "explicit"
    assert updated["value_sources"]["target"] == "explicit_update"


def test_merge_updates_dispatch_from_v1_module(tmp_path):
    """既存request.merge_updatesがschema="pipeline-request.v2"のrequestを
    request_v2.merge_updatesへ委譲することを確認する。"""
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    updated = request_v1.merge_updates(resolved, {"target": "exploratory"}, profile=profile)
    assert updated["schema"] == request_v2.SCHEMA
    assert updated["target"] == "exploratory"


def test_resolve_request_dispatch_from_v1_module(tmp_path):
    """既存request.resolve_requestがexplicit["schema"]=="pipeline-request.v2"を
    検出してrequest_v2.resolveへ委譲することを確認する。"""
    root = tmp_path / "source"
    root.mkdir()
    profile = _profile()
    resolved = request_v1.resolve_request(
        root, {"schema": "pipeline-request.v2", "profile_file": "profile.json",
               "statistics": [_pca()]},
        profile=profile,
    )
    assert resolved["schema"] == "pipeline-request.v2"
    assert resolved["profile_file"] == str(root / "profile.json")


# ---------- resolveは入力を書き換えない ----------

def test_resolve_does_not_mutate_input_profile_or_data():
    profile = _profile(statistics=[_welch()])
    profile_copy = copy.deepcopy(profile)
    data = {"profile_file": "profile.json"}
    data_copy = copy.deepcopy(data)
    resolve(data, profile)
    assert profile == profile_copy
    assert data == data_copy


# ---------- Concern 1: analysis-request.json の4段優先順位 ----------
# spec §6「値の優先順位はMCP明示値 > analysis-request.json > profile既定値 >
# v2既定値」。metabolomix.pipeline.request.resolve_request が v1 と同じ場所
# （request層自身）でファイルを読むので、v2でもここで読む。

def _write_request_file(root, payload):
    (root / "analysis-request.json").write_text(json.dumps(payload), encoding="utf-8")


def test_request_file_fills_gaps_left_by_explicit(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _write_request_file(root, {
        "schema": "pipeline-request.v2",
        "profile_file": "profile-from-file.json",
        "sample_manifest": "manifest-from-file.tsv",
    })
    profile = _profile()
    resolved = request_v1.resolve_request(root, {"statistics": [_pca()]}, profile=profile)
    assert resolved["schema"] == request_v2.SCHEMA
    assert resolved["profile_file"] == str(root / "profile-from-file.json")
    assert resolved["value_sources"]["profile_file"] == "request_file"
    assert resolved["sample_manifest"] == "manifest-from-file.tsv"
    assert resolved["value_sources"]["sample_manifest"] == "request_file"
    # explicitで渡したstatisticsはfileより優先(そもそもfile側に無い)
    assert resolved["value_sources"]["statistics"] == "explicit"


def test_explicit_value_beats_request_file_for_v2(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _write_request_file(root, {
        "schema": "pipeline-request.v2", "profile_file": "profile-from-file.json",
    })
    profile = _profile()
    resolved = request_v1.resolve_request(
        root,
        {"schema": "pipeline-request.v2", "profile_file": "profile-from-explicit.json",
         "statistics": [_pca()]},
        profile=profile,
    )
    assert resolved["profile_file"] == str(root / "profile-from-explicit.json")
    assert resolved["value_sources"]["profile_file"] == "explicit"


def test_request_file_schema_v2_dispatches_without_explicit_schema_key(tmp_path):
    """schemaをexplicit側で省略しても、analysis-request.json自身がv2を宣言して
    いればv2として解決する——explicit/file双方が省略したときだけがA01の
    「schema省略はv1」に該当する。"""
    root = tmp_path / "source"
    root.mkdir()
    _write_request_file(root, {
        "schema": "pipeline-request.v2", "profile_file": "profile-from-file.json",
    })
    profile = _profile()
    resolved = request_v1.resolve_request(root, {"statistics": [_pca()]}, profile=profile)
    assert resolved["schema"] == request_v2.SCHEMA


def test_request_file_v2_shape_rejects_v1_only_keys(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    _write_request_file(root, {
        "schema": "pipeline-request.v2", "profile_file": "profile.json",
        "method_file": "method.txt",
    })
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        request_v1.resolve_request(root, {"statistics": [_pca()]}, profile=_profile())


def test_no_request_file_v2_behaves_as_before(tmp_path):
    """ファイルが無い場合の既存の3段（explicit > profile既定 > v2既定）は
    そのまま——回帰していないことを確認する。"""
    root = tmp_path / "source"
    root.mkdir()
    profile = _profile()
    resolved = request_v1.resolve_request(
        root,
        {"schema": "pipeline-request.v2", "profile_file": "profile.json",
         "statistics": [_pca()]},
        profile=profile,
    )
    assert resolved["value_sources"]["sample_manifest"] == "default"


def test_broken_request_file_still_rejected_when_schema_omitted_everywhere(tmp_path):
    """壊れたJSONを黙ってv1既定へフォールバックしない
    （v1の既存挙動: test_a_broken_request_file_is_rejected_not_ignored と同じ
    不変条件をv2追加後も保つことを確認する）。"""
    root = tmp_path / "source"
    root.mkdir()
    (root / "analysis-request.json").write_text("{ not json", encoding="utf-8")
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        request_v1.resolve_request(root)


# ---------- Concern 3: feature_bindings（spec §6.2） ----------
# 「更新payloadはdataset hash、target_idごとのfeature_id、選択理由」の構造検証
# だけがこのモジュールの責務。候補解決・許容規則との突合（実データ照合）は
# 後続task（resolve_feature_bindings stage）が担う——ここでは実装しない。

_DATASET_HASH = "a" * 64


def _feature_bindings(target_id="gaba", feature_id="F00123", reason="single candidate"):
    return {
        "dataset_hash": _DATASET_HASH,
        "selections": {target_id: {"feature_id": feature_id, "reason": reason}},
    }


def test_feature_bindings_valid_structure_accepted_on_initial_resolve():
    profile = _profile(feature_targets={"gaba": {}})
    resolved = resolve(
        {"profile_file": "profile.json", "statistics": [_pca()],
         "feature_bindings": _feature_bindings()},
        profile,
    )
    assert resolved["feature_bindings"]["dataset_hash"] == _DATASET_HASH
    assert resolved["feature_bindings"]["selections"]["gaba"]["feature_id"] == "F00123"


def test_feature_bindings_omitted_defaults_to_none():
    profile = _profile()
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    assert resolved["feature_bindings"] is None


def test_feature_bindings_explicit_null_rejected():
    profile = _profile()
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "feature_bindings": None},
            profile,
        )


def test_feature_bindings_unknown_target_id_rejected_on_initial_resolve():
    profile = _profile(feature_targets={"gaba": {}})
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "feature_bindings": _feature_bindings(target_id="not_gaba")},
            profile,
        )


def test_feature_bindings_bad_dataset_hash_rejected():
    profile = _profile(feature_targets={"gaba": {}})
    bindings = _feature_bindings()
    bindings["dataset_hash"] = "not-a-sha256"
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "feature_bindings": bindings},
            profile,
        )


def test_feature_bindings_missing_reason_rejected():
    profile = _profile(feature_targets={"gaba": {}})
    bindings = {"dataset_hash": _DATASET_HASH, "selections": {"gaba": {"feature_id": "F1"}}}
    with pytest.raises(DomainError, match="PIPELINE_REQUEST_INVALID"):
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "feature_bindings": bindings},
            profile,
        )


def test_feature_bindings_is_updatable_via_resume():
    profile = _profile(feature_targets={"gaba": {}})
    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    assert resolved["feature_bindings"] is None
    updated = merge_updates(resolved, {"feature_bindings": _feature_bindings()}, profile)
    assert updated["feature_bindings"]["selections"]["gaba"]["feature_id"] == "F00123"
    assert updated["value_sources"]["feature_bindings"] == "explicit_update"


def test_feature_bindings_initial_and_resume_go_through_identical_validation():
    """初回指定とresumeが同じ検査関数を通ることをピン留めする
    （brief「feature_bindingsは初回指定もresumeも同じ検査を通す」）。同じ
    不正payload（存在しないtarget_id参照）が、初回resolveでもresumeの
    merge_updatesでも同一のDomainError codeで拒否されることを確認する。
    """
    profile = _profile(feature_targets={"gaba": {}})
    bad_bindings = _feature_bindings(target_id="not_gaba")

    with pytest.raises(DomainError) as initial_exc:
        resolve(
            {"profile_file": "profile.json", "statistics": [_pca()],
             "feature_bindings": bad_bindings},
            profile,
        )

    resolved = resolve({"profile_file": "profile.json", "statistics": [_pca()]}, profile)
    with pytest.raises(DomainError) as resume_exc:
        merge_updates(resolved, {"feature_bindings": bad_bindings}, profile)

    assert initial_exc.value.code == resume_exc.value.code == "PIPELINE_REQUEST_INVALID"
