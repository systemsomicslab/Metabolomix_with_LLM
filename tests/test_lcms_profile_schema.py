"""`lcms-profile.v1` と `lcms-profile-validation.v1` の厳密な契約検証（spec §5, §5.1）。

`metabolomix.console.profile_schema` はファイルを読まない純ロジックなので、fixtureは
すべてこのテスト自身が組み立てるPython dict（実データ・実profileに依存しない）。
"""
import copy
import math

import pytest

from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.console.profile_schema import (
    CERTIFICATE_SCHEMA,
    SCHEMA,
    profile_content_hash,
    validate_certificate,
    validate_profile,
)


# ---------- fixture: 合成の最小合法profile ----------

def _base_profile() -> dict:
    return {
        "schema": SCHEMA,
        "profile_id": "kanzo-lcms-metabolomics",
        "revision": 1,
        "omics": "metabolomics",
        "acquisition": {
            "separation": "lc",
            "acquisition_type": "dda",
            "polarity": "positive",
            "instrument": "ExionLC AD / ZenoTOF 7600",
            "lc": {
                "column": "ACQUITY UPLC BEH 1.7um 100x2.1mm",
                "mobile_phase_a": "0.1% formic acid in water",
                "mobile_phase_b": "methanol",
                "flow_rate_ul_min": 150.0,
                "column_temperature_c": 30.0,
                "gradient_profile": None,
            },
            "ms_range": {
                "ms1_low_mz": 50.0,
                "ms1_high_mz": 1500.0,
                "ms2_low_mz": None,
                "ms2_high_mz": None,
            },
            "sample_matrix": "licorice extract",
            "scope": "single LC method, positive polarity, DDA acquisition, peak height only",
        },
        "software": {
            "msdial_version": "5.x",
            "executable_path": "C:/tools/MSDIAL5/MSDIALCUI.exe",
            "executable_sha256": "a" * 64,
            "adapter_version": "1.0.0",
        },
        "processing": {
            "method_path": "method/kanzo_pos_param.txt",
            "method_sha256": "b" * 64,
            "measure": "peak_height",
            "dependencies": [
                {
                    "dependency_id": "msp-lib-1",
                    "kind": "msp",
                    "method_key": "MSP file",
                    "path": "library/lib1.msp",
                    "sha256": "c" * 64,
                    "required": True,
                },
            ],
            "effective_settings": {"Minimum peak height": "1000"},
        },
        "analysis_recipe": {
            "statistics": [
                {
                    "statistic_id": "pca",
                    "kind": "pca",
                    "matrix_recipe_id": "default",
                    "transform": "none",
                    "feature_scope": {"mode": "all_eligible"},
                    "scaling": "autoscale",
                    "n_components": 2,
                },
            ],
            "internal_standards": [
                {"target_id": "gaba", "standard_target_id": "gaba_d6"},
            ],
        },
        "feature_targets": {
            "gaba": {
                "compound_identifiers": {"name": "GABA", "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
                "adduct": "[M+H]+",
                "charge": 1,
                "expected_mz": 104.0706,
                "mz_tolerance_ppm": 10.0,
                "expected_rt_min": 1.2,
                "rt_tolerance_min": 0.1,
                "required_evidence": [{"kind": "mass_rt"}],
            },
            "gaba_d6": {
                "compound_identifiers": {"name": "GABA-d6", "inchikey": "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"},
                "adduct": "[M+H]+",
                "charge": 1,
                "expected_mz": 110.11,
                "mz_tolerance_ppm": 10.0,
                "expected_rt_min": 1.2,
                "rt_tolerance_min": 0.1,
                "required_evidence": [{"kind": "mass_rt"}],
            },
        },
        "matrix_recipes": {
            "default": {
                "base": "peak_height",
                "normalize": "none",
                "drift_correct": False,
                "filter": None,
                "impute": "none",
            },
        },
        "qc_policy": {
            "pooled_qc_rsd": {
                "scope": "all_features",
                "target_ids": None,
                "required": True,
                "threshold": {"operator": "<=", "value": 30.0},
                "evidence_requirement": False,
                "minimum_pass_fraction": 0.8,
            },
        },
        "evidence": {
            "acquisition.instrument": {
                "value": "ExionLC AD / ZenoTOF 7600",
                "reason": None,
                "tier": "paper_explicit",
                "source_uri": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12738721/",
                "source_hash": None,
                "location": "Methods section",
            },
            "acquisition.lc.gradient_profile": {
                "value": None,
                "reason": "論文本文にLC時間プログラムの記載がない。",
                "tier": "proposed",
                "source_uri": None,
                "source_hash": None,
                "location": None,
            },
            "acquisition.ms_range.ms2_low_mz": {
                "value": None,
                "reason": "甘草抽出物methodの本文にMS2質量範囲の明示的な記載がない。",
                "tier": "proposed",
                "source_uri": None,
                "source_hash": None,
                "location": None,
            },
            "acquisition.ms_range.ms2_high_mz": {
                "value": None,
                "reason": "甘草抽出物methodの本文にMS2質量範囲の明示的な記載がない。",
                "tier": "proposed",
                "source_uri": None,
                "source_hash": None,
                "location": None,
            },
        },
        "validation": {
            "status": "draft",
            "scope": "single LC method, positive polarity, DDA acquisition, peak height only",
            "certificate_path": None,
            "certificate_sha256": None,
        },
    }


def _validated_profile(certificate_overrides: dict | None = None) -> dict:
    """draftを土台に、hash-pinされた証明書を持つvalidated profileと対応証明書を作る。

    ``certificate_overrides``は証明書へ追加/上書きする任意キー（例:
    ``routine_overrides``）——hashは上書き後の内容から計算するので、常に
    整合したcertificate_sha256になる。
    """
    draft = _base_profile()
    content_hash = profile_content_hash(draft)

    certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "profile_content_sha256": content_hash,
        "dependency_hashes": {"msp-lib-1": "c" * 64},
        "fixed_input_hashes": {"validation-run-1": "d" * 64},
        "reference_file_hashes": {"reference-1": "e" * 64},
        "validation_output_hashes": {"output-1": "f" * 64},
        "criteria": [
            {"criterion_id": "mass_rt_pass", "status": "pass", "required": True, "detail": None},
            {"criterion_id": "optional_check", "status": "not_evaluable", "required": False, "detail": None},
        ],
        "performed_by": "researcher@example.org",
        "performed_at": "2026-09-15T00:00:00Z",
        "scope": draft["acquisition"]["scope"],
    }
    if certificate_overrides:
        certificate.update(certificate_overrides)
    certificate_sha256 = canonical_hash(certificate)

    profile = copy.deepcopy(draft)
    profile["validation"] = {
        "status": "validated",
        "scope": draft["acquisition"]["scope"],
        "certificate_path": "certificates/validation-1.json",
        "certificate_sha256": certificate_sha256,
    }

    observed_hashes = {
        "dependencies": {"msp-lib-1": "c" * 64},
        "fixed_inputs": {"validation-run-1": "d" * 64},
        "reference_files": {"reference-1": "e" * 64},
        "validation_outputs": {"output-1": "f" * 64},
    }
    return profile, certificate, observed_hashes


# ---------- brief記載のRED ----------

def test_certificate_metadata_does_not_change_content_identity():
    p = {"schema": "lcms-profile.v1", "validation": {"status": "draft"}}
    q = {**p, "validation": {"status": "validated", "certificate_sha256": "a" * 64}}
    assert profile_content_hash(p) == profile_content_hash(q)


# ---------- 正常系: 最小合法profile ----------

def test_valid_draft_profile_round_trips():
    normalized = validate_profile(_base_profile())
    assert normalized["schema"] == SCHEMA
    assert normalized["validation"]["status"] == "draft"
    # 未知キーが混入していない・正規化後も全トップレベルキーが揃っている
    assert set(normalized) == set(_base_profile())


def test_valid_validated_profile_and_certificate_pass():
    profile, certificate, observed_hashes = _validated_profile()
    normalized = validate_profile(profile)
    assert normalized["validation"]["status"] == "validated"
    # 例外を投げなければ合格
    validate_certificate(normalized, certificate, observed_hashes)


# ---------- 未知キー ----------

def test_unknown_top_level_key_rejected():
    data = _base_profile()
    data["unexpected_field"] = "x"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_unknown_nested_key_rejected():
    data = _base_profile()
    data["acquisition"]["unexpected"] = "x"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_unknown_feature_target_key_rejected():
    data = _base_profile()
    data["feature_targets"]["gaba"]["unexpected"] = "x"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- NaN ----------

def test_nan_in_numeric_field_rejected():
    data = _base_profile()
    data["feature_targets"]["gaba"]["expected_mz"] = math.nan
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_infinity_rejected():
    data = _base_profile()
    data["feature_targets"]["gaba"]["mz_tolerance_ppm"] = math.inf
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 非整数revision ----------

def test_non_integer_revision_rejected():
    data = _base_profile()
    data["revision"] = 1.5
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_string_revision_rejected():
    data = _base_profile()
    data["revision"] = "1"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_zero_revision_rejected():
    data = _base_profile()
    data["revision"] = 0
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_bool_revision_rejected():
    data = _base_profile()
    data["revision"] = True
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 負ppm ----------

def test_negative_mz_tolerance_ppm_rejected():
    data = _base_profile()
    data["feature_targets"]["gaba"]["mz_tolerance_ppm"] = -10.0
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_zero_mz_tolerance_ppm_rejected():
    data = _base_profile()
    data["feature_targets"]["gaba"]["mz_tolerance_ppm"] = 0.0
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 重複target/recipe ----------

def test_duplicate_internal_standard_target_mapping_rejected():
    data = _base_profile()
    # 3つ目のtargetを追加し、gabaがgaba_d6とgaba_alt2の両方に対応付けられる状態を作る
    data["feature_targets"]["gaba_alt2"] = copy.deepcopy(data["feature_targets"]["gaba_d6"])
    data["analysis_recipe"]["internal_standards"].append(
        {"target_id": "gaba", "standard_target_id": "gaba_alt2"})
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_duplicate_statistic_id_rejected():
    data = _base_profile()
    data["analysis_recipe"]["statistics"].append(
        copy.deepcopy(data["analysis_recipe"]["statistics"][0]))
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 循環内部標準 ----------

def test_circular_internal_standard_mapping_rejected():
    data = _base_profile()
    # gaba -> gaba_d6 (既存) に加え gaba_d6 -> gaba を足して循環させる
    data["analysis_recipe"]["internal_standards"].append(
        {"target_id": "gaba_d6", "standard_target_id": "gaba"})
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_self_referential_internal_standard_rejected():
    data = _base_profile()
    data["analysis_recipe"]["internal_standards"] = [
        {"target_id": "gaba", "standard_target_id": "gaba"},
    ]
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 不存在recipe参照 ----------

def test_statistic_referencing_missing_matrix_recipe_rejected():
    data = _base_profile()
    data["analysis_recipe"]["statistics"][0]["matrix_recipe_id"] = "does-not-exist"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_internal_standard_referencing_missing_target_rejected():
    data = _base_profile()
    data["analysis_recipe"]["internal_standards"] = [
        {"target_id": "gaba", "standard_target_id": "no-such-target"},
    ]
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 既定recipe欠落 ----------

def test_missing_default_matrix_recipe_rejected():
    data = _base_profile()
    data["matrix_recipes"] = {
        "custom": data["matrix_recipes"]["default"],
    }
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


# ---------- 証明書: hash一致・必須判定・適用範囲・改変検知 ----------

def test_certificate_missing_is_rejected():
    profile, _certificate, observed_hashes = _validated_profile()
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, None, observed_hashes)


def test_hand_written_validated_status_without_matching_certificate_is_rejected():
    """profileを手書きvalidatedにしても証明書不足なら拒否する。"""
    data = _base_profile()
    data["validation"] = {
        "status": "validated",
        "scope": data["acquisition"]["scope"],
        "certificate_path": "certificates/fake.json",
        # 実在しない証明書と対応しないでたらめなhash
        "certificate_sha256": "0" * 64,
    }
    normalized = validate_profile(data)
    fake_certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "profile_content_sha256": profile_content_hash(normalized),
        "dependency_hashes": {},
        "fixed_input_hashes": {},
        "reference_file_hashes": {},
        "validation_output_hashes": {},
        "criteria": [{"criterion_id": "x", "status": "pass", "required": True, "detail": None}],
        "performed_by": "someone",
        "performed_at": "2026-09-15T00:00:00Z",
        "scope": data["acquisition"]["scope"],
    }
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, fake_certificate, {})


@pytest.mark.parametrize("mutate", [
    lambda cert: cert.__setitem__("performed_by", "someone-else"),
    lambda cert: cert["dependency_hashes"].__setitem__("msp-lib-1", "9" * 64),
    lambda cert: cert["criteria"].append(
        {"criterion_id": "extra", "status": "fail", "required": False, "detail": None}),
    lambda cert: cert.__setitem__("scope", "a different scope"),
])
def test_certificate_single_element_tampering_is_rejected(mutate):
    """証明書の一要素改変も拒否する（GREEN確認の核心）。"""
    profile, certificate, observed_hashes = _validated_profile()
    normalized = validate_profile(profile)
    tampered = copy.deepcopy(certificate)
    mutate(tampered)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, tampered, observed_hashes)


def test_certificate_hash_mismatch_with_observed_hashes_rejected():
    profile, certificate, observed_hashes = _validated_profile()
    normalized = validate_profile(profile)
    tampered_observed = copy.deepcopy(observed_hashes)
    tampered_observed["dependencies"]["msp-lib-1"] = "1" * 64
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, tampered_observed)


def test_certificate_required_criterion_not_passing_rejected():
    profile, certificate, observed_hashes = _validated_profile()
    draft_for_hash = copy.deepcopy(profile)

    failing_certificate = copy.deepcopy(certificate)
    failing_certificate["criteria"][0]["status"] = "fail"
    new_sha = canonical_hash(failing_certificate)

    profile["validation"]["certificate_sha256"] = new_sha
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, failing_certificate, observed_hashes)


def test_certificate_scope_mismatch_with_profile_scope_rejected():
    profile, certificate, observed_hashes = _validated_profile()
    mismatched_certificate = copy.deepcopy(certificate)
    mismatched_certificate["scope"] = "a totally different scope"
    new_sha = canonical_hash(mismatched_certificate)
    profile["validation"]["certificate_sha256"] = new_sha
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, mismatched_certificate, observed_hashes)


def test_draft_profile_rejects_certificate_check():
    data = _base_profile()
    normalized = validate_profile(data)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, {}, {})


# ---------- fix round 1: evidence must actually back null conditions and real values ----------

def test_null_condition_without_matching_evidence_entry_rejected():
    """測定条件の未記載値には理由付きevidenceが必須（nullを装置既定値と解釈しない）。

    reviewerの反例そのもの: nullなacquisition条件がありながらevidence={}でも
    以前は通ってしまっていた。
    """
    data = _base_profile()
    data["evidence"] = {}
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_null_condition_missing_just_one_evidence_entry_rejected():
    """複数のnull条件のうち1つだけevidenceを外しても検出される。"""
    data = _base_profile()
    del data["evidence"]["acquisition.ms_range.ms2_low_mz"]
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_evidence_value_contradicting_actual_field_value_rejected():
    """evidence[path]["value"]は実際のprofile値と一致しなければならない。

    reviewerの反例: 非nullなacquisition.instrumentに矛盾するevidence.valueを
    与えても以前は通ってしまっていた。
    """
    data = _base_profile()
    data["evidence"]["acquisition.instrument"]["value"] = "a different instrument entirely"
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_evidence_referencing_nonexistent_field_path_rejected():
    data = _base_profile()
    data["evidence"]["acquisition.nonexistent_field"] = {
        "value": None,
        "reason": "でっちあげのfield path。",
        "tier": "proposed",
        "source_uri": None,
        "source_hash": None,
        "location": None,
    }
    with pytest.raises(DomainError, match="PROFILE_INVALID"):
        validate_profile(data)


def test_evidence_entry_for_present_condition_with_matching_value_accepted():
    """value一致・reason省略（nullでないので任意）の正常系は引き続き通る。"""
    data = _base_profile()
    normalized = validate_profile(data)
    assert normalized["evidence"]["acquisition.instrument"]["value"] == \
        normalized["acquisition"]["instrument"]


# ---------- routine_overrides（証明書の明示許容集合。spec §6） ----------
# controller裁定: routineが上書きできる範囲は証明書自身が明示的に宣言する。
# 省略・空はfail-closed（何も上書きを許可しない）。値からの推測は一切しない。

def test_routine_overrides_absent_keeps_existing_fixture_valid():
    """既存fixture（routine_overridesを持たない証明書）はマイグレーション不要で
    引き続き有効——欠落は必須キー扱いにしない。"""
    profile, certificate, observed_hashes = _validated_profile()
    assert "routine_overrides" not in certificate
    normalized = validate_profile(profile)
    result = validate_certificate(normalized, certificate, observed_hashes)
    assert result["routine_overrides"] == {}


def test_routine_overrides_valid_structure_with_any_and_values_modes_accepted():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {
            "preprocess": {
                "default": {
                    "drift_correct": {"mode": "any"},
                    "normalize": {"mode": "values", "values": ["none", "tic"]},
                },
            },
        },
    })
    normalized = validate_profile(profile)
    result = validate_certificate(normalized, certificate, observed_hashes)
    assert result["routine_overrides"]["preprocess"]["default"]["drift_correct"] == {"mode": "any"}
    assert result["routine_overrides"]["preprocess"]["default"]["normalize"] == {
        "mode": "values", "values": ["none", "tic"],
    }


def test_routine_overrides_unknown_category_rejected():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {"statistics": {}},
    })
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, observed_hashes)


def test_routine_overrides_unknown_preprocess_field_rejected():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {
            "preprocess": {"default": {"policy": {"mode": "any"}}},
        },
    })
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, observed_hashes)


def test_routine_overrides_invalid_mode_rejected():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {
            "preprocess": {"default": {"normalize": {"mode": "everything"}}},
        },
    })
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, observed_hashes)


def test_routine_overrides_values_mode_with_empty_values_rejected():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {
            "preprocess": {"default": {"normalize": {"mode": "values", "values": []}}},
        },
    })
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, observed_hashes)


def test_routine_overrides_any_mode_rejects_extra_keys():
    profile, certificate, observed_hashes = _validated_profile({
        "routine_overrides": {
            "preprocess": {"default": {"normalize": {"mode": "any", "values": ["none"]}}},
        },
    })
    normalized = validate_profile(profile)
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, certificate, observed_hashes)


def test_hand_written_validated_with_mismatched_certificate_rejected_before_allowance_consulted():
    """spec §6裁定の核心: `validation.status='validated'`ラベルだけでは何も
    証明しない。証明書がpermissiveなroutine_overrides（何でも上書き可）を
    持っていても、証明書自体がprofileのpinされたhashと一致しなければ、その
    routine_overridesの中身は一切参照されない（PROFILE_VALIDATION_INVALIDが
    先に発生する）。"""
    data = _base_profile()
    data["validation"] = {
        "status": "validated",
        "scope": data["acquisition"]["scope"],
        "certificate_path": "certificates/fake.json",
        "certificate_sha256": "0" * 64,  # 実在しない証明書と対応しないでたらめなhash
    }
    normalized = validate_profile(data)
    permissive_but_unpinned_certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "profile_content_sha256": profile_content_hash(normalized),
        "dependency_hashes": {},
        "fixed_input_hashes": {},
        "reference_file_hashes": {},
        "validation_output_hashes": {},
        "criteria": [{"criterion_id": "x", "status": "pass", "required": True, "detail": None}],
        "performed_by": "someone",
        "performed_at": "2026-09-15T00:00:00Z",
        "scope": data["acquisition"]["scope"],
        # 全recipe・全fieldを無制限に許可する——中身だけ見れば理想的に見える
        "routine_overrides": {
            "preprocess": {"default": {
                "normalize": {"mode": "any"}, "drift_correct": {"mode": "any"},
                "filter": {"mode": "any"}, "impute": {"mode": "any"},
            }},
        },
    }
    with pytest.raises(DomainError, match="PROFILE_VALIDATION_INVALID"):
        validate_certificate(normalized, permissive_but_unpinned_certificate, {})
