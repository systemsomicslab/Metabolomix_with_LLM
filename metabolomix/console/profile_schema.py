"""`lcms-profile.v1` 契約と `lcms-profile-validation.v1` 証明書の検証（spec §5, §5.1）。

このモジュールは依存グラフの leaf（stdlib と `metabolomix.core.atomic_io` だけに依存する）。
`metabolomix.core.session_state` / `metabolomix.core.mcp_core` / `metabolomix.pipeline.*` /
`metabolomix.tools.*` を import してはいけない——後続タスク（profile読込・request v2・
pipeline stage）がここから import する側であり、逆向きの依存を作ると循環になる。

**このモジュールはファイルを読まない**。method/依存ファイルの実体解決・hash計算・
証明書ファイルの読込は呼び出し側（後続タスク）の責務で、ここでは渡された
dict をそのまま検証するだけ。

## 契約の要点

- `validate_profile(data)`: `lcms-profile.v1` を構造的に検証し、正規化した新しい
  dict を返す。未知キー・型違反・NaN/Infinity・参照整合性違反はすべて
  `DomainError("PROFILE_INVALID", ...)`。
- `profile_content_hash(data)`: `validation` を除いた profile 本体の canonical hash。
  証明書のメタデータ（status を draft→validated に変える等）は本体の内容identityを
  変えない——これが `validation` を除く理由（spec §5「証明書はprofile本体
  （validation部分を除く）のcanonical hash…を必須とする」）。
- `validate_certificate(profile, certificate, observed_hashes)`: 証明書が
  「この profile を正当に validated と名乗らせる」ものかを検証し、全検証に
  通った場合のみ正規化済み証明書（dict）を返す。署名基盤は作らない
  ——検証は (1) 証明書自体の hash が profile に pin された
  `certificate_sha256` と一致するか（証明書の一要素改変で必ず変わる）、
  (2) profile 本体・全依存・固定入力・基準ファイル・検証出力の hash が
  `observed_hashes`（実測値、呼び出し側が計算して渡す）と一致するか、
  (3) 必須基準が全て pass か、(4) 適用範囲（scope）が profile の宣言と
  一致するか、の4点に限る。不一致・不足はすべて
  `DomainError("PROFILE_VALIDATION_INVALID", ...)`。
- 証明書の`routine_overrides`（任意キー。省略時は`{}`＝fail-closed）は
  `execution_purpose="routine"`のpreprocess上書きが許される明示範囲を
  宣言する——`metabolomix.pipeline.request_v2`がこの戻り値を消費して判定する
  （spec §6「証明書の明示許容集合」）。詳細は`docs/schema/lcms-profile-v1.md`
  「証明書 `routine_overrides`」節を参照。

## 「未記載の測定条件」の表現

spec は「測定条件の未記載値は値nullと理由を持たせる」と言うが、その
value/reason の対を実際に保持する場所は `evidence`（field pathキーのdict）で
ある。`acquisition`/`processing`/`software` 自身のフィールドは素の値
（nullable な項目は素の `None`）を持ち、「なぜ null なのか」「どの根拠
tier（`paper_explicit` 等）に基づくか」は `evidence["<field path>"]` 側に
一元化する。この分離により、条件値そのものと出典証跡を別々に更新でき、
未記載理由を書き忘れたまま値だけ埋める事故を「evidenceの`reason`必須」で
機械的に防げる。
"""
from __future__ import annotations

import json
import re

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "SCHEMA",
    "CERTIFICATE_SCHEMA",
    "validate_profile",
    "validate_matrix_recipe",
    "profile_content_hash",
    "validate_certificate",
]

SCHEMA = "lcms-profile.v1"
CERTIFICATE_SCHEMA = "lcms-profile-validation.v1"

_PROFILE_INVALID = "PROFILE_INVALID"
_CERT_INVALID = "PROFILE_VALIDATION_INVALID"

#: 安全なslug。パス脱出・区切り文字を含まない（pipeline.request._SAFE_ID_REと同じ規約）。
_SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
#: evidenceのfield pathキー。ドット区切りの識別子列だけを許す（例:
#: "acquisition.lc.column"）。任意文字列を許すとエラーメッセージ・将来のツール化
#: （field pathからprofileへ書き戻す等）で解釈が割れるため制限する。
_FIELD_PATH_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*$")

_OMICS_VALUES = frozenset({"metabolomics"})
_POLARITY_VALUES = frozenset({"positive", "negative"})
_DEPENDENCY_KINDS = frozenset({"msp", "text_identification", "rt_reference", "lbm"})
_EVIDENCE_TIERS = frozenset({
    "paper_explicit", "paper_cross_reference", "supplied",
    "locally_validated", "proposed",
})
_REQUIRED_EVIDENCE_KINDS = frozenset({"mass_rt", "library_match", "authentic_standard_match"})
_COMPOUND_ID_KEYS = frozenset({
    "name", "inchikey", "formula", "cas", "hmdb_id", "kegg_id", "pubchem_cid",
})
#: nameだけでは同定証拠として不十分（spec §8.1「化合物の識別子（名前だけは不可）」）。
_COMPOUND_ID_NON_NAME_KEYS = _COMPOUND_ID_KEYS - {"name"}
_STATISTIC_KINDS = frozenset({"pca", "welch", "anova_tukey"})
_TRANSFORM_VALUES = frozenset({"none", "log2"})
_QC_METRICS = frozenset({
    "pooled_qc_rsd", "blank_fold", "detection_rate",
    "standard_rt_error", "standard_mass_error_ppm", "internal_standard_valid_fraction",
})
_QC_SCOPE_VALUES = frozenset({"targets", "all_features"})
_QC_OPERATORS = frozenset({"<=", "<", ">=", ">", "=="})
_VALIDATION_STATUS_VALUES = frozenset({"draft", "validated"})
_CRITERION_STATUS_VALUES = frozenset({"pass", "fail", "not_evaluable"})

_TOP_LEVEL_KEYS = frozenset({
    "schema", "profile_id", "revision", "omics", "acquisition", "software",
    "processing", "analysis_recipe", "feature_targets", "matrix_recipes",
    "qc_policy", "evidence", "validation",
})
_ACQUISITION_KEYS = frozenset({
    "separation", "acquisition_type", "polarity", "instrument", "lc",
    "ms_range", "sample_matrix", "scope",
})
_LC_KEYS = frozenset({
    "column", "mobile_phase_a", "mobile_phase_b", "flow_rate_ul_min",
    "column_temperature_c", "gradient_profile",
})
_LC_STR_KEYS = frozenset({"column", "mobile_phase_a", "mobile_phase_b", "gradient_profile"})
_LC_NUM_KEYS = frozenset({"flow_rate_ul_min", "column_temperature_c"})
_MS_RANGE_KEYS = frozenset({"ms1_low_mz", "ms1_high_mz", "ms2_low_mz", "ms2_high_mz"})
_SOFTWARE_KEYS = frozenset({
    "msdial_version", "executable_path", "executable_sha256", "adapter_version",
})
_PROCESSING_KEYS = frozenset({
    "method_path", "method_sha256", "measure", "dependencies", "effective_settings",
})
_DEPENDENCY_KEYS = frozenset({"dependency_id", "kind", "method_key", "path", "sha256", "required"})
_ANALYSIS_RECIPE_KEYS = frozenset({"statistics", "internal_standards"})
_INTERNAL_STANDARD_KEYS = frozenset({"target_id", "standard_target_id"})
_STATISTIC_COMMON_KEYS = frozenset({
    "statistic_id", "kind", "matrix_recipe_id", "transform", "feature_scope",
})
_STATISTIC_EXTRA_KEYS = {
    "pca": frozenset({"scaling", "n_components"}),
    "welch": frozenset({"reference_group", "test_group", "q_threshold", "log2fc_threshold"}),
    "anova_tukey": frozenset({"groups", "alpha"}),
}
_FEATURE_TARGET_KEYS = frozenset({
    "compound_identifiers", "adduct", "charge", "expected_mz", "mz_tolerance_ppm",
    "expected_rt_min", "rt_tolerance_min", "required_evidence",
})
_EVIDENCE_REQUIREMENT_EXTRA_KEYS = {
    "mass_rt": frozenset(),
    "library_match": frozenset({"library_id", "library_sha256", "score_field", "score_threshold"}),
    "authentic_standard_match": frozenset({"mz_tolerance_ppm", "rt_tolerance_min"}),
}
_MATRIX_RECIPE_KEYS = frozenset({"base", "normalize", "drift_correct", "filter", "impute"})
_MATRIX_RECIPE_BASE_VALUES = frozenset({"peak_height", "internal_standard_ratio"})
_MATRIX_RECIPE_NORMALIZE_VALUES = frozenset({"none", "tic", "median", "pqn"})
_MATRIX_RECIPE_IMPUTE_VALUES = frozenset({"none", "half_min", "knn", "column_mean"})
_MATRIX_RECIPE_FILTER_KEYS = frozenset({"min_detection_rate"})
_DEFAULT_RECIPE_ID = "default"
_QC_POLICY_KEYS = frozenset({
    "scope", "target_ids", "required", "threshold", "evidence_requirement",
    "minimum_pass_fraction",
})
_QC_THRESHOLD_KEYS = frozenset({"operator", "value"})
_EVIDENCE_ENTRY_KEYS = frozenset({
    "value", "reason", "tier", "source_uri", "source_hash", "location",
})
_VALIDATION_KEYS = frozenset({"status", "scope", "certificate_path", "certificate_sha256"})

#: 「測定条件の未記載値は値nullと理由を持たせる」（spec §5）の対象となる、
#: nullを取りうる測定条件のfield path一覧。ここに挙げたpathがnullの場合、
#: 対応する`evidence[<path>]`エントリ（reason必須）が無ければ`PROFILE_INVALID`
#: にする——nullを「装置既定値」や「単なる未指定」と解釈させないための
#: 機械的な歯止め（controller裁定: fix round 1 finding 1）。
_NULLABLE_CONDITION_FIELD_PATHS = (
    "acquisition.instrument",
    "acquisition.sample_matrix",
    "acquisition.lc.column",
    "acquisition.lc.mobile_phase_a",
    "acquisition.lc.mobile_phase_b",
    "acquisition.lc.flow_rate_ul_min",
    "acquisition.lc.column_temperature_c",
    "acquisition.lc.gradient_profile",
    "acquisition.ms_range.ms2_low_mz",
    "acquisition.ms_range.ms2_high_mz",
)

_CERTIFICATE_REQUIRED_KEYS = frozenset({
    "schema", "profile_content_sha256", "dependency_hashes", "fixed_input_hashes",
    "reference_file_hashes", "validation_output_hashes", "criteria",
    "performed_by", "performed_at", "scope",
})
#: `routine_overrides`だけが任意（controller裁定: spec §6の「証明書の明示許容
#: 集合」をここに置く。省略はfail-closedで「何も上書きを許可しない」——
#: Task 1が書いた既存fixtureをこのキー無しのまま有効に保つための設計）。
_CERTIFICATE_OPTIONAL_KEYS = frozenset({"routine_overrides"})
_CERTIFICATE_KEYS = _CERTIFICATE_REQUIRED_KEYS | _CERTIFICATE_OPTIONAL_KEYS
_CRITERION_KEYS = frozenset({"criterion_id", "status", "required", "detail"})
_HASH_MAP_CATEGORIES = (
    ("dependency_hashes", "dependencies"),
    ("fixed_input_hashes", "fixed_inputs"),
    ("reference_file_hashes", "reference_files"),
    ("validation_output_hashes", "validation_outputs"),
)

#: `routine_overrides`が宣言できるカテゴリ。現状は`preprocess`（matrix recipe
#: の上書き）だけ——spec §6が例示するもう1つの動機（比較群の変更）は、v2の
#: `statistics`にはまだ「上書き」という概念自体が無い（profileのstatisticsを
#: requestが部分的に上書きするのではなく、request自身が全体を宣言するか
#: profile既定を丸ごと使うかの二択）ため、今回は対象にしない。将来
#: 「statistics」のような別カテゴリを追加する場合もこの集合へ足すだけで済む
#: よう、`routine_overrides`のtop-levelはカテゴリ名をキーにした開いた構造に
#: してある。
_ROUTINE_OVERRIDE_CATEGORIES = frozenset({"preprocess"})
#: `preprocess`カテゴリで許容を宣言できるフィールド
#: （`metabolomix.pipeline.request_v2._PREPROCESS_OVERRIDE_KEYS`と同じ集合。
#: 循環importを避けるためここでも独立に持つ）。
_ROUTINE_OVERRIDE_PREPROCESS_FIELDS = frozenset({"normalize", "drift_correct", "filter", "impute"})
_ROUTINE_OVERRIDE_ALLOWANCE_MODES = frozenset({"any", "values"})


def _fail(code: str, message: str, **details) -> None:
    raise DomainError(code, message, details)


def _fail_profile(message: str, **details) -> None:
    _fail(_PROFILE_INVALID, message, **details)


def _fail_certificate(message: str, **details) -> None:
    _fail(_CERT_INVALID, message, **details)


def _require_json_safe(value: object, label: str) -> None:
    """NaN/Infinity・非JSON型を拒否する。

    `canonical_hash`（`json.dumps(..., allow_nan=False)`）が同じ入力に対して
    素の`ValueError`を投げてしまう前に、ここで`DomainError(PROFILE_INVALID)`
    として検出する。
    """
    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError) as exc:
        _fail_profile(
            f"{label}はJSON直列化可能である必要があります"
            f"（NaN/Infinity・非対応型は不可）: {exc}",
        )


def _require_dict(value: object, label: str) -> dict:
    if not isinstance(value, dict):
        _fail_profile(f"{label}はオブジェクトである必要があります: {value!r}", value=value)
    return value


def _require_list(value: object, label: str) -> list:
    if not isinstance(value, list):
        _fail_profile(f"{label}は配列である必要があります: {value!r}", value=value)
    return value


def _reject_unknown_keys(value: dict, allowed: frozenset, label: str) -> None:
    unknown = set(value) - allowed
    if unknown:
        _fail_profile(f"{label}に未知のキーがあります: {sorted(unknown)}",
                      unknown_keys=sorted(unknown))


def _require_exact_keys(value: dict, allowed: frozenset, label: str) -> None:
    """value のキー集合が allowed とちょうど一致することを要求する。"""
    _reject_unknown_keys(value, allowed, label)
    missing = allowed - set(value)
    if missing:
        _fail_profile(f"{label}に必須キーが不足しています: {sorted(missing)}",
                      missing_keys=sorted(missing))


def _require_slug(value: object, label: str, *, code: str = _PROFILE_INVALID) -> str:
    if not isinstance(value, str) or not _SAFE_ID_RE.fullmatch(value):
        _fail(
            code,
            f"{label}は安全な文字だけのID（英数字・アンダースコア・ハイフン）"
            f"である必要があります: {value!r}",
            value=value,
        )
    return value


def _require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or not _SHA256_RE.fullmatch(value):
        _fail_profile(f"{label}は64桁小文字16進のSHA-256文字列である必要があります: {value!r}",
                      value=value)
    return value


def _require_nonempty_str(value: object, label: str, *, code: str = _PROFILE_INVALID) -> str:
    if not isinstance(value, str) or not value:
        _fail(code, f"{label}は空でない文字列である必要があります: {value!r}", value=value)
    return value


def _require_nullable_str(value: object, label: str) -> str | None:
    if value is not None and (not isinstance(value, str) or value == ""):
        _fail_profile(f"{label}は空でない文字列またはnullである必要があります: {value!r}",
                      value=value)
    return value


def _require_bool(value: object, label: str, *, code: str = _PROFILE_INVALID) -> bool:
    if not isinstance(value, bool):
        _fail(code, f"{label}はboolである必要があります: {value!r}", value=value)
    return value


def _is_finite_number(value: object) -> bool:
    """boolを除いた真の有限数値だけを真にする（NaN/Infinityは既に上流で拒否済みだが、
    ここでは型そのもの——bool・非数値——を拒否する）。
    """
    if isinstance(value, bool):
        return False
    return isinstance(value, (int, float))


def _require_finite_number(value: object, label: str) -> float:
    if not _is_finite_number(value):
        _fail_profile(f"{label}は数値である必要があります: {value!r}", value=value)
    return float(value)


def _require_positive_number(value: object, label: str) -> float:
    number = _require_finite_number(value, label)
    if number <= 0:
        _fail_profile(f"{label}は正の数である必要があります: {value!r}", value=value)
    return number


def _require_positive_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        _fail_profile(f"{label}は正の整数である必要があります（boolは不可）: {value!r}",
                      value=value)
    return value


def _require_probability(value: object, label: str) -> float:
    """0 < x < 1（確率閾値。spec §6.2「確率閾値は0より大きく1未満」）。"""
    number = _require_finite_number(value, label)
    if not (0 < number < 1):
        _fail_profile(f"{label}は0より大きく1未満である必要があります: {value!r}", value=value)
    return number


def _require_unit_fraction(value: object, label: str) -> float:
    """0 <= x <= 1。"""
    number = _require_finite_number(value, label)
    if not (0.0 <= number <= 1.0):
        _fail_profile(f"{label}は0以上1以下である必要があります: {value!r}", value=value)
    return number


# ---------- acquisition ----------

def _validate_lc(value: object) -> dict:
    filled = dict.fromkeys(_LC_KEYS)
    if value is not None:
        _require_dict(value, "acquisition.lc")
        _reject_unknown_keys(value, _LC_KEYS, "acquisition.lc")
        filled.update(value)
    out: dict = {}
    for key in _LC_STR_KEYS:
        out[key] = _require_nullable_str(filled[key], f"acquisition.lc.{key}")
    for key in _LC_NUM_KEYS:
        raw = filled[key]
        out[key] = None if raw is None else _require_positive_number(raw, f"acquisition.lc.{key}")
    return out


def _validate_ms_range(value: object) -> dict:
    _require_dict(value, "acquisition.ms_range")
    _require_exact_keys(value, _MS_RANGE_KEYS, "acquisition.ms_range")

    ms1_low = _require_positive_number(value["ms1_low_mz"], "acquisition.ms_range.ms1_low_mz")
    ms1_high = _require_positive_number(value["ms1_high_mz"], "acquisition.ms_range.ms1_high_mz")
    if ms1_high <= ms1_low:
        _fail_profile("acquisition.ms_range.ms1_high_mzはms1_low_mzより大きい必要があります。",
                      ms1_low_mz=ms1_low, ms1_high_mz=ms1_high)

    ms2_low, ms2_high = value["ms2_low_mz"], value["ms2_high_mz"]
    if (ms2_low is None) != (ms2_high is None):
        _fail_profile("acquisition.ms_range.ms2_low_mz/ms2_high_mzはnullを片方だけには"
                      "できません。", ms2_low_mz=ms2_low, ms2_high_mz=ms2_high)
    if ms2_low is not None:
        ms2_low = _require_positive_number(ms2_low, "acquisition.ms_range.ms2_low_mz")
        ms2_high = _require_positive_number(ms2_high, "acquisition.ms_range.ms2_high_mz")
        if ms2_high <= ms2_low:
            _fail_profile(
                "acquisition.ms_range.ms2_high_mzはms2_low_mzより大きい必要があります。",
                ms2_low_mz=ms2_low, ms2_high_mz=ms2_high)

    return {
        "ms1_low_mz": ms1_low, "ms1_high_mz": ms1_high,
        "ms2_low_mz": ms2_low, "ms2_high_mz": ms2_high,
    }


def _validate_acquisition(value: object) -> dict:
    _require_dict(value, "acquisition")
    _require_exact_keys(value, _ACQUISITION_KEYS, "acquisition")

    if value["separation"] != "lc":
        _fail_profile(f"acquisition.separationは'lc'固定です: {value['separation']!r}")
    if value["acquisition_type"] != "dda":
        _fail_profile(
            f"acquisition.acquisition_typeは'dda'固定です: {value['acquisition_type']!r}")
    polarity = value["polarity"]
    if polarity not in _POLARITY_VALUES:
        _fail_profile(f"acquisition.polarityが不正です: {polarity!r}", polarity=polarity)

    return {
        "separation": "lc",
        "acquisition_type": "dda",
        "polarity": polarity,
        "instrument": _require_nullable_str(value["instrument"], "acquisition.instrument"),
        "lc": _validate_lc(value["lc"]),
        "ms_range": _validate_ms_range(value["ms_range"]),
        "sample_matrix": _require_nullable_str(value["sample_matrix"], "acquisition.sample_matrix"),
        "scope": _require_nonempty_str(value["scope"], "acquisition.scope"),
    }


# ---------- software / processing ----------

def _validate_software(value: object) -> dict:
    _require_dict(value, "software")
    _require_exact_keys(value, _SOFTWARE_KEYS, "software")
    return {
        "msdial_version": _require_nonempty_str(value["msdial_version"], "software.msdial_version"),
        "executable_path": _require_nonempty_str(value["executable_path"], "software.executable_path"),
        "executable_sha256": _require_sha256(value["executable_sha256"], "software.executable_sha256"),
        "adapter_version": _require_nonempty_str(value["adapter_version"], "software.adapter_version"),
    }


def _validate_dependency(value: object, index: int, seen_ids: set) -> dict:
    _require_dict(value, f"processing.dependencies[{index}]")
    _require_exact_keys(value, _DEPENDENCY_KEYS, f"processing.dependencies[{index}]")

    dependency_id = _require_slug(value["dependency_id"], f"processing.dependencies[{index}].dependency_id")
    if dependency_id in seen_ids:
        _fail_profile(f"dependency_idが重複しています: {dependency_id!r}",
                      dependency_id=dependency_id)
    seen_ids.add(dependency_id)

    kind = value["kind"]
    if kind not in _DEPENDENCY_KINDS:
        _fail_profile(f"processing.dependencies[{index}].kindが不正です: {kind!r}",
                      index=index, kind=kind)

    return {
        "dependency_id": dependency_id,
        "kind": kind,
        "method_key": _require_nonempty_str(value["method_key"], f"processing.dependencies[{index}].method_key"),
        "path": _require_nonempty_str(value["path"], f"processing.dependencies[{index}].path"),
        "sha256": _require_sha256(value["sha256"], f"processing.dependencies[{index}].sha256"),
        "required": _require_bool(value["required"], f"processing.dependencies[{index}].required"),
    }


def _validate_processing(value: object) -> dict:
    _require_dict(value, "processing")
    _require_exact_keys(value, _PROCESSING_KEYS, "processing")

    if value["measure"] != "peak_height":
        _fail_profile(f"processing.measureは'peak_height'固定です: {value['measure']!r}")

    dependencies_raw = _require_list(value["dependencies"], "processing.dependencies")
    seen_ids: set = set()
    dependencies = [_validate_dependency(item, i, seen_ids) for i, item in enumerate(dependencies_raw)]

    effective_settings = _require_dict(value["effective_settings"], "processing.effective_settings")

    return {
        "method_path": _require_nonempty_str(value["method_path"], "processing.method_path"),
        "method_sha256": _require_sha256(value["method_sha256"], "processing.method_sha256"),
        "measure": "peak_height",
        "dependencies": dependencies,
        # 実効設定snapshotは開いた辞書（MS-DIALの実パラメータキー集合を固定しない）。
        # NaN/Infinity・非JSON型の拒否は validate_profile 冒頭の _require_json_safe が
        # profile全体に対してすでに行っている。
        "effective_settings": dict(effective_settings),
    }


# ---------- feature_targets ----------

def _require_nonzero_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value == 0:
        _fail_profile(f"{label}は0でない整数である必要があります（boolは不可）: {value!r}",
                      value=value)
    return value


def _validate_compound_identifiers(value: object, label: str) -> dict:
    _require_dict(value, label)
    _reject_unknown_keys(value, _COMPOUND_ID_KEYS, label)
    out: dict = {}
    has_non_name = False
    for key in _COMPOUND_ID_KEYS:
        if key not in value:
            continue
        text = _require_nonempty_str(value[key], f"{label}.{key}")
        out[key] = text
        if key in _COMPOUND_ID_NON_NAME_KEYS:
            has_non_name = True
    if not has_non_name:
        _fail_profile(
            f"{label}は名前以外の識別子（{sorted(_COMPOUND_ID_NON_NAME_KEYS)}のいずれか）"
            "を少なくとも1つ持つ必要があります（名前だけは不可）。",
            provided=sorted(out),
        )
    return out


def _validate_required_evidence_item(value: object, index: int, label: str) -> dict:
    item_label = f"{label}[{index}]"
    _require_dict(value, item_label)
    kind = value.get("kind")
    if kind not in _REQUIRED_EVIDENCE_KINDS:
        _fail_profile(f"{item_label}.kindが不正です: {kind!r}", index=index, kind=kind)

    allowed = frozenset({"kind"}) | _EVIDENCE_REQUIREMENT_EXTRA_KEYS[kind]
    _require_exact_keys(value, allowed, item_label)

    if kind == "mass_rt":
        return {"kind": "mass_rt"}
    if kind == "library_match":
        return {
            "kind": "library_match",
            "library_id": _require_nonempty_str(value["library_id"], f"{item_label}.library_id"),
            "library_sha256": _require_sha256(value["library_sha256"], f"{item_label}.library_sha256"),
            "score_field": _require_nonempty_str(value["score_field"], f"{item_label}.score_field"),
            "score_threshold": _require_finite_number(value["score_threshold"], f"{item_label}.score_threshold"),
        }
    # authentic_standard_match
    return {
        "kind": "authentic_standard_match",
        "mz_tolerance_ppm": _require_positive_number(value["mz_tolerance_ppm"], f"{item_label}.mz_tolerance_ppm"),
        "rt_tolerance_min": _require_positive_number(value["rt_tolerance_min"], f"{item_label}.rt_tolerance_min"),
    }


def _validate_feature_target(value: object, target_id: str) -> dict:
    label = f"feature_targets[{target_id!r}]"
    _require_dict(value, label)
    _require_exact_keys(value, _FEATURE_TARGET_KEYS, label)

    required_evidence_raw = _require_list(value["required_evidence"], f"{label}.required_evidence")
    if not required_evidence_raw:
        _fail_profile(f"{label}.required_evidenceは空にできません。")
    required_evidence = [
        _validate_required_evidence_item(item, i, f"{label}.required_evidence")
        for i, item in enumerate(required_evidence_raw)
    ]
    seen_kinds = {item["kind"] for item in required_evidence}
    if len(seen_kinds) != len(required_evidence):
        _fail_profile(f"{label}.required_evidenceに同じkindが重複しています。")

    return {
        "compound_identifiers": _validate_compound_identifiers(
            value["compound_identifiers"], f"{label}.compound_identifiers"),
        "adduct": _require_nonempty_str(value["adduct"], f"{label}.adduct"),
        "charge": _require_nonzero_int(value["charge"], f"{label}.charge"),
        "expected_mz": _require_positive_number(value["expected_mz"], f"{label}.expected_mz"),
        "mz_tolerance_ppm": _require_positive_number(value["mz_tolerance_ppm"], f"{label}.mz_tolerance_ppm"),
        "expected_rt_min": _require_positive_number(value["expected_rt_min"], f"{label}.expected_rt_min"),
        "rt_tolerance_min": _require_positive_number(value["rt_tolerance_min"], f"{label}.rt_tolerance_min"),
        "required_evidence": required_evidence,
    }


def _validate_feature_targets(value: object) -> dict:
    _require_dict(value, "feature_targets")
    out: dict = {}
    for target_id, target_value in value.items():
        _require_slug(target_id, "feature_targetsのキー(target_id)")
        out[target_id] = _validate_feature_target(target_value, target_id)
    return out


# ---------- matrix_recipes ----------

def _validate_matrix_recipe_filter(value: object, label: str) -> dict | None:
    if value is None:
        return None
    _require_dict(value, label)
    _require_exact_keys(value, _MATRIX_RECIPE_FILTER_KEYS, label)
    return {
        "min_detection_rate": _require_unit_fraction(
            value["min_detection_rate"], f"{label}.min_detection_rate"),
    }


def _validate_matrix_recipe(value: object, recipe_id: str) -> dict:
    label = f"matrix_recipes[{recipe_id!r}]"
    _require_dict(value, label)
    _require_exact_keys(value, _MATRIX_RECIPE_KEYS, label)

    base = value["base"]
    if base not in _MATRIX_RECIPE_BASE_VALUES:
        _fail_profile(f"{label}.baseが不正です: {base!r}", base=base)

    normalize = value["normalize"]
    if normalize not in _MATRIX_RECIPE_NORMALIZE_VALUES:
        _fail_profile(f"{label}.normalizeが不正です: {normalize!r}", normalize=normalize)
    if base == "internal_standard_ratio" and normalize != "none":
        # spec §8「内部標準比viewへのTIC/median/PQNは初期版では拒否する。
        # 二重正規化を暗黙に実行しない。」
        _fail_profile(
            f"{label}: base='internal_standard_ratio'にnormalize={normalize!r}は"
            "許可されません（二重正規化の禁止）。", base=base, normalize=normalize)

    impute = value["impute"]
    if impute not in _MATRIX_RECIPE_IMPUTE_VALUES:
        _fail_profile(f"{label}.imputeが不正です: {impute!r}", impute=impute)

    return {
        "base": base,
        "normalize": normalize,
        "drift_correct": _require_bool(value["drift_correct"], f"{label}.drift_correct"),
        "filter": _validate_matrix_recipe_filter(value["filter"], f"{label}.filter"),
        "impute": impute,
    }


def validate_matrix_recipe(recipe: object, recipe_id: str = _DEFAULT_RECIPE_ID) -> dict:
    """matrix recipe 1件を検証して正規化する（profileの外からも使う公開API）。

    profile の `matrix_recipes` を検証するのと**同じ規則**をそのまま使う——
    単体ツールから recipe を渡す経路（`dataset_build_matrix`）が、profile 経由
    より緩い recipe を通してしまうと、同じ recipe_id の行列が経路によって
    別物になる。戻り値に `recipe_id` は含めない（呼出側が付ける）。
    """
    return _validate_matrix_recipe(recipe, recipe_id)


def _validate_matrix_recipes(value: object) -> dict:
    _require_dict(value, "matrix_recipes")
    out: dict = {}
    for recipe_id, recipe_value in value.items():
        _require_slug(recipe_id, "matrix_recipesのキー(recipe_id)")
        out[recipe_id] = _validate_matrix_recipe(recipe_value, recipe_id)
    if _DEFAULT_RECIPE_ID not in out:
        _fail_profile(f"matrix_recipesに既定recipe {_DEFAULT_RECIPE_ID!r} がありません。")
    return out


# ---------- analysis_recipe ----------

def _validate_feature_scope(value: object, label: str, feature_target_ids: frozenset) -> dict:
    _require_dict(value, label)
    mode = value.get("mode")
    if mode == "all_eligible":
        _require_exact_keys(value, frozenset({"mode"}), label)
        return {"mode": "all_eligible"}
    if mode == "targets":
        _require_exact_keys(value, frozenset({"mode", "target_ids"}), label)
        target_ids = _require_list(value["target_ids"], f"{label}.target_ids")
        if not target_ids:
            _fail_profile(f"{label}.target_idsは空にできません。")
        if len(set(target_ids)) != len(target_ids):
            _fail_profile(f"{label}.target_idsに重複があります。", target_ids=target_ids)
        for target_id in target_ids:
            if target_id not in feature_target_ids:
                _fail_profile(f"{label}.target_idsが存在しないtarget_idを参照しています: "
                              f"{target_id!r}", target_id=target_id)
        return {"mode": "targets", "target_ids": list(target_ids)}
    _fail_profile(f"{label}.modeが不正です: {mode!r}", mode=mode)


def _validate_statistic(
    value: object, index: int, matrix_recipe_ids: frozenset, feature_target_ids: frozenset,
) -> dict:
    label = f"analysis_recipe.statistics[{index}]"
    _require_dict(value, label)
    missing_common = _STATISTIC_COMMON_KEYS - set(value)
    if missing_common:
        _fail_profile(f"{label}に必須キーが不足しています: {sorted(missing_common)}",
                      missing_keys=sorted(missing_common))
    kind = value["kind"]
    if kind not in _STATISTIC_KINDS:
        _fail_profile(f"{label}.kindが不正です: {kind!r}", index=index, kind=kind)

    allowed = _STATISTIC_COMMON_KEYS | _STATISTIC_EXTRA_KEYS[kind]
    _require_exact_keys(value, allowed, label)

    statistic_id = _require_slug(value["statistic_id"], f"{label}.statistic_id")

    matrix_recipe_id = value["matrix_recipe_id"]
    if matrix_recipe_id not in matrix_recipe_ids:
        _fail_profile(f"{label}.matrix_recipe_idが存在しないrecipeを参照しています: "
                      f"{matrix_recipe_id!r}", matrix_recipe_id=matrix_recipe_id)

    transform = value["transform"]
    if transform not in _TRANSFORM_VALUES:
        _fail_profile(f"{label}.transformが不正です: {transform!r}", transform=transform)

    feature_scope = _validate_feature_scope(value["feature_scope"], f"{label}.feature_scope",
                                            feature_target_ids)

    out = {
        "statistic_id": statistic_id,
        "kind": kind,
        "matrix_recipe_id": matrix_recipe_id,
        "transform": transform,
        "feature_scope": feature_scope,
    }

    if kind == "pca":
        scaling = value["scaling"]
        if scaling not in ("none", "autoscale"):
            _fail_profile(f"{label}.scalingが不正です: {scaling!r}", scaling=scaling)
        out["scaling"] = scaling
        out["n_components"] = _require_positive_int(value["n_components"], f"{label}.n_components")
    elif kind == "welch":
        reference_group = _require_nonempty_str(value["reference_group"], f"{label}.reference_group")
        test_group = _require_nonempty_str(value["test_group"], f"{label}.test_group")
        if reference_group == test_group:
            _fail_profile(f"{label}: reference_groupとtest_groupが同一です: {reference_group!r}")
        out["reference_group"] = reference_group
        out["test_group"] = test_group
        out["q_threshold"] = _require_probability(value["q_threshold"], f"{label}.q_threshold")
        out["log2fc_threshold"] = _require_finite_number(
            value["log2fc_threshold"], f"{label}.log2fc_threshold")
        if out["log2fc_threshold"] < 0:
            _fail_profile(f"{label}.log2fc_thresholdは非負である必要があります: "
                          f"{value['log2fc_threshold']!r}")
    else:  # anova_tukey
        groups = _require_list(value["groups"], f"{label}.groups")
        if len(groups) < 3:
            _fail_profile(f"{label}.groupsは3群以上である必要があります: {groups!r}")
        if len(set(groups)) != len(groups):
            _fail_profile(f"{label}.groupsに重複があります: {groups!r}")
        for group in groups:
            _require_nonempty_str(group, f"{label}.groups[]")
        out["groups"] = list(groups)
        out["alpha"] = _require_probability(value["alpha"], f"{label}.alpha")

    return out


def _validate_statistics(value: object, matrix_recipe_ids: frozenset, feature_target_ids: frozenset) -> list:
    items = _require_list(value, "analysis_recipe.statistics")
    seen_ids: set = set()
    out = []
    for index, item in enumerate(items):
        normalized = _validate_statistic(item, index, matrix_recipe_ids, feature_target_ids)
        if normalized["statistic_id"] in seen_ids:
            _fail_profile(
                f"analysis_recipe.statistics[{index}].statistic_idが重複しています: "
                f"{normalized['statistic_id']!r}", statistic_id=normalized["statistic_id"])
        seen_ids.add(normalized["statistic_id"])
        out.append(normalized)
    return out


def _validate_internal_standards(value: object, feature_target_ids: frozenset) -> list:
    items = _require_list(value, "analysis_recipe.internal_standards")
    normalized = []
    edges: dict[str, str] = {}
    for index, item in enumerate(items):
        label = f"analysis_recipe.internal_standards[{index}]"
        _require_dict(item, label)
        _require_exact_keys(item, _INTERNAL_STANDARD_KEYS, label)

        target_id = item["target_id"]
        standard_target_id = item["standard_target_id"]
        for role, tid in (("target_id", target_id), ("standard_target_id", standard_target_id)):
            if tid not in feature_target_ids:
                _fail_profile(f"{label}.{role}が存在しないtarget_idを参照しています: {tid!r}",
                              target_id=tid)
        if target_id == standard_target_id:
            _fail_profile(f"{label}: target_idとstandard_target_idが同一です（自己参照）: "
                          f"{target_id!r}")
        if target_id in edges:
            _fail_profile(f"{label}: target_id {target_id!r} の内部標準対応が重複しています。",
                          target_id=target_id)
        edges[target_id] = standard_target_id
        normalized.append({"target_id": target_id, "standard_target_id": standard_target_id})

    for start in edges:
        visited = set()
        current = start
        while current in edges:
            if current in visited:
                _fail_profile(f"内部標準対応が循環しています: {start!r} から到達する対応が"
                              "自身に戻ります。", start=start)
            visited.add(current)
            current = edges[current]

    return normalized


def _validate_analysis_recipe(value: object, matrix_recipe_ids: frozenset, feature_target_ids: frozenset) -> dict:
    _require_dict(value, "analysis_recipe")
    _require_exact_keys(value, _ANALYSIS_RECIPE_KEYS, "analysis_recipe")
    return {
        "statistics": _validate_statistics(value["statistics"], matrix_recipe_ids, feature_target_ids),
        "internal_standards": _validate_internal_standards(value["internal_standards"], feature_target_ids),
    }


# ---------- qc_policy ----------

def _validate_qc_threshold(value: object, label: str) -> dict:
    _require_dict(value, label)
    _require_exact_keys(value, _QC_THRESHOLD_KEYS, label)
    operator = value["operator"]
    if operator not in _QC_OPERATORS:
        _fail_profile(f"{label}.operatorが不正です: {operator!r}", operator=operator)
    return {
        "operator": operator,
        "value": _require_finite_number(value["value"], f"{label}.value"),
    }


def _validate_qc_entry(value: object, metric: str, feature_target_ids: frozenset) -> dict:
    label = f"qc_policy[{metric!r}]"
    _require_dict(value, label)
    _require_exact_keys(value, _QC_POLICY_KEYS, label)

    scope = value["scope"]
    if scope not in _QC_SCOPE_VALUES:
        _fail_profile(f"{label}.scopeが不正です: {scope!r}", scope=scope)

    target_ids = value["target_ids"]
    if scope == "targets":
        target_ids = _require_list(target_ids, f"{label}.target_ids")
        if not target_ids:
            _fail_profile(f"{label}.target_idsは空にできません（scope='targets'）。")
        if len(set(target_ids)) != len(target_ids):
            _fail_profile(f"{label}.target_idsに重複があります。", target_ids=target_ids)
        for target_id in target_ids:
            if target_id not in feature_target_ids:
                _fail_profile(f"{label}.target_idsが存在しないtarget_idを参照しています: "
                              f"{target_id!r}", target_id=target_id)
        target_ids = list(target_ids)
    else:
        if target_ids is not None:
            _fail_profile(f"{label}.target_idsはscope='all_features'ではnullである必要が"
                          f"あります: {target_ids!r}")

    minimum_pass_fraction = value["minimum_pass_fraction"]
    if minimum_pass_fraction is not None:
        minimum_pass_fraction = _require_unit_fraction(
            minimum_pass_fraction, f"{label}.minimum_pass_fraction")

    return {
        "scope": scope,
        "target_ids": target_ids,
        "required": _require_bool(value["required"], f"{label}.required"),
        "threshold": _validate_qc_threshold(value["threshold"], f"{label}.threshold"),
        "evidence_requirement": _require_bool(value["evidence_requirement"], f"{label}.evidence_requirement"),
        "minimum_pass_fraction": minimum_pass_fraction,
    }


def _validate_qc_policy(value: object, feature_target_ids: frozenset) -> dict:
    _require_dict(value, "qc_policy")
    out: dict = {}
    for metric, entry in value.items():
        if metric not in _QC_METRICS:
            _fail_profile(f"qc_policyに未知のmetricがあります: {metric!r}", metric=metric)
        out[metric] = _validate_qc_entry(entry, metric, feature_target_ids)
    return out


# ---------- evidence ----------

_MISSING = object()  # 「keyが無い」と「値がNone」を区別するための番兵


def _resolve_field_path(profile_so_far: dict, field_path: str):
    """ドット区切りのfield pathをprofile内の実際の値へ解決する。

    辞書のキーだけを辿る（`_FIELD_PATH_RE`が角括弧・添字を許さないため、
    リスト要素へは到達できない——それでよい。evidenceが記述する対象は
    固定schemaのスカラー値であり、`processing.dependencies[]`のような
    リストの特定要素を指す用途は想定しない）。

    戻り値は解決できた場合は実際の値、できなければ`_MISSING`番兵。
    """
    current: object = profile_so_far
    for part in field_path.split("."):
        if not isinstance(current, dict) or part not in current:
            return _MISSING
        current = current[part]
    return current


def _validate_evidence_entry(value: object, field_path: str) -> dict:
    label = f"evidence[{field_path!r}]"
    _require_dict(value, label)
    _require_exact_keys(value, _EVIDENCE_ENTRY_KEYS, label)

    tier = value["tier"]
    if tier not in _EVIDENCE_TIERS:
        _fail_profile(f"{label}.tierが不正です: {tier!r}", tier=tier)

    # NaN/Infinity・非JSON型は validate_profile 冒頭の _require_json_safe が
    # profile全体（このentry_valueを含む）に対してすでに拒否している。
    entry_value = value["value"]
    reason = value["reason"]
    if entry_value is None:
        reason = _require_nonempty_str(reason, f"{label}.reason")
    else:
        reason = _require_nullable_str(reason, f"{label}.reason")

    source_hash = value["source_hash"]
    if source_hash is not None:
        source_hash = _require_sha256(source_hash, f"{label}.source_hash")

    return {
        "value": entry_value,
        "reason": reason,
        "tier": tier,
        "source_uri": _require_nullable_str(value["source_uri"], f"{label}.source_uri"),
        "source_hash": source_hash,
        "location": _require_nullable_str(value["location"], f"{label}.location"),
    }


def _validate_evidence(value: object, profile_so_far: dict) -> dict:
    """evidenceを検証する。

    各キーは`profile_so_far`（`evidence`/`validation`を除く、検証済みの
    profile本体）内の実在するfield pathへ解決できなければならず（controller
    裁定: fix round 1 finding 2）、解決できた場合は`entry["value"]`がそこに
    実際にある値と一致しなければならない——evidenceが「実際のprofile内容とは
    無関係な自己申告」になることを防ぎ、証明書が守る監査証跡としての意味を
    保つ。
    """
    _require_dict(value, "evidence")
    out: dict = {}
    for field_path, entry in value.items():
        if not isinstance(field_path, str) or not _FIELD_PATH_RE.fullmatch(field_path):
            _fail_profile(f"evidenceのキー(field path)が不正です: {field_path!r}",
                          field_path=field_path)
        actual_value = _resolve_field_path(profile_so_far, field_path)
        if actual_value is _MISSING:
            _fail_profile(
                f"evidence[{field_path!r}]が参照するfield pathがprofile内に"
                "存在しません。", field_path=field_path)
        normalized_entry = _validate_evidence_entry(entry, field_path)
        if normalized_entry["value"] != actual_value:
            _fail_profile(
                f"evidence[{field_path!r}].valueが実際の値と一致しません: "
                f"evidence={normalized_entry['value']!r} actual={actual_value!r}",
                field_path=field_path, evidence_value=normalized_entry["value"],
                actual_value=actual_value)
        out[field_path] = normalized_entry
    return out


def _require_evidence_for_null_conditions(profile_so_far: dict, evidence: dict) -> None:
    """nullな測定条件には対応するevidenceエントリの存在を必須にする。

    `_validate_evidence_entry`は「evidenceエントリが存在する場合」その
    value/reasonの整合を検証するが、そもそもエントリ自体が無い場合を
    ここで検出する（spec §5「測定条件の未記載値は値nullと理由を持たせる。
    nullを装置既定値と解釈しない」・controller裁定: fix round 1 finding 1）。
    """
    for field_path in _NULLABLE_CONDITION_FIELD_PATHS:
        if _resolve_field_path(profile_so_far, field_path) is None and field_path not in evidence:
            _fail_profile(
                f"測定条件 {field_path!r} がnullですが、対応する"
                f"evidence[{field_path!r}]がありません"
                "（未記載値には理由付きのevidenceが必須です）。",
                field_path=field_path)


# ---------- validation ----------

def _validate_validation(value: object) -> dict:
    _require_dict(value, "validation")
    _require_exact_keys(value, _VALIDATION_KEYS, "validation")

    status = value["status"]
    if status not in _VALIDATION_STATUS_VALUES:
        _fail_profile(f"validation.statusが不正です: {status!r}", status=status)

    scope = _require_nonempty_str(value["scope"], "validation.scope")

    certificate_path = value["certificate_path"]
    certificate_sha256 = value["certificate_sha256"]
    if status == "draft":
        if certificate_path is not None or certificate_sha256 is not None:
            _fail_profile(
                "validation.status='draft'ではcertificate_path/certificate_sha256は"
                "nullである必要があります（未検証のprofileに証明書ポインタを持たせない）。",
                certificate_path=certificate_path, certificate_sha256=certificate_sha256)
    else:  # validated
        certificate_path = _require_nonempty_str(certificate_path, "validation.certificate_path")
        certificate_sha256 = _require_sha256(certificate_sha256, "validation.certificate_sha256")

    return {
        "status": status,
        "scope": scope,
        "certificate_path": certificate_path,
        "certificate_sha256": certificate_sha256,
    }


# ---------- 公開API ----------

def validate_profile(data: dict) -> dict:
    """`lcms-profile.v1` を検証し、正規化した新しいdictを返す。

    未知キー・型違反・NaN/Infinity・参照整合性違反（存在しないtarget_id/recipe_id
    参照、内部標準対応の重複・循環、既定matrix recipeの欠落等）はすべて
    `DomainError("PROFILE_INVALID", ...)`。
    """
    _require_dict(data, "profile")
    _require_json_safe(data, "profile")
    _require_exact_keys(data, _TOP_LEVEL_KEYS, "profile")

    if data["schema"] != SCHEMA:
        _fail_profile(f"schemaは{SCHEMA!r}のみ許可されます: {data['schema']!r}")

    profile_id = _require_slug(data["profile_id"], "profile_id")
    revision = _require_positive_int(data["revision"], "revision")

    omics = data["omics"]
    if omics not in _OMICS_VALUES:
        _fail_profile(f"omicsが不正です: {omics!r}", omics=omics)

    acquisition = _validate_acquisition(data["acquisition"])
    software = _validate_software(data["software"])
    processing = _validate_processing(data["processing"])
    feature_targets = _validate_feature_targets(data["feature_targets"])
    matrix_recipes = _validate_matrix_recipes(data["matrix_recipes"])
    feature_target_ids = frozenset(feature_targets)
    matrix_recipe_ids = frozenset(matrix_recipes)
    analysis_recipe = _validate_analysis_recipe(data["analysis_recipe"], matrix_recipe_ids, feature_target_ids)
    qc_policy = _validate_qc_policy(data["qc_policy"], feature_target_ids)

    # evidenceは「evidenceを除くprofile本体」に対して検証する
    # （参照先の実在・実際の値との一致、そしてnull条件に対応エントリが
    # あることをここまでに検証済みの本体を土台に確認する）。
    profile_so_far = {
        "schema": SCHEMA, "profile_id": profile_id, "revision": revision, "omics": omics,
        "acquisition": acquisition, "software": software, "processing": processing,
        "analysis_recipe": analysis_recipe, "feature_targets": feature_targets,
        "matrix_recipes": matrix_recipes, "qc_policy": qc_policy,
    }
    evidence = _validate_evidence(data["evidence"], profile_so_far)
    _require_evidence_for_null_conditions(profile_so_far, evidence)
    validation = _validate_validation(data["validation"])

    return {
        "schema": SCHEMA,
        "profile_id": profile_id,
        "revision": revision,
        "omics": omics,
        "acquisition": acquisition,
        "software": software,
        "processing": processing,
        "analysis_recipe": analysis_recipe,
        "feature_targets": feature_targets,
        "matrix_recipes": matrix_recipes,
        "qc_policy": qc_policy,
        "evidence": evidence,
        "validation": validation,
    }


def profile_content_hash(data: dict) -> str:
    """`validation`を除いたprofile本体のcanonical hash。

    証明書の付与・status変更（draft→validated等）はprofile本体の内容identityを
    変えない——これが`validation`を除く理由（spec §5）。
    """
    return canonical_hash({k: v for k, v in data.items() if k != "validation"})


def _validate_hash_map(value: object, label: str) -> dict:
    _require_dict(value, label)
    out = {}
    for key, digest in value.items():
        if not isinstance(key, str) or not key:
            _fail_certificate(f"{label}のキーは空でない文字列である必要があります: {key!r}")
        if not isinstance(digest, str) or not _SHA256_RE.fullmatch(digest):
            _fail_certificate(f"{label}[{key!r}]は64桁小文字16進のSHA-256文字列である必要が"
                              f"あります: {digest!r}", key=key, value=digest)
        out[key] = digest
    return out


def _validate_criterion(value: object, index: int) -> dict:
    label = f"criteria[{index}]"
    _require_dict(value, label)
    _require_exact_keys(value, _CRITERION_KEYS, label)
    status = value["status"]
    if status not in _CRITERION_STATUS_VALUES:
        _fail_certificate(f"{label}.statusが不正です: {status!r}", status=status)
    required = _require_bool(value["required"], f"{label}.required", code=_CERT_INVALID)
    detail = value["detail"]
    if detail is not None and not isinstance(detail, str):
        _fail_certificate(f"{label}.detailは文字列またはnullである必要があります: {detail!r}")
    criterion_id = _require_slug(value["criterion_id"], f"{label}.criterion_id", code=_CERT_INVALID)
    return {
        "criterion_id": criterion_id,
        "status": status,
        "required": required,
        "detail": detail,
    }


# ---------- routine_overrides（spec §6の「証明書の明示許容集合」） ----------

def _validate_routine_override_allowance(value: object, label: str) -> dict:
    """1つのrecipe_id×fieldに対する許容を検証する。

    `{"mode": "any"}`（任意の値を許可）または`{"mode": "values", "values": [...]}`
    （列挙した値だけを許可）のどちらか——構造から範囲を推測する第3の道は無い
    （controller裁定: 「explicit only」）。

    ``_require_dict``/``_require_list``/``_require_exact_keys``（profile本体
    検証用の共有helper）はここでは使わない——それらは無条件で
    ``PROFILE_INVALID``を送出するため、証明書の一部であるこの検証が
    ``PROFILE_VALIDATION_INVALID``を返すべき契約と食い違う
    （``_validate_certificate_shape``自身が同じ理由でこれらのhelperを
    避け、手書きの検査＋``_fail_certificate``で統一しているのと同じ規約）。
    """
    if not isinstance(value, dict):
        _fail_certificate(f"{label}はオブジェクトである必要があります: {value!r}", value=value)
    mode = value.get("mode")
    if mode == "any":
        unknown = set(value) - {"mode"}
        if unknown:
            _fail_certificate(f"{label}に未知のキーがあります: {sorted(unknown)}",
                              unknown_keys=sorted(unknown))
        return {"mode": "any"}
    if mode == "values":
        unknown = set(value) - {"mode", "values"}
        if unknown:
            _fail_certificate(f"{label}に未知のキーがあります: {sorted(unknown)}",
                              unknown_keys=sorted(unknown))
        if "values" not in value:
            _fail_certificate(f"{label}に必須キーが不足しています: ['values']",
                              missing_keys=["values"])
        values = value["values"]
        if not isinstance(values, list):
            _fail_certificate(f"{label}.valuesは配列である必要があります: {values!r}", value=values)
        if not values:
            _fail_certificate(f"{label}.valuesは空にできません。")
        return {"mode": "values", "values": list(values)}
    _fail_certificate(f"{label}.modeが不正です（{sorted(_ROUTINE_OVERRIDE_ALLOWANCE_MODES)}の"
                      f"いずれか）: {mode!r}", mode=mode)


def _validate_routine_overrides_preprocess(value: object, label: str) -> dict:
    if not isinstance(value, dict):
        _fail_certificate(f"{label}はオブジェクトである必要があります: {value!r}", value=value)
    out: dict = {}
    for recipe_id, fields in value.items():
        _require_slug(recipe_id, f"{label}のキー(recipe_id)", code=_CERT_INVALID)
        recipe_label = f"{label}[{recipe_id!r}]"
        if not isinstance(fields, dict):
            _fail_certificate(f"{recipe_label}はオブジェクトである必要があります: {fields!r}",
                              value=fields)
        unknown_fields = set(fields) - _ROUTINE_OVERRIDE_PREPROCESS_FIELDS
        if unknown_fields:
            _fail_certificate(
                f"{recipe_label}に未知のフィールドがあります: {sorted(unknown_fields)}",
                unknown_keys=sorted(unknown_fields),
            )
        out[recipe_id] = {
            field: _validate_routine_override_allowance(allowance, f"{recipe_label}.{field}")
            for field, allowance in fields.items()
        }
    return out


def _validate_routine_overrides(value: object) -> dict:
    """spec §6「routineで許される...変更等は証明書の適用範囲内に限定する」・
    このbriefの裁定「証明書が明示許容集合を持つ」を実装する。

    キー自体が無ければ`{}`——**fail-closed**: `execution_purpose="routine"`
    では何も上書きを許可しない、という既定になる（Task 1が書いた既存の
    証明書fixtureに、このキーを足すマイグレーションを要求しない）。宣言
    されたrecipe_id×field以外はrequest_v2側で常に拒否される——ここでの
    検証は構造だけを保証し、「profile本体に既にある値だから許容する」的な
    推測は一切行わない（controller裁定: 前回のrequest_v2独自ルールは
    このため差し戻された）。
    """
    if value is None:
        return {}
    if not isinstance(value, dict):
        _fail_certificate(f"certificate.routine_overridesはオブジェクトである必要があります: "
                          f"{value!r}", value=value)
    unknown = set(value) - _ROUTINE_OVERRIDE_CATEGORIES
    if unknown:
        _fail_certificate(
            f"certificate.routine_overridesに未知のカテゴリがあります: {sorted(unknown)}",
            unknown_keys=sorted(unknown),
        )
    out: dict = {}
    if "preprocess" in value:
        out["preprocess"] = _validate_routine_overrides_preprocess(
            value["preprocess"], "certificate.routine_overrides.preprocess")
    return out


def _validate_certificate_shape(certificate: object) -> dict:
    """証明書自体の構造（型・enum・必須キー）だけを検証する（内容の突合はしない）。"""
    if not isinstance(certificate, dict):
        _fail_certificate(f"certificateはオブジェクトである必要があります: {certificate!r}")
    try:
        json.dumps(certificate, allow_nan=False)
    except (TypeError, ValueError) as exc:
        _fail_certificate(f"certificateがJSON直列化できません（NaN/Infinity・非対応型は"
                          f"不可）: {exc}")

    unknown = set(certificate) - _CERTIFICATE_KEYS
    if unknown:
        _fail_certificate(f"certificateに未知のキーがあります: {sorted(unknown)}",
                          unknown_keys=sorted(unknown))
    # routine_overridesだけは任意（_CERTIFICATE_OPTIONAL_KEYS）——欠落を許すのは
    # このキーだけで、他の10キーはこれまでどおり必須のまま。
    missing = _CERTIFICATE_REQUIRED_KEYS - set(certificate)
    if missing:
        _fail_certificate(f"certificateに必須キーが不足しています: {sorted(missing)}",
                          missing_keys=sorted(missing))

    if certificate["schema"] != CERTIFICATE_SCHEMA:
        _fail_certificate(f"certificate.schemaは{CERTIFICATE_SCHEMA!r}のみ許可されます: "
                          f"{certificate['schema']!r}")

    if not isinstance(certificate["profile_content_sha256"], str) or not _SHA256_RE.fullmatch(
            certificate["profile_content_sha256"]):
        _fail_certificate("certificate.profile_content_sha256は64桁小文字16進のSHA-256"
                          f"文字列である必要があります: {certificate['profile_content_sha256']!r}")

    hash_maps = {}
    for cert_key, _obs_key in _HASH_MAP_CATEGORIES:
        hash_maps[cert_key] = _validate_hash_map(certificate[cert_key], f"certificate.{cert_key}")

    criteria_raw = certificate["criteria"]
    if not isinstance(criteria_raw, list) or not criteria_raw:
        _fail_certificate("certificate.criteriaは空でない配列である必要があります: "
                          f"{criteria_raw!r}")
    criteria = [_validate_criterion(item, i) for i, item in enumerate(criteria_raw)]
    seen_criterion_ids: set = set()
    for item in criteria:
        if item["criterion_id"] in seen_criterion_ids:
            _fail_certificate(f"certificate.criteria: criterion_idが重複しています: "
                              f"{item['criterion_id']!r}")
        seen_criterion_ids.add(item["criterion_id"])

    performed_by = certificate["performed_by"]
    if not isinstance(performed_by, str) or not performed_by:
        _fail_certificate(f"certificate.performed_byは空でない文字列である必要があります: "
                          f"{performed_by!r}")
    performed_at = certificate["performed_at"]
    if not isinstance(performed_at, str) or not performed_at:
        _fail_certificate(f"certificate.performed_atは空でない文字列である必要があります: "
                          f"{performed_at!r}")
    scope = certificate["scope"]
    if not isinstance(scope, str) or not scope:
        _fail_certificate(f"certificate.scopeは空でない文字列である必要があります: {scope!r}")

    routine_overrides = _validate_routine_overrides(certificate.get("routine_overrides"))

    return {
        "schema": certificate["schema"],
        "profile_content_sha256": certificate["profile_content_sha256"],
        **hash_maps,
        "criteria": criteria,
        "performed_by": performed_by,
        "performed_at": performed_at,
        "scope": scope,
        "routine_overrides": routine_overrides,
    }


def validate_certificate(profile: dict, certificate: dict, observed_hashes: dict) -> dict:
    """証明書がprofileを正当にvalidatedと名乗らせるものかを検証する（副作用なし）。

    `observed_hashes`は呼び出し側が実測した参照値で、次の4キーを持つ想定:
    `dependencies` / `fixed_inputs` / `reference_files` / `validation_outputs`
    （それぞれ `{id: sha256}` の辞書）。ここではハッシュの一致・適用範囲一致・
    必須基準のpassだけを検証し、署名や実測自体の真正性は検証しない
    （spec §5「署名認証基盤は作らない」）。

    不一致・不足はすべて`DomainError("PROFILE_VALIDATION_INVALID", ...)`。

    全検証に通った場合のみ、正規化済み証明書（`routine_overrides`を含む）を
    返す——呼び出し側（`metabolomix.pipeline.request_v2`）はこの戻り値の
    `routine_overrides`だけを使い、`execution_purpose="routine"`のpreprocess
    上書きをこの明示許容集合とだけ照合する（spec §6）。`routine_overrides`の
    値は、この関数が真正性（hash pin・profile内容・適用範囲・必須基準）を
    確認し終えた**後**でしか呼び出し側に渡らない——手書きで`validated`かつ
    それらしい`routine_overrides`を持つ証明書を用意しても、hash不一致等は
    このチェック列のどこよりも先に検出され、`routine_overrides`が「参照
    できる」段階にすら到達しない（下記の検証順序を参照）。
    """
    profile_validation = _require_dict(profile, "profile").get("validation")
    if not isinstance(profile_validation, dict):
        _fail_certificate("profile.validationがオブジェクトではありません（validate_profile"
                          "を通した正規化済みprofileを渡してください）。")

    if profile_validation.get("status") != "validated":
        _fail_certificate(
            "profile.validation.statusが'validated'ではありません（draft profileに"
            "証明書は適用できません）。", status=profile_validation.get("status"))

    pinned_hash = profile_validation.get("certificate_sha256")
    if not isinstance(pinned_hash, str) or not _SHA256_RE.fullmatch(pinned_hash):
        _fail_certificate("profile.validation.certificate_sha256が不正です。",
                          certificate_sha256=pinned_hash)

    normalized_certificate = _validate_certificate_shape(certificate)

    # 証明書のどの一要素を改変しても、pinされたhashと不一致になる
    # （証明書自体のcanonical hashをここで計算する。profile_content_hashが
    # profile本体の内容identityを守るのと対になる仕組み）。
    actual_hash = canonical_hash(certificate)
    if actual_hash != pinned_hash:
        _fail_certificate(
            "certificateのhashがprofile.validation.certificate_sha256と一致しません"
            "（証明書が改変されたか、pinされた証明書と異なります）。",
            expected=pinned_hash, actual=actual_hash)

    expected_profile_hash = profile_content_hash(profile)
    if normalized_certificate["profile_content_sha256"] != expected_profile_hash:
        _fail_certificate(
            "certificate.profile_content_sha256がprofile本体のcanonical hashと一致しません。",
            expected=expected_profile_hash,
            actual=normalized_certificate["profile_content_sha256"])

    if normalized_certificate["scope"] != profile_validation.get("scope"):
        _fail_certificate(
            "certificate.scopeがprofile.validation.scopeと一致しません（適用範囲の不一致）。",
            certificate_scope=normalized_certificate["scope"],
            profile_scope=profile_validation.get("scope"))

    observed = _require_dict(observed_hashes, "observed_hashes")
    for cert_key, obs_key in _HASH_MAP_CATEGORIES:
        expected = observed.get(obs_key, {})
        if not isinstance(expected, dict):
            expected = {}
        actual = normalized_certificate[cert_key]
        if actual != expected:
            _fail_certificate(
                f"certificate.{cert_key}が実測hash（observed_hashes[{obs_key!r}]）と"
                "一致しません。", category=cert_key, expected=expected, actual=actual)

    for item in normalized_certificate["criteria"]:
        if item["required"] and item["status"] != "pass":
            _fail_certificate(
                f"必須基準がpassではありません: criterion_id={item['criterion_id']!r} "
                f"status={item['status']!r}",
                criterion_id=item["criterion_id"], status=item["status"])

    return normalized_certificate
