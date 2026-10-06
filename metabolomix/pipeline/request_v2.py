"""pipeline-request.v2 の厳密な検証・解決・更新契約（spec §6, §6.2）。

v1（`metabolomix.pipeline.request`）が生データ由来の`method_file`/`lbm_file`を
直接受け取るのに対し、v2は検証済みLC–MSプロファイル（`lcms-profile.v1`、
`metabolomix.console.profile_schema`）を上流条件の唯一の出所とする。上流条件
（method・依存・極性・measure）はprofile経由でのみ解決し、v2ではこれらを
requestへ直接指定できない。

**このモジュールはファイルを読まない**（`metabolomix.console.profile_schema`と
同じ規約）——ただし1つだけ例外がある。`analysis-request.json`（v1と同じ
`REQUEST_FILE_NAME`）の実体読込・JSON解析・v2としての構造検証は、v1の
`read_request_file`と対になる`read_request_file`関数としてここに持つ
（v1がファイル読込を「request層自身」で行っているため、同じ層に置く——
`metabolomix/pipeline/request.py`のdispatch側docstring参照）。`profile_file`が
指す実体（profile本体）の読込・検証・hash計算はそれとは別で、呼び出し側の
責務のまま（後続タスクの`pipeline_plan`/`pipeline_run`側が担う）。

## 4段の優先順位

spec §6「値の優先順位はMCP明示値 > analysis-request.json > profile既定値 >
v2既定値」。`statistics`だけがこの4段すべてに乗る（他フィールドはprofile
既定値を持たないので実質3段: 明示値 > ファイル > v2既定値）。
`value_sources`の各キーには採用段を記録する——`statistics`は
`"explicit"` / `"request_file"` / `"profile_default"` / `"v2_default"`の
4値、他は`"explicit"` / `"request_file"` / `"default"`の3値。

## routine許容範囲: 証明書の`routine_overrides`とだけ照合する

spec §6「routineで証明書の許容範囲外になる場合は`PROFILE_SCOPE_MISMATCH`」を、
証明書（`lcms-profile-validation.v1`）が持つ`routine_overrides`という
明示的・構造化されたフィールドと照合することで実装する
（`profile_schema.validate_certificate`の戻り値。詳細は
`docs/schema/lcms-profile-v1.md`「証明書 `routine_overrides`」節）。
`execution_purpose="routine"`のpreprocess上書きは、この
`routine_overrides["preprocess"][recipe_id][field]`が
`{"mode": "any"}`（任意の値を許可）または`{"mode": "values", "values": [...]}`
（列挙した値だけ許可）と宣言していない限り、**すべて**`PROFILE_SCOPE_MISMATCH`
で拒否する（**fail-closed・explicit only**——`routine_overrides`省略時・
対象recipe_id/field未宣言時は無条件で拒否。profile本体の他の場所（他の
recipeの値等）から範囲を推測することは一切しない）。

このモジュールは以前の実装で「profile自身が持ついずれかのmatrix_recipes
エントリに既に現れる値」という代替規則を独自に発明していたが、controller
裁定によりこれは差し戻された——存在しない契約を肩代わりして実装すると、
後から見て「これが仕様だ」と誤読される。証明書側に構造化フィールドを
追加する形で解決済み（詳細はTask 3のfix report「Concern 2」参照）。

`execution_purpose="validation"`ではこの照合を一切行わない——validation実行は
まだ証明書が検証していない値を試すためのものだから。`feature_bindings`
選択の許容判定（§6.2「選択はprofileの許容規則内に限定し、外れれば
PROFILE_SCOPE_MISMATCH」）は実データ（dataset）へのアクセスを要するため
引き続き未実装——後続の`resolve_feature_bindings`stageの責務（下記参照）。

## feature_bindings（spec §6.2）の扱い

「statistics、feature_bindings、standard_assaysはresumeの更新対象とし、
元のrequestと同じ厳密検証を行う」——これはrequest **validation**の話であり、
実際に候補を解決する`resolve_feature_bindings`stage（後続task）とは別。
このモジュールが持つのは「更新payload（dataset hash・target_idごとの
feature_id・選択理由）」の**構造検証**だけ——候補一覧との突合・許容誤差・
証拠条件の判定は実データ（dataset）へのアクセスを要し、この層の責務外
（後続taskが担う）。初回`resolve`とresumeの`merge_updates`の両方が同じ
`_validate_feature_bindings`を通ることで、brief「feature_bindingsは初回指定も
resumeも同じ検査を通す」を満たす。
"""
from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path

from metabolomix.core.atomic_io import DomainError

__all__ = [
    "REQUEST_FILE_NAME",
    "SCHEMA",
    "UPDATABLE",
    "merge_updates",
    "read_request_file",
    "resolve",
    "validate_statistics",
]

SCHEMA = "pipeline-request.v2"

#: v1と同じファイル名・同じ場所（source_root直下）。v1の
#: `metabolomix.pipeline.request.REQUEST_FILE_NAME`と値は同一だが、循環import
#: を避けるためここでも定数として独立に持つ（`request.py`がこちらをimportする
#: 側であり、逆向きの依存は作らない）。
REQUEST_FILE_NAME = "analysis-request.json"

#: resumeで変更可能なトップレベルキー（spec §6.1「sample manifest・許可された
#: 前処理・統計定義・§8.1のバッチ固有対応付けの修正は新request revisionで
#: 再開する」、spec §6.2「statistics、feature_bindings、standard_assaysは
#: resumeの更新対象とし、元のrequestと同じ厳密検証を行う」）。
UPDATABLE = {
    "target", "sample_manifest", "preprocess", "statistics", "standard_assays",
    "feature_bindings",
}

_TOP_LEVEL_KEYS = frozenset({
    "schema", "omics", "profile_file", "execution_purpose", "target",
    "sample_manifest", "standard_assays", "preprocess", "statistics",
    "feature_bindings", "timeout_s", "save_project", "output_root", "keep_extension",
})

#: 内部専用（merge_updates内部でのみ現れる。外部入力には許可しない）。
_INTERNAL_KEYS = frozenset({"effective_target", "value_sources"})

_OMICS_VALUES = frozenset({"metabolomics"})
_EXECUTION_PURPOSE_VALUES = frozenset({"routine", "validation"})
_TARGET_VALUES = frozenset({"auto", "exploratory", "differential"})

_STATISTIC_KINDS = frozenset({"pca", "welch", "anova_tukey"})
_TRANSFORM_VALUES = frozenset({"none", "log2"})
_STATISTIC_COMMON_KEYS = frozenset({
    "statistic_id", "kind", "matrix_recipe_id", "transform", "feature_scope",
})
_STATISTIC_EXTRA_KEYS = {
    "pca": frozenset({"scaling", "n_components"}),
    "welch": frozenset({"reference_group", "test_group", "q_threshold", "log2fc_threshold"}),
    "anova_tukey": frozenset({"groups", "alpha"}),
}

_PREPROCESS_OVERRIDE_KEYS = frozenset({"normalize", "drift_correct", "filter", "impute"})
_NORMALIZE_VALUES = frozenset({"none", "tic", "median", "pqn"})
_IMPUTE_VALUES = frozenset({"none", "half_min", "knn", "column_mean"})
_FILTER_KEYS = frozenset({"min_detection_rate"})

_SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

_DEFAULT_TIMEOUT_S = 21600

#: statisticsが省略され、かつprofileにも既定統計が無い場合の唯一のfallback
#: （spec §6.2「既定PCAはstatistic_id=pca, matrix_recipe_id=default,
#: transform=none, scaling=autoscale, n_components=2」）。
_DEFAULT_PCA_STATISTIC = {
    "statistic_id": "pca", "kind": "pca", "matrix_recipe_id": "default",
    "transform": "none", "feature_scope": {"mode": "all_eligible"},
    "scaling": "autoscale", "n_components": 2,
}

#: feature_bindings更新payloadが受け付ける2キーちょうど（spec §6.2「更新payload
#: はdataset hash、target_idごとのfeature_id、選択理由」）。
_FEATURE_BINDINGS_KEYS = frozenset({"dataset_hash", "selections"})
_FEATURE_BINDING_SELECTION_KEYS = frozenset({"feature_id", "reason"})

#: 明示nullを拒否するトップレベルキー（v1の`_NULL_REJECTED_TOP_LEVEL_KEYS`と
#: 同じ規約——「無効化」を意味しない設定に明示nullを許すと未指定と区別が
#: つかなくなる）。`sample_manifest`だけがv1同様に例外（明示解除の語彙）。
#: `preprocess`はここに含めない——全体をnullにすることは「上書き無し」と
#: 同義であり、無効化ではなく単なる省略の言い換えとして扱う。
_NULL_REJECTED_TOP_LEVEL_KEYS = frozenset({
    "schema", "omics", "profile_file", "execution_purpose", "target",
    "standard_assays", "statistics", "feature_bindings", "output_root", "keep_extension",
})


def _fail(message: str, **details) -> None:
    raise DomainError("PIPELINE_REQUEST_INVALID", message, details)


def _fail_scope(message: str, **details) -> None:
    raise DomainError("PROFILE_SCOPE_MISMATCH", message, details)


def _reject_disallowed_explicit_null(source: dict, keys) -> None:
    for key in keys:
        if key in source and source[key] is None:
            _fail(f"{key}に明示nullは許可されません（無効化を意味しない設定のため、"
                  f"未指定にしてください）。", **{key: None})


def _is_finite_number(value: object) -> bool:
    if isinstance(value, bool):
        return False
    if not isinstance(value, (int, float)):
        return False
    return math.isfinite(float(value))


def _require_positive_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        _fail(f"{label}は正の整数である必要があります（boolは不可）: {value!r}",
              **{label: value})
    return value


def _require_probability(value: object, label: str) -> float:
    if not (_is_finite_number(value) and 0 < float(value) < 1):
        _fail(f"{label}は0より大きく1未満である必要があります: {value!r}",
              **{label: value})
    return float(value)


def _require_optional_str(data: dict, field: str) -> None:
    value = data.get(field)
    if value is not None and (not isinstance(value, str) or value == ""):
        _fail(f"{field}は空でない文字列またはnullである必要があります: {value!r}",
              **{field: value})


# ---------- statistics（spec §6.2） ----------

def _validate_feature_scope(value: object, label: str, feature_target_ids: frozenset) -> dict:
    if not isinstance(value, dict):
        _fail(f"{label}はオブジェクトである必要があります: {value!r}")
    mode = value.get("mode")
    if mode == "all_eligible":
        unknown = set(value) - {"mode"}
        if unknown:
            _fail(f"{label}に未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))
        return {"mode": "all_eligible"}
    if mode == "targets":
        unknown = set(value) - {"mode", "target_ids"}
        if unknown:
            _fail(f"{label}に未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))
        target_ids = value.get("target_ids")
        if not isinstance(target_ids, list) or not target_ids:
            _fail(f"{label}.target_idsは空でない配列である必要があります: {target_ids!r}")
        if len(set(target_ids)) != len(target_ids):
            _fail(f"{label}.target_idsに重複があります。", target_ids=target_ids)
        for target_id in target_ids:
            if target_id not in feature_target_ids:
                _fail(f"{label}.target_idsが存在しないtarget_idを参照しています: "
                      f"{target_id!r}", target_id=target_id)
        return {"mode": "targets", "target_ids": list(target_ids)}
    _fail(f"{label}.modeが不正です: {mode!r}", mode=mode)


def _validate_statistic_item(
    value: object, index: int, matrix_recipe_ids: frozenset, feature_target_ids: frozenset,
) -> dict:
    label = f"statistics[{index}]"
    if not isinstance(value, dict):
        _fail(f"{label}はオブジェクトである必要があります: {value!r}")

    missing_common = _STATISTIC_COMMON_KEYS - set(value)
    if missing_common:
        _fail(f"{label}に必須キーが不足しています: {sorted(missing_common)}",
              missing_keys=sorted(missing_common))

    kind = value["kind"]
    if kind not in _STATISTIC_KINDS:
        _fail(f"{label}.kindが不正です: {kind!r}", kind=kind)

    allowed = _STATISTIC_COMMON_KEYS | _STATISTIC_EXTRA_KEYS[kind]
    unknown = set(value) - allowed
    if unknown:
        _fail(f"{label}に未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))

    statistic_id = value["statistic_id"]
    if not isinstance(statistic_id, str) or not _SAFE_ID_RE.fullmatch(statistic_id):
        _fail(f"{label}.statistic_idは安全な文字だけのID（英数字・アンダースコア・"
              f"ハイフン）である必要があります: {statistic_id!r}", statistic_id=statistic_id)

    matrix_recipe_id = value["matrix_recipe_id"]
    if matrix_recipe_id not in matrix_recipe_ids:
        _fail(f"{label}.matrix_recipe_idが存在しないrecipeを参照しています: "
              f"{matrix_recipe_id!r}", matrix_recipe_id=matrix_recipe_id)

    transform = value["transform"]
    if transform not in _TRANSFORM_VALUES:
        _fail(f"{label}.transformが不正です: {transform!r}", transform=transform)

    feature_scope = _validate_feature_scope(
        value["feature_scope"], f"{label}.feature_scope", feature_target_ids)

    out = {
        "statistic_id": statistic_id, "kind": kind, "matrix_recipe_id": matrix_recipe_id,
        "transform": transform, "feature_scope": feature_scope,
    }

    if kind == "pca":
        scaling = value.get("scaling", "autoscale")
        if scaling not in ("none", "autoscale"):
            _fail(f"{label}.scalingが不正です: {scaling!r}", scaling=scaling)
        n_components = value.get("n_components", 2)
        out["scaling"] = scaling
        out["n_components"] = _require_positive_int(n_components, f"{label}.n_components")
    elif kind == "welch":
        reference_group = value["reference_group"]
        test_group = value["test_group"]
        if not isinstance(reference_group, str) or not reference_group:
            _fail(f"{label}.reference_groupが不正です: {reference_group!r}")
        if not isinstance(test_group, str) or not test_group:
            _fail(f"{label}.test_groupが不正です: {test_group!r}")
        if reference_group == test_group:
            _fail(f"{label}: reference_groupとtest_groupが同一です: {reference_group!r}",
                  group=reference_group)
        q_threshold = value.get("q_threshold", 0.05)
        log2fc_threshold = value.get("log2fc_threshold", 1.0)
        out["reference_group"] = reference_group
        out["test_group"] = test_group
        out["q_threshold"] = _require_probability(q_threshold, f"{label}.q_threshold")
        if not (_is_finite_number(log2fc_threshold) and float(log2fc_threshold) >= 0):
            _fail(f"{label}.log2fc_thresholdは非負である必要があります: "
                  f"{log2fc_threshold!r}", log2fc_threshold=log2fc_threshold)
        out["log2fc_threshold"] = float(log2fc_threshold)
    else:  # anova_tukey
        groups = value["groups"]
        if not isinstance(groups, list) or len(groups) < 3:
            _fail(f"{label}.groupsは3群以上である必要があります: {groups!r}")
        if len(set(groups)) != len(groups):
            _fail(f"{label}.groupsに重複があります: {groups!r}", groups=groups)
        for group in groups:
            if not isinstance(group, str) or not group:
                _fail(f"{label}.groups[]は空でない文字列である必要があります: {group!r}")
        alpha = value.get("alpha", 0.05)
        out["groups"] = list(groups)
        out["alpha"] = _require_probability(alpha, f"{label}.alpha")

    return out


def validate_statistics(items: list, profile: dict) -> list:
    """`statistics`配列を検証し、既定値を埋めた正規化済みリストを返す（spec §6.2）。

    ``profile``は``matrix_recipes``（recipe_id存在確認）と``feature_targets``
    （target_id存在確認）だけを参照する——``lcms-profile.v1``の完全な形で
    ある必要はない（呼び出し側のテストは最小限のdictでよい）。
    """
    if not isinstance(profile, dict):
        _fail("profileはオブジェクトである必要があります。", value=profile)
    matrix_recipe_ids = frozenset((profile.get("matrix_recipes") or {}).keys())
    feature_target_ids = frozenset((profile.get("feature_targets") or {}).keys())

    if not isinstance(items, list):
        _fail("statisticsは配列である必要があります。", value=items)
    if not items:
        _fail("statisticsは非空の配列である必要があります。")

    normalized: list = []
    seen_ids: set = set()
    for index, item in enumerate(items):
        normalized_item = _validate_statistic_item(
            item, index, matrix_recipe_ids, feature_target_ids)
        statistic_id = normalized_item["statistic_id"]
        if statistic_id in seen_ids:
            _fail(f"statistic_idが重複しています: {statistic_id!r}", statistic_id=statistic_id)
        seen_ids.add(statistic_id)
        normalized.append(normalized_item)
    return normalized


def _effective_target(statistics: list) -> str:
    """spec §6.2「target=autoは検定を含めばdifferential、PCAだけなら
    exploratory」。"""
    return "differential" if any(s["kind"] != "pca" for s in statistics) else "exploratory"


def _compute_effective_target(target: str, statistics: list) -> str:
    if target != "auto":
        return target
    return _effective_target(statistics)


def _validate_target_consistency(target: str, statistics: list) -> None:
    """spec §6.2「exploratoryに検定を指定、またはdifferentialに検定がなければ
    request不正とする」。"""
    has_test = any(s["kind"] != "pca" for s in statistics)
    if target == "exploratory" and has_test:
        _fail(
            "target='exploratory'に検定(welch/anova_tukey)を含むstatisticsは指定"
            "できません（検定を実行するならtarget='differential'または'auto'に"
            "してください）。",
            target=target,
            statistic_ids=[s["statistic_id"] for s in statistics if s["kind"] != "pca"],
        )
    if target == "differential" and not has_test:
        _fail(
            "target='differential'には検定(welch/anova_tukey)を含むstatisticsが"
            "少なくとも1つ必要です。",
            target=target,
        )


# ---------- standard_assays ----------

def _validate_standard_assays(value: object, feature_target_ids: frozenset) -> dict:
    if value is None:
        value = {}
    if not isinstance(value, dict):
        _fail("standard_assaysはオブジェクトである必要があります。", value=value)
    out: dict = {}
    for target_id, sample_ids in value.items():
        if target_id not in feature_target_ids:
            _fail(f"standard_assaysが存在しないtarget_idを参照しています: "
                  f"{target_id!r}", target_id=target_id)
        if not isinstance(sample_ids, list) or not sample_ids:
            _fail(f"standard_assays[{target_id!r}]は空でない配列である必要があります: "
                  f"{sample_ids!r}", target_id=target_id)
        normalized_ids = []
        for sample_id in sample_ids:
            if not isinstance(sample_id, str) or not sample_id:
                _fail(f"standard_assays[{target_id!r}]の要素は空でない文字列である"
                      f"必要があります: {sample_id!r}", target_id=target_id)
            normalized_ids.append(sample_id)
        out[target_id] = normalized_ids
    return out


# ---------- feature_bindings（spec §6.2 resume payload、構造検証のみ） ----------

def _validate_feature_bindings(value: object, feature_target_ids: frozenset) -> dict | None:
    """resume/初回共通の``feature_bindings``payloadを構造検証する。

    spec §6.2「更新payloadはdataset hash、target_idごとのfeature_id、選択
    理由」。ここで確認するのはこの3つの**形**だけ——選択されたfeature_idが
    実際にそのdataset上で許容誤差・証拠条件を満たすかどうかの判定
    （§6.2「選択はprofileの許容規則内に限定し、外れればPROFILE_SCOPE_MISMATCH」）
    は実データ（当該datasetの`assay-feature-evidence.v1`等）へのアクセスを
    要する``resolve_feature_bindings``stage（後続task）の責務であり、ここでは
    実装しない（このモジュールはprofileだけを受け取り、datasetを受け取らない）。
    """
    if value is None:
        return None
    if not isinstance(value, dict):
        _fail("feature_bindingsはオブジェクトである必要があります。", value=value)
    unknown = set(value) - _FEATURE_BINDINGS_KEYS
    if unknown:
        _fail(f"feature_bindingsに未知のキーがあります: {sorted(unknown)}",
              unknown_keys=sorted(unknown))
    missing = _FEATURE_BINDINGS_KEYS - set(value)
    if missing:
        _fail(f"feature_bindingsに必須キーが不足しています: {sorted(missing)}",
              missing_keys=sorted(missing))

    dataset_hash = value["dataset_hash"]
    if not isinstance(dataset_hash, str) or not _SHA256_RE.fullmatch(dataset_hash):
        _fail(f"feature_bindings.dataset_hashは64桁小文字16進のSHA-256文字列で"
              f"ある必要があります: {dataset_hash!r}", dataset_hash=dataset_hash)

    selections = value["selections"]
    if not isinstance(selections, dict) or not selections:
        _fail(f"feature_bindings.selectionsは空でないオブジェクトである必要が"
              f"あります: {selections!r}", selections=selections)

    normalized_selections: dict = {}
    for target_id, selection in selections.items():
        if target_id not in feature_target_ids:
            _fail(f"feature_bindings.selectionsが存在しないtarget_idを参照して"
                  f"います: {target_id!r}", target_id=target_id)
        label = f"feature_bindings.selections[{target_id!r}]"
        if not isinstance(selection, dict):
            _fail(f"{label}はオブジェクトである必要があります: {selection!r}")
        unknown_sel = set(selection) - _FEATURE_BINDING_SELECTION_KEYS
        if unknown_sel:
            _fail(f"{label}に未知のキーがあります: {sorted(unknown_sel)}",
                  unknown_keys=sorted(unknown_sel))
        missing_sel = _FEATURE_BINDING_SELECTION_KEYS - set(selection)
        if missing_sel:
            _fail(f"{label}に必須キーが不足しています: {sorted(missing_sel)}",
                  missing_keys=sorted(missing_sel))

        feature_id = selection["feature_id"]
        if not isinstance(feature_id, str) or not feature_id:
            _fail(f"{label}.feature_idは空でない文字列である必要があります: "
                  f"{feature_id!r}")
        reason = selection["reason"]
        if not isinstance(reason, str) or not reason:
            _fail(f"{label}.reasonは空でない文字列である必要があります: {reason!r}")
        normalized_selections[target_id] = {"feature_id": feature_id, "reason": reason}

    return {"dataset_hash": dataset_hash, "selections": normalized_selections}


# ---------- preprocess（matrix_recipe override） ----------

def _validate_matrix_recipe_override(value: object, label: str, recipe: dict) -> dict:
    if not isinstance(value, dict):
        _fail(f"{label}はオブジェクトである必要があります: {value!r}")
    unknown = set(value) - _PREPROCESS_OVERRIDE_KEYS
    if unknown:
        _fail(f"{label}に未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))

    out: dict = {}
    if "normalize" in value:
        normalize = value["normalize"]
        if normalize not in _NORMALIZE_VALUES:
            _fail(f"{label}.normalizeが不正です: {normalize!r}", normalize=normalize)
        out["normalize"] = normalize
    if "drift_correct" in value:
        drift_correct = value["drift_correct"]
        if not isinstance(drift_correct, bool):
            _fail(f"{label}.drift_correctはboolである必要があります: {drift_correct!r}")
        out["drift_correct"] = drift_correct
    if "impute" in value:
        impute = value["impute"]
        if impute not in _IMPUTE_VALUES:
            _fail(f"{label}.imputeが不正です: {impute!r}", impute=impute)
        out["impute"] = impute
    if "filter" in value:
        filter_value = value["filter"]
        if filter_value is not None:
            if not isinstance(filter_value, dict) or set(filter_value) != _FILTER_KEYS:
                _fail(f"{label}.filterはnullまたは{{min_detection_rate}}である"
                      f"必要があります: {filter_value!r}")
            rate = filter_value["min_detection_rate"]
            if not (_is_finite_number(rate) and 0.0 <= float(rate) <= 1.0):
                _fail(f"{label}.filter.min_detection_rateは0以上1以下である"
                      f"必要があります: {rate!r}")
            filter_value = {"min_detection_rate": float(rate)}
        out["filter"] = filter_value

    # base='internal_standard_ratio'のrecipeに対するnormalize上書きは、
    # 上書き後も二重正規化禁止（spec §8）を満たす必要がある
    # （profile_schema._validate_matrix_recipeと同じ規則をoverride後の
    # 実効値に対しても適用する）。この規則はexecution_purposeに関わらず
    # 常に適用する——routine/validationの区別が要る話ではなく、profile自身の
    # 契約（base×normalizeの組合せ）を破らないという構造的な制約だから。
    effective_normalize = out.get("normalize", recipe.get("normalize"))
    effective_base = recipe.get("base")
    if effective_base == "internal_standard_ratio" and effective_normalize != "none":
        _fail(
            f"{label}: base='internal_standard_ratio'のrecipeにnormalize="
            f"{effective_normalize!r}への上書きは許可されません（二重正規化の禁止）。",
            base=effective_base, normalize=effective_normalize,
        )
    return out


def _check_routine_scope(
    recipe_id: str, field: str, value: object, routine_overrides: dict | None,
) -> None:
    """`execution_purpose='routine'`のpreprocess上書きを、証明書の明示許容集合
    （`profile_schema.validate_certificate`の戻り値の`routine_overrides`）と
    だけ照合する（spec §6・controller裁定）。

    **fail-closed・explicit only**: ``routine_overrides``が``None``・
    対象のcategoryが無い・対象の``recipe_id``が無い・対象の``field``が
    宣言されていない、のいずれでも拒否する。profile本体（`matrix_recipes`
    の他エントリ等）から範囲を推測することは一切しない——以前の実装
    （「profile自身が使っている値なら許可」）はこの理由でcontrollerに
    差し戻された。
    """
    preprocess_allowances = (routine_overrides or {}).get("preprocess") or {}
    recipe_allowances = preprocess_allowances.get(recipe_id) or {}
    allowance = recipe_allowances.get(field)
    if not isinstance(allowance, dict):
        _fail_scope(
            f"execution_purpose='routine'ではpreprocess[{recipe_id!r}].{field}の"
            "上書きは証明書の明示許容集合(routine_overrides)に含まれていません。",
            recipe_id=recipe_id, field=field,
        )
        return
    mode = allowance.get("mode")
    if mode == "any":
        return
    if mode == "values":
        allowed_values = allowance.get("values") or []
        if value not in allowed_values:
            _fail_scope(
                f"execution_purpose='routine'ではpreprocess[{recipe_id!r}].{field}を"
                f"{value!r}へ上書きできません（証明書が許可する値: {allowed_values!r}）。",
                recipe_id=recipe_id, field=field, value=value, allowed_values=allowed_values,
            )
        return
    # routine_overridesはprofile_schema.validate_certificateが検証済みの
    # 前提で受け取るため、通常ここへは来ない（防御的な扱い）。
    _fail_scope(
        f"証明書のroutine_overrides[{recipe_id!r}][{field!r}]が不正です（mode不明）: "
        f"{mode!r}", recipe_id=recipe_id, field=field,
    )


def _validate_preprocess(
    value: object, matrix_recipes: dict, execution_purpose: str, routine_overrides: dict | None,
) -> dict:
    """spec §6「preprocess: profileの既定値を明示指定で上書き可能」。

    ``execution_purpose="routine"``では、上書きした各フィールドを
    ``routine_overrides``（証明書由来の明示許容集合）と照合する
    （`_check_routine_scope`）。``"validation"``ではこの照合をしない——
    validation実行はまだ証明書が検証していない値を試すためのものだから。
    """
    if value is None:
        return {}
    if not isinstance(value, dict):
        _fail("preprocessはオブジェクトである必要があります。", value=value)

    out: dict = {}
    for recipe_id, override in value.items():
        if recipe_id not in matrix_recipes:
            _fail(f"preprocessが存在しないrecipeを参照しています: {recipe_id!r}",
                  recipe_id=recipe_id)
        normalized_override = _validate_matrix_recipe_override(
            override, f"preprocess[{recipe_id!r}]", matrix_recipes[recipe_id])
        if execution_purpose == "routine":
            for field, field_value in normalized_override.items():
                _check_routine_scope(recipe_id, field, field_value, routine_overrides)
        out[recipe_id] = normalized_override
    return out


# ---------- analysis-request.json（v1と同じ場所・同じファイル名） ----------

def read_request_file(source_root: Path) -> dict:
    """元フォルダ直下の``analysis-request.json``をv2 shapeとして読む（spec §6）。

    無ければ空dict。schemaキーの整合性はここでは検証しない（呼び出し側の
    ``metabolomix.pipeline.request.resolve_request``が、explicit/fileどちらの
    schemaを採用するかを既に決定した上でこの関数を呼ぶ——「v2として読む」と
    決まった後の構造検証だけがここの責務）。不正JSON・非オブジェクト・v2に
    存在しないキー・許可されない明示nullは、v1の``read_request_file``と同じ
    考え方で即座に``PIPELINE_REQUEST_INVALID``にする——「置いてあるのに
    黙って無視された」を避ける。
    """
    path = Path(source_root) / REQUEST_FILE_NAME
    if not path.is_file():
        return {}
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError as exc:
        _fail(f"{REQUEST_FILE_NAME}を読めません: {exc}", path=str(path))
    try:
        data = json.loads(raw)
    except ValueError as exc:
        _fail(f"{REQUEST_FILE_NAME}がJSONとして壊れています: {exc}", path=str(path))
    if not isinstance(data, dict):
        _fail(f"{REQUEST_FILE_NAME}はオブジェクトである必要があります。",
              path=str(path), value=data)
    _check_known_keys(data)
    _reject_disallowed_explicit_null(data, _NULL_REJECTED_TOP_LEVEL_KEYS)
    return data


def _check_known_keys(source: dict) -> None:
    """``_TOP_LEVEL_KEYS``にないキーを、v1からの移行者向けの案内つきで拒否する。"""
    unknown = set(source) - _TOP_LEVEL_KEYS
    if not unknown:
        return
    if {"method_file", "lbm_file"} & unknown:
        _fail(
            "pipeline-request.v2ではmethod_file/lbm_fileを直接指定できません"
            "（profile_fileから解決してください）。", unknown_keys=sorted(unknown),
        )
    if "comparisons" in unknown:
        _fail(
            "pipeline-request.v2ではcomparisonsではなくstatisticsを指定して"
            "ください。", unknown_keys=sorted(unknown),
        )
    _fail(f"未知のキーがあります: {sorted(unknown)}", unknown_keys=sorted(unknown))


def _layer_top_level(explicit: dict, from_file: dict) -> tuple[dict, dict]:
    """トップレベル各キーを「explicit > from_file > (欠落のまま)」で重ね、
    (merged, sources) を返す。``sources``はここで決着済み——``merged``に
    キーが無ければ、``_resolve_core``が既定値で埋める（そのケースの出所は
    ``"default"``、または``statistics``だけ``_resolve_core``がprofile既定値/
    v2既定値のどちらかへさらに解決する）。
    """
    merged: dict = {}
    sources: dict = {}
    for key in _TOP_LEVEL_KEYS:
        if key in explicit:
            merged[key] = explicit[key]
            sources[key] = "explicit"
        elif key in from_file:
            merged[key] = from_file[key]
            sources[key] = "request_file"
        else:
            sources[key] = "default"
    return merged, sources


# ---------- 公開API ----------

def resolve(
    data: dict, profile: dict, *, from_file: dict | None = None,
    routine_overrides: dict | None = None,
) -> dict:
    """`pipeline-request.v2`を検証・解決する（spec §6, §6.2）。

    ``data``は「MCP明示値」の層そのもの。``from_file``（省略可）は
    ``source_root``直下の``analysis-request.json``の内容（``read_request_file``
    で既にv2 shapeとして検証済み）——渡すと優先順位は「``data`` >
    ``from_file`` > profile既定値（statisticsのみ） > v2既定値」の4段
    （spec §6）になる。省略時（``None``、既定）は「``data`` > profile既定値 >
    v2既定値」の3段——既存の呼び出し・テストはこちらの経路のまま。

    ``profile``は既に検証済みの``lcms-profile.v1``相当のdict（このモジュール
    が実際に参照するのは``matrix_recipes`` / ``feature_targets`` /
    ``analysis_recipe.statistics``の3つだけ）。``routine_overrides``（省略可）
    は``profile_schema.validate_certificate``の戻り値の同名キー——
    ``execution_purpose="routine"``のpreprocess上書きをこの明示許容集合とだけ
    照合する（spec §6）。**省略時（``None``）はfail-closed**: routineでは
    何も上書きできない。

    未指定と明示nullの区別（``_NULL_REJECTED_TOP_LEVEL_KEYS``）は**新規の
    明示入力に対してだけ**行う——``data``と``from_file``それぞれに対して
    個別に検査する。``merge_updates``が``_resolve_core``を直接呼ぶのは
    このため：既に解決済みの要求（``keep_extension``等の未指定フィールドが
    素の``None``で埋まっている）を丸ごとここへ再度通すと、「省略していた」
    という事実が失われ、すべてのNoneが「明示null」に見えてしまう（v1の
    ``resolve_request`` vs ``validate_request``と同じ区別）。
    """
    if not isinstance(data, dict):
        _fail("requestはオブジェクトである必要があります。", value=data)
    _check_known_keys(data)
    _reject_disallowed_explicit_null(data, _NULL_REJECTED_TOP_LEVEL_KEYS)

    if from_file is None:
        from_file = {}
    elif not isinstance(from_file, dict):
        _fail("analysis-request.jsonの内容はオブジェクトである必要があります。",
              value=from_file)

    merged, sources = _layer_top_level(data, from_file)
    resolved, statistics_fallback_source = _resolve_core(merged, profile, routine_overrides)
    resolved["value_sources"] = sources
    if "statistics" not in merged:
        resolved["value_sources"]["statistics"] = statistics_fallback_source
    return resolved


def _resolve_core(
    data: dict, profile: dict, routine_overrides: dict | None = None,
) -> tuple[dict, str | None]:
    """``resolve``/``merge_updates``が共有する検証本体（明示null拒否を含まない）。

    戻り値は``(解決済みdict（value_sourcesを含まない）, statisticsが
    dataに無かった場合のfallback出所（"profile_default"/"v2_default"）
    ——dataに"statistics"があれば``None``)``の2要素タプル。呼び出し側
    （``resolve``/``merge_updates``）がそれぞれの流儀で``value_sources``を
    組み立てる。
    """
    if not isinstance(profile, dict):
        _fail("profileはオブジェクトである必要があります"
              "（v2はprofileの解決結果が必須です）。", value=profile)
    _check_known_keys(data)

    matrix_recipes = profile.get("matrix_recipes") or {}
    feature_target_ids = frozenset((profile.get("feature_targets") or {}).keys())

    schema = data.get("schema", SCHEMA)
    if schema != SCHEMA:
        _fail(f"schemaは{SCHEMA!r}のみ許可されます: {schema!r}", schema=schema)

    omics = data.get("omics", "metabolomics")
    if omics not in _OMICS_VALUES:
        _fail(f"omicsが不正です: {omics!r}", omics=omics)

    profile_file = data.get("profile_file")
    if not isinstance(profile_file, str) or not profile_file:
        _fail(f"profile_fileは空でない文字列である必要があります: {profile_file!r}",
              profile_file=profile_file)

    execution_purpose = data.get("execution_purpose", "routine")
    if execution_purpose not in _EXECUTION_PURPOSE_VALUES:
        _fail(f"execution_purposeが不正です: {execution_purpose!r}",
              execution_purpose=execution_purpose)

    target = data.get("target", "auto")
    if target not in _TARGET_VALUES:
        _fail(f"targetが不正です: {target!r}", target=target)

    _require_optional_str(data, "sample_manifest")
    sample_manifest = data.get("sample_manifest")

    _require_optional_str(data, "output_root")
    output_root = data.get("output_root")

    _require_optional_str(data, "keep_extension")
    keep_extension = data.get("keep_extension")

    timeout_s = data.get("timeout_s", _DEFAULT_TIMEOUT_S)
    if isinstance(timeout_s, bool) or not isinstance(timeout_s, int) or timeout_s <= 0:
        _fail(f"timeout_sは正の整数である必要があります（boolは不可）: {timeout_s!r}",
              timeout_s=timeout_s)

    save_project = data.get("save_project", True)
    if not isinstance(save_project, bool):
        _fail(f"save_projectはboolのみ許可されます: {save_project!r}",
              save_project=save_project)

    statistics_fallback_source = None
    if "statistics" in data:
        statistics_input = data["statistics"]
    else:
        profile_statistics = (profile.get("analysis_recipe") or {}).get("statistics") or []
        if profile_statistics:
            statistics_fallback_source = "profile_default"
            statistics_input = profile_statistics
        else:
            statistics_fallback_source = "v2_default"
            statistics_input = [copy.deepcopy(_DEFAULT_PCA_STATISTIC)]
    statistics = validate_statistics(statistics_input, profile)
    _validate_target_consistency(target, statistics)
    effective_target = _compute_effective_target(target, statistics)

    standard_assays = _validate_standard_assays(data.get("standard_assays"), feature_target_ids)
    feature_bindings = _validate_feature_bindings(
        data.get("feature_bindings"), feature_target_ids)
    preprocess = _validate_preprocess(
        data.get("preprocess"), matrix_recipes, execution_purpose, routine_overrides)

    resolved = {
        "schema": SCHEMA,
        "omics": omics,
        "profile_file": profile_file,
        "execution_purpose": execution_purpose,
        "target": target,
        "sample_manifest": sample_manifest,
        "standard_assays": standard_assays,
        "feature_bindings": feature_bindings,
        "preprocess": preprocess,
        "statistics": statistics,
        "timeout_s": timeout_s,
        "save_project": save_project,
        "output_root": output_root,
        "keep_extension": keep_extension,
        "effective_target": effective_target,
    }
    return resolved, statistics_fallback_source


def merge_updates(
    current: dict, updates: dict, profile: dict, *, routine_overrides: dict | None = None,
) -> dict:
    """resumeの入力訂正を反映し、再検証した v2 request を返す（spec §6.1, §6.2）。

    ``UPDATABLE``（target・sample_manifest・preprocess・statistics・
    standard_assays・feature_bindings）以外のキーを変えようとした場合は
    ``NEW_PIPELINE_REQUIRED``にする——profile_file・execution_purpose・omics等
    の変更は上流条件そのものの変更であり、新しいpipelineが要る（spec §6.1
    「profile・method・library・raw・実行環境・極性・measureの変更は
    NEW_PIPELINE_REQUIRED」）。

    ``feature_bindings``の更新は初回``resolve``と全く同じ``_validate_feature_bindings``
    を通る（``_resolve_core``経由）——brief「feature_bindingsは初回指定も
    resumeも同じ検査を通す」はこの共有によって満たされる。``preprocess``の
    更新も同様に、初回``resolve``と全く同じ``_check_routine_scope``を通る
    ——``routine_overrides``は呼び出し側が都度渡す（resumeでも証明書の
    再検証・再取得は呼び出し側の責務のまま）。
    """
    if not isinstance(current, dict):
        _fail("requestはオブジェクトである必要があります。", value=current)
    if not isinstance(updates, dict):
        _fail("updatesはオブジェクトである必要があります。", value=updates)

    unknown = set(updates) - UPDATABLE
    if unknown:
        raise DomainError(
            "NEW_PIPELINE_REQUIRED", "上流条件の変更には新しい解析が必要です",
            {"unknown_keys": sorted(unknown)},
        )
    _reject_disallowed_explicit_null(updates, _NULL_REJECTED_TOP_LEVEL_KEYS & UPDATABLE)

    merged = {key: copy.deepcopy(value) for key, value in current.items()
              if key not in _INTERNAL_KEYS}
    for key, value in updates.items():
        merged[key] = copy.deepcopy(value)

    resolved, _fallback_source = _resolve_core(merged, profile, routine_overrides)

    # 出所は「保存済みの既存値」を土台に、updatesへ挙がったキーだけを
    # explicit_updateへ格上げする——merged自体は常に全キー揃っているため、
    # _resolve_coreの出力だけからは「更新していないフィールドの元の出所」が
    # 再現できない。
    resolved["value_sources"] = copy.deepcopy(current.get("value_sources", {}))
    for key in updates:
        resolved["value_sources"][key] = "explicit_update"
    return resolved
