"""v2 統計: 変換・効果量・検定（spec §10）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

v1 の `differential.two_group_test` は `log2(x + pseudo_count)` で変換し、効果量も
log 空間の平均差から作る。それは v1 の契約として**そのまま変えない**（過去の
結果と比較できなくなる）。v2 は別の定義を持ち、両者は共存する:

**変換。** `log2` は有限の正値にだけ適用する。0以下を1へ clip したり pseudocount を
足して「変換できた」ことにしない——clip は 0 と 1 を同じ値に潰し、pseudocount は
比の大きさを「足した定数の大きさ」にすり替える。変換できない値は欠損にする。

**効果量。** `effect_size_definition = log2_arithmetic_mean_ratio`。統計変換**前**の
test/reference の**算術**平均比の log2 で、log 空間の平均差（＝幾何平均比）とは
別物。どちらも「log2FC」と呼ばれるので、定義を結果に書いて持ち回る。

**BH の母集団。** 同じ statistic_id で**検定できた** feature 全体。検定不能を母集団へ
入れると q が不当に緩み、しかも「何件入れたか」を書かないと再現できない。

**母集団の検査は検定の前。** 群が n<2 なら止める。n 不足の群を落として別の検定に
すり替えない（spec §10）——落とした結果は、指定した比較とは違う比較の答えになる。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis import differential, multigroup
from metabolomix.analysis.sample_manifest import (
    select_statistical_samples,
    validate_independent_samples,
)
from metabolomix.core.atomic_io import DomainError

__all__ = [
    "EFFECT_SIZE_DEFINITION",
    "STATISTIC_SCHEMA",
    "arithmetic_log2fc",
    "run_statistic",
    "transform_values",
]

STATISTIC_SCHEMA = "statistic-result.v1"

#: 結果に必ず載せる効果量の定義。v1 の log 空間の平均差と取り違えないための印。
EFFECT_SIZE_DEFINITION = "log2_arithmetic_mean_ratio"

_SPEC_INVALID = "STATISTIC_SPECIFICATION_INVALID"

_TRANSFORMS = ("none", "log2")
_KINDS = ("pca", "welch", "anova_tukey")

#: 各群に必要な独立試料数。1点の群は群内分散を持たない。
_MIN_GROUP_N = 2


def transform_values(values: np.ndarray, transform: str) -> np.ndarray:
    """統計変換（元の配列は変更しない）。

    `log2` は有限の正値だけに適用し、0以下・非有限は NaN にする。clip も
    pseudocount も入れない（spec §10「log2は有限正値のみ許可し、0以下は欠損化
    して理由を保存する。1へのclipや暗黙pseudocountは使わない」）。
    """
    if transform not in _TRANSFORMS:
        raise DomainError(
            _SPEC_INVALID,
            f"未知の変換です: {transform!r}（対応: {', '.join(_TRANSFORMS)}）",
            {"transform": transform, "supported": list(_TRANSFORMS)})
    array = np.array(values, dtype=float, copy=True)
    if transform == "none":
        return array
    out = np.full(array.shape, np.nan, dtype=float)
    usable = np.isfinite(array) & (array > 0)
    out[usable] = np.log2(array[usable])
    return out


def arithmetic_log2fc(a, b) -> float:
    """`log2(mean(b) / mean(a))`。a=reference、b=test。

    平均は**算術**平均で、統計変換前の値に対して取る。どちらかの平均が 0 以下
    なら NA——微小定数を足して「計算できた」ことにしない。
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if not len(a) or not len(b) or a.mean() <= 0 or b.mean() <= 0:
        return float("nan")
    return float(np.log2(b.mean() / a.mean()))


# ---------- 母集団 ----------

def _row_index(metadata: list[dict], matrix: dict) -> dict:
    """sample_id → 解析行列の行番号（metadata は assay 列と同じ並び）。"""
    assay_ids = list(matrix["assay_ids"])
    if len(metadata) != len(assay_ids):
        raise DomainError(
            _SPEC_INVALID,
            f"metadataの件数が解析行列のassay数と一致しません"
            f"（metadata={len(metadata)}, assays={len(assay_ids)}）。",
            {"metadata": len(metadata), "assays": len(assay_ids)})
    return {row["sample_id"]: i for i, row in enumerate(metadata)}


def _groups_for(metadata: list[dict], matrix: dict,
                group_labels: list[str]) -> dict[str, list[int]]:
    """指定群ごとの行番号。選択規則は sample_manifest が唯一の定義。"""
    index_of = _row_index(metadata, matrix)
    selected = select_statistical_samples(metadata, group_labels)
    # 独立性（反復注入）の検査は検定の直前。ここを通らない限り n を数えない。
    validate_independent_samples(metadata, group_labels)

    rows: dict[str, list[int]] = {label: [] for label in group_labels}
    for row in selected:
        rows[row["group"]].append(index_of[row["sample_id"]])

    too_small = {label: len(idx) for label, idx in rows.items()
                 if len(idx) < _MIN_GROUP_N}
    if too_small:
        raise DomainError(
            "STATISTIC_GROUP_TOO_SMALL",
            f"n<{_MIN_GROUP_N} の群があります: "
            + ", ".join(f"{k}={v}" for k, v in sorted(too_small.items()))
            + "。n不足の群を落として別の検定にはしません。",
            {"group_n": {k: len(v) for k, v in rows.items()},
             "minimum": _MIN_GROUP_N})
    return rows


def _feature_indices(matrix: dict, specification: dict,
                     target_features: dict | None) -> list[tuple[int, str]]:
    """解析対象 feature の (列番号, feature_id)。eligibility と scope の AND。"""
    eligibility = np.asarray(matrix["eligibility_mask"], dtype=bool)
    feature_ids = list(matrix["feature_ids"])
    scope = specification.get("feature_scope") or {"mode": "all_eligible"}
    mode = scope.get("mode")

    if mode == "all_eligible":
        wanted = set(feature_ids)
    elif mode == "targets":
        resolved = target_features or {}
        missing = [t for t in (scope.get("target_ids") or []) if t not in resolved]
        if missing:
            raise DomainError(
                _SPEC_INVALID,
                f"feature_scope が未解決の target を指しています: {', '.join(missing)}"
                "（binding を先に解決してください）。",
                {"unresolved_target_ids": missing})
        wanted = {resolved[t] for t in (scope.get("target_ids") or [])}
    else:
        raise DomainError(
            _SPEC_INVALID, f"未知の feature_scope.mode です: {mode!r}",
            {"mode": mode})

    return [(i, fid) for i, fid in enumerate(feature_ids)
            if eligibility[i] and fid in wanted]


def _envelope(specification: dict, matrix: dict, **extra) -> dict:
    return {
        "schema": STATISTIC_SCHEMA,
        "statistic_id": specification.get("statistic_id"),
        "kind": specification.get("kind"),
        "matrix_id": matrix.get("matrix_id"),
        "matrix_recipe_id": (specification.get("matrix_recipe_id")
                             or matrix.get("recipe_id")),
        "transform": specification.get("transform", "none"),
        **extra,
    }


# ---------- 検定 ----------

def _run_welch(matrix, specification, metadata, indices) -> dict:
    reference = specification.get("reference_group")
    test = specification.get("test_group")
    if not reference or not test or reference == test:
        raise DomainError(
            _SPEC_INVALID,
            "reference_group と test_group を異なる値で明示してください。",
            {"reference_group": reference, "test_group": test})

    rows = _groups_for(metadata, matrix, [reference, test])
    values = np.asarray(matrix["values"], dtype=float)
    transformed = transform_values(values, specification.get("transform", "none"))

    features = []
    for index, feature_id in indices:
        raw_a = values[rows[reference], index]
        raw_b = values[rows[test], index]
        a = transformed[rows[reference], index]
        b = transformed[rows[test], index]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]

        entry = {"feature_id": feature_id, "status": "not_testable",
                 "reason": None, "p_value": None, "q_value": None,
                 "t_statistic": None,
                 "log2fc": arithmetic_log2fc(raw_a, raw_b),
                 "n_reference": int(a.size), "n_test": int(b.size)}
        if a.size < _MIN_GROUP_N or b.size < _MIN_GROUP_N:
            entry["reason"] = "insufficient_finite_values"
        elif a.var(ddof=1) == 0 and b.var(ddof=1) == 0:
            entry["reason"] = "zero_variance"
        else:
            t, p = differential.welch_t(a, b)
            if not math.isfinite(p):
                entry["reason"] = "test_statistic_not_finite"
            else:
                entry.update(status="tested", t_statistic=float(t),
                             p_value=float(p))
        if math.isnan(entry["log2fc"]):
            entry["log2fc"] = None
            entry["log2fc_reason"] = "non_positive_arithmetic_mean"
        features.append(entry)

    tested = [f for f in features if f["status"] == "tested"]
    for feature, q in zip(tested, differential.bh_fdr([f["p_value"] for f in tested])):
        feature["q_value"] = None if q is None or not math.isfinite(q) else float(q)

    status = "completed" if tested else "not_evaluable"
    return _envelope(
        specification, matrix,
        status=status,
        reason=None if tested else "no_testable_feature",
        effect_size_definition=EFFECT_SIZE_DEFINITION,
        reference_group=reference, test_group=test,
        groups={reference: len(rows[reference]), test: len(rows[test])},
        bh_population=len(tested),
        q_threshold=specification.get("q_threshold"),
        log2fc_threshold=specification.get("log2fc_threshold"),
        features=features)


def _run_anova(matrix, specification, metadata, indices) -> dict:
    labels = list(specification.get("groups") or [])
    if len(labels) < 3 or len(set(labels)) != len(labels):
        raise DomainError(
            _SPEC_INVALID,
            f"anova_tukey には重複のない3群以上が必要です: {labels}",
            {"groups": labels})

    rows = _groups_for(metadata, matrix, labels)
    values = np.asarray(matrix["values"], dtype=float)
    transformed = transform_values(values, specification.get("transform", "none"))
    alpha = specification.get("alpha", 0.05)

    features = []
    for index, feature_id in indices:
        columns = [transformed[rows[label], index] for label in labels]
        result = multigroup.test_feature(columns, alpha, labels=labels)
        features.append({"feature_id": feature_id, "q_value": None, **result})

    tested = [f for f in features if f["status"] == "tested"]
    for feature, q in zip(tested, differential.bh_fdr([f["p_value"] for f in tested])):
        feature["q_value"] = None if q is None or not math.isfinite(q) else float(q)

    status = "completed" if tested else "not_evaluable"
    return _envelope(
        specification, matrix,
        status=status,
        reason=None if tested else "no_testable_feature",
        groups={label: len(rows[label]) for label in labels},
        alpha=alpha,
        bh_population=len(tested),
        # Tukey の補正は feature 内の群対に閉じている。feature 横断の FDR
        # （ANOVA p の q 値）と同じ欄に並べない。
        tukey_correction_scope="within_feature_group_pairs",
        features=features)


def _run_pca(matrix, specification, metadata, indices) -> dict:
    from metabolomix.analysis.pca import run_pca

    # 学習対象は included の role=sample だけ。QC・blank・標準品は初期版では
    # 学習に入れない（spec §10）——測定の道具の分散を生物の分散として読ませない。
    keep = [i for i, row in enumerate(metadata)
            if row.get("include", True) and row.get("role") == "sample"]
    values = np.asarray(matrix["values"], dtype=float)
    transformed = transform_values(values, specification.get("transform", "none"))
    columns = [index for index, _ in indices]
    subset = transformed[np.ix_(keep, columns)] if keep and columns else np.empty((0, 0))

    finite_columns = [j for j in range(subset.shape[1])
                      if np.isfinite(subset[:, j]).all()
                      and float(np.std(subset[:, j])) > 0]
    n_excluded = subset.shape[1] - len(finite_columns)
    used = subset[:, finite_columns] if finite_columns else np.empty((len(keep), 0))

    if used.shape[0] < 2 or used.shape[1] < 2:
        return _envelope(
            specification, matrix, status="not_evaluable",
            reason="insufficient_samples_or_features",
            n_samples=int(used.shape[0]), n_features_used=int(used.shape[1]),
            n_features_excluded=int(n_excluded), features=[])

    requested = specification.get("n_components") or 5
    scaling = specification.get("scaling", "autoscale")
    if scaling not in ("none", "autoscale"):
        raise DomainError(
            _SPEC_INVALID, f"未知の scaling です: {scaling!r}",
            {"scaling": scaling, "supported": ["none", "autoscale"]})

    # PCA が学習する中心化・スケーリング後の行列で rank を制限する。
    from sklearn.preprocessing import StandardScaler

    fitted_values = StandardScaler(with_std=(scaling == "autoscale")).fit_transform(used)
    n_components = int(min(requested, np.linalg.matrix_rank(fitted_values),
                           used.shape[0] - 1, used.shape[1]))

    result = run_pca(used, n_components=n_components, log_transform=False,
                     scaling=scaling)
    return _envelope(
        specification, matrix, status="completed", reason=None,
        scaling=scaling,
        n_components=n_components,
        n_samples=int(used.shape[0]),
        n_features_used=int(used.shape[1]),
        n_features_excluded=int(n_excluded),
        sample_ids=[metadata[i]["sample_id"] for i in keep],
        explained_variance_ratio=[round(float(v), 6) for v in
                                  result["explained_variance_ratio"][:n_components]],
        scores=result["components"],
        features=[])


def run_statistic(matrix: dict, specification: dict, metadata: list[dict],
                  target_features: dict | None = None) -> dict:
    """1つの統計指定を、指定された解析行列に対して実行する。

    `metadata` は解析行列の assay 列と同じ順・同じ件数の sample-manifest 行。
    `target_features` は `feature_scope.mode="targets"` のときに要る
    `target_id -> feature_id`（binding の解決結果）。

    行列は `matrix_id` で名指しして結果に記録する。「直近の前処理」ではなく
    どの行列で出した数字かを、結果自身が持つ。
    """
    kind = specification.get("kind")
    if kind not in _KINDS:
        raise DomainError(
            _SPEC_INVALID,
            f"未知の統計種別です: {kind!r}（対応: {', '.join(_KINDS)}）",
            {"kind": kind, "supported": list(_KINDS)})

    indices = _feature_indices(matrix, specification, target_features)
    if kind == "welch":
        return _run_welch(matrix, specification, metadata, indices)
    if kind == "anova_tukey":
        return _run_anova(matrix, specification, metadata, indices)
    return _run_pca(matrix, specification, metadata, indices)
