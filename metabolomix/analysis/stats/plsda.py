"""PLS-DA と VIP（統計エンジン。analysis-spec v2 の解析の種類 `plsda`）。

前処理済みの個体 × 特徴の log2（`PreparedMatrix.log2`）で、`design.factors` の水準の組（cell）を
クラスにした PLS 判別分析を当てはめ、次を**必ず**返す（過学習を数字で見えるようにするため）:

- 成分数（`n_components`。個体数と特徴数から決まる上限で頭打ちにしたら警告）
- 交差検証の Q2（層化 k 分割。`q2_cum` は成分 1〜A のそれぞれ。`Q2 = 1 − PRESS / TSS`、TSS は全個体の
  Y の平均の周りの平方和）。R2X・R2Y（当てはめた全個体で）も並べる
- 並べ替え検定の p（クラスの札を `n_permutations` 回入れ替え、同じ分割の交差検証で Q2 を出し直した
  帰無分布で `(1 + #{Q2_perm ≥ Q2}) / (1 + n_permutations)`）
- VIP（全個体の当てはめの変数重要度。`VIP_j = √(p·Σ_a SSY_a·(w_ja/‖w_a‖)² / Σ_a SSY_a)`、SSY_a は成分 a が
  説明する Y の平方和）

Y はクラスの指示変数（2 クラスなら後のクラスを 1 とした 1 列、3 クラス以上は 1 クラス 1 列）。X は
`scaling`（`autoscale`＝平均 0・SD 1 / `pareto`＝SD の平方根で割る）で尺度をそろえ、交差検証では
**学習の個体だけ**で中心と尺度を決める。欠測は PCA と同じ規則（`impute_for_pca`）で、検定前に 1 回だけ補う
（交差検証の外。`imputation` に件数）。全個体で値が変わらない特徴は外す（`n_features_dropped_constant`）。
QC は入れない（`PreparedMatrix.log2` は生体の個体だけ）。乱数（分割と並べ替え）は `seed` で決まる。

MCP・セッションには依存しない。数字は丸めない（非有限は None）。
"""
from __future__ import annotations

import math
import warnings as _warnings

import numpy as np

from metabolomix.analysis.stats.preprocess import PreparedMatrix, impute_for_pca
from metabolomix.core.atomic_io import DomainError

__all__ = ["SCHEMA", "SCALINGS", "IMPUTES", "plsda"]

SCHEMA = "stats-plsda.v1"
SCALINGS = ("autoscale", "pareto")
IMPUTES = ("half_min", "knn", "none")
_MIN_UNITS = 2
_SD_FLOOR = 1e-12


def _invalid(message: str, **details) -> DomainError:
    return DomainError("PLSDA_SPECIFICATION_INVALID", message, details)


def _num(x) -> float | None:
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def _label(unit, factors) -> str | None:
    values = []
    for f in factors:
        if f in ("group", "batch"):
            raw = getattr(unit, f, None)
        else:
            raw = (getattr(unit, "attributes", None) or {}).get(f)
        if raw is None or str(raw) == "":
            return None
        values.append(str(raw))
    return ":".join(values)


def _scale_params(x: np.ndarray, scaling: str) -> tuple[np.ndarray, np.ndarray]:
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=1) if x.shape[0] > 1 else np.ones(x.shape[1])
    sd = np.where(sd > _SD_FLOOR, sd, 1.0)
    return mean, (sd if scaling == "autoscale" else np.sqrt(sd))


def _fit(x: np.ndarray, y: np.ndarray, n_components: int):
    from sklearn.cross_decomposition import PLSRegression

    model = PLSRegression(n_components=n_components, scale=False)
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        model.fit(x, y)
    return model


def _predict_by_components(model, x_centered: np.ndarray, y_mean: np.ndarray, n_components: int):
    """成分 1〜A それぞれの予測（A 成分の当てはめの先頭 a 成分＝a 成分の当てはめ。NIPALS の P'W は上三角）。"""
    rot, q = model.x_rotations_, model.y_loadings_
    scores = x_centered @ rot
    return [y_mean + scores[:, :a] @ q[:, :a].T for a in range(1, n_components + 1)]


def _folds(labels: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """層化 k 分割の割り当て（クラスごとに並べ替えて順に配る）。"""
    fold = np.empty(labels.size, dtype=int)
    start = 0
    for cls in sorted(set(labels.tolist())):
        idx = np.flatnonzero(labels == cls)
        idx = idx[rng.permutation(idx.size)]
        fold[idx] = (np.arange(idx.size) + start) % k
        start += idx.size
    return fold


def _cv_q2(x: np.ndarray, y: np.ndarray, fold: np.ndarray, k: int, n_components: int,
           scaling: str) -> list[float]:
    """交差検証の Q2（成分 1〜A）。分割ごとに学習の個体だけで中心と尺度を決めて当てはめる。"""
    press = np.zeros(n_components)
    for f in range(k):
        test = fold == f
        if not test.any():
            continue
        train = ~test
        mean, scale = _scale_params(x[train], scaling)
        xs_train = (x[train] - mean) / scale
        xs_test = (x[test] - mean) / scale
        model = _fit(xs_train, y[train], n_components)
        x_center = xs_train.mean(axis=0)
        preds = _predict_by_components(model, xs_test - x_center, y[train].mean(axis=0), n_components)
        for a, pred in enumerate(preds):
            press[a] += float(((y[test] - pred) ** 2).sum())
    tss = float(((y - y.mean(axis=0)) ** 2).sum())
    return [1.0 - p / tss if tss > 0 else math.nan for p in press]


def _vip(model) -> np.ndarray:
    w = model.x_weights_                       # p × A（列は単位長）
    t = model.x_scores_                        # n × A
    q = model.y_loadings_                      # m × A
    ssy = (t ** 2).sum(axis=0) * (q ** 2).sum(axis=0)
    total = ssy.sum()
    if not total > 0:
        return np.full(w.shape[0], np.nan)
    wn = w / np.linalg.norm(w, axis=0, keepdims=True)
    return np.sqrt(w.shape[0] * ((wn ** 2) @ ssy) / total)


def plsda(prepared: PreparedMatrix, *, factors: list[str], classes: list[str] | None = None,
          n_components: int = 2, scaling: str = "autoscale", impute: str = "half_min",
          cv_folds: int = 7, n_permutations: int = 200, seed: int = 0) -> dict:
    """PLS-DA を当てはめ、成分数・交差検証の Q2・並べ替え検定の p・VIP を返す。"""
    if scaling not in SCALINGS:
        raise _invalid(f"scaling は {SCALINGS} のどれか", field="scaling", value=repr(scaling))
    if impute not in IMPUTES:
        raise _invalid(f"impute は {IMPUTES} のどれか", field="impute", value=repr(impute))
    for name, value, low in (("n_components", n_components, 1), ("cv_folds", cv_folds, 2),
                             ("n_permutations", n_permutations, 1), ("seed", seed, 0)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < low:
            raise _invalid(f"{name} は {low} 以上の整数", field=name, value=repr(value))
    factors = list(factors)
    warnings: list[dict] = []

    labels_all = [_label(u, factors) for u in prepared.units]
    lacking = [u.unit_id for u, lab in zip(prepared.units, labels_all) if lab is None]
    if lacking:
        warnings.append({"code": "units_missing_model_values",
                         "message": f"{', '.join(factors)} の値が欠けた {len(lacking)} 個体を外しました。",
                         "unit_ids": lacking})
    present = sorted({lab for lab in labels_all if lab is not None})
    if classes is not None:
        unknown = [c for c in classes if c not in present]
        if unknown:
            raise _invalid(f"classes に無い cell がある: {unknown}", field="classes",
                           unknown=unknown, available=present)
        order = list(dict.fromkeys(classes))
    else:
        order = present
    keep = [i for i, lab in enumerate(labels_all) if lab is not None and lab in order]
    labels = np.array([labels_all[i] for i in keep])
    counts = {c: int((labels == c).sum()) for c in order}
    if len(order) < 2 or any(n < _MIN_UNITS for n in counts.values()):
        raise DomainError("STATISTIC_GROUP_TOO_SMALL",
                          f"PLS-DA には {_MIN_UNITS} 個体以上のクラスが 2 つ以上要る",
                          {"class_n": counts})

    x = np.asarray(prepared.log2, dtype=float)[keep]
    n_missing = int((~np.isfinite(x)).sum())
    if n_missing:
        if impute == "none":
            raise DomainError("PLSDA_MISSING_VALUES", "欠測があるため PLS-DA を計算できません（補完なし）",
                              {"missing_cells": n_missing})
        x, rep = impute_for_pca(x, method=impute)
        imputation = {"method": rep["method"], "imputed_cells": rep["imputed_cells"], "used": True}
    else:
        imputation = {"method": impute, "imputed_cells": 0, "used": False}
    finite_cols = np.isfinite(x).all(axis=0)
    with np.errstate(invalid="ignore"):
        sd = np.where(finite_cols, np.nanstd(np.where(np.isfinite(x), x, 0.0), axis=0), 0.0)
        mean = np.where(finite_cols, np.nanmean(np.where(np.isfinite(x), x, 0.0), axis=0), 0.0)
    kept = finite_cols & (sd > _SD_FLOOR * np.maximum(1.0, np.abs(mean)))
    if int(kept.sum()) < 2:
        raise DomainError("PLSDA_TOO_FEW_FEATURES", "値が変わる特徴が 2 未満のため PLS-DA を計算できません",
                          {"kept": int(kept.sum())})
    x = x[:, kept]
    kept_ids = [fid for fid, k in zip(prepared.feature_ids, kept) if k]
    n, p = x.shape

    y = (np.column_stack([(labels == order[1]).astype(float)]) if len(order) == 2
         else np.column_stack([(labels == c).astype(float) for c in order]))
    k = min(int(cv_folds), n)
    if k < cv_folds:
        warnings.append({"code": "cv_folds_capped",
                         "message": f"分割数を個体数 {n} に合わせて {k} にしました（1 個抜き）。"})
    smallest_train = n - int(math.ceil(n / k))
    cap = max(1, min(smallest_train - 1, p))
    a = min(int(n_components), cap)
    if a < n_components:
        warnings.append({"code": "n_components_capped",
                         "message": f"成分数を {n_components} から {a} にしました（交差検証の学習の個体数と"
                                    "特徴数から決まる上限）。", "requested": int(n_components),
                         "used": a})

    rng = np.random.default_rng(seed)
    fold = _folds(labels, k, rng)
    q2_cum = _cv_q2(x, y, fold, k, a, scaling)

    mean_all, scale_all = _scale_params(x, scaling)
    xs = (x - mean_all) / scale_all
    model = _fit(xs, y, a)
    preds = _predict_by_components(model, xs - xs.mean(axis=0), y.mean(axis=0), a)
    yss = float(((y - y.mean(axis=0)) ** 2).sum())
    r2y = [1.0 - float(((y - pr) ** 2).sum()) / yss for pr in preds]
    xss = float((xs ** 2).sum())
    t, pl = model.x_scores_, model.x_loadings_
    r2x = [1.0 - float(((xs - xs.mean(axis=0) - t[:, :c] @ pl[:, :c].T) ** 2).sum()) / xss
           for c in range(1, a + 1)]
    vip = _vip(model)

    q2 = q2_cum[-1]
    null = np.empty(int(n_permutations))
    for i in range(int(n_permutations)):
        perm = rng.permutation(n)
        null[i] = _cv_q2(x, y[perm], fold, k, a, scaling)[-1]
    finite_null = null[np.isfinite(null)]
    exceed = int((finite_null >= q2).sum()) if math.isfinite(q2) else len(finite_null)
    perm_p = (1 + exceed) / (1 + len(finite_null)) if len(finite_null) else None
    if not (math.isfinite(q2) and q2 > 0):
        warnings.append({"code": "plsda_q2_not_positive",
                         "message": f"交差検証の Q2 が {q2:.3g} で正でありません（クラスを予測できていない）。"
                                    "VIP は当てはめの重要度で、予測できない模型の VIP は判別の根拠になりません。"})

    order_vip = np.argsort(-np.nan_to_num(vip, nan=-1.0), kind="stable")
    return {
        "schema": SCHEMA,
        "matrix_id": prepared.matrix_id,
        "factors": factors,
        "classes": counts,
        "y_coding": ("one_column_second_class" if len(order) == 2 else "indicator_per_class"),
        "n_units": int(n),
        "unit_ids": [prepared.units[i].unit_id for i in keep],
        "scaling": scaling,
        "imputation": imputation,
        "n_features": int(p),
        "n_features_dropped_constant": int(len(prepared.feature_ids) - p),
        "n_components_requested": int(n_components),
        "n_components": int(a),
        "r2x_cum": [_num(v) for v in r2x],
        "r2y_cum": [_num(v) for v in r2y],
        "q2_cum": [_num(v) for v in q2_cum],
        "q2": _num(q2),
        "cv": {"method": "stratified_kfold", "folds": int(k), "seed": int(seed),
               "scaling_within_fold": True, "q2_definition": "1 - PRESS/TSS"},
        "permutation": {"n": int(n_permutations), "n_finite": int(len(finite_null)),
                        "p": _num(perm_p), "statistic": "q2",
                        "q2_null_mean": _num(finite_null.mean()) if len(finite_null) else None,
                        "q2_null_q95": (_num(np.quantile(finite_null, 0.95))
                                        if len(finite_null) else None)},
        "vip": [{"feature_id": kept_ids[j], "vip": _num(vip[j])} for j in order_vip],
        "scores": {"unit_ids": [prepared.units[i].unit_id for i in keep],
                   "values": [[_num(v) for v in row] for row in t]},
        "warnings": warnings,
    }
