"""前処理後の行列での PCA（統計エンジン spec §10）。

軸は生体の個体だけで決め、QC は同じ軸へ投影する。欠測は PCA のためだけに補う
（`impute_for_pca`）。`half_min` は基準行を生体の個体にし、QC の欠測も生体の個体から決めた
特徴ごとの値（最小の log2 − 1）で埋める。`knn` は個体と QC を積んだ行列の上で、各行（個体・QC）の
欠測をその行の近い特徴の値の平均で埋める（近さは全行で測り、基準行は使わない）。
計算本体は `analysis.pca.run_pca_fit_subset`（決定的にするため `svd_solver="full"`）。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis.pca import run_pca_fit_subset
from metabolomix.analysis.stats.preprocess import PreparedMatrix, impute_for_pca
from metabolomix.core.atomic_io import DomainError


def pca(prepared: PreparedMatrix, *, scaling: str = "autoscale",
        n_components: int | None = 5, impute: str = "half_min") -> dict:
    """生体の個体で軸を決めた PCA。`loadings` は成分 × 残した特徴（run_pca_fit_subset と同じ向き）。"""
    if scaling not in ("autoscale", "pareto"):
        raise DomainError("PCA_SPECIFICATION_INVALID", f"未知の scaling: {scaling!r}",
                          {"field": "scaling", "value": repr(scaling)})
    if impute not in ("half_min", "knn", "none"):
        raise DomainError("PCA_SPECIFICATION_INVALID", f"未知の補完法: {impute!r}",
                          {"field": "impute", "value": repr(impute)})
    if n_components is not None and (isinstance(n_components, bool)
                                     or not isinstance(n_components, (int, np.integer))
                                     or n_components < 1):
        raise DomainError("PCA_SPECIFICATION_INVALID", "n_components は 1 以上の整数か None にしてください",
                          {"field": "n_components", "value": repr(n_components)})
    units = np.asarray(prepared.log2, dtype=float)
    n_units = units.shape[0]
    qc = np.asarray(prepared.qc_log2, dtype=float)
    if qc.size == 0:
        qc = np.empty((0, units.shape[1]))
    stacked = np.vstack([units, qc])
    n_nan = int((~np.isfinite(stacked)).sum())
    if n_nan:
        if impute == "none":
            raise DomainError("PCA_MISSING_VALUES", "欠測があるため PCA を計算できません（補完なし）",
                              {"missing_cells": n_nan})
        try:
            stacked, rep = impute_for_pca(stacked, method=impute, reference_rows=np.arange(n_units))
        except ValueError as exc:
            raise DomainError("PCA_SPECIFICATION_INVALID", str(exc), {"field": "impute"}) from exc
        imputation = {"method": rep["method"], "imputed_cells": rep["imputed_cells"], "used": True}
    else:
        imputation = {"method": impute, "imputed_cells": 0, "used": False}
    still = np.flatnonzero(~np.isfinite(stacked).all(axis=0))
    if still.size:
        raise DomainError("PCA_MISSING_VALUES", "補完後も欠測が残る特徴があります",
                          {"columns": [prepared.feature_ids[j] for j in still]})
    if n_units < 3:
        raise DomainError("STATISTIC_GROUP_TOO_SMALL", "PCA の計算に使う個体が 3 未満です",
                          {"n_fit": n_units})
    fit_rows = stacked[:n_units]
    mean, sd = fit_rows.mean(axis=0), fit_rows.std(axis=0)
    kept = int((sd > 1e-12 * np.maximum(1.0, np.abs(mean))).sum())
    if kept < 2:
        raise DomainError("PCA_TOO_FEW_FEATURES", "値が変わる特徴が 2 未満のため PCA を計算できません",
                          {"kept": kept})
    fit = np.zeros(stacked.shape[0], dtype=bool)
    fit[:n_units] = True
    try:
        res = run_pca_fit_subset(stacked, fit, n_components, scaling=scaling, svd_solver="full")
    except ValueError as exc:
        raise DomainError("PCA_SPECIFICATION_INVALID", str(exc), {}) from exc
    comp = res["components"]
    return {
        "schema": "stats-pca.v1",
        "matrix_id": prepared.matrix_id,
        "scaling": scaling,
        "log2": True,
        "imputation": imputation,
        "n_fit": res["n_fit"],
        "unit_ids": [u.unit_id for u in prepared.units],
        "qc_ids": list(prepared.qc_ids),
        "scores": {"units": comp[:n_units], "qc": comp[n_units:]},
        "explained_variance_ratio": res["explained_variance_ratio"],
        "loadings": res["loadings"],
        "kept_features": [prepared.feature_ids[j] for j in res["kept_features"]],
    }
