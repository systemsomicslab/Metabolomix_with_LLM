"""前処理後の行列での PCA（統計エンジン spec §10）。

軸は生体の個体だけで決め、QC は同じ軸へ投影する。欠測は PCA のためだけに補う
（`impute_for_pca`。基準行は生体の個体で、QC も同じ特徴ごとの値で埋める）。
計算本体は `analysis.pca.run_pca_fit_subset`。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis.pca import run_pca_fit_subset
from metabolomix.analysis.stats.preprocess import PreparedMatrix, impute_for_pca
from metabolomix.core.atomic_io import DomainError


def pca(prepared: PreparedMatrix, *, scaling: str = "autoscale",
        n_components: int | None = 5, impute: str = "half_min") -> dict:
    """生体の個体で軸を決めた PCA。`loadings` は成分 × 残した特徴（run_pca_fit_subset と同じ向き）。"""
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
        stacked, rep = impute_for_pca(stacked, method=impute, reference_rows=np.arange(n_units))
        imputation = {"method": rep["method"], "imputed_cells": rep["imputed_cells"], "used": True}
    else:
        imputation = {"method": impute, "imputed_cells": 0, "used": False}
    fit = np.zeros(stacked.shape[0], dtype=bool)
    fit[:n_units] = True
    try:
        res = run_pca_fit_subset(stacked, fit, n_components, scaling=scaling)
    except ValueError as exc:
        if n_units < 3:
            raise DomainError("STATISTIC_GROUP_TOO_SMALL", "PCA の計算に使う個体が 3 未満です",
                              {"n_fit": n_units}) from exc
        raise
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
