"""PCA（入力形式に依存しない数値処理）。

計算に使う試料だけで軸を決めて全試料を投影する `run_pca_fit_subset` と、ローディングを相関 r に
換算する `loading_correlations`。統計エンジンの PCA（`analysis/stats/pca.py`）と、選んだ分子種の PCA
（`tools/plot_species_kinds.py`）が使う。

依存は numpy + sklearn のみ。metabolomix.* を import しない leaf。
"""
from __future__ import annotations

import math

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca_fit_subset(matrix, fit_mask, n_components: int | None = None,
                       scaling: str = "autoscale", svd_solver: str = "auto") -> dict:
    """計算に使う試料（fit_mask=True）だけで autoscale と PCA の軸を決め、全試料をその軸へ投影する。

    低信頼の試料を主成分の計算から外しつつ、どこに落ちるかは見たいときに使う。autoscale は
    計算に使う試料の平均と**母標準偏差**（ddof=0。StandardScaler と同じ）。計算に使う試料で
    分散が 0 の列は外し、`kept_features` に残した列の番号を返す。

    `scaling="pareto"` は中心化のあと SD の**平方根**で割る（SD・定数列の扱いは autoscale と同じ）。

    `svd_solver` は sklearn の PCA にそのまま渡す。既定の "auto" は大きな行列で乱数を使う
    randomized SVD を選び、実行ごとに結果が揺れうる。決定的な結果が要る呼び出しは "full" を渡す。
    """
    if scaling not in ("autoscale", "pareto"):
        raise ValueError(f"unknown scaling: {scaling!r}")
    matrix = np.asarray(matrix, dtype=float)
    fit = np.asarray(fit_mask, dtype=bool)
    if matrix.ndim != 2 or fit.shape != (matrix.shape[0],):
        raise ValueError("matrix は 2 次元、fit_mask は行数と同じ長さにしてください。")
    n_fit = int(fit.sum())
    if n_fit < 3:
        raise ValueError(f"PCA の計算に使う試料は 3 以上必要です（現在: {n_fit}）。")
    if n_components is not None and (
            isinstance(n_components, bool)
            or not isinstance(n_components, (int, np.integer))
            or n_components < 1):
        raise ValueError(f"n_components は 1 以上の整数にしてください（現在: {n_components!r}）。")
    fit_rows = matrix[fit]
    mean_all = fit_rows.mean(axis=0)
    sd_all = fit_rows.std(axis=0)
    # 定数列の母標準偏差は浮動小数点の誤差で 0 にならないことがあるので、相対閾値で外す。
    kept = np.flatnonzero(sd_all > 1e-12 * np.maximum(1.0, np.abs(mean_all)))
    if kept.size < 2:
        raise ValueError(f"計算に使う試料で値が変わる特徴量が 2 未満のため PCA を計算できません（現在: {kept.size}）。")
    mean = mean_all[kept]
    sd = sd_all[kept]
    scaled = (matrix[:, kept] - mean) / (sd if scaling == "autoscale" else np.sqrt(sd))
    max_components = min(n_fit, kept.size)
    target = max_components if n_components is None else min(n_components, max_components)
    pca = PCA(n_components=target, svd_solver=svd_solver).fit(scaled[fit])
    return {
        "components": pca.transform(scaled).tolist(),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "singular_values": pca.singular_values_.tolist(),
        "loadings": pca.components_.tolist(),
        "kept_features": kept.tolist(),
        "n_fit": n_fit,
    }


def loading_correlations(loadings, singular_values, n_fit: int) -> list[list[float]]:
    """各特徴量と主成分スコアの相関 r（主成分 × 特徴量）。

    r_jk = 固有ベクトルの成分_jk × 特異値_k / √n。母標準偏差で autoscale した PCA でだけ
    成り立つ（n は主成分の計算に使った試料数）。中心化だけの PCA には使わない。
    """
    coef = np.asarray(loadings, dtype=float)
    s = np.asarray(singular_values, dtype=float)[: coef.shape[0]]
    return (coef * s[:, None] / math.sqrt(n_fit)).tolist()
