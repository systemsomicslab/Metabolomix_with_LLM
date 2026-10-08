"""PCA（入力形式に依存しない数値処理）。

もとは metabolomix/arf/reader.py にあった run_pca をここへ移した。ARF と
DatasetState（mzTab-M）の双方から**同一実装**を呼ぶため、形式非依存の
analysis/ 層に置く。analysis/ が arf/ を import する層序の逆転を作らない。

metabolomix/arf/reader.py はトップレベルで同名を再エクスポートしている。
既存 ARF テストが patch.object(server.arf_reader, "run_pca", ...) で
モジュール属性を差し替えるため、その束縛は消してはならない。

依存は numpy + sklearn のみ。metabolomix.* を import しない leaf。
"""
from __future__ import annotations

import math

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(matrix: np.ndarray, n_components: int | None = None,
            log_transform: bool = False, *, scaling: str = "autoscale") -> dict:
    """
    標準化（スケーリング）を行った上で多次元PCAを実行する

    引数:
      log_transform: True のとき標準化の前に log10 変換を適用する（既定 False=従来動作）。
        ピーク強度は右に大きく歪むため、log変換で正規性が改善し条件分離が向上しやすい。
        0以下/極小値対策として下限1.0でクリップしてから log10 を取る。
      scaling: `"autoscale"`（既定・従来動作＝平均0分散1へ標準化）または `"none"`
        （中心化だけ行い分散で割らない）。v2（spec §10）は profile が
        `scaling` を明示するため、既定値へ黙って倒さない経路が要る。
        `"none"` でも中心化はする——中心化しない PCA の第1成分は「平均」を
        表すだけで、群構造を見るという目的に合わない。
    """
    n_samples, n_features = matrix.shape
    max_components = min(n_samples, n_features)

    if max_components < 2:
        raise ValueError(f"PCAには2つ以上のサンプルと特徴量が必要です。（現在: サンプル={n_samples}, 特徴量={n_features}）")

    # n_componentsが指定されていない場合は最大次元数まで計算
    target_components = max_components if n_components is None else min(n_components, max_components)

    work_matrix = matrix
    # 任意: log変換（既定オフ）。強度の歪みを抑える。
    if log_transform:
        work_matrix = np.log10(np.clip(matrix, 1.0, None))

    # 【重要】スケールの異なる多変量（Height, RTなど）を扱うための標準化
    if scaling not in ("autoscale", "none"):
        raise ValueError(f"unknown scaling: {scaling!r}")
    scaler = StandardScaler(with_std=(scaling == "autoscale"))
    scaled_matrix = scaler.fit_transform(work_matrix)

    pca = PCA(n_components=target_components)
    transformed = pca.fit_transform(scaled_matrix)

    return {
        "components": transformed.tolist(),
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "singular_values": pca.singular_values_.tolist(),
        "loadings": pca.components_.tolist()
    }


def run_pca_fit_subset(matrix, fit_mask, n_components: int | None = None) -> dict:
    """計算に使う試料（fit_mask=True）だけで autoscale と PCA の軸を決め、全試料をその軸へ投影する。

    低信頼の試料を主成分の計算から外しつつ、どこに落ちるかは見たいときに使う。autoscale は
    計算に使う試料の平均と**母標準偏差**（ddof=0。StandardScaler と同じ）。計算に使う試料で
    分散が 0 の列は外し、`kept_features` に残した列の番号を返す。
    """
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
    scaled = (matrix[:, kept] - mean) / sd
    max_components = min(n_fit, kept.size)
    target = max_components if n_components is None else min(n_components, max_components)
    pca = PCA(n_components=target).fit(scaled[fit])
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
