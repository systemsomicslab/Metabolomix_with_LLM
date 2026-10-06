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
