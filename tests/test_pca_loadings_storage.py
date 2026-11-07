"""既存の PCA（ARF / mzTab）がローディング図に要る情報をセッションに残すこと。"""
from __future__ import annotations

import numpy as np



def test_run_dataset_pca_returns_singular_values():
    from metabolomix.analysis.dataset_analysis import run_dataset_pca

    class DS:
        pp_matrix = np.array([[1.0, 2.0, 3.0], [2.0, 1.0, 0.5], [3.0, 5.0, 1.0], [0.5, 0.2, 2.0]])
        pp_sample_names = ["s1", "s2", "s3", "s4"]
        pp_feature_names = ["f1", "f2", "f3"]
        roles = {}
        sample_meta = {}

    out = run_dataset_pca(DS(), n_components=2)
    assert len(out["singular_values"]) == 2 and len(out["loadings"]) == 2
