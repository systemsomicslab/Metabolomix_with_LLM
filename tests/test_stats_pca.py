"""前処理後の行列での PCA（統計エンジン spec §10）と run_pca_fit_subset の scaling 引数。"""
from __future__ import annotations

import numpy as np
import pytest
from sklearn.decomposition import PCA

from metabolomix.analysis.pca import run_pca_fit_subset
from metabolomix.analysis.stats.pca import pca
from metabolomix.analysis.stats.preprocess import PreparedMatrix, UnitInfo
from metabolomix.core.atomic_io import DomainError


def _prep(log2, qc=None) -> PreparedMatrix:
    log2 = np.asarray(log2, dtype=float)
    n, g = log2.shape
    qc = np.empty((0, g)) if qc is None else np.asarray(qc, dtype=float)
    units = [UnitInfo(unit_id=f"u{i}", group="A" if i % 2 else "B", batch="b1", subject_id=None,
                      injection_order=None, n_injections=1, sample_ids=[f"u{i}_inj0"])
             for i in range(n)]
    return PreparedMatrix(
        matrix_id="mx_pca", log2=log2, units=units, feature_ids=[f"f{j}" for j in range(g)],
        detected_units=None, qc_log2=qc, qc_ids=[f"qc{i}" for i in range(qc.shape[0])],
        steps=[], config={}, quant=None, dataset_hash=None, subject_source="injection",
        warnings=[], feature_classes={}, annotated_features=frozenset())


def _data(n=8, g=6, seed=0):
    return np.random.default_rng(seed).normal(20.0, 1.0, (n, g)) * np.arange(1, g + 1)


def test_run_pca_fit_subset_default_unchanged():
    x = _data()
    mask = np.ones(len(x), bool)
    a = run_pca_fit_subset(x, mask, 3)
    b = run_pca_fit_subset(x, mask, 3, scaling="autoscale")
    assert a == b
    sc = (x - x.mean(0)) / x.std(0)
    ref = PCA(n_components=3).fit(sc)
    assert a["components"] == ref.transform(sc).tolist()


def test_pareto_scaling_divides_by_sqrt_sd():
    x = _data()
    mask = np.ones(len(x), bool)
    out = run_pca_fit_subset(x, mask, 3, scaling="pareto")
    sc = (x - x.mean(0)) / np.sqrt(x.std(0))
    ref = PCA(n_components=3).fit(sc)
    np.testing.assert_allclose(out["components"], ref.transform(sc), atol=1e-12)
    np.testing.assert_allclose(out["explained_variance_ratio"], ref.explained_variance_ratio_)
    with pytest.raises(ValueError):
        run_pca_fit_subset(x, mask, 3, scaling="bogus")


def test_axes_fit_on_biological_units_and_qc_projected():
    x = _data()
    r1 = pca(_prep(x, qc=x[:2]), n_components=3)
    r2 = pca(_prep(x, qc=x[:2] + 100.0), n_components=3)
    assert r1["loadings"] == r2["loadings"]
    assert r1["scores"]["units"] == r2["scores"]["units"]
    assert r1["scores"]["qc"] != r2["scores"]["qc"]
    assert len(r1["scores"]["qc"]) == 2 and len(r1["scores"]["units"]) == 8
    assert r1["n_fit"] == 8 and r1["qc_ids"] == ["qc0", "qc1"]
    assert r1["kept_features"] == [f"f{j}" for j in range(6)]
    assert r1["schema"] == "stats-pca.v1" and r1["log2"] is True
    assert pca(_prep(x))["scores"]["qc"] == []


def test_imputation_recorded_when_used():
    x = _data()
    x[1, 2] = np.nan
    qc = x[:2].copy()  # 行 1 の欠測も写る（計 3 セル）
    qc[0, 3] = np.nan
    r = pca(_prep(x, qc=qc), n_components=2)
    assert r["imputation"] == {"method": "half_min", "imputed_cells": 3, "used": True}
    assert np.isfinite(np.asarray(r["scores"]["qc"])).all()


def test_no_imputation_when_complete():
    r = pca(_prep(_data()), n_components=2)
    assert r["imputation"] == {"method": "half_min", "imputed_cells": 0, "used": False}


def test_impute_none_with_missing_raises():
    x = _data()
    x[0, 0] = np.nan
    with pytest.raises(DomainError) as e:
        pca(_prep(x), impute="none")
    assert e.value.code == "PCA_MISSING_VALUES"


def test_too_few_units():
    with pytest.raises(DomainError) as e:
        pca(_prep(_data(n=2)))
    assert e.value.code == "STATISTIC_GROUP_TOO_SMALL"
    assert e.value.details["n_fit"] == 2


def test_invalid_scaling_and_n_components():
    for kw in ({"scaling": "bogus"}, {"n_components": 0}, {"n_components": True}, {"n_components": 1.5}):
        with pytest.raises(DomainError) as e:
            pca(_prep(_data()), **kw)
        assert e.value.code == "PCA_SPECIFICATION_INVALID"


def test_too_few_varying_features():
    x = np.ones((6, 4))
    x[:, 0] = np.arange(6)
    with pytest.raises(DomainError) as e:
        pca(_prep(x))
    assert e.value.code == "PCA_TOO_FEW_FEATURES" and e.value.details["kept"] == 1


def test_all_missing_column_after_imputation():
    x = _data()
    x[:, 2] = np.nan
    with pytest.raises(DomainError) as e:
        pca(_prep(x))
    assert e.value.code == "PCA_MISSING_VALUES" and e.value.details["columns"] == ["f2"]
