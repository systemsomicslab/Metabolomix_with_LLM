import numpy as np
import pytest


def test_arf_reader_does_not_reexport_run_pca():
    """v1 の `run_pca` は D 段で消えた。arf/reader.py に後方互換の束縛を残さない
    （ARF のツールとそのモックは B-2b で消えた）。"""
    from metabolomix.arf import reader
    assert not hasattr(reader, "run_pca")


def test_run_pca_fit_subset_projects_unfitted_rows_and_drops_constant_columns():
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    rng = np.random.default_rng(0)
    fit_rows = rng.normal(size=(6, 4))
    fit_rows[:, 3] = 5.0                                  # 計算に使う試料で一定 → 外す
    outlier = np.array([[10.0, -10.0, 3.0, 99.0]])
    matrix = np.vstack([fit_rows, outlier])
    fit = np.array([True] * 6 + [False])
    out = run_pca_fit_subset(matrix, fit, n_components=2)
    assert out["kept_features"] == [0, 1, 2] and out["n_fit"] == 6
    assert len(out["components"]) == 7                    # 投影した試料にもスコアがある
    # 計算に使った試料だけで決めた軸: 外れ値を入れても fit 試料のスコアは変わらない
    alone = run_pca_fit_subset(fit_rows, np.ones(6, bool), n_components=2)
    assert np.allclose(np.abs(out["components"][:6]), np.abs(alone["components"]), atol=1e-9)


def test_loading_correlations_equal_pearson_r_on_fitted_rows():
    import numpy as np
    from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset

    rng = np.random.default_rng(1)
    matrix = rng.normal(size=(8, 5))
    fit = np.array([True] * 7 + [False])
    out = run_pca_fit_subset(matrix, fit, n_components=2)
    r = np.asarray(loading_correlations(out["loadings"], out["singular_values"], out["n_fit"]))
    scores = np.asarray(out["components"])[fit]
    for k in range(2):
        for j in range(5):
            expected = np.corrcoef(matrix[fit][:, j], scores[:, k])[0, 1]
            assert r[k, j] == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("fit", [[True, True, False, False], [True, True, True, False]])
def test_run_pca_fit_subset_needs_three_fitted_rows_and_two_varying_columns(fit):
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    matrix = np.array([[1.0, 2.0], [2.0, 2.0], [3.0, 2.0], [4.0, 2.0]])   # 列 1 は一定
    with pytest.raises(ValueError):
        run_pca_fit_subset(matrix, np.array(fit))


def test_run_pca_fit_subset_drops_float_noise_constant_columns():
    """定数列の母標準偏差は浮動小数点の誤差で 0 にならないことがある。相対閾値で外す。"""
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    rng = np.random.default_rng(2)
    varying = rng.normal(size=(7, 2))
    matrix = np.hstack([varying[:, :1], np.full((7, 1), 0.1),
                        varying[:, 1:], np.full((7, 1), 1234.567)])
    out = run_pca_fit_subset(matrix, np.ones(7, bool), n_components=2)
    assert out["kept_features"] == [0, 2]


def test_run_pca_fit_subset_rejects_non_positive_n_components():
    import numpy as np
    from metabolomix.analysis.pca import run_pca_fit_subset

    matrix = np.random.default_rng(3).normal(size=(6, 4))
    with pytest.raises(ValueError):
        run_pca_fit_subset(matrix, np.ones(6, bool), n_components=0)
