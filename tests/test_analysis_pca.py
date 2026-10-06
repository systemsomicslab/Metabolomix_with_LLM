import numpy as np
import pytest


def test_run_pca_importable_from_analysis():
    from metabolomix.analysis.pca import run_pca
    rng = np.random.default_rng(0)
    result = run_pca(rng.random((6, 10)), n_components=2)
    assert set(result) == {"components", "explained_variance_ratio",
                           "singular_values", "loadings"}
    assert len(result["components"]) == 6
    assert len(result["explained_variance_ratio"]) == 2


def test_arf_reader_reexports_the_same_object():
    """arf/reader.py の束縛が analysis/pca.py の関数と同一であること。

    既存 ARF テストが patch.object(server.arf_reader, "run_pca", ...) で
    差し替えるため、この束縛が消えるとモックが効かなくなる。
    """
    from metabolomix.analysis.pca import run_pca as canonical
    from metabolomix.arf.reader import run_pca as reexported
    assert reexported is canonical


def test_run_pca_rejects_degenerate_matrix():
    from metabolomix.analysis.pca import run_pca
    with pytest.raises(ValueError):
        run_pca(np.array([[1.0, 2.0, 3.0]]))  # n_samples=1
