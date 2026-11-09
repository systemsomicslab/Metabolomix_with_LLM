import numpy as np
import pytest

from metabolomix.mztab.dataset_state import DatasetState


def _make_ds(n_features=20, n_samples=8, with_blank=False):
    """最小限の DatasetState。ロールは dataset_analysis がサンプル名から推定する。"""
    rng = np.random.default_rng(42)
    ds = DatasetState()
    names = [f"ctrl_{i}" for i in range(n_samples // 2)]
    names += [f"treat_{i}" for i in range(n_samples - len(names))]
    if with_blank:
        names[-1] = "blank_1"
    ds.feature_matrix = rng.random((n_features, n_samples)) * 1000.0
    ds.sample_names = names
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    return ds


def _ds_with_assay_meta(names, assay_meta):
    """sample_names と assay_metadata を持つ最小の DatasetState を作る。"""
    ds = DatasetState()
    ds.feature_matrix = np.arange(3 * len(names), dtype=float).reshape(3, len(names)) + 1.0
    ds.sample_names = list(names)
    ds.feature_ids = ["1", "2", "3"]
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(len(names))]
    ds.assay_metadata = assay_meta
    return ds


# ---------- build_dataset_pp_inputs ----------

def test_build_dataset_pp_inputs_transposes():
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _make_ds(n_features=20, n_samples=8)
    matrix, sample_names, feature_names, roles, sample_meta = build_dataset_pp_inputs(ds)
    assert matrix.shape == (8, 20)          # (n_samples, n_features)
    assert sample_names == ds.sample_names
    assert feature_names == ds.feature_ids
    assert set(roles) == set(ds.sample_names)
    assert set(sample_meta) == set(ds.sample_names)


def test_build_dataset_pp_inputs_detects_blank_role():
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _make_ds(with_blank=True)
    _, _, _, roles, _ = build_dataset_pp_inputs(ds)
    assert roles["blank_1"] == "blank"


def test_build_dataset_pp_inputs_uses_mztab_injection_order():
    """assay 順序を全件 None にすると QC ドリフト補正が常に skipped になる。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["s1", "s2"], {
        "assay[1]": {"run_order": 7, "batch": "1"},
        "assay[2]": {"run_order": 3, "batch": "1"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["s1"]["run_order"] == 7
    assert sample_meta["s2"]["run_order"] == 3
    assert sample_meta["s1"]["run_order_source"] == "mztab_injection_sequence"


def test_build_dataset_pp_inputs_prefers_filename_date_when_batch_is_constant():
    """既定の定数バッチを採ると日付由来の交絡を見落とす。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["20220901_a", "20220902_b"], {
        "assay[1]": {"run_order": 1, "batch": "1"},
        "assay[2]": {"run_order": 2, "batch": "1"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["20220901_a"]["batch"] == "20220901"
    assert sample_meta["20220901_a"]["batch_source"] == "filename_date"


def test_build_dataset_pp_inputs_uses_mztab_batch_when_it_varies():
    """変化する assay バッチを捨てると真の実験バッチを使えない。"""
    from metabolomix.analysis.dataset_analysis import build_dataset_pp_inputs
    ds = _ds_with_assay_meta(["20220901_a", "20220901_b"], {
        "assay[1]": {"run_order": 1, "batch": "B1"},
        "assay[2]": {"run_order": 2, "batch": "B2"},
    })
    _, _, _, _, sample_meta = build_dataset_pp_inputs(ds)
    assert sample_meta["20220901_a"]["batch"] == "B1"
    assert sample_meta["20220901_b"]["batch"] == "B2"
    assert sample_meta["20220901_a"]["batch_source"] == "mztab_batch_label"


def test_build_dataset_pp_inputs_rejects_empty_matrix():
    from metabolomix.analysis.dataset_analysis import (
        PreconditionError, build_dataset_pp_inputs,
    )
    ds = DatasetState()
    with pytest.raises(PreconditionError):
        build_dataset_pp_inputs(ds)
