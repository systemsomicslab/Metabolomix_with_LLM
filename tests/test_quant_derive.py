"""定量の派生データセット: quant_state・derive・ハンドルへの登録。"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.quant import derive as qd


def _ds() -> DatasetState:
    ds = DatasetState()
    ds.source_format = "msdial_folder"
    ds.feature_ids = ["f1", "f2"]
    ds.sample_names = ["s1", "s2", "s3"]
    ds.sample_assay_ids = ["a1", "a2", "a3"]
    ds.feature_matrix = np.arange(6, dtype=float).reshape(2, 3)
    ds.detected_mask = np.ones((2, 3), dtype=bool)
    ds.feature_metadata = {"f1": {"name": "x"}, "f2": {"name": "y"}}
    ds.pp_matrix = np.zeros((3, 2))
    ds.pp_sample_names = ["s1"]
    ds.pp_feature_names = ["f1"]
    ds.results = {"r": {"v": 1}}
    ds.analysis_matrices = {"m": {}}
    ds.last_pca = {"p": 1}
    ds.last_differential = {"d": 1}
    ds.preprocess_id = "pp"
    ds.preprocess_metadata_hash = "h"
    ds.preprocessing_recipe = {"a": 1}
    ds.handoff_report = {"restored": []}
    ds.roles = {"s1": "sample", "s2": "qc", "s3": "blank"}
    ds.sample_meta = {"s1": {"batch": 1}, "s2": {"batch": 2}, "s3": {"batch": 3}}
    ds.sample_ids = ["i1", "i2", "i3"]
    ds.sample_raw_paths = ["r1", "r2", "r3"]
    ds.msdial_samples = [{"name": "s1"}, {"name": "s2"}, {"name": "s3"}]
    ds.per_sample_mz = {"f1": {"s1": 1.0, "s3": 3.0}, "f2": {"s2": 2.0}}
    ds.feature_candidates = {"f1": [{"a": 1}], "f2": [{"a": 2}]}
    ds.feature_annotations = {"f1": {"l": 1}, "f2": {"l": 2}}
    ds.feature_tags = {"f1": ["t"], "f2": ["u"]}
    ds.sml_rows = [{"x": 1}]
    ds.sme_rows = [{"x": 1}]
    ds.inchikey_coverage = {"total_features": 2}
    ds.feature_qc = {"source": "arf", "n_cells": 6, "n_detected": 6, "gap_filled_rate": 0.0}
    return ds


def _state():
    return qd.quant_state(normalized=True, method="internal_standard", unit="nmol/mL",
                          record_id="rec1", steps=["normalize"], skipped=[])


def test_quant_state_shape_and_extra():
    s = qd.quant_state(normalized=False, method=None, unit=None, record_id="r",
                       steps=[], skipped=[], internal_standard_ids=["a"])
    assert s["schema"] == "lipidmix.quant-state.v1" == qd.QUANT_STATE_SCHEMA
    assert s["normalized"] is False and s["method"] is None
    assert s["internal_standard_ids"] == ["a"]


def test_derive_resets_results_and_ids():
    ds = _ds()
    new = qd.derive(ds, matrix=ds.feature_matrix * 2, detected=ds.detected_mask.copy(),
                    state=_state())
    assert new is not ds and new.dataset_id != ds.dataset_id
    assert new.dataset_id.startswith("ds_")
    assert new.results == {} and new.analysis_matrices == {}
    assert new.pp_matrix is None and new.pp_sample_names == [] and new.pp_feature_names == []
    assert new.last_pca is None and new.last_differential is None
    assert new.preprocess_id is None and new.preprocess_metadata_hash is None
    assert new.preprocessing_recipe == {} and new.handoff_report is None
    assert new.source_format == "msdial_folder"
    assert new.quant_state["normalized"] is True
    assert float(new.feature_matrix[1, 2]) == 10.0


def test_derive_substitutes_fields_and_validates_shape():
    ds = _ds()
    new = qd.derive(ds, matrix=np.ones((1, 2)), detected=np.ones((1, 2), dtype=bool),
                    feature_ids=["f2"], feature_metadata={"f2": {"name": "y"}},
                    sample_names=["s1", "s2"], state=_state())
    assert new.feature_ids == ["f2"] and new.sample_names == ["s1", "s2"]
    with pytest.raises(ValueError):
        qd.derive(ds, matrix=np.ones((1, 2)), detected=np.ones((1, 2), dtype=bool),
                  state=_state())
    with pytest.raises(ValueError):
        qd.derive(ds, matrix=np.ones((2, 3)), detected=np.ones((2, 2), dtype=bool),
                  state=_state())


def test_derive_keeps_original_untouched():
    ds = _ds()
    before = ds.feature_matrix.copy()
    new = qd.derive(ds, matrix=ds.feature_matrix * 0, detected=ds.detected_mask,
                    state=_state())
    new.feature_metadata["f1"]["name"] = "changed"
    assert np.array_equal(ds.feature_matrix, before)
    assert ds.feature_metadata["f1"]["name"] == "x"
    assert ds.results == {"r": {"v": 1}} and ds.quant_state is None


def test_register_new_handle_becomes_current():
    ds = _ds()
    session_state.session.dataset = ds
    src_id = session_state.session.current_handle_id
    src = handles.registry.get(src_id)
    src.arf2_path, src.dcl_path = "a.arf2", "a.dcl"
    new = qd.derive(ds, matrix=ds.feature_matrix, detected=ds.detected_mask, state=_state())
    hid = qd.register(src, new, step="normalize", parent_handle_ids=[src_id])
    assert hid != src_id and session_state.session.current_handle_id == hid
    assert session_state.session.dataset is new
    assert handles.registry.get(src_id).dataset is ds  # 元は残る
    h = handles.registry.get(hid)
    assert h.polarity_source == "quant"
    assert h.source == {"kind": "quant", "step": "normalize", "parents": [src_id]}
    assert h.arf2_path == "a.arf2" and h.dcl_path == "a.dcl"
    assert h.library is not src.library


def test_register_merge_leaves_paths_none():
    ds = _ds()
    session_state.session.dataset = ds
    src = handles.registry.get(session_state.session.current_handle_id)
    src.arf2_path = "a.arf2"
    new = qd.derive(ds, matrix=ds.feature_matrix, detected=ds.detected_mask, state=_state())
    hid = qd.register(src, new, step="merge", parent_handle_ids=["x", "y"])
    assert handles.registry.get(hid).arf2_path is None


def test_register_respects_limit(monkeypatch):
    monkeypatch.setenv(handles.MAX_HANDLES_ENV, "1")
    ds = _ds()
    session_state.session.dataset = ds
    src_id = session_state.session.current_handle_id
    src = handles.registry.get(src_id)
    new = qd.derive(ds, matrix=ds.feature_matrix, detected=ds.detected_mask, state=_state())
    with pytest.raises(DomainError) as exc:
        qd.register(src, new, step="normalize", parent_handle_ids=[src_id])
    assert exc.value.code == "DATASET_HANDLE_LIMIT"
    assert session_state.session.current_handle_id == src_id
    assert handles.registry.peek(src_id) is src


def test_derive_realigns_reordered_subset_by_name():
    ds = _ds()
    new = qd.derive(ds, matrix=np.ones((2, 2)), detected=np.ones((2, 2), dtype=bool),
                    sample_names=["s3", "s1"], state=_state())
    assert new.sample_assay_ids == ["a3", "a1"]
    assert new.sample_ids == ["i3", "i1"] and new.sample_raw_paths == ["r3", "r1"]
    assert new.roles == {"s3": "blank", "s1": "sample"}
    assert new.sample_meta == {"s3": {"batch": 3}, "s1": {"batch": 1}}
    assert [m["name"] for m in new.msdial_samples] == ["s3", "s1"]
    assert new.per_sample_mz == {"f1": {"s3": 3.0, "s1": 1.0}, "f2": {}}
    with pytest.raises(ValueError):
        qd.derive(ds, matrix=np.ones((2, 1)), detected=np.ones((2, 1), dtype=bool),
                  sample_names=["nope"], state=_state())


def test_derive_realigns_renamed_by_index():
    ds = _ds()
    new = qd.derive(ds, matrix=np.ones((2, 2)), detected=np.ones((2, 2), dtype=bool),
                    sample_names=["M1", "M2"], sample_index=[2, 0], state=_state())
    assert new.sample_names == ["M1", "M2"]
    assert new.sample_assay_ids == ["a3", "a1"] and new.sample_raw_paths == ["r3", "r1"]
    assert new.roles == {"M1": "blank", "M2": "sample"}
    assert new.per_sample_mz["f1"] == {"M1": 3.0, "M2": 1.0}
    with pytest.raises(ValueError):
        qd.derive(ds, matrix=np.ones((2, 2)), detected=np.ones((2, 2), dtype=bool),
                  sample_index=[0], state=_state())


def test_derive_narrowing_prunes_feature_side_data():
    ds = _ds()
    new = qd.derive(ds, matrix=np.ones((1, 3)), detected=np.ones((1, 3), dtype=bool),
                    feature_ids=["f2"], state=_state())
    assert set(new.feature_metadata) == {"f2"} and set(new.feature_candidates) == {"f2"}
    assert set(new.feature_annotations) == {"f2"} and set(new.feature_tags) == {"f2"}
    assert set(new.per_sample_mz) == {"f2"}
    assert new.sml_rows is None and new.sme_rows is None and new.inchikey_coverage == {}
    assert new.feature_qc["n_cells"] == 3 and new.feature_qc["source"] == "arf"
    # 元は触らない
    assert set(ds.feature_candidates) == {"f1", "f2"}


def test_derive_does_not_alias_inputs():
    ds = _ds()
    m = np.ones((2, 3))
    d = np.ones((2, 3), dtype=bool)
    st = _state()
    new = qd.derive(ds, matrix=m, detected=d, state=st)
    assert not np.shares_memory(new.feature_matrix, m)
    assert not np.shares_memory(new.detected_mask, d)
    assert new.quant_state is not st and new.quant_state == st
    assert new.feature_matrix.dtype == float and new.detected_mask.dtype == bool


def test_set_current_sets_handle_id():
    session_state.session.set_current("dh_x")
    assert session_state.session.current_handle_id == "dh_x"
    session_state.session.set_current(None)
