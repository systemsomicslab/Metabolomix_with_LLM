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
