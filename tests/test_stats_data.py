"""stats.data: DatasetState からの入口と 3 状態。"""
from types import SimpleNamespace

import numpy as np

from metabolomix.analysis.stats import data as D


def _ds(**kw):
    base = dict(feature_matrix=np.array([[1.0, 2.0, 0.0], [4.0, 5.0, 6.0]]),
                detected_mask=np.array([[True, False, False], [True, True, True]]),
                feature_ids=["f1", "f2"], sample_names=["a", "b", "c"])
    base.update(kw)
    return SimpleNamespace(**base)


def test_from_dataset_reads_metadata_rows_and_subject_id():
    rows = [{"sample_id": "a", "role": "Sample", "group": "G1", "batch": "B1",
             "injection_order": "3", "subject_id": "S1", "biological_sample_id": "X", "age": "40"},
            {"sample_id": "b", "role": "qc", "biological_sample_id": "BIO2", "include": "false"},
            {"sample_id": "c", "group": "G2"}]
    sd = D.from_dataset(_ds(sample_metadata_rows=rows, sample_ids=["a", "b", "c"]))
    a, b, c = sd.samples
    assert (a.subject_id, a.role, a.group, a.batch, a.injection_order) == ("S1", "sample", "G1", "B1", 3.0)
    assert a.attributes == {"age": "40"}
    assert b.subject_id == "BIO2" and b.role == "qc" and b.include is False
    assert c.subject_id is None and c.role == "sample" and c.include is True
    assert sd.values.shape == (3, 2) and sd.detected.shape == (3, 2)


def test_from_dataset_falls_back_to_msdial_samples_then_roles():
    ms = [{"name": "a", "group": "G1", "role": "qc", "batch": "1", "injection_order": 7},
          {"name": "b", "group": "G2", "role": "sample", "batch": None, "injection_order": None}]
    sd = D.from_dataset(_ds(msdial_samples=ms))
    assert sd.samples[0].role == "qc" and sd.samples[0].injection_order == 7.0
    assert sd.samples[1].group == "G2"
    sd2 = D.from_dataset(_ds(roles={"a": "blank"}, sample_meta={"b": {"run_order": 5, "batch": "Z"}}))
    assert sd2.samples[0].role == "blank" and sd2.samples[1].injection_order == 5.0
    assert sd2.samples[1].batch == "Z" and sd2.samples[2].role == "sample"


def test_quant_state_port_absent_and_present():
    assert D.quant_state_of(_ds()) is None
    assert D.quant_state_of(_ds(quant_state=None)) is None
    q = D.quant_state_of(_ds(quant_state={"schema": "lipidmix.quant-state.v1", "normalized": True,
                                          "method": "internal_standard", "unit": "pmol", "steps": []}))
    assert q == D.QuantState(True, "internal_standard", "pmol")
    assert D.from_dataset(_ds(quant_state={"normalized": False})).quant.normalized is False


def test_cell_states_three_way_and_unknown():
    v = np.array([[1.0, 2.0, 0.0, -1.0, np.nan]])
    d = np.array([[True, False, True, True, True]])
    assert D.cell_states(v, d).tolist() == [[1, 2, 0, 0, 0]]
    assert D.cell_states(v, None).tolist() == [[3, 3, 0, 0, 0]]
    assert D.cell_states(v, d).dtype == np.int8


def test_replicate_name_warnings_flags_two_level_numbers_only():
    assert D.replicate_name_warnings(["HC_1", "HC_2", "RA_1"]) == []
    w = D.replicate_name_warnings(["M01_1", "M01_2", "M02_1", "S3_inj1", "S3_inj2"])
    assert {x["stem"] for x in w} == {"M01", "S3"}
    assert all(x["code"] == "possible_technical_replicates" for x in w)


def test_from_dataset_features_and_hash_fallback():
    ds = _ds(feature_metadata={"f1": {"ontology": "PC"}, "f2": {"ontology": ""}},
             feature_candidates={"f1": [{"rank": 1}], "f2": []})
    sd = D.from_dataset(ds)
    assert sd.feature_classes == {"f1": "PC"} and sd.annotated_features == frozenset({"f1"})
    assert sd.dataset_hash is None or isinstance(sd.dataset_hash, str)
    assert D.from_dataset(_ds(detected_mask=None)).detected is None


def test_from_dataset_does_not_mutate_dataset():
    ds = _ds()
    fm, dm = ds.feature_matrix.copy(), ds.detected_mask.copy()
    sd = D.from_dataset(ds)
    sd.values[:] = 99
    sd.detected[:] = False
    assert np.array_equal(ds.feature_matrix, fm) and np.array_equal(ds.detected_mask, dm)
