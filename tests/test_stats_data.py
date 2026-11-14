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
            {"sample_id": "b", "role": "qc", "biological_sample_id": "BIO2", "include": False},
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
    sd2 = D.from_dataset(_ds(roles={"a": "blank"}, sample_meta={"a": {"batch": "Y"},
                                                    "b": {"run_order": 5, "batch": "Z"}}))
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


def test_role_inferred_from_names_when_no_source_gives_one():
    ds = _ds(sample_names=["QC_01", "Blank_1", "S01"])
    assert [s.role for s in D.from_dataset(ds).samples] == ["qc", "blank", "sample"]


def test_detected_mask_shape_mismatch_raises():
    import pytest
    from metabolomix.core.atomic_io import DomainError
    with pytest.raises(DomainError) as e:
        D.from_dataset(_ds(detected_mask=np.ones((2, 2), dtype=bool)))
    assert e.value.code == "STUDY_DATA_INVALID"


def test_fallbacks_found_by_display_name_when_ids_differ():
    ds = _ds(sample_ids=["id1", "id2", "id3"], sample_names=["a", "b", "c"],
             sample_metadata_rows=[{}, {}, {}],
             roles={"a": "QC"}, sample_meta={"a": {"batch": "Y"},
                                  "b": {"role": "blank", "batch": "Z", "run_order": 4}},
             msdial_samples=[{"name": "a"}, {"name": "b"}, {"name": "c", "group": "G9"}])
    s = D.from_dataset(ds).samples
    assert s[0].sample_id == "id1" and s[0].role == "qc"
    assert s[1].role == "blank" and s[1].batch == "Z" and s[1].injection_order == 4.0
    assert s[2].group == "G9"


def test_include_follows_analysis_spec_truthiness():
    ds = _ds(sample_metadata_rows=[{"include": False}, {"include": "0"}, {}], sample_ids=["a", "b", "c"])
    assert [s.include for s in D.from_dataset(ds).samples] == [False, True, True]


def test_role_inferred_from_display_names_when_ids_are_opaque():
    ds = _ds(sample_ids=["s001", "s002", "s003"], sample_names=["QC_01", "Blank_1", "S01"],
             sample_metadata_rows=[{}, {}, {}])
    assert [s.role for s in D.from_dataset(ds).samples] == ["qc", "blank", "sample"]


def test_role_inference_uses_msdial_class_id():
    ds = _ds(sample_names=["x1", "x2", "x3"],
             msdial_samples=[{"name": "x1", "group": "QC"}, {"name": "x2", "group": "Blank"},
                             {"name": "x3", "group": "KO"}])
    assert [s.role for s in D.from_dataset(ds).samples] == ["qc", "blank", "sample"]


# --- 実際の試料の対応表の行（帳簿の欄は因子にしない） -----------------------------------

def _manifest_ds():
    """`sample_manifest` の実際の行（注入ごとに違う provenance・conflicts・qc_pool 付き）を持つ DatasetState。"""
    from tests.test_stats_results import _ds
    ds = _ds()
    for i, row in enumerate(ds.sample_metadata_rows):
        row["provenance"] = {"role": {"value": row["role"], "source": f"file_{i}", "confidence": "inferred"},
                             "injection_order": {"value": i, "source": "mztab", "confidence": "unverified"}}
        row["conflicts"] = [{"field": "batch", "value": f"B{i}"}] if i % 3 == 0 else []
        row["qc_pool"] = f"pool{i}" if row["role"] == "qc" else None
    return ds


def test_bookkeeping_keys_do_not_become_attributes():
    ds = _manifest_ds()
    study = D.from_dataset(ds)
    for sample in study.samples:
        assert not {"provenance", "conflicts", "qc_pool"} & set(sample.attributes)
        assert all(isinstance(v, str) for v in sample.attributes.values())


def test_non_scalar_row_values_are_not_attributes():
    ds = _manifest_ds()
    ds.sample_metadata_rows[4]["tags"] = ["a", "b"]
    ds.sample_metadata_rows[4]["nested"] = {"x": 1}
    ds.sample_metadata_rows[4]["sex"] = "F"
    sample = D.from_dataset(ds).samples[4]
    assert sample.attributes == {"sex": "F"}


def test_replicates_collapse_with_real_manifest_rows():
    """技術反復（同じ subject_id の 2 注入）は、provenance が注入ごとに違っても 1 個体にまとまる。"""
    from metabolomix.analysis.stats import preprocess as P
    ds = _manifest_ds()
    study = D.from_dataset(ds)
    off = {"drift": {"mode": "off"},
           "filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
           "normalize": {"method": "none"}}
    prepared = P.preprocess(study, off)
    assert len(prepared.units) == 8                    # 8 個体 × 2 注入（以前は 16 に割れた）
    assert all(u.n_injections == 2 for u in prepared.units)
    assert not any(w.get("code") == "replicates_split" for w in prepared.warnings)


# --- バッチの読み方（v1・解析フォルダ経路と同じ規則。定義は sample_manifest の 1 つ） ---

def _batches(**kw):
    return [x.batch for x in D.from_dataset(_ds(**kw)).samples]


def test_constant_declared_batch_is_ignored_and_filename_dates_are_only_candidates():
    """A7: 試料名の 8 桁はバッチに当てない（候補として `batch_candidates` に持つだけ）。"""
    names = ["smp_20260101_a", "smp_20260101_b", "smp_20260102_c"]
    ms = [{"name": n, "batch": "1"} for n in names]
    assert _batches(sample_names=names, msdial_samples=ms) == [None, None, None]
    meta = {n: {"batch": "1"} for n in names}
    study = D.from_dataset(_ds(sample_names=names, sample_meta=meta))
    assert [x.batch for x in study.samples] == [None, None, None]
    assert study.batch_candidates == {
        "smp_20260101_a": {"value": "20260101", "source": "filename_date"},
        "smp_20260101_b": {"value": "20260101", "source": "filename_date"},
        "smp_20260102_c": {"value": "20260102", "source": "filename_date"}}


def test_filename_date_batch_from_an_old_row_or_sample_meta_is_not_used():
    names = ["smp_20260101_a", "smp_20260102_b", "c"]
    rows = [{"sample_id": n, "batch": "20260101" if i == 0 else None,
             "provenance": {"batch": {"value": "20260101", "source": "filename_date",
                                      "confidence": "inferred"}} if i == 0 else {}}
            for i, n in enumerate(names)]
    assert _batches(sample_names=names, sample_ids=names, sample_metadata_rows=rows) == [None, None, None]
    meta = {n: {"batch": d, "batch_source": "filename_date"}
            for n, d in zip(names, ["20260101", "20260102", None])}
    assert _batches(sample_names=names, sample_meta=meta) == [None, None, None]


def test_informative_declared_batches_win_over_dates():
    names = ["smp_20260101_a", "smp_20260101_b", "smp_20260102_c"]
    ms = [{"name": n, "batch": b} for n, b in zip(names, ["X", "Y", "Y"])]
    assert _batches(sample_names=names, msdial_samples=ms) == ["X", "Y", "Y"]


def test_no_declared_batch_and_no_date_stays_unknown_and_row_batch_wins():
    assert _batches() == [None, None, None]
    names = ["smp_20260101_a", "b", "c"]
    rows = [{"sample_id": n, "batch": "B7" if n == "b" else None} for n in names]
    study = D.from_dataset(_ds(sample_names=names, sample_ids=names, sample_metadata_rows=rows))
    assert [x.batch for x in study.samples] == [None, "B7", None]
    assert study.batch_candidates == {"smp_20260101_a": {"value": "20260101", "source": "filename_date"}}


# --- 注入順の出所とファイルの順（A6） -------------------------------------------------

def test_injection_order_source_by_origin():
    names = ["a", "b", "c"]
    rows = [{"sample_id": "a", "injection_order": 3,
             "provenance": {"injection_order": {"source": "user_manifest"}}},
            {"sample_id": "b", "injection_order": 1,
             "provenance": {"injection_order": {"source": "msdial_mddata"}}},
            {"sample_id": "c", "injection_order": 2,
             "provenance": {"injection_order": {"source": "mztab"}}}]
    study = D.from_dataset(_ds(sample_names=names, sample_ids=names, sample_metadata_rows=rows))
    assert [x.injection_order_source for x in study.samples] == ["sheet", "mddata", "mzTab"]
    ms = [{"name": n, "injection_order": i} for i, n in enumerate(names, start=1)]
    assert [x.injection_order_source for x in D.from_dataset(_ds(msdial_samples=ms)).samples] == [
        "mddata"] * 3
    assert [x.injection_order_source for x in D.from_dataset(_ds()).samples] == [None] * 3


def test_msdial_order_matching_file_names_is_flagged():
    """MS-DIAL の既定（ファイル名の順に 1, 2, …）と一致する注入順には印を付ける。"""
    names = ["s_B", "s_A", "s_C"]
    file_order = [{"name": n, "injection_order": o} for n, o in zip(names, [2, 1, 3])]
    flagged = D.from_dataset(_ds(sample_names=names, msdial_samples=file_order)).samples
    assert all(x.injection_order_may_be_file_order for x in flagged)
    shuffled = [{"name": n, "injection_order": o} for n, o in zip(names, [1, 3, 2])]
    assert not any(x.injection_order_may_be_file_order
                   for x in D.from_dataset(_ds(sample_names=names, msdial_samples=shuffled)).samples)
    # 利用者のシートの値は（ファイルの順と一致しても）印を付けない
    rows = [{"sample_id": n, "injection_order": o, "provenance": {"injection_order": {"source": "user_manifest"}}}
            for n, o in zip(names, [2, 1, 3])]
    sheet = D.from_dataset(_ds(sample_names=names, sample_ids=names, sample_metadata_rows=rows)).samples
    assert not any(x.injection_order_may_be_file_order for x in sheet)


def test_file_order_uses_raw_file_names_when_known():
    """鍵は raw のファイル名（全試料で分かるとき）。試料名の順とは違ってよい。"""
    names = ["zeta", "alpha", "mid"]
    ms = [{"name": n, "injection_order": o, "file_path": f"C:/raw/{f}"}
          for n, o, f in zip(names, [1, 2, 3], ["r01.wiff", "r02.wiff", "r03.wiff"])]
    samples = D.from_dataset(_ds(sample_names=names, msdial_samples=ms)).samples
    assert all(x.injection_order_may_be_file_order for x in samples)
