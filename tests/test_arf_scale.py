""".arf の大きさの見積もりと `ARF_TOO_LARGE`、基のキャッシュの大きさの上限（B-2a Task 5）。

fixture はすべてこのテストが tmp に作る。
"""
from __future__ import annotations

import tracemalloc

import numpy as np
import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.mztab import evidence
from metabolomix.mztab.dataset_state import DatasetState
from tests.curation_fixtures import arf_row, write_arf
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.test_arf_streaming import _old_expectation


def _small_folder(tmp_path):
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(3)]
    spots = [{"name": "Unknown", "mz": 500.0 + k, "heights": [10.0, 20.0, 30.0]}
             for k in range(4)]
    return write_msdial_folder(tmp_path, spots=spots, samples=samples)


def test_too_large_refused_before_read(tmp_path, monkeypatch):
    paths = _small_folder(tmp_path)
    monkeypatch.setattr(handles, "available_memory_bytes", lambda: 1000)

    def no_read(*args, **kwargs):
        raise AssertionError("読む前に断るはずが、読んだ")

    monkeypatch.setattr(arf_reader, "iter_groups", no_read)
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ARF_TOO_LARGE"
    details = exc.value.details
    assert details["path_name"] == paths["arf"].name
    assert details["available_bytes"] == 1000
    assert details["estimated_bytes"] > 800
    assert "mztab" in details["remedy"]
    assert evidence.estimate_scan_bytes(paths["arf"]) == details["estimated_bytes"]


def test_estimate_scales_with_uncompressed_size(tmp_path):
    small = write_arf(tmp_path / "small.arf", [[arf_row(file_id=i, mz=300.0, rt=1.0,
                                                        height=1.0, gap_filled=False)
                                                for i in range(4)]] * 2)
    large = write_arf(tmp_path / "large.arf", [[arf_row(file_id=i, mz=300.0, rt=1.0,
                                                        height=1.0, gap_filled=False)
                                                for i in range(4)]] * 20)
    with open(large, "rb") as handle:
        raw = arf_reader.uncompressed_size(handle)
    assert evidence.estimate_scan_bytes(large) > raw
    assert evidence.estimate_scan_bytes(large) > 5 * evidence.estimate_scan_bytes(small)


def test_unknown_available_memory_does_not_refuse(tmp_path, monkeypatch):
    paths = _small_folder(tmp_path)
    monkeypatch.setattr(handles, "available_memory_bytes", lambda: None)
    ds = load_msdial_folder(select_alignment(tmp_path))
    assert ds.feature_matrix.shape == (4, 3)
    assert paths["arf"].is_file()


def test_available_memory_bytes_is_positive_or_none():
    value = handles.available_memory_bytes()
    assert value is None or value > 0


def test_memory_error_becomes_too_large(tmp_path, monkeypatch):
    _small_folder(tmp_path)
    real = arf_reader.iter_groups

    def failing(file_like):
        for index, group in enumerate(real(file_like)):
            if index == 2:
                raise arf_reader.ArfMemoryError(spots=2, rows=6)
            yield group

    monkeypatch.setattr(arf_reader, "iter_groups", failing)
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ARF_TOO_LARGE"
    assert exc.value.details["path_name"].endswith(".arf")
    assert exc.value.details["estimated_bytes"] > 0
    assert "remedy" in exc.value.details


def test_plain_memory_error_also_becomes_too_large(tmp_path, monkeypatch):
    _small_folder(tmp_path)

    def failing(file_like):
        raise MemoryError()
        yield  # pragma: no cover

    monkeypatch.setattr(arf_reader, "iter_groups", failing)
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ARF_TOO_LARGE"


def test_mztab_attach_reports_too_large_not_unreadable(tmp_path, monkeypatch):
    write_arf(tmp_path / "x_PeakProperties.arf", [[arf_row(file_id=i, mz=700.5, rt=1.0,
                                                           height=1.0, gap_filled=False)
                                                   for i in range(2)]])
    (tmp_path / "a.mzTab").write_text("x", encoding="utf-8")
    ds = DatasetState()
    ds.feature_ids = ["0"]
    ds.sample_names = ["s0", "s1"]
    ds.feature_matrix = np.ones((1, 2))
    ds.artifact_paths = {}
    monkeypatch.setattr(handles, "available_memory_bytes", lambda: 100)  # 小さな .arf でも超える

    evidence.attach_to_dataset(ds, tmp_path / "a.mzTab")

    assert ds.detected_mask is None
    assert ds.feature_qc["source"] is None
    assert ds.feature_qc["reason"] == "arf_too_large"
    detail = ds.feature_qc["detail"]
    assert detail["estimated_bytes"] > 80 and detail["available_bytes"] == 100
    assert "remedy" in detail
    assert any("arf_too_large" in w for w in ds.validation_result["warnings"])


def test_matching_arfs_lets_too_large_propagate(tmp_path, monkeypatch):
    path = write_arf(tmp_path / "a.arf", [[arf_row(file_id=0, mz=700.5, rt=1.0, height=1.0,
                                                   gap_filled=False)]])
    monkeypatch.setattr(handles, "available_memory_bytes", lambda: 100)
    with pytest.raises(DomainError) as exc:
        evidence.matching_arfs([path], feature_mz=[700.5], sample_names=["sample_0"])
    assert exc.value.code == "ARF_TOO_LARGE"


def test_matching_arfs_does_not_keep_rejected_scan_alive(tmp_path, monkeypatch):
    import weakref

    def spot(mz):  # 群は 3 行以上（1〜2 試料の群は Task 6 まで読まない）
        return [[arf_row(file_id=i, mz=mz, rt=1.0, height=1.0, gap_filled=False)
                 for i in range(3)]]

    good = write_arf(tmp_path / "a.arf", spot(700.5))
    bad = write_arf(tmp_path / "b.arf", spot(100.0))
    refs: list = []
    real = evidence.scan_arf

    def tracking(path, **kwargs):
        scan = real(path, **kwargs)
        refs.append(weakref.ref(scan))
        return scan

    monkeypatch.setattr(evidence, "scan_arf", tracking)
    seen_alive: list[bool] = []
    real_scan_evidence = evidence.scan_evidence

    def probe(scan, **kwargs):
        # 2 本目を接合する時点で、1 本目（却下）が生きていないこと。
        seen_alive.append(len(refs) > 1 and refs[0]() is not None)
        return real_scan_evidence(scan, **kwargs)

    monkeypatch.setattr(evidence, "scan_evidence", probe)
    passed, rejected, first = evidence.matching_arfs(
        [bad, good], feature_mz=[700.5], sample_names=["sample_0", "sample_1", "sample_2"])
    assert passed == [good] and "b.arf" in rejected and first is not None
    assert seen_alive == [False, False]


def test_400_samples_stream_within_bound(tmp_path):
    n_spots, n_samples = 300, 400
    groups = [[arf_row(file_id=i, mz=300.0 + k * 0.5, rt=1.0 + k * 0.01,
                       height=100.0 + 7 * i + k, gap_filled=(i + k) % 5 == 0)
               for i in range(n_samples)] for k in range(n_spots)]
    path = write_arf(tmp_path / "big.arf", groups)
    del groups
    matrix_bytes = n_spots * n_samples * 9
    with open(path, "rb") as handle:
        uncompressed = arf_reader.uncompressed_size(handle)
    # 1 ブロックの伸長（python-lz4 は出力の 2 倍を持つ）と圧縮ペイロードの写しは避けられない。
    # それに行列の 20 倍までなら許す（行ごとの dict を持つ実装は 150 MB 前後になり超える）。
    bound = 2 * uncompressed + 2 * path.stat().st_size + 20 * matrix_bytes

    tracemalloc.start()
    try:
        scan = evidence.scan_arf(path, measures=("height", "area"))
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    print(f"\ntracemalloc peak = {peak} bytes, bound = {bound} bytes, "
          f"matrix = {matrix_bytes} bytes, estimate = {evidence.estimate_scan_bytes(path)}")
    assert peak < bound
    assert evidence.estimate_scan_bytes(path) >= peak  # 見積もりはピークを下回らない

    old = _old_expectation(path)
    assert scan.file_ids == old["file_ids"]
    assert scan.spot_mz == old["spot_mz"]
    np.testing.assert_array_equal(scan.values["height"], old["height"])
    np.testing.assert_array_equal(scan.values["area"], old["area"])
    np.testing.assert_array_equal(scan.detected, old["detected"])
    np.testing.assert_array_equal(scan.per_sample_mz, old["mz"])
    np.testing.assert_array_equal(scan.msms_assigned, old["msms"])


def test_memory_error_while_building_dataset_becomes_too_large(tmp_path, monkeypatch):
    from metabolomix.dataset import msdial_adapter

    _small_folder(tmp_path)

    def failing(*args, **kwargs):
        raise MemoryError()

    monkeypatch.setattr(msdial_adapter, "load_spot_fields", failing)
    with pytest.raises(DomainError) as exc:
        load_msdial_folder(select_alignment(tmp_path))
    assert exc.value.code == "ARF_TOO_LARGE"
    details = exc.value.details
    assert details["path_name"].endswith(".arf") and details["estimated_bytes"] > 0
    assert "available_bytes" in details and "remedy" in details


# ---------- 基のキャッシュの大きさ ----------

def _state(n_features: int, n_samples: int = 10) -> DatasetState:
    ds = DatasetState()
    ds.feature_matrix = np.zeros((n_features, n_samples))
    return ds


def test_estimate_state_bytes_counts_arrays_and_dicts():
    ds = _state(100, 10)  # 8000
    ds.detected_mask = np.zeros((100, 10), dtype=bool)  # 1000
    ds.pp_matrix = np.zeros((10, 100))  # 8000
    ds.per_sample_mz = {"0": {"s0": 1.0, "s1": 2.0}, "1": {"s0": 1.0}}  # 3 要素
    ds.feature_metadata = {"0": {}, "1": {}}  # 2 件
    assert handles.estimate_state_bytes(ds) == 8000 + 1000 + 8000 + 3 * 100 + 2 * 1500
    assert handles.estimate_state_bytes(DatasetState()) == 0


def test_estimate_state_bytes_counts_stat_matrices():
    import numpy as np
    from types import SimpleNamespace
    prepared = SimpleNamespace(log2=np.zeros((4, 5)), qc_log2=np.zeros((2, 5)),
                               sample_values=np.zeros((6, 5)), detected_units=np.zeros((4, 5), dtype=bool),
                               sample_detected=None)
    ds = DatasetState()
    ds.stat_matrices = {"mx_1": {"prepared": prepared, "matrix_id": "mx_1"}}
    assert handles.estimate_state_bytes(ds) == 160 + 80 + 240 + 20


def test_base_cache_size_cap(monkeypatch):
    registry = handles.HandleRegistry()
    monkeypatch.setattr(handles, "MAX_BASE_BYTES", 20_000)
    a, b, c = _state(100), _state(100), _state(100)  # 8000 バイトずつ
    registry.put_base("a", a, "sig")
    registry.put_base("b", b, "sig")
    assert registry.cached_base("a") is not None  # a を直近にする
    registry.put_base("c", c, "sig")  # 24000 > 20000: 最も使われていない b が追い出される
    assert registry.cached_base("b") is None
    assert registry.cached_base("a") is not None and registry.cached_base("c") is not None

    huge = _state(1000)  # 80000 バイト: 1 つで上限を超える基は入れない
    registry.put_base("huge", huge, "sig")
    assert registry.cached_base("huge") is None
    assert registry.cached_base("a") is not None  # 既存の基は追い出さない
    # 同じキーの置き換えで古い大きい基を残さない。
    registry.put_base("a", huge, "sig")
    assert registry.cached_base("a") is None


def test_base_cache_count_cap_still_applies(monkeypatch):
    registry = handles.HandleRegistry()
    monkeypatch.setattr(handles, "MAX_BASES", 2)
    for key in "abc":
        registry.put_base(key, _state(1), "sig")
    assert registry.cached_base("a") is None
    assert registry.cached_base("b") is not None and registry.cached_base("c") is not None
