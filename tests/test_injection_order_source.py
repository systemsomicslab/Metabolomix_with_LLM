"""注入順の出所と「MS-DIAL の既定（ファイルの順）かもしれない」印（A6）、シートの空欄の .mddata での補い。

対象: `sample_manifest` の出所の部品と `_merge_row`、交絡の判定（`confounding_report`）、
質量精度の注入順の傾向（`mass_accuracy.injection_trend`）。fixture はテスト自身が作る。
"""
from __future__ import annotations

from types import SimpleNamespace

from metabolomix.analysis import mass_accuracy, sample_manifest
from metabolomix.analysis.stats.confounding import confounding_report


def test_matches_file_order_rules():
    assert sample_manifest.matches_file_order([2, 1, 3], ["b", "a", "c"])
    assert sample_manifest.matches_file_order([1, 2, 3], ["A1", "a2", "B3"])          # 大文字小文字を無視
    assert not sample_manifest.matches_file_order([1, 3, 2], ["a", "b", "c"])
    assert not sample_manifest.matches_file_order([1, 2], ["a", "b"])                  # 3 試料未満は判定しない
    assert sample_manifest.matches_file_order([1, None, 2, 3], ["a", "z", "b", "c"])   # 値の無い試料は数えない


def test_source_labels():
    label = sample_manifest.injection_order_source_label
    assert [label(s) for s in (None, "user_manifest", "msdial_mddata", "mztab",
                               "mztab_injection_sequence", "raw_acquisition_time")] == [
        "sheet", "sheet", "mddata", "mzTab", "mzTab", "raw_acquisition_time"]


def _msdial_ds(orders, *, batches=("1", "1", "1"), names=("s_c", "s_a", "s_b")):
    samples = [{"name": n, "role": "sample", "group": "G", "batch": b, "injection_order": o,
                "file_path": f"C:/raw/{n}.wiff"} for n, b, o in zip(names, batches, orders)]
    return SimpleNamespace(sample_names=list(names), msdial_samples=samples, sample_assay_ids=[],
                           assay_metadata={}, assay_sources={})


def test_file_order_warning_from_default_rows():
    ds = _msdial_ds([3, 1, 2])                      # s_a=1, s_b=2, s_c=3 → ファイル名の順
    warn = sample_manifest.injection_order_file_order_warning(ds)
    assert warn["code"] == "injection_order_may_be_file_order" and warn["n_samples"] == 3
    assert sample_manifest.injection_order_file_order_warning(_msdial_ds([1, 2, 3])) is None


def _sheet_row(name, *, batch=None, order=None):
    explicit = {"source": "user_manifest", "confidence": "confirmed"}
    empty = {"value": None, "source": None, "confidence": None}
    return {"sample_id": name, "source_file": f"C:/raw/{name}.wiff", "role": "sample", "group": "G",
            "batch": batch, "injection_order": order, "qc_pool": None, "include": True,
            "biological_sample_id": None,
            "provenance": {"sample_id": {"value": name, **explicit},
                           "source_file": {"value": f"C:/raw/{name}.wiff", **explicit},
                           "role": {"value": "sample", **explicit}, "group": {"value": "G", **explicit},
                           "batch": {"value": batch, **explicit} if batch else dict(empty),
                           "injection_order": {"value": order, **explicit} if order else dict(empty),
                           "qc_pool": dict(empty), "include": {"value": True, **explicit},
                           "biological_sample_id": dict(empty)},
            "conflicts": []}


def test_sheet_blanks_are_filled_from_mddata_with_source():
    """A6: シートの空欄の注入順・バッチは .mddata の宣言で埋め、値と出所を行に残す（統計と同じ値）。"""
    ds = _msdial_ds([3, 1, 2], batches=("X", "Y", "Y"))
    merged = sample_manifest._merge_row(ds, 0, "s_c", _sheet_row("s_c"))
    assert merged["injection_order"] == 3 and merged["batch"] == "X"
    assert merged["provenance"]["injection_order"] == {"value": 3, "source": "msdial_mddata",
                                                       "confidence": "declared"}
    assert merged["provenance"]["batch"]["source"] == "msdial_mddata"
    # シートの値は上書きせず、食い違う .mddata の値は conflicts に残す
    kept = sample_manifest._merge_row(ds, 0, "s_c", _sheet_row("s_c", batch="B9", order=7))
    assert kept["injection_order"] == 7 and kept["batch"] == "B9"
    assert {c["field"] for c in kept["conflicts"]} == {"batch", "injection_order"}
    # 定数の .mddata のバッチ（GUI の既定 1）は情報が無いので埋めない
    const = sample_manifest._merge_row(_msdial_ds([3, 1, 2]), 0, "s_c", _sheet_row("s_c"))
    assert const["batch"] is None and const["injection_order"] == 3


def test_default_row_batch_keeps_filename_date_only_as_candidate():
    ds = _msdial_ds([1, 2, 3], names=("smp_20260101_a", "smp_20260101_b", "smp_20260102_c"))
    rows = sample_manifest.resolve_metadata(ds, None)
    assert [r["batch"] for r in rows] == [None, None, None]
    assert rows[2]["provenance"]["batch"]["candidate"] == {"value": "20260102", "source": "filename_date",
                                                           "confidence": "inferred"}


def test_confounding_report_flags_file_order():
    units = [SimpleNamespace(group=g, batch="b1", injection_order=float(o),
                             injection_order_may_be_file_order=True)
             for g, o in (("A", 1), ("A", 2), ("A", 3), ("B", 4), ("B", 5), ("B", 6))]
    rep = confounding_report(units, ["A", "B"])
    assert rep["injection_order"]["may_be_file_order"] is True
    warn = next(w for w in rep["warnings"] if w["code"] == "injection_order_may_be_file_order")
    assert warn["n_units"] == 6
    plain = [SimpleNamespace(group=u.group, batch=u.batch, injection_order=u.injection_order) for u in units]
    rep = confounding_report(plain, ["A", "B"])
    assert "may_be_file_order" not in rep["injection_order"]
    assert not any(w["code"] == "injection_order_may_be_file_order" for w in rep["warnings"])


def test_mass_accuracy_rho_carries_the_file_order_warning():
    names = [f"s{i}" for i in range(6)]
    ds = _msdial_ds(list(range(1, 7)), batches=("1",) * 6, names=tuple(names))
    ds.sample_metadata_rows = sample_manifest.resolve_metadata(ds, None)
    ds.per_sample_mz = {"f1": {n: 100.0 + 1e-4 * i for i, n in enumerate(names)}}
    rows = [{"feature_id": "f1", "theoretical_mz": 100.0}]
    trend = mass_accuracy.injection_trend(ds, rows)
    assert trend["injection_order_source"] == "metadata" and trend["spearman_rho"] is not None
    assert trend["warnings"][0]["code"] == "injection_order_may_be_file_order"
    shuffled = _msdial_ds([2, 1, 4, 3, 6, 5], batches=("1",) * 6, names=tuple(names))
    shuffled.sample_metadata_rows = sample_manifest.resolve_metadata(shuffled, None)
    shuffled.per_sample_mz = ds.per_sample_mz
    assert "warnings" not in mass_accuracy.injection_trend(shuffled, rows)
