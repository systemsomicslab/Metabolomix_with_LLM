"""結果 ID・来歴・派生結果の無効化（spec §6）。

入力が変わったのに古い PCA / 差次的結果がセッションに残っていると、図と TSV が
**別々の入力から出た数字**を並べる。しかもどちらも「直近の結果」を名乗るので、
読み手には見分けがつかない。ここで縛るのは 2 つ:

1. **どの入力から出た結果かを結果自身が持つ**（来歴の `dataset_hash` と `metadata_revision`）。
2. **入力が変わったら派生状態を消す**。残す（＝古い数字を返し続ける）くらいなら、
   無いことにして計算し直させるほうが安全。
"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.analysis.result_state import (
    array_fingerprint,
    assert_current,
    dataset_fingerprint,
    invalidate_results,
    is_current,
    metadata_fingerprints,
)
from metabolomix.core.atomic_io import DomainError
from tests.engine_chain import register_two_group
from tests.pipeline_fixtures import make_dataset


# ---------- fingerprint ----------

def test_array_fingerprint_is_stable_for_equal_arrays():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[1.0, 2.0], [3.0, 4.0]])
    assert array_fingerprint(a) == array_fingerprint(b)


def test_array_fingerprint_separates_nan_from_zero():
    """NaN を JSON 経由で扱うと null になり 0 と混ざる。バイト列で見る。"""
    a = np.array([[np.nan, 1.0]])
    b = np.array([[0.0, 1.0]])
    assert array_fingerprint(a) != array_fingerprint(b)


def test_array_fingerprint_is_equal_for_two_nans():
    """NaN != NaN だが、同じ配列は同じ指紋でなければ再計算の判定が壊れる。"""
    assert array_fingerprint(np.array([np.nan])) == array_fingerprint(np.array([np.nan]))


def test_array_fingerprint_separates_shape_and_orientation():
    """転置は別の行列。要素の並びだけを見ると (n, m) と (m, n) が同一になる。"""
    a = np.arange(6.0).reshape(2, 3)
    assert array_fingerprint(a) != array_fingerprint(a.T)
    assert array_fingerprint(a) != array_fingerprint(a.reshape(3, 2))


def test_array_fingerprint_separates_dtype():
    assert (array_fingerprint(np.array([1, 2], dtype=np.int64))
            != array_fingerprint(np.array([1.0, 2.0], dtype=np.float64)))


def test_dataset_fingerprint_changes_with_the_matrix():
    ds = make_dataset()
    before = dataset_fingerprint(ds)
    ds.feature_matrix = ds.feature_matrix + 1.0
    assert dataset_fingerprint(ds) != before


def test_dataset_fingerprint_changes_with_the_detection_mask():
    """検出マスクは min_detection_rate の入力。変われば前処理の入力が変わる。"""
    ds = make_dataset()
    before = dataset_fingerprint(ds)
    ds.detected_mask = np.ones(ds.feature_matrix.shape, dtype=bool)
    assert dataset_fingerprint(ds) != before


def test_dataset_fingerprint_ignores_result_slots():
    """派生結果は入力ではない。ここが変わるたびに前処理をやり直させない。"""
    ds = make_dataset()
    before = dataset_fingerprint(ds)
    ds.last_pca = {"scores": [[1.0]]}
    assert dataset_fingerprint(ds) == before


def test_metadata_fingerprints_are_per_field():
    rows = [{"sample_id": "S0", "group": "ctrl", "role": "sample", "batch": "1"},
            {"sample_id": "S1", "group": "treat", "role": "sample", "batch": "1"}]
    before = metadata_fingerprints(rows)
    rows[1]["group"] = "ctrl"
    after = metadata_fingerprints(rows)
    assert after["group"] != before["group"]
    assert after["batch"] == before["batch"]
    assert after["role"] == before["role"]


# ---------- 無効化 ----------

def _with_results(ds):
    ds.last_pca = {"scores": []}
    ds.last_differential = {"kind": "compare"}
    return ds


def test_preprocess_relevant_metadata_change_invalidates_everything():
    ds = _with_results(make_dataset())
    invalidate_results(ds, {"batch"})
    assert ds.last_pca is None
    assert ds.last_differential is None


def test_group_change_invalidates_only_the_comparison():
    """群の付け替えは PCA の入力ではない。無関係な結果まで捨てない。"""
    ds = _with_results(make_dataset())
    invalidate_results(ds, {"group"})
    assert ds.last_pca is not None
    assert ds.last_differential is None


def test_pca_settings_change_invalidates_only_pca():
    ds = _with_results(make_dataset())
    invalidate_results(ds, {"pca_settings"})
    assert ds.last_pca is None
    assert ds.last_differential is not None


def test_unrelated_change_invalidates_nothing():
    ds = _with_results(make_dataset())
    invalidate_results(ds, {"note"})
    assert ds.last_pca is not None
    assert ds.last_differential is not None


def test_preprocess_relevant_change_also_drops_the_registered_matrices():
    """行列を残したまま結果だけ消すと、次の PCA が古い行列で走る。"""
    ds = make_dataset()
    ds.stat_matrices = {"mx_0": {}}
    invalidate_results(ds, {"batch"})
    assert ds.stat_matrices == {}


def test_changing_only_the_group_membership_keeps_the_matrices():
    """群は前処理の入力ではない。付け替えのたびに行列を作り直させない。"""
    ds = _with_results(make_dataset())
    ds.stat_matrices = {"mx_0": {}}
    invalidate_results(ds, {"group"})
    assert ds.stat_matrices == {"mx_0": {}}
    assert ds.last_pca is not None


def test_detection_mask_change_invalidates_the_matrices():
    ds = _with_results(make_dataset())
    ds.stat_matrices = {"mx_0": {}}
    invalidate_results(ds, {"detection"})
    assert ds.stat_matrices == {}
    assert ds.last_pca is None


# ---------- 有効性の検査 ----------

def test_a_result_without_provenance_is_not_called_stale():
    """来歴が無いのは「古い」の証拠ではない。手で組んだ結果を一律に殺さない。"""
    ds = make_dataset()
    assert is_current(ds, {"scores": []}) is True
    assert_current(ds, {"scores": []})


def test_a_species_pca_from_another_dataset_is_not_current():
    ds = make_dataset()
    foreign = {"provenance": {"result_id": "res_x", "kind": "species_pca", "dataset_id": "ds_other"}}
    assert is_current(ds, foreign) is False
    with pytest.raises(DomainError) as exc:
        assert_current(ds, foreign)
    assert exc.value.code == "STALE_ANALYSIS_RESULT"


# ---------- 群の訂正で差次的結果が古くなる（レビュー指摘3） ----------

def _dataset_with_groups(first_half: str, second_half: str):
    """8 検体すべて role=sample・include=true で、前半後半を別群に割った ds。"""
    from metabolomix.analysis.sample_manifest import apply_metadata
    from tests.pipeline_fixtures import metadata_rows

    ds = make_dataset()
    rows = metadata_rows(ds, n_qc=0)
    for i, row in enumerate(rows):
        row["group"] = first_half if i < 4 else second_half
    apply_metadata(ds, rows)
    return ds, rows


def _regroup(ds, rows, groups: list[str]):
    from metabolomix.analysis.sample_manifest import apply_metadata

    updated = [dict(row) for row in rows]
    for row, group in zip(updated, groups):
        row["group"] = group
    apply_metadata(ds, updated)
    return updated


def test_a_fresh_comparison_is_current_and_records_the_dataset_fingerprint():
    ds, _rows = _dataset_with_groups("control", "treated")
    record = register_two_group(ds, "control", "treated")
    assert ds.results[record["result_id"]] is record
    assert record["provenance"]["dataset_hash"] == dataset_fingerprint(ds)
    assert_current(ds, record)                      # 例外が出ないこと


def test_a_differential_result_goes_stale_when_the_groups_are_corrected():
    """群だけを直しても行列の入力は変わらない。だからこそ、旧群の結果は自分で古くなる。"""
    ds, rows = _dataset_with_groups("control", "treated")
    stale = register_two_group(ds, "control", "treated")

    # 群だけの訂正（1 検体の割当を直す）。
    _regroup(ds, rows, ["control"] * 3 + ["treated"] * 5)

    assert is_current(ds, stale) is False
    with pytest.raises(DomainError) as excinfo:
        assert_current(ds, stale)
    assert excinfo.value.code == "STALE_ANALYSIS_RESULT"


def test_exporting_a_result_from_the_old_groups_is_refused(tmp_path):
    """`dataset_export_differential`が約束している拒否を、実際の書き出しで確かめる。"""
    from metabolomix.analysis.dataset_export import export_dataset_result

    ds, rows = _dataset_with_groups("control", "treated")
    stale = register_two_group(ds, "control", "treated")
    _regroup(ds, rows, ["control"] * 3 + ["treated"] * 5)

    out = tmp_path / "differential.tsv"
    with pytest.raises(DomainError) as excinfo:
        export_dataset_result(ds, stale, out)
    assert excinfo.value.code == "STALE_ANALYSIS_RESULT"
    assert not out.exists()


def test_a_freshly_recomputed_comparison_is_current_again():
    """訂正後に計算し直した結果は当然通る（拒否が広すぎないことの確認）。"""
    ds, rows = _dataset_with_groups("control", "treated")
    register_two_group(ds, "control", "treated")
    _regroup(ds, rows, ["control"] * 3 + ["treated"] * 5)

    fresh = register_two_group(ds, "control", "treated")
    assert is_current(ds, fresh) is True
