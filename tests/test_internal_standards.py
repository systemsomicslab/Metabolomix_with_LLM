"""内部標準比（spec §8.1・§8）。

比は同一注入の `target_height / standard_height` で、単位は
`internal_standard_ratio`。絶対濃度ではない。ここで縛るのは、**分母が信用
できないセルを「小さい値」に化けさせない**こと:

- 分母が0・非有限・負、または profile が要求する検出条件を満たさないセルは欠損。
  微小定数を足して割らない（足せば比は分母の大きさを反映しない巨大値になる）。
- その欠損は**後段の補完も禁止**する。補完で埋めると「分母が無かった」という
  事実が消え、平均・分散に他試料の値が混ざる。
- 元の行列は破壊しない。生値（未補正 height）は QC と再計算の基準として残る。
"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.analysis.internal_standards import (
    apply_internal_standards,
    ratio,
    validate_pairs,
)
from metabolomix.core.atomic_io import DomainError


# ---------- brief記載のRED ----------


def test_different_target_ids_cannot_assign_two_standards_to_one_feature():
    pairs = [
        {"target_id": "a", "target_feature_id": "f", "standard_feature_id": "s1"},
        {"target_id": "b", "target_feature_id": "f", "standard_feature_id": "s2"},
    ]
    with pytest.raises(DomainError) as caught:
        apply_internal_standards([[12., 2., 3.]], ["f", "s1", "s2"], pairs, None)
    assert caught.value.code == "INTERNAL_STANDARD_MAP_INVALID"

def test_invalid_denominator_is_locked_missing():
    values, locked = ratio(np.array([10., 10.]), np.array([2., 0.]), None)
    assert values[0] == 5 and np.isnan(values[1])
    assert locked.tolist() == [False, True]


# ---------- 分母の条件 ----------

@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_every_invalid_denominator_locks_the_cell(bad):
    values, locked = ratio(np.array([10., 10.]), np.array([2., bad]), None)
    assert np.isnan(values[1])
    assert locked.tolist() == [False, True]


def test_undetected_standard_locks_even_with_a_usable_number():
    """gap-fill された分母は補間値。数値が入っていても比の根拠にしない。"""
    values, locked = ratio(np.array([10., 10.]), np.array([2., 2.]),
                           np.array([True, False]))
    assert values[0] == 5 and np.isnan(values[1])
    assert locked.tolist() == [False, True]


def test_missing_target_is_not_locked():
    """分子の欠損は通常の欠測。分母不正とは別で、後段の補完対象になりうる。"""
    values, locked = ratio(np.array([10., np.nan]), np.array([2., 2.]), None)
    assert np.isnan(values[1])
    assert locked.tolist() == [False, False]


def test_inputs_are_not_modified():
    target = np.array([10., 10.])
    standard = np.array([2., 0.])
    ratio(target, standard, None)
    assert target.tolist() == [10., 10.]
    assert standard.tolist() == [2., 0.]


def test_shape_mismatch_is_rejected():
    with pytest.raises(ValueError):
        ratio(np.array([1., 2., 3.]), np.array([1., 2.]), None)


# ---------- mapの検証 ----------

def _pairs(*items) -> list[dict]:
    return [{"target_id": t, "target_feature_id": tf,
             "standard_target_id": s, "standard_feature_id": sf}
            for t, tf, s, sf in items]


def test_one_target_with_two_standards_is_invalid():
    pairs = _pairs(("gaba", "1", "d6", "2"), ("gaba", "1", "d10", "3"))
    with pytest.raises(DomainError) as caught:
        validate_pairs(pairs, ["1", "2", "3"])
    assert caught.value.code == "INTERNAL_STANDARD_MAP_INVALID"


def test_self_reference_is_invalid():
    pairs = _pairs(("gaba", "1", "gaba", "1"))
    with pytest.raises(DomainError) as caught:
        validate_pairs(pairs, ["1"])
    assert caught.value.code == "INTERNAL_STANDARD_MAP_INVALID"


def test_cycle_is_invalid():
    pairs = _pairs(("a", "1", "b", "2"), ("b", "2", "a", "1"))
    with pytest.raises(DomainError) as caught:
        validate_pairs(pairs, ["1", "2"])
    assert caught.value.code == "INTERNAL_STANDARD_MAP_INVALID"


def test_unknown_feature_id_is_invalid():
    pairs = _pairs(("gaba", "1", "d6", "99"))
    with pytest.raises(DomainError) as caught:
        validate_pairs(pairs, ["1", "2"])
    assert caught.value.code == "INTERNAL_STANDARD_MAP_INVALID"


def test_valid_pairs_pass():
    validate_pairs(_pairs(("gaba", "1", "d6", "2")), ["1", "2", "3"])


# ---------- 行列への適用 ----------

def _matrix() -> np.ndarray:
    """assay 3 × feature 3。列0=対象、列1=内部標準、列2=未対応。"""
    return np.array([[10., 2., 100.],
                     [20., 4., 200.],
                     [30., 0., 300.]])


def test_units_distinguish_corrected_and_uncorrected_features():
    out = apply_internal_standards(
        _matrix(), ["1", "2", "3"], _pairs(("gaba", "1", "d6", "2")), None)
    assert out["units"] == ["internal_standard_ratio", "peak_height", "peak_height"]


def test_mapped_column_becomes_a_ratio_and_others_keep_raw_values():
    out = apply_internal_standards(
        _matrix(), ["1", "2", "3"], _pairs(("gaba", "1", "d6", "2")), None)
    assert out["values"][:, 0].tolist()[:2] == [5.0, 5.0]
    assert np.isnan(out["values"][2, 0])            # 分母0
    assert out["values"][:, 2].tolist() == [100., 200., 300.]


def test_the_standard_column_is_kept_as_support():
    """内部標準そのものは比にせず、support として元の値のまま残す。"""
    out = apply_internal_standards(
        _matrix(), ["1", "2", "3"], _pairs(("gaba", "1", "d6", "2")), None)
    assert out["values"][:, 1].tolist() == [2., 4., 0.]
    assert out["support_feature_ids"] == ["2"]


def test_locked_cells_record_why_they_are_missing():
    out = apply_internal_standards(
        _matrix(), ["1", "2", "3"], _pairs(("gaba", "1", "d6", "2")), None)
    assert out["locked_mask"][2, 0]
    assert not out["locked_mask"][0, 0]
    reason = out["missing_reasons"]["1"]["2"]
    assert reason["code"] == "standard_denominator_invalid"
    assert reason["standard_feature_id"] == "2"


def test_detection_requirement_uses_the_assay_mask():
    detected = np.array([[True, True, True],
                         [True, False, True],
                         [True, True, True]])
    out = apply_internal_standards(
        _matrix(), ["1", "2", "3"], _pairs(("gaba", "1", "d6", "2")), detected)
    assert np.isnan(out["values"][1, 0])
    assert out["locked_mask"][1, 0]


def test_the_source_matrix_is_never_modified():
    matrix = _matrix()
    before = matrix.copy()
    apply_internal_standards(matrix, ["1", "2", "3"],
                             _pairs(("gaba", "1", "d6", "2")), None)
    assert np.array_equal(matrix, before)


def test_no_pairs_leaves_every_column_uncorrected():
    out = apply_internal_standards(_matrix(), ["1", "2", "3"], [], None)
    assert out["units"] == ["peak_height"] * 3
    assert np.array_equal(out["values"], _matrix())
    assert not out["locked_mask"].any()
