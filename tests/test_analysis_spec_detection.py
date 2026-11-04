"""統計の行の群ごとの検出の例数・Fisher の p・detection_driven（reanalysis v9 #4）。"""
from __future__ import annotations

import numpy as np
import pytest

from metabolomix.analysis import analysis_spec as A
from metabolomix.analysis.spec_statistics import detection_driven, detection_fisher_p
from tests.test_alignment_qc import ADDUCT, ALL, FIRST_TWO, _make_ds
from tests.test_analysis_spec import _rows as _result_rows
from tests.test_analysis_spec import files  # noqa: F401  外部の行列の fixture
from tests.test_analysis_spec import _spec, MAIN


def test_fisher_presence_absence():
    assert detection_fisher_p(0, 6, 6, 6) == pytest.approx(2 / 924)


@pytest.mark.parametrize("direction,ref,test,expected", [
    ("test_higher", 0, 2, True),     # 高い側 2/6 < 3
    ("test_higher", 0, 6, True),     # 低い側の検出 0
    ("test_higher", 3, 6, False),
    ("test_higher", 1, 3, False),    # 3 は 3 未満ではない
    ("ref_higher", 2, 0, True),      # 高い側は ref
    (None, 0, 6, None),
    ("test_higher", None, 6, None),
])
def test_detection_driven_rule(direction, ref, test, expected):
    assert detection_driven(direction, ref, 6, test, 6) is expected


def _rows(ds, spec):
    result = A.run_spec(A.normalize_spec(spec), ds)
    return {r["feature_id"]: r for r in result["analyses"][0]["rows"]}


def _spec_for(merge_split):
    return {"comparison": {"reference": "A", "test": "B"},
            "features": {"all": True, "merge_split": merge_split},
            "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]}


def test_rows_carry_detection_counts_from_the_mask():
    ds = _make_ds([("1", 300.0, 2.0, ALL, "X", ADDUCT),
                   ("2", 400.0, 3.0, set(range(21, 42)), "Y", ADDUCT)], detected=True)
    rows = _rows(ds, _spec_for(False))
    assert (rows["1"]["detected_reference"], rows["1"]["detected_test"]) == (21, 21)
    assert rows["1"]["detection_driven"] is False
    assert (rows["2"]["detected_reference"], rows["2"]["detected_test"]) == (0, 21)
    assert rows["2"]["direction"] == "test_higher" and rows["2"]["detection_driven"] is True
    assert rows["2"]["detection_fisher_p"] < 1e-6
    assert type(rows["2"]["detected_test"]) is int


def test_no_detection_state_gives_nulls(files):  # noqa: F811
    result = A.run_spec(A.normalize_spec(_spec(files, [MAIN])))
    row = _result_rows(result)["F1"]
    assert [row[k] for k in ("detected_reference", "detected_test",
                             "detection_fisher_p", "detection_driven")] == [None] * 4


def test_dataset_without_detected_mask_gives_nulls():
    ds = _make_ds([("1", 300.0, 2.0, ALL, "X", ADDUCT)], detected=False)
    row = _rows(ds, _spec_for(False))["1"]
    assert row["detected_reference"] is None and row["detection_driven"] is None


def test_merged_column_counts_either_spot():
    # 1265 は最初の 2 注入と 10 番目を受け持たず、1266 は最初の 2 注入だけを受け持つ。
    # まとめた列は「どちらかで検出」なので、検出されないのは 10 番目の 1 試料だけ（A 群 20 / B 群 21）。
    ds = _make_ds([("1265", 650.6000, 20.00, ALL - FIRST_TWO - {10}, "CE 18:2", ADDUCT),
                   ("1266", 650.6050, 20.61, FIRST_TWO, None, ADDUCT),
                   ("2000", 900.5000, 12.00, ALL, "PC 34:1", ADDUCT)], detected=True)
    merged = _rows(ds, _spec_for(True))
    assert sorted(merged) == ["1265", "2000"] and merged["1265"]["merged_from"] == ["1265", "1266"]
    assert (merged["1265"]["detected_reference"], merged["1265"]["detected_test"]) == (20, 21)
    split = _rows(ds, _spec_for(False))
    assert split["1265"]["detected_reference"] == 18 and split["1266"]["detected_reference"] == 2
    assert np.all([r["detected_test"] is not None for r in merged.values()])
