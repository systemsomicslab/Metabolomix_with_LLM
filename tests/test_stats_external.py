"""`stats.external.prepared_from_injections`（注入ごとの log2 → 個体の `PreparedMatrix`）。"""
import numpy as np
import pytest

from metabolomix.analysis.stats.compare import compare
from metabolomix.analysis.stats.external import prepared_from_injections

NAN = np.nan


def test_collapses_replicates_by_subject_and_group_with_nanmean():
    log2 = [[1.0, 4.0], [3.0, NAN], [5.0, NAN], [7.0, 8.0]]
    p = prepared_from_injections(log2, sample_ids=["a1", "a2", "b1", "b2"],
                                 groups=["A", "A", "B", "B"],
                                 subject_ids=["s1", "s1", "s2", None], feature_ids=["f1", "f2"])
    assert [u.unit_id for u in p.units] == ["s1@A", "s2@B", "b2"]
    assert [u.n_injections for u in p.units] == [2, 1, 1]
    assert p.units[0].sample_ids == ["a1", "a2"]
    assert p.log2[0].tolist() == [2.0, 4.0]          # f2 は有限な 1 注入だけの平均
    assert p.log2[1].tolist()[0] == 5.0 and np.isnan(p.log2[1][1])
    assert p.subject_source == "mixed"


def test_same_subject_in_two_groups_stays_two_units():
    p = prepared_from_injections([[1.0], [2.0]], sample_ids=["x", "y"], groups=["A", "B"],
                                 subject_ids=["s", "s"], feature_ids=["f"])
    assert [u.group for u in p.units] == ["A", "B"]
    assert p.subject_source == "subject_id"


def test_shape_mismatch_is_value_error():
    with pytest.raises(ValueError):
        prepared_from_injections([[1.0]], sample_ids=["x", "y"], groups=["A"], subject_ids=[None],
                                 feature_ids=["f"])


def test_result_runs_through_compare():
    rng = np.random.default_rng(0)
    log2 = np.vstack([rng.normal(0, .1, (3, 4)), rng.normal(1, .1, (3, 4))])
    p = prepared_from_injections(log2, sample_ids=[f"s{i}" for i in range(6)],
                                 groups=["A"] * 3 + ["B"] * 3, subject_ids=[None] * 6,
                                 feature_ids=list("wxyz"))
    out = compare(p, factors=["group"], comparisons=[
        {"id": "c", "factor": "group", "reference": "A", "test": "B"}])
    assert len(out["rows"]) == 4 and all(r["status"] == "tested" for r in out["rows"])
