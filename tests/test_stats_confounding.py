from types import SimpleNamespace as U

import pytest

from metabolomix.analysis.stats.confounding import assert_not_confounded, confounding_report
from metabolomix.core.atomic_io import DomainError


def mk(spec):
    """spec: [(group, batch, order), ...]"""
    return [U(group=g, batch=b, injection_order=o) for g, b, o in spec]


def test_complete_confounding_blocks_unless_allowed():
    units = mk([("A", "b1", None)] * 3 + [("B", "b2", None)] * 3)
    rep = confounding_report(units, ["A", "B"])
    assert rep["batch_status"] == "complete"
    assert any(w["code"] == "batch_group_complete_confounding" for w in rep["warnings"])
    with pytest.raises(DomainError) as e:
        assert_not_confounded(rep, allow_confounded=False)
    assert e.value.code == "STATISTIC_CONFOUNDED"
    assert e.value.details["batch_table"] == {"A": {"b1": 3}, "B": {"b2": 3}}
    assert_not_confounded(rep, allow_confounded=True)


def test_partial_bias_single_group_batch_warns():
    units = mk([("A", "b1", None)] * 2 + [("A", "b2", None)] * 3 + [("B", "b2", None)] * 3)
    rep = confounding_report(units, ["A", "B"])
    assert rep["batch_status"] == "partial"
    assert any(w["code"] == "batch_group_partial_confounding" for w in rep["warnings"])
    assert_not_confounded(rep, allow_confounded=False)  # 完全でなければ止めない


def test_cramers_v_threshold():
    # [[8,2],[2,8]]: n=20, 期待値は全部 5, chi2 = 4*(3^2/5) = 7.2, V = sqrt(7.2/20) = 0.6
    biased = mk([("A", "b1", None)] * 8 + [("A", "b2", None)] * 2
                + [("B", "b1", None)] * 2 + [("B", "b2", None)] * 8)
    rep = confounding_report(biased, ["A", "B"])
    assert rep["cramers_v"] == pytest.approx(0.6)
    assert rep["batch_status"] == "partial"
    even = mk([(g, b, None) for g in "AB" for b in ("b1", "b2") for _ in range(5)])
    rep = confounding_report(even, ["A", "B"])
    assert rep["cramers_v"] == pytest.approx(0.0)
    assert rep["batch_status"] == "none"
    assert rep["batch_table"] == {"A": {"b1": 5, "b2": 5}, "B": {"b1": 5, "b2": 5}}


def test_missing_batch_is_not_evaluable_not_none():
    units = mk([("A", "b1", None), ("A", None, None), ("B", "b1", None), ("B", "b2", None)])
    rep = confounding_report(units, ["A", "B"])
    assert rep["batch_status"] == "not_evaluable"
    assert any(w["code"] == "batch_not_evaluable" for w in rep["warnings"])


def test_single_batch_and_uncounted_groups_ignored():
    units = mk([("A", "b1", None)] * 3 + [("B", "b1", None)] * 3 + [("C", "b2", None)] * 3)
    rep = confounding_report(units, ["A", "B"])
    assert rep["batch_status"] == "single_batch"
    assert "C" not in rep["batch_table"]


def test_injection_order_grouped_runs_warn():
    grouped = mk([("A", "b1", float(i)) for i in range(6)] + [("B", "b1", float(i)) for i in range(6, 12)])
    rep = confounding_report(grouped, ["A", "B"])
    io = rep["injection_order"]
    assert io["status"] == "biased" and io["kruskal_p"] < 0.05
    assert io["range_by_group"] == {"A": [0.0, 5.0], "B": [6.0, 11.0]}
    assert any(w["code"] == "injection_order_biased" for w in rep["warnings"])
    inter = mk([(("A", "B")[i % 2], "b1", float(i)) for i in range(12)])
    rep = confounding_report(inter, ["A", "B"])
    assert rep["injection_order"]["status"] == "ok"
    assert not any(w["code"].startswith("injection_order") for w in rep["warnings"])


def test_injection_order_not_evaluable_cases():
    none = confounding_report(mk([("A", "b1", None)] * 3 + [("B", "b1", None)] * 3), ["A", "B"])
    assert none["injection_order"]["status"] == "not_evaluable"
    assert not none["warnings"]
    partial = mk([("A", "b1", 1.0), ("A", "b1", 2.0), ("B", "b1", 3.0), ("B", "b1", None)])
    rep = confounding_report(partial, ["A", "B"])
    assert rep["injection_order"]["status"] == "not_evaluable"
    assert any(w["code"] == "injection_order_not_evaluable" for w in rep["warnings"])
    same = confounding_report(mk([("A", "b1", 1.0)] * 3 + [("B", "b1", 1.0)] * 3), ["A", "B"])
    assert same["injection_order"]["status"] == "not_evaluable"
    assert any(w["code"] == "injection_order_not_evaluable" for w in same["warnings"])
