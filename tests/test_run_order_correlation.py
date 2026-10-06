"""注入順と主成分の相関 — QC が無いバッチでドリフトを検出する唯一の手段。

QC 試料が無いと `qc_drift_correct` も QC-RSD フィルタも動かない（実データの
POS 60 サンプルがまさにそれで、60/60 が role=sample だった）。しかし
mzTab-M の `assay[N]-custom[...]` には injection sequence label が入っており、
主成分と注入順の相関を見れば、補正はできなくても**ドリフトの存在は言える**。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis.preprocessing import run_order_correlation


def _names(n):
    return [f"s{i}" for i in range(n)]


def test_perfect_positive_drift_gives_r_one():
    names = _names(6)
    order = {n: i + 1 for i, n in enumerate(names)}
    components = np.array([[float(i)] for i in range(6)])
    assert run_order_correlation(components, names, order) == [1.0]


def test_perfect_negative_drift_gives_r_minus_one():
    names = _names(6)
    order = {n: i + 1 for i, n in enumerate(names)}
    components = np.array([[float(-i)] for i in range(6)])
    assert run_order_correlation(components, names, order) == [-1.0]


def test_reports_one_value_per_component():
    names = _names(5)
    order = {n: i + 1 for i, n in enumerate(names)}
    components = np.column_stack([np.arange(5.0), np.zeros(5), np.arange(5.0)[::-1]])
    result = run_order_correlation(components, names, order)
    assert len(result) == 3
    assert result[0] == 1.0
    assert result[2] == -1.0


def test_missing_run_order_gives_none():
    """注入順を持たない mzTab-M では相関を主張してはいけない。"""
    names = _names(5)
    components = np.column_stack([np.arange(5.0), np.arange(5.0)])
    assert run_order_correlation(components, names, {n: None for n in names}) == [None, None]


def test_too_few_ordered_samples_gives_none():
    """2 点なら必ず |r|=1 になる。相関を語れる下限を 3 点に置く。"""
    names = _names(5)
    order = {"s0": 1, "s1": 2, "s2": None, "s3": None, "s4": None}
    components = np.array([[float(i)] for i in range(5)])
    assert run_order_correlation(components, names, order) == [None]


def test_constant_component_gives_none_not_nan():
    """分散ゼロの主成分は相関が定義できない。nan を返すと JSON が壊れる。"""
    names = _names(5)
    order = {n: i + 1 for i, n in enumerate(names)}
    components = np.zeros((5, 1))
    assert run_order_correlation(components, names, order) == [None]


def test_uses_only_samples_that_have_run_order():
    """注入順が欠けたサンプルは相関計算から外す（0 で埋めない）。"""
    names = _names(6)
    order = {"s0": 1, "s1": 2, "s2": 3, "s3": 4, "s4": None, "s5": None}
    # 注入順を持つ 4 点は完全な直線、持たない 2 点は外れ値。
    components = np.array([[0.0], [1.0], [2.0], [3.0], [999.0], [-999.0]])
    assert run_order_correlation(components, names, order) == [1.0]


def test_result_is_rounded_for_payload():
    names = _names(4)
    order = {n: i + 1 for i, n in enumerate(names)}
    components = np.array([[0.0], [1.0], [1.9], [3.3]])
    (r,) = run_order_correlation(components, names, order)
    assert isinstance(r, float)
    assert r == round(r, 3)
