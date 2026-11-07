from types import SimpleNamespace as U

import numpy as np
import pytest

from metabolomix.analysis.stats.design import build_design
from metabolomix.core.atomic_io import DomainError


def u(group, subject=None, batch=None, **attrs):
    return U(group=group, batch=batch, subject_id=subject, attributes=dict(attrs))


def test_two_groups_contrast_is_test_minus_reference():
    units = [u("HC"), u("HC"), u("RA"), u("RA")]
    d = build_design(units, factors=["group"],
                     comparisons=[{"id": "RA_vs_HC", "factor": "group",
                                   "reference": "HC", "test": "RA"}])
    assert d.cells == ["HC", "RA"] and d.columns == ["HC", "RA"]
    assert d.matrix.shape == (4, 2)
    assert list(d.contrasts["RA_vs_HC"]) == [-1.0, 1.0]
    assert d.contrast_info["RA_vs_HC"]["test_cells"] == ["RA"]


def _two_by_two():
    units = []
    for a in ("a1", "a2"):
        for b in ("b1", "b2"):
            units += [u(a, B=b), u(a, B=b)]
    return units


def test_two_by_two_simple_effect_and_interaction():
    d = build_design(_two_by_two(), factors=["group", "B"], comparisons=[
        {"id": "s", "factor": "group", "reference": "a1", "test": "a2", "at": {"B": "b2"}},
        {"id": "i", "kind": "interaction", "factors": ["group", "B"],
         "levels": {"group": ["a1", "a2"], "B": ["b1", "b2"]}}])
    ix = {c: i for i, c in enumerate(d.columns)}
    assert d.cells == ["a1:b1", "a1:b2", "a2:b1", "a2:b2"]
    s = d.contrasts["s"]
    assert s[ix["a2:b2"]] == 1 and s[ix["a1:b2"]] == -1 and s.sum() == 0
    i = d.contrasts["i"]
    assert i[ix["a2:b2"]] == 1 and i[ix["a1:b2"]] == -1
    assert i[ix["a2:b1"]] == -1 and i[ix["a1:b1"]] == 1


def test_simple_effect_requires_at_for_two_factors():
    with pytest.raises(DomainError) as e:
        build_design(_two_by_two(), factors=["group", "B"], comparisons=[
            {"id": "s", "factor": "group", "reference": "a1", "test": "a2"}])
    assert e.value.code == "DESIGN_INVALID"


def test_unknown_level_is_invalid():
    with pytest.raises(DomainError) as e:
        build_design([u("HC"), u("RA")], factors=["group"], comparisons=[
            {"id": "x", "factor": "group", "reference": "HC", "test": "ZZ"}])
    assert e.value.code == "DESIGN_INVALID"


def test_block_nested_in_group_is_rank_deficient():
    units = [u("HC", "s1"), u("HC", "s2"), u("RA", "s3"), u("RA", "s4")]
    with pytest.raises(DomainError) as e:
        build_design(units, factors=["group"], comparisons=[], block="subject_id")
    assert e.value.code == "DESIGN_RANK_DEFICIENT"
    assert "columns" in e.value.details


def test_paired_block_is_full_rank():
    units = [u(g, s) for s in ("s1", "s2", "s3") for g in ("HC", "RA")]
    d = build_design(units, factors=["group"], block="subject_id", comparisons=[
        {"id": "c", "factor": "group", "reference": "HC", "test": "RA"}])
    assert d.columns == ["HC", "RA", "subject_id:s2", "subject_id:s3"]
    assert np.linalg.matrix_rank(d.matrix) == 4


def test_categorical_covariate_dummies_and_reference():
    units = [u("HC", batch="b1"), u("HC", batch="b2"), u("RA", batch="b3"),
             u("RA", batch="b1"), u("HC", batch="b3"), u("RA", batch="b2")]
    d = build_design(units, factors=["group"], comparisons=[],
                     covariates=[{"name": "batch", "type": "categorical", "reference": "b1"}])
    assert d.columns == ["HC", "RA", "batch:b2", "batch:b3"]


def test_numeric_covariate_and_missing_value():
    units = [u("HC", age="30"), u("HC", age="40"), u("RA", age="35"), u("RA", age="50")]
    d = build_design(units, factors=["group"], comparisons=[],
                     covariates=[{"name": "age", "type": "numeric"}])
    assert d.columns[-1] == "age" and list(d.matrix[:, -1]) == [30, 40, 35, 50]
    units[0].attributes.pop("age")
    with pytest.raises(DomainError) as e:
        build_design(units, factors=["group"], comparisons=[],
                     covariates=[{"name": "age", "type": "numeric"}])
    assert e.value.code == "DESIGN_INVALID"


def test_omnibus_has_cells_minus_one_columns():
    units = [u("A"), u("A"), u("B"), u("B"), u("C"), u("C")]
    d = build_design(units, factors=["group"], comparisons=[])
    assert d.omnibus.shape == (3, 2)
    assert list(d.omnibus[:, 1]) == [-1, 0, 1]
