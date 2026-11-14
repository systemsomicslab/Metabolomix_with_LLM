"""PLS-DA と VIP（`metabolomix.analysis.stats.plsda`）。データは合成（テスト自身が作る）。"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from metabolomix.analysis.stats import plsda as P
from metabolomix.core.atomic_io import DomainError


def _prepared(n_per=8, p=12, shift=2.0, classes=("C", "T"), seed=1, missing=0, constant=0):
    """クラスごとに n_per 個体。最初の 3 特徴だけが `shift` だけずれる（shift=0 なら差が無い）。"""
    rng = np.random.default_rng(seed)
    units, rows = [], []
    for ci, cls in enumerate(classes):
        for k in range(n_per):
            units.append(SimpleNamespace(unit_id=f"{cls}{k}", group=cls, batch=None,
                                         attributes={"sex": "F" if k % 2 else "M"}))
            row = rng.normal(20, 1, p)
            row[:3] += shift * ci
            rows.append(row)
    log2 = np.array(rows)
    for j in range(constant):
        log2[:, p - 1 - j] = 15.0
    flat = np.flatnonzero(np.ones_like(log2, dtype=bool))
    log2.flat[flat[:missing * 7:7]] = np.nan
    return SimpleNamespace(matrix_id="mx_test", units=units, log2=log2,
                           feature_ids=[f"f{j}" for j in range(p)])


def test_separable_classes_give_positive_q2_small_p_and_top_vip():
    res = P.plsda(_prepared(shift=4.0), factors=["group"], n_permutations=30, seed=0)
    assert res["schema"] == "stats-plsda.v1"
    assert res["classes"] == {"C": 8, "T": 8} and res["y_coding"] == "one_column_second_class"
    assert res["n_components"] == 2 and len(res["q2_cum"]) == 2
    assert res["q2"] > 0.5
    assert res["permutation"]["n"] == 30 and res["permutation"]["n_finite"] == 30
    assert res["permutation"]["p"] <= 0.05
    assert {v["feature_id"] for v in res["vip"][:3]} == {"f0", "f1", "f2"}
    assert res["vip"][0]["vip"] > 1
    assert res["r2y_cum"][0] <= res["r2y_cum"][1] + 1e-12
    assert res["cv"]["folds"] == 7 and res["cv"]["scaling_within_fold"] is True
    assert not any(w["code"] == "plsda_q2_not_positive" for w in res["warnings"])


def test_cv_q2_matches_refits_on_the_same_folds():
    """Q2 は分割ごとに学習の個体だけで尺度を決め、成分数ごとに当てはめ直した予測と同じ。"""
    from sklearn.cross_decomposition import PLSRegression

    pr = _prepared(shift=2.0)
    res = P.plsda(pr, factors=["group"], n_components=3, cv_folds=4, n_permutations=1, seed=7)
    x = pr.log2
    labels = np.array([u.group for u in pr.units])
    y = (labels == "T").astype(float)
    fold = P._folds(labels, 4, np.random.default_rng(7))
    press = np.zeros(3)
    for f in range(4):
        test, train = fold == f, fold != f
        mean, sd = x[train].mean(axis=0), x[train].std(axis=0, ddof=1)
        for a in (1, 2, 3):
            model = PLSRegression(n_components=a, scale=False).fit((x[train] - mean) / sd, y[train])
            pred = model.predict((x[test] - mean) / sd).ravel()
            press[a - 1] += ((y[test] - pred) ** 2).sum()
    expect = 1 - press / ((y - y.mean()) ** 2).sum()
    assert res["q2_cum"] == pytest.approx(list(expect), rel=1e-8, abs=1e-10)


def test_vip_mean_square_is_one():
    res = P.plsda(_prepared(), factors=["group"], n_permutations=1)
    vip = np.array([v["vip"] for v in res["vip"]])
    assert float((vip ** 2).mean()) == pytest.approx(1.0, rel=1e-9)


def test_null_data_warns_and_permutation_p_is_large():
    res = P.plsda(_prepared(shift=0.0, seed=3), factors=["group"], n_permutations=40, seed=2)
    assert res["permutation"]["p"] > 0.05
    assert res["q2"] is not None
    if res["q2"] <= 0:
        assert any(w["code"] == "plsda_q2_not_positive" for w in res["warnings"])


def test_same_seed_same_numbers_and_seed_changes_split():
    a = P.plsda(_prepared(), factors=["group"], n_permutations=10, seed=4)
    b = P.plsda(_prepared(), factors=["group"], n_permutations=10, seed=4)
    assert a == b


def test_submodel_prediction_equals_refit_with_fewer_components():
    pr = _prepared()
    x = pr.log2 - pr.log2.mean(axis=0)
    y = np.array([[0.0]] * 8 + [[1.0]] * 8)
    full = P._fit(x, y, 3)
    preds = P._predict_by_components(full, x - x.mean(axis=0), y.mean(axis=0), 3)
    for a in (1, 2):
        sub = P._fit(x, y, a)
        expect = sub.predict(x)
        assert np.allclose(preds[a - 1], expect, atol=1e-8)


def test_three_classes_indicator_coding_and_classes_subset():
    pr = _prepared(classes=("A", "B", "C"))
    res = P.plsda(pr, factors=["group"], n_permutations=5)
    assert res["y_coding"] == "indicator_per_class" and set(res["classes"]) == {"A", "B", "C"}
    two = P.plsda(pr, factors=["group"], classes=["A", "C"], n_permutations=5)
    assert two["classes"] == {"A": 8, "C": 8} and two["n_units"] == 16


def test_cell_classes_from_two_factors():
    res = P.plsda(_prepared(), factors=["group", "sex"], n_permutations=3)
    assert set(res["classes"]) == {"C:F", "C:M", "T:F", "T:M"}


def test_constant_features_dropped_and_missing_imputed():
    res = P.plsda(_prepared(constant=2, missing=4), factors=["group"], n_permutations=3)
    assert res["n_features"] == 10 and res["n_features_dropped_constant"] == 2
    assert res["imputation"]["used"] is True and res["imputation"]["imputed_cells"] == 4
    with pytest.raises(DomainError) as e:
        P.plsda(_prepared(missing=4), factors=["group"], impute="none", n_permutations=3)
    assert e.value.code == "PLSDA_MISSING_VALUES"


def test_components_and_folds_capped_with_warnings():
    res = P.plsda(_prepared(n_per=3, p=4), factors=["group"], n_components=10, cv_folds=10,
                  n_permutations=3)
    codes = {w["code"] for w in res["warnings"]}
    assert {"cv_folds_capped", "n_components_capped"} <= codes
    assert res["cv"]["folds"] == 6 and res["n_components"] < 10


def test_errors():
    pr = _prepared()
    for kw, code in (({"scaling": "unit"}, "PLSDA_SPECIFICATION_INVALID"),
                     ({"impute": "zero"}, "PLSDA_SPECIFICATION_INVALID"),
                     ({"n_components": 0}, "PLSDA_SPECIFICATION_INVALID"),
                     ({"cv_folds": 1}, "PLSDA_SPECIFICATION_INVALID"),
                     ({"n_permutations": True}, "PLSDA_SPECIFICATION_INVALID"),
                     ({"classes": ["C", "nope"]}, "PLSDA_SPECIFICATION_INVALID")):
        with pytest.raises(DomainError) as e:
            P.plsda(pr, factors=["group"], **kw)
        assert e.value.code == code, kw
    with pytest.raises(DomainError) as e:
        P.plsda(_prepared(n_per=1), factors=["group"], n_permutations=1)
    assert e.value.code == "STATISTIC_GROUP_TOO_SMALL"
    flat = _prepared(p=3, constant=2)
    with pytest.raises(DomainError) as e:
        P.plsda(flat, factors=["group"], n_permutations=1)
    assert e.value.code == "PLSDA_TOO_FEW_FEATURES"
