"""PLS-DA と VIP（`metabolomix.analysis.stats.plsda`）。データは合成（テスト自身が作る）。"""
from __future__ import annotations

import math
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


def test_half_min_q2_and_null_match_training_only_imputation():
    from sklearn.cross_decomposition import PLSRegression

    pr = _prepared(n_per=8, p=6, seed=19)
    x = pr.log2
    x[0, 0] = -10
    x[[1, 2, 5, 10, 12, 15], 0] = np.nan
    res = P.plsda(pr, factors=["group"], n_components=1, cv_folds=4, n_permutations=1, seed=4)
    labels = np.array([u.group for u in pr.units])
    groups = np.array(["u:" + u.unit_id for u in pr.units])
    rng = np.random.default_rng(4)
    fold = P._folds(labels, 4, rng, groups)
    plan = P._permutation_plan(labels, groups, ["C", "T"])
    for outcome, expected in ((labels, res["q2"]),
                              (P._permuted(labels, plan, rng), res["permutation"]["q2_null_mean"])):
        y = (outcome == "T").astype(float)
        pred = np.empty_like(y)
        for i in range(4):
            train = fold != i
            fill = np.nanmin(x[train], axis=0) - 1.0
            filled = np.where(np.isfinite(x), x, fill)
            model = PLSRegression(n_components=1).fit(filled[train], y[train])
            pred[~train] = model.predict(filled[~train]).ravel()
        assert expected == pytest.approx(1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def test_knn_neighbours_are_learned_only_from_training_rows():
    rng = np.random.default_rng(42)
    raw = rng.normal(size=(12, 18))
    train = np.arange(12) < 10
    raw[2, 0] = raw[10, 0] = np.nan
    expected, _ = P.impute_for_pca(raw[train], method="knn")
    mu, sd = P._scale_params(expected, "autoscale")
    fit, held = P._fold_inputs(raw, train, "knn", "autoscale")
    np.testing.assert_allclose(fit, (expected - mu) / sd)
    # 評価試料の欠測は学習時の近傍から、その試料自身が持つ値を選ぶ。
    distances = np.nanmean((raw[train] - raw[train, 0, None]) ** 2, axis=0)
    neighbours = np.argsort(distances[1:])[:10] + 1
    assert held[0, 0] == pytest.approx((raw[10, neighbours].mean() - mu[0]) / sd[0])
    raw[11] *= 1000
    fit2, held2 = P._fold_inputs(raw, train, "knn", "autoscale")
    np.testing.assert_allclose(fit, fit2)
    np.testing.assert_allclose(held[0], held2[0])


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


def test_plsda_summary_in_analysis_run_has_no_rows(tmp_path):
    """dataset_analysis_run（spec v2）の要約: plsda は行の表を持たないので `rows` を付けない（pca と同じ）。"""
    import json

    from metabolomix.core import session_state
    from metabolomix.tools.reanalysis_tools import dataset_analysis_run
    from tests.test_stats_results import OFF, T_VS_C, _ds

    state = _ds()
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    spec = {"schema": "analysis-spec.v2", "preprocess": OFF, "comparisons": [T_VS_C],
            "analyses": [{"analysis_id": "cmp", "kind": "compare"},
                         {"analysis_id": "pl", "kind": "plsda", "n_permutations": 5}]}
    out = json.loads(dataset_analysis_run(spec))
    cmp_, pl = out["analyses"]
    assert cmp_["kind"] == "compare" and cmp_["rows"].startswith("contrast_id\t")
    assert pl["kind"] == "plsda" and "rows" not in pl
    assert {"q2", "permutation", "vip_top"} <= pl.keys()


# --- 独立な単位（同じ個体の複数の単位）で分割・並べ替えをする -------------------------------

def _nested(seed=100, n_subj=6, per=3, p=20, subject_sd=2.0, subject_ids=True):
    """2 群 × n_subj 個体 × 個体ごとに per 単位（時点の属性で分かれた単位）。個体の差だけが大きく、群の差は無い。"""
    rng = np.random.default_rng(seed)
    units, rows = [], []
    for cls in ("A", "B"):
        for s in range(n_subj):
            base = rng.normal(0.0, subject_sd, p)
            for t in range(per):
                sid = f"{cls}{s}"
                units.append(SimpleNamespace(unit_id=f"{sid}|t{t}", group=cls, batch=None,
                                             subject_id=sid if subject_ids else None,
                                             attributes={"time": f"t{t}"}))
                rows.append(20.0 + base + rng.normal(0.0, 0.3, p))
    return SimpleNamespace(matrix_id="mx_nested", units=units, log2=np.array(rows),
                           feature_ids=[f"f{j}" for j in range(p)])


def _paired(n_subj=8, effect=0.0, seed=0, p=20):
    """個体ごとに pre / post の 2 単位（対応のある計画）。"""
    rng = np.random.default_rng(seed)
    units, rows = [], []
    for s in range(n_subj):
        base = rng.normal(0.0, 2.0, p)
        for cls in ("pre", "post"):
            units.append(SimpleNamespace(unit_id=f"S{s}|{cls}", group=cls, batch=None, subject_id=f"S{s}",
                                         attributes={"pair": f"S{s}"}))
            row = 20.0 + base + rng.normal(0.0, 0.3, p)
            if cls == "post":
                row[:3] += effect
            rows.append(row)
    return SimpleNamespace(matrix_id="mx_paired", units=units, log2=np.array(rows),
                           feature_ids=[f"f{j}" for j in range(p)])


def test_units_of_one_subject_stay_in_one_fold_and_classes_are_spread():
    pr = _nested()
    labels = np.array([u.group for u in pr.units])
    groups = np.array([u.subject_id for u in pr.units])
    fold = P._folds(labels, 4, np.random.default_rng(0), groups)
    for sid in set(groups.tolist()):
        assert len(set(fold[groups == sid].tolist())) == 1          # 個体の単位は同じ分割
    for cls in ("A", "B"):
        per_fold = np.bincount([fold[groups == sid][0] for sid in sorted(set(groups[labels == cls].tolist()))],
                               minlength=4)
        assert per_fold.max() - per_fold.min() <= 1                 # クラスの個体を分割へ均等に配る


def test_permutation_moves_whole_subjects_or_stays_within_pairs():
    pr = _nested()
    labels = np.array([u.group for u in pr.units])
    groups = np.array([u.subject_id for u in pr.units])
    plan = P._permutation_plan(labels, groups, ["A", "B"])
    assert plan["scheme"] == "between_groups"
    assert plan["n_labelings"] == math.comb(12, 6) and plan["min_attainable_p"] == 2 / math.comb(12, 6)
    rng = np.random.default_rng(1)
    for _ in range(20):
        out = P._permuted(labels, plan, rng)
        for sid in set(groups.tolist()):
            assert len(set(out[groups == sid].tolist())) == 1      # 個体の単位は同じ札のまま
        assert {c: int((out == c).sum()) for c in ("A", "B")} == {"A": 18, "B": 18}
    pr = _paired()
    labels = np.array([u.group for u in pr.units])
    groups = np.array([u.subject_id for u in pr.units])
    plan = P._permutation_plan(labels, groups, ["post", "pre"])
    assert plan["scheme"] == "within_groups"
    assert plan["n_labelings"] == 2 ** 8 and plan["min_attainable_p"] == 2 / 2 ** 8
    out = P._permuted(labels, plan, np.random.default_rng(2))
    for sid in set(groups.tolist()):
        assert sorted(out[groups == sid].tolist()) == ["post", "pre"]  # 個体の中で入れ替える


def test_null_data_with_several_units_per_subject_is_not_validated():
    """群差の無いデータで、同じ個体の単位が学習と予測にまたがると Q2 と p が個体の違いで良く見えた。"""
    for seed in (100, 101, 102):
        res = P.plsda(_nested(seed=seed), factors=["group"], n_permutations=60, seed=seed)
        assert res["cv"]["method"] == "stratified_group_kfold"
        assert res["cv"]["group_by"] == "subject_id" and res["cv"]["n_groups"] == 12
        assert res["cv"]["class_groups"] == {"A": 6, "B": 6} and res["n_units"] == 36
        assert res["permutation"]["scheme"] == "between_groups"
        assert res["permutation"]["p"] > 0.1 and res["q2"] < 0.3


def test_real_difference_between_subjects_is_still_found():
    pr = _nested(seed=7, subject_sd=0.5)
    pr.log2[:18, :4] += 3.0                                         # A の個体だけ 4 特徴が高い
    res = P.plsda(pr, factors=["group"], n_permutations=60, seed=0)
    assert res["q2"] > 0.5 and res["permutation"]["p"] <= 0.05


def test_paired_design_permutes_within_subjects_and_block_names_the_group():
    res = P.plsda(_paired(effect=2.0), factors=["group"], n_permutations=60, seed=0)
    assert res["cv"]["group_by"] == "subject_id" and res["cv"]["n_groups"] == 8
    assert res["permutation"]["scheme"] == "within_groups"
    assert res["q2"] > 0.3 and res["permutation"]["p"] <= 0.05
    blocked = P.plsda(_paired(effect=2.0), factors=["group"], n_permutations=5, seed=0, block="pair")
    assert blocked["cv"]["group_by"] == "pair" and blocked["cv"]["n_groups"] == 8


def test_few_independent_units_warn_that_the_permutation_p_cannot_be_small():
    res = P.plsda(_paired(n_subj=3, effect=3.0), factors=["group"], n_permutations=20, seed=0)
    w = next(w for w in res["warnings"] if w["code"] == "permutation_few_labelings")
    assert w["n_labelings"] == 8 and w["min_attainable_p"] == 0.25
    assert res["permutation"]["min_attainable_p"] == 0.25 and res["permutation"]["p"] >= 0.2
    # 3 個体 × 2 クラス（入れ子）: C(6, 3) = 20 通り、名前の入れ替えで同じ Q2 なので p は 0.1 を下回れない
    res = P.plsda(_nested(n_subj=3, per=2), factors=["group"], n_permutations=40, seed=0)
    assert res["permutation"]["min_attainable_p"] == pytest.approx(0.1)
    assert any(w["code"] == "permutation_few_labelings" for w in res["warnings"])
    assert res["permutation"]["p"] > 0.05


def test_fold_cap_counts_subjects_and_one_subject_class_is_refused():
    res = P.plsda(_nested(n_subj=3), factors=["group"], cv_folds=10, n_permutations=3)
    assert res["cv"]["folds"] == 6 and res["n_units"] == 18
    capped = next(w for w in res["warnings"] if w["code"] == "cv_folds_capped")
    assert "subject_id" in capped["message"]
    with pytest.raises(DomainError) as e:
        P.plsda(_nested(n_subj=1, per=3), factors=["group"], n_permutations=3)
    assert e.value.code == "STATISTIC_GROUP_TOO_SMALL"
    assert e.value.details["class_groups"] == {"A": 1, "B": 1}


def test_units_without_subject_are_each_a_group_and_replicate_names_warn():
    pr = _nested(subject_ids=False)
    pr.warnings = [{"code": "possible_technical_replicates", "stem": "A0"}]
    res = P.plsda(pr, factors=["group"], n_permutations=3)
    assert res["cv"]["group_by"] == "unit" and res["cv"]["n_groups"] == 36
    assert any(w["code"] == "plsda_units_may_not_be_independent" for w in res["warnings"])


def test_spec_v2_passes_the_design_block_and_preprocess_split_units_are_grouped():
    """spec v2 の design.block が PLS-DA の独立な単位になり、前処理で属性に分かれた個体の単位はまとめて分割する。"""
    from metabolomix.analysis.stats.data import SampleInfo, StudyData
    from metabolomix.analysis.stats.spec_v2 import normalize_spec_v2, run_spec_v2

    rng = np.random.default_rng(11)
    samples, rows = [], []
    for cls in ("A", "B"):
        for s in range(6):
            base = rng.normal(0.0, 1.0, 20)
            for t in range(3):
                samples.append(SampleInfo(sample_id=f"{cls}{s}_t{t}", role="sample", group=cls, batch="b1",
                                          injection_order=None, subject_id=f"{cls}{s}", include=True,
                                          attributes={"time": f"t{t}"}))
                rows.append(2.0 ** (20.0 + base + rng.normal(0.0, 0.1, 20)))
    study = StudyData(values=np.array(rows), detected=None, sample_ids=[s.sample_id for s in samples],
                      feature_ids=[f"f{j}" for j in range(20)], samples=samples, quant=None, dataset_hash="h")
    off = {"drift": {"mode": "off"}, "normalize": {"method": "none"},
           "filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None}}
    spec = {"schema": "analysis-spec.v2", "preprocess": off,
            "comparisons": [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"}],
            "analyses": [{"analysis_id": "pl", "kind": "plsda", "n_permutations": 40}]}
    out = run_spec_v2(normalize_spec_v2(spec), study)
    res = out["analyses"][0]["result"]
    assert res["n_units"] == 36                                      # 時点で分かれた単位
    assert res["cv"]["group_by"] == "subject_id" and res["cv"]["n_groups"] == 12
    assert res["permutation"]["p"] > 0.1
    blocked = run_spec_v2(normalize_spec_v2({**spec, "design": {"factors": ["group"], "block": "subject_id"},
                                             "analyses": [{"analysis_id": "pl", "kind": "plsda",
                                                           "n_permutations": 3}]}), study)
    assert blocked["analyses"][0]["result"]["cv"]["group_by"] == "subject_id"
