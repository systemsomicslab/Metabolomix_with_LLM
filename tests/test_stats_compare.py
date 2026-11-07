"""比較のエンジン `compare`（統計エンジン spec §5〜§9）。

`PreparedMatrix` はテストの中で直接組み立てる（`preprocess()` を通さない）。
limma の参照値は `tests/test_stats_limma.py` と同じ TSV / JSON を読む。
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest
from pytest import approx
from scipy import stats

from metabolomix.analysis import spec_statistics
from metabolomix.analysis.stats.compare import BH_POPULATIONS, EFFECT_SIZES, TESTS, compare
from metabolomix.analysis.stats.preprocess import PreparedMatrix, UnitInfo
from metabolomix.core.atomic_io import DomainError

REF = Path(__file__).parent / "stats_reference"
EXPECTED = json.loads((REF / "limma_expected.json").read_text(encoding="utf-8"))
TOL = dict(rel=1e-6, abs=1e-10)


# --- fixture の組み立て ---------------------------------------------------------

def _prepared(log2, groups, *, batches=None, subjects=None, detected=None, feature_ids=None,
              attributes=None, annotated=(), injection_order=None, n_injections=None,
              unit_ids=None) -> PreparedMatrix:
    """個体 × 特徴の log2 を UnitInfo で包む。"""
    log2 = np.asarray(log2, dtype=float)
    n, g = log2.shape
    feature_ids = feature_ids or [f"f{j:02d}" for j in range(g)]
    units = []
    for i in range(n):
        uid = unit_ids[i] if unit_ids else f"u{i:02d}"
        k = n_injections[i] if n_injections else 1
        units.append(UnitInfo(
            unit_id=uid, group=groups[i], batch=batches[i] if batches else "b1",
            subject_id=subjects[i] if subjects else None,
            injection_order=injection_order[i] if injection_order else None,
            n_injections=k, sample_ids=[f"{uid}_inj{r}" for r in range(k)],
            attributes=dict(attributes[i]) if attributes else {}))
    return PreparedMatrix(
        matrix_id="mx_test", log2=log2, units=units, feature_ids=list(feature_ids),
        detected_units=None if detected is None else np.asarray(detected, bool),
        qc_log2=np.empty((0, g)), qc_ids=[], steps=[], config={}, quant=None,
        dataset_hash=None, subject_source="injection", warnings=[], feature_classes={},
        annotated_features=frozenset(annotated))


def _two_groups(n_a, n_b, g, *, sd=0.3, shift=None, seed=0, **kw):
    rng = np.random.default_rng(seed)
    x = rng.normal(20.0, sd, (n_a + n_b, g))
    if shift:
        for j, d in shift.items():
            x[n_a:, j] += d
    return _prepared(x, ["A"] * n_a + ["B"] * n_b, **kw), x


BVA = [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"}]


def _rows(res, cid=None):
    return {r["feature_id"]: r for r in res["rows"] if cid is None or r["contrast_id"] == cid}


def _limma_inputs():
    with (REF / "limma_input.tsv").open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    samples = rows[0][1:]
    features = [r[0] for r in rows[1:]]
    values = np.array([[np.nan if v == "NA" else float(v) for v in r[1:]] for r in rows[1:]])
    with (REF / "limma_design.tsv").open(encoding="utf-8", newline="") as fh:
        info = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return features, samples, values, info


def _limma_prepared(design_id):
    features, samples, values, info = _limma_inputs()
    exp = EXPECTED[design_id]
    fi = [features.index(f) for f in exp["features"]]
    si = [samples.index(s) for s in exp["samples"]]
    log2 = values[np.ix_(fi, si)].T
    ss = exp["samples"]
    prep = _prepared(log2, [info[s]["group"] for s in ss], batches=[info[s]["batch"] for s in ss],
                     subjects=[info[s]["subject"] for s in ss], feature_ids=exp["features"],
                     unit_ids=ss)
    return prep, exp


def _check_against_limma(res, exp, mapping):
    """mapping: compare の対比 id -> limma の対比名。"""
    for cid, name in mapping.items():
        ref = exp["contrast"][name]
        rows = [r for r in res["rows"] if r["contrast_id"] == cid]
        assert [r["feature_id"] for r in rows] == exp["features"]
        assert all(r["status"] == "tested" for r in rows)
        assert np.array([r["log2fc"] for r in rows]) == approx(np.array(ref["logFC"]), **TOL)
        assert np.array([r["statistic"] for r in rows]) == approx(np.array(ref["t"]), **TOL)
        assert np.array([r["p"] for r in rows]) == approx(np.array(ref["P.Value"]), **TOL)
        assert np.array([r["ci_low"] for r in rows]) == approx(np.array(ref["CI.L"]), **TOL)
        assert np.array([r["ci_high"] for r in rows]) == approx(np.array(ref["CI.R"]), **TOL)
        assert np.array([r["df"] for r in rows]) == approx(np.array(exp["df_total"]), **TOL)
    assert res["prior"]["s2_prior"] == approx(float(exp["s2_prior"]), rel=1e-6)
    df_prior = math.inf if exp["df_prior"] == "Inf" else float(exp["df_prior"])
    assert res["prior"]["df_prior"] == approx(df_prior, rel=1e-6)


# --- limma の参照値 -------------------------------------------------------------

def test_moderated_t_two_groups_matches_limma_reference():
    prep, exp = _limma_prepared("cell_means")
    comps = [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"},
             {"id": "CvsA", "factor": "group", "reference": "A", "test": "C"}]
    res = compare(prep, factors=["group"], comparisons=comps)
    _check_against_limma(res, exp, {"BvsA": "BvsA", "CvsA": "CvsA"})
    assert res["group_n"] == {"A": 4, "B": 4, "C": 4}


def test_three_groups_omnibus_f_matches_limma():
    prep, exp = _limma_prepared("cell_means")
    comps = [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"},
             {"id": "CvsA", "factor": "group", "reference": "A", "test": "C"}]
    res = compare(prep, factors=["group"], comparisons=comps)
    om = res["omnibus"]
    assert [o["feature_id"] for o in om] == exp["features"]
    assert np.array([o["f"] for o in om]) == approx(np.array(exp["F"]), rel=1e-6)
    assert np.array([o["f_p"] for o in om]) == approx(np.array(exp["F_p_value"]), rel=1e-6)
    assert [o["f_q"] for o in om] == approx(spec_statistics.bh([o["f_p"] for o in om]))

    # 3 cell なら対比が 1 つでも omnibus を出す
    one = compare(prep, factors=["group"], comparisons=comps[:1])
    assert np.array([o["f"] for o in one["omnibus"]]) == approx(np.array(exp["F"]), rel=1e-6)
    # 2 cell・対比 1 つなら出さない
    two, _ = _two_groups(4, 4, 5)
    assert compare(two, factors=["group"], comparisons=BVA)["omnibus"] is None


def test_batch_covariate_matches_limma_reference():
    prep, exp = _limma_prepared("batch_covariate")
    res = compare(prep, factors=["group"],
                  comparisons=[{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"}],
                  covariates=[{"name": "batch", "type": "categorical"}])
    _check_against_limma(res, exp, {"BvsA": "BvsA"})
    assert res["design"]["columns"] == ["A", "B", "C", "batch:b2"]


def test_paired_block_matches_limma_reference():
    prep, exp = _limma_prepared("paired_block")
    res = compare(prep, factors=["group"],
                  comparisons=[{"id": "CvsB", "factor": "group", "reference": "B", "test": "C"}],
                  block="subject_id")
    _check_against_limma(res, exp, {"CvsB": "CvsB"})
    assert res["design"]["block"] == "subject_id"


# --- 判定 -----------------------------------------------------------------------

def test_verdicts_changed_unchanged_inconclusive():
    rng = np.random.default_rng(1)
    big = rng.normal(20.0, 0.1, (40, 30))
    big[20:, 0] += 2.0
    big[20:, 1] += big[:20, 1].mean() - big[20:, 1].mean()   # 差をちょうど 0 に
    res = compare(_prepared(big, ["A"] * 20 + ["B"] * 20), factors=["group"], comparisons=BVA)
    rows = _rows(res)
    assert rows["f00"]["verdict"] == "changed"
    assert rows["f00"]["verdict_reason"] == "q_below_threshold"
    assert rows["f01"]["verdict"] == "unchanged"
    assert rows["f01"]["verdict_reason"] == "ci_within_equivalence"
    eq = math.log2(1.2)
    assert -eq < rows["f01"]["ci_low"] and rows["f01"]["ci_high"] < eq

    small = rng.normal(20.0, 0.5, (6, 30))
    small[3:, 1] += small[:3, 1].mean() - small[3:, 1].mean()
    res = compare(_prepared(small, ["A"] * 3 + ["B"] * 3), factors=["group"], comparisons=BVA)
    r = _rows(res)["f01"]
    assert r["log2fc"] == approx(0.0, abs=1e-9)
    assert r["q"] > 0.05
    assert r["verdict"] == "inconclusive"           # 有意でない ≠ 変わっていない
    assert r["verdict_reason"] == "neither"


def test_quantity_flag_forces_inconclusive_even_if_significant():
    prep, x = _two_groups(4, 4, 20, sd=0.1, shift={0: 3.0})
    det = np.ones_like(x, dtype=bool)
    det[4:7, 0] = False          # B では 4 個体中 1 個体だけが実検出（残りは gap-fill）
    prep.detected_units = det
    res = compare(prep, factors=["group"], comparisons=BVA)
    r = _rows(res)["f00"]
    assert r["status"] == "tested" and r["q"] <= 0.05
    assert r["quantity_judgeable"] is False
    assert r["verdict"] == "inconclusive"
    assert r["verdict_reason"] == "quantity_not_judgeable"
    assert r["detected_reference"] == 4 and r["detected_test"] == 1
    assert _rows(res)["f01"]["quantity_judgeable"] is True
    # 閾値は引数で変えられる
    res = compare(prep, factors=["group"], comparisons=BVA, min_detected_fraction=0.25)
    assert _rows(res)["f00"]["verdict"] == "changed"


def test_feature_absent_in_one_group_is_presence_only():
    prep, x = _two_groups(6, 6, 10)
    prep.log2[6:, 0] = np.nan
    det = np.ones_like(x, dtype=bool)
    det[6:, 0] = False
    prep.detected_units = det
    res = compare(prep, factors=["group"], comparisons=BVA)
    r = _rows(res)["f00"]
    assert r["status"] == "not_tested" and r["reason"] == "insufficient_units"
    assert r["n_reference"] == 6 and r["n_test"] == 0
    assert r["log2fc"] is None and r["p"] is None and r["q"] is None
    assert r["presence_p"] == approx(spec_statistics.detection_fisher_p(6, 6, 0, 6))
    assert r["presence_q"] is not None
    assert r["verdict"] == "inconclusive" and r["verdict_reason"] == "not_tested"
    assert res["bh_population"]["n_by_contrast"] == {"BvsA": 9}


def test_no_detected_mask_runs_quantity_only():
    prep, _ = _two_groups(4, 4, 10, shift={0: 3.0})
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["presence"]["detection_available"] is False
    assert res["presence"]["reason"] == "detected_mask_unavailable"
    r = _rows(res)["f00"]
    assert r["status"] == "tested"
    assert r["quantity_judgeable"] is None
    assert r["presence_p"] is None and r["presence_q"] is None
    assert r["detected_reference"] is None and r["detected_test"] is None
    assert r["verdict"] == "changed"


def test_replicates_collapse_before_group_size_check():
    # 1 匹 3 注入 × 各群 1 匹: 個体でまとめた後は各群 n=1（n=3 と数えない）
    prep = _prepared(np.array([[20.0, 21.0], [22.0, 23.0]]), ["A", "B"], n_injections=[3, 3],
                     subjects=["m1", "m2"])
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA)
    assert e.value.code == "STATISTIC_GROUP_TOO_SMALL"
    assert e.value.details["cell_n"] == {"A": 1, "B": 1}


def test_presence_bh_is_separate_from_quantity_bh():
    prep, x = _two_groups(5, 5, 10)
    det = np.ones_like(x, dtype=bool)
    for j in (0, 1):
        prep.log2[5:, j] = np.nan
        det[5:, j] = False
    prep.detected_units = det
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["presence"]["bh_n_by_contrast"] == {"BvsA": 10}
    assert res["bh_population"]["n_by_contrast"] == {"BvsA": 8}
    rows = [r for r in res["rows"]]
    assert [r["presence_q"] for r in rows] == approx(spec_statistics.bh([r["presence_p"] for r in rows]))
    tested = [r for r in rows if r["status"] == "tested"]
    assert [r["q"] for r in tested] == approx(spec_statistics.bh([r["p"] for r in tested]))


def test_bh_population_reported_marks_reproduction_only():
    prep, _ = _two_groups(4, 4, 10, shift={0: 2.0})
    reported = ["f00", "f03", "f07", "nope"]
    res = compare(prep, factors=["group"], comparisons=BVA, bh_population="reported",
                  reported_feature_ids=reported)
    assert res["bh_population"]["name"] == "reported"
    assert res["bh_population"]["q_interpretation"] == "reproduction_only"
    assert res["bh_population"]["n_by_contrast"] == {"BvsA": 3}
    rows = _rows(res)
    inside = [rows[f] for f in ("f00", "f03", "f07")]
    assert [r["q"] for r in inside] == approx(spec_statistics.bh([r["p"] for r in inside]))
    assert all(rows[f]["q"] is None for f in rows if f not in reported)
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA, bh_population="reported")
    assert e.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_bh_population_annotated_uses_annotated_features():
    prep, _ = _two_groups(4, 4, 10, annotated={"f01", "f02", "f05", "f09"})
    res = compare(prep, factors=["group"], comparisons=BVA, bh_population="annotated")
    assert res["bh_population"]["name"] == "annotated"
    assert res["bh_population"]["q_interpretation"] == "discovery_fdr"
    assert res["bh_population"]["n_by_contrast"] == {"BvsA": 4}
    rows = _rows(res)
    assert rows["f00"]["q"] is None and rows["f01"]["q"] is not None


# --- 2 群の検定の選択肢 ----------------------------------------------------------

def test_welch_option_uses_spec_statistics_numbers():
    prep, x = _two_groups(4, 6, 5, shift={0: 1.0}, seed=3)
    res = compare(prep, factors=["group"], comparisons=BVA, test="welch")
    assert res["test"] == "welch" and res["prior"] is None and res["omnibus"] is None
    for j, r in enumerate(res["rows"]):
        a, b = x[:4, j], x[4:, j]
        ref = spec_statistics.two_group(a, b, kind="welch", transform="none")
        assert r["p"] == ref["p"] and r["statistic"] == ref["statistic"]
        va, vb = a.var(ddof=1) / 4, b.var(ddof=1) / 6
        df = (va + vb) ** 2 / (va ** 2 / 3 + vb ** 2 / 5)
        half = stats.t.ppf(0.975, df) * math.sqrt(va + vb)
        d = b.mean() - a.mean()
        assert r["log2fc"] == approx(d)
        assert r["df"] == approx(df)
        assert (r["ci_low"], r["ci_high"]) == approx((d - half, d + half))

    res = compare(prep, factors=["group"], comparisons=BVA, test="student")
    r = res["rows"][0]
    a, b = x[:4, 0], x[4:, 0]
    sp2 = (3 * a.var(ddof=1) + 5 * b.var(ddof=1)) / 8
    half = stats.t.ppf(0.975, 8) * math.sqrt(sp2 * (1 / 4 + 1 / 6))
    assert r["p"] == spec_statistics.two_group(a, b, kind="student", transform="none")["p"]
    assert r["df"] == 8
    assert (r["ci_low"], r["ci_high"]) == approx((b.mean() - a.mean() - half, b.mean() - a.mean() + half))

    # 2 群の検定は 1 因子・共変量なし・ブロックなしだけ
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA, test="welch",
                covariates=[{"name": "batch", "type": "categorical"}])
    assert e.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_mannwhitney_and_arithmetic_ratio_have_no_ci_and_never_unchanged():
    rng = np.random.default_rng(5)
    x = rng.normal(20.0, 0.05, (40, 6))
    x[20:, 1] += x[:20, 1].mean() - x[20:, 1].mean()
    prep = _prepared(x, ["A"] * 20 + ["B"] * 20)
    res = compare(prep, factors=["group"], comparisons=BVA, test="mannwhitney")
    for r in res["rows"]:
        assert r["ci_low"] is None and r["ci_high"] is None
        assert r["ci_unavailable_reason"] == "rank_test"
        assert r["df"] is None
        assert r["verdict"] != "unchanged"
    r = _rows(res)["f01"]
    assert r["p"] == spec_statistics.two_group(x[:20, 1], x[20:, 1], kind="mannwhitney",
                                               transform="none")["p"]
    assert r["verdict_reason"] == "no_ci"

    res = compare(prep, factors=["group"], comparisons=BVA, effect_size="arithmetic_mean_ratio")
    assert res["effect_size_definition"]["name"] == "arithmetic_mean_ratio"
    for j, r in enumerate(res["rows"]):
        expected = math.log2(np.mean(2 ** x[20:, j]) / np.mean(2 ** x[:20, j]))
        assert r["log2fc"] == approx(expected)
        assert r["ci_low"] is None and r["ci_unavailable_reason"] == "arithmetic_mean_ratio"
        assert r["verdict"] != "unchanged"


# --- 2 因子 -------------------------------------------------------------------

def test_two_by_two_interaction_contrast_runs():
    rng = np.random.default_rng(11)
    attrs, rows = [], []
    for geno in ("wt", "ko"):
        for treat in ("ctl", "drug"):
            for _ in range(4):
                attrs.append({"geno": geno, "treat": treat})
    x = rng.normal(20.0, 0.1, (16, 8))
    for i, a in enumerate(attrs):
        if a["geno"] == "ko" and a["treat"] == "drug":
            x[i, 0] += 2.0                                  # 交互作用 = 2
    prep = _prepared(x, ["g"] * 16, attributes=attrs)
    comps = [{"id": "drug_in_wt", "factor": "treat", "reference": "ctl", "test": "drug",
              "at": {"geno": "wt"}},
             {"id": "inter", "kind": "interaction", "factors": ["geno", "treat"],
              "levels": {"geno": ["wt", "ko"], "treat": ["ctl", "drug"]}}]
    res = compare(prep, factors=["geno", "treat"], comparisons=comps)
    mean = {(g, t): x[[i for i, a in enumerate(attrs) if a == {"geno": g, "treat": t}], 0].mean()
            for g in ("wt", "ko") for t in ("ctl", "drug")}
    r = _rows(res, "inter")["f00"]
    want = mean["ko", "drug"] - mean["wt", "drug"] - mean["ko", "ctl"] + mean["wt", "ctl"]
    assert r["log2fc"] == approx(want)
    assert r["verdict"] == "changed"
    assert r["presence_p"] is None and r["presence_q"] is None
    assert r["n_reference"] == 8 and r["n_test"] == 8
    assert set(res["presence"]["bh_n_by_contrast"]) <= {"drug_in_wt"}
    assert res["omnibus"] is not None and len(res["omnibus"]) == 8
    assert res["group_n"] == {"ko:ctl": 4, "ko:drug": 4, "wt:ctl": 4, "wt:drug": 4}


# --- 交絡 ----------------------------------------------------------------------

def test_complete_confounding_raises_unless_allowed():
    prep, _ = _two_groups(4, 4, 6, batches=["b1"] * 4 + ["b2"] * 4)
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA)
    assert e.value.code == "STATISTIC_CONFOUNDED"
    res = compare(prep, factors=["group"], comparisons=BVA, allow_confounded=True)
    assert res["confounding"]["batch_status"] == "complete"
    assert any(w["code"] == "batch_group_complete_confounding" for w in res["warnings"])


def test_result_always_has_batch_table():
    prep, _ = _two_groups(4, 4, 6, batches=["b1", "b2"] * 4)
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["confounding"]["batch_table"] == {"A": {"b1": 2, "b2": 2}, "B": {"b1": 2, "b2": 2}}
    prep, _ = _two_groups(4, 4, 6, batches=[None] * 8)
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert "batch_table" in res["confounding"]
    assert res["confounding"]["batch_status"] == "not_evaluable"


# --- 結果の定義と縁の場合 ---------------------------------------------------------

def test_result_states_definitions():
    prep, _ = _two_groups(4, 4, 6)
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["schema"] == "stats-compare.v1"
    assert res["matrix_id"] == "mx_test"
    assert res["test"] == "moderated_t"
    assert res["effect_size_definition"]["name"] == "log2_mean_difference"
    assert res["ci_level"] == 0.95
    assert res["bh_population"]["name"] == "matrix"
    assert res["bh_population"]["q_interpretation"] == "discovery_fdr"
    assert res["verdict_rule"] == {"q_threshold": 0.05, "equivalence_fc": 1.2,
                                   "equivalence_log2": approx(math.log2(1.2))}
    assert res["subject_source"] == "injection"
    assert res["design"] == {"columns": ["A", "B"], "factors": ["group"], "covariates": [],
                             "block": None}
    assert TESTS == ("moderated_t", "welch", "student", "mannwhitney")
    assert EFFECT_SIZES == ("log2_mean_difference", "arithmetic_mean_ratio")
    assert BH_POPULATIONS == ("matrix", "annotated", "reported")
    keys = {"feature_id", "contrast_id", "reference", "test", "status", "reason", "n_reference",
            "n_test", "detected_reference", "detected_test", "quantity_judgeable", "log2fc",
            "ci_low", "ci_high", "ci_unavailable_reason", "statistic", "df", "p", "q",
            "presence_p", "presence_q", "verdict", "verdict_reason"}
    assert all(set(r) == keys for r in res["rows"])
    assert res["rows"][0]["reference"] == "A" and res["rows"][0]["test"] == "B"


@pytest.mark.parametrize("kw", [{"test": "ttest"}, {"effect_size": "median"},
                                {"bh_population": "all"}, {"q_threshold": 0.0},
                                {"equivalence_fc": 1.0}, {"min_detected_fraction": 1.5}])
def test_invalid_specification(kw):
    prep, _ = _two_groups(4, 4, 6)
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA, **kw)
    assert e.value.code == "STATISTIC_SPECIFICATION_INVALID"


def test_units_missing_factor_values_are_excluded_with_warning():
    prep, _ = _two_groups(4, 4, 6)
    prep.units[0].group = None
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["group_n"] == {"A": 3, "B": 4}
    assert any(w["code"] == "units_missing_model_values" for w in res["warnings"])


def test_no_eligible_feature_returns_all_not_tested():
    prep, _ = _two_groups(3, 3, 3)
    prep.log2[3:5, :] = np.nan            # B は各特徴 1 個体だけ
    res = compare(prep, factors=["group"], comparisons=BVA)
    assert res["prior"] is None
    assert all(r["status"] == "not_tested" and r["reason"] == "insufficient_units"
               for r in res["rows"])


def test_zero_variance_feature_is_not_tested_not_infinite():
    x = np.full((6, 2), 20.0)
    x[3:, :] = 21.0                       # どの特徴も群の中で分散 0
    res = compare(_prepared(x, ["A"] * 3 + ["B"] * 3), factors=["group"], comparisons=BVA)
    for r in res["rows"]:
        assert r["status"] == "not_tested" and r["reason"] == "degenerate_variance"
        assert r["statistic"] is None and r["p"] is None
