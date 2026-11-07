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
    for f, cols in (exp.get("na_cells") or {}).items():    # make_limma_expected.R が NA にした値
        log2[[ss.index(c) for c in cols], exp["features"].index(f)] = np.nan
    prep = _prepared(log2, [info[s]["group"] for s in ss], batches=[info[s]["batch"] for s in ss],
                     subjects=[info[s]["subject"] for s in ss], feature_ids=exp["features"],
                     unit_ids=ss)
    return prep, exp


def _check_against_limma(res, exp, mapping, skip=None):
    """mapping: compare の対比 id -> limma の対比名。skip: {対比 id: 突き合わせない特徴}。"""
    for cid, name in mapping.items():
        ref = exp["contrast"][name]
        rows = [r for r in res["rows"] if r["contrast_id"] == cid]
        assert [r["feature_id"] for r in rows] == exp["features"]
        keep = [i for i, r in enumerate(rows) if r["feature_id"] not in (skip or {}).get(cid, ())]
        rows = [rows[i] for i in keep]
        assert all(r["status"] == "tested" for r in rows)

        def want(values):
            return approx(np.array([values[i] for i in keep], dtype=float), **TOL)
        assert np.array([r["log2fc"] for r in rows]) == want(ref["logFC"])
        assert np.array([r["statistic"] for r in rows]) == want(ref["t"])
        assert np.array([r["p"] for r in rows]) == want(ref["P.Value"])
        assert np.array([r["ci_low"] for r in rows]) == want(ref["CI.L"])
        assert np.array([r["ci_high"] for r in rows]) == want(ref["CI.R"])
        assert np.array([r["df"] for r in rows]) == want(exp["df_total"])
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


def test_cell_without_values_matches_limma_na_coefficients():
    """C が全部欠測の特徴（f05）と C が 1 試料だけの特徴（f12）: 当てはめと事前分布には入り（limma と同じ）、
    B_vs_A は検定し、C_vs_A は insufficient_units、omnibus は出さない。"""
    prep, exp = _limma_prepared("missing_cell")
    assert np.isnan(prep.log2[8:, exp["features"].index("f05")]).all()
    comps = [{"id": "BvsA", "factor": "group", "reference": "A", "test": "B"},
             {"id": "CvsA", "factor": "group", "reference": "A", "test": "C"}]
    res = compare(prep, factors=["group"], comparisons=comps)
    _check_against_limma(res, exp, {"BvsA": "BvsA", "CvsA": "CvsA"},
                         skip={"CvsA": {"f05", "f12"}})
    rows = _rows(res, "CvsA")
    for f, n_test in (("f05", 0), ("f12", 1)):
        assert rows[f]["status"] == "not_tested" and rows[f]["reason"] == "insufficient_units"
        assert rows[f]["n_reference"] == 4 and rows[f]["n_test"] == n_test
    om = {o["feature_id"]: o for o in res["omnibus"]}
    assert om["f05"] == {"feature_id": "f05", "f": None, "f_p": None, "f_q": None}
    assert om["f12"]["f"] is None
    others = [i for i, f in enumerate(exp["features"]) if f not in ("f05", "f12")]
    assert [om[exp["features"][i]]["f"] for i in others] == approx([exp["F"][i] for i in others], rel=1e-6)
    f_p = [om[exp["features"][i]]["f_p"] for i in others]
    assert [om[exp["features"][i]]["f_q"] for i in others] == approx(spec_statistics.bh(f_p))


def _three_groups_with_empty_c():
    rng = np.random.default_rng(21)
    x = rng.normal(20.0, 0.1, (12, 6))
    x[4:8, 0] += 2.0          # B − A = +2
    x[8:, 0] = np.nan         # C は全部欠測
    return _prepared(x, ["A"] * 4 + ["B"] * 4 + ["C"] * 4), x


THREE = [{"id": "B_vs_A", "factor": "group", "reference": "A", "test": "B"},
         {"id": "C_vs_A", "factor": "group", "reference": "A", "test": "C"}]


@pytest.mark.parametrize("test", ["moderated_t", "welch"])
def test_eligibility_is_per_contrast(test):
    prep, x = _three_groups_with_empty_c()
    res = compare(prep, factors=["group"], comparisons=THREE, test=test)
    r = _rows(res, "B_vs_A")["f00"]
    assert r["status"] == "tested" and r["verdict"] == "changed"
    assert r["log2fc"] == approx(x[4:8, 0].mean() - x[:4, 0].mean())
    c = _rows(res, "C_vs_A")["f00"]
    assert c["status"] == "not_tested" and c["reason"] == "insufficient_units"
    assert (c["n_reference"], c["n_test"]) == (4, 0)
    if test == "moderated_t":
        om = {o["feature_id"]: o for o in res["omnibus"]}
        assert om["f00"]["f"] is None and om["f00"]["f_q"] is None
        assert om["f01"]["f"] is not None
        # 独立の計算: t = (B − A) / (√s2_post · stdev_unscaled)、s2_post は結果の事前分布で縮めた分散
        a, b = x[:4, 0], x[4:8, 0]
        s2 = (((a - a.mean()) ** 2).sum() + ((b - b.mean()) ** 2).sum()) / 6
        d0, s0 = res["prior"]["df_prior"], res["prior"]["s2_prior"]
        s2_post = (d0 * s0 + 6 * s2) / (d0 + 6)
        assert r["statistic"] == approx((b.mean() - a.mean()) / math.sqrt(s2_post * (1 / 4 + 1 / 4)),
                                        rel=1e-9)
    else:
        assert r["p"] == spec_statistics.two_group(x[:4, 0], x[4:8, 0], kind="welch",
                                                   transform="none")["p"]


def test_inestimable_contrast_is_not_tested(monkeypatch):
    """検定できる対比でも、当てはめで係数が推定できなければ inestimable_contrast（防御の分岐）。"""
    from metabolomix.analysis.stats import compare as module

    real = module.lm_fit

    def broken(y, x, names):
        fit = real(y, x, names)
        fit.coefficients[0, 1] = np.nan
        fit.stdev_unscaled[0, 1] = np.nan
        return fit

    monkeypatch.setattr(module, "lm_fit", broken)
    prep, _ = _two_groups(4, 4, 5)
    r = _rows(compare(prep, factors=["group"], comparisons=BVA))["f00"]
    assert r["status"] == "not_tested" and r["reason"] == "inestimable_contrast"


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


def test_replicates_collapse_through_preprocess_then_group_too_small():
    """preprocess() → compare(): 1 匹 3 注入 × 各群 1 匹は、まとめた後の n=1 で止まる。"""
    from metabolomix.analysis.stats.data import SampleInfo, StudyData
    from metabolomix.analysis.stats.preprocess import preprocess

    samples = [SampleInfo(sample_id=f"{g}{k}", role="sample", group=g, batch="b1",
                          injection_order=None, subject_id=f"m_{g}", include=True)
               for g in ("A", "B") for k in range(3)]
    rng = np.random.default_rng(4)
    study = StudyData(values=rng.uniform(1e5, 2e5, (6, 4)), detected=None,
                      sample_ids=[s.sample_id for s in samples],
                      feature_ids=[f"f{j}" for j in range(4)], samples=samples)
    prep = preprocess(study, {"filters": {"blank_ratio": None, "detection": None, "qc_rsd": None,
                                          "d_ratio": None},
                              "normalize": {"method": "none"}, "drift": {"mode": "off"}})
    assert [u.n_injections for u in prep.units] == [3, 3]
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
    inside = ["f01", "f02", "f05", "f09"]
    assert [rows[f]["q"] for f in inside] == approx(spec_statistics.bh([rows[f]["p"] for f in inside]))
    assert all(rows[f]["q"] is None for f in rows if f not in inside)


def test_outside_bh_population_is_never_changed_or_unchanged():
    rng = np.random.default_rng(2)
    x = rng.normal(20.0, 0.1, (40, 10))
    x[20:, 0] += 2.0                                         # 大きな差
    x[20:, 1] += x[:20, 1].mean() - x[20:, 1].mean()         # 差 0・区間は狭い
    prep = _prepared(x, ["A"] * 20 + ["B"] * 20)
    res = compare(prep, factors=["group"], comparisons=BVA, bh_population="reported",
                  reported_feature_ids=["f02", "f03"])
    rows = _rows(res)
    for f in ("f00", "f01"):
        assert rows[f]["status"] == "tested" and rows[f]["q"] is None and rows[f]["p"] is not None
        assert rows[f]["verdict"] == "inconclusive"
        assert rows[f]["verdict_reason"] == "outside_bh_population"
    # 同じ特徴は母集団の中なら changed / unchanged になる
    full = _rows(compare(prep, factors=["group"], comparisons=BVA))
    assert full["f00"]["verdict"] == "changed" and full["f01"]["verdict"] == "unchanged"


def test_verdict_boundaries_are_strict_for_equivalence_and_inclusive_for_q():
    from metabolomix.analysis.stats.compare import _verdict

    eq = math.log2(1.2)
    base = {"status": "tested", "quantity_judgeable": True, "q": 0.5}
    assert _verdict({**base, "ci_low": -0.1, "ci_high": eq}, 0.05, eq) == ("inconclusive", "neither")
    assert _verdict({**base, "ci_low": -eq, "ci_high": 0.1}, 0.05, eq) == ("inconclusive", "neither")
    assert _verdict({**base, "ci_low": -0.1, "ci_high": 0.1}, 0.05, eq) ==         ("unchanged", "ci_within_equivalence")
    # q が閾値ちょうどなら changed（compare を通して: 閾値を行の q そのものにする）
    prep, _ = _two_groups(4, 4, 10, shift={0: 1.0})
    q = _rows(compare(prep, factors=["group"], comparisons=BVA))["f00"]["q"]
    r = _rows(compare(prep, factors=["group"], comparisons=BVA, q_threshold=q))["f00"]
    assert r["q"] == q and r["verdict"] == "changed"
    r = _rows(compare(prep, factors=["group"], comparisons=BVA, q_threshold=q * (1 - 1e-12)))["f00"]
    assert r["verdict"] != "changed"


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
    sizes = {("wt", "ctl"): 3, ("wt", "drug"): 4, ("ko", "ctl"): 5, ("ko", "drug"): 3}
    attrs = [{"geno": g, "treat": t} for (g, t), n in sizes.items() for _ in range(n)]
    x = rng.normal(20.0, 0.1, (len(attrs), 8))
    for i, a in enumerate(attrs):
        if a["geno"] == "ko" and a["treat"] == "drug":
            x[i, 0] += 2.0                                  # 交互作用 = 2
    prep = _prepared(x, ["g"] * len(attrs), attributes=attrs)
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
    # 交互作用の reference は係数が負の cell、test は正の cell
    assert r["reference"] == "wt:drug + ko:ctl" and r["test"] == "ko:drug + wt:ctl"
    assert r["n_reference"] == 4 + 5 and r["n_test"] == 3 + 3
    s = _rows(res, "drug_in_wt")["f00"]
    assert (s["reference"], s["test"], s["n_reference"], s["n_test"]) == ("ctl", "drug", 3, 4)
    assert set(res["presence"]["bh_n_by_contrast"]) <= {"drug_in_wt"}
    assert res["omnibus"] is not None and len(res["omnibus"]) == 8
    assert res["group_n"] == {"ko:ctl": 5, "ko:drug": 3, "wt:ctl": 3, "wt:drug": 4}


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


# --- 数値の共変量の尺度（F2） -----------------------------------------------------

def _protein_case():
    rng = np.random.default_rng(1)
    x = rng.normal(20, 0.3, (10, 20))
    cov = rng.normal(5e7, 1e6, 10)
    return x, cov


def test_large_scale_numeric_covariate_runs_and_matches_prestandardized():
    """桁の大きい共変量でも compare が通り、手で標準化した共変量と同じ B−A の log2fc・t になる。"""
    x, cov = _protein_case()
    groups = ["A"] * 5 + ["B"] * 5
    raw = _prepared(x, groups, attributes=[{"protein": repr(float(c))} for c in cov])
    z = (cov - cov.mean()) / cov.std()
    std = _prepared(x, groups, attributes=[{"protein": repr(float(c))} for c in z])
    covs = [{"name": "protein", "type": "numeric"}]
    a = _rows(compare(raw, factors=["group"], comparisons=BVA, covariates=covs))
    b = _rows(compare(std, factors=["group"], comparisons=BVA, covariates=covs))
    assert all(r["status"] == "tested" for r in a.values())
    for fid in a:
        assert a[fid]["log2fc"] == approx(b[fid]["log2fc"], rel=1e-9, abs=1e-12)
        assert a[fid]["statistic"] == approx(b[fid]["statistic"], rel=1e-9, abs=1e-12)
        assert a[fid]["p"] == approx(b[fid]["p"], rel=1e-9, abs=1e-15)


def test_lm_fit_value_error_becomes_design_rank_deficient(monkeypatch):
    """lm_fit の ValueError（階数の判定）は DESIGN_RANK_DEFICIENT で返す（素の ValueError を漏らさない）。"""
    from metabolomix.analysis.stats import compare as module

    def broken(y, x, names):
        raise ValueError("design が列フルランクでない（推定できない係数がある）")

    monkeypatch.setattr(module, "lm_fit", broken)
    prep, _ = _two_groups(4, 4, 5)
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA)
    assert e.value.code == "DESIGN_RANK_DEFICIENT"


def test_group_spanning_batches_is_complete_confounding_in_compare():
    """A は b1 だけ、B は b2・b3: compare は allow_confounded なしで STATISTIC_CONFOUNDED。"""
    prep, _ = _two_groups(4, 4, 5, batches=["b1"] * 4 + ["b2", "b2", "b3", "b3"])
    with pytest.raises(DomainError) as e:
        compare(prep, factors=["group"], comparisons=BVA)
    assert e.value.code == "STATISTIC_CONFOUNDED"
    res = compare(prep, factors=["group"], comparisons=BVA, allow_confounded=True)
    assert res["confounding"]["batch_status"] == "complete"



# --- 対応のある計画が preprocess を通る（F5）・警告（F5・F6） --------------------------

def test_paired_design_through_preprocess_matches_hand_built_units():
    """同じ 4 個体 × 群 A・B × 2 注入 → 8 単位。block=subject_id の比較が、手で組んだ行列と一致する。"""
    from metabolomix.analysis.stats.data import SampleInfo, StudyData
    from metabolomix.analysis.stats.preprocess import preprocess

    subjects, g = ("m0", "m1", "m2", "m3"), 6
    samples = [SampleInfo(sample_id=f"{m}_{grp}_{r}", role="sample", group=grp, batch="b1",
                          injection_order=None, subject_id=m, include=True)
               for m in subjects for grp in ("A", "B") for r in (0, 1)]
    rng = np.random.default_rng(7)
    subj_eff = {m: rng.normal(0, 1.0, g) for m in subjects}
    log_vals = np.array([20 + subj_eff[s.subject_id] + (0.8 if s.group == "B" else 0.0)
                         + rng.normal(0, 0.1, g) for s in samples])
    study = StudyData(values=np.exp2(log_vals), detected=None,
                      sample_ids=[s.sample_id for s in samples],
                      feature_ids=[f"f{j:02d}" for j in range(g)], samples=samples)
    prep = preprocess(study, {"filters": {"blank_ratio": None, "detection": None, "qc_rsd": None,
                                          "d_ratio": None},
                              "normalize": {"method": "none"}, "drift": {"mode": "off"}})
    assert len(prep.units) == 8
    res = compare(prep, factors=["group"], comparisons=BVA, block="subject_id")
    hand_log2 = log_vals.reshape(8, 2, g).mean(axis=1)
    hand = _prepared(hand_log2, ["A", "B"] * 4, subjects=[m for m in subjects for _ in "AB"])
    ref = compare(hand, factors=["group"], comparisons=BVA, block="subject_id")
    assert res["design"]["columns"] == ref["design"]["columns"]
    a, b = _rows(res), _rows(ref)
    for fid in a:
        assert a[fid]["status"] == "tested"
        for key in ("log2fc", "statistic", "p", "ci_low", "ci_high"):
            assert a[fid][key] == approx(b[fid][key], rel=1e-9, abs=1e-12)


def test_units_missing_model_values_names_variables_and_units():
    prep, _ = _two_groups(3, 3, 4, attributes=[{"age": "30"}, {}, {"age": "40"}, {"age": "35"},
                                               {"age": "50"}, {"age": "45"}])
    res = compare(prep, factors=["group"], comparisons=BVA,
                  covariates=[{"name": "age", "type": "numeric"}])
    w = next(w for w in res["warnings"] if w["code"] == "units_missing_model_values")
    assert w["unit_ids"] == ["u01"] and w["variables"] == ["age"]
    assert "age" in w["message"] and "u01" in w["message"]


def test_annotated_population_empty_warns():
    prep, _ = _two_groups(4, 4, 5, annotated={"not_in_matrix"})
    res = compare(prep, factors=["group"], comparisons=BVA, bh_population="annotated")
    w = next(w for w in res["warnings"] if w["code"] == "bh_population_empty")
    assert w["contrasts"] == ["BvsA"] and w["n_annotated"] == 1 and w["n_annotated_in_matrix"] == 0


def test_reported_ids_not_in_matrix_warn():
    prep, _ = _two_groups(4, 4, 5)
    res = compare(prep, factors=["group"], comparisons=BVA, bh_population="reported",
                  reported_feature_ids=["f00", "nope", "zzz"])
    w = next(w for w in res["warnings"] if w["code"] == "reported_features_not_in_matrix")
    assert w["n_reported"] == 3 and w["n_not_in_matrix"] == 2 and w["examples"] == ["nope", "zzz"]
