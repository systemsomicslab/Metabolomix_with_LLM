"""analysis-spec.v1 の外部の表（MAF・mwTab・xlsx）・`correlate`・和の寄与・群の中で一定の印・
効かない設定の notes・aggregate と features.feature_ids の併用の誤り。fixture は合成のファイル。"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from metabolomix.analysis import analysis_spec as A
from metabolomix.analysis import spec_statistics as S
from metabolomix.core.atomic_io import DomainError
from tests.external_table_fixtures import write_maf, write_mwtab, write_xlsx

SAMPLES = ["HC1", "HC2", "HC3", "RA1", "RA2", "RA3", "QC1", "QC2"]
VALUES = {
    "F1": [10, 12, 11, 20, 22, 21, 15, 16],
    "F2": [5, 6, 5, 5, 6, 5, 5, 5],
    "F3": [0, 0, 0, 1, 0, 0, 0, 0],          # HC の 3 試料がすべて 0（群の中で一定）
    "F4": [100, 110, 90, 50, 55, 45, 75, 76],
}
#: 外部の表の値（試料は SAMPLES の順）。
GLUCOSE = [1, 3, 2, 5, 8, 6, 4, 7]
LACTATE = [9, 8, 7, 6, 5, 4, 3, 2]


def _write(path, header, rows):
    path.write_text("\n".join(["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]) + "\n",
                    encoding="utf-8")
    return path


@pytest.fixture
def files(tmp_path):
    matrix = _write(tmp_path / "matrix.tsv", ["feature_id", *SAMPLES],
                    [[f, *v] for f, v in VALUES.items()])
    groups = _write(tmp_path / "groups.tsv", ["sample", "group", "role"],
                    [[s, s[:2] if not s.startswith("QC") else "", "qc" if s.startswith("QC") else "sample"]
                     for s in SAMPLES])
    mwtab = write_mwtab(tmp_path / "paper.txt", SAMPLES,
                        [("Glucose", GLUCOSE), ("Lactate", LACTATE)])
    return {"matrix": matrix, "groups": groups, "mwtab": mwtab, "dir": tmp_path}


def _spec(files, analyses=None, **extra):
    spec = {"input": {"source": "tsv", "matrix_path": str(files["matrix"])},
            "groups": {"source": "tsv", "path": str(files["groups"]), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": analyses or [{"analysis_id": "main", "test": {"kind": "student"}}]}
    spec.update(extra)
    return spec


def _run(spec):
    return A.run_spec(A.normalize_spec(spec))


def _rows(result, analysis=0):
    return {r["feature_id"]: r for r in result["analyses"][analysis]["rows"]}


# ---------- Spearman と群の中で一定の印（spec_statistics） ----------

def test_spearman_matches_scipy_with_ties_and_drops_non_finite_pairs():
    x = [1.0, 2.0, 2.0, 4.0, 5.0, 7.0, math.nan]
    y = [2.0, 1.0, 4.0, 3.0, 6.0, 5.0, 1.0]
    out = S.spearman(x, y)
    expected = stats.spearmanr(x[:6], y[:6])
    assert out["n"] == 6 and out["reason"] is None
    assert out["rho"] == pytest.approx(expected.statistic)
    assert out["p"] == pytest.approx(expected.pvalue)
    perfect = S.spearman([1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6])
    assert perfect["rho"] == 1.0 and perfect["p"] == 0.0     # scipy と同じ（|ρ|=1 は p=0）
    reverse = S.spearman([1, 2, 3, 4, 5], [10, 9, 8, 7, 6])
    expected = stats.spearmanr([1, 2, 3, 4, 5], [10, 9, 8, 7, 6])
    assert reverse["rho"] == pytest.approx(-1.0)
    assert reverse["p"] == pytest.approx(expected.pvalue, abs=1e-12)


def test_spearman_too_few_pairs_and_constant_are_null_with_reason():
    assert S.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == {
        "n": 4, "rho": None, "p": None, "reason": "too_few_pairs"}
    assert S.spearman([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])["reason"] == "constant"


def test_constant_groups():
    assert S.constant_groups({"A": [0, 0, 0], "B": [1, 0, 0]}) == ["A"]
    assert S.constant_groups({"A": [5, math.nan, 5], "B": [math.nan, 2, math.nan]}) == ["A"]
    assert S.constant_groups({"A": [1, 2], "B": [3, 4]}) is None


def test_constant_in_group_marks_rows_without_changing_the_test(files):
    result = _run(_spec(files))
    rows = _rows(result)
    assert rows["F3"]["constant_in_group"] == ["HC"]
    assert rows["F1"]["constant_in_group"] is None
    assert list(rows["F1"])[-1] == "constant_in_group"        # 行の最後の欄


# ---------- aggregate（和）の寄与・E9 ----------

def test_sum_aggregate_returns_member_contributions(files):
    spec = _spec(files, aggregate={"kind": "sum", "label": "total", "feature_ids": ["F1", "F2", "F4"]})
    info = _run(spec)["analyses"][0]["aggregate"]
    assert info["n_samples_for_shares"] == 6                 # 比べる 2 群の試料だけ（QC は入れない）
    assert [c["feature_id"] for c in info["contributions"]] == ["F1", "F2", "F4"]
    block = np.array([VALUES[f][:6] for f in ("F1", "F2", "F4")], dtype=float)
    shares = block / block.sum(axis=0)
    by_id = {c["feature_id"]: c for c in info["contributions"]}
    for k, fid in enumerate(("F1", "F2", "F4")):
        assert by_id[fid]["share_median"] == pytest.approx(np.median(shares[k]))
        assert by_id[fid]["share_median_ref"] == pytest.approx(np.median(shares[k][:3]))
        assert by_id[fid]["share_median_test"] == pytest.approx(np.median(shares[k][3:]))
        assert by_id[fid]["share_max"] == pytest.approx(shares[k].max())
    assert "_block" not in info
    # 試料ごとの割合の和は 1（寄与は和を分けたもの）
    assert np.allclose(shares.sum(axis=0), 1.0)


def test_ratio_aggregate_has_no_contributions(files):
    spec = _spec(files, aggregate={"kind": "ratio", "label": "F1/F4", "numerator": "F1",
                                   "denominator": "F4"})
    info = _run(spec)["analyses"][0]["aggregate"]
    assert "contributions" not in info and "_block" not in info


def test_aggregate_with_feature_ids_is_an_error(files):
    spec = _spec(files, aggregate={"kind": "sum", "label": "total", "feature_ids": ["F1", "F2"]},
                 features={"feature_ids": ["F1"]})
    with pytest.raises(DomainError) as info:
        A.normalize_spec(spec)
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == "features"
    assert "aggregate" in info.value.message


# ---------- 効かない設定は notes ----------

def test_covariates_on_a_non_ols_test_are_dropped_with_a_note(files):
    covariates = _write(files["dir"] / "cov.tsv", ["sample", "age"],
                        [[s, 30 + i] for i, s in enumerate(SAMPLES)])
    cov = {"path": str(covariates), "columns": [{"name": "age", "type": "numeric"}]}
    with_cov = _spec(files, [{"analysis_id": "main", "test": {"kind": "welch", "covariates": ["age"]}}],
                     covariates=cov)
    notes: list[str] = []
    normalized = A.normalize_spec(with_cov, notes=notes)
    assert normalized["analyses"][0]["test"]["covariates"] == []
    assert len(notes) == 1 and "analyses[0].test.covariates" in notes[0] and "ols" in notes[0]
    plain = A.normalize_spec(_spec(files, [{"analysis_id": "main", "test": {"kind": "welch"}}],
                                   covariates=cov))
    assert A.rule_hash(normalized) == A.rule_hash(plain)
    # 知らない共変量はこれまでどおり誤り（効かない設定でも名前は確かめる、ではなく ols でだけ確かめる）
    ols = _spec(files, [{"analysis_id": "main", "test": {"kind": "ols", "covariates": ["sex"]}}],
                covariates=cov)
    with pytest.raises(DomainError):
        A.normalize_spec(ols)


def test_recipe_settings_not_in_steps_are_notes_not_errors(files):
    notes: list[str] = []
    spec = _spec(files, [{"analysis_id": "main", "test": {"kind": "student"},
                          "matrix": {"steps": [], "normalize": {"method": "pqn"}}}])
    normalized = A.normalize_spec(spec, notes=notes)
    assert normalized["analyses"][0]["matrix"] == {"schema": "matrix-recipe.v2", "steps": []}
    assert notes and notes[0].startswith("analyses[0].matrix:")


def test_alignment_qc_on_external_input_is_a_note(files):
    notes: list[str] = []
    normalized = A.normalize_spec(_spec(files, features={"all": True, "alignment_qc": "x.json"}),
                                  notes=notes)
    assert "alignment_qc" not in normalized["features"]
    assert any("外部の行列" in n for n in notes)


# ---------- 外部の入力（MAF・mwTab・xlsx） ----------

def test_maf_input_reads_like_the_tsv_and_reports_unreadable_cells(files):
    maf_rows = [(f, [*v[:2], "ND" if f == "F2" else v[2], *v[3:]]) for f, v in VALUES.items()]
    maf = write_maf(files["dir"] / "m_paper.tsv", SAMPLES, maf_rows)
    result = _run(_spec(files, input={"source": "maf", "path": str(maf)}))
    assert result["input_read"] == {"source": "maf", "n_features": 4, "n_samples": 8,
                                    "n_unreadable_cells": 1, "unreadable_examples": ["ND"]}
    rows = _rows(result)
    assert sorted(rows) == ["row_1:F1", "row_2:F2", "row_3:F3", "row_4:F4"]
    tsv_rows = _rows(_run(_spec(files)))
    assert rows["row_1:F1"]["p"] == pytest.approx(tsv_rows["F1"]["p"])
    assert rows["row_2:F2"]["n_ref"] == 2                    # ND の 1 セルは欠測
    assert len(result["input_fingerprint"]["files"]["matrix"]) == 64
    # TSV の入力には input_read を作らない
    assert "input_read" not in _run(_spec(files))


def test_mwtab_and_xlsx_inputs(files):
    mwtab = write_mwtab(files["dir"] / "mw.txt", SAMPLES, list(VALUES.items()))
    via_mwtab = _rows(_run(_spec(files, input={"source": "mwtab", "path": str(mwtab)})))
    sheet = [["Supplementary Table"], ["Metabolite", *SAMPLES],
             *[[f, *v] for f, v in VALUES.items()]]
    xlsx = write_xlsx(files["dir"] / "supp.xlsx", {"S1": sheet})
    result = _run(_spec(files, input={"source": "xlsx", "path": str(xlsx)}))
    via_xlsx = _rows(result)
    assert result["spec"]["input"] == {"source": "xlsx", "path": str(xlsx), "sheet": None,
                                       "feature_column": "A", "first_sample_column": "B"}
    tsv = _rows(_run(_spec(files)))
    for fid in VALUES:
        assert via_mwtab[fid]["p"] == pytest.approx(tsv[fid]["p"], nan_ok=True)
        assert via_xlsx[fid]["p"] == pytest.approx(tsv[fid]["p"], nan_ok=True)


def test_external_input_rule_hash_ignores_the_path(files, tmp_path):
    a = write_mwtab(files["dir"] / "a.txt", SAMPLES, list(VALUES.items()))
    b = write_mwtab(tmp_path / "b.txt", SAMPLES, list(VALUES.items()))
    ha = A.rule_hash(A.normalize_spec(_spec(files, input={"source": "mwtab", "path": str(a)})))
    hb = A.rule_hash(A.normalize_spec(_spec(files, input={"source": "mwtab", "path": str(b)})))
    assert ha == hb


@pytest.mark.parametrize("raw, field", [
    ({"source": "maf"}, "input"),
    ({"source": "xlsx", "path": "x.xlsx", "feature_column": "C", "first_sample_column": "B"}, "input"),
    ({"source": "xlsx", "path": "x.xlsx", "feature_column": "1"}, "input.feature_column"),
    ({"source": "csv", "path": "x.csv"}, "input.source"),
])
def test_external_input_spec_errors(files, raw, field):
    with pytest.raises(DomainError) as info:
        A.normalize_spec(_spec(files, input=raw))
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == field


def test_external_input_needs_groups_tsv(files):
    spec = _spec(files, input={"source": "maf", "path": "m.tsv"})
    del spec["groups"]
    with pytest.raises(DomainError) as info:
        A.normalize_spec(spec)
    assert "input.source=maf" in info.value.message


# ---------- correlate ----------

def _corr(files, **extra):
    return {"external": {"source": "mwtab", "path": str(files["mwtab"])},
            "pairs": [{"feature_id": "F1", "external_id": "Glucose"},
                      {"feature_id": "F4", "external_id": "Lactate"}], **extra}


def test_correlate_input_stage_is_spearman_against_the_external_values(files):
    result = _run(_spec(files, correlate=_corr(files)))
    block = result["correlate"]
    assert block["stage"] == "input" and block["n_samples_matched"] == 8
    assert block["unmatched_samples"] == [] and block["unmatched_external_samples"] == []
    assert block["external"] == {"source": "mwtab", "n_features": 2, "n_samples": 8,
                                 "n_unreadable_cells": 0, "unreadable_examples": []}
    first, second = block["rows"]
    expected = stats.spearmanr(VALUES["F1"], GLUCOSE)
    assert first["analysis_id"] is None and first["merged_into"] is None
    assert (first["feature_id"], first["external_id"], first["n"]) == ("F1", "Glucose", 8)
    assert first["rho"] == pytest.approx(expected.statistic)
    assert first["p"] == pytest.approx(expected.pvalue)
    assert second["rho"] == pytest.approx(stats.spearmanr(VALUES["F4"], LACTATE).statistic)
    assert len(result["input_fingerprint"]["files"]["correlate.external"]) == 64
    # 報告の数字は変えない・rule_hash は correlate の有無で分かれる
    plain = _run(_spec(files))
    assert result["analyses"] == plain["analyses"]
    assert result["rule_hash"] != plain["rule_hash"]
    assert "correlate" not in plain


def test_correlate_sample_map_and_unmatched(files):
    external = [f"P{i}" for i in range(1, 10)]
    mwtab = write_mwtab(files["dir"] / "renamed.txt", external,
                        [("Glucose", [*GLUCOSE, 99]), ("Lactate", [*LACTATE, 99])])
    sample_map = _write(files["dir"] / "map.tsv", ["sample", "external_sample"],
                        [[s, f"P{i}"] for i, s in enumerate(SAMPLES[:7], start=1)])
    spec = _spec(files, correlate={"external": {"source": "mwtab", "path": str(mwtab)},
                                   "pairs": [{"feature_id": "F1", "external_id": "Glucose"}],
                                   "sample_map": {"path": str(sample_map)}})
    block = _run(spec)["correlate"]
    assert block["n_samples_matched"] == 7
    assert block["unmatched_samples"] == ["QC2"]            # 対応表に無い試料は使わない
    assert block["unmatched_external_samples"] == ["P8", "P9"]
    [row] = block["rows"]
    assert row["n"] == 7
    assert row["rho"] == pytest.approx(stats.spearmanr(VALUES["F1"][:7], GLUCOSE[:7]).statistic)


def test_correlate_matrix_stage_reads_each_analysis_matrix(files):
    analyses = [{"analysis_id": "raw", "test": {"kind": "student"}},
                {"analysis_id": "filtered", "test": {"kind": "student"},
                 "matrix": {"steps": ["detection_filter"],
                            "detection_filter": {"basis": "positive", "min_rate": {"sample": 0.5}}}}]
    corr = _corr(files, stage="matrix",
                 pairs=[{"feature_id": "F1", "external_id": "Glucose"},
                        {"feature_id": "F3", "external_id": "Lactate"}])
    block = _run(_spec(files, analyses, correlate=corr))["correlate"]
    assert block["stage"] == "matrix"
    by_key = {(r["analysis_id"], r["feature_id"]): r for r in block["rows"]}
    assert set(by_key) == {("raw", "F1"), ("raw", "F3"), ("filtered", "F1"), ("filtered", "F3")}
    assert by_key[("raw", "F1")]["rho"] == pytest.approx(
        stats.spearmanr(VALUES["F1"], GLUCOSE).statistic)
    # recipe が落とした特徴（F3）は not_in_matrix
    assert by_key[("filtered", "F3")] == {
        "analysis_id": "filtered", "feature_id": "F3", "external_id": "Lactate",
        "merged_into": None, "n": 0, "rho": None, "p": None, "reason": "not_in_matrix"}


def test_correlate_errors(files):
    unknown_ext = _corr(files, pairs=[{"feature_id": "F1", "external_id": "Citrate"}])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=unknown_ext))
    assert info.value.code == "FEATURE_NOT_FOUND" and info.value.details["feature_ids"] == ["Citrate"]
    unknown = _corr(files, pairs=[{"feature_id": "F9", "external_id": "Glucose"}])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate=unknown))
    assert info.value.code == "FEATURE_NOT_FOUND" and info.value.details["feature_ids"] == ["F9"]
    other = write_mwtab(files["dir"] / "other.txt", ["X1", "X2"], [("Glucose", [1, 2])])
    with pytest.raises(DomainError) as info:
        _run(_spec(files, correlate={"external": {"source": "mwtab", "path": str(other)},
                                     "pairs": [{"feature_id": "F1", "external_id": "Glucose"}]}))
    assert info.value.code == "SAMPLE_MISMATCH"


def test_correlate_too_few_matched_samples_is_null(files):
    sample_map = _write(files["dir"] / "few.tsv", ["sample", "external_sample"],
                        [[s, s] for s in SAMPLES[:4]])
    block = _run(_spec(files, correlate=_corr(files, sample_map={"path": str(sample_map)})))["correlate"]
    assert block["n_samples_matched"] == 4
    assert all(r["rho"] is None and r["reason"] == "too_few_pairs" for r in block["rows"])


@pytest.mark.parametrize("raw, field", [
    ({"external": {"source": "mwtab", "path": "x"}, "pairs": []}, "correlate.pairs"),
    ({"external": {"source": "mwtab", "path": "x"}, "pairs": [{"feature_id": "F1"}]},
     "correlate.pairs[0]"),
    ({"external": {"source": "mwtab", "path": "x"},
      "pairs": [{"feature_id": "F1", "external_id": "G"}], "stage": "after"}, "correlate.stage"),
    ({"external": {"source": "dataset"}, "pairs": [{"feature_id": "F1", "external_id": "G"}]},
     "correlate.external.source"),
])
def test_correlate_spec_errors(files, raw, field):
    with pytest.raises(DomainError) as info:
        A.normalize_spec(_spec(files, correlate=raw))
    assert info.value.code == "ANALYSIS_SPEC_INVALID"
    assert info.value.details["field"] == field
