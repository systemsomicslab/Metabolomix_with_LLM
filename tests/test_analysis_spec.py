"""analysis-spec.v1（spec 2026-10-22 §4 の 3・4・5・10・11・12）。外部の TSV で組む。"""
from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy import stats

from metabolomix.analysis import analysis_spec as A
from metabolomix.core.atomic_io import DomainError

SAMPLES = ["HC1", "HC2", "HC3", "RA1", "RA2", "RA3", "QC1", "QC2"]
VALUES = {
    "F1": [10, 12, 11, 20, 22, 21, 15, 16],
    "F2": [5, 6, 5, 5, 6, 5, 5, 5],
    "F3": [0, 0, 0, 1, 0, 0, 0, 0],          # 検定できない（log2 で値がほぼ無い）
    "F4": [100, 110, 90, 50, 55, 45, 75, 76],
}


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
    covariates = _write(tmp_path / "cov.tsv", ["sample", "sex", "age"],
                        [["HC1", "F", 30], ["HC2", "M", 40], ["HC3", "F", 35],
                         ["RA1", "M", 50], ["RA2", "F", ""], ["RA3", "M", 41],
                         ["OTHER", "F", 99]])
    order = _write(tmp_path / "order.tsv", ["sample", "order"],
                   [[s, i + 1] for i, s in enumerate(SAMPLES)])
    return {"matrix": matrix, "groups": groups, "covariates": covariates, "order": order}


def _spec(files, analyses, **extra):
    spec = {"input": {"source": "tsv", "matrix_path": str(files["matrix"])},
            "groups": {"source": "tsv", "path": str(files["groups"]), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": analyses}
    spec.update(extra)
    return spec


MAIN = {"analysis_id": "main", "test": {"kind": "student"},
        "thresholds": {"p": 0.05, "fc": 1.5}}


def _run(spec):
    return A.run_spec(A.normalize_spec(spec))


def _rows(result, analysis=0):
    return {r["feature_id"]: r for r in result["analyses"][analysis]["rows"]}


def test_defaults_filled_and_unknown_key_rejected(files):
    normalized = A.normalize_spec(_spec(files, [MAIN]))
    analysis = normalized["analyses"][0]
    assert analysis["test"] == {"kind": "student", "transform": "log2", "covariates": []}
    assert analysis["fold_change"] == "arithmetic_mean_ratio"
    assert analysis["correction"] == "bh" and analysis["bh_population"] == "matrix"
    assert analysis["matrix"] == {"schema": "matrix-recipe.v2", "steps": []}
    assert normalized["features"] == {"all": True}
    with pytest.raises(DomainError):
        A.normalize_spec(_spec(files, [MAIN], extra_key=1))
    with pytest.raises(DomainError):
        A.normalize_spec(_spec(files, [{**MAIN, "test": {"kind": "student", "tails": 1}}]))


def test_group_names_direction_and_numbers_on_every_row(files):
    result = _run(_spec(files, [MAIN]))
    assert result["comparison"] == "RA_vs_HC"
    row = _rows(result)["F1"]
    assert row["group_ref"] == "HC" and row["group_test"] == "RA"
    assert row["fc_test_over_ref"] == pytest.approx(21 / 11)
    assert row["log2fc_test_over_ref"] == pytest.approx(math.log2(21 / 11))
    assert row["p"] == pytest.approx(stats.ttest_ind(np.log2([20, 22, 21]), np.log2([10, 12, 11])).pvalue)
    assert row["direction"] == "test_higher" and row["meets_thresholds"] is True
    assert _rows(result)["F4"]["direction"] == "ref_higher"


def test_not_testable_excluded_from_bh(files):
    result = _run(_spec(files, [MAIN]))
    rows = _rows(result)
    assert rows["F3"]["status"] == "not_testable" and rows["F3"]["q"] is None
    assert result["analyses"][0]["bh_population_n"] == 3


def test_bh_population_matrix_vs_reported(files):
    spec = _spec(files, [MAIN, {**MAIN, "analysis_id": "rep", "bh_population": "reported"}],
                 features={"feature_ids": ["F1"]})
    result = _run(spec)
    whole, reported = _rows(result, 0)["F1"], _rows(result, 1)["F1"]
    assert reported["q"] == pytest.approx(reported["p"])
    assert whole["q"] > whole["p"]


def test_sensitivity_analyses_return_each_numbers(files):
    analyses = [
        {"analysis_id": "none", "test": {"kind": "student"}},
        {"analysis_id": "half_min", "matrix": {"steps": ["impute"]}, "test": {"kind": "student"}},
        {"analysis_id": "welch", "test": {"kind": "welch"}},
        {"analysis_id": "mw", "test": {"kind": "mannwhitney"}},
        {"analysis_id": "ols", "test": {"kind": "ols", "covariates": ["sex", "age"]}},
    ]
    spec = _spec(files, analyses, features={"feature_ids": ["F1", "F3"]},
                 covariates={"path": str(files["covariates"]),
                             "columns": [{"name": "sex", "type": "categorical", "reference": "F"},
                                         {"name": "age", "type": "numeric"}]})
    result = _run(spec)
    assert [a["analysis_id"] for a in result["analyses"]] == ["none", "half_min", "welch", "mw", "ols"]
    assert _rows(result, 0)["F3"]["status"] == "not_testable"
    assert _rows(result, 1)["F3"]["status"] == "tested"          # 半最小値で埋まる
    ols = result["analyses"][4]
    assert ols["covariate_excluded_samples"] == ["RA2"]          # 年齢の欠測を除く
    assert _rows(result, 4)["F1"]["n_test"] == 2
    assert _rows(result, 4)["F1"]["beta_log2"] is not None


def test_same_spec_same_json_and_rule_hash(files):
    spec = _spec(files, [MAIN])
    first = json.dumps(_run(spec), sort_keys=True, ensure_ascii=False)
    second = json.dumps(_run(spec), sort_keys=True, ensure_ascii=False)
    assert first == second


def test_rule_hash_ignores_paths_but_fingerprint_tracks_files(files, tmp_path):
    moved = tmp_path / "elsewhere"
    moved.mkdir()
    copy = moved / "matrix.tsv"
    copy.write_bytes(files["matrix"].read_bytes())
    a = _run(_spec(files, [MAIN]))
    b = _run({**_spec(files, [MAIN]), "input": {"source": "tsv", "matrix_path": str(copy)}})
    assert a["rule_hash"] == b["rule_hash"]
    assert a["input_fingerprint"] == b["input_fingerprint"]
    copy.write_text(files["matrix"].read_text(encoding="utf-8").replace("100", "101"),
                    encoding="utf-8")
    c = _run({**_spec(files, [MAIN]), "input": {"source": "tsv", "matrix_path": str(copy)}})
    assert c["input_fingerprint"]["files"]["matrix"] != a["input_fingerprint"]["files"]["matrix"]


def test_sample_mismatch_names_listed(files, tmp_path):
    groups = _write(tmp_path / "g2.tsv", ["sample", "group"],
                    [["HC1", "HC"], ["HC2", "HC"], ["HC3", "HC"], ["RA1", "RA"], ["RA2", "RA"],
                     ["XX9", "RA"]])
    spec = _spec(files, [MAIN], groups={"source": "tsv", "path": str(groups)})
    with pytest.raises(DomainError) as info:
        _run(spec)
    assert info.value.code == "SAMPLE_MISMATCH"
    assert set(info.value.details["missing"]) == {"RA3", "QC1", "QC2"}
    assert info.value.details["extra"] == ["XX9"]


def test_aggregate_sum_lists_features_and_dropped(files):
    spec = _spec(files, [{"analysis_id": "sum", "test": {"kind": "student"},
                          "matrix": {"steps": ["detection_filter"],
                                     "detection_filter": {"basis": "positive",
                                                          "min_rate": {"sample": 0.5}}}}],
                 aggregate={"kind": "sum", "label": "total", "feature_ids": ["F1", "F2", "F3"]})
    result = _run(spec)
    analysis = result["analyses"][0]
    assert analysis["aggregate"]["summed"] == ["F1", "F2"]
    assert analysis["aggregate"]["dropped"] == ["F3"]
    row = analysis["rows"][0]
    assert row["feature_id"] == "total"
    assert row["mean_ref"] == pytest.approx((15 + 18 + 16) / 3)


def test_aggregate_ratio(files):
    spec = _spec(files, [{"analysis_id": "ratio", "test": {"kind": "student"}}],
                 aggregate={"kind": "ratio", "label": "F1/F4", "numerator": "F1",
                            "denominator": "F4"})
    row = _run(spec)["analyses"][0]["rows"][0]
    assert row["mean_test"] == pytest.approx(np.mean([20 / 50, 22 / 55, 21 / 45]))


def test_drift_with_order_tsv(files):
    spec = _spec(files, [{"analysis_id": "d", "test": {"kind": "student"},
                          "matrix": {"steps": ["drift_correct"],
                                     "drift_correct": {"min_qc": 2}}}],
                 injection_order={"source": "tsv", "path": str(files["order"])})
    result = _run(spec)
    assert result["analyses"][0]["matrix"]["steps"][0]["step"] == "drift_correct"


def test_unknown_reported_feature_rejected(files):
    with pytest.raises(DomainError) as info:
        _run(_spec(files, [MAIN], features={"feature_ids": ["nope"]}))
    assert info.value.code == "FEATURE_NOT_FOUND"


def test_small_group_rejected(files, tmp_path):
    groups = _write(tmp_path / "g3.tsv", ["sample", "group", "role"],
                    [[s, "HC" if s == "HC1" else ("RA" if s.startswith("RA") else "X"),
                      "qc" if s.startswith("QC") else "sample"] for s in SAMPLES])
    spec = _spec(files, [MAIN], groups={"source": "tsv", "path": str(groups), "role_column": "role"})
    with pytest.raises(DomainError) as info:
        _run(spec)
    assert info.value.code == "STATISTIC_GROUP_TOO_SMALL"


def _metadata_dataset():
    from types import SimpleNamespace

    names = ["HC1-pos", "HC2-pos", "HC3-pos", "RA1-pos", "RA2-pos", "RA3-pos", "QC1-pos", "X-pos"]
    ids = ["HC1", "HC2", "HC3", "RA1", "RA2", "RA3", "QC1", "X"]
    rows = [{"sample_id": i, "role": "qc" if i.startswith("QC") else "sample",
             "group": "" if i.startswith("QC") else i[:2], "include": i != "X",
             "injection_order": k + 1} for k, i in enumerate(ids)]
    matrix = np.array([[10, 12, 11, 20, 22, 21, 15, 999],
                       [5, 6, 5, 5, 6, 5, 5, 999]], dtype=float)
    return SimpleNamespace(sample_names=names, sample_ids=ids, sample_metadata_rows=rows,
                           feature_matrix=matrix, feature_ids=["1", "2"], detected_mask=None,
                           feature_candidates={}, feature_annotations={}, source_files={},
                           dataset_id="ds_test")


def test_dataset_input_uses_metadata_groups_roles_include_and_order(monkeypatch):
    from metabolomix.analysis import result_state

    monkeypatch.setattr(result_state, "dataset_fingerprint", lambda ds: "fp")
    ds = _metadata_dataset()
    spec = {"comparison": {"reference": "HC", "test": "RA"},
            "injection_order": {"source": "metadata"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"}}]}
    result = A.run_spec(A.normalize_spec(spec), ds)
    row = {r["feature_id"]: r for r in result["analyses"][0]["rows"]}["1"]
    assert row["n_ref"] == 3 and row["n_test"] == 3            # QC と include=false の X は群に入らない
    assert row["mean_test"] == pytest.approx(21.0)
    assert result["input_fingerprint"]["dataset_hash"] == "fp"


def test_dataset_without_metadata_needs_groups_tsv():
    ds = _metadata_dataset()
    ds.sample_metadata_rows = None
    spec = {"comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"}}]}
    with pytest.raises(DomainError) as info:
        A.run_spec(A.normalize_spec(spec), ds)
    assert info.value.code == "GROUPS_UNAVAILABLE"
