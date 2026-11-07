"""analysis-spec.v1 の数字を黄金値で固定する（C §12）。

v1 の実行経路（`analysis_spec.normalize_spec` / `run_spec`）が出す数を、保存した黄金値と
突き合わせる。以降の配線で v1 の数字が動いたら、ここが落ちる。黄金値の作り直しは
`UPDATE_V1_GOLDEN=1` を付けて実行したときだけ（意図した変更のときに限る）。
入力はテスト自身が tmp_path に作る（乱数は legacy RandomState。版をまたいで同じ列を返す）。
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from metabolomix.analysis import analysis_spec as A
from metabolomix.analysis import recipe_matrix, spec_statistics

GOLDEN = Path(__file__).parent / "stats_reference" / "analysis_spec_v1_golden.json"

GROUP_OF = {**{f"HC{i}": "HC" for i in range(1, 5)}, **{f"RA{i}": "RA" for i in range(1, 5)},
            **{f"QC{i}": "" for i in range(1, 5)}}
SAMPLES = list(GROUP_OF)
FEATURES = [f"F{i}" for i in range(1, 9)]


def _values() -> np.ndarray:
    """(特徴 8 × 試料 12)。0・欠測・gap 的な小さな値を含む。"""
    rng = np.random.RandomState(20261107)
    base = np.array([100.0, 50.0, 200.0, 30.0, 80.0, 10.0, 5.0, 60.0])
    effect = np.array([1.0, 1.6, 0.7, 1.0, 2.2, 1.3, 1.0, 1.0])
    out = np.zeros((8, 12))
    for j, sample in enumerate(SAMPLES):
        for i in range(8):
            level = base[i] * (effect[i] if sample.startswith("RA") else 1.0)
            out[i, j] = round(level * float(np.exp(rng.normal(0, 0.15))) * (1 + 0.02 * j), 2)
    out[2, [0, 5]] = 0.0              # 0
    out[3, [1, 9]] = np.nan           # 空欄
    out[5, [2, 6, 8]] = 0.3           # gap 的な小さな値
    out[6, [0, 1, 2, 4, 5, 8, 10]] = 0.0   # 欠測が多い（KNN の重い行・検出率の篩で落ちる）
    out[7, 4] = 0.0
    return out


def _fmt(x: float) -> str:
    return "" if np.isnan(x) else repr(float(x))


def _write(path: Path, header, rows) -> Path:
    text = "\n".join(["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]) + "\n"
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return path


def _tsv_files(tmp_path: Path) -> dict:
    values = _values()
    matrix = _write(tmp_path / "matrix.tsv", ["feature_id", *SAMPLES],
                    [[f, *[_fmt(v) for v in values[i]]] for i, f in enumerate(FEATURES)])
    groups = _write(tmp_path / "groups.tsv", ["sample", "group", "role"],
                    [[s, GROUP_OF[s], "qc" if s.startswith("QC") else "sample"] for s in SAMPLES])
    order = _write(tmp_path / "order.tsv", ["sample", "order"],
                   [[s, 12 - i if i % 2 else i + 1] for i, s in enumerate(SAMPLES)])
    sex = ["F", "M"]
    ages = [30, 41, 35, 52, 47, 38, 44, 60, 33, 36, 39, 40]
    rows = [[s, sex[(i * 5 + i // 3) % 2], "" if s == "RA2" else ages[i]]
            for i, s in enumerate(SAMPLES)]
    covariates = _write(tmp_path / "cov.tsv", ["sample", "sex", "age"], rows)
    return {"matrix": matrix, "groups": groups, "order": order, "covariates": covariates}


def _tsv_spec(files: dict, analyses: list, **extra) -> dict:
    spec = {"input": {"source": "tsv", "matrix_path": str(files["matrix"])},
            "groups": {"source": "tsv", "path": str(files["groups"]), "role_column": "role"},
            "injection_order": {"source": "tsv", "path": str(files["order"])},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": analyses}
    spec.update(extra)
    return spec


def _dataset() -> SimpleNamespace:
    values = _values()
    ids = list(SAMPLES) + ["BLK1"]
    matrix = np.column_stack([values, np.full(8, 1.0)])
    rows = [{"sample_id": s, "role": "qc" if s.startswith("QC") else "sample",
             "group": GROUP_OF.get(s, ""), "include": s != "BLK1",
             "injection_order": k + 1} for k, s in enumerate(ids)]
    mask = np.isfinite(matrix) & (matrix > 0.5)
    mask[7, [3, 4]] = False
    return SimpleNamespace(
        sample_names=[f"{s}-pos" for s in ids], sample_ids=ids, sample_metadata_rows=rows,
        feature_matrix=np.nan_to_num(matrix, nan=0.0), feature_ids=[str(i) for i in range(1, 9)],
        detected_mask=mask,
        feature_candidates={"1": [{"chemical_name": "PC 16:0_18:1"}],
                            "2": [{"chemical_name": "PE 18:0_20:4"}],
                            "5": [{"chemical_name": "TG 50:2"}], "8": []},
        feature_annotations={}, quantification_measure="height", source_files={},
        dataset_id="ds_golden")


def _recipe_full() -> dict:
    return {"steps": ["detection_filter", "impute", "normalize", "drift_correct", "qc_rsd_filter"],
            "detection_filter": {"basis": "positive", "min_rate": {"sample": 0.6, "qc": 0.5}},
            "impute": {"method": "half_min"},
            "normalize": {"method": "pqn"},
            "drift_correct": {"span": 0.9, "min_qc": 4},
            "qc_rsd_filter": {"max_rsd": 0.5}}


def _case_tsv_student_recipe(tmp_path):
    files = _tsv_files(tmp_path)
    return _tsv_spec(files, [
        {"analysis_id": "student_full", "test": {"kind": "student"}, "matrix": _recipe_full(),
         "fold_change": "geometric_mean_ratio", "thresholds": {"p": 0.05, "q": 0.1, "fc": 1.5}},
        {"analysis_id": "student_none_corr", "test": {"kind": "student", "transform": "none"},
         "correction": "none"}]), None


def _case_tsv_welch_knn_tic(tmp_path):
    files = _tsv_files(tmp_path)
    recipe = {"steps": ["impute", "normalize"],
              "impute": {"method": "knn", "k": 3, "scale": "log2", "rowmax": 0.5},
              "normalize": {"method": "tic"}}
    return _tsv_spec(files, [
        {"analysis_id": "welch_knn", "test": {"kind": "welch"}, "matrix": recipe,
         "bh_population": "reported"}],
        features={"feature_ids": ["F1", "F2", "F5", "F7", "F8"]}), None


def _case_tsv_mannwhitney_median(tmp_path):
    files = _tsv_files(tmp_path)
    recipe = {"steps": ["detection_filter", "impute", "normalize", "qc_rsd_filter"],
              "detection_filter": {"basis": "positive", "min_rate": {"sample": 0.5}},
              "impute": {"method": "half_min", "missing": "nan"},
              "normalize": {"method": "median"},
              "qc_rsd_filter": {"max_rsd": 0.4}}
    return _tsv_spec(files, [
        {"analysis_id": "mw_median", "test": {"kind": "mannwhitney"}, "matrix": recipe,
         "bh_population": "matrix", "thresholds": {"q": 0.2}},
        {"analysis_id": "mw_reported", "test": {"kind": "mannwhitney"}, "bh_population": "reported",
         "fold_change": "geometric_mean_ratio"}],
        features={"feature_ids": ["F2", "F5", "F7"]}), None


def _case_tsv_ols(tmp_path):
    files = _tsv_files(tmp_path)
    cov = {"path": str(files["covariates"]),
           "columns": [{"name": "sex", "type": "categorical", "reference": "F"},
                       {"name": "age", "type": "numeric"}]}
    recipe = {"steps": ["normalize", "drift_correct"], "normalize": {"method": "pqn"},
              "drift_correct": {"span": 0.9, "min_qc": 4}}
    return _tsv_spec(files, [
        {"analysis_id": "ols_cov", "test": {"kind": "ols", "covariates": ["sex", "age"]},
         "matrix": recipe, "thresholds": {"p": 0.05}},
        {"analysis_id": "ols_plain", "test": {"kind": "ols", "transform": "none"}}],
        covariates=cov), None


def _case_tsv_aggregate_sum(tmp_path):
    files = _tsv_files(tmp_path)
    return _tsv_spec(files, [
        {"analysis_id": "sum", "test": {"kind": "welch"},
         "matrix": {"steps": ["impute"], "impute": {"method": "half_min"}}}],
        aggregate={"kind": "sum", "label": "F1_F2_sum", "feature_ids": ["F1", "F2"]}), None


def _case_tsv_aggregate_ratio(tmp_path):
    files = _tsv_files(tmp_path)
    return _tsv_spec(files, [
        {"analysis_id": "ratio", "test": {"kind": "student"}, "fold_change": "geometric_mean_ratio",
         "thresholds": {"p": 0.05}}],
        aggregate={"kind": "ratio", "label": "F5_over_F1", "numerator": "F5",
                   "denominator": "F1"}), None


def _dataset_spec(analyses, **extra) -> dict:
    spec = {"comparison": {"reference": "HC", "test": "RA"},
            "features": {"merge_split": False},
            "injection_order": {"source": "metadata"}, "analyses": analyses}
    spec.update(extra)
    return spec


def _case_dataset_annotated(tmp_path):
    recipe = {"steps": ["detection_filter", "impute", "normalize", "drift_correct"],
              "detection_filter": {"basis": "detected", "min_rate": {"sample": 0.5}},
              "impute": {"method": "knn", "k": 2, "scale": "none"},
              "normalize": {"method": "pqn", "reference": "sample"},
              "drift_correct": {"span": 0.9, "min_qc": 4}}
    return _dataset_spec([
        {"analysis_id": "annotated", "test": {"kind": "welch"}, "matrix": recipe,
         "bh_population": "annotated", "thresholds": {"p": 0.05, "q": 0.1, "fc": 1.2}},
        {"analysis_id": "matrix_pop", "test": {"kind": "student"}, "matrix": recipe,
         "bh_population": "matrix"},
        {"analysis_id": "ols_ds", "test": {"kind": "ols"}, "bh_population": "reported"}],
    ), _dataset()


def _case_dataset_aggregate_select(tmp_path):
    return _dataset_spec(
        [{"analysis_id": "pc_sum", "test": {"kind": "welch"}, "bh_population": "reported"}],
        aggregate={"kind": "sum", "label": "PC_PE",
                   "select": {"name_regex": "^P[CE] ", "evidence": ["sme"]}}), _dataset()


CASES = {
    "tsv_student_recipe": _case_tsv_student_recipe,
    "tsv_welch_knn_tic": _case_tsv_welch_knn_tic,
    "tsv_mannwhitney_median": _case_tsv_mannwhitney_median,
    "tsv_ols": _case_tsv_ols,
    "tsv_aggregate_sum": _case_tsv_aggregate_sum,
    "tsv_aggregate_ratio": _case_tsv_aggregate_ratio,
    "dataset_annotated": _case_dataset_annotated,
    "dataset_aggregate_select": _case_dataset_aggregate_select,
}


def _run_case(case: str, tmp_path: Path) -> dict:
    spec, ds = CASES[case](tmp_path)
    return A.run_spec(A.normalize_spec(spec), ds)


def _numbers(result: dict) -> dict:
    """数だけを取り出す。パス（input_fingerprint.files）・dataset_hash は入れない。"""
    out = {"rule_hash": result["rule_hash"], "analyses": []}
    for analysis in result["analyses"]:
        item = {"analysis_id": analysis["analysis_id"], "rows": analysis["rows"],
                "matrix": {"n_features_in": analysis["matrix"]["n_features_in"],
                           "n_features_out": analysis["matrix"]["n_features_out"],
                           "steps": analysis["matrix"]["steps"]},
                "bh_population_n": analysis["bh_population_n"]}
        for key in ("covariate_excluded_samples", "aggregate"):
            if key in analysis:
                item[key] = analysis[key]
        out["analyses"].append(item)
    return _jsonable(out)


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return None if math.isnan(value) else value
    return value


def _assert_same(actual, expected, where="") -> None:
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and list(actual) == list(expected), where
        for key in expected:
            _assert_same(actual[key], expected[key], f"{where}/{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), where
        for i, (a, e) in enumerate(zip(actual, expected)):
            _assert_same(a, e, f"{where}[{i}]")
    elif isinstance(expected, float) and not isinstance(expected, bool):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), where
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12), (where, actual, expected)
    else:
        assert actual == expected, (where, actual, expected)


@pytest.mark.parametrize("case", sorted(CASES))
def test_analysis_spec_v1_numbers_are_frozen(case, tmp_path):
    actual = _numbers(_run_case(case, tmp_path))
    if os.environ.get("UPDATE_V1_GOLDEN") == "1":
        golden = json.loads(GOLDEN.read_text(encoding="utf-8")) if GOLDEN.exists() else {}
        golden[case] = actual
        GOLDEN.parent.mkdir(parents=True, exist_ok=True)
        GOLDEN.write_text(json.dumps(dict(sorted(golden.items())), ensure_ascii=False, indent=1)
                          + "\n", encoding="utf-8")
        return
    golden = json.loads(GOLDEN.read_text(encoding="utf-8"))
    _assert_same(actual, golden[case], case)


def test_golden_has_exactly_the_cases():
    if os.environ.get("UPDATE_V1_GOLDEN") == "1":
        return
    assert sorted(json.loads(GOLDEN.read_text(encoding="utf-8"))) == sorted(CASES)


def test_golden_covers_every_v1_option(tmp_path):
    kinds, fcs, pops, steps, imputes, norms = set(), set(), set(), set(), set(), set()
    aggregates = set()
    for builder in CASES.values():
        spec, _ = builder(tmp_path)
        normalized = A.normalize_spec(spec)
        if normalized["aggregate"]:
            aggregates.add(normalized["aggregate"]["kind"])
        for analysis in normalized["analyses"]:
            kinds.add(analysis["test"]["kind"])
            fcs.add(analysis["fold_change"])
            pops.add(analysis["bh_population"])
            recipe = analysis["matrix"]
            steps.update(recipe["steps"])
            if "impute" in recipe:
                imputes.add(recipe["impute"]["method"])
            if "normalize" in recipe:
                norms.add(recipe["normalize"]["method"])
    assert kinds == set(spec_statistics.KINDS)
    assert fcs == set(spec_statistics.FC_DEFINITIONS)
    assert pops == {"matrix", "reported", "annotated"}
    assert steps == set(recipe_matrix.STEPS)
    assert {"half_min", "knn"} <= imputes and {"pqn", "tic", "median"} <= norms
    assert aggregates == {"sum", "ratio"}
