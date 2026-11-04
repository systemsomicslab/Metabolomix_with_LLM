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
    assert normalized["features"] == {"all": True, "merge_split": True}
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


# ---------- 割れたスポットのまとめ（merge_split）・BH の母集団 annotated ----------

def _split_ds(**kwargs):
    """CE 18:2 が最初の 2 注入だけ別スポット（1266）に割れている（injection_block）+ 割れない 2000。"""
    from tests.test_alignment_qc import ADDUCT, ALL, FIRST_TWO, _make_ds

    return _make_ds([("1265", 650.6000, 20.00, ALL - FIRST_TWO, "CE 18:2", ADDUCT),
                     ("1266", 650.6050, 20.61, FIRST_TWO, None, ADDUCT),
                     ("2000", 900.5000, 12.00, ALL, "PC 34:1", ADDUCT)], **kwargs)


def _split_spec(**extra):
    spec = {"comparison": {"reference": "A", "test": "B"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]}
    spec.update(extra)
    return spec


def _ds_rows(ds, spec, analysis=0):
    result = A.run_spec(A.normalize_spec(spec), ds)
    return result, {r["feature_id"]: r for r in result["analyses"][analysis]["rows"]}


def test_merge_split_takes_per_sample_max():
    ds = _split_ds()
    merged = np.maximum(ds.feature_matrix[0], ds.feature_matrix[1])     # 試料ごとに大きい方
    result, rows = _ds_rows(ds, _split_spec(features={"all": True, "merge_split": True}))
    assert sorted(rows) == ["1265", "2000"]                    # 最も小さい ID を残す
    assert rows["1265"]["mean_ref"] == pytest.approx(merged[:21].mean())
    assert rows["1265"]["mean_test"] == pytest.approx(merged[21:].mean())
    assert rows["1265"]["p"] == pytest.approx(
        stats.ttest_ind(np.log2(merged[21:]), np.log2(merged[:21]), equal_var=False).pvalue)
    assert "merge_split" not in result["analyses"][0] and result["merge_split"]["skipped"] == []
    assert result["merge_split"]["enabled"] and result["merge_split"]["source"] == "find_split_pairs"
    assert result["merge_split"]["groups"] == [["1265", "1266"]]


def test_merge_split_false_keeps_both():
    ds = _split_ds()
    result, rows = _ds_rows(ds, _split_spec(features={"all": True, "merge_split": False}))
    assert sorted(rows) == ["1265", "1266", "2000"]
    assert rows["1265"]["mean_ref"] == pytest.approx(ds.feature_matrix[0][:21].mean())
    assert all(r["merged_from"] is None for r in rows.values())
    assert result["merge_split"]["enabled"] is False and result["merge_split"]["groups"] == []


def test_merge_split_default_is_true():
    assert A.normalize_spec(_split_spec())["features"] == {"all": True, "merge_split": True}
    assert A.normalize_spec(_split_spec(features={"feature_ids": ["1265"]}))["features"] == {
        "feature_ids": ["1265"], "merge_split": True}
    _, rows = _ds_rows(_split_ds(), _split_spec())                  # features を書かなくても
    assert sorted(rows) == ["1265", "2000"]
    with pytest.raises(DomainError):                                # 型の違う値・未知のキーは拒否
        A.normalize_spec(_split_spec(features={"all": True, "merge_split": "yes"}))
    with pytest.raises(DomainError):
        A.normalize_spec(_split_spec(features={"all": True, "merge": True}))
    with pytest.raises(DomainError):                                # 効かない設定は書かせない
        A.normalize_spec(_split_spec(features={"all": True, "merge_split": False,
                                               "alignment_qc": "qc.json"}))


def test_merged_from_in_rows():
    ds = _split_ds()
    _, rows = _ds_rows(ds, _split_spec())
    assert rows["1265"]["merged_from"] == ["1265", "1266"]
    assert rows["2000"]["merged_from"] is None
    # 消えた側の ID で報告を頼んでも、残した 1 行（1 回だけ）に当たる
    for ids in (["1266"], ["1266", "1265"]):
        _, rows = _ds_rows(ds, _split_spec(features={"feature_ids": ids}))
        assert list(rows) == ["1265"] and rows["1265"]["merged_from"] == ["1265", "1266"]
    # 行の TSV に merged_from の列が出る
    from metabolomix.tools import reanalysis_tools as R

    table = R.rows_table(list(rows.values())).splitlines()
    header, first = table[0].split("	"), table[1].split("	")
    assert first[header.index("merged_from")] == "1265;1266"


def test_aggregate_counts_split_group_once():
    # C6 の型: 和に割れた 2 スポット（1185/1186）の両方が選ばれても 1 回だけ数える
    from tests.test_alignment_qc import ADDUCT, ALL, FIRST_TWO, _make_ds

    ds = _make_ds([("1185", 650.6000, 20.00, ALL - FIRST_TWO, "CE 18:2", ADDUCT),
                   ("1186", 650.6050, 20.61, FIRST_TWO, "CE 18:2", ADDUCT),
                   ("1187", 700.5000, 22.00, ALL, "CE 20:4", ADDUCT)])
    merged = np.maximum(ds.feature_matrix[0], ds.feature_matrix[1])
    by_name = {"kind": "sum", "label": "CE", "select": {"name_regex": "^CE ", "evidence": ["sme"]}}
    by_ids = {"kind": "sum", "label": "CE", "feature_ids": ["1185", "1186", "1187"]}
    for aggregate in (by_name, by_ids):
        result, _ = _ds_rows(ds, _split_spec(aggregate=aggregate))
        info = result["analyses"][0]["aggregate"]
        assert info["summed"] == ["1185", "1187"]
        assert info["merged_from"] == {"1185": ["1185", "1186"]}
        assert info["members"] == ["1185", "1186", "1187"] and info["dropped"] == []
        row = result["analyses"][0]["rows"][0]
        assert row["mean_ref"] == pytest.approx((merged[:21] + ds.feature_matrix[2][:21]).mean())
        off = A.run_spec(A.normalize_spec(_split_spec(
            aggregate=aggregate, features={"all": True, "merge_split": False})), ds)
        assert off["analyses"][0]["aggregate"]["summed"] == ["1185", "1186", "1187"]
        assert off["analyses"][0]["rows"][0]["mean_ref"] == pytest.approx(
            ds.feature_matrix[:, :21].sum(axis=0).mean())
    # 比の分子・分母が同じ物質にまとまるなら、比にならない
    with pytest.raises(DomainError):
        _ds_rows(ds, _split_spec(aggregate={"kind": "ratio", "label": "r", "numerator": "1185",
                                            "denominator": "1186"}))


def test_only_alignment_qc_merge_groups_are_merged():
    # 群ごとに別の異性体の可能性がある group_aligned は、まとめない
    from metabolomix.analysis import alignment_qc
    from tests.test_alignment_qc import ADDUCT, ALL, GROUPS_ALTERNATE, _make_ds

    owned = {i for i, g in enumerate(GROUPS_ALTERNATE) if g == "B"}     # 群が注入順で交互
    ds = _make_ds([("1", 500.3000, 8.0, ALL - owned, None, ADDUCT),
                   ("2", 500.3010, 8.3, owned, None, ADDUCT)], groups=GROUPS_ALTERNATE)
    (pair,) = alignment_qc.find_split_pairs(ds)["pairs"]
    assert pair["type"] == "group_aligned" and not pair["mergeable"]
    result, rows = _ds_rows(ds, _split_spec())
    assert sorted(rows) == ["1", "2"] and result["merge_split"]["groups"] == []
    assert all(r["merged_from"] is None for r in rows.values())


def test_alignment_qc_json_is_used_and_fingerprinted(tmp_path):
    from metabolomix.analysis import alignment_qc

    ds = _split_ds()
    path = tmp_path / "alignment_qc_pos.json"
    path.write_text(json.dumps(alignment_qc.find_split_pairs(ds)), encoding="utf-8")
    spec = _split_spec(features={"all": True, "alignment_qc": str(path)})
    result, rows = _ds_rows(ds, spec)
    assert rows["1265"]["merged_from"] == ["1265", "1266"]
    assert result["merge_split"]["source"] == "alignment_qc"
    assert len(result["input_fingerprint"]["files"]["alignment_qc"]) == 64
    live = A.run_spec(A.normalize_spec(_split_spec()), ds)
    assert result["merge_split"]["groups_hash"] == live["merge_split"]["groups_hash"]
    assert result["rule_hash"] == live["rule_hash"]               # 置き場所・出どころでは変わらない
    qc = alignment_qc.find_split_pairs(ds)
    # 空の merge_groups の JSON ならまとめない（まとめの判断は JSON の側）
    path.write_text(json.dumps({**qc, "merge_groups": []}), encoding="utf-8")
    _, rows = _ds_rows(ds, spec)
    assert sorted(rows) == ["1265", "1266", "2000"]
    # 行列に 2 つ揃わない群はまとめず、skipped に挙げる
    extra = {"merge_group_id": "M9", "features": ["1265", "9999"], "mergeable": True}
    path.write_text(json.dumps({**qc, "merge_groups": [extra]}), encoding="utf-8")
    result, rows = _ds_rows(ds, spec)
    assert sorted(rows) == ["1265", "1266", "2000"]
    assert result["merge_split"]["groups"] == [] and result["merge_split"]["skipped"] == [
        ["1265", "9999"]]
    path.write_text("{}", encoding="utf-8")
    with pytest.raises(DomainError) as info:
        _ds_rows(ds, spec)
    assert info.value.code == "ANALYSIS_INPUT_INVALID"


def test_alignment_qc_json_from_another_dataset_is_rejected(tmp_path):
    # 特徴 ID は小さい整数で別の run にも同じ ID がある。指紋の違う JSON は黙ってまとめず拒む
    from metabolomix.analysis import alignment_qc

    ds = _split_ds()
    qc = alignment_qc.find_split_pairs(ds)
    assert isinstance(qc["dataset_fingerprint"], str) and len(qc["dataset_fingerprint"]) == 64
    spec = _split_spec(features={"all": True, "alignment_qc": str(tmp_path / "qc.json")})
    other = _split_ds()
    other.feature_matrix = other.feature_matrix * 2.0                # 同じ ID・別の値 = 別の run
    assert alignment_qc.find_split_pairs(other)["dataset_fingerprint"] != qc["dataset_fingerprint"]
    for bad in ({**qc, "dataset_fingerprint": alignment_qc.find_split_pairs(other)["dataset_fingerprint"]},
                {k: v for k, v in qc.items() if k != "dataset_fingerprint"}):
        (tmp_path / "qc.json").write_text(json.dumps(bad), encoding="utf-8")
        with pytest.raises(DomainError) as info:
            _ds_rows(ds, spec)
        assert info.value.code == "ANALYSIS_INPUT_INVALID"
        assert info.value.details["field"] == "features.alignment_qc"


def test_merge_happens_before_the_recipe_so_a_detection_filter_keeps_the_group():
    # 割れの少数側（1266）は検出が 2 試料だけで、単独なら検出率の篩に落ちる。recipe の前にまとめ、
    # まとめた列の検出状態を「どちらかで検出」にするので、群は篩を通る
    ds = _split_ds(detected=True)
    recipe = {"steps": ["detection_filter"],
              "detection_filter": {"basis": "detected", "min_rate": {"sample": 0.5}}}
    analyses = [{"analysis_id": "main", "test": {"kind": "welch"}, "matrix": recipe}]
    result, rows = _ds_rows(ds, _split_spec(analyses=analyses))
    assert sorted(rows) == ["1265", "2000"] and rows["1265"]["merged_from"] == ["1265", "1266"]
    assert rows["1265"]["status"] == "tested"
    assert result["analyses"][0]["matrix"]["n_features_in"] == 2           # まとめた後の特徴数
    assert result["analyses"][0]["matrix"]["steps"][0]["removed"] == []
    off, _ = _ds_rows(ds, _split_spec(analyses=analyses, features={"all": True,
                                                                  "merge_split": False}))
    assert off["analyses"][0]["matrix"]["steps"][0]["removed"] == ["1266"]


def test_rule_hash_changes_with_merge_split():
    ds = _split_ds()
    on = A.run_spec(A.normalize_spec(_split_spec()), ds)
    off = A.run_spec(A.normalize_spec(_split_spec(features={"all": True, "merge_split": False})), ds)
    assert on["rule_hash"] != off["rule_hash"]
    normalized = A.normalize_spec(_split_spec())
    # 使った群の一覧が変われば（同じ spec でも）別の規則
    hashes = {A.rule_hash(normalized, [["1", "2"]]), A.rule_hash(normalized, [["1", "3"]]),
              A.rule_hash(normalized, [])}
    assert len(hashes) == 3
    assert A.rule_hash(normalized, [["2", "1"]]) == A.rule_hash(normalized, [["1", "2"]])
    assert on["merge_split"]["groups_hash"] == A.merge_groups_hash([["1266", "1265"]])


def test_merge_split_is_a_no_op_for_tsv_input(files):
    # 外部の行列にはアライメントが無い（較正用）。merge_split は何もせず、hash にも入らない
    on = _run(_spec(files, [MAIN], features={"all": True, "merge_split": True}))
    off = _run(_spec(files, [MAIN], features={"all": True, "merge_split": False}))
    assert on["merge_split"]["enabled"] is False and on["merge_split"]["groups"] == []
    assert on["rule_hash"] == off["rule_hash"]
    assert on["analyses"] == off["analyses"] and "merge_split" not in on["analyses"][0]
    assert all(r["merged_from"] is None for r in on["analyses"][0]["rows"])
    with pytest.raises(DomainError):
        A.normalize_spec(_spec(files, [MAIN], features={"all": True, "alignment_qc": "x.json"}))


def _annotated_ds():
    from tests.test_alignment_qc import ADDUCT, ALL, _make_ds

    names = ["PC 34:1", None, "PC 36:2", None, "TG 52:3"]
    ds = _make_ds([(str(k), 300.0 + 100 * k, 2.0 + k, ALL, name, ADDUCT)
                   for k, name in enumerate(names, start=1)])
    rng = np.random.default_rng(7)
    shift = np.outer([0.0, 0.1, 0.5, 0.0, 0.9], np.arange(42) >= 21)
    ds.feature_matrix = np.exp(rng.normal(8.0, 0.5, size=(5, 42)) + shift)
    return ds


def test_bh_population_annotated():
    from metabolomix.analysis import spec_statistics

    ds = _annotated_ds()
    analyses = [{"analysis_id": "all", "test": {"kind": "welch"}},
                {"analysis_id": "ann", "test": {"kind": "welch"}, "bh_population": "annotated"}]
    result = A.run_spec(A.normalize_spec(_split_spec(analyses=analyses)), ds)
    whole = {r["feature_id"]: r for r in result["analyses"][0]["rows"]}
    ann = {r["feature_id"]: r for r in result["analyses"][1]["rows"]}
    assert result["analyses"][0]["bh_population_n"] == 5
    assert result["analyses"][1]["bh_population_n"] == 3
    assert ann["2"]["q"] is None and ann["4"]["q"] is None        # SME の無い特徴は母集団の外
    expected = spec_statistics.bh([whole[f]["p"] for f in ("1", "3", "5")])
    assert [ann[f]["q"] for f in ("1", "3", "5")] == pytest.approx(expected)
    assert ann["1"]["p"] == whole["1"]["p"]                       # p は同じ。q だけが母集団で変わる
    assert expected != pytest.approx([whole[f]["q"] for f in ("1", "3", "5")])


def test_bh_population_annotated_counts_a_merged_group_by_either_member():
    # 1265（SME あり）と 1266（SME なし）をまとめた 1 本は annotated、割れない 2000（SME あり）も
    ds = _split_ds()
    analyses = [{"analysis_id": "ann", "test": {"kind": "welch"}, "bh_population": "annotated"}]
    result = A.run_spec(A.normalize_spec(_split_spec(analyses=analyses)), ds)
    assert result["analyses"][0]["bh_population_n"] == 2
    # 残る ID（小さい方）が SME なしでも、相手に SME があれば母集団に入る
    ds.feature_candidates["1265"] = []
    ds.feature_candidates["1266"] = [{"rank": 1, "chemical_name": "CE 18:2", "adduct": "[M+H]1+"}]
    result = A.run_spec(A.normalize_spec(_split_spec(analyses=analyses)), ds)
    assert result["analyses"][0]["bh_population_n"] == 2
    # どちらにも SME が無ければ入らない
    ds.feature_candidates["1266"] = []
    result = A.run_spec(A.normalize_spec(_split_spec(analyses=analyses)), ds)
    assert result["analyses"][0]["bh_population_n"] == 1


def test_bh_population_annotated_needs_a_dataset(files):
    with pytest.raises(DomainError):
        A.normalize_spec(_spec(files, [{**MAIN, "bh_population": "annotated"}]))
