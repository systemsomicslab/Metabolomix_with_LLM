"""全feature出力と三軸レポート（spec §11）。

v2 は「終わったか」を3つの独立した軸で持つ:

- `execution_status` — pipeline が回りきったか（既存）
- `qc_status` — バッチ品質（pass/fail/not_evaluable）
- `analysis_status` — 指定した統計が計算できたか（ready/limited/not_evaluable）

QC が fail でも統計が計算できるなら診断用に続行し、結果に不合格を付記する。
逆に、統計が全部計算できなくても QC は独立に語れる。1つの「成功/失敗」へ
畳むと、どちらの理由で止まったのかが読み手から消える。

出力側で縛るのは:

- 全feature表は**注釈が無い行も捨てない**。捨てると「同定できたものだけの世界」
  になり、未同定の大半が存在しないことになる。
- 元SMEの複数候補を保持し、best hit を確定同定と呼ばない。
- 欠損・±inf を有意値として書かない。
- ANOVA の結果を既存15列 volcano TSV へ押し込まない（列の意味が変わる）。
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from metabolomix.analysis.feature_export import export_features, export_statistic
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.pipeline.report import analysis_status, required_outputs_v2


# ---------- brief記載のRED ----------

def test_computable_and_uncomputable_are_limited():
    assert analysis_status([{"status": "ready"}, {"status": "not_evaluable"}]) == "limited"


def test_all_ready_is_ready_and_none_ready_is_not_evaluable():
    assert analysis_status([{"status": "ready"}, {"status": "ready"}]) == "ready"
    assert analysis_status([{"status": "not_evaluable"}]) == "not_evaluable"
    assert analysis_status([]) == "not_evaluable"


# ---------- fixture ----------

def _dataset() -> DatasetState:
    ds = DatasetState()
    ds.feature_ids = ["1", "2"]
    ds.sample_assay_ids = ["assay[1]", "assay[2]"]
    ds.sample_names = ["S1", "S2"]
    ds.feature_metadata = {
        "1": {"name": "GABA", "mz": 104.0706, "rt": 1.2,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
        # 2 は完全に未同定（名前もInChIKeyも無い）。
        "2": {"name": None, "mz": 500.1, "rt": 5.0, "inchikey": None},
    }
    ds.feature_candidates = {
        "1": [{"sme_id": "1", "rank": 1, "chemical_name": "GABA",
               "database_identifier": "HMDB:HMDB0000112", "adduct": "[M+H]1+",
               "charge": 1, "confidence_value": 0.92,
               "confidence_measure": "MS-DIAL total score"},
              {"sme_id": "2", "rank": 2, "chemical_name": "other",
               "database_identifier": "HMDB:HMDB0000999", "adduct": "[M+Na]1+",
               "charge": 1, "confidence_value": 0.41,
               "confidence_measure": "MS-DIAL total score"}],
        "2": [],
    }
    return ds


def _matrix(matrix_id="mat_a", *, units=None, eligibility=None, imputed=None,
            values=None, detected=None, stage="finalized") -> dict:
    values = np.asarray(values if values is not None
                        else [[10., 20.], [11., 21.]], dtype=float)
    return {
        "schema": "analysis-matrix.v1", "matrix_id": matrix_id,
        "parent_id": None, "stage": stage, "recipe_id": "default",
        "dataset_id": "ds_test",
        "assay_ids": ["assay[1]", "assay[2]"], "feature_ids": ["1", "2"],
        "values": values,
        "units": units or ["peak_height", "peak_height"],
        "eligibility_mask": (np.ones(2, dtype=bool) if eligibility is None
                             else np.asarray(eligibility, dtype=bool)),
        "detected_mask": (None if detected is None
                          else np.asarray(detected, dtype=bool)),
        "imputed_mask": (np.zeros(values.shape, dtype=bool) if imputed is None
                         else np.asarray(imputed, dtype=bool)),
        "locked_mask": np.zeros(values.shape, dtype=bool),
        "missing_reasons": {}, "support_feature_ids": [],
        "correction_history": [], "caveats": [],
    }


def _read(path):
    """TSV本体だけを返す（先頭の `# key = value` メタ行は既存exportと同じ流儀）。"""
    return [line.split("\t") for line in
            path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


# ---------- 全feature表 ----------

def test_every_feature_is_exported_even_without_annotation(tmp_path):
    ds = _dataset()
    out = export_features(ds, [_matrix()], tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header, body = rows[0], rows[1:]
    assert header[:5] == ["dataset_id", "feature_id", "assay_id",
                          "matrix_result_id", "value"]
    assert len(body) == 4                     # feature 2 × assay 2
    assert {r[1] for r in body} == {"1", "2"}
    assert out["n_features"] == 2


def test_units_distinguish_raw_from_corrected(tmp_path):
    ds = _dataset()
    matrices = [_matrix("mat_raw"),
                _matrix("mat_ratio",
                        units=["internal_standard_ratio", "peak_height"])]
    export_features(ds, matrices, tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header = rows[0]
    unit_col = header.index("unit")
    matrix_col = header.index("matrix_result_id")
    units = {(r[matrix_col], r[1]): r[unit_col] for r in rows[1:]}
    assert units[("mat_raw", "1")] == "peak_height"
    assert units[("mat_ratio", "1")] == "internal_standard_ratio"


def test_detection_and_imputation_state_are_columns(tmp_path):
    ds = _dataset()
    matrix = _matrix(detected=[[True, False], [True, True]],
                     imputed=[[False, True], [False, False]])
    export_features(ds, [matrix], tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header = rows[0]
    by_key = {(r[1], r[2]): r for r in rows[1:]}
    detected_col = header.index("detected")
    gap_col = header.index("gap_filled")
    imputed_col = header.index("imputed")
    assert by_key[("2", "assay[1]")][detected_col] == "false"
    assert by_key[("2", "assay[1]")][gap_col] == "true"
    assert by_key[("2", "assay[1]")][imputed_col] == "true"


def test_unknown_detection_is_not_written_as_detected(tmp_path):
    ds = _dataset()
    export_features(ds, [_matrix(detected=None)], tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header = rows[0]
    detected_col = header.index("detected")
    assert {r[detected_col] for r in rows[1:]} == {"unknown"}


def test_excluded_features_say_why(tmp_path):
    ds = _dataset()
    matrix = _matrix(eligibility=[True, False])
    export_features(ds, [matrix], tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header = rows[0]
    reason_col = header.index("exclusion_reason")
    by_feature = {r[1]: r[reason_col] for r in rows[1:]}
    assert by_feature["1"] == ""
    assert by_feature["2"] == "not_eligible"


def test_annotations_go_to_a_side_table_keyed_by_feature_id(tmp_path):
    ds = _dataset()
    out = export_features(ds, [_matrix()], tmp_path / "features.tsv")
    annotations = _read(tmp_path / out["annotation_file"])
    header, body = annotations[0], annotations[1:]
    assert header[0] == "feature_id"
    by_feature = {}
    for row in body:
        by_feature.setdefault(row[0], []).append(row)
    # 未同定の feature も1行残す（候補0件であることを書く）。
    assert set(by_feature) == {"1", "2"}
    assert len(by_feature["1"]) == 2          # SME候補を全部残す
    assert len(by_feature["2"]) == 1


def test_the_best_hit_is_not_called_a_confirmed_identification(tmp_path):
    ds = _dataset()
    out = export_features(ds, [_matrix()], tmp_path / "features.tsv")
    annotations = _read(tmp_path / out["annotation_file"])
    header = annotations[0]
    assert "identification_status" in header
    status_col = header.index("identification_status")
    statuses = {row[status_col] for row in annotations[1:]}
    assert "confirmed" not in statuses
    assert "candidate" in statuses


def test_non_finite_values_are_written_as_empty_not_as_numbers(tmp_path):
    ds = _dataset()
    matrix = _matrix(values=[[np.nan, np.inf], [1., 2.]])
    export_features(ds, [matrix], tmp_path / "features.tsv")
    rows = _read(tmp_path / "features.tsv")
    header = rows[0]
    value_col = header.index("value")
    values = {(r[1], r[2]): r[value_col] for r in rows[1:]}
    assert values[("1", "assay[1]")] == ""
    assert values[("2", "assay[1]")] == ""
    assert "inf" not in {v.lower() for v in values.values()}


# ---------- 統計結果の出力 ----------

def _welch_result(**overrides) -> dict:
    result = {
        "schema": "statistic-result.v1", "statistic_id": "w1", "kind": "welch",
        "status": "completed", "reason": None, "matrix_id": "mat_a",
        "matrix_recipe_id": "default", "transform": "none",
        "effect_size_definition": "log2_arithmetic_mean_ratio",
        "reference_group": "control", "test_group": "treated",
        "groups": {"control": 3, "treated": 3}, "bh_population": 1,
        "features": [
            {"feature_id": "1", "status": "tested", "reason": None,
             "p_value": 0.001, "q_value": 0.001, "t_statistic": 12.0,
             "log2fc": 2.0, "n_reference": 3, "n_test": 3},
            {"feature_id": "2", "status": "not_testable", "reason": "zero_variance",
             "p_value": None, "q_value": None, "t_statistic": None,
             "log2fc": None, "n_reference": 3, "n_test": 3},
        ],
    }
    result.update(overrides)
    return result


def test_welch_export_keeps_untestable_rows_with_their_reason(tmp_path):
    out = export_statistic(_welch_result(), tmp_path / "w1.tsv")
    rows = _read(tmp_path / "w1.tsv")
    header, body = rows[0], rows[1:]
    assert len(body) == 2
    reason_col = header.index("reason")
    status_col = header.index("status")
    by_feature = {r[header.index("feature_id")]: r for r in body}
    assert by_feature["2"][status_col] == "not_testable"
    assert by_feature["2"][reason_col] == "zero_variance"
    assert out["n_tested"] == 1


def test_missing_p_values_are_blank_not_zero(tmp_path):
    export_statistic(_welch_result(), tmp_path / "w1.tsv")
    rows = _read(tmp_path / "w1.tsv")
    header = rows[0]
    p_col = header.index("p_value")
    by_feature = {r[header.index("feature_id")]: r for r in rows[1:]}
    assert by_feature["2"][p_col] == ""


def test_export_preserves_small_nonzero_p_and_q_values(tmp_path):
    result = _welch_result()
    result["features"][0].update(p_value=1e-9, q_value=2.123456789e-8)
    export_statistic(result, tmp_path / "small.tsv")
    header, row, *_ = _read(tmp_path / "small.tsv")
    assert float(row[header.index("p_value")]) == 1e-9
    assert float(row[header.index("q_value")]) == 2.123456789e-8


def test_effect_size_definition_travels_with_the_table(tmp_path):
    out = export_statistic(_welch_result(), tmp_path / "w1.tsv")
    text = (tmp_path / "w1.tsv").read_text(encoding="utf-8")
    assert "log2_arithmetic_mean_ratio" in text
    assert out["effect_size_definition"] == "log2_arithmetic_mean_ratio"


def test_infinite_effect_size_is_not_written_as_a_number(tmp_path):
    result = _welch_result()
    result["features"][0]["log2fc"] = math.inf
    export_statistic(result, tmp_path / "w1.tsv")
    rows = _read(tmp_path / "w1.tsv")
    header = rows[0]
    fc_col = header.index("log2fc")
    assert rows[1][fc_col] == ""


# ---------- ANOVA 専用出力 ----------

def _anova_result() -> dict:
    return {
        "schema": "statistic-result.v1", "statistic_id": "a1",
        "kind": "anova_tukey", "status": "completed", "reason": None,
        "matrix_id": "mat_a", "matrix_recipe_id": "default", "transform": "none",
        "groups": {"g1": 3, "g2": 3, "g3": 3}, "alpha": 0.05,
        "bh_population": 1,
        "tukey_correction_scope": "within_feature_group_pairs",
        "features": [
            {"feature_id": "1", "status": "tested", "reason": None,
             "f_statistic": 7.0, "p_value": 0.027, "q_value": 0.027,
             "df_between": 2, "df_within": 6,
             "group_n": {"g1": 3, "g2": 3, "g3": 3},
             "tukey": [{"reference_group": "g1", "test_group": "g3",
                        "mean_difference": 3.0, "ci_low": 0.49, "ci_high": 5.51,
                        "p_adjusted": 0.0242}]},
        ],
    }


def test_anova_writes_its_own_table_not_the_fifteen_column_contract(tmp_path):
    out = export_statistic(_anova_result(), tmp_path / "a1.tsv")
    header = _read(tmp_path / "a1.tsv")[0]
    assert "f_statistic" in header and "df_between" in header
    assert "log2fc" not in header          # volcano契約の列を流用しない
    assert out["kind"] == "anova_tukey"


def test_tukey_pairs_go_to_a_separate_table(tmp_path):
    out = export_statistic(_anova_result(), tmp_path / "a1.tsv")
    pairs = _read(tmp_path / out["tukey_file"])
    header, body = pairs[0], pairs[1:]
    assert header[:4] == ["feature_id", "reference_group", "test_group",
                          "mean_difference"]
    assert len(body) == 1
    assert "p_adjusted" in header


def test_tukey_table_states_that_it_is_not_a_cross_feature_fdr(tmp_path):
    out = export_statistic(_anova_result(), tmp_path / "a1.tsv")
    text = (tmp_path / out["tukey_file"]).read_text(encoding="utf-8")
    assert "within_feature_group_pairs" in text


def test_a_not_evaluable_statistic_still_writes_its_reason(tmp_path):
    result = _welch_result(status="not_evaluable", reason="no_testable_feature",
                           features=[])
    out = export_statistic(result, tmp_path / "w1.tsv")
    text = (tmp_path / "w1.tsv").read_text(encoding="utf-8")
    assert "no_testable_feature" in text
    assert out["status"] == "not_evaluable"


# ---------- v2 必須成果物 ----------

def _v2_request(**overrides) -> dict:
    request = {
        "schema": "pipeline-request.v2",
        "statistics": [{"statistic_id": "w1", "kind": "welch"},
                       {"statistic_id": "p1", "kind": "pca"}],
        "standard_assays": {"gaba": ["s1"]},
    }
    request.update(overrides)
    return request


def test_v2_required_outputs_cover_the_spec_list():
    names = required_outputs_v2(_v2_request())
    for expected in ("profile", "execution_manifest", "sample_manifest",
                     "assay_evidence", "feature_bindings", "matrix",
                     "qc_population", "qc", "feature_table", "quality_report"):
        assert expected in names


def test_every_statistic_gets_its_own_required_output():
    names = required_outputs_v2(_v2_request())
    assert "statistic:w1" in names
    assert "statistic:p1" in names


def test_without_internal_standards_binding_is_still_required():
    """内部標準未使用でも binding の not_applicable を記録する（spec §11）。"""
    names = required_outputs_v2(_v2_request(standard_assays={}))
    assert "feature_bindings" in names


def test_pathway_tsv_is_not_a_required_v2_output():
    """既存15列TSVは任意。対象ゼロでpipelineを未完了にしない。"""
    names = required_outputs_v2(_v2_request())
    assert not any(name.startswith("pathway") for name in names)


def test_v1_required_outputs_are_unchanged():
    from metabolomix.pipeline.report import required_outputs

    v1 = required_outputs({"save_project": False, "effective_target": "exploratory",
                           "comparisons": []})
    assert v1 == ["preprocess", "pca", "pca_figure", "quality_report"]


# ---------- 15列TSVへの任意出力 ----------

def test_v1_export_meta_is_unchanged_without_an_effect_size_definition():
    """v1の結果にはこの欄が無い。無い結果へ行を足さない。"""
    from metabolomix.analysis.dataset_export import _effect_size_lines

    assert _effect_size_lines({"provenance": {}}) == []


def test_v2_effect_size_definition_is_inherited_into_the_meta():
    from metabolomix.analysis.dataset_export import _effect_size_lines

    lines = _effect_size_lines(
        {"effect_size_definition": "log2_arithmetic_mean_ratio"})
    assert lines == ["# effect_size_definition = log2_arithmetic_mean_ratio"]



# ---------- SML 由来の注釈（spec 2026-09-17-mztab-sml-annotation-design §7.1）----------

def _ms1_annotation(**over) -> dict:
    annotation = {"sml_id": "1", "ambiguous": False, "name": "GABA",
                  "database_identifier": "TextDB:GABA", "chemical_formula": None,
                  "smiles": None, "adduct": "[M+H]1+", "reliability": None,
                  "confidence_measure": "MS-DIAL algorithm matching score",
                  "confidence_value": 0.99,
                  "inchikey": None, "inchikey_source": "none"}
    annotation.update(over)
    return annotation


def _annotation_rows(tmp_path, ds) -> list[dict]:
    """`.annotations.tsv` を dict の一覧で返す。"""
    out = export_features(ds, [_matrix()], tmp_path / "features.tsv")
    table = _read(tmp_path / out["annotation_file"])
    header = table[0]
    return [dict(zip(header, row)) for row in table[1:]]


def test_an_ms1_annotation_gets_its_own_identification_status(tmp_path):
    """SME 候補が無く SML 注釈だけがある特徴は `ms1_annotation`。

    `unidentified` に寄せると「同定できなかった」と読まれ、`candidate` に
    寄せると MS/MS 裏付けがあるように読まれる。
    """
    ds = _dataset()
    ds.feature_annotations = {"2": _ms1_annotation()}

    rows = _annotation_rows(tmp_path, ds)

    row = next(r for r in rows if r["feature_id"] == "2")
    assert row["identification_status"] == "ms1_annotation"
    assert row["name"] == "GABA"
    assert row["database_identifier"] == "TextDB:GABA"
    assert row["adduct"] == "[M+H]1+"
    assert row["candidate_rank"] == ""


def test_an_sme_backed_feature_is_unaffected_by_an_annotation(tmp_path):
    """SME 候補がある特徴は従来どおり candidate。注釈があっても変わらない。"""
    ds = _dataset()
    ds.feature_annotations = {"1": _ms1_annotation(name="別の名前"),
                              "2": _ms1_annotation()}

    rows = _annotation_rows(tmp_path, ds)

    for row in (r for r in rows if r["feature_id"] == "1"):
        assert row["identification_status"] == "candidate"
        assert row["name"] != "別の名前"


def test_a_feature_without_any_annotation_stays_unidentified(tmp_path):
    ds = _dataset()
    ds.feature_annotations = {}

    rows = _annotation_rows(tmp_path, ds)

    row = next(r for r in rows if r["feature_id"] == "2")
    assert row["identification_status"] == "unidentified"


def test_an_ambiguous_annotation_stays_unidentified(tmp_path):
    """どの名前かを決められない注釈で名前を出すと嘘になる。"""
    ds = _dataset()
    ds.feature_annotations = {
        "2": {"ambiguous": True, "sml_ids": ["1", "2"], "name": None}}

    rows = _annotation_rows(tmp_path, ds)

    row = next(r for r in rows if r["feature_id"] == "2")
    assert row["identification_status"] == "unidentified"
