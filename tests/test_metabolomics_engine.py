"""v2 stage を worker へ接続する（spec §6.2・§8〜§11）。

ここで縛るのは、工程を繋ぐときに一番壊れやすい4点:

1. **どの schema の要求かで handler が変わる。** `preprocess` / `export` / `report`
   は v1 と v2 で同じ stage 名だが中身が違う。1つの辞書へ後勝ちで入れると、
   片方の pipeline が黙ってもう片方の計算を走らせる。
2. **needs_input はそのまま止める。** binding が未解決のまま自動で先へ進まないし、
   同じ未解決条件で自動再試行もしない（再試行しても同じ結果で、記録だけが増える）。
3. **recipe ごとに別の行列。** 2つの recipe が同じスロットを奪い合わない。
4. **復元は hash 照合の上で。** 再開時に再計算した値で古い result ID を上書きしない。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from metabolomix.pipeline import metabolomics_handlers, stage_plan
from metabolomix.pipeline.service import build_handlers


# ---------- brief記載のRED ----------

def test_binding_change_does_not_restart_console():
    dirty = stage_plan.invalidated_v2({"feature_bindings"},
                                      {"statistics": [{"statistic_id": "s"}]})
    assert "resolve_feature_bindings" in dirty
    assert "statistics:s" in dirty
    assert "execute_console" not in dirty


def test_standard_assays_change_starts_at_binding_not_at_console():
    dirty = stage_plan.invalidated_v2({"standard_assays"},
                                      {"statistics": [{"statistic_id": "s"}]})
    assert "resolve_feature_bindings" in dirty
    assert "execute_console" not in dirty
    assert "prepare_inputs" not in dirty


def test_adding_or_removing_one_statistic_touches_only_that_statistic():
    request = {"statistics": [{"statistic_id": "a"}, {"statistic_id": "b"}]}
    dirty = stage_plan.invalidated_v2({"statistics:b"}, request)
    assert {"statistics:b", "export:b", "report"} <= dirty
    assert "statistics:a" not in dirty
    assert "preprocess" not in dirty


def test_recipe_change_starts_at_preprocess():
    dirty = stage_plan.invalidated_v2({"preprocess"},
                                      {"statistics": [{"statistic_id": "s"}]})
    assert "preprocess" in dirty and "qc_processed" in dirty
    assert "qc_raw" not in dirty
    assert "resolve_feature_bindings" not in dirty


def test_metadata_change_that_moves_the_qc_population_also_redoes_qc_raw():
    dirty = stage_plan.invalidated_v2({"qc_raw_scope"},
                                      {"statistics": [{"statistic_id": "s"}]})
    assert "qc_raw" in dirty and "preprocess" in dirty


# ---------- handler の網羅と dispatch ----------

def _v2_request(**overrides) -> dict:
    request = {"schema": stage_plan.REQUEST_SCHEMA_V2,
               "statistics": [{"statistic_id": "w1", "kind": "welch"}]}
    request.update(overrides)
    return request


def test_every_v2_stage_has_a_handler():
    handlers = build_handlers()
    plan = stage_plan.build_v2(_v2_request())
    missing = [entry["handler"] for entry in plan
               if entry["handler"] not in handlers]
    assert missing == []


def test_v1_handler_keys_are_still_present():
    handlers = build_handlers()
    for key in ("prepare_input", "upstream", "pca", "resolve_comparisons",
                "differential"):
        assert key in handlers


def test_shared_stage_names_dispatch_on_the_request_schema():
    """`preprocess` は v1 と v2 で別の計算。要求の schema で選ぶ。"""
    handlers = build_handlers()
    seen = []

    def spy(name):
        def handler(context):
            seen.append(name)
            return {"status": "succeeded", "result_refs": [], "warnings": [],
                    "error": None}
        return handler

    dispatch = metabolomics_handlers.by_schema(spy("v1"), spy("v2"))
    dispatch({"request": {"schema": "pipeline-request.v1"}})
    dispatch({"request": _v2_request()})
    assert seen == ["v1", "v2"]
    assert "preprocess" in handlers


# ---------- 合成 context で handler の挙動を見る ----------

def _profile(tmp_path: Path) -> Path:
    """draft の最小 profile（`execution_purpose="validation"` で読める）。"""
    from metabolomix.console.profile_schema import validate_profile

    profile = {
        "schema": "lcms-profile.v1", "profile_id": "p1", "revision": 1,
        "omics": "metabolomics",
        "acquisition": {
            "separation": "lc", "acquisition_type": "dda", "polarity": "positive",
            "instrument": None,
            "lc": {"column": None, "mobile_phase_a": None, "mobile_phase_b": None,
                   "flow_rate_ul_min": None, "column_temperature_c": None,
                   "gradient_profile": None},
            "ms_range": {"ms1_low_mz": 50.0, "ms1_high_mz": 1500.0,
                         "ms2_low_mz": None, "ms2_high_mz": None},
            "sample_matrix": None, "scope": "synthetic test profile"},
        "software": {"msdial_version": "5.x", "executable_path": "C:/x/MSDIALCUI.exe",
                     "executable_sha256": "a" * 64, "adapter_version": "1.0.0"},
        "processing": {
            "method_path": "method/param.txt", "method_sha256": "b" * 64,
            "measure": "peak_height",
            "dependencies": [{"dependency_id": "msp-lib-1", "kind": "msp",
                              "method_key": "Msp file path", "path": "lib/l.msp",
                              "sha256": "c" * 64, "required": True}],
            "effective_settings": {}},
        "analysis_recipe": {
            "statistics": [], "internal_standards": []},
        "feature_targets": {
            "gaba": {"compound_identifiers": {
                        "name": "GABA",
                        "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
                     "adduct": "[M+H]+", "charge": 1, "expected_mz": 104.0706,
                     "mz_tolerance_ppm": 10.0, "expected_rt_min": 1.2,
                     "rt_tolerance_min": 0.1,
                     "required_evidence": [{"kind": "mass_rt"}]}},
        "matrix_recipes": {
            "default": {"base": "peak_height", "normalize": "none",
                        "drift_correct": False, "filter": None, "impute": "none"},
            "ratio": {"base": "internal_standard_ratio", "normalize": "none",
                      "drift_correct": False, "filter": None, "impute": "none"}},
        "qc_policy": {
            "detection_rate": {"scope": "all_features", "target_ids": None,
                               "required": True,
                               "threshold": {"operator": ">=", "value": 0.5},
                               "evidence_requirement": False,
                               "minimum_pass_fraction": None}},
        # null の測定条件には理由付き evidence が要る（profile_schema の契約）。
        "evidence": {path: {"value": None, "reason": "合成テスト用のため未記載",
                            "tier": "proposed", "source_uri": None,
                            "source_hash": None, "location": None}
                     for path in (
                         "acquisition.instrument", "acquisition.sample_matrix",
                         "acquisition.lc.column", "acquisition.lc.mobile_phase_a",
                         "acquisition.lc.mobile_phase_b",
                         "acquisition.lc.flow_rate_ul_min",
                         "acquisition.lc.column_temperature_c",
                         "acquisition.lc.gradient_profile",
                         "acquisition.ms_range.ms2_low_mz",
                         "acquisition.ms_range.ms2_high_mz")},
        "validation": {"status": "draft", "scope": "synthetic test profile",
                       "certificate_path": None,
                       "certificate_sha256": None},
    }
    validate_profile(profile)
    path = tmp_path / "profile.json"
    path.write_text(json.dumps(profile), encoding="utf-8")
    return path


def _dataset():
    from metabolomix.mztab.dataset_state import DatasetState

    ds = DatasetState()
    ds.feature_matrix = np.array([[10., 12., 40., 44.],
                                  [2., 2., 2., 2.]])
    ds.feature_ids = ["1", "2"]
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(4)]
    ds.sample_names = [f"S{i + 1}" for i in range(4)]
    ds.feature_metadata = {"1": {"mz": 104.0706, "rt": 1.2,
                                 "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                                 "name": "GABA"},
                           "2": {"mz": 110.11, "rt": 1.2, "inchikey": None,
                                 "name": None}}
    ds.feature_candidates = {"1": [], "2": []}
    return ds


def _metadata():
    return [{"sample_id": f"S{i + 1}", "source_file": f"S{i + 1}.wiff",
             "role": "sample", "group": "control" if i < 2 else "treated",
             "batch": "B1", "injection_order": i + 1, "qc_pool": None,
             "include": True, "biological_sample_id": f"bio{i}"}
            for i in range(4)]


def _context(tmp_path: Path, **overrides) -> dict:
    ds = _dataset()
    context = {
        "pipeline_root": tmp_path,
        "pipeline_id": "pipe-1",
        "identity": {"pipeline_root": str(tmp_path), "source_root": str(tmp_path)},
        "stage_id": "preprocess",
        "statistic_id": None,
        "comparison_id": None,
        "attempt": 0,
        "request": _v2_request(profile_file=str(_profile(tmp_path)),
                               execution_purpose="validation",
                               preprocess={}, standard_assays={},
                               feature_bindings={}),
        "request_meta": {"revision": 1},
        "inputs": {},
        "upstream": {},
        "results": [],
        "runtime": {"dataset": ds, "metadata": _metadata()},
    }
    context.update(overrides)
    return context


def _evidence(ds, *, availability=True):
    rows = [{"feature_id": fid, "assay_id": assay_id, "observed_rt_min": 1.2,
             "observed_mz": 104.0706, "detection_status": "detected"}
            for fid in ds.feature_ids for assay_id in ds.sample_assay_ids]
    return {"schema": "assay-feature-evidence.v1", "availability": availability,
            "reasons": [], "rows": rows if availability else [],
            "source_artifact_hash": "d" * 64}


# ---------- 各handler ----------

def test_unavailable_evidence_still_persists_a_reason(tmp_path):
    """証拠が取れなくても理由JSONは必須成果物（spec §11）。黙って欠かさない。"""
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["dataset"].artifact_paths = {}
    outcome = handlers["load_assay_evidence"](context)

    assert outcome["status"] == "succeeded"
    assert [r["output_name"] for r in outcome["result_refs"]] == ["assay_evidence"]
    assert context["runtime"]["assay_evidence"]["availability"] is False


def test_unresolved_binding_stops_and_does_not_retry_itself(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    # feature 1 は m/z も RT も合うが、feature_targets は GABA を1件要求する。
    context["runtime"]["dataset"].feature_metadata["2"]["mz"] = 104.0706
    context["runtime"]["dataset"].feature_metadata["2"]["inchikey"] = \
        "BTCSSZJGUNDROE-UHFFFAOYSA-N"

    first = handlers["resolve_feature_bindings"](context)
    assert first["status"] == "needs_input"
    assert first["error"]["code"] == "FEATURE_BINDING_UNRESOLVED"

    second = handlers["resolve_feature_bindings"](context)
    assert second["status"] == "needs_input"
    # 同じ未解決条件では結果も同じ。新しい成果物を積み増さない。
    assert second["error"]["details"]["unresolved"] == \
        first["error"]["details"]["unresolved"]


def test_resolved_binding_records_the_result(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    outcome = handlers["resolve_feature_bindings"](context)
    assert outcome["status"] == "succeeded"
    assert [r["output_name"] for r in outcome["result_refs"]] == ["feature_bindings"]
    assert context["runtime"]["bindings"]["status"] == "resolved"


def test_two_recipes_produce_two_matrices(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    handlers["resolve_feature_bindings"](context)
    handlers["qc_raw"](context)
    outcome = handlers["preprocess"](context)

    assert outcome["status"] == "succeeded"
    matrices = context["runtime"]["matrices"]
    assert set(matrices) == {"default", "ratio"}
    assert matrices["default"]["matrix_id"] != matrices["ratio"]["matrix_id"]
    # ds の単一スロットを奪い合わない。
    assert context["runtime"]["dataset"].pp_matrix is None


def test_qc_processed_commits_qc_and_the_final_matrix_together(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    handlers["resolve_feature_bindings"](context)
    handlers["qc_raw"](context)
    handlers["preprocess"](context)
    outcome = handlers["qc_processed"](context)

    names = [r["output_name"] for r in outcome["result_refs"]]
    assert "qc" in names and "matrix" in names
    final = context["runtime"]["final_matrices"]["default"]
    assert final["stage"] == "finalized"
    assert final["parent_id"] == context["runtime"]["matrices"]["default"]["matrix_id"]


def test_qc_population_is_fixed_at_qc_raw_and_reused(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    handlers["resolve_feature_bindings"](context)
    handlers["qc_raw"](context)
    fixed = context["runtime"]["qc_population"]
    handlers["preprocess"](context)
    handlers["qc_processed"](context)
    assert context["runtime"]["qc_population"] == fixed


def test_statistics_use_the_final_matrix_only(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    handlers["resolve_feature_bindings"](context)
    handlers["qc_raw"](context)
    handlers["preprocess"](context)
    handlers["qc_processed"](context)

    context["statistic_id"] = "w1"
    context["stage_id"] = "statistics:w1"
    context["request"]["statistics"] = [{
        "statistic_id": "w1", "kind": "welch", "matrix_recipe_id": "default",
        "transform": "none", "feature_scope": {"mode": "all_eligible"},
        "reference_group": "control", "test_group": "treated",
        "q_threshold": 0.05, "log2fc_threshold": 1.0}]
    outcome = handlers["statistics"](context)

    assert outcome["status"] == "succeeded"
    result = context["runtime"]["statistics"]["w1"]
    final_id = context["runtime"]["final_matrices"]["default"]["matrix_id"]
    assert result["matrix_id"] == final_id


def test_a_statistic_that_cannot_run_is_recorded_not_skipped(tmp_path):
    handlers = metabolomics_handlers.build_handlers()
    context = _context(tmp_path)
    context["runtime"]["assay_evidence"] = _evidence(context["runtime"]["dataset"])
    handlers["resolve_feature_bindings"](context)
    handlers["qc_raw"](context)
    handlers["preprocess"](context)
    handlers["qc_processed"](context)

    context["statistic_id"] = "p1"
    context["stage_id"] = "statistics:p1"
    context["request"]["statistics"] = [{
        "statistic_id": "p1", "kind": "pca", "matrix_recipe_id": "default",
        "transform": "none", "feature_scope": {"mode": "all_eligible"},
        "scaling": "autoscale", "n_components": 5}]
    outcome = handlers["statistics"](context)
    # 2 feature しか無いので PCA は成立しうる/しないのどちらでもよいが、
    # どちらでも「結果JSON」は残る（評価不能も成果物）。
    assert [r["output_name"] for r in outcome["result_refs"]] == ["statistic:p1"]


# ---------- restore_results ----------

def test_restore_results_verifies_hashes(tmp_path):
    from metabolomix.pipeline import report as report_mod

    ref = report_mod.persist_result(tmp_path, {
        "output_name": "feature_bindings", "kind": "feature_bindings",
        "result_id": "res_binding_1",
        "data": {"schema": "feature-bindings.v1", "status": "resolved"}})
    record = {"identity": {"pipeline_root": str(tmp_path)}, "results": [ref]}

    restored = metabolomics_handlers.restore_results(record, tmp_path, _dataset())
    assert restored["bindings"]["status"] == "resolved"


def test_restore_results_refuses_a_tampered_result(tmp_path):
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import report as report_mod

    ref = report_mod.persist_result(tmp_path, {
        "output_name": "feature_bindings", "kind": "feature_bindings",
        "result_id": "res_binding_1",
        "data": {"schema": "feature-bindings.v1", "status": "resolved"}})
    record = {"identity": {"pipeline_root": str(tmp_path)}, "results": [ref]}

    target = tmp_path / ref["relative_path"]
    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["data"]["status"] = "needs_input"
    target.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(DomainError) as caught:
        metabolomics_handlers.restore_results(record, tmp_path, _dataset())
    assert caught.value.code == "RESULT_INTEGRITY_MISMATCH"


def test_restore_results_ignores_missing_results(tmp_path):
    record = {"identity": {"pipeline_root": str(tmp_path)}, "results": []}
    restored = metabolomics_handlers.restore_results(record, tmp_path, _dataset())
    assert restored == {}
