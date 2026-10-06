"""v2 stage契約・pipeline-run.v2・analysis-job.v3（spec §6.1, §6.2）。

ここで縛るのは4つ。

1. **v2 stage列は`metabolomix.pipeline.stage_plan`が唯一の正準**。store（`record["stages"]`
   を作る側）とengine（実行順を決める側）が別々に同じ順序を持つと、片方だけ直したときに
   「recordに無いstageを計画する」または「計画に無いstageを実行しない」で静かにずれる。
   両者が同じbuilderを呼ぶことをテストで固定する。
2. **依存無効化表も同じモジュールに1つだけ置く**（spec §6.2の「metadata変更は
   resolve_metadata以降」等）。特にbinding訂正でConsoleを再実行しないこと。
3. **schema dispatchはreaderにもwriterにも要る**。v1のrecordはv1のまま読み、v1のまま
   書く（インプレース変換しない）。未知schemaは推測せず拒否する。
4. **job v3はprofile snapshotを持つときだけv3で書く**。持たないジョブの書き出しは
   v2のまま1バイトも変えない（v1経路の不変性）。

profileはrequest_v2が実際に参照するキーだけを持つ最小dictで足りる
（tests/test_metabolomics_request.py と同じ流儀）。
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import engine as engine_mod
from metabolomix.pipeline import recovery as recovery_mod
from metabolomix.pipeline import request_v2
from metabolomix.pipeline import store as store_mod
from metabolomix.pipeline.stage_plan import build_v2, invalidated_v2

REPO_ROOT = Path(__file__).resolve().parents[1]

#: spec §6.2のv2 stage列（角括弧展開前の固定部分）。実装の定数を参照せず、
#: specの本文からそのまま書き写す——実装側の並びが変わったらここが落ちる。
_BASE_V2_FROM_SPEC = [
    "prepare_inputs", "execute_console", "validate_outputs", "load_dataset",
    "resolve_metadata", "load_assay_evidence", "resolve_feature_bindings",
    "qc_raw", "preprocess", "qc_processed",
]


# ---------- 最小fixture ----------

def _profile():
    return {
        "matrix_recipes": {
            "default": {"base": "peak_height", "normalize": "none",
                        "drift_correct": False, "filter": None, "impute": "none"},
        },
        "feature_targets": {},
        "analysis_recipe": {"statistics": [], "internal_standards": []},
    }


def _pca(statistic_id="pca", **extra):
    base = {"statistic_id": statistic_id, "kind": "pca", "matrix_recipe_id": "default",
            "transform": "none", "feature_scope": {"mode": "all_eligible"}}
    base.update(extra)
    return base


def _welch(statistic_id="welch1", **extra):
    base = {"statistic_id": statistic_id, "kind": "welch", "matrix_recipe_id": "default",
            "transform": "log2", "feature_scope": {"mode": "all_eligible"},
            "reference_group": "control", "test_group": "treated"}
    base.update(extra)
    return base


def _anova(statistic_id="anova1", **extra):
    base = {"statistic_id": statistic_id, "kind": "anova_tukey",
            "matrix_recipe_id": "default", "transform": "none",
            "feature_scope": {"mode": "all_eligible"},
            "groups": ["a", "b", "c"]}
    base.update(extra)
    return base


def _resolved_v2(statistics=None, **extra):
    """実物の`request_v2.resolve`を通したv2 request（検証済みの形）。"""
    data = {"profile_file": "profile.json", "target": "auto"}
    if statistics is not None:
        data["statistics"] = statistics
    data.update(extra)
    return request_v2.resolve(data, _profile())


def _profile_with_target():
    """`feature_targets`を1件持つprofile（standard_assays/feature_bindings検証用）。"""
    profile = _profile()
    profile["feature_targets"] = {"IS1": {"role": "internal_standard"}}
    return profile


#: routineでのpreprocess上書きを1フィールドだけ許す証明書相当の許容集合
#: （`profile_schema.validate_certificate`の戻り値の同名キー）。
_ROUTINE_OVERRIDES = {"preprocess": {"default": {"impute": {"mode": "any"}}}}


def _resolved_v2_with_target(statistics=None, **extra):
    data = {"profile_file": "profile.json", "target": "auto"}
    if statistics is not None:
        data["statistics"] = statistics
    data.update(extra)
    return request_v2.resolve(data, _profile_with_target(),
                              routine_overrides=_ROUTINE_OVERRIDES)


def _plain_v2(*statistic_ids):
    """stage計画だけを見たいときの最小v2 request。"""
    return {"schema": request_v2.SCHEMA,
            "statistics": [{"statistic_id": sid} for sid in statistic_ids]}


def _source(tmp_path, name="source") -> Path:
    root = tmp_path / name
    root.mkdir(parents=True, exist_ok=True)
    return root


def _inputs(root: Path) -> dict:
    return {"source_root": str(root), "fingerprint": "a" * 64, "raw_inventory": []}


# ---------- brief記載のRED（原文のまま） ----------

def test_statistic_ids_expand_once():
    plan = build_v2({"statistics": [{"statistic_id": "a", "kind": "anova_tukey"}]})
    ids = [s["stage_id"] for s in plan]
    assert ids[-3:] == ["statistics:a", "export:a", "report"]
    assert len(ids) == len(set(ids))


# ---------- build_v2: 混在統計のstage ----------

def test_mixed_statistics_expand_in_request_order():
    request = _resolved_v2([_pca("explore"), _welch("w1"), _anova("a1")],
                           target="differential")
    ids = [entry["stage_id"] for entry in build_v2(request)]

    assert ids == _BASE_V2_FROM_SPEC + [
        "statistics:explore", "export:explore",
        "statistics:w1", "export:w1",
        "statistics:a1", "export:a1",
        "report",
    ]


def test_base_stages_carry_their_own_handler_and_no_statistic_id():
    plan = {entry["stage_id"]: entry for entry in build_v2(_plain_v2("s1"))}
    for stage_id in _BASE_V2_FROM_SPEC + ["report"]:
        assert plan[stage_id]["handler"] == stage_id
        assert plan[stage_id]["statistic_id"] is None


def test_per_statistic_stages_share_two_handlers():
    plan = {entry["stage_id"]: entry for entry in build_v2(_plain_v2("s1", "s2"))}
    assert plan["statistics:s1"]["handler"] == "statistics"
    assert plan["statistics:s1"]["statistic_id"] == "s1"
    assert plan["export:s2"]["handler"] == "export"
    assert plan["export:s2"]["statistic_id"] == "s2"


def test_report_is_always_the_last_stage():
    plan = build_v2(_plain_v2("z", "a"))
    assert plan[-1]["stage_id"] == "report"


def test_duplicate_statistic_ids_are_rejected_not_collapsed():
    """dictキーとして潰れると、stageが黙って1つ消える。"""
    with pytest.raises(DomainError) as caught:
        build_v2(_plain_v2("dup", "dup"))
    assert caught.value.code == "PIPELINE_REQUEST_INVALID"


def test_empty_statistics_is_rejected():
    with pytest.raises(DomainError) as caught:
        build_v2({"schema": request_v2.SCHEMA, "statistics": []})
    assert caught.value.code == "PIPELINE_REQUEST_INVALID"


def test_store_and_engine_share_one_builder():
    """順序の写しが2つあると片方だけ腐る。同じbuilderの出力であることを縛る。"""
    request = _resolved_v2([_pca("p"), _welch("w")], target="differential")
    assert store_mod.stage_ids_for(request) == [
        entry["stage_id"] for entry in engine_mod.build_stages(request)]


def test_v1_stage_ids_are_untouched_by_v2():
    """旧BASE_STAGE_IDSは変更しない（v1 goldenの不変性）。"""
    assert store_mod.STAGE_IDS == (
        "prepare_input", "upstream", "validate_outputs", "load_dataset",
        "resolve_metadata", "preprocess", "pca", "resolve_comparisons", "report")


# ---------- invalidated_v2: 依存無効化表 ----------

def test_binding_change_does_not_restart_console():
    dirty = invalidated_v2({"feature_bindings"}, _plain_v2("s"))
    assert "resolve_feature_bindings" in dirty
    assert "statistics:s" in dirty
    assert "execute_console" not in dirty
    assert "prepare_inputs" not in dirty
    assert "load_dataset" not in dirty


def test_standard_assays_change_starts_at_feature_bindings():
    dirty = invalidated_v2({"standard_assays"}, _plain_v2("s"))
    assert dirty == {
        "resolve_feature_bindings", "qc_raw", "preprocess", "qc_processed",
        "statistics:s", "export:s", "report"}


def test_metadata_change_starts_at_resolve_metadata():
    dirty = invalidated_v2({"metadata"}, _plain_v2("s"))
    assert dirty == {
        "resolve_metadata", "load_assay_evidence", "resolve_feature_bindings",
        "qc_raw", "preprocess", "qc_processed", "statistics:s", "export:s", "report"}
    assert invalidated_v2({"sample_manifest"}, _plain_v2("s")) == dirty


def test_recipe_change_starts_at_preprocess_and_leaves_raw_qc_alone():
    dirty = invalidated_v2({"preprocess"}, _plain_v2("s"))
    assert dirty == {"preprocess", "qc_processed", "statistics:s", "export:s", "report"}
    assert "qc_raw" not in dirty


def test_recipe_change_that_moves_the_raw_qc_scope_starts_at_qc_raw():
    dirty = invalidated_v2({"qc_raw_scope"}, _plain_v2("s"))
    assert dirty == {"qc_raw", "preprocess", "qc_processed",
                     "statistics:s", "export:s", "report"}


def test_one_statistic_change_leaves_the_others_alone():
    dirty = invalidated_v2({"statistics:w"}, _plain_v2("p", "w"))
    assert dirty == {"statistics:w", "export:w", "report"}


def test_a_wholesale_statistics_change_touches_every_statistic():
    dirty = invalidated_v2({"statistics"}, _plain_v2("p", "w"))
    assert dirty == {"statistics:p", "export:p", "statistics:w", "export:w", "report"}
    assert invalidated_v2({"target"}, _plain_v2("p", "w")) == dirty


def test_a_removed_statistic_still_invalidates_its_stages_and_the_report():
    """削除済み統計の成果物をcurrentのまま残さない（spec §6.2）。"""
    dirty = invalidated_v2({"statistics:gone"}, _plain_v2("kept"))
    assert dirty == {"statistics:gone", "export:gone", "report"}
    assert "statistics:kept" not in dirty


def test_unknown_change_token_is_rejected_not_ignored():
    """黙って何も無効化しないのが一番危ない。"""
    with pytest.raises(DomainError) as caught:
        invalidated_v2({"bogus_aspect"}, _plain_v2("s"))
    assert caught.value.code == "PIPELINE_STAGE_PLAN_INVALID"


def test_no_change_invalidates_nothing():
    assert invalidated_v2(set(), _plain_v2("s")) == set()


# ---------- store: schemaごとの読込・書込版 ----------

def test_v2_request_creates_a_pipeline_run_v2_with_v2_stages(tmp_path):
    root = _source(tmp_path)
    request = _resolved_v2([_pca("p")])
    pipeline_root = store_mod.create_run(root, request, _inputs(root))

    record = store_mod.load_run(pipeline_root)
    assert record["schema"] == store_mod.SCHEMA_V2 == "pipeline-run.v2"
    assert list(record["stages"]) == _BASE_V2_FROM_SPEC + [
        "statistics:p", "export:p", "report"]


def test_v1_request_still_creates_a_pipeline_run_v1(tmp_path):
    from metabolomix.pipeline.request import resolve_request

    root = _source(tmp_path)
    request = resolve_request(root)
    pipeline_root = store_mod.create_run(root, request, _inputs(root))

    record = store_mod.load_run(pipeline_root)
    assert record["schema"] == store_mod.SCHEMA == "pipeline-run.v1"
    assert set(record["stages"]) == {
        "prepare_input", "upstream", "validate_outputs", "load_dataset",
        "resolve_metadata", "preprocess", "pca", "resolve_comparisons", "report"}


def test_load_run_still_reads_an_existing_v1_record(tmp_path):
    """v1 recordはv1のまま読める（v2へ自動昇格しない）。"""
    from metabolomix.pipeline.request import resolve_request

    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, resolve_request(root), _inputs(root))
    raw = json.loads((pipeline_root / store_mod.RUN_FILENAME).read_text(encoding="utf-8"))
    assert raw["schema"] == "pipeline-run.v1"
    assert store_mod.load_run(pipeline_root)["schema"] == "pipeline-run.v1"


def test_load_run_rejects_an_unknown_schema(tmp_path):
    from metabolomix.pipeline.request import resolve_request

    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, resolve_request(root), _inputs(root))
    run_path = pipeline_root / store_mod.RUN_FILENAME
    raw = json.loads(run_path.read_text(encoding="utf-8"))
    raw["schema"] = "pipeline-run.v9"
    run_path.write_text(json.dumps(raw), encoding="utf-8")

    with pytest.raises(DomainError) as caught:
        store_mod.load_run(pipeline_root)
    assert caught.value.code == "PIPELINE_RUN_INVALID"


def test_create_run_rejects_an_unknown_request_schema(tmp_path):
    root = _source(tmp_path)
    request = {"schema": "pipeline-request.v9", "target": "auto"}
    with pytest.raises(DomainError) as caught:
        store_mod.create_run(root, request, _inputs(root))
    assert caught.value.code == "PIPELINE_REQUEST_INVALID"


def test_save_run_refuses_to_rewrite_a_v1_record_as_v2(tmp_path):
    """過去の記録をインプレース変換しない（spec §6, CONSTRAINTS）。"""
    from metabolomix.pipeline.request import resolve_request

    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, resolve_request(root), _inputs(root))
    record = store_mod.load_run(pipeline_root)
    record["schema"] = store_mod.SCHEMA_V2

    with pytest.raises(DomainError) as caught:
        store_mod.save_run(pipeline_root, record, expected_revision=record["state_revision"])
    assert caught.value.code == "PIPELINE_RUN_SCHEMA_IMMUTABLE"
    assert store_mod.load_run(pipeline_root)["schema"] == "pipeline-run.v1"


def test_save_run_keeps_v2_and_still_enforces_append_only_results(tmp_path):
    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, _resolved_v2([_pca("p")]), _inputs(root))

    record = store_mod.load_run(pipeline_root)
    record["results"] = [{"result_id": "r1"}]
    store_mod.save_run(pipeline_root, record, expected_revision=record["state_revision"])
    assert store_mod.load_run(pipeline_root)["schema"] == store_mod.SCHEMA_V2

    record = store_mod.load_run(pipeline_root)
    record["results"] = [{"result_id": "rewritten"}]
    with pytest.raises(DomainError) as caught:
        store_mod.save_run(pipeline_root, record, expected_revision=record["state_revision"])
    assert caught.value.code == "PIPELINE_RESULTS_IMMUTABLE"


def test_save_run_still_detects_a_state_revision_conflict_on_v2(tmp_path):
    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, _resolved_v2([_pca("p")]), _inputs(root))
    record = store_mod.load_run(pipeline_root)
    with pytest.raises(DomainError) as caught:
        store_mod.save_run(pipeline_root, record, expected_revision=99)
    assert caught.value.code == "STATE_REVISION_CONFLICT"


# ---------- engine: v2 planの文脈 ----------

def test_engine_context_carries_statistic_id_for_v2(tmp_path):
    root = _source(tmp_path)
    request = _resolved_v2([_welch("w")], target="differential")
    pipeline_root = store_mod.create_run(root, request, _inputs(root))
    record = store_mod.load_run(pipeline_root)

    context = engine_mod.make_context(
        record, record["stages"]["export:w"], {}, request)
    assert context["statistic_id"] == "w"
    assert context["comparison_id"] is None


def test_engine_never_re_runs_console_on_its_own_for_v2():
    """`execute_console`は永続状態を信頼する側（自動再試行しない）。"""
    stage = {"stage_id": "execute_console", "status": "succeeded", "result_refs": []}
    assert engine_mod.stage_inputs_unchanged(stage, pipeline_root=Path("."))


# ---------- recovery: 同じbuilder・同じ依存表 ----------

def test_recovery_reset_for_v2_uses_the_shared_invalidation_table():
    current = _resolved_v2([_pca("p"), _welch("w")], target="differential")
    merged = _resolved_v2([_pca("p"), _welch("w", q_threshold=0.01)],
                          target="differential")
    stage_ids = set(store_mod.stage_ids_for(merged))

    reset = recovery_mod._stages_to_reset(
        current, merged, rerun_upstream=False, stage_ids=stage_ids)

    assert reset == {"statistics:w", "export:w", "report"}


def test_recovery_reset_for_v2_handles_added_and_removed_statistics():
    current = _resolved_v2([_pca("p"), _welch("gone")], target="differential")
    merged = _resolved_v2([_pca("p"), _welch("added")], target="differential")
    stage_ids = set(store_mod.stage_ids_for(merged)) | {"statistics:gone", "export:gone"}

    reset = recovery_mod._stages_to_reset(
        current, merged, rerun_upstream=False, stage_ids=stage_ids)

    assert reset == {"statistics:added", "export:added",
                     "statistics:gone", "export:gone", "report"}


def test_recovery_reset_for_v2_treats_a_rewritten_sheet_as_metadata_change():
    request = _resolved_v2([_pca("p")])
    stage_ids = set(store_mod.stage_ids_for(request))

    reset = recovery_mod._stages_to_reset(
        request, request, rerun_upstream=False, stage_ids=stage_ids,
        manifest_content_changed=True)

    assert "resolve_metadata" in reset
    assert "execute_console" not in reset


def test_recovery_reset_for_v1_is_unchanged(tmp_path):
    """v1の依存区分はそのまま（golden）。"""
    from metabolomix.pipeline.request import resolve_request

    root = _source(tmp_path)
    comparisons = [{"comparison_id": "c1", "reference_group": "a", "test_group": "b"}]
    current = resolve_request(root, {"target": "differential", "comparisons": comparisons})
    merged = resolve_request(root, {
        "target": "differential",
        "comparisons": [{**comparisons[0], "q_threshold": 0.01}]})
    stage_ids = set(store_mod.stage_ids_for(merged))

    reset = recovery_mod._stages_to_reset(
        current, merged, rerun_upstream=False, stage_ids=stage_ids)

    assert reset == {"resolve_comparisons", "differential:c1", "export:c1", "report"}


def test_resume_of_a_v2_run_looks_at_execute_console_not_upstream(tmp_path, monkeypatch):
    """v2 recordに`upstream` stageは無い。上流の再実行判定はexecute_consoleで行う。"""
    monkeypatch.setenv("LIPIDMIX_PIPELINE_INDEX_DIR", str(tmp_path / "_index"))
    root = _source(tmp_path)
    pipeline_root = store_mod.create_run(root, _resolved_v2([_pca("p")]), _inputs(root))

    with pytest.raises(DomainError) as caught:
        recovery_mod.prepare_resume(pipeline_root)
    assert caught.value.code == "UPSTREAM_RERUN_REQUIRED"

    record = store_mod.load_run(pipeline_root)
    record["stages"]["execute_console"]["status"] = "succeeded"
    store_mod.save_run(pipeline_root, record, expected_revision=record["state_revision"])

    receipt = recovery_mod.prepare_resume(pipeline_root)
    assert receipt["reset_stage_ids"] == []


# ---------- handoff: analysis-job.v3 ----------

def _profile_snapshot(run_dir: Path) -> dict:
    return {
        "profile_id": "lcms-neg-v1",
        "profile_revision": 1,
        "profile_content_hash": "b" * 64,
        "adapter": {"adapter_id": "msdial-5", "dependency_keys": {}},
        "polarity": "negative",
        "measure": "peak_height",
        "method": {"source_path": str(run_dir / "params.txt"), "sha256": "c" * 64},
        "dependencies": [{"method_key": "Msp file path", "kind": "msp",
                          "source_path": str(run_dir / "lib.msp"), "sha256": "d" * 64,
                          "present": True}],
        "execution_environment": {"adapter_id": "msdial-5", "exe_sha256": "e" * 64},
        "effective_method_relative_path": "profile/effective-method.txt",
        "effective_method_sha256": "f" * 64,
        "plan_identity_hash": "0" * 64,
    }


def test_a_job_without_a_profile_snapshot_is_still_written_as_v2(tmp_path):
    from metabolomix.handoff.schema import SCHEMA_VERSION, AnalysisJob

    job = AnalysisJob(
        schema=SCHEMA_VERSION, job_id="j", status="planned", created_at="t",
        updated_at="t", dataset_root=str(tmp_path), input_count=0,
        software_name="MS-DIAL", software_version="", execution_mode="console",
        method_file="m.txt", omics="lipidomics", polarity="negative",
        measure="peak_height", run_dir=str(tmp_path))
    path = tmp_path / "analysis-job.json"
    job.save(path)

    written = json.loads(path.read_text(encoding="utf-8"))
    assert written["schema"] == "analysis-job.v2"
    assert "profile" not in written


def test_a_job_with_a_profile_snapshot_is_written_as_v3_and_round_trips(tmp_path):
    from metabolomix.handoff.schema import SCHEMA_VERSION_V3, AnalysisJob

    snapshot = _profile_snapshot(tmp_path)
    job = AnalysisJob(
        schema="analysis-job.v2", job_id="j", status="planned", created_at="t",
        updated_at="t", dataset_root=str(tmp_path), input_count=0,
        software_name="MS-DIAL", software_version="", execution_mode="console",
        method_file="m.txt", omics="metabolomics", polarity="negative",
        measure="peak_height", run_dir=str(tmp_path), profile_snapshot=snapshot)
    path = tmp_path / "analysis-job.json"
    job.save(path)

    written = json.loads(path.read_text(encoding="utf-8"))
    assert written["schema"] == SCHEMA_VERSION_V3 == "analysis-job.v3"
    assert written["profile"] == snapshot

    loaded = AnalysisJob.load(path)
    assert loaded.schema == "analysis-job.v3"
    assert loaded.profile_snapshot == snapshot
    assert loaded.dependencies == snapshot["dependencies"]
    assert loaded.environment == snapshot["execution_environment"]


def test_job_dependency_and_environment_manifests_come_from_the_snapshot(tmp_path):
    """写しを2つ持たない——manifestはsnapshotの中の1箇所だけが正準。"""
    from metabolomix.handoff.schema import AnalysisJob

    job = AnalysisJob(
        schema="analysis-job.v2", job_id="j", status="planned", created_at="t",
        updated_at="t", dataset_root=str(tmp_path), input_count=0,
        software_name="MS-DIAL", software_version="", execution_mode="console",
        method_file="m.txt", omics="metabolomics", polarity="negative",
        measure="peak_height", run_dir=str(tmp_path))
    assert job.dependencies == []
    assert job.environment == {}


def test_job_load_rejects_an_unknown_schema(tmp_path):
    from metabolomix.handoff.schema import AnalysisJob

    path = tmp_path / "analysis-job.json"
    path.write_text(json.dumps({"schema": "analysis-job.v9", "job_id": "j",
                                "status": "planned", "created_at": "t", "updated_at": "t",
                                "run_dir": str(tmp_path)}), encoding="utf-8")
    with pytest.raises(ValueError):
        AnalysisJob.load(path)


# ---------- mztab/loading: v3の完了検証 ----------

def _v3_job_with_outputs(tmp_path, *, effective_method_text="Ion mode: Positive\n"):
    """v3ジョブ（profile snapshot付き）と、その実効メソッド写しを作る。"""
    from metabolomix.console.execution import receipt_path, write_supervision_inputs
    from metabolomix.console.job_manager import create_job, save_job
    from metabolomix.handoff.schema import Artifact, MztabEntry, sha256_file
    from tests.pipeline_fixtures import execution_record, write_mztab

    source = _source(tmp_path)
    raws = [source / "S1.abf", source / "S2.abf"]
    for index, raw in enumerate(raws, start=1):
        raw.write_bytes(b"\x00" * (32 + index))
    method = source / "params.txt"
    method.write_text("Ion mode: Positive\n", encoding="ascii")

    job, job_path = create_job(source, method, "positive", "peak_height",
                               omics="metabolomics", input_count=len(raws))
    run_dir = Path(job.run_dir)
    mztab = run_dir / "msdial" / "Height_AlignmentResult_1.mzTab"
    mztab.parent.mkdir(parents=True, exist_ok=True)
    write_mztab(mztab, raws)

    effective = run_dir / "profile" / "effective-method.txt"
    effective.parent.mkdir(parents=True, exist_ok=True)
    effective.write_text(effective_method_text, encoding="ascii")

    snapshot = _profile_snapshot(run_dir)
    snapshot["effective_method_sha256"] = sha256_file(effective)

    job.primary_mztab_files = [MztabEntry(
        path="msdial/Height_AlignmentResult_1.mzTab", polarity="positive",
        measure="peak_height", sha256=sha256_file(mztab), root="run_dir")]
    job.artifacts = [Artifact(path="S1_1.pai2", role="sample_peaks", format="pai2",
                              sha256="x", root="dataset_root")]
    job.status = "completed"
    job.profile_snapshot = snapshot
    save_job(job, job_path)

    write_supervision_inputs(run_dir, {"raw_inventory": [str(p.resolve()) for p in raws]})
    receipt_path(run_dir).write_text(
        json.dumps(execution_record(job_id=job.job_id)), encoding="utf-8")
    return job_path, effective


def test_a_v3_job_with_an_intact_effective_method_is_verified(tmp_path):
    from metabolomix.mztab.loading import load_dataset_state

    job_path, _effective = _v3_job_with_outputs(tmp_path)
    ds = load_dataset_state(job_path=job_path)
    assert ds.source_verification == "verified"


def test_a_v3_job_whose_effective_method_was_altered_is_not_verified(tmp_path):
    """実際に走ったメソッドを裏取りできない出力をverifiedと名乗らせない（spec §6.1）。"""
    from metabolomix.mztab.loading import load_dataset_state

    job_path, effective = _v3_job_with_outputs(tmp_path)
    effective.write_text("Ion mode: Negative\n", encoding="ascii")

    ds = load_dataset_state(job_path=job_path)
    assert ds.source_verification == "legacy_unverified"
    assert any("実効メソッド" in w for w in ds.validation_result.get("warnings", []))


def test_a_v3_job_whose_effective_method_is_missing_is_not_verified(tmp_path):
    from metabolomix.mztab.loading import load_dataset_state

    job_path, effective = _v3_job_with_outputs(tmp_path)
    effective.unlink()

    ds = load_dataset_state(job_path=job_path)
    assert ds.source_verification == "legacy_unverified"


# ---------- AST: stage_planも禁止importを持たない ----------

_FORBIDDEN_IMPORT_PREFIXES = (
    "metabolomix.core.session_state",
    "metabolomix.core.mcp_core",
    "metabolomix.tools",
)


def test_stage_plan_does_not_import_the_global_session():
    """store/engine/recoveryがimportする以上、同じ規則の下に置く。"""
    source = (REPO_ROOT / "metabolomix/pipeline/stage_plan.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    for name in names:
        assert not any(name == prefix or name.startswith(prefix + ".")
                       for prefix in _FORBIDDEN_IMPORT_PREFIXES), \
            f"stage_plan.py が禁止importを含む: {name}"


# ---------- request_fingerprint: schemaごとの内容hash ----------

#: v1要求の内容hashの実測値（このタスクの変更前に`request_fingerprint`が返していた
#: 値そのもの）。v1の対象キー集合・正規化を1文字でも動かすとここが落ちる
#: ——受付冪等性（`find_or_create_run`）と再開時の`request_changed`判定は
#: どちらもこの値に乗っているため、v1 runの同一性が静かに変わってはいけない。
_V1_DEFAULT_FINGERPRINT = "ea87c236da35811be6148b234f7fe3987e833376d97ac886991706bbff1c7385"
_V1_DIFFERENTIAL_FINGERPRINT = (
    "0c15056a9c6f79577ad4c9e285eb2d86eee96a2fe75c7e658dbf924a3d2fa2bb")


def test_v1_request_fingerprint_is_unchanged(tmp_path):
    """v1 goldenの同一性: 既知のv1要求のcontent_hashは1バイトも変わらない。"""
    from metabolomix.pipeline.request import request_fingerprint, resolve_request

    root = _source(tmp_path)
    assert request_fingerprint(resolve_request(root)) == _V1_DEFAULT_FINGERPRINT
    differential = resolve_request(root, {
        "target": "differential",
        "comparisons": [{"comparison_id": "c1", "reference_group": "a",
                         "test_group": "b"}]})
    assert request_fingerprint(differential) == _V1_DIFFERENTIAL_FINGERPRINT


def test_v2_fingerprint_separates_requests_that_differ_only_in_statistics():
    """**統計定義だけ**が違う2つのv2要求を、同一内容と誤判定しない（spec §6.1）。

    targetも群も同じで、違うのは`statistics`の中身だけ。ここが衝突すると
    `find_or_create_run`が別の解析を1本のrunへ畳み、`prepare_resume`は
    統計の訂正を「何も変わっていない」と読む。
    """
    from metabolomix.pipeline.request import request_fingerprint

    two_components = _resolved_v2([_pca("p", n_components=2)])
    three_components = _resolved_v2([_pca("p", n_components=3)])
    assert two_components["target"] == three_components["target"]
    assert two_components["effective_target"] == three_components["effective_target"]
    assert request_fingerprint(two_components) != request_fingerprint(three_components)


def test_v2_fingerprint_separates_an_added_statistic_with_the_same_target():
    from metabolomix.pipeline.request import request_fingerprint

    one = _resolved_v2([_welch("w")], target="differential")
    two = _resolved_v2([_welch("w"), _anova("a")], target="differential")
    assert request_fingerprint(one) != request_fingerprint(two)


def test_v2_fingerprint_separates_every_resumable_v2_field():
    """v2固有フィールドはどれ一つとしてhashから落ちていない。"""
    from metabolomix.pipeline.request import request_fingerprint

    base = _resolved_v2_with_target([_pca("p")])
    base_hash = request_fingerprint(base)

    variants = {
        "statistics": _resolved_v2_with_target([_pca("p", n_components=3)]),
        "profile_file": _resolved_v2_with_target([_pca("p")], profile_file="other.json"),
        "execution_purpose": _resolved_v2_with_target(
            [_pca("p")], execution_purpose="validation"),
        "sample_manifest": _resolved_v2_with_target([_pca("p")],
                                                    sample_manifest="sheet.tsv"),
        "standard_assays": _resolved_v2_with_target(
            [_pca("p")], standard_assays={"IS1": ["std-01"]}),
        "preprocess": _resolved_v2_with_target(
            [_pca("p")], preprocess={"default": {"impute": "half_min"}}),
        "feature_bindings": _resolved_v2_with_target(
            [_pca("p")],
            feature_bindings={"dataset_hash": "a" * 64,
                              "selections": {"IS1": {"feature_id": "F1",
                                                     "reason": "手動選択"}}}),
    }
    for field, variant in variants.items():
        assert request_fingerprint(variant) != base_hash,             f"{field}の違いがcontent_hashに反映されていません"


def test_v2_fingerprint_ignores_provenance_only_differences():
    """value_sources/effective_targetは由来情報。同じ実効値なら同一内容。"""
    from metabolomix.pipeline.request import request_fingerprint

    request = _resolved_v2([_pca("p")])
    restated = dict(request)
    restated["value_sources"] = {"statistics": "explicit_update"}
    assert request_fingerprint(restated) == request_fingerprint(request)


def test_v2_fingerprint_is_stable_for_the_same_content():
    from metabolomix.pipeline.request import request_fingerprint

    assert (request_fingerprint(_resolved_v2([_pca("p"), _welch("w")], target="differential"))
            == request_fingerprint(_resolved_v2([_pca("p"), _welch("w")],
                                                target="differential")))


def test_a_v2_statistics_change_bumps_the_request_revision(tmp_path):
    """content_hashのschema dispatchが、store/recoveryの同一性判定まで効いている。"""
    root = _source(tmp_path)
    one = _resolved_v2([_pca("p", n_components=2)])
    two = _resolved_v2([_pca("p", n_components=3)])

    pipeline_root = store_mod.create_run(root, one, _inputs(root))
    recorded = store_mod.load_run(pipeline_root)["request"]["content_hash"]

    from metabolomix.pipeline.request import request_fingerprint
    assert recorded != request_fingerprint(two)


# ---------- docs/schema: 記録形式の正準文書 ----------

_SCHEMA_DOCS = REPO_ROOT / "docs" / "schema"


def _first_table_column(text: str, heading: str) -> list[str]:
    """`heading`節のmarkdown表から、第1列のバッククォート付き値を順に取り出す。"""
    import re

    section = text.split(heading, 1)[1].split("\n## ", 1)[0]
    return [match.group(1) for match in
            re.finditer(r"^\|\s*`([^`]+)`\s*\|", section, flags=re.MULTILINE)]


def test_pipeline_run_v2_doc_lists_the_canonical_stage_order():
    """stage列の写しが腐ったら落ちる（docs/workflowの腐敗防止と同じ流儀）。"""
    from metabolomix.pipeline import stage_plan

    doc = (_SCHEMA_DOCS / "pipeline-run-v2.md").read_text(encoding="utf-8")
    expected = (list(stage_plan.BASE_V2_STAGE_IDS)
                + [f"{handler}:<statistic_id>"
                   for handler in stage_plan.PER_STATISTIC_HANDLERS]
                + [stage_plan.FINAL_V2_STAGE_ID])
    assert _first_table_column(doc, "## stage列") == expected


def test_pipeline_run_v2_doc_lists_every_invalidation_token():
    from metabolomix.pipeline import stage_plan

    doc = (_SCHEMA_DOCS / "pipeline-run-v2.md").read_text(encoding="utf-8")
    documented = set(_first_table_column(doc, "## 依存無効化"))
    assert documented == set(stage_plan.CHANGE_ENTRY_POINTS) | {
        "statistics", "target", "statistics:<statistic_id>"}


def test_pipeline_run_v2_doc_names_both_schema_constants():
    doc = (_SCHEMA_DOCS / "pipeline-run-v2.md").read_text(encoding="utf-8")
    assert store_mod.SCHEMA_V2 in doc
    assert store_mod.SCHEMA in doc


def test_analysis_job_v3_doc_documents_the_v3_only_fields():
    from metabolomix.handoff.schema import SCHEMA_VERSION, SCHEMA_VERSION_V3

    doc = (_SCHEMA_DOCS / "analysis-job-v3.md").read_text(encoding="utf-8")
    assert SCHEMA_VERSION_V3 in doc
    assert SCHEMA_VERSION in doc
    documented = set(_first_table_column(doc, "## v3で追加されるフィールド"))
    assert documented == {"profile", "dependencies", "environment"}
