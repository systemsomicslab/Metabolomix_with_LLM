"""再開・取消・状態の読取表示のテスト（pipeline-run.v3）。

engine は注入 handler（`tests.pipeline_worker_harness.synthetic_handlers`）で回し、
run は `tests.pipeline_fixtures.create_synthetic_run`（入力計画を最小にした v3 の run）。
固定入力の再検証（D05）だけは、本物の `inspect_inputs` / `stage_inputs` で入力を配置する。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.process_control import process_identity
from metabolomix.pipeline.engine import cancel_request_path, run_engine
from metabolomix.pipeline.recovery import prepare_resume, read_status, request_cancel
from metabolomix.pipeline.store import create_run, load_run, save_run
import metabolomix.pipeline.recovery as recovery_module
import metabolomix.pipeline.store as store_module
import metabolomix.pipeline.receipt as receipt_mod

from tests.pipeline_fixtures import (create_synthetic_run, make_source,
                                     write_sample_manifest)
from tests.pipeline_worker_harness import synthetic_handlers

_UNIT_STAGES = ("prepare_input", "upstream", "validate_outputs")


def _unit_stage_ids(polarity: str) -> list[str]:
    return [f"{handler}:{polarity}" for handler in _UNIT_STAGES]


def _process_identity_of_self():
    return process_identity(os.getpid())


def _set(pipeline_path: Path, **fields) -> None:
    record = load_run(pipeline_path)
    record.update(fields)
    save_run(pipeline_path, record, expected_revision=record["state_revision"])


def _completed(tmp_path, polarities=("negative",)):
    path = create_synthetic_run(tmp_path, polarities=polarities)
    assert run_engine(path, synthetic_handlers())["status"] == "completed"
    return path


# ---------- read_status: 読むだけ・判定不能は dead と読まない ----------

@pytest.fixture
def run_with_lost_worker(tmp_path):
    """status=running・worker identity 設定済みの run（自分自身の生きた identity を使う）。"""
    path = create_synthetic_run(tmp_path)
    _set(path, status="running",
         worker={"identity": _process_identity_of_self(),
                 "started_at": "2026-09-05T00:00:00+00:00"})
    return path / store_module.RUN_FILENAME


def test_status_does_not_finalize_or_rewrite(run_with_lost_worker, monkeypatch):
    path = run_with_lost_worker
    before = path.read_bytes()
    monkeypatch.setattr("metabolomix.pipeline.recovery.same_process", lambda identity: False)
    status = read_status(path)
    assert status["observed_health"] == "worker_missing"
    assert path.read_bytes() == before
    assert status["status"] != "completed"


def test_unknown_worker_identity_is_unknown_not_dead(tmp_path):
    path = create_synthetic_run(tmp_path)
    _set(path, status="running", worker={"identity": None, "started_at": None})

    status = read_status(path)
    assert status["observed_health"] == "unknown"
    assert status["status"] == "running"
    assert status["recovery_hint"]["code"] == "SUPERVISION_UNKNOWN"


def test_real_dead_pid_is_worker_missing_without_mocking_same_process(tmp_path):
    path = create_synthetic_run(tmp_path)
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    identity = process_identity(proc.pid)
    proc.wait(timeout=30)
    _set(path, status="running",
         worker={"identity": identity, "started_at": "2026-09-05T00:00:00+00:00"})

    assert read_status(path)["observed_health"] == "worker_missing"


def test_status_of_a_completed_run_shows_units_and_receipt(tmp_path):
    path = _completed(tmp_path, ("negative", "positive"))

    status = read_status(path)
    assert status["status"] == "completed"
    assert status["observed_health"] == "ok"
    assert status["units"] == {"negative": "validated", "positive": "validated"}
    assert status["stage_statuses"]["receipt"] == "succeeded"
    receipt = status["receipt"]
    assert receipt["relative_path"] == "results/receipt.json"
    assert set(receipt["analysis_dirs"]) == {"negative", "positive"}
    assert Path(receipt["analysis_dirs"]["negative"]) == path / "negative" / "input"
    # 次の呼び出しは実際の pipeline の場所で示す（置き換え前の目印を残さない）。
    assert "<pipeline_path>" not in receipt["next"]
    assert "load_dataset" in receipt["next"]


def test_status_of_a_legacy_run_is_readable_and_marked(tmp_path):
    path = _legacy_v1_run(tmp_path)
    status = read_status(path)
    assert status["legacy"] is True
    assert "load_dataset" in status["note"]
    assert status["units"] == {"negative": "validated"}
    assert "receipt" not in status


# ---------- request_cancel ----------

def test_request_cancel_accepts_without_killing_and_is_idempotent(tmp_path):
    path = create_synthetic_run(tmp_path)

    first = request_cancel(path)
    assert first["accepted"] is True
    assert first["status"] != "cancelled"  # 受理≠確定
    assert request_cancel(path)["accepted"] is True

    flag = json.loads(cancel_request_path(path).read_text(encoding="utf-8"))
    assert flag["cancel_requested"] is True


def test_request_cancel_then_engine_stops_between_stages(tmp_path):
    path = create_synthetic_run(tmp_path)
    calls: list[str] = []
    handlers = synthetic_handlers(calls)
    original = handlers["upstream"]

    def upstream_then_cancel(context):
        outcome = original(context)
        request_cancel(path)
        return outcome

    handlers["upstream"] = upstream_then_cancel
    assert run_engine(path, handlers)["status"] == "cancelled"
    assert "validate_outputs:negative" not in calls
    assert "receipt" not in calls


# ---------- prepare_resume: 取消フラグ・running のまま残った run ----------

def _completed_then(tmp_path: Path, status: str, *, identity=None) -> Path:
    path = _completed(tmp_path)
    fields = {"status": status}
    if identity is not None:
        fields["worker"] = {"identity": identity, "started_at": "2026-09-05T00:00:00+00:00"}
    _set(path, **fields)
    return path


def test_resume_after_cancel_withdraws_the_cancel_flag(tmp_path):
    path = _completed_then(tmp_path, "cancelled")
    request_cancel(path)
    assert cancel_request_path(path).is_file()

    assert prepare_resume(path)["status"] == "planned"
    assert not cancel_request_path(path).exists()


def test_resume_of_a_lost_running_run_returns_to_planned(tmp_path, monkeypatch):
    path = _completed_then(tmp_path, "running", identity=_process_identity_of_self())
    monkeypatch.setattr(recovery_module, "same_process", lambda identity: False)
    assert read_status(path)["recovery_hint"]["code"] == "PIPELINE_INTERRUPTED"

    assert prepare_resume(path)["status"] == "planned"


def test_resume_never_disturbs_a_running_run_whose_worker_is_alive(tmp_path):
    path = _completed_then(tmp_path, "running", identity=_process_identity_of_self())
    request_cancel(path)
    assert read_status(path)["observed_health"] == "ok"

    receipt = prepare_resume(path)

    assert receipt["status"] == "running", "生きているworkerのrunを再開扱いにした"
    assert cancel_request_path(path).is_file(), "出された取消を握り潰した"
    assert load_run(path)["status"] == "running"


def test_resume_of_a_running_run_with_an_unresolvable_worker_does_nothing(tmp_path,
                                                                          monkeypatch):
    path = _completed_then(tmp_path, "running", identity=_process_identity_of_self())
    request_cancel(path)

    def _probe_fails(identity):
        raise OSError("probeそのものが失敗した")

    monkeypatch.setattr(recovery_module, "same_process", _probe_fails)
    assert read_status(path)["observed_health"] == "unknown"

    assert prepare_resume(path)["status"] == "running"
    assert cancel_request_path(path).is_file()


def test_resume_keeps_the_cancel_of_a_planned_run_whose_worker_has_not_registered_yet(tmp_path):
    """起動直後（identity をまだ刻んでいない）の窓で、出された取消を黙って消さない。"""
    path = _completed_then(tmp_path, "planned")
    _set(path, worker={"identity": None, "started_at": None})
    request_cancel(path)

    prepare_resume(path)
    assert cancel_request_path(path).is_file(), "起動直後の窓で、出された取消を握り潰した"

    _set(path, status="cancelled")
    prepare_resume(path)
    assert not cancel_request_path(path).exists()


# ---------- prepare_resume: 上流の再実行は明示だけ ----------

def test_upstream_not_succeeded_requires_explicit_rerun_upstream(tmp_path):
    path = create_synthetic_run(tmp_path)
    result = run_engine(path, synthetic_handlers(fail_stages={"upstream:negative"}))
    assert result["status"] == "failed"  # 解析フォルダが 1 つも無い

    with pytest.raises(DomainError, match="UPSTREAM_RERUN_REQUIRED"):
        prepare_resume(path)

    accepted = prepare_resume(path, rerun_upstream=True)
    assert accepted["reused"] is False
    assert set(accepted["reset_stage_ids"]) == {*_unit_stage_ids("negative"), "receipt"}
    record = load_run(path)
    assert record["stages"]["upstream:negative"]["status"] == "pending"
    assert record["status"] == "planned"


def test_rerun_upstream_resets_only_failed_units(tmp_path):
    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    first = run_engine(path, synthetic_handlers(fail_stages={"upstream:negative"}))
    assert first["status"] == "partial"
    positive_before = {sid: dict(first["stages"][sid]) for sid in _unit_stage_ids("positive")}

    # 明示しない再開は、失敗した極性の Console をやり直さない（上流には触れない）。
    without = prepare_resume(path)
    assert without["reset_stage_ids"] == []
    assert without["not_retried_units"] == ["negative"]
    assert "rerun_upstream" in without["hint"]
    record = load_run(path)
    assert record["stages"]["upstream:negative"]["status"] == "failed"
    assert record["stages"]["upstream:negative"]["attempt"] == 1

    # 明示すれば、失敗した極性の 3 工程と受領だけを戻す。
    accepted = prepare_resume(path, rerun_upstream=True)
    assert set(accepted["reset_stage_ids"]) == {*_unit_stage_ids("negative"), "receipt"}
    assert "not_retried_units" not in accepted
    record = load_run(path)
    for sid in _unit_stage_ids("positive"):
        assert record["stages"][sid]["status"] == "succeeded"
        assert record["stages"][sid]["attempt"] == positive_before[sid]["attempt"]

    calls: list[str] = []
    second = run_engine(path, synthetic_handlers(calls))
    assert second["status"] == "completed"
    assert calls == [*_unit_stage_ids("negative"), "receipt"]


def test_rerun_upstream_leaves_a_unit_running_under_a_live_worker(tmp_path):
    """worker が生きていて upstream が running の単位は、rerun_upstream でも戻さない。"""
    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    record = load_run(path)
    for sid, status in (("prepare_input:negative", "succeeded"), ("upstream:negative", "failed"),
                        ("prepare_input:positive", "succeeded"),
                        ("upstream:positive", "running")):
        record["stages"][sid]["status"] = status
        record["stages"][sid]["attempt"] = 1
    record["status"] = "running"
    record["worker"] = {"identity": _process_identity_of_self(),
                        "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(path, record, expected_revision=record["state_revision"])

    accepted = prepare_resume(path, rerun_upstream=True)

    assert set(accepted["reset_stage_ids"]) == {*_unit_stage_ids("negative"), "receipt"}
    after = load_run(path)
    assert after["stages"]["upstream:positive"]["status"] == "running"
    assert after["stages"]["upstream:positive"]["attempt"] == 1
    assert after["status"] == "running"

    # 対応表だけの訂正も、稼働中の Console を理由に止めない。
    manifest = prepare_resume(path, updates={"units": {
        "positive": {"sample_manifest": "positive/pos.tsv"}}})
    assert manifest["reset_stage_ids"] == ["receipt"]
    assert load_run(path)["stages"]["upstream:positive"]["status"] == "running"


def test_a_unit_running_under_a_lost_worker_still_requires_rerun(tmp_path, monkeypatch):
    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    record = load_run(path)
    record["stages"]["upstream:positive"]["status"] = "running"
    record["status"] = "running"
    record["worker"] = {"identity": _process_identity_of_self(),
                        "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(path, record, expected_revision=record["state_revision"])
    monkeypatch.setattr(recovery_module, "same_process", lambda identity: False)

    with pytest.raises(DomainError, match="UPSTREAM_RERUN_REQUIRED"):
        prepare_resume(path)
    accepted = prepare_resume(path, rerun_upstream=True)
    assert set(_unit_stage_ids("positive")) <= set(accepted["reset_stage_ids"])


def test_execution_unresolved_when_console_alive_but_unverified(tmp_path, monkeypatch):
    path = create_synthetic_run(tmp_path)
    run_dir = path / "negative" / "console" / "attempt-0001"
    run_dir.mkdir(parents=True)
    (run_dir / "execution-result.json").write_text(json.dumps({
        "schema": "console-execution.v1",
        "process_identity": {"pid": 4242, "creation_time": 7}}), encoding="utf-8")
    (run_dir / "analysis-job.json").write_text("{}", encoding="utf-8")

    record = load_run(path)
    record["units"]["negative"]["upstream"]["console_job_path"] = \
        "negative/console/attempt-0001/analysis-job.json"
    save_run(path, record, expected_revision=record["state_revision"])

    monkeypatch.setattr(recovery_module, "same_process", lambda identity: True)
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(path)
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(path, rerun_upstream=True)


def test_mismatched_console_job_path_convention_is_not_silently_verified(tmp_path):
    """`console_job_path` が job.json でなく run_dir を指していて終了証跡を読めないとき、
    未検証の上流を黙って素通りさせない（`rerun_upstream` でも回避できない）。"""
    path = create_synthetic_run(tmp_path)
    run_dir = path / "negative" / "console" / "attempt-0001"
    run_dir.mkdir(parents=True)
    (run_dir / "execution-result.json").write_text(json.dumps({
        "schema": "console-execution.v1",
        "process_identity": {"pid": 4242, "creation_time": 7}}), encoding="utf-8")

    record = load_run(path)
    assert record["stages"]["upstream:negative"]["status"] != "succeeded"
    record["units"]["negative"]["upstream"]["console_job_path"] = \
        "negative/console/attempt-0001"
    save_run(path, record, expected_revision=record["state_revision"])

    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(path)
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(path, rerun_upstream=True)


def test_verified_upstream_without_receipt_file_is_not_unresolved(tmp_path):
    """検証済み（upstream succeeded）の単位は、終了証跡が読めなくても監視対象に数えない。"""
    path = create_synthetic_run(tmp_path)
    record = load_run(path)
    record["stages"]["upstream:negative"]["status"] = "succeeded"
    record["units"]["negative"]["upstream"]["console_job_path"] = \
        "negative/console/attempt-0001/analysis-job.json"
    assert recovery_module._console_supervision_states(record, path) == []


# ---------- prepare_resume: 試料の対応表の訂正は受領だけを戻す ----------

def _manifest_text(group: str) -> str:
    return f"# schema = sample-manifest.v1\nsample_id\tgroup\nS0\t{group}\n"


def test_manifest_update_resets_only_receipt(tmp_path):
    path = _completed(tmp_path)
    source = Path(load_run(path)["identity"]["source_root"])
    (source / "corrected.tsv").write_text(_manifest_text("treated"), encoding="utf-8",
                                          newline="\n")

    accepted = prepare_resume(path, updates={"sample_manifest": "corrected.tsv"})

    assert accepted["reset_stage_ids"] == ["receipt"]
    assert accepted["request_revision"] == 2
    assert accepted["status"] == "planned"
    record = load_run(path)
    for sid in _unit_stage_ids("negative"):
        assert record["stages"][sid]["status"] == "succeeded"

    calls: list[str] = []
    assert run_engine(path, synthetic_handlers(calls))["status"] == "completed"
    assert calls == ["receipt"]
    receipt = json.loads((path / "results" / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["request_revision"] == 2
    assert receipt["units"]["negative"]["samples"]["manifest"]["path"] == "corrected.tsv"


def test_manifest_update_for_one_polarity_of_two(tmp_path):
    path = _completed(tmp_path, ("negative", "positive"))
    record = load_run(path)
    source = Path(record["identity"]["source_root"])
    # 受領は Console が読んだ raw（`<unit_root>/input`）でシートを解す。
    unit_input = path / record["units"]["positive"]["unit_root"] / "input"
    unit_input.mkdir(parents=True, exist_ok=True)
    (unit_input / "S0.mzML").write_bytes(b"raw")
    write_sample_manifest(source / "positive" / "pos.tsv", source_names=["S0.mzML"],
                          groups=["control"])

    # 単位の `sample_manifest` も `dataset_root` 相対（単位の raw_dir 相対ではない）。
    accepted = prepare_resume(
        path, updates={"units": {"positive": {"sample_manifest": "positive/pos.tsv"}}})
    assert accepted["reset_stage_ids"] == ["receipt"]

    assert run_engine(path, synthetic_handlers())["status"] == "completed"
    data = receipt_mod.read_receipt(path, load_run(path))
    assert data["units"]["positive"]["samples"]["manifest"]["path"] == "positive/pos.tsv"
    assert data["units"]["negative"]["samples"]["manifest"] is None
    assert "SAMPLE_MANIFEST_INVALID" not in [w["code"] for w in data["warnings"]]


def test_downstream_updates_are_refused_on_resume(tmp_path):
    path = _completed(tmp_path)
    with pytest.raises(DomainError) as caught:
        prepare_resume(path, updates={"comparisons": []})
    assert caught.value.code == "PIPELINE_DOWNSTREAM_REMOVED"
    with pytest.raises(DomainError) as caught:
        prepare_resume(path, updates={"method": "other.txt"})
    assert caught.value.code == "NEW_PIPELINE_REQUIRED"
    assert load_run(path)["status"] == "completed"


def _completed_run_with_manifest(tmp_path):
    path = create_synthetic_run(tmp_path)
    source = Path(load_run(path)["identity"]["source_root"])
    manifest = source / "sample-manifest.tsv"
    manifest.write_text(_manifest_text("control"), encoding="utf-8", newline="\n")
    assert run_engine(path, synthetic_handlers())["status"] == "completed"
    return path, manifest


def test_rewriting_the_same_manifest_resets_only_the_receipt(tmp_path):
    """パスは同じまま中身だけ直したシートを、無視せず受領へ反映させる。"""
    path, manifest = _completed_run_with_manifest(tmp_path)
    manifest.write_text(_manifest_text("treated"), encoding="utf-8", newline="\n")

    accepted = prepare_resume(path)

    assert accepted["status"] == "planned"
    assert accepted["reset_stage_ids"] == ["receipt"]
    calls: list[str] = []
    assert run_engine(path, synthetic_handlers(calls))["status"] == "completed"
    assert calls == ["receipt"]


def test_an_untouched_manifest_keeps_the_completed_run_as_is(tmp_path):
    path, _manifest = _completed_run_with_manifest(tmp_path)

    accepted = prepare_resume(path)

    assert accepted["status"] == "completed"
    assert accepted["reset_stage_ids"] == []


# ---------- 旧 v1 の run ----------

def _legacy_v1_run(tmp_path: Path) -> Path:
    pipeline_root = tmp_path / "legacy" / "pipeline_legacy"
    pipeline_root.mkdir(parents=True)
    record = {
        "schema": "pipeline-run.v1",
        "identity": {"pipeline_id": "legacy", "source_root": str(tmp_path / "legacy"),
                     "pipeline_root": str(pipeline_root)},
        "request": {"revision": 1, "request_id": None, "content_hash": "h",
                    "saved_path": "requests/revision-0001.json",
                    "effective_target": "exploratory"},
        "status": "completed", "state_revision": 3,
        "stages": {sid: {"stage_id": sid, "status": "succeeded"} for sid in (
            "prepare_input", "upstream", "validate_outputs", "load_dataset",
            "resolve_metadata", "preprocess", "pca", "resolve_comparisons", "report")},
        "upstream": {"console_job_path": None, "execution_id": "e",
                     "verification": {"status": "completed"}},
        "inputs": {"polarity": {"value": "negative", "source": "request_explicit"}},
        "results": [], "needs_input": None, "warnings": [], "worker": {},
    }
    (pipeline_root / "pipeline-run.json").write_text(json.dumps(record), encoding="utf-8")
    return pipeline_root


def test_legacy_run_resume_is_refused(tmp_path):
    path = _legacy_v1_run(tmp_path)
    before = (path / "pipeline-run.json").read_bytes()

    with pytest.raises(DomainError) as caught:
        prepare_resume(path, rerun_upstream=True)
    assert caught.value.code == "PIPELINE_RUN_LEGACY"
    assert "load_dataset" in caught.value.details["replacements"]["read"]
    with pytest.raises(DomainError) as caught:
        request_cancel(path)
    assert caught.value.code == "PIPELINE_RUN_LEGACY"
    assert not cancel_request_path(path).exists()
    assert (path / "pipeline-run.json").read_bytes() == before


# ---------- 冪等性 ----------

def test_same_request_id_resend_does_not_create_new_revision(tmp_path):
    path = _completed(tmp_path)
    first = prepare_resume(path, updates={"sample_manifest": "a.tsv"}, request_id="resume-1")
    assert first["reused"] is False
    second = prepare_resume(path, updates={"sample_manifest": "a.tsv"}, request_id="resume-1")
    assert second["reused"] is True
    assert second["request_revision"] == first["request_revision"]
    assert load_run(path)["request"]["revision"] == first["request_revision"]


def test_same_request_id_resend_with_different_updates_conflicts(tmp_path):
    path = _completed(tmp_path)
    prepare_resume(path, updates={"sample_manifest": "a.tsv"}, request_id="resume-1")
    with pytest.raises(DomainError, match="IDEMPOTENCY_CONFLICT"):
        prepare_resume(path, updates={"sample_manifest": "b.tsv"}, request_id="resume-1")


def test_same_updates_resent_without_request_id_is_a_no_op(tmp_path):
    path = _completed(tmp_path)
    first = prepare_resume(path, updates={"sample_manifest": "a.tsv"})
    second = prepare_resume(path, updates={"sample_manifest": "a.tsv"})
    assert second["request_revision"] == first["request_revision"]
    assert load_run(path)["request"]["revision"] == first["request_revision"]


def test_no_op_resume_does_not_disturb_a_completed_run(tmp_path):
    path = _completed(tmp_path)
    accepted = prepare_resume(path)
    assert accepted["reused"] is False
    assert accepted["reset_stage_ids"] == []
    assert load_run(path)["status"] == "completed"


def test_prepare_resume_absorbs_concurrent_state_bump(tmp_path, monkeypatch):
    path = _completed(tmp_path)
    real_load_run = store_module.load_run
    calls_n = {"n": 0}

    def flaky_load_run(p):
        record = real_load_run(p)
        calls_n["n"] += 1
        if calls_n["n"] == 1:
            concurrent = real_load_run(p)
            concurrent["status"] = "planned"
            store_module.save_run(p, concurrent, expected_revision=concurrent["state_revision"])
        return record

    monkeypatch.setattr(store_module, "load_run", flaky_load_run)

    accepted = prepare_resume(path, updates={"sample_manifest": "a.tsv"})
    assert accepted["reused"] is False
    assert calls_n["n"] >= 2


# ---------- D05: 固定入力・実効メソッドの改変を検出する ----------

def _run_with_real_inputs(tmp_path, monkeypatch):
    """本物の `inspect_inputs` / `stage_inputs` で入力を配置した、上流まで済んだ run。"""
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    from metabolomix.pipeline import request as request_mod
    from metabolomix.pipeline.inputs import inspect_inputs, stage_inputs

    source = make_source(tmp_path / "source")
    root = source["root"]
    request = request_mod.resolve_request(root, {"polarity": "negative", "save_project": False})
    plan = inspect_inputs(root, request_mod.unit_request(request, 0, root),
                          exe_path=source["exe"])
    path = create_run(root, request, {"fingerprint": "0" * 64, "units": {"negative": plan}})

    unit_root = path / "negative"
    record = load_run(path)
    staged = stage_inputs(
        store_module.absolutize_unit_inputs(record["units"]["negative"]["inputs"], unit_root),
        unit_root)
    record["units"]["negative"]["inputs"] = store_module.relativize_unit_inputs(staged, unit_root)
    record["stages"]["prepare_input:negative"]["status"] = "succeeded"
    record["stages"]["upstream:negative"]["status"] = "succeeded"
    save_run(path, record, expected_revision=record["state_revision"])
    return path, source


def test_resume_detects_tampered_effective_method(tmp_path, monkeypatch):
    path, _source = _run_with_real_inputs(tmp_path, monkeypatch)
    effective = path / "negative" / "inputs" / "effective-method.txt"
    effective.write_text(effective.read_text(encoding="ascii") + "\ntampered", encoding="ascii")

    with pytest.raises(DomainError, match="STAGED_INPUT_MISMATCH"):
        prepare_resume(path)


def test_resume_detects_changed_raw_input(tmp_path, monkeypatch):
    path, source = _run_with_real_inputs(tmp_path, monkeypatch)
    (source["root"] / "S0.wiff").write_text("raw-0-changed", encoding="ascii")

    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        prepare_resume(path)


def test_resume_with_rerun_upstream_skips_staged_input_check_for_rerun_units(tmp_path,
                                                                             monkeypatch):
    """やり直す単位は再検証しない（上流が自分で再検証・再配置する）。"""
    path, source = _run_with_real_inputs(tmp_path, monkeypatch)
    _set(path, stages={**load_run(path)["stages"],
                       "validate_outputs:negative": {"stage_id": "validate_outputs:negative",
                                                     "status": "failed"}})
    (source["root"] / "S0.wiff").write_text("raw-0-changed", encoding="ascii")

    accepted = prepare_resume(path, rerun_upstream=True)
    assert accepted["reused"] is False
    assert "prepare_input:negative" in accepted["reset_stage_ids"]


# ---------- 受領の成果物の改変 ----------

def test_tampered_receipt_is_rewritten_without_touching_the_units(tmp_path):
    path = _completed(tmp_path)
    (path / "results" / "receipt.json").write_text("{}", encoding="utf-8")

    calls: list[str] = []
    assert run_engine(path, synthetic_handlers(calls))["status"] == "completed"
    assert calls == ["receipt"]


def test_untampered_run_is_skipped_on_second_engine_pass(tmp_path):
    path = _completed(tmp_path)
    calls: list[str] = []
    assert run_engine(path, synthetic_handlers(calls))["status"] == "completed"
    assert calls == []
