"""pipeline の Console の前後の守り（2026-11-14 の棚卸しの直し）。

- B7-T63: 受付と prepare_input で切れた mzML を名指しで断る（`PIPELINE_INPUT_INVALID`）。
  Console が完了しなかったら見立て（`details.diagnosis`）を返す。
- B8-T66: 機械全体で Console は 1 本（待つ間は工程の `waiting`）、起こす前の commit の空きの
  警告（`LOW_COMMIT_MEMORY`）、終了コード 1450 を資源不足（`details.cause`）と読む。
- E2-T02: worker がまだ Windows の Job に入っていれば run の記録と戻り値に `WORKER_IN_JOB`。

計測ファイル・Console・worker の起動はすべて合成（偽の supervise と起動）。
"""
from __future__ import annotations

import copy
import threading
import time
from pathlib import Path

import pytest

from metabolomix.console import machine_slot
from metabolomix.core import process_control
from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import engine, recovery, service
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline.inputs import inspect_inputs, stage_inputs
from metabolomix.pipeline.store import load_run, save_run
from tests.config_fixtures import use_config
from tests.pipeline_fixtures import create_synthetic_run
from tests.pipeline_worker_harness import synthetic_handlers

_GB = 1024 ** 3
_PAD = "<spectrum/>\n" * 4000


# ---------- B7: 切れた mzML を Console の前で断る ----------

def _write_mzml(path: Path, *, closed: bool = True) -> Path:
    body = '<?xml version="1.0"?>\n<indexedmzML>\n<mzML>\n' + _PAD
    if closed:
        body += "</mzML>\n</indexedmzML>\n"
    path.write_text(body, encoding="utf-8")
    return path


def _mzml_source(root: Path, *, broken=()) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    for i in range(3):
        name = f"S{i}.mzML"
        _write_mzml(root / name, closed=name not in broken)
    (root / "fake.lbm2").write_bytes(b"fake-lbm-library")
    (root / "lab_param_202609050001.txt").write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nLbm file path: fake.lbm2\n",
        encoding="ascii", newline="\n")
    exe = root / "fake.exe"
    exe.write_text("fake console executable placeholder", encoding="ascii")
    return {"root": root, "exe": exe}


def _plan(source: dict, monkeypatch) -> dict:
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    root = source["root"]
    request = request_mod.unit_request(request_mod.resolve_request(root, None), 0, root)
    return inspect_inputs(root, request, exe_path=source["exe"])


def test_intake_names_a_truncated_mzml(tmp_path, monkeypatch):
    source = _mzml_source(tmp_path / "raw", broken={"S1.mzML"})
    with pytest.raises(DomainError) as caught:
        _plan(source, monkeypatch)
    assert caught.value.code == "PIPELINE_INPUT_INVALID"
    assert "S1.mzML" in caught.value.message
    assert caught.value.details["files"] == [
        {"name": "S1.mzML", "size": (source["root"] / "S1.mzML").stat().st_size,
         "reason": "truncated"}]
    assert caught.value.details["n_files"] == 1


def test_intake_accepts_closed_mzml(tmp_path, monkeypatch):
    plan = _plan(_mzml_source(tmp_path / "raw"), monkeypatch)
    assert plan["selected_format"] == "mzml"


def test_prepare_input_refuses_a_file_cut_after_intake(tmp_path, monkeypatch):
    source = _mzml_source(tmp_path / "raw")
    plan = _plan(source, monkeypatch)
    _write_mzml(source["root"] / "S2.mzML", closed=False)  # 受付の後に切れた

    pipeline_root = tmp_path / "pipeline_run"
    with pytest.raises(DomainError) as caught:
        stage_inputs(plan, pipeline_root)
    assert caught.value.code == "PIPELINE_INPUT_INVALID"
    assert [f["name"] for f in caught.value.details["files"]] == ["S2.mzML"]
    assert not (pipeline_root / "input").exists()  # 何も置かない


# ---------- upstream handler（偽の supervise） ----------

class _FakeJob:
    def __init__(self, status="completed", error=None):
        self.status = status
        self.warnings: list[str] = []
        self.error = error


def _upstream_context(tmp_path: Path, **extra) -> dict:
    pipeline_root = tmp_path / "pipeline"
    unit_root = pipeline_root / "negative"
    (unit_root / "input").mkdir(parents=True, exist_ok=True)
    return {
        "pipeline_root": pipeline_root, "pipeline_id": "pipe1",
        "stage_id": "upstream:negative", "attempt": 1, "unit": "negative",
        "unit_root": unit_root,
        "unit_request": {"measure": "peak_height", "save_project": False, "timeout_s": 60},
        "unit_state": {"inputs": {
            "method": {"effective_relative_path": None, "source_path": "C:/m.txt"},
            "polarity": {"value": "negative"}, "exe": {"path": "C:/fake.exe"}}},
        **extra,
    }


def _fake_console(monkeypatch, *, job=None, supervised=None, log_text=None, calls=None):
    job = job or _FakeJob()
    supervised = supervised or {"execution_id": "exec-1", "termination": "exited",
                                "exit_code": 0}

    def fake_supervise(job_path, *, cancel_path):
        if calls is not None:
            calls.append(job_path)
        if log_text is not None:
            (Path(job_path).parent / "msdial.log").write_text(log_text, encoding="utf-8")
        return dict(supervised)

    monkeypatch.setattr(service.console_execution, "supervise", fake_supervise)
    monkeypatch.setattr(service.job_manager, "load_job", lambda job_path: job)
    monkeypatch.setattr(service.job_manager, "count_raw_inputs", lambda root: 1)
    monkeypatch.setattr(service.job_manager, "list_raw_inputs", lambda root: [])
    monkeypatch.setattr(service.console_execution, "write_supervision_inputs",
                        lambda run_dir, payload: None)


def _commit(monkeypatch, available_gb: float) -> None:
    monkeypatch.setattr(process_control, "commit_memory_status", lambda: {
        "available_commit_bytes": int(available_gb * _GB),
        "commit_limit_bytes": 64 * _GB, "available_physical_bytes": int(available_gb * _GB)})


def test_low_commit_is_warned_and_console_still_runs(tmp_path, monkeypatch):
    calls: list = []
    _fake_console(monkeypatch, calls=calls)
    _commit(monkeypatch, 2.0)

    outcome = service._handle_upstream(_upstream_context(tmp_path))

    assert outcome["status"] == "succeeded"
    assert len(calls) == 1  # 止めない
    [warning] = [w for w in outcome["warnings"] if w["code"] == "LOW_COMMIT_MEMORY"]
    assert warning["details"]["available_commit_gb"] == 2.0
    assert outcome["unit_updates"]["upstream"]["memory_before"]["min_free_commit_gb"] == 8


def test_enough_commit_leaves_no_warning(tmp_path, monkeypatch):
    _fake_console(monkeypatch)
    _commit(monkeypatch, 48.0)
    outcome = service._handle_upstream(_upstream_context(tmp_path))
    assert [w["code"] for w in outcome["warnings"]] == []
    assert "console_slot_waited_s" not in outcome["unit_updates"]["upstream"]


def test_exit_code_1450_is_classified_as_resource_shortage(tmp_path, monkeypatch):
    _fake_console(monkeypatch, job=_FakeJob("failed", error="termination=exited exit_code=1450"),
                  supervised={"execution_id": "exec-1", "termination": "exited",
                              "exit_code": 1450})
    _commit(monkeypatch, 48.0)

    outcome = service._handle_upstream(_upstream_context(tmp_path))

    assert outcome["status"] == "failed"
    error = outcome["error"]
    assert error["code"] == "MSDIAL_EXECUTION_FAILED"
    assert error["details"]["cause"] == "system_resources"
    assert error["details"]["exit_code"] == 1450
    assert "1450" in error["message"]


def test_failed_console_names_the_file_from_its_log(tmp_path, monkeypatch):
    context = _upstream_context(tmp_path)
    raw = _write_mzml(context["unit_root"] / "input" / "S7.mzML")
    _fake_console(monkeypatch, job=_FakeJob("failed", error="termination=exited exit_code=1"),
                  supervised={"execution_id": "exec-1", "termination": "exited",
                              "exit_code": 1},
                  log_text=f"Unhandled Exception: System.Xml.XmlException in {raw}\n")
    monkeypatch.setattr(service.job_manager, "list_raw_inputs", lambda root: [raw])
    _commit(monkeypatch, 48.0)

    outcome = service._handle_upstream(context)

    details = outcome["error"]["details"]
    assert details["cause"] == "input_suspected"
    assert details["diagnosis"]["suspect_files"] == [
        {"name": "S7.mzML", "reason": "named_in_error"}]


def test_cancelled_console_carries_no_diagnosis(tmp_path, monkeypatch):
    _fake_console(monkeypatch, job=_FakeJob("failed", error="cancelled"),
                  supervised={"execution_id": "exec-1", "termination": "cancelled",
                              "exit_code": 1})
    _commit(monkeypatch, 48.0)
    details = service._handle_upstream(_upstream_context(tmp_path))["error"]["details"]
    assert "diagnosis" not in details and "cause" not in details


def test_upstream_waits_for_the_machine_wide_console_slot(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {})
    calls: list = []
    _fake_console(monkeypatch, calls=calls)
    _commit(monkeypatch, 48.0)
    monkeypatch.setattr(machine_slot, "POLL_INTERVAL_S", 0.02)
    reports: list = []
    context = _upstream_context(tmp_path, report_waiting=reports.append)

    inside = threading.Event()
    release = threading.Event()

    def other_run():
        with machine_slot.console_slot(owner={"pipeline_id": "other", "unit": "positive"}):
            inside.set()
            release.wait(10)

    holder = threading.Thread(target=other_run)
    holder.start()
    try:
        assert inside.wait(10)
        threading.Timer(0.3, release.set).start()
        outcome = service._handle_upstream(context)
    finally:
        release.set()
        holder.join(10)

    assert outcome["status"] == "succeeded"
    assert len(calls) == 1
    assert reports[0]["reason"] == "console_slot"
    assert reports[0]["holder"]["pipeline_id"] == "other"
    assert reports[-1] is None  # 順番が来たら待ちの表示を消す
    assert outcome["unit_updates"]["upstream"]["console_slot_waited_s"] >= 0


def test_cancel_while_waiting_does_not_start_console(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {})
    calls: list = []
    _fake_console(monkeypatch, calls=calls)
    monkeypatch.setattr(machine_slot, "POLL_INTERVAL_S", 0.02)
    context = _upstream_context(tmp_path)
    cancel = engine.cancel_request_path(context["pipeline_root"])

    def request_cancel(_waiting):
        cancel.parent.mkdir(parents=True, exist_ok=True)
        cancel.write_text("{}", encoding="utf-8")

    context["report_waiting"] = request_cancel
    with process_control.file_lock(machine_slot.console_lock_path()):
        outcome = service._handle_upstream(context)

    assert outcome["status"] == "failed"
    assert outcome["error"]["code"] == service.CONSOLE_WAIT_CANCELLED
    assert calls == []


# ---------- engine: 待ちの表示と worker の Job ----------

def test_waiting_is_shown_while_the_stage_runs_and_cleared_after(tmp_path):
    path = create_synthetic_run(tmp_path)
    seen: dict = {}
    handlers = synthetic_handlers()
    upstream = handlers["upstream"]

    def waiting_upstream(context):
        context["report_waiting"]({"reason": "console_slot", "message": "待っている"})
        seen["status"] = recovery.read_status(path)
        context["report_waiting"](None)
        seen["cleared"] = recovery.read_status(path).get("waiting")
        context["report_waiting"]({"reason": "console_slot", "message": "もう一度"})
        return upstream(context)

    handlers["upstream"] = waiting_upstream
    result = engine.run_engine(path, handlers)

    assert result["status"] == "completed"
    assert seen["status"]["waiting"] == {
        "upstream:negative": {"reason": "console_slot", "message": "待っている"}}
    assert seen["cleared"] is None
    # 工程の結果を書くときに待ちの表示を消す（残った表示を完了後に見せない）。
    assert "waiting" not in load_run(path)["stages"]["upstream:negative"]
    assert "waiting" not in recovery.read_status(path)


def test_waiting_is_not_written_on_a_stage_that_is_not_running(tmp_path):
    path = create_synthetic_run(tmp_path)
    engine.report_stage_waiting(path, "upstream:negative", {"reason": "console_slot"})
    assert "waiting" not in load_run(path)["stages"]["upstream:negative"]


def test_worker_in_a_job_is_warned_once_on_the_run(tmp_path):
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, True)
    engine._record_worker_identity(path, True)
    record = load_run(path)
    assert record["worker"]["in_job"] is True
    jobs = [w for w in record["warnings"] if w["code"] == "WORKER_IN_JOB"]
    assert len(jobs) == 1 and jobs[0]["stage_id"] is None
    assert "pipeline_resume" in jobs[0]["message"]


def test_worker_outside_any_job_is_not_warned(tmp_path):
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, False)
    record = load_run(path)
    assert record["worker"]["in_job"] is False
    assert "WORKER_IN_JOB" not in [w["code"] for w in record.get("warnings") or []]


def test_engine_without_the_check_records_nothing_about_jobs(tmp_path):
    """試験の harness（`run_engine` を直に呼ぶ）は確かめないので、記録は今までどおり。"""
    path = create_synthetic_run(tmp_path)
    assert engine.run_engine(path, synthetic_handlers())["status"] == "completed"
    record = load_run(path)
    assert "in_job" not in record["worker"]
    assert "WORKER_IN_JOB" not in [w["code"] for w in record.get("warnings") or []]


def test_production_worker_passes_its_own_job_check(tmp_path, monkeypatch):
    from metabolomix.pipeline import worker

    captured = {}
    monkeypatch.setattr(process_control, "current_process_in_job", lambda: True)
    monkeypatch.setattr(engine, "run_engine",
                        lambda path, handlers, *, worker_in_job=None:
                        captured.update(worker_in_job=worker_in_job) or {"status": "completed"})
    worker.run_worker(tmp_path)
    assert captured == {"worker_in_job": True}


# ---------- 起動の戻り値: WORKER_IN_JOB ----------

def _upstream_done(pipeline_root: Path) -> None:
    record = copy.deepcopy(load_run(pipeline_root))
    record["stages"]["upstream:negative"]["status"] = "succeeded"
    save_run(pipeline_root, record, expected_revision=record["state_revision"])


@pytest.mark.parametrize("in_job, expected", [(True, ["WORKER_IN_JOB"]), (False, []),
                                              (None, [])])
def test_resume_returns_a_warning_when_the_worker_stays_in_a_job(tmp_path, monkeypatch,
                                                                 in_job, expected):
    pipeline_root = create_synthetic_run(tmp_path)
    _upstream_done(pipeline_root)
    monkeypatch.setattr(service, "launch_detached", lambda command, *, cwd, log_path: {
        "pid": 999, "identity": {"pid": 999, "creation_time": 1}, "breakaway": "inherited_job",
        "relay_breakaway": "inherited_job", "worker_in_job": in_job, "log_path": str(log_path)})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.05)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.01)

    summary = service.resume_pipeline(pipeline_root)

    assert summary["launch"]["worker_in_job"] is in_job
    assert [w["code"] for w in summary.get("warnings") or []
            if w["code"] == "WORKER_IN_JOB"] == expected


def test_start_returns_a_warning_when_the_worker_stays_in_a_job(tmp_path, monkeypatch):
    from tests.pipeline_fixtures import make_source

    source = make_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    monkeypatch.setattr(service, "launch_detached", lambda command, *, cwd, log_path: {
        "pid": 999, "identity": {"pid": 999, "creation_time": 1}, "breakaway": "inherited_job",
        "relay_breakaway": "inherited_job", "worker_in_job": True, "log_path": str(log_path)})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.05)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.01)

    summary = service.start_pipeline(source["root"])

    assert summary["launched"] is True
    assert [w["code"] for w in summary["warnings"]] == ["WORKER_IN_JOB"]
