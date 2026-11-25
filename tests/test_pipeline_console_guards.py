"""pipeline の Console の前後の守り（2026-11-14 の棚卸しの直し）。

- B7-T63: 受付と prepare_input で切れた mzML を名指しで断る（`PIPELINE_INPUT_INVALID`）。
  Console が完了しなかったら見立て（`details.diagnosis`）を返す。
- B8-T66: 機械全体で Console は 1 本（待つ間は工程の `waiting`）、起こす前の commit の空きの
  警告（`LOW_COMMIT_MEMORY`）、終了コード 1450 を資源不足（`details.cause`）と読む。
- E2-T02: worker が、閉じたら止めるかの分からない Windows の Job に残っていれば run の記録と戻り値に
  `WORKER_IN_JOB`（直近の Job が閉じても止めない `inherited_job` なら警告しない）。
- 見立ては `pipeline_status` の既定の戻り値（`errors`）と受領の単位の `error` にも出る。

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


_HEADER_EMPTY = ('<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n<indexedmzML>\n<mzML>\n'
                 '<referenceableParamGroupList count="1">\n'
                 '<referenceableParamGroup id="CommonInstrumentParams"/>\n'
                 '</referenceableParamGroupList>\n'
                 '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">'
                 '<referenceableParamGroupRef ref="CommonInstrumentParams"/></instrumentConfiguration>'
                 '</instrumentConfigurationList>\n<run id="r">\n')


def _write_header_mzml(path: Path, header: str, *, encoding: str = "utf-8") -> Path:
    path.write_bytes((header + _PAD + "</run>\n</mzML>\n</indexedmzML>\n").encode(encoding))
    return path


def test_intake_names_an_mzml_with_an_empty_self_closing_param_group(tmp_path, monkeypatch):
    # 登録の前に書き直された mzML（ST004581）: 閉じタグは正しいが、Console はピーク 0 で終わる
    source = _mzml_source(tmp_path / "raw")
    _write_header_mzml(source["root"] / "S1.mzML", _HEADER_EMPTY)
    with pytest.raises(DomainError) as caught:
        _plan(source, monkeypatch)
    assert caught.value.code == "PIPELINE_INPUT_INVALID"
    assert "S1.mzML" in caught.value.message and "referenceableParamGroup" in caught.value.message
    assert caught.value.details["files"] == [
        {"name": "S1.mzML", "size": (source["root"] / "S1.mzML").stat().st_size,
         "reason": "empty_header_element", "element": "referenceableParamGroup"}]


def test_intake_accepts_the_open_close_param_group_and_the_self_closing_ref(tmp_path, monkeypatch):
    source = _mzml_source(tmp_path / "raw")
    fixed = _HEADER_EMPTY.replace('<referenceableParamGroup id="CommonInstrumentParams"/>',
                                  '<referenceableParamGroup id="CommonInstrumentParams"></referenceableParamGroup>')
    _write_header_mzml(source["root"] / "S1.mzML", fixed)
    assert _plan(source, monkeypatch)["selected_format"] == "mzml"


def test_header_check_reads_utf16_and_ignores_the_spectra(tmp_path):
    from metabolomix.console import raw_integrity
    utf16 = _write_header_mzml(tmp_path / "U.mzML", _HEADER_EMPTY, encoding="utf-16-le")
    assert raw_integrity.open_format_problem(utf16)["reason"] == "empty_header_element"
    # 見出しの後（スペクトルの中）の自己終了の要素は見ない
    late = tmp_path / "L.mzML"
    late.write_text('<?xml version="1.0"?>\n<indexedmzML>\n<mzML>\n<run id="r">\n'
                    '<referenceableParamGroup id="X"/>\n' + _PAD + "</run>\n</mzML>\n</indexedmzML>\n",
                    encoding="utf-8")
    assert raw_integrity.open_format_problem(late) is None


@pytest.mark.parametrize("element,dangerous", [
    ('<sampleList count="0"/>', "sampleList"),
    ('<fileContent/>', "fileContent"),
    ('<referenceableParamGroupList count="0"/>', "referenceableParamGroupList"),
    ('<softwareList count="0"/>', None),            # 読み込みが閉じタグを待たない要素は害が無い
    ('<sourceFileList count="0"/>', None),
])
def test_header_check_names_only_the_elements_msdial_reads_to_the_end_tag(tmp_path, element, dangerous):
    from metabolomix.console import raw_integrity
    path = tmp_path / "H.mzML"
    path.write_text('<?xml version="1.0"?>\n<indexedmzML>\n<mzML>\n' + element + '\n<run id="r">\n' + _PAD
                    + "</run>\n</mzML>\n</indexedmzML>\n", encoding="utf-8")
    problem = raw_integrity.open_format_problem(path)
    assert (problem or {}).get("element") == dangerous


@pytest.mark.parametrize("encoding", ["utf-8", "utf-16-le", "utf-16-be"])
def test_header_ignores_comments_cdata_and_quoted_delimiters(tmp_path, encoding):
    from metabolomix.console.raw_integrity import open_format_problem

    # コメントは読み込みの境界をまたぎ、中に終了の目印と自己終了タグの両方がある。
    comment = "<!--" + " " * 66000 + "<run><sampleList/> -->"
    prefix = ("<indexedmzML><mzML>" + comment + "<note><![CDATA[<sampleList/>]]></note>"
              + '<sampleList description="text /> end"></sampleList>')
    path = tmp_path / "header.mzML"
    _write_header_mzml(path, prefix + "<run>", encoding=encoding)
    assert open_format_problem(path) is None
    _write_header_mzml(path, prefix + '<fileContent attr="漢字 a > b"/><run>', encoding=encoding)
    assert open_format_problem(path)["element"] == "fileContent"


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


def _failing_upstream(handlers: dict, error: dict, *, unit: str = "negative") -> dict:
    upstream = handlers["upstream"]

    def failing(context):
        if context["unit"] == unit:
            return {"status": "failed", "result_refs": [], "warnings": [], "error": error}
        return upstream(context)

    handlers["upstream"] = failing
    return handlers


_CONSOLE_FAILURE = {
    "code": "MSDIAL_EXECUTION_FAILED",
    "message": "termination=exited exit_code=1 途中で切れた計測ファイルがある（suspect_files）。",
    "details": {"termination": "exited", "exit_code": 1, "unit": "negative",
                "cause": "broken_input",
                "diagnosis": {"cause": "broken_input", "message": "…",
                              "suspect_files": [{"name": "S1.mzML", "reason": "truncated"}],
                              "log_lines": ["Exception: S1.mzML"]}},
}


def test_status_returns_the_console_diagnosis_without_details(tmp_path):
    """1 極性の run が Console で落ちても（受領は無い）、既定の pipeline_status で原因が読める。"""
    path = create_synthetic_run(tmp_path)
    handlers = _failing_upstream(synthetic_handlers(), _CONSOLE_FAILURE)
    assert engine.run_engine(path, handlers)["status"] == "failed"

    status = recovery.read_status(path)

    assert "record" not in status and "receipt" not in status
    assert status["errors"] == {"upstream:negative": {
        "code": "MSDIAL_EXECUTION_FAILED", "message": _CONSOLE_FAILURE["message"],
        "cause": "broken_input", "suspect_files": [{"name": "S1.mzML", "reason": "truncated"}]}}


def test_status_names_the_files_prepare_input_refused(tmp_path):
    path = create_synthetic_run(tmp_path)
    handlers = synthetic_handlers()
    files = [{"name": f"S{i}.mzML", "size": 10, "reason": "truncated"} for i in range(12)]

    def refusing(context):
        raise DomainError("PIPELINE_INPUT_INVALID", "計測ファイルが途中で切れている",
                          {"files": files, "n_files": 25, "min_bytes": 32768})

    handlers["prepare_input"] = refusing
    engine.run_engine(path, handlers)

    error = recovery.read_status(path)["errors"]["prepare_input:negative"]
    assert error["code"] == "PIPELINE_INPUT_INVALID"
    assert error["suspect_files"] == [{"name": f"S{i}.mzML", "reason": "truncated"}
                                      for i in range(10)]
    assert error["n_suspect_files"] == 25


def test_status_has_no_errors_when_nothing_failed(tmp_path):
    path = create_synthetic_run(tmp_path)
    assert engine.run_engine(path, synthetic_handlers())["status"] == "completed"
    assert "errors" not in recovery.read_status(path)


def test_receipt_unit_error_carries_the_diagnosis(tmp_path):
    from metabolomix.pipeline import receipt as receipt_mod

    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    handlers = _failing_upstream(synthetic_handlers(), _CONSOLE_FAILURE)
    assert engine.run_engine(path, handlers)["status"] == "partial"

    error = receipt_mod.read_receipt(path, load_run(path))["units"]["negative"]["error"]
    assert error["cause"] == "broken_input"
    assert error["suspect_files"] == [{"name": "S1.mzML", "reason": "truncated"}]
    assert error["stage_id"] == "upstream:negative"
    assert "log_lines" not in error


def test_error_summary_clips_a_long_message():
    from metabolomix.pipeline import receipt as receipt_mod

    summary = receipt_mod.stage_error_summary({"code": "X", "message": "あ" * 2000, "details": {}})
    assert summary == {"code": "X", "message": "あ" * 600 + "…"}


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
    job = {"in_job": True, "nearest_job": None}
    engine._record_worker_identity(path, job)
    engine._record_worker_identity(path, job)
    record = load_run(path)
    assert record["worker"]["in_job"] is True
    jobs = [w for w in record["warnings"] if w["code"] == "WORKER_IN_JOB"]
    assert len(jobs) == 1 and jobs[0]["stage_id"] is None
    assert jobs[0]["details"] == {"nearest_job": None}
    assert "pipeline_resume" in jobs[0]["message"]


@pytest.mark.parametrize("nearest_job", [None, "unsupported", "silent_breakaway"])
def test_worker_in_a_job_that_may_kill_it_is_warned(tmp_path, nearest_job):
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, {"in_job": True, "nearest_job": nearest_job})
    [warning] = [w for w in load_run(path)["warnings"] if w["code"] == "WORKER_IN_JOB"]
    assert warning["details"] == {"nearest_job": nearest_job}


def test_worker_in_a_job_known_not_to_kill_on_close_is_not_warned(tmp_path):
    """直近の Job が閉じても止めない（`inherited_job`。Windows 11 のシェルの配下でいつもこうなる）なら
    `WORKER_IN_JOB` は出さず、その外側の段（kill-on-close の Job があるか）を確かめていないことだけを
    run の記録の `worker` に注記する（Job の問い合わせは直近の Job しか読めない）。"""
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, {"in_job": True, "nearest_job": "inherited_job"})
    record = load_run(path)
    assert record["worker"]["in_job"] is True
    assert record["worker"]["nearest_job"] == "inherited_job"
    assert "WORKER_IN_JOB" not in [w["code"] for w in record.get("warnings") or []]
    assert record["worker"]["outer_jobs"] == "unchecked"
    note = record["worker"]["outer_jobs_note"]
    assert "外側" in note and "kill-on-close" in note and "pipeline_resume" in note


@pytest.mark.parametrize("job", [{"in_job": True, "nearest_job": None},
                                 {"in_job": True, "nearest_job": "unsupported"},
                                 {"in_job": False, "nearest_job": None}])
def test_outer_jobs_note_only_when_the_nearest_job_is_inherited(tmp_path, job):
    """警告する場合（直近の Job が分からない・閉じたら止める）と Job の外では、外側の段の注記は付けない。"""
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, job)
    worker = load_run(path)["worker"]
    assert "outer_jobs" not in worker and "outer_jobs_note" not in worker


def test_worker_outside_any_job_is_not_warned(tmp_path):
    path = create_synthetic_run(tmp_path)
    engine._record_worker_identity(path, {"in_job": False, "nearest_job": None})
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
    state = {"in_job": True, "nearest_job": "inherited_job"}
    monkeypatch.setattr(process_control, "current_job_state", lambda: dict(state))
    monkeypatch.setattr(engine, "run_engine",
                        lambda path, handlers, *, worker_job=None:
                        captured.update(worker_job=worker_job) or {"status": "completed"})
    worker.run_worker(tmp_path)
    assert captured == {"worker_job": state}


# ---------- 起動の戻り値: WORKER_IN_JOB ----------

def _upstream_done(pipeline_root: Path) -> None:
    record = copy.deepcopy(load_run(pipeline_root))
    record["stages"]["upstream:negative"]["status"] = "succeeded"
    save_run(pipeline_root, record, expected_revision=record["state_revision"])


@pytest.mark.parametrize("in_job, relay_breakaway, expected", [
    # 中継は直近の Job を抜けさせたのに worker が Job にいる＝見えない外側の段（閉じたら止めるか分からない）。
    (True, "silent_breakaway", ["WORKER_IN_JOB"]),
    (True, "breakaway_ok", ["WORKER_IN_JOB"]),
    (True, None, ["WORKER_IN_JOB"]),
    # 中継の直近の Job が閉じても止めない Job（抜けられないので worker もそこにいる）。毎回の雑音にしない。
    (True, "inherited_job", []),
    (False, "silent_breakaway", []),
    (None, "silent_breakaway", []),
])
def test_resume_returns_a_warning_when_the_worker_stays_in_a_job(tmp_path, monkeypatch, in_job,
                                                                 relay_breakaway, expected):
    pipeline_root = create_synthetic_run(tmp_path)
    _upstream_done(pipeline_root)
    monkeypatch.setattr(service, "launch_detached", lambda command, *, cwd, log_path: {
        "pid": 999, "identity": {"pid": 999, "creation_time": 1}, "breakaway": "silent_breakaway",
        "relay_breakaway": relay_breakaway, "worker_in_job": in_job, "log_path": str(log_path)})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.05)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.01)

    summary = service.resume_pipeline(pipeline_root)

    assert summary["launch"]["worker_in_job"] is in_job
    warnings = [w for w in summary.get("warnings") or [] if w["code"] == "WORKER_IN_JOB"]
    assert [w["code"] for w in warnings] == expected
    for warning in warnings:
        assert warning["details"] == {"breakaway": "silent_breakaway",
                                      "relay_breakaway": relay_breakaway}
    # 警告しない inherited_job でも、外側の段を確かめていないことは起動の結果に注記する
    if in_job is True and relay_breakaway == "inherited_job":
        assert summary["launch"]["outer_jobs"] == "unchecked"
        assert "外側" in summary["launch"]["outer_jobs_note"]
    else:
        assert "outer_jobs" not in summary["launch"]


def test_start_returns_a_warning_when_the_worker_stays_in_a_job(tmp_path, monkeypatch):
    from tests.pipeline_fixtures import make_source

    source = make_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    monkeypatch.setattr(service, "launch_detached", lambda command, *, cwd, log_path: {
        "pid": 999, "identity": {"pid": 999, "creation_time": 1}, "breakaway": "silent_breakaway",
        "relay_breakaway": "silent_breakaway", "worker_in_job": True, "log_path": str(log_path)})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.05)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.01)

    summary = service.start_pipeline(source["root"])

    assert summary["launched"] is True
    assert [w["code"] for w in summary["warnings"]] == ["WORKER_IN_JOB"]
