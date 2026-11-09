"""工程エンジン（`metabolomix.pipeline.engine`）とsessionを持たないworkerのテスト。

engineは注入handler（`tests.pipeline_worker_harness.synthetic_handlers`）だけで検証する。
run は `tests.pipeline_fixtures.create_synthetic_run`（入力計画を最小にした v3 の run）。
"""
from __future__ import annotations

import ast
import json
import subprocess
import sys
from pathlib import Path

import pytest

import metabolomix.pipeline.store as store_module
from metabolomix.core.atomic_io import DomainError, atomic_write_json
from metabolomix.pipeline.engine import build_stages, cancel_request_path, run_engine
from metabolomix.pipeline.store import load_run
from tests.pipeline_fixtures import create_synthetic_run
from tests.pipeline_worker_harness import synthetic_handlers

REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(tmp_path, polarities=("negative",)):
    path = create_synthetic_run(tmp_path, polarities=polarities)
    calls: list[str] = []
    return path, calls


# ---------- 工程の計画 ----------

def test_stage_plan_is_per_unit_then_receipt():
    request = {"units": [{"polarity": "negative"}, {"polarity": "positive"}]}
    plan = build_stages(request)
    assert [s["stage_id"] for s in plan] == [
        "prepare_input:negative", "upstream:negative", "validate_outputs:negative",
        "prepare_input:positive", "upstream:positive", "validate_outputs:positive",
        "receipt"]
    assert plan[1] == {"stage_id": "upstream:negative", "handler": "upstream",
                       "unit": "negative"}
    assert plan[-1] == {"stage_id": "receipt", "handler": "receipt", "unit": None}


def test_record_stages_follow_the_plan(tmp_path):
    path, _calls = _run(tmp_path, ("negative", "positive"))
    record = load_run(path)
    assert list(record["stages"]) == [
        "prepare_input:negative", "upstream:negative", "validate_outputs:negative",
        "prepare_input:positive", "upstream:positive", "validate_outputs:positive",
        "receipt"]
    assert record["schema"] == "pipeline-run.v3"
    assert set(record["units"]) == {"negative", "positive"}


# ---------- 単位の成否と最終の status ----------

def test_all_units_validated_is_completed_and_receipt_runs_last(tmp_path):
    path, calls = _run(tmp_path, ("negative", "positive"))
    result = run_engine(path, synthetic_handlers(calls))
    assert result["status"] == "completed"
    assert calls[-1] == "receipt"
    assert result["units"]["negative"]["outputs"]["alignment"]["arf"]
    assert result["units"]["positive"]["upstream"]["execution_id"] == "exec-fake-positive"
    assert (path / "results" / "receipt.json").is_file()


def test_failed_unit_does_not_stop_the_other_unit(tmp_path):
    path, calls = _run(tmp_path, ("negative", "positive"))
    result = run_engine(path, synthetic_handlers(calls, fail_stages={"upstream:negative"}))

    assert result["status"] == "partial"
    # 失敗した単位の残りの工程は呼ばれず pending のまま。
    assert "validate_outputs:negative" not in calls
    assert result["stages"]["upstream:negative"]["status"] == "failed"
    assert result["stages"]["validate_outputs:negative"]["status"] == "pending"
    # もう片方は最後まで通り、受領も走る。
    assert calls[-4:] == ["prepare_input:positive", "upstream:positive",
                          "validate_outputs:positive", "receipt"]
    receipt = json.loads((path / "results" / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["units"]["negative"]["status"] == "failed"
    assert receipt["units"]["negative"]["error"]["code"] == "SYNTHETIC_FAILURE"
    assert receipt["units"]["positive"]["status"] == "validated"
    assert 'polarity="positive"' in receipt["next"]


def test_all_units_failed_is_failed_without_receipt(tmp_path):
    path, calls = _run(tmp_path, ("negative", "positive"))
    handlers = synthetic_handlers(calls, fail_stages={"upstream:negative", "upstream:positive"})
    result = run_engine(path, handlers)
    assert result["status"] == "failed"
    assert "receipt" not in calls
    assert result["stages"]["receipt"]["status"] == "pending"
    assert not (path / "results" / "receipt.json").exists()


def test_unmatched_analysis_folder_is_partial(tmp_path):
    path, calls = _run(tmp_path)
    result = run_engine(path, synthetic_handlers(calls, unmatched_units={"negative"}))
    assert result["status"] == "partial"
    assert result["stages"]["receipt"]["status"] == "succeeded"


def test_a_failed_unit_is_not_retried_by_a_second_engine_pass(tmp_path):
    """失敗した工程を engine が自分で選び直さない（Console の自動再試行はしない）。"""
    path, calls = _run(tmp_path, ("negative", "positive"))
    run_engine(path, synthetic_handlers(calls, fail_stages={"upstream:negative"}))
    calls.clear()
    result = run_engine(path, synthetic_handlers(calls))
    assert "upstream:negative" not in calls
    assert result["stages"]["upstream:negative"]["status"] == "failed"
    assert result["status"] == "partial"


def test_a_succeeded_receipt_whose_file_vanished_is_rewritten(tmp_path):
    path, calls = _run(tmp_path)
    run_engine(path, synthetic_handlers(calls))
    (path / "results" / "receipt.json").unlink()
    calls.clear()
    result = run_engine(path, synthetic_handlers(calls))
    assert calls == ["receipt"]          # 上流は信頼して飛ばし、受領だけ書き直す
    assert result["status"] == "completed"


def test_engine_refuses_a_legacy_v1_run(tmp_path):
    path, _calls = _run(tmp_path)
    record = json.loads((path / "pipeline-run.json").read_text(encoding="utf-8"))
    record["schema"] = "pipeline-run.v1"
    atomic_write_json(path / "pipeline-run.json", record)
    with pytest.raises(DomainError) as excinfo:
        run_engine(path, synthetic_handlers())
    assert excinfo.value.code == "PIPELINE_RUN_LEGACY"


def test_an_exception_in_a_handler_fails_only_that_unit(tmp_path):
    path, calls = _run(tmp_path, ("negative", "positive"))
    handlers = synthetic_handlers(calls)

    def boom(context):
        if context["unit"] == "negative":
            raise DomainError("STAGED_INPUT_MISMATCH", "壊れた")
        return synthetic_handlers()["prepare_input"](context)

    handlers["prepare_input"] = boom
    result = run_engine(path, handlers)
    assert result["stages"]["prepare_input:negative"]["error"]["code"] == "STAGED_INPUT_MISMATCH"
    assert result["status"] == "partial"


def test_context_carries_the_unit(tmp_path):
    path, _calls = _run(tmp_path, ("negative", "positive"))
    seen = []
    handlers = synthetic_handlers()
    real = handlers["upstream"]

    def capture(context):
        seen.append((context["unit"], context["unit_index"], context["unit_root"].name,
                     context["unit_request"]["polarity"], context["unit_request"]["raw_dir"]))
        return real(context)

    handlers["upstream"] = capture
    run_engine(path, handlers)
    assert seen == [("negative", 0, "negative", "negative", "negative"),
                    ("positive", 1, "positive", "positive", "positive")]


# ---------- 取消 ----------

def test_cancel_requested_stops_before_next_stage(tmp_path):
    path, calls = _run(tmp_path, ("negative", "positive"))
    handlers = synthetic_handlers(calls)
    real = handlers["validate_outputs"]

    def then_cancel(context):
        outcome = real(context)
        atomic_write_json(cancel_request_path(path), {"cancel_requested": True})
        return outcome

    handlers["validate_outputs"] = then_cancel
    result = run_engine(path, handlers)
    assert result["status"] == "cancelled"
    assert "prepare_input:positive" not in calls


def test_a_stage_failing_under_cancellation_settles_as_cancelled(tmp_path):
    """工程の最中（Console）に取り消された run は `cancelled` に確定する。"""
    path, calls = _run(tmp_path, ("negative", "positive"))
    handlers = synthetic_handlers(calls)

    def cancelled_mid_stage(context):
        atomic_write_json(cancel_request_path(path), {"cancel_requested": True})
        return {"status": "failed", "result_refs": [], "warnings": [],
                "error": {"code": "MSDIAL_EXECUTION_FAILED",
                          "message": "termination=cancelled", "details": {}}}

    handlers["upstream"] = cancelled_mid_stage
    result = run_engine(path, handlers)
    assert result["status"] == "cancelled"
    assert result["stages"]["upstream:negative"]["error"]["code"] == "MSDIAL_EXECUTION_FAILED"
    assert "prepare_input:positive" not in calls


# ---------- state保存 ----------

def test_state_save_failure_propagates_and_keeps_prior_results(tmp_path, monkeypatch):
    path, calls = _run(tmp_path)
    real_save_run = store_module.save_run

    def flaky_save_run(pipeline_path, record, *, expected_revision):
        if record["stages"]["upstream:negative"]["status"] == "running":
            raise OSError("simulated disk failure")
        return real_save_run(pipeline_path, record, expected_revision=expected_revision)

    monkeypatch.setattr(store_module, "save_run", flaky_save_run)
    with pytest.raises(OSError, match="simulated disk failure"):
        run_engine(path, synthetic_handlers(calls))

    record = load_run(path)
    assert record["stages"]["prepare_input:negative"]["status"] == "succeeded"
    assert record["stages"]["upstream:negative"]["status"] == "pending"
    assert "upstream:negative" not in calls


def test_a_concurrent_state_revision_bump_does_not_kill_the_worker(tmp_path):
    """他アクターが`state_revision`を進めても、engineはworkerを落とさない。"""
    path, calls = _run(tmp_path)
    handlers = synthetic_handlers(calls)
    real = handlers["upstream"]

    def with_bump(context):
        outcome = real(context)
        latest = load_run(path)
        store_module.save_run(path, latest, expected_revision=latest["state_revision"])
        return outcome

    handlers["upstream"] = with_bump
    result = run_engine(path, handlers)
    assert result["status"] == "completed"
    assert result["units"]["negative"]["upstream"]["execution_id"] == "exec-fake-negative"


def test_second_engine_start_is_refused(tmp_path):
    path, calls = _run(tmp_path)
    handlers = synthetic_handlers(calls)
    real = handlers["upstream"]

    def reentrant(context):
        with pytest.raises(DomainError) as excinfo:
            run_engine(path, handlers)
        assert excinfo.value.code == "PIPELINE_ALREADY_RUNNING"
        return real(context)

    handlers["upstream"] = reentrant
    assert run_engine(path, handlers)["status"] == "completed"


# ---------- 稼働中workerが訂正・差し戻しを取り込む ----------

def test_a_manifest_update_during_the_run_rewrites_the_receipt(tmp_path):
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path)
    handlers = synthetic_handlers(calls)
    real = handlers["receipt"]
    done = {"count": 0}

    def receipt_then_update(context):
        outcome = real(context)
        if done["count"] == 0:
            done["count"] += 1
            prepare_resume(path, updates={"sample_manifest": "other.tsv"})
        return outcome

    handlers["receipt"] = receipt_then_update
    result = run_engine(path, handlers)
    assert calls.count("receipt") == 2
    assert result["request"]["revision"] == 2
    assert result["status"] == "completed"


def test_a_rerun_during_the_run_re_executes_the_failed_unit(tmp_path):
    """稼働中に `rerun_upstream=True` が失敗した単位を戻したら、このworkerがやり直す。"""
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path, ("negative", "positive"))
    handlers = synthetic_handlers(calls)
    failing = synthetic_handlers(fail_stages={"upstream:negative"})["upstream"]
    real_upstream = handlers["upstream"]
    state = {"failed_once": False, "reset": False}

    def upstream(context):
        if context["unit"] == "negative" and not state["failed_once"]:
            state["failed_once"] = True
            calls.append(context["stage_id"])
            return failing(context)
        return real_upstream(context)

    real_validate = handlers["validate_outputs"]

    def validate(context):
        outcome = real_validate(context)
        if context["unit"] == "positive" and not state["reset"]:
            state["reset"] = True
            prepare_resume(path, rerun_upstream=True)
        return outcome

    handlers["upstream"] = upstream
    handlers["validate_outputs"] = validate
    result = run_engine(path, handlers)
    assert calls.count("upstream:negative") == 2
    assert calls.count("upstream:positive") == 1        # 成功した単位はやり直さない
    assert result["status"] == "completed"


def test_an_update_right_before_a_stage_starts_restarts_the_pass(tmp_path, monkeypatch):
    """版の確認と `mark_stage_running` の間に届いた訂正で、古い要求のまま工程を回さない。"""
    from metabolomix.pipeline import engine
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path)
    handlers = synthetic_handlers(calls)
    real_receipt = handlers["receipt"]
    seen_manifests: list = []

    def receipt(context):
        seen_manifests.append(context["request"]["units"][0].get("sample_manifest"))
        return real_receipt(context)

    handlers["receipt"] = receipt
    real_cancel_requested = engine.cancel_requested
    count = {"n": 0}

    def cancel_requested_then_update(pipeline_path):
        count["n"] += 1
        if count["n"] == 4:  # receipt の直前（版の確認の後、running にする前）
            prepare_resume(path, updates={"sample_manifest": "other.tsv"})
        return real_cancel_requested(pipeline_path)

    monkeypatch.setattr(engine, "cancel_requested", cancel_requested_then_update)
    result = run_engine(path, handlers)

    assert seen_manifests == ["other.tsv"]
    assert result["request"]["revision"] == 2
    assert result["status"] == "completed"


def _failing_negative_then_resume_inside_positive_upstream(path, calls, resume):
    """negative の upstream を 1 回失敗させ、positive の upstream の**最中**に `resume()` を呼ぶ。"""
    handlers = synthetic_handlers(calls)
    failing = synthetic_handlers(fail_stages={"upstream:negative"})["upstream"]
    real_upstream = handlers["upstream"]
    state = {"failed_once": False, "resumed": False}

    def upstream(context):
        if context["unit"] == "negative" and not state["failed_once"]:
            state["failed_once"] = True
            calls.append(context["stage_id"])
            return failing(context)
        if context["unit"] == "positive" and not state["resumed"]:
            state["resumed"] = True
            resume()
        return real_upstream(context)

    handlers["upstream"] = upstream
    return handlers


def test_a_rerun_sent_while_a_console_runs_does_not_restart_that_console(tmp_path):
    """稼働中の worker が回している Console（upstream running）は rerun_upstream でも捨てない。"""
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path, ("negative", "positive"))
    accepted = {}

    def resume():
        accepted.update(prepare_resume(path, rerun_upstream=True))

    result = run_engine(path, _failing_negative_then_resume_inside_positive_upstream(
        path, calls, resume))
    assert set(accepted["reset_stage_ids"]) == {
        "prepare_input:negative", "upstream:negative", "validate_outputs:negative", "receipt"}
    assert calls.count("upstream:positive") == 1        # 走っていた Console をやり直さない
    assert calls.count("upstream:negative") == 2
    assert result["stages"]["upstream:positive"]["attempt"] == 1
    assert result["status"] == "completed"


def test_a_manifest_update_sent_while_a_console_runs_is_accepted(tmp_path):
    """稼働中の Console があっても、対応表の訂正は UPSTREAM_RERUN_REQUIRED にしない。"""
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path, ("negative", "positive"))
    accepted = {}

    def resume():
        accepted.update(prepare_resume(path, updates={"units": {
            "positive": {"sample_manifest": "positive/pos.tsv"}}}))

    result = run_engine(path, _failing_negative_then_resume_inside_positive_upstream(
        path, calls, resume))
    assert accepted["reset_stage_ids"] == ["receipt"]
    assert calls.count("upstream:positive") == 1
    assert result["request"]["revision"] == 2
    assert result["status"] == "partial"     # negative は明示が無いのでやり直さない


def test_endless_request_updates_stop_the_worker_instead_of_spinning(tmp_path):
    from metabolomix.pipeline.engine import _REQUEST_RELOAD_MAX_PASSES
    from metabolomix.pipeline.recovery import prepare_resume

    path, calls = _run(tmp_path)
    handlers = synthetic_handlers(calls)
    real = handlers["receipt"]
    bumps = {"count": 0}

    def receipt_then_always_correct(context):
        outcome = real(context)
        bumps["count"] += 1
        prepare_resume(path, updates={"sample_manifest": f"s{bumps['count']}.tsv"})
        return outcome

    handlers["receipt"] = receipt_then_always_correct
    with pytest.raises(DomainError) as excinfo:
        run_engine(path, handlers)
    assert excinfo.value.code == "REQUEST_REVISION_CHURN"
    assert bumps["count"] <= _REQUEST_RELOAD_MAX_PASSES


# ---------- AST: 禁止import ----------

_FORBIDDEN_IMPORT_PREFIXES = (
    "metabolomix.core.session_state",
    "metabolomix.core.mcp_core",
    "metabolomix.tools",
)


def _imported_module_names(source: str) -> set[str]:
    tree = ast.parse(source)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module)
    return names


@pytest.mark.parametrize("relative", [
    "metabolomix/pipeline/engine.py", "metabolomix/pipeline/worker.py",
    "metabolomix/pipeline/service.py", "metabolomix/pipeline/recovery.py",
    "metabolomix/pipeline/store.py", "metabolomix/pipeline/receipt.py",
])
def test_pipeline_modules_do_not_import_global_session(relative):
    source = (REPO_ROOT / relative).read_text(encoding="utf-8")
    for name in _imported_module_names(source):
        forbidden = any(name == prefix or name.startswith(prefix + ".")
                        for prefix in _FORBIDDEN_IMPORT_PREFIXES)
        assert not forbidden, f"{relative} が禁止importを含む: {name}"


# ---------- worker.py のCLI配線 ----------

def test_worker_main_invokes_run_worker_with_parsed_path(monkeypatch, tmp_path, capsys):
    from metabolomix.pipeline import worker

    seen = {}

    def fake_run_worker(path):
        seen["path"] = Path(path)
        return {"status": "completed"}

    monkeypatch.setattr(worker, "run_worker", fake_run_worker)
    exit_code = worker.main(["--pipeline", str(tmp_path)])
    assert exit_code == 0
    assert seen["path"] == tmp_path
    assert "status=completed" in capsys.readouterr().out


def test_worker_harness_runs_pipeline_as_real_subprocess(tmp_path):
    path, _calls = _run(tmp_path)
    harness = REPO_ROOT / "tests" / "pipeline_worker_harness.py"
    completed = subprocess.run(
        [sys.executable, str(harness), "--pipeline", str(path)],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    assert payload["status"] == "completed"


def test_handler_keys_cover_every_planned_stage():
    from metabolomix.pipeline.service import build_handlers

    handlers = build_handlers()
    for entry in build_stages({"units": [{"polarity": "negative"}]}):
        assert entry["handler"] in handlers
    assert set(handlers) == {"prepare_input", "upstream", "validate_outputs", "receipt"}
