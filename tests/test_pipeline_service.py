"""metabolomix.pipeline.service の単体テスト（受付・起動・実 handler）。

`tests/test_pipeline_engine.py` / `test_pipeline_recovery.py` は合成 handler で engine 自体を
見る。ここでは受付（`plan_pipeline` / `start_pipeline` / `resume_pipeline`）の戻り値と起動の
規律、`build_handlers()` が組み立てる実 handler を見る。
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

from metabolomix.pipeline import service
from metabolomix.pipeline.store import load_run, save_run
from tests.config_fixtures import use_config
from tests.pipeline_fixtures import create_synthetic_run, make_two_polarity_source

_FAKE_CONSOLE = Path(__file__).resolve().parent / "fixtures" / "fake_console.py"


# ---------- 実 handler の単体 ----------

def test_handlers_are_the_four_upstream_handlers():
    assert set(service.build_handlers()) == {
        "prepare_input", "upstream", "validate_outputs", "receipt"}


def test_upstream_handler_calls_supervise_directly(tmp_path, monkeypatch):
    """`upstream` handler は共通監視経路 `execution.supervise` を直接呼び、
    attempt ごとの run_dir を `<極性>/console/attempt-NNNN` に作る。"""
    supervise_calls = []

    def fake_supervise(job_path, *, cancel_path):
        supervise_calls.append((job_path, cancel_path))
        return {"execution_id": "exec-1", "termination": "exited", "exit_code": 0}

    class _FakeJob:
        status = "completed"
        warnings: list[str] = []
        error = None

    monkeypatch.setattr(service.console_execution, "supervise", fake_supervise)
    monkeypatch.setattr(service.job_manager, "load_job", lambda job_path: _FakeJob())
    monkeypatch.setattr(service.job_manager, "count_raw_inputs", lambda root: 1)
    monkeypatch.setattr(service.job_manager, "list_raw_inputs", lambda root: [])
    monkeypatch.setattr(service.console_execution, "write_supervision_inputs",
                        lambda run_dir, payload: None)

    pipeline_root = tmp_path / "pipeline"
    unit_root = pipeline_root / "negative"
    (unit_root / "input").mkdir(parents=True)
    context = {
        "pipeline_root": pipeline_root,
        "pipeline_id": "pipe1",
        "stage_id": "upstream:negative",
        "attempt": 1,
        "unit": "negative",
        "unit_root": unit_root,
        "unit_request": {"measure": "peak_height", "save_project": False, "timeout_s": 60},
        "unit_state": {"inputs": {
            "method": {"effective_relative_path": None, "source_path": "C:/m.txt"},
            "polarity": {"value": "negative"}, "exe": {"path": "C:/fake.exe"}}},
    }

    outcome = service._handle_upstream(context)

    assert len(supervise_calls) == 1
    assert supervise_calls[0][0] == unit_root / "console" / "attempt-0001" / "analysis-job.json"
    assert outcome["status"] == "succeeded"
    upstream = outcome["unit_updates"]["upstream"]
    assert upstream["console_job_path"] == "negative/console/attempt-0001/analysis-job.json"
    assert upstream["execution_id"] == "exec-1"


def test_launch_pipeline_worker_uses_sys_executable_and_repo_root_cwd(tmp_path, monkeypatch):
    captured = {}

    def fake_launch_detached(command, *, cwd, log_path):
        captured.update(command=command, cwd=cwd, log_path=log_path)
        return {"launched": True, "pid": 4321}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    pipeline_path = tmp_path / "pipeline"
    (pipeline_path / "control").mkdir(parents=True)

    assert service.launch_pipeline_worker(pipeline_path) == {"launched": True, "pid": 4321}
    assert captured["command"][0] == sys.executable
    assert captured["command"][1:4] == ["-m", "metabolomix.pipeline.worker", "--pipeline"]
    assert captured["command"][4] == str(pipeline_path)
    assert captured["cwd"] == service._REPO_ROOT


# ---------- resume_pipeline: pipeline-run.json 自身のパスも受け付ける ----------

def _upstream_done(pipeline_root: Path, *, worker: dict | None = None) -> None:
    record = copy.deepcopy(load_run(pipeline_root))
    record["stages"]["upstream:negative"]["status"] = "succeeded"
    if worker is not None:
        record["worker"] = worker
    save_run(pipeline_root, record, expected_revision=record["state_revision"])


def test_resume_pipeline_normalizes_the_run_json_path_before_launching(tmp_path, monkeypatch):
    from metabolomix.pipeline import store as store_mod

    pipeline_root = create_synthetic_run(tmp_path)
    _upstream_done(pipeline_root)
    captured = {}

    def fake_launch_detached(command, *, cwd, log_path):
        captured.update(command=command, log_path=log_path)
        return {"launched": True, "pid": 555}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    receipt = service.resume_pipeline(pipeline_root / store_mod.RUN_FILENAME)

    assert receipt["pipeline_path"] == str(pipeline_root)
    assert captured["command"][-1] == str(pipeline_root)
    assert captured["log_path"] == pipeline_root / store_mod.CONTROL_SUBDIR / "worker-launch.log"


# ---------- 受付: 再利用した活動中・終端の run の上に二重起動しない ----------

def _prepare_source(tmp_path, monkeypatch):
    """`start_pipeline` の実処理（受付冪等性込み）を通すための生データフォルダ。"""
    from tests.pipeline_fixtures import make_source
    source = make_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")
    return source


def _counting_launch(monkeypatch, *, pid_base: int = 1000) -> list:
    calls: list = []

    def fake_launch_detached(command, *, cwd, log_path):
        calls.append(command)
        return {"launched": True, "pid": pid_base + len(calls)}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)
    return calls


def test_start_pipeline_does_not_relaunch_over_an_active_reused_run(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls = _counting_launch(monkeypatch)

    first = service.start_pipeline(source["root"])
    assert first["launched"] is True
    assert first["units"] == ["negative"]

    pipeline_path = Path(first["pipeline_path"])
    record = copy.deepcopy(load_run(pipeline_path))
    record["status"] = "running"
    record["worker"] = {"identity": {"pid": 424242, "creation_time": 1},
                        "started_at": "2026-09-06T00:00:00+00:00"}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    second = service.start_pipeline(source["root"])
    assert second["pipeline_path"] == first["pipeline_path"]
    assert second["launched"] is False
    assert len(launch_calls) == 1


def test_start_pipeline_does_not_relaunch_over_a_reused_completed_run(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls = _counting_launch(monkeypatch, pid_base=2000)

    first = service.start_pipeline(source["root"])
    pipeline_path = Path(first["pipeline_path"])
    record = copy.deepcopy(load_run(pipeline_path))
    record["status"] = "completed"
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    assert service.start_pipeline(source["root"])["launched"] is False
    assert len(launch_calls) == 1


# ---------- handshake は「新しい identity」だけを confirmed とする ----------

def test_resume_handshake_is_not_confirmed_when_worker_identity_is_unchanged(tmp_path, monkeypatch):
    pipeline_root = create_synthetic_run(tmp_path)
    _upstream_done(pipeline_root, worker={"identity": {"pid": 111111, "creation_time": 1},
                                          "started_at": "2026-09-05T00:00:00+00:00"})
    monkeypatch.setattr(service, "launch_detached",
                        lambda command, *, cwd, log_path: {"launched": True, "pid": 999})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.2)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    assert service.resume_pipeline(pipeline_root)["launch"]["handshake"] == "not_confirmed"


def test_resume_handshake_confirms_on_a_genuinely_new_identity(tmp_path, monkeypatch):
    pipeline_root = create_synthetic_run(tmp_path)
    _upstream_done(pipeline_root, worker={"identity": {"pid": 111111, "creation_time": 1},
                                          "started_at": "2026-09-05T00:00:00+00:00"})

    def fake_launch_detached(command, *, cwd, log_path):
        current = copy.deepcopy(load_run(pipeline_root))
        current["worker"] = {"identity": {"pid": 222222, "creation_time": 2},
                             "started_at": "2026-09-06T00:00:00+00:00"}
        save_run(pipeline_root, current, expected_revision=current["state_revision"])
        return {"launched": True, "pid": 222222}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.5)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    assert service.resume_pipeline(pipeline_root)["launch"]["handshake"] == "confirmed"


# ---------- plan_pipeline / start_pipeline の戻り値 ----------

def test_plan_pipeline_never_launches_a_worker(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)

    def poison(*a, **kw):
        raise AssertionError("plan_pipelineがworkerを起動した")

    monkeypatch.setattr(service, "launch_detached", poison)

    receipt = service.plan_pipeline(source["root"])
    assert receipt["launched"] is False
    assert receipt["status"] == "planned"
    assert load_run(Path(receipt["pipeline_path"]))["status"] == "planned"


def test_plan_receipt_shows_the_resolved_settings_per_unit(tmp_path, monkeypatch):
    """`pipeline_plan` は単位ごとの解決結果（メソッド・LBM・極性…）を起動前に見せる。
    大きな中間データ（raw 一覧・stat・entries・overrides）は載せない。"""
    source = _prepare_source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source["root"])

    assert receipt["units"] == ["negative"]
    resolved = receipt["resolved"]["units"]["negative"]
    assert resolved["method"]["source_path"] == str(Path(source["method"]).resolve())
    assert len(resolved["method"]["sha256"]) == 64
    assert resolved["lbm"]["path"] == str(Path(source["lbm"]).resolve())
    assert resolved["polarity"]["value"] == "negative"
    assert resolved["polarity"]["source"]
    assert set(resolved) == {"method", "lbm", "polarity", "library_mode", "msp",
                             "alignment_reference"}
    assert resolved["library_mode"] == "auto"
    assert resolved["msp"] == {"file": None, "sha256": None, "source": "not_used"}
    text = json.dumps(receipt, ensure_ascii=False)
    assert "raw_stat" not in text and "entries" not in text
    assert len(text) < 1500, len(text)

    # 起動する側（pipeline_run）の要約は最小のまま。
    assert "resolved" not in service._dispatch_summary(
        Path(receipt["pipeline_path"]), launched=False)


def test_plan_of_two_polarities_resolves_each_unit(tmp_path, monkeypatch):
    source = make_two_polarity_source(tmp_path / "source")
    use_config(monkeypatch, tmp_path, {"msdial.exe": source["negative"]["exe"]})
    monkeypatch.setattr("metabolomix.console.runner.probe_console_exe", lambda *a, **k: "console")

    receipt = service.plan_pipeline(source["root"], {"units": [
        {"polarity": "negative", "raw_dir": "neg"}, {"polarity": "positive", "raw_dir": "pos"}]})

    assert receipt["units"] == ["negative", "positive"]
    resolved = receipt["resolved"]["units"]
    assert resolved["negative"]["polarity"]["value"] == "negative"
    assert resolved["positive"]["polarity"]["value"] == "positive"
    record = load_run(Path(receipt["pipeline_path"]))
    assert record["units"]["positive"]["raw_dir"] == "pos"
    assert list(record["stages"])[-1] == "receipt"


def test_request_file_downstream_keys_appear_as_a_warning_in_the_plan(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    (source["root"] / "analysis-request.json").write_text(
        json.dumps({"target": "differential", "comparisons": []}), encoding="utf-8")

    receipt = service.plan_pipeline(source["root"])

    assert [w["code"] for w in receipt["warnings"]] == ["REQUEST_FILE_DOWNSTREAM_IGNORED"]


def test_start_pipeline_saves_needs_input_without_launching_on_bad_manifest(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    (source["root"] / "sample-manifest.tsv").write_text(
        "not-a-valid-manifest-header\n", encoding="utf-8")

    def poison(*a, **kw):
        raise AssertionError("不正な入力を検出したのにworkerを起動した")

    monkeypatch.setattr(service, "launch_detached", poison)

    receipt = service.start_pipeline(source["root"])

    assert receipt["launched"] is False
    assert receipt["status"] == "needs_input"
    assert receipt["needs_input"]["code"] == "SAMPLE_MANIFEST_INVALID"
    record = load_run(Path(receipt["pipeline_path"]))
    assert record["status"] == "needs_input"
    assert record["needs_input"]["stage_id"] == "prepare_input:negative"


def test_resume_of_a_needs_input_run_rechecks_the_manifest_before_launching(tmp_path, monkeypatch):
    """受付の事前検査で止まった run は、resume でもシートを検査し直す。直っていなければ起動しない。"""
    source = _prepare_source(tmp_path, monkeypatch)
    manifest = source["root"] / "sample-manifest.tsv"
    manifest.write_text("not-a-valid-manifest-header\n", encoding="utf-8")

    def poison(*a, **kw):
        raise AssertionError("壊れたシートのまま worker を起動した")

    monkeypatch.setattr(service, "launch_detached", poison)
    pipeline_path = Path(service.start_pipeline(source["root"])["pipeline_path"])

    still_broken = service.resume_pipeline(pipeline_path)

    assert still_broken["status"] == "needs_input"
    assert still_broken["launch"] == {"launched": False}
    assert still_broken["needs_input"]["code"] == "SAMPLE_MANIFEST_INVALID"
    record = load_run(pipeline_path)
    assert record["status"] == "needs_input"
    assert record["needs_input"]["stage_id"] == "prepare_input:negative"

    manifest.unlink()  # シートを外せば検査は通る
    launch_calls = _counting_launch(monkeypatch)
    fixed = service.resume_pipeline(pipeline_path)
    assert fixed["status"] == "planned"
    assert len(launch_calls) == 1
    assert load_run(pipeline_path)["needs_input"] is None


def test_start_pipeline_handshake_timeout_does_not_relaunch(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls = _counting_launch(monkeypatch, pid_base=700)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.15)

    receipt = service.start_pipeline(source["root"])

    assert receipt["launch"]["handshake"] == "not_confirmed"
    assert len(launch_calls) == 1
    assert receipt["status"] == "planned"


# ---------- handler 一式が本物として連鎖すること ----------

def _use_fake_console_script(monkeypatch, scenario: str = "success") -> None:
    """Console のコマンドライン組み立てだけを偽 Console（tests/fixtures/fake_console.py）へ。"""
    def build_msdial_cmd(exe, dataset_root, msdial_out_dir, method_file, save_project=False):
        return [sys.executable, str(_FAKE_CONSOLE), "--scenario", scenario,
                "-i", str(dataset_root), "-o", str(msdial_out_dir), "-m", str(method_file)]

    monkeypatch.setattr("metabolomix.console.runner.build_msdial_cmd", build_msdial_cmd)


def test_worker_run_worker_completes_a_real_pipeline(tmp_path, monkeypatch):
    """`worker.run_worker` を実際に呼び、prepare_input の実配置・upstream の実監視・
    validate_outputs（mzTab-M と解析フォルダの `.arf` の照合）・受領を本物で通す。"""
    from metabolomix.pipeline import worker

    source = _prepare_source(tmp_path, monkeypatch)
    plan = service.plan_pipeline(source["root"], {"save_project": False})
    pipeline_path = Path(plan["pipeline_path"])
    _use_fake_console_script(monkeypatch)

    result = worker.run_worker(pipeline_path)

    assert result["status"] == "completed", result.get("warnings")
    assert {sid: s["status"] for sid, s in result["stages"].items()} == {
        "prepare_input:negative": "succeeded", "upstream:negative": "succeeded",
        "validate_outputs:negative": "succeeded", "receipt": "succeeded"}
    unit = result["units"]["negative"]
    assert unit["outputs"]["alignment"]["arf"].endswith("PeakProperties.arf")
    assert (pipeline_path / "negative" / "input" / unit["outputs"]["alignment"]["arf"]).is_file()
    receipt = json.loads((pipeline_path / "results" / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["units"]["negative"]["status"] == "validated"


# ---------- 研究室ライブラリの置き場所を戻り値に出さない ----------

def test_plan_receipt_shows_the_lab_msp_without_its_location(tmp_path, monkeypatch):
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_only"})
    resolved = receipt["resolved"]["units"]["negative"]
    assert resolved["library_mode"] == "msp_only"
    assert resolved["lbm"] == {"path": None, "sha256": None}
    assert resolved["msp"]["file"] == "lab_neg.msp"
    assert len(resolved["msp"]["sha256"]) == 64
    assert resolved["msp"]["source"] == "config_file"
    assert resolved["dot_product_cutoffs"]["scale"] == "square_root"
    assert "secret_share" not in json.dumps(receipt, ensure_ascii=False)


def test_plan_receipt_shows_lbm_without_location(tmp_path, monkeypatch):
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share", lbm=True)
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_then_lbm"})
    resolved = receipt["resolved"]["units"]["negative"]
    assert resolved["library_mode"] == "msp_then_lbm"
    assert set(resolved["lbm"]) == {"file", "sha256", "source"}
    assert resolved["lbm"]["file"] == "lab_lipids.lbm2"
    assert resolved["lbm"]["source"] == "config_file"
    assert resolved["msp"]["file"] == "lab_neg.msp"
    assert resolved["msp"]["priority"] == 2
    cutoffs = resolved["dot_product_cutoffs"]
    assert (cutoffs["weighted"], cutoffs["simple"], cutoffs["reverse"]) == (0.6, 0.6, 0.8)
    assert "secret_share" not in json.dumps(receipt, ensure_ascii=False)


def test_plan_fingerprint_follows_the_lbm_of_msp_then_lbm():
    plan = {
        "raw_stat": [{"relative_path": "a.wiff", "size": 1, "mtime_ns": 2}],
        "selected_format": "wiff",
        "method": {"sha256": "m" * 64}, "lbm": {"sha256": "l" * 64},
        "exe": {"sha256": "e" * 64}, "polarity": {"value": "negative"},
        "msp": {"sha256": "s" * 64},
    }
    other = {**plan, "lbm": {"sha256": "x" * 64}}
    assert service._plan_fingerprint(plan) != service._plan_fingerprint(other)


def test_plan_fingerprint_is_unchanged_without_an_msp():
    from metabolomix.core.atomic_io import canonical_hash
    plan = {
        "raw_stat": [{"relative_path": "a.wiff", "size": 1, "mtime_ns": 2}],
        "selected_format": "wiff",
        "method": {"sha256": "m" * 64}, "lbm": {"sha256": "l" * 64},
        "exe": {"sha256": "e" * 64}, "polarity": {"value": "negative"},
    }
    legacy = canonical_hash({
        "raw_stat": plan["raw_stat"], "selected_format": "wiff",
        "method_sha256": "m" * 64, "lbm_sha256": "l" * 64, "exe_sha256": "e" * 64,
        "polarity": "negative"})
    assert service._plan_fingerprint(plan) == legacy
    plan["msp"] = {"path": None, "file": None, "sha256": None, "source": "not_used"}
    assert service._plan_fingerprint(plan) == legacy
    plan["msp"] = {"path": "x", "file": "x.msp", "sha256": "s" * 64, "source": "config_file"}
    assert service._plan_fingerprint(plan) != legacy


def test_status_details_do_not_reveal_the_lab_msp_location(tmp_path, monkeypatch):
    from metabolomix.pipeline import recovery, store
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_only"})
    status = recovery.read_status(Path(receipt["pipeline_path"]), include_details=True)
    assert "secret_share" not in json.dumps(status, ensure_ascii=False)
    inputs = status["record"]["units"]["negative"]["inputs"]
    assert inputs["msp"]["file"] == "lab_neg.msp"
    assert inputs["method"]["overrides"]["Msp file path"] == "lab_neg.msp"
    # データ側の記録そのものは絶対パスのまま（resume の再検査に要る）。
    saved = store.load_run(Path(receipt["pipeline_path"]))
    assert "secret_share" in saved["units"]["negative"]["inputs"]["msp"]["path"]


def test_status_details_do_not_reveal_the_lab_lbm_location(tmp_path, monkeypatch):
    from metabolomix.pipeline import recovery, store
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share", lbm=True)
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_then_lbm"})
    saved = store.load_run(Path(receipt["pipeline_path"]))
    assert "secret_share" in saved["units"]["negative"]["inputs"]["lbm"]["path"]

    status = recovery.read_status(Path(receipt["pipeline_path"]), include_details=True)
    assert "secret_share" not in json.dumps(status, ensure_ascii=False)
    inputs = status["record"]["units"]["negative"]["inputs"]
    assert inputs["lbm"]["path"] == "lab_lipids.lbm2"
    assert inputs["method"]["overrides"]["Lbm file path"] == "lab_lipids.lbm2"
    assert inputs["method"]["overrides"]["Msp file path"] == ""


def _auto_source_with_library_lbm(tmp_path, monkeypatch):
    """auto で LBM が `[library] lbm`（出どころ config_file）から解ける生データフォルダ。"""
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    source["method"].write_text(
        "Ion mode: negative\nTarget omics: Lipidomics\nLbm file path: \n",
        encoding="ascii", newline="\n")
    paths = write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share", lbm=True)
    return source, paths


def test_auto_plan_shows_only_the_name_of_a_library_lbm(tmp_path, monkeypatch):
    from metabolomix.pipeline import recovery, store
    source, paths = _auto_source_with_library_lbm(tmp_path, monkeypatch)
    receipt = service.plan_pipeline(source["root"], {"library_mode": "auto"})
    assert "secret_share" not in json.dumps(receipt, ensure_ascii=False)
    resolved = receipt["resolved"]["units"]["negative"]
    assert resolved["lbm"]["path"] == "lab_lipids.lbm2"
    assert len(resolved["lbm"]["sha256"]) == 64

    saved = store.load_run(Path(receipt["pipeline_path"]))
    inputs = saved["units"]["negative"]["inputs"]
    assert inputs["lbm"]["path"] == str(paths["lbm"].resolve())
    assert inputs["lbm"]["source"] == "config_file"

    status = recovery.read_status(Path(receipt["pipeline_path"]), include_details=True)
    assert "secret_share" not in json.dumps(status, ensure_ascii=False)
    assert status["record"]["units"]["negative"]["inputs"]["lbm"]["path"] == "lab_lipids.lbm2"


def test_auto_plan_keeps_the_full_path_of_a_declared_lbm(tmp_path, monkeypatch):
    source = _prepare_source(tmp_path, monkeypatch)
    receipt = service.plan_pipeline(source["root"], {"library_mode": "auto"})
    assert receipt["resolved"]["units"]["negative"]["lbm"]["path"] == \
        str(source["lbm"].resolve())


def test_verify_inputs_names_only_the_file_of_a_changed_library_lbm(tmp_path):
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import inputs as inputs_mod
    lib = tmp_path / "secret_share" / "lab_lipids.lbm2"
    lib.parent.mkdir()
    lib.write_bytes(b"changed")
    snapshot = {"source_root": str(tmp_path), "raw_stat": [], "library_mode": "auto",
                "lbm": {"path": str(lib), "sha256": "0" * 64, "source": "config_file"}}
    with pytest.raises(DomainError) as exc:
        inputs_mod.verify_inputs(snapshot)
    assert exc.value.code == "INPUT_CHANGED"
    assert "secret_share" not in exc.value.message
    assert "secret_share" not in json.dumps(exc.value.details, ensure_ascii=False)
    assert exc.value.details["lbm_file"] == "lab_lipids.lbm2"


def test_without_library_location_redacts_absolute_annotator_tables():
    from metabolomix.pipeline import recovery
    inputs = {
        "library_mode": "msp_then_lbm",
        "lbm": {"path": "//nas/secret_share/lab.lbm2", "file": "lab.lbm2"},
        "msp": {"path": "//nas/secret_share/lab_neg.msp", "file": "lab_neg.msp"},
        "method": {"overrides": {
            "Lbm file path": "//nas/secret_share/lab.lbm2",
            "MSP annotator settings file path": "//nas/secret_share/msp-annotators.tsv",
            "MSP search settings file path": "inputs/msp-annotators.tsv",
        }},
    }
    record = {"schema": "pipeline-run.v3", "units": {"negative": {"inputs": inputs}}}
    redacted = recovery._without_library_location(record)
    assert "secret_share" not in json.dumps(redacted)
    overrides = redacted["units"]["negative"]["inputs"]["method"]["overrides"]
    assert overrides["MSP annotator settings file path"] == "msp-annotators.tsv"
    assert overrides["MSP search settings file path"] == "inputs/msp-annotators.tsv"
    assert inputs["lbm"]["path"] == "//nas/secret_share/lab.lbm2"  # 元は変えない
    # auto の LBM は従来どおり path を残す。
    auto = {"schema": "pipeline-run.v3", "units": {"negative": {"inputs": {
        "library_mode": "auto", "lbm": {"path": "C:/MS-DIAL/x.lbm2"}, "msp": {"path": None},
        "method": {"overrides": {"Lbm file path": "C:/MS-DIAL/x.lbm2"}}}}}}
    assert recovery._without_library_location(auto)["units"]["negative"]["inputs"]["lbm"][
        "path"] == "C:/MS-DIAL/x.lbm2"


def test_plan_receipt_shows_alignment_reference(tmp_path, monkeypatch):
    from tests.pipeline_fixtures import write_sample_manifest
    source = _prepare_source(tmp_path, monkeypatch)
    names = [f"S{i}.wiff" for i in range(8)]
    sheet = write_sample_manifest(source["root"] / "sample-manifest.tsv", source_names=names)
    # 注入順をファイル名の逆にする（S7 が 1 番目）。
    text = sheet.read_text(encoding="utf-8").splitlines()
    for i in range(8):
        cells = text[2 + i].split("\t")
        cells[5] = str(8 - i)
        text[2 + i] = "\t".join(cells)
    sheet.write_text("\n".join(text) + "\n", encoding="utf-8", newline="\n")

    receipt = service.plan_pipeline(source["root"])
    reference = receipt["resolved"]["units"]["negative"]["alignment_reference"]
    assert reference["injection_order_source"] == "manifest"
    assert (reference["file_id"], reference["sample_file"]) == (4, "S4.wiff")

    method_receipt = service.plan_pipeline(source["root"], {"alignment_reference": "method"})
    assert method_receipt["resolved"]["units"]["negative"]["alignment_reference"][
        "file_id"] is None
