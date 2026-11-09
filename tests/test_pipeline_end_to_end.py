"""偽 Console による一気通貫の受入検証（pipeline-run.v3: 上流だけ）。

ここでの偽物は**worker が起動する Console のコマンドだけ**。受付・engine・入力の配置・
Console の監視・mzTab-M と解析フォルダ（`.arf`）の照合・受領はすべて本物を通す。
偽 Console（`tests/fixtures/fake_console.py`）は `-o` へ mzTab-M を、`-i`（解析フォルダ
＝`<極性>/input/`）へ同じ特徴・試料名の `.arf` / `.arf2` を書く。
"""
from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from metabolomix.pipeline import recovery, store
from tests.pipeline_fixtures import PipelineHarness, WAIT_TIMEOUT_S


@pytest.fixture
def pipeline_harness(tmp_path, monkeypatch):
    harness = PipelineHarness(tmp_path, monkeypatch)
    try:
        yield harness
    finally:
        harness.close()


def _arf_files(folder: Path) -> list[str]:
    return sorted(p.name for p in folder.glob("*.arf"))


# ---------- 受領で終わる ----------

def test_single_polarity_run_ends_with_receipt(pipeline_harness):
    run = pipeline_harness.start()
    record = pipeline_harness.wait(run, expected="completed")
    root = pipeline_harness.pipeline_root(run)

    assert record["schema"] == "pipeline-run.v3"
    assert list(record["stages"]) == [
        "prepare_input:negative", "upstream:negative", "validate_outputs:negative", "receipt"]
    assert all(s["status"] == "succeeded" for s in record["stages"].values())
    assert pipeline_harness.console_start_count(run) == 1

    unit = record["units"]["negative"]
    assert unit["upstream"]["verification"]["status"] == "completed"
    assert unit["upstream"]["console_job_path"].startswith("negative/console/attempt-0001/")
    analysis_dir = root / "negative" / "input"
    assert unit["outputs"]["alignment"]["arf"] in _arf_files(analysis_dir)
    assert (analysis_dir / unit["outputs"]["alignment"]["arf2"]).is_file()
    assert (root / unit["outputs"]["mztab"]).is_file()

    receipt = pipeline_harness.receipt(run)
    assert receipt["units"]["negative"]["status"] == "validated"
    assert receipt["units"]["negative"]["analysis_dir"] == "negative/input"
    assert receipt["units"]["negative"]["alignment"] == unit["outputs"]["alignment"]
    assert receipt["units"]["negative"]["samples"]["n_raw"] == 8
    assert receipt["warnings"] == []
    assert receipt["next"] == \
        'load_dataset(pipeline_path="<pipeline_path>", polarity="negative")'

    status = recovery.read_status(root)
    assert status["units"] == {"negative": "validated"}
    assert status["receipt"]["next"] == \
        f'load_dataset(pipeline_path="{str(root).replace(chr(92), "/")}", polarity="negative")'


def test_two_polarity_run_has_two_analysis_folders(pipeline_harness):
    run = pipeline_harness.start(two_polarity=True)
    record = pipeline_harness.wait(run, expected="completed")
    root = pipeline_harness.pipeline_root(run)

    # 極性ごとに Console を 1 回ずつ、要求の順に回す。
    assert pipeline_harness.console_start_count(run) == 2
    assert list(record["stages"]) == [
        "prepare_input:negative", "upstream:negative", "validate_outputs:negative",
        "prepare_input:positive", "upstream:positive", "validate_outputs:positive",
        "receipt"]
    for polarity, raw_dir in (("negative", "neg"), ("positive", "pos")):
        unit = record["units"][polarity]
        assert unit["raw_dir"] == raw_dir
        analysis_dir = root / polarity / "input"
        assert sorted(p.name for p in analysis_dir.glob("S*.wiff")) == \
            [f"S{i}.wiff" for i in range(8)]
        assert unit["outputs"]["alignment"]["arf"] in _arf_files(analysis_dir)
        assert pipeline_harness.console_receipt(run, polarity=polarity)["termination"] == \
            "exited"

    receipt = pipeline_harness.receipt(run)
    assert {p: u["analysis_dir"] for p, u in receipt["units"].items()} == {
        "negative": "negative/input", "positive": "positive/input"}
    assert 'polarity="negative" または "positive"' in receipt["next"]


def test_unmatched_arf_is_partial_with_warning(pipeline_harness):
    """解析フォルダの `.arf` が mzTab-M と合わなければ、単位は validated のまま
    `alignment: null` と警告 `ANALYSIS_FOLDER_UNMATCHED`、run は `partial`。"""
    run = pipeline_harness.start(scenario="unmatched_arf")
    record = pipeline_harness.wait(run, expected="partial")

    unit = record["units"]["negative"]
    assert record["stages"]["validate_outputs:negative"]["status"] == "succeeded"
    assert unit["outputs"]["alignment"] is None
    assert unit["outputs"]["mztab"]  # mzTab-M（source="mztab"）なら読める
    codes = [w["code"] for w in record["stages"]["validate_outputs:negative"]["warnings"]]
    assert codes == ["ANALYSIS_FOLDER_UNMATCHED"]
    warning = record["stages"]["validate_outputs:negative"]["warnings"][0]
    assert "source=\"mztab\"" in warning["message"]
    assert "ANALYSIS_FOLDER_UNMATCHED" in {w["code"] for w in record["warnings"]}

    receipt = pipeline_harness.receipt(run)
    assert receipt["units"]["negative"]["status"] == "validated"
    assert receipt["units"]["negative"]["alignment"] is None
    assert "ANALYSIS_FOLDER_UNMATCHED" in {
        w["code"] for w in receipt["units"]["negative"]["validation"]["warnings"]}


def test_receipt_names_library_files_only(pipeline_harness, tmp_path, monkeypatch):
    """研究室の MSP と `[library] lbm` の LBM は、受領にも `pipeline_status` にも
    ファイル名と sha256 だけを示す（保存する run の記録は絶対パスのまま）。"""
    from metabolomix.tools.pipeline_tools import pipeline_status
    from tests.lab_msp_fixtures import write_lab_msp_config

    pipeline_harness._configure()
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share", lbm=True)
    run = pipeline_harness.start(extra_request={"library_mode": "msp_then_lbm"})
    record = pipeline_harness.wait(run, expected="completed")
    root = pipeline_harness.pipeline_root(run)

    receipt_text = (root / "results" / "receipt.json").read_text(encoding="utf-8")
    assert "secret_share" not in receipt_text
    library = json.loads(receipt_text)["units"]["negative"]["library"]
    assert library["mode"] == "msp_then_lbm"
    assert library["msp"]["file"] == "lab_neg.msp"
    assert library["lbm"]["file"] == "lab_lipids.lbm2"
    assert len(library["lbm"]["sha256"]) == 64

    for include_details in (False, True):
        text = pipeline_status(str(root), include_details=include_details)
        assert "secret_share" not in text, include_details
    # 保存する記録は resume の再検査のために絶対パスを持つ。
    assert "secret_share" in record["units"]["negative"]["inputs"]["lbm"]["path"]


# ---------- Console の終了と検証 ----------

@pytest.mark.parametrize("scenario, termination, exit_code, error_codes, retained", [
    ("nonzero", "exited", 1, [], "intermediate.pai2"),
    ("invalid", "exited", 0, ["MZTAB_STRUCTURE_INVALID"], None),
    ("missing_sample", "exited", 0, ["SAMPLE_MAPPING_MISSING"], None),
    ("hang", "timeout", None, [], "intermediate.pai2"),
])
def test_execution_scenarios(pipeline_harness, scenario, termination, exit_code,
                             error_codes, retained):
    """異常終了・不正 mzTab・assay 欠落・timeout は単位の失敗で、受領は走らず run は failed。"""
    run = pipeline_harness.start_scenario(scenario, timeout_s=3)
    record = pipeline_harness.wait(run, expected="failed")

    receipt = pipeline_harness.console_receipt(run)
    assert receipt["termination"] == termination
    if exit_code is None:
        assert receipt["exit_code"] != 0, "timeoutを0で補完してはいけない"
    else:
        assert receipt["exit_code"] == exit_code
    for code in error_codes:
        assert code in receipt["validation"]["errors"], receipt["validation"]

    # 「終了検証が失敗し、成果物は保持された」——ファイルの有無を成功と読み替えない。
    assert record["stages"]["upstream:negative"]["status"] == "failed"
    assert record["stages"]["upstream:negative"]["error"]["code"] == "MSDIAL_EXECUTION_FAILED"
    assert record["units"]["negative"]["upstream"]["verification"]["status"] != "completed"
    assert record["stages"]["validate_outputs:negative"]["status"] == "pending"
    assert record["stages"]["receipt"]["status"] == "pending"
    assert receipt["collection"]["status"] == "succeeded"
    if retained:
        artifacts = pipeline_harness.console_job(run).artifacts
        assert [a for a in artifacts if a.path.endswith(retained)], \
            f"{retained} が保持されていません: {[a.path for a in artifacts]}"


def test_rerun_upstream_creates_a_new_attempt_and_keeps_the_old_evidence(pipeline_harness):
    """Console の再試行は明示（`rerun_upstream=true`）だけで、新しい attempt を作る。
    前回 attempt の終了証跡・ログ・実行前スナップショットは 1 バイトも変えない。"""
    from metabolomix.pipeline import service

    run = pipeline_harness.start_scenario("hang", timeout_s=3)
    pipeline_harness.wait(run, expected="failed")

    first_dir = pipeline_harness.console_run_dir(run)
    first_receipt = json.loads((first_dir / "execution-result.json").read_text(encoding="utf-8"))
    assert first_receipt["termination"] == "timeout"
    first_log = (first_dir / "msdial.log").read_bytes()
    first_state = json.loads((first_dir / "worker.json").read_text(encoding="utf-8"))
    assert first_state["recovery"]["befores"], first_state

    pipeline_harness.set_launch_options(run, scenario="success")
    service.resume_pipeline(pipeline_harness.pipeline_root(run), rerun_upstream=True)
    completed = pipeline_harness.wait(run, expected="completed")

    assert pipeline_harness.console_start_count(run) == 2
    second_dir = pipeline_harness.console_run_dir(run)
    assert second_dir != first_dir, "再実行が前回と同じrun_dirを使っています"
    job_rel = completed["units"]["negative"]["upstream"]["console_job_path"]
    assert (pipeline_harness.pipeline_root(run) / job_rel).parent == second_dir

    assert json.loads(
        (first_dir / "execution-result.json").read_text(encoding="utf-8")) == first_receipt
    assert (first_dir / "msdial.log").read_bytes() == first_log
    assert json.loads((first_dir / "worker.json").read_text(encoding="utf-8")) == first_state
    assert json.loads(
        (second_dir / "execution-result.json").read_text(encoding="utf-8"))["termination"] \
        == "exited"


def test_rerun_after_an_alignment_was_written_picks_the_new_batch(pipeline_harness):
    """前の attempt のバッチが解析フォルダに残っても、受領は今回のジョブのバッチを選ぶ
    （`dataset/opening` と同じ規則: ジョブの `peak_matrix_source` が先）。"""
    from metabolomix.pipeline import service

    run = pipeline_harness.start_scenario("arf_then_fail")
    pipeline_harness.wait(run, expected="failed")
    analysis_dir = pipeline_harness.pipeline_root(run) / "negative" / "input"
    first = _arf_files(analysis_dir)
    assert len(first) == 1

    pipeline_harness.set_launch_options(run, scenario="success")
    service.resume_pipeline(pipeline_harness.pipeline_root(run), rerun_upstream=True)
    record = pipeline_harness.wait(run, expected="completed")

    batches = _arf_files(analysis_dir)
    assert len(batches) == 2
    newer = max(batches)
    assert newer != first[0]
    assert record["units"]["negative"]["outputs"]["alignment"]["arf"] == newer
    assert "ANALYSIS_FOLDER_UNMATCHED" not in {w["code"] for w in record["warnings"]}
    assert pipeline_harness.receipt(run)["units"]["negative"]["alignment"]["arf"] == newer


def test_missing_gui_project_is_a_warning_in_the_receipt(pipeline_harness):
    """save_project=true なのに `.mdproject` が無ければ警告 `GUI_PROJECT_UNAVAILABLE`。
    解析フォルダは確かめられているので単位は validated のまま。"""
    run = pipeline_harness.start(save_project=True)
    record = pipeline_harness.wait(run, expected="completed")

    assert "GUI_PROJECT_UNAVAILABLE" in {w["code"] for w in record["warnings"]}
    assert record["units"]["negative"]["outputs"]["gui_project"] is None
    receipt = pipeline_harness.receipt(run)
    validation = receipt["units"]["negative"]["validation"]
    assert validation["gui_project"] is None
    assert "GUI_PROJECT_UNAVAILABLE" in {w["code"] for w in validation["warnings"]}


# ---------- ディレクトリ形式 raw ----------

def test_directory_style_raw_runs_end_to_end_without_touching_the_source(pipeline_harness):
    """1 測定単位がフォルダ（Agilent `.d` 等）でも規定どおり準備して進む。"""
    source = pipeline_harness.source("dirsource", extension=".d", kind="dir")

    def snapshot() -> dict:
        return {str(path.relative_to(source["root"])): path.stat().st_mtime_ns
                for path in sorted(Path(source["root"]).rglob("*"))
                if "runs" not in path.relative_to(source["root"]).parts}

    before = snapshot()
    run = pipeline_harness.start(source_name="dirsource")
    record = pipeline_harness.wait(run, expected="completed")

    entries = [e for e in record["units"]["negative"]["inputs"]["entries"]
               if e["role"] == "primary"]
    assert [e["kind"] for e in entries] == ["dir"] * 8
    staged = pipeline_harness.pipeline_root(run) / "negative" / "input"
    assert sorted(p.name for p in staged.iterdir() if p.suffix == ".d") == \
        [f"S{i}.d" for i in range(8)]
    assert all((staged / f"S{i}.d" / "AcqData.bin").is_file() for i in range(8))
    assert snapshot() == before


# ---------- 同一要求の再送・同時送信 ----------

def test_identical_request_resend_reuses_the_run_and_starts_console_once(pipeline_harness):
    first = pipeline_harness.start(request_id="req-1")
    record = pipeline_harness.wait(first, expected="completed")

    second = pipeline_harness.start(request_id="req-1")
    third = pipeline_harness.start()  # request_id無しの再送

    assert second == first and third == first
    assert pipeline_harness.console_start_count(first) == 1
    assert pipeline_harness.record(first)["status"] == "completed"
    runs_parent = Path(record["identity"]["pipeline_root"]).parent
    assert len(list(runs_parent.iterdir())) == 1


def test_two_concurrent_starts_of_the_same_request_launch_console_once(pipeline_harness):
    """同一要求の `start_pipeline` が同時に来ても、run は 1 つ・Console は 1 回だけ。

    受付直後にバリアを挟み、2 つの呼び出しが**どちらも受付を終えてから**起動判定へ
    進む瞬間を強制する。二重起動を止めるのは worker 側の owner lock。
    """
    barrier = threading.Barrier(2, timeout=WAIT_TIMEOUT_S)
    real_find_or_create_run = store.find_or_create_run

    def _synced_find_or_create_run(*args, **kwargs):
        result = real_find_or_create_run(*args, **kwargs)
        barrier.wait()
        return result

    pipeline_harness.monkeypatch.setattr(
        store, "find_or_create_run", _synced_find_or_create_run)
    pipeline_harness.source("source")

    results: list = []
    errors: list = []

    def _call():
        try:
            results.append(pipeline_harness.start(request_id="req-simul"))
        except Exception as exc:  # pragma: no cover - 失敗時の診断用
            errors.append(exc)

    threads = [threading.Thread(target=_call) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=WAIT_TIMEOUT_S)

    assert not errors, f"同時起動でエラー（バリアが揃わなかった可能性）: {errors}"
    first, second = results
    assert first == second, "2つのstart_pipeline呼び出しが別のrunを指しました"
    workers = pipeline_harness.worker_processes(first)
    assert len(workers) == 2, "2つの呼び出しが2つのworker起動を試みていません"

    record = pipeline_harness.wait(first, expected="completed")
    assert pipeline_harness.console_start_count(first) == 1
    runs_parent = Path(record["identity"]["pipeline_root"]).parent
    assert len(list(runs_parent.iterdir())) == 1
    assert {w.returncode for w in workers} <= {0, 1}


def test_two_concurrent_resumes_of_the_same_run_launch_console_once(pipeline_harness):
    """完了した run への試料の対応表の訂正が同時に来ても、Console は 1 回のまま。"""
    run = pipeline_harness.start()
    pipeline_harness.wait(run, expected="completed")
    manifest = pipeline_harness.write_manifest()

    barrier = threading.Barrier(2, timeout=WAIT_TIMEOUT_S)
    real_prepare_resume = recovery.prepare_resume

    def _synced_prepare_resume(*args, **kwargs):
        result = real_prepare_resume(*args, **kwargs)
        barrier.wait()
        return result

    pipeline_harness.monkeypatch.setattr(recovery, "prepare_resume", _synced_prepare_resume)

    pipeline_root = pipeline_harness.pipeline_root(run)
    results: list = []
    errors: list = []

    def _call():
        from metabolomix.pipeline import service
        try:
            results.append(service.resume_pipeline(
                pipeline_root, updates={"sample_manifest": manifest.name}))
        except Exception as exc:  # pragma: no cover - 失敗時の診断用
            errors.append(exc)

    threads = [threading.Thread(target=_call) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=WAIT_TIMEOUT_S)

    assert not errors, f"同時resumeでエラー: {errors}"
    assert len(results) == 2
    record = pipeline_harness.wait(run, expected="completed")
    assert pipeline_harness.console_start_count(run) == 1
    assert record["request"]["revision"] == 2
    receipt = pipeline_harness.receipt(run)
    assert receipt["request_revision"] == 2
    assert receipt["units"]["negative"]["samples"]["manifest"]["path"] == manifest.name
    assert {w.returncode for w in pipeline_harness.worker_processes(run)[1:]} <= {0, 1}


# ---------- 入力・メソッドの変化 ----------

def test_changed_raw_is_detected_on_resume(pipeline_harness):
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import service

    run = pipeline_harness.start()
    pipeline_harness.wait(run, expected="completed")
    source_root = Path(pipeline_harness.source("source")["root"])
    (source_root / "S3.wiff").write_text("tampered-raw-content", encoding="ascii")

    with pytest.raises(DomainError) as excinfo:
        service.resume_pipeline(pipeline_harness.pipeline_root(run),
                                updates={"sample_manifest": "sample-manifest.tsv"})
    assert excinfo.value.code == "INPUT_CHANGED"
    assert "S3.wiff" in excinfo.value.details["changed"]


def test_changed_effective_method_is_detected_on_resume(pipeline_harness):
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.pipeline import service

    run = pipeline_harness.start()
    record = pipeline_harness.wait(run, expected="completed")

    unit_root = pipeline_harness.pipeline_root(run) / "negative"
    effective = unit_root / record["units"]["negative"]["inputs"]["method"][
        "effective_relative_path"]
    effective.write_text(effective.read_text(encoding="ascii") + "Extra: 1\n", encoding="ascii")

    with pytest.raises(DomainError) as excinfo:
        service.resume_pipeline(pipeline_harness.pipeline_root(run),
                                updates={"sample_manifest": "sample-manifest.tsv"})
    assert excinfo.value.code == "STAGED_INPUT_MISMATCH"
