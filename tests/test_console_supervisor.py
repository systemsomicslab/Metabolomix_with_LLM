"""Console 監視（metabolomix.console.execution.supervise）の全終了経路の検証。

旧実装は非ゼロ終了で即 return し、生成物の収集を丸ごと飛ばしていた。**失敗した実行が
何の証拠も残さない**のがこのテストが防ぐ欠陥なので、成功・非ゼロ・timeout・取消・
起動失敗のどれでも終了証跡が残ることを、モックではなく実際の子プロセス
（`tests/fixtures/fake_console.py`）で確かめる。

`hang` シナリオは孫の実 Python まで起こす。親 pid だけを殺す実装でも「親は死んだ」
テストは通ってしまうため、Job Object が一族を停止したことを孫の identity で確認する。
"""
from __future__ import annotations

import contextlib
import json
import os
import signal
import sys
import threading
import time
from pathlib import Path

import pytest

from metabolomix.console.execution import supervise, write_supervision_inputs
from metabolomix.console.job_manager import create_job, load_job
from metabolomix.core.process_control import process_identity, same_process
from metabolomix.handoff.schema import sha256_file

pytestmark = pytest.mark.skipif(
    os.name != "nt", reason="Job Object による実プロセス監視は Windows 専用")

_FAKE_CONSOLE = Path(__file__).parent / "fixtures" / "fake_console.py"

#: 合成 raw の本数（brief 指定）。
_RAW_COUNT = 4


# ---------- fixture ----------

@pytest.fixture
def planned_fake_job(tmp_path):
    """planned 状態の実ジョブと、偽 Console の起動回数を数えるカウンタを返す。

    実 raw も既存の解析成果物も使わない。`tmp_path/source` に合成 raw を 4 本作り、
    既存の `create_job` でジョブを作ってから、監視入力（input inventory・method /
    exe の hash）を run_dir へ固定する。
    """
    source = tmp_path / "source"
    source.mkdir()
    raws = [source / f"S{i}.abf" for i in range(1, _RAW_COUNT + 1)]
    for i, raw in enumerate(raws, start=1):
        raw.write_bytes(b"\x00" * (32 + i))
    method = source / "method.txt"
    method.write_text("Ion mode: Positive\n", encoding="ascii")
    exe = tmp_path / "fake-msdialcui.exe"
    exe.write_bytes(b"fake console executable")

    job, job_path = create_job(source, method, "positive", "peak_height",
                               input_count=len(raws))
    job.timeout_s = 60
    job.save(job_path)

    write_supervision_inputs(Path(job.run_dir), {
        "raw_inventory": [str(raw.resolve()) for raw in raws],
        "method_sha256": sha256_file(method),
        "exe_path": str(exe),
        "exe_sha256": sha256_file(exe),
    })
    return job_path, tmp_path / "launches.tsv"


# ---------- helper ----------

def _fake_command(job_path: Path, counter: Path, scenario: str) -> list[str]:
    return [sys.executable, str(_FAKE_CONSOLE), "--scenario", scenario,
            "--job", str(job_path), "--counter", str(counter)]


def _set_timeout(job_path: Path, timeout_s: int) -> None:
    job = load_job(job_path)
    job.timeout_s = timeout_s
    job.save(job_path)


def _raw_files(job_path: Path) -> list[Path]:
    dataset_root = Path(load_job(job_path).dataset_root)
    return sorted(p for p in dataset_root.iterdir() if p.suffix == ".abf")


def _receipt_on_disk(job_path: Path) -> dict:
    path = Path(load_job(job_path).run_dir) / "execution-result.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _supervision_state(job_path: Path) -> dict:
    path = Path(load_job(job_path).run_dir) / "worker.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _launch_lines(counter: Path) -> list[list[str]]:
    if not counter.exists():
        return []
    return [line.split("\t") for line in
            counter.read_text(encoding="utf-8").splitlines() if line.strip()]


def _wait_until_gone(identity: dict | None, timeout: float = 10.0) -> bool:
    if identity is None:
        return True
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not same_process(identity):
            return True
        time.sleep(0.02)
    return not same_process(identity)


def _kill_pid_if_same(identity: dict | None) -> None:
    """後始末。pid 再利用を踏まないよう identity 一致時だけ終了させる。"""
    if not identity or not same_process(identity):
        return
    try:
        os.kill(int(identity["pid"]), signal.SIGTERM)
    except OSError:
        pass


@contextlib.contextmanager
def _watch_grandchild(counter: Path, *, on_seen=None):
    """supervise が走っている間に、孫プロセスの identity を捕まえておく見張り。

    supervise から戻ってから孫を探しても、既に停止済みで identity が取れない。
    それでは「孫が残っていない」の検査が **孫が起動しなかった場合にも通って
    しまう**（本当に一族を停止したのか分からない）。生きているうちに捕まえて
    おき、停止したことを identity 一致で確かめる。

    `on_seen` は孫の起動を確認した直後に呼ぶ。取消要求を「孫が確実に居る状態」で
    置くために使う。
    """
    seen: dict = {}
    stop = threading.Event()

    def _run() -> None:
        while not stop.is_set():
            for row in _launch_lines(counter):
                if len(row) >= 3 and row[2] not in ("", "-"):
                    identity = process_identity(int(row[2]))
                    if identity is not None:
                        seen["identity"] = identity
                        if on_seen is not None:
                            on_seen()
                        return
            stop.wait(0.02)

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    try:
        yield seen
    finally:
        stop.set()
        thread.join(timeout=5)
        _kill_pid_if_same(seen.get("identity"))


# ---------- 代表テスト（Task 4 brief） ----------

def test_nonzero_execution_collects_before_return(planned_fake_job):
    """非ゼロ終了でも収集まで進む。旧経路はここで即 return して成果物を捨てていた。"""
    job_path, counter = planned_fake_job
    command = _fake_command(job_path, counter, "nonzero")

    receipt = supervise(job_path, command=command)

    assert receipt["exit_code"] == 1
    assert receipt["collection"]["status"] == "succeeded"
    assert receipt["validation"]["status"] != "succeeded"


# ---------- 終了経路ごとの証跡 ----------

def test_success_execution_completes_job(planned_fake_job):
    job_path, counter = planned_fake_job

    receipt = supervise(job_path, command=_fake_command(job_path, counter, "success"))

    assert receipt["termination"] == "exited"
    assert receipt["exit_code"] == 0
    assert receipt["validation"]["status"] == "succeeded"
    assert receipt["validation"]["errors"] == []
    assert load_job(job_path).status == "completed"
    assert len(_launch_lines(counter)) == 1, "偽 Console は 1 回だけ起動するはず"


def test_nonzero_execution_keeps_partial_artifacts_on_job(planned_fake_job):
    job_path, counter = planned_fake_job

    supervise(job_path, command=_fake_command(job_path, counter, "nonzero"))

    job = load_job(job_path)
    assert job.status == "partial"
    assert [a.path for a in job.artifacts if a.path.endswith("intermediate.pai2")]


def test_timeout_terminates_process_tree_and_records_receipt(planned_fake_job):
    """timeout で親も孫も止め、証跡を残す。親だけ殺す実装ならここで落ちる。"""
    job_path, counter = planned_fake_job
    _set_timeout(job_path, 2)

    with _watch_grandchild(counter) as seen:
        started = time.monotonic()
        receipt = supervise(job_path, command=_fake_command(job_path, counter, "hang"))
        elapsed = time.monotonic() - started

        assert receipt["termination"] == "timeout"
        assert receipt["exit_code"] != 0
        assert receipt["timeout_s"] == 2
        assert 2.0 <= elapsed < 60.0, f"timeout の効き方が想定外です: {elapsed}s"
        assert not same_process(receipt["process_identity"]), "偽 Console が残っている"
        assert seen.get("identity") is not None, "孫プロセスを捕まえられなかった"
        assert _wait_until_gone(seen["identity"]), "孫プロセスが残っている"
        assert _receipt_on_disk(job_path)["termination"] == "timeout"
        assert "record_error" not in receipt, receipt.get("record_error")


def test_cancel_request_stops_the_run(planned_fake_job, tmp_path):
    """取消要求で一族を止める。取消は孫が確実に起動してから置く。"""
    job_path, counter = planned_fake_job
    cancel_path = tmp_path / "cancel.json"

    def _request_cancel() -> None:
        cancel_path.write_text('{"reason": "test"}', encoding="utf-8")

    with _watch_grandchild(counter, on_seen=_request_cancel) as seen:
        receipt = supervise(job_path, command=_fake_command(job_path, counter, "hang"),
                            cancel_path=cancel_path)

        assert receipt["termination"] == "cancelled"
        assert not same_process(receipt["process_identity"])
        assert seen.get("identity") is not None, "孫プロセスを捕まえられなかった"
        assert _wait_until_gone(seen["identity"]), "孫プロセスが残っている"
        assert load_job(job_path).status in ("partial", "failed")


def test_launch_failure_persists_launch_failed_receipt(planned_fake_job):
    """起動できなかった実行も証跡を残す。何も残らないと再開判断ができない。"""
    job_path, _counter = planned_fake_job
    missing = Path(load_job(job_path).run_dir) / "does-not-exist.exe"

    receipt = supervise(job_path, command=[str(missing)])

    assert receipt["termination"] == "launch_failed"
    assert receipt["exit_code"] is None
    assert receipt["pid"] is None
    assert receipt["process_identity"] is None
    assert _receipt_on_disk(job_path)["termination"] == "launch_failed"
    assert load_job(job_path).status == "failed"
    assert "record_error" not in receipt, receipt.get("record_error")


# ---------- 出力検証 ----------

def test_invalid_mztab_is_not_completed(planned_fake_job):
    job_path, counter = planned_fake_job

    receipt = supervise(job_path, command=_fake_command(job_path, counter, "invalid"))

    assert receipt["exit_code"] == 0
    assert "MZTAB_STRUCTURE_INVALID" in receipt["validation"]["errors"]
    assert load_job(job_path).status != "completed"


def test_missing_sample_is_detected(planned_fake_job):
    job_path, counter = planned_fake_job

    receipt = supervise(job_path,
                        command=_fake_command(job_path, counter, "missing_sample"))

    assert receipt["exit_code"] == 0
    assert "SAMPLE_MAPPING_MISSING" in receipt["validation"]["errors"]
    assert load_job(job_path).status != "completed"


def test_no_output_run_is_failed(planned_fake_job):
    """0 で終わっても何も出さない実行は failed。exit code だけでは完了にしない。"""
    job_path, _counter = planned_fake_job

    receipt = supervise(job_path, command=[sys.executable, "-c", "raise SystemExit(0)"])

    assert receipt["exit_code"] == 0
    assert receipt["collection"]["status"] == "succeeded"
    assert receipt["collection"]["artifacts"] == 0
    assert load_job(job_path).status == "failed"


def test_raw_input_change_during_run_is_detected(planned_fake_job):
    """実行中に元 raw が変われば INPUT_CHANGED。fingerprint を比べない実装なら通らない。"""
    job_path, _counter = planned_fake_job
    victim = _raw_files(job_path)[0]
    script = (f"from pathlib import Path;"
              f"Path(r'{victim}').open('ab').write(b'appended-during-run')")

    receipt = supervise(job_path, command=[sys.executable, "-c", script])

    assert receipt["inputs"]["changed"] == ["INPUT_CHANGED"]
    assert "INPUT_CHANGED" in receipt["validation"]["errors"]
    assert load_job(job_path).status != "completed"


def test_method_change_during_run_is_detected(planned_fake_job):
    job_path, _counter = planned_fake_job
    method = Path(load_job(job_path).method_file)
    script = (f"from pathlib import Path;"
              f"Path(r'{method}').open('ab').write(b'Ion mode: Negative\\n')")

    receipt = supervise(job_path, command=[sys.executable, "-c", script])

    assert "METHOD_CHANGED" in receipt["inputs"]["changed"]
    assert "METHOD_CHANGED" in receipt["validation"]["errors"]


# ---------- 収集・保存の失敗 ----------

def test_collection_failure_keeps_before_snapshot(planned_fake_job, monkeypatch):
    """収集に失敗しても証跡と再収集用の before スナップショットを残す。"""
    job_path, counter = planned_fake_job

    def _boom(*args, **kwargs):
        raise OSError("収集中のディスクエラー")

    monkeypatch.setattr("metabolomix.console.execution.collect_artifacts", _boom)

    receipt = supervise(job_path, command=_fake_command(job_path, counter, "success"))

    assert receipt["termination"] == "exited"
    assert receipt["collection"]["status"] == "failed"
    assert receipt["validation"]["status"] == "skipped"
    state = _supervision_state(job_path)
    assert state["recovery"]["befores"]["run_dir"] is not None
    assert _receipt_on_disk(job_path)["collection"]["status"] == "failed"


def test_job_save_failure_returns_recovery_information(planned_fake_job, monkeypatch):
    """job を保存できなくても、証跡と復旧に必要な情報を応答から失わない。"""
    job_path, counter = planned_fake_job

    def _boom(self, path):
        raise OSError("analysis-job.json を書けない")

    monkeypatch.setattr("metabolomix.handoff.schema.AnalysisJob.save", _boom)

    receipt = supervise(job_path, command=_fake_command(job_path, counter, "success"))

    assert receipt["job_save"]["status"] == "failed"
    recovery = receipt["job_save"]["recovery"]
    assert recovery["intended_status"] == "completed"
    assert Path(recovery["receipt_path"]).is_file()
    assert Path(recovery["supervision_state_path"]).is_file()
    assert _receipt_on_disk(job_path)["job_save"]["status"] == "failed"


def test_receipt_is_persisted_before_stage_results(planned_fake_job, monkeypatch):
    """証跡は収集より先に確定する。収集の途中で落ちても終了の事実は残る。"""
    job_path, counter = planned_fake_job
    seen: dict = {}

    def _capture(*args, **kwargs):
        seen["receipt"] = _receipt_on_disk(job_path)
        raise OSError("収集中のディスクエラー")

    monkeypatch.setattr("metabolomix.console.execution.collect_artifacts", _capture)

    supervise(job_path, command=_fake_command(job_path, counter, "success"))

    assert seen["receipt"]["termination"] == "exited"
    assert seen["receipt"]["exit_code"] == 0
    assert seen["receipt"]["collection"]["status"] == "pending"


def test_receipt_passes_execution_record_validation(planned_fake_job):
    from metabolomix.console.execution import validate_execution_record
    job_path, counter = planned_fake_job

    receipt = supervise(job_path, command=_fake_command(job_path, counter, "success"))

    assert validate_execution_record(_receipt_on_disk(job_path))["execution_id"] \
        == receipt["execution_id"]
