"""切り離し実行 — 40 分級の Console を 1 ツール呼び出しで抱えず、監視も失わない。

実データ 60 サンプルの実行は 44 分かかる。これを 1 つの MCP ツール呼び出しで待つと、
呼び出し元の都合でプロセスツリーごと落とされたときに生成物ごと失う（実測 2 回）。

旧経路は「起動して pid を返し、完了は `console_status` が後から拾う」形だった。
これには誰も監視していない時間があり、`console_status` が呼ばれなければ生成物は
永久にジョブへ載らず、**終了コードも終了理由も回収できない**。ファイルが増えた
という事実だけを根拠に completed へ進めていたため、途中で落ちた実行が「完了」に
化けた。

いまは `console_run(detach=True)` が切り離しワーカー
（`python -m metabolomix.console.worker --job <path>`）を起こし、そのワーカーが
job 単位の OS ロックを持ったまま `supervise` を最後まで回す。`console_status` は
保存済みの状態を読むだけになり、完了処理を担わない。
"""
from __future__ import annotations

import json as _json
import os
import sys
import time
from pathlib import Path

import pytest

from metabolomix.console.job_manager import create_job, load_job

_WINDOWS_ONLY = pytest.mark.skipif(
    os.name != "nt", reason="切り離しワーカーの実プロセス起動は Windows 専用")

_FAKE_CONSOLE = Path(__file__).parent / "fixtures" / "fake_console.py"

#: 実ワーカーの完了を待つ上限（秒）。子は Python の起動と numpy の import を
#: 挟むので、fake Console 自体の実行時間よりずっと長く見ておく。
_WORKER_TIMEOUT_S = 120


# ---------- helper ----------

def _planned(tmp_path, monkeypatch, *, raw_count: int = 2):
    """planned 状態のジョブを 1 件作る（合成 raw と method ファイル付き）。"""
    monkeypatch.setenv("MSDIAL_EXE", str(tmp_path / "fake-msdialcui.exe"))
    (tmp_path / "fake-msdialcui.exe").write_bytes(b"fake console executable")
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    source = tmp_path / "source"
    source.mkdir(exist_ok=True)
    for i in range(1, raw_count + 1):
        (source / f"S{i}.abf").write_bytes(b"\x00" * (32 + i))
    method = source / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    _, job_path = create_job(dataset_root=source, method_file=method,
                             polarity="negative", measure="peak_height")
    return job_path


def _run_dir(job_path: Path) -> Path:
    return Path(load_job(job_path).run_dir)


def _wait_for_status(job_path: Path, statuses: set[str],
                     timeout: float = _WORKER_TIMEOUT_S) -> str:
    """ジョブが終端状態になるまで待つ（console_status は呼ばない）。"""
    deadline = time.monotonic() + timeout
    status = load_job(job_path).status
    while status not in statuses and time.monotonic() < deadline:
        time.sleep(0.1)
        status = load_job(job_path).status
    return status


def _fake_console_args(job_path: Path, counter: Path, scenario: str) -> list[str]:
    return [sys.executable, str(_FAKE_CONSOLE), "--scenario", scenario,
            "--job", str(job_path), "--counter", str(counter)]


def _prime_supervision(job_path: Path) -> None:
    """監視入力（予定した raw と method / exe の hash）を固定する。

    通常は `console_run` がこれを書いてからワーカーを起こす。ワーカーを直接
    起動するテストでは、同じ前提をここで作る（入力目録が無いと
    `INPUT_INVENTORY_MISSING` で completed へ進めない）。
    """
    from metabolomix.console.execution import write_supervision_inputs
    from metabolomix.handoff.schema import sha256_file
    job = load_job(job_path)
    source = Path(job.dataset_root)
    raws = sorted(p.resolve() for p in source.iterdir() if p.suffix == ".abf")
    exe = Path(os.environ["MSDIAL_EXE"])
    write_supervision_inputs(Path(job.run_dir), {
        "raw_inventory": [str(p) for p in raws],
        "method_sha256": sha256_file(Path(job.method_file)),
        "exe_path": str(exe),
        "exe_sha256": sha256_file(exe),
    })


# ---------- 旧 detach 状態（.detached-state.json）の扱い ----------

def test_legacy_detached_state_is_not_promoted_by_files(tmp_path, monkeypatch):
    """ファイルが増えたという事実だけで completed にしてはいけない。

    旧経路が残した `.detached-state.json` には終了コードも process identity も
    無い。この状態のジョブは「どう終わったか分からない」のであって、
    「成功した」ではない。
    """
    from metabolomix.tools.console_tools import console_status
    root = tmp_path / "source"
    root.mkdir()
    method = root / "method.txt"
    method.write_text("Ion mode: Positive\n", encoding="ascii")
    job, path = create_job(root, method, "positive", "peak_height")
    out = Path(job.run_dir)
    (out / ".detached-state.json").write_text('{"pid":123,"befores":{}}')
    (out / "intermediate.pai2").write_bytes(b"intermediate")
    monkeypatch.setattr("metabolomix.console.runner.is_process_running", lambda pid: False)
    console_status(str(path))
    assert load_job(path).status != "completed"


def test_legacy_detached_state_reports_execution_unresolved(tmp_path, monkeypatch):
    """終了コード不明は成功でも失敗でもない。機械可読に「未解決」と言う。"""
    from metabolomix.tools.console_tools import console_status
    job_path = _planned(tmp_path, monkeypatch)
    from metabolomix.console.detached import write_detached_state
    write_detached_state(_run_dir(job_path), 123, {})
    monkeypatch.setattr("metabolomix.console.runner.is_process_running", lambda pid: False)

    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["error"]["code"] == "EXECUTION_UNRESOLVED"
    assert parsed["error"]["details"]["pid"] == 123


def test_legacy_detached_state_still_alive_is_reported_as_running(tmp_path, monkeypatch):
    """まだ生きている旧実行は、未解決ではなく実行中として読める。"""
    from metabolomix.tools.console_tools import console_status
    job_path = _planned(tmp_path, monkeypatch)
    from metabolomix.console.detached import write_detached_state
    write_detached_state(_run_dir(job_path), 123, {})
    monkeypatch.setattr("metabolomix.console.runner.is_process_running", lambda pid: True)

    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["detached"] == {"pid": 123, "alive": True, "legacy": True}


def test_is_process_running_true_for_self():
    """旧 state の安全な検出のため互換 wrapper は残す。"""
    from metabolomix.console.runner import is_process_running
    assert is_process_running(os.getpid()) is True


def test_is_process_running_false_for_absent_pid():
    from metabolomix.console.runner import is_process_running
    # 実在しない可能性が極めて高い pid。誤検出しても False を返すのが安全側。
    assert is_process_running(999_999_999) is False


# ---------- ワーカーの起動 ----------

def test_worker_command_runs_the_module_from_the_repo_root(tmp_path, monkeypatch):
    """`-m metabolomix.console.worker` は cwd がリポジトリルートでなければ解決しない。"""
    from metabolomix.console import worker
    seen: dict = {}

    def _fake_launch(command, *, cwd, log_path):
        seen.update(command=command, cwd=cwd, log_path=log_path)
        return {"pid": 4242, "identity": {"pid": 4242, "creation_time": 7},
                "breakaway": "not_in_job", "log_path": str(log_path)}

    monkeypatch.setattr("metabolomix.console.worker.launch_detached", _fake_launch)
    job_path = _planned(tmp_path, monkeypatch)

    info = worker.launch_console_worker(job_path)

    assert info["pid"] == 4242
    assert seen["command"] == [sys.executable, "-m", "metabolomix.console.worker",
                               "--job", str(job_path)]
    assert (seen["cwd"] / "metabolomix" / "console" / "worker.py").is_file()
    assert Path(seen["log_path"]).name == "worker.log"


def test_console_run_detach_returns_the_worker_pid(tmp_path, monkeypatch):
    """受理応答の pid はワーカー。Console の pid は起動後の証跡に載る。"""
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setattr(
        "metabolomix.console.worker.launch_detached",
        lambda command, *, cwd, log_path: {
            "pid": 31337, "identity": {"pid": 31337, "creation_time": 11},
            "breakaway": "not_in_job", "log_path": str(log_path)})
    job_path = _planned(tmp_path, monkeypatch)

    parsed = _json.loads(console_run(str(job_path), detach=True))

    assert parsed["status"] == "running"
    assert parsed["pid"] == 31337
    assert parsed["pid_of"] == "worker"
    assert load_job(job_path).status == "running"


def test_console_run_detach_fixes_the_inputs_before_launching(tmp_path, monkeypatch):
    """何を入力として実行するかは、起動より前に固定しなければ検証できない。"""
    from metabolomix.console.execution import read_supervision_state
    from metabolomix.tools.console_tools import console_run
    seen: dict = {}

    def _fake_launch(command, *, cwd, log_path):
        seen["state"] = read_supervision_state(Path(log_path).parent)
        return {"pid": 4242, "identity": {"pid": 4242, "creation_time": 7},
                "breakaway": "not_in_job", "log_path": str(log_path)}

    monkeypatch.setattr("metabolomix.console.worker.launch_detached", _fake_launch)
    job_path = _planned(tmp_path, monkeypatch, raw_count=3)

    console_run(str(job_path), detach=True)

    inputs = seen["state"]["inputs"]
    assert len(inputs["raw_inventory"]) == 3
    assert all(p.endswith(".abf") for p in inputs["raw_inventory"])
    assert len(inputs["method_sha256"]) == 64
    assert inputs["exe_path"].endswith("fake-msdialcui.exe")


def test_console_run_detach_writes_no_legacy_state(tmp_path, monkeypatch):
    """旧サイドカーを新規に書かない（読み手が 2 つの真実を持たされる）。"""
    from metabolomix.console.detached import read_detached_state
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setattr(
        "metabolomix.console.worker.launch_detached",
        lambda command, *, cwd, log_path: {
            "pid": 31337, "identity": {"pid": 31337, "creation_time": 11},
            "breakaway": "not_in_job", "log_path": str(log_path)})
    job_path = _planned(tmp_path, monkeypatch)

    console_run(str(job_path), detach=True)

    assert read_detached_state(_run_dir(job_path)) is None


def test_console_run_detach_reports_the_worker_even_if_the_owner_write_fails(
        tmp_path, monkeypatch):
    """起動後に所有記録を書けなくても、走り出したワーカーを黙って捨てない。"""
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setattr(
        "metabolomix.console.worker.launch_detached",
        lambda command, *, cwd, log_path: {
            "pid": 31337, "identity": {"pid": 31337, "creation_time": 11},
            "breakaway": "not_in_job", "log_path": str(log_path)})

    def _boom(run_dir, owner):
        raise OSError("worker.json を書けない")

    monkeypatch.setattr("metabolomix.console.worker.write_owner", _boom)
    job_path = _planned(tmp_path, monkeypatch)

    parsed = _json.loads(console_run(str(job_path), detach=True))

    assert parsed["status"] == "running"
    assert parsed["pid"] == 31337
    assert any("worker.json" in w for w in parsed["warnings"])


def test_console_run_detach_marks_failed_when_the_worker_cannot_start(tmp_path, monkeypatch):
    """起動できなかったなら running のまま放置しない。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.tools.console_tools import console_run

    def _boom(command, *, cwd, log_path):
        raise DomainError("DETACH_UNSUPPORTED", "親 Job から抜けられない")

    monkeypatch.setattr("metabolomix.console.worker.launch_detached", _boom)
    job_path = _planned(tmp_path, monkeypatch)

    parsed = _json.loads(console_run(str(job_path), detach=True))

    assert parsed["error"]["code"] == "WORKER_LAUNCH_FAILED"
    assert load_job(job_path).status == "failed"


# ---------- job ID ----------

def test_job_ids_are_unique_within_the_same_second(tmp_path):
    """同じ秒に 2 件計画しても run_dir を共有させない（UUID を含める）。"""
    source = tmp_path / "source"
    source.mkdir()
    method = source / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")

    first, first_path = create_job(source, method, "negative", "peak_height")
    second, second_path = create_job(source, method, "negative", "peak_height")

    assert first.job_id != second.job_id
    assert first.run_dir != second.run_dir
    assert first_path.is_file() and second_path.is_file()


# ---------- 実ワーカー（Windows のみ） ----------

@_WINDOWS_ONLY
def test_worker_finalizes_the_job_without_any_console_status_call(tmp_path, monkeypatch):
    """console_status を一度も呼ばなくても、ワーカーがジョブを確定させる。"""
    from metabolomix.console import worker
    job_path = _planned(tmp_path, monkeypatch, raw_count=2)
    _prime_supervision(job_path)
    counter = tmp_path / "launches.tsv"

    worker.launch_console_worker(
        job_path, console_args=_fake_console_args(job_path, counter, "success"))

    status = _wait_for_status(job_path, {"completed", "partial", "failed"})
    assert status == "completed"

    receipt = _json.loads(
        (_run_dir(job_path) / "execution-result.json").read_text(encoding="utf-8"))
    assert (receipt["termination"], receipt["exit_code"]) == ("exited", 0)
    assert counter.read_text(encoding="utf-8").count("\n") == 1


@_WINDOWS_ONLY
def test_worker_records_a_nonzero_exit_instead_of_completing(tmp_path, monkeypatch):
    """非ゼロ終了は成果物が残っていても completed にしない。"""
    from metabolomix.console import worker
    job_path = _planned(tmp_path, monkeypatch, raw_count=2)
    _prime_supervision(job_path)
    counter = tmp_path / "launches.tsv"

    worker.launch_console_worker(
        job_path, console_args=_fake_console_args(job_path, counter, "nonzero"))

    status = _wait_for_status(job_path, {"completed", "partial", "failed"})
    assert status == "partial"
    receipt = _json.loads(
        (_run_dir(job_path) / "execution-result.json").read_text(encoding="utf-8"))
    assert receipt["exit_code"] == 1


@_WINDOWS_ONLY
def test_a_second_execution_is_refused_while_the_job_is_locked(tmp_path, monkeypatch):
    """同じジョブを 2 つのプロセスが同時に走らせない（OS ロックで排他する）。"""
    from metabolomix.console.worker import lock_path
    from metabolomix.core.process_control import file_lock
    from metabolomix.tools.console_tools import console_run
    job_path = _planned(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.worker.LOCK_TIMEOUT_S", 0.2)

    with file_lock(lock_path(_run_dir(job_path))):
        parsed = _json.loads(console_run(str(job_path)))

    assert parsed["error"]["code"] == "JOB_BUSY"


# ---------- 旧スキーマのジョブ ----------

def test_console_status_reads_a_legacy_v1_job_without_rewriting_it(tmp_path):
    """analysis-job.v1 のジョブも読める。状態照会は読むだけで書き換えない。

    v1 には run_dir 側の証跡も detach 状態も無い。「終わったことになっている」
    古い記録をそのまま表示し、勝手に再解釈しない。
    """
    from metabolomix.tools.console_tools import console_status
    run_dir = tmp_path / "runs" / "old"
    run_dir.mkdir(parents=True)
    job_path = run_dir / "analysis-job.json"
    original = {
        "schema": "analysis-job.v1", "job_id": "old", "status": "completed",
        "created_at": "t", "updated_at": "t",
        "source": {"dataset_root": str(tmp_path), "input_count": 1},
        "software": {"name": "MS-DIAL", "version": "", "execution_mode": "console",
                     "method_file": "m.txt"},
        "project": {"omics": "lipidomics", "polarity": "negative",
                    "measure": "peak_height"},
        "run_dir": str(run_dir),
        "primary_mztab_files": [{"path": "A.mzTab", "polarity": "negative",
                                 "measure": "peak_height", "sha256": "x",
                                 "validation": {}}],
        "artifacts": [{"path": "S1.pai2", "role": "sample_peaks",
                       "format": "pai2", "sha256": "y"}],
        "warnings": [], "error": None,
    }
    job_path.write_text(_json.dumps(original), encoding="utf-8")

    parsed = _json.loads(console_status(str(job_path)))

    assert parsed["status"] == "completed"
    assert parsed["mztab_files"][0]["root"] == "run_dir"
    assert "execution_receipt" not in parsed  # 証跡が無いものを作り出さない
    assert _json.loads(job_path.read_text(encoding="utf-8")) == original


# ---------- 先に終わったワーカーの結果を親が潰さない（レビュー指摘P2） ----------

def test_console_run_detach_does_not_overwrite_a_worker_that_already_finished(
        tmp_path, monkeypatch):
    """起動直後に走り切ったワーカーの completed を、親が running で上書きしない。

    `console_status` は読取専用なので、ここで潰された状態は誰も直せない
    ——終わっている実行が永久に running のまま残る。
    """
    from metabolomix.console.job_manager import update_status
    from metabolomix.tools.console_tools import console_run

    job_path = _planned(tmp_path, monkeypatch)

    def _fast_worker(command, *, cwd, log_path):
        # 起動したワーカーが、親が戻るより先に走り切った状況。
        update_status(job_path, "running")
        update_status(job_path, "completed")
        return {"pid": 31337, "identity": {"pid": 31337, "creation_time": 11},
                "breakaway": "not_in_job", "log_path": str(log_path)}

    monkeypatch.setattr("metabolomix.console.worker.launch_detached", _fast_worker)

    parsed = _json.loads(console_run(str(job_path), detach=True))

    assert load_job(job_path).status == "completed"
    assert parsed["status"] == "completed"
    assert parsed["pid"] == 31337          # 起動した事実は変わらず返す


def test_console_run_detach_does_not_overwrite_the_owner_written_by_the_worker(
        tmp_path, monkeypatch):
    """ワーカーが自分の identity を書いていたら、親の pid 記録で上書きしない。"""
    from metabolomix.console import worker as console_worker
    from metabolomix.tools.console_tools import console_run

    job_path = _planned(tmp_path, monkeypatch)
    run_dir = _run_dir(job_path)

    def _worker_that_claims_ownership(command, *, cwd, log_path):
        console_worker.write_owner(run_dir, {
            "kind": console_worker.OWNER_KIND_CONSOLE, "pid": 999,
            "identity": {"pid": 999, "creation_time": 5},
            "job_path": str(job_path), "status": "running"})
        return {"pid": 31337, "identity": {"pid": 31337, "creation_time": 11},
                "breakaway": "not_in_job", "log_path": str(log_path)}

    monkeypatch.setattr("metabolomix.console.worker.launch_detached",
                        _worker_that_claims_ownership)

    console_run(str(job_path), detach=True)

    assert console_worker.read_owner(run_dir)["pid"] == 999


# ---------- 二重実行と所有の解除（レビュー指摘P2） ----------

def test_run_job_refuses_a_job_that_another_process_already_finished(tmp_path, monkeypatch):
    """ロック待ちの間に先行プロセスが走り切っていたら、MS-DIAL を二度起動しない。"""
    from metabolomix.console import worker as console_worker
    from metabolomix.console.job_manager import update_status
    from metabolomix.core.atomic_io import DomainError

    job_path = _planned(tmp_path, monkeypatch)
    _prime_supervision(job_path)
    update_status(job_path, "completed")

    def _never(*args, **kwargs):
        raise AssertionError("完了済みジョブで Console を起動してはいけない")

    monkeypatch.setattr("metabolomix.console.worker.supervise", _never)

    with pytest.raises(DomainError) as excinfo:
        console_worker.run_job(job_path)
    assert excinfo.value.code == "JOB_ALREADY_FINISHED"
    assert load_job(job_path).status == "completed"   # 勝者の状態を壊さない


def test_console_run_reports_an_already_finished_job_without_touching_its_state(
        tmp_path, monkeypatch):
    """同時に届いた2本目の console_run は、勝者の終端状態を failed で塗り替えない。"""
    from metabolomix.console.job_manager import update_status
    from metabolomix.tools.console_tools import console_run

    job_path = _planned(tmp_path, monkeypatch)
    _prime_supervision(job_path)

    def _finish_first(job_path_arg, *, command=None):
        # ロックを取るまでの間に先行プロセスが走り切っていた、という筋書き。
        from metabolomix.core.atomic_io import DomainError
        update_status(job_path, "completed")
        raise DomainError("JOB_ALREADY_FINISHED", "既に実行を終えています",
                          {"status": "completed"})

    monkeypatch.setattr("metabolomix.console.worker.run_job", _finish_first)

    parsed = _json.loads(console_run(str(job_path)))

    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"
    assert load_job(job_path).status == "completed"


def test_owner_is_released_when_the_run_finishes_in_this_process(tmp_path, monkeypatch):
    """同期実行の所有者は MCP サーバ自身。終わったあとも生きているのは当然で、
    それを「まだ実行中」と読むと console_cleanup が永久に JOB_BUSY になる。
    """
    import os

    from metabolomix.console import worker as console_worker
    from metabolomix.core.process_control import process_identity

    job_path = _planned(tmp_path, monkeypatch)
    run_dir = _run_dir(job_path)
    identity = process_identity(os.getpid())
    console_worker.write_owner(run_dir, {
        "kind": console_worker.OWNER_KIND_CONSOLE, "pid": identity["pid"],
        "identity": identity, "job_path": str(job_path), "status": "running"})
    assert console_worker.owner_is_active(run_dir) is True

    console_worker.clear_owner(run_dir)

    assert console_worker.owner_is_active(run_dir) is False
    assert console_worker.owner_summary(run_dir)["status"] == "finished"
