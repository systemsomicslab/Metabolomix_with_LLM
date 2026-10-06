"""切断・喪失・取消・並行実行を実プロセスで検証する（spec §13 A04〜A06・D06・D09）。

このファイルの全ケースは**本物のプロセス**を起こす。Job Objectによる一族停止・
pid再利用・切り離し起動はWindows固有のAPIで実装されているため、モックで
置き換えると「停止したつもり」「生きているつもり」が緑のまま通ってしまう。

R22（Task16からの持ち越し）の実証もここで行う。「所有者が分からないrunの取消が、
無関係なプロセスを殺さないこと」は、終了APIを呼ぶ場所（＝実際にプロセスを
止める経路）でしか確かめられない。`recovery.request_cancel`が終了APIを一切
呼ばないという実装上の事実だけでは、**取消の全経路が無害である**ことの証明に
ならない——巻き添えの当事者になりうる無関係プロセスを実際に立てて、取消を
一通り通した後もそれが生きていることを見る。
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from lipidmix.core.process_control import process_identity, resolve_python_launch, same_process
from lipidmix.pipeline import engine, recovery
from tests.pipeline_fixtures import DEFAULT_COMPARISON, PipelineHarness, read_contract_tsv

pytestmark = pytest.mark.skipif(
    os.name != "nt",
    reason="Job Object・pid再利用の見分け・切り離し起動はWindows専用の実装")

#: 巻き添え検査用の無関係プロセスの生存時間（秒）。テストが必ず自分で始末する。
_BYSTANDER_LIFETIME = 300


@pytest.fixture
def pipeline_harness(tmp_path, monkeypatch):
    harness = PipelineHarness(tmp_path, monkeypatch)
    try:
        yield harness
        # 子プロセスを漏らしていないこと（worker・偽Console・その孫まで）。
        assert harness.survivors() == [], f"プロセスが残っています: {harness.survivors()}"
    finally:
        harness.close()


class Bystander:
    """このpipelineとは何の関係も無い、生きている実プロセス。

    取消・喪失判定・resumeを通したあとで**まだ生きていること**を確かめるための
    巻き添え検査用。identityを持つので、pidが再利用されても取り違えない。
    """

    def __init__(self):
        # venvのリダイレクタを挟まない（挟むとidentityがリダイレクタのものになる）。
        command, env = resolve_python_launch(
            [sys.executable, "-c", f"import time; time.sleep({_BYSTANDER_LIFETIME})"])
        self.proc = subprocess.Popen(
            command, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.identity = process_identity(self.proc.pid)
        assert self.identity is not None, "巻き添え検査用プロセスを起こせませんでした"

    @property
    def pid(self) -> int:
        return int(self.proc.pid)

    def alive(self) -> bool:
        return same_process(self.identity)

    def close(self) -> None:
        if self.proc.poll() is None:
            self.proc.kill()
        self.proc.wait(timeout=30)


@pytest.fixture
def bystander():
    process = Bystander()
    try:
        yield process
    finally:
        process.close()


def _kill_process_tree(proc) -> None:
    """workerを「落ちた」状態にする（猶予なしの強制終了）。"""
    subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                   stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, timeout=60)
    proc.wait(timeout=60)


# ---------- A05: MCP相当の親が消えても、statusを一度も呼ばずに完走する ----------

def test_run_completes_after_the_launcher_process_exits_without_any_status_call(
        pipeline_harness):
    """spec A05/D06: 起動役が消えたあとも、workerが終了・成果物確定・下流進行を行う。

    `pipeline_status`/`read_status`はこのテストの中で一度も呼ばない
    （harnessの`wait`も状態ツールを通さず、workerの生死だけを見て最後に1回読む）。
    """
    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"},
        launch=False)  # plan止め。起動は切り離し起動役に任せる。
    info = pipeline_harness.launch_detached(run)

    # 起動役はもう居ない（`launch_detached`は起動役プロセスの終了を待ってから返る）。
    assert info["pid"] > 0 and info["identity"]["pid"] == info["pid"]
    assert info["breakaway"] != "unsupported"

    record = pipeline_harness.wait(run, expected="completed")

    assert pipeline_harness.console_start_count(run) == 1
    names = {ref["output_name"] for ref in record["results"]}
    assert {"preprocess", "pca", "pca_figure", "differential:treated_vs_control",
            "volcano:treated_vs_control", "tsv:treated_vs_control",
            "quality_report"} <= names
    assert Path(pipeline_harness.output_path(record, "tsv:treated_vs_control")).is_file()


# ---------- A04: timeoutで一族が止まる ----------

def test_timeout_stops_the_whole_console_family_and_leaves_a_receipt(pipeline_harness):
    """spec A04: 非同期でtimeout到達、子プロセスあり → 一族停止＋timeout証跡。"""
    run = pipeline_harness.start_scenario("hang", timeout_s=3)
    record = pipeline_harness.wait(run, expected="failed")

    rows = pipeline_harness.console_launch_rows(run)
    assert len(rows) == 1
    console_pid, grandchild_pid = int(rows[0][1]), int(rows[0][2])
    assert grandchild_pid > 0, "孫プロセスが起動していません（一族停止の検証にならない）"
    assert process_identity(console_pid) is None, "偽Consoleが残っています"
    assert process_identity(grandchild_pid) is None, "孫プロセスが残っています"

    receipt = pipeline_harness.console_receipt(run)
    assert receipt["termination"] == "timeout"
    assert receipt["timeout_s"] == 3
    assert receipt["exit_code"] != 0
    assert record["upstream"]["verification"]["termination"] == "timeout"


# ---------- D09/R22: 取消 ----------

def test_cancel_during_upstream_stops_only_the_owned_console_tree(
        pipeline_harness, bystander):
    """R22: 取消はpipelineが所有するConsole一族だけを止め、無関係なPIDに触らない。"""
    run = pipeline_harness.start_scenario("hang", timeout_s=600)
    pipeline_harness.wait_for(
        lambda: len(pipeline_harness.console_launch_rows(run)) == 1,
        what="偽Consoleと孫プロセスの起動")
    rows = pipeline_harness.console_launch_rows(run)
    console_pid, grandchild_pid = int(rows[0][1]), int(rows[0][2])
    assert process_identity(grandchild_pid) is not None
    assert bystander.alive()

    receipt = recovery.request_cancel(pipeline_harness.pipeline_root(run))
    assert receipt["accepted"] is True

    # 工程の**最中**に取り消しても `cancelled` に確定する。以前はここだけ
    # `failed` になっていた——`cancel_requested` の判定が工程の境界でしか
    # 走らず、Console を殺された handler の failed が格下げ規則に落ちていたため。
    # 境界での取消（`test_cancel_at_a_stage_boundary_keeps_earlier_results`）と
    # 語彙が食い違い、`pipeline_cancel` が契約する「`cancelled` への確定」を
    # 待つクライアントが永久に確認できなかった。
    record = pipeline_harness.wait(run, expected="cancelled")

    # 所有する一族は止まった。
    assert process_identity(console_pid) is None
    assert process_identity(grandchild_pid) is None
    assert pipeline_harness.console_receipt(run)["termination"] == "cancelled"
    assert record["stages"]["upstream"]["error"]["code"] == "MSDIAL_EXECUTION_FAILED"
    # 無関係なプロセスには触れていない。
    assert bystander.alive(), "取消が無関係なプロセスを巻き込みました"
    # 再実行を捏造しない（Consoleは1回しか起動していない）。
    assert pipeline_harness.console_start_count(run) == 1


def test_cancel_at_a_stage_boundary_keeps_earlier_results(pipeline_harness):
    """spec D09: 下流は工程境界で止まり、そこまでの有効な成果物は残る。"""
    run = pipeline_harness.start(target="exploratory",
                                 sleep_stage="pca", sleep_seconds=5)
    pipeline_harness.wait_for_stage(run, "pca")

    recovery.request_cancel(pipeline_harness.pipeline_root(run))
    record = pipeline_harness.wait(run, expected="cancelled")

    assert record["status"] == "cancelled"
    names = {ref["output_name"] for ref in record["results"]}
    assert {"preprocess", "pca", "pca_figure"} <= names
    assert record["stages"]["report"]["status"] == "pending"
    # 取消要求は受理された事実として残る（受理と停止完了は別状態）。
    assert engine.cancel_requested(pipeline_harness.pipeline_root(run)) is True


def test_cancelled_run_can_be_resumed(pipeline_harness):
    """spec D09: 取消後のresumeは、過去attemptを保持したまま続きから進められること。

    resumeで渡すのはシートだけで、比較は足さない（`comparisons=[]`）
    ——このrunは`target="exploratory"`であり、探索目標へ比較を足す要求は
    矛盾として受付で拒否されるようになった（最終レビュー指摘5）。
    このテストが見るのは取消後の再開そのものなので、比較の有無は本題では
    ない。assertは一切緩めていない。
    """
    run = pipeline_harness.start(target="exploratory",
                                 sleep_stage="pca", sleep_seconds=5)
    pipeline_harness.wait_for_stage(run, "pca")
    recovery.request_cancel(pipeline_harness.pipeline_root(run))
    cancelled = pipeline_harness.wait(run, expected="cancelled")
    before_attempt = cancelled["stages"]["pca"]["attempt"]

    pipeline_harness.resume_with_groups(run, comparisons=[], sleep_seconds=0.0)
    resumed = pipeline_harness.wait(run, expected="completed")

    assert resumed["stages"]["pca"]["attempt"] > before_attempt
    assert pipeline_harness.console_start_count(run) == 1  # 上流は再実行しない


@pytest.mark.parametrize("identity_factory, expected_health, hint_code", [
    (lambda pid: None, "unknown", "SUPERVISION_UNKNOWN"),
    (lambda pid: {"pid": pid, "creation_time": 1}, "worker_missing", "PIPELINE_INTERRUPTED"),
])
def test_cancel_with_unknown_owner_never_touches_another_process(
        pipeline_harness, bystander, identity_factory, expected_health, hint_code):
    """R22: 所有者identityが不明／pid再利用の状態で取消しても、別PIDを止めない。

    2つの「所有者不明」を作る。

    * identityがそもそも無い（起動直後にworkerを失った等）。
    * 記録されたpidが**無関係な生きているプロセス**に再利用されている
      （creation_timeが違うので`same_process`はFalseになるべき）。

    どちらでも、取消の受理・状態読取・再開準備を通した後で、その無関係な
    プロセスが生きていることを確かめる。
    """
    from lipidmix.pipeline import store

    run = pipeline_harness.start(target="exploratory", launch=False)
    pipeline_root = pipeline_harness.pipeline_root(run)

    record = store.load_run(pipeline_root)
    record["status"] = "running"
    record["worker"] = {"identity": identity_factory(bystander.pid),
                        "started_at": "2026-09-05T00:00:00+00:00"}
    store.save_run(pipeline_root, record, expected_revision=record["state_revision"])

    cancel_receipt = recovery.request_cancel(pipeline_root)
    status = recovery.read_status(pipeline_root)

    assert cancel_receipt["accepted"] is True
    assert engine.cancel_requested(pipeline_root) is True
    # 所有者を確定できないので、完了も停止完了も捏造しない。
    assert status["status"] == "running"
    assert status["observed_health"] == expected_health
    assert status["recovery_hint"]["code"] == hint_code
    # 無関係なプロセスは生きたまま（pid再利用を根拠に停止していない）。
    assert bystander.alive(), "所有者不明の取消が無関係なPIDを止めました"
    assert process_identity(bystander.pid) is not None


# ---------- A06: workerを失う ----------

def _lose_the_worker(pipeline_harness):
    """pca工程の途中でworkerを強制終了し、(pipeline_root, run)を返す。"""
    run = pipeline_harness.start(target="exploratory",
                                 sleep_stage="pca", sleep_seconds=60)
    pipeline_harness.wait_for_stage(run, "pca")
    worker = pipeline_harness.worker_processes(run)[0]
    worker_identity = process_identity(worker.pid)
    assert worker_identity is not None

    _kill_process_tree(worker)
    assert not same_process(worker_identity)
    return pipeline_harness.pipeline_root(run), run


def test_lost_worker_is_reported_and_never_fabricated(pipeline_harness):
    """spec A06: workerが落ちても、完了も再起動も捏造せず「中断」として報告する。"""
    pipeline_root, run = _lose_the_worker(pipeline_harness)

    status = recovery.read_status(pipeline_root)

    assert status["status"] == "running"          # 永続状態は書き換えない
    assert status["observed_health"] == "worker_missing"
    assert status["recovery_hint"]["code"] == "PIPELINE_INTERRUPTED"
    record = pipeline_harness.record(run)
    assert record["stages"]["pca"]["status"] == "running"   # 成功を捏造しない
    assert record["status"] == "running"
    assert record["stages"]["report"]["status"] == "pending"
    # 勝手に再起動していない（Consoleは最初の1回だけ）。
    assert pipeline_harness.console_start_count(run) == 1
    assert "quality_report" not in {ref["output_name"] for ref in record["results"]}


def test_lost_run_can_be_resumed_as_read_status_advertises(pipeline_harness):
    """spec A06: 中断したrunは、案内どおりのresumeで最後まで進められること。"""
    from lipidmix.pipeline import service

    pipeline_root, run = _lose_the_worker(pipeline_harness)
    assert recovery.read_status(pipeline_root)["recovery_hint"]["code"] == "PIPELINE_INTERRUPTED"

    pipeline_harness.set_launch_options(run, sleep_stage=None)
    service.resume_pipeline(pipeline_root)
    resumed = pipeline_harness.wait(run, expected="completed")

    assert pipeline_harness.console_start_count(run) == 1  # 上流は再実行しない
    assert {"preprocess", "pca", "pca_figure", "quality_report"} <= \
        {ref["output_name"] for ref in resumed["results"]}


def test_resumed_export_alone_still_names_a_recorded_differential_result(pipeline_harness):
    """spec §6.1/§9.1: `export`だけが動く再開passでも、出所の連結が切れないこと。

    `export`の途中でworkerを失うと、`differential:<cid>`はsucceededのまま残る。
    再開したworkerはそれをskipするので、`export`は差次的結果を持たないまま
    走り出す——ここで結果を計算し直すだけで永続化を忘れると、TSVの
    `# result_id`とvolcanoの`parent_ids`が`record["results"]`に存在しない
    結果を指す（Task19が見つけた欠陥の、再開経路側の顔）。
    """
    from lipidmix.pipeline import service

    cid = DEFAULT_COMPARISON["comparison_id"]
    pipeline_harness.write_manifest()
    run = pipeline_harness.start(
        target="differential", comparisons=[DEFAULT_COMPARISON],
        extra_request={"sample_manifest": "sample-manifest.tsv"},
        sleep_stage=f"export:{cid}", sleep_seconds=60)
    pipeline_harness.wait_for_stage(run, f"export:{cid}")
    worker = pipeline_harness.worker_processes(run)[0]
    _kill_process_tree(worker)
    interrupted = pipeline_harness.record(run)
    assert interrupted["stages"][f"differential:{cid}"]["status"] == "succeeded"

    pipeline_harness.set_launch_options(run, sleep_stage=None)
    service.resume_pipeline(pipeline_harness.pipeline_root(run))
    resumed = pipeline_harness.wait(run, expected="completed")

    # differentialは再実行されず（skip）、exportだけがこのpassで動いた。
    assert pipeline_harness.console_start_count(run) == 1
    known_ids = {ref["result_id"] for ref in resumed["results"]}
    dangling = {ref["output_name"]: [pid for pid in ref.get("parent_ids") or []
                                     if pid not in known_ids]
                for ref in resumed["results"]}
    assert {name: ids for name, ids in dangling.items() if ids} == {}
    meta, _rows, _fields = read_contract_tsv(
        pipeline_harness.output_path(resumed, f"tsv:{cid}"))
    assert meta["result_id"] in known_ids


def test_second_worker_on_the_same_run_refuses_instead_of_running_twice(pipeline_harness):
    """owner lock: 同じrunに2つ目のworkerが来ても、二重に進めずに拒否する。"""
    run = pipeline_harness.start(target="exploratory",
                                 sleep_stage="pca", sleep_seconds=5)
    pipeline_harness.wait_for_stage(run, "pca")

    pipeline_harness.relaunch(run)
    record = pipeline_harness.wait(run, expected="completed")

    first, second = pipeline_harness.worker_processes(run)
    assert first.returncode == 0
    assert second.returncode == 1, "2つ目のworkerが拒否されずに走りました"
    logs = "\n".join(path.read_text(encoding="utf-8", errors="replace")
                     for path in sorted(pipeline_harness.logs_dir.glob("*.log")))
    assert "PIPELINE_ALREADY_RUNNING" in logs
    assert pipeline_harness.console_start_count(run) == 1
    assert record["status"] == "completed"


# ---------- D06: 2つの異なるpipelineの並行実行 ----------

def test_two_concurrent_runs_do_not_mix_samples_matrices_or_figures(pipeline_harness):
    """spec D06: 群割当を逆にした2つのrunを同時に走らせても混ざらない。"""
    reversed_groups = ("treated",) * 4 + ("control",) * 4
    pipeline_harness.write_manifest(source_name="alpha")
    pipeline_harness.write_manifest(source_name="beta", groups=reversed_groups)

    # pca工程で待たせることで、「両方が同時に走っている瞬間」を作って観測する
    # （先に終わった片方を後から見る、という順次実行では並行の検証にならない）。
    def _start(source_name):
        return pipeline_harness.start(
            target="differential", comparisons=[DEFAULT_COMPARISON],
            source_name=source_name,
            extra_request={"sample_manifest": "sample-manifest.tsv"},
            sleep_stage="pca", sleep_seconds=4)

    alpha, beta = _start("alpha"), _start("beta")
    alpha_worker = pipeline_harness.worker_processes(alpha)[0]
    beta_worker = pipeline_harness.worker_processes(beta)[0]
    assert alpha_worker.pid != beta_worker.pid

    pipeline_harness.wait_for_stage(alpha, "pca")
    pipeline_harness.wait_for_stage(beta, "pca")
    assert alpha_worker.poll() is None and beta_worker.poll() is None, \
        "2つのrunが同時に走っていません（並行実行の検証になっていない）"

    alpha_record = pipeline_harness.wait(alpha, expected="completed")
    beta_record = pipeline_harness.wait(beta, expected="completed")

    assert (alpha_record["identity"]["pipeline_id"]
            != beta_record["identity"]["pipeline_id"])
    assert pipeline_harness.console_start_count(alpha) == 1
    assert pipeline_harness.console_start_count(beta) == 1

    alpha_meta, alpha_rows, _ = read_contract_tsv(
        pipeline_harness.output_path(alpha_record, "tsv:treated_vs_control"))
    beta_meta, beta_rows, _ = read_contract_tsv(
        pipeline_harness.output_path(beta_record, "tsv:treated_vs_control"))

    # 検体・行列・比較が混ざっていない: 群割当が逆なので log2FC の符号も逆になる。
    assert all(float(row["log2fc"]) > 0 for row in alpha_rows)
    assert all(float(row["log2fc"]) < 0 for row in beta_rows)
    assert alpha_meta["result_id"] != beta_meta["result_id"]
    # 成果物refはすべて自分のpipeline_root配下の実ファイルを指している。
    for record in (alpha_record, beta_record):
        root = Path(record["identity"]["pipeline_root"])
        for ref in record["results"]:
            assert (root / ref["relative_path"]).is_file()
    # 各runのmzTabは自分のinput配下だけを見ている。
    assert str(Path(alpha_record["identity"]["pipeline_root"])) in alpha_meta["source_mztab"]
    assert str(Path(beta_record["identity"]["pipeline_root"])) in beta_meta["source_mztab"]


# ---------- 後始末 ----------

def test_stopping_the_worker_takes_its_console_family_and_nothing_else(
        pipeline_harness, bystander):
    """漏れ検査そのものの検証: 検出器が空でないことを示してから、空になるのを見る。

    「無関係なpython.exeの総数」ではなく**harnessが起こしたプロセス**を数える
    （他のプロセスが同時に増減する環境で意味のある判定にするため）。まず
    worker・偽Console・その孫の3種が生きているのを検出し、次にworkerだけを
    止めて、Job ObjectのKILL_ON_JOB_CLOSEで一族が道連れになり、巻き添え検査用
    プロセスだけが残ることを見る。
    """
    run = pipeline_harness.start_scenario("hang", timeout_s=600)
    pipeline_harness.wait_for(
        lambda: len(pipeline_harness.console_launch_rows(run)) == 1,
        what="偽Consoleと孫プロセスの起動")

    kinds = {entry["kind"] for entry in pipeline_harness.survivors()}
    assert {"worker", "fake_console", "fake_console_grandchild"} <= kinds, kinds

    pipeline_harness.stop_workers()  # 所有workerだけを止める

    pipeline_harness.wait_for(lambda: pipeline_harness.survivors() == [],
                              what="worker一族が全部消えること")
    assert bystander.alive(), "workerの停止が無関係なプロセスを巻き込みました"
    os.kill(bystander.pid, signal.SIGTERM)
    pipeline_harness.wait_for(lambda: not bystander.alive(),
                              what="巻き添え検査用プロセスの終了")
