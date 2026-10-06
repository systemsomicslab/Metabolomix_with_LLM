"""再開・再構築・取消・中断の読取表示のテスト（Task16）。

engineの`stage_inputs_unchanged`/`make_context`の拡張はここで、注入handler
（test_pipeline_engine.pyと同じ流儀）を使って検証する——本番handler一式
（`build_handlers`）はTask18が実装する。fixtureは合成dictビルダーとして
本ファイル内に閉じる（controller裁定: 既存`tests.pipeline_fixtures`の
契約は壊さない範囲でしか触らない）。
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError, atomic_write_json
from metabolomix.core.process_control import process_identity
from metabolomix.pipeline.engine import cancel_request_path, run_engine
from metabolomix.pipeline.recovery import prepare_resume, read_status, request_cancel
from metabolomix.pipeline.request import resolve_request
from metabolomix.pipeline.store import create_run, load_run, save_run
import metabolomix.pipeline.recovery as recovery_module
import metabolomix.pipeline.store as store_module

from tests.pipeline_fixtures import make_source

REPO_ROOT = Path(__file__).resolve().parents[1]


# ---------- 合成fixture ----------

def _minimal_inputs(root: Path, *, fingerprint: str = "f" * 64) -> dict:
    return {"source_root": str(root), "fingerprint": fingerprint, "raw_inventory": []}


def _persist(pipeline_root, name: str):
    """`output_name=name`の実ファイル付きref（Task18: `evaluate_target`は
    hash照合込みで`output_name`付きrefだけを「達成」と数えるため、文字列の
    ダミーrefでは`finish_success`が常にfailed/partialへ落ちる）。"""
    from metabolomix.pipeline import report as report_mod
    return report_mod.persist_result(pipeline_root, {
        "output_name": name, "kind": "synthetic",
        "result_id": f"{name.replace(':', '_')}-result",
        "data": {"synthetic": True, "name": name},
    })


def _recorder(calls: list[str], *, status: str = "succeeded",
              output_names: list[str] | None = None, error: dict | None = None):
    def handler(context: dict) -> dict:
        calls.append(context["stage_id"])
        refs = [_persist(context["pipeline_root"], name) for name in (output_names or [])]
        return {"status": status, "result_refs": refs, "warnings": [], "error": error}
    return handler


def _upstream_recorder(calls: list[str]):
    """`record["upstream"]["verification"]`をcompletedへ書く合成upstream
    （`evaluate_target`はこれが無いとどのoutputも達成と認めない）。"""
    def handler(context: dict) -> dict:
        calls.append(context["stage_id"])
        return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None,
                "record_updates": {"upstream": {
                    "console_job_path": None, "execution_id": "exec-fake",
                    "verification": {"status": "completed"}}}}
    return handler


def _differential_recorder(calls: list[str]):
    def handler(context: dict) -> dict:
        calls.append(context["stage_id"])
        cid = context["comparison_id"]
        ref = _persist(context["pipeline_root"], f"differential:{cid}")
        return {"status": "succeeded", "result_refs": [ref], "warnings": [], "error": None}
    return handler


def _export_recorder(calls: list[str]):
    def handler(context: dict) -> dict:
        calls.append(context["stage_id"])
        cid = context["comparison_id"]
        refs = [_persist(context["pipeline_root"], f"volcano:{cid}"),
                _persist(context["pipeline_root"], f"tsv:{cid}")]
        return {"status": "succeeded", "result_refs": refs, "warnings": [], "error": None}
    return handler


_HANDLER_KEYS = (
    "prepare_input", "upstream", "validate_outputs", "load_dataset",
    "resolve_metadata", "preprocess", "pca", "resolve_comparisons",
    "differential", "export", "report",
)


def _build_pipeline(tmp_path: Path, *, target: str, comparisons: list | None = None):
    """`save_project=False`を明示する（`evaluate_target`の`gui_project`必須化を
    避ける——本ファイルのfixtureはvalidate_outputsでGUI projectを作らない）。"""
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": target, "comparisons": comparisons or [],
                                            "save_project": False})
    inputs = _minimal_inputs(source_root)
    pipeline_path = create_run(source_root, request, inputs)

    calls: list[str] = []
    handlers = {key: _recorder(calls) for key in _HANDLER_KEYS}
    handlers["upstream"] = _upstream_recorder(calls)
    handlers["preprocess"] = _recorder(calls, output_names=["preprocess"])
    handlers["pca"] = _recorder(calls, output_names=["pca", "pca_figure"])
    handlers["differential"] = _differential_recorder(calls)
    handlers["export"] = _export_recorder(calls)
    handlers["report"] = _recorder(calls, output_names=["quality_report"])
    return pipeline_path, handlers, calls


@pytest.fixture
def run_with_lost_worker(tmp_path):
    """`create_run`/`load_run`/`save_run`だけでstatus=running・worker identity
    設定済みのrunを作る（実PIDを終了させない——自分自身の生きたidentityを使う）。
    """
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    record = load_run(pipeline_path)
    record["status"] = "running"
    record["worker"] = {"identity": process_identity_of_self(), "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return pipeline_path / store_module.RUN_FILENAME


def process_identity_of_self():
    import os
    return process_identity(os.getpid())


# ---------- brief記載のRED（そのまま） ----------

def test_status_does_not_finalize_or_rewrite(run_with_lost_worker, monkeypatch):
    path = run_with_lost_worker
    before = path.read_bytes()
    monkeypatch.setattr("metabolomix.pipeline.recovery.same_process", lambda identity: False)
    status = read_status(path)
    assert status["observed_health"] == "worker_missing"
    assert path.read_bytes() == before
    assert status["status"] != "completed"


# ---------- read_status: 判定不能はdeadと解釈しない / 生存確認 ----------

def test_unknown_worker_identity_is_unknown_not_dead(tmp_path):
    """worker identity自体が記録されていないrunning ——判定不能はdeadと断定しない。"""
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    record = load_run(pipeline_path)
    record["status"] = "running"
    record["worker"] = {"identity": None, "started_at": None}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    status = read_status(pipeline_path)
    assert status["observed_health"] == "unknown"
    assert status["status"] == "running"
    assert status["recovery_hint"]["code"] == "SUPERVISION_UNKNOWN"


def test_real_dead_pid_is_worker_missing_without_mocking_same_process(tmp_path):
    """PID再利用/実プロセス終了: モックせず本物のsame_processで判定する。"""
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    completed = subprocess.run([sys.executable, "-c", "pass"], timeout=30)
    assert completed.returncode == 0
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    identity = process_identity(proc.pid)
    proc.wait(timeout=30)

    record = load_run(pipeline_path)
    record["status"] = "running"
    record["worker"] = {"identity": identity, "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    status = read_status(pipeline_path)
    assert status["observed_health"] == "worker_missing"


def test_status_readable_on_old_completed_run_with_only_results(tmp_path):
    """旧completed job: control/等が無くてもresultsとstatusだけで読める（E04）。"""
    pipeline_path, handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    result = run_engine(pipeline_path, handlers)
    assert result["status"] == "completed"

    status = read_status(pipeline_path)
    assert status["status"] == "completed"
    assert status["observed_health"] == "ok"
    assert status["stage_statuses"]["pca"] == "succeeded"


# ---------- request_cancel: 受理のみ・再送は冪等 ----------

def test_request_cancel_accepts_without_killing_and_is_idempotent(tmp_path):
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")

    first = request_cancel(pipeline_path)
    assert first["accepted"] is True
    assert first["status"] != "cancelled"  # 受理≠確定

    second = request_cancel(pipeline_path)
    assert second["accepted"] is True

    flag = json.loads(cancel_request_path(pipeline_path).read_text(encoding="utf-8"))
    assert flag["cancel_requested"] is True


def test_request_cancel_then_engine_stops_between_stages(tmp_path):
    pipeline_path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    original_preprocess = handlers["preprocess"]

    def preprocess_then_cancel(context):
        outcome = original_preprocess(context)
        request_cancel(pipeline_path)
        return outcome

    handlers["preprocess"] = preprocess_then_cancel
    result = run_engine(pipeline_path, handlers)
    assert result["status"] == "cancelled"
    assert "pca" not in calls


# ---------- prepare_resume: 取消フラグの取り下げ / running のまま残ったrun ----------

def _completed_then(tmp_path: Path, status: str, *, identity=None):
    """一度完走させてから、`status`（と必要ならworker identity）を上書きしたrun。

    上流をsucceededにしておかないと`prepare_resume`が
    `UPSTREAM_RERUN_REQUIRED`で先に止まり、取消フラグの扱いまで到達しない。
    """
    pipeline_path, handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    assert run_engine(pipeline_path, handlers)["status"] == "completed"
    record = load_run(pipeline_path)
    record["status"] = status
    if identity is not None:
        record["worker"] = {"identity": identity, "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return pipeline_path


def test_resume_after_cancel_withdraws_the_cancel_flag(tmp_path):
    """D09: 取消フラグを消さないと、再開したworkerが最初のstage境界で即止まる。"""
    pipeline_path = _completed_then(tmp_path, "cancelled")
    request_cancel(pipeline_path)
    assert cancel_request_path(pipeline_path).is_file()

    receipt = prepare_resume(pipeline_path)

    assert receipt["status"] == "planned"
    assert not cancel_request_path(pipeline_path).exists()


def test_resume_of_a_lost_running_run_returns_to_planned(tmp_path, monkeypatch):
    """A06/§9.3: `read_status`が案内するとおり、worker喪失runはresumeで動く。"""
    pipeline_path = _completed_then(tmp_path, "running", identity=process_identity_of_self())
    monkeypatch.setattr(recovery_module, "same_process", lambda identity: False)
    assert read_status(pipeline_path)["recovery_hint"]["code"] == "PIPELINE_INTERRUPTED"

    receipt = prepare_resume(pipeline_path)

    assert receipt["status"] == "planned"


def test_resume_never_disturbs_a_running_run_whose_worker_is_alive(tmp_path):
    """生きているworkerのrunは`planned`へ戻さず、その取消要求も握り潰さない。

    identityは**このpytestプロセス自身**（確実に生きている）。`same_process`は
    モックしない——pidと生成時刻の両方を見る本物の判定を通す。
    """
    pipeline_path = _completed_then(tmp_path, "running", identity=process_identity_of_self())
    request_cancel(pipeline_path)
    assert read_status(pipeline_path)["observed_health"] == "ok"

    receipt = prepare_resume(pipeline_path)

    assert receipt["status"] == "running", "生きているworkerのrunを再開扱いにした"
    assert cancel_request_path(pipeline_path).is_file(), "出された取消を握り潰した"
    assert load_run(pipeline_path)["status"] == "running"


def test_resume_of_a_running_run_with_an_unresolvable_worker_does_nothing(tmp_path,
                                                                          monkeypatch):
    """判定不能（probe失敗）は「死んでいる」ではない——runにも取消にも触れない。"""
    pipeline_path = _completed_then(tmp_path, "running", identity=process_identity_of_self())
    request_cancel(pipeline_path)

    def _probe_fails(identity):
        raise OSError("probeそのものが失敗した")

    monkeypatch.setattr(recovery_module, "same_process", _probe_fails)
    assert read_status(pipeline_path)["observed_health"] == "unknown"

    receipt = prepare_resume(pipeline_path)

    assert receipt["status"] == "running"
    assert cancel_request_path(pipeline_path).is_file()


def test_resume_keeps_the_cancel_of_a_planned_run_whose_worker_has_not_registered_yet(tmp_path):
    """最終レビュー保留1: 起動直後の窓で、出された取消を黙って消さない。

    `engine.run_engine`がworker identityを刻む（`_record_worker_identity`）のは
    owner lockを取ったあと——起動から数秒間、runは`planned`のままidentityが
    無い。この窓で`_worker_health`は`"unknown"`（identity無し）を返すが、
    それは「死んでいる証拠」ではない（`_RESUMABLE_RUNNING_HEALTH`が
    `"unknown"`を入れないのと同じ理由）。ここでフラグを消すと、これから
    stage loopへ入るworkerは取消を一度も見ずMS-DIALを起動する
    ——利用者が出した取消が黙って無効になる唯一の経路。

    握り潰さない代償は、resumeがもう1往復要ること（起きたworkerが最初の
    stage境界で`cancelled`を確定させ、そのあとのresumeでフラグが消える）
    だけで、恒久的に再開不能にはならない。
    """
    pipeline_path = _completed_then(tmp_path, "planned")
    record = load_run(pipeline_path)
    record["worker"] = {"identity": None, "started_at": None}  # まだ刻まれていない
    save_run(pipeline_path, record, expected_revision=record["state_revision"])
    assert read_status(pipeline_path)["status"] == "planned"

    request_cancel(pipeline_path)
    assert cancel_request_path(pipeline_path).is_file()

    prepare_resume(pipeline_path)

    assert cancel_request_path(pipeline_path).is_file(), \
        "起動直後の窓で、出された取消を握り潰した"

    # 停止が確定（cancelled）したあとのresumeでは、従来どおり取り下げる。
    record = load_run(pipeline_path)
    record["status"] = "cancelled"
    save_run(pipeline_path, record, expected_revision=record["state_revision"])
    prepare_resume(pipeline_path)
    assert not cancel_request_path(pipeline_path).exists()


# ---------- prepare_resume: EXECUTION_UNRESOLVED / UPSTREAM_RERUN_REQUIRED ----------

def test_upstream_not_succeeded_requires_explicit_rerun_upstream(tmp_path):
    """D09: 取消後・timeout後のresumeは上流再実行の明示指定が無いと進めない。"""
    pipeline_path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")

    def failing_upstream(context):
        calls.append(context["stage_id"])
        raise RuntimeError("boom-during-upstream")

    handlers["upstream"] = failing_upstream
    result = run_engine(pipeline_path, handlers)
    assert result["status"] == "failed"  # 有効な出力が無いのでfailed

    with pytest.raises(DomainError, match="UPSTREAM_RERUN_REQUIRED"):
        prepare_resume(pipeline_path)

    # 明示すれば進める(upstream stageがpendingへ戻る)。
    accepted = prepare_resume(pipeline_path, rerun_upstream=True)
    assert accepted["reused"] is False
    assert "upstream" in accepted["reset_stage_ids"]
    record = load_run(pipeline_path)
    assert record["stages"]["upstream"]["status"] == "pending"
    assert record["status"] == "planned"


def test_execution_unresolved_when_console_alive_but_unverified(tmp_path, monkeypatch):
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")
    job_dir = tmp_path / "console_job"
    job_dir.mkdir()
    receipt = {
        "schema": "console-execution.v1", "process_identity": {"pid": 4242, "creation_time": 7},
    }
    (job_dir / "execution-result.json").write_text(json.dumps(receipt), encoding="utf-8")
    # console_job_pathはanalysis-job.jsonそのものへのパス(run_dirはその親)という
    # コードベース全体の流儀(store.register_job_owner等)に合わせる。
    job_json_path = job_dir / "analysis-job.json"
    job_json_path.write_text("{}", encoding="utf-8")

    record = load_run(pipeline_path)
    record["upstream"]["console_job_path"] = str(job_json_path)
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    monkeypatch.setattr(recovery_module, "same_process", lambda identity: True)
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(pipeline_path)


def test_mismatched_console_job_path_convention_is_not_silently_verified(tmp_path):
    """Finding2: console_job_pathが記録されているのに終了証跡が読めない
    （Task18が本タスクの推測と逆の規約を選んだ場合を想定）状況では、
    「console_job_path未記録」と同じ(None, verified)を返して素通りしては
    ならない——upstreamが未verifiedのままEXECUTION_UNRESOLVEDで騒ぐべきで、
    静かにUPSTREAM_RERUN_REQUIRED（または成功）へ落ちてはいけない。
    """
    pipeline_path, _handlers, _calls = _build_pipeline(tmp_path, target="exploratory")

    # 「逆の規約」を再現する: console_job_pathがrun_dirそのものを指しており、
    # 本実装(job_path.parent=run_dir)が読みに行く場所には何も無い。
    job_dir = tmp_path / "console_job"
    job_dir.mkdir()
    (job_dir / "execution-result.json").write_text(
        json.dumps({"schema": "console-execution.v1",
                    "process_identity": {"pid": 4242, "creation_time": 7}}),
        encoding="utf-8")

    record = load_run(pipeline_path)
    # upstreamはまだ未verified(pending)のまま。
    assert record["stages"]["upstream"]["status"] != "succeeded"
    record["upstream"]["console_job_path"] = str(job_dir)  # job.jsonではなくrun_dir自体
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    # 静かに(None, verified=False)へ丸め込まれてUPSTREAM_RERUN_REQUIREDだけが
    # 出るのではなく、規約不一致・判定不能がEXECUTION_UNRESOLVEDとして
    # 見えることを確認する。rerun_upstream=Trueを渡しても素通りしない
    # (既存のEXECUTION_UNRESOLVEDチェック同様、rerun_upstreamでは回避できない)。
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(pipeline_path)
    with pytest.raises(DomainError, match="EXECUTION_UNRESOLVED"):
        prepare_resume(pipeline_path, rerun_upstream=True)


# ---------- D03: 上流を再実行せずcomparisonを解決して完走する ----------

def test_resume_reuses_upstream_and_completes_after_adding_comparisons(tmp_path):
    comparisons = [{"comparison_id": "treated_vs_control",
                    "reference_group": "control", "test_group": "treated"}]
    path, handlers, calls = _build_pipeline(tmp_path, target="differential", comparisons=[])
    handlers["resolve_comparisons"] = _recorder(
        calls, status="needs_input",
        error={"code": "COMPARISON_REQUIRED", "message": "比較を指定してください。", "details": {}})

    waiting = run_engine(path, handlers)
    assert waiting["status"] == "needs_input"

    accepted = prepare_resume(path, updates={"comparisons": comparisons})
    assert accepted["reused"] is False
    assert "resolve_comparisons" in accepted["reset_stage_ids"]

    # resolve_comparisonsは今度こそ本来の(needs_inputを返さない)handlerへ戻す。
    handlers["resolve_comparisons"] = _recorder(calls)

    calls.clear()
    completed = run_engine(path, handlers)
    assert completed["status"] == "completed"
    # upstream/prepare_input/validate_outputsは再実行されない(R20)。
    assert "upstream" not in calls
    assert "prepare_input" not in calls
    assert "validate_outputs" not in calls
    # load_dataset~pcaは再構築のため再度呼ばれる。
    assert "load_dataset" in calls
    assert "preprocess" in calls
    assert "pca" in calls
    assert calls.count("differential:treated_vs_control") == 1
    assert calls.count("export:treated_vs_control") == 1

    record = load_run(path)
    # 探索段階の旧resultは消えていない(append-only)かつ重複追記されていない
    # (R20のcommit_stage_outcome重複防止)。pcaは_ALWAYS_RECONSTRUCT_STAGE_IDSに
    # 属し2回目のrun_engineでも必ずhandlerを呼び直すが、前回と完全に同じ
    # result_refsを返す限り"results"へは1回しか積まれない。件数チェックでない
    # 存在確認だけでは、重複防止を後退させて"pca"のrefが2回積まれても
    # 素通りしてしまい、この回帰を検出できない。
    pca_entries = [r for r in record["results"]
                   if isinstance(r, dict) and r.get("output_name") == "pca"]
    assert len(pca_entries) == 1


# ---------- group-only更新: 変更した比較だけreset ----------

def test_group_only_update_resets_only_the_changed_comparison(tmp_path):
    comparisons = [
        {"comparison_id": "cmp1", "reference_group": "control", "test_group": "treated"},
        {"comparison_id": "cmp2", "reference_group": "control", "test_group": "other"},
    ]
    path, handlers, calls = _build_pipeline(tmp_path, target="differential", comparisons=comparisons)
    result = run_engine(path, handlers)
    assert result["status"] == "completed"

    changed = [
        {"comparison_id": "cmp1", "reference_group": "control", "test_group": "recovered"},
        {"comparison_id": "cmp2", "reference_group": "control", "test_group": "other"},
    ]
    accepted = prepare_resume(path, updates={"comparisons": changed})
    assert set(accepted["reset_stage_ids"]) == {
        "resolve_comparisons", "report", "differential:cmp1", "export:cmp1",
    }
    record = load_run(path)
    assert record["stages"]["differential:cmp2"]["status"] == "succeeded"  # 無関係は不変
    assert record["stages"]["differential:cmp1"]["status"] == "pending"


# ---------- 冪等性: 同一request_id再送 / request_id無し同一更新再送 ----------

def test_same_request_id_resend_does_not_create_new_revision(tmp_path):
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    result = run_engine(path, handlers)
    assert result["status"] == "completed"

    first = prepare_resume(path, updates={"preprocess": {"normalize": "tic"}},
                           request_id="resume-1")
    assert first["reused"] is False
    second = prepare_resume(path, updates={"preprocess": {"normalize": "tic"}},
                            request_id="resume-1")
    assert second["reused"] is True
    assert second["request_revision"] == first["request_revision"]

    record = load_run(path)
    assert record["request"]["revision"] == first["request_revision"]


def test_same_request_id_resend_with_different_updates_conflicts(tmp_path):
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    run_engine(path, handlers)

    prepare_resume(path, updates={"preprocess": {"normalize": "tic"}}, request_id="resume-1")
    with pytest.raises(DomainError, match="IDEMPOTENCY_CONFLICT"):
        prepare_resume(path, updates={"preprocess": {"normalize": "pqn"}}, request_id="resume-1")


def test_same_updates_resent_without_request_id_is_a_no_op(tmp_path):
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    run_engine(path, handlers)

    first = prepare_resume(path, updates={"preprocess": {"normalize": "tic"}})
    second = prepare_resume(path, updates={"preprocess": {"normalize": "tic"}})
    assert second["request_revision"] == first["request_revision"]
    record = load_run(path)
    assert record["request"]["revision"] == first["request_revision"]


def test_no_op_resume_does_not_disturb_a_completed_run(tmp_path):
    """完全に無変化なresume呼出しはcompletedを崩さない(no-opはno-opのまま)。"""
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    result = run_engine(path, handlers)
    assert result["status"] == "completed"

    accepted = prepare_resume(path)  # updates無し・rerun_upstream無し
    assert accepted["reused"] is False
    assert accepted["reset_stage_ids"] == []
    record = load_run(path)
    assert record["status"] == "completed"


# ---------- 二重呼出し: 状態競合を有限回で吸収する ----------

def test_prepare_resume_absorbs_concurrent_state_bump(tmp_path, monkeypatch):
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    run_engine(path, handlers)

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

    accepted = prepare_resume(path, updates={"preprocess": {"normalize": "tic"}})
    assert accepted["reused"] is False
    assert calls_n["n"] >= 2


# ---------- D05: 固定入力・実効メソッドの改変を検出する ----------

def _build_pipeline_with_real_inputs(tmp_path, monkeypatch, *, target="exploratory"):
    """Task13の`inspect_inputs`/`stage_inputs`で実際に入力を固定したrunを作る。

    `stage_inputs`は確定済みのpipeline_rootへしか書けないため、まず仮の
    fingerprintで`create_run`し、得られた実際のpipeline_rootへ配置してから
    `record["inputs"]`を実snapshotへ差し替える。upstream stageはこのfixtureの
    時点で（Task18のhandlerを経由せず）直接`succeeded`にする——D05は「上流は
    既に検証済み」という前提のもとでの再検証なので、ここでは
    UPSTREAM_RERUN_REQUIREDゲートを通過させるためだけに必要。
    """
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    from metabolomix.pipeline.inputs import inspect_inputs, stage_inputs
    from metabolomix.core.atomic_io import canonical_hash

    source = make_source(tmp_path / "source")
    request = resolve_request(source["root"], {"target": target})
    plan = inspect_inputs(source["root"], request, exe_path=source["exe"])

    placeholder_inputs = _minimal_inputs(source["root"], fingerprint="0" * 64)
    pipeline_path = create_run(source["root"], request, placeholder_inputs)

    snapshot = stage_inputs(plan, pipeline_path)
    snapshot["fingerprint"] = canonical_hash({"raw_stat": snapshot["raw_stat"]})

    record = load_run(pipeline_path)
    record["inputs"] = snapshot
    record["stages"]["upstream"]["status"] = "succeeded"
    save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return pipeline_path, source


def test_resume_detects_tampered_effective_method(tmp_path, monkeypatch):
    pipeline_path, _source = _build_pipeline_with_real_inputs(tmp_path, monkeypatch)
    effective = pipeline_path / "inputs" / "effective-method.txt"
    effective.write_text(effective.read_text(encoding="ascii") + "\ntampered", encoding="ascii")

    with pytest.raises(DomainError, match="STAGED_INPUT_MISMATCH"):
        prepare_resume(pipeline_path)


def test_resume_detects_changed_raw_input(tmp_path, monkeypatch):
    pipeline_path, source = _build_pipeline_with_real_inputs(tmp_path, monkeypatch)
    (source["root"] / "S0.wiff").write_text("raw-0-changed", encoding="ascii")

    with pytest.raises(DomainError, match="INPUT_CHANGED"):
        prepare_resume(pipeline_path)


def test_resume_with_rerun_upstream_skips_staged_input_check(tmp_path, monkeypatch):
    """rerun_upstream=Trueのときは再検証しない(上流が自分で再検証・再配置する)。"""
    pipeline_path, source = _build_pipeline_with_real_inputs(tmp_path, monkeypatch)
    (source["root"] / "S0.wiff").write_text("raw-0-changed", encoding="ascii")

    accepted = prepare_resume(pipeline_path, rerun_upstream=True)
    assert accepted["reused"] is False


# ---------- hash不一致: 成果物改変は無条件再利用しない(engine側) ----------

def test_tampered_result_artifact_forces_recompute(tmp_path):
    comparisons = [{"comparison_id": "cmp1", "reference_group": "control", "test_group": "treated"}]
    path, handlers, calls = _build_pipeline(tmp_path, target="differential", comparisons=comparisons)

    output_file = path / "report.txt"

    def report_handler(context):
        calls.append(context["stage_id"])
        output_file.write_text("report-v1", encoding="utf-8")
        import hashlib
        digest = hashlib.sha256(output_file.read_bytes()).hexdigest()
        return {"status": "succeeded",
                "result_refs": [{"output_name": "quality_report", "result_id": "report-result",
                                 "relative_path": "report.txt", "hash": digest}],
                "warnings": [], "error": None}

    handlers["report"] = report_handler
    result = run_engine(path, handlers)
    assert result["status"] == "completed"

    # 成果物を改変(hash不一致)してから再度run_engineを呼ぶ。
    output_file.write_text("tampered", encoding="utf-8")
    calls.clear()
    result2 = run_engine(path, handlers)
    assert result2["status"] == "completed"
    assert "report" in calls  # hash不一致を検出し再実行された


def test_untampered_result_artifact_is_skipped_on_second_run(tmp_path):
    comparisons = [{"comparison_id": "cmp1", "reference_group": "control", "test_group": "treated"}]
    path, handlers, calls = _build_pipeline(tmp_path, target="differential", comparisons=comparisons)
    result = run_engine(path, handlers)
    assert result["status"] == "completed"

    calls.clear()
    # 何も変えずにもう一度run_engineを呼ぶ(resumeの二重呼出し相当)。
    result2 = run_engine(path, handlers)
    assert result2["status"] == "completed"
    assert "upstream" not in calls
    assert "differential:cmp1" not in calls
    assert "export:cmp1" not in calls
    assert "report" not in calls
    # runtime再構築対象は再度呼ばれる。
    assert "load_dataset" in calls
    assert "pca" in calls



# ---------- 同じシートを書き直した訂正を検出する（レビュー指摘2） ----------

def _completed_run_with_manifest(tmp_path):
    """比較まで完走したrunと、`resolve_metadata`が記録するシート指紋を用意する。

    本番では`service._handle_resolve_metadata`が
    `record["inputs"]["manifest_source"]`を書く。ここではその形だけを合成する
    （このファイルの流儀どおり、handlerは注入fixtureのまま）。
    """
    from metabolomix.pipeline import inputs as inputs_mod

    comparisons = [{"comparison_id": "cmp1", "reference_group": "control",
                    "test_group": "treated"}]
    path, handlers, calls = _build_pipeline(tmp_path, target="differential",
                                            comparisons=comparisons)
    source_root = tmp_path / "source"
    manifest = source_root / inputs_mod.DEFAULT_MANIFEST_NAME
    manifest.write_text("# schema = sample-manifest.v1\nsample_id\tgroup\nS0\tcontrol\n",
                        encoding="utf-8", newline="\n")

    assert run_engine(path, handlers)["status"] == "completed"

    record = load_run(path)
    record["inputs"] = {**(record.get("inputs") or {}),
                        "manifest_source": inputs_mod.manifest_source_record(
                            source_root, {"sample_manifest": None})}
    save_run(path, record, expected_revision=record["state_revision"])
    calls.clear()
    return path, handlers, calls, manifest


def test_rewriting_the_same_manifest_resets_the_downstream_stages(tmp_path):
    """パスは同じまま中身だけ直したシートを、無視せず差し戻す。

    要求の`sample_manifest`はパス文字列でしかないので、パスの比較だけでは
    「書き直した」が変更として現れない——旧群割当のまま`completed`が保たれ、
    訂正が無視されたことすら利用者に伝わらない。
    """
    path, handlers, calls, manifest = _completed_run_with_manifest(tmp_path)

    manifest.write_text("# schema = sample-manifest.v1\nsample_id\tgroup\nS0\ttreated\n",
                        encoding="utf-8", newline="\n")
    receipt = prepare_resume(path)

    assert receipt["status"] == "planned"          # 再開待ちへ戻っている
    assert "differential:cmp1" in receipt["reset_stage_ids"]
    assert "export:cmp1" in receipt["reset_stage_ids"]
    assert "report" in receipt["reset_stage_ids"]
    # 実際に再開すると、その3工程がやり直される。
    assert run_engine(path, handlers)["status"] == "completed"
    assert "differential:cmp1" in calls and "export:cmp1" in calls


def test_an_untouched_manifest_keeps_the_completed_run_as_is(tmp_path):
    """内容が変わっていないシートでは、no-op resumeをcompletedのままにする。"""
    path, _handlers, _calls, _manifest = _completed_run_with_manifest(tmp_path)

    receipt = prepare_resume(path)

    assert receipt["status"] == "completed"
    assert receipt["reset_stage_ids"] == []


# ---------- make_context (R19) ----------

def test_make_context_carries_full_request_and_upstream_identity(tmp_path):
    path, handlers, calls = _build_pipeline(tmp_path, target="exploratory")
    seen = {}

    def capturing_upstream(context):
        seen["context"] = context
        calls.append(context["stage_id"])
        return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None}

    handlers["upstream"] = capturing_upstream
    run_engine(path, handlers)

    context = seen["context"]
    assert context["request"]["schema"] == "pipeline-request.v1"
    assert "value_sources" in context["request"]
    assert "inputs" in context and "upstream" in context and "identity" in context
    assert context["identity"]["pipeline_root"]
    assert context["request_meta"]["revision"] == 1


# ---------- v2 の上流stage ID（既存バグの再現） ----------

def test_v2_upstream_verified_state_uses_the_v2_stage_id(tmp_path):
    """v2 record の上流stageは`execute_console`で、`upstream`ではない。

    `verified`を`stages["upstream"]`から読むと、v2では常にFalseになる。その状態で
    終了証跡が読めないと、`execute_console`が succeeded として永続化されていても
    `EXECUTION_UNRESOLVED`を投げ、`prepare_resume`が止まる——同モジュールの
    `prepare_resume`は既に`_upstream_stage_id`を使っており、ここだけが
    取り残されていた。
    """
    from metabolomix.pipeline import recovery

    job_path = tmp_path / "console" / "attempt-0001" / "analysis-job.json"
    job_path.parent.mkdir(parents=True)
    job_path.write_text("{}", encoding="utf-8")
    record = {
        "schema": "pipeline-run.v2",
        "stages": {"execute_console": {"status": "succeeded"}},
        "upstream": {"console_job_path": str(job_path)},
    }

    # 終了証跡は置かない。verifiedを読み違えるとここで例外になる。
    identity, verified = recovery._console_supervision_state(record)
    assert verified is True
    assert identity is None


def test_v1_upstream_verified_state_is_unchanged(tmp_path):
    from metabolomix.pipeline import recovery

    job_path = tmp_path / "console" / "attempt-0001" / "analysis-job.json"
    job_path.parent.mkdir(parents=True)
    job_path.write_text("{}", encoding="utf-8")
    record = {
        "schema": "pipeline-run.v1",
        "stages": {"upstream": {"status": "succeeded"}},
        "upstream": {"console_job_path": str(job_path)},
    }
    identity, verified = recovery._console_supervision_state(record)
    assert verified is True
    assert identity is None
