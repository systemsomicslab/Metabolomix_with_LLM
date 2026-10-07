"""metabolomix.pipeline.service の単体テスト。

`tests/test_pipeline_engine.py` / `test_pipeline_recovery.py` は合成handlerで
engine自体（stage順序・owner lock・冪等性）を検証しており、`build_handlers()`が
組み立てる実handler一式そのものは別の場所で検証されていなかった。ここでは
そのうち、carry-forwardルーリングとして明示された振る舞いに絞って検証する。

- R18: PCA不成立(`PreconditionError`)がhandler境界で`needs_input`へ変換されること。
- Console用workerを二重起動しないこと（`upstream` handlerが `console.worker` を
  一切呼ばず `console.execution.supervise` だけを呼ぶこと）。
- `launch_pipeline_worker` が `sys.executable` を使い、cwdをこのcheckout自身へ
  固定すること（Task2: 利用者のcwdや別checkoutを暗黙に使わない）。
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

from metabolomix.analysis.dataset_analysis import PreconditionError
from metabolomix.pipeline import engine, service
from metabolomix.pipeline.request import resolve_request
from metabolomix.pipeline.store import create_run, load_run, save_run


def _pca_context(pipeline_root: Path) -> dict:
    return {
        "runtime": {"dataset": object()},  # pca_dataset自体をpatchするので中身は使わない
        "pipeline_root": pipeline_root,
        "request_meta": {"revision": 1},
    }


def test_pca_precondition_error_becomes_needs_input_stage_result(tmp_path, monkeypatch):
    """R18: PreconditionErrorはhandler境界でneeds_inputへ変換され、failedにならない。"""
    def raise_precondition(*a, **kw):
        raise PreconditionError("missing_state", "検出状態が無くPCAを実行できません。",
                                {"detail": "no-detection-state"}, state="preprocessed_matrix")

    monkeypatch.setattr(service.dataset_service, "pca_dataset", raise_precondition)

    outcome = service._handle_pca(_pca_context(tmp_path))

    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] == "PCA_PRECONDITION_MISSING"
    assert outcome["result_refs"] == []


def test_pca_success_persists_data_and_figure_refs(tmp_path, monkeypatch):
    """成功時はpca/pca_figureの2件のoutput_name付きrefをresult_refsへ積む。"""
    fake_result = {
        "provenance": {"result_id": "res_pca_fake", "input_fingerprint": "f" * 40,
                      "parent_ids": [], "warnings": []},
        "scores": [[0.1, 0.2]], "loadings": [[0.3, 0.4]], "explained_variance": [0.9],
    }
    def fake_save_figure(ds, result, figure_path, *, kind):
        Path(figure_path).parent.mkdir(parents=True, exist_ok=True)
        Path(figure_path).write_bytes(b"fake-png")

    monkeypatch.setattr(service.dataset_service, "pca_dataset", lambda *a, **kw: fake_result)
    monkeypatch.setattr(service.result_output, "save_result_figure", fake_save_figure)

    outcome = service._handle_pca(_pca_context(tmp_path))

    assert outcome["status"] == "succeeded"
    names = {ref["output_name"] for ref in outcome["result_refs"]}
    assert names == {"pca", "pca_figure"}


def test_upstream_handler_never_calls_console_worker_module(tmp_path, monkeypatch):
    """Console用workerを二重起動しない: `metabolomix.console.worker`を一切呼ばない。"""
    import metabolomix.console.worker as console_worker_mod

    def _forbidden(*a, **kw):
        raise AssertionError("upstream handlerがmetabolomix.console.workerを呼び出した"
                            "（Console用workerの二重起動）")

    for name in ("run_job", "launch_console_worker"):
        if hasattr(console_worker_mod, name):
            monkeypatch.setattr(console_worker_mod, name, _forbidden)

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
    monkeypatch.setattr(service.store, "register_job_owner", lambda job_path, root: None)
    monkeypatch.setattr(service.console_execution, "write_supervision_inputs",
                        lambda run_dir, payload: None)

    pipeline_root = tmp_path / "pipeline"
    pipeline_root.mkdir()
    (pipeline_root / "input").mkdir()
    context = {
        "pipeline_root": pipeline_root,
        "pipeline_id": "pipe1",
        "request": {"measure": "peak_height", "save_project": False, "timeout_s": 60},
        "inputs": {"method": {"effective_relative_path": None, "source_path": "C:/m.txt"},
                  "polarity": {"value": "negative"}, "exe": {"path": "C:/fake.exe"}},
    }

    outcome = service._handle_upstream(context)

    assert len(supervise_calls) == 1
    assert outcome["status"] == "succeeded"


def test_launch_pipeline_worker_uses_sys_executable_and_repo_root_cwd(tmp_path, monkeypatch):
    """Task2: launchはsys.executableで、cwdは呼び出し元のcwdでなくこのcheckout自身。"""
    captured = {}

    def fake_launch_detached(command, *, cwd, log_path):
        captured["command"] = command
        captured["cwd"] = cwd
        captured["log_path"] = log_path
        return {"launched": True, "pid": 4321}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    pipeline_path = tmp_path / "pipeline"
    (pipeline_path / "control").mkdir(parents=True)

    result = service.launch_pipeline_worker(pipeline_path)

    assert result == {"launched": True, "pid": 4321}
    assert captured["command"][0] == sys.executable
    assert captured["command"][1:4] == ["-m", "metabolomix.pipeline.worker", "--pipeline"]
    assert captured["command"][4] == str(pipeline_path)
    assert captured["cwd"] == service._REPO_ROOT


# ---------- レビュー指摘5: CONFOUNDED_COMPARISONはneeds_input（failedにしない）----------

def test_confounded_comparison_becomes_needs_input_not_failed(monkeypatch):
    """群とバッチが完全交絡した比較は、allow_confounded=trueという利用者入力を
    待つneeds_inputであって、failedではない（spec §7.4「入力を待つ」、R18と同じ
    defect class）。`_handle_resolve_comparisons`自身はDomainErrorを投げるだけで、
    needs_inputへの変換はengineの`_invoke_handler`が`_NEEDS_INPUT_CODES`経由で行う
    ——両方を実際に呼んで確認する。
    """
    comparisons = [{"comparison_id": "treated_vs_control",
                   "reference_group": "control", "test_group": "treated"}]

    def fake_resolve_comparison(ds, comparison, metadata):
        return {"reference_group": comparison["reference_group"],
                "test_group": comparison["test_group"],
                "confounding": {"confounded": True}, "allow_confounded": False}

    monkeypatch.setattr(service.dataset_service, "resolve_comparison", fake_resolve_comparison)
    context = {"runtime": {"dataset": object(), "metadata": []},
              "request": {"comparisons": comparisons}}

    outcome = engine._invoke_handler(service._handle_resolve_comparisons, context)

    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] == "CONFOUNDED_COMPARISON"


# ---------- 最終レビュー指摘2: resolve_metadataが解決した行をrecordへ残す ----------

def test_resolve_metadata_records_the_resolved_manifest_into_inputs():
    """spec §7.3/§11: 解決したrole/group/batch/orderの出所をrecordへ残す。

    `record["inputs"]["manifest"]`は`metabolomix.pipeline.report._section_sample_
    provenance`が読み、`store`が相対化し、`recovery`が絶対化する——**読む側は
    3つあるのに書く側が居なかった**。`_handle_resolve_metadata`が
    `record_updates`を返さない限り、解決結果は`context["runtime"]`にしか
    残らず、品質レポートのサンプル来歴は毎回
    `(no sample manifest recorded)`になる。
    """
    from tests.pipeline_fixtures import make_dataset

    ds = make_dataset()
    context = {
        "runtime": {"dataset": ds},
        "request": {"sample_manifest": None},
        "identity": {"source_root": str(Path.cwd())},
        "inputs": {"source_root": "src", "method": {"source_path": "inputs/method.txt"}},
    }

    outcome = service._handle_resolve_metadata(context)

    assert outcome["status"] == "succeeded"
    manifest = outcome["record_updates"]["inputs"]["manifest"]
    assert [row["sample_id"] for row in manifest] == ds.sample_names
    # 既存の入力スナップショットを取りこぼさない（manifestを足すだけ）。
    assert outcome["record_updates"]["inputs"]["method"] == {"source_path": "inputs/method.txt"}
    # レポートが読む8列がそのまま揃っている。
    for column in ("sample_id", "source_file", "role", "group", "batch",
                   "injection_order", "qc_pool", "include"):
        assert column in manifest[0], column


# ---------- レビュー指摘4: quality_report自身の未生成を理由にした偽陽性警告 ----------

def _minimal_inputs_for(root: Path) -> dict:
    return {"source_root": str(root), "fingerprint": "f" * 64, "raw_inventory": []}


def _build_exploratory_run(tmp_path: Path) -> Path:
    """探索目標のrunを作り、upstream検証済み・preprocess/pca/pca_figureの3件を
    実ref（`persist_result`が返す実file・実hash付き）で先に達成させておく。
    quality_report以外は既に揃っている状態を作るのが目的。
    """
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": "exploratory", "save_project": False})
    pipeline_root = create_run(source_root, request, _minimal_inputs_for(source_root))

    record = load_run(pipeline_root)
    record = copy.deepcopy(record)
    record["upstream"] = {"console_job_path": None, "execution_id": "exec-1",
                          "verification": {"status": "completed"}}
    refs = [service.report_mod.persist_result(pipeline_root, {
        "output_name": name, "kind": "synthetic", "result_id": f"res_{name}",
        "data": {"synthetic": True, "name": name},
    }) for name in ("preprocess", "pca", "pca_figure")]
    record["results"] = refs
    save_run(pipeline_root, record, expected_revision=record["state_revision"])
    return pipeline_root


def test_report_handler_does_not_warn_when_only_its_own_ref_was_missing(tmp_path):
    """`quality_report`以外はすべて達成済みのrunでは、reportハンドラの警告
    `REPORT_INCOMPLETE_AT_WRITE_TIME`は出ない——quality_report自身がまだ
    `results`に無い時点の一度きりのevaluate_targetだけを見ると常にTrueになって
    しまう（レビュー指摘4）。
    """
    pipeline_root = _build_exploratory_run(tmp_path)
    record = load_run(pipeline_root)
    context = {"pipeline_root": pipeline_root,
              "identity": {"pipeline_id": record["identity"]["pipeline_id"]}}

    outcome = service._handle_report(context)

    assert outcome["warnings"] == []


def test_report_handler_still_warns_when_something_else_is_genuinely_missing(tmp_path):
    """quality_report以外にも本当に未達成の出力が残っている場合は、引き続き
    警告が出る（Finding4の修正が「常にFalse」へ倒れていないことの反証）。
    """
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": "exploratory", "save_project": False})
    pipeline_root = create_run(source_root, request, _minimal_inputs_for(source_root))
    record = load_run(pipeline_root)
    record = copy.deepcopy(record)
    record["upstream"] = {"console_job_path": None, "execution_id": "exec-1",
                          "verification": {"status": "completed"}}
    # pcaを欠いたまま(preprocessだけ達成)。
    ref = service.report_mod.persist_result(pipeline_root, {
        "output_name": "preprocess", "kind": "synthetic", "result_id": "res_preprocess",
        "data": {"synthetic": True},
    })
    record["results"] = [ref]
    save_run(pipeline_root, record, expected_revision=record["state_revision"])
    context = {"pipeline_root": pipeline_root,
              "identity": {"pipeline_id": record["identity"]["pipeline_id"]}}

    outcome = service._handle_report(context)

    codes = {w["code"] for w in outcome["warnings"]}
    assert "REPORT_INCOMPLETE_AT_WRITE_TIME" in codes


# ---------- レビュー指摘3: resume_pipelineはpipeline-run.json自身も受け付ける ----------

def test_resume_pipeline_normalizes_the_run_json_path_before_launching(tmp_path, monkeypatch):
    """`pipeline_status`/`pipeline_cancel`と同じく、`resume_pipeline`も
    `pipeline-run.json`そのものへのパスを受け付ける（`pipeline_tools.py`が
    アドバタイズする形）。`prepare_resume`は内部で正規化するが、`resume_pipeline`
    自身がその後の`launch_pipeline_worker`/receiptで生のパスを使うと、
    二重に`pipeline-run.json`を連結した不正なパスへ`--pipeline`を渡し、
    `_await_launch_handshake`の`store.load_run`が`PIPELINE_RUN_NOT_FOUND`で
    落ちる（レビュー指摘3）。
    """
    from metabolomix.pipeline import store as store_mod

    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": "exploratory", "save_project": False})
    pipeline_root = store_mod.create_run(source_root, request, _minimal_inputs_for(source_root))
    run_json_path = pipeline_root / store_mod.RUN_FILENAME

    # 上流は既に検証済みという前提にする(D05と同じ流儀)——本テストの対象は
    # 「パス正規化」だけで、UPSTREAM_RERUN_REQUIREDの判定は無関係。
    record = load_run(pipeline_root)
    record = copy.deepcopy(record)
    record["stages"]["upstream"]["status"] = "succeeded"
    save_run(pipeline_root, record, expected_revision=record["state_revision"])

    captured = {}

    def fake_launch_detached(command, *, cwd, log_path):
        captured["command"] = command
        captured["log_path"] = log_path
        return {"launched": True, "pid": 555}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    receipt = service.resume_pipeline(run_json_path)

    assert receipt["pipeline_path"] == str(pipeline_root)
    assert captured["command"][-1] == str(pipeline_root)
    assert captured["log_path"] == pipeline_root / store_mod.CONTROL_SUBDIR / "worker-launch.log"


# ---------- レビュー指摘1: 再利用した活動中/終端runの上に二重起動しない ----------

def _prepare_source(tmp_path, monkeypatch):
    """`start_pipeline`の実処理(受付冪等性込み)を通すための実生データフォルダ。

    `MSDIAL_EXE`環境変数と`is_console_exe`だけを差し替え、他はすべて
    `inspect_inputs`/`find_or_create_run`本体を実際に通す
    （`tests/test_pipeline_recovery.py::_build_pipeline_with_real_inputs`と
    同じ流儀）。
    """
    from tests.pipeline_fixtures import make_source
    source = make_source(tmp_path / "source")
    monkeypatch.setenv("MSDIAL_EXE", str(source["exe"]))
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    monkeypatch.setenv("LIPIDMIX_PIPELINE_INDEX_DIR", str(tmp_path / "_pipeline_index_base"))
    return source


def test_start_pipeline_does_not_relaunch_over_an_active_reused_run(tmp_path, monkeypatch):
    """`find_or_create_run`が同一fingerprint/要求の既存runを再利用したとき、
    そのrunが既に活動中(running)なら二度目のworkerを起動しない
    （レビュー指摘1: 二重起動は片方が`PIPELINE_ALREADY_RUNNING`で無駄死にし、
    ログとプロセスが1つずつ漏れる）。
    """
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls: list = []

    def fake_launch_detached(command, *, cwd, log_path):
        launch_calls.append(command)
        return {"launched": True, "pid": 1000 + len(launch_calls)}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    first = service.start_pipeline(source["root"])
    assert first["launched"] is True
    assert len(launch_calls) == 1

    # 実際に別workerがこのrunを引き受けた状態を模す(statusをrunningへ進める)。
    pipeline_path = Path(first["pipeline_path"])
    record = load_run(pipeline_path)
    record = copy.deepcopy(record)
    record["status"] = "running"
    record["worker"] = {"identity": {"pid": 424242, "creation_time": 1},
                        "started_at": "2026-09-06T00:00:00+00:00"}
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    second = service.start_pipeline(source["root"])

    assert second["pipeline_path"] == first["pipeline_path"]
    assert second["launched"] is False
    assert len(launch_calls) == 1  # 再launchしていない


def test_start_pipeline_does_not_relaunch_over_a_reused_completed_run(tmp_path, monkeypatch):
    """completed runの再利用でも同様に起動しない(`_ALWAYS_RECONSTRUCT_STAGE_IDS`
    の無条件再計算・figureの上書きを避ける)。"""
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls: list = []

    def fake_launch_detached(command, *, cwd, log_path):
        launch_calls.append(command)
        return {"launched": True, "pid": 2000 + len(launch_calls)}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.1)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    first = service.start_pipeline(source["root"])
    pipeline_path = Path(first["pipeline_path"])
    record = load_run(pipeline_path)
    record = copy.deepcopy(record)
    record["status"] = "completed"
    save_run(pipeline_path, record, expected_revision=record["state_revision"])

    second = service.start_pipeline(source["root"])

    assert second["launched"] is False
    assert len(launch_calls) == 1


# ---------- レビュー指摘2: handshakeは「新しいidentity」だけをconfirmedとする ----------

def test_resume_handshake_is_not_confirmed_when_worker_identity_is_unchanged(tmp_path, monkeypatch):
    """resumeで前回workerの残骸identityが既に記録されている場合、今回のlaunchが
    新しいidentityを書かない(=実質何も起動していない)限りconfirmedにならない
    ——`identity is not None`だけを見ると即座にconfirmedへ誤判定する
    （レビュー指摘2）。
    """
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": "exploratory", "save_project": False})
    pipeline_root = create_run(source_root, request, _minimal_inputs_for(source_root))

    stale_identity = {"pid": 111111, "creation_time": 1}
    record = load_run(pipeline_root)
    record = copy.deepcopy(record)
    record["stages"]["upstream"]["status"] = "succeeded"
    record["worker"] = {"identity": stale_identity, "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(pipeline_root, record, expected_revision=record["state_revision"])

    # launchは起きるが、新しいidentityは一切書かない(=起動が実質失敗した状況を模す)。
    monkeypatch.setattr(service, "launch_detached",
                        lambda command, *, cwd, log_path: {"launched": True, "pid": 999})
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.2)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    receipt = service.resume_pipeline(pipeline_root)

    assert receipt["launch"]["handshake"] == "not_confirmed"


def test_resume_handshake_confirms_on_a_genuinely_new_identity(tmp_path, monkeypatch):
    """対照: 新しいidentityが実際に書かれればconfirmedになる(過剰検出でないこと)。"""
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {"target": "exploratory", "save_project": False})
    pipeline_root = create_run(source_root, request, _minimal_inputs_for(source_root))

    stale_identity = {"pid": 111111, "creation_time": 1}
    record = load_run(pipeline_root)
    record = copy.deepcopy(record)
    record["stages"]["upstream"]["status"] = "succeeded"
    record["worker"] = {"identity": stale_identity, "started_at": "2026-09-05T00:00:00+00:00"}
    save_run(pipeline_root, record, expected_revision=record["state_revision"])

    def fake_launch_detached(command, *, cwd, log_path):
        current = load_run(pipeline_root)
        current = copy.deepcopy(current)
        current["worker"] = {"identity": {"pid": 222222, "creation_time": 2},
                             "started_at": "2026-09-06T00:00:00+00:00"}
        save_run(pipeline_root, current, expected_revision=current["state_revision"])
        return {"launched": True, "pid": 222222}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.5)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    receipt = service.resume_pipeline(pipeline_root)

    assert receipt["launch"]["handshake"] == "confirmed"


# ---------- レビュー指摘6: 受付層のさらなる被覆 ----------

def test_plan_pipeline_never_launches_a_worker(tmp_path, monkeypatch):
    """`plan_pipeline`は検査・保存のみ——`launch_detached`を一切呼ばない。"""
    source = _prepare_source(tmp_path, monkeypatch)

    def poison(*a, **kw):
        raise AssertionError("plan_pipelineがworkerを起動した")

    monkeypatch.setattr(service, "launch_detached", poison)

    receipt = service.plan_pipeline(source["root"])

    assert receipt["launched"] is False
    assert receipt["status"] == "planned"
    record = load_run(Path(receipt["pipeline_path"]))
    assert record["status"] == "planned"


def test_plan_receipt_shows_the_resolved_settings_the_instructions_promise(tmp_path, monkeypatch):
    """最終レビュー指摘6: `MCP_INSTRUCTIONS`は`pipeline_plan`を
    「method file / LBM / polarityの解決結果を起動前に確認する」入口として
    案内している。receiptがそれを持っていなければ、案内された確認は
    `pipeline_status(include_details=True)`でrecord全体を引っぱるしかない
    ——コンパクトなreceiptという制約と正面から衝突する。

    同時に**肥大させない**ことも縛る（CLAUDE.md「戻り値を肥大させない」）:
    raw一覧・stat・entries・overridesのような大きな中間データは載せない。
    """
    source = _prepare_source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source["root"])

    resolved = receipt["resolved"]
    assert resolved["method"]["source_path"] == str(Path(source["method"]).resolve())
    assert len(resolved["method"]["sha256"]) == 64
    assert resolved["lbm"]["path"] == str(Path(source["lbm"]).resolve())
    assert resolved["polarity"]["value"] in ("positive", "negative")
    assert resolved["polarity"]["source"]

    # コンパクトさ: 大きな中間データを持ち込んでいない。
    assert set(resolved) == {"method", "lbm", "polarity", "library_mode", "msp"}
    assert resolved["library_mode"] == "auto"
    assert resolved["msp"] == {"file": None, "sha256": None, "source": "not_used"}
    assert set(resolved["method"]) == {"source_path", "sha256"}
    text = json.dumps(receipt, ensure_ascii=False)
    assert "raw_stat" not in text and "entries" not in text
    assert len(text) < 1500, len(text)

    # 起動する側（pipeline_run）のreceiptは従来どおり最小のまま。
    assert "resolved" not in service._dispatch_receipt(
        Path(receipt["pipeline_path"]), launched=False)


def test_start_pipeline_saves_needs_input_without_launching_on_bad_manifest(tmp_path, monkeypatch):
    """`_mark_needs_input_without_launch`の実経路: 壊れたsample_manifestは
    Console起動前に検出され、runは`needs_input`のまま保存されるが
    workerは一切起動されない(spec §9.2)。
    """
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


def test_start_pipeline_handshake_timeout_does_not_relaunch(tmp_path, monkeypatch):
    """handshakeが上限まで確認できなくても、再launchはしない(brief拘束)。"""
    source = _prepare_source(tmp_path, monkeypatch)
    launch_calls: list = []

    def fake_launch_detached(command, *, cwd, log_path):
        # 現実のworkerが実際に起動する場合と違い、ここではrecordに一切触れない
        # (=identityが書かれない)——起動が実質確認できない状況を模す。
        launch_calls.append(command)
        return {"launched": True, "pid": 777}

    monkeypatch.setattr(service, "launch_detached", fake_launch_detached)
    monkeypatch.setattr(service, "_HANDSHAKE_TIMEOUT_S", 0.15)
    monkeypatch.setattr(service, "_HANDSHAKE_POLL_S", 0.02)

    receipt = service.start_pipeline(source["root"])

    assert receipt["launch"]["handshake"] == "not_confirmed"
    assert len(launch_calls) == 1  # timeoutしても再launchしない
    assert receipt["status"] == "planned"  # 起動失敗と決め付けてstatusを書き換えない


# ---------- レビュー指摘6: handler一式が本物として連鎖すること ----------

def test_worker_run_worker_completes_a_real_exploratory_pipeline(tmp_path, monkeypatch):
    """`build_handlers()`を`run_engine`へ通す経路を、モックした本体ではなく
    `worker.run_worker`を実際に呼んで検証する（レビュー指摘6）。

    Consoleの実バイナリだけは起動不可能なので、`console_runner.build_msdial_cmd`
    （既存のConsole層テストが使う唯一の注入口。`_handle_upstream`自体や
    `build_handlers`はモックしない）を実データを書く偽コマンドへ差し替える。
    それ以外——`prepare_input`の実配置、`upstream`の実`supervise`監視、
    `validate_outputs`、実mzTabからの`load_dataset`、`resolve_metadata`の
    自動推定、実`preprocess_auto`、実`pca_dataset`、`report`の
    `evaluate_target`——はすべて本物を通す。
    """
    from metabolomix.pipeline import worker
    from tests.pipeline_fixtures import fake_console_command, make_source, mztab_text, use_fake_console

    source = _prepare_source(tmp_path, monkeypatch)
    assert source is not None

    plan_receipt = service.plan_pipeline(
        source["root"], {"target": "exploratory", "save_project": False})
    pipeline_path = Path(plan_receipt["pipeline_path"])

    # stage_inputsが配置する先(pipeline_root/input)のファイル名は元名を保つ
    # (metabolomix.pipeline.inputs.stage_inputs)。まだ配置されていない時点でも、
    # この規約からpost-stage時点の絶対パスを先に計算できる。
    staged_sources = [pipeline_path / "input" / f"S{i}.wiff" for i in range(8)]
    mztab_content = mztab_text(tmp_path, staged_sources)

    # Console実行はattemptごとに別run_dirを持つ（最終レビュー指摘3）。この
    # runは初回なのでattempt-0001。
    run_dir = pipeline_path / "console" / "attempt-0001"
    mztab_out = run_dir / "msdial" / "Height_AlignmentResult_1.mzTab"
    use_fake_console(monkeypatch, fake_console_command({mztab_out: mztab_content}))

    result = worker.run_worker(pipeline_path)

    assert result["status"] == "completed"
    stage_statuses = {sid: s["status"] for sid, s in result["stages"].items()}
    assert stage_statuses["prepare_input"] == "succeeded"
    assert stage_statuses["upstream"] == "succeeded"
    assert stage_statuses["validate_outputs"] == "succeeded"
    assert stage_statuses["load_dataset"] == "succeeded"
    assert stage_statuses["resolve_metadata"] == "succeeded"
    assert stage_statuses["preprocess"] == "succeeded"
    assert stage_statuses["pca"] == "succeeded"
    assert stage_statuses["report"] == "succeeded"
    # Finding4の回帰確認: 本物の完走でもREPORT_INCOMPLETE_AT_WRITE_TIMEの
    # 偽陽性警告が出ない。
    codes = {w["code"] for w in result.get("warnings") or []}
    assert "REPORT_INCOMPLETE_AT_WRITE_TIME" not in codes


# ---------- PreconditionErrorはneeds_inputへ（レビュー指摘P2） ----------

def test_a_precondition_error_from_any_handler_becomes_needs_input(monkeypatch):
    """数値層の`PreconditionError`は`DomainError`ではないので、素通りさせると
    engineの汎用例外分岐に落ち、`error.code`がPythonのクラス名のまま`failed`に
    なる——入力を直せば再開できる停止なのに「回復不能な失敗」に見える。
    """
    from metabolomix.analysis.dataset_analysis import PreconditionError
    from metabolomix.pipeline import service as service_mod

    def _boom(context):
        raise PreconditionError(
            "bad_request",
            "この DatasetState は検出状態を持たないため min_detection_rate を適用できません。",
            {"min_detection_rate": 1.0})

    monkeypatch.setattr(service_mod, "_handle_preprocess", _boom)
    outcome = service_mod.build_handlers()["preprocess"]({"stage_id": "preprocess"})

    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] == "ANALYSIS_PRECONDITION_INVALID"
    assert outcome["error"]["details"]["min_detection_rate"] == 1.0


def test_a_missing_state_precondition_carries_the_missing_state_name(monkeypatch):
    from metabolomix.analysis.dataset_analysis import PreconditionError
    from metabolomix.pipeline import service as service_mod

    def _boom(context):
        raise PreconditionError("missing_state", "前処理済み行列がありません。",
                                state="dataset_preprocessed")

    monkeypatch.setattr(service_mod, "_handle_pca", _boom)
    outcome = service_mod.build_handlers()["pca"]({"stage_id": "pca"})

    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] == "ANALYSIS_PRECONDITION_MISSING"
    assert outcome["error"]["details"]["state"] == "dataset_preprocessed"


def test_a_precondition_error_stops_the_engine_as_needs_input_not_failed(tmp_path, monkeypatch):
    """engine経由でも`failed`にならないこと（`_invoke_handler`の汎用分岐を通さない）。"""
    from metabolomix.analysis.dataset_analysis import PreconditionError
    from metabolomix.pipeline import service as service_mod
    from metabolomix.pipeline.engine import _invoke_handler

    def _boom(context):
        raise PreconditionError("bad_request", "レシピが不正です。", {"recipe": {}})

    monkeypatch.setattr(service_mod, "_handle_preprocess", _boom)
    outcome = _invoke_handler(service_mod.build_handlers()["preprocess"],
                              {"stage_id": "preprocess"})

    assert outcome["status"] == "needs_input"
    assert outcome["error"]["code"] != "PreconditionError"


def test_plan_receipt_shows_the_lab_msp_without_its_location(tmp_path, monkeypatch):
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_only"})
    resolved = receipt["resolved"]
    assert resolved["library_mode"] == "msp_only"
    assert resolved["lbm"] == {"path": None, "sha256": None}
    assert resolved["msp"]["file"] == "lab_neg.msp"
    assert len(resolved["msp"]["sha256"]) == 64
    assert resolved["msp"]["source"] == "config_file"
    assert "secret_share" not in json.dumps(receipt, ensure_ascii=False)


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
    """最終レビュー #3: pipeline_status(include_details=True) の record から研究室 MSP の置き場所を除く。"""
    from metabolomix.pipeline import recovery
    from tests.lab_msp_fixtures import write_lab_msp_config
    source = _prepare_source(tmp_path, monkeypatch)
    write_lab_msp_config(tmp_path, monkeypatch, directory_name="secret_share")
    receipt = service.plan_pipeline(source["root"], {"library_mode": "msp_only"})
    status = recovery.read_status(Path(receipt["pipeline_path"]), include_details=True)
    text = json.dumps(status, ensure_ascii=False)
    assert "secret_share" not in text
    inputs = status["record"]["inputs"]
    assert inputs["msp"]["file"] == "lab_neg.msp"
    assert inputs["method"]["overrides"]["Msp file path"] == "lab_neg.msp"
    # データ側の記録そのものは絶対パスのまま（resume の再検査に要る）。
    from metabolomix.pipeline import store
    saved = store.load_run(recovery.normalize_pipeline_root(Path(receipt["pipeline_path"])))
    assert saved["inputs"]["msp"]["path"].endswith("lab_neg.msp")
    assert "secret_share" in saved["inputs"]["msp"]["path"]
