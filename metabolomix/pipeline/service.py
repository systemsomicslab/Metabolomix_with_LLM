"""生データフォルダpipelineの公開サービス層（`pipeline_*` ツールの実体。plan 2026-11-09）。

pipeline は上流だけを回す: 極性の単位ごとに raw を配置し（`prepare_input`）、MS-DIAL
Console を回し（`upstream`）、出来た mzTab-M と解析フォルダ（`.arf`）を確かめ
（`validate_outputs`）、最後に受領（`receipt`）を書く。続き（データセットを開く・前処理・
比較）は対話の側で `load_dataset(pipeline_path=…, polarity=…)` から。

このモジュールが解く問題は2つ。

1. **受付**（`plan_pipeline` / `start_pipeline`）: `resolve_request` → 単位ごとに
   `inspect_inputs`（読取専用の入力計画）と `sample_manifest` の事前検査 →
   `find_or_create_run` → （`start_pipeline`だけ）workerの起動、という順で進める。
   `plan_pipeline`は検査・保存のみで起動しない。壊れたシート等の既知の不正入力は、
   Console起動前に検出し、runは作るがworkerは起動せず`needs_input`のまま保存する
   （spec §9.2「既知の不正入力はConsole起動前に拒否する」）。
2. **工程handler一式**（`build_handlers`）: `engine` が要求するstage契約
   （`context: dict -> StageResult`）へ薄く配線する。ワーカー（`pipeline.worker`）だけが
   これを呼ぶ——MCP接続やグローバルsessionへは一切触れない（`metabolomix.core.session_state`
   / `metabolomix.core.mcp_core` / `metabolomix.tools.*` をimportしない）。

**Console用workerを二重起動しない**: `upstream` handlerは共通監視経路
`metabolomix.console.execution.supervise`を直接呼ぶ——pipeline用workerが既に
「起動して見張る」役を担っているため。単位は要求の順に 1 つずつ回す（1 台の PC で
Console を 2 本並べない）。

**launchはsys.executableで、cwdは常にこのcheckout自身**（`_REPO_ROOT`）。
"""
from __future__ import annotations

import copy
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.analysis.sample_manifest import parse_manifest
from metabolomix.console import execution as console_execution
from metabolomix.console import job_manager
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.core.process_control import launch_detached
from metabolomix.handoff.schema import SCHEMA_VERSION, AnalysisJob
from metabolomix.dataset import alignment_files
from metabolomix.mztab import evidence as mztab_evidence
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import engine
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import recovery
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import store

__all__ = [
    "build_handlers",
    "launch_pipeline_worker",
    "plan_pipeline",
    "resume_pipeline",
    "start_pipeline",
]

#: このcheckoutのルート（`metabolomix/pipeline/`の2階層上）。launchのcwdを
#: 呼び出し元のcwdや別checkoutへ暗黙に依存させないため、ここで固定する。
_REPO_ROOT = Path(__file__).resolve().parents[2]

#: 単位の根（`<pipeline_root>/<極性>`）の下の書込み先。
_CONSOLE_RUN_SUBDIR = "console"
_INPUT_SUBDIR = "input"  # metabolomix.pipeline.inputs._INPUT_SUBDIRと同じ値
#: `pipeline_plan` / `pipeline_run` の戻り値に出す受付の警告（`request.resolve_request` が出す）。
_INTAKE_WARNING_CODES = frozenset({"REQUEST_FILE_DOWNSTREAM_IGNORED",
                                   "REQUEST_FILE_UNIT_KEYS_IGNORED"})

#: 起動受理を待つ上限（秒）。「短い上限」——超えても失敗と決め付けず、
#: pipeline_pathと（未確認の）launch状態を返す。
_HANDSHAKE_TIMEOUT_S = 3.0
_HANDSHAKE_POLL_S = 0.05


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# 受付: plan_pipeline / start_pipeline
# ---------------------------------------------------------------------------

def _resolve_exe_path() -> Path:
    from metabolomix.console import runner as console_runner
    try:
        exe = console_runner.get_exe_path()
    except EnvironmentError as exc:
        # 未設定は MSDIAL_EXE_NOT_FOUND、設定ファイルが読めなければ CONFIG_INVALID。
        # どちらも利用者が直す手掛かり（設定ファイルの場所・キー・行）を details に持つ。
        code = getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND")
        details = dict(getattr(exc, "details", {}) or {})
        if code == "MSDIAL_EXE_NOT_FOUND":
            # 設定ファイルの編集は MCP クライアントからはできない。人間の作業であることと
            # コピペできる手順を機械可読で添える。
            details = {**console_runner.msdial_exe_setup_help(), **details}
        raise DomainError(code, str(exc), details) from exc
    return Path(exe)


def _plan_fingerprint(plan: dict) -> str:
    """単位の入力計画の内容hash（受付冪等性の根拠。spec §9.3）。

    採用する形式・raw stat・メソッド/LBM/実行体のhash・極性だけを対象にする
    ——`overrides`（絶対パスを含みうる）や`unverified`（判断の説明文）は
    「同一入力」の判定に無関係なので混ぜない。
    """
    content = {
        "raw_stat": plan["raw_stat"],
        "selected_format": plan["selected_format"],
        "method_sha256": plan["method"]["sha256"],
        "lbm_sha256": plan["lbm"]["sha256"],
        "exe_sha256": plan["exe"]["sha256"],
        "polarity": plan["polarity"]["value"],
    }
    # 研究室 MSP は msp_only / msp_then_lbm のときだけ入れる（auto の指紋は変えない）。
    msp_sha256 = (plan.get("msp") or {}).get("sha256")
    if msp_sha256:
        content["msp_sha256"] = msp_sha256
    return canonical_hash(content)


def _precheck_manifest(raw_root: Path, unit_req: dict, plan: dict) -> DomainError | None:
    """既知の不正入力（壊れたシート等）をConsole起動前に検出する（spec §9.2）。

    `expected_sources`は`plan["entries"]`のトップレベル名（primaryのみ）——
    ディレクトリ形式raw（`.d`等）でも1測定単位=1entryになるのはここだけ。
    """
    manifest_path = inputs_mod.resolve_manifest_path(raw_root, unit_req)
    if manifest_path is None:
        return None
    if not manifest_path.is_file():
        return DomainError(
            "SAMPLE_MANIFEST_NOT_FOUND",
            f"sample_manifestが指すファイルが見つかりません: {manifest_path}",
            {"path": str(manifest_path)})
    expected_sources = [e["name"] for e in plan["entries"] if e.get("role") == "primary"]
    try:
        parse_manifest(manifest_path, source_root=raw_root, expected_sources=expected_sources)
    except DomainError as exc:
        return exc
    return None


def _manifest_rows_for_planning(raw_root: Path, unit_req: dict) -> list[dict] | None:
    """計画（アライメントの参照ファイルの選択）に使う実験情報シートの行を、先に読む。

    シートが無い・読めない・予定の raw と合わない場合は None（ファイル名順へ戻す）。
    不正なシートの報告は`_precheck_manifest`の役目で、ここでは止めない。
    """
    manifest_path = inputs_mod.resolve_manifest_path(raw_root, unit_req)
    if manifest_path is None or not manifest_path.is_file():
        return None
    try:
        expected_sources = inputs_mod.list_primary_names(raw_root, unit_req)
        return parse_manifest(manifest_path, source_root=raw_root,
                              expected_sources=expected_sources)
    except DomainError:
        return None


def _mark_needs_input_without_launch(pipeline_path: Path, exc: DomainError, *, stage_id: str) -> None:
    """事前検査で見つかった不正入力を、workerを起動せずneeds_inputとして保存する。

    既に処理済み（find_or_create_runの再利用等でplanned以外）のrunは書き換えない。
    """
    record = store.load_run(pipeline_path)
    if record["status"] != "planned":
        return
    record = copy.deepcopy(record)
    record["status"] = "needs_input"
    record["needs_input"] = {
        "code": exc.code, "stage_id": stage_id, "message": exc.message,
        "details": dict(exc.details),
    }
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])


def _resolved_settings(inputs: dict) -> dict:
    """`pipeline_plan`の戻り値へ載せる、単位の「解決した実行条件」だけを抜き出す。

    載せるのは method / lbm / polarity / library_mode / msp / alignment_reference
    （msp_only / msp_then_lbm では実効の足切り dot_product_cutoffs も）だけ。研究室の
    ライブラリは置き場所を載せずファイル名・sha256・出どころ（spec 2026-10-07 §3.3）。
    `raw_stat`・`entries`・`overrides`のような大きな中間データは載せない。
    """
    method = inputs.get("method") or {}
    lbm = inputs.get("lbm") or {}
    polarity = inputs.get("polarity") or {}
    msp = inputs.get("msp") or {}
    library_mode = inputs.get("library_mode", "auto")
    resolved_msp = {"file": msp.get("file"), "sha256": msp.get("sha256"),
                    "source": msp.get("source", "not_used")}
    if library_mode == "msp_then_lbm":
        resolved_lbm = {"file": lbm.get("file"), "sha256": lbm.get("sha256"),
                        "source": lbm.get("source")}
        resolved_msp["priority"] = msp.get("priority")
    elif inputs_mod.is_lab_lbm(inputs):
        # auto でも [library] lbm から解いた LBM は研究室の資産。path の中身をファイル名にする。
        resolved_lbm = {"path": Path(lbm["path"]).name if lbm.get("path") else None,
                        "sha256": lbm.get("sha256")}
    else:
        resolved_lbm = {"path": lbm.get("path"), "sha256": lbm.get("sha256")}
    resolved = {
        "method": {"source_path": method.get("source_path"), "sha256": method.get("sha256")},
        "lbm": resolved_lbm,
        "polarity": {"value": polarity.get("value"), "source": polarity.get("source")},
        "library_mode": library_mode,
        "msp": resolved_msp,
        "alignment_reference": inputs.get("alignment_reference"),
    }
    if library_mode in ("msp_only", "msp_then_lbm"):
        resolved["dot_product_cutoffs"] = inputs.get("dot_product_cutoffs")
    return resolved


def _dispatch_summary(pipeline_path: Path, *, launched: bool, launch: dict | None = None,
                      include_resolved: bool = False) -> dict:
    """コンパクトな発送の要約（CLAUDE.md「戻り値を肥大させない」）。

    status・pipeline_id/path・単位（極性）・needs_input（あれば）・受付の警告（あれば）・
    launch状態だけを返す。`include_resolved=True`（`plan_pipeline`だけが渡す）のときに
    限り、単位ごとの解決済み条件を`resolved.units.<極性>`として足す。
    """
    record = store.load_run(pipeline_path)
    units = store.run_units(record)
    summary: dict = {
        "status": record["status"],
        "pipeline_id": record["identity"]["pipeline_id"],
        "pipeline_path": str(pipeline_path),
        "units": list(units),
        "launched": launched,
    }
    if include_resolved:
        summary["resolved"] = {"units": {p: _resolved_settings(u.get("inputs") or {})
                                         for p, u in units.items()}}
    if record.get("needs_input"):
        summary["needs_input"] = record["needs_input"]
    # 受付の警告（要求ファイルから読み飛ばしたもの）。読み飛ばした値こそ受付の時点で見せる。
    intake = [w for w in record.get("warnings") or []
              if w.get("code") in _INTAKE_WARNING_CODES]
    if intake:
        summary["warnings"] = intake
    if launch is not None:
        summary["launch"] = launch
    return summary


def _recheck_manifests(pipeline_path: Path) -> tuple[str, DomainError] | None:
    """保存済みの要求と入力計画で、単位ごとに `_precheck_manifest` をやり直す。

    最初に見つかった不正を `(stage_id, error)` で返す。どれも通れば None。
    """
    record = store.load_run(pipeline_path)
    request = json.loads(
        (pipeline_path / record["request"]["saved_path"]).read_text(encoding="utf-8"))
    source_root = Path(record["identity"]["source_root"])
    for index, unit in enumerate(request["units"]):
        polarity = unit["polarity"]
        unit_state = record["units"][polarity]
        unit_root = pipeline_path / (unit_state.get("unit_root") or polarity)
        plan = store.absolutize_unit_inputs(unit_state["inputs"], unit_root)
        error = _precheck_manifest(
            request_mod.unit_source_root(source_root, unit),
            request_mod.unit_request(request, index, source_root), plan)
        if error is not None:
            return store.unit_stage_id("prepare_input", polarity), error
    return None


def _prepare_run(dataset_root: Path, request: dict | None,
                 request_id: str | None) -> tuple[Path, DomainError | None]:
    """resolve_request→単位ごとの inspect_inputs と事前検査→find_or_create_run。

    戻り値は`(pipeline_path, manifest_error)`——`manifest_error`が非Noneなら
    workerを起動してはいけない。
    """
    source_root = Path(dataset_root).expanduser()
    resolved = request_mod.resolve_request(source_root, request)
    exe_path = _resolve_exe_path()

    plans: dict[str, dict] = {}
    fingerprints: dict[str, str] = {}
    manifest_error: tuple[str, DomainError] | None = None
    for index in range(len(resolved["units"])):
        unit = resolved["units"][index]
        raw_root = request_mod.unit_source_root(source_root, unit)
        unit_req = request_mod.unit_request(resolved, index, source_root)
        # アライメントの参照ファイルを注入順の中央にするため、シートを先に解して渡す。
        plan = inputs_mod.inspect_inputs(
            raw_root, unit_req, exe_path=exe_path,
            manifest_rows=_manifest_rows_for_planning(raw_root, unit_req))
        polarity = plan["polarity"]["value"]
        if polarity in plans:
            raise DomainError(
                "PIPELINE_REQUEST_INVALID",
                f"2 つの単位が同じ極性（{polarity}）に解決されました。polarity を明示してください。",
                {"polarity": polarity})
        resolved = request_mod.pin_polarity(resolved, index, polarity)
        plans[polarity] = plan
        fingerprints[polarity] = _plan_fingerprint(plan)
        if manifest_error is None:
            error = _precheck_manifest(raw_root, unit_req, plan)
            if error is not None:
                manifest_error = (polarity, error)

    inputs = {"fingerprint": canonical_hash({"units": fingerprints}), "units": plans}
    pipeline_path = store.find_or_create_run(
        source_root, resolved, inputs, request_id=request_id)

    if manifest_error is not None:
        polarity, error = manifest_error
        _mark_needs_input_without_launch(
            pipeline_path, error, stage_id=store.unit_stage_id("prepare_input", polarity))
        return pipeline_path, error
    return pipeline_path, None


def plan_pipeline(dataset_root: Path, request: dict | None = None,
                  request_id: str | None = None) -> dict:
    """入力検査・不足情報・固定要求の保存のみを行う。Consoleは起動しない。

    戻り値には単位ごとの解決済み条件を`resolved.units.<極性>`として載せる
    ——起動前に解決結果を確認する入口として案内しているため（`_resolved_settings`）。
    """
    pipeline_path, _manifest_error = _prepare_run(dataset_root, request, request_id)
    return _dispatch_summary(pipeline_path, launched=False, include_resolved=True)


def start_pipeline(dataset_root: Path, request: dict | None = None,
                   request_id: str | None = None) -> dict:
    """計画と起動を一括実行する。

    起動受理は、workerのidentity保存と起動handshakeを確認した後に返す
    （`_await_launch_handshake`）。応答待ちには短い上限を設け、上限時は起動失敗と
    決め付けず再起動もしない。

    `find_or_create_run`が既存runを再利用した場合は起動しない。新規に作られたrunだけが
    `status="planned"`のまま返るので、それだけを起動条件とする。既存runを本当に
    進めたい利用者は`pipeline_resume`を使う。
    """
    pipeline_path, manifest_error = _prepare_run(dataset_root, request, request_id)
    if manifest_error is not None:
        return _dispatch_summary(pipeline_path, launched=False)

    baseline_record = store.load_run(pipeline_path)
    if baseline_record["status"] != "planned":
        return _dispatch_summary(pipeline_path, launched=False)

    launch_info = _launch_and_await(pipeline_path, baseline_record=baseline_record)
    return _dispatch_summary(pipeline_path, launched=True, launch=launch_info)


def resume_pipeline(path: Path, updates: dict | None = None,
                    request_id: str | None = None, rerun_upstream: bool = False) -> dict:
    """試料の対応表の訂正・停止した工程からの再開。`prepare_resume`の後、結果が
    `planned`のときだけworkerを起動する（no-op resumeでは起動しない）。

    `path`は`pipeline-run.json`そのものでもその親（pipeline_root）でもよい。正規化した
    `pipeline_path`を戻り値・`launch_pipeline_worker`・handshakeの全経路で使い回す。
    """
    pipeline_path = recovery.normalize_pipeline_root(path)
    prior_status = store.load_run(pipeline_path)["status"]
    result = recovery.prepare_resume(
        pipeline_path, updates=updates, request_id=request_id, rerun_upstream=rerun_upstream)

    summary = dict(result)
    summary["pipeline_path"] = str(pipeline_path)
    if result["status"] == "planned" and prior_status == "needs_input":
        # 受付の事前検査で止まった run は、起動する前に同じ検査をやり直す（シートが
        # 直っていないまま worker を起動しない）。
        found = _recheck_manifests(pipeline_path)
        if found is not None:
            stage_id, error = found
            _mark_needs_input_without_launch(pipeline_path, error, stage_id=stage_id)
            record = store.load_run(pipeline_path)
            summary["status"] = record["status"]
            summary["needs_input"] = record.get("needs_input")
            summary["launch"] = {"launched": False}
            return summary
    if result["status"] == "planned":
        baseline_record = store.load_run(pipeline_path)
        launch_info = _launch_and_await(pipeline_path, baseline_record=baseline_record)
        summary["launch"] = launch_info
    else:
        summary["launch"] = {"launched": False}
    return summary


# ---------------------------------------------------------------------------
# worker起動
# ---------------------------------------------------------------------------

def launch_pipeline_worker(path: Path) -> dict:
    """`metabolomix.pipeline.worker`を切り離して起動する（待たない）。

    `sys.executable`で起動し、cwdは`_REPO_ROOT`へ固定する。本番workerは試験専用の
    環境変数・任意import名を一切受け付けない（引数は`--pipeline`のみ）。
    stdout/stderrは`launch_detached`がログファイルへ結ぶ。
    """
    pipeline_path = Path(path)
    command = [sys.executable, "-m", "metabolomix.pipeline.worker",
               "--pipeline", str(pipeline_path)]
    log_path = pipeline_path / store.CONTROL_SUBDIR / "worker-launch.log"
    return launch_detached(command, cwd=_REPO_ROOT, log_path=log_path)


def _launch_and_await(pipeline_path: Path, *, baseline_record: dict) -> dict:
    """launch直前の状態を基準にworkerを起動し、handshakeを待つ。"""
    launch_info = launch_pipeline_worker(pipeline_path)
    baseline_identity = (baseline_record.get("worker") or {}).get("identity")
    baseline_status = baseline_record["status"]
    launch_info["handshake"] = _await_launch_handshake(
        pipeline_path, baseline_identity=baseline_identity, baseline_status=baseline_status)
    return launch_info


def _await_launch_handshake(pipeline_path: Path, *, baseline_identity, baseline_status: str) -> str:
    """workerのidentity保存を短時間だけ待つ。超えても失敗と決め付けない。

    今回のlaunchで**新しい**identityが書かれた（baselineと異なる）か、statusが
    baselineから動いた場合だけを「本当にworkerが起動した」証拠として扱う。
    """
    deadline = time.monotonic() + _HANDSHAKE_TIMEOUT_S
    while time.monotonic() < deadline:
        record = store.load_run(pipeline_path)
        identity = (record.get("worker") or {}).get("identity")
        if identity is not None and identity != baseline_identity:
            return "confirmed"
        if record.get("status") != baseline_status:
            return "confirmed"
        time.sleep(_HANDSHAKE_POLL_S)
    return "not_confirmed"


# ---------------------------------------------------------------------------
# build_handlers: 工程handler一式
# ---------------------------------------------------------------------------

def build_handlers() -> dict:
    """`engine.run_engine`へ渡すhandler一式（キーはhandlerキー。stage_idではない）。"""
    return {
        "prepare_input": _handle_prepare_input,
        "upstream": _handle_upstream,
        "validate_outputs": _handle_validate_outputs,
        "receipt": _handle_receipt,
    }


def _ok(warnings: list | None = None, **extra) -> dict:
    return {"status": "succeeded", "result_refs": [], "warnings": list(warnings or []),
            "error": None, **extra}


def _pipeline_relative(path: Path, pipeline_root: Path) -> str:
    """pipeline_root 相対（記録に残すパスの規約）。外を指すなら絶対のまま。"""
    try:
        return str(Path(path).resolve().relative_to(Path(pipeline_root).resolve())).replace("\\", "/")
    except ValueError:
        return str(path)


# ---------- prepare_input ----------

def _handle_prepare_input(context: dict) -> dict:
    """単位の入力計画（`unit_state["inputs"]`＝`inspect_inputs`の出力）を
    `<pipeline_root>/<極性>/input/` へ実配置する。

    `rerun_upstream=True`がこのstageをpendingへ戻した場合は、既存の配置済みファイルと
    突き合わせて再配置する（`stage_inputs`自身の冪等性）。
    """
    unit_root = context["unit_root"]
    plan = store.absolutize_unit_inputs(context["unit_state"]["inputs"], unit_root)
    snapshot = inputs_mod.stage_inputs(plan, unit_root)
    snapshot = store.relativize_unit_inputs(snapshot, unit_root)
    warnings = [{"code": "INPUT_UNVERIFIED", "message": note}
                for note in plan.get("unverified", [])]
    return _ok(warnings, unit_updates={"inputs": snapshot})


# ---------- upstream ----------

def _handle_upstream(context: dict) -> dict:
    """`execution.supervise`を直接呼ぶ（Console用workerを二重起動しない）。

    `console_job_path`は`analysis-job.json`そのものへの、pipeline_root 相対のパスとして
    記録する（`run_dir = job_path.parent` の規約）。

    **run_dirはattemptごとに分ける**（`<極性>/console/attempt-NNNN/`。spec §9.3
    「Consoleの再試行だけは`rerun_upstream=true`の明示を必要とし、新しいjob/attemptを
    作る」）。固定パスにすると、再実行が前回の終了証跡・ログ・実行前スナップショットを
    上書きし、timeoutで打ち切った事実が後から証明できなくなる。stageの`attempt`は
    `engine.mark_stage_running`が単調に増やすので、これが「何回目のConsole実行か」になる。
    """
    pipeline_root = context["pipeline_root"]
    unit_root = context["unit_root"]
    request = context["unit_request"]
    inputs_snapshot = store.absolutize_unit_inputs(context["unit_state"]["inputs"], unit_root)
    attempt = max(1, int(context.get("attempt") or 1))
    run_dir = unit_root / _CONSOLE_RUN_SUBDIR / f"attempt-{attempt:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    job_path = run_dir / job_manager.JOB_FILENAME
    dataset_root = unit_root / _INPUT_SUBDIR

    if not job_path.is_file():
        method_rel = inputs_snapshot["method"].get("effective_relative_path")
        method_abs = ((unit_root / method_rel) if method_rel
                      else Path(inputs_snapshot["method"]["source_path"]))
        now = _now_iso()
        job = AnalysisJob(
            schema=SCHEMA_VERSION, job_id=f"pipeline_{context['pipeline_id']}_{context['unit']}",
            status="planned", created_at=now, updated_at=now,
            dataset_root=str(dataset_root),
            input_count=job_manager.count_raw_inputs(dataset_root),
            software_name="MS-DIAL", software_version="", execution_mode="console",
            method_file=str(method_abs), omics="lipidomics",
            polarity=inputs_snapshot["polarity"]["value"], measure=request.get("measure"),
            run_dir=str(run_dir), save_project=request["save_project"],
            timeout_s=request["timeout_s"],
        )
        job.save(job_path)
        raw_inventory = [str(p.resolve()) for p in job_manager.list_raw_inputs(dataset_root)]
        console_execution.write_supervision_inputs(run_dir, {
            "raw_inventory": raw_inventory,
            "method_sha256": inputs_snapshot["method"].get("effective_sha256"),
            "exe_path": inputs_snapshot["exe"]["path"],
            "exe_sha256": inputs_snapshot["exe"].get("sha256"),
        })

    supervised = console_execution.supervise(
        job_path, cancel_path=engine.cancel_request_path(pipeline_root))
    job = job_manager.load_job(job_path)

    upstream_update = {
        "console_job_path": _pipeline_relative(job_path, pipeline_root),
        "execution_id": supervised.get("execution_id"),
        "verification": {"status": job.status, "termination": supervised.get("termination"),
                         "exit_code": supervised.get("exit_code")},
    }
    warnings = [{"code": "CONSOLE_WARNING", "message": w} for w in job.warnings[:5]]

    if job.status == "completed":
        return _ok(warnings, unit_updates={"upstream": upstream_update})

    message = job.error or (
        f"MS-DIAL Consoleの実行が完了しませんでした"
        f"（termination={supervised.get('termination')}）。")
    return {"status": "failed", "result_refs": [], "warnings": warnings,
            "error": {"code": "MSDIAL_EXECUTION_FAILED", "message": message,
                      "details": {"termination": supervised.get("termination"),
                                  "exit_code": supervised.get("exit_code"),
                                  "unit": context["unit"]}},
            "unit_updates": {"upstream": upstream_update}}


# ---------- validate_outputs ----------

def _find_gui_project(run_dir: Path) -> Path | None:
    candidates = sorted(run_dir.rglob("*.mdproject"))
    return candidates[0] if candidates else None


def _matching_alignment(gate, input_dir: Path) -> tuple[dict | None, dict | None]:
    """解析フォルダ（`<極性>/input/`）の `.arf` のうち、正準 mzTab-M と同じアラインメントだと
    確かめられたものがちょうど 1 本あれば `{"arf", "arf2"}`（名前）を返す。

    選び方は `dataset/opening` と同じ（`mztab.evidence.first_matching_arfs`）: まずこの
    attempt のジョブが集めた `peak_matrix_source`、どれも通らなければ解析フォルダのバッチを
    新しい順に 1 本ずつ。照合は `matching_arfs`（特徴数・m/z・試料名）。`rerun_upstream` の
    やり直しで前の attempt のバッチが残っていても、それを 2 本目に数えない。
    合わなければ (None, 警告)。
    """
    details: dict = {"analysis_dir": f"{input_dir.parent.name}/{input_dir.name}"}

    axes: list = []

    def verify(paths):
        # 照合の軸（特徴の m/z・試料名）は mzTab-M を 1 回だけ読んで使い回す（フォールバックの
        # バッチごとに呼ばれるので、毎回読むとバッチの数だけ同じ mzTab-M を解析する）。
        if not axes:
            state, _parsed = mztab_loading.read_mztab_parsed(gate.mztab_path)
            axes.extend(mztab_evidence.join_axes(state))
        return mztab_evidence.matching_arfs(paths, feature_mz=axes[0], sample_names=axes[1])

    try:
        passed, rejected, _scan, candidates = mztab_evidence.first_matching_arfs(
            gate.mztab_path, mztab_loading.job_artifact_paths(gate.job), verify,
            alignment_files.batch_arfs_newest_first(input_dir))
    except DomainError as exc:
        return None, {"code": "ANALYSIS_FOLDER_UNMATCHED",
                      "message": f"解析フォルダの .arf を mzTab-M と照合できませんでした"
                                 f"（{exc.code}）。mzTab-M（source=\"mztab\"）なら読めます。",
                      "details": {**details, "error_code": exc.code}}
    details["candidates"] = [p.name for p in candidates]
    if len(passed) == 1:
        # 同じバッチの `.arf2`（`<stem>_PeakProperties.arf` → `<stem>.arf2`。開く側と同じ規則）。
        try:
            arf2 = alignment_files.select_alignment(passed[0]).arf2.name
        except DomainError:
            arf2 = None
        return {"arf": passed[0].name, "arf2": arf2}, None
    if not candidates:
        return None, {"code": "ANALYSIS_FOLDER_UNMATCHED",
                      "message": "解析フォルダに MS-DIAL のアラインメント（.arf）がありません。"
                                 "mzTab-M（source=\"mztab\"）なら読めます。",
                      "details": details}
    reason = ("複数あります: " + ", ".join(p.name for p in passed)) if passed else "ありません"
    return None, {"code": "ANALYSIS_FOLDER_UNMATCHED",
                  "message": f"mzTab-M と同じアラインメントだと確かめられた .arf が{reason}。"
                             "mzTab-M（source=\"mztab\"）なら読めます。",
                  "details": {**details, "matching": [p.name for p in passed],
                              "rejected": rejected}}


def _handle_validate_outputs(context: dict) -> dict:
    """Console の出力を確かめる: 正準 mzTab-M（ジョブの門）・解析フォルダの `.arf`・
    `save_project` なら GUI project（`*.mdproject`）。

    mzTab-M が無い・ジョブが完了していないなら失敗（その単位は読めない）。解析フォルダの
    `.arf` が mzTab-M と合わないのは失敗にしない（警告 `ANALYSIS_FOLDER_UNMATCHED`。
    `source="mztab"` なら読めるので、run は `partial`）。
    """
    pipeline_root = context["pipeline_root"]
    unit_root = context["unit_root"]
    job_rel = (context["unit_state"].get("upstream") or {}).get("console_job_path")
    if not job_rel:
        raise DomainError("PIPELINE_RUN_INVALID", "upstream の記録（console_job_path）がありません。",
                          {"unit": context["unit"]})
    job_path = pipeline_root / job_rel
    gate = mztab_loading.job_gate(job_path, allow_incomplete=False)

    warnings: list[dict] = []
    result_refs: list[dict] = []
    alignment, warning = _matching_alignment(gate, unit_root / _INPUT_SUBDIR)
    if warning is not None:
        warnings.append(warning)

    gui_project = None
    if context["unit_request"].get("save_project"):
        project_path = _find_gui_project(job_path.parent)
        if project_path is None:
            warnings.append({"code": "GUI_PROJECT_UNAVAILABLE",
                             "message": "save_project=trueですがGUIプロジェクト(.mdproject)が"
                                        "見つかりません。"})
        else:
            ref = receipt_mod.persist_result(pipeline_root, {
                "output_name": f"gui_project:{context['unit']}", "kind": "gui_project",
                "result_id": f"res_gui_project_{context['pipeline_id']}_{context['unit']}",
                "path": str(project_path),
            })
            result_refs.append(ref)
            gui_project = ref["relative_path"]

    outputs = {"mztab": _pipeline_relative(gate.mztab_path, pipeline_root),
               "alignment": alignment, "gui_project": gui_project}
    return {"status": "succeeded", "result_refs": result_refs, "warnings": warnings,
            "error": None, "unit_updates": {"outputs": outputs}}


# ---------- receipt ----------

def _handle_receipt(context: dict) -> dict:
    """受領（`results/receipt.json`）を書いて登録する。受領の警告は工程の警告にする。"""
    pipeline_root = context["pipeline_root"]
    record = store.load_run(pipeline_root)
    receipt, ref = receipt_mod.write_receipt(record, pipeline_root)
    warnings = [{"code": w["code"], "message": w["message"]}
                for w in receipt.get("warnings") or []]
    return {"status": "succeeded", "result_refs": [ref], "warnings": warnings, "error": None}
