"""生データフォルダpipelineの公開サービス層（Task18: 5件のMCPツールの実体）。

このモジュールが解く問題は2つ。

1. **受付**（`plan_pipeline` / `start_pipeline`）: `resolve_request` →
   `inspect_inputs`（読取専用の入力計画）→ `sample_manifest` の事前検査 →
   `find_or_create_run` → （`start_pipeline`だけ）workerの起動、という順で
   進める。`plan_pipeline`は検査・保存のみで起動しない。壊れたシート等の
   既知の不正入力は、Console起動前に検出し、runは作るがworkerは起動せず
   `needs_input`のまま保存する（spec §9.2「既知の不正入力はConsole起動前に
   拒否する」）。target=differential・comparisons未指定はここでは特別扱いしない
   ——`resolve_comparisons` stage（`build_handlers`）がworker内で
   `COMPARISON_REQUIRED`のneeds_inputとして検出する（「下流の群不足だけは
   workerを起動できる」）。
2. **工程handler一式**（`build_handlers`）: Task 3/4/6/7/8/9/11/12/13/17が
   実装した純関数を、`metabolomix.pipeline.engine`が要求するstage契約
   （`context: dict -> StageResult`）へ薄く配線する。ワーカー
   （`metabolomix.pipeline.worker`）だけがこれを呼ぶ——MCP接続やグローバル
   sessionへは一切触れない（`metabolomix.core.session_state` /
   `metabolomix.core.mcp_core` / `metabolomix.tools.*` をimportしない。
   `tests/test_pipeline_engine.py`のASTテストと同じ制約を本モジュールにも課す）。

**Console用workerを二重起動しない**: `upstream` handlerは
`metabolomix.console.worker.run_job`/`launch_console_worker`を一切呼ばず、
共通監視経路`metabolomix.console.execution.supervise`を直接呼ぶ——pipeline用
workerが既に「起動して見張る」役を担っているため。

**launchはsys.executableで、cwdは常にこのcheckout自身**（`_REPO_ROOT`）。
利用者のcwdや別checkoutを暗黙に使わない。起動時のstdout/stderrは
`process_control.launch_detached`がログファイルへ結ぶため、launcher（この
プロセス）がworkerのstdoutパイプを継承して`communicate()`相当が詰まる
経路は無い。
"""
from __future__ import annotations

import copy
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.analysis import dataset_export, dataset_service
from metabolomix.analysis.dataset_analysis import PreconditionError
from metabolomix.analysis.sample_manifest import parse_manifest
from metabolomix.analysis.sample_manifest import resolve_metadata as resolve_sample_metadata
from metabolomix.console import execution as console_execution
from metabolomix.console import job_manager
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.core.process_control import launch_detached
from metabolomix.handoff.schema import SCHEMA_VERSION, AnalysisJob
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import engine
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import recovery
from metabolomix.pipeline import report as report_mod
from metabolomix.console import profiles as profiles_mod
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import request_v2 as request_v2_mod
from metabolomix.pipeline import stage_plan
from metabolomix.pipeline import store
from metabolomix.plots import result_output

__all__ = [
    "build_handlers",
    "launch_pipeline_worker",
    "plan_pipeline",
    "resume_pipeline",
    "start_pipeline",
]

#: このcheckoutのルート（`metabolomix/pipeline/`の2階層上）。launchのcwdを
#: 呼び出し元のcwdや別checkoutへ暗黙に依存させないため、ここで固定する
#: （`metabolomix/console/worker.py::_repo_root`と同じ流儀）。
_REPO_ROOT = Path(__file__).resolve().parents[2]

#: request.sample_manifest省略時に探す既定ファイル名（spec §7.1/§10.1「既定名を
#: 探索。存在しなければ自動一覧生成」）。実体は`metabolomix.pipeline.inputs`にあり、
#: 受付の事前検査・`resolve_metadata` handler・再開時のシート内容照合
#: （`recovery.prepare_resume`）が同じ1つの定数と解決規則を共有する。
_DEFAULT_MANIFEST_NAME = inputs_mod.DEFAULT_MANIFEST_NAME

#: pipeline_root配下の各種書込み先。
_CONSOLE_RUN_SUBDIR = "console"
_INPUT_SUBDIR = "input"  # metabolomix.pipeline.inputs._INPUT_SUBDIRと同じ値
_RESULTS_SUBDIR = "results"

#: 起動受理を待つ上限（秒）。「短い上限」——超えても失敗と決め付けず、
#: pipeline_pathと（未確認の）launch状態を返す（brief拘束）。
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
        raise DomainError(getattr(exc, "code", "MSDIAL_EXE_NOT_FOUND"), str(exc),
                          dict(getattr(exc, "details", {}) or {})) from exc
    return Path(exe)


def _plan_fingerprint(plan: dict) -> str:
    """入力計画の内容hash（受付冪等性の根拠。spec §9.3）。

    採用する形式・raw stat・メソッド/LBM/実行体のhash・極性だけを対象にする
    ——`overrides`（絶対パスを含みうる）や`unverified`（判断の説明文）は
    「同一入力」の判定に無関係なので混ぜない。
    """
    return canonical_hash({
        "raw_stat": plan["raw_stat"],
        "selected_format": plan["selected_format"],
        "method_sha256": plan["method"]["sha256"],
        "lbm_sha256": plan["lbm"]["sha256"],
        "exe_sha256": plan["exe"]["sha256"],
        "polarity": plan["polarity"]["value"],
    })


def _resolve_manifest_path(source_root: Path, request: dict) -> Path | None:
    """`inputs.resolve_manifest_path`への薄い委譲（解決規則を二重に持たない）。"""
    return inputs_mod.resolve_manifest_path(source_root, request)


def _precheck_manifest(source_root: Path, request: dict, plan: dict) -> DomainError | None:
    """既知の不正入力（壊れたシート等）をConsole起動前に検出する（spec §9.2）。

    `expected_sources`は`plan["entries"]`のトップレベル名（primaryのみ）——
    ディレクトリ形式raw（`.d`等）でも1測定単位=1entryになるのはここだけで、
    `raw_stat`は内部ファイルへ展開済みのため代用できない。
    """
    manifest_path = _resolve_manifest_path(source_root, request)
    if manifest_path is None:
        return None
    if not manifest_path.is_file():
        return DomainError(
            "SAMPLE_MANIFEST_NOT_FOUND",
            f"sample_manifestが指すファイルが見つかりません: {manifest_path}",
            {"path": str(manifest_path)})
    expected_sources = [e["name"] for e in plan["entries"] if e.get("role") == "primary"]
    try:
        parse_manifest(manifest_path, source_root=source_root, expected_sources=expected_sources)
    except DomainError as exc:
        return exc
    return None


def _mark_needs_input_without_launch(pipeline_path: Path, exc: DomainError, *, stage_id: str) -> None:
    """事前検査で見つかった不正入力を、workerを起動せずneeds_inputとして保存する。

    既に処理済み（find_or_create_runの再利用等でplanned以外）のrunは書き換えない
    ——同一fingerprintの活動中/completed runを再利用した場合に、その状態を
    ここで壊してはいけない。
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


def _resolved_settings(record: dict) -> dict:
    """`pipeline_plan`のreceiptへ載せる「解決した実行条件」だけを抜き出す。

    `MCP_INSTRUCTIONS`のENTRY POINTは`pipeline_plan`を「method file・LBM・
    polarityの解決結果を、何も起動しないうちに確認する」入口として案内して
    いる。receiptがそれを持たないと、案内された確認は
    `pipeline_status(include_details=True)`でrecord全体を引くしかない
    ——コンパクトなreceiptという拘束と正面から衝突する。

    載せるのは3項目だけ。`raw_stat`・`entries`・`companions`・`overrides`の
    ような大きな中間データは`record["inputs"]`に残したまま、receiptへは
    出さない（CLAUDE.md「戻り値を肥大させない」）。
    """
    inputs = record.get("inputs") or {}
    method = inputs.get("method") or {}
    lbm = inputs.get("lbm") or {}
    polarity = inputs.get("polarity") or {}
    return {
        "method": {"source_path": method.get("source_path"), "sha256": method.get("sha256")},
        "lbm": {"path": lbm.get("path"), "sha256": lbm.get("sha256")},
        "polarity": {"value": polarity.get("value"), "source": polarity.get("source")},
    }


def _dispatch_receipt(pipeline_path: Path, *, launched: bool, launch: dict | None = None,
                      include_resolved: bool = False) -> dict:
    """コンパクトな発送receipt（CLAUDE.md「戻り値を肥大させない」）。

    行列・スコア・volcano点列は一切含めない。status・pipeline_id/path・
    effective_target・needs_input（あれば）・launch状態だけを返す。

    `include_resolved=True`（`plan_pipeline`だけが渡す）のときに限り、
    解決済みのmethod/LBM/polarityを`resolved`として足す——起動する側
    （`pipeline_run`）のreceiptは従来どおり最小のままにする。
    """
    record = store.load_run(pipeline_path)
    receipt: dict = {
        "status": record["status"],
        "pipeline_id": record["identity"]["pipeline_id"],
        "pipeline_path": str(pipeline_path),
        "effective_target": record["request"].get("effective_target"),
        "launched": launched,
    }
    if include_resolved:
        receipt["resolved"] = _resolved_settings(record)
    if record.get("needs_input"):
        receipt["needs_input"] = record["needs_input"]
    if launch is not None:
        receipt["launch"] = launch
    return receipt


def _profile_arguments(source_root: Path, request: dict | None) -> dict:
    """v2要求のときだけ、`profile_file`を読んで`resolve_request`へ渡す引数を作る。

    profileの**読み込み**はここ（受付層）の責務で、`resolve_request`は検証済みの
    dictしか受け取らない（request.pyのdocstringの契約）。v1要求では何も返さない
    ——v1にprofileの概念は無く、空のprofileを渡すと「未指定」と「指定したが空」が
    区別できなくなる。

    `execution_purpose="routine"`のときは証明書の`routine_overrides`も渡す
    ——省略すると fail-closed（何も上書きできない）になり、profileが明示的に
    許した範囲の上書きまで拒否されてしまう。
    """
    explicit = request if isinstance(request, dict) else {}
    schema = (explicit.get("schema")
              if "schema" in explicit
              else request_mod._peek_schema_declared_by_request_file(source_root))
    if schema != stage_plan.REQUEST_SCHEMA_V2:
        return {}

    from_file = request_v2_mod.read_request_file(source_root)
    profile_file = explicit.get("profile_file", from_file.get("profile_file"))
    if not isinstance(profile_file, str) or not profile_file:
        raise DomainError(
            "PIPELINE_REQUEST_INVALID",
            "v2要求にはprofile_fileが必要です（検証済みprofileへのパス）。",
            {"schema": schema})

    profile_path = Path(profile_file).expanduser()
    if not profile_path.is_absolute():
        profile_path = source_root / profile_path
    purpose = explicit.get("execution_purpose", from_file.get("execution_purpose", "routine"))
    if not isinstance(purpose, str) or purpose not in {"routine", "validation"}:
        raise DomainError("PIPELINE_REQUEST_INVALID", "execution_purposeが不正です。",
                          {"execution_purpose": purpose})
    profile = profiles_mod.load_profile(profile_path, purpose)
    routine_overrides = None
    if purpose == "routine":
        routine_overrides = profiles_mod.certificate_routine_overrides(
            profile, profile_path)
    return {"profile": profile, "routine_overrides": routine_overrides}


def _prepare_run(dataset_root: Path, request: dict | None,
                 request_id: str | None) -> tuple[Path, DomainError | None]:
    """resolve_request→inspect_inputs→manifest事前検査→find_or_create_run。

    `plan_pipeline`/`start_pipeline`の共通前半。戻り値は
    `(pipeline_path, manifest_error)`——`manifest_error`が非Noneなら
    workerを起動してはいけない。
    """
    source_root = Path(dataset_root).expanduser()
    request_resolved = request_mod.resolve_request(
        source_root, request, **_profile_arguments(source_root, request))
    if stage_plan.is_v2_request(request_resolved):
        # v1受付へ落とすとprofile以外のmethod/LBM/実行体を採用してしまう。
        # v2はprofileが唯一の情報源なので、計画を作る経路そのものを分ける。
        profile = _profile_arguments(source_root, request).get("profile")
        if not profile:
            raise DomainError(
                "PIPELINE_REQUEST_INVALID",
                "v2要求にはprofileが必要です（profile_fileが解決できませんでした）。",
                {"schema": request_resolved["schema"]})
        plan = inputs_mod.plan_from_profile(source_root, request_resolved, profile)
    else:
        exe_path = _resolve_exe_path()
        plan = inputs_mod.inspect_inputs(source_root, request_resolved, exe_path=exe_path)
    plan["fingerprint"] = _plan_fingerprint(plan)

    manifest_error = _precheck_manifest(source_root, request_resolved, plan)
    pipeline_path = store.find_or_create_run(
        source_root, request_resolved, plan, request_id=request_id)

    if manifest_error is not None:
        _mark_needs_input_without_launch(pipeline_path, manifest_error, stage_id="resolve_metadata")
    return pipeline_path, manifest_error


def plan_pipeline(dataset_root: Path, request: dict | None = None,
                  request_id: str | None = None) -> dict:
    """入力検査・不足情報・固定要求の保存のみを行う。Consoleは起動しない。

    receiptには解決済みのmethod/LBM/polarityを`resolved`として載せる
    ——`MCP_INSTRUCTIONS`がこのツールを「起動前に解決結果を確認する」入口
    として案内しているため（`_resolved_settings`）。
    """
    pipeline_path, _manifest_error = _prepare_run(dataset_root, request, request_id)
    return _dispatch_receipt(pipeline_path, launched=False, include_resolved=True)


def start_pipeline(dataset_root: Path, request: dict | None = None,
                   request_id: str | None = None) -> dict:
    """計画と起動を一括実行する。フォルダだけの通常入口。

    起動受理は、workerのidentity保存と起動handshakeを確認した後に返す
    （`_await_launch_handshake`）。応答待ちには短い上限を設け、上限時は
    起動失敗と決め付けず再起動もしない——`pipeline_path`とlaunch状態を返す。

    `find_or_create_run`が既存run（活動中／completed／失敗・取消・部分完了
    済みの終端状態）を再利用した場合は起動しない（レビュー指摘1）。新規に
    作られたrunだけが`status="planned"`のまま返るので、それだけを起動条件と
    する——さもないと、resendが活動中runの上で`PIPELINE_ALREADY_RUNNING`を
    抱えたworkerを1個ずつ増やしたり、completed runの`_ALWAYS_RECONSTRUCT_
    STAGE_IDS`（load_dataset/resolve_metadata/preprocess/pca）を無条件に
    再計算したりする。既存runを本当に進めたい利用者は`pipeline_resume`を使う。
    """
    pipeline_path, manifest_error = _prepare_run(dataset_root, request, request_id)
    if manifest_error is not None:
        return _dispatch_receipt(pipeline_path, launched=False)

    baseline_record = store.load_run(pipeline_path)
    if baseline_record["status"] != "planned":
        # 新規作成直後のrunは常にplanned。それ以外なら`find_or_create_run`が
        # 既存runを再利用したということ——ここでは起動せず、見つかった状態を
        # そのまま発送receiptで伝える。
        return _dispatch_receipt(pipeline_path, launched=False)

    launch_info = _launch_and_await(pipeline_path, baseline_record=baseline_record)
    return _dispatch_receipt(pipeline_path, launched=True, launch=launch_info)


def resume_pipeline(path: Path, updates: dict | None = None,
                    request_id: str | None = None, rerun_upstream: bool = False) -> dict:
    """入力訂正・下流revision作成・停止工程からの再開。`prepare_resume`の後、
    結果が`planned`のときだけworkerを起動する（no-op resumeでは起動しない）。

    `path`は`pipeline_status`/`pipeline_cancel`と同じく、`pipeline-run.json`
    そのものでもその親（pipeline_root）でもよい（レビュー指摘3）。
    `recovery.prepare_resume`は内部でも正規化するが、ここで先に正規化した
    `pipeline_path`をreceipt・`launch_pipeline_worker`・handshakeの全経路で
    使い回す——正規化前の生パスをどこかに残すと、`--pipeline`引数やlogの
    置き場所が二重に`pipeline-run.json`を連結した不正なパスになる。
    """
    pipeline_path = recovery.normalize_pipeline_root(path)
    result = recovery.prepare_resume(
        pipeline_path, updates=updates, request_id=request_id, rerun_upstream=rerun_upstream)

    receipt = dict(result)
    receipt["pipeline_path"] = str(pipeline_path)
    if result["status"] == "planned":
        # handshakeの基準（レビュー指摘2）は「launch直前」の状態でなければ
        # ならない——`prepare_resume`の戻り値はworker識別子を持たないため、
        # ここで改めて読み直す（前回workerの残骸identityが既に残っている
        # ケースを正しく基準にするため）。
        baseline_record = store.load_run(pipeline_path)
        launch_info = _launch_and_await(pipeline_path, baseline_record=baseline_record)
        receipt["launch"] = launch_info
    else:
        receipt["launch"] = {"launched": False}
    return receipt


# ---------------------------------------------------------------------------
# worker起動
# ---------------------------------------------------------------------------

def launch_pipeline_worker(path: Path) -> dict:
    """`metabolomix.pipeline.worker`を切り離して起動する（待たない）。

    `sys.executable`で起動し、cwdは`_REPO_ROOT`（このcheckout自身）へ固定する
    ——呼び出し元のcwdや別checkoutを暗黙に使わない。本番workerは試験専用の
    環境変数・任意import名を一切受け付けない（`metabolomix/pipeline/worker.py`の
    引数は`--pipeline`のみ）。

    stdout/stderrは`launch_detached`がログファイルへ結ぶ——launcher（この
    プロセス）がworkerのstdoutパイプを継承して詰まる経路（bInheritHandles=TRUE
    かつハンドル指定無し）を作らない。
    """
    pipeline_path = Path(path)
    command = [sys.executable, "-m", "metabolomix.pipeline.worker",
              "--pipeline", str(pipeline_path)]
    log_path = pipeline_path / store.CONTROL_SUBDIR / "worker-launch.log"
    return launch_detached(command, cwd=_REPO_ROOT, log_path=log_path)


def _launch_and_await(pipeline_path: Path, *, baseline_record: dict) -> dict:
    """launch直前の状態を基準にworkerを起動し、handshakeを待つ。

    呼び出し前提: `baseline_record`は`launch_pipeline_worker`を呼ぶ**直前**に
    読んだrecordであること（`start_pipeline`/`resume_pipeline`とも、launch
    条件を判定するために既にrecordを読んでいるので、それをそのまま渡す）。
    """
    launch_info = launch_pipeline_worker(pipeline_path)
    baseline_identity = (baseline_record.get("worker") or {}).get("identity")
    baseline_status = baseline_record["status"]
    launch_info["handshake"] = _await_launch_handshake(
        pipeline_path, baseline_identity=baseline_identity, baseline_status=baseline_status)
    return launch_info


def _await_launch_handshake(pipeline_path: Path, *, baseline_identity, baseline_status: str) -> str:
    """workerのidentity保存を短時間だけ待つ。超えても失敗と決め付けない。

    `baseline_identity`/`baseline_status`はlaunch直前の状態（レビュー指摘2）。
    resumeでは前回workerの残骸identityが既に記録されていることがあり、
    `identity is not None`だけを見ると「何も起きていない」のに即
    `"confirmed"`を返してしまう——今回のlaunchで**新しい**identityが書かれた
    (baselineと異なる)か、statusがbaselineから動いた場合だけを「本当に
    workerが起動した」証拠として扱う。
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

#: `PreconditionError.kind` → StageResultのerror code。どちらも「入力を
#: 与え直せば先へ進める」種類の停止で、needs_inputとして扱う
#: （`missing_state`は前提状態の不足、`bad_request`はその入力では成立しない
#: 指定——`min_detection_rate>0`を検出状態の無いmzTabへ指定する、等）。
_PRECONDITION_CODES = {
    "missing_state": "ANALYSIS_PRECONDITION_MISSING",
    "bad_request": "ANALYSIS_PRECONDITION_INVALID",
}


def _as_needs_input(handler):
    """handlerの`PreconditionError`をneeds_inputのStageResultへ変換して包む。

    `PreconditionError`は`DomainError`ではないので、素通りさせると
    `engine._invoke_handler`の汎用`except Exception`分岐に落ち、
    `error.code`がPythonのクラス名（`"PreconditionError"`）のまま`failed`に
    なる——実際には入力を直せば再開できる停止なのに、機械可読なcodeを持たない
    「回復不能な失敗」に見える。R18が`pca`だけに個別に施した回避を、handler
    境界の共通の規則へ引き上げる（新しいhandlerが同じ穴を開け直さないように、
    包む場所は`build_handlers`の1か所だけにする）。

    `_handle_pca`のように自分でPreconditionErrorを捕らえて固有のcodeを返す
    handlerはそのまま——ここへは到達しない。
    """
    def wrapped(context: dict) -> dict:
        try:
            return handler(context)
        except PreconditionError as exc:
            details = dict(exc.details or {})
            if exc.state:
                details.setdefault("state", exc.state)
            return {"status": "needs_input", "result_refs": [], "warnings": [],
                    "error": {
                        "code": _PRECONDITION_CODES.get(
                            exc.kind, "ANALYSIS_PRECONDITION_INVALID"),
                        "message": exc.message,
                        "details": details}}
    wrapped.__name__ = getattr(handler, "__name__", "handler")
    wrapped.__doc__ = handler.__doc__
    return wrapped


def _handle_prepare_inputs_v2(context: dict) -> dict:
    """v2の`prepare_inputs`: profile/実行環境の固定（＋rawの実配置）。

    v2の入力を決めるのは**profile**（method・依存・実行体・rawの選択と hash）で、
    v1のようにフォルダを読んでmethodとLBMを推定しない——v1のLBM必須アクセスは
    profile adapterの`dependency_keys`へ移した（plan Task 12）。

    rawの実配置は、受付層がv1形式の配置計画（`inputs.inspect_inputs`の出力）を
    作っていたときだけ行う。v2の受付が profile 由来の配置計画を作る経路は
    実Console接続（Task 15）で入れる——それまでは、profile plan の `raw`
    （実測hash付き）が「どのrawを固定したか」の記録になる。
    """
    from metabolomix.pipeline import metabolomics_handlers

    outcome = {"status": "succeeded", "result_refs": [], "warnings": [],
               "error": None}
    if "raw_stat" in (context.get("inputs") or {}):
        outcome = _handle_prepare_input(context)
        if outcome["status"] != "succeeded":
            return outcome
    return metabolomics_handlers.snapshot_profile_outcome(context, outcome)


def _handle_resolve_metadata_v2(context: dict) -> dict:
    """v2の`resolve_metadata`: v1と同じ解決に、試料対応表の成果物固定を足す。"""
    from metabolomix.pipeline import metabolomics_handlers

    outcome = _handle_resolve_metadata(context)
    if outcome["status"] != "succeeded":
        return outcome
    return metabolomics_handlers.sample_manifest_outcome(context, outcome)


def build_handlers() -> dict:
    """`metabolomix.pipeline.engine.run_engine`へ渡すhandler一式を組み立てる。

    キーはhandlerキー（stage_idではない。`differential`/`export`は
    `comparison_id`ごとの複数stageで1つのhandlerキーを共有する）。

    全handlerを`_as_needs_input`で包む——数値層が投げる`PreconditionError`
    （`DomainError`ではない）を、engineの汎用例外分岐へ落とさずneeds_inputへ
    変換するため。
    """
    from metabolomix.pipeline import metabolomics_handlers

    v2 = metabolomics_handlers.build_handlers()
    handlers = {
        # v1 の stage 名
        "prepare_input": _handle_prepare_input,
        "upstream": _handle_upstream,
        "pca": _handle_pca,
        "resolve_comparisons": _handle_resolve_comparisons,
        "differential": _handle_differential,
        # v1/v2 で同じ計算をする工程（v2 の stage 名も同じ）
        "validate_outputs": _handle_validate_outputs,
        "load_dataset": _handle_load_dataset,
        "resolve_metadata": metabolomics_handlers.by_schema(
            _handle_resolve_metadata, _handle_resolve_metadata_v2),
        # v2 だけの stage 名（v1 の計画には現れない）
        "prepare_inputs": _handle_prepare_inputs_v2,
        "execute_console": _handle_upstream,
        "load_assay_evidence": v2["load_assay_evidence"],
        "resolve_feature_bindings": v2["resolve_feature_bindings"],
        "qc_raw": v2["qc_raw"],
        "qc_processed": v2["qc_processed"],
        "statistics": v2["statistics"],
        # 名前は同じだが計算が違う工程。要求の schema で選ぶ——後勝ちで1つに
        # すると、片方の pipeline がもう片方の計算を黙って走らせる。
        "preprocess": metabolomics_handlers.by_schema(
            _handle_preprocess, v2["preprocess"]),
        "export": metabolomics_handlers.by_schema(
            _handle_export, v2["export"]),
        "report": metabolomics_handlers.by_schema(
            _handle_report, v2["report"]),
    }
    return {key: _as_needs_input(handler) for key, handler in handlers.items()}


# ---------- prepare_input（Task13） ----------

def _handle_prepare_input(context: dict) -> dict:
    """入力計画（`context["inputs"]`＝`inspect_inputs`の出力）を実配置する。

    `_TRUST_PERSISTED_STAGE_IDS`に属するため、通常は初回の1回しか呼ばれない
    （2回目以降はengineが永続statusを信頼してskipする）。`rerun_upstream=True`が
    このstageをpendingへ戻した場合は、既存の配置済みファイルと突き合わせて
    再配置する（`stage_inputs`自身の冪等性）。
    """
    plan = context["inputs"]
    pipeline_root = context["pipeline_root"]
    snapshot = inputs_mod.stage_inputs(plan, pipeline_root)
    warnings = [{"code": "INPUT_UNVERIFIED", "message": note}
               for note in plan.get("unverified", [])]
    return {"status": "succeeded", "result_refs": [], "warnings": warnings, "error": None,
            "record_updates": {"inputs": snapshot}}


# ---------- upstream（Task4） ----------

def _handle_upstream(context: dict) -> dict:
    """`execution.supervise`を直接呼ぶ（Console用workerを二重起動しない）。

    `console_job_path`は`analysis-job.json`そのものへのパスとして記録する
    （`run_dir = Path(job_path).parent`という、コードベース全体の規約——
    `metabolomix.pipeline.store.register_job_owner`・`metabolomix.console.job_manager`と
    同じ。`metabolomix.pipeline.recovery._console_supervision_state`もこの規約を
    前提にしている）。

    **run_dirはattemptごとに分ける**（`console/attempt-NNNN/`。spec §9.3
    「Consoleの再試行だけは`rerun_upstream=true`の明示を必要とし、**新しい
    job/attemptを作る**」）。固定パスにすると、再実行が前回の終了証跡
    （`execution-result.json`）を上書きし、`msdial.log`を切り詰め
    （`process_control`は`_CREATE_ALWAYS`で開く）、`supervise`が
    `worker.json`の`recovery.befores`（再収集に要る実行前スナップショット）を
    置き換える——timeoutで打ち切った事実そのものが後から証明できなくなる
    （spec §9.1「既存の失敗記録を消さない」）。stageの`attempt`は
    `engine.mark_stage_running`が単調に増やすので、これがそのまま
    「何回目のConsole実行か」になる。
    """
    pipeline_root = context["pipeline_root"]
    request = context["request"]
    inputs_snapshot = context["inputs"]
    attempt = max(1, int(context.get("attempt") or 1))
    run_dir = pipeline_root / _CONSOLE_RUN_SUBDIR / f"attempt-{attempt:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    job_path = run_dir / job_manager.JOB_FILENAME
    dataset_root = pipeline_root / _INPUT_SUBDIR

    # schemaで分ける値は3つだけ。supervise・attemptディレクトリ・
    # write_supervision_inputsの規約は v1/v2 で共通のまま。
    omics, measure, profile_snapshot = "lipidomics", request.get("measure"), None
    if stage_plan.is_v2_request(request):
        # snapshotはruntimeではなく保存済み成果物から読む——再開でprepare_inputsが
        # skipされるとruntimeは空で、そのとき書くjobだけがsnapshotを失う。
        saved = store.read_result_data(
            Path(pipeline_root), context.get("results") or [], "profile")
        if not saved:
            return {"status": "needs_input", "result_refs": [], "warnings": [],
                    "error": {"code": "PROFILE_SNAPSHOT_MISSING",
                              "message": "profile成果物がまだ固定されていません"
                                         "（prepare_inputsを先に通してください）。",
                              "details": {"pipeline_root": str(pipeline_root)}},
                    "record_updates": {}}
        omics = request["omics"]
        measure = saved[-1]["profile"]["processing"]["measure"]
        profile_snapshot = saved[-1]["snapshot"]

    if not job_path.is_file():
        method_rel = inputs_snapshot["method"].get("effective_relative_path")
        method_abs = ((pipeline_root / method_rel) if method_rel
                     else Path(inputs_snapshot["method"]["source_path"]))
        now = _now_iso()
        job = AnalysisJob(
            schema=SCHEMA_VERSION, job_id=f"pipeline_{context['pipeline_id']}",
            status="planned", created_at=now, updated_at=now,
            dataset_root=str(dataset_root),
            input_count=job_manager.count_raw_inputs(dataset_root),
            software_name="MS-DIAL", software_version="", execution_mode="console",
            method_file=str(method_abs), omics=omics,
            polarity=inputs_snapshot["polarity"]["value"], measure=measure,
            run_dir=str(run_dir), save_project=request["save_project"],
            timeout_s=request["timeout_s"],
            # snapshotを持つjobだけがv3として書かれる（handoff/schema.py）。
            profile_snapshot=profile_snapshot,
        )
        job.save(job_path)
        store.register_job_owner(job_path, pipeline_root)
        raw_inventory = [str(p.resolve()) for p in job_manager.list_raw_inputs(dataset_root)]
        console_execution.write_supervision_inputs(run_dir, {
            "raw_inventory": raw_inventory,
            "method_sha256": inputs_snapshot["method"].get("effective_sha256"),
            "exe_path": inputs_snapshot["exe"]["path"],
            "exe_sha256": inputs_snapshot["exe"].get("sha256"),
        })

    receipt = console_execution.supervise(
        job_path, cancel_path=engine.cancel_request_path(pipeline_root))
    job = job_manager.load_job(job_path)

    upstream_update = {
        "console_job_path": str(job_path),
        "execution_id": receipt.get("execution_id"),
        "verification": {"status": job.status, "termination": receipt.get("termination"),
                        "exit_code": receipt.get("exit_code")},
    }
    warnings = [{"code": "CONSOLE_WARNING", "message": w} for w in job.warnings[:5]]

    if job.status == "completed":
        return {"status": "succeeded", "result_refs": [], "warnings": warnings, "error": None,
                "record_updates": {"upstream": upstream_update}}

    message = job.error or (
        f"MS-DIAL Consoleの実行が完了しませんでした"
        f"（termination={receipt.get('termination')}）。")
    return {"status": "failed", "result_refs": [], "warnings": warnings,
            "error": {"code": "MSDIAL_EXECUTION_FAILED", "message": message,
                     "details": {"termination": receipt.get("termination"),
                                 "exit_code": receipt.get("exit_code")}},
            "record_updates": {"upstream": upstream_update}}


# ---------- validate_outputs（Task3） ----------

def _find_gui_project(run_dir: Path) -> Path | None:
    candidates = sorted(run_dir.rglob("*.mdproject"))
    return candidates[0] if candidates else None


def _handle_validate_outputs(context: dict) -> dict:
    """`upstream`が既に完了判定を確定させているので、ここでは`save_project`要求
    に対するGUI project(`*.mdproject`)の有無だけを追加確認する。

    `upstream` stageが失敗した場合、engineはこのstageを呼ばずに止まる
    （`_run_stage_loop`は`needs_input`/`failed`直後に`finish_interrupted`へ
    抜ける）ため、ここに到達する時点でjob.status=="completed"は保証済み。
    """
    request = context["request"]
    if not request.get("save_project"):
        return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None}

    job_path = Path(context["upstream"]["console_job_path"])
    job = job_manager.load_job(job_path)
    project_path = _find_gui_project(Path(job.run_dir))
    if project_path is None:
        return {"status": "succeeded", "result_refs": [], "warnings": [
            {"code": "GUI_PROJECT_UNAVAILABLE",
             "message": "save_project=trueですがGUIプロジェクト(.mdproject)が見つかりません。"}],
            "error": None}

    ref = report_mod.persist_result(context["pipeline_root"], {
        "output_name": "gui_project", "kind": "gui_project",
        "result_id": f"res_gui_project_{context['pipeline_id']}", "path": str(project_path),
    })
    return {"status": "succeeded", "result_refs": [ref], "warnings": [], "error": None}


# ---------- load_dataset（Task7） ----------

def _handle_load_dataset(context: dict) -> dict:
    job_path = Path(context["upstream"]["console_job_path"])
    ds = mztab_loading.load_dataset_state(job_path=job_path)
    context["runtime"]["dataset"] = ds
    return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None}


# ---------- resolve_metadata（Task9） ----------

def _raw_manifest_layout(ds) -> tuple[Path, list[str]]:
    """`ds.assay_sources`（実行時に予定したrawの絶対パス）から、
    sample-manifest.v1検証に使う`source_root`/相対パス一覧を逆算する。

    `metabolomix.analysis.dataset_service._raw_manifest_layout`と同じ規則
    （private関数を跨いで参照しない、既存コードベースの複製方針にならう）。
    """
    paths = [Path(p) for p in (getattr(ds, "assay_sources", None) or {}).values() if p]
    if not paths:
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            "raw参照(assay_sources)が無いため実験情報シートを検証できません。", {})
    root = paths[0].parent if len(paths) == 1 else Path(os.path.commonpath([str(p) for p in paths]))
    expected = [str(p.relative_to(root)) for p in paths]
    return root, expected


def _handle_resolve_metadata(context: dict) -> dict:
    """解決したrole/group/batch/orderを`record["inputs"]["manifest"]`へ残す。

    解決結果を`context["runtime"]`だけに置くと、このworkerプロセスが終わった
    時点で消える。`record["inputs"]["manifest"]`は品質レポートのサンプル来歴
    セクション（`metabolomix.pipeline.report._section_sample_provenance`、spec §11
    が要求し §7.3 が出所の記録を義務付ける）が読む場所そのもので、
    `recovery._absolutize_inputs_paths`（再開時に絶対化。既に絶対なら何もしない
    no-op）も既にこのキーを前提にしている——**読む側が2つあるのに書く側が
    居なかった**。明示シートでも自動生成でも同じ形の行が出るので、両方ここで
    記録する。
    **`store._relativize_inputs_paths`はここでは効かない**——`create_run`
    （このstageより前）でしか呼ばれず、`save_run`自身は`record`を
    deepcopyしてそのまま書くだけ（相対化を挟まない）。つまりこの`manifest`が
    ここで初めて記録される`source_file`は保存後も絶対パスのまま残る。
    結果としてrun recordはpipeline_rootを移設しても`source_file`が指す実体と
    ずれ、品質レポートの来歴テーブルにホスト側の絶対パスがそのまま出る
    （最終レビューで判明・controller裁定によりそのまま出荷、
    `docs/superpowers/notes/2026-09-05-raw-folder-pipeline-validation.md`
    「既知の制限」参照）。
    """
    ds = context["runtime"]["dataset"]
    request = context["request"]
    source_root = Path(context["identity"]["source_root"])
    manifest_path = _resolve_manifest_path(source_root, request)

    rows = None
    if manifest_path is not None:
        raw_root, expected_sources = _raw_manifest_layout(ds)
        rows = parse_manifest(manifest_path, source_root=raw_root, expected_sources=expected_sources)
    metadata = resolve_sample_metadata(ds, rows)
    context["runtime"]["metadata"] = metadata
    # シートの**中身**の指紋も残す。`recovery.prepare_resume`はこれを今の
    # ファイルと突き合わせて「同じパスのシートを書き直した訂正」を検出する
    # ——パス文字列の比較だけでは、書き直しは変更として現れない。
    inputs_update = {
        **(context.get("inputs") or {}),
        "manifest": metadata,
        "manifest_source": inputs_mod.manifest_source_record(source_root, request),
    }
    return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None,
            "record_updates": {"inputs": inputs_update}}


# ---------- preprocess（Task11） ----------
#
# **計算結果のrefは`provenance.result_id`そのもので登録する**（spec §6.1/§9.1）。
# 合成id（`res_preprocess_<fingerprint>`等）で登録すると、下流のrefが
# `parent_ids`へ入れる**計算側の実ID**（`result_state.new_provenance`が発行する
# UUID）と別の名前空間になり、`record["results"]`の中で親を解決できなくなる。
# 図やTSVのような「成果物ファイル」だけは、どの必須出力かで引ける合成id
# （`res_pca_figure_*` / `res_volcano_<cid>` / `res_tsv_<cid>`）のまま——
# これらは計算結果ではなく、`parent_ids`で親の計算結果を名指しする側なので、
# 誰かの`parent_ids`に現れることがない。

def _handle_preprocess(context: dict) -> dict:
    ds = context["runtime"]["dataset"]
    metadata = context["runtime"]["metadata"]
    requested = context["request"].get("preprocess") or {}
    request_revision = context["request_meta"].get("revision")

    plan = dataset_service.preprocess_auto(ds, requested, metadata, request_revision=request_revision)

    provenance = plan["result"]["provenance"]
    data = {k: v for k, v in plan.items() if k != "result"}
    data["result"] = {k: v for k, v in plan["result"].items() if k != "provenance"}
    ref = report_mod.persist_result(context["pipeline_root"], {
        "output_name": "preprocess", "kind": "preprocess",
        "result_id": provenance["result_id"],
        "parent_ids": list(provenance.get("parent_ids") or []),
        "data": data, "request_revision": request_revision,
    })
    warnings = [{"code": "PREPROCESS_CAVEAT", "message": w}
               for w in plan["result"].get("caveats", [])]
    return {"status": "succeeded", "result_refs": [ref], "warnings": warnings, "error": None}


# ---------- pca（Task6/Task8, R18） ----------

def _handle_pca(context: dict) -> dict:
    """PreconditionError（PCA不成立）は例外のままengineへ渡さず、ここで
    直接needs_inputのStageResultへ変換する（R18）。

    engineの`_invoke_handler`は`DomainError`のうち`_NEEDS_INPUT_CODES`
    しかneeds_inputへ変換しない——`PreconditionError`は別クラスなので、
    捕らえずに投げると汎用`except Exception`分岐に落ちて`failed`になる。
    ここでhandlerが直接dictを返すことで、その分岐を経由させない。
    """
    ds = context["runtime"]["dataset"]
    pipeline_root = context["pipeline_root"]
    request_revision = context["request_meta"].get("revision")

    try:
        result = dataset_service.pca_dataset(ds, request_revision=request_revision)
    except PreconditionError as exc:
        return {"status": "needs_input", "result_refs": [], "warnings": [],
                "error": {"code": "PCA_PRECONDITION_MISSING", "message": exc.message,
                         "details": dict(exc.details)}}

    fingerprint = result["provenance"]["input_fingerprint"]
    # `scores`（検体ごとのPC座標）は残す——これが図の中身そのもので、これを落とすと
    # 永続化されたPCA結果から図を描き直せない（`result_output._render`は
    # `points`が無ければ`scores`から射影する）。`loadings`だけは特徴量×成分の
    # 大きな行列なので落とす。
    summary = {k: v for k, v in result.items() if k not in ("provenance", "loadings")}
    data_ref = report_mod.persist_result(pipeline_root, {
        "output_name": "pca", "kind": "pca", "result_id": result["provenance"]["result_id"],
        "data": summary, "parent_ids": list(result["provenance"].get("parent_ids") or []),
        "request_revision": request_revision,
    })

    figure_path = pipeline_root / _RESULTS_SUBDIR / "pca.png"
    result_output.save_result_figure(ds, result, figure_path, kind="pca")
    figure_ref = report_mod.persist_result(pipeline_root, {
        "output_name": "pca_figure", "kind": "pca_figure",
        "result_id": f"res_pca_figure_{fingerprint[:24]}", "path": str(figure_path),
        "parent_ids": [result["provenance"]["result_id"]], "request_revision": request_revision,
    })

    warnings = [{"code": "PCA_CAVEAT", "message": w}
               for w in result["provenance"].get("warnings", [])]
    return {"status": "succeeded", "result_refs": [data_ref, figure_ref],
            "warnings": warnings, "error": None}


# ---------- resolve_comparisons（Task12ゲート） ----------

def _handle_resolve_comparisons(context: dict) -> dict:
    """空配列を`COMPARISON_REQUIRED`とし、各比較の方向・群・交絡を検証してから
    比較別stageへ進む（`_NEEDS_INPUT_CODES`が`COMPARISON_REQUIRED`をneeds_input
    へ変換する）。実際の統計計算（`compare_dataset`）はここでは行わない
    ——`differential`/`export` handlerがそれぞれ独立に`run_comparison`を
    呼ぶため、ここでの検証は「先に全件止めるためのゲート」に留める。
    """
    ds = context["runtime"]["dataset"]
    metadata = context["runtime"]["metadata"]
    comparisons = context["request"].get("comparisons") or []
    if not comparisons:
        raise DomainError(
            "COMPARISON_REQUIRED",
            "comparisonsが指定されていません。差次的解析には最低1件の比較定義が必要です。", {})

    for comparison in comparisons:
        resolved = dataset_service.resolve_comparison(ds, comparison, metadata)
        confounding = resolved["confounding"]
        if confounding["confounded"] and not resolved["allow_confounded"]:
            raise DomainError(
                "CONFOUNDED_COMPARISON",
                f"{resolved['reference_group']!r}と{resolved['test_group']!r}は"
                "群とバッチが完全に交絡しています。allow_confounded=trueを"
                "明示しない限り差次的解析は実行しません。",
                {"comparison_id": comparison["comparison_id"]})
    return {"status": "succeeded", "result_refs": [], "warnings": [], "error": None}


def _find_comparison(request: dict, comparison_id: str) -> dict:
    for comparison in request.get("comparisons") or []:
        if comparison["comparison_id"] == comparison_id:
            return comparison
    raise DomainError(
        "COMPARISON_REQUIRED",
        f"comparison_id={comparison_id!r}が要求のcomparisonsに見つかりません。",
        {"comparison_id": comparison_id})


# ---------- differential（Task12） ----------

def _persist_differential(context: dict, cid: str, result: dict) -> dict:
    """差次的結果を永続化し、同じworker内の`export`が読めるようruntimeへ残す。

    refのresult_idは`provenance.result_id`そのもの——volcano/TSVの`parent_ids`と
    TSVメタ行の`# result_id`が名指しするのはこのIDで、`record["results"]`の中で
    解決できなければ出所を辿れない（spec §6.1/§9.1）。
    """
    ref = report_mod.persist_result(context["pipeline_root"], {
        "output_name": f"differential:{cid}", "kind": "differential",
        "result_id": result["provenance"]["result_id"], "data": result,
        "parent_ids": list(result["provenance"].get("parent_ids") or []),
        "request_revision": context["request_meta"].get("revision"),
    })
    context["runtime"].setdefault("differential", {})[cid] = result
    return ref


def _handle_differential(context: dict) -> dict:
    ds = context["runtime"]["dataset"]
    metadata = context["runtime"]["metadata"]
    cid = context["comparison_id"]
    comparison = _find_comparison(context["request"], cid)

    result = dataset_service.run_comparison(ds, comparison, metadata)
    ref = _persist_differential(context, cid, result)
    warnings = [{"code": "DIFFERENTIAL_CAVEAT", "message": w}
               for w in result.get("caveats", [])]
    return {"status": "succeeded", "result_refs": [ref], "warnings": warnings, "error": None}


# ---------- export（Task8） ----------

def _handle_export(context: dict) -> dict:
    """volcano図とcontract TSVを書く。

    使う差次的結果は、**同じworker passで`differential` handlerが計算し永続化
    したそのもの**（`runtime["differential"][cid]`）。`run_comparison`を独立に
    呼び直すと`result_state.new_provenance`が新しいUUIDを発行するため、
    volcano/TSVの`parent_ids`とTSVメタ行の`# result_id`が
    `record["results"]`に載っていない2つ目の結果を指してしまい、出所の連結が
    切れる（spec §6.1/§9.1、E02「result_idの整合」）。

    worker再起動をまたいで`differential`が既にsucceededでskipされ、`export`
    だけがこのpassで動く場合はruntimeに結果が無い。永続化済みJSONは前回の
    `ds.preprocess_id`を親に持ち`assert_current`を通らない（`ds`は
    `_ALWAYS_RECONSTRUCT_STAGE_IDS`により作り直され、preprocessのIDも新しい）
    ので、その場合だけ再計算する——ただし**再計算した結果も
    `differential:<cid>`として永続化する**（append-only）。そうしないと、
    今書くvolcano/TSVが名指しするIDだけがrecordに存在しない状態に戻る。
    """
    ds = context["runtime"]["dataset"]
    metadata = context["runtime"]["metadata"]
    cid = context["comparison_id"]
    comparison = _find_comparison(context["request"], cid)
    pipeline_root = context["pipeline_root"]
    request_revision = context["request_meta"].get("revision")

    result_refs: list[dict] = []
    result = (context["runtime"].get("differential") or {}).get(cid)
    if result is None:
        result = dataset_service.run_comparison(ds, comparison, metadata)
        result_refs.append(_persist_differential(context, cid, result))

    volcano_path = pipeline_root / _RESULTS_SUBDIR / f"volcano_{cid}.png"
    result_output.save_result_figure(ds, result, volcano_path, kind="volcano")
    volcano_ref = report_mod.persist_result(pipeline_root, {
        "output_name": f"volcano:{cid}", "kind": "volcano",
        "result_id": f"res_volcano_{cid}", "path": str(volcano_path),
        "parent_ids": [result["provenance"]["result_id"]], "request_revision": request_revision,
    })

    result_refs.append(volcano_ref)
    warnings: list[dict] = []
    tsv_path = pipeline_root / _RESULTS_SUBDIR / f"differential_{cid}.tsv"
    try:
        dataset_export.export_dataset_result(ds, result, tsv_path)
    except DomainError as exc:
        if exc.code != "NO_ANNOTATED_FEATURES":
            raise
        warnings.append({"code": "EXPORT_BACKGROUND_EMPTY", "message": exc.message})
    else:
        tsv_ref = report_mod.persist_result(pipeline_root, {
            "output_name": f"tsv:{cid}", "kind": "tsv",
            "result_id": f"res_tsv_{cid}", "path": str(tsv_path),
            "parent_ids": [result["provenance"]["result_id"]], "request_revision": request_revision,
        })
        result_refs.append(tsv_ref)

    return {"status": "succeeded", "result_refs": result_refs, "warnings": warnings, "error": None}


# ---------- report（Task17） ----------

def _handle_report(context: dict) -> dict:
    """レポートは完全な`record`（`warnings`/`needs_input`を含む）を必要とする
    ため、`context`から組み立てず`store.load_run`で読み直す（読取専用）。
    """
    pipeline_root = context["pipeline_root"]
    record = store.load_run(pipeline_root)
    report_path = pipeline_root / "reports" / "pipeline-quality-report.md"
    summary = report_mod.write_pipeline_report(record, report_path)

    ref = report_mod.persist_result(pipeline_root, {
        "output_name": "quality_report", "kind": "quality_report",
        "result_id": f"res_quality_report_{context['identity']['pipeline_id']}",
        "path": str(report_path),
    })

    # `summary["status"]`は`write_pipeline_report`内部の`evaluate_target`が
    # このrefをまだ`record["results"]`へ追記する前に呼んだ判定であり、
    # quality_report自身は必ず未達成として出る（意図通り。markdownの
    # `status_at_report_time`はこの生成時点のstatusを正直に書くための場所）。
    # ここでの`REPORT_INCOMPLETE_AT_WRITE_TIME`警告は「quality_report自身が
    # まだ無い」以外の理由で本当に不完全な場合だけ意味を持つ——refを積んだ
    # 仮想recordでもう一度evaluate_targetを呼び直し、quality_report自身の
    # 未生成を理由にした偽陽性（レビュー指摘4）を除く。
    record_with_report = copy.deepcopy(record)
    record_with_report["results"] = list(record_with_report.get("results") or []) + [ref]
    final_evaluation = report_mod.evaluate_target(record_with_report)

    warnings: list[dict] = []
    if final_evaluation["status"] != "completed":
        warnings.append({
            "code": "REPORT_INCOMPLETE_AT_WRITE_TIME",
            "message": f"レポート生成時点でstatus={final_evaluation['status']}でした。",
        })
    return {"status": "succeeded", "result_refs": [ref], "warnings": warnings, "error": None}
