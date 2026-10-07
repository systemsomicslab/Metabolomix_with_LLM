"""再開・再構築・取消・中断の読取表示（spec §9.2/9.3）。

このモジュールが解く問題は3つ、公開関数もその3つに対応する。

1. **`read_status`**: pipeline-run.jsonを読むだけの、絶対に書き換えないstatus取得。
   `status="running"`のときだけ、記録済みの`worker`（`metabolomix.pipeline.engine.
   run_engine`が刻む自分自身のidentity）が今も生きているかを`same_process`で見る。
   このprobeがどう転んでも（生存・消失・判定不能）ファイルには一切触れない
   ——observed_healthという別軸で返すだけで、statusフィールド自体は変えない。
2. **`request_cancel`**: 協調的な取消フラグ（`metabolomix.pipeline.engine.
   cancel_request_path`が指す小さなJSON）を書いて受理を返すだけ。ここでは
   一切のプロセスを直接殺さない——上流のConsole停止は`supervise`側の
   `cancel_path`監視（既存機構）に、下流の停止はengineのstage境界チェックに
   委ねる。受理（このファイルが書けたこと）と、実際にworkerが止まったこと
   （`cancelled`への確定）は別状態として扱う。
3. **`prepare_resume`**: 新しいattempt/revisionを用意し、どのstageを`pending`へ
   戻すかを決めて保存するだけ——**launchしない**。以下の順で安全装置を通す:

   a. 監視を失ったまま稼働中のConsoleがあれば`EXECUTION_UNRESOLVED`（終了を
      確認できないものを勝手にverified completedへ昇格しない）。
   b. 上流が未完了（succeeded以外）なのに`rerun_upstream`が無ければ
      `UPSTREAM_RERUN_REQUIRED`（自動再試行はしない）。
   c. `rerun_upstream=False`のときだけ、固定済み入力（raw stat・実効method
      ハッシュを含む）を再検証する——Task13の`inputs.verify_inputs`に加え、
      Consoleが実際に読む`inputs/effective-method.txt`のhashもここで照合する
      （Task13の申し送り。`verify_inputs`自体は原本methodの原本ハッシュしか
      見ない）。
   d. `updates`を`metabolomix.pipeline.request.merge_updates`で検証・統合し、
      内容が変わっていなければ新しいrequest revisionを作らない。
   e. Task 6の依存区分に倣い、更新内容からどのstageを`pending`へ戻すかを決める
      （`_stages_to_reset`）。
   f. `planned`へ戻したなら、協調取消フラグ（`control/cancel-request.json`）を
      取り下げる（`_clear_cancel_request`）。消し忘れると、再開したworkerが
      最初のstage境界で`cancel_requested()`を見て即座に止まる。生きている
      workerを観測できたときは消さない——出された取消を黙って無効にしない。

   永続statusが`running`のまま固まったrun（worker喪失）も、そのworkerが
   もう存在しないと観測できたときだけ`planned`へ戻す（`_RESUMABLE_RUNNING_HEALTH`。
   `read_status`が`PIPELINE_INTERRUPTED`で案内する再開経路の実体）。

   `request_id`を指定した再送は、直前と同じ`updates`/`rerun_upstream`なら
   新しいrevision/attemptを作らず直前の結果を返す（`resume_log`。source_root単位の
   受付索引＝`store.find_or_create_run`とは別物で、pipeline_root内に閉じる）。
   `request_id`を指定しない再送でも、統合後の内容が現行のrequestと同一なら
   同様に新revisionを作らない。

このモジュールはsessionを一切importしない（`metabolomix.core.session_state` /
`metabolomix.core.mcp_core` / `metabolomix.tools.*`）——read_status/request_cancel/
prepare_resumeはすべてMCP接続やグローバル状態と無関係に、ファイルだけで完結する。
"""
from __future__ import annotations

import copy
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.console import execution as console_execution
from metabolomix.core.atomic_io import DomainError, atomic_write_json, canonical_hash
from metabolomix.core.process_control import same_process
from metabolomix.pipeline import engine
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import stage_plan
from metabolomix.pipeline import store

__all__ = ["normalize_pipeline_root", "prepare_resume", "read_status", "request_cancel"]

#: prepare_resumeの状態競合吸収リトライ上限（store.pyの
#: _REQUEST_ID_PATCH_MAX_ATTEMPTSと同じ考え方——有限回の読み直しで十分安全に
#: 吸収できる。無関係な同時更新のせいでresume受付全体を失敗させない）。
_RESUME_RETRY_MAX_ATTEMPTS = 5

#: 中断由来の終端状態。`completed`と違い、何も変えなくてもresumeが呼ばれた
#: 時点で常に`planned`へ戻してよい（D09「取消後・timeout後のresume」）。
_INTERRUPTED_STATUSES = frozenset({"partial", "failed", "cancelled", "needs_input"})

#: 永続statusが`running`のまま残ったrunを「中断」と見なしてよい観測結果。
#:
#: workerを失うと（強制終了・OSごと落ちる）statusは`running`のまま固まる。
#: `read_status`はこれを`observed_health="worker_missing"`＋
#: `PIPELINE_INTERRUPTED`（「prepare_resumeで再開してください」）として報告する
#: のだから、`prepare_resume`はその案内どおりに`planned`へ戻せなければならない
#: （spec A06/§9.3）。
#:
#: **`"unknown"`は入れない。** identityが記録されていない・生存確認そのものが
#: 失敗した、のいずれも「死んでいる証拠」ではない（brief「liveness probeの失敗を
#: deadと解釈しない」）。生きているworkerのrunを`planned`へ戻して別workerを
#: 起こしにいくのは、記録の書き換えと二重実行の両方を招く。
_RESUMABLE_RUNNING_HEALTH = frozenset({"worker_missing"})


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_pipeline_root(path: Path) -> Path:
    """`path`がpipeline-run.jsonそのものでも、その親（pipeline_root）でも
    受け付ける。

    Task14/15の他モジュール（store.load_run/save_run、engine.run_engine）は
    一貫してpipeline_root（ディレクトリ）を受け取る。このモジュールの3公開
    関数も同じ流儀に合わせるが、`read_status`をfile-levelで検証したい
    呼び出し側（このファイルの`.read_bytes()`で「一切書き換えていない」ことを
    直接証明できる）にも応えられるよう、run.jsonファイルのパスを渡された
    場合はその親を使う。

    公開関数にしているのは、`pipeline_status`/`pipeline_cancel`が既に
    アドバタイズしているこの2形式受け入れを、`metabolomix.pipeline.service`の
    `resume_pipeline`でも同じ規則で揃えるため（レビュー指摘3: 正規化を
    `prepare_resume`の内側だけで完結させると、呼び出し元がその後の
    `launch_pipeline_worker`/receiptで生のパスを使ってしまい、二重に
    `pipeline-run.json`を連結した不正なパスを作ってしまう）。
    """
    p = Path(path)
    if p.is_file():
        return p.parent
    return p


# ---------- read_status ----------

def _stage_statuses(record: dict) -> dict:
    return {sid: stage.get("status") for sid, stage in record.get("stages", {}).items()}


def _worker_health(record: dict) -> str:
    """記録済みworker identityの生死を`same_process`で観測する（読取専用）。

    `"ok"`（生きている）/ `"worker_missing"`（記録されたidentityはもう存在
    しない）/ `"unknown"`（identityが無い、または生存確認自体が失敗した）。
    probeがどう転んでもファイルには一切触れない——`read_status`（表示）と
    `prepare_resume`（再開判断）が同じ観測を共有するためだけの純関数。
    """
    identity = (record.get("worker") or {}).get("identity")
    if not identity:
        return "unknown"
    try:
        alive = same_process(identity)
    except Exception:  # noqa: BLE001 - probe自体の失敗をdeadと解釈しない
        return "unknown"
    return "ok" if alive else "worker_missing"


def read_status(path: Path, *, include_details: bool = False) -> dict:
    """pipeline状態を読むだけで、一切書き換えない（brief「read-onlyとする」）。

    `status`が"running"のときだけ、記録済みworker identityの生死を`same_process`
    で見る。この生存確認自体が失敗しても（例外・判定不能）"dead"とは読まず、
    `observed_health="unknown"`として呼び出し側に復旧判断を委ねる——
    「liveness probeの失敗をdeadと解釈しない」というbriefの制約そのもの。
    """
    pipeline_root = normalize_pipeline_root(path)
    record = store.load_run(pipeline_root)  # 読取専用。ここでは一切保存しない。

    status = record.get("status")
    observed_health = _worker_health(record) if status == "running" else "ok"

    result = {
        "pipeline_id": record.get("identity", {}).get("pipeline_id"),
        "status": status,
        "effective_target": record.get("request", {}).get("effective_target"),
        "observed_health": observed_health,
        "needs_input": record.get("needs_input"),
        "warnings": list(record.get("warnings") or []),
        "stage_statuses": _stage_statuses(record),
    }
    if status == "running" and observed_health in {"worker_missing", "unknown"}:
        result["recovery_hint"] = {
            "code": "PIPELINE_INTERRUPTED" if observed_health == "worker_missing"
                    else "SUPERVISION_UNKNOWN",
            "message": "workerの生存を確認できません。prepare_resumeで再開してください。",
        }
    if include_details:
        result["record"] = _without_library_location(record)
    return result


def _without_library_location(record: dict) -> dict:
    """戻り値用の record の写しから、研究室 MSP の置き場所を除く（ファイル名だけ残す）。

    戻り値は LLM の文脈に入る（spec 2026-10-07 §3.3）。保存された record は変えない
    ——resume の再検査（`inputs.verify_inputs`）が絶対パスを要る。
    """
    inputs = record.get("inputs")
    if not isinstance(inputs, dict) or not isinstance(inputs.get("msp"), dict):
        return record
    redacted = copy.deepcopy(record)
    msp = redacted["inputs"]["msp"]
    if msp.get("path"):
        msp["path"] = Path(msp["path"]).name
    overrides = (redacted["inputs"].get("method") or {}).get("overrides")
    if isinstance(overrides, dict):
        for key, value in overrides.items():
            if key.lower() == "msp file path" and value:
                overrides[key] = Path(value).name
    return redacted


# ---------- request_cancel ----------

def request_cancel(path: Path) -> dict:
    """協調的な取消フラグを保存し、受理だけを返す（brief「cancelledは実際に
    workerが停止を確認した後に確定する」）。

    ここでは一切のプロセスを直接殺さない。上流Consoleの停止はTask18の
    upstream handlerが`supervise`へ`cancel_path=engine.cancel_request_path(...)`
    を渡すことで、既存の`_monitor`が拾う（同じファイルの存在を見る）。
    Job Objectのowner processが不明な状況でこの関数が別PIDを終了させることは
    無い——それ自体をしないという設計でその要件を満たす。
    """
    pipeline_root = normalize_pipeline_root(path)
    record = store.load_run(pipeline_root)  # 存在・schema確認（読取専用）。

    requested_at = _now_iso()
    atomic_write_json(engine.cancel_request_path(pipeline_root),
                      {"cancel_requested": True, "requested_at": requested_at})

    status = read_status(pipeline_root)
    return {
        "pipeline_id": record.get("identity", {}).get("pipeline_id"),
        "accepted": True,
        "cancel_requested_at": requested_at,
        "status": status["status"],
        "observed_health": status["observed_health"],
    }


# ---------- prepare_resume: 補助関数 ----------

def _load_full_request(pipeline_root: Path, record: dict) -> dict:
    saved_path = record["request"]["saved_path"]
    return json.loads((Path(pipeline_root) / saved_path).read_text(encoding="utf-8"))


def _console_supervision_state(record: dict) -> tuple[dict | None, bool]:
    """(console process identity, 検証済みか)を返す。

    `verified`はpipelineの`upstream`stage自体がsucceededかどうか
    （＝validate_outputsまで含め検証済みかどうか）で判定する。まだ
    Console jobを起動していない（`console_job_path`が無い）ならidentityは
    Noneで返す——これは「そもそも確認すべき対象が無い」正当な状態。

    `console_job_path`は本コードベース全体の流儀（`metabolomix.pipeline.store.
    register_job_owner`・`metabolomix.console.job_manager.load_job`等）に合わせ、
    analysis-job.jsonそのものへのパスとして扱う——終了証跡
    （`execution-result.json`）はその親（run_dir）に置かれる。**ただしTask18
    （build_handlers）はまだ無く、この規約は本タスクの推測にすぎない
    （report「懸念3」参照）。**

    Task18が逆の規約（`console_job_path`をrun_dirそのものとして書く等）を
    選んだ場合、ここで`console_job_path`は記録されているのに終了証跡が
    読めない・想定の形をしていない、という状況が起こる。これを
    「`console_job_path`が最初から記録されていない」場合と同じ
    `(None, verified)`で返すと、`prepare_resume`のEXECUTION_UNRESOLVEDゲート
    （呼び出し元）を素通りしてしまい、監視を失って稼働中かもしれない
    Consoleが黙ってverified扱いされかねない（spec §9.3「その実行を
    verified completedへ昇格しない」）。

    そこで、`upstream`stageがまだ`verified`でない（＝succeededと永続化
    されていない）のに`console_job_path`は記録されている、という組合せで
    終了証跡が読めない・process_identityを持たない場合は、ここで
    `DomainError("EXECUTION_UNRESOLVED", ...)`を送出し、呼び出し元
    （`prepare_resume`）へ「規約不一致・判定不能」を騒がしく伝える
    （`verified`が既にTrueなら、record自体が既にsucceededと確定させている
    ので、受信できない終了証跡があっても実害は無く例外にしない）。
    """
    # 上流stageはv1が`upstream`、v2が`execute_console`。ここを直書きにすると
    # v2では常にverified=Falseになり、succeededと永続化済みの実行まで
    # EXECUTION_UNRESOLVEDで止まる（prepare_resumeは既に_upstream_stage_idを使う）。
    verified = record.get("stages", {}).get(
        _upstream_stage_id(record), {}).get("status") == "succeeded"
    job_path = (record.get("upstream") or {}).get("console_job_path")
    if not job_path:
        return None, verified
    run_dir = Path(job_path).parent
    try:
        raw = console_execution.receipt_path(run_dir).read_text(encoding="utf-8")
        data = json.loads(raw)
    except (OSError, ValueError) as exc:
        if verified:
            return None, verified
        raise DomainError(
            "EXECUTION_UNRESOLVED",
            "console_job_pathは記録されていますが終了証跡を読めません"
            f"（規約不一致の疑い、Task18確認要）: {run_dir}",
            {"console_job_path": str(job_path), "run_dir": str(run_dir)},
        ) from exc
    identity = data.get("process_identity")
    if not isinstance(identity, dict):
        if verified:
            return None, verified
        raise DomainError(
            "EXECUTION_UNRESOLVED",
            f"終了証跡がprocess_identityを持たない不正な形式です: {run_dir}",
            {"console_job_path": str(job_path), "run_dir": str(run_dir)},
        )
    return identity, verified


def _absolutize_inputs_paths(inputs: dict, pipeline_root: Path) -> dict:
    """`store._relativize_inputs_paths`の逆操作。

    method/lbm/exe・manifest[].source_fileがpipeline_root相対に固定されている
    場合だけ、再び絶対パスへ戻す（`inputs.verify_inputs`はこれらを絶対パスと
    して扱う）。
    """
    pipeline_root = Path(pipeline_root).resolve()

    def _fix(value: str) -> str:
        if isinstance(value, str) and not Path(value).is_absolute():
            return str((pipeline_root / value).resolve())
        return value

    result = copy.deepcopy(inputs)
    for top, sub in (("method", "source_path"), ("lbm", "path"), ("exe", "path")):
        section = result.get(top)
        if isinstance(section, dict) and isinstance(section.get(sub), str):
            section[sub] = _fix(section[sub])
    manifest_rows = result.get("manifest")
    if isinstance(manifest_rows, list):
        for row in manifest_rows:
            if isinstance(row, dict) and isinstance(row.get("source_file"), str):
                row["source_file"] = _fix(row["source_file"])
    result["pipeline_root"] = str(pipeline_root)
    return result


def _verify_effective_method(snapshot: dict, pipeline_root: Path) -> None:
    """Task13の申し送り: `inputs.verify_inputs`は原本methodの原本ハッシュしか
    見ず、実際にConsoleが読む`inputs/effective-method.txt`のhashを再検証しない。
    ここで`method.effective_relative_path`/`effective_sha256`（Task13
    `stage_inputs`が固定した値）を実ファイルと照合する。
    """
    method = snapshot.get("method") or {}
    rel = method.get("effective_relative_path")
    expected = method.get("effective_sha256")
    if not rel or not expected:
        return  # 最小fixture等、実効メソッドを固定していないrunは対象外
    target = Path(pipeline_root) / rel
    try:
        actual = _sha256_file(target)
    except OSError as exc:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"実効メソッドファイルが読めません（再開時の照合）: {target}",
            {"path": str(target)},
        ) from exc
    if actual != expected:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"実効メソッドファイルが元の固定内容と一致しません（再開時の照合）: {target}",
            {"path": str(target)},
        )


def _verify_staged_inputs(record: dict) -> None:
    """再開時の入力再検証（D05）。`rerun_upstream=True`のときは呼ばない
    ——上流をやり直す以上、`prepare_input`が自分で再検証・再配置する。
    """
    pipeline_root = Path(record["identity"]["pipeline_root"])
    snapshot = _absolutize_inputs_paths(record.get("inputs") or {}, pipeline_root)
    if not snapshot.get("source_root") or not snapshot.get("raw_stat"):
        return  # 検証対象を持たない最小fixture（単体test）向けの安全側スキップ
    inputs_mod.verify_inputs(snapshot)
    _verify_effective_method(snapshot, pipeline_root)


def _diff_comparisons(old: list, new: list) -> set:
    """変更・新規追加されたcomparison_idの集合を返す（内容が完全一致するものは
    含めない——B04「groupだけ訂正」で、無関係なcomparisonまで巻き込まないため）。
    """
    old_by_id = {c["comparison_id"]: c for c in (old or [])}
    changed = set()
    for c in new or []:
        cid = c["comparison_id"]
        if old_by_id.get(cid) != c:
            changed.add(cid)
    return changed


def _manifest_content_changed(record: dict, request: dict) -> bool:
    """実験情報シートの**中身**が、前回解決したときから変わったかを返す。

    要求の`sample_manifest`はパス文字列でしかないため、同じ`sample-manifest.tsv`を
    書き直して`pipeline_resume`しても、要求の比較では「変更なし」に見える
    ——群の訂正が下流へ一切波及せず、旧群割当のまま`completed`になり、
    訂正が無視されたことすら利用者に伝わらない。

    比較の基準は`resolve_metadata` handlerが記録した
    `record["inputs"]["manifest_source"]["sha256"]`（そのとき実際に読んだ
    シートの内容hash）。記録が無い版のrun・自動一覧生成（シート無し）では
    比較材料が無いので`False`——「判定できない」を「変わった」と読み替えて
    下流を無条件にやり直させない。
    """
    recorded = (record.get("inputs") or {}).get("manifest_source") or {}
    previous = recorded.get("sha256")
    if not previous:
        return False
    source_root = Path(record["identity"]["source_root"])
    current = inputs_mod.manifest_source_record(source_root, request)
    return current.get("sha256") != previous


def _diff_statistics(old: list, new: list) -> set:
    """変更・追加・**削除**されたstatistic_idの`statistics:<id>`token集合を返す。

    削除も含めるのがcomparisonとの違い（spec §6.2「統計の追加・削除でもstage集合を
    再構築し、削除済み統計の成果物を必須出力に残さない」）——削除された統計の
    stageはrecordには残り続けるので、そこをpendingへ戻してcurrentから外す。
    """
    old_by_id = {item["statistic_id"]: item for item in (old or [])}
    new_by_id = {item["statistic_id"]: item for item in (new or [])}
    changed = {sid for sid, item in new_by_id.items() if old_by_id.get(sid) != item}
    changed |= set(old_by_id) - set(new_by_id)
    return {f"statistics:{sid}" for sid in changed}


def _profile_arguments(request: dict) -> dict:
    """v2要求のときだけ、profileを読んで`merge_updates`へ渡す引数を作る。

    v2の更新検証（統計の形・target参照・routineの上書き範囲）はprofileが無いと
    成立しない。読み込みはここ（再開の受付）の責務で、`merge_updates`自身は
    検証済みのdictしか受け取らない——`service._profile_arguments`と同じ規約を
    再開経路でも守る（片方だけがprofileを見ると、初回とresumeで検証の厳しさが
    食い違う）。
    """
    if not stage_plan.is_v2_request(request):
        return {}
    from metabolomix.console import profiles as profiles_mod

    profile_file = request.get("profile_file")
    if not profile_file:
        raise DomainError(
            "PIPELINE_REQUEST_INVALID",
            "v2要求にprofile_fileがありません（再開時の検証にprofileが要ります）。",
            {"schema": request.get("schema")})
    purpose = request.get("execution_purpose", "routine")
    profile = profiles_mod.load_profile(Path(profile_file), purpose)
    routine_overrides = (
        profiles_mod.certificate_routine_overrides(profile, Path(profile_file))
        if purpose == "routine" else None)
    return {"profile": profile, "routine_overrides": routine_overrides}


def _changed_aspects_v2(current_request: dict, merged_request: dict, *,
                        manifest_content_changed: bool) -> set:
    """v2要求の差分を`stage_plan.invalidated_v2`が読むtoken集合へ畳む。

    `qc_raw_scope`（元データQCの評価集合を動かすrecipe変更）はここでは立てない
    ——要求の差分だけでは「そのrecipe変更がQC対象集合を動かすか」を判定できない。
    判定できるのは実際に集合を組む工程側（Task 12のhandler）で、そこが必要なら
    同じtokenを渡す。
    """
    changed: set = set()
    if manifest_content_changed:
        changed.add("metadata")
    for key in ("sample_manifest", "preprocess", "standard_assays",
                "feature_bindings", "target"):
        if current_request.get(key) != merged_request.get(key):
            changed.add(key)
    changed |= _diff_statistics(current_request.get("statistics"),
                                merged_request.get("statistics"))
    return changed


def _stages_to_reset(current_request: dict, merged_request: dict, *,
                     rerun_upstream: bool, stage_ids: set,
                     manifest_content_changed: bool = False) -> set:
    """Task 6の依存区分に倣い、更新内容からどのstageを`pending`へ戻すかを返す。

    `load_dataset`/`resolve_metadata`/`preprocess`/`pca`は対象にしない——
    `engine._ALWAYS_RECONSTRUCT_STAGE_IDS`により、resumeのたびに無条件で
    handlerを呼び直すため、明示的にpendingへ戻す必要が無い。

    `manifest_content_changed`は「同じパスのシートを書き直した訂正」
    （`_manifest_content_changed`）。要求の値としては何も変わっていないが、
    群・バッチ・include が変わりうる以上、`sample_manifest`を差し替えたときと
    同じ範囲を差し戻す。

    `pipeline-request.v2`（メタボロミクス）の依存区分はspec §6.2が定めており、
    その表は`metabolomix.pipeline.stage_plan.invalidated_v2`が唯一の正準
    （store/engineが使うstage builderと同じモジュール）。ここでは要求の差分を
    そのtoken集合へ畳んで渡すだけで、区分そのものを持たない。
    """
    if rerun_upstream:
        # 上流のやり直しは下流すべてに波及する（Consoleの出力自体が変わりうる）。
        return set(stage_ids)

    if stage_plan.is_v2_request(merged_request):
        # v2は依存表も`stage_plan`が正準（store/engineと同じモジュール）。
        # 戻り値はrecordが実際に持つstageへ絞る——削除済み統計のstage IDのように
        # 新しい計画には無いがrecordには残っているものを含みうるため。
        changed = _changed_aspects_v2(
            current_request, merged_request,
            manifest_content_changed=manifest_content_changed)
        return stage_plan.invalidated_v2(changed, merged_request) & set(stage_ids)

    reset: set = set()
    metadata_or_preprocess_changed = (
        manifest_content_changed
        or current_request.get("sample_manifest") != merged_request.get("sample_manifest")
        or current_request.get("preprocess") != merged_request.get("preprocess")
    )
    target_changed = current_request.get("target") != merged_request.get("target")
    changed_comparison_ids = _diff_comparisons(
        current_request.get("comparisons"), merged_request.get("comparisons"))

    if metadata_or_preprocess_changed or target_changed:
        for sid in stage_ids:
            if sid == "resolve_comparisons" or sid == "report" \
                    or sid.startswith(("differential:", "export:")):
                reset.add(sid)
    elif changed_comparison_ids:
        reset.add("resolve_comparisons")
        reset.add("report")
        for cid in changed_comparison_ids:
            reset.add(f"differential:{cid}")
            reset.add(f"export:{cid}")
    return reset


def _upstream_stage_id(record: dict) -> str:
    """このrunでConsoleを起動するstageのID（v1 `upstream` / v2 `execute_console`）。

    recordのschemaで決める——`prepare_resume`が上流の状態を見るのは要求本体を
    読み込む前で、そこで使える版情報はrecord側にしかない。
    """
    if record.get("schema") == store.SCHEMA_V2:
        return stage_plan.UPSTREAM_V2_STAGE_ID
    return "upstream"


def _reset_stage_for_resume(stage: dict) -> dict:
    reset = dict(stage)
    reset["status"] = "pending"
    reset["input_fingerprint"] = None
    reset["result_refs"] = []
    reset["error"] = None
    return reset


def _clear_cancel_request(pipeline_root: Path, *, prior_status: str,
                          worker_health: str) -> None:
    """協調取消フラグ（`control/cancel-request.json`）を取り下げる。

    これを消さないと、`planned`へ戻して新しいworkerを起こしても
    `engine._run_stage_loop`が最初のstage境界で`cancel_requested()`を見て即座に
    `finish_cancelled`するため、**取消したrunは二度と再開できない**（D09
    「取消後・timeout後のresumeは何も変えなくても続行できて当然」）。

    ただし**停止処理がまだ進行中かもしれない取消は握り潰さない**。
    `prior_status`が`running`/`planned`（＝workerが今このrunを進めている最中で
    ありうる状態）なら、そのworkerが**もう居ないと観測できたとき**
    （`worker_health == "worker_missing"`）以外は消さずに残す
    ——そこで消すと、利用者が出した取消がworkerに届く前に黙って無効になる。

    `"unknown"`（identityがまだ記録されていない／生存確認そのものが失敗した）を
    「居ない」と読まないのが要点で、`_RESUMABLE_RUNNING_HEALTH`が`"unknown"`を
    入れないのと同じ理由。`engine.run_engine`がworker identityを刻むのは
    owner lockを取ったあと——起動から数秒間、runは`planned`のままidentityが
    無く`"unknown"`になる。その窓でフラグを消すと、これからstage loopへ入る
    workerは取消を一度も見ずにMS-DIALを起動してしまう。

    握り潰さない代償は、resumeがもう1往復要ること（起きたworkerが最初の
    stage境界で`cancelled`を確定させ、そのあとのresumeでフラグが消える）だけ
    で、恒久的に再開不能にはならない。

    逆に`cancelled`（`engine.finish_cancelled`が停止完了を確定させた）や
    その他の終端状態から再開するときは、識別子として記録されているworkerが
    まだ生きていても（終了処理中・pid再利用）取消はもう完了しているので消す。
    """
    if prior_status in {"running", "planned"} \
            and worker_health not in _RESUMABLE_RUNNING_HEALTH:
        return
    engine.cancel_request_path(pipeline_root).unlink(missing_ok=True)


def _resume_receipt(record: dict, *, reused: bool, reset_stage_ids: list | None = None) -> dict:
    return {
        "pipeline_id": record["identity"]["pipeline_id"],
        "status": record["status"],
        "request_revision": record["request"]["revision"],
        "effective_target": record["request"].get("effective_target"),
        "reused": reused,
        "reset_stage_ids": sorted(reset_stage_ids or []),
    }


# ---------- prepare_resume ----------

def prepare_resume(path: Path, *, updates: dict | None = None,
                   request_id: str | None = None, rerun_upstream: bool = False) -> dict:
    """新しいrevision/attemptと再利用判断を保存する。**launchしない**。

    手順は本モジュールdocstringのa〜eのとおり。`STATE_REVISION_CONFLICT`
    （他アクターがこのrunを同時に更新した）は有限回まで読み直して吸収する
    （`store._patch_request_id_with_retry`と同じ考え方）。
    """
    pipeline_root = normalize_pipeline_root(path)
    updates = dict(updates) if updates else {}
    updates_hash = canonical_hash({"updates": updates, "rerun_upstream": rerun_upstream})

    for _ in range(_RESUME_RETRY_MAX_ATTEMPTS):
        record = store.load_run(pipeline_root)

        # --- 冪等性: 同じrequest_idの再送は直前の結果をそのまま返す ---
        resume_log = dict(record.get("resume_log") or {})
        if request_id is not None and request_id in resume_log:
            prior = resume_log[request_id]
            if prior["updates_hash"] != updates_hash:
                raise DomainError(
                    "IDEMPOTENCY_CONFLICT",
                    f"request_id={request_id!r}は既存の異なる内容のresumeと衝突しています。",
                    {"request_id": request_id})
            return _resume_receipt(record, reused=True,
                                   reset_stage_ids=prior.get("reset_stage_ids"))

        # --- a. 監視を失った稼働中Console ---
        console_identity, verified_receipt = _console_supervision_state(record)
        if console_identity and same_process(console_identity) and not verified_receipt:
            raise DomainError("EXECUTION_UNRESOLVED", "監視を失ったConsoleが稼働しています")

        # --- b. 上流再実行の明示要求 ---
        upstream_stage_id = _upstream_stage_id(record)
        upstream_status = record["stages"].get(upstream_stage_id, {}).get("status")
        upstream_needs_rerun = upstream_status != "succeeded"
        if upstream_needs_rerun and not rerun_upstream:
            raise DomainError("UPSTREAM_RERUN_REQUIRED", "上流の再実行を明示してください")

        # --- c. 固定済み入力の再検証（上流をやり直さない場合だけ） ---
        if not rerun_upstream:
            _verify_staged_inputs(record)

        # --- d. updatesの検証・統合 ---
        current_request = _load_full_request(pipeline_root, record)
        merged_request = (
            request_mod.merge_updates(current_request, updates,
                                      **_profile_arguments(current_request))
            if updates else current_request)
        new_content_hash = request_mod.request_fingerprint(merged_request)
        request_changed = new_content_hash != record["request"]["content_hash"]

        record = copy.deepcopy(record)

        if request_changed:
            next_revision = int(record["request"]["revision"]) + 1
            saved_rel = f"{store.REQUESTS_SUBDIR}/revision-{next_revision:04d}.json"
            atomic_write_json(pipeline_root / saved_rel, merged_request)
            record["request"] = {
                "revision": next_revision,
                "request_id": record["request"].get("request_id"),
                "content_hash": new_content_hash,
                "saved_path": saved_rel,
                "effective_target": merged_request.get("effective_target"),
            }

        # --- e. 依存区分に基づくstageのpending化 ---
        for sid in store.stage_ids_for(merged_request):
            if sid not in record["stages"]:
                record["stages"][sid] = store.initial_stage(sid)

        manifest_changed = _manifest_content_changed(record, merged_request)
        reset_ids = _stages_to_reset(
            current_request, merged_request, rerun_upstream=rerun_upstream,
            stage_ids=set(record["stages"].keys()),
            manifest_content_changed=manifest_changed)
        for sid in reset_ids:
            record["stages"][sid] = _reset_stage_for_resume(record["stages"][sid])

        # 中断由来の終端状態（partial/failed/cancelled）とneeds_inputは、resumeが
        # 呼ばれた時点で常に`planned`へ戻す——「取消後・timeout後のresume」は
        # 何も変えなくても続行できて当然（D09）。一方`completed`は、実際に何か
        # 変わった場合（request内容・stageのpending化・上流再実行の明示）だけ
        # `planned`へ戻す——変化の無いno-op resumeでcompletedを崩さない。
        #
        # `running`のまま固まったrun（worker喪失）は、**そのworkerがもう居ない
        # ことを観測できたときだけ**同じ扱いにする（A06/§9.3。`read_status`が
        # `PIPELINE_INTERRUPTED`で案内する再開経路の実体）。判定不能
        # （`observed_health="unknown"`）は「死んでいる」ではないので動かさない。
        made_changes = (request_changed or bool(reset_ids) or rerun_upstream
                        or manifest_changed)
        worker_health = _worker_health(record)
        prior_status = record["status"]
        interrupted = (prior_status in _INTERRUPTED_STATUSES
                       or (prior_status == "running"
                           and worker_health in _RESUMABLE_RUNNING_HEALTH))
        if interrupted:
            record["status"] = "planned"
            record["needs_input"] = None
        elif prior_status == "completed" and made_changes:
            record["status"] = "planned"

        if request_id is not None:
            resume_log[request_id] = {
                "updates_hash": updates_hash,
                "revision": record["request"]["revision"],
                "reset_stage_ids": sorted(reset_ids),
            }
            record["resume_log"] = resume_log

        try:
            store.save_run(pipeline_root, record, expected_revision=record["state_revision"])
        except DomainError as exc:
            if exc.code == "STATE_REVISION_CONFLICT":
                continue  # 他アクターが進めた最新を読み直して再試行
            raise
        if record["status"] == "planned":
            _clear_cancel_request(pipeline_root, prior_status=prior_status,
                                  worker_health=worker_health)
        saved = store.load_run(pipeline_root)
        return _resume_receipt(saved, reused=False, reset_stage_ids=reset_ids)

    raise DomainError(
        "STATE_REVISION_CONFLICT",
        f"resumeの状態競合の再試行上限（{_RESUME_RETRY_MAX_ATTEMPTS}回）に達しました: "
        f"{pipeline_root}",
        {"pipeline_root": str(pipeline_root)})
