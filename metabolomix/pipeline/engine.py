"""pipeline-run.v3 の工程を1回分進めるエンジン（spec §9・plan 2026-11-09）。

**sessionを一切importしない。** `metabolomix.core.session_state` / `metabolomix.core.mcp_core` /
`metabolomix.tools.*` はここからimportしてはいけない（`tests/test_pipeline_engine.py` が
ASTで検査する）。解析の進行状況は `pipeline-run.json`（`metabolomix.pipeline.store`）に住む。

**工程は上流だけ。** 極性の単位ごとに `prepare_input:<p>` → `upstream:<p>` →
`validate_outputs:<p>` を要求の順に回し、最後に単位を持たない `receipt` を 1 回だけ回す
（`build_stages` が唯一の正準）。データセットを開く・前処理・比較は対話の側
（`load_dataset(pipeline_path=…)` → `dataset_preprocess` / `dataset_compare` / `dataset_pca`）。

**単位の失敗は他の単位を止めない。** ある単位の工程が `failed` になったら、その単位の
残りの工程は `pending` のまま止め、次の単位へ進む（spec §5「片方が失敗しても、もう片方の
結果は残す」）。失敗した工程を engine が自分で選び直すことは無い——Console の自動再試行は
しない。やり直しは `recovery.prepare_resume(rerun_upstream=True)` が明示的に `pending` へ
戻したときだけ。`receipt` は検証を通った単位が 1 つでもあれば走る。

**stage_inputs_unchanged は2種類に分ける（handlerキーで判定）。**

1. `prepare_input` / `upstream` / `validate_outputs`（`_TRUST_PERSISTED_HANDLERS`）:
   永続状態が succeeded/skipped ならそのまま `True`（信頼して skip）。hash の再検証は
   ここでは行わない——`upstream` をこの関数の都合で選び直すと、Console を黙って
   再起動しかねない。壊れ・改変の検出は `recovery.prepare_resume` の責務。
2. `receipt`: 自身の `result_refs` を `store.verify_result_refs` で実ファイルと照合し、
   一致すれば `True`、消えた・改変されていれば `False`（Console を起動しない工程なので、
   自動的にやり直しても安全）。

**例外はStageResultの中でだけ表現する。** handler呼出しは必ずtry/exceptで囲み、
`DomainError` も他の例外も traceback をログへ残した上で `failed` の `StageResult` に
変換する——「例外が成功結果に化ける」ことは絶対にない。worker の中の `needs_input` は
無い（受付の事前検査だけが `needs_input` を保存する）。

**owner lockとstate更新lockは別ファイル。** `metabolomix.core.process_control.file_lock`
は同一ファイルへの新しい `os.open` ごとに独立したロック要求として扱われ、
同一プロセス内の入れ子取得でも自己デッドロックする。そこで:

- **owner lock**（`control/worker.lock`）: `run_engine` がstage loop全体の間
  ずっと保持する「このpipelineを今実行しているのは自分だけだ」という主張。
  取得に失敗したら即座に `PIPELINE_ALREADY_RUNNING`。
- **state更新lock**（`control/pipeline.lock`、`store.save_run` が呼出しごとに取る）:
  stage境界の1回の保存だけを直列化する短いロック。
"""
from __future__ import annotations

import copy
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.process_control import file_lock, process_identity
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import store

__all__ = [
    "RECEIPT_STAGE_ID",
    "build_stages",
    "cancel_request_path",
    "cancel_requested",
    "commit_stage_outcome",
    "finish_cancelled",
    "final_status",
    "finish_success",
    "make_context",
    "mark_stage_running",
    "run_engine",
    "stage_inputs_unchanged",
]

_logger = logging.getLogger(__name__)

#: 単位を持たない最後の工程（受領）。
RECEIPT_STAGE_ID = "receipt"

#: owner lockの取得タイムアウト（秒）。長時間待つ理由がない——待っても取れない
#: ということは、本当に別workerがこのpipelineを実行中だということ。
_OWNER_LOCK_TIMEOUT_S = 1.0
_OWNER_LOCK_FILENAME = "worker.lock"

#: 取消の協調フラグ。`recovery.request_cancel` が書き、ここではその置き場所の契約と
#: 読み手だけを持つ。owner lockの解放を待たずに書ける小さな独立ファイルにする。
_CANCEL_REQUEST_FILENAME = "cancel-request.json"


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _owner_lock_path(pipeline_path: Path) -> Path:
    return Path(pipeline_path) / store.CONTROL_SUBDIR / _OWNER_LOCK_FILENAME


def cancel_request_path(pipeline_path: Path) -> Path:
    """取消フラグの置き場所（`request_cancel` と読み手を1か所で揃える）。"""
    return Path(pipeline_path) / store.CONTROL_SUBDIR / _CANCEL_REQUEST_FILENAME


def cancel_requested(pipeline_path: Path) -> bool:
    """協調的な取消要求が保存済みかを返す（owner lockを一切取らない読取り専用）。"""
    path = cancel_request_path(pipeline_path)
    try:
        raw = path.read_text(encoding="utf-8")
    except OSError:
        return False
    try:
        data = json.loads(raw)
    except ValueError:
        return False
    return bool(data.get("cancel_requested"))


# ---------- stage計画 ----------

def build_stages(request: dict) -> list[dict]:
    """要求から、このrunが実行すべきstageの計画を組み立てる（実行順の唯一の正準）。

    各要素は `{"stage_id", "handler", "unit"}`。単位（要求の `units` の順）ごとに
    `prepare_input:<p>` / `upstream:<p>` / `validate_outputs:<p>`、最後に
    `{"stage_id": "receipt", "handler": "receipt", "unit": None}`。`store.create_run` が
    `record["stages"]` を作るときも、この関数の stage_id を使う。
    """
    plan: list[dict] = []
    for unit in request.get("units") or []:
        polarity = unit.get("polarity")
        if not polarity:
            raise DomainError(
                "PIPELINE_REQUEST_INVALID",
                "極性の確定していない単位からは工程を組めません。",
                {"unit": unit})
        for handler in store.UNIT_HANDLERS:
            plan.append({"stage_id": store.unit_stage_id(handler, polarity),
                         "handler": handler, "unit": polarity})
    plan.append({"stage_id": RECEIPT_STAGE_ID, "handler": RECEIPT_STAGE_ID, "unit": None})
    return plan


# ---------- ループ内helper ----------

#: 永続状態がsucceeded/skippedのままなら無条件に信頼してskipする工程（handlerキー）。
#: `upstream`のConsole再起動は`rerun_upstream=True`を通じてrecoveryが明示的に
#: stageを`pending`へ戻さない限り、このengineが自発的に選び直すことは絶対にない。
_TRUST_PERSISTED_HANDLERS = frozenset(store.UNIT_HANDLERS)


def _handler_of(stage_id: str) -> str:
    return str(stage_id).split(":", 1)[0]


def stage_inputs_unchanged(stage: dict, *, pipeline_root: Path) -> bool:
    """既存の成功/skip結果をそのまま使ってよいか（＝handlerを呼ばず飛ばせるか）。

    呼び出し元は、対象stageの永続statusが既に `succeeded`/`skipped` である場合にしか
    本関数を呼ばない。分類はモジュール docstring の2分類のとおり。
    """
    if _handler_of(stage.get("stage_id")) in _TRUST_PERSISTED_HANDLERS:
        return True
    failures = store.verify_result_refs(Path(pipeline_root), stage.get("result_refs") or [])
    return not failures


def make_context(record: dict, stage: dict, entry: dict, request: dict) -> dict:
    """handlerへ渡すstage contextを組み立てる。

    `request` は保存済みの完全な要求（v3 の正規形）、`request_meta` は
    `record["request"]` の要約（どの版か）。単位の工程には `unit`（極性）・
    `unit_index`（要求の `units` の何番目か）・`unit_root`（`<pipeline_root>/<p>`）・
    `unit_state`（`record["units"][p]` の写し）・`unit_request`（共有のキーと単位のキーを
    平らにした、入力計画が読む形）も渡す。`receipt` ではこれらは None。
    """
    pipeline_root = Path(record["identity"]["pipeline_root"])
    unit = entry.get("unit")
    context = {
        "pipeline_root": pipeline_root,
        "pipeline_id": record["identity"]["pipeline_id"],
        "identity": dict(record["identity"]),
        "stage_id": stage["stage_id"],
        "attempt": stage.get("attempt", 0),
        "request": copy.deepcopy(request),
        "request_meta": dict(record["request"]),
        "results": list(record.get("results") or []),
        "unit": unit,
        "unit_index": None,
        "unit_root": None,
        "unit_state": None,
        "unit_request": None,
    }
    if unit is not None:
        index = request_mod.unit_index(request, unit)
        unit_state = copy.deepcopy((record.get("units") or {}).get(unit) or {})
        context.update({
            "unit_index": index,
            "unit_root": pipeline_root / (unit_state.get("unit_root") or unit),
            "unit_state": unit_state,
            "unit_request": request_mod.unit_request(
                request, index, Path(record["identity"]["source_root"])),
        })
    return context


#: engineの状態保存が`STATE_REVISION_CONFLICT`を読み直して再適用する上限。
#: `store._REQUEST_ID_PATCH_MAX_ATTEMPTS`・`recovery._RESUME_RETRY_MAX_ATTEMPTS`と
#: 同じ考え方（有限回で十分吸収できる。無関係な同時更新でworkerを殺さない）。
_STATE_SAVE_MAX_ATTEMPTS = 5

#: 要求版の変化を検知してstage計画を組み直す上限（`_run_stage_loop`）。
#: ここに達するのは、要求の更新が途切れずに届き続けてworkerが1passも完走できない
#: 状態だけで、そのときは黙って回り続けるより止めて知らせるほうが正しい。
_REQUEST_RELOAD_MAX_PASSES = 20


def _save_with_retry(pipeline_path: Path, record: dict, apply_change) -> dict:
    """`apply_change(record) -> record`を保存し、状態競合は読み直して再適用する。

    長時間のstage（upstream）の最中に`pipeline_resume`が届くと——`prepare_resume`は
    workerが生きている限りstatusを動かさないが`state_revision`は必ず1つ進める——
    stage開始時に読んだrecordでの保存が`STATE_REVISION_CONFLICT`になる。吸収しないと
    workerが死に、既に成功していたConsole実行を捨てさせる。

    `apply_change`は「最新recordへの純粋な再適用」でなければならない
    （変更内容をrecordの現在値から計算し、`state_revision`には触れない）。
    """
    for _ in range(_STATE_SAVE_MAX_ATTEMPTS):
        updated = apply_change(copy.deepcopy(record))
        try:
            store.save_run(pipeline_path, updated,
                          expected_revision=updated["state_revision"])
        except DomainError as exc:
            if exc.code == "STATE_REVISION_CONFLICT":
                record = store.load_run(pipeline_path)  # 最新を読み直して再適用
                continue
            raise
        return store.load_run(pipeline_path)
    raise DomainError(
        "STATE_REVISION_CONFLICT",
        f"工程状態の保存が状態競合の再試行上限（{_STATE_SAVE_MAX_ATTEMPTS}回）に"
        f"達しました: {pipeline_path}",
        {"pipeline_root": str(pipeline_path)})


class _StaleRequest(Exception):
    """running にする直前に要求の版が変わっていた（pass をやり直す合図。外へは出ない）。"""


def mark_stage_running(pipeline_path: Path, record: dict, stage: dict, *,
                       identity: tuple | None = None) -> dict:
    """stageを`running`にして原子的に保存し、再読込した最新recordを返す。

    状態競合は`_save_with_retry`が読み直して再適用する——`attempt`は毎回
    「そのとき読んだrecordの値+1」なので、再適用しても二重加算にならない。
    `identity`（pass が実行している要求の版）を渡すと、保存する記録の版がそれと違う
    とき running にせず `_StaleRequest` を送出する——版の確認と保存の間に届いた訂正を
    古い要求のまま回し、その結果を以後信頼してしまわないため。
    """
    stage_id = stage["stage_id"]

    def _apply(current: dict) -> dict:
        if identity is not None and _request_identity(current) != identity:
            raise _StaleRequest()
        target = current["stages"][stage_id]
        now = _now_iso()
        target["status"] = "running"
        target["attempt"] = int(target.get("attempt") or 0) + 1
        target["started_at"] = now
        target["updated_at"] = now
        target["error"] = None
        current["status"] = "running"
        return current

    return _save_with_retry(pipeline_path, record, _apply)


def commit_stage_outcome(pipeline_path: Path, record: dict, stage: dict, outcome: dict, *,
                         unit: str | None = None) -> dict:
    """StageResultをstageと`results`へ反映して原子的に保存する。

    `outcome["result_refs"]` は `record["results"]`（append-only）へ追記する。ただし
    そのstageの直前の`result_refs`と**完全に同じ内容**なら追記しない（同じ成果物を
    指す outcome で履歴を無限に膨らませない）。

    handler は自分で `store.save_run` を呼ばない（ここでの1回の保存と競合するため）。
    record の更新は戻り値で運ぶ:
      - `record_updates`: `record` のトップレベルキーを丸ごと置き換える。
      - `unit_updates`: `record["units"][unit]` のキーごとの置き換え（単位の
        `inputs`・`upstream`・`outputs`）。`unit` はこの工程の単位。

    状態競合は`_save_with_retry`が読み直して再適用する。

    handler の実行中に `prepare_resume` がこの工程を `pending` へ戻していたら（訂正・
    `rerun_upstream`）、outcome は反映しない——戻した側が新しい条件でのやり直しを
    求めているので、古い条件の結果で `succeeded` に上書きしない（次の pass でやり直す）。
    """
    stage_id = stage["stage_id"]

    def _apply(current: dict) -> dict:
        target = current["stages"][stage_id]
        if target.get("status") == "pending":
            return current
        previous_result_refs = list(target.get("result_refs") or [])
        target["status"] = outcome.get("status")
        target["updated_at"] = _now_iso()
        target["result_refs"] = list(outcome.get("result_refs") or [])
        target["warnings"] = list(outcome.get("warnings") or [])
        target["error"] = outcome.get("error")

        result_refs = outcome.get("result_refs") or []
        if result_refs and result_refs != previous_result_refs:
            current["results"] = list(current.get("results") or []) + list(result_refs)

        for warning in outcome.get("warnings") or []:
            entry = dict(warning)
            entry.setdefault("stage_id", stage_id)
            current.setdefault("warnings", [])
            current["warnings"].append(entry)

        for key, value in (outcome.get("record_updates") or {}).items():
            current[key] = copy.deepcopy(value)
        unit_updates = outcome.get("unit_updates") or {}
        if unit_updates:
            if unit is None or unit not in (current.get("units") or {}):
                raise DomainError(
                    "PIPELINE_RUN_INVALID",
                    f"単位の無い工程が単位の更新を返しました: {stage_id}",
                    {"stage_id": stage_id, "unit": unit})
            for key, value in unit_updates.items():
                current["units"][unit][key] = copy.deepcopy(value)
        return current

    return _save_with_retry(pipeline_path, record, _apply)


def _unit_statuses(record: dict) -> dict[str, str]:
    return {p: store.unit_status(record, p) for p in store.run_units(record)}


def final_status(record: dict) -> str:
    """単位の状態と受領から run の最終 status を決める（取消は別）。

    全単位が validated・受領が成功・全単位の解析フォルダ（`.arf`）が mzTab-M と
    照合できた → `completed`。validated が 1 つ以上で、どれかが欠ける → `partial`。
    validated が 0 → `failed`。
    """
    statuses = _unit_statuses(record)
    validated = [p for p, s in statuses.items() if s == "validated"]
    if not validated:
        return "failed"
    units = record.get("units") or {}
    receipt_ok = ((record.get("stages") or {}).get(RECEIPT_STAGE_ID) or {}).get("status") \
        == "succeeded"
    aligned = all(((units.get(p) or {}).get("outputs") or {}).get("alignment")
                  for p in validated)
    if len(validated) == len(statuses) and receipt_ok and aligned:
        return "completed"
    return "partial"


def finish_success(pipeline_path: Path, record: dict) -> dict:
    """計画を最後まで通した（取消以外で止まらずに終えた）時点の最終判定。"""

    def _apply(current: dict) -> dict:
        current["status"] = final_status(current)
        current["needs_input"] = None
        return current

    return _save_with_retry(pipeline_path, record, _apply)


def finish_cancelled(pipeline_path: Path, record: dict) -> dict:
    def _apply(current: dict) -> dict:
        current["status"] = "cancelled"
        return current

    return _save_with_retry(pipeline_path, record, _apply)


# ---------- handler呼出し ----------

def _domain_error_payload(exc: DomainError) -> dict:
    return {"code": exc.code, "message": exc.message, "details": dict(exc.details)}


def _invoke_handler(handler: Callable[[dict], dict], context: dict) -> dict:
    """handlerを呼び、例外を成功結果へ絶対に変換しないStageResultへ落とす。"""
    try:
        outcome = handler(context)
    except DomainError as exc:
        _logger.exception("pipeline stage %s: DomainError", context.get("stage_id"))
        return {"status": "failed", "result_refs": [], "warnings": [],
                "error": _domain_error_payload(exc)}
    except Exception as exc:  # noqa: BLE001 - 未知例外を成功へ変換しないための意図的な全捕捉
        _logger.exception("pipeline stage %s: 未知の例外", context.get("stage_id"))
        return {"status": "failed", "result_refs": [], "warnings": [],
                "error": {"code": exc.__class__.__name__, "message": str(exc), "details": {}}}

    outcome = dict(outcome)
    outcome.setdefault("result_refs", [])
    outcome.setdefault("warnings", [])
    outcome.setdefault("error", None)
    return outcome


# ---------- 本体 ----------

def _load_request(pipeline_path: Path, record: dict) -> dict:
    """`record["request"]["saved_path"]` の完全な要求（v3 の正規形）を読む。"""
    saved_path = record["request"]["saved_path"]
    full_path = Path(pipeline_path) / saved_path
    return json.loads(full_path.read_text(encoding="utf-8"))


def _request_identity(record: dict) -> tuple:
    """「今どの要求版を実行しているか」を表す最小の識別子。

    `state_revision`は使わない——無関係な同時更新（取消フラグの保存など）でも進むため。
    """
    request_meta = record.get("request") or {}
    return (request_meta.get("revision"), request_meta.get("content_hash"),
            request_meta.get("saved_path"))


def _resume_reset_seen(record: dict, visited: dict[str, str]) -> bool:
    """このpassで既に通過したstageが、外から`pending`へ戻されたかを返す。

    `visited` は通過したときに見た status。最初から `pending` だったもの（失敗した単位の
    残りの工程）は差し戻しではない。`rerun_upstream=True` が失敗した工程を `pending` へ
    戻した場合は、要求の版が変わらなくても計画をやり直す（そうしないと稼働中の worker が
    戻された単位を見ないまま終わり、誰も再開しない）。
    """
    stages = record.get("stages") or {}
    return any((stages.get(sid) or {}).get("status") == "pending" and seen != "pending"
               for sid, seen in visited.items())


def _run_stage_loop(pipeline_path: Path, handlers: dict) -> dict:
    """要求版が変わっていないことを**stage境界ごとに**確かめながら1回分進める。

    workerが稼働している間に届いた`pipeline_resume`は、`pipeline-run.json`のrevisionと
    保存済み要求だけを差し替えて戻ってくる（稼働中workerがいる限り新しいworkerを
    起こさない）。そこで各stageの直前に最新recordを読み、要求版か、通過済みstageの
    差し戻しを観測したら、その場でstage計画を組み直してpassをやり直す。
    """
    record = store.load_run(pipeline_path)
    request = _load_request(pipeline_path, record)
    identity = _request_identity(record)

    for _ in range(_REQUEST_RELOAD_MAX_PASSES):
        visited: dict[str, str] = {}
        restarted = False

        for entry in build_stages(request):
            stage_id = entry["stage_id"]
            latest = store.load_run(pipeline_path)
            if _request_identity(latest) != identity or _resume_reset_seen(latest, visited):
                record = latest
                request = _load_request(pipeline_path, latest)
                identity = _request_identity(latest)
                restarted = True
                break
            record = latest
            stage = record["stages"].get(stage_id)
            if stage is None:
                raise DomainError(
                    "PIPELINE_RUN_INVALID",
                    f"計画の工程が run の記録にありません: {stage_id}",
                    {"pipeline_root": str(pipeline_path), "stage_id": stage_id})

            unit = entry["unit"]
            if unit is not None and store.unit_status(record, unit) == "failed":
                # 失敗した単位の残りの工程は pending のまま止める（自動でやり直さない）。
                visited[stage_id] = stage["status"]
                continue
            if unit is None and "validated" not in _unit_statuses(record).values():
                visited[stage_id] = stage["status"]
                continue  # 検証を通った単位が無ければ受領は書かない

            if (stage["status"] in {"succeeded", "skipped"}
                    and stage_inputs_unchanged(stage, pipeline_root=pipeline_path)):
                visited[stage_id] = stage["status"]
                continue

            if cancel_requested(pipeline_path):
                return finish_cancelled(pipeline_path, record)

            try:
                record = mark_stage_running(pipeline_path, record, stage, identity=identity)
            except _StaleRequest:
                record = store.load_run(pipeline_path)
                request = _load_request(pipeline_path, record)
                identity = _request_identity(record)
                restarted = True
                break
            stage = record["stages"][stage_id]
            context = make_context(record, stage, entry, request)
            outcome = _invoke_handler(handlers[entry["handler"]], context)
            record = commit_stage_outcome(pipeline_path, record, stage, outcome, unit=unit)
            visited[stage_id] = outcome.get("status")

            # 長い工程（Console）の最中の取消は、監視経路が Console を止めて handler が
            # failed を返す形で届く。`pipeline_cancel` の契約どおり cancelled で確定する。
            if outcome.get("status") == "failed" and cancel_requested(pipeline_path):
                return finish_cancelled(pipeline_path, record)

        if not restarted:
            # 最後の工程の最中に届いた訂正・差し戻しも、終える前に取り込む。
            latest = store.load_run(pipeline_path)
            if _request_identity(latest) != identity or _resume_reset_seen(latest, visited):
                record = latest
                request = _load_request(pipeline_path, latest)
                identity = _request_identity(latest)
                continue
            return finish_success(pipeline_path, latest)

    raise DomainError(
        "REQUEST_REVISION_CHURN",
        f"要求の更新が続いて工程計画が{_REQUEST_RELOAD_MAX_PASSES}回組み直されました。"
        f"更新を止めてから再開してください: {pipeline_path}",
        {"pipeline_root": str(pipeline_path)})


def _record_worker_identity(pipeline_path: Path) -> None:
    """自分（現在プロセス）のidentityをrunへ刻む（spec §9.3「ロックにはowner
    identityを記録し」の実体）。`recovery.read_status`/`prepare_resume`がこれを読み、
    `same_process`でこのworkerがまだ生きているかを判定する。
    """
    record = store.load_run(pipeline_path)

    def _apply(current: dict) -> dict:
        current["worker"] = {"identity": process_identity(os.getpid()), "started_at": _now_iso()}
        return current

    _save_with_retry(pipeline_path, record, _apply)


def run_engine(pipeline_path: Path, handlers: dict[str, Callable[[dict], dict]]) -> dict:
    """永続化済みrunを、計画の最後まで（または取消まで）1回分進める。

    `pipeline_path` は `store.create_run`/`find_or_create_run` が返す
    pipeline_root。戻り値は保存後に読み直した最新の `pipeline-run.json`。
    旧 v1 の run は `PIPELINE_RUN_LEGACY`（読むだけ）。

    owner lock（`control/worker.lock`）をstage loop全体の間ずっと保持する。
    取得できなければ即 `PIPELINE_ALREADY_RUNNING`。ロック取得の成否だけを判定するために
    `__enter__`/`__exit__` を手で呼ぶ（stage loop 本体の中の state 更新lockの
    `LOCK_TIMEOUT` まで `PIPELINE_ALREADY_RUNNING` に化けさせないため）。
    """
    pipeline_path = Path(pipeline_path)
    if store.is_legacy(store.load_run(pipeline_path)):
        raise store.legacy_refusal(pipeline_path, "実行")
    owner_lock = file_lock(_owner_lock_path(pipeline_path), timeout=_OWNER_LOCK_TIMEOUT_S)
    try:
        owner_lock.__enter__()
    except DomainError as exc:
        if exc.code == "LOCK_TIMEOUT":
            raise DomainError(
                "PIPELINE_ALREADY_RUNNING",
                f"このpipelineは既に別のworkerが実行中です: {pipeline_path}",
                {"pipeline_root": str(pipeline_path)},
            ) from exc
        raise
    try:
        _record_worker_identity(pipeline_path)
        return _run_stage_loop(pipeline_path, handlers)
    finally:
        owner_lock.__exit__(None, None, None)
