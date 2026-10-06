"""pipeline-run.v1 の工程を1回分進めるエンジン（spec §9）。

**sessionを一切importしない。** `metabolomix.core.session_state` / `metabolomix.core.mcp_core` /
`metabolomix.tools.*` はここからimportしてはいけない（`tests/test_pipeline_engine.py` が
ASTで検査する）。解析の進行状況は `pipeline-run.json`（Task14 `metabolomix.pipeline.store`）と、
このプロセスだけが持つ `runtime`（worker固有の `DatasetState` を積む素のdict）に住む。
`runtime` は絶対に永続化しない——再開後の新しいworkerプロセスは空の `runtime` から
始まる。したがって「前回成功した」という記録だけでstageを飛ばすと、後続stageが
参照する `runtime` の中身が無いまま呼ばれる事故になる。

**stage_inputs_unchanged（Task16で実装）は3種類に分ける。**

1. `load_dataset` / `resolve_metadata` / `preprocess` / `pca`
   （`_ALWAYS_RECONSTRUCT_STAGE_IDS`）: 常に `False`（絶対にskipしない）。
   `runtime` を実際に組み立てる工程だからで、再開のたびに——たとえ前回成功していても
   ——handlerを呼び直して `runtime` を作り直す（brief「loading→metadata→
   必要ならpreprocess/PCAを再構築する」）。これらはConsoleのような外部プロセスを
   起動しない、安価な再計算という前提に立つ。
2. `prepare_input` / `upstream` / `validate_outputs`
   （`_TRUST_PERSISTED_STAGE_IDS`）: 永続状態がsucceeded/skippedならそのまま `True`
   （信頼してskip）。**hashの再検証はここでは行わない**——`upstream`をこの関数の
   都合で自動的に選び直すと、壊れたreceiptを検出した瞬間にConsoleを黙って
   再起動しかねない（「自動再試行は行わない」に反する）。壊れ・改変の検出は
   `metabolomix.pipeline.recovery.prepare_resume` が resume の前段で行う責務とし、
   本当に再実行が要るときは `rerun_upstream=True` で明示的にこれらのstageを
   `pending` へ戻す（recoveryの責務）。エンジン自身は「まだ`succeeded`のまま
   残っているなら、それはrecoveryが既に安全と判断した結果」として信頼する。
3. それ以外（`resolve_comparisons` / `differential:*` / `export:*` / `report`）:
   自身の `result_refs` を `store.verify_result_refs` で実ファイルと照合し、
   一致すれば `True`、成果物が消えた・改変されていれば `False`（再実行——
   Consoleを起動しない工程なので、自動的にやり直しても安全）。

`stage_inputs_unchanged` 自体は「このstageは既にrecoveryの判断を経て
succeeded/skippedのままだ」という前提の上でしか呼ばれない
（`_run_stage_loop` が先にstatusを見てから渡す）。REQUESTの内容が変わった
かどうかの判定・どのstageを`pending`へ戻すかの決定は本モジュールの外
（`metabolomix.pipeline.recovery.prepare_resume`）にある。

**handlerキーとstage_idは別物。** `differential:<comparison_id>` / `export:<comparison_id>`
という複数のstage_idは、`differential` / `export` という**1つのhandlerキー**を共有し、
`comparison_id` はstage_idから機械的に取り出して context へ渡す（`build_stages` が
そのマッピングを組み立てる）。

**例外はStageResultの中でだけ表現する。** handler呼出しは必ずtry/exceptで囲み、
`DomainError` のうち「入力で解消する」codeだけを `needs_input` へ変換する
（`_NEEDS_INPUT_CODES`）。それ以外の `DomainError` および他の例外は
tracebackをログへ残した上で `failed` の `StageResult` に変換する——
どちらの経路でも「例外が成功結果に化ける」ことは絶対にない。

**owner lockとstate更新lockは別ファイル。** `metabolomix.core.process_control.file_lock`
は同一ファイルへの新しい `os.open` ごとに独立したロック要求として扱われ、
同一プロセス内の入れ子取得でも自己デッドロック（timeoutまで待って失敗）する
（`metabolomix/pipeline/store.py` 冒頭の説明と同じ注意）。そこで:

- **owner lock**（`control/worker.lock`）: `run_engine` がstage loop全体の間
  ずっと保持する「このpipelineを今実行しているのは自分だけだ」という主張。
  取得に失敗したら即座に `PIPELINE_ALREADY_RUNNING`（長時間待たない——
  待つ理由がない。取れなければ本当に別workerが動いている）。
- **state更新lock**（`control/pipeline.lock`、Task14 `store.save_run` が
  呼出しごとに取る）: stage境界の1回の保存だけを直列化する短いロック。
  ワーカーがowner lockを握っている間も、`pipeline_status`（読取専用）や
  取消要求の保存はこちらを使うだけで進められる——owner lockの解放を待たない。
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
from metabolomix.pipeline import report as report_mod
from metabolomix.pipeline import stage_plan
from metabolomix.pipeline import store

__all__ = [
    "build_stages",
    "cancel_request_path",
    "cancel_requested",
    "commit_stage_outcome",
    "finish_cancelled",
    "finish_interrupted",
    "finish_success",
    "make_context",
    "mark_stage_running",
    "run_engine",
    "stage_inputs_unchanged",
]

_logger = logging.getLogger(__name__)

#: brief記載の4コード＋レビュー指摘5の1コード。共通context「Task11/Task12」で
#: 名指しされた、入力を与え直せば先へ進める（＝ユーザー操作で解消しうる）
#: DomainErrorだけをここに置く。
#: PCA不成立（"an uncomputable PCA"）はbrief/spec文言上は同格に挙げられているが、
#: 現行実装（`metabolomix/analysis/dataset_analysis.py::run_dataset_pca`）は
#: `DomainError` ではなく別クラスの `PreconditionError` を送出しており、この
#: whitelistでは検出できない。**Task18のR18でこのwhitelistは変更していない**
#: ——`metabolomix.pipeline.service`の`pca` handler自身が`PreconditionError`を
#: 捕らえ、`needs_input`のStageResultを直接返すことで解決する（`_invoke_handler`
#: の汎用`except Exception`分岐に落ちて`failed`になるのを handler境界で防ぐ）。
_NEEDS_INPUT_CODES = frozenset({
    "PREPROCESS_PREREQUISITE_MISSING",  # Task11: 前提を欠く明示的な前処理要求
    "NORMALIZATION_DEGENERATE",          # Task11: 正規化係数が0/非有限の試料が残る
    "PREPROCESS_FEATURES_EXHAUSTED",     # Task11: 前処理後に特徴量が0件
    "COMPARISON_REQUIRED",               # Task12: 比較の対照群/比較群が未指定
    "CONFOUNDED_COMPARISON",             # Task12/レビュー指摘5: 群とバッチが完全交絡
                                         # （spec §7.4「探索結果を残してneeds_inputで
                                         # 入力を待つ」。allow_confounded=trueは
                                         # pipeline_resumeで与える利用者入力そのもの
                                         # ——R18がPCA不成立に適用したのと同じ理由で、
                                         # `resolve_comparisons`handlerと
                                         # `dataset_service.run_comparison`の両方の
                                         # 送出元をここ1箇所で救う）。
})

#: owner lockの取得タイムアウト（秒）。長時間待つ理由がない——待っても取れない
#: ということは、本当に別workerがこのpipelineを実行中だということ。
_OWNER_LOCK_TIMEOUT_S = 1.0
_OWNER_LOCK_FILENAME = "worker.lock"

#: 取消の協調フラグ。`request_cancel`（Task16 `metabolomix.pipeline.recovery`）が
#: 書き、ここではその置き場所の契約と読み手だけを持つ。state更新lockと同じく、
#: owner lockの解放を待たずに書ける小さな独立ファイルにする。
_CANCEL_REQUEST_FILENAME = "cancel-request.json"


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _owner_lock_path(pipeline_path: Path) -> Path:
    return Path(pipeline_path) / store.CONTROL_SUBDIR / _OWNER_LOCK_FILENAME


def cancel_request_path(pipeline_path: Path) -> Path:
    """取消フラグの置き場所（Task16 `request_cancel` と読み手を1か所で揃える）。"""
    return Path(pipeline_path) / store.CONTROL_SUBDIR / _CANCEL_REQUEST_FILENAME


def cancel_requested(pipeline_path: Path) -> bool:
    """協調的な取消要求が保存済みかを返す（owner lockを一切取らない読取り専用）。

    `atomic_write_json` は置換が原子的なので、ロック無しで読んでも「書きかけの
    中途半端な内容」を見ることはない（見えるのは直前の内容か、書き終わった
    新しい内容のどちらか）。
    """
    path = cancel_request_path(pipeline_path)
    try:
        raw = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return False
    except OSError:
        return False
    try:
        data = json.loads(raw)
    except ValueError:
        return False
    return bool(data.get("cancel_requested"))


# ---------- stage計画 ----------

def _comparison_id_from_stage_id(stage_id: str) -> str | None:
    for prefix in ("differential:", "export:"):
        if stage_id.startswith(prefix):
            return stage_id[len(prefix):]
    return None


def build_stages(request: dict) -> list[dict]:
    """要求から、このrunが実行すべきstageの計画を組み立てる。

    `pipeline-request.v2`（メタボロミクス）は`metabolomix.pipeline.stage_plan.build_v2`
    へ丸ごと委譲する——**v2のstage列の正準はあちら1箇所だけ**で、`store`も同じ
    builderを呼ぶ（順序の写しを2つ持たない）。以下はv1の規則。

    各要素は `{"stage_id", "handler", "comparison_id"}`。`differential:<id>` /
    `export:<id>` はどちらも1件のcomparisonに対応するstageで、handlerキーは
    それぞれ `"differential"` / `"export"` の1つだけを共有する
    （`comparison_id` で呼出し先の区別を渡す）。

    `effective_target == "exploratory"` のときは `resolve_comparisons` を計画から
    除く（spec §9.2「比較群を要求しない」）——`store.create_run` が作る
    `record["stages"]` にはtargetを問わず常にこのstage_idが存在するため
    （`store._stage_ids_for` はtargetを見ない）、計画に無いstage_idは
    `run_engine` 側で「対象外」として扱う（handlerを呼ばずskip）。

    stage_id集合そのものの計算（`_BASE_STAGE_IDS` の並び・
    `differential:`/`export:` の追加条件）は `store._stage_ids_for` と
    意図的に同じ規則にしてある——両者がズレると `record["stages"]` に無い
    stage_idを計画してしまい、`run_engine` が `KeyError` で落ちる。
    """
    if stage_plan.is_v2_request(request):
        return stage_plan.build_v2(request)

    effective_target = request.get("effective_target") or request.get("target")
    comparisons = request.get("comparisons") or []

    plan: list[dict] = []
    for stage_id in store.STAGE_IDS[:-1]:  # "report"を除いた基本8種
        if stage_id == "resolve_comparisons" and effective_target == "exploratory":
            continue
        plan.append({"stage_id": stage_id, "handler": stage_id, "comparison_id": None})

    if effective_target != "exploratory":
        for comparison in comparisons:
            cid = comparison["comparison_id"]
            plan.append({"stage_id": f"differential:{cid}", "handler": "differential",
                        "comparison_id": cid})
            plan.append({"stage_id": f"export:{cid}", "handler": "export",
                        "comparison_id": cid})

    final_stage_id = store.STAGE_IDS[-1]  # "report"
    plan.append({"stage_id": final_stage_id, "handler": final_stage_id, "comparison_id": None})
    return plan


# ---------- ループ内helper（brief step3） ----------

#: `runtime`を実際に組み立てる工程。再開のたびに必ずhandlerを呼び直す
#: （R20/spec §9.3「下流の再開は…DatasetStateと必要結果を再構築できることを必須とする」）。
_ALWAYS_RECONSTRUCT_STAGE_IDS = frozenset({
    "load_dataset", "resolve_metadata", "preprocess", "pca",
    # v2（spec §6.2）の解析鎖。成果物JSONには行列の**値**を入れていない
    # （数十MBになるうえ、記録の意味は「どの入力から出たか」であって数値の
    # 保管庫ではない）ので、runtimeは再計算でしか復元できない。skipすると
    # 下流のstatisticsが「行列が無い」で止まる。どれも決定論的で、Consoleは
    # 起動しない。
    "load_assay_evidence", "resolve_feature_bindings", "qc_raw", "qc_processed",
})

#: 永続状態がsucceeded/skippedのままなら無条件に信頼してskipするstage。
#: `upstream`のConsole再起動は`rerun_upstream=True`を通じてrecoveryが明示的に
#: stageを`pending`へ戻さない限り、このengineが自発的に選び直すことは絶対にない
#: （「自動再試行は行わない」）。
_TRUST_PERSISTED_STAGE_IDS = frozenset({
    "prepare_input", "upstream", "validate_outputs",
    # v2（spec §6.2）の同じ位置の工程。`validate_outputs`は綴りが同じなので上の
    # 行が兼ねる。ここに載せないと`execute_console`が成果物hashの照合で
    # 「変わった」と読まれ、engineが自分の判断でMS-DIALを起動し直しうる
    # ——「Consoleの自動再試行は行わない」に反する。
    "prepare_inputs", "execute_console",
})


def stage_inputs_unchanged(stage: dict, *, pipeline_root: Path) -> bool:
    """既存の成功/skip結果をそのまま使ってよいか（＝handlerを呼ばず飛ばせるか）。

    呼び出し元（`_run_stage_loop`）は、対象stageの永続statusが既に
    `succeeded`/`skipped`である場合にしか本関数を呼ばない。つまりここでの
    「True」は「まだ信頼してよい」、「False」は「（安全に）handlerを呼び直す」
    という意味であり、要求内容が変わったかどうかの判定・どのstageを`pending`へ
    戻すかの決定は本関数の外（`metabolomix.pipeline.recovery.prepare_resume`）が
    既に済ませている前提に立つ。

    分類はモジュール docstring の3分類のとおり:
      - `_ALWAYS_RECONSTRUCT_STAGE_IDS`: 常に`False`（`runtime`復元のため）。
      - `_TRUST_PERSISTED_STAGE_IDS`: 常に`True`（`upstream`の自動再起動を防ぐ）。
      - それ以外: 自身の`result_refs`のhashを実ファイルと照合し、一致すれば
        `True`、成果物が消失・改変されていれば`False`（安全に再実行できる工程
        のみがここに属する——Consoleを起動しない）。
    """
    stage_id = stage.get("stage_id")
    if stage_id in _ALWAYS_RECONSTRUCT_STAGE_IDS:
        return False
    if stage_id in _TRUST_PERSISTED_STAGE_IDS:
        return True
    failures = store.verify_result_refs(Path(pipeline_root), stage.get("result_refs") or [])
    return not failures


def make_context(record: dict, stage: dict, runtime: dict, request: dict) -> dict:
    """handlerへ渡すstage contextを組み立てる（R19: run_record全体を見渡せるよう拡張）。

    Task15時点では`request`（`record["request"]`の要約——revision/request_id/
    content_hash/saved_path/effective_targetだけ）しか渡しておらず、handlerが
    `target`/`comparisons`/`preprocess`/`sample_manifest`本体へ触るには
    `store.load_run`を自分で呼び直す必要があった。ここでは呼び出し元
    （`_run_stage_loop`。既に`_load_request`で読み込み済みの完全な
    pipeline-request.v1）から丸ごと受け取り、`context["request"]`として渡す
    ——`record["request"]`の要約は`context["request_meta"]`として別に残す
    （revision/saved_pathなど、要求本体ではなく「この要求のどの版か」を指す
    情報はこちらにしかない）。

    合わせて`identity`（pipeline_id/source_root/pipeline_root/作成・更新時刻）・
    `inputs`（Task13が固定した入力スナップショット）・`upstream`（Console job
    path・execution_id・検証記録）も渡す——いずれも`record`が最初から持つ
    フィールドで、handlerが`store.load_run`を再呼出ししなくても読めるようにする。

    `runtime` はこのworkerプロセスだけが持つ素のdict（例: `runtime["dataset"]`に
    `DatasetState` を積む）への参照そのものを渡す——handler間で状態を共有するのは
    同一worker内のこの1呼出し系列だけ、という契約をそのまま体現する。
    """
    stage_id = stage["stage_id"]
    # `export:<id>`はv1（comparison）とv2（statistic）で綴りが重なる。要求のschemaで
    # どちらの意味かを決め、両方に値が入った曖昧なcontextを作らない。
    is_v2 = stage_plan.is_v2_request(request)
    return {
        "pipeline_root": Path(record["identity"]["pipeline_root"]),
        "pipeline_id": record["identity"]["pipeline_id"],
        "identity": dict(record["identity"]),
        "stage_id": stage_id,
        "comparison_id": None if is_v2 else _comparison_id_from_stage_id(stage_id),
        "statistic_id": stage_plan.statistic_id_from_stage_id(stage_id) if is_v2 else None,
        "attempt": stage.get("attempt", 0),
        "request": dict(request),
        "request_meta": dict(record["request"]),
        "inputs": copy.deepcopy(record.get("inputs") or {}),
        "upstream": copy.deepcopy(record.get("upstream") or {}),
        "results": list(record.get("results") or []),
        "runtime": runtime,
    }


#: engineの状態保存が`STATE_REVISION_CONFLICT`を読み直して再適用する上限。
#: `store._REQUEST_ID_PATCH_MAX_ATTEMPTS`・`recovery._RESUME_RETRY_MAX_ATTEMPTS`と
#: 同じ考え方（有限回で十分吸収できる。無関係な同時更新でworkerを殺さない）。
_STATE_SAVE_MAX_ATTEMPTS = 5

#: 要求版の変化を検知してstage計画を組み直す上限（`_run_stage_loop`）。
#: 実運用で人が訂正を重ねても数回で収まる。ここに達するのは、要求の更新が
#: 途切れずに届き続けてworkerが1passも完走できない状態だけで、そのときは
#: 黙って回り続けるより止めて知らせるほうが正しい。
_REQUEST_RELOAD_MAX_PASSES = 20


def _save_with_retry(pipeline_path: Path, record: dict, apply_change) -> dict:
    """`apply_change(record) -> record`を保存し、状態競合は読み直して再適用する。

    engineの保存だけが、他の状態書込み側（`store._patch_request_id_with_retry`・
    `recovery.prepare_resume`はどちらも有限回のリトライを持つ）と違って
    競合を吸収していなかった。長時間のstage（upstream）の最中に
    `pipeline_resume`が届くと——`prepare_resume`はworkerが生きている限り
    statusを動かさないが`state_revision`は必ず1つ進める——stage開始時に
    読んだrecordでの保存が`STATE_REVISION_CONFLICT`になり、その送出は
    `_invoke_handler`のtryの**外**なので`run_engine`を貫いてworkerが死ぬ。
    ディスク上はstageが`running`のまま残り、次のresumeは
    `UPSTREAM_RERUN_REQUIRED`——既に成功していたConsole実行を捨てさせる。

    `apply_change`は「最新recordへの純粋な再適用」でなければならない
    （変更内容をrecordの現在値から計算し、`state_revision`には触れない）。
    そうであれば、読み直してもう一度適用するだけで正しい結果になる。
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


def mark_stage_running(pipeline_path: Path, record: dict, stage: dict) -> dict:
    """stageを`running`にして原子的に保存し、再読込した最新recordを返す。

    毎回 `save_run` 後に `load_run` で読み直す（`save_run` は渡した引数の
    dictを書き換えず、内部で深複製した別dictへ`state_revision`/`updated_at`を
    刻んで書くため、呼び出し側が手で追随するより読み直す方が確実）。

    状態競合は`_save_with_retry`が読み直して再適用する——`attempt`は毎回
    「そのとき読んだrecordの値+1」なので、再適用しても二重加算にならない。
    """
    stage_id = stage["stage_id"]

    def _apply(current: dict) -> dict:
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


def commit_stage_outcome(pipeline_path: Path, record: dict, stage: dict, outcome: dict) -> dict:
    """StageResultをstageと`results`へ反映して原子的に保存する。

    `outcome["result_refs"]` は無条件に `record["results"]`（append-only、
    `store.save_run` が既存要素の書換えを拒否する）へ追記する——failed/needs_input
    でも、そのstageが実際に有効な出力を残していたなら「有効な解析出力があるか」
    （`finish_interrupted` のpartial/failed判定）の材料になる。

    ただし、そのstageの直前の`result_refs`と**完全に同じ内容**を返した場合は
    追記しない（Task16: `_ALWAYS_RECONSTRUCT_STAGE_IDS` は再開のたびに
    handlerを呼び直すため、決定論的なhandlerが前回と同じ成果物を指すoutcomeを
    返すと、何もしなければ`results`が無限に重複して膨らむ。「旧結果を書き換え
    ない」はここでは満たしたまま——新しい内容の追記だけが対象で、既存要素は
    一切触らない）。

    `outcome["record_updates"]`（Task18追加、既定なしで後方互換）は
    `record`のトップレベルキーを丸ごと置き換える唯一の経路。`prepare_input`
    handlerが`record["inputs"]`を実配置後のsnapshotへ、`upstream` handlerが
    `record["upstream"]`をConsole終了証跡の参照へ更新するのに使う。handlerが
    自分で`store.save_run`を呼ばないのは、ここでの1回の保存と競合し
    `STATE_REVISION_CONFLICT`になるため（`record`は`mark_stage_running`が
    返した、handler呼出し**前**のスナップショット）——`record_updates`は
    handlerの戻り値として運び、ここで初めて`record`へ反映してから1回だけ
    保存する。

    状態競合は`_save_with_retry`が読み直して再適用する。反映内容は
    `outcome`と「そのとき読んだrecord」だけから計算するため、通常は
    再適用してもresultsが二重に積まれることはない（`previous_result_refs`
    と一致すれば追記しない、上記の判定による）。

    ただし一つ例外がある——再適用のために読み直したstageが、その間に
    並行する`prepare_resume`から**同じstage**のリセット（`result_refs`を
    `[]`にする）を受けていた場合、`previous_result_refs`は`[]`に見えて
    しまい、既に`results`へ入っている同じrefsをもう一度追記しうる。
    `store._assert_results_append_only`は追記であることに変わりはないため
    これを拒否せず、`report._achieved_ref`は一致する最後の要素を採るため
    現時点で実害はないが、履歴には重複が残る（詳細は
    `docs/superpowers/notes/2026-09-05-raw-folder-pipeline-validation.md`の
    既知の制限「リトライ経路でresult_refsが重複しうる」を参照）。
    """
    stage_id = stage["stage_id"]

    def _apply(current: dict) -> dict:
        target = current["stages"][stage_id]
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

        record_updates = outcome.get("record_updates")
        if record_updates:
            for key, value in record_updates.items():
                current[key] = copy.deepcopy(value)
        return current

    return _save_with_retry(pipeline_path, record, _apply)


def _mark_out_of_scope(pipeline_path: Path, record: dict, stage_id: str) -> dict:
    """計画に無いstage（exploratory目標の`resolve_comparisons`）を対象外として記録する。

    「stageが無い＝completedへ数える」を避けるため、`skipped`のまま明示的に
    理由をwarningsへ残す（spec「stageがないからcompletedとしない」の裏返しで、
    こちらは「対象外だからhandlerを呼ばずに進める」ことを可視化する）。
    """
    record = copy.deepcopy(record)
    stage = record["stages"][stage_id]
    if stage["status"] == "skipped":
        return record  # 既に反映済み（run_engineを繰り返し呼んでも増殖させない）
    stage["status"] = "skipped"
    stage["updated_at"] = _now_iso()
    stage["warnings"] = list(stage.get("warnings") or []) + [{
        "code": "STAGE_OUT_OF_SCOPE_FOR_TARGET",
        "stage_id": stage_id,
        "message": f"effective_targetにより対象外のためskipしました: {stage_id}",
    }]
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return store.load_run(pipeline_path)


def finish_success(pipeline_path: Path, record: dict) -> dict:
    """全stageがsucceeded/skippedで止まらずに終えた時点の最終判定（Task18）。

    Task17時点はここを`record["warnings"]`の既知codeだけで`partial`へ
    格下げする狭い判定（`_PARTIAL_ON_WARNING_CODES`、消去済み）に留めていた。
    `metabolomix.pipeline.report.evaluate_target`が持つ、目標別必須出力の
    hash/ID照合まで含む厳密な判定と統合していなかったため、「全stageが
    succeededでも必須出力が足りない」ケース（例: exportがwarningだけ残して
    result_refsを空で返した場合）を狭い判定でしか捉えられなかった。

    build_handlersが各stageの`result_refs`へ`output_name`付きの正しいrefを
    積むようになった（Task18）ので、ここで`evaluate_target(record)`を1回
    呼び、その`status`をそのまま採用する——`_PARTIAL_ON_WARNING_CODES`が
    捉えていた「EXPORT_BACKGROUND_EMPTY（tsv:<cid>が無い）」
    「GUI_PROJECT_UNAVAILABLE（gui_projectが無い）」は、どちらも
    `evaluate_target`の`missing_outputs`判定に自然に含まれるため、狭い
    判定を残す必要が無くなった。

    `evaluate_target`が"needs_input"を返すのは`target=differential`かつ
    comparisonsが空の場合だけだが、その組合せは実handler
    （`resolve_comparisons`）がCOMPARISON_REQUIREDを`_NEEDS_INPUT_CODES`経由で
    先に検出し`finish_interrupted`へ抜けるため、正常経路ではここへ到達しない
    ——防御的に残すだけの分岐。
    """
    record = copy.deepcopy(record)
    evaluation = report_mod.evaluate_target(record)
    if evaluation["status"] == "completed":
        record["status"] = "completed"
    elif evaluation["status"] == "needs_input":
        record["status"] = "needs_input"
        record["needs_input"] = {
            "code": evaluation["reason_codes"][0] if evaluation["reason_codes"]
                    else "REQUIRED_OUTPUT_MISSING",
            "stage_id": None,
            "message": "必須出力の最終判定でneeds_inputが検出されました。",
            "details": {"missing_outputs": evaluation["missing_outputs"]},
        }
    else:
        record["status"] = "partial" if evaluation["achieved_outputs"] else "failed"
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return store.load_run(pipeline_path)


def finish_interrupted(pipeline_path: Path, record: dict, outcome: dict) -> dict:
    """needs_input/failed/取消で停止する。

    needs_inputはそのままneeds_input。failedは「有効な解析出力
    （`record['results']`が非空）」があればpartial、無ければfailedへ格下げする
    （brief「有効な解析出力があればpartial、なければfailed」）。

    ただし**取消が要求済みなら`cancelled`**にする。`cancel_requested`の判定は
    工程の境界でしか走らないので、長い工程（実測13分半のConsole）の最中に
    取り消すと、監視経路がConsoleを殺してhandlerがfailedを返し、ここへ来る。
    格下げ規則だけを当てるとpartialになり、`pipeline_cancel`が契約している
    「`cancelled`への確定を`pipeline_status`で確認する」が永久に成立しない。
    工程側の停止理由（`stages[...]["error"]`）は既に
    `commit_stage_outcome`が記録済みなので、ここで失われるものは無い。
    """
    record = copy.deepcopy(record)
    if outcome.get("status") == "failed" and cancel_requested(pipeline_path):
        record["status"] = "cancelled"
        store.save_run(pipeline_path, record, expected_revision=record["state_revision"])
        return store.load_run(pipeline_path)
    if outcome.get("status") == "needs_input":
        record["status"] = "needs_input"
        error = outcome.get("error") or {}
        record["needs_input"] = {
            "code": error.get("code"),
            "stage_id": outcome.get("stage_id"),
            "message": error.get("message"),
            "details": error.get("details") or {},
        }
    else:
        record["status"] = "partial" if record.get("results") else "failed"
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return store.load_run(pipeline_path)


def finish_cancelled(pipeline_path: Path, record: dict) -> dict:
    record = copy.deepcopy(record)
    record["status"] = "cancelled"
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])
    return store.load_run(pipeline_path)


# ---------- handler呼出し ----------

def _domain_error_payload(exc: DomainError) -> dict:
    return {"code": exc.code, "message": exc.message, "details": dict(exc.details)}


def _invoke_handler(handler: Callable[[dict], dict], context: dict) -> dict:
    """handlerを呼び、例外を成功結果へ絶対に変換しないStageResultへ落とす。"""
    try:
        outcome = handler(context)
    except DomainError as exc:
        if exc.code in _NEEDS_INPUT_CODES:
            return {"status": "needs_input", "result_refs": [], "warnings": [],
                    "error": _domain_error_payload(exc)}
        _logger.exception("pipeline stage %s: 入力解消不可のDomainError",
                          context.get("stage_id"))
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
    """`record["request"]["saved_path"]` の完全なpipeline-request.v1を読む。

    `record["request"]` 自体はrevision/request_id/content_hash/saved_path/
    effective_targetだけを持つ要約（`store.create_run` 参照）で、
    `target`/`comparisons` 等の本体は別ファイルにある。
    """
    saved_path = record["request"]["saved_path"]
    full_path = Path(pipeline_path) / saved_path
    return json.loads(full_path.read_text(encoding="utf-8"))


def _stage_order(record: dict, plan: dict) -> list[str]:
    """実行順を決める。**永続dictの挿入順ではなく`build_stages`の計画順**が正準。

    `record["stages"]`の挿入順は、`recovery.prepare_resume`がresume時に
    「recordに無いstage_id」を末尾へ**追記**するせいで計画順と一致しなくなる
    ——比較を後から足したrunでは`differential:*`/`export:*`が`report`の
    **後ろ**に並ぶ。その順で回すと、比較を1件も実行していない時点でレポートを
    書き、そのレポート自身が必須出力の欠落を告げているのにrunはcompletedに
    なる（`docs/workflow/pipeline.md`の手順21はreportを差次的解析・exportの
    後に置いている）。

    計画に無いstage_id（exploratory目標の`resolve_comparisons`など、record
    だけが持つもの）は計画順のあとへ回す——`_mark_out_of_scope`がhandlerを
    呼ばずにskipとして畳むだけなので、どこで処理しても順序上の意味は無い。
    計画にあってrecordに無いstage_idは`record["stages"][...]`を引けないため
    ここでは扱わない（`build_stages`と`store._stage_ids_for`は同じ規則で
    組み立てられており、通常は起こらない）。
    """
    stages = record["stages"]
    ordered = [sid for sid in plan if sid in stages]
    planned = set(ordered)
    ordered += [sid for sid in stages if sid not in planned]
    return ordered


def _request_identity(record: dict) -> tuple:
    """「今どの要求版を実行しているか」を表す最小の識別子。

    `state_revision`は使わない——無関係な同時更新（`pipeline_status`の
    probe・取消フラグの保存）でも進むため、これで再計画すると要求が
    変わっていないのにstage loopを何度も組み直すことになる。
    """
    request_meta = record.get("request") or {}
    return (request_meta.get("revision"), request_meta.get("content_hash"),
            request_meta.get("saved_path"))


def _resume_reset_seen(record: dict, visited: list[str]) -> bool:
    """このpassで既に通過したstageが、外から`pending`へ戻されたかを返す。

    `prepare_resume`は`rerun_upstream=True`や比較の訂正で該当stageを
    `pending`へ戻すが、要求内容そのものは変わらないことがある
    （`rerun_upstream`だけの再開）。要求版の変化だけを見ていると、
    その差し戻しを見ないまま先へ進み、古い前提のまま後続stageを
    `succeeded`にしてしまう。
    """
    stages = record.get("stages") or {}
    return any((stages.get(sid) or {}).get("status") == "pending" for sid in visited)


def _run_stage_loop(pipeline_path: Path, handlers: dict) -> dict:
    """要求版が変わっていないことを**stage境界ごとに**確かめながら1回分進める。

    要求（`pipeline-request.v1`）はloop開始時に1回読むだけでは足りない。
    workerが稼働している間に届いた`pipeline_resume`は——`prepare_resume`が
    「稼働中workerがいるならstatusを動かさない」規則を持つため——
    `pipeline_run.json`のrevisionと保存済み要求だけを差し替えて戻ってくる
    （`resume_pipeline`はstatusが`planned`でないので新しいworkerを起こさない）。
    このとき走り続けているworkerが**旧revisionの条件で計算した結果**を
    新revisionの結果として`succeeded`にすると、群の訂正が結果へ反映されないまま
    completedになる——planの中心要件（訂正・再開で結果を取り違えない）に
    正面から反する。

    そこで各stageの直前に最新recordを読み、要求版（`_request_identity`）か
    既に通過したstageの差し戻し（`_resume_reset_seen`）を観測したら、
    **その場でstage計画を組み直してpassをやり直す**。停止して次のresumeを
    待つのではなくこのworkerが続けるのは、稼働中workerがいる限り
    `resume_pipeline`が新しいworkerを起こさない（＝誰も再開しない）ため。

    `runtime`はpassをやり直すとき必ず捨てる。`DatasetState`も
    `runtime["differential"]`のキャッシュも旧要求の条件（旧シート・旧群割当）で
    組み立てられており、これを引き継ぐと「新しい要求で計算した」と称して
    古い数字を書き出すことになる。捨てても`_ALWAYS_RECONSTRUCT_STAGE_IDS`が
    load_dataset〜pcaを必ず呼び直すので、次のpassで作り直される。
    """
    record = store.load_run(pipeline_path)
    request = _load_request(pipeline_path, record)
    identity = _request_identity(record)
    runtime: dict = {}

    for _ in range(_REQUEST_RELOAD_MAX_PASSES):
        plan = {entry["stage_id"]: entry for entry in build_stages(request)}
        visited: list[str] = []
        restarted = False

        for stage_id in _stage_order(record, plan):
            latest = store.load_run(pipeline_path)
            if _request_identity(latest) != identity or _resume_reset_seen(latest, visited):
                record = latest
                request = _load_request(pipeline_path, latest)
                identity = _request_identity(latest)
                runtime.clear()
                restarted = True
                break
            record = latest
            visited.append(stage_id)

            entry = plan.get(stage_id)
            if entry is None:
                record = _mark_out_of_scope(pipeline_path, record, stage_id)
                continue

            stage = record["stages"][stage_id]
            if (stage["status"] in {"succeeded", "skipped"}
                    and stage_inputs_unchanged(stage, pipeline_root=pipeline_path)):
                continue

            if cancel_requested(pipeline_path):
                return finish_cancelled(pipeline_path, record)

            record = mark_stage_running(pipeline_path, record, stage)
            stage = record["stages"][stage_id]
            context = make_context(record, stage, runtime, request)
            outcome = _invoke_handler(handlers[entry["handler"]], context)
            record = commit_stage_outcome(pipeline_path, record, stage, outcome)

            if outcome["status"] in {"needs_input", "failed"}:
                return finish_interrupted(pipeline_path, record,
                                          {**outcome, "stage_id": stage_id})

        if not restarted:
            return finish_success(pipeline_path, record)

    raise DomainError(
        "REQUEST_REVISION_CHURN",
        f"要求の更新が続いて工程計画が{_REQUEST_RELOAD_MAX_PASSES}回組み直されました。"
        f"更新を止めてから再開してください: {pipeline_path}",
        {"pipeline_root": str(pipeline_path)})


def _record_worker_identity(pipeline_path: Path) -> None:
    """自分（現在プロセス）のidentityをrunへ刻む（spec §9.3「ロックにはowner
    identityを記録し」の実体）。owner lock取得直後、stage loop開始前に1回だけ
    呼ぶ。Task16 `metabolomix.pipeline.recovery.read_status`/`prepare_resume`が
    これを読み、`same_process`でこのworkerがまだ生きているかを判定する
    （pidだけでは再利用と見分けられないため、creation_timeと組みで記録する）。
    """
    record = store.load_run(pipeline_path)
    record = copy.deepcopy(record)
    record["worker"] = {"identity": process_identity(os.getpid()), "started_at": _now_iso()}
    store.save_run(pipeline_path, record, expected_revision=record["state_revision"])


def run_engine(pipeline_path: Path, handlers: dict[str, Callable[[dict], dict]]) -> dict:
    """永続化済みrunを、必要な工程がなくなるかstageが止まるまで1回分進める。

    `pipeline_path` は `store.create_run`/`find_or_create_run` が返す
    pipeline_root。戻り値は保存後に読み直した最新の `pipeline-run.json`。

    owner lock（`control/worker.lock`）をstage loop全体の間ずっと保持する。
    取得できなければ即 `PIPELINE_ALREADY_RUNNING`（別workerが実行中、または
    このプロセス自身の入れ子呼出し——`file_lock` は同一プロセスからの
    入れ子取得も別ロック要求として扱うため区別できない。区別する必要も無い:
    どちらの場合も「今は進めてはいけない」という結論は同じ）。

    ロック取得の成否だけを判定するために `__enter__`/`__exit__` を手で呼ぶ
    （`with`一つで両方を包むと、stage loop本体の中で別の場所——例えば
    `store.save_run` 自身の `control/pipeline.lock` 待ち——が`LOCK_TIMEOUT`に
    なった場合まで`PIPELINE_ALREADY_RUNNING`に化けてしまう。それは誤り:
    state更新lockの競合は「別workerが実行中」ではなく単なる一時的な混雑）。
    """
    pipeline_path = Path(pipeline_path)
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
