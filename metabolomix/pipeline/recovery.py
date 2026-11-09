"""再開・取消・状態の読取表示（spec §9.2/9.3・plan 2026-11-09）。

このモジュールが解く問題は3つ、公開関数もその3つに対応する。

1. **`read_status`**: pipeline-run.jsonを読むだけの、絶対に書き換えないstatus取得。
   `status="running"`のときだけ、記録済みの`worker`（`engine.run_engine`が刻む
   自分自身のidentity）が今も生きているかを`same_process`で見る。このprobeが
   どう転んでも（生存・消失・判定不能）ファイルには一切触れない。極性ごとの単位の
   状態（`units`）と、受領があればその要約（`receipt`: 解析フォルダの場所と `next`）を
   添える。旧 v1 の run も読め、`legacy: true` と注意 1 行を添える。
2. **`request_cancel`**: 協調的な取消フラグ（`engine.cancel_request_path`が指す小さな
   JSON）を書いて受理を返すだけ。プロセスを直接殺さない——Console の停止は
   `supervise` 側の `cancel_path` 監視に、工程の停止は engine の stage 境界チェックに
   委ねる。旧 v1 の run は `PIPELINE_RUN_LEGACY`。
3. **`prepare_resume`**: どの工程を`pending`へ戻すかを決めて保存するだけ——
   **launchしない**。旧 v1 の run は `PIPELINE_RUN_LEGACY`（記録を書き換えない）。
   以下の順で安全装置を通す:

   a. `updates` を `request.merge_updates` で検証する（下流のキーは
      `PIPELINE_DOWNSTREAM_REMOVED`、上流の条件は `NEW_PIPELINE_REQUIRED`）。
   b. 監視を失ったまま稼働中のConsoleがどれかの単位にあれば`EXECUTION_UNRESOLVED`。
   c. `rerun_upstream`が無いのに、Console の状態が分からない単位（upstream が
      `running` のまま）がある、または進められる単位が 1 つも無い（全単位が失敗）なら
      `UPSTREAM_RERUN_REQUIRED`（自動再試行はしない）。失敗した単位があっても他の単位を
      進められるなら止めず、失敗した単位はそのまま（`not_retried_units` で知らせる）。
   d. `rerun_upstream=True` は upstream が `succeeded` でない単位の 3 工程を `pending` へ
      戻す。成功した単位は戻さない（Console を再起動しない）。
   e. 戻さない単位の固定済み入力（raw stat・実効メソッド・注釈器の表）を再検証する。
   f. 試料の対応表の更新・シートの中身の変化は `receipt` だけを戻す（アライメントの
      参照ファイルは計画の時点で固定済み）。どれかを戻したら `receipt` も戻す。
   g. `planned`へ戻したなら、協調取消フラグを取り下げる（`_clear_cancel_request`）。

   永続statusが`running`のまま固まったrun（worker喪失）も、そのworkerがもう存在しないと
   観測できたときだけ`planned`へ戻す（`_RESUMABLE_RUNNING_HEALTH`）。

   `request_id`を指定した再送は、直前と同じ`updates`/`rerun_upstream`なら新しい
   revisionを作らず直前の結果を返す（`resume_log`）。

このモジュールはsessionを一切importしない（`metabolomix.core.session_state` /
`metabolomix.core.mcp_core` / `metabolomix.tools.*`）。
"""
from __future__ import annotations

import copy
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.console import execution as console_execution
from metabolomix.console import method_file as method_file_mod
from metabolomix.core.atomic_io import DomainError, atomic_write_json, canonical_hash
from metabolomix.core.process_control import same_process
from metabolomix.pipeline import engine
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import store

__all__ = ["normalize_pipeline_root", "prepare_resume", "read_status", "request_cancel"]

#: prepare_resumeの状態競合吸収リトライ上限（有限回の読み直しで十分安全に吸収できる）。
_RESUME_RETRY_MAX_ATTEMPTS = 5

#: 中断由来の終端状態。`completed`と違い、何も変えなくてもresumeが呼ばれた
#: 時点で常に`planned`へ戻してよい（D09「取消後・timeout後のresume」）。
_INTERRUPTED_STATUSES = frozenset({"partial", "failed", "cancelled", "needs_input"})

#: 永続statusが`running`のまま残ったrunを「中断」と見なしてよい観測結果。
#: **`"unknown"`は入れない。** identityが記録されていない・生存確認そのものが
#: 失敗した、のいずれも「死んでいる証拠」ではない。
_RESUMABLE_RUNNING_HEALTH = frozenset({"worker_missing"})

_LEGACY_NOTE = ("旧形式（pipeline-run.v1）の run です。読むだけで、再開・取消はできません"
                "（読むなら load_dataset(pipeline_path=…)、Console をやり直すなら新しい pipeline_run）。")


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_pipeline_root(path: Path) -> Path:
    """`path`がpipeline-run.jsonそのものでも、その親（pipeline_root）でも受け付ける。

    `pipeline_status`/`pipeline_cancel`/`pipeline_resume` と `load_dataset(pipeline_path=…)`
    が同じ規則で揃える。
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

    `"ok"` / `"worker_missing"` / `"unknown"`（identityが無い、または生存確認自体が失敗した）。
    """
    identity = (record.get("worker") or {}).get("identity")
    if not identity:
        return "unknown"
    try:
        alive = same_process(identity)
    except Exception:  # noqa: BLE001 - probe自体の失敗をdeadと解釈しない
        return "unknown"
    return "ok" if alive else "worker_missing"


def _receipt_summary(pipeline_root: Path, record: dict) -> dict | None:
    """受領の要約（解析フォルダの場所・照合した `.arf`・`next`）。受領が無ければ None。"""
    data = receipt_mod.read_receipt(pipeline_root, record)
    if data is None:
        return None
    units = data.get("units") or {}
    validated = [p for p, u in units.items() if u.get("status") == "validated"]
    next_call = data.get("next")
    if next_call:
        next_call = next_call.replace(receipt_mod.PIPELINE_PATH_PLACEHOLDER,
                                      str(pipeline_root).replace("\\", "/"))
    return {
        "relative_path": receipt_mod.RECEIPT_RELATIVE_PATH,
        "analysis_dirs": {p: str(pipeline_root / units[p]["analysis_dir"]) for p in validated},
        "alignment": {p: units[p].get("alignment") for p in validated},
        "warnings": [w.get("code") for w in data.get("warnings") or []],
        "next": next_call,
    }


def read_status(path: Path, *, include_details: bool = False) -> dict:
    """pipeline状態を読むだけで、一切書き換えない。

    `status`が"running"のときだけ、記録済みworker identityの生死を`same_process`
    で見る。この生存確認自体が失敗しても"dead"とは読まず、`observed_health="unknown"`
    として呼び出し側に復旧判断を委ねる。
    """
    pipeline_root = normalize_pipeline_root(path)
    record = store.load_run(pipeline_root)  # 読取専用。ここでは一切保存しない。

    status = record.get("status")
    observed_health = _worker_health(record) if status == "running" else "ok"

    result = {
        "pipeline_id": record.get("identity", {}).get("pipeline_id"),
        "status": status,
        "observed_health": observed_health,
        "needs_input": record.get("needs_input"),
        "warnings": list(record.get("warnings") or []),
        "stage_statuses": _stage_statuses(record),
        "units": {p: store.unit_status(record, p) for p in store.run_units(record)},
    }
    if store.is_legacy(record):
        result["legacy"] = True
        result["note"] = _LEGACY_NOTE
    else:
        summary = _receipt_summary(pipeline_root, record)
        if summary is not None:
            result["receipt"] = summary
    if status == "running" and observed_health in {"worker_missing", "unknown"}:
        result["recovery_hint"] = {
            "code": "PIPELINE_INTERRUPTED" if observed_health == "worker_missing"
                    else "SUPERVISION_UNKNOWN",
            "message": "workerの生存を確認できません。prepare_resumeで再開してください。",
        }
    if include_details:
        result["record"] = _without_library_location(record)
    return result


#: 研究室 MSP の注釈器設定表を指すキー（`method_file.MSP_ANNOTATOR_SETTINGS_KEYS` の小文字）。
_MSP_TABLE_KEYS = frozenset(key.lower() for key in method_file_mod.MSP_ANNOTATOR_SETTINGS_KEYS)


def _redact_inputs(inputs: dict) -> None:
    """入力計画（の写し）から研究室ライブラリの置き場所を除く（ファイル名だけ残す）。

    研究室 MSP（`msp.path` と `Msp file path` の上書き）、`[library] lbm` 由来の LBM
    （`lbm.path` と `Lbm file path` の上書き。`inputs.is_lab_lbm`）、絶対パスで書かれた
    MSP 注釈器設定表の上書きが対象。
    """
    lab_lbm = inputs_mod.is_lab_lbm(inputs)
    sections = [inputs["msp"]] if isinstance(inputs.get("msp"), dict) else []
    if lab_lbm and isinstance(inputs.get("lbm"), dict):
        sections.append(inputs["lbm"])
    for section in sections:
        if section.get("path"):
            section["path"] = Path(section["path"]).name
    overrides = (inputs.get("method") or {}).get("overrides")
    if isinstance(overrides, dict):
        for key, value in overrides.items():
            if not value or not isinstance(value, str):
                continue
            lowered = key.lower()
            if (lowered == "msp file path"
                    or (lab_lbm and lowered == "lbm file path")
                    or (lowered in _MSP_TABLE_KEYS and Path(value).is_absolute())):
                overrides[key] = Path(value).name


def _without_library_location(record: dict) -> dict:
    """戻り値用の record の写しから、研究室ライブラリの置き場所を除く。

    戻り値は LLM の文脈に入る（spec 2026-10-07 §3.3）。保存された record は変えない
    ——resume の再検査（`inputs.verify_inputs`）が絶対パスを要る。
    """
    redacted = copy.deepcopy(record)
    if store.is_legacy(redacted):
        if isinstance(redacted.get("inputs"), dict):
            _redact_inputs(redacted["inputs"])
        return redacted
    for unit in (redacted.get("units") or {}).values():
        if isinstance(unit, dict) and isinstance(unit.get("inputs"), dict):
            _redact_inputs(unit["inputs"])
    return redacted


# ---------- request_cancel ----------

def request_cancel(path: Path) -> dict:
    """協調的な取消フラグを保存し、受理だけを返す（`cancelled`は実際にworkerが停止を
    確認した後に確定する）。プロセスを直接殺さない。
    """
    pipeline_root = normalize_pipeline_root(path)
    record = store.load_run(pipeline_root)  # 存在・schema確認（読取専用）。
    if store.is_legacy(record):
        raise store.legacy_refusal(pipeline_root, "取消")

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


def _console_supervision_states(record: dict, pipeline_root: Path) -> list[tuple[dict, bool]]:
    """単位ごとの (Console process identity, 検証済みか)。

    `verified` はその単位の `upstream` 工程が succeeded かどうか。Console job をまだ
    起動していない単位は数えない。`console_job_path`（pipeline_root 相対の
    `analysis-job.json`）は記録されているのに、未検証のまま終了証跡が読めない・
    process_identity を持たない場合は `EXECUTION_UNRESOLVED`（監視を失って稼働中かも
    しれない Console を黙って verified 扱いしない）。
    """
    states: list[tuple[dict, bool]] = []
    stages = record.get("stages") or {}
    for polarity, unit in store.run_units(record).items():
        verified = (stages.get(store.unit_stage_id("upstream", polarity)) or {}).get(
            "status") == "succeeded"
        job_rel = (unit.get("upstream") or {}).get("console_job_path")
        if not job_rel:
            continue
        run_dir = (Path(pipeline_root) / job_rel).parent
        try:
            data = json.loads(console_execution.receipt_path(run_dir).read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            if verified:
                continue
            raise DomainError(
                "EXECUTION_UNRESOLVED",
                f"console_job_pathは記録されていますが終了証跡を読めません: {run_dir}",
                {"console_job_path": str(job_rel), "run_dir": str(run_dir), "unit": polarity},
            ) from exc
        identity = data.get("process_identity")
        if not isinstance(identity, dict):
            if verified:
                continue
            raise DomainError(
                "EXECUTION_UNRESOLVED",
                f"終了証跡がprocess_identityを持たない不正な形式です: {run_dir}",
                {"console_job_path": str(job_rel), "run_dir": str(run_dir), "unit": polarity})
        states.append((identity, verified))
    return states


def _verify_effective_method(snapshot: dict, unit_root: Path) -> None:
    """Console が実際に読む `inputs/effective-method.txt` の hash を照合する
    （`inputs.verify_inputs` は原本メソッドの hash しか見ない）。"""
    method = snapshot.get("method") or {}
    rel = method.get("effective_relative_path")
    expected = method.get("effective_sha256")
    if not rel or not expected:
        return  # 実効メソッドを固定していない run（配置前）は対象外
    target = Path(unit_root) / rel
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


def _verify_msp_annotator_table(snapshot: dict, unit_root: Path) -> None:
    """library_mode="msp_then_lbm" の MSP 注釈器設定表（`stage_inputs` が書いて sha256 を
    固定した `inputs/msp-annotators.tsv`）を実ファイルと照合する。

    文面と details は単位の根の相対のパスだけ（表の中身は研究室 MSP の置き場所を持つ）。
    """
    msp = snapshot.get("msp") or {}
    rel = msp.get("annotator_table_relative_path")
    expected = msp.get("annotator_table_sha256")
    if not rel or not expected:
        return
    try:
        actual = _sha256_file(Path(unit_root) / rel)
    except OSError as exc:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"MSP の注釈器設定表が読めません（再開時の照合）: {rel}",
            {"path": rel},
        ) from exc
    if actual != expected:
        raise DomainError(
            "STAGED_INPUT_MISMATCH",
            f"MSP の注釈器設定表が元の固定内容と一致しません（再開時の照合）: {rel}",
            {"path": rel},
        )


def _verify_staged_inputs(record: dict, polarity: str) -> None:
    """単位の固定済み入力を再検証する（D05）。やり直す単位には呼ばない
    ——上流をやり直す以上、`prepare_input` が自分で再検証・再配置する。
    """
    pipeline_root = Path(record["identity"]["pipeline_root"])
    unit = store.run_units(record)[polarity]
    unit_root = pipeline_root / (unit.get("unit_root") or polarity)
    snapshot = store.absolutize_unit_inputs(unit.get("inputs") or {}, unit_root)
    if not snapshot.get("source_root") or not snapshot.get("raw_stat"):
        return  # 検証対象を持たない最小fixture（単体test）向けの安全側スキップ
    inputs_mod.verify_inputs(snapshot)
    _verify_effective_method(snapshot, unit_root)
    _verify_msp_annotator_table(snapshot, unit_root)


def _manifest_content_changed(pipeline_root: Path, record: dict, request: dict) -> bool:
    """実験情報シートの**中身**が、受領を書いたときから変わったかを返す。

    要求の`sample_manifest`はパス文字列でしかないため、同じシートを書き直して
    `pipeline_resume`しても、要求の比較では「変更なし」に見える。比較の基準は受領が
    記録した `samples.manifest.sha256`（そのとき読んだシートの内容hash）。受領が無い・
    シートが無かった単位は比較材料が無いので数えない（「判定できない」を「変わった」と
    読み替えない）。
    """
    data = receipt_mod.read_receipt(pipeline_root, record)
    if data is None:
        return False
    dataset_root = Path(record["identity"]["source_root"])
    for index, unit in enumerate(request.get("units") or []):
        recorded = (((data.get("units") or {}).get(unit.get("polarity")) or {})
                    .get("samples") or {}).get("manifest") or {}
        previous = recorded.get("sha256")
        if not previous:
            continue
        current = inputs_mod.manifest_source_record(
            request_mod.unit_source_root(dataset_root, unit),
            request_mod.unit_request(request, index, dataset_root))
        if current.get("sha256") != previous:
            return True
    return False


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

    消さないと、`planned`へ戻して新しいworkerを起こしても最初のstage境界で
    `cancel_requested()`を見て即座に`finish_cancelled`するため、取消したrunは
    二度と再開できない。

    ただし停止処理がまだ進行中かもしれない取消は握り潰さない。`prior_status`が
    `running`/`planned`なら、そのworkerがもう居ないと観測できたとき
    （`worker_health == "worker_missing"`）以外は消さずに残す。`"unknown"`を
    「居ない」と読まない（起動直後の worker はまだ identity を刻んでいない）。
    """
    if prior_status in {"running", "planned"} \
            and worker_health not in _RESUMABLE_RUNNING_HEALTH:
        return
    engine.cancel_request_path(pipeline_root).unlink(missing_ok=True)


def _resume_receipt(record: dict, *, reused: bool, reset_stage_ids: list | None = None,
                    not_retried_units: list | None = None,
                    held_units: list | None = None) -> dict:
    out = {
        "pipeline_id": record["identity"]["pipeline_id"],
        "status": record["status"],
        "request_revision": record["request"]["revision"],
        "reused": reused,
        "reset_stage_ids": sorted(reset_stage_ids or []),
    }
    if not_retried_units:
        out["not_retried_units"] = sorted(not_retried_units)
        out["hint"] = ("失敗した極性の Console は自動ではやり直しません。"
                       "やり直すなら rerun_upstream=true で再開してください。")
    if held_units:
        # rerun_upstream=true でも、生きている worker が回している単位は戻さない。黙って
        # reset_stage_ids から消えるだけだと「やり直しが効かなかった」と読まれる。
        out["held_units"] = sorted(held_units)
        out["hint"] = (f"極性 {', '.join(sorted(held_units))} は Console が稼働中なので"
                       "やり直していません（worker がそのまま最後まで進めます）。"
                       "pipeline_status で結果を見てから、なお必要なら rerun_upstream=true で"
                       "再開してください。")
    return out


# ---------- prepare_resume ----------

def prepare_resume(path: Path, *, updates: dict | None = None,
                   request_id: str | None = None, rerun_upstream: bool = False) -> dict:
    """再開の判断（どの工程を戻すか・新しい要求の版）を保存する。**launchしない**。

    手順は本モジュールdocstringのa〜gのとおり。`STATE_REVISION_CONFLICT`
    （他アクターがこのrunを同時に更新した）は有限回まで読み直して吸収する。
    """
    pipeline_root = normalize_pipeline_root(path)
    updates = dict(updates) if updates else {}
    updates_hash = canonical_hash({"updates": updates, "rerun_upstream": rerun_upstream})

    for _ in range(_RESUME_RETRY_MAX_ATTEMPTS):
        record = store.load_run(pipeline_root)
        if store.is_legacy(record):
            raise store.legacy_refusal(pipeline_root, "再開")

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

        # --- a. updatesの検証・統合（下流のキー・上流の条件の変更はここで止める） ---
        current_request = _load_full_request(pipeline_root, record)
        merged_request = (request_mod.merge_updates(current_request, updates)
                          if updates else current_request)

        # --- b. 監視を失った稼働中Console ---
        for identity, verified in _console_supervision_states(record, pipeline_root):
            if not verified and same_process(identity):
                raise DomainError("EXECUTION_UNRESOLVED", "監視を失ったConsoleが稼働しています")

        # --- c. 上流再実行の明示要求 ---
        units = store.run_units(record)
        upstream_status = {
            p: (record["stages"].get(store.unit_stage_id("upstream", p)) or {}).get("status")
            for p in units}
        unit_states = {p: store.unit_status(record, p) for p in units}
        # 生きている worker が今まさに回している単位（工程が running）は、監視が生きて
        # いるので「状態が分からない」ではない。やり直しの対象にも、未解決にも数えない
        # ——戻すと worker はその Console の結果を捨て、Console をもう 1 回起動する。
        worker_health = _worker_health(record)
        live_worker = record["status"] in ("running", "planned") and worker_health == "ok"
        held = ({p for p in units
                 if any((record["stages"].get(store.unit_stage_id(h, p)) or {}).get("status")
                        == "running" for h in store.UNIT_HANDLERS)}
                if live_worker else set())
        if not rerun_upstream:
            unresolved = sorted(p for p, s in upstream_status.items()
                                if s == "running" and p not in held)
            if unresolved or all(s == "failed" for s in unit_states.values()):
                raise DomainError(
                    "UPSTREAM_RERUN_REQUIRED", "上流の再実行を明示してください",
                    {"units": unresolved or sorted(unit_states),
                     "unit_statuses": unit_states})

        # --- d. やり直す単位（upstream が succeeded でない・失敗した単位。稼働中は除く） ---
        rerun_units = (sorted(p for p, s in upstream_status.items()
                              if p not in held
                              and (s != "succeeded" or unit_states[p] == "failed"))
                       if rerun_upstream else [])
        not_retried = ([] if rerun_upstream
                       else sorted(p for p, s in unit_states.items() if s == "failed"))

        # --- e. 戻さない単位の固定済み入力の再検証 ---
        for polarity in units:
            if polarity in rerun_units or unit_states[polarity] == "failed":
                continue
            prepare = record["stages"].get(store.unit_stage_id("prepare_input", polarity)) or {}
            if prepare.get("status") == "succeeded":
                _verify_staged_inputs(record, polarity)

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
            }

        # --- f. 工程のpending化 ---
        for sid in store.stage_ids_for(merged_request):
            if sid not in record["stages"]:
                record["stages"][sid] = store.initial_stage(sid)
        manifest_changed = _manifest_content_changed(pipeline_root, record, merged_request)
        reset_ids: set = set()
        for polarity in rerun_units:
            for handler in store.UNIT_HANDLERS:
                reset_ids.add(store.unit_stage_id(handler, polarity))
        if reset_ids or request_changed or manifest_changed:
            reset_ids.add(engine.RECEIPT_STAGE_ID)
        for sid in reset_ids:
            record["stages"][sid] = _reset_stage_for_resume(record["stages"][sid])

        # 中断由来の終端状態（partial/failed/cancelled）とneeds_inputは、resumeが
        # 呼ばれた時点で常に`planned`へ戻す。`completed`は、実際に何か変わった場合
        # だけ`planned`へ戻す（no-op resumeでcompletedを崩さない）。`running`のまま
        # 固まったrunは、そのworkerがもう居ないと観測できたときだけ同じ扱いにする。
        made_changes = request_changed or bool(reset_ids) or manifest_changed
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
        return _resume_receipt(saved, reused=False, reset_stage_ids=sorted(reset_ids),
                               not_retried_units=not_retried,
                               held_units=sorted(held) if rerun_upstream else None)

    raise DomainError(
        "STATE_REVISION_CONFLICT",
        f"resumeの状態競合の再試行上限（{_RESUME_RETRY_MAX_ATTEMPTS}回）に達しました: "
        f"{pipeline_root}",
        {"pipeline_root": str(pipeline_root)})
