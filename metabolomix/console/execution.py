"""Console終了証跡（console-execution.v1）の検証と、同期・非同期共通の監視。

前半は `execution-result.json` の契約。`validate_execution_record`は検証済みdictだけを
消費者へ渡すための唯一の入口で、妥当性を確認しないまま `pipeline_status` 等へ
流用してはいけない。

後半の `supervise()` が **Console実行の唯一の監視経路**。同期呼出し（console_run）も
切り離しワーカー（Task 5 の `metabolomix.console.worker`）も、pipelineの上流工程も
ここを通る。旧経路は非ゼロ終了で即returnして収集を飛ばしていたため、失敗した実行が
何の証拠も残さなかった。`supervise`は成功・非ゼロ・timeout・取消・起動不能の
**どの終了経路でも証跡を保存する**。

段階の順序には理由がある。証跡（終了の事実）を先に確定し、そのあとで収集 →
出力検証 → job保存を順に記録する。逆順にすると、収集や保存が落ちたときに
「Consoleがどう終わったか」まで失われる。例外はどの段階でも成功へ変換しない。

このモジュールが run_dir へ書く運用ファイルは 2 つだけで、どちらも
`metabolomix.console.output_collector._OPERATIONAL_FILES` に載っている
（載っていないと解析成果物として誤収集され、完了判定が静かに劣化する）:

- `execution-result.json` … 終了証跡（console-execution.v1）
- `worker.json` … 監視入力・起動情報・復旧用スナップショット（console-supervision.v1）
"""
from __future__ import annotations

import json
import os
import re
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

from metabolomix.console.job_manager import load_job, save_job
from metabolomix.console.output_collector import RUNS_SUBDIR, collect_artifacts, snapshot
from metabolomix.console.validation import completion_status, validate_outputs
from metabolomix.core.atomic_io import DomainError, atomic_write_json, canonical_hash
from metabolomix.core.process_control import start_owned_process
from metabolomix.handoff.schema import AnalysisJob, sha256_file

SCHEMA = "console-execution.v1"

# spec 5.1: 終了証跡の termination 種別。
TERMINATIONS = frozenset({"exited", "timeout", "cancelled", "worker_lost", "launch_failed"})

# 起動自体に失敗した／監視を見失ったケースは、pid・process_identityの不明を
# 起動情報（command/method/exe の各hash）とは別に許可する（spec 5.1）。
_UNKNOWN_PID_TERMINATIONS = frozenset({"worker_lost", "launch_failed"})

_HASH_FIELDS = ("command_sha256", "method_sha256", "exe_sha256")
_ID_FIELDS = ("execution_id", "job_id")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _fail(message: str, details: dict | None = None) -> None:
    raise DomainError("EXECUTION_RECORD_INVALID", message, details)


def validate_exit_fields(data: dict) -> None:
    """exit_codeの型とterminationとの整合を検査する。

    exit_codeキー自体が無い入力も、値がnullの入力と同様に「未回収」として
    扱い、辞書直接indexingでKeyErrorを送出しない（このモジュールの他の
    フィールド検査はすべて`.get()`経由であり、契約はDomainErrorのみ）。
    """
    rc = data.get("exit_code")
    if rc is not None and type(rc) is not int:
        raise DomainError("EXECUTION_RECORD_INVALID", "終了コードが不正です")
    if data.get("termination") == "exited" and rc is None:
        raise DomainError("EXECUTION_RECORD_INVALID", "終了コードがありません")


def _parse_utc(value: object, field_name: str) -> datetime:
    if not isinstance(value, str) or not value:
        _fail(f"{field_name} が空です", {field_name: value})
    text = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        dt = datetime.fromisoformat(text)
    except ValueError:
        _fail(f"{field_name} の日時形式が不正です", {field_name: value})
    if dt.tzinfo is None or dt.utcoffset() != timedelta(0):
        _fail(f"{field_name} はUTC（オフセット+00:00）である必要があります", {field_name: value})
    return dt.astimezone(timezone.utc)


def _validate_pid_and_identity(data: dict) -> None:
    pid = data.get("pid")
    identity = data.get("process_identity")
    termination = data.get("termination")

    if pid is None:
        if termination not in _UNKNOWN_PID_TERMINATIONS:
            _fail("pid が不正です（起動失敗・監視喪失以外でnullは許可されません）",
                 {"pid": pid, "termination": termination})
        if identity is not None:
            _fail("pidがnullのprocess_identityはnullでなければなりません", {"process_identity": identity})
        return

    if type(pid) is not int or pid <= 0:
        _fail("pid が不正です", {"pid": pid})
    if not isinstance(identity, dict) or identity.get("pid") != pid:
        _fail("process_identity がpidと整合しません", {"pid": pid, "process_identity": identity})
    if "creation_time" not in identity:
        _fail("process_identity に creation_time がありません", {"process_identity": identity})


def validate_execution_record(data: dict) -> dict:
    """console-execution.v1 の必須フィールドと整合性を検証する。

    妥当なら検証済みdict（浅いコピー）を返す。不正ならDomainError
    （code="EXECUTION_RECORD_INVALID"）を送出する。呼び出し側は検証済み
    dictだけを保存・伝播させ、生の入力をそのまま流用しない。
    """
    if data.get("schema") != SCHEMA:
        _fail(f"schema が {SCHEMA} ではありません", {"schema": data.get("schema")})

    for field_name in _ID_FIELDS:
        value = data.get(field_name)
        if not isinstance(value, str) or not value:
            _fail(f"{field_name} が空です", {field_name: value})

    termination = data.get("termination")
    if termination not in TERMINATIONS:
        _fail("termination が不正です", {"termination": termination})

    timeout_s = data.get("timeout_s")
    if type(timeout_s) is not int or timeout_s <= 0:
        _fail("timeout_s は正の整数である必要があります", {"timeout_s": timeout_s})

    for field_name in _HASH_FIELDS:
        value = data.get(field_name)
        if not isinstance(value, str) or not _HEX64.match(value):
            _fail(f"{field_name} が64桁の16進数ではありません", {field_name: value})

    started_at = _parse_utc(data.get("started_at"), "started_at")
    ended_at = _parse_utc(data.get("ended_at"), "ended_at")
    if ended_at < started_at:
        _fail("ended_at が started_at より前です",
             {"started_at": data.get("started_at"), "ended_at": data.get("ended_at")})

    validate_exit_fields(data)
    _validate_pid_and_identity(data)

    return dict(data)


# =====================================================================
# 監視（supervise）: 同期・非同期に共通する唯一のConsole実行経路
# =====================================================================

#: run_dirのサイドカー（監視入力・起動情報・復旧用スナップショット）のスキーマ。
SUPERVISION_SCHEMA = "console-supervision.v1"

#: run_dirへ書く運用ファイル。どちらも output_collector._OPERATIONAL_FILES に
#: 載っている名前を使う（新しい名前を増やすとそちらの除外集合の更新が要る）。
RECEIPT_FILENAME = "execution-result.json"
SUPERVISION_FILENAME = "worker.json"

#: Consoleのstdout/stderrを流し込むログ。runner._open_log と同じ名前。
CONSOLE_LOG_FILENAME = "msdial.log"

#: Consoleの -o に渡すサブディレクトリ。runner.build_msdial_cmd の呼び出し規約。
MSDIAL_OUT_SUBDIR = "msdial"

#: 監視ループの見張り間隔（秒）。取消要求とtimeoutの検出粒度。
POLL_INTERVAL_S = 0.1

#: terminate_tree後に終了コードを回収するまでの猶予（秒）。
REAP_TIMEOUT_S = 5

#: 証跡へ載せる可変長リストの上限。戻り値がそのままLLMの文脈を占めるため、
#: 検証の詳細は件数を絞る（全文はvalidate_outputsの呼び出し側で読める）。
_MAX_RECORDED = 5


def receipt_path(run_dir: Path) -> Path:
    """終了証跡（console-execution.v1）の保存先。"""
    return Path(run_dir) / RECEIPT_FILENAME


def supervision_state_path(run_dir: Path) -> Path:
    """監視サイドカー（console-supervision.v1）の保存先。"""
    return Path(run_dir) / SUPERVISION_FILENAME


def read_supervision_state(run_dir: Path) -> dict:
    """監視サイドカーを読む。無い・壊れている場合は空dict。

    読めないことを例外にしない。サイドカーは補助情報であり、これが無いせいで
    「終了証跡すら残さない」状態にしてはいけない。
    """
    try:
        data = json.loads(supervision_state_path(run_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def write_supervision_inputs(run_dir: Path, inputs: dict) -> Path:
    """監視入力（input inventory・method / exe のhash）をrun_dirへ固定する。

    `supervise`を呼ぶ側が「何を入力として実行するか」を先に宣言する場所。

    `inputs`のキー:

    ``raw_inventory``
        Consoleが実際に読む raw のパス一覧（pipelineでは`input/`配下の準備済み
        ファイル。元の`source_root`ではない）。`validate_outputs`の
        `expected_sources`にそのまま渡り、mzTabの`ms_run[N]-location`と1対1で
        突き合わされる。**無い場合は`INPUT_INVENTORY_MISSING`を立てて
        completedへ進ませない**（分からないものを「問題なし」と読み替えない）。
    ``method_sha256`` / ``exe_path`` / ``exe_sha256``
        起動時に固定する情報。省略時は`job.method_file`と`exe_path`から算出し、
        読めなければ`_unavailable_hash`の標識を置いて`inputs.unavailable`に残す。

    既存のサイドカーへ併合して書く。起動情報や復旧用スナップショットを
    後から消さないため、丸ごと置換はしない。
    """
    state = read_supervision_state(run_dir)
    state["schema"] = SUPERVISION_SCHEMA
    state["inputs"] = {**state.get("inputs", {}), **inputs}
    path = supervision_state_path(run_dir)
    atomic_write_json(path, state)
    return path


def _utc_now() -> str:
    """UTCの現在時刻（オフセット+00:00付き）。証跡の時刻はすべてこれで作る。"""
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def _unavailable_hash(field: str) -> str:
    """内容を読めなかったhash欄に入れる決定的な標識。

    証跡のhash欄は64桁hexという契約なので空文字を入れられない。かといって
    0埋めは「実在するhash」と見分けがつかない。読めなかったという事実そのものを
    hashして、`inputs.unavailable`にどの欄が標識かを併記する。
    """
    return canonical_hash({"unavailable": field})


def _file_sha256(path: object) -> str | None:
    """ファイルのSHA-256。読めなければNone（例外にしない）。"""
    if not path:
        return None
    try:
        return sha256_file(Path(str(path)))
    except OSError:
        return None


def raw_fingerprint(paths: list[str]) -> str:
    """生入力のstat fingerprint（パス・サイズ・mtime）を返す。

    spec 4.2の区別に従い、これは**内容検証ではない**。rawの全量SHA-256は初期版の
    必須条件にしないと決めてあるので、実行中に入力が差し替わったことを安く
    検出するためのstat比較に留める。
    """
    entries = []
    for raw in sorted((str(p) for p in paths), key=str.lower):
        try:
            st = os.stat(raw)
            entries.append({"path": raw, "size": st.st_size, "mtime_ns": st.st_mtime_ns})
        except OSError:
            entries.append({"path": raw, "size": None, "mtime_ns": None})
    return canonical_hash(entries)


def _monitor(proc, *, timeout_s: int, cancel_path: Path | None) -> tuple[str, int | None]:
    """終了・timeout・取消のいずれかまで見張り、(termination, exit_code)を返す。

    期限は`time.monotonic()`で測る。壁時計だとNTP補正や夏時間で期限が前後し、
    「44分かかる実行が10分で打ち切られる」「timeoutが永久に来ない」が起きる。

    停止させたのに終了コードを回収できなかった場合は0で補わず`worker_lost`へ
    落とす（回収不能を成功と読み替えない。spec 5.1）。
    """
    deadline = time.monotonic() + timeout_s
    termination = "exited"
    while proc.poll() is None:
        if cancel_path is not None and cancel_path.exists():
            termination = "cancelled"
            proc.terminate_tree()
            break
        if time.monotonic() >= deadline:
            termination = "timeout"
            proc.terminate_tree()
            break
        time.sleep(POLL_INTERVAL_S)
    try:
        exit_code = proc.wait(timeout=REAP_TIMEOUT_S)
    except DomainError:
        return ("worker_lost" if termination == "exited" else termination), None
    return termination, exit_code


def _default_command(job: AnalysisJob, dataset_root: Path, msdial_out_dir: Path,
                     inputs: dict) -> list[str]:
    """jobから実Consoleのコマンドラインを組む。

    `metabolomix.console.runner`をモジュール先頭でimportしない。runner側は共通経路
    （このモジュール）を指す立場なので、双方をmodule levelでimportすると循環になる。
    """
    from metabolomix.console.runner import build_msdial_cmd, get_exe_path

    exe = inputs.get("exe_path") or get_exe_path()
    return build_msdial_cmd(exe, dataset_root, msdial_out_dir,
                            Path(job.method_file), job.save_project)


def _write_receipt(run_dir: Path, receipt: dict) -> Path:
    """証跡を検証してから原子的に保存する。

    自分が組んだ記録が契約に反していても保存自体は行い、理由を`record_error`へ
    残す。保存を諦めると「Consoleがどう終わったか」ごと失われるため。読み手は
    `validate_execution_record`を通してから使う契約なので、不正な記録が黙って
    妥当なものとして流通することはない。
    """
    path = receipt_path(run_dir)
    try:
        validate_execution_record(receipt)
        receipt.pop("record_error", None)
    except DomainError as exc:
        receipt["record_error"] = str(exc)
    atomic_write_json(path, receipt)
    return path


def _validation_envelope(job: AnalysisJob, receipt: dict,
                         expected_sources: list[str]) -> dict:
    """`validate_outputs`を呼び、例外は失敗の封筒へ畳む（成功へ変換しない）。"""
    try:
        return validate_outputs(job, receipt, expected_sources)
    except Exception as exc:  # 出力が壊れている可能性そのものを検証している
        return {"ok": False, "errors": ["OUTPUT_VALIDATION_FAILED"], "warnings": [],
                "primary_path": None, "sample_map": {}, "structure_errors": [repr(exc)]}


def _job_error_summary(receipt: dict, errors: list[str]) -> str:
    """completedにならなかった理由を1行にまとめる（jobのerror欄用）。"""
    head = f"termination={receipt['termination']} exit_code={receipt['exit_code']}"
    if not errors:
        return head
    return f"{head}: {','.join(errors[:_MAX_RECORDED])}"


def supervise(job_path: Path, *, command: list[str] | None = None,
              cancel_path: Path | None = None) -> dict:
    """Consoleを起動して終わりまで見張り、全終了経路で証跡と成果物を残す。

    同期呼出しも切り離しワーカーもこの関数を通る。段階の順序は
    **証跡の確定 → 収集 → 出力検証 → jobの保存**で固定で、各段階の成否を証跡へ
    書き足していく。逆順にすると収集や保存が落ちたときに終了の事実まで失われる。
    どの段階の例外も成功結果には変換しない。

    Parameters
    ----------
    job_path:
        analysis-job.json のパス。
    command:
        **内部試験用の注入口**。偽Consoleを実プロセスとして走らせるためだけに
        あり、MCPの公開引数にしてはいけない（任意コマンドの実行口になる）。
        省略時はjobのmethod / exe / save_projectから組む。
    cancel_path:
        存在したら取消要求とみなすファイル。

    Returns
    -------
    dict
        終了証跡（console-execution.v1）に段階記録`collection` `validation`
        `job_save`と入力比較`inputs`を足したもの。job保存に失敗した場合は
        `job_save.recovery`に復旧に必要なパスと意図した状態が入る。
    """
    job_path = Path(job_path)
    job = load_job(job_path)
    run_dir = Path(job.run_dir)
    dataset_root = Path(job.dataset_root)
    run_dir.mkdir(parents=True, exist_ok=True)
    msdial_out_dir = run_dir / MSDIAL_OUT_SUBDIR
    msdial_out_dir.mkdir(parents=True, exist_ok=True)

    state = read_supervision_state(run_dir)
    inputs = dict(state.get("inputs", {}))
    expected_sources = [str(p) for p in inputs.get("raw_inventory", [])]

    unavailable: list[str] = []
    method_sha256 = inputs.get("method_sha256") or _file_sha256(job.method_file)
    if not method_sha256:
        method_sha256 = _unavailable_hash("method")
        unavailable.append("method_sha256")
    exe_sha256 = inputs.get("exe_sha256") or _file_sha256(inputs.get("exe_path"))
    if not exe_sha256:
        exe_sha256 = _unavailable_hash("exe")
        unavailable.append("exe_sha256")

    # 二重ルートの実行前スナップショット。Consoleは -o にエクスポートを、
    # 生データフォルダに .pai2 / .dcl / .arf 等を出すので両方撮る。
    roots = {"run_dir": run_dir, "dataset_root": dataset_root}
    befores = {
        "run_dir": snapshot(run_dir),
        "dataset_root": snapshot(dataset_root, exclude_dir_names={RUNS_SUBDIR}),
    }
    start_fingerprint = raw_fingerprint(expected_sources)
    started_at = _utc_now()

    state["schema"] = SUPERVISION_SCHEMA
    state["inputs"] = inputs
    state["job_path"] = str(job_path)
    state["status"] = "running"
    state["started_at"] = started_at
    # 収集に失敗しても再収集できるよう、実行前スナップショットを残す。
    # 最終化前にこれを消してはいけない（消すと差分の起点が永久に失われる）。
    state["recovery"] = {"roots": {k: str(v) for k, v in roots.items()},
                         "befores": befores}
    atomic_write_json(supervision_state_path(run_dir), state)

    proc = None
    failure_note: str | None = None
    try:
        if command is None:
            command = _default_command(job, dataset_root, msdial_out_dir, inputs)
        command = [str(part) for part in command]
        proc = start_owned_process(command, cwd=run_dir,
                                   log_path=run_dir / CONSOLE_LOG_FILENAME)
    except (DomainError, OSError) as exc:
        failure_note = f"起動に失敗しました: {exc!r}"

    if proc is None:
        termination, exit_code, pid, identity = "launch_failed", None, None, None
    else:
        pid, identity = proc.pid, proc.identity
        try:
            # 起動情報の記録も try の内側に置く。start_owned_process が返った後・
            # try に入る前に例外が出ると、Job ハンドルを閉じられず孫が孤児になる。
            state["launch"] = {"command": command, "pid": pid,
                               "process_identity": identity}
            atomic_write_json(supervision_state_path(run_dir), state)
            termination, exit_code = _monitor(proc, timeout_s=job.timeout_s,
                                              cancel_path=cancel_path)
        except (DomainError, OSError) as exc:
            # 停止も回収もできなかった。成功にはせず「見失った」と記録する。
            termination, exit_code = "worker_lost", None
            failure_note = f"監視に失敗しました: {exc!r}"
        finally:
            # Jobハンドルを閉じるまでが所有権。ここを飛ばすと孫が孤児になる。
            proc.close()

    ended_at = _utc_now()
    end_fingerprint = raw_fingerprint(expected_sources)
    method_sha256_end = _file_sha256(job.method_file) or _unavailable_hash("method")
    changed: list[str] = []
    if expected_sources and end_fingerprint != start_fingerprint:
        changed.append("INPUT_CHANGED")
    if method_sha256_end != method_sha256:
        changed.append("METHOD_CHANGED")

    receipt: dict = {
        "schema": SCHEMA,
        "execution_id": f"exec_{uuid.uuid4().hex}",
        "job_id": job.job_id,
        "started_at": started_at,
        "ended_at": ended_at,
        "pid": pid,
        "process_identity": identity,
        "command_sha256": canonical_hash(command or []),
        "method_sha256": method_sha256,
        "exe_sha256": exe_sha256,
        "exit_code": exit_code,
        "termination": termination,
        "timeout_s": int(job.timeout_s),
        "collection": {"status": "pending"},
        "validation": {"status": "pending"},
        "job_save": {"status": "pending"},
        "inputs": {
            "raw_count": len(expected_sources),
            "fingerprint_start": start_fingerprint,
            "fingerprint_end": end_fingerprint,
            "method_sha256_end": method_sha256_end,
            "changed": changed,
            "unavailable": unavailable,
        },
    }
    if failure_note:
        receipt["error"] = failure_note
    # 終了の事実を先に確定させる。以降の段階が落ちてもここまでは残る。
    _write_receipt(run_dir, receipt)

    mztab_entries: list = []
    other_artifacts: list = []
    try:
        mztab_entries, other_artifacts = collect_artifacts(
            roots, befores,
            declared_polarity=job.polarity, declared_measure=job.measure)
        receipt["collection"] = {"status": "succeeded",
                                 "mztab_files": len(mztab_entries),
                                 "artifacts": len(other_artifacts)}
    except Exception as exc:
        receipt["collection"] = {
            "status": "failed", "error": repr(exc),
            "recovery": {"supervision_state_path": str(supervision_state_path(run_dir))},
        }
    _write_receipt(run_dir, receipt)

    if receipt["collection"]["status"] == "succeeded":
        job.primary_mztab_files = mztab_entries
        job.artifacts = other_artifacts
        validation = _validation_envelope(job, receipt, expected_sources)
        if not expected_sources:
            validation["errors"].append("INPUT_INVENTORY_MISSING")
        validation["errors"].extend(changed)
        validation["ok"] = not validation["errors"]
        receipt["validation"] = {
            "status": "succeeded" if validation["ok"] else "failed",
            "errors": validation["errors"],
            "warnings": validation["warnings"][:_MAX_RECORDED],
            "structure_errors": validation["structure_errors"][:_MAX_RECORDED],
            "primary_path": validation["primary_path"],
            "sample_count": len(validation["sample_map"]),
        }
        job.status = completion_status(
            receipt, validation, bool(mztab_entries or other_artifacts))
        job.error = (None if job.status == "completed"
                     else _job_error_summary(receipt, validation["errors"]))
    else:
        # 収集できていない以上、成果物の有無を語れない。completedへは進めない。
        receipt["validation"] = {"status": "skipped", "reason": "COLLECTION_FAILED",
                                 "errors": ["COLLECTION_FAILED"]}
        job.status = "failed"
        job.error = _job_error_summary(receipt, ["COLLECTION_FAILED"])

    try:
        save_job(job, job_path)
        receipt["job_save"] = {"status": "succeeded", "job_status": job.status}
    except Exception as exc:
        # 保存できなくても復旧できるだけの情報を応答へ残す。
        receipt["job_save"] = {
            "status": "failed", "error": repr(exc),
            "recovery": {
                "job_path": str(job_path),
                "intended_status": job.status,
                "receipt_path": str(receipt_path(run_dir)),
                "supervision_state_path": str(supervision_state_path(run_dir)),
            },
        }
    _write_receipt(run_dir, receipt)

    # 最終化。before スナップショットは消さずに残す（再収集の起点）。
    state["status"] = "finalized"
    state["execution_id"] = receipt["execution_id"]
    state["ended_at"] = ended_at
    atomic_write_json(supervision_state_path(run_dir), state)
    return receipt
