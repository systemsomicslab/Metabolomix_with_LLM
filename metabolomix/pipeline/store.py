"""pipeline-run（v3。旧 v1 は読むだけ）の永続化と受付冪等性（spec §9, §6・plan 2026-11-09）。

このモジュールが解く問題は4つある。

1. **schema**: 新しい run は `pipeline-run.v3`（極性ごとの単位 `units` を持つ）を書く。
   旧 `pipeline-run.v1`（下流の工程まで回した平らな形）は**読むだけ**で、保存は
   `PIPELINE_RUN_LEGACY` で拒む（読むなら `load_dataset(pipeline_path=…)`、Console を
   やり直すなら新しい `pipeline_run`）。未知の schema（実物の無い v2 を含む）は
   推測せず拒否し、既存 record の schema を書き換える保存も拒否する。
2. **永続化**: `pipeline-run.json` を原子的に保存・読込する。`state_revision` は
   楽観的並行制御用のカウンタで、要求自身が持つ `revision`（訂正ごとに増える別の
   カウンタ）とは意味も増分タイミングも異なる。
3. **受付冪等性**: 同一 source_root へ複数プロセス・複数 output_root から同時に
   `find_or_create_run` が呼ばれても、同じ入力・同じ要求なら二重に run を作らない。
   排他は `source_root` 単位（正規化した絶対パスの hash）に限定する。
4. **単位の読み替え**: 読む側（status・受領・データセットの入口）は `run_units` を通して
   極性ごとの `inputs` / `upstream` を読む。v1 の平らな record も 1 単位として読める。

依存グラフ上は `metabolomix.pipeline.inputs` と同じ層。`metabolomix.core.mcp_core` は
import しない（読取専用の `pipeline_status` が import だけで状態ディレクトリを作らない）。
索引の置き場所（既定 `%LOCALAPPDATA%`）はテストが `LIPIDMIX_PIPELINE_INDEX_DIR` で挿し替える。
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.core.atomic_io import (
    DomainError,
    atomic_write_json,
    canonical_hash,
    read_text_stable,
)
from metabolomix.core.process_control import file_lock, process_identity
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline.request import request_fingerprint

__all__ = [
    "ACTIVE_STATUSES",
    "LEGACY_SCHEMA",
    "SCHEMA",
    "SUPPORTED_SCHEMAS",
    "RUN_FILENAME",
    "read_result_data",
    "REQUESTS_SUBDIR",
    "TERMINAL_STATUSES",
    "UNIT_HANDLERS",
    "absolutize_unit_inputs",
    "relativize_unit_inputs",
    "create_run",
    "current_result_refs",
    "find_or_create_run",
    "initial_stage",
    "is_legacy",
    "legacy_refusal",
    "load_run",
    "run_schema_for",
    "run_units",
    "save_run",
    "stage_ids_for",
    "unit_stage_id",
    "unit_status",
    "verify_result_refs",
]

SCHEMA = "pipeline-run.v3"
LEGACY_SCHEMA = "pipeline-run.v1"

#: reader が受ける版。ここに無い schema は推測せず拒否する——未知の版を既知のものとして
#: 読むと、知らないフィールドを落としたまま保存し直して記録を壊す。v2 は実物が無いので読まない。
SUPPORTED_SCHEMAS = frozenset({LEGACY_SCHEMA, SCHEMA})

RUN_FILENAME = "pipeline-run.json"
REQUESTS_SUBDIR = "requests"
CONTROL_SUBDIR = "control"
RUNS_PARENT_SUBDIR = "runs"  # <source_root>/runs/pipeline_<UUID>/ （spec §4.3）

_INDEX_BASE_ENV = "LIPIDMIX_PIPELINE_INDEX_DIR"
_INDEX_DIRNAME = "pipeline-index"
_INDEX_FILENAME = "index.json"
_INDEX_LOCK_FILENAME = "index.lock"

#: spec §9.1: run の状態。active は running 側からしか動けない中間状態。
ACTIVE_STATUSES = frozenset({"planned", "running", "needs_input"})
TERMINAL_STATUSES = frozenset({"completed", "partial", "failed", "cancelled"})

#: 単位ごとに回る工程（この順）。最後に単位を持たない `receipt` が 1 つ付く
#: （計画の正準は `engine.build_stages`）。
UNIT_HANDLERS = ("prepare_input", "upstream", "validate_outputs")

#: run の記録に絶対パスのまま保存してよいのは、実際に単位のフォルダの外を指す場合だけ
#: （原本メソッド・LBM・実行体は通常その外にある）。中を指すときだけ相対化する対象。
_RELATIVIZABLE_INPUT_PATHS = (("method", "source_path"), ("lbm", "path"), ("exe", "path"))

_VANISHED = object()  # 索引entryはあるがpipeline_root/pipeline-run.jsonが読めない標識

#: find_or_create_runがrequest_idを既存runへ書き足す際、他アクターによる状態競合
#: (STATE_REVISION_CONFLICT)を有限回まで読み直して吸収する上限。
_REQUEST_ID_PATCH_MAX_ATTEMPTS = 5

#: 旧 v1 の run の代わりの道（`PIPELINE_RUN_LEGACY` の details）。
_LEGACY_REPLACEMENTS = {
    "read": "load_dataset(pipeline_path=…)",
    "rerun_console": "pipeline_run（新しい run）",
}


class _CorruptedReuse:
    """明示request_id無しで見つかった、成果物hashの検証に失敗したcompleted run。

    新runは作るが、この情報を新runのwarningsへ記録して破損を可視化する。
    """

    __slots__ = ("pipeline_root", "failures")

    def __init__(self, pipeline_root: Path, failures: list[dict]) -> None:
        self.pipeline_root = pipeline_root
        self.failures = failures


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------- 旧 v1 の扱い・単位の読み替え ----------

def is_legacy(record: dict) -> bool:
    """旧 `pipeline-run.v1`（読むだけの記録）か。"""
    return isinstance(record, dict) and record.get("schema") == LEGACY_SCHEMA


def legacy_refusal(pipeline_root: Path | str, action: str) -> DomainError:
    """旧 v1 の run を書き換える操作を断る `PIPELINE_RUN_LEGACY`（代わりの道つき）。"""
    return DomainError(
        "PIPELINE_RUN_LEGACY",
        f"旧形式（{LEGACY_SCHEMA}）の run は読むだけです（{action} はできません）。"
        "読むなら load_dataset(pipeline_path=…)、Console をやり直すなら新しい pipeline_run を使ってください。",
        {"pipeline_root": str(pipeline_root), "schema": LEGACY_SCHEMA, "action": action,
         "replacements": dict(_LEGACY_REPLACEMENTS)})


def _legacy_polarity(record: dict) -> str:
    polarity = ((record.get("inputs") or {}).get("polarity") or {})
    value = polarity.get("value") if isinstance(polarity, dict) else polarity
    return value if isinstance(value, str) and value else "unknown"


def run_units(record: dict) -> dict[str, dict]:
    """極性ごとの単位（`raw_dir` `unit_root` `inputs` `upstream` …）を返す。

    v3 は `record["units"]` そのもの。v1 は平らな `inputs` / `upstream` を、run の
    フォルダ直下を単位の根（`unit_root="."`）とする 1 単位として読み替える。
    読む側（status・受領・データセットの入口）はすべてここを通す。
    """
    if is_legacy(record):
        return {_legacy_polarity(record): {
            "raw_dir": ".", "unit_root": ".",
            "inputs": record.get("inputs") or {},
            "upstream": record.get("upstream") or {}}}
    return dict(record.get("units") or {})


def unit_stage_id(handler: str, polarity: str | None) -> str:
    """単位の工程の stage_id（`<handler>:<極性>`。単位を持たない工程は handler のまま）。"""
    return handler if polarity is None else f"{handler}:{polarity}"


def unit_status(record: dict, polarity: str) -> str:
    """単位の状態: `validated` / `failed` / `running` / `pending`。

    3 工程のどれかが failed・cancelled なら failed、全部 succeeded（か skipped）なら
    validated、どれかが running なら running、それ以外は pending。旧 v1 は接尾辞の無い
    stage_id を読む。
    """
    stages = record.get("stages") or {}
    legacy = is_legacy(record)
    statuses = []
    for handler in UNIT_HANDLERS:
        sid = handler if legacy else unit_stage_id(handler, polarity)
        statuses.append((stages.get(sid) or {}).get("status", "pending"))
    if any(s in ("failed", "cancelled") for s in statuses):
        return "failed"
    if all(s in ("succeeded", "skipped") for s in statuses):
        return "validated"
    if any(s == "running" for s in statuses):
        return "running"
    return "pending"


# ---------- stage初期化 ----------

def run_schema_for(request: dict) -> str:
    """要求のschemaから、このrunが書くpipeline-runのschemaを決める。

    新しい run は v3 の正規形の要求（`resolve_request` の戻り値）からしか作らない。
    それ以外は推測せず拒否する。
    """
    schema = request.get("schema") if isinstance(request, dict) else None
    if schema == request_mod.SCHEMA:
        return SCHEMA
    raise DomainError(
        "PIPELINE_REQUEST_INVALID",
        f"未知の要求schemaです: {schema!r}（既知: {request_mod.SCHEMA!r}）",
        {"schema": schema})


def stage_ids_for(request: dict) -> list[str]:
    """要求から record の stage ID 列を作る（計画の正準 `engine.build_stages` の stage_id）。"""
    # engine は store を import するので、ここは呼び出し時に引く（循環を避ける）。
    from metabolomix.pipeline import engine
    return [stage["stage_id"] for stage in engine.build_stages(request)]


def _initial_stage(stage_id: str) -> dict:
    """spec §9.1「stageの状態はpending/running/succeeded/skipped/needs_input/
    failed/cancelled」の初期値。"""
    return {
        "stage_id": stage_id,
        "status": "pending",
        "attempt": 0,
        "input_fingerprint": None,
        "result_refs": [],
        "warnings": [],
        "started_at": None,
        "updated_at": None,
        "error": None,
    }


def initial_stage(stage_id: str) -> dict:
    """`_initial_stage`の公開版（再開時に stage を初期値へ戻すのに使う）。"""
    return _initial_stage(stage_id)


# ---------- 絶対パスの相対化 ----------

def _relativize_if_under(value: str, pipeline_root: Path) -> str:
    """valueが絶対パスで、かつ実際にpipeline_root配下を指すときだけ相対化する。

    pipeline_root外を指す絶対パス（原本メソッド・LBM・実行体は通常ここに該当する
    ——source_rootや共有ライブラリフォルダはpipeline_rootの外にあるのが普通）は、
    無理に相対化すると`..`だらけの読みにくい文字列になるうえ意味も無いので、
    そのまま残す。
    """
    try:
        candidate = Path(value)
        if not candidate.is_absolute():
            return value
        rel = candidate.resolve().relative_to(Path(pipeline_root).resolve())
    except (OSError, ValueError):
        return value
    return str(rel).replace(os.sep, "/")


def _relativize_inputs_paths(inputs: dict, pipeline_root: Path) -> dict:
    """Task 9/13が絶対パスのまま持つフィールドのうち、実際にpipeline_root配下を
    指すものだけpipeline_root相対へ書き換えて返す（引数のdictは変更しない）。

    対象は次の4種類:
      - `method.source_path` / `lbm.path` / `exe.path`（Task13 `inspect_inputs`の
        出力。通常はpipeline_root外を指すので、その場合は絶対のまま残る）。
      - `manifest`（sample-manifest.v1の行リスト、Task 9 `parse_manifest`が返す）
        の各行の`source_file`。これはConsoleが実際に読んだ場所
        （mzTabの`ms_run[N]-location`と対応する`ds.assay_sources`）を指すため、
        Consoleの`dataset_root`が`pipeline_root/input`である以上、実際には
        pipeline_root配下を指すのが通常のケース。
    """
    result = copy.deepcopy(inputs)
    for top, sub in _RELATIVIZABLE_INPUT_PATHS:
        section = result.get(top)
        if isinstance(section, dict) and isinstance(section.get(sub), str):
            section[sub] = _relativize_if_under(section[sub], pipeline_root)
    manifest_rows = result.get("manifest")
    if isinstance(manifest_rows, list):
        for row in manifest_rows:
            if isinstance(row, dict) and isinstance(row.get("source_file"), str):
                row["source_file"] = _relativize_if_under(row["source_file"], pipeline_root)
    return result


def relativize_unit_inputs(inputs: dict, unit_root: Path) -> dict:
    """入力計画（配置後の snapshot を含む）の、単位の根の下を指すパスを相対にした写し。"""
    return _relativize_inputs_paths(inputs, unit_root)


def absolutize_unit_inputs(inputs: dict, unit_root: Path) -> dict:
    """`_relativize_inputs_paths` の逆。単位の根の相対に固定したパスを絶対パスへ戻す。

    `inputs.verify_inputs` などはこれらを絶対パスとして扱う。戻り値の
    `pipeline_root` は単位の根（`stage_inputs` に渡した場所）。
    """
    unit_root = Path(unit_root).resolve()

    def _fix(value: str) -> str:
        if isinstance(value, str) and not Path(value).is_absolute():
            return str((unit_root / value).resolve())
        return value

    result = copy.deepcopy(inputs)
    for top, sub in _RELATIVIZABLE_INPUT_PATHS:
        section = result.get(top)
        if isinstance(section, dict) and isinstance(section.get(sub), str):
            section[sub] = _fix(section[sub])
    manifest_rows = result.get("manifest")
    if isinstance(manifest_rows, list):
        for row in manifest_rows:
            if isinstance(row, dict) and isinstance(row.get("source_file"), str):
                row["source_file"] = _fix(row["source_file"])
    result["pipeline_root"] = str(unit_root)
    return result


# ---------- create_run / load_run / save_run ----------

def create_run(source_root: Path, request: dict, inputs: dict, *,
               warnings: list[dict] | None = None) -> Path:
    """新規pipeline_rootを採番して作成し、初期pipeline-run.jsonを書く。

    書込先は既定で`<source_root>/runs/pipeline_<UUID>/`、要求に`output_root`が
    あればその配下（spec §4.3「書込先は要求でリポジトリ外の別フォルダへ変更
    できる」）。pipeline_idはUUID4（秒精度の時刻だけに頼らない、spec §9.3）。
    ディレクトリ作成が衝突した場合（極めて起こりにくいが同一秒UUIDの理論上の
    再採番要求に備える）は新しいUUIDで採番し直す。

    `inputs` は `{"fingerprint": <受付の指紋>, "units": {<極性>: <入力計画>}}`。
    単位の順は要求の `units` の順で、各単位の根は `<pipeline_root>/<極性>`
    （入力計画の中の、その下を指す絶対パスは単位の根の相対にする）。

    `warnings` は、`find_or_create_run` が破損した completed run を黙って上書きした
    事実と、受付の警告（`request["intake_warnings"]`。読み飛ばした下流のキーなど）を
    新 run へ最初から記録するために使う。
    """
    source_root = Path(source_root).expanduser()
    if not source_root.is_dir():
        raise DomainError("DATASET_ROOT_NOT_FOUND",
                          f"source_rootが存在しません: {source_root}",
                          {"source_root": str(source_root)})
    if "fingerprint" not in inputs:
        raise DomainError("PIPELINE_INPUTS_INVALID",
                          "inputsにfingerprintがありません（受付冪等性の根拠に必須）。",
                          {"keys": sorted(inputs)})

    # 未知の要求schemaはここで止める。ディレクトリを掘る前に判定するので、
    # 拒否した要求が空のpipeline_rootを残さない。
    run_schema = run_schema_for(request)
    unit_plans = inputs.get("units") or {}
    polarities = [unit.get("polarity") for unit in request.get("units") or []]
    if not polarities or any(not p for p in polarities) or set(polarities) != set(unit_plans):
        raise DomainError(
            "PIPELINE_INPUTS_INVALID",
            "要求の単位（極性の確定したもの）と入力計画の単位が一致しません。",
            {"request_units": polarities, "planned_units": sorted(unit_plans)})
    stage_ids = stage_ids_for(request)

    output_root = request.get("output_root")
    base_dir = Path(output_root).expanduser() if output_root else source_root / RUNS_PARENT_SUBDIR

    pipeline_root: Path | None = None
    for _ in range(8):
        pipeline_id = uuid.uuid4().hex
        candidate = base_dir / f"pipeline_{pipeline_id}"
        try:
            candidate.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            continue
        pipeline_root = candidate
        break
    if pipeline_root is None:
        raise DomainError("PIPELINE_ID_COLLISION",
                          "pipeline_idの採番に失敗しました（再試行上限に達しました）。",
                          {"base_dir": str(base_dir)})

    now = _now_iso()
    request_rel = f"{REQUESTS_SUBDIR}/revision-0001.json"
    # 要求の保存は書換え不可（append-only）。以降のrevisionは新しいファイル名で
    # 追加する（Task 16の責務）。ここで書くのは最初の1件だけ。
    atomic_write_json(pipeline_root / request_rel, request)

    units = {}
    for unit in request["units"]:
        polarity = unit["polarity"]
        units[polarity] = {
            "raw_dir": unit.get("raw_dir") or ".",
            "unit_root": polarity,
            "inputs": _relativize_inputs_paths(unit_plans[polarity], pipeline_root / polarity),
            "upstream": {"console_job_path": None, "execution_id": None, "verification": None},
            "outputs": None,
        }
    seeded = list(warnings or []) + [dict(w) for w in request.get("intake_warnings") or []]

    record = {
        "schema": run_schema,
        "identity": {
            "pipeline_id": pipeline_id,
            "source_root": str(source_root.resolve()),
            "pipeline_root": str(pipeline_root.resolve()),
            "created_at": now,
            "updated_at": now,
        },
        "request": {
            "revision": 1,
            "request_id": None,  # find_or_create_runが必要なら後で埋める
            "content_hash": request_fingerprint(request),
            "saved_path": request_rel,
        },
        "status": "planned",
        "state_revision": 0,
        "stages": {sid: _initial_stage(sid) for sid in stage_ids},
        "units": units,
        "inputs": {"fingerprint": inputs["fingerprint"]},
        "results": [],
        "needs_input": None,
        "warnings": seeded,
        # spec §9.3「ロックにはowner identityを記録し」の実体。run_engine が owner lock
        # 取得直後にここへ自分のprocess identityを刻む。pidだけでは再利用を見分けられない
        # ため、read_status/prepare_resumeはcreation_timeと組みで`same_process`により
        # 生死を判定する（PIDだけを根拠にしない）。
        "worker": {"identity": None, "started_at": None},
    }
    atomic_write_json(pipeline_root / RUN_FILENAME, record)
    return pipeline_root


def load_run(path: Path) -> dict:
    """pipeline_root（`create_run`/`find_or_create_run`の戻り値）を受け取り、
    その配下の`pipeline-run.json`を読んで返す。

    読取りは`read_text_stable`——この1ファイルだけは、workerが`save_run`で
    置換している最中に監視（`pipeline_status`）が読みにくる。Windowsでは
    その一瞬に`PermissionError`（共有違反）が出るので、有界に再試行して
    吸収する（恒久的な権限エラーは窓を過ぎればそのまま伝わる）。
    """
    pipeline_root = Path(path)
    run_path = pipeline_root / RUN_FILENAME
    try:
        raw = read_text_stable(run_path)
    except FileNotFoundError as exc:
        raise DomainError("PIPELINE_RUN_NOT_FOUND",
                          f"pipeline-run.jsonが見つかりません: {run_path}",
                          {"pipeline_root": str(pipeline_root)}) from exc
    except OSError as exc:
        raise DomainError("PIPELINE_RUN_INVALID", f"pipeline-run.jsonを読めません: {exc}",
                          {"pipeline_root": str(pipeline_root)}) from exc
    try:
        data = json.loads(raw)
    except ValueError as exc:
        raise DomainError("PIPELINE_RUN_INVALID", f"pipeline-run.jsonが不正なJSONです: {exc}",
                          {"pipeline_root": str(pipeline_root)}) from exc
    # readerは既知の版だけを受け、読んだrecordのschemaのまま返す（spec §6）。
    if data.get("schema") not in SUPPORTED_SCHEMAS:
        raise DomainError("PIPELINE_RUN_INVALID",
                          f"schemaが不正です（既知 {sorted(SUPPORTED_SCHEMAS)}）: "
                          f"{data.get('schema')!r}",
                          {"pipeline_root": str(pipeline_root),
                           "schema": data.get("schema")})
    return data


def _assert_schema_unchanged(current_schema: object, new_schema: object,
                             pipeline_root: Path) -> None:
    """既存recordのschemaを書き換える保存を拒む（writer側のschema dispatch）。

    「runを別のschemaへ読み替えて保存し直す」はspec §6が明示的に禁じる
    （「v1 pipelineは旧契約・旧schemaで新規作成・再開する。過去の記録を
    インプレース変換しない」）。そのrunのstage列・結果IDの意味はv1の契約で
    確定しており、器だけ新しい版に付け替えても中身は付いてこない。
    """
    if new_schema == current_schema:
        return
    raise DomainError(
        "PIPELINE_RUN_SCHEMA_IMMUTABLE",
        f"既存recordのschemaを変更できません（現在 {current_schema!r} → "
        f"{new_schema!r}）。過去の記録はインプレース変換しません。",
        {"current_schema": current_schema, "new_schema": new_schema,
         "pipeline_root": str(pipeline_root)})


def _assert_results_append_only(current_results: list, new_results: list) -> None:
    """resultsは追記のみ許す（spec/brief「結果履歴は書換え不可」）。

    既存の各エントリと、新しいリストの同じ位置のエントリが完全一致することを
    要求する。件数が減る・途中の内容が変わるのはどちらも拒否する。
    """
    if len(new_results) < len(current_results):
        raise DomainError("PIPELINE_RESULTS_IMMUTABLE",
                          "resultsの件数が既存より減っています（履歴は書換え不可）。",
                          {"current": len(current_results), "new": len(new_results)})
    for index, existing in enumerate(current_results):
        if new_results[index] != existing:
            raise DomainError("PIPELINE_RESULTS_IMMUTABLE",
                              f"results[{index}]が既存の記録と一致しません（履歴は書換え不可）。",
                              {"index": index})


def save_run(path: Path, record: dict, *, expected_revision: int) -> None:
    """pipeline-run.jsonを楽観的並行制御つきで原子的に保存する。

    pipeline_root配下の`control/pipeline.lock`（OSの排他ロック）の中で、
    ディスク上の現在の`state_revision`が`expected_revision`と一致するかを
    照合してから書く。一致しなければ`STATE_REVISION_CONFLICT`（呼び出し側は
    最新を読み直してからやり直す）。書込み自体が失敗した場合（ディスク障害等）は
    例外がそのまま伝播し、ロックは解放されるが`pipeline-run.json`は直前の
    有効な内容のまま変化しない（`atomic_write_json`の性質）。
    """
    pipeline_root = Path(path)
    lock_path = pipeline_root / CONTROL_SUBDIR / "pipeline.lock"
    with file_lock(lock_path):
        current = load_run(pipeline_root)
        # 旧 v1 の記録は読むだけ（engine が 2 つの形を持ち続けないため）。
        if is_legacy(current):
            raise legacy_refusal(pipeline_root, "保存")
        if current.get("state_revision") != expected_revision:
            raise DomainError(
                "STATE_REVISION_CONFLICT",
                f"state_revisionが競合しています（期待 {expected_revision}、"
                f"現在 {current.get('state_revision')}）。最新を読み直してください。",
                {"expected_revision": expected_revision,
                 "current_revision": current.get("state_revision"),
                 "pipeline_root": str(pipeline_root)})
        _assert_schema_unchanged(current.get("schema"), record.get("schema"), pipeline_root)
        _assert_results_append_only(current.get("results", []), record.get("results", []))

        to_write = copy.deepcopy(record)
        to_write["state_revision"] = expected_revision + 1
        identity = dict(to_write.get("identity") or {})
        identity["updated_at"] = _now_iso()
        to_write["identity"] = identity
        atomic_write_json(pipeline_root / RUN_FILENAME, to_write)


# ---------- 受付索引（source_root単位の冪等性） ----------

def _index_base_dir() -> Path:
    override = os.environ.get(_INDEX_BASE_ENV)
    if override:
        return Path(override)
    local_appdata = os.environ.get("LOCALAPPDATA") or tempfile.gettempdir()
    return Path(local_appdata) / "Lipidmix" / _INDEX_DIRNAME


def _normalized_source_key(source_root: Path) -> str:
    """source_rootの比較用正規化文字列。Windowsは大小文字を同一視する。"""
    resolved = str(Path(source_root).resolve())
    if os.name == "nt":
        resolved = resolved.lower()
    return resolved


def _index_dir_for_source(source_root: Path) -> Path:
    """このsource_rootの受付索引を置くディレクトリ（output_rootには依存しない）。"""
    source_hash = canonical_hash({"source_root": _normalized_source_key(source_root)})
    return _index_base_dir() / source_hash


def _read_index(index_dir: Path) -> dict:
    path = index_dir / _INDEX_FILENAME
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {"entries": []}
    except (OSError, ValueError):
        # 索引はキャッシュ（真実の所在は各pipeline_rootのpipeline-run.json）。
        # 壊れていても受付自体を止めない——空として扱い、以後の書込みで復旧する。
        return {"entries": []}


def _write_index(index_dir: Path, data: dict) -> None:
    atomic_write_json(index_dir / _INDEX_FILENAME, data)


def _current_owner_identity() -> dict:
    identity = process_identity(os.getpid())
    return identity if identity is not None else {"pid": os.getpid(), "creation_time": None}


def _write_index_entry(index_dir: Path, *, pipeline_root: Path,
                       request_id: str | None, request_hash: str, input_hash: str) -> None:
    """索引へ1件をupsertする（run path/request_id/hash/ownerだけを持つ）。

    request_id指定時はrequest_idごとに高々1件（同じrequest_idの古い記録を
    置き換える）。request_id未指定時はpipeline_rootごとに高々1件（無名エントリの
    無限増殖を防ぐ）。解析ログ・行列はここに一切置かない。
    """
    data = _read_index(index_dir)
    entries = list(data.get("entries", []))

    def _is_same_slot(entry: dict) -> bool:
        if request_id is not None:
            return entry.get("request_id") == request_id
        return entry.get("request_id") is None and entry.get("pipeline_root") == str(pipeline_root)

    entries = [e for e in entries if not _is_same_slot(e)]
    entries.append({
        "pipeline_root": str(pipeline_root),
        "request_id": request_id,
        "request_hash": request_hash,
        "input_hash": input_hash,
        "owner": _current_owner_identity(),
        "created_at": _now_iso(),
    })
    _write_index(index_dir, {"entries": entries})


def read_result_data(pipeline_root: Path, results: list, output_name: str) -> list[dict]:
    """`output_name`の成果物を、参照hashを照合してから古い順に読む。

    worker の工程 handler と、run を対話セッションへ
    引き渡す経路（`metabolomix.pipeline.session_handoff`）が共有する。record の
    `results` がどう並び、data がどのキーに入るかを知る場所を1つに保つ。
    """
    pipeline_root = Path(pipeline_root)
    matches = [r for r in (results or [])
               if isinstance(r, dict) and r.get("output_name") == output_name]
    failures = verify_result_refs(pipeline_root, matches)
    if failures:
        raise DomainError(
            "RESULT_INTEGRITY_MISMATCH",
            f"保存済みの成果物が変更されています: {output_name}",
            {"output_name": output_name, "failures": failures})
    payloads = []
    for ref in matches:
        raw = (pipeline_root / ref["relative_path"]).read_text(encoding="utf-8")
        payloads.append(json.loads(raw).get("data"))
    return payloads


def verify_result_refs(pipeline_root: Path, result_refs: list) -> list[dict]:
    """`result_refs`各件の成果物hashを実ファイルと照合し、不一致だけを返す。

    戻り値は不一致（またはファイル自体を読めない）entryのリスト。空リストは全件一致
    （result_refsが空の場合も、検証対象が無いという意味で空虚に真とする）。
    `_artifacts_verify`（run全体のresults検証）とTask16
    `metabolomix.pipeline.engine.stage_inputs_unchanged`（stage単位の再利用可否判定）の
    両方がこれを共有する——「resultsの一部だけを検証したい」呼び出し元向けに、
    record全体を要求しない形で公開する。
    """
    failures: list[dict] = []
    for result in result_refs or []:
        rel = result.get("relative_path")
        expected = result.get("hash")
        if not rel or not expected:
            continue
        target = Path(pipeline_root) / rel
        try:
            actual = _sha256_file(target)
        except OSError:
            actual = None
        if actual != expected:
            failures.append({"relative_path": rel, "expected_hash": expected, "actual_hash": actual})
    return failures


def current_result_refs(record: dict) -> list[dict]:
    """`record["results"]`から「今そのrunが指している成果物」だけを取り出す。

    `results`は追記専用の履歴で、同じ`output_name`が複数revisionにまたがって
    並ぶ（正当な訂正→再開のたびに新しいrefが積まれる）。一方で成果物の
    書き出し先パスはoutput_nameごとに固定なので、再計算は同じファイルを
    上書きする——**古いrefのhashはもうどのファイルとも一致しない**。

    ここでの選び方（同じoutput_nameの最後の1件を「現在の成果物」とする）は、
    成果物を読む側と同じでなければならない。片方が最後の1件だけを見て、
    もう片方が履歴を全件検証すると、正当な訂正のあとで
    `RESULT_INTEGRITY_MISMATCH`——実際には壊れていないのに壊れたという報告——
    が出る（`find_or_create_run`の同一request_id再送がこれを踏んだ）。

    `output_name`を持たない古いrefは、代わりに`relative_path`で畳む
    （同じ出力先を指す最後の1件が現在の成果物、という規則は同じ）。
    """
    current: dict = {}
    for result in record.get("results") or []:
        if not isinstance(result, dict):
            continue
        key = result.get("output_name") or result.get("relative_path") or result.get("result_id")
        current[key] = result
    return list(current.values())


def _artifacts_verify(record: dict, pipeline_root: Path) -> list[dict]:
    """record["results"]が今指している成果物のhashを実ファイルと照合する。

    戻り値は不一致（またはファイル自体を読めない）entryのリスト。空リストは
    全件一致（resultsが空——探索・比較のいずれもまだ結果を持たないcompleted、
    通常は起こらないが安全側で許す——場合も同様に空虚に真とする）。

    検証対象は履歴全件ではなく`current_result_refs`が選ぶ現行ref
    （上書きされた旧revisionの成果物を「壊れている」と読まないため）。
    """
    return verify_result_refs(pipeline_root, current_result_refs(record))


def _resolve_reusable(pipeline_root_str: str, *, explicit_request_id: bool):
    """索引entryが指すpipeline_rootを、そのまま返してよいか判定する。

    戻り値:
        Path             そのまま返してよい（活動中／completedで成果物検証OK／
                         失敗・取消・部分完了などそのままの状態で返してよい終端状態）。
        _VANISHED        pipeline_root自体（pipeline-run.json）が読めない
                         （索引はあるがrunが消失している）。呼び出し側はこの索引
                         entryを無視して次の候補・新規作成へフォールバックしてよい。
        _CorruptedReuse  明示request_id無しで見つかったcompletedだが成果物hashが
                         壊れているrun。呼び出し側は消失と同様フォールバックする
                         が、破損の事実を新runの警告として残さなければならない
                         （controller ruling R17）。

    例外:
        DomainError("RESULT_INTEGRITY_MISMATCH")
            明示request_idが一致する既存runがcompletedだが、記録済み成果物hashが
            実ファイルと一致しない（壊れている）場合。新runを黙って作らない。
    """
    pipeline_root = Path(pipeline_root_str)
    try:
        record = load_run(pipeline_root)
    except DomainError:
        return _VANISHED
    if is_legacy(record):
        return _VANISHED  # 旧 v1 の run は新しい受付で再利用しない（指紋も別）
    status = record.get("status")
    if status == "completed":
        failures = _artifacts_verify(record, pipeline_root)
        if failures:
            if explicit_request_id:
                raise DomainError(
                    "RESULT_INTEGRITY_MISMATCH",
                    f"request_idに一致する既存runの成果物を検証できませんでした: {pipeline_root}",
                    {"pipeline_root": str(pipeline_root), "failures": failures})
            # 内容一致だけのフォールバックでは黙って再利用しない。ただし静かに
            # 上書きもしない——呼び出し側（find_or_create_run）が新runの
            # warningsへ書き足せるよう、破損情報を運ぶ。
            return _CorruptedReuse(pipeline_root, failures)
    return pipeline_root


def _lookup_request(index_dir: Path, *, request_id: str | None,
                    request_hash: str, input_hash: str,
                    corrupted: list) -> Path | None:
    """索引から再利用可能なpipeline_pathを探す。無ければNone（create_runへ進む）。

    `corrupted`は呼び出し側（find_or_create_run）が用意する出力用リスト。
    明示request_id無しのフォールバック探索中に見つかった破損completed run
    （`_CorruptedReuse`）をここへ積む——新runを作ることになった場合、
    呼び出し側がこれを警告として書き足す。
    """
    data = _read_index(index_dir)
    entries = data.get("entries", [])

    if request_id is not None:
        for entry in entries:
            if entry.get("request_id") != request_id:
                continue
            if entry.get("request_hash") != request_hash or entry.get("input_hash") != input_hash:
                raise DomainError(
                    "IDEMPOTENCY_CONFLICT",
                    f"request_id={request_id!r}は既存の異なる内容の要求と衝突しています。",
                    {"request_id": request_id, "pipeline_root": entry.get("pipeline_root")})
            resolved = _resolve_reusable(entry["pipeline_root"], explicit_request_id=True)
            if resolved is _VANISHED:
                break  # このrequest_idの記録は消失。素通りしてcontentフォールバックへ
            # explicit_request_id=Trueは_CorruptedReuseを返さない（例外で止まる）。
            return resolved

    for entry in entries:
        if entry.get("request_hash") != request_hash or entry.get("input_hash") != input_hash:
            continue
        resolved = _resolve_reusable(entry["pipeline_root"], explicit_request_id=False)
        if resolved is _VANISHED:
            continue
        if isinstance(resolved, _CorruptedReuse):
            corrupted.append(resolved)
            continue
        return resolved
    return None


def _corruption_warning(corrupted: _CorruptedReuse) -> dict:
    """R17: 破損した既存completed runを黙って上書きした事実を新runの
    warningsへ記録する（spec §9.1「warnings: 機械可読code、対象工程、説明」）。
    stage非依存の事象なので対象工程はNone。
    """
    broken = ", ".join(f["relative_path"] for f in corrupted.failures)
    return {
        "code": "STALE_RUN_ARTIFACT_CORRUPTED",
        "stage_id": None,
        "message": (
            f"既存run {corrupted.pipeline_root} は成果物ハッシュの検証に失敗した"
            f"ため再利用せず、新しいrunを作成しました。破損した成果物: {broken}"
        ),
    }


def _patch_request_id_with_retry(pipeline_path: Path, request_id: str) -> None:
    """既存run（新規作成直後含む）の`request.request_id`を書き足す。

    ここは`source_root`単位の`index.lock`では保護されるが、その run自身の
    `control/pipeline.lock`では保護されない一瞬の読取り（`load_run`）を挟む。
    その間にConsoleワーカーなど別アクターがこのrunの`state_revision`を進める
    （例えばstages更新の`save_run`）と、続く`save_run`の楽観チェックが
    `STATE_REVISION_CONFLICT`を投げる——受付自体は正当（runは既に存在し
    再利用してよい）のに、無関係な同時更新のせいで受付全体を失敗させては
    いけない（レビュー finding 1）。

    pipeline自身のlockを先に取ってから読み書きする案は採らなかった:
    `save_run`が内部で同じ`control/pipeline.lock`を取るため、同一プロセスから
    二重に取得することになり、`file_lock`は同一ファイルへの新しい`os.open`ごとに
    別のロック要求として扱われる（POSIXの`flock`もWindowsの`msvcrt.locking`も
    ファイル記述子単位）——自己デッドロックの危険がある。`save_run`を
    「ロック取得済み」モード対応に改修する手もあるが、ロック解放漏れの経路を
    増やすだけの割に、有限回のリトライで十分安全に吸収できる。
    """
    for _ in range(_REQUEST_ID_PATCH_MAX_ATTEMPTS):
        record = load_run(pipeline_path)
        if record["request"].get("request_id") == request_id:
            return
        record["request"]["request_id"] = request_id
        try:
            save_run(pipeline_path, record, expected_revision=record["state_revision"])
        except DomainError as exc:
            if exc.code == "STATE_REVISION_CONFLICT":
                continue  # 他アクターが進めた最新を読み直して再試行
            raise
        else:
            return
    raise DomainError(
        "STATE_REVISION_CONFLICT",
        f"request_idの記録が状態競合の再試行上限（{_REQUEST_ID_PATCH_MAX_ATTEMPTS}回）"
        f"に達しました: {pipeline_path}",
        {"pipeline_root": str(pipeline_path), "request_id": request_id})


def find_or_create_run(source_root: Path, request: dict, inputs: dict, *,
                       request_id: str | None = None) -> Path:
    """同一source_root・同一要求・同一入力の受付を冪等にする（spec §9.3）。

    排他はsource_root単位（`output_root`には依存しない）。同一input fingerprint・
    同一要求の活動中runがあればそれを返す。成果物hashの一致するcompleted runも
    再利用する。request_idを指定した場合、既存の異なる内容の要求と衝突すれば
    `IDEMPOTENCY_CONFLICT`、一致した既存completed runの成果物が壊れていれば
    `RESULT_INTEGRITY_MISMATCH`（新runを黙って作らない）。
    """
    source_root = Path(source_root).expanduser()
    if "fingerprint" not in inputs:
        raise DomainError("PIPELINE_INPUTS_INVALID",
                          "inputsにfingerprintがありません（受付冪等性の根拠に必須）。",
                          {"keys": sorted(inputs)})

    request_hash = request_fingerprint(request)
    input_hash = inputs["fingerprint"]

    index_dir = _index_dir_for_source(source_root)
    lock_path = index_dir / _INDEX_LOCK_FILENAME
    with file_lock(lock_path):
        corrupted: list[_CorruptedReuse] = []
        pipeline_path = _lookup_request(index_dir, request_id=request_id,
                                        request_hash=request_hash, input_hash=input_hash,
                                        corrupted=corrupted)
        if pipeline_path is None:
            seed_warnings = [_corruption_warning(c) for c in corrupted]
            pipeline_path = create_run(source_root, request, inputs, warnings=seed_warnings)

        if request_id is not None:
            _patch_request_id_with_retry(pipeline_path, request_id)

        _write_index_entry(index_dir, pipeline_root=pipeline_path, request_id=request_id,
                           request_hash=request_hash, input_hash=input_hash)
        return pipeline_path
