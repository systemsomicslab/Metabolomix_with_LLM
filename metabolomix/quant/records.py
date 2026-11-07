"""定量の工程の記録（1 工程 1 ファイル、解析フォルダ配下の `reports/quant/`）。

記録は追記専用: 上書きせず、入力が変わったら後ろの記録を「古い」と示すだけで消さない。
やり直しは新しい記録になる。リポジトリ側には書かない。
"""
from __future__ import annotations

import json
import os
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.core.atomic_io import DomainError, read_text_stable

RECORD_SCHEMA = "lipidmix.quant-record.v1"
STEPS = ("normalize", "merge", "select")
AUDIT_SUFFIX = ".audit.json"


def record_dir(folder: Path) -> Path:
    return Path(folder) / "reports" / "quant"


def _new_record_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"q_{stamp}_{uuid.uuid4().hex[:8]}"


def _create_exclusive(path: Path, payload: dict) -> None:
    """一時ファイルへ書いてから `os.link` で本名にする（既存があれば上書きせず失敗）。

    `os.replace` は上書きしてしまうので使えない。`os.link` は宛先が在れば
    FileExistsError になり、書き途中の本名ファイルが読まれることも無い。
    """
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.",
            suffix=".tmp", delete=False,
        ) as tmp:
            tmp_path = Path(tmp.name)
            json.dump(payload, tmp, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
            tmp.flush()
            os.fsync(tmp.fileno())
        try:
            os.link(tmp_path, path)
        except FileExistsError:
            raise DomainError("QUANT_RECORD_EXISTS",
                              f"同名の記録が既にある（上書きしない）: {path.name}") from None
    finally:
        if tmp_path is not None and tmp_path.exists():
            os.unlink(tmp_path)


def write_record(
    folder: Path, step: str, *, inputs: dict[str, str], parent_ids: list[str],
    settings: dict, summary: dict, skipped: bool = False,
) -> dict:
    """記録を 1 件書いて、その内容を返す。`inputs` はパス → sha256。"""
    if step not in STEPS:
        raise DomainError("QUANT_RECORD_STEP_INVALID", f"未知の工程: {step}")
    record_id = _new_record_id()
    record = {
        "schema": RECORD_SCHEMA,
        "record_id": record_id,
        "step": step,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": dict(inputs),
        "parent_ids": list(parent_ids),
        "settings": settings,
        "summary": summary,
        "skipped": skipped,
        "rules_version": settings.get("rules_version") if isinstance(settings, dict) else None,
    }
    target_dir = record_dir(folder)
    target_dir.mkdir(parents=True, exist_ok=True)
    _create_exclusive(target_dir / f"{record_id}.json", record)
    return record


def write_audit(folder: Path, record_id: str, payload: dict) -> Path:
    """記録の横の `<record_id>.audit.json`（上書きしない。`read_records` は記録として読まない）。"""
    target = record_dir(folder)
    target.mkdir(parents=True, exist_ok=True)
    path = target / f"{record_id}{AUDIT_SUFFIX}"
    _create_exclusive(path, payload)
    return path


def read_records(folder: Path) -> tuple[list[dict], list[str]]:
    """(作った順の記録, 除いたファイルを挙げる警告)。schema 違い・読めないものは除く。"""
    directory = record_dir(folder)
    if not directory.is_dir():
        return [], []
    found: list[dict] = []
    warnings: list[str] = []
    for path in sorted(directory.glob("*.json")):
        if path.name.endswith(AUDIT_SUFFIX):  # 記録の横の audit（記録ではない）
            continue
        try:
            data = json.loads(read_text_stable(path))
        except (OSError, ValueError) as exc:
            warnings.append(f"{path.name}: 読めないので除いた（{type(exc).__name__}）")
            continue
        if not isinstance(data, dict) or data.get("schema") != RECORD_SCHEMA:
            warnings.append(f"{path.name}: schema が {RECORD_SCHEMA} ではないので除いた")
            continue
        found.append(data)
    found.sort(key=lambda r: (str(r.get("created_at", "")), str(r.get("record_id", ""))))
    return found, warnings


def is_stale(
    record: dict, current_inputs: dict[str, str], known_ids: set[str],
    *, stale_ids: set[str] | frozenset[str] = frozenset(),
) -> str | None:
    """古ければ理由、新しければ None。

    `known_ids` は「在って古くない」記録 ID、`stale_ids` は「在るが古い」記録 ID。
    呼び出し側が連鎖を作った順にたどって渡す。記録にある入力だけを比べ、
    現在の入力に無いパスは変わったものと数える。
    """
    for path, sha in record.get("inputs", {}).items():
        if current_inputs.get(path) != sha:
            return f"input_changed:{Path(path).name}"
    for parent in record.get("parent_ids", []):
        if parent in known_ids:
            continue
        if parent in stale_ids:
            return f"parent_stale:{parent}"
        return f"parent_missing:{parent}"
    return None


def merged_record_folder(folder_a: Path, folder_b: Path) -> Path:
    """統合後の記録の置き場所＝2 つの解析フォルダの共通の親。

    ドライブが違う・同じ場所・どちらかが祖先なら置き場所を決められないので止める。
    """
    a = Path(folder_a).resolve()
    b = Path(folder_b).resolve()
    if os.path.normcase(a.drive) != os.path.normcase(b.drive):
        raise DomainError("QUANT_RECORD_LOCATION_REQUIRED",
                          "2 つの解析フォルダのドライブが違う。記録の置き場所を指定してほしい。")
    if _is_within(a, b) or _is_within(b, a):
        raise DomainError("QUANT_RECORD_LOCATION_REQUIRED",
                          "一方がもう一方の中にある。記録の置き場所を指定してほしい。")
    return Path(os.path.commonpath([a, b]))


def _is_within(child: Path, parent: Path) -> bool:
    c = os.path.normcase(str(child))
    p = os.path.normcase(str(parent))
    return c == p or c.startswith(p.rstrip("\\/") + os.sep)
