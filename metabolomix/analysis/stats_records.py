"""統計の工程の記録（1 工程 1 ファイル、解析フォルダ配下の `reports/stats/`）。

前処理・比較・PCA の対話ツールが、何をどの設定で計算したかを残す。記録は追記専用
（`atomic_io.write_new_json`。同名があれば `STATS_RECORD_EXISTS`）で、リポジトリ側には書かない。
数値は `json_safe` を通してから `allow_nan=False` で書く（無限大は `"inf"`、NaN は null。
厳密な JSON の読み手が止まらない）。定量の記録（`quant/records.py`）と同じ形の別系統。
"""
from __future__ import annotations

import uuid
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.core.atomic_io import DomainError, write_new_json
from metabolomix.core.serialization import json_safe

__all__ = ["RECORD_SCHEMA", "STEPS", "record_dir", "write_stats_record"]

RECORD_SCHEMA = "lipidmix.stats-record.v1"
STEPS = ("preprocess", "compare", "pca")


def record_dir(reports_dir: Path) -> Path:
    """記録の置き場所。`reports_dir` は記録の置き場所の基準（例 `<解析フォルダ>/reports`）。"""
    return Path(reports_dir) / "stats"


def _new_record_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"s_{stamp}_{uuid.uuid4().hex[:8]}"


def write_stats_record(reports_dir: Path, step: str, *, inputs: dict, parent_ids: list[str],
                       settings: dict, summary: dict) -> dict:
    """記録を 1 件 `<reports_dir>/stats/<record_id>.json` に書き、その内容（と `path`）を返す。

    `reports_dir` は解析フォルダ配下の `reports/`（解析フォルダそのものではない。
    `mcp_core._resolve_report_dir(holder)` の戻り値を渡す）。
    `inputs` は入力の識別（データの指紋・行列 ID など）、`parent_ids` は前の工程の記録 ID、
    `settings` は効いた設定、`summary` は要約。`path` は返すだけでファイルには書かない
    （解析フォルダを動かしても記録の中身は変わらない）。
    """
    if step not in STEPS:
        raise DomainError("STATS_RECORD_STEP_INVALID", f"未知の工程: {step}",
                          {"step": step, "steps": list(STEPS)})
    record_id = _new_record_id()
    record = json_safe({
        "schema": RECORD_SCHEMA,
        "record_id": record_id,
        "step": step,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": dict(inputs),
        "parent_ids": list(parent_ids),
        "settings": settings,
        "summary": summary,
    })
    path = record_dir(reports_dir) / f"{record_id}.json"
    try:
        write_new_json(path, record)
    except FileExistsError:
        raise DomainError("STATS_RECORD_EXISTS",
                          f"同名の記録が既にある（上書きしない）: {path.name}",
                          {"record_id": record_id}) from None
    return {**record, "path": str(path)}
