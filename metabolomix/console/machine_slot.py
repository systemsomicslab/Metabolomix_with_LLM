"""MS-DIAL Console を機械全体で同時に 1 本にする排他と、起こす前の commit の空きの確認。

Console は研究室の参照ライブラリ（MSP・LBM）を読む脂質の run で 20 GB 以上の commit を確保する。
極性ごとに別の run を同時に起こすと Console が 2 本並び、どちらも遅くなるうえ、commit の上限で
Windows のエラー 1450（`ERROR_NO_SYSTEM_RESOURCES`）で落ちうる（reanalysis-study の記録。
1 本の run の中では極性を 1 つずつ回している）。そこで run を跨いで Console を 1 本にする:

- **排他**は OS のファイルロック（`process_control.wait_for_lock`）。持っている worker が
  死ねば OS が外す（その worker の Console も Job ごと止まる）ので、取り残しの掃除は要らない。
  置き場所は設定ファイルのフォルダの `data/.msdial-console.lock`（参照ライブラリのキャッシュの
  既定と同じ規則）。この機械のサーバ・worktree・worker は同じ設定ファイル（main ツリーの
  `lipidmix.local.toml`）を指すので、同じロックを取り合う。テストは設定ファイルをテストごとの
  tmp へ向けるので、本物の run と取り合わない。
- **commit の空き**は止めずに測るだけ。下限（設定 `[msdial] min_free_commit_gb`、既定 8。
  0 で測らない）を下回れば警告を返す。

依存は stdlib と `metabolomix.core` の leaf（`atomic_io`・`process_control`・`user_config`）だけ
（pipeline の worker が読む）。
"""
from __future__ import annotations

import contextlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.core import process_control, user_config
from metabolomix.core.atomic_io import atomic_write_json

__all__ = [
    "DEFAULT_MIN_FREE_COMMIT_GB",
    "LOW_COMMIT_WARNING_CODE",
    "MIN_FREE_COMMIT_SETTING_KEY",
    "check_commit_memory",
    "console_lock_path",
    "console_slot",
    "read_owner",
]

LOCK_FILENAME = ".msdial-console.lock"
#: ロックを持っている run の申告（参考。ロックの正体は OS が握る状態で、このファイルではない）。
OWNER_FILENAME = ".msdial-console.owner.json"

MIN_FREE_COMMIT_SETTING_KEY = "msdial.min_free_commit_gb"
DEFAULT_MIN_FREE_COMMIT_GB = 8
LOW_COMMIT_WARNING_CODE = "LOW_COMMIT_MEMORY"

#: 順番待ちの間に取消を見る間隔（秒）。
POLL_INTERVAL_S = 2.0

_GB = 1024 ** 3


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def console_lock_path() -> Path:
    """機械全体の Console の排他のロックファイル（設定ファイルのフォルダの `data/` の下）。"""
    return user_config.config_file_path().parent / "data" / LOCK_FILENAME


def _owner_path() -> Path:
    return console_lock_path().with_name(OWNER_FILENAME)


def read_owner() -> dict | None:
    """いまロックを持っている run の申告（読めなければ None）。参考の情報で、古いことがある。"""
    try:
        data = json.loads(_owner_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


@contextlib.contextmanager
def console_slot(*, owner: dict, on_wait=None, should_abort=None,
                 poll_interval: float = POLL_INTERVAL_S):
    """機械全体で 1 つの Console の順番を取る（取れるまで待つ）。yield するのは待った秒数。

    `on_wait()` は順番待ちが始まったときに 1 回だけ呼ぶ。`should_abort()` が真になったら
    `LOCK_WAIT_ABORTED`（`process_control.wait_for_lock`）。取れたら `owner`（run と単位）を
    申告ファイルに書く。
    """
    with process_control.wait_for_lock(console_lock_path(), poll_interval=poll_interval,
                                       on_wait=on_wait, should_abort=should_abort) as waited:
        with contextlib.suppress(OSError):
            atomic_write_json(_owner_path(), {**owner, "pid": os.getpid(), "since": _now_iso()})
        yield waited


def _min_free_commit_gb() -> tuple[int, str | None]:
    """下限（GB）と、設定を読めなかったときの注記。"""
    try:
        setting = user_config.get_setting(MIN_FREE_COMMIT_SETTING_KEY)
    except user_config.ConfigInvalidError as exc:
        return DEFAULT_MIN_FREE_COMMIT_GB, f"設定ファイルを読めないので既定の {DEFAULT_MIN_FREE_COMMIT_GB} GB で測った（{exc.message}）"
    if setting is None:
        return DEFAULT_MIN_FREE_COMMIT_GB, None
    return max(0, int(setting.value)), None


def check_commit_memory() -> tuple[dict | None, dict | None]:
    """commit の空きを測り、`(測った値, 警告)` を返す。測れなければ `(None, None)`。

    下限が 0 なら測らない。警告は工程の警告（`{"code", "message", "details"}`）の形。
    """
    threshold_gb, note = _min_free_commit_gb()
    if threshold_gb <= 0:
        return None, None
    status = process_control.commit_memory_status()
    if status is None:
        return None, None
    available_gb = round(status["available_commit_bytes"] / _GB, 1)
    measured = {"available_commit_gb": available_gb,
                "commit_limit_gb": round(status["commit_limit_bytes"] / _GB, 1),
                "min_free_commit_gb": threshold_gb}
    if note:
        measured["note"] = note
    if status["available_commit_bytes"] >= threshold_gb * _GB:
        return measured, None
    warning = {
        "code": LOW_COMMIT_WARNING_CODE,
        "message": (f"Console を起こす前の commit の空き（物理メモリとページファイルのうち、まだ確保されて"
                    f"いない量）が {available_gb} GB で、下限 {threshold_gb} GB"
                    f"（{user_config.config_key(MIN_FREE_COMMIT_SETTING_KEY)}）を下回っていた。Console は"
                    "そのまま起こした。研究室の MSP・LBM を読む脂質の run は 20 GB 以上を確保し、足りないと"
                    "Windows のエラー 1450 で落ちることがある。落ちたら、ほかの重いプロセスを止めるか"
                    "ページファイルを広げてから rerun_upstream=true で再開する。"),
        "details": measured,
    }
    return measured, warning
