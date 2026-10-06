"""**旧方式**の切り離し実行が残したサイドカー（読むためだけに残す）。

かつて `console_run(detach=True)` は Console を起動して即座に戻り、完了時の
生成物収集を別の呼び出し（`console_status`）に引き継がせていた。その引き継ぎに
使ったのがこのファイル（run_dir の `.detached-state.json`）で、pid と実行前
スナップショットだけを持つ。

**新しい実行はこれを書かない**。いまは監視ワーカー（`metabolomix.console.worker`）が
最後まで見張り、終了証跡（`execution-result.json`）と成果物を自分で残す。旧 state に
無いのは終了コードと process identity で、それが無いと「もう走っていない」ことしか
分からない——成功したかどうかは分からない。ファイルが増えた事実だけで完了と
判定すると、途中で落ちた実行が completed に化ける。

そのため残っている旧 state は `console_status` が `EXECUTION_UNRESOLVED` として
報告するだけで、成果物の収集も状態の昇格も行わない。`write_detached_state` は
その状況を再現するテストのためだけに残している。

置き場所が `analysis-job.json` ではなく run_dir のサイドカーなのは、
`analysis-job.v2` が別リポジトリ（massbank-context）との契約でもあるため。
実行方式という内部事情でスキーマを増やさない。
"""
from __future__ import annotations

import json
from pathlib import Path

STATE_FILENAME = ".detached-state.json"


def state_path(run_dir: Path) -> Path:
    return Path(run_dir) / STATE_FILENAME


def write_detached_state(run_dir: Path, pid: int, befores: dict[str, dict]) -> Path:
    """pid と実行前スナップショットを run_dir に残す。"""
    path = state_path(run_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"pid": int(pid), "befores": befores}, ensure_ascii=False),
        encoding="utf-8")
    return path


def read_detached_state(run_dir: Path) -> dict | None:
    """引き継ぎ状態を返す。無い・壊れているなら None。"""
    path = state_path(run_dir)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or "pid" not in data:
        return None
    data.setdefault("befores", {})
    return data


def clear_detached_state(run_dir: Path) -> None:
    """収集を終えたら消す。残すと 2 回目の console_status がもう一度収集し、
    確定済みのジョブを上書きしてしまう。"""
    try:
        state_path(run_dir).unlink()
    except OSError:
        pass
