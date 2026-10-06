"""単体 Console 実行の切り離しワーカーと、実行の所有権。

`console_run` の同期呼出しも `detach=True` も、**このモジュールの `run_job` を
通る**。run_job は job 単位の OS ロックを取り、`execution.supervise` を最後まで
回し、証跡から読めない人向けの所見（極性の裏取り・未対応 mzTab・サンプル別
ファイルの欠落）をジョブへ書き足すところまでを 1 か所で行う。同期と非同期で
結果が違ってはいけないので、経路を 2 本持たない。

旧経路は「起動して pid を返し、完了は `console_status` が後から拾う」形だった。
誰も監視していない時間があり、`console_status` が呼ばれなければ生成物は永久に
ジョブへ載らず、終了コードも終了理由も回収できなかった。いまは切り離した
ワーカー自身が最後まで見張るので、`console_status` は保存済みの状態を読むだけの
読取専用ツールになる。

**所有権**（`worker.json` の `owner`）は「このジョブを今まさに実行している
プロセス」の記録。生きているかどうかは pid ではなく process identity で見る
（pid は再利用される）。`console_cleanup` はこれが生きている間、生成物を消さない
——監視ワーカーが書いている最中のファイルを足元から抜くことになるため。

CLI:

    python -m metabolomix.console.worker --job <analysis-job.json>

`--console-arg` は**試験用の注入口**で、偽 Console を実プロセスとして走らせる
ためだけにある（`supervise(command=...)` と同じ位置付け）。MCP からは到達
できない——`console_run` はこれを渡さない。
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.console.execution import (
    read_supervision_state,
    supervise,
    supervision_state_path,
)
from metabolomix.console.job_manager import load_job, save_job, update_status
from metabolomix.console.output_collector import (
    meta_conflict_warnings,
    missing_per_sample_output_warnings,
    polarity_crosscheck_warnings,
    record_mztab_provenance,
    unsupported_mztab_warnings,
)
from metabolomix.core.atomic_io import DomainError, atomic_write_json
from metabolomix.core.process_control import file_lock, launch_detached, same_process

__all__ = ["LOCK_FILENAME", "WORKER_LOG_FILENAME", "annotate_job",
           "clear_owner", "launch_console_worker", "lock_path", "owner_is_active",
           "owner_summary", "read_owner", "run_job", "write_owner"]

#: job 単位の排他ロック。実体は OS が握る byte-range lock で、このファイルの
#: 存在や中身ではない（消してもロックは解けないし、残っていても解けている）。
LOCK_FILENAME = "job.lock"

#: 切り離しワーカー自身の stdout/stderr。Console のログ（msdial.log）とは別物で、
#: ワーカーが起動直後に落ちた場合の唯一の手掛かりになる。
WORKER_LOG_FILENAME = "worker.log"

#: 所有者の種別。pipeline が所有する場合の値は Task 14 で足す。
OWNER_KIND_CONSOLE = "console_worker"

#: ロック取得の待ち時間（秒）。ここを超えるなら別のプロセスが本当に走っている。
LOCK_TIMEOUT_S = 2.0

#: 「もう実行し終えた」と読むジョブ状態。ロック取得後の再確認（`run_job`）で
#: これに当たれば MS-DIAL は起動しない。`planned`（未実行）と`running`
#: （切り離し起動の親が先に印を付けた状態）だけが実行してよい入口。
_FINISHED_STATUSES = frozenset({"completed", "partial", "failed", "cleaned"})


def lock_path(run_dir: Path) -> Path:
    return Path(run_dir) / LOCK_FILENAME


def worker_log_path(run_dir: Path) -> Path:
    return Path(run_dir) / WORKER_LOG_FILENAME


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------- 所有権 ----------

def read_owner(run_dir: Path) -> dict | None:
    """このジョブの実行を所有しているプロセスの記録を返す。無ければ None。"""
    owner = read_supervision_state(run_dir).get("owner")
    return owner if isinstance(owner, dict) else None


def write_owner(run_dir: Path, owner: dict) -> None:
    """所有者を `worker.json` へ併合して書く（他のキーは残す）。"""
    state = read_supervision_state(run_dir)
    state["owner"] = owner
    atomic_write_json(supervision_state_path(run_dir), state)


def clear_owner(run_dir: Path) -> None:
    """所有を解いたことを記録する。

    キーごと消さずに `status` を落とすだけにする。「誰が最後に走らせたか」は
    失敗の追跡に要る情報で、実行が終わったこととは別の事実だから。
    """
    owner = read_owner(run_dir)
    if owner is None:
        return
    owner = {**owner, "status": "finished", "finished_at": _utc_now()}
    write_owner(run_dir, owner)


def owner_summary(run_dir: Path) -> dict | None:
    """所有者の要点（誰が・いつから・生きているか）を返す。無ければ None。

    拒否の封筒へ入れるためのもの。「別の何かが使っています」だけでは、呼び出し側は
    待てばよいのか手で直すのか判断できない。所有者の種別（単体実行か pipeline か）
    まで返す。
    """
    owner = read_owner(run_dir)
    if owner is None:
        return None
    return {"kind": owner.get("kind"), "pid": owner.get("pid"),
            "since": owner.get("since"), "status": owner.get("status"),
            "active": owner_is_active(run_dir)}


def owner_is_active(run_dir: Path) -> bool:
    """所有者のプロセスが**今も**生きているかを返す。

    判定は process identity（pid ＋ 生成時刻）で行う。pid だけを見ると、
    再利用された pid を掴んで無関係なプロセスを「まだ走っている」と誤読する。

    identity を確認できない環境（Windows 以外）では、記録上 running のままの
    ものを active として扱う——分からないものを「終わっている」と読み替えて
    実行中の生成物を消させないため。

    **記録された status を identity より先に見る。** `clear_owner` は所有を
    解いたことを `status="finished"` で表すだけで、identity は「誰が最後に
    走らせたか」として残す。同期実行（`console_run` の detach=False）の所有者は
    MCP サーバ自身なので、identity だけで判定すると `same_process` は
    サーバが生きている限り真を返し続ける——実行はとっくに終わっているのに
    `console_cleanup` が永久に JOB_BUSY になる。
    """
    owner = read_owner(run_dir)
    if owner is None:
        return False
    if owner.get("status") != "running":
        return False
    identity = owner.get("identity")
    if not isinstance(identity, dict):
        return owner.get("status") == "running"
    try:
        return same_process(identity)
    except DomainError:
        return owner.get("status") == "running"


# ---------- 実行 ----------

def annotate_job(job_path: Path) -> list[str]:
    """収集済みの成果物から、証跡だけでは読めない所見をジョブへ書き足す。

    `supervise` が確定させるのは機械が検証できる契約（終了証跡・完了ゲート）で、
    ここで足すのは人が読む所見（mzTab から採る software 版・極性の裏取り・
    ファイル名と宣言の食い違い・未対応 mzTab・サンプル別ファイルの欠落）。
    同期実行と切り離し実行で結果が変わらないよう、両方がこの関数を通る。

    追加した warning を返す。
    """
    job = load_job(job_path)
    record_mztab_provenance(job, job.primary_mztab_files)
    added = [
        *meta_conflict_warnings(job.primary_mztab_files),
        *polarity_crosscheck_warnings(job.primary_mztab_files),
        *unsupported_mztab_warnings(job.artifacts),
        *missing_per_sample_output_warnings(job.artifacts),
    ]
    job.warnings.extend(added)
    save_job(job, job_path)
    return added


def run_job(job_path: Path, *, console_command: list[str] | None = None) -> dict:
    """ジョブを排他的に実行し、終了証跡を返す。

    ロックを保持したまま監視・収集・保存・所見付与までを終える。ロックを
    取れなければ `DomainError("LOCK_TIMEOUT")`——別のプロセスが同じジョブを
    走らせている状態で二重に起動すると、同じ出力先へ 2 つの Console が書く。
    ロックを取れても、待っている間に別プロセスが走り切っていれば
    `DomainError("JOB_ALREADY_FINISHED")`（先行の結果を上書きしない）。

    `console_command` は試験用の注入口（`supervise` の `command` と同じ）。
    """
    job_path = Path(job_path)
    run_dir = Path(load_job(job_path).run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    with file_lock(lock_path(run_dir), timeout=LOCK_TIMEOUT_S):
        # ロックを取ってから、もう一度ジョブの状態を見る。ロック待ちの間に
        # 先行プロセスが同じジョブを走らせ切っていることがあり（同時に届いた
        # 2 つの console_run は、どちらも status=planned を見てからここへ来る）、
        # そのまま進むと MS-DIAL が同じ出力先へ二度走る。終端状態に達した
        # ジョブは実行し直さない——やり直すなら計画からやり直す。
        current = load_job(job_path)
        if current.status in _FINISHED_STATUSES:
            raise DomainError(
                "JOB_ALREADY_FINISHED",
                f"このジョブは既に実行を終えています（status={current.status}）。"
                "再実行が必要なら console_plan からやり直してください。",
                {"job_path": str(job_path), "status": current.status})

        identity = _current_identity()
        write_owner(run_dir, {
            "kind": OWNER_KIND_CONSOLE,
            "pid": identity["pid"],
            "identity": identity,
            "job_path": str(job_path),
            "status": "running",
            "since": _utc_now(),
        })
        # running はロックを取ってから書く。取る前に書くと、ロックを取れなかった
        # 側が「誰も走っていないのに running」という状態を残す。
        update_status(job_path, "running")
        try:
            receipt = supervise(job_path, command=console_command)
            if receipt["job_save"]["status"] == "succeeded":
                # 保存できなかったジョブへ所見だけ足すと、状態と所見の出所が
                # ずれる。証跡の復旧情報のほうが先。
                try:
                    annotate_job(job_path)
                except Exception as exc:  # 所見の失敗で証跡を失わせない
                    receipt["annotation"] = {"status": "failed", "error": repr(exc)}
        finally:
            clear_owner(run_dir)
    return receipt


def _current_identity() -> dict:
    """このプロセスの identity。取得できない環境では pid だけを記録する。"""
    import os

    from metabolomix.core.process_control import process_identity
    try:
        identity = process_identity(os.getpid())
    except DomainError:
        identity = None
    return identity or {"pid": os.getpid(), "creation_time": None}


# ---------- 切り離し起動 ----------

def _repo_root() -> Path:
    """`-m metabolomix.console.worker` を解決できる作業ディレクトリ。

    `-m` は cwd を sys.path の先頭に置く。ここを間違えるとワーカーは
    ModuleNotFoundError で即死し、その事実は worker.log にしか残らない。
    """
    return Path(__file__).resolve().parents[2]


def launch_console_worker(job_path: Path, *,
                          console_args: list[str] | None = None) -> dict:
    """監視ワーカーを親から切り離して起動し、pid と identity を返す（待たない）。

    返る pid は**ワーカーの pid** で、Console のものではない。Console の pid は
    起動後に `worker.json` の `launch` と終了証跡へ載る。ここを取り違えると、
    「まだ走っているか」の判定が別のプロセスを見に行く。

    `console_args` は試験用の注入口（`run_job(console_command=...)` に渡る）。
    """
    job_path = Path(job_path)
    run_dir = Path(load_job(job_path).run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    command = [sys.executable, "-m", "metabolomix.console.worker", "--job", str(job_path)]
    for arg in console_args or []:
        # `--console-arg=<値>` の形で渡す。値は `-i` のようにハイフンで始まるので、
        # 空白区切りだと argparse が次のオプション名と読み違える。
        command.append(f"--console-arg={arg}")

    return launch_detached(command, cwd=_repo_root(),
                           log_path=worker_log_path(run_dir))


# ---------- CLI ----------

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="metabolomix.console.worker",
        description="MS-DIAL Console を監視し、終了証跡と成果物を残す。")
    parser.add_argument("--job", required=True, help="analysis-job.json のパス")
    parser.add_argument("--console-arg", dest="console_arg", action="append",
                        default=None,
                        help="試験用: 実 Console の代わりに起動するコマンドの引数")
    args = parser.parse_args(argv)

    receipt = run_job(Path(args.job), console_command=args.console_arg)
    # 終了コードは「監視できたか」であって Console の成否ではない。Console の
    # 終了コードは証跡（execution-result.json）にある。
    print(f"execution_id={receipt['execution_id']} "
          f"termination={receipt['termination']} exit_code={receipt['exit_code']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
