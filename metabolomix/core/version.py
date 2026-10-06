"""稼働中のサーバ版数と、origin/main からの遅れを 1 行で出すためのモジュール。

**なぜ要るか**: MCP サーバはクライアント（Claude Desktop 等）が起動したまま
生き続けるので、`git pull` でコードを更新しても**次の会話でも古いプロセスが
動き続ける**。実際に、削除済みの caveat 文面が返ってきて初めて気づいた
（docs/HISTRY.md 2026-09-04(4)）。まず版数を状態ツールの戻り値に刻み、
ズレていたらサーバを再起動する、というのが元の対処。

配布先（DEPLOY.md のメンバー各自のクローン）では、そもそも更新が出たことに
気づけない。そこで起動時に 1 回だけ `git fetch` し、`origin/main` より遅れて
いれば通知する。**自動 pull はしない** —— 走っているプロセスの裏でコードが
差し替わると、import 済みモジュールだけが古いまま残る混在状態になり、
上の事故より厄介な壊れ方をする。更新の適用はユーザが明示的に行う。

判定できないときは**すべて黙る**（オフライン・認証失敗・git 不在・main 以外の
ブランチ）。通知機能が解析の邪魔をしてはいけないので、不明は通知なしに倒す。
"""
from __future__ import annotations

import subprocess
import sys
import threading
from functools import lru_cache
from pathlib import Path

VERSION = "0.1.0"

#: 配布の正準ブランチ。ユーザのクローンはここに乗っている前提。
UPSTREAM_BRANCH = "main"
#: ローカルの git 操作（ネットワークを伴わない）の上限。
GIT_TIMEOUT_SEC = 10.0
#: fetch の上限。daemon スレッドで走るので、切れても誰も待たされない。
FETCH_TIMEOUT_SEC = 15.0
#: pip install の上限。依存が増えた更新ではホイールの取得に分単位かかりうる。
PIP_TIMEOUT_SEC = 600.0
#: 依存表。pull でこれが変わったときだけ pip を回す。
REQUIREMENTS_FILE = "requirements.txt"
#: pip の出力を戻り値へ載せる上限。失敗の理由は末尾に出るので後ろから取る。
PIP_LOG_TAIL_CHARS = 500

_UNSET = object()
#: fetch が**成功した**ことを示す。未完了と失敗はどちらも「黙る」なので同一視してよい。
_fetch_ok = threading.Event()
_update_cache: object = _UNSET
_check_started = False
_check_lock = threading.Lock()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git(*args: str, timeout: float = GIT_TIMEOUT_SEC) -> str | None:
    """git の標準出力（strip 済み）を返す。失敗・git 不在なら None。

    呼び出し側は None を「判定できなかった」として黙る側に倒すこと。
    """
    try:
        out = subprocess.run(
            ["git", "-C", str(_repo_root()), *args],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            stdin=subprocess.DEVNULL, timeout=timeout, text=True, errors="replace")
    except (OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    return (out.stdout or "").strip()


@lru_cache(maxsize=1)
def server_version() -> str:
    """`<version>+<git の短縮 SHA>` を返す。git が無ければ版数だけ。

    プロセス生存中は変わらないのでキャッシュする（起動時のコードが何かを
    示す値なので、途中で作業ツリーが変わっても更新しないのが正しい）。
    """
    sha = _git("rev-parse", "--short", "HEAD")
    return f"{VERSION}+{sha}" if sha else VERSION


def _run_update_fetch() -> None:
    """`origin/<UPSTREAM_BRANCH>` の ref を 1 回だけ更新する。

    成功したときだけフラグを立てる。失敗を握り潰すのは意図的で、
    オフラインのメンバーに通知の失敗を見せる意味がないため。
    """
    if _git("fetch", "--quiet", "origin", UPSTREAM_BRANCH,
            timeout=FETCH_TIMEOUT_SEC) is not None:
        _fetch_ok.set()


def start_update_check() -> None:
    """起動時に 1 回だけ、裏で fetch を投げる（プロセスにつき 1 回）。

    ツール呼び出しの中で同期実行すると、最初の 1 呼び出しがネットワーク待ちで
    止まる。daemon スレッドにしてあるので、終わっていなければ通知しないだけで
    サーバの終了を妨げない。
    """
    global _check_started
    with _check_lock:
        if _check_started:
            return
        _check_started = True
    threading.Thread(target=_run_update_fetch,
                     name="metabolomix-update-check", daemon=True).start()


def _compute_update_status() -> dict | None:
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    if branch != UPSTREAM_BRANCH:
        # 開発者は feature ブランチに、巻き戻し検証は detached HEAD にいる。
        # どちらも「配布物が遅れている」話ではないので鳴らさない。
        return None
    count = _git("rev-list", "--count", f"HEAD..origin/{UPSTREAM_BRANCH}")
    if count is None:
        return None
    try:
        behind = int(count)
    except ValueError:
        return None
    if behind <= 0:
        return None
    return {
        "behind": behind,
        "message": (f"origin/{UPSTREAM_BRANCH} より {behind} コミット遅れています。"
                    "`server_update` ツールを呼ぶか、手元で git pull してください。"
                    "どちらの場合も反映には MCP クライアントの再起動が要ります"
                    "（pull だけでは起動中の古いプロセスが動き続けます）。"),
    }


def update_status() -> dict | None:
    """遅れていれば `{"behind": N, "message": ...}`、それ以外は None。

    fetch が終わるまでは None を返すが**キャッシュしない**。起動直後の
    1 回目で黙ったことを覚えてしまうと、その後 fetch が届いても永久に
    通知できなくなるため。
    """
    global _update_cache
    if not _fetch_ok.is_set():
        return None
    if _update_cache is _UNSET:
        _update_cache = _compute_update_status()
    return _update_cache  # type: ignore[return-value]


# --- 更新の適用 -------------------------------------------------------------
#
# 通知（update_status）が「遅れている」と言った後、利用者がその場で適用できる
# ようにするための実体。`metabolomix/tools/maintenance.py` の `server_update` が
# 薄く包んで MCP に出す。
#
# **プロセスは触らない**。MCP サーバは Claude Desktop の子プロセスなので、
# ここで親を終了させるとツールの応答を返す前に自分が道連れで死ぬ。再起動は
# 利用者が手で行い、ここは「再起動してください」と言うところまでを担う。
#
# 拒否を厚くしてあるのは、失敗の仕方が両方とも静かだから: pull が手元の
# 未コミット変更を巻き込むと復元手段がなく、pip の失敗を見逃すとコードだけ
# 新しく依存が古い状態で次回起動時に ImportError になる。


def _pip_install() -> tuple[bool, str]:
    """`requirements.txt` を現在の Python へ入れ直す。`(成功したか, 出力)`。

    `sys.executable` を使うのは、PATH 上の `pip` が別の Python を指しうるため
    （このプロジェクトは `C:/Python314/python.exe` 固定で、`.venv` を持たない）。
    """
    requirements = _repo_root() / REQUIREMENTS_FILE
    if not requirements.is_file():
        return True, ""
    try:
        out = subprocess.run(
            [sys.executable, "-m", "pip", "install", "-r", str(requirements)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL, timeout=PIP_TIMEOUT_SEC,
            text=True, errors="replace")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{type(exc).__name__}: {exc}"
    return out.returncode == 0, (out.stdout or "").strip()


def _refused(reason: str, message: str) -> dict:
    return {"status": "refused", "reason": reason, "message": message}


def apply_update() -> dict:
    """クローンを `origin/<UPSTREAM_BRANCH>` へ早送りし、必要なら依存も入れ直す。

    戻り値の `status` は `up_to_date` / `updated` / `updated_with_warning` /
    `refused` のいずれか。`refused` のときは `reason` に機械可読な理由が入る。
    例外は投げない——更新は解析の本筋ではないので、呼び出し側に分岐を強いない。
    """
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    if branch is None:
        return _refused(
            "git_unavailable",
            "git が使えないか、このサーバが git クローンではありません。"
            "配布物を展開しただけの場合は、クローンし直してください。")
    if branch != UPSTREAM_BRANCH:
        return _refused(
            "not_on_main",
            f"現在のブランチは {branch} です。配布物の更新は "
            f"{UPSTREAM_BRANCH} 上でのみ行います（開発中の作業を早送りしません）。")

    if _git("fetch", "--quiet", "origin", UPSTREAM_BRANCH,
            timeout=FETCH_TIMEOUT_SEC) is None:
        return _refused(
            "fetch_failed",
            "origin から取得できませんでした（オフライン・認証切れなど）。"
            "古い情報のまま「最新です」とは言えないので、何もしていません。")

    count = _git("rev-list", "--count", f"HEAD..origin/{UPSTREAM_BRANCH}")
    try:
        behind = int(count)
    except (TypeError, ValueError):
        return _refused(
            "distance_unknown",
            f"origin/{UPSTREAM_BRANCH} との差分を数えられませんでした。何もしていません。")

    running = server_version()
    if behind <= 0:
        return {"status": "up_to_date", "head": _git("rev-parse", "--short", "HEAD"),
                "message": f"すでに最新です（{running}）。更新は不要です。"}

    # 以降は実際に作業ツリーへ書き込む。遅れていないときに手元の汚れを
    # 咎めないよう、この検査は「何もしない」分岐を通り過ぎてから行う。
    dirty = _git("status", "--porcelain")
    if dirty is None:
        return _refused("git_unavailable", "作業ツリーの状態を確認できませんでした。")
    if dirty:
        return _refused(
            "dirty_worktree",
            "未コミットの変更があるため更新しませんでした。"
            "巻き込むと復元できないためです。`git status` で確認し、"
            "退避（コミット）してから呼び直してください。")

    changed = _git("diff", "--name-only", f"HEAD..origin/{UPSTREAM_BRANCH}") or ""
    requirements_changed = any(
        line.strip() == REQUIREMENTS_FILE for line in changed.splitlines())

    if _git("merge", "--ff-only", f"origin/{UPSTREAM_BRANCH}") is None:
        return _refused(
            "not_fast_forward",
            "早送りできませんでした（このクローンに独自のコミットがあります）。"
            "merge も rebase も自動では行いません。手元で解決してください。")

    head = _git("rev-parse", "--short", "HEAD")
    restart = ("反映するには Claude Desktop（MCP クライアント）を再起動してください。"
               "起動中のプロセスは古いコードのまま動き続けます。")
    if not requirements_changed:
        return {
            "status": "updated", "behind": behind, "head": head,
            "dependencies": "unchanged",
            "message": f"{behind} コミット分を更新しました（{running} → {head}）。{restart}",
        }

    installed, log = _pip_install()
    if installed:
        return {
            "status": "updated", "behind": behind, "head": head,
            "dependencies": "installed",
            "message": (f"{behind} コミット分を更新し、依存も入れ直しました"
                        f"（{running} → {head}）。{restart}"),
        }
    return {
        "status": "updated_with_warning", "behind": behind, "head": head,
        "dependencies": "failed",
        "message": (f"コードは {behind} コミット分更新しましたが、依存の更新に失敗しました"
                    f"（{running} → {head}）。このまま再起動すると起動しない可能性があります。"
                    f"手元で `pip install -r {REQUIREMENTS_FILE}` を実行してください。"
                    f" pip の出力（末尾）: {log[-PIP_LOG_TAIL_CHARS:]}"),
    }
