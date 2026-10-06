"""数値・ドメイン層で共有する原子的JSON保存とhash。

MCPに依存しない依存グラフのleaf（stdlibのみ）。metabolomix.core.mcp_core /
metabolomix.<形式>.tools / metabolomix.tools.* をimportしてはいけない
（`metabolomix/core/mcp_core.py` はこのモジュールを含む「leaf」からさらに
上位を組み立てる側なので、循環を避けるためここからは何も引かない）。

`DomainError`は本plan（生データフォルダ起点pipeline）の数値・ドメイン層が
共通で使う唯一の例外。console終了証跡・pipeline状態・結果fingerprintなど、
後続タスクの多くがここからimportする。
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from pathlib import Path

#: 共有違反を吸収する有界リトライの上限（秒）と待ち時間。
#:
#: **Windowsでは、誰かがそのファイルを読むために開いているだけで`os.replace`が
#: `PermissionError`（WinError 5）になる。** 読み手側に`FILE_SHARE_DELETE`を
#: 付けても`MoveFileExW(REPLACE_EXISTING)`は通らないので、読み手側の書き方では
#: 直せない（Task19が実測）。つまり`pipeline-run.json`をポーリングして眺めるだけで
#: **動いているworkerの保存を失敗させられる**。監視のために解析を落とすのは
#: 本末転倒なので、置換する側が短時間だけ譲って待つ。
#:
#: 逆向き（置換の最中に読もうとして開けない）も同じ性質で起こるため、
#: `read_text_stable`が読み手側にも同じ有界リトライを与える。
#:
#: 「有界」であることが要点——待つのはOSレベルの一瞬の共有違反だけで、
#: 権限不足のような恒久的な失敗はこの窓を超えた時点で元の例外のまま送出する
#: （握り潰して静かに書けていないことにしない）。実測では最悪10回程度の
#: 再試行（≒0.3秒）で通る。
_SHARING_RETRY_DEADLINE_S = 2.0
_SHARING_RETRY_FIRST_S = 0.001
_SHARING_RETRY_MAX_S = 0.05


class DomainError(Exception):
    """数値・ドメイン層の共通例外。

    `code`は機械可読な種別（例 "EXECUTION_RECORD_INVALID"）、`message`は
    人間可読な日本語説明。`str(exc)`は`f"{code}: {message}"`を返す。
    """

    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message
        self.details = details if details is not None else {}

    def __str__(self) -> str:
        return f"{self.code}: {self.message}"


def canonical_hash(value: object) -> str:
    """値を正規化JSON（key昇順・改行なし）にしてSHA-256を返す。

    呼び出し側の責務: 日時・UUIDなど実行のたびに変わる値をvalueへ混ぜない
    （fingerprintに再現性がなくなるため）。
    """
    raw = json.dumps(value, ensure_ascii=False, sort_keys=True,
                     separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _retry_on_sharing_violation(operation):
    """共有違反（`PermissionError`）だけを有界に再試行する。

    `operation`は1回分の試行。成功したらその戻り値を返す。期限を過ぎたら
    **元の例外をそのまま**送出する（別のcodeに包み直さない——呼び出し側の
    既存のエラー処理を変えないため）。
    """
    deadline = time.monotonic() + _SHARING_RETRY_DEADLINE_S
    delay = _SHARING_RETRY_FIRST_S
    while True:
        try:
            return operation()
        except PermissionError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(delay)
            delay = min(delay * 2, _SHARING_RETRY_MAX_S)


def read_text_stable(path: Path, *, encoding: str = "utf-8") -> str:
    """`atomic_write_json`で置換されうるファイルを読む（共有違反だけ有界に再試行）。

    置換の瞬間に開こうとすると`PermissionError`になりうる。ここで吸収しないと、
    読取専用のはずの監視（`pipeline_status`）が書き込み中にランダムで失敗する。
    `FileNotFoundError`はそのまま送出する——「まだ無い」は競合ではなく事実。
    """
    path = Path(path)
    return _retry_on_sharing_violation(lambda: path.read_text(encoding=encoding))


def atomic_write_bytes(path: Path, data: bytes) -> None:
    """バイト列を原子的に保存する（手順と失敗時の扱いは `atomic_write_json` と同じ）。"""
    path = Path(path)
    parent = path.parent
    parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=parent, prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as tmp_file:
            tmp_path = Path(tmp_file.name)
            tmp_file.write(data)
            tmp_file.flush()
            os.fsync(tmp_file.fileno())
        _retry_on_sharing_violation(lambda: os.replace(tmp_path, path))
    except Exception:
        if tmp_path is not None and tmp_path.exists():
            os.unlink(tmp_path)
        raise


def atomic_write_json(path: Path, data: dict) -> None:
    """JSONを原子的に保存する。

    同じ親ディレクトリへNamedTemporaryFileで書き、UTF-8・allow_nan=False で
    直列化し、flush・os.fsync してから os.replace で置換する。この順序を
    崩すと、置換直前にプロセスが落ちた際に中途半端な内容を確定状態として
    読ませてしまう。

    置換は`_retry_on_sharing_violation`越しに行う。Windowsでは**誰かがそのファイルを
    読むために開いているだけで`os.replace`が失敗する**ため、これが無いと
    `pipeline_status`のポーリングが動いているworkerの保存を落とせてしまう
    （Task19実測。読み手側の共有モードでは直せない）。

    失敗時は未確定の一時ファイルだけを片付け、既存の path には一切触れない
    （path.unlink はしない）。
    """
    path = Path(path)
    parent = path.parent
    parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as tmp_file:
            tmp_path = Path(tmp_file.name)
            json.dump(data, tmp_file, ensure_ascii=False,
                     separators=(",", ":"), allow_nan=False)
            tmp_file.flush()
            os.fsync(tmp_file.fileno())
        _retry_on_sharing_violation(lambda: os.replace(tmp_path, path))
    except Exception:
        if tmp_path is not None and tmp_path.exists():
            os.unlink(tmp_path)
        raise
