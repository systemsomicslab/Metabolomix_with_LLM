"""切り離し起動の中継（`process_control.launch_detached`だけが起こす）。

MCPサーバが直接ワーカーを起こすと、ワーカーの親はサーバになる。デスクトップアプリが
呼び出しを打ち切ってサーバをプロセスの木ごと消すと（`taskkill /T`と同じく親子の記録を
たどる）、切り離しのフラグを付けていてもワーカーが巻き込まれる。そこでサーバはこの中継を
起こし、中継がワーカーを起こしてそのpidとidentityを書き、すぐ終わる。ワーカーの親（中継）は
もう居ないので、サーバの木をたどってもワーカーに届かない。

    python -m metabolomix.core.detach_relay <要求のJSON>

要求は`{"command", "cwd", "log_path", "launcher", "result_path"}`（`command`は
`resolve_python_launch`を通した後のもの。`launcher`はvenvの印`__PYVENV_LAUNCHER__`の値か
null）。結果は`result_path`へJSONで書く: 起こせたら`{"pid", "identity", "breakaway"}`
（pidとidentityはワーカーのもの）、`DomainError`なら`{"error": {"code", "message",
"details"}}`。書いたらすぐ終わる。それ以外の失敗はstderr（中継のログ）へ出して0以外で終わる。

依存はstdlibと`metabolomix.core`のleaf（`atomic_io`・`process_control`）だけ。
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from metabolomix.core import process_control
from metabolomix.core.atomic_io import DomainError, atomic_write_json

__all__ = ["main"]


def _worker_environment(launcher: str | None) -> dict | None:
    """ワーカーへ渡す環境。Noneならこの中継の環境をそのまま受け継がせる。

    中継の環境は呼び出し側の環境そのもの（解釈系は起動時に`__PYVENV_LAUNCHER__`を
    自分の環境から消す）。呼び出し側の`resolve_python_launch`がvenvの印を付けていた
    （`launcher`がある）ときだけ、同じ値を足し直す。
    """
    name = process_control._PYVENV_LAUNCHER_ENV
    if launcher:
        env = dict(os.environ)
        env[name] = launcher
        return env
    if name in os.environ:
        env = dict(os.environ)
        env.pop(name)
        return env
    return None


def _json_safe(details: dict) -> dict:
    return json.loads(json.dumps(details or {}, default=str))


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print("usage: python -m metabolomix.core.detach_relay <要求のJSON>", file=sys.stderr)
        return 2
    # 中継自身のstdout/stderr（中継のログ）をワーカーへ受け継がせない。
    for fd in (0, 1, 2):
        try:
            os.set_inheritable(fd, False)
        except OSError:
            pass
    request = json.loads(argv[0])
    result_path = Path(request["result_path"])
    try:
        launched = process_control._spawn_detached(
            [str(part) for part in request["command"]], Path(request["cwd"]),
            Path(request["log_path"]), _worker_environment(request.get("launcher")))
    except DomainError as exc:
        atomic_write_json(result_path, {"error": {
            "code": exc.code, "message": exc.message, "details": _json_safe(exc.details)}})
        return 1
    atomic_write_json(result_path, launched)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
