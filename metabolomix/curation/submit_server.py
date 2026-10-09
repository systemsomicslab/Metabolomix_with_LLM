"""キュレーションのビューアが直接送信する受け口（127.0.0.1 の HTTP）。spec 2026-10-09。

プロセスに 1 つを共有し、curation_review / curation_suggest のたびにトークンを登録する。
登録は `{kind, review_id, arf2_path, last_seen}`。書き込み先は登録側の値だけを使い、
ページからは flags しか受け取らない（`/v1/ping` `/v1/submit` `/v1/finish`）。Host は
127.0.0.1:<port> のみ許し、本文は MAX_BODY_BYTES まで。最後の通信から `IDLE_TIMEOUT_MIN` 分で登録を外し、
0 件になれば待受けを止める。daemon スレッドなので MCP サーバが終われば一緒に消える。

**stdout に何も書かない**（stdio の MCP ではプロトコルの通り道）。`log_message` も止める。
deps: curation.submission、標準ライブラリ。session を import しない（HTTP スレッドから呼ぶ）。
"""
from __future__ import annotations

import hmac
import json
import secrets
import socketserver
import threading
import time
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from metabolomix.curation import flags as flag_log
from metabolomix.curation import submission

IDLE_TIMEOUT_MIN = 30
SWEEP_INTERVAL_S = 60
MAX_BODY_BYTES = 8 * 1024 * 1024        # MAX_SPOTS（3000 件）にメモが付いても収まる
KINDS = ("review", "suggest")

_clock = time.monotonic                 # テストが差し替える
_lock = threading.Lock()
_entries: dict[str, dict] = {}          # token -> {kind, review_id, arf2_path, last_seen}
_server: "_Server | None" = None
_sweeper_stop: threading.Event | None = None


def _error_body(message: str, **details) -> dict:
    return {"status": "error", "message": message, **details}


class _Server(ThreadingHTTPServer):
    allow_reuse_address = False         # Windows の SO_REUSEADDR は同じポートの横取りを許す

    def server_bind(self):
        # HTTPServer.server_bind は socket.getfqdn() を呼び、Windows では逆引きで数秒止まることがある。
        socketserver.TCPServer.server_bind(self)
        self.server_name, self.server_port = self.server_address[:2]


class _Handler(BaseHTTPRequestHandler):
    server_version = "ms-data-parser-curation"
    sys_version = ""
    timeout = 30                    # 止まったクライアントがハンドラスレッドを握り続けないように

    def log_message(self, format, *args):       # 既定のアクセスログ（stderr）を出さない
        pass

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Private-Network", "true")

    def _reply(self, code: int, body: dict):
        data = json.dumps(body, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self._cors()
        self.end_headers()
        self.wfile.write(data)

    def _host_ok(self) -> bool:      # DNS rebinding 対策: 127.0.0.1:<port> 以外の Host は拒む
        return self.headers.get("Host") == f"127.0.0.1:{self.server.server_address[1]}"

    def do_OPTIONS(self):
        if not self._host_ok():
            return self._reply(403, _error_body("Host が不正です。"))
        self.send_response(204)
        self._cors()
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_GET(self):               # GET を返さないので、他のサイトから中身を覗けない
        self._reply(405, _error_body("POST だけを受け付けます。"))

    def _read_json(self):
        """本文を dict として読む。読めなければ (None, (code, payload))。"""
        try:
            length = int(self.headers.get("Content-Length") or "")
        except ValueError:
            length = -1
        if length < 0:      # 負の値で rfile.read すると EOF までハンドラスレッドが止まる
            return None, (400, _error_body("Content-Length がないか不正です。"))
        if length > MAX_BODY_BYTES:
            # 読まずに閉じると Windows では RST になり、ページに 413 が届かない。上限の 4 倍までは
            # 読み捨ててから答え、それを超える申告は読まずに閉じる。
            self.close_connection = True
            if length <= 4 * MAX_BODY_BYTES:
                remaining = length
                while remaining > 0:
                    chunk = self.rfile.read(min(remaining, 65536))
                    if not chunk:
                        break
                    remaining -= len(chunk)
            return None, (413, _error_body(f"本文が上限 {MAX_BODY_BYTES} バイトを超えています。"))
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            return None, (400, _error_body("本文が JSON ではありません。"))
        if not isinstance(body, dict):
            return None, (400, _error_body("本文は JSON オブジェクトにしてください。"))
        return body, None

    def do_POST(self):
        if not self._host_ok():
            return self._reply(403, _error_body("Host が不正です。"))
        route = _ROUTES.get(self.path)
        if route is None:
            return self._reply(404, _error_body(f"不明なパスです: {self.path}"))
        body, problem = self._read_json()
        if problem:
            return self._reply(*problem)
        token = body.get("token")
        entry = _touch(token) if isinstance(token, str) else None
        if entry is None:
            return self._reply(401, _error_body("トークンが不明か期限切れです。送信用テキストをコピーしてチャットに貼ってください。"))
        try:
            code, payload = route(self, token, entry, body)
        except Exception as exc:     # noqa: BLE001 — 接続を落とさず理由をページへ返す
            code, payload = 500, _error_body(f"{type(exc).__name__}: {exc}")
        self._reply(code, payload)
        if route is _finish:         # 応答を返した後で外す（最後の 1 件なら別スレッドで止まる）
            unregister(token, wait=False)


def _ping(handler, token, entry, body):
    return 200, {"status": "ok", "kind": entry["kind"], "review_id": entry["review_id"],
                 "expires_at": _expires_at(entry)}


_SUBMISSION_STATUS = {"invalid": 400, "alignment_changed": 409, "flag_file": 500}


def _submit(handler, token, entry, body):
    if body.get("review_id") != entry["review_id"]:
        return 403, _error_body("review_id がこのビューアの登録と違います（別の HTML の取り違え）。")
    try:
        saved = submission.load_saved(flag_log.curation_dir(entry["arf2_path"]), entry["review_id"])
    except FileNotFoundError:
        return 410, _error_body(f"review_id={entry['review_id']} の保存済みレビューが見つかりません。"
                                "curation_review（候補付けなら curation_suggest）をやり直してください。")
    try:
        return 200, submission.submit_flags(saved, body.get("flags"), review_id=entry["review_id"],
                                            source="user")
    except submission.SubmissionError as exc:
        return _SUBMISSION_STATUS[exc.kind], _error_body(str(exc), **exc.details)


def _finish(handler, token, entry, body):
    return 200, {"status": "ok"}


_ROUTES = {"/v1/ping": _ping, "/v1/submit": _submit, "/v1/finish": _finish}


def _expired(entry: dict, now: float) -> bool:
    return now - entry["last_seen"] >= IDLE_TIMEOUT_MIN * 60


def _expires_at(entry: dict) -> str:
    remaining = IDLE_TIMEOUT_MIN * 60 - (_clock() - entry["last_seen"])
    return (datetime.now(timezone.utc) + timedelta(seconds=remaining)).isoformat(timespec="seconds")


def _touch(token: str) -> dict | None:
    """トークンを定数時間で照合し、生きていれば期限を延ばして登録の写しを返す。"""
    now = _clock()
    given = token.encode("utf-8")
    with _lock:
        for key, entry in _entries.items():
            if hmac.compare_digest(key.encode("utf-8"), given):
                if _expired(entry, now):
                    return None
                entry["last_seen"] = now
                return dict(entry)
    return None


def _ensure_running() -> None:
    global _server, _sweeper_stop
    if _server is not None:
        return
    server = _Server(("127.0.0.1", 0), _Handler)
    threading.Thread(target=server.serve_forever, name="curation-submit", daemon=True).start()
    stop = threading.Event()
    threading.Thread(target=_sweep_loop, args=(stop,), name="curation-submit-sweep", daemon=True).start()
    _server, _sweeper_stop = server, stop


def _sweep_loop(stop: threading.Event) -> None:
    while not stop.wait(SWEEP_INTERVAL_S):
        sweep()


def _detach_if_idle():
    """（_lock の中で呼ぶ）登録が 0 件なら待受けを切り離して返す。止めるのは呼び出し側が鍵の外で。"""
    global _server, _sweeper_stop
    if _entries or _server is None:
        return None
    detached = (_server, _sweeper_stop)
    _server = _sweeper_stop = None
    return detached


def _shutdown(detached, *, wait: bool) -> None:
    server, stop = detached
    stop.set()

    def run():
        server.shutdown()
        server.server_close()
    if wait:
        run()
    else:       # 要求を処理しているスレッドからは、応答を返した後に別スレッドで止める
        threading.Thread(target=run, name="curation-submit-stop", daemon=True).start()


def register(kind: str, review_id: str, arf2_path) -> dict:
    if kind not in KINDS:
        raise ValueError(f"kind は {KINDS} のいずれか: {kind!r}")
    with _lock:
        _ensure_running()
        token = secrets.token_urlsafe(32)
        _entries[token] = {"kind": kind, "review_id": review_id, "arf2_path": str(arf2_path),
                           "last_seen": _clock()}
        port = _server.server_address[1]
    return {"port": port, "token": token, "idle_timeout_min": IDLE_TIMEOUT_MIN}


def unregister(token: str, *, wait: bool = True) -> None:
    with _lock:
        _entries.pop(token, None)
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=wait)


def sweep() -> None:
    now = _clock()
    with _lock:
        for token in [t for t, e in _entries.items() if _expired(e, now)]:
            del _entries[token]
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=True)


def is_running() -> bool:
    with _lock:
        return _server is not None


def current_port() -> int | None:
    with _lock:
        return None if _server is None else _server.server_address[1]


def reset() -> None:
    """テスト用: 登録を全部外して待受けを止める。"""
    with _lock:
        _entries.clear()
        detached = _detach_if_idle()
    if detached:
        _shutdown(detached, wait=True)
