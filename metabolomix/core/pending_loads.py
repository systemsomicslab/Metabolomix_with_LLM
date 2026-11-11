"""時間のかかる `load_dataset` の受付（受付番号 `ld_…`）と、その結果の控え。

Claude のデスクトップアプリの中継は MCP の呼び出しを約 80 秒で打ち切る。読み込みがワーカー
スレッドで回っていても（`handle_dispatch.ALWAYS_THREADED`）、応答が遅れれば中継が諦める。
そこで `load_dataset` は読み込みをこの受付のスレッドで始め、期限（`wait_seconds`）の中で
終われば戻り値をそのまま返し、終わらなければ受付番号を返す。クライアントは
`dataset_load_status(load_id)` で結果を引く。

- 受付が持つのは読み込みの戻り値（JSON の文字列。成功の戻り値か誤りの封筒）だけ。ハンドルの
  登録と直近の切り替えは、読み込み自身（`handle_tools._load_now` の `adopt_handle`）が終わる
  ときに 1 回だけする。結果を引いても状態は変わらない。
- 読み込みがどう終わっても受付は必ず終わる。`DomainError` は `load_dataset` と同じ封筒
  （`mztab_error`。`details` が JSON にならなければ `details` を落とす）、それ以外は
  `DATASET_LOAD_FAILED`（例外の文言は置き場所を含みうるので、種類だけを書く。詳細はサーバの
  ログ）。封筒づくり自体が落ちても、最後の封筒で終える。
- `submit` は受付そのもの（`PendingLoad`）を返す。受付を出した呼び出しはそれを直に待つので、
  受付簿から捨てられても結果を失わない。
- 終わった受付は終わってから `TTL_SECONDS`（1 時間）覚え、何度引いても同じ文字列を返す。
  件数が `MAX_ENTRIES`（64）に達したら、終わったものを終わった順に古いものから捨てる。
  読み込み中のものは捨てない。すべてが読み込み中なら `DATASET_LOAD_LIMIT` で受け付けない。

依存: atomic_io・mcp_errors・serialization（どれも leaf）だけ。tools・session は import しない。
"""
from __future__ import annotations

import logging
import threading
import time
import uuid
from collections.abc import Callable

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload

__all__ = ["MAX_ENTRIES", "TTL_SECONDS", "PendingLoad", "PendingLoads", "failure_envelope",
           "registry"]

#: 終わった受付を覚えておく秒数（終わったときから数える）。
TTL_SECONDS = 3600.0
#: 覚えておく受付の件数の上限（読み込み中のものを含む）。
MAX_ENTRIES = 64

_logger = logging.getLogger(__name__)

#: 封筒づくりがどれも落ちたときの最後の封筒（ここは落ちない）。
_LAST_RESORT = json_payload({"error": {
    "code": "DATASET_LOAD_FAILED",
    "message": "読み込みの結果を封筒にできませんでした。詳細はサーバのログにあります。"
               "dataset_list で開いたハンドルを確かめ、load_dataset を呼び直してください。"}})


class PendingLoad:
    """1 つの受付。`payload` と `finished_at` は受付簿のロックの下で一緒に書く。"""

    __slots__ = ("load_id", "started_at", "finished_at", "payload", "done")

    def __init__(self, load_id: str, started_at: float):
        self.load_id = load_id
        self.started_at = started_at
        self.finished_at: float | None = None
        self.payload: str | None = None
        self.done = threading.Event()

    def wait(self, seconds: float) -> str | None:
        """最長 `seconds` 秒待つ。終われば戻り値の JSON 文字列、終わらなければ None。

        受付簿から捨てられた後でも、この受付そのものを待つので結果を失わない。
        """
        timeout = min(max(0.0, float(seconds)), threading.TIMEOUT_MAX)
        if self.done.wait(timeout):
            return self.payload
        return None


class PendingLoads:
    """受付簿。受付・参照・捨てるのはプロセス内の 1 本のロックで守る。

    `clock` は経過と期限を測る時計（既定 `time.monotonic`。テストが差し替える）。
    """

    def __init__(self, *, clock: Callable[[], float] = time.monotonic,
                 ttl_seconds: float = TTL_SECONDS, max_entries: int = MAX_ENTRIES):
        self._clock = clock
        self._ttl = float(ttl_seconds)
        self._max = max(1, int(max_entries))
        self._lock = threading.Lock()
        self._entries: dict[str, PendingLoad] = {}

    def submit(self, fn: Callable[[], str]) -> PendingLoad:
        """`fn` を daemon スレッドで回し、受付（`load_id` は `ld_` + 12 桁）を返す。

        上限なら終わった受付を古いものから捨てて空きを作る。すべてが読み込み中なら
        `DATASET_LOAD_LIMIT`（DomainError。スレッドは起こさない）。
        """
        with self._lock:
            self._prune()
            self._make_room()
            load_id = f"ld_{uuid.uuid4().hex[:12]}"
            entry = PendingLoad(load_id, self._clock())
            self._entries[load_id] = entry
        thread = threading.Thread(target=self._run, args=(entry, fn),
                                  name=f"load_dataset-{load_id}", daemon=True)
        try:
            thread.start()
        except BaseException:
            with self._lock:  # 起こせなかった受付を読み込み中のまま残さない
                self._entries.pop(load_id, None)
            raise
        return entry

    def wait(self, load_id: str, seconds: float) -> str | None:
        """受付番号で引いて最長 `seconds` 秒待つ（`PendingLoad.wait`）。知らない番号は
        `DATASET_LOAD_NOT_FOUND`。0 秒なら待たずに今の状態を見る（既に終わっていれば返す）。"""
        with self._lock:
            self._prune()
            entry = self._lookup(load_id)
        return entry.wait(seconds)

    def status(self, load_id: str) -> tuple[str, str | None]:
        """`("loading", None)` か `("done", payload)`。知らない受付番号は `DATASET_LOAD_NOT_FOUND`。"""
        state, payload, _elapsed = self.snapshot(load_id)
        return state, payload

    def snapshot(self, load_id: str) -> tuple[str, str | None, float]:
        """`status` に受け付けてからの秒（終わったものは読み込みにかかった秒）を足したもの。

        3 つは同じロックの下で読む（状態と経過が食い違わない）。
        """
        with self._lock:
            self._prune()
            entry = self._lookup(load_id)
            if entry.finished_at is None:
                return "loading", None, self._clock() - entry.started_at
            return "done", entry.payload, entry.finished_at - entry.started_at

    def discard(self, load_id: str) -> None:
        """受付を忘れる（無ければ何もしない）。読み込み中のものも忘れるが、読み込みは止めない。"""
        with self._lock:
            self._entries.pop(load_id, None)

    # ---------- 内部 ----------

    def _run(self, entry: PendingLoad, fn: Callable[[], str]) -> None:
        """読み込みのスレッドの本体。どの出口でも戻り値を置いて受付を終える。"""
        payload: object = None
        try:
            try:
                payload = fn()
            except DomainError as exc:  # 本体の try の外で出たもの（本体の中のものは封筒で返る）
                payload = _domain_envelope(exc)
            except BaseException as exc:  # noqa: BLE001 - 受付を読み込み中のまま残さない
                _logger.exception("load_dataset の読み込みが予期しない誤りで止まった（%s）",
                                  entry.load_id)
                payload = failure_envelope(exc)
            if not isinstance(payload, str):  # `wait` の None（終わっていない）と取り違えない
                payload = failure_envelope(
                    TypeError(f"戻り値が文字列ではない: {type(payload).__name__}"))
        except BaseException:  # noqa: BLE001 - 封筒づくりが落ちても受付は終える
            _logger.exception("load_dataset の結果を封筒にできなかった（%s）", entry.load_id)
        finally:
            self._finish(entry, payload if isinstance(payload, str) else _LAST_RESORT)

    def _finish(self, entry: PendingLoad, payload: str) -> None:
        """戻り値を置いて受付を終える（`done` は必ず立てる）。"""
        try:
            with self._lock:
                entry.payload = payload
                entry.finished_at = self._clock()
        finally:
            entry.done.set()

    def _lookup(self, load_id: str) -> PendingLoad:
        """`self._lock` を持って呼ぶ。"""
        entry = self._entries.get(load_id)
        if entry is None:
            raise DomainError(
                "DATASET_LOAD_NOT_FOUND",
                f"受付番号 {load_id!r} を知りません（終わってから 1 時間を過ぎた、件数の上限で"
                "古いものから捨てた、サーバが再起動した、のどれか）。読み込みが成功していたなら"
                "ハンドルは開いているので dataset_list で確かめ、無ければ load_dataset で"
                "読み直してください。",
                {"load_id": load_id})
        return entry

    def _prune(self) -> None:
        """終わってから期限を過ぎた受付を捨てる。`self._lock` を持って呼ぶ。"""
        now = self._clock()
        expired = [load_id for load_id, entry in self._entries.items()
                   if entry.finished_at is not None and now - entry.finished_at >= self._ttl]
        for load_id in expired:
            del self._entries[load_id]

    def _make_room(self) -> None:
        """1 件足せるように、終わった受付を終わった順に古いものから捨てる。`self._lock` を持って呼ぶ。"""
        while len(self._entries) >= self._max:
            finished = [entry for entry in self._entries.values() if entry.finished_at is not None]
            if not finished:
                raise DomainError(
                    "DATASET_LOAD_LIMIT",
                    f"読み込み中の受付が上限（{self._max}）に達しています。dataset_load_status で"
                    "どれかが終わるのを待ってから load_dataset を呼び直してください。",
                    {"loading": sorted(self._entries), "max_entries": self._max})
            oldest = min(finished, key=lambda entry: entry.finished_at)
            del self._entries[oldest.load_id]


def _domain_envelope(exc: DomainError) -> str:
    """`DomainError` を `load_dataset` と同じ封筒にする。

    `details` が JSON にならない（`Path` など）・`code` / `message` が空のときも落ちない:
    `details` を落とした封筒、それも作れなければ `DATASET_LOAD_FAILED`。
    """
    try:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except Exception:  # noqa: BLE001 - 受付を読み込み中のまま残さない
        _logger.exception("DomainError（%s）を封筒にできなかった", getattr(exc, "code", None))
    try:
        return mztab_error(exc.code, exc.message, {"details_dropped": True})
    except Exception:  # noqa: BLE001
        return failure_envelope(exc)


def failure_envelope(exc: BaseException) -> str:
    """読み込みが `DomainError` 以外で止まったときの封筒（例外の文言は出さず、種類だけ）。"""
    return mztab_error(
        "DATASET_LOAD_FAILED",
        f"読み込みが予期しない誤り（{type(exc).__name__}）で止まりました。詳細はサーバのログに"
        "あります。dataset_list で開いたハンドルを確かめ、入力を見直して load_dataset を"
        "呼び直してください。",
        {"exception": type(exc).__name__})


registry = PendingLoads()
