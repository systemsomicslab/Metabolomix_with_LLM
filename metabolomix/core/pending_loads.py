"""時間のかかる呼び出しの受付（`load_dataset` の受付番号 `ld_…` と、長い道具の受付番号 `job_…`）と、その結果の控え。

Claude のデスクトップアプリの中継は MCP の呼び出しを 36〜60 秒ほどで打ち切る。読み込みがワーカー
スレッドで回っていても（`handle_dispatch.ALWAYS_THREADED`）、応答が遅れれば中継が諦める。
そこで `load_dataset` は読み込みをこの受付のスレッドで始め、期限（`wait_seconds`）の中で
終われば戻り値をそのまま返し、終わらなければ受付番号を返す。クライアントは
`dataset_load_status(load_id)` で結果を引く。長くかかりうる道具（`handle_dispatch.LONG_RUNNING`。
`pathway_project` / `dataset_homolog_check`）も同じ形で、受付番号（`job_…`）を `job_status(job_id)` で引く。

- 受付が持つのは呼び出しの戻り値だけ（`load_dataset` は JSON の文字列。成功の戻り値か誤りの封筒）。
  ハンドルの登録と直近の切り替えは、読み込み自身（`handle_tools._load_now` の `adopt_handle`）が終わる
  ときに 1 回だけする。結果を引いても状態は変わらない。
- 呼び出しがどう終わっても受付は必ず終わる。`DomainError` は `load_dataset` と同じ封筒
  （`mztab_error`。`details` が JSON にならなければ `details` を落とす）、それ以外は
  `DATASET_LOAD_FAILED` / `JOB_FAILED`（例外の文言は置き場所を含みうるので、種類だけを書く。
  詳細はサーバのログ）。封筒づくり自体が落ちても、最後の封筒で終える。
- `submit` は受付そのもの（`PendingLoad`）を返す。受付を出した呼び出しはそれを直に待つので、
  受付簿から捨てられても結果を失わない。
- 終わった受付は終わってから `TTL_SECONDS`（1 時間）覚え、何度引いても同じ結果を返す。
  件数が `MAX_ENTRIES`（64）に達したら、終わったものを終わった順に古いものから捨てる。
  走っているものは捨てない。すべてが走っていれば `DATASET_LOAD_LIMIT` / `JOB_LIMIT` で受け付けない。
- 受付番号を省いた問い合わせのために、まだ結果を渡していない受付（`unclaimed`）を並べられる。
  結果を渡した受付（`claim`）はそこから外れる。受付簿はプロセスに 1 つで、クライアントを区別しない
  （HTTP 待受では別の接続の受付も並ぶ）。

依存: atomic_io・mcp_errors・serialization（どれも leaf）だけ。tools・session は import しない。
"""
from __future__ import annotations

import json
import logging
import math
import threading
import time
import uuid
from collections.abc import Callable

from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import mztab_error
from metabolomix.core.serialization import json_payload

__all__ = ["MAX_ENTRIES", "MAX_WAIT_SECONDS", "TTL_SECONDS", "PendingCalls", "PendingJobs",
           "PendingLoad", "PendingLoads", "clamp_wait", "failure_envelope", "jobs", "registry",
           "with_wait_note"]

#: 終わった受付を覚えておく秒数（終わったときから数える）。
TTL_SECONDS = 3600.0
#: 覚えておく受付の件数の上限（走っているものを含む）。
MAX_ENTRIES = 64
#: `wait_seconds` の上限（秒）。デスクトップアプリの中継が呼び出しを打ち切る 36〜60 秒ほどより短く
#: （試し走りで 36 秒の呼び出しが打ち切られた）。
MAX_WAIT_SECONDS = 30.0

_logger = logging.getLogger(__name__)


def clamp_wait(value, *, maximum: float = MAX_WAIT_SECONDS) -> tuple[float, str | None]:
    """`wait_seconds` を 0〜`maximum` 秒に丸める。(使う秒数, 丸めたときの注記 or None)。

    範囲の外の数（負・上限超え・±無限）は断らずに丸め、丸めたことを注記で返す。数でないもの
    （bool・文字列・NaN）は丸めようがないので `DATASET_BAD_REQUEST`（DomainError）。
    """
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or (isinstance(value, float) and math.isnan(value))):
        raise DomainError("DATASET_BAD_REQUEST",
                          f"wait_seconds は秒数（0〜{maximum:g}）です: {value!r}",
                          {"field": "wait_seconds"})
    try:
        number = float(value)
    except OverflowError:   # 浮動小数に収まらない巨大な整数（JSON の整数は桁数に上限が無い）
        number = math.inf if value > 0 else -math.inf
    used = min(max(number, 0.0), float(maximum))
    if used == value:
        return used, None
    shown = ("-inf" if number < 0 else "inf") if math.isinf(number) else f"{number:g}"
    return used, (f"wait_seconds={shown} は 0〜{maximum:g} 秒の外なので {used:g} 秒に丸めた"
                  "（デスクトップアプリの中継が 36〜60 秒ほどで呼び出しを打ち切るため）。")


def with_wait_note(payload, note: str | None):
    """`clamp_wait` の注記を、JSON の dict の戻り値に `wait_seconds_note` として足す。

    注記が無い・JSON の dict でない戻り値（画像を含む一覧など）はそのまま返す。
    """
    if not note or not isinstance(payload, str):
        return payload
    try:
        value = json.loads(payload)
    except ValueError:
        return payload
    if not isinstance(value, dict):
        return payload
    value["wait_seconds_note"] = note
    return json_payload(value)


class PendingLoad:
    """1 つの受付。`payload` と `finished_at` は受付簿のロックの下で一緒に書く。

    `tool` は受け付けた道具の名前、`label` は一覧に見せる短い説明（入力のフォルダの名前など）。
    `claimed` は結果を呼び出し側へ渡したか（`unclaimed` から外す）。
    """

    __slots__ = ("load_id", "started_at", "finished_at", "payload", "done", "tool", "label",
                 "claimed")

    def __init__(self, load_id: str, started_at: float, *, tool: str | None = None,
                 label: str | None = None):
        self.load_id = load_id
        self.started_at = started_at
        self.finished_at: float | None = None
        self.payload = None
        self.done = threading.Event()
        self.tool = tool
        self.label = label
        self.claimed = False

    @property
    def job_id(self) -> str:
        """受付番号（`load_id` と同じ。長い道具の受付では `job_…`）。"""
        return self.load_id

    def wait(self, seconds: float):
        """最長 `seconds` 秒待つ。終われば戻り値、終わらなければ None。

        受付簿から捨てられた後でも、この受付そのものを待つので結果を失わない。
        """
        timeout = min(max(0.0, float(seconds)), threading.TIMEOUT_MAX)
        if self.done.wait(timeout):
            return self.payload
        return None


class PendingCalls:
    """受付簿の共通部。受付・参照・捨てるのはプロセス内の 1 本のロックで守る。

    `clock` は経過と期限を測る時計（既定 `time.monotonic`。テストが差し替える）。
    受付番号の頭・誤りのコード・失敗の封筒はサブクラスが決める（`PendingLoads` / `PendingJobs`）。
    """

    #: 受付番号の頭（`<prefix>_` + 12 桁）。
    prefix = "call"
    #: 受付のスレッドの名前の頭（ログ用）。
    thread_label = "call"

    def __init__(self, *, clock: Callable[[], float] = time.monotonic,
                 ttl_seconds: float = TTL_SECONDS, max_entries: int = MAX_ENTRIES):
        self._clock = clock
        self._ttl = float(ttl_seconds)
        self._max = max(1, int(max_entries))
        self._lock = threading.Lock()
        self._entries: dict[str, PendingLoad] = {}

    def submit(self, fn: Callable[[], object], *, tool: str | None = None,
               label: str | None = None) -> PendingLoad:
        """`fn` を daemon スレッドで回し、受付を返す。

        上限なら終わった受付を古いものから捨てて空きを作る。すべてが走っていれば上限の誤り
        （DomainError。スレッドは起こさない）。
        """
        with self._lock:
            self._prune()
            self._make_room()
            load_id = f"{self.prefix}_{uuid.uuid4().hex[:12]}"
            entry = PendingLoad(load_id, self._clock(), tool=tool, label=label)
            self._entries[load_id] = entry
        thread = threading.Thread(target=self._run, args=(entry, fn),
                                  name=f"{self.thread_label}-{load_id}", daemon=True)
        try:
            thread.start()
        except BaseException:
            with self._lock:  # 起こせなかった受付を走っているまま残さない
                self._entries.pop(load_id, None)
            raise
        return entry

    def wait(self, load_id: str, seconds: float):
        """受付番号で引いて最長 `seconds` 秒待つ（`PendingLoad.wait`）。知らない番号は
        `_not_found`。0 秒なら待たずに今の状態を見る（既に終わっていれば返す）。"""
        with self._lock:
            self._prune()
            entry = self._lookup(load_id)
        return entry.wait(seconds)

    def status(self, load_id: str) -> tuple[str, object]:
        """`("loading", None)` か `("done", payload)`。知らない受付番号は `_not_found`。"""
        state, payload, _elapsed = self.snapshot(load_id)
        return state, payload

    def snapshot(self, load_id: str) -> tuple[str, object, float]:
        """`status` に受け付けてからの秒（終わったものはかかった秒）を足したもの。

        3 つは同じロックの下で読む（状態と経過が食い違わない）。
        """
        with self._lock:
            self._prune()
            entry = self._lookup(load_id)
            if entry.finished_at is None:
                return "loading", None, self._clock() - entry.started_at
            return "done", entry.payload, entry.finished_at - entry.started_at

    def discard(self, load_id: str) -> None:
        """受付を忘れる（無ければ何もしない）。走っているものも忘れるが、呼び出しは止めない。"""
        with self._lock:
            self._entries.pop(load_id, None)

    def claim(self, load_id: str) -> None:
        """結果を呼び出し側へ渡した印を付ける（`unclaimed` から外す。無ければ何もしない）。"""
        with self._lock:
            entry = self._entries.get(load_id)
            if entry is not None:
                entry.claimed = True

    def tool_of(self, load_id: str) -> str | None:
        """受け付けた道具の名前（知らない番号・名前の無い受付は None）。"""
        with self._lock:
            entry = self._entries.get(load_id)
            return entry.tool if entry is not None else None

    def unclaimed(self) -> list[dict]:
        """まだ結果を渡していない受付（走っているものと、終わって引かれていないもの）を受け付けた順に。

        各要素は `{"id", "state": "loading"|"done", "elapsed_s", "tool", "label"}`。
        """
        with self._lock:
            self._prune()
            now = self._clock()
            entries = sorted((e for e in self._entries.values() if not e.claimed),
                             key=lambda e: e.started_at)
            return [{"id": e.load_id,
                     "state": "loading" if e.finished_at is None else "done",
                     "elapsed_s": round((e.finished_at if e.finished_at is not None else now)
                                        - e.started_at, 1),
                     "tool": e.tool, "label": e.label} for e in entries]

    # ---------- サブクラスが決めるもの ----------

    def _not_found(self, load_id: str) -> DomainError:
        return DomainError("CALL_NOT_FOUND", f"受付番号 {load_id!r} を知りません。",
                           {"id": load_id})

    def _limit(self) -> DomainError:
        return DomainError("CALL_LIMIT", f"走っている受付が上限（{self._max}）に達しています。",
                           {"running": sorted(self._entries), "max_entries": self._max})

    def _failure(self, exc: BaseException):
        return failure_envelope(exc)

    def _last_resort(self):
        return _LAST_RESORT

    def _accepts(self, payload) -> bool:
        """受付に置ける戻り値か（`wait` の None と取り違えないよう、None は置かない）。"""
        return payload is not None

    # ---------- 内部 ----------

    def _run(self, entry: PendingLoad, fn: Callable[[], object]) -> None:
        """受付のスレッドの本体。どの出口でも戻り値を置いて受付を終える。"""
        payload: object = None
        try:
            try:
                payload = fn()
            except DomainError as exc:  # 本体の try の外で出たもの（本体の中のものは封筒で返る）
                payload = _domain_envelope(exc, self._failure)
            except BaseException as exc:  # noqa: BLE001 - 受付を走っているまま残さない
                _logger.exception("%s が予期しない誤りで止まった（%s）",
                                  entry.tool or self.thread_label, entry.load_id)
                payload = self._failure(exc)
            if not self._accepts(payload):  # `wait` の None（終わっていない）と取り違えない
                payload = self._failure(
                    TypeError(f"戻り値の型が受付に置けない: {type(payload).__name__}"))
        except BaseException:  # noqa: BLE001 - 封筒づくりが落ちても受付は終える
            _logger.exception("%s の結果を封筒にできなかった（%s）",
                              entry.tool or self.thread_label, entry.load_id)
        finally:
            self._finish(entry, payload if self._safe_accepts(payload) else self._last_resort())

    def _safe_accepts(self, payload) -> bool:
        try:
            return bool(self._accepts(payload))
        except Exception:  # noqa: BLE001
            return False

    def _finish(self, entry: PendingLoad, payload) -> None:
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
            raise self._not_found(load_id)
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
                raise self._limit()
            oldest = min(finished, key=lambda entry: entry.finished_at)
            del self._entries[oldest.load_id]


class PendingLoads(PendingCalls):
    """`load_dataset` の受付簿（受付番号 `ld_…`）。戻り値は JSON の文字列だけを置く。"""

    prefix = "ld"
    thread_label = "load_dataset"

    def _not_found(self, load_id: str) -> DomainError:
        return DomainError(
            "DATASET_LOAD_NOT_FOUND",
            f"受付番号 {load_id!r} を知りません（終わってから 1 時間を過ぎた、件数の上限で"
            "古いものから捨てた、サーバが再起動した、のどれか）。読み込みが成功していたなら"
            "ハンドルは開いているので dataset_list で確かめ、無ければ load_dataset で"
            "読み直してください。",
            {"load_id": load_id})

    def _limit(self) -> DomainError:
        return DomainError(
            "DATASET_LOAD_LIMIT",
            f"読み込み中の受付が上限（{self._max}）に達しています。dataset_load_status で"
            "どれかが終わるのを待ってから load_dataset を呼び直してください。",
            {"loading": sorted(self._entries), "max_entries": self._max})

    def _accepts(self, payload) -> bool:
        return isinstance(payload, str)


class PendingJobs(PendingCalls):
    """長い道具（`handle_dispatch.LONG_RUNNING`）の受付簿（受付番号 `job_…`）。戻り値は道具のものをそのまま置く。"""

    prefix = "job"
    thread_label = "job"

    def _not_found(self, load_id: str) -> DomainError:
        return DomainError(
            "JOB_NOT_FOUND",
            f"受付番号 {load_id!r} を知りません（終わってから 1 時間を過ぎた、件数の上限で"
            "古いものから捨てた、サーバが再起動した、のどれか）。道具を呼び直してください。",
            {"job_id": load_id})

    def _limit(self) -> DomainError:
        return DomainError(
            "JOB_LIMIT",
            f"走っている受付が上限（{self._max}）に達しています。job_status でどれかが終わるのを"
            "待ってから呼び直してください。",
            {"running": sorted(self._entries), "max_entries": self._max})

    def _failure(self, exc: BaseException):
        return job_failure_envelope(exc)

    def _last_resort(self):
        return _JOB_LAST_RESORT


#: 封筒づくりがどれも落ちたときの最後の封筒（ここは落ちない）。
_LAST_RESORT = json_payload({"error": {
    "code": "DATASET_LOAD_FAILED",
    "message": "読み込みの結果を封筒にできませんでした。詳細はサーバのログにあります。"
               "dataset_list で開いたハンドルを確かめ、load_dataset を呼び直してください。"}})
_JOB_LAST_RESORT = json_payload({"error": {
    "code": "JOB_FAILED",
    "message": "道具の結果を封筒にできませんでした。詳細はサーバのログにあります。"}})


def _domain_envelope(exc: DomainError, failure=None) -> str:
    """`DomainError` を `load_dataset` と同じ封筒にする。

    `details` が JSON にならない（`Path` など）・`code` / `message` が空のときも落ちない:
    `details` を落とした封筒、それも作れなければ失敗の封筒（既定 `DATASET_LOAD_FAILED`）。
    """
    try:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except Exception:  # noqa: BLE001 - 受付を走っているまま残さない
        _logger.exception("DomainError（%s）を封筒にできなかった", getattr(exc, "code", None))
    try:
        return mztab_error(exc.code, exc.message, {"details_dropped": True})
    except Exception:  # noqa: BLE001
        return (failure or failure_envelope)(exc)


def failure_envelope(exc: BaseException) -> str:
    """読み込みが `DomainError` 以外で止まったときの封筒（例外の文言は出さず、種類だけ）。"""
    return mztab_error(
        "DATASET_LOAD_FAILED",
        f"読み込みが予期しない誤り（{type(exc).__name__}）で止まりました。詳細はサーバのログに"
        "あります。dataset_list で開いたハンドルを確かめ、入力を見直して load_dataset を"
        "呼び直してください。",
        {"exception": type(exc).__name__})


def job_failure_envelope(exc: BaseException) -> str:
    """長い道具が `DomainError` 以外で止まったときの封筒（例外の文言は出さず、種類だけ）。"""
    return mztab_error(
        "JOB_FAILED",
        f"道具が予期しない誤り（{type(exc).__name__}）で止まりました。詳細はサーバのログに"
        "あります。入力を見直して呼び直してください。",
        {"exception": type(exc).__name__})


registry = PendingLoads()
jobs = PendingJobs()
