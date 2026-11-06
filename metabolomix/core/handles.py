"""データセットのハンドル（spec 2026-10-22 §4 の 1）。

対話のセッション（`session_state.session`）はデータセットとライブラリを 1 つずつしか
持てない。並行の再解析のエージェントが別の極性・データセットを同時に扱うと上書きし合う
ので、1 つのサーバのプロセスにハンドル（ID）ごとの状態を並べて持つ。

`DatasetHandle` は `dataset` と `library` を `AnalysisSession` と**同じ名前**で持つ。
ツールは `holder_for(dataset_handle)` が返したものの `.dataset` / `.library` を読むだけで、
ハンドルを渡さない呼び出しは `session_state.session` に効く。

**読み込んだデータセットはすべてここに住む**（spec 2026-11-06 §5）。`dataset_load` も
`session.dataset = ds` でハンドルを 1 本足し、対話のセッションはその直近のハンドル
（`session.current_handle_id`）を指すだけになる。`dataset_handle_open` で開いたハンドルは
直近にならないので、そのハンドルの状態は対話の呼び出しから見えない。同じ入力を開き直した
ときは読み込み済みの基（DatasetState）を使い回すが、ハンドルには複製を渡す——片方の
前処理がもう片方に見えないようにするため。

同じハンドルへの呼び出しは `lock` で直列化する（`core.handle_dispatch` が取る。ハンドルを
省いた対話の呼び出しも、直近のハンドルのロックを取る）。

依存: session_state（LibraryState）と atomic_io だけ。tools は import しない。
"""
from __future__ import annotations

import itertools
import os
import threading
import uuid
from collections import OrderedDict

from metabolomix.core import session_state
from metabolomix.core.atomic_io import DomainError

__all__ = [
    "DEFAULT_MAX_HANDLES",
    "DatasetHandle",
    "HandleRegistry",
    "MAX_HANDLES_ENV",
    "holder_for",
    "registry",
]

MAX_HANDLES_ENV = "LIPIDMIX_MAX_DATASET_HANDLES"
DEFAULT_MAX_HANDLES = 8
#: 読み込み済みの基を何本まで覚えておくか（LRU）。mzTab-M の基は 1 本数十 MB 程度だが、解析フォルダの
#: 大きな基はそれでは済まない（`per_sample_mz` だけで特徴 × 試料の数百万要素。`dataset_load` は
#: 基と開くたびの複製を持つ）。大きさを見た上限は B-2 で入れる（ここは本数だけ）。
MAX_BASES = 4


def max_handles() -> int:
    raw = (os.environ.get(MAX_HANDLES_ENV) or "").strip()
    if not raw:
        return DEFAULT_MAX_HANDLES
    try:
        value = int(raw)
    except ValueError:
        return DEFAULT_MAX_HANDLES
    return max(1, value)


class DatasetHandle:
    """1 つのハンドルの状態。`dataset` / `library` は AnalysisSession と同じ名前。"""

    def __init__(self, dataset, *, polarity, polarity_source, source: dict,
                 arf2_path=None, dcl_path=None, library=None):
        self.handle_id = f"dh_{uuid.uuid4().hex[:12]}"
        self.dataset = dataset
        self.library = library if library is not None else session_state.LibraryState()
        self.lock = threading.RLock()
        self.polarity = polarity
        self.polarity_source = polarity_source
        self.source = dict(source)
        self.arf2_path = arf2_path
        self.dcl_path = dcl_path
        self.created_seq = 0
        self.last_used_seq = 0

    def close(self) -> None:
        store = self.library.store
        if store is not None:
            store.close()
        self.library.store = None
        self.library.last_match = None


class HandleRegistry:
    """ハンドルの登録簿。登録・参照・解放はプロセス内の 1 本のロックで守る。"""

    def __init__(self):
        self._lock = threading.Lock()
        self._handles: dict[str, DatasetHandle] = {}
        self._bases: OrderedDict = OrderedDict()
        self._seq = itertools.count(1)

    def _open_ids(self) -> list[str]:
        return sorted(self._handles)

    def _refuse_if_full(self) -> None:
        """上限なら `DATASET_HANDLE_LIMIT`。`self._lock` を持って呼ぶ。"""
        limit = max_handles()
        if len(self._handles) >= limit:
            raise DomainError(
                "DATASET_HANDLE_LIMIT",
                f"開いているハンドルが上限（{limit}）に達しています。使い終えたハンドルを "
                f"dataset_handle_close で閉じてください（上限は {MAX_HANDLES_ENV} で変えられます）。",
                {"open_handles": self._open_ids(), "max_handles": limit})

    def check_capacity(self) -> None:
        """もう 1 本足せるかを読み込みの**前に**確かめる（上限なら `DATASET_HANDLE_LIMIT`）。

        長い読み込みを済ませてから断らないための前置き。読み込みの間に別の呼び出しが
        足しうるので、`add` も同じ検査をもう一度する。
        """
        with self._lock:
            self._refuse_if_full()

    def add(self, handle: DatasetHandle) -> None:
        with self._lock:
            self._refuse_if_full()
            seq = next(self._seq)
            handle.created_seq = handle.last_used_seq = seq
            self._handles[handle.handle_id] = handle

    def get(self, handle_id: str) -> DatasetHandle:
        with self._lock:
            handle = self._handles.get(handle_id)
            if handle is None:
                raise DomainError(
                    "DATASET_HANDLE_NOT_FOUND",
                    f"ハンドル {handle_id!r} は開いていません（閉じたか、サーバが再起動した）。"
                    "dataset_handle_open で開き直してください。対話のセッションの状態には"
                    "落としていません。",
                    {"dataset_handle": handle_id, "open_handles": self._open_ids()})
            handle.last_used_seq = next(self._seq)
            return handle

    def peek(self, handle_id) -> DatasetHandle | None:
        """使用の記録を動かさずに引く（無ければ None）。"""
        with self._lock:
            return self._handles.get(handle_id)

    def close(self, handle_id: str) -> None:
        with self._lock:
            handle = self._handles.pop(handle_id, None)
        if handle is None:
            raise DomainError(
                "DATASET_HANDLE_NOT_FOUND",
                f"ハンドル {handle_id!r} は開いていません。",
                {"dataset_handle": handle_id, "open_handles": self._open_ids()})
        with handle.lock:
            handle.close()

    def list(self) -> list[DatasetHandle]:
        with self._lock:
            return sorted(self._handles.values(), key=lambda h: h.created_seq)

    # ---------- 読み込み済みの基 ----------

    def cached_base(self, key: str):
        """`(dataset, signature)` を返す（無ければ None）。指紋の照合は呼び出し側。"""
        with self._lock:
            entry = self._bases.get(key)
            if entry is not None:
                self._bases.move_to_end(key)
            return entry

    def put_base(self, key: str, dataset, signature) -> None:
        with self._lock:
            self._bases[key] = (dataset, signature)
            self._bases.move_to_end(key)
            while len(self._bases) > MAX_BASES:
                self._bases.popitem(last=False)


registry = HandleRegistry()


def holder_for(dataset_handle: str | None):
    """`None` なら対話のセッション、そうでなければそのハンドル（無ければ DomainError）。

    `registry` は module 修飾で読む（テストが差し替える）。
    """
    if dataset_handle is None:
        return session_state.session
    return registry.get(dataset_handle)
