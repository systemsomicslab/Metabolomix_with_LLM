"""データセットのハンドル（spec 2026-10-22 §4 の 1）。

対話のセッション（`session_state.session`）はデータセットとライブラリを 1 つずつしか
持てない。並行の再解析のエージェントが別の極性・データセットを同時に扱うと上書きし合う
ので、1 つのサーバのプロセスにハンドル（ID）ごとの状態を並べて持つ。

`DatasetHandle` は `dataset` と `library` を `AnalysisSession` と**同じ名前**で持つ。
ツールは `holder_for(dataset_handle)` が返したものの `.dataset` / `.library` を読むだけで、
ハンドルを渡さない呼び出しは `session_state.session` に効く。

**読み込んだデータセットはすべてここに住む**（spec 2026-11-06 §5）。`load_dataset` も
`session.dataset = ds` でハンドルを 1 本足し、対話のセッションはその直近のハンドル
（`session.current_handle_id`）を指すだけになる。`load_dataset` で開いたハンドルは
常に直近になる（`adopt_handle` が登録と直近の切り替えを 1 つのロックの下でやる）。同じ入力を開き直した
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

from metabolomix.core import session_state, user_config
from metabolomix.core.atomic_io import DomainError

__all__ = [
    "DEFAULT_MAX_HANDLES",
    "MAX_BASES",
    "MAX_BASE_BYTES",
    "available_memory_bytes",
    "estimate_state_bytes",
    "DatasetHandle",
    "HandleRegistry",
    "holder_for",
    "registry",
]

#: ハンドルの上限の設定キー（`[server] max_dataset_handles`）と、それが無いときの上限。
MAX_HANDLES_SETTING = "server.max_dataset_handles"
DEFAULT_MAX_HANDLES = 8
#: 読み込み済みの基を何本まで覚えておくか（LRU）。mzTab-M の基は 1 本数十 MB 程度だが、解析フォルダの
#: 大きな基はそれでは済まない（`per_sample_mz` だけで特徴 × 試料の数百万要素。`load_dataset` は
#: 基と開くたびの複製を持つ）。本数の上限がこれで、大きさの上限は `MAX_BASE_BYTES`（両方で追い出す）。
MAX_BASES = 4
#: 覚えておく基の見積もりの合計の上限（バイト。`estimate_state_bytes` の和）。1 つでこれを超える基は
#: キャッシュしない（毎回読み直す）。ハンドルの複製は基と別に持つので、実メモリはこの数倍になりうる。
MAX_BASE_BYTES = 2 * 1024 ** 3

#: `per_sample_mz` の 1 要素（試料名 → m/z の dict の 1 項目）と `feature_metadata` の 1 件の
#: 見積もりの大きさ（バイト。CPython の dict と float / 小さな dict の実測の桁）。
_BYTES_PER_MZ_CELL = 100
_BYTES_PER_FEATURE_META = 1500


def available_memory_bytes() -> int | None:
    """今、使える物理メモリ（バイト）。取れなければ None（呼び出し側は「分からない」として進む）。

    Windows は `GlobalMemoryStatusEx` の `ullAvailPhys`、POSIX は `SC_AVPHYS_PAGES × SC_PAGE_SIZE`。
    """
    try:
        if os.name == "nt":
            import ctypes

            class _MemoryStatus(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            status = _MemoryStatus()
            status.dwLength = ctypes.sizeof(_MemoryStatus)
            if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return None
            value = int(status.ullAvailPhys)
        else:
            value = int(os.sysconf("SC_AVPHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return None
    return value if value > 0 else None


def _array_bytes(value) -> int:
    """配列（と、配列を値に持つ dict / list）の `nbytes` の和。"""
    nbytes = getattr(value, "nbytes", None)
    if isinstance(nbytes, int):
        return nbytes
    if isinstance(value, dict):
        return sum(_array_bytes(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return sum(_array_bytes(v) for v in value)
    return 0


def estimate_state_bytes(ds) -> int:
    """`DatasetState` が抱える大きさの見積もり（バイト）。

    行列（`feature_matrix` `detected_mask` `pp_matrix` と `analysis_matrices` の配列）の `nbytes` に、
    検出された全セルの m/z（`per_sample_mz` の要素数 × 100）と特徴のメタデータ（件数 × 1500）を足す。
    基のキャッシュの上限の判定だけに使う概算で、他の入れ物（同定の候補など）は数えない。
    """
    total = sum(_array_bytes(getattr(ds, name, None))
                for name in ("feature_matrix", "detected_mask", "pp_matrix", "analysis_matrices"))
    per_sample_mz = getattr(ds, "per_sample_mz", None) or {}
    total += sum(len(cells) for cells in per_sample_mz.values()) * _BYTES_PER_MZ_CELL
    total += len(getattr(ds, "feature_metadata", None) or {}) * _BYTES_PER_FEATURE_META
    return total


def max_handles() -> int:
    """開いておけるハンドルの上限（設定ファイルの `[server] max_dataset_handles`。1 未満は 1）。

    設定ファイルが読めなければ `CONFIG_INVALID`（DomainError。`DATASET_HANDLE_LIMIT` と同じく
    読み込みの入口が封筒にする）。黙って既定の上限に落とさない。
    """
    try:
        setting = user_config.get_setting(MAX_HANDLES_SETTING)
    except user_config.ConfigInvalidError as exc:
        raise DomainError(exc.code, exc.message, exc.details()) from exc
    if setting is None:
        return DEFAULT_MAX_HANDLES
    return max(1, int(setting.value))


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
        self._base_sizes: dict[str, int] = {}  # key -> `estimate_state_bytes`（置いた時点）
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
                "dataset_close で閉じてください（上限は "
                f"{user_config.config_file_path().name} の {user_config.config_key(MAX_HANDLES_SETTING)} "
                "で変えられます）。",
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
                    "load_dataset で開き直してください。対話のセッションの状態には"
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

    def handles_for_arf2(self, sha256: str) -> list[DatasetHandle]:
        """開いているハンドルのうち、データセットの `.arf2` の sha256 が `sha256` のもの（開いた順）。

        ハンドルの `.arf2`（`load_dataset(arf2_path=…)` の明示、無ければデータセットの記録）の
        ハッシュは、読み込みが `source_files` に残した値を使い、無ければ `file_digest.sha256_of`
        で測る。読めない・無い `.arf2` のハンドルは数えない。ハンドルのロックは取らない
        （呼び出し側が、ハンドルごとに自分のロックの下で触る）。
        """
        return [h for h in self.list() if _arf2_sha256(h) == sha256]

    # ---------- 読み込み済みの基 ----------

    def cached_base(self, key: str):
        """`(dataset, signature)` を返す（無ければ None）。指紋の照合は呼び出し側。"""
        with self._lock:
            entry = self._bases.get(key)
            if entry is not None:
                self._bases.move_to_end(key)
            return entry

    def put_base(self, key: str, dataset, signature) -> None:
        """基を覚える。件数（`MAX_BASES`）と見積もりの合計（`MAX_BASE_BYTES`）の両方で LRU を追い出す。

        1 つで `MAX_BASE_BYTES` を超える基は入れない（同じキーの古い基は捨てる。既存の他の基は残す）。
        """
        size = estimate_state_bytes(dataset)
        with self._lock:
            self._bases.pop(key, None)
            self._base_sizes.pop(key, None)
            if size > MAX_BASE_BYTES:
                return
            self._bases[key] = (dataset, signature)
            self._base_sizes[key] = size
            while len(self._bases) > MAX_BASES or (
                    len(self._bases) > 1 and sum(self._base_sizes.values()) > MAX_BASE_BYTES):
                evicted, _entry = self._bases.popitem(last=False)
                self._base_sizes.pop(evicted, None)


def _arf2_sha256(handle: DatasetHandle) -> str | None:
    """ハンドルの `.arf2` の sha256（無い・読めなければ None）。"""
    from metabolomix.core import file_digest
    from metabolomix.dataset.artifacts import artifact_path

    ds = handle.dataset
    path = handle.arf2_path or artifact_path(ds, "arf2")
    if not path:
        return None
    recorded = (getattr(ds, "source_files", None) or {}).get(str(path))
    if recorded:
        return recorded
    try:
        return file_digest.sha256_of(path)
    except OSError:
        return None


registry = HandleRegistry()


def holder_for(dataset_handle: str | None):
    """`None` なら対話のセッション、そうでなければそのハンドル（無ければ DomainError）。

    `registry` は module 修飾で読む（テストが差し替える）。
    """
    if dataset_handle is None:
        return session_state.session
    return registry.get(dataset_handle)
