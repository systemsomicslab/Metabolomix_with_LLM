"""ファイルの sha256 を (パス・大きさ・更新時刻) で覚える（stdlib だけの leaf）。

解析フォルダの読み込みは同じファイルの sha256 を何度も取る（読み込み・ライブラリの記憶・
ジョブの出力の集計）。大きなファイルを毎回読み直さないよう、キーが同じ間は覚えた値を返す。
キーのどれかが変われば計算し直す。上限 256 件の LRU。
"""
from __future__ import annotations

import hashlib
import os
import threading
from collections import OrderedDict
from pathlib import Path

__all__ = ["sha256_of", "clear"]

_CHUNK = 1 << 20  # 1 MiB
_LIMIT = 256
_MEMO: "OrderedDict[tuple, str]" = OrderedDict()
_LOCK = threading.Lock()


def sha256_of(path: str | Path) -> str:
    """`path` の sha256（16 進）。読めなければ `OSError`（呼び出し側が従来どおり受ける）。"""
    resolved = str(Path(path).resolve())
    stat = os.stat(resolved)
    key = (resolved, stat.st_size, stat.st_mtime_ns)
    with _LOCK:
        hit = _MEMO.get(key)
        if hit is not None:
            _MEMO.move_to_end(key)
            return hit
    digest = hashlib.sha256()
    with open(resolved, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            digest.update(chunk)
    # 読んでいる間に変わったファイルは覚えない（読み始め前の stat と違えば捨てる）。
    after = os.stat(resolved)
    value = digest.hexdigest()
    if (after.st_size, after.st_mtime_ns) == (stat.st_size, stat.st_mtime_ns):
        with _LOCK:
            _MEMO[key] = value
            _MEMO.move_to_end(key)
            while len(_MEMO) > _LIMIT:
                _MEMO.popitem(last=False)
    return value


def clear() -> None:
    """覚えた値をすべて捨てる（テスト用）。"""
    with _LOCK:
        _MEMO.clear()
