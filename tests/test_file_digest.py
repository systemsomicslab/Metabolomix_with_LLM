"""`core/file_digest`（sha256 の記憶）。"""
from __future__ import annotations

import ast
import builtins
import hashlib
import os
import sys
from pathlib import Path

import pytest

from metabolomix.core import file_digest


@pytest.fixture(autouse=True)
def _clean():
    file_digest.clear()
    yield
    file_digest.clear()


def _count_opens(monkeypatch, target: Path) -> list:
    opened = []
    real = builtins.open

    def counting(file, *args, **kwargs):
        if Path(str(file)).resolve() == target.resolve():
            opened.append(file)
        return real(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", counting)
    return opened


def test_same_file_hashed_once(tmp_path, monkeypatch):
    path = tmp_path / "a.bin"
    path.write_bytes(b"hello")
    opened = _count_opens(monkeypatch, path)
    expected = hashlib.sha256(b"hello").hexdigest()
    assert file_digest.sha256_of(path) == expected
    assert file_digest.sha256_of(str(path)) == expected
    assert len(opened) == 1


def test_changed_size_or_mtime_rehashes(tmp_path):
    path = tmp_path / "a.bin"
    path.write_bytes(b"hello")
    first = file_digest.sha256_of(path)
    path.write_bytes(b"hello!")
    assert file_digest.sha256_of(path) == hashlib.sha256(b"hello!").hexdigest() != first
    # 大きさが同じでも更新時刻が違えば計算し直す。
    stat = path.stat()
    path.write_bytes(b"HELLO!")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
    assert file_digest.sha256_of(path) == hashlib.sha256(b"HELLO!").hexdigest()


def test_unreadable_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        file_digest.sha256_of(tmp_path / "missing")


def test_memo_is_bounded(tmp_path):
    for i in range(file_digest._LIMIT + 5):
        p = tmp_path / f"f{i}"
        p.write_bytes(bytes([i % 256]))
        file_digest.sha256_of(p)
    assert len(file_digest._MEMO) == file_digest._LIMIT


def test_core_file_digest_is_stdlib_only():
    source = Path(file_digest.__file__).read_text(encoding="utf-8")
    roots = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            roots.add(node.module.split(".")[0])
    assert roots - {"__future__"} <= set(sys.stdlib_module_names)
