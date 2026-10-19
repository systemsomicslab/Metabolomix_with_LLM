"""`require_pymupdf` が「入っていない」を skip で隠さないこと。"""
from __future__ import annotations

import importlib.util
import sys

import pytest

from tests.pymupdf_guard import require_pymupdf


def test_missing_pymupdf_fails_instead_of_skipping(monkeypatch):
    # 宣言済みの依存が無いのは環境の壊れ。skip にすると黙って緑になる。
    monkeypatch.setattr(importlib.util, "find_spec", lambda name, *a, **k: None)
    monkeypatch.setitem(sys.modules, "fitz", None)
    with pytest.raises(ImportError):
        require_pymupdf()


def test_installed_but_unloadable_pymupdf_skips(monkeypatch):
    real_find_spec = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec",
                        lambda name, *a, **k: object() if name == "fitz" else real_find_spec(name, *a, **k))
    monkeypatch.setitem(sys.modules, "fitz", None)
    with pytest.raises(pytest.skip.Exception, match="PyMuPDF"):
        require_pymupdf()
