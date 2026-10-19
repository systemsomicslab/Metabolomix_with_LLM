"""PyMuPDF が「入っているのに読み込めない」ときだけ、要るテストを理由つきで飛ばす。

Windows の Smart App Control（アプリケーション制御）は、評判の無いネイティブ
拡張をファイル単位で拒否する。PyMuPDF では `_mupdf.pyd` が拒否され、
`import fitz` が `ModuleNotFoundError: No module named 'mupdf'` で落ちる
（2026-10-19 実測。Code Integrity のイベント 3077）。機械の設定の問題で、
コードの退行ではないので、PDF を作る・描くテストだけを skip にする。

**入っていない**（`find_spec` が None）ときは skip しない。PyMuPDF は宣言済みの
依存なので、無いのは環境の壊れで、黙って緑にしてはいけない——そのまま
import させて落とす。
"""
from __future__ import annotations

import importlib.util

import pytest


def require_pymupdf():
    """`fitz` を返す。入っているのに読み込めなければ、呼んだテストを skip する。"""
    if importlib.util.find_spec("fitz") is None:
        import fitz  # noqa: F401 - 無いなら ModuleNotFoundError で落とす（skip しない）
    try:
        import fitz
    except ImportError as exc:
        cause = exc.__cause__ or exc.__context__ or exc
        pytest.skip(f"PyMuPDF は入っているがネイティブ拡張を読み込めない"
                    f"（Windows のアプリケーション制御の可能性）: {cause!r}")
    return fitz
