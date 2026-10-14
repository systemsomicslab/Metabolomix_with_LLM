"""PDF の 1 ページ目を PNG にする（チャット表示用）。

画像トークンは 幅×高さ/750。長辺を上限で揃え、A4 横 1600px ≒ 2,400 トークンに抑える。
"""
from __future__ import annotations

from pathlib import Path

MAX_SIDE_PX = 1600


def pdf_first_page_png(pdf_path: Path, max_side: int = MAX_SIDE_PX) -> tuple[bytes, int, int]:
    import fitz  # PyMuPDF。遅延 import: 無くても server.py は起動できる（図を描くときだけ要る）

    with fitz.open(pdf_path) as doc:
        page = doc[0]
        rect = page.rect
        zoom = max_side / max(rect.width, rect.height)
        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        return pix.tobytes("png"), pix.width, pix.height
