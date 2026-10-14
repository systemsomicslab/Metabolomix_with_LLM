import fitz

from metabolomix.plots import pdf_raster


def _pdf(tmp_path, w=842, h=595):          # A4 横（pt）
    doc = fitz.open()
    page = doc.new_page(width=w, height=h)
    page.insert_text((72, 72), "PC")
    path = tmp_path / "x.pdf"
    doc.save(path)
    return path


def test_long_side_is_capped(tmp_path):
    png, w, h = pdf_raster.pdf_first_page_png(_pdf(tmp_path))
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert max(w, h) == 1600
    assert abs(w / h - 842 / 595) < 0.01


def test_small_cap(tmp_path):
    _, w, h = pdf_raster.pdf_first_page_png(_pdf(tmp_path), max_side=400)
    assert max(w, h) == 400
