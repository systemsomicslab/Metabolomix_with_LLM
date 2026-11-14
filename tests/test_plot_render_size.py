"""返す PNG（`render_png`）の画素数の上限（B12）。ファイルに書く図（`render_files`）は対象外。"""
import io

import matplotlib

matplotlib.use("Agg")

from PIL import Image as PILImage

from metabolomix.plots import render


def _figure(width: float, height: float):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(width, height))
    ax.plot([0, 1], [0, 1])
    ax.set_title("size")
    return fig


def _size(png: bytes) -> tuple[int, int]:
    with PILImage.open(io.BytesIO(png)) as image:
        return image.size


def test_small_figure_keeps_the_default_dpi():
    fig = _figure(6, 5)
    try:
        assert render.image_dpi(fig) == render.DEFAULT_IMAGE_DPI
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)


def test_large_figure_is_drawn_at_a_lower_dpi_within_the_pixel_cap():
    width, height = _size(render.render_png(_figure, 12, 6))       # 重ね描きの EIC の大きさ
    assert width * height <= render.MAX_IMAGE_PIXELS * 1.05
    assert width < 12 * render.DEFAULT_IMAGE_DPI


def test_dpi_never_drops_below_the_floor():
    import matplotlib.pyplot as plt

    fig = _figure(20, 20)
    try:
        assert render.image_dpi(fig) == render.MIN_IMAGE_DPI
    finally:
        plt.close(fig)


def test_saved_files_are_not_capped(tmp_path):
    out = tmp_path / "big.png"
    render.render_files(_figure, {"png": out}, 12, 6, dpi=120)
    width, height = _size(out.read_bytes())
    assert width * height > render.MAX_IMAGE_PIXELS
