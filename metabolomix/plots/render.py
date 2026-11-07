"""matplotlib の Figure を PNG バイト列にする共通処理。

MCP は画像を content ブロック（`ImageContent`）として返せる。座標点列を LLM の
文脈へ流し込むより桁違いに安く、実測で volcano は 183,578 字（≒数万トークン）に
対し PNG は約 470 画像トークンだった。画像トークンはおおよそ `幅×高さ/750` なので、
**dpi を上げるとトークンが二乗で増える**点に注意（既定 100 を保つこと）。

図の保存（`save_figure`）と画像返し（`arf_plot_volcano` 等）の両方がここを通る。
"""
from __future__ import annotations

import io
import os
import tempfile
import threading
from pathlib import Path

import matplotlib

# ヘッドレスで描く。pyplot を読む前に 1 度だけ指定する（描画ごとには指定しない）。
matplotlib.use("Agg")

# 描画系ツールの戻り値の形。
IMAGE = "image"
PAYLOAD = "payload"

# クライアント単位の既定（起動ごとに変わるので設定ファイルでなく環境変数）。チャット系クライアント
# （Claude Desktop / Claude Code）は画像を表示できるので既定は image。Plotly で
# 対話的に描く Use-LLLM 側は、MCP サーバ起動時の env に
# `LIPIDMIX_PLOT_OUTPUT=payload` を置いて点列を受け取る。
PLOT_OUTPUT_ENV = "LIPIDMIX_PLOT_OUTPUT"


def resolve_plot_output(requested: str | None) -> str:
    """ツール引数 > 環境変数 > 既定(image) の順で戻り値の形を決める。"""
    value = requested if requested is not None else os.getenv(PLOT_OUTPUT_ENV, IMAGE)
    value = str(value).strip().lower()
    if value not in (IMAGE, PAYLOAD):
        raise ValueError(
            f"output は {IMAGE!r} か {PAYLOAD!r} を指定してください（受け取った値: {requested!r}）。"
        )
    return value

# 画像トークンは画素数に比例する。1000x750 で約 1,000 トークン。ツール戻り値として
# 妥当な上限をこのあたりに置き、既定 dpi は上げない。
DEFAULT_IMAGE_DPI = 100


# matplotlib の pyplot は図の登録簿（gcf）などの大域状態を持ち、スレッドセーフでない。
# ツールはワーカースレッドから呼ばれるので、図の生成から破棄までを 1 本のロックで直列にする。
# RLock なのは、描画関数の中から別の描画の部品が入れ子で呼ばれてもデッドロックしないため。
PLOT_LOCK = threading.RLock()


def figure_to_png(fig, dpi: int = DEFAULT_IMAGE_DPI) -> bytes:
    """Figure を PNG バイト列にして Figure を閉じる。このモジュールの中だけで使う（外からは `render_png`）。"""
    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    try:
        fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    finally:
        plt.close(fig)
    return buffer.getvalue()


def render_png(build, /, *args, dpi: int = DEFAULT_IMAGE_DPI, **kwargs) -> bytes:
    """`build(*args, **kwargs)` が返す Figure を PNG にする。生成から破棄までを `PLOT_LOCK` の下で行う。"""
    with PLOT_LOCK:
        fig = build(*args, **kwargs)
        return figure_to_png(fig, dpi)


def render_files(build, paths: dict[str, Path], /, *args, dpi: int = DEFAULT_IMAGE_DPI,
                 atomic: bool = False, **kwargs) -> None:
    """`build(*args, **kwargs)` が返す Figure を、形式（"png" / "svg" ...）ごとに `paths` へ保存する。

    `save_figure` と報告書用の図が使う。`atomic=True` は同じ親ディレクトリの一時ファイルへ
    書いてから置換する（描画の途中で落ちても、既にある図を壊れたファイルで置き換えない）。
    """
    import matplotlib.pyplot as plt

    with PLOT_LOCK:
        fig = build(*args, **kwargs)
        try:
            for fmt, path in paths.items():
                path = Path(path)
                if not atomic:
                    fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight")
                    continue
                fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp",
                                                dir=str(path.parent))
                os.close(fd)
                tmp = Path(tmp_name)
                try:
                    fig.savefig(tmp, format=fmt, dpi=dpi, bbox_inches="tight")
                    os.replace(tmp, path)
                except BaseException:
                    tmp.unlink(missing_ok=True)
                    raise
        finally:
            plt.close(fig)
