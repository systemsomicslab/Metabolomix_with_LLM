"""matplotlib の Figure を PNG バイト列にする共通処理。

MCP は画像を content ブロック（`ImageContent`）として返せる。座標点列を LLM の
文脈へ流し込むより桁違いに安く、実測で volcano は 183,578 字（≒数万トークン）に
対し PNG は約 470 画像トークンだった。画像トークンはおおよそ `幅×高さ/750` なので、
**dpi を上げるとトークンが二乗で増える**点に注意（既定 100 を保つこと）。

図の保存（`save_*_figure`）と画像返し（`arf_plot_volcano` 等）の両方がここを通る。
"""
from __future__ import annotations

import io
import os

# 描画系ツールの戻り値の形。
IMAGE = "image"
PAYLOAD = "payload"

# デプロイ単位の既定（`LIPIDMIX_CAVEAT_MODE` と同じ流儀）。チャット系クライアント
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


def figure_to_png(fig, dpi: int = DEFAULT_IMAGE_DPI) -> bytes:
    """Figure を PNG バイト列にして Figure を閉じる。"""
    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    try:
        fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    finally:
        plt.close(fig)
    return buffer.getvalue()
