"""Client-neutral volcano plot payload built from a two-group differential result.

`metabolomix.plots.eic` と同じ役割で、MCP に依存しない。`metabolomix.analysis.differential`
が作った点列を、クライアントがそのまま描ける形へ整形する純ロジック層＋
（`render_eic_plot` と同じく）明示的に図が要るときだけ使う matplotlib 描画。
matplotlib は関数内 import に留め、payload 組み立てだけを使う経路には載せない。

点数は特徴量数ぶん（数千〜数万件）になりうるため `ns` 点だけを決定的に間引く。
`up`/`down` を切らないのは、有意点を落とすと図の意味が壊れるため。
"""
from __future__ import annotations

import math
from typing import TypedDict

VOLCANO_PLOT_SCHEMA = "lipidmix.volcano.v1"

DEFAULT_Q_THRESHOLD = 0.05
DEFAULT_LOG2FC_THRESHOLD = 1.0
_SIGNIFICANT = ("up", "down")

# 点座標の丸め桁数。log2FC も -log10(p) も図と解釈に効くのは小数 2〜3 桁までで、
# 4 桁あれば十分。
POINT_DIGITS = 4

# payload モードの既定点数。特徴量数は実データで 714〜19,388 件あり、旧既定の
# 3000 は「ほぼ全点」＝間引きが働かない値だった。up/down は常に全件残るので、
# ns をこの範囲に収めても有意点の情報は落ちない。
DEFAULT_MAX_POINTS = 800


class PlotAxis(TypedDict):
    label: str
    unit: str | None
    scale: str


class PlotAxes(TypedDict):
    x: PlotAxis
    y: PlotAxis


class VolcanoComparison(TypedDict):
    group_a: str
    group_b: str
    n_a: int | None
    n_b: int | None


class VolcanoThresholds(TypedDict):
    q: float
    log2fc: float


class VolcanoPoint(TypedDict):
    feature: str
    log2fc: float
    neg_log10_p: float
    sig: str


class VolcanoSelection(TypedDict):
    total: int
    plotted: int
    significant_total: int
    significant_plotted: int
    ns_total: int
    ns_plotted: int
    max_points: int
    dropped_nonfinite: int


class VolcanoGuides(TypedDict):
    x: list[float]
    y: list[float]


class VolcanoRenderHints(TypedDict):
    mode: str
    color_by: str
    show_legend: bool
    guides: VolcanoGuides
    hover_fields: list[str]


class VolcanoPlotPayload(TypedDict):
    plot_schema: str
    plot_type: str
    title: str
    comparison: VolcanoComparison
    axes: PlotAxes
    thresholds: VolcanoThresholds
    points: list[VolcanoPoint]
    selection: VolcanoSelection
    render_hints: VolcanoRenderHints
    caveats: list[str]


def _is_drawable(point: dict) -> bool:
    """log2fc と -log10 p の両方が有限なら描画できる。"""
    values = (point.get("log2fc"), point.get("neg_log10_p"))
    return all(
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        for value in values
    )


def _subsample(points: list[dict], quota: int) -> list[dict]:
    """先頭から等間隔に最大 `quota` 件を抽出する（決定的・乱数なし）。

    `(index * total) // quota` は `total >= quota` のとき厳密に単調増加するため、
    重複なくちょうど `quota` 件を返す。x/y の分布形を保ったまま点数だけ落とす。
    """
    if quota <= 0 or not points:
        return []
    total = len(points)
    if total <= quota:
        return list(points)
    return [points[(index * total) // quota] for index in range(quota)]


def build_volcano_plot_payload(
    last_differential: dict,
    *,
    max_points: int = DEFAULT_MAX_POINTS,
    title: str | None = None,
) -> VolcanoPlotPayload:
    """`session.arf.last_differential` から描画中立な volcano ペイロードを組み立てる。"""
    if isinstance(max_points, bool) or not isinstance(max_points, int) or max_points < 1:
        raise ValueError("max_points must be a positive integer")

    raw = list(last_differential.get("volcano") or [])
    q_threshold = float(last_differential.get("q_threshold") or DEFAULT_Q_THRESHOLD)
    log2fc_threshold = float(
        last_differential.get("log2fc_threshold") or DEFAULT_LOG2FC_THRESHOLD
    )
    group_a = str(last_differential.get("a") or "")
    group_b = str(last_differential.get("b") or "")

    caveats: list[str] = []
    drawable = [point for point in raw if _is_drawable(point)]
    dropped_nonfinite = len(raw) - len(drawable)
    if dropped_nonfinite:
        caveats.append(
            f"log2fc または p が有限値でない {dropped_nonfinite} 件は描画対象から外しました"
            "（分散0・欠損・片群のみ検出などで検定不能。有意でないという意味ではありません）。"
        )

    significant = [p for p in drawable if p.get("sig") in _SIGNIFICANT]
    ns = [p for p in drawable if p.get("sig") not in _SIGNIFICANT]
    quota = max_points - len(significant)
    if quota <= 0:
        ns_plotted: list[dict] = []
        if ns:
            caveats.append(
                f"有意点 {len(significant)} 件だけで max_points={max_points} に達したため、"
                f"ns 点 {len(ns)} 件は全て省略しました（有意点は間引いていません）。"
            )
    else:
        ns_plotted = _subsample(ns, quota)
        if len(ns_plotted) < len(ns):
            caveats.append(
                f"ns 点は {len(ns)} 件から {len(ns_plotted)} 件へ等間隔で間引きました"
                "（up/down は全件保持）。件数の判断には selection を参照してください。"
            )

    caveats.append(
        "横破線は -log10(q しきい値) の目安であり、各点の y は -log10(p) です"
        "（q と p は別量なので破線位置と有意判定は厳密には一致しません）。"
    )

    # spec §7.4(R16): allow_confounded=true で継続した比較は、その旨をこの図にも
    # 残す。run_comparison が `provenance.comparison.unadjusted_confounded` へ
    # 記録する事実を読むだけで、ここで独自に交絡を判定はしない。
    comparison_prov = ((last_differential.get("provenance") or {}).get("comparison") or {})
    if comparison_prov.get("unadjusted_confounded"):
        caveats.append(
            "allow_confounded=true が明示されたため、群とバッチの完全交絡を"
            "未調整のまま解析を継続した比較です（この図はその補正前の結果です）。"
        )

    # 座標は 4 桁で丸める。float の既定 repr は 17 桁まで出す（-0.17092926812859943）
    # ため、数百〜数千点の payload では桁の大半が図にも解釈にも寄与しない尾数で
    # 埋まる。実測で 1 点あたり 104 字 → 約 70 字。
    points: list[VolcanoPoint] = [
        {
            "feature": str(point.get("feature") or ""),
            "log2fc": round(float(point["log2fc"]), POINT_DIGITS),
            "neg_log10_p": round(float(point["neg_log10_p"]), POINT_DIGITS),
            "sig": str(point.get("sig") or "ns"),
        }
        for point in significant + ns_plotted
    ]

    guide_y = -math.log10(q_threshold) if q_threshold > 0 else 0.0
    return {
        "plot_schema": VOLCANO_PLOT_SCHEMA,
        "plot_type": "scatter",
        "title": title or f"Volcano ({group_a} vs {group_b})",
        "comparison": {
            "group_a": group_a,
            "group_b": group_b,
            "n_a": last_differential.get("n_a"),
            "n_b": last_differential.get("n_b"),
        },
        "axes": {
            "x": {"label": "log2 fold change", "unit": None, "scale": "linear"},
            "y": {"label": "-log10 p", "unit": None, "scale": "linear"},
        },
        "thresholds": {"q": q_threshold, "log2fc": log2fc_threshold},
        "points": points,
        "selection": {
            "total": len(raw),
            "plotted": len(points),
            "significant_total": len(significant),
            "significant_plotted": len(significant),
            "ns_total": len(ns),
            "ns_plotted": len(ns_plotted),
            "max_points": max_points,
            "dropped_nonfinite": dropped_nonfinite,
        },
        "render_hints": {
            "mode": "markers",
            "color_by": "sig",
            "show_legend": True,
            "guides": {
                "x": [-log2fc_threshold, log2fc_threshold],
                "y": [round(guide_y, 4)],
            },
            "hover_fields": ["feature", "log2fc", "neg_log10_p", "sig"],
        },
        "caveats": caveats,
    }


# 配色は Use-LLLM の volcano-plot.js と一致させてある（PNG と画面で色が入れ替わる
# と、同じ図の話をしているのか判別できなくなる）。
SIG_COLORS = {"up": "#c0392b", "down": "#2471a3", "ns": "#95a5a6"}


def render_volcano_plot(last_differential: dict, title: str | None = None):
    """直近の2群差次的解析を volcano の matplotlib Figure にする。

    点列は payload（間引き済み）ではなく `session.arf.last_differential` の**全量**
    から描く。間引きは payload のトークン対策であって、図には不要なため。
    `save_figure(kind="volcano")`（PNG 保存）と `arf_plot_volcano`（画像返し）の共通描画。
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    for sig in ("ns", "up", "down"):
        pts = [p for p in last_differential.get("volcano") or []
               if p["sig"] == sig and _is_drawable(p)]
        if pts:
            ax.scatter([p["log2fc"] for p in pts], [p["neg_log10_p"] for p in pts],
                       s=12, c=SIG_COLORS[sig], label=sig, alpha=0.7)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p")
    ax.set_title(title or f"Volcano ({last_differential.get('a')} vs {last_differential.get('b')})")
    ax.legend()
    return fig
