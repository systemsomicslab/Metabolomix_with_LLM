"""PCA スコア図（MCP 非依存）。分子種 PCA の画像返しと save_figure(kind="pca") が共有する。

points の `group` があれば色分けと凡例、`fitted=False`（主成分の計算に使わず投影だけした試料）は
白抜き。軸の範囲は計算に使った点に合わせ、そこから外れた投影点は枠の端に矢印と座標で示す
（外れ値 1 点のために他の点が潰れないように）。
"""
from __future__ import annotations

from metabolomix.plots.group_intensity import GROUP_COLORS


def render_pca_scores(plot: dict, *, title: str | None = None):
    import matplotlib.pyplot as plt

    points = plot.get("points") or []
    if not points:
        raise ValueError("PCA 図に描ける点がありません。")
    fitted = [p for p in points if p.get("fitted", True)] or points
    xs = [float(p["x"]) for p in fitted]
    ys = [float(p["y"]) for p in fitted]
    pad_x = 0.15 * (max(xs) - min(xs)) + 0.5
    pad_y = 0.15 * (max(ys) - min(ys)) + 0.5
    lo_x, hi_x, lo_y, hi_y = min(xs) - pad_x, max(xs) + pad_x, min(ys) - pad_y, max(ys) + pad_y
    groups = []
    for p in points:
        g = p.get("group")
        if g is not None and g not in groups:
            groups.append(g)
    color_of = {g: GROUP_COLORS[(i + 1) % len(GROUP_COLORS)] for i, g in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(lo_y, hi_y)
    labelled = set()
    has_projected = False
    for p in points:
        x, y = float(p["x"]), float(p["y"])
        projected = not p.get("fitted", True)
        has_projected = has_projected or projected
        if projected and not (lo_x <= x <= hi_x and lo_y <= y <= hi_y):
            cx, cy = min(max(x, lo_x), hi_x), min(max(y, lo_y), hi_y)
            mx, my = (lo_x + hi_x) / 2, (lo_y + hi_y) / 2
            ax.annotate(f"{p.get('label', '')} (off-scale: {x:.0f}, {y:.0f})", xy=(cx, cy),
                        xytext=(cx - 0.25 * (cx - mx), cy - 0.25 * (cy - my)), fontsize=7, ha="center",
                        arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "#555555"})
            continue
        g = p.get("group")
        label = g if (g is not None and g not in labelled and not projected) else None
        if label:
            labelled.add(g)
        ax.scatter(x, y, s=60, facecolor="white" if projected else color_of.get(g, "#3182bd"),
                   edgecolor="black", linewidth=0.7, zorder=3, label=label)
        if p.get("label"):
            ax.annotate(str(p["label"]), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.axhline(0, color="#eeeeee", lw=0.6, zorder=1)
    ax.axvline(0, color="#eeeeee", lw=0.6, zorder=1)
    ax.set_xlabel(plot.get("x_label", "PC1"))
    ax.set_ylabel(plot.get("y_label", "PC2"))
    ax.set_title(title or plot.get("title", "PCA"))
    if groups or has_projected:
        handles, labels = ax.get_legend_handles_labels()
        if has_projected:
            handles.append(plt.Line2D([], [], marker="o", ls="", mfc="white", mec="black"))
            labels.append("projected (not used to fit)")
        ax.legend(handles, labels, fontsize=7, loc="best")
    fig.tight_layout()
    return fig
