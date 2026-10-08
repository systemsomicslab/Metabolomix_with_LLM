"""分子種ごとの図: 1 パネル = 1 アラインメントスポット、群ごとに試料別の値を並べる（MCP 非依存）。

value="share": 試料ごとに「スポットの高さ ÷ 分母（share_basis のスポットの高さの合計）× 100」。
  棒 = 群平均、ひげ = SD（n ≥ 2）、点 = 各試料。
value="height": PeakHeight。log10 の平均 ± SD（群別強度と同じ規則）、0 は軸の下端。
縦軸はパネルごと（分子種ごとに桁が違い、群間の変動を見るのが目的のため）。検定はしない。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §4。
"""
from __future__ import annotations

import math
import statistics

from metabolomix.plots.group_intensity import GROUP_COLORS

SPECIES_SCHEMA = "lipidmix.species_intensity.v1"
MAX_PANELS = 40
VALUES = ("share", "height")


def build_species_payload(spots, groups, rows_by_spot, *, value, basis_spots=None, basis_items=None,
                          low_reliability=frozenset(), excluded=None, caveats=None):
    """スポット × 群 × 試料の値を組み立てる（描画しない）。"""
    if value not in VALUES:
        raise ValueError(f"value は {' / '.join(VALUES)} のどちらかを指定してください（受け取った値: {value!r}）。")
    index = {sid: {r.get("file_name"): r for r in rows} for sid, rows in rows_by_spot.items()}
    samples = []
    for group in groups:
        for name in group["samples"]:
            if name not in samples:
                samples.append(name)

    def height(sid, name):
        row = index.get(sid, {}).get(name)
        return float(row.get("height") or 0.0) if row else 0.0

    denominators = {}
    basis_ids = []
    zero = []
    if value == "share":
        basis_ids = [int(s["spot_id"]) for s in (basis_spots if basis_spots is not None else spots)]
        denominators = {name: sum(height(sid, name) for sid in basis_ids) for name in samples}
        zero = [name for name in samples if denominators[name] <= 0]

    out = []
    for spot in spots:
        sid = int(spot["spot_id"])
        out_groups = []
        for group in groups:
            entries = []
            for name in group["samples"]:
                row = index.get(sid, {}).get(name)
                h = height(sid, name)
                share = None
                if value == "share" and denominators.get(name, 0.0) > 0:
                    share = round(h / denominators[name] * 100.0, 4)
                entries.append({"sample": name, "height": round(h, 1), "share": share,
                                "gap_filled": bool(row.get("is_gap_filled")) if row else False,
                                "low_reliability": name in low_reliability})
            if value == "share":
                stats = [e["share"] for e in entries if e["share"] is not None and not e["low_reliability"]]
            else:
                stats = [math.log10(e["height"]) for e in entries if e["height"] > 0 and not e["low_reliability"]]
            out_groups.append({
                "label": group["label"], "samples": entries,
                "mean": round(statistics.mean(stats), 4) if stats else None,
                "sd": round(statistics.stdev(stats), 4) if len(stats) > 1 else None,
                "n_in_stats": len(stats),
            })
        out.append({**spot, "groups": out_groups})
    return {
        "plot_schema": SPECIES_SCHEMA,
        "value": value,
        "stat": "mean ± SD of share (%)" if value == "share" else "mean ± SD of log10(peak height)",
        "share_basis": ({"items": list(basis_items) if basis_items else None, "spot_ids": basis_ids}
                        if value == "share" else None),
        "groups": [{"label": g["label"], "samples": list(g["samples"])} for g in groups],
        "low_reliability_samples": sorted(low_reliability),
        "zero_denominator_samples": zero,
        "spots": out,
        "excluded": {key: list(values) for key, values in (excluded or {}).items()},
        "caveats": list(caveats or []),
    }


def render_species_plot(payload, *, title=None, ncols=None):
    """payload を 1 スポット 1 パネルの図にする。"""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    spots = payload["spots"]
    if not spots:
        raise ValueError("描く分子種がありません。")
    labels = [g["label"] for g in payload["groups"]]
    share = payload["value"] == "share"
    ncols = ncols or min(5, len(spots))
    nrows = math.ceil(len(spots) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.3 * ncols, 2.25 * nrows + 0.7), squeeze=False)
    for ax, spot in zip(axes.flat, spots):
        has_zero = False
        for i, group in enumerate(spot["groups"]):
            color = GROUP_COLORS[i % len(GROUP_COLORS)]
            if group["mean"] is not None:
                if share:
                    ax.bar(i, group["mean"], width=0.62, color=color, alpha=0.55, edgecolor=color, zorder=2)
                else:
                    ax.hlines(group["mean"], i - 0.28, i + 0.28, color="black", linewidth=1.6, zorder=3)
                if group["sd"] is not None:
                    ax.errorbar(i, group["mean"], yerr=group["sd"], color="black", capsize=3, linewidth=0.9, zorder=3)
            n = len(group["samples"])
            for j, s in enumerate(group["samples"]):
                x = i + (j - (n - 1) / 2) * 0.12
                if share:
                    if s["share"] is None:
                        continue
                    y = s["share"]
                elif s["height"] <= 0:
                    has_zero = True
                    ax.scatter(x, 0.0, marker="v", s=18, facecolor="white", edgecolor="#555555", zorder=4)
                    continue
                else:
                    y = math.log10(s["height"])
                ax.scatter(x, y, s=18, marker="D" if s["gap_filled"] else "o",
                           facecolor="white" if s["low_reliability"] else color,
                           edgecolor="black", linewidth=0.5, zorder=4)
        ax.set_title(f"{spot['label']}\n{spot['adduct']} (ID {spot['spot_id']})", fontsize=7, loc="left")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=6)
        ax.set_xlim(-0.6, len(labels) - 0.4)
        ax.tick_params(axis="y", labelsize=6)
        if share:
            ax.set_ylim(bottom=0)
        elif has_zero:
            ax.set_ylim(bottom=-0.3)
        ax.set_ylabel("% of basis" if share else "log10(peak height)", fontsize=6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in list(axes.flat)[len(spots):]:
        ax.axis("off")
    legend = [Patch(facecolor=GROUP_COLORS[i % len(GROUP_COLORS)], alpha=0.7, label=label)
              for i, label in enumerate(labels)]
    legend += [
        Line2D([], [], marker="o", ls="", mfc="#888888", mec="black", label="detected peak"),
        Line2D([], [], marker="D", ls="", mfc="#888888", mec="black", label="gap-filled value"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="black", label="low-reliability (not in mean/SD)"),
    ]
    if not share:
        legend.append(Line2D([], [], marker="v", ls="", mfc="white", mec="#555555", label="zero (axis floor)"))
    fig.legend(handles=legend, loc="lower center", ncol=min(len(legend), 6), fontsize=6.5, frameon=False)
    default = ("Share of each species (bar = mean, whisker = SD)" if share
               else "Peak height of each species (line = mean of log10, whisker = SD)")
    fig.suptitle(title or default, fontsize=10)
    fig.tight_layout(rect=(0, 0.6 / fig.get_figheight(), 1, 0.97))
    return fig
