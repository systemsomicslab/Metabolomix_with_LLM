"""PCA ローディング図（MCP 非依存）。ARF・分子種 PCA・mzTab の 3 経路が同じ形に変換して渡す。

1 パネル = 1 主成分の横棒。値は相関 r（既定。−1〜1、autoscale した PCA でだけ出せる）か
固有ベクトルの成分（coefficient）。top_n は主成分ごとに正の上位 N と負の上位 N（重複は 1 回）。
top_n=None は全件で、全パネルの行を最初の主成分の値の順にそろえる（行の並びで見比べられるように）。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §6。
"""
from __future__ import annotations

LOADINGS_SCHEMA = "lipidmix.pca_loadings.v1"
MAX_ALL_FEATURES = 60
MAX_PCS = 3
_CLASS_PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#7e57c2",
                  "#00897b", "#8d6e63", "#5c6bc0", "#c0ca33"]


def _check_pcs(pcs, n_available):
    if not pcs or len(pcs) > MAX_PCS or len(set(pcs)) != len(pcs):
        raise ValueError(f"pcs は重複のない 1〜{MAX_PCS} 個の主成分番号で指定してください（受け取った値: {pcs!r}）。")
    for pc in pcs:
        if isinstance(pc, bool) or not isinstance(pc, int) or not 1 <= pc <= n_available:
            raise ValueError(f"pcs の {pc!r} はありません（この PCA の主成分は 1〜{n_available}）。")


def _row(feature, k, use_r):
    coef = float(feature["coefficient"][k])
    r = float(feature["r"][k]) if feature.get("r") is not None else None
    return {"feature_id": feature["feature_id"], "label": feature["label"], "ontology": feature.get("ontology"),
            "adduct": feature.get("adduct") or "", "mz": feature.get("mz"), "rt": feature.get("rt"),
            "coefficient": round(coef, 4), "r": round(r, 4) if r is not None else None,
            "value": round(r if use_r else coef, 4)}


def build_loadings_payload(features, *, source, result_id, explained_variance_ratio, pcs, top_n, value,
                           r_available, caveats=None):
    """特徴量の一覧から、主成分ごとに描く行を選んで payload にする。"""
    if value not in ("r", "coefficient"):
        raise ValueError(f"value は r / coefficient のどちらかを指定してください（受け取った値: {value!r}）。")
    if top_n is not None and (isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1):
        raise ValueError(f"top_n は 1 以上の整数か None で指定してください（受け取った値: {top_n!r}）。")
    n_available = min(len(explained_variance_ratio), min((len(f["coefficient"]) for f in features), default=0))
    _check_pcs(pcs, n_available)
    caveats = list(caveats or [])
    use_r = value == "r" and r_available
    if value == "r" and not r_available:
        caveats.append("この PCA は相関 r を出せない（autoscale していない、または特異値が無い）ため、"
                       "固有ベクトルの成分（coefficient）で描きました。")
    aligned = top_n is None
    if aligned and len(features) > MAX_ALL_FEATURES:
        raise ValueError(f"特徴量が {len(features)} あり、全件は描けません（{MAX_ALL_FEATURES} まで）。"
                         "top_n で絞ってください。")
    panels = []
    order = None
    for pc in pcs:
        k = pc - 1
        rows = [_row(f, k, use_r) for f in features]
        if aligned:
            if order is None:
                order = [row["feature_id"] for row in sorted(rows, key=lambda r: r["value"], reverse=True)]
            by_id = {row["feature_id"]: row for row in rows}
            chosen = [by_id[fid] for fid in order]
        else:
            ranked = sorted(rows, key=lambda r: r["value"], reverse=True)
            picked = ranked[:top_n] + ranked[-top_n:]
            seen = set()
            chosen = []
            for row in sorted(picked, key=lambda r: r["value"], reverse=True):
                if row["feature_id"] not in seen:
                    seen.add(row["feature_id"])
                    chosen.append(row)
        panels.append({"pc": f"PC{pc}", "explained_pct": round(float(explained_variance_ratio[k]) * 100, 1),
                       "rows": chosen})
    return {"plot_schema": LOADINGS_SCHEMA, "source": source, "result_id": result_id,
            "value": "r" if use_r else "coefficient", "value_requested": value, "top_n": top_n,
            "aligned": aligned, "panels": panels, "caveats": caveats}


def render_loadings_plot(payload, *, title=None):
    """payload を主成分ごとの横棒にする（色はクラス）。"""
    import matplotlib.pyplot as plt

    panels = payload["panels"]
    if not panels or not any(p["rows"] for p in panels):
        raise ValueError("描くローディングがありません。")
    classes = sorted({row["ontology"] for p in panels for row in p["rows"] if row.get("ontology")})
    color_of = {c: _CLASS_PALETTE[i % len(_CLASS_PALETTE)] for i, c in enumerate(classes)}
    n_rows = max(len(p["rows"]) for p in panels)
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 0.19 * n_rows + 1.6),
                             sharey=payload["aligned"], squeeze=False)
    is_r = payload["value"] == "r"
    for ax, panel in zip(axes[0], panels):
        rows = list(reversed(panel["rows"]))                 # 上ほど値が大きい
        ys = range(len(rows))
        ax.barh(list(ys), [row["value"] for row in rows], height=0.72,
                color=[color_of.get(row.get("ontology"), "#9e9e9e") for row in rows])
        ax.axvline(0, color="#6b6a64", lw=0.8)
        if is_r:
            ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel(f"{'loading r' if is_r else 'coefficient'} ({panel['pc']}, {panel['explained_pct']}%)",
                      fontsize=8)
        if not payload["aligned"] or ax is axes[0][0]:
            ax.set_yticks(list(ys))
            ax.set_yticklabels([f"{row['label']} {row['adduct']}".strip() for row in rows], fontsize=6.5)
        ax.grid(axis="x", color="#e4e3dc", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    if classes:
        handles = [plt.Rectangle((0, 0), 1, 1, color=color_of[c]) for c in classes]
        fig.legend(handles, classes, loc="lower center", ncol=min(len(classes), 8), frameon=False, fontsize=7)
    default = f"PCA loadings ({payload['source']}; {'r = correlation with the PC score' if is_r else 'eigenvector coefficient'})"
    fig.suptitle(title or default, fontsize=9)
    fig.tight_layout(rect=(0, 0.45 / fig.get_figheight() if classes else 0, 1, 0.97))
    return fig
