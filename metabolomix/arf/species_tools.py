"""ARF の分子種ごとの図（arf_plot_species）と分子種 PCA（arf_pca_species）の MCP ツール。

項目の解決は arf/selection.py、描画は plots/species.py / plots/pca_scores.py。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §4・§5。
"""
from __future__ import annotations

import numpy as np
from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset
from metabolomix.analysis.result_state import array_fingerprint, new_provenance
from metabolomix.arf import selection as arf_selection
from metabolomix.arf.tools import _check_ncols, _require_arf_with_arf2
from metabolomix.core import session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.curation import flags as curation_flags
from metabolomix.plots import render as plot_render
from metabolomix.plots.pca_scores import render_pca_scores
from metabolomix.plots import species as species_plot

__all__ = ["arf_plot_species", "arf_pca_species"]

_EXCLUDED_LABELS = (("internal_standard", "内部標準"), ("manual", "手動除外"), ("curation", "判断"),
                    ("auto_likely_wrong", "自動判定"), ("standard_only", "標準液"), ("no_msms", "MS/MS なし"))


def _excluded_summary(excluded: dict) -> str:
    return "・".join(f"{label} {len(excluded.get(key, []))}" for key, label in _EXCLUDED_LABELS)


def _species_caption(payload: dict) -> str:
    groups = "、".join(f"{g['label']} n={len(g['samples'])}" for g in payload["groups"])
    kind = "割合 %" if payload["value"] == "share" else "PeakHeight"
    caption = (f"分子種ごとの{kind}: {len(payload['spots'])} パネル（{groups}）。"
               f"除外 — {_excluded_summary(payload['excluded'])}。")
    if payload["zero_denominator_samples"]:
        caption += f" 分母が 0 で描かなかった試料: {'、'.join(payload['zero_denominator_samples'])}。"
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_plot_species(
    items: list[str],
    groups: list[str],
    value: str = "share",
    share_basis: list[str] | None = None,
    low_reliability_samples: list[str] | None = None,
    apply_curation: bool = True,
    exclude_auto_likely_wrong: bool = False,
    standard_samples: list[str] | None = None,
    require_msms: bool = False,
    title: str | None = None,
    ncols: int | None = None,
    output: str | None = None,
) -> list | str:
    """選んだクラスを分子種（アラインメントスポット）に展開し、1 分子種 1 パネルで群ごとの試料別の値を並べる。

    - items: クラス名（.arf2 の Ontology に完全一致）は配下のスポットに展開、それ以外は分子種名。
      同じ名前でも付加イオンが違えば別パネル（見出しに付加イオンと ID）。1 回 40 パネルまで。
    - groups: arf_plot_group_intensity と同じトークン規則（QC・ブランクは qc / blank を書いたときだけ）。
    - value: "share"（既定。試料ごとに 高さ ÷ 分母の合計 × 100。棒 = 群平均、ひげ = SD）/ "height"
      （PeakHeight を log10 で。横線 = log10 の平均、ひげ = SD、0 は軸の下端）。検定はしない。
    - share_basis: 割合の分母にする項目（省略時は items）。例: PG のページでも分母は
      ["PG", "DGDG", "MGDG", "CL"]。分母にも同じ除外を適用する。
    - 除外は arf_plot_group_intensity と同じ（内部標準・standard_samples だけにある分子種・
      apply_curation の判断・exclude_auto_likely_wrong・arf_exclude したスポット）。require_msms=True で
      MS/MS の裏付けがない分子種も除く。除いたものは payload の excluded。
    - low_reliability_samples: 白抜きで描き、平均と SD から外す。菱形は gap-filled の値。
    - output: "image"（既定）/ "payload"（lipidmix.species_intensity.v1）。PNG ファイルは
      save_figure(kind="species")。
    - 入力の誤りは {"status": "error", "message": ...}。失敗した呼び出しの後は前回の図を破棄する。
    """
    arf2_path = _require_arf_with_arf2()
    if isinstance(arf2_path, str):
        return arf2_path
    arf_state = session_state.session.arf
    arf_state.last_species_plot = None
    filters = dict(low_reliability_samples=low_reliability_samples, apply_curation=apply_curation,
                   exclude_auto_likely_wrong=exclude_auto_likely_wrong, standard_samples=standard_samples)
    try:
        mode = plot_render.resolve_plot_output(output)
        if value not in species_plot.VALUES:
            raise ValueError(f"value は share / height のどちらかを指定してください（受け取った値: {value!r}）。")
        _check_ncols(ncols)
        selected = arf_selection.build_selection(arf_state, arf2_path, items=items, groups=groups, **filters)
        spots, no_msms = arf_selection.expand_spots(selected, require_msms=require_msms)
        if not spots:
            raise ValueError("描く分子種がありません（項目が当たらないか、すべて除外されました）。")
        if len(spots) > species_plot.MAX_PANELS:
            raise ValueError(f"パネルが {len(spots)} 枚になります（1 回 {species_plot.MAX_PANELS} 枚まで）。"
                             "items を分けて呼んでください（share_basis を同じにすれば割合はそろいます）。")
        basis_spots = None
        caveats = list(selected.caveats)
        if value == "share" and share_basis:
            basis_sel = arf_selection.build_selection(arf_state, arf2_path, items=share_basis, groups=groups, **filters)
            basis_spots, _ = arf_selection.expand_spots(basis_sel, require_msms=require_msms)
            outside = sorted({s["spot_id"] for s in spots} - {s["spot_id"] for s in basis_spots})
            if outside:
                caveats.append(f"分母に含まれない分子種が {len(outside)} 個あります（share_basis の外。ID {outside}）。"
                               "割合が 100 % を超えることがあります。")
    except curation_flags.FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    excluded = {**selected.excluded, "no_msms": no_msms}
    payload = species_plot.build_species_payload(
        spots, selected.groups, selected.rows_by_spot, value=value, basis_spots=basis_spots,
        basis_items=share_basis if value == "share" else None, low_reliability=selected.low_reliability,
        excluded=excluded, caveats=caveats)
    if mode == plot_render.PAYLOAD:
        arf_state.last_species_plot = {"payload": payload, "title": title, "ncols": ncols}
        return json_payload(payload)
    png = plot_render.figure_to_png(species_plot.render_species_plot(payload, title=title, ncols=ncols))
    arf_state.last_species_plot = {"payload": payload, "title": title, "ncols": ncols}
    return [_species_caption(payload), Image(data=png, format="png")]


_PCA_COMPONENTS = 5


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_pca_species(
    items: list[str],
    groups: list[str],
    low_reliability_samples: list[str] | None = None,
    apply_curation: bool = True,
    exclude_auto_likely_wrong: bool = False,
    standard_samples: list[str] | None = None,
    require_msms: bool = False,
    normalize: str = "none",
    log_transform: bool = True,
    orient_by: str | None = None,
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """選んだ分子種（arf_plot_species と同じ項目指定・除外）で PCA を回し、スコア図を返す。

    - 対象の試料は groups に入る試料（2 群に入る試料は最初の群の色で 1 回だけ）。
    - low_reliability_samples は主成分の計算に使わず、投影だけする（図では白抜き）。
    - normalize: "none"（既定）/ "total"（試料ごとに選んだ分子種の合計で割り、全試料の合計の中央値を掛ける。
      試料量の差を除いて組成を比べるとき）。log_transform（既定 True）で log10(x + 1)。
      その後、計算に使う試料の平均と母標準偏差で autoscale。分散 0 の分子種は外す。
    - orient_by: 群名。その群の平均スコアが正になるよう主成分の符号をそろえる（符号は本来任意）。
    - ローディング（固有ベクトルの成分と相関 r）はセッションに残し、plot_pca_loadings(source="species")
      が描く。スコア図のファイルは save_figure(kind="pca", source="species")。
    - output: "image"（既定）/ "payload"（点列・寄与率・PC1 の r の上位と下位）。
    - 計算に使う試料が 3 未満などの入力の誤りは {"status": "error", "message": ...}。
    """
    arf2_path = _require_arf_with_arf2()
    if isinstance(arf2_path, str):
        return arf2_path
    arf_state = session_state.session.arf
    arf_state.last_species_pca = None
    settings = {"items": list(items), "groups": list(groups), "normalize": normalize,
                "log_transform": log_transform, "require_msms": require_msms, "orient_by": orient_by,
                "low_reliability_samples": list(low_reliability_samples or [])}
    try:
        mode = plot_render.resolve_plot_output(output)
        if normalize not in ("none", "total"):
            raise ValueError(f"normalize は none / total のどちらかを指定してください（受け取った値: {normalize!r}）。")
        selected = arf_selection.build_selection(
            arf_state, arf2_path, items=items, groups=groups, low_reliability_samples=low_reliability_samples,
            apply_curation=apply_curation, exclude_auto_likely_wrong=exclude_auto_likely_wrong,
            standard_samples=standard_samples)
        labels = [g["label"] for g in selected.groups]
        if orient_by is not None and orient_by not in labels:
            raise ValueError(f"orient_by は groups のどれかを指定してください（{labels}）。")
        spots, no_msms = arf_selection.expand_spots(selected, require_msms=require_msms)
        if len(spots) < 2:
            raise ValueError(f"PCA には分子種が 2 つ以上必要です（現在: {len(spots)}）。")
        samples, group_of = [], {}
        for group in selected.groups:
            for name in group["samples"]:
                if name not in group_of:
                    group_of[name] = group["label"]
                    samples.append(name)
        index = {s["spot_id"]: {r.get("file_name"): r for r in selected.rows_by_spot.get(s["spot_id"], [])}
                 for s in spots}
        X = np.array([[float((index[s["spot_id"]].get(n) or {}).get("height") or 0.0) for s in spots]
                      for n in samples])
        if normalize == "total":
            totals = X.sum(axis=1)
            zero = [n for n, t in zip(samples, totals) if t <= 0]
            if zero:
                raise ValueError(f"選んだ分子種の合計が 0 の試料があり、合計で割れません: {'、'.join(zero)}")
            X = X / totals[:, None] * float(np.median(totals))
        if log_transform:
            X = np.log10(X + 1.0)
        fit = np.array([n not in selected.low_reliability for n in samples])
        result = run_pca_fit_subset(X, fit, n_components=_PCA_COMPONENTS)
    except curation_flags.FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})

    scores = np.asarray(result["components"])
    coef = np.asarray(result["loadings"])
    if orient_by is not None:
        mask = np.array([group_of[n] == orient_by for n in samples]) & fit
        for k in range(coef.shape[0]):
            if scores[mask, k].mean() < 0:
                scores[:, k] *= -1
                coef[k] *= -1
    r = np.asarray(loading_correlations(coef, result["singular_values"], result["n_fit"]))
    kept = [spots[i] for i in result["kept_features"]]
    evr = result["explained_variance_ratio"]
    loadings_rows = [{
        "feature_id": str(s["spot_id"]), "label": s["label"], "ontology": s["ontology"], "adduct": s["adduct"],
        "mz": s["mz"], "rt": s["rt"],
        "coefficient": [round(float(coef[k][i]), 4) for k in range(coef.shape[0])],
        "r": [round(float(r[k][i]), 4) for k in range(coef.shape[0])],
    } for i, s in enumerate(kept)]
    plot = {
        "title": title or "PCA of selected species",
        "x_label": f"PC1 ({evr[0] * 100:.1f}%)",
        "y_label": f"PC2 ({evr[1] * 100:.1f}%)" if len(evr) > 1 else "PC2",
        "points": [{"x": round(float(scores[i][0]), 4),
                    "y": round(float(scores[i][1]), 4) if scores.shape[1] > 1 else 0.0,
                    "label": n, "group": group_of[n], "fitted": bool(fit[i])} for i, n in enumerate(samples)],
        "explained_variance_ratio": [round(float(v), 4) for v in evr],
        "loadings_rows": loadings_rows, "scaling": "autoscale", "n_fit": result["n_fit"],
        "n_species": len(kept), "dropped_zero_variance": len(spots) - len(kept),
        "excluded": {**selected.excluded, "no_msms": no_msms}, "caveats": list(selected.caveats),
        "settings": settings,
        "provenance": new_provenance(arf_state, kind="pca", input_fingerprint=array_fingerprint(X),
                                     effective_parameters=settings),
    }
    ranked = sorted(loadings_rows, key=lambda row: row["r"][0], reverse=True)
    top = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked[:5]]
    bottom = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked[-5:][::-1]]
    caption = (f"分子種 PCA: {len(kept)} 分子種、計算に使った試料 {result['n_fit']}"
               f"（投影 {int((~fit).sum())}）。PC1 {evr[0] * 100:.1f}%"
               + (f"、PC2 {evr[1] * 100:.1f}%" if len(evr) > 1 else "") + "。"
               f" PC1 の r 上位: {', '.join(top)}。下位: {', '.join(bottom)}。"
               + (f" 分散 0 で外した分子種 {plot['dropped_zero_variance']}。" if plot["dropped_zero_variance"] else "")
               + " ローディング図は plot_pca_loadings(source=\"species\")。")
    if mode == plot_render.PAYLOAD:
        arf_state.last_species_pca = plot
        return json_payload({k: plot[k] for k in ("title", "x_label", "y_label", "points",
                                                  "explained_variance_ratio", "n_species", "n_fit",
                                                  "dropped_zero_variance", "excluded", "caveats")}
                            | {"pc1_r_top": top, "pc1_r_bottom": bottom,
                               "result_id": plot["provenance"].get("result_id")})
    png = plot_render.figure_to_png(render_pca_scores(plot, title=plot["title"]))
    arf_state.last_species_pca = plot
    return [caption, Image(data=png, format="png")]
