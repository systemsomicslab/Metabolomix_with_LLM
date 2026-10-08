"""PCA のローディング図（plot_pca_loadings）。ARF・分子種 PCA・mzTab の 3 経路を同じ図にする。

どの PCA を描くかは save_figure(kind="pca") と同じ規則で選ぶ（どれも優先せず、曖昧なら止める）。
ARF は上位 N の選び方とスポット情報に既存の arf/reader.py get_pca_loading_features() を使う。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §6。
"""
from __future__ import annotations

import math

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.core import mcp_errors, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.plots import pca_loadings as loadings_plot
from metabolomix.plots import render as plot_render

__all__ = ["plot_pca_loadings"]


def _candidates() -> list[dict]:
    from metabolomix.analysis.result_state import is_current

    out = []
    arf = session_state.session.arf
    plot = getattr(arf, "last_pca_plot", None)
    if plot and plot.get("loadings") and plot.get("feature_names"):
        prov = plot.get("provenance") or {}
        out.append({"source": "arf", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": True, "result": plot})
    species = getattr(arf, "last_species_pca", None)
    if species and species.get("loadings_rows"):
        prov = species.get("provenance") or {}
        out.append({"source": "species", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": True, "result": species})
    ds = getattr(session_state.session, "dataset", None)
    ds_pca = getattr(ds, "last_pca", None) if ds is not None else None
    if ds_pca and ds_pca.get("loadings"):
        prov = ds_pca.get("provenance") or {}
        out.append({"source": "mztab", "result_id": prov.get("result_id"), "dataset_id": prov.get("dataset_id"),
                    "valid": is_current(ds, ds_pca), "result": ds_pca, "ds": ds})
    return out


def _arf_features(plot: dict, pcs: list[int], top_n: int | None) -> tuple[list[dict], bool]:
    from metabolomix.arf.reader import get_pca_loading_features
    from metabolomix.arf.selection import display_name
    from metabolomix.arf.tools import _sibling_arf2_path
    from metabolomix.arf2 import reader as arf2_reader

    names = plot["feature_names"]
    loadings = plot["loadings"]
    if max(pcs) > len(loadings):
        raise ValueError(f"pcs の {max(pcs)} はありません（この PCA の主成分は 1〜{len(loadings)}）。")
    if top_n is None and len(names) > loadings_plot.MAX_ALL_FEATURES:
        raise ValueError(f"特徴量が {len(names)} あり、全件は描けません（{loadings_plot.MAX_ALL_FEATURES} まで）。"
                         "top_n で絞ってください。")
    picked = get_pca_loading_features(
        {"loadings": loadings, "explained_variance_ratio": plot.get("explained_variance_ratio") or []},
        session_state.session.arf.features or [], names, top_n=top_n or len(names), n_pcs=max(pcs))
    arf2_path = _sibling_arf2_path()
    catalog = ({int(r["MasterAlignmentID"]): r for r in arf2_reader.load_catalog(str(arf2_path))}
               if arf2_path else {})
    column_of = {}
    for j, name in enumerate(names):
        parts = name.split("_")
        if len(parts) >= 3 and parts[1].isdigit():
            column_of.setdefault(int(parts[1]), j)
    ids = []
    for pc in picked:
        for item in pc["positive"] + pc["negative"]:
            if item["id"] not in ids:
                ids.append(item["id"])
    r_ok = plot.get("scaling") == "autoscale" and plot.get("singular_values") and plot.get("n_fit")
    factors = ([float(s) / math.sqrt(plot["n_fit"]) for s in plot["singular_values"]] if r_ok else None)
    features = []
    meta = {item["id"]: item for pc in picked for item in pc["positive"] + pc["negative"]}
    for sid in ids:
        j = column_of[sid]
        coef = [float(loadings[k][j]) for k in range(len(loadings))]
        row = catalog.get(sid, {})
        features.append({
            "feature_id": str(sid),
            "label": display_name(row.get("Name") or "") or meta[sid].get("annotation") or f"Spot {sid}",
            "ontology": (row.get("Ontology") or "").strip() or None, "adduct": row.get("AdductType") or "",
            "mz": meta[sid].get("m_z"), "rt": meta[sid].get("rt"), "coefficient": coef,
            "r": [c * factors[k] for k, c in enumerate(coef)] if factors else None})
    return features, bool(factors)


def _mztab_features(chosen: dict) -> tuple[list[dict], bool]:
    last = chosen["result"]
    ds = chosen["ds"]
    loadings = last["loadings"]
    s = last.get("singular_values")
    n = last.get("n_samples")
    factors = [float(v) / math.sqrt(n) for v in s] if s and n else None
    annotations = getattr(ds, "feature_annotations", {}) or {}
    metadata = getattr(ds, "feature_metadata", {}) or {}
    features = []
    for j, fid in enumerate(ds.pp_feature_names):
        coef = [float(loadings[k][j]) for k in range(len(loadings))]
        ann = annotations.get(fid) or {}
        meta = metadata.get(fid) or {}
        features.append({"feature_id": str(fid), "label": ann.get("name") or meta.get("name") or str(fid),
                         "ontology": None, "adduct": "", "mz": meta.get("mz"), "rt": meta.get("rt"),
                         "coefficient": coef,
                         "r": [c * factors[k] for k, c in enumerate(coef)] if factors else None})
    return features, bool(factors)


def _caption(payload: dict) -> str:
    parts = "、".join(f"{p['pc']} {p['explained_pct']}%（{len(p['rows'])} 行）" for p in payload["panels"])
    caption = (f"PCA ローディング（source={payload['source']} result_id={payload['result_id'] or '(なし)'}、"
               f"値={payload['value']}）: {parts}。")
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def plot_pca_loadings(
    source: str = "auto",
    result_id: str | None = None,
    pcs: list[int] | None = None,
    top_n: int | None = 15,
    value: str = "r",
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """保存済みの PCA（arf_parser / arf_pca_preprocessed / arf_pca_species / dataset_pca）のローディングを横棒で描く。

    - source: "auto"（既定）/ "arf" / "species" / "mztab"。auto はどれも優先せず、候補が 2 つ以上あると
      AMBIGUOUS_RESULT_SOURCE で止まる。result_id で特定の結果を名指しできる。
    - pcs: 描く主成分（既定 [1, 2]、1〜3 個）。
    - top_n: 主成分ごとに正の上位 N と負の上位 N（既定 15）。None は全件（60 特徴量まで）で、
      全パネルの行を最初の主成分の順にそろえる。
    - value: "r"（既定。各特徴量と主成分スコアの相関、−1〜1）/ "coefficient"（固有ベクトルの成分）。
      r は autoscale した PCA でだけ出せ、出せないときは coefficient で描いて説明文に書く。
    - 棒の色はクラス（ARF・分子種 PCA は .arf2 の Ontology）。
    - output: "image"（既定）/ "payload"（lipidmix.pca_loadings.v1）。ファイルは save_figure(kind="pca_loadings")。
    """
    session = session_state.session
    session.last_loadings_plot = None
    try:
        mode = plot_render.resolve_plot_output(output)
        pcs = list(pcs) if pcs is not None else [1, 2]
        if value not in ("r", "coefficient"):
            raise ValueError(f"value は r / coefficient のどちらかを指定してください（受け取った値: {value!r}）。")
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    candidates = _candidates()
    if not candidates:
        return mcp_errors.missing_state(
            "pca_result", ["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "dataset_pca"],
            "先に PCA を実行してください（ローディングを持つ PCA 結果がありません）。")
    from metabolomix.plots.result_output import select_result
    try:
        chosen = select_result(candidates, source=source, result_id=result_id)
    except DomainError as exc:
        return mcp_errors.mztab_error(exc.code, exc.message, exc.details or None)
    try:
        result = chosen["result"]
        if chosen["source"] == "arf":
            features, r_ok = _arf_features(result, pcs, top_n)
        elif chosen["source"] == "species":
            features, r_ok = list(result["loadings_rows"]), True
        else:
            features, r_ok = _mztab_features(chosen)
        payload = loadings_plot.build_loadings_payload(
            features, source=chosen["source"], result_id=chosen.get("result_id"),
            explained_variance_ratio=result.get("explained_variance_ratio") or [], pcs=pcs, top_n=top_n,
            value=value, r_available=r_ok)
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    if mode == plot_render.PAYLOAD:
        session.last_loadings_plot = {"payload": payload, "title": title}
        return json_payload(payload)
    png = plot_render.figure_to_png(loadings_plot.render_loadings_plot(payload, title=title))
    session.last_loadings_plot = {"payload": payload, "title": title}
    return [_caption(payload), Image(data=png, format="png")]
