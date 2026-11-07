"""`dataset_plot` の kind のうち、項目（クラス・分子種）を選んで群ごとに描く 3 つ（B-2b Task 7）。

    group_intensity  lipidmix.group_intensity.v1   項目ごとの群別の試料別強度（旧 arf_plot_group_intensity）
    species          lipidmix.species_intensity.v1 分子種（スポット）ごとの群別の値（旧 arf_plot_species）
    species_pca      lipidmix.species_pca.v1       選んだ分子種の PCA のスコア図（旧 arf_pca_species）

`plot_tools.py` の `_KINDS` へ登録される（`plot_tools` が末尾でこの関数群を引く）。`plot_tools` を読み込みの時に
import しない（循環するため）——記録の組み立てと `_Request` は呼ばれたときに引く。項目の解決・行・除外は
`plots/dataset_rows.py`、数値の組み立ては `plots/group_intensity.py` `plots/species.py`、PCA は
`analysis/pca.run_pca_fit_subset`。これらの図は群を試料の対応表から引いて図に焼き込むので、記録は
`metadata_revision` を持つ（試料の情報が更新されると `save_figure` が古い図として拒む）。
`species_pca` は結果（`ds.results`、来歴の kind は `species_pca`）も登録し、結果から作った図として
`invalidate_results` と `save_figure` の古さの確認に載る。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis import result_state
from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.artifacts import artifact_path
from metabolomix.msdial import analysis_params
from metabolomix.plots import dataset_rows as rows_mod
from metabolomix.plots import group_intensity as gi
from metabolomix.plots import species as species_plot

#: `species_pca` が回す主成分の数。
_PCA_COMPONENTS = 5

#: `value`（species）の値。`intensity` はピークの強度（payload の `value` は形式名の都合で "height"）。
_SPECIES_VALUES = {"share": "share", "intensity": "height"}

_EXCLUDED_LABELS = (("internal_standard", "内部標準"), ("manual", "手動除外"), ("curation", "判断"),
                    ("auto_likely_wrong", "自動判定"), ("standard_only", "標準液"), ("no_msms", "MS/MS なし"))


def _helpers():
    from metabolomix.tools import plot_tools
    return plot_tools


def _bad(message: str, **details) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details)


def _excluded_summary(excluded: dict, keys=None) -> str:
    return "・".join(f"{label} {len(excluded.get(key, []))}" for key, label in _EXCLUDED_LABELS
                     if keys is None or key in keys)


def _common(ds, req):
    """3 つの kind に共通の解決: 項目・群・低信頼・標準液。`(selection, resolved, rows_caveats)`。"""
    items = list(req.items or [])
    if not items:
        raise _bad("items を指定してください（クラス名・分子種名・特徴 ID。`+` で合算）。")
    if not req.groups:
        raise _bad("groups を 1 群以上指定してください（試料名・Class ID・因子のトークン）。")
    samples = rows_mod.resolve_samples(ds, req.groups, req.low_reliability_samples, req.standard_samples)
    plotted = {n for g in samples.groups for n in g["samples"]}
    standard_only = rows_mod.standard_only_spots(ds, samples.standard, plotted)
    caveats = list(samples.caveats)
    measure = getattr(ds, "quantification_measure", None)
    if measure and measure != "peak_height":
        caveats.append(f"値は定量値（{measure}）で、ピークの高さ（PeakHeight）ではありません。")
    resolved = rows_mod.resolve_selection(
        ds, items, apply_curation=True if req.apply_curation is None else req.apply_curation,
        exclude_auto_likely_wrong=bool(req.exclude_auto_likely_wrong),
        standard_only=standard_only, caveats=caveats)
    if ds.detected_mask is None:
        resolved.caveats.append("検出の印（detected_mask）が無いため、gap-filled の値に印を付けていません。")
    return samples, resolved


def _rows(ds, req, resolved, spot_ids):
    fids = [str(ds.feature_ids[resolved.index[sid]]) for sid in spot_ids]
    return rows_mod.rows_by_spot(ds, fids, groups=list(req.groups),
                                 low_reliability_samples=req.low_reliability_samples,
                                 standard_samples=req.standard_samples)


def _revision(ds) -> int:
    return getattr(ds, "metadata_revision", 0)


def _unmatched_caveat(resolved) -> list[str]:
    unmatched = rows_mod.unmatched_parts(resolved)
    return ([f"items のうち当たらなかった（またはすべて除外された）ものがあります: {'、'.join(unmatched)}。"]
            if unmatched else [])


# ---------- group_intensity ----------

def _detection_limit(ds, requested):
    """検出下限の破線。引数が優先、無ければ param ファイルの Minimum peak height。`(値, 出所)`。"""
    if requested is not None:
        if (isinstance(requested, bool) or not isinstance(requested, (int, float))
                or not math.isfinite(requested) or requested <= 0):
            raise _bad(f"detection_limit は 0 より大きい有限の数で指定してください（受け取った値: {requested!r}）。")
        return float(requested), "argument"
    arf2 = artifact_path(ds, "arf2")
    param = analysis_params.find_param_file(arf2) if arf2 else None
    value = analysis_params.read_analysis_params(param)["min_peak_height"] if param else None
    return (value, "param_file") if value else (None, None)


def group_intensity(ds, req) -> dict:
    limit, source = _detection_limit(ds, req.detection_limit)
    samples, resolved = _common(ds, req)
    spot_ids = sorted({s["spot_id"] for item in resolved.items for s in item["spots"]})
    rows = _rows(ds, req, resolved, spot_ids)
    payload = gi.build_group_intensity_payload(
        resolved.items, samples.groups, rows, msms=resolved.msms, low_reliability=samples.low_reliability,
        excluded=resolved.excluded, detection_limit=limit, detection_limit_source=source,
        caveats=resolved.caveats)
    return _helpers()._record("group_intensity", payload, title=req.title,
                              caption=_group_intensity_caption(payload), metadata_revision=_revision(ds))


def _group_intensity_caption(payload: dict) -> str:
    """画像に添える 1 行。図から読めない内訳（試料数・除外数・N.D.）を言葉で残す。"""
    groups = "、".join(f"{g['label']} n={len(g['samples'])}" for g in payload["groups"])
    nd = [it["item"] for it in payload["items"] if not it["detected"]]
    per_item = "、".join(
        f"{it['item']} {it['n_spots']}"
        + (f" (MS/MS {it['n_spots_msms']})" if it.get("n_spots_msms") is not None else "")
        for it in payload["items"])
    caption = (f"群別強度: {len(payload['items'])} 項目（{groups}）。スポット数 — {per_item}。"
               f"除外 — {_excluded_summary(payload['excluded'], gi._EXCLUDED_KEYS)}。")
    if nd:
        caption += f" N.D.: {'、'.join(nd)}。"
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption


# ---------- species ----------

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


def species(ds, req) -> dict:
    value_name = req.value or "share"
    if value_name not in _SPECIES_VALUES:
        raise _bad(f"value は {' / '.join(_SPECIES_VALUES)} のどちらかを指定してください"
                   f"（受け取った値: {value_name!r}）。", value=value_name)
    value = _SPECIES_VALUES[value_name]
    samples, resolved = _common(ds, req)
    spots, no_msms = rows_mod.expand_spots(ds, resolved, require_msms=bool(req.require_msms))
    caveats = list(resolved.caveats) + _unmatched_caveat(resolved)
    if not spots:
        raise _bad("描く分子種がありません（項目が当たらないか、すべて除外されました）。")
    if len(spots) > species_plot.MAX_PANELS:
        raise _bad(f"パネルが {len(spots)} 枚になります（1 回 {species_plot.MAX_PANELS} 枚まで）。"
                   "items を分けて呼んでください（share_basis を同じにすれば割合はそろいます）。")
    basis_spots, spot_ids = None, [s["spot_id"] for s in spots]
    if value == "height" and req.share_basis:
        caveats.append("value=\"intensity\" では share_basis は無視されます（割合を描くときだけ使います）。")
    if value == "share" and req.share_basis:
        basis = rows_mod.resolve_selection(
            ds, list(req.share_basis), apply_curation=True if req.apply_curation is None else req.apply_curation,
            exclude_auto_likely_wrong=bool(req.exclude_auto_likely_wrong),
            standard_only=rows_mod.standard_only_spots(
                ds, samples.standard, {n for g in samples.groups for n in g["samples"]}))
        basis_unmatched = rows_mod.unmatched_parts(basis)
        if basis_unmatched:
            raise _bad("share_basis のうち当たらなかった（またはすべて除外された）ものがあります: "
                       f"{'、'.join(basis_unmatched)}。分母が変わるため描きません。"
                       "名前を直すか share_basis から外してください。")
        basis_spots, _ = rows_mod.expand_spots(ds, basis, require_msms=bool(req.require_msms))
        outside = sorted({s["spot_id"] for s in spots} - {s["spot_id"] for s in basis_spots})
        if outside:
            caveats.append(f"分母に含まれない分子種が {len(outside)} 個あります（share_basis の外。ID {outside}）。"
                           "割合が 100 % を超えることがあります。")
        spot_ids = sorted(set(spot_ids) | {s["spot_id"] for s in basis_spots})
    rows = _rows(ds, req, resolved, spot_ids)
    payload = species_plot.build_species_payload(
        spots, samples.groups, rows, value=value, basis_spots=basis_spots,
        basis_items=req.share_basis if value == "share" else None, low_reliability=samples.low_reliability,
        excluded={**resolved.excluded, "no_msms": no_msms}, caveats=caveats)
    return _helpers()._record("species", payload, title=req.title, caption=_species_caption(payload),
                              metadata_revision=_revision(ds))


# ---------- species_pca ----------

def _rank_text(rows: list[dict]) -> tuple[list[str], list[str]]:
    ranked = sorted(rows, key=lambda row: row["r"][0], reverse=True)
    top = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked if row["r"][0] > 0][:5]
    bottom = [f"{row['label']} {row['r'][0]:+.2f}" for row in ranked[::-1] if row["r"][0] < 0][:5]
    return top, bottom


def species_pca(ds, req) -> dict:
    normalize = req.normalize or "none"
    if normalize not in ("none", "total"):
        raise _bad(f"normalize は none / total のどちらかを指定してください（受け取った値: {normalize!r}）。")
    log_transform = True if req.log_transform is None else bool(req.log_transform)
    samples, resolved = _common(ds, req)
    labels = [g["label"] for g in samples.groups]
    if req.orient_by is not None and req.orient_by not in labels:
        raise _bad(f"orient_by は groups のどれかを指定してください（{labels}）。")
    spots, no_msms = rows_mod.expand_spots(ds, resolved, require_msms=bool(req.require_msms))
    caveats = list(resolved.caveats) + _unmatched_caveat(resolved)
    if len(spots) < 2:
        raise _bad(f"PCA には分子種が 2 つ以上必要です（現在: {len(spots)}）。")
    names, group_of = [], {}
    for group in samples.groups:
        for name in group["samples"]:
            if name not in group_of:
                group_of[name] = group["label"]
                names.append(name)
    column = {n: j for j, n in enumerate(ds.sample_names)}
    matrix = np.asarray(ds.feature_matrix, dtype=float)
    X = np.array([[float(np.nan_to_num(matrix[resolved.index[s["spot_id"]], column[n]], nan=0.0))
                   for s in spots] for n in names])
    if normalize == "total":
        totals = X.sum(axis=1)
        zero = [n for n, t in zip(names, totals) if t <= 0]
        if zero:
            raise _bad(f"選んだ分子種の合計が 0 の試料があり、合計で割れません: {'、'.join(zero)}")
        X = X / totals[:, None] * float(np.median(totals))
    if log_transform:
        X = np.log10(X + 1.0)
    fit = np.array([n not in samples.low_reliability for n in names])
    if req.orient_by is not None and not any(group_of[n] == req.orient_by and f for n, f in zip(names, fit)):
        raise _bad("orient_by の群に主成分の計算に使った試料がありません"
                   f"（{req.orient_by}。低信頼の試料は計算に使わず投影だけします）。")
    try:
        result = run_pca_fit_subset(X, fit, n_components=_PCA_COMPONENTS)
    except ValueError as exc:
        raise _bad(str(exc)) from exc

    scores = np.asarray(result["components"])
    coef = np.asarray(result["loadings"])
    if req.orient_by is not None:
        mask = np.array([group_of[n] == req.orient_by for n in names]) & fit
        for k in range(coef.shape[0]):
            if scores[mask, k].mean() < 0:
                scores[:, k] *= -1
                coef[k] *= -1
    r = np.asarray(loading_correlations(coef, result["singular_values"], result["n_fit"]))
    kept = [spots[i] for i in result["kept_features"]]
    evr = result["explained_variance_ratio"]
    loadings_rows = [{
        "feature_id": str(ds.feature_ids[resolved.index[s["spot_id"]]]), "label": s["label"],
        "ontology": s["ontology"], "adduct": s["adduct"], "mz": s["mz"], "rt": s["rt"],
        "coefficient": [round(float(coef[k][i]), 4) for k in range(coef.shape[0])],
        "r": [round(float(r[k][i]), 4) for k in range(coef.shape[0])],
    } for i, s in enumerate(kept)]
    top, bottom = _rank_text(loadings_rows)
    settings = {"items": list(req.items), "groups": list(req.groups), "normalize": normalize,
                "log_transform": log_transform, "require_msms": bool(req.require_msms),
                "orient_by": req.orient_by,
                "low_reliability_samples": list(req.low_reliability_samples or []),
                "apply_curation": True if req.apply_curation is None else req.apply_curation,
                "exclude_auto_likely_wrong": bool(req.exclude_auto_likely_wrong),
                "standard_samples": list(req.standard_samples or [])}
    provenance = result_state.new_provenance(
        ds, kind="species_pca", input_fingerprint=result_state.array_fingerprint(X),
        effective_parameters=settings, warnings=caveats)
    result_id = provenance["result_id"]
    dropped = len(spots) - len(kept)
    payload = {
        "plot_schema": "lipidmix.species_pca.v1", "plot_type": "scatter",
        "title": req.title or "PCA of selected species",
        "x_label": f"PC1 ({evr[0] * 100:.1f}%)",
        "y_label": f"PC2 ({evr[1] * 100:.1f}%)" if len(evr) > 1 else "PC2",
        "result_id": result_id,
        "points": [{"x": round(float(scores[i][0]), 4),
                    "y": round(float(scores[i][1]), 4) if scores.shape[1] > 1 else 0.0,
                    "label": n, "group": group_of[n], "fitted": bool(fit[i])} for i, n in enumerate(names)],
        "explained_variance_ratio": [round(float(v), 4) for v in evr],
        "scaling": "autoscale", "n_fit": result["n_fit"], "n_species": len(kept),
        "dropped_zero_variance": dropped,
        "excluded": {**resolved.excluded, "no_msms": no_msms},
        "pc1_r_top": top, "pc1_r_bottom": bottom,
        "caveats": caveats, "settings": settings,
    }
    caption = (f"分子種 PCA: {len(kept)} 分子種、計算に使った試料 {result['n_fit']}"
               f"（投影 {int((~fit).sum())}）。PC1 {evr[0] * 100:.1f}%"
               + (f"、PC2 {evr[1] * 100:.1f}%" if len(evr) > 1 else "") + "。"
               f" PC1 の r 上位: {', '.join(top) or 'なし'}。下位: {', '.join(bottom) or 'なし'}。"
               + (f" 分散 0 で外した分子種 {dropped}。" if dropped else "")
               + f" 除外 — {_excluded_summary(payload['excluded'])}。"
               + (" 注意: " + " / ".join(caveats) if caveats else "")
               + f" 結果は result_id={result_id}（ローディングの全行は ds.results に保持）。")
    # 結果に登録するのは成功したときだけ（ここまでの検査で落ちた呼び出しは何も残さない）。
    result_state.register_result(ds, {
        "kind": "species_pca", "scores": payload["points"], "loadings_rows": loadings_rows,
        "explained_variance_ratio": payload["explained_variance_ratio"], "n_fit": result["n_fit"],
        "scaling": "autoscale", "settings": settings, "provenance": provenance})
    return _helpers()._record("species_pca", payload, title=req.title, caption=caption,
                              result_id=result_id, result=ds.results[result_id],
                              metadata_revision=_revision(ds))
