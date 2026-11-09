"""描画の MCP ツール `dataset_plot`（spec §8）。`kind` ごとの関数へ振り分ける。

入力は `DatasetState`（ハンドル）。数値・描画の組み立ては `plots/*` の純粋な部品を使い、ここは
結果の選び方（古い結果は描かない）・引数の取り違えの拒否・`.EIC.aef` の出所・戻り値の形だけを持つ。
作った図の記録は `ds.figures[kind]` に置き、`save_figure(kind=...)` がそこから描く。

`kind` と形式名（Use-LLLM が読む）:
    volcano        lipidmix.volcano.v1        直近の差次的結果（2 群）か `result_id`
    pca            lipidmix.pca_scores.v1     直近の `dataset_pca` か `result_id`。`color_by`
    loadings       lipidmix.pca_loadings.v1   同じ PCA。`pcs` `top_n` `value`
    group_intensity  lipidmix.group_intensity.v1   項目（クラス・分子種）ごとの群別の試料別強度
    species        lipidmix.species_intensity.v1 分子種（スポット）ごとの群別の値（割合か強度）
    species_pca    lipidmix.species_pca.v1    選んだ分子種の PCA のスコア図（結果も登録する）
    eic            lipidmix.eic.v1            `feature_id`（= MasterAlignmentID）の試料別クロマトグラム
    eic_compounds  lipidmix.eic.multi.v1      1 試料に同定済みの物質を重ねる

group_intensity / species / species_pca は `plot_species_kinds.py`（項目の解決・行は `plots/dataset_rows.py`）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.analysis import result_state, stats_results
from metabolomix.analysis.dataset_analysis import included_sample_indices
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.core.tool_helpers import is_stat_pca, stat_pca_rows
from metabolomix.dataset.artifacts import artifact_path
from metabolomix.plots import eic as eic_plot
from metabolomix.plots import pca_loadings as loadings_plot
from metabolomix.plots import group_intensity as group_intensity_plot
from metabolomix.plots import pca_scores as scores_plot
from metabolomix.plots import render as plot_render
from metabolomix.plots import species as species_plot
from metabolomix.plots import volcano as volcano_plot
from metabolomix.plots.result_output import build_annotated, figure_annotations
from metabolomix.tools import plot_species_kinds
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_plot"]

#: 結果（`ds.last_*` / `ds.results`）から描く kind。`result_id` を受ける。
RESULT_KINDS = ("volcano", "pca", "loadings")

#: 項目を選んで群ごとに描く kind（`plot_species_kinds`）が共有する引数。
_ITEM_ARGS = ("items", "groups", "low_reliability_samples", "standard_samples", "apply_curation",
              "exclude_auto_likely_wrong")

#: kind ごとの引数。ここに無い引数が渡されたら、どの kind の引数かを添えて拒む。
_KIND_ARGS: dict[str, tuple[str, ...]] = {
    "volcano": ("contrast_id",),
    "pca": ("color_by",),
    "loadings": ("pcs", "top_n", "value"),
    "group_intensity": _ITEM_ARGS + ("detection_limit",),
    "species": _ITEM_ARGS + ("value", "share_basis", "require_msms"),
    "species_pca": _ITEM_ARGS + ("require_msms", "normalize", "log_transform", "orient_by"),
    "eic": ("feature_id", "file_ids"),
    "eic_compounds": ("file_id", "names", "ontologies"),
}

#: `top_n=0`（loadings）は全特徴量。既定は主成分ごとに正負の上位 15。
_DEFAULT_TOP_N = 15


class _Missing(Exception):
    """前提の状態が無い（`missing_state` の封筒にする）。"""

    def __init__(self, state: str, required_tools: list[str], message: str):
        super().__init__(message)
        self.state, self.required_tools, self.message = state, required_tools, message


def _bad(message: str, **details) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details)


@dataclass
class _Request:
    result_id: str | None = None
    contrast_id: str | None = None
    title: str | None = None
    color_by: str | None = None
    pcs: list[int] | None = None
    top_n: int | None = None
    value: str | None = None
    feature_id: str | None = None
    file_ids: list[int] | None = None
    file_id: int | None = None
    names: list[str] | None = None
    ontologies: list[str] | None = None
    items: list[str] | None = None
    groups: list[str] | None = None
    low_reliability_samples: list[str] | None = None
    standard_samples: list[str] | None = None
    apply_curation: bool | None = None
    exclude_auto_likely_wrong: bool | None = None
    detection_limit: float | None = None
    share_basis: list[str] | None = None
    require_msms: bool | None = None
    normalize: str | None = None
    log_transform: bool | None = None
    orient_by: str | None = None
    given: set = field(default_factory=set)


def _record(kind: str, payload: dict, *, title, caption: str, result_id=None, result=None,
            metadata_revision: int | None = None) -> dict:
    """`ds.figures[kind]` に置く図の記録（`save_figure` が読む）。

    `metadata_revision` は、試料の情報（群・バッチ・役割）を色や題に焼き込んだ図だけが持つ。
    情報が更新されると（結果は現在のままでも）図は古い——`save_figure` が版の違いで拒む。
    """
    return {"kind": kind, "payload": payload, "title": title, "caption": caption,
            "result_id": result_id, "result": result, "metadata_revision": metadata_revision}


# ---------- 結果の選び方（古い結果は描かない） ----------

def _pick_result(ds, family: str, result_id: str | None, *, also_kinds: tuple = ()) -> dict:
    """PCA（family="pca"）か差次的結果（"differential"）を 1 つに決める。

    `result_id` があれば `ds.results` から引き、現在のものでなければ `STALE_ANALYSIS_RESULT`
    （古い結果を現在のものとして描かない）。無ければ直近（`ds.last_pca` / `ds.last_differential`）で、
    無ければ作り直すツールを案内する。`also_kinds` は `result_id` で名指しされたときに family のほかに
    受ける来歴の kind（loadings が `species_pca` の結果を受ける）。
    """
    tool = "dataset_pca" if family == "pca" else "dataset_compare"
    if result_id is None:
        last = ds.last_pca if family == "pca" else ds.last_differential
        if not last:
            what = "PCA" if family == "pca" else "2 群の差次的解析"
            raise _Missing("pca_result" if family == "pca" else "dataset_compare_result", [tool],
                           f"直近の{what}の結果がありません。先に {tool} を実行してください"
                           "（前処理・除外・群の変更で古くなった結果は捨てられます）。")
        result_state.assert_current(ds, last)
        return last
    result = (ds.results or {}).get(result_id)
    if result is None:
        raise _bad(f"result_id={result_id} の結果がありません。{tool} の戻り値の result_id を渡すか、"
                   "省いて直近の結果を使ってください。", result_id=result_id)
    accepted = (*also_kinds, "stat_pca" if family == "pca" else "compare")
    if (result.get("provenance") or {}).get("kind") not in accepted:
        raise _bad(f"result_id={result_id} は{'PCA' if family == 'pca' else '差次的解析'}の結果ではありません"
                   f"（{tool} の result_id を渡してください）。", result_id=result_id)
    result_state.assert_current(ds, result)
    return result


def _result_id(result: dict) -> str | None:
    return (result.get("provenance") or {}).get("result_id")


# ---------- volcano / pca / loadings ----------

def _volcano(ds, req: _Request) -> dict:
    result_id, contrast_id = req.result_id, req.contrast_id
    if result_id is not None:                       # `result_id#contrast_id`（特徴検索が出す形）も受ける
        result_id, contrast_id = stats_results.split_result_ref(result_id, contrast_id)
    record = _pick_result(ds, "differential", result_id)
    # 対比 1 つを 2 群の形に読む（複数の対比で contrast_id が無ければ CONTRAST_REQUIRED）。
    result = stats_results.as_two_group(record, contrast_id)
    if not result.get("volcano"):
        raise _bad("この差次的結果は 2 群比較の volcano 点列を持っていません（volcano は 2 群比較だけ）。",
                   result_id=_result_id(result))
    payload = volcano_plot.build_volcano_plot_payload(result, title=req.title)
    points = [p for p in result["volcano"] if volcano_plot._is_drawable(p)]
    up = sum(1 for p in points if p.get("sig") == "up")
    down = sum(1 for p in points if p.get("sig") == "down")
    caption = (f"Volcano: {result.get('a')} (n={result.get('n_a')}) vs {result.get('b')} (n={result.get('n_b')})"
               f" — 全 {len(result['volcano'])} 特徴のうち up={up} / down={down} / ns={len(points) - up - down}"
               f"（検定不能 {len(result['volcano']) - len(points)} 件は非描画）。"
               f" しきい値 q<{result.get('q_threshold')}, |log2FC|>={result.get('log2fc_threshold')}。"
               " 個々の特徴量名と統計値は dataset_compare の結果を参照。")
    return _record("volcano", payload, title=req.title, caption=caption,
                   result_id=_result_id(result), result=result)


def _explained(result: dict) -> list:
    """PCA の記録の寄与率（統計エンジンの結果は `engine` の中、分子種 PCA は記録の直下）。"""
    source = result.get("engine") if is_stat_pca(result) else result
    return list((source or {}).get("explained_variance_ratio") or [])


def _pca(ds, req: _Request) -> dict:
    color_by = req.color_by or "group"
    if color_by not in scores_plot.COLOR_BY:
        raise _bad(f"color_by は {' / '.join(scores_plot.COLOR_BY)} のどれかを指定してください"
                   f"（受け取った値: {color_by!r}）。", color_by=color_by)
    result = _pick_result(ds, "pca", req.result_id)
    caveats = []
    scores = stat_pca_rows(result)
    projected = [row["name"] for row in scores if not row["fitted"]]
    if color_by == "group":
        colors = {row["name"]: row["group"] for row in scores}
    elif color_by == "role":
        colors = {row["name"]: row["role"] for row in scores}
    else:
        entry = stats_results.get_matrix(ds, result.get("matrix_id"))
        batch = dict(stats_results.sample_batches(ds))           # QC の注入（試料 ID）
        batch.update({u.unit_id: u.batch for u in entry["prepared"].units} if entry else {})
        colors = {row["name"]: batch.get(row["id"]) for row in scores}
    if projected:
        caveats.append("QC は主成分の軸の計算に使わず、同じ軸へ投影した点（白抜き）。")
    try:
        payload = scores_plot.build_pca_scores_payload(
            scores, _explained(result), colors, color_by=color_by,
            result_id=_result_id(result), title=req.title, caveats=caveats, projected=projected)
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    evr = _explained(result)
    caption = (f"PCA スコア（{len(payload['points'])} 試料、色={color_by}、result_id={_result_id(result)}）: "
               + "、".join(f"PC{i + 1} {v * 100:.1f}%" for i, v in enumerate(evr[:2])) + "。")
    return _record("pca", payload, title=req.title, caption=caption, result_id=_result_id(result),
                   metadata_revision=getattr(ds, "metadata_revision", 0))


def _loadings(ds, req: _Request) -> dict:
    pcs = list(req.pcs) if req.pcs is not None else [1, 2]
    value = req.value or "r"
    if req.top_n is not None and (isinstance(req.top_n, bool) or not isinstance(req.top_n, int)
                                  or req.top_n < 0):
        raise _bad(f"top_n は 0 以上の整数を指定してください（0 は全特徴量。受け取った値: {req.top_n!r}）。")
    top_n = _DEFAULT_TOP_N if req.top_n is None else (None if req.top_n == 0 else req.top_n)
    result = _pick_result(ds, "pca", req.result_id, also_kinds=("species_pca",))
    try:
        if (result.get("provenance") or {}).get("kind") == "species_pca":
            # 分子種 PCA: 特徴の行は結果が持つ（相関 r は autoscale した PCA で計算済み）
            features, r_ok = result["loadings_rows"], True
        else:
            features, r_ok = loadings_plot.dataset_features(ds, result)
        payload = loadings_plot.build_loadings_payload(
            features, source="dataset", result_id=_result_id(result),
            explained_variance_ratio=_explained(result), pcs=pcs,
            top_n=top_n, value=value, r_available=r_ok)
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    parts = "、".join(f"{p['pc']} {p['explained_pct']}%（{len(p['rows'])} 行）" for p in payload["panels"])
    caption = (f"PCA ローディング（result_id={payload['result_id'] or '(なし)'}、値={payload['value']}）: {parts}。")
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return _record("loadings", payload, title=req.title, caption=caption, result_id=_result_id(result))


# ---------- eic / eic_compounds ----------

def _eic_file(ds):
    path = artifact_path(ds, "eic")
    if path is None:
        raise DomainError(
            "COMPANION_ARTIFACT_MISSING",
            "このデータセットには .EIC.aef が無いため、EIC は描けません（mzTab-M だけから読んだデータセット、"
            "または .EIC.aef を持たない解析フォルダ）。.EIC.aef を含む MS-DIAL の解析フォルダを "
            "load_dataset(directory=...) で開き直してください。",
            {"role": "eic", "artifact_paths": sorted(getattr(ds, "artifact_paths", None) or {})})
    return path


def _catalog_records(ds) -> list[dict]:
    """データセットの特徴表を、EIC の照合部品（`eic/identity_map` `eic/mz_backfill`）が食う
    `.arf2` 風のレコードにする。特徴 ID は MasterAlignmentID（整数）のものだけ。"""
    heights = np.zeros(len(ds.feature_ids))
    matrix = getattr(ds, "feature_matrix", None)
    if matrix is not None and len(ds.feature_ids):
        values = np.asarray(matrix, dtype=float)
        kept = included_sample_indices(ds)
        if kept is not None:
            values = values[:, kept]
        finite = np.isfinite(values)
        heights = np.where(finite, values, 0.0).sum(axis=1) / np.maximum(finite.sum(axis=1), 1)
    records = []
    for index, fid in enumerate(ds.feature_ids):
        try:
            alignment_id = int(fid)
        except (TypeError, ValueError):
            continue
        meta = (ds.feature_metadata or {}).get(fid) or {}
        records.append({"AlignmentID": alignment_id, "Name": meta.get("name"),
                        "Ontology": meta.get("ontology"), "AdductType": meta.get("adduct"),
                        "RT": meta.get("rt"), "MassCenter": meta.get("mz"),
                        "HeightAverage": float(heights[index])})
    return records


def _backfill_notes(ds, spots: list[dict], records: list[dict]) -> list[str]:
    """欠損 m/z（古い MS-DIAL の -1）を特徴表の m/z で補い、注記を返す（`eic/mz_backfill`）。"""
    from metabolomix.eic import mz_backfill

    if not mz_backfill.has_missing_mz(spots):
        return []
    arf2 = artifact_path(ds, "arf2")
    note = mz_backfill.backfill_caveat(mz_backfill.backfill_missing_mz(spots, records),
                                       arf2.name if arf2 else "データセットの特徴表")
    return [note] if note else []


def _contiguity_note(ds) -> list[str]:
    ids = []
    for fid in ds.feature_ids:
        try:
            ids.append(int(fid))
        except (TypeError, ValueError):
            return []
    if sorted(ids) == list(range(len(ids))):
        return []
    return ["特徴 ID は .arf2 の MasterAlignmentID で、連番になっていない（イオンモビリティのデータなど）。"
            ".EIC.aef のスポット番号も同じ MasterAlignmentID の空間なので、特徴 ID をそのまま番号として引いた"
            "（RT と m/z を特徴表と照合している）。"]


def _default_file_ids(ds) -> list[int] | None:
    """`file_ids` を省いたときの試料。解析に残さない試料（`included_sample_indices`）があるときだけ、
    残す試料の FileID に絞る（全試料ならリーダーの既定＝すべて）。FileID を決められなければ None。"""
    kept = included_sample_indices(ds)
    if kept is None:
        return None
    assays = list(getattr(ds, "sample_assay_ids", None) or [])
    if len(assays) != len(ds.sample_names) or not all(str(a).isdigit() for a in assays):
        return None
    return [int(assays[i]) for i in kept] or None


def _eic(ds, req: _Request) -> dict:
    from metabolomix.eic.identity_map import verify_spot_match
    from metabolomix.eic.reader import read_eic_spot_css1

    path = _eic_file(ds)
    if req.feature_id is None:
        raise _bad("kind=eic には feature_id が要ります（dataset_feature_search の feature_id）。")
    fid = str(req.feature_id)
    if fid not in set(ds.feature_ids):
        raise _bad(f"feature_id={fid!r} はこのデータセットの特徴にありません"
                   "（dataset_feature_search で探せます）。", feature_id=fid)
    try:
        spot_id = int(fid)
    except ValueError:
        raise _bad(f"feature_id={fid!r} は整数ではなく、.EIC.aef のスポット番号（MasterAlignmentID）に"
                   "対応づけられません。", feature_id=fid) from None
    file_ids = req.file_ids if req.file_ids is not None else _default_file_ids(ds)
    try:
        spot = read_eic_spot_css1(str(path), spot_id, file_ids=file_ids)
    except ValueError as exc:
        raise _bad(f"EIC を読めませんでした: {exc}（file_ids で 12 試料以内に絞ってください）。",
                   feature_id=fid) from exc
    notes = _backfill_notes(ds, [spot], _catalog_records(ds)) + _contiguity_note(ds)
    meta = (ds.feature_metadata or {}).get(fid) or {}
    if meta.get("rt") is not None and meta.get("mz") is not None and float(spot["mz"]) > 0:
        reason = verify_spot_match({"rt": meta["rt"], "mz": meta["mz"]}, spot)
        if reason:
            notes.append(f"feature {fid} の RT/m/z（{float(meta['rt']):.3f} min / {float(meta['mz']):.4f}）と "
                         f"EIC スポット {spot_id} が合わない（{reason}）。番号の対応が崩れている可能性があり、"
                         "別の特徴のクロマトグラムを見ているかもしれません。")
    payload = eic_plot.build_eic_plot_payload(spot, path, title=req.title)
    payload["caveats"].extend(notes)
    name = meta.get("name") or fid
    caption = (f"EIC（feature {fid} {name}、m/z {spot['mz']:.4f}、RT {spot['rt']:.2f} min）: "
               f"{len(payload['series'])} 試料のクロマトグラム。")
    return _record("eic", payload, title=req.title, caption=caption)


def _eic_compounds(ds, req: _Request) -> dict:
    from metabolomix.eic.identity_map import select_identity_candidates
    from metabolomix.eic.reader import read_eic_spots_css1

    path = _eic_file(ds)
    file_id = req.file_id
    if isinstance(file_id, bool) or not isinstance(file_id, int):
        raise _bad("kind=eic_compounds には整数の file_id が要ります（1 回の呼び出しで 1 試料。"
                   "試料を比べるときは試料ごとに呼びます）。")
    records = _catalog_records(ds)
    try:
        candidates, caveats, total_matched = select_identity_candidates(
            records, names=req.names, ontologies=req.ontologies)
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    if not candidates:
        raise _bad(f"names={req.names} / ontologies={req.ontologies} に一致する特徴がありません"
                   "（dataset_feature_search で名前・クラスの表記を確認してください）。")
    try:
        spots = read_eic_spots_css1(str(path), [c["spot_id"] for c in candidates], file_ids=[file_id])
    except ValueError as exc:
        raise _bad(f"EIC を読めませんでした: {exc}") from exc
    caveats = caveats + _backfill_notes(ds, spots, records) + _contiguity_note(ds)
    arf2 = artifact_path(ds, "arf2")
    try:
        payload = eic_plot.build_multi_compound_plot_payload(
            candidates, spots, file_id=file_id, file_path=path, arf2_path=arf2,
            top_n=eic_plot.DEFAULT_COMPOUND_TOP_N, title=req.title, queries=req.names,
            ontologies=req.ontologies, caveats=caveats, total_matched=total_matched)
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    selection = payload["selection"]
    drawn = "、".join(item["label"] for item in payload["series"]) or "なし"
    sample = payload["sample"].get("sample_name") or f"FileID {file_id}"
    caption = (f"EIC 重ね描き（{sample}）: 候補 {selection['candidates']} 件のうち "
               f"{selection['plotted']} 件を描画 — {drawn}。")
    if selection["dropped"]:
        caption += (f" 除外 {len(selection['dropped'])} 件（理由は output='payload' の "
                    "selection.dropped、または caveats を参照）。")
    return _record("eic_compounds", payload, title=req.title, caption=caption)


#: kind -> 図の記録を作る関数 `(ds, _Request) -> record`。
_KINDS: dict[str, Callable[[object, _Request], dict]] = {
    "volcano": _volcano,
    "pca": _pca,
    "loadings": _loadings,
    "group_intensity": plot_species_kinds.group_intensity,
    "species": plot_species_kinds.species,
    "species_pca": plot_species_kinds.species_pca,
    "eic": _eic,
    "eic_compounds": _eic_compounds,
}


def build_figure(ds, kind: str, **args) -> dict:
    """`kind` の図の記録を作る（`ds.figures` へは置かない。`save_figure(result_id=...)` も使う）。"""
    request = _Request(**args)
    return _KINDS[kind](ds, request)


# ---------- 描画 ----------

def _png(ds, record: dict) -> bytes:
    kind, title = record["kind"], record["title"]
    if kind == "volcano":
        notes = figure_annotations(ds, record["result"])
        return plot_render.render_png(build_annotated, volcano_plot.render_volcano_plot, notes,
                                      record["result"], title=title)
    if kind == "pca":
        notes = figure_annotations(ds, None)
        return plot_render.render_png(build_annotated, scores_plot.render_pca_scores, notes,
                                      record["payload"], title=title)
    if kind == "loadings":
        return plot_render.render_png(loadings_plot.render_loadings_plot, record["payload"], title=title)
    if kind == "group_intensity":
        return plot_render.render_png(group_intensity_plot.render_group_intensity_plot,
                                      record["payload"], title=title)
    if kind == "species":
        return plot_render.render_png(species_plot.render_species_plot, record["payload"], title=title)
    if kind == "species_pca":
        return plot_render.render_png(scores_plot.render_pca_scores, record["payload"],
                                      title=title or record["payload"]["title"])
    return plot_render.render_png(eic_plot.render_eic_plot, record["payload"], title=title)


def _check_arguments(kind: str, request: _Request) -> None:
    allowed = set(_KIND_ARGS[kind]) | ({"result_id"} if kind in RESULT_KINDS else set())
    for name in sorted(request.given - allowed):
        owners = ([k for k, args in _KIND_ARGS.items() if name in args]
                  or (list(RESULT_KINDS) if name == "result_id" else []))
        raise _bad(f"{name} は kind={' / '.join(owners)} の引数です（kind={kind!r} では使いません）。",
                   argument=name, kind=kind, belongs_to=owners)


@mcp.tool(annotations=ToolAnnotations(
    # 図を ds.figures へ、species_pca は結果を ds.results へ登録する（メモリ上。ファイルは save_figure）。
    readOnlyHint=False, destructiveHint=False, idempotentHint=False),
    structured_output=False)
def dataset_plot(
    kind: str,
    dataset_handle: str | None = None,
    result_id: str | None = None,
    output: str | None = None,
    title: str | None = None,
    color_by: str | None = None,
    pcs: list[int] | None = None,
    top_n: int | None = None,
    value: str | None = None,
    feature_id: str | None = None,
    file_ids: list[int] | None = None,
    file_id: int | None = None,
    names: list[str] | None = None,
    ontologies: list[str] | None = None,
    items: list[str] | None = None,
    groups: list[str] | None = None,
    low_reliability_samples: list[str] | None = None,
    standard_samples: list[str] | None = None,
    apply_curation: bool | None = None,
    exclude_auto_likely_wrong: bool | None = None,
    detection_limit: float | None = None,
    share_basis: list[str] | None = None,
    require_msms: bool | None = None,
    normalize: str | None = None,
    log_transform: bool | None = None,
    orient_by: str | None = None,
    contrast_id: str | None = None,
) -> list | str:
    """データセットの図を描く — 既定は PNG 画像と 1 行の説明。ファイルは書かない（保存は `save_figure`）。

    kind と、その kind だけの引数（別の kind の引数を渡すと DATASET_BAD_REQUEST）:
      - `volcano`: 直近の 2 群の差次的結果（`dataset_compare`）。`result_id` で名指しもできる。対比が複数の結果は
        `contrast_id`（`T_vs_C` など）で 1 つ選ぶ（省くと `CONTRAST_REQUIRED`、交互作用は `CONTRAST_NOT_TWO_GROUP`）。
      - `pca`: 直近の `dataset_pca`（`result_id` で名指し）。`color_by="group"`（既定）/ `"batch"` / `"role"`。
      - `loadings`: 同じ PCA のローディング（`result_id` に `species_pca` の結果も渡せる。横棒、色は特徴のクラス＝`ontology`）。`pcs`（既定 [1, 2]、
        1〜3 個）、`top_n`（主成分ごとに正負の上位 N。既定 15、0 は全特徴量で 60 まで）、
        `value="r"`（既定。特徴量と主成分スコアの相関）/ `"coefficient"`。
      - `eic`: `feature_id`（`dataset_feature_search` の ID）の試料別クロマトグラム。`file_ids` で試料を絞る
        （省けば全試料。12 を超えるときは要指定。除外した試料は省いた呼び出しから外れる）。
      - `eic_compounds`: 1 試料（`file_id`）に、同定済みの物質（`names` は名前の部分一致、`ontologies` は
        クラスの完全一致。どちらかが要る）を重ねる。RT と m/z が特徴表と合わない物質は描かず、理由は
        payload の `selection.dropped`。
      - `group_intensity`: 項目（`items` = クラス名〔`ontology` に完全一致。配下の全特徴の合計〕・分子種名〔`|` で
        分けた候補名に完全一致〕・特徴 ID。`"PE+EtherPE"` のように `+` で合算。当たらない項目は N.D.、1〜30 件）
        ごとに、群（`groups` = 試料名・Class ID・因子のトークン。1 群以上。QC・ブランクは `qc` / `blank` を
        書いたときだけ入る）の試料別の強度合計を 1 パネルに並べる。縦軸は log10 で全パネル共通、群は log10 の
        平均 ± SD。検定はしない（`dataset_compare`）。`detection_limit` は検出下限の破線（省略時は
        param ファイルの Minimum peak height）。
      - `species`: `items` を分子種（特徴）に展開し、1 分子種 1 パネル（1 回 40 まで）。`value="share"`（既定。
        試料ごとの 値 ÷ 分母の合計 × 100。棒 = 群平均、ひげ = SD）/ `"intensity"`（log10 の平均 ± SD。payload の `value` は形式を保つため `"height"`）。
        `share_basis` は割合の分母にする項目（省略時は `items`。当たらない部品があると分母が変わるので描かずに
        エラー）。`require_msms=True` で MS/MS の裏付けの無い分子種も（分母からも）除く。
      - `species_pca`: 選んだ分子種（`species` と同じ `items` `groups` `require_msms`）の PCA のスコア図。
        `normalize="none"`（既定）/ `"total"`、`log_transform`（既定 True）、`orient_by`（群名。その群の平均
        スコアが正になる向きにそろえる）。`low_reliability_samples` は計算に使わず投影だけ。結果は
        `result_id`（`lipidmix.species_pca.v1`）で登録する。
      3 つの共通: `low_reliability_samples`（白抜き。平均・SD から外す）、`standard_samples`（その試料にだけ
        実質的な強度がある分子種をクラスの合計から外す）、`apply_curation`（既定 True。`.arf2` の人の判断
        〔wrong / redundant 除外・assign の付け替え〕）、`exclude_auto_likely_wrong`（最新の `curation_review` の
        likely_wrong も除く）。除いたものは payload の `excluded`。内部標準（`(d7)` 等）はクラスの合計に入れない。
        `dataset_exclude` した試料は描かず、除外した特徴はクラスにも名前指定にも入れない。
        群・試料の情報は図に焼き込む（試料の情報を更新したあとの図は `save_figure` が古いとして拒む）。
    `.EIC.aef` の無いデータセット（mzTab-M だけ）の `eic` 系は COMPANION_ARTIFACT_MISSING。

    結果は**現在のものだけ**描く: 前処理・除外・群の変更で古くなった結果の `result_id` は
    STALE_ANALYSIS_RESULT で拒む。`result_id` を省いて直近の結果が無ければ missing_state（作り直すツールを案内）。

    output: `"image"`（既定）= PNG と説明。`"payload"` = 描画中立の点列 JSON（形式名 `lipidmix.volcano.v1`
      `lipidmix.pca_scores.v1` `lipidmix.pca_loadings.v1` `lipidmix.group_intensity.v1`
      `lipidmix.species_intensity.v1` `lipidmix.species_pca.v1` `lipidmix.eic.v1` `lipidmix.eic.multi.v1`。
      volcano は有意な点を全件、ns の点を間引いて返す）。環境変数 `LIPIDMIX_PLOT_OUTPUT` で既定を変えられる。
      座標を LLM に渡すより画像のほうが桁違いに安いので、画面に出すだけなら image。
    title: 図の題。dataset_handle: `load_dataset` が返したハンドル（省略時は直近）。

    作った図の記録は `ds.figures[kind]` に置き（kind ごとに最後の 1 つ）、`save_figure(kind=...)` が描く。
    """
    request = _Request(
        result_id=result_id, contrast_id=contrast_id, title=title, color_by=color_by, pcs=pcs, top_n=top_n, value=value,
        feature_id=feature_id, file_ids=file_ids, file_id=file_id, names=names, ontologies=ontologies,
        items=items, groups=groups, low_reliability_samples=low_reliability_samples,
        standard_samples=standard_samples, apply_curation=apply_curation,
        exclude_auto_likely_wrong=exclude_auto_likely_wrong, detection_limit=detection_limit,
        share_basis=share_basis, require_msms=require_msms, normalize=normalize,
        log_transform=log_transform, orient_by=orient_by)
    request.given = {name for name, given in (
        ("result_id", result_id), ("contrast_id", contrast_id), ("color_by", color_by), ("pcs", pcs), ("top_n", top_n),
        ("value", value), ("feature_id", feature_id), ("file_ids", file_ids),
        ("file_id", file_id), ("names", names), ("ontologies", ontologies), ("items", items),
        ("groups", groups), ("low_reliability_samples", low_reliability_samples),
        ("standard_samples", standard_samples), ("apply_curation", apply_curation),
        ("exclude_auto_likely_wrong", exclude_auto_likely_wrong), ("detection_limit", detection_limit),
        ("share_basis", share_basis), ("require_msms", require_msms), ("normalize", normalize),
        ("log_transform", log_transform), ("orient_by", orient_by)) if given is not None}
    if kind not in _KINDS:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"kind は {', '.join(_KINDS)} のどれかを指定してください（受け取った値: {kind!r}）。",
                           {"kinds": list(_KINDS)})
    try:
        _check_arguments(kind, request)
        mode = plot_render.resolve_plot_output(output)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except ValueError as exc:
        return mztab_error("DATASET_BAD_REQUEST", str(exc), {"output": output})
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    try:
        record = _KINDS[kind](ds, request)
        if mode == plot_render.PAYLOAD:
            body = json_payload(record["payload"])
        else:
            png = _png(ds, record)
            body = [record["caption"], Image(data=png, format="png")]
    except _Missing as exc:
        return missing_state(exc.state, exc.required_tools, exc.message)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except ValueError as exc:                      # 描くものが無い等（描画部品の検査）
        return mztab_error("DATASET_BAD_REQUEST", str(exc), {"kind": kind})
    ds.figures[kind] = record
    return body
