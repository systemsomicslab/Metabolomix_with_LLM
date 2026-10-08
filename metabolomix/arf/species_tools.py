"""ARF の分子種ごとの図（arf_plot_species）と分子種 PCA（arf_pca_species）の MCP ツール。

項目の解決は arf/selection.py、描画は plots/species.py / plots/pca_scores.py。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §4・§5。
"""
from __future__ import annotations

from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations

from metabolomix.arf import selection as arf_selection
from metabolomix.arf.tools import _check_ncols, _require_arf_with_arf2
from metabolomix.core import session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.curation import flags as curation_flags
from metabolomix.plots import render as plot_render
from metabolomix.plots import species as species_plot

__all__ = ["arf_plot_species"]

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
