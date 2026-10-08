"""ARF の項目指定（クラス・分子種名）を、描画・PCA が使うスポットの集合に解決する。

群別強度（arf_plot_group_intensity）・分子種ごとの図（arf_plot_species）・分子種 PCA
（arf_pca_species）が同じ規則で分子種を選ぶための共通部品。項目の文字列の解決そのものは
`plots/group_intensity.py resolve_items()`。ここは ARF の行・`.arf2`・キュレーションの判断と
最新レビュー・標準液の判定を集めて渡す層。MCP には依存しない（session の ArfState は引数で受ける）。
spec: docs/superpowers/specs/2026-10-08-species-plot-and-pca-loadings-design.md §3。
"""
from __future__ import annotations

from dataclasses import dataclass, field

from metabolomix.arf import reader as arf_reader
from metabolomix.arf2 import reader as arf2_reader
from metabolomix.arf2.match_results import load_spot_annotations, name_prefix
from metabolomix.curation import apply as curation_apply
from metabolomix.curation import flags as curation_flags
from metabolomix.curation import suggest as curation_suggest
from metabolomix.msdial import sample_factors
from metabolomix.plots import group_intensity as gi


@dataclass
class Selection:
    items: list
    excluded: dict
    groups: list
    low_reliability: frozenset
    rows_by_spot: dict
    msms: dict
    catalog_by_id: dict
    caveats: list = field(default_factory=list)


def display_name(name: str) -> str:
    """`PG 34:1|PG 16:0_18:1` → `PG 16:0_18:1`（鎖まで決まった最後の候補）。接頭辞は残す。"""
    return (name or "").split("|")[-1].strip()


def build_selection(arf_state, arf2_path, *, items, groups, low_reliability_samples=None,
                    apply_curation=True, exclude_auto_likely_wrong=False, standard_samples=None):
    """項目・群・除外を解決する。入力の誤りは ValueError、記録ファイルの破損は FlagFileError。"""
    excluded_samples = set(arf_state.excluded_samples or ())
    rows_by_spot = {}
    all_rows_by_spot = {}    # 標準液の判定用: arf_exclude 済みの試料の行も含める
    for feature in arf_state.features:
        rows = [r for r in (arf_reader.alignment_feature_row(raw) for raw in feature["AlignedPeakProperties"]) if r]
        all_rows_by_spot[feature["MasterAlignmentID"]] = rows
        rows_by_spot[feature["MasterAlignmentID"]] = [
            r for r in rows if r.get("file_name") not in excluded_samples]
    names = [n for n in sample_factors.arf_sample_names(arf_state.features) if n not in excluded_samples]
    facets = sample_factors.build_sample_facets(names, arf_state.class_index)
    resolved_groups, caveats = gi.resolve_groups(groups, facets)
    plotted = {n for g in resolved_groups for n in g["samples"]}
    low = gi.resolve_sample_specs(low_reliability_samples or [], facets)

    catalog = arf2_reader.load_catalog(str(arf2_path))
    curation = {}
    if apply_curation:
        flag_set = curation_apply.flags_for_arf2(arf2_path)
        curation = {
            "wrong": {int(s) for s in flag_set["wrong"]},
            "redundant": {int(s) for s in flag_set["redundant"]},
            "assign": {int(s): {"name": row.get("name") or "", "ontology": row.get("ontology") or ""}
                       for s, row in flag_set["assign"].items()},
        }
        if flag_set["orphaned"]:
            caveats.append(curation_flags.orphaned_warning(flag_set["orphaned"]))
    if exclude_auto_likely_wrong:
        review = curation_suggest.latest_review(
            arf2_path, curation_flags.alignment_key(arf2_path)["alignment_sha256"])
        if review is None:
            caveats.append("このアラインメントの curation_review のレビューが無いため、自動判定 likely_wrong は除いていません。")
        else:
            curation["auto_likely_wrong"] = {
                int(s["spot_id"]) for s in review["spots"] if s.get("verdict") == "likely_wrong"}

    standard = set()
    if standard_samples:
        standard = gi.resolve_sample_specs(standard_samples, sample_factors.build_sample_facets(
            sample_factors.arf_sample_names(arf_state.features), arf_state.class_index))
    standard_only = gi.standard_only_spots(all_rows_by_spot, standard, plotted) if standard else frozenset()
    manual = {int(s) for s in (arf_state.excluded_spots or ())}    # arf_exclude で除いたスポット
    resolved_items, excluded = gi.resolve_items(
        items, catalog, curation=curation, standard_only=standard_only, manual=manual)

    # MS/MS の裏付け: 照合結果に MS/MS があり、かつ .arf2 の Name が `no MS2:` / `w/o MS2:` でない
    # （接頭辞は照合結果ではなく .arf2 の Name に付く）。assign 済みでも元のスペクトルの有無で判断する。
    annotations = load_spot_annotations(str(arf2_path))
    name_by_id = {int(r["MasterAlignmentID"]): r.get("Name") or "" for r in catalog}
    msms = {}
    for sid, ann in annotations.items():
        rep = ann.get("representative") or {}
        msms[sid] = (bool(rep.get("has_msms"))
                     and name_prefix(name_by_id.get(sid, "")) not in ("no MS2", "w/o MS2"))
    return Selection(items=resolved_items, excluded=excluded, groups=resolved_groups,
                     low_reliability=frozenset(low), rows_by_spot=rows_by_spot, msms=msms,
                     catalog_by_id={int(r["MasterAlignmentID"]): r for r in catalog}, caveats=caveats)


def expand_spots(selection: Selection, *, require_msms: bool = False):
    """解決済みの項目をスポット単位に展開する（項目の順 → 表示名の順。同じスポットは 1 回）。"""
    seen = set()
    spots = []
    dropped = []
    for item in selection.items:
        ordered = sorted(item["spots"], key=lambda s: (display_name(s["name"]).casefold(), s["spot_id"]))
        for spot in ordered:
            sid = int(spot["spot_id"])
            if sid in seen:
                continue
            seen.add(sid)
            has_msms = bool(selection.msms.get(sid))
            if require_msms and not has_msms:
                dropped.append(f"#{sid} {spot['name']}")
                continue
            row = selection.catalog_by_id.get(sid, {})
            spots.append({
                "spot_id": sid, "name": spot["name"], "label": display_name(spot["name"]),
                "ontology": spot["ontology"], "item": item["item"],
                "adduct": row.get("AdductType") or "",
                "mz": round(float(row.get("MassCenter") or 0.0), 4),
                "rt": round(float(row.get("RT") or 0.0), 3),
                "msms": has_msms,
            })
    return spots, sorted(dropped)
