"""`DatasetState` から、群別強度・分子種の図・分子種 PCA が使う行（特徴 × 試料）を選ぶ（MCP 非依存）。

旧 ARF 経路（`arf/selection.py`）が ARF の行と `.arf2` から作っていた「項目 → スポット」「スポット × 試料の行」を、
データセットの正準の欄から作り直したもの。項目の文字列の解決そのものは `plots/group_intensity.py resolve_items()`
（クラス名・分子種名・`+` の合算）を使い、ここは次の欄を渡す層になる:

    値            `feature_matrix`（特徴 × 試料。NaN は行が無いのと同じ）
    gap-fill     `detected_mask`（False が gap-filled。無ければ印を付けない）
    群・役割      試料の対応表（`sample_meta` の元の行。`_metadata_rows`）の `group` / `role`
    除外          試料は `included_sample_indices`（include=false ∪ `dataset_exclude`）、特徴は `excluded_features`
    キュレーション `ds.curation_flags`（`load_dataset` が `.arf2` から読んだもの）と、最新レビューの likely_wrong
    MS/MS        `feature_metadata[fid]["msms_assigned"]`（`.arf` のどれかのセルが MS/MS を持つか）と、
                 `.arf2` の Name の `no MS2:` / `w/o MS2:` の接頭辞
    標準液        `standard_samples` にだけ実質的な強度があるスポット（`gi.STANDARD_ONLY_FOLD`）

特徴 ID は整数（`.arf2` の MasterAlignmentID、mzTab-M の SMF_ID）のものだけを扱う（整数でないものは数えて注記する）。
旧ツールの引数・除外の意味を保つ。入力の誤りは `DomainError("DATASET_BAD_REQUEST")`。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from metabolomix.analysis import sample_manifest
from metabolomix.analysis.dataset_analysis import included_sample_indices
from metabolomix.arf2.match_results import name_prefix
from metabolomix.core.atomic_io import DomainError
from metabolomix.curation import flags as curation_flags
from metabolomix.curation import suggest as curation_suggest
from metabolomix.curation.apply import confirmed_spots
from metabolomix.dataset.artifacts import artifact_path
from metabolomix.msdial import sample_factors
from metabolomix.plots import group_intensity as gi

__all__ = ["Resolved", "SampleSelection", "display_name", "expand_spots", "resolve_samples",
           "resolve_selection", "rows_by_spot", "select_items", "standard_only_spots", "unmatched_parts"]


def _bad(message: str, **details) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details)


def display_name(name: str) -> str:
    """`PG 34:1|PG 16:0_18:1` → `PG 16:0_18:1`（鎖まで決まった最後の候補）。接頭辞は残す。"""
    return (name or "").split("|")[-1].strip()


# ---------- 試料（群・低信頼・標準液） ----------

@dataclass
class SampleSelection:
    groups: list                       # [{"label", "samples"}]（指定順。除外した試料は含まない）
    low_reliability: frozenset
    standard: frozenset                # 標準液の試料（除外した試料も含む。標準液だけにあるスポットの判定用）
    caveats: list = field(default_factory=list)


def _facets(ds, names: list[str]) -> dict:
    # 試料の対応表の行（適用済みの `sample_metadata_rows`。無ければ読み込みが決めた既定の行）
    table = getattr(ds, "sample_metadata_rows", None) or sample_manifest.resolve_metadata(ds, None)
    rows = dict(zip(ds.sample_names, table))
    return sample_factors.build_sample_facets(
        names, None, {name: {"role": (rows.get(name) or {}).get("role"),
                             "group": (rows.get(name) or {}).get("group")} for name in names})


def _included_names(ds) -> list[str]:
    keep = included_sample_indices(ds)
    names = list(ds.sample_names)
    return names if keep is None else [names[i] for i in keep]


def resolve_samples(ds, groups, low_reliability_samples=None, standard_samples=None) -> SampleSelection:
    """群・低信頼・標準液の指定を試料名にする。群は残す試料だけから解く（除外した試料は描かない）。"""
    included = _included_names(ds)
    try:
        resolved, caveats = gi.resolve_groups(groups, _facets(ds, included))
        low = gi.resolve_sample_specs(low_reliability_samples or [], _facets(ds, included))
        standard = (gi.resolve_sample_specs(standard_samples, _facets(ds, list(ds.sample_names)))
                    if standard_samples else set())
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    return SampleSelection(resolved, frozenset(low), frozenset(standard), list(caveats))


# ---------- 特徴の表 ----------

def _spot_ids(ds) -> tuple[dict[int, int], int]:
    """{整数の特徴 ID: ds.feature_ids の位置}と、整数でなくて外した特徴の数。"""
    index: dict[int, int] = {}
    skipped = 0
    for position, fid in enumerate(ds.feature_ids):
        try:
            index[int(fid)] = position
        except (TypeError, ValueError):
            skipped += 1
    return index, skipped


def _meta(ds, spot_id: int, index: dict[int, int]) -> dict:
    return (ds.feature_metadata or {}).get(ds.feature_ids[index[spot_id]]) or {}


def _catalog(ds, index: dict[int, int]) -> list[dict]:
    return [{"MasterAlignmentID": sid, "Name": _meta(ds, sid, index).get("name") or "",
             "Ontology": _meta(ds, sid, index).get("ontology") or ""} for sid in sorted(index)]


def _msms(ds, index: dict[int, int]) -> dict[int, bool] | None:
    """スポットごとの MS/MS の裏付け。証拠の欄（`msms_assigned`）を持たないデータセットは None。

    `.arf` のどれかのセルが MS/MS を持ち、かつ Name が `no MS2:` / `w/o MS2:` で始まらない
    （接頭辞は照合結果ではなく `.arf2` の Name に付く）。assign 済みでも元のスペクトルの有無で判断する。
    """
    if not any("msms_assigned" in (ds.feature_metadata or {}).get(fid, {}) for fid in ds.feature_ids):
        return None
    out = {}
    for sid in index:
        meta = _meta(ds, sid, index)
        out[sid] = bool(meta.get("msms_assigned")) and name_prefix(meta.get("name") or "") not in (
            "no MS2", "w/o MS2")
    return out


def _curation(ds, *, apply_curation: bool, exclude_auto_likely_wrong: bool, caveats: list) -> dict:
    """キュレーションの判断（`resolve_items` の `curation`）。人の判断は `ds.curation_flags`、
    自動判定は最新レビューの likely_wrong（注釈を確かめたスポットには当てない）。"""
    curation: dict = {}
    if apply_curation:
        flag_set = ds.curation_flags
        if flag_set is None:
            caveats.append("このデータセットはキュレーションの記録（.arf2 の判断）を持たないため、"
                           "キュレーションの当て方は何もしていません。")
        else:
            curation = {
                "wrong": {int(s) for s in flag_set["wrong"]},
                "redundant": {int(s) for s in flag_set["redundant"]},
                "assign": {int(s): {"name": row.get("name") or "", "ontology": row.get("ontology") or ""}
                           for s, row in flag_set["assign"].items()},
            }
            if flag_set.get("orphaned"):
                caveats.append(curation_flags.orphaned_warning(flag_set["orphaned"]))
    if exclude_auto_likely_wrong:
        arf2 = artifact_path(ds, "arf2")
        if arf2 is None:
            caveats.append("このデータセットには .arf2 が無いため、自動判定 likely_wrong は除いていません。")
            return curation
        try:
            sha = curation_flags.alignment_key(arf2)["alignment_sha256"]
            review = curation_suggest.latest_review(arf2, sha)
            if review is None:
                caveats.append("このアラインメントの curation_review のレビューが無いため、"
                               "自動判定 likely_wrong は除いていません。")
            else:
                # 注釈を確かめたスポット（confirmed の記録・_tags.xml の Confirmed）には自動判定を当てない
                effective = curation_flags.FlagStore(curation_flags.curation_dir(arf2)).effective(sha)
                confirmed = confirmed_spots(arf2, effective)
                curation["auto_likely_wrong"] = {
                    int(s["spot_id"]) for s in review["spots"]
                    if s.get("verdict") == "likely_wrong" and int(s["spot_id"]) not in confirmed}
        except curation_flags.FlagFileError as exc:
            raise DomainError("DATASET_BAD_REQUEST", str(exc), exc.details()) from exc
    return curation


def standard_only_spots(ds, standard: frozenset, plotted: set) -> frozenset[int]:
    """`standard` の試料にだけ実質的な強度があるスポット（`gi.standard_only_spots` と同じ規則を行列で）。

    標準液での最大が、描く試料での最大の `gi.STANDARD_ONLY_FOLD` 倍以上（かつ正）。標準液の試料が
    `dataset_exclude` で除外されていても判定に使う（標準液は描かないために指定するもの）。
    """
    if not standard or ds.feature_matrix is None:
        return frozenset()
    names = list(ds.sample_names)
    matrix = np.nan_to_num(np.asarray(ds.feature_matrix, dtype=float), nan=0.0)
    std_cols = [j for j, n in enumerate(names) if n in standard]
    grp_cols = [j for j, n in enumerate(names) if n in plotted]
    if not std_cols:
        return frozenset()
    std = matrix[:, std_cols].max(axis=1)
    grp = matrix[:, grp_cols].max(axis=1) if grp_cols else np.zeros(len(matrix))
    hit = (std > 0) & (std >= gi.STANDARD_ONLY_FOLD * np.maximum(grp, 1.0))
    index, _ = _spot_ids(ds)
    return frozenset(sid for sid, position in index.items() if hit[position])


# ---------- 項目の解決 ----------

@dataclass
class Resolved:
    items: list                        # `resolve_items` の戻り（item / parts / spots）
    excluded: dict
    caveats: list
    msms: dict | None
    index: dict                        # 整数の特徴 ID -> ds.feature_ids の位置


def resolve_selection(ds, items, *, apply_curation: bool = True, exclude_auto_likely_wrong: bool = False,
                      standard_only=frozenset(), caveats=None) -> Resolved:
    """項目（クラス名・分子種名・特徴 ID、`+` で合算）を特徴の集合に解く。

    除外した特徴（`ds.excluded_features`）はクラスにも名前指定にも入れない（excluded の `manual`）。
    クラス名・分子種名に当たらなかった部品が特徴 ID に一致すればその特徴 1 つに当てる（キュレーションの
    判断と手動除外は同じに効く）。
    """
    caveats = list(caveats or [])
    index, skipped = _spot_ids(ds)
    if not index:
        raise _bad("整数の ID を持つ特徴がありません（群別強度・分子種の図は MasterAlignmentID / SMF_ID が整数の"
                   "データセットだけ）。")
    if skipped:
        caveats.append(f"ID が整数でない特徴 {skipped} 件は選べないため除きました。")
    curation = _curation(ds, apply_curation=apply_curation,
                         exclude_auto_likely_wrong=exclude_auto_likely_wrong, caveats=caveats)
    manual = {sid for sid in index if ds.feature_ids[index[sid]] in (ds.excluded_features or {})}
    catalog = _catalog(ds, index)
    try:
        resolved, excluded = gi.resolve_items(
            items, catalog, curation=curation, standard_only=standard_only, manual=manual)
    except ValueError as exc:
        raise _bad(str(exc)) from exc
    _attach_feature_ids(resolved, excluded, catalog, curation, manual, index, ds)
    return Resolved(resolved, excluded, caveats, _msms(ds, index), index)


def _attach_feature_ids(resolved, excluded, catalog, curation, manual, index, ds) -> None:
    """どれにも当たらなかった部品が特徴 ID なら、その特徴を項目に足す（`resolve_items` は名前とクラスだけ解く）。"""
    by_id = {str(ds.feature_ids[position]): sid for sid, position in index.items()}
    row_of = {r["MasterAlignmentID"]: r for r in catalog}
    wrong = set(curation.get("wrong") or ()) | set(curation.get("redundant") or ())
    auto = set(curation.get("auto_likely_wrong") or ())
    assign = curation.get("assign") or {}
    extra = {key: set(values) for key, values in excluded.items()}
    for item in resolved:
        for part in item["parts"]:
            sid = by_id.get(part["part"]) if part["kind"] == "none" else None
            if sid is None:
                continue
            row = row_of[sid]
            name = (assign.get(sid) or {}).get("name") or row["Name"]
            ontology = ((assign.get(sid) or {}).get("ontology") or row["Ontology"]).strip()
            label = f"#{sid} {name}"
            reason = ("manual" if sid in manual else "curation" if sid in wrong
                      else "auto_likely_wrong" if sid in auto else None)
            if reason:
                extra[reason].add(label)
                continue
            part["kind"], part["n_spots"] = "feature", 1
            if all(s["spot_id"] != sid for s in item["spots"]):
                item["spots"].append({"spot_id": sid, "name": name, "ontology": ontology})
        item["spots"].sort(key=lambda s: s["spot_id"])
    for key, values in extra.items():
        excluded[key] = sorted(values)


def unmatched_parts(resolved: Resolved) -> list[str]:
    """どのスポットにも当たらなかった項目の部品（クラスでも分子種名でも特徴 ID でもない、または除外で空になった名前）。"""
    return [p["part"] for item in resolved.items for p in item["parts"] if p.get("kind") == "none"]


def expand_spots(ds, resolved: Resolved, *, require_msms: bool = False):
    """解決済みの項目をスポット単位に展開する（項目の順 → 表示名の順。同じスポットは 1 回）。

    `require_msms` で MS/MS の裏付けが無いスポットを落とす（落とした表示は戻り値の 2 つ目）。
    MS/MS の証拠の欄を持たないデータセット（mzTab-M だけ）では要求できない。
    """
    if require_msms and resolved.msms is None:
        raise _bad("このデータセットは MS/MS の裏付けの欄（.arf のセルの MS/MS の有無）を持たないため、"
                   "require_msms は使えません（MS-DIAL の解析フォルダを load_dataset で開いたときだけ）。")
    seen: set[int] = set()
    spots, dropped = [], []
    for item in resolved.items:
        ordered = sorted(item["spots"], key=lambda s: (display_name(s["name"]).casefold(), s["spot_id"]))
        for spot in ordered:
            sid = int(spot["spot_id"])
            if sid in seen:
                continue
            seen.add(sid)
            has_msms = bool((resolved.msms or {}).get(sid))
            if require_msms and not has_msms:
                dropped.append(f"#{sid} {spot['name']}")
                continue
            meta = (ds.feature_metadata or {}).get(ds.feature_ids[resolved.index[sid]]) or {}
            spots.append({
                "spot_id": sid, "name": spot["name"], "label": display_name(spot["name"]),
                "ontology": spot["ontology"], "item": item["item"],
                "adduct": meta.get("adduct") or "",
                "mz": round(float(meta.get("mz") or 0.0), 4),
                "rt": round(float(meta.get("rt") or 0.0), 3),
                "msms": has_msms,
            })
    return spots, sorted(dropped)


def select_items(ds, items: list[str], *, apply_curation: bool, exclude_auto_likely_wrong: bool,
                 require_msms: bool, standard_only=frozenset()) -> list[str]:
    """項目（クラス名・分子種名・特徴 ID）を特徴 ID（`ds.feature_ids` の表記）の並びにする。

    クラス名は `feature_metadata["ontology"]` で配下の特徴に展開する。除外した特徴
    （`ds.excluded_features`）・キュレーションの判断（`apply_curation` / `exclude_auto_likely_wrong`）・
    MS/MS の要求（`require_msms`）を当てたあとの特徴だけを、項目の順 → 表示名の順で、重複なしで返す。
    `standard_only` は標準液だけにある特徴（`standard_only_spots`）で、クラスの展開から外す。
    """
    resolved = resolve_selection(ds, items, apply_curation=apply_curation,
                                 exclude_auto_likely_wrong=exclude_auto_likely_wrong,
                                 standard_only=standard_only)
    spots, _ = expand_spots(ds, resolved, require_msms=require_msms)
    return [str(ds.feature_ids[resolved.index[s["spot_id"]]]) for s in spots]


# ---------- 行（特徴 × 試料） ----------

def rows_by_spot(ds, feature_ids, *, groups: list[str], low_reliability_samples=None,
                 standard_samples=None) -> dict[int, list[dict]]:
    """選んだ特徴の、描く試料の行 `{spot_id: [{"file_name", "height", "is_gap_filled"}]}`。

    `build_group_intensity_payload` / `build_species_payload` が食う形。値は `feature_matrix`（NaN のセルは
    行にしない）、gap-fill の印は `detected_mask`（False が gap-filled。マスクが無ければ False のまま）。
    行を作る試料は、`groups` に入る試料と `low_reliability_samples` の試料（どちらも除外した試料を除く）。
    `standard_samples` は描かないので行にしない（標準液だけにあるスポットは `standard_only_spots` が判定する）。
    """
    selection = resolve_samples(ds, groups, low_reliability_samples, standard_samples)
    wanted = {n for g in selection.groups for n in g["samples"]} | set(selection.low_reliability)
    names = list(ds.sample_names)
    columns = [(j, n) for j, n in enumerate(names) if n in wanted]
    index, _ = _spot_ids(ds)
    matrix = np.asarray(ds.feature_matrix, dtype=float)
    mask = ds.detected_mask
    out: dict[int, list[dict]] = {}
    for fid in feature_ids:
        sid = int(fid)
        i = index[sid]
        rows = []
        for j, name in columns:
            value = matrix[i, j]
            if not math.isfinite(value):
                continue
            rows.append({"file_name": name, "height": float(value),
                         "is_gap_filled": bool(mask is not None and not mask[i, j])})
        out[sid] = rows
    return out
