"""候補付けレビューの組み立て・保存・要約・送信内容の展開。spec §3〜§8。

1 回の候補付け = 対象（wrong フラグ・likely_wrong・未注釈）のスポットごとの ①②（candidates）と
④（relations）。`<arf2 のフォルダ>/curation/` に `suggest-<id>.json`（正準）と `.html` を書く。
送信用テキストは候補 ID だけを運び、記録行の中身はここで保存済みの候補から展開する。
deps: curation.{evidence,candidates,relations,flags,review,trend,viewer}、arf2.{reader,match_results,
ion_features}、msdial.analysis_params。tools_* は import しない。
"""
from __future__ import annotations

import json
import re
import secrets
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.arf2.ion_features import load_ion_features
from metabolomix.arf2.match_results import load_spot_annotations, load_spot_candidates
from metabolomix.arf2.reader import load_catalog
from metabolomix.core.atomic_io import atomic_write_json
from metabolomix.curation import candidates, evidence, flags, relations, review, trend, viewer
from metabolomix.curation.apply import confirmed_spots
from metabolomix.eic.reader import read_eic_spot_css1
from metabolomix.library.defaults import pick_tol
from metabolomix.msdial.analysis_params import resolve_analysis_params

SUGGESTION_ID_RE = re.compile(r"^cs-\d{8}-\d{6}-[0-9a-f]{4}$")
#: 相手 Y の EIC を持たせる関係の数（スポットごと、並び順の上から）。JSON を膨らませないため。
PARTNER_EIC_RELATIONS = 3
WRONG_MODES = ("flagged_or_likely", "flagged")
LEVELS = ("sum", "species")
DEFAULTS = {"wrong": "flagged_or_likely", "unannotated": True, "include_decided": False, "top_n": 5,
            "relation_mz_tol": 0.010, "relation_min_r": 0.8}
TSV_COLUMNS = ["spot_id", "name", "target", "top_candidate", "source", "total_score", "reasons",
               "relation", "strong"]


def new_suggestion_id() -> str:
    return f"cs-{datetime.now():%Y%m%d-%H%M%S}-{secrets.token_hex(2)}"


def is_valid_suggestion_id(value) -> bool:
    return isinstance(value, str) and SUGGESTION_ID_RE.fullmatch(value) is not None


def latest_review(arf2_path, alignment_sha256: str) -> dict | None:
    """そのアラインメント（sha256 一致）の最新のレビュー。ID の時刻順に新しいものから読む。"""
    directory = flags.curation_dir(arf2_path)
    for path in sorted(directory.glob("review-cr-*.json"), key=lambda p: p.name, reverse=True):
        review_id = path.stem[len("review-"):]
        if not review.is_valid_review_id(review_id):
            continue
        try:
            saved = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue    # 壊れた・読めないレビュー 1 件で候補付け全体を止めない
        if isinstance(saved, dict) and (saved.get("alignment") or {}).get("alignment_sha256") == alignment_sha256:
            return saved
    return None


def select_targets(catalog, annotations, base_review, effective, *, wrong, unannotated, include_decided,
                   confirmed=frozenset()):
    """`confirmed`（注釈を確かめたスポット）は likely_wrong の対象にしない。"""
    decisions = flags.split_decisions(effective)
    decided = set(decisions["assign"]) | set(decisions["redundant"])
    verdicts = {s["spot_id"]: s for s in base_review["spots"]}
    targets = []
    for spot in catalog:
        spot_id = spot["MasterAlignmentID"]
        if spot_id in decided and not include_decided:
            continue
        annotated = (annotations.get(spot_id) or {}).get("representative") is not None
        judged = verdicts.get(spot_id)
        if annotated and (spot_id in decisions["wrong"] or spot_id in decided):
            # include_decided=True のとき、判断済み（assign / redundant）の注釈付きスポットも "flagged" で戻す
            targets.append({"spot": spot, "target_kind": "flagged",
                            "target_reasons": list((judged or {}).get("reasons") or [])})
        elif (annotated and wrong == "flagged_or_likely" and judged and judged["verdict"] == "likely_wrong"
              and spot_id not in confirmed):
            targets.append({"spot": spot, "target_kind": "likely_wrong",
                            "target_reasons": list(judged.get("reasons") or [])})
        elif not annotated and unannotated:
            targets.append({"spot": spot, "target_kind": "unannotated", "target_reasons": []})
    return targets


def _heights(rows: list[dict]) -> list:
    return [row.get("height") for row in sorted(rows, key=lambda r: r.get("file_id") or 0)]


def _partner(spot: dict, annotation: dict, rows, assigned: dict | None = None) -> dict:
    """`assigned` は有効な assign 行。あればその name / ontology / formula / adduct（空なら catalog の値）を使う
    （置き換えられた MS-DIAL の名前や組成式で同位体の期待比などを計算しない）。"""
    rep = annotation.get("representative") or {}
    row = assigned or {}
    return {"spot_id": spot["MasterAlignmentID"],
            "name": row.get("name") or spot.get("Name") or rep.get("name"),
            "ontology": row.get("ontology") or spot.get("Ontology") or "",
            "mz": spot.get("MassCenter"), "rt": spot.get("RT"),
            "adduct": row.get("adduct") or spot.get("AdductType"),
            "formula": row.get("formula") or spot.get("Formula"), "heights": _heights(rows)}


def _partner_trace(eic_path, spot_id: int, annotations: dict) -> dict | None:
    """相手 Y の代表試料の EIC 1 本（ビューアの重ね描き用。spec §7.2）。読めなければ None。"""
    rep_file = (annotations.get(spot_id) or {}).get("representative_file_id")
    if eic_path is None or rep_file is None:
        return None
    try:
        eic = read_eic_spot_css1(str(eic_path), spot_id, [rep_file], max_traces=1, max_total_points=5000)
    except (ValueError, OSError):
        return None
    if not eic["samples"]:
        return None
    sample = eic["samples"][0]
    left, right = sample["peak_left"], sample["peak_right"]
    points = evidence._trim(sample["chromatogram"], left, right)
    return {"file_id": sample["file_id"], "left": round(left, 4), "top": round(sample["peak_top"], 4),
            "right": round(right, 4), "points": evidence._downsample_points(points, left, right)}


def _scoring(store) -> dict:
    sp = store.summary().get("search_params") or {}
    return {"mz_tol": pick_tol(None, sp, "ms1_tolerance", 0.01),
            "ms2_tol": pick_tol(None, sp, "ms2_tolerance", 0.025),
            "rt_tol": pick_tol(None, sp, "rt_tolerance", 0.2),
            "mass_begin": pick_tol(None, sp, "mass_range_begin", 0.0),
            "mass_end": pick_tol(None, sp, "mass_range_end", 2000.0),
            "relative_amp_cutoff": pick_tol(None, sp, "relative_amp_cutoff", 0.0),
            "absolute_amp_cutoff": pick_tol(None, sp, "absolute_amp_cutoff", 0.0),
            "use_rt": bool(sp.get("use_time_for_annotation_scoring", False))}


def run_suggestion(arf2_path, *, base_review, store, th, options) -> dict:
    arf2_path = Path(arf2_path).resolve()
    opts = {**DEFAULTS, **{k: v for k, v in (options or {}).items() if v is not None}}
    if opts["wrong"] not in WRONG_MODES:
        raise ValueError(f"wrong は {WRONG_MODES} のいずれかにしてください。")
    alignment = flags.alignment_key(arf2_path)
    if base_review["alignment"]["alignment_sha256"] != alignment["alignment_sha256"]:
        raise ValueError("元のレビューは別の版のアラインメント（.arf2）に対して作られています。"
                         "curation_review をやり直してください。")
    effective = flags.FlagStore(flags.curation_dir(arf2_path)).effective(alignment["alignment_sha256"])
    catalog = load_catalog(arf2_path)
    annotations = load_spot_annotations(arf2_path)
    all_candidates = load_spot_candidates(arf2_path)
    ion_features = load_ion_features(arf2_path)
    links = {spot_id: f["links"] for spot_id, f in ion_features.items()}
    ion_mode = next((s.get("IonMode") for s in catalog if s.get("IonMode")), None)
    params = resolve_analysis_params(arf2_path, ion_mode)
    rt_window = opts.get("rt_window") or params["rt_window"]
    warnings = []
    if params["source"] == "default":
        warnings.append("解析の param ファイル（*_param_*.txt）が見つからないので、既定の検索アダクトと "
                        f"RT 窓 {rt_window} 分を使いました。")

    targets = select_targets(catalog, annotations, base_review, effective, wrong=opts["wrong"],
                             unannotated=opts["unannotated"], include_decided=opts["include_decided"],
                             confirmed=confirmed_spots(arf2_path, effective))
    if len(targets) > evidence.MAX_SPOTS:
        raise ValueError(f"対象が {len(targets)} 件あり、上限 {evidence.MAX_SPOTS} を超えています。")
    target_ids = {t["spot"]["MasterAlignmentID"] for t in targets}
    decisions = flags.split_decisions(effective)
    excluded_partners = target_ids | decisions["wrong"] | set(decisions["redundant"])

    files = evidence.sibling_files(arf2_path)
    rows_by_spot = evidence.arf_rows_by_spot(files["arf"])
    if files["arf"] is None:
        warnings.append("PeakProperties.arf が無いので試料間の相関を計算できません（MS-DIAL のリンクだけで裏付けを判定）。")
    partners, points = [], []
    for spot in catalog:
        spot_id = spot["MasterAlignmentID"]
        annotation = annotations.get(spot_id) or {}
        if annotation.get("representative") is None or spot_id in excluded_partners:
            continue
        partner = _partner(spot, annotation, rows_by_spot.get(spot_id, []), decisions["assign"].get(spot_id))
        partners.append(partner)
        comp = trend.composition(partner["name"])
        if comp is not None and spot.get("RT") is not None:
            points.append({"spot_id": spot_id, "ontology": partner["ontology"], "rt": spot["RT"],
                           "mz": spot.get("MassCenter"), "carbon": comp[0], "db": comp[1]})
    trends = trend.fit_trends(points, th)

    scoring = _scoring(store)
    evs, stats = evidence.collect(arf2_path, [t["spot"] for t in targets], store=store,
                                  ms2_tol=scoring["ms2_tol"], th=th, keep_measured=True)
    by_id = {t["spot"]["MasterAlignmentID"]: t for t in targets}
    spots_out, n_hard = [], 0
    for ev in evs:
        target = by_id[ev["spot_id"]]
        annotation = annotations.get(ev["spot_id"]) or {}
        rep = annotation.get("representative")
        msdial = [c for c in all_candidates.get(ev["spot_id"], []) if c is not rep and c != rep]
        lib = candidates.build_library_candidates(
            ev, msdial_matches=msdial, representative=rep, store=store, scoring=scoring, trends=trends,
            th=th, exclude_current=target["target_kind"] != "unannotated", top_n=opts["top_n"])
        n_hard += lib["n_hard_removed"]
        me = _partner(target["spot"], annotation, rows_by_spot.get(ev["spot_id"], []))
        ion = relations.find_relations(me, partners, adducts=params["searched_adducts"], rt_window=rt_window,
                                       mz_tol=opts["relation_mz_tol"], min_r=opts["relation_min_r"],
                                       links=links)
        for relation in ion[:PARTNER_EIC_RELATIONS]:
            relation["partner_eic"] = _partner_trace(files["eic"], relation["of"], annotations)
        strong = next((r for r in ion if r["strong"]), None)
        ev.pop("_measured", None)
        spots_out.append({
            "spot_id": ev["spot_id"], "name": ev["name"], "ontology": ev["ontology"], "adduct": ev["adduct"],
            "ion_mode": ev["ion_mode"], "mz": ev["mz"], "rt": ev["rt"], "rep_mz": ev["rep_mz"],
            "rep_rt": ev["rep_rt"], "target_kind": target["target_kind"],
            "target_reasons": target["target_reasons"],
            "current": None if rep is None else {"name": rep.get("name"), "inchikey": rep.get("inchikey")},
            "eic": ev["eic"], "measured": lib["measured"], "candidates": lib["candidates"],
            "relations": ion, "strong": strong is not None,
            "preset": strong["candidate_id"] if strong else None})

    spots_out.sort(key=_spot_order)      # ビューアのカードも summary_tsv と同じ並びにする
    counts = {"targets": {k: sum(1 for s in spots_out if s["target_kind"] == k)
                          for k in ("flagged", "likely_wrong", "unannotated")},
              "with_candidates": sum(1 for s in spots_out if s["candidates"] or
                                     any(not r["informational"] for r in s["relations"])),
              "with_strong_relation": sum(1 for s in spots_out if s["strong"]),
              "hard_removed": n_hard}
    if stats["missing_files"]:
        warnings.append(f"兄弟ファイルがありません: {', '.join(stats['missing_files'])}。")
    warnings.append("likely_wrong は元のレビュー（" + base_review["review_id"] + "）の対象範囲についてだけ分かります。")
    return {"suggestion_id": new_suggestion_id(), "created_at": datetime.now(timezone.utc).isoformat(),
            "arf2_path": str(arf2_path), "alignment": alignment, "base_review_id": base_review["review_id"],
            "library": {"sha256": store.summary().get("source_sha256")},
            "thresholds": th, "options": {**opts, "rt_window": rt_window}, "analysis_params": params,
            "scoring": scoring, "counts": counts, "warnings": warnings, "spots": spots_out}


def save_suggestion(s: dict, submit_endpoint: dict | None = None) -> dict:
    directory = flags.curation_dir(s["arf2_path"])
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / f"suggest-{s['suggestion_id']}.json"
    html_path = directory / f"suggest-{s['suggestion_id']}.html"
    atomic_write_json(json_path, s)
    html_path.write_text(viewer.render_suggest_html(s, submit_endpoint), encoding="utf-8")
    return {"json": json_path, "html": html_path}


def load_suggestion(arf2_path_or_dir, suggestion_id: str) -> dict:
    if not is_valid_suggestion_id(suggestion_id):
        raise ValueError(f"suggestion_id={suggestion_id!r} の形が不正です（cs-YYYYMMDD-HHMMSS-xxxx）。")
    base = Path(arf2_path_or_dir)
    directory = base if base.is_dir() else flags.curation_dir(base)
    path = directory / f"suggest-{suggestion_id}.json"
    if not path.is_file():
        raise FileNotFoundError(f"suggestion_id={suggestion_id} の候補付けがありません: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


_KIND_ORDER = {"flagged": 0, "likely_wrong": 1, "unannotated": 2}


def _spot_order(sp: dict) -> tuple:
    """強い説明が先、次に wrong フラグ → likely_wrong → 未注釈、最後に spot_id。"""
    return (not sp["strong"], _KIND_ORDER[sp["target_kind"]], sp["spot_id"])


def summary_tsv(s: dict, max_rows: int | None = None) -> str:
    spots = sorted(s["spots"], key=_spot_order)
    if max_rows is not None:
        spots = spots[:max(0, max_rows)]
    lines = ["\t".join(TSV_COLUMNS)]
    for sp in spots:
        top = sp["candidates"][0] if sp["candidates"] else None
        rel = sp["relations"][0] if sp["relations"] else None
        cells = [sp["spot_id"], sp["name"], sp["target_kind"],
                 top and (top["sum_name"] or top["name"]), top and top["source"],
                 top and top["scores"]["total_score"], top and ",".join(top["soft"] + top["info"]),
                 rel and f"{rel['relation']}(#{rel['of']})", sp["strong"]]
        lines.append("\t".join(review._cell(c) for c in cells))
    return "\n".join(lines)


def _note(entry) -> str:
    note = str(entry.get("note") or "")
    for control in ("\r", "\n", "\t"):
        note = note.replace(control, " ")
    return note[:500]


def expand_entries(entries, suggestion: dict) -> list[dict]:
    if not isinstance(entries, list) or not entries:
        raise ValueError("flags は 1 件以上のリストで渡してください。")
    spots = {sp["spot_id"]: sp for sp in suggestion["spots"]}
    out = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"flags[{index}] が dict ではありません。")
        spot_id, flag = entry.get("spot_id"), entry.get("flag")
        if isinstance(spot_id, bool) or not isinstance(spot_id, int):
            raise ValueError(f"flags[{index}].spot_id が整数ではありません: {spot_id!r}")
        if spot_id not in spots:
            raise ValueError(f"flags[{index}].spot_id={spot_id} はこの候補付けの対象外です。")
        if flag not in flags.SUGGEST_FLAG_VALUES:
            raise ValueError(f"flags[{index}].flag={flag!r} は {flags.SUGGEST_FLAG_VALUES} のいずれかにしてください。")
        base = {"spot_id": spot_id, "flag": flag, "note": _note(entry)}
        if flag == "clear":
            out.append(base)
            continue
        pool = spots[spot_id]["candidates"] if flag == "assign" else spots[spot_id]["relations"]
        chosen = next((c for c in pool if c["candidate_id"] == entry.get("candidate")), None)
        if chosen is None or (flag == "redundant" and chosen.get("informational")):
            raise ValueError(f"flags[{index}].candidate={entry.get('candidate')!r} はスポット {spot_id} の"
                             f"{'候補' if flag == 'assign' else 'イオン関係'}にありません。")
        if flag == "assign":
            level = entry.get("level") or "sum"
            if level not in LEVELS:
                raise ValueError(f"flags[{index}].level={level!r} は {LEVELS} のいずれかにしてください。")
            name = (chosen["sum_name"] or chosen["name"]) if level == "sum" else chosen["name"]
            out.append({**base, "name": name, "level": level, "species_name": chosen["name"],
                        "ontology": chosen.get("ontology"), "adduct": chosen.get("adduct"),
                        "formula": chosen.get("formula"), "inchikey": chosen.get("inchikey"),
                        "candidate_source": chosen["source"],
                        "library": {"sha256": suggestion["library"].get("sha256"),
                                    "library_id": chosen.get("library_id"),
                                    "record_index": chosen.get("record_index")},
                        "scores": chosen["scores"], "suggestion_id": suggestion["suggestion_id"]})
        else:
            out.append({**base, "of": chosen["of"], "relation": chosen["relation"],
                        "evidence": {k: chosen.get(k) for k in ("dmz_mda", "drt", "r", "msdial_links")},
                        "suggestion_id": suggestion["suggestion_id"]})
    return out
