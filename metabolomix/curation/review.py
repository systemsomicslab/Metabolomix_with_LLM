"""レビューの生成・保存・要約。spec §3・§5。

1 回のレビュー = 対象スポットの証拠＋判別＋傾向＋既存フラグ。`<arf2 のフォルダ>/curation/`
に `review-<id>.json`（正準）と `review-<id>.html`（ビューア）を書く。LLM へ返すのは
`summary_tsv` の要約だけで、EIC 系列やスペクトル座標は返さない。
"""
from __future__ import annotations

import json
import re
import secrets
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.arf2.match_results import load_spot_annotations
from metabolomix.arf2.reader import load_catalog
from metabolomix.core.atomic_io import atomic_write_json
from metabolomix.curation import adduct_isomer, evidence, flags, trend, viewer
from metabolomix.curation.judge import auto_note, judge_spot, lipid_rules_active

VERDICT_RANK = {"ok": 0, "suspect": 1, "likely_wrong": 2}
REFERENCE_WARN_FRACTION = 0.5
TSV_COLUMNS = ["spot_id", "name", "ontology", "adduct", "verdict", "reasons",
               "ppm", "dmz_mda", "drt", "wdot", "mpp", "eic_good", "trend_z", "flag"]


#: `new_review_id()` の形。review_id はファイル名（`review-<id>.json`）に使うので、
#: 外から来た値はパスに使う前に必ずこの形かを検める（`../` などを通さない）。
REVIEW_ID_RE = re.compile(r"^cr-\d{8}-\d{6}-[0-9a-f]{4}$")


def new_review_id() -> str:
    return f"cr-{datetime.now():%Y%m%d-%H%M%S}-{secrets.token_hex(2)}"


def is_valid_review_id(review_id) -> bool:
    return isinstance(review_id, str) and REVIEW_ID_RE.fullmatch(review_id) is not None


def run_review(arf2_path, spots, *, store, ms2_tol, th, file_ids, max_traces, selection) -> dict:
    """1 回のレビューを組み立てる。既存フラグは重い証拠収集の**前**に読む——`flags.jsonl`
    が壊れていれば（`flags.FlagFileError`）数分の計算を捨てずに先に止めるため。"""
    arf2_path = Path(arf2_path).resolve()
    alignment = flags.alignment_key(arf2_path)
    flag_rows = flags.FlagStore(flags.curation_dir(arf2_path)).rows()
    existing = flags.effective_flags(flag_rows, alignment["alignment_sha256"])
    orphaned = flags.orphaned_count(flag_rows, alignment)
    cleared = flags.cleared_spots(flag_rows, alignment["alignment_sha256"])
    evs, stats = evidence.collect(arf2_path, spots, store=store, ms2_tol=ms2_tol, th=th,
                                  file_ids=file_ids, max_traces=max_traces)
    points = []
    for ev in evs:
        comp = trend.composition(ev["name"])
        if comp is not None and ev["rt"] is not None and ev["mz"] is not None:
            points.append({"spot_id": ev["spot_id"], "ontology": ev["ontology"] or "",
                           "rt": ev["rt"], "mz": ev["mz"], "carbon": comp[0], "db": comp[1]})
        ev["composition"] = list(comp) if comp else None
    trends = trend.fit_trends(points, th)
    pool = _adduct_isomer_pool(arf2_path)
    wrong = flags.split_decisions(existing)["wrong"]
    partners = [p for spot_id, p in pool.items() if spot_id not in wrong]

    counts = {"ok": 0, "suspect": 0, "likely_wrong": 0}
    lipid_rules = lipid_rules_active(evs)
    for ev in evs:
        target = pool.get(ev["spot_id"])
        ev["adduct_isomer"] = (adduct_isomer.find_adduct_isomer(target, partners, th)
                               if target is not None else None)
        ev.update(judge_spot(ev, trends["spots"].get(ev["spot_id"]), th, lipid_rules=lipid_rules))
        ev["trend"] = trends["spots"].get(ev["spot_id"])
        row = existing.get(ev["spot_id"]) or {}
        is_flag = row.get("flag") in ("wrong", "suspect")
        ev["flag"] = row.get("flag") if is_flag else None
        ev["flag_note"] = row.get("note") if is_flag else None
        ev["decision"] = ({"flag": "assign", "name": row.get("name")} if row.get("flag") == "assign"
                          else {"flag": "redundant", "of": row.get("of")} if row.get("flag") == "redundant"
                          else None)
        ev["flag_cleared"] = ev["spot_id"] in cleared
        ev["auto_note"] = auto_note(ev, th)
        counts[ev["verdict"]] += 1

    warnings = []
    if stats["n_with_match"] and stats["n_reference_resolved"] / stats["n_with_match"] < REFERENCE_WARN_FRACTION:
        warnings.append(
            f"照合結果を持つ {stats['n_with_match']} 件のうち参照を引けたのは "
            f"{stats['n_reference_resolved']} 件です。アラインメントに使われたものと別のライブラリを"
            "読んでいる可能性があります（同じフォルダの *_Loaded.msp2.dbs を library_load してください）。")
    if stats["missing_files"]:
        warnings.append(f"兄弟ファイルがありません: {', '.join(stats['missing_files'])}"
                        "（その系統の判別は UNKNOWN になります）。")
    if orphaned:
        warnings.append(flags.orphaned_warning(orphaned))
    if stats["n_with_match"] and not lipid_rules:
        warnings.append(
            "MS-DIAL の脂質規則フラグ（IsLipidClassMatch / IsLipidChainsMatch / IsOtherLipidMatch）が"
            "どのスポットにも立っていません（脂質以外の採点器の出力とみなしました）。"
            "規則による判定（class_rule_rejected・class_rules_not_run・chains_unsupported）は使っていません。")

    return {
        "review_id": new_review_id(),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "arf2_path": str(arf2_path), "alignment": alignment, "selection": selection,
        "thresholds": th, "ms2_tol": ms2_tol, "counts": counts, "stats": stats,
        "warnings": warnings, "trend": {"classes": trends["classes"],
                                        "groups": {k: {str(db): g for db, g in v.items()}
                                                   for k, v in trends["groups"].items()}},
        "spots": evs,
    }


def _adduct_isomer_pool(arf2_path) -> dict[int, dict]:
    """別アダクトの取り違えを比べる相手の候補。レビューの対象に絞らず、アラインメントの注釈付き
    スポット全部（`ontology=["PI"]` のレビューでも DGDG の相手を見つけるため）。有効な `wrong`
    フラグのあるスポットを相手から外すのは呼び出し側（対象自身はここから引く）。"""
    annotations = load_spot_annotations(arf2_path)
    pool = {}
    for spot in load_catalog(arf2_path):
        entry = adduct_isomer.pool_entry(spot, annotations.get(spot["MasterAlignmentID"]))
        if entry is not None:
            pool[entry["spot_id"]] = entry
    return pool


def save_review(review: dict) -> dict:
    directory = flags.curation_dir(review["arf2_path"])
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / f"review-{review['review_id']}.json"
    html_path = directory / f"review-{review['review_id']}.html"
    atomic_write_json(json_path, review)
    html_path.write_text(viewer.render_html(review), encoding="utf-8")
    return {"json": json_path, "html": html_path}


def load_review(arf2_path_or_dir, review_id: str) -> dict:
    """保存済みレビュー。review_id が `REVIEW_ID_RE` の形でなければ ValueError（パスに使わない）。"""
    if not is_valid_review_id(review_id):
        raise ValueError(f"review_id={review_id!r} の形が不正です（cr-YYYYMMDD-HHMMSS-xxxx）。")
    base = Path(arf2_path_or_dir)
    directory = base if base.is_dir() else flags.curation_dir(base)
    path = directory / f"review-{review_id}.json"
    if not path.is_file():
        raise FileNotFoundError(f"review_id={review_id} のレビューがありません: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return str(value).replace("\t", " ").replace("\n", " ")


def _summary_spots(review: dict, min_verdict: str) -> list[dict]:
    """表に載せるスポット（`min_verdict` 以上とフラグ済み）を判定の重い順に。"""
    floor = VERDICT_RANK[min_verdict]
    return [s for s in sorted(review["spots"],
                              key=lambda s: (-VERDICT_RANK[s["verdict"]], s["spot_id"]))
            if VERDICT_RANK[s["verdict"]] >= floor or s.get("flag")]


def n_summary_rows(review: dict, *, min_verdict: str = "suspect") -> int:
    return len(_summary_spots(review, min_verdict))


def summary_tsv(review: dict, *, min_verdict: str = "suspect", max_rows: int | None = None) -> str:
    """`min_verdict` 以上とフラグ済みのスポットの TSV。`max_rows` があれば重い順の先頭だけ
    （実データでは約 1000 行になり LLM の文脈を占めるため。全件は HTML ビューアにある）。"""
    lines = ["\t".join(TSV_COLUMNS)]
    spots = _summary_spots(review, min_verdict)
    if max_rows is not None:
        spots = spots[:max(0, max_rows)]
    for s in spots:
        match = s.get("match") or {}
        weighted = match.get("squared_weighted_dot_product")
        lines.append("\t".join(_cell(v) for v in [
            s["spot_id"], s["name"], s["ontology"], s["adduct"], s["verdict"],
            ",".join(s["reasons"]), s["ppm"], s.get("dmz_mda"), s["drt"],
            round(weighted ** 0.5, 3) if weighted is not None and weighted >= 0 else None,
            match.get("matched_peaks_percentage"),
            (s.get("eic_shape") or {}).get("good_fraction"),
            (s.get("trend") or {}).get("z"), s.get("flag")]))
    return "\n".join(lines)


def trend_summary(review: dict) -> dict:
    """クラス別傾向の要約（点数・外れ値を除いた R²・外れ数。spec §5.1）。"""
    return {name: {"n": c["n"], "r2": c["r2"], "n_outliers": c.get("n_outliers", 0)}
            for name, c in review["trend"]["classes"].items()}


def page(review: dict, page_index: int, page_size: int = 50) -> dict:
    spots = review["spots"]
    n_pages = max(1, -(-len(spots) // page_size))
    start = page_index * page_size
    return {"spots": spots[start:start + page_size], "page": page_index, "n_pages": n_pages}
