"""アラインメントのキュレーション: curation_review / curation_submit / curation_flags /
curation_view_data と、MCP Apps 用の ui:// ビューア。spec 2026-09-28。

deps: mcp_core / session_state / mcp_errors / serialization / path_resolvers / curation.*。
lipidmix.tools.* の他モジュールと server は import しない。
"""
from __future__ import annotations

from pathlib import Path

from mcp.types import ToolAnnotations

from lipidmix.arf2.reader import load_catalog
from lipidmix.core import mcp_errors, session_state
from lipidmix.core.mcp_core import mcp
from lipidmix.core.path_resolvers import resolve_arf2_file_path
from lipidmix.core.serialization import json_payload, round_floats
# モジュール名を flag_log にするのは、curation_submit の引数 `flags`（公開 API の名前）が
# モジュールを隠すため。
from lipidmix.curation import evidence, judge, review, viewer
from lipidmix.curation import flags as flag_log
from lipidmix.library.defaults import DEFAULT_MS2_TOL, pick_tol

__all__ = ["curation_review", "curation_submit", "curation_flags", "curation_view_data"]

VIEWER_URI = "ui://ms-data-parser/curation-viewer"
_UI_META = {"ui": {"resourceUri": VIEWER_URI}, "ui/resourceUri": VIEWER_URI}
_APP_ONLY_META = {"ui": {"resourceUri": VIEWER_URI, "visibility": ["app"]}}
_LOCAL_WRITE = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True)
_LOCAL_WRITE_APPEND = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False)


def _error(message: str, **details) -> str:
    return json_payload({"status": "error", "message": message, **details})


def _review_dir(review_id: str) -> Path | None:
    known = session_state.session.curation.review_dirs.get(review_id)
    if known:
        return Path(known)
    arf2 = resolve_arf2_file_path(None)
    return flag_log.curation_dir(arf2) if arf2 else None


@mcp.resource(VIEWER_URI, name="curation-viewer", mime_type="text/html;profile=mcp-app",
              description="キュレーションの一覧ビューア（MCP Apps）。データは curation_view_data から取る。")
def curation_viewer_resource() -> str:
    return viewer.render_html(None)


@mcp.tool(annotations=_LOCAL_WRITE, structured_output=False, meta=_UI_META)
def curation_review(ontology: list[str] | None = None, name_contains: str | None = None,
                    file_ids: list[int] | None = None, max_traces: int = 12,
                    thresholds: dict | None = None, file_path: str | None = None) -> str:
    """注釈付きスポットを一覧で確かめるレビューを作る（EIC・対向プロット・Δppm・ΔRT・
    RT–m/z 傾向と機械判別）。**先に library_load**（アラインメントに使われた
    `*_Loaded.msp2.dbs` を推奨）。

    対象: 既定は注釈付き全部。`ontology=["PG"]` でクラス、`name_contains` で名前の部分一致。
    戻り値は suspect 以上(とフラグ済み)のスポットだけの TSV、判定の件数、クラス別の
    傾向要約、HTML ビューアのパス。**EIC 系列やスペクトルは返さない**——ユーザーには
    `html_path` をブラウザで開いてもらい、ビューアで付けたフラグを「送信用テキストを
    コピー」→ チャットに貼ってもらう。貼られたら curation_submit(submission_text=...) に渡す。

    判定: likely_wrong（強い不一致: ppm_out / polarity_mismatch / precursor_unmatched）、
    suspect（low_score / drt_out / eic_poor、または弱い兆候の重なり）、ok。
    RT–m/z 傾向（trend_outlier）は補強にしか使わない。UNKNOWN は不一致に数えない。
    `thresholds` で既定のしきい値（ppm_pass=5, ppm_borderline=10, drt_pass=0.5 分 など）を上書きできる。
    """
    arf2_path = resolve_arf2_file_path(file_path)
    if not arf2_path:
        return mcp_errors.missing_state(
            "arf2_file", ["load_dataset", "arf2_parser"],
            ".arf2 が見つかりません。先に load_dataset で MS-DIAL の出力フォルダを指定してください。")
    store = session_state.session.library.store
    if store is None:
        return mcp_errors.missing_state(
            "library", ["library_load"],
            "参照ライブラリが読み込まれていません。先に library_load を実行してください"
            "（アラインメントと同じフォルダの *_Loaded.msp2.dbs を推奨）。")
    try:
        th = judge.resolve_thresholds(thresholds)
    except ValueError as exc:
        return _error(str(exc))

    spots = evidence.select_spots(load_catalog(arf2_path), ontology=ontology,
                                  name_contains=name_contains)
    if not spots:
        return _error("条件に合う注釈付きスポットがありません。",
                      selection={"ontology": ontology, "name_contains": name_contains})
    if len(spots) > evidence.MAX_SPOTS:
        return _error(f"対象が {len(spots)} 件あり、1 回のレビューの上限 {evidence.MAX_SPOTS} を"
                      "超えています。ontology か name_contains で絞ってください。")

    search_params = store.summary().get("search_params") or {}
    ms2_tol = pick_tol(None, search_params, "ms2_tolerance", DEFAULT_MS2_TOL)
    result = review.run_review(arf2_path, spots, store=store, ms2_tol=ms2_tol, th=th,
                               file_ids=file_ids, max_traces=max_traces,
                               selection={"ontology": ontology, "name_contains": name_contains})
    saved = review.save_review(result)
    session_state.session.curation.last_review_id = result["review_id"]
    session_state.session.curation.review_dirs[result["review_id"]] = str(saved["json"].parent)

    return json_payload(round_floats({
        "review_id": result["review_id"],
        "n_spots": len(result["spots"]),
        "counts": result["counts"],
        "warnings": result["warnings"],
        "trend": {name: {"n": c["n"], "r2": c["r2"]} for name, c in result["trend"]["classes"].items()},
        "table": review.summary_tsv(result),
        "html_path": str(saved["html"]),
        "thresholds": th, "ms2_tol": ms2_tol,
    }, 4))


@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False)
def curation_submit(submission_text: str | None = None, review_id: str | None = None,
                    flags: list[dict] | None = None, source: str = "user") -> str:
    """キュレーションのフラグを記録する（追記。フラグの無いスポットは「間違っていない」で何も書かない）。

    ユーザーがビューアの「送信用テキストをコピー」で貼った文をそのまま `submission_text` に渡す
    （`CURATION_SUBMIT ` で始まる行だけを読む。書き写さないこと）。直接渡すなら `review_id` と
    `flags=[{"spot_id": 12, "flag": "wrong" | "suspect" | "clear", "note": "..."}]`。
    `source` は誰の判断か: ユーザー自身の判断は "user"、LLM の提案にユーザーがチャットで同意した
    ものは "llm"。**ユーザーの同意なしに呼ばない。**
    不正な要素が 1 つでもあれば何も書かずにエラーを返す（同じ spot_id を 1 回の送信で
    2 回以上名指しした場合も含む——どちらを採るか決められないため）。
    """
    if source not in ("user", "llm"):
        return _error("source は 'user' か 'llm' にしてください。")
    if submission_text:
        try:
            parsed = flag_log.parse_submission_text(submission_text)
        except ValueError as exc:
            return _error(str(exc))
        review_id, entries = parsed["review_id"], parsed["flags"]
    else:
        entries = flags
    if not review_id:
        return _error("review_id がありません。")
    directory = _review_dir(review_id)
    try:
        saved = review.load_review(directory, review_id) if directory else None
    except FileNotFoundError:
        saved = None
    if saved is None:
        return _error(f"review_id={review_id} のレビューが見つかりません。curation_review をやり直してください。")
    try:
        cleaned = flag_log.validate_entries(
            entries, allowed_spot_ids={s["spot_id"] for s in saved["spots"]})
    except ValueError as exc:
        return _error(str(exc))
    seen_spot_ids = set()
    duplicated = sorted({e["spot_id"] for e in cleaned if e["spot_id"] in seen_spot_ids
                         or seen_spot_ids.add(e["spot_id"])})
    if duplicated:
        return _error(f"同じ spot_id を 1 回の送信で複数回指定しています: {duplicated}")
    current = flag_log.alignment_key(saved["arf2_path"])
    if current["alignment_sha256"] != saved["alignment"]["alignment_sha256"]:
        return _error("レビューの後でアラインメント（.arf2）が変わっています。curation_review をやり直してください。")
    store = flag_log.FlagStore(flag_log.curation_dir(saved["arf2_path"]))
    n = store.append(cleaned, alignment=current, review_id=review_id, source=source)
    effective = store.effective(current["alignment_sha256"])
    return json_payload({"status": "ok", "recorded": n, "review_id": review_id,
                         "n_wrong": sum(1 for r in effective.values() if r["flag"] == "wrong"),
                         "n_suspect": sum(1 for r in effective.values() if r["flag"] == "suspect")})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def curation_flags(file_path: str | None = None) -> str:
    """現在のアラインメントで有効なフラグ（スポットごとの最新 1 行、clear 済みは除く）を TSV で返す。"""
    arf2_path = resolve_arf2_file_path(file_path)
    if not arf2_path:
        return mcp_errors.missing_state("arf2_file", ["load_dataset", "arf2_parser"],
                                        ".arf2 が見つかりません。先に load_dataset を実行してください。")
    key = flag_log.alignment_key(arf2_path)
    effective = flag_log.FlagStore(flag_log.curation_dir(arf2_path)).effective(key["alignment_sha256"])
    lines = ["spot_id\tflag\tnote\tsource\tts"]
    for spot_id in sorted(effective):
        row = effective[spot_id]
        lines.append("\t".join([str(spot_id), row["flag"], (row.get("note") or "").replace("\t", " "),
                                row.get("source", ""), row.get("ts", "")]))
    return json_payload({"alignment_file": key["alignment_file"], "n_flags": len(effective),
                         "table": "\n".join(lines)})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False, meta=_APP_ONLY_META)
def curation_view_data(review_id: str, page: int = 0) -> str:
    """ビューア（MCP Apps）専用。保存済みレビューのスポットをページ単位で返す。LLM は呼ばない。"""
    directory = _review_dir(review_id)
    try:
        saved = review.load_review(directory, review_id) if directory else None
    except FileNotFoundError:
        saved = None
    if saved is None:
        return _error(f"review_id={review_id} のレビューが見つかりません。")
    chunk = review.page(saved, page)
    head = {k: v for k, v in saved.items() if k != "spots"}
    return json_payload({"review": head, **chunk})
