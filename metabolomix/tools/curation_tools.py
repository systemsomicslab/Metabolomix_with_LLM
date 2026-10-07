"""アラインメントのキュレーション: curation_review / curation_suggest / curation_submit /
curation_flags / curation_view_data と、MCP Apps 用の ui:// ビューア。spec 2026-09-28。

deps: mcp_core / session_state / mcp_errors / serialization / path_resolvers / curation.*。
metabolomix.tools.* の他モジュールと server は import しない。
"""
from __future__ import annotations

from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.arf2.reader import load_catalog
from metabolomix.core import mcp_errors, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.path_resolvers import resolve_arf2_file_path
from metabolomix.core.serialization import json_payload, round_floats
# モジュール名を flag_log にするのは、curation_submit の引数 `flags`（公開 API の名前）が
# モジュールを隠すため。
from metabolomix.curation import evidence, judge, msdial_writeback, review, suggest, viewer
from metabolomix.curation import flags as flag_log
from metabolomix.library.defaults import DEFAULT_MS2_TOL, pick_tol

__all__ = ["curation_review", "curation_suggest", "curation_submit", "curation_flags",
           "curation_view_data"]

VIEWER_URI = "ui://ms-data-parser/curation-viewer"
_UI_META = {"ui": {"resourceUri": VIEWER_URI}, "ui/resourceUri": VIEWER_URI}
_APP_ONLY_META = {"ui": {"resourceUri": VIEWER_URI, "visibility": ["app"]}}
_LOCAL_WRITE_APPEND = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False)
# curation_submit は MS-DIAL の _tags.xml を書き換える（clear で Misannotation を外す）ので destructive。
_SUBMIT_WRITE = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=False)


def _error(message: str, **details) -> str:
    return json_payload({"status": "error", "message": message, **details})


def _flag_file_error(exc: flag_log.FlagFileError) -> str:
    return _error(str(exc), **exc.details())


def _bad_review_id(review_id) -> str | None:
    if review.is_valid_review_id(review_id) or suggest.is_valid_suggestion_id(review_id):
        return None
    return _error(f"review_id={review_id!r} の形が不正です（cr-… か cs-…）。"
                  "ビューアの送信用テキストをそのまま貼ってください。")


def _candidate_dirs(review_id: str, *arf2_paths) -> list[Path]:
    """レビューを探すフォルダの候補を優先順に: セッションの記録 → 渡された `.arf2`
    （送信用テキストの `arf2_path`、引数 `file_path` の順）→ 既定の解決。"""
    candidates = []
    known = session_state.session.curation.review_dirs.get(review_id)
    if known:
        candidates.append(Path(known))
    for arf2 in arf2_paths:
        if arf2:
            path = Path(arf2)
            candidates.append(path if path.is_dir() else flag_log.curation_dir(path))
    default = resolve_arf2_file_path(None)
    if default:
        candidates.append(flag_log.curation_dir(default))
    unique = []
    for path in candidates:
        if path not in unique:
            unique.append(path)
    return unique


def _load_any(directory, review_id: str) -> dict:
    if suggest.is_valid_suggestion_id(review_id):
        return suggest.load_suggestion(directory, review_id)
    return review.load_review(directory, review_id)


def _find_review(review_id: str, *arf2_paths) -> tuple[dict | None, list[Path]]:
    searched = _candidate_dirs(review_id, *arf2_paths)
    for directory in searched:
        try:
            return _load_any(directory, review_id), searched
        except FileNotFoundError:
            continue
    return None, searched


def _not_found(review_id: str, searched: list[Path]) -> str:
    where = "、".join(str(p) for p in searched) or "（候補なし: .arf2 を特定できませんでした）"
    return _error(f"review_id={review_id} のレビューが見つかりません。探したフォルダ: {where}。"
                  "ビューアの「送信用テキストをコピー」で作った文をもう一度そのまま貼るか、"
                  "レビューを作った .arf2 を file_path で指定してください。",
                  searched=[str(p) for p in searched])


@mcp.resource(VIEWER_URI, name="curation-viewer", mime_type="text/html;profile=mcp-app",
              description="キュレーションの一覧ビューア（MCP Apps）。データは curation_view_data から取る。")
def curation_viewer_resource() -> str:
    return viewer.render_html(None)


@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False, meta=_UI_META)
def curation_review(ontology: list[str] | None = None, name_contains: str | None = None,
                    file_ids: list[int] | None = None, max_traces: int = 12,
                    thresholds: dict | None = None, file_path: str | None = None,
                    max_rows: int = 100) -> str:
    """注釈付きスポットを一覧で確かめるレビューを作る（EIC・対向プロット・Δppm・ΔRT・
    RT–m/z 傾向と機械判別）。**先に library_load**（アラインメントに使われた
    `*_Loaded.msp2.dbs` を推奨）。

    対象: 既定は注釈付き全部。`ontology=["PG"]` でクラス、`name_contains` で名前の部分一致。
    戻り値は suspect 以上(とフラグ済み)のスポットだけの TSV（判定の重い順に先頭 `max_rows`
    行、既定 100。総数は `n_table_rows_total`、載せた数は `n_table_rows_shown`。全件は
    HTML ビューアにある）、判定の件数、クラス別の傾向要約（点数・R²・外れ数）、HTML
    ビューアのパス。`file_ids` は `.arf` の行にある試料 ID だけ（無い ID はエラー）。
    **EIC 系列やスペクトルは返さない**——ユーザーには
    `html_path` をブラウザで開いてもらい、ビューアで付けたフラグを「送信用テキストを
    コピー」→ チャットに貼ってもらう。貼られたら curation_submit(submission_text=...) に渡す。

    判定: likely_wrong（強い不一致: polarity_mismatch / precursor_unmatched / dmz_out /
    class_rule_rejected（MS-DIAL の脂質クラス規則による棄却で、脂質規則が走ったデータに限る）/
    adduct_isomer_of:<spot>（m/z が同時溶出する証拠の強い別物質スポットの別アダクトで説明できる。
    相手はレビューの対象に絞らずアラインメント全体から探す））、
    suspect（ppm_out / low_score / drt_out / eic_poor / adduct_isomer_minor_of:<spot>、
    または弱い兆候の重なり）、ok。
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
    try:
        result = review.run_review(arf2_path, spots, store=store, ms2_tol=ms2_tol, th=th,
                                   file_ids=file_ids, max_traces=max_traces,
                                   selection={"ontology": ontology, "name_contains": name_contains})
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    except evidence.UnknownFileIdsError as exc:
        return _error(str(exc), missing_file_ids=exc.missing)
    saved = review.save_review(result)
    session_state.session.curation.last_review_id = result["review_id"]
    session_state.session.curation.review_dirs[result["review_id"]] = str(saved["json"].parent)

    n_total = review.n_summary_rows(result)
    table = review.summary_tsv(result, max_rows=max_rows)
    return json_payload(round_floats({
        "review_id": result["review_id"],
        "n_spots": len(result["spots"]),
        "counts": result["counts"],
        "warnings": result["warnings"],
        "trend": review.trend_summary(result),
        "table": table,
        "n_table_rows_total": n_total,
        "n_table_rows_shown": len(table.splitlines()) - 1,
        "table_note": ("table は判定の重い順の先頭だけです。全件は html_path のビューアにあります。"
                       if len(table.splitlines()) - 1 < n_total else
                       "table は該当する全件です。カードは html_path のビューアで見られます。"),
        "html_path": str(saved["html"]),
        "thresholds": th, "ms2_tol": ms2_tol,
    }, 4))


@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False)
def curation_suggest(review_id: str | None = None, wrong: str = "flagged_or_likely",
                     unannotated: bool = True, include_decided: bool = False, top_n: int = 5,
                     rt_window: float | None = None, relation_mz_tol: float | None = None,
                     relation_min_r: float | None = None, thresholds: dict | None = None,
                     file_path: str | None = None, max_rows: int = 100) -> str:
    """キュレーションで「間違い」になったスポットと未注釈スポットに、注釈の候補を並べる
    （候補付けレビュー）。**先に library_load と curation_review**。

    対象: `wrong="flagged_or_likely"`（既定）= wrong フラグ ＋ 元レビューの likely_wrong、
    `"flagged"` = wrong フラグだけ。`unannotated=True` で未注釈スポットも。判断済み
    （assign / redundant）のスポットは `include_decided=True` のときだけ含める。
    元レビューは `review_id`（省略時はこのアラインメントの最新のレビュー）。
    候補: ① MS-DIAL の下位候補、② 閾値を緩めた再検索（スコアの足切りなし）、
    ④ 注釈付きの別スポットの同位体・アダクト・インソース断片としての説明。
    スペクトル類似度で並べ、極性矛盾・|Δm/z| ≥ 10 mDa は削り、非典型アダクト・RT–m/z 傾向の外れ・
    一致ピーク 0 は順位を下げる。MS-DIAL の脂質規則は評価していない（鎖組成は保証しない）ので、
    既定では和組成で記録する。
    戻り値は 1 スポット 1 行の TSV（強い説明のあるものが先、先頭 `max_rows` 行）・件数・`html_path`。
    ユーザーには `html_path` をブラウザで開いてもらい、選んだ内容を「送信用テキストをコピー」で
    チャットに貼ってもらう。貼られたら curation_submit(submission_text=...) に渡す。
    **ユーザーの同意なしに curation_submit を呼ばない。**
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
    if review_id is not None:
        if not review.is_valid_review_id(review_id):
            return _error(f"review_id={review_id!r} の形が不正です（cr-YYYYMMDD-HHMMSS-xxxx）。")
        base, searched = _find_review(review_id, arf2_path)
        if base is None:
            return _not_found(review_id, searched)
    else:
        base = suggest.latest_review(arf2_path, flag_log.alignment_key(arf2_path)["alignment_sha256"])
        if base is None:
            return mcp_errors.missing_state(
                "curation_review", ["curation_review"],
                "このアラインメントのレビューがありません。先に curation_review を実行してください"
                "（likely_wrong の判定をそこから読みます）。")
    options = {"wrong": wrong, "unannotated": unannotated, "include_decided": include_decided,
               "top_n": top_n, "rt_window": rt_window, "relation_mz_tol": relation_mz_tol,
               "relation_min_r": relation_min_r}
    try:
        result = suggest.run_suggestion(arf2_path, base_review=base, store=store, th=th, options=options)
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    except ValueError as exc:
        return _error(str(exc))
    saved = suggest.save_suggestion(result)
    session_state.session.curation.review_dirs[result["suggestion_id"]] = str(saved["json"].parent)
    table = suggest.summary_tsv(result, max_rows=max_rows)
    return json_payload(round_floats({
        "suggestion_id": result["suggestion_id"], "base_review_id": result["base_review_id"],
        "n_spots": len(result["spots"]), "counts": result["counts"], "warnings": result["warnings"],
        "table": table, "n_table_rows_shown": len(table.splitlines()) - 1,
        "html_path": str(saved["html"]),
        "library": {**result["library"], "path": session_state.session.library.source_path},
        "analysis_params": {k: result["analysis_params"][k] for k in ("source", "path", "rt_window")},
        "options": result["options"], "thresholds": th}, 4))


@mcp.tool(annotations=_SUBMIT_WRITE, structured_output=False)
def curation_submit(submission_text: str | None = None, review_id: str | None = None,
                    flags: list[dict] | None = None, source: str = "user",
                    file_path: str | None = None) -> str:
    """キュレーションのフラグを記録する（追記。フラグの無いスポットは「間違っていない」で何も書かない）。

    ユーザーがビューアの「送信用テキストをコピー」で貼った文をそのまま `submission_text` に渡す
    （`CURATION_SUBMIT ` で始まる行だけを読む。書き写さないこと）。直接渡すなら `review_id` と
    `flags=[{"spot_id": 12, "flag": "wrong" | "suspect" | "clear", "note": "..."}]`。
    候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign|redundant|clear, candidate: "L1"|"R1",
    level: sum|species, note}]`。候補の中身（名前・InChIKey）は保存済みの候補付けから展開するので、
    送信側で名前を書かない。assign / redundant は `_tags.xml` を変えない。
    `source` は誰の判断か: ユーザー自身の判断は "user"、LLM の提案にユーザーがチャットで同意した
    ものは "llm"。**ユーザーの同意なしに呼ばない。**
    不正な要素が 1 つでもあれば何も書かずにエラーを返す（同じ spot_id を 1 回の送信で
    2 回以上名指しした場合も含む——どちらを採るか決められないため）。
    レビューはセッションの記録 → 送信用テキストの `arf2_path` → `file_path`（レビューを
    作った `.arf2`）→ 既定の `.arf2` の順に探すので、サーバ再起動の後でも貼った文で送れる。
    記録の後、アラインメントの `_tags.xml` の Misannotation に反映する（wrong → 付ける、
    clear → 外す、suspect → 触らない。控えは `curation/tags-backup/`）。結果は `tags_xml`。
    反映に失敗しても記録は残る。MS-DIAL でプロジェクトを開いたままだと GUI の保存で
    上書きされるので、`tags_xml.note` をユーザーに伝える。
    """
    if source not in ("user", "llm"):
        return _error("source は 'user' か 'llm' にしてください。")
    text_arf2 = None
    if submission_text:
        try:
            parsed = flag_log.parse_submission_text(submission_text)
        except ValueError as exc:
            return _error(str(exc))
        review_id, entries, text_arf2 = parsed["review_id"], parsed["flags"], parsed["arf2_path"]
    else:
        entries = flags
    if not review_id:
        return _error("review_id がありません。")
    bad = _bad_review_id(review_id)
    if bad:
        return bad
    saved, searched = _find_review(review_id, text_arf2, file_path)
    if saved is None:
        return _not_found(review_id, searched)
    try:
        if suggest.is_valid_suggestion_id(review_id):
            cleaned = suggest.expand_entries(entries, saved)
        else:
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
        return _error("レビューの後でアラインメント（.arf2）が変わっています。curation_review（候補付けなら curation_suggest）をやり直してください。")
    store = flag_log.FlagStore(flag_log.curation_dir(saved["arf2_path"]))
    try:
        store.rows()                       # 壊れた記録に追記しない（先に読めるか確かめる）
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    n = store.append(cleaned, alignment=current, review_id=review_id, source=source)
    effective = store.effective(current["alignment_sha256"])
    tags_xml = msdial_writeback.sync_misannotation(saved["arf2_path"], cleaned)
    return json_payload({"status": "ok", "recorded": n, "review_id": review_id,
                         "n_wrong": sum(1 for r in effective.values() if r["flag"] == "wrong"),
                         "n_suspect": sum(1 for r in effective.values() if r["flag"] == "suspect"),
                         "n_assign": sum(1 for r in effective.values() if r["flag"] == "assign"),
                         "n_redundant": sum(1 for r in effective.values() if r["flag"] == "redundant"),
                         "tags_xml": tags_xml})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def curation_flags(file_path: str | None = None) -> str:
    """現在のアラインメントで有効な判断（スポットごとの最新 1 行、clear 済みは除く）を TSV で返す。
    列は spot_id・flag・name（assign の記録名）・of（redundant の相手スポット）・note・source・ts。"""
    arf2_path = resolve_arf2_file_path(file_path)
    if not arf2_path:
        return mcp_errors.missing_state("arf2_file", ["load_dataset", "arf2_parser"],
                                        ".arf2 が見つかりません。先に load_dataset を実行してください。")
    key = flag_log.alignment_key(arf2_path)
    try:
        effective = flag_log.FlagStore(flag_log.curation_dir(arf2_path)).effective(
            key["alignment_sha256"])
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    lines = ["spot_id\tflag\tname\tof\tnote\tsource\tts"]
    for spot_id in sorted(effective):
        row = effective[spot_id]
        lines.append("\t".join([str(spot_id), row["flag"], str(row.get("name") or ""),
                                "" if row.get("of") is None else str(row["of"]),
                                (row.get("note") or "").replace("\t", " "),
                                row.get("source", ""), row.get("ts", "")]))
    return json_payload({"alignment_file": key["alignment_file"], "n_flags": len(effective),
                         "table": "\n".join(lines)})


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False, meta=_APP_ONLY_META)
def curation_view_data(review_id: str, page: int = 0) -> str:
    """ビューア（MCP Apps）専用。保存済みレビューのスポットをページ単位で返す。LLM は呼ばない。"""
    if suggest.is_valid_suggestion_id(review_id):
        return _error("curation_view_data はレビュー（cr-…）専用です。")
    bad = _bad_review_id(review_id)
    if bad:
        return bad
    saved, searched = _find_review(review_id)
    if saved is None:
        return _error(f"review_id={review_id} のレビューが見つかりません。",
                      searched=[str(p) for p in searched])
    chunk = review.page(saved, page)
    head = {k: v for k, v in saved.items() if k != "spots"}
    return json_payload({"review": head, **chunk})
