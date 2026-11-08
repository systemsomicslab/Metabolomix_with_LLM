"""アラインメントのキュレーション: curation_review / curation_suggest / curation_submit /
curation_flags / curation_view_data と、MCP Apps 用の ui:// ビューア。spec 2026-09-28。

キュレーションは `DatasetState`（ハンドル）の上で動く。`.arf2` はデータセットが覚えているもの
（`dataset.artifacts.artifact_path(ds, "arf2")`）で、`data/` は探さない。ライブラリは
データセットのもの（`holder.library`）。`curation_submit` が書いたら、同じ `.arf2` を開いている
ハンドルすべての `curation_flags` を読み直す（`submission.after_write`）。

deps: mcp_core / session_state / handles / mcp_errors / serialization / dataset.* / curation.*。
metabolomix.tools.* の `_holder` 以外と server は import しない。
"""
from __future__ import annotations

from pathlib import Path, PureWindowsPath

from mcp.types import ToolAnnotations

from metabolomix.arf2.reader import load_catalog
from metabolomix.core import file_digest, handles, mcp_errors, session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.dataset.artifacts import artifact_path
from metabolomix.dataset.msdial_adapter import reload_curation
from metabolomix.tools._holder import resolve_holder
# モジュール名を flag_log にするのは、curation_submit の引数 `flags`（公開 API の名前）が
# モジュールを隠すため。
from metabolomix.curation import evidence, judge, review, submission, submit_server, suggest, viewer
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

#: 他のハンドルのロックを待つ上限（秒）。自分のハンドルのロックを持ったまま別のハンドルを待つので、
#: 無限に待つと 2 つのハンドルから互いに提出したときにロックの順で止まる。
REFRESH_LOCK_TIMEOUT = 5.0


def _file_name(path) -> str | None:
    """ライブラリの置き場所は戻り値に出さない（ファイル名だけ。区切りは `/` でも バックスラッシュでも）。"""
    return PureWindowsPath(str(path)).name if path else None


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
    （送信用テキストの `arf2_path`、データセットの `.arf2` の順）。`data/` は探さない。"""
    candidates = []
    known = session_state.session.curation.review_dirs.get(review_id)
    if known:
        candidates.append(Path(known))
    for arf2 in arf2_paths:
        if arf2:
            path = Path(arf2)
            candidates.append(path if path.is_dir() else flag_log.curation_dir(path))
    unique = []
    for path in candidates:
        if path not in unique:
            unique.append(path)
    return unique


def _load_any(directory, review_id: str) -> dict:
    return submission.load_saved(directory, review_id)


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
                  "レビューを作ったデータセットを dataset_handle で指定してください。",
                  searched=[str(p) for p in searched])


def _no_arf2(ds) -> str:
    return _error(
        "このデータセットには .arf2 が無いので、キュレーションできません（キュレーションは MS-DIAL の "
        "アラインメントの .arf2 に対して行う）。MS-DIAL の解析フォルダを load_dataset(directory=...) で"
        "開き直すか、load_dataset の arf2_path で .arf2 を指定してください。",
        source_format=getattr(ds, "source_format", None))


def _dataset_arf2(holder):
    """データセットの `.arf2`（`load_dataset(arf2_path=...)` の明示 → 記録の `.arf2`。無ければ None）。"""
    ds = holder.dataset
    if ds is None:
        return None
    return holder.arf2_path or artifact_path(ds, "arf2")


def _dataset_context(dataset_handle):
    """`(holder, arf2_path)`、足りなければ `(封筒, None)`（str）。

    ハンドルが無い・閉じていれば `missing_state("dataset", ["load_dataset"])` か
    `DATASET_HANDLE_NOT_FOUND`（別のハンドルや `data/` へ落ちない）。`.arf2` を持たない
    データセット（mzTab-M だけ）は理由つきのエラー。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder, None
    ds = holder.dataset
    if ds is None:
        return mcp_errors.missing_state(
            "dataset", ["load_dataset"],
            "データセットが開いていません。先に load_dataset で MS-DIAL の解析フォルダを開いてください。"), None
    arf2 = _dataset_arf2(holder)
    if not arf2:
        return _no_arf2(ds), None
    return holder, Path(arf2)


def _library_store(holder):
    """データセットのライブラリの store か、`library_load` を促す封筒（str）。"""
    store = holder.library.store
    if store is None:
        return None, mcp_errors.missing_state(
            "library", ["library_load"],
            "このデータセットに参照ライブラリが読み込まれていません。先に library_load を実行してください"
            "（アラインメントと同じフォルダの *_Loaded.msp2.dbs を推奨）。")
    return store, None


def _refresh_open_handles(arf2_path) -> dict | None:
    """`arf2_path` を開いているハンドルすべての `curation_flags` を読み直す（`submission.after_write`）。

    同じ `.arf2` は sha256 で照合する（ハンドルの `.arf2` の記録。`handles.handles_for_arf2`）。
    ハンドルごとに**そのハンドルのロックの下で**、**1 つずつ**（入れ子にしない）読み直す。
    呼び出し元の `curation_submit` が自分のハンドルのロックを持っていても、RLock なので
    そのハンドルは待たずに通る。別のスレッドが持っているハンドルは `REFRESH_LOCK_TIMEOUT` 秒
    だけ待ち、取れなければ読み直さず `not_refreshed` に挙げる（ロックの順で止まらないため）。
    """
    try:
        sha = file_digest.sha256_of(arf2_path)
    except OSError:
        return None
    targets = handles.registry.handles_for_arf2(sha)
    if not targets:
        return None
    refreshed, skipped = [], []
    for handle in targets:
        if not handle.lock.acquire(timeout=REFRESH_LOCK_TIMEOUT):
            skipped.append(handle.handle_id)
            continue
        try:
            if handles.registry.peek(handle.handle_id) is not handle:
                continue                                     # 待つ間に閉じられた
            reload_curation(handle.dataset, handle.arf2_path)
            refreshed.append(handle.handle_id)
        finally:
            handle.lock.release()
    return {"refreshed": refreshed, "not_refreshed": skipped}


submission.after_write = _refresh_open_handles


@mcp.resource(VIEWER_URI, name="curation-viewer", mime_type="text/html;profile=mcp-app",
              description="キュレーションの一覧ビューア（MCP Apps）。データは curation_view_data から取る。")
def curation_viewer_resource() -> str:
    return viewer.render_html(None)


def _open_submit(kind: str, review_id: str, arf2_path) -> tuple[dict | None, dict]:
    """ビューアの受け口に登録する。立てられなければ Copy の経路だけ（記録は curation_submit）。"""
    try:
        endpoint = submit_server.register(kind, review_id, arf2_path)
    except OSError as exc:
        return None, {"via": "copy", "reason": f"受け口を立てられませんでした: {exc}"}
    return endpoint, {"via": "viewer", "idle_timeout_min": endpoint["idle_timeout_min"]}


@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False, meta=_UI_META)
def curation_review(ontology: list[str] | None = None, name_contains: str | None = None,
                    file_ids: list[int] | None = None, max_traces: int = 12,
                    thresholds: dict | None = None, max_rows: int = 100,
                    dataset_handle: str | None = None) -> str:
    """注釈付きスポットを一覧で確かめるレビューを作る（EIC・対向プロット・Δppm・ΔRT・
    RT–m/z 傾向と機械判別）。**先に load_dataset と library_load**（アラインメントに使われた
    `*_Loaded.msp2.dbs` を推奨）。

    対象のデータセットは `dataset_handle`（省けば直近）。`.arf2` と、EIC・`.dcl` はそのデータセットが
    覚えているものを使い、`data/` は探さない。ライブラリもそのデータセットのもの。`.arf2` の無い
    データセット（mzTab-M だけ）は理由つきのエラー。

    対象: 既定は注釈付き全部。`ontology=["PG"]` でクラス、`name_contains` で名前の部分一致。
    戻り値は suspect 以上(とフラグ済み)のスポットだけの TSV（判定の重い順に先頭 `max_rows`
    行、既定 100。総数は `n_table_rows_total`、載せた数は `n_table_rows_shown`。全件は
    HTML ビューアにある）、判定の件数、クラス別の傾向要約（点数・R²・外れ数）、HTML
    ビューアのパス。`file_ids` は `.arf` の行にある試料 ID だけ（無い ID はエラー）。
    **EIC 系列やスペクトルは返さない**——ユーザーには
    `html_path` をブラウザで開いてもらい、ビューアの **Submit** で送ってもらう（記録と
    `_tags.xml` の反映はビューアが直接行う。`submit.via == "viewer"`）。送信用テキストを貼るよう
    頼まない。送った結果は curation_flags で確かめられる。Copy のテキストが貼られたとき
    （受け口の時間切れ・`submit.via == "copy"`）は従来どおり curation_submit(submission_text=...) に渡す。

    判定: likely_wrong（強い不一致: polarity_mismatch / precursor_unmatched / dmz_out /
    class_rule_rejected（MS-DIAL の脂質クラス規則による棄却で、脂質規則が走ったデータに限る）/
    adduct_isomer_of:<spot>（m/z が同時溶出する証拠の強い別物質スポットの別アダクトで説明できる。
    相手はレビューの対象に絞らずアラインメント全体から探す））、
    suspect（ppm_out / low_score / drt_out / eic_poor / adduct_isomer_minor_of:<spot>、
    または弱い兆候の重なり）、ok。
    RT–m/z 傾向（trend_outlier）は補強にしか使わない。UNKNOWN は不一致に数えない。
    `thresholds` で既定のしきい値（ppm_pass=5, ppm_borderline=10, drt_pass=0.5 分 など）を上書きできる。
    """
    holder, arf2_path = _dataset_context(dataset_handle)
    if arf2_path is None:
        return holder
    store, missing = _library_store(holder)
    if missing:
        return missing
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
                                   selection={"ontology": ontology, "name_contains": name_contains},
                                   files=evidence.dataset_files(holder.dataset))
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    except evidence.UnknownFileIdsError as exc:
        return _error(str(exc), missing_file_ids=exc.missing)
    endpoint, submit_info = _open_submit("review", result["review_id"], result["arf2_path"])
    saved = review.save_review(result, submit_endpoint=endpoint)
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
        "submit": submit_info,
        "thresholds": th, "ms2_tol": ms2_tol,
    }, 4))


@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False)
def curation_suggest(review_id: str | None = None, wrong: str = "flagged_or_likely",
                     unannotated: bool = True, include_decided: bool = False, top_n: int = 5,
                     rt_window: float | None = None, relation_mz_tol: float | None = None,
                     relation_min_r: float | None = None, thresholds: dict | None = None,
                     max_rows: int = 100, dataset_handle: str | None = None) -> str:
    """キュレーションで「間違い」になったスポットと未注釈スポットに、注釈の候補を並べる
    （候補付けレビュー）。**先に load_dataset・library_load・curation_review**。
    対象のデータセットは `dataset_handle`（省けば直近。`curation_review` と同じ規則）。

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
    ユーザーには `html_path` をブラウザで開いてもらい、選んだ内容をビューアの **送信** で送ってもらう
    （記録と `_tags.xml` の反映はビューアが直接行う。`submit.via == "viewer"`）。送信用テキストを
    貼るよう頼まない。Copy のテキストが貼られたとき（受け口の時間切れ・`submit.via == "copy"`）は
    curation_submit(submission_text=...) に渡す。
    **ユーザーの同意なしに curation_submit を呼ばない。**
    """
    holder, arf2_path = _dataset_context(dataset_handle)
    if arf2_path is None:
        return holder
    store, missing = _library_store(holder)
    if missing:
        return missing
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
        result = suggest.run_suggestion(arf2_path, base_review=base, store=store, th=th, options=options,
                                        files=evidence.dataset_files(holder.dataset))
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    except ValueError as exc:
        return _error(str(exc))
    endpoint, submit_info = _open_submit("suggest", result["suggestion_id"], result["arf2_path"])
    saved = suggest.save_suggestion(result, submit_endpoint=endpoint)
    session_state.session.curation.review_dirs[result["suggestion_id"]] = str(saved["json"].parent)
    table = suggest.summary_tsv(result, max_rows=max_rows)
    return json_payload(round_floats({
        "suggestion_id": result["suggestion_id"], "base_review_id": result["base_review_id"],
        "n_spots": len(result["spots"]), "counts": result["counts"], "warnings": result["warnings"],
        "table": table, "n_table_rows_shown": len(table.splitlines()) - 1,
        "html_path": str(saved["html"]),
        "submit": submit_info,
        "library": {**result["library"], "path": _file_name(holder.library.source_path)},
        "analysis_params": {k: result["analysis_params"][k] for k in ("source", "path", "rt_window")},
        "options": result["options"], "thresholds": th}, 4))


@mcp.tool(annotations=_SUBMIT_WRITE, structured_output=False)
def curation_submit(submission_text: str | None = None, review_id: str | None = None,
                    flags: list[dict] | None = None, source: str = "user",
                    dataset_handle: str | None = None) -> str:
    """キュレーションのフラグを記録する（追記。フラグの無いスポットは「間違っていない」で何も書かない）。

    ユーザーがビューアの「送信用テキストをコピー」で貼った文をそのまま `submission_text` に渡す
    （`CURATION_SUBMIT ` で始まる行だけを読む。書き写さないこと）。直接渡すなら `review_id` と
    `flags=[{"spot_id": 12, "flag": "wrong" | "suspect" | "confirmed" | "clear", "note": "..."}]`
    （confirmed = 注釈が正しいと確かめた。wrong / suspect と排他）。
    候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign|redundant|clear, candidate: "L1"|"R1",
    level: sum|species, note}]`。候補の中身（名前・InChIKey）は保存済みの候補付けから展開するので、
    送信側で名前を書かない。assign / redundant は `_tags.xml` を変えない。
    `source` は誰の判断か: ユーザー自身の判断は "user"、LLM の提案にユーザーがチャットで同意した
    ものは "llm"。**ユーザーの同意なしに呼ばない。**
    不正な要素が 1 つでもあれば何も書かずにエラーを返す（同じ spot_id を 1 回の送信で
    2 回以上名指しした場合も含む——どちらを採るか決められないため）。
    レビューはセッションの記録 → 送信用テキストの `arf2_path` → データセット（`dataset_handle`、
    省けば直近）の `.arf2` のフォルダの順に探す（`data/` は探さない）ので、サーバ再起動の後でも、
    貼った文なら送れる。データセットも貼った文の `arf2_path` も無ければ `missing_state`。
    書いた後、同じ `.arf2`（sha256 で照合）を開いているハンドルすべての `curation_flags` を読み直す
    （結果の `open_handles`: `refreshed` と、別の呼び出しがロックを持っていて読み直せなかった
    `not_refreshed`）。書き出しの `apply_curation` は新しい判断をすぐ使う。
    記録の後、アラインメントの `_tags.xml` に反映する（wrong → Misannotation を付け Confirmed を外す、
    confirmed → Confirmed を付け Misannotation を外す、suspect → Confirmed を外す、clear → 両方外す。
    控えは `curation/tags-backup/`）。結果は `tags_xml`（`added` / `removed` は Misannotation、
    `confirmed` は Confirmed の変化）。
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
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    saved, searched = _find_review(review_id, text_arf2, _dataset_arf2(holder))
    if saved is None and not searched:
        return mcp_errors.missing_state(
            "dataset", ["load_dataset"],
            "レビューのあるデータセットが分かりません。先に load_dataset でレビューを作った解析フォルダを"
            "開いてください（送信用テキストの arf2_path があればそれも使えます）。")
    if saved is None:
        return _not_found(review_id, searched)
    try:
        result = submission.submit_flags(saved, entries, review_id=review_id, source=source)
    except submission.SubmissionError as exc:
        return _error(str(exc), **exc.details)
    return json_payload(result)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def curation_flags(dataset_handle: str | None = None) -> str:
    """データセット（`dataset_handle`、省けば直近）の `.arf2` で有効な判断（スポットごとの最新 1 行、
    clear 済みは除く）を TSV で返す。
    列は spot_id・flag・name（assign の記録名）・of（redundant の相手スポット）・note・source・ts。"""
    holder, arf2_path = _dataset_context(dataset_handle)
    if arf2_path is None:
        return holder
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
def curation_view_data(review_id: str, page: int = 0, dataset_handle: str | None = None) -> str:
    """ビューア（MCP Apps）専用。保存済みレビューのスポットをページ単位で返す。LLM は呼ばない。
    レビューはセッションの記録か、データセット（`dataset_handle`、省けば直近）の `.arf2` のフォルダから探す。"""
    if suggest.is_valid_suggestion_id(review_id):
        return _error("curation_view_data はレビュー（cr-…）専用です。")
    bad = _bad_review_id(review_id)
    if bad:
        return bad
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    saved, searched = _find_review(review_id, _dataset_arf2(holder))
    if saved is None:
        return _error(f"review_id={review_id} のレビューが見つかりません。",
                      searched=[str(p) for p in searched])
    chunk = review.page(saved, page)
    head = {k: v for k, v in saved.items() if k != "spots"}
    return json_payload({"review": head, **chunk})
