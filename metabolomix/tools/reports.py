"""解析・解釈レポートと図の保存ツール群。

write_report / read_report（ID を省くと一覧）/ save_figure。
レポート先は `[paths] reports_dir`（あれば先）→ データセットの解析フォルダ配下 `reports/`（mcp_core）。
各ツールは `dataset_handle` で読み書きするハンドルを選ぶ（省略時は直近）。save_figure は `dataset_plot` が置いた
`ds.figures` から描き、ファイルへの書き方は `plot_tools.write_figure_files`（`dataset_plot(output_path=...)` と共有）。
deps: mcp_core / plot_tools（遅延 import）/ notes。
server は import しない。
"""

from metabolomix.tools import notes
from metabolomix.core import handles, mcp_errors
from metabolomix.core.mcp_errors import config_invalid_as_envelope
from mcp.types import ToolAnnotations
from metabolomix.core.serialization import json_payload
from metabolomix.core.mcp_core import mcp, _resolve_report_dir, _report_dir_candidates, _build_report_meta
from metabolomix.tools._holder import resolve_holder as _resolve_holder

__all__ = [
    "write_report",
    "read_report",
    "save_figure",
]


# --- 解析・解釈レポート（reports/<analysis_id>.md） ---
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
@config_invalid_as_envelope
def write_report(
    analysis_id: str,
    dataset: str,
    body: str,
    status: str = "draft",
    dataset_handle: str | None = None,
) -> str:
    """解析・解釈レポートを reports/<analysis_id>.md に上書き保存する（成果物＋記録）。

    body は frontmatter を含まない Markdown 本文（見出しは自由）。所見が増えたら本文を作り直して再度呼ぶ
    （ファイルは毎回上書き）。
    lipidmix.local.toml の [paths] reports_dir の明示先を優先し、書込不可なら解析フォルダ配下reports/。
    未指定時は解析フォルダ配下reports/だけを使う（リポジトリ側へは退避しない）。
    解析フォルダは直近のデータセット（または dataset_handle）のもの。設定の置き場所も
    データセットも無いときは missing_state（load_dataset）。

    analysis_id はファイル名 slug の元になるため ASCII で一意に。別IDが同一slugに潰れて既存
    レポートを上書きしそうな場合は保存を中止して通知する（objectiveレコードと同じIDを推奨）。
    """
    holder = _resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    slug = notes.make_slug(analysis_id)
    # slug衝突ガード: 別の analysis_id が同一ファイル名に潰れる場合は、既存レポートを
    # 黙って上書きせず中止する（同一IDの上書きは意図どおり許可）。非ASCII/記号違いのIDで起きうる。
    for directory in _report_dir_candidates(holder):
        existing = directory / f"{slug}.md"
        if existing.is_file():
            ex_meta, _ = notes.parse_frontmatter(existing.read_text(encoding="utf-8"))
            ex_id = str(ex_meta.get("analysis_id", ""))
            if ex_id and ex_id != analysis_id:
                return (
                    f"slug衝突のため中止: analysis_id '{analysis_id}' はファイル名 '{slug}.md' に潰れますが、"
                    f"そこには別の '{ex_id}' のレポートが既にあります（{existing}）。"
                    "上書きを避けました。ASCIIで一意な analysis_id を指定してください。"
                )
            break  # 同一ID → 上書きしてよい
    reports_dir = _resolve_report_dir(holder)
    meta = _build_report_meta(analysis_id, dataset, status)
    path = notes.write_note(reports_dir, slug, meta, body)
    return f"レポートを保存: {path}（status={status}）。read_report('{analysis_id}') で読み戻せます。"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
@config_invalid_as_envelope
def read_report(analysis_id: str | None = None, dataset_handle: str | None = None) -> str:
    """過去レポートを読み戻す。analysis_id を省くと既存レポートの 1 行索引を返す。

    dataset_handle: 解析フォルダを決めるハンドル（省略時は直近）。
    analysis_id を渡すと、候補ディレクトリ横断で最新更新のものの本文を返す（セッション継続用）。
    解析フォルダ配下と退避先の両方に同名レポートが残る場合（書き込み可否が途中で変化した等）、
    古い方を返さないよう mtime が最新のファイルを採用する。
    省いた場合の索引は 1 行ずつ `analysis_id / date / status`。候補ディレクトリ横断で同一
    analysis_id が重複するときは mtime が最新の 1 件を採用する。
    """
    holder = _resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    if analysis_id is None:
        return _list_reports(holder)
    slug = notes.make_slug(analysis_id)
    matches = [d / f"{slug}.md" for d in _report_dir_candidates(holder)]
    matches = [p for p in matches if p.is_file()]
    if not matches:
        return f"レポートが見つかりません: {analysis_id}（write_report で作成してください）"
    newest = max(matches, key=lambda p: p.stat().st_mtime)
    return newest.read_text(encoding="utf-8")


def _list_reports(holder) -> str:
    """既存レポートの 1 行索引（analysis_id / date / status）。"""
    # analysis_id -> (mtime, 表示行)。最新更新の行を残す。
    best: dict[str, tuple[float, str]] = {}
    for directory in _report_dir_candidates(holder):
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.md")):
            meta, _body = notes.parse_frontmatter(path.read_text(encoding="utf-8"))
            if meta.get("type") != "report":
                continue
            aid = str(meta.get("analysis_id", path.stem))
            line = f"- {aid} | date={meta.get('date', '?')} | status={meta.get('status', '?')}"
            mtime = path.stat().st_mtime
            if aid not in best or mtime > best[aid][0]:
                best[aid] = (mtime, line)
    lines = ["# レポート一覧"]
    for aid, (_mtime, line) in sorted(best.items()):
        lines.append(line)
    if len(lines) == 1:
        lines.append("（レポートはまだありません）")
    return "\n".join(lines)


#: save_figure が受け付ける図の種類（`dataset_plot` の kind）。
FIGURE_KINDS = ("pca", "volcano", "loadings", "group_intensity", "species", "species_pca",
                "eic", "eic_compounds")
#: `result_id` で登録された結果から描ける kind（`ds.results`）。
_RESULT_KINDS = ("pca", "volcano", "loadings")
#: 結果から作った図（保存の直前に、その結果が今のものか確かめる）。`species_pca` は項目の指定から作る
#: 結果なので `result_id` では描き直せない（`_RESULT_KINDS` には入れない）。`save_figure` は名指しされた
#: `result_id` と置き場の図の結果が同じかだけを確かめる（`_other_species_pca`）。
_RESULT_DERIVED_KINDS = _RESULT_KINDS + ("species_pca",)
_LABELS = {"pca": "PCA", "volcano": "volcano", "loadings": "ローディング",
           "group_intensity": "群別強度", "species": "分子種", "species_pca": "分子種の PCA",
           "eic": "EIC", "eic_compounds": "EIC 重ね描き"}


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
@config_invalid_as_envelope
def save_figure(kind: str, analysis_id: str, title: str | None = None,
                result_id: str | None = None, dataset_handle: str | None = None,
                contrast_id: str | None = None) -> str:
    """データセットの図を reports/figures/<analysis_id>_<kind>.png に保存する（報告書に埋め込む図）。

    kind は `dataset_plot` の kind: "pca" / "volcano"（間引き前の全点）/ "loadings" / "group_intensity" /
    "species" / "species_pca" / "eic" / "eic_compounds"。loadings / group_intensity / species / species_pca は
    dpi 300 の PNG と同名 .svg、pca / volcano / eic / eic_compounds は dpi 120 の PNG だけ。`dataset_plot` が最後に作った図（`ds.figures[kind]`）を描く。
    まだ一度も描いていなければ missing_state（`dataset_plot` を案内）。結果から作った図（pca / volcano / loadings /
    species_pca）は、その結果が現在のものでなければ保存せず `STALE_FIGURE`（理由と `required_tools=["dataset_plot"]`。
    前処理・解析をやり直してから描き直す）。除外（`dataset_exclude`）と試料の表の更新
    （`dataset_set_sample_metadata`。前処理の入力の列なら pca / loadings / species_pca / volcano、群・足した因子の列
    だけなら volcano）は図ごと捨てるが、捨てた図も `STALE_FIGURE`（`details.reason` に実際に変わったもの: `exclusion`
    や列名）で、
    描いたことの無い図（missing_state）と見分けられる。group_intensity / species / species_pca は群を図に焼き込むので、
    試料の情報（群・役割）が更新されたあとも `STALE_FIGURE`（species_pca は結果から作った図でもある）。
    eic / eic_compounds は結果から作った図ではないので、描いたとおりに保存する（あとで `dataset_exclude` しても変わらない）。
    result_id: kind が pca / volcano / loadings で効く。登録された結果を名指しして描く
      （`dataset_plot` を呼ばなくてよい。古い結果は STALE_ANALYSIS_RESULT）。volcano では
      `"result_id#contrast_id"`（`dataset_feature_search` が出す形）でも対比を選べる。species_pca では描き直しに
      使わず（項目の指定から作る結果）、置き場の図がその結果から描いたものかを確かめる: 同じなら保存し（同じ指定で
      描き直しただけの別の結果なら、置き場の図を保存して `notes:` 行にそう書く）、違えば保存せず `STALE_FIGURE`
      （`details.held_result_id` に置き場の図の結果、`details.mismatch` に違いの種類）。`settings`（項目の指定が違う）は
      `details.settings` の指定で `dataset_plot(kind="species_pca", ...)` を描き直してから保存する。`stale`（名指しした
      結果が除外・試料の表の更新より前のもの）と `input`（同じ指定だが、キュレーションの判断の更新などで入力の行列が
      違う）はその結果の図をもう描けないので、置き場の図を保存するなら `held_result_id` を渡すか省く（`stale` で指定も
      違えば `details.settings` も付く）。分子種の PCA の結果でなければ DATASET_BAD_REQUEST。そのハンドルの結果に
      無い ID も DATASET_BAD_REQUEST で、開いているほかのハンドルの分子種の PCA の結果なら `details.result_handle` に
      そのハンドル（その `dataset_handle` で呼び直す）。ほかのハンドルの別の種類の結果（比較・PCA・概観）なら呼び直しを
      案内せず、分子種の PCA の結果ではないと言ってそのハンドルを文で名指す（`details.result_kind` にその種類）。
      ほかの kind（group_intensity / species / eic / eic_compounds）では使わず、置き場の図を保存して戻り値の
      `notes:` 行にそう書く。
    contrast_id: kind が volcano で効く。`dataset_compare` の結果の対比（`T_vs_C` など）を選んで描き直す
      （`result_id` を省けば直近の差次的結果）。`result_id#contrast_id` と食い違えば DATASET_BAD_REQUEST。
      2 群の対比が複数ある結果で対比を選ばなければ `needs_input`（選択肢は対比の一覧）。ほかの kind では使わず、
      戻り値の `notes:` 行にそう書く。
    dataset_handle: `load_dataset` が返したハンドル（省略時は直近）。
    画面に出すだけなら `dataset_plot`（画像を返す。`output_path` で任意の場所にも書ける）。返り値の相対パスは
    write_report の本文に `![...](figures/<analysis_id>_<kind>.png)` として埋め込める。
    """
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.tools import plot_tools
    from metabolomix.tools._contrast_choice import contrast_envelope
    from metabolomix.tools._holder import resolve_holder

    if kind not in FIGURE_KINDS:
        return json_payload({"status": "error", "message":
                             f"kind は {', '.join(FIGURE_KINDS)} のどれかを指定してください（受け取った値: {kind!r}）。"})
    # 効かない引数は断らずに使わず、戻り値の notes に書く（その kind の図は引数によらず同じものを保存する）。
    notes_out = []
    # species_pca の result_id は描き直しには使えないが、置き場の図がその結果から描いたものかを確かめる
    # （別の項目の指定の図を、名指しした結果の図として保存しない）。
    species_pca_id = None
    if result_id is not None and kind == "species_pca":
        species_pca_id, result_id = result_id, None
    elif result_id is not None and kind not in _RESULT_KINDS:
        notes_out.append(f"result_id={result_id!r} は kind が {' / '.join(_RESULT_KINDS)} のときだけ効くので"
                         f"使っていない（kind={kind!r} は dataset_plot が置いた図を保存した）。")
        result_id = None
    if contrast_id is not None and kind != "volcano":
        notes_out.append(f"contrast_id={contrast_id!r} は kind が volcano のときだけ効くので使っていない"
                         f"（kind={kind!r}）。")
        contrast_id = None
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return mcp_errors.missing_state("dataset", ["load_dataset"],
                                        "DatasetState がありません。先に load_dataset を実行してください。")
    try:
        if result_id is not None or contrast_id is not None:
            record = plot_tools.build_figure(ds, kind, result_id=result_id, contrast_id=contrast_id)
        else:
            record = ds.figures.get(kind)
            if record is None:
                discarded = (getattr(ds, "discarded_figures", None) or {}).get(kind)
                if discarded is not None:
                    return mcp_errors.mztab_error("STALE_FIGURE", _discarded_message(kind, discarded), {
                        "kind": kind, "reason": list(discarded), "required_tools": ["dataset_plot"]})
                return mcp_errors.missing_state(
                    f"{kind}_figure", ["dataset_plot"],
                    f"先に dataset_plot(kind=\"{kind}\") で図を作ってください（このデータセットで {kind} の図は"
                    "まだ描かれていません）。")
            # 置き場の図そのものが古ければ、名指しした結果と照らす前に描き直しを求める（照らすのは今の図とだけ）。
            stale = _stale_figure(ds, record)
            if stale:
                return mcp_errors.mztab_error("STALE_FIGURE", stale, {
                    "kind": kind, "result_id": record.get("result_id"),
                    "required_tools": ["dataset_plot"]})
            if species_pca_id is not None and record.get("result_id") != species_pca_id:
                refusal = _other_species_pca(ds, record, species_pca_id, holder=holder)
                if refusal is not None:
                    return refusal
                notes_out.append(f"result_id={species_pca_id!r} は置き場の図の結果 {record.get('result_id')} と同じ"
                                 "入力・指定の分子種の PCA なので、置き場の図を保存した。")
        return _save_record(ds, record, analysis_id, title, holder, notes_out)
    except plot_tools._Missing as exc:
        return mcp_errors.missing_state(exc.state, exc.required_tools, exc.message)
    except DomainError as exc:
        return contrast_envelope(exc, "save_figure")


def _other_species_pca(ds, record: dict, requested: str, *, holder=None) -> str | None:
    """置き場の分子種の PCA の図が、`save_figure(result_id=...)` が名指しした結果と別の結果から描いたもののときの封筒。

    置き場の図は今のもの（呼ぶ前に `_stale_figure` が確かめる）。名指しした結果が置き場の図の結果と同じ計算
    （同じデータセット・入力の行列の指紋・項目の指定・試料の情報の版。同じ指定で描き直しただけ）なら None
    （置き場の図はその結果の図と同じなので保存してよい）。分子種の PCA の結果でなければ `DATASET_BAD_REQUEST`。
    それ以外は `STALE_FIGURE` で、違いの種類（`details.mismatch`）を先に決めてから案内を選ぶ:

    - `stale`: 名指しした結果が今のデータセットの状態から出たものでない（`is_current` が偽。そのあとの除外・
      試料の表の更新か、別のデータセットの結果）。その図はもう描けない（同じ指定で描き直しても今の状態の新しい
      結果になる）。置き場の図を保存するなら `details.held_result_id` を渡すか省く。
    - `input`: 今の結果で項目の指定も置き場の図と同じだが、入力の行列（選んだ分子種の値。キュレーションの判断の
      更新で変わり、`is_current` はそれを古さと読まない）が違う。これももう描けない。案内は `stale` と同じ。
    - `settings`: 今の結果で、項目の指定が置き場の図と違う。その指定で描き直してから保存する。

    名指しした結果の項目の指定が置き場の図と違うときだけ `details.settings`（名指しした結果の指定。そのまま
    `dataset_plot(kind="species_pca", ...)` に渡せる）と `details.required_tools=["dataset_plot"]` を付ける。同じ指定の
    描き直しは置き場の図と同じ図を作るだけで、名指しした結果には戻らないので案内しない。
    """
    from metabolomix.analysis.result_state import is_current

    held = record.get("result_id")
    result = (ds.results or {}).get(requested)
    held_result = (ds.results or {}).get(held) or record.get("result") or {}
    if result is not None and _same_species_pca(result, held_result):
        return None
    if result is None:
        # このハンドルの結果でない（別のハンドルの結果か、無い ID）。開いているほかのハンドルが分子種の PCA の
        # 結果として持っていればそのハンドルを名指す（読むだけ。食い違っても文面だけ）。ほかのハンドルの別の種類の
        # 結果（比較・PCA・概観）はそのハンドルで呼び直しても分子種の PCA の図にならないので、呼び直しを案内せず
        # 「分子種の PCA の結果ではない」と言い、持っているハンドルを文で名指す。
        open_handles = handles.registry.list()
        here = getattr(holder, "handle_id", None) or next(
            (h.handle_id for h in open_handles if h.dataset is ds), None)
        found = [(h.handle_id, ((getattr(h.dataset, "results", None) or {}).get(requested) or {}))
                 for h in open_handles
                 if h.handle_id != here and h.dataset is not None and h.dataset is not ds
                 and requested in (getattr(h.dataset, "results", None) or {})]
        owner = next((hid for hid, res in found if (res.get("provenance") or {}).get("kind") == "species_pca"),
                     None)
        other = None if owner else next(iter(found), None)
        details = {"kind": "species_pca", "result_id": requested, "held_result_id": held,
                   **({"dataset_handle": here} if here else {}),
                   **({"result_handle": owner} if owner else {})}
        where = f"ハンドル {here} " if here else "このデータセット"
        if owner:
            message = (f"result_id={requested} は{where}の結果ではなく、ハンドル {owner} の結果です。その図を保存するなら "
                       f"save_figure(..., result_id={requested!r}, dataset_handle={owner!r}) で呼んでください"
                       f"（{where}の置き場の分子種の PCA の図は結果 {held} から描いたもの）。")
        elif other is not None:
            kind = (other[1].get("provenance") or {}).get("kind") or other[1].get("kind")
            details["result_kind"] = kind
            message = (f"result_id={requested} は分子種の PCA の結果ではありません（ハンドル {other[0]} の "
                       f"{kind} の結果。{where}の置き場の分子種の PCA の図は結果 {held} から描いたもの）。"
                       "dataset_plot(kind=\"species_pca\") の戻り値の result_id を渡すか、省いて置き場の図を保存して"
                       "ください。")
        else:
            message = (f"result_id={requested} は{where}の結果にありません（{where}の置き場の分子種の PCA の図は結果 "
                       f"{held} から描いたもの）。別のハンドルで描いた図なら、そのハンドルを dataset_handle に渡してください。")
        return mcp_errors.mztab_error("DATASET_BAD_REQUEST", message, details)
    if (result.get("provenance") or {}).get("kind") != "species_pca":
        return mcp_errors.mztab_error(
            "DATASET_BAD_REQUEST",
            f"result_id={requested} は分子種の PCA の結果ではありません（置き場の分子種の PCA の図は結果 {held} から"
            "描いたもの）。dataset_plot(kind=\"species_pca\") の戻り値の result_id を渡すか、省いて置き場の図を"
            "保存してください。",
            {"kind": "species_pca", "result_id": requested, "held_result_id": held})
    settings = result.get("settings") or {}
    other_settings = settings != (held_result.get("settings") or {})
    details = {"kind": "species_pca", "result_id": requested, "held_result_id": held}
    redraw = (f"{requested} と同じ項目の指定を今のデータで描くなら、details.settings の指定で "
              "dataset_plot(kind=\"species_pca\", ...) を描き直し、その戻り値の result_id で保存してください"
              "（result_id を省いても描き直した図を保存します）。")
    keep_held = (f"置き場の図（結果 {held}。今のデータセットの状態から描いたもの）を保存するなら、result_id に "
                 "details.held_result_id を渡すか省いてください。")
    if not is_current(ds, result):
        details["mismatch"] = "stale"
        message = (f"名指しした結果 {requested} は今のデータセットの状態から出たものではありません（そのあとに試料の"
                   "情報が除外・試料の表の更新で変わったか、別のデータセットの結果）。その図はもう描けません"
                   f"（同じ指定で描き直しても今の状態の新しい結果になり、{requested} にはなりません）。" + keep_held)
        if other_settings:
            message += " " + redraw
    elif not other_settings:
        details["mismatch"] = "input"
        message = (f"名指しした結果 {requested} と置き場の図の結果 {held} は同じ項目の指定ですが、入力の行列（選んだ"
                   "分子種の値。キュレーションの判断の更新などで変わる）が違います。"
                   f"{requested} の図はもう描けません（同じ指定で描き直しても今の入力の結果になります）。" + keep_held)
    else:
        details["mismatch"] = "settings"
        message = (f"置き場の分子種の PCA の図は結果 {held} から描いたもので、名指しした結果 {requested} の図では"
                   "ありません（項目の指定が違う）。details.settings の指定で dataset_plot(kind=\"species_pca\", ...) を"
                   f"描き直してから保存してください（{requested} のあとに入力の行列〔キュレーションの判断など〕が"
                   f"変わっていれば、描き直した図は {requested} と別の結果になります。そのときは描き直した結果の "
                   "result_id を渡すか省いてください）。")
    if other_settings:
        details.update({"settings": settings, "required_tools": ["dataset_plot"]})
    return mcp_errors.mztab_error("STALE_FIGURE", message, details)


def _same_species_pca(a: dict, b: dict) -> bool:
    """2 つの分子種の PCA の結果が同じ計算か（同じデータセット・入力の行列・項目の指定・試料の情報の版）。"""
    pa, pb = a.get("provenance") or {}, b.get("provenance") or {}
    keys = ("kind", "dataset_id", "input_fingerprint", "metadata_revision")
    return (pa.get("kind") == "species_pca" and all(pa.get(k) == pb.get(k) for k in keys)
            and a.get("settings") == b.get("settings"))


#: 図を捨てた理由のうち、試料の表の列でないもの（`result_state.invalidate_results` の `changed`）。
_DISCARD_CAUSES = {"exclusion": "除外（dataset_exclude）", "dataset": "データの入れ替え",
                   "detection": "検出の印の変更", "comparison": "比較の指定の変更",
                   "pca_settings": "PCA の設定の変更"}


def _discarded_message(kind: str, changed) -> str:
    """捨てた図（`ds.discarded_figures[kind]`）の STALE_FIGURE の文言。"""
    causes = [_DISCARD_CAUSES[c] for c in changed if c in _DISCARD_CAUSES]
    columns = [str(c) for c in changed if c not in _DISCARD_CAUSES]
    if columns:
        causes.append(f"試料の表の更新（変わった列: {', '.join(columns)}）")
    cause = "・".join(causes) or "解析の入力の変更"
    return (f"{_LABELS[kind]}の図は、{cause}で元の結果ごと古くなったので捨てました。前処理・解析をやり直してから "
            f"dataset_plot(kind=\"{kind}\") で描き直してください。")


def _stale_figure(ds, record: dict) -> str | None:
    """結果から作った図が、今の結果から出たものでなければ理由を返す（`invalidate_results` が拾えない
    再前処理などの取りこぼしを、保存の直前に塞ぐ）。"""
    from metabolomix.analysis.result_state import is_current

    revision = record.get("metadata_revision")
    if revision is not None and revision != getattr(ds, "metadata_revision", 0):
        return (f"この図は試料の情報（群・バッチ・役割）を色に焼き込んだものですが、その後に試料の情報が更新されました"
                f"（実験情報シートの適用・除外など）。dataset_plot(kind=\"{record['kind']}\") で作り直してください。")
    result_id = record.get("result_id")
    if record["kind"] not in _RESULT_DERIVED_KINDS or result_id is None:
        return None
    result = (ds.results or {}).get(result_id)
    if result is None or not is_current(ds, result):
        return (f"この図は結果 {result_id} から作ったものですが、その結果はいまのデータセットの状態から出たものでは"
                f"ありません。dataset_plot(kind=\"{record['kind']}\") で作り直してください。")
    return None


def _save_record(ds, record: dict, analysis_id: str, title: str | None, holder,
                 notes_out: list[str] | None = None) -> str:
    """図の記録を `reports/figures/` に保存し、どの結果から描いたかを添えて返す（`notes_out` は使わなかった
    引数の注記。最後の `notes:` 行になる）。"""
    from metabolomix.tools.plot_tools import write_figure_files

    kind = record["kind"]
    out_path = _figures_dir(holder) / f"{notes.make_slug(analysis_id)}_{kind}.png"
    written = write_figure_files(ds, record, out_path, title=title)
    both = "（同名の .svg も保存）" if len(written) > 1 else ""
    label = _LABELS[kind]
    origin = f" result_id={record['result_id']}" if record.get("result_id") else ""
    text = (f"{label}図を保存: {out_path}{both}{origin}\n"
            f"本文に ![{label}](figures/{out_path.name}) で埋め込めます。")
    if notes_out:
        text += "\nnotes: " + " ".join(notes_out)
    return text


def _figures_dir(holder):
    figures_dir = _resolve_report_dir(holder) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    return figures_dir
