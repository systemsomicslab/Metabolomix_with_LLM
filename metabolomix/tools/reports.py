"""解析・解釈レポートと図の保存ツール群。

write_report / read_report（ID を省くと一覧）/ save_figure。
レポート先は `[paths] reports_dir`（あれば先）→ データセットの解析フォルダ配下 `reports/`（mcp_core）。
各ツールは `dataset_handle` で読み書きするハンドルを選ぶ（省略時は直近）。save_figure は `dataset_plot` が置いた
`ds.figures` から描き、ファイルへの書き方は `plot_tools.write_figure_files`（`dataset_plot(output_path=...)` と共有）。
deps: mcp_core / plot_tools（遅延 import）/ notes。
server は import しない。
"""

from metabolomix.tools import notes
from metabolomix.core import mcp_errors
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
#: 結果なので `result_id` では描き直せない（`_RESULT_KINDS` には入れない）。
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
    （`dataset_set_sample_metadata`。前処理の入力の列なら pca / loadings / species_pca / volcano、群だけなら volcano）は
    図ごと捨てるが、捨てた図も `STALE_FIGURE`（`details.reason` に変わったもの: `exclusion` や列名）で、
    描いたことの無い図（missing_state）と見分けられる。group_intensity / species / species_pca は群を図に焼き込むので、
    試料の情報（群・役割）が更新されたあとも `STALE_FIGURE`（species_pca は結果から作った図でもある）。
    eic / eic_compounds は結果から作った図ではないので、描いたとおりに保存する（あとで `dataset_exclude` しても変わらない）。
    result_id: kind が pca / volcano / loadings で効く。登録された結果を名指しして描く
      （`dataset_plot` を呼ばなくてよい。古い結果は STALE_ANALYSIS_RESULT）。volcano では
      `"result_id#contrast_id"`（`dataset_feature_search` が出す形）でも対比を選べる。ほかの kind では使わず、
      置き場の図を保存して戻り値の `notes:` 行にそう書く。
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
    if result_id is not None and kind not in _RESULT_KINDS:
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
            stale = _stale_figure(ds, record)
            if stale:
                return mcp_errors.mztab_error("STALE_FIGURE", stale, {
                    "kind": kind, "result_id": record.get("result_id"),
                    "required_tools": ["dataset_plot"]})
        return _save_record(ds, record, analysis_id, title, holder, notes_out)
    except plot_tools._Missing as exc:
        return mcp_errors.missing_state(exc.state, exc.required_tools, exc.message)
    except DomainError as exc:
        return contrast_envelope(exc, "save_figure")


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
