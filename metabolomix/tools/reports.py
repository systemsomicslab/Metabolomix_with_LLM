"""解析・解釈レポートと図の保存ツール群。

write/read/list_reports, save_figure。
レポート先の解決は mcp_core（DATA_DIR を動的参照）に委ねる。deps: mcp_core /
session_state / tool_helpers / knowledge_store / matplotlib。
tools_* / server は import しない。
"""
import matplotlib.pyplot as plt

from metabolomix.corpus import knowledge_store
from metabolomix.core import mcp_errors
from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.serialization import json_payload
from metabolomix.core.mcp_core import mcp, _resolve_report_dir, _report_dir_candidates, _build_report_meta
from metabolomix.core.tool_helpers import _pca_scatter_arrays, dataset_pca_plot
from metabolomix.plots.eic import render_eic_plot
from metabolomix.plots.volcano import render_volcano_plot

__all__ = [
    "write_report",
    "read_report",
    "list_reports",
    "save_figure",
]


# --- 図の入力元の選択（ARF 経路 / mzTab-M 経路） ---
#
# 図の描画は 2 つのセッションスロットから来る。ARF 経路は session.arf に、
# mzTab-M 経路は session.dataset（DatasetState）に結果を置く。
#
# **どちらかを優先しない**。旧実装は両方あるとき黙って ARF を採っていたので、
# 新しく読み込んだ mzTab の解析をしているつもりでも前のデータの図が保存され得た。
# 候補を並べて select_result に決めさせ、決まらなければ描かずに止める。

def _pca_candidates() -> list[dict]:
    from metabolomix.analysis.result_state import is_current

    candidates: list[dict] = []
    plot = getattr(session_state.session.arf, "last_pca_plot", None)
    if plot and plot.get("points"):
        prov = plot.get("provenance") or {}
        candidates.append({
            "source": "arf", "result_id": prov.get("result_id"),
            "dataset_id": prov.get("dataset_id"), "valid": True,
            "result": plot, "ds": None})
    species = getattr(session_state.session.arf, "last_species_pca", None)
    if species and species.get("points"):
        prov = species.get("provenance") or {}
        candidates.append({
            "source": "species", "result_id": prov.get("result_id"),
            "dataset_id": prov.get("dataset_id"), "valid": True,
            "result": species, "ds": None})
    ds = getattr(session_state.session, "dataset", None)
    ds_pca = getattr(ds, "last_pca", None) if ds is not None else None
    if ds_pca:
        projected = dataset_pca_plot(ds_pca)
        if projected:
            prov = ds_pca.get("provenance") or {}
            candidates.append({
                "source": "mztab", "result_id": prov.get("result_id"),
                "dataset_id": prov.get("dataset_id"),
                "valid": is_current(ds, ds_pca), "result": projected, "ds": ds})
    return candidates


def _differential_candidates() -> list[dict]:
    from metabolomix.analysis.result_state import is_current

    candidates: list[dict] = []
    last = getattr(session_state.session.arf, "last_differential", None)
    if last and last.get("volcano"):
        prov = last.get("provenance") or {}
        candidates.append({
            "source": "arf", "result_id": prov.get("result_id"),
            "dataset_id": prov.get("dataset_id"), "valid": True,
            "result": last, "ds": None})
    ds = getattr(session_state.session, "dataset", None)
    ds_diff = getattr(ds, "last_differential", None) if ds is not None else None
    if ds_diff and ds_diff.get("volcano"):
        prov = ds_diff.get("provenance") or {}
        candidates.append({
            "source": "mztab", "result_id": prov.get("result_id"),
            "dataset_id": prov.get("dataset_id"),
            "valid": is_current(ds, ds_diff), "result": ds_diff, "ds": ds})
    return candidates


def _figure_missing_state(kind: str) -> str:
    if kind == "pca":
        return mcp_errors.missing_state(
            "pca_result",
            ["arf_parser", "arf_pca_preprocessed", "arf_pca_species", "load_dataset", "dataset_pca"],
            "先に arf_parser / arf_pca_preprocessed / arf_pca_species / load_dataset 等でPCAを実行してください"
            "（PCA結果がありません）。")
    return mcp_errors.missing_state(
        "differential_result", ["arf_differential", "dataset_differential"],
        "[error] 直近の差次的解析（volcano データ）がありません。"
        "先に arf_differential または dataset_differential を実行してください。")


def _choose_figure_result(kind: str, source: str, result_id: str | None):
    """候補から 1 件選ぶ。選べないときはエラーエンベロープ文字列を返す。"""
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.plots.result_output import select_result

    candidates = (_pca_candidates() if kind == "pca" else _differential_candidates())
    if not candidates:
        return _figure_missing_state(kind)
    try:
        return select_result(candidates, source=source, result_id=result_id)
    except DomainError as exc:
        return mcp_errors.mztab_error(exc.code, exc.message, exc.details or None)


# --- 解析・解釈レポート（reports/<analysis_id>.md） ---
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
def write_report(
    analysis_id: str,
    dataset: str,
    body: str,
    status: str = "draft",
    knowledge_refs: list[str] | None = None,
) -> str:
    """解析・解釈レポートを reports/<analysis_id>.md に上書き保存する（成果物＋記録）。

    body は frontmatter を含まない Markdown 本文。推奨セクション見出し:
    `## 目的` / `## 実施した解析` / `## 主要な所見` / `## 解釈` /
    `## 注意点・コンフリクト` / `## 結論`。所見が増えたら本文を作り直して再度呼ぶ
    （ファイルは毎回上書き）。引用した knowledge/playbook の slug を knowledge_refs に渡す。
    LIPIDMIX_REPORTS_DIRの明示先を優先し、書込不可なら解析フォルダ配下reports/。
    未指定時は解析フォルダ配下reports/を優先し、不可なら<project>/reports。

    analysis_id はファイル名 slug の元になるため ASCII で一意に。別IDが同一slugに潰れて既存
    レポートを上書きしそうな場合は保存を中止して通知する（objectiveレコードと同じIDを推奨）。
    """
    slug = knowledge_store.make_slug(analysis_id)
    # slug衝突ガード: 別の analysis_id が同一ファイル名に潰れる場合は、既存レポートを
    # 黙って上書きせず中止する（同一IDの上書きは意図どおり許可）。非ASCII/記号違いのIDで起きうる。
    for directory in _report_dir_candidates():
        existing = directory / f"{slug}.md"
        if existing.is_file():
            ex_meta, _ = knowledge_store.parse_frontmatter(existing.read_text(encoding="utf-8"))
            ex_id = str(ex_meta.get("analysis_id", ""))
            if ex_id and ex_id != analysis_id:
                return (
                    f"slug衝突のため中止: analysis_id '{analysis_id}' はファイル名 '{slug}.md' に潰れますが、"
                    f"そこには別の '{ex_id}' のレポートが既にあります（{existing}）。"
                    "上書きを避けました。ASCIIで一意な analysis_id を指定してください。"
                )
            break  # 同一ID → 上書きしてよい
    reports_dir = _resolve_report_dir()
    meta = _build_report_meta(analysis_id, dataset, status, knowledge_refs)
    path = knowledge_store.write_note(reports_dir, slug, meta, body)
    return f"レポートを保存: {path}（status={status}）。read_report('{analysis_id}') で読み戻せます。"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def read_report(analysis_id: str) -> str:
    """過去レポートを読み戻す（候補ディレクトリ横断で最新更新のものを返す）。セッション継続用。

    解析フォルダ配下と退避先の両方に同名レポートが残る場合（書き込み可否が途中で変化した等）、
    古い方を返さないよう mtime が最新のファイルを採用する。
    """
    slug = knowledge_store.make_slug(analysis_id)
    matches = [d / f"{slug}.md" for d in _report_dir_candidates()]
    matches = [p for p in matches if p.is_file()]
    if not matches:
        return f"レポートが見つかりません: {analysis_id}（write_report で作成してください）"
    newest = max(matches, key=lambda p: p.stat().st_mtime)
    return newest.read_text(encoding="utf-8")


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def list_reports() -> str:
    """既存レポートの1行索引（analysis_id / date / status）を返す。

    候補ディレクトリ横断で同一 analysis_id が重複する場合は mtime が最新の1件を採用する。
    """
    # analysis_id -> (mtime, 表示行)。最新更新の行を残す。
    best: dict[str, tuple[float, str]] = {}
    for directory in _report_dir_candidates():
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.md")):
            meta, _body = knowledge_store.parse_frontmatter(path.read_text(encoding="utf-8"))
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


#: save_figure が受け付ける図の種類。
FIGURE_KINDS = ("pca", "volcano", "eic", "group_intensity", "species", "pca_loadings")
#: ARF と mzTab のどちらの結果を描くかを選ぶ必要がある種類（source / result_id を使う）。
_SOURCE_KINDS = ("pca", "volcano")


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
def save_figure(kind: str, analysis_id: str, title: str | None = None,
                source: str = "auto", result_id: str | None = None) -> str:
    """明示的なユーザー要求時だけ、直近の図を reports/figures/<analysis_id>_<kind>.png に保存する。

    kind: "pca"（arf_parser / arf_pca_preprocessed / arf_pca_species / dataset_pca のスコア。
      arf_pca_species の結果は source="species"）/ "volcano"
      （arf_differential / dataset_differential。間引き前の全点）/ "eic"（eic_plot_chromatograms /
      eic_plot_compounds）/ "group_intensity"（arf_plot_group_intensity。dpi 300 の PNG と同名 .svg）/
      "species"（arf_plot_species。dpi 300 の PNG と同名 .svg）/
      "pca_loadings"（plot_pca_loadings。dpi 300 の PNG と同名 .svg）。
    source / result_id: kind が pca / volcano のときだけ使う。"auto"（既定）はどちらも優先せず、
      有効な結果が 2 つ以上あると AMBIGUOUS_RESULT_SOURCE で止まる（"arf" / "species"（arf_pca_species）/
      "mztab" で指定する）。
    通常の対話描画では呼ばない（各描画ツールが画像を返す）。返り値の相対パスは write_report の
    本文に `![...](figures/<analysis_id>_<kind>.png)` として埋め込める。
    """
    if kind not in FIGURE_KINDS:
        return json_payload({"status": "error", "message":
                             f"kind は {', '.join(FIGURE_KINDS)} のどれかを指定してください（受け取った値: {kind!r}）。"})
    if kind not in _SOURCE_KINDS and (source != "auto" or result_id is not None):
        return json_payload({"status": "error", "message":
                             f"source / result_id は kind が {' / '.join(_SOURCE_KINDS)} のときだけ指定できます"
                             f"（kind={kind!r}）。"})
    if kind == "pca":
        chosen = _choose_figure_result("pca", source, result_id)
        if isinstance(chosen, str):
            return chosen
        return _save_figure(chosen, analysis_id, title, kind="pca", label="PCA")
    if kind == "volcano":
        chosen = _choose_figure_result("differential", source, result_id)
        if isinstance(chosen, str):
            return chosen
        return _save_figure(chosen, analysis_id, title, kind="volcano", label="volcano")
    if kind == "eic":
        return _save_eic(analysis_id, title)
    if kind == "group_intensity":
        return _save_group_intensity(analysis_id, title)
    if kind == "species":
        return _save_species(analysis_id, title)
    if kind == "pca_loadings":
        return _save_pca_loadings(analysis_id, title)
    raise AssertionError(f"unhandled figure kind: {kind!r}")  # FIGURE_KINDS に足したら分岐も足す


def _save_figure(chosen: dict, analysis_id: str, title: str | None, *,
                 kind: str, label: str) -> str:
    """選ばれた結果を PNG にして、どの結果から描いたかを添えて返す。"""
    from metabolomix.plots.result_output import save_result_figure

    slug = knowledge_store.make_slug(analysis_id)
    figures_dir = _resolve_report_dir() / "figures"
    out_path = save_result_figure(
        chosen.get("ds"), chosen["result"],
        figures_dir / f"{slug}_{kind if kind != 'volcano' else 'volcano'}.png",
        kind=kind, title=title)

    rel = f"figures/{out_path.name}"
    result_id = chosen.get("result_id") or "(なし)"
    return (f"{label}図を保存: {out_path}"
            f"（source={chosen['source']} result_id={result_id}）\n"
            f"本文に ![{label}]({rel}) で埋め込めます。")


def _save_eic(analysis_id: str, title: str | None) -> str:
    """直近EICプロット情報（eic_plot_chromatograms / eic_plot_compounds）をPNGとして保存する。"""
    plot = getattr(session_state.session.eic, "last_plot", None)
    if not plot or not plot.get("series"):
        return mcp_errors.missing_state(
            "eic_plot", ["eic_plot_chromatograms", "eic_plot_compounds"],
            "先に eic_plot_chromatograms または eic_plot_compounds を実行してください"
            "（EICプロット情報がありません）。")

    slug = knowledge_store.make_slug(analysis_id)
    reports_dir = _resolve_report_dir()
    figures_dir = reports_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / f"{slug}_eic.png"
    fig = render_eic_plot(plot, title=title)
    try:
        fig.savefig(out_path, dpi=120, format="png", bbox_inches="tight")
    finally:
        plt.close(fig)

    rel = f"figures/{out_path.name}"
    return f"EIC図を保存: {out_path}\n本文に ![EIC]({rel}) で埋め込めます。"


def _save_group_intensity(analysis_id: str, title: str | None) -> str:
    """直前の arf_plot_group_intensity の図を PNG（dpi 300）と同名 .svg に保存する。"""
    from metabolomix.plots.group_intensity import render_group_intensity_plot

    last = getattr(session_state.session.arf, "last_group_intensity", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "group_intensity_plot", ["arf_plot_group_intensity"],
            "先に arf_plot_group_intensity を実行してください（群別強度の図がありません）。")
    out_path = _figures_dir() / f"{knowledge_store.make_slug(analysis_id)}_group_intensity.png"
    fig = render_group_intensity_plot(last["payload"], title=title or last.get("title"), ncols=last.get("ncols"))
    _write_png_and_svg(fig, out_path)
    rel = f"figures/{out_path.name}"
    return (f"群別強度の図を保存: {out_path}（同名の .svg も保存）\n"
            f"本文に ![group intensity]({rel}) で埋め込めます。")


def _save_species(analysis_id: str, title: str | None) -> str:
    """直前の arf_plot_species の図を PNG（dpi 300）と同名 .svg に保存する。"""
    from metabolomix.plots.species import render_species_plot

    last = getattr(session_state.session.arf, "last_species_plot", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "species_plot", ["arf_plot_species"],
            "先に arf_plot_species を実行してください（分子種ごとの図がありません）。")
    out_path = _figures_dir() / f"{knowledge_store.make_slug(analysis_id)}_species.png"
    fig = render_species_plot(last["payload"], title=title or last.get("title"), ncols=last.get("ncols"))
    _write_png_and_svg(fig, out_path)
    return (f"分子種ごとの図を保存: {out_path}（同名の .svg も保存）\n"
            f"本文に ![species](figures/{out_path.name}) で埋め込めます。")


def _save_pca_loadings(analysis_id: str, title: str | None) -> str:
    """直前の plot_pca_loadings の図を PNG（dpi 300）と同名 .svg に保存する。"""
    from metabolomix.plots.pca_loadings import render_loadings_plot

    last = getattr(session_state.session, "last_loadings_plot", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "pca_loadings_plot", ["plot_pca_loadings"],
            "先に plot_pca_loadings を実行してください（ローディング図がありません）。")
    out_path = _figures_dir() / f"{knowledge_store.make_slug(analysis_id)}_pca_loadings.png"
    _write_png_and_svg(render_loadings_plot(last["payload"], title=title or last.get("title")), out_path)
    payload = last["payload"]
    return (f"ローディング図を保存: {out_path}（同名の .svg も保存。source={payload['source']} "
            f"result_id={payload['result_id'] or '(なし)'}）\n"
            f"本文に ![PCA loadings](figures/{out_path.name}) で埋め込めます。")


def _figures_dir():
    figures_dir = _resolve_report_dir() / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    return figures_dir


def _write_png_and_svg(fig, out_path) -> None:
    """報告書用の図を dpi 300 の PNG と同名 SVG で書く（Figure は閉じる）。"""
    try:
        fig.savefig(out_path, dpi=300, format="png", bbox_inches="tight")
        fig.savefig(out_path.with_suffix(".svg"), format="svg", bbox_inches="tight")
    finally:
        plt.close(fig)
