"""EIC/AEF ツール群: eic_parser, eic_plot_chromatograms, eic_plot_compounds,
eic_rank_by_max_intensity, m/z 範囲検索, RT 範囲検索。

deps: mcp_core / session_state / path_resolvers / metabolomix.eic.reader /
metabolomix.eic.identity_map / metabolomix.plots.eic。tools_* / server は import しない。
"""
from pathlib import Path

from metabolomix.core import session_state
from mcp.types import ToolAnnotations
from metabolomix.core.mcp_core import mcp
from metabolomix.core.serialization import json_payload
from metabolomix.core.path_resolvers import resolve_arf2_file_path, resolve_eicaef_file_path
from metabolomix.eic.reader import (
    read_eic_spot_css1,
    read_eic_spots_css1,
    summarize_eic_data,
    top_eic_spots_by_max_intensity,
    search_eic_by_mz_range,
    search_eic_by_rt_range,
)
from metabolomix.eic.identity_map import load_arf2_records, select_identity_candidates
from metabolomix.plots import render as plot_render
from metabolomix.plots.eic import (
    EICMultiPlotPayload,
    build_eic_plot_payload,
    build_multi_compound_plot_payload,
    render_eic_plot,
)
from mcp.server.fastmcp import Image

# 1枚に重ねて読める本数の上限。旧既定 24 は payload を実測 119,768 字まで膨らませ、
# 図としても凡例で埋まって判読できなかった。
DEFAULT_COMPOUND_TOP_N = 8

# EIC 検索の既定 m/z 上限。リピドミクスでは TG / DGDG / AHexCer / CL など中性・複合脂質
# が m/z 1000 を超える。旧既定 1000 は実データ（liver POS）で同定済み 1,994 スポットを
# 含む 13% を黙って範囲外にしていた。実測の最大 m/z は約 1,250。
DEFAULT_MAX_MZ = 1500.0


def _spot_line(spot: dict, with_intensity: bool = False) -> str:
    """検索結果 1 行。float をそのまま埋めると 17 桁出る（rt=12.376221656799316）ため、
    RT は 3 桁（0.001 分 = 60 ms）、m/z は慣例どおり 4 桁、強度は整数に丸める。"""
    line = (
        f"spot_id={spot['spot_id']} rt={spot['rt']:.3f} mz={spot['mz']:.4f}"
    )
    if with_intensity:
        line += f" max_intensity={spot['max_intensity']:.0f}"
    return line + f" num_samples={spot['num_samples']}"

__all__ = [
    "eic_parser",
    "eic_plot_chromatograms",
    "eic_plot_compounds",
    "eic_rank_by_max_intensity",
    "eic_search_by_mz_range",
    "eic_search_by_rt_range",
]


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_parser(file_path: str | None = None) -> str:
    """
    .EIC.aef ファイルに対応する解析用関数
    ファイルを解析し、テキスト要約を返します。
    """
    file_path = resolve_eicaef_file_path(file_path)
    if not file_path:
        return "データディレクトリに .aef ファイルが見つかりませんでした。"

    try:
        parsed = session_state.session.eic.load_data(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        summary = summarize_eic_data(parsed)
        output_text = [
            f"EIC解析完了: {Path(file_path).name}",
            "== 基本要約 ==",
            json_payload(summary),
            "== コメント ==",
            "この結果をもとに、m/z範囲検索やRT範囲検索、強度上位抽出（eic_rank_by_max_intensity）を実行できます。",
        ]
        return session_state.session.maybe_prepend_caveat(
            "\n".join(output_text), topic="eic",
        )
    except Exception as e:
        return f"EIC解析中にエラーが発生しました: {str(e)}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_plot_chromatograms(
    spot_id: int,
    file_path: str | None = None,
    file_ids: list[int] | None = None,
    normalize: str = "none",
    title: str | None = None,
) -> str:
    """Return client-neutral plot data for selected EIC chromatograms.

    This read-only tool returns structured ``lipidmix.eic.v1`` JSON only. It does
    not render an image and does not write files: Use-LLLM may render the series
    with Plotly, while Claude Desktop or another MCP client may use its own UI.
    Call ``save_figure(kind="eic")`` only after the user explicitly requests PNG output.

    ``spot_id`` is exact and should normally come from the existing m/z or RT
    search tools. ``file_ids`` selects at most 12 sample traces; when omitted all
    traces are returned only if the spot contains 12 samples or fewer.
    """
    resolved = resolve_eicaef_file_path(file_path)
    if not resolved:
        raise FileNotFoundError("データディレクトリに .aef ファイルが見つかりませんでした。")
    spot = read_eic_spot_css1(resolved, spot_id, file_ids=file_ids)
    payload = build_eic_plot_payload(
        spot, resolved, normalize=normalize, title=title,
    )
    session_state.session.eic.last_plot = payload
    return json_payload(payload)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_plot_compounds(
    file_id: int,
    names: list[str] | None = None,
    ontologies: list[str] | None = None,
    file_path: str | None = None,
    arf2_path: str | None = None,
    normalize: str = "none",
    top_n: int = DEFAULT_COMPOUND_TOP_N,
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """Overlay several identified compounds' EIC traces for ONE sample — as a PNG image by default.

    Read-only; writes no file.

    ``output="image"`` (the default) renders the overlay server-side and returns it
    as an image block plus a one-line caption, so the figure appears directly in the
    chat. ``output="payload"`` instead returns the client-neutral
    ``lipidmix.eic.multi.v1`` JSON for a client that draws its own interactive chart
    (Use-LLLM's Plotly view); the deployment default can be flipped with the
    ``LIPIDMIX_PLOT_OUTPUT`` environment variable. Prefer the image — one overlay's
    coordinate payload costs tens of thousands of tokens and an LLM cannot read a
    trace usefully.

    One call plots one sample: ``file_id`` is required. To compare samples, call
    this tool once per sample and show the figures side by side. ``top_n`` caps how
    many compounds are overlaid (default 8: more than that is unreadable in one
    figure and multiplies the payload).

    Compounds are selected from the ARF2 annotation: ``names`` matches ARF2
    ``Name`` as a case-insensitive substring and ``ontologies`` matches ARF2
    ``Ontology`` exactly (case-insensitive); the two results are unioned. At
    least one of them is required.

    ``spot_id = AlignmentID`` is verified against the EIC spot RT and m/z, so a
    compound whose identity cannot be confirmed is not drawn. Everything that was
    excluded — mismatch, missing trace, or falling outside ``top_n`` — is listed
    with its reason in ``selection.dropped``; read it before concluding that a
    compound is absent from the sample.
    """
    if isinstance(file_id, bool) or not isinstance(file_id, int):
        raise ValueError("file_id must be an integer")
    if normalize not in {"none", "per_trace_max"}:
        raise ValueError("normalize must be 'none' or 'per_trace_max'")
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    mode = plot_render.resolve_plot_output(output)

    resolved = resolve_eicaef_file_path(file_path)
    if not resolved:
        raise FileNotFoundError("データディレクトリに .aef ファイルが見つかりませんでした。")
    resolved_arf2 = resolve_arf2_file_path(arf2_path)
    if not resolved_arf2:
        raise FileNotFoundError("データディレクトリに .arf2 ファイルが見つかりませんでした。")

    records = load_arf2_records(resolved_arf2)
    candidates, caveats, total_matched = select_identity_candidates(
        records, names=names, ontologies=ontologies,
    )
    if not candidates:
        raise ValueError(
            f"names={names} / ontologies={ontologies} に一致する ARF2 スポットがありません。"
            "arf2_parser で脂質名・オントロジーの表記を確認してください。"
        )

    spots = read_eic_spots_css1(
        resolved,
        [candidate["spot_id"] for candidate in candidates],
        file_ids=[file_id],
    )
    payload = build_multi_compound_plot_payload(
        candidates,
        spots,
        file_id=file_id,
        file_path=resolved,
        arf2_path=resolved_arf2,
        normalize=normalize,
        top_n=top_n,
        title=title,
        queries=names,
        ontologies=ontologies,
        caveats=caveats,
        total_matched=total_matched,
    )
    session_state.session.eic.last_plot = payload
    if mode == plot_render.PAYLOAD:
        # dict をそのまま返すと FastMCP が indent=2 で整形する。点列 payload では
        # 空白だけで数千〜数万字になるため、最小形の JSON 文字列で返す。
        return json_payload(payload)
    png = plot_render.figure_to_png(render_eic_plot(payload, title=title))
    return [_multi_plot_caption(payload), Image(data=png, format="png")]


def _multi_plot_caption(payload: EICMultiPlotPayload) -> str:
    """画像に添える1行要約。図から読み取れない選択の内訳を言葉で残す。"""
    selection = payload["selection"]
    drawn = "、".join(item["label"] for item in payload["series"]) or "なし"
    sample = payload["sample"].get("sample_name") or f"FileID {payload['sample']['file_id']}"
    caption = (
        f"EIC 重ね描き（{sample}）: 候補 {selection['candidates']} 件のうち "
        f"{selection['plotted']} 件を描画 — {drawn}。"
    )
    if selection["dropped"]:
        caption += (
            f" 除外 {len(selection['dropped'])} 件（理由は output='payload' の "
            "selection.dropped、または caveats を参照）。"
        )
    return caption


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_rank_by_max_intensity(file_path: str | None = None, top_n: int = 20) -> str:
    """EICスポットを強度（各サンプルのクロマトグラム最大強度の最大値）で降順に返します。

    強度上位のスポットを抽出します。旧 eicaef_top_peak_tops は peak_top（ピーク頂点の
    横軸=RT座標。強度ではない）で並べていたため強度上位にならず、本ツールで置き換えました。
    """
    file_path = resolve_eicaef_file_path(file_path)
    if not file_path:
        return "データディレクトリに .aef ファイルが見つかりませんでした。"

    try:
        parsed = session_state.session.eic.load_data(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        top_spots = top_eic_spots_by_max_intensity(parsed, top_n=top_n)
        output_text = [
            f"EIC強度上位スポット: {Path(file_path).name}",
            f"上位{top_n}件（max_intensity 降順）:",
        ]
        for spot in top_spots:
            output_text.append(_spot_line(spot, with_intensity=True))
        return "\n".join(output_text)
    except Exception as e:
        return f"EIC強度上位抽出中にエラーが発生しました: {str(e)}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_search_by_mz_range(file_path: str | None = None, min_mz: float = 0.0,
                           max_mz: float = DEFAULT_MAX_MZ, max_results: int = 20) -> str:
    """EICデータのm/z範囲でスポットを検索します。

    既定の上限は 1500。リピドミクスでは TG / DGDG / CL など m/z 1000 超の脂質が
    普通にあるため（実データで全スポットの 5〜13%、その大半が同定済み）、既定で
    そこを切り落とさない。
    """
    file_path = resolve_eicaef_file_path(file_path)
    if not file_path:
        return "データディレクトリに .aef ファイルが見つかりませんでした。"

    try:
        parsed = session_state.session.eic.load_data(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        matches = search_eic_by_mz_range(parsed, min_mz, max_mz)
        output_text = [
            f"EIC m/z範囲検索: {min_mz} - {max_mz}",
            f"一致件数: {len(matches)}",
            "上位結果:",
        ]
        for spot in matches[:max_results]:
            output_text.append(_spot_line(spot))
        if len(matches) > max_results:
            output_text.append(f"(表示上限: {max_results} 件)")
        return "\n".join(output_text)
    except Exception as e:
        return f"EIC m/z検索中にエラーが発生しました: {str(e)}"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_search_by_rt_range(file_path: str | None = None, min_rt: float = 0.0, max_rt: float = 20.0, max_results: int = 20) -> str:
    """
    EICデータのRT範囲でスポットを検索します。
    """
    file_path = resolve_eicaef_file_path(file_path)
    if not file_path:
        return "データディレクトリに .aef ファイルが見つかりませんでした。"

    try:
        parsed = session_state.session.eic.load_data(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        matches = search_eic_by_rt_range(parsed, min_rt, max_rt)
        output_text = [
            f"EIC RT範囲検索: {min_rt} - {max_rt}",
            f"一致件数: {len(matches)}",
            "上位結果:",
        ]
        for spot in matches[:max_results]:
            output_text.append(_spot_line(spot))
        if len(matches) > max_results:
            output_text.append(f"(表示上限: {max_results} 件)")
        return "\n".join(output_text)
    except Exception as e:
        return f"EIC RT検索中にエラーが発生しました: {str(e)}"
