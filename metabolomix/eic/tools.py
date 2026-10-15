"""EIC/AEF ツール群: eic_parser, eic_search（spot ID・m/z・RT の AND 検索と強度順）,
eic_plot_chromatograms, eic_plot_compounds。

deps: mcp_core / session_state / path_resolvers / metabolomix.eic.reader /
metabolomix.eic.identity_map / metabolomix.eic.mz_backfill / metabolomix.msdial.alignment_spot /
metabolomix.arf2.reader（TSV 整形）/ metabolomix.plots.eic。tools_* / server は import しない。
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
    search_eic_spots,
    summarize_eic_data,
)
from metabolomix.eic.identity_map import load_arf2_records, select_identity_candidates
from metabolomix.eic.mz_backfill import (
    backfill_caveat,
    backfill_missing_mz,
    has_missing_mz,
    sibling_arf2_path,
)
from metabolomix.msdial.alignment_spot import (
    alignment_sibling,
    msms_peak_counts,
    precursor_agrees,
)
from metabolomix.arf2.reader import format_spots_as_table
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

def _backfill_mz(
    spots: list[dict], eic_path: str, records: list[dict] | None = None,
    arf2_name: str | None = None,
) -> list[str]:
    """欠損 m/z（古い MS-DIAL の -1）を ARF2 で補い、注記を返す（mz_backfill 参照）。

    `records` を渡さなければ同じバッチの隣の `.arf2` を読む。
    """
    if not has_missing_mz(spots):
        return []
    if records is None:
        arf2 = sibling_arf2_path(eic_path)
        records = load_arf2_records(arf2) if arf2 else []
        arf2_name = arf2.name if arf2 else None
    note = backfill_caveat(backfill_missing_mz(spots, records), arf2_name)
    return [note] if note else []


def _load_session_spots(file_path: str) -> tuple[list, list[str]]:
    """セッションの EIC を読み、欠損 m/z を 1 度だけ補う。注記はキャッシュと一緒に残す。"""
    state = session_state.session.eic
    parsed = state.load_data(file_path)
    if isinstance(parsed, list) and state.mz_caveats is None:
        state.mz_caveats = _backfill_mz(parsed, file_path)
    return parsed, list(state.mz_caveats or [])


def _with_notes(lines: list[str], notes: list[str]) -> str:
    if notes:
        lines = lines + ["== 注意 =="] + notes
    return "\n".join(lines)


__all__ = [
    "eic_parser",
    "eic_plot_chromatograms",
    "eic_plot_compounds",
    "eic_search",
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
        parsed, notes = _load_session_spots(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        summary = summarize_eic_data(parsed)
        output_text = [
            f"EIC解析完了: {Path(file_path).name}",
            "== 基本要約 ==",
            json_payload(summary),
            "== コメント ==",
            "この結果をもとに、eic_search で spot ID・m/z 範囲・RT 範囲の検索や強度上位抽出（sort=\"max_intensity\"）を実行できます。",
        ]
        return session_state.session.maybe_prepend_caveat(
            _with_notes(output_text, notes), topic="eic",
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
    notes = _backfill_mz([spot], resolved)
    payload = build_eic_plot_payload(
        spot, resolved, normalize=normalize, title=title,
    )
    payload["caveats"].extend(notes)
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
    compound whose identity cannot be confirmed is not drawn. EIC files written
    by MS-DIAL before 2023-11-20 carry no spot m/z (-1); those spots borrow the
    ARF2 ``MassCenter`` after an RT check, and a caveat says so. Everything that was
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
    # 欠損 m/z を照合前に補う。補うときに同じ AlignmentID の RT 一致を確かめるので、
    # 以後の m/z 照合が自明に通っても RT による取り違え検出は残る。
    caveats = caveats + _backfill_mz(
        spots, resolved, records=records, arf2_name=Path(resolved_arf2).name,
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


_SORT_KEYS = ("spot_id", "max_intensity")
_SEARCH_COLUMNS = [
    "spot_id", "rt", "mz", "max_intensity", "num_samples",
    "name", "ontology", "adduct", "n_msms",
]


def _bound(low, high, unit: str) -> str:
    low_text = "-∞" if low is None else f"{low}"
    high_text = "∞" if high is None else f"{high}"
    return f"{unit} {low_text}–{high_text}"


def _search_condition(spot_ids, min_mz, max_mz, min_rt, max_rt, sort) -> str:
    parts = []
    if spot_ids is not None:
        parts.append(f"spot_ids={sorted(set(spot_ids))}")
    if min_mz is not None or max_mz is not None:
        parts.append(_bound(min_mz, max_mz, "m/z"))
    if min_rt is not None or max_rt is not None:
        parts.append(_bound(min_rt, max_rt, "RT") + " min")
    return (" かつ ".join(parts) or "条件なし（全件）") + f" / 並び: {sort}"


def _annotate_rows(spots: list[dict], eic_path: str) -> tuple[list[dict], list[str]]:
    """同じバッチの ARF2（同定名）とアラインメント .dcl（MS/MS ピーク数）を行へ付ける。"""
    notes: list[str] = []
    arf2 = sibling_arf2_path(eic_path)
    records = {}
    if arf2 is None:
        notes.append("同じバッチの .arf2 が無いため name / ontology / adduct は空。")
    else:
        for record in load_arf2_records(arf2):
            spot_id = record.get("AlignmentID")
            if isinstance(spot_id, int) and not isinstance(spot_id, bool):
                records.setdefault(spot_id, record)

    dcl = alignment_sibling(eic_path, ".dcl")
    counts = msms_peak_counts(dcl) if dcl is not None else None
    if counts is None:
        notes.append("同じバッチのアライメント .dcl が無いため n_msms は不明（-）。")

    rows, has_dash = [], False
    for spot in spots:
        record = records.get(spot["spot_id"], {})
        mz = spot["mz"] if spot["mz"] >= 0 else record.get("MassCenter")
        n_msms = "-"
        if counts is not None and spot["spot_id"] < len(counts):
            entry = counts[spot["spot_id"]]
            if entry["n_peaks"] > 0 and precursor_agrees(
                entry["precursor_mz"], mz, entry["rt"], spot["rt"],
            ):
                n_msms = entry["n_peaks"]
        has_dash = has_dash or n_msms == "-"
        rows.append({
            "spot_id": spot["spot_id"],
            "rt": round(float(spot["rt"]), 3),
            "mz": round(float(spot["mz"]), 4),
            "max_intensity": int(round(float(spot["max_intensity"]))),
            "num_samples": spot["num_samples"],
            "name": str(record.get("Name") or ""),
            "ontology": str(record.get("Ontology") or ""),
            "adduct": str(record.get("AdductType") or ""),
            "n_msms": n_msms,
        })
    if has_dash and counts is not None:
        notes.append(
            "n_msms の「-」は MS/MS 未取得か precursor 不一致（索引対応の崩れ）で、"
            "0 本という意味ではない。"
        )
    return rows, notes


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def eic_search(
    spot_ids: list[int] | None = None,
    min_mz: float | None = None,
    max_mz: float | None = None,
    min_rt: float | None = None,
    max_rt: float | None = None,
    sort: str = "spot_id",
    max_results: int = 20,
    file_path: str | None = None,
) -> str:
    """EIC スポットを spot ID・m/z 範囲・RT 範囲で検索し、TSV で返す。

    指定した条件はすべて AND で掛かる（何も指定しなければ全件）。境界値は含む。
    ``spot_ids`` は MasterAlignmentID（ARF / ARF2 / ``.dcl`` と共通）の完全一致で、
    範囲外の ID は注記に列挙する。``sort="max_intensity"`` は各試料のクロマトグラム
    最大強度の、そのまた最大値で降順に並べる（強度上位ランキング）。既定は spot_id 昇順。

    各行は ``spot_id rt mz max_intensity num_samples name ontology adduct n_msms``。
    ``name`` / ``ontology`` / ``adduct`` は同じバッチの ``.arf2``、``n_msms`` は
    アラインメント ``.dcl`` の代表試料 MS/MS のピーク数。``-`` は未取得または
    precursor 不一致で、0 本ではない。未注釈で ``n_msms`` がある行は
    ``library_match_feature(spot_id=...)`` → ``library_plot_mirror`` で照合・確認できる。
    """
    if sort not in _SORT_KEYS:
        return f"sort は {' / '.join(_SORT_KEYS)} のいずれかを指定してください（指定: {sort!r}）。"
    if min_mz is not None and max_mz is not None and min_mz > max_mz:
        return f"min_mz（{min_mz}）が max_mz（{max_mz}）より大きいため一致しません。"
    if min_rt is not None and max_rt is not None and min_rt > max_rt:
        return f"min_rt（{min_rt}）が max_rt（{max_rt}）より大きいため一致しません。"
    if isinstance(max_results, bool) or not isinstance(max_results, int) or max_results < 1:
        return "max_results は 1 以上の整数を指定してください。"

    file_path = resolve_eicaef_file_path(file_path)
    if not file_path:
        return "データディレクトリに .aef ファイルが見つかりませんでした。"

    try:
        parsed, notes = _load_session_spots(file_path)
        if not isinstance(parsed, list):
            return "EIC解析結果がリストではありません。"

        matches, unknown = search_eic_spots(
            parsed, spot_ids=spot_ids, min_mz=min_mz, max_mz=max_mz,
            min_rt=min_rt, max_rt=max_rt,
        )
        if unknown:
            notes.append(
                f"spot_ids のうち {unknown} はこの EIC に無い（範囲外。"
                f"スポットは 0〜{len(parsed) - 1}）。"
            )
        if sort == "max_intensity":
            matches.sort(key=lambda item: item["max_intensity"], reverse=True)
        shown = matches[:max_results]

        output_text = [
            f"EIC 検索: {Path(file_path).name} — "
            f"{_search_condition(spot_ids, min_mz, max_mz, min_rt, max_rt, sort)}",
            f"一致件数: {len(matches)}",
        ]
        if shown:
            rows, row_notes = _annotate_rows(shown, file_path)
            output_text.append(format_spots_as_table(rows, columns=_SEARCH_COLUMNS))
            notes.extend(row_notes)
        if len(matches) > max_results:
            output_text.append(f"(表示上限: {max_results} 件。条件か sort で絞り込める)")
        return _with_notes(output_text, notes)
    except Exception as e:
        return f"EIC検索中にエラーが発生しました: {str(e)}"
