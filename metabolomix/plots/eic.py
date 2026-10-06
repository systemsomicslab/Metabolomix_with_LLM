"""Client-neutral EIC plot payload and optional matplotlib rendering."""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

from metabolomix.core.serialization import round_floats
from metabolomix.eic.identity_map import IdentityCandidate, verify_spot_match
from metabolomix.msdial.classes import parse_analysis_file_classes, resolve_mddata_path


# 座標の丸め桁数。RT は分単位で 4 桁（0.0001 分 = 6 ms）あれば図も解釈も足り、
# m/z も 4 桁が慣例。float の既定 repr は 17 桁まで出すため、1 トレース 200 点
# × 60 サンプルでは桁の大半が無意味な尾数になる（実測 -36%）。
POINT_DIGITS = 4


class PlotAxis(TypedDict):
    label: str
    unit: str | None
    scale: str


class PlotAxes(TypedDict):
    x: PlotAxis
    y: PlotAxis


class PlotSource(TypedDict):
    file: str
    file_name: str


class MultiPlotSource(TypedDict):
    file: str
    file_name: str
    arf2_file: str


class PlotSpot(TypedDict):
    spot_id: int
    rt: float
    ri: float
    mz: float
    drift: float
    main_type: int
    total_samples: int
    selected_samples: int


class PlotSeries(TypedDict):
    id: str
    label: str
    file_id: int
    sample_name: str | None
    class_id: str | None
    x: list[float]
    y: list[float]
    peak_left: float
    peak_top: float
    peak_right: float
    max_intensity: float
    mean_intensity: float
    point_count: int


class RenderHints(TypedDict):
    mode: str
    connect_points: bool
    show_legend: bool
    hover_fields: list[str]


class EICPlotPayload(TypedDict):
    plot_schema: str
    plot_type: str
    title: str
    source: PlotSource
    axes: PlotAxes
    spot: PlotSpot
    normalization: str
    series: list[PlotSeries]
    render_hints: RenderHints
    caveats: list[str]


MULTI_PLOT_SCHEMA = "lipidmix.eic.multi.v1"


class PlotAnnotation(TypedDict):
    text: str
    x: float
    y: float


class PlotSample(TypedDict):
    file_id: int
    sample_name: str | None
    class_id: str | None


class MultiPlotSeries(TypedDict):
    id: str
    label: str
    spot_id: int
    name: str
    ontology: str
    adduct: str
    mz: float
    rt: float
    x: list[float]
    y: list[float]
    peak_left: float
    peak_top: float
    peak_right: float
    max_intensity: float
    mean_intensity: float
    point_count: int
    annotation: PlotAnnotation


class DroppedCompound(TypedDict):
    spot_id: int
    name: str
    reason: str


class PlotSelection(TypedDict):
    queries: list[str]
    ontologies: list[str]
    candidates: int
    candidates_evaluated: int
    plotted: int
    top_n: int
    dropped: list[DroppedCompound]


class MultiRenderHints(TypedDict):
    mode: str
    connect_points: bool
    show_legend: bool
    show_annotations: bool
    hover_fields: list[str]


class EICMultiPlotPayload(TypedDict):
    plot_schema: str
    plot_type: str
    title: str
    source: MultiPlotSource
    axes: PlotAxes
    sample: PlotSample
    normalization: str
    series: list[MultiPlotSeries]
    selection: PlotSelection
    render_hints: MultiRenderHints
    caveats: list[str]


def _axis_definition(main_type: int) -> tuple[str, str | None]:
    if main_type == 0:
        return "RT", "min"
    return f"Coordinate (main_type={main_type})", None


def _sample_metadata(file_path: str | Path) -> tuple[dict[int, dict], list[str]]:
    caveats = []
    try:
        mddata_path = resolve_mddata_path(file_path)
        if mddata_path is None:
            return {}, ["Adjacent .mddata was not found; series labels use FileID."]
        records = parse_analysis_file_classes(mddata_path)
        return {
            record["file_id"]: record
            for record in records
            if isinstance(record.get("file_id"), int)
        }, caveats
    except Exception as exc:
        return {}, [f"Sample metadata could not be loaded; series labels use FileID: {exc}"]


def _series_xy(sample: dict, normalize: str) -> tuple[list[float], list[float]]:
    """クロマトグラム点列を x/y 配列へ分解し、必要なら trace 内最大で正規化する。"""
    x_values = [float(point[0]) for point in sample["chromatogram"]]
    raw_y = [float(point[1]) for point in sample["chromatogram"]]
    if normalize == "per_trace_max":
        denominator = max(raw_y, default=0.0)
        return x_values, [
            value / denominator if denominator > 0 else 0.0 for value in raw_y
        ]
    return x_values, raw_y


def _apex_point(
    x_values: list[float], y_values: list[float], peak_top: float,
) -> tuple[float, float]:
    """peak_top に最も近いデータ点の (x, y) を返す。点列が空なら (peak_top, 0.0)。"""
    if not x_values:
        return float(peak_top), 0.0
    index = min(
        range(len(x_values)), key=lambda position: abs(x_values[position] - peak_top)
    )
    return x_values[index], y_values[index]


def build_eic_plot_payload(
    spot: dict,
    file_path: str | Path,
    *,
    normalize: str = "none",
    title: str | None = None,
) -> EICPlotPayload:
    """Build a renderer-neutral line-plot payload from selected EIC traces."""
    if normalize not in {"none", "per_trace_max"}:
        raise ValueError("normalize must be 'none' or 'per_trace_max'")

    path = Path(file_path).resolve()
    metadata, caveats = _sample_metadata(path)
    x_label, x_unit = _axis_definition(int(spot["main_type"]))
    series = []
    for sample in spot["samples"]:
        file_id = int(sample["file_id"])
        record = metadata.get(file_id, {})
        sample_name = record.get("file_name")
        class_id = record.get("class_id")
        label = sample_name or f"FileID {file_id}"
        x_values, y_values = _series_xy(sample, normalize)
        series.append({
            "id": f"file-{file_id}",
            "label": str(label),
            "file_id": file_id,
            "sample_name": str(sample_name) if sample_name else None,
            "class_id": str(class_id) if class_id else None,
            "x": x_values,
            "y": y_values,
            "peak_left": float(sample["peak_left"]),
            "peak_top": float(sample["peak_top"]),
            "peak_right": float(sample["peak_right"]),
            "max_intensity": float(sample["max_intensity"]),
            "mean_intensity": float(sample["mean_intensity"]),
            "point_count": int(sample["num_points"]),
        })

    y_label = "Relative intensity" if normalize == "per_trace_max" else "Intensity"
    plot_title = title or (
        f"EIC spot {spot['spot_id']} | m/z {spot['mz']:.4f} | RT {spot['rt']:.2f} min"
    )
    return round_floats({
        "plot_schema": "lipidmix.eic.v1",
        "plot_type": "line",
        "title": plot_title,
        "source": {"file": str(path), "file_name": path.name},
        "axes": {
            "x": {"label": x_label, "unit": x_unit, "scale": "linear"},
            "y": {"label": y_label, "unit": None, "scale": "linear"},
        },
        "spot": {
            "spot_id": int(spot["spot_id"]),
            "rt": float(spot["rt"]),
            "ri": float(spot["ri"]),
            "mz": float(spot["mz"]),
            "drift": float(spot["drift"]),
            "main_type": int(spot["main_type"]),
            "total_samples": int(spot["num_samples"]),
            "selected_samples": int(spot["selected_samples"]),
        },
        "normalization": normalize,
        "series": series,
        "render_hints": {
            "mode": "lines",
            "connect_points": True,
            "show_legend": True,
            "hover_fields": [
                "label", "file_id", "class_id", "peak_left", "peak_top",
                "peak_right", "max_intensity",
            ],
        },
        "caveats": caveats,
    }, POINT_DIGITS)


def _dropped(candidate: IdentityCandidate, reason: str) -> DroppedCompound:
    return {
        "spot_id": int(candidate["spot_id"]),
        "name": str(candidate["name"]),
        "reason": reason,
    }


_NAMED_DROP_REASONS = ("rt_mismatch", "mz_mismatch")
_AGGREGATE_ONLY_DROP_REASONS = ("spot_out_of_range", "file_id_absent", "below_top_n")
_NAMED_DROP_CAP = 10


def _dropped_caveats(dropped: list[DroppedCompound]) -> list[str]:
    """`dropped[]` を理由ごとに1行へ集約した caveat 文を返す。

    `rt_mismatch` / `mz_mismatch` は個々の物質名が診断に有用なので、最大
    `_NAMED_DROP_CAP` 件まで名指しし、残りは件数のみ添える。それ以外の理由
    （`spot_out_of_range` / `file_id_absent` / `below_top_n`）は件数のみの
    1行に集約する。`selection.dropped` 自体は完全なまま変更しない。
    """
    if not dropped:
        return []

    by_reason: dict[str, list[DroppedCompound]] = {}
    for item in dropped:
        by_reason.setdefault(item["reason"], []).append(item)

    notes: list[str] = []
    for reason in _NAMED_DROP_REASONS:
        items = by_reason.get(reason)
        if not items:
            continue
        names = [
            f"spot_id={item['spot_id']} ({item['name'] or 'Unknown'})"
            for item in items[:_NAMED_DROP_CAP]
        ]
        note = f"{reason} で {len(items)} 件を除外しました: " + ", ".join(names)
        if len(items) > _NAMED_DROP_CAP:
            note += f" ほか {len(items) - _NAMED_DROP_CAP} 件（詳細は selection.dropped）。"
        else:
            note += "。"
        notes.append(note)

    for reason in _AGGREGATE_ONLY_DROP_REASONS:
        items = by_reason.get(reason)
        if not items:
            continue
        notes.append(
            f"{reason} で {len(items)} 件を除外しました（詳細は selection.dropped）。"
        )

    return notes


def build_multi_compound_plot_payload(
    candidates: list[IdentityCandidate],
    spots: list[dict],
    *,
    file_id: int,
    file_path: str | Path,
    arf2_path: str | Path,
    normalize: str = "none",
    top_n: int = 24,
    title: str | None = None,
    queries: list[str] | None = None,
    ontologies: list[str] | None = None,
    caveats: list[str] | None = None,
    total_matched: int | None = None,
) -> EICMultiPlotPayload:
    """1 サンプル分の複数物質オーバーレイ用ペイロードを組み立てる。

    `candidates` は ARF2 由来の同定候補、`spots` は同じ `spot_id` を要求して読んだ
    EIC スポット。rt/mz 検証、`top_n` での強度打ち切り、RT 昇順の並べ替えを行い、
    除外された物質は理由付きで `selection.dropped` に残す。

    `total_matched` は `select_identity_candidates` の `max_candidates` 予備選抜が
    行われる前のクエリ一致件数。省略時は `len(candidates)`（予備選抜なし）とみなす。
    `selection.candidates` にはこの値を、実際に読み出し・検証した件数は
    `selection.candidates_evaluated`（= `len(candidates)`）に入れる。
    """
    if normalize not in {"none", "per_trace_max"}:
        raise ValueError("normalize must be 'none' or 'per_trace_max'")
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")

    path = Path(file_path).resolve()
    arf2 = Path(arf2_path).resolve()
    metadata, meta_caveats = _sample_metadata(path)
    notes = list(caveats or []) + meta_caveats
    dropped: list[DroppedCompound] = []

    spot_by_id = {int(spot["spot_id"]): spot for spot in spots}
    accepted = []
    for candidate in candidates:
        spot = spot_by_id.get(int(candidate["spot_id"]))
        if spot is None:
            dropped.append(_dropped(candidate, "spot_out_of_range"))
            continue
        sample = next(
            (item for item in spot["samples"] if int(item["file_id"]) == int(file_id)),
            None,
        )
        if sample is None:
            dropped.append(_dropped(candidate, "file_id_absent"))
            continue
        reason = verify_spot_match(candidate, spot)
        if reason:
            dropped.append(_dropped(candidate, reason))
            continue
        accepted.append((candidate, spot, sample))

    accepted.sort(key=lambda item: float(item[2]["max_intensity"]), reverse=True)
    for candidate, _spot, _sample in accepted[top_n:]:
        dropped.append(_dropped(candidate, "below_top_n"))
    accepted = accepted[:top_n]
    accepted.sort(key=lambda item: float(item[1]["rt"]))

    if not accepted:
        reason_counts: dict[str, int] = {}
        for item in dropped:
            reason_counts[item["reason"]] = reason_counts.get(item["reason"], 0) + 1
        breakdown = ", ".join(
            f"{reason}={count}" for reason, count in reason_counts.items()
        )
        examples = ", ".join(
            f"{item['name'] or item['spot_id']}={item['reason']}" for item in dropped[:5]
        )
        example_note = f"（例: {examples}）" if examples else ""
        raise ValueError(
            f"検証を通過した物質がありません（候補 {len(candidates)} 件、除外内訳: "
            f"{breakdown}{example_note}）。クエリ、file_id、または対象バッチを見直してください。"
        )

    series: list[MultiPlotSeries] = []
    main_types: list[int] = []
    for candidate, spot, sample in accepted:
        x_values, y_values = _series_xy(sample, normalize)
        label = str(candidate["name"] or "")
        if not label or label.casefold() == "unknown":
            label = f"spot {int(candidate['spot_id'])} (m/z {float(spot['mz']):.4f})"
        apex_x, apex_y = _apex_point(x_values, y_values, float(sample["peak_top"]))
        main_types.append(int(spot["main_type"]))
        series.append({
            "id": f"spot-{int(candidate['spot_id'])}",
            "label": label,
            "spot_id": int(candidate["spot_id"]),
            "name": str(candidate["name"]),
            "ontology": str(candidate["ontology"]),
            "adduct": str(candidate["adduct"]),
            "mz": float(spot["mz"]),
            "rt": float(spot["rt"]),
            "x": x_values,
            "y": y_values,
            "peak_left": float(sample["peak_left"]),
            "peak_top": float(sample["peak_top"]),
            "peak_right": float(sample["peak_right"]),
            "max_intensity": float(sample["max_intensity"]),
            "mean_intensity": float(sample["mean_intensity"]),
            "point_count": int(sample["num_points"]),
            "annotation": {
                "text": f"{label} / {float(sample['peak_top']):.3f}",
                "x": apex_x,
                "y": apex_y,
            },
        })

    notes.extend(_dropped_caveats(dropped))

    main_type = max(set(main_types), key=main_types.count)
    if len(set(main_types)) > 1:
        notes.append(
            f"選択スポットの main_type が混在しています（採用: {main_type}）。"
        )
    x_label, x_unit = _axis_definition(main_type)
    y_label = "Relative intensity" if normalize == "per_trace_max" else "Intensity"
    record = metadata.get(int(file_id), {})
    sample_name = record.get("file_name")
    class_id = record.get("class_id")
    sample_label = str(sample_name) if sample_name else f"FileID {int(file_id)}"
    plot_title = title or f"EIC overlay | {len(series)} compounds | {sample_label}"

    return round_floats({
        "plot_schema": MULTI_PLOT_SCHEMA,
        "plot_type": "line",
        "title": plot_title,
        "source": {
            "file": str(path), "file_name": path.name, "arf2_file": str(arf2),
        },
        "axes": {
            "x": {"label": x_label, "unit": x_unit, "scale": "linear"},
            "y": {"label": y_label, "unit": None, "scale": "linear"},
        },
        "sample": {
            "file_id": int(file_id),
            "sample_name": str(sample_name) if sample_name else None,
            "class_id": str(class_id) if class_id else None,
        },
        "normalization": normalize,
        "series": series,
        "selection": {
            "queries": [str(item) for item in (queries or [])],
            "ontologies": [str(item) for item in (ontologies or [])],
            "candidates": (
                total_matched if total_matched is not None else len(candidates)
            ),
            "candidates_evaluated": len(candidates),
            "plotted": len(series),
            "top_n": top_n,
            "dropped": dropped,
        },
        "render_hints": {
            "mode": "lines",
            "connect_points": True,
            "show_legend": True,
            "show_annotations": True,
            "hover_fields": [
                "label", "spot_id", "ontology", "adduct", "mz", "rt",
                "peak_left", "peak_top", "peak_right", "max_intensity",
            ],
        },
        "caveats": notes,
    }, POINT_DIGITS)


def render_eic_plot(payload, title: str | None = None):
    """Render an approved EIC payload to a matplotlib figure for explicit saving."""
    schema = payload.get("plot_schema")
    if schema == MULTI_PLOT_SCHEMA:
        return _render_multi_compound(payload, title)
    if schema == "lipidmix.eic.v1":
        return _render_single_spot(payload, title)
    raise ValueError(f"Unsupported EIC plot schema: {schema!r}")


def _render_single_spot(payload: EICPlotPayload, title: str | None = None):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 6))
    for item in payload["series"]:
        line = ax.plot(item["x"], item["y"], label=item["label"], linewidth=1.4)[0]
        if item["x"]:
            apex_x, apex_y = _apex_point(item["x"], item["y"], item["peak_top"])
            ax.scatter([apex_x], [apex_y], color=line.get_color(), s=20, zorder=3)
    if len(payload["series"]) == 1:
        item = payload["series"][0]
        ax.axvspan(item["peak_left"], item["peak_right"], color="#999999", alpha=0.12)
    _apply_axes(ax, payload, title)
    ax.legend(fontsize=8, loc="best")
    return fig


def _render_multi_compound(payload: EICMultiPlotPayload, title: str | None = None):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 6))
    show_annotations = payload["render_hints"].get("show_annotations", True)
    for item in payload["series"]:
        line = ax.plot(item["x"], item["y"], label=item["label"], linewidth=1.2)[0]
        annotation = item.get("annotation")
        if not annotation:
            continue
        ax.scatter(
            [annotation["x"]], [annotation["y"]],
            color=line.get_color(), s=16, zorder=3,
        )
        if show_annotations:
            ax.annotate(
                annotation["text"],
                xy=(annotation["x"], annotation["y"]),
                xytext=(0, 6),
                textcoords="offset points",
                fontsize=6,
                rotation=45,
                ha="left",
                va="bottom",
            )
    _apply_axes(ax, payload, title)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    return fig


def _apply_axes(ax, payload, title: str | None) -> None:
    x_axis = payload["axes"]["x"]
    y_axis = payload["axes"]["y"]
    x_label = x_axis["label"] + (f" ({x_axis['unit']})" if x_axis["unit"] else "")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_axis["label"])
    ax.set_title(title or payload["title"])
    ax.grid(alpha=0.2)
