"""群別強度プロット: 選んだクラス・分子種ごとに、試料群の試料別強度を並べる（MCP 非依存）。

1 パネル = 1 項目、1 点 = 1 試料（項目に当たったスポットの PeakHeight の合計）、群ごとに
log10 の平均 ± SD。検定はしない（`arf_differential` の役割）。spec:
docs/superpowers/specs/2026-10-06-group-intensity-plot-design.md。
"""
from __future__ import annotations

import math
import re
import statistics

from metabolomix.msdial import sample_factors

GROUP_INTENSITY_SCHEMA = "lipidmix.group_intensity.v1"
MAX_ITEMS = 30
#: 重水素などで標識した内部標準（`PC 33:1(d7)`、`FA 16:0(d3)`）。クラスの合計に入れない。
INTERNAL_STANDARD = re.compile(r"\(d\d+\)")
#: 標準液での最大高さが、描く試料での最大高さのこの倍数以上なら「標準液にだけある」。
STANDARD_ONLY_FOLD = 10.0
#: 合計のうち gap-fill の値がこの割合を超える点は形を変える。
GAPFILL_MAJORITY = 0.5
_PREFIX = re.compile(r"^\s*(?:low score|no MS2|w/o MS2|unsettled)\s*:\s*", re.IGNORECASE)
_EXCLUDED_KEYS = ("internal_standard", "curation", "auto_likely_wrong", "standard_only")


def _candidate_names(name: str) -> set[str]:
    return {_PREFIX.sub("", part).strip().casefold() for part in (name or "").split("|") if part.strip()}


def resolve_items(items, catalog, *, curation=None, standard_only=frozenset()):
    """項目（クラス・名前、`+` で合算）をスポット集合に解決する。spec §3.1 / §5。"""
    items = [str(item) for item in (items or [])]
    if not 1 <= len(items) <= MAX_ITEMS:
        raise ValueError(f"items は 1〜{MAX_ITEMS} 件で指定してください（受け取った数: {len(items)}）。")
    curation = curation or {}
    wrong = set(curation.get("wrong") or ()) | set(curation.get("redundant") or ())
    auto = set(curation.get("auto_likely_wrong") or ())
    assign = curation.get("assign") or {}

    spots = []
    for row in catalog:
        sid = int(row["MasterAlignmentID"])
        name, ontology = row.get("Name") or "", (row.get("Ontology") or "").strip()
        if sid in assign:
            name = assign[sid].get("name") or name
            ontology = (assign[sid].get("ontology") or ontology).strip()
        spots.append({"spot_id": sid, "name": name, "ontology": ontology})
    ontologies = {s["ontology"].casefold() for s in spots if s["ontology"]}

    excluded = {key: set() for key in _EXCLUDED_KEYS}

    def keep(spot, *, as_class):
        label = f"#{spot['spot_id']} {spot['name']}"
        if spot["spot_id"] in wrong:
            excluded["curation"].add(label)
            return False
        if spot["spot_id"] in auto:
            excluded["auto_likely_wrong"].add(label)
            return False
        if as_class and INTERNAL_STANDARD.search(spot["name"]):
            excluded["internal_standard"].add(label)
            return False
        if as_class and spot["spot_id"] in standard_only:
            excluded["standard_only"].add(label)
            return False
        return True

    resolved = []
    for item in items:
        parts = [p.strip() for p in item.split("+") if p.strip()]
        if not parts:
            raise ValueError(f"空の項目は指定できません: {item!r}")
        chosen: dict[int, dict] = {}
        out_parts = []
        for part in parts:
            key = part.casefold()
            if key in ontologies:
                hits = [s for s in spots if s["ontology"].casefold() == key and keep(s, as_class=True)]
                kind = "class"
            else:
                hits = [s for s in spots if key in _candidate_names(s["name"]) and keep(s, as_class=False)]
                kind = "name" if hits else "none"
            for s in hits:
                chosen[s["spot_id"]] = s
            out_parts.append({"part": part, "kind": kind, "n_spots": len(hits)})
        resolved.append({"item": item, "parts": out_parts,
                         "spots": [chosen[k] for k in sorted(chosen)]})
    return resolved, {key: sorted(values) for key, values in excluded.items()}


def standard_only_spots(rows_by_spot, standard_samples, plotted_samples, fold=STANDARD_ONLY_FOLD):
    """標準液の試料にだけ実質的な強度があるスポット。spec §5-4。"""
    found = set()
    for sid, rows in rows_by_spot.items():
        std = max((float(r.get("height") or 0.0) for r in rows if r.get("file_name") in standard_samples),
                  default=0.0)
        grp = max((float(r.get("height") or 0.0) for r in rows if r.get("file_name") in plotted_samples),
                  default=0.0)
        if std > 0 and std >= fold * max(grp, 1.0):
            found.add(sid)
    return frozenset(found)


def resolve_groups(group_specs, facets):
    """群の指定を試料名に解決する（指定順）。spec §3.2。"""
    specs = [str(s) for s in (group_specs or [])]
    if not specs:
        raise ValueError("groups を 1 群以上指定してください。")
    groups = []
    for spec in specs:
        tokens = sample_factors.split_tokens(spec)
        if not tokens:
            raise ValueError(f"空の群指定は使えません: {spec!r}")
        roles = ["sample"] + [role for role in ("blank", "qc") if role in tokens]
        matches, _ = sample_factors.expand_sample_specs([spec], facets, include_roles=tuple(roles))
        groups.append({"label": spec, "samples": list(matches[spec])})
    caveats = []
    seen: dict[str, str] = {}
    for group in groups:
        for name in group["samples"]:
            if name in seen and seen[name] != group["label"]:
                caveats.append(f"試料 {name} は群 {seen[name]} と {group['label']} の両方に当たり、両方に描きました。")
            seen.setdefault(name, group["label"])
    return groups, caveats


def resolve_sample_specs(specs, facets):
    """試料名（完全一致）またはトークン指定を試料名の集合にする（役割で絞らない）。"""
    found = set()
    for spec in specs or []:
        if spec in facets:
            found.add(spec)
            continue
        matches, _ = sample_factors.expand_sample_specs([spec], facets, include_roles=None)
        found.update(matches.get(spec, []))
    return found


def build_group_intensity_payload(resolved_items, groups, rows_by_spot, *, msms=None,
                                  low_reliability=frozenset(), excluded=None,
                                  detection_limit=None, detection_limit_source=None, caveats=None):
    """項目 × 群 × 試料のクラス合計強度を組み立てる（描画しない）。spec §4 / §7。"""
    index = {sid: {r.get("file_name"): r for r in rows} for sid, rows in rows_by_spot.items()}
    out_items = []
    for item in resolved_items:
        sids = [s["spot_id"] for s in item["spots"]]
        out_groups = []
        for group in groups:
            samples = []
            for name in group["samples"]:
                total = gap = 0.0
                for sid in sids:
                    row = index.get(sid, {}).get(name)
                    if row is None:
                        continue
                    h = float(row.get("height") or 0.0)
                    total += h
                    if row.get("is_gap_filled"):
                        gap += h
                samples.append({"sample": name, "value": round(total, 1),
                                "gap_filled_fraction": round(gap / total, 3) if total > 0 else None,
                                "low_reliability": name in low_reliability})
            logs = [math.log10(s["value"]) for s in samples if s["value"] > 0 and not s["low_reliability"]]
            out_groups.append({
                "label": group["label"], "samples": samples,
                "log10_mean": round(statistics.mean(logs), 4) if logs else None,
                "log10_sd": round(statistics.stdev(logs), 4) if len(logs) > 1 else None,
                "n_in_stats": len(logs),
            })
        out_items.append({
            "item": item["item"], "parts": item["parts"], "spots": item["spots"],
            "n_spots": len(sids),
            "n_spots_msms": sum(1 for sid in sids if (msms or {}).get(sid)) if msms is not None else None,
            "detected": any(s["value"] > 0 for g in out_groups for s in g["samples"]),
            "groups": out_groups,
        })
    return {
        "plot_schema": GROUP_INTENSITY_SCHEMA,
        "value": "sum of PeakHeight (log10 on the plot)",
        "detection_limit": detection_limit, "detection_limit_source": detection_limit_source,
        "groups": [{"label": g["label"], "samples": list(g["samples"])} for g in groups],
        "low_reliability_samples": sorted(low_reliability),
        "items": out_items,
        "excluded": {key: list((excluded or {}).get(key, [])) for key in _EXCLUDED_KEYS},
        "caveats": list(caveats or []),
    }


GROUP_COLORS = ["#bdbdbd", "#9ecae1", "#3182bd", "#de2d26", "#fd8d3c", "#31a354", "#756bb1", "#8c564b"]


def render_group_intensity_plot(payload, *, title=None, ncols=None):
    """payload を 1 項目 1 パネルの点 + 平均 ± SD 図にする。spec §4。

    縦軸は全パネル共通（log10）。パネルごとに拡大すると強度 10 前後のノイズが信号に見え、
    「見つかるか」の判断を誤らせるため。0 の試料は軸の底に ▽。
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    items = payload["items"]
    if not items:
        raise ValueError("描く項目がありません。")
    labels = [g["label"] for g in payload["groups"]]
    ncols = ncols or min(5, len(items))
    nrows = math.ceil(len(items) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 2.7 * nrows), squeeze=False)
    positive = [s["value"] for it in items for g in it["groups"] for s in g["samples"] if s["value"] > 0]
    floor = 0.0
    top = math.ceil(math.log10(max(positive)) + 0.2) if positive else 1.0
    limit = payload.get("detection_limit")
    for ax, item in zip(axes.flat, items):
        ms1_only = item["n_spots"] > 0 and item.get("n_spots_msms") == 0
        msms_note = "" if item.get("n_spots_msms") is None else f", MS/MS {item['n_spots_msms']}"
        ax.set_title(f"{item['item']}  (spots {item['n_spots']}{msms_note})", fontsize=9,
                     color="#777777" if ms1_only else "black", style="italic" if ms1_only else "normal")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
        ax.set_xlim(-0.6, len(labels) - 0.4)
        ax.set_ylim(floor - 0.3, top)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("log10(Intensity)", fontsize=7)
        if limit:
            ax.axhline(math.log10(limit), color="#888888", linestyle="--", linewidth=0.8, zorder=1)
        if ms1_only:
            ax.set_facecolor("#f4f4f4")
            ax.text(0.03, 0.95, "MS1-only (unconfirmed)", transform=ax.transAxes, ha="left", va="top",
                    fontsize=7, color="#777777", style="italic")
        if not item["detected"]:
            reason = "no annotated species" if item["n_spots"] == 0 else "intensity 0 in all samples"
            ax.text(0.5, 0.5, f"N.D.\n({reason})", transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="#555555",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 2}, zorder=5)
            continue
        for i, group in enumerate(item["groups"]):
            color = GROUP_COLORS[i % len(GROUP_COLORS)]
            n = len(group["samples"])
            for j, s in enumerate(group["samples"]):
                x = i + (j - (n - 1) / 2) * 0.12
                if s["value"] <= 0:
                    ax.scatter(x, floor, marker="v", s=22, facecolor="white", edgecolor="#555555", zorder=3)
                    continue
                gapfilled = (s["gap_filled_fraction"] or 0) > GAPFILL_MAJORITY
                ax.scatter(x, math.log10(s["value"]), s=30, marker="D" if gapfilled else "o",
                           facecolor="white" if s["low_reliability"] else color,
                           edgecolor="black", linewidth=0.6, zorder=3)
            if group["log10_mean"] is not None:
                m, sd = group["log10_mean"], group["log10_sd"]
                ax.hlines(m, i - 0.28, i + 0.28, color="black", linewidth=1.6, zorder=4)
                if sd is not None:
                    ax.errorbar(i, m, yerr=sd, color="black", capsize=4, linewidth=1.0, zorder=4)
        ax.axhline(floor, color="#dddddd", linewidth=0.6, zorder=1)
    for ax in list(axes.flat)[len(items):]:
        ax.axis("off")
    legend = [
        Line2D([], [], marker="o", ls="", mfc="#888888", mec="black", label="sample (mostly detected peaks)"),
        Line2D([], [], marker="D", ls="", mfc="#888888", mec="black", label="sample (mostly gap-filled values)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="black", label="low-reliability sample (not in mean/SD)"),
        Line2D([], [], marker="v", ls="", mfc="white", mec="#555555", label="zero (at the axis floor)"),
        Line2D([], [], color="black", lw=1.6, label="mean ± SD of log10"),
    ]
    if limit:
        legend.append(Line2D([], [], color="#888888", ls="--", lw=0.8,
                             label=f"detection limit ({limit:,.0f})"))
    fig.legend(handles=legend, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.suptitle(title or "Intensity by sample group", fontsize=11)
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    return fig
