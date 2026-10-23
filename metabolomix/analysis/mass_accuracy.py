"""m/z 帯ごとの質量精度と、新しい知見の候補の一覧（spec 2026-10-22 §4 の 9(c)・11）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

- 質量精度: SME（MS/MS の照合で付いた同定）の実測 m/z と理論 m/z の差（ppm）の分布。
  許容幅の決定はしない（再解析の規則で決める）。四分位は `statistics.quantiles(n=4)`
  （exclusive 法）。
- 新しい知見の候補: SME の同定・m/z 帯ごとの許容幅・q の閾値で絞り、除く対象を外した
  特徴の一覧。統計は 2 群比較の結果（`result_id`）か analysis-spec の解析から引く。
"""
from __future__ import annotations

import math
import statistics

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["NOVEL_SCHEMA", "band_summary", "normalize_novel_spec", "novel_candidates",
           "sme_ppm", "summarize"]

NOVEL_SCHEMA = "novel-candidates.v1"
_INVALID = "NOVEL_SPEC_INVALID"


def _number(value) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return float(value)


def sme_ppm(ds, rank: int = 1) -> list[dict]:
    """特徴ごとに、その rank の SME 候補の (実測 − 理論)/理論 × 1e6。並びは ds.feature_ids。"""
    out = []
    for fid in ds.feature_ids:
        for candidate in (ds.feature_candidates or {}).get(fid) or []:
            if candidate.get("rank") != rank:
                continue
            exp = _number(candidate.get("exp_mass_to_charge"))
            theo = _number(candidate.get("theoretical_mass_to_charge"))
            if exp is None or theo is None or theo <= 0:
                continue
            out.append({"feature_id": fid, "theoretical_mz": theo, "mz": exp,
                        "ppm": (exp - theo) / theo * 1e6, "name": candidate.get("chemical_name"),
                        "inchikey": candidate.get("inchikey"), "adduct": candidate.get("adduct")})
            break
    return out


def summarize(values: list[float]) -> dict:
    values = sorted(values)
    if len(values) < 2:
        return {"n": len(values)}
    q1, _, q3 = statistics.quantiles(values, n=4)
    return {"n": len(values), "median": statistics.median(values), "q1": q1, "q3": q3,
            "within_5ppm": sum(abs(v) <= 5 for v in values) / len(values),
            "within_10ppm": sum(abs(v) <= 10 for v in values) / len(values)}


def band_summary(rows: list[dict], band_edges: list[float]) -> dict:
    edges = sorted(band_edges)
    bounds = list(zip([None, *edges], [*edges, None]))
    bands = []
    for low, high in bounds:
        chosen = [r["ppm"] for r in rows
                  if (low is None or r["theoretical_mz"] >= low)
                  and (high is None or r["theoretical_mz"] < high)]
        bands.append({"min_mz": low, "max_mz": high, **summarize(chosen)})
    return {"overall": summarize([r["ppm"] for r in rows]), "bands": bands}


# ---------- 新しい知見の候補 ----------

def normalize_novel_spec(spec) -> dict:
    if not isinstance(spec, dict) or set(spec) - {"schema", "statistics", "ppm", "q_max", "exclude"}:
        raise DomainError(_INVALID, "spec は {statistics, ppm, q_max, exclude}", {})
    stats_spec = spec.get("statistics")
    if not isinstance(stats_spec, dict):
        raise DomainError(_INVALID, "statistics は {result_id} か {analysis_spec, analysis_id}",
                          {"field": "statistics"})
    if set(stats_spec) == {"result_id"}:
        statistics_out = {"result_id": str(stats_spec["result_id"])}
    elif set(stats_spec) == {"analysis_spec", "analysis_id"}:
        from metabolomix.analysis import analysis_spec
        statistics_out = {"analysis_spec": analysis_spec.normalize_spec(stats_spec["analysis_spec"]),
                          "analysis_id": str(stats_spec["analysis_id"])}
        ids = [a["analysis_id"] for a in statistics_out["analysis_spec"]["analyses"]]
        if statistics_out["analysis_id"] not in ids:
            raise DomainError(_INVALID, f"analysis_id が spec に無い（{ids}）",
                              {"field": "statistics.analysis_id"})
    else:
        raise DomainError(_INVALID, "statistics は {result_id} か {analysis_spec, analysis_id}",
                          {"field": "statistics"})
    ppm = spec.get("ppm", {"default": 10.0, "bands": []})
    if (not isinstance(ppm, dict) or set(ppm) - {"default", "bands"}
            or _number(ppm.get("default", 10.0)) is None):
        raise DomainError(_INVALID, "ppm は {default, bands: [{max_mz, ppm}]}", {"field": "ppm"})
    bands = []
    for i, band in enumerate(ppm.get("bands") or []):
        if (not isinstance(band, dict) or set(band) != {"max_mz", "ppm"}
                or _number(band["max_mz"]) is None or _number(band["ppm"]) is None):
            raise DomainError(_INVALID, f"ppm.bands[{i}] は {{max_mz, ppm}}", {"field": "ppm"})
        bands.append({"max_mz": float(band["max_mz"]), "ppm": float(band["ppm"])})
    q_max = _number(spec.get("q_max", 0.05))
    if q_max is None or not 0 < q_max <= 1:
        raise DomainError(_INVALID, "q_max は 0 < q_max ≤ 1", {"field": "q_max"})
    exclude = spec.get("exclude") or {}
    if not isinstance(exclude, dict) or set(exclude) - {"feature_ids", "inchikey_blocks"}:
        raise DomainError(_INVALID, "exclude は {feature_ids, inchikey_blocks}", {"field": "exclude"})
    return {"schema": NOVEL_SCHEMA, "statistics": statistics_out,
            "ppm": {"default": float(ppm.get("default", 10.0)),
                    "bands": sorted(bands, key=lambda b: b["max_mz"])},
            "q_max": q_max,
            "exclude": {"feature_ids": sorted(set(map(str, exclude.get("feature_ids") or []))),
                        "inchikey_blocks": sorted({str(b).strip()[:14].upper()
                                                   for b in exclude.get("inchikey_blocks") or []})}}


def rule_hash(normalized: dict) -> str:
    from metabolomix.analysis.analysis_spec import _strip_paths

    return canonical_hash(_strip_paths(normalized))


def _tolerance(mz: float, ppm_spec: dict) -> float:
    for band in ppm_spec["bands"]:
        if mz < band["max_mz"]:
            return band["ppm"]
    return ppm_spec["default"]


def novel_candidates(ds, normalized: dict, stats_rows: dict, *, group_ref, group_test) -> list[dict]:
    """`stats_rows` は feature_id → {log2fc_test_over_ref, p, q}（検定できた行）。"""
    excluded_ids = set(normalized["exclude"]["feature_ids"])
    excluded_blocks = set(normalized["exclude"]["inchikey_blocks"])
    out = []
    for row in sme_ppm(ds, rank=1):
        fid = row["feature_id"]
        stats = stats_rows.get(fid)
        if stats is None or stats.get("q") is None or not stats["q"] < normalized["q_max"]:
            continue
        if abs(row["ppm"]) > _tolerance(row["theoretical_mz"], normalized["ppm"]):
            continue
        block = (row["inchikey"] or "").strip()[:14].upper()
        if fid in excluded_ids or (block and block in excluded_blocks):
            continue
        log2fc = stats.get("log2fc_test_over_ref")
        out.append({**row, "rt": (ds.feature_metadata.get(fid) or {}).get("rt"),
                    "group_ref": group_ref, "group_test": group_test,
                    "log2fc_test_over_ref": log2fc, "p": stats.get("p"), "q": stats["q"],
                    "direction": None if not log2fc else
                    ("test_higher" if log2fc > 0 else "ref_higher")})
    out.sort(key=lambda r: (r["q"], r["feature_id"]))
    return out
