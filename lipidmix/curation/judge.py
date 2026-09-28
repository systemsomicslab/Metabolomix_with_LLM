"""機械判別（純関数）。spec §4。

各系統を PASS / BORDERLINE / FAIL / UNKNOWN の帯と理由コードにし、総合判定を出す。
- UNKNOWN は FAIL に数えない（MS/MS 未取得・参照 RT なしは「合わなかった」ではない）。
- MS-DIAL 自身がその解析の許容幅で出した判定（IsPrecursorMzMatch / IsReferenceMatched）を
  そのまま使う。ppm と ΔRT の帯はこちらの固定既定値で、引数で上書きできる。
- ⑤傾向は単独で総合判定を上げない。①〜④の BORDERLINE と重なったときだけ補強する。
"""
from __future__ import annotations

import math

DEFAULT_THRESHOLDS: dict[str, float] = {
    "ppm_pass": 5.0, "ppm_borderline": 10.0,
    "drt_pass": 0.5, "drt_borderline": 1.0,
    "eic_min_points": 5, "eic_min_r2": 0.8, "eic_max_maxima": 2,
    "eic_pass_frac": 0.5, "eic_borderline_frac": 0.2, "eic_rt_scatter_sd": 0.1,
    "trend_min_points": 5, "trend_outlier_z": 3.0, "trend_min_r2": 0.7,
    "rescore_tolerance": 0.1,
}

# ppm_out は adduct 非依存(実測: 全 adduct で中央値 約 -0.8 ppm)で、precursor_unmatched=0
# の実データでも likely_wrong が ppm_out 単独からしか出ていなかった。MS-DIAL 自身が
# その解析の許容幅で判定した precursor_unmatched / polarity_mismatch とは信頼度が違うため、
# 弱い理由へ格下げする(ユーザー決定 2026-09-29)。
STRONG_REASONS = frozenset({"polarity_mismatch", "precursor_unmatched"})
WEAK_REASONS = frozenset({"ppm_out", "low_score", "drt_out", "eic_poor"})
_ORDER = ["polarity_mismatch", "precursor_unmatched",
          "ppm_out", "low_score", "drt_out", "eic_poor",
          "ppm_borderline", "drt_borderline", "eic_borderline", "rt_scatter", "trend_outlier"]


def resolve_thresholds(overrides: dict | None) -> dict:
    unknown = set(overrides or {}) - set(DEFAULT_THRESHOLDS)
    if unknown:
        raise ValueError(f"未知のしきい値キー: {sorted(unknown)}（有効: {sorted(DEFAULT_THRESHOLDS)}）")
    return {**DEFAULT_THRESHOLDS, **(overrides or {})}


def _check(band: str, reasons: list[str] | None = None) -> dict:
    return {"band": band, "reasons": list(reasons or [])}


def _msms(ev: dict, info: list[str], th: dict) -> dict:
    m = ev.get("match")
    prefix = ev.get("name_prefix")
    if prefix == "unsettled":
        info.append("manually_unsettled")
    if m is None:
        info.append("no_match_result")
        return _check("UNKNOWN")
    if m.get("is_manually_modified"):
        info.append("manually_modified")
    rescore = ev.get("rescore")
    msdial_weighted = m.get("squared_weighted_dot_product")
    if rescore and msdial_weighted is not None and msdial_weighted >= 0:
        ours = rescore.get("weighted_dot_product")
        if ours is not None and ours >= 0 and abs(ours - math.sqrt(msdial_weighted)) > th["rescore_tolerance"]:
            info.append("rescore_discrepancy")
    if not m.get("has_msms") or prefix in ("no MS2", "w/o MS2"):
        info.append("msms_absent")
        return _check("UNKNOWN")
    if m.get("is_reference_matched"):
        return _check("PASS")
    return _check("FAIL", ["low_score"])


def _mz(ev: dict, th: dict) -> dict:
    reasons = []
    if ev.get("adduct_band") == "FAIL":
        reasons.append("polarity_mismatch")
    m = ev.get("match")
    if m is not None and m.get("is_precursor_mz_match") is False:
        reasons.append("precursor_unmatched")
    ppm = ev.get("ppm")
    if ppm is not None:
        if abs(ppm) > th["ppm_borderline"]:
            reasons.append("ppm_out")
        elif abs(ppm) > th["ppm_pass"]:
            reasons.append("ppm_borderline")
    if any(r in STRONG_REASONS or r in WEAK_REASONS for r in reasons):
        return _check("FAIL", reasons)     # ppm_out は弱いが、mz 系統自体は FAIL にする
    if reasons:
        return _check("BORDERLINE", reasons)
    return _check("UNKNOWN" if ppm is None else "PASS")


def _rt(ev: dict, info: list[str], th: dict) -> dict:
    if ev.get("reference") is None:
        info.append("reference_not_found")
        return _check("UNKNOWN")
    drt = ev.get("drt")
    if drt is None:
        info.append("reference_rt_absent")
        return _check("UNKNOWN")
    if abs(drt) > th["drt_borderline"]:
        return _check("FAIL", ["drt_out"])
    if abs(drt) > th["drt_pass"]:
        return _check("BORDERLINE", ["drt_borderline"])
    return _check("PASS")


def _eic(ev: dict, th: dict) -> dict:
    shape = ev.get("eic_shape") or {"band": "UNKNOWN"}
    band = shape.get("band", "UNKNOWN")
    reasons = {"FAIL": ["eic_poor"], "BORDERLINE": ["eic_borderline"]}.get(band, [])
    sd = shape.get("apex_rt_sd")
    if sd is not None and sd > th["eic_rt_scatter_sd"]:
        reasons.append("rt_scatter")
        if band == "PASS":
            band = "BORDERLINE"
    return _check(band, reasons)


def _trend(entry: dict | None, info: list[str]) -> dict:
    if entry is None:
        return _check("UNKNOWN")
    if entry.get("outlier") and entry.get("reliable"):
        return _check("BORDERLINE", ["trend_outlier"])
    if entry.get("outlier"):
        info.append("trend_outlier_unreliable")
    return _check("PASS")


def judge_spot(ev: dict, trend_entry: dict | None, th: dict) -> dict:
    info: list[str] = list(ev.get("notes") or [])
    ref_adduct = ev.get("reference_adduct")
    if ref_adduct and ev.get("adduct") and ref_adduct != ev.get("adduct"):
        info.append("adduct_differs_from_reference")
    checks = {"msms": _msms(ev, info, th), "mz": _mz(ev, th), "rt": _rt(ev, info, th),
              "eic": _eic(ev, th), "trend": _trend(trend_entry, info)}
    core = [checks[k] for k in ("msms", "mz", "rt", "eic")]
    reasons = [r for c in checks.values() for r in c["reasons"]]
    reasons.sort(key=lambda r: _ORDER.index(r) if r in _ORDER else len(_ORDER))

    if any(r in STRONG_REASONS for r in reasons):
        verdict = "likely_wrong"
    else:
        weak_fail = any(c["band"] == "FAIL" for c in core)
        borderlines = sum(1 for c in core if c["band"] == "BORDERLINE")
        trend_flag = checks["trend"]["band"] == "BORDERLINE"
        if weak_fail or borderlines >= 2 or (borderlines >= 1 and trend_flag):
            verdict = "suspect"
        else:
            verdict = "ok"
    return {"checks": checks, "verdict": verdict, "reasons": reasons,
            "info": sorted(set(info))}
