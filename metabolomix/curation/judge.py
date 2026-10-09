"""機械判別（純関数）。spec §4。

各系統を PASS / BORDERLINE / FAIL / UNKNOWN の帯と理由コードにし、総合判定を出す。
- UNKNOWN は FAIL に数えない（MS/MS 未取得・参照 RT なしは「合わなかった」ではない）。
- MS-DIAL 自身がその解析の許容幅で出した判定（IsPrecursorMzMatch / IsReferenceMatched）を
  そのまま使う。ppm と ΔRT の帯はこちらの固定既定値で、引数で上書きできる。
- ⑤傾向は単独で総合判定を上げない。①〜④の BORDERLINE と重なったときだけ補強する。
- MS-DIAL の脂質規則フラグ（IsLipidClassMatch / IsLipidChainsMatch / IsOtherLipidMatch）は
  Lipidomics 採点器でしか立たない。メタボロミクス採点器では全部 False なので、そのレビューで
  規則が走った証拠（`lipid_rules_active`）があるときだけ読む（`lipid_rules=True`）。
"""
from __future__ import annotations

import math
import re

DEFAULT_THRESHOLDS: dict[str, float] = {
    "ppm_pass": 5.0, "ppm_borderline": 10.0, "dmz_fail_mda": 10.0,
    "drt_pass": 0.5, "drt_borderline": 1.0,
    "eic_min_points": 5, "eic_min_r2": 0.8, "eic_max_maxima": 2,
    "eic_pass_frac": 0.5, "eic_borderline_frac": 0.2, "eic_rt_scatter_sd": 0.1,
    "trend_min_points": 5, "trend_outlier_z": 3.0, "trend_min_r2": 0.7,
    "rescore_tolerance": 0.1,
    "adduct_isomer_ppm": 5.0, "adduct_isomer_drt": 0.05, "adduct_isomer_height_ratio": 3.0,
}

# ppm_out は adduct 非依存(実測: 全 adduct で中央値 約 -0.8 ppm)で、precursor_unmatched=0
# の実データでも likely_wrong が ppm_out 単独からしか出ていなかった。MS-DIAL 自身が
# その解析の許容幅で判定した precursor_unmatched / polarity_mismatch とは信頼度が違うため、
# 弱い理由へ格下げする(ユーザー決定 2026-09-29)。
# class_rule_rejected: MS/MS ありで脂質クラス規則を評価して棄却(class=F ∧ other=F)。
# 診断イオン規則が注釈のクラスを否定したので強い理由にする(ユーザー承認 2026-09-29)。
# dmz_out: |Δm/z| が dmz_fail_mda（既定 10 mDa）以上。ppm は m/z に比例して緩むので、絶対差で
# 強い理由にする（ユーザー決定 2026-09-29）。Δppm>10 の ppm_out は弱いまま。
# adduct_isomer_of / adduct_isomer_minor_of: m/z が同時溶出する別物質のスポットの別アダクトで説明でき、
# 相手の証拠の段階が上なら強い・同じ段階で強度が 3 倍以上なら弱い（ユーザー承認 2026-10-07。
# 判定は adduct_isomer.py）。理由には相手のスポット番号を `:<spot>` で付ける。
STRONG_REASONS = frozenset({"polarity_mismatch", "precursor_unmatched", "class_rule_rejected",
                            "dmz_out", "adduct_isomer_of"})
WEAK_REASONS = frozenset({"ppm_out", "low_score", "drt_out", "eic_poor", "adduct_isomer_minor_of"})
_ORDER = ["polarity_mismatch", "precursor_unmatched", "dmz_out", "class_rule_rejected",
          "adduct_isomer_of", "ppm_out", "low_score", "drt_out", "eic_poor", "adduct_isomer_minor_of",
          "ppm_borderline", "drt_borderline", "eic_borderline", "rt_scatter", "trend_outlier"]


_RULE_FLAGS = ("is_lipid_class_match", "is_lipid_chains_match", "is_other_lipid_match")
#: 鎖を区切る `_` / `/` が鎖表記（`16:0` など）の直後に来る名前。`PC 34:1|PC 16:0_18:1` の
#: `|` 以降も拾う。単鎖（`LPC 16:0`）は種名と同じなので鎖レベルに数えない。
_CHAIN_LEVEL_RE = re.compile(r"\d+:\d+[^\s_/|]*[_/]")


def lipid_rules_active(evs) -> bool:
    """規則フラグが 1 件でも True なら、このデータは脂質規則で採点されている。"""
    return any(any((ev.get("match") or {}).get(k) for k in _RULE_FLAGS) for ev in evs)


def is_chain_level_name(name) -> bool:
    return bool(name) and _CHAIN_LEVEL_RE.search(str(name)) is not None


def resolve_thresholds(overrides: dict | None) -> dict:
    unknown = set(overrides or {}) - set(DEFAULT_THRESHOLDS)
    if unknown:
        raise ValueError(f"未知のしきい値キー: {sorted(unknown)}（有効: {sorted(DEFAULT_THRESHOLDS)}）")
    return {**DEFAULT_THRESHOLDS, **(overrides or {})}


def reason_code(reason: str) -> str:
    """理由から `:<spot>` の添え字を除いたコード（`adduct_isomer_of:189` → `adduct_isomer_of`）。"""
    return reason.split(":", 1)[0]


def _check(band: str, reasons: list[str] | None = None) -> dict:
    return {"band": band, "reasons": list(reasons or [])}


def _msms(ev: dict, info: list[str], th: dict, lipid_rules: bool) -> dict:
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
    class_rejected = False
    if lipid_rules:
        other = m.get("is_other_lipid_match")
        if m.get("is_lipid_class_match") is False:
            if other is True:
                info.append("class_rules_not_run")
            elif other is False:
                class_rejected = True
        if m.get("is_lipid_chains_match") is False and is_chain_level_name(ev.get("name")):
            info.append("chains_unsupported")
    if m.get("is_reference_matched"):
        return _check("PASS")
    return _check("FAIL", (["class_rule_rejected"] if class_rejected else []) + ["low_score"])


def _mz(ev: dict, th: dict) -> dict:
    reasons = []
    if ev.get("adduct_band") == "FAIL":
        reasons.append("polarity_mismatch")
    m = ev.get("match")
    if m is not None and m.get("is_precursor_mz_match") is False:
        reasons.append("precursor_unmatched")
    dmz = ev.get("dmz_mda")
    if dmz is not None and abs(dmz) >= th["dmz_fail_mda"]:
        reasons.append("dmz_out")
    ppm = ev.get("ppm")
    if ppm is not None:
        if abs(ppm) > th["ppm_borderline"]:
            reasons.append("ppm_out")
        elif abs(ppm) > th["ppm_pass"]:
            reasons.append("ppm_borderline")
    if any(r in STRONG_REASONS for r in reasons) or "ppm_out" in reasons:
        return _check("FAIL", reasons)     # ppm_out は弱いが、mz 系統自体は FAIL にする
    if reasons:
        return _check("BORDERLINE", reasons)
    return _check("UNKNOWN" if ppm is None else "PASS")


def _rt(ev: dict, info: list[str], th: dict) -> dict:
    if ev.get("reference") is None:
        info.append("reference_not_found")
        return _check("UNKNOWN")
    if ev.get("rt_used_by_annotation") is False:
        # 参照 RT は別のクロマトグラフィーの予測値でありうる（LBM 等）。MS-DIAL が
        # 使っていない RT のずれで誤りとは言えない。
        info.append("rt_not_used_by_annotation")
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


def _ion(ev: dict) -> dict:
    """別アダクトの取り違え。`adduct_isomer` が無ければ未評価（UNKNOWN）、None なら説明する相手が無い。"""
    if "adduct_isomer" not in ev:
        return _check("UNKNOWN")
    hit = ev["adduct_isomer"]
    if hit is None:
        return _check("PASS")
    code = "adduct_isomer_of" if hit["severity"] == "strong" else "adduct_isomer_minor_of"
    return _check("FAIL", [f"{code}:{hit['of']}"])


def _trend(entry: dict | None, info: list[str]) -> dict:
    if entry is None:
        return _check("UNKNOWN")
    if entry.get("outlier") and entry.get("reliable"):
        return _check("BORDERLINE", ["trend_outlier"])
    if entry.get("outlier"):
        info.append("trend_outlier_unreliable")
    return _check("PASS")


#: 理由コード → 判定根拠の文（ビューアのメモ欄の既定値）。値は ev と th から差し込む。
REASON_TEXT = {
    "polarity_mismatch": lambda ev, th: "Adduct — charge sign contradicts the ion mode",
    "precursor_unmatched": lambda ev, th: "m/z — MS-DIAL precursor match failed",
    "dmz_out": lambda ev, th: f"m/z — Δm/z {_num(ev.get('dmz_mda'))} mDa (≥{_num(th['dmz_fail_mda'])} mDa)",
    "class_rule_rejected": lambda ev, th: "MS2 — rejected by the lipid class rules (diagnostic ions)",
    "ppm_out": lambda ev, th: f"m/z — Δppm {_num(ev.get('ppm'))} (>{_num(th['ppm_borderline'])})",
    "low_score": lambda ev, th: "MS2 — no match to the reference (low score)",
    "drt_out": lambda ev, th: f"RT — ΔRT {_num(ev.get('drt'), 2)} min (>{_num(th['drt_borderline'])} min)",
    "eic_poor": lambda ev, th: "EIC — poor peak shape",
    "ppm_borderline": lambda ev, th: f"m/z — Δppm {_num(ev.get('ppm'))} (borderline)",
    "drt_borderline": lambda ev, th: f"RT — ΔRT {_num(ev.get('drt'), 2)} min (borderline)",
    "eic_borderline": lambda ev, th: "EIC — borderline peak shape",
    "rt_scatter": lambda ev, th: "EIC — apex RT scatters across samples",
    "trend_outlier": lambda ev, th: "RT–m/z trend — off the class trend",
    "adduct_isomer_of": lambda ev, th: "Adduct — " + _isomer_text(ev, "better-evidenced"),
    "adduct_isomer_minor_of": lambda ev, th: "Adduct — " + _isomer_text(ev, "more intense same-tier"),
}


def _isomer_text(ev: dict, kind: str) -> str:
    hit = ev.get("adduct_isomer") or {}
    return (f"explained as {hit.get('as_adduct')} of the co-eluting {kind} #{hit.get('of')} {hit.get('of_name')}"
            f" ({_num(hit.get('ppm'))} ppm, ΔRT {_num(hit.get('drt'), 3)} min)")


def _num(value, digits: int | None = None) -> str:
    if value is None:
        return "–"
    if digits is not None and isinstance(value, (int, float)):
        value = round(float(value), digits)
    return f"{value:g}" if isinstance(value, float) else str(value)


def auto_note(spot: dict, th: dict) -> str | None:
    """判定済みスポット（`judge_spot` の結果を持つ）の判定根拠。`ok` なら None。"""
    if spot.get("verdict") in (None, "ok") or not spot.get("reasons"):
        return None
    parts = [REASON_TEXT[reason_code(r)](spot, th) if reason_code(r) in REASON_TEXT else r
             for r in spot["reasons"]]
    return "Auto: " + " / ".join(parts)


def judge_spot(ev: dict, trend_entry: dict | None, th: dict, *, lipid_rules: bool = False) -> dict:
    info: list[str] = list(ev.get("notes") or [])
    ref_adduct = ev.get("reference_adduct")
    if ref_adduct and ev.get("adduct") and ref_adduct != ev.get("adduct"):
        info.append("adduct_differs_from_reference")
    checks = {"msms": _msms(ev, info, th, lipid_rules), "mz": _mz(ev, th), "rt": _rt(ev, info, th),
              "eic": _eic(ev, th), "ion": _ion(ev), "trend": _trend(trend_entry, info)}
    core = [checks[k] for k in ("msms", "mz", "rt", "eic", "ion")]
    reasons = [r for c in checks.values() for r in c["reasons"]]
    reasons.sort(key=lambda r: _ORDER.index(reason_code(r)) if reason_code(r) in _ORDER else len(_ORDER))

    if any(reason_code(r) in STRONG_REASONS for r in reasons):
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
