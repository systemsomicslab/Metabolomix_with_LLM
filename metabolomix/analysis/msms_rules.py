"""MS/MS の合否（msms-rules.v1）・診断イオン・中性損失（spec 2026-10-22 §4 の 6）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。照合の結果
（`library.matching.score_candidates` の `scored`）と測定スペクトルを受け、spec に書いた
規則を機械的に当てる。**規則に無い判断はしない**（「明らかに別物」かどうかは再解析の
計画の規則で決める。道具は材料と、書かれた規則への合否だけを返す）。

設定していない検査は行を出さない。`verdict` は、すべて pass → `qualified`、
どれか fail → `not_qualified`、fail が無く判定できない検査がある → `not_evaluable`。
"""
from __future__ import annotations

import copy
import math
import re

from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["RULES_SCHEMA", "find_ions", "normalize_ions", "normalize_rules", "qualify",
           "rules_hash"]

RULES_SCHEMA = "msms-rules.v1"
_INVALID = "MSMS_RULES_INVALID"
_BLOCK_RE = re.compile(r"^[A-Z]{14}$")
DEFAULT_ION_TOL = 0.01

_DEFAULTS = {
    "top_rank": 5,
    "min_fragments_excl_precursor": None,
    "precursor_exclusion_tol": None,
    "same_formula_margin": None,
    "lipid_class": None,
    "diagnostic_ions": [],
    "neutral_losses": [],
    "ions_match": "all",
}


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


def _positive(value, label: str, *, allow_zero=False) -> float:
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            or value < 0 or (value == 0 and not allow_zero)):
        raise _fail(f"{label} は正の数で指定する: {value!r}", field=label)
    return float(value)


def inchikey_block(value) -> str:
    block = value.strip()[:14].upper() if isinstance(value, str) else ""
    if not _BLOCK_RE.match(block):
        raise _fail("target.inchikey は InChIKey（先頭 14 文字の英大文字ブロックを使う）",
                    field="target.inchikey")
    return block


def normalize_ions(items, *, kind: str, label: str) -> list[dict]:
    """診断イオン（`mz`）か中性損失（`delta`）の一覧。"""
    key = "mz" if kind == "ion" else "delta"
    if not isinstance(items, list):
        raise _fail(f"{label} はリスト", field=label)
    out = []
    for i, item in enumerate(items):
        allowed = {key, "label", "tol"} | ({"min_relative_intensity"} if kind == "ion" else set())
        if not isinstance(item, dict) or key not in item or set(item) - allowed:
            raise _fail(f"{label}[{i}] は {{{key}, label, tol}} の辞書", field=f"{label}[{i}]")
        out.append({
            key: _positive(item[key], f"{label}[{i}].{key}"),
            "label": str(item.get("label") or f"{key}={item[key]}"),
            "tol": _positive(item.get("tol", DEFAULT_ION_TOL), f"{label}[{i}].tol"),
            **({"min_relative_intensity": _positive(item.get("min_relative_intensity", 0.0),
                                                    f"{label}[{i}].min_relative_intensity",
                                                    allow_zero=True)} if kind == "ion" else {}),
        })
    return out


def normalize_rules(rules, *, target: dict | None = None) -> dict:
    """規則を検証して既定を埋める。`target` を渡すと規則の target はそれで決まる（証拠の束）。"""
    if not isinstance(rules, dict):
        raise _fail("rules は辞書で指定する")
    allowed = {"schema", "target", *_DEFAULTS}
    unknown = sorted(set(rules) - allowed)
    if unknown:
        raise _fail(f"rules の未知のキー: {unknown}", unknown=unknown)
    if rules.get("schema", RULES_SCHEMA) != RULES_SCHEMA:
        raise _fail(f"schema は {RULES_SCHEMA}")
    raw_target = target if target is not None else rules.get("target")
    if not isinstance(raw_target, dict) or set(raw_target) - {"inchikey", "formula"}:
        raise _fail("target は {inchikey, formula} の辞書", field="target")
    formula = raw_target.get("formula")
    if formula is not None and (not isinstance(formula, str) or not formula.strip()):
        raise _fail("target.formula は文字列", field="target.formula")
    out = copy.deepcopy(_DEFAULTS)
    out.update({k: copy.deepcopy(v) for k, v in rules.items() if k in _DEFAULTS})
    top_rank = out["top_rank"]
    if isinstance(top_rank, bool) or not isinstance(top_rank, int) or top_rank < 1:
        raise _fail("top_rank は 1 以上の整数", field="top_rank")
    if out["min_fragments_excl_precursor"] is not None:
        value = out["min_fragments_excl_precursor"]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise _fail("min_fragments_excl_precursor は 0 以上の整数か null",
                        field="min_fragments_excl_precursor")
    if out["precursor_exclusion_tol"] is not None:
        out["precursor_exclusion_tol"] = _positive(out["precursor_exclusion_tol"],
                                                   "precursor_exclusion_tol")
    margin = out["same_formula_margin"]
    if margin is not None:
        if not isinstance(margin, dict) or set(margin) - {"score", "min_delta"}:
            raise _fail("same_formula_margin は {score, min_delta}", field="same_formula_margin")
        score = margin.get("score", "total_score")
        if score not in ("total_score", "weighted_dot_product", "reverse_dot_product",
                         "simple_dot_product", "entropy_similarity"):
            raise _fail(f"same_formula_margin.score を読めない: {score!r}",
                        field="same_formula_margin.score")
        out["same_formula_margin"] = {"score": score,
                                      "min_delta": float(margin.get("min_delta", 0.0))}
    if out["lipid_class"] is not None and (not isinstance(out["lipid_class"], str)
                                           or not out["lipid_class"].strip()):
        raise _fail("lipid_class は文字列か null", field="lipid_class")
    if out["ions_match"] not in ("all", "any"):
        raise _fail("ions_match は all / any", field="ions_match")
    out["diagnostic_ions"] = normalize_ions(out["diagnostic_ions"], kind="ion",
                                            label="diagnostic_ions")
    out["neutral_losses"] = normalize_ions(out["neutral_losses"], kind="loss",
                                           label="neutral_losses")
    return {"schema": RULES_SCHEMA,
            "target": {"inchikey_block": inchikey_block(raw_target.get("inchikey")),
                       "formula": formula.strip() if formula else None},
            **out}


def rules_hash(normalized: dict) -> str:
    return canonical_hash(normalized)


# ---------- 照会 ----------

def find_ions(measured, ions: list[dict], *, precursor_mz: float | None) -> list[dict]:
    """測定スペクトルで診断イオン（`mz`）・中性損失（`delta`）を探す。

    相対強度は測定スペクトルの最大ピークに対する比。±tol に複数あれば最も強いもの。
    """
    peaks = [(float(mz), float(intensity)) for mz, intensity in (measured or [])]
    top = max((i for _, i in peaks), default=0.0)
    out = []
    for ion in ions:
        if "mz" in ion:
            kind, target = "ion", ion["mz"]
        else:
            kind = "neutral_loss"
            target = None if precursor_mz is None else precursor_mz - ion["delta"]
        hits = [] if target is None else [(mz, i) for mz, i in peaks
                                          if abs(mz - target) <= ion["tol"] + 1e-9]
        best = max(hits, key=lambda p: p[1], default=None)
        relative = best[1] / top if best and top > 0 else None
        minimum = ion.get("min_relative_intensity", 0.0)
        out.append({
            "label": ion["label"], "kind": kind, "target_mz": target,
            "found": bool(best) and (relative or 0.0) >= minimum,
            "observed_mz": best[0] if best else None,
            "relative_intensity": relative,
        })
    return out


def _block(record) -> str | None:
    key = record.get("inchikey")
    return key.strip()[:14].upper() if isinstance(key, str) and key.strip() else None


def _lipid_class(name) -> str | None:
    from metabolomix.msdial.lipid_identity import normalize_lipid_name

    if not isinstance(name, str) or not name.strip():
        return None
    parsed = normalize_lipid_name(name)
    if not parsed.get("parse_ok") or not parsed.get("normalized"):
        return None
    return parsed["normalized"].split(" ", 1)[0]


def _check(check: str, status: str, value=None, detail=None) -> dict:
    return {"check": check, "status": status, "value": value, "detail": detail}


def qualify(scored: list, measured, precursor_mz: float, rules: dict, *, ms2_tol: float) -> dict:
    """照合の結果（`score_candidates` の `scored`。総合スコアの降順）に規則を当てる。"""
    block = rules["target"]["inchikey_block"]
    checks = []
    target_index = next((i for i, (record, _) in enumerate(scored) if _block(record) == block), None)
    in_top = target_index is not None and target_index < rules["top_rank"]
    checks.append(_check("target_in_top_rank", "pass" if in_top else "fail",
                         None if target_index is None else target_index + 1,
                         {"top_rank": rules["top_rank"]}))
    target = scored[target_index] if target_index is not None else None

    if rules["min_fragments_excl_precursor"] is not None:
        tol = rules["precursor_exclusion_tol"] or ms2_tol
        if target is None:
            checks.append(_check("min_fragments", "not_evaluable", None,
                                 {"reason": "target_not_among_candidates"}))
        else:
            count = fragment_matches_excluding_precursor(target[1].get("alignment") or [],
                                                         precursor_mz, tol)
            status = ("not_evaluable" if count is None
                      else "pass" if count >= rules["min_fragments_excl_precursor"] else "fail")
            checks.append(_check("min_fragments", status, count,
                                 {"minimum": rules["min_fragments_excl_precursor"],
                                  "precursor_exclusion_tol": tol}))

    margin = rules["same_formula_margin"]
    if margin is not None:
        if target is None:
            checks.append(_check("same_formula_margin", "not_evaluable", None,
                                 {"reason": "target_not_among_candidates"}))
        else:
            formula = rules["target"]["formula"] or target[0].get("formula")
            alternatives = [(r, s) for r, s in scored
                            if formula and (r.get("formula") or "").strip() == formula
                            and _block(r) not in (None, block)]
            if not formula:
                checks.append(_check("same_formula_margin", "not_evaluable", None,
                                     {"reason": "formula_unknown"}))
            elif not alternatives:
                checks.append(_check("same_formula_margin", "pass", None,
                                     {"reason": "no_same_formula_alternative"}))
            else:
                best = max(alternatives, key=lambda item: item[1][margin["score"]])
                delta = target[1][margin["score"]] - best[1][margin["score"]]
                checks.append(_check("same_formula_margin",
                                     "pass" if delta >= margin["min_delta"] else "fail", delta,
                                     {"score": margin["score"], "min_delta": margin["min_delta"],
                                      "best_alternative": best[0].get("name"),
                                      "best_alternative_inchikey": best[0].get("inchikey")}))

    if rules["lipid_class"] is not None:
        top_name = scored[0][0].get("name") if scored else None
        found = _lipid_class(top_name)
        status = ("not_evaluable" if found is None
                  else "pass" if found == rules["lipid_class"] else "fail")
        checks.append(_check("lipid_class", status, found,
                             {"expected": rules["lipid_class"], "rank1_name": top_name}))

    for check, ions in (("diagnostic_ions", rules["diagnostic_ions"]),
                        ("neutral_losses", rules["neutral_losses"])):
        if not ions:
            continue
        found = find_ions(measured, ions, precursor_mz=precursor_mz)
        hits = [f["found"] for f in found]
        passed = all(hits) if rules["ions_match"] == "all" else any(hits)
        checks.append(_check(check, "pass" if passed else "fail", sum(hits),
                             {"ions": found, "match": rules["ions_match"]}))

    statuses = {c["status"] for c in checks}
    verdict = ("not_qualified" if "fail" in statuses
               else "not_evaluable" if "not_evaluable" in statuses else "qualified")
    return {"verdict": verdict, "checks": checks,
            "target_rank": None if target_index is None else target_index + 1}
