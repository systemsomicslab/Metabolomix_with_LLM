"""同族列の RT の単調性の検査（spec 2026-10-22 §4 の 6）。

逆相の LC では、同じクラス・同じ二重結合数なら炭素数が増えるほど RT が遅れ、同じ炭素数なら
二重結合が増えるほど RT が早まる。同定された種がこの並びに反していれば、どこかの同定が
疑わしい。数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

並びの向きと許す外れの数は spec の `rule` で決める（既定は逆相の向き・外れ 0）。
同じ（炭素数, 二重結合数）の特徴が複数あれば RT の中央値で並べ、`duplicates` に挙げる。
隣どうし（同じ二重結合数の中で炭素数の昇順、同じ炭素数の中で二重結合数の昇順）で、
規則の向きに反する対を違反として数える。RT が同じ対は違反にしない。
"""
from __future__ import annotations

import re
from statistics import median

from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["SPEC_SCHEMA", "check", "members_from_names", "normalize_spec", "species_of"]

SPEC_SCHEMA = "homolog-check.v1"
_INVALID = "HOMOLOG_SPEC_INVALID"
_SPECIES_RE = re.compile(r"^(?P<cls>\S+) (?P<c>\d+):(?P<db>\d+)(?:;.*)?$")
_DIRECTIONS = ("increasing", "decreasing")


def _fail(message, **details):
    return DomainError(_INVALID, message, details)


def normalize_spec(spec) -> dict:
    if not isinstance(spec, dict) or set(spec) - {"schema", "members", "auto", "rule"}:
        raise _fail("spec は {members | auto, rule} の辞書")
    if ("members" in spec and spec["members"] is not None) == (
            "auto" in spec and spec["auto"] is not None):
        raise _fail("members と auto のどちらか一方を指定する")
    rule = spec.get("rule") or {}
    if not isinstance(rule, dict) or set(rule) - {"carbon_direction", "double_bond_direction",
                                                   "max_violations"}:
        raise _fail("rule は {carbon_direction, double_bond_direction, max_violations}", field="rule")
    out_rule = {"carbon_direction": rule.get("carbon_direction", "increasing"),
                "double_bond_direction": rule.get("double_bond_direction", "decreasing"),
                "max_violations": rule.get("max_violations", 0)}
    for key in ("carbon_direction", "double_bond_direction"):
        if out_rule[key] not in _DIRECTIONS:
            raise _fail(f"rule.{key} は increasing / decreasing", field=f"rule.{key}")
    if (isinstance(out_rule["max_violations"], bool) or not isinstance(out_rule["max_violations"], int)
            or out_rule["max_violations"] < 0):
        raise _fail("rule.max_violations は 0 以上の整数", field="rule.max_violations")
    out = {"schema": SPEC_SCHEMA, "rule": out_rule, "members": None, "auto": None}
    if spec.get("members") is not None:
        members = spec["members"]
        if not isinstance(members, list) or len(members) < 2:
            raise _fail("members は 2 件以上", field="members")
        normalized = []
        for i, m in enumerate(members):
            if (not isinstance(m, dict) or set(m) != {"feature_id", "carbons", "double_bonds"}
                    or not isinstance(m["feature_id"], str)
                    or not all(isinstance(m[k], int) and not isinstance(m[k], bool) and m[k] >= 0
                               for k in ("carbons", "double_bonds"))):
                raise _fail(f"members[{i}] は {{feature_id, carbons, double_bonds}}", field=f"members[{i}]")
            normalized.append(dict(m))
        out["members"] = sorted(normalized, key=lambda m: (m["carbons"], m["double_bonds"],
                                                           m["feature_id"]))
    else:
        auto = spec["auto"]
        if not isinstance(auto, dict) or set(auto) - {"lipid_class", "evidence", "adduct"} \
                or not isinstance(auto.get("lipid_class"), str):
            raise _fail("auto は {lipid_class, evidence, adduct}", field="auto")
        evidence = auto.get("evidence", ["sme"])
        if not isinstance(evidence, list) or not evidence or set(evidence) - {"sme", "sml"}:
            raise _fail("auto.evidence は sme / sml のリスト", field="auto.evidence")
        out["auto"] = {"lipid_class": auto["lipid_class"],
                       "evidence": [e for e in ("sme", "sml") if e in evidence],
                       "adduct": auto.get("adduct")}
    return out


def species_of(name) -> tuple[str, int, int] | None:
    """名前を GOSLIN で種のレベル（`CLASS C:DB`）に読み、(class, C, DB) を返す。"""
    from metabolomix.msdial.lipid_identity import normalize_lipid_name

    if not isinstance(name, str) or not name.strip():
        return None
    parsed = normalize_lipid_name(name)
    text = parsed.get("normalized") if parsed.get("parse_ok") else None
    if not text:
        return None
    species = text.split("|", 1)[0]
    # 分子種（16:0_18:1）は鎖を足して種のレベルに直す
    head, _, tail = species.partition(" ")
    # エーテル（O- / P-）は別のクラスとして扱う（PC O-34:1 は PC 34:1 の同族ではない）
    ether = re.match(r"([OP])-", tail)
    if ether:
        head = f"{head} {ether.group(1)}-"
    chains = re.findall(r"(\d+):(\d+)", tail)
    if not chains:
        return None
    carbons = sum(int(c) for c, _ in chains)
    double_bonds = sum(int(d) for _, d in chains)
    return head, carbons, double_bonds


def members_from_names(features: list[dict], auto: dict) -> list[dict]:
    """`features` の各要素 {feature_id, names: {sme, sml}, adduct} から同族の候補を集める。"""
    out = []
    for feature in features:
        if auto["adduct"] and feature.get("adduct") != auto["adduct"]:
            continue
        for source in auto["evidence"]:
            species = species_of((feature.get("names") or {}).get(source))
            if species and species[0] == auto["lipid_class"]:
                out.append({"feature_id": feature["feature_id"], "carbons": species[1],
                            "double_bonds": species[2], "name_source": source})
                break
    return out


def _violates(direction: str, rt_low: float, rt_high: float) -> bool:
    """値が小さい側の RT と大きい側の RT が、規則の向きに反しているか。"""
    if rt_low == rt_high:
        return False
    return rt_high < rt_low if direction == "increasing" else rt_high > rt_low


def check(members: list[dict], rt_of: dict, rule: dict) -> dict:
    """members（{feature_id, carbons, double_bonds}）の RT（rt_of: feature_id → 分）を検査する。"""
    groups: dict[tuple[int, int], list[str]] = {}
    for m in members:
        if rt_of.get(m["feature_id"]) is not None:
            groups.setdefault((m["carbons"], m["double_bonds"]), []).append(m["feature_id"])
    rt_by_species = {key: median(rt_of[f] for f in fids) for key, fids in groups.items()}
    duplicates = {f"{c}:{d}": sorted(fids) for (c, d), fids in groups.items() if len(fids) > 1}
    violations = []
    for axis, direction in (("carbons", rule["carbon_direction"]),
                            ("double_bonds", rule["double_bond_direction"])):
        other = "double_bonds" if axis == "carbons" else "carbons"
        series: dict[int, list[tuple[int, int]]] = {}
        for key in rt_by_species:
            species = {"carbons": key[0], "double_bonds": key[1]}
            series.setdefault(species[other], []).append(key)
        for keys in series.values():
            keys.sort(key=lambda k: k[0] if axis == "carbons" else k[1])
            for low, high in zip(keys, keys[1:]):
                if _violates(direction, rt_by_species[low], rt_by_species[high]):
                    violations.append({"axis": axis, "low": f"{low[0]}:{low[1]}",
                                       "high": f"{high[0]}:{high[1]}",
                                       "rt_low": rt_by_species[low], "rt_high": rt_by_species[high],
                                       "features_low": sorted(groups[low]),
                                       "features_high": sorted(groups[high])})
    violations.sort(key=lambda v: (v["axis"], v["low"], v["high"]))
    flagged = {f for v in violations for f in v["features_low"] + v["features_high"]}
    rows = [{"feature_id": m["feature_id"], "carbons": m["carbons"],
             "double_bonds": m["double_bonds"], "rt": rt_of.get(m["feature_id"]),
             "status": ("no_rt" if rt_of.get(m["feature_id"]) is None
                        else "violation" if m["feature_id"] in flagged else "ok")}
            for m in sorted(members, key=lambda m: (m["carbons"], m["double_bonds"],
                                                    m["feature_id"]))]
    return {"verdict": "pass" if len(violations) <= rule["max_violations"] else "fail",
            "n_violations": len(violations), "violations": violations, "members": rows,
            "duplicates": duplicates}


def rule_hash(normalized: dict) -> str:
    return canonical_hash(normalized)
