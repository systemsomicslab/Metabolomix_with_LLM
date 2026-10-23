"""特徴どうしの関係と、1 主張の証拠の束（spec 2026-10-22 §4 の 7・8）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない（DatasetState・
ライブラリの store・`.arf2` / `.dcl` のパスを引数で受ける）。

**証拠の束は決めない。** 候補の網羅の走査（付加体・許容幅を広げた走査を含む）・各候補の
照合と MS/MS の合否・関係・**規則で選んだ既定の特徴**（`dataset_feature_search` の目標照合と
同じ規則。決まらなければ `needs_input`）を 1 つにまとめるだけ。特徴を最終的に決めるのは
再解析の実行役で、その特徴での統計は analysis-spec で別に呼ぶ。

**同じ spec・同じ入力なら同じ JSON。** 時刻を入れず、並びは特徴 ID・規則の順で固定する。
"""
from __future__ import annotations

import copy
import math

import numpy as np

from metabolomix.analysis import feature_search, msms_rules
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["CLAIM_SCHEMA", "build", "feature_relations", "normalize_claim_spec"]

CLAIM_SCHEMA = "claim-evidence.v1"
SPEC_SCHEMA = "claim-evidence-spec.v1"
_INVALID = "CLAIM_SPEC_INVALID"
_DEFAULT_RELATIONS = {"rt_window": 0.05, "mz_tol": 0.01, "min_r": 0.7}


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


# ---------- 特徴どうしの関係 ----------

def _heights(ds, index: int) -> list:
    row = np.asarray(ds.feature_matrix[index], dtype=float)
    return [float(v) if math.isfinite(v) else None for v in row]


def _top_candidate(ds, fid) -> dict:
    candidates = (ds.feature_candidates or {}).get(fid) or []
    return candidates[0] if candidates else {}


def formula_index(ds) -> dict:
    """SME_ID → 分子式（SME 行の `chemical_formula`。DatasetState の候補は持たないので行から引く）。"""
    out = {}
    for row in getattr(ds, "sme_rows", None) or []:
        formula = row.get("chemical_formula")
        if formula not in (None, "", "null"):
            out[str(row.get("SME_ID"))] = formula
    return out


def _partner(ds, index: int, formulas: dict) -> dict:
    fid = ds.feature_ids[index]
    meta = ds.feature_metadata.get(fid) or {}
    top = _top_candidate(ds, fid)
    return {"spot_id": fid, "mz": meta.get("mz"), "rt": meta.get("rt"),
            "name": top.get("chemical_name"), "adduct": top.get("adduct"),
            "formula": formulas.get(top.get("sme_id")), "heights": _heights(ds, index)}


def feature_relations(ds, feature_id: str, *, rt_window: float, mz_tol: float, min_r: float,
                      adducts: list[str], formulas: dict | None = None) -> dict:
    """同時溶出の特徴と、付加体・同位体・源内断片の関係（`curation.relations` を流用）。

    `as_ion_of`: この特徴が、相手の（付加体・同位体・源内断片の）イオンとして説明できる関係。
    `ions_of_this`: 相手が、この特徴のイオンとして説明できる関係。
    MS-DIAL のリンク（`.arf2`）は使わないので、裏付けは試料間の強度の相関だけ。
    """
    from metabolomix.curation.relations import find_relations, profile_correlation

    formulas = formula_index(ds) if formulas is None else formulas
    index_of = {fid: i for i, fid in enumerate(ds.feature_ids)}
    if feature_id not in index_of:
        raise DomainError("FEATURE_NOT_FOUND", f"特徴が無い: {feature_id!r}",
                          {"feature_ids": [feature_id]})
    target = _partner(ds, index_of[feature_id], formulas)
    if target["rt"] is None or target["mz"] is None:
        raise DomainError("FEATURE_NOT_FOUND", f"特徴 {feature_id} に m/z か RT が無い",
                          {"feature_ids": [feature_id]})
    partners = [_partner(ds, i, formulas) for i, fid in enumerate(ds.feature_ids)
                if fid != feature_id and (ds.feature_metadata.get(fid) or {}).get("rt") is not None
                and abs(ds.feature_metadata[fid]["rt"] - target["rt"]) <= rt_window + 1e-9]
    partners = [p for p in partners if p["mz"] is not None]
    coeluting = sorted(
        ({"feature_id": p["spot_id"], "mz": p["mz"], "rt": p["rt"],
          "drt": p["rt"] - target["rt"],
          "r": profile_correlation(target["heights"], p["heights"]), "name": p["name"]}
         for p in partners),
        key=lambda row: (abs(row["drt"]), row["feature_id"]))
    as_ion_of = find_relations(target, partners, adducts=adducts, rt_window=rt_window,
                               mz_tol=mz_tol, min_r=min_r, links={})
    ions_of_this = []
    for partner in partners:
        for entry in find_relations(partner, [target], adducts=adducts, rt_window=rt_window,
                                    mz_tol=mz_tol, min_r=min_r, links={}):
            ions_of_this.append({**entry, "ion": partner["spot_id"], "of": feature_id})
    ions_of_this.sort(key=lambda e: (str(e["ion"]), e["relation"]))
    for entry in as_ion_of:
        entry.pop("candidate_id", None)
    for entry in ions_of_this:
        entry.pop("candidate_id", None)
    return {"feature_id": feature_id, "coeluting": coeluting, "as_ion_of": as_ion_of,
            "ions_of_this": ions_of_this}


# ---------- 証拠の束 ----------

def normalize_claim_spec(spec, *, dataset_polarity: str | None) -> dict:
    allowed = {"schema", "claim_id", "target", "scan", "msms", "relations", "max_candidates"}
    if not isinstance(spec, dict) or set(spec) - allowed:
        raise _fail(f"spec の未知のキー: {sorted(set(spec or {}) - allowed)}")
    if spec.get("schema", SPEC_SCHEMA) != SPEC_SCHEMA:
        raise _fail(f"schema は {SPEC_SCHEMA}")
    claim_id = spec.get("claim_id")
    if not isinstance(claim_id, str) or not claim_id.strip():
        raise _fail("claim_id は空でない文字列", field="claim_id")
    raw_target = spec.get("target")
    if not isinstance(raw_target, dict):
        raise _fail("target は dataset_feature_search の目標と同じ形の辞書", field="target")
    raw_target = {"target_id": claim_id.strip(), **raw_target}
    target = feature_search.normalize_targets([raw_target], dataset_polarity=dataset_polarity)[0]

    scan = spec.get("scan") or {}
    if not isinstance(scan, dict) or set(scan) - {"ppm", "rt_tol_min", "adduct_tiers"}:
        raise _fail("scan は {ppm, rt_tol_min, adduct_tiers}", field="scan")

    def _levels(key, default):
        values = scan.get(key, default)
        if (not isinstance(values, list) or not values
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0
                           for v in values)):
            raise _fail(f"scan.{key} は正の数のリスト", field=f"scan.{key}")
        return sorted({float(v) for v in values})

    out_scan = {"ppm": _levels("ppm", [target["ppm"]]),
                "rt_tol_min": _levels("rt_tol_min", [target["rt_tol_min"]]),
                "adduct_tiers": scan.get("adduct_tiers")}
    if out_scan["adduct_tiers"] is not None:
        if target["expected_mz"] is not None:
            raise _fail("expected_mz の目標は付加体が固定なので scan.adduct_tiers は使えない",
                        field="scan.adduct_tiers")
        # 検証は目標の正規化に任せる（読めない付加体・極性の食い違いを同じ規則で落とす）
        feature_search.normalize_targets(
            [{**raw_target, "adduct_tiers": out_scan["adduct_tiers"]}],
            dataset_polarity=dataset_polarity)

    msms = spec.get("msms")
    out_msms = None
    if msms is not None:
        if not isinstance(msms, dict) or set(msms) - {"top_n", "rules"}:
            raise _fail("msms は {top_n, rules}", field="msms")
        top_n = msms.get("top_n", 10)
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
            raise _fail("msms.top_n は 1 以上の整数", field="msms.top_n")
        rules = msms.get("rules") or {}
        if "target" in rules:
            raise _fail("msms.rules.target は書かない（主張の target から決まる）",
                        field="msms.rules")
        # 種の主張は種を、InChIKey の主張は InChIKey を規則の target に渡す
        identifier = ({"lipid_species": target["lipid_species"]} if target.get("lipid_species")
                      else {"inchikey": raw_target.get("inchikey")})
        out_msms = {"top_n": top_n, "rules": msms_rules.normalize_rules(
            rules, target={**identifier, "formula": target["formula"]})}

    relations = spec.get("relations")
    out_relations = None
    if relations is not None:
        if not isinstance(relations, dict) or set(relations) - set(_DEFAULT_RELATIONS):
            raise _fail("relations は {rt_window, mz_tol, min_r}", field="relations")
        out_relations = {k: float(relations.get(k, v)) for k, v in _DEFAULT_RELATIONS.items()}

    max_candidates = spec.get("max_candidates", 50)
    if isinstance(max_candidates, bool) or not isinstance(max_candidates, int) or max_candidates < 1:
        raise _fail("max_candidates は 1 以上の整数", field="max_candidates")
    return {"schema": SPEC_SCHEMA, "claim_id": claim_id.strip(), "target": target,
            "raw_target": raw_target, "scan": out_scan, "msms": out_msms,
            "relations": out_relations, "max_candidates": max_candidates}


def rule_hash(normalized: dict) -> str:
    rule = {k: v for k, v in normalized.items() if k != "raw_target"}
    if "lipid_species_input" in rule["target"]:
        # 種の綴り（PC 18:3_16:1 と PC 34:4）は同じ規則（feature_search.rule_hash と同じ扱い）
        rule["target"] = {k: v for k, v in rule["target"].items() if k != "lipid_species_input"}
    return canonical_hash(rule)


def _selection(match: dict) -> dict:
    return {k: match[k] for k in ("status", "reason", "selected_feature_id", "selected_adduct",
                                  "selected_evidence", "adduct_tier_used", "counts")}


def _msms_for(fid, *, store, spot_query, arf2_path, dcl_path, msms_spec) -> dict:
    from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor
    from metabolomix.library import matching

    if store is None:
        return {"status": "no_library"}
    try:
        spot_id = int(fid)
    except (TypeError, ValueError):
        return {"status": "spot_id_not_integer"}
    spot = spot_query(spot_id, arf2_path=arf2_path, dcl_file=dcl_path, rt=None, ion_mode=None)
    if spot["status"] != "ok":
        return {"status": spot["status"], "message": spot.get("message")}
    matched = matching.score_candidates(store, precursor_mz=spot["precursor_mz"], rt=spot["rt"],
                                        ion_mode=spot["ion_mode"], measured=spot["measured"])
    ms2_tol = matched["query"]["ms2_tol"]
    tol = msms_spec["rules"]["precursor_exclusion_tol"] or ms2_tol
    top = []
    for rank, (record, result) in enumerate(matched["scored"][: msms_spec["top_n"]], start=1):
        top.append({"rank": rank, "name": record.get("name"), "formula": record.get("formula"),
                    "inchikey": record.get("inchikey"), "adduct": record.get("adduct"),
                    "total_score": result["total_score"],
                    "matched_peaks_count": result["matched_peaks_count"],
                    "fragment_matches_excl_precursor": fragment_matches_excluding_precursor(
                        result["alignment"], spot["precursor_mz"], tol)})
    if matched["status"] == "no_candidates":
        qualify = {"verdict": "not_evaluable", "checks": [], "target_rank": None,
                   "reason": "no_library_candidates"}
    else:
        qualify = msms_rules.qualify(matched["scored"], spot["measured"], spot["precursor_mz"],
                                     msms_spec["rules"], ms2_tol=ms2_tol)
    return {"status": matched["status"], "precursor_mz": spot["precursor_mz"],
            "n_measured_peaks": len(spot["measured"]), "n_candidates": matched["n_candidates"],
            "top": top, "qualify": qualify}


def build(ds, normalized: dict, *, store, spot_query, arf2_path, dcl_path,
          dataset_polarity: str | None, library_sha256: str | None) -> dict:
    """証拠の束（claim-evidence.v1）を作る。`spot_query` は `library.tools.spot_msms_query`。"""
    from metabolomix.analysis.result_state import dataset_fingerprint

    mz_array = feature_search._feature_mz_array(ds)
    unlimited = len(ds.feature_ids) * 16 + 1
    default = feature_search.match_target(ds, normalized["target"], max_results=unlimited,
                                          mz_array=mz_array)
    seen: dict[str, dict] = {}

    def _collect(match, scan_label):
        for row in match["candidates"]:
            entry = seen.setdefault(row["feature_id"], {
                "feature_id": row["feature_id"], "mz": row["mz"], "rt": row["rt"],
                "evidence": set(), "adducts": set(), "annotated_name": row["annotated_name"],
                "annotated_inchikey": row["annotated_inchikey"], "scans": set()})
            entry["evidence"].add(row["evidence"])
            entry["adducts"].add(row["adduct"])
            entry["scans"].add(scan_label)

    _collect(default, "default")
    scans = []
    raw = normalized["raw_target"]
    for ppm in normalized["scan"]["ppm"]:
        for rt_tol in normalized["scan"]["rt_tol_min"]:
            variant = {**raw, "ppm": ppm, "rt_tol_min": rt_tol}
            if normalized["scan"]["adduct_tiers"] is not None:
                variant["adduct_tiers"] = normalized["scan"]["adduct_tiers"]
            target = feature_search.normalize_targets([variant],
                                                      dataset_polarity=dataset_polarity)[0]
            match = feature_search.match_target(ds, target, max_results=unlimited,
                                                mz_array=mz_array)
            label = f"ppm={ppm:g},rt_tol_min={rt_tol:g}"
            _collect(match, label)
            scans.append({"ppm": ppm, "rt_tol_min": rt_tol, **_selection(match),
                          "candidate_feature_ids": sorted({r["feature_id"] for r in match["candidates"]},
                                                          key=_feature_order(ds))})

    order = _feature_order(ds)
    ordered = sorted(seen, key=order)
    truncated = max(0, len(ordered) - normalized["max_candidates"])
    ordered = ordered[: normalized["max_candidates"]]
    polarity = dataset_polarity
    adducts = [a for tier in feature_search.DEFAULT_ADDUCT_TIERS.get(polarity or "", [])
               for a in tier]
    formulas = formula_index(ds) if normalized["relations"] is not None else {}
    candidates = []
    for fid in ordered:
        entry = seen[fid]
        item = {"feature_id": fid, "mz": entry["mz"], "rt": entry["rt"],
                "evidence": sorted(entry["evidence"]), "adducts": sorted(entry["adducts"]),
                "annotated_name": entry["annotated_name"],
                "annotated_inchikey": entry["annotated_inchikey"],
                "scans": sorted(entry["scans"])}
        if normalized["msms"] is not None:
            item["msms"] = _msms_for(fid, store=store, spot_query=spot_query,
                                     arf2_path=arf2_path, dcl_path=dcl_path,
                                     msms_spec=normalized["msms"])
        if normalized["relations"] is not None:
            relations = feature_relations(ds, fid, adducts=adducts, formulas=formulas,
                                          **normalized["relations"])
            item["relations"] = {"as_ion_of": relations["as_ion_of"],
                                 "ions_of_this": relations["ions_of_this"],
                                 "n_coeluting": len(relations["coeluting"])}
        candidates.append(item)

    rule = {k: copy.deepcopy(v) for k, v in normalized.items() if k != "raw_target"}
    return {
        "schema": CLAIM_SCHEMA,
        "claim_id": normalized["claim_id"],
        "rule_hash": rule_hash(normalized),
        "input_fingerprint": {"dataset_hash": dataset_fingerprint(ds),
                              "library_sha256": library_sha256},
        "rule": rule,
        "dataset_polarity": polarity,
        "default_selection": _selection(default),
        "scans": scans,
        "candidates": candidates,
        "truncated": truncated,
    }


def _feature_order(ds):
    position = {fid: i for i, fid in enumerate(ds.feature_ids)}
    return lambda fid: (position.get(fid, len(position)), fid)
