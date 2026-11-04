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

from metabolomix.analysis import feature_search, mass_accuracy, msms_rules
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
    allowed = {"schema", "claim_id", "target", "scan", "msms", "relations", "max_candidates",
               "merge_split", "ppm_offset"}
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
    # 質量の窓の中心のずれ（任意）。既定の目標と走査の全変種に同じ spec を渡す（幅は変えない）
    offset = mass_accuracy.normalize_ppm_offset(spec.get("ppm_offset"), code="CLAIM_SPEC_INVALID")
    target = feature_search.normalize_targets([raw_target], dataset_polarity=dataset_polarity,
                                              ppm_offset=offset)[0]

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
            dataset_polarity=dataset_polarity, ppm_offset=offset)

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
    merge = spec.get("merge_split", True)
    if not isinstance(merge, bool):
        raise _fail("merge_split は true / false", field="merge_split")
    out = {"schema": SPEC_SCHEMA, "claim_id": claim_id.strip(), "target": target,
           "raw_target": raw_target, "scan": out_scan, "msms": out_msms,
           "relations": out_relations, "max_candidates": max_candidates}
    if merge:
        # 既定（true）だけを規則に載せる。false は v7 と同じ規則（同じ rule_hash）になる
        out["merge_split"] = True
    if offset is not None:
        # 渡したときだけ規則に載せる（省けば v8 と同じ spec・同じ rule_hash）
        out["ppm_offset"] = offset
    return out


def rule_hash(normalized: dict, merge_groups_hash: str | None = None) -> str:
    """規則の hash。`merge_split` が真なら、判定でまとめた割れの組の hash（`merge_groups_hash`）も
    入れる（同じ spec でも組が変われば別の規則。`analysis_spec.rule_hash` と同じ扱い）。偽なら v7 と同じ。"""
    rule = {k: v for k, v in normalized.items() if k != "raw_target"}
    if "lipid_species_input" in rule["target"]:
        # 種の綴り（PC 18:3_16:1 と PC 34:4）は同じ規則（feature_search.rule_hash と同じ扱い）
        rule["target"] = {k: v for k, v in rule["target"].items() if k != "lipid_species_input"}
    if normalized.get("merge_split") and merge_groups_hash is not None:
        rule["merge_groups_hash"] = merge_groups_hash
    return canonical_hash(rule)


def _selection(match: dict, verdicts: dict[str, str | None] | None) -> dict:
    """目標照合の判定と段ごとの件数（`counts`）。割れの組は 1 つの候補（`match_target` の
    `merged_from`。`selected_feature_id` は残す ID、`merged_from` にその群）。

    `counts.qualified` は**規則で qualified になった特徴（質量・RT・注釈）のうち、MS/MS の合否が
    `qualified` の数**（`verdicts` は特徴 ID → `msms.qualify.verdict`。MS/MS を頼んでいない・
    ライブラリが無いときは None で、数も None。割れの組はどちらかのスポット——窓の外の
    `merge_partner` も含む——が `qualified` なら 1）。
    `counts.msms_available` は**旧 `qualified` の意味の数**——規則（質量の窓・RT の窓・注釈）で
    qualified になった特徴の数で、MS/MS があるかどうかの数ではない（名前は計画のとおり。
    `dataset_feature_search` の `counts.qualified` と同じ数）。
    """
    out = {k: match[k] for k in ("status", "reason", "selected_feature_id", "selected_adduct",
                                 "selected_evidence", "adduct_tier_used")}
    if "merged_from" in match:
        out["merged_from"] = match["merged_from"]
    counts = match["counts"]
    units: dict[str, set[str]] = {}
    for row in match["candidates"]:
        if row["qualified"]:
            # 割れの組は群の両スポット（窓の外の merge_partner も含む）を見る
            spots = row.get("merged_from") or [row["feature_id"]]
            units.setdefault(spots[0], set()).update(spots)
    out["counts"] = {
        "mass": counts["mass"], "rt": counts["rt"],
        "qualified": (None if verdicts is None
                      else sum(any(verdicts.get(fid) == "qualified" for fid in spots)
                               for spots in units.values())),
        "msms_available": counts["qualified"],
        "by_evidence": dict(counts["by_evidence"]),
    }
    return out


def _split(info: dict | None) -> dict | None:
    """候補の割れの相手（`alignment_qc.split_partners` の 1 件）。組に入らなければ None。

    `partners` はまとめる相手（`merge_groups`。互いに最も近い mergeable な組）で、まとめない
    特徴（`group_aligned` など）では空。`component` はすべての組の連結成分のほかの特徴（報告用）。
    """
    if info is None:
        return None
    return {"group_id": info["group_id"], "partners": list(info["partners"]),
            "type": info["type"], "mergeable": info["mergeable"],
            "component": list(info["component"])}


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
    """証拠の束（claim-evidence.v1）を作る。`spot_query` は `library.tools.spot_msms_query`。

    候補ごとの `split` は割れたアライメントスポットの相手（`alignment_qc.find_split_pairs` を
    既定の閾値で回した結果。`dataset_alignment_qc` の既定と同じ）。`merge_split`（spec。既定 true）
    なら、その `merge_groups` の組を既定の特徴と走査の判定で 1 つの候補として数え（候補の一覧は
    スポットごとのまま）、組の片方だけが窓に入ったときは相手も候補に足す（`scans` に
    `merge_partner`。`split`・`msms`・`relations` も相手について作る）。束の `merge_split` と
    `merge_groups_hash`（判定でまとめた組の hash）を規則の隣に置き、真のときだけ `rule_hash` に入れる。
    """
    from metabolomix.analysis.analysis_spec import merge_groups_hash
    from metabolomix.analysis import alignment_qc
    from metabolomix.analysis.result_state import dataset_fingerprint

    dataset_hash = dataset_fingerprint(ds)
    qc = alignment_qc.find_split_pairs(ds, dataset_hash=dataset_hash)
    partners = alignment_qc.split_partners(qc)
    # 判定（既定の特徴・走査）でも割れの組は 1 つの候補（dataset_feature_search の merge_split と同じ）
    merge = bool(normalized.get("merge_split"))
    merged_map = (feature_search.merge_map(alignment_qc.usable_merge_groups(qc, ds.feature_ids)[0])
                  if merge else None)
    groups_hash = merge_groups_hash(list(merged_map.values())) if merge else None
    mz_array = feature_search._feature_mz_array(ds)
    unlimited = len(ds.feature_ids) * 16 + 1
    default = feature_search.match_target(ds, normalized["target"], max_results=unlimited,
                                          mz_array=mz_array, merged_from=merged_map)
    seen: dict[str, dict] = {}

    def _collect(match, scan_label):
        for row in match["candidates"]:
            entry = seen.setdefault(row["feature_id"], {
                "feature_id": row["feature_id"], "mz": row["mz"], "rt": row["rt"],
                "evidence": set(), "adducts": set(), "annotated_name": row["annotated_name"],
                "annotated_inchikey": row["annotated_inchikey"], "scans": set()})
            entry["evidence"].add(row["evidence"])
            entry["adducts"].add(row["adduct"])
            # 割れの相手として足した行（質量の窓の外）は、その走査で見えたものと区別する
            partner = "merge_partner" in row["reasons"].split(";")
            entry["scans"].add("merge_partner" if partner else scan_label)

    _collect(default, "default")
    scan_matches = []
    raw = normalized["raw_target"]
    for ppm in normalized["scan"]["ppm"]:
        for rt_tol in normalized["scan"]["rt_tol_min"]:
            variant = {**raw, "ppm": ppm, "rt_tol_min": rt_tol}
            if normalized["scan"]["adduct_tiers"] is not None:
                variant["adduct_tiers"] = normalized["scan"]["adduct_tiers"]
            target = feature_search.normalize_targets([variant], dataset_polarity=dataset_polarity,
                                                      ppm_offset=normalized.get("ppm_offset"))[0]
            match = feature_search.match_target(ds, target, max_results=unlimited,
                                                mz_array=mz_array, merged_from=merged_map)
            label = f"ppm={ppm:g},rt_tol_min={rt_tol:g}"
            _collect(match, label)
            scan_matches.append((ppm, rt_tol, match))

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
                "scans": sorted(entry["scans"]),
                "split": _split(partners.get(fid))}
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

    # MS/MS の合否（counts.qualified）。頼んでいない・ライブラリが無ければ数えない（None）。
    verdicts = None
    if normalized["msms"] is not None and store is not None:
        verdicts = {c["feature_id"]: ((c["msms"].get("qualify") or {}).get("verdict"))
                    for c in candidates}
    scans = [{"ppm": ppm, "rt_tol_min": rt_tol, **_selection(match, verdicts),
              "candidate_feature_ids": sorted({r["feature_id"] for r in match["candidates"]},
                                              key=order)}
             for ppm, rt_tol, match in scan_matches]
    rule = {k: copy.deepcopy(v) for k, v in normalized.items() if k != "raw_target"}
    return {
        "schema": CLAIM_SCHEMA,
        "claim_id": normalized["claim_id"],
        "rule_hash": rule_hash(normalized, groups_hash),
        "input_fingerprint": {"dataset_hash": dataset_hash,
                              "library_sha256": library_sha256},
        "rule": rule,
        "merge_split": merge,
        "merge_groups_hash": groups_hash,
        "dataset_polarity": polarity,
        "default_selection": _selection(default, verdicts),
        "scans": scans,
        "candidates": candidates,
        "truncated": truncated,
    }


def _feature_order(ds):
    position = {fid: i for i, fid in enumerate(ds.feature_ids)}
    return lambda fid: (position.get(fid, len(position)), fid)
