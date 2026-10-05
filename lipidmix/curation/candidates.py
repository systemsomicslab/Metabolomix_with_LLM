# lipidmix/curation/candidates.py
"""候補付けの ① MS-DIAL の下位候補と ② 閾値を緩めた再検索。spec §5。

両方を同じ物差し（`match_spectrum` + `total_score`、採点の前処理は store の search_params）で
採点し、同じレコードは 1 つにまとめる。制約（P）: ハードは削る、ソフトは順位を下げる、情報は
何もしない。並び順は（ソフト理由の数 昇順, total_score 降順）。
deps: analysis.spectral_match、library.store、msdial.adducts、msdial.peak_verification、
plots.mirror、curation.evidence（payload の間引き）、curation.trend。tools_* は import しない。
"""
from __future__ import annotations

from lipidmix.analysis.spectral_match import match_spectrum, total_score
from lipidmix.curation import trend
from lipidmix.curation.evidence import MIRROR_MAX_PEAKS, REFERENCE_MZ_WINDOW, _keep_most_intense_peaks
from lipidmix.msdial.adducts import parse_adduct
from lipidmix.msdial.peak_verification import adduct_consistency
from lipidmix.plots.mirror import build_mirror_payload

SCORE_KEYS = ("total_score", "weighted_dot_product", "simple_dot_product", "reverse_dot_product",
              "matched_peaks_percentage", "matched_peaks_count")


def record_key(record: dict) -> tuple:
    return (record.get("library_id"), int(record["record_index"]))


def _polarity(ion_mode) -> str | None:
    mode = str(ion_mode or "").strip().lower()
    return "+" if mode.startswith("pos") else "-" if mode.startswith("neg") else None


def constraint_reasons(identity: dict, *, rep_mz, ion_mode, measured, scores, trend_entry, th) -> dict:
    hard, soft, info = [], [], []
    adduct = parse_adduct(identity.get("adduct"))
    polarity = _polarity(ion_mode)
    if adduct is not None and polarity is not None and adduct.polarity != polarity:
        hard.append("polarity_mismatch")
    ref_mz = identity.get("precursor_mz")
    if ref_mz is not None and rep_mz is not None and abs(float(rep_mz) - float(ref_mz)) * 1000 >= th["dmz_fail_mda"]:
        hard.append("dmz_out")
    if adduct_consistency(identity.get("adduct"), ion_mode, identity.get("ontology")).get("class_typical") is False:
        soft.append("adduct_atypical")
    if not measured:
        soft.append("msms_absent")
    elif (scores or {}).get("weighted_dot_product", 0.0) >= 0 and (scores or {}).get("matched_peaks_count", 0) <= 0:
        # 参照スペクトルが無く比較していない候補（スコアが -1 の印）には付けない（reference_unresolved の領分）
        soft.append("no_matched_peaks")
    if trend_entry is None:
        info.append("trend_unknown")
    elif abs(trend_entry["z"]) > th["trend_outlier_z"]:
        soft.append("trend_outlier")
    return {"hard": hard, "soft": soft, "info": info}


def rank(items: list[dict]) -> list[dict]:
    return sorted(items, key=lambda c: (len(c["soft"]), -(c["scores"].get("total_score") or 0.0)))


def _trend_entry(identity: dict, rt, trends: dict) -> dict | None:
    comp = trend.composition(identity.get("name"))
    if comp is None or rt is None:
        return None
    predicted = trend.predict_rt(trends.get("classes") or {}, identity.get("ontology") or "", *comp)
    if predicted is None:
        return None
    residual = float(rt) - predicted["predicted_rt"]
    return {"predicted_rt": round(predicted["predicted_rt"], 3), "residual": round(residual, 3),
            "z": round(residual / predicted["scale"], 2)}


def _score(measured, reference_spectrum, *, rep_mz, rt, ref_mz, ref_rt, scoring) -> dict:
    if not measured or not reference_spectrum:
        return {**{k: -1.0 for k in SCORE_KEYS}, "matched_peaks_count": 0, "alignment": []}
    result = match_spectrum(measured, reference_spectrum, ms2_tol=scoring["ms2_tol"],
                            mass_begin=scoring["mass_begin"], mass_end=scoring["mass_end"],
                            relative_amp_cutoff=scoring["relative_amp_cutoff"],
                            absolute_amp_cutoff=scoring["absolute_amp_cutoff"])
    result.update(total_score(result, precursor_mz=rep_mz, reference_precursor_mz=ref_mz,
                              ms1_tol=scoring["mz_tol"], rt=rt, reference_rt=ref_rt,
                              rt_tol=scoring.get("rt_tol"), use_rt=scoring["use_rt"]))
    return result


def _from_record(record: dict, source: str) -> dict:
    return {"source": source, "name": record["name"], "ontology": record.get("ontology") or record.get("compound_class"),
            "adduct": record.get("adduct"), "formula": record.get("formula"), "inchikey": record.get("inchikey"),
            "precursor_mz": record.get("precursor_mz"), "ref_rt": record.get("rt"),
            "library_id": record.get("library_id"), "record_index": record.get("record_index"),
            "_spectrum": record.get("spectrum") or [], "msdial_total_score": None}


def _from_match(match: dict, store) -> dict:
    return {"source": "msdial", "name": match.get("name"), "ontology": None, "adduct": None,
            "formula": None, "inchikey": match.get("inchikey"), "precursor_mz": None, "ref_rt": None,
            "library_id": store.library_id_for(match.get("annotator_id")),
            "record_index": match.get("library_id"), "_spectrum": [],
            "msdial_total_score": match.get("total_score"), "_unresolved": True}


def _resolve_match(store, match: dict, rep_mz) -> dict | None:
    if match.get("library_id") is None or rep_mz is None:
        return None
    record = store.record_by_scan_id(int(match["library_id"]),
                                     library_id=store.library_id_for(match.get("annotator_id")),
                                     precursor_mz=float(rep_mz), mz_tol=REFERENCE_MZ_WINDOW)
    if record is None:
        return None
    if match.get("inchikey") and record.get("inchikey") and match["inchikey"] != record["inchikey"]:
        return None
    return record


def build_library_candidates(ev: dict, *, msdial_matches, representative, store, scoring, trends, th,
                             exclude_current, top_n) -> dict:
    measured = ev.get("_measured") or []
    rep_mz, rt = ev.get("rep_mz"), ev.get("rep_rt")
    pool: dict[tuple, dict] = {}
    unresolved: list[dict] = []
    current = _resolve_match(store, representative, rep_mz) if representative else None
    current_key = record_key(current) if current else None

    for match in msdial_matches:
        record = _resolve_match(store, match, rep_mz)
        if record is None:
            unresolved.append(_from_match(match, store))
            continue
        item = pool.setdefault(record_key(record), _from_record(record, "msdial"))
        item["msdial_total_score"] = match.get("total_score")
    if measured and rep_mz is not None:
        for record in store.candidates(float(rep_mz), mz_tol=scoring["mz_tol"], ion_mode=ev.get("ion_mode")):
            key = record_key(record)
            if key in pool:
                if pool[key]["source"] == "msdial":
                    pool[key]["source"] = "msdial+research"
            else:
                pool[key] = _from_record(record, "research")
    if exclude_current and current_key is not None:
        pool.pop(current_key, None)
    if exclude_current and representative is not None:
        unresolved = [u for u in unresolved if u["name"] != representative.get("name")]

    kept, n_hard, matched_measured = [], 0, set()
    for item in [*pool.values(), *unresolved]:
        result = _score(measured, item["_spectrum"], rep_mz=rep_mz, rt=rt, ref_mz=item["precursor_mz"],
                        ref_rt=item["ref_rt"], scoring=scoring)
        trend_entry = _trend_entry(item, rt, trends)
        reasons = constraint_reasons(item, rep_mz=rep_mz, ion_mode=ev.get("ion_mode"), measured=measured,
                                     scores=result, trend_entry=trend_entry, th=th)
        if item.pop("_unresolved", False):
            reasons["info"].append("reference_unresolved")
        if reasons["hard"]:
            n_hard += 1
            continue
        mirror = None
        if measured and item["_spectrum"]:
            payload = build_mirror_payload(measured, item["_spectrum"], result["alignment"],
                                           title=item["name"] or "", ms2_tol=scoring["ms2_tol"])
            matched_measured.update(payload.get("matched_measured_mz") or ())
            mirror = {"reference": _keep_most_intense_peaks(
                          payload["reference"], lambda p, m=set(payload.get("matched_mz") or ()): p[0] in m,
                          MIRROR_MAX_PEAKS),
                      "matched_mz": payload.get("matched_mz") or [],
                      "matched_measured_mz": payload.get("matched_measured_mz") or []}
        item.pop("_spectrum")
        ppm = dmz = None
        if item["precursor_mz"] and rep_mz is not None:
            dmz = round((float(rep_mz) - item["precursor_mz"]) * 1000, 2)
            ppm = round((float(rep_mz) - item["precursor_mz"]) / item["precursor_mz"] * 1e6, 2)
        kept.append({**item, "kind": "library", "sum_name": trend.sum_composition(item["name"]),
                     "dmz_mda": dmz, "ppm": ppm,
                     "scores": {k: (round(result[k], 4) if isinstance(result[k], float) else result[k])
                                for k in SCORE_KEYS},
                     "trend": trend_entry, **reasons, "mirror": mirror})
    ranked = rank(kept)[:top_n]
    for index, candidate in enumerate(ranked, 1):
        candidate["candidate_id"] = f"L{index}"
    cut_measured = _keep_most_intense_peaks(measured, lambda p: p[0] in matched_measured, MIRROR_MAX_PEAKS)
    return {"candidates": ranked, "n_hard_removed": n_hard, "measured": cut_measured}
