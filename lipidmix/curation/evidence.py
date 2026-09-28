"""スポット単位の証拠収集(I/O)。spec §3・§4。

ファイルの対応は上流の実装と実測で確定している(plan「確定済みの事実」):
- アラインメント `.dcl` は `dcl[MasterAlignmentID]` が代表試料のスキャン。
- `.EIC.aef` のスポット番号 = MasterAlignmentID、試料の並びは `.arf` の行と同順。
- 参照は `record_index == LibraryID` かつ `library_id == AnnotatorID - "_<n>"`。
- Δ は代表試料の行(FileID == RepresentativeFileID)の Mass / RT と参照の差。
deps: arf / arf2 / dcl / eic reader、library.store、analysis.spectral_match、
plots.mirror、msdial.peak_verification、curation.eic_shape。tools_* は import しない。
"""
from __future__ import annotations

from pathlib import Path

from lipidmix.analysis.spectral_match import match_spectrum
from lipidmix.arf import reader as arf_reader
from lipidmix.arf2.match_results import load_spot_annotations, name_prefix
from lipidmix.curation.eic_shape import spot_shape
from lipidmix.dcl.reader import deserialize_dcl
from lipidmix.eic.reader import read_eic_spot_css1
from lipidmix.library.store import library_id_from_annotator
from lipidmix.msdial.peak_verification import adduct_consistency, mass_error_ppm
from lipidmix.plots.mirror import build_mirror_payload

MAX_SPOTS = 3000
REFERENCE_MZ_WINDOW = 0.05     # 参照 precursor の検索窓(実測の最大差 0.023 Da の倍以上)
EIC_WINDOW_FACTOR = 1.5        # 積分範囲の外側に残す幅(範囲の幅に対する倍率)
_UNANNOTATED = {"", "unknown"}


def sibling_files(arf2_path) -> dict:
    arf2_path = Path(arf2_path)
    stem = arf2_path.name[: -len(".arf2")]
    folder = arf2_path.parent

    def existing(name):
        path = folder / name
        return path if path.is_file() else None

    return {"dcl": existing(f"{stem}.dcl"), "eic": existing(f"{stem}.EIC.aef"),
            "arf": existing(f"{stem}_PeakProperties.arf")}


def select_spots(catalog, *, ontology, name_contains) -> list[dict]:
    wanted = {o.strip().lower() for o in ontology} if ontology else None
    needle = name_contains.lower() if name_contains else None
    selected = []
    for spot in catalog:
        name = (spot.get("Name") or "").strip()
        if name.lower() in _UNANNOTATED:
            continue
        if wanted is not None and (spot.get("Ontology") or "").strip().lower() not in wanted:
            continue
        if needle is not None and needle not in name.lower():
            continue
        selected.append(spot)
    return selected


def _arf_rows(arf_path) -> dict[int, list[dict]]:
    if arf_path is None:
        return {}
    with open(arf_path, "rb") as handle:
        spots = arf_reader.deserialize(handle)
    return {spot["MasterAlignmentID"]: [arf_reader.alignment_feature_row(row)
                                        for row in spot["AlignedPeakProperties"]]
            for spot in spots}


def _trim(points, left, right):
    margin = (right - left) * EIC_WINDOW_FACTOR
    return [[round(x, 3), round(y, 1)] for x, y in points
            if left - margin <= x <= right + margin]


def _choose_file_ids(rows, representative_file_id, file_ids, max_traces):
    if file_ids is not None:
        return list(file_ids)[:max_traces]
    detected = sorted((r for r in rows if not r.get("is_gap_filled")),
                      key=lambda r: -(r.get("height") or 0.0))
    chosen = [representative_file_id] if representative_file_id is not None else []
    for row in detected:
        if len(chosen) >= max_traces:
            break
        if row["file_id"] not in chosen:
            chosen.append(row["file_id"])
    for row in rows:                       # 検出が足りなければ gap-fill で埋めて比較対象を残す
        if len(chosen) >= max_traces:
            break
        if row["file_id"] not in chosen:
            chosen.append(row["file_id"])
    return chosen


def _reference(store, match, rep_mz):
    if store is None or match is None or match.get("library_id") is None or rep_mz is None:
        return None
    record = store.record_by_scan_id(
        int(match["library_id"]),
        library_id=library_id_from_annotator(match.get("annotator_id")),
        precursor_mz=float(rep_mz), mz_tol=REFERENCE_MZ_WINDOW)
    if record is None:
        return None
    if match.get("inchikey") and record.get("inchikey") and match["inchikey"] != record["inchikey"]:
        return None                        # ScanID が別の行に当たった(別ライブラリ)
    return record


def collect(arf2_path, spots, *, store, ms2_tol, th, file_ids=None, max_traces=12):
    files = sibling_files(arf2_path)
    annotations = load_spot_annotations(arf2_path)
    dcl = (deserialize_dcl(str(files["dcl"]), include_spectrum=True, top_n_peaks=None)
           if files["dcl"] else [])
    rows_by_spot = _arf_rows(files["arf"])
    stats = {"n_reference_resolved": 0, "n_with_match": 0,
             "missing_files": sorted(k for k, v in files.items() if v is None)}

    results = []
    for spot in spots:
        spot_id = spot["MasterAlignmentID"]
        annotation = annotations.get(spot_id) or {}
        match = annotation.get("representative")
        rep_file = annotation.get("representative_file_id")
        rows = rows_by_spot.get(spot_id, [])
        rep_row = next((r for r in rows if r.get("file_id") == rep_file), None)
        rep_mz = rep_row.get("m_z") if rep_row else spot.get("MassCenter")
        rep_rt = rep_row.get("rt") if rep_row else spot.get("RT")
        notes = []
        if match is not None:
            stats["n_with_match"] += 1

        reference = _reference(store, match, rep_mz)
        if reference is not None:
            stats["n_reference_resolved"] += 1
            ppm = (float(rep_mz) - reference["precursor_mz"]) / reference["precursor_mz"] * 1e6
            ppm_basis = "reference"
            ref_rt = reference.get("rt")
            drt = float(rep_rt) - float(ref_rt) if ref_rt and ref_rt > 0 else None
        else:
            computed = mass_error_ppm(rep_mz, spot.get("Formula"), spot.get("AdductType"))
            ppm = computed["ppm"]
            ppm_basis = "formula" if ppm is not None else None
            drt = None

        measured = []
        if spot_id < len(dcl):
            entry = dcl[spot_id]
            if rep_mz is not None and abs(entry["precursor_mz"] - float(rep_mz)) > 0.01:
                notes.append("dcl_precursor_mismatch")
            else:
                measured = entry["msms_spectrum"]
        rescore = mirror = None
        if measured and reference is not None and reference.get("spectrum"):
            scores = match_spectrum(measured, reference["spectrum"], ms2_tol=ms2_tol)
            rescore = {k: round(scores[k], 4) for k in (
                "weighted_dot_product", "simple_dot_product", "reverse_dot_product",
                "matched_peaks_percentage")}
            mirror = build_mirror_payload(measured, reference["spectrum"], scores["alignment"],
                                          title=spot.get("Name") or "", ms2_tol=ms2_tol)

        samples = []
        if files["eic"] is not None and rows:
            chosen = _choose_file_ids(rows, rep_file, file_ids, max_traces)
            detected = {r["file_id"]: not r.get("is_gap_filled") for r in rows}
            eic = read_eic_spot_css1(str(files["eic"]), spot_id, chosen,
                                     max_traces=max_traces, max_total_points=max_traces * 5000)
            for sample in eic["samples"]:
                samples.append({
                    "file_id": sample["file_id"],
                    "detected": detected.get(sample["file_id"], False),
                    "representative": sample["file_id"] == rep_file,
                    "left": round(sample["peak_left"], 4), "top": round(sample["peak_top"], 4),
                    "right": round(sample["peak_right"], 4),
                    "points": _trim(sample["chromatogram"], sample["peak_left"], sample["peak_right"]),
                })
        shape = spot_shape([{**s, "chromatogram": s["points"]} for s in samples], th)
        shape.pop("per_sample", None)

        results.append({
            "spot_id": spot_id, "name": spot.get("Name"), "ontology": spot.get("Ontology"),
            "adduct": spot.get("AdductType"), "ion_mode": spot.get("IonMode"),
            "mz": spot.get("MassCenter"), "rt": spot.get("RT"),
            "rep_mz": rep_mz, "rep_rt": rep_rt,
            "name_prefix": name_prefix(spot.get("Name")),
            "match": match, "n_candidates": annotation.get("n_candidates", 0),
            "reference": None if reference is None else {
                k: reference.get(k) for k in ("name", "precursor_mz", "rt", "adduct",
                                              "inchikey", "record_index", "library_id")},
            "reference_adduct": None if reference is None else reference.get("adduct"),
            "ppm": None if ppm is None else round(ppm, 2), "ppm_basis": ppm_basis,
            "adduct_band": adduct_consistency(spot.get("AdductType"), spot.get("IonMode"),
                                              spot.get("Ontology"))["band"],
            "drt": None if drt is None else round(drt, 4),
            "rescore": rescore, "mirror": mirror,
            "eic": {"samples": samples}, "eic_shape": shape, "notes": notes,
        })
    return results, stats
