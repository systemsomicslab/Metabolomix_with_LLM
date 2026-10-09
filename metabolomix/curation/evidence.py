"""スポット単位の証拠収集(I/O)。spec §3・§4。

ファイルの対応は上流の実装と実測で確定している(plan「確定済みの事実」):
- アラインメント `.dcl` は `dcl[MasterAlignmentID]` が代表試料のスキャン。
- `.EIC.aef` のスポット番号 = MasterAlignmentID、試料の並びは `.arf` の行と同順。
- 参照は `record_index == LibraryID` かつ `library_id == AnnotatorID - "_<n>"`。
- Δ は代表試料の行(FileID == RepresentativeFileID)の Mass / RT と参照の差。
deps: arf / arf2 / dcl / eic reader、library.store、analysis.spectral_match、
plots.mirror、msdial.peak_verification、curation.eic_shape / isotope。tools_* は import しない。

**EIC と mirror の座標列は payload だけ間引く(判定は変えない)**。実データ check
(kidney pos, 2196 spots)で JSON が 41.7 MB になり、内訳は EIC 31.0 MB(12 トレース
× 中央値 70 点)・mirror 5.5 MB(測定ピーク中央値 45、最大 3835)だった。
`EIC_MAX_POINTS`(既定 40)は `_downsample_points()` が `spot_shape()` 計算の**後**に
1 サンプルずつ間引く(先頭・末尾・頂点・left/right 最近傍を必ず残し、残りはトレース全体に等間隔に散らす)。
`MIRROR_MAX_PEAKS`(既定 150)は `_cut_mirror_for_payload()` が `build_mirror_payload()`
で満スペクトルから `matched_mz`/`matched_measured_mz`/`labels` を確定させた**後**に
measured/reference それぞれを間引く(一致ピークは必ず残し、残りは強度降順)。
`rescore` はどちらも満スペクトルから計算済みの値をそのまま使う。
"""
from __future__ import annotations

from pathlib import Path

from metabolomix.analysis.spectral_match import match_spectrum
from metabolomix.arf import reader as arf_reader
from metabolomix.arf2.match_results import load_spot_annotations, name_prefix
from metabolomix.arf2.reader import load_isotopic_peaks
from metabolomix.curation import isotope
from metabolomix.curation.eic_shape import spot_shape
from metabolomix.dcl.reader import deserialize_dcl
from metabolomix.eic.reader import read_eic_spot_css1
from metabolomix.msdial.peak_verification import adduct_consistency, mass_error_ppm
from metabolomix.plots.mirror import build_mirror_payload

MAX_SPOTS = 3000
REFERENCE_MZ_WINDOW = 0.05     # 参照 precursor の検索窓(実測の最大差 0.023 Da の倍以上)
EIC_WINDOW_FACTOR = 1.5        # 積分範囲の外側に残す幅(範囲の幅に対する倍率)
_UNANNOTATED = {"", "unknown"}

# --- payload だけを間引く定数(実データ check: pos の JSON が 41.7 MB、内訳は EIC 31.0 MB /
# mirror 5.5 MB。判定(spot_shape・rescore・matched_mz など)は必ず全点/全スペクトルで計算した
# 後に、返す座標列だけを削る。判定を変えたら別のバグになる。---

#: EIC 1 サンプルあたりの payload 点数の上限。spot_shape は間引き前の全点で計算する。
EIC_MAX_POINTS = 40
#: mirror の measured/reference それぞれの payload 点数の上限。matched_mz / matched_measured_mz /
#: labels / rescore は間引き前の全スペクトルで計算する。
MIRROR_MAX_PEAKS = 150


def _spaced_indices(n: int, k: int) -> set[int]:
    """0..n-1 に等間隔な k 点（両端を含む）。`numpy.linspace(0, n-1, k)` の四捨五入の整数版。"""
    if k <= 1:
        return {0}
    return {(i * (n - 1) + (k - 1) // 2) // (k - 1) for i in range(k)}


def _downsample_points(points: list, left: float, right: float,
                       max_points: int = EIC_MAX_POINTS) -> list:
    """EIC の 1 サンプル分の点列を payload 用に間引く(`spot_shape` を計算した**後**に
    呼ぶこと——形状指標は間引き前の全点で確定済みでなければならない)。

    先頭・末尾・頂点(最大強度)・`left`/`right` に最も近い点は必ず残し、残りは
    トレース全体に等間隔に散らした位置で埋める(先頭から詰めると、点数が上限の
    2 倍未満のトレースでピークの右半分が丸ごと落ちる)。等間隔の点数 k は、
    アンカーとの和集合が `max_points` 以下になるまで減らす。強度は整数に丸める
    (x は `_trim()` で既に小数 3 桁に丸め済みなのでそのまま)。点数が上限以下なら
    そのまま(強度の丸めだけ行う)。
    """
    n = len(points)
    if n <= max_points:
        return [[x, round(y)] for x, y in points]
    apex = max(range(n), key=lambda i: points[i][1])
    nearest_left = min(range(n), key=lambda i: abs(points[i][0] - left))
    nearest_right = min(range(n), key=lambda i: abs(points[i][0] - right))
    anchors = {0, n - 1, apex, nearest_left, nearest_right}
    k = max_points
    keep = anchors | _spaced_indices(n, k)
    while len(keep) > max_points and k > 1:
        k -= 1
        keep = anchors | _spaced_indices(n, k)
    return [[points[i][0], round(points[i][1])] for i in sorted(keep)]


def _keep_most_intense_peaks(points: list, is_matched, max_points: int) -> list:
    """一致した点は必ず残し、残りは強度降順で `max_points` まで埋める(m/z 昇順に戻して返す)。
    一致した点だけで `max_points` を超える場合はそれでも全部残す(「一致ピークは必ず残す」を優先)。
    """
    if len(points) <= max_points:
        return points
    matched = [p for p in points if is_matched(p)]
    rest = sorted((p for p in points if not is_matched(p)), key=lambda p: -p[1])
    budget = max(0, max_points - len(matched))
    kept = matched + rest[:budget]
    kept.sort(key=lambda p: p[0])
    return kept


def _cut_mirror_for_payload(mirror: dict | None) -> dict | None:
    """mirror payload の `measured`/`reference` を payload 用に間引く。

    呼ぶ前提: `mirror` は満スペクトルから `build_mirror_payload` で組み立て済みで、
    `matched_mz` / `matched_measured_mz` / `labels` は既に確定している(このまま変えない)。
    一致した測定ピークの判定は `build_mirror_payload` が `ms2_tol` で計算済みの
    `matched_measured_mz` をそのまま使う(ここで許容幅を取り直さない)。
    `scored_peak_count` / `unscored_peak_count` も**間引き前の満スペクトルの値のまま**
    返す——「採点対象外だった点数」という満スペクトル上の意味を保つため、間引き後の
    表示点数には合わせない(payload only の決めごと)。
    """
    if mirror is None:
        return None
    matched_measured = set(mirror.get("matched_measured_mz") or ())
    matched_reference = set(mirror.get("matched_mz") or ())
    mirror = dict(mirror)
    mirror["measured"] = _keep_most_intense_peaks(
        mirror["measured"], lambda p: p[0] in matched_measured, MIRROR_MAX_PEAKS)
    mirror["reference"] = _keep_most_intense_peaks(
        mirror["reference"], lambda p: p[0] in matched_reference, MIRROR_MAX_PEAKS)
    return mirror


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


def arf_rows_by_spot(arf_path) -> dict[int, list[dict]]:
    if arf_path is None:
        return {}
    with open(arf_path, "rb") as handle:
        spots = arf_reader.deserialize(handle)
    return {spot["MasterAlignmentID"]: [arf_reader.alignment_feature_row(row)
                                        for row in spot["AlignedPeakProperties"]]
            for spot in spots}


_arf_rows = arf_rows_by_spot


def _trim(points, left, right):
    margin = (right - left) * EIC_WINDOW_FACTOR
    return [[round(x, 3), round(y, 1)] for x, y in points
            if left - margin <= x <= right + margin]


class UnknownFileIdsError(ValueError):
    """`file_ids` に `.arf` の行に無い試料 ID が含まれる。`missing` はその昇順リスト。"""

    def __init__(self, missing: list[int]):
        self.missing = missing
        super().__init__(f"file_ids のうち {missing} はこのアラインメントの .arf の行にありません。")


def _check_file_ids(rows_by_spot: dict, file_ids) -> None:
    if file_ids is None or not rows_by_spot:
        return
    known = {r.get("file_id") for rows in rows_by_spot.values() for r in rows}
    missing = sorted(set(file_ids) - known)
    if missing:
        raise UnknownFileIdsError(missing)


def _choose_file_ids(rows, representative_file_id, file_ids, max_traces):
    if file_ids is not None:
        return list(file_ids)[:max_traces]
    rows = [r for r in rows if r.get("file_id") is not None]   # 壊れた(空の)行は数えない
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
        library_id=store.library_id_for(match.get("annotator_id")),
        precursor_mz=float(rep_mz), mz_tol=REFERENCE_MZ_WINDOW)
    if record is None:
        return None
    if match.get("inchikey") and record.get("inchikey") and match["inchikey"] != record["inchikey"]:
        return None                        # ScanID が別の行に当たった(別ライブラリ)
    return record


def collect(arf2_path, spots, *, store, ms2_tol, th, file_ids=None, max_traces=12,
            keep_measured: bool = False):
    """スポットごとの証拠を集める。`file_ids` に `.arf` の行に無い ID があれば、重い読み込み
    (`.dcl`・EIC)の前に `UnknownFileIdsError` を投げる。"""
    files = sibling_files(arf2_path)
    rows_by_spot = _arf_rows(files["arf"])
    _check_file_ids(rows_by_spot, file_ids)
    annotations = load_spot_annotations(arf2_path)
    isotopic_peaks = load_isotopic_peaks(arf2_path)
    dcl = (deserialize_dcl(str(files["dcl"]), include_spectrum=True, top_n_peaks=None)
           if files["dcl"] else [])
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
            dmz = float(rep_mz) - reference["precursor_mz"]
            ppm_basis = "reference"
            ref_rt = reference.get("rt")
            drt = float(rep_rt) - float(ref_rt) if ref_rt and ref_rt > 0 else None
        else:
            computed = mass_error_ppm(rep_mz, spot.get("Formula"), spot.get("AdductType"))
            ppm = computed["ppm"]
            theoretical = computed.get("theoretical_mz")
            dmz = float(rep_mz) - theoretical if ppm is not None and theoretical else None
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
            mirror = _cut_mirror_for_payload(mirror)

        samples = []
        if files["eic"] is not None and rows:
            chosen = _choose_file_ids(rows, rep_file, file_ids, max_traces)
            detected = {r["file_id"]: not r.get("is_gap_filled") for r in rows
                        if r.get("file_id") is not None}
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
        for s in samples:                  # payload only: 形状計算の後に間引く(判定は変えない)
            s["points"] = _downsample_points(s["points"], s["left"], s["right"])

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
            "dmz_mda": None if dmz is None else round(dmz * 1000, 2),
            "adduct_band": adduct_consistency(spot.get("AdductType"), spot.get("IonMode"),
                                              spot.get("Ontology"))["band"],
            "drt": None if drt is None else round(drt, 4),
            # MS-DIAL がこの照合で RT を使ったか（False なら RT 判定をしない。None は不明）。
            "rt_used_by_annotation": (store.rt_used_for(match.get("annotator_id"))
                                      if store is not None and match is not None else None),
            "rescore": rescore, "mirror": mirror,
            "eic": {"samples": samples}, "eic_shape": shape, "notes": notes,
            "isotopes": {"measured": isotope.measured_envelope(isotopic_peaks.get(spot_id)),
                         "theoretical": isotope.theoretical_envelope(spot.get("Formula"),
                                                                     spot.get("AdductType"))},
            **({"_measured": measured} if keep_measured else {}),
        })
    return results, stats
