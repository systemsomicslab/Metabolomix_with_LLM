"""MsScanMatchResult（`.arf2` Key 56 の中身）の復号と、GUI が表示する代表の選択。

Key 番号の正準は docs/schema/MsScanMatchResult.md。代表の選択は上流
`MsScanMatchResultContainer.Representative`（`ResultOrder` の argmax）の移植で、
シリアライズ済みの値だけから決まるので忠実に再現できる。
deps: arf2.reader / msdial.lipid_identity（接頭辞の正規表現）。tools_* は import しない。
"""
from __future__ import annotations

import math

from lipidmix.arf2.reader import load_raw_spots
from lipidmix.msdial.lipid_identity import _MSDIAL_QUALIFIER_RE

MATCH_KEYS = {
    0: "name", 1: "inchikey", 2: "total_score",
    3: "squared_weighted_dot_product", 4: "squared_simple_dot_product",
    5: "squared_reverse_dot_product", 6: "matched_peaks_count",
    7: "matched_peaks_percentage", 8: "essential_fragment_matched_score",
    9: "rt_similarity", 10: "ri_similarity", 11: "ccs_similarity",
    12: "isotope_similarity", 13: "accurate_mass_similarity", 14: "library_id",
    15: "is_precursor_mz_match", 16: "is_spectrum_match", 17: "is_rt_match",
    18: "is_ccs_match", 19: "is_lipid_class_match", 20: "is_lipid_chains_match",
    21: "is_lipid_position_match", 22: "is_other_lipid_match", 23: "is_ri_match",
    24: "library_id_when_ordered", 26: "source", 27: "annotator_id",
    28: "spectrum_id", 29: "andromeda_score", 30: "is_decoy", 31: "priority",
    32: "pep_score", 33: "is_reference_matched", 34: "is_annotation_suggested",
    35: "is_lipid_double_bond_position_match", 36: "collision_energy",
    37: "enhanced_dot_product", 38: "spectral_entropy",
}

SOURCE_UNKNOWN = 1
SOURCE_MANUAL = 64

# Mapping from normalized qualifier text to label.
_QUALIFIER_LABELS = {
    "noms2": "no MS2",
    "lowscore": "low score",
    "unsettled": "unsettled",
    "woms2": "w/o MS2",
    "woms1": "w/o MS2",  # Map w/o MS1 to w/o MS2 (same family)
}


def _clean(value):
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def decode_match_result(raw) -> dict | None:
    """MsScanMatchResult の生配列を名前付き dict にする。読めなければ None。

    追加で 2 つの導出値を持つ:
    - `is_manually_modified`: `(Source & Manual) != 0`（上流の定義）
    - `has_msms`: 重み付きドット積が -1（比較不能の番兵）でない
    """
    if not isinstance(raw, list) or len(raw) < 35:
        return None
    decoded = {name: _clean(raw[key]) if key < len(raw) else None
               for key, name in MATCH_KEYS.items()}
    source = decoded.get("source") or 0
    decoded["is_manually_modified"] = bool(int(source) & SOURCE_MANUAL)
    weighted = decoded.get("squared_weighted_dot_product")
    decoded["has_msms"] = weighted is not None and float(weighted) >= 0
    return decoded


def _order_key(result: dict) -> tuple:
    return (
        bool(result.get("is_manually_modified")),
        bool(result.get("is_reference_matched")),
        bool(result.get("is_annotation_suggested")),
        int(result.get("priority") or 0),
        float(result.get("total_score") or 0.0),
    )


def _candidates(container) -> list[dict]:
    if not isinstance(container, list) or not container or not isinstance(container[0], list):
        return []
    decoded = [decode_match_result(raw) for raw in container[0]]
    return [d for d in decoded if d is not None]


def representative(container) -> dict | None:
    """GUI が表示する代表（非 decoy・Source != Unknown の ResultOrder 最大）。無ければ None。"""
    usable = [c for c in _candidates(container)
              if not c.get("is_decoy") and (c.get("source") or 0) != SOURCE_UNKNOWN]
    if not usable:
        return None
    return max(usable, key=_order_key)


def name_prefix(name) -> str | None:
    """Extract and normalize MS-DIAL qualifier prefix from name.

    Uses the canonical regex from lipid_identity._MSDIAL_QUALIFIER_RE to detect
    qualifiers like "no MS2:", "low score:", etc. Returns normalized label or None.
    """
    if not name:
        return None
    s = str(name).strip()
    m = _MSDIAL_QUALIFIER_RE.match(s)
    if not m:
        return None
    # Extract the matched text (including colon) and normalize it
    matched = m.group(0)
    # Remove leading/trailing whitespace and colon, then normalize by removing spaces and slashes
    normalized = matched.strip().rstrip(":").replace(" ", "").replace("/", "").lower()
    # Map to canonical label
    return _QUALIFIER_LABELS.get(normalized)


def spot_annotation(raw_spot: list) -> dict:
    container = raw_spot[56] if len(raw_spot) > 56 else None
    return {
        "spot_id": int(raw_spot[0]),
        "representative_file_id": int(raw_spot[3]) if raw_spot[3] is not None else None,
        "representative": representative(container),
        "n_candidates": len(_candidates(container)),
    }


def load_spot_annotations(file_path) -> dict[int, dict]:
    """`.arf2` の全スポットについて代表の照合結果を MasterAlignmentID で引ける形にする。"""
    return {entry["spot_id"]: entry
            for entry in (spot_annotation(raw) for raw in load_raw_spots(file_path))}
