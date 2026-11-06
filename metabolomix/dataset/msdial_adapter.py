"""MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）から `DatasetState` を作る。

値・検出状態・試料・極性までを担う。同定の写しは別の Task が同じモジュールへ足す。
"""
from __future__ import annotations

import io
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from metabolomix.analysis import sample_manifest
from metabolomix.arf import reader as arf_reader
from metabolomix.arf2.spot_fields import load_spot_fields
from metabolomix.core.atomic_io import DomainError
from metabolomix.curation.apply import flags_for_arf2
from metabolomix.curation.flags import FlagFileError
from metabolomix.dataset import msdial_identity
from metabolomix.dataset.alignment_files import AlignmentSelection
from metabolomix.msdial.classes import parse_analysis_file_classes
from metabolomix.msdial.tags import parse_tag_file
from metabolomix.mztab import evidence
from metabolomix.mztab.dataset_state import DatasetState, _sha256

__all__ = ["load_msdial_folder"]

_MEASURES = {"peak_height": "height", "peak_area_above_zero": "area"}
_ROLE_BY_FILE_TYPE = {0: "sample", 1: "standard", 2: "qc", 3: "blank"}


def _read_groups(arf_path: Path) -> list[dict]:
    with open(arf_path, "rb") as handle:
        return arf_reader.deserialize(io.BytesIO(handle.read()))


def _cells(group: dict) -> dict[int, dict]:
    cells: dict[int, dict] = {}
    for raw in group.get("AlignedPeakProperties") or []:
        feature = arf_reader.alignment_feature_row(raw)
        if feature and isinstance(feature.get("file_id"), int):
            cells[feature["file_id"]] = feature
    return cells


def _number(value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return math.nan
    return float(value)


def _check_alignment(groups: list[dict], fields: list[dict]) -> None:
    if len(groups) != len(fields):
        raise DomainError(
            "ALIGNMENT_MISMATCH",
            f".arf のスポット数（{len(groups)}）と .arf2 のスポット数（{len(fields)}）が違う。"
            "別のバッチの出力が混ざっている。",
            {"arf_spots": len(groups), "arf2_spots": len(fields)})
    ok, detail = evidence.compare_feature_mz(
        [g.get("MassCenter") for g in groups], [f["mz"] for f in fields])
    if not ok:
        raise DomainError(
            "ALIGNMENT_MISMATCH",
            f".arf と .arf2 のスポットの m/z が合わない（最悪の差 {detail['worst_mz_delta']:.4f}、"
            f"{detail['at_feature_index']} 番目）。別のバッチの出力が混ざっている。",
            detail)


def _ion_mode_majority(fields: list[dict]) -> str | None:
    positive = sum(1 for f in fields if f["ion_mode"] == "Positive")
    negative = sum(1 for f in fields if f["ion_mode"] == "Negative")
    if positive == negative:
        return None
    return "positive" if positive > negative else "negative"


def _clean(name: str | None) -> str | None:
    """表示名の前後の空白を落とす（mzTab-M の MTD 値と同じ扱い）。"""
    return name.strip() if isinstance(name, str) else name


def _msdial_samples(file_ids: list[int], arf_names: dict[int, str | None],
                    mddata: Path | None, notes: list[str]) -> tuple[list[str], list[dict] | None]:
    """列の試料名と、.mddata 由来の試料の行（無ければ None）。名前は .mddata を優先する。"""
    if mddata is None:
        notes.append(".mddata が無いので、試料の役割・群・バッチは宣言から読めず、"
                     "ファイル名からの推定になる。")
        return [_clean(arf_names[i]) or f"Sample_{i}" for i in file_ids], None
    by_id = {r["file_id"]: r for r in parse_analysis_file_classes(mddata)
             if isinstance(r["file_id"], int)}
    names, samples = [], []
    for file_id in file_ids:
        record = by_id.get(file_id)
        if record is None:
            name = _clean(arf_names[file_id]) or f"Sample_{file_id}"
            notes.append(f".mddata に FileID {file_id}（{name}）が無いので、その試料の役割は sample とした。")
            record = {}
        else:
            name = _clean(record["file_name"]) or _clean(arf_names[file_id]) or f"Sample_{file_id}"
            if arf_names[file_id] is not None and arf_names[file_id] != record["file_name"]:
                notes.append(f"FileID {file_id}: .mddata の名前 {record['file_name']!r} と "
                             f".arf の名前 {arf_names[file_id]!r} が違う。FileID で結び、"
                             ".mddata の名前を使った。")
        file_type = record.get("analysis_file_type")
        batch = record.get("analysis_batch")
        names.append(name)
        samples.append({
            "file_id": file_id, "name": name, "arf_file_name": arf_names[file_id],
            "group": record.get("class_id") or None,
            "role": _ROLE_BY_FILE_TYPE.get(file_type, "unknown") if record else "sample",
            "file_type": file_type, "batch": None if batch is None else str(batch),
            "injection_order": record.get("analytical_order"),
            "file_path": record.get("file_path") or None,
        })
    return names, samples


def _fill_identity(ds, fid: str, spot: dict, meta: dict, cache) -> None:
    rep = spot["representative"]
    msms = meta["msms_assigned"]
    inchikey = msdial_identity.resolve_inchikey(spot)
    meta.update({
        "inchikey": inchikey, "inchikey_source": "arf2" if inchikey else None,
        "smiles": spot["smiles"], "inchi": None,
        "annotator_id": (rep or {}).get("annotator_id"),
        "annotation_source": msdial_identity.annotation_source(spot),
        "spectrum_matched": msdial_identity.spectrum_matched(spot),
        "total_score": (rep or {}).get("total_score"),
        "identity": cache.columns(spot, msms_assigned=msms),
    })
    if msdial_identity.is_sme_like(spot, msms_assigned=msms):
        ds.feature_candidates[fid] = [msdial_identity.sme_candidate(spot, fid)]
    annotation = msdial_identity.sml_annotation(spot, msms_assigned=msms)
    if annotation is not None:
        ds.feature_annotations[fid] = annotation


def _load_tags(ds, selection: AlignmentSelection, fields: list[dict]) -> None:
    if selection.tags is None:
        return
    try:
        parsed = parse_tag_file(selection.tags)
    except (ET.ParseError, ValueError, OSError) as exc:
        ds.load_notes.append(f"タグファイルを読めなかったので無視した: {exc}")
        return
    # Peak の Id は MasterAlignmentID（既存の attach_tags_to_spots と同じ）。
    for spot in fields:
        tag_ids = parsed["peaks"].get(spot["master_id"])
        if tag_ids:
            ds.feature_tags[str(spot["master_id"])] = [
                parsed["definitions"].get(t, str(t)) for t in sorted(tag_ids)]


def _load_curation(ds, selection: AlignmentSelection) -> None:
    try:
        ds.curation_flags = flags_for_arf2(selection.arf2)
    except FlagFileError as exc:
        ds.load_notes.append(f"キュレーションのフラグを読めなかったので無視した: {exc}")


def load_msdial_folder(selection: AlignmentSelection, *, measure: str = "peak_height") -> DatasetState:
    if measure not in _MEASURES:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            f"measure は {sorted(_MEASURES)} のどちらか: {measure!r}", {"measure": measure})
    cell_key = _MEASURES[measure]

    groups = _read_groups(selection.arf)
    fields = load_spot_fields(selection.arf2)
    _check_alignment(groups, fields)

    cells_by_spot = [_cells(g) for g in groups]
    file_ids = sorted({fid for cells in cells_by_spot for fid in cells})
    arf_names: dict[int, str | None] = {}
    for cells in cells_by_spot:
        for fid, cell in cells.items():
            arf_names.setdefault(fid, cell.get("file_name"))

    ds = DatasetState()
    ds.source_format = "msdial_folder"
    ds.quantification_measure = measure
    names, samples = _msdial_samples(file_ids, arf_names, selection.mddata, ds.load_notes)
    ds.sample_names = names
    ds.sample_assay_ids = [str(i) for i in file_ids]
    ds.msdial_samples = samples
    ds.ion_mode_majority = _ion_mode_majority(fields)

    n_spots, n_samples = len(groups), len(file_ids)
    matrix = np.full((n_spots, n_samples), np.nan)
    detected = np.zeros((n_spots, n_samples), dtype=bool)
    per_sample_mz: dict[str, dict[str, float]] = {}
    feature_ids, metadata = [], {}
    identity_cache = msdial_identity.IdentityCache()
    # 特徴 ID は .arf2 の MasterID（Key0）。mzTab-M の SMF_ID と同じ番号空間で、
    # イオンモビリティでは連番にならない。
    for i, (cells, spot) in enumerate(zip(cells_by_spot, fields)):
        fid = str(spot["master_id"])
        feature_ids.append(fid)
        per_sample_mz[fid] = {}
        for j, file_id in enumerate(file_ids):
            cell = cells.get(file_id)
            if cell is None:
                continue
            matrix[i, j] = _number(cell.get(cell_key))
            if not cell.get("is_gap_filled"):
                detected[i, j] = True
                mz = _number(cell.get("m_z"))
                if not math.isnan(mz):
                    per_sample_mz[fid][names[j]] = mz
        istd = spot["internal_standard_id"]
        metadata[fid] = {
            "name": spot["name"], "mz": spot["mz"], "rt": spot["rt"],
            "ontology": spot["ontology"], "adduct": spot["adduct"], "formula": spot["formula"],
            "comment": spot["comment"], "fill_percentage": spot["fill_percentage"],
            # MS-DIAL の -1 は「内部標準なし」。
            "internal_standard_id": None if istd is None or istd < 0 else istd,
            "msms_assigned": any(c.get("is_msms") for c in cells.values()),
        }
        _fill_identity(ds, fid, spot, metadata[fid], identity_cache)
    ds.feature_ids = feature_ids
    ds.feature_matrix = matrix
    ds.detected_mask = detected
    ds.per_sample_mz = per_sample_mz
    ds.feature_metadata = metadata
    _load_tags(ds, selection, fields)
    _load_curation(ds, selection)

    for path in (selection.arf, selection.arf2, selection.mddata):
        if path is not None:
            ds.source_files[str(path)] = _sha256(path)
    for path in (selection.dcl, selection.eic):
        if path is not None:
            ds.artifact_sha256[str(path)] = _sha256(path)
    ds.artifact_paths = {
        "dcl": [str(selection.dcl)] if selection.dcl else [],
        "eic": [str(selection.eic)] if selection.eic else [],
        "tags": [str(selection.tags)] if selection.tags else [],
    }

    rows = sample_manifest.resolve_metadata(ds, None)
    sample_manifest.apply_metadata(ds, rows)
    ds.roles = {row["sample_id"]: row["role"] for row in rows}
    ds.sample_meta = {
        row["sample_id"]: {
            "role": row["role"], "batch": row["batch"],
            "batch_source": row["provenance"]["batch"]["source"],
            "run_order": row["injection_order"],
            "run_order_source": row["provenance"]["injection_order"]["source"],
        } for row in rows}
    return ds
