"""MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）から `DatasetState` を作る。

値・検出状態・試料・極性までを担う。同定の写しは別の Task が同じモジュールへ足す。
"""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path, PureWindowsPath

from metabolomix.analysis import sample_manifest
from metabolomix.arf2.spot_fields import load_spot_fields
from metabolomix.core.atomic_io import DomainError
from metabolomix.curation.apply import flags_for_arf2
from metabolomix.curation.flags import FlagFileError
from metabolomix.dataset import msdial_identity
from metabolomix.dataset.alignment_files import AlignmentSelection
from metabolomix.msdial.classes import MDDATA_READ_ERRORS, parse_analysis_file_classes
from metabolomix.msdial.tags import parse_tag_file
from metabolomix.mztab import evidence
from metabolomix.mztab.dataset_state import DatasetState, _sha256

__all__ = ["MEASURES", "load_msdial_folder"]

MEASURES = {"peak_height": "height", "peak_area_above_zero": "area"}
_ROLE_BY_FILE_TYPE = {0: "sample", 1: "standard", 2: "qc", 3: "blank"}


def _check_alignment(scan: evidence.ArfScan, fields: list[dict]) -> None:
    """`.arf` のスポット（`scan`）と `.arf2` のスポットが同じアラインメントかを確かめる。"""
    if scan.n_spots != len(fields):
        raise DomainError(
            "ALIGNMENT_MISMATCH",
            f".arf のスポット数（{scan.n_spots}）と .arf2 のスポット数（{len(fields)}）が違う。"
            "同じアラインメントの出力でないか、.arf を読み切れていない。",
            {"arf_spots": scan.n_spots, "arf2_spots": len(fields)})
    ok, detail = evidence.compare_feature_mz(scan.spot_mz, [f["mz"] for f in fields])
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
    records = None
    if mddata is not None:
        try:
            records = parse_analysis_file_classes(mddata)
        except MDDATA_READ_ERRORS as exc:
            # 読めない・解けない .mddata は、無いときと同じに続ける（宣言とは名乗らない）。
            notes.append(f".mddata（{mddata.name}）を読めなかった（{type(exc).__name__}: {exc}）。"
                         ".mddata が無いときと同じに、試料の役割・群・バッチはファイル名からの"
                         "推定にした。")
    if records is None:
        if mddata is None:
            notes.append(".mddata が無いので、試料の役割・群・バッチは宣言から読めず、"
                         "ファイル名からの推定になる。")
        return [_clean(arf_names[i]) or f"Sample_{i}" for i in file_ids], None
    by_id = {r["file_id"]: r for r in records if isinstance(r["file_id"], int)}
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
            "in_mddata": bool(record),
        })
    return names, samples


def _annotator_file_name(annotator_id):
    """注釈器 ID はライブラリのフルパスのことがある（研究室ライブラリの置き場所）。名前だけ残す。"""
    if not annotator_id:
        return annotator_id
    return PureWindowsPath(str(annotator_id)).name or annotator_id


def _fill_identity(ds, fid: str, spot: dict, meta: dict, cache) -> None:
    rep = spot["representative"]
    msms = meta["msms_assigned"]
    inchikey = msdial_identity.resolve_inchikey(spot)
    meta.update({
        "inchikey": inchikey, "inchikey_source": "arf2" if inchikey else None,
        "smiles": spot["smiles"], "inchi": None,
        "annotator_id": _annotator_file_name((rep or {}).get("annotator_id")),
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


def _fill_summaries(ds, selection: AlignmentSelection, detected, metadata: dict) -> None:
    """検出状態と InChIKey 被覆の要約（mzTab-M のアダプタの `feature_qc` / `inchikey_coverage` と同じキー）。

    検出状態は `.arf`（`is_gap_filled`）由来なので `source` は "arf"（mzTab-M のアダプタが
    隣の `.arf` から取り込んだときと同じ値。キュレーションの適用の「接合済み」判定がそのまま
    成り立つ）。特徴の軸は `.arf2` の MasterAlignmentID で、フラグの ID 空間と同じ。
    """
    from metabolomix.mztab import identity as mztab_identity

    n_cells, n_detected = int(detected.size), int(detected.sum())
    ds.feature_qc = {
        "source": "arf", "tried": [str(selection.arf)],
        "n_cells": n_cells, "n_detected": n_detected,
        "gap_filled_rate": round(1.0 - n_detected / n_cells, 4) if n_cells else None,
        "join": {"path": str(selection.arf), "feature_axis": "arf2_master_id"},
    }
    with_ik = sum(1 for m in metadata.values() if m.get("inchikey"))
    sme = len(ds.feature_candidates)
    sml_only = sum(1 for fid in ds.feature_annotations if fid not in ds.feature_candidates)
    ds.inchikey_coverage = {
        "total_features": len(metadata), "with_inchikey": with_ik,
        "by_source": {"arf2": with_ik, "none": len(metadata) - with_ik},
        "rdkit_available": mztab_identity.rdkit_available(),
        "msp_lookup": {"status": "not_applicable"},
        "identified_by": {"sme": sme, "sml_only": sml_only,
                          "none": len(metadata) - sme - sml_only},
    }


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


def _usable_scan(scan: evidence.ArfScan | None, selection: AlignmentSelection,
                  cell_key: str) -> evidence.ArfScan | None:
    """渡された `ArfScan` が、選んだ `.arf` を要る定量値つきで読んだものなら使う。"""
    if scan is None or cell_key not in scan.values:
        return None
    try:
        same = Path(scan.path).resolve() == Path(selection.arf).resolve()
    except OSError:
        same = False
    return scan if same else None


def load_msdial_folder(selection: AlignmentSelection, *, measure: str = "peak_height",
                       scan: evidence.ArfScan | None = None) -> DatasetState:
    """解析フォルダの `.arf` / `.arf2`（と `.mddata` ほか）から `DatasetState` を作る。

    `scan` は `selection.arf` を読んだ `evidence.scan_arf` の結果（job の経路が照合に読んだもの）。
    渡れば `.arf` を読まない（別のファイルや定量値の欠けたものなら読み直す）。
    """
    if measure not in MEASURES:
        raise DomainError(
            "DATASET_BAD_REQUEST",
            f"measure は {sorted(MEASURES)} のどちらか: {measure!r}", {"measure": measure})
    cell_key = MEASURES[measure]

    scan = _usable_scan(scan, selection, cell_key)
    if scan is None:
        scan = evidence.scan_arf(selection.arf, measures=(cell_key,))
    try:
        return _build_dataset(selection, scan, measure, cell_key)
    except MemoryError as exc:
        # 読み取りのあと、行列・per_sample_mz（セルあたり約 100 バイト）を組む間に尽きた。
        raise evidence.too_large_error(
            Path(selection.arf), "DatasetState を組む間にメモリが尽きました") from exc


def _build_dataset(selection: AlignmentSelection, scan: evidence.ArfScan, measure: str,
                   cell_key: str) -> DatasetState:
    """`scan`（`.arf` の読み取り結果）と `.arf2` ほかから `DatasetState` を組む。"""
    fields = load_spot_fields(selection.arf2)
    _check_alignment(scan, fields)

    file_ids = scan.file_ids
    arf_names = scan.file_names

    ds = DatasetState()
    ds.source_format = "msdial_folder"
    if selection.mddata_error:
        ds.load_notes.append(
            f".mddata を決められなかった（{selection.mddata_error}）。.mddata が無いときと同じに続ける。")
    ds.quantification_measure = measure
    names, samples = _msdial_samples(file_ids, arf_names, selection.mddata, ds.load_notes)
    ds.sample_names = names
    ds.sample_assay_ids = [str(i) for i in file_ids]
    ds.msdial_samples = samples
    ds.ion_mode_majority = _ion_mode_majority(fields)

    matrix = scan.values[cell_key]
    detected = scan.detected
    per_sample_mz: dict[str, dict[str, float]] = {}
    feature_ids, metadata = [], {}
    identity_cache = msdial_identity.IdentityCache()
    # 特徴 ID は .arf2 の MasterID（Key0）。mzTab-M の SMF_ID と同じ番号空間で、
    # イオンモビリティでは連番にならない。
    for i, spot in enumerate(fields):
        fid = str(spot["master_id"])
        feature_ids.append(fid)
        # 検出されたセルの実測 m/z（数でないものと gap-fill は scan の時点で NaN）。
        per_sample_mz[fid] = {names[j]: float(mz) for j, mz in enumerate(scan.per_sample_mz[i])
                              if not math.isnan(mz)}
        istd = spot["internal_standard_id"]
        metadata[fid] = {
            "name": spot["name"], "mz": spot["mz"], "rt": spot["rt"],
            "ontology": spot["ontology"], "adduct": spot["adduct"], "formula": spot["formula"],
            "comment": spot["comment"], "fill_percentage": spot["fill_percentage"],
            # MS-DIAL の -1 は「内部標準なし」。
            "internal_standard_id": None if istd is None or istd < 0 else istd,
            "msms_assigned": bool(scan.msms_assigned[i]),
        }
        _fill_identity(ds, fid, spot, metadata[fid], identity_cache)
    ds.feature_ids = feature_ids
    ds.feature_matrix = matrix
    ds.detected_mask = detected
    ds.per_sample_mz = per_sample_mz
    ds.feature_metadata = metadata
    _fill_summaries(ds, selection, detected, metadata)
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

    sample_manifest.apply_defaults_and_inputs(ds, sample_manifest.resolve_metadata(ds, None))
    return ds
