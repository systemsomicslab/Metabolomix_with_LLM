"""Read MS-DIAL AnalysisFileClass values and join them to ARF samples."""

from __future__ import annotations

from collections import Counter
import io
from pathlib import Path
import re
import zipfile

import lz4.block
import msgpack

from metabolomix.msdial.tags import normalize_sample_name
from metabolomix.msdial.sample_factors import (
    arf_sample_names,
    assign_factor_groups,
    build_sample_facets,
    expand_sample_specs,
)


MSGPACK_LZ4_BLOCK_TYPE = 99
_ALIGNMENT_TIMESTAMP_RE = re.compile(r"(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})")
_COMPACT_TIMESTAMP_RE = re.compile(r"(\d{12,14})")

# MsdialDataStorageBase MessagePack keys.
ANALYSIS_FILES_INDEX = 0

# AnalysisFileBean MessagePack keys.
ANALYSIS_FILE_PATH_INDEX = 0
ANALYSIS_FILE_NAME_INDEX = 1
ANALYSIS_FILE_TYPE_INDEX = 2
ANALYSIS_FILE_CLASS_INDEX = 3
ANALYTICAL_ORDER_INDEX = 4
ANALYSIS_FILE_ID_INDEX = 5
ANALYSIS_FILE_INCLUDED_INDEX = 6

# ProjectBaseParameter MessagePack keys.
PROJECT_FOLDER_INDEX = 3
PROJECT_FILE_INDEX = 4

# AlignmentChromPeakFeature MessagePack keys.
ARF_FILE_ID_INDEX = 0
ARF_FILE_NAME_INDEX = 1


def decode_msdial_messagepack(data: bytes):
    """Decode a MessagePack-CSharp Lz4Block payload used by MS-DIAL."""
    packed = msgpack.unpackb(data, raw=False, strict_map_key=False)
    if not isinstance(packed, msgpack.ExtType):
        return packed
    if packed.code != MSGPACK_LZ4_BLOCK_TYPE:
        raise ValueError(f"Unsupported MS-DIAL MessagePack extension type: {packed.code}")

    unpacker = msgpack.Unpacker(io.BytesIO(packed.data), raw=False)
    uncompressed_size = next(unpacker)
    compressed = packed.data[unpacker.tell():]
    unpacked = lz4.block.decompress(
        compressed,
        uncompressed_size=uncompressed_size,
    )
    return msgpack.unpackb(unpacked, raw=False, strict_map_key=False)


def resolve_mddata_path(
    source_path: str | Path,
    mddata_path: str | Path | None = None,
) -> Path | None:
    """Resolve an mddata file explicitly, from mdproject, or beside an ARF."""
    if mddata_path is not None:
        explicit = Path(mddata_path)
        if not explicit.is_file():
            raise FileNotFoundError(f"MS-DIAL dataset file was not found: {explicit}")
        return explicit

    source = Path(source_path)
    if source.suffix.casefold() == ".mddata":
        return source if source.is_file() else None
    if source.suffix.casefold() == ".mdproject":
        candidates = _mddata_paths_from_project(source)
        return _select_single_dataset(candidates, source)

    directory = source.parent if source.suffix else source
    direct = sorted(directory.glob("*.mddata"))
    if direct:
        return _select_single_dataset(direct, source)

    projects = sorted(directory.glob("*.mdproject"))
    project_candidates = []
    for project in projects:
        project_candidates.extend(_mddata_paths_from_project(project))
    return _select_single_dataset(project_candidates, source)


def parse_analysis_file_classes(mddata_path: str | Path) -> list[dict]:
    """Extract AnalysisFileBean metadata, including user-defined Class ID."""
    path = Path(mddata_path)
    storage = decode_msdial_messagepack(path.read_bytes())
    if not isinstance(storage, list) or len(storage) <= ANALYSIS_FILES_INDEX:
        raise ValueError(f"Invalid MS-DIAL mddata structure: {path}")
    analysis_files = storage[ANALYSIS_FILES_INDEX]
    if not isinstance(analysis_files, list):
        raise ValueError(f"AnalysisFiles is not a list in: {path}")

    records = []
    for row in analysis_files:
        if not isinstance(row, list) or len(row) <= ANALYSIS_FILE_ID_INDEX:
            continue
        records.append({
            "file_id": _at(row, ANALYSIS_FILE_ID_INDEX),
            "file_name": _decode_text(_at(row, ANALYSIS_FILE_NAME_INDEX)),
            "file_path": _decode_text(_at(row, ANALYSIS_FILE_PATH_INDEX)),
            "class_id": _decode_text(_at(row, ANALYSIS_FILE_CLASS_INDEX)),
            "analysis_file_type": _at(row, ANALYSIS_FILE_TYPE_INDEX),
            "analytical_order": _at(row, ANALYTICAL_ORDER_INDEX),
            "included": bool(_at(row, ANALYSIS_FILE_INCLUDED_INDEX, True)),
        })
    return records


def discover_arf_class_index(
    arf_path: str | Path,
    mddata_path: str | Path | None = None,
) -> dict | None:
    """Resolve mddata and build FileID/FileName lookup tables for an ARF."""
    resolved = resolve_mddata_path(arf_path, mddata_path)
    if resolved is None:
        return None
    records = parse_analysis_file_classes(resolved)

    by_file_id: dict[int, dict] = {}
    by_file_name: dict[str, dict] = {}
    for record in records:
        file_id = record["file_id"]
        if isinstance(file_id, int):
            if file_id in by_file_id:
                raise ValueError(f"Duplicate AnalysisFileId in {resolved}: {file_id}")
            by_file_id[file_id] = record

        name_key = normalize_sample_name(
            record["file_name"], strip_processing_timestamp=False,
        )
        if name_key:
            if name_key in by_file_name:
                raise ValueError(
                    f"Duplicate AnalysisFileName after normalization in {resolved}: "
                    f"{record['file_name']}"
                )
            by_file_name[name_key] = record

    class_counts = Counter(record["class_id"] for record in records)
    return {
        "mddata_path": str(resolved),
        "records": records,
        "by_file_id": by_file_id,
        "by_file_name": by_file_name,
        "class_counts": dict(sorted(class_counts.items(), key=lambda item: str(item[0]))),
    }


def attach_class_ids_to_spots(features: list[dict], class_index: dict | None) -> list[dict]:
    """Attach per-sample Class ID metadata to every ARF alignment spot."""
    if class_index is None:
        return features
    for spot in features:
        sample_classes = {}
        for row in spot.get("AlignedPeakProperties") or []:
            file_id, file_name = _arf_sample_identity(row)
            record = resolve_sample_class(class_index, file_id, file_name)
            if record is None:
                continue
            sample_classes[_sample_key(file_id, file_name)] = record
        spot["SampleClasses"] = sample_classes
    return features


def _class_tokens(class_id) -> set[str]:
    """Split a Class ID into its underscore-delimited factor tokens (casefold)."""
    return {token for token in str(class_id).casefold().split("_") if token}


def expand_class_specs(
    specs: list[str],
    available_class_ids: list[str],
) -> dict[str, list[str]]:
    """Expand each (possibly partial) Class spec to the matching full Class IDs.

    A spec is split into underscore-delimited tokens; a Class ID matches when it
    contains ALL of the spec's tokens (order-independent AND). A full exact Class
    ID therefore matches only itself. Each spec must match at least one Class ID
    or a ValueError is raised.

    Returns a mapping {original_spec: [matched_class_id, ...]} (original casing
    preserved for both keys and values).
    """
    result: dict[str, list[str]] = {}
    for spec in specs:
        spec_str = str(spec).strip()
        if not spec_str:
            continue
        spec_tokens = _class_tokens(spec_str)
        matched = [
            class_id
            for class_id in available_class_ids
            if spec_tokens <= _class_tokens(class_id)
        ]
        if not matched:
            raise ValueError(
                f"No Class ID matched spec '{spec_str}'. "
                f"Available: {', '.join(str(c) for c in available_class_ids)}"
            )
        result[spec] = matched
    return result


def assign_sample_groups(
    sample_names: list[str],
    class_index: dict | None,
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
) -> dict[str, str | None]:
    """各 PCA サンプル名を、色分け用の群ラベルへ対応づける。

    - 因子指定なし: 群はサンプルの完全 Class ID（class_index が無ければ None）。
    - ``group_levels``（1因子の値トークン、例 ["gf", "spf"]）: その因子だけで統合する。
      どの値にも該当しなければ "other"、2つ以上に該当すれば ValueError。
    - ``group_factors``（因子軸のリスト、例 [["control","ILG"], ["0h","6h"]]）: 軸ごとの
      値を解決して直積ラベル（"control|0h"）にする。group_levels より優先。

    値トークンは Class ID だけでなくサンプル名からも解決される。時点や複製のように
    Class ID に入っていない因子で色分けできるようにするため（詳細は sample_factors）。
    """
    facets = build_sample_facets(sample_names, class_index)
    return assign_factor_groups(facets, group_factors=group_factors, group_levels=group_levels)


def filter_arf_by_class_ids(
    features: list[dict],
    class_index: dict | None,
    class_ids: list[str] | None,
    *,
    missing_sample_policy: str = "error",
    include_roles=("sample",),
) -> tuple[list[dict], dict]:
    """サンプル名 ∪ Class ID の因子トークンでサンプル別ピーク行を絞り込む。

    ``class_ids`` の各要素は `_` 区切りの部分指定（要素内 AND・順不同、要素間 OR）。
    Class ID だけでなくサンプル名のトークンにも一致するため、Class ID に入っていない
    因子（時点・複製・測定日）でも絞り込める。``class_index`` が None（.mddata 未検出）
    でもサンプル名だけで成立する。

    ``include_roles`` は既定 ("sample",) で QC/blank を落とし、内訳を stats の
    ``excluded_by_role`` に残す。行の照合キーは AlignedPeakProperties の FileName。
    """
    specs = [str(value) for value in class_ids or [] if str(value).strip()]
    before_rows = _count_rows(features)
    if not specs:
        return features, {
            "requested_class_ids": [],
            "matched_class_ids": [],
            "before_spots": len(features),
            "after_spots": len(features),
            "before_sample_peaks": before_rows,
            "after_sample_peaks": before_rows,
        }
    if missing_sample_policy not in {"error", "exclude"}:
        raise ValueError("missing_sample_policy must be 'error' or 'exclude'")

    facets = build_sample_facets(arf_sample_names(features), class_index)
    matches, excluded = expand_sample_specs(specs, facets, include_roles=include_roles)
    selected = {name for names in matches.values() for name in names}
    excluded_by_role = sorted({name for names in excluded.values() for name in names})

    # Class メタデータを持つはずなのに解決できないサンプルは、FileID/FileName の
    # 食い違いを示す。従来どおり既定でエラーにする（class_index が無いときは
    # そもそも名前トークンだけで動く設計なので「未対応」ではない）。
    missing_samples: set[str] = set()
    missing_excluded: list[str] = []
    if class_index is not None:
        missing_samples = {name for name, facet in facets.items() if facet.class_id is None}
        if missing_samples and missing_sample_policy == "error":
            raise ValueError(
                "No Class ID metadata matched these ARF samples: "
                + ", ".join(sorted(missing_samples))
            )
        # policy="exclude" は「Class メタデータを解決できないサンプルを解析から外す」の
        # 意。統合トークン空間ではそういうサンプルもサンプル名トークンで spec に一致
        # するため、ここで明示的に差し引かないと「未対応 N 件」と開示しながら実際には
        # 解析に入ってしまう（旧 Class ID 専用実装は行ループで読み飛ばしていた）。
        if missing_samples:
            missing_excluded = sorted(selected & missing_samples)
            selected -= missing_samples

    filtered = []
    for spot in features:
        kept_rows = []
        for row in spot.get("AlignedPeakProperties") or []:
            _, file_name = _arf_sample_identity(row)
            if file_name in selected:
                kept_rows.append(row)
        if kept_rows:
            copied = spot.copy()
            copied["AlignedPeakProperties"] = kept_rows
            filtered.append(copied)

    matched_class_ids = sorted({
        facets[name].class_id for name in selected if facets[name].class_id
    })
    return filtered, {
        "requested_class_ids": specs,
        "matched_class_ids": matched_class_ids,
        "matched_samples": sorted(selected),
        "excluded_by_role": excluded_by_role,
        "before_spots": len(features),
        "after_spots": len(filtered),
        "before_sample_peaks": before_rows,
        "after_sample_peaks": _count_rows(filtered),
        "missing_samples": len(missing_samples),
        "missing_samples_excluded": missing_excluded,
    }


def resolve_sample_class(
    class_index: dict,
    file_id: int | None,
    file_name: str | None,
) -> dict | None:
    """Resolve class metadata by FileID first and FileName second."""
    by_id = class_index.get("by_file_id", {})
    by_name = class_index.get("by_file_name", {})
    id_record = by_id.get(file_id) if file_id is not None else None
    name_key = normalize_sample_name(file_name, strip_processing_timestamp=False)
    name_record = by_name.get(name_key) if name_key else None
    if id_record is not None and name_record is not None and id_record is not name_record:
        raise ValueError(
            f"ARF FileID/FileName resolve to different mddata samples: "
            f"FileID={file_id}, FileName={file_name}"
        )
    return id_record or name_record


def get_sample_class_id(
    spot: dict,
    file_id: int | None,
    file_name: str | None,
) -> str | None:
    record = (spot.get("SampleClasses") or {}).get(_sample_key(file_id, file_name))
    return record.get("class_id") if record else None


def _mddata_paths_from_project(project_path: Path) -> list[Path]:
    if not project_path.is_file():
        return []
    with zipfile.ZipFile(project_path) as archive:
        project_data = decode_msdial_messagepack(archive.read("Project"))
    parameters = project_data.get("ProjectParameters", []) if isinstance(project_data, dict) else []
    candidates = []
    for parameter in parameters:
        if not isinstance(parameter, list) or len(parameter) <= PROJECT_FILE_INDEX:
            continue
        folder = Path(_decode_text(parameter[PROJECT_FOLDER_INDEX]) or project_path.parent)
        filename = Path(_decode_text(parameter[PROJECT_FILE_INDEX]) or "")
        local_candidate = project_path.parent / filename.name
        saved_candidate = folder / filename
        if local_candidate.is_file():
            candidates.append(local_candidate)
        elif saved_candidate.is_file():
            candidates.append(saved_candidate)
    return candidates


def _select_single_dataset(candidates, source: Path) -> Path | None:
    unique = sorted({Path(candidate).resolve() for candidate in candidates if Path(candidate).is_file()})
    if not unique:
        return None
    source_timestamp = _metadata_timestamp_key(source)
    if source_timestamp:
        not_newer_than_source = [
            path for path in unique
            if (candidate_timestamp := _metadata_timestamp_key(path))
            and candidate_timestamp <= source_timestamp
        ]
        if not_newer_than_source:
            return max(not_newer_than_source, key=_metadata_recency_key)
    return max(unique, key=_metadata_recency_key)


def _metadata_recency_key(path: Path) -> tuple[str, float]:
    timestamp = _metadata_timestamp_key(path)
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0
    return timestamp, mtime


def _metadata_timestamp_key(path: Path) -> str:
    match = _ALIGNMENT_TIMESTAMP_RE.search(path.name)
    if match:
        return match.group(1).replace("_", "")
    compact = _COMPACT_TIMESTAMP_RE.search(path.name)
    return compact.group(1) if compact else ""


def _arf_sample_identity(row: list) -> tuple[int | None, str | None]:
    if not isinstance(row, list):
        return None, None
    file_id = row[ARF_FILE_ID_INDEX] if len(row) > ARF_FILE_ID_INDEX else None
    file_name = _decode_text(row[ARF_FILE_NAME_INDEX]) if len(row) > ARF_FILE_NAME_INDEX else None
    return file_id if isinstance(file_id, int) else None, file_name


def _sample_key(file_id: int | None, file_name: str | None) -> int | str:
    if file_id is not None:
        return file_id
    return normalize_sample_name(file_name, strip_processing_timestamp=False)


def _decode_text(value) -> str | None:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore")
    if value is None:
        return None
    return str(value)


def _at(row: list, index: int, default=None):
    return row[index] if len(row) > index else default


def _count_rows(features: list[dict]) -> int:
    return sum(len(spot.get("AlignedPeakProperties") or []) for spot in features)
