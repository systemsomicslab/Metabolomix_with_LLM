"""Read MS-DIAL ``*_tags.xml`` sidecars and apply them to ARF data."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape, quoteattr

from lipidmix.core.atomic_io import atomic_write_bytes


TAG_MODES = {"any", "all", "none", "not_all"}
TAG_SCOPES = {"sample_peak", "alignment_spot"}
MISSING_SAMPLE_POLICIES = {"error", "exclude", "untagged"}

TAG_XML_SUFFIX = "_tags.xml"
TAG_SUFFIX = "_tags"
PEAK_PROPERTIES_SUFFIX = "_PeakProperties"
PROCESSING_TIMESTAMP_DIGITS = 12
_TIMESTAMP_SUFFIX = re.compile(rf"_\d{{{PROCESSING_TIMESTAMP_DIGITS}}}$")
_COMPACT_TIMESTAMP_RE = re.compile(r"(\d{12,14})")

# AlignmentChromPeakFeature MessagePack keys from the MS-DIAL schema.
FILE_ID_INDEX = 0
FILE_NAME_INDEX = 1
MASTER_PEAK_ID_INDEX = 2
# MS-DIAL は `GetFileNameWithoutExtension` でサンプル名を作るので、剥がす集合は
# `SupportMsRawDataExtension` 全 12 形式を覆う必要がある（`.d` を残すと、`QC01.d`
# 由来の名前が `qc01.d` のままになり `qc01` を名乗るサイドカーと結合できない）。
# 先に長い複合サフィックス（`.wiff.scan` 等）を置く。前方から順に 1 つだけ剥がす。
_KNOWN_SUFFIXES = (
    ".wiff.scan",
    ".timeseries.data",
    ".pai2",
    ".wiff2",
    ".wiff",
    ".mzml",
    ".mzxml",
    ".imzml",
    ".cdf",
    ".abf",
    ".ibf",
    ".raw",
    ".qgd",
    ".lcd",
    ".lrp",
    ".d",
)


def normalize_sample_name(
    value: str | Path | None,
    *,
    strip_processing_timestamp: bool = True,
) -> str:
    """Normalize an ARF FileName or tag-sidecar stem for matching.

    **メモ化されている**。`.arf` では spot × 注入の直積ぶん呼ばれる（393MB の実測で
    881 万回）のに、異なる値は注入数しかない。実処理は `Path(...).name` を通るので
    素で呼ぶと pathlib のパス分解だけで数十秒になる。

    ハッシュ不能な値（壊れた行から来た list など）はキャッシュに載せられないので、
    その場合だけ素の計算へ落とす——例外にはしない（従来は `str(value)` で素通り
    していた）。
    """
    if value is None:
        return ""
    try:
        return _normalize_sample_name_cached(value, strip_processing_timestamp)
    except TypeError:
        return _normalize_sample_name_impl(value, strip_processing_timestamp)


@lru_cache(maxsize=8192)
def _normalize_sample_name_cached(value, strip_processing_timestamp: bool) -> str:
    return _normalize_sample_name_impl(value, strip_processing_timestamp)


def _normalize_sample_name_impl(value, strip_processing_timestamp: bool) -> str:
    name = Path(str(value)).name.strip()
    lower = name.casefold()
    if lower.endswith(TAG_XML_SUFFIX):
        name = name[:-len(TAG_XML_SUFFIX)]
    elif lower.endswith(TAG_SUFFIX):
        name = name[:-len(TAG_SUFFIX)]

    lower = name.casefold()
    for suffix in _KNOWN_SUFFIXES:
        if lower.endswith(suffix):
            name = name[:-len(suffix)]
            break

    if strip_processing_timestamp:
        name = _TIMESTAMP_SUFFIX.sub("", name)
    return name.casefold()


#: メモ化の管理面（テストと、長時間走るセッションでの明示的な解放用）。
normalize_sample_name.cache_info = _normalize_sample_name_cached.cache_info
normalize_sample_name.cache_clear = _normalize_sample_name_cached.cache_clear


def parse_tag_file(path: str | Path) -> dict:
    """Parse one MS-DIAL PeakSpotTags XML sidecar."""
    tag_path = Path(path)
    root = ET.parse(tag_path).getroot()
    if root.tag != "PeakSpotTags":
        raise ValueError(f"Unsupported tag XML root '{root.tag}': {tag_path}")

    definitions: dict[int, str] = {}
    for element in root.findall("./Definitions/Tag"):
        id_text = element.findtext("Id")
        label = (element.findtext("Label") or "").strip()
        if id_text is None:
            continue
        try:
            tag_id = int(id_text)
        except ValueError:
            continue
        definitions[tag_id] = label or str(tag_id)

    peaks: dict[int, frozenset[int]] = {}
    for element in root.findall("./Peaks/Peak"):
        id_text = element.get("Id")
        if id_text is None:
            continue
        try:
            peak_id = int(id_text)
        except ValueError:
            continue
        tag_ids = set()
        for tag_element in element.findall("Tag"):
            try:
                tag_ids.add(int(tag_element.text or ""))
            except ValueError:
                continue
        peaks[peak_id] = frozenset(tag_ids)

    return {
        "path": str(tag_path),
        "definitions": definitions,
        "peaks": peaks,
    }


#: MS-DIAL 5 の `PeakSpotTag.AllTypes()`（上流 `MsdialCore/DataObj/PeakSpotTag.cs`）。
#: `_tags.xml` が無いとき、MS-DIAL と同じ定義でファイルを作る。
MSDIAL_TAG_DEFINITIONS = ((1, "Confirmed"), (2, "Low quality spectrum"), (3, "Misannotation"),
                          (4, "Coelution (mixed spectra)"), (5, "Overannotation"))
MISANNOTATION_TAG_ID = 3


def alignment_tag_path(arf2_path: str | Path) -> Path:
    """アラインメントのタグファイル。上流 `AlignmentResultContainer.Save` と同じく
    アラインメントファイルと同じフォルダの `<stem>_tags.xml`。"""
    path = Path(arf2_path)
    return path.parent / f"{path.stem}{TAG_XML_SUFFIX}"


def _serialize_tag_file(definitions: dict[int, str], peaks: dict[int, frozenset[int]]) -> bytes:
    """上流 `XElement.Save` と同じ見た目（UTF-8 BOM・2 空白字下げ・CRLF・末尾改行なし）。
    MS-DIAL の読み取り（`AlignmentResultContainer.Load`）は `Peak` の `Id` 属性と子 `Tag`
    の整数しか見ないが、GUI が書くものと差分が出ないようにそろえる。"""
    lines = ['<?xml version="1.0" encoding="utf-8"?>', "<PeakSpotTags>", "  <Definitions>"]
    for tag_id, label in sorted(definitions.items()):
        lines += ["    <Tag>", f"      <Id>{tag_id}</Id>", f"      <Label>{escape(label)}</Label>", "    </Tag>"]
    lines.append("  </Definitions>")
    kept = {peak_id: ids for peak_id, ids in peaks.items() if ids}
    if not kept:
        lines.append("  <Peaks />")
    else:
        lines.append("  <Peaks>")
        for peak_id in sorted(kept):
            lines.append(f"    <Peak Id={quoteattr(str(peak_id))}>")
            lines += [f"      <Tag>{tag_id}</Tag>" for tag_id in sorted(kept[peak_id])]
            lines.append("    </Peak>")
        lines.append("  </Peaks>")
    lines.append("</PeakSpotTags>")
    return "\ufeff".encode("utf-8") + "\r\n".join(lines).encode("utf-8")


def update_alignment_tag(path: str | Path, *, tag_id: int, add, remove) -> dict:
    """アラインメントの `_tags.xml` で、`add` のスポットに `tag_id` を付け `remove` から外す。

    他のタグ・定義はそのまま残し、タグが 1 つも残らない `Peak` は書かない（上流と同じ）。
    ファイルが無ければ MS-DIAL と同じ定義で作る。書き込みは原子的。戻り値の
    `added` / `removed` は実際に変わったスポットだけ（昇順）。
    """
    tag_path = Path(path)
    created = not tag_path.exists()
    if created:
        definitions, peaks = dict(MSDIAL_TAG_DEFINITIONS), {}
    else:
        parsed = parse_tag_file(tag_path)
        definitions, peaks = dict(parsed["definitions"]), dict(parsed["peaks"])
    if tag_id not in definitions:
        definitions[tag_id] = dict(MSDIAL_TAG_DEFINITIONS).get(tag_id, str(tag_id))
    added, removed = [], []
    for spot in sorted({int(s) for s in add}):
        if tag_id not in peaks.get(spot, frozenset()):
            peaks[spot] = peaks.get(spot, frozenset()) | {tag_id}
            added.append(spot)
    for spot in sorted({int(s) for s in remove}):
        if tag_id in peaks.get(spot, frozenset()):
            peaks[spot] = peaks[spot] - {tag_id}
            removed.append(spot)
    if created or added or removed:
        atomic_write_bytes(tag_path, _serialize_tag_file(definitions, peaks))
    return {"path": str(tag_path), "added": added, "removed": removed, "created": created}


def resolve_alignment_tag_file(
    arf_path: str | Path,
    tag_directory: str | Path | None = None,
) -> Path | None:
    """Resolve the alignment-level sidecar for a PeakProperties ARF file."""
    path = Path(arf_path)
    stem = path.stem
    if stem.endswith(PEAK_PROPERTIES_SUFFIX):
        stem = stem[:-len(PEAK_PROPERTIES_SUFFIX)]
    directory = Path(tag_directory) if tag_directory else path.parent
    candidates = [directory / f"{stem}{TAG_XML_SUFFIX}", directory / f"{stem}{TAG_SUFFIX}"]
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def discover_arf_tag_index(
    arf_path: str | Path,
    features: list[dict] | None = None,
    tag_directory: str | Path | None = None,
) -> dict:
    """Discover alignment and per-sample tag files adjacent to an ARF file."""
    path = Path(arf_path)
    directory = Path(tag_directory) if tag_directory else path.parent
    alignment_path = resolve_alignment_tag_file(path, directory)

    definitions: dict[int, str] = {}
    alignment_peaks: dict[int, frozenset[int]] = {}
    if alignment_path is not None:
        parsed = parse_tag_file(alignment_path)
        definitions.update(parsed["definitions"])
        alignment_peaks = parsed["peaks"]

    sample_names = _collect_arf_sample_names(features or [])
    exact_sample_names = _group_sample_names(sample_names, strip_processing_timestamp=False)
    relaxed_sample_names = _group_sample_names(sample_names, strip_processing_timestamp=True)
    sample_files: dict[str, dict] = {}
    sample_tag_candidates: dict[str, list[Path]] = {}
    ignored_files: list[str] = []
    tag_paths = set(directory.glob(f"*{TAG_XML_SUFFIX}")) | set(directory.glob(f"*{TAG_SUFFIX}"))
    for tag_path in sorted(tag_paths):
        if alignment_path is not None and tag_path.resolve() == alignment_path.resolve():
            continue
        sample_key = _resolve_tag_sample_key(
            tag_path,
            exact_sample_names,
            relaxed_sample_names,
        )
        if sample_names and sample_key is None:
            ignored_files.append(str(tag_path))
            continue
        if sample_key is None:
            sample_key = normalize_sample_name(
                tag_path.name, strip_processing_timestamp=False,
            )
        sample_tag_candidates.setdefault(sample_key, []).append(tag_path)

    duplicate_sample_tag_files: dict[str, list[str]] = {}
    for sample_key, candidates in sorted(sample_tag_candidates.items()):
        selected = _pick_latest_tag_path(candidates)
        if len(candidates) > 1:
            duplicate_sample_tag_files[sample_key] = [str(path) for path in sorted(candidates)]
        parsed = parse_tag_file(selected)
        definitions.update(parsed["definitions"])
        sample_files[sample_key] = parsed

    index = {
        "arf_path": str(path),
        "tag_directory": str(directory),
        "definitions": definitions,
        "alignment_file": str(alignment_path) if alignment_path else None,
        "alignment_peaks": alignment_peaks,
        "sample_files": sample_files,
        "duplicate_sample_tag_files": duplicate_sample_tag_files,
        "arf_sample_names": sample_names,
        "ignored_files": ignored_files,
        "unmatched_arf_samples": sorted(set(sample_names) - set(sample_files)),
    }
    index["summary"] = summarize_tag_index(index)
    return index


def attach_tags_to_spots(features: list[dict], tag_index: dict) -> list[dict]:
    """Attach alignment tags and non-empty sample-peak tags to ARF spots."""
    definitions = tag_index.get("definitions", {})
    alignment_peaks = tag_index.get("alignment_peaks", {})
    sample_files = tag_index.get("sample_files", {})

    for spot in features:
        alignment_id = spot.get("MasterAlignmentID")
        alignment_ids = alignment_peaks.get(alignment_id, frozenset())
        spot["TagIds"] = sorted(alignment_ids)
        spot["Tags"] = _labels_for_ids(alignment_ids, definitions)

        sample_tags: dict[int | str, dict] = {}
        for row in spot.get("AlignedPeakProperties") or []:
            file_id, file_name, master_peak_id = _sample_row_identity(row)
            parsed = sample_files.get(normalize_sample_name(
                file_name, strip_processing_timestamp=False,
            ))
            if parsed is None or master_peak_id is None:
                continue
            tag_ids = parsed["peaks"].get(master_peak_id, frozenset())
            if not tag_ids:
                continue
            key = _sample_tag_key(file_id, file_name)
            sample_tags[key] = {
                "FileID": file_id,
                "FileName": file_name,
                "MasterPeakID": master_peak_id,
                "TagIds": sorted(tag_ids),
                "Tags": _labels_for_ids(tag_ids, definitions),
                "TagFile": parsed["path"],
            }
        spot["SamplePeakTags"] = sample_tags
    return features


def filter_arf_by_tags(
    features: list[dict],
    tag_index: dict,
    tag_labels: list[str] | None,
    mode: str = "any",
    scope: str = "sample_peak",
    missing_sample_policy: str = "error",
) -> tuple[list[dict], dict]:
    """Filter ARF alignment spots or their sample-level peak rows by tags."""
    if mode not in TAG_MODES:
        raise ValueError(f"tag_mode must be one of {sorted(TAG_MODES)}")
    if scope not in TAG_SCOPES:
        raise ValueError(f"tag_scope must be one of {sorted(TAG_SCOPES)}")
    if missing_sample_policy not in MISSING_SAMPLE_POLICIES:
        raise ValueError(
            f"missing_sample_policy must be one of {sorted(MISSING_SAMPLE_POLICIES)}"
        )

    requested_ids = resolve_tag_ids(tag_index, tag_labels or [])
    if not requested_ids:
        return features, {
            "scope": scope,
            "mode": mode,
            "requested_tags": [],
            "before_spots": len(features),
            "after_spots": len(features),
            "before_sample_peaks": _count_sample_rows(features),
            "after_sample_peaks": _count_sample_rows(features),
            "missing_sample_policy": missing_sample_policy,
        }

    if scope == "alignment_spot":
        if not tag_index.get("alignment_file"):
            raise ValueError(
                "No alignment-level MS-DIAL tag file was found for alignment_spot filtering."
            )
        alignment_peaks = tag_index.get("alignment_peaks", {})
        filtered = [
            spot for spot in features
            if _tag_match(alignment_peaks.get(spot.get("MasterAlignmentID"), frozenset()), requested_ids, mode)
        ]
    else:
        sample_files = tag_index.get("sample_files", {})
        unmatched_samples = tag_index.get("unmatched_arf_samples", [])
        if unmatched_samples and missing_sample_policy == "error":
            names = tag_index.get("arf_sample_names", {})
            display_names = [names.get(key, key) for key in unmatched_samples]
            raise ValueError(
                "No MS-DIAL tag file matched these ARF samples: "
                + ", ".join(display_names)
                + ". Use missing_sample_policy='exclude' or 'untagged' to continue explicitly."
            )
        filtered = []
        for spot in features:
            kept_rows = []
            for row in spot.get("AlignedPeakProperties") or []:
                _, file_name, master_peak_id = _sample_row_identity(row)
                parsed = sample_files.get(normalize_sample_name(
                    file_name, strip_processing_timestamp=False,
                ))
                if parsed is None:
                    if missing_sample_policy == "exclude":
                        continue
                    actual_ids = frozenset()
                else:
                    actual_ids = (
                        parsed["peaks"].get(master_peak_id, frozenset())
                        if master_peak_id is not None else frozenset()
                    )
                if _tag_match(actual_ids, requested_ids, mode):
                    kept_rows.append(row)
            if kept_rows:
                copied = spot.copy()
                copied["AlignedPeakProperties"] = kept_rows
                filtered.append(copied)

    definitions = tag_index.get("definitions", {})
    return filtered, {
        "scope": scope,
        "mode": mode,
        "requested_tags": _labels_for_ids(requested_ids, definitions),
        "before_spots": len(features),
        "after_spots": len(filtered),
        "before_sample_peaks": _count_sample_rows(features),
        "after_sample_peaks": _count_sample_rows(filtered),
        "missing_sample_policy": missing_sample_policy,
    }


def resolve_tag_ids(tag_index: dict, labels: list[str]) -> frozenset[int]:
    definitions = tag_index.get("definitions", {})
    by_label = {label.casefold(): tag_id for tag_id, label in definitions.items()}
    resolved = set()
    unknown = []
    for label in labels:
        text = str(label).strip()
        if not text:
            continue
        if text.isdigit() and int(text) in definitions:
            resolved.add(int(text))
        elif text.casefold() in by_label:
            resolved.add(by_label[text.casefold()])
        else:
            unknown.append(text)
    if unknown:
        available = ", ".join(definitions.values()) or "none"
        raise ValueError(f"Unknown MS-DIAL tag(s): {', '.join(unknown)}. Available: {available}")
    return frozenset(resolved)


def summarize_tag_index(tag_index: dict) -> dict:
    definitions = tag_index.get("definitions", {})
    alignment_peaks = tag_index.get("alignment_peaks", {})
    sample_files = tag_index.get("sample_files", {})
    alignment_counts = Counter()
    sample_counts = Counter()
    for tag_ids in alignment_peaks.values():
        alignment_counts.update(tag_ids)
    for parsed in sample_files.values():
        for tag_ids in parsed["peaks"].values():
            sample_counts.update(tag_ids)
    return {
        "definitions": [
            {
                "id": tag_id,
                "label": label,
                "alignment_spots": alignment_counts[tag_id],
                "sample_peaks": sample_counts[tag_id],
            }
            for tag_id, label in sorted(definitions.items())
        ],
        "sample_tag_files": len(sample_files),
        "matched_arf_samples": len(set(sample_files) & set(tag_index.get("arf_sample_names", {}))),
        "arf_samples": len(tag_index.get("arf_sample_names", {})),
        "alignment_tag_file": tag_index.get("alignment_file"),
        "tagged_alignment_spots": len(alignment_peaks),
        "tagged_sample_peaks": sum(len(parsed["peaks"]) for parsed in sample_files.values()),
        "ignored_tag_files": len(tag_index.get("ignored_files", [])),
        "duplicate_sample_tag_files": len(tag_index.get("duplicate_sample_tag_files", {})),
        "unmatched_arf_samples": len(tag_index.get("unmatched_arf_samples", [])),
    }


def _tag_recency_key(path: Path) -> tuple[str, float]:
    match = _COMPACT_TIMESTAMP_RE.search(path.name)
    timestamp = match.group(1) if match else ""
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0
    return timestamp, mtime


def _pick_latest_tag_path(paths: list[Path]) -> Path:
    return max(paths, key=_tag_recency_key)


def _collect_arf_sample_names(features: list[dict]) -> dict[str, str]:
    names: dict[str, str] = {}
    for spot in features:
        for row in spot.get("AlignedPeakProperties") or []:
            _, file_name, _ = _sample_row_identity(row)
            key = normalize_sample_name(file_name, strip_processing_timestamp=False)
            if key:
                if key in names and names[key] != str(file_name):
                    raise ValueError(
                        f"ARF contains sample names that collide after normalization: "
                        f"{names[key]!r}, {file_name!r}"
                    )
                names[key] = str(file_name)
    return names


def _group_sample_names(
    sample_names: dict[str, str],
    *,
    strip_processing_timestamp: bool,
) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for exact_key, display_name in sample_names.items():
        key = normalize_sample_name(
            display_name,
            strip_processing_timestamp=strip_processing_timestamp,
        )
        groups.setdefault(key, []).append(exact_key)
    return groups


def _resolve_tag_sample_key(
    tag_path: Path,
    exact_sample_names: dict[str, list[str]],
    relaxed_sample_names: dict[str, list[str]],
) -> str | None:
    exact_key = normalize_sample_name(
        tag_path.name, strip_processing_timestamp=False,
    )
    exact_matches = exact_sample_names.get(exact_key, [])
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise ValueError(f"Ambiguous exact sample match for tag file: {tag_path}")

    relaxed_key = normalize_sample_name(tag_path.name, strip_processing_timestamp=True)
    relaxed_matches = relaxed_sample_names.get(relaxed_key, [])
    if len(relaxed_matches) == 1:
        return relaxed_matches[0]
    if len(relaxed_matches) > 1:
        raise ValueError(
            f"Ambiguous timestamp-stripped sample match for tag file {tag_path}: "
            + ", ".join(relaxed_matches)
        )
    return None


def _sample_row_identity(row: list) -> tuple[int | None, str | None, int | None]:
    if not isinstance(row, list):
        return None, None, None
    file_id = (
        row[FILE_ID_INDEX]
        if len(row) > FILE_ID_INDEX and isinstance(row[FILE_ID_INDEX], int)
        else None
    )
    file_name = row[FILE_NAME_INDEX] if len(row) > FILE_NAME_INDEX else None
    if isinstance(file_name, bytes):
        file_name = file_name.decode("utf-8", errors="ignore")
    elif not isinstance(file_name, str):
        file_name = None
    master_peak_id = (
        row[MASTER_PEAK_ID_INDEX]
        if len(row) > MASTER_PEAK_ID_INDEX and isinstance(row[MASTER_PEAK_ID_INDEX], int)
        else None
    )
    return file_id, file_name, master_peak_id


def _sample_tag_key(file_id: int | None, file_name: str | None) -> int | str:
    if file_id is not None:
        return file_id
    return normalize_sample_name(file_name, strip_processing_timestamp=False)


def get_sample_peak_tag_info(
    spot: dict,
    file_id: int | None,
    file_name: str | None,
) -> dict:
    """Return attached sample-peak tag metadata using the same fallback key."""
    return (spot.get("SamplePeakTags") or {}).get(
        _sample_tag_key(file_id, file_name), {}
    )


def _labels_for_ids(tag_ids, definitions: dict[int, str]) -> list[str]:
    return [definitions.get(tag_id, str(tag_id)) for tag_id in sorted(tag_ids)]


def _tag_match(actual_ids, requested_ids: frozenset[int], mode: str) -> bool:
    actual = set(actual_ids)
    requested = set(requested_ids)
    if mode == "any":
        return bool(actual & requested)
    if mode == "all":
        return requested <= actual
    if mode == "none":
        return not bool(actual & requested)
    return not requested <= actual


def _count_sample_rows(features: list[dict]) -> int:
    return sum(len(spot.get("AlignedPeakProperties") or []) for spot in features)
