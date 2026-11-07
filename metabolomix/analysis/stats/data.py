"""DatasetState からの統計エンジンの入口と、検出・gap-fill・欠測の 3 状態（spec §5・§6）。

このモジュールは DatasetState を**読むだけ**で書き換えない。MCP・セッションには依存しない。
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

from metabolomix.core.atomic_io import DomainError

#: セルの状態。DETECTED = 実検出、GAP_FILLED = 値はあるが検出は偽（gap-fill）、
#: MISSING = 値が無い（0 以下・非有限）、PRESENT_UNKNOWN = 値はあるが検出の別が分からない。
DETECTED, GAP_FILLED, MISSING, PRESENT_UNKNOWN = 1, 2, 0, 3

_CORE_KEYS = frozenset({"sample_id", "role", "group", "batch", "injection_order",
                        "include", "subject_id", "biological_sample_id", "source_file"})


@dataclass(frozen=True)
class SampleInfo:
    """試料（注入）1 本の情報。`attributes` は因子・共変量に使う追加の欄。"""
    sample_id: str
    role: str
    group: str | None
    batch: str | None
    injection_order: float | None
    subject_id: str | None
    include: bool
    attributes: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class QuantState:
    """定量値が内部標準などで換算済みか（B-3 の `DatasetState.quant_state` から読む）。"""
    normalized: bool
    method: str | None
    unit: str | None


@dataclass
class StudyData:
    """n_samples × n_features の生の値と、試料の表・検出の別をまとめた入口。"""
    values: np.ndarray
    detected: np.ndarray | None
    sample_ids: list[str]
    feature_ids: list[str]
    samples: list[SampleInfo]
    quant: QuantState | None = None
    dataset_hash: str | None = None
    feature_classes: dict[str, str] = field(default_factory=dict)
    annotated_features: frozenset[str] = frozenset()


def quant_state_of(ds) -> QuantState | None:
    """`ds.quant_state`（lipidmix.quant-state.v1 の dict）から読む。無い・None なら None。"""
    raw = getattr(ds, "quant_state", None)
    if not isinstance(raw, dict):
        return None
    return QuantState(normalized=bool(raw.get("normalized", False)),
                      method=raw.get("method"), unit=raw.get("unit"))


def cell_states(values, detected) -> np.ndarray:
    """セルごとの状態（int8）。有限の正値だけが「値あり」で、それ以外は MISSING。"""
    v = np.asarray(values, dtype=float)
    present = np.isfinite(v) & (v > 0)
    if detected is None:
        return np.where(present, PRESENT_UNKNOWN, MISSING).astype(np.int8)
    d = np.asarray(detected, dtype=bool)
    out = np.full(v.shape, MISSING, dtype=np.int8)
    out[present & d] = DETECTED
    out[present & ~d] = GAP_FILLED
    return out


# 反復注入らしい試料名（plan §15）。茎が数字で終わる 2 段の番号、または inj/rep/tech 接尾辞。
_REPLICATE_A = re.compile(r"^(.*\d)[_\-.](?:inj|rep|tech|r)?\d{1,2}$")
_REPLICATE_B = re.compile(r"(?i)[_\-.](inj|rep|tech)\d+$")


def replicate_name_warnings(sample_ids: list[str]) -> list[dict]:
    """同じ茎の名前が 2 つ以上ある反復注入らしい試料を注意として返す（まとめはしない）。"""
    stems: dict[str, list[str]] = {}
    for name in sample_ids:
        m = _REPLICATE_A.match(name)
        if m:
            stem = m.group(1)
        else:
            m = _REPLICATE_B.search(name)
            if not m:
                continue
            stem = name[:m.start()]
        stems.setdefault(stem, []).append(name)
    return [{"code": "possible_technical_replicates", "stem": stem, "samples": names}
            for stem, names in stems.items() if len(names) >= 2]


def _str_or_none(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _float_or_none(value) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _include(row: dict) -> bool:
    """`include` の読み方は analysis_spec._load と同じ: 欄が無ければ真、あれば Python の真偽値。

    （文字列 "false" や "0" も真になる。sample_metadata_rows の include は取り込み時に bool へ
    正規化済みという前提で、ここで独自に解釈し直して _load と食い違わせない。）
    """
    return bool(row.get("include", True))


def _role(value) -> str | None:
    text = _str_or_none(value)
    return text.lower() if text else None


def _samples_table(ds, names: list[str]) -> list[SampleInfo]:
    """試料の表。sample_metadata_rows → msdial_samples → roles / sample_meta の順に欄ごと埋める。"""
    rows = list(getattr(ds, "sample_metadata_rows", None) or [])
    msdial = list(getattr(ds, "msdial_samples", None) or [])
    roles = getattr(ds, "roles", None) or {}
    meta = getattr(ds, "sample_meta", None) or {}
    by_name = {s.get("name"): s for s in msdial}
    display = list(getattr(ds, "sample_names", None) or [])
    detected_roles = None
    out = []
    for i, name in enumerate(names):
        # 安定 ID を先に、次に表示名（sample_names[i]）で引く。ID ≠ 表示名でも取りこぼさない。
        keys = [name] + ([display[i]] if i < len(display) and display[i] != name else [])
        row = rows[i] if i < len(rows) else {}
        ms = (next((by_name[k] for k in keys if k in by_name), None)
              or (msdial[i] if len(msdial) == len(names) else {}))
        sm = next((meta[k] for k in keys if meta.get(k)), {})
        role = (_role(row.get("role")) or _role(ms.get("role"))
                or next((r for k in keys if (r := _role(roles.get(k)))), None)
                or _role(sm.get("role")))
        if role is None:
            # どこにも役割が無い試料は名前から推定する（msdial/sample_factors と同じ規則）。
            if detected_roles is None:
                from metabolomix.analysis.preprocessing import detect_sample_roles
                detected_roles = detect_sample_roles(list(names))
            role = _role(detected_roles.get(name)) or "sample"
        group = _str_or_none(row.get("group")) or _str_or_none(ms.get("group"))
        batch = (_str_or_none(row.get("batch")) or _str_or_none(ms.get("batch"))
                 or _str_or_none(sm.get("batch")))
        order = next((o for o in (_float_or_none(row.get("injection_order")),
                                  _float_or_none(ms.get("injection_order")),
                                  _float_or_none(sm.get("run_order"))) if o is not None), None)
        subject = _str_or_none(row.get("subject_id")) or _str_or_none(row.get("biological_sample_id"))
        attrs = {k: s for k, v in row.items() if k not in _CORE_KEYS
                 and (s := _str_or_none(v)) is not None}
        out.append(SampleInfo(sample_id=name, role=role, group=group, batch=batch,
                              injection_order=order, subject_id=subject,
                              include=_include(row), attributes=attrs))
    return out


def _fingerprint(ds) -> str | None:
    try:
        from metabolomix.analysis.result_state import dataset_fingerprint
        return dataset_fingerprint(ds)
    except Exception:
        return None


def from_dataset(ds) -> StudyData:
    """DatasetState から StudyData を作る。ds は書き換えない（配列も写す）。"""
    values = np.array(np.asarray(ds.feature_matrix, dtype=float).T, copy=True)
    mask = getattr(ds, "detected_mask", None)
    if mask is not None and np.shape(mask) != values.T.shape:
        raise DomainError("STUDY_DATA_INVALID", "detected_mask の形が feature_matrix と合わない",
                          {"feature_matrix_shape": list(values.T.shape),
                           "detected_mask_shape": list(np.shape(mask))})
    detected = None if mask is None else np.array(np.asarray(mask, dtype=bool).T, copy=True)
    rows = getattr(ds, "sample_metadata_rows", None)
    ids = list(getattr(ds, "sample_ids", None) or []) if rows else []
    names = ids or list(ds.sample_names)
    classes = {}
    for fid, meta in (getattr(ds, "feature_metadata", None) or {}).items():
        ontology = (meta or {}).get("ontology")
        if ontology:
            classes[fid] = str(ontology)
    annotated = frozenset(fid for fid, cands in (getattr(ds, "feature_candidates", None) or {}).items()
                          if cands)
    return StudyData(values=values, detected=detected, sample_ids=names,
                     feature_ids=list(ds.feature_ids), samples=_samples_table(ds, names),
                     quant=quant_state_of(ds), dataset_hash=_fingerprint(ds),
                     feature_classes=classes, annotated_features=annotated)
