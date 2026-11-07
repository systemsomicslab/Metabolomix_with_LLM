"""`DatasetState` → 規則の入力行、とキュレーションの適用。

行の変換は msdial-workflow の `project_reader`（GUI の表を作る側）に合わせる:
名前が空・None → "Unknown"、ontology・adduct・comment の None → ""、
msms_assigned は "True" / "False"、NaN の高さ → None。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from metabolomix.core.atomic_io import DomainError

__all__ = ["FeatureRow", "dataset_polarity", "feature_rows", "apply_curation"]

_POLARITY_ALIASES = {"positive": "positive", "pos": "positive",
                     "negative": "negative", "neg": "negative"}


@dataclass(frozen=True)
class FeatureRow:
    feature_id: str
    name: str
    ontology: str
    adduct: str
    mz: float
    rt: float
    comment: str
    msms_assigned: str  # "True" / "False"
    heights: tuple[float | None, ...]
    detected: tuple[bool, ...]


def _norm_polarity(value) -> str | None:
    if not isinstance(value, str):
        return None
    return _POLARITY_ALIASES.get(value.strip().lower())


def dataset_polarity(ds, holder, polarity: str | None) -> str:
    """引数 → ハンドルの `polarity` → `ds.ion_mode_majority` の順。

    解釈できない値は次の候補へ落ちる。決まらなければ QUANT_POLARITY_REQUIRED。
    `holder` は `.polarity` を持たない対話のセッションでもよい。
    """
    for candidate in (polarity, getattr(holder, "polarity", None),
                      getattr(ds, "ion_mode_majority", None)):
        found = _norm_polarity(candidate)
        if found:
            return found
    raise DomainError("QUANT_POLARITY_REQUIRED",
                      'データセットの極性を決められません。polarity="positive" か '
                      'polarity="negative" を渡してください。')


def _height(value) -> float | None:
    if value is None:
        return None
    value = float(value)
    return None if math.isnan(value) else value


def feature_rows(ds) -> list[FeatureRow]:
    """特徴ごとに 1 行（`ds.feature_ids` の順）。gap-fill を区別できない入力では定量しない。"""
    if getattr(ds, "detected_mask", None) is None or getattr(ds, "feature_matrix", None) is None:
        raise DomainError("QUANT_NEEDS_DETECTION",
                          "検出の印（detected_mask）と値の行列が無い入力では、gap-fill を区別できないため"
                          "定量しません（MS-DIAL の解析フォルダを読んでください）。")
    rows = []
    for i, fid in enumerate(ds.feature_ids):
        meta = ds.feature_metadata.get(fid, {}) or {}
        name = meta.get("name")
        rows.append(FeatureRow(
            feature_id=str(fid),
            name=str(name).strip() if name and str(name).strip() else "Unknown",
            ontology=meta.get("ontology") or "",
            adduct=meta.get("adduct") or "",
            mz=float(meta.get("mz")),
            rt=float(meta.get("rt")),
            comment=meta.get("comment") or "",
            msms_assigned=str(bool(meta.get("msms_assigned"))),
            heights=tuple(_height(v) for v in ds.feature_matrix[i]),
            detected=tuple(bool(v) for v in ds.detected_mask[i]),
        ))
    return rows


def _spot(feature_id: str):
    try:
        return int(feature_id)
    except (TypeError, ValueError):
        return feature_id


def apply_curation(rows: list[FeatureRow], flags: dict | None) -> tuple[list[FeatureRow], list[dict]]:
    """人の判断を当てる。`wrong` / `redundant` は除外、`assign` は name・ontology・adduct の置換。

    `wrong` と `assign` が同じスポットに当たったら `wrong` が勝つ（除外。置換の audit は書かない）。
    `assign` の判断の行に無い欄は今の値を残す。`suspect` は当てない。
    audit は {feature_id, action: excluded|replaced, reason, (replaced なら from / to)}。
    """
    flags = flags or {}
    wrong = set(flags.get("wrong") or ())
    redundant = set(flags.get("redundant") or ())
    assign = flags.get("assign") or {}
    out, audit = [], []
    for row in rows:
        spot = _spot(row.feature_id)
        if spot in wrong:
            audit.append({"feature_id": row.feature_id, "action": "excluded",
                          "reason": "curation_wrong"})
            continue
        if spot in redundant:
            audit.append({"feature_id": row.feature_id, "action": "excluded",
                          "reason": "curation_redundant"})
            continue
        decision = assign.get(spot)
        if decision is not None:
            new = {k: decision.get(k) or getattr(row, k) for k in ("name", "ontology", "adduct")}
            old = {k: getattr(row, k) for k in new}
            audit.append({"feature_id": row.feature_id, "action": "replaced",
                          "reason": "curation_assign", "from": old, "to": new})
            row = replace(row, **new)
        out.append(row)
    return out, audit
