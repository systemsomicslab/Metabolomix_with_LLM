"""`DatasetState` → 規則の入力行、とキュレーションの適用。

行の変換は msdial-workflow の `project_reader`（GUI の表を作る側）に合わせる:
名前が空・None → "Unknown"、ontology・adduct・comment の None → ""、
msms_assigned は "True" / "False"、NaN の高さ → None。
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, replace

from metabolomix.core.atomic_io import DomainError
from metabolomix.curation.apply import identity_for

__all__ = ["FeatureRow", "resolve_quant_polarity", "require_msdial_source", "feature_rows",
           "apply_curation", "has_curation_decisions", "curation_digest"]

#: 定量を掛けられる入力の形式（`DatasetState.source_format`）。
MSDIAL_SOURCE = "msdial_folder"

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


def resolve_quant_polarity(ds, holder, polarity: str | None) -> str:
    """引数 → ハンドルの `polarity` → データセット由来（`mztab.polarity.dataset_polarity`）の順。

    明示された値（引数・ハンドル）が解釈できなければ QUANT_POLARITY_INVALID（打ち間違いで
    別の極性を黙って選ばない）。データセット由来が決まらなければ QUANT_POLARITY_REQUIRED。
    `holder` は `.polarity` を持たない対話のセッションでもよい。
    """
    for candidate in (polarity, getattr(holder, "polarity", None)):
        if candidate is None or (isinstance(candidate, str) and not candidate.strip()):
            continue
        found = _norm_polarity(candidate)
        if found is None:
            raise DomainError("QUANT_POLARITY_INVALID",
                              'polarity は "positive" / "negative"（"pos" / "neg" も可）のいずれかです。',
                              {"polarity": candidate})
        return found
    from metabolomix.mztab.polarity import dataset_polarity
    found = _norm_polarity(dataset_polarity(ds)[0])
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


def require_msdial_source(ds) -> None:
    """MS-DIAL の解析フォルダから読んだデータセットでなければ QUANT_SOURCE_UNSUPPORTED。

    規則の入力（クラス・付加体・コメント・MS/MS の有無）は解析フォルダのアダプタが
    `feature_metadata` に置く欄だけを読む。ほかの形式はこれらを同じ意味で持たないため、
    黙って空欄のまま定量しない。
    """
    source = getattr(ds, "source_format", None)
    if source != MSDIAL_SOURCE:
        raise DomainError("QUANT_SOURCE_UNSUPPORTED",
                          "定量は MS-DIAL の解析フォルダ（.arf / .arf2）から読んだデータセットにだけ"
                          "掛けます（クラス・付加体・コメント・MS/MS の有無を同じ意味で読めるのは"
                          "解析フォルダだけ）。解析フォルダを dataset_load / dataset_handle_open で"
                          "開いてください。", {"source_format": source})


def feature_rows(ds) -> list[FeatureRow]:
    """特徴ごとに 1 行（`ds.feature_ids` の順）。解析フォルダ以外の入力と、gap-fill を区別
    できない入力では定量しない。"""
    require_msdial_source(ds)
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
    """人の判断を当てる。`wrong` / `redundant` は除外、`assign` は name・ontology の置換。

    `assign` は今の書き出し（`curation.apply.identity_for`）と同じ意味: 判断に無い欄は空に
    なる。置換後の名前が空なら MS-DIAL の空名と同じく "Unknown"。adduct は変えない。
    `wrong` と `assign` が同じスポットに当たったら `wrong` が勝つ（除外。置換の audit は書かない）。
    `suspect` は当てない。
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
        if spot in assign:
            ident = identity_for(spot, flags)
            new = {"name": ident["name"].strip() or "Unknown", "ontology": ident["ontology"]}
            old = {k: getattr(row, k) for k in new}
            audit.append({"feature_id": row.feature_id, "action": "replaced",
                          "reason": "curation_assign", "from": old, "to": new})
            row = replace(row, **new)
        out.append(row)
    return out, audit


def has_curation_decisions(flags: dict | None) -> bool:
    """当てる判断（wrong / redundant / assign）が 1 件でもあるか。"""
    flags = flags or {}
    return bool(flags.get("wrong") or flags.get("redundant") or flags.get("assign"))


def curation_digest(flags: dict | None) -> str:
    """`apply_curation` が当てる判断だけの安定した sha256（suspect などは含めない）。

    wrong / redundant のスポットと、assign のスポットと置換後の name・ontology。
    判断が無ければ空の集合の digest（後から判断が付けば変わる）。
    """
    flags = flags or {}
    canonical = {
        "wrong": sorted(str(s) for s in flags.get("wrong") or ()),
        "redundant": sorted(str(s) for s in flags.get("redundant") or ()),
        "assign": sorted([str(s), (identity_for(s, flags) or {}).get("name", ""),
                          (identity_for(s, flags) or {}).get("ontology", "")]
                         for s in (flags.get("assign") or {})),
    }
    text = json.dumps(canonical, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
