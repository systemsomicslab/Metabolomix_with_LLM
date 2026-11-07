"""`.arf2` のスポットから、同定の写しに要る欄を読む。

Key 番号の正準は docs/schema/AlignmentSpotProperty.md と MsScanMatchResult.md。
マッチ結果のコンテナは Key0 結果の一覧・Key1 MSP の辞書・Key2 Text DB の一覧。
deps: arf2.reader / arf2.match_results のみ。
"""
from __future__ import annotations

from metabolomix.arf2.match_results import (
    _candidates, decode_match_result, representative, spot_master_id)
from metabolomix.arf2.reader import IonMode, _decode, _extract_times, _nested_at, load_raw_spots

SOURCE_TEXT_DB = 16  # MsScanMatchResult Key26 の TextDB ビット


def _get(raw: list, index: int):
    return raw[index] if len(raw) > index else None


def _text(value) -> str | None:
    return _decode(value) or None


def _float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_text_db_representative(container, rep: dict | None) -> bool:
    """上流 `IsTextDbBasedRepresentative` の移植。`Source & 16`、または
    Key2 の一覧に代表と同じ復号結果がある。"""
    if rep is None:
        return False
    if int(rep.get("source") or 0) & SOURCE_TEXT_DB:
        return True
    if isinstance(container, list) and len(container) > 2 and isinstance(container[2], list):
        # NaN は復号時に None へ畳まれるので、dict の等値比較で足りる
        return any(decode_match_result(raw) == rep for raw in container[2])
    return False


def spot_fields(raw_spot: list) -> dict:
    master_id = spot_master_id(raw_spot)
    chromx = _get(raw_spot, 4)
    times = _extract_times(chromx) if isinstance(chromx, list) else {}
    container = _get(raw_spot, 56)
    rep = representative(container)
    comment = _get(raw_spot, 25)
    istd = _get(raw_spot, 7)
    return {
        "master_id": master_id,
        "name": _text(_get(raw_spot, 12)),
        "formula": _text(_nested_at(_get(raw_spot, 13), 0)),
        "ontology": _text(_get(raw_spot, 14)),
        "smiles": _text(_get(raw_spot, 15)),
        "inchikey": _text(_get(raw_spot, 16)),
        "adduct": _text(_nested_at(_get(raw_spot, 54), 2)),
        "mz": _float(_get(raw_spot, 5)),
        "rt": times.get("rt"),
        "ion_mode": IonMode.from_int(_get(raw_spot, 11)),
        "comment": _text(comment),
        "internal_standard_id": int(istd) if isinstance(istd, (int, float)) and not isinstance(istd, bool) else None,
        "blank_filtered": bool(_get(raw_spot, 58)),
        "fill_percentage": _float(_get(raw_spot, 49)),
        "representative": rep,
        "is_text_db_representative": _is_text_db_representative(container, rep),
        "is_manually_modified": bool(rep and rep.get("is_manually_modified")),
        "n_candidates": len(_candidates(container)),
    }


def load_spot_fields(arf2_path) -> list[dict]:
    """`.arf2` の全スポットの欄をファイル順に返す。"""
    return [spot_fields(raw) for raw in load_raw_spots(arf2_path)]
