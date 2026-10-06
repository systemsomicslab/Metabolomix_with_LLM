"""内部標準比（spec §8・§8.1）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

補正値は**同一注入**の `target_height / standard_height`。単位は
`internal_standard_ratio` であり、絶対濃度ではない（校正線も回収率も入っていない）。

## 分母を信用できないセルをどう扱うか

0・非有限・負、または検出条件を満たさない分母で割った値は、比として意味を持たない。
微小定数（pseudocount）を足して割るのは**最悪の選択**で、比は分母の大きさではなく
足した定数の大きさを反映した巨大値になり、しかも数値として通ってしまう。ここでは
欠損にし、その欠損を `locked` として印を付けて**後段の補完も禁止**する
（spec §8「分母不正で生じた欠損は後段補完も禁止し、sample/featureの理由を保存する」）。

分子（対象）側の欠損は locked にしない。これは通常の欠測で、補完してよいかは
recipe の `impute` が決める——分母が無かったことと、対象が検出されなかったことは
別の事実で、同じ扱いにすると前者が後者の顔をして埋まる。

## 内部標準の列そのものは残す

内部標準は解析対象ではなく **support**。比を作った後も元の値のまま列に残す
（比にしてしまうと自分自身で割って常に1になる）。解析対象かどうかは
eligibility mask が持つ別の軸で、列の有無で表現しない——列を消すと、後から
分母を再計算できなくなる。
"""
from __future__ import annotations

import numpy as np

from metabolomix.core.atomic_io import DomainError

__all__ = [
    "RATIO_UNIT",
    "RAW_UNIT",
    "apply_internal_standards",
    "ratio",
    "validate_pairs",
]

RATIO_UNIT = "internal_standard_ratio"
RAW_UNIT = "peak_height"

_MAP_INVALID = "INTERNAL_STANDARD_MAP_INVALID"


def ratio(target, standard, standard_detected):
    """`target / standard` と、補完禁止 mask を返す（元の配列は変更しない）。

    `standard_detected` は分母側の検出 mask（`None` なら検出条件を課さない）。
    gap-fill された分母は補間値であって観測ではないので、値が入っていても
    比の根拠にしない。

    戻り値の `locked` は「分母が不正だったセル」。分子側の欠損で NaN になった
    セルは locked にしない（後段で補完してよいかは別の判断）。
    """
    target = np.asarray(target, dtype=float)
    standard = np.asarray(standard, dtype=float)
    if target.shape != standard.shape:
        raise ValueError(
            f"target と standard の形が違います: {target.shape} != {standard.shape}")

    valid = np.isfinite(standard) & (standard > 0)
    if standard_detected is not None:
        detected = np.asarray(standard_detected, dtype=bool)
        if detected.shape != standard.shape:
            raise ValueError(
                f"standard_detected の形が違います: {detected.shape} "
                f"!= {standard.shape}")
        valid = valid & detected

    out = np.full(target.shape, np.nan, dtype=float)
    np.divide(target, standard, out=out, where=valid & np.isfinite(target))
    return out, ~valid


def validate_pairs(pairs: list[dict], feature_ids: list[str]) -> None:
    """`internal-standard-map.v1` の対応が数値計算に載る形かを確認する。

    一対象に複数標準・自己参照・循環・存在しないIDは `INTERNAL_STANDARD_MAP_INVALID`
    （spec §8）。循環を通すと「補正済みの値で補正する」入れ子ができ、単位が
    何なのか誰にも言えなくなる。
    """
    known = set(feature_ids)
    seen_targets: set[str] = set()
    edges: dict[str, str] = {}

    for pair in pairs:
        target_fid = pair.get("target_feature_id")
        standard_fid = pair.get("standard_feature_id")
        target_id = pair.get("target_id")

        if target_id in seen_targets or target_fid in edges:
            raise DomainError(
                _MAP_INVALID,
                f"1つの対象に複数の内部標準が対応しています: target_id={target_id!r}",
                {"target_id": target_id, "target_feature_id": target_fid})
        seen_targets.add(target_id)

        missing = [fid for fid in (target_fid, standard_fid) if fid not in known]
        if missing:
            raise DomainError(
                _MAP_INVALID,
                f"存在しないfeature_idを参照しています: {', '.join(map(repr, missing))}",
                {"target_id": target_id, "missing_feature_ids": missing})

        if target_fid == standard_fid:
            raise DomainError(
                _MAP_INVALID,
                f"自分自身を内部標準にしています: feature_id={target_fid!r}",
                {"target_id": target_id, "feature_id": target_fid})

        edges[target_fid] = standard_fid

    for start in edges:
        seen = {start}
        node = edges[start]
        while node in edges:
            if node in seen:
                raise DomainError(
                    _MAP_INVALID,
                    f"内部標準の対応が循環しています（feature_id={node!r}）。",
                    {"feature_id": node, "chain": sorted(seen)})
            seen.add(node)
            node = edges[node]


def apply_internal_standards(matrix, feature_ids: list[str], pairs: list[dict],
                             detected_mask) -> dict:
    """(assay × feature) 行列へ内部標準比を適用する（元の行列は変更しない）。

    `detected_mask` は同じ形の bool 行列（`None` なら検出条件を課さない）。

    対応の無い feature は元の単位のまま残し、`units` で区別する——補正済みと
    未補正を「同じ数値」として1つの統計へ流さないための印で、列を分けるのでは
    なく単位で持つ（列を分けると feature 軸が recipe ごとに変わる）。
    """
    values = np.array(matrix, dtype=float, copy=True)
    validate_pairs(pairs, feature_ids)

    index_of = {fid: i for i, fid in enumerate(feature_ids)}
    units = [RAW_UNIT] * len(feature_ids)
    locked = np.zeros(values.shape, dtype=bool)
    missing_reasons: dict[str, dict] = {}
    support: list[str] = []

    source = np.asarray(matrix, dtype=float)
    for pair in pairs:
        target_index = index_of[pair["target_feature_id"]]
        standard_index = index_of[pair["standard_feature_id"]]
        standard_detected = (None if detected_mask is None
                             else np.asarray(detected_mask, dtype=bool)[:, standard_index])
        column, column_locked = ratio(source[:, target_index],
                                      source[:, standard_index],
                                      standard_detected)
        values[:, target_index] = column
        locked[:, target_index] = column_locked
        units[target_index] = RATIO_UNIT
        if pair["standard_feature_id"] not in support:
            support.append(pair["standard_feature_id"])
        missing_reasons[pair["target_feature_id"]] = {
            str(row): {"code": "standard_denominator_invalid",
                       "standard_feature_id": pair["standard_feature_id"],
                       "target_id": pair.get("target_id"),
                       "standard_target_id": pair.get("standard_target_id")}
            for row in np.flatnonzero(column_locked)
        }

    return {
        "values": values,
        "units": units,
        "locked_mask": locked,
        "missing_reasons": missing_reasons,
        "support_feature_ids": support,
    }
