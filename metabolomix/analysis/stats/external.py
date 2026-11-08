"""外から渡された注入ごとの log2 値から `PreparedMatrix` を作る（前処理は通さない）。

前処理済みの行列を持たない利用者（例: MSLipidMapper が返すスコープの値を比べる
`pathway/compare.py`）がエンジンの `compare` を使うための入口。注入を個体にまとめる規則は
前処理の `_collapse` と同じ: 個体（subject）× 群が同じ注入を log2 の上で、有限な値だけの平均に
する（全注入が欠測なら欠測）。個体の指定が無い注入は、それぞれが 1 個体。
MCP・セッションには依存しない。
"""
from __future__ import annotations

import warnings

import numpy as np

from metabolomix.analysis.stats.preprocess import PreparedMatrix, UnitInfo

__all__ = ["prepared_from_injections"]


def prepared_from_injections(log2, *, sample_ids: list[str], groups: list[str],
                             subject_ids: list[str | None], feature_ids: list[str],
                             matrix_id: str = "external") -> PreparedMatrix:
    """`log2`（注入 × 特徴。NaN = 欠測）を個体 × 特徴の `PreparedMatrix` にまとめる。

    `groups` は注入ごとの群、`subject_ids` は注入ごとの個体 ID（None はその注入だけで 1 個体）。
    同じ個体でも群が違えば別の個体にする（前処理と同じ）。個体の並びは注入の初出の順。
    """
    log2 = np.asarray(log2, dtype=float)
    n = len(sample_ids)
    if log2.ndim != 2 or log2.shape != (n, len(feature_ids)) or len(groups) != n or len(subject_ids) != n:
        raise ValueError("log2 の形（注入 × 特徴）と sample_ids / groups / subject_ids / feature_ids の長さが合わない")
    members: dict[tuple, list[int]] = {}
    for k in range(n):
        key = (("s", subject_ids[k], groups[k]) if subject_ids[k] is not None else ("i", sample_ids[k]))
        members.setdefault(key, []).append(k)
    units, rows, taken = [], [], set()
    for key, idx in members.items():
        base = sample_ids[idx[0]] if key[0] == "i" else f"{key[1]}@{key[2]}"
        uid, num = base, 2
        while uid in taken:
            uid, num = f"{base}#{num}", num + 1
        taken.add(uid)
        units.append(UnitInfo(unit_id=uid, group=groups[idx[0]], batch=None,
                              subject_id=None if key[0] == "i" else key[1],
                              injection_order=None, n_injections=len(idx),
                              sample_ids=[sample_ids[k] for k in idx], attributes={}))
        if len(idx) > 1:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                rows.append(np.nanmean(log2[idx], axis=0))
        else:
            rows.append(log2[idx[0]])
    with_subject = sum(s is not None for s in subject_ids)
    source = ("subject_id" if with_subject == n and n else "injection" if with_subject == 0 else "mixed")
    matrix = np.vstack(rows) if rows else np.empty((0, len(feature_ids)))
    return PreparedMatrix(
        matrix_id=matrix_id, log2=matrix, units=units, feature_ids=[str(f) for f in feature_ids],
        detected_units=None, qc_log2=np.empty((0, len(feature_ids))), qc_ids=[], steps=[],
        config={}, quant=None, dataset_hash=None, subject_source=source, warnings=[],
        feature_classes={}, annotated_features=frozenset(str(f) for f in feature_ids))
