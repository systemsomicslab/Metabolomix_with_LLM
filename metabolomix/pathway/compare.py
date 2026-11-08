# metabolomix/pathway/compare.py
"""JSON-LD のスコープ値で 2 群比較し、隣接ノードの逆向きを拾う（純ロジック）。

検定は統計エンジン（`analysis.stats.compare`）の moderated t: スコープ種別（whole / 鎖ごと）ごとに
試料 × スコープの log2 行列（正の有限値の log2、それ以外は欠測）を作り、`subject_of` で同じ個体の
注入を log2 の上で平均してから比べる。BH はパスウェイ（＝1 グラフ）× スコープ種別の単位
（エンジンの母集団 `matrix` ＝その種別のスコープ全部）でかかる。判定（verdict）は q ≤ 閾値と CI で決まる。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis.stats.compare import compare as engine_compare
from metabolomix.analysis.stats.external import prepared_from_injections
from metabolomix.pathway.jsonld import PathwayGraph, adjacent_lipid_pairs

_CONTRAST = {"id": "c", "factor": "group"}


def _scope_key(scope) -> str:
    return "whole" if scope.scope_type == "subclass" else "+".join(scope.chains)


def _num(v):
    return None if v is None or not math.isfinite(float(v)) else float(v)


def _log2(v) -> float:
    """スコープの値 → log2。正の有限値だけが値で、それ以外（0・負・欠測）は欠測（NaN）。"""
    if v is None:
        return math.nan
    v = float(v)
    return math.log2(v) if math.isfinite(v) and v > 0 else math.nan


def _call(row: dict) -> str:
    """verdict が changed なら符号で up / down、検定していなければ NA、それ以外は ns。"""
    if row["status"] != "tested":
        return "NA"
    if row["verdict"] == "changed":
        return "up" if row["log2fc"] > 0 else "down"
    return "ns"


def compare_scopes(graph: PathwayGraph, group_a: str, group_b: str, q_threshold: float,
                   *, subject_of: dict[str, str] | None = None) -> list[dict]:
    """スコープ種別ごとに 1 回、エンジンの moderated t で group_a（基準）と group_b を比べる。

    `subject_of`: 試料（JSON-LD の注入）→ 個体。同じ個体の注入は log2 の有限な値の平均で 1 個体にまとまる
    （前処理と同じ規則）。無い・載っていない試料はその注入だけで 1 個体。log2fc は正なら group_b が高い。
    各群の個体数が 2 未満なら `DomainError(STATISTIC_GROUP_TOO_SMALL)`。
    """
    subject_of = subject_of or {}
    samples = [s for s, g in graph.sample_groups.items() if g in (group_a, group_b)]
    labels = [graph.sample_groups[s] for s in samples]
    by_kind: dict[str, list] = {}
    for scope in graph.scopes.values():
        if scope.n_features == 0:
            # 上流はデータに無いサブクラスのノード・データに無い鎖にも空スコープ
            # （memberFeature: []、値は全部 null）を出す。値の付いたノードではない（spec §3）
            continue
        by_kind.setdefault(_scope_key(scope), []).append(scope)

    comparison = {**_CONTRAST, "reference": group_a, "test": group_b}
    out: list[dict] = []
    kinds = ["whole"] + sorted(k for k in by_kind if k != "whole")
    for kind in kinds:
        scopes = sorted(by_kind.get(kind, []), key=lambda s: graph.node_labels.get(s.node_id, ""))
        if not scopes:
            continue
        log2 = np.array([[_log2(graph.values.get(s.id, {}).get(name)) for s in scopes]
                         for name in samples], dtype=float).reshape(len(samples), len(scopes))
        prepared = prepared_from_injections(
            log2, sample_ids=samples, groups=labels,
            subject_ids=[subject_of.get(name) for name in samples],
            feature_ids=[s.id for s in scopes], matrix_id=f"pathway:{kind}")
        result = engine_compare(prepared, factors=["group"], comparisons=[comparison],
                                test="moderated_t", q_threshold=q_threshold)
        for scope, r in zip(scopes, result["rows"]):
            out.append({"node": graph.node_labels.get(scope.node_id, ""),
                        "scope": kind, "node_id": scope.node_id,
                        "n_features": scope.n_features,
                        "mean_a": _num(r["mean_reference"]), "mean_b": _num(r["mean_test"]),
                        "log2fc": _num(r["log2fc"]), "p": _num(r["p"]),
                        "ci_low": _num(r["ci_low"]), "ci_high": _num(r["ci_high"]),
                        "q": _num(r["q"]), "verdict": r["verdict"], "call": _call(r)})
    return out


def opposite_pairs(graph: PathwayGraph, rows: list[dict]) -> tuple[list[dict], int]:
    whole = {r["node_id"]: r for r in rows if r["scope"] == "whole"}
    pairs, undeterminable = [], 0
    for a, b, via in adjacent_lipid_pairs(graph):
        ra, rb = whole.get(a), whole.get(b)
        if ra is None or rb is None:
            continue                     # 値の付かないノード同士は比較の対象外
        calls = {ra["call"], rb["call"]}
        if calls & {"ns", "NA"}:
            undeterminable += 1
            continue
        if ra["call"] != rb["call"]:
            first, second = sorted((ra, rb), key=lambda r: r["node"])
            pairs.append({"node_1": first["node"], "call_1": first["call"],
                          "log2fc_1": first["log2fc"],
                          "node_2": second["node"], "call_2": second["call"],
                          "log2fc_2": second["log2fc"], "via": via})
    pairs.sort(key=lambda p: (p["node_1"], p["node_2"]))
    return pairs, undeterminable


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def rows_to_tsv(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    lines += ["\t".join(_fmt(r.get(c)) for c in columns) for r in rows]
    return "\n".join(lines)
