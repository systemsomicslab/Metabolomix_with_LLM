# metabolomix/pathway/compare.py
"""JSON-LD のスコープ値で 2 群比較し、隣接ノードの逆向きを拾う（純ロジック）。

検定は既存の v1 経路と同じ: log2(x+1) 空間の Welch（`differential.two_group_test`）。
BH はパスウェイ（＝1 グラフ）× スコープ種別（whole / 鎖ごと）の単位でかける。
"""
from __future__ import annotations

import math

import numpy as np

from metabolomix.analysis import differential
from metabolomix.pathway.jsonld import PathwayGraph, adjacent_lipid_pairs


def _scope_key(scope) -> str:
    return "whole" if scope.scope_type == "subclass" else "+".join(scope.chains)


def _num(v):
    return None if v is None or not math.isfinite(float(v)) else float(v)


def compare_scopes(graph: PathwayGraph, group_a: str, group_b: str,
                   q_threshold: float) -> list[dict]:
    samples = [s for s, g in graph.sample_groups.items() if g in (group_a, group_b)]
    labels = [graph.sample_groups[s] for s in samples]
    by_kind: dict[str, list] = {}
    for scope in graph.scopes.values():
        if scope.n_features == 0:
            # 上流はデータに無いサブクラスのノード・データに無い鎖にも空スコープ
            # （memberFeature: []、値は全部 null）を出す。値の付いたノードではない（spec §3）
            continue
        by_kind.setdefault(_scope_key(scope), []).append(scope)

    out: list[dict] = []
    kinds = ["whole"] + sorted(k for k in by_kind if k != "whole")
    for kind in kinds:
        scopes = sorted(by_kind.get(kind, []), key=lambda s: graph.node_labels.get(s.node_id, ""))
        if not scopes:
            continue
        matrix = np.array([[np.nan if graph.values.get(s.id, {}).get(name) is None
                            else graph.values[s.id][name] for s in scopes]
                           for name in samples], dtype=float)
        results = differential.two_group_test(matrix, [s.id for s in scopes], labels,
                                              group_a, group_b, log_transform=True)
        qs = differential.bh_fdr([r["p"] for r in results])
        for scope, r, q in zip(scopes, results, qs):
            q = _num(q)
            fc = _num(r["log2fc"])
            if q is None or fc is None:
                call = "NA"
            elif q < q_threshold:
                call = "up" if fc > 0 else "down"
            else:
                call = "ns"
            out.append({"node": graph.node_labels.get(scope.node_id, ""),
                        "scope": kind, "node_id": scope.node_id,
                        "n_features": scope.n_features,
                        "mean_a": _num(r["mean_a"]), "mean_b": _num(r["mean_b"]),
                        "log2fc": fc, "q": q, "call": call})
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
