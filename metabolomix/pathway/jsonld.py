# metabolomix/pathway/jsonld.py
"""MSLipidMapper の pathway measurements JSON-LD（schema 2.0.0）を読む（純ロジック）。

上流文書: docs/pathway-measurements.md（MSLipidMapper リポジトリ）。型は `mslm:` 接頭辞付きの
compact 形で来る。参照は {"@id": ...}。値の欠損は JSON null（MissingValue）。
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA_VERSION = "2.0.0"
_KIND = {"mslm:PathwayLipidNode": "lipid", "mslm:PathwayEntity": "entity",
         "mslm:PathwayJunction": "junction"}


class JsonldError(Exception):
    code = "MSLM_OUTPUT_INVALID"

    def __init__(self, message: str):
        super().__init__(message)
        self.message = message


@dataclass
class Scope:
    id: str
    node_id: str
    scope_type: str
    chains: tuple[str, ...]
    n_features: int


@dataclass
class PathwayGraph:
    node_labels: dict[str, str] = field(default_factory=dict)
    node_kinds: dict[str, str] = field(default_factory=dict)
    edges: list[tuple[str, str]] = field(default_factory=list)
    scopes: dict[str, Scope] = field(default_factory=dict)
    sample_groups: dict[str, str] = field(default_factory=dict)
    values: dict[str, dict[str, float | None]] = field(default_factory=dict)
    unmatched_subclasses: list[str] = field(default_factory=list)


def _ref(obj, key) -> str | None:
    value = obj.get(key)
    return value.get("@id") if isinstance(value, dict) else None


def load_graph(path: Path) -> PathwayGraph:
    try:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        items = doc["@graph"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise JsonldError(f"JSON-LD を読めません: {path} ({exc})") from exc
    datasets = [x for x in items if x.get("@type") == "mslm:Dataset"]
    version = datasets[0].get("schemaVersion") if datasets else None
    if version != SCHEMA_VERSION:
        raise JsonldError(f"schemaVersion が {version!r} です（想定 {SCHEMA_VERSION}）。"
                          "上流 MSLipidMapper の出力仕様が変わった可能性があります。")

    g = PathwayGraph()
    sample_uri_to_id: dict[str, str] = {}
    subclass_labels: dict[str, str] = {}
    represented: set[str] = set()
    measurements = []
    for x in items:
        t = x.get("@type")
        if t in _KIND:
            g.node_kinds[x["@id"]] = _KIND[t]
            g.node_labels[x["@id"]] = str(x.get("label") or "")
            rep = _ref(x, "representsEntity")
            if rep:
                represented.add(rep)
        elif t == "mslm:PathwayConnection":
            src, tgt = _ref(x, "sourceEntity"), _ref(x, "targetEntity")
            if src and tgt:
                g.edges.append((src, tgt))
        elif t == "mslm:LipidMeasurementScope":
            g.scopes[x["@id"]] = Scope(
                id=x["@id"], node_id=_ref(x, "onPathwayNode") or "",
                scope_type=str(x.get("scopeType") or ""),
                chains=tuple(x.get("acylChainSignature") or ()),
                n_features=len(x.get("memberFeature") or ()))
        elif t == "mslm:Sample":
            sample_uri_to_id[x["@id"]] = str(x.get("identifier"))
            g.sample_groups[str(x.get("identifier"))] = str(x.get("sampleGroup") or "")
        elif t == "mslm:MSDialMappableLipidEntity":
            subclass_labels[x["@id"]] = str(x.get("label") or "")
        elif t == "mslm:SampleMeasurement":
            measurements.append(x)

    values: dict[str, dict[str, float | None]] = defaultdict(dict)
    for x in measurements:
        entity = _ref(x, "measuredEntity")
        if entity not in g.scopes:
            continue                       # 分子種（feature）の値はスコープ値に混ぜない
        sample = sample_uri_to_id.get(_ref(x, "inSample") or "")
        if sample is None:
            continue
        v = x.get("value")
        values[entity][sample] = None if v is None else float(v)
    g.values = dict(values)
    g.unmatched_subclasses = sorted(label for uri, label in subclass_labels.items()
                                    if uri not in represented)
    return g


def adjacent_lipid_pairs(graph: PathwayGraph) -> list[tuple[str, str, str]]:
    """接合ノードだけを挟んでつながる脂質ノードの組（向きは問わない）。"""
    neighbors: dict[str, set[str]] = defaultdict(set)
    direct: set[frozenset] = set()
    for a, b in graph.edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
        if graph.node_kinds.get(a) == "lipid" and graph.node_kinds.get(b) == "lipid":
            direct.add(frozenset((a, b)))
    via_junction: set[frozenset] = set()
    for start, kind in graph.node_kinds.items():
        if kind != "lipid":
            continue
        seen, stack = {start}, [n for n in neighbors[start]
                                if graph.node_kinds.get(n) == "junction"]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            for nxt in neighbors[node]:
                k = graph.node_kinds.get(nxt)
                if k == "junction" and nxt not in seen:
                    stack.append(nxt)
                elif k == "lipid" and nxt != start:
                    via_junction.add(frozenset((start, nxt)))
    out = []
    for pair in direct | via_junction:
        a, b = sorted(pair)
        out.append((a, b, "direct" if pair in direct else "junction"))
    return sorted(out)
