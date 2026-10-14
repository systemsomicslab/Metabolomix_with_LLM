# tests/test_pathway_jsonld.py
import json

import pytest

from metabolomix.pathway import jsonld

B = "https://x/d/"
V = "https://systemomicslab.github.io/MSLipidMapper/ontology/"


def _ref(i):
    return {"@id": B + i}


def _graph(extra=(), schema="2.0.0"):
    g = [
        {"@id": B + "ds", "@type": "mslm:Dataset", "schemaVersion": schema,
         "identifier": "t", "pathway": _ref("pw")},
        {"@id": B + "pw", "@type": "mslm:Pathway"},
        {"@id": B + "sub/PC", "@type": "mslm:MSDialMappableLipidEntity", "label": "PC"},
        {"@id": B + "sub/PE", "@type": "mslm:MSDialMappableLipidEntity", "label": "PE"},
        {"@id": B + "sub/LPC", "@type": "mslm:MSDialMappableLipidEntity", "label": "LPC"},
        {"@id": B + "sub/TG", "@type": "mslm:MSDialMappableLipidEntity", "label": "TG"},
        {"@id": B + "n/PC", "@type": "mslm:PathwayLipidNode", "label": "PC", "representsEntity": _ref("sub/PC")},
        {"@id": B + "n/PE", "@type": "mslm:PathwayLipidNode", "label": "PE", "representsEntity": _ref("sub/PE")},
        {"@id": B + "n/LPC", "@type": "mslm:PathwayLipidNode", "label": "LPC", "representsEntity": _ref("sub/LPC")},
        {"@id": B + "n/PA", "@type": "mslm:PathwayEntity", "label": "PA"},
        {"@id": B + "n/DG", "@type": "mslm:PathwayLipidNode", "label": "DG"},
        {"@id": B + "n/j1", "@type": "mslm:PathwayJunction"},
        # PC -> j1 -> PE, j1 -> LPC（分岐）／ PE -> PA -> DG（PA で隣接が切れる）／ PC -> LPC 直結
        {"@id": B + "c1", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PC"), "targetEntity": _ref("n/j1")},
        {"@id": B + "c2", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/j1"), "targetEntity": _ref("n/PE")},
        {"@id": B + "c3", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/j1"), "targetEntity": _ref("n/LPC")},
        {"@id": B + "c4", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PE"), "targetEntity": _ref("n/PA")},
        {"@id": B + "c5", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PA"), "targetEntity": _ref("n/DG")},
        {"@id": B + "c6", "@type": "mslm:PathwayConnection", "sourceEntity": _ref("n/PC"), "targetEntity": _ref("n/LPC")},
        {"@id": B + "s/1", "@type": "mslm:Sample", "identifier": "S1", "sampleGroup": "C"},
        {"@id": B + "s/2", "@type": "mslm:Sample", "identifier": "S2", "sampleGroup": "P"},
        {"@id": B + "n/PC/scope/whole", "@type": "mslm:LipidMeasurementScope", "onPathwayNode": _ref("n/PC"),
         "scopeType": "subclass", "acylChainSignature": [], "memberFeature": [_ref("f/1"), _ref("f/2")]},
        {"@id": B + "n/PC/scope/acyl-18%3A0", "@type": "mslm:LipidMeasurementScope", "onPathwayNode": _ref("n/PC"),
         "scopeType": "acyl_subset", "acylChainSignature": ["18:0"], "memberFeature": [_ref("f/2")]},
        {"@id": B + "m1", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("n/PC/scope/whole"),
         "inSample": _ref("s/1"), "value": 15, "valueStatus": {"@id": V + "AvailableValue"}},
        {"@id": B + "m2", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("n/PC/scope/whole"),
         "inSample": _ref("s/2"), "value": None, "valueStatus": {"@id": V + "MissingValue"}},
        {"@id": B + "m3", "@type": "mslm:SampleMeasurement", "measuredEntity": _ref("f/1"),
         "inSample": _ref("s/1"), "value": 10},
        *extra,
    ]
    return {"@context": {}, "@graph": g}


def _write(tmp_path, doc):
    p = tmp_path / "remodeling-measurements.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return p


def test_nodes_scopes_samples_and_values(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    assert g.node_kinds[B + "n/PC"] == "lipid"
    assert g.node_kinds[B + "n/PA"] == "entity"
    assert g.node_kinds[B + "n/j1"] == "junction"
    whole = g.scopes[B + "n/PC/scope/whole"]
    assert (whole.scope_type, whole.chains, whole.n_features) == ("subclass", (), 2)
    assert g.scopes[B + "n/PC/scope/acyl-18%3A0"].chains == ("18:0",)
    assert g.sample_groups == {"S1": "C", "S2": "P"}
    assert g.values[B + "n/PC/scope/whole"] == {"S1": 15.0, "S2": None}
    # 分子種（feature）の測定値はスコープ値に混ぜない
    assert B + "f/1" not in g.values


def test_unmatched_subclasses(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    assert g.unmatched_subclasses == ["TG"]


def test_adjacency_through_junctions_and_not_through_entities(tmp_path):
    g = jsonld.load_graph(_write(tmp_path, _graph()))
    pairs = {(g.node_labels[a], g.node_labels[b], via)
             for a, b, via in jsonld.adjacent_lipid_pairs(g)}
    assert ("LPC", "PC", "direct") in pairs or ("PC", "LPC", "direct") in pairs
    labels = {frozenset((a, b)) for a, b, _ in pairs}
    assert frozenset(("PC", "PE")) in labels
    assert frozenset(("LPC", "PE")) in labels      # 分岐した接合ノードの枝同士
    assert frozenset(("PE", "DG")) not in labels   # PA（entity）で切れる
    # PC–LPC は direct と junction の両方で結ばれているが 1 回だけ
    assert sum(1 for a, b, _ in pairs if {a, b} == {"PC", "LPC"}) == 1


def test_wrong_schema_version_is_rejected(tmp_path):
    with pytest.raises(jsonld.JsonldError) as exc:
        jsonld.load_graph(_write(tmp_path, _graph(schema="3.0.0")))
    assert exc.value.code == "MSLM_OUTPUT_INVALID"


def test_broken_json_is_rejected(tmp_path):
    p = tmp_path / "x.json"
    p.write_text("{", encoding="utf-8")
    with pytest.raises(jsonld.JsonldError):
        jsonld.load_graph(p)
