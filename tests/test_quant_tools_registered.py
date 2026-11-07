"""登録された定量ツール 4 本を、server 経由（handle_dispatch の包みを通して）呼ぶ。"""
from __future__ import annotations

import asyncio
import json

import pytest

import server
from metabolomix.core import handles
from tests.test_quant_normalize import REF, _ambiguous, _point, _simple

TOOLS = ("dataset_quant_preview", "dataset_quant_apply", "dataset_merge_polarities",
         "dataset_select_features")


def _call(name, **args):
    out = asyncio.run(server.mcp.call_tool(name, args))
    content = out[0] if isinstance(out, tuple) else out
    return json.loads(content[0].text)


def _open(ds, polarity="positive"):
    handle = handles.DatasetHandle(ds, polarity=polarity, polarity_source="test",
                                   source={"kind": "test"})
    handles.registry.add(handle)
    return handle.handle_id


def test_four_tools_registered_with_annotations():
    tm = server.mcp._tool_manager._tools
    for name in TOOLS:
        assert name in tm and tm[name].annotations is not None
        assert tm[name].output_schema is None
    for name in ("dataset_quant_preview", "dataset_quant_apply", "dataset_select_features"):
        assert tm[name].is_async  # dataset_handle を持つので handle_dispatch が包む
    assert "dataset_handle" not in tm["dataset_merge_polarities"].parameters["properties"]


def test_preview_and_apply_via_dataset_handle(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ref = tmp_path / "ref.csv"
    ref.write_text(REF, encoding="utf-8")
    hid = _open(_simple(tmp_path))
    prev = _call("dataset_quant_preview", reference_path=str(ref), dataset_handle=hid)
    assert prev["status"] == "ready"
    assert str(tmp_path) not in json.dumps(prev)
    applied = _call("dataset_quant_apply", reference_path=str(ref), dataset_handle=hid)
    assert applied["dataset_handle"] != hid
    assert handles.registry.get(applied["dataset_handle"]).dataset is not None


def test_preview_needs_input_envelope_for_ambiguous(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ref = tmp_path / "ref.csv"
    ref.write_text(REF, encoding="utf-8")
    hid = _open(_ambiguous(tmp_path))
    out = _call("dataset_quant_apply", reference_path=str(ref), dataset_handle=hid)
    assert out["status"] == "needs_input"


@pytest.mark.parametrize("name,args", [
    ("dataset_quant_preview", {"reference_path": "x.csv"}),
    ("dataset_quant_apply", {"skip": True}),
    ("dataset_select_features", {"skip_adduct_filter": True}),
])
def test_no_dataset_is_missing_state(name, args):
    out = _call(name, **args)
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["state"] == "dataset"


def test_unknown_handle_is_envelope_not_exception():
    out = _call("dataset_quant_apply", skip=True, dataset_handle="h-nope")
    assert out["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"
    out = _call("dataset_merge_polarities", positive_handle="h-nope", negative_handle="h-no2")
    assert out["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"


def test_rules_unavailable_envelopes_do_not_leak_paths(tmp_path, monkeypatch):
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    ref = tmp_path / "ref.csv"
    ref.write_text(REF, encoding="utf-8")
    hid = _open(_simple(tmp_path))
    hid2 = _open(_simple(tmp_path / "second"), polarity="negative")
    calls = [
        ("dataset_quant_preview", {"reference_path": str(ref), "dataset_handle": hid}),
        ("dataset_quant_apply", {"reference_path": str(ref), "dataset_handle": hid}),
        ("dataset_merge_polarities", {"positive_handle": hid, "negative_handle": hid2}),
        ("dataset_select_features", {"skip_adduct_filter": True, "dataset_handle": hid}),
    ]
    for name, args in calls:
        out = _call(name, **args)
        assert out["error"]["code"] == "quant_rules_unavailable", (name, out)
        assert str(tmp_path) not in json.dumps(out), name
