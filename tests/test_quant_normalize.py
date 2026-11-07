"""内部標準の確認と正規化（preview / apply / skip）。人工の解析フォルダと偽の規則で。"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys

import numpy as np
import pytest

from metabolomix.core import handles, session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant import normalize as qn
from metabolomix.quant.inputs import reference_input
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.quant_fake_rules import write_fake_rules

REF = ("TargetClass,StandardName,Coefficient,Unit\nPC,IS-A,2.0,nmol\n"
       "PE,IS-B,3.0,nmol\nAny others,IS-A,1.0,nmol\n")
SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(3)]


def _point(tmp_path, monkeypatch):
    folder = write_fake_rules(tmp_path / "rules")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nrules = '{folder}'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))


def _ref(tmp_path, text=REF, **kw):
    path = tmp_path / "ref_table.csv"
    path.write_text(text, encoding="utf-8")
    return reference_input(str(path), **kw)


def _spot(name, ontology, heights, gap=None, mz=500.0):
    return {"name": name, "ontology": ontology, "heights": heights, "gap_filled": gap, "mz": mz}


def _load(tmp_path, spots):
    folder = tmp_path / "analysis"
    write_msdial_folder(folder, spots=spots, samples=SAMPLES)
    return load_msdial_folder(select_alignment(folder))


def _simple(tmp_path):
    return _load(tmp_path, [
        _spot("IS-A", "PC", [500.0, 500.0, 0.0], gap=[False, True, False]),
        _spot("PC 34:1", "PC", [1000.0, 1000.0, 1000.0]),
        _spot("TG 50:1", "TG", [100.0, 100.0, 100.0], mz=800.0)])


def _ambiguous(tmp_path):
    return _load(tmp_path, [
        _spot("IS-A", "PC", [500.0] * 3),
        _spot("IS-B", "PE", [400.0] * 3, mz=600.0),
        _spot("IS-B", "PE", [200.0] * 3, mz=610.0),
        _spot("PE 36:2", "PE", [1200.0] * 3, mz=700.0)])


def test_preview_needs_input_for_ambiguous_standard(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _ambiguous(tmp_path)
    ref = _ref(tmp_path)
    out = qn.preview(ds, ref, None, "positive")
    assert out["status"] == "needs_input" and out["polarity"] == "positive"
    assert out["unresolved"][0]["standard_name"] == "IS-B"
    assert out["unresolved"][0]["status"] == "standard_ambiguous"
    assert out["unresolved"][0]["target_classes"] == ["PE"]
    assert len(out["unresolved"][0]["candidates"]) == 2
    assert {c["feature_id"] for c in out["unresolved"][0]["candidates"]} == {"1", "2"}
    assert out["reference"]["label"] == "ref_table.csv"
    text = qn.preview_tool(ds, ref, None, "positive")
    assert str(tmp_path) not in text and "ref_table.csv" in text
    ready = qn.preview(ds, ref, {"PE": 2}, "positive")
    assert ready["status"] == "ready" and ready["unresolved"] == []
    pe = next(c for c in ready["classes"] if c["target_class"] == "PE")
    assert pe["feature_id"] == "2" and pe["coefficient"] == 3.0 and pe["via_any_others"] is False


def test_apply_with_bindings_values(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    new, info = qn.apply(ds, _ref(tmp_path), None, "positive")
    row = ds.feature_ids.index("1")
    assert new.feature_matrix[row, 0] == pytest.approx(4.0)
    assert np.isnan(new.feature_matrix[row, 2])
    assert info["status_counts"]["standard_invalid_height"] >= 1
    assert new.quant_state["normalized"] is True and new.quant_state["method"] == "internal_standard"
    assert new.quant_state["unit"] == "nmol" and info["unit"] == "nmol"
    assert new.quant_state["internal_standard_ids"] == ["0"]
    assert ds.quant_state is None and ds.feature_matrix[row, 0] == 1000.0  # 元は不変


def test_apply_with_ambiguous_binding(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _ambiguous(tmp_path)
    new, info = qn.apply(ds, _ref(tmp_path), {"PE": 2}, "positive")
    row = ds.feature_ids.index("3")
    assert new.feature_matrix[row, 0] == pytest.approx(1200 / 200 * 3.0)
    none, preview = qn.apply(ds, _ref(tmp_path), None, "positive")
    assert none is None and preview["status"] == "needs_input"


def test_gapfill_mark_carried(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    new, _ = qn.apply(ds, _ref(tmp_path), None, "positive")
    row = ds.feature_ids.index("1")
    assert new.feature_matrix[row, 1] == pytest.approx(4.0)  # 内部標準が gap-fill でも計算する
    assert bool(new.detected_mask[row, 0]) is True
    assert bool(new.detected_mask[row, 1]) is False  # 印を引き継ぐ


def test_any_others_row_used(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    new, _ = qn.apply(ds, _ref(tmp_path), None, "positive")
    row = ds.feature_ids.index("2")
    assert new.feature_matrix[row, 0] == pytest.approx(100 / 500 * 1.0)
    pv = qn.preview(ds, _ref(tmp_path), None, "positive")
    tg = next(c for c in pv["classes"] if c["target_class"] == "TG")
    assert tg["via_any_others"] is True


def test_apply_refuses_already_normalized(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    new, _ = qn.apply(ds, _ref(tmp_path), None, "positive")
    for fn in (qn.apply, qn.preview):
        with pytest.raises(DomainError) as exc:
            fn(new, _ref(tmp_path), None, "positive")
        assert exc.value.code == "QUANT_ALREADY_NORMALIZED"


def test_apply_writes_record_and_new_handle(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    session_state.session.dataset = ds
    parent = session_state.session.current_handle_id
    out = json.loads(qn.apply_tool(ds, _ref(tmp_path), None, "positive",
                                   holder=session_state.session, parent_handle_id=parent))
    assert out["status"] == "applied" and out["dataset_handle"] != parent
    assert session_state.session.current_handle_id == out["dataset_handle"]
    assert str(tmp_path) not in json.dumps(out)
    rec = json.loads((tmp_path / "analysis" / "reports" / "quant"
                      / f"{out['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["step"] == "normalize" and rec["skipped"] is False
    assert any(p.endswith(".arf") for p in rec["inputs"])
    assert any(p.endswith(".arf2") for p in rec["inputs"])
    assert rec["settings"]["reference"]["sha256"]
    assert rec["settings"]["rules"]["digest"]
    new_ds = handles.registry.get(out["dataset_handle"]).dataset
    assert new_ds.quant_state["record_id"] == out["record_id"]


def test_concentration_form_saves_table(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ref = _ref(tmp_path, concentrations={"IS-A": 10.0, "IS-B": 20.0},
               solvent_volume=100.0, cells_millions=2.0, unit="pmol / 10^6 cells")
    new, info = qn.apply(ds, ref, None, "positive")
    csv = tmp_path / "analysis" / "reports" / "quant" / f"{info['record_id']}.reference.csv"
    assert csv.is_file()
    assert info["reference"]["sha256"] == hashlib.sha256(csv.read_bytes()).hexdigest()
    assert info["reference"]["sha256"] != ref.sha256
    assert new.quant_state["unit"] == "pmol / 10^6 cells"
    row = ds.feature_ids.index("1")
    assert new.feature_matrix[row, 0] == pytest.approx(1000 / 500 * (10.0 * 100.0 / 2.0))


def test_skip_records_and_marks_state(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    new, info = qn.skip(ds, "positive")
    assert np.array_equal(new.feature_matrix, ds.feature_matrix)
    assert np.array_equal(new.detected_mask, ds.detected_mask)
    assert new.quant_state["normalized"] is False and new.quant_state["skipped"] == ["normalize"]
    rec = json.loads((tmp_path / "analysis" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["skipped"] is True and rec["step"] == "normalize"


def test_rules_unavailable_envelope(tmp_path):
    ds = _simple(tmp_path)
    ref = _ref(tmp_path)
    for text in (qn.preview_tool(ds, ref, None, "positive"),
                 qn.apply_tool(ds, ref, None, "positive", holder=None, parent_handle_id=None)):
        assert json.loads(text)["error"]["code"] == "quant_rules_unavailable"


def test_non_integer_feature_id_refused(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ds.feature_ids[1] = "positive:1"
    ds.feature_metadata["positive:1"] = ds.feature_metadata.pop("1")
    with pytest.raises(DomainError) as exc:
        qn.preview(ds, _ref(tmp_path), None, "positive")
    assert exc.value.code == "QUANT_FEATURE_ID_NOT_INTEGER"


def test_no_analysis_folder_is_coded(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ds.source_files = {}
    with pytest.raises(DomainError) as exc:
        qn.apply(ds, _ref(tmp_path), None, "positive")
    assert exc.value.code == "QUANT_RECORD_LOCATION_REQUIRED"


def test_tools_register_in_subprocess():
    code = (
        "from metabolomix.core.mcp_core import mcp\n"
        "import metabolomix.quant.tools as t\n"
        "tm = mcp._tool_manager._tools\n"
        "assert set(t.__all__) == {'dataset_quant_preview','dataset_quant_apply'}\n"
        "for n in t.__all__:\n"
        "    assert n in tm and tm[n].annotations is not None and tm[n].output_schema is None\n"
        "    assert 'dataset_handle' in tm[n].parameters['properties']\n"
        "print('ok')\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "ok" in out.stdout, out.stderr
