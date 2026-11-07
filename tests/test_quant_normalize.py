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
from metabolomix.quant import derive as qd
from metabolomix.quant import normalize as qn
from metabolomix.quant.inputs import reference_input
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.quant_fake_rules import write_fake_rules

UNIT = "pmol / 10^6 cells"  # 本物の規則も受け付ける書き方
REF = (f"TargetClass,StandardName,Coefficient,Unit\nPC,IS-A,2.0,{UNIT}\n"
       f"PE,IS-B,3.0,{UNIT}\nAny others,IS-A,1.0,{UNIT}\n")
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
    ready = qn.preview(ds, ref, {"IS-B": 2}, "positive")
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
    assert new.quant_state["unit"] == UNIT and info["unit"] == UNIT
    assert new.quant_state["internal_standard_ids"] == ["0"]
    assert ds.quant_state is None and ds.feature_matrix[row, 0] == 1000.0  # 元は不変


def test_apply_with_ambiguous_binding(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _ambiguous(tmp_path)
    new, info = qn.apply(ds, _ref(tmp_path), {"IS-B": 2}, "positive")
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
    other = next(c for c in pv["classes"] if c["target_class"] == "Any others")
    assert other["via_any_others"] is True and other["dataset_classes"] == ["TG"]


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


def test_record_next_to_arf2_even_if_mddata_elsewhere(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ds.source_files[str(tmp_path / "elsewhere" / "proj.mddata")] = "0" * 64
    _new, info = qn.apply(ds, _ref(tmp_path), None, "positive")
    assert (tmp_path / "analysis" / "reports" / "quant" / f"{info['record_id']}.json").is_file()
    assert not (tmp_path / "elsewhere").exists()


def test_analysis_folder_rules(tmp_path):
    ds = _simple(tmp_path)
    arf2 = next(p for p in ds.source_files if p.endswith(".arf2"))
    arf = next(p for p in ds.source_files if p.endswith(".arf"))
    assert qn.analysis_folder(ds) == tmp_path / "analysis"
    ds.source_files = {arf: "h", str(tmp_path / "o" / "a.mddata"): "h"}  # arf2 無し → arf
    assert qn.analysis_folder(ds) == tmp_path / "analysis"
    ds.source_files = {str(tmp_path / "z" / "x.mzTab"): "h", str(tmp_path / "o" / "a.mddata"): "h"}
    with pytest.raises(DomainError) as exc:  # mzTab-M は解析フォルダの主ファイルではない
        qn.analysis_folder(ds)
    assert exc.value.code == "QUANT_RECORD_LOCATION_REQUIRED"
    assert "mzTab" not in exc.value.message
    ds.source_files = {str(tmp_path / "o" / "a.mddata"): "h"}  # 主ファイル無し
    with pytest.raises(DomainError) as exc:
        qn.analysis_folder(ds)
    assert exc.value.code == "QUANT_RECORD_LOCATION_REQUIRED"
    ds.source_files = {arf2: "h", str(tmp_path / "other" / "b.arf2"): "h"}  # 主ファイルが別フォルダ
    with pytest.raises(DomainError) as exc:
        qn.analysis_folder(ds)
    assert exc.value.code == "QUANT_RECORD_LOCATION_REQUIRED"


def test_mask_follows_resolved_class_alias_and_any_others(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _load(tmp_path, [
        _spot("IS-A", "PC", [500.0, 500.0, 500.0], gap=[False, True, False]),
        _spot("Cer alias", "Cer_NS", [1000.0] * 3, mz=600.0),
        _spot("TG 50:1", "TG", [100.0] * 3, mz=800.0)])
    ref = _ref(tmp_path, text=REF + f"Cer-NS,IS-A,2.0,{UNIT}\n")  # Cer_NS は Cer-NS の別名
    new, _ = qn.apply(ds, ref, None, "positive")
    for fid, val in (("1", 4.0), ("2", 0.2)):
        row = ds.feature_ids.index(fid)
        assert new.feature_matrix[row, 0] == pytest.approx(val)
        assert [bool(x) for x in new.detected_mask[row]] == [True, False, True]


def test_record_inputs_hold_no_reference_path(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ref = _ref(tmp_path)
    _new, info = qn.apply(ds, ref, None, "positive")
    rec = json.loads((tmp_path / "analysis" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["inputs"]["reference:ref_table.csv"] == ref.sha256
    assert not any(str(tmp_path / "ref_table.csv") == k for k in rec["inputs"])


def test_reference_read_error_carries_label_only(tmp_path, monkeypatch):
    from pathlib import Path
    path = tmp_path / "ref_table.csv"
    path.write_text(REF, encoding="utf-8")

    def boom(self):
        raise PermissionError(13, "denied", str(self))
    monkeypatch.setattr(Path, "read_bytes", boom)
    with pytest.raises(DomainError) as exc:
        reference_input(str(path))
    assert str(tmp_path) not in exc.value.message and str(tmp_path) not in str(exc.value.details)
    assert exc.value.details.get("label") == "ref_table.csv"


def test_apply_without_reference_needs_skip(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    out = qn.apply_tool(ds, None, None, "positive", holder=None, parent_handle_id=None)
    assert json.loads(out)["error"]["code"] == "QUANT_BAD_REQUEST"


def test_tools_register_in_subprocess():
    code = (
        "from metabolomix.core.mcp_core import mcp\n"
        "import metabolomix.quant.tools as t\n"
        "tm = mcp._tool_manager._tools\n"
        "assert {'dataset_quant_preview','dataset_quant_apply'} <= set(t.__all__)\n"
        "for n in ('dataset_quant_preview','dataset_quant_apply'):\n"
        "    assert n in tm and tm[n].annotations is not None and tm[n].output_schema is None\n"
        "    assert 'dataset_handle' in tm[n].parameters['properties']\n"
        "print('ok')\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "ok" in out.stdout, out.stderr


def test_reference_rows_for_absent_classes_do_not_block(tmp_path, monkeypatch):
    """参照表に載る内部標準が、データセットに無いクラスのぶんは needs_input にしない（本物の規則は
    参照表の全行の assignments を返すため、現れたクラスだけに絞る）。"""
    _point(tmp_path, monkeypatch)
    ds = _load(tmp_path, [_spot("IS-A", "PC", [500.0] * 3), _spot("PC 34:1", "PC", [1000.0] * 3)])
    ref = _ref(tmp_path, text=("TargetClass,StandardName,Coefficient,Unit\n"
                               f"PC,IS-A,2.0,{UNIT}\nPE,IS-B,3.0,{UNIT}\n"))
    out = qn.preview(ds, ref, None, "positive")
    assert out["status"] == "ready" and out["unresolved"] == []
    assert [c["target_class"] for c in out["classes"]] == ["PC"]
    new, info = qn.apply(ds, ref, None, "positive")
    assert new is not None and info["classes"][0]["target_class"] == "PC"


# ---- 入力の形式・工程の順序・極性の引き継ぎ・ハンドルの上限 ----

def test_non_msdial_source_refused(tmp_path, monkeypatch):
    """mzTab-M 由来はクラス・付加体・コメント・MS/MS の欄の意味が違うので定量しない。"""
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ds.source_format = "mztab"
    for call in (lambda: qn.preview(ds, _ref(tmp_path), None, "positive"),
                 lambda: qn.apply(ds, _ref(tmp_path), None, "positive"),
                 lambda: qn.skip(ds, "positive")):
        with pytest.raises(DomainError) as exc:
            call()
        assert exc.value.code == "QUANT_SOURCE_UNSUPPORTED"
    out = json.loads(qn.apply_tool(ds, _ref(tmp_path), None, "positive", holder=None,
                                   parent_handle_id=None))
    assert out["error"]["code"] == "QUANT_SOURCE_UNSUPPORTED"
    assert not (tmp_path / "analysis" / "reports").exists()


def test_apply_refused_after_select(tmp_path, monkeypatch):
    """選別で行を落としたデータセットは正規化しない（内部標準の行がもう無いかもしれない）。"""
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    ds.quant_state = qd.quant_state(normalized=False, method=None, unit=None, record_id="q_s",
                                    steps=["select"], skipped=[])
    for call in (lambda: qn.preview(ds, _ref(tmp_path), None, "positive"),
                 lambda: qn.apply(ds, _ref(tmp_path), None, "positive"),
                 lambda: qn.skip(ds, "positive")):
        with pytest.raises(DomainError) as exc:
            call()
        assert exc.value.code == "QUANT_ALREADY_SELECTED"


def _unknown_polarity(monkeypatch):
    for target in ("metabolomix.mztab.polarity.dataset_polarity",
                   "metabolomix.quant.derive.dataset_polarity"):
        monkeypatch.setattr(target, lambda d: (None, "unknown", []))


def test_new_handle_keeps_resolved_polarity(tmp_path, monkeypatch):
    """極性の決まらないデータセットに polarity を渡して正規化したら、新しいハンドルもその極性。"""
    _point(tmp_path, monkeypatch)
    _unknown_polarity(monkeypatch)
    for skip_step, root in ((False, "a"), (True, "b")):
        ds = _simple(tmp_path / root)
        handle = handles.DatasetHandle(ds, polarity=None, polarity_source="test",
                                       source={"kind": "test"})
        handles.registry.add(handle)
        out = json.loads(qn.apply_tool(ds, None if skip_step else _ref(tmp_path), None,
                                       "positive", holder=handle,
                                       parent_handle_id=handle.handle_id, skip_step=skip_step))
        assert handles.registry.get(out["dataset_handle"]).polarity == "positive", skip_step


def test_handle_capacity_checked_before_record(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _simple(tmp_path)
    session_state.session.dataset = ds
    parent = session_state.session.current_handle_id
    monkeypatch.setenv("LIPIDMIX_MAX_DATASET_HANDLES", "1")
    for skip_step in (False, True):
        out = json.loads(qn.apply_tool(ds, None if skip_step else _ref(tmp_path), None,
                                       "positive", holder=session_state.session,
                                       parent_handle_id=parent, skip_step=skip_step))
        assert out["error"]["code"] == "DATASET_HANDLE_LIMIT", skip_step
    assert not (tmp_path / "analysis" / "reports").exists()  # 登録できない記録を残さない
