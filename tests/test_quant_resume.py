"""定量の鎖の復元（resume.restore_chain）。人工の解析フォルダと、偽（`--quant-rules` 付きなら本物）の規則で。"""
from __future__ import annotations

import json

import numpy as np

from metabolomix.core import handles, session_state
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant import inputs as qi
from metabolomix.quant import normalize as qn
from metabolomix.quant import records as qr
from metabolomix.quant import resume
from metabolomix.quant import select as qs
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.quant_fake_rules import write_fake_rules

UNIT = "pmol / 10^6 cells"
REF = f"TargetClass,StandardName,Coefficient,Unit\nPC,IS-A,2.0,{UNIT}\n"
ADDUCTS = "Class,Adduct,Ion mode,IsSelected\nPC,[M+H]+,Positive,True\nPC,[M+Na]+,Positive,False\n"
SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A"} for i in range(3)]


def _point(tmp_path, *, rules=True, adducts=True):
    lines = []
    if rules:
        lines.append(f"rules = '{write_fake_rules(tmp_path / 'rules')}'")
    if adducts:
        table = tmp_path / "adduct-table.csv"
        table.write_text(ADDUCTS, encoding="utf-8")
        lines.append(f"adducts = '{table}'")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[quant]\n" + "\n".join(lines) + "\n", encoding="utf-8")
    return cfg


def _env(tmp_path, monkeypatch, **kw):
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(_point(tmp_path, **kw)))


def _base(tmp_path):
    spots = [{"name": "IS-A", "ontology": "PC", "heights": [500.0, 250.0, 0.0],
              "gap_filled": [False, True, False]},
             {"name": "PC 34:1", "ontology": "PC", "heights": [1000.0, 800.0, 300.0], "mz": 600.0},
             {"name": "PC 36:2", "ontology": "PC", "adduct": "[M+Na]+",
              "heights": [10.0, 20.0, 30.0], "mz": 700.0},
             {"name": "PC 38:4", "ontology": "PC", "heights": [5.0, 6.0, 7.0], "mz": 800.0}]
    write_msdial_folder(tmp_path / "analysis", spots=spots, samples=SAMPLES)
    ds = load_msdial_folder(select_alignment(tmp_path / "analysis"))
    for meta in ds.feature_metadata.values():
        meta["msms_assigned"] = True
        meta["comment"] = ""
    ds.curation_flags = {"wrong": {3}}
    return ds


def _holder(ds):
    session_state.session.dataset = ds
    return handles.registry.get(session_state.session.current_handle_id)


def _ref_file(tmp_path, text=REF):
    path = tmp_path / "lab_reference.csv"
    path.write_text(text, encoding="utf-8")
    return path


def _chain(tmp_path, ds):
    """normalize → select を本当に走らせて記録を書く。"""
    ref = qi.reference_input(str(_ref_file(tmp_path)))
    n1, _ = qn.apply(ds, ref, None, "positive")
    n2, info = qs.select(n1, adduct=qi.adduct_table(), polarity="positive")
    return n1, n2, info


def _folder(tmp_path):
    return tmp_path / "analysis"


def _derived(entry):
    return handles.registry.get(entry["dataset_handle"]).dataset


def _refs(tmp_path):
    return [tmp_path / "lab_reference.csv"]


def test_restore_rebuilds_same_values(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    n1, n2, _info = _chain(tmp_path, ds)
    holder = _holder(ds)
    before = len(qr.read_records(_folder(tmp_path))[0])
    out = resume.restore_chain(_folder(tmp_path), ds, holder, reference_paths=_refs(tmp_path))
    assert [e["step"] for e in out["restored"]] == ["normalize", "select"]
    assert out["stale"] == [] and out["skipped"] == []
    r1, r2 = (_derived(e) for e in out["restored"])
    assert np.array_equal(r1.feature_matrix, n1.feature_matrix, equal_nan=True)
    assert np.array_equal(r1.detected_mask, n1.detected_mask)
    assert np.array_equal(r2.feature_matrix, n2.feature_matrix, equal_nan=True)
    assert np.array_equal(r2.detected_mask, n2.detected_mask)
    assert r2.feature_ids == n2.feature_ids and r2.feature_ids  # 空でない選別結果
    assert r2.quant_state == n2.quant_state and r1.quant_state == n1.quant_state
    assert len(qr.read_records(_folder(tmp_path))[0]) == before  # 記録を増やさない
    assert session_state.session.current_handle_id == out["restored"][-1]["dataset_handle"]
    again = resume.restore_chain(_folder(tmp_path), ds, holder,  # 決定論的
                                 reference_paths=_refs(tmp_path))
    assert np.array_equal(_derived(again["restored"][-1]).feature_matrix, r2.feature_matrix,
                          equal_nan=True)
    assert ds.quant_state is None  # 元は不変


def test_restore_skips_stale_chain(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    holder = _holder(ds)
    key = next(k for k in ds.source_files if k.lower().endswith(".arf2"))
    ds.source_files[key] = "0" * 64  # .arf2 が変わった
    out = resume.restore_chain(_folder(tmp_path), ds, holder, reference_paths=_refs(tmp_path))
    assert out["restored"] == []
    assert [s["step"] for s in out["stale"]] == ["normalize", "select"]
    assert out["stale"][0]["reason"].startswith("input_changed:")
    assert out["stale"][1]["reason"].startswith("input_changed:")  # 選別の入力も同じ .arf2


def test_fresh_ignores_records(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    out = resume.restore_chain(_folder(tmp_path), ds, _holder(ds), fresh=True,
                               reference_paths=_refs(tmp_path))
    assert out == {"restored": [], "stale": [], "skipped": []}


def test_no_records_restores_nothing(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    assert resume.restore_chain(_folder(tmp_path), ds, _holder(ds)) == {
        "restored": [], "stale": [], "skipped": []}


def test_reference_missing_or_changed_stops_chain(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    holder = _holder(ds)
    out = resume.restore_chain(_folder(tmp_path), ds, holder)  # 参照表の場所を渡さない
    assert out["restored"] == [] and out["stale"] == []
    assert out["skipped"][0]["step"] == "normalize"
    assert out["skipped"][0]["reason"] == "reference_missing"
    assert str(tmp_path) not in json.dumps(out)
    _ref_file(tmp_path, REF.replace("2.0", "3.0"))  # 同じ名前で中身が違う
    out = resume.restore_chain(_folder(tmp_path), ds, holder, reference_paths=_refs(tmp_path))
    assert out["restored"] == [] and out["skipped"][0]["reason"] == "reference_missing"


def test_concentration_form_restores_from_saved_table(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    std = "TargetClass,StandardName,Unit\nPC,IS-A,pmol\n"
    ref = qi.reference_input(str(_ref_file(tmp_path, std)), concentrations={"IS-A": 2.0},
                             solvent_volume=10.0, cells_millions=4.0)
    n1, _ = qn.apply(ds, ref, None, "positive")
    (tmp_path / "lab_reference.csv").unlink()  # 保存した係数表だけで戻る（元の表は要らない）
    out = resume.restore_chain(_folder(tmp_path), ds, _holder(ds))
    assert [e["step"] for e in out["restored"]] == ["normalize"]
    assert np.array_equal(_derived(out["restored"][0]).feature_matrix, n1.feature_matrix,
                          equal_nan=True)


def test_rules_digest_mismatch_marks_stale(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    for path in (_folder(tmp_path) / "reports" / "quant").glob("*.json"):
        if path.name.endswith(".audit.json"):
            continue
        rec = json.loads(path.read_text(encoding="utf-8"))
        rec["settings"]["rules"]["digest"] = "0" * 64  # 記録が使った規則と違う
        path.write_text(json.dumps(rec), encoding="utf-8")
    out = resume.restore_chain(_folder(tmp_path), ds, _holder(ds), reference_paths=_refs(tmp_path))
    assert out["restored"] == []
    assert out["stale"][0]["step"] == "normalize"
    assert out["stale"][0]["reason"] == "rules_digest_changed"


def test_rules_unset_stops_chain(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    holder = _holder(ds)
    _env(tmp_path, monkeypatch, rules=False)
    out = resume.restore_chain(_folder(tmp_path), ds, holder, reference_paths=_refs(tmp_path))
    assert out["restored"] == [] and out["skipped"][0]["reason"] == "rules_unavailable"


def test_skipped_normalize_then_select_restores(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    n1, _ = qn.skip(ds, "positive")
    n2, _ = qs.select(n1, adduct=qi.adduct_table(), polarity="positive")
    out = resume.restore_chain(_folder(tmp_path), ds, _holder(ds))
    assert [e["step"] for e in out["restored"]] == ["normalize", "select"]
    assert _derived(out["restored"][1]).quant_state == n2.quant_state
    assert np.array_equal(_derived(out["restored"][1]).feature_matrix, n2.feature_matrix,
                          equal_nan=True)


def test_adduct_table_change_marks_select_stale(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    _chain(tmp_path, ds)
    (tmp_path / "adduct-table.csv").write_text(ADDUCTS.replace("False", "True"), encoding="utf-8")
    out = resume.restore_chain(_folder(tmp_path), ds, _holder(ds), reference_paths=_refs(tmp_path))
    assert [e["step"] for e in out["restored"]] == ["normalize"]
    assert out["stale"][0]["step"] == "select"
    assert out["stale"][0]["reason"].startswith("input_changed:")


def test_merge_records_are_skipped_with_reason(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    ds = _base(tmp_path)
    folder = _folder(tmp_path)
    merge = qr.write_record(folder, "merge", inputs={}, parent_ids=["q_a", "q_b"],
                            settings={}, summary={})
    post = qr.write_record(folder, "select", inputs={}, parent_ids=[merge["record_id"]],
                           settings={}, summary={})
    out = resume.restore_chain(folder, ds, _holder(ds))
    assert out["restored"] == [] and out["stale"] == []
    assert {(s["record_id"], s["reason"]) for s in out["skipped"]} == {
        (merge["record_id"], "merge_restore_needs_both_polarities"),
        (post["record_id"], "merge_restore_needs_both_polarities")}
