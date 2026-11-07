"""付加体の選別と除外（select / dataset_select_features）。人工の解析フォルダと偽の規則で。"""
from __future__ import annotations

import json
import subprocess
import sys

from metabolomix.core import handles, session_state
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant import derive as qd
from metabolomix.quant import inputs as qi
from metabolomix.quant import merge as qm
from metabolomix.quant import records as qr
from metabolomix.quant import select as qs
from tests.msdial_folder_fixtures import write_msdial_folder
from tests.quant_fake_rules import write_fake_rules

ADDUCTS = ("Class,Adduct,Ion mode,IsSelected\n"
           "PC,[M+H]+,Positive,True\nPC,[M+Na]+,Positive,False\nPE,[M-H]-,Negative,True\n")
NAMES = ("S1", "S2", "S3")


def _point(tmp_path, monkeypatch, adducts=True):
    folder = write_fake_rules(tmp_path / "rules")
    lines = [f"rules = '{folder}'"]
    if adducts:
        table = tmp_path / "adduct-table.csv"
        table.write_text(ADDUCTS, encoding="utf-8")
        lines.append(f"adducts = '{table}'")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text("[quant]\n" + "\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))


def _spot(name, ontology="PC", adduct="[M+H]+", mz=500.0):
    return {"name": name, "ontology": ontology, "adduct": adduct, "mz": mz,
            "heights": [10.0, 20.0, 30.0]}


def _load(folder, spots, names=NAMES):
    samples = [{"file_id": i, "name": n, "class_id": "A"} for i, n in enumerate(names)]
    write_msdial_folder(folder, spots=spots, samples=samples)
    ds = load_msdial_folder(select_alignment(folder))
    for meta in ds.feature_metadata.values():
        meta["msms_assigned"] = True
        meta["comment"] = ""
    return ds


def _normalized(ds, is_ids=("0",), unit="nmol"):
    ds.quant_state = qd.quant_state(
        normalized=True, method="internal_standard", unit=unit, record_id="q_norm",
        steps=["normalize"], skipped=[], internal_standard_ids=list(is_ids))
    return ds


def _pos(tmp_path, root="run"):
    ds = _load(tmp_path / root, [
        _spot("IS-A"),                       # 0 internal standard
        _spot("PC 34:1"),                    # 1 selected adduct
        _spot("PC 36:2", adduct="[M+Na]+"),  # 2 adduct not selected
        _spot("Unknown"),                    # 3 name
        _spot("PC 40:6"),                    # 4 no MS/MS
        _spot("PC 30:0"),                    # 5 comment
    ])
    ds.feature_metadata["4"]["msms_assigned"] = False
    ds.feature_metadata["5"]["comment"] = "x bad peak"
    return _normalized(ds)


def _adduct():
    return qi.adduct_table()


def _audit(tmp_path, info, root="run"):
    path = tmp_path / root / "reports" / "quant" / f"{info['record_id']}.audit.json"
    return json.loads(path.read_text(encoding="utf-8"))


def test_select_reasons_counted(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    assert info["reasons"] == {"internal_standard": 1, "reference_selected": 1,
                               "reference_not_selected": 1, "name_unknown": 1,
                               "msms_not_assigned": 1, "comment_excluded": 1}
    assert info["n_retained"] == 1 and info["n_excluded"] == 5
    assert new.feature_ids == ["1"] and new.feature_matrix.shape == (1, 3)
    assert ds.feature_ids == ["0", "1", "2", "3", "4", "5"]  # 元は不変
    assert new.quant_state["steps"] == ["normalize", "select"]
    assert new.quant_state["unit"] == "nmol" and new.quant_state["normalized"] is True


def test_internal_standard_rows_removed(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    _new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    rows = {d["feature_id"]: d for d in _audit(tmp_path, info)["decisions"]}
    assert rows["0"]["action"] == "excluded" and rows["0"]["reason"] == "internal_standard"
    assert rows["0"]["mode"] == "positive"


def test_not_normalized_role_unassigned(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    ds.quant_state = None  # 正規化していない: 内部標準の行を知らない
    new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    assert "internal_standard" not in info["reasons"]
    assert new.feature_ids == ["0", "1"]  # 0 は PC [M+H]+ として通る
    assert new.quant_state["normalized"] is False and new.quant_state["method"] is None
    assert new.quant_state["steps"] == ["select"] and new.quant_state["skipped"] == []
    ds2 = _pos(tmp_path, "run2")
    ds2.quant_state = qd.quant_state(normalized=False, method=None, unit=None, record_id="q_s",
                                     steps=[], skipped=["normalize"], internal_standard_ids=["0"])
    new2, info2 = qs.select(ds2, adduct=_adduct(), polarity="positive")
    assert "internal_standard" not in info2["reasons"] and new2.feature_ids == ["0", "1"]
    assert new2.quant_state["skipped"] == ["normalize"]


def test_curation_applied_before_rules(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    ds.curation_flags = {"wrong": {1}, "redundant": set(),
                         "assign": {2: {"name": "PC 36:2", "ontology": "PE"},
                                    5: {"name": "PC 30:0", "ontology": "PC"}}}
    new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    audit = _audit(tmp_path, info)
    dec = {d["feature_id"]: d for d in audit["decisions"]}
    assert dec["1"]["reason"] == "curation_excluded" and dec["1"]["action"] == "excluded"
    # assign は name / ontology を差し替えてから規則へ（PE の [M+H]+ は採用表に無い）
    assert dec["2"]["original_class"] == "PE" and dec["2"]["reason"] == "reference_not_selected"
    assert dec["5"]["reason"] == "comment_excluded"  # 置換でコメントは消えない
    assert new.feature_ids == []
    assert info["curation"] == {"excluded": 1, "replaced": 2}
    assert {(a["feature_id"], a["action"]) for a in audit["curation"]} == {
        ("1", "excluded"), ("2", "replaced"), ("5", "replaced")}
    # apply_curation=False: 判断を当てない
    ds.curation_flags = {"wrong": {1}}
    new2, info2 = qs.select(ds, adduct=_adduct(), polarity="positive", apply_curation=False)
    assert "curation_excluded" not in info2["reasons"] and "1" in new2.feature_ids
    assert info2["curation"] is None


def test_curation_keeps_internal_standard_reason(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    ds.curation_flags = {"wrong": {0}}
    _new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    assert info["reasons"]["internal_standard"] == 1  # 既に除外の行の理由は書き換えない
    assert "curation_excluded" not in info["reasons"]


def test_replaced_identity_carried_to_retained_rows(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    ds.curation_flags = {"assign": {1: {"name": "PC 35:1", "ontology": "PC"}}}
    new, _info = qs.select(ds, adduct=_adduct(), polarity="positive")
    assert new.feature_metadata["1"]["name"] == "PC 35:1"
    assert ds.feature_metadata["1"]["name"] == "PC 34:1"


def test_skip_filter_records_skip(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, adducts=False)
    ds = _pos(tmp_path)
    new, info = qs.select(ds, adduct=None, polarity="positive")
    assert info["adducts"] is None
    assert info["reasons"]["filter_disabled"] == 2  # PC 34:1 と PC 36:2。付加体は見ない
    assert info["reasons"]["internal_standard"] == 1
    assert "reference_not_selected" not in info["reasons"]
    assert "adduct_filter" in new.quant_state["skipped"]
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["step"] == "select" and rec["settings"]["adducts"] is None
    assert rec["settings"]["apply_curation"] is True
    assert rec["settings"]["rules"]["rules_version"] and rec["settings"]["rules"]["digest"]
    assert not any(k.startswith("adducts:") for k in rec["inputs"])
    assert rec["parent_ids"] == ["q_norm"]


def test_record_with_adduct_table(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    ref = _adduct()
    _new, info = qs.select(ds, adduct=ref, polarity="positive")
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert rec["inputs"]["adducts:adduct-table.csv"] == ref.sha256
    assert rec["settings"]["adducts"] == {"sha256": ref.sha256, "label": "adduct-table.csv"}
    assert info["adducts"] == {"sha256": ref.sha256, "label": "adduct-table.csv"}
    assert str(tmp_path) not in json.dumps(info)


def test_audit_file_written_not_returned(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    spots = [_spot("IS-A")] + [_spot(f"PC 3{i}:1", adduct="[M+Na]+") for i in range(9)]
    ds = _load(tmp_path / "run", spots)
    _normalized(ds)
    _new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    audit = _audit(tmp_path, info)
    assert len(audit["decisions"]) == 10 and audit["record_id"] == info["record_id"]
    assert len(info["examples"]) == 5
    assert "decisions" not in info
    assert json.dumps(info).count('"feature_id"') <= 5
    assert info["audit_file"] == f"{info['record_id']}.audit.json"
    assert str(tmp_path) not in json.dumps(info)


def test_read_records_ignores_audit(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    _new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    recs, warnings = qr.read_records(tmp_path / "run")
    assert [r["record_id"] for r in recs] == [info["record_id"]] and warnings == []


# ---- 統合後のデータセット ----

def _merged(tmp_path):
    pos = _load(tmp_path / "run" / "pos", [_spot("IS-A"), _spot("PC 34:1"),
                                           _spot("PC 36:2", adduct="[M+Na]+")],
                names=("S1_pos", "S2_pos", "S3_pos"))
    neg = _load(tmp_path / "run" / "neg", [_spot("IS-B", "PE", "[M-H]-"),
                                           _spot("PE 36:2", "PE", "[M-H]-"),
                                           _spot("PE 38:4", "PE", "[M-H]-")],
                names=("S1_neg", "S2_neg", "S3_neg"))
    _normalized(pos, ("0",))
    _normalized(neg, ("0",))
    neg.quant_state["record_id"] = "q_neg"
    pos.curation_flags = {"wrong": {1}}
    neg.curation_flags = {"wrong": {2}, "assign": {1: {"name": "PE 40:6", "ontology": "PE"}}}
    merged, _ = qm.merge(pos, neg, polarities=("positive", "negative"))
    return merged


def test_merged_dataset_curation_per_polarity(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    merged = _merged(tmp_path)
    assert merged.curation_flags is None
    assert merged.quant_state["internal_standard_ids"] == ["positive:0", "negative:0"]
    new, info = qs.select(merged, adduct=_adduct())
    audit = _audit(tmp_path, info)
    dec = {d["feature_id"]: d for d in audit["decisions"]}
    assert dec["positive:0"]["reason"] == "internal_standard"
    assert dec["negative:0"]["reason"] == "internal_standard"
    # 極性ごとの判断: positive の 1 は除外。negative の 1 は assign、negative の 2 は除外
    assert dec["positive:1"]["reason"] == "curation_excluded"
    assert dec["negative:1"]["reason"] == "reference_selected"
    assert dec["negative:1"]["mode"] == "negative"
    assert dec["negative:2"]["reason"] == "curation_excluded"
    assert dec["positive:2"]["reason"] == "reference_not_selected"
    assert new.feature_ids == ["negative:1"]
    assert new.feature_metadata["negative:1"]["name"] == "PE 40:6"
    assert new.feature_metadata["negative:1"]["origin"]["feature_id"] == "1"
    assert new.quant_state["steps"][-2:] == ["merge", "select"]
    assert new.quant_state["internal_standard_ids"] == []
    assert new.quant_state["curation_by_polarity"]["positive"] == {"wrong": {1}}
    assert {(a["feature_id"], a["action"]) for a in audit["curation"]} == {
        ("positive:1", "excluded"), ("negative:2", "excluded"), ("negative:1", "replaced")}
    recs, warnings = qr.read_records(tmp_path / "run")  # 統合後の記録は共通の親
    assert warnings == [] and recs[-1]["record_id"] == info["record_id"]
    assert recs[-1]["parent_ids"] == [merged.quant_state["record_id"]]


def test_merged_polarity_without_flags_reports_no_flags(tmp_path, monkeypatch):
    """統合後で apply_curation=True なのに、ある極性の判断が無いときは黙って飛ばさず理由を返す。"""
    _point(tmp_path, monkeypatch)
    merged = _merged(tmp_path)
    merged.quant_state["curation_by_polarity"].pop("negative")
    _new, info = qs.select(merged, adduct=_adduct())
    assert info["curation_by_polarity"]["negative"] == {"applied": False, "reason": "no_flags"}
    assert info["curation_by_polarity"]["positive"]["applied"] is True
    _new, off = qs.select(merged, adduct=_adduct(), apply_curation=False)
    assert off["curation_by_polarity"] is None  # 当てない指定なら報告するものが無い


def test_adduct_key_matches_record_inputs(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    assert qi.adduct_key("adduct-table.csv") == "adducts:adduct-table.csv"
    ds = _pos(tmp_path)
    _new, info = qs.select(ds, adduct=_adduct(), polarity="positive")
    rec = json.loads((tmp_path / "run" / "reports" / "quant"
                      / f"{info['record_id']}.json").read_text(encoding="utf-8"))
    assert qi.adduct_key(_adduct().label) in rec["inputs"]
    assert rec["settings"]["polarity"] == "positive"


# ---- ツール ----

def _open(ds, polarity):
    session_state.session.dataset = ds
    holder = handles.registry.get(session_state.session.current_handle_id)
    holder.polarity = polarity
    return holder.handle_id


def _tool(handle, **kw):
    return qs.select_tool(handles.registry.get(handle), parent_handle_id=handle, **kw)


def test_missing_adduct_table_needs_input(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, adducts=False)
    h = _open(_pos(tmp_path), "positive")
    out = json.loads(_tool(h))
    assert out["status"] == "needs_input" and out["reason"] == "adduct_table_required"
    assert "[quant] adducts" in out["how_to_answer"]
    assert "skip_adduct_filter=True" in out["how_to_answer"]
    assert not (tmp_path / "run" / "reports").exists()
    assert str(tmp_path) not in json.dumps(out)


def test_tool_selects_into_new_handle(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    h = _open(ds, "positive")
    out = json.loads(_tool(h))
    assert out["status"] == "selected" and out["dataset_handle"] != h
    assert out["n_retained"] == 1 and out["adducts"]["label"] == "adduct-table.csv"
    new = handles.registry.get(out["dataset_handle"])
    assert new.dataset.feature_ids == ["1"] and new.source["parents"] == [h]
    assert session_state.session.current_handle_id == out["dataset_handle"]
    assert handles.registry.get(h).dataset is ds
    assert str(tmp_path) not in json.dumps(out) and len(json.dumps(out)) < 2500
    assert (tmp_path / "run" / "reports" / "quant" / out["audit_file"]).is_file()


def test_tool_merged_dataset_needs_no_polarity(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    merged = _merged(tmp_path)
    h = _open(merged, None)
    out = json.loads(_tool(h))
    assert out["status"] == "selected" and out["n_retained"] == 1


def test_tool_skip_filter(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch, adducts=False)
    h = _open(_pos(tmp_path), "positive")
    out = json.loads(_tool(h, skip_adduct_filter=True))
    assert out["status"] == "selected" and out["adducts"] is None
    assert "adduct_filter" in out["skipped"]


def test_tool_polarity_required_for_unknown_single(tmp_path, monkeypatch):
    _point(tmp_path, monkeypatch)
    ds = _pos(tmp_path)
    h = _open(ds, None)
    monkeypatch.setattr("metabolomix.mztab.polarity.dataset_polarity",
                        lambda d: (None, None, []))
    out = json.loads(_tool(h))
    assert out["error"]["code"] == "QUANT_POLARITY_REQUIRED"
    ok = json.loads(_tool(h, polarity="positive"))
    assert ok["status"] == "selected"


def test_tool_missing_dataset_and_rules(tmp_path, monkeypatch):
    class Empty:
        dataset = None
    out = json.loads(qs.select_tool(Empty(), parent_handle_id=None))
    assert out["error"]["code"] == "missing_state"
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(tmp_path / "absent.toml"))
    h = _open(_pos(tmp_path), "positive")
    out = json.loads(_tool(h, skip_adduct_filter=True))
    assert out["error"]["code"] == "quant_rules_unavailable"


def test_select_tool_registers_in_subprocess():
    code = (
        "from metabolomix.core.mcp_core import mcp\n"
        "import metabolomix.quant.tools as t\n"
        "tm = mcp._tool_manager._tools\n"
        "n = 'dataset_select_features'\n"
        "assert n in t.__all__ and n in tm\n"
        "assert tm[n].annotations is not None and tm[n].output_schema is None\n"
        "assert {'skip_adduct_filter','apply_curation','polarity','dataset_handle'} "
        "<= set(tm[n].parameters['properties'])\n"
        "print('ok')\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "ok" in out.stdout, out.stderr
