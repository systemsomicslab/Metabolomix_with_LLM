"""再解析のツール（spec 2026-10-22 §4 の 3〜11）を MCP のツールとして呼ぶ。

ハンドルの上でも対話のセッションでも動き、どちらの状態も書き換えないこと・同じ spec なら
同じ JSON になることを確かめる。fixture はすべてテスト自身が作る。
"""
from __future__ import annotations

import json
import textwrap

import pytest

from metabolomix.core import handles, mcp_core, session_state
from metabolomix.dcl import tools as dcl_tools
from metabolomix.msdial import alignment_spot
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools, mztab_tools
from metabolomix.tools import reanalysis_tools as R
from tests.dcl_fixture import build_dcl_bytes
from tests.reanalysis_fixtures import SAMPLES, write_mztab

GABA_KEY = "BTCSSZJGUNDROE-UHFFFAOYSA-N"
BASE = [10, 12, 11, 20, 22, 21, 15, 16]


def _scaled(factor):
    return [round(v * factor, 4) for v in BASE]


#: GABA の [M+H]+ と、同時に溶出する [M+Na]+・M+1 同位体・-H2O の源内断片・無関係の特徴。
RELATION_FEATURES = (
    {"id": "20", "mz": 104.0706, "rt": 2.00, "values": _scaled(100),
     "sme": [{"name": "GABA", "inchikey": GABA_KEY, "formula": "C4H9NO2", "adduct": "[M+H]1+",
              "theo": 104.07061, "exp": 104.0707}]},
    {"id": "21", "mz": 126.0526, "rt": 2.01, "values": _scaled(30), "sme": []},
    {"id": "22", "mz": 105.0740, "rt": 2.00, "values": _scaled(5), "sme": []},
    {"id": "23", "mz": 86.0600, "rt": 2.02, "values": _scaled(8), "sme": []},
    {"id": "24", "mz": 300.0000, "rt": 2.00, "values": [5, 1, 9, 2, 8, 3, 4, 6], "sme": []},
    {"id": "25", "mz": 126.0526, "rt": 5.00, "values": _scaled(2), "sme": []},
    {"id": "26", "mz": 250.1000, "rt": 6.00, "values": [9, 10, 8, 1, 2, 1, 5, 5],
     "sme": [{"name": "Other", "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "formula": "C10H15NO5",
              "adduct": "[M+H]1+", "theo": 250.09997, "exp": 250.1000}]},
)

_MSP = textwrap.dedent("""\
    NAME: GABA
    PRECURSORMZ: 104.0706
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: C4H9NO2
    INCHIKEY: BTCSSZJGUNDROE-UHFFFAOYSA-N
    Num Peaks: 3
    87.0441 999
    69.0335 500
    104.0706 100
""")


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def data(tmp_path):
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=RELATION_FEATURES)
    groups = tmp_path / "groups.tsv"
    groups.write_text("sample\tgroup\trole\n" + "".join(
        f"{s}\t{s[:2] if not s.startswith('QC') else ''}\t{'qc' if s.startswith('QC') else 'sample'}\n"
        for s in SAMPLES), encoding="utf-8")
    return {"mztab": mztab, "groups": groups}


@pytest.fixture
def library(tmp_path, monkeypatch, data):
    """研究室のライブラリの代わりの小さな .msp と、spot_id = SMF_ID の .arf2 / .dcl。"""
    msp = tmp_path / "lab" / "pos.msp"
    msp.parent.mkdir()
    msp.write_text(_MSP, encoding="ascii")
    monkeypatch.setenv("MSDIAL_MSP_POS", str(msp))
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    run = data["mztab"].parent
    (run / "AlignResult-1.arf2").write_bytes(b"placeholder")
    entries = [{"precursor_mz": 999.0, "rt": 9.0, "spectrum": []} for _ in range(27)]
    entries[20] = {"precursor_mz": 104.0706, "rt": 2.0,
                   "spectrum": [(87.0441, 999.0), (69.0335, 400.0), (104.0706, 300.0)]}
    entries[21] = {"precursor_mz": 126.0526, "rt": 2.01, "spectrum": [(104.07, 50.0)]}
    (run / "AlignResult-1.dcl").write_bytes(build_dcl_bytes(entries))
    records = [{"AlignmentID": int(f["id"]), "MassCenter": f["mz"], "RT": f["rt"],
                "IonMode": "Positive", "Name": "x"} for f in RELATION_FEATURES]
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: records)
    return msp


def _open(data, **kwargs):
    kwargs.setdefault("library", "none")
    out = json.loads(handle_tools.dataset_handle_open(mztab_path=str(data["mztab"]), **kwargs))
    assert out["status"] == "success", out
    return out["dataset_handle"]


def _spec(data, **extra):
    spec = {"groups": {"source": "tsv", "path": str(data["groups"]), "role_column": "role"},
            "comparison": {"reference": "HC", "test": "RA"},
            "features": {"feature_ids": ["20", "26"]},
            "analyses": [{"analysis_id": "main", "test": {"kind": "student"},
                          "thresholds": {"p": 0.05, "fc": 1.5}}]}
    spec.update(extra)
    return spec


# ---------- dataset_analysis_run ----------

def test_analysis_run_on_handle_returns_rows_and_leaves_states(data):
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    assert out["status"] == "success" and out["comparison"] == "RA_vs_HC"
    rows = out["analyses"][0]["rows"].splitlines()
    assert rows[0].startswith("feature_id\tstatus\treason\tgroup_ref\tgroup_test")
    assert rows[1].split("\t")[:5] == ["20", "tested", "", "HC", "RA"]
    assert session_state.session.dataset is None
    assert handles.registry.get(handle).dataset.pp_matrix is None


def test_analysis_run_on_session_dataset(data):
    mztab_tools.dataset_load(mztab_path=str(data["mztab"]))
    out = json.loads(R.dataset_analysis_run(_spec(data)))
    assert out["status"] == "success"
    assert out["input_fingerprint"]["dataset_hash"]


def test_analysis_run_without_dataset_is_missing_state(data):
    out = json.loads(R.dataset_analysis_run(_spec(data)))
    assert out["error"]["code"] == "missing_state"


def test_analysis_run_output_path_is_deterministic(data, tmp_path):
    handle = _open(data)
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    first = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle, output_path=str(a)))
    R.dataset_analysis_run(_spec(data), dataset_handle=handle, output_path=str(b))
    assert a.read_bytes() == b.read_bytes()
    assert "rows" not in first["analyses"][0]
    saved = json.loads(a.read_text(encoding="utf-8"))
    assert saved["rule_hash"] == first["rule_hash"] and saved["schema"] == "analysis-result.v1"


def test_analysis_run_requires_output_path_for_many_rows(tmp_path):
    matrix = tmp_path / "m.tsv"
    lines = ["feature_id\t" + "\t".join(SAMPLES)]
    lines += [f"F{i}\t" + "\t".join(str(v + i) for v in BASE) for i in range(51)]
    matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")
    groups = tmp_path / "g.tsv"
    groups.write_text("sample\tgroup\n" + "".join(f"{s}\t{s[:2]}\n" for s in SAMPLES),
                      encoding="utf-8")
    spec = {"input": {"source": "tsv", "matrix_path": str(matrix)},
            "groups": {"source": "tsv", "path": str(groups)},
            "comparison": {"reference": "HC", "test": "RA"},
            "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]}
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["error"]["code"] == "ANALYSIS_OUTPUT_REQUIRED"
    ok = json.loads(R.dataset_analysis_run(spec, output_path=str(tmp_path / "o.json")))
    assert ok["analyses"][0]["n_rows"] == 51


def test_analysis_run_bad_spec_is_an_error_envelope(data):
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run({"comparison": {}}, dataset_handle=handle))
    assert out["error"]["code"] == "ANALYSIS_SPEC_INVALID"


# ---------- 質量精度・新しい知見の候補 ----------

def test_mass_accuracy_by_band(data):
    handle = _open(data)
    out = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))
    assert out["overall"]["n"] == 2
    low, high = out["bands"]
    assert (low["min_mz"], low["max_mz"], low["n"]) == (None, 200.0, 1)
    assert high["n"] == 1


def test_novel_candidates_from_result_id_and_from_spec(data):
    handle = _open(data)
    T.dataset_preprocess(dataset_handle=handle)
    diff = json.loads(T.dataset_differential(group_a=["HC1", "HC2", "HC3"],
                                             group_b=["RA1", "RA2", "RA3"], dataset_handle=handle))
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    out = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert out["status"] == "success" and out["n_rows"] == 2
    assert out["group_ref"] == "group_a"
    excluded = json.loads(R.dataset_novel_candidates(
        {**spec, "exclude": {"inchikey_blocks": [GABA_KEY]}}, dataset_handle=handle))
    assert excluded["n_rows"] == 1
    from_spec = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": _spec(data, features={"all": True}),
                        "analysis_id": "main"}, "q_max": 1.0}, dataset_handle=handle))
    assert from_spec["group_test"] == "RA" and from_spec["n_rows"] == 2
    assert out["rule_hash"] != from_spec["rule_hash"]


# ---------- 関係 ----------

def test_feature_relations_name_adduct_isotope_and_insource(data):
    handle = _open(data)
    out = json.loads(R.dataset_feature_relations("20", dataset_handle=handle))
    coeluting = out["coeluting"].splitlines()
    assert {line.split("\t")[0] for line in coeluting[1:]} == {"21", "22", "23", "24"}
    ions = {line.split("\t")[0]: line.split("\t")[1] for line in out["ions_of_this"].splitlines()[1:]}
    assert ions["21"].startswith("adduct:[M+Na]+")
    assert ions["22"] == "isotope_M+1"
    assert ions["23"] == "insource:-H2O"
    as_ion = json.loads(R.dataset_feature_relations("21", dataset_handle=handle))["as_ion_of"]
    assert "adduct:[M+Na]+/[M+H]1+" in as_ion          # 相手の付加体は MS-DIAL の表記のまま


def test_feature_relations_unknown_feature(data):
    handle = _open(data)
    out = json.loads(R.dataset_feature_relations("999", dataset_handle=handle))
    assert out["error"]["code"] == "FEATURE_NOT_FOUND"


# ---------- 証拠の束 ----------

CLAIM = {"claim_id": "C1",
         "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0},
         "scan": {"ppm": [10, 20], "rt_tol_min": [0.3]},
         "relations": {"rt_window": 0.05}}


def test_claim_evidence_default_matches_feature_search(data):
    handle = _open(data)
    out = json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=handle))
    assert out["status"] == "success"
    assert out["default_selection"]["status"] == "resolved"
    assert out["default_selection"]["selected_feature_id"] == "20"
    search = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", "inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0}],
        dataset_handle=handle))
    assert search["targets"][0]["selected_feature_id"] == "20"
    assert [s["ppm"] for s in out["scans"]] == [10.0, 20.0]
    by_id = {c["feature_id"]: c for c in out["candidates"]}
    assert "20" in by_id and "21" in by_id            # [M+Na]+ は第 2 層の候補
    assert "msms" not in by_id["20"]                    # msms を頼んでいない
    assert by_id["20"]["relations"]["n_coeluting"] == 4


def test_claim_evidence_is_deterministic(data, tmp_path):
    handle = _open(data)
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    R.dataset_claim_evidence(CLAIM, dataset_handle=handle, output_path=str(a))
    R.dataset_claim_evidence(CLAIM, dataset_handle=handle, output_path=str(b))
    assert a.read_bytes() == b.read_bytes()
    saved = json.loads(a.read_text(encoding="utf-8"))
    assert saved["schema"] == "claim-evidence.v1" and len(saved["rule_hash"]) == 64


def test_claim_evidence_with_library_runs_msms_rules(data, library):
    handle = _open(data, library="auto")
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2,
                                                     "diagnostic_ions": [{"mz": 87.0441}]}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    msms = {c["feature_id"]: c["msms"] for c in out["candidates"]}
    assert msms["20"]["qualify"]["verdict"] == "qualified"
    assert msms["20"]["top"][0]["fragment_matches_excl_precursor"] == 2
    assert out["input_fingerprint"]["library_sha256"]


def test_claim_evidence_without_library_says_so(data):
    handle = _open(data)
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert {c["msms"]["status"] for c in out["candidates"]} == {"no_library"}


def test_claim_evidence_unresolved_is_needs_input(data):
    handle = _open(data)
    claim = {"claim_id": "C2", "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2"}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["default_selection"]["status"] in ("resolved", "needs_input")
    claim = {"claim_id": "C3", "target": {"inchikey": "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",
                                          "formula": "C4H9NO2", "rt_min": 2.0}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["default_selection"]["status"] == "needs_input"


# ---------- MS/MS の合否・診断イオン・同族列 ----------

def test_msms_qualify_tool_on_handle(data, library):
    handle = _open(data, library="auto")
    out = json.loads(R.library_msms_qualify(
        {"target": {"inchikey": GABA_KEY}, "min_fragments_excl_precursor": 2},
        spot_id=20, dataset_handle=handle))
    assert out["verdict"] == "qualified" and len(out["rule_hash"]) == 64
    assert out["candidates"].splitlines()[1].split("\t")[1] == "GABA"
    missing = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=24,
                                                dataset_handle=handle))
    assert missing["verdict"] == "not_evaluable"


def test_diagnostic_ions_tool(data, library):
    handle = _open(data, library="auto")
    out = json.loads(dcl_tools.dcl_diagnostic_ions(
        [20, 21], [{"mz": 87.0441, "label": "a"}, {"neutral_loss": 17.0265, "label": "-NH3"}],
        dataset_handle=handle))
    rows = [line.split("\t") for line in out["ions"].splitlines()[1:]]
    found = {(r[0], r[1]): r[4] for r in rows}
    assert found[("20", "a")] == "true" and found[("21", "a")] == "false"
    assert found[("20", "-NH3")] == "true"               # 104.0706 − 17.0265 = 87.0441


def test_homolog_check_tool_with_members(data):
    handle = _open(data)
    out = json.loads(R.dataset_homolog_check(
        {"members": [{"feature_id": "20", "carbons": 34, "double_bonds": 1},
                     {"feature_id": "26", "carbons": 36, "double_bonds": 1}]},
        dataset_handle=handle))
    assert out["verdict"] == "pass" and out["n_violations"] == 0
    bad = json.loads(R.dataset_homolog_check(
        {"members": [{"feature_id": "26", "carbons": 34, "double_bonds": 1},
                     {"feature_id": "20", "carbons": 36, "double_bonds": 1}]},
        dataset_handle=handle))
    assert bad["verdict"] == "fail"


# ---------- dataset_alignment_qc ----------

#: 最初の 2 注入（ファイル名の順で HC1・HC2）だけが RT +0.6 の別スポットに割れた CE。
SPLIT_FEATURES = (
    {"id": "30", "mz": 650.6000, "rt": 20.00, "values": [20, 21, 1000, 1010, 1020, 1030, 990, 995],
     "sme": [{"name": "CE 18:2", "inchikey": "RJECHNNFRHZQKU-RMUVNZEASA-N", "formula": "C45H76O2",
              "adduct": "[M+H]1+", "theo": 650.6000}]},
    {"id": "31", "mz": 650.6040, "rt": 20.60, "values": [1005, 1015, 22, 20, 21, 23, 20, 22],
     "sme": []},
    {"id": "32", "mz": 300.0000, "rt": 3.00, "values": [5, 6, 5, 5, 6, 5, 5, 5], "sme": []},
)


def test_tool_writes_output_and_returns_summary(tmp_path):
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=SPLIT_FEATURES)
    handle = _open({"mztab": mztab})
    target = tmp_path / "results" / "alignment_qc_pos.json"
    out = json.loads(R.dataset_alignment_qc(output_path=str(target), dataset_handle=handle))
    assert out["status"] == "success", out
    assert out["injection_order_source"] == "filename"
    assert out["summary"]["by_type"]["injection_block"] == 1
    assert out["summary"]["n_features_affected"] == 2
    assert out["output_path"] == str(target.resolve())
    assert "pairs" not in out                             # 全体はファイルにだけ書く
    written = json.loads(target.read_text(encoding="utf-8"))
    assert written["schema"] == "alignment-qc.v1"
    (pair,) = written["pairs"]
    assert pair["features"] == ["30", "31"] and pair["mergeable"] is True
    assert pair["block"] == {"spot": "31", "injections": [1, 2]}
    assert written["params"]["ppm"] == 15.0 and written["params"]["min_block_size"] == 2
    (merge,) = written["merge_groups"]
    assert merge["features"] == ["30", "31"] and merge["mergeable"] is True
    assert pair["merged"] is True and out["summary"]["n_merged_pairs"] == 1
    # 同じ入力なら同じファイル・出力先なしなら要約だけ・ハンドルの状態は変えない
    again = tmp_path / "again.json"
    R.dataset_alignment_qc(output_path=str(again), dataset_handle=handle)
    assert again.read_bytes() == target.read_bytes()
    inline = json.loads(R.dataset_alignment_qc(dataset_handle=handle))
    assert inline["status"] == "success" and "output_path" not in inline
    assert handles.registry.get(handle).dataset.pp_matrix is None
    bad = json.loads(R.dataset_alignment_qc(ppm=-1, dataset_handle=handle))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


def test_alignment_qc_without_dataset_is_missing_state():
    out = json.loads(R.dataset_alignment_qc())
    assert out["error"]["code"] == "missing_state"
