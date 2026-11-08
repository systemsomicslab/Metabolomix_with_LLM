"""脂質の種の目標（`lipid_species`）を MCP のツールとして通す（目標照合・証拠の束・MS/MS の合否）。

MTBLS3233 の主張 C1〜C4 型（「PC 18:3_16:1」。InChIKey が一意に決まらない）。fixture はテスト
自身が作る（合成の mzTab-M・.msp・.dcl）。
"""
from __future__ import annotations

import json
import textwrap

import pytest

from metabolomix.analysis.feature_search import theoretical_mz
from metabolomix.core import handles, session_state
from metabolomix.library import tools as library_tools
from metabolomix.msdial import alignment_spot
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools
from metabolomix.tools import reanalysis_tools as R
from tests.dcl_fixture import build_dcl_bytes
from tests.reanalysis_fixtures import write_mztab
from tests.config_fixtures import use_config

PC_FORMULA = "C42H76NO8P"
PC_MH = round(theoretical_mz(PC_FORMULA, "[M+H]+"), 4)
BASE = [10, 12, 11, 20, 22, 21, 15, 16]

FEATURES = (
    {"id": "30", "mz": PC_MH, "rt": 8.40, "values": BASE,
     "sme": [{"name": "PC 16:1_18:3", "inchikey": "null", "formula": PC_FORMULA,
              "adduct": "[M+H]1+", "theo": PC_MH}]},
    {"id": "31", "mz": PC_MH, "rt": 8.45, "values": BASE,
     "sme": [{"name": "PE 37:4", "inchikey": "null", "formula": PC_FORMULA,
              "adduct": "[M+H]1+", "theo": PC_MH}]},
    {"id": "32", "mz": 300.0, "rt": 2.0, "values": BASE, "sme": []},
)

_MSP = textwrap.dedent(f"""\
    NAME: PE 37:4
    PRECURSORMZ: {PC_MH}
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: {PC_FORMULA}
    Num Peaks: 2
    613.5 999
    {PC_MH} 100

    NAME: PC 16:1_18:3
    PRECURSORMZ: {PC_MH}
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: {PC_FORMULA}
    Num Peaks: 3
    184.0733 999
    478.33 200
    {PC_MH} 100
""")

CLAIM = {"claim_id": "C1",
         "target": {"lipid_species": "PC 18:3_16:1", "formula": PC_FORMULA, "rt_min": 8.4},
         "scan": {"ppm": [10, 20], "rt_tol_min": [0.3]}}


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


@pytest.fixture
def mztab(tmp_path):
    return write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=FEATURES)


@pytest.fixture
def library(tmp_path, monkeypatch, mztab):
    msp = tmp_path / "lab" / "pos.msp"
    msp.parent.mkdir()
    msp.write_text(_MSP, encoding="ascii")
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    run = mztab.parent
    (run / "AlignResult-1.arf2").write_bytes(b"placeholder")
    entries = [{"precursor_mz": 999.0, "rt": 9.0, "spectrum": []} for _ in range(33)]
    entries[30] = {"precursor_mz": PC_MH, "rt": 8.4,
                   "spectrum": [(184.0733, 999.0), (478.33, 150.0), (PC_MH, 300.0)]}
    (run / "AlignResult-1.dcl").write_bytes(build_dcl_bytes(entries))
    records = [{"AlignmentID": int(f["id"]), "MassCenter": f["mz"], "RT": f["rt"],
                "IonMode": "Positive", "Name": "x"} for f in FEATURES]
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: records)
    return msp


def _open(mztab, msp=None):
    out = json.loads(handle_tools.load_dataset(mztab_path=str(mztab), library="none"))
    assert out["status"] == "success", out
    if msp is not None:
        # mzTab-M の記憶は `unknown` なので、library_load には file_path で明示する（mzTab-M の本番の経路）。
        loaded = json.loads(library_tools.library_load(file_path=str(msp),
                                                       dataset_handle=out["dataset_handle"]))
        assert loaded["status"] == "success", loaded
    return out["dataset_handle"]


def test_feature_search_resolves_the_species_and_contradicts_the_other_class(mztab):
    handle = _open(mztab)
    out = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", "lipid_species": "PC 18:3_16:1", "formula": PC_FORMULA,
                  "rt_min": 8.4}], dataset_handle=handle))
    [target] = out["targets"]
    assert target["status"] == "resolved" and target["selected_feature_id"] == "30"
    assert target["selected_evidence"] == "sme_match"
    header, *rows = out["candidates"].splitlines()
    cols = header.split("\t")
    by_feature = {r.split("\t")[cols.index("feature_id")]: dict(zip(cols, r.split("\t")))
                  for r in rows if r.split("\t")[cols.index("adduct")] == "[M+H]+"}
    assert by_feature["31"]["evidence"] == "contradicted"
    assert by_feature["31"]["reasons"] == "identifier_mismatch"


def test_feature_search_species_target_without_formula_is_invalid(mztab):
    handle = _open(mztab)
    out = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", "lipid_species": "PC"}], dataset_handle=handle))
    assert out["error"]["code"] == "FEATURE_TARGET_INVALID"
    assert out["error"]["details"]["field"] == "lipid_species"


def test_claim_evidence_for_a_species_claim(mztab):
    handle = _open(mztab)
    out = json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=handle))
    assert out["status"] == "success"
    assert out["default_selection"]["status"] == "resolved"
    assert out["default_selection"]["selected_feature_id"] == "30"
    assert out["rule"]["target"]["lipid_species"] == "PC 34:4"
    assert out["rule"]["target"]["lipid_species_input"] == "PC 18:3_16:1"
    by_id = {c["feature_id"]: c for c in out["candidates"]}
    assert by_id["30"]["evidence"] == ["sme_match"]
    assert "contradicted" in by_id["31"]["evidence"]
    same = json.loads(R.dataset_claim_evidence(
        {**CLAIM, "target": {**CLAIM["target"], "lipid_species": "PC 34:4"}},
        dataset_handle=handle))
    assert same["rule_hash"] == out["rule_hash"]


def test_claim_evidence_species_msms_rules(mztab, library):
    handle = _open(mztab, msp=library)
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {
        "lipid_class": "PC", "min_fragments_excl_precursor": 2,
        "diagnostic_ions": [{"mz": 184.0733, "label": "phosphocholine"}]}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["rule"]["msms"]["rules"]["target"] == {"lipid_species": "PC 34:4",
                                                      "formula": PC_FORMULA}
    msms = {c["feature_id"]: c["msms"] for c in out["candidates"]}["30"]
    assert "qualify" in msms, msms
    qualify = msms["qualify"]
    assert qualify["verdict"] == "qualified" and qualify["target_rank"] == 1, qualify


def test_msms_qualify_tool_takes_a_species_target(mztab, library):
    handle = _open(mztab, msp=library)
    out = json.loads(R.library_msms_qualify(
        {"target": {"lipid_species": "PC 18:3_16:1"}, "min_fragments_excl_precursor": 2},
        spot_id=30, dataset_handle=handle))
    assert out["verdict"] == "qualified", out
    assert out["rules"]["target"] == {"lipid_species": "PC 34:4", "formula": None}
