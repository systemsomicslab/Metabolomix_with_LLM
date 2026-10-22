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
    {"id": "26", "mz": 250.1000, "rt": 6.00, "values": [9, 9, 9, 1, 1, 1, 5, 5],
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
