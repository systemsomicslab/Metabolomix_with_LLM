"""照合の詳細（前駆体を除いた断片の本数）・ライブラリの検索・spot の代表スペクトル
（spec 2026-10-22 §4 の 2）。"""
from __future__ import annotations

import json
import sqlite3
import textwrap

import pytest

from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor
from metabolomix.core import mcp_core, session_state
from metabolomix.dcl import tools as dcl_tools
from metabolomix.library import tools
from metabolomix.library.store import open_store
from metabolomix.msdial import alignment_spot
from tests.dcl_fixture import build_dcl_bytes
from tests.config_fixtures import use_config

_MSP = textwrap.dedent("""\
    NAME: PC 34:1
    PRECURSORMZ: 760.5851
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: C42H82NO8P
    INCHIKEY: JLPULHDHAOZNQI-ZTIMHPMXSA-N
    Num Peaks: 2
    184.0733 999
    760.5851 300

    NAME: PC 16:0_18:1 isomer
    PRECURSORMZ: 760.5851
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: C42H82NO8P
    INCHIKEY: WTJKGGKOPKCXLL-VYOBOKEXSA-N
    Num Peaks: 1
    760.5851 999

    NAME: GABA
    PRECURSORMZ: 104.0706
    PRECURSORTYPE: [M+H]+
    IONMODE: Positive
    FORMULA: C4H9NO2
    INCHIKEY: BTCSSZJGUNDROE-UHFFFAOYSA-N
    Num Peaks: 2
    87.0441 999
    69.0335 500
""")

_RECORDS = [
    {"AlignmentID": 0, "Name": "PC 34:1", "RT": 3.0, "MassCenter": 760.5851,
     "IonMode": "Positive"},
    {"AlignmentID": 1, "Name": "Unknown", "RT": 1.2, "MassCenter": 104.0706,
     "IonMode": "Positive"},
]


@pytest.fixture(autouse=True)
def workspace(tmp_path, monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    (tmp_path / "lib.msp").write_text(_MSP, encoding="utf-8")
    (tmp_path / "AlignmentResult_1.arf2").write_bytes(b"placeholder")
    (tmp_path / "AlignmentResult_1.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 3.0,
         "spectrum": [(184.0733, 999.0), (760.5851, 400.0)]},
        {"precursor_mz": 104.0706, "rt": 1.2, "spectrum": [(87.0441, 999.0), (69.0335, 500.0)]},
    ]))
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: _RECORDS)
    tools.library_load()
    return tmp_path


def _match(**kwargs):
    return json.loads(tools.library_match_feature(**kwargs))


def test_fragment_count_excludes_precursor():
    alignment = [{"mz": 184.07, "measured": 1, "reference": 1, "matched": True},
                 {"mz": 500.0, "measured": 0, "reference": 1, "matched": False},
                 {"mz": 760.585, "measured": 1, "reference": 1, "matched": True}]
    assert fragment_matches_excluding_precursor(alignment, 760.5851, 0.01) == 1
    assert fragment_matches_excluding_precursor([], 760.5851, 0.01) is None


def test_detail_false_output_has_no_new_keys():
    out = _match(spot_id=0)
    assert "alignments" not in out and "measured_peaks" not in out
    assert out["candidates_table"].splitlines()[0].endswith("entropy_similarity")


def test_detail_true_adds_columns_and_alignments():
    out = _match(spot_id=0, detail=True)
    lines = out["candidates_table"].splitlines()
    assert lines[0].endswith("inchikey\tfragment_matches_excl_precursor")
    by_name = {line.split("\t")[1]: line.split("\t") for line in lines[1:]}
    assert by_name["PC 34:1"][-2] == "JLPULHDHAOZNQI-ZTIMHPMXSA-N"
    assert by_name["PC 34:1"][-1] == "1"                       # 184 だけ（前駆体を除く）
    assert by_name["PC 16:0_18:1 isomer"][-1] == "0"           # 前駆体だけの一致
    assert out["precursor_exclusion_tol"] == out["query"]["ms2_tol"]
    assert out["measured_peaks"][0][0] == pytest.approx(184.0733)
    ranks = [a["rank"] for a in out["alignments"]]
    assert ranks == list(range(1, len(ranks) + 1))
    assert all(len(pair) == 3 for a in out["alignments"] for pair in a["matched"])


def test_detail_does_not_change_existing_columns():
    plain = _match(spot_id=0)["candidates_table"].splitlines()
    detail = _match(spot_id=0, detail=True)["candidates_table"].splitlines()
    assert [line.rsplit("\t", 2)[0] for line in detail] == plain


def test_search_by_name_formula_inchikey():
    def search(**kwargs):
        return json.loads(tools.library_search(**kwargs))

    by_name = search(name="pc 34")
    assert by_name["n_hits"] == 1 and "PC 34:1" in by_name["records"]
    by_formula = search(formula="C42H82NO8P")
    assert by_formula["n_hits"] == 2
    by_key = search(inchikey="BTCSSZJGUNDROE-XXXXXXXXXX-N")       # 先頭ブロックで照合
    assert by_key["n_hits"] == 1 and "GABA" in by_key["records"]
    assert search()["status"] == "invalid_arguments"
    limited = search(formula="C42H82NO8P", limit=1)
    assert limited["n_hits"] == 1 and limited["truncated"] == 1


def test_search_without_library_is_missing_state():
    session_state.session = session_state.AnalysisSession()
    out = json.loads(tools.library_search(name="x"))
    assert out["error"]["code"] == "missing_state"


def test_read_only_store_rejects_writes(workspace):
    store = open_store(workspace / "lib.msp", read_only=True)
    with pytest.raises(sqlite3.OperationalError):
        store._conn.execute("DELETE FROM record")
    store.close()


def test_spot_spectrum_returns_peaks(workspace):
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=1))
    assert out["status"] == "ok" and out["n_peaks"] == 2
    assert out["peaks"][0] == [87.0441, 999.0]
    assert out["precursor_mz"] == pytest.approx(104.0706)
    missing = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=9))
    assert missing["status"] == "spot_not_found"
