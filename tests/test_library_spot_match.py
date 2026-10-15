"""`library_match_feature(spot_id=...)`: アラインメントのスポットを ID で照合する。

測定スペクトルは precursor m/z で `.dcl` を探さず、アラインメント `.dcl` の
`dcl[spot_id]` をそのまま使う（同じ m/z の別スポットを拾わないため）。
"""
import json
import textwrap

import pytest

from metabolomix.core import mcp_core, session_state
from metabolomix.library import tools
from metabolomix.msdial import alignment_spot
from tests.dcl_fixture import build_dcl_bytes

_MSP = textwrap.dedent("""\
    NAME: GABA
    PRECURSORMZ: 104.0706
    IONMODE: Positive
    Num Peaks: 2
    87.0441 999
    69.0335 500
""")

_RECORDS = [
    {"AlignmentID": 0, "Name": "PC 34:1", "RT": 3.0, "MassCenter": 760.5851,
     "IonMode": "Positive"},
    {"AlignmentID": 1, "Name": "Unknown", "RT": 1.2, "MassCenter": 104.0706,
     "IonMode": "Positive"},
    {"AlignmentID": 2, "Name": "Unknown", "RT": 2.0, "MassCenter": 300.0,
     "IonMode": "Positive"},
]


@pytest.fixture(autouse=True)
def workspace(tmp_path, monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    monkeypatch.setenv("LIPIDMIX_LIBRARY_CACHE_DIR", str(tmp_path / "cache"))
    (tmp_path / "lib.msp").write_text(_MSP, encoding="utf-8")
    (tmp_path / "AlignmentResult_1.arf2").write_bytes(b"placeholder")
    (tmp_path / "AlignmentResult_1.dcl").write_bytes(build_dcl_bytes([
        # 0: 同じ precursor の別スポットとして、m/z 検索なら拾われうる
        {"precursor_mz": 760.5851, "rt": 3.0, "spectrum": [(184.07, 999.0)]},
        {"precursor_mz": 104.0706, "rt": 1.2,
         "spectrum": [(87.0441, 999.0), (69.0335, 500.0)]},
        {"precursor_mz": 999.0, "rt": 2.0, "spectrum": [(100.0, 10.0)]},
    ]))
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: _RECORDS)
    tools.library_load()
    return tmp_path


def _match(**kwargs):
    return json.loads(tools.library_match_feature(**kwargs))


def test_spot_id_takes_precursor_rt_and_polarity_from_arf2():
    payload = _match(spot_id=1)
    assert payload["status"] == "success"
    query = payload["query"]
    assert query["spot_id"] == 1
    assert query["precursor_mz"] == pytest.approx(104.0706)
    assert query["rt"] == pytest.approx(1.2)
    assert query["ion_mode"] == "positive"
    assert "GABA" in payload["candidates_table"]
    assert session_state.session.library.last_match["spot_id"] == 1
    assert session_state.session.library.last_match["measured"] == [
        [87.0441, 999.0], [69.0335, 500.0]]


def test_spot_id_uses_the_spot_entry_not_a_precursor_search(monkeypatch):
    def fail(*_a, **_k):
        raise AssertionError("spot_id 経由で precursor 検索を使ってはいけない")

    monkeypatch.setattr(tools, "_measured_spectrum", fail)
    assert _match(spot_id=1)["status"] == "success"


def test_precursor_mismatch_refuses_to_score():
    payload = _match(spot_id=2)
    assert payload["status"] == "dcl_precursor_mismatch"
    assert session_state.session.library.last_match is None


def test_unknown_spot_is_reported():
    payload = _match(spot_id=42)
    assert payload["status"] == "spot_not_found"
    assert "42" in payload["message"]


def test_exactly_one_of_precursor_or_spot_is_required():
    assert _match()["status"] == "invalid_arguments"
    assert _match(precursor_mz=104.07, spot_id=1)["status"] == "invalid_arguments"


def test_spot_without_msms_is_not_found(workspace):
    (workspace / "AlignmentResult_1.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 3.0, "spectrum": []},
        {"precursor_mz": 104.0706, "rt": 1.2, "spectrum": []},
    ]))
    payload = _match(spot_id=1)
    assert payload["status"] == "not_found"
    assert "未取得" in payload["message"]


def test_mirror_plot_follows_a_spot_match():
    _match(spot_id=1)
    text = tools.library_plot_mirror(output="payload")
    assert "GABA" in text
