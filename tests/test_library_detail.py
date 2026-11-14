"""照合の詳細（前駆体を除いた断片の本数）・ライブラリの検索・spot の代表スペクトル
（spec 2026-10-22 §4 の 2）。"""
from __future__ import annotations

import json
import sqlite3
import textwrap

import pytest

from metabolomix.analysis.spectral_match import fragment_matches_excluding_precursor
from metabolomix.core import session_state
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
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
    (tmp_path / "lib.msp").write_text(_MSP, encoding="utf-8")
    (tmp_path / "AlignmentResult_1.arf2").write_bytes(b"placeholder")
    (tmp_path / "AlignmentResult_1.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 3.0,
         "spectrum": [(184.0733, 999.0), (760.5851, 400.0)]},
        {"precursor_mz": 104.0706, "rt": 1.2, "spectrum": [(87.0441, 999.0), (69.0335, 500.0)]},
    ]))
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: _RECORDS)
    assert json.loads(tools.library_load(file_path=str(tmp_path / "lib.msp")))["status"] == "success"
    _FILES.update(arf2_path=str(tmp_path / "AlignmentResult_1.arf2"),
                  dcl_file=str(tmp_path / "AlignmentResult_1.dcl"))
    return tmp_path


#: workspace が書いた `.arf2` / `.dcl`（データセットを開かずに、明示のパスで照合する）。
_FILES: dict = {}


def _match(**kwargs):
    for key, value in _FILES.items():
        kwargs.setdefault(key, value)
    return json.loads(tools.library_match_feature(**kwargs))


def test_fragment_count_excludes_precursor():
    alignment = [{"mz": 184.07, "measured": 1, "reference": 1, "matched": True},
                 {"mz": 500.0, "measured": 0, "reference": 1, "matched": False},
                 {"mz": 760.585, "measured": 1, "reference": 1, "matched": True}]
    assert fragment_matches_excluding_precursor(alignment, 760.5851, 0.01) == 1
    assert fragment_matches_excluding_precursor([], 760.5851, 0.01) is None


def test_fragment_count_floor_is_relative_to_the_base_peak():
    # measured は基準ピーク 100 の強度。下限 0.01 は 1.0 未満の窓を数えない
    alignment = [{"mz": 184.07, "measured": 100.0, "reference": 1, "matched": True},
                 {"mz": 300.0, "measured": 0.4, "reference": 1, "matched": True},
                 {"mz": 760.585, "measured": 40.0, "reference": 1, "matched": True}]
    assert fragment_matches_excluding_precursor(alignment, 760.5851, 0.01) == 2
    assert fragment_matches_excluding_precursor(alignment, 760.5851, 0.01,
                                                min_relative_intensity=0.01) == 1


def test_detail_false_output_has_no_new_keys():
    out = _match(spot_id=0)
    assert "alignments" not in out and "measured_peaks" not in out
    assert out["candidates_table"].splitlines()[0].endswith("entropy_similarity")


def test_default_table_has_inchikey_and_its_standard_flag():
    # InChIKey は既定の列（同じ名前の立体異性体を見分ける。事実確認 B6）
    lines = _match(spot_id=0)["candidates_table"].splitlines()
    header = lines[0].split("\t")
    assert header.index("inchikey") == header.index("formula") + 1
    assert header.index("inchikey_standard") == header.index("inchikey") + 1
    rows = {line.split("\t")[1]: dict(zip(header, line.split("\t"))) for line in lines[1:]}
    assert rows["PC 34:1"]["inchikey"] == "JLPULHDHAOZNQI-ZTIMHPMXSA-N"
    assert rows["PC 34:1"]["inchikey_standard"] == "true"


def test_inchikey_standard_reads_the_second_block_flag():
    from metabolomix.library.record import inchikey_standard

    assert inchikey_standard("JLPULHDHAOZNQI-ZTIMHPMXSA-N") is True
    assert inchikey_standard("JLPULHDHAOZNQI-ZTIMHPMXNA-N") is False
    assert inchikey_standard("JLPULHDHAOZNQI") is None
    assert inchikey_standard(None) is None


def test_detail_true_adds_columns_and_alignments():
    out = _match(spot_id=0, detail=True)
    lines = out["candidates_table"].splitlines()
    assert lines[0].endswith("entropy_similarity\tfragment_matches_excl_precursor")
    by_name = {line.split("\t")[1]: line.split("\t") for line in lines[1:]}
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
    assert [line.rsplit("\t", 1)[0] for line in detail] == plain


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
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=1, **_FILES))
    assert out["status"] == "ok" and out["n_peaks"] == 2
    assert out["peaks"][0] == [87.0441, 999.0]
    assert out["precursor_mz"] == pytest.approx(104.0706)
    missing = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=9, **_FILES))
    assert missing["status"] == "spot_not_found"


# ---------- MS-DIAL の注釈との食い違い（E6）と spot_id の対向プロット（C16） ----------

def test_spot_match_reports_msdial_agreement():
    out = _match(spot_id=0)
    assert out["msdial_name"] == "PC 34:1" and out["msdial_rank"] == 1
    assert out["top_differs_from_msdial"] is False
    unknown = _match(spot_id=1)                     # 未注釈（Unknown）のスポットには足さない
    assert "msdial_name" not in unknown and "top_differs_from_msdial" not in unknown


def test_msdial_agreement_compares_inchikey_species_then_name():
    from metabolomix.library.matching import msdial_agreement

    scored = [({"name": "PC 16:0_18:1", "inchikey": "WTJKGGKOPKCXLL-VYOBOKEXSA-N"}, {}),
              ({"name": "PC 34:1", "inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N"}, {})]
    by_key = msdial_agreement(scored, {"name": "PC 34:1",
                                       "inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N"})
    assert by_key == {"msdial_name": "PC 34:1", "msdial_rank": 2, "top_differs_from_msdial": True}
    # InChIKey が片方に無ければ脂質の種で比べる（同じ種の別の分子種も同じ物質に数える）
    by_species = msdial_agreement(scored, {"name": "PC 18:1_16:0", "inchikey": None})
    assert by_species["msdial_rank"] == 1 and by_species["top_differs_from_msdial"] is False
    by_name = msdial_agreement([({"name": "GABA"}, {})], {"name": "w/o MS2: gaba"})
    assert by_name["msdial_rank"] == 1
    assert msdial_agreement(scored, {"name": "Unknown"}) == {}
    assert msdial_agreement(scored, None) == {}
    assert msdial_agreement([], {"name": "PC 34:1"}) == {}


def test_mirror_by_spot_id_on_a_dataset(monkeypatch):
    from metabolomix.tools import _holder

    real = _holder.alignment_paths

    def explicit(holder, dataset_handle, arf2_path, dcl_file, **kwargs):
        return real(holder, dataset_handle, arf2_path or _FILES["arf2_path"],
                    dcl_file or _FILES["dcl_file"], **kwargs)

    monkeypatch.setattr(_holder, "alignment_paths", explicit)
    _match(precursor_mz=104.0706, ion_mode="positive")       # 別の照合が直前にある
    payload = json.loads(tools.library_plot_mirror(output="payload", spot_id=0))
    assert payload["match"]["spot_id"] == 0
    assert payload["match"]["precursor_mz"] == pytest.approx(760.5851)
    assert payload["match"]["name"] == "PC 34:1"
    isomer = json.loads(tools.library_plot_mirror(output="payload", spot_id=0,
                                                  candidate="WTJKGGKOPKCXLL"))
    assert isomer["match"]["name"] == "PC 16:0_18:1 isomer" and isomer["match"]["rank"] == 2
    caption = tools.library_plot_mirror(spot_id=0)[0]
    assert caption.startswith("Mirror: spot 0 (precursor m/z=760.5851)")
    assert "InChIKey=JLPULHDHAOZNQI-ZTIMHPMXSA-N" in caption
    missing = json.loads(tools.library_plot_mirror(output="payload", spot_id=99))
    assert missing["status"] == "spot_not_found"
    assert session_state.session.library.last_match is None   # 古い照合を黙って描かない


# ---------- 対向プロットの候補の選び方（全体の InChIKey は立体異性体まで見分ける） ----------

def test_pick_candidate_prefers_the_full_inchikey_over_its_first_block():
    candidates = [{"name": "D-Glucose", "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N"},
                  {"name": "D-Galactose", "inchikey": "WQZGKKKJIJFFOK-SVZMEOIVSA-N"}]
    # 同じ骨格（先頭ブロック）の 2 位を全体の InChIKey で選ぶ（大小を問わない）
    assert tools._pick_candidate(candidates, "WQZGKKKJIJFFOK-SVZMEOIVSA-N") == (1, "inchikey")
    assert tools._pick_candidate(candidates, " wqzgkkkjijffok-svzmeoivsa-n ") == (1, "inchikey")
    # 先頭 14 文字だけなら、先頭ブロックが同じ最上位
    assert tools._pick_candidate(candidates, "WQZGKKKJIJFFOK") == (0, "inchikey_block")
    # 全体で合う候補が無いときだけ先頭ブロックで選ぶ（選び方で分かる）
    assert tools._pick_candidate(candidates, "WQZGKKKJIJFFOK-UHFFFAOYSA-N") == (0, "inchikey_block")
    assert tools._pick_candidate(candidates, "d-galactose") == (1, "name")
    assert tools._pick_candidate(candidates, "BTCSSZJGUNDROE-UHFFFAOYSA-N") is None
    assert tools._pick_candidate(candidates, "Fructose") is None


def test_mirror_full_inchikey_draws_the_stereoisomer_asked_for():
    _match(spot_id=0)
    candidates = session_state.session.library.last_match["candidates"]
    # 2 位の候補を 1 位と同じ骨格の別の立体異性体にする（合成の InChIKey）
    candidates[1]["inchikey"] = "JLPULHDHAOZNQI-UHFFFAOYSA-N"
    exact = json.loads(tools.library_plot_mirror(output="payload",
                                                 candidate="JLPULHDHAOZNQI-UHFFFAOYSA-N"))
    assert exact["match"]["rank"] == 2 and exact["match"]["picked_by"] == "inchikey"
    assert exact["match"]["inchikey"] == "JLPULHDHAOZNQI-UHFFFAOYSA-N"
    by_rank = json.loads(tools.library_plot_mirror(output="payload"))
    assert by_rank["match"]["rank"] == 1 and by_rank["match"]["picked_by"] == "rank"
    # 全体で合う候補が無ければ先頭ブロックの最上位を描き、caption で知らせる
    fallback = json.loads(tools.library_plot_mirror(output="payload",
                                                    candidate="JLPULHDHAOZNQI-XXXXXXXXXX-N"))
    assert fallback["match"]["rank"] == 1 and fallback["match"]["picked_by"] == "inchikey_block"
    caption = tools.library_plot_mirror(candidate="JLPULHDHAOZNQI-XXXXXXXXXX-N")[0]
    assert "No candidate has the full InChIKey given" in caption
    assert "No candidate has the full InChIKey" not in tools.library_plot_mirror(
        candidate="JLPULHDHAOZNQI")[0]


# ---------- dcl_spot_spectrum の極性（照合の絞り込みとは別に .arf2 の値を残す） ----------

@pytest.mark.parametrize("ion_mode, expected", [
    ("Positive", "positive"), ("Negative", "negative"), ("Both", "both"), ("Unknown", "unknown")])
def test_spot_spectrum_reports_the_arf2_polarity(monkeypatch, ion_mode, expected):
    records = [dict(_RECORDS[0]), dict(_RECORDS[1], IonMode=ion_mode)]
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: records)
    out = json.loads(dcl_tools.dcl_spot_spectrum(spot_id=1, **_FILES))
    assert out["status"] == "ok" and out["ion_mode"] == expected
    spot = tools.spot_msms_query(1, rt=None, ion_mode=None, **_FILES)
    # 照合の絞り込みは Both / Unknown を None（絞らない）にする
    assert spot["ion_mode"] == (expected if expected in ("positive", "negative") else None)
    assert spot["arf2_ion_mode"] == expected


# ---------- detail の断片の本数は合否の規則と同じ下限で数える ----------

def _weak_fragment_dcl(workspace):
    # spot 0 の 184.0733 を基準ピーク（760.5851）の 0.5 % にする（既定の下限 1 % 未満）
    (workspace / "AlignmentResult_1.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": 760.5851, "rt": 3.0,
         "spectrum": [(184.0733, 5.0), (760.5851, 1000.0)]},
        {"precursor_mz": 104.0706, "rt": 1.2, "spectrum": [(87.0441, 999.0), (69.0335, 500.0)]},
    ]))


def _fragment_cells(out):
    lines = out["candidates_table"].splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return {row["name"]: row["fragment_matches_excl_precursor"] for row in rows}


def test_detail_fragment_count_uses_the_rule_floor_by_default(workspace):
    _weak_fragment_dcl(workspace)
    out = _match(spot_id=0, detail=True)
    assert out["min_fragment_relative_intensity"] == 0.01
    assert _fragment_cells(out)["PC 34:1"] == "0"           # 0.5 % の 184 は数えない
    unfloored = _match(spot_id=0, detail=True, min_fragment_relative_intensity=0)
    assert unfloored["min_fragment_relative_intensity"] == 0
    assert _fragment_cells(unfloored)["PC 34:1"] == "1"
    # 同じ名前の列が、既定のまま library_msms_qualify と同じ本数になる
    from metabolomix.tools import reanalysis_tools

    qualify = json.loads(reanalysis_tools.library_msms_qualify(
        {"target": {"inchikey": "JLPULHDHAOZNQI-ZTIMHPMXSA-N"}}, spot_id=0, **_FILES))
    assert _fragment_cells({"candidates_table": qualify["candidates"]})["PC 34:1"] == "0"


@pytest.mark.parametrize("bad", [-0.1, 1.5, float("nan"), True, "0.01"])
def test_detail_floor_rejects_values_outside_zero_to_one(bad):
    out = _match(spot_id=0, detail=True, min_fragment_relative_intensity=bad)
    assert out["status"] == "invalid_arguments"
    assert "min_fragment_relative_intensity" in out["message"]
