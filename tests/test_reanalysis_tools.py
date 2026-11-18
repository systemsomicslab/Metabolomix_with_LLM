"""再解析のツール（spec 2026-10-22 §4 の 3〜11）を MCP のツールとして呼ぶ。

ハンドルの上でも対話のセッションでも動き、どちらの状態も書き換えないこと・同じ spec なら
同じ JSON になることを確かめる。fixture はすべてテスト自身が作る。
"""
from __future__ import annotations

import json
import textwrap

import pytest

from metabolomix.core import handles, session_state
from metabolomix.dcl import tools as dcl_tools
from metabolomix.library import tools as library_tools
from metabolomix.msdial import alignment_spot
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools, mztab_tools
from metabolomix.tools import reanalysis_tools as R
from tests.dcl_fixture import build_dcl_bytes
from tests.reanalysis_fixtures import SAMPLES, write_mztab
from tests.config_fixtures import use_config
from tests.engine_chain import declare_groups, engine_compare, hc_ra

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
    use_config(monkeypatch, tmp_path, {"paths.library_cache": tmp_path / "cache"})
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


def _open(data, msp=None):
    kwargs = {"library": "none"}
    # 旧 `dataset_handle_open` は直近を動かさなかった。この検査は「ハンドルを渡す呼び出しは対話の状態に
    # 触れない」ことを縛るので、開いたあとで直近を元に戻す（`load_dataset` はいつも直近にする）。
    previous = session_state.session.current_handle_id
    out = json.loads(handle_tools.load_dataset(mztab_path=str(data["mztab"]), **kwargs))
    assert out["status"] == "success", out
    session_state.session.current_handle_id = previous
    if msp is not None:
        # mzTab-M の記憶は `unknown` なので、library_load には file_path で明示する（mzTab-M の本番の経路）。
        loaded = json.loads(library_tools.library_load(file_path=str(msp),
                                                       dataset_handle=out["dataset_handle"]))
        assert loaded["status"] == "success", loaded
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

def _compare_on(handle) -> dict:
    """ハンドルのデータセットで HC 対 RA を統計エンジンで比較する（群は MS-DIAL の試料の宣言で与える）。"""
    declare_groups(handles.registry.get(handle).dataset, hc_ra)
    return engine_compare("HC", "RA", dataset_handle=handle)


def _analysis_state(handle):
    ds = handles.registry.get(handle).dataset
    return set(ds.results), set(ds.stat_matrices)


def test_analysis_run_on_handle_returns_rows_and_leaves_states(data):
    handle = _open(data)
    before = _analysis_state(handle)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    assert out["status"] == "success" and out["comparison"] == "RA_vs_HC"
    rows = out["analyses"][0]["rows"].splitlines()
    assert rows[0].startswith("feature_id\tstatus\treason\tgroup_ref\tgroup_test")
    assert rows[1].split("\t")[:5] == ["20", "tested", "", "HC", "RA"]
    assert session_state.session.dataset is None
    assert _analysis_state(handle) == before


def test_analysis_run_reports_merge_split_state_and_merged_from_column(data):
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    assert out["merge_split"]["enabled"] is True and out["merge_split"]["source"] == "find_split_pairs"
    assert "merged_from" in out["analyses"][0]["rows"].splitlines()[0].split("\t")
    off = json.loads(R.dataset_analysis_run(
        _spec(data, features={"feature_ids": ["20", "26"], "merge_split": False}),
        dataset_handle=handle))
    assert off["merge_split"]["enabled"] is False and off["rule_hash"] != out["rule_hash"]


def test_analysis_run_rows_have_detection_columns(data):
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    analysis = out["analyses"][0]
    header = analysis["rows"].splitlines()[0].split("	")
    assert header[-5:] == ["detected_reference", "detected_test", "detection_fisher_p",
                           "detection_driven", "constant_in_group"]
    assert analysis["detection_available"] is False       # 合成の mzTab は .arf が無い
    assert analysis["n_detection_driven"] == 0


def test_analysis_run_on_session_dataset(data):
    handle_tools.load_dataset(mztab_path=str(data["mztab"]))
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


def test_analysis_run_many_rows_are_truncated_not_refused(tmp_path):
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
    item = out["analyses"][0]
    assert "error" not in out and item["n_rows"] == 51 and item["truncated"] is True
    assert item["rows_order"] == "p_ascending" and any("output_path" in n for n in out["notes"])
    table = item["rows"].splitlines()
    header = table[0].split("\t")
    shown = [dict(zip(header, line.split("\t"))) for line in table[1:]]
    assert len(shown) == R.MAX_INLINE_ROWS
    ps = [float(r["p"]) for r in shown]
    assert ps == sorted(ps)                                   # p の小さい順の先頭 50 行
    ok = json.loads(R.dataset_analysis_run(spec, output_path=str(tmp_path / "o.json")))
    assert ok["analyses"][0]["n_rows"] == 51 and "truncated" not in ok["analyses"][0]
    assert "notes" not in ok
    saved = json.loads((tmp_path / "o.json").read_text(encoding="utf-8"))
    assert len(saved["analyses"][0]["rows"]) == 51


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
    # 勧める許容幅は 10 件から。2 件では出さず、理由を添える（A3-T12）
    assert low["recommended_mz_tol_ppm"] is None and high["recommended_mz_tol_ppm"] is None
    assert out["recommended_tolerance"]["mz_tol_ppm"] is None
    assert out["recommended_tolerance"]["reason"] == "too_few"


def test_mass_accuracy_recommends_a_tolerance(tmp_path):
    # ppm が -4.0 … +1.5（0.5 刻み）の 12 個の同定 → 中央値 -1.25、MAD 1.5 → ceil(1.25 + 4.5) = 6
    features = tuple(
        {"id": str(i + 1), "mz": 300.0 + i, "rt": 1.0 + i, "values": [10] * 8, "matched": True,
         "sme": [{"name": f"N{i}", "inchikey": f"{chr(65 + i) * 14}-UHFFFAOYSA-N", "formula": "C1",
                  "adduct": "[M+H]1+", "theo": 300.0 + i,
                  "exp": (300.0 + i) * (1 + (-4.0 + 0.5 * i) * 1e-6)}]}
        for i in range(12))
    handle = _open({"mztab": write_mztab(tmp_path / "t" / "AlignResult-1.mzTab", features=features)})
    out = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))
    rec = out["recommended_tolerance"]
    assert rec["n"] == 12 and rec["mz_tol_ppm"] == 6.0
    assert rec["median_ppm"] == pytest.approx(-1.25, abs=1e-3)
    assert rec["mad_ppm"] == pytest.approx(1.5, abs=1e-3)
    low, high = out["bands"]
    assert low["n"] == 0 and low["recommended_mz_tol_ppm"] is None
    assert high["n"] == 12 and high["recommended_mz_tol_ppm"] == 6.0


#: SML の `opt_global_spectrum_matched` を持つ合成の特徴（ppm はそれぞれ +2 / +8 / -3、1 つは同定なし）。
MATCH_FEATURES = (
    {"id": "1", "mz": 200.0, "rt": 1.0, "values": [10] * 8, "matched": True,
     "sme": [{"name": "A", "inchikey": "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "formula": "C1", "adduct": "[M+H]1+",
              "theo": 200.0, "exp": 200.0004}]},
    {"id": "2", "mz": 300.0, "rt": 2.0, "values": [10] * 8, "matched": False,
     "sme": [{"name": "B", "inchikey": "BBBBBBBBBBBBBB-UHFFFAOYSA-N", "formula": "C2", "adduct": "[M+H]1+",
              "theo": 300.0, "exp": 300.0024}]},
    {"id": "3", "mz": 400.0, "rt": 3.0, "values": [10] * 8, "matched": True,
     "sme": [{"name": "C", "inchikey": "CCCCCCCCCCCCCC-UHFFFAOYSA-N", "formula": "C3", "adduct": "[M+H]1+",
              "theo": 400.0, "exp": 399.9988}]},
    {"id": "4", "mz": 500.0, "rt": 4.0, "values": [10] * 8, "matched": True, "sme": []},
)


@pytest.fixture
def match_data(tmp_path):
    return {"mztab": write_mztab(tmp_path / "m" / "AlignResult-1.mzTab", features=MATCH_FEATURES)}


def test_mass_accuracy_msms_matched_subset(match_data):
    handle = _open(match_data)
    full = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))
    assert full["subset"] == "rank1" and full["overall"]["n"] == 3
    matched = json.loads(R.dataset_mass_accuracy(subset="msms_matched", dataset_handle=handle))
    assert matched["subset"] == "msms_matched"
    # 特徴 1 と 3（SML の spectrum_matched が真）だけ。特徴 2（偽）と 4（SME なし）は入らない
    assert matched["overall"]["n"] == 2
    assert matched["overall"]["median"] == pytest.approx((2.0 + -3.0) / 2)
    assert matched["overall"]["within_5ppm"] == 1.0


def test_mass_accuracy_msms_matched_needs_the_sml_columns(data):
    handle = _open(data)           # この mzTab の SML には opt_global_spectrum_matched が無い
    out = json.loads(R.dataset_mass_accuracy(subset="msms_matched", dataset_handle=handle))
    # エラーにしない（手順の後退を止めない）。n=0 と、選べなかった理由を返す
    assert out["status"] == "success" and out["subset"] == "msms_matched"
    assert out["subset_available"] is False and "opt_global_spectrum_matched" in out["note"]
    assert out["overall"] == {"n": 0} and all(b["n"] == 0 for b in out["bands"])
    # 既定の rank1 は今までどおり動く（旗の有無の注記は付かない）
    plain = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))
    assert plain["subset"] == "rank1" and "subset_available" not in plain
    assert plain["theoretical_mz_source"] == "mztab"
    bad = json.loads(R.dataset_mass_accuracy(subset="nope", dataset_handle=handle))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


def test_mass_accuracy_msms_matched_is_rank1_only(match_data):
    handle = _open(match_data)
    out = json.loads(R.dataset_mass_accuracy(rank=2, subset="msms_matched", dataset_handle=handle))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert out["error"]["details"]["field"] == "rank"
    assert json.loads(R.dataset_mass_accuracy(rank=2, dataset_handle=handle))["subset"] == "rank1"


def test_injection_trend_unavailable_without_per_sample_mz(match_data):
    handle = _open(match_data)
    assert handles.registry.get(handle).dataset.per_sample_mz is None
    out = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))
    assert out["injection_trend"] == {"status": "unavailable"}


def _with_trend(match_data, ppm_of_sample):
    """特徴 1・3 の検出されたセルに、試料名の並び（注入順の代用）に比例する誤差を載せる。"""
    handle = _open(match_data)
    ds = handles.registry.get(handle).dataset
    theo = {"1": 200.0, "3": 400.0}
    ds.per_sample_mz = {fid: {name: mz * (1 + ppm_of_sample(i) * 1e-6)
                              for i, name in enumerate(sorted(ds.sample_names)) if i != 7}
                        for fid, mz in theo.items()}
    # 照合されていない特徴 2（msms_matched に入らない）には、まったく別の傾向（一定 +50 ppm）
    ds.per_sample_mz["2"] = {name: 300.0 * (1 + 50e-6) for name in ds.sample_names}
    return handle


def _with_sheet_order(handle):
    """実験情報シートに注入順（試料名の並びと同じ）を載せる。出所が filename でなくなる。"""
    ds = handles.registry.get(handle).dataset
    ds.sample_metadata_rows = [{"sample_id": name, "injection_order": i + 1}
                               for i, name in enumerate(sorted(ds.sample_names))]
    return handle


def test_injection_trend_spearman(match_data):
    handle = _with_sheet_order(_with_trend(match_data, lambda i: 1.0 + 0.5 * i))
    out = json.loads(R.dataset_mass_accuracy(subset="msms_matched", dataset_handle=handle))
    trend = out["injection_trend"]
    assert trend["spearman_rho"] == pytest.approx(1.0)
    assert "rho_unavailable_reason" not in trend
    assert trend["n_samples"] == 7                      # 8 試料のうち 1 試料はどの特徴も検出なし
    assert trend["injection_order_source"] == "metadata"
    per = trend["per_injection"]
    assert [p["injection_order"] for p in per] == sorted(p["injection_order"] for p in per)
    assert per[0]["sample"] == "HC1" and per[0]["n"] == 2
    assert per[0]["median_ppm"] == pytest.approx(1.0, abs=1e-3)
    # 特徴 2 の +50 ppm は subset の行に入らないので、傾向に混ざらない（試料 8 も載らないまま）
    assert max(p["median_ppm"] for p in per) < 10
    # rank1 では特徴 2 も数える（8 試料すべてに値が付き、1 試料目は 3 特徴）
    rank1 = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))["injection_trend"]
    assert rank1["n_samples"] == 8 and rank1["per_injection"][0]["n"] == 3
    flipped = _with_sheet_order(_with_trend(match_data, lambda i: 10.0 - i))
    assert json.loads(R.dataset_mass_accuracy(subset="msms_matched", dataset_handle=flipped)
                      )["injection_trend"]["spearman_rho"] == pytest.approx(-1.0)


def test_injection_trend_rho_is_unavailable_under_filename_order(match_data):
    """採取時刻（注入順）が読めずファイル名の順に落ちたときは、順位相関を出さない（設計 4.3(c)）。
    試料ごとの中央値と試料数は、傾向を見る材料として残す。"""
    handle = _with_trend(match_data, lambda i: 1.0 + 0.5 * i)   # 試料名の順に完全な単調増加
    trend = json.loads(R.dataset_mass_accuracy(subset="msms_matched", dataset_handle=handle)
                       )["injection_trend"]
    assert trend["injection_order_source"] == "filename"
    assert trend["spearman_rho"] is None
    assert trend["rho_unavailable_reason"] == "filename_order"
    assert trend["n_samples"] == 7 and len(trend["per_injection"]) == 7
    assert trend["per_injection"][0]["sample"] == "HC1"


def test_injection_trend_rho_is_none_below_five_samples(match_data):
    handle = _open(match_data)
    ds = handles.registry.get(handle).dataset
    ds.per_sample_mz = {"1": {name: 200.0 for name in sorted(ds.sample_names)[:4]}}
    trend = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))["injection_trend"]
    assert trend["n_samples"] == 4 and trend["spearman_rho"] is None
    # 出所が filename なので、その理由が先に付く（試料数が足りても足りなくても rho は出ない）
    assert trend["rho_unavailable_reason"] == "filename_order"
    _with_sheet_order(handle)
    trend = json.loads(R.dataset_mass_accuracy(dataset_handle=handle))["injection_trend"]
    assert trend["spearman_rho"] is None and trend["rho_unavailable_reason"] == "too_few_samples"


def test_novel_candidates_from_result_id_and_from_spec(data):
    handle = _open(data)
    diff = _compare_on(handle)
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    out = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert out["status"] == "success" and out["n_rows"] == 2
    assert out["group_ref"] == "HC"
    assert out["theoretical_mz_source"] == "mztab"      # mzTab-M の SME の値
    excluded = json.loads(R.dataset_novel_candidates(
        {**spec, "exclude": {"inchikey_blocks": [GABA_KEY]}}, dataset_handle=handle))
    assert excluded["n_rows"] == 1
    from_spec = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": _spec(data, features={"all": True}),
                        "analysis_id": "main"}, "q_max": 1.0}, dataset_handle=handle))
    assert from_spec["group_test"] == "RA" and from_spec["n_rows"] == 2
    assert out["rule_hash"] != from_spec["rule_hash"]


def test_novel_candidates_from_an_uncorrected_analysis_say_why_the_list_is_empty(data):
    """C7: 補正しない解析（correction=none）は q を持たないので、q_max で絞る候補は出ない。空の一覧が
    「該当なし」に読まれないよう notes に理由を書く（同じ spec の bh なら候補は出る）。"""
    handle = _open(data)
    uncorrected = _spec(data, features={"all": True})
    uncorrected["analyses"][0]["correction"] = "none"
    out = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": uncorrected, "analysis_id": "main"}, "q_max": 1.0},
        dataset_handle=handle))
    assert out["status"] == "success" and out["n_rows"] == 0
    assert any("correction=none" in n and "q_max" in n and "statistics.analysis_id=main" in n
               for n in out["notes"])
    corrected = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": _spec(data, features={"all": True}), "analysis_id": "main"},
         "q_max": 1.0}, dataset_handle=handle))
    assert corrected["n_rows"] == 2 and not any("correction=none" in n for n in corrected.get("notes", []))


def test_analysis_run_v1_notes_ineffective_settings_kept_in_the_spec(data, tmp_path):
    """C7: v1 の効かない設定（補正なしの q の閾値・drift_correct の無い注入順・ols の無い共変量）は断らず notes に書き、
    spec と rule_hash には残す。使わない注入順の表は読まない（欠けた試料・無いファイルで断らない）。"""
    handle = _open(data)
    order = tmp_path / "order.tsv"
    order.write_text("sample\torder\nHC1\t1\n", encoding="utf-8")      # ほかの試料は欠ける
    spec = _spec(data, injection_order={"source": "tsv", "path": str(order)},
                 covariates={"path": str(tmp_path / "missing.tsv"),
                             "columns": [{"name": "age", "type": "numeric"}]})
    spec["analyses"][0].update(correction="none", thresholds={"p": 0.05, "q": 0.05})
    out = json.loads(R.dataset_analysis_run(spec, dataset_handle=handle))
    assert out["status"] == "success", out
    heads = [n.split(" ")[0] for n in out["notes"]]
    assert heads == ["analyses[0].thresholds.q", "injection_order", "covariates"]
    assert out["analyses"][0]["n_meets_thresholds"] == 0
    assert "injection_order" not in out["input_fingerprint"]["files"]
    # 書いた設定は spec に残る（捨てない）ので、書かない spec とは前と同じく別の rule_hash
    bare = _spec(data)
    bare["analyses"][0].update(correction="none", thresholds={"p": 0.05, "q": 0.05})
    plain = json.loads(R.dataset_analysis_run(bare, dataset_handle=handle))
    assert plain["status"] == "success" and "notes" in plain
    assert out["rule_hash"] != plain["rule_hash"]
    assert out["analyses"][0]["rows"] == plain["analyses"][0]["rows"]


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


def _tsv(text):
    header, *lines = text.splitlines()
    names = header.split("\t")
    return names, [dict(zip(names, line.split("\t"))) for line in lines]


def test_feature_relations_columns_carry_strong_ambiguous_and_informational(data):
    handle = _open(data)
    names, ions = _tsv(json.loads(R.dataset_feature_relations("20", dataset_handle=handle))["ions_of_this"])
    assert names == ["ion", "relation", "dmz_mda", "drt", "r", "isotope_consistent", "strong",
                     "ambiguous", "informational"]
    by_relation = {(row["ion"], row["relation"]): row for row in ions}
    # MS-DIAL のリンクが無い mzTab でも、同時溶出・Δm/z・試料間の相関 r ≥ min_r がそろえば強い（B1-T17）
    na = by_relation[("21", "adduct:[M+Na]+/[M+H]1+")]
    assert (na["strong"], na["ambiguous"], na["informational"]) == ("true", "false", "false")
    iso = by_relation[("22", "isotope_M+1")]
    assert (iso["isotope_consistent"], iso["strong"], iso["ambiguous"]) == ("true", "true", "false")
    assert by_relation[("23", "insource:-H2O")]["strong"] == "true"
    names, as_ion = _tsv(json.loads(R.dataset_feature_relations("22", dataset_handle=handle))["as_ion_of"])
    assert names == ["of", "of_name", "relation", "dmz_mda", "drt", "r", "isotope_ratio",
                     "expected_ratio", "isotope_consistent", "strong", "ambiguous", "informational"]
    [row] = as_ion
    assert (row["of"], row["of_name"], row["relation"]) == ("20", "GABA", "isotope_M+1")
    # イオン C4H10NO2+ の M+1 の和: 13C 4·q_C + 2H 10·q_H + 15N 1·q_N + 17O 2·q_O（q = a / (1 − a)）
    # = 0.043263 + 0.001150 + 0.003653 + 0.000762 ≈ 0.0488。実測 5 / 100 = 0.05（1.5 倍以内）
    hand = (4 * 0.0107 / 0.9893 + 10 * 0.000115 / 0.999885 + 0.00364 / 0.99636
            + 2 * 0.00038 / 0.99757)
    assert hand == pytest.approx(0.0488, abs=1e-4)
    assert float(row["expected_ratio"]) == pytest.approx(hand, abs=1e-4)
    assert float(row["isotope_ratio"]) == pytest.approx(0.05, abs=1e-3)


def test_feature_relations_do_not_call_a_coeluting_isobar_an_adduct_of_itself(tmp_path):
    """mzTab-M の付加体は `[M+H]1+` と書かれる。既定の付加体 `[M+H]+` と綴りが違っても Y 自身の付加体なので、
    同時に溶出して強度が相関する同じ m/z の特徴を「Y の [M+H]+」（強い説明）と言わない。"""
    isobar = {"id": "27", "mz": 104.0716, "rt": 2.00, "values": _scaled(40), "sme": []}
    mztab = write_mztab(tmp_path / "iso" / "AlignResult-1.mzTab", features=(RELATION_FEATURES[0], isobar))
    out = json.loads(handle_tools.load_dataset(mztab_path=str(mztab), library="none"))
    assert out["status"] == "success", out
    handle = out["dataset_handle"]
    as_ion = json.loads(R.dataset_feature_relations("27", dataset_handle=handle))
    assert as_ion["n_coeluting"] == 1 and as_ion["as_ion_of"] == ""
    ions = json.loads(R.dataset_feature_relations("20", dataset_handle=handle))["ions_of_this"]
    assert "adduct:" not in ions


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
    handle = _open(data, msp=library)
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
    handle = _open(data, msp=library)
    out = json.loads(R.library_msms_qualify(
        {"target": {"inchikey": GABA_KEY}, "min_fragments_excl_precursor": 2},
        spot_id=20, dataset_handle=handle))
    assert out["verdict"] == "qualified" and len(out["rule_hash"]) == 64
    assert out["candidates"].splitlines()[1].split("\t")[1] == "GABA"
    missing = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=24,
                                                dataset_handle=handle))
    assert missing["verdict"] == "not_evaluable"


def test_diagnostic_ions_tool(data, library):
    handle = _open(data, msp=library)
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
    before = _analysis_state(handle)
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
    # 同じ入力なら同じファイル・出力先なしなら要約だけ・ハンドルの解析の状態（結果・行列）は変えない
    again = tmp_path / "again.json"
    R.dataset_alignment_qc(output_path=str(again), dataset_handle=handle)
    assert again.read_bytes() == target.read_bytes()
    inline = json.loads(R.dataset_alignment_qc(dataset_handle=handle))
    assert inline["status"] == "success" and "output_path" not in inline
    assert _analysis_state(handle) == before
    bad = json.loads(R.dataset_alignment_qc(ppm=-1, dataset_handle=handle))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


def test_alignment_qc_without_dataset_is_missing_state():
    out = json.loads(R.dataset_alignment_qc())
    assert out["error"]["code"] == "missing_state"


# ---------- 割れたスポット: 証拠の束・counts・新しい知見の候補（spec 2026-10-28 §4.3 (c)） ----------

CE_KEY = "RJECHNNFRHZQKU-RMUVNZEASA-N"
CE_CLAIM = {"claim_id": "C9",
            "target": {"inchikey": CE_KEY, "expected_mz": 650.6, "adduct": "[M+H]+", "rt_min": 20.0},
            "scan": {"ppm": [10], "rt_tol_min": [1.0]}}


def _split_handle(tmp_path, features=SPLIT_FEATURES):
    mztab = write_mztab(tmp_path / "split" / "AlignResult-1.mzTab", features=features)
    return _open({"mztab": mztab})


def test_claim_evidence_candidate_has_split(tmp_path, data):
    handle = _split_handle(tmp_path)
    out = json.loads(R.dataset_claim_evidence(CE_CLAIM, dataset_handle=handle))
    assert out["status"] == "success", out
    by_id = {c["feature_id"]: c for c in out["candidates"]}
    keys = ("group_id", "partners", "type", "mergeable")
    assert {k: by_id["30"]["split"][k] for k in keys} == {
        "group_id": "G1", "partners": ["31"], "type": "injection_block", "mergeable": True}
    assert by_id["31"]["split"]["partners"] == ["30"]
    # 割れの組に入らない候補は null
    plain = json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=_open(data)))
    assert all(c["split"] is None for c in plain["candidates"])


def test_counts_qualified_is_msms_verdict(data, library):
    handle = _open(data, msp=library)
    target = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0,
              "adduct_tiers": [["[M+H]+"]]}

    def counts(min_fragments):
        claim = {"claim_id": "C1", "target": target,
                 "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": min_fragments}}}
        out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
        assert out["default_selection"]["selected_feature_id"] == "20"
        verdict = {c["feature_id"]: c["msms"]["qualify"]["verdict"] for c in out["candidates"]}
        return out["default_selection"]["counts"], verdict["20"]

    strict, verdict = counts(3)                        # 断片は前駆体を除いて 2 本
    assert verdict == "not_qualified"
    assert strict["qualified"] == 0 and strict["msms_available"] == 1
    assert strict["rule_qualified"] == 1 and strict["msms_spectrum_available"] == 1
    loose, verdict = counts(2)
    assert verdict == "qualified"
    assert loose["qualified"] == 1 and loose["msms_available"] == 1
    # MS/MS を頼まなければ合否の数は無い
    bare = json.loads(R.dataset_claim_evidence({"claim_id": "C1", "target": target},
                                               dataset_handle=handle))
    assert bare["default_selection"]["counts"]["qualified"] is None
    assert bare["default_selection"]["counts"]["msms_available"] == 1
    assert bare["default_selection"]["counts"]["msms_spectrum_available"] is None


def test_spectrum_available_counts_mass_window_even_when_none_rule_qualified(data, library):
    # 目標の RT を特徴 20（RT 2.0）から外す: 規則で qualified は 0、質量の窓の 20 の代表 MS/MS は数える
    target = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 9.0, "adduct_tiers": [["[M+H]+"]]}
    claim = {"claim_id": "C1", "target": target, "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=_open(data, msp=library)))
    counts = out["default_selection"]["counts"]
    assert counts["rule_qualified"] == 0 and counts["msms_spectrum_available"] == 1


def test_spectrum_available_without_library(data, library):
    # ライブラリを開かない（library="none"）が .arf2 / .dcl はある（library フィクスチャが置く）
    target = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0, "adduct_tiers": [["[M+H]+"]]}
    claim = {"claim_id": "C1", "target": target, "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=_open(data)))
    c = {x["feature_id"]: x for x in out["candidates"]}["20"]
    assert c["msms"]["status"] == "no_library" and c["msms"]["spectrum_available"] is True


def test_spectrum_available_count_is_a_lower_bound_when_some_are_unknown():
    from metabolomix.analysis import claim_evidence

    def row(fid):
        return {"feature_id": fid, "qualified": False, "merged_from": None}

    match = {"status": "resolved", "reason": None, "selected_feature_id": "1",
             "selected_adduct": "[M+H]+", "selected_evidence": "sme_match",
             "adduct_tiers_used": 1, "adduct_tier_used": 1,
             "counts": {"mass": 3, "rt": 3, "qualified": 0, "by_evidence": {}},
             "candidates": [row("1"), row("2"), row("3")]}

    def count(spectra):
        return claim_evidence._selection(match, None, spectra)["counts"]["msms_spectrum_available"]

    assert count({"1": None, "2": None, "3": None}) is None
    assert count({"1": True, "2": None, "3": False}) == 1       # 下限
    assert count({"1": False, "2": False, "3": False}) == 0


def test_spectrum_available_is_unknown_when_the_files_are_missing():
    from metabolomix.analysis import claim_evidence

    missing = {"status": "not_found", "reason": "no_arf2"}
    assert claim_evidence._spectrum_available(missing) is None
    assert claim_evidence._spectrum_available({"status": "not_found", "reason": "no_dcl"}) is None
    assert claim_evidence._spectrum_available({"status": "not_found", "reason": "no_spectrum"}) is False
    assert claim_evidence._spectrum_available({"status": "not_found"}) is False


def test_counts_qualified_counts_a_merged_pair_by_either_spot():
    """割れの組は、どちらのスポットが MS/MS で qualified でも 1。相手が窓の外（規則で qualified で
    ない merge_partner）でもよい。"""
    from metabolomix.analysis import claim_evidence

    def row(fid, qualified):
        return {"feature_id": fid, "qualified": qualified, "merged_from": ["30", "31"]}

    match = {"status": "resolved", "reason": None, "selected_feature_id": "30",
             "selected_adduct": "[M+H]+", "selected_evidence": "sme_match",
             "adduct_tier_used": 1, "merged_from": ["30", "31"],
             "counts": {"mass": 1, "rt": 1, "qualified": 1, "by_evidence": {}},
             # 窓に入って規則で qualified なのは 31 だけ。30 は merge_partner（qualified でない）
             "candidates": [row("31", True), row("30", False)]}
    got = claim_evidence._selection(match, {"30": "qualified", "31": "not_qualified"})
    assert got["counts"]["qualified"] == 1
    assert claim_evidence._selection(match, {"30": "not_qualified",
                                             "31": "not_qualified"})["counts"]["qualified"] == 0
    assert claim_evidence._selection(match, None)["counts"]["qualified"] is None
    # 組のどちらかに代表 MS/MS があれば 1（相手が merge_partner でも）。頼んでいなければ None
    assert claim_evidence._selection(match, None, {"30": True, "31": False}
                                     )["counts"]["msms_spectrum_available"] == 1
    assert claim_evidence._selection(match, None, {"30": False, "31": False}
                                     )["counts"]["msms_spectrum_available"] == 0
    # True が無く、分からない（None）が 1 つでもあれば 0 ではなく None（#12 の再発を防ぐ）
    assert claim_evidence._selection(match, None, {"30": False, "31": None}
                                     )["counts"]["msms_spectrum_available"] is None
    assert claim_evidence._selection(match, None, {"30": None, "31": None}
                                     )["counts"]["msms_spectrum_available"] is None
    assert claim_evidence._selection(match, None)["counts"]["msms_spectrum_available"] is None
    # 組でない候補は従来どおり自分のスポットだけで数える
    single = {**match, "merged_from": None,
              "candidates": [{"feature_id": "31", "qualified": True, "merged_from": None}]}
    assert claim_evidence._selection(single, {"30": "qualified", "31": "not_qualified"}
                                     )["counts"]["qualified"] == 0


def _tsv_rows(text):
    lines = text.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _groups_tsv(tmp_path):
    groups = tmp_path / "groups_split.tsv"
    groups.write_text("sample\tgroup\trole\n" + "".join(
        f"{s}\t{s[:2] if not s.startswith('QC') else ''}\t{'qc' if s.startswith('QC') else 'sample'}\n"
        for s in SAMPLES), encoding="utf-8")
    return groups


#: SPLIT_FEATURES と同じ割れ方（最初の 2 注入だけが 31 に割れる）で、まとめた値が HC で低く RA で高い
#: （統計エンジンの判定が changed になり、新しい知見の候補に載る）。
SPLIT_CHANGED = (
    {**SPLIT_FEATURES[0], "values": [20, 21, 100, 1010, 1020, 1030, 500, 505]},
    {**SPLIT_FEATURES[1], "values": [105, 98, 22, 20, 21, 23, 20, 22]},
    SPLIT_FEATURES[2],
)


def test_novel_candidates_count_split_once(tmp_path):
    handle = _split_handle(tmp_path, SPLIT_CHANGED)
    diff = _compare_on(handle)
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    out = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    assert out["status"] == "success", out
    [row] = _tsv_rows(out["rows"])                     # 30 と 31 は 1 つの物質として 1 行
    assert row["feature_id"] == "30" and row["merged_from"] == "30;31"
    # q はまとめた行列で検定し直した値（dataset_feature_search の stats と同じ）
    search = json.loads(T.dataset_feature_search(feature_ids=["30"], dataset_handle=handle))
    [stat] = _tsv_rows(search["stats"])
    [unmerged] = _tsv_rows(search["stats_unmerged"])
    assert abs(float(row["q"]) - float(stat["q"])) <= 1e-3 * float(stat["q"])
    assert stat["q"] != unmerged["q"]
    off = json.loads(R.dataset_novel_candidates(
        {**spec, "statistics": {"result_id": diff["result_id"], "merge_split": False}},
        dataset_handle=handle))
    assert off["rule_hash"] != out["rule_hash"]
    # merge_split=false は v7 と同じ形・同じ規則（merged_from 列も merge_split も無い）
    assert "merge_split" not in off
    assert off["rows"].splitlines()[0].split("\t")[-1] == "direction"
    from metabolomix.analysis import mass_accuracy
    v7 = mass_accuracy.normalize_novel_spec({**spec, "statistics": {"result_id": diff["result_id"],
                                                                    "merge_split": False}})
    assert v7["statistics"] == {"result_id": diff["result_id"]}
    assert off["rule_hash"] == mass_accuracy.rule_hash(v7)

    # exclude.feature_ids に消える側の ID（31）だけを書いても、組ごと外れる
    gone = json.loads(R.dataset_novel_candidates({**spec, "exclude": {"feature_ids": ["31"]}},
                                                 dataset_handle=handle))
    assert gone["n_rows"] == 0

    # SME（rank 1）が消える側（31）にしか無くても、残す ID の 1 候補として数える
    moved = (dict(SPLIT_FEATURES[0], sme=[]), dict(SPLIT_FEATURES[1], sme=SPLIT_FEATURES[0]["sme"]),
             SPLIT_FEATURES[2])
    other = tmp_path / "moved"
    other.mkdir()
    handle = _split_handle(other, moved)
    analysis = {"groups": {"source": "tsv", "path": str(_groups_tsv(tmp_path)), "role_column": "role"},
                "comparison": {"reference": "HC", "test": "RA"}, "features": {"all": True},
                "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]}
    written = tmp_path / "novel.json"
    from_spec = json.loads(R.dataset_novel_candidates(
        {"statistics": {"analysis_spec": analysis, "analysis_id": "main"}, "q_max": 1.0},
        dataset_handle=handle, output_path=str(written)))
    assert from_spec["status"] == "success", from_spec
    [saved] = json.loads(written.read_text(encoding="utf-8"))["rows"]
    assert saved["feature_id"] == "30" and saved["sme_feature_id"] == "31"
    assert saved["merged_from"] == ["30", "31"] and saved["name"] == "CE 18:2"


def test_novel_candidates_rule_hash_includes_merge_groups(tmp_path):
    """まとめた割れの組の hash を規則に入れる（まとめるときだけ。false は v7 と同じ hash）。"""
    from metabolomix.analysis import mass_accuracy

    handle = _split_handle(tmp_path)
    diff = _compare_on(handle)
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    out = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    groups_hash = out["merge_split"]["groups_hash"]
    assert len(groups_hash) == 64
    normalized = mass_accuracy.normalize_novel_spec(spec)
    assert out["rule_hash"] == mass_accuracy.rule_hash(normalized, groups_hash)
    assert out["rule_hash"] != mass_accuracy.rule_hash(normalized)
    # 組が変われば別の規則（同じ spec でも）
    assert mass_accuracy.rule_hash(normalized, "0" * 64) != out["rule_hash"]
    # merge_split=false は groups_hash を載せず、v7 と同じ hash
    off = json.loads(R.dataset_novel_candidates(
        {**spec, "statistics": {"result_id": diff["result_id"], "merge_split": False}},
        dataset_handle=handle))
    assert "merge_split" not in off
    v7 = mass_accuracy.normalize_novel_spec({**spec, "statistics": {
        "result_id": diff["result_id"], "merge_split": False}})
    assert off["rule_hash"] == mass_accuracy.rule_hash(v7)


def test_claim_evidence_default_selection_counts_a_merged_pair_once(tmp_path):
    """組の 2 スポットが窓に入れば既定の特徴は残す ID。無関係な 2 スポットは needs_input のまま。"""
    features = ({**SPLIT_FEATURES[0], "sme": []}, SPLIT_FEATURES[1], SPLIT_FEATURES[2],
                {"id": "40", "mz": 500.0000, "rt": 8.00, "values": [100, 110, 90, 120, 105, 95, 100, 101],
                 "sme": []},
                {"id": "41", "mz": 500.0010, "rt": 8.20, "values": [98, 112, 91, 118, 104, 97, 99, 100],
                 "sme": []})
    handle = _split_handle(tmp_path, features)
    claim = {"claim_id": "C9", "target": {"inchikey": CE_KEY, "expected_mz": 650.6,
                                          "adduct": "[M+H]+", "rt_min": 20.3, "rt_tol_min": 0.5}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    selection = out["default_selection"]
    assert selection["status"] == "resolved" and selection["selected_feature_id"] == "30"
    assert selection["merged_from"] == ["30", "31"]
    assert selection["counts"]["msms_available"] == 1
    assert {c["feature_id"] for c in out["candidates"]} == {"30", "31"}   # 候補はスポットごと
    unrelated = {"claim_id": "C10", "target": {"inchikey": CE_KEY, "expected_mz": 500.0,
                                               "adduct": "[M+H]+", "rt_min": 8.1}}
    other = json.loads(R.dataset_claim_evidence(unrelated, dataset_handle=handle))
    assert other["default_selection"]["status"] == "needs_input"
    assert other["default_selection"]["merged_from"] is None


def test_claim_evidence_adds_the_partner_outside_the_mass_window(tmp_path):
    """組の片方（31）だけが質量の窓に入っても、選んだ残す ID（30）が候補に入る（レビューの再現）。"""
    features = ({**SPLIT_FEATURES[0], "sme": []}, SPLIT_FEATURES[1], SPLIT_FEATURES[2])
    handle = _split_handle(tmp_path, features)
    claim = {"claim_id": "C", "target": {"inchikey": CE_KEY, "expected_mz": 650.608,
                                         "adduct": "[M+H]+", "rt_min": 20.6, "rt_tol_min": 0.3}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["default_selection"]["selected_feature_id"] == "30"
    by_id = {c["feature_id"]: c for c in out["candidates"]}
    assert set(by_id) == {"30", "31"}
    assert by_id["30"]["scans"] == ["merge_partner"]
    assert "default" in by_id["31"]["scans"]
    assert by_id["30"]["split"]["partners"] == ["31"]
    assert out["merge_split"] is True and len(out["merge_groups_hash"]) == 64

    # merge_split=false は v7 と同じ規則（同じ rule_hash）で、組を 1 つにしない・相手を足さない
    from metabolomix.analysis import claim_evidence
    off = json.loads(R.dataset_claim_evidence({**claim, "merge_split": False},
                                              dataset_handle=handle))
    assert off["merge_split"] is False and off["merge_groups_hash"] is None
    assert "merge_split" not in off["rule"]
    assert off["default_selection"]["selected_feature_id"] == "31"
    assert {c["feature_id"] for c in off["candidates"]} == {"31"}
    normalized = claim_evidence.normalize_claim_spec(claim, dataset_polarity="positive")
    v7 = {k: v for k, v in normalized.items() if k != "merge_split"}
    assert off["rule_hash"] == claim_evidence.rule_hash(v7)
    assert out["rule_hash"] != off["rule_hash"]


def test_novel_candidates_labels_column_and_counts(data, tmp_path):
    handle = _open(data)
    diff = _compare_on(handle)
    spec = {"statistics": {"result_id": diff["result_id"]}, "q_max": 1.0}
    plain = json.loads(R.dataset_novel_candidates(spec, dataset_handle=handle))
    labels = tmp_path / "labels.tsv"
    first = _tsv_rows(plain["rows"])[0]["feature_id"]
    labels.write_text(f"feature_id\tlevel\n{first}\tA\n999999\tB\n", encoding="utf-8")
    out = json.loads(R.dataset_novel_candidates(
        {**spec, "labels": {"path": str(labels), "name": "xenobiotic"}}, dataset_handle=handle))
    assert out["n_rows"] == plain["n_rows"] and out["rule_hash"] == plain["rule_hash"]
    assert out["rows"].splitlines()[0].split("\t")[-1] == "labels.xenobiotic"
    assert ([r["feature_id"] for r in _tsv_rows(out["rows"])]
            == [r["feature_id"] for r in _tsv_rows(plain["rows"])])
    assert out["labels"]["counts"]["all"]["A"] == 1
    assert out["labels"]["n_unmatched"] == 1
    assert out["labels"]["name"] == "xenobiotic" and out["labels"]["n_rows"] == 2
    assert out["labels"]["order"] == ["A", "B", "C", "D"] and len(out["labels"]["sha256"]) == 64


# ---------- ppm_offset（質量の窓の中心のずれ） ----------

#: ppm_offset を足す前（v8）の実装で CLAIM を処理した rule_hash。陽イオンの既定の付加体の
#: 第 2 段に [M+CH3OH+H]+ と [M+H-2H2O]+ を足したとき（A12）に段の並びが変わったので
#: 書き直した（旧 41f68e94… は旧い段で再現することを確かめた）。ppm_offset を省いた
#: ときに規則が変わらないことを縛るのは、この値と ppm_offset の有無の比較
_CLAIM_RULE_HASH_V8 = "10003a763b0bc7c19f996242b458b66567e2a550fcbfde582ecd43e3b92d7292"


def test_claim_evidence_ppm_offset_matches_feature_search(data):
    handle = _open(data)
    offset = {"bands": [{"max_mz": 200, "offset": -8.0}]}
    claim = {**CLAIM, "ppm_offset": offset}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    search = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", **CLAIM["target"]}], ppm_offset=offset, dataset_handle=handle))
    assert out["default_selection"]["status"] == search["targets"][0]["status"]
    assert out["default_selection"]["selected_feature_id"] == search["targets"][0]["selected_feature_id"]
    assert out["rule"]["ppm_offset"]["bands"][0]["offset"] == -8.0


def test_claim_evidence_without_offset_keeps_rule_hash(data):
    handle = _open(data)
    before = _CLAIM_RULE_HASH_V8
    assert json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=handle))["rule_hash"] == before
    shifted = json.loads(R.dataset_claim_evidence({**CLAIM, "ppm_offset": {"default": -3}},
                                                  dataset_handle=handle))
    assert shifted["rule_hash"] != before


def test_claim_evidence_bad_offset_is_spec_invalid(data):
    out = json.loads(R.dataset_claim_evidence({**CLAIM, "ppm_offset": {"default": "x"}},
                                              dataset_handle=_open(data)))
    assert out["error"]["code"] == "CLAIM_SPEC_INVALID"


# ---------- rt_offset_min（RT の窓の中心のずれ） ----------

#: rt_min を特徴 20（RT 2.00）から 0.5 分ずらした目標。rt_tol_min 0.3 では窓 [1.2, 1.8] に入らない。
SHIFTED_TARGET = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 1.5}


def test_claim_evidence_rt_offset_matches_feature_search(data):
    handle = _open(data)
    claim = {"claim_id": "C1", "target": SHIFTED_TARGET, "scan": {"ppm": [10, 20], "rt_tol_min": [0.3]},
             "rt_offset_min": 0.5}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    search = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", **SHIFTED_TARGET}], rt_offset_min=0.5, dataset_handle=handle))
    assert out["default_selection"]["status"] == search["targets"][0]["status"] == "resolved"
    assert out["default_selection"]["selected_feature_id"] == search["targets"][0]["selected_feature_id"] == "20"
    # 走査の全変種にも同じずれが効く
    assert [s["selected_feature_id"] for s in out["scans"]] == ["20", "20"]
    assert out["rule"]["rt_offset_min"] == 0.5 and out["rule"]["target"]["rt_offset_min"] == 0.5
    assert "notes" not in out
    # ずらさなければ RT の窓の外（判定が変わる）
    plain = json.loads(R.dataset_claim_evidence({k: v for k, v in claim.items() if k != "rt_offset_min"},
                                                dataset_handle=handle))
    assert plain["default_selection"]["selected_feature_id"] != "20"
    assert "rt_offset_min" not in plain["rule"] and "rt_offset_min" not in plain["rule"]["target"]


def test_claim_evidence_rt_offset_enters_the_rule_hash_only_when_set(data):
    handle = _open(data)
    # 省けば従来と同じ hash（v8 の値）。渡せば（0 でも）規則の一部
    assert json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=handle))["rule_hash"] == _CLAIM_RULE_HASH_V8
    for value in (0, -0.2):
        shifted = json.loads(R.dataset_claim_evidence({**CLAIM, "rt_offset_min": value},
                                                      dataset_handle=handle))
        assert shifted["rule_hash"] != _CLAIM_RULE_HASH_V8
    assert (json.loads(R.dataset_claim_evidence({**CLAIM, "rt_offset_min": 0.1}, dataset_handle=handle))["rule_hash"]
            != json.loads(R.dataset_claim_evidence({**CLAIM, "rt_offset_min": 0.2}, dataset_handle=handle))["rule_hash"])


def test_claim_evidence_rt_offset_without_rt_min_is_ignored_with_a_note(data, tmp_path):
    handle = _open(data)
    bare = {"claim_id": "C2", "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2"}}
    plain = json.loads(R.dataset_claim_evidence(bare, dataset_handle=handle))
    out = json.loads(R.dataset_claim_evidence({**bare, "rt_offset_min": 0.5}, dataset_handle=handle))
    assert "error" not in out
    assert out["rule_hash"] == plain["rule_hash"] and "rt_offset_min" not in out["rule"]
    assert any("rt_offset_min" in n and "rt_min" in n for n in out["notes"])
    assert "notes" not in plain
    # 理由は書いた JSON と要約にも残る
    path = tmp_path / "c2.json"
    summary = json.loads(R.dataset_claim_evidence({**bare, "rt_offset_min": 0.5}, dataset_handle=handle,
                                                  output_path=str(path)))
    assert summary["notes"] == out["notes"]
    assert json.loads(path.read_text(encoding="utf-8"))["notes"] == out["notes"]


@pytest.mark.parametrize("value", ["x", True, float("nan"), float("inf")])
def test_claim_evidence_bad_rt_offset_is_spec_invalid(data, value):
    out = json.loads(R.dataset_claim_evidence({**CLAIM, "rt_offset_min": value},
                                              dataset_handle=_open(data)))
    assert out["error"]["code"] == "CLAIM_SPEC_INVALID"
    assert out["error"]["details"]["field"] == "rt_offset_min"


# ---------- 候補が多いときの戻り値（断らず先頭だけ） ----------

def test_claim_evidence_inline_truncates_with_the_selected_feature_first(data, tmp_path, monkeypatch):
    handle = _open(data)
    # [M+Na]+ だけの目標で、RT 5.0 の 25 を既定の特徴にする（特徴の並びでは 21 が先）
    claim = {"claim_id": "C5", "target": {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 5.0,
                                          "adduct_tiers": [["[M+Na]+"]]}}
    full = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert full["default_selection"]["selected_feature_id"] == "25"
    assert [c["feature_id"] for c in full["candidates"]] == ["21", "25"]
    assert full["truncated"] == 0 and "n_candidates_shown" not in full     # 上限内は束のまま

    monkeypatch.setattr(R, "MAX_INLINE_CANDIDATES", 1)
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["status"] == "success" and "error" not in out
    assert [c["feature_id"] for c in out["candidates"]] == ["25"]
    assert out["truncated"] is True
    assert out["n_candidates"] == 2 and out["n_candidates_shown"] == 1
    assert out["candidates_order"] == "default_selection_first"
    assert out["max_candidates_truncated"] == 0
    assert any("output_path" in n for n in out["notes"])
    # 判定・規則・走査は束と同じ（候補の一覧だけを切る）
    for key in ("default_selection", "scans", "rule_hash", "rule", "input_fingerprint"):
        assert out[key] == full[key]
    # output_path は全件を書く（束の truncated は max_candidates の打ち切りの数のまま）
    path = tmp_path / "c5.json"
    summary = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle, output_path=str(path)))
    assert summary["n_candidates"] == 2 and summary["truncated"] == 0
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert [c["feature_id"] for c in saved["candidates"]] == ["21", "25"] and saved["truncated"] == 0
    assert "notes" not in saved


def test_claim_evidence_inline_cut_reports_the_max_candidates_cut_too(data, monkeypatch):
    handle = _open(data)
    monkeypatch.setattr(R, "MAX_INLINE_CANDIDATES", 1)
    claim = {**CLAIM, "max_candidates": 2}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["truncated"] is True and out["n_candidates"] == 2
    assert out["max_candidates_truncated"] > 0
    assert any("max_candidates=2" in n for n in out["notes"])
    assert out["candidates"][0]["feature_id"] == out["default_selection"]["selected_feature_id"]


def test_claim_evidence_on_an_lbm_only_dataset_is_not_evaluable_not_no_library(data, library):
    """同定が LBM だけのデータセット（MSP が無い）。no_library は library_load を促すので使わない。"""
    from metabolomix.core import handles

    handle = _open(data)
    handles.registry.get(handle).dataset.library_record = {
        "status": "lbm_only", "kind": "lbm", "file": "lab.lbm2", "path": None}
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    msms = {c["feature_id"]: c["msms"] for c in out["candidates"]}
    assert {m["status"] for m in msms.values()} == {"not_evaluable"}
    assert all(m["reason"] == "lbm_only" for m in msms.values())
    assert "LBM" in out["msms_not_evaluable_reason"]
    assert msms["20"]["spectrum_available"] is True


def test_not_evaluable_reason_is_last_key(data, library):
    """理由の文は束の先頭に置かず、末尾に置く（読み手が最初に見るのは status / schema）。"""
    from metabolomix.core import handles

    handle = _open(data)
    handles.registry.get(handle).dataset.library_record = {
        "status": "lbm_only", "kind": "lbm", "file": "lab.lbm2", "path": None}
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    keys = list(out)
    assert "msms_not_evaluable_reason" in keys
    assert keys[-1] == "msms_not_evaluable_reason" and keys[0] != "msms_not_evaluable_reason"
    assert "schema" in keys


# ---------- dataset_analysis_run: v1 / v2 の振り分け（C の配線 Task 7） ----------

def _v2_state(tmp_path):
    from tests.test_stats_results import _ds
    state = _ds()
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


_V2_OFF = {"drift": {"mode": "off"},
           "filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
           "normalize": {"method": "none"}}


def _v2_spec(**extra):
    spec = {"schema": "analysis-spec.v2", "preprocess": _V2_OFF,
            "comparisons": [{"id": "T_vs_C", "factor": "group", "reference": "C", "test": "T"}],
            "analyses": [{"analysis_id": "cmp", "kind": "compare"}]}
    spec.update(extra)
    return spec


def _strict_loads(text):
    return json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))


def test_analysis_run_dispatches_v2_by_schema(tmp_path):
    _v2_state(tmp_path)
    out = json.loads(R.dataset_analysis_run(_v2_spec(analyses=[
        {"analysis_id": "cmp", "kind": "compare"}, {"analysis_id": "pc", "kind": "pca"},
        {"analysis_id": "comp", "kind": "composition"}], primary="cmp")))
    assert out["status"] == "success" and out["schema"] == "analysis-result.v2"
    assert {"rule_hash", "matrix_id", "input_fingerprint", "warnings"} <= out.keys()
    cmp_, pca_, comp = out["analyses"]
    assert (cmp_["analysis_id"], cmp_["kind"], cmp_["n_rows"], cmp_["n_tested"]) == ("cmp", "compare", 10, 10)
    assert sum(cmp_["verdicts"].values()) == 10 and cmp_["verdicts"]["changed"] >= 3
    assert "n_q_below" in cmp_ and cmp_["rows"].splitlines()[0].split("\t")[:2] == ["contrast_id", "feature_id"]
    assert pca_["kind"] == "pca" and pca_["explained_variance_ratio"] and "rows" not in pca_
    assert comp["kind"] == "composition" and "n_rows" in comp
    # v1 の spec（schema なし）は v1 のまま
    v1 = json.loads(R.dataset_analysis_run({"comparison": {}}))
    assert v1["error"]["code"] in ("ANALYSIS_SPEC_INVALID", "missing_state")


def test_analysis_run_v1_unchanged_golden(data, tmp_path):
    from metabolomix.analysis import analysis_spec
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    assert list(out) == ["status", "schema", "rule_hash", "merge_split", "input_fingerprint",
                         "comparison", "group_ref", "group_test", "analyses", "note"]
    assert out["schema"] == "analysis-result.v1"
    assert list(out["analyses"][0]) == [
        "analysis_id", "n_features_in", "n_features_out", "bh_population_n", "n_rows", "n_tested",
        "n_meets_thresholds", "detection_available", "n_detection_driven", "n_constant_in_group",
        "rows", "truncated"]
    assert out["analyses"][0]["truncated"] is False
    path = tmp_path / "v1.json"
    R.dataset_analysis_run(_spec(data), dataset_handle=handle, output_path=str(path))
    ds = handles.registry.get(handle).dataset
    expected = R.write_result(str(tmp_path / "exp.json"), analysis_spec.run_spec(
        analysis_spec.normalize_spec(_spec(data)), ds))
    assert path.read_bytes() == (tmp_path / "exp.json").read_bytes() and expected


def test_analysis_run_v2_output_path_strict_json(tmp_path):
    _v2_state(tmp_path)
    path = tmp_path / "out" / "r.json"
    out = json.loads(R.dataset_analysis_run(_v2_spec(), output_path=str(path)))
    assert out["output_path"] == str(path.resolve()) and "rows" not in out["analyses"][0]
    saved = _strict_loads(path.read_text(encoding="utf-8"))
    assert saved["schema"] == "analysis-result.v2" and saved["rule_hash"] == out["rule_hash"]
    assert len(saved["analyses"][0]["result"]["rows"]) == 10
    again = tmp_path / "out" / "r2.json"
    R.dataset_analysis_run(_v2_spec(), output_path=str(again))
    assert again.read_bytes() == path.read_bytes()


def test_analysis_run_v2_composition_summary_carries_the_inner_confounding(tmp_path):
    """C4: 既定の allow_confounded=true は完全交絡でも組成を比べる。止める代わりの印（交絡）は要約に残る:
    組成の要約に内側の比較の `confounding` と、その注意（完全交絡ほか）をコードで重ねずに載せる。"""
    state = _v2_state(tmp_path)
    for k, fid in enumerate(state.feature_ids):
        state.feature_metadata[fid]["ontology"] = "PC" if k < 3 else ("PE" if k < 6 else "TG")
    for row in state.sample_metadata_rows:
        if row.get("group") == "T":
            row["batch"] = "B2"          # C は全部 B1・T は全部 B2: どのバッチも 1 群だけ＝完全交絡
    spec = _v2_spec(analyses=[{"analysis_id": "cmp", "kind": "compare"},
                              {"analysis_id": "comp", "kind": "composition"}], primary="comp")
    path = tmp_path / "full.json"
    full = json.loads(R.dataset_analysis_run(spec, output_path=str(path)))
    out = json.loads(R.dataset_analysis_run(spec))
    assert out["status"] == "success" and full["analyses"] == [
        {k: v for k, v in a.items() if k not in ("rows", "truncated", "rows_order")} for a in out["analyses"]]
    cmp_, comp = out["analyses"]
    saved = _strict_loads(path.read_text(encoding="utf-8"))
    inner = saved["analyses"][1]["result"]
    assert inner["absolute"]["confounding"]["batch_status"] == "complete"
    expected = {"batch_status": "complete", "by_contrast": {"T_vs_C": {"batch_status": "complete",
                                                                       "confounded": True}}}
    assert comp["confounding"] == expected and cmp_["confounding"] == expected
    codes = [w["code"] for w in comp["warnings"]]
    assert codes.count("batch_group_complete_confounding") == 1          # 2 つの表で重ねない
    warn = next(w for w in comp["warnings"] if w["code"] == "batch_group_complete_confounding")
    assert warn["tables"] == ["absolute", "clr"]
    inner_codes = {w["code"] for name in ("absolute", "clr") for w in inner[name]["warnings"]}
    assert inner_codes <= set(codes) and len(codes) == len(set(codes))
    assert "batch_group_complete_confounding" in {w["code"] for w in cmp_["warnings"]}


def _lipid_mixed_state(tmp_path, assay_kind):
    """脂質の名前の特徴と、研究室 MSP の天然物の名前の特徴（f6・f7）と、研究室 MSP の名前が ID の TG（f9。
    ClassyFire のクラスと組成式から TG に救われる）が混ざったデータセット。"""
    state = _v2_state(tmp_path)
    names = {"f0": "PC 34:1", "f1": "PC 36:2", "f2": "PC 38:4", "f3": "PE 36:2", "f4": "PE 38:4",
             "f5": "PE 40:6", "f6": "Ginsenoside Rb1", "f7": "no MS2: Saikosaponin A", "f8": "TG 52:2",
             "f9": "MassSpecGymID0000001"}
    for fid, name in names.items():
        onto = ("Triterpene saponins" if fid in ("f6", "f7") else
                "Triacylglycerols" if fid == "f9" else name.split()[0])
        state.feature_metadata[fid].update(name=name, ontology=onto)
    state.feature_metadata["f9"]["formula"] = "C57H104O6"
    state.assay_kind = assay_kind
    return state


def test_analysis_run_v2_lipid_composition_keeps_non_lipid_names_out_of_the_classes(tmp_path):
    spec = _v2_spec(analyses=[{"analysis_id": "comp", "kind": "composition"}])
    _lipid_mixed_state(tmp_path, "lipid")
    path = tmp_path / "lipid.json"
    out = json.loads(R.dataset_analysis_run(spec, output_path=str(path)))
    assert out["status"] == "success" and out["assay_kind"] == "lipid"
    comp = out["analyses"][0]
    assert comp["class_label_kind"] == "lipid_class" and comp["n_classes"] == 3
    assert comp["non_lipid"] == {"label": "non-lipid", "n_features": 2,
                                 "by_ontology": {"Triterpene saponins": 2}}
    assert comp["lipid_class_relabel"] == {
        "n_features": 1, "n_rescued": 1,
        "moves": [{"from": "Triacylglycerols", "to": "TG", "route": "class_formula", "n_features": 1,
                   "n_rescued": 1}]}
    inner = _strict_loads(path.read_text(encoding="utf-8"))["analyses"][0]["result"]
    assert set(inner["classes"]) == {"PC", "PE", "TG"} and inner["classes"]["TG"] == ["f8", "f9"]
    assert inner["non_lipid"]["feature_ids"] == ["f6", "f7"]
    assert {r["feature_id"] for r in inner["absolute"]["rows"]} == {"PC", "PE", "TG"}
    # 一般代謝物はこれまでどおり: クラスはどの ontology も数え、non_lipid・クラスの移し替えは無い
    _lipid_mixed_state(tmp_path, "metabolite")
    out_m = json.loads(R.dataset_analysis_run(spec))
    comp_m = out_m["analyses"][0]
    assert out_m["assay_kind"] == "metabolite" and comp_m["n_classes"] == 5
    assert "non_lipid" not in comp_m and "lipid_class_relabel" not in comp_m


def test_analysis_run_v2_inline_row_limit(tmp_path):
    import numpy as np
    state = _v2_state(tmp_path)
    nf = 60
    rng = np.random.default_rng(1)
    state.feature_matrix = np.exp2(rng.normal(20, 1, (nf, len(state.sample_names))))
    state.detected_mask = np.ones_like(state.feature_matrix, dtype=bool)
    state.feature_ids = [f"f{k}" for k in range(nf)]
    state.feature_metadata = {f"f{k}": {"name": f"F{k}"} for k in range(nf)}
    out = json.loads(R.dataset_analysis_run(_v2_spec()))
    item = out["analyses"][0]
    assert "error" not in out and item["n_rows"] == 60
    assert item["truncated"] is True and item["rows_order"] == "p_ascending"
    assert any("['cmp']" in n for n in out["notes"])
    lines = item["rows"].splitlines()
    header = lines[0].split("\t")
    ps = [float(dict(zip(header, line.split("\t")))["p"]) for line in lines[1:]]
    assert len(ps) == R.MAX_INLINE_ROWS and ps == sorted(ps)
    ok = json.loads(R.dataset_analysis_run(_v2_spec(), output_path=str(tmp_path / "o.json")))
    assert ok["analyses"][0]["n_rows"] == 60 and "truncated" not in ok["analyses"][0]


def test_analysis_run_v2_notes_inert_settings_instead_of_refusing(tmp_path):
    """dataset_compare と同じく、効かない設定は断らない: reported 以外の母集団の reported_feature_ids・
    normalize=none の係数の特徴・drift=off の loess_degree は捨てて notes に書き、結果と rule_hash は書かない spec と同じ。"""
    _v2_state(tmp_path)
    plain = json.loads(R.dataset_analysis_run(_v2_spec()))
    inert = _v2_spec(
        preprocess={**_V2_OFF, "drift": {"mode": "off", "loess_degree": 2},
                    "normalize": {"method": "none", "features": ["f0"]}},
        analyses=[{"analysis_id": "cmp", "kind": "compare", "reported_feature_ids": ["f0", "f1"]}])
    out = json.loads(R.dataset_analysis_run(inert))
    assert "error" not in out and out["status"] == "success" and "notes" not in plain
    assert (out["rule_hash"], out["matrix_id"]) == (plain["rule_hash"], plain["matrix_id"])
    assert out["analyses"][0]["rows"] == plain["analyses"][0]["rows"]
    assert len(out["notes"]) == 3
    assert any("analyses[0].reported_feature_ids" in n for n in out["notes"])
    assert any("loess_degree" in n for n in out["notes"])
    assert any("normalization_features" in n for n in out["notes"])


def test_analysis_run_v2_normalized_with_qc_does_not_ask(tmp_path):
    """drift=auto は換算済み・QC ありでも尋ねない（needs_input を返さずに実行する。2026-11-14 の C3）。"""
    state = _v2_state(tmp_path)
    state.quant_state = {"normalized": True, "method": "is_ratio", "unit": "ratio",
                         "record_id": "q_test_1"}
    spec = _v2_spec(preprocess={"filters": _V2_OFF["filters"], "normalize": {"method": "none"}})
    out = json.loads(R.dataset_analysis_run(spec, output_path=str(tmp_path / "o.json")))
    assert "needs_input" not in out and "error" not in out, out
    assert out["analyses"]


def test_analysis_run_v2_dataset_exclusions_not_in_spec_warn(tmp_path):
    state = _v2_state(tmp_path)
    state.excluded_samples = {"C0_r0": "x"}
    state.excluded_features = {"f1": "y"}
    out = json.loads(R.dataset_analysis_run(_v2_spec(), output_path=str(tmp_path / "o.json")))
    warn = [w for w in out["warnings"] if w["code"] == "dataset_exclusions_not_in_spec"]
    assert warn == [{"code": "dataset_exclusions_not_in_spec", "n_samples": 1, "n_features": 1}]
    assert out["analyses"][0]["n_rows"] == 10            # ds の除外は使わない
    listed = json.loads(R.dataset_analysis_run(_v2_spec(input={
        "exclude_samples": ["C0_r0"], "exclude_features": ["f1"]}), output_path=str(tmp_path / "p.json")))
    assert not [w for w in listed["warnings"] if w["code"] == "dataset_exclusions_not_in_spec"]
    assert listed["analyses"][0]["n_rows"] == 9


def test_analysis_run_v2_exclusions_compare_names_and_ids_on_one_key(tmp_path):
    """spec は試料を表示名でも試料 ID でも書ける（名前が先）。`dataset_exclude` の除外（表示名）と同じ鍵で
    比べるので、ID で書いた同じ試料を「spec に無い除外」と数えない。"""
    state = _v2_state(tmp_path)
    assert state.sample_ids[state.sample_names.index("C0_r0")] == "id_C0_r0"
    state.excluded_samples = {"C0_r0": "x", "T1_r1": "y"}
    by_id = json.loads(R.dataset_analysis_run(_v2_spec(input={"exclude_samples": ["id_C0_r0", "T1_r1"]}),
                                              output_path=str(tmp_path / "o.json")))
    assert by_id["status"] == "success"
    assert not [w for w in by_id["warnings"] if w["code"] == "dataset_exclusions_not_in_spec"]
    partial = json.loads(R.dataset_analysis_run(_v2_spec(input={"exclude_samples": ["id_C0_r0"]}),
                                                output_path=str(tmp_path / "p.json")))
    warn = [w for w in partial["warnings"] if w["code"] == "dataset_exclusions_not_in_spec"]
    assert warn == [{"code": "dataset_exclusions_not_in_spec", "n_samples": 1, "n_features": 0}]


def test_analysis_run_does_not_touch_session_or_reports(tmp_path):
    import copy
    state = _v2_state(tmp_path)
    before = (dict(state.__dict__.get("results", {}) or {}), copy.deepcopy(state.stat_matrices),
              state.sample_names[:], state.feature_matrix.copy())
    handle_before = session_state.session.current_handle_id
    R.dataset_analysis_run(_v2_spec(analyses=[{"analysis_id": "cmp", "kind": "compare"},
                                              {"analysis_id": "pc", "kind": "pca"}], primary="cmp"))
    after = (dict(state.__dict__.get("results", {}) or {}), copy.deepcopy(state.stat_matrices),
             state.sample_names[:], state.feature_matrix)
    assert after[:3] == before[:3] and (after[3] == before[3]).all()
    assert session_state.session.dataset is state
    assert session_state.session.current_handle_id == handle_before
    assert not (tmp_path / "analysis" / "reports").exists()


# ---------- 2026-11-13 の事実確認（担当 1: MS/MS とライブラリ） ----------

def _annotate_spot(monkeypatch, spot_id, name, inchikey=None):
    """`.arf2` の代表の注釈（MS-DIAL の名前・InChIKey）を差し替える（library フィクスチャの後に呼ぶ）。"""
    records = [{"AlignmentID": int(f["id"]), "MassCenter": f["mz"], "RT": f["rt"],
                "IonMode": "Positive", "Name": "Unknown"} for f in RELATION_FEATURES]
    records[[r["AlignmentID"] for r in records].index(spot_id)].update(Name=name, InChIKey=inchikey)
    monkeypatch.setattr(alignment_spot, "load_catalog", lambda _path: records)


def test_msms_qualify_reports_window_delta_ppm_and_msdial_agreement(data, library, monkeypatch):
    _annotate_spot(monkeypatch, 20, "GABA", GABA_KEY)
    handle = _open(data, msp=library)
    out = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=20,
                                            dataset_handle=handle, mz_tol_ppm=5))
    assert out["verdict"] == "qualified"      # 既定でも断片 1 本（87・69 の 2 本）を求める
    assert out["rules"]["min_fragments_excl_precursor"] == 1
    assert out["target_rank"] == 1 and out["top1_is_target"] is True
    assert out["query"]["mz_tol"] == pytest.approx(0.01)
    assert out["query"]["mz_window"] == pytest.approx(0.01)   # m/z 500 以下は Da のまま
    assert out["query"]["mz_tol_ppm"] == 5
    assert out["msdial_name"] == "GABA" and out["msdial_rank"] == 1
    assert out["top_differs_from_msdial"] is False
    lines = out["candidates"].splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert {"inchikey", "delta_ppm"} <= set(header) and "ppm_error" not in header
    assert row["inchikey"] == GABA_KEY and float(row["delta_ppm"]) == pytest.approx(0.0)
    assert row["fragment_matches_excl_precursor"] == "2"
    assert "passes_dbs_cutoffs" not in header             # .msp には注釈の足切りが無い


def test_msms_qualify_flags_a_top_hit_that_differs_from_msdial(data, library, monkeypatch):
    _annotate_spot(monkeypatch, 20, "Glutamine")
    handle = _open(data, msp=library)
    out = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=20,
                                            dataset_handle=handle))
    assert out["msdial_name"] == "Glutamine" and out["msdial_rank"] is None
    assert out["top_differs_from_msdial"] is True
    _annotate_spot(monkeypatch, 20, "Unknown")
    plain = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=20,
                                              dataset_handle=handle))
    assert "msdial_name" not in plain and "top_differs_from_msdial" not in plain


def test_msms_qualify_reads_short_ion_mode_and_rejects_unreadable(data, library):
    handle = _open(data, msp=library)
    rules = {"target": {"inchikey": GABA_KEY}}
    ok = json.loads(R.library_msms_qualify(rules, precursor_mz=104.0706, rt=2.0, ion_mode="pos",
                                           dataset_handle=handle))
    assert ok["query"]["ion_mode"] == "positive" and ok["target_rank"] == 1
    neg = json.loads(R.library_msms_qualify(rules, precursor_mz=104.0706, rt=2.0, ion_mode="-",
                                            dataset_handle=handle))
    assert neg["status"] == "no_candidates"               # 陽イオンの記録は陰イオンの照合に出ない
    bad = json.loads(R.library_msms_qualify(rules, precursor_mz=104.0706, ion_mode="sideways",
                                            dataset_handle=handle))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"
    assert bad["error"]["details"]["field"] == "ion_mode"


@pytest.mark.parametrize("bad", [0, -1, float("nan"), True, "5"])
def test_msms_qualify_rejects_a_bad_mz_tol_ppm(data, library, bad):
    out = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=20,
                                            dataset_handle=_open(data, msp=library),
                                            mz_tol_ppm=bad))
    assert out["error"]["details"]["field"] == "mz_tol_ppm"


def test_claim_evidence_msms_takes_mz_tol_ppm_and_shows_the_window(data, library, monkeypatch):
    _annotate_spot(monkeypatch, 20, "GABA", GABA_KEY)
    handle = _open(data, msp=library)
    claim = {**CLAIM, "msms": {"top_n": 3, "rules": {}, "mz_tol_ppm": 5}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    assert out["rule"]["msms"]["mz_tol_ppm"] == 5
    msms = {c["feature_id"]: c["msms"] for c in out["candidates"]}["20"]
    assert msms["window"] == {"mz_tol": 0.01, "mz_window": 0.01, "mz_tol_ppm": 5.0}
    assert msms["qualify"]["top1_is_target"] is True
    assert msms["msdial_name"] == "GABA" and msms["top_differs_from_msdial"] is False
    top = msms["top"][0]
    assert top["inchikey"] == GABA_KEY and top["delta_ppm"] == pytest.approx(0.0)
    # 省けば規則に載らない（rule_hash を変えない）
    plain = json.loads(R.dataset_claim_evidence({**CLAIM, "msms": {"top_n": 3, "rules": {}}},
                                                dataset_handle=handle))
    assert "mz_tol_ppm" not in plain["rule"]["msms"]
    assert {c["feature_id"]: c["msms"] for c in plain["candidates"]}["20"]["window"][
        "mz_tol_ppm"] is None
    bad = json.loads(R.dataset_claim_evidence(
        {**CLAIM, "msms": {"rules": {}, "mz_tol_ppm": -1}}, dataset_handle=handle))
    assert bad["error"]["code"] == "CLAIM_SPEC_INVALID"


def test_diagnostic_ions_neutral_loss_floor_and_key_named_errors(data, library):
    handle = _open(data, msp=library)
    out = json.loads(dcl_tools.dcl_diagnostic_ions(
        [20], [{"neutral_loss": 17.0265, "label": "-NH3", "min_relative_intensity": 0.5},
               {"neutral_loss": 35.0371, "label": "-NH3-H2O", "min_relative_intensity": 0.5}],
        dataset_handle=handle))
    rows = {r.split("\t")[1]: r.split("\t") for r in out["ions"].splitlines()[1:]}
    assert rows["-NH3"][4] == "true"                       # 87.0441 は基準ピーク
    assert rows["-NH3-H2O"][4] == "false"                  # 69.0335 は 0.4 しかない
    bad = json.loads(dcl_tools.dcl_diagnostic_ions(
        [20], [{"neutral_loss": 17.0265, "min_intensity": 0.5}], dataset_handle=handle))
    assert bad["error"]["code"] == "MSMS_RULES_INVALID"
    assert "min_intensity" in bad["error"]["message"] and "[0][0]" not in bad["error"]["message"]
    assert bad["error"]["details"]["field"] == "ions[0]"


def test_msms_qualify_shows_whether_each_candidate_passes_the_dbs_cutoffs(data, library, tmp_path):
    from tests.dbs_fixture import record, write_dbs

    dbs = write_dbs(tmp_path / "lab" / "lib.msp2.dbs", [("lib", [
        record(ion_mode=0), record(name="Decoy", inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N",
                                   ion_mode=0, peaks=((30.0, 999.0),))], ["k"])])
    handle = _open(data, msp=dbs)
    out = json.loads(R.library_msms_qualify({"target": {"inchikey": GABA_KEY}}, spot_id=20,
                                            dataset_handle=handle))
    assert out["dbs_cutoffs"]["total_score_cutoff"] == pytest.approx(0.8)
    assert out["dbs_cutoffs"]["minimum_spectrum_match"] == pytest.approx(1.0)
    lines = out["candidates"].splitlines()
    header = lines[0].split("\t")
    rows = {r.split("\t")[1]: dict(zip(header, r.split("\t"))) for r in lines[1:]}
    assert header[-1] == "passes_dbs_cutoffs"
    assert rows["GABA"]["passes_dbs_cutoffs"] == "true"
    assert rows["Decoy"]["passes_dbs_cutoffs"] == "false"
    assert out["verdict"] == "qualified"                  # 足切りは材料。合否には当てない


def test_dbs_cutoffs_compare_squared_dots_and_skip_the_sentinel():
    from metabolomix.library import matching

    assert matching.dbs_cutoffs({"ms1_tolerance": 0.01}) is None
    cutoffs = matching.dbs_cutoffs({"squared_weighted_dot_cutoff": 0.25, "total_score_cutoff": 1.0,
                                    "minimum_spectrum_match": 2.0})
    assert set(cutoffs) == {"squared_weighted_dot_cutoff", "total_score_cutoff",
                            "minimum_spectrum_match"}
    good = {"weighted_dot_product": 0.6, "total_score": 1.5, "matched_peaks_count": 3}
    assert matching.passes_dbs_cutoffs(good, cutoffs) is True
    # weighted dot 0.45 の二乗は 0.2025 < 0.25
    assert matching.passes_dbs_cutoffs({**good, "weighted_dot_product": 0.45}, cutoffs) is False
    assert matching.passes_dbs_cutoffs({**good, "weighted_dot_product": -1}, cutoffs) is False
    assert matching.passes_dbs_cutoffs({**good, "matched_peaks_count": 1}, cutoffs) is False
    assert matching.passes_dbs_cutoffs(good, None) is None


# ---------- dataset_homolog_check の auto は付加体を組成で比べる（最終検証 F26） ----------

#: mzTab-M の SME の付加体は `[M+H]1+`。RT は炭素数が増えるほど早い（逆相の向きに反する）。
HOMOLOG_FEATURES = tuple(
    {"id": str(30 + i), "mz": mz, "rt": rt, "values": _scaled(10),
     "sme": [{"name": name, "inchikey": key, "formula": formula, "adduct": "[M+H]1+", "theo": mz}]}
    for i, (name, key, formula, mz, rt) in enumerate((
        ("PC 32:1", "AAAAAAAAAAAAAA-UHFFFAOYSA-N", "C40H78NO8P", 732.5538, 5.0),
        ("PC 34:1", "BBBBBBBBBBBBBB-UHFFFAOYSA-N", "C42H82NO8P", 760.5851, 4.0),
        ("PC 36:1", "CCCCCCCCCCCCCC-UHFFFAOYSA-N", "C44H86NO8P", 788.6164, 3.0),
    )))


def test_homolog_check_auto_adduct_matches_the_mztab_spelling(tmp_path):
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=HOMOLOG_FEATURES)
    handle = _open({"mztab": mztab})
    for adduct in ("[M+H]1+", "[M+H]+"):
        out = json.loads(R.dataset_homolog_check(
            {"auto": {"lipid_class": "PC", "adduct": adduct}}, dataset_handle=handle))
        assert out["verdict"] == "fail" and out["n_violations"] == 2, out
        assert out["n_members"] == 3 and out["n_pairs"] == 2 and "notes" not in out
    absent = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PC", "adduct": "[M+Na]+"}}, dataset_handle=handle))
    assert absent["verdict"] == "not_evaluable" and absent["n_members"] == 0 and absent["n_pairs"] == 0
    [note] = absent["notes"]
    assert "not_evaluable" in note and "0 件" in note


def test_homolog_check_auto_reads_noted_sme_names_and_skips_internal_standards(tmp_path):
    """Console / mzTab-M の SME の名前は MSP の記録名（注記付き）のまま。同族列も目標照合・MS/MS の合否と同じ読み手で
    注記を外して読み、RT の規則違反を見つける。同位体の標識のある内部標準は天然の種の列に混ぜない（レビュー V2）。"""
    noted = tuple({**f, "sme": [{**f["sme"][0], "name": f"{f['sme'][0]['name']}; ExampleDB ID-000{i}"}]}
                  for i, f in enumerate(HOMOLOG_FEATURES, 1))
    standard = {"id": "40", "mz": 767.6290, "rt": 9.0, "values": _scaled(10),
                "sme": [{"name": "PC 15:0_18:1(d7)", "inchikey": "DDDDDDDDDDDDDD-UHFFFAOYSA-N",
                         "formula": "C41H73D7NO8P", "adduct": "[M+H]1+", "theo": 767.6290}]}
    mztab = write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=noted + (standard,))
    handle = _open({"mztab": mztab})
    out = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PC", "adduct": "[M+H]+"}}, dataset_handle=handle))
    assert out["verdict"] == "fail" and out["n_violations"] == 2, out
    assert out["n_members"] == 3 and out["n_pairs"] == 2
    header, *rows = out["members"].splitlines()
    assert sorted(r.split("	")[0] for r in rows) == ["30", "31", "32"]     # 内部標準（40）は集めない
    # P- は種のレベルで O- に直るので、P- のクラスは集まらず、notes がそう書く（断らない）
    plasmalogen = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PC P-"}}, dataset_handle=handle))
    assert plasmalogen["verdict"] == "not_evaluable" and plasmalogen["n_members"] == 0
    assert any("'PC O-'" in note for note in plasmalogen["notes"])
    # 修飾を書いた P- の指定（"PC P-;O"）も同じ案内で、案内のクラスに修飾を付け直す（最終の確かめ: 案内が消えていた）
    oxidized = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PC P-;O"}}, dataset_handle=handle))
    assert oxidized["verdict"] == "not_evaluable" and oxidized["n_members"] == 0
    assert any("'PC O-;O'" in note for note in oxidized["notes"]), oxidized["notes"]


def test_homolog_check_auto_keeps_oxidized_species_in_their_own_series(tmp_path):
    """SME の `PE 36:2;O (oxPE)` は PE 36:2 の重複ではない（最終検証 WB）。混ぜると 36:2 の RT の中央値が
    (12.0 + 9.0) / 2 = 10.5 < 36:1 の 11.0 に下がって違反が消えていた。修飾ごとの列で比べて違反を残し、
    集めた修飾の内訳を notes に、`members` に `modification` の列を書く。rule_hash は spec だけで決まる（変わらない）。"""
    rows = (("PE 18:0_18:2", "EEEEEEEEEEEEEE-UHFFFAOYSA-N", "C41H78NO8P", 744.5538, 12.0),
            ("PE 36:2;O (oxPE)", "FFFFFFFFFFFFFF-UHFFFAOYSA-N", "C41H78NO9P", 760.5487, 9.0),
            ("PE 18:0_18:1", "GGGGGGGGGGGGGG-UHFFFAOYSA-N", "C41H80NO8P", 746.5694, 11.0))
    features = tuple({"id": str(50 + i), "mz": mz, "rt": rt, "values": _scaled(10),
                      "sme": [{"name": name, "inchikey": key, "formula": formula, "adduct": "[M+H]1+",
                               "theo": mz}]}
                     for i, (name, key, formula, mz, rt) in enumerate(rows))
    handle = _open({"mztab": write_mztab(tmp_path / "run" / "AlignResult-1.mzTab", features=features)})
    spec = {"auto": {"lipid_class": "PE", "adduct": "[M+H]+"}}
    out = json.loads(R.dataset_homolog_check(spec, dataset_handle=handle))
    assert out["verdict"] == "fail" and out["n_violations"] == 1 and out["n_pairs"] == 1, out
    assert out["duplicates"] == {}
    [violation] = _tsv_rows(out["violations"])
    assert (violation["low"], violation["high"]) == ("36:1", "36:2")
    members = {r["feature_id"]: r for r in _tsv_rows(out["members"])}
    assert members["51"]["modification"] == "O" and members["51"]["status"] == "ok"
    assert members["50"]["modification"] == "" and members["50"]["status"] == "violation"
    [note] = out["notes"]
    assert "無修飾 2 件" in note and ";O 1 件" in note and "'PE;O'" in note
    from metabolomix.analysis import homolog
    assert out["rule_hash"] == homolog.rule_hash(homolog.normalize_spec(spec))
    # 修飾を書けばその列だけ: 1 種なので not_evaluable（修飾の内訳の注記は出ない）
    oxidized = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PE;O"}}, dataset_handle=handle))
    assert oxidized["verdict"] == "not_evaluable" and oxidized["n_members"] == 1
    assert len(oxidized["notes"]) == 1 and "not_evaluable" in oxidized["notes"][0]
    # 修飾の無い列は前と同じ形（modification の列を足さない）
    plain = json.loads(R.dataset_homolog_check(
        {"auto": {"lipid_class": "PE;"}}, dataset_handle=handle))
    assert plain["members"].splitlines()[0].split("\t") == ["feature_id", "carbons", "double_bonds", "rt",
                                                             "status"]
