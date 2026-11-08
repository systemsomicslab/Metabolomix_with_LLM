"""再解析のツール（spec 2026-10-22 §4 の 3〜11）を MCP のツールとして呼ぶ。

ハンドルの上でも対話のセッションでも動き、どちらの状態も書き換えないこと・同じ spec なら
同じ JSON になることを確かめる。fixture はすべてテスト自身が作る。
"""
from __future__ import annotations
from tests.v1_chain import v1_preprocess

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

def test_analysis_run_on_handle_returns_rows_and_leaves_states(data):
    handle = _open(data)
    out = json.loads(R.dataset_analysis_run(_spec(data), dataset_handle=handle))
    assert out["status"] == "success" and out["comparison"] == "RA_vs_HC"
    rows = out["analyses"][0]["rows"].splitlines()
    assert rows[0].startswith("feature_id\tstatus\treason\tgroup_ref\tgroup_test")
    assert rows[1].split("\t")[:5] == ["20", "tested", "", "HC", "RA"]
    assert session_state.session.dataset is None
    assert handles.registry.get(handle).dataset.pp_matrix is None


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
    assert header[-4:] == ["detected_reference", "detected_test", "detection_fisher_p",
                           "detection_driven"]
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
    v1_preprocess(dataset_handle=handle)
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


def test_novel_candidates_count_split_once(tmp_path):
    handle = _split_handle(tmp_path)
    v1_preprocess(dataset_handle=handle)
    diff = json.loads(T.dataset_differential(group_a=["HC1", "HC2", "HC3"],
                                             group_b=["RA1", "RA2", "RA3"], dataset_handle=handle))
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
    v1_preprocess(dataset_handle=handle)
    diff = json.loads(T.dataset_differential(group_a=["HC1", "HC2", "HC3"],
                                             group_b=["RA1", "RA2", "RA3"], dataset_handle=handle))
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
    v1_preprocess(dataset_handle=handle)
    diff = json.loads(T.dataset_differential(group_a=["HC1", "HC2", "HC3"],
                                             group_b=["RA1", "RA2", "RA3"], dataset_handle=handle))
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

#: ppm_offset を足す前（v8）の実装で CLAIM を処理した rule_hash
_CLAIM_RULE_HASH_V8 = "41f68e9459154e44cd3fe48101cdb088bc4969e917940bcfa3fb25b7410b524c"


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
