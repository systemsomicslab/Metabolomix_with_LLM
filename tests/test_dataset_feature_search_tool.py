"""dataset_feature_search（MCP ツール層）のテスト。

純ロジックは tests/test_feature_search.py が固定している。ここでは
ツール層の仕事——前提状態の封筒・引数の取り違えの拒否・極性の出所・
差次的結果の選び方（名指しの古い結果は missing_state）・丸めと TSV・
注記——を見る。fixture はテスト自身が作る。
"""
import asyncio
import json
from types import SimpleNamespace

import numpy as np
import pytest

from metabolomix.analysis.feature_search import CANDIDATE_COLUMNS, STATS_COLUMNS, WINDOW_COLUMNS
from metabolomix.analysis.result_state import dataset_fingerprint
from metabolomix.core import session_state
from metabolomix.core.mcp_errors import MISSING_STATE
from metabolomix.mztab.dataset_state import DatasetState

GABA_KEY = "BTCSSZJGUNDROE-UHFFFAOYSA-N"
OTHER_KEY = "QNAYBMKLOCPYGJ-REOHCLBHSA-N"


@pytest.fixture(autouse=True)
def reset_session():
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def _load_ds(n_features=20, n_samples=8):
    rng = np.random.default_rng(0)
    ds = DatasetState()
    ds.feature_matrix = rng.random((n_features, n_samples)) * 1000.0
    ds.sample_names = ([f"ctrl_{i}" for i in range(n_samples // 2)]
                       + [f"treat_{i}" for i in range(n_samples - n_samples // 2)])
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    ds.feature_metadata = {
        f"f{i}": {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                  "inchikey": None, "inchikey_source": None}
        for i in range(n_features)
    }
    # f5 だけ SME 候補（MS/MS 同定）を持つ
    ds.feature_candidates = {"f5": [{"rank": 1, "chemical_name": "GABA", "inchikey": GABA_KEY,
                                     "adduct": "[M+H]1+", "charge": 1}]}
    ds.validation_result = {"ok": True, "errors": [], "warnings": []}
    session_state.session.dataset = ds
    return ds


def _differential(normalize="none"):
    """前処理 → 2 群比較をツール経由で回し、result_id を返す。

    同じレシピの前処理は同じ前処理として扱われるので、古い結果を作るときは
    `normalize` を変えて前処理をやり直す。
    """
    from metabolomix.tools.dataset_analysis_tools import dataset_differential, dataset_preprocess
    dataset_preprocess(normalize=normalize)
    ds = session_state.session.dataset
    a = [n for n in ds.pp_sample_names if n.startswith("ctrl")]
    b = [n for n in ds.pp_sample_names if n.startswith("treat")]
    return json.loads(dataset_differential(group_a=a, group_b=b))["result_id"]


def _target(**extra):
    target = {"target_id": "c1", "inchikey": GABA_KEY, "expected_mz": 105.0,
              "adduct": "[M+H]+", "rt_min": 1.5}
    target.update(extra)
    return target


def _search(**kw):
    from metabolomix.tools.dataset_analysis_tools import dataset_feature_search
    return json.loads(dataset_feature_search(**kw))


def _tsv(text):
    lines = text.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def test_without_dataset_returns_missing_state():
    parsed = _search(min_mz=100.0)
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["state"] == "dataset"
    assert parsed["error"]["required_tools"] == ["dataset_load"]


def test_window_and_targets_together_is_bad_request():
    _load_ds()
    parsed = _search(min_mz=100.0, targets=[_target()])
    assert parsed["error"]["code"] == "DATASET_BAD_REQUEST"
    # 引数の取り違えは別のツールを先に呼んでも直らない
    assert parsed["error"]["code"] != MISSING_STATE


def test_target_mode_payload_shape_and_rule_hash():
    _load_ds()
    rid = _differential()
    parsed = _search(targets=[_target()])
    assert set(parsed) == {"targets", "candidates", "stats", "stats_unmerged", "merge_split",
                           "results", "rule", "stats_source", "notes"}

    [target] = parsed["targets"]
    assert target["status"] == "resolved" and target["reason"] is None
    assert target["selected_feature_id"] == "f5"
    assert target["selected_adduct"] == "[M+H]+"
    assert target["selected_evidence"] == "sme_match"
    assert target["truncated"] == 0
    assert "candidates" not in target                     # 候補は TSV の側に 1 回だけ

    assert parsed["candidates"].splitlines()[0] == "\t".join(CANDIDATE_COLUMNS)
    [row] = _tsv(parsed["candidates"])
    assert row["feature_id"] == "f5" and row["evidence"] == "sme_match"
    assert row["ppm"] == "0.0" and row["mz"] == "105.0"

    # merge_split（既定 true）は末尾に merged_from を足す。割れの無い特徴では空で、値は従来どおり
    assert parsed["stats"].splitlines()[0] == "\t".join(STATS_COLUMNS + ("merged_from",))
    [stat] = _tsv(parsed["stats"])
    assert stat["feature_id"] == "f5" and stat["result_id"] == rid
    assert stat["status"] == "tested" and stat["p"] != "" and stat["merged_from"] == ""
    assert _tsv(parsed["stats_unmerged"]) == [{k: v for k, v in stat.items() if k != "merged_from"}]
    assert parsed["merge_split"]["groups"] == [] and parsed["merge_split"]["recomputed"] == []

    rule = parsed["rule"]
    assert len(rule["rule_hash"]) == 64
    assert rule["polarity"] is None and rule["polarity_source"] == "unknown"
    assert rule["targets"][0]["target_id"] == "c1"
    assert parsed["stats_source"] == {
        "result_ids": [rid], "dataset_hash": dataset_fingerprint(session_state.session.dataset)}
    [summary] = parsed["results"]
    assert summary["result_id"] == rid
    assert summary["bh_population"]["n_features_total"] == 20

    # 既定値を明示しても、別の呼び出しでも、規則の hash は同じ（result_id を混ぜない）
    again = _search(targets=[_target(ppm=10, rt_tol_min=0.3)])
    assert again["rule"]["rule_hash"] == rule["rule_hash"]


def test_window_mode_payload_shape_and_rounding():
    ds = _load_ds()
    ds.feature_metadata["f1"]["mz"] = 101.123456789
    ds.feature_metadata["f1"]["rt"] = 1.123456789
    _differential()
    parsed = _search(min_mz=100.0, max_mz=102.0)
    assert set(parsed) == {"features", "stats", "stats_unmerged", "merge_split", "results",
                           "stats_source", "notes"}
    assert parsed["features"].splitlines()[0] == "\t".join(WINDOW_COLUMNS)
    rows = _tsv(parsed["features"])
    assert [r["feature_id"] for r in rows] == ["f0", "f1", "f2"]
    assert rows[1]["mz"] == "101.1235" and rows[1]["rt"] == "1.123"
    stats = _tsv(parsed["stats"])
    assert [s["feature_id"] for s in stats] == ["f0", "f1", "f2"]
    for s in stats:
        # p・q・平均は有効 4 桁
        for key in ("p", "q", "mean_a", "mean_b"):
            mantissa = s[key].split("e")[0].replace("-", "").replace(".", "").lstrip("0")
            assert 1 <= len(mantissa) <= 4, (key, s[key])
    # 検出マスクが無いので検出例数は空で、その理由が注記される
    assert all(s["detected_a"] == "" for s in stats)
    assert any("検出" in note for note in parsed["notes"])


def test_window_truncation_is_noted():
    _load_ds()
    parsed = _search(min_mz=100.0, max_results=2)
    assert len(_tsv(parsed["features"])) == 2
    assert any("18" in note for note in parsed["notes"])


def test_empty_feature_ids_is_a_condition_with_zero_hits_and_noted():
    _load_ds()
    parsed = _search(feature_ids=[])
    assert parsed["features"] == ""
    assert any("feature_ids" in note for note in parsed["notes"])


def test_named_stale_result_returns_missing_state():
    _load_ds()
    old = _differential()
    _differential(normalize="median")                 # 前処理のやり直しで old は古くなる
    parsed = _search(min_mz=100.0, result_ids=[old])
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["state"] == "dataset_differential_result"
    assert parsed["error"]["required_tools"] == ["dataset_differential"]
    assert old in parsed["error"]["message"]

    missing = _search(min_mz=100.0, result_ids=["res_nope"])
    assert missing["error"]["code"] == MISSING_STATE


def test_omitted_result_ids_skip_stale_results_and_note_them():
    _load_ds()
    old = _differential()
    new = _differential(normalize="median")
    parsed = _search(min_mz=100.0, max_mz=101.0)
    assert parsed["stats_source"]["result_ids"] == [new]
    assert any(old in note for note in parsed["notes"])


def test_named_result_ids_are_deduplicated():
    _load_ds()
    rid = _differential()
    parsed = _search(min_mz=100.0, max_mz=100.5, result_ids=[rid, rid])
    assert parsed["stats_source"]["result_ids"] == [rid]
    assert len(_tsv(parsed["stats"])) == 1


def test_no_differential_result_returns_without_stats_and_notes_it():
    _load_ds()
    parsed = _search(min_mz=100.0, max_mz=102.0)
    assert len(_tsv(parsed["features"])) == 3
    assert parsed["stats"] == "" and parsed["results"] == []
    assert parsed["stats_source"]["result_ids"] == []
    assert any("差次的" in note for note in parsed["notes"])


def _write_negative_mztab(path):
    path.write_text(
        "MTD\tmzTab-version\t2.0.0-M\n"
        "SMH\tSML_ID\tadduct_ions\n"
        "SML\t1\t[M-H]1-\n"
        "SML\t2\t[M-H]1-\n",
        encoding="utf-8")


def test_polarity_from_job_declaration_is_preferred(tmp_path, monkeypatch):
    from metabolomix.console import job_manager
    ds = _load_ds()
    mztab = tmp_path / "align.mzTab"
    _write_negative_mztab(mztab)
    ds.source_files = {str(mztab): "sha"}
    ds.job_path = str(tmp_path / "analysis-job.json")
    monkeypatch.setattr(job_manager, "load_job", lambda path: SimpleNamespace(polarity="positive"))

    parsed = _search(targets=[{"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2"}])
    assert parsed["rule"]["polarity"] == "positive"
    assert parsed["rule"]["polarity_source"] == "job_declared"
    assert parsed["rule"]["targets"][0]["adduct_tiers"][0] == ["[M+H]+", "[M]+"]
    # 宣言と中身の多数決が食い違うことは注記する
    assert any("negative" in note for note in parsed["notes"])


def test_polarity_falls_back_to_adduct_majority(tmp_path):
    ds = _load_ds()
    mztab = tmp_path / "align.mzTab"
    _write_negative_mztab(mztab)
    ds.source_files = {str(mztab): "sha"}
    parsed = _search(targets=[{"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2"}])
    assert parsed["rule"]["polarity"] == "negative"
    assert parsed["rule"]["polarity_source"] == "adduct_majority"


def test_unknown_polarity_without_tiers_is_an_error_envelope():
    _load_ds()
    parsed = _search(targets=[{"target_id": "c1", "inchikey": GABA_KEY, "formula": "C4H9NO2"}])
    assert parsed["error"]["code"] == "POLARITY_UNKNOWN"


def test_invalid_target_is_an_error_envelope_not_an_exception():
    _load_ds()
    parsed = _search(targets=[{"target_id": "c1", "formula": "C4H9NO2"}])   # inchikey が無い
    assert parsed["error"]["code"] == "FEATURE_TARGET_INVALID"
    assert parsed["error"]["details"]["field"] == "inchikey"


def test_empty_window_query_is_an_error_envelope():
    _load_ds()
    parsed = _search()
    assert parsed["error"]["code"] == "FEATURE_SEARCH_EMPTY_QUERY"


def test_multiple_candidates_need_input_and_name_who_checks_next():
    ds = _load_ds()
    ds.feature_metadata["f6"]["mz"] = 105.0       # 同じ m/z・RT 窓に 2 つの未注釈の特徴
    ds.feature_metadata["f6"]["rt"] = 1.5
    ds.feature_candidates = {}
    parsed = _search(targets=[_target()])
    [target] = parsed["targets"]
    assert target["status"] == "needs_input"
    assert target["reason"] == "multiple_qualified_candidates"
    assert target["selected_feature_id"] is None
    assert target["selected_evidence"] == "unannotated"
    assert len(_tsv(parsed["candidates"])) == 2


def test_curation_flags_are_noted_but_not_applied(tmp_path, monkeypatch):
    from metabolomix.curation import apply as curation_apply
    ds = _load_ds()
    ds.source_files = {str(tmp_path / "align.mzTab"): "sha"}
    monkeypatch.setattr(curation_apply, "arf2_for_mztab", lambda path: tmp_path / "align.arf2")
    monkeypatch.setattr(curation_apply, "flags_for_arf2", lambda path: {"n": 3})
    parsed = _search(min_mz=100.0, max_mz=101.0)
    assert any("キュレーション" in note for note in parsed["notes"])


def test_unreadable_curation_flags_do_not_fail_the_search(tmp_path, monkeypatch):
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError
    ds = _load_ds()
    ds.source_files = {str(tmp_path / "align.mzTab"): "sha"}
    monkeypatch.setattr(curation_apply, "arf2_for_mztab", lambda path: tmp_path / "align.arf2")

    def _broken(path):
        raise FlagFileError(path, 3, "JSON として読めない")
    monkeypatch.setattr(curation_apply, "flags_for_arf2", _broken)
    parsed = _search(min_mz=100.0, max_mz=101.0)
    assert "error" not in parsed
    assert any("キュレーション" in note for note in parsed["notes"])


def test_payload_has_no_structured_duplicate():
    import server
    tools = asyncio.run(server.mcp.list_tools())
    tool = next(item for item in tools if item.name == "dataset_feature_search")
    assert tool.outputSchema is None              # content と structuredContent の二重送出をしない
    assert tool.annotations.readOnlyHint is True


def test_dataset_fingerprint_is_computed_once_per_call(monkeypatch):
    """データの指紋は行列全体の hash。引き継いだ結果ごとに計算し直さず、1 回で使い回す。"""
    import copy

    from metabolomix.analysis import result_state

    _load_ds()
    session_rid = _differential()
    ds = session_state.session.dataset
    handed = []
    for i in range(3):
        result = copy.deepcopy(ds.results[session_rid])
        rid = f"res_handed_{i}"
        result["provenance"].update(result_id=rid, dataset_id="ds_run",
                                    handoff={"pipeline_root": "x", "comparison_id": f"c{i}"})
        ds.results[rid] = result
        handed.append(rid)

    calls = []
    original = result_state.dataset_fingerprint

    def counting(target):
        calls.append(target)
        return original(target)

    monkeypatch.setattr(result_state, "dataset_fingerprint", counting)
    parsed = _search(min_mz=100.0, max_mz=101.0)

    assert len(calls) == 1
    assert parsed["stats_source"]["result_ids"] == [session_rid, *handed]
    assert parsed["stats_source"]["dataset_hash"] == original(ds)


# ---------- 割れたスポットのまとめ（merge_split。spec 2026-10-28 §4.3 (c)） ----------

CE_KEY = "RJECHNNFRHZQKU-RMUVNZEASA-N"

#: 最初の 2 注入（ファイル名の順で HC1・HC2）だけが RT +0.6 の別スポット 31 に割れた CE。
#: 32〜37 は割れの無い特徴（BH の母集団を作る）。
SPLIT = (
    {"id": "30", "mz": 650.6000, "rt": 20.00, "values": [20, 21, 1000, 1010, 1020, 1030, 990, 995],
     "sme": [{"name": "CE 18:2", "inchikey": CE_KEY, "formula": "C45H76O2", "adduct": "[M+H]1+",
              "theo": 650.6000}]},
    {"id": "31", "mz": 650.6040, "rt": 20.60, "values": [1005, 1015, 22, 20, 21, 23, 20, 22],
     "sme": []},
) + tuple(
    {"id": str(32 + i), "mz": 300.0 + 37.0 * i, "rt": 3.0 + i,
     "values": [50 + i, 61 + 2 * i, 55, 70 + 3 * i, 66, 90 - i, 60, 61], "sme": []}
    for i in range(6))


def _merged_values():
    by_id = {f["id"]: f for f in SPLIT}
    return [max(a, b) for a, b in zip(by_id["30"]["values"], by_id["31"]["values"])]


def _load_mztab(path, features):
    from metabolomix.tools import mztab_tools
    from tests.reanalysis_fixtures import write_mztab

    mztab_tools.dataset_load(mztab_path=str(write_mztab(path, features=features)))
    assert session_state.session.dataset is not None


def _differential_on_session():
    from metabolomix.tools.dataset_analysis_tools import dataset_differential, dataset_preprocess
    dataset_preprocess()
    out = json.loads(dataset_differential(group_a=["HC1", "HC2", "HC3"],
                                          group_b=["RA1", "RA2", "RA3"]))
    return out["result_id"]


def _oracle_stats(tmp_path):
    """30 の列を試料ごとの大きい方の値にし、31 を除いた行列を最初から読み込んで検定した 30 の行。"""
    merged = ({**SPLIT[0], "values": _merged_values()},) + SPLIT[2:]
    _load_mztab(tmp_path / "oracle" / "AlignResult-1.mzTab", merged)
    _differential_on_session()
    [row] = _tsv(_search(feature_ids=["30"], merge_split=False)["stats"])
    return row


def _split_session(tmp_path):
    _load_mztab(tmp_path / "split" / "AlignResult-1.mzTab", SPLIT)
    return _differential_on_session()


_NUMBERS = ("log2fc", "p", "q", "mean_a", "mean_b")


def test_feature_search_merge_split_recomputes_q(tmp_path):
    oracle = _oracle_stats(tmp_path)
    rid = _split_session(tmp_path)
    unmerged = {r["feature_id"]: r for r in _tsv(_search(feature_ids=["30", "31"],
                                                         merge_split=False)["stats"])}

    parsed = _search(feature_ids=["30", "31"])
    [row] = _tsv(parsed["stats"])                      # 2 スポットは 1 つの物質として 1 行
    assert parsed["stats"].splitlines()[0] == "\t".join(STATS_COLUMNS + ("merged_from",))
    assert row["feature_id"] == "30" and row["merged_from"] == "30;31"
    assert row["result_id"] == rid and row["status"] == "tested"
    # まとめた行列で最初から回した検定（BH はまとめた後の全特徴）と同じ数字
    assert {k: row[k] for k in _NUMBERS} == {k: oracle[k] for k in _NUMBERS}
    assert row["q"] != unmerged["30"]["q"]
    # まとめなかった場合の stats（従来の値）も同じ行の並びで残す
    assert {r["feature_id"]: r for r in _tsv(parsed["stats_unmerged"])} == unmerged
    assert parsed["stats_unmerged"].splitlines()[0] == "\t".join(STATS_COLUMNS)
    block = parsed["merge_split"]
    assert block["enabled"] is True and block["source"] == "find_split_pairs"
    assert block["groups"] == [["30", "31"]] and block["recomputed"] == [rid]
    assert block["bh_population"][rid]["n_features_total"] == len(SPLIT) - 1

    # 消える側の ID（31）を名指しても、残す ID の行に解く
    [only] = _tsv(_search(feature_ids=["31"])["stats"])
    assert only["feature_id"] == "30" and only["q"] == oracle["q"]

    # 目標照合でも、照合した特徴の stats はまとめた後の値
    target = {"target_id": "ce", "inchikey": CE_KEY, "expected_mz": 650.6, "adduct": "[M+H]+",
              "rt_min": 20.0}
    matched = _search(targets=[target])
    assert matched["targets"][0]["selected_feature_id"] == "30"
    [stat] = _tsv(matched["stats"])
    assert stat["merged_from"] == "30;31" and stat["q"] == oracle["q"]
    assert _tsv(matched["stats_unmerged"])[0]["q"] == unmerged["30"]["q"]


def test_feature_search_no_split_stats_unchanged(tmp_path, monkeypatch):
    from metabolomix.analysis import feature_search

    _split_session(tmp_path)
    plain = _tsv(_search(feature_ids=["32", "33"], merge_split=False)["stats"])

    def _must_not_run(*args, **kwargs):
        raise AssertionError("割れの組を含まない特徴では検定をやり直さない")

    monkeypatch.setattr(feature_search, "retest_merged", _must_not_run)
    parsed = _search(feature_ids=["32", "33"])
    rows = _tsv(parsed["stats"])
    assert [{k: v for k, v in r.items() if k != "merged_from"} for r in rows] == plain
    assert all(r["merged_from"] == "" for r in rows)
    assert _tsv(parsed["stats_unmerged"]) == plain
    assert parsed["merge_split"]["recomputed"] == [] and parsed["merge_split"]["groups"] == []


def test_feature_search_merge_split_false(tmp_path):
    _split_session(tmp_path)
    parsed = _search(feature_ids=["30", "31"], merge_split=False)
    assert "stats_unmerged" not in parsed
    assert parsed["merge_split"] == {"enabled": False}
    assert parsed["stats"].splitlines()[0] == "\t".join(STATS_COLUMNS)
    assert [r["feature_id"] for r in _tsv(parsed["stats"])] == ["30", "31"]

    # 効かない設定は書かせない・別のデータセットの JSON は拒む
    qc_path = tmp_path / "qc.json"
    qc_path.write_text(json.dumps({"schema": "alignment-qc.v1", "dataset_fingerprint": "other",
                                   "merge_groups": [{"features": ["30", "31"], "mergeable": True}]}),
                       encoding="utf-8")
    bad = _search(feature_ids=["30"], merge_split=False, alignment_qc_path=str(qc_path))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"
    foreign = _search(feature_ids=["30"], alignment_qc_path=str(qc_path))
    assert foreign["error"]["code"] == "ANALYSIS_INPUT_INVALID"

    # dataset_alignment_qc が書いた JSON を渡すと、その merge_groups でまとめる
    from metabolomix.tools import reanalysis_tools as R
    written = tmp_path / "alignment_qc_pos.json"
    R.dataset_alignment_qc(output_path=str(written))
    via_file = _search(feature_ids=["31"], alignment_qc_path=str(written))
    assert via_file["merge_split"]["source"] == "alignment_qc"
    assert _tsv(via_file["stats"])[0]["merged_from"] == "30;31"
