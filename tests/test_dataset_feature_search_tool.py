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
    assert set(parsed) == {"targets", "candidates", "stats", "results", "rule",
                           "stats_source", "notes"}

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

    assert parsed["stats"].splitlines()[0] == "\t".join(STATS_COLUMNS)
    [stat] = _tsv(parsed["stats"])
    assert stat["feature_id"] == "f5" and stat["result_id"] == rid
    assert stat["status"] == "tested" and stat["p"] != ""

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
    assert set(parsed) == {"features", "stats", "results", "stats_source", "notes"}
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
