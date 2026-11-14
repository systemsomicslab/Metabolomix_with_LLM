"""`dataset_compare`（統計エンジンの比較ツール）。fixture はテスト自身が tmp に作る。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from metabolomix.analysis import stats_results
from metabolomix.analysis.dataset_exclusions import apply_exclusions
from metabolomix.analysis.stats import spec_v2
from metabolomix.analysis.stats.spec_v2 import normalize_spec_v2, run_spec_v2
from metabolomix.core import session_state
from metabolomix.tools import stats_tools
from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess
from tests.test_stats_results import _ds

FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")


def _install(state, tmp_path):
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


@pytest.fixture
def ds(tmp_path):
    return _install(_ds(), tmp_path)


def _pre(**kw) -> dict:
    return json.loads(dataset_preprocess(**{**FLAT_OFF, **kw}))


def _cmp(**kw) -> dict:
    return json.loads(dataset_compare(**kw))


def _strict(text: str) -> dict:
    return json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))


def _contrast(out, cid):
    return next(c for c in out["contrasts"] if c["contrast_id"] == cid)


# --- 群の指定 ---------------------------------------------------------------------------

def test_two_groups_default_moderated_t_summary_shape(ds):
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    assert out["status"] == "success"
    assert out["matrix_id"] == mx and out["result_id"].startswith("cmp_")
    assert out["test"] == "moderated_t"
    assert out["effect_size_definition"]["name"] == "log2_mean_difference"
    assert out["ci_level"] == pytest.approx(0.95)
    assert out["verdict_rule"]["q_threshold"] == 0.05
    assert out["bh_population"]["name"] == "matrix"
    assert out["bh_population"]["n_by_contrast"] == {"T_vs_C": 10}
    assert out["group_n"] == {"C": 4, "T": 4}
    assert out["prior"].keys() == {"s2_prior", "df_prior"}
    assert out["omnibus"] is None
    assert len(out["contrasts"]) == 1
    c = _contrast(out, "T_vs_C")
    assert (c["reference"], c["test"]) == ("C", "T")
    assert c["n_tested"] == 10 and c["n_not_tested"] == 0
    assert c["verdicts"]["changed"] >= 3
    assert sum(c["verdicts"].values()) == 10
    assert c["large_changes"] >= 3
    assert c["presence"]["available"] is True
    lines = out["top"].splitlines()
    assert lines[0].split("\t") == ["contrast_id", "feature_id", "name", "id_level", "log2fc",
                                    "ci_low", "ci_high", "q", "verdict"]
    assert 1 < len(lines) <= 11
    first = lines[1].split("\t")
    assert first[0] == "T_vs_C" and first[-1] == "changed" and first[2].startswith("F")
    assert first[3] == "none"                                  # 候補も注釈も無い合成データ
    assert out["record"]["written"] is True and "next" in out
    assert list(out)[:2] == ["status", "confounding"]          # 交絡は先頭（読み落とさせない）
    assert out["confounding"].keys() >= {"batch_table", "by_contrast", "warnings"}
    assert out["matrix_source"] == "matrix_id"
    assert out["p_adjust_method"] == "bh"
    assert c["confounded"] is False
    assert out["factorial"] is None and "f_top" not in out and "output" not in out


def test_two_groups_without_matrix_id_registers_implicit_matrix(ds):
    out = _cmp(reference_group="C", test_group="T")
    assert out["status"] == "success"
    assert out["matrix_source"] == "default"
    entry = stats_results.get_matrix(ds, out["matrix_id"])
    assert entry is not None and entry["record_id"]
    assert out["record"]["preprocess_record"]["written"] is True
    again = _cmp(reference_group="C", test_group="T")
    assert again["matrix_id"] == out["matrix_id"]
    assert "preprocess_record" not in again["record"]        # 既存の行列は再利用


def test_three_groups_reference_first_and_omnibus(tmp_path):
    state = _install(_ds(groups=("C", "T", "U")), tmp_path)
    mx = _pre()["matrix_id"]
    out = _cmp(groups=["C", "T", "U"], matrix_id=mx)
    assert [c["contrast_id"] for c in out["contrasts"]] == ["T_vs_C", "U_vs_C"]
    assert all(c["reference"] == "C" for c in out["contrasts"])
    assert out["omnibus"]["n_tested"] == 10
    assert out["omnibus"]["n_q_below"] >= 3
    assert out["group_n"] == {"C": 4, "T": 4, "U": 4}
    assert out["bh_population"]["n_by_contrast"] == {"T_vs_C": 10, "U_vs_C": 10}
    heads = {line.split("\t")[0] for line in out["top"].splitlines()[1:]}
    assert heads == {"T_vs_C", "U_vs_C"}
    assert state.last_differential["result_id"] == out["result_id"]


def test_two_by_two_comparisons_with_interaction(tmp_path):
    _install(_ds(per_group=6, sex=True), tmp_path)
    mx = _pre()["matrix_id"]
    comparisons = [
        {"id": "T_vs_C_F", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "F"}},
        {"id": "T_vs_C_M", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "M"}},
        {"id": "inter", "kind": "interaction", "factors": ["group", "sex"],
         "levels": {"group": ["C", "T"], "sex": ["F", "M"]}},
    ]
    out = _cmp(comparisons=comparisons, factors=["group", "sex"], matrix_id=mx)
    assert out["status"] == "success"
    assert [c["contrast_id"] for c in out["contrasts"]] == ["T_vs_C_F", "T_vs_C_M", "inter"]
    inter = _contrast(out, "inter")
    assert inter["presence"]["available"] is False
    assert _contrast(out, "T_vs_C_F")["verdicts"]["changed"] >= 3
    assert set(out["group_n"]) == {"C:F", "C:M", "T:F", "T:M"}


def test_exactly_one_group_specification_required(ds):
    mx = _pre()["matrix_id"]
    comp = [{"id": "x", "factor": "group", "reference": "C", "test": "T"}]
    bad = [
        {},                                                                     # 何も無い
        {"reference_group": "C"},                                               # 片方だけ
        {"test_group": "T"},
        {"reference_group": "C", "test_group": "T", "groups": ["C", "T"]},      # 2 つ
        {"groups": ["C", "T"], "comparisons": comp, "factors": ["group"]},
        {"groups": ["C"]},                                                      # 2 群未満
        {"comparisons": comp},                                                  # factors が無い
        {"reference_group": "C", "test_group": "T", "factors": ["group"]},      # factors は comparisons だけ
    ]
    for kw in bad:
        out = _cmp(matrix_id=mx, **kw)
        assert out["error"]["code"] == "DATASET_BAD_REQUEST", kw
    assert _cmp(matrix_id=mx, comparisons=comp, factors=["group"])["status"] == "success"


def test_unknown_group_is_an_error_envelope(ds):
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="NOPE", matrix_id=mx)
    assert "error" in out and "needs_input" not in out


def test_comparison_ids_are_sanitized(tmp_path):
    state = _ds()
    for row in state.sample_metadata_rows:
        if row["group"] == "T":
            row["group"] = "high dose"
    _install(state, tmp_path)
    out = _cmp(reference_group="C", test_group="high dose", matrix_id=_pre()["matrix_id"])
    assert out["status"] == "success"
    assert out["contrasts"][0]["contrast_id"] == "high_dose_vs_C"
    assert out["contrasts"][0]["test"] == "high dose"


# --- 共変量・ブロック ---------------------------------------------------------------------

def _with_columns(tmp_path):
    state = _ds(sex=True)
    for i, row in enumerate(state.sample_metadata_rows):
        if row["role"] == "sample":
            k = int(row["subject_id"][1:])
            row["age"] = str(20 + k * 3 + (0 if row["group"] == "C" else 1))
    return _install(state, tmp_path)


def test_covariate_type_inferred_numeric_vs_categorical(tmp_path):
    state = _with_columns(tmp_path)
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", covariates=["age", "sex"], matrix_id=mx)
    assert out["status"] == "success"
    spec = state.last_differential["spec"]
    assert spec["design"]["covariates"] == [{"name": "age", "type": "numeric"},
                                            {"name": "sex", "type": "categorical"}]
    missing = _cmp(reference_group="C", test_group="T", covariates=["nope"], matrix_id=mx)
    assert missing["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "age" in missing["error"]["details"]["available"]


def test_block_only_subject_id(tmp_path):
    state = _ds()
    for row in state.sample_metadata_rows:                 # 同じ個体が両方の群に入る（対応のある計画）
        if row["role"] == "sample":
            row["subject_id"] = "S" + row["subject_id"][1:]
    ds = _install(state, tmp_path)
    mx = _pre(collapse_replicates=False)["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", block="subject_id", matrix_id=mx)
    assert out["status"] == "success"
    assert ds.last_differential["spec"]["design"]["block"] == "subject_id"
    bad = _cmp(reference_group="C", test_group="T", block="batch", matrix_id=mx)
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


# --- 行列 ID ----------------------------------------------------------------------------

def test_matrix_id_from_preprocess_is_reused(ds):
    pre = _pre(normalize="tic")
    out = _cmp(reference_group="C", test_group="T", matrix_id=pre["matrix_id"])
    assert out["matrix_id"] == pre["matrix_id"]
    assert set(ds.stat_matrices) == {pre["matrix_id"]}               # 暗黙の前処理を足さない
    entry = stats_results.get_matrix(ds, pre["matrix_id"])
    spec = ds.last_differential["spec"]
    assert spec["preprocess"] == entry["config"]
    assert spec["preprocess"]["normalize"]["method"] == "tic"


def test_matrix_id_carries_its_groups_and_rejects_group_arguments(ds):
    mx = _pre(group_levels=["r0", "r1"])["matrix_id"]
    out = _cmp(reference_group="r0", test_group="r1", matrix_id=mx)
    assert out["status"] == "success"
    assert ds.last_differential["spec"]["input"]["group_overrides"]
    clash = _cmp(reference_group="C", test_group="T", matrix_id=mx, group_levels=["C", "T"])
    assert clash["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "dataset_preprocess" in clash["error"]["message"]


def test_unknown_or_stale_matrix_id_is_missing_state(ds):
    out = _cmp(reference_group="C", test_group="T", matrix_id="mx_nope")
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["state"] == "stat_matrix"
    assert out["error"]["required_tools"] == ["dataset_preprocess"]
    mx = _pre()["matrix_id"]
    apply_exclusions(ds, samples=["C0_r0"], features=[], mode="add", reason="bad")
    stale = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    assert stale["error"]["code"] == "missing_state" and stale["error"]["state"] == "stat_matrix"


def test_compare_without_dataset_is_missing_state():
    session_state.session.dataset = None
    out = _cmp(reference_group="C", test_group="T")
    assert out["error"]["code"] == "missing_state" and out["error"]["state"] == "dataset"


# --- needs_input ------------------------------------------------------------------------

def test_compare_needs_input_then_retry(ds):
    ds.quant_state = {"normalized": True, "method": "is_ratio", "unit": "ratio", "record_id": "q1"}
    out = _cmp(reference_group="C", test_group="T")
    body = out["needs_input"]
    assert "error" not in out
    assert body["field"] == "preprocess.drift.mode" and body["argument"] == "drift_mode"
    assert "dataset_preprocess" in body["reason"] and "matrix_id" in body["reason"]
    assert body["retry_tools"] == ["dataset_preprocess"]
    assert not ds.stat_matrices
    mx = _pre(drift_mode="off")["matrix_id"]
    ok = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    assert ok["status"] == "success" and ok["matrix_id"] == mx


# --- 戻り値 ----------------------------------------------------------------------------

def test_compare_payload_strict_json(ds):
    mx = _pre()["matrix_id"]
    text = dataset_compare(reference_group="C", test_group="T", matrix_id=mx)
    out = _strict(text)
    assert out["prior"]["df_prior"] == "inf"                          # 分散が揃うと無限大
    assert len(text) < 6000
    # 記録も厳密な JSON
    path = Path(ds.analysis_dir) / "reports" / "stats" / f"{out['record']['record_id']}.json"
    _strict(path.read_text(encoding="utf-8"))


def test_long_unit_lists_truncated_in_payload(tmp_path, monkeypatch):
    _install(_ds(per_group=4), tmp_path)
    mx = _pre()["matrix_id"]
    long_ids = [f"u{i}" for i in range(60)]

    real = stats_tools.spec_v2.run_spec_v2

    def noisy(*a, **k):
        res = real(*a, **k)
        res["analyses"][0]["result"]["warnings"].append(
            {"code": "units_missing_model_values", "message": "m", "unit_ids": long_ids})
        return res
    monkeypatch.setattr(stats_tools.spec_v2, "run_spec_v2", noisy)
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    w = next(w for w in out["warnings"] if w["code"] == "units_missing_model_values")
    assert len(w["unit_ids"]) == 20 and w["unit_ids_n_more"] == 40


def test_result_registered_and_last_differential_points_to_it(ds):
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    rec = ds.results[out["result_id"]]
    assert rec["kind"] == "compare" and ds.last_differential is rec
    assert rec["matrix_id"] == mx and rec["rule_hash"] == out["rule_hash"]
    # 同じ呼び出しは同じ ID、閾値違いは別の記録
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx)["result_id"] == out["result_id"]
    other = _cmp(reference_group="C", test_group="T", matrix_id=mx, log2fc_threshold=3.0)
    assert other["result_id"] != out["result_id"] and other["rule_hash"] == out["rule_hash"]
    assert ds.last_differential["result_id"] == other["result_id"]


def test_spec_in_result_reproduces_numbers_via_run_spec_v2(tmp_path):
    state = _install(_ds(), tmp_path)
    apply_exclusions(state, samples=["C0_r0"], features=["f9"], mode="add", reason="bad")
    out = _cmp(reference_group="C", test_group="T")                  # 暗黙の前処理・除外つき
    rec = state.last_differential
    spec = normalize_spec_v2(json.loads(json.dumps(rec["spec"])))
    assert spec["input"]["exclude_samples"] == ["C0_r0"]
    assert spec["input"]["exclude_features"] == ["f9"]
    rerun = run_spec_v2(spec, state)
    assert rerun["matrix_id"] == out["matrix_id"] and rerun["rule_hash"] == out["rule_hash"]
    old = {(r["feature_id"], r["contrast_id"]): (r["verdict"], r["q"])
           for r in rec["engine"]["rows"]}
    new = {(r["feature_id"], r["contrast_id"]): (r["verdict"], r["q"])
           for r in rerun["analyses"][0]["result"]["rows"]}
    assert old == new


def test_log2fc_threshold_does_not_change_verdicts(ds):
    mx = _pre()["matrix_id"]
    low = _cmp(reference_group="C", test_group="T", matrix_id=mx, log2fc_threshold=0.0)
    high = _cmp(reference_group="C", test_group="T", matrix_id=mx, log2fc_threshold=100.0)
    assert _contrast(low, "T_vs_C")["verdicts"] == _contrast(high, "T_vs_C")["verdicts"]
    assert _contrast(low, "T_vs_C")["large_changes"] >= 3
    assert _contrast(high, "T_vs_C")["large_changes"] == 0
    assert low["rule_hash"] == high["rule_hash"]
    bad = _cmp(reference_group="C", test_group="T", matrix_id=mx, log2fc_threshold=-1)
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


def test_top_n_limits_rows_per_contrast(ds):
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx, top_n=2)
    assert len(out["top"].splitlines()) == 3
    none = _cmp(reference_group="C", test_group="T", matrix_id=mx, top_n=0)
    assert none["top"] == ""
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx, top_n=-1)["error"]["code"] \
        == "DATASET_BAD_REQUEST"


def _confounded(tmp_path):
    state = _ds()
    for row in state.sample_metadata_rows:
        if row["role"] == "sample":
            row["batch"] = "B1" if row["group"] == "C" else "B2"
    return _install(state, tmp_path)


def test_confounded_compared_with_marks_by_default_and_refused_when_disallowed(tmp_path):
    _confounded(tmp_path)
    mx = _pre()["matrix_id"]
    ok = _cmp(reference_group="C", test_group="T", matrix_id=mx)      # 既定 allow_confounded=True
    assert ok["status"] == "success"
    assert list(ok)[:2] == ["status", "confounding"]
    assert ok["confounding"]["batch_status"] == "complete"
    assert ok["confounding"]["by_contrast"] == {"T_vs_C": {"batch_status": "complete",
                                                           "confounded": True}}
    assert any(w["code"] == "batch_group_complete_confounding" for w in ok["confounding"]["warnings"])
    assert _contrast(ok, "T_vs_C")["confounded"] is True
    assert _contrast(ok, "T_vs_C")["batch_status"] == "complete"
    assert any(w["code"] == "batch_group_complete_confounding" for w in ok["warnings"])
    refused = _cmp(reference_group="C", test_group="T", matrix_id=mx, allow_confounded=False)
    assert refused["error"]["code"] == "STATISTIC_CONFOUNDED"
    assert refused["error"]["details"]["contrasts"] == ["T_vs_C"]


def test_compare_writes_record_with_parent_matrix_record(ds):
    pre = _pre()
    out = _cmp(reference_group="C", test_group="T", matrix_id=pre["matrix_id"])
    rec_info = out["record"]
    assert rec_info["written"] is True
    path = Path(ds.analysis_dir) / "reports" / "stats" / f"{rec_info['record_id']}.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["step"] == "compare"
    assert data["parent_ids"] == [pre["record"]["record_id"]]
    assert data["inputs"]["matrix_id"] == pre["matrix_id"]
    assert data["inputs"]["result_id"] == out["result_id"]
    assert data["settings"]["spec"]["schema"] == "analysis-spec.v2"
    assert data["settings"]["log2fc_threshold"] == 1.0
    assert "top" not in data["summary"] and data["summary"]["result_id"] == out["result_id"]


def test_compare_implicit_preprocess_chain(ds):
    out = _cmp(reference_group="C", test_group="T")
    stats_dir = Path(ds.analysis_dir) / "reports" / "stats"
    cmp_rec = json.loads((stats_dir / f"{out['record']['record_id']}.json").read_text(encoding="utf-8"))
    pre_id = out["record"]["preprocess_record"]["record_id"]
    assert cmp_rec["parent_ids"] == [pre_id]
    assert json.loads((stats_dir / f"{pre_id}.json").read_text(encoding="utf-8"))["step"] == "preprocess"


def test_compare_without_report_dir_still_succeeds():
    state = _ds()
    session_state.session.dataset = state
    mx = _pre()["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    assert out["status"] == "success"
    assert out["record"]["written"] is False and "NO_REPORT_DIR" in out["record"]["reason"]


def test_handle_isolation(tmp_path):
    a, b = _ds(), _ds()
    a.analysis_dir = str(tmp_path / "a")
    b.analysis_dir = str(tmp_path / "b")
    session_state.session.dataset = a
    handle_a = session_state.session.current_handle_id
    session_state.session.dataset = b
    out = _cmp(reference_group="C", test_group="T", dataset_handle=handle_a)
    assert out["status"] == "success"
    assert a.last_differential is not None and b.last_differential is None
    assert not b.stat_matrices



# --- 修正ラウンド 1 ----------------------------------------------------------------------

def test_failed_implicit_call_keeps_matrix_and_reports_its_id(tmp_path):
    state = _confounded(tmp_path)
    out = _cmp(reference_group="C", test_group="T", allow_confounded=False)
    assert out["error"]["code"] == "STATISTIC_CONFOUNDED"
    mx = out["error"]["details"]["matrix_id"]
    assert mx in state.stat_matrices
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx)["status"] == "success"


def test_n_q_below_counts_q_equal_to_threshold():
    engine = {"verdict_rule": {"q_threshold": 0.05},
              "presence": {"detection_available": True, "bh_n_by_contrast": {"c": 3}},
              "rows": [{"contrast_id": "c", "reference": "a", "test": "b", "verdict": "changed",
                        "status": "tested", "log2fc": 1.0, "presence_q": q}
                       for q in (0.05, 0.0501, 0.01)],
              "omnibus": [{"f_p": 0.1, "f_q": 0.05}, {"f_p": 0.1, "f_q": 0.06}]}
    spec = {"comparisons": [{"id": "c"}]}
    assert stats_tools._contrast_summaries(engine, spec, 1.0)[0]["presence"]["n_q_below"] == 2
    assert stats_tools._omnibus_summary(engine)["n_q_below"] == 1


def test_top_n_capped(ds):
    mx = _pre()["matrix_id"]
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx, top_n=50)["status"] == "success"
    big = _cmp(reference_group="C", test_group="T", matrix_id=mx, top_n=51)
    assert big["error"]["code"] == "DATASET_BAD_REQUEST" and "result_id" in big["error"]["message"]


def test_reported_population_needs_its_ids_and_is_written_into_spec(ds):
    mx = _pre()["matrix_id"]
    no_ids = _cmp(reference_group="C", test_group="T", matrix_id=mx, bh_population="reported")
    assert no_ids["error"]["code"] == "DATASET_BAD_REQUEST"
    # reported 以外で ID を渡しても誤りにしない（効かない引数の組は notes に書き、母集団は行列のまま）
    stray = _cmp(reference_group="C", test_group="T", matrix_id=mx, reported_feature_ids=["f0"])
    assert stray["status"] == "success"
    assert stray["bh_population"]["name"] == "matrix"
    assert stray["bh_population"]["n_by_contrast"] == {"T_vs_C": 10}
    assert len(stray["notes"]) == 1 and 'bh_population="reported"' in stray["notes"][0]
    assert "reported_feature_ids" not in ds.last_differential["spec"]["analyses"][0]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx, bh_population="reported",
               reported_feature_ids=["f0", "f1", "f2", "f1"])
    assert out["status"] == "success" and "notes" not in out
    assert out["bh_population"]["name"] == "reported"
    assert out["bh_population"]["n_by_contrast"] == {"T_vs_C": 3}
    analysis = ds.last_differential["spec"]["analyses"][0]
    assert analysis["bh_population"] == "reported"
    assert analysis["reported_feature_ids"] == ["f0", "f1", "f2"]
    path = Path(ds.analysis_dir) / "reports" / "stats" / f"{out['record']['record_id']}.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["settings"]["spec"]["analyses"][0]["reported_feature_ids"] == ["f0", "f1", "f2"]
    # q は報告した 3 特徴だけで補正する（行列の母集団とは別の数字）
    rows = {r["feature_id"]: r for r in ds.last_differential["engine"]["rows"]}
    assert rows["f5"]["q"] is None and rows["f0"]["q"] is not None


def test_p_adjust_method_reaches_payload_record_and_q(ds):
    mx = _pre()["matrix_id"]
    bh = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    bh_rows = {r["feature_id"]: r for r in ds.last_differential["engine"]["rows"]}
    holm = _cmp(reference_group="C", test_group="T", matrix_id=mx, p_adjust="holm")
    rec = ds.last_differential
    assert holm["p_adjust_method"] == "holm" and rec["p_adjust_method"] == "holm"
    assert holm["bh_population"]["p_adjust_method"] == "holm"
    assert holm["rule_hash"] != bh["rule_hash"]
    assert rec["spec"]["analyses"][0]["p_adjust"] == "holm"
    for r in rec["engine"]["rows"]:                            # FWER は FDR 以上
        if r["q"] is not None:
            assert r["q"] >= bh_rows[r["feature_id"]]["q"] - 1e-12
    none = _cmp(reference_group="C", test_group="T", matrix_id=mx, p_adjust="none")
    assert all(r["q"] == r["p"] for r in ds.last_differential["engine"]["rows"])
    assert none["p_adjust_method"] == "none"
    assert ds.results[bh["result_id"]]["p_adjust_method"] == "bh"
    bad = _cmp(reference_group="C", test_group="T", matrix_id=mx, p_adjust="fdr")
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


# --- matrix_id の省略（直近の dataset_preprocess の行列） ------------------------------------

def test_omitted_matrix_id_uses_latest_preprocess_of_the_handle(ds):
    first = _pre(normalize="tic")["matrix_id"]
    second = _pre(normalize="median")["matrix_id"]
    assert first != second
    out = _cmp(reference_group="C", test_group="T")
    assert out["matrix_id"] == second and out["matrix_source"] == "latest_preprocess"
    assert "preprocess_record" not in out["record"]               # 既定の前処理を足さない
    assert set(ds.stat_matrices) == {first, second}
    assert ds.last_differential["spec"]["preprocess"]["normalize"]["method"] == "median"
    named = _cmp(reference_group="C", test_group="T", matrix_id=first)
    assert named["matrix_id"] == first and named["matrix_source"] == "matrix_id"
    # 前処理を作り直すと、それが直近になる（同じ ID でも順番は新しくなる）
    _pre(normalize="tic")
    assert _cmp(reference_group="C", test_group="T")["matrix_id"] == first
    # 除外で古くなれば既定の前処理へ
    apply_exclusions(ds, samples=["C0_r0"], features=[], mode="add", reason="bad")
    stale = _cmp(reference_group="C", test_group="T")
    assert stale["matrix_source"] == "default" and stale["matrix_id"] not in (first, second)


def test_omitted_matrix_id_with_group_arguments_matches_groups(ds):
    mx = _pre(group_levels=["r0", "r1"])["matrix_id"]
    same = _cmp(reference_group="r0", test_group="r1", group_levels=["r0", "r1"])
    assert same["matrix_id"] == mx and same["matrix_source"] == "latest_preprocess"
    plain = _cmp(reference_group="r0", test_group="r1")             # 群の指定なし → 直近の行列の群
    assert plain["matrix_id"] == mx and plain["matrix_source"] == "latest_preprocess"
    other = _cmp(reference_group="C", test_group="T", group_levels=["C", "T"])
    assert other["status"] == "success"
    assert other["matrix_source"] == "default" and other["matrix_id"] != mx


def test_default_matrix_is_not_taken_as_latest_preprocess(ds):
    implicit = _cmp(reference_group="C", test_group="T")
    assert implicit["matrix_source"] == "default"
    again = _cmp(reference_group="C", test_group="T")
    assert again["matrix_source"] == "default" and again["matrix_id"] == implicit["matrix_id"]
    assert stats_results.latest_preprocess_matrix(ds) is None


def test_latest_preprocess_is_per_handle(tmp_path):
    a, b = _ds(), _ds()
    a.analysis_dir = str(tmp_path / "a")
    b.analysis_dir = str(tmp_path / "b")
    session_state.session.dataset = a
    handle_a = session_state.session.current_handle_id
    mx_a = _pre(normalize="tic")["matrix_id"]
    session_state.session.dataset = b
    out = _cmp(reference_group="C", test_group="T")
    assert out["matrix_source"] == "default" and out["matrix_id"] != mx_a
    on_a = _cmp(reference_group="C", test_group="T", dataset_handle=handle_a)
    assert on_a["matrix_id"] == mx_a and on_a["matrix_source"] == "latest_preprocess"


# --- 多群・要因の F ----------------------------------------------------------------------

def test_factorial_effects_returned_beside_contrasts(tmp_path):
    _install(_ds(per_group=6, sex=True), tmp_path)
    mx = _pre()["matrix_id"]
    comparisons = [{"id": "T_vs_C_F", "factor": "group", "reference": "C", "test": "T",
                    "at": {"sex": "F"}}]
    out = _cmp(comparisons=comparisons, factors=["group", "sex"], matrix_id=mx)
    assert out.get("status") == "success", out
    assert set(out["factorial"]) == {"group", "sex", "group:sex"}
    assert out["factorial"]["group"]["df1"] == 1 and out["factorial"]["group"]["n_tested"] == 10
    assert out["factorial"]["group"]["n_q_below"] >= 3
    lines = out["f_top"].splitlines()
    assert lines[0].split("\t") == ["effect", "feature_id", "name", "id_level", "f", "df1", "df2", "f_q"]
    assert {line.split("\t")[0] for line in lines[1:]} >= {"group"}
    engine = session_state.session.dataset.last_differential["engine"]
    assert engine["factorial"]["effects"]["group:sex"]["factors"] == ["group", "sex"]
    assert engine["omnibus"] is not None                       # cell 全体の F は別の欄のまま


def test_three_groups_f_top_lists_omnibus(tmp_path):
    _install(_ds(groups=("C", "T", "U")), tmp_path)
    out = _cmp(groups=["C", "T", "U"], matrix_id=_pre()["matrix_id"])
    lines = out["f_top"].splitlines()
    assert {line.split("\t")[0] for line in lines[1:]} == {"omnibus"}
    assert out["factorial"] is None


# --- output_path（全行の TSV） --------------------------------------------------------------

def _read_tsv(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def test_output_path_writes_all_rows_and_columns(ds, tmp_path):
    from metabolomix.analysis import compare_table

    mx = _pre()["matrix_id"]
    target = tmp_path / "out" / "compare.tsv"
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx, output_path=str(target), top_n=1)
    assert out["output"] == {"written": True, "path": str(target), "n_rows": 10,
                             "ftests_path": None, "n_ftest_rows": 0}
    head = target.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert head == list(compare_table.COLUMNS)
    rows = _read_tsv(target)
    assert len(rows) == 10 and {r["feature_id"] for r in rows} == {f"f{k}" for k in range(10)}
    engine = {r["feature_id"]: r for r in ds.last_differential["engine"]["rows"]}
    for r in rows:
        e = engine[r["feature_id"]]
        assert r["contrast_id"] == "T_vs_C" and r["name"] == f"F{r['feature_id'][1:]}"
        assert r["id_level"] == "none" and r["verdict"] == e["verdict"]
        assert float(r["log2fc"]) == pytest.approx(e["log2fc"], rel=1e-9)
        assert float(r["se"]) == pytest.approx(e["se"], rel=1e-9)
        assert float(r["df"]) == pytest.approx(e["df"], rel=1e-9)
        assert int(r["n_reference"]) == e["n_reference"] == 4
        assert r["detected_test"] == str(e["detected_test"])
        assert r["confounded"] == "false"
        assert r["constant_in_group"] == ""
    # 戻り値の top は切っても、ファイルは全行
    assert len(out["top"].splitlines()) == 2


def test_output_path_writes_ftests_for_multi_group(tmp_path):
    _install(_ds(groups=("C", "T", "U")), tmp_path)
    target = tmp_path / "multi.tsv"
    out = _cmp(groups=["C", "T", "U"], matrix_id=_pre()["matrix_id"], output_path=str(target))
    assert out["output"]["n_rows"] == 20
    ftests = tmp_path / "multi.ftests.tsv"
    assert out["output"]["ftests_path"] == str(ftests) and out["output"]["n_ftest_rows"] == 10
    rows = _read_tsv(ftests)
    assert {r["kind"] for r in rows} == {"omnibus"} and {r["effect"] for r in rows} == {""}
    assert all(r["df1"] == "2" for r in rows)


def test_output_path_validation_and_write_failure(ds, tmp_path, monkeypatch):
    from metabolomix.analysis import compare_table

    mx = _pre()["matrix_id"]
    folder = _cmp(reference_group="C", test_group="T", matrix_id=mx, output_path=str(tmp_path))
    assert folder["error"]["code"] == "DATASET_BAD_REQUEST"
    empty = _cmp(reference_group="C", test_group="T", matrix_id=mx, output_path=" ")
    assert empty["error"]["code"] == "DATASET_BAD_REQUEST"

    def boom(*a, **k):
        raise OSError("disk full")
    monkeypatch.setattr(compare_table, "atomic_write_bytes", boom)
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx, output_path=str(tmp_path / "x.tsv"))
    assert out["status"] == "success"                         # 比較は成功し、書けなかった理由を返す
    assert out["output"]["written"] is False and "disk full" in out["output"]["reason"]


def test_japanese_group_names_give_unique_stable_ids(tmp_path):
    state = _ds(groups=("C", "T", "U"))
    names = {"C": "対照", "T": "高用量", "U": "低用量"}
    for row in state.sample_metadata_rows:
        if row["role"] == "sample":
            row["group"] = names[row["group"]]
    _install(state, tmp_path)
    mx = _pre()["matrix_id"]
    out = _cmp(groups=["対照", "高用量", "低用量"], matrix_id=mx)
    ids = [c["contrast_id"] for c in out["contrasts"]]
    assert ids == ["g1_vs_g0", "g2_vs_g0"]
    assert [c["test"] for c in out["contrasts"]] == ["高用量", "低用量"]
    assert [c["reference"] for c in out["contrasts"]] == ["対照", "対照"]
    again = _cmp(groups=["対照", "高用量", "低用量"], matrix_id=mx)
    assert [c["contrast_id"] for c in again["contrasts"]] == ids
    two = _cmp(reference_group="対照", test_group="高用量", matrix_id=mx)
    assert two["contrasts"][0]["contrast_id"] == "g1_vs_g0"


def test_matrix_id_path_does_not_rerun_preprocess(ds, monkeypatch):
    mx = _pre()["matrix_id"]

    def boom(*a, **k):
        raise AssertionError("preprocess must not run")
    monkeypatch.setattr(spec_v2, "preprocess", boom)               # spec_v2 が呼ぶ場所
    monkeypatch.setattr(stats_tools.stats_preprocess, "preprocess", boom)
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx)["status"] == "success"


def test_matrix_id_spec_reproduces_numbers(tmp_path):
    state = _install(_ds(), tmp_path)
    apply_exclusions(state, samples=["C0_r0"], features=["f9"], mode="add", reason="bad")
    mx = _pre(group_levels=["C", "T"], normalize="tic")["matrix_id"]
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    rec = state.last_differential
    rerun = run_spec_v2(normalize_spec_v2(json.loads(json.dumps(rec["spec"]))), state)
    assert rerun["matrix_id"] == mx == out["matrix_id"] and rerun["rule_hash"] == rec["rule_hash"]
    old = {(r["feature_id"], r["contrast_id"]): (r["verdict"], r["q"]) for r in rec["engine"]["rows"]}
    new = {(r["feature_id"], r["contrast_id"]): (r["verdict"], r["q"])
           for r in rerun["analyses"][0]["result"]["rows"]}
    assert old == new


# --- レビューの指摘（交絡の 3 値・対比ごとの印・行列の出所・F の df2） -------------------------

def _split_batches(state):
    """C → B1、T → B2、U は個体で 2 つに割る: まとめた表は partial、T_vs_C だけが完全交絡。"""
    for row in state.sample_metadata_rows:
        if row["role"] != "sample":
            continue
        g, subject = row["group"], row["subject_id"]
        row["batch"] = ("B1" if g == "C" else "B2" if g == "T"
                        else "B1" if subject in ("U0", "U1") else "B2")
    return state


def test_one_confounded_pair_among_three_groups_is_marked_not_stopped(tmp_path):
    ds = _install(_split_batches(_ds(groups=("C", "T", "U"))), tmp_path)
    mx = _pre()["matrix_id"]
    out = _cmp(groups=["C", "T", "U"], matrix_id=mx)
    assert out["status"] == "success"
    assert out["confounding"]["batch_status"] == "partial"
    assert out["confounding"]["by_contrast"]["T_vs_C"] == {"batch_status": "complete", "confounded": True}
    assert out["confounding"]["by_contrast"]["U_vs_C"]["confounded"] is False
    assert _contrast(out, "T_vs_C")["confounded"] is True
    assert any(w["code"] == "contrast_complete_confounding" for w in out["warnings"])
    # 消費側: 対比ごとの完全交絡にも未調整の印が付き、部分交絡の対比には付かない
    rec = ds.last_differential
    tc = stats_results.as_two_group(rec, "T_vs_C")["provenance"]
    assert tc["comparison"]["unadjusted_confounded"] is True
    assert tc["comparison"]["confounding"]["batch_status"] == "complete"
    assert tc["comparison"]["confounding"]["pooled_batch_status"] == "partial"
    assert set(tc["comparison"]["confounding"]["batch_table"]) == {"C", "T"}
    assert stats_results.UNADJUSTED_CONFOUNDED_EXPLANATION in tc["warnings"]
    assert "comparison" not in stats_results.as_two_group(rec, "U_vs_C")["provenance"]
    # allow_confounded=false の明示だけが止める（どの対比が完全交絡かを details に）
    refused = _cmp(groups=["C", "T", "U"], matrix_id=mx, allow_confounded=False)
    assert refused["error"]["code"] == "STATISTIC_CONFOUNDED"
    assert refused["error"]["details"]["contrasts"] == ["T_vs_C"]


def test_default_call_on_complete_confounding_does_not_claim_an_explicit_choice(tmp_path):
    ds = _confounded(tmp_path)
    out = _cmp(reference_group="C", test_group="T")             # allow_confounded を渡さない
    assert out["status"] == "success"
    warnings = stats_results.as_two_group(ds.last_differential)["provenance"]["warnings"]
    text = next(w for w in warnings if "未調整" in w)
    assert "明示" not in text and "既定" in text


def test_unknown_batches_leave_confounded_unknown_not_false(tmp_path):
    state = _ds()
    for row in state.sample_metadata_rows:
        if row["role"] == "sample" and row["group"] == "T":
            row["batch"] = None
    ds = _install(state, tmp_path)
    mx = _pre()["matrix_id"]
    path = tmp_path / "rows.tsv"
    out = _cmp(reference_group="C", test_group="T", matrix_id=mx, output_path=str(path))
    assert out["confounding"]["batch_status"] == "not_evaluable"
    assert out["confounding"]["by_contrast"]["T_vs_C"] == {"batch_status": "not_evaluable",
                                                           "confounded": None}
    assert _contrast(out, "T_vs_C")["confounded"] is None
    assert all(r["confounded"] is None for r in ds.last_differential["engine"]["rows"])
    assert {row["confounded"] for row in _read_tsv(path)} == {""}       # 空欄（false ではない）
    assert "comparison" not in stats_results.as_two_group(ds.last_differential)["provenance"]
    # 評価できない交絡は allow_confounded=false でも止めない（完全交絡ではない）
    assert _cmp(reference_group="C", test_group="T", matrix_id=mx,
                allow_confounded=False)["status"] == "success"


def test_error_after_latest_preprocess_names_the_matrix_and_its_source(ds):
    mx = _pre(group_factors=[["C"], ["T"]])["matrix_id"]
    out = _cmp(reference_group="C", test_group="T")
    assert out["error"]["code"] == "DESIGN_INVALID"
    details = out["error"]["details"]
    assert details["matrix_id"] == mx and details["matrix_source"] == "latest_preprocess"
    assert "dataset_preprocess" in details["matrix_note"]
    # 名指しした行列のエラーには足さない（呼び手が知っている）
    named = _cmp(reference_group="C", test_group="T", matrix_id=mx)
    assert "matrix_source" not in (named["error"].get("details") or {})


def test_error_after_implicit_default_matrix_reports_source(tmp_path):
    state = _confounded(tmp_path)
    out = _cmp(reference_group="C", test_group="T", allow_confounded=False)
    assert out["error"]["details"]["matrix_source"] == "default"
    assert out["error"]["details"]["matrix_id"] in state.stat_matrices
    assert "matrix_note" not in out["error"]["details"]


def test_infinite_df2_is_written_as_inf_not_blank(tmp_path):
    import math
    ds = _install(_ds(groups=("C", "T", "U")), tmp_path)
    mx = _pre()["matrix_id"]
    path = tmp_path / "rows.tsv"
    text = dataset_compare(groups=["C", "T", "U"], matrix_id=mx, output_path=str(path), top_n=50)
    out = _strict(text)
    engine = ds.last_differential["engine"]
    assert math.isinf(engine["prior"]["df_prior"])                  # この合成データは事前分布の自由度が ∞
    assert all(o["df2"] == math.inf for o in engine["omnibus"] if o["f"] is not None)
    f_lines = out["f_top"].splitlines()
    assert len(f_lines) > 1 and all(line.split("\t")[6] == "inf" for line in f_lines[1:])
    ftests = _read_tsv(Path(out["output"]["ftests_path"]))
    assert ftests and {row["df2"] for row in ftests} == {"inf"}


def test_overview_compares_three_groups_when_only_one_pair_is_confounded(tmp_path):
    from metabolomix.tools.overview_tools import dataset_overview
    ds = _install(_split_batches(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _strict(dataset_overview())
    assert out["status"] == "success"
    assert any(f.startswith("G3 moderated F") for f in out["global"]["facts"]), out["global"]
    rec = ds.results[out["result_ids"]["compare"]]
    assert rec["p_adjust_method"] == "bh"                     # どの比較の記録にも補正の方法が残る
    assert rec["engine"]["confounding"]["by_contrast"]["T_vs_C"]["confounded"] is True
    assert all(r["confounded"] is True for r in rec["engine"]["rows"] if r["contrast_id"] == "T_vs_C")
