"""`dataset_preprocess`（統計エンジンの前処理ツール）。fixture はテスト自身が tmp に作る。"""
from __future__ import annotations

import json

import pytest

from metabolomix.analysis import stats_results
from metabolomix.analysis.dataset_exclusions import apply_exclusions
from metabolomix.core import session_state
from metabolomix.tools import stats_tools
from metabolomix.tools.stats_tools import dataset_preprocess
from tests.config_fixtures import use_config
from tests.test_stats_results import _ds

# 小さな fixture で足切りに落とされないよう、足切りは使い道のあるものだけ残す。
FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")


@pytest.fixture
def ds(tmp_path):
    state = _ds()
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


def _call(**kw) -> dict:
    return json.loads(dataset_preprocess(**kw))


def test_preprocess_returns_matrix_id_and_registers_matrix(ds):
    out = _call(**FLAT_OFF)
    assert out["status"] == "success"
    entry = stats_results.get_matrix(ds, out["matrix_id"])
    assert entry is not None
    assert entry["prepared"].matrix_id == out["matrix_id"]
    assert out["n_units"] == len(entry["prepared"].units) == 8      # 8 個体 × 2 注入をまとめた
    assert out["n_qc"] == 4
    assert out["n_features"] == out["n_features_before"] == 10
    assert out["subject_source"]
    assert {s["step"] for s in out["steps"]} >= {"drift", "log2"}
    assert "next" in out


def test_preprocess_same_args_same_matrix_id(ds):
    a = _call(**FLAT_OFF)
    b = _call(**FLAT_OFF)
    c = _call(**{**FLAT_OFF, "normalize": "tic"})
    assert a["matrix_id"] == b["matrix_id"] != c["matrix_id"]
    assert set(ds.stat_matrices) == {a["matrix_id"], c["matrix_id"]}


def test_preprocess_none_disables_filter(ds):
    on = _call(**{**FLAT_OFF, "max_qc_rsd": 0.1})
    off = _call(**FLAT_OFF)
    assert on["n_features"] == 3 and off["n_features"] == 10
    cfg = stats_results.get_matrix(ds, off["matrix_id"])["config"]
    assert cfg["filters"]["qc_rsd"] is None and cfg["filters"]["detection"] is None


def test_preprocess_bad_argument_is_config_envelope(ds):
    out = _call(**{**FLAT_OFF, "normalize": "zscore"})
    assert out["error"]["code"] == "PREPROCESS_CONFIG_INVALID"
    assert out["error"]["details"]["field"] == "normalize.method"


def _normalized(ds):
    ds.quant_state = {"normalized": True, "method": "is_ratio", "unit": "ratio",
                      "record_id": "q_test_1"}


def test_preprocess_needs_input_for_normalized_with_qc(ds):
    _normalized(ds)
    out = _call(**{**FLAT_OFF, "drift_mode": "auto"})
    body = out["needs_input"]
    assert "error" not in out
    assert body["field"] == "preprocess.drift.mode"
    assert body["argument"] == "drift_mode"
    assert body["retry_tools"] == ["dataset_preprocess"]
    assert {o["value"] for o in body["options"]} == {"on", "off"}
    assert not ds.stat_matrices                                   # 行列は作らない


def test_preprocess_drift_off_after_needs_input(ds):
    _normalized(ds)
    out = _call(**{**FLAT_OFF, "drift_mode": "off"})
    assert out["status"] == "success"
    assert out["quant"] == {"normalized": True, "method": "is_ratio", "unit": "ratio"}


def test_preprocess_writes_record_under_reports_stats(ds, tmp_path):
    _normalized(ds)
    out = _call(**FLAT_OFF)
    record = out["record"]
    assert record["written"] is True
    path = tmp_path / "analysis" / "reports" / "stats" / f"{record['record_id']}.json"
    assert path.is_file()
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["step"] == "preprocess"
    assert data["inputs"]["matrix_id"] == out["matrix_id"]
    assert data["parent_ids"] == ["q_test_1"]
    assert data["settings"]["config"]["normalize"]["method"] == "none"
    assert data["summary"]["matrix_id"] == out["matrix_id"]
    # 2 回目は別の記録（追記専用）。行列 ID は同じ。
    again = _call(**FLAT_OFF)
    assert again["record"]["record_id"] != record["record_id"]
    assert again["matrix_id"] == out["matrix_id"]


def test_preprocess_without_report_dir_still_succeeds():
    state = _ds()                                    # analysis_dir 無し・設定の reports_dir も無し
    session_state.session.dataset = state
    out = _call(**FLAT_OFF)
    assert out["status"] == "success"
    assert out["record"]["written"] is False
    assert "NO_REPORT_DIR" in out["record"]["reason"]
    assert stats_results.get_matrix(state, out["matrix_id"]) is not None


def test_preprocess_record_write_failure_still_succeeds(ds, monkeypatch):
    def boom(*a, **k):
        raise OSError("disk full")
    monkeypatch.setattr(stats_tools.stats_records, "write_stats_record", boom)
    out = _call(**FLAT_OFF)
    assert out["status"] == "success"
    assert out["record"]["written"] is False and "disk full" in out["record"]["reason"]


def test_preprocess_uses_configured_reports_dir(ds, monkeypatch, tmp_path):
    target = tmp_path / "elsewhere"
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": str(target)})
    out = _call(**FLAT_OFF)
    assert (target / "stats" / f"{out['record']['record_id']}.json").is_file()


def test_preprocess_applies_dataset_exclusions(ds):
    base = _call(**FLAT_OFF)
    apply_exclusions(ds, samples=["C0_r0"], features=["f9"], mode="add", reason="bad injection")
    out = _call(**FLAT_OFF)
    assert out["matrix_id"] != base["matrix_id"]
    assert out["n_features"] == 9
    entry = stats_results.get_matrix(ds, out["matrix_id"])
    assert entry["exclusions"] == {"samples": ["C0_r0"], "features": ["f9"]}
    assert "f9" not in entry["prepared"].feature_ids
    assert out["n_units"] == 8                       # C0 は 1 注入が残る（個体は残る）


def test_preprocess_group_levels_from_sample_names(ds):
    out = _call(**FLAT_OFF, group_levels=["r0", "r1"])
    entry = stats_results.get_matrix(ds, out["matrix_id"])
    overrides = entry["group_overrides"]
    assert overrides["C0_r0"] == "r0" and overrides["T3_r1"] == "r1"
    assert not any(name.startswith("QC") for name in overrides)   # QC は上書きしない
    # 上書き違いは別の行列
    assert out["matrix_id"] != _call(**FLAT_OFF)["matrix_id"]


def test_preprocess_group_levels_not_exclusive_is_bad_request(ds):
    out = _call(**FLAT_OFF, group_levels=["c0", "r0"])           # C0_r0 が 2 つの値に当たる
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"


def test_preprocess_does_not_set_the_latest_results(ds):
    _call(**FLAT_OFF)
    assert ds.last_pca is None and ds.last_differential is None


def test_preprocess_without_dataset_is_missing_state():
    session_state.session.dataset = None
    out = _call()
    assert out["error"]["code"] == "missing_state"
    assert out["error"]["state"] == "dataset"
    assert out["error"]["required_tools"] == ["load_dataset"]


def test_preprocess_unknown_handle_is_not_found():
    out = _call(dataset_handle="dh_nope")
    assert out["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"


def test_preprocess_payload_is_strict_json_and_small(ds):
    text = dataset_preprocess(**FLAT_OFF)
    json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))
    assert len(text) < 6000


def test_preprocess_invalid_injection_columns_is_config_envelope(ds):
    out = _call(**FLAT_OFF, injection_columns=["a", ""])
    assert out["error"]["code"] == "PREPROCESS_CONFIG_INVALID"


def test_stats_tool_module_public_surface():
    assert stats_tools.__all__ == ["dataset_preprocess", "dataset_compare", "dataset_pca"]


# --- 工程の要約（足切りごとの内訳・群ごとの個体数） -------------------------------------------

def _step(out, name):
    return next(s for s in out["steps"] if s["step"] == name)


def test_filter_step_reports_each_filter(ds):
    out = _call(**{**FLAT_OFF, "max_qc_rsd": 0.1})
    flt = _step(out, "filter")
    assert flt["dropped"] == 7
    assert flt["filters"]["qc_rsd"]["status"] == "applied"
    assert flt["filters"]["qc_rsd"]["dropped"] == 7
    assert flt["filters"]["blank_ratio"] == {"status": "skipped", "reason": "disabled"}
    assert "dropped_ids" not in flt["filters"]["qc_rsd"]


def test_filter_step_reports_detection_skipped_without_mask(ds):
    ds.detected_mask = None
    out = _call(**{**FLAT_OFF, "detection_min_rate": 0.5})
    det = _step(out, "filter")["filters"]["detection"]
    assert det == {"status": "skipped", "reason": "no_detected_mask"}


def test_select_step_reports_per_role_counts(ds):
    ds.sample_metadata_rows[0]["role"] = "standard"      # QC0 を役割が対象外の試料にする
    out = _call(**FLAT_OFF)
    sel = _step(out, "select")
    assert sel["n_sample"] == 16 and sel["n_qc"] == 3
    assert sel["excluded_by_role"] == {"standard": 1}


def test_payload_has_group_counts_and_single_group_warning(ds):
    out = _call(**FLAT_OFF)
    assert out["groups"] == {"C": 4, "T": 4}
    assert not any(w.get("code") == "single_group" for w in out["warnings"])
    typo = _call(**FLAT_OFF, group_levels=["nomatch"])
    assert typo["groups"] == {"other": 8}
    assert any(w.get("code") == "single_group" and w["group"] == "other" for w in typo["warnings"])


# --- 詳細の切り詰め ---------------------------------------------------------------------

def test_error_details_are_truncated_at_every_depth(ds, monkeypatch):
    from metabolomix.core.atomic_io import DomainError

    def boom(*a, **k):
        raise DomainError("STUDY_DATA_INVALID", "bad",
                          {"unknown_samples": [f"s{i}" for i in range(50)],
                           "excluded_by_role": {"std": [f"x{i}" for i in range(50)], "few": ["a"]}})
    monkeypatch.setattr(stats_tools, "from_dataset", boom)
    err = _call(**FLAT_OFF)["error"]["details"]
    assert len(err["unknown_samples"]) == 20 and err["unknown_samples_n_more"] == 30
    assert len(err["excluded_by_role"]["std"]) == 20 and err["excluded_by_role"]["std_n_more"] == 30
    assert err["excluded_by_role"]["few"] == ["a"]


# --- 公開面 ----------------------------------------------------------------------------

def test_handle_isolation(tmp_path):
    a, b = _ds(), _ds()
    a.analysis_dir = str(tmp_path / "a")
    b.analysis_dir = str(tmp_path / "b")
    session_state.session.dataset = a
    handle_a = session_state.session.current_handle_id
    session_state.session.dataset = b                  # 直近は B
    out = _call(**FLAT_OFF, dataset_handle=handle_a)
    assert out["status"] == "success"
    assert out["matrix_id"] in a.stat_matrices
    assert not b.stat_matrices and session_state.session.dataset is b
    assert (tmp_path / "a" / "reports" / "stats" / f"{out['record']['record_id']}.json").is_file()
    assert not (tmp_path / "b").exists()


def test_invalid_config_file_does_not_fail_tool(ds, monkeypatch, tmp_path):
    cfg = tmp_path / "bad.toml"
    cfg.write_text('[msdial]\nlbm = "x.lbm2"\n', encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))
    out = _call(**FLAT_OFF)
    assert out["status"] == "success"
    assert out["record"]["written"] is False and "CONFIG_INVALID" in out["record"]["reason"]


def test_compact_truncates_long_strings_with_an_ellipsis():
    """要約の長い文字列（注意・理由）は黙って落とさず、上限の長さで「…」を付けて切る。"""
    from metabolomix.tools._matrix_store import _compact
    long_text = "あ" * 200
    out = _compact({"reason": long_text, "short": "ok", "edge": "x" * 80, "n": 3})
    assert out["short"] == "ok" and out["n"] == 3
    assert out["edge"] == "x" * 80                                # 上限ちょうどは切らない
    assert len(out["reason"]) == 80 and out["reason"].endswith("…")
    assert out["reason"][:-1] == long_text[:79]
    nested = _compact({"filters": {"note": "y" * 81}})
    assert nested["filters"]["note"] == "y" * 79 + "…"


def test_default_arguments_are_the_engine_defaults():
    import inspect

    from metabolomix.analysis.stats.preprocess import PREPROCESS_DEFAULTS, normalize_preprocess_config
    params = inspect.signature(dataset_preprocess).parameters
    flat = {name: p.default for name, p in params.items()
            if name in inspect.signature(stats_tools._config).parameters}
    assert normalize_preprocess_config(stats_tools._config(**flat)) == normalize_preprocess_config(None)
    assert normalize_preprocess_config(None) == PREPROCESS_DEFAULTS
