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
    assert body["retry_tools"] == ["dataset_preprocess", "dataset_compare", "dataset_pca"]
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


def test_preprocess_does_not_write_v1_pp_slot(ds):
    _call(**FLAT_OFF)
    assert ds.pp_matrix is None
    assert not getattr(ds, "pp_sample_names", None)
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
    assert stats_tools.__all__ == ["dataset_preprocess"]
