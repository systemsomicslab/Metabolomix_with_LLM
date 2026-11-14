"""`dataset_status(view="results")`（棚卸し B6-G05）と、要約の試料の役割が試料の表に従うこと（A10-T23）。"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import session_state
from metabolomix.tools import mztab_tools
from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess
from tests.test_stats_results import _ds

FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")


@pytest.fixture
def ds(tmp_path, monkeypatch):
    monkeypatch.setattr(session_state, "session", session_state.AnalysisSession())
    state = _ds()
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


def _table(text: str) -> list[dict]:
    lines = text.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _results() -> dict:
    return json.loads(mztab_tools.dataset_status(view="results"))


def test_results_view_is_empty_before_any_analysis(ds):
    out = _results()
    assert out["view"] == "results"
    assert out["n_matrices"] == 0 and out["matrices"] == ""
    assert out["n_results"] == 0 and out["results"] == ""


def test_results_view_lists_matrices_and_results(ds):
    mx = json.loads(dataset_preprocess(**FLAT_OFF))["matrix_id"]
    cmp = json.loads(dataset_compare(reference_group="C", test_group="T", matrix_id=mx))
    out = _results()
    assert out["n_matrices"] == 1 and out["n_results"] == 1
    matrix = _table(out["matrices"])[0]
    assert list(matrix) == ["matrix_id", "current", "drift", "filters", "normalize", "n_features",
                            "n_units", "excluded_samples", "excluded_features", "group_overrides"]
    assert matrix["matrix_id"] == mx and matrix["current"] == "True"
    assert int(matrix["n_features"]) == 10 and int(matrix["n_units"]) > 0
    result = _table(out["results"])[0]
    assert result["result_id"] == cmp["result_id"]
    assert result["kind"] == "compare" and result["matrix_id"] == mx
    assert result["contrasts"] == cmp["contrasts"][0]["contrast_id"]
    assert result["current"] == "True" and result["last"] == "True"


def test_results_view_marks_stale_results_after_metadata_change(ds):
    mx = json.loads(dataset_preprocess(**FLAT_OFF))["matrix_id"]
    json.loads(dataset_compare(reference_group="C", test_group="T", matrix_id=mx))
    ds.metadata_revision += 1      # 試料の表を直した（行列と結果は古くなる）
    out = _results()
    assert _table(out["matrices"])[0]["current"] == "False"
    assert _table(out["results"])[0]["current"] == "False"


def test_results_view_ignores_query_with_a_note(ds):
    """query は samples / tags でだけ効く: results では断らず、一覧をそのまま返して notes に書く。"""
    mx = json.loads(dataset_preprocess(**FLAT_OFF))["matrix_id"]
    json.loads(dataset_compare(reference_group="C", test_group="T", matrix_id=mx))
    plain = _results()
    out = json.loads(mztab_tools.dataset_status(view="results", query="x"))
    assert "error" not in out and "notes" not in plain
    assert (out["matrices"], out["results"]) == (plain["matrices"], plain["results"])
    assert len(out["notes"]) == 1 and "query" in out["notes"][0] and "view=results" in out["notes"][0]


def test_summary_sample_roles_follow_the_sample_table(ds):
    """要約の `samples`（name / role）は試料の表の役割（`_metadata_rows`）を読む（名前の推定ではない）。"""
    rows = _table(json.loads(mztab_tools.dataset_status())["samples"])
    by_name = {row["name"]: row["role"] for row in rows}
    assert by_name["QC0"] == "qc"
    assert by_name["C0_r0"] == "sample"
    # 名前からは QC と読めない試料でも、試料の表が qc と言えば qc。
    target = ds.sample_names.index("C0_r0")
    ds.sample_metadata_rows[target]["role"] = "qc"
    rows = _table(json.loads(mztab_tools.dataset_status())["samples"])
    assert {row["name"]: row["role"] for row in rows}["C0_r0"] == "qc"
