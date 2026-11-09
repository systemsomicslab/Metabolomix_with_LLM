"""dataset_feature_search の with_group_labels（spec 2026-10-22 §4 の 5）。"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import session_state
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools, mztab_tools
from tests.engine_chain import declare_groups, engine_compare, hc_ra
from tests.reanalysis_fixtures import write_mztab


@pytest.fixture(autouse=True)
def loaded(tmp_path):
    session_state.session = session_state.AnalysisSession()
    handle_tools.load_dataset(mztab_path=str(write_mztab(tmp_path / "AlignResult-1.mzTab")))
    declare_groups(session_state.session.dataset, hc_ra)
    engine_compare("HC", "RA")
    yield
    session_state.session = session_state.AnalysisSession()


def _stats(**kwargs):
    out = json.loads(T.dataset_feature_search(feature_ids=["1"], **kwargs))
    lines = out["stats"].splitlines()
    return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]


def test_default_has_no_label_columns():
    plain = ["feature_id", "result_id", "status", "log2fc", "p", "q", "mean_a", "mean_b",
             "detected_a", "detected_b"]
    header, _ = _stats(merge_split=False)
    assert header == plain
    # merge_split（既定 true）は末尾に merged_from だけを足す
    header, _ = _stats()
    assert header == plain + ["merged_from"]


def test_labels_add_comparison_groups_and_direction():
    header, rows = _stats(with_group_labels=True)
    assert header[-5:] == ["comparison", "group_a", "group_b", "log2fc_test_over_ref", "merged_from"]
    row = dict(zip(header, rows[0]))
    assert row["comparison"] == "RA_vs_HC"
    assert (row["group_a"], row["group_b"]) == ("HC", "RA")
    assert row["log2fc_test_over_ref"] == row["log2fc"]
    assert float(row["log2fc"]) > 0
