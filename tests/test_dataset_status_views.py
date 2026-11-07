"""`dataset_status` の views（B-2b Task 3）。

`summary` は群・役割・タグの件数を持ち、`samples` は試料の TSV を因子トークンで絞り、
`tags` はタグごとの特徴を返す。因子トークンの文法は旧 `sample_search` の `specs` と同じ
（`msdial/sample_factors.py`。要素内 AND・順不同、要素間 OR）。fixture はテスト自身が作る。
"""
import json

import numpy as np
import pytest

from metabolomix.analysis import sample_manifest
from metabolomix.core import session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools.mztab_tools import dataset_status

SAMPLES = [
    "20220901_RAW_control_6h_1_NEG",
    "20220902_RAW_ILG_6h_1_NEG",
    "20220902_RAW_ILG_6h_2_NEG",
    "20220901_QC_RAW_NEG_1",
]


@pytest.fixture(autouse=True)
def reset_session():
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def _dataset(groups: dict[str, str] | None = None) -> DatasetState:
    ds = DatasetState()
    ds.feature_matrix = np.random.default_rng(3).uniform(100, 200, (3, len(SAMPLES)))
    ds.sample_names = list(SAMPLES)
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(len(SAMPLES))]
    ds.feature_ids = ["F0", "F1", "F2"]
    ds.quantification_measure = "peak_height"
    ds.validation_result = {"ok": True, "errors": [], "warnings": []}
    if groups:
        rows = sample_manifest.resolve_metadata(ds, None)
        for row in rows:
            row["group"] = groups.get(row["sample_id"])
        sample_manifest.apply_defaults_and_inputs(ds, rows)
    return ds


def _status(**kwargs) -> dict:
    return json.loads(dataset_status(**kwargs))


def test_summary_counts_groups_roles_tags():
    ds = _dataset({SAMPLES[0]: "control", SAMPLES[1]: "ILG", SAMPLES[2]: "ILG"})
    ds.feature_tags = {"F0": ["isotope"], "F1": ["isotope", "adduct"], "F2": []}
    session_state.session.dataset = ds

    status = _status()

    assert status["groups"] == {"control": 1, "ILG": 2}
    assert status["roles"] == {"sample": 3, "qc": 1}
    assert status["n_ungrouped"] == 1
    assert status["tags"] == {"isotope": 2, "adduct": 1}
    assert status["excluded"] == {"samples": 0, "features": 0}
    # 既存の鍵は変えない
    assert status["samples"].splitlines()[0] == "name\trole"
    assert status["n_samples"] == 4


def test_samples_view_lists_columns_and_excluded_reason():
    ds = _dataset({SAMPLES[0]: "control"})
    ds.excluded_samples = {SAMPLES[1]: "injection failed"}
    session_state.session.dataset = ds

    status = _status(view="samples")

    lines = status["samples"].splitlines()
    assert lines[0] == "name\tgroup\trole\tbatch\tinjection_order\tsubject_id\texcluded"
    by_name = {line.split("\t")[0]: line.split("\t") for line in lines[1:]}
    assert len(by_name) == 4
    assert by_name[SAMPLES[0]][1] == "control"
    assert by_name[SAMPLES[3]][2] == "qc"
    assert by_name[SAMPLES[1]][6] == "injection failed"
    assert by_name[SAMPLES[0]][6] == ""
    assert status["n_matched"] == 4


def test_samples_view_filters_by_factor_tokens():
    ds = _dataset()
    session_state.session.dataset = ds

    one = _status(view="samples", query="ILG_6h")
    names = [line.split("\t")[0] for line in one["samples"].splitlines()[1:]]
    assert names == [SAMPLES[1], SAMPLES[2]]
    assert one["n_matched"] == 2 and one["n_total"] == 4

    # 順不同の AND、要素間は OR
    either = _status(view="samples", query="6h_ILG,qc")
    names = [line.split("\t")[0] for line in either["samples"].splitlines()[1:]]
    assert names == [SAMPLES[1], SAMPLES[2], SAMPLES[3]]

    # QC も落とさない（role は列で判断する）
    raw = _status(view="samples", query="raw")
    assert raw["n_matched"] == 4

    # 完全なサンプル名でも引ける
    exact = _status(view="samples", query=SAMPLES[1])
    assert exact["n_matched"] == 1


def test_samples_view_zero_match_explains_with_vocabulary():
    session_state.session.dataset = _dataset()

    out = _status(view="samples", query="24h")

    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "24h" in out["error"]["message"]
    assert "6h" in out["error"]["details"]["token_vocabulary"]["tokens"]


def test_samples_view_uses_group_as_a_token():
    session_state.session.dataset = _dataset({SAMPLES[0]: "ctrlgrp"})

    out = _status(view="samples", query="ctrlgrp")

    assert [line.split("\t")[0] for line in out["samples"].splitlines()[1:]] == [SAMPLES[0]]


def test_tags_view_lists_features_per_tag_and_caps_at_fifty():
    ds = _dataset()
    ds.feature_ids = [f"F{i}" for i in range(60)]
    ds.feature_tags = {f"F{i}": ["isotope"] for i in range(60)}
    ds.feature_tags["F3"] = ["isotope", "adduct"]
    session_state.session.dataset = ds

    out = _status(view="tags")

    assert out["tags"]["isotope"]["n_features"] == 60
    assert len(out["tags"]["isotope"]["feature_ids"]) == 50
    assert out["tags"]["isotope"]["truncated"] is True
    assert out["tags"]["adduct"] == {"n_features": 1, "feature_ids": ["F3"], "truncated": False}

    only = _status(view="tags", query="adduct")
    assert list(only["tags"]) == ["adduct"]


def test_unknown_view_is_refused_and_missing_state_still_wins():
    out = json.loads(dataset_status(view="nope"))
    assert out["error"]["code"] == "missing_state"

    session_state.session.dataset = _dataset()
    bad = _status(view="nope")
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"
    assert _status(view="summary", query="x")["error"]["code"] == "DATASET_BAD_REQUEST"
