"""`dataset_sample_peaks`: 特徴 × 試料のピークの RT・m/z・高さ・検出の印（読み取り）。

fixture はこのテストが tmp に作る。RT は試料ごとに違う値を `.arf` に書き、読んだ値が合うことを見る。
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from metabolomix.core import handles, session_state
from metabolomix.dataset.opening import open_dataset
from metabolomix.quant.derive import derive
from metabolomix.tools import reanalysis_tools as R
from tests.msdial_folder_fixtures import write_matching_mztab, write_msdial_folder


@pytest.fixture(autouse=True)
def fresh_state():
    session_state.session = session_state.AnalysisSession()
    handles.registry = handles.HandleRegistry()
    yield
    for handle in handles.registry.list():
        handles.registry.close(handle.handle_id)
    handles.registry = handles.HandleRegistry()
    session_state.session = session_state.AnalysisSession()


def _spec() -> dict:
    samples = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"}
               for i in range(3)]
    spots = [
        {"name": "PC 34:1", "mz": 760.5851, "heights": [100.0, 0.0, 50.0],
         "rts": [12.01, 12.50, 12.07], "gap_filled": [False, True, False]},
        {"name": "Unknown", "mz": 500.1, "heights": [10.0, 20.0, 30.0],
         "rts": [3.30, 3.31, 3.35]},
    ]
    return {"spots": spots, "samples": samples}


def _open_folder(tmp_path):
    write_msdial_folder(tmp_path, **_spec())
    ds = open_dataset(directory=str(tmp_path)).dataset
    session_state.session.dataset = ds
    return ds


def _rows(out: str) -> list[dict]:
    payload = json.loads(out)
    assert "error" not in payload, payload
    return payload["rows"]


def test_rt_matches_arf_and_gap_fill_is_not_detected(tmp_path):
    _open_folder(tmp_path)
    rows = _rows(R.dataset_sample_peaks(["0"]))
    assert [r["sample"] for r in rows] == ["s0", "s1", "s2"]
    assert {r["feature_id"] for r in rows} == {"0"}
    assert [r["rt"] for r in rows] == pytest.approx([12.01, 12.50, 12.07])
    assert [r["detected"] for r in rows] == [True, False, True]
    assert [r["height"] for r in rows] == pytest.approx([100.0, 0.0, 50.0])
    assert rows[0]["mz"] == pytest.approx(760.5851)
    assert rows[1]["mz"] is None   # gap-fill のセルに実測 m/z は無い


def test_int_and_str_ids_and_order_of_request(tmp_path):
    _open_folder(tmp_path)
    rows = _rows(R.dataset_sample_peaks([1, "0"]))
    assert [r["feature_id"] for r in rows] == ["1"] * 3 + ["0"] * 3
    assert [r["rt"] for r in rows[:3]] == pytest.approx([3.30, 3.31, 3.35])


def test_unknown_feature_is_an_error(tmp_path):
    _open_folder(tmp_path)
    out = json.loads(R.dataset_sample_peaks(["0", "999"]))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "999" in json.dumps(out["error"]["details"])


def test_more_than_200_features_is_an_error(tmp_path):
    _open_folder(tmp_path)
    out = json.loads(R.dataset_sample_peaks([str(i) for i in range(201)]))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "200" in out["error"]["message"]
    assert len(_rows(R.dataset_sample_peaks(["0"] * 1))) == 3


def test_empty_list_is_an_error(tmp_path):
    _open_folder(tmp_path)
    assert json.loads(R.dataset_sample_peaks([]))["error"]["code"] == "DATASET_BAD_REQUEST"


def test_mztab_dataset_reports_per_sample_rt_unavailable(tmp_path):
    write_matching_mztab(tmp_path / "x.mzTab", _spec())
    ds = open_dataset(mztab_path=str(tmp_path / "x.mzTab")).dataset
    session_state.session.dataset = ds
    out = json.loads(R.dataset_sample_peaks(["0"]))
    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert "per_sample_rt" in out["error"]["message"]


def test_no_dataset_is_a_missing_state_envelope():
    out = json.loads(R.dataset_sample_peaks(["0"]))
    assert "error" in out or "missing_state" in json.dumps(out)


def test_rt_array_is_compact_and_follows_derive(tmp_path):
    ds = _open_folder(tmp_path)
    assert isinstance(ds.per_sample_rt, np.ndarray)
    assert ds.per_sample_rt.shape == ds.feature_matrix.shape
    sub = derive(ds, matrix=ds.feature_matrix[[1]][:, [2, 0]],
                 detected=ds.detected_mask[[1]][:, [2, 0]],
                 feature_ids=["1"], sample_names=["s2", "s0"], state={})
    assert sub.per_sample_rt[0].tolist() == pytest.approx([3.35, 3.30])
    session_state.session.dataset = sub
    rows = _rows(R.dataset_sample_peaks(["1"]))
    assert [(r["sample"], r["rt"]) for r in rows] == [("s2", pytest.approx(3.35)),
                                                       ("s0", pytest.approx(3.30))]
