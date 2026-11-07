"""`dataset_exclude`（B-2b Task 3）: 試料・特徴の除外を可逆に持ち、前処理が当てる。

除外は解析の状態（`DatasetState.excluded_samples` / `excluded_features`。値は理由）で、
入力の指紋には入らない。変えたら `metadata_revision` を進め、前処理済み行列から先の結果を古くする。
"""
import json

import pytest

from metabolomix.analysis import dataset_service, result_state
from metabolomix.analysis.sample_manifest import apply_metadata
from metabolomix.core import session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.tools.sample_tools import dataset_exclude
from tests.pipeline_fixtures import make_dataset, metadata_rows


@pytest.fixture(autouse=True)
def reset_session():
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def _loaded(with_rows: bool = False):
    ds = make_dataset()
    if with_rows:
        apply_metadata(ds, metadata_rows(ds, n_qc=0))
    session_state.session.dataset = ds
    return ds


def _call(**kwargs) -> dict:
    return json.loads(dataset_exclude(**kwargs))


def test_exclude_is_reversible():
    ds = _loaded()

    added = _call(samples=["S0", "S1"], features=["F2"], reason="bad injection")
    assert added["status"] == "success"
    assert ds.excluded_samples == {"S0": "bad injection", "S1": "bad injection"}
    assert ds.excluded_features == {"F2": "bad injection"}
    assert added["n_excluded_samples"] == 2 and added["n_excluded_features"] == 1

    removed = _call(samples=["S0"], mode="remove")
    assert ds.excluded_samples == {"S1": "bad injection"}
    assert removed["n_excluded_samples"] == 1

    cleared_one_kind = _call(samples=[], mode="clear")
    assert ds.excluded_samples == {} and ds.excluded_features == {"F2": "bad injection"}
    assert cleared_one_kind["n_excluded_features"] == 1

    _call(samples=["S3"], reason="x")
    cleared = _call(mode="clear")
    assert ds.excluded_samples == {} and ds.excluded_features == {}
    assert cleared["n_excluded_samples"] == 0 and cleared["n_excluded_features"] == 0


def test_exclude_accepts_sample_id_and_stores_the_sample_name():
    ds = _loaded(with_rows=True)
    rows = ds.sample_metadata_rows
    rows[2]["sample_id"] = "ID-2"
    apply_metadata(ds, rows)

    _call(samples=["ID-2"], reason="r")

    assert ds.excluded_samples == {"S2": "r"}


def test_exclude_unknown_name_refuses_whole_request():
    ds = _loaded()
    before_revision = ds.metadata_revision

    out = _call(samples=["S0", "NOPE"], features=["F1"], reason="r")

    assert out["error"]["code"] == "DATASET_BAD_REQUEST"
    assert out["error"]["details"]["unknown_samples"] == ["NOPE"]
    assert ds.excluded_samples == {} and ds.excluded_features == {}
    assert ds.metadata_revision == before_revision

    bad_feature = _call(samples=["S0"], features=["Z9"], reason="r")
    assert bad_feature["error"]["details"]["unknown_features"] == ["Z9"]
    assert ds.excluded_samples == {}

    # remove も存在しない名前は全体を拒む
    _call(samples=["S0"], reason="r")
    refused = _call(samples=["S0", "NOPE"], mode="remove")
    assert refused["error"]["code"] == "DATASET_BAD_REQUEST"
    assert ds.excluded_samples == {"S0": "r"}


def test_exclude_requires_reason_for_add_and_something_to_act_on():
    ds = _loaded()

    assert _call(samples=["S0"])["error"]["code"] == "DATASET_BAD_REQUEST"
    assert _call(samples=["S0"], reason="  ")["error"]["code"] == "DATASET_BAD_REQUEST"
    assert _call(reason="r")["error"]["code"] == "DATASET_BAD_REQUEST"
    assert _call(samples=["S0"], mode="weird", reason="r")["error"]["code"] == "DATASET_BAD_REQUEST"
    assert ds.excluded_samples == {}


def test_exclude_without_dataset_is_missing_state():
    out = json.loads(dataset_exclude(samples=["S0"], reason="r"))
    assert out["error"]["code"] == "missing_state"
    assert "load_dataset" in out["error"]["required_tools"]


def test_exclude_add_is_idempotent_and_does_not_invalidate_when_nothing_changes():
    ds = _loaded()
    _call(samples=["S0"], reason="r")
    dataset_service.preprocess_dataset(ds, {"impute": "half_min"})
    revision = ds.metadata_revision
    preprocess_id = ds.preprocess_id

    again = _call(samples=["S0"], reason="r")

    assert again["changed"] is False
    assert ds.metadata_revision == revision
    assert ds.preprocess_id == preprocess_id and ds.pp_matrix is not None


def test_exclude_applied_by_preprocess():
    ds = _loaded()
    _call(samples=["S0", "S5"], features=["F2"], reason="r")

    dataset_service.preprocess_dataset(ds, {"impute": "half_min"})

    assert ds.pp_matrix.shape == (6, 5)
    assert "S0" not in ds.pp_sample_names and "S5" not in ds.pp_sample_names
    assert "F2" not in ds.pp_feature_names
    assert ds.sample_names[0] == "S0"  # 入力そのものは変えない
    assert len(ds.feature_ids) == 6


def test_exclude_applied_by_preprocess_with_explicit_metadata_and_detection_mask():
    import numpy as np

    ds = _loaded(with_rows=True)
    mask = np.ones((6, 8), dtype=bool)
    mask[0, 7] = False  # 除外する試料でだけ未検出
    ds.detected_mask = mask
    ds.feature_qc = {"source": "arf"}
    _call(samples=["S7"], features=["F4"], reason="r")

    dataset_service.preprocess_dataset(
        ds, {"impute": "none", "min_detection_rate": 1.0})

    assert "S7" not in ds.pp_sample_names
    assert "F0" in ds.pp_feature_names  # 除外した試料の未検出は分母に入らない
    assert "F4" not in ds.pp_feature_names


def test_exclude_marks_downstream_results_stale():
    """除外の後に、除外前の結果を「現在の結果」として通さない（Review Focus 3）。"""
    ds = _loaded()
    pre = dataset_service.preprocess_dataset(ds, {"impute": "half_min"})
    pca = dataset_service.pca_dataset(ds, 2)
    old_pca_id = pca["provenance"]["result_id"]
    assert result_state.is_current(ds, ds.results[old_pca_id])
    assert ds.last_pca is not None
    revision = ds.metadata_revision

    out = _call(samples=["S0"], reason="r")

    assert out["changed"] is True
    assert ds.metadata_revision == revision + 1
    assert ds.pp_matrix is None and ds.last_pca is None and ds.last_differential is None
    assert ds.preprocess_id is None
    # ds.results には残るが、古い結果として拒む
    assert old_pca_id in ds.results
    assert not result_state.is_current(ds, ds.results[old_pca_id])
    with pytest.raises(DomainError) as info:
        result_state.assert_current(ds, ds.results[old_pca_id])
    assert info.value.code == "STALE_ANALYSIS_RESULT"
    assert not result_state.is_current(ds, ds.results[pre["provenance"]["result_id"]])
    assert "dataset_preprocess" in out["next"]


def test_exclusions_do_not_change_the_input_fingerprint():
    ds = _loaded()
    before = result_state.dataset_fingerprint(ds)

    _call(samples=["S0"], features=["F1"], reason="r")

    assert result_state.dataset_fingerprint(ds) == before


def test_preprocess_fingerprint_separates_exclusion_states():
    ds = _loaded()
    recipe = {"impute": "half_min"}
    plain = result_state.preprocess_fingerprint(ds, recipe, None)
    ds.excluded_samples = {"S0": "r"}
    assert result_state.preprocess_fingerprint(ds, recipe, None) != plain
    ds.excluded_samples = {}
    assert result_state.preprocess_fingerprint(ds, recipe, None) == plain


def test_feature_search_says_when_a_named_feature_is_excluded():
    from metabolomix.tools.dataset_analysis_tools import dataset_feature_search

    ds = _loaded()
    _call(features=["F2"], reason="co-eluting artifact")

    out = json.loads(dataset_feature_search(feature_ids=["F2"]))

    notes = " ".join(out.get("notes") or [])
    assert "F2" in notes and "co-eluting artifact" in notes, out
