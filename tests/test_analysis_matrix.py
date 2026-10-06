"""版管理した解析行列（spec §8.2 `analysis-matrix.v1`）。

統計は「直近の前処理結果」ではなく**matrix result ID**を参照する。同じ
データセットから recipe 違いの行列が2本出るのが前提だからで、`ds.pp_matrix`
のような単一スロットへ上書きすると、図と表が別 recipe の数字を並べても
どちらも「最新」を名乗れてしまう。

ここで縛るのは4点:

1. `make_matrix` は**補完しない**。補完は QC 後 filter を確定させてからでないと、
   QC で落ちる予定の試料・featureの値を混ぜ込む。補完は `finalize_matrix` だけ。
2. 解析対象から外れた feature（検出率filter）も**列としては残す**。内部標準は
   解析対象外でも support として必要で、削除すると比の分母が消える。
   「対象かどうか」は eligibility mask が持ち、値の有無とは別の軸。
3. 分母不正で欠損したセル（locked）は補完後も欠損のまま。
4. 保存した行列は、値hash・軸順・mask形状を load 時に確認する。1バイト変われば
   読み込みを拒否する——黙って別の数字を統計へ渡さない。
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from metabolomix.analysis.matrix_state import (
    MATRIX_SCHEMA,
    finalize_matrix,
    load_matrix,
    make_matrix,
    save_matrix,
)
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState


def _recipe(**overrides) -> dict:
    recipe = {"base": "peak_height", "normalize": "none", "drift_correct": False,
              "filter": None, "impute": "none"}
    recipe.update(overrides)
    return recipe


def _dataset() -> DatasetState:
    """assay 3 × feature 3。列0=対象、列1=内部標準、列2=無関係。"""
    ds = DatasetState()
    # feature_matrix は (features, assays)。解析行列は (assays, features)。
    ds.feature_matrix = np.array([[10., 20., 30.],
                                  [2., 4., 5.],
                                  [100., 200., np.nan]])
    ds.feature_ids = ["1", "2", "3"]
    ds.sample_assay_ids = ["assay[1]", "assay[2]", "assay[3]"]
    ds.sample_names = ["S1", "S2", "S3"]
    ds.quantification_measure = "peak_height"
    return ds


def _bindings(ds, pairs=()) -> dict:
    return {"schema": "feature-bindings.v1", "dataset_id": ds.dataset_id,
            "dataset_hash": "h" * 64, "rule_hash": "r" * 64, "status": "resolved",
            "bindings": {}, "unresolved": [],
            "internal_standard_map": {"schema": "internal-standard-map.v1",
                                      "dataset_id": ds.dataset_id,
                                      "dataset_hash": "h" * 64,
                                      "rule_hash": "r" * 64,
                                      "pairs": list(pairs)}}


def _pair() -> dict:
    return {"target_id": "gaba", "target_feature_id": "1",
            "standard_target_id": "d6", "standard_feature_id": "2"}


def _evidence(ds, *, availability=True, undetected=()) -> dict:
    rows = []
    for fid in ds.feature_ids:
        for assay_id in ds.sample_assay_ids:
            status = ("gap_filled" if (fid, assay_id) in undetected else "detected")
            rows.append({"feature_id": fid, "assay_id": assay_id,
                         "observed_rt_min": 1.0, "observed_mz": 100.0,
                         "detection_status": status})
    return {"schema": "assay-feature-evidence.v1", "dataset_id": ds.dataset_id,
            "availability": availability, "reasons": [],
            "rows": rows if availability else []}


def _all_eligible(ds) -> np.ndarray:
    return np.ones(len(ds.feature_ids), dtype=bool)


# ---------- make_matrix ----------

def test_peak_height_matrix_is_the_transposed_raw_matrix():
    ds = _dataset()
    result = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds),
                         _all_eligible(ds))
    assert result["schema"] == MATRIX_SCHEMA
    assert result["assay_ids"] == ds.sample_assay_ids
    assert result["feature_ids"] == ds.feature_ids
    assert result["values"][:, 0].tolist() == [10., 20., 30.]
    assert result["units"] == ["peak_height"] * 3
    assert result["stage"] == "preprocessed"
    assert result["parent_id"] is None


def test_make_matrix_does_not_touch_the_dataset():
    """recipeごとに別の行列を作る。ds.pp_matrix を上書きしない。"""
    ds = _dataset()
    make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    assert ds.pp_matrix is None
    assert ds.preprocessing_recipe == {}


def test_make_matrix_never_imputes():
    ds = _dataset()
    result = make_matrix(ds, _recipe(impute="half_min"), _bindings(ds),
                         _evidence(ds), _all_eligible(ds))
    assert np.isnan(result["values"][2, 2])
    assert not result["imputed_mask"].any()


def test_two_recipes_give_two_matrix_ids():
    ds = _dataset()
    a = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    b = make_matrix(ds, _recipe(normalize="median"), _bindings(ds),
                    _evidence(ds), _all_eligible(ds))
    assert a["matrix_id"] != b["matrix_id"]


def test_internal_standard_ratio_is_a_different_matrix_than_raw():
    ds = _dataset()
    raw = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    ratio = make_matrix(ds, _recipe(base="internal_standard_ratio"),
                        _bindings(ds, [_pair()]), _evidence(ds), _all_eligible(ds))
    assert raw["matrix_id"] != ratio["matrix_id"]
    assert ratio["units"][0] == "internal_standard_ratio"
    assert raw["units"][0] == "peak_height"
    assert ratio["values"][:, 0].tolist() == [5., 5., 6.]


def test_undetected_standard_locks_the_ratio_cell():
    ds = _dataset()
    result = make_matrix(ds, _recipe(base="internal_standard_ratio"),
                         _bindings(ds, [_pair()]),
                         _evidence(ds, undetected=[("2", "assay[2]")]),
                         _all_eligible(ds))
    assert np.isnan(result["values"][1, 0])
    assert result["locked_mask"][1, 0]


def test_ratio_without_evidence_records_that_detection_is_unknown():
    ds = _dataset()
    result = make_matrix(ds, _recipe(base="internal_standard_ratio"),
                         _bindings(ds, [_pair()]),
                         _evidence(ds, availability=False), _all_eligible(ds))
    assert "detection_unknown" in result["caveats"]


def test_normalizing_a_ratio_matrix_is_refused():
    """内部標準比への TIC/median/PQN は二重正規化。初期版では拒否する。"""
    ds = _dataset()
    with pytest.raises(DomainError) as caught:
        make_matrix(ds, _recipe(base="internal_standard_ratio", normalize="pqn"),
                    _bindings(ds, [_pair()]), _evidence(ds), _all_eligible(ds))
    assert caught.value.code == "MATRIX_RECIPE_INVALID"


def test_drift_correction_without_qc_prerequisites_stops():
    ds = _dataset()
    with pytest.raises(DomainError) as caught:
        make_matrix(ds, _recipe(drift_correct=True), _bindings(ds),
                    _evidence(ds), _all_eligible(ds))
    assert caught.value.code == "QC_PREREQUISITE_MISSING"


# ---------- filter は列を消さない ----------

def test_detection_filter_narrows_eligibility_but_keeps_the_column():
    ds = _dataset()
    undetected = [("3", a) for a in ds.sample_assay_ids]
    result = make_matrix(ds, _recipe(filter={"min_detection_rate": 0.8}),
                         _bindings(ds), _evidence(ds, undetected=undetected),
                         _all_eligible(ds))
    assert result["values"].shape == (3, 3)
    assert result["feature_ids"] == ["1", "2", "3"]
    assert result["eligibility_mask"].tolist() == [True, True, False]


def test_a_filtered_internal_standard_stays_available_as_support():
    """検出filterで解析対象から外れた内部標準も、support行列からは削除しない。"""
    ds = _dataset()
    undetected = [("2", a) for a in ds.sample_assay_ids]
    result = make_matrix(ds, _recipe(base="internal_standard_ratio",
                                     filter={"min_detection_rate": 0.8}),
                         _bindings(ds, [_pair()]),
                         _evidence(ds, undetected=undetected), _all_eligible(ds))
    index = result["feature_ids"].index("2")
    assert result["eligibility_mask"][index] is np.False_
    assert result["values"][:, index].tolist() == [2., 4., 5.]
    assert "2" in result["support_feature_ids"]


# ---------- finalize_matrix ----------

def test_finalize_imputes_and_links_to_its_parent():
    ds = _dataset()
    base = make_matrix(ds, _recipe(impute="half_min"), _bindings(ds),
                       _evidence(ds), _all_eligible(ds))
    final = finalize_matrix(base, _all_eligible(ds), "half_min")
    assert final["parent_id"] == base["matrix_id"]
    assert final["matrix_id"] != base["matrix_id"]
    assert final["stage"] == "finalized"
    assert not np.isnan(final["values"][2, 2])
    assert final["imputed_mask"][2, 2]


def test_locked_cells_stay_missing_after_imputation():
    """分母が使えず欠損したセルは、補完器を通しても欠損のまま残す。"""
    ds = _dataset()
    base = make_matrix(ds, _recipe(base="internal_standard_ratio"),
                       _bindings(ds, [_pair()]),
                       _evidence(ds, undetected=[("2", "assay[3]")]),
                       _all_eligible(ds))
    assert base["locked_mask"][2, 0]
    assert np.isnan(base["values"][2, 0])

    final = finalize_matrix(base, _all_eligible(ds), "half_min")
    assert np.isnan(final["values"][2, 0])
    assert not final["imputed_mask"][2, 0]
    assert final["missing_reasons"]["1"]["assay[3]"]["code"] \
        == "standard_denominator_invalid"


def test_finalize_narrows_eligibility_with_the_qc_mask():
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    qc_mask = np.array([True, False, True])
    final = finalize_matrix(base, qc_mask, "none")
    assert final["eligibility_mask"].tolist() == [True, False, True]
    assert final["values"].shape == (3, 3)      # 列は消さない


def test_finalize_does_not_impute_ineligible_columns():
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    final = finalize_matrix(base, np.array([True, True, False]), "half_min")
    assert np.isnan(final["values"][2, 2])
    assert not final["imputed_mask"].any()


# ---------- 保存と読み戻し ----------

def test_save_and_load_round_trip(tmp_path):
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    reference = save_matrix(base, tmp_path)
    loaded = load_matrix(reference, tmp_path)

    assert loaded["matrix_id"] == base["matrix_id"]
    assert loaded["feature_ids"] == base["feature_ids"]
    assert loaded["assay_ids"] == base["assay_ids"]
    assert np.allclose(loaded["values"], base["values"], equal_nan=True)
    assert loaded["eligibility_mask"].tolist() == base["eligibility_mask"].tolist()


def test_tampered_values_are_refused(tmp_path):
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    reference = save_matrix(base, tmp_path)

    arrays_path = tmp_path / reference["arrays_file"]
    with np.load(arrays_path, allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    arrays["values"][0, 0] = 999.0
    np.savez(arrays_path, **arrays)

    with pytest.raises(DomainError) as caught:
        load_matrix(reference, tmp_path)
    assert caught.value.code == "MATRIX_INTEGRITY_MISMATCH"


def test_tampered_axis_order_is_refused(tmp_path):
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    reference = save_matrix(base, tmp_path)

    meta_path = tmp_path / reference["meta_file"]
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["feature_ids"] = list(reversed(meta["feature_ids"]))
    meta_path.write_text(json.dumps(meta), encoding="utf-8")

    with pytest.raises(DomainError) as caught:
        load_matrix(reference, tmp_path)
    assert caught.value.code == "MATRIX_INTEGRITY_MISMATCH"


def test_missing_file_is_reported_not_crashed(tmp_path):
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    reference = save_matrix(base, tmp_path)
    (tmp_path / reference["arrays_file"]).unlink()
    with pytest.raises(DomainError) as caught:
        load_matrix(reference, tmp_path)
    assert caught.value.code == "MATRIX_NOT_FOUND"


def test_saved_arrays_do_not_allow_pickle(tmp_path):
    ds = _dataset()
    base = make_matrix(ds, _recipe(), _bindings(ds), _evidence(ds), _all_eligible(ds))
    reference = save_matrix(base, tmp_path)
    with np.load(tmp_path / reference["arrays_file"], allow_pickle=False) as data:
        assert set(data.files) >= {"values", "eligibility_mask", "locked_mask",
                                   "imputed_mask"}


# ---------- 無効化 ----------

def test_rebinding_invalidates_the_preprocessed_matrix():
    """どのfeatureが内部標準かが変われば、比の行列は作り直しになる。"""
    from metabolomix.analysis.result_state import invalidate_results

    ds = _dataset()
    ds.pp_matrix = object()
    ds.last_pca = {"x": 1}
    invalidate_results(ds, {"bindings"})
    assert ds.pp_matrix is None
    assert ds.last_pca is None


def test_matrix_recipe_change_invalidates_the_preprocessed_matrix():
    from metabolomix.analysis.result_state import invalidate_results

    ds = _dataset()
    ds.pp_matrix = object()
    invalidate_results(ds, {"matrix_recipe"})
    assert ds.pp_matrix is None


@pytest.mark.parametrize("role,include", [("standard", True), ("blank", True), ("sample", False), ("qc", False)])
def test_pqn_and_detection_filter_ignore_nonstudy_reference_rows(role, include):
    ds = _dataset()
    ds.feature_matrix = np.array([[1., 2., 100.], [2., 4., 1.], [3., 6., 50.]])
    ds.sample_metadata_rows = [
        {"sample_id": "a", "role": "sample", "include": True},
        {"sample_id": "b", "role": "sample", "include": True},
        {"sample_id": "c", "role": role, "include": include}]
    out = make_matrix(ds, _recipe(normalize="pqn", filter={"min_detection_rate": 1.}),
                      _bindings(ds), _evidence(ds, undetected=[("1", "assay[3]")]), _all_eligible(ds))
    np.testing.assert_allclose(out["values"][:2], [[1.5, 3., 4.5], [1.5, 3., 4.5]])
    assert out["eligibility_mask"].all()
    assert out["values"].shape == (3, 3)
    assert out["units"] == ["normalized_height"] * 3


@pytest.mark.parametrize("defect", ["excluded", "batch", "pool", "missing_pool", "missing_order", "outside"])
def test_drift_prerequisites_fail_closed(defect):
    ds = _dataset()
    ds.feature_matrix = np.ones((3, 5))
    ds.sample_assay_ids = [f"assay[{i}]" for i in range(5)]
    ds.sample_metadata_rows = [dict(sample_id=str(i), role="qc" if i < 4 else "sample",
                                   include=True, batch="B", qc_pool="P" if i < 4 else None,
                                   injection_order=i * 2 + 1 if i < 4 else 4) for i in range(5)]
    if defect == "excluded": ds.sample_metadata_rows[0]["include"] = False
    if defect == "batch": ds.sample_metadata_rows[0]["batch"] = "other"
    if defect == "pool": ds.sample_metadata_rows[0]["qc_pool"] = "other"
    if defect == "missing_pool": ds.sample_metadata_rows[0]["qc_pool"] = None
    if defect == "missing_order": ds.sample_metadata_rows[4]["injection_order"] = None
    if defect == "outside": ds.sample_metadata_rows[4]["injection_order"] = 99
    with pytest.raises(DomainError) as caught:
        make_matrix(ds, _recipe(drift_correct=True), _bindings(ds), _evidence(ds), _all_eligible(ds))
    assert caught.value.code == "QC_PREREQUISITE_MISSING"
