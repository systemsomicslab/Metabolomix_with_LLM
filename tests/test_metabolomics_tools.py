"""pipeline の v2 の入口（spec §12）が使う部品の境界。

公開ツール（旧 `dataset_statistic` / `dataset_build_matrix` と、v1 の 2 群比較ツール）は 2026-11 に消えた。
v2 の解析行列と統計は pipeline の工程だけが作る（工程のテストは `test_metabolomics_end_to_end.py`）。
ここに残すのは、ツールを経ずに今も効く 2 点——worker の handler が揃っていること、
行列の identity と「揃っていない前提を揃ったことにしない」規則（`matrix_state.make_matrix`）。
旧ツール固有だった前提欠落の封筒・サマリの形・未知の kind の拒否は、`test_statistics_v2.py`
（`run_statistic` の止まり方）と `test_analysis_matrix.py` が持つ。
"""
from __future__ import annotations

import numpy as np

from metabolomix.core import session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.pipeline.metabolomics_handlers import build_handlers


# ---------- brief記載のRED ----------

def test_worker_handlers_cover_new_stages():
    assert {"load_assay_evidence", "resolve_feature_bindings", "statistics"} \
        <= set(build_handlers())


# ---------- fixture ----------

def _rows():
    groups = ["control"] * 3 + ["treated"] * 3
    return [{"sample_id": f"s{i}", "source_file": f"S{i}.wiff", "role": "sample",
             "group": group, "batch": "B1", "injection_order": i + 1,
             "qc_pool": None, "include": True, "biological_sample_id": f"bio{i}"}
            for i, group in enumerate(groups)]


def _dataset_for_matrix() -> DatasetState:
    ds = DatasetState()
    ds.feature_ids = ["1", "2"]
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(6)]
    ds.sample_names = [f"S{i}" for i in range(6)]
    ds.sample_metadata_rows = _rows()
    # DatasetState は (feature × assay)。解析行列側で転置される。
    ds.feature_matrix = np.array([[10., 11., 9., 40., 42., 38.],
                                  [5., 6., 4., 20., 21., 19.]])
    session_state.session.dataset = ds
    return ds


def _recipe(**overrides) -> dict:
    recipe = {"base": "peak_height", "normalize": "none", "drift_correct": False,
              "filter": None, "impute": "none"}
    recipe.update(overrides)
    return recipe


def _make(ds, recipe):
    from metabolomix.analysis.matrix_state import make_matrix

    return make_matrix(ds, recipe, {}, {}, np.ones(len(ds.feature_ids), dtype=bool))


# ---------- v1 の既定値を変えない ----------

def test_v1_compare_dataset_defaults_are_unchanged():
    """pipeline が今も呼ぶ v1 の 2 群比較の既定値（対話ツールは統計エンジンへ移った）。"""
    import inspect

    from metabolomix.analysis.dataset_service import compare_dataset

    signature = inspect.signature(compare_dataset)
    assert signature.parameters["q_threshold"].default == 0.05
    assert signature.parameters["log2fc_threshold"].default == 1.0


# ---------- 解析行列（pipeline の preprocess 工程が作る形） ----------

def test_the_same_recipe_and_dataset_give_the_same_matrix_id():
    ds = _dataset_for_matrix()
    assert _make(ds, _recipe())["matrix_id"] == _make(ds, _recipe())["matrix_id"]


def test_detection_filter_is_not_evaluable_when_detection_is_unknown():
    ds = _dataset_for_matrix()
    matrix = _make(ds, _recipe(filter={"min_detection_rate": 0.8}))
    steps = {step["step"]: step for step in matrix["correction_history"]}
    assert steps["detection_filter"]["status"] == "not_evaluable"
    # 検出状態が分からないことを理由に feature を落とさない。
    assert int(sum(bool(x) for x in matrix["eligibility_mask"])) == 2
    assert "detection_unknown" in matrix["caveats"]
