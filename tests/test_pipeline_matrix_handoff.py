"""pipeline が作った解析行列を、対話セッションへ引き渡す経路（S2）。

worker は別プロセスで走るので、行列は `context["runtime"]` にしか無く、
`_persist` は values を落として保存していた——run が終わると数値は残らない。
ここで縛るのは2つ:

- **保存側**: 行列は `save_matrix`（npz + meta hash + 配列 fingerprint）で
  results/matrices/ へ残り、参照が run record 側に載る。
- **読み側**: `dataset_load(pipeline_path=...)` が、その run の dataset・
  解決済み試料対応表・binding・行列を**まとめて**セッションへ載せる。
  行列だけを別ツールで載せられると、run A の dataset に run B の行列を
  混ぜられてしまう（群の対応が黙ってずれる）。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from metabolomix.core import session_state
from metabolomix.pipeline import store
from tests.metabolomics_fixtures import MetabolomicsHarness
from tests.test_mztab_tools import _make_job_json


def _payload(text: str) -> dict:
    return json.loads(text)


def _matrix_refs(record: dict) -> list[dict]:
    return [r for r in (record.get("results") or [])
            if r.get("output_name") == "matrix"]


def _result_data(record: dict, ref: dict) -> dict:
    path = Path(record["identity"]["pipeline_root"]) / ref["relative_path"]
    return json.loads(path.read_text(encoding="utf-8"))["data"]


def _completed_run(tmp_path) -> tuple[MetabolomicsHarness, Path]:
    """harness で v2 run を1本回し、`console_job_path` を実 job.json へ向ける。

    harness の fake Console は console_job_path に mzTab を入れる。読み側は
    本番と同じく `analysis-job.json` として読むので、ここで実物へ差し替える。
    """
    harness = MetabolomicsHarness(tmp_path)
    harness.run()
    job_path, _ = _make_job_json(tmp_path, "../../source/Height_synthetic.mzTab")
    record = store.load_run(harness.pipeline_path)
    record["upstream"]["console_job_path"] = str(job_path)
    store.save_run(harness.pipeline_path, record,
                   expected_revision=record["state_revision"])
    return harness, job_path


# ---------- 保存側 ----------

def test_matrices_are_saved_with_their_values_not_only_a_summary(tmp_path):
    harness = MetabolomicsHarness(tmp_path)
    harness.run()
    record = harness.record()
    refs = _matrix_refs(record)
    assert refs, "matrix 成果物が1つも無い"

    root = Path(record["identity"]["pipeline_root"])
    for ref in refs:
        storage = _result_data(record, ref).get("storage")
        assert storage is not None, "matrix 成果物に storage 参照が無い"
        assert (root / "results" / "matrices" / storage["arrays_file"]).is_file()
        assert (root / "results" / "matrices" / storage["meta_file"]).is_file()


def test_saved_matrices_can_be_read_back_with_the_recorded_reference(tmp_path):
    from metabolomix.analysis.matrix_state import load_matrix

    harness = MetabolomicsHarness(tmp_path)
    harness.run()
    record = harness.record()
    root = Path(record["identity"]["pipeline_root"])

    ref = _matrix_refs(record)[0]
    storage = _result_data(record, ref)["storage"]
    matrix = load_matrix(storage, root / "results" / "matrices")
    assert matrix["values"].shape == (len(matrix["assay_ids"]),
                                      len(matrix["feature_ids"]))


# ---------- 読み側 ----------

def test_loading_a_run_brings_the_dataset_and_its_matrices_together(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    session_state.session.dataset = None

    # 成功時は要約テキスト（エラー時だけ JSON 封筒）。
    summary = dataset_load(pipeline_path=str(harness.pipeline_path))
    assert "解析行列" in summary, summary

    ds = session_state.session.dataset
    assert ds is not None
    recorded = {_result_data(harness.record(), ref)["matrix_id"]
                for ref in _matrix_refs(harness.record())}
    assert recorded <= set(ds.analysis_matrices)


def test_loading_a_run_restores_the_resolved_sample_metadata(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    session_state.session.dataset = None
    dataset_load(pipeline_path=str(harness.pipeline_path))

    rows = session_state.session.dataset.sample_metadata_rows
    expected = harness.result_data("sample_manifest")["rows"]
    assert rows == expected
    # 群が入っていなければ、載せた行列で統計は組めない。
    assert {row.get("group") for row in rows if row.get("role") == "sample"} != {None}


def test_loading_a_run_restores_the_resolved_binding_targets(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    session_state.session.dataset = None
    dataset_load(pipeline_path=str(harness.pipeline_path))

    bindings = harness.result_data("feature_bindings") or {}
    expected = {tid: b["selected_feature_id"]
                for tid, b in (bindings.get("bindings") or {}).items()
                if b.get("status") == "resolved"}
    assert session_state.session.dataset.feature_binding_targets == expected


def test_a_restored_matrix_can_be_used_by_dataset_statistic(tmp_path):
    from metabolomix.tools.dataset_analysis_tools import dataset_statistic
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    session_state.session.dataset = None
    dataset_load(pipeline_path=str(harness.pipeline_path))

    matrix_id = sorted(session_state.session.dataset.analysis_matrices)[0]
    result = _payload(dataset_statistic(
        {"statistic_id": "check", "kind": "pca", "matrix_recipe_id": "default",
         "transform": "none", "feature_scope": {"mode": "all_eligible"},
         "scaling": "autoscale", "n_components": 2},
        matrix_id))
    assert result.get("status") == "success", result


def test_pipeline_path_cannot_be_combined_with_the_other_entry_points(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    payload = _payload(dataset_load(mztab_path=str(harness.mztab_path),
                                    pipeline_path=str(harness.pipeline_path)))
    assert payload["error"]["code"] == "DATASET_BAD_REQUEST"


def test_a_tampered_matrix_file_is_refused_and_leaves_the_session_alone(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _completed_run(tmp_path)
    record = harness.record()
    root = Path(record["identity"]["pipeline_root"])
    storage = _result_data(record, _matrix_refs(record)[0])["storage"]
    arrays_path = root / "results" / "matrices" / storage["arrays_file"]
    with np.load(arrays_path, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    arrays["values"] = arrays["values"] * 2.0
    np.savez(arrays_path, **arrays)

    session_state.session.dataset = None
    payload = _payload(dataset_load(pipeline_path=str(harness.pipeline_path)))
    assert payload["error"]["code"] == "MATRIX_INTEGRITY_MISMATCH"
    assert session_state.session.dataset is None
