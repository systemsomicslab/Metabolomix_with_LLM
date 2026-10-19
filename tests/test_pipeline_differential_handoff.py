"""pipeline が作った差次的結果を、対話セッションへ引き渡す経路（spec 2026-10-19 §11）。

`dataset_load(pipeline_path=...)` は run の dataset を作り直すので、`dataset_id` は
run 側と必ず食い違い、session の前処理も通っていない。ここで縛るのは:

- 指紋（データ・群）が一致する差次的結果だけを、**result_id を付け替えずに**
  `ds.results` へ戻す（pipeline の TSV のメタ行と同じ ID で引けるように）。
- 戻せなかったもの（旧 run・データ違い・群違い）は理由つきで報告する。
- 同じ比較が何度も永続化されていても、**今の ref だけ**を読む
  （上書きされた古い ref の hash を「壊れた」と読まない）。
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from metabolomix.analysis import dataset_service, result_state
from metabolomix.core import session_state
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import report, session_handoff, store
from tests.test_pipeline_matrix_handoff import _completed_run

_RECIPE = {"normalize": "none", "impute": "half_min"}
_COMPARISON = {"comparison_id": "c1", "reference_group": "g1", "test_group": "g2"}


def _run_rows(harness) -> list[dict]:
    return list(harness.result_data("sample_manifest")["rows"])


def _differential_from_same_job(job_path: Path, rows: list[dict]) -> dict:
    """run と同じ job から ds を作り、pipeline と同じ関数で差次的結果を作る。"""
    ds = mztab_loading.load_dataset_state(job_path=job_path)
    ds.sample_metadata_rows = list(rows)
    dataset_service.preprocess_dataset(ds, _RECIPE)
    return dataset_service.run_comparison(ds, _COMPARISON, rows)


def _append_refs(harness, *refs: dict) -> None:
    record = store.load_run(harness.pipeline_path)
    record["results"] = list(record.get("results") or []) + list(refs)
    store.save_run(harness.pipeline_path, record,
                   expected_revision=record["state_revision"])


def _persist_differential(harness, result: dict, cid: str = "c1") -> dict:
    return report.persist_result(harness.pipeline_path, {
        "output_name": f"differential:{cid}", "kind": "differential",
        "result_id": result["provenance"]["result_id"], "data": result,
        "parent_ids": list(result["provenance"].get("parent_ids") or []),
    })


def _run_with_differential(tmp_path, *, mutate=None):
    harness, job_path = _completed_run(tmp_path)
    result = _differential_from_same_job(job_path, _run_rows(harness))
    if mutate is not None:
        result = mutate(copy.deepcopy(result))
    _append_refs(harness, _persist_differential(harness, result))
    return harness, result


def test_loading_a_run_restores_its_differential_results(tmp_path):
    harness, result = _run_with_differential(tmp_path)
    result_id = result["provenance"]["result_id"]

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    restored = ds.results[result_id]
    assert restored["provenance"]["result_id"] == result_id
    assert restored["provenance"]["handoff"]["comparison_id"] == "c1"
    assert restored["provenance"]["handoff"]["pipeline_root"] == str(harness.pipeline_path)
    assert result_state.is_current(ds, restored) is True
    # 「直近の結果」は session で計算したものに限る。
    assert ds.last_differential is None
    assert ds.handoff_report["restored"] == [
        {"comparison_id": "c1", "result_id": result_id,
         "group_a": "g1", "group_b": "g2"}]
    assert ds.handoff_report["rejected"] == []


def test_legacy_result_without_dataset_hash_is_reported_not_restored(tmp_path):
    def drop_hash(result):
        result["provenance"].pop("dataset_hash")
        return result

    harness, result = _run_with_differential(tmp_path, mutate=drop_hash)
    result_id = result["provenance"]["result_id"]

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    assert result_id not in ds.results
    assert ds.handoff_report["restored"] == []
    assert ds.handoff_report["rejected"] == [
        {"comparison_id": "c1", "result_id": result_id,
         "reason": "legacy_result_without_dataset_hash"}]


def test_result_from_other_data_is_reported_not_restored(tmp_path):
    def other_data(result):
        result["provenance"]["dataset_hash"] = "0" * 64
        return result

    harness, result = _run_with_differential(tmp_path, mutate=other_data)

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    assert result["provenance"]["result_id"] not in ds.results
    assert [r["reason"] for r in ds.handoff_report["rejected"]] == ["dataset_mismatch"]


def test_result_from_other_groups_is_reported_not_restored(tmp_path):
    harness, result = _run_with_differential(tmp_path)

    # 結果を保存したあとで、run の試料対応表を群違いで再保存する。
    rows = _run_rows(harness)
    swapped = [dict(row, group={"g1": "g2", "g2": "g1"}.get(row.get("group"),
                                                           row.get("group")))
               for row in rows]
    manifest_ref = report.persist_result(harness.pipeline_path, {
        "output_name": "sample_manifest", "kind": "sample_manifest",
        "result_id": "res_sample_manifest_regrouped", "data": {"rows": swapped},
    })
    _append_refs(harness, manifest_ref)

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    assert result["provenance"]["result_id"] not in ds.results
    assert ds.handoff_report["rejected"] == [
        {"comparison_id": "c1", "result_id": result["provenance"]["result_id"],
         "reason": "group_mismatch"}]


def test_handoff_uses_only_current_refs_after_reexport(tmp_path):
    harness, job_path = _completed_run(tmp_path)
    rows = _run_rows(harness)
    first = _differential_from_same_job(job_path, rows)
    first_ref = _persist_differential(harness, first)
    # export が同じ比較を再計算して、もう一度 `differential:c1` として永続化した。
    second = _differential_from_same_job(job_path, rows)
    second_ref = _persist_differential(harness, second)
    _append_refs(harness, first_ref, second_ref)
    # 同じ出力先を再利用していた版では、古い ref の hash はもうファイルと一致しない。
    stale = Path(harness.pipeline_path) / first_ref["relative_path"]
    stale.write_text(json.dumps({"overwritten": True}), encoding="utf-8")

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    assert second["provenance"]["result_id"] in ds.results
    assert first["provenance"]["result_id"] not in ds.results
    assert [r["result_id"] for r in ds.handoff_report["restored"]] == [
        second["provenance"]["result_id"]]


def test_a_tampered_current_differential_is_refused(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, _ = _run_with_differential(tmp_path)
    ref = [r for r in store.load_run(harness.pipeline_path)["results"]
           if r.get("output_name") == "differential:c1"][-1]
    (Path(harness.pipeline_path) / ref["relative_path"]).write_text(
        json.dumps({"tampered": True}), encoding="utf-8")

    session_state.session.dataset = None
    payload = json.loads(dataset_load(pipeline_path=str(harness.pipeline_path)))
    assert payload["error"]["code"] == "RESULT_INTEGRITY_MISMATCH"
    assert session_state.session.dataset is None


def test_a_run_without_differentials_reports_nothing(tmp_path):
    harness, _ = _completed_run(tmp_path)

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    assert ds.handoff_report == {"restored": [], "rejected": []}


def test_dataset_load_summary_lists_restored_and_rejected(tmp_path):
    from metabolomix.tools.mztab_tools import dataset_load

    harness, job_path = _completed_run(tmp_path)
    rows = _run_rows(harness)
    kept = _differential_from_same_job(job_path, rows)
    legacy = _differential_from_same_job(job_path, rows)
    legacy["provenance"].pop("dataset_hash")
    _append_refs(harness, _persist_differential(harness, kept, "c1"),
                 _persist_differential(harness, legacy, "c2"))

    session_state.session.dataset = None
    summary = dataset_load(pipeline_path=str(harness.pipeline_path))

    assert "comparison_id\tresult_id\tgroup_a\tgroup_b" in summary, summary
    assert f"c1\t{kept['provenance']['result_id']}\tg1\tg2" in summary
    assert (f"c2\t{legacy['provenance']['result_id']}"
            "\tlegacy_result_without_dataset_hash") in summary
    # 指紋を持たない旧 run の結果は、session で作り直す手順を案内する。
    assert "dataset_differential" in summary
    ds = session_state.session.dataset
    assert kept["provenance"]["result_id"] in ds.results
