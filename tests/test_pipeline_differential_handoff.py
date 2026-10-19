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
from metabolomix.pipeline import report, service, session_handoff, store
from tests.test_pipeline_matrix_handoff import _completed_run

_COMPARISON = {"comparison_id": "c1", "reference_group": "g1", "test_group": "g2"}
#: session で後からやり直す、run（conservative-v1 の既定 = pqn）とは別のレシピ。
_SESSION_RECIPE = {"normalize": "none", "impute": "half_min"}


def _run_rows(harness) -> list[dict]:
    return list(harness.result_data("sample_manifest")["rows"])


def _pipeline_preprocess_and_compare(harness, job_path: Path) -> tuple[dict, dict]:
    """run と同じ job から ds を作り、pipeline の前処理 handler と比較関数を回す。

    前処理は `service._handle_preprocess` そのもので永続化する（run が残す
    成果物の形をテスト側で真似しない）。戻り値は (前処理の ref, 差次的結果)。
    """
    rows = _run_rows(harness)
    ds = mztab_loading.load_dataset_state(job_path=job_path)
    context = {"runtime": {"dataset": ds, "metadata": rows},
               "request": {"preprocess": {}}, "request_meta": {"revision": 1},
               "pipeline_root": Path(harness.pipeline_path)}
    preprocess_ref = service._handle_preprocess(context)["result_refs"][0]
    return preprocess_ref, dataset_service.run_comparison(ds, _COMPARISON, rows)


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


def _run_with_differential(tmp_path, *, mutate=None, with_preprocess=True):
    harness, job_path = _completed_run(tmp_path)
    preprocess_ref, result = _pipeline_preprocess_and_compare(harness, job_path)
    if mutate is not None:
        result = mutate(copy.deepcopy(result))
    refs = [preprocess_ref] if with_preprocess else []
    _append_refs(harness, *refs, _persist_differential(harness, result))
    return harness, result


def _exported_preprocess_line(tmp_path, result_id: str) -> str:
    from metabolomix.tools.dataset_analysis_tools import dataset_export_differential

    out = tmp_path / "restored.tsv"
    payload = json.loads(dataset_export_differential(str(out), result_id=result_id))
    assert payload.get("status") == "success", payload
    return next(line for line in out.read_text(encoding="utf-8").splitlines()
                if line.startswith("# preprocess = "))


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
    first_pre, first = _pipeline_preprocess_and_compare(harness, job_path)
    first_ref = _persist_differential(harness, first)
    # export が同じ比較を再計算して、もう一度 `differential:c1` として永続化した。
    second_pre, second = _pipeline_preprocess_and_compare(harness, job_path)
    second_ref = _persist_differential(harness, second)
    _append_refs(harness, first_pre, first_ref, second_pre, second_ref)
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
    preprocess_ref, kept = _pipeline_preprocess_and_compare(harness, job_path)
    _, legacy = _pipeline_preprocess_and_compare(harness, job_path)
    legacy["provenance"].pop("dataset_hash")
    _append_refs(harness, preprocess_ref, _persist_differential(harness, kept, "c1"),
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


# ---------- 書き出しのメタ行 `# preprocess = ...`（run の条件を名乗る） ----------

def _run_recipe(harness) -> dict:
    refs = [r for r in store.load_run(harness.pipeline_path)["results"]
            if r.get("output_name") == "preprocess"]
    return store.read_result_data(
        harness.pipeline_path, refs[-1:], "preprocess")[0]["resolved_recipe"]


def test_exporting_a_restored_result_states_the_runs_recipe(tmp_path):
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    from metabolomix.tools.mztab_tools import dataset_load

    harness, result = _run_with_differential(tmp_path)
    run_recipe = _run_recipe(harness)
    assert run_recipe != _SESSION_RECIPE

    session_state.session.dataset = None
    dataset_load(pipeline_path=str(harness.pipeline_path))
    # session で別のレシピの前処理をやり直しても、引き継いだ結果の条件は run のもの。
    dataset_preprocess(**_SESSION_RECIPE)
    session_recipe = session_state.session.dataset.preprocessing_recipe
    assert session_recipe["normalize"] == "none" != run_recipe["normalize"]

    line = _exported_preprocess_line(tmp_path, result["provenance"]["result_id"])
    assert line == f"# preprocess = {run_recipe}"


def test_exporting_a_restored_result_without_the_runs_preprocess_says_unknown(tmp_path):
    from metabolomix.analysis.dataset_export import HANDOFF_PREPROCESS_UNKNOWN
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    from metabolomix.tools.mztab_tools import dataset_load

    harness, result = _run_with_differential(tmp_path, with_preprocess=False)
    result_id = result["provenance"]["result_id"]

    session_state.session.dataset = None
    dataset_load(pipeline_path=str(harness.pipeline_path))
    assert result_id in session_state.session.dataset.results
    dataset_preprocess(**_SESSION_RECIPE)

    line = _exported_preprocess_line(tmp_path, result_id)
    assert line == f"# preprocess = {HANDOFF_PREPROCESS_UNKNOWN}"
    assert line.startswith("# preprocess = unknown")


def test_the_runs_preprocess_result_is_not_restored_into_session_results(tmp_path):
    harness, result = _run_with_differential(tmp_path)

    ds = session_handoff.load_run_dataset(harness.pipeline_path)

    # 戻すのはレシピだけ。前処理結果そのものを戻すと「今の前処理」に見える。
    assert set(ds.results) == {result["provenance"]["result_id"]}
    assert ds.preprocess_id is None
    assert ds.results[result["provenance"]["result_id"]]["provenance"][
        "handoff"]["preprocess"] == _run_recipe(harness)
