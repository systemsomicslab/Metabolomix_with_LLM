"""v2 の一気通貫（plan Task 14・Review D）。

合成入力だが、通るのは本物の経路: 実ファイルの mzTab-M と `.arf` を
`parse_mztab` / `build_dataset_state` / `arf_reader.deserialize` で読み、
`run_engine` が v2 handler を順に呼び、`pipeline-run.json` に保存された成果物を
**ファイルから読み戻して**検証する（モックの成功JSONは見ない）。

ここで縛る不変条件:

1. binding を直して再開しても **Console を起こし直さない**（spec §6.2）。
2. QC が fail でも、統計が計算できるなら続行して QC 不合格を付記する（spec §11）。
3. 全 feature が未同定でも出力は完成する（注釈の無い行を捨てない）。
4. 統計の追加・削除は該当 stage だけを動かす。
5. 必須成果物が揃って初めて `completed`。
"""
from __future__ import annotations

import json

import pytest

from tests.metabolomics_fixtures import DEFAULT_STATISTICS, MetabolomicsHarness


# ---------- brief記載のRED ----------

def test_new_binding_resumes_without_second_console(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    first = h.run(binding_mode="ambiguous", annotated=False)
    assert first["status"] == "needs_input"
    h.resume({"feature_bindings": first["allowed_binding_update"]})
    assert h.launch_count == 1


# ---------- 正常系 ----------

def test_a_complete_run_produces_every_required_output(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    result = h.run()
    assert result["status"] == "completed"
    produced = set(h.output_names())
    for required in ("profile", "execution_manifest", "sample_manifest",
                     "assay_evidence", "feature_bindings", "matrix",
                     "qc_population", "qc", "feature_table", "quality_report",
                     "statistic:anova", "statistic:welch"):
        assert required in produced, required


def test_the_evidence_table_is_per_injection(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    evidence = h.result_data("assay_evidence")
    assert evidence["availability"] is True
    # 15注入 × 4 feature（代表値の複製なら注入数ぶんの行は出ない）。
    assert evidence["coverage"]["rows"] == 15 * 4
    assert evidence["rt_unit"] == {"source": "minute", "normalized": "minute"}


def test_bindings_resolve_both_targets_and_build_the_standard_map(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    bindings = h.result_data("feature_bindings")
    assert bindings["status"] == "resolved"
    assert bindings["bindings"]["gaba"]["selected_feature_id"] == "1"
    assert bindings["bindings"]["gaba_d6"]["selected_feature_id"] == "2"
    assert bindings["internal_standard_map"]["pairs"] == [
        {"target_id": "gaba", "target_feature_id": "1",
         "standard_target_id": "gaba_d6", "standard_feature_id": "2"}]


def test_two_recipes_give_two_persisted_matrices(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    record = h.record()
    matrix_ids = {r["result_id"] for r in record["results"]
                  if r["output_name"] == "matrix"}
    # default / ratio × (preprocessed, finalized) = 4本。どれも別ID。
    assert len(matrix_ids) == 4


def test_the_anova_matches_the_hand_computed_example(tmp_path):
    """合成値 [1,2,3]/[2,3,4]/[4,5,6]（×100）は F=7・df=(2,6)。"""
    h = MetabolomicsHarness(tmp_path)
    h.run()
    result = h.result_data("statistic:anova")
    target = next(f for f in result["features"] if f["feature_id"] == "1")
    assert target["f_statistic"] == pytest.approx(7.0)
    assert (target["df_between"], target["df_within"]) == (2, 6)
    assert len(target["tukey"]) == 3


def test_the_welch_statistic_uses_the_internal_standard_ratio_matrix(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    result = h.result_data("statistic:welch")
    assert result["effect_size_definition"] == "log2_arithmetic_mean_ratio"
    matrices = [json.loads("{}")]           # 明示: 行列IDは結果自身が持つ
    assert result["matrix_id"]
    assert result["matrix_recipe_id"] == "ratio"


def test_the_feature_table_keeps_every_feature_and_both_units(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    text = h.result_data("feature_table")
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    header, body = lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]
    feature_ids = {row[header.index("feature_id")] for row in body}
    assert feature_ids == {"1", "2", "3", "4"}
    units = {row[header.index("unit")] for row in body}
    assert {"peak_height", "internal_standard_ratio"} <= units


# ---------- binding が決まらない ----------

def test_ambiguous_binding_stops_before_any_statistic(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    result = h.run(binding_mode="ambiguous")
    assert result["status"] == "needs_input"
    assert h.stage_status("resolve_feature_bindings") == "needs_input"
    assert h.stage_status("qc_raw") == "pending"
    assert "statistic:anova" not in h.output_names()


def test_the_binding_correction_completes_the_run(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    first = h.run(binding_mode="ambiguous")
    assert first["status"] == "needs_input"

    second = h.resume({"feature_bindings": first["allowed_binding_update"]})
    assert second["status"] == "completed"
    bindings = h.result_data("feature_bindings")
    assert bindings["bindings"]["gaba"]["selection"] == "manual"
    assert bindings["bindings"]["gaba"]["selection_reason"]


def test_a_binding_payload_from_another_dataset_is_refused(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    first = h.run(binding_mode="ambiguous")
    payload = dict(first["allowed_binding_update"])
    payload["dataset_hash"] = "f" * 64
    h.resume({"feature_bindings": payload})
    # run 全体の状態ではなく、止まった stage の理由を見る（「どこで・なぜ」が
    # 分からない停止は、利用者が次の一手を選べない）。
    stage = h.record()["stages"]["resolve_feature_bindings"]
    assert stage["status"] == "failed"
    assert stage["error"]["code"] == "FEATURE_BINDING_OVERRIDE_INVALID"
    assert "statistic:anova" not in h.output_names()


# ---------- QC fail でも統計は続行 ----------

def test_a_failing_qc_does_not_stop_the_statistics(tmp_path):
    """QC が落ちても統計は診断用に計算し、結果へ不合格を付記する（spec §11）。"""
    h = MetabolomicsHarness(tmp_path)
    result = h.run(qc_fail=True)
    assert result["status"] == "completed"

    qc = h.result_data("qc")
    rsd = qc["metrics"]["pooled_qc_rsd"]
    # 対象 feature の QC 値は [1,10,100]。RSD は閾値30%を大きく超える。
    failing = [e for e in rsd["elements"] if e["status"] == "fail"]
    assert failing, "QC値を荒らしてもfail要素が出ていない"
    assert h.result_data("statistic:anova")["status"] == "completed"


def test_the_report_separates_qc_from_analysis(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run(qc_fail=True)
    report = h.result_data("quality_report")
    # 3軸が別々に出ること。1つの成否へ畳まない。
    assert "qc_status:" in report
    assert "analysis_status: ready" in report


# ---------- 未同定・分母欠損 ----------

def test_every_feature_unannotated_still_completes(tmp_path):
    """全 feature が未同定でも、対象は m/z と RT（mass_rt）で対応付く。"""
    h = MetabolomicsHarness(tmp_path)
    result = h.run(annotated=False)
    assert result["status"] == "completed"
    text = h.result_data("feature_table")
    assert text.count("\n") > 1


def test_a_missing_standard_denominator_is_locked_not_imputed(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run(standard_missing=True)
    text = h.result_data("feature_table")
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    unit = header.index("unit")
    value = header.index("value")
    imputed = header.index("imputed")
    reason = header.index("exclusion_reason")
    locked = [r for r in rows if r[unit] == "internal_standard_ratio"
              and r[value] == ""]
    assert locked, "分母が使えないセルが欠損として残っていない"
    assert all(r[imputed] == "false" for r in locked)
    assert all("standard_denominator_invalid" in r[reason] for r in locked)


# ---------- 統計の追加・削除 ----------

def test_removing_a_statistic_drops_its_requirement(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run()
    assert "statistic:welch" in h.output_names()

    result = h.resume({"statistics": [DEFAULT_STATISTICS[0]]})
    assert result["status"] == "completed"
    record = h.record()
    assert "statistics:welch" not in [
        s for s, v in record["stages"].items() if v["status"] == "succeeded"]


def test_adding_a_statistic_does_not_redo_the_matrices(tmp_path):
    h = MetabolomicsHarness(tmp_path)
    h.run(statistics=[DEFAULT_STATISTICS[0]])
    before = {r["result_id"] for r in h.record()["results"]
              if r["output_name"] == "matrix"}

    result = h.resume({"statistics": DEFAULT_STATISTICS})
    assert result["status"] == "completed"
    assert "statistic:welch" in h.output_names()
    # 行列は内容hash由来のIDなので、同じ入力なら同じID（作り直しても増えない）。
    after = {r["result_id"] for r in h.record()["results"]
             if r["output_name"] == "matrix"}
    assert before <= after
    assert h.launch_count == 1


# ---------- 記録の健全性 ----------

def test_results_are_append_only_and_verifiable(tmp_path):
    from metabolomix.pipeline import store

    h = MetabolomicsHarness(tmp_path)
    h.run()
    record = h.record()
    failures = store.verify_result_refs(h.pipeline_path, record["results"])
    assert failures == []


def test_a_tampered_artifact_is_detected(tmp_path):
    from pathlib import Path

    from metabolomix.pipeline import store

    h = MetabolomicsHarness(tmp_path)
    h.run()
    record = h.record()
    ref = next(r for r in record["results"] if r["output_name"] == "qc")
    target = Path(record["identity"]["pipeline_root"]) / ref["relative_path"]
    target.write_text("{}", encoding="utf-8")
    assert store.verify_result_refs(h.pipeline_path, [ref])
