"""必須成果物の確定・決定的な品質レポート・部分完了判定（spec §7.4, §11）。

``required_outputs`` / ``evaluate_target`` / ``write_pipeline_report`` /
``persist_result`` の4関数を、ケースごとに小さくTDDする。fixtureは
``tests.pipeline_fixtures`` の合成dictビルダーを再利用し、実rawや既存成果物は
使わない（CLAUDE.md方針）。``evaluate_target`` の部分完了ケースだけは
``metabolomix.pipeline.store.create_run`` と本モジュールの ``persist_result`` で
実ファイル・実hash・実IDを ``tmp_path`` へ作る——path/hashを持たないダミーの
result_refで有効成果物を代用しない。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline.report import (
    evaluate_target,
    persist_result,
    required_outputs,
    write_pipeline_report,
)
from metabolomix.pipeline.request import resolve_request
from metabolomix.pipeline.store import create_run, load_run


def _minimal_inputs(root: Path) -> dict:
    return {"source_root": str(root), "fingerprint": "f" * 64, "raw_inventory": []}


# ---------- required_outputs ----------

def test_differential_requires_background_tsv_for_each_comparison():
    request = {
        "effective_target": "differential", "save_project": False,
        "comparisons": [{"comparison_id": "treated_vs_control"},
                        {"comparison_id": "recovery_vs_control"}],
    }
    outputs = required_outputs(request)
    assert "pca" in outputs
    assert "tsv:treated_vs_control" in outputs
    assert "tsv:recovery_vs_control" in outputs
    assert "gui_project" not in outputs


def test_exploratory_target_has_no_comparison_outputs():
    request = {"effective_target": "exploratory", "save_project": False, "comparisons": []}
    outputs = required_outputs(request)
    assert outputs == ["preprocess", "pca", "pca_figure", "quality_report"]


def test_save_project_adds_the_gui_project_output():
    request = {"effective_target": "exploratory", "save_project": True, "comparisons": []}
    outputs = required_outputs(request)
    assert "gui_project" in outputs


# ---------- persist_result ----------

def test_persist_result_registers_an_already_written_file_with_a_real_hash(tmp_path):
    import hashlib

    pca_png = tmp_path / "figures" / "pca.png"
    pca_png.parent.mkdir(parents=True)
    pca_png.write_bytes(b"fake-png-bytes")
    expected_hash = hashlib.sha256(b"fake-png-bytes").hexdigest()

    ref = persist_result(tmp_path, {
        "output_name": "pca_figure", "kind": "pca_figure", "result_id": "res_abc",
        "path": pca_png, "parent_ids": ["res_pca"], "request_revision": 1,
    })

    assert ref["output_name"] == "pca_figure"
    assert ref["result_id"] == "res_abc"
    assert ref["kind"] == "pca_figure"
    assert ref["relative_path"] == "figures/pca.png"
    assert ref["hash"] == expected_hash
    assert ref["parent_ids"] == ["res_pca"]
    assert ref["request_revision"] == 1


def test_persist_result_rejects_a_path_outside_the_pipeline_root(tmp_path):
    outside = tmp_path.parent / "outside.png"
    outside.write_bytes(b"x")
    try:
        with pytest.raises(DomainError) as exc:
            persist_result(tmp_path, {"output_name": "pca_figure", "kind": "pca_figure",
                                      "result_id": "res_abc", "path": outside})
        assert exc.value.code == "OUTPUT_OUTSIDE_PIPELINE_ROOT"
    finally:
        outside.unlink(missing_ok=True)


def test_persist_result_rejects_a_missing_file():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.raises(DomainError) as exc:
            persist_result(Path(tmp), {"output_name": "pca_figure", "kind": "pca_figure",
                                       "result_id": "res_abc", "path": Path(tmp) / "none.png"})
        assert exc.value.code == "OUTPUT_ARTIFACT_MISSING"


def test_persist_result_writes_inline_numeric_data_converting_nonfinite_to_null(tmp_path):
    import math

    ref = persist_result(tmp_path, {
        "output_name": "pca", "kind": "pca", "result_id": "res_pca",
        "data": {"explained_variance_ratio": [0.4, math.nan], "n_samples": 8},
    })

    written = json.loads((tmp_path / ref["relative_path"]).read_text(encoding="utf-8"))
    assert written["data"]["explained_variance_ratio"][0] == 0.4
    assert written["data"]["explained_variance_ratio"][1] is None
    assert written["nonfinite_reasons"]
    assert ref["nonfinite_reasons"]


# ---------- 再実行での成果物の上書き（2026-09-17 Stage B の #2）----------

def test_rewriting_a_result_with_new_content_does_not_clobber_the_old_file(tmp_path):
    """同じ result_id で内容が変わったら、別のファイルへ書く。

    `rerun_upstream=true` の再実行は同じ result_id を使う（pipeline_id 由来で
    revision も上がらない）。同じパスへ上書きすると、先に記録された ref の
    (relative_path, hash) が永久に不一致になり、以後の `verify_result_refs` が
    必ず落ちる。実際 Stage B の実データ run はこれで止まった。
    """
    first = persist_result(tmp_path, {
        "output_name": "profile", "kind": "profile", "result_id": "res_profile_x",
        "data": {"profile": {"revision": 1}},
    })
    second = persist_result(tmp_path, {
        "output_name": "profile", "kind": "profile", "result_id": "res_profile_x",
        "data": {"profile": {"revision": 2}},
    })

    assert first["relative_path"] != second["relative_path"]
    assert first["hash"] != second["hash"]
    # 先の ref が指すファイルは、先の hash のまま検証できること。
    import hashlib
    actual = hashlib.sha256((tmp_path / first["relative_path"]).read_bytes()).hexdigest()
    assert actual == first["hash"]


def test_rewriting_a_result_with_identical_content_reuses_the_same_file(tmp_path):
    """内容が同じ再実行はファイルを増やさない（冪等）。"""
    payload = {"output_name": "profile", "kind": "profile",
               "result_id": "res_profile_y", "data": {"profile": {"revision": 1}}}
    first = persist_result(tmp_path, dict(payload))
    second = persist_result(tmp_path, dict(payload))

    assert first["relative_path"] == second["relative_path"]
    assert first["hash"] == second["hash"]
    assert len(list((tmp_path / "results").iterdir())) == 1


# ---------- evaluate_target ----------

def _build_record(tmp_path, *, target, comparisons=None, save_project=False):
    """create_runでpipeline_rootとpipeline-run.v1レコードを1本作る。

    `upstream.verification`は既定で完了検証済み（Task 3の`completion_status`
    語彙で"completed"）にする——`evaluate_target`がこれを検査する
    （レビュー指摘3）ため、ここは「一通り正しく進んだpipeline」を表す既定値で、
    上流が未検証であることを試したいテストだけが明示的に上書きする。
    """
    source_root = tmp_path / "source"
    source_root.mkdir()
    request = resolve_request(source_root, {
        "target": target, "comparisons": comparisons or [], "save_project": save_project})
    inputs = _minimal_inputs(source_root)
    pipeline_root = create_run(source_root, request, inputs)
    record = load_run(pipeline_root)
    record["upstream"]["verification"] = {"status": "completed"}
    return pipeline_root, record


def _write_dummy_file(pipeline_root, relative, content=b"x") -> Path:
    path = Path(pipeline_root) / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def test_evaluate_target_is_needs_input_when_differential_has_no_comparisons(tmp_path):
    pipeline_root, record = _build_record(tmp_path, target="differential", comparisons=[])
    evaluation = evaluate_target(record)
    assert evaluation["status"] == "needs_input"
    assert evaluation["reason_codes"] == ["COMPARISON_REQUIRED"]


def test_evaluate_target_is_completed_when_all_exploratory_outputs_are_achieved(tmp_path):
    pipeline_root, record = _build_record(tmp_path, target="exploratory")
    results = []
    for name in ("preprocess", "pca", "pca_figure", "quality_report"):
        path = _write_dummy_file(pipeline_root, f"results/{name}.bin", name.encode("ascii"))
        results.append(persist_result(pipeline_root, {
            "output_name": name, "kind": name, "result_id": f"res_{name}", "path": path}))
    record["results"] = results

    evaluation = evaluate_target(record)

    assert evaluation["status"] == "completed"
    assert evaluation["reason_codes"] == []
    assert evaluation["missing_outputs"] == []
    assert set(evaluation["achieved_outputs"]) == {"preprocess", "pca", "pca_figure", "quality_report"}


def _differential_pipeline(tmp_path):
    """差次的目標の1本を、実ファイル・実hash・実IDで一通り「達成」させて作る。

    brief「evaluate_targetの部分完了テスト」の共有fixture——TSVを生成した
    正常ケース・EXPORT_BACKGROUND_EMPTYによる部分完了ケース・保存後に改変した
    hash不一致ケースを、すべてここから枝分かれさせる。
    """
    from metabolomix.analysis.dataset_export import export_dataset_result
    from metabolomix.analysis.dataset_service import compare_dataset, pca_dataset, preprocess_dataset
    from tests.pipeline_fixtures import make_dataset

    comparisons = [{"comparison_id": "t_vs_c", "reference_group": "control", "test_group": "treated"}]
    pipeline_root, record = _build_record(tmp_path, target="differential", comparisons=comparisons)

    ds = make_dataset()
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    pca_result = pca_dataset(ds, n_components=2)
    differential_result = compare_dataset(ds, ds.sample_names[:4], ds.sample_names[4:])

    results = []
    results.append(persist_result(pipeline_root, {
        "output_name": "preprocess", "kind": "preprocess", "result_id": ds.preprocess_id,
        "data": {"recipe": ds.preprocessing_recipe}}))
    results.append(persist_result(pipeline_root, {
        "output_name": "pca", "kind": "pca",
        "result_id": pca_result["provenance"]["result_id"],
        "data": {"explained_variance_ratio": pca_result["explained_variance_ratio"]}}))
    pca_png = _write_dummy_file(pipeline_root, "figures/pca.png", b"\x89PNG-pca")
    results.append(persist_result(pipeline_root, {
        "output_name": "pca_figure", "kind": "pca_figure", "result_id": "res_pca_fig",
        "path": pca_png}))
    qr_path = _write_dummy_file(pipeline_root, "report.md", b"# placeholder")
    results.append(persist_result(pipeline_root, {
        "output_name": "quality_report", "kind": "quality_report", "result_id": "res_qr",
        "path": qr_path}))
    results.append(persist_result(pipeline_root, {
        "output_name": "differential:t_vs_c", "kind": "differential",
        "result_id": differential_result["provenance"]["result_id"],
        "data": {"a": differential_result["a"], "b": differential_result["b"]}}))
    volcano_png = _write_dummy_file(pipeline_root, "figures/volcano_t_vs_c.png", b"\x89PNG-volcano")
    results.append(persist_result(pipeline_root, {
        "output_name": "volcano:t_vs_c", "kind": "volcano", "result_id": "res_volcano",
        "path": volcano_png}))

    tsv_path = Path(pipeline_root) / "exports" / "t_vs_c.tsv"
    tsv_path.parent.mkdir(parents=True)
    export_dataset_result(ds, differential_result, tsv_path)
    results.append(persist_result(pipeline_root, {
        "output_name": "tsv:t_vs_c", "kind": "tsv",
        "result_id": differential_result["provenance"]["result_id"] + "_tsv",
        "path": tsv_path}))

    record["results"] = results
    return pipeline_root, record


def test_evaluate_target_is_completed_with_a_real_generated_tsv(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)
    evaluation = evaluate_target(record)
    assert evaluation["status"] == "completed"
    assert evaluation["reason_codes"] == []
    assert "tsv:t_vs_c" in evaluation["achieved_outputs"]


def test_evaluate_target_is_partial_when_the_tsv_is_missing_for_an_empty_background(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)
    record["results"] = [r for r in record["results"] if r["output_name"] != "tsv:t_vs_c"]
    record["output_failures"] = {"tsv:t_vs_c": "EXPORT_BACKGROUND_EMPTY"}

    evaluation = evaluate_target(record)

    assert evaluation["status"] == "partial"
    assert "EXPORT_BACKGROUND_EMPTY" in evaluation["reason_codes"]
    assert "tsv:t_vs_c" in evaluation["missing_outputs"]
    assert "tsv:t_vs_c" not in evaluation["achieved_outputs"]


def test_evaluate_target_treats_a_tampered_output_as_not_achieved(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)
    tsv_ref = next(r for r in record["results"] if r["output_name"] == "tsv:t_vs_c")
    tampered = Path(pipeline_root) / tsv_ref["relative_path"]
    tampered.write_text("tampered contents", encoding="utf-8")

    evaluation = evaluate_target(record)

    assert evaluation["status"] == "partial"
    assert "RESULT_INTEGRITY_MISMATCH" in evaluation["reason_codes"]
    assert "tsv:t_vs_c" in evaluation["missing_outputs"]
    assert "tsv:t_vs_c" not in evaluation["achieved_outputs"]


@pytest.mark.parametrize("verification", [None, {"status": "failed"}])
def test_evaluate_target_does_not_achieve_outputs_when_upstream_is_not_verified(tmp_path, verification):
    """レビュー指摘3: brief「evaluate_targetは上流verifiedと全必須outputのhash/ID
    を検査し」にもかかわらず、`record["upstream"]["verification"]`が一切
    読まれておらず、上流が未解決(None)・failedでも成果物refのhash/IDさえ
    有効なら達成扱いになっていた。Task 3が確立した`completion_status`語彙
    ("completed"/"partial"/"failed")を再利用し、"completed"以外
    （未解決のNoneも含む）は上流未検証として扱う。"""
    pipeline_root, record = _differential_pipeline(tmp_path)
    record["upstream"]["verification"] = verification

    evaluation = evaluate_target(record)

    assert evaluation["status"] != "completed"
    assert "UPSTREAM_NOT_VERIFIED" in evaluation["reason_codes"]
    assert "tsv:t_vs_c" not in evaluation["achieved_outputs"]


# ---------- write_pipeline_report ----------

_SECTION_HEADERS_IN_ORDER = (
    "## Source", "## Method", "## Execution Receipt",
    "## Role / Group / Batch / Order Provenance", "## Applied / Skipped",
    "## PCA", "## Comparisons", "## InChIKey Coverage", "## Unverified Conditions",
)


def test_report_sections_appear_in_the_fixed_order(tmp_path):
    pipeline_root, record = _build_record(tmp_path, target="exploratory")

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    positions = [text.index(header) for header in _SECTION_HEADERS_IN_ORDER]
    assert positions == sorted(positions)


def test_report_escapes_newlines_and_pipes_in_table_cells(tmp_path):
    pipeline_root, record = _build_record(tmp_path, target="exploratory")
    record["inputs"]["manifest"] = [{
        "sample_id": "s1|weird\nname", "source_file": "S1.wiff", "role": "sample",
        "group": "control", "batch": "B1", "injection_order": 1, "qc_pool": None,
        "include": True,
    }]

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    assert "s1|weird\nname" not in text
    assert "s1\\|weird<br>name" in text


def test_report_makes_no_unverified_biological_interpretation(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    forbidden = ("生物学的", "biologically significant", "示唆", "統計的に有意差がある")
    for phrase in forbidden:
        assert phrase not in text


def test_report_does_not_leak_an_older_revisions_artifact(tmp_path):
    """append-onlyのresultsに古いrevisionの成果物が残っていても、レポートは
    最後に追記された（現行revisionの）成果物だけを指す。"""
    pipeline_root, record = _build_record(tmp_path, target="exploratory")
    stale_pca = _write_dummy_file(pipeline_root, "results/pca_old.bin", b"stale")
    fresh_pca = _write_dummy_file(pipeline_root, "results/pca_new.bin", b"fresh")
    stale_ref = persist_result(pipeline_root, {
        "output_name": "pca", "kind": "pca", "result_id": "res_old", "path": stale_pca,
        "request_revision": 1})
    fresh_ref = persist_result(pipeline_root, {
        "output_name": "pca", "kind": "pca", "result_id": "res_new", "path": fresh_pca,
        "request_revision": 2})
    record["results"] = [stale_ref, fresh_ref]

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    assert "res_new" in text
    assert "res_old" not in text


def test_report_marks_unadjusted_confounded_comparisons(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)
    differential_ref = next(r for r in record["results"] if r["output_name"] == "differential:t_vs_c")
    payload = json.loads((pipeline_root / differential_ref["relative_path"]).read_text(encoding="utf-8"))
    payload["data"]["provenance"] = {"comparison": {"unadjusted_confounded": True}}
    (pipeline_root / differential_ref["relative_path"]).write_text(
        json.dumps(payload), encoding="utf-8")
    # ハッシュを付け替えないと改変検知に引っかかって「達成」扱いされない。
    import hashlib
    differential_ref["hash"] = hashlib.sha256(
        (pipeline_root / differential_ref["relative_path"]).read_bytes()).hexdigest()

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    assert "unadjusted" in text.lower()


def test_report_states_real_inchikey_coverage_from_the_generated_tsv(tmp_path):
    """レビュー指摘1: `tsv_summary:<cid>`は`required_outputs`が作らず誰も永続化
    しない存在しない出力名なので、このセクションは常に(n/a)だった（死んだ
    コード）。`tsv:<cid>`として登録済みのTSV自身のメタ行
    （`export_contract.build_meta`が書く`# n_features_total = ...`行）から
    読む——別の出力名を新設しない。make_dataset()は6特徴・全件InChIKey付き
    なので6/6/0のはず。"""
    pipeline_root, record = _differential_pipeline(tmp_path)

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    coverage_section = text.split("## InChIKey Coverage", 1)[1].split(
        "## Unverified Conditions", 1)[0]
    assert "(n/a)" not in coverage_section
    assert "| t_vs_c | 6 | 6 | 0 |" in coverage_section


def test_report_comparisons_table_does_not_say_achieved_when_the_result_is_tampered(tmp_path):
    """レビュー指摘2: `_section_comparisons`はJSONとして読めるかどうかだけで
    achieved/missingを決めており、`evaluate_target`が同じ`record`から出す
    hash整合性検証を無視していた。改変後も構文的にはJSONとして読める内容
    （フィールド値だけ書き換え、hashは付け替えない）で、ヘッダの
    status_at_report_timeと矛盾しないことを見る——これが本来落ちるべきgapで、
    `test_evaluate_target_treats_a_tampered_output_as_not_achieved`は
    `evaluate_target`の戻り値しか見ておらずレポート本文の矛盾を検知できない。"""
    pipeline_root, record = _differential_pipeline(tmp_path)
    differential_ref = next(r for r in record["results"] if r["output_name"] == "differential:t_vs_c")
    tampered_path = Path(pipeline_root) / differential_ref["relative_path"]
    payload = json.loads(tampered_path.read_text(encoding="utf-8"))
    payload["data"]["a"] = "tampered-group-name"
    tampered_path.write_text(json.dumps(payload), encoding="utf-8")
    # hashは付け替えない——ハッシュ不一致こそが今回検知したい整合性違反。

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    assert "status_at_report_time: completed" not in text
    assert "RESULT_INTEGRITY_MISMATCH" in text
    comparisons_section = text.split("## Comparisons", 1)[1].split(
        "## InChIKey Coverage", 1)[0]
    assert "achieved" not in comparisons_section


def test_write_pipeline_report_returns_the_status_snapshot(tmp_path):
    pipeline_root, record = _build_record(tmp_path, target="exploratory")

    summary = write_pipeline_report(record, pipeline_root / "report.md")

    assert summary["path"] == str(pipeline_root / "report.md")
    assert summary["status"] in {"completed", "partial", "failed", "needs_input"}
    assert "quality_report" in summary["missing_outputs"]


def test_report_links_stay_relative_and_inside_the_pipeline(tmp_path):
    pipeline_root, record = _differential_pipeline(tmp_path)

    write_pipeline_report(record, pipeline_root / "report.md")

    text = (pipeline_root / "report.md").read_text(encoding="utf-8")
    links = re.findall(r"\]\(([^)]+)\)", text)
    assert links  # PCAセクションで少なくとも1件は張られる
    for link in links:
        assert ".." not in link
        assert not Path(link).is_absolute()


def test_method_section_shows_the_msp_file_but_not_its_location():
    from metabolomix.pipeline.report import _section_method
    record = {"inputs": {
        "library_mode": "msp_only",
        "method": {"source_path": "m.txt", "sha256": "a" * 64},
        "lbm": {"path": None, "sha256": None},
        "msp": {"path": "//nas/secret_share/lab_neg.msp", "file": "lab_neg.msp",
                "sha256": "b" * 64, "source": "config_file"},
        "exe": {"version": None},
    }}
    text = _section_method(record)
    assert "- library_mode: msp_only" in text
    assert "- msp_file: lab_neg.msp" in text
    assert f"- msp_sha256: {'b' * 64}" in text
    assert "- msp_source: config_file" in text
    assert "secret_share" not in text
