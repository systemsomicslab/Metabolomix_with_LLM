"""受領（`metabolomix.pipeline.receipt`）と成果物の登録（`persist_result`）のテスト。

受領は run の最後の工程が書く `results/receipt.json`（`lipidmix.pipeline-receipt.v1`）。
極性ごとの解析フォルダの場所・照合した `.arf`・使ったメソッドとライブラリ（研究室の資産は
ファイル名と sha256 だけ）・試料の対応表・検証の結果と、次の呼び出しを載せる。
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline.engine import run_engine
from metabolomix.pipeline.receipt import build_receipt, persist_result, read_receipt
from metabolomix.pipeline.store import load_run, save_run
from tests.pipeline_fixtures import create_synthetic_run
from tests.pipeline_worker_harness import synthetic_handlers


# ---------- persist_result ----------

def test_persist_result_registers_an_already_written_file_with_a_real_hash(tmp_path):
    png = tmp_path / "figures" / "pca.png"
    png.parent.mkdir(parents=True)
    png.write_bytes(b"fake-png-bytes")

    ref = persist_result(tmp_path, {
        "output_name": "pca_figure", "kind": "pca_figure", "result_id": "res_abc",
        "path": png, "parent_ids": ["res_pca"], "request_revision": 1,
    })

    assert ref["output_name"] == "pca_figure"
    assert ref["result_id"] == "res_abc"
    assert ref["relative_path"] == "figures/pca.png"
    assert ref["hash"] == hashlib.sha256(b"fake-png-bytes").hexdigest()
    assert ref["parent_ids"] == ["res_pca"]
    assert ref["request_revision"] == 1


def test_persist_result_rejects_a_path_outside_the_pipeline_root(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "outside.png"
    outside.write_bytes(b"x")
    with pytest.raises(DomainError) as exc:
        persist_result(root, {"output_name": "x", "kind": "x", "result_id": "r",
                              "path": outside})
    assert exc.value.code == "OUTPUT_OUTSIDE_PIPELINE_ROOT"


def test_persist_result_rejects_a_missing_file(tmp_path):
    with pytest.raises(DomainError) as exc:
        persist_result(tmp_path, {"output_name": "x", "kind": "x", "result_id": "r",
                                  "path": tmp_path / "none.png"})
    assert exc.value.code == "OUTPUT_ARTIFACT_MISSING"


def test_persist_result_writes_inline_numeric_data_converting_nonfinite_to_null(tmp_path):
    ref = persist_result(tmp_path, {
        "output_name": "x", "kind": "x", "result_id": "res_x",
        "data": {"values": [0.4, math.nan], "n": 8},
    })
    written = json.loads((tmp_path / ref["relative_path"]).read_text(encoding="utf-8"))
    assert written["data"]["values"] == [0.4, None]
    assert written["nonfinite_reasons"]
    assert ref["nonfinite_reasons"]


def test_rewriting_a_result_with_new_content_does_not_clobber_the_old_file(tmp_path):
    """同じ result_id で内容が変わったら別のファイルへ書く（先の ref の hash を壊さない）。"""
    first = persist_result(tmp_path, {"output_name": "x", "kind": "x", "result_id": "res_y",
                                      "data": {"revision": 1}})
    second = persist_result(tmp_path, {"output_name": "x", "kind": "x", "result_id": "res_y",
                                       "data": {"revision": 2}})
    assert first["relative_path"] != second["relative_path"]
    actual = hashlib.sha256((tmp_path / first["relative_path"]).read_bytes()).hexdigest()
    assert actual == first["hash"]


def test_rewriting_a_result_with_identical_content_reuses_the_same_file(tmp_path):
    payload = {"output_name": "x", "kind": "x", "result_id": "res_z", "data": {"revision": 1}}
    first = persist_result(tmp_path, dict(payload))
    second = persist_result(tmp_path, dict(payload))
    assert first["relative_path"] == second["relative_path"]
    assert len(list((tmp_path / "results").iterdir())) == 1


# ---------- 受領 ----------

def test_receipt_of_a_two_polarity_run(tmp_path):
    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    assert run_engine(path, synthetic_handlers())["status"] == "completed"

    receipt = json.loads((path / receipt_mod.RECEIPT_RELATIVE_PATH).read_text(encoding="utf-8"))
    assert receipt["schema"] == "lipidmix.pipeline-receipt.v1"
    assert receipt["request_revision"] == 1
    assert list(receipt["units"]) == ["negative", "positive"]
    negative = receipt["units"]["negative"]
    assert negative["status"] == "validated"
    assert negative["analysis_dir"] == "negative/input"
    assert negative["alignment"] == {"arf": "negative_PeakProperties.arf",
                                     "arf2": "negative.arf2"}
    assert negative["validation"]["job_status"] == "completed"
    assert receipt["next"] == ('load_dataset(pipeline_path="<pipeline_path>", '
                               'polarity="negative" または "positive")')
    assert read_receipt(path, load_run(path)) == receipt


def test_receipt_lists_a_failed_unit_with_its_error_and_leaves_it_out_of_next(tmp_path):
    path = create_synthetic_run(tmp_path, polarities=("negative", "positive"))
    result = run_engine(path, synthetic_handlers(fail_stages={"upstream:negative"}))
    assert result["status"] == "partial"

    receipt = read_receipt(path, load_run(path))
    assert receipt["units"]["negative"]["status"] == "failed"
    assert receipt["units"]["negative"]["error"]["code"] == "SYNTHETIC_FAILURE"
    assert receipt["units"]["negative"]["error"]["stage_id"] == "upstream:negative"
    assert '"negative"' not in receipt["next"]
    assert 'polarity="positive"' in receipt["next"]


def test_a_tampered_receipt_is_not_read(tmp_path):
    path = create_synthetic_run(tmp_path)
    run_engine(path, synthetic_handlers())
    (path / receipt_mod.RECEIPT_RELATIVE_PATH).write_text("{}", encoding="utf-8")
    assert read_receipt(path, load_run(path)) is None


def test_receipt_names_library_files_only(tmp_path):
    """研究室の MSP と `[library] lbm` の LBM は、受領にファイル名と sha256 だけを書く。"""
    path = create_synthetic_run(tmp_path)
    record = load_run(path)
    record["units"]["negative"]["inputs"].update({
        "library_mode": "msp_then_lbm",
        "method": {"source_path": "C:/elsewhere/secret_share/method.txt", "sha256": "a" * 64},
        "lbm": {"path": "//nas/secret_share/lab.lbm2", "sha256": "c" * 64,
                "source": "config_file"},
        "msp": {"path": "//nas/secret_share/lab_neg.msp", "sha256": "b" * 64,
                "source": "config_file"},
        "exe": {"path": "C:/elsewhere/secret_share/MsdialConsoleApp.exe", "sha256": "d" * 64},
    })
    save_run(path, record, expected_revision=record["state_revision"])

    receipt = build_receipt(load_run(path), path)
    library = receipt["units"]["negative"]["library"]
    assert library == {
        "mode": "msp_then_lbm",
        "msp": {"file": "lab_neg.msp", "sha256": "b" * 64, "source": "config_file"},
        "lbm": {"file": "lab.lbm2", "sha256": "c" * 64, "source": "config_file"},
    }
    assert receipt["units"]["negative"]["method"]["file"] == "method.txt"
    assert receipt["units"]["negative"]["console"]["exe_file"] == "MsdialConsoleApp.exe"
    assert "secret_share" not in json.dumps(receipt)


def test_receipt_records_the_sample_manifest_it_read(tmp_path):
    path = create_synthetic_run(tmp_path)
    source = Path(load_run(path)["identity"]["source_root"])
    manifest = source / "sample-manifest.tsv"
    manifest.write_text("# schema = sample-manifest.v1\nsample_id\tgroup\nS0\tcontrol\n",
                        encoding="utf-8", newline="\n")
    run_engine(path, synthetic_handlers())

    samples = read_receipt(path, load_run(path))["units"]["negative"]["samples"]
    assert samples["manifest"] == {
        "path": "sample-manifest.tsv",
        "sha256": hashlib.sha256(manifest.read_bytes()).hexdigest()}


def test_a_manifest_that_vanished_is_a_receipt_warning_not_a_crash(tmp_path):
    path = create_synthetic_run(tmp_path, explicit={"sample_manifest": "gone.tsv"})
    assert run_engine(path, synthetic_handlers())["status"] == "completed"
    receipt = read_receipt(path, load_run(path))
    assert [w["code"] for w in receipt["warnings"]] == ["SAMPLE_MANIFEST_INVALID"]
