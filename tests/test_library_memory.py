"""データセットが同定に使ったライブラリを覚える（Task 6）。

ライブラリのファイル名は作り物（研究室の実ファイルに依存しない・書かない）。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.library_memory import (
    is_lbm_annotator,
    public_view,
    remember_from_folder,
    remember_from_job,
    remember_from_run,
    unknown_record,
    verify_record,
)
from metabolomix.handoff.schema import SCHEMA_VERSION, AnalysisJob
from tests.dbs_fixture import record, write_dbs
from tests.msdial_folder_fixtures import write_msdial_folder

_SAMPLES = [{"file_id": i, "name": f"s{i}", "class_id": "A" if i < 2 else "B"} for i in range(3)]


def _folder(tmp_path):
    write_msdial_folder(tmp_path, spots=[{"name": "PC 34:1"}], samples=_SAMPLES)
    return tmp_path


def _dbs(path):
    return write_dbs(path, [("lib", [record()], ["k"])])


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------- 解析フォルダ ----------

def test_folder_with_one_loaded_dbs_is_known(tmp_path):
    _folder(tmp_path)
    _dbs(tmp_path / "Proj_Loaded.msp2.dbs")
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "known" and rec["kind"] == "dbs" and rec["origin"] == "loaded_dbs"
    assert rec["file"] == "Proj_Loaded.msp2.dbs"
    assert rec["sha256"] == _sha(tmp_path / "Proj_Loaded.msp2.dbs")
    assert rec["lbm_also"] is False
    assert "path" in rec and "path" not in public_view(rec)


def test_folder_with_two_loaded_dbs_is_unknown_with_candidates(tmp_path):
    _folder(tmp_path)
    _dbs(tmp_path / "Proj_Loaded.msp2.dbs")
    _dbs(tmp_path / "Proj_Loaded.text2.dbs")
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "unknown"
    assert rec["candidates"] == ["Proj_Loaded.msp2.dbs", "Proj_Loaded.text2.dbs"]
    assert rec["file"] is None and rec["kind"] is None


def test_folder_with_only_lbm_dbs_is_lbm_only(tmp_path):
    _folder(tmp_path)
    _dbs(tmp_path / "Proj_Loaded.lbm2.dbs")
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "lbm_only" and rec["kind"] == "lbm"


def test_folder_with_one_msp_dbs_and_lbm_dbs_is_known_with_lbm_also(tmp_path):
    _folder(tmp_path)
    _dbs(tmp_path / "Proj_Loaded.msp2.dbs")
    _dbs(tmp_path / "Proj_Loaded.lbm2.dbs")
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "known" and rec["kind"] == "dbs"
    assert rec["file"] == "Proj_Loaded.msp2.dbs" and rec["lbm_also"] is True


def test_folder_without_loaded_dbs_is_unknown(tmp_path):
    _folder(tmp_path)
    rec = remember_from_folder(select_alignment(tmp_path))
    assert rec["status"] == "unknown" and rec["candidates"] == []


# ---------- ジョブ ----------

def _job(tmp_path, method_text: str) -> Path:
    method = tmp_path / "method.txt"
    method.write_text(method_text, encoding="ascii")
    job = AnalysisJob(
        schema=SCHEMA_VERSION, job_id="j1", status="completed", created_at="t", updated_at="t",
        dataset_root=str(tmp_path), input_count=1, software_name="MS-DIAL",
        software_version="5", execution_mode="console", method_file=str(method),
        omics="lipidomics", polarity="positive", measure="peak_height",
        run_dir=str(tmp_path / "run"))
    path = tmp_path / "analysis-job.json"
    job.save(path)
    return path


def test_job_records_effective_method_msp(tmp_path):
    msp = tmp_path / "lab.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    rec = remember_from_job(_job(tmp_path, f"Msp file path: {msp}\n"))
    assert rec["status"] == "known" and rec["kind"] == "msp"
    assert rec["origin"] == "effective_method" and rec["file"] == "lab.msp"
    assert rec["sha256"] == _sha(msp)
    assert "path" not in public_view(rec)


def test_job_with_only_lbm_declared_is_lbm_only(tmp_path):
    lbm = tmp_path / "lab.lbm2"
    lbm.write_text("x", encoding="ascii")
    rec = remember_from_job(_job(tmp_path, f"Msp file path: \nLbm file path: {lbm}\n"))
    assert rec["status"] == "lbm_only" and rec["kind"] == "lbm" and rec["file"] == "lab.lbm2"


def test_job_without_library_is_unknown(tmp_path):
    assert remember_from_job(_job(tmp_path, "Ion mode: Positive\n"))["status"] == "unknown"


# ---------- パイプラインの run ----------

def _run(tmp_path, *, msp, library_mode, lbm=None) -> Path:
    inputs = {"library_mode": library_mode, "msp": msp}
    if lbm is not None:
        inputs["lbm"] = lbm
    root = tmp_path / "pipeline_x"
    root.mkdir()
    (root / "pipeline-run.json").write_text(
        json.dumps({"schema": "pipeline-run.v1", "inputs": inputs}), encoding="utf-8")
    return root


def _msp_input(path: Path, sha: str | None = None) -> dict:
    return {"path": str(path), "file": path.name,
            "sha256": sha or _sha(path), "source": "config_file"}


def test_run_msp_only_is_known_from_run_inputs(tmp_path):
    msp = tmp_path / "lab.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    rec = remember_from_run(_run(tmp_path, msp=_msp_input(msp), library_mode="msp_only"))
    assert rec["status"] == "known" and rec["kind"] == "msp" and rec["origin"] == "run_inputs"
    assert rec["file"] == "lab.msp" and rec["sha256"] == _sha(msp) and rec["lbm_also"] is False


def test_run_pinned_sha_mismatch(tmp_path):
    msp = tmp_path / "lab.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    root = _run(tmp_path, msp=_msp_input(msp, sha="0" * 64), library_mode="msp_only")
    rec = remember_from_run(root)
    assert rec["status"] == "mismatch" and rec["file"] == "lab.msp"


def test_run_with_missing_msp_file_is_mismatch(tmp_path):
    msp = tmp_path / "gone.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    root = _run(tmp_path, msp=_msp_input(msp), library_mode="msp_only")
    msp.unlink()
    assert remember_from_run(root)["status"] == "mismatch"


def test_msp_then_lbm_records_msp_with_lbm_also(tmp_path):
    msp = tmp_path / "lab.msp"
    msp.write_text("NAME: x\n", encoding="ascii")
    lbm = {"path": str(tmp_path / "lab.lbm2"), "file": "lab.lbm2", "sha256": "1" * 64,
           "source": "config_file"}
    rec = remember_from_run(_run(tmp_path, msp=_msp_input(msp), library_mode="msp_then_lbm",
                                 lbm=lbm))
    assert rec["status"] == "known" and rec["kind"] == "msp" and rec["lbm_also"] is True


def test_auto_run_with_lbm_only_is_lbm_only(tmp_path):
    lbm = tmp_path / "lab.lbm2"
    lbm.write_text("x", encoding="ascii")
    rec = remember_from_run(_run(
        tmp_path, library_mode="auto",
        msp={"path": None, "file": None, "sha256": None, "source": "not_used"},
        lbm={"path": str(lbm), "sha256": _sha(lbm)}))
    assert rec["status"] == "lbm_only" and rec["kind"] == "lbm" and rec["file"] == "lab.lbm2"


def test_auto_run_without_any_library_is_unknown(tmp_path):
    rec = remember_from_run(_run(
        tmp_path, library_mode="auto",
        msp={"path": None, "file": None, "sha256": None, "source": "not_used"}))
    assert rec["status"] == "unknown"


# ---------- 共通 ----------

def test_mztab_direct_is_unknown():
    rec = unknown_record()
    assert rec["status"] == "unknown"
    assert rec["kind"] is None and rec["file"] is None and rec["candidates"] == []


def test_verify_record_remeasures_and_flags_a_changed_file(tmp_path):
    _folder(tmp_path)
    dbs = _dbs(tmp_path / "Proj_Loaded.msp2.dbs")
    rec = remember_from_folder(select_alignment(tmp_path))
    assert verify_record(rec)["status"] == "known"
    dbs.write_bytes(dbs.read_bytes() + b"x")
    assert verify_record(rec)["status"] == "mismatch"
    assert rec["status"] == "known"  # 元の記録は書き換えない
    dbs.unlink()
    assert verify_record(rec)["status"] == "mismatch"


def test_public_view_passes_none_and_hides_only_path():
    assert public_view(None) is None
    view = public_view({"status": "known", "path": "/x/y.msp", "file": "y.msp"})
    assert view == {"status": "known", "file": "y.msp"}


def test_is_lbm_annotator_follows_the_console_key_rule():
    # Console は LBM の注釈器キーに LBM のファイルのパス（拡張子 .lbm + 数字）を使う。
    assert is_lbm_annotator("D:/x/some.lbm2") is True
    assert is_lbm_annotator("\\\\host\\share\\some.LBM") is True
    # MSP の注釈器キーは利用者が付けた名前。名前だけでは LBM でないとは言えない（GUI は自由名）。
    assert is_lbm_annotator("LabMsp") is None
    assert is_lbm_annotator(None) is None
    assert is_lbm_annotator("") is None
