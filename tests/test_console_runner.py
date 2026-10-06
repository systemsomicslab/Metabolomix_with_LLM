# tests/test_console_runner.py
#
# Console の実行は `metabolomix.console.execution.supervise` が唯一の経路で、
# `console_run` は監視ワーカー（metabolomix.console.worker）を通してそこへ入る。
# ここでの実行系テストは **subprocess をモックしない** —— 監視・停止・収集の
# どれが壊れても緑のままになるため、実際に子プロセスを起こし、その終了コードと
# 生成物で判定する（`_fake_console`）。
import json
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from metabolomix.console.runner import (
    MsdialExeNotFoundError,
    build_msdial_cmd,
    get_exe_path,
)
from metabolomix.console.output_collector import (
    collect_artifacts,
    snapshot,
    _assign_role,
    _infer_mztab_meta,
)
from metabolomix.console.job_manager import (
    create_job,
    count_raw_inputs,
    list_jobs,
    load_job,
    update_status,
)


# ---------- runner ----------

def test_get_exe_path_missing(monkeypatch):
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    with pytest.raises(MsdialExeNotFoundError):
        get_exe_path()


def test_get_exe_path_set(monkeypatch):
    monkeypatch.setenv("MSDIAL_EXE", "C:/MsDial/MsdialConsoleApp.exe")
    assert get_exe_path() == "C:/MsDial/MsdialConsoleApp.exe"


def test_is_console_exe_accepts_output_with_lcms():
    from metabolomix.console.runner import is_console_exe
    completed = MagicMock()
    completed.stdout = "MSDIAL Console Application 5.5\n  lcms   Run LC-MS data processing\n"
    with patch("subprocess.run", return_value=completed):
        assert is_console_exe("fake.exe") is True


def test_is_console_exe_rejects_gui():
    """GUI はコンソール出力を持たず、--help でウィンドウを開いて返らない。"""
    from metabolomix.console.runner import is_console_exe
    with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("x", 15)):
        assert is_console_exe("gui.exe") is False


def test_is_console_exe_rejects_missing_file():
    from metabolomix.console.runner import is_console_exe
    with patch("subprocess.run", side_effect=OSError("not found")):
        assert is_console_exe("nope.exe") is False


def test_build_msdial_cmd_passes_project_flag(tmp_path):
    """save_project=True のとき -p を渡す（GUI で開ける .mdproject を出させる）。"""
    cmd = build_msdial_cmd("fake.exe", tmp_path, tmp_path / "out",
                           tmp_path / "params.txt", save_project=True)
    assert cmd[-1] == "-p"


def test_build_msdial_cmd_omits_project_flag_by_default(tmp_path):
    cmd = build_msdial_cmd("fake.exe", tmp_path, tmp_path / "out",
                           tmp_path / "params.txt")
    assert "-p" not in cmd


def test_build_msdial_cmd_is_the_only_argument_layout(tmp_path):
    """-i / -o / -m の並びはここ 1 か所でしか組まない（経路で分裂させない）。"""
    cmd = build_msdial_cmd("fake.exe", tmp_path / "raw", tmp_path / "out",
                           tmp_path / "params.txt")
    assert cmd[:2] == ["fake.exe", "lcms"]
    assert cmd[2::2] == ["-i", "-o", "-m"]


# ---------- output_collector ----------

def test_snapshot_empty(tmp_path):
    assert snapshot(tmp_path) == {}


def test_snapshot_captures_files(tmp_path):
    (tmp_path / "a.txt").write_text("hello")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b.txt").write_text("world")
    s = snapshot(tmp_path)
    assert "a.txt" in s
    assert str(Path("sub") / "b.txt") in s


def test_snapshot_excludes_named_dirs(tmp_path):
    (tmp_path / "runs" / "job1").mkdir(parents=True)
    (tmp_path / "runs" / "job1" / "analysis-job.json").write_text("{}", encoding="utf-8")
    (tmp_path / "S1.pai2").write_bytes(b"x")
    snap = snapshot(tmp_path, exclude_dir_names={"runs"})
    assert "S1.pai2" in snap
    assert not any("runs" in k for k in snap)


def test_collect_artifacts_tags_root_per_source(tmp_path):
    """-o のエクスポートと生データフォルダの生成物を 1 回で集め、出所を刻む。"""
    run_dir = tmp_path / "runs" / "job1"
    (run_dir / "msdial").mkdir(parents=True)
    (run_dir / "msdial" / "AlignResult-1.mzTab").write_text("MTD\n", encoding="utf-8")
    (run_dir / "msdial" / "S1.mdpeak").write_bytes(b"a")
    (tmp_path / "S1_1.pai2").write_bytes(b"b")
    (tmp_path / "S1_1.dcl").write_bytes(b"c")

    roots = {"run_dir": run_dir, "dataset_root": tmp_path}
    befores = {"run_dir": {}, "dataset_root": {}}
    entries, artifacts = collect_artifacts(
        roots, befores, declared_polarity="negative", declared_measure="peak_height")

    assert [(e.path, e.root) for e in entries] == [
        (str(Path("msdial") / "AlignResult-1.mzTab"), "run_dir")]
    by_role = {a.role: a for a in artifacts}
    assert by_role["sample_peak_table"].root == "run_dir"
    assert by_role["sample_peaks"].root == "dataset_root"
    assert by_role["msms_evidence"].root == "dataset_root"


def test_collect_artifacts_does_not_double_count_run_dir(tmp_path):
    """run_dir は dataset_root の配下にある。同じファイルを 2 回集めない。"""
    run_dir = tmp_path / "runs" / "job1"
    (run_dir / "msdial").mkdir(parents=True)
    (run_dir / "msdial" / "S1.mdpeak").write_bytes(b"a")
    roots = {"run_dir": run_dir, "dataset_root": tmp_path}
    befores = {"run_dir": {}, "dataset_root": snapshot(tmp_path, exclude_dir_names={"runs"})}
    _, artifacts = collect_artifacts(roots, befores)
    assert [a.path for a in artifacts] == [str(Path("msdial") / "S1.mdpeak")]


@pytest.mark.parametrize("filename,expected_role,expected_fmt", [
    ("msdial/S1.mdpeak", "sample_peak_table", "mdpeak"),
    ("msdial/AlignResult-2026931617.mdalign", "alignment_table", "mdalign"),
    ("msdial/S1.mdmsp", "msms_spectra", "mdmsp"),
    ("msdial/AlignResult-2026931617.qa.tsv", "quality_matrix", "qatsv"),
    ("msdial/Project-2609030417.mdproject", "gui_project", "mdproject"),
    ("Project-2609030417.mddata", "project_data", "mddata"),
    ("S1_2026931617_tags.xml", "peak_tags", "tagsxml"),
    ("S1_2026931617.pai2", "sample_peaks", "pai2"),
])
def test_assign_role_console_outputs(filename, expected_role, expected_fmt):
    role, fmt = _assign_role(filename)
    assert (role, fmt) == (expected_role, expected_fmt)


def test_collect_artifacts_skips_hash_for_unknown_role(tmp_path):
    """role の付かないファイルはハッシュしない。"""
    (tmp_path / "mystery.bin").write_bytes(b"x" * 1024)
    (tmp_path / "S1.pai2").write_bytes(b"y" * 16)
    _, artifacts = collect_artifacts({"run_dir": tmp_path}, {"run_dir": {}})
    by_path = {a.path: a for a in artifacts}
    assert by_path["mystery.bin"].role == "unknown"
    assert by_path["mystery.bin"].sha256 == ""
    assert by_path["S1.pai2"].sha256 != ""


def test_collect_artifacts_detects_new_files(tmp_path):
    before = snapshot(tmp_path)
    mzdial_out = tmp_path / "msdial"
    mzdial_out.mkdir()
    (mzdial_out / "Height_AlignmentResult_001.mzTab").write_text("MTD\t")
    (mzdial_out / "AlignmentResult_001.arf").write_bytes(b"\x00" * 10)

    mztabs, others = collect_artifacts({"run_dir": tmp_path}, {"run_dir": before})
    assert len(mztabs) == 1
    assert len(others) == 1
    assert others[0].role == "peak_matrix_source"
    assert others[0].format == "arf"


def test_collect_artifacts_ignores_unchanged(tmp_path):
    (tmp_path / "existing.txt").write_text("old")
    before = snapshot(tmp_path)
    mztabs, others = collect_artifacts({"run_dir": tmp_path}, {"run_dir": before})
    assert mztabs == []
    assert others == []


@pytest.mark.parametrize("filename,expected_role,expected_fmt", [
    ("Height_AlignmentResult.mzTab", "primary_mztab", "mztab"),
    ("AlignmentResult.arf2",         "spot_catalog",   "arf2"),
    ("AlignmentResult.arf",          "peak_matrix_source", "arf"),
    ("sample01.pai2",                "sample_peaks",   "pai2"),
    ("sample01.dcl",                 "msms_evidence",  "dcl"),
    ("sample01.EIC.aef",             "chromatogram",   "eicaef"),
    ("unknown.bin",                  "unknown",        "bin"),
])
def test_assign_role(filename, expected_role, expected_fmt):
    role, fmt = _assign_role(filename)
    assert role == expected_role
    assert fmt == expected_fmt


@pytest.mark.parametrize("filename,expected_polarity,expected_measure", [
    ("Height_AlignmentResult_Neg.mzTab", "negative", "peak_height"),
    ("Area_AlignmentResult_Pos.mzTab",   "positive", "peak_area_above_zero"),
    # 信号ゼロのファイル名。既定値へ落とすのは呼び出し側（collect_artifacts）の
    # 判断で、推定関数は「何も語っていない」を None で返す。
    ("AlignmentResult.mzTab",            None,       None),
])
def test_infer_mztab_meta(filename, expected_polarity, expected_measure):
    polarity, measure = _infer_mztab_meta(filename)
    assert polarity == expected_polarity
    assert measure == expected_measure


# ---------- job_manager ----------

def test_create_job_writes_json(tmp_path):
    method = tmp_path / "params.msdial"
    method.touch()
    job, job_path = create_job(
        dataset_root=tmp_path,
        method_file=method,
        polarity="negative",
        measure="peak_height",
    )
    assert job_path.is_file()
    assert job.status == "planned"
    assert job.polarity == "negative"


def test_create_job_run_dir_inside_dataset_root(tmp_path):
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(
        dataset_root=tmp_path,
        method_file=method,
        polarity="positive",
        measure="peak_height",
    )
    assert str(tmp_path) in str(job_path)


def test_update_status(tmp_path):
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(
        dataset_root=tmp_path, method_file=method,
        polarity="positive", measure="peak_height",
    )
    job = update_status(job_path, "completed")
    assert job.status == "completed"
    assert load_job(job_path).status == "completed"


def test_update_status_with_error(tmp_path):
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(
        dataset_root=tmp_path, method_file=method,
        polarity="positive", measure="peak_height",
    )
    job = update_status(job_path, "failed", error="timeout")
    assert job.error == "timeout"
    assert load_job(job_path).error == "timeout"


def test_list_jobs_empty(tmp_path):
    assert list_jobs(tmp_path) == []


def test_list_jobs_finds_jobs(tmp_path):
    method = tmp_path / "params.msdial"
    method.touch()
    _, p1 = create_job(dataset_root=tmp_path, method_file=method,
                       polarity="positive", measure="peak_height")
    _, p2 = create_job(dataset_root=tmp_path, method_file=method,
                       polarity="negative", measure="peak_height")
    found = list_jobs(tmp_path)
    assert len(found) == 2
    assert p1 in found or p2 in found


def test_raw_input_summary_counts_by_extension(tmp_path):
    from metabolomix.console.job_manager import raw_input_summary
    (tmp_path / "a.wiff").touch()
    (tmp_path / "b.wiff").touch()
    (tmp_path / "a.wiff2").touch()
    (tmp_path / "note.txt").touch()
    (tmp_path / "c.d").mkdir()
    assert raw_input_summary(tmp_path) == {"wiff": 2, "wiff2": 1, "d": 1}


def test_count_raw_inputs_totals_summary(tmp_path):
    (tmp_path / "a.wiff").touch()
    (tmp_path / "a.abf").touch()
    assert count_raw_inputs(tmp_path) == 2


def test_vendor_directories_count_as_measurements(tmp_path):
    """`.d` / `.raw` の**フォルダ**は 1 検体。MS-DIAL がそう読む。

    上流 `AnalysisFilesParser.ReadFolderContents` の `isVendorDirectory`
    （`Directory.Exists(path) && (extension == ".raw" || extension == ".d")`）。
    Agilent/Bruker の `.d` と Waters の `.raw` がこれに当たる。
    """
    from metabolomix.console.job_manager import list_raw_inputs, raw_input_summary
    (tmp_path / "a.d").mkdir()
    (tmp_path / "b.d").mkdir()
    (tmp_path / "waters.raw").mkdir()
    assert raw_input_summary(tmp_path) == {"d": 2, "raw": 1}
    assert [p.name for p in list_raw_inputs(tmp_path)] == ["a.d", "b.d", "waters.raw"]


def test_non_vendor_directories_are_not_measurements(tmp_path):
    """`.d` / `.raw` 以外はフォルダだと計測データにならない。

    `DataAccess.IsDataFormatSupported` は `File.Exists` を要求する。名前だけ
    `.mzml` のフォルダを 1 検体と数えると、`input_count` と `raw_inventory` に
    実在しないサンプルが入り、`ms_run[N]-location` との 1 対 1 照合が
    実行の最後になって落ちる。
    """
    from metabolomix.console.job_manager import raw_input_summary
    (tmp_path / "backup.mzml").mkdir()
    (tmp_path / "old.wiff").mkdir()
    (tmp_path / "real.mzml").touch()
    assert raw_input_summary(tmp_path) == {"mzml": 1}


def test_all_msdial_raw_extensions_are_recognised(tmp_path):
    """MS-DIAL の `SupportMsRawDataExtension` 12 形式をそのまま数える。

    出典: `src/MSDIAL5/MsdialCore/Enum/SupportFormat.cs`
    `enum SupportMsRawDataExtension { abf, ibf, cdf, mzml, wiff, raw, d, wiff2, qgd, lcd, lrp, imzml }`
    """
    from metabolomix.console.job_manager import raw_input_summary
    extensions = ["abf", "ibf", "cdf", "mzml", "wiff", "raw",
                  "d", "wiff2", "qgd", "lcd", "lrp", "imzml"]
    for ext in extensions:
        (tmp_path / f"sample.{ext}").touch()
    assert raw_input_summary(tmp_path) == {ext: 1 for ext in extensions}



def _stub_lbm(tmp_path, monkeypatch):
    """脂質ライブラリを 1 件用意する。

    console_plan は lipidomics で LBM を解決できなければ停止する（GUI と同じ規則）。
    LBM そのものを見ていないテストは、ここで最小の 1 件を与えて本題に集中する。
    LBM 解決の検証は tests/test_console_plan_method.py。
    """
    lbm = tmp_path / "stub.lbm2"
    lbm.touch()
    monkeypatch.setenv("MSDIAL_LBM", str(lbm))
    return lbm


# 偽 Console のヘルパの正準は tests/pipeline_fixtures.py（他のテストと共有する）。
from tests.pipeline_fixtures import (  # noqa: E402
    fake_console_command as _fake_console,
    mztab_text as _mztab_text,
    use_fake_console as _use_console,
)


# ---------- console_tools ----------

def test_console_plan_unsupported_area(tmp_path):
    import json as _json
    from metabolomix.tools.console_tools import console_plan
    result = console_plan(
        dataset_root=str(tmp_path),
        method_file=str(tmp_path / "m.msdial"),
        polarity="positive",
        measure="peak_area_above_zero",
    )
    parsed = _json.loads(result)
    assert parsed["error"]["code"] == "UNSUPPORTED_AREA_CONSOLE"


def test_console_plan_missing_exe(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.delenv("MSDIAL_EXE", raising=False)
    method = tmp_path / "params.msdial"
    method.write_text("Ion mode: Positive\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    result = console_plan(
        dataset_root=str(tmp_path),
        method_file=str(method),
        polarity="positive",
        measure="peak_height",
    )
    parsed = _json.loads(result)
    assert parsed["error"]["code"] == "MSDIAL_EXE_NOT_FOUND"


def test_console_plan_missing_method_file(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    from metabolomix.tools.console_tools import console_plan
    result = console_plan(
        dataset_root=str(tmp_path),
        method_file=str(tmp_path / "nonexistent.msdial"),
        polarity="positive",
        measure="peak_height",
    )
    parsed = _json.loads(result)
    assert parsed["error"]["code"] == "METHOD_FILE_NOT_FOUND"


def test_console_plan_rejects_non_console_exe(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "gui.exe")
    _stub_lbm(tmp_path, monkeypatch)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    (tmp_path / "a.wiff").touch()
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: False)
    from metabolomix.tools.console_tools import console_plan
    result = console_plan(dataset_root=str(tmp_path), method_file=str(method),
                          polarity="negative", measure="peak_height")
    assert _json.loads(result)["error"]["code"] == "MSDIAL_EXE_NOT_CONSOLE"


def test_console_plan_success(tmp_path, monkeypatch):
    import json as _json
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.msdial"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    (tmp_path / "a.wiff").touch()
    from metabolomix.tools.console_tools import console_plan
    result = console_plan(
        dataset_root=str(tmp_path),
        method_file=str(method),
        polarity="negative",
        measure="peak_height",
    )
    parsed = _json.loads(result)
    assert parsed["status"] == "planned"
    assert parsed["polarity"] == "negative"
    assert session_state.session.current_job_path is not None


def test_console_plan_rejects_mixed_raw_formats(tmp_path, monkeypatch):
    """.wiff と .wiff2 の併存は MS-DIAL が対話プロンプトを出す条件。"""
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    (tmp_path / "a.wiff").touch()
    (tmp_path / "a.wiff2").touch()
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "MIXED_RAW_FORMATS"
    assert parsed["error"]["details"]["formats"] == {"wiff": 1, "wiff2": 1}


def test_console_plan_rejects_empty_dataset_root(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "MIXED_RAW_FORMATS"


def test_console_plan_warns_existing_alignment_results(tmp_path, monkeypatch):
    """実行のたびに dataset_root へ別タイムスタンプのアライメント一式が積まれる。"""
    import json as _json
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    (tmp_path / "a.wiff").touch()
    (tmp_path / "AlignResult-2026931617_PeakProperties.arf").touch()
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["status"] == "planned"
    assert any("既存のアライメント結果" in w for w in parsed["warnings"])


def test_console_plan_rejects_binary_method_file(tmp_path, monkeypatch):
    """.mdproject は ZIP。渡すと MS-DIAL は全パラメータ既定値で走ってしまう。"""
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "a.wiff").touch()
    method = tmp_path / "project.mdproject"
    method.write_bytes(b"PK\x03\x04\x00\x00binary")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_FILE_NOT_TEXT"


def test_console_plan_rejects_text_without_key_value(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "a.wiff").touch()
    method = tmp_path / "empty.txt"
    method.write_text("# comment only\n\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["error"]["code"] == "METHOD_FILE_NOT_TEXT"


def test_console_plan_accepts_key_value_method_file(tmp_path, monkeypatch):
    import json as _json
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "a.wiff").touch()
    method = tmp_path / "params.txt"
    method.write_text("# MS-DIAL param\nIon mode: Negative\n", encoding="ascii")
    from metabolomix.tools.console_tools import console_plan
    parsed = _json.loads(console_plan(dataset_root=str(tmp_path), method_file=str(method),
                                      polarity="negative", measure="peak_height"))
    assert parsed["status"] == "planned"


def test_console_status_no_job():
    import json as _json
    from metabolomix.core import session_state
    session_state.session = session_state.AnalysisSession()
    from metabolomix.tools.console_tools import console_status
    result = console_status()
    parsed = _json.loads(result)
    assert parsed["error"]["code"] == "JOB_NOT_FOUND"


def test_console_run_nonplanned_job(tmp_path, monkeypatch):
    import json as _json
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(
        dataset_root=tmp_path, method_file=method,
        polarity="positive", measure="peak_height",
    )
    update_status(job_path, "completed")
    from metabolomix.tools.console_tools import console_run
    result = console_run(str(job_path))
    parsed = _json.loads(result)
    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"


def test_console_run_rejects_changed_non_console_exe_before_launch(tmp_path, monkeypatch):
    """計画後に MSDIAL_EXE が GUI へ変わっても、ジョブを planned のまま止める。"""
    import json as _json
    from metabolomix.core import session_state
    from metabolomix.console import runner
    from metabolomix.tools.console_tools import console_plan, console_run

    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "console.exe")
    _stub_lbm(tmp_path, monkeypatch)
    real_is_console_exe = runner.is_console_exe
    monkeypatch.setattr(runner, "is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    (tmp_path / "a.wiff").touch()
    planned = _json.loads(console_plan(
        dataset_root=str(tmp_path), method_file=str(method),
        polarity="negative", measure="peak_height"))
    job_path = planned["job_path"]

    monkeypatch.setenv("MSDIAL_EXE", "MSDIAL.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr(runner, "is_console_exe", real_is_console_exe)
    completed = MagicMock()
    completed.stdout = "MS-DIAL GUI output without console commands\n"
    with patch("subprocess.run", return_value=completed) as mock_run:
        result = _json.loads(console_run(job_path))

    assert result["error"]["code"] == "MSDIAL_EXE_NOT_CONSOLE"
    assert load_job(Path(job_path)).status == "planned"
    mock_run.assert_called_once_with(
        ["MSDIAL.exe", "--help"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL,
        timeout=15,
        text=True,
        errors="replace",
    )


def test_console_run_marks_missing_changed_exe_failed(tmp_path, monkeypatch):
    """計画後に MSDIAL_EXE が消えた場合は従来どおり NOT_FOUND/failed にする。"""
    import json as _json
    from metabolomix.console.job_manager import create_job
    from metabolomix.tools.console_tools import console_run

    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    job, job_path = create_job(
        dataset_root=tmp_path, method_file=method,
        polarity="negative", measure="peak_height",
    )
    monkeypatch.setenv("MSDIAL_EXE", "removed-console.exe")
    _stub_lbm(tmp_path, monkeypatch)
    with patch("subprocess.run", side_effect=OSError("exe not found")) as mock_run:
        result = _json.loads(console_run(str(job_path)))

    assert result["error"]["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert load_job(job_path).status == "failed"
    mock_run.assert_called_once()
    assert mock_run.call_args.args[0][1] == "--help"


# ---------- Task 0: ブロッカー回帰 ----------

def test_create_job_allows_default_data_dir(tmp_path, monkeypatch):
    """既定データディレクトリ（<repo>/data 配下）を dataset_root にできる。"""
    from metabolomix.core import mcp_core
    from metabolomix.console.job_manager import create_job
    # リポジトリ直下の data/ を模す: BASE_DIR 配下だが DATA_DIR 配下でもある
    fake_repo = tmp_path / "repo"
    data_root = fake_repo / "data" / "study-001"
    data_root.mkdir(parents=True)
    method = tmp_path / "params.msdial"
    method.touch()
    monkeypatch.setattr(mcp_core, "BASE_DIR", fake_repo)
    monkeypatch.setenv("LIPIDMIX_DATA_DIR", str(fake_repo / "data"))

    job, job_path = create_job(dataset_root=data_root, method_file=method,
                              polarity="positive", measure="peak_height")
    assert job_path.is_file()


def test_create_job_still_rejects_source_tree(tmp_path, monkeypatch):
    """データディレクトリ外のリポジトリ内パスは従来どおり拒否する。"""
    from metabolomix.core import mcp_core
    from metabolomix.console.job_manager import create_job
    fake_repo = tmp_path / "repo"
    (fake_repo / "metabolomix").mkdir(parents=True)
    (fake_repo / "data").mkdir(parents=True)
    method = tmp_path / "params.msdial"
    method.touch()
    monkeypatch.setattr(mcp_core, "BASE_DIR", fake_repo)
    monkeypatch.setenv("LIPIDMIX_DATA_DIR", str(fake_repo / "data"))

    with pytest.raises(ValueError):
        create_job(dataset_root=fake_repo / "metabolomix", method_file=method,
                   polarity="positive", measure="peak_height")


def test_collect_artifacts_excludes_operational_files(tmp_path):
    """msdial.log / analysis-job.json は生成物として数えない。"""
    from metabolomix.console.output_collector import collect_artifacts, snapshot
    before = snapshot(tmp_path)
    (tmp_path / "msdial.log").write_text("CMD: fake\n")
    (tmp_path / "analysis-job.json").write_text("{}")
    mztabs, others = collect_artifacts({"run_dir": tmp_path}, {"run_dir": before})
    assert mztabs == []
    assert others == []


def test_console_run_reports_no_output(tmp_path, monkeypatch):
    """終了コード 0 でも生成物が 1 件も無ければ NO_JOB_OUTPUT になる。"""
    import json as _json
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="positive", measure="peak_height")
    _use_console(monkeypatch, _fake_console())

    parsed = _json.loads(console_run(str(job_path)))

    assert parsed["error"]["code"] == "NO_JOB_OUTPUT"
    assert load_job(job_path).status == "failed"


def test_console_run_marks_failed_when_the_executable_cannot_be_started(tmp_path, monkeypatch):
    """実在しない exe で起動に失敗しても封筒を返し、running に固着させない。

    `is_console_exe` の事前確認をすり抜けた場合（計画後に exe が消える等）の
    最後の砦。起動できなかった事実は終了証跡に termination=launch_failed として残る。
    """
    import json as _json
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", str(tmp_path / "definitely_not_here.exe"))
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="positive", measure="peak_height")

    parsed = _json.loads(console_run(str(job_path)))

    assert parsed["error"]["code"] == "MSDIAL_EXE_NOT_FOUND"
    assert parsed["error"]["details"]["termination"] == "launch_failed"
    assert load_job(job_path).status == "failed"


# ---------- F1: 実行後処理（collect_artifacts 以降）のガード ----------

def test_console_run_post_run_failure_marks_job_failed_not_running(tmp_path, monkeypatch):
    """Console 自体は成功しても、収集で例外が出たら running に固着させない。

    Windows でウイルススキャナ等が生成直後のファイルをロックしていると
    `sha256_file` が PermissionError を投げる実運用の再現。ここが無防備だと
    analysis-job.json は running のまま残り、以降の console_run は全部
    JOB_NOT_PLANNED で拒否されて誰も直せなくなる。
    """
    import json as _json
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.msdial"
    method.touch()
    job, job_path = create_job(dataset_root=tmp_path, method_file=method,
                               polarity="positive", measure="peak_height")
    out = Path(job.run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console(
        {out / "Height_AlignmentResult_ctrl_1.mzTab": "MTD\t"}))
    monkeypatch.setattr("metabolomix.console.output_collector.sha256_file",
                        lambda *a, **k: (_ for _ in ()).throw(
                            PermissionError("locked by AV scanner")))

    parsed = _json.loads(console_run(str(job_path)))

    assert parsed["error"]["code"] == "JOB_POST_RUN_FAILED"
    reloaded = load_job(job_path)
    assert reloaded.status == "failed"
    assert reloaded.error


# ---------- 成果物メタの正準化（job の宣言値 vs ファイル名の推定） ----------
# 背景: MS-DIAL のアライメント出力名（Height_AlignmentResult_<timestamp>.mzTab）は
# 極性トークンを持たない。旧実装は「信号なし」を positive と区別せず既定に落として
# いたため、negative で計画したジョブの mzTab エントリが全て positive と記録された。

def test_infer_mztab_meta_returns_none_when_filename_has_no_polarity_token():
    polarity, measure = _infer_mztab_meta("Height_AlignmentResult_2026_05_15.mzTab")
    assert polarity is None       # 「信号なし」。positive と断定してはいけない
    assert measure == "peak_height"


def test_infer_mztab_meta_ignores_polarity_substring_inside_a_word():
    polarity, _ = _infer_mztab_meta("Height_Negev_cohort_AlignmentResult.mzTab")
    assert polarity is None       # "Negev" の neg は極性トークンではない


def test_infer_mztab_meta_returns_none_measure_for_normalized_prefix():
    # spec §8.1: normalized value は未対応。peak_height と偽ってはいけない。
    _, measure = _infer_mztab_meta("NormalizedHeight_AlignmentResult_2026.mzTab")
    assert measure is None


def test_collect_artifacts_fills_polarity_from_job_declaration(tmp_path):
    (tmp_path / "Height_AlignmentResult_2026.mzTab").write_text("x", encoding="utf-8")
    mztabs, _ = collect_artifacts(
        {"run_dir": tmp_path}, {"run_dir": {}},
        declared_polarity="negative", declared_measure="peak_height")
    assert len(mztabs) == 1
    assert mztabs[0].polarity == "negative"
    assert mztabs[0].validation["polarity_source"] == "job_declared"


def test_collect_artifacts_prefers_filename_polarity_over_declaration_and_records_conflict(tmp_path):
    (tmp_path / "Height_AlignmentResult_Neg.mzTab").write_text("x", encoding="utf-8")
    mztabs, _ = collect_artifacts(
        {"run_dir": tmp_path}, {"run_dir": {}},
        declared_polarity="positive", declared_measure="peak_height")
    # ファイルは実物の性質を語る。宣言は意図でしかないので、食い違いは
    # ファイル側を採ったうえで記録する（黙って上書きしない）。
    assert mztabs[0].polarity == "negative"
    assert mztabs[0].validation["polarity_source"] == "filename"
    assert "polarity" in mztabs[0].validation["conflicts"]


def test_collect_artifacts_excludes_normalized_mztab_from_primary_candidates(tmp_path):
    (tmp_path / "NormalizedHeight_AlignmentResult.mzTab").write_text("x", encoding="utf-8")
    mztabs, others = collect_artifacts(
        {"run_dir": tmp_path}, {"run_dir": {}},
        declared_polarity="negative", declared_measure="peak_height")
    assert mztabs == []
    assert [a.role for a in others] == ["unsupported_mztab"]


def test_console_run_records_declared_polarity_on_mztab_entries(tmp_path, monkeypatch):
    """極性トークンを持たない出力名でも、宣言した negative が記録される。"""
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.abf").write_bytes(b"raw")
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignmentResult_2026.mzTab":
            _mztab_text(tmp_path, [tmp_path / "S1.abf"])}))

    console_run(str(job_path))

    saved = load_job(job_path)
    assert saved.status == "completed"
    assert [e.polarity for e in saved.primary_mztab_files] == ["negative"]


def test_console_run_warns_when_filename_contradicts_declared_polarity(tmp_path, monkeypatch):
    """ファイル名の極性が宣言と食い違うなら、黙らずに warning へ残す。"""
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.abf").write_bytes(b"raw")
    method = tmp_path / "params.msdial"
    method.touch()
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignmentResult_Pos.mzTab":
            _mztab_text(tmp_path, [tmp_path / "S1.abf"])}))

    console_run(str(job_path))

    saved = load_job(job_path)
    assert [e.polarity for e in saved.primary_mztab_files] == ["positive"]
    assert any("polarity" in w for w in saved.warnings)


# ---------- サイドカー廃止（案 c）: 消費者のいない生成物を出さない ----------
# feature-qc.tsv は「単なる記録」として作られたが、リポジトリ内外に読み手が
# 一人もおらず（別リポ massbank-context / Use-LLLM も参照ゼロ）、書いていた
# 内容（role / batch）は build_dataset_pp_inputs がサンプル名から独立に
# 再導出している。spec §10.1 が要求する feature-qc（SMF_ID × assay_id の
# is_gap_filled / detected_peak）とは別物なので、廃止しても仕様上の穴は
# 増えない（その穴は最初から埋まっていない）。

def test_console_run_writes_no_sidecar_files(tmp_path, monkeypatch):
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.abf").write_bytes(b"raw")
    method = tmp_path / "params.msdial"
    method.touch()
    job, job_path = create_job(dataset_root=tmp_path, method_file=method,
                               polarity="positive", measure="peak_height")
    run_dir = Path(job.run_dir)
    _use_console(monkeypatch, _fake_console({
        run_dir / "msdial" / "Height_AlignmentResult_1.mzTab":
            _mztab_text(tmp_path, [tmp_path / "S1.abf"]),
        tmp_path / "ctrl_1.pai2": "peaks"}))

    assert json.loads(console_run(str(job_path)))["status"] == "completed"

    assert not (run_dir / "sidecars").exists()
    assert "sample_qc_sidecar" not in {a.role for a in load_job(job_path).artifacts}


def test_console_run_warns_when_no_per_sample_pai2(tmp_path, monkeypatch):
    """.pai2 が 1 つも出なければ warning を残す（MS/MS 経路が空になる合図）。

    サイドカー生成の副産物として出ていた警告だが、価値はサイドカーではなく
    「サンプル別ファイルが無い＝pai2_parser / dcl_find_msms が読むものが無い」
    という事実の側にある。廃止後も独立した検査として残す。
    """
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.abf").write_bytes(b"raw")
    method = tmp_path / "params.msdial"
    method.touch()
    job, job_path = create_job(dataset_root=tmp_path, method_file=method,
                               polarity="positive", measure="peak_height")
    _use_console(monkeypatch, _fake_console({
        Path(job.run_dir) / "msdial" / "Height_Alignment.mzTab":
            _mztab_text(tmp_path, [tmp_path / "S1.abf"])}))

    assert json.loads(console_run(str(job_path)))["status"] == "completed"

    warnings = load_job(job_path).warnings
    assert any(".pai2" in w for w in warnings)
    assert not any("feature-qc" in w for w in warnings)


def test_console_run_no_pai2_warning_when_pai2_present(tmp_path, monkeypatch):
    """.pai2 が出ているときに上の警告を出してはいけない（過剰警告の防止）。"""
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.abf").write_bytes(b"raw")
    method = tmp_path / "params.msdial"
    method.touch()
    job, job_path = create_job(dataset_root=tmp_path, method_file=method,
                               polarity="positive", measure="peak_height")
    _use_console(monkeypatch, _fake_console({
        Path(job.run_dir) / "msdial" / "Height_Alignment.mzTab":
            _mztab_text(tmp_path, [tmp_path / "S1.abf"]),
        tmp_path / "ctrl_1.pai2": "peaks"}))

    console_run(str(job_path))

    assert not any(".pai2" in w for w in load_job(job_path).warnings)


# ---------- Task 8: Console 実行時の二重ルート収集と実行オプション ----------

def _planned_task8_job(tmp_path, monkeypatch, **kwargs):
    """console_plan を通して Task 8 のジョブを 1 件作り、パスを返す。"""
    import json as _json
    from metabolomix.core import session_state
    from metabolomix.tools.console_tools import console_plan

    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    (tmp_path / "S1.wiff").touch()
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    parsed = _json.loads(console_plan(
        dataset_root=str(tmp_path), method_file=str(method), polarity="negative",
        measure="peak_height", **kwargs))
    assert parsed["status"] == "planned"
    return Path(parsed["job_path"])


def test_console_run_collects_dataset_root_outputs(tmp_path, monkeypatch):
    """生データフォルダの .pai2/.dcl も収集し、誤った不在警告を出さない。"""
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"]),
        out / "S1.mdpeak": "a",
        tmp_path / "S1_1.pai2": "b",
        tmp_path / "S1_1.dcl": "c"}))

    parsed = _json.loads(console_run(str(job_path)))

    assert parsed["status"] == "completed"
    assert not any(".pai2" in warning for warning in parsed["warnings"])
    saved = load_job(job_path)
    assert (saved.save_project, saved.timeout_s) == (True, 21600)
    roots = {artifact.role: artifact.root for artifact in saved.artifacts}
    assert roots["sample_peaks"] == "dataset_root"
    assert roots["msms_evidence"] == "dataset_root"
    assert roots["sample_peak_table"] == "run_dir"


def test_console_run_warns_when_no_dataset_root_sample_files(tmp_path, monkeypatch):
    """両ルートに .pai2 がなければ初めて MS/MS 根拠不足を警告する。"""
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"])}))

    parsed = _json.loads(console_run(str(job_path)))
    assert any(".pai2" in warning for warning in parsed["warnings"])


def test_console_plan_persists_and_console_run_passes_execution_options(tmp_path, monkeypatch):
    """計画時の timeout/save_project が永続化され、実行へそのまま渡る。

    timeout は Console のコマンドラインではなく監視側の期限なので、渡ったことは
    終了証跡の `timeout_s` で確かめる（証跡は監視が実際に使った値を書く）。
    """
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch, save_project=True, timeout_s=1234)
    out = Path(load_job(job_path).run_dir) / "msdial"
    command = _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"])})
    seen = {}

    def _capture(exe, dataset_root, msdial_out_dir, method_file, save_project=False):
        seen.update(exe=exe, save_project=save_project,
                    dataset_root=str(dataset_root), out=str(msdial_out_dir))
        return command

    monkeypatch.setattr("metabolomix.console.runner.build_msdial_cmd", _capture)
    console_run(str(job_path))

    saved = load_job(job_path)
    assert (saved.save_project, saved.timeout_s) == (True, 1234)
    assert seen["exe"] == "fake.exe"
    assert seen["save_project"] is True
    assert seen["dataset_root"] == str(tmp_path)
    receipt = _json.loads(
        (Path(saved.run_dir) / "execution-result.json").read_text(encoding="utf-8"))
    assert receipt["timeout_s"] == 1234


@pytest.mark.parametrize("save_project, timeout_s", [
    ("true", 1), (1, 1), (True, 0), (True, -1), (True, 1.5), (True, True),
])
def test_console_plan_rejects_invalid_execution_options_before_job_side_effects(
        tmp_path, monkeypatch, save_project, timeout_s):
    """不正な実行オプションは runs 作成も current_job_path 更新も起こさせない。"""
    import json as _json
    from metabolomix.core import session_state
    from metabolomix.tools.console_tools import console_plan

    session_state.session = session_state.AnalysisSession()
    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    _stub_lbm(tmp_path, monkeypatch)
    (tmp_path / "S1.wiff").touch()
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")

    parsed = _json.loads(console_plan(
        dataset_root=str(tmp_path), method_file=str(method), polarity="negative",
        measure="peak_height", save_project=save_project, timeout_s=timeout_s))

    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"
    assert session_state.session.current_job_path is None
    assert not (tmp_path / "runs").exists()


def test_console_run_rejects_invalid_loaded_execution_options_before_running(tmp_path, monkeypatch):
    """手編集された不正ジョブは Console を起動せず planned のまま拒否する。"""
    import json as _json
    from metabolomix.console.job_manager import load_job, save_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    job = load_job(job_path)
    job.timeout_s = 0
    save_job(job, job_path)
    monkeypatch.setattr(
        "metabolomix.console.execution.start_owned_process",
        lambda *args, **kwargs: pytest.fail("不正な実行オプションで Console を起動した"),
    )

    parsed = _json.loads(console_run(str(job_path)))
    assert parsed["error"]["code"] == "JOB_NOT_PLANNED"
    assert load_job(job_path).status == "planned"


def test_console_run_timeout_persists_partial_outputs_from_both_roots(tmp_path, monkeypatch):
    """タイムアウト後も既出力は partial ジョブに保存し、根を失わない。"""
    import json as _json
    from metabolomix.console.job_manager import load_job, save_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    job = load_job(job_path)
    job.timeout_s = 1
    save_job(job, job_path)
    out = Path(job.run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console(
        {out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"]),
         tmp_path / "S1_1.pai2": "partial"},
        sleep_s=30))

    parsed = _json.loads(console_run(str(job_path)))
    saved = load_job(job_path)

    assert parsed["error"]["code"] == "MSDIAL_TIMEOUT"
    assert parsed["error"]["details"]["status"] == "partial"
    assert saved.status == "partial"
    assert "termination=timeout" in (saved.error or "")
    assert {entry.root for entry in saved.primary_mztab_files} == {"run_dir"}
    assert {artifact.root for artifact in saved.artifacts} == {"dataset_root"}


def test_console_run_timeout_without_outputs_fails(tmp_path, monkeypatch):
    """タイムアウト時に出力ゼロなら partial と偽らず failed にする。"""
    import json as _json
    from metabolomix.console.job_manager import load_job, save_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    job = load_job(job_path)
    job.timeout_s = 1
    save_job(job, job_path)
    _use_console(monkeypatch, _fake_console(sleep_s=30))

    parsed = _json.loads(console_run(str(job_path)))
    assert parsed["error"]["code"] == "MSDIAL_TIMEOUT"
    assert parsed["error"]["details"]["status"] == "failed"
    assert load_job(job_path).status == "failed"


def test_console_run_timeout_collection_failure_marks_failed_with_timeout_context(tmp_path, monkeypatch):
    """タイムアウト後の収集例外でも running を残さず、両方の原因を保存する。

    封筒は timeout（Console がどう終わったか）を名乗り、収集が落ちた事実は
    終了証跡の collection と job の error に残る。**どちらか一方だけを記録して
    もう一方を消さない**。
    """
    import json as _json
    from metabolomix.console.job_manager import load_job, save_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    job = load_job(job_path)
    job.timeout_s = 1
    save_job(job, job_path)
    _use_console(monkeypatch, _fake_console(sleep_s=30))
    monkeypatch.setattr(
        "metabolomix.console.execution.collect_artifacts",
        lambda *args, **kwargs: (_ for _ in ()).throw(PermissionError("locked")),
    )

    parsed = _json.loads(console_run(str(job_path)))
    saved = load_job(job_path)
    receipt = _json.loads(
        (Path(saved.run_dir) / "execution-result.json").read_text(encoding="utf-8"))

    assert parsed["error"]["code"] == "MSDIAL_TIMEOUT"
    assert parsed["error"]["details"]["status"] == "failed"
    assert saved.status == "failed"
    assert "termination=timeout" in (saved.error or "")
    assert "COLLECTION_FAILED" in (saved.error or "")
    assert "locked" in receipt["collection"]["error"]


def test_console_run_persistence_failure_marks_successful_run_failed(tmp_path, monkeypatch):
    """収集済み成果物の保存に失敗しても running を残さず封筒を返す。"""
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"])}))
    monkeypatch.setattr(
        "metabolomix.console.execution.save_job",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("save locked")),
    )

    parsed = _json.loads(console_run(str(job_path)))
    saved = load_job(job_path)
    assert parsed["error"]["code"] == "JOB_POST_RUN_FAILED"
    assert saved.status == "failed"
    assert "save locked" in (saved.error or "")
    # 保存できなくても復旧できるだけの情報を封筒に残す。
    assert parsed["error"]["details"]["recovery"]["intended_status"] == "completed"


def test_console_run_timeout_persistence_failure_preserves_timeout_context(tmp_path, monkeypatch):
    """timeout 後の保存失敗も timeout 封筒と failed 状態に収束させる。"""
    import json as _json
    from metabolomix.console.job_manager import load_job, save_job
    from metabolomix.tools.console_tools import console_run

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    job = load_job(job_path)
    job.timeout_s = 1
    save_job(job, job_path)
    out = Path(job.run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console(
        {out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"])},
        sleep_s=30))
    monkeypatch.setattr(
        "metabolomix.console.execution.save_job",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("save locked")),
    )

    parsed = _json.loads(console_run(str(job_path)))
    saved = load_job(job_path)
    assert parsed["error"]["code"] == "MSDIAL_TIMEOUT"
    assert parsed["error"]["details"]["status"] == "failed"
    assert saved.status == "failed"
    assert "save locked" in (saved.error or "")


def test_console_status_exposes_roots_artifacts_and_execution_options(tmp_path, monkeypatch):
    """状態照会だけで生成物の由来と実行設定を追跡できる。"""
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run, console_status

    job_path = _planned_task8_job(tmp_path, monkeypatch, save_project=True, timeout_s=456)
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"]),
        tmp_path / "S1_1.pai2": "x"}))
    console_run(str(job_path))

    parsed = _json.loads(console_status(str(job_path), include_artifacts=True))

    assert parsed["dataset_root"] == str(tmp_path)
    assert parsed["execution"] == {"save_project": True, "timeout_s": 456}
    assert parsed["mztab_files"][0]["root"] == "run_dir"
    assert parsed["artifacts"] == (
        "path\trole\tformat\troot\n"
        "S1_1.pai2\tsample_peaks\tpai2\tdataset_root"
    )
    # 状態照会は保存済みの証跡を読むだけ。終わり方は execution_receipt に出る。
    assert parsed["execution_receipt"]["termination"] == "exited"
    assert parsed["execution_receipt"]["exit_code"] == 0


def test_method_text_check_reads_only_the_head(tmp_path, monkeypatch):
    """判定は先頭だけを読む。弾く対象の .mddata は GB 級になり得る。"""
    from metabolomix.tools import console_tools
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n" + "x" * 50000, encoding="ascii")

    def _forbidden(self):
        raise AssertionError("ファイル全体を read_bytes してはいけない")

    monkeypatch.setattr(Path, "read_bytes", _forbidden)
    assert console_tools._looks_like_method_text(method) is True


def test_console_status_returns_artifacts_as_tsv(tmp_path, monkeypatch):
    """生成物の全文は行が並ぶ一覧なので TSV（列名 1 回）で返す（include_artifacts=True 時）。"""
    import json as _json
    from metabolomix.console.job_manager import load_job
    from metabolomix.tools.console_tools import console_run, console_status

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    out = Path(load_job(job_path).run_dir) / "msdial"
    _use_console(monkeypatch, _fake_console({
        out / "Height_AlignResult-1.mzTab": _mztab_text(tmp_path, [tmp_path / "S1.wiff"]),
        out / "S1.mdpeak": "a",
        tmp_path / "S1_1.pai2": "x"}))
    console_run(str(job_path))

    parsed = _json.loads(console_status(str(job_path), include_artifacts=True))

    tsv = parsed["artifacts"]
    assert isinstance(tsv, str)
    lines = tsv.splitlines()
    assert lines[0] == "path\trole\tformat\troot"
    assert len(lines) == 3  # 列名 1 行 + 生成物 2 件
    assert "S1_1.pai2\tsample_peaks\tpai2\tdataset_root" in lines
    # run_dir 側の相対パスは -o の msdial/ を含む。root 列と合わせて出所が読める。
    assert any(l.endswith("S1.mdpeak\tsample_peak_table\tmdpeak\trun_dir") for l in lines)


def test_console_status_artifacts_tsv_is_empty_string_when_none(tmp_path, monkeypatch):
    """生成物ゼロなら空文字。列名だけの行を返して件数を誤読させない。"""
    import json as _json
    from metabolomix.tools.console_tools import console_status

    job_path = _planned_task8_job(tmp_path, monkeypatch)
    parsed = _json.loads(console_status(str(job_path), include_artifacts=True))
    assert parsed["artifacts"] == ""
