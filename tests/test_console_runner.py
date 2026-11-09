# tests/test_console_runner.py
#
# Console の実行は `metabolomix.console.execution.supervise` が唯一の経路で、
# pipeline の worker もここへ入る。
# ここは runner / output_collector / job_manager の単体テスト。Console の実際の起動と監視は
# test_console_supervisor.py と pipeline のテストが見る。
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
from tests.config_fixtures import use_config


def _use_exe(tmp_path, monkeypatch, name: str) -> Path:
    """`[msdial] exe` に置き場所だけの偽の実行体（`tmp_path/_app/<name>`）を書く。

    設定ファイルの値は実在を確かめる（`get_exe_path`）ので、ファイルは作る。中身は
    実行できないので、起動まで進むテストは probe か build_msdial_cmd を差し替える。
    """
    exe = tmp_path / "_app" / name
    exe.parent.mkdir(parents=True, exist_ok=True)
    exe.write_text("fake console executable placeholder", encoding="ascii")
    use_config(monkeypatch, tmp_path, {"msdial.exe": exe})
    return exe


# ---------- runner ----------

def test_get_exe_path_missing(monkeypatch):
    with pytest.raises(MsdialExeNotFoundError):
        get_exe_path()


def test_get_exe_path_set(tmp_path, monkeypatch):
    exe = _use_exe(tmp_path, monkeypatch, "MsdialConsoleApp.exe")
    assert get_exe_path() == str(exe)


def test_get_exe_path_rejects_a_configured_path_that_does_not_exist(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"msdial.exe": tmp_path / "gone" / "MSDIALCUI.exe"})
    with pytest.raises(MsdialExeNotFoundError) as exc:
        get_exe_path()
    assert "[msdial] exe" in str(exc.value)


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


def test_probe_returns_timeout_after_retry():
    """--help が時間切れを繰り返せば "timeout"（NOT_CONSOLE と区別する）。1 回やり直す。"""
    from metabolomix.console.runner import probe_console_exe
    with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("x", 60)) as run:
        assert probe_console_exe("slow.exe") == "timeout"
    assert run.call_count == 2
    assert run.call_args.kwargs["timeout"] == 60


def test_probe_retry_succeeds():
    """1 回目は冷間起動で時間切れ、2 回目に lcms を含む出力が返れば "console"。"""
    from metabolomix.console.runner import probe_console_exe
    completed = MagicMock()
    completed.stdout = "MSDIAL Console Application\n  lcms   Run LC-MS data processing\n"
    with patch("subprocess.run",
               side_effect=[subprocess.TimeoutExpired("x", 60), completed]) as run:
        assert probe_console_exe("cold.exe") == "console"
    assert run.call_count == 2


def test_probe_not_console_does_not_retry():
    """出力が返って lcms が無いのは確定判定。やり直さない。"""
    from metabolomix.console.runner import probe_console_exe
    completed = MagicMock()
    completed.stdout = "GUI banner\n"
    with patch("subprocess.run", return_value=completed) as run:
        assert probe_console_exe("gui.exe") == "not_console"
    assert run.call_count == 1


def test_probe_os_error_is_not_console_unless_raised():
    from metabolomix.console.runner import probe_console_exe
    with patch("subprocess.run", side_effect=OSError("not found")):
        assert probe_console_exe("nope.exe") == "not_console"
        with pytest.raises(OSError):
            probe_console_exe("nope.exe", raise_on_os_error=True)


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



# ---------- Task 0: ブロッカー回帰 ----------

def test_create_job_allows_repo_data_dir(tmp_path, monkeypatch):
    """<repo>/data 配下（BASE_DIR / "data"）を dataset_root にできる。"""
    from metabolomix.core import mcp_core
    from metabolomix.console.job_manager import create_job
    # リポジトリ直下の data/ を模す: BASE_DIR 配下だが BASE_DIR / "data" 配下でもある
    fake_repo = tmp_path / "repo"
    data_root = fake_repo / "data" / "study-001"
    data_root.mkdir(parents=True)
    method = tmp_path / "params.msdial"
    method.touch()
    monkeypatch.setattr(mcp_core, "BASE_DIR", fake_repo)

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


# ---------- F1: 実行後処理（collect_artifacts 以降）のガード ----------

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


def test_method_text_check_reads_only_the_head(tmp_path, monkeypatch):
    """判定は先頭だけを読む。弾く対象の .mddata は GB 級になり得る。"""
    from metabolomix.tools import console_tools
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n" + "x" * 50000, encoding="ascii")

    def _forbidden(self):
        raise AssertionError("ファイル全体を read_bytes してはいけない")

    monkeypatch.setattr(Path, "read_bytes", _forbidden)
    assert console_tools._looks_like_method_text(method) is True
