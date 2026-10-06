"""出所の記録 — 版数・極性の裏取り・サーバ版数。

実データの境界（analysis-job.v2）で残っていた穴:
  - `software.version` が空文字のまま。mzTab の `MTD software[1]` には
    `Msdial console 5.5.241113` が入っているのに拾っていなかった。
  - `polarity_source` が `job_declared` にしかならないため、宣言を間違えても
    検出できない（Console 出力のファイル名に極性トークンが無いのは仕様なので
    ここは変えられない。代わりにアダクトで裏取りする）。
  - サーバ版数がどの戻り値にも出ないため、古いプロセスが動いていても気づけない。
"""
from __future__ import annotations

import json as _json

from metabolomix.console.output_collector import (
    read_adduct_polarity,
    read_software_version,
)

_MTD = (
    "MTD\tmzTab-version\t2.0.0-M\n"
    "MTD\tsoftware[1]\t[MS, MS:1003082, MS-DIAL, Msdial console 5.5.241113]\n"
)


def _mztab(tmp_path, body=""):
    p = tmp_path / "AlignResult-1.mzTab"
    p.write_text(_MTD + body, encoding="utf-8")
    return p


# ---------- software version ----------

def test_reads_software_version_from_mtd(tmp_path):
    assert read_software_version(_mztab(tmp_path)) == "Msdial console 5.5.241113"


def test_software_version_is_none_when_absent(tmp_path):
    p = tmp_path / "x.mzTab"
    p.write_text("MTD\tmzTab-version\t2.0.0-M\n", encoding="utf-8")
    assert read_software_version(p) is None


def test_software_version_is_none_for_unreadable_file(tmp_path):
    assert read_software_version(tmp_path / "missing.mzTab") is None


# ---------- polarity crosscheck ----------

def _sml_rows(adducts):
    head = "SMH\tSML_ID\tSMF_ID_REFS\tdatabase_identifier\tadduct_ions\n"
    return head + "".join(f"SML\t{i}\t{i}\tnull\t{a}\n" for i, a in enumerate(adducts))


def test_adduct_majority_reads_positive(tmp_path):
    p = _mztab(tmp_path, _sml_rows(["[M+H]1+", "[M+NH4]1+", "[M-H]1-"]))
    result = read_adduct_polarity(p)
    assert result["adduct_majority"] == "positive"
    assert result["n_positive"] == 2
    assert result["n_negative"] == 1


def test_adduct_majority_reads_negative(tmp_path):
    p = _mztab(tmp_path, _sml_rows(["[M-H]1-", "[M+CH3COO]1-", "[M+H]1+"]))
    assert read_adduct_polarity(p)["adduct_majority"] == "negative"


def test_adduct_majority_is_none_without_usable_rows(tmp_path):
    p = _mztab(tmp_path, _sml_rows(["null", "null"]))
    assert read_adduct_polarity(p)["adduct_majority"] is None


def test_adduct_majority_handles_missing_column(tmp_path):
    p = _mztab(tmp_path, "SMH\tSML_ID\tdatabase_identifier\nSML\t0\tnull\n")
    assert read_adduct_polarity(p)["adduct_majority"] is None


# ---------- 結線 ----------

def test_console_run_records_software_version_from_mztab(tmp_path, monkeypatch):
    """Console 実行からは分からない版数を、成果物の mzTab から採る。"""
    from pathlib import Path
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    from tests.pipeline_fixtures import fake_console_command, use_fake_console

    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")
    out = Path(load_job(job_path).run_dir) / "msdial"
    use_fake_console(monkeypatch, fake_console_command({
        out / "AlignResult-1.mzTab": _MTD + _sml_rows(["[M-H]1-", "[M-H]1-"])}))

    console_run(str(job_path))

    assert load_job(job_path).software_version == "Msdial console 5.5.241113"


def test_console_run_crosschecks_polarity_without_overwriting_the_source(tmp_path, monkeypatch):
    """裏取りは別フィールドに置く。polarity_source は job_declared のまま嘘をつかない。

    宣言（negative）と中身（positive のアダクトばかり）が食い違うので、この実行は
    completed にならない。それでも出所と裏取りの記録は残す——「完了しなかった」と
    「何も分からない」は別のことなので、原因を読める形にしてから止める。
    """
    from pathlib import Path
    from metabolomix.console.job_manager import create_job, load_job
    from metabolomix.tools.console_tools import console_run
    from tests.pipeline_fixtures import fake_console_command, use_fake_console

    monkeypatch.setenv("MSDIAL_EXE", "fake.exe")
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")
    out = Path(load_job(job_path).run_dir) / "msdial"
    use_fake_console(monkeypatch, fake_console_command({
        out / "AlignResult-1.mzTab":
            _MTD + _sml_rows(["[M+H]1+", "[M+NH4]1+", "[M+Na]1+"])}))

    parsed = _json.loads(console_run(str(job_path)))
    job = load_job(job_path)

    entry = job.primary_mztab_files[0]
    assert entry.validation["polarity_source"] == "job_declared"
    assert entry.validation["polarity_crosscheck"]["adduct_majority"] == "positive"
    assert entry.validation["polarity_crosscheck"]["agrees"] is False
    assert "POLARITY_MISMATCH" in parsed["error"]["details"]["errors"]
    assert any("アダクト" in w for w in parsed["error"]["details"]["warnings"])


# ---------- サーバ版数 ----------

def test_console_status_reports_server_version(tmp_path, monkeypatch):
    from pathlib import Path
    from metabolomix.console.job_manager import create_job
    from metabolomix.tools.console_tools import console_status
    method = tmp_path / "params.txt"
    method.write_text("Ion mode: Negative\n", encoding="ascii")
    _, job_path = create_job(dataset_root=tmp_path, method_file=method,
                             polarity="negative", measure="peak_height")
    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["server_version"]


def test_dataset_status_reports_server_version(monkeypatch):
    """更新後に再起動を忘れると、古いプロセスが黙って動き続ける。"""
    import numpy as np
    from metabolomix.core import session_state
    from metabolomix.mztab.dataset_state import DatasetState
    from metabolomix.tools.mztab_tools import dataset_status
    ds = DatasetState()
    ds.feature_matrix = np.ones((2, 2))
    ds.sample_names = ["s1", "s2"]
    ds.feature_ids = ["1", "2"]
    session_state.session = session_state.AnalysisSession()
    session_state.session.dataset = ds
    parsed = _json.loads(dataset_status())
    assert parsed["server_version"]
