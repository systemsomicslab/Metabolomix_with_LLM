"""戻り値の肥大を防ぐ縛り。

`console_status` は生成物を 1 サンプルにつき 5 件出す（.pai2 / .dcl / _tags.xml /
.mdpeak / .mdmsp）。60 サンプルの実走で 313 件になり、全文 TSV は数万字に達して
後ろの warnings / error を埋没させた（このテストのフィクスチャは同じ形を
`_job_with_artifacts` の既定 60 サンプルで再現し、300 件になる）。件数と
役割別内訳が既定で、全文は明示要求時のみ。
"""
from __future__ import annotations

import json as _json

from metabolomix.console.job_manager import create_job, load_job, save_job
from metabolomix.handoff.schema import Artifact
from metabolomix.tools.console_tools import console_status


def _job_with_artifacts(tmp_path, n_samples=60):
    method = tmp_path / "param.txt"
    method.write_text("Ion mode: Positive\n", encoding="ascii")
    job, job_path = create_job(dataset_root=tmp_path, method_file=method,
                               polarity="positive", measure="peak_height")
    roles = ("sample_peaks", "msms_evidence", "peak_tags",
             "sample_peak_table", "msms_spectra")
    job.artifacts = [
        Artifact(path=f"s{i}.{role}", role=role, format=role, sha256="")
        for i in range(n_samples) for role in roles
    ]
    job.status = "completed"
    save_job(job, job_path)
    return job_path


def test_console_status_omits_the_artifact_table_by_default(tmp_path):
    job_path = _job_with_artifacts(tmp_path)
    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["artifact_count"] == 300
    assert "artifacts" not in parsed


def test_console_status_reports_counts_per_role_instead(tmp_path):
    job_path = _job_with_artifacts(tmp_path)
    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["artifacts_by_role"] == {
        "msms_evidence": 60, "msms_spectra": 60, "peak_tags": 60,
        "sample_peak_table": 60, "sample_peaks": 60,
    }


def test_console_status_returns_the_full_table_on_request(tmp_path):
    job_path = _job_with_artifacts(tmp_path)
    parsed = _json.loads(console_status(str(job_path), include_artifacts=True))
    lines = parsed["artifacts"].splitlines()
    assert lines[0] == "path\trole\tformat\troot"
    assert len(lines) == 301


def test_console_status_keeps_reporting_zero_artifacts(tmp_path):
    """生成物ゼロは「まだ何も出ていない」で、要約が消えてよい理由にはならない。"""
    job_path = _job_with_artifacts(tmp_path, n_samples=0)
    parsed = _json.loads(console_status(str(job_path)))
    assert parsed["artifact_count"] == 0
    assert parsed["artifacts_by_role"] == {}
