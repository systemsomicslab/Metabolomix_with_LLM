# tests/test_handoff_schema.py
import json
import pytest
from metabolomix.handoff.schema import (
    SCHEMA_VERSION,
    AnalysisJob,
    Artifact,
    MztabEntry,
    SampleManifest,
    sha256_file,
)


def _make_job(**kwargs) -> AnalysisJob:
    defaults = dict(
        schema=SCHEMA_VERSION,
        job_id="job_20260902_120000_pos_h",
        status="planned",
        created_at="2026-09-02T12:00:00+09:00",
        updated_at="2026-09-02T12:00:00+09:00",
        dataset_root="/data/study-001",
        input_count=56,
        software_name="MS-DIAL",
        software_version="5.5.260820",
        execution_mode="console",
        method_file="/data/params/lc_lipidomics.msdial",
        omics="lipidomics",
        polarity="positive",
        measure="peak_height",
        run_dir="/data/study-001/runs/job_20260902_120000_pos_h",
    )
    defaults.update(kwargs)
    return AnalysisJob(**defaults)


def test_round_trip(tmp_path):
    job = _make_job()
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    loaded = AnalysisJob.load(job_path)
    assert loaded.job_id == job.job_id
    assert loaded.status == "planned"
    assert loaded.polarity == "positive"
    assert loaded.measure == "peak_height"


def test_serialized_json_is_compact(tmp_path):
    job = _make_job()
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    text = job_path.read_text(encoding="utf-8")
    assert "\n" not in text, "JSON に改行が含まれています（compact 形式でない）"
    assert "  " not in text, "JSON にインデントが含まれています"


def test_schema_version_mismatch(tmp_path):
    job_path = tmp_path / "analysis-job.json"
    job_path.write_text(
        json.dumps({"schema": "analysis-job.v0", "job_id": "x"}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="schema mismatch"):
        AnalysisJob.load(job_path)


def test_mztab_entries_round_trip(tmp_path):
    job = _make_job(
        primary_mztab_files=[
            MztabEntry(
                path="mztab/neg-height.mzTab",
                polarity="negative",
                measure="peak_height",
                sha256="abc123",
                validation={"structure": "passed_with_warnings"},
            )
        ]
    )
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    loaded = AnalysisJob.load(job_path)
    assert len(loaded.primary_mztab_files) == 1
    e = loaded.primary_mztab_files[0]
    assert e.polarity == "negative"
    assert e.sha256 == "abc123"
    assert e.validation == {"structure": "passed_with_warnings"}


def test_artifacts_round_trip(tmp_path):
    job = _make_job(
        artifacts=[
            Artifact(path="msdial/AlignmentResult_001.arf", role="peak_matrix_source",
                     format="arf", sha256="def456"),
        ]
    )
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    loaded = AnalysisJob.load(job_path)
    assert loaded.artifacts[0].role == "peak_matrix_source"


def test_sample_manifest_round_trip(tmp_path):
    job = _make_job(
        sample_manifest=SampleManifest(
            path="sidecars/sample-manifest.tsv",
            sha256="aabbcc",
            status="approved",
        )
    )
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    loaded = AnalysisJob.load(job_path)
    assert loaded.sample_manifest is not None
    assert loaded.sample_manifest.status == "approved"


def test_no_sample_manifest(tmp_path):
    job = _make_job()
    job_path = tmp_path / "analysis-job.json"
    job.save(job_path)
    loaded = AnalysisJob.load(job_path)
    assert loaded.sample_manifest is None


def test_sha256_file(tmp_path):
    import hashlib
    data = b"hello world"
    f = tmp_path / "dummy.txt"
    f.write_bytes(data)
    h = sha256_file(f)
    assert len(h) == 64
    assert h == hashlib.sha256(data).hexdigest()


def test_status_transitions(tmp_path):
    job_path = tmp_path / "analysis-job.json"
    job = _make_job()
    job.save(job_path)

    job.status = "running"
    job.save(job_path)
    assert AnalysisJob.load(job_path).status == "running"

    job.status = "completed"
    job.save(job_path)
    assert AnalysisJob.load(job_path).status == "completed"


def test_schema_version_is_v2():
    assert SCHEMA_VERSION == "analysis-job.v2"


def test_artifact_root_defaults_to_run_dir():
    a = Artifact(path="msdial/x.mdpeak", role="sample_peak_table", format="mdpeak", sha256="ab")
    assert a.root == "run_dir"


def test_roundtrip_preserves_root_and_execution_fields(tmp_path):
    job = _make_job(
        status="completed",
        primary_mztab_files=[MztabEntry(path="msdial/A.mzTab", polarity="negative",
                                        measure="peak_height", sha256="c1")],
        artifacts=[Artifact(path="S1.pai2", role="sample_peaks", format="pai2",
                            sha256="c2", root="dataset_root")],
        save_project=True, timeout_s=21600,
    )
    p = tmp_path / "analysis-job.json"
    job.save(p)
    loaded = AnalysisJob.load(p)
    assert loaded.artifacts[0].root == "dataset_root"
    assert loaded.primary_mztab_files[0].root == "run_dir"
    assert loaded.save_project is True
    assert loaded.timeout_s == 21600


def test_load_accepts_v1_and_defaults_root(tmp_path):
    p = tmp_path / "analysis-job.json"
    p.write_text(json.dumps({
        "schema": "analysis-job.v1", "job_id": "old", "status": "completed",
        "created_at": "t", "updated_at": "t",
        "source": {"dataset_root": str(tmp_path), "input_count": 1},
        "software": {"name": "MS-DIAL", "version": "", "execution_mode": "console",
                     "method_file": "m.txt"},
        "project": {"omics": "lipidomics", "polarity": "negative", "measure": "peak_height"},
        "run_dir": str(tmp_path),
        "primary_mztab_files": [{"path": "A.mzTab", "polarity": "negative",
                                 "measure": "peak_height", "sha256": "x", "validation": {}}],
        "artifacts": [{"path": "S1.pai2", "role": "sample_peaks", "format": "pai2", "sha256": "y"}],
        "warnings": [], "error": None,
    }), encoding="utf-8")
    loaded = AnalysisJob.load(p)
    assert loaded.artifacts[0].root == "run_dir"
    assert loaded.primary_mztab_files[0].root == "run_dir"
    assert loaded.save_project is False
    assert loaded.timeout_s == 3600


def test_loaded_v1_job_is_saved_as_v2(tmp_path):
    p = tmp_path / "analysis-job.json"
    p.write_text(json.dumps({"schema": "analysis-job.v1", "job_id": "old",
        "status": "completed", "created_at": "t", "updated_at": "t",
        "source": {}, "software": {}, "project": {}, "run_dir": "r",
        "primary_mztab_files": [], "artifacts": [], "warnings": [], "error": None}), encoding="utf-8")
    loaded = AnalysisJob.load(p)
    loaded.save(p)
    assert json.loads(p.read_text(encoding="utf-8"))["schema"] == "analysis-job.v2"
