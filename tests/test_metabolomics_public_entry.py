"""公開入口から起こした v2 run を、実 handler 一式で通す。

`MetabolomicsHarness` は `run_engine` を直接回し、上流4工程を test double へ
差し替える（工程単体の試験）。ここはそれと目的が違う——**受付を公開入口に
通し、上流も実 handler で回す**。「工程は通るが入口から起動できない」という
状態を検出できるのはこちらだけで、実際それが 2026-09-16 の監査で見つかった
穴だった。

Console の実バイナリだけは起動できないので、既存の Console 層テストと同じ
唯一の注入口（`console_runner.build_msdial_cmd`）を偽コマンドへ差し替える。
受付・入力配置・supervise 監視・出力収集・mzTab 読込・以降の全工程は本物。
"""
from __future__ import annotations

import json
from pathlib import Path

from metabolomix.pipeline import service, store, worker
from tests.metabolomics_fixtures import (
    DEFAULT_STATISTICS, INJECTIONS, write_arf, write_manifest_v2, write_mztab_v2,
    write_profile,
)
from tests.pipeline_fixtures import fake_console_command, use_fake_console


def _source(tmp_path, monkeypatch) -> tuple[Path, dict, list[Path]]:
    monkeypatch.setattr("metabolomix.console.runner.is_console_exe", lambda *a, **k: True)
    source_root = tmp_path / "source"
    source_root.mkdir()
    sources = []
    for row in INJECTIONS:
        path = source_root / f"{row['sample_id']}.wiff"
        path.write_bytes(b"synthetic raw")
        sources.append(path)
    manifest = write_manifest_v2(source_root / "sample-manifest.tsv", sources)
    profile = write_profile(source_root / "profile.json", statistics=DEFAULT_STATISTICS)
    request = {
        "schema": "pipeline-request.v2",
        "profile_file": str(profile),
        "execution_purpose": "validation",
        "sample_manifest": str(manifest),
        "statistics": DEFAULT_STATISTICS,
    }
    return source_root, request, sources


def _console_outputs(tmp_path, pipeline_path: Path) -> tuple[dict, dict]:
    """偽 Console が書く成果物（staged 入力に対応する mzTab と ARF）。

    `stage_inputs` は元の名前を保って `pipeline_root/input` へ配置するので、
    配置前でも staged 絶対パスを先に計算できる（v1 の E2E と同じ規約）。
    ARF は MessagePack なのでバイト列のまま渡す——テキストとして書くと壊れ、
    注入証拠が読めないまま `load_assay_evidence` が成功してしまう。
    """
    staged = [pipeline_path / "input" / f"{row['sample_id']}.wiff" for row in INJECTIONS]
    scratch = tmp_path / "_console_template"
    scratch.mkdir(exist_ok=True)
    mztab_text = write_mztab_v2(scratch / "Height.mzTab", staged).read_text(encoding="utf-8")
    arf_path = write_arf(scratch / "AlignmentResult_PeakProperties.arf", staged)

    run_dir = pipeline_path / "console" / "attempt-0001"
    text_files = {run_dir / "msdial" / "Height_AlignmentResult_1.mzTab": mztab_text}
    binary_files = {run_dir / "msdial" / "AlignmentResult_PeakProperties.arf":
                    arf_path.read_bytes()}
    return text_files, binary_files


def _fake_console(monkeypatch, tmp_path, pipeline_path: Path) -> None:
    text_files, binary_files = _console_outputs(tmp_path, pipeline_path)
    use_fake_console(monkeypatch,
                     fake_console_command(text_files, binary=binary_files))


def test_the_public_entry_creates_a_v2_run(tmp_path, monkeypatch):
    source_root, request, _sources = _source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source_root, request)

    record = store.load_run(Path(receipt["pipeline_path"]))
    assert record["schema"] == "pipeline-run.v2"
    assert "execute_console" in record["stages"]
    assert record["inputs"]["polarity"]["source"] == "profile"


def test_the_public_entry_runs_every_stage_with_the_real_handlers(tmp_path, monkeypatch):
    source_root, request, _sources = _source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source_root, request)
    pipeline_path = Path(receipt["pipeline_path"])
    _fake_console(monkeypatch, tmp_path, pipeline_path)

    result = worker.run_worker(pipeline_path)

    unfinished = {sid: stage["status"] for sid, stage in result["stages"].items()
                  if stage["status"] != "succeeded"}
    assert not unfinished, unfinished
    assert result["status"] == "completed", result
    # 上流4工程が test double でないことを、この試験だけが示す。
    assert {"prepare_inputs", "execute_console", "validate_outputs",
            "load_dataset"} <= set(result["stages"])

    produced = {ref["output_name"] for ref in (result.get("results") or [])}
    assert {"profile", "execution_manifest", "sample_manifest", "assay_evidence",
            "feature_bindings", "matrix", "qc_population", "qc", "feature_table",
            "quality_report"} <= produced, sorted(produced)


def test_the_started_job_is_v3_and_carries_the_profile_snapshot(tmp_path, monkeypatch):
    source_root, request, _sources = _source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source_root, request)
    pipeline_path = Path(receipt["pipeline_path"])
    _fake_console(monkeypatch, tmp_path, pipeline_path)
    worker.run_worker(pipeline_path)

    job_path = next(pipeline_path.rglob("analysis-job.json"))
    data = json.loads(job_path.read_text(encoding="utf-8"))
    assert data["schema"] == "analysis-job.v3"
    assert data["profile"]["profile_id"] == "synthetic-metabolomics"
    assert data["project"]["omics"] == "metabolomics"


def test_a_console_that_produces_nothing_does_not_complete_the_run(tmp_path, monkeypatch):
    """完走の判定が exit 0 ではなく成果物に基づくことを示す。

    これが無いと上の完走試験は空虚になる——「偽 Console が何を書いても
    completed になる」経路でも通ってしまう。
    """
    source_root, request, _sources = _source(tmp_path, monkeypatch)

    receipt = service.plan_pipeline(source_root, request)
    pipeline_path = Path(receipt["pipeline_path"])
    use_fake_console(monkeypatch, fake_console_command({}))  # exit 0 だが何も書かない

    result = worker.run_worker(pipeline_path)

    assert result["status"] != "completed"
    assert result["stages"]["execute_console"]["status"] != "succeeded" \
        or result["stages"]["validate_outputs"]["status"] != "succeeded"


def test_a_completed_run_can_be_opened_for_interactive_work(tmp_path, monkeypatch):
    """公開入口で回した run を、そのまま対話セッションへ載せられること。

    S1（dataset_build_matrix）・S2（dataset_load の pipeline_path モード）・
    S3（この公開入口）が1本に繋がっていることを、この試験だけが示す。
    """
    from metabolomix.core import session_state
    from metabolomix.tools.dataset_analysis_tools import dataset_statistic
    from metabolomix.tools.mztab_tools import dataset_load

    source_root, request, _sources = _source(tmp_path, monkeypatch)
    receipt = service.plan_pipeline(source_root, request)
    pipeline_path = Path(receipt["pipeline_path"])
    _fake_console(monkeypatch, tmp_path, pipeline_path)
    assert worker.run_worker(pipeline_path)["status"] == "completed"

    session_state.session.dataset = None
    summary = dataset_load(pipeline_path=str(pipeline_path))
    assert "解析行列" in summary, summary

    ds = session_state.session.dataset
    assert ds.analysis_matrices
    assert ds.sample_metadata_rows

    matrix_id = sorted(ds.analysis_matrices)[0]
    result = json.loads(dataset_statistic(
        {"statistic_id": "check", "kind": "pca", "matrix_recipe_id": "default",
         "transform": "none", "feature_scope": {"mode": "all_eligible"},
         "scaling": "autoscale", "n_components": 2},
        matrix_id))
    assert result.get("status") == "success", result
