"""session に依存しない mzTab 読み込みと、出所の信用度（spec §7）。

読み込みが答えるべきことは 2 つある。**読めたか**と、**どれだけ信用してよいか**。
旧実装は前者しか答えず、`status == "completed"` という文字列だけを後者の代わりに
していた。旧経路の completed は「実行後にファイルが増えた」以上の意味を持たず、
途中で落ちた実行にも付いた——つまり信用の根拠になっていなかった。

ここで縛るのは:

- 完了していない実行の出力は**既定で読まない**。明示して読んだ場合は探索専用と
  記録し、2 群比較と差次的エクスポートを拒否する。
- `verified` は終了証跡（exit code・identity）とファイル hash が揃って初めて名乗る。
- 読み込みは `session_state` も `mcp_core` も触らない（ワーカーから呼べる）。
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from metabolomix.console.job_manager import create_job, save_job
from metabolomix.core.atomic_io import DomainError
from metabolomix.handoff.schema import Artifact, MztabEntry, sha256_file
from metabolomix.mztab.loading import load_dataset_state
from tests.pipeline_fixtures import execution_record, write_mztab


def _job_with_outputs(tmp_path, *, status="completed", raw_count=2,
                      with_receipt=True, receipt_overrides=None,
                      corrupt_hash=False):
    """合成 raw と主 mzTab を持つジョブを、実行せずに組み立てる。

    実 MS-DIAL も実行成果物も使わない。証跡（execution-result.json）と
    監視サイドカー（worker.json）も、この関数が意図した内容で書く。
    """
    source = tmp_path / "source"
    source.mkdir(exist_ok=True)
    raws = [source / f"S{i}.abf" for i in range(1, raw_count + 1)]
    for i, raw in enumerate(raws, start=1):
        raw.write_bytes(b"\x00" * (32 + i))
    method = source / "params.txt"
    method.write_text("Ion mode: Positive\n", encoding="ascii")

    job, job_path = create_job(source, method, "positive", "peak_height",
                               input_count=len(raws))
    run_dir = Path(job.run_dir)
    mztab = run_dir / "msdial" / "Height_AlignmentResult_1.mzTab"
    mztab.parent.mkdir(parents=True, exist_ok=True)
    write_mztab(mztab, raws)

    recorded_hash = sha256_file(mztab)
    if corrupt_hash:
        recorded_hash = "0" * 64
    job.primary_mztab_files = [MztabEntry(
        path="msdial/Height_AlignmentResult_1.mzTab", polarity="positive",
        measure="peak_height", sha256=recorded_hash, root="run_dir")]
    job.artifacts = [Artifact(path="S1_1.pai2", role="sample_peaks",
                              format="pai2", sha256="x", root="dataset_root")]
    job.status = status
    save_job(job, job_path)

    from metabolomix.console.execution import (
        receipt_path, supervision_state_path, write_supervision_inputs)
    write_supervision_inputs(run_dir, {
        "raw_inventory": [str(p.resolve()) for p in raws]})
    assert supervision_state_path(run_dir).is_file()
    if with_receipt:
        record = execution_record(job_id=job.job_id, **(receipt_overrides or {}))
        receipt_path(run_dir).write_text(json.dumps(record), encoding="utf-8")
    return job_path, raws


# ---------- 引数 ----------

def test_requires_exactly_one_source(tmp_path):
    with pytest.raises(DomainError) as exc:
        load_dataset_state()
    assert exc.value.code == "DATASET_BAD_REQUEST"
    with pytest.raises(DomainError):
        load_dataset_state(mztab_path=tmp_path / "a.mzTab",
                           job_path=tmp_path / "analysis-job.json")


def test_allow_incomplete_must_be_a_real_bool(tmp_path):
    """"false" は真だが False ではない。型で弾く。"""
    with pytest.raises(DomainError) as exc:
        load_dataset_state(job_path=tmp_path / "analysis-job.json",
                           allow_incomplete="false")
    assert exc.value.code == "DATASET_BAD_REQUEST"


# ---------- 直接読み ----------

def test_direct_read_is_unverified(tmp_path):
    """直接読んだファイルは、どの実行から出たかを語れない。"""
    raw = tmp_path / "S1.abf"
    raw.write_bytes(b"\x00")
    path = write_mztab(tmp_path / "direct.mzTab", [raw])

    ds = load_dataset_state(mztab_path=path)

    assert ds.source_verification == "direct_unverified"
    assert ds.exploratory_only is False
    assert len(ds.sample_names) == 1


def test_direct_read_of_a_missing_file(tmp_path):
    with pytest.raises(DomainError) as exc:
        load_dataset_state(mztab_path=tmp_path / "nope.mzTab")
    assert exc.value.code == "MZTAB_NOT_FOUND"


# ---------- ジョブ経由 ----------

def test_completed_job_with_a_matching_receipt_is_verified(tmp_path):
    job_path, raws = _job_with_outputs(tmp_path)

    ds = load_dataset_state(job_path=job_path)

    assert ds.source_verification == "verified"
    assert ds.exploratory_only is False
    assert ds.job_path == str(job_path)
    assert len(ds.sample_names) == len(raws)
    assert ds.artifact_paths["sample_peaks"]


def test_assay_sources_map_to_the_planned_raw_files(tmp_path):
    """どの列がどの生データから来たかは、完了ゲートと同じ対応で持つ。"""
    job_path, raws = _job_with_outputs(tmp_path)

    ds = load_dataset_state(job_path=job_path)

    assert sorted(ds.assay_sources.values()) == sorted(str(r.resolve()) for r in raws)


@pytest.mark.parametrize("status", ["partial", "failed", "running"])
def test_incomplete_jobs_are_refused_by_default(tmp_path, status):
    job_path, _ = _job_with_outputs(tmp_path, status=status)

    with pytest.raises(DomainError) as exc:
        load_dataset_state(job_path=job_path)
    assert exc.value.code == "INCOMPLETE_ANALYSIS_JOB"
    assert exc.value.details["status"] == status


def test_partial_requires_explicit_exploratory_mode(tmp_path):
    job_path, _ = _job_with_outputs(tmp_path, status="partial")

    with pytest.raises(DomainError, match="INCOMPLETE_ANALYSIS_JOB"):
        load_dataset_state(job_path=job_path)

    ds = load_dataset_state(job_path=job_path, allow_incomplete=True)
    assert ds.exploratory_only is True
    assert "partial" in ds.validation_result["warnings"][0]


def test_a_completed_job_without_a_receipt_is_legacy_unverified(tmp_path):
    """completed という文字列は裏取りではない。旧実行はそう名乗らせる。"""
    job_path, _ = _job_with_outputs(tmp_path, with_receipt=False)

    ds = load_dataset_state(job_path=job_path)

    assert ds.source_verification == "legacy_unverified"
    assert any("裏取り" in w for w in ds.validation_result["warnings"])


def test_a_nonzero_exit_receipt_is_not_verified(tmp_path):
    job_path, _ = _job_with_outputs(
        tmp_path, receipt_overrides={"exit_code": 1})

    assert load_dataset_state(job_path=job_path).source_verification == "legacy_unverified"


def test_a_receipt_from_another_job_is_not_verified(tmp_path):
    job_path, _ = _job_with_outputs(tmp_path)
    from metabolomix.console.execution import receipt_path
    from metabolomix.console.job_manager import load_job
    run_dir = Path(load_job(job_path).run_dir)
    record = json.loads(receipt_path(run_dir).read_text(encoding="utf-8"))
    record["job_id"] = "someone_elses_job"
    receipt_path(run_dir).write_text(json.dumps(record), encoding="utf-8")

    assert load_dataset_state(job_path=job_path).source_verification == "legacy_unverified"


def test_a_mztab_that_no_longer_matches_its_recorded_hash_is_not_verified(tmp_path):
    """記録した hash と中身が違うなら、それは記録された実行の出力ではない。"""
    job_path, _ = _job_with_outputs(tmp_path, corrupt_hash=True)

    ds = load_dataset_state(job_path=job_path)

    assert ds.source_verification == "legacy_unverified"
    assert any("hash" in w for w in ds.validation_result["warnings"])


def test_job_without_primary_mztab_files(tmp_path):
    job_path, _ = _job_with_outputs(tmp_path)
    from metabolomix.console.job_manager import load_job, save_job as _save
    job = load_job(job_path)
    job.primary_mztab_files = []
    _save(job, job_path)

    with pytest.raises(DomainError) as exc:
        load_dataset_state(job_path=job_path)
    assert exc.value.code == "MZTAB_NOT_FOUND"


def test_legacy_v1_job_is_loadable(tmp_path):
    """analysis-job.v1 の completed も読める（読込互換は維持する）。"""
    job_path, raws = _job_with_outputs(tmp_path, with_receipt=False)
    data = json.loads(job_path.read_text(encoding="utf-8"))
    data["schema"] = "analysis-job.v1"
    for entry in data["primary_mztab_files"]:
        entry.pop("root", None)
    job_path.write_text(json.dumps(data), encoding="utf-8")

    ds = load_dataset_state(job_path=job_path)

    assert len(ds.sample_names) == len(raws)
    assert ds.source_verification == "legacy_unverified"


# ---------- 層の境界 ----------

def test_loading_does_not_import_session_or_mcp_core():
    """ワーカーから呼べる読み込みであること。グローバル状態を触らせない。"""
    import metabolomix.mztab.loading as loading
    tree = ast.parse(Path(loading.__file__).read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not any(m.startswith("metabolomix.core.session_state")
                   or m.startswith("metabolomix.core.mcp_core")
                   or m.startswith("metabolomix.tools") for m in imported), imported
