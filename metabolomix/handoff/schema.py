"""analysis-job.json のスキーマ定義と読み書き。

AnalysisJob はジョブの永続スナップショット。MCP サーバ再起動後の再開・GUI 表示・
hash 検証に使う。ファイルはデータフォルダ内のランディレクトリに置き、
リポジトリには書き込まない（metabolomix/console/job_manager.py が保証する）。

このモジュールは依存グラフの leaf（stdlib のみ）。tools_* / server を import しない。
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from metabolomix.core.atomic_io import atomic_write_json

# v2: 生成物のルートと実行オプションを明示する。
SCHEMA_VERSION = "analysis-job.v2"

# v3: 検証済みプロファイルで実行したジョブ（spec §6, §6.1）。profile snapshot
# ——プロファイル同一性・method 原本と実効コピーの hash・依存ファイル一覧・
# 実行環境 manifest——を丸ごと持つ。
#
# 書き出す版は**ジョブの内容で決める**。profile snapshot を持つジョブだけが v3 で、
# 持たないジョブは今までどおり v2 で書く（v1 経路の出力は 1 バイトも変わらない）。
SCHEMA_VERSION_V3 = "analysis-job.v3"

# 読み込みは旧版も受け続ける（spec §6「readerは旧版を継続して受ける」）。
SUPPORTED_SCHEMA_VERSIONS = frozenset(
    {"analysis-job.v1", SCHEMA_VERSION, SCHEMA_VERSION_V3})

ArtifactRoot = Literal["run_dir", "dataset_root"]

JobStatus = Literal["planned", "running", "needs_input", "partial", "completed", "failed", "cleaned"]
Polarity = Literal["positive", "negative"]
MeasureType = Literal["peak_height", "peak_area_above_zero"]
OmicsType = Literal["lipidomics", "metabolomics"]


@dataclass
class MztabEntry:
    path: str
    polarity: Polarity
    measure: MeasureType
    sha256: str
    validation: dict = field(default_factory=dict)
    root: ArtifactRoot = "run_dir"


@dataclass
class Artifact:
    path: str
    role: str
    format: str
    sha256: str
    root: ArtifactRoot = "run_dir"


@dataclass
class SampleManifest:
    path: str
    sha256: str
    status: Literal["pending", "approved", "auto_generated"] = "pending"


@dataclass
class AnalysisJob:
    schema: str
    job_id: str
    status: JobStatus
    created_at: str
    updated_at: str
    dataset_root: str
    input_count: int
    software_name: str
    software_version: str
    execution_mode: str
    method_file: str
    omics: OmicsType
    polarity: Polarity
    measure: MeasureType
    run_dir: str
    primary_mztab_files: list[MztabEntry] = field(default_factory=list)
    artifacts: list[Artifact] = field(default_factory=list)
    sample_manifest: SampleManifest | None = None
    warnings: list[str] = field(default_factory=list)
    error: str | None = None
    save_project: bool = False
    timeout_s: int = 3600
    #: v3 のみ。metabolomix.console.profiles.snapshot_profile() の戻り値そのまま
    #: （profile 同一性・method・依存一覧・実行環境 manifest・実効メソッドの hash）。
    #: None なら v1/v2 のジョブで、書き出しも従来どおり v2。
    profile_snapshot: dict | None = None

    @property
    def dependencies(self) -> list[dict]:
        """依存ファイル manifest（method_key / kind / path / sha256）。

        snapshot の中の 1 箇所だけを正準にする——同じ一覧をジョブ側にも複製すると、
        どちらかだけ更新された記録が生まれて「hash は合うのに中身が違う」になる。
        """
        return list((self.profile_snapshot or {}).get("dependencies") or [])

    @property
    def environment(self) -> dict:
        """実行環境 manifest（adapter 版・実行体 hash 等）。出所は同じく snapshot。"""
        return dict((self.profile_snapshot or {}).get("execution_environment") or {})

    def save(self, path: Path) -> None:
        """analysis-job.json を原子的に保存する（内容に応じて v2 / v3）。

        書込途中でプロセスが落ちても、既存ファイルは壊れた内容に置き換わらない
        （metabolomix.core.atomic_io.atomic_write_json）。
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_json(path, _to_dict(self))

    @staticmethod
    def load(path: Path) -> "AnalysisJob":
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("schema") not in SUPPORTED_SCHEMA_VERSIONS:
            raise ValueError(
                f"analysis-job schema mismatch: expected one of "
                f"{sorted(SUPPORTED_SCHEMA_VERSIONS)}, got {data.get('schema')!r}"
            )
        return _from_dict(data)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------- シリアライズ補助 ----------

def _schema_version_for(job: AnalysisJob) -> str:
    """書き出す schema を**ジョブの内容から**決める（writer 側の schema dispatch）。

    profile snapshot を持つジョブだけが v3。持たないジョブは v2 のまま——読み込んだ
    v1 のジョブを保存し直すと v2 になる従来の挙動も、そのまま変わらない。
    """
    return SCHEMA_VERSION_V3 if job.profile_snapshot is not None else SCHEMA_VERSION


def _to_dict(job: AnalysisJob) -> dict:
    data = {
        "schema": _schema_version_for(job),
        "job_id": job.job_id,
        "status": job.status,
        "created_at": job.created_at,
        "updated_at": job.updated_at,
        "source": {
            "dataset_root": job.dataset_root,
            "input_count": job.input_count,
        },
        "software": {
            "name": job.software_name,
            "version": job.software_version,
            "execution_mode": job.execution_mode,
            "method_file": job.method_file,
        },
        "project": {
            "omics": job.omics,
            "polarity": job.polarity,
            "measure": job.measure,
        },
        "run_dir": job.run_dir,
        "primary_mztab_files": [
            {
                "path": e.path,
                "polarity": e.polarity,
                "measure": e.measure,
                "sha256": e.sha256,
                "validation": e.validation,
                "root": e.root,
            }
            for e in job.primary_mztab_files
        ],
        "artifacts": [
            {"path": a.path, "role": a.role, "format": a.format,
             "sha256": a.sha256, "root": a.root}
            for a in job.artifacts
        ],
        "sample_manifest": (
            {
                "path": job.sample_manifest.path,
                "sha256": job.sample_manifest.sha256,
                "status": job.sample_manifest.status,
            }
            if job.sample_manifest
            else None
        ),
        "warnings": job.warnings,
        "error": job.error,
        "execution": {
            "save_project": job.save_project,
            "timeout_s": job.timeout_s,
        },
    }
    if job.profile_snapshot is not None:
        data["profile"] = job.profile_snapshot
    return data


def _from_dict(d: dict) -> AnalysisJob:
    src = d.get("source", {})
    sw = d.get("software", {})
    proj = d.get("project", {})

    mztab = [
        MztabEntry(
            path=e["path"],
            polarity=e["polarity"],
            measure=e["measure"],
            sha256=e.get("sha256", ""),
            validation=e.get("validation", {}),
            root=e.get("root", "run_dir"),
        )
        for e in d.get("primary_mztab_files", [])
    ]
    artifacts = [
        Artifact(path=a["path"], role=a["role"], format=a["format"],
                 sha256=a.get("sha256", ""), root=a.get("root", "run_dir"))
        for a in d.get("artifacts", [])
    ]
    sm_raw = d.get("sample_manifest")
    sample_manifest = (
        SampleManifest(
            path=sm_raw["path"],
            sha256=sm_raw.get("sha256", ""),
            status=sm_raw.get("status", "pending"),
        )
        if sm_raw
        else None
    )
    execution = d.get("execution", {})
    return AnalysisJob(
        schema=d["schema"],
        job_id=d["job_id"],
        status=d["status"],
        created_at=d["created_at"],
        updated_at=d["updated_at"],
        dataset_root=src.get("dataset_root", ""),
        input_count=src.get("input_count", 0),
        software_name=sw.get("name", "MS-DIAL"),
        software_version=sw.get("version", ""),
        execution_mode=sw.get("execution_mode", "console"),
        method_file=sw.get("method_file", ""),
        omics=proj.get("omics", "lipidomics"),
        polarity=proj.get("polarity", "positive"),
        measure=proj.get("measure", "peak_height"),
        run_dir=d.get("run_dir", ""),
        primary_mztab_files=mztab,
        artifacts=artifacts,
        sample_manifest=sample_manifest,
        warnings=d.get("warnings", []),
        error=d.get("error"),
        save_project=execution.get("save_project", False),
        timeout_s=execution.get("timeout_s", 3600),
        profile_snapshot=d.get("profile"),
    )
