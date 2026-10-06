"""AnalysisJob の CRUD とランディレクトリ管理。

ランディレクトリはデータフォルダ内にのみ作成する（リポジトリ外を強制）。
"""
from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.handoff.schema import (
    SCHEMA_VERSION,
    AnalysisJob,
    MeasureType,
    OmicsType,
    Polarity,
)

JOB_FILENAME = "analysis-job.json"
RUNS_SUBDIR = "runs"

#: pipelineがConsole jobを所有していることを記す小さなsidecar（Task 14）。
#: `metabolomix/console/output_collector.py`の`_OPERATIONAL_FILES`に登録済みなので
#: 生成物としては収集されない。中身は`{"pipeline_path": "<絶対str>"}`のみ
#: （所有の可否は毎回そのpipeline-run.jsonの現在状態を動的に見て判定するため、
# sidecar自体には状態を持たせない）。
PIPELINE_OWNER_FILENAME = "pipeline-owner.json"

# MS-DIAL の SupportMsRawDataExtension と同じ集合
# （src/MSDIAL5/MsdialCore/Enum/SupportFormat.cs）。
_RAW_EXTENSIONS = frozenset({
    "abf", "ibf", "cdf", "mzml", "wiff", "raw", "d", "wiff2", "qgd", "lcd", "lrp", "imzml",
})

#: フォルダそのものが 1 検体になれる拡張子。Agilent / Bruker の `.d` と
#: Waters の `.raw` がこれに当たる（Thermo の `.raw` は同じ拡張子のファイル）。
#: 上流 `AnalysisFilesParser.ReadFolderContents` の `isVendorDirectory`
#: （`Directory.Exists(path) && (extension == ".raw" || extension == ".d")`）と
#: 同じ規則。`_RAW_EXTENSIONS` の残りは `DataAccess.IsDataFormatSupported` が
#: `File.Exists` を要求するため、同名のフォルダがあっても入力にならない。
_VENDOR_DIR_EXTENSIONS = frozenset({"d", "raw"})


def is_raw_input(entry: Path) -> bool:
    """`entry` が MS-DIAL の計測データ 1 件として読まれるかを返す。

    拡張子だけで判定すると、名前が `backup.mzml` のフォルダまで 1 検体として
    数えてしまう。それは `input_count` と `raw_inventory` に実在しない
    サンプルを混ぜ、`ms_run[N]-location` との 1 対 1 照合が実行の最後に
    なって落ちる形で表面化する。
    """
    ext = entry.suffix.lower().lstrip(".")
    if ext not in _RAW_EXTENSIONS:
        return False
    if entry.is_dir():
        return ext in _VENDOR_DIR_EXTENSIONS
    return True


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _job_id(polarity: str, measure: str) -> str:
    """人が読める接頭辞に UUID を足したジョブ ID を作る。

    秒精度のタイムスタンプだけでは、同じ秒に 2 件計画すると ID が衝突して
    **同じランディレクトリを 2 つのジョブが共有する**（後から計画したほうが
    先のジョブの analysis-job.json を上書きし、実行中の証跡も混ざる）。
    人が一覧で読める部分は残したまま、一意性は UUID 側に持たせる。
    """
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    m_short = "h" if measure == "peak_height" else "a"
    p_short = polarity[:3]
    return f"job_{ts}_{p_short}_{m_short}_{uuid.uuid4().hex[:8]}"


def create_job(
    dataset_root: Path,
    method_file: Path,
    polarity: Polarity,
    measure: MeasureType,
    omics: OmicsType = "lipidomics",
    software_version: str = "",
    input_count: int = 0,
) -> tuple[AnalysisJob, Path]:
    """新規ジョブを作成してランディレクトリと analysis-job.json を返す。

    ランディレクトリが BASE_DIR（このリポジトリ）の配下にないことを確認する。
    """
    _assert_not_in_repo(dataset_root)

    job_id = _job_id(polarity, measure)
    run_dir = dataset_root / RUNS_SUBDIR / job_id
    run_dir.mkdir(parents=True, exist_ok=True)

    now = _now_iso()
    job = AnalysisJob(
        schema=SCHEMA_VERSION,
        job_id=job_id,
        status="planned",
        created_at=now,
        updated_at=now,
        dataset_root=str(dataset_root),
        input_count=input_count,
        software_name="MS-DIAL",
        software_version=software_version,
        execution_mode="console",
        method_file=str(method_file),
        omics=omics,
        polarity=polarity,
        measure=measure,
        run_dir=str(run_dir),
    )
    job_path = run_dir / JOB_FILENAME
    job.save(job_path)
    return job, job_path


def load_job(job_path: Path) -> AnalysisJob:
    if not job_path.is_file():
        raise FileNotFoundError(f"analysis-job.json が見つかりません: {job_path}")
    return AnalysisJob.load(job_path)


def update_status(job_path: Path, status: str, error: str | None = None) -> AnalysisJob:
    job = load_job(job_path)
    job.status = status  # type: ignore[assignment]
    job.updated_at = _now_iso()
    if error is not None:
        job.error = error
    job.save(job_path)
    return job


def save_job(job: AnalysisJob, job_path: Path) -> None:
    job.updated_at = _now_iso()
    job.save(job_path)


def list_jobs(dataset_root: Path) -> list[Path]:
    """dataset_root/runs/ 以下の analysis-job.json を新しい順に返す。"""
    runs_dir = dataset_root / RUNS_SUBDIR
    if not runs_dir.is_dir():
        return []
    found = sorted(
        runs_dir.rglob(JOB_FILENAME),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return found


def list_raw_inputs(dataset_root: Path) -> list[Path]:
    """データフォルダ直下の計測ファイルを安定順で返す。

    「Console が実際に読む入力はどれか」の唯一の答え。監視入力の目録
    （`write_supervision_inputs` の `raw_inventory`）はこの一覧で固定し、
    実行後に mzTab の `ms_run[N]-location` と 1 対 1 で突き合わせる。
    数えるだけの `raw_input_summary` もここを通す（数と中身が食い違わない）。
    """
    # is_file() で絞らない。Agilent の `.d` と Waters の `.raw` はフォルダ
    # そのものが 1 検体の計測データで、MS-DIAL もフォルダを入力として受ける。
    # どの拡張子でフォルダが認められるかは `is_raw_input` が上流と揃える。
    entries = [entry for entry in Path(dataset_root).iterdir() if is_raw_input(entry)]
    return sorted(entries, key=lambda p: str(p).lower())


def raw_input_summary(dataset_root: Path) -> dict[str, int]:
    """データフォルダ直下の計測ファイルを拡張子ごとに数える。"""
    counts: dict[str, int] = {}
    for entry in list_raw_inputs(dataset_root):
        ext = entry.suffix.lower().lstrip(".")
        counts[ext] = counts.get(ext, 0) + 1
    return counts


def count_raw_inputs(dataset_root: Path) -> int:
    """データフォルダ内の計測ファイル数を返す（raw_input_summary の合計）。"""
    return sum(raw_input_summary(dataset_root).values())


def _assert_not_in_repo(path: Path) -> None:
    """path が「リポジトリ内かつデータディレクトリ外」でないことを確認する。

    ランディレクトリを版管理下のソースツリーに掘らせないためのガード。
    ただし既定のデータディレクトリは <repo>/data（data_config.DEFAULT_DATA_DIR）で
    あり、これはリポジトリ配下だが .gitignore 済みの運用領域なので許可する。
    「リポジトリ配下すべて禁止」にすると既定構成が丸ごと使えなくなる。
    """
    from metabolomix.core import mcp_core
    from metabolomix.core.data_config import get_data_dir

    target = path.resolve()
    repo_root = mcp_core.BASE_DIR.resolve()
    if not _is_relative_to(target, repo_root):
        return  # リポジトリ外。何も言わない
    data_dir = get_data_dir().resolve()
    if _is_relative_to(target, data_dir):
        return  # <repo>/data 配下は運用領域として許可
    raise ValueError(
        f"ランディレクトリをリポジトリのソースツリー内 ({repo_root}) に作成しようとしました。"
        f"dataset_root はリポジトリ外か、データディレクトリ ({data_dir}) 配下を"
        f"指定してください: {path}"
    )


def pipeline_owner_path(run_dir: Path) -> Path:
    """このrun_dir（Console jobのランディレクトリ）向けの所有権sidecarのパス。"""
    return Path(run_dir) / PIPELINE_OWNER_FILENAME


def read_pipeline_owner(run_dir: Path) -> dict | None:
    """所有権sidecarを読む。

    「ファイルが無い」（＝一度も所有登録されていない）と「ファイルはあるが
    読めない／壊れている」（＝判定不能）を区別する。前者は None（呼び出し側は
    「所有記録なし」として続行してよい）。後者は空dict（`pipeline_path`を
    持たない）を返し、呼び出し側（`pipeline_owner_block_reason`）に
    「判定不能なので拒否」を選ばせる——記録が消えた／壊れただけで単体
    console_run/console_cleanup の保護が抜けるのは安全側の設計として誤り。
    """
    path = pipeline_owner_path(run_dir)
    try:
        raw = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    except OSError:
        return {}
    try:
        data = json.loads(raw)
    except ValueError:
        return {}
    return data if isinstance(data, dict) else {}


def _is_relative_to(child: Path, parent: Path) -> bool:
    """child が parent 配下かを bool で返す（例外を制御フローに使わない）。"""
    try:
        child.relative_to(parent)
    except ValueError:
        return False
    return True
