"""完了した run の成果物を、対話セッションが使える DatasetState へ戻す層。

`load_dataset(pipeline_path=..., polarity=...)`（`dataset.opening`）が `resolve_unit` で
極性の単位を 1 つ選び、その単位の解析フォルダ・Console のジョブ・試料の対応表を引く。
v3 の run は最後の工程が書いた受領（`receipt.read_receipt`）から読み、受領が無い・照合できない
ときだけ run の記録の単位から読む。旧 v1 の run は読むだけの記録（`store.is_legacy`）で、
1 単位として今までどおり開く（試料の対応表は載せるが、run の差次的結果は戻さない）。

グローバル session は import しない（この層は DatasetState を作って返すだけで、
どこへ置くかは呼び出し側＝MCP 公開層が決める）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import receipt as receipt_mod
from metabolomix.pipeline import recovery, store

__all__ = ["UnitRun", "load_run_dataset", "resolve_unit", "run_manifest_rows"]

_INPUT_SUBDIR = "input"


@dataclass(frozen=True)
class UnitRun:
    """run の極性の単位 1 つを開くのに要るもの（パスはすべて絶対）。"""

    pipeline_root: Path
    record: dict
    polarity: str
    legacy: bool
    unit_root: Path
    #: MS-DIAL が `.arf` 等を書いた解析フォルダ（受領の `analysis_dir`。無ければ `<単位>/input`）。
    analysis_dir: Path
    job_path: Path
    #: 受領の時点で解いた試料の対応表の行（無ければ None）。
    manifest_rows: list[dict] | None
    notes: list[str] = field(default_factory=list)


def _polarity_options(record: dict, units: dict) -> list[dict]:
    return [{"value": p, "label": f"{p}（{store.unit_status(record, p)}）"}
            for p in units if store.unit_status(record, p) == "validated"]


def resolve_unit(pipeline_path, polarity: str | None = None, *,
                 allow_incomplete: bool = False) -> UnitRun:
    """run の極性の単位を 1 つ選び、その解析フォルダ・ジョブ・試料の対応表を返す。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい。`polarity` を省いて
    単位が 1 つならそれ、2 つ以上なら `PIPELINE_POLARITY_REQUIRED`（`details.options` が選べる
    極性）。v3 の単位は検証を通ったもの（`validated`）だけ開き、それ以外は
    `PIPELINE_UNIT_NOT_VALIDATED`（`allow_incomplete` なら探索用に読む）。
    """
    root = recovery.normalize_pipeline_root(Path(pipeline_path).expanduser())
    record = store.load_run(root)
    legacy = store.is_legacy(record)
    units = store.run_units(record)
    if not units:
        raise DomainError(
            "PIPELINE_UPSTREAM_INCOMPLETE",
            "この run はまだ上流（Console）の成果物を持っていません。"
            "pipeline_status で状態を確認してください。", {"pipeline_root": str(root)})

    if polarity is None:
        if len(units) > 1:
            options = _polarity_options(record, units)
            if not options and not allow_incomplete:
                raise _not_validated(root, record, units, None)
            raise DomainError(
                "PIPELINE_POLARITY_REQUIRED",
                "この run には極性の単位が複数あります。どれを開くか polarity で選んでください。",
                {"options": options or [{"value": p, "label": p} for p in units]})
        polarity = next(iter(units))
    elif polarity not in units:
        raise DomainError(
            "PIPELINE_POLARITY_NOT_FOUND",
            f"この run に極性 {polarity!r} の単位はありません（{', '.join(units)}）。",
            {"polarity": polarity, "available": list(units)})

    if not legacy and store.unit_status(record, polarity) != "validated" and not allow_incomplete:
        raise _not_validated(root, record, units, polarity)

    unit = units[polarity]
    receipt = None if legacy else receipt_mod.read_receipt(root, record)
    unit_receipt = ((receipt or {}).get("units") or {}).get(polarity)
    notes: list[str] = []
    if not legacy and unit_receipt is None:
        notes.append("run の受領（results/receipt.json）が無い・照合できないので、run の記録の"
                     "単位から解析フォルダを引いた（試料の対応表は載らない）。群を使うなら"
                     "dataset_set_sample_metadata で実験情報シートを載せ直すか、pipeline_status で"
                     "受領を確かめてから開き直す。")

    unit_root = root / (unit.get("unit_root") or ".")
    job_rel = ((unit_receipt or {}).get("console_job")
               or (unit.get("upstream") or {}).get("console_job_path"))
    if not job_rel:
        raise DomainError(
            "PIPELINE_UPSTREAM_INCOMPLETE",
            "この run はまだ上流（Console）の成果物を持っていません。"
            "pipeline_status で状態を確認してください。",
            {"pipeline_root": str(root), "polarity": polarity})
    if unit_receipt is not None and unit_receipt.get("analysis_dir"):
        analysis_dir = root / unit_receipt["analysis_dir"]
    else:
        analysis_dir = unit_root / _INPUT_SUBDIR

    if legacy:
        rows = run_manifest_rows(root, record)
    else:
        rows = _receipt_rows(root, unit_receipt)
    return UnitRun(pipeline_root=root, record=record, polarity=polarity, legacy=legacy,
                   unit_root=unit_root, analysis_dir=analysis_dir, job_path=root / job_rel,
                   manifest_rows=rows, notes=notes)


def _not_validated(root: Path, record: dict, units: dict, polarity: str | None) -> DomainError:
    statuses = {p: store.unit_status(record, p) for p in units}
    which = f"極性 {polarity!r} の単位は" if polarity else "どの単位も"
    return DomainError(
        "PIPELINE_UNIT_NOT_VALIDATED",
        f"{which}検証を通っていません（{statuses}）。pipeline_status で確かめてください。"
        "中断・失敗した実行の出力を探索で読むなら allow_incomplete=True。",
        {"pipeline_root": str(root), "polarity": polarity, "units": statuses,
         "validated": [p for p, s in statuses.items() if s == "validated"]})


def _receipt_rows(root: Path, unit_receipt: dict | None) -> list[dict] | None:
    """受領の単位の試料の対応表の行。`source_file` は pipeline_root 相対なので絶対にして返す。"""
    rows = ((unit_receipt or {}).get("samples") or {}).get("rows")
    if not isinstance(rows, list) or not rows:
        return None
    out = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        row = dict(row)
        source = row.get("source_file")
        if isinstance(source, str) and source and not Path(source).is_absolute():
            row["source_file"] = str(root / source)
        out.append(row)
    return out or None


def load_run_dataset(pipeline_path, *, polarity: str | None = None,
                     allow_incomplete: bool = False):
    """run の成果物から DatasetState を組み立てて返す（保存も session 更新もしない）。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい
    （`pipeline_status` と同じ規約）。選んだ単位のジョブと試料の対応表だけを載せる
    （旧 v1 の run が持つ差次的結果は戻さない）。
    """
    unit = resolve_unit(pipeline_path, polarity, allow_incomplete=allow_incomplete)
    ds = mztab_loading.load_dataset_state(job_path=unit.job_path,
                                          allow_incomplete=allow_incomplete)
    if unit.manifest_rows is not None:
        ds.sample_metadata_rows = unit.manifest_rows
    return ds


def run_manifest_rows(pipeline_root: Path, record: dict) -> list[dict] | None:
    """旧 v1 の run が解決した試料対応表の行を返す。どこにも無ければ None。

    旧 v1 の run の `resolve_metadata` が `record["inputs"]["manifest"]` にだけ書いた行
    （`resolve_sample_metadata` の行）を返す。これを載せないと、開いたデータセットの群が空になる。
    """
    rows = (record.get("inputs") or {}).get("manifest")
    if isinstance(rows, list) and rows:
        return [dict(row) for row in rows if isinstance(row, dict)]
    return None
