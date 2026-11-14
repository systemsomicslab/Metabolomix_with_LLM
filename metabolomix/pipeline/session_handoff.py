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

__all__ = ["UnitRun", "bind_manifest_rows", "load_run_dataset", "resolve_unit",
           "run_manifest_rows"]

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
    （旧 v1 の run が持つ差次的結果は戻さない）。試料の対応表の行は試料へ結んでから当てる
    （`bind_manifest_rows`。位置では当てない。結べなければ当てず、注記を `ds.load_notes` に残す）。
    """
    unit = resolve_unit(pipeline_path, polarity, allow_incomplete=allow_incomplete)
    ds = mztab_loading.load_dataset_state(job_path=unit.job_path,
                                          allow_incomplete=allow_incomplete)
    if unit.manifest_rows is not None:
        bind_manifest_rows(ds, unit.manifest_rows)
    return ds


def _rows_by_name(rows: list[dict], sample_names: list[str]) -> list[dict | None]:
    """試料名ごとの行。行の `source_file`（raw）の名前（拡張子を除く）か `sample_id` で結ぶ
    （解析フォルダの経路 `dataset.opening._apply_run_rows` と同じ鍵）。"""
    by_key: dict[str, dict] = {}
    for row in rows:
        source = row.get("source_file")
        for key in (Path(str(source)).stem if source else None, row.get("sample_id")):
            if key:
                by_key.setdefault(str(key).strip().casefold(), row)
    return [by_key.get(str(name or "").strip().casefold()) for name in sample_names]


def bind_manifest_rows(ds, rows: list[dict]) -> str | None:
    """run の試料の対応表の行（受領、旧 v1 は記録）をデータセットの試料へ結んで当てる。

    `pipeline_path` で開く両方のアダプタ（mzTab-M は `load_run_dataset`、解析フォルダは
    `dataset.opening._folder_via_run`）が通る。当てられなければ当てず、注記を `ds.load_notes` に
    残して返す（`load_dataset` の `summary` の「読み込みの注記」と `dataset_status` の `load_notes` に出る。
    読み込みは止めない）。当てられれば None。

    行の並びはシートの並び（`pipeline_sample_sheet` なら取得順）で、`ds.sample_names`（mzTab-M の
    assay の順＝MS-DIAL のファイルの順）とは限らない。試料の表を位置で読む消費側（統計エンジンの
    試料の表・`included_sample_indices`）に順の違う行を渡すと、群と include=false が別の試料に付く。
    そこでまず raw のパス（mzTab-M の `ms_run[N]-location` → `ds.assay_sources`。解析フォルダのアダプタは
    `_folder_via_job` が同じ値を `ds.sample_raw_paths` に写す。`sample_manifest.resolve_metadata`）で結び、
    結べなければ（raw を動かした run など）行の raw の名前か `sample_id` を試料名と結ぶ（run の中の
    解析フォルダを `directory` で開く経路 `dataset.opening._apply_run_rows` と同じ鍵）。どちらでも全試料を
    結べなければ当てない（群は `dataset_set_sample_metadata` で与える）。読み込みの一部なので
    `metadata_revision` は進めない。
    """
    note = _bind_manifest_rows(ds, rows)
    if note:
        ds.load_notes.append(note)
    return note


def _bind_manifest_rows(ds, rows: list[dict]) -> str | None:
    """`bind_manifest_rows` の本体。当てられなければ注記を返す。"""
    from metabolomix.analysis import sample_manifest

    def resolve():
        try:
            return sample_manifest.resolve_metadata(ds, rows)
        except (DomainError, TypeError):    # TypeError: source_file の無い行（raw のパスで結べない）
            return None

    previous = getattr(ds, "sample_raw_paths", None)
    resolved = resolve()
    if resolved is None:
        matched = _rows_by_name(rows, list(ds.sample_names))
        missing = [name for name, row in zip(ds.sample_names, matched) if row is None]
        if missing:
            return (f"試料 {missing[:5]} を run の試料の対応表と raw のパスでも名前でも結べないので、"
                    "対応表は当てていない（群は dataset_set_sample_metadata で与える）。")
        ds.sample_raw_paths = [row.get("source_file") for row in matched]
        resolved = resolve()
        if resolved is None:
            ds.sample_raw_paths = previous
            return ("run の試料の対応表の行に raw のパスが無く試料へ結べないので、対応表は"
                    "当てていない（群は dataset_set_sample_metadata で与える）。")
    revision = ds.metadata_revision
    try:
        sample_manifest.apply_defaults_and_inputs(ds, resolved)
    except DomainError as exc:
        ds.sample_raw_paths = previous
        return (f"run の試料の対応表を当てられなかった（{exc.code}: {exc.message}）。"
                "群は dataset_set_sample_metadata で与える。")
    ds.metadata_revision = revision
    return None


def run_manifest_rows(pipeline_root: Path, record: dict) -> list[dict] | None:
    """旧 v1 の run が解決した試料対応表の行を返す。どこにも無ければ None。

    旧 v1 の run の `resolve_metadata` が `record["inputs"]["manifest"]` にだけ書いた行
    （`resolve_sample_metadata` の行）を返す。これを載せないと、開いたデータセットの群が空になる。
    """
    rows = (record.get("inputs") or {}).get("manifest")
    if isinstance(rows, list) and rows:
        return [dict(row) for row in rows if isinstance(row, dict)]
    return None
