"""データセットが同定に使った参照ライブラリを覚える。

後段の MS/MS 照合（Task 10）が、利用者にライブラリを名指しさせずに同じものを引けるようにする。
記録（LibraryRecord）は dict で、`path`（絶対パス）はデータセットの内側だけで持つ。
**戻り値に出すときは `public_view` を通す**——研究室ライブラリの置き場所は LLM の文脈に出さない。

status: `known`（1 本に決まった）/ `unknown`（決められない）/ `mismatch`（記録した sha256 と
今のファイルが違う、またはファイルが無い）/ `lbm_only`（LBM だけで MSP は使っていない）。
"""
from __future__ import annotations

import re
from pathlib import Path

from metabolomix.console.method_file import LBM_KEY, read_method_keys
from metabolomix.core.atomic_io import DomainError
from metabolomix.handoff.schema import sha256_file

# 解析フォルダの `*_Loaded.<kind>.dbs` のうち LBM のもの（`.lbm2.dbs` など）。
_LBM_DBS_RE = re.compile(r"\.lbm\d*\.dbs$", re.IGNORECASE)
# Console は LBM の注釈器キーに LBM ファイルのパスを使う（`.lbm` + 数字で終わる）。
_LBM_KEY_RE = re.compile(r"\.lbm\d*$", re.IGNORECASE)


def _record(**fields) -> dict:
    base = {"status": "unknown", "kind": None, "file": None, "path": None, "sha256": None,
            "origin": None, "candidates": [], "lbm_also": False}
    base.update(fields)
    return base


def unknown_record() -> dict:
    """mzTab-M を直接読んだとき（同定に使ったライブラリは mzTab-M から分からない）。"""
    return _record()


def _measure(path: Path) -> str | None:
    try:
        return sha256_file(path)
    except OSError:
        return None


def remember_from_folder(selection) -> dict:
    """解析フォルダの `*_Loaded.*.dbs`（MS-DIAL が同定に使った実体）から覚える。"""
    lbm = [p for p in selection.loaded_libraries if _LBM_DBS_RE.search(p.name)]
    other = [p for p in selection.loaded_libraries if not _LBM_DBS_RE.search(p.name)]
    if len(other) == 1:
        path = other[0]
        return _record(status="known", kind="dbs", file=path.name, path=str(path.resolve()),
                       sha256=_measure(path), origin="loaded_dbs", lbm_also=bool(lbm))
    if len(other) >= 2:
        return _record(candidates=[p.name for p in other])
    if lbm:
        path = lbm[0]
        return _record(status="lbm_only", kind="lbm", file=path.name, path=str(path.resolve()),
                       sha256=_measure(path), origin="loaded_dbs")
    return _record()


def remember_from_job(job_path: Path) -> dict:
    """ジョブの実効メソッドが宣言するライブラリから覚える（sha256 は今測る）。"""
    from metabolomix.handoff.schema import AnalysisJob
    from metabolomix.mztab.loading import effective_method_msp

    job = AnalysisJob.load(Path(job_path))
    msp = effective_method_msp(job)
    if msp is not None:
        return _record(status="known", kind="msp", file=msp.name, path=str(msp),
                       sha256=_measure(msp), origin="effective_method",
                       lbm_also=_declared_lbm(job) is not None)
    lbm = _declared_lbm(job)
    if lbm is not None:
        return _record(status="lbm_only", kind="lbm", file=lbm.name, path=str(lbm),
                       sha256=_measure(lbm), origin="effective_method")
    return _record()


def _declared_lbm(job) -> Path | None:
    method_text = (getattr(job, "method_file", "") or "").strip()
    if not method_text:
        return None
    method = Path(method_text)
    try:
        keys = read_method_keys(method)
    except (OSError, DomainError):
        return None
    declared = (keys.get(LBM_KEY.lower()) or "").strip()
    if not declared:
        return None
    lbm = Path(declared)
    return lbm if lbm.is_absolute() else method.parent / lbm


def _abs_in_run(value, pipeline_root: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else pipeline_root / path


def remember_from_run(pipeline_root: Path) -> dict:
    """パイプラインの run が固定した入力（`inputs.msp` / `inputs.lbm`）から覚える。

    run が固定した sha256 と今のファイルを突き合わせ、違う・無いなら `mismatch`。
    """
    from metabolomix.pipeline.store import load_run

    pipeline_root = Path(pipeline_root)
    inputs = load_run(pipeline_root).get("inputs") or {}
    msp = inputs.get("msp") or {}
    lbm = inputs.get("lbm") or {}
    mode = inputs.get("library_mode")
    if msp.get("source") != "not_used" and msp.get("path"):
        path = _abs_in_run(msp["path"], pipeline_root)
        rec = _record(status="known", kind="msp", file=msp.get("file") or path.name,
                      path=str(path), sha256=msp.get("sha256"), origin="run_inputs",
                      lbm_also=mode == "msp_then_lbm")
    elif lbm.get("path"):
        path = _abs_in_run(lbm["path"], pipeline_root)
        rec = _record(status="lbm_only", kind="lbm", file=lbm.get("file") or path.name,
                      path=str(path), sha256=lbm.get("sha256"), origin="run_inputs")
    else:
        return _record()
    return verify_record(rec)


def verify_record(record: dict) -> dict:
    """今のファイルの sha256 を測り直す。違う・ファイルが無いなら `mismatch` の写しを返す。

    記録に sha256 が無いとき（ファイルが読めなかった記録）は、ファイルの有無だけを見る。
    元の記録は書き換えない。
    """
    result = dict(record)
    if record.get("status") not in ("known", "lbm_only") or not record.get("path"):
        return result
    current = _measure(Path(record["path"]))
    pinned = record.get("sha256")
    if current is None or (pinned is not None and current != pinned):
        result["status"] = "mismatch"
    return result


def public_view(record: dict | None) -> dict | None:
    """戻り値に出す写し。置き場所（`path`）を外す。"""
    if record is None:
        return None
    return {k: v for k, v in record.items() if k != "path"}


def is_lbm_annotator(annotator_id: str | None) -> bool | None:
    """代表の `annotator_id` が LBM 由来なら True、判定できなければ None（False は返さない）。

    Console は LBM の注釈器キーに LBM ファイルのパス（末尾が `.lbm` + 数字）を使い、MSP の
    注釈器キーは利用者の付けた名前になる（実データで確認。docs/output_format/library.md）。
    GUI の注釈器 ID は自由名なので、規則に合わないことは「LBM でない」の証拠にならない。
    注釈器キーは研究室ライブラリの置き場所を含みうるので、この関数の呼び出し側も ID を出さない。
    """
    if not annotator_id:
        return None
    return True if _LBM_KEY_RE.search(annotator_id.strip()) else None
