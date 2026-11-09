"""受領（`results/receipt.json`、`lipidmix.pipeline-receipt.v1`）と成果物の登録（plan 2026-11-09）。

pipeline は上流（Console → mzTab-M → 解析フォルダ）までを回し、最後の `receipt` 工程で
「極性ごとに、どの解析フォルダが・どの条件で・どれだけ確かめられて出来たか」を 1 つの
JSON に残す。続きは対話の側で `load_dataset(pipeline_path=…, polarity=…)` から。

- `build_receipt(record, pipeline_root)`: run の記録から受領を組み立てる（書かない）。
- `write_receipt(record, pipeline_root)`: 組み立てて書き、成果物として登録した ref を返す。
- `read_receipt(pipeline_root, record)`: 今の受領を、登録した hash と照合してから読む。
- `persist_result`: 書き終えた成果物（または、この場で書く数値結果）を ref にする。

研究室の参照ライブラリ（MSP・`[library] lbm` の LBM）は、ファイル名と sha256 と出どころ
だけを示す（run の記録の `inputs` は絶対パスのまま）。受領の中のパスは `pipeline_root`
相対（`next` の `<pipeline_path>` は `pipeline_status` が実際のパスで埋める）。

**sessionを一切importしない**（`tests/test_pipeline_engine.py` の AST の検査）。
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.analysis.sample_manifest import parse_manifest
from metabolomix.console import job_manager
from metabolomix.core.atomic_io import DomainError, atomic_write_json
from metabolomix.pipeline import inputs as inputs_mod
from metabolomix.pipeline import request as request_mod
from metabolomix.pipeline import store

__all__ = [
    "PIPELINE_PATH_PLACEHOLDER",
    "RECEIPT_OUTPUT_NAME",
    "RECEIPT_RELATIVE_PATH",
    "RECEIPT_SCHEMA",
    "build_receipt",
    "persist_result",
    "read_receipt",
    "write_receipt",
]

RECEIPT_SCHEMA = "lipidmix.pipeline-receipt.v1"
RECEIPT_RELATIVE_PATH = "results/receipt.json"
RECEIPT_OUTPUT_NAME = "receipt"
#: 受領の `next` の中で、run の場所を指す置き場（読む側が実際のパスで埋める）。
PIPELINE_PATH_PLACEHOLDER = "<pipeline_path>"

_INPUT_SUBDIR = "input"


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------- 成果物の登録 ----------

def _sanitize_nonfinite(value, *, reasons: list, path_prefix: str = ""):
    """非有限float（NaN/Infinity）をnullへ変換し、理由をreasonsへ記録する。

    `json.dump(..., allow_nan=False)`（`atomic_write_json`の既定）は非有限floatに
    出会うと例外で落ちる。呼ぶ前にここで潰す（pickleは使わない）。
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        reason = "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")
        reasons.append({"path": path_prefix or "$", "reason": reason})
        return None
    if isinstance(value, dict):
        return {k: _sanitize_nonfinite(v, reasons=reasons,
                                       path_prefix=f"{path_prefix}.{k}" if path_prefix else str(k))
                for k, v in value.items()}
    if isinstance(value, list):
        return [_sanitize_nonfinite(v, reasons=reasons, path_prefix=f"{path_prefix}[{i}]")
                for i, v in enumerate(value)]
    return value


def _result_path_for(pipeline_root: Path, relative: str, payload: dict) -> Path:
    """`relative` へ書く。既存と内容が違うなら内容 hash を足した別名を返す。

    既存が無い、または既存の中身が完全に同じなら `relative` をそのまま使う。
    別名は `<stem>__<先頭8桁>.json` で、同じ内容なら同じ名前に落ちる（決定論的）。
    """
    path = pipeline_root / relative
    # 直列化は atomic_write_json と同一にする（ensure_ascii=False・区切り詰め・
    # allow_nan=False）。ここがずれると、同じ内容でも別名に落ちて冪等でなくなる。
    encoded = json.dumps(payload, ensure_ascii=False,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    if not path.is_file():
        return path
    if path.read_bytes() == encoded:
        return path
    digest = hashlib.sha256(encoded).hexdigest()[:8]
    return path.with_name(f"{path.stem}__{digest}{path.suffix}")


def persist_result(pipeline_root, result: dict) -> dict:
    """成果物への参照(ref)を確定する。

    `result` は次のいずれかの形:

    - 既に書き終えたファイルを登録するだけ:
      ``{"output_name", "kind", "result_id", "path", parent_ids?, request_revision?}``。
      hashを計算し、pipeline_root相対へ書き換えるだけで、ファイルには書き込まない。
    - この場でJSONとして書く数値結果:
      ``{"output_name", "kind", "result_id", "data", relative_path?, parent_ids?,
      request_revision?}``。非有限値はnull+理由へ変換してから書く。

    戻り値は `output_name/result_id/kind/relative_path/hash/parent_ids/request_revision`
    を持つref（`store.verify_result_refs` がそのまま検証できる形）。
    """
    pipeline_root = Path(pipeline_root)
    output_name = result["output_name"]
    kind = result["kind"]
    result_id = result["result_id"]
    nonfinite_reasons: list = []

    if "path" in result:
        path = Path(result["path"])
        if not path.is_absolute():
            path = pipeline_root / path
        if not path.is_file():
            raise DomainError(
                "OUTPUT_ARTIFACT_MISSING",
                f"登録しようとした成果物ファイルがありません: {path}",
                {"path": str(path)})
    elif "data" in result:
        relative = result.get("relative_path") or f"results/{result_id}.json"
        sanitized = _sanitize_nonfinite(result["data"], reasons=nonfinite_reasons)
        payload = {
            "schema": "pipeline-result-data.v1",
            "kind": kind,
            "data": sanitized,
            "nonfinite_reasons": nonfinite_reasons,
        }
        # 既に記録された ref の (relative_path, hash) は検証可能なまま残す。内容が
        # 変わったときだけ別名へ逃がす——内容が同じ再実行ではファイルを増やさない。
        path = _result_path_for(pipeline_root, relative, payload)
        atomic_write_json(path, payload)
    else:
        raise DomainError(
            "OUTPUT_RESULT_INVALID",
            "resultには書き込み済みファイルのpath、またはこの場で書くdataの"
            "いずれかが必要です。",
            {"keys": sorted(result)})

    try:
        relative_path = str(path.resolve().relative_to(pipeline_root.resolve()))
    except ValueError as exc:
        raise DomainError(
            "OUTPUT_OUTSIDE_PIPELINE_ROOT",
            f"成果物がpipeline_root外にあります: {path}",
            {"path": str(path), "pipeline_root": str(pipeline_root)}) from exc
    relative_path = relative_path.replace(os.sep, "/")

    ref = {
        "output_name": output_name,
        "result_id": result_id,
        "kind": kind,
        "relative_path": relative_path,
        "hash": _sha256_file(path),
        "parent_ids": list(result.get("parent_ids") or []),
        "request_revision": result.get("request_revision"),
    }
    if nonfinite_reasons:
        ref["nonfinite_reasons"] = nonfinite_reasons
    return ref


# ---------- 受領の組み立て ----------

def _file_name(path: str | None) -> str | None:
    return Path(path).name if path else None


def _relative_to(path: Path, root: Path) -> str | None:
    try:
        return str(Path(path).resolve().relative_to(Path(root).resolve())).replace(os.sep, "/")
    except (OSError, ValueError):
        return None


def _load_request(pipeline_root: Path, record: dict) -> dict | None:
    saved = (record.get("request") or {}).get("saved_path")
    if not saved:
        return None
    try:
        return json.loads((Path(pipeline_root) / saved).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _library(inputs: dict) -> dict:
    """ライブラリの欄（研究室の資産はファイル名と sha256 と出どころだけ）。"""
    msp = inputs.get("msp") or {}
    lbm = inputs.get("lbm") or {}
    msp_out = None
    if msp.get("sha256"):
        msp_out = {"file": msp.get("file") or _file_name(msp.get("path")),
                   "sha256": msp.get("sha256"), "source": msp.get("source")}
    lbm_out = None
    if lbm.get("path") or lbm.get("sha256"):
        lbm_out = {"file": lbm.get("file") or _file_name(lbm.get("path")),
                   "sha256": lbm.get("sha256"), "source": lbm.get("source")}
    return {"mode": inputs.get("library_mode", "auto"), "msp": msp_out, "lbm": lbm_out}


def _job_or_none(job_path: Path | None):
    if job_path is None or not job_path.is_file():
        return None
    try:
        return job_manager.load_job(job_path)
    except Exception:  # noqa: BLE001 - 受領は読める範囲で書く（壊れたジョブで受領を止めない）
        return None


def _samples(pipeline_root: Path, record: dict, request: dict | None, polarity: str,
             unit_root: Path, warnings: list[dict]) -> dict:
    """試料の欄。今の実験情報シートを、Console が読んだ raw の一覧で解して載せる。"""
    input_dir = unit_root / _INPUT_SUBDIR
    raw_names = ([p.name for p in job_manager.list_raw_inputs(input_dir)]
                 if input_dir.is_dir() else [])
    out = {"manifest": None, "rows": None, "n_raw": len(raw_names)}
    if request is None:
        return out
    dataset_root = Path(record["identity"]["source_root"])
    try:
        index = request_mod.unit_index(request, polarity)
    except DomainError:
        return out
    unit_req = request_mod.unit_request(request, index, dataset_root)
    raw_root = request_mod.unit_source_root(dataset_root, request["units"][index])
    source = inputs_mod.manifest_source_record(raw_root, unit_req)
    if not source.get("path"):
        return out
    manifest_path = Path(source["path"])
    out["manifest"] = {"path": _relative_to(manifest_path, dataset_root) or str(manifest_path),
                       "sha256": source.get("sha256")}
    if not manifest_path.is_file():
        warnings.append({"code": "SAMPLE_MANIFEST_INVALID", "unit": polarity,
                         "message": f"実験情報シートが受領の時点で見つかりません: "
                                    f"{out['manifest']['path']}"})
        return out
    try:
        rows = parse_manifest(manifest_path, source_root=input_dir, expected_sources=raw_names)
    except DomainError as exc:
        warnings.append({"code": "SAMPLE_MANIFEST_INVALID", "unit": polarity,
                         "message": f"実験情報シートを受領の時点で読めませんでした（{exc.code}）: "
                                    f"{exc.message}"})
        return out
    for row in rows:
        if isinstance(row.get("source_file"), str):
            row["source_file"] = (_relative_to(Path(row["source_file"]), pipeline_root)
                                  or row["source_file"])
    out["rows"] = rows
    return out


def _unit_receipt(pipeline_root: Path, record: dict, request: dict | None,
                  polarity: str, unit: dict, warnings: list[dict]) -> dict:
    status = store.unit_status(record, polarity)
    unit_root_rel = unit.get("unit_root") or polarity
    unit_root = pipeline_root / unit_root_rel
    inputs = unit.get("inputs") or {}
    upstream = unit.get("upstream") or {}
    outputs = unit.get("outputs") or {}
    method = inputs.get("method") or {}
    exe = inputs.get("exe") or {}
    job_rel = upstream.get("console_job_path")
    job = _job_or_none(pipeline_root / job_rel if job_rel else None)
    verification = upstream.get("verification") or {}

    stages = record.get("stages") or {}
    stage_warnings: list[dict] = []
    error = None
    for handler in store.UNIT_HANDLERS:
        sid = store.unit_stage_id(handler, polarity)
        stage = stages.get(sid) or {}
        for w in stage.get("warnings") or []:
            stage_warnings.append({"code": w.get("code"), "stage_id": sid,
                                   "message": w.get("message")})
        if error is None and stage.get("status") in ("failed", "cancelled") and stage.get("error"):
            error = {"code": stage["error"].get("code"), "message": stage["error"].get("message"),
                     "stage_id": sid}

    return {
        "status": {"validated": "validated", "failed": "failed"}.get(status, "not_run"),
        "raw_dir": unit.get("raw_dir") or ".",
        "analysis_dir": f"{unit_root_rel}/{_INPUT_SUBDIR}",
        "console_job": job_rel,
        "mztab": outputs.get("mztab"),
        "alignment": outputs.get("alignment"),
        "method": {"file": _file_name(method.get("source_path")),
                   "sha256": method.get("sha256"),
                   "effective_sha256": method.get("effective_sha256")},
        "library": _library(inputs),
        "console": {"exe_file": _file_name(exe.get("path")), "exe_sha256": exe.get("sha256"),
                    "software_version": (getattr(job, "software_version", None) or None)},
        "alignment_reference": inputs.get("alignment_reference"),
        "samples": _samples(pipeline_root, record, request, polarity, unit_root, warnings),
        "validation": {"job_status": verification.get("status"),
                       "termination": verification.get("termination"),
                       "exit_code": verification.get("exit_code"),
                       "gui_project": outputs.get("gui_project"),
                       "warnings": stage_warnings},
        "error": error,
    }


def _next_call(validated: list[str]) -> str | None:
    if not validated:
        return None
    choice = " または ".join(f'"{p}"' for p in validated)
    return f'load_dataset(pipeline_path="{PIPELINE_PATH_PLACEHOLDER}", polarity={choice})'


def build_receipt(record: dict, pipeline_root) -> dict:
    """run の記録から受領を組み立てる（ファイルには書かない）。

    単位の順は run の単位の順（要求の順）。`warnings` は受領の時点で見つかった問題
    （事前検査の後に書き換えられて読めなくなった実験情報シートなど）。
    """
    pipeline_root = Path(pipeline_root)
    request = _load_request(pipeline_root, record)
    warnings: list[dict] = []
    units = {polarity: _unit_receipt(pipeline_root, record, request, polarity, unit, warnings)
             for polarity, unit in store.run_units(record).items()}
    validated = [p for p, u in units.items() if u["status"] == "validated"]
    return {
        "schema": RECEIPT_SCHEMA,
        "pipeline_id": (record.get("identity") or {}).get("pipeline_id"),
        "created_at": _now_iso(),
        "request_revision": (record.get("request") or {}).get("revision"),
        "units": units,
        "warnings": warnings,
        "next": _next_call(validated),
    }


def write_receipt(record: dict, pipeline_root) -> tuple[dict, dict]:
    """受領を組み立てて `results/receipt.json` へ書き、(受領, 登録した ref) を返す。

    同じ場所へ上書きする（受領は run の今の姿を 1 つだけ指す）。前の ref は
    `store.current_result_refs` が最後の 1 件に畳むので、上書きで壊れたとは読まれない。
    """
    pipeline_root = Path(pipeline_root)
    receipt = build_receipt(record, pipeline_root)
    path = pipeline_root / RECEIPT_RELATIVE_PATH
    atomic_write_json(path, receipt)
    ref = persist_result(pipeline_root, {
        "output_name": RECEIPT_OUTPUT_NAME, "kind": "receipt",
        "result_id": f"res_receipt_{receipt['pipeline_id']}", "path": str(path),
        "request_revision": receipt["request_revision"],
    })
    return receipt, ref


def read_receipt(pipeline_root, record: dict) -> dict | None:
    """今の受領を読む。登録が無い・hash が合わない・読めないなら None。"""
    pipeline_root = Path(pipeline_root)
    refs = [r for r in store.current_result_refs(record)
            if r.get("output_name") == RECEIPT_OUTPUT_NAME]
    if not refs:
        return None
    ref = refs[-1]
    if store.verify_result_refs(pipeline_root, [ref]):
        return None
    try:
        return json.loads((pipeline_root / ref["relative_path"]).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
