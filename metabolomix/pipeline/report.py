"""必須成果物の確定判定と決定的な品質レポート（spec §7.4, §11, brief task-17）。

4つの公開関数で「runは何を達成すれば完了と言えるか」「今どこまで達成したか」
「それを人が読める形でどう記録するか」を分ける:

- `required_outputs`: 要求（target・save_project・comparisons）から、目標別の
  必須出力名の一覧を機械的に組み立てる（spec §9.2）。
- `persist_result`: 書き終えた成果物（PNG/TSV）を登録する、または数値結果を
  この場でJSONへ書く（非有限値はnull+理由へ変換。pickle不使用）。戻り値は
  `metabolomix.pipeline.store.verify_result_refs`がそのまま検証できるref。
- `evaluate_target`: `record["results"]`の各refを`required_outputs`と照合し、
  `completed`/`partial`/`failed`/`needs_input`を判定する。report生成前の
  判定にも、report保存後の最終確定にも同じ関数を使う——「レポートが自分の
  未生成を理由に永久失敗する」循環を、呼び出しタイミングの違いだけで解く
  （関数自身は分岐を持たない）。
- `write_pipeline_report`: 決定的なMarkdownテンプレートへ落とす。LLM生成文・
  自動生物学的解釈には一切依存しない。
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from metabolomix.core.atomic_io import DomainError, atomic_write_json
from metabolomix.pipeline import store

REPORT_SCHEMA = "pipeline-quality-report.v1"

__all__ = [
    "analysis_status",
    "evaluate_target",
    "persist_result",
    "required_outputs",
    "required_outputs_v2",
    "write_pipeline_report",
]


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _sanitize_nonfinite(value, *, reasons: list, path_prefix: str = ""):
    """非有限float（NaN/Infinity）をnullへ変換し、理由をreasonsへ記録する。

    `json.dump(..., allow_nan=False)`（`metabolomix.core.atomic_io.atomic_write_json`
    の既定）は非有限floatに出会うと例外で落ちる。呼ぶ前にここで潰す
    （brief拘束「非有限値をnull+理由に変換し、pickleを使わない」）。
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


#: exploratoryでも常に要る基本4種（spec §9.2）。gui_projectはsave_project時のみ、
#: differential系の3種（differential/volcano/tsv）はcomparisonごとに増える。
_BASE_REQUIRED_OUTPUTS = ("preprocess", "pca", "pca_figure", "quality_report")


def required_outputs(request: dict) -> list[str]:
    """spec §9.2の目標別必須出力を、要求から機械的に組み立てる（brief step3）。"""
    names = list(_BASE_REQUIRED_OUTPUTS)
    if request["save_project"]:
        names.append("gui_project")
    if request["effective_target"] == "differential":
        for comparison in request["comparisons"]:
            cid = comparison["comparison_id"]
            names.extend([f"differential:{cid}", f"volcano:{cid}", f"tsv:{cid}"])
    return names


#: v2（LC–MSメタボロミクス）の必須成果物（spec §11）。統計ごとの成果物は
#: 要求の`statistics`から機械的に足す。pathway 15列TSVは**任意**で、対象ゼロの
#: `NO_ANNOTATED_FEATURES`だけでpipelineを未完了にしない。
_V2_REQUIRED_OUTPUTS = (
    "profile",              # 解決済みprofile
    "execution_manifest",   # 実行条件/依存manifest
    "sample_manifest",      # 試料対応表
    "assay_evidence",       # 注入単位証拠表（取得不能なら理由JSON）
    "feature_bindings",     # binding結果（内部標準未使用ならnot_applicableを記録）
    "matrix",               # 前処理matrixと履歴
    "qc_population",        # QC評価集合
    "qc",                   # QC JSON
    "feature_table",        # 全feature定量表
    "quality_report",       # 最終Markdownレポート
)


def required_outputs_v2(request: dict) -> list[str]:
    """v2要求の必須成果物を列挙する（spec §11）。

    v1の`required_outputs`とは別関数にする——v1は「目標(target)」で分岐する
    契約で、v2は「指定した統計の集合」で決まる。1つの関数に両方を詰めると、
    片方の分岐を触ったときにもう片方が黙って変わる。

    統計は`statistic:<statistic_id>`として1件ずつ必須にする。評価不能でも
    「評価不能である」という結果JSONは成果物なので、欠けていれば未完了。
    """
    names = list(_V2_REQUIRED_OUTPUTS)
    for statistic in (request.get("statistics") or []):
        names.append(f"statistic:{statistic['statistic_id']}")
    return names


def analysis_status(results: list[dict]) -> str:
    """指定した統計がどこまで計算できたか（spec §11）。

    execution（回りきったか）とも qc（品質）とも独立の軸。全部計算できれば
    `ready`、一部だけなら `limited`、1つも計算できなければ `not_evaluable`。
    1つでも欠けた状態を `ready` と呼ばない——「解析できた」と読まれる。
    """
    ready = sum(r.get("status") == "ready" for r in results)
    return ("ready" if results and ready == len(results)
            else "limited" if ready else "not_evaluable")


def _load_full_request(record: dict) -> dict:
    """`record["request"]`は要約(revision/request_id/content_hash/saved_path/
    effective_target)しか持たないため、`required_outputs`が要る本体
    (target/comparisons/save_project)は保存済みのpipeline-request.v1から読む
    （`metabolomix.pipeline.store.create_run`が書いた`requests/revision-NNNN.json`）。
    """
    pipeline_root = Path(record["identity"]["pipeline_root"])
    saved_path = record["request"]["saved_path"]
    return json.loads((pipeline_root / saved_path).read_text(encoding="utf-8"))


def _achieved_ref(record: dict, name: str) -> dict | None:
    """`name`に対応する最新のresult refを返す（無ければNone）。

    `results`はappend-only（`metabolomix.pipeline.store`が既存要素の書換えを拒否
    する）ため、resumeで同じoutput_nameが複数revisionにまたがって残ることが
    ある。最後に追記された1件だけを「現在のレポートが指すべき成果物」として
    扱い、古いrevisionの成果物が新しいレポートへ紛れ込まないようにする。
    """
    matches = [r for r in record.get("results", [])
              if isinstance(r, dict) and r.get("output_name") == name]
    return matches[-1] if matches else None


def _upstream_verified(record: dict) -> bool:
    """`record["upstream"]["verification"]`が上流の完了検証済みを示すか判定する。

    Task 3が確立した`completion_status`語彙（"completed"/"partial"/"failed"、
    `metabolomix/console/validation.py`）をそのまま再利用する——新しい語彙は作らない。
    "completed"だけを検証済みとし、未解決（キー自体が無い/None）・"partial"・
    "failed"はいずれも「上流が検証済みでない」として扱う（brief「evaluate_target
    は上流verifiedと全必須outputのhash/IDを検査し」、レビュー指摘3）。
    `_section_execution`と同じ形（dictの"status"、または生の文字列）を読む。
    """
    upstream = record.get("upstream") or {}
    verification = upstream.get("verification")
    status = verification.get("status") if isinstance(verification, dict) else verification
    return status == "completed"


def evaluate_target(record: dict) -> dict:
    """目標別の必須出力が揃っているかを判定する（spec §9.2、brief step3）。

    quality_report作成前の判定にも、report保存後の最終確定にも同じ関数を使う
    ——呼び出し時点の`record`が持つ`results`だけを見るので、report自身が
    まだ`results`に無い時点で呼べば「quality_reportは未達成」と正直に出る。
    これは意図的な仕様（brief「レポートには生成時点のstatusを載せる」）で、
    最終的な`completed`確定は、reportのrefが`results`へ追記された**後**に
    もう一度この関数を呼ぶ側（`report`ステージハンドラ、Task18）の責務。

    差次的目標で比較定義が空の場合は、上流・探索解析の達成状況に関わらず
    `needs_input`/`COMPARISON_REQUIRED`を返す（spec §9.2「比較群を要求しない」
    exploratory目標とは違い、differentialは比較の向きが要る）。

    各出力refのhash/ID検証に加え、`_upstream_verified(record)`
    （`record["upstream"]["verification"]`）も検査する（brief「上流verifiedと
    全必須outputのhash/IDを検査し」）。ref自体が有効でも上流が検証済み
    （"completed"）でなければ「達成」と認めない——`evaluate_target`は単体で
    呼べる関数であり、engineのstage順序が上流成功を保証してくれるとは限らない。
    """
    from metabolomix.pipeline import stage_plan

    request = _load_full_request(record)
    pipeline_root = Path(record["identity"]["pipeline_root"])
    is_v2 = stage_plan.is_v2_request(request)
    target = request.get("effective_target")
    comparisons = request.get("comparisons") or []

    if not is_v2 and target == "differential" and not comparisons:
        return {
            "status": "needs_input",
            "reason_codes": ["COMPARISON_REQUIRED"],
            "required_outputs": required_outputs(request),
            "achieved_outputs": [],
            "missing_outputs": [],
            "output_failures": {},
        }

    # v2は「目標(target)」ではなく「指定した統計の集合」で必須成果物が決まる。
    # v1の分岐はそのまま残す（1つの関数に両方を詰めると、片方を触ったときに
    # もう片方が黙って変わる）。
    required = (required_outputs_v2(request) if is_v2
                else required_outputs(request))
    output_failures = dict(record.get("output_failures") or {})
    upstream_verified = _upstream_verified(record)

    achieved: list[str] = []
    missing: list[str] = []
    reason_codes: list[str] = []

    for name in required:
        ref = _achieved_ref(record, name)
        if ref is not None:
            failures = store.verify_result_refs(pipeline_root, [ref])
            if failures:
                reason_codes.append("RESULT_INTEGRITY_MISMATCH")
                missing.append(name)
                continue
            if not upstream_verified:
                # ref自体は有効でも、上流が検証済みでなければ「達成」と認めない
                # （レビュー指摘3）。hash/IDが揃っていることは、その計算の
                # 前提（Console実行が実際に検証済み完了した）を保証しない。
                reason_codes.append("UPSTREAM_NOT_VERIFIED")
                missing.append(name)
                continue
            achieved.append(name)
            continue
        missing.append(name)
        code = output_failures.get(name)
        reason_codes.append(code if code else "REQUIRED_OUTPUT_MISSING")

    if not missing:
        status = "completed"
    elif achieved:
        status = "partial"
    else:
        status = "failed"

    seen: set[str] = set()
    ordered_reasons = []
    for code in reason_codes:
        if code not in seen:
            seen.add(code)
            ordered_reasons.append(code)

    return {
        "status": status,
        "reason_codes": ordered_reasons,
        "required_outputs": required,
        "achieved_outputs": achieved,
        "missing_outputs": missing,
        "output_failures": {k: v for k, v in output_failures.items() if k in required},
    }


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _atomic_write_text(path: Path, text: str) -> None:
    """同じ親の一時ファイルへ書いてから置換する（`dataset_export.py`と同じ流儀）。

    途中で落ちた出力を「新しい完成品」として残さない。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp",
                                    dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _escape_cell(value) -> str:
    """Markdownテーブルのセル内で改行・パイプ文字が列を壊さないようにする。

    改行は`<br>`（GFMテーブル内で最も広く解釈される改行表現）、パイプは
    `\\|`へ変換する。バックスラッシュ自体も先に変換しないと二重変換になる。
    """
    text = "" if value is None else str(value)
    text = text.replace("\\", "\\\\").replace("|", "\\|")
    text = text.replace("\r\n", "<br>").replace("\n", "<br>").replace("\r", "<br>")
    return text


def _md_table(headers: list[str], rows: list[list], *, empty_note: str = "(no data)") -> str:
    header_line = "| " + " | ".join(headers) + " |"
    sep_line = "| " + " | ".join("---" for _ in headers) + " |"
    if not rows:
        return "\n".join([header_line, sep_line, "| " + " | ".join([empty_note] * len(headers)) + " |"])
    body_lines = ["| " + " | ".join(_escape_cell(cell) for cell in row) + " |" for row in rows]
    return "\n".join([header_line, sep_line, *body_lines])


def _read_result_data(pipeline_root: Path, ref: dict | None):
    """persist_resultが``data``モードで書いたJSONの中身（``data``キーの下）を読む。

    ``path``モード（PNG/TSV等）で登録されたrefや、読めない/形の違うファイルは
    Noneを返す——レポートは「取得できなかった」ことを正直に書くだけで、無い
    情報を埋めない。
    """
    if not ref:
        return None
    try:
        payload = json.loads((Path(pipeline_root) / ref["relative_path"]).read_text(encoding="utf-8"))
    except (OSError, ValueError, KeyError):
        return None
    if isinstance(payload, dict) and "data" in payload:
        return payload["data"]
    return None


_MANIFEST_COLUMNS = ("sample_id", "source_file", "role", "group", "batch",
                    "injection_order", "qc_pool", "include")


def _read_tsv_coverage(pipeline_root: Path, ref: dict | None) -> dict | None:
    """``tsv:<cid>``として登録済みのTSVファイル自身から、InChIKey被覆3件を読む。

    `tsv_summary:<cid>`のような別出力名は`required_outputs`が作らず、
    誰も永続化しない（レビュー指摘1: 死んだコード）。件数は既に
    `export_contract.build_meta`がTSVのメタ行へ書いている
    （`# n_features_total = ...\\tn_with_inchikey = ...\\tn_unannotated = ...`）
    ので、そこを読むだけで足りる——`_read_result_data`は``data``モード
    （JSON）専用でTSVの中身をJSONとして読もうとして失敗するため、別に用意する。
    """
    if not ref:
        return None
    try:
        text = (Path(pipeline_root) / ref["relative_path"]).read_text(encoding="utf-8")
    except (OSError, KeyError):
        return None
    for line in text.splitlines():
        if not line.startswith("# n_features_total"):
            continue
        fields: dict[str, str] = {}
        for part in line[2:].split("\t"):
            key, sep, value = part.partition("=")
            if sep:
                fields[key.strip()] = value.strip()
        try:
            return {
                "n_features_total": int(fields["n_features_total"]),
                "n_with_inchikey": int(fields["n_with_inchikey"]),
                "n_unannotated": int(fields["n_unannotated"]),
            }
        except (KeyError, ValueError):
            return None
    return None


def _section_source(record: dict) -> str:
    identity = record.get("identity") or {}
    lines = ["## Source", "",
             f"- pipeline_id: {identity.get('pipeline_id')}",
             f"- source_root: {identity.get('source_root')}",
             f"- pipeline_root: {identity.get('pipeline_root')}",
             f"- created_at: {identity.get('created_at')}",
             f"- updated_at: {identity.get('updated_at')}"]
    return "\n".join(lines)


def _section_method(record: dict) -> str:
    inputs = record.get("inputs") or {}
    method = inputs.get("method") or {}
    lbm = inputs.get("lbm") or {}
    exe = inputs.get("exe") or {}
    # 研究室 MSP はファイル名だけ（置き場所は載せない。spec 2026-10-07 §3.3）。
    msp = inputs.get("msp") or {}
    lines = [
        "## Method / LBM / Version", "",
        f"- method_source_path: {method.get('source_path') or '(unknown)'}",
        f"- method_sha256: {method.get('sha256') or '(unknown)'}",
        f"- lbm_path: {lbm.get('path') or '(unknown)'}",
        f"- lbm_sha256: {lbm.get('sha256') or '(unknown)'}",
        f"- library_mode: {inputs.get('library_mode') or 'auto'}",
        f"- msp_file: {msp.get('file') or '(not used)'}",
        f"- msp_sha256: {msp.get('sha256') or '(not used)'}",
        f"- msp_source: {msp.get('source') or 'not_used'}",
        f"- exe_version: {exe.get('version') or '(unknown)'}",
    ]
    return "\n".join(lines)


def _section_execution(record: dict) -> str:
    upstream = record.get("upstream") or {}
    verification = upstream.get("verification")
    status = verification.get("status") if isinstance(verification, dict) else verification
    lines = [
        "## Execution Receipt", "",
        f"- console_job_path: {upstream.get('console_job_path') or '(unknown)'}",
        f"- execution_id: {upstream.get('execution_id') or '(unknown)'}",
        f"- verification_status: {status if status is not None else '(unknown)'}",
    ]
    return "\n".join(lines)


def _section_sample_provenance(record: dict) -> str:
    manifest = (record.get("inputs") or {}).get("manifest") or []
    rows = [[row.get(column) for column in _MANIFEST_COLUMNS] for row in manifest]
    table = _md_table(list(_MANIFEST_COLUMNS), rows, empty_note="(no sample manifest recorded)")
    return "\n".join(["## Role / Group / Batch / Order Provenance", "", table])


def _section_applied_skipped(record: dict, pipeline_root: Path, evaluation: dict) -> str:
    """`evaluation["achieved_outputs"]`だけを可否の根拠にする（レビュー指摘2）。

    以前はファイルが読めるかどうかだけを見ており、`evaluate_target`が同じ
    `record`から出すhash整合性検証を無視していた——改変されたファイルでも
    構文的に読めれば「applied_steps」を出してしまい、レポート冒頭の
    `status_at_report_time`（partial/RESULT_INTEGRITY_MISMATCH）と矛盾する。
    """
    lines = ["## Applied / Skipped", ""]
    if "preprocess" not in evaluation["achieved_outputs"]:
        lines.append("(preprocess result not available)")
        return "\n".join(lines)
    data = _read_result_data(pipeline_root, _achieved_ref(record, "preprocess")) or {}
    applied = data.get("applied_steps") or []
    skipped = data.get("skipped_steps") or []
    lines.append(f"- applied_steps: {', '.join(applied) if applied else '(none)'}")
    lines.append(f"- skipped_steps: {', '.join(skipped) if skipped else '(none)'}")
    reasons = data.get("reasons") or {}
    if reasons:
        lines.append("")
        lines.append(_md_table(["step", "reason"], [[k, v] for k, v in sorted(reasons.items())]))
    return "\n".join(lines)


def _section_pca(record: dict, pipeline_root: Path, evaluation: dict) -> str:
    """PCAも同じ理由（レビュー指摘2）で`evaluation["achieved_outputs"]`を見る。

    以前は`ref`の有無だけを見ており、hashが一致しない（改変された）refでも
    result_id・リンク・数値をそのまま出していた。
    """
    lines = ["## PCA", ""]
    if "pca" not in evaluation["achieved_outputs"]:
        lines.append("(PCA result not available)")
        return "\n".join(lines)
    ref = _achieved_ref(record, "pca")
    lines.append(f"- result_id: {ref.get('result_id')}")
    lines.append(f"- source: [{ref['relative_path']}]({ref['relative_path']})")
    data = _read_result_data(pipeline_root, ref)
    if not data:
        return "\n".join(lines)
    evr = data.get("explained_variance_ratio") or []
    lines.append(f"- n_samples: {data.get('n_samples', '(unknown)')}")
    lines.append(f"- n_features: {data.get('n_features', '(unknown)')}")
    lines.append("- explained_variance_ratio: "
                + (", ".join(f"{v:.4f}" for v in evr) if evr else "(unknown)"))
    for caveat in data.get("caveats") or []:
        lines.append(f"- caveat: {caveat}")
    return "\n".join(lines)


def _section_comparisons(record: dict, pipeline_root: Path, request: dict,
                         evaluation: dict) -> str:
    """`status`列は`evaluation["achieved_outputs"]`から取る（レビュー指摘2）。

    以前は`_read_result_data`がJSONとして読めたかどうかだけで
    achieved/missingを決めており、`write_pipeline_report`が同じ`record`から
    既に呼んでいる`evaluate_target`のhash整合性検証と無関係だった。改変されて
    もJSONとして読める限り「achieved」と出てしまい、レポート冒頭の
    `status_at_report_time`（partial/RESULT_INTEGRITY_MISMATCH）と矛盾する
    ——一つの`record`から一つの整合性判定だけを使う。
    """
    comparisons = request.get("comparisons") or []
    lines = ["## Comparisons", ""]
    if not comparisons:
        lines.append("(no comparisons requested)")
        return "\n".join(lines)
    headers = ["comparison_id", "reference_group", "test_group", "status", "unadjusted_confounded"]
    rows = []
    for comparison in comparisons:
        cid = comparison["comparison_id"]
        name = f"differential:{cid}"
        achieved = name in evaluation["achieved_outputs"]
        data = _read_result_data(pipeline_root, _achieved_ref(record, name)) if achieved else None
        prov_comparison = ((data or {}).get("provenance") or {}).get("comparison") or {}
        rows.append([
            cid, comparison.get("reference_group"), comparison.get("test_group"),
            "achieved" if achieved else "missing",
            "yes" if prov_comparison.get("unadjusted_confounded") else "no",
        ])
    lines.append(_md_table(headers, rows))
    return "\n".join(lines)


def _section_inchikey_coverage(record: dict, pipeline_root: Path, request: dict) -> str:
    comparisons = request.get("comparisons") or []
    lines = ["## InChIKey Coverage", ""]
    if not comparisons:
        lines.append("(no comparisons requested)")
        return "\n".join(lines)
    headers = ["comparison_id", "n_features_total", "n_with_inchikey", "n_unannotated"]
    rows = []
    for comparison in comparisons:
        cid = comparison["comparison_id"]
        coverage = _read_tsv_coverage(pipeline_root, _achieved_ref(record, f"tsv:{cid}"))
        if coverage:
            rows.append([cid, coverage["n_features_total"], coverage["n_with_inchikey"],
                        coverage["n_unannotated"]])
        else:
            rows.append([cid, "(n/a)", "(n/a)", "(n/a)"])
    lines.append(_md_table(headers, rows))
    return "\n".join(lines)


def _section_unverified(record: dict) -> str:
    lines = ["## Unverified Conditions", ""]
    warnings = record.get("warnings") or []
    if warnings:
        rows = [[w.get("code"), w.get("stage_id"), w.get("message")] for w in warnings]
        lines.append(_md_table(["code", "stage_id", "message"], rows))
    else:
        lines.append("(no recorded warnings)")
    needs_input = record.get("needs_input")
    if needs_input:
        lines.append("")
        lines.append(f"- needs_input: {needs_input.get('code')} ({needs_input.get('message')})")
    return "\n".join(lines)


def write_pipeline_report(record: dict, path) -> dict:
    """決定的なテンプレートで品質レポートを書く（spec §11、brief step3）。

    セクション順は固定（source→method/LBM/版→実行証跡→role/group/batch/order
    出所→applied/skip→PCA→比較→InChIKey被覆→未検証条件）。LLM生成文には一切
    依存せず、`record`が持つ機械可読な値をそのまま並べるだけ——生物学的解釈や
    MS/MS確認の自動判定はしない。

    `evaluate_target(record)`をこの関数の内側で**一度だけ**呼ぶ。このレポート
    自身のrefはまだ`record["results"]`に無い時点の呼び出しなので、返る
    `status`は「report作成前の時点の状態」であり、`quality_report`自体は
    まだ未達成として出る——これは意図的（brief「レポートには生成時点の
    statusを載せる」「最終completedはreport保存後のstate更新で確定する」）。
    最終的な`completed`確定は、このrefが`results`へ追記された後にもう一度
    `evaluate_target`を呼ぶ側の責務（呼び出し元がpersist_resultで登録し、
    改めて判定する）。
    """
    path = Path(path)
    pipeline_root = Path(record["identity"]["pipeline_root"])
    request = _load_full_request(record)
    evaluation = evaluate_target(record)

    header = [
        "# Pipeline Quality Report", "",
        f"- schema: {REPORT_SCHEMA}",
        f"- generated_at: {_now_iso()}",
        f"- status_at_report_time: {evaluation['status']}",
        f"- reason_codes: {', '.join(evaluation['reason_codes']) or '(none)'}",
        "",
        "This report is generated from a fixed deterministic template. It does not "
        "contain automated biological interpretation or MS/MS identification "
        "confirmation.",
    ]

    sections = [
        _section_source(record),
        _section_method(record),
        _section_execution(record),
        _section_sample_provenance(record),
        _section_applied_skipped(record, pipeline_root, evaluation),
        _section_pca(record, pipeline_root, evaluation),
        _section_comparisons(record, pipeline_root, request, evaluation),
        _section_inchikey_coverage(record, pipeline_root, request),
        _section_unverified(record),
    ]

    text = "\n\n".join(["\n".join(header), *sections]) + "\n"
    _atomic_write_text(path, text)

    return {
        "path": str(path),
        "status": evaluation["status"],
        "reason_codes": evaluation["reason_codes"],
        "required_outputs": evaluation["required_outputs"],
        "achieved_outputs": evaluation["achieved_outputs"],
        "missing_outputs": evaluation["missing_outputs"],
    }


def _result_path_for(pipeline_root: Path, relative: str, payload: dict) -> Path:
    """`relative` へ書く。既存と内容が違うなら内容 hash を足した別名を返す。

    既存が無い、または既存の中身が完全に同じなら `relative` をそのまま使う。
    別名は `<stem>__<先頭8桁>.json` で、同じ内容なら同じ名前に落ちる
    （決定論的——再実行のたびに名前が増えることはない）。
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
    """成果物への参照(ref)を確定する（Task 8が残した「manifestへの登録」の穴）。

    `result` は次のいずれかの形:

    - 既に書き終えたファイル（PNG/TSV等）を登録するだけ:
      ``{"output_name", "kind", "result_id", "path", parent_ids?, request_revision?}``。
      ここではhashを計算し、pipeline_root相対へ書き換えるだけで、ファイルには
      一切書き込まない。
    - この場でJSONとして書く数値結果（PCA/差次的解析等）:
      ``{"output_name", "kind", "result_id", "data", relative_path?,
      parent_ids?, request_revision?}``。非有限値はnull+理由へ変換してから
      `atomic_write_json`（pickle不使用）で書く。

    戻り値は `result_id/kind/relative_path/hash/parent_ids/request_revision` を
    持つref（`metabolomix.pipeline.store.verify_result_refs` がそのまま検証できる
    形——store側は`relative_path`/`hash`というキー名で読む）。加えて
    `output_name`（このrefがどの必須出力に対応するかの照合キー、
    `evaluate_target` が使う）を持つ。
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
        # 既に記録された ref の (relative_path, hash) は永久に検証可能でなければ
        # ならない。`rerun_upstream=true` の再実行は同じ result_id を使い
        # （pipeline_id 由来で request_revision も上がらない）、素朴に同じパスへ
        # 書くと先の ref が恒久的に hash 不一致になり `verify_result_refs` が
        # 必ず落ちる。内容が変わったときだけ別名へ逃がす——内容が同じ再実行では
        # ファイルを増やさない（冪等）。
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
