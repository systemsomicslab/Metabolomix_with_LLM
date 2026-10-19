"""完了した run の成果物を、対話セッションが使える DatasetState へ戻す層。

`dataset_load(pipeline_path=...)` だけがここを呼ぶ。行列だけを別ツールで
載せられる形にしないのは、run A の dataset に run B の行列を混ぜられると
群の対応が黙ってずれるため——dataset・解決済み試料対応表・binding・解析行列は
同じ run から**まとめて**載せる。

グローバル session は import しない（この層は DatasetState を作って返すだけで、
どこへ置くかは呼び出し側＝MCP 公開層が決める）。
"""
from __future__ import annotations

from pathlib import Path

from metabolomix.analysis import feature_bindings, matrix_state, result_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import recovery, store

__all__ = ["load_run_dataset"]

#: `metabolomics_handlers._MATRIX_SUBDIR` と同じ置き場（保存側と読み側の対）。
_MATRIX_SUBDIR = ("results", "matrices")

#: `service._persist_differential` が差次的結果に付ける output_name の接頭辞。
_DIFFERENTIAL_PREFIX = "differential:"


def load_run_dataset(pipeline_path, *, allow_incomplete: bool = False):
    """run の成果物から DatasetState を組み立てて返す（保存も session 更新もしない）。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい
    （`pipeline_status` と同じ規約）。行列は `save_matrix` が残した参照から
    読むので、meta・配列のどちらが書き換わっても `load_matrix` が拒否する。
    差次的結果は指紋（データ・群）が一致するものだけを `ds.results` へ戻し、
    戻したもの・戻さなかったものを `ds.handoff_report` に残す
    （`restore_differentials`）。
    """
    pipeline_root = recovery.normalize_pipeline_root(Path(pipeline_path).expanduser())
    record = store.load_run(pipeline_root)

    job_path = (record.get("upstream") or {}).get("console_job_path")
    if not job_path:
        raise DomainError(
            "PIPELINE_UPSTREAM_INCOMPLETE",
            "この run はまだ上流（Console）の成果物を持っていません。"
            "pipeline_status で状態を確認してください。",
            {"pipeline_root": str(pipeline_root)})

    ds = mztab_loading.load_dataset_state(job_path=Path(job_path),
                                          allow_incomplete=allow_incomplete)
    results = record.get("results") or []

    manifests = store.read_result_data(pipeline_root, results, "sample_manifest")
    if manifests:
        # 群・batch・注入順が無ければ、載せた行列で統計を組めない。
        ds.sample_metadata_rows = list((manifests[-1] or {}).get("rows") or [])

    bindings = store.read_result_data(pipeline_root, results, "feature_bindings")
    if bindings:
        ds.feature_binding_targets = feature_bindings.resolved_targets(bindings[-1])

    directory = pipeline_root.joinpath(*_MATRIX_SUBDIR)
    for data in store.read_result_data(pipeline_root, results, "matrix"):
        storage = (data or {}).get("storage")
        if not storage:
            # 値を保存していない古い run。要約しか無いので復元しない
            # （空の行列を置くと、統計が「行列はある」と誤読する）。
            continue
        matrix = matrix_state.load_matrix(storage, directory)
        ds.analysis_matrices[matrix["matrix_id"]] = matrix

    # 群の指紋は試料対応表を戻した後でなければ run と同じ規則で計算できない。
    ds.handoff_report = restore_differentials(ds, pipeline_root, record)
    return ds


def _handoff_rejection(ds, provenance: dict) -> str | None:
    """引き継げない理由を返す（引き継げるなら None）。

    判定の規則そのものは `result_state.is_current` の handoff 分岐と同じ。
    ここでは「なぜ戻さなかったか」を利用者に言えるよう、理由を分けて返す。
    """
    recorded_data = provenance.get("dataset_hash")
    if recorded_data is None:
        return "legacy_result_without_dataset_hash"
    if recorded_data != result_state.dataset_fingerprint(ds):
        return "dataset_mismatch"
    recorded_groups = provenance.get("group_fingerprint")
    if recorded_groups is not None and recorded_groups != result_state.group_fingerprint(ds):
        return "group_mismatch"
    return None


def restore_differentials(ds, pipeline_root: Path, record: dict) -> dict:
    """run の差次的結果のうち、指紋が今の ds と一致するものを `ds.results` へ戻す。

    読むのは `store.current_result_refs` が選ぶ**今の ref だけ**——同じ比較が
    export で再計算されると `differential:<cid>` が履歴に何本も並び、上書きされた
    古い ref の hash はもうファイルと一致しない（それを「壊れた」と読まない）。
    今の ref の成果物が改変されていれば、行列と同じく `RESULT_INTEGRITY_MISMATCH`
    で読み込みごと止める。

    result_id は付け替えない（pipeline の TSV のメタ行 `# result_id` と同じ ID で
    引けるように）。`ds.last_differential` は触らない——「直近の結果」は session で
    計算したものに限り、戻した結果は result_id で名指しして使う。
    """
    restored: list[dict] = []
    rejected: list[dict] = []
    for ref in store.current_result_refs(record):
        output_name = ref.get("output_name") or ""
        if not output_name.startswith(_DIFFERENTIAL_PREFIX):
            continue
        comparison_id = output_name.split(":", 1)[1]
        result = store.read_result_data(pipeline_root, [ref], output_name)[0] or {}
        provenance = dict(result.get("provenance") or {})
        result_id = provenance.get("result_id") or ref.get("result_id")
        reason = _handoff_rejection(ds, provenance)
        if reason is not None:
            rejected.append({"comparison_id": comparison_id, "result_id": result_id,
                             "reason": reason})
            continue
        provenance["handoff"] = {"pipeline_root": str(pipeline_root),
                                 "comparison_id": comparison_id}
        result["provenance"] = provenance
        ds.results[result_id] = result
        restored.append({"comparison_id": comparison_id, "result_id": result_id,
                         "group_a": result.get("a"), "group_b": result.get("b")})
    return {"restored": restored, "rejected": rejected}
