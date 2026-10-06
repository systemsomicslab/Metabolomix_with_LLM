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

from metabolomix.analysis import feature_bindings, matrix_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading as mztab_loading
from metabolomix.pipeline import recovery, store

__all__ = ["load_run_dataset"]

#: `metabolomics_handlers._MATRIX_SUBDIR` と同じ置き場（保存側と読み側の対）。
_MATRIX_SUBDIR = ("results", "matrices")


def load_run_dataset(pipeline_path, *, allow_incomplete: bool = False):
    """run の成果物から DatasetState を組み立てて返す（保存も session 更新もしない）。

    `pipeline_path` は `pipeline-run.json` でもその親ディレクトリでもよい
    （`pipeline_status` と同じ規約）。行列は `save_matrix` が残した参照から
    読むので、meta・配列のどちらが書き換わっても `load_matrix` が拒否する。
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
    return ds
