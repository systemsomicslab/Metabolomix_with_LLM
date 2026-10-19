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

#: `service._handle_preprocess` が前処理の成果物に付ける output_name。
_PREPROCESS_OUTPUT = "preprocess"


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

    # 群・batch・注入順が無ければ、載せた行列で統計を組めない。
    rows = _run_manifest_rows(pipeline_root, record)
    if rows is not None:
        ds.sample_metadata_rows = rows

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


def _run_manifest_rows(pipeline_root: Path, record: dict) -> list[dict] | None:
    """run が解決した試料対応表の行を返す。どこにも無ければ None。

    v2 は `sample_manifest` を hash 付きの result ref として残すので、そちらを先に読む。
    v1 はそれを作らず、`resolve_metadata` が `record["inputs"]["manifest"]` にだけ
    書く（`metabolomics_handlers.sample_manifest_outcome`）——どちらも
    `resolve_sample_metadata` の同じ行なので、無ければそこへフォールバックする。
    フォールバックしないと群の指紋が「群なし」になり、v1 の差次的結果が全部
    `group_mismatch` で戻らない。`source_file` の絶対/相対の違いは群の指紋
    （`group` と `biological_sample_id` だけを見る）に効かない。
    """
    manifests = store.read_result_data(pipeline_root, record.get("results") or [],
                                       "sample_manifest")
    if manifests:
        return list((manifests[-1] or {}).get("rows") or [])
    rows = (record.get("inputs") or {}).get("manifest")
    if isinstance(rows, list) and rows:
        return [dict(row) for row in rows if isinstance(row, dict)]
    return None


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


def _run_preprocess_parameters(pipeline_root: Path, record: dict,
                               parent_ids: list) -> dict | None:
    """差次的結果の親（run の前処理）が実際に使ったレシピを返す。読めなければ None。

    書き出しのメタ行 `# preprocess = ...` は本来 `ds.results` の親結果から引くが、
    run の前処理結果は session へ戻さない（`is_current` の handoff 分岐は kind を
    見ないので、戻すと「今の前処理」に見えてしまう）。そこでレシピだけを
    引き継ぎ結果の来歴へ写す。写さないと、書き出しは session の
    `preprocessing_recipe`（空、または session で後からやり直した別のレシピ）を
    書き、正しい数字に違う条件のラベルが付く。

    ここは `current_result_refs` ではなく**履歴を result_id で引く**。前処理は
    resume のたびに作り直されて新しい result_id で永続化される一方、成果物で
    検証が通った差次的工程は skip されるので、差次的結果の親は最後の前処理とは
    限らない。前処理の成果物は result_id ごとに別ファイル（`results/<id>.json`）
    なので、古い ref も hash 照合できる。照合できなければ None（不明）を返し、
    読み込み自体は止めない——数字は差次的結果側で検証済みで、欠けるのはラベルだけ。
    """
    wanted = set(parent_ids or [])
    for ref in reversed(record.get("results") or []):
        if not isinstance(ref, dict) or ref.get("output_name") != _PREPROCESS_OUTPUT:
            continue
        if ref.get("result_id") not in wanted:
            continue
        try:
            data = store.read_result_data(pipeline_root, [ref], _PREPROCESS_OUTPUT)[0]
        except (DomainError, OSError, ValueError):
            return None
        # `preprocess_auto` が `_commit_preprocess` に渡したレシピ＝来歴の
        # `effective_parameters` と同じもの（`_handle_preprocess` は provenance を
        # 落として plan を保存するので、ここから読む）。
        recipe = (data or {}).get("resolved_recipe")
        return dict(recipe) if isinstance(recipe, dict) else None
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
        provenance["handoff"] = {
            "pipeline_root": str(pipeline_root),
            "comparison_id": comparison_id,
            # None は「run の前処理結果を読めなかった」。書き出しは不明と書く。
            "preprocess": _run_preprocess_parameters(
                pipeline_root, record, provenance.get("parent_ids")),
        }
        result["provenance"] = provenance
        ds.results[result_id] = result
        restored.append({"comparison_id": comparison_id, "result_id": result_id,
                         "group_a": result.get("a"), "group_b": result.get("b")})
    return {"restored": restored, "rejected": rejected}
