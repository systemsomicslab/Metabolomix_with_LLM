"""v1 の前処理（`ds.pp_matrix` を書く経路）をテストから呼ぶための足場。

`dataset_preprocess` ツールは統計エンジンの前処理（`ds.stat_matrices`）になった。まだ v1 の
`pp_matrix` を読む利用者（旧 2 群比較 / `dataset_pca`・`pathway_project`・書き出し・
根拠バンドル・特徴検索）のテストは、ここから `dataset_service.preprocess_dataset`（旧ツールが呼んでいた
v1 のサービス）を同じ引数で呼んで、その状態を作る。D 段で v1 を消すときにこのファイルごと消える。
"""
from __future__ import annotations

from metabolomix.analysis.dataset_analysis import PreconditionError
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.tools._holder import resolve_holder


def _count_roles(roles: dict, sample_names: list[str]) -> dict:
    counts: dict[str, int] = {}
    for name in sample_names:
        role = roles.get(name, "sample")
        counts[role] = counts.get(role, 0) + 1
    return counts


def v1_preprocess(normalize: str = "none", blank_min_fold: float | None = None,
                  drift_correct: bool = False, max_qc_rsd: float | None = None,
                  impute: str = "half_min", min_detection_rate: float = 0.0,
                  dataset_handle: str | None = None) -> str:
    """旧 `dataset_preprocess` ツールと同じ引数・同じ戻り値（JSON 文字列）で v1 の前処理を回す。"""
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    from metabolomix.analysis.dataset_service import preprocess_dataset

    recipe = {"normalize": normalize, "blank_min_fold": blank_min_fold,
              "drift_correct": drift_correct, "max_qc_rsd": max_qc_rsd, "impute": impute,
              "min_detection_rate": min_detection_rate}
    try:
        report = preprocess_dataset(ds, recipe)
    except PreconditionError as exc:
        if exc.kind == "missing_state":
            return missing_state(exc.state, ["dataset_preprocess"], exc.message)
        return mztab_error("DATASET_BAD_REQUEST", exc.message, exc.details or None)
    return json_payload({
        "status": "success",
        "result_id": report["provenance"]["result_id"],
        "n_samples": len(ds.pp_sample_names),
        "n_features": len(ds.pp_feature_names),
        "features_before": report.get("features_before"),
        "features_removed_total": report.get("features_removed_total"),
        "recipe_applied": report.get("recipe_applied", []),
        "excluded_from_matrix": report.get("excluded_from_matrix", {}),
        "role_counts": _count_roles(ds.roles, ds.pp_sample_names),
        "steps": report.get("steps", {}),
        "caveats": report.get("caveats", []),
        "next": "dataset_pca または dataset_compare を実行してください",
        **({"detection": report["detection"]} if "detection" in report else {}),
        **({"detection_filter": report["detection_filter"]}
           if "detection_filter" in report else {}),
    })


def v1_pca(n_components: int = 5, dataset_handle: str | None = None) -> dict:
    """v1 の PCA（`ds.last_pca` に scores / loadings / singular_values を置く経路）を回す。

    pipeline の引き渡しが今も作る形。`v1_preprocess` の後に呼ぶ。D 段で v1 を消すときにこのファイルごと消える。
    """
    holder = resolve_holder(dataset_handle)
    ds = holder.dataset
    from metabolomix.analysis.dataset_service import pca_dataset

    return pca_dataset(ds, n_components=n_components)


def v1_differential(group_a: list[str], group_b: list[str], q_threshold: float = 0.05,
                    log2fc_threshold: float = 1.0, log_transform: bool = True,
                    group_a_label: str = "group_a", group_b_label: str = "group_b",
                    dataset_handle: str | None = None, handoff: bool = True) -> str:
    """pipeline の引き渡し（`session_handoff.restore_differentials`）が session へ戻す v1 の 2 群比較を再現する。

    v1 の `two_group` 結果が MCP の呼び出しだけで session に載る道は、この削除のあと無い（対話のツールは
    統計エンジンの `dataset_compare`）。載るのは pipeline run の引き渡しだけで、その結果は
    `provenance.handoff`（run の前処理レシピ・比較 ID）と引き継いだ時点の `metadata_revision` を持つ。
    ここでは `dataset_service.compare_dataset` で結果を作り、同じ印を付けて旧ツールと同じ引数・同じ戻り値
    （要約の JSON 文字列）で返す。結果の読まれ方（書き出し・根拠バンドル・特徴検索・volcano）を確かめる
    テスト専用で、MCP から作れる状態を確かめるテストは `dataset_compare` を使うこと。
    D 段で v1 を消すときにこのファイルごと消える。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    from metabolomix.analysis.dataset_service import compare_dataset

    try:
        result = compare_dataset(
            ds, group_a, group_b, q_threshold=q_threshold,
            log2fc_threshold=log2fc_threshold, log_transform=log_transform,
            group_a_label=group_a_label, group_b_label=group_b_label)
    except PreconditionError as exc:
        if exc.kind == "missing_state":
            return missing_state(exc.state, ["dataset_preprocess"], exc.message)
        return mztab_error("DATASET_BAD_REQUEST", exc.message, exc.details or None)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    provenance = result["provenance"]
    if not handoff:
        # 同じ session で計算した v1 の結果（親の前処理が変われば古くなる）。session で前処理を
        # やり直すと古くなる規則そのものを見るテストだけが使う。
        return _summary(result)
    parent = (ds.results or {}).get((provenance.get("parent_ids") or [None])[0]) or {}
    provenance["handoff"] = {
        "pipeline_root": "(test)", "comparison_id": "c1",
        "preprocess": dict((parent.get("provenance") or {}).get("effective_parameters") or {}) or None}
    provenance["metadata_revision"] = getattr(ds, "metadata_revision", 0)
    return _summary(result)


def _summary(result: dict) -> str:
    payload = {k: v for k, v in result.items()
               if k not in ("results", "volcano", "provenance")}
    payload["status"] = "success"
    payload["result_id"] = result["provenance"]["result_id"]
    payload["differential_contract_version"] = result["contract_version"]
    payload["results_note"] = (
        "全特徴の結果と volcano 点列は本要約に非同梱（セッション保持）。"
        "InChIKey 付きの全行が必要なら dataset_export_differential を実行してください。")
    return json_payload(payload)
