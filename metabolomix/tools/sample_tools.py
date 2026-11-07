"""試料の MCP ツール: `dataset_exclude`（試料・特徴の除外）。

試料の確認は `dataset_status`（`view="samples"`）、実験情報シートによる訂正は
`dataset_set_sample_metadata`。ここは**解析から外す**ことだけを持つ。数値・状態の
規則は `analysis/dataset_exclusions.py`。
"""
from __future__ import annotations

from mcp.types import ToolAnnotations

from metabolomix.analysis import dataset_exclusions
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_exclude"]

#: 返す除外の一覧の上限（特徴は数が多くなりうる）。
_LIST_CAP = 50


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def dataset_exclude(
    samples: list[str] | None = None,
    features: list[str] | None = None,
    mode: str = "add",
    reason: str | None = None,
    dataset_handle: str | None = None,
) -> str:
    """試料・特徴を解析から外す（可逆）。`dataset_preprocess` が行列から外して適用する。

    samples: 試料名（`dataset_status(view="samples")` の name）または実験情報シートの sample_id。
    features: 特徴 ID（`dataset_feature_search` の feature_id）。
    mode: `"add"`（既定。除外に足す。**`reason` が要る**。除外済みは理由を置き換える）/
        `"remove"`（除外から戻す）/ `"clear"`（`samples` / `features` のうち渡した種類を空にする。
        どちらも渡さなければ両方）。
    reason: 除外の理由（`add` で必須。後で戻す・説明するための記録で、`dataset_status` の
        `samples` の `excluded` 列に出る）。

    **存在しない試料名・特徴 ID が 1 つでもあれば全体を拒み、何も変更しない**。全試料の除外も拒む。
    入力（定量行列・特徴 ID・試料名）は変えず、`dataset_fingerprint` も動かさない。

    実際に除外が変わったときだけ、前処理済み行列から先の結果（PCA・差次的解析）を古いとして
    捨てる（`dataset_preprocess` からやり直す。古い結果は `ds.results` に残るが現在の結果として
    は通らない）。同じ引数の再実行は何も動かさない（冪等）。結果の件数を返す。

    dataset_handle: `load_dataset` が返したハンドル（再解析の並行のエージェント向け）。
        省略時は対話のセッションのデータセット。ハンドルの状態だけを読み書きする。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    try:
        outcome = dataset_exclusions.apply_exclusions(
            ds, samples=samples, features=features, mode=mode, reason=reason)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    payload = {"status": "success", "mode": mode, **outcome}
    for kind, excluded in (("samples", ds.excluded_samples), ("features", ds.excluded_features)):
        if excluded:
            shown = dict(list(excluded.items())[:_LIST_CAP])
            payload[f"excluded_{kind}"] = shown
            if len(excluded) > len(shown):
                payload[f"excluded_{kind}_truncated"] = True
    if outcome["changed"]:
        payload["next"] = ("前処理済み行列から先の結果を古いとして捨てました。"
                           "dataset_preprocess からやり直してください。")
    return json_payload(payload)
