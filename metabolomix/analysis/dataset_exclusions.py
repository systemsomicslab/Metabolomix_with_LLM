"""データセットの試料・特徴の除外（`dataset_exclude` の数値・状態層。MCP 非依存）。

除外は `DatasetState.excluded_samples` / `excluded_features`（名前または特徴 ID -> 理由）に
持つ。**可逆**で、入力の指紋（`result_state.dataset_fingerprint`）には入れない
（`sample_meta` の編集と同じ解析の状態）。前処理（`dataset_analysis`）が行列から外す。

約束は 2 つ:

1. **全体を検証してから反映する。** 存在しない試料名・特徴 ID が 1 つでもあれば全体を拒み、
   一部だけ当てない。
2. **実際に変わったときだけ**、`metadata_revision` を進めて前処理済み行列から先の結果を古くする
   （`result_state.invalidate_results`。`sample_manifest.apply_metadata` と同じ経路）。同じ引数の
   再実行や、何も外れない `remove` は状態を動かさない（冪等）。
"""
from __future__ import annotations

from metabolomix.analysis.result_state import invalidate_results
from metabolomix.core.atomic_io import DomainError

__all__ = ["MODES", "apply_exclusions"]

MODES = ("add", "remove", "clear")


def _bad(message: str, details: dict | None = None) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details or {})


def _resolve_samples(ds, requested: list[str]) -> tuple[list[str], list[str]]:
    """試料名または sample_id -> 試料名。(解決した名前, 見つからなかったもの)。

    名前が優先で、無ければ `ds.sample_ids`（実験情報シートの sample_id）で引く。保存は
    いつも試料名（前処理が名前で外す）。
    """
    names = list(getattr(ds, "sample_names", []) or [])
    known = set(names)
    by_id = {sid: name for sid, name in zip(getattr(ds, "sample_ids", []) or [], names)}
    resolved: list[str] = []
    unknown: list[str] = []
    for item in requested:
        name = item if item in known else by_id.get(item)
        if name is None:
            unknown.append(item)
        elif name not in resolved:
            resolved.append(name)
    return resolved, unknown


def _resolve_features(ds, requested: list[str]) -> tuple[list[str], list[str]]:
    known = set(getattr(ds, "feature_ids", []) or [])
    resolved: list[str] = []
    unknown: list[str] = []
    for item in requested:
        if item not in known:
            unknown.append(item)
        elif item not in resolved:
            resolved.append(item)
    return resolved, unknown


def apply_exclusions(ds, *, samples: list[str] | None, features: list[str] | None,
                     mode: str, reason: str | None) -> dict:
    """除外を足す・戻す・空にする。結果の件数と、変わったかを返す。

    - `add`: `reason` が要る（空は拒む。記録のない除外は、後で戻す・説明する材料が無い）。
      既に除外済みの名前は理由を新しいものに置き換える。
    - `remove`: 指定した名前・ID を除外から戻す（除外されていないものは何もしない）。
    - `clear`: 指定した種類（`samples` / `features` を渡したほう）を空にする。どちらも
      渡さなければ両方。名前は検証しない。

    全部の試料を除外する指定は拒む（解析に残る試料が無い）。
    """
    if mode not in MODES:
        raise _bad(f"mode は {list(MODES)} のどれかです: {mode!r}", {"mode": mode})

    new_samples = dict(getattr(ds, "excluded_samples", None) or {})
    new_features = dict(getattr(ds, "excluded_features", None) or {})

    if mode == "clear":
        clear_samples = samples is not None or features is None
        clear_features = features is not None or samples is None
        if clear_samples:
            new_samples = {}
        if clear_features:
            new_features = {}
    else:
        if not samples and not features:
            raise _bad("samples か features に、対象の試料名（または sample_id）・特徴 ID を"
                       "1 つ以上渡してください。")
        text = (reason or "").strip()
        if mode == "add" and not text:
            raise _bad("mode=\"add\" には reason（除外の理由）が要ります。"
                       "後で戻す・説明するための記録です。")
        resolved_samples, unknown_samples = _resolve_samples(ds, list(samples or []))
        resolved_features, unknown_features = _resolve_features(ds, list(features or []))
        if unknown_samples or unknown_features:
            raise _bad(
                "存在しない試料名・特徴 ID があるため、何も変更していません。",
                {**({"unknown_samples": unknown_samples} if unknown_samples else {}),
                 **({"unknown_features": unknown_features} if unknown_features else {}),
                 "hint": "試料は dataset_status(view=\"samples\")、特徴は dataset_feature_search で確かめる。"})
        for target, items in ((new_samples, resolved_samples), (new_features, resolved_features)):
            for item in items:
                if mode == "add":
                    target[item] = text
                else:
                    target.pop(item, None)

    n_samples = len(getattr(ds, "sample_names", []) or [])
    if n_samples and len(new_samples) >= n_samples:
        raise _bad("全ての試料を除外することはできません（解析に残る試料が無くなります）。"
                   "何も変更していません。", {"n_samples": n_samples})

    changed = (new_samples != (getattr(ds, "excluded_samples", None) or {})
               or new_features != (getattr(ds, "excluded_features", None) or {}))
    if changed:
        # 前処理済み行列から先の結果を古くする（`apply_metadata` の無効化と同じ経路）。
        invalidate_results(ds, {"exclusion"})
        ds.excluded_samples = new_samples
        ds.excluded_features = new_features
        ds.metadata_revision = getattr(ds, "metadata_revision", 0) + 1
    return {
        "changed": changed,
        "n_excluded_samples": len(new_samples),
        "n_excluded_features": len(new_features),
        "n_samples_remaining": n_samples - len(new_samples),
        "n_features_remaining": len(getattr(ds, "feature_ids", []) or []) - len(new_features),
    }
