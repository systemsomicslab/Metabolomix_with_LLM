"""統計エンジンの結果の登録と、2 群の結果への読み替え（構成見直し C の配線の共通部品）。

対話のツール（前処理・比較・PCA）と、差次的結果を読む既存の利用者（書き出し・根拠バンドル・
特徴検索・volcano）の間に立つ。MCP・セッションには依存しない（`ds` を受け取って書くだけ）。

- 行列: `register_matrix` が前処理済みの行列（`PreparedMatrix`）を設定・群の上書き・除外・指紋と一緒に
  `ds.stat_matrices[matrix_id]` へ置き、`get_matrix` が古ければ None を返す。
- 結果: `register_compare` / `register_pca` が `ds.results[result_id]` と `ds.last_differential` /
  `ds.last_pca` に置く。`result_id` は決定的（同じ計算は同じ鍵を上書きする）。
  来歴は `{result_id, kind: "compare"|"stat_pca", dataset_hash, metadata_revision, parent_ids: [matrix_id]}`
  で、古さは `is_stats_current`（`result_state.is_current` もここへ回す）。
- 読み替え: `as_two_group` が `compare` の記録の 1 つの対比を、v1 の 2 群の結果と同じ欄に写す。

試料は**表示名**（`ds.sample_names`）で持つ。エンジンの試料 ID（試料の表があれば `ds.sample_ids`）は
位置で表示名へ結ぶ（既存の利用者はすべて表示名で試料を引く）。
"""
from __future__ import annotations

import copy
import math

from metabolomix.analysis import differential, export_contract
from metabolomix.analysis.result_state import dataset_fingerprint, is_stats_current
from metabolomix.analysis.stats.compare import _model_value
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "as_two_group",
    "contrast_ids",
    "get_matrix",
    "is_stats_current",
    "register_compare",
    "register_matrix",
    "register_pca",
    "study_inputs",
    "unit_group_labels",
]


# --- 入力（除外と群の上書き） ------------------------------------------------------

def _facets(ds) -> dict:
    """試料の表示名 → 因子トークンの束。役割と群は試料の対応表の行から、無ければ読み込みの既定の行から
    取る（`plots/dataset_rows._facets` と同じ作り）。"""
    from metabolomix.analysis import sample_manifest
    from metabolomix.msdial import sample_factors

    names = list(ds.sample_names)
    table = getattr(ds, "sample_metadata_rows", None) or sample_manifest.resolve_metadata(ds, None)
    meta = {name: {"role": (row or {}).get("role"), "group": (row or {}).get("group")}
            for name, row in zip(names, table)}
    return sample_factors.build_sample_facets(names, None, meta)


def study_inputs(ds, *, group_levels=None, group_factors=None) -> dict:
    """エンジンの入力に写す 3 つ: `{"exclude_samples", "exclude_features", "group_overrides"}`。

    除外は `ds.excluded_samples` / `ds.excluded_features`（`dataset_exclude` の状態）を名前の昇順で。
    群の上書きは `group_levels`（1 因子の値）か `group_factors`（因子軸のリスト。優先）を
    `msdial.sample_factors.assign_factor_groups` で試料の表示名に当てた結果（試料名 → 群）。どちらも
    無ければ `{}`（試料の表の群をそのまま使う）。役割が sample でない試料（QC・ブランク）は上書きしない。
    値が相互排他でない（1 つの試料が 1 軸の複数の値に当たる）指定は `DATASET_BAD_REQUEST`。
    これで作った dict は spec v2 の `input.exclude_samples` / `exclude_features` / `group_overrides`
    にそのまま入る（spec だけで同じ数字を出し直せる）。
    """
    from metabolomix.msdial import sample_factors

    overrides: dict[str, str] = {}
    if group_levels or group_factors:
        facets = _facets(ds)
        try:
            groups = sample_factors.assign_factor_groups(
                facets, group_factors=group_factors, group_levels=group_levels)
        except ValueError as exc:
            raise DomainError("DATASET_BAD_REQUEST", str(exc),
                              {"group_levels": group_levels, "group_factors": group_factors}) from exc
        overrides = {name: label for name, label in groups.items()
                     if label and facets[name].role == "sample"}
    return {
        "exclude_samples": sorted(getattr(ds, "excluded_samples", None) or {}),
        "exclude_features": sorted(getattr(ds, "excluded_features", None) or {}),
        "group_overrides": dict(sorted(overrides.items())),
    }


# --- 行列 ---------------------------------------------------------------------

def _entry_current(entry: dict, dataset_hash: str, revision) -> bool:
    return (entry.get("metadata_revision") == revision
            and entry.get("dataset_hash") == dataset_hash)


def register_matrix(ds, prepared, config, group_overrides, *, exclusions: dict | None = None) -> dict:
    """前処理済みの行列を `ds.stat_matrices[prepared.matrix_id]` に置き、その登録を返す。

    登録は `{matrix_id, prepared, config, group_overrides, exclusions: {samples, features},
    dataset_hash, metadata_revision}`。`config` は前処理の設定（None なら `prepared.config`）、
    `group_overrides` は試料名 → 群（None は空）。`exclusions` を渡さなければ、登録のときの
    `ds.excluded_samples` / `ds.excluded_features`（`dataset_exclude` の除外）を名前の並びで記録する
    ——行列は今の除外で作る前提で、除外が変われば `metadata_revision` が進んで登録ごと古くなる。
    同じ行列 ID は上書きし、古くなった登録はこのとき捨てる。
    別のデータから作った行列（`prepared.dataset_hash` が今の指紋と違う）は ValueError。
    """
    current = dataset_fingerprint(ds)
    recorded = getattr(prepared, "dataset_hash", None)
    if recorded is not None and recorded != current:
        raise ValueError("行列は今のデータセットとは別のデータから作られている（指紋が違う）")
    revision = getattr(ds, "metadata_revision", 0)
    store = getattr(ds, "stat_matrices", None)
    if store is None:
        store = ds.stat_matrices = {}
    for stale in [mid for mid, entry in store.items() if not _entry_current(entry, current, revision)]:
        del store[stale]
    if exclusions is None:
        exclusions = {"samples": sorted(getattr(ds, "excluded_samples", None) or {}),
                      "features": sorted(getattr(ds, "excluded_features", None) or {})}
    else:
        exclusions = {"samples": list(exclusions.get("samples") or []),
                      "features": list(exclusions.get("features") or [])}
    entry = {
        "matrix_id": prepared.matrix_id,
        "prepared": prepared,
        "config": copy.deepcopy(config if config is not None else prepared.config),
        "group_overrides": dict(sorted((group_overrides or {}).items())),
        "exclusions": exclusions,
        "dataset_hash": current,
        "metadata_revision": revision,
    }
    store[prepared.matrix_id] = entry
    return entry


def get_matrix(ds, matrix_id: str) -> dict | None:
    """登録した行列。無い・古い（データの指紋か `metadata_revision` が登録のときと違う）なら None。"""
    entry = (getattr(ds, "stat_matrices", None) or {}).get(matrix_id)
    if entry is None:
        return None
    if not _entry_current(entry, dataset_fingerprint(ds), getattr(ds, "metadata_revision", 0)):
        return None
    return entry


# --- 結果の登録 -------------------------------------------------------------------

def _display_names(ds) -> dict[str, str]:
    """エンジンの試料 ID → 表示名。ID の決め方は `stats.data.from_dataset` と同じで、位置で結ぶ。"""
    names = list(getattr(ds, "sample_names", None) or [])
    rows = getattr(ds, "sample_metadata_rows", None)
    ids = (list(getattr(ds, "sample_ids", None) or []) if rows else []) or names
    out: dict[str, str] = {}
    for sample_id, name in zip(ids, names):
        out.setdefault(sample_id, name)
    return out


def _unit_cell(unit, design: dict) -> str | None:
    """比較に入る個体の cell（因子の水準を ":" で結んだもの。`stats.design` と同じ）。

    モデルの変数（因子・共変量・ブロック）の値が欠けた個体はエンジンが比較から外すので None。
    """
    values = [_model_value(unit, f) for f in design["factors"]]
    if any(v is None for v in values):
        return None
    if any(_model_value(unit, c["name"]) is None for c in design.get("covariates") or []):
        return None
    block = design.get("block")
    if block is not None and _model_value(unit, block, block=True) is None:
        return None
    return ":".join(values)


def _provenance(ds, result_id: str, kind: str, matrix_id: str, dataset_hash: str) -> dict:
    return {"result_id": result_id, "kind": kind, "dataset_hash": dataset_hash,
            "metadata_revision": getattr(ds, "metadata_revision", 0), "parent_ids": [matrix_id]}


def register_compare(ds, run_result: dict, *, log2fc_threshold, units) -> dict:
    """`run_spec_v2` の結果（compare の解析がちょうど 1 つ）を登録し、その記録を返す。

    `units` は同じ行列の `PreparedMatrix.units`。記録の `units` は 個体 ID → {group, subject_id,
    cell（比較に入らない個体は None）, samples（注入の表示名）}。`log2fc_threshold` は判定に効かない
    （図の色分け・書き出しの significant 列だけ）が、ID には含める（同じ結果でも閾値違いは別の記録）。
    compare が 1 つでない・別のデータから出た結果は ValueError。
    """
    analyses = [a for a in (run_result or {}).get("analyses") or [] if a.get("kind") == "compare"]
    if len(analyses) != 1:
        raise ValueError(f"compare の解析がちょうど 1 つの結果だけを登録できる（{len(analyses)} 個）")
    threshold = float(log2fc_threshold)
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError(f"log2fc_threshold は 0 以上の有限の数: {log2fc_threshold!r}")
    analysis = analyses[0]
    engine = analysis["result"]
    matrix_id = run_result["matrix_id"]
    rule_hash = run_result["rule_hash"]
    current = dataset_fingerprint(ds)
    recorded = (run_result.get("input_fingerprint") or {}).get("dataset_hash")
    if recorded is not None and recorded != current:
        raise ValueError("結果は今のデータセットとは別のデータから出ている（指紋が違う）")
    spec = run_result["spec"]
    names = _display_names(ds)
    unit_map = {u.unit_id: {"group": u.group, "subject_id": u.subject_id,
                            "cell": _unit_cell(u, spec["design"]),
                            "samples": [names.get(s, s) for s in u.sample_ids]}
                for u in units}
    result_id = "cmp_" + canonical_hash({"rule_hash": rule_hash, "matrix_id": matrix_id,
                                         "log2fc_threshold": threshold})[:12]
    record = {
        "kind": "compare",
        "result_id": result_id,
        "analysis_id": analysis["analysis_id"],
        "spec": spec,
        "rule_hash": rule_hash,
        "matrix_id": matrix_id,
        "input_fingerprint": run_result.get("input_fingerprint"),
        "preprocess": run_result.get("preprocess"),
        "engine": engine,
        "units": unit_map,
        "log2fc_threshold": threshold,
        "provenance": _provenance(ds, result_id, "compare", matrix_id, current),
    }
    ds.results[result_id] = record
    ds.last_differential = record
    return record


def unit_group_labels(ds, prepared, group_overrides: dict | None = None) -> dict[str, str]:
    """個体 ID → 色分けの群（PCA の表示だけに使う。数字には効かない）。

    `group_overrides`（試料の表示名 → 群。`study_inputs` が解いたもの）があれば、その個体の試料に
    当たった群（1 つなら、複数にまたがれば `(mixed)`）、当たる試料が無い個体と上書きの無いときは
    行列の個体の群（`unit.group`）。どれも無ければ `(none)`。
    """
    overrides = group_overrides or {}
    names = _display_names(ds)
    out: dict[str, str] = {}
    for unit in prepared.units:
        hit = {overrides[names.get(sid, sid)] for sid in unit.sample_ids
               if names.get(sid, sid) in overrides}
        if len(hit) == 1:
            out[unit.unit_id] = next(iter(hit))
        elif hit:
            out[unit.unit_id] = "(mixed)"
        else:
            out[unit.unit_id] = unit.group if unit.group is not None else "(none)"
    return out


def register_pca(ds, pca_result: dict, groups: dict | None) -> dict:
    """`stats.pca.pca` の結果（stats-pca.v1）を登録し、その記録を返す。

    `groups` は 個体 ID → 色分けの群（無ければ None）。PCA の数字には効かないが ID には含める。
    ID の `n_components` と `impute` は結果から読む（成分の数・補完の方法）。
    """
    matrix_id = pca_result["matrix_id"]
    scaling = pca_result["scaling"]
    n_components = len(pca_result.get("explained_variance_ratio") or [])
    impute = (pca_result.get("imputation") or {}).get("method")
    groups = None if groups is None else dict(groups)
    result_id = "pca_" + canonical_hash({"matrix_id": matrix_id, "scaling": scaling,
                                         "n_components": n_components, "impute": impute,
                                         "groups": groups})[:12]
    record = {
        "kind": "stat_pca",
        "result_id": result_id,
        "matrix_id": matrix_id,
        "scaling": scaling,
        "n_components": n_components,
        "impute": impute,
        "groups": groups,
        "engine": pca_result,
        "provenance": _provenance(ds, result_id, "stat_pca", matrix_id, dataset_fingerprint(ds)),
    }
    ds.results[result_id] = record
    ds.last_pca = record
    return record


# --- 2 群への読み替え ---------------------------------------------------------------

def contrast_ids(record: dict, *, two_group_only: bool = False) -> list[str]:
    """`compare` の記録の対比 ID（spec の順）。`two_group_only` なら交互作用を除く。他の記録は []。"""
    if (record or {}).get("kind") != "compare":
        return []
    return [c["id"] for c in record["spec"]["comparisons"]
            if not (two_group_only and c.get("kind") == "interaction")]


def _cell(assign: dict, factors: list[str]) -> str:
    return ":".join(str(assign.get(f, "")) for f in factors)


def as_two_group(record: dict, contrast_id: str | None = None) -> dict:
    """`compare` の記録の 1 つの対比を、v1 の 2 群の結果（`kind="two_group"`）と同じ欄に写す。

    `a` は対比の基準（reference）、`b` は比較（test）で、log2fc は正なら b が高い
    （`export_contract.LOG2FC_SIGN`）。`samples_a` / `samples_b` は比較に入った個体に属する注入の
    表示名、`n_a` / `n_b` は注入の数、`n_units_a` / `n_units_b` は個体の数。`results` は対比の全特徴
    （検定していない行は p・q・log2fc が None）、`volcano` は p が有限の行だけ。`provenance` は元の
    記録のもの（`kind` は `"compare"`）。v1 の `two_group` の記録はそのまま返す。

    対比が 1 つでなく `contrast_id` が無い・未知の ID なら `CONTRAST_REQUIRED`、交互作用なら
    `CONTRAST_NOT_TWO_GROUP`（どちらも details に候補）。比較の記録でなければ `ANALYSIS_RESULT_NOT_FOUND`。
    """
    kind = (record or {}).get("kind")
    if kind == "two_group":
        return record
    if kind != "compare":
        raise DomainError("ANALYSIS_RESULT_NOT_FOUND",
                          "2 群の結果に読み替えられる比較の結果ではありません。", {"kind": kind})
    ids = contrast_ids(record)
    simple = contrast_ids(record, two_group_only=True)
    candidates = {"contrast_ids": ids, "two_group_contrast_ids": simple}
    if contrast_id is None:
        if len(ids) != 1:
            raise DomainError("CONTRAST_REQUIRED",
                              f"この比較には対比が {len(ids)} つあります。contrast_id を指定してください。",
                              candidates)
        contrast_id = ids[0]
    elif contrast_id not in ids:
        raise DomainError("CONTRAST_REQUIRED", f"この比較に対比 {contrast_id!r} はありません。",
                          {"contrast_id": contrast_id, **candidates})
    comp = next(c for c in record["spec"]["comparisons"] if c["id"] == contrast_id)
    if comp.get("kind") == "interaction":
        raise DomainError("CONTRAST_NOT_TWO_GROUP",
                          f"対比 {contrast_id!r} は交互作用で、2 群の結果に読み替えられません。",
                          {"contrast_id": contrast_id, **candidates})

    factors = record["spec"]["design"]["factors"]
    at = dict(comp.get("at") or {})
    ref_cell = _cell({**at, comp["factor"]: comp["reference"]}, factors)
    test_cell = _cell({**at, comp["factor"]: comp["test"]}, factors)
    units = list(record["units"].values())
    units_a = [u for u in units if u.get("cell") == ref_cell]
    units_b = [u for u in units if u.get("cell") == test_cell]
    samples_a = [s for u in units_a for s in u["samples"]]
    samples_b = [s for u in units_b for s in u["samples"]]

    engine = record["engine"]
    q_threshold = engine["verdict_rule"]["q_threshold"]
    log2fc_threshold = record["log2fc_threshold"]
    rows = [{"feature": r["feature_id"], "mean_a": r["mean_reference"], "mean_b": r["mean_test"],
             "log2fc": r["log2fc"], "t": r["statistic"], "p": r["p"], "q": r["q"],
             "ci_low": r["ci_low"], "ci_high": r["ci_high"], "verdict": r["verdict"]}
            for r in engine["rows"] if r["contrast_id"] == contrast_id]
    finite = [r for r in rows if r["p"] is not None and math.isfinite(r["p"])]
    bh = engine["bh_population"]
    return {
        "kind": "two_group",
        "result_id": record["result_id"],
        "contrast_id": contrast_id,
        "matrix_id": record["matrix_id"],
        "rule_hash": record["rule_hash"],
        "a": comp["reference"],
        "b": comp["test"],
        "samples_a": samples_a,
        "samples_b": samples_b,
        "n_a": len(samples_a),
        "n_b": len(samples_b),
        "n_units_a": len(units_a),
        "n_units_b": len(units_b),
        "q_threshold": q_threshold,
        "log2fc_threshold": log2fc_threshold,
        "log_transform": True,
        "contract_version": export_contract.CONTRACT_VERSION,
        "log2fc_sign": export_contract.LOG2FC_SIGN,
        "results": rows,
        "volcano": differential.volcano_data(finite, q_threshold, log2fc_threshold),
        "caveats": [w.get("message") or str(w.get("code"))
                    for w in engine.get("warnings") or [] if w.get("message") or w.get("code")],
        "effect_size_definition": dict(engine["effect_size_definition"]),
        "statistics": {
            "test": engine["test"],
            "ci_level": engine["ci_level"],
            "verdict_rule": dict(engine["verdict_rule"]),
            "bh_population": {"name": bh["name"], "n": (bh.get("n_by_contrast") or {}).get(contrast_id),
                              "q_interpretation": bh["q_interpretation"]},
            "subject_source": engine["subject_source"],
            "n_units_a": len(units_a),
            "n_units_b": len(units_b),
        },
        "provenance": record["provenance"],
    }
