"""統計エンジンの結果の登録と、2 群の結果への読み替え（構成見直し C の配線の共通部品）。

対話のツール（前処理・比較・PCA）と、差次的結果を読む既存の利用者（書き出し・根拠バンドル・
特徴検索・volcano）の間に立つ。MCP・セッションには依存しない（`ds` を受け取って書くだけ）。

- 行列: `register_matrix` が前処理済みの行列（`PreparedMatrix`）を設定・群の上書き・除外・指紋と一緒に
  `ds.stat_matrices[matrix_id]` へ置き、`get_matrix` が古ければ None を返す。
- 結果: `register_compare` / `register_pca` が `ds.results[result_id]` と `ds.last_differential` /
  `ds.last_pca` に置く。`result_id` は決定的（同じ計算は同じ鍵を上書きする）。
  来歴は `{result_id, kind: "compare"|"stat_pca", dataset_hash, metadata_revision, parent_ids: [matrix_id]}`
  で、古さは `is_stats_current`（`result_state.is_current` もここへ回す）。
- 読み替え: `as_two_group` が `compare` の記録の 1 つの対比を、利用者が読む 2 群の形（`kind="two_group"`）に写す。
  2×2 の交互作用（差の差）は `allow_interaction=True` を渡した利用者（差次的の書き出し）だけが同じ形で読む。

試料は**表示名**（`ds.sample_names`）で持つ。エンジンの試料 ID（試料の表があれば `ds.sample_ids`）は
位置で表示名へ結ぶ（既存の利用者はすべて表示名で試料を引く）。
"""
from __future__ import annotations

import copy
import itertools
import math
import re

from metabolomix.analysis import differential, export_contract
from metabolomix.analysis.result_state import dataset_fingerprint, is_stats_current
from metabolomix.analysis.stats import preprocess as stats_preprocess
from metabolomix.analysis.stats.compare import _model_value
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "MatrixNotFound",
    "as_two_group",
    "contrast_choices",
    "contrast_ids",
    "p_adjust_method",
    "q_interpretation",
    "get_matrix",
    "sample_ids_by_name",
    "split_result_ref",
    "compare_spec",
    "contrast_id_for",
    "is_stats_current",
    "latest_preprocess_matrix",
    "matrix_inputs",
    "prepare_default",
    "register_compare",
    "register_matrix",
    "register_pca",
    "sample_batches",
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


#: 登録の順番（プロセス内で単調に増える）。`latest_preprocess_matrix` が「最後に作った」行列を選ぶのに使う。
_REGISTRATION_SEQ = itertools.count(1)
#: `dataset_preprocess` が作った行列の `origin`（`matrix_id` を省いた比較・PCA が使う行列の印）。
ORIGIN_PREPROCESS = "dataset_preprocess"


def register_matrix(ds, prepared, config, group_overrides, *, exclusions: dict | None = None,
                    origin: str | None = None) -> dict:
    """前処理済みの行列を `ds.stat_matrices[prepared.matrix_id]` に置き、その登録を返す。

    登録は `{matrix_id, prepared, config, group_overrides, exclusions: {samples, features},
    dataset_hash, metadata_revision, origin, seq}`。`origin` は作り手（`dataset_preprocess` が作った
    行列は `ORIGIN_PREPROCESS`、既定の前処理は ``"default"``）、`seq` は登録の順番。`config` は前処理の設定（None なら `prepared.config`）、
    `group_overrides` は試料名 → 群（None は空）。`exclusions` を渡さなければ、登録のときの
    `ds.excluded_samples` / `ds.excluded_features`（`dataset_exclude` の除外）を名前の並びで記録する
    ——行列は今の除外で作る前提で、除外が変われば `metadata_revision` が進んで登録ごと古くなる。
    同じ行列 ID は上書きし、古くなった登録はこのとき捨てる。
    別のデータから作った行列（`prepared.dataset_hash` が今の指紋と違う）と、データの指紋を持たない行列
    （`stats.external.prepared_from_injections` が作る外部の行列など）は ValueError。
    """
    current = dataset_fingerprint(ds)
    recorded = getattr(prepared, "dataset_hash", None)
    if recorded is None:
        raise ValueError("行列はデータセットの指紋を持たない（外部の値から作った行列は登録できない）")
    if recorded != current:
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
        "origin": origin,
        "seq": next(_REGISTRATION_SEQ),
    }
    store[prepared.matrix_id] = entry
    return entry


def latest_preprocess_matrix(ds) -> dict | None:
    """このデータセット（ハンドル）で `dataset_preprocess` が最後に作った、まだ古くない行列の登録。無ければ None。"""
    current = dataset_fingerprint(ds)
    revision = getattr(ds, "metadata_revision", 0)
    entries = [e for e in (getattr(ds, "stat_matrices", None) or {}).values()
               if e.get("origin") == ORIGIN_PREPROCESS and _entry_current(e, current, revision)]
    return max(entries, key=lambda e: e.get("seq") or 0) if entries else None


def get_matrix(ds, matrix_id: str) -> dict | None:
    """登録した行列。無い・古い（データの指紋か `metadata_revision` が登録のときと違う）なら None。"""
    entry = (getattr(ds, "stat_matrices", None) or {}).get(matrix_id)
    if entry is None:
        return None
    if not _entry_current(entry, dataset_fingerprint(ds), getattr(ds, "metadata_revision", 0)):
        return None
    return entry


class MatrixNotFound(LookupError):
    """名指しされた行列が無い・古い（`get_matrix` が None）。ツールは missing_state にする。"""

    def __init__(self, matrix_id: str):
        super().__init__(matrix_id)
        self.matrix_id = matrix_id


def _entry_inputs(entry: dict) -> dict:
    return {"exclude_samples": list(entry["exclusions"]["samples"]),
            "exclude_features": list(entry["exclusions"]["features"]),
            "group_overrides": dict(entry["group_overrides"])}


def matrix_inputs(ds, matrix_id: str | None = None, *, group_levels=None, group_factors=None,
                  prefer_latest: bool = False):
    """「`matrix_id` か直近の前処理か既定の行列」の入力側（純粋。登録も記録もしない）:
    `(登録 or None, inputs, 前処理の設定)`。

    `matrix_id` があれば登録済みの行列（`entry`）とその除外・群の上書き・設定をそのまま返す
    （無い・古いなら `MatrixNotFound`）。無ければ、`prefer_latest` のときは `latest_preprocess_matrix`
    （このハンドルで `dataset_preprocess` が最後に作った行列）を返す——ただし `group_levels` /
    `group_factors` を渡したときは、その群の上書きがその行列のものと同じときだけ（違う群で比べたいなら
    既定の前処理へ回す）。どちらでもなければ `entry=None` で、今の除外と `group_levels` /
    `group_factors` から `study_inputs` を作り、設定は既定（`PREPROCESS_DEFAULTS` を正規化したもの）。
    `inputs` は `from_dataset(ds, **inputs)` にそのまま渡せる。
    """
    if matrix_id is not None:
        entry = get_matrix(ds, matrix_id)
        if entry is None:
            raise MatrixNotFound(matrix_id)
        return entry, _entry_inputs(entry), entry["config"]
    inputs = study_inputs(ds, group_levels=group_levels, group_factors=group_factors)
    if prefer_latest:
        latest = latest_preprocess_matrix(ds)
        if latest is not None and (not (group_levels or group_factors)
                                   or dict(latest["group_overrides"])
                                   == dict(sorted(inputs["group_overrides"].items()))):
            return latest, _entry_inputs(latest), latest["config"]
    return None, inputs, stats_preprocess.normalize_preprocess_config(None)


def prepare_default(study, config):
    """既定の前処理（純粋。登録も記録もしない）。`PreparedMatrix` を返し、作れなければ `DomainError`。

    利用者には尋ねない: drift=auto が補正を当てなかった行列は `warnings` に `drift_not_applied_by_default`
    （理由と drift_mode="on" での呼び直し）を持つ。"""
    return stats_preprocess.preprocess(study, config)


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


def sample_ids_by_name(ds) -> dict[str, str]:
    """試料の表示名 → エンジンの試料 ID（`register_compare` が `units[*].samples` に写す向きの逆）。"""
    return {name: sample_id for sample_id, name in _display_names(ds).items()}


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


def register_compare(ds, run_result: dict, *, log2fc_threshold, units,
                     make_current: bool = True, origin: str | None = None) -> dict:
    """`run_spec_v2` の結果（compare の解析がちょうど 1 つ）を登録し、その記録を返す。

    `units` は同じ行列の `PreparedMatrix.units`。記録の `units` は 個体 ID → {group, subject_id,
    cell（比較に入らない個体は None）, samples（注入の表示名）}。`log2fc_threshold` は判定に効かない
    （図の色分け・書き出しの significant 列だけ）が、ID には含める（同じ結果でも閾値違いは別の記録）。
    compare が 1 つでない・別のデータから出た結果は ValueError。
    `make_current=False` は `ds.last_differential`（図・書き出しの既定）を動かさない（最初の概観が使う。
    利用者が明示に回した比較を中立の概観で上書きしない）。
    `origin`（例 ``"overview"``）は記録の ``origin`` に残し、特徴検索の既定の選択（`result_ids` 省略）から外す
    印になる。同じ ID の記録が印なしで既にあれば（利用者が同じ比較を先に回した）、印を付けない。
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
        # 多重比較の補正の方法（bh / bonferroni / holm / none）。書き出しのメタ行と根拠バンドルが読む。
        "p_adjust_method": engine.get("p_adjust_method", "bh"),
        "provenance": _provenance(ds, result_id, "compare", matrix_id, current),
    }
    existing = ds.results.get(result_id)
    if origin is not None and not (existing is not None and existing.get("origin") is None):
        record["origin"] = origin
    ds.results[result_id] = record
    if make_current:
        ds.last_differential = record
    return record


def sample_batches(ds) -> dict[str, str | None]:
    """エンジンの試料 ID → バッチ（個体にまとめる前の注入ごと。QC も含む。`from_dataset` と同じ試料の表から）。"""
    from metabolomix.analysis.stats.data import from_dataset

    study = from_dataset(ds, exclude_samples=[], exclude_features=[], group_overrides={})
    return {s.sample_id: s.batch for s in study.samples}


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


def register_pca(ds, pca_result: dict, groups: dict | None, *, make_current: bool = True) -> dict:
    """`stats.pca.pca` の結果（stats-pca.v1）を登録し、その記録を返す。

    `groups` は 個体 ID → 色分けの群（無ければ None）。PCA の数字には効かないが ID には含める。
    ID の `n_components` と `impute` は結果から読む（成分の数・補完の方法）。
    `make_current=False` は `ds.last_pca` を動かさない（最初の概観が使う）。
    """
    matrix_id = pca_result["matrix_id"]
    scaling = pca_result["scaling"]
    n_components = len(pca_result.get("explained_variance_ratio") or [])
    impute = (pca_result.get("imputation") or {}).get("method")
    groups = None if groups is None else dict(groups)
    names = _display_names(ds)
    labels = {i: names[i] for i in [*pca_result.get("unit_ids", []), *pca_result.get("qc_ids", [])]
              if names.get(i) not in (None, i)}              # 図の表示名（ID と違うものだけ）
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
        "labels": labels,
        "engine": pca_result,
        "provenance": _provenance(ds, result_id, "stat_pca", matrix_id, dataset_fingerprint(ds)),
    }
    ds.results[result_id] = record
    if make_current:
        ds.last_pca = record
    return record


# --- 比較の spec の組み立て（dataset_compare と最初の概観が共有する） --------------------------

_ID_UNSAFE = re.compile(r"[^A-Za-z0-9_.\-]")


def _id_part(name: str, position: int) -> str:
    """対比 ID の群の部分。英数字が 1 つも残らない名前（日本語など）は群の並びの位置（`g<n>`）にする。"""
    part = _ID_UNSAFE.sub("_", name)
    return part if part.strip("_") else f"g{position}"


def contrast_id_for(test: str, reference: str, taken: set[str], *, test_pos: int, ref_pos: int) -> str:
    """群名から対比 ID を作る（spec の id は英数字と `_.-` だけ。重なれば番号を足す）。
    名前が ID にできない群は並びの位置で呼ぶ（呼び出しをまたいで安定。実名は対比の要約の
    reference / test に残る）。"""
    base = f"{_id_part(test, test_pos)}_vs_{_id_part(reference, ref_pos)}"
    cid, n = base, 2
    while cid in taken:
        cid, n = f"{base}_{n}", n + 1
    taken.add(cid)
    return cid


def compare_spec(*, factors, comparisons, covariates, block, method, effect_size, bh_population,
                 q_threshold, equivalence_fc, allow_confounded, preprocess, inputs,
                 p_adjust: str = "bh", reported_feature_ids=None) -> dict:
    """引数から組み立てた spec v2（正規化前）。`analyses` は compare 1 つ。結果の `spec` はこれを
    `normalize_spec_v2` したもの（`dataset_analysis_run` に渡せば同じ数字を出し直せる）。
    `reported_feature_ids` は渡されたときだけ書く（`bh_population="reported"` の補正の母集団）。"""
    from metabolomix.analysis.stats import spec_v2

    analysis = {"analysis_id": "compare", "kind": "compare", "test": method,
                "effect_size": effect_size, "bh_population": bh_population, "p_adjust": p_adjust,
                "verdict": {"q": q_threshold, "equivalence_fc": equivalence_fc},
                "allow_confounded": allow_confounded}
    if reported_feature_ids is not None:
        analysis["reported_feature_ids"] = reported_feature_ids
    return {
        "schema": spec_v2.ANALYSIS_SCHEMA_V2,
        "input": {"source": "dataset", "exclude_samples": list(inputs["exclude_samples"]),
                  "exclude_features": list(inputs["exclude_features"]),
                  "group_overrides": dict(inputs["group_overrides"])},
        "preprocess": preprocess,
        "design": {"factors": factors, "covariates": covariates, "block": block},
        "comparisons": comparisons,
        "analyses": [analysis],
        "primary": "compare",
    }


# --- 2 群への読み替え ---------------------------------------------------------------

def contrast_ids(record: dict, *, two_group_only: bool = False) -> list[str]:
    """`compare` の記録の対比 ID（spec の順）。`two_group_only` なら交互作用を除く。他の記録は []。"""
    if (record or {}).get("kind") != "compare":
        return []
    return [c["id"] for c in record["spec"]["comparisons"]
            if not (two_group_only and c.get("kind") == "interaction")]


def _cell(assign: dict, factors: list[str]) -> str:
    return ":".join(str(assign.get(f, "")) for f in factors)


def _interaction_cells(comp: dict, factors: list[str]) -> tuple[list[str], list[str]]:
    """2×2 の交互作用の対比の (基準の cell, 比較の cell)。`stats.design` と同じ向き:
    (a1,b1) と (a0,b0) が +1（比較）、(a0,b1) と (a1,b0) が −1（基準）。"""
    a, b = comp["factors"]
    (a0, a1), (b0, b1) = comp["levels"][a], comp["levels"][b]
    test = [_cell({a: a1, b: b1}, factors), _cell({a: a0, b: b0}, factors)]
    reference = [_cell({a: a0, b: b1}, factors), _cell({a: a1, b: b0}, factors)]
    return reference, test


def _contrast_label(comp: dict) -> str:
    """対比の 1 行の説明（`needs_input` の選択肢に出す。向きは test vs reference）。"""
    if comp.get("kind") == "interaction":
        a, b = comp["factors"]
        (a0, a1), (b0, b1) = comp["levels"][a], comp["levels"][b]
        return (f"interaction {a} x {b}: ({a1} vs {a0} at {b}={b1}) minus "
                f"({a1} vs {a0} at {b}={b0})")
    at = comp.get("at") or {}
    where = f" at {', '.join(f'{k}={v}' for k, v in sorted(at.items()))}" if at else ""
    return f"{comp['test']} vs {comp['reference']} ({comp['factor']}{where})"


def contrast_choices(record: dict, *, allow_interaction: bool = False) -> list[dict]:
    """`compare` の記録の対比の選択肢 `[{"id", "kind", "label"}]`（spec の順）。
    `allow_interaction` が偽なら交互作用を除く（2 群の形に読む利用者の選択肢）。"""
    if (record or {}).get("kind") != "compare":
        return []
    return [{"id": c["id"], "kind": "interaction" if c.get("kind") == "interaction" else "simple",
             "label": _contrast_label(c)}
            for c in record["spec"]["comparisons"]
            if allow_interaction or c.get("kind") != "interaction"]


#: 多重比較の補正の方法（記録の `p_adjust_method`）。無ければ BH（補正の方法を選べる前の記録）。
P_ADJUST_METHODS = ("bh", "bonferroni", "holm", "none")
_DEFAULT_P_ADJUST = "bh"


def _declared_p_adjust(engine: dict, spec: dict, analysis_id) -> str | None:
    """エンジンの結果か spec の解析が宣言した補正の方法（どちらも無ければ None）。"""
    value = (engine or {}).get("p_adjust_method")
    if value is None:
        value = (engine or {}).get("p_adjust")
    if value is None:
        analysis = next((a for a in (spec or {}).get("analyses") or []
                         if a.get("analysis_id") == analysis_id), {})
        value = analysis.get("p_adjust_method", analysis.get("p_adjust"))
    if isinstance(value, dict):                     # {"method": ...} の形も受ける
        value = value.get("method")
    return str(value).lower() if value is not None else None


def p_adjust_method(record: dict) -> str:
    """記録の多重比較の補正の方法（`bh` / `bonferroni` / `holm` / `none`）。宣言が無ければ `bh`。"""
    value = (record or {}).get("p_adjust_method")
    if value is None:
        value = _declared_p_adjust((record or {}).get("engine") or {}, (record or {}).get("spec") or {},
                                   (record or {}).get("analysis_id"))
    return str(value).lower() if value is not None else _DEFAULT_P_ADJUST


#: 補正の方法 → `q` の読み方（報告した特徴だけの母集団 `reported` は方法によらず再現の確認）。
_Q_INTERPRETATION = {"bonferroni": "fwer", "holm": "fwer", "none": "unadjusted"}


def q_interpretation(bh: dict, method: str) -> str:
    """`q` の読み方。BH はエンジンの値（`discovery_fdr` / `reproduction_only`）のまま、ほかの補正は方法から
    決める: `bonferroni` / `holm` は `fwer`、`none` は `unadjusted`（補正しない p）。母集団が `reported` なら
    方法によらず `reproduction_only`。補正の方法を記録する前のエンジンが BH の読み方を書いたままの記録でも、
    メタ行が `q_interpretation = discovery_fdr` と `p_adjust_method = holm` を並べて食い違わないようにする。"""
    if bh.get("name") == "reported":
        return "reproduction_only"
    if method in _Q_INTERPRETATION:
        return _Q_INTERPRETATION[method]
    return bh.get("q_interpretation") or "discovery_fdr"


def split_result_ref(result_id: str, contrast_id: str | None = None) -> tuple[str, str | None]:
    """`result_id` か `result_id#contrast_id`（特徴検索が出す形）を `(記録の ID, 対比 ID)` に割る。

    `#` の後ろと `contrast_id` が両方あって食い違えば `DATASET_BAD_REQUEST`（名指しの対比と引数が衝突）。
    どちらか一方なら、その対比。どちらも無ければ対比は None（呼び出し側が `as_two_group` に任せる）。
    """
    base, sep, named = str(result_id).partition("#")
    if not sep:
        return base, contrast_id
    if contrast_id is not None and contrast_id != named:
        raise DomainError("DATASET_BAD_REQUEST",
                          f"result_id が名指しした対比 {named!r} と contrast_id={contrast_id!r} が食い違っています。",
                          {"result_id": result_id, "contrast_id": contrast_id})
    return base, named


def _detection_driven(row: dict, n_units_ref: int, n_units_test: int) -> bool | None:
    """差が検出の有無から来ているか（`spec_statistics.detection_driven`。個体単位の検出数と個体数）。

    検定できた行で、向き（log2fc の符号）と検出数があるときだけ。ほかは None。
    """
    from metabolomix.analysis import spec_statistics

    log2fc = row.get("log2fc")
    if (row.get("status") != "tested" or not isinstance(log2fc, (int, float))
            or not math.isfinite(log2fc) or log2fc == 0):
        return None
    direction = "test_higher" if log2fc > 0 else "ref_higher"
    return spec_statistics.detection_driven(direction, row.get("detected_reference"), n_units_ref,
                                            row.get("detected_test"), n_units_test)


def as_two_group(record: dict, contrast_id: str | None = None, *,
                 allow_interaction: bool = False) -> dict:
    """`compare` の記録の 1 つの対比を、利用者が読む 2 群の形（`kind="two_group"`）に写す。

    `a` は対比の基準（reference）、`b` は比較（test）で、log2fc は正なら b が高い
    （`export_contract.LOG2FC_SIGN`）。`samples_a` / `samples_b` は比較に入った個体に属する注入の
    表示名、`n_a` / `n_b` は注入の数、`n_units_a` / `n_units_b` は個体の数。`results` は対比の全特徴
    （検定していない行は p・q・log2fc が None）、`volcano` は p が有限の行だけ。各行は記録の検出の有無の
    検定 `presence_p` / `presence_q`（Fisher の正確検定。個体単位）と個体単位の検出数 `detected_units_a` /
    `detected_units_b` を写し、`detection_driven`（差が値の大小でなく検出の有無から来ているか。
    `spec_statistics.detection_driven`、検定できた行だけ）を足す。各行の `significant` は
    `export_contract.row_significant`（記録の閾値に加えて `verdict == "changed"`。R16）で、ここで 1 回だけ
    付ける（書き出し・根拠バンドル・volcano はこの印を読む）。`provenance` は元の記録のもの（`kind` は
    `"compare"`）。`allow_confounded=true`（既定）で完全交絡（全対比をまとめた表か、この対比の表）を越えた比較だけ、
    その写しに `comparison.unadjusted_confounded: true` と説明（`warnings`。`caveats` にも）を足す（記録は変えない）。
    `contrast_kind` は `simple` か `interaction`、`statistics.p_adjust_method` は補正の方法（`p_adjust_method`）、
    `statistics.bh_population.q_interpretation` はその方法での `q` の読み方（`q_interpretation`）。

    交互作用（2×2 の差の差）は `allow_interaction=True` のときだけ読む: `a` / `b` は基準の 2 cell・比較の
    2 cell（`"C:M + T:F"` の形）、log2fc は (a1 vs a0 at b1) − (a1 vs a0 at b0)、`mean_a` / `mean_b` は None
    （cell の平均の和に意味が無い）。`interaction` に因子・水準・定義が入る。

    `contrast_id` を省いたとき、読める対比（`allow_interaction` が偽なら交互作用を除く）が 1 つならそれを使う。
    2 つ以上なら `CONTRAST_REQUIRED`（details の `ambiguous: true` と選択肢 `choices`。ツールはこれを
    `needs_input` にする）、交互作用しか無ければ `CONTRAST_NOT_TWO_GROUP`。未知の ID は `CONTRAST_REQUIRED`、
    交互作用を名指しして `allow_interaction` が偽なら `CONTRAST_NOT_TWO_GROUP`（どれも details に候補）。
    比較の記録でなければ `ANALYSIS_RESULT_NOT_FOUND`。
    """
    kind = (record or {}).get("kind")
    if kind != "compare":
        raise DomainError("ANALYSIS_RESULT_NOT_FOUND",
                          "2 群の結果に読み替えられる比較の結果ではありません。", {"kind": kind})
    ids = contrast_ids(record)
    simple = contrast_ids(record, two_group_only=True)
    candidates = {"contrast_ids": ids, "two_group_contrast_ids": simple}
    usable = ids if allow_interaction else simple
    if contrast_id is None:
        if len(usable) == 1:
            contrast_id = usable[0]
        elif not usable:
            raise DomainError("CONTRAST_NOT_TWO_GROUP",
                              "この比較の対比は交互作用だけで、2 群の結果に読み替えられません"
                              "（交互作用は dataset_export_differential で書き出せます）。",
                              {"contrast_id": ids[0] if len(ids) == 1 else None, **candidates})
        else:
            rid = record.get("result_id")
            raise DomainError("CONTRAST_REQUIRED",
                              f"この比較（result_id={rid}）には対比が {len(usable)} つあります。"
                              "どれを使うか contrast_id で選んでください"
                              "（result_id#contrast_id の形でも渡せます）。",
                              {**candidates, "ambiguous": True, "result_id": rid,
                               "choices": contrast_choices(record, allow_interaction=allow_interaction)})
    elif contrast_id not in ids:
        raise DomainError("CONTRAST_REQUIRED", f"この比較に対比 {contrast_id!r} はありません。",
                          {"contrast_id": contrast_id, **candidates})
    comp = next(c for c in record["spec"]["comparisons"] if c["id"] == contrast_id)
    interaction = comp.get("kind") == "interaction"
    if interaction and not allow_interaction:
        raise DomainError("CONTRAST_NOT_TWO_GROUP",
                          f"対比 {contrast_id!r} は交互作用で、2 群の結果に読み替えられません"
                          "（交互作用は dataset_export_differential で書き出せます）。",
                          {"contrast_id": contrast_id, **candidates})

    factors = record["spec"]["design"]["factors"]
    if interaction:
        ref_cells, test_cells = _interaction_cells(comp, factors)
        label_a, label_b = " + ".join(ref_cells), " + ".join(test_cells)
    else:
        at = dict(comp.get("at") or {})
        ref_cells = [_cell({**at, comp["factor"]: comp["reference"]}, factors)]
        test_cells = [_cell({**at, comp["factor"]: comp["test"]}, factors)]
        label_a, label_b = comp["reference"], comp["test"]
    units = list(record["units"].values())
    units_a = [u for u in units if u.get("cell") in ref_cells]
    units_b = [u for u in units if u.get("cell") in test_cells]
    samples_a = [s for u in units_a for s in u["samples"]]
    samples_b = [s for u in units_b for s in u["samples"]]

    engine = record["engine"]
    q_threshold = engine["verdict_rule"]["q_threshold"]
    log2fc_threshold = record["log2fc_threshold"]
    rows = [{"feature": r["feature_id"], "mean_a": r["mean_reference"], "mean_b": r["mean_test"],
             "log2fc": r["log2fc"], "t": r["statistic"], "p": r["p"], "q": r["q"],
             "ci_low": r["ci_low"], "ci_high": r["ci_high"], "verdict": r["verdict"],
             # 検出の有無（個体単位）。記録の値をそのまま写す（検出状態が無ければ None）
             "presence_p": r.get("presence_p"), "presence_q": r.get("presence_q"),
             "detected_units_a": r.get("detected_reference"),
             "detected_units_b": r.get("detected_test"),
             "detection_driven": _detection_driven(r, len(units_a), len(units_b))}
            for r in engine["rows"] if r["contrast_id"] == contrast_id]
    for row in rows:
        # 有意の印はここで 1 回だけ付ける（R16: 閾値に加えて判定が changed）。書き出し・根拠バンドル・
        # volcano・新しい知見の候補はこの印を読み、閾値から計算し直さない。
        row["significant"] = export_contract.row_significant(
            row, q_threshold=q_threshold, log2fc_threshold=log2fc_threshold)
    finite = [r for r in rows if r["p"] is not None and math.isfinite(r["p"])]
    bh = engine["bh_population"]
    method = p_adjust_method(record)
    caveats = [w.get("message") or str(w.get("code"))
               for w in engine.get("warnings") or [] if w.get("message") or w.get("code")]
    provenance = record["provenance"]
    trace = _unadjusted_confounded(record, label_a, label_b, contrast_id, (*ref_cells, *test_cells))
    if trace is not None:
        # 書き出しのメタ行・volcano・図の但し書きはここを読む。
        # 記録の来歴は書き換えず、読み替えた結果の来歴にだけ足す。
        provenance = {**provenance, "comparison": trace,
                      "warnings": [*(provenance.get("warnings") or []), UNADJUSTED_CONFOUNDED_EXPLANATION]}
        caveats.append(UNADJUSTED_CONFOUNDED_EXPLANATION)
    view = {
        "kind": "two_group",
        "result_id": record["result_id"],
        "contrast_id": contrast_id,
        "contrast_kind": "interaction" if interaction else "simple",
        "matrix_id": record["matrix_id"],
        "rule_hash": record["rule_hash"],
        "a": label_a,
        "b": label_b,
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
        "caveats": caveats,
        "effect_size_definition": dict(engine["effect_size_definition"]),
        "statistics": {
            "test": engine["test"],
            "ci_level": engine["ci_level"],
            "verdict_rule": dict(engine["verdict_rule"]),
            "bh_population": {"name": bh["name"], "n": (bh.get("n_by_contrast") or {}).get(contrast_id),
                              "q_interpretation": q_interpretation(bh, method)},
            "p_adjust_method": method,
            "subject_source": engine["subject_source"],
            "n_units_a": len(units_a),
            "n_units_b": len(units_b),
        },
        "provenance": provenance,
    }
    if interaction:
        a, b = comp["factors"]
        view["interaction"] = {
            "factors": [a, b], "levels": {a: list(comp["levels"][a]), b: list(comp["levels"][b])},
            "reference_cells": ref_cells, "test_cells": test_cells,
            "definition": INTERACTION_DEFINITION.format(
                a=a, b=b, a0=comp["levels"][a][0], a1=comp["levels"][a][1],
                b0=comp["levels"][b][0], b1=comp["levels"][b][1])}
    return view


#: 交互作用の対比の log2fc の定義（書き出しのメタ行 `interaction_definition`。`=` とタブを含めない）。
INTERACTION_DEFINITION = ("difference of differences: log2fc is ({a1} vs {a0} at {b} {b1}) minus "
                          "({a1} vs {a0} at {b} {b0}); mean_a and mean_b are not defined (blank)")


#: 完全交絡を未調整のまま越えた比較の来歴の `warnings` と `caveats` に書く説明
#: （書き出しの `unadjusted_confounded_note` は来歴の warnings から「未調整」を含む行を拾う）。
UNADJUSTED_CONFOUNDED_EXPLANATION = (
    "群とバッチが完全交絡の比較を、allow_confounded=true（既定。false なら止める）のまま未調整で"
    "解析しました（群差とバッチ差を区別できない。図・TSV・レポートにこの旨を残すこと）。")


def _unadjusted_confounded(record: dict, reference: str, test: str, contrast_id: str,
                           cells=()) -> dict | None:
    """`allow_confounded=true`（既定）で完全交絡を越えて走った対比なら、来歴の `comparison` に足す
    痕跡（`unadjusted_confounded: true`）。そうでなければ None（何も足さない）。

    完全交絡はエンジンの交絡の表で、全対比をまとめた表（`engine.confounding.batch_status == "complete"`）か、
    この対比の群だけで見た表（`engine.confounding.by_contrast[contrast_id].confounded is True`。多群でまとめた表が
    部分交絡でも、一部の対比だけが完全交絡になりうる）のどちらか。許可は記録の spec のその比較の解析の
    `allow_confounded`（spec v2 の既定は true）。部分的な交絡・評価できない交絡には印を付けない（警告は `caveats`）。
    痕跡の `confounding` は、この対比の状態（`batch_status`）と、その cell（2 群なら 2 つ、交互作用なら 4 つ）に
    絞った群 × バッチの表、全対比をまとめた表の状態（`pooled_batch_status`）。`reference` / `test` は
    写しの `a` / `b` の表示（交互作用なら cell を `+` でつないだもの）。
    """
    confounding = record["engine"].get("confounding") or {}
    analysis = next((a for a in record["spec"].get("analyses") or []
                     if a.get("analysis_id") == record.get("analysis_id")), {})
    allowed = bool(analysis.get("allow_confounded", True))
    pooled = confounding.get("batch_status")
    own = (confounding.get("by_contrast") or {}).get(contrast_id) or {}
    if not (allowed and (pooled == "complete" or own.get("confounded") is True)):
        return None
    table = confounding.get("batch_table") or {}
    if cells and all(c in table for c in cells):
        table = {c: table[c] for c in cells}
    return {"comparison_id": contrast_id,
            "reference_group": reference,
            "test_group": test,
            "confounding": {"batch_status": own.get("batch_status") or pooled,
                            "batch_table": copy.deepcopy(table),
                            "pooled_batch_status": pooled},
            "allow_confounded": allowed,
            "unadjusted_confounded": True}
