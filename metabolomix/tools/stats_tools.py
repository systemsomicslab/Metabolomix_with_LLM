"""統計エンジン（`metabolomix.analysis.stats`）の対話ツール。

`dataset_preprocess` が前処理（spec §4）を回して行列を作り、`matrix_id` を返す。行列は
`ds.stat_matrices[matrix_id]` に設定・群の上書き・除外・データの指紋と一緒に置く（v1 の `pp_*` には
書かない）。同じ入力と設定は同じ `matrix_id`。後続の比較・PCA はこの ID を名指しできる。

戻り値は要約だけ（工程ごとの名前・状態・件数）。全量は `ds.stat_matrices` に保持する。
1 回の呼び出しは `reports/stats/` に記録を 1 件追記する（書けないときも前処理は成功し、
`record.written=false` で理由を返す）。
"""
from __future__ import annotations

import math
import re

from mcp.types import ToolAnnotations

from metabolomix.analysis import stats_records, stats_results
from metabolomix.analysis.stats import preprocess as stats_preprocess
from metabolomix.analysis.stats import pca as stats_pca
from metabolomix.analysis.stats import spec_v2
from metabolomix.analysis.stats.data import from_dataset
from metabolomix.core import mcp_core
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error, needs_input
from metabolomix.core.serialization import json_payload, json_safe, round_floats, truncate_lists
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_preprocess", "dataset_compare", "dataset_pca"]

#: `argument`（drift_mode）を受けて呼び直せるツール（`needs_input` の retry_tools の意味）。
#: 受けるのは dataset_preprocess だけ。比較・PCA は返った matrix_id を渡して呼び直す。
_RETRY_TOOLS = ["dataset_preprocess"]
#: 問いの欄 → 呼び直すときの引数名。
_DRIFT_ARGUMENT = "drift_mode"
#: details に載せる一覧の上限（長いものは切る）。
_DETAIL_LIMIT = 20
_STEP_LIST_LIMIT = 20


def _config(*, drift_mode="auto", drift_span=0.75, blank_min_ratio=3.0, detection_min_rate=0.5,
            detection_by="group", max_qc_rsd=0.30, max_d_ratio=0.5, normalize="auto",
            collapse_replicates=True, injection_columns=None) -> dict:
    """平らな引数を前処理の設定（`PREPROCESS_DEFAULTS` の形）へ写す。

    足切りの引数が None ならその足切りを行わない（`filters.<name> = None`）。検証はここでは
    しない（`normalize_preprocess_config` が `PREPROCESS_CONFIG_INVALID` で止める）。
    """
    return {
        "drift": {"mode": drift_mode, "span": drift_span},
        "filters": {
            "blank_ratio": None if blank_min_ratio is None else {"min": blank_min_ratio},
            "detection": (None if detection_min_rate is None
                          else {"by": detection_by, "min_rate": detection_min_rate}),
            "qc_rsd": None if max_qc_rsd is None else {"max": max_qc_rsd},
            "d_ratio": None if max_d_ratio is None else {"max": max_d_ratio},
        },
        "normalize": {"method": normalize},
        "collapse_replicates": collapse_replicates,
        "injection_columns": [] if injection_columns is None else injection_columns,
    }


def _cut(obj):
    """どの深さでも `_DETAIL_LIMIT` を超える list を先頭で切る（非破壊）。dict の値の list は隣に
    `<key>_n_more`（切った件数）を置く。"""
    if isinstance(obj, dict):
        out = {}
        for key, value in obj.items():
            if isinstance(value, (list, tuple)) and len(value) > _DETAIL_LIMIT:
                out[key] = [_cut(v) for v in list(value)[:_DETAIL_LIMIT]]
                out[f"{key}_n_more"] = len(value) - _DETAIL_LIMIT
            else:
                out[key] = _cut(value)
        return out
    if isinstance(obj, (list, tuple)):
        return [_cut(v) for v in list(obj)[:_DETAIL_LIMIT]]
    return obj


def _truncate_details(details: dict | None) -> dict | None:
    """エラーの details の長い一覧を切る（入れ子の中の一覧も。件数は `<key>_n_more` に残す）。"""
    if not details:
        return None
    return json_safe(_cut(details))


def _domain_envelope(exc: DomainError) -> str:
    return mztab_error(exc.code, exc.message, _truncate_details(exc.details))


def _compact(value: dict) -> dict:
    """工程の dict を小さな要約にする（入れ子も）。スカラーはそのまま、一覧は件数（`n_<key>`）、
    値がすべて一覧の dict（`excluded_by_role` など）は `{キー: 件数}`、入れ子の dict（`filters` など）は
    再帰する。特徴 ID の一覧（`dropped_ids`）は件数が `dropped` にあるので載せない。"""
    if value and all(isinstance(v, (list, tuple, set, frozenset)) for v in value.values()):
        return {k: len(v) for k, v in value.items()}
    out: dict = {}
    for key, item in value.items():
        if key == "dropped_ids":
            continue
        if isinstance(item, bool) or item is None or isinstance(item, (int, float)):
            out[key] = item
        elif isinstance(item, str):
            if len(item) <= 80:
                out[key] = item
        elif isinstance(item, (list, tuple, set, frozenset)):
            out[f"n_{key}"] = len(item)
        elif isinstance(item, dict):
            if len(item) > 8 and all(isinstance(v, (int, float)) for v in item.values()):
                out[f"n_{key}"] = len(item)
                out[f"sum_{key}"] = sum(item.values())
            else:
                out[key] = _compact(item)
    return out


def _step_summary(step: dict) -> dict:
    """工程 1 つの要約: 名前・状態・理由を先に、あとに `_compact` した件数と足切りごとの内訳。"""
    head = {"step": step.get("step"), "status": step.get("status")}
    if step.get("reason") is not None:
        head["reason"] = step["reason"]
    rest = {k: v for k, v in step.items() if k not in ("step", "status", "reason")}
    return {**head, **_compact(rest)}


def _group_counts(prepared) -> dict:
    counts: dict[str, int] = {}
    for unit in prepared.units:
        label = unit.group if unit.group is not None else "(none)"
        counts[label] = counts.get(label, 0) + 1
    return dict(sorted(counts.items()))


def _warning_summary(warnings: list[dict]) -> list[dict]:
    """注意の要約（コードと、短い欄だけ。試料名の長い一覧は切る）。"""
    return truncate_lists(json_safe(warnings), _STEP_LIST_LIMIT, keys=("samples", "sample_ids"))


def _write_record(holder, ds, prepared, config, inputs, summary) -> dict:
    """記録を書く。書き先が決まらない・書けないときは `{"written": False, "reason": ...}`。"""
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        parents = []
        record_id = (getattr(ds, "quant_state", None) or {}).get("record_id")
        if record_id:
            parents.append(record_id)
        written = stats_records.write_stats_record(
            reports_dir, "preprocess",
            inputs={"dataset_hash": prepared.dataset_hash, "matrix_id": prepared.matrix_id},
            parent_ids=parents,
            settings={"config": config, "group_overrides": inputs["group_overrides"],
                      "exclusions": {"samples": inputs["exclude_samples"],
                                     "features": inputs["exclude_features"]}},
            summary=summary)
    except DomainError as exc:
        return {"written": False, "reason": f"{exc.code}: {exc.message}"}
    except OSError as exc:
        return {"written": False, "reason": f"OSError: {exc}"}
    return {"written": True, "record_id": written["record_id"], "path": written["path"]}


def _preprocess_summary(study, prepared) -> dict:
    """前処理の戻り値の要約（工程・群ごとの個体数・注意。記録の summary にもなる）。"""
    quant = prepared.quant
    groups = _group_counts(prepared)
    warnings = list(prepared.warnings)
    if len(groups) == 1:
        warnings.append({"code": "single_group", "group": next(iter(groups)),
                         "message": "すべての個体が 1 つの群に入っている（比較できない）。"
                                    "group_levels / group_factors の綴りか試料の表の群を確かめてください。"})
    summary = {
        "matrix_id": prepared.matrix_id,
        "n_units": len(prepared.units),
        "n_qc": len(prepared.qc_ids),
        "n_features_before": len(study.feature_ids),
        "n_features": len(prepared.feature_ids),
        "subject_source": prepared.subject_source,
        "groups": groups,
        "quant": ({"normalized": quant.normalized, "method": quant.method, "unit": quant.unit}
                  if quant is not None else None),
        "steps": [_step_summary(s) for s in prepared.steps],
        "warnings": _warning_summary(warnings),
    }
    return truncate_lists(json_safe(round_floats(summary)), _STEP_LIST_LIMIT)


def _store_matrix(holder, ds, study, prepared, config, inputs):
    """行列を登録し、前処理の記録を 1 件書く。`(登録, 要約, record 欄)` を返す。

    記録が書けたときは登録に `record_id` を持たせる（後続の比較の記録がここを親にする）。
    別のデータから作った行列は DomainError(STUDY_DATA_INVALID)。
    """
    try:
        entry = stats_results.register_matrix(
            ds, prepared, config, inputs["group_overrides"],
            exclusions={"samples": inputs["exclude_samples"],
                        "features": inputs["exclude_features"]})
    except ValueError as exc:
        raise DomainError("STUDY_DATA_INVALID", str(exc)) from exc
    summary = _preprocess_summary(study, prepared)
    record = _write_record(holder, ds, prepared, config, inputs, summary)
    if record.get("written"):
        entry["record_id"] = record["record_id"]
    return entry, summary, record


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_preprocess(
    drift_mode: str = "auto",
    drift_span: float = 0.75,
    blank_min_ratio: float | None = 3.0,
    detection_min_rate: float | None = 0.5,
    detection_by: str = "group",
    max_qc_rsd: float | None = 0.30,
    max_d_ratio: float | None = 0.5,
    normalize: str = "auto",
    collapse_replicates: bool = True,
    injection_columns: list[str] | None = None,
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
    dataset_handle: str | None = None,
) -> str:
    """統計エンジンの前処理を回して、前処理済みの行列を作る（行列 ID を返す）。

    工程は 選択 → ドリフト補正 → 品質の足切り → 正規化 → log2 → 反復のまとめ（個体ごとの行列）。
    同じデータ・除外・群・設定は同じ `matrix_id`。比較（dataset_compare）・PCA（dataset_pca）は
    この ID を受ける（省略時は既定の設定）。`dataset_exclude` の除外は自動で反映される。
    `reports/stats/` に記録を 1 件追記する（書けないときは record.written=false。前処理は成功）。

    drift_mode: QC の LOESS ドリフト補正。"auto"（既定。条件が揃えば行う）/ "on" / "off"。
        内部標準などで換算済みの値に QC がある場合は auto だと尋ねる（needs_input の封筒。
        答えを drift_mode に入れて呼び直す）。drift_span: LOESS の窓（0〜1、既定 0.75）。
    blank_min_ratio: 生体試料の平均がブランク平均のこの倍数未満の特徴を除く（既定 3.0）。
    detection_min_rate: 検出率（gap-fill を除く）の下限（既定 0.5）。detection_by: "group"（群ごと）/
        "all"（全体）。検出の別が無いデータでは評価できない旨を steps に出す。
    max_qc_rsd: QC の RSD がこれを超える特徴を除く（既定 0.30）。max_d_ratio: D-ratio の上限（既定 0.5）。
    上の 4 つの足切りは None でその足切りを行わない。
    normalize: "auto"（既定。換算済みなら何もしない）/ "none" / "tic" / "median" / "pqn"。
    collapse_replicates: 同じ個体（subject_id）の反復注入を 1 個体にまとめる（既定 True）。
    injection_columns: 個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う属性）。
    group_levels / group_factors: 試料名（と試料の表の群）の因子トークンで群を解く（group_factors は因子軸の
        リストで優先）。解いた群は行列の設定に記録される。無ければ試料の表の群。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    try:
        config = stats_preprocess.normalize_preprocess_config(_config(
            drift_mode=drift_mode, drift_span=drift_span, blank_min_ratio=blank_min_ratio,
            detection_min_rate=detection_min_rate, detection_by=detection_by,
            max_qc_rsd=max_qc_rsd, max_d_ratio=max_d_ratio, normalize=normalize,
            collapse_replicates=collapse_replicates, injection_columns=injection_columns))
        inputs = stats_results.study_inputs(ds, group_levels=group_levels,
                                            group_factors=group_factors)
        study = from_dataset(ds, **inputs)
        prepared = stats_preprocess.preprocess(study, config)
    except DomainError as exc:
        return _domain_envelope(exc)

    if isinstance(prepared, stats_preprocess.NeedsInput):
        return needs_input(prepared.field, prepared.reason, prepared.options,
                           argument=_DRIFT_ARGUMENT, retry_tools=_RETRY_TOOLS)

    try:
        _entry, summary, record = _store_matrix(holder, ds, study, prepared, config, inputs)
    except DomainError as exc:
        return _domain_envelope(exc)
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "record": record,
        "next": ("dataset_compare（2 群以上の比較）または dataset_pca を、この matrix_id を"
                 "渡して（または省略して既定の前処理で）実行してください"),
    }))


# --- dataset_compare ---------------------------------------------------------------------

#: 戻り値・記録の警告で、長い一覧を切るキー（エンジンの注意が持つ一覧）。
_WARNING_LIST_KEYS = ("samples", "sample_ids", "unit_ids", "examples", "subjects")
_TOP_COLUMNS = ["contrast_id", "feature_id", "name", "log2fc", "ci_low", "ci_high", "q", "verdict"]
#: top の件数の上限（戻り値を肥大させない）。
_TOP_N_MAX = 50
_ID_UNSAFE = re.compile(r"[^A-Za-z0-9_.\-]")


def _bad_request(message: str, **details) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details or None)


def _id_part(name: str, position: int) -> str:
    """対比 ID の群の部分。英数字が 1 つも残らない名前（日本語など）は群の並びの位置（`g<n>`）にする。"""
    part = _ID_UNSAFE.sub("_", name)
    return part if part.strip("_") else f"g{position}"


def _contrast_id(test: str, reference: str, taken: set[str], *, test_pos: int, ref_pos: int) -> str:
    """群名から対比 ID を作る（spec の id は英数字と `_.-` だけ。重なれば番号を足す）。
    名前が ID にできない群は並びの位置で呼ぶ（呼び出しをまたいで安定。実名は対比の要約の
    reference / test に残る）。"""
    base = f"{_id_part(test, test_pos)}_vs_{_id_part(reference, ref_pos)}"
    cid, n = base, 2
    while cid in taken:
        cid, n = f"{base}_{n}", n + 1
    taken.add(cid)
    return cid


def _design_and_comparisons(reference_group, test_group, groups, comparisons, factors):
    """群の指定（ちょうど 1 つ）を spec v2 の `(design.factors, comparisons)` に写す。

    2 群は `reference_group`+`test_group`、3 群以上は `groups`（先頭が基準で各群 − 基準）、
    2 因子などは `comparisons`（spec の形そのまま。`factors` が要る）。それ以外は DATASET_BAD_REQUEST。
    """
    pair = (reference_group is not None, test_group is not None)
    if any(pair) and not all(pair):
        raise _bad_request("reference_group と test_group は対で指定してください。")
    chosen = [name for name, on in (("reference_group+test_group", all(pair)),
                                    ("groups", groups is not None),
                                    ("comparisons", comparisons is not None)) if on]
    if len(chosen) != 1:
        raise _bad_request("群の指定はちょうど 1 つです（reference_group+test_group / groups / "
                           f"comparisons のどれか）。指定された: {chosen or 'なし'}",
                           specified=chosen)
    if comparisons is None and factors is not None:
        raise _bad_request("factors は comparisons（2 因子などの対比）を指定するときだけ使います。")
    if comparisons is not None:
        if not factors:
            raise _bad_request("comparisons を指定するときは factors（モデルの因子）も指定してください。")
        return list(factors), list(comparisons)
    if groups is not None:
        if (not isinstance(groups, list) or len(groups) < 2 or len(set(groups)) != len(groups)
                or any(not isinstance(g, str) or not g for g in groups)):
            raise _bad_request("groups は重複のない 2 つ以上の群名のリスト（先頭が基準）です。",
                               groups=groups)
        ref, tests = groups[0], groups[1:]
    else:
        ref, tests = reference_group, [test_group]
    taken: set[str] = set()
    return ["group"], [{"id": _contrast_id(t, ref, taken, test_pos=i, ref_pos=0), "factor": "group",
                        "reference": ref, "test": t} for i, t in enumerate(tests, start=1)]


def _is_number(text: str) -> bool:
    try:
        return math.isfinite(float(text))
    except (TypeError, ValueError):
        return False


def _covariate_specs(study, covariates) -> list[dict]:
    """共変量の列名 → `{"name", "type"}`。試料の表の列（生体の試料の値）がすべて数値なら numeric、
    さもなくば categorical。どの試料にも値が無い列は DATASET_BAD_REQUEST（候補の列名を details に）。"""
    out = []
    samples = [s for s in study.samples if s.role == "sample" and s.include]
    for name in covariates or []:
        if name in ("group", "batch"):
            values = [getattr(s, name) for s in samples if getattr(s, name)]
        else:
            values = [s.attributes[name] for s in samples if name in s.attributes]
        if not values:
            available = sorted({k for s in samples for k in s.attributes})
            raise _bad_request(f"共変量 {name!r} は試料の表に値がありません。",
                               covariate=name, available=available[:_DETAIL_LIMIT])
        out.append({"name": name,
                    "type": "numeric" if all(_is_number(v) for v in values) else "categorical"})
    return out


def _spec_for(*, factors, comparisons, covariates, block, method, effect_size, bh_population,
              q_threshold, equivalence_fc, allow_confounded, preprocess, inputs) -> dict:
    """引数から組み立てた spec v2（正規化前）。`analyses` は compare 1 つ。結果の `spec` はこれを
    `normalize_spec_v2` したもの（`dataset_analysis_run` に渡せば同じ数字を出し直せる）。"""
    return {
        "schema": spec_v2.ANALYSIS_SCHEMA_V2,
        "input": {"source": "dataset", "exclude_samples": list(inputs["exclude_samples"]),
                  "exclude_features": list(inputs["exclude_features"]),
                  "group_overrides": dict(inputs["group_overrides"])},
        "preprocess": preprocess,
        "design": {"factors": factors, "covariates": covariates, "block": block},
        "comparisons": comparisons,
        "analyses": [{"analysis_id": "compare", "kind": "compare", "test": method,
                      "effect_size": effect_size, "bh_population": bh_population,
                      "verdict": {"q": q_threshold, "equivalence_fc": equivalence_fc},
                      "allow_confounded": allow_confounded}],
        "primary": "compare",
    }


def _fmt(value, spec: str = ".4g") -> str:
    return "" if value is None else format(value, spec)


def _top_table(ds, rows, cids, limit: int) -> str:
    """対比ごとに changed の上位（q の小さい順）を TSV（列名は 1 回）にする。無ければ空文字。"""
    lines = []
    meta = getattr(ds, "feature_metadata", None) or {}
    for cid in cids:
        changed = [r for r in rows if r["contrast_id"] == cid and r["verdict"] == "changed"
                   and r["q"] is not None]
        changed.sort(key=lambda r: (r["q"], r["feature_id"]))
        for r in changed[:limit]:
            name = str((meta.get(r["feature_id"]) or {}).get("name") or "")
            lines.append("\t".join([cid, str(r["feature_id"]),
                                    name.replace("\t", " ").replace("\n", " "),
                                    _fmt(r["log2fc"]), _fmt(r["ci_low"]), _fmt(r["ci_high"]),
                                    _fmt(r["q"], ".3g"), r["verdict"]]))
    return "\n".join(["\t".join(_TOP_COLUMNS), *lines]) if lines else ""


def _contrast_summaries(engine: dict, spec: dict, log2fc_threshold: float) -> list[dict]:
    """対比ごとの要約。`large_changes` は判定が changed で |log2fc| が閾値以上の件数。"""
    q_thr = engine["verdict_rule"]["q_threshold"]
    presence = engine["presence"]
    out = []
    for comp in spec["comparisons"]:
        cid = comp["id"]
        rows = [r for r in engine["rows"] if r["contrast_id"] == cid]
        verdicts = {"changed": 0, "unchanged": 0, "inconclusive": 0}
        for r in rows:
            verdicts[r["verdict"]] += 1
        n_tested = sum(1 for r in rows if r["status"] == "tested")
        large = sum(1 for r in rows if r["verdict"] == "changed" and r["log2fc"] is not None
                    and abs(r["log2fc"]) >= log2fc_threshold)
        has_presence = bool(presence["detection_available"]) and cid in presence["bh_n_by_contrast"]
        out.append({
            "contrast_id": cid,
            "reference": rows[0]["reference"] if rows else None,
            "test": rows[0]["test"] if rows else None,
            "n_tested": n_tested, "n_not_tested": len(rows) - n_tested,
            "verdicts": verdicts, "large_changes": large,
            "presence": {"available": has_presence,
                         "n_q_below": (sum(1 for r in rows if r["presence_q"] is not None
                                           and r["presence_q"] <= q_thr) if has_presence else None)},
        })
    return out


def _omnibus_summary(engine: dict):
    if engine["omnibus"] is None:
        return None
    q_thr = engine["verdict_rule"]["q_threshold"]
    tested = [o for o in engine["omnibus"] if o["f_p"] is not None]
    return {"n_tested": len(tested),
            "n_q_below": sum(1 for o in tested if o["f_q"] is not None and o["f_q"] <= q_thr)}


def _prepare_default(study, config, tool: str):
    """既定の前処理（純粋。登録も記録もしない）。利用者への確認が要るときは `needs_input` の封筒（str）、
    さもなくば `PreparedMatrix`。`tool` は理由に書く呼び直し先（drift_mode の引数を持たないツール）。"""
    prepared = stats_preprocess.preprocess(study, config)
    if isinstance(prepared, stats_preprocess.NeedsInput):
        return needs_input(
            prepared.field,
            f"{prepared.reason} {tool} には drift_mode 引数がありません。"
            f"dataset_preprocess(drift_mode=...) で行列を作り、その matrix_id を "
            f"{tool} に渡して呼び直してください。",
            prepared.options, argument=_DRIFT_ARGUMENT, retry_tools=_RETRY_TOOLS)
    return prepared


def _register_default(holder, ds, study, prepared, config, inputs):
    """既定の前処理の行列を登録する（書き込みを伴う部分）。同じ行列が登録済みなら再利用する。
    `(登録, 前処理の記録 or None)` を返す（記録を書いたのは新しく登録したときだけ）。"""
    entry = stats_results.get_matrix(ds, prepared.matrix_id)
    if entry is not None:
        return entry, None
    entry, _summary, record = _store_matrix(holder, ds, study, prepared, config, inputs)
    return entry, record


def _write_compare_record(holder, entry, record, spec, summary) -> dict:
    """比較の記録を書く（親は行列の前処理の記録）。書けないときは `{"written": False, "reason"}`。"""
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        parents = [entry["record_id"]] if entry.get("record_id") else []
        written = stats_records.write_stats_record(
            reports_dir, "compare",
            inputs={"dataset_hash": entry["dataset_hash"], "matrix_id": record["matrix_id"],
                    "result_id": record["result_id"]},
            parent_ids=parents,
            settings={"spec": spec, "log2fc_threshold": record["log2fc_threshold"]},
            summary=summary)
    except DomainError as exc:
        return {"written": False, "reason": f"{exc.code}: {exc.message}"}
    except OSError as exc:
        return {"written": False, "reason": f"OSError: {exc}"}
    return {"written": True, "record_id": written["record_id"], "path": written["path"]}


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_compare(
    reference_group: str | None = None,
    test_group: str | None = None,
    groups: list[str] | None = None,
    comparisons: list[dict] | None = None,
    factors: list[str] | None = None,
    covariates: list[str] | None = None,
    block: str | None = None,
    matrix_id: str | None = None,
    method: str = "moderated_t",
    effect_size: str = "log2_mean_difference",
    bh_population: str = "matrix",
    q_threshold: float = 0.05,
    equivalence_fc: float = 1.2,
    log2fc_threshold: float = 1.0,
    allow_confounded: bool = False,
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
    top_n: int = 10,
    dataset_handle: str | None = None,
) -> str:
    """2 群・多群・2 因子の比較を統計エンジン（analysis-spec v2）で行い、要約を返す。

    群の指定はちょうど 1 つ: reference_group+test_group（2 群）/ groups（3 群以上。先頭が基準で
    各群 − 基準の対比と omnibus F）/ comparisons（spec v2 の comparisons の形。2 因子・`at`・
    交互作用。factors が要る）。判定（changed / unchanged / inconclusive）は q と CI で決まり、
    log2fc_threshold は判定に使わない（large_changes の数え方・図の色分け・書き出しの significant 列だけ）。
    全行は result_id でセッションに保持し、直近の結果（図・書き出しの既定）にもなる。
    `reports/stats/` に記録を 1 件追記する（書けないときは record.written=false。比較は成功）。

    covariates: 試料の表の列名。値がすべて数値なら numeric、さもなくば categorical。
    block: "subject_id" だけ（同じ個体の反復を対応づける。collapse_replicates=False の行列で意味を持つ）。
    matrix_id: dataset_preprocess が返した行列 ID（その前処理・群・除外のまま比較する）。省略時は既定の
        前処理で行列を作って登録する（QC のある換算済みデータでは drift_mode を尋ねる needs_input。
        dataset_preprocess で行列を作ってその ID を渡す）。古い・無い ID は missing_state。
        ID と group_levels / group_factors の併用は不可（群は行列が決めている）。
    method: "moderated_t"（既定）など（エンジンの検定）。effect_size: "log2_mean_difference"（既定）/
        "arithmetic_mean_ratio"。bh_population: "matrix"（既定）/ "annotated"（"reported" は報告した特徴の
        一覧が要るので dataset_analysis_run の spec v2 だけで使える。ここでは DATASET_BAD_REQUEST）。
    q_threshold / equivalence_fc: 判定の q の閾値（既定 0.05）と同等とみなす倍率（既定 1.2）。
    allow_confounded: 群とバッチが完全交絡でも進める（既定は STATISTIC_CONFOUNDED で止める）。
    top_n: 対比ごとに changed の上位（q の小さい順）を top に何件載せるか（既定 10、上限 50。
        全行は result_id・書き出しで取り出す）。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    implicit_matrix = None
    try:
        if isinstance(top_n, bool) or not isinstance(top_n, int) or not 0 <= top_n <= _TOP_N_MAX:
            raise _bad_request(f"top_n は 0〜{_TOP_N_MAX} の整数です。全行は結果（result_id）に保持して"
                               "あり、書き出し（dataset_export_differential）で取り出せます。",
                               top_n=top_n)
        if bh_population == "reported":
            raise _bad_request('bh_population="reported" はこのツールでは使えません（報告した特徴の一覧を'
                               "渡せない）。dataset_analysis_run の spec v2（reported_feature_ids）で"
                               "指定してください。", bh_population=bh_population)
        if (isinstance(log2fc_threshold, bool) or not isinstance(log2fc_threshold, (int, float))
                or not math.isfinite(log2fc_threshold) or log2fc_threshold < 0):
            raise _bad_request("log2fc_threshold は 0 以上の有限の数です。",
                               log2fc_threshold=log2fc_threshold)
        if block is not None and block != "subject_id":
            raise _bad_request('block は "subject_id" だけを受けます。', block=block)
        design_factors, spec_comparisons = _design_and_comparisons(
            reference_group, test_group, groups, comparisons, factors)

        entry = None
        if matrix_id is not None:
            if group_levels or group_factors:
                raise _bad_request("matrix_id を渡すと群は行列が決めています。group_levels / "
                                   "group_factors は dataset_preprocess に渡してください。")
            entry = stats_results.get_matrix(ds, matrix_id)
            if entry is None:
                return missing_state(
                    "stat_matrix", ["dataset_preprocess"],
                    f"行列 {matrix_id!r} がありません（無い・または除外や試料の表の変更で古くなりました）。"
                    "dataset_preprocess を実行し直して、その matrix_id を渡してください。")
            inputs = {"exclude_samples": list(entry["exclusions"]["samples"]),
                      "exclude_features": list(entry["exclusions"]["features"]),
                      "group_overrides": dict(entry["group_overrides"])}
            preprocess = entry["config"]
        else:
            inputs = stats_results.study_inputs(ds, group_levels=group_levels,
                                                group_factors=group_factors)
            preprocess = stats_preprocess.normalize_preprocess_config(_config())
        study = from_dataset(ds, **inputs)
        spec = spec_v2.normalize_spec_v2(_spec_for(
            factors=design_factors, comparisons=spec_comparisons,
            covariates=_covariate_specs(study, covariates), block=block, method=method,
            effect_size=effect_size, bh_population=bh_population, q_threshold=q_threshold,
            equivalence_fc=equivalence_fc, allow_confounded=allow_confounded,
            preprocess=preprocess, inputs=inputs))

        prepared = (entry["prepared"] if entry is not None
                    else _prepare_default(study, spec["preprocess"], "dataset_compare"))
        if isinstance(prepared, str):
            return prepared

        preprocess_record = None
        if entry is None:
            entry, preprocess_record = _register_default(holder, ds, study, prepared,
                                                         preprocess, inputs)
            implicit_matrix = prepared.matrix_id
        run = spec_v2.run_spec_v2(spec, study, prepared=prepared)
        try:
            record = stats_results.register_compare(ds, run, log2fc_threshold=log2fc_threshold,
                                                    units=prepared.units)
        except ValueError as exc:
            raise DomainError("STUDY_DATA_INVALID", str(exc)) from exc
    except DomainError as exc:
        if implicit_matrix is not None:           # 行列は登録済み。呼び直しにそのまま使える
            exc = DomainError(exc.code, exc.message, {**(exc.details or {}),
                                                      "matrix_id": implicit_matrix})
        return _domain_envelope(exc)

    engine = record["engine"]
    conf = engine["confounding"]
    summary = {
        "result_id": record["result_id"],
        "matrix_id": record["matrix_id"],
        "rule_hash": record["rule_hash"],
        "test": engine["test"],
        "effect_size_definition": {"name": engine["effect_size_definition"]["name"]},
        "ci_level": engine["ci_level"],
        "verdict_rule": engine["verdict_rule"],
        "bh_population": engine["bh_population"],
        "subject_source": engine["subject_source"],
        "group_n": engine["group_n"],
        "prior": engine["prior"],
        "contrasts": _contrast_summaries(engine, spec, record["log2fc_threshold"]),
        "omnibus": _omnibus_summary(engine),
        "confounding": {"batch_status": conf.get("batch_status"),
                        "batch_table": conf.get("batch_table"),
                        "warnings": conf.get("warnings")},
        "warnings": engine["warnings"],
    }
    summary = truncate_lists(json_safe(round_floats(summary)), _STEP_LIST_LIMIT,
                             keys=_WARNING_LIST_KEYS)
    written = _write_compare_record(holder, entry, record, spec, summary)
    if preprocess_record is not None:
        written = {**written, "preprocess_record": preprocess_record}
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "top": _top_table(ds, engine["rows"], [c["id"] for c in spec["comparisons"]], top_n),
        "record": written,
        "next": ("全行は result_id で保持した（直近の結果として書き出し・図の既定になる）。対比が複数の結果は"
                 "書き出し・図で contrast_id を指定する（交互作用の対比は 2 群の形にならない）"),
    }))


# --- dataset_pca -------------------------------------------------------------------------

def _write_pca_record(holder, entry, record, settings, summary) -> dict:
    """PCA の記録を書く（親は行列の前処理の記録）。書けないときは `{"written": False, "reason"}`。"""
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        parents = [entry["record_id"]] if entry.get("record_id") else []
        written = stats_records.write_stats_record(
            reports_dir, "pca",
            inputs={"dataset_hash": entry["dataset_hash"], "matrix_id": record["matrix_id"],
                    "result_id": record["result_id"]},
            parent_ids=parents, settings=settings, summary=summary)
    except DomainError as exc:
        return {"written": False, "reason": f"{exc.code}: {exc.message}"}
    except OSError as exc:
        return {"written": False, "reason": f"OSError: {exc}"}
    return {"written": True, "record_id": written["record_id"], "path": written["path"]}


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=False),
          structured_output=False)
def dataset_pca(
    matrix_id: str | None = None,
    scaling: str = "autoscale",
    n_components: int = 5,
    impute: str = "half_min",
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
    dataset_handle: str | None = None,
) -> str:
    """前処理済みの行列で PCA を行う（統計エンジンの PCA。stats-pca.v1）。

    軸は生体の個体だけで決め、QC は同じ軸へ投影する（軸の計算には使わない）。欠測は PCA のためだけに
    補う（行列そのものは変えない）。値は log2 のまま（log_transform 引数は無い）。全量（スコア・
    ローディング）は戻り値に載せず result_id でセッションに保持し、直近の PCA（図の既定）にもなる。
    `reports/stats/` に記録を 1 件追記する（書けないときは record.written=false。PCA は成功）。

    matrix_id: dataset_preprocess が返した行列 ID。省略時は既定の前処理で行列を作って登録する
        （QC のある換算済みデータでは drift_mode を尋ねる needs_input。dataset_preprocess で行列を
        作ってその ID を渡す）。古い・無い ID は missing_state。
    scaling: "autoscale"（既定。平均 0・SD 1）/ "pareto"（SD の平方根で割る）。
    n_components: 主成分数（既定 5。個体数・特徴数の小さい方が上限）。
    impute: 欠測の補い方 "half_min"（既定）/ "knn" / "none"（欠測があれば PCA_MISSING_VALUES）。
    group_levels / group_factors: **色分けの表示と結果の記録にだけ**使う群（試料名の因子トークン）。
        PCA の数字・行列の選択には効かない（matrix_id と併用できる）。無ければ行列の個体の群。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    implicit_matrix = None
    preprocess_record = None
    try:
        # 群の指定の検証は行列を登録する前に済ませる（不正な指定が行列・記録を残さない）
        label_overrides = stats_results.study_inputs(
            ds, group_levels=group_levels, group_factors=group_factors)["group_overrides"]
        if matrix_id is not None:
            entry = stats_results.get_matrix(ds, matrix_id)
            if entry is None:
                return missing_state(
                    "stat_matrix", ["dataset_preprocess"],
                    f"行列 {matrix_id!r} がありません（無い・または除外や試料の表の変更で古くなりました）。"
                    "dataset_preprocess を実行し直して、その matrix_id を渡してください。")
            prepared = entry["prepared"]
        else:
            config = stats_preprocess.normalize_preprocess_config(_config())
            inputs = stats_results.study_inputs(ds)           # 群は行列を選ばない
            study = from_dataset(ds, **inputs)
            prepared = _prepare_default(study, config, "dataset_pca")
            if isinstance(prepared, str):
                return prepared
            entry, preprocess_record = _register_default(holder, ds, study, prepared, config, inputs)
            implicit_matrix = prepared.matrix_id
        result = stats_pca.pca(prepared, scaling=scaling, n_components=n_components, impute=impute)
        labels = stats_results.unit_group_labels(ds, prepared, label_overrides)
        record = stats_results.register_pca(ds, result, labels)
    except DomainError as exc:
        if implicit_matrix is not None:           # 行列は登録済み。呼び直しにそのまま使える
            exc = DomainError(exc.code, exc.message, {**(exc.details or {}),
                                                      "matrix_id": implicit_matrix})
        return _domain_envelope(exc)

    groups: dict[str, int] = {}
    for label in labels.values():
        groups[label] = groups.get(label, 0) + 1
    summary = truncate_lists(json_safe(round_floats({
        "result_id": record["result_id"],
        "matrix_id": record["matrix_id"],
        "scaling": result["scaling"],
        "log2": result["log2"],
        "imputation": result["imputation"],
        "n_fit": result["n_fit"],
        "n_qc_projected": len(result["qc_ids"]),
        "explained_variance_ratio": result["explained_variance_ratio"],
        "groups": dict(sorted(groups.items())),
    })), _STEP_LIST_LIMIT)
    written = _write_pca_record(
        holder, entry, record,
        {"scaling": scaling, "n_components": n_components, "impute": impute,
         "group_labels": label_overrides}, summary)
    if preprocess_record is not None:
        written = {**written, "preprocess_record": preprocess_record}
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "record": written,
        "loadings_note": ("スコアとローディングの全量（特徴数 × 主成分数）は本要約に非同梱。result_id で"
                          "セッションに保持しており、図は dataset_plot(kind=\"pca\" / \"loadings\") で出す。"
                          "QC は軸の計算に使わず投影だけ（n_qc_projected）。"),
    }))
