"""統計エンジン（`metabolomix.analysis.stats`）の対話ツール。

`dataset_preprocess` が前処理（spec §4）を回して行列を作り、`matrix_id` を返す。行列は
`ds.stat_matrices[matrix_id]` に設定・群の上書き・除外・データの指紋と一緒に置く（v1 の `pp_*` には
書かない）。同じ入力と設定は同じ `matrix_id`。後続の比較・PCA はこの ID を名指しできる。

戻り値は要約だけ（工程ごとの名前・状態・件数）。全量は `ds.stat_matrices` に保持する。
1 回の呼び出しは `reports/stats/` に記録を 1 件追記する（書けないときも前処理は成功し、
`record.written=false` で理由を返す）。
"""
from __future__ import annotations

from mcp.types import ToolAnnotations

from metabolomix.analysis import stats_records, stats_results
from metabolomix.analysis.stats import preprocess as stats_preprocess
from metabolomix.analysis.stats.data import from_dataset
from metabolomix.core import mcp_core
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error, needs_input
from metabolomix.core.serialization import json_payload, json_safe, round_floats, truncate_lists
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_preprocess"]

#: 答えを入れて呼び直せるツール（`preprocess.drift.mode` の問いを受ける）。
_RETRY_TOOLS = ["dataset_preprocess", "dataset_compare", "dataset_pca"]
#: 問いの欄 → 呼び直すときの引数名。
_DRIFT_ARGUMENT = "drift_mode"
#: details に載せる一覧の上限（長いものは切る）。
_DETAIL_LIMIT = 20
_STEP_LIST_LIMIT = 10


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


def _truncate_details(details: dict | None) -> dict | None:
    """エラーの details の長い一覧を切る（件数は `<key>_n_more` に残す）。"""
    if not details:
        return None
    keys = {k for k, v in details.items() if isinstance(v, list)}
    return json_safe(truncate_lists(details, _DETAIL_LIMIT, keys=keys or ("__none__",)))


def _domain_envelope(exc: DomainError) -> str:
    return mztab_error(exc.code, exc.message, _truncate_details(exc.details))


def _step_summary(step: dict) -> dict:
    """工程 1 つの要約: 名前・状態・理由と、件数（スカラー値と、一覧の長さ）。特徴ごとの一覧は載せない。"""
    out: dict = {"step": step.get("step"), "status": step.get("status")}
    if step.get("reason") is not None:
        out["reason"] = step["reason"]
    for key, value in step.items():
        if key in ("step", "status", "reason"):
            continue
        if isinstance(value, bool) or value is None or isinstance(value, (int, float, str)):
            if not isinstance(value, str) or len(value) <= 80:
                out[key] = value
        elif isinstance(value, (list, tuple, set, frozenset)):
            out[f"n_{key}"] = len(value)
        elif isinstance(value, dict) and all(isinstance(v, (int, float)) for v in value.values()):
            if len(value) <= 8:
                out[key] = value                           # 小さな件数の表（欠測セルの数など）はそのまま
            else:
                out[f"n_{key}"] = len(value)
                out[f"sum_{key}"] = sum(value.values())
    return out


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
    group_levels / group_factors: 試料名・Class ID の因子トークンで群を解く（group_factors は因子軸の
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
        stats_results.register_matrix(ds, prepared, config, inputs["group_overrides"],
                                      exclusions={"samples": inputs["exclude_samples"],
                                                  "features": inputs["exclude_features"]})
    except ValueError as exc:
        return mztab_error("STUDY_DATA_INVALID", str(exc))

    quant = prepared.quant
    summary = {
        "matrix_id": prepared.matrix_id,
        "n_units": len(prepared.units),
        "n_qc": len(prepared.qc_ids),
        "n_features_before": len(study.feature_ids),
        "n_features": len(prepared.feature_ids),
        "subject_source": prepared.subject_source,
        "quant": ({"normalized": quant.normalized, "method": quant.method, "unit": quant.unit}
                  if quant is not None else None),
        "steps": [_step_summary(s) for s in prepared.steps],
        "warnings": _warning_summary(prepared.warnings),
    }
    summary = truncate_lists(json_safe(round_floats(summary)), _STEP_LIST_LIMIT)
    record = _write_record(holder, ds, prepared, config, inputs, summary)
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "record": record,
        "next": ("dataset_compare（2 群以上の比較）または dataset_pca を、この matrix_id を"
                 "渡して（または省略して既定の前処理で）実行してください"),
    }))
