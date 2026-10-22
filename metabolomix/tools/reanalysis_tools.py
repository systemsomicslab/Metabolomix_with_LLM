"""再解析の spec で呼ぶ MCP ツール（spec 2026-10-22 §4）。

**再解析で使う**（判定の数字を、その場で書くコードではなく版のあるテスト済みの実装から
出すため）。数値の実装は `metabolomix/analysis/` の純ロジックにあり、ここは spec を
正規化して呼び、戻り値を小さくするだけ。どのツールも `dataset_handle` を受け、省略時は
対話のセッションのデータセットを読む（書き換えない）。

大きな結果は `output_path` に JSON（キーの順を固定・改行なし）で書き、要約だけを返す。
同じ spec・同じ入力なら同じファイルになる。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.core.atomic_io import DomainError, atomic_write_bytes
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload, round_floats
from metabolomix.tools._holder import resolve_holder

__all__ = ["dataset_analysis_run"]

#: output_path なしで行を返す上限（解析ごと）。超えれば output_path を求める。
MAX_INLINE_ROWS = 50

ROW_COLUMNS = ("feature_id", "status", "reason", "group_ref", "group_test", "n_ref", "n_test",
               "mean_ref", "mean_test", "fc_test_over_ref", "log2fc_test_over_ref",
               "statistic", "p", "q", "beta_log2", "direction", "meets_thresholds")
_SIGNIFICANT = frozenset({"mean_ref", "mean_test", "statistic", "p", "q", "beta_log2"})
_DECIMAL4 = frozenset({"fc_test_over_ref", "log2fc_test_over_ref"})


def _cell(column: str, value):
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and math.isfinite(value):
        if column in _SIGNIFICANT:
            return str(float(f"{value:.4g}"))
        if column in _DECIMAL4:
            rounded = round(value, 4)
            return str(0.0 if rounded == 0 else rounded)
        # そのほかの float は有効 6 桁（既定 repr の 17 桁を返さない）
        return str(float(f"{value:.6g}"))
    return str(value)


def rows_table(rows: list[dict], columns=ROW_COLUMNS) -> str:
    from metabolomix.arf2.reader import format_spots_as_table

    cells = [{c: _cell(c, row.get(c)) for c in columns} for row in rows]
    return format_spots_as_table(cells, columns=list(columns)) if cells else ""


def write_result(output_path: str, payload: dict) -> str:
    """結果の JSON を書く（キーの順を固定・改行なし・NaN を許さない）。書いた絶対パスを返す。"""
    path = Path(output_path)
    data = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                      allow_nan=False)
    atomic_write_bytes(path, data.encode("utf-8"))
    return str(path.resolve())


def _dataset_or_error(holder, needed: bool):
    """(ds, error)。データセットが要るのに無ければ missing_state の封筒。"""
    ds = getattr(holder, "dataset", None)
    if needed and ds is None:
        return None, missing_state(
            "dataset", ["dataset_load", "dataset_handle_open"],
            "データセットがありません。dataset_handle_open で開いたハンドルを渡すか、"
            "dataset_load を実行してください。")
    return ds, None


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                      idempotentHint=True),
          structured_output=False)
def dataset_analysis_run(spec: dict, dataset_handle: str | None = None,
                         output_path: str | None = None) -> str:
    """論文の方法を写した spec（analysis-spec.v1）で、行列を作って 2 群を検定する（再解析で使う）。

    spec の形と欄の意味は `lipidmix://docs/output-format/reanalysis` §2。要点:
    `comparison`（`reference` / `test` の群の名前。比較名は `<test>_vs_<ref>`）、
    `analyses`（名前つきの解析を並べる——感度解析は 1 つの spec に並べる。各解析は
    `matrix`（matrix-recipe.v2: 段の順序・検出率・KNN・PQN・LOESS・RSD）と
    `test`（`student` / `welch` / `mannwhitney` / `ols`）・`fold_change`・`correction`・
    `bh_population`・`thresholds`）、入力は `input`（データセットか外部の行列 TSV）、
    群は `groups`（実験情報シートか TSV）、`covariates`（OLS の共変量の TSV）、
    `aggregate`（クラスの和・2 特徴の比を 1 本の特徴として検定）。

    返すのは数字と向き（`direction`）と閾値を機械的に当てた `meets_thresholds` だけで、
    再現したかの判定は返さない。結果には正規化した spec・`rule_hash`・`input_fingerprint`
    が入り、同じ spec・同じ入力なら同じ JSON になる。

    output_path: 結果の JSON をここへ書き、要約だけを返す。省略時は解析ごとの行が
        50 以下のときだけ行を TSV で返す（超えれば `ANALYSIS_OUTPUT_REQUIRED`）。
    対話のセッション・ハンドルの状態は変えない（読むだけ）。
    """
    from metabolomix.analysis import analysis_spec

    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    try:
        normalized = analysis_spec.normalize_spec(spec)
        ds, error = _dataset_or_error(holder, normalized["input"]["source"] == "dataset"
                                      or normalized["aggregate"] is not None
                                      and "select" in (normalized["aggregate"] or {}))
        if error:
            return error
        result = analysis_spec.run_spec(normalized, ds)
        too_many = [a["analysis_id"] for a in result["analyses"]
                    if len(a["rows"]) > MAX_INLINE_ROWS]
        if output_path is None and too_many:
            raise DomainError(
                "ANALYSIS_OUTPUT_REQUIRED",
                f"報告する行が {MAX_INLINE_ROWS} を超える解析がある（{too_many}）。output_path を"
                "指定して JSON に書くか、features.feature_ids で絞ってください。",
                {"analyses": too_many})
        written = write_result(output_path, result) if output_path else None
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    summary = []
    for analysis in result["analyses"]:
        rows = analysis["rows"]
        item = {
            "analysis_id": analysis["analysis_id"],
            "n_features_in": analysis["matrix"]["n_features_in"],
            "n_features_out": analysis["matrix"]["n_features_out"],
            "bh_population_n": analysis["bh_population_n"],
            "n_rows": len(rows),
            "n_tested": sum(1 for r in rows if r["status"] == "tested"),
            "n_meets_thresholds": sum(1 for r in rows if r["meets_thresholds"]),
        }
        if "aggregate" in analysis:
            item["aggregate"] = analysis["aggregate"]
        if "covariate_excluded_samples" in analysis:
            item["covariate_excluded_samples"] = analysis["covariate_excluded_samples"]
        if written is None:
            item["rows"] = rows_table(rows)
        summary.append(item)
    return json_payload({
        "status": "success",
        "schema": result["schema"],
        "rule_hash": result["rule_hash"],
        "input_fingerprint": result["input_fingerprint"],
        "comparison": result["comparison"],
        "group_ref": result["group_ref"],
        "group_test": result["group_test"],
        **({"output_path": written} if written else {}),
        "analyses": summary,
        "note": "判定（再現したか）は返さない。direction は group_test から見た向き。",
    })
