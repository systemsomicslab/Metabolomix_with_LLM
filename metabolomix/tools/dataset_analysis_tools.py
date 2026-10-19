"""DatasetState に対する解析 MCP ツール。

dataset_load → dataset_preprocess → dataset_pca / dataset_differential
→ dataset_export_differential の順に実行する。ARF 経路
（arf_preprocess → arf_pca_preprocessed → arf_differential →
arf_export_differential）と同じ純関数を共有しており、同じ入力からは
同じ数字が出る。dataset_feature_search は読み取り専用で、特徴の検索・目的の物質の
判定に差次的結果の統計を添えて返す（判定の純ロジックは analysis/feature_search.py）。

戻り値には要約だけを載せ、全量（PCA の loadings・差次的の results/volcano）は
session.dataset に保持する（CLAUDE.md の戻り値肥大禁止）。
"""
from __future__ import annotations

import math
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.analysis import export_contract
from metabolomix.analysis.dataset_analysis import PreconditionError
from metabolomix.core import session_state
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error
from metabolomix.core.serialization import json_payload

__all__ = [
    "dataset_preprocess",
    "dataset_pca",
    "dataset_differential",
    "dataset_export_differential",
    "dataset_set_sample_metadata",
    "dataset_statistic",
    "dataset_build_matrix",
    "dataset_feature_search",
]

# 欠けている状態 → それを作れるツール。missing_state の required_tools になる。
# **自分自身は入れない**（クライアントが同じ呼び出しを繰り返すループになる）。
_RECOVERY_TOOLS = {
    "dataset": ["dataset_load"],
    "dataset_preprocessed": ["dataset_preprocess"],
    "dataset_differential_result": ["dataset_differential"],
    # v2 の解析行列は pipeline（preprocess / qc_processed）が作る。単体ツールの
    # dataset_build_matrix でも作れるが、そちらは binding と証拠を持たない分
    # できることが狭い（内部標準比は作れず、検出率filterは評価不能になる）。
    # 先頭ほど優先＝上流の証拠ごと揃う pipeline を先に案内する。
    "analysis_matrix": ["pipeline_run", "pipeline_status", "dataset_load",
                        "dataset_build_matrix"],
    # 対象feature↔内部標準の対応付け。pipeline の resolve_feature_bindings が作る。
    "feature_bindings": ["pipeline_run", "pipeline_status"],
}


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_preprocess(
    normalize: str = "none",
    blank_min_fold: float | None = None,
    drift_correct: bool = False,
    max_qc_rsd: float | None = None,
    impute: str = "half_min",
    min_detection_rate: float = 0.0,
) -> str:
    """DatasetState の定量行列に前処理レシピを適用する。

    引数は arf_preprocess と同一。以降の dataset_pca / dataset_differential は
    ここで作った前処理済み行列を消費する。

    normalize: "none"（既定）/ "tic"（行総和）/ "median"（行中央値）/ "pqn"。
    blank_min_fold: 生体試料平均がブランク平均のこの倍数未満の特徴量を背景として
        除去する（例 3.0）。None（既定）でブランク除去なし。
    drift_correct: QC 注入順ドリフト補正。注入順は mzTab-M の
        `assay[N]-custom[...]` の injection sequence label（MS:4000089）から読む。
        それを持たない mzTab-M では実施できず、caveat で報告する。
    max_qc_rsd: QC 群の RSD がこの値を超える特徴量を除去する（例 0.30）。
    impute: "half_min"（既定）/ "knn" / "column_mean" / "none"。
    min_detection_rate: 実検出率（gap-fill を除く）による特徴量の足切り 0.0-1.0
        （既定 0.0=無効）。**検出状態は mzTab-M 単体には無い**ので、隣接する `.arf`
        から取り込めた場合にだけ使える。取り込めていない状態で 0 より大きい値を
        渡すと引数エラーを返す（黙って未検出 0 件として通さない）。
        取り込み状況は dataset_status の `detection` を見る。

    成功すると session.dataset.pp_matrix に前処理済み行列が設定される。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    from metabolomix.analysis.dataset_service import preprocess_dataset

    recipe = {
        "normalize": normalize,
        "blank_min_fold": blank_min_fold,
        "drift_correct": drift_correct,
        "max_qc_rsd": max_qc_rsd,
        "impute": impute,
        "min_detection_rate": min_detection_rate,
    }
    try:
        # 状態の差し替え（成功時のみの反映・派生結果の無効化・来歴の付与）は
        # サービス側が持つ。ツールはそれを呼んで要約を返すだけにする——
        # ここで直接 ds を書き換えると、pipeline ワーカー経由の実行と規則がずれる。
        report = preprocess_dataset(ds, recipe)
    except PreconditionError as exc:
        return _from_precondition(exc)

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
        "next": "dataset_pca または dataset_differential を実行してください",
        **({"detection": report["detection"]} if "detection" in report else {}),
        **({"detection_filter": report["detection_filter"]}
           if "detection_filter" in report else {}),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_pca(n_components: int = 5, log_transform: bool = False) -> str:
    """前処理済み DatasetState に PCA を実行する。

    dataset_preprocess を先に実行しておくこと。
    n_components: 主成分数（既定 5。サンプル数・特徴量数の小さい方で上限が決まる）。
    log_transform: 標準化の前に log10 変換を適用する（既定 False）。

    ローディング全量は戻り値に載せず session.dataset.last_pca に保持する。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    from metabolomix.analysis.dataset_service import pca_dataset

    try:
        result = pca_dataset(ds, n_components=n_components, log_transform=log_transform)
    except PreconditionError as exc:
        return _from_precondition(exc)

    payload = {k: v for k, v in result.items()
               if k not in ("loadings", "singular_values", "provenance")}
    payload["status"] = "success"
    payload["result_id"] = result["provenance"]["result_id"]
    payload["loadings_note"] = (
        "ローディング全量（特徴量数 × 主成分数）は本要約に非同梱。"
        "セッションに保持しており、寄与特徴量が必要になったら別途取得します。")
    return json_payload(payload)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_differential(
    group_a: list[str],
    group_b: list[str],
    q_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    log_transform: bool = True,
    group_a_label: str = "group_a",
    group_b_label: str = "group_b",
) -> str:
    """前処理済み DatasetState で 2 群の差次的解析（Welch t 検定 + BH-FDR）を実行する。

    group_a / group_b: サンプル名のリスト。使える名前と役割は dataset_status の
        `samples`（name/role の TSV）で確認する。前処理済み行列に無い
        名前、QC/ブランクは除外し、caveat で名指しする。
    q_threshold: BH-FDR 補正後の有意水準（既定 0.05）。
    log2fc_threshold: この絶対値以上の log2FC を有意として数える（既定 1.0）。
    log_transform: log2(x + 1) 空間で検定する（既定 True。MS 強度は対数正規に近い）。

    **log2FC は正なら group_b が高い（上昇）。** group_a が基準（対照）。
    全特徴量の結果と volcano 点列は戻り値に載せず session.dataset.last_differential
    に保持する。エクスポートは dataset_export_differential を使う。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    from metabolomix.analysis.dataset_service import compare_dataset

    from metabolomix.core.atomic_io import DomainError

    try:
        result = compare_dataset(
            ds,
            group_a,
            group_b,
            q_threshold=q_threshold,
            log2fc_threshold=log2fc_threshold,
            log_transform=log_transform,
            group_a_label=group_a_label,
            group_b_label=group_b_label,
        )
    except PreconditionError as exc:
        return _from_precondition(exc)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    payload = {k: v for k, v in result.items()
               if k not in ("results", "volcano", "provenance")}
    payload["status"] = "success"
    payload["result_id"] = result["provenance"]["result_id"]
    payload["differential_contract_version"] = result["contract_version"]
    payload["results_note"] = (
        "全特徴の結果と volcano 点列は本要約に非同梱（セッション保持）。"
        "InChIKey 付きの全行が必要なら dataset_export_differential を実行してください。")
    return json_payload(payload)


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def dataset_export_differential(output_path: str,
                                result_id: str | None = None,
                                apply_curation: bool = True) -> str:
    """指定した差次的結果を InChIKey 付きの 1 ファイルへ書き出す。

    先に dataset_preprocess → dataset_differential を実行しておくこと。
    出力は arf_export_differential と**同一の契約**（15 列 + contract_version
    メタ行）なので、下流のパスウェイ解析にそのまま渡せる。

    apply_curation（既定 true）: curation_submit で wrong を付けたスポットを
    同定なしとして扱い、出力から外す。フラグがあればメタ行 # curation = ... で
    適用状況を宣言する。フラグが無ければ出力は変わらない。成功 payload の
    `curation` は `{state, wrong_excluded, suspect, orphaned}`（arf_export_differential と
    同じ形）。assign（curation_submit の候補付け）は同定を置き換え（name_source /
    inchikey_source = curation）、redundant は除外する。どちらかがあれば `curation` に
    `assigned` と `redundant_excluded` が足される。
    フラグ記録に読めない行があれば書き出さずに `CURATION_FLAGS_INVALID`。

    result_id: 書き出す結果を名指しする（省略時は直近の差次的結果）。
        **前処理をやり直した後の古い結果は書き出しません**（`STALE_ANALYSIS_RESULT`）。
        古い数字に現在の前処理条件のラベルが付いた TSV は、どちらも正しく見えて
        ずれが分からなくなります。

    InChIKey は DatasetState.feature_metadata から取る（mzTab-M の
    database_identifier / InChI / SMILES 由来。.arf2 との結合は不要）。
    濃縮解析の背景を保つため、有意な行だけでなく InChIKey が付いた全行を出す。
    ontology と msi_level は mzTab-M に対応物が無いため空欄で、その旨をメタ行に
    書く（空欄を「該当なし」と読み違えさせない）。メタ行には result_id と
    preprocess_id も入るので、あとから「どの計算の出力か」を辿れる。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    last = ds.last_differential
    if result_id is not None:
        last = (ds.results or {}).get(result_id)
    if not last or last.get("kind") != "two_group":
        return _missing(
            "dataset_differential_result",
            "先に dataset_preprocess → dataset_differential（2群）を実行してください。"
            if result_id is None else
            f"result_id={result_id} の 2 群比較の結果がありません。")

    from metabolomix.analysis.dataset_export import export_dataset_result
    from metabolomix.core.atomic_io import DomainError
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError, orphaned_warning

    curation = None
    flag_set = None
    arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p) for p in ds.source_files)), None)
    if arf2 is not None:
        try:
            flag_set = curation_apply.flags_for_arf2(arf2)
        except FlagFileError as exc:
            return mztab_error("CURATION_FLAGS_INVALID", str(exc), exc.details())
        if flag_set["n"]:
            # SMF_ID と MasterAlignmentID が同じ空間だと言えるのは、隣接 .arf との
            # 対応が検証済み（feature_axis = "arf_spot_index == mztab_smf_id"）のときだけ。
            mapped = (getattr(ds, "feature_qc", None) or {}).get("source") == "arf"
            state = ("not_applied" if not apply_curation
                     else "applied" if mapped else "unmapped")
            curation = {"state": state, "flag_set": flag_set}

    try:
        info = export_dataset_result(ds, last, Path(output_path), curation=curation)
    except DomainError as exc:
        if exc.code in ("STALE_ANALYSIS_RESULT", "ANALYSIS_RESULT_NOT_FOUND"):
            # 再計算すれば直る。どのツールを呼べばよいかを機械可読に返す。
            return _missing("dataset_differential_result", exc.message)
        return mztab_error(exc.code, exc.message, exc.details or None)

    curation_stats = info.pop("curation_stats", None)
    info["curation"] = curation_apply.payload_summary(info["curation"], flag_set, curation_stats)
    orphaned = (flag_set or {}).get("orphaned", 0)
    return json_payload({
        "status": "success",
        **info,
        **({"warnings": [orphaned_warning(orphaned)]} if orphaned else {}),
        "log2fc_sign": "log2fc は正なら group_b が高い（上昇）。",
        "note": ("n_unannotated は InChIKey が付かず書き出さなかった行数です。"
                 "「変化が無かった」ではなく「調べていない」行です。"),
    })


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def dataset_set_sample_metadata(manifest_path: str) -> str:
    """実験情報シート(sample-manifest.v1)を読み込み、検証してから DatasetState へ
    一括適用する（spec §7.3）。

    上流のConsole実行をやり直さずに、群・バッチ・注入順・qc_pool などの誤記入を
    このツール単体で訂正できる。列の意味・拒否条件は sample-manifest.v1 の仕様
    （`docs/superpowers/specs/2026-09-05-raw-folder-pipeline-integrity-and-metadata-design.md`
    §7.2）を参照。

    検証は全件そろってから一括反映する（**一部だけ適用される状態は作らない**）。
    シートに不備があれば `session.dataset` を一切変更せずエラーを返すので、
    そのまま再送してよい。

    前処理入力（role/batch/injection_order/qc_pool/include/sample_id/source_file）が
    変われば前処理済み行列ごと・PCA・差次的解析まで無効化する。`group` だけの変更なら
    差次的解析だけを無効化し、PCA は生かしたまま返す（`changed_fields` で確認できる）。

    `dataset_differential` を直接呼ぶ経路（比較する2群のサンプル名を都度自分で
    組み立てる）とは別に、比較を明示するのは pipeline 側の役割になる——群名だけで
    対照/処置の向きを決めない・QC/blank/unknown/include=false を混ぜない・完全交絡を
    止める、という前提検証は `resolve_comparison`/`run_comparison`（内部関数）が持つ。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    from metabolomix.analysis.dataset_service import apply_sample_manifest
    from metabolomix.core.atomic_io import DomainError

    try:
        summary = apply_sample_manifest(ds, manifest_path)
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)
    except OSError as exc:
        return mztab_error(
            "SAMPLE_MANIFEST_NOT_FOUND", str(exc), {"manifest_path": manifest_path})

    return json_payload({
        "status": "success",
        **summary,
        "next": ("changed_fields に前処理入力が含まれる場合は dataset_preprocess から"
                 "やり直してください。group のみの変更なら dataset_differential を"
                 "再実行するだけで構いません。"),
    })


# ---------- 内部ヘルパ ----------

def _missing(state: str, message: str) -> str:
    return missing_state(state, _RECOVERY_TOOLS[state], message)


def _from_precondition(exc: PreconditionError) -> str:
    """PreconditionError を kind に応じた封筒へ振り分ける。

    missing_state は「先に別のツールを呼べば直る」場合だけに使う。引数エラーを
    missing_state にすると、契約に従うクライアントが同じ呼び出しを繰り返す。
    """
    if exc.kind == "missing_state":
        return _missing(exc.state, exc.message)
    return mztab_error("DATASET_BAD_REQUEST", exc.message, exc.details or None)


def _count_roles(roles: dict, sample_names: list[str]) -> dict:
    counts: dict[str, int] = {}
    for name in sample_names:
        role = roles.get(name, "sample")
        counts[role] = counts.get(role, 0) + 1
    return counts


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_statistic(specification: dict, matrix_result_id: str) -> str:
    """v2 の統計（welch / anova_tukey / pca）を、名指しした解析行列に実行する。

    specification: `analysis_recipe.statistics` の1要素と同じ形
      （`statistic_id` / `kind` / `transform` / `feature_scope` と、kind別の
      `reference_group`+`test_group`、`groups`+`alpha`、`scaling`+`n_components`）。
    matrix_result_id: `analysis-matrix.v1` の ID。pipeline の preprocess /
      qc_processed が作った行列を名指しする。**直近の前処理を暗黙に使わない**
      ——recipe 違いの行列が2本ある前提の設計で、どちらの数字かを結果に残す。

    v1 の `dataset_differential` とは別物: 効果量は統計変換前の算術平均比の log2
    （`effect_size_definition`）、log2 変換に pseudocount を足さない、BH の母集団は
    検定できた feature だけ。v1 の既定値・数値契約は変更していない。

    全量（feature ごとの結果）は戻り値に載せず session に保持する。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")
    if not getattr(ds, "analysis_matrices", None):
        return _missing(
            "analysis_matrix",
            "解析行列（analysis-matrix.v1）がありません。"
            "pipeline_run で v2 要求を実行し、preprocess/qc_processed が作った"
            "matrix_result_id を指定してください。")

    from metabolomix.analysis.dataset_service import statistic_dataset
    from metabolomix.core.atomic_io import DomainError

    try:
        result = statistic_dataset(ds, specification, matrix_result_id)
    except DomainError as exc:
        return json_payload({"error": {"code": exc.code, "message": str(exc.message),
                                       "details": exc.details}})

    session_state.session.dataset.results[
        f"stat_{result.get('statistic_id')}"] = result
    payload = {k: v for k, v in result.items() if k != "features"}
    payload["status"] = "success"
    payload["n_features"] = len(result.get("features") or [])
    payload["features_note"] = (
        "feature ごとの結果全量は本要約に非同梱（セッションに保持）。"
        "TSV が必要なら pipeline の export 工程を使ってください。")
    return json_payload(payload)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_build_matrix(recipe: dict, recipe_id: str = "default") -> str:
    """v2 の解析行列（analysis-matrix.v1）を、セッションの DatasetState から1本作る。

    recipe: profile の `matrix_recipes` の1要素と同じ形
      （`base` / `normalize` / `drift_correct` / `filter` / `impute` の5キー。
      過不足はそのまま拒否する）。検証は profile と同じ規則を共有しているので、
      pipeline 経由より緩い recipe がここから入ることはない。
    recipe_id: 行列の識別に使う名前。`matrix_id` は recipe 内容と dataset・
      metadata から決まるので、**同じ dataset に同じ recipe を与えれば同じ ID**
      になる（pipeline が作った行列と突き合わせられる）。

    揃っていない前提を揃ったことにしない:

    - `base="internal_standard_ratio"` は対象feature↔内部標準の対応付けが要る。
      この単体ツールは対応付けを作らないので `missing_state` で止める
      （空の対応付けで「比」を名乗る行列を作らない）。
    - 検出状態が分からない dataset に `filter.min_detection_rate` を指定しても
      feature は落ちない。`history` に `not_evaluable` が残る。

    値そのものは戻り値に載せず session に保持する。`dataset_statistic` へは
    戻り値の `matrix_id` を渡す。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset",
                        "DatasetState がありません。先に dataset_load を実行してください。")

    from metabolomix.analysis.dataset_service import build_analysis_matrix
    from metabolomix.console.profile_schema import validate_matrix_recipe
    from metabolomix.core.atomic_io import DomainError

    try:
        normalized = validate_matrix_recipe(recipe, recipe_id)
    except DomainError as exc:
        # profile 経由なら PROFILE_INVALID だが、ここでの入力は profile ではなく
        # 引数そのもの。直す場所を取り違えさせない。
        return mztab_error("MATRIX_RECIPE_INVALID", exc.message, exc.details or None)

    if normalized["base"] == "internal_standard_ratio":
        return _missing(
            "feature_bindings",
            "内部標準比の行列には対象feature↔内部標準の対応付けが必要です。"
            "この単体ツールは対応付けを作らないため、pipeline の経路で作成してください。")

    try:
        matrix = build_analysis_matrix(ds, {**normalized, "recipe_id": recipe_id})
    except DomainError as exc:
        return json_payload({"error": {"code": exc.code, "message": str(exc.message),
                                       "details": exc.details}})

    return json_payload({
        "status": "success",
        "matrix_id": matrix["matrix_id"],
        "recipe_id": recipe_id,
        "stage": matrix["stage"],
        "n_assays": len(matrix["assay_ids"]),
        "n_features": len(matrix["feature_ids"]),
        "eligible_features": int(sum(bool(x) for x in matrix["eligibility_mask"])),
        "units": sorted(set(matrix["units"])),
        "history": matrix["correction_history"],
        "caveats": matrix["caveats"],
        "values_note": "値は戻り値に非同梱（セッションに保持）。dataset_statistic に matrix_id を渡す。",
    })


# ---------- dataset_feature_search（spec 2026-10-19） ----------

_POLARITIES = ("positive", "negative")

# 丸めの桁（spec §4.2）。m/z 4 桁・ppm 2 桁・RT 3 桁・log2FC 4 桁。p・q・平均は有効 4 桁。
_CANDIDATE_DIGITS = {"theoretical_mz": 4, "mz": 4, "ppm": 2, "rt": 3, "delta_rt": 3}
_WINDOW_DIGITS = {"mz": 4, "rt": 3}
_STATS_DIGITS = {"log2fc": 4}
_STATS_SIGNIFICANT = ("p", "q", "mean_a", "mean_b")


def _round(value, digits: int):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    value = float(value)
    if not math.isfinite(value):
        return value
    rounded = round(value, digits)
    return 0.0 if rounded == 0 else rounded          # "-0.0" を出さない


def _significant(value, digits: int = 4):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    value = float(value)
    if not math.isfinite(value):
        return value
    return float(f"{value:.{digits}g}")


def _cell(value):
    """TSV の 1 セル。数値は丸め済みのものを repr のまま出す（format_spots_as_table の
    固定小数 4 桁に任せると p = 1e-6 が 0 になる）。bool は JSON と同じ小文字。"""
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _table(rows: list[dict], columns, *, digits: dict | None = None,
           significant: tuple = ()) -> str:
    """行を TSV（列名 1 回）にする。行が無ければ空文字。"""
    from metabolomix.arf2.reader import format_spots_as_table

    digits = digits or {}
    cells = []
    for row in rows:
        out = {}
        for column in columns:
            value = row.get(column)
            if column in digits:
                value = _round(value, digits[column])
            elif column in significant:
                value = _significant(value)
            out[column] = _cell(value)
        cells.append(out)
    return format_spots_as_table(cells, columns=list(columns)) if cells else ""


def _dataset_polarity(ds) -> tuple[str | None, str, list[str]]:
    """データセットの極性・その出所（`job_declared` / `adduct_majority` / `unknown`）・注記。

    ジョブの宣言極性を優先し、mzTab-M の SML 付加体の多数決は裏取りに使う（spec §7）。
    両方あって食い違えば宣言を使い、食い違いを注記する。
    """
    notes: list[str] = []
    declared = None
    if getattr(ds, "job_path", None):
        from metabolomix.console import job_manager
        try:
            declared = job_manager.load_job(Path(ds.job_path)).polarity
        except (OSError, ValueError, KeyError, TypeError) as exc:
            notes.append(f"ジョブ記録を読めず、宣言極性を使っていない（{exc}）。")
        if declared is not None and declared not in _POLARITIES:
            notes.append(f"ジョブの宣言極性 {declared!r} を解釈できず、使っていない。")
            declared = None

    from metabolomix.console.output_collector import read_adduct_polarity
    majorities = {read_adduct_polarity(path)["adduct_majority"]
                  for path in (getattr(ds, "source_files", None) or {})}
    majorities.discard(None)
    majority = next(iter(majorities)) if len(majorities) == 1 else None
    if len(majorities) > 1:
        notes.append("読み込んだ mzTab-M の付加体の多数決が陽陰に割れたので、多数決は使っていない。")

    if declared is not None:
        if majority is not None and majority != declared:
            notes.append(
                f"ジョブの宣言極性（{declared}）と mzTab-M の付加体の多数決（{majority}）が"
                "食い違う。宣言を使った。")
        return declared, "job_declared", notes
    if majority is not None:
        return majority, "adduct_majority", notes
    notes.append(
        "データセットの極性を決められなかった（ジョブの宣言も mzTab-M の付加体の多数決も無い）。"
        "目標の polarity による適用判定（not_applicable）はしていない。")
    return None, "unknown", notes


def _curation_note(ds) -> str | None:
    """キュレーションのフラグがあれば「反映していない」注記。読めなくても検索は止めない。

    フラグの見つけ方は dataset_export_differential と同じ（兄弟の `.arf2`）。
    """
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation.flags import FlagFileError

    arf2 = next(filter(None, (curation_apply.arf2_for_mztab(p)
                              for p in (getattr(ds, "source_files", None) or {}))), None)
    if arf2 is None:
        return None
    try:
        n_flags = curation_apply.flags_for_arf2(arf2)["n"]
    except (FlagFileError, OSError) as exc:
        return (f"キュレーションのフラグ記録を読めなかった（{exc}）。この検索はキュレーションを"
                "反映していない（SME / SML の生の同定だけを見ている）。")
    if not n_flags:
        return None
    return (f"キュレーションのフラグが {n_flags} 件あるが、この検索は反映していない"
            "（SME / SML の生の同定だけを見ている）。dataset_export_differential の同定とは"
            "食い違いうる。")


def _target_summary(match: dict) -> dict:
    """目標 1 件の判定（候補の行は TSV の側に 1 回だけ載せるので外す）。"""
    keys = ("target_id", "status", "reason", "selected_feature_id", "selected_adduct",
            "selected_evidence", "adduct_tier_used", "counts", "truncated")
    out = {key: match[key] for key in keys}
    out["theoretical_mz"] = {adduct: _round(mz, 4)
                             for adduct, mz in match["theoretical_mz"].items()}
    return out


def _rule_targets(normalized: list[dict]) -> list[dict]:
    """戻り値に写す正規化後の目標（理論 m/z だけ 4 桁に丸める。hash は丸める前から作る）。"""
    shown = []
    for target in normalized:
        item = dict(target)
        item["theoretical_mz"] = {adduct: _round(mz, 4)
                                  for adduct, mz in target["theoretical_mz"].items()}
        shown.append(item)
    return shown


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def dataset_feature_search(
    feature_ids: list[str] | None = None,
    min_mz: float | None = None,
    max_mz: float | None = None,
    min_rt: float | None = None,
    max_rt: float | None = None,
    targets: list[dict] | None = None,
    result_ids: list[str] | None = None,
    max_results: int = 20,
) -> str:
    """DatasetState の特徴を探し、目的の物質に当たる特徴を規則で判定する（読み取り専用）。

    2 つの使い方があり、`targets` の有無で切り替わる（同時指定は `DATASET_BAD_REQUEST`）。

    **窓検索**（`targets` なし）: `feature_ids`（SMF_ID の完全一致）・`min_mz`/`max_mz`・
    `min_rt`/`max_rt`（分）の AND。境界を含む。何も指定しなければ `FEATURE_SEARCH_EMPTY_QUERY`
    （全件を返さない）。`feature_ids=[]` は「どの ID でもない」条件で 0 件。並びは m/z 昇順。
    戻り値 `features` は TSV（`feature_id mz rt evidence name inchikey adduct n_sme`。
    `evidence` は `sme`/`sml`/`none`）。注釈の無い特徴にも統計を付ける。

    **目標照合**（`targets` あり、1 回 50 件まで）: 各要素は `target_id`・`inchikey`（必須。
    先頭 14 文字で照合）と、`formula` か（`expected_mz` と `adduct`）。任意で `name` `ppm`
    （既定 10）`rt_min` `rt_tol_min`（既定 0.3）`adduct_tiers` `polarity`。付加体の層を上から、
    証拠の水準 `sme_match` → `sml_match` → `unannotated` の順に見て、qualified がちょうど 1 つの
    ときだけ `resolved`。2 つ以上なら `needs_input`（`multiple_qualified_candidates`）、無ければ
    `needs_input`（`no_qualified_candidate`）。**強度・ppm・RT の近さでは選ばない。**
    `unannotated` は `rt_min` があるときだけ qualified。SME（MS/MS 同定）が別物質なら落とす。
    極性はジョブの宣言 → mzTab-M の付加体の多数決の順に決め、`rule.polarity_source` に書く。
    戻り値 `targets`（目標ごとの判定と段ごとの件数）・`candidates`（TSV。目標ごとに
    `max_results` 件まで）・`rule`（正規化した規則と `rule_hash`。result_id を含まないので、
    同じ目標・同じ極性なら別の session でも同じ hash）。hash は表示用の丸め（理論 m/z を
    4 桁に）より前の値から作るので、返された `rule.targets` から計算し直すことはできない
    ——照合は `rule_hash` の値どうしで行う。

    result_ids: 統計を付ける 2 群比較（`dataset_differential`）の結果。省略時は今も有効な
        2 群比較の全部（古いものは使わず `notes` に挙げる）。名指しした ID が無い・古いときは
        `missing_state`。差次的結果が 1 件も無ければ統計なしで返す。
    max_results: 窓検索の行数・目標ごとの候補数の上限（既定 20）。打ち切りは `notes` に書く。

    `stats` は TSV（`feature_id result_id status log2fc p q mean_a mean_b detected_a detected_b`）。
    `status` は `tested` / `not_in_matrix`（前処理で落ちた）/ `not_testable`（p が無い）で、
    後の 2 つは数値欄が空（「差が無かった」ではない）。`detected_*` は群ごとの実検出例数
    （gap-fill を数えない。検出状態が無ければ空）。`results` に BH の母集団
    （`bh_population`）と caveat。キュレーションのフラグは反映しない（あれば `notes`）。
    欄の意味は `lipidmix://docs/output-format/mztab`。
    """
    ds = session_state.session.dataset
    if ds is None:
        return _missing("dataset", "DatasetState がありません。先に dataset_load を実行してください。")

    window = {"feature_ids": feature_ids, "min_mz": min_mz, "max_mz": max_mz,
              "min_rt": min_rt, "max_rt": max_rt}
    given_window = [name for name, value in window.items() if value is not None]
    if targets is not None and given_window:
        return mztab_error(
            "DATASET_BAD_REQUEST",
            "targets（目標照合）と窓の引数（feature_ids・m/z・RT）は同時に指定できません"
            "（どちらの意味か曖昧になる）。どちらか一方で呼び直してください。",
            {"window_arguments": given_window})
    if isinstance(max_results, bool) or not isinstance(max_results, int) or max_results < 1:
        return mztab_error("DATASET_BAD_REQUEST",
                           f"max_results は 1 以上の整数で指定してください: {max_results!r}",
                           {"field": "max_results"})

    from metabolomix.analysis import feature_search
    from metabolomix.analysis.result_state import dataset_fingerprint
    from metabolomix.core.atomic_io import DomainError

    notes: list[str] = []
    payload: dict = {}
    try:
        if targets is None:
            rows, truncated = feature_search.search_window(
                ds, **window, max_results=max_results)
            if feature_ids is not None and not feature_ids:
                notes.append("feature_ids=[] は「どの ID でもない」条件として扱った（0 件）。")
            if truncated:
                notes.append(f"max_results={max_results} で打ち切った（ほか {truncated} 件）。"
                             "窓を狭めるか max_results を増やしてください。")
            payload["features"] = _table(rows, feature_search.WINDOW_COLUMNS,
                                         digits=_WINDOW_DIGITS)
            stat_features = [row["feature_id"] for row in rows]
        else:
            if isinstance(targets, list) and not targets:
                return mztab_error("FEATURE_TARGET_INVALID", "targets が空です（1 件以上渡す）。",
                                   {"target_id": None, "field": "targets"})
            polarity, polarity_source, polarity_notes = _dataset_polarity(ds)
            notes.extend(polarity_notes)
            normalized = feature_search.normalize_targets(targets, dataset_polarity=polarity)
            # 特徴の m/z 配列は 1 回だけ作り、全目標で使い回す
            mz_array = feature_search._feature_mz_array(ds)
            matches = [feature_search.match_target(ds, target, max_results=max_results,
                                                   mz_array=mz_array)
                       for target in normalized]
            cut = [f"{m['target_id']}: {m['truncated']} 件" for m in matches if m["truncated"]]
            if cut:
                notes.append(f"候補を目標ごとに max_results={max_results} 件で打ち切った"
                             f"（{', '.join(cut)}）。")
            candidates = [row for m in matches for row in m["candidates"]]
            payload["targets"] = [_target_summary(m) for m in matches]
            payload["candidates"] = _table(candidates, feature_search.CANDIDATE_COLUMNS,
                                           digits=_CANDIDATE_DIGITS)
            stat_features = list(dict.fromkeys(row["feature_id"] for row in candidates))
    except DomainError as exc:
        return mztab_error(exc.code, exc.message, exc.details or None)

    if result_ids is not None:
        chosen = feature_search.select_results(ds, list(dict.fromkeys(result_ids)))
        if chosen["missing"] or chosen["stale"]:
            # 名指しの結果を黙って落とすと、呼び出し側は統計が付いたつもりで読む
            return _missing(
                "dataset_differential_result",
                "result_ids に指定した差次的結果が使えません"
                f"（無い・2 群比較でない: {chosen['missing']} / 古い: {chosen['stale']}）。"
                "dataset_differential で計算し直すか、result_ids を省略して今も有効な結果だけを"
                "使ってください。")
    else:
        chosen = feature_search.select_results(ds, None)
        if chosen["stale"]:
            notes.append("古い差次的結果（前処理・群の訂正前）は統計に使っていない: "
                         + ", ".join(chosen["stale"]))
    results = chosen["results"]
    if not results:
        notes.append("差次的結果（2 群比較）が 1 件も無いので統計を付けていない"
                     "（dataset_differential を実行すると付く）。")

    stats = feature_search.feature_stats(ds, stat_features, results)
    if any(row["detected_a"] is None for row in stats):
        if getattr(ds, "detected_mask", None) is None:
            notes.append("検出状態（gap-fill の区別）を取り込めていないので detected_a / "
                         "detected_b は空（非ゼロ値は実測と補間を区別しない）。")
        else:
            notes.append("結果の試料と検出状態の列の対応を保証できず、detected_a / detected_b を"
                         "数えていない結果がある。")
    payload["stats"] = _table(stats, feature_search.STATS_COLUMNS, digits=_STATS_DIGITS,
                              significant=_STATS_SIGNIFICANT)
    payload["results"] = feature_search.result_summaries(ds, results)
    if targets is not None:
        payload["rule"] = {"targets": _rule_targets(normalized), "polarity": polarity,
                           "polarity_source": polarity_source,
                           "rule_hash": feature_search.rule_hash(normalized)}
    payload["stats_source"] = {
        "result_ids": [(r.get("provenance") or {}).get("result_id") for r in results],
        "dataset_hash": dataset_fingerprint(ds),
    }
    curation = _curation_note(ds)
    if curation:
        notes.append(curation)
    payload["notes"] = notes
    return json_payload(payload)
