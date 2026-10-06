"""DatasetState を明示的に受け取る解析サービス（spec §6）。

MCP ツールも pipeline ワーカーもここを通る。ツール層はグローバルな
`session.dataset` を渡すだけ、ワーカーは自分が持つ DatasetState を渡すだけで、
**同じ関数が同じ更新規則を実行する**。実行方式によって状態の更新順が変わると、
「単体ツールでは無効化されるが pipeline では残る」ような差が生まれる。

この層の約束は 3 つ:

1. **成功してからまとめて反映する。** 前処理は途中で失敗しうる。行列だけ差し替えて
   レシピは古いまま、という半端な状態を作らない（`preprocess_dataset`）。
2. **入力が変わったら派生結果を捨てる。** 判断は `result_state.invalidate_results`
   に一本化する。
3. **結果に来歴を付ける。** 既存の戻り値のキーは一切改名せず、`provenance` を
   足すだけ（`dataset_export_differential` や図の保存が読むキーを動かさない）。
"""
from __future__ import annotations

import os
from pathlib import Path

from metabolomix.analysis import differential, result_state
from metabolomix.analysis.dataset_analysis import (
    run_dataset_differential,
    run_dataset_pca,
    run_dataset_preprocess,
)
from metabolomix.analysis.preprocess_policy import check_applied_policy, resolve_policy
from metabolomix.analysis.sample_manifest import (
    apply_metadata,
    parse_manifest,
    resolve_metadata,
    select_statistical_samples,
)
from metabolomix.core.atomic_io import DomainError

__all__ = [
    "apply_sample_manifest",
    "build_analysis_matrix",
    "compare_dataset",
    "pca_dataset",
    "preprocess_auto",
    "preprocess_dataset",
    "resolve_comparison",
    "run_comparison",
    "statistic_dataset",
]


def _reusable_preprocess(ds, fingerprint: str) -> dict | None:
    """同じ入力・同じレシピの前処理結果が既にあるならそれを返す（無ければ None）。

    同値のレシピを同じ入力へ再適用した場合は**計算し直さない**。同じ計算に別の
    ID を付けると来歴が「別の実行」に見えるし、無関係な再計算を強いる。
    """
    current = ds.results.get(ds.preprocess_id) if ds.preprocess_id else None
    if (current is not None and ds.pp_matrix is not None
            and current["provenance"]["input_fingerprint"] == fingerprint):
        return current
    return None


def _commit_preprocess(ds, recipe: dict, computed: tuple, *, metadata_hash: str | None,
                       fingerprint: str, request_revision: int | None) -> dict:
    """計算済みの前処理結果を ds へ反映して登録する（whole-or-nothing の後半）。

    `computed` は `run_dataset_preprocess` の戻り値そのもの。計算と反映を
    分けてあるのは、呼び出し側（`preprocess_auto`）が反映の**前に**結果を
    検査できるようにするため——検査で落ちれば ds は一切書き換わらない。
    """
    (pp_matrix, pp_sample_names, pp_feature_names, roles, sample_meta, report) = computed

    # 前処理をやり直した以上、前の行列から出た結果は全部古い。
    result_state.invalidate_results(ds, {"recipe"})

    ds.pp_matrix = pp_matrix
    ds.pp_sample_names = pp_sample_names
    ds.pp_feature_names = pp_feature_names
    ds.roles = roles
    ds.sample_meta = sample_meta
    ds.preprocessing_recipe = dict(recipe)
    ds.preprocess_metadata_hash = metadata_hash

    result = dict(report)
    result["provenance"] = result_state.new_provenance(
        ds, kind="preprocess", input_fingerprint=fingerprint,
        effective_parameters=dict(recipe),
        warnings=list(report.get("caveats", [])),
        request_revision=request_revision)
    ds.preprocess_id = result["provenance"]["result_id"]
    return result_state.register_result(ds, result)


def preprocess_dataset(ds, recipe: dict, *, request_revision: int | None = None) -> dict:
    """前処理を適用し、成功したときだけ状態を差し替える。

    同値のレシピを同じ入力へ再適用した場合は**計算し直さず**、前回の結果と
    result_id をそのまま返す（派生結果も生かしたまま）。

    Returns
    -------
    dict
        `run_dataset_preprocess` の report に `provenance` を足したもの。
    """
    metadata_hash = _metadata_hash(ds)
    fingerprint = result_state.preprocess_fingerprint(ds, recipe, metadata_hash)

    reusable = _reusable_preprocess(ds, fingerprint)
    if reusable is not None:
        return reusable

    # ここから先が「新しい前処理」。計算は一時変数へ受け、全部成功してから反映する。
    computed = run_dataset_preprocess(ds, recipe)
    return _commit_preprocess(ds, recipe, computed, metadata_hash=metadata_hash,
                              fingerprint=fingerprint, request_revision=request_revision)


def preprocess_auto(ds, requested: dict, metadata: list[dict], *,
                    request_revision: int | None = None) -> dict:
    """conservative-v1でレシピを解決し、検査してから前処理をcommitする（spec §8）。

    順序は「一時計算 → 検査 → commit」を守る:

    1. `resolve_policy` でレシピを決める（``ds`` は書き換えない）。
    2. 実験情報シートを確定する（`apply_metadata`。Task 9 の一括適用契約——
       検証で1件でも落ちれば ``ds`` には触れない）。
    3. `run_dataset_preprocess` で**一時的に**計算する（``ds.pp_matrix`` 等はまだ
       書き換わらない）。
    4. `check_applied_policy` でその結果を検査する。ここで例外なら ``ds`` の
       前処理系フィールドは一切書き換わっていない（前段の metadata 確定だけが
       残る——実験情報の確定自体は前処理の成否と独立に有効な情報のため）。
    5. 検査を通ってはじめて `_commit_preprocess` で本commitする。

    **計算は1回だけ**。`preprocess_dataset` をそのまま呼ぶと、検査用に一度
    計算した同じ行列を commit のためにもう一度計算することになる（間に ``ds``
    の入力は何も変わらないので、二度目は結果まで同じ）——検体数×特徴量数の
    行列に対する全量の重複計算で、得られるものは無い。ここでは
    `preprocess_dataset` と同じ再利用判定・同じ commit 手順を、計算を1回に
    畳んだ順序で組み立てる。
    """
    plan = resolve_policy(ds, requested, metadata)
    apply_metadata(ds, metadata)
    recipe = plan["resolved_recipe"]

    metadata_hash = _metadata_hash(ds)
    fingerprint = result_state.preprocess_fingerprint(ds, recipe, metadata_hash)
    reusable = _reusable_preprocess(ds, fingerprint)
    if reusable is not None:
        # 同じ入力・同じレシピの結果がもう手元にある。検査だけはやり直す
        # ——policy の解決が前回と同じ結論に達したことを確かめる材料は
        # レポート（＝この結果）に全部入っている。
        check_applied_policy(plan, reusable)
        result = reusable
    else:
        computed = run_dataset_preprocess(ds, recipe)
        check_applied_policy(plan, computed[-1])  # 6要素目が report
        result = _commit_preprocess(ds, recipe, computed, metadata_hash=metadata_hash,
                                    fingerprint=fingerprint,
                                    request_revision=request_revision)

    plan["applied_steps"] = list(result.get("recipe_applied", []))
    plan["result"] = result
    return plan


def pca_dataset(ds, n_components: int = 5, log_transform: bool = False, *,
                request_revision: int | None = None) -> dict:
    """PCA を実行し、成功したときだけ `last_pca` を差し替える。"""
    settings = {"n_components": n_components, "log_transform": log_transform}
    result = run_dataset_pca(ds, n_components=n_components,
                             log_transform=log_transform)
    result["provenance"] = result_state.new_provenance(
        ds, kind="pca",
        input_fingerprint=_derived_fingerprint(ds, settings),
        effective_parameters=settings,
        parent_ids=[ds.preprocess_id] if ds.preprocess_id else [],
        warnings=list(result.get("caveats", [])),
        request_revision=request_revision)
    ds.last_pca = result
    return result_state.register_result(ds, result)


def compare_dataset(ds, group_a: list[str], group_b: list[str], *,
                    q_threshold: float = 0.05, log2fc_threshold: float = 1.0,
                    log_transform: bool = True,
                    group_a_label: str = "group_a",
                    group_b_label: str = "group_b",
                    request_revision: int | None = None) -> dict:
    """2 群比較を実行し、成功したときだけ `last_differential` を差し替える。

    探索専用のデータセット（中断された実行の出力）では拒否する。欠けた検体は
    mzTab 上「その群に無い」ようにしか見えず、比較はその欠落ごと結論にする。
    """
    _refuse_exploratory(ds, "2 群比較")
    settings = {
        "group_a": list(group_a), "group_b": list(group_b),
        "q_threshold": q_threshold, "log2fc_threshold": log2fc_threshold,
        "log_transform": log_transform,
        "group_a_label": group_a_label, "group_b_label": group_b_label,
    }
    result = run_dataset_differential(
        ds, group_a_samples=group_a, group_b_samples=group_b,
        q_threshold=q_threshold, log2fc_threshold=log2fc_threshold,
        log_transform=log_transform,
        group_a_label=group_a_label, group_b_label=group_b_label)
    result["provenance"] = result_state.new_provenance(
        ds, kind="differential",
        input_fingerprint=_derived_fingerprint(ds, settings),
        effective_parameters=settings,
        parent_ids=[ds.preprocess_id] if ds.preprocess_id else [],
        warnings=list(result.get("caveats", [])),
        request_revision=request_revision)
    # 群割当の指紋を来歴へ刻む。差次的結果だけは前処理・データセットが同じでも
    # 「どの群割当で計算したか」で古くなる——群を訂正したあとにこの結果を
    # result_id で名指しされたとき、`result_state.is_current` がここを読んで
    # 拒否する（訂正前の数字が現在の条件のラベル付きで書き出されるのを防ぐ）。
    result["provenance"]["group_fingerprint"] = result_state.group_fingerprint(ds)
    ds.last_differential = result
    return result_state.register_result(ds, result)


def resolve_comparison(ds, comparison: dict, metadata: list[dict]) -> dict:
    """比較定義とサンプルメタデータから、実際に比較へ回す2群を解決する（spec §7.4）。

    群名から向きを推測しない——``comparison`` に ``reference_group``/``test_group``が
    無ければ、群ラベルがどれだけ「対照/処置」らしく見えても ``COMPARISON_REQUIRED`` にする。

    ``metadata`` は ``ds.sample_names`` と同じ順・同じ件数の sample-manifest.v1 行
    （Task 9 ``resolve_metadata``/``metadata_rows`` と同じ契約）。``include=true`` かつ
    ``role=="sample"`` の行だけを比較対象にする——QC・blank・unknown（Task 9 が
    ``build_dataset_pp_inputs`` で素通りさせた3値目）は前処理入力には残るが、
    比較には混ぜない（`docs/superpowers/.../task-12-brief.md` step 3）。

    行の ``sample_id`` は前処理へ渡す実際のサンプル名（``ds.sample_names``）とは別の
    ユーザー定義IDになりうるため、位置対応（``zip(ds.sample_names, metadata)``）で
    変換してから返す——差次的解析（`compare_dataset`）は ``ds.pp_sample_names`` に
    含まれる名前しか受け付けない。

    Returns
    -------
    dict
        ``reference_group``/``test_group``（実際に使った群ラベル）、
        ``reference_sample_ids``/``test_sample_ids``（sample-manifest.v1 の sample_id）、
        ``reference_samples``/``test_samples``（`compare_dataset` にそのまま渡せる
        サンプル名）、``confounding``（``differential.check_confounding`` の結果）、
        ``allow_confounded``（明示されたoverride）を持つ。
    """
    reference_group = comparison.get("reference_group")
    test_group = comparison.get("test_group")
    if not reference_group or not test_group or reference_group == test_group:
        raise DomainError(
            "COMPARISON_REQUIRED",
            "対照群(reference_group)と比較群(test_group)を異なる値で明示してください"
            "（群名だけでは比較の向きを決めません）。",
            {"reference_group": reference_group, "test_group": test_group},
        )

    sample_names = list(getattr(ds, "sample_names", []) or [])
    if len(metadata) != len(sample_names):
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"metadataの行数がサンプル数と一致しません(rows={len(metadata)}, "
            f"samples={len(sample_names)})。",
            {"rows": len(metadata), "samples": len(sample_names)},
        )

    id_to_name: dict[str, str] = {}
    for name, row in zip(sample_names, metadata):
        sample_id = row.get("sample_id")
        if not sample_id:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"sample_idが空です（sample={name!r}）。",
                {"sample_name": name},
            )
        if sample_id in id_to_name:
            raise DomainError(
                "SAMPLE_MANIFEST_INVALID",
                f"sample_idが重複しています: {sample_id!r}",
                {"sample_id": sample_id},
            )
        id_to_name[sample_id] = name

    # 比較対象の選択規則（include=true・role=="sample"・指定群）は
    # `sample_manifest.select_statistical_samples` が唯一の定義。ここで書き直すと、
    # v2で足した role="standard"（標準品注入）が片方の経路にだけ混ざる。
    selected = select_statistical_samples(metadata, [reference_group, test_group])

    def _group_members(label: str) -> tuple[list[str], list[str], list]:
        members = [row for row in selected if row.get("group") == label]
        ids = [row["sample_id"] for row in members]
        names = [id_to_name[sid] for sid in ids]
        batches = [row.get("batch") for row in members]
        return ids, names, batches

    ref_ids, ref_names, ref_batches = _group_members(reference_group)
    test_ids, test_names, test_batches = _group_members(test_group)

    if len(ref_ids) < 2:
        raise DomainError(
            "COMPARISON_GROUP_TOO_SMALL",
            f"reference_group={reference_group!r} に該当する試料が2件未満です"
            f"(n={len(ref_ids)})。群ラベルの誤り、またはinclude/role除外の可能性があります。",
            {"group": reference_group, "n": len(ref_ids)},
        )
    if len(test_ids) < 2:
        raise DomainError(
            "COMPARISON_GROUP_TOO_SMALL",
            f"test_group={test_group!r} に該当する試料が2件未満です"
            f"(n={len(test_ids)})。群ラベルの誤り、またはinclude/role除外の可能性があります。",
            {"group": test_group, "n": len(test_ids)},
        )

    group_labels = [reference_group] * len(ref_ids) + [test_group] * len(test_ids)
    batch_labels = ref_batches + test_batches
    confounding = differential.check_confounding(group_labels, batch_labels)

    return {
        "comparison_id": comparison.get("comparison_id"),
        "reference_group": reference_group,
        "test_group": test_group,
        "reference_sample_ids": ref_ids,
        "test_sample_ids": test_ids,
        "reference_samples": ref_names,
        "test_samples": test_names,
        "confounding": confounding,
        "allow_confounded": bool(comparison.get("allow_confounded", False)),
    }


def run_comparison(ds, comparison: dict, metadata: list[dict]) -> dict:
    """比較前提を検証してから `compare_dataset` を実行する（spec §7.4）。

    完全交絡（群⟂バッチが分離不能）は ``allow_confounded=true`` の明示が無い限り
    ``CONFOUNDED_COMPARISON`` で止める。バッチ情報不足で交絡を評価できない場合
    （``confounding.assessable=False``）は「評価不可」であって「交絡あり」ではない
    ——ここで止めない（spec §7.4、common-context 参照）。

    正の log2FC は ``test_group`` が高い方向に固定する（`compare_dataset` の
    ``group_b`` に ``test_group`` を渡すことで実現。`export_contract.LOG2FC_SIGN`
    と整合）。
    """
    resolved = resolve_comparison(ds, comparison, metadata)
    confounding = resolved["confounding"]
    if confounding["confounded"] and not resolved["allow_confounded"]:
        raise DomainError(
            "CONFOUNDED_COMPARISON",
            f"{resolved['reference_group']!r} と {resolved['test_group']!r} は群と"
            "バッチが完全に交絡しており、処理効果と測定バッチを分離できません。"
            "allow_confounded=true を明示しない限り差次的解析は実行しません。",
            {"reference_group": resolved["reference_group"],
             "test_group": resolved["test_group"],
             "detail": confounding["detail"]},
        )

    result = compare_dataset(
        ds, resolved["reference_samples"], resolved["test_samples"],
        q_threshold=comparison.get("q_threshold", 0.05),
        log2fc_threshold=comparison.get("log2fc_threshold", 1.0),
        log_transform=comparison.get("log_transform", True),
        group_a_label=resolved["reference_group"],
        group_b_label=resolved["test_group"],
    )

    unadjusted = bool(confounding["confounded"] and resolved["allow_confounded"])
    result["provenance"]["comparison"] = {
        "comparison_id": resolved["comparison_id"],
        "reference_group": resolved["reference_group"],
        "test_group": resolved["test_group"],
        "reference_sample_ids": resolved["reference_sample_ids"],
        "test_sample_ids": resolved["test_sample_ids"],
        "confounding": confounding,
        "allow_confounded": resolved["allow_confounded"],
        "unadjusted_confounded": unadjusted,
    }
    if unadjusted:
        result["caveats"].append(
            "allow_confounded=true が明示されたため、群とバッチの完全交絡を未調整の"
            "まま解析を継続しました（図・TSV・レポートにこの旨を残すこと）。")
    # compare_dataset は provenance.warnings を「compare_dataset を呼んだ時点の
    # result['caveats']」でスナップショットする（result_state.new_provenance に
    # list(...) で渡すコピー）。run_comparison がこの後に caveats へ追記した注記
    # （上記の未調整注記を含む）はそのコピーには反映されない——spec §7.4が求める
    # 「未調整であることを...レポートへ記録して継続」の記録先は永続化された
    # provenance（図/TSV/レポートが読む来歴）であって result['caveats'] という
    # その場限りの戻り値ではないため、ここで両者を同期し直す（レビュー Finding 2）。
    result["provenance"]["warnings"] = list(result["caveats"])
    return result


def apply_sample_manifest(ds, manifest_path: str) -> dict:
    """実験情報シート(sample-manifest.v1)を読み込み、検証してから一括適用する（spec §7.3）。

    Task 9 の ``parse_manifest`` → ``resolve_metadata`` → ``apply_metadata`` を
    この順に呼ぶだけの薄い合成——検証で1件でも落ちれば ``apply_metadata`` 自身の
    契約により ``ds`` には一切触れない（対話的なDatasetState経路にも同じ検証・出所・
    無効化規則を適用する、spec §7.3）。

    ``DatasetState`` は ``source_root`` 自体を保持しないため、実行時に予定した raw の
    絶対パス（``ds.assay_sources``）から逆算する。
    """
    source_root, expected_sources = _raw_manifest_layout(ds)
    rows = parse_manifest(Path(manifest_path), source_root=source_root,
                          expected_sources=expected_sources)
    resolved = resolve_metadata(ds, rows)
    return apply_metadata(ds, resolved)


def _raw_manifest_layout(ds) -> tuple[Path, list[str]]:
    """``ds.assay_sources``（実行時に予定したrawの絶対パス）から、
    sample-manifest.v1 検証に使う ``source_root``/相対パス一覧を逆算する。
    """
    paths = [Path(p) for p in (getattr(ds, "assay_sources", None) or {}).values() if p]
    if not paths:
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            "raw参照(assay_sources)が無いため実験情報シートを検証できません"
            "（console_run経由のjob_pathからdataset_loadしたデータセットが必要です）。",
            {},
        )
    root = paths[0].parent if len(paths) == 1 else Path(os.path.commonpath([str(p) for p in paths]))
    expected = [str(p.relative_to(root)) for p in paths]
    return root, expected


def _refuse_exploratory(ds, what: str) -> None:
    """探索専用のデータセットで結論を出す操作を止める。"""
    from metabolomix.core.atomic_io import DomainError

    if getattr(ds, "exploratory_only", False):
        raise DomainError(
            "EXPLORATORY_ONLY_DATASET",
            f"このデータセットは中断された解析の出力（探索専用）なので、{what}には"
            "使えません。Console 実行を完了させてから読み込み直してください。",
            {"job_path": getattr(ds, "job_path", None),
             "source_verification": getattr(ds, "source_verification", None)})


# ---------- 内部ヘルパ ----------

def _metadata_hash(ds) -> str | None:
    """前処理に効くサンプルメタデータの指紋（明示メタデータが無ければ None）。

    明示メタデータ（実験情報シート由来）は Task 9 で `ds.sample_metadata_rows` へ
    入る。無いうちは None で、指紋はデータセットとレシピだけで決まる。
    """
    rows = getattr(ds, "sample_metadata_rows", None)
    if not rows:
        return None
    from metabolomix.core.atomic_io import canonical_hash

    fingerprints = result_state.metadata_fingerprints(rows)
    # 前処理に効く列だけを混ぜる。group を混ぜると、群を付け替えただけで
    # 前処理の指紋が変わり、無関係な再前処理を強いる。
    return canonical_hash({k: v for k, v in fingerprints.items()
                           if k in result_state.PP_FIELDS})


def _derived_fingerprint(ds, settings: dict) -> str:
    """前処理済み行列と設定から、派生結果の入力指紋を作る。"""
    from metabolomix.core.atomic_io import canonical_hash

    return canonical_hash({
        "preprocess_id": ds.preprocess_id,
        "pp_matrix": result_state.array_fingerprint(ds.pp_matrix),
        "samples": list(ds.pp_sample_names),
        "features": list(ds.pp_feature_names),
        "settings": settings,
    })


# ---------- v2: 解析行列と統計（spec §8.2・§10） ----------

def build_analysis_matrix(ds, recipe: dict, bindings: dict | None = None,
                          evidence: dict | None = None,
                          eligibility=None) -> dict:
    """`analysis-matrix.v1` を1本作り、`ds.analysis_matrices` へ登録する。

    `ds.pp_matrix`（v1 の単一スロット）へは**書かない**。v2 は recipe ごとに
    行列を持つので、単一スロットに置くと2本目が1本目を黙って上書きする
    ——統計が「どの行列の数字か」を言えなくなる。
    """
    import numpy as np

    from metabolomix.analysis import matrix_state

    if eligibility is None:
        eligibility = np.ones(len(getattr(ds, "feature_ids", []) or []), dtype=bool)
    matrix = matrix_state.make_matrix(ds, recipe, bindings or {}, evidence or {},
                                      eligibility)
    ds.analysis_matrices[matrix["matrix_id"]] = matrix
    return matrix


def statistic_dataset(ds, specification: dict, matrix_result_id: str) -> dict:
    """指定した解析行列に対して v2 統計を1件実行する（spec §10）。

    行列は `matrix_result_id` で**名指し**する。「直近の前処理」を暗黙に使うと、
    recipe 違いの行列が2本あるときにどちらの数字か言えなくなる。該当IDが
    無ければ `ANALYSIS_RESULT_NOT_FOUND` で止める（近い行列で代用しない）。

    数値の実装はここには無い——`analysis/statistics_v2.py` が唯一の実装で、
    この層は「どの行列・どのメタデータで呼ぶか」を解決するだけ。
    """
    from metabolomix.analysis import statistics_v2

    matrices = getattr(ds, "analysis_matrices", None) or {}
    matrix = matrices.get(matrix_result_id)
    if matrix is None:
        raise DomainError(
            "ANALYSIS_RESULT_NOT_FOUND",
            f"指定されたmatrix_result_idの解析行列がありません: {matrix_result_id!r}"
            "（pipelineのpreprocess/qc_processedが作る行列IDを指定してください。"
            "近い行列で代用はしません）。",
            {"matrix_result_id": matrix_result_id,
             "available_matrix_ids": sorted(matrices)})

    metadata = getattr(ds, "sample_metadata_rows", None) or []
    target_features = getattr(ds, "feature_binding_targets", None) or {}
    return statistics_v2.run_statistic(matrix, specification, metadata,
                                       target_features)
