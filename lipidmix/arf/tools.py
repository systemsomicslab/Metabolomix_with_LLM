"""ARF（PeakProperties 等・サンプル別強度）ツール群。

arf_list_tags/classes/sample_roles, arf_exclude, arf_preprocess,
arf_pca_preprocessed, arf_parser, arf_differential, arf_export_differential。
deps: mcp_core / session_state /
path_resolvers / tool_helpers / msdial_* / preprocessing / differential /
arf_reader（arf_reader は関数内 import で、テストの patch.object(server.arf_reader)
が共有 module 経由で効くようにする）。

_filter_arf_spots / _pp_build_matrix はテストが差し替える対象なので、正準定義元を
module 修飾（path_resolvers.* / tool_helpers.*）で参照し patch が確実に効くようにする。
"""
import math
import re
from pathlib import Path

from lipidmix.analysis import differential
from lipidmix.analysis import export_contract
from lipidmix.arf import exclusions
from lipidmix.arf import identity_join
from lipidmix.core import mcp_errors
from lipidmix.core import path_resolvers
from lipidmix.analysis import preprocessing
from lipidmix.msdial import lipid_identity, sample_factors
from lipidmix.core import session_state
from lipidmix.core import tool_helpers
from lipidmix.plots import render as plot_render
from lipidmix.plots import volcano as volcano_plot
from mcp.server.fastmcp import Image
from mcp.types import ToolAnnotations
from pydantic import StrictInt
from lipidmix.core.mcp_core import mcp
from lipidmix.core.serialization import json_payload
from lipidmix.arf2.reader import format_spots_as_table
from lipidmix.msdial.classes import assign_sample_groups, filter_arf_by_class_ids
from lipidmix.msdial.tags import filter_arf_by_tags
from lipidmix.core.path_resolvers import resolve_arf_file_path
from lipidmix.core.session_state import _build_sample_meta
from lipidmix.core.tool_helpers import (
    _pp_has_preprocessed,
    _class_factors_by_position,
    _remember_arf_pca_plot,
    _format_pca_plot_block,
    _format_pca_loadings_md,
    _format_arf_tag_summary,
    _format_arf_class_summary,
    _format_arf_parse_summary,
    _format_arf_class_filter,
    _format_arf_tag_filter,
)

__all__ = [
    "arf_list_tags",
    "arf_list_classes",
    "arf_list_sample_roles",
    "arf_exclude",
    "arf_preprocess",
    "arf_pca_preprocessed",
    "arf_parser",
    "arf_differential",
    "arf_plot_volcano",
    "arf_export_differential",
]


# readOnlyHint の意味: 「サーバの外に副作用が無い」＝ファイルとネットワークを変更
# しない。ARF ツールは session_state.session.arf を更新するが、これはサーバ自身の
# 解析セッション状態であり副作用に数えない（その依存は missing_state エンベロープで
# 伝えるため、annotations で二重に表現しない）。
@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_list_tags() -> str:
    """List MS-DIAL tags discovered for the currently loaded ARF dataset."""
    if (
        session_state.session.arf.features is None
        or session_state.session.arf.tag_index is None
        or not str(session_state.session.arf.current_file_path or "").lower().endswith(".arf")
    ):
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"],
            "先に arf_parser を実行してARFデータとタグファイルを読み込んでください。")
    return json_payload(session_state.session.arf.tag_index.get("summary", {}))


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_list_classes() -> str:
    """絞り込み・群分けに使える因子トークンの一覧を返す（Class ID とサンプル名の両方）。

    class_ids / group_levels / group_factors / group_a / group_b に何を書けるかを
    発見するための入口。`.mddata` が無くてもサンプル名の語彙は返る。
    `sample_token_vocabulary.tokens[<token>]` は出現サンプル数・role 内訳・
    サンプル名内の出現位置を持つ。位置は因子の並びを推測する手掛かりだが、値が
    2トークンに割れる場合（G_uralensis）にずれるため指定には使わない。
    """
    if (
        session_state.session.arf.features is None
        or not str(session_state.session.arf.current_file_path or "").lower().endswith(".arf")
    ):
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"],
            "先に arf_parser を実行してARFデータを読み込んでください。")
    class_index = session_state.session.arf.class_index or {}
    class_counts = class_index.get("class_counts", {})
    # 語彙は**常にデータセット全体**（features）から作る。filtered_features を優先すると
    # 直前の arf_parser(class_ids=...) の絞り込みで語彙が黙って痩せ、「何で絞れるか」を
    # 尋ねる本ツールが「control サンプルは存在しない」と誤答する（class_counts は
    # mddata 全体を見ているため、同じ payload の中で2つのスコープが混在してしまう）。
    names = sample_factors.arf_sample_names(session_state.session.arf.features)
    facets = sample_factors.build_sample_facets(
        names, session_state.session.arf.class_index)
    return json_payload({
        "mddata_path": class_index.get("mddata_path"),
        "class_counts": class_counts,
        "factors_by_position": _class_factors_by_position(class_counts.keys()),
        "sample_token_vocabulary": sample_factors.token_vocabulary(facets),
    })




@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_list_sample_roles() -> str:
    """ロード済み ARF のサンプルを sample/qc/blank に分類し、TSV 表で返す（前処理の適用前確認）。

    列は sample / role / group / batch / run_order / excluded。`group` は
    group_levels 未指定時の既定＝完全 Class ID。`batch` は MS-DIAL に専用項目が
    無いためファイル名中の8桁日付からの推定で、その出所（`batch_source`）は
    全行共通なのでヘッダ行にまとめてある。
    """
    if session_state.session.arf.filtered_features is None:
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"],
            "先に arf_parser で ARF を読み込んでください。")
    _, sample_names, _ = tool_helpers._pp_build_matrix(session_state.session.arf.filtered_features, ["height"])
    meta = _build_sample_meta(sample_names, session_state.session.arf.class_index)
    counts = {"sample": 0, "qc": 0, "blank": 0}
    for m in meta.values():
        counts[m["role"]] = counts.get(m["role"], 0) + 1
    # 同じ 6 キーをサンプル数ぶん繰り返す JSON より、列名 1 回の TSV のほうが安い
    # （実データ 60 サンプルで 12,418 字 → 約 2,600 字）。
    rows = [
        {
            "sample": name,
            "role": m["role"],
            "group": m["group"],
            "batch": m["batch"],
            "run_order": m["run_order"],
            "excluded": name in session_state.session.arf.excluded_samples,
        }
        for name, m in meta.items()
    ]
    sources = {m["batch_source"] for m in meta.values() if m["batch_source"]}
    header = (
        f"# サンプル役割一覧: sample={counts['sample']} / qc={counts['qc']} / blank={counts['blank']}"
        f"（計 {len(rows)} 件）\n"
        f"# batch の出所: {', '.join(sorted(sources)) or '不明（ファイル名に日付なし）'}"
    )
    return f"{header}\n{format_spots_as_table(rows)}"


def _resolve_exclude_specs(req_samples, avail_samples):
    """除外指定を実サンプル名へ解決する（完全一致優先、次に因子トークン spec）。

    完全一致を先に見るのは、サンプル名が他のサンプル名のトークン部分集合に
    なっている場合でも「その1件だけ」という従来の意図を保つため。
    一致ゼロの指定は例外にせず unmatched へ回す（既存方針: 除外指定の取りこぼし
    だけで解析全体を止めない）。
    """
    if not req_samples:
        return [], {}, []
    facets = sample_factors.build_sample_facets(
        sorted(avail_samples), session_state.session.arf.class_index)
    matched: list[str] = []
    resolved: dict[str, list[str]] = {}
    unmatched: list[str] = []
    for spec in req_samples:
        if spec in avail_samples:
            resolved[spec] = [spec]
            matched.append(spec)
            continue
        try:
            hits, _ = sample_factors.expand_sample_specs(
                [spec], facets, include_roles=None)
        except ValueError:
            unmatched.append(spec)
            continue
        if spec not in hits:
            # トークンが空集合になる spec（空文字・空白・アンダースコアのみ）は
            # expand_sample_specs 内で continue されるため matches に入らない。
            # ValueError と同じ扱いで unmatched へ回す（KeyError で落とさない）。
            unmatched.append(spec)
            continue
        resolved[spec] = hits[spec]
        matched.extend(hits[spec])
    return matched, resolved, unmatched


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_exclude(
    exclude_samples: list[str] | None = None,
    exclude_spots: list[int] | None = None,
    mode: str = "add",
) -> str:
    """PCA 外れサンプルや特定ピークを名前/ID で手動除外・再包含する（可逆・非破壊）。

    先に arf_parser で ARF を読み込んでおくこと。除外は session に保持され、以降の
    arf_parser の再実行 / arf_preprocess（→ arf_pca_preprocessed / arf_differential）へ反映される。
    filtered_features 自体は変更しないため、mode="remove"/"clear" で元に戻せる。

    引数:
    - exclude_samples: 除外するサンプル。完全なサンプル名に加え、`_` 区切りの因子
      トークン指定（`"ILG_6h"` = ILG かつ 6h の全サンプル）を受ける。トークンは
      Class ID とサンプル名の両方から解決される。解決結果は payload の
      resolved_samples に必ず出る。一致ゼロの指定は unmatched_samples 行きで、
      除外指定の取りこぼしだけで解析を止めない。
    - exclude_spots: 除外するスポットの MasterAlignmentID（int）のリスト。
    - mode: add（既定・追加）/ remove（再包含）/ clear（全消去）/ list（現状表示のみ）。
    """
    spots = session_state.session.arf.filtered_features
    if spots is None:
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"],
            "先に arf_parser で ARF を読み込んでください。")

    avail_samples, avail_ids = exclusions.roster(spots)
    es = session_state.session.arf.excluded_samples
    esp = session_state.session.arf.excluded_spots
    req_samples = list(exclude_samples or [])
    req_spots = list(exclude_spots or [])
    resolved_samples: dict[str, list[str]] = {}
    unmatched_samples: list[str] = []
    unmatched_spots: list[int] = []
    caveats: list[str] = []

    if mode == "clear":
        es.clear()
        esp.clear()
    elif mode == "list":
        pass
    elif mode in ("add", "remove"):
        matched_samples, resolved_samples, unmatched_samples = _resolve_exclude_specs(
            req_samples, avail_samples)
        matched_spots = [i for i in req_spots if i in avail_ids]
        unmatched_spots = [i for i in req_spots if i not in avail_ids]
        if mode == "add":
            es.update(matched_samples)
            esp.update(matched_spots)
        else:  # remove
            es.difference_update(matched_samples)
            esp.difference_update(req_spots)
        if unmatched_samples:
            preview = ", ".join(sorted(avail_samples)[:10])
            caveats.append(
                f"未一致サンプル {unmatched_samples} は現データに存在しません（無視）。"
                f"利用可能サンプル例: {preview}")
        if unmatched_spots:
            caveats.append(
                f"未一致スポット {unmatched_spots} は現データに存在しません（無視）。")
    else:
        return json_payload({"status": "error",
                           "message": f"unknown mode: {mode!r}（add/remove/clear/list）"})

    pruned = exclusions.prune_spots(spots, es, esp)
    pruned_names, pruned_ids = exclusions.roster(pruned)
    payload = {
        "status": "success",
        "mode": mode,
        "excluded_samples": sorted(es),
        "excluded_spots": sorted(esp),
        "resolved_samples": resolved_samples,
        "samples_before": len(avail_samples),
        "samples_after": len(pruned_names),
        "spots_before": len(avail_ids),
        "spots_after": len(pruned_ids),
        "unmatched_samples": unmatched_samples,
        "unmatched_spots": unmatched_spots,
        "caveats": caveats,
    }
    if pruned_names == set() or pruned_ids == set():
        payload["caveats"].append(
            "除外の結果、残サンプルまたは残スポットが 0 件です。PCA/差次的解析は実行できません。")
    return json_payload(payload)


def _manual_exclusion_caveat(n_excl_samples: int, n_excl_spots: int) -> str:
    """手動除外（arf_exclude）の状態を常に明示するcaveat文を返す。

    非ゼロのときだけ発話すると、arf_parser/load_dataset の再実行で
    excluded_samples/excluded_spots が reset_analysis() によって黙って
    ゼロ化されたとき、開示自体も一緒に消えてしまう。ユーザーが外れQC等を
    除外したつもりのまま、外れ値込みの統計が「同じ解析」として出てくる
    事故を防ぐため、ゼロ件のときも「現在なし」を積極的に述べる。
    """
    if n_excl_samples or n_excl_spots:
        return f"ユーザ手動除外: サンプル {n_excl_samples} 件 / スポット {n_excl_spots} 件を除外済み。"
    return "ユーザ手動除外: 現在なし（サンプル0件・スポット0件。arf_exclude による除外は適用されていません）。"


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_preprocess(
    normalize: str = "none",
    blank_min_fold: float | None = None,
    drift_correct: bool = False,
    max_qc_rsd: float | None = None,
    impute: str = "half_min",
    props: list[str] | None = None,
) -> str:
    """ロード済み ARF 行列に前処理レシピを適用し、session を更新して報告を返す。

    以降の PCA/差次的解析は session_state.session.feature_matrix（前処理後）を消費する。
    """
    if session_state.session.arf.filtered_features is None:
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"],
            "先に arf_parser で ARF を読み込んでください。")
    props = props or ["height"]
    # 手動除外（PCA 外れサンプル / 特定ピーク）を行列構築前に適用（非破壊）
    active = exclusions.prune_spots(
        session_state.session.arf.filtered_features,
        session_state.session.arf.excluded_samples,
        session_state.session.arf.excluded_spots,
    )
    matrix, sample_names, feature_names = tool_helpers._pp_build_matrix(active, props)
    meta = _build_sample_meta(sample_names, session_state.session.arf.class_index)
    roles = {n: meta[n]["role"] for n in sample_names}
    run_order = {n: meta[n]["run_order"] for n in sample_names}

    # プールQC が層別（複数 QC サブグループ）かの簡易警告材料
    qc_batches = {meta[n]["batch"] for n in sample_names if meta[n]["role"] == "qc"}
    # バッチ(日付)だけでなく QC 試料名の層別（部位別 QC 等）も検出する。
    qc_strata = preprocessing.detect_qc_strata(sample_names, roles)

    recipe = {
        "normalize": normalize,
        "blank_min_fold": blank_min_fold,
        "drift_correct": drift_correct,
        "max_qc_rsd": max_qc_rsd,
        "impute": impute,
        "props": props,
    }
    matrix2, kept_idx, report = preprocessing.preprocess(
        matrix, sample_names, roles, run_order, recipe,
    )
    kept_feature_names = [feature_names[i] for i in kept_idx]

    # ブランクは背景除去（blank_filter）の参照として使い終えたので、ここで解析行列から
    # 外す。残すと生体試料と桁違いに低い総強度が PCA の PC1 を支配し、群分離の解釈が
    # 壊れる。QC は残す——QC クラスタの締まり具合を PCA で見るのは分析の定番手段。
    # 群に混ざると困る差次的解析側は arf_differential が別途 QC を群から外す。
    matrix2, pp_sample_names, dropped = preprocessing.drop_samples_by_role(
        matrix2, sample_names, roles, drop_roles=("blank",),
    )
    report["excluded_from_matrix"] = dropped
    if dropped.get("blank"):
        report.setdefault("caveats", []).append(
            f"ブランク {len(dropped['blank'])} 件（{', '.join(dropped['blank'])}）は背景除去に"
            "使用後、解析行列（PCA/差次的解析）から除外しました。QC は PCA での品質確認の"
            "ため残しています。"
        )

    session_state.session.arf.feature_matrix = matrix2
    session_state.session.arf.pp_sample_names = pp_sample_names
    session_state.session.arf.pp_feature_names = kept_feature_names
    session_state.session.arf.sample_meta = meta
    session_state.session.arf.preprocessing_recipe = recipe
    if len(qc_batches) > 1:
        report.setdefault("caveats", []).append(
            "プールQC が複数バッチ/層に分かれています。全体一律のドリフト補正は近似です。"
        )
    if len(qc_strata) > 1:
        labels = ", ".join(sorted(s for s in qc_strata if s))
        report.setdefault("caveats", []).append(
            f"プールQC が層別（{len(qc_strata)} サブグループ{f': {labels}' if labels else ''}）"
            "と検出されました。全 QC を1系列として扱うドリフト補正/RSD フィルタは近似です。"
        )
    n_excl_s = len(session_state.session.arf.excluded_samples)
    n_excl_p = len(session_state.session.arf.excluded_spots)
    report.setdefault("caveats", []).append(_manual_exclusion_caveat(n_excl_s, n_excl_p))
    # 除外が過度で行列が空（残サンプル0 または 残特徴量0）になった場合を前景化する。
    if matrix2.size == 0:
        report.setdefault("caveats", []).append(
            "前処理後の行列が空です（残サンプルまたは残特徴量が 0 件）。手動除外が過度な"
            "可能性があります。PCA/差次的解析は実行できません（arf_exclude の mode=remove/clear で復帰）。")
    report["status"] = "success"
    report["matrix_shape"] = list(matrix2.shape)
    report["recipe"] = recipe
    return json_payload(report)




@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_pca_preprocessed(
    components: int | None = None,
    top_features: int = 10,
    log_transform: bool = False,
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
) -> str:
    """arf_preprocess で用意した前処理後行列で PCA を実行する（前処理後経路の PCA 入口）。

    正規化・QC フィルタ・欠損補完を経た行列に対する PCA。生スポットへ直接フィルタして
    やり直す PCA は arf_parser。手順は arf_preprocess → 本ツール。色分け(group_levels)や
    log 変換だけを変えて再実行しても、前処理はやり直さない。
    群分けは group_levels（1因子）または group_factors（因子軸の直積、例
    [["control","ILG"],["0h","6h"]] → "control|0h"）で指定する。値トークンは Class ID
    とサンプル名の両方から解決されるため、Class ID に無い時点などでも色分けできる。
    """
    if not _pp_has_preprocessed():
        return mcp_errors.missing_state(
            "preprocessed_matrix", ["arf_preprocess"],
            "前処理後の行列がありません。先に arf_preprocess を実行してください。")
    from lipidmix.arf.reader import run_pca, get_pca_loading_features
    matrix = session_state.session.arf.feature_matrix
    sample_names = session_state.session.arf.pp_sample_names
    feature_names = session_state.session.arf.pp_feature_names
    try:
        pca_result = run_pca(matrix, n_components=components, log_transform=log_transform)
    except Exception as exc:
        return f"[ERROR] PCA 実行に失敗しました: {exc}"

    sample_groups = assign_sample_groups(
        sample_names, session_state.session.arf.class_index, group_levels,
        group_factors=group_factors,
    )
    plot_block = _format_pca_plot_block(
        pca_result, sample_names,
        title="PCA (preprocessed ARF)",
        intro="\n#### 📊 PCA スコアプロット用データ（前処理後）\n",
        groups=sample_groups,
    )
    _remember_arf_pca_plot(pca_result, sample_names,
                           title="PCA (preprocessed ARF)", groups=sample_groups)
    loading_features = get_pca_loading_features(
        pca_result, session_state.session.arf.features or [], feature_names, top_n=top_features,
    )
    loadings_block = _format_pca_loadings_md(
        loading_features, header="#### 📊 PCA Loadings 寄与度分析（前処理後）\n",
    )
    text = (
        f"### 📈 前処理後 ARF PCA 解析\n"
        f"- **前処理レシピ**: {session_state.session.arf.preprocessing_recipe}\n"
        f"- **PCA入力行列の形状**: {tuple(matrix.shape)} (サンプル数 x 特徴量数)\n"
        f"- **PC1 説明分散比**: {pca_result['explained_variance_ratio'][0]*100:.2f}%\n"
        f"- **PC2 説明分散比**: {pca_result['explained_variance_ratio'][1]*100:.2f}%\n"
        f"{loadings_block}{plot_block}"
    )
    return text


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_parser(
    file_path: str | None = None,
    props: list[str] | None = None,
    components: int | None = None,
    top_features: int = 10,
    log_transform: bool = False,
    min_detection_rate: float = 0.0,
    min_intensity: float = 0.0,
    annotation_keyword: str | None = None,
    tag_labels: list[str] | None = None,
    tag_mode: str = "any",
    tag_scope: str = "sample_peak",
    tag_directory: str | None = None,
    missing_sample_policy: str = "error",
    class_ids: list[str] | None = None,
    class_missing_sample_policy: str = "error",
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
    include_roles: list[str] | None = None,
    spot_ids: list[StrictInt] | None = None,
    ontologies: list[str] | None = None,
) -> str:
    """.arf（サンプル別強度）を読み込み、フィルタを適用して PCA を実行する。

    生スポット行列に対する PCA の唯一の入口。**フィルタ条件を変えて PCA をやり直したい
    ときは、引数を変えて本ツールを再呼び出しする**（ファイルはセッションキャッシュされ
    再パースは走らない）。手動除外（arf_exclude）も PCA 前に反映される。
    正規化・QC・欠損補完を経た「前処理後」行列での PCA は arf_preprocess → arf_pca_preprocessed。

    正負の Loading 上位とスコアプロット用 JSON を含む Markdown 要約を返す。

    引数:
    - file_path: 解析する .arf ファイルのパス (省略時は自動検索)
    - props: PCAに使用するプロパティのリスト (デフォルト: ["height"])
    - components: 計算する主成分の数
    - top_features: 各主成分から抽出する正・負の寄与トップ件数 (デフォルト: 10)
    - log_transform: [任意] PCA前に log10 変換を適用する（強度の歪みを抑え条件分離が向上しやすい。既定 False）
    - min_detection_rate: [任意] 特徴量の実検出率(非ギャップフィル)による足切り 0.0-1.0（既定 0.0=無効）
    - min_intensity: [任意] スポット平均強度(HeightAverage)の最小閾値（既定 0.0=無効）
    - annotation_keyword: [任意] 同一バッチARF2の名前を優先した部分一致検索。
      ARF2が無い/未注釈ならARF名を使う。PCはLPCにも部分一致するため、クラス限定はontologiesを使う。
    - spot_ids: [任意] MasterAlignmentIDの完全一致リスト（0以上、空リスト不可）。
    - ontologies: [任意] 脂質クラスの完全一致リスト（大文字小文字を区別しない）。
      spot_ids、キーワード、他のフィルタとはAND。元のピーク値は変更しない。
      これは生スポットの選択であり、選択変更後はarf_preprocessを再実行する。
      部分集合でTIC正規化すると分母が変わる。全体の正規化/FDRを維持する用途では
      全体で解析・エクスポートしてから対象を参照し、部分集合で再計算しない。
    - tag_labels: [任意] MS-DIALタグ名またはタグIDのリスト
    - tag_mode: any/all/none/not_all のいずれか
    - tag_scope: sample_peak（サンプル別Peak ID）または alignment_spot（MasterAlignmentID）
    - tag_directory: [任意] *_tags.xml の探索先。既定はARFと同じディレクトリ
    - missing_sample_policy: タグファイル未対応サンプルの扱い。error/exclude/untagged（既定 error）
    - class_ids: [任意] 因子トークンによるサンプル絞り込み。各要素は `_` 区切りの部分指定
      で、指定した全トークンを含むサンプルに一致する（要素内AND・順不同）。要素間はOR。
      トークンは **Class ID とサンプル名の両方**から解決されるため、Class ID に入って
      いない因子（時点・複製・測定日）でも絞れる。例: `["6h"]`, `["ILG_6h","control_6h"]`。
    - class_missing_sample_policy: Class IDメタデータ未対応サンプルの扱い。error/exclude（既定 error）
    - group_levels: [任意] PCA点の色分け因子の値トークン（例: `["gf","spf"]`）。
      未指定なら各サンプルの完全Class IDで色分け。指定するとその因子だけで統合し、
      該当しないサンプルは "other" 群になる。
    - group_factors: [任意] 因子軸のリスト。各軸は値トークンのリストで、直積ラベルを
      作る（例 `[["control","ILG"],["0h","6h"]]` → `"control|0h"`）。group_levels より優先。
      処置×時点のような多群を手書き列挙せずに色分けできる。
    - include_roles: [任意] class_ids のフィルタに含める role（既定 `["sample"]`）。
      統合トークン空間ではサンプル名経由で QC/blank を掴み得るため既定で落とす。
      QC も含めたいときだけ `["sample","qc"]` のように明示する。
    """
    props = props or ["height"]

    file_path = resolve_arf_file_path(file_path)
    if not file_path:
        return "データディレクトリに .arf ファイルが見つかりませんでした。"

    # 外部モジュールからのインポート
    from lipidmix.arf.reader import count_peak_property_rows, build_pca_matrix, run_pca, get_pca_loading_features

    try:
        # 不正な選択引数でファイルや既存の解析状態を切り替えない。
        path_resolvers._filter_arf_spots([], spot_ids=spot_ids, ontologies=ontologies)
        deserialized_and_formatted_data = session_state.session.arf.load_data(file_path, tag_directory=tag_directory)
        if not isinstance(deserialized_and_formatted_data, list):
            return "デシリアライズ結果がリストではありません。"

        analysis_data, tag_filter_stats = filter_arf_by_tags(
            deserialized_and_formatted_data,
            session_state.session.arf.tag_index or {},
            tag_labels,
            mode=tag_mode,
            scope=tag_scope,
            missing_sample_policy=missing_sample_policy,
        )
        if not analysis_data:
            return "指定されたMS-DIALタグ条件に一致するARFピークが見つかりませんでした。arf_list_tags で利用可能タグと件数を確認してください。"

        analysis_data, class_filter_stats = filter_arf_by_class_ids(
            analysis_data,
            session_state.session.arf.class_index,
            class_ids,
            missing_sample_policy=class_missing_sample_policy,
            include_roles=tuple(include_roles) if include_roles else ("sample",),
        )
        if not analysis_data:
            return "指定されたClass IDに一致するARFサンプルが見つかりませんでした。arf_list_classes で利用可能なClass IDと件数を確認してください。"

        # 統合注釈を用いる際も他バッチのARF2や別セッションのカタログは参照しない。
        selection_requested = bool(annotation_keyword) or spot_ids is not None or ontologies is not None
        catalog = None
        sibling = _sibling_arf2_path() if selection_requested else None
        if sibling:
            from lipidmix.arf2.reader import load_catalog
            catalog = load_catalog(sibling)
        analysis_data = path_resolvers._filter_arf_spots(
            analysis_data, min_intensity, annotation_keyword,
            catalog=catalog, spot_ids=spot_ids, ontologies=ontologies)
        if not analysis_data:
            return (
                f"指定された条件（強度 >= {min_intensity}, キーワード: '{annotation_keyword or '指定なし'}', "
                f"spot_ids={spot_ids}, ontologies={ontologies}）"
                "に一致するARFピークが見つかりませんでした。"
            )

        # 手動除外（PCA 外れサンプル / 特定ピーク）を PCA 前に適用（非破壊）
        active_data = exclusions.prune_spots(
            analysis_data,
            session_state.session.arf.excluded_samples,
            session_state.session.arf.excluded_spots,
        )

        # 使うのは行数だけなので DataFrame は組まない（spot × 注入ぶんの dict を
        # pandas へ積むと、393MB の `.arf` で 881 万行になり数十秒を捨てる）。
        peak_row_count = count_peak_property_rows(active_data)
        avg_samples = 0
        if len(active_data) > 0 and peak_row_count > 0:
            avg_samples = peak_row_count / len(active_data)

        # PCA行列構築（min_detection_rate は任意の検出率フィルタ）
        matrix, sample_names, feature_names = build_pca_matrix(
            active_data, use_properties=props,
            min_detection_rate=min_detection_rate,
        )

        if matrix.size == 0 and not (selection_requested and sample_names):
            return "[ERROR] PCA 用データを構築できませんでした（フィルタ・除外が過度な可能性があります）。"

        selection_note = _arf_selection_note(analysis_data, sibling, spot_ids, ontologies) if selection_requested else ""
        if selection_requested and min(matrix.shape) < 2:
            _commit_arf_selection(analysis_data, selection_requested)
            return session_state.session.maybe_prepend_caveat(
                f"### ARF脂質選択完了: {Path(file_path).name}\n{selection_note}"
                f"- **適用フィルタ条件**: 強度最小値=`{min_intensity}`, アノテーションキーワード=`'{annotation_keyword or '指定なし'}'`\n"
                f"{_format_arf_class_filter(class_filter_stats)}"
                f"{_format_arf_tag_filter(tag_filter_stats)}"
                f"- **PCA入力行列の形状**: {matrix.shape} (サンプル数 x 特徴量数)\n"
                f"- PCAはスキップ: 分散・検出率フィルタ後に2サンプル・2特徴以上が必要です（行列形状: {matrix.shape}）。\n"
                "- 選択は保存済みです。必要ならarf_preprocessで後続処理へ進めます。\n", topic="arf")

        # PCA実行（log_transform は任意のlog10変換）
        pca_result = run_pca(matrix, n_components=components, log_transform=log_transform)

        _commit_arf_selection(analysis_data, selection_requested)

        # サンプル別の群ラベル（既定=完全Class ID、group_levels 指定時はその因子で統合）
        sample_groups = assign_sample_groups(
            sample_names, session_state.session.arf.class_index, group_levels,
            group_factors=group_factors,
        )

        # PCAスコアプロット用データ（共通ヘルパー）
        plot_instruction_text = _format_pca_plot_block(
            pca_result, sample_names,
            title=f"PCA Score Plot ({Path(file_path).name})",
            intro=(
                "\n#### 📊 PCA スコアプロット用データ\n"
                "以下のJSONデータを用いて、見やすい散布図（Scatter Plot）を描画してください。\n"
                "各点には `sample` の名前をラベルとして表示するか、ホバー時に確認できるようにしてください。\n"
                "`group` フィールドがある場合は、群ごとに色分け（凡例付き）して群間比較が分かるようにしてください。\n"
            ),
            groups=sample_groups,
        )
        _remember_arf_pca_plot(
            pca_result, sample_names,
            title=f"PCA Score Plot ({Path(file_path).name})",
            groups=sample_groups,
        )

        # Loadings 寄与上位（arf_reader の構造化関数 + 共通整形ヘルパー）
        loading_features = get_pca_loading_features(
            pca_result, analysis_data, feature_names, top_n=top_features,
        )
        loadings_summary_text = _format_pca_loadings_md(
            loading_features, header="#### 📊 PCA Loadings 寄与度分析 (各極値トップ件数)\n",
        )

        n_excl_s = len(session_state.session.arf.excluded_samples)
        n_excl_p = len(session_state.session.arf.excluded_spots)
        exclude_note = (
            f"- **ユーザ手動除外**: サンプル {n_excl_s} 件 / スポット {n_excl_p} 件\n"
            if (n_excl_s or n_excl_p) else ""
        )
        filter_note = (
            f"- **適用フィルタ条件**: 強度最小値=`{min_intensity}`, アノテーションキーワード=`'{annotation_keyword or '指定なし'}'`\n"
            if (min_intensity or annotation_keyword) else ""
        )
        filter_note += selection_note

        # 基本的な要約テキストの作成
        output_text = (
            f"### 📈 ARF 多変量PCA解析完了: {Path(file_path).name}\n"
            f"- **読み込んだ総スポット数**: {len(deserialized_and_formatted_data)}\n"
            f"{filter_note}"
            f"{exclude_note}"
            f"{_format_arf_parse_summary(deserialized_and_formatted_data)}"
            f"{_format_arf_class_summary(session_state.session.arf.class_index)}"
            f"{_format_arf_class_filter(class_filter_stats)}"
            f"{_format_arf_tag_summary(session_state.session.arf.tag_index)}"
            f"{_format_arf_tag_filter(tag_filter_stats)}"
            f"- **抽出された総ピークレコード数**: {peak_row_count}\n"
            f"- **平均サンプル数/スポット**: {avg_samples:.2f}\n"
            f"- **PCA入力行列の形状**: {matrix.shape} (サンプル数 x 特徴量数)\n"
            f"- **PC1 説明分散比**: {pca_result['explained_variance_ratio'][0]*100:.2f}%\n"
            f"- **PC2 説明分散比**: {pca_result['explained_variance_ratio'][1]*100:.2f}%\n"
            f"{loadings_summary_text}"
            f"{plot_instruction_text}"  # ← 座標ブロックは末尾（loadings の後）へ
        )

        return session_state.session.maybe_prepend_caveat(output_text, topic="arf")

    except MemoryError as e:
        return _format_arf_memory_error(
            file_path, e, getattr(session_state.session.arf, "features", None),
        )
    except Exception as e:
        return f"[ERROR] ARF解析に失敗しました: {str(e)}"


# 実測モデル: RSS ≈ 8.0 KB × 行数、行数 = スポット × 注入（docs/HISTRY.md 2026-09-20）。
# 123k〜2.94M 行で線形。ファイルのバイト数とは相関しない。
_ARF_BYTES_PER_ROW = 8 * 1024


def _count_arf_rows(features) -> int | None:
    """読み込み済みスポットが抱えているピーク行数を数える（DataFrame を組まない）。"""
    if not isinstance(features, list) or not features:
        return None
    total = 0
    for spot in features:
        aligned = spot.get("AlignedPeakProperties") if isinstance(spot, dict) else None
        if isinstance(aligned, list):
            total += len(aligned)
    return total


def _format_arf_memory_error(file_path: str, exc: BaseException, features) -> str:
    """メモリ不足を、規模と対処つきで返す。

    `str(MemoryError())` は空文字なので、そのまま返すと理由なしのエラーになる。
    400 検体規模で最初に当たる失敗がこれ。
    """
    spots = getattr(exc, "spots", None)
    rows = getattr(exc, "rows", None)
    if rows is not None:
        stage = "読み込み中に尽きた"
    elif (counted := _count_arf_rows(features)) is not None:
        spots, rows = len(features), counted
        stage = "読み込みは終わり、解析中に尽きた"
    else:
        stage = "どこで尽きたかを特定できなかった"

    if rows is None:
        scale = "- **読めた範囲**: 不明\n"
        needed = ("- **必要な RAM**: 行数が取れないため算出できない"
                  "（実測モデルは 8.0 KB × スポット × 注入）\n")
    else:
        scale = f"- **読めた範囲**: {spots:,} スポット / {rows:,} 行（行数 = スポット × 注入）\n"
        needed = (f"- **必要な RAM の目安**: 約 {rows * _ARF_BYTES_PER_ROW / 1e9:.1f} GB"
                  "（実測モデル 8.0 KB/行）\n")

    return (
        f"[ERROR] ARF解析に失敗しました: メモリ不足（MemoryError）。{Path(file_path).name}\n"
        f"- **打ち切り**: {stage}\n"
        f"{scale}"
        f"{needed}"
        "- **対処**: `class_ids` / `tag_labels` / `spot_ids` で対象を絞るか、"
        "RAM の大きい機械で読む。ファイルのバイト数は指標にならない"
        "——効くのはスポット × 注入。\n"
    )


def _commit_arf_selection(spots: list[dict], selection_requested: bool) -> None:
    """選択確定時に、異なるサンプル・注釈に由来する解析結果を無効化する。"""
    def signature(items):
        return [(s.get("MasterAlignmentID"), exclusions.roster([s])[0],
                 s.get("Name"), s.get("Ontology"), s.get("annotation_source")) for s in items]

    state = session_state.session.arf
    if selection_requested or signature(state.filtered_features or []) != signature(spots):
        state.feature_matrix = None
        state.pp_sample_names = None
        state.pp_feature_names = None
        state.sample_meta = {}
        state.preprocessing_recipe = {}
        state.last_differential = None
    state.last_pca_plot = None
    state.pca_result = None
    state.filtered_features = spots


def _arf_selection_note(spots, sibling, spot_ids, ontologies) -> str:
    """選択に使った統合注釈の出所と正規化の注意を表示する。"""
    sources = {s.get("annotation_source", "arf") for s in spots}
    conflicts = [s.get("MasterAlignmentID") for s in spots if s.get("annotation_conflict")]
    return (
        f"- **脂質選択**: spot_ids={spot_ids}, ontologies={ontologies}; 採用 {len(spots)} 件\n"
        f"- **注釈の出所**: {', '.join(sorted(sources))}; 同一バッチARF2={sibling or 'なし'}\n"
        f"- **ARF/ARF2名前の不一致**: {len(conflicts)} 件（ID先頭10件: {conflicts[:10]}）\n"
        "- 生スポットを選択しました。前処理・差次解析は再実行が必要です。"
        "部分集合でのTIC正規化やFDRは全体解析と分母・対象数が変わります。\n")


_SPOT_ID_RE = re.compile(r"^Spot_(\d+)")


def _pool_group_labels(sample_names, group_labels, group_a, group_b):
    """完全 Class ID / 群ラベルだけでなく、サンプル名の因子トークンでもプール群を作る。

    Class ID は MS-DIAL 上で入力された1文字列にすぎず、時点や複製のような因子は
    サンプル名にしか無いことがある。トークン空間を tokens(サンプル名) ∪ tokens(群ラベル)
    に統合し、`group_a="ILG_6h"` のような多因子指定を通す。完全な群ラベルを渡した
    場合は（そのトークン集合を含む他のサンプルが無い限り）従来どおりの2群比較になる。

    group_labels が None のサンプル（QC/blank・群未解決）は候補から外す。
    戻り値 (relabeled, resolved, resolved_samples)。一致ゼロ・両群の重複は ValueError。
    """
    if str(group_a) == str(group_b):
        raise ValueError(f"group_a と group_b が同一です: {group_a!r}")
    available = sorted({str(g) for g in group_labels if g is not None})
    if not available:
        raise ValueError(
            "Class ID メタデータが解決できていないため群を特定できません（.mddata 未検出）。")

    meta = {name: {"group": label, "role": "sample"}
            for name, label in zip(sample_names, group_labels) if label is not None}
    facets = sample_factors.build_sample_facets(list(meta), None, sample_meta=meta)
    try:
        # role は呼び出し側が group_labels=None で既に落としているので、ここでは絞らない。
        matches, _ = sample_factors.expand_sample_specs(
            [group_a, group_b], facets, include_roles=None)
    except ValueError as exc:
        raise ValueError(f"{exc} 利用可能な群ラベル: {', '.join(available)}") from exc

    # トークンが空集合になる spec（空文字・空白のみ・アンダースコアのみ）は
    # expand_sample_specs 内で読み飛ばされ matches に入らない。ここで ValueError へ
    # 変換しないと KeyError が MCP ツールの外まで抜ける（他の失敗経路はすべて構造化
    # エラー JSON を返す）。_resolve_exclude_specs の同名ガードと同じ失敗モード。
    blank_specs = [str(spec) for spec in (group_a, group_b) if spec not in matches]
    if blank_specs:
        raise ValueError(
            f"群指定 {blank_specs} に因子トークンがありません"
            "（空文字・空白のみ・アンダースコアのみは群指定として無効です）。"
            f"利用可能な群ラベル: {', '.join(available)}")

    a_names, b_names = set(matches[group_a]), set(matches[group_b])
    overlap = sorted({str(meta[n]["group"]) for n in (a_names & b_names)})
    if overlap:
        raise ValueError(
            f"group_a='{group_a}' と group_b='{group_b}' が同じサンプルを含みます "
            f"({', '.join(overlap)})。群は排他である必要があります。")

    relabeled = []
    for name, label in zip(sample_names, group_labels):
        if label is None:
            relabeled.append(None)
        elif name in a_names:
            relabeled.append(group_a)
        elif name in b_names:
            relabeled.append(group_b)
        else:
            relabeled.append(None)
    resolved = {
        "group_a": sorted({str(meta[n]["group"]) for n in a_names}),
        "group_b": sorted({str(meta[n]["group"]) for n in b_names}),
    }
    resolved_samples = {"group_a": sorted(a_names), "group_b": sorted(b_names)}
    return relabeled, resolved, resolved_samples


def _spot_id_of(feature_name) -> int | None:
    match = _SPOT_ID_RE.match(str(feature_name))
    return int(match.group(1)) if match else None


def _sibling_arf2_path() -> Path | None:
    """読み込み中の ARF と同一アラインメント実行の .arf2 を返す。

    MasterAlignmentID はアラインメント実行ごとに振り直されるため、別バッチの
    .arf2 を引くと ID 対応が黙って崩れる。GUI/Consoleの既知の語幹が
    完全一致する兄弟だけを許し、無ければ None（注釈は諦める）。
    """
    current = getattr(session_state.session.arf, "current_file_path", None)
    if not current:
        return None
    path = Path(current)
    match = re.fullmatch(
        r"(AlignmentResult_\d{4}(?:_\d{2}){5}|AlignResult-\d+)"
        r"(?:_(?:PeakProperties|DriftSpots|DriftSopts))?\.arf", path.name)
    if not match:
        return None
    sibling = path.with_name(f"{match.group(1)}.arf2")
    return sibling if sibling.is_file() else None


def _annotate_with_names(rows: list[dict]) -> dict:
    """差次的解析の上位ヒットに脂質名/Ontology を付す。ARF が Unknown なら ARF2 を引く。

    同一アラインメントの ARF と ARF2 は同じ MasterAlignmentID を指すが、代表 Name は
    食い違うことがある（ARF 側 Unknown・ARF2 側は注釈あり）。上位ヒットが
    `Spot_474_height` のままだと解釈に到達できないため橋渡しし、出所を name_source で
    開示する。ARF2 は大きいので、ARF で埋まらない ID が残るときだけ読む。
    """
    spots = {s.get("MasterAlignmentID"): s for s in (session_state.session.arf.features or [])}
    unresolved: list[dict] = []
    for row in rows:
        sid = _spot_id_of(row.get("feature"))
        row["spot_id"] = sid
        name = (spots.get(sid) or {}).get("Name") or ""
        if name and name.strip().lower() != "unknown":
            row["name"] = name
            row["name_source"] = "arf"
        else:
            row["name"] = None
            row["name_source"] = None
            if sid is not None:
                unresolved.append(row)
    report = {"annotated_from_arf": len(rows) - len(unresolved)}
    if not unresolved:
        return report

    arf2_path = _sibling_arf2_path()
    if not arf2_path:
        report["arf2_lookup"] = "skipped: 同一アラインメントの .arf2 が隣接していません"
        return report
    from lipidmix.arf2.reader import load_catalog
    catalog = {s.get("MasterAlignmentID"): s for s in load_catalog(arf2_path)}
    filled = 0
    for row in unresolved:
        spot = catalog.get(row["spot_id"])
        if not spot:
            continue
        name = spot.get("Name") or ""
        if name and name.strip().lower() != "unknown":
            row["name"] = name
            row["ontology"] = spot.get("Ontology") or None
            row["name_source"] = "arf2"
            filled += 1
    report["annotated_from_arf2"] = filled
    report["arf2_path"] = str(arf2_path)
    if filled:
        report["note"] = (
            f"{filled} 件は ARF 側 Name が Unknown で、ARF2 カタログの注釈を採用しました"
            "（name_source=arf2）。ARF と ARF2 で代表 Name は食い違い得ます。")
    return report


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_differential(
    group_a: str | None = None,
    group_b: str | None = None,
    q_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    log_transform: bool = True,
) -> str:
    """前処理後行列で2群差次的解析（Welch t 検定 + log2FC）を行う。

    group_a / group_b の両方を指定する。多群 ANOVA は「因子（加齢/菌叢等）→水準」の
    対応が MS-DIAL メタに無く安全に導けないため現状は非対応（誤って全 Class ID を
    水準にした結果を返さないよう封鎖している）。3群以上を比べたいときは、下記の
    因子トークン・プール指定で関心のある2群を切り出して呼ぶこと。

    先に arf_preprocess を実行して session_state.session.arf.feature_matrix を用意すること
    （未実行ならエラーを返す。生行列への暗黙フォールバックはしない）。

    - group_a / group_b: 完全な Class ID（`24M_GF_F`）に加え、**因子トークンによる
      プール群指定**を受け付ける。`group_a="24M", group_b="9w"` のように書くと、その
      トークンを含む全 Class ID がプールされる（`24M_GF_F` + `24M_SPF_M` + …）。
      複数トークンの AND 指定も可（`"24M_GF"`）。両群が同じ Class ID を掴むとエラー。
    - log_transform: [既定 True] log2(x+1) 空間で検定する。MS 強度は対数正規に近く、
      生強度での t 検定/ANOVA は正規性仮定を外れやすいため既定で有効。2群では log2FC も
      log2 空間の群平均差（＝幾何平均比）になる。生スケールで検定したい場合のみ False。
    """
    matrix = getattr(session_state.session.arf, "feature_matrix", None)
    if matrix is None:
        return mcp_errors.missing_state(
            "preprocessed_matrix", ["arf_preprocess"],
            "先に arf_preprocess を実行してください（前処理後行列が必要）。")
    sample_names = session_state.session.arf.pp_sample_names
    feature_names = session_state.session.arf.pp_feature_names
    meta = session_state.session.arf.sample_meta or {}
    batch_labels = [(meta.get(n) or {}).get("batch") for n in sample_names]

    # 生体試料以外（QC・ブランク）は比較対象から外す。Class ID が group_a/group_b の
    # 因子トークンを含むと（例 QC_24M と group_a="24M"）プールに紛れ込み、群平均を
    # 汚染するため、ラベルを None にして _pool_group_labels のどちらにも寄らせない。
    # ブランクは arf_preprocess で行ごと落ちている想定だが、ここでも二重に守る。
    _NON_SAMPLE_ROLES = {"qc", "blank"}
    excluded_roles: dict[str, list[str]] = {}
    group_labels: list = []
    for name in sample_names:
        entry = meta.get(name) or {}
        role = entry.get("role", "sample")
        if role in _NON_SAMPLE_ROLES:
            excluded_roles.setdefault(role, []).append(name)
            group_labels.append(None)
        else:
            group_labels.append(entry.get("group"))

    caveats: list[str] = []
    caveats.append(
        _manual_exclusion_caveat(
            len(session_state.session.arf.excluded_samples),
            len(session_state.session.arf.excluded_spots),
        )
    )
    recipe = session_state.session.arf.preprocessing_recipe or {}
    if recipe.get("normalize", "none") == "none":
        caveats.append("正規化が未適用のため log2FC は測定量差を含み得ます（arf_preprocess の normalize を検討）。")
    if log_transform:
        caveats.append("log2(x+1) 変換後に検定を実施（強度の歪みを補正）。log2FC は群平均の log2 差＝幾何平均比です。")
    if excluded_roles:
        detail = "; ".join(
            f"{role}={len(names)}件（{', '.join(names)}）"
            for role, names in sorted(excluded_roles.items())
        )
        caveats.append(f"比較対象から除外（生体試料でないため）: {detail}。")

    batch_source = next((m.get("batch_source") for m in meta.values() if m.get("batch_source")), None)
    src_note = "（バッチはファイル名の日付から推定。実バッチ設計と異なる場合あり）" \
        if batch_source == "filename_date" else ""

    if group_a is not None and group_b is not None:
        try:
            group_labels, resolved, resolved_samples = _pool_group_labels(
                sample_names, group_labels, group_a, group_b)
        except ValueError as exc:
            return json_payload({"status": "error", "message": str(exc)})

        # 交絡は「実際に比較した2群」に対して見る。プール前の Class ID 単位で判定すると、
        # 細粒度ラベルほど各群が単一バッチになりやすく偽の交絡警告を出す（プールすれば
        # 両群ともバッチ混在、という設計を交絡と誤報していた）。
        paired = [(g, b) for g, b in zip(group_labels, batch_labels) if g is not None]
        conf = differential.check_confounding(
            [g for g, _ in paired], [b for _, b in paired],
        )
        if conf["confounded"]:
            caveats.append("交絡: " + conf["detail"] + src_note)
        elif not conf.get("assessable", True):
            caveats.append("交絡評価不可: " + conf["detail"] + src_note)
        # プールしているかどうかは Class ID の数ではなくサンプルの数で決まる。
        # 統合トークン空間では 1 つの Class ID の中を時点等で切り分けられるため、
        # len(resolved[...]) > 1 で判定すると「ILG_6h = 3サンプル」のような
        # プールが発火せず、平均化されている事実が隠れる。
        if any(len(names) > 1 for names in resolved_samples.values()):
            caveats.append(
                "プール群として解決: "
                + "; ".join(
                    f"{spec} = {len(names)} サンプル"
                    + (f"（Class ID: {' + '.join(ids)}）" if ids else "")
                    for spec, names, ids in (
                        (group_a, resolved_samples["group_a"], resolved["group_a"]),
                        (group_b, resolved_samples["group_b"], resolved["group_b"]))
                )
                + "。因子内の他要因（性・菌叢等）はプール内で平均化されます。")
        # resolved_class_ids は「選択されたサンプルが持つ Class ID」であって群の定義
        # ではない。両群が同じ Class ID を持つ構成（ILG_6h vs ILG_0h）では、そのまま
        # 読むとツール自身が拒否すべき縮退比較に見える。群を分けている実体は
        # サンプル名の因子トークンなので、その旨を明示する。
        shared_class_ids = sorted(set(resolved["group_a"]) & set(resolved["group_b"]))
        if shared_class_ids:
            caveats.append(
                f"resolved_class_ids は選択サンプルの Class ID であり群の定義では"
                f"ありません。group_a='{group_a}' と group_b='{group_b}' は "
                f"Class ID では区別されず（共通: {', '.join(shared_class_ids)}）、"
                "サンプル名の因子トークンで分離されています。実際の群構成は "
                "resolved_samples を参照してください。")
        caveats.append(
            "比較サンプル: "
            + "; ".join(
                f"{spec} = {', '.join(names)}"
                for spec, names in ((group_a, resolved_samples["group_a"]),
                                    (group_b, resolved_samples["group_b"]))
            )
            + "。指定トークンは Class ID とサンプル名の両方から解決されます。"
        )
        caveats.append(
            f"log2FC の向き: 正なら {group_b} が高い（上昇）、負なら {group_a} が高い（低下）。"
            "2026-08-31 に慣習へ合わせて反転したため、それ以前の出力とは符号が逆である。")
        results = differential.two_group_test(matrix, feature_names, group_labels,
                                              group_a, group_b, log_transform=log_transform)
        results = differential.add_fdr(results)
        summary = differential.summarize_two_group(results, q_threshold, log2fc_threshold)
        annotation = _annotate_with_names(summary.get("top") or [])
        if annotation.get("note"):
            caveats.append(annotation["note"])
        volcano = differential.volcano_data(results, q_threshold, log2fc_threshold)
        n_a = group_labels.count(group_a)
        n_b = group_labels.count(group_b)
        if n_a < 2 or n_b < 2:
            caveats.append(
                f"群サイズ不足（{group_a}={n_a}, {group_b}={n_b}）: 各群 n>=2 が必要です。"
                "群名の誤り、または前処理での試料脱落の可能性があります。")
        elif min(n_a, n_b) < 4:
            caveats.append(f"小n（{group_a}={n_a}, {group_b}={n_b}）につき検出力が限られます。")
        n_tested = summary["n_tested"]
        if n_tested == 0:
            caveats.append(
                "検定可能な特徴が0件（全特徴で p=NaN）。群が空・分散0・または正規化で試料が"
                "NaN化した可能性があります。『有意0件』を『群間差なし』と解釈しないでください。")
        elif n_tested < 0.2 * len(feature_names):
            caveats.append(
                f"検定できた特徴は {n_tested}/{len(feature_names)} 件のみ（多くが p=NaN）。"
                "群内 n 不足・分散0・欠損が多い可能性があります（前処理の見直しを検討）。")
        # contract_version / log2fc_sign は export_contract を単一情報源として直接参照する
        # （リテラルで複製すると CONTRACT_VERSION を上げたときここが追随せず、
        # arf_export_differential の互換性チェックに永遠に落ち続ける再現ループになる。
        # dataset_analysis.run_dataset_differential と同じ理由・同じ直し方）。
        from lipidmix.analysis.result_state import array_fingerprint, new_provenance
        _arf = session_state.session.arf
        _provenance = new_provenance(
            _arf, kind="differential",
            input_fingerprint=array_fingerprint(getattr(_arf, "feature_matrix", None)),
            effective_parameters={"group_a": group_a, "group_b": group_b,
                                  "q_threshold": q_threshold,
                                  "log2fc_threshold": log2fc_threshold})
        session_state.session.arf.last_differential = {"provenance": _provenance,
                                       "kind": "two_group", "a": group_a, "b": group_b,
                                     "n_a": n_a, "n_b": n_b,
                                     "q_threshold": q_threshold,
                                     "log2fc_threshold": log2fc_threshold,
                                     "contract_version": export_contract.CONTRACT_VERSION,
                                     "log2fc_sign": export_contract.LOG2FC_SIGN,
                                     "log_transform": log_transform,
                                     "results": results, "volcano": volcano}
        # 全量 volcano（~特徴数）は上の last_differential に保持し、arf_plot_volcano
        # （構造化点列）と save_volcano_figure（PNG）から使う。payload には載せない
        # ——先頭の summary が巨大 volcano 配列＋文脈切り詰めで埋没し、解釈モデルが
        # 有意件数を読めず「全て ns」と誤読する退行を避けるため。
        payload = {"status": "success", "kind": "two_group",
                   "group_a": group_a, "group_b": group_b,
                   "resolved_class_ids": resolved,
                   "resolved_samples": resolved_samples,
                   "n_a": n_a, "n_b": n_b,
                   "summary": summary, "caveats": caveats,
                   "differential_contract_version": export_contract.CONTRACT_VERSION,
                   "log2fc_sign": ("log2fc は正なら group_b が高い（上昇）。"
                                   "2026-08-31 以前の出力とは符号が逆である。"),
                   "volcano_note": "全特徴の volcano 点列は本要約に非同梱。"
                                   "arf_plot_volcano で構造化した点列を取得し、"
                                   "クライアント側で散布図を描画してください。"
                                   "PNG が必要だとユーザーが明示した場合のみ "
                                   "save_volcano_figure を実行します。"}
    else:
        return json_payload({"status": "error",
                           "message": "group_a と group_b の両方を指定してください（2群比較）。"
                                      "完全 Class ID か因子トークン（例 group_a='24M', group_b='9w'）で"
                                      "関心のある2群を切り出せます。多群 ANOVA は現状非対応です。"})
    return json_payload(payload)


# 契約の実体は lipidmix.analysis.export_contract（leaf・別リポとの契約）。
# ここでは後方互換の別名束縛のみ行う。tests/test_export_contract.py が `is` で
# 同一性を見るため、値をコピーせず同じオブジェクトを束縛する。
_EXPORT_COLUMNS = export_contract.EXPORT_COLUMNS
_DIFFERENTIAL_CONTRACT_VERSION = export_contract.CONTRACT_VERSION
_LOG2FC_SIGN = export_contract.LOG2FC_SIGN
_format_export_number = export_contract.format_number


@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True),
    structured_output=False)
def arf_export_differential(output_path: str) -> str:
    """差次的結果を InChIKey 付きの 1 ファイルへ書き出す（spec §6.1）。

    先に arf_preprocess → arf_differential を実行しておくこと。同一
    アラインメントの兄弟 .arf2 から同定情報を MasterAlignmentID で結合する。
    濃縮解析の背景を保つため、有意な行だけでなく InChIKey が付いた全行を出す。
    """
    last = getattr(session_state.session.arf, "last_differential", None)
    if not last or last.get("kind") != "two_group":
        return mcp_errors.missing_state(
            "two_group_differential", ["arf_differential"],
            "先に arf_preprocess → arf_differential（2群）を実行してください。")
    if (last.get("contract_version") != _DIFFERENTIAL_CONTRACT_VERSION
            or last.get("log2fc_sign") != _LOG2FC_SIGN):
        return mcp_errors.missing_state(
            "compatible_two_group_differential", ["arf_differential"],
            "直近の差次的結果は現行エクスポート契約と互換性がありません。"
            "arf_differential（2群）を再実行してください。")

    arf2_path = _sibling_arf2_path()
    if not arf2_path:
        return mcp_errors.missing_state(
            "sibling_arf2", ["arf_parser"],
            "同一アラインメントの .arf2 が隣接していません。InChIKey を補えないため"
            "書き出しません（空欄のまま出すと、下流でパスウェイ無しと区別できません）。")

    from lipidmix.arf2.reader import load_catalog
    catalog = {spot.get("MasterAlignmentID"): spot
               for spot in load_catalog(arf2_path)}
    rows, report = identity_join.join_identity(last.get("results") or [], catalog)
    if not rows:
        return json_payload({
            "status": "error",
            "message": ("InChIKey が付いた特徴が 0 件のため書き出しません。"
                        "下流のパスウェイ解析に使える背景集合がありません。"),
            "n_features_total": report["n_features_total"],
            "n_with_inchikey": 0,
            "n_unannotated": report["n_unannotated"],
        })

    q_threshold = last.get("q_threshold")
    log2fc_threshold = last.get("log2fc_threshold")

    meta = export_contract.build_meta(
        group_a=last["a"], n_a=last["n_a"],
        group_b=last["b"], n_b=last["n_b"],
        q_threshold=q_threshold, log2fc_threshold=log2fc_threshold,
        log_transform=last.get("log_transform"),
        n_features_total=report["n_features_total"],
        n_with_inchikey=report["n_with_inchikey"],
        n_unannotated=report["n_unannotated"],
        msi_note="# msi_level は .arf2 由来の注釈確度。MS/MS の有無ではない",
        source_lines=[
            f"# source_arf = {getattr(session_state.session.arf, 'current_file_path', '')}",
            f"# source_arf2 = {arf2_path}",
        ],
        preprocess_line=f"# preprocess = {getattr(session_state.session.arf, 'preprocessing_recipe', None)}",
    )

    lines = [*meta, "\t".join(export_contract.EXPORT_COLUMNS)]
    identity_tables = tool_helpers._identity_tables()
    for row in rows:
        identity_name = row["name"]
        if identity_name.strip().lower() == "unknown":
            identity_name = ""
        identity = lipid_identity.build_identity_block(
            {"name": identity_name, "ontology": row["ontology"], "has_msms": False},
            identity_tables,
            mass_error_band="UNKNOWN",
            adduct_band="UNKNOWN",
        )
        lines.append(export_contract.format_row({
            "spot_id": row["spot_id"],
            "name": row["name"],
            "name_source": "arf2",
            "ontology": row["ontology"],
            "inchikey": row["inchikey"],
            "inchikey_source": "arf2",
            "msi_level": identity["msi"]["level"],
            "mz": row["mz"], "rt": row["rt"],
            "log2fc": row["log2fc"],
            "p_value": row["p_value"], "q_value": row["q_value"],
            "mean_a": row["mean_a"], "mean_b": row["mean_b"],
            "significant": export_contract.is_significant(
                q=row["q_value"], log2fc=row["log2fc"],
                q_threshold=q_threshold, log2fc_threshold=log2fc_threshold),
        }))

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return json_payload({
        "status": "success",
        "output_path": str(out),
        "contract_version": last["contract_version"],
        "group_a": last["a"],
        "group_b": last["b"],
        "n_features_total": report["n_features_total"],
        "n_with_inchikey": report["n_with_inchikey"],
        "n_unannotated": report["n_unannotated"],
        "log2fc_sign": "log2fc は正なら group_b が高い（上昇）。",
        "note": ("n_unannotated は注釈が付かず書き出さなかった行数です。"
                 "「変化が無かった」ではなく「調べていない」行です。"),
    })


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_plot_volcano(
    max_points: int = volcano_plot.DEFAULT_MAX_POINTS,
    title: str | None = None,
    output: str | None = None,
) -> list | str:
    """Show the latest two-group volcano — as a PNG image by default.

    Run ``arf_differential`` (two-group) first. Read-only; writes no file.

    ``output="image"`` (the default) renders the plot server-side and returns it as
    an image block plus a one-line caption, so the figure appears directly in the
    chat. ``output="payload"`` instead returns the client-neutral
    ``lipidmix.volcano.v1`` JSON for a client that draws its own interactive chart
    (Use-LLLM's Plotly view); the deployment default can be flipped with the
    ``LIPIDMIX_PLOT_OUTPUT`` environment variable.

    Prefer the image: the coordinate payload costs tens of thousands of tokens and
    an LLM cannot read a point cloud usefully — the numbers you should reason about
    (counts, top hits, caveats) are in ``arf_differential``'s own result. In payload
    mode, ``up``/``down`` points are always kept in full and only ``ns`` points are
    thinned to fit ``max_points``; every count is reported in ``selection``. The
    image always draws every feature.

    Call ``save_volcano_figure`` only when the user wants the PNG written to disk.
    """
    last = getattr(session_state.session.arf, "last_differential", None)
    if not last or last.get("kind") != "two_group" or not last.get("volcano"):
        return mcp_errors.missing_state(
            "differential_result", ["arf_differential"],
            "先に arf_differential（2群比較）を実行してください"
            "（直近の2群差次的解析の volcano データがありません）。")
    try:
        mode = plot_render.resolve_plot_output(output)
    except ValueError as exc:
        return json_payload({"status": "error", "message": str(exc)})
    if mode == plot_render.PAYLOAD:
        # 最小形の JSON 文字列で返す。FastMCP は dict の戻り値を indent=2 で整形する
        # ため、数百点の payload では実測 86,403 字 → 約 5 万字の差が出る。
        # クライアント（Use-LLLM の volcano-plot.js）は content のテキストを parse する。
        return json_payload(volcano_plot.build_volcano_plot_payload(
            last, max_points=max_points, title=title,
        ))

    counts = _volcano_counts(last)
    png = plot_render.figure_to_png(volcano_plot.render_volcano_plot(last, title=title))
    caption = (
        f"Volcano: {last.get('a')} (n={last.get('n_a')}) vs {last.get('b')} (n={last.get('n_b')})"
        f" — 全 {counts['total']} 特徴のうち up={counts['up']} / down={counts['down']} /"
        f" ns={counts['ns']}（検定不能 {counts['nonfinite']} 件は非描画）。"
        f" しきい値 q<{last.get('q_threshold')}, |log2FC|>={last.get('log2fc_threshold')}。"
        " 個々の特徴量名と統計値は arf_differential の結果を参照。"
    )
    return [caption, Image(data=png, format="png")]


def _volcano_counts(last: dict) -> dict:
    """画像のキャプション用に有意/非有意の件数を数える（図から読み取らせない）。"""
    points = list(last.get("volcano") or [])
    drawable = [p for p in points if volcano_plot._is_drawable(p)]
    return {
        "total": len(points),
        "up": sum(1 for p in drawable if p.get("sig") == "up"),
        "down": sum(1 for p in drawable if p.get("sig") == "down"),
        "ns": sum(1 for p in drawable if p.get("sig") not in ("up", "down")),
        "nonfinite": len(points) - len(drawable),
    }
