"""統計エンジン（`metabolomix.analysis.stats`）の対話ツール。

`dataset_preprocess` が前処理（spec §4）を回して行列を作り、`matrix_id` を返す。行列は
`ds.stat_matrices[matrix_id]` に設定・群の上書き・除外・データの指紋と一緒に置く。同じ入力と設定は同じ `matrix_id`。後続の比較・PCA はこの ID を名指しできる。

戻り値は要約だけ（工程ごとの名前・状態・件数）。全量は `ds.stat_matrices` に保持する。
1 回の呼び出しは `reports/stats/` に記録を 1 件追記する（書けないときも前処理は成功し、
`record.written=false` で理由を返す）。
"""
from __future__ import annotations

import itertools
import math
import re
from pathlib import Path

from mcp.types import ToolAnnotations

from metabolomix.analysis import compare_table, stats_records, stats_results
from metabolomix.analysis.stats import preprocess as stats_preprocess
from metabolomix.analysis.stats import pca as stats_pca
from metabolomix.analysis.stats import spec_v2
from metabolomix.analysis.stats.compare import P_ADJUST_METHODS
from metabolomix.analysis.stats.data import from_dataset, sample_table
from metabolomix.analysis.stats_run_summary import factorial_summary
from metabolomix.core import mcp_core
from metabolomix.core.atomic_io import DomainError
from metabolomix.core.mcp_core import mcp
from metabolomix.core.mcp_errors import missing_state, mztab_error, needs_input
from metabolomix.core.serialization import json_payload, json_safe, round_floats, truncate_lists
from metabolomix.tools._holder import resolve_holder
from metabolomix.tools._matrix_store import register_default, store_matrix

__all__ = ["dataset_preprocess", "dataset_compare", "dataset_pca"]

#: details に載せる一覧の上限（長いものは切る）。
_DETAIL_LIMIT = 20
_STEP_LIST_LIMIT = 20


def _config(*, drift_mode="auto", drift_span=0.75, blank_min_ratio=3.0, detection_min_rate=0.5,
            detection_by="group", max_qc_rsd=0.30, max_d_ratio=0.5, normalize="auto",
            collapse_replicates=True, injection_columns=None, loess_degree=1,
            normalization_features=None, normalization_exclude_features=None) -> dict:
    """平らな引数を前処理の設定（`PREPROCESS_DEFAULTS` の形）へ写す。

    足切りの引数が None ならその足切りを行わない（`filters.<name> = None`）。検証はここでは
    しない（`normalize_preprocess_config` が `PREPROCESS_CONFIG_INVALID` で止める）。任意のキー
    （`drift.loess_degree`・`normalize.features` / `exclude_features`）は既定値なら載せない（既定の
    行列の `matrix_id` を、このキーを足す前と同じに保つ）。
    """
    drift = {"mode": drift_mode, "span": drift_span}
    if loess_degree != 1:
        drift["loess_degree"] = loess_degree
    normalize_cfg = {"method": normalize}
    if normalization_features is not None:
        normalize_cfg["features"] = normalization_features
    if normalization_exclude_features is not None:
        normalize_cfg["exclude_features"] = normalization_exclude_features
    return {
        "drift": drift,
        "filters": {
            "blank_ratio": None if blank_min_ratio is None else {"min": blank_min_ratio},
            "detection": (None if detection_min_rate is None
                          else {"by": detection_by, "min_rate": detection_min_rate}),
            "qc_rsd": None if max_qc_rsd is None else {"max": max_qc_rsd},
            "d_ratio": None if max_d_ratio is None else {"max": max_d_ratio},
        },
        "normalize": normalize_cfg,
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
    loess_degree: int = 1,
    normalization_features: list[str] | None = None,
    normalization_exclude_features: list[str] | None = None,
) -> str:
    """統計エンジンの前処理を回して、前処理済みの行列を作る（行列 ID を返す）。

    工程は 選択 → ドリフト補正 → 品質の足切り → 正規化 → log2 → 反復のまとめ（個体ごとの行列）。同じデータ・除外・群・設定は
    同じ `matrix_id`。比較（dataset_compare）・PCA（dataset_pca）・パスウェイ（pathway_project）はこの ID を受ける（省略時は、
    このハンドルで最後に作ったこの行列。無ければ既定の設定）。`dataset_exclude` の除外は自動で反映される。`reports/stats/` に記録を
    1 件追記する（書けないときは record.written=false。前処理は成功）。陽陰を統合したデータセット（dataset_merge_polarities）は、
    ドリフト補正をその極性の注入順・バッチで、正規化の係数をその極性の特徴だけから、極性ごとに当てる。段の `by_polarity` に極性
    ごとの記録（ドリフトの段は一部の極性だけ補正すれば status=partial）、陽陰で注入順・バッチが違えば warnings の
    `polarity_runs_differ`。戻り値・段・注意の意味は `docs_read(topic="dataset")` の 10.2。

    drift_mode: QC の LOESS ドリフト補正。"auto"（既定。条件が揃えば行う）/ "on" / "off"。auto は条件が揃っても、換算済みのとき・
        注入順が MS-DIAL の既定（ファイルの順）かもしれないときは補正せずに計算し、warnings の `drift_not_applied_by_default` に理由と
        `drift_mode="on"` での呼び直しを返す（尋ねない）。"on" は条件が揃えば必ず補正し、揃わなければ PREPROCESS_PRECONDITION。
    drift_span: LOESS の窓（0〜1、既定 0.75）。
    loess_degree: LOESS の局所多項式の次数 1（既定。局所線形）/ 2（局所 2 次。要る QC が 1 つ多い）。drift_mode="off" では使わず、
        warnings の `loess_degree_ignored` に書く（matrix_id は渡さないときと同じ）。
    blank_min_ratio: 生体試料の平均がブランク平均のこの倍数未満の特徴を除く（既定 3.0）。
    detection_min_rate: 検出率（gap-fill を除く）の下限（既定 0.5）。detection_by: "group"（群ごと）/ "all"（全体）。
    max_qc_rsd: QC の RSD がこれを超える特徴を除く（既定 0.30）。max_d_ratio: D-ratio の上限（既定 0.5）。
        上の 4 つの足切りは None でその足切りを行わない。QC の RSD の足切りは QC が 2 本未満なら飛ばす。検出率の足切りは gap-fill を
        数えないので、少数の試料にだけある特徴（服薬など）は群ごとでも全体でも落ちうる。
    normalize: "auto"（既定。換算済みか QC が無ければ何もしない。それ以外は PQN）/ "none" / "tic" / "median" / "pqn"。QC が無い
        データで "pqn" を明示すると、参照は生体試料（全群）の特徴ごとの中央値（warnings の `pqn_reference_is_samples`）。
        auto の PQN は QC が 1 本でもあれば当て、参照はその QC の中央値（QC が少ないと参照が揺れる）。TIC・中央値・PQN の係数は
        試料の全特徴から求めるので、一部の試料だけに強い信号（薬剤など）があると、その試料のほかの特徴が見かけ上下がる
        （`normalization_exclude_features` でその信号を係数から外せる）。
        ドリフト補正は、バッチごとに QC が 4 本以上あり、QC が試料の注入の両端を挟むときだけ当たる。
    normalization_features: 正規化の係数をこの特徴 ID だけから計算する（足切りで落ちた特徴も使える。内部標準など）。
    normalization_exclude_features: 係数の計算から外す特徴 ID。
        どちらも係数の計算だけで、外した特徴も行列には残り、係数で割られる。正規化しないとき（"none" と、"none" に決まった "auto"）は
        使わず、warnings の `normalization_features_ignored` に書く（matrix_id は渡さないときと同じ）。空のリストは渡さないのと同じ
        （`normalization_features_empty`）。統合したデータセットのある極性だけ係数の特徴が残らなければ断らず、その極性は残った特徴
        すべてから求めて注意（`normalization_features_not_in_polarity` / `normalization_exclude_features_cover_polarity`）を返す。
    collapse_replicates: 同じ個体（subject_id）の反復注入を 1 個体にまとめる（既定 True）。subject_id が無ければ注入が単位で、
        同じ個体の反復注入は独立な試料として検定に入る（`dataset_overview` が反復らしいまとまりを示す）。
    injection_columns: 個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う属性）。
    group_levels / group_factors: 試料名（と試料の表の群）の因子トークンで群を解く（group_factors は因子軸のリストで優先）。解いた群は
        行列の設定に記録される。無ければ試料の表の群。
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
        raw_config = _config(
            drift_mode=drift_mode, drift_span=drift_span, blank_min_ratio=blank_min_ratio,
            detection_min_rate=detection_min_rate, detection_by=detection_by,
            max_qc_rsd=max_qc_rsd, max_d_ratio=max_d_ratio, normalize=normalize,
            collapse_replicates=collapse_replicates, injection_columns=injection_columns,
            loess_degree=loess_degree, normalization_features=normalization_features,
            normalization_exclude_features=normalization_exclude_features)
        # 設定の誤りは群の解決より先に返す。効かない引数を捨てた注記は前処理が warnings に載せる
        # （前処理は生の設定を正規化し直すので、その注記を拾うために生の設定を渡す）
        config = stats_preprocess.normalize_preprocess_config(raw_config)
        inputs = stats_results.study_inputs(ds, group_levels=group_levels,
                                            group_factors=group_factors)
        study = from_dataset(ds, **inputs)
        prepared = stats_preprocess.preprocess(study, raw_config)
    except DomainError as exc:
        return _domain_envelope(exc)

    try:
        _entry, summary, record = store_matrix(holder, ds, study, prepared, config, inputs)
    except DomainError as exc:
        return _domain_envelope(exc)
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "record": record,
        "next": ("dataset_compare（2 群以上の比較）または dataset_pca を、この matrix_id を"
                 "渡して実行してください（省略しても、このハンドルで最後に作ったこの行列を使う）"),
    }))


# --- dataset_compare ---------------------------------------------------------------------

#: 戻り値・記録の警告で、長い一覧を切るキー（エンジンの注意が持つ一覧）。
_WARNING_LIST_KEYS = ("samples", "sample_ids", "unit_ids", "examples", "subjects")
_TOP_COLUMNS = ["contrast_id", "feature_id", "name", "id_level", "log2fc", "ci_low", "ci_high", "q",
                "verdict"]
#: F 検定の上位（omnibus・要因の効果）の列。
_F_TOP_COLUMNS = ["effect", "feature_id", "name", "id_level", "f", "df1", "df2", "f_q"]
#: top の件数の上限（戻り値を肥大させない）。
_TOP_N_MAX = 50


def _bad_request(message: str, **details) -> DomainError:
    return DomainError("DATASET_BAD_REQUEST", message, details or None)


def _reported_ids(bh_population, reported_feature_ids) -> tuple[list[str] | None, list[str]]:
    """`bh_population="reported"` の母集団（特徴 ID の文字列のリスト。整数の ID も文字列にする）と注記。

    `reported` なのに ID が無い（空・不正）ときだけ DATASET_BAD_REQUEST（補正の母集団が決まらない）。
    `reported` 以外で ID を渡したときは効かない引数の組なので、ID を使わずに注記を返す（誤りにしない）。"""
    if bh_population != "reported":
        if reported_feature_ids is None:
            return None, []
        return None, [f'reported_feature_ids は bh_population="reported" のときだけ補正の母集団に使う。'
                      f'bh_population="{bh_population}" なので使っていない（報告した特徴だけで補正するなら'
                      ' bh_population="reported" を渡して呼び直す）。']
    if (not isinstance(reported_feature_ids, list) or not reported_feature_ids
            or any(isinstance(f, bool) or not isinstance(f, (str, int)) or f == ""
                   for f in reported_feature_ids)):
        raise _bad_request('bh_population="reported" には reported_feature_ids（報告する特徴 ID の'
                           "空でないリスト）が要ります。")
    return list(dict.fromkeys(str(f) for f in reported_feature_ids)), []


def _with_matrix(exc: DomainError, matrix_id, matrix_source) -> DomainError:
    """`matrix_id` を省いた呼び出しのエラーに、使った（作った）行列とその出所を足す。

    既定の前処理で作った行列は呼び直しにそのまま使え、直近の `dataset_preprocess` の行列なら、群がその前処理の
    `group_levels` / `group_factors` で決まっていたことが分かる（既定の群が消えた理由を読めるように）。"""
    if matrix_id is None or matrix_source in (None, "matrix_id"):
        return exc
    details = {**(exc.details or {}), "matrix_id": matrix_id, "matrix_source": matrix_source}
    if matrix_source == "latest_preprocess":
        details["matrix_note"] = ("matrix_id を省いたので、このハンドルで dataset_preprocess が最後に作った行列を"
                                  "使った（群・除外・設定はその前処理のもの）。別の行列は matrix_id で名指しする。")
    return DomainError(exc.code, exc.message, details)


def _output_path(output_path):
    """`output_path` の検証（比較を回す前に済ませる）。無ければ None。"""
    if output_path is None:
        return None
    if not isinstance(output_path, str) or not output_path.strip():
        raise _bad_request("output_path は書き出し先のファイルのパスです。", output_path=output_path)
    path = Path(output_path)
    if path.is_dir():
        raise _bad_request("output_path がフォルダです。ファイルのパスを渡してください。",
                           output_path=output_path)
    return path


def _design_and_comparisons(reference_group, test_group, groups, comparisons, factors):
    """群の指定（ちょうど 1 つ）を spec v2 の `(design.factors, comparisons)` に写す。

    2 群は `reference_group`+`test_group`、3 群以上は `groups`（先頭が基準で各群 − 基準）、
    2 因子などは `comparisons`（spec の形そのまま。`factors` が要る）。それ以外は DATASET_BAD_REQUEST。
    `comparisons` の無い `factors` はここでは見ない（`_factors_without_comparisons` が注記か問いにする）。
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
    return ["group"], [{"id": stats_results.contrast_id_for(t, ref, taken, test_pos=i, ref_pos=0), "factor": "group",
                        "reference": ref, "test": t} for i, t in enumerate(tests, start=1)]


#: needs_input の選択肢 1 つに並べる対比の上限。超える読みは選択肢にせず理由の文に書く（ほかに選択肢が無ければ並べる）。
_MAX_OPTION_CONTRASTS = 50


def _sample_level(sample, factor: str):
    """試料の因子の値（`stats.design._value` と同じ読み方。group / batch は欄、ほかは試料の表の列）。"""
    if factor in ("group", "batch"):
        value = getattr(sample, factor, None)
    else:
        value = (getattr(sample, "attributes", None) or {}).get(factor)
    return None if value is None or str(value) == "" else str(value)


def _id_suffix(level: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", level) or "_"


def _model_samples(study) -> list:
    """比較のモデルに入る試料（生体試料で、除外していないもの）。"""
    return [s for s in study.samples if s.role == "sample" and s.include]


def _named_levels(spec_comparisons) -> set[str]:
    return {c["reference"] for c in spec_comparisons} | {c["test"] for c in spec_comparisons}


def _factor_design_problem(study, factors: list[str]) -> str | None:
    """`factors` のモデル（spec v2 の cell means の設計）が組めない理由。組めるなら None。

    `stats.design.build_design` と同じ規則を試料で先に見る: 比較のモデルに入る試料にその因子の値が無ければ断られ、
    因子が 2 つ以上なら水準の直積のどれかの cell に試料が無いだけで断られる（DESIGN_INVALID）。
    """
    samples = _model_samples(study)
    missing = sorted({f for f in factors for s in samples if _sample_level(s, f) is None})
    if missing:
        return f"因子 {missing} の値が無い試料がある"
    if len(factors) >= 2:
        cells = {tuple(_sample_level(s, f) for f in factors) for s in samples}
        full = set(itertools.product(*(sorted({c[i] for c in cells}) for i in range(len(factors)))))
        empty = sorted(":".join(c) for c in full - cells)
        if empty:
            return f"試料の無い cell がある（{empty[:5]}{' ほか' if len(empty) > 5 else ''}）"
    return None


def _factor_model_comparisons(study, factors: list[str], factor: str, spec_comparisons) -> list[dict] | None:
    """群の指定の対比を「`factors` のモデルで、因子 `factor` の水準として比べる」comparisons に写す。

    群の名前がすべて `factor` の水準でなければ None（その読みは成り立たない）。因子が 2 つ以上なら、ほかの因子の
    水準の組ごとに `at` を付けた対比にする（spec v2 は複数因子の単純な対比に `at` を要る）。モデルが組めること
    （値の欠けも空の cell も無いこと。`_factor_design_problem`）は呼び手が先に確かめる。
    """
    samples = _model_samples(study)
    if not _named_levels(spec_comparisons) <= {_sample_level(s, factor) for s in samples}:
        return None
    others = [f for f in factors if f != factor]
    combos = sorted({tuple(_sample_level(s, f) for f in others) for s in samples})
    out, taken = [], set()
    for c in spec_comparisons:
        for combo in combos:
            if not others:
                item = {**c, "factor": factor}
            else:
                cid = f"{c['id']}_at_{'_'.join(_id_suffix(v) for v in combo)}"
                n = 2
                while cid in taken:
                    cid, n = f"{cid}_{n}", n + 1
                item = {**c, "factor": factor, "id": cid, "at": dict(zip(others, combo))}
            if item["id"] not in taken:
                taken.add(item["id"])
                out.append(item)
    return out or None


def _factors_without_comparisons(study, factors, spec_comparisons,
                                 group_args: list[str]) -> tuple[str | None, str | None]:
    """`comparisons` の無い `factors` を `(注記, needs_input の封筒)` のどちらかにする（渡されていなければ両方 None）。

    群の指定（reference_group+test_group / groups）は `group` の 1 因子のモデルで比べる。`factors` が空（因子を
    挙げていない）か、ちょうど `["group"]`（その設計と同じ）なら効かない（重複した）引数なので注記して進める。
    ほかの因子を挙げたとき、群の名前がどれかの因子の水準なら「`factors` のモデルでその因子の水準として比べる
    （ほかの因子の水準ごとの `at`）」と読め、群の名前が group 列の水準でもあれば「因子を外して group 列で 1 因子の
    モデルで比べる」とも読めて数字が変わるので尋ねる。答えは 1 つの引数に収まらない（前者は `comparisons` を足して
    群の指定 `group_args` を外し、後者は `factors` を外す）ので、選択肢の値は替える引数の dict（`patch` の封筒）。
    因子の水準でないか、`factors` のモデルが組めない（因子の値が無い試料・2 因子以上で試料の無い cell。
    `_factor_design_problem`）なら前者は成り立たないので、group 列で比べて注記する（group 列の水準でもなければ、
    比較そのものが群の誤りを返す）。
    """
    if factors is None:
        return None, None
    if not isinstance(factors, list) or any(not isinstance(f, str) or not f for f in factors):
        raise _bad_request("factors は因子名（文字列）のリストです。", factors=factors)
    if not factors:
        return ("factors=[] は因子を挙げていない（comparisons を渡したときだけ使う引数）ので使わず、群の指定の"
                '設計 design.factors=["group"] で比べた。'), None
    unique = list(dict.fromkeys(factors))
    if unique == ["group"]:
        return ('factors は comparisons を渡したときだけ使う引数で、この比較は design.factors=["group"]（群の指定の'
                "設計）と同じなので、そのまま比べた。"), None
    named = ", ".join(f"{c['reference']} → {c['test']}" for c in spec_comparisons)
    group_reading = _named_levels(spec_comparisons) <= {
        v for s in _model_samples(study) if (v := _sample_level(s, "group")) is not None}
    problem = _factor_design_problem(study, unique)
    readings = [] if problem else [
        (factor, comps) for factor in unique
        if (comps := _factor_model_comparisons(study, unique, factor, spec_comparisons)) is not None]
    if not readings:
        if not group_reading:
            return None, None            # どちらの読みも成り立たない: 比較が群の誤り（使える群の一覧つき）を返す
        why = (f"このモデルは組めない（{problem}）" if problem
               else f"群の名前（{named}）はそのどの因子の水準でもない")
        return (f"factors={unique} は comparisons を渡したときだけ使う引数で、{why}ので、試料の表の群（group 列）の"
                "1 因子のモデルで比べた（ほかの因子で調整するなら covariates に渡す）。"), None
    fits = [(f, c) for f, c in readings if len(c) <= _MAX_OPTION_CONTRASTS]
    shown = fits if (fits or group_reading) else readings
    left_out = [(f, len(c)) for f, c in readings if (f, c) not in shown]
    dropped = " / ".join(group_args)
    options = [{"value": {"comparisons": comps, **{name: None for name in group_args}},
                "label": (f"factors={unique} のモデルで、因子 {factor} の水準として比べる"
                          + (f"（ほかの因子の水準ごとに at を付けた {len(comps)} 対比）" if len(unique) > 1 else "")
                          + f"。comparisons にこの対比を渡し、{dropped} を外す（factors はそのまま）")}
               for factor, comps in shown]
    if group_reading:
        options.append({"value": {"factors": None},
                        "label": "factors を外し、試料の表の群（group 列）の水準として 1 因子のモデルで比べる"
                                 "（ほかの因子で調整するなら covariates に渡す）"})
    reason = (f"factors={unique} は comparisons（対比）を渡したときのモデルの因子です。群の指定（{named}）だけでは、"
              "その因子をモデルに入れてどの水準の組で比べるのかが決まらず、読み方で数字が変わります。選択肢の value の"
              "引数を置き換えて（null の引数は外して）呼び直してください。")
    if left_out:
        reason += "".join(f"因子 {f} の水準として比べる読みは対比が {n} 個になるので選択肢に載せていません"
                          "（その読みなら comparisons を組んで渡す）。" for f, n in left_out)
    return None, needs_input("comparisons", reason, options, patch=True, retry_tools=["dataset_compare"])


def _group_names(spec_comparisons) -> set[str]:
    """対比が group 因子の水準として名指す群の名前（単純な対比の reference / test、交互作用の levels.group、
    at.group）。spec の正規化の前に読むので、形の崩れた対比は飛ばす（誤りは正規化が返す）。"""
    out: set[str] = set()
    for c in spec_comparisons or []:
        if not isinstance(c, dict):
            continue
        if c.get("factor") == "group":
            out.update(v for v in (c.get("reference"), c.get("test")) if isinstance(v, str))
        levels = c.get("levels")
        if isinstance(levels, dict) and isinstance(levels.get("group"), list):
            out.update(v for v in levels["group"] if isinstance(v, str))
        at = c.get("at")
        if isinstance(at, dict) and isinstance(at.get("group"), str):
            out.add(at["group"])
    return out


def _group_levels_under(ds, inputs, group_overrides) -> set[str]:
    """行列の入力 `inputs` の群の上書きを `group_overrides` にしたとき、比較のモデルに入る試料の群の水準。"""
    table = sample_table(ds, exclude_samples=inputs.get("exclude_samples") or (),
                         group_overrides=group_overrides or None)
    return {v for s in table if s.role == "sample" and s.include
            and (v := _sample_level(s, "group")) is not None}


def _matrix_group_conflict(ds, inputs, matrix_id, group_levels, group_factors,
                           spec_comparisons) -> tuple[str | None, str | None]:
    """`matrix_id` と `group_levels` / `group_factors` を一緒に渡したときの `(注記, needs_input の封筒)`。

    行列の群は前処理のときに決まっている。渡した群の引数が、前処理に入る試料をこの行列（`inputs` は行列の入力）と
    同じ群に分けるなら効かない（重複した）引数なので注記して進める——上書きの dict ではなく試料ごとの群で比べる
    （試料の表の群と同じ分け方の `group_levels` は上書きの無い行列と同じ行列になる。`stats_results.groups_match`）。
    違えば「行列の群で比べる（群の引数を外す）」と「渡した群で比べる（`matrix_id` を外す。行列が変わる）」の 2 通りに
    読め、数字が変わるので尋ねる（選択肢の値は替える引数の dict。`patch` の封筒）。選択肢は対比の群の名前がその群の
    水準にある読みだけ: 行列の群でしか比べられなければ群の引数は効かないので注記して進め、どちらでも比べられなければ
    比較が群の誤りを返す（両方 None）。群の引数を渡していなければ両方 None。
    """
    if not (group_levels or group_factors):
        return None, None
    wanted = stats_results.study_inputs(ds, group_levels=group_levels,
                                        group_factors=group_factors)["group_overrides"]
    if stats_results.groups_match(ds, inputs, wanted):
        return ("matrix_id を渡すと群は行列が決める。渡した group_levels / group_factors は試料をこの行列と同じ群に"
                "分けるので、そのまま比べた。"), None
    names = _group_names(spec_comparisons)
    on_matrix = names <= _group_levels_under(ds, inputs, inputs.get("group_overrides"))
    on_passed = names <= _group_levels_under(ds, inputs, wanted)
    if not (on_matrix or on_passed):
        return None, None                # どちらの群でも比べられない: 比較が群の誤り（行列の群の一覧つき）を返す
    if not on_passed:
        return (f"matrix_id を渡すと群は行列が決める。渡した group_levels / group_factors の群には {sorted(names)} が"
                "無く、その群では比べられないので使わず、この行列の群で比べた。"), None
    group_args = {name: None for name, value in (("group_levels", group_levels), ("group_factors", group_factors))
                  if value is not None}
    options = []
    if on_matrix:
        options.append({"value": group_args,
                        "label": f"この行列の群（前処理で決めた群）で比べる（{' / '.join(group_args)} を外す）"})
    options.append({"value": {"matrix_id": None},
                    "label": ("渡した群で比べる（matrix_id を外す。群が同じ直近の dataset_preprocess の行列か、無ければ"
                              "既定の前処理で行列を作る。この行列と同じ前処理の設定で比べるなら、dataset_preprocess に"
                              "その設定と群を渡して作った matrix_id を渡す）")})
    reason = (f"matrix_id={matrix_id} の行列の群と、渡した group_levels / group_factors の群が違います。行列の群は"
              "前処理のときに決まっているので、どちらの群で比べるかで数字が変わります。")
    if not on_matrix:
        reason += (f"この行列の群には {sorted(names)} が無いので、この行列のままでは比べられません（渡した群で比べるには"
                   "行列が変わり、数字はその行列の前処理の設定で決まります）。")
    reason += "選択肢の value の引数を置き換えて（null の引数は外して）呼び直してください。"
    return None, needs_input("matrix_id", reason, options, patch=True, retry_tools=["dataset_compare"])


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


def _fmt(value, spec: str = ".4g") -> str:
    return "" if value is None else format(value, spec)


def _clean(text) -> str:
    return str(text or "").replace("\t", " ").replace("\n", " ")


def _top_table(ds, rows, cids, limit: int) -> tuple[str, dict[str, int]]:
    """対比ごとに changed の上位（q の小さい順）を TSV（列名は 1 回）にする。無ければ空文字。
    名前は `feature_metadata`、無ければ `feature_annotations`（MS1 注釈）から。`id_level` は msms / ms1 / none。
    2 つ目の戻り値は `limit` で載せきれなかった対比 → 行数（載せきれたら空）。"""
    lines, more = [], {}
    for cid in cids:
        changed = [r for r in rows if r["contrast_id"] == cid and r["verdict"] == "changed"
                   and r["q"] is not None]
        changed.sort(key=lambda r: (r["q"], r["feature_id"]))
        if len(changed) > limit:
            more[cid] = len(changed) - limit
        for r in changed[:limit]:
            ident = compare_table.identity_of(ds, r["feature_id"])
            lines.append("\t".join([cid, str(r["feature_id"]), _clean(ident["name"]), ident["id_level"],
                                    _fmt(r["log2fc"]), _fmt(r["ci_low"]), _fmt(r["ci_high"]),
                                    _fmt(r["q"], ".3g"), r["verdict"]]))
    return ("\n".join(["\t".join(_TOP_COLUMNS), *lines]) if lines else ""), more


def _f_top_table(ds, engine: dict, limit: int) -> tuple[str, dict[str, int]]:
    """omnibus F と要因の F の効果ごとに、q が閾値以下の上位（q の小さい順）を TSV にする。無ければ空文字。
    effect は omnibus なら `omnibus`、要因なら効果の名前（`group` / `group:B` など）。
    2 つ目の戻り値は `limit` で載せきれなかった効果 → 行数（載せきれたら空）。"""
    q_thr = engine["verdict_rule"]["q_threshold"]
    blocks = []
    if engine.get("omnibus"):
        blocks.append(("omnibus", engine["omnibus"]))
    factorial = engine.get("factorial") or {}
    for name in factorial.get("effects") or {}:
        blocks.append((name, [r for r in factorial["rows"] if r["effect"] == name]))
    lines, more = [], {}
    for effect, rows in blocks:
        hits = [r for r in rows if r["f_q"] is not None and r["f_q"] <= q_thr]
        hits.sort(key=lambda r: (r["f_q"], r["feature_id"]))
        if len(hits) > limit:
            more[effect] = len(hits) - limit
        for r in hits[:limit]:
            ident = compare_table.identity_of(ds, r["feature_id"])
            lines.append("\t".join([effect, str(r["feature_id"]), _clean(ident["name"]), ident["id_level"],
                                    _fmt(r["f"]), _fmt(r["df1"]), _fmt(r["df2"]), _fmt(r["f_q"], ".3g")]))
    return ("\n".join(["\t".join(_F_TOP_COLUMNS), *lines]) if lines else ""), more


def _contrast_summaries(engine: dict, spec: dict, log2fc_threshold: float) -> list[dict]:
    """対比ごとの要約。`large_changes` は判定が changed で |log2fc| が閾値以上の件数。"""
    q_thr = engine["verdict_rule"]["q_threshold"]
    presence = engine["presence"]
    by_contrast = (engine.get("confounding") or {}).get("by_contrast") or {}
    out = []
    for comp in spec["comparisons"]:
        cid = comp["id"]
        conf = by_contrast.get(cid) or {}
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
            "confounded": conf.get("confounded"),          # true / false / null（評価できない）
            "batch_status": conf.get("batch_status"),
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


def _matrix_not_found(matrix_id) -> str:
    return missing_state(
        "stat_matrix", ["dataset_preprocess"],
        f"行列 {matrix_id!r} がありません（無い・または除外や試料の表の変更で古くなりました）。"
        "dataset_preprocess を実行し直して、その matrix_id を渡してください。")


def _write_compare_record(holder, entry, record, spec, summary, matrix_source: str) -> dict:
    """比較の記録を書く（親は行列の前処理の記録）。書けないときは `{"written": False, "reason"}`。"""
    try:
        reports_dir = mcp_core._resolve_report_dir(holder)
        parents = [entry["record_id"]] if entry.get("record_id") else []
        written = stats_records.write_stats_record(
            reports_dir, "compare",
            inputs={"dataset_hash": entry["dataset_hash"], "matrix_id": record["matrix_id"],
                    "result_id": record["result_id"]},
            parent_ids=parents,
            settings={"spec": spec, "log2fc_threshold": record["log2fc_threshold"],
                      "matrix_source": matrix_source},
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
    reported_feature_ids: list[str] | None = None,
    p_adjust: str = "bh",
    q_threshold: float = 0.05,
    equivalence_fc: float = 1.2,
    log2fc_threshold: float = 1.0,
    allow_confounded: bool = True,
    group_levels: list[str] | None = None,
    group_factors: list[list[str]] | None = None,
    top_n: int = 10,
    output_path: str | None = None,
    dataset_handle: str | None = None,
) -> str:
    """2 群・多群・2 因子の比較を統計エンジン（analysis-spec v2）で行い、要約を返す。

    群の指定はちょうど 1 つ: reference_group+test_group（2 群）/ groups（3 群以上。先頭が基準で各群 − 基準の対比と omnibus F）/
    comparisons（spec v2 の comparisons の形。2 因子・`at`・交互作用。factors が要る。因子が 2 つ以上で cell がそろえば主効果・
    交互作用の F も `factorial` に返す）。comparisons の無い factors は、空（[]）か ["group"] なら notes に書いてそのまま比べる。
    ほかの因子で、群の名前がその因子の水準でそのモデルが組めるなら数字が変わるので needs_input（patch の封筒。選択肢の value は
    置き換える引数の dict）、組めなければ group 列で比べて notes に書く。
    判定（changed / unchanged / inconclusive）は q と CI で決まり、log2fc_threshold は判定に使わない（large_changes の数え方・図の
    色分け・書き出しの significant 列だけ）。全行は result_id でセッションに保持し、直近の結果（図・書き出しの既定）にもなる。
    `reports/stats/` に記録を 1 件追記する（書けないときは record.written=false。比較は成功）。戻り値の先頭は `confounding`（群と
    バッチの交絡。対比ごとの `by_contrast`）。戻り値と needs_input の選択肢は `docs_read(topic="dataset")` の 11.2。

    covariates: 試料の表の列名。値がすべて数値なら numeric、さもなくば categorical。
    block: "subject_id" だけ（個体の固定効果。同じ個体が複数の群・cell に出る対応のある計画で個体差を除く）。既定のまとめた行列
        （dataset_preprocess の collapse_replicates=True）でも効く。個体が群に入れ子（どの個体も 1 つの群にしか出ない）なら計画行列が
        階数不足になり DESIGN_RANK_DEFICIENT（入れ子の技術反復は block ではなく既定のまとめで扱う）。
    matrix_id: dataset_preprocess が返した行列 ID（その前処理・群・除外のまま比較する）。省略時は、このハンドルで dataset_preprocess が
        最後に作った行列（group_levels / group_factors を渡したときは群が同じものだけ）、無ければ既定の前処理で行列を作って登録する
        （drift=auto が補正を当てなければ警告 drift_not_applied_by_default。補正するなら dataset_preprocess(drift_mode="on") の ID を
        渡す）。使った行列は戻り値の matrix_id と matrix_source（"matrix_id" / "latest_preprocess" / "default"。省いた呼び出しの
        エラーにも details.matrix_id / matrix_source）。古い・無い ID は missing_state。ID と一緒に渡した group_levels /
        group_factors は、試料を行列と同じ群に分けるなら notes に書いてそのまま比べ、違えば needs_input（patch の封筒。群の引数を
        外すか matrix_id を外すか）。
    method: "moderated_t"（既定）など（エンジンの検定）。welch / student / mannwhitney は 1 因子・共変量なし・ブロックなしだけで、
        2 因子の主効果・交互作用の F は moderated_t のときだけ出る。effect_size: "log2_mean_difference"（既定）/
        "arithmetic_mean_ratio"。log2fc は test − reference（正なら test が高い）。
        割れたアライメントスポット（`dataset_alignment_qc`）はまとめずに比べるので、q はまとめる前の全スポットの母集団の値。
        まとめた行列で比べるなら `dataset_analysis_run`（`features.merge_split` 既定 true）。
    bh_population: 補正の母集団 "matrix"（既定）/ "annotated" / "reported"（reported_feature_ids の特徴だけ。ID が無ければ
        DATASET_BAD_REQUEST。reported 以外で渡した ID は使わず notes に書く）。
    p_adjust: 多重比較の補正 "bh"（既定。FDR）/ "bonferroni" / "holm"（FWER）/ "none"。q 列はこの方法の値で、戻り値と記録に
        p_adjust_method として残る。dataset_export_differential・dataset_export_evidence_bundle の q_value は方法によらず BH の q。
    q_threshold / equivalence_fc: 判定の q の閾値（既定 0.05）と同等とみなす倍率（既定 1.2）。
    allow_confounded: 既定 true（群とバッチが完全交絡でも比較し、confounding と対比の confounded に印を付ける）。false なら完全交絡
        （全体か、どれかの対比）で STATISTIC_CONFOUNDED にして止める。
    top_n: 対比ごとに changed の上位（q の小さい順）を top に載せる件数（既定 10。F 検定があれば効果ごとの上位を f_top に）。50 を
        超える値は断らず 50 に下げて notes に書き、50 行で切った対比・効果があれば top_truncated: true と top_n_more / f_top_n_more。
        整数でない・負の値は DATASET_BAD_REQUEST。
    output_path: 全特徴 × 全対比の全列（統計量・df・CI・SE・群ごとの n と検出数・同定の列）を TSV で書く先。F 検定があれば
        `<stem>.ftests<拡張子>` にも書く。dataset_export_differential の 15 列とは別物。書けなくても比較は成功し output.written=false。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    used_matrix = matrix_source = None
    try:
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 0:
            raise _bad_request("top_n は 0 以上の整数です。", top_n=top_n)
        reported, notes = _reported_ids(bh_population, reported_feature_ids)
        top_capped = top_n > _TOP_N_MAX
        if top_capped:
            # 長いだけの指定は断らない: 上限に下げて比べ、全行の取り出し方を書く
            notes.append(f"top_n={top_n} は上限の {_TOP_N_MAX} に下げた（戻り値を肥大させない）。全行は result_id・"
                         "output_path（全列の TSV）・dataset_export_differential で取り出せる。")
            top_n = _TOP_N_MAX
        if p_adjust not in P_ADJUST_METHODS:
            raise _bad_request(f"p_adjust は {list(P_ADJUST_METHODS)} のどれかです。", p_adjust=p_adjust)
        out_path = _output_path(output_path)
        if (isinstance(log2fc_threshold, bool) or not isinstance(log2fc_threshold, (int, float))
                or not math.isfinite(log2fc_threshold) or log2fc_threshold < 0):
            raise _bad_request("log2fc_threshold は 0 以上の有限の数です。",
                               log2fc_threshold=log2fc_threshold)
        if block is not None and block != "subject_id":
            raise _bad_request('block は "subject_id" だけを受けます。', block=block)
        design_factors, spec_comparisons = _design_and_comparisons(
            reference_group, test_group, groups, comparisons, factors)

        try:
            entry, inputs, preprocess = stats_results.matrix_inputs(
                ds, matrix_id, group_levels=group_levels, group_factors=group_factors,
                prefer_latest=True)
        except stats_results.MatrixNotFound:
            return _matrix_not_found(matrix_id)
        if matrix_id is not None:
            note, asked = _matrix_group_conflict(ds, inputs, matrix_id, group_levels, group_factors,
                                                 spec_comparisons)
            if asked is not None:
                return asked
            if note is not None:
                notes.append(note)
        matrix_source = ("matrix_id" if matrix_id is not None
                         else "latest_preprocess" if entry is not None else "default")
        used_matrix = entry["matrix_id"] if entry is not None else None
        study = from_dataset(ds, **inputs)
        if comparisons is None:
            group_args = (["groups"] if groups is not None else ["reference_group", "test_group"])
            note, asked = _factors_without_comparisons(study, factors, spec_comparisons, group_args)
            if asked is not None:
                return asked
            if note is not None:
                notes.append(note)
        spec = spec_v2.normalize_spec_v2(stats_results.compare_spec(
            factors=design_factors, comparisons=spec_comparisons,
            covariates=_covariate_specs(study, covariates), block=block, method=method,
            effect_size=effect_size, bh_population=bh_population, q_threshold=q_threshold,
            equivalence_fc=equivalence_fc, allow_confounded=allow_confounded,
            preprocess=preprocess, inputs=inputs, p_adjust=p_adjust, reported_feature_ids=reported))

        prepared = (entry["prepared"] if entry is not None
                    else stats_results.prepare_default(study, spec["preprocess"]))

        preprocess_record = None
        if entry is None:
            entry, preprocess_record = register_default(holder, ds, study, prepared,
                                                         preprocess, inputs)
            used_matrix = prepared.matrix_id      # 行列は登録済み。呼び直しにそのまま使える
        run = spec_v2.run_spec_v2(spec, study, prepared=prepared)
        try:
            record = stats_results.register_compare(ds, run, log2fc_threshold=log2fc_threshold,
                                                    units=prepared.units)
        except ValueError as exc:
            raise DomainError("STUDY_DATA_INVALID", str(exc)) from exc
    except DomainError as exc:
        return _domain_envelope(_with_matrix(exc, used_matrix, matrix_source))

    engine = record["engine"]
    conf = engine["confounding"]
    q_thr = engine["verdict_rule"]["q_threshold"]
    summary = {
        # 交絡は先頭に置く（既定では止めないので、読み落とさせない）
        "confounding": {"batch_status": conf.get("batch_status"),
                        "batch_table": conf.get("batch_table"),
                        "by_contrast": conf.get("by_contrast"),
                        "warnings": conf.get("warnings")},
        "result_id": record["result_id"],
        "matrix_id": record["matrix_id"],
        "matrix_source": matrix_source,
        "rule_hash": record["rule_hash"],
        "test": engine["test"],
        "p_adjust_method": engine["p_adjust_method"],
        "effect_size_definition": {"name": engine["effect_size_definition"]["name"]},
        "ci_level": engine["ci_level"],
        "verdict_rule": engine["verdict_rule"],
        "bh_population": engine["bh_population"],
        "subject_source": engine["subject_source"],
        "group_n": engine["group_n"],
        "prior": engine["prior"],
        "contrasts": _contrast_summaries(engine, spec, record["log2fc_threshold"]),
        "omnibus": _omnibus_summary(engine),
        "factorial": factorial_summary(engine, q_thr),
        "warnings": engine["warnings"],
    }
    summary = truncate_lists(json_safe(round_floats(summary)), _STEP_LIST_LIMIT,
                             keys=_WARNING_LIST_KEYS)
    output = None
    if out_path is not None:
        try:
            output = {"written": True, **compare_table.write_compare_tables(out_path, ds, record)}
        except OSError as exc:
            output = {"written": False, "path": str(out_path), "reason": f"OSError: {exc}"}
    written = _write_compare_record(holder, entry, record, spec, summary, matrix_source)
    if preprocess_record is not None:
        written = {**written, "preprocess_record": preprocess_record}
    cids = [c["id"] for c in spec["comparisons"]]
    top, top_more = _top_table(ds, engine["rows"], cids, top_n)
    f_top, f_top_more = _f_top_table(ds, engine, top_n)
    # 上限に下げた呼び出しで切った対比・効果があれば、切った印と行数（下げていなければ頼んだ件数どおり）
    cut = {} if not top_capped else {
        **({"top_truncated": True} if top_more or f_top_more else {}),
        **({"top_n_more": top_more} if top_more else {}),
        **({"f_top_n_more": f_top_more} if f_top_more else {})}
    return json_payload(json_safe({
        "status": "success",
        **summary,
        "top": top,
        **({"f_top": f_top} if f_top else {}),
        **cut,
        **({"output": output} if output is not None else {}),
        **({"notes": notes} if notes else {}),
        "record": written,
        "next": ("全行は result_id で保持した（直近の結果として書き出し・図の既定になる）。対比が複数の結果は"
                 "dataset_export_differential・dataset_export_evidence_bundle・dataset_plot(kind=volcano)・"
                 "dataset_feature_search には contrast_id を渡して対比を選ぶ（交互作用の対比は 2 群の形にならない）"),
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

    軸は生体の個体だけで決め、QC は同じ軸へ投影する。欠測は PCA のためだけに補う（行列は変えない）。値は log2 のまま。スコア・
    ローディングの全量は戻り値に載せず result_id でセッションに保持し、直近の PCA（図の既定）にもなる。`reports/stats/` に記録を
    1 件追記する（書けないときは record.written=false。PCA は成功）。

    matrix_id: dataset_preprocess が返した行列 ID。省略時は、このハンドルで dataset_preprocess が最後に作った行列、無ければ既定の
        前処理で行列を作って登録する（drift=auto が補正を当てなかったら警告 drift_not_applied_by_default。補正するなら
        dataset_preprocess(drift_mode="on") の ID を渡す）。使った行列は戻り値の matrix_id と matrix_source（"matrix_id" /
        "latest_preprocess" / "default"。省いた呼び出しのエラーにも details.matrix_id / matrix_source）。古い・無い ID は missing_state。
    scaling: "autoscale"（既定。平均 0・SD 1）/ "pareto"（SD の平方根で割る）。
    n_components: 主成分数（既定 5。個体数・特徴数の小さい方が上限）。
    impute: 欠測の補い方 "half_min"（既定）/ "knn" / "none"（欠測があれば PCA_MISSING_VALUES）。
    group_levels / group_factors: 色分けの表示と結果の記録にだけ使う群（試料名の因子トークン）。PCA の数字・行列の選択には効かない
        （matrix_id と併用できる）。無ければ行列の個体の群。
    dataset_handle: `load_dataset` が返したハンドル。省略時は対話のセッションのデータセット。
    """
    holder = resolve_holder(dataset_handle)
    if isinstance(holder, str):
        return holder
    ds = holder.dataset
    if ds is None:
        return missing_state("dataset", ["load_dataset"],
                             "DatasetState がありません。先に load_dataset を実行してください。")
    used_matrix = matrix_source = None
    preprocess_record = None
    try:
        # 群の指定の検証は行列を登録する前に済ませる（不正な指定が行列・記録を残さない）
        label_overrides = stats_results.study_inputs(
            ds, group_levels=group_levels, group_factors=group_factors)["group_overrides"]
        try:
            # 群は行列を選ばない（色分けだけ）。省略時は直近の dataset_preprocess の行列
            entry, inputs, config = stats_results.matrix_inputs(ds, matrix_id, prefer_latest=True)
        except stats_results.MatrixNotFound:
            return _matrix_not_found(matrix_id)
        matrix_source = ("matrix_id" if matrix_id is not None
                         else "latest_preprocess" if entry is not None else "default")
        used_matrix = entry["matrix_id"] if entry is not None else None
        if entry is not None:
            prepared = entry["prepared"]
        else:
            study = from_dataset(ds, **inputs)
            prepared = stats_results.prepare_default(study, config)
            entry, preprocess_record = register_default(holder, ds, study, prepared, config, inputs)
            used_matrix = prepared.matrix_id      # 行列は登録済み。呼び直しにそのまま使える
        result = stats_pca.pca(prepared, scaling=scaling, n_components=n_components, impute=impute)
        labels = stats_results.unit_group_labels(ds, prepared, label_overrides)
        record = stats_results.register_pca(ds, result, labels)
    except DomainError as exc:
        return _domain_envelope(_with_matrix(exc, used_matrix, matrix_source))

    groups: dict[str, int] = {}
    for label in labels.values():
        groups[label] = groups.get(label, 0) + 1
    summary = truncate_lists(json_safe(round_floats({
        "result_id": record["result_id"],
        "matrix_id": record["matrix_id"],
        "matrix_source": matrix_source,
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
         "group_labels": label_overrides, "matrix_source": matrix_source}, summary)
    if preprocess_record is not None:
        written = {**written, "preprocess_record": preprocess_record}
    # drift=auto が補正を当てなかった行列なら、その理由と呼び直し方を添える（数字は補正なしで出ている）
    held = [w for w in prepared.warnings if w.get("code") == "drift_not_applied_by_default"]
    return json_payload(json_safe({
        "status": "success",
        **summary,
        **({"warnings": held} if held else {}),
        "record": written,
        "loadings_note": ("スコアとローディングの全量（特徴数 × 主成分数）は本要約に非同梱。result_id で"
                          "セッションに保持しており、図は dataset_plot(kind=\"pca\" / \"loadings\") で出す。"
                          "QC は軸の計算に使わず投影だけ（n_qc_projected）。"),
    }))
