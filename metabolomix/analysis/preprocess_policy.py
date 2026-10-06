"""conservative-v1: 自動前処理レシピの解決と適用状況の監査（spec §8）。

自動パイプラインは「証拠が無いのに補正した」と「補正したと言いながら実は
補正していない」の両方をやってはいけない。ここが持つのは2つの役目:

``resolve_policy``
    include=true のサンプルメタデータ（``confidence`` を含む）とQC健全性
    （既存 :func:`metabolomix.analysis.preprocessing.detect_failed_qc` の生強度評価）
    から、conservative-v1 が採用してよいレシピを決める。**未確認の batch /
    injection_order / qc_pool は auto のドリフト補正・正規化・QC-RSDフィルタを
    一切有効化しない**（値があるだけでは足りない。confidence を見る）。
    明示要求（"auto" 以外の値）は最後に検査して反映するが、前提を欠く明示要求は
    ``PREPROCESS_PREREQUISITE_MISSING`` にする——auto での同じ状況は
    ``skipped_steps`` への記録に留める。

``check_applied_policy``
    resolve_policy の判断が、実際に数値層（:func:`preprocessing.preprocess`）を
    走らせた結果（report）と一致しているかを再照合する。数値層の caveat だけに
    委ねると、「要求した」と「実施できた」の区別が呼び出し側から見えなくなる
    （arf 側で既に踏んだ罠、CLAUDE.md 参照）。

数値式の変更は一切行わない。ここは「既存実装を条件付きで組み合わせる」だけの層。
"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis import preprocessing
from metabolomix.core.atomic_io import DomainError

__all__ = ["EXPLICIT_POLICY_VERSION", "POLICY_VERSION", "check_applied_policy",
           "explicit_policy", "resolve_policy"]

#: 数値閾値を変えるときはここを上げる（spec §8.1: 「数値閾値の変更はpolicy版を上げる」）。
POLICY_VERSION = "conservative-v1"

#: conservative-v1 が解決する7フィールドの既定値（"auto" は本モジュールが判断する）。
_PP_DEFAULTS = {
    "policy": POLICY_VERSION,
    "normalize": "auto",
    "blank_min_fold": "auto",
    "drift_correct": "auto",
    "max_qc_rsd": "auto",
    "impute": "half_min",
    "min_detection_rate": 0.0,
}

#: resolved_recipe のキー → preprocessing.preprocess() が報告するステップ名。
_STEP_NAMES = {
    "normalize": "normalize",
    "blank_min_fold": "blank_filter",
    "drift_correct": "drift_correct",
    "max_qc_rsd": "qc_rsd_filter",
}

#: 正規化・ドリフト補正・QC-RSDフィルタの「健全なQCが何件以上必要か」（spec §8.1 の表）。
_MIN_QC_FOR_NORMALIZE = 3
_MIN_QC_FOR_DRIFT = 4
_MIN_QC_FOR_RSD = 3
_MAX_QC_RSD_DEFAULT = 0.30


#: v2（LC–MSメタボロミクス profile）が使う policy 名。conservative-v1 とは
#: 別の名前にする——v2 では profile が全ステップを明示するので、ここが
#: 「自動で決めた」判断を1つも持たないことを名前で言い切る。
EXPLICIT_POLICY_VERSION = "profile-explicit-v2"


def explicit_policy(recipe: dict) -> dict:
    """profile の `matrix_recipes` 1件を、そのまま plan の形にして返す。

    conservative-v1 は「証拠が無ければ自動で無効化する」層で、未確認の
    batch/injection_order/qc_pool を見て勝手にステップを落とす。v2 の profile は
    normalize・drift_correct・filter・impute を**明示**しており、そこへ自動判断を
    重ねると「profile にそう書いてあるのに実行されていない」状態が黙って生まれる
    （前提不足は `matrix_state` が `QC_PREREQUISITE_MISSING` で止める役目）。

    だから v2 はこの関数を通す。推測は1つも足さず、`skipped_steps` は常に空で、
    plan の形だけを conservative-v1 と揃える（`check_applied_policy` を共用する
    ため）。
    """
    resolved = dict(recipe)
    resolved["policy"] = EXPLICIT_POLICY_VERSION
    return {
        "requested_recipe": dict(recipe),
        "resolved_recipe": resolved,
        "applied_steps": [],
        "skipped_steps": [],
        "reasons": {},
        "policy_version": EXPLICIT_POLICY_VERSION,
        "assumptions": [],
    }


def _prerequisite_missing(step: str, message: str, **details) -> None:
    raise DomainError("PREPROCESS_PREREQUISITE_MISSING", message,
                      {"step": step, **details})


def _row_field(row: dict, field: str) -> tuple[object, object]:
    """行のvalueとprovenance.confidenceを返す（provenance欠落は None 扱い）。"""
    value = row.get(field)
    prov = (row.get("provenance") or {}).get(field) or {}
    return value, prov.get("confidence")


def _merge_requested(requested: dict | None) -> dict:
    merged = dict(_PP_DEFAULTS)
    merged.update(requested or {})
    return merged


def _gate_reason(*, batch_ok: bool, same_pool: bool, pool_ambiguous: bool,
                 n_qc: int, min_qc: int, failed_qc: list[str],
                 need_order: bool = False, order_confirmed_all: bool = True,
                 need_enclosure: bool = False, all_samples_enclosed: bool = True) -> str:
    """auto でこの手順を見送った理由を日本語で組み立てる。"""
    parts: list[str] = []
    if not batch_ok:
        parts.append("include対象のbatchが単一の確認済み値ではありません")
    if pool_ambiguous:
        parts.append("QCのqc_poolが複数の確認済み値に分かれています")
    elif not same_pool:
        parts.append("QCのqc_poolが単一の確認済み値ではありません")
    if failed_qc:
        parts.append(f"失敗疑いのQCがあります: {', '.join(failed_qc)}")
    if n_qc < min_qc:
        parts.append(f"健全なQCが{min_qc}件未満です(n_qc={n_qc})")
    if need_order and not order_confirmed_all:
        parts.append("include対象の注入順が全件確認済みではありません")
    if need_enclosure and not all_samples_enclosed:
        parts.append("QCの注入順区間に入らないsampleがあります")
    return "、".join(parts) if parts else "条件不足のため見送りました"


def resolve_policy(ds, requested: dict, metadata: list[dict]) -> dict:
    """conservative-v1でrequestedを解決し、根拠付きのplanを返す。

    `metadata` は ``ds.sample_names`` と同じ順・同じ件数のsample-manifest.v1行
    （Task 9 ``resolve_metadata``/``metadata_rows`` と同じ契約）。この関数は
    ``ds`` を一切書き換えない（判断材料として ``ds.feature_matrix`` を読むだけ）。
    """
    requested_recipe = _merge_requested(requested)
    policy = requested_recipe.get("policy", POLICY_VERSION)
    if policy != POLICY_VERSION:
        raise DomainError(
            "PREPROCESS_POLICY_UNKNOWN",
            f"未知のpolicyです: {policy!r}（対応済みは{POLICY_VERSION!r}のみ）。",
            {"policy": policy},
        )

    sample_names = list(getattr(ds, "sample_names", []) or [])
    if len(metadata) != len(sample_names):
        raise DomainError(
            "SAMPLE_MANIFEST_INVALID",
            f"metadataの行数がサンプル数と一致しません(rows={len(metadata)}, "
            f"samples={len(sample_names)})。",
            {"rows": len(metadata), "samples": len(sample_names)},
        )
    if getattr(ds, "feature_matrix", None) is None:
        _prerequisite_missing(
            "resolve_policy", "DatasetStateに定量行列がありません。dataset_loadを先に実行してください。")

    # ---------- include=true 対象の集計（constraint: 全集計はinclude=trueのみ） ----------
    included = [(name, row) for name, row in zip(sample_names, metadata)
                if row.get("include", True)]
    n_qc = sum(1 for _, row in included if row.get("role") == "qc")
    n_blank = sum(1 for _, row in included if row.get("role") == "blank")
    n_sample = sum(1 for _, row in included if row.get("role") == "sample")

    # batch: include対象が単一の確認済みbatchに属するか。
    batch_values: set = set()
    batch_ok = True
    for _, row in included:
        value, confidence = _row_field(row, "batch")
        if confidence != "confirmed" or value is None:
            batch_ok = False
        batch_values.add(value)
    if len(batch_values) > 1:
        batch_ok = False

    # pool: QC行のqc_poolが単一の確認済み値か。
    qc_pool_values: set = set()
    pool_missing = False
    for _, row in included:
        if row.get("role") != "qc":
            continue
        value, confidence = _row_field(row, "qc_pool")
        if confidence != "confirmed" or value is None:
            pool_missing = True
        else:
            qc_pool_values.add(value)
    pool_ambiguous = len(qc_pool_values) > 1
    same_pool = (not pool_missing) and len(qc_pool_values) == 1

    # QC健全性: 既存detect_failed_qcの生強度評価（include対象のみに絞った部分行列で判定）。
    matrix = np.asarray(ds.feature_matrix, dtype=float).T  # samples x features
    name_index = {name: i for i, name in enumerate(sample_names)}
    included_idx = [name_index[name] for name, _ in included]
    sub_matrix = matrix[included_idx] if included_idx else matrix[:0]
    sub_names = [name for name, _ in included]
    sub_roles = {name: (row.get("role") or "sample") for name, row in included}
    qc_health = preprocessing.detect_failed_qc(sub_matrix, sub_roles, sub_names)
    failed_qc = list(qc_health.get("failed", []))

    qc_eligible = batch_ok and same_pool and not failed_qc

    # 注入順: 解析対象(include=true)の全件が確認済みか。全sampleがQC区間内か。
    order_confirmed_all = True
    orders_by_name: dict[str, object] = {}
    for name, row in included:
        value, confidence = _row_field(row, "injection_order")
        orders_by_name[name] = value
        if confidence != "confirmed" or value is None:
            order_confirmed_all = False

    qc_orders = [orders_by_name[name] for name, row in included
                if row.get("role") == "qc" and orders_by_name.get(name) is not None]
    sample_orders = [orders_by_name[name] for name, row in included
                     if row.get("role") == "sample" and orders_by_name.get(name) is not None]
    if sample_orders and qc_orders:
        lo, hi = min(qc_orders), max(qc_orders)
        all_samples_enclosed = all(lo <= o <= hi for o in sample_orders)
    else:
        # sampleが無い（比較対象が無い）なら「区間外」も起こりようがない。
        all_samples_enclosed = not sample_orders

    reasons: dict[str, str] = {}
    skipped_steps: list[str] = []

    def _skip(step: str, reason: str) -> None:
        reasons[step] = reason
        skipped_steps.append(step)

    resolved = dict(requested_recipe)
    assumptions_normalize_weak: dict | None = None

    # ---------- normalize ----------
    normalize_req = requested_recipe["normalize"]
    normalize_auto_ok = qc_eligible and n_qc >= _MIN_QC_FOR_NORMALIZE
    if normalize_req == "auto":
        if normalize_auto_ok:
            resolved["normalize"] = "pqn"
        else:
            resolved["normalize"] = "none"
            _skip("normalize", _gate_reason(
                batch_ok=batch_ok, same_pool=same_pool, pool_ambiguous=pool_ambiguous,
                n_qc=n_qc, min_qc=_MIN_QC_FOR_NORMALIZE, failed_qc=failed_qc))
    elif normalize_req == "pqn":
        # 明示PQNは可能——拒否するのは次の2つだけ:
        #   (1) 複数poolを単一QC参照にする指定（spec §8.1「異poolを単一QC参照に
        #       する要求は拒否する」）。
        #   (2) 参照に使えるQCがinclude対象に1件も無い（preprocessing.normalize()
        #       の pqn 分岐は QC 行が0件だと全サンプル中央値へ無言でフォール
        #       バックする——「PQN(QC参照)を要求したのに実際は全サンプル参照で
        #       計算した」という乖離を防ぐため、これは実施不能として止める。
        #       spec §8.2「明示的に要求した処理を実施できなければneeds_input」）。
        # batch未確認・QC不健全・QC件数がauto閾値未満は拒否しない
        # （spec §8.1本文は明示指定の拒否条件をpool曖昧性だけに絞っている。
        # brief step3の「QC不足の明示補正はPREPREQUISITE_MISSING」という一般則は
        # drift_correct/max_qc_rsd等の「補正」系に対する記述で、normalizeには
        # 及ばないとcontroller裁定R14で確認済み——docs/task-11-report.md参照）。
        if pool_ambiguous:
            _prerequisite_missing(
                "normalize",
                "複数のQC poolを単一のQC参照として扱うPQN正規化は要求できません。",
                qc_pools=sorted(qc_pool_values))
        if n_qc == 0:
            _prerequisite_missing(
                "normalize",
                "PQN正規化の参照に使える健全なQCがinclude対象に1件もありません。",
                n_qc=n_qc)
        resolved["normalize"] = "pqn"
        # auto方針が要求する健全QC件数(3件以上)未満、または失敗疑いのQCを含んだ
        # ままの明示PQNは、実施はするが「QC補正済み」という誤った含意を避けるため
        # 対象を記録する（spec §8.1「失敗疑いのQCがある場合…対象sample_idを報告
        # する」、spec §8.2「『QC補正済み』等の誤った表現をしない」）。
        if n_qc < _MIN_QC_FOR_NORMALIZE or failed_qc:
            qc_sample_ids = [name for name, row in included if row.get("role") == "qc"]
            weak_parts = []
            if n_qc < _MIN_QC_FOR_NORMALIZE:
                weak_parts.append(f"健全なQCが{_MIN_QC_FOR_NORMALIZE}件未満です(n_qc={n_qc})")
            if failed_qc:
                weak_parts.append(f"失敗疑いのQCを含んだまま続行します: {', '.join(failed_qc)}")
            reasons["normalize"] = (
                "明示PQN正規化は条件不足のまま続行しました（"
                + "、".join(weak_parts)
                + "）。QC補正済みとはみなせません。"
            )
            assumptions_normalize_weak = {
                "n_qc": n_qc,
                "min_qc_for_auto": _MIN_QC_FOR_NORMALIZE,
                "failed_qc": list(failed_qc),
                "sample_ids": qc_sample_ids,
            }
    else:
        # tic/median/none の明示はQC前提を問わずそのまま反映する。
        resolved["normalize"] = normalize_req

    # ---------- blank_min_fold ----------
    blank_req = requested_recipe["blank_min_fold"]
    blank_prereq_ok = n_blank >= 1 and n_sample >= 1
    if blank_req == "auto":
        if blank_prereq_ok:
            resolved["blank_min_fold"] = 3.0
        else:
            resolved["blank_min_fold"] = None
            _skip("blank_filter",
                 f"ブランクまたは生体試料がinclude対象にありません(n_blank={n_blank}, "
                 f"n_sample={n_sample})。")
    elif blank_req is None:
        resolved["blank_min_fold"] = None
    else:
        if not blank_prereq_ok:
            _prerequisite_missing(
                "blank_filter",
                "ブランク背景除去にはブランクと生体試料の両方がinclude対象に必要です。",
                n_blank=n_blank, n_sample=n_sample)
        resolved["blank_min_fold"] = blank_req

    # ---------- drift_correct ----------
    drift_req = requested_recipe["drift_correct"]
    drift_prereq_ok = (qc_eligible and n_qc >= _MIN_QC_FOR_DRIFT
                      and order_confirmed_all and all_samples_enclosed)
    if drift_req == "auto":
        resolved["drift_correct"] = bool(drift_prereq_ok)
        if not drift_prereq_ok:
            _skip("drift_correct", _gate_reason(
                batch_ok=batch_ok, same_pool=same_pool, pool_ambiguous=pool_ambiguous,
                n_qc=n_qc, min_qc=_MIN_QC_FOR_DRIFT, failed_qc=failed_qc,
                need_order=True, order_confirmed_all=order_confirmed_all,
                need_enclosure=True, all_samples_enclosed=all_samples_enclosed))
    elif drift_req is False:
        resolved["drift_correct"] = False
    else:  # True
        if not drift_prereq_ok:
            _prerequisite_missing(
                "drift_correct",
                "ドリフト補正の前提（単一の確認済みbatch・単一の確認済みpool・健全なQC"
                f"{_MIN_QC_FOR_DRIFT}件以上・注入順が全件確認済み・全sampleがQC区間内）"
                "を満たしません。",
                n_qc=n_qc, batch_ok=batch_ok, same_pool=same_pool, failed_qc=failed_qc,
                order_confirmed_all=order_confirmed_all,
                all_samples_enclosed=all_samples_enclosed)
        resolved["drift_correct"] = True

    # ---------- max_qc_rsd ----------
    rsd_req = requested_recipe["max_qc_rsd"]
    rsd_prereq_ok = qc_eligible and n_qc >= _MIN_QC_FOR_RSD
    if rsd_req == "auto":
        if rsd_prereq_ok:
            resolved["max_qc_rsd"] = _MAX_QC_RSD_DEFAULT
        else:
            resolved["max_qc_rsd"] = None
            _skip("qc_rsd_filter", _gate_reason(
                batch_ok=batch_ok, same_pool=same_pool, pool_ambiguous=pool_ambiguous,
                n_qc=n_qc, min_qc=_MIN_QC_FOR_RSD, failed_qc=failed_qc))
    elif rsd_req is None:
        resolved["max_qc_rsd"] = None
    else:
        if not rsd_prereq_ok:
            _prerequisite_missing(
                "qc_rsd_filter",
                f"QC-RSDフィルタの前提（単一の確認済みbatch・単一の確認済みpool・健全な"
                f"QC{_MIN_QC_FOR_RSD}件以上）を満たしません。",
                n_qc=n_qc, batch_ok=batch_ok, same_pool=same_pool, failed_qc=failed_qc)
        resolved["max_qc_rsd"] = rsd_req

    # ---------- min_detection_rate ----------
    detection_req = requested_recipe["min_detection_rate"]
    if detection_req and detection_req > 0.0 and getattr(ds, "detected_mask", None) is None:
        _prerequisite_missing(
            "detection_filter",
            "検出状態（gap-fillの区別、detected_mask）が無いためmin_detection_rateを"
            "適用できません。",
            min_detection_rate=detection_req)
    resolved["min_detection_rate"] = detection_req

    # ---------- impute ----------
    # imputeに"auto"は無い（request.pyのスキーマに従う）。conservative-v1は常にそのまま
    # 反映する——half_minが既定で、NaN以外（ゼロ・gap-fill）を新たに欠測へ置換しない
    # という既定挙動はpreprocessing.impute自体の契約（spec §8.1表・数値式は変えない）。
    resolved["impute"] = requested_recipe["impute"]

    assumptions = {
        "n_qc": n_qc, "n_blank": n_blank, "n_sample": n_sample,
        "batch_ok": batch_ok, "same_pool": same_pool, "pool_ambiguous": pool_ambiguous,
        "qc_eligible": qc_eligible, "failed_qc": failed_qc,
        "order_confirmed_all": order_confirmed_all,
        "all_samples_enclosed": all_samples_enclosed,
    }
    if assumptions_normalize_weak is not None:
        # 明示PQNがauto閾値未満のQC件数・失敗疑いQCを抱えたまま続行した記録
        # （spec §8.1「失敗疑いのQCがある場合…対象sample_idを報告する」、
        # spec §8.2「『QC補正済み』等の誤った表現をしない」）。
        assumptions["normalize_weak_qc_reference"] = assumptions_normalize_weak

    return {
        "requested_recipe": requested_recipe,
        "resolved_recipe": resolved,
        "applied_steps": [],  # 計算完了まで空のまま（本taskの中心不変条件）。
        "skipped_steps": skipped_steps,
        "reasons": reasons,
        "policy_version": POLICY_VERSION,
        "assumptions": assumptions,
    }


def _step_should_run(field: str, value: object) -> bool:
    if field == "normalize":
        return value not in (None, "none")
    if field == "drift_correct":
        return bool(value)
    return value is not None  # blank_min_fold / max_qc_rsd


def check_applied_policy(plan: dict, report: dict) -> None:
    """resolved_recipeが要求した各ステップが、実際に適用済みとしてreportへ
    現れているかを照合する。数値層のcaveatだけに任せると、「要求した」と
    「(前提不足で)実施できなかった」の区別が消えるため、ここで二重に確認する。

    合わせて、正規化係数が0・非有限（unscaled_samples）の試料が残っている場合と、
    前処理後に特徴量が0件になった場合を ``DomainError`` にする（spec §8.2）。
    """
    resolved = plan.get("resolved_recipe", {})
    applied = set(report.get("recipe_applied", []))
    skipped = set(report.get("recipe_skipped", []))

    for field, step_name in _STEP_NAMES.items():
        if not _step_should_run(field, resolved.get(field)):
            continue
        if step_name in skipped or step_name not in applied:
            raise DomainError(
                "PREPROCESS_STEP_NOT_APPLIED",
                f"resolved_recipeは{step_name}の適用を要求しましたが、実行結果では"
                "適用されませんでした（数値層のcaveat/status を確認してください）。",
                {"step": step_name, "resolved_recipe": dict(resolved),
                 "recipe_applied": sorted(applied), "recipe_skipped": sorted(skipped),
                 "step_report": (report.get("steps") or {}).get(step_name)},
            )

    normalize_step = (report.get("steps") or {}).get("normalize") or {}
    unscaled = normalize_step.get("unscaled_samples", 0)
    if unscaled:
        raise DomainError(
            "NORMALIZATION_DEGENERATE",
            f"{unscaled}件の試料で正規化係数が0または非有限のため未正規化のまま残って"
            "います。未正規化試料を混ぜたまま比較へは自動進行できません。レシピの見直し"
            "か試料除外が必要です。",
            {"unscaled_samples": unscaled, "n_samples": normalize_step.get("n_samples"),
             "method": normalize_step.get("method")},
        )

    if report.get("features_after") == 0:
        raise DomainError(
            "PREPROCESS_FEATURES_EXHAUSTED",
            "前処理後に特徴量が0件になりました。フィルタ閾値（max_qc_rsd/blank_min_fold）"
            "を見直すか、レシピを変更してください。",
            {"features_before": report.get("features_before"),
             "features_removed_total": report.get("features_removed_total")},
        )
