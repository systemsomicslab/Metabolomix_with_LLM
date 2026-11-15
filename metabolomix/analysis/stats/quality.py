"""品質の純計算（最初の概観の品質の段。spec 2026-11-06-dataset-overview-design §4.1）。

StudyData（前処理の前の値）だけを読み、QC の RSD の分布・失敗した QC 注入・生体試料の実検出率と
gap-fill の割合・ドリフトの兆候を数え、交絡の報告（`confounding.confounding_report`）と合わせて
品質の判定（ok / caution / poor / not_evaluable）を出す。値は変えない（検出だけ）。

セッションにも MCP にも依存しない（`metabolomix/analysis/stats/` の規則）。失敗した QC 注入の
規則（旧 `analysis.preprocessing.detect_failed_qc` と同じもの）はここだけが持つ。

陽陰を統合したデータセット（`StudyData.feature_polarity`）は、注入ごとの総量を極性ごとに数える
（失敗した QC 注入・ドリフトの兆候）。陽イオンと陰イオンは別の注入なので、強度の大きい極性の総量が
もう一方の失敗やドリフトを隠さないように。ドリフトの兆候はその極性の測定の注入順・バッチで見る。
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats as sps

from metabolomix.analysis.stats.data import (DETECTED, GAP_FILLED, MISSING, StudyData, cell_states,
                                             order_may_be_file_order, polarity_samples)

__all__ = ["QUALITY_THRESHOLDS", "qc_rsd_summary", "failed_qc", "detection_summary",
           "drift_signs", "quality_verdict"]

#: 判定の閾値（概観の記録の settings に書く）。
QUALITY_THRESHOLDS = {
    "qc_rsd_ok": 0.30,          # RSD がこれ以下の特徴を「良い」と数え、中央値がこれを超えたら poor
    "failed_qc_ratio": 0.2,     # QC の総量が QC の中央値のこの倍未満なら失敗した注入
    "drift_min_qc": 5,          # ドリフトの兆候を見るバッチの QC の最小数
    "drift_rho": 0.6,           # QC の総量の log2 と注入順の Spearman の |ρ|
    "drift_p": 0.05,
}
_NO_BATCH = "(none)"            # preprocess._batch_key と同じ鍵


def _rows(study: StudyData, role: str) -> list[int]:
    return [i for i, s in enumerate(study.samples) if s.role == role and s.include]


def _positive(block: np.ndarray) -> np.ndarray:
    b = np.asarray(block, dtype=float)
    return np.where(np.isfinite(b) & (b > 0), b, np.nan)


def _polarity_columns(study: StudyData) -> list[tuple[str, np.ndarray]]:
    """統合したデータセットの (極性, 列の位置)（特徴の順に出てきた極性の順）。統合でなければ空。"""
    if not study.feature_polarity:
        return []
    labels = [study.feature_polarity.get(f, "(unknown)") for f in study.feature_ids]
    return [(label, np.array([j for j, lab in enumerate(labels) if lab == label], dtype=int))
            for label in dict.fromkeys(labels)]


def qc_rsd_summary(study: StudyData) -> dict:
    """QC の特徴ごとの RSD（有限正値、ddof=1。前処理の qc_rsd の足切りと同じ定義）の分布。"""
    qc = _rows(study, "qc")
    if len(qc) < 2:
        return {"status": "not_evaluable", "reason": "too_few_qc", "n_qc": len(qc)}
    pos = _positive(study.values[qc])
    n = np.isfinite(pos).sum(axis=0)
    with np.errstate(all="ignore"):
        mean = np.nanmean(np.where(n >= 2, pos, np.nan), axis=0)
        sd = np.nanstd(np.where(n >= 2, pos, np.nan), axis=0, ddof=1)
        rsd = sd / mean
    ok = np.isfinite(rsd)
    if not ok.any():
        return {"status": "not_evaluable", "reason": "no_computable_feature", "n_qc": len(qc)}
    vals = rsd[ok]
    return {"status": "ok", "n_qc": len(qc), "n_computable": int(ok.sum()),
            "median": float(np.median(vals)),
            "frac_le_030": float((vals <= QUALITY_THRESHOLDS["qc_rsd_ok"]).mean())}


def _failed_in(study: StudyData, qc: list[int], block: np.ndarray) -> list[str]:
    totals = np.nansum(_positive(block), axis=1)
    median = float(np.median(totals))
    if not math.isfinite(median) or median <= 0:
        return []
    cut = QUALITY_THRESHOLDS["failed_qc_ratio"] * median
    return sorted(study.samples[i].sample_id for i, t in zip(qc, totals) if t < cut)


def failed_qc(study: StudyData) -> dict:
    """総量が QC の中央値の `failed_qc_ratio` 倍未満の QC 注入（試料 ID）。QC が 3 本未満は見ない。

    統合したデータセットは極性ごとに総量を数え、どれかの極性で失敗した注入を `failed` に、極性ごとの内訳を
    `by_polarity` に返す（統合の試料 ID は陽陰の 2 本の注入を表す）。
    """
    qc = _rows(study, "qc")
    if len(qc) < 3:
        return {"checked": len(qc), "failed": []}
    blocks = _polarity_columns(study)
    if not blocks:
        return {"checked": len(qc), "failed": _failed_in(study, qc, study.values[qc])}
    by = {label: _failed_in(study, qc, study.values[np.ix_(qc, cols)]) for label, cols in blocks}
    return {"checked": len(qc), "failed": sorted(set().union(*by.values())), "by_polarity": by}


def detection_summary(study: StudyData) -> dict:
    """生体試料のセルの内訳（実検出・gap-fill・欠測の割合）。検出の別が無ければ評価できない。"""
    if study.detected is None:
        return {"status": "not_evaluable", "reason": "no_detected_mask"}
    bio = _rows(study, "sample")
    if not bio or not study.values.shape[1]:
        return {"status": "not_evaluable", "reason": "no_biological_cells"}
    states = cell_states(study.values[bio], np.asarray(study.detected, bool)[bio])
    total = states.size
    return {"status": "ok", "n_cells": int(total),
            "detected_frac": float((states == DETECTED).sum() / total),
            "gap_filled_frac": float((states == GAP_FILLED).sum() / total),
            "missing_frac": float((states == MISSING).sum() / total)}


def _drift_signs_in(study: StudyData, samples, cols: np.ndarray | None) -> list[dict]:
    by_batch: dict[str, list[int]] = {}
    for i in _rows(study, "qc"):
        s = samples[i]
        if s.injection_order is None:
            continue
        by_batch.setdefault(s.batch if s.batch is not None else _NO_BATCH, []).append(i)
    out = []
    for batch in sorted(by_batch):
        rows = by_batch[batch]
        if len(rows) < QUALITY_THRESHOLDS["drift_min_qc"]:
            continue
        block = study.values[rows] if cols is None else study.values[np.ix_(rows, cols)]
        totals = np.nansum(_positive(block), axis=1)
        keep = totals > 0
        if keep.sum() < QUALITY_THRESHOLDS["drift_min_qc"]:
            continue
        order = np.array([samples[i].injection_order for i in rows], float)[keep]
        rho, p = sps.spearmanr(order, np.log2(totals[keep]))
        rho, p = float(rho), float(p)
        if not (math.isfinite(rho) and math.isfinite(p)):
            continue
        flag = abs(rho) >= QUALITY_THRESHOLDS["drift_rho"] and p < QUALITY_THRESHOLDS["drift_p"]
        # ρ に使った QC の注入順が実際の測定順とは限らない（MS-DIAL の既定・シートの宣言でファイル名の順）なら印を
        # 付ける（ρ はファイル名の並びとの相関かもしれない）。確かな順の行には鍵ごと付けない
        unsure = any(order_may_be_file_order(samples[i]) for i, k in zip(rows, keep) if k)
        out.append({"batch": batch, "n_qc": int(keep.sum()), "rho": rho, "p": p, "flag": bool(flag),
                    **({"injection_order_may_be_file_order": True} if unsure else {})})
    return out


def drift_signs(study: StudyData) -> list[dict]:
    """バッチごと（QC が `drift_min_qc` 本以上）の、QC の総量の log2 と注入順の Spearman。

    統合したデータセットは極性ごとに、その極性の特徴の総量とその極性の測定の注入順・バッチで見て、
    各行に `polarity` を足す。ρ に使った QC の注入順が実際の測定順とは限らない（`data.order_may_be_file_order`）
    行には `injection_order_may_be_file_order: True` を足す（ρ は出す。判定はそのまま）。
    """
    blocks = _polarity_columns(study)
    if not blocks:
        return _drift_signs_in(study, study.samples, None)
    return [{"polarity": label, **row} for label, cols in blocks
            for row in _drift_signs_in(study, polarity_samples(study, label), cols)]


def quality_verdict(*, rsd: dict, failed: dict, drift: list[dict], confounding: dict,
                    replicate_warnings: list) -> tuple[str, list[str]]:
    """品質の判定と理由のコード（計画書 2026-11-09-a-overview「品質の判定」）。"""
    poor, caution = [], []
    status = confounding.get("batch_status")
    inj = (confounding.get("injection_order") or {}).get("status")
    if status == "complete":
        poor.append("complete_confounding")
    if failed.get("failed"):
        poor.append("failed_qc")
    if rsd.get("status") == "ok" and rsd["median"] > QUALITY_THRESHOLDS["qc_rsd_ok"]:
        poor.append("qc_rsd_high")
    if status == "partial":
        caution.append("partial_confounding")
    if inj == "biased":
        caution.append("injection_order_biased")
    if any(d.get("flag") for d in drift):
        caution.append("drift_sign")
    if replicate_warnings:
        caution.append("possible_technical_replicates")
    qc_missing = rsd.get("status") != "ok"
    tail = ["qc_not_evaluable"] if qc_missing else []
    if poor:
        return "poor", poor + caution + tail
    if caution:
        return "caution", caution + tail
    if qc_missing and status == "not_evaluable":
        return "not_evaluable", ["no_qc_and_batch_unknown"]
    if qc_missing:
        # QC が無い（足りない）ときは ok と言わない（spec §4.1: 評価できないことをそのまま書く）
        return "not_evaluable", ["qc_not_evaluable"]
    return "ok", []
