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

from metabolomix.analysis.stats.data import (DETECTED, GAP_FILLED, MISSING, PRESENT_UNKNOWN, StudyData,
                                             cell_states, order_may_be_file_order, polarity_samples)

__all__ = ["QUALITY_THRESHOLDS", "REPLICATE_BLOCK_THRESHOLDS", "qc_rsd_summary", "failed_qc",
           "detection_summary", "drift_signs", "replicate_blocks", "quality_verdict"]

#: 判定の閾値（概観の記録の settings に書く）。
QUALITY_THRESHOLDS = {
    "qc_rsd_ok": 0.30,          # RSD がこれ以下の特徴を「良い」と数え、中央値がこれを超えたら poor
    "failed_qc_ratio": 0.2,     # QC の総量が QC の中央値のこの倍未満なら失敗した注入
    "drift_min_qc": 5,          # ドリフトの兆候を見るバッチの QC の最小数
    "drift_rho": 0.6,           # QC の総量の log2 と注入順の Spearman の |ρ|
    "drift_p": 0.05,
}
_NO_BATCH = "(none)"            # preprocess._batch_key と同じ鍵
#: 技術反復らしいまとまり（`replicate_blocks`）の閾値。独立な試料で鳴らないことを優先して決めた（合成データで
#: 群の中の性差・注入順のドリフト・偶然の近さを試し、どれも 1 つ以上の条件で落ちる値）。
REPLICATE_BLOCK_THRESHOLDS = {
    "min_blocks": 3,            # まとまりの数（n / k）の下限。2 つだと前半・後半のドリフトと見分けられない
    "min_features": 30,         # 群のすべての注入で実検出の特徴がこれ未満なら見ない
    "max_features": 2000,       # 使う特徴の上限（強度の中央値の高い順）
    "margin": 0.3,              # まとまりの中の相関の中央値 − まとまりの間の相関の中央値
    "purity": 0.8,              # 各注入の相関の高い k−1 本のうちまとまりの仲間の割合（注入の平均）
    "median_icc": 0.4,          # 特徴ごとの級内相関 ICC(1) の中央値（大半の特徴でまとまりが分散を説明する）
    "adjacent_excess": 0.5,     # 並びで k 本未満しか離れていない別のまとまりの対の相関が、間の中央値を margin の何割まで超えてよいか
    "single_group_untrended": 0.25,  # 見られる群が 1 つだけのとき、まとまりの平均の差のうち並びの直線で説明できない割合の中央値の下限
}


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


def _group_order(study: StudyData, rows: list[int]) -> tuple[list[int], str]:
    """群の注入の並び: すべてに注入順があればその順（同じ値は試料の表の順）、無ければ試料の表の順。"""
    orders = [study.samples[i].injection_order for i in rows]
    if all(o is not None for o in orders):
        return [i for _o, i in sorted(zip(orders, rows))], "injection_order"
    return list(rows), "table"


def _block_profiles(study: StudyData, rows: list[int]) -> tuple[np.ndarray | None, int, str]:
    """群の注入 × 特徴の標準化した log2 と、使った特徴の数・セルの扱い（`detected` / `present`）。

    特徴は群のすべての注入で実検出（gap-fill を含まない。検出の別が無ければ値あり）のものを、強度の中央値の高い順に
    `max_features` まで。注入ごとに全特徴が同じだけずれる分（希釈・総量）を中央値で外してから、特徴ごとに平均 0・SD 1 に
    する（強度の桁の違いで相関が 1 に張り付かず、少数の大きくばらつく特徴が相関を決めないように）。足りなければ None。
    """
    th = REPLICATE_BLOCK_THRESHOLDS
    vals = np.asarray(study.values[rows], dtype=float)
    det = None if study.detected is None else np.asarray(study.detected, dtype=bool)[rows]
    states = cell_states(vals, det)
    basis = "detected" if det is not None else "present"
    usable = np.flatnonzero((states == (DETECTED if det is not None else PRESENT_UNKNOWN)).all(axis=0))
    if usable.size < th["min_features"]:
        return None, int(usable.size), basis
    x = np.log2(vals[:, usable])
    med = np.median(x, axis=0)
    if usable.size > th["max_features"]:
        top = np.argsort(-med, kind="stable")[:th["max_features"]]
        x, med = x[:, top], med[top]
    x = x - np.median(x - med, axis=1, keepdims=True)
    sd = x.std(axis=0, ddof=1)
    keep = np.isfinite(sd) & (sd > 1e-12)
    if keep.sum() < th["min_features"]:
        return None, int(keep.sum()), basis
    x = x[:, keep]
    return (x - x.mean(axis=0)) / sd[keep], int(keep.sum()), basis


def _sample_corr(z: np.ndarray) -> np.ndarray:
    """注入どうしの Pearson 相関（標準化した特徴の上で）。"""
    zc = z - z.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(zc, axis=1, keepdims=True)
    zc = zc / np.where(norm > 0, norm, 1.0)
    return zc @ zc.T


def _greedy_blocks(c: np.ndarray, k: int) -> list[list[int]]:
    """相関の高い注入を k 本ずつ束ねる（種は、残りの中で近い k−1 本との相関の平均が最も高い注入）。"""
    left = list(range(len(c)))
    blocks = []
    while len(left) > k:
        sub = c[np.ix_(left, left)].copy()
        np.fill_diagonal(sub, -np.inf)
        near = np.argsort(-sub, axis=1, kind="stable")[:, :k - 1]
        seed = int(np.argmax(np.take_along_axis(sub, near, axis=1).mean(axis=1)))
        chosen = sorted([left[seed], *(left[j] for j in near[seed])])
        blocks.append(chosen)
        left = [i for i in left if i not in chosen]
    blocks.append(left)
    return sorted(blocks)


def _layouts(c: np.ndarray, k: int):
    """試すまとまりの分け方（並びの位置で）: 並びで連続・並びで巡回（1, 2, 3, 1, 2, 3, …）・相関からの束ね。"""
    n = len(c)
    m = n // k
    yield "consecutive", [list(range(b * k, (b + 1) * k)) for b in range(m)]
    yield "interleaved", [list(range(b, n, m)) for b in range(m)]
    yield "data", _greedy_blocks(c, k)


def _block_stats(z: np.ndarray, c: np.ndarray, blocks: list[list[int]], k: int) -> dict:
    """分け方の数: 中と間の相関の中央値・差（margin）・近い k−1 本の仲間の割合（purity）・特徴ごとの ICC(1) の中央値・
    並びで k 本未満しか離れていない別のまとまりの対の相関の超過（`adjacent_excess` = `(その中央値 − 間の中央値) / margin`。
    注入順のドリフトが作った連続のまとまりなら、境目をまたぐ近い対もまとまりの中と同じくらい似るので 1 に近い）・
    まとまりの平均の差のうち並びの位置の直線で説明できない割合の中央値（`untrended`。直線のドリフトなら 0 に近い）。"""
    n, m = len(c), len(blocks)
    label = np.empty(n, dtype=int)
    for b, members in enumerate(blocks):
        label[members] = b
    same = label[:, None] == label[None, :]
    off = ~np.eye(n, dtype=bool)
    within = float(np.median(c[same & off]))
    between = float(np.median(c[~same]))
    margin = within - between
    near = np.argsort(-np.where(off, c, -np.inf), axis=1, kind="stable")[:, :k - 1]
    purity = float(same[np.arange(n)[:, None], near].mean())
    arr = np.stack([z[members] for members in blocks])                    # m × k × 特徴
    means = arr.mean(axis=1)
    msb = k * ((means - arr.mean(axis=(0, 1))) ** 2).sum(axis=0) / (m - 1)
    msw = ((arr - means[:, None, :]) ** 2).sum(axis=(0, 1)) / (m * (k - 1))
    den = msb + (k - 1) * msw
    icc = np.divide(msb - msw, den, out=np.zeros_like(den), where=den > 0)
    near_pairs = [c[p, q] for p in range(n) for q in range(p + 1, min(n, p + k)) if label[p] != label[q]]
    excess = (float(np.median(near_pairs)) - between) / margin if near_pairs and margin > 0 else None
    # まとまりの平均の並びのうち、並びの位置の直線で説明できない割合（特徴ごと）。ドリフトなら 0 に近い
    fitted = means[label]                                                   # n × 特徴（その注入のまとまりの平均）
    pos = np.arange(n, dtype=float) - (n - 1) / 2
    ssb = (fitted ** 2).sum(axis=0)
    lin = (pos @ fitted) ** 2 / float(pos @ pos)
    untrended = np.divide(ssb - lin, ssb, out=np.zeros_like(ssb), where=ssb > 0)
    return {"within": within, "between": between, "margin": margin, "purity": purity,
            "median_icc": float(np.median(icc)), "adjacent_excess": excess,
            "untrended": float(np.median(untrended))}


def _passes(stats: dict, *, single: bool) -> bool:
    """採る条件。`single`（見られる群が 1 つだけ）なら、まとまりの平均の差が並びの直線（ドリフト）で説明しきれないことも求める。"""
    th = REPLICATE_BLOCK_THRESHOLDS
    return (stats["margin"] >= th["margin"] and stats["purity"] >= th["purity"]
            and stats["median_icc"] >= th["median_icc"]
            and stats["adjacent_excess"] is not None and stats["adjacent_excess"] <= th["adjacent_excess"]
            and (not single or stats["untrended"] >= th["single_group_untrended"]))


def replicate_blocks(study: StudyData) -> dict:
    """subject_id が無いとき、群の中の注入が同じ個体の技術反復らしい k 本ずつのまとまりに分かれるか（データから。まとめはしない）。

    試料名に反復の印が無くても（`data.replicate_name_warnings` が見つけない）、群の注入が「個体 m 匹 × 注入 k 本」なら、
    同じ個体の注入どうしはほぼすべての特徴で似る。含める生体試料（`role == "sample"`）を群ごとに注入順（無ければ試料の表の
    順）に並べ、`_block_profiles` の標準化した log2 の相関で、k が群の注入の数を割り切り、まとまりが `min_blocks` 以上になる
    k ごとに、分け方（並びで連続 → 並びで巡回 → 相関からの束ね）を順に試し、最初に次のすべてを満たしたものを採る:
    中と間の相関の中央値の差 ≥ `margin`、各注入の近い k−1 本のうちまとまりの仲間の割合の平均 ≥ `purity`（性差など群の中の
    大きな部分集団は割合が下がる）、特徴ごとの ICC(1) の中央値 ≥ `median_icc`（一部の特徴だけが動く部分集団は中央値が
    0 近くに残る）、並びで k 本未満しか離れていない別のまとまりの対が間の対より似ていない（`adjacent_excess`。注入順の
    ドリフトは連続のまとまりを作るが、境目をまたぐ近い対も似る）。

    群をまたいで同じ k が、その k を当てられるすべての群（2 群以上評価できたなら 2 群以上）で通ったときだけ `found`。
    見られる群が 1 つだけなら群をまたいだ一貫性が使えないので、まとまりの平均の差の大半が並びの直線（ドリフト）で
    説明できないこと（`untrended` ≥ `single_group_untrended`）も求める（強いドリフトは連続の両端のまとまりを締まって見せる）。
    複数の k が通れば margin の平均の大きい方（同じなら小さい k）。

    戻り値の `status`: `found`（`k`・`groups`〔群ごとの `blocks`＝試料 ID の組・`layout`・`order`・数〕・`skipped`）/
    `not_found` / `not_evaluable`（見られる群が無い）/ `subject_id_set`（含める生体試料のどれかに subject_id がある。見ない）。
    """
    th = REPLICATE_BLOCK_THRESHOLDS
    bio = [i for i, s in enumerate(study.samples) if s.include and s.role == "sample"]
    if any(study.samples[i].subject_id is not None for i in bio):
        return {"status": "subject_id_set"}
    by_group: dict = {}
    for i in bio:
        by_group.setdefault(study.samples[i].group, []).append(i)
    ready, skipped = [], []
    for group, rows in by_group.items():
        n = len(rows)
        ks = [k for k in range(2, n // th["min_blocks"] + 1) if n % k == 0]
        if not ks:
            skipped.append({"group": group, "n": n, "reason": "no_block_size"})
            continue
        ordered, order_basis = _group_order(study, rows)
        z, n_features, cells = _block_profiles(study, ordered)
        if z is None:
            skipped.append({"group": group, "n": n, "reason": "too_few_features", "n_features": n_features})
            continue
        ready.append((group, n, ks, ordered, order_basis, cells, n_features, z))
    single = len(ready) == 1                 # 群をまたいだ一貫性が使えない: ドリフトの直線も見る
    evaluated = []
    for group, n, ks, ordered, order_basis, cells, n_features, z in ready:
        c = _sample_corr(z)
        passing = {}
        for k in ks:
            for layout, blocks in _layouts(c, k):
                stats = _block_stats(z, c, blocks, k)
                if _passes(stats, single=single):
                    passing[k] = {"layout": layout, **stats,
                                  "blocks": [[study.samples[ordered[p]].sample_id for p in b] for b in blocks]}
                    break
        evaluated.append({"group": group, "n": n, "order": order_basis, "cells": cells,
                          "n_features": n_features, "ks": ks, "passing": passing})
    if not evaluated:
        return {"status": "not_evaluable", "skipped": skipped, "thresholds": th}
    found = []
    for k in sorted({k for g in evaluated for k in g["ks"]}):
        applicable = [g for g in evaluated if k in g["ks"]]
        if not all(k in g["passing"] for g in applicable):
            continue
        if len(evaluated) >= 2 and len(applicable) < 2:
            continue
        found.append((float(np.mean([g["passing"][k]["margin"] for g in applicable])), -k, applicable))
    if not found:
        return {"status": "not_found", "n_groups": len(evaluated), "skipped": skipped, "thresholds": th}
    _margin, neg_k, applicable = max(found, key=lambda x: (x[0], x[1]))
    k = -neg_k
    groups = [{"group": g["group"], "n": g["n"], "order": g["order"], "cells": g["cells"],
               "n_features": g["n_features"], **g["passing"][k]} for g in applicable]
    others = [{"group": g["group"], "n": g["n"], "reason": "block_size_does_not_divide"}
              for g in evaluated if k not in g["ks"]]
    return {"status": "found", "k": k, "groups": groups, "skipped": skipped + others, "thresholds": th}


def quality_verdict(*, rsd: dict, failed: dict, drift: list[dict], confounding: dict,
                    replicate_warnings: list, replicate_blocks: dict | None = None) -> tuple[str, list[str]]:
    """品質の判定と理由のコード（計画書 2026-11-09-a-overview「品質の判定」）。

    `replicate_blocks`（`replicate_blocks` の戻り値）が `found` なら caution の理由 `possible_replicate_blocks`。"""
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
    if (replicate_blocks or {}).get("status") == "found":
        caution.append("possible_replicate_blocks")
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
