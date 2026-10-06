"""前処理後のサンプル×特徴量行列に対する差次的解析（純ロジック層、MCP 非依存）。

2群 Welch t 検定・一元配置 ANOVA・log2 fold change・BH-FDR・volcano・交絡検出を提供する。
行列は行=サンプル、列=特徴量。
"""

from __future__ import annotations

import math

import numpy as np


# --- scipy 不在時のフォールバック用の自前分布関数（Numerical Recipes 系、正確） ---
def _betacf(a, b, x):
    """正則化不完全ベータ関数の連分数（Lentz 法）。"""
    MAXIT, EPS, FPMIN = 200, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


def _betai(a, b, x):
    """正則化不完全ベータ関数 I_x(a, b)（t/F 分布の裾確率に使用）。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_beta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(ln_beta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _t_sf_two_sided(t, df):
    """自由度 df の t 分布の両側裾確率 P(|T| >= |t|)。scipy と一致（近似ではない）。"""
    if not math.isfinite(t) or df <= 0:
        return math.nan
    return _betai(df / 2.0, 0.5, df / (df + t * t))


def _f_sf(f, df1, df2):
    """F 分布 F(df1, df2) の上側裾確率 P(F >= f)。scipy と一致（近似ではない）。"""
    if not math.isfinite(f) or f <= 0 or df1 <= 0 or df2 <= 0:
        return math.nan
    return _betai(df2 / 2.0, df1 / 2.0, df2 / (df2 + df1 * f))


def _safe_mean(x):
    """有限値のみで平均。全欠損なら NaN（空スライス警告を避ける）。"""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(x.mean()) if x.size else math.nan


def _log2_transform(arr, pseudo_count):
    """log2(max(x, 0) + pseudo_count)。負値・欠損に頑健。"""
    a = np.asarray(arr, dtype=float)
    return np.log2(np.clip(a, 0.0, None) + pseudo_count)


def bh_fdr(pvalues):
    """Benjamini-Hochberg で p 値を q 値に補正する（NaN は位置保持で除外）。"""
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    finite_idx = np.where(np.isfinite(p))[0]
    if finite_idx.size == 0:
        return q.tolist()
    pv = p[finite_idx]
    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    adj = ranked * m / (np.arange(1, m + 1))
    # 単調化（後ろから最小を累積）
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty(m)
    out[order] = adj
    q[finite_idx] = out
    return q.tolist()


def welch_t(a, b):
    """Welch t 統計量と両側 p 値（純粋な数値部品）。`a` が基準、`b` が比較で、
    t は `b − a` の向き（`_welch_t` 参照）。

    v1 の `two_group_test` と v2 の `statistics_v2` が**同じ実装**を使うための
    公開名。v1 側の変換・効果量の定義（pseudo_count と log 空間の平均差）は
    v2 とは別物なので共有しない——共有するのはここ（検定統計量）と
    `bh_fdr` だけで、v1 の数値契約は変えない。
    """
    return _welch_t(a, b)


def _welch_t(a, b):
    """Welch t 統計量と両側 p 値。

    `a` が基準（対照）、`b` が比較対象。**t = (mean(b) − mean(a)) / SE** で、
    log2FC（正 = b が高い）・Tukey の mean_difference（test − reference）と同じ向き。
    p 値は scipy があれば使い、無ければ自前の t 分布裾確率（こちらも正確）。
    """
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return math.nan, math.nan
    va, vb = a.var(ddof=1), b.var(ddof=1)
    na, nb = a.size, b.size
    denom = va / na + vb / nb
    if denom <= 0:
        return math.nan, math.nan
    t = (b.mean() - a.mean()) / math.sqrt(denom)
    df = denom ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    try:
        from scipy import stats
        p = 2.0 * stats.t.sf(abs(t), df)
    except Exception:
        # scipy 不在時は自前の t 分布両側裾確率（正確。旧実装の正規近似は小 df で
        # 非保守＝偽陽性寄りだった）。
        p = _t_sf_two_sided(t, df)
    return t, p


def _log2fc(mean_a, mean_b, pseudo_count):
    """log2 fold change。**正なら group_b が高い（＝上昇）**。

    慣習（log2FC = log2(比較対象 / 基準)）に合わせる。group_a が基準（対照）、
    group_b が比較対象である。2026-08-31 に向きを反転した（旧実装は
    正 = group_a が高い、で慣習と逆だった）。過去の解析結果とは符号が逆になる。
    """
    num = max(mean_b, 0.0) + pseudo_count
    den = max(mean_a, 0.0) + pseudo_count
    return math.log2(num / den)


def two_group_test(matrix, feature_names, group_labels, group_a, group_b,
                   *, log2=True, pseudo_count=1.0, log_transform=False):
    """群 a/b について特徴量ごとに Welch t 検定と log2 fold change を計算する。

    **log2fc は正なら group_b が高い（上昇）。** group_a が基準（対照）、
    group_b が比較対象である。log_transform=True のときは log2(x+pseudo_count)
    空間で検定し、log2FC も log2 空間の群平均差（＝幾何平均比の log2）とする。
    MS 強度は対数正規に近く、生強度での t 検定は正規性仮定を外れやすいため、
    こちらが推奨経路。mean_a/mean_b は解釈用に常に生強度平均を返す。
    """
    matrix = np.asarray(matrix, dtype=float)
    labels = np.asarray(group_labels)
    idx_a = np.where(labels == group_a)[0]
    idx_b = np.where(labels == group_b)[0]
    results = []
    for j, name in enumerate(feature_names):
        col = matrix[:, j]
        a, b = col[idx_a], col[idx_b]
        mean_a = _safe_mean(a)
        mean_b = _safe_mean(b)
        if log_transform:
            la, lb = _log2_transform(a, pseudo_count), _log2_transform(b, pseudo_count)
            t, p = _welch_t(la, lb)
            lm_a, lm_b = _safe_mean(la), _safe_mean(lb)
            fc = (lm_b - lm_a) if (log2 and math.isfinite(lm_a) and math.isfinite(lm_b)) else math.nan
        else:
            t, p = _welch_t(a, b)
            fc = _log2fc(mean_a, mean_b, pseudo_count) if (log2 and math.isfinite(mean_a) and math.isfinite(mean_b)) else math.nan
        results.append({
            "feature": name, "mean_a": mean_a, "mean_b": mean_b,
            "log2fc": fc, "t": t, "p": p,
        })
    return results


def _one_way_f(groups):
    """一元配置 ANOVA の F 統計量と p 値（scipy があれば使用）。"""
    groups = [g[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if g.size >= 2]
    k = len(groups)
    if k < 2:
        return math.nan, math.nan, k
    grand = np.concatenate(groups)
    n = grand.size
    grand_mean = grand.mean()
    ss_between = sum(g.size * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_within = sum(((g - g.mean()) ** 2).sum() for g in groups)
    df_b, df_w = k - 1, n - k
    if df_w <= 0 or ss_within == 0:
        return math.nan, math.nan, k
    f = (ss_between / df_b) / (ss_within / df_w)
    try:
        from scipy import stats
        p = stats.f.sf(f, df_b, df_w)
    except Exception:
        # scipy 不在時は自前の F 分布上側裾確率（不完全ベータ、任意 df で正確）。
        # 旧実装の chi2 級数近似は df_b が奇数だと誤りだった。
        p = _f_sf(f, df_b, df_w)
    return f, p, k


def one_way_anova(matrix, feature_names, group_labels, *,
                  log_transform=False, pseudo_count=1.0):
    """factor の全水準について特徴量ごとに一元配置 ANOVA を実行する。

    log_transform=True のときは log2(x+pseudo_count) 空間で検定する（two_group_test と同様、
    MS 強度の対数正規性に配慮）。
    """
    matrix = np.asarray(matrix, dtype=float)
    labels = np.asarray(group_labels)
    levels = [lv for lv in sorted(set(labels.tolist())) if lv is not None]
    results = []
    for j, name in enumerate(feature_names):
        col = matrix[:, j]
        groups = [col[labels == lv] for lv in levels]
        if log_transform:
            groups = [_log2_transform(g, pseudo_count) for g in groups]
        f, p, k = _one_way_f(groups)
        results.append({"feature": name, "F": f, "p": p, "n_groups": k})
    return results


def add_fdr(results):
    """results の p 列に BH-FDR を適用し、各要素に q を付与する。"""
    qs = bh_fdr([r["p"] for r in results])
    for r, q in zip(results, qs):
        r["q"] = q
    return results


def volcano_data(results, q_thr=0.05, log2fc_thr=1.0):
    """volcano 用の点列（log2fc, -log10 p, 有意フラグ up/down/ns）を返す。"""
    points = []
    for r in results:
        p = r.get("p")
        q = r.get("q", p)
        fc = r.get("log2fc")
        neg_log10_p = (-math.log10(p)) if (p is not None and math.isfinite(p) and p > 0) else math.nan
        sig = "ns"
        if q is not None and math.isfinite(q) and q <= q_thr and fc is not None and math.isfinite(fc):
            if fc >= log2fc_thr:
                sig = "up"
            elif fc <= -log2fc_thr:
                sig = "down"
        points.append({"feature": r["feature"], "log2fc": fc,
                       "neg_log10_p": neg_log10_p, "sig": sig})
    return points


def check_confounding(group_labels, batch_labels):
    """各群が単一バッチに偏るか（群 ⟂ バッチの交絡）を判定する。

    ``assessable`` は交絡の有無を判定できたかを示す。バッチ情報が無い、または
    全試料が単一バッチのときは交絡を評価できない（``assessable=False``）。
    「バッチが1つしかない＝交絡なし」と誤読しないよう明示的に分ける。
    """
    if not batch_labels or all(b is None for b in batch_labels):
        return {"confounded": False, "assessable": False,
                "detail": "バッチ情報が無いため交絡判定不可。"}
    known_batches = {b for b in batch_labels if b is not None}
    if len(known_batches) < 2:
        return {"confounded": False, "assessable": False,
                "detail": "全試料が単一バッチ（判明分）のため、群⟂バッチ交絡は評価できません。"}
    by_group: dict[str, set] = {}
    for g, b in zip(group_labels, batch_labels):
        if b is None:
            continue
        by_group.setdefault(g, set()).add(b)
    single = {g: next(iter(bs)) for g, bs in by_group.items() if len(bs) == 1}
    confounded = len(single) == len(by_group) and len(set(single.values())) > 1
    if confounded:
        detail = "各群が単一バッチに対応し、処理効果と測定バッチを分離できません: " + \
                 ", ".join(f"{g}->{b}" for g, b in single.items())
    else:
        detail = "群とバッチは交絡していません（または部分的）。"
    return {"confounded": confounded, "assessable": True, "detail": detail}


def summarize_two_group(results, q_thr=0.05, log2fc_thr=1.0, top_n=15):
    """2群比較結果の有意 up/down 件数と上位特徴量を要約する。"""
    tested = [r for r in results if r.get("p") is not None and math.isfinite(r["p"])]

    def is_sig(r):
        q = r.get("q", r.get("p"))
        return q is not None and math.isfinite(q) and q <= q_thr and \
            r.get("log2fc") is not None and abs(r["log2fc"]) >= log2fc_thr

    sig = [r for r in tested if is_sig(r)]
    n_up = sum(1 for r in sig if r["log2fc"] >= log2fc_thr)
    n_down = sum(1 for r in sig if r["log2fc"] <= -log2fc_thr)
    top = sorted(sig, key=lambda r: r.get("q", r.get("p", 1.0)))[:top_n]
    return {"n_tested": len(tested), "n_significant": len(sig),
            "n_up": n_up, "n_down": n_down, "top": top}
