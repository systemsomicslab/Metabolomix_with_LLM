"""MS-DIAL ロード後のサンプル×特徴量行列に対する前処理・QC（純ロジック層、MCP 非依存）。

役割検出・ブランク処理・正規化・QCドリフト補正・特徴量フィルタ・欠損補完を提供する。
lipidmix.corpus.knowledge_store / lipidmix.msdial.peak_verification と同じく MCP に依存しない純関数群。
行列は行=サンプル、列=特徴量（test_arf.build_pca_matrix の向き）。
"""

from __future__ import annotations

import re

import numpy as np

DEFAULT_ROLE_TOKENS: dict[str, set[str]] = {
    "qc": {"qc"},
    "blank": {"blank"},
}

# QC 層別キーの算出で除くトークン（極性）。日付/純数字/QC 自体は別途除外。
_POLARITY_TOKENS = {"neg", "pos"}


def _segments(text: str) -> set[str]:
    return {seg for seg in str(text).casefold().split("_") if seg}


def detect_sample_roles(
    sample_names: list[str],
    class_ids: dict[str, str] | None = None,
    config: dict | None = None,
) -> dict[str, str]:
    """各サンプルを sample/qc/blank に分類する。

    ファイル名と Class ID の `_` 区切りセグメントを大小無視でトークン照合する。
    blank を qc より優先評価する（曖昧語衝突を避けるため決定的順序）。
    """
    cfg = config or {}
    tokens = {
        "qc": {t.casefold() for t in cfg.get("qc_tokens", DEFAULT_ROLE_TOKENS["qc"])},
        "blank": {t.casefold() for t in cfg.get("blank_tokens", DEFAULT_ROLE_TOKENS["blank"])},
    }
    class_ids = class_ids or {}
    roles: dict[str, str] = {}
    for name in sample_names:
        segs = _segments(name) | _segments(class_ids.get(name, ""))
        if segs & tokens["blank"]:
            roles[name] = "blank"
        elif segs & tokens["qc"]:
            roles[name] = "qc"
        else:
            roles[name] = "sample"
    return roles


def detect_qc_strata(sample_names: list[str], roles: dict[str, str]) -> set[str]:
    """QC 試料名から層別サブグループのキー集合を返す。

    日付(8桁)・`qc`・極性(neg/pos)・純数字(反復や連番タイムスタンプ)トークンを
    除いた残りを層識別子とする。例: `20240311_QC_Cerebellum_ICR_NEG_1` → `cerebellum_icr`。
    層が2つ以上なら「プールQC が層別」と判断でき、全 QC を1系列扱いするドリフト補正/
    RSD が近似になる旨の caveat を出す材料になる（設計 §3.4）。
    """
    strata: set[str] = set()
    for name in sample_names:
        if roles.get(name) != "qc":
            continue
        toks = []
        for tok in str(name).split("_"):
            if not tok:
                continue
            low = tok.lower()
            if low == "qc" or low in _POLARITY_TOKENS:
                continue
            if re.fullmatch(r"\d{8}", tok) or tok.isdigit():  # 日付 / 反復・連番
                continue
            toks.append(low)
        strata.add("_".join(toks))
    return strata


def _reference_rows(matrix, roles, sample_names):
    """QC 行が特定できれば QC のみ、無ければ全行を参照集合として返す。"""
    if roles and sample_names:
        qc_idx = [i for i, n in enumerate(sample_names) if roles.get(n) == "qc"]
        if qc_idx:
            return matrix[qc_idx], True
    return matrix, False


def _rows_for_role(matrix, roles, sample_names, role):
    """指定ロールに該当する行のインデックスを抽出し、対応する行部分行列を返す。"""
    idx = [i for i, n in enumerate(sample_names) if roles.get(n) == role]
    return matrix[idx] if idx else None


def blank_filter(matrix, roles, sample_names, min_fold=3.0):
    """生体試料平均 < min_fold × ブランク平均 の特徴量を背景として除去する。"""
    matrix = np.asarray(matrix, dtype=float)
    n_features = matrix.shape[1]
    blanks = _rows_for_role(matrix, roles, sample_names, "blank")
    samples = _rows_for_role(matrix, roles, sample_names, "sample")
    if blanks is None or samples is None:
        return np.ones(n_features, dtype=bool), {
            "removed": 0, "status": "skipped",
            "caveat": "ブランクまたは生体試料が無いため背景除去は未実施。",
        }
    blank_mean = np.nanmean(blanks, axis=0)
    sample_mean = np.nanmean(samples, axis=0)
    # ブランクがゼロの特徴量は常に keep（除算回避）
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(blank_mean > 0, sample_mean / blank_mean, np.inf)
    keep = ratio >= min_fold
    return keep, {"removed": int((~keep).sum()), "min_fold": min_fold}


def normalize(
    matrix,
    method: str,
    roles: dict[str, str] | None = None,
    sample_names: list[str] | None = None,
):
    """行（サンプル）ごとにスケーリングして測定量の系統差を補正する。

    tic=行総和, median=行中央値, pqn=Probabilistic Quotient Normalization
    （参照は QC 中央値、無ければ全サンプル中央値）。none は恒等。
    """
    matrix = np.asarray(matrix, dtype=float)
    report = {"method": method}
    n = matrix.shape[0]
    if method == "none":
        return matrix, np.ones(n), report

    if method == "tic":
        factors = np.nansum(matrix, axis=1)
    elif method == "median":
        factors = np.nanmedian(matrix, axis=1)
    elif method == "pqn":
        ref_rows, used_qc = _reference_rows(matrix, roles, sample_names)
        reference = np.nanmedian(ref_rows, axis=0)  # 参照スペクトル
        safe_ref = np.where(reference == 0, np.nan, reference)
        quotients = matrix / safe_ref
        factors = np.nanmedian(quotients, axis=1)
        report["pqn_reference"] = "qc_median" if used_qc else "all_sample_median"
    else:
        raise ValueError(f"unknown normalization method: {method!r}")

    # 正規化係数が0/非有限のサンプル（検出特徴が疎で行中央値=0 等）を検出する。
    degenerate = (factors == 0) | ~np.isfinite(factors)
    # スケール後の全体水準を保つため、TIC/median は除数そのものを参照係数の
    # 中央値で正規化する（有効な係数のみで参照を計算する）。
    if method in ("tic", "median"):
        valid = factors[~degenerate]
        ref = np.nanmedian(valid) if valid.size else 1.0
        if not np.isfinite(ref) or ref == 0:
            ref = 1.0
        factors = factors / ref
    # 退化サンプルは行全体を NaN 化して破棄せず、未正規化のまま残す（factor=1.0）。
    # 旧実装は factor=NaN で行を全消去し、疎データで多数の試料を無言で失っていた。
    factors = np.where(degenerate, 1.0, factors)
    scaled = matrix / factors[:, None]
    n = matrix.shape[0]
    report["factors_finite"] = int((~degenerate).sum())
    report["n_samples"] = int(n)
    if degenerate.any():
        n_deg = int(degenerate.sum())
        report["unscaled_samples"] = n_deg
        report["caveat"] = (
            f"{n_deg}/{n} 試料は正規化係数が0または非有限（検出特徴が疎で行中央値=0 等）の"
            f"ため未正規化のまま残置しました（{method}）。該当試料の定量比較は測定量差を含み得ます。"
        )
    return scaled, factors, report


def _moving_median(values, window=5):
    """奇数窓の移動中央値（端は縮小窓）。numpy のみで外れ値に頑健な平滑化。

    LOESS（局所重み付き回帰）ではない。QC-RLSC の「QC で系統ドリフトを推定して
    補正する」枠組みを、平滑化器を移動中央値に置き換えて実装している。
    """
    n = len(values)
    half = window // 2
    smoothed = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        smoothed[i] = np.nanmedian(values[lo:hi])
    return smoothed


#: 相関を語るのに要る最小の点数。2 点なら必ず |r|=1 になり意味が無い。
_MIN_ORDERED_SAMPLES = 3


def run_order_correlation(components, sample_names, run_order):
    """主成分ごとに注入順との Pearson 相関を返す（語れないものは None）。

    QC が無いバッチでは `qc_drift_correct` も QC-RSD フィルタも動かない。
    その状態でも注入順さえあれば「その主成分が分析ドリフトを写しているか」は
    言える。**補正の代わりではなく、補正できないことの可視化**に使う。

    注入順を持たないサンプルは計算から外す（欠測を 0 で埋めない）。分散ゼロの
    主成分は相関が定義できないので None を返す —— nan は JSON にできない。
    """
    components = np.asarray(components, dtype=float)
    if components.ndim != 2:
        raise ValueError(f"components は 2 次元である必要があります: {components.shape}")

    orders = [run_order.get(n) for n in sample_names]
    keep = [i for i, o in enumerate(orders) if o is not None]
    n_pc = components.shape[1]
    if len(keep) < _MIN_ORDERED_SAMPLES:
        return [None] * n_pc

    x = np.array([float(orders[i]) for i in keep], dtype=float)
    if np.std(x) == 0:
        return [None] * n_pc

    out = []
    for pc in range(n_pc):
        y = components[keep, pc]
        if not np.all(np.isfinite(y)) or np.std(y) == 0:
            out.append(None)
            continue
        out.append(round(float(np.corrcoef(x, y)[0, 1]), 3))
    return out


def qc_drift_correct(matrix, roles, sample_names, run_order, min_qc=4, window=5):
    """QC を注入順に平滑化した系統ドリフトで、特徴量ごとに全サンプルを補正する。

    平滑化は QC 値の移動中央値（`_moving_median`）で、それを `np.interp` で全注入へ
    区分線形に補間し、`median(QC) / trend` を掛ける。全ランを 1 系列として扱い、
    バッチ別には補正しない。
    注入順が取れない、または QC が min_qc 未満なら未実施（skipped）。
    プールQC が層別と疑われる場合の警告は呼び出し側（各ツール）で付す。
    """
    matrix = np.asarray(matrix, dtype=float)
    orders = np.array([run_order.get(n) for n in sample_names], dtype=object)
    have_order = np.array([o is not None for o in orders])
    qc_mask = np.array([roles.get(n) == "qc" and run_order.get(n) is not None
                        for n in sample_names])
    if not have_order.all() or qc_mask.sum() < min_qc:
        return matrix, {
            "status": "skipped",
            "caveat": "注入順が欠落、または QC が不足のためドリフト補正は未実施。",
            "qc_used": int(qc_mask.sum()),
        }

    order_int = np.array([int(o) for o in orders])
    span = _qc_interspersion(order_int, qc_mask, roles, sample_names)
    if span["covered"] == 0:
        return matrix, {
            "status": "skipped",
            "caveat": (
                f"QC が試料列に挿入されていません（注入順 QC={span['qc_range']}, "
                f"試料={span['sample_range']}）。QC-RLSC は QC が試料の前後に散在することを"
                "前提とするため、この設計では補正が成立しません（内挿が全て外挿になる）。"),
            "qc_used": int(qc_mask.sum()),
            "qc_interspersion": span,
        }
    qc_order = order_int[qc_mask]
    sort = np.argsort(qc_order)
    qc_order_sorted = qc_order[sort]
    n_qc = qc_mask.sum()
    # 窓サイズをQCサンプル数以下に制限し、奇数を保証
    eff_window = min(window, int(n_qc))
    if eff_window % 2 == 0:
        eff_window = max(1, eff_window - 1)
    corrected = matrix.copy()
    for j in range(matrix.shape[1]):
        qc_vals = matrix[qc_mask, j][sort]
        trend = _moving_median(qc_vals, eff_window)
        global_level = np.nanmedian(qc_vals)
        if not np.isfinite(global_level) or global_level == 0:
            continue
        # 全サンプル注入順に対しトレンドを内挿し、補正係数=global_level/trend を適用
        interp_trend = np.interp(order_int, qc_order_sorted, trend)
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = np.where(interp_trend > 0, global_level / interp_trend, 1.0)
        corrected[:, j] = matrix[:, j] * factor
    report = {"status": "applied", "qc_used": int(qc_mask.sum()), "qc_interspersion": span}
    if span["covered"] < 0.5 * span["n_samples"]:
        report["caveat"] = (
            f"試料 {span['n_samples']} 件のうち QC 注入順区間に入るのは {span['covered']} 件のみ。"
            "区間外の試料は外挿補正となり、ドリフト補正の信頼性は限定的です。")
    return corrected, report


def _qc_interspersion(order_int, qc_mask, roles, sample_names) -> dict:
    """QC の注入順が試料の注入順を挟んでいるか（QC-RLSC の前提）を測る。

    QC を全試料の後ろにまとめて流す設計では np.interp が端値で頭打ちになり、
    「補正した」外見だけが残る。実際に内挿できた試料数を covered として返す。
    """
    sample_mask = np.array([roles.get(n) == "sample" for n in sample_names])
    qc_orders = order_int[qc_mask]
    sample_orders = order_int[sample_mask]
    if sample_orders.size == 0:
        return {"covered": 0, "n_samples": 0, "qc_range": None, "sample_range": None}
    lo, hi = int(qc_orders.min()), int(qc_orders.max())
    covered = int(((sample_orders >= lo) & (sample_orders <= hi)).sum())
    return {
        "covered": covered,
        "n_samples": int(sample_orders.size),
        "qc_range": [lo, hi],
        "sample_range": [int(sample_orders.min()), int(sample_orders.max())],
    }


def qc_rsd_filter(matrix, roles, sample_names, max_rsd=0.30):
    """QC 群での相対標準偏差 (SD/mean) が max_rsd を超える特徴量を除去する。"""
    matrix = np.asarray(matrix, dtype=float)
    n_features = matrix.shape[1]
    qc = _rows_for_role(matrix, roles, sample_names, "qc")
    if qc is None or qc.shape[0] < 2:
        return np.ones(n_features, dtype=bool), {
            "removed": 0, "status": "skipped",
            "caveat": "QC が無い/不足のため RSD フィルタは未実施。",
        }
    mean = np.nanmean(qc, axis=0)
    sd = np.nanstd(qc, axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rsd = np.where(mean > 0, sd / mean, np.inf)
    keep = rsd <= max_rsd
    return keep, {"removed": int((~keep).sum()), "max_rsd": max_rsd}


def detect_failed_qc(matrix, roles, sample_names, min_ratio=0.2) -> dict:
    """総強度が QC 中央値の min_ratio 未満に落ちた QC 注入（失敗注入）を検出する。

    失敗 QC を残したまま qc_rsd_filter を掛けると QC の RSD が全特徴で跳ね上がり、
    ほぼ全特徴が除去される（実測: 1345 → 51）。フィルタ側の閾値問題に見えるため、
    原因である QC 側を名指しで前景化する。行列は非破壊（検出のみ）。
    """
    matrix = np.asarray(matrix, dtype=float)
    qc_names = [n for n in sample_names if roles.get(n) == "qc"]
    if len(qc_names) < 3:
        return {"failed": [], "checked": len(qc_names)}
    index = {n: i for i, n in enumerate(sample_names)}
    tic = {n: float(np.nansum(matrix[index[n], :])) for n in qc_names}
    median = float(np.median(list(tic.values())))
    if not np.isfinite(median) or median <= 0:
        return {"failed": [], "checked": len(qc_names)}
    failed = sorted(n for n in qc_names if tic[n] < min_ratio * median)
    report = {
        "failed": failed,
        "checked": len(qc_names),
        "min_ratio": min_ratio,
        "qc_total_intensity": {n: tic[n] for n in qc_names},
    }
    if failed:
        report["caveat"] = (
            f"QC {len(failed)}/{len(qc_names)} 件の総強度が QC 中央値の {min_ratio:.0%} 未満です"
            f"（{', '.join(failed)}）。失敗注入の可能性が高く、残したまま max_qc_rsd を適用すると"
            "ほぼ全特徴が除去されます。arf_exclude で除外してから前処理し直してください。")
    return report


def impute(matrix, method="half_min"):
    """行列生成後に残る欠損 (NaN) を補完する。

    half_min=特徴量最小値の半分（既定）/ knn=sklearn KNNImputer /
    column_mean=列平均（現行互換）/ none=補完しない。
    """
    matrix = np.asarray(matrix, dtype=float)
    report = {"method": method, "missing": int(np.isnan(matrix).sum())}
    if method == "none" or report["missing"] == 0:
        return matrix, report
    out = matrix.copy()
    if method == "half_min":
        col_min = np.nanmin(np.where(np.isnan(out), np.inf, out), axis=0)
        col_min = np.where(np.isfinite(col_min), col_min, 0.0)
        fill = col_min / 2.0
        idx = np.where(np.isnan(out))
        out[idx] = np.take(fill, idx[1])
    elif method == "column_mean":
        col_mean = np.nanmean(out, axis=0)
        col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
        idx = np.where(np.isnan(out))
        out[idx] = np.take(col_mean, idx[1])
    elif method == "knn":
        from sklearn.impute import KNNImputer
        out = KNNImputer(n_neighbors=min(5, max(1, out.shape[0] - 1))).fit_transform(out)
    else:
        raise ValueError(f"unknown imputation method: {method!r}")
    return out, report


def drop_samples_by_role(matrix, sample_names, roles, drop_roles=("blank",)):
    """指定ロールの行を解析行列から外す。戻り値 ``(matrix, sample_names, dropped)``。

    ブランクは背景除去（blank_filter）の参照として使い終えたら外す。残したまま PCA に
    かけると、生体試料とは桁違いに低い総強度が PC1 を支配し、群分離の解釈が壊れる。
    一方 **QC は既定では外さない**——QC クラスタの締まり具合を PCA で見るのは分析の
    定番手段であり、除いてしまうとその確認手段を失うため。差次的解析のように群に
    混ざると困る側は、呼び出し側が drop_roles で明示する。

    ``dropped`` は ``{role: [sample_name, ...]}``（該当なしのロールは載らない）。
    どのロールも該当しなければ行列は素通しする。
    """
    matrix = np.asarray(matrix, dtype=float)
    targets = set(drop_roles)
    dropped: dict[str, list[str]] = {}
    keep_idx: list[int] = []
    kept_names: list[str] = []
    for i, name in enumerate(sample_names):
        role = roles.get(name, "sample")
        if role in targets:
            dropped.setdefault(role, []).append(name)
        else:
            keep_idx.append(i)
            kept_names.append(name)
    if not dropped:
        return matrix, list(sample_names), {}
    return matrix[keep_idx], kept_names, dropped


def detection_rates(detected_mask) -> np.ndarray:
    """特徴ごとの実検出率（gap-fill を除いた割合）。

    `detected_mask` は (特徴 × サンプル) の bool 行列。形式に依存しない純関数なので
    ここに置く。mzTab-M 経路は `.arf` 由来のマスクを、ARF 経路は自前の検出数を
    使うが、率の定義は 1 つに保つ。
    """
    mask = np.asarray(detected_mask, dtype=bool)
    if mask.ndim != 2 or mask.shape[1] == 0:
        return np.zeros(mask.shape[0] if mask.ndim == 2 else 0, dtype=float)
    return mask.sum(axis=1) / mask.shape[1]


def detection_rate_filter(detected_mask, min_rate: float):
    """検出率が `min_rate` 未満の特徴を落とす keep mask を返す。

    `detected_mask` は `detection_rates` と同じ (特徴 × サンプル)。率の定義を
    1 か所に保つため、閾値判定もここに置く——呼び出し側（解析行列の recipe
    filter、QC）が各々 `mask.mean(axis=...)` を書くと、軸を取り違えた側だけが
    静かに別の集合を落とす。

    gap-fill 後の非ゼロ値で検出を代用しない（それは検出率ではなく「値の有無」）。
    検出状態が不明なら呼び出し側が filter 自体を not_evaluable にする。
    """
    rates = detection_rates(detected_mask)
    keep = rates >= float(min_rate)
    return keep, {"method": "detection_rate", "min_detection_rate": float(min_rate),
                  "removed": int((~keep).sum()), "n_features": int(keep.size)}


def preprocess(matrix, sample_names, roles, run_order, recipe):
    """順序: ブランク除去 → 正規化 → ドリフト補正 → QC RSD フィルタ → 補完。

    列（特徴量）を落とすステップはマスクを蓄積し、最後にまとめて適用する。
    """
    matrix = np.asarray(matrix, dtype=float)
    n_features = matrix.shape[1]
    keep = np.ones(n_features, dtype=bool)
    applied: list[str] = []
    skipped: list[str] = []
    caveats: list[str] = []
    steps: dict = {}

    def _record(name: str, rep: dict) -> None:
        """要求されたステップを、実際に適用できたかで振り分ける。

        `recipe_applied` に載せるのは本当に走ったものだけ。skipped を混ぜると
        「やった」と「できなかった」が呼び出し側から区別できなくなる
        （QC 不在のバッチで drift_correct を要求した場合が実例）。
        """
        steps[name] = rep
        if rep.get("status") == "skipped":
            skipped.append(name)
        else:
            applied.append(name)
        if "caveat" in rep:
            caveats.append(rep["caveat"])

    # QC もブランクも無いバッチでは、前処理レシピの大半が原理的に効かない。
    # それでも preprocess は成功を返すので、要求の有無によらず先に明示する
    # （要求していないだけなのか、できないのかを呼び出し側が区別できるように）。
    # present_roles は `sample_names`（この関数が実際に評価する試料の集合）
    # を起点に roles.get() で引く——他の全消費者（_rows_for_role/
    # detect_failed_qc/_qc_interspersion 等）と同じ組み立て方に揃える。
    # `set(roles.values())` のように roles 辞書そのものを走査すると、
    # 呼び出し側（dataset_analysis.build_dataset_pp_inputs）が監査・表示用に
    # 維持している「include=false で除外済みの試料も含む全件分」の roles から
    # stale な role を拾ってしまい、実際には評価していない QC/blank の存在を
    # 誤認する（controller裁定 round2 のregression: include=false のQCが
    # matrix/sample_names からは正しく落ちているのに、この注意書きだけは
    # 「QC がある」と誤判定して出なくなっていた）。
    present_roles = {roles.get(n) for n in sample_names} if roles else set()
    if not present_roles & {"qc", "blank"}:
        caveats.append(
            "この試料群には QC もブランクも含まれていないため、"
            "ブランク背景除去・QC-RSD フィルタ・QC ドリフト補正は"
            "いずれも適用できません（要求の有無によらず）。"
            "注入順があれば、PCA の run_order_correlation でドリフトの有無だけは確認できます。")

    # 前処理の前に QC 自体の健全性を見る。失敗 QC は下流の全フィルタを汚染するため、
    # 正規化でスケールが動く前の生強度で判定する。
    qc_health = detect_failed_qc(matrix, roles, sample_names)
    steps["qc_health"] = qc_health
    if "caveat" in qc_health:
        caveats.append(qc_health["caveat"])

    if recipe.get("blank_min_fold") is not None:
        mask, rep = blank_filter(matrix, roles, sample_names, recipe["blank_min_fold"])
        keep &= mask
        _record("blank_filter", rep)

    method = recipe.get("normalize", "none")
    if method != "none":
        matrix, _, rep = normalize(matrix, method, roles, sample_names)
        _record("normalize", rep)

    if recipe.get("drift_correct"):
        matrix, rep = qc_drift_correct(matrix, roles, sample_names, run_order)
        _record("drift_correct", rep)

    if recipe.get("max_qc_rsd") is not None:
        mask, rep = qc_rsd_filter(matrix, roles, sample_names, recipe["max_qc_rsd"])
        keep &= mask
        _record("qc_rsd_filter", rep)

    features_after = int(keep.sum())
    # 特徴量フィルタが過度（閾値が厳しすぎる等）で解析不能になる場合を前景化する。
    # なお steps[*]["removed"] は各フィルタ独立のマスク件数で重複し得るため加算不可。
    # 実際に落ちた総数は features_removed_total（before-after）を参照すること。
    if n_features and features_after == 0:
        caveats.append(
            "すべての特徴量がフィルタで除去されました（残存0件）。閾値（max_qc_rsd 等）が"
            "厳しすぎる可能性があります。差次的解析/PCA は実行できません。")
    elif n_features and features_after < 0.1 * n_features:
        pct = 100.0 * (1.0 - features_after / n_features)
        caveats.append(
            f"特徴量の {pct:.0f}% が除去され {features_after}/{n_features} 件のみ残存しました。"
            "フィルタ閾値（max_qc_rsd/blank_min_fold）の見直しを検討してください。")

    matrix = matrix[:, keep]
    kept_idx = list(np.where(keep)[0])

    impute_method = recipe.get("impute", "half_min")
    matrix, rep = impute(matrix, impute_method)
    applied.append("impute")
    steps["impute"] = rep

    report = {
        "recipe_applied": applied,
        "recipe_skipped": skipped,
        "steps": steps,
        "caveats": caveats,
        "features_before": n_features,
        "features_after": features_after,
        "features_removed_total": n_features - features_after,
    }
    return matrix, kept_idx, report
