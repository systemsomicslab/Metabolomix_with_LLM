"""RT–m/z 傾向（純関数）。判断材料の 1 つにとどめる（spec §4 ⑤）。

判定はクラスごとに 1 本の加法モデル RT = a + b·C + c·DB を Huber の IRLS で頑健回帰し、
残差を頑健尺度（1.4826·MAD）で割った z で外れを見る。不飽和度が上がると線形性が落ちる
ので、群（クラス×DB）ごとの RT–m/z 直線の R² と点数も返し、ビューアが当てはまりの
悪い群を弱めて表示できるようにする。
"""
from __future__ import annotations

import numpy as np

from lipidmix.msdial.lipid_identity import _clean_msdial_name

# 線形にきれいに揃ったクラスでは Huber スケール（1.4826·MAD）がほぼ 0 に潰れ、
# 素の残差/スケールでは全点が外れ値（scale=0 で除算不可）にも無外れ値（brief 案の
# else 0.0）にもなり得る。分単位の実務上無意味な微小スケールを下駄で底上げし、
# 「事実上完全に線形」なクラスでも 1 点だけの外れを検出できるようにする（Ruling 2）。
MIN_SCALE = 0.02


def composition(name) -> tuple[int, int] | None:
    if not name or not str(name).strip():
        return None
    clean, _ = _clean_msdial_name(name)
    if not clean:
        return None
    try:
        from pygoslin.parser.Parser import LipidParser
        lipid = LipidParser().parse(clean)
        info = lipid.lipid.info
        return int(info.num_carbon), int(info.double_bonds)
    except Exception:  # noqa: BLE001 - 脂質名として読めない注釈は傾向の対象外
        return None


def _huber(X: np.ndarray, y: np.ndarray, iterations: int = 30):
    weights = np.ones(len(y))
    beta = np.zeros(X.shape[1])
    scale = 0.0
    for _ in range(iterations):
        sw = np.sqrt(weights)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        residual = y - X @ beta
        mad = np.median(np.abs(residual - np.median(residual)))
        scale = 1.4826 * mad
        if scale <= 1e-9:
            break
        k = 1.345 * scale
        weights = np.where(np.abs(residual) <= k, 1.0, k / np.maximum(np.abs(residual), 1e-12))
    return beta, y - X @ beta, scale


def _r2(y, residual) -> float | None:
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return round(1.0 - float(np.sum(residual ** 2)) / ss_tot, 4) if ss_tot > 0 else None


def fit_trends(points: list[dict], th: dict) -> dict:
    classes: dict = {}
    groups: dict = {}
    spots: dict = {}
    by_class: dict[str, list[dict]] = {}
    for p in points:
        by_class.setdefault(p["ontology"] or "", []).append(p)

    for ontology, members in by_class.items():
        if len(members) < th["trend_min_points"] or len({p["carbon"] for p in members}) < 2:
            continue
        y = np.array([p["rt"] for p in members], dtype=float)
        columns = [np.ones(len(members)), np.array([p["carbon"] for p in members], dtype=float)]
        if len({p["db"] for p in members}) >= 2:
            columns.append(np.array([p["db"] for p in members], dtype=float))
        X = np.column_stack(columns)
        beta, residual, scale = _huber(X, y)
        effective_scale = max(scale, MIN_SCALE)
        z_all = residual / effective_scale
        outlier_all = np.abs(z_all) > th["trend_outlier_z"]

        inlier_mask = ~outlier_all
        r2 = _r2(y[inlier_mask], residual[inlier_mask]) if inlier_mask.any() else None
        reliable = r2 is not None and r2 >= th["trend_min_r2"]
        classes[ontology] = {"n": len(members), "r2": r2,
                             "coef": [round(float(b), 5) for b in beta],
                             "scale": round(float(effective_scale), 5)}
        for p, r, z, outlier in zip(members, residual, z_all, outlier_all):
            spots[p["spot_id"]] = {"residual": round(float(r), 4), "z": round(float(z), 3),
                                   "outlier": bool(outlier),
                                   "reliable": reliable}

        per_db: dict[int, dict] = {}
        for db in sorted({p["db"] for p in members}):
            subset = [p for p in members if p["db"] == db]
            entry = {"n": len(subset), "slope": None, "intercept": None, "r2": None}
            if len(subset) >= 3 and len({p["mz"] for p in subset}) >= 2:
                xs = np.array([p["mz"] for p in subset], dtype=float)
                ys = np.array([p["rt"] for p in subset], dtype=float)
                slope, intercept = np.polyfit(xs, ys, 1)
                entry.update(slope=round(float(slope), 6), intercept=round(float(intercept), 4),
                             r2=_r2(ys, ys - (slope * xs + intercept)))
            per_db[db] = entry
        groups[ontology] = per_db
    return {"classes": classes, "groups": groups, "spots": spots}
