"""EIC の形状指標（純関数）。画像ではなく系列の数値に対して計算する。

最適化を伴うガウス当てはめは収束の失敗が判定を揺らすので使わない。頂点と半値幅から
決まる理想ガウスとの R² を「ガウスらしさ」とする（決定的で、同じ入力に同じ値）。
"""
from __future__ import annotations

import math
import statistics

_FWHM_TO_SIGMA = 2.354820045


def _half_height_x(xs, ys, apex, half, step):
    i = apex
    while 0 <= i + step < len(xs) and ys[i + step] > half:
        i += step
    j = i + step
    if not 0 <= j < len(xs):
        return xs[i]
    y0, y1 = ys[i], ys[j]
    if y0 == y1:
        return xs[j]
    return xs[i] + (xs[j] - xs[i]) * (y0 - half) / (y0 - y1)


def sample_shape(chromatogram, *, left, top, right) -> dict:
    points = [(float(x), float(y)) for x, y in chromatogram if left <= x <= right]
    if not points:
        return {"n_points": 0, "apex_rt": None, "apex_in_window": False,
                "gauss_r2": None, "n_maxima": 0, "symmetry": None}
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    apex = max(range(len(ys)), key=ys.__getitem__)
    height = ys[apex]
    apex_in_window = 0 < apex < len(ys) - 1

    n_maxima = sum(
        1 for i in range(1, len(ys) - 1)
        if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1] and ys[i] >= 0.1 * height)
    if apex_in_window and n_maxima == 0:
        n_maxima = 1

    gauss_r2 = None
    symmetry = None
    if height > 0 and len(ys) >= 3:
        half = height / 2.0
        left_x = _half_height_x(xs, ys, apex, half, -1)
        right_x = _half_height_x(xs, ys, apex, half, +1)
        sigma = (right_x - left_x) / _FWHM_TO_SIGMA
        if sigma > 0:
            model = [height * math.exp(-((x - xs[apex]) ** 2) / (2 * sigma ** 2)) for x in xs]
            mean = sum(ys) / len(ys)
            ss_tot = sum((y - mean) ** 2 for y in ys)
            ss_res = sum((y - m) ** 2 for y, m in zip(ys, model))
            gauss_r2 = round(1.0 - ss_res / ss_tot, 4) if ss_tot > 0 else None
        lw, rw = xs[apex] - left_x, right_x - xs[apex]
        symmetry = round(lw / rw, 3) if lw > 0 and rw > 0 else None

    return {"n_points": len(points), "apex_rt": round(xs[apex], 4),
            "apex_in_window": apex_in_window, "gauss_r2": gauss_r2,
            "n_maxima": n_maxima, "symmetry": symmetry}


def sample_is_good(shape: dict, th: dict) -> bool:
    return (shape["apex_in_window"]
            and shape["n_points"] >= th["eic_min_points"]
            and shape["gauss_r2"] is not None and shape["gauss_r2"] >= th["eic_min_r2"]
            and shape["n_maxima"] <= th["eic_max_maxima"])


def spot_shape(samples: list[dict], th: dict) -> dict:
    per_sample = {}
    detected = []
    for sample in samples:
        shape = sample_shape(sample["chromatogram"], left=sample["left"],
                             top=sample["top"], right=sample["right"])
        per_sample[sample["file_id"]] = shape
        if sample["detected"]:
            detected.append((sample, shape))

    n_detected = len(detected)
    n_good = sum(1 for _, shape in detected if sample_is_good(shape, th))
    fraction = round(n_good / n_detected, 3) if n_detected else None
    tops = [s["top"] for s, _ in detected]
    apex_rt_sd = round(statistics.pstdev(tops), 4) if len(tops) >= 2 else None

    if fraction is None:
        band = "UNKNOWN"
    elif fraction >= th["eic_pass_frac"]:
        band = "PASS"
    elif fraction >= th["eic_borderline_frac"]:
        band = "BORDERLINE"
    else:
        band = "FAIL"
    return {"n_detected": n_detected, "n_good": n_good, "good_fraction": fraction,
            "apex_rt_sd": apex_rt_sd, "per_sample": per_sample, "band": band}
