"""MS1 の同位体パターン（純関数）。ビューアの MS1 パネル用で、判定には使わない。

実測は `.arf2` の Key 53 `IsotopicPeaks`（代表試料の M, M+1, M+2。要素の Key 0 が
相対強度 %、Key 1 が m/z、Key 4 が絶対強度。docs/schema/AlignmentSpotProperty.md・molecule_ms_reference.md）。
理論は組成式＋アダクトの原子から、整数質量の分解能（M+k の k ごと）で畳み込む。
アダクトが読めなければ中性の組成式だけで計算し、`basis` で区別する。
deps: msdial.adducts / msdial.peak_verification。
"""
from __future__ import annotations

from metabolomix.msdial.adducts import adduct_composition
from metabolomix.msdial.peak_verification import parse_formula

#: 元素ごとの天然存在比（整数質量のずれ k → 存在比）。IUPAC 代表値。
ABUNDANCES: dict[str, dict[int, float]] = {
    "H": {0: 0.999885, 1: 0.000115},
    "C": {0: 0.9893, 1: 0.0107},
    "N": {0: 0.99636, 1: 0.00364},
    "O": {0: 0.99757, 1: 0.00038, 2: 0.00205},
    "P": {0: 1.0},
    "S": {0: 0.9499, 1: 0.0075, 2: 0.0425, 4: 0.0001},
    "Na": {0: 1.0},
    "K": {0: 0.932581, 1: 0.000117, 2: 0.067302},
    "Cl": {0: 0.7576, 2: 0.2424},
    "F": {0: 1.0},
    "Br": {0: 0.5069, 2: 0.4931},
    "D": {0: 1.0},
}
N_PEAKS = 3          # M, M+1, M+2（MS-DIAL の IsotopicPeaks と同じ本数）


def _convolve(a: list[float], b: dict[int, float], n: int) -> list[float]:
    out = [0.0] * n
    for i, x in enumerate(a):
        if not x:
            continue
        for k, p in b.items():
            if i + k < n:
                out[i + k] += x * p
    return out


def _power(dist: dict[int, float], count: int, n: int) -> list[float]:
    """元素 1 種 `count` 個ぶんの分布（二乗法）。"""
    result = [1.0] + [0.0] * (n - 1)
    base = [dist.get(k, 0.0) for k in range(n)]
    while count:
        if count & 1:
            result = _convolve(result, dict(enumerate(base)), n)
        count >>= 1
        if count:
            base = _convolve(base, dict(enumerate(base)), n)
    return result


def theoretical_envelope(formula, adduct, n: int = N_PEAKS) -> dict | None:
    """M=100 とした M..M+(n-1) の相対強度。組成式が読めなければ None。"""
    if not formula or str(formula).strip() in ("", "Unknown"):
        return None
    try:
        counts = parse_formula(str(formula))
    except ValueError:
        return None
    composition = adduct_composition(adduct)
    basis = "formula"
    if composition is not None:
        n_mol, _charge, delta = composition
        merged = {k: v * n_mol for k, v in counts.items()}
        for k, v in delta.items():
            merged[k] = merged.get(k, 0) + v
        if all(v >= 0 for v in merged.values()):
            counts, basis = merged, "formula+adduct"
    if any(element not in ABUNDANCES for element in counts):
        return None
    dist = [1.0] + [0.0] * (n - 1)
    for element, count in counts.items():
        if count:
            dist = _convolve(dist, dict(enumerate(_power(ABUNDANCES[element], count, n))), n)
    if dist[0] <= 0:
        return None
    return {"relative": [round(x / dist[0] * 100, 2) for x in dist], "basis": basis}


def measured_envelope(raw) -> list[list]:
    """`.arf2` Key 53 の生配列から [[m/z, 相対強度 %, 絶対強度], ...]。絶対強度は要素の Key 4
    `AbsoluteAbundance`（整数に丸める。無い・読めなければ None）。m/z と相対強度が読めない要素が
    あれば全体を捨てる。"""
    if not isinstance(raw, list):
        return []
    peaks = []
    for item in raw:
        try:
            peaks.append([round(float(item[1]), 4), round(float(item[0]), 2), _absolute(item)])
        except (TypeError, ValueError, IndexError):
            return []
    return peaks


def _absolute(item) -> int | None:
    try:
        return round(float(item[4]))
    except (TypeError, ValueError, IndexError, OverflowError):
        return None
