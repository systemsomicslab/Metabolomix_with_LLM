"""MS1 の同位体パターン（純関数）。ビューアの MS1 パネル（`theoretical_envelope`）と、イオン関係の
同位体の期待比（`curation.relations`。`isotopologues`）が使う。

実測は `.arf2` の Key 53 `IsotopicPeaks`（代表試料の M, M+1, M+2。要素の Key 0 が
相対強度 %、Key 1 が m/z、Key 4 が絶対強度。docs/schema/AlignmentSpotProperty.md・molecule_ms_reference.md）。
理論は組成式＋アダクトの原子（`msdial.adducts.ion_counts`: n_mol 倍＋付加・脱離する原子）から、
整数質量の分解能（M+k の k ごと）で畳み込む。M+2 は 13C2・18O・34S・37Cl・81Br・41K と交差項（13C+15N など）の和。
アダクトが読めなければ中性の組成式だけで計算し、`basis` で区別する。
`isotopologues` は同じ M+1・M+2 を同位体の組み合わせごとに、正確な質量差つきで分けて返す（k ごとの和は
`envelope_ratios` と同じ）。イオン関係は、そのうち X の m/z の窓に入るものだけを足す。
deps: msdial.adducts / msdial.peak_verification。
"""
from __future__ import annotations

from metabolomix.msdial.adducts import ion_counts
from metabolomix.msdial.peak_verification import parse_formula, strip_charge

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
    # 単一同位体の元素（同位体のピークを作らない）。甲状腺ホルモン・ヨウ素造影剤、ヒ素化合物、コバラミン
    "I": {0: 1.0},
    "As": {0: 1.0},
    "Co": {0: 1.0},
}
#: 重い同位体の、最も軽い同位体からの質量差 [Da]（整数質量のずれ k → 質量差）。NIST の原子質量
#: （Atomic Weights and Isotopic Compositions）の差。`ABUNDANCES` の k ≥ 1 の同位体ごとに 1 つ。
ISOTOPE_SHIFTS: dict[str, dict[int, float]] = {
    "H": {1: 1.00627674589},                                     # 2H − 1H
    "C": {1: 1.00335483507},                                     # 13C − 12C
    "N": {1: 0.99703489445},                                     # 15N − 14N
    "O": {1: 1.00421713693, 2: 2.00424499329},                   # 17O・18O − 16O
    "S": {1: 0.9993877354, 2: 1.9957958296, 4: 3.9950095356},    # 33S・34S・36S − 32S
    "K": {1: 1.0002916796, 2: 1.9981187715},                     # 40K・41K − 39K
    "Cl": {2: 1.99704992},                                       # 37Cl − 35Cl
    "Br": {2: 1.9979521},                                        # 81Br − 79Br
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


def envelope_ratios(counts: dict[str, int], n: int = N_PEAKS) -> list[float] | None:
    """元素の数（イオンの組成）から、M..M+(n-1) の強度の M に対する比（M が 1。丸めない）。

    表（`ABUNDANCES`）に無い元素・負の個数があれば None。
    """
    if any(element not in ABUNDANCES or count < 0 for element, count in counts.items()):
        return None
    dist = [1.0] + [0.0] * (n - 1)
    for element, count in counts.items():
        if count:
            dist = _convolve(dist, dict(enumerate(_power(ABUNDANCES[element], count, n))), n)
    if dist[0] <= 0:
        return None
    return [x / dist[0] for x in dist]


def isotopologues(counts: dict[str, int]) -> list[tuple[int, float, float]] | None:
    """元素の数（イオンの組成）から、M+1 と M+2 の同位体の組み合わせごとの (k, 質量差 [Da], M に対する比)。

    M+1 は k = 1 の重い同位体 1 個（13C・15N・2H・17O・33S・40K）、M+2 は k = 2 の重い同位体 1 個
    （18O・34S・37Cl・81Br・41K）か k = 1 の重い同位体 2 個（13C2・13C+15N など）。比は q = 重い同位体の存在比 /
    最も軽い同位体の存在比として、1 個なら n·q、同じ元素の 2 個なら C(n,2)·q²、別の元素の 2 個なら
    n_a·q_a·n_b·q_b（多項分布の比で、近似ではない）。k ごとの和は `envelope_ratios` の M+k と一致する。
    表（`ABUNDANCES`）に無い元素・負の個数があれば None。
    """
    if any(element not in ABUNDANCES or count < 0 for element, count in counts.items()):
        return None
    singles = []
    for element, count in counts.items():
        if not count:
            continue
        light = ABUNDANCES[element][0]
        for k, abundance in ABUNDANCES[element].items():
            if k in (1, 2) and abundance:
                singles.append((element, k, ISOTOPE_SHIFTS[element][k], count * abundance / light))
    terms = [(k, shift, ratio) for _, k, shift, ratio in singles]
    ones = [single for single in singles if single[1] == 1]
    for index, (element, _, shift, ratio) in enumerate(ones):
        count = counts[element]
        if count >= 2:              # 同じ元素の 2 個: C(n,2)·q² = (n·q)²·(n−1)/(2n)
            terms.append((2, shift + shift, ratio * ratio * (count - 1) / (2 * count)))
        for _, _, other_shift, other_ratio in ones[index + 1:]:
            terms.append((2, shift + other_shift, ratio * other_ratio))
    return terms


def theoretical_envelope(formula, adduct, n: int = N_PEAKS) -> dict | None:
    """M=100 とした M..M+(n-1) の相対強度。組成式が読めなければ None。"""
    if not formula or str(formula).strip() in ("", "Unknown"):
        return None
    try:
        counts = parse_formula(strip_charge(str(formula)))
    except ValueError:
        return None
    basis = "formula"
    ion = ion_counts(counts, adduct)
    if ion is not None:
        counts, basis = ion, "formula+adduct"
    ratios = envelope_ratios(counts, n)
    if ratios is None:
        return None
    return {"relative": [round(x * 100, 2) for x in ratios], "basis": basis}


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
