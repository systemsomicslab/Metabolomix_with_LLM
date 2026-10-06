"""④ 注釈付きの別スポット Y のイオンとしての説明（純関数）。spec §6。

X と Y が同時に溶出し（|ΔRT| ≤ rt_window）、Δm/z が既知の関係（アダクトの組・同位体・インソース断片）
で説明できるか、MS-DIAL の FoundInUpperMsMs リンクがあるときに候補にする。裏付けは MS-DIAL の
リンク（correl_similar / chrom_similar / found_in_upper_msms）か、試料間の log 強度の相関 r ≥ min_r。
M+1/M+2 は強度比が炭素数から期待される比の 1.5 倍以下のときだけ「同位体で説明できる」とし、
それを超えれば X は実在の別物質とみなして情報として添えるだけにする（Type II の重なり）。

FoundInUpperMsMs のリンクは無向の対（実データで相互 100%、相手が高 m/z 側になるのは半々）。
断片は前駆体より軽いので、リンクだけで関係にするのも、リンクを強い理由にするのも、
X が Y より軽い（X の m/z < Y の m/z）ときに限る。upper（前駆体）側は m/z の大小で決める。
リンクは裏付け（supported）には向きによらず数える。

deps: msdial.adducts、msdial.peak_verification（組成式）。
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from statistics import median

from metabolomix.msdial.adducts import mz_from_neutral, neutral_from_mz, parse_adduct
from metabolomix.msdial.peak_verification import parse_formula

C13_ABUNDANCE = 0.0107
ISOTOPE_RATIO_FACTOR = 1.5
SUPPORT_LINKS = frozenset({"correl_similar", "chrom_similar", "found_in_upper_msms"})
ISOTOPE_DELTAS = {"isotope_M+1": (1, 1.003355), "isotope_M+2": (2, 2.006710)}
#: 質量は各中性損失の組成式から monoisotopic_mass で計算した値（Task 5 で照合済み）。想定は spec §6 の表。
INSOURCE_LOSSES = [
    ("insource:-H2O", 18.010565),          # H2O。脱水（Cer・DG など）
    ("insource:-2H2O", 36.021129),         # H4O2。二重の脱水（Cer）
    ("insource:-NH3", 17.026549),          # NH3。[M+NH4]+ からの脱アンモニア
    ("insource:-HCOOCH3", 60.021129),      # C2H4O2。PC・SM の [M+HCOO]- → [M-CH3]-
    ("insource:-CH3COOCH3", 74.036779),    # C3H6O2。PC・SM の [M+CH3COO]- → [M-CH3]-
    ("insource:-C3H5NO2", 87.032028),      # C3H5NO2。PS の頭部基（セリン）脱離（neg）
    ("insource:-C2H8NO4P", 141.019094),    # C2H8NO4P。PE の頭部基脱離（pos）
    ("insource:-C3H8NO6P", 185.008923),    # C3H8NO6P。PS の頭部基脱離（pos）
    ("insource:-C6H10O5", 162.052824),     # C6H10O5。ヘキソース脱離（HexCer・PI など）
]


def expected_isotope_ratio(n_carbon: int, k: int) -> float:
    q = C13_ABUNDANCE / (1 - C13_ABUNDANCE)
    return math.comb(n_carbon, k) * q ** k


def _carbon(formula) -> int | None:
    try:
        return parse_formula(formula).get("C") if formula else None
    except ValueError:
        return None


def _pairs(x, y):
    return [(a, b) for a, b in zip(x, y) if a is not None and b is not None and a > 0 and b > 0]


def profile_correlation(x, y, min_n: int = 5) -> float | None:
    pairs = _pairs(x, y)
    if len(pairs) < min_n:
        return None
    lx = [math.log1p(a) for a, _ in pairs]
    ly = [math.log1p(b) for _, b in pairs]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    sx = math.sqrt(sum((v - mx) ** 2 for v in lx))
    sy = math.sqrt(sum((v - my) ** 2 for v in ly))
    if sx == 0 or sy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / (sx * sy)


def median_ratio(x, y) -> float | None:
    pairs = _pairs(x, y)
    return median(a / b for a, b in pairs) if pairs else None


def _link_kinds(links: dict, a: int, b: int) -> list[str]:
    kinds = {l["kind"] for l in links.get(a, []) if l["spot_id"] == b}
    kinds |= {l["kind"] for l in links.get(b, []) if l["spot_id"] == a}
    return sorted(kinds)


def _explanations(x: dict, y: dict, adducts: list[str], mz_tol: float) -> list[tuple[str, float]]:
    """(relation, Δm/z [Da]) の一覧。"""
    found = []
    y_adduct = parse_adduct(y.get("adduct"))
    if y_adduct is not None:
        neutral = neutral_from_mz(y["mz"], y_adduct)
        for name in adducts:
            other = parse_adduct(name)
            if other is None or other.name == y_adduct.name:
                continue
            delta = x["mz"] - mz_from_neutral(neutral, other)
            if abs(delta) <= mz_tol:
                found.append((f"adduct:{other.name}/{y_adduct.name}", delta))
    charge = y_adduct.charge if y_adduct else 1
    for relation, (_, shift) in ISOTOPE_DELTAS.items():
        delta = x["mz"] - (y["mz"] + shift / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    for relation, loss in INSOURCE_LOSSES:
        delta = x["mz"] - (y["mz"] - loss / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    return found


def find_relations(target: dict, partners: list[dict], *, adducts, rt_window, mz_tol, min_r, links) -> list[dict]:
    ordered = sorted(partners, key=lambda p: p["rt"])
    rts = [p["rt"] for p in ordered]
    lo, hi = bisect_left(rts, target["rt"] - rt_window), bisect_right(rts, target["rt"] + rt_window)
    out = []
    for y in ordered[lo:hi]:
        if y["spot_id"] == target["spot_id"]:
            continue
        kinds = _link_kinds(links, target["spot_id"], y["spot_id"])
        r = profile_correlation(target.get("heights") or [], y.get("heights") or [])
        supported = bool(SUPPORT_LINKS & set(kinds)) or (r is not None and r >= min_r)
        # 無向のリンクの向きは m/z で決める: X が断片（軽い側）、Y が前駆体（重い側）のときだけ。
        y_is_precursor = target["mz"] < y["mz"]
        upper_link = "found_in_upper_msms" in kinds and y_is_precursor
        explained = _explanations(target, y, adducts, mz_tol)
        if not explained and upper_link:
            explained = [("found_in_upper_msms", None)]
        for relation, delta in explained:
            entry = {"kind": "ion", "of": y["spot_id"], "of_name": y.get("name"), "relation": relation,
                     "dmz_mda": None if delta is None else round(delta * 1000, 2),
                     "drt": round(target["rt"] - y["rt"], 3), "r": None if r is None else round(r, 3),
                     "msdial_links": kinds, "isotope_ratio": None, "expected_ratio": None,
                     "isotope_consistent": None, "informational": False}
            if relation in ISOTOPE_DELTAS:
                k = ISOTOPE_DELTAS[relation][0]
                carbon = _carbon(y.get("formula"))
                ratio = median_ratio(target.get("heights") or [], y.get("heights") or [])
                expected = expected_isotope_ratio(carbon, k) if carbon else None
                entry.update(isotope_ratio=None if ratio is None else round(ratio, 4),
                             expected_ratio=None if expected is None else round(expected, 4))
                if ratio is not None and expected is not None:
                    entry["isotope_consistent"] = ratio <= expected * ISOTOPE_RATIO_FACTOR
                    entry["informational"] = not entry["isotope_consistent"]
            entry["strong"] = bool(supported and not entry["informational"] and (
                relation == "found_in_upper_msms" or upper_link
                or entry["isotope_consistent"] is True))
            entry["soft"] = [] if supported else ["no_support"]
            out.append(entry)
    out.sort(key=lambda e: (not e["strong"], e["informational"], bool(e["soft"]),
                            abs(e["dmz_mda"]) if e["dmz_mda"] is not None else 0.0))
    for index, entry in enumerate(out, 1):
        entry["candidate_id"] = f"R{index}"
    return out
