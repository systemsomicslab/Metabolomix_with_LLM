"""④ 注釈付きの別スポット Y のイオンとしての説明（純関数）。spec §6。

X と Y が同時に溶出し（|ΔRT| ≤ rt_window）、Δm/z が既知の関係（アダクトの組・同位体・インソース断片）
で説明できるか、MS-DIAL の FoundInUpperMsMs リンクがあるときに候補にする。裏付けは MS-DIAL の
リンク（correl_similar / chrom_similar / found_in_upper_msms）か、試料間の log 強度の相関 r ≥ min_r。
同位体は Y の組成式から期待される強度比の 1.5 倍以下のときだけ「同位体で説明できる」とし、
それを超えれば X は実在の別物質とみなして情報として添えるだけにする（Type II の重なり）。

強い説明（`strong`）は、(a) 前駆体側の FoundInUpperMsMs リンク（裏付けつき）、(b) 同位体で強度比が
期待比に合う（裏付けつき・M+2 の取り違えの恐れが無い）、(c) アダクトの組・インソース断片で、Δm/z と
同時溶出に加えて試料間の相関 r ≥ min_r がそろう、のどれか。Δm/z が合うだけでは断定しない。

M+2 は 13C2（+2.006710）・37Cl（+1.997050）・81Br（+1.997952）の 3 つを見る。Y の組成式が分かっていて
その元素（Cl・Br）を持たなければ、そのハロゲンの説明は出さない。2 つ以上の M+2 が同じ Δm/z の窓に
入るときは一方に決めず、どれにも `ambiguous` を付けて並べ（強い説明にしない）、すべての期待比に
合わないときだけ情報（別物質）にする——13C2 の期待比に合わないことだけで「別物質」と言わない。

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
#: 重い同位体の存在比（IUPAC の代表値）。期待比は C(n,k)·(a/(1−a))^k。
ISOTOPE_ABUNDANCE = {"C": C13_ABUNDANCE, "Cl": 0.2424, "Br": 0.4931}
ISOTOPE_RATIO_FACTOR = 1.5
SUPPORT_LINKS = frozenset({"correl_similar", "chrom_similar", "found_in_upper_msms"})
#: 関係コード → (置き換わる重い同位体の数 k, 質量差 [Da], 元素)。質量差は同位体の質量の差から計算した値。
ISOTOPE_DELTAS = {
    "isotope_M+1": (1, 1.003355, "C"),
    "isotope_M+2": (2, 2.006710, "C"),
    "isotope_M+2_37Cl": (1, 1.997050, "Cl"),     # 37Cl − 35Cl
    "isotope_M+2_81Br": (1, 1.997952, "Br"),     # 81Br − 79Br
}
#: M+2 の同位体（Δm/z の窓が重なりうる。既定の 10 mDa では 13C2 と 37Cl の差 9.66 mDa が窓の内側）。
M2_ISOTOPES = frozenset({"isotope_M+2", "isotope_M+2_37Cl", "isotope_M+2_81Br"})
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


def expected_isotope_ratio(n_atoms: int, k: int, abundance: float = C13_ABUNDANCE) -> float:
    """n 個の原子のうち k 個が重い同位体に置き換わったイオンの、単一同位体イオンに対する強度比。"""
    q = abundance / (1 - abundance)
    return math.comb(n_atoms, k) * q ** k


def _element_counts(formula) -> dict | None:
    """組成式の元素の数。組成式が無い・読めなければ None（「その元素を持たない」とは区別する）。"""
    try:
        return parse_formula(formula) if formula else None
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
    counts = _element_counts(y.get("formula"))
    for relation, (_, shift, element) in ISOTOPE_DELTAS.items():
        if element != "C" and counts is not None and not counts.get(element):
            continue                    # Y の組成式にその元素が無い（その同位体では説明できない）
        delta = x["mz"] - (y["mz"] + shift / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    for relation, loss in INSOURCE_LOSSES:
        delta = x["mz"] - (y["mz"] - loss / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    return found


def _isotope_check(entry: dict, relation: str, x: dict, y: dict) -> None:
    """同位体の行に強度比・期待比・`isotope_consistent` を書く（期待比が求まらなければ consistent は None）。"""
    k, _, element = ISOTOPE_DELTAS[relation]
    n_atoms = (_element_counts(y.get("formula")) or {}).get(element)
    ratio = median_ratio(x.get("heights") or [], y.get("heights") or [])
    expected = (expected_isotope_ratio(n_atoms, k, ISOTOPE_ABUNDANCE[element])
                if n_atoms and n_atoms >= k else None)
    entry.update(isotope_ratio=None if ratio is None else round(ratio, 4),
                 expected_ratio=None if expected is None else round(expected, 4))
    if ratio is not None and expected is not None:
        entry["isotope_consistent"] = ratio <= expected * ISOTOPE_RATIO_FACTOR


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
        correlated = r is not None and r >= min_r
        supported = bool(SUPPORT_LINKS & set(kinds)) or correlated
        # 無向のリンクの向きは m/z で決める: X が断片（軽い側）、Y が前駆体（重い側）のときだけ。
        y_is_precursor = target["mz"] < y["mz"]
        upper_link = "found_in_upper_msms" in kinds and y_is_precursor
        explained = _explanations(target, y, adducts, mz_tol)
        if not explained and upper_link:
            explained = [("found_in_upper_msms", None)]
        entries = []
        for relation, delta in explained:
            entry = {"kind": "ion", "of": y["spot_id"], "of_name": y.get("name"), "relation": relation,
                     "dmz_mda": None if delta is None else round(delta * 1000, 2),
                     "drt": round(target["rt"] - y["rt"], 3), "r": None if r is None else round(r, 3),
                     "msdial_links": kinds, "isotope_ratio": None, "expected_ratio": None,
                     "isotope_consistent": None, "informational": False, "ambiguous": False}
            if relation in ISOTOPE_DELTAS:
                _isotope_check(entry, relation, target, y)
                entry["informational"] = entry["isotope_consistent"] is False
            entries.append(entry)
        m2 = [e for e in entries if e["relation"] in M2_ISOTOPES]
        if len(m2) > 1:
            # 窓が重なる M+2 は一方に決めない。どの期待比にも合わないときだけ別物質（情報）とする。
            all_inconsistent = all(e["isotope_consistent"] is False for e in m2)
            for e in m2:
                e["ambiguous"] = True
                e["informational"] = all_inconsistent
        for entry in entries:
            relation = entry["relation"]
            # (a) 前駆体側の FoundInUpperMsMs（裏付けつき）
            by_upper = supported and (relation == "found_in_upper_msms" or upper_link)
            if relation in ISOTOPE_DELTAS:
                # (b) 同位体: 裏付けつきで期待比に合い、M+2 の取り違えの恐れが無い
                by_own = supported and entry["isotope_consistent"] is True and not entry["ambiguous"]
            elif relation == "found_in_upper_msms":
                by_own = False
            else:
                # (c) アダクトの組・インソース断片: Δm/z と同時溶出（rt_window）に加えて相関 r ≥ min_r
                by_own = correlated
            entry["strong"] = bool(not entry["informational"] and (by_upper or by_own))
            entry["soft"] = [] if supported else ["no_support"]
            out.append(entry)
    out.sort(key=lambda e: (not e["strong"], e["informational"], bool(e["soft"]),
                            abs(e["dmz_mda"]) if e["dmz_mda"] is not None else 0.0))
    for index, entry in enumerate(out, 1):
        entry["candidate_id"] = f"R{index}"
    return out
