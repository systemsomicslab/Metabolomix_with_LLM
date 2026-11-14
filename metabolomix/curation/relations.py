"""④ 注釈付きの別スポット Y のイオンとしての説明（純関数）。spec §6。

X と Y が同時に溶出し（|ΔRT| ≤ rt_window）、Δm/z が既知の関係（アダクトの組・同位体・インソース断片）
で説明できるか、MS-DIAL の FoundInUpperMsMs リンクがあるときに候補にする。裏付けは MS-DIAL の
リンク（correl_similar / chrom_similar / found_in_upper_msms）か、試料間の log 強度の相関 r ≥ min_r。
同位体は Y のイオンの組成から期待される強度比の 1.5 倍以下のときだけ「同位体で説明できる」とし、
それを超えれば X は実在の別物質とみなして情報として添えるだけにする（Type II の重なり）。

アダクトの組は Y の中性質量に `adducts` の各アダクトを付けた m/z で探す。Y 自身の付加体は綴りでなく組成で
見分けて除く（`msdial.adducts.same_adduct`。mzTab-M の `[M+H]1+` と `[M+H]+`、`[M+HCOO]-` と `[M+FA-H]-` は
同じ）。Y と同じ m/z になる付加体（`[M+H]+` に対する `[2M+2H]2+` など）も、Δm/z が何も区別しないので除く。

同位体の期待比は Y の**イオン**の組成（中性の組成式の n_mol 倍に付加体の原子を足したもの。`[2M+H]+` は炭素も
2 倍、`[M+Cl]-` は Cl を 1 個持つ）の M+k の同位体の組み合わせ（`curation.isotope.isotopologues`。M+1 は 13C・15N・
2H・17O・33S・40K、M+2 は 13C2・18O・34S・37Cl・81Br・41K と 13C+15N などの交差項。それぞれ正確な質量差を持つ）の
うち、**m/z が X から `mz_tol` 以内に来るものの比の和**。X のピークに入りうるのはその窓の中の同位体だけなので、
分離された 1 つの同位体（13C2 だけ・34S だけ）はその位置の近くの同位体だけと比べる。分離された位置の強度を窓の外の
同位体（Cl を持つイオンの 13C2 に 37Cl など）まで足した和と比べると、その位置の同位体の何倍もの別物質を「同位体で
説明できる」と言ってしまう。分離できない M+2 の束（重心が 13C2 と 34S の間に来る）が束全体の和と比べられるのは、
X から束のどの同位体までの差も `mz_tol` 以内のときだけ——重心が真ん中なら 13C2 と 34S の差 10.91 mDa の半分（約 5.5 mDa
を電荷で割った値）以上の窓が要る（既定の 10 mDa は覆う）。それより狭い窓では、重心は近い方の位置の行だけに当たり、
期待比は遠い方の同位体を欠くので、Y 自身の束が `informational`（別物質）になりうる。重心から両端のどちらまでの差より
狭い窓（タウロコール酸 [M-H]- で約 5.2 mDa 未満）では 13C2 と 34S の行は出ず、窓に 18O の位置（重心から約 2.7 mDa）が
入れば 18O の行だけが立って、やはり `informational` になる。束を分離しない装置では `mz_tol` を
約 6 mDa / 電荷以上にする。窓は X の Δm/z で決まるので、同じ X と Y の M+2 の行は（どの位置の行でも）同じ期待比を持つ。

強い説明（`strong`）は、(a) 前駆体側の FoundInUpperMsMs リンク（裏付けつき）、(b) 同位体で強度比が
期待比に合う（裏付けつき・M+2 の取り違えの恐れが無い）、(c) アダクトの組・インソース断片で、Δm/z と
同時溶出に加えて試料間の相関 r ≥ min_r がそろう、のどれか。Δm/z が合うだけでは断定しない。

M+2 の位置は 13C2（+2.006710）・18O（+2.004245）・34S（+1.995796）・37Cl（+1.997050）・81Br（+1.997952）・
41K（+1.998119。`[M+K]+` など K を持つイオン）を見る。
18O は 13C2 から 2.46 mDa（を電荷で割った値）しか離れないので、窓（`mz_tol`）がその差以上なら 13C2 と分けられず、
X が 13C2 の窓に入るかぎり 18O の行（`isotope_M+2_18O`）は出さずに 13C2 の行（`isotope_M+2`）が 18O も含む。
窓がその差より狭い（分離された 18O を見る）ときは 18O を別の説明として見る。41K も同じく、13C2 から 8.59 mDa
（を電荷で割った値）以上の窓（既定の 10 mDa）では X が 13C2 の窓に入るかぎり 13C2 の行が 41K も含み、それより狭い
窓では 41K の行（`isotope_M+2_41K`）を出す（K を持つイオンの M+2 がどれも 13C2 と 41K の 2 行に当たって `ambiguous` に
なるのを避ける）。X が 13C2 の窓の外なら、18O・41K の行は X が当たった M+2 の位置のうちいちばん近いときだけ出す
（窓を広げて説明を失わない。`_SUBSUMED`）。
41K は 37Cl から 1.07 mDa、81Br から 0.17 mDa なので、K と Cl・Br を両方持つイオンではそれらと重なる。Y のイオンの組成が
分かっていてその元素（O・S・Cl・Br・K）を持たなければ、その説明は出さない（同位体の表に無い元素を含む組成でも、組成式が
読めればこの門には使い、期待比だけを求めない）。2 つ以上の M+2 が同じ Δm/z の窓に入るときは一方に決めず、どれにも
`ambiguous` を付けて並べ（強い説明にしない）、期待比（X の窓に入る M+2 の和。行によらず同じ）に合わないときだけ
情報（別物質）にする。

FoundInUpperMsMs のリンクは無向の対（実データで相互 100%、相手が高 m/z 側になるのは半々）。
断片は前駆体より軽いので、リンクだけで関係にするのも、リンクを強い理由にするのも、
X が Y より軽い（X の m/z < Y の m/z）ときに限る。upper（前駆体）側は m/z の大小で決める。
リンクは裏付け（supported）には向きによらず数える。

deps: msdial.adducts、msdial.peak_verification（組成式）、curation.isotope（同位体の組み合わせと質量差）。
"""
from __future__ import annotations

import math
import re
from bisect import bisect_left, bisect_right
from functools import lru_cache
from statistics import median

from metabolomix.curation.isotope import ISOTOPE_SHIFTS, isotopologues
from metabolomix.msdial.adducts import (SAME_MZ_TOL, ion_counts, mz_from_neutral, neutral_from_mz, parse_adduct,
                                        same_adduct)
from metabolomix.msdial.peak_verification import parse_formula

ISOTOPE_RATIO_FACTOR = 1.5
SUPPORT_LINKS = frozenset({"correl_similar", "chrom_similar", "found_in_upper_msms"})
#: 関係コード → (整数質量のずれ k（M+k）, 質量差 [Da], その位置の重い同位体の元素)。質量差は同位体の質量の差
#: （`curation.isotope.ISOTOPE_SHIFTS`）。期待比は X の窓に入る M+k の同位体の和（`expected_isotope_ratio`）。
ISOTOPE_DELTAS = {
    "isotope_M+1": (1, ISOTOPE_SHIFTS["C"][1], "C"),            # 13C（+1.003355）
    "isotope_M+2": (2, 2 * ISOTOPE_SHIFTS["C"][1], "C"),        # 13C2（+2.006710。18O と分けられない窓では 18O も含む）
    "isotope_M+2_18O": (2, ISOTOPE_SHIFTS["O"][2], "O"),        # 18O（+2.004245）。窓が 13C2 との差より狭いときだけ（`_SUBSUMED`）
    "isotope_M+2_34S": (2, ISOTOPE_SHIFTS["S"][2], "S"),        # 34S（+1.995796）
    "isotope_M+2_37Cl": (2, ISOTOPE_SHIFTS["Cl"][2], "Cl"),     # 37Cl（+1.997050）
    "isotope_M+2_81Br": (2, ISOTOPE_SHIFTS["Br"][2], "Br"),     # 81Br（+1.997952）
    "isotope_M+2_41K": (2, ISOTOPE_SHIFTS["K"][2], "K"),        # 41K（+1.998119）。窓が 13C2 との差より狭いときだけ（`_SUBSUMED`）
}
#: M+2 の同位体（Δm/z の窓が重なりうる。既定の 10 mDa では 13C2 と 37Cl の差 9.66 mDa、34S と 37Cl の差 1.25 mDa が
#: 窓の内側。13C2 と 34S の差 10.91 mDa は外側だが、分離できない M+2 の束の重心はその間に来る）。
M2_ISOTOPES = frozenset(code for code, (k, _, _) in ISOTOPE_DELTAS.items() if k == 2)
#: 窓（mz_tol）が 2 つの位置の差（を電荷で割った値）以上で、X が右の行にも当たるなら、左の行を出さずに右の行に含める
#: （窓では分けられない）。18O は 13C2 から 2.46 mDa で、ほぼすべての有機物のイオンが C と O を持つ——分けられない窓で
#: 両方を並べると、どの M+2 も `ambiguous`（強い説明にならない）になってしまう。41K（13C2 から 8.59 mDa）も同じで、
#: 既定の 10 mDa では K を持つイオン（`[M+K]+`）の M+2 がどれも 13C2 と 41K の両方に当たる。期待比は X の窓で決まり
#: 行によらないので、含めても比の判定は変わらない。X が右の行に当たらない（右の位置は X の窓の外。41K の位置から
#: 13C2 と逆の側へ「窓 − 8.59 mDa」より離れた X）なら、左の行は X が当たった M+2 の位置のうちいちばん近いときだけ
#: 右の行の代わりに出す——窓の幅だけで消すと、その X は狭い窓では説明され、広い窓では説明が付かない。いちばん近い
#: ときに限るのは、分離された 34S（18O から 8.45 mDa）のように、より近い位置で説明される X を遠い 18O と並べて
#: `ambiguous` にしないため。
_SUBSUMED = {"isotope_M+2_18O": "isotope_M+2", "isotope_M+2_41K": "isotope_M+2"}
#: 同位体の位置が X の窓に入るかの比較で、浮動小数の誤差だけを吸う許容 [Da]（行の窓と期待比の窓を一致させる）。
_WINDOW_EPS = 1e-9
#: 元素記号（と標識の D）。組成式の中にこれ以外の記号があれば読めない組成式とみなす（`C10Xq3` の `Xq`）。
_ELEMENTS = frozenset("""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr
Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt
Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc
Lv Ts Og D
""".split())
#: 組成式は元素記号と個数の並びだけ（`Unknown` や `None` のような語は組成式でない）。末尾の電荷（`C5H14NO+`）は
#: 読み飛ばす（`parse_formula` と同じ）。
_FORMULA_RE = re.compile(r"(?:[A-Z][a-z]?\d*)+(?:\d*[+-]|[+-]\d*)?")
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


@lru_cache(maxsize=4096)
def _ion_composition(formula: str | None, adduct: str | None) -> tuple | None:
    if not formula or not _FORMULA_RE.fullmatch(formula.strip()):
        return None
    counts = parse_formula(formula)
    if any(element not in _ELEMENTS for element in counts):
        return None
    ion = ion_counts(counts, adduct)
    counts = ion if ion is not None else counts
    return tuple(sorted(counts.items()))


def _text(value) -> str | None:
    return value if isinstance(value, str) else None


def _charge(adduct) -> int:
    parsed = parse_adduct(_text(adduct))
    return parsed.charge if parsed else 1


def ion_element_counts(formula, adduct) -> dict | None:
    """Y のイオンの元素の数（中性の組成式の n_mol 倍に付加体の原子を足したもの。付加体が読めない・差し引いて
    負になるなら中性の組成式のまま）。組成式が無い・読めない（元素記号でない記号を含む）なら None
    （「その元素を持たない」とは区別する）。同位体の表に無い元素（Se・Fe など）を含んでも数は返す。"""
    items = _ion_composition(_text(formula), _text(adduct))
    return None if items is None else dict(items)


@lru_cache(maxsize=4096)
def _isotopologues(items: tuple) -> tuple | None:
    terms = isotopologues(dict(items))
    return None if terms is None else tuple(terms)


def expected_isotope_ratio(formula, adduct, k: int, dmz: float, mz_tol: float) -> float | None:
    """X（Y からの Δm/z が `dmz`）が Y のイオンの M+k の同位体ピークだとしたときの、X の強度の Y に対する期待比。

    Y のイオンの M+k の同位体の組み合わせ（`curation.isotope.isotopologues`）のうち、m/z（Y + 質量差 / 価数）が
    X から `mz_tol` 以内に来るものの比の和。組成式が無い・読めない・同位体の表に無い元素を含むなら None。
    """
    items = _ion_composition(_text(formula), _text(adduct))
    terms = None if items is None else _isotopologues(items)
    if terms is None:
        return None
    charge = _charge(adduct)
    return sum(ratio for term_k, shift, ratio in terms
               if term_k == k and abs(dmz - shift / charge) <= mz_tol + _WINDOW_EPS)


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
            # Y 自身の付加体（綴りでなく組成で比べる）は Y そのもの
            if other is None or same_adduct(other.name, y_adduct.name):
                continue
            theoretical = mz_from_neutral(neutral, other)
            if abs(theoretical - y["mz"]) <= SAME_MZ_TOL:
                continue                # Y と同じ m/z（[2M+2H]2+ と [M+H]+ など）。Δm/z が何も区別しない
            delta = x["mz"] - theoretical
            if abs(delta) <= mz_tol:
                found.append((f"adduct:{other.name}/{y_adduct.name}", delta))
    charge = y_adduct.charge if y_adduct else 1
    counts = ion_element_counts(y.get("formula"), y.get("adduct"))
    isotopes = {}
    for relation, (_, shift, element) in ISOTOPE_DELTAS.items():
        if element != "C" and counts is not None and not counts.get(element):
            continue                    # Y のイオンにその元素が無い（その同位体では説明できない）
        delta = x["mz"] - (y["mz"] + shift / charge)
        if abs(delta) <= mz_tol:
            isotopes[relation] = delta
    for relation, delta in isotopes.items():
        into = _SUBSUMED.get(relation)
        if into is not None and mz_tol >= abs(ISOTOPE_DELTAS[into][1] - ISOTOPE_DELTAS[relation][1]) / charge:
            # 窓では 13C2 と分けられない 18O・41K: X が 13C2 の行にも当たれば 13C2 の行に含める。当たらなければ
            # （X は 13C2 の窓の外）、当たった M+2 の位置のうち X にいちばん近いときだけ 13C2 の代わりに立てる
            # （窓の幅だけで消すと、狭い窓では説明される X が広い窓で説明を失う）
            if into in isotopes or any(abs(other) < abs(delta) for code, other in isotopes.items()
                                       if code != relation and code in M2_ISOTOPES):
                continue
        found.append((relation, delta))
    for relation, loss in INSOURCE_LOSSES:
        delta = x["mz"] - (y["mz"] - loss / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    return found


def _isotope_check(entry: dict, relation: str, x: dict, y: dict, mz_tol: float) -> None:
    """同位体の行に強度比・期待比・`isotope_consistent` を書く（期待比が求まらなければ consistent は None）。

    期待比は Y のイオンの M+k の同位体のうち X の窓（`mz_tol`）に入るものの和（`expected_isotope_ratio`）。"""
    k = ISOTOPE_DELTAS[relation][0]
    expected = expected_isotope_ratio(y.get("formula"), y.get("adduct"), k, x["mz"] - y["mz"], mz_tol)
    ratio = median_ratio(x.get("heights") or [], y.get("heights") or [])
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
                _isotope_check(entry, relation, target, y, mz_tol)
                entry["informational"] = entry["isotope_consistent"] is False
            entries.append(entry)
        m2 = [e for e in entries if e["relation"] in M2_ISOTOPES]
        if len(m2) > 1:
            # 窓が重なる M+2 は一方に決めない（強い説明にしない）。期待比は X の窓に入る M+2 の和で行によらず
            # 同じなので、期待比に合わないとき（別物質の情報）もそろう。
            for e in m2:
                e["ambiguous"] = True
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
