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
2 倍、`[M+Cl]-` は Cl を 1 個持つ）の整数質量の同位体分布（`curation.isotope.envelope_ratios`）から求める。
M+1 はその M+1 の和（13C・15N・2H・17O・33S）、M+2 はその M+2 の和（13C2・18O・34S・37Cl・81Br・41K と
交差項）。M+2 の行はどの同位体の位置で見つけても同じ和を期待比にする——分離できない M+2 の束の強度は和で、
分離された 1 つの同位体の強度もその和を超えないので、同位体を別物質と言わない側の上限になる。

強い説明（`strong`）は、(a) 前駆体側の FoundInUpperMsMs リンク（裏付けつき）、(b) 同位体で強度比が
期待比に合う（裏付けつき・M+2 の取り違えの恐れが無い）、(c) アダクトの組・インソース断片で、Δm/z と
同時溶出に加えて試料間の相関 r ≥ min_r がそろう、のどれか。Δm/z が合うだけでは断定しない。

M+2 の位置は 13C2（+2.006710）・18O（+2.004246）・34S（+1.995796）・37Cl（+1.997050）・81Br（+1.997952）を見る。
18O は 13C2 から 2.46 mDa（を電荷で割った値）しか離れないので、窓（`mz_tol`）がその差以上なら 13C2 と分けられず、
18O の行（`isotope_M+2_18O`）は出さずに 13C2 の行（`isotope_M+2`）が 18O も含む——既定の 10 mDa ではいつもこちら。
窓がその差より狭い（分離された 18O を見る）ときだけ 18O を別の説明として見る。Y のイオンの組成が分かっていて
その元素（O・S・Cl・Br）を持たなければ、その説明は出さない。2 つ以上の M+2 が同じ Δm/z の窓に入るときは
一方に決めず、どれにも `ambiguous` を付けて並べ（強い説明にしない）、期待比（M+2 の和）に合わないときだけ
情報（別物質）にする。

FoundInUpperMsMs のリンクは無向の対（実データで相互 100%、相手が高 m/z 側になるのは半々）。
断片は前駆体より軽いので、リンクだけで関係にするのも、リンクを強い理由にするのも、
X が Y より軽い（X の m/z < Y の m/z）ときに限る。upper（前駆体）側は m/z の大小で決める。
リンクは裏付け（supported）には向きによらず数える。

deps: msdial.adducts、msdial.peak_verification（組成式）、curation.isotope（同位体分布）。
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from functools import lru_cache
from statistics import median

from metabolomix.curation.isotope import ABUNDANCES, envelope_ratios
from metabolomix.msdial.adducts import ion_counts, mz_from_neutral, neutral_from_mz, parse_adduct, same_adduct
from metabolomix.msdial.peak_verification import parse_formula

ISOTOPE_RATIO_FACTOR = 1.5
SUPPORT_LINKS = frozenset({"correl_similar", "chrom_similar", "found_in_upper_msms"})
#: 関係コード → (整数質量のずれ k（M+k）, 質量差 [Da], その位置の重い同位体の元素)。質量差は同位体の質量の差から
#: 計算した値。期待比は位置によらず M+k の和（`expected_isotope_ratios`）。
ISOTOPE_DELTAS = {
    "isotope_M+1": (1, 1.003355, "C"),
    "isotope_M+2": (2, 2.006710, "C"),           # 13C2（18O と分けられない窓では 18O も含む）
    "isotope_M+2_18O": (2, 2.004246, "O"),       # 18O − 16O。窓が 13C2 との差より狭いときだけ（`_SUBSUMED`）
    "isotope_M+2_34S": (2, 1.995796, "S"),       # 34S − 32S
    "isotope_M+2_37Cl": (2, 1.997050, "Cl"),     # 37Cl − 35Cl
    "isotope_M+2_81Br": (2, 1.997952, "Br"),     # 81Br − 79Br
}
#: M+2 の同位体（Δm/z の窓が重なりうる。既定の 10 mDa では 13C2 と 37Cl の差 9.66 mDa、34S と 37Cl の差 1.25 mDa が
#: 窓の内側。13C2 と 34S の差 10.91 mDa は外側だが、分離できない M+2 の束の重心はその間に来る）。
M2_ISOTOPES = frozenset(code for code, (k, _, _) in ISOTOPE_DELTAS.items() if k == 2)
#: 窓（mz_tol）が 2 つの位置の差（を電荷で割った値）以上なら、左の位置を見ずに右の行に含める（窓では分けられない）。
#: 18O は 13C2 から 2.46 mDa で、ほぼすべての有機物のイオンが C と O を持つ——分けられない窓で両方を並べると、
#: どの M+2 も `ambiguous`（強い説明にならない）になってしまう。期待比はどちらの行も M+2 の和なので、含めても
#: 比の判定は変わらない。
_SUBSUMED = {"isotope_M+2_18O": "isotope_M+2"}
#: Y と同じ m/z になる付加体（Δm/z が何も区別しない）を除くための許容 [Da]。浮動小数の誤差だけを吸う。
SAME_MZ_TOL = 1e-6
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
    if not formula:
        return None
    try:
        counts = parse_formula(formula)
    except ValueError:
        return None
    ion = ion_counts(counts, adduct)
    counts = ion if ion is not None else counts
    if any(element not in ABUNDANCES for element in counts):
        return None
    return tuple(sorted(counts.items()))


def _text(value) -> str | None:
    return value if isinstance(value, str) else None


def ion_element_counts(formula, adduct) -> dict | None:
    """Y のイオンの元素の数（中性の組成式の n_mol 倍に付加体の原子を足したもの。付加体が読めない・差し引いて
    負になるなら中性の組成式のまま）。組成式が無い・読めない・同位体の表に無い元素を含むなら None
    （「その元素を持たない」とは区別する）。"""
    items = _ion_composition(_text(formula), _text(adduct))
    return None if items is None else dict(items)


@lru_cache(maxsize=4096)
def _ratios(items: tuple) -> tuple | None:
    ratios = envelope_ratios(dict(items))
    return None if ratios is None else tuple(ratios)


def expected_isotope_ratios(formula, adduct) -> tuple[float, float] | None:
    """Y のイオンの (M+1 の和 / M, M+2 の和 / M)。組成が分からなければ None。"""
    items = _ion_composition(_text(formula), _text(adduct))
    ratios = None if items is None else _ratios(items)
    return None if ratios is None else (ratios[1], ratios[2])


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
    for relation, (_, shift, element) in ISOTOPE_DELTAS.items():
        if element != "C" and counts is not None and not counts.get(element):
            continue                    # Y のイオンにその元素が無い（その同位体では説明できない）
        into = _SUBSUMED.get(relation)
        if into is not None and mz_tol >= abs(ISOTOPE_DELTAS[into][1] - shift) / charge:
            continue                    # 窓では 13C2 と分けられない 18O は 13C2 の行に含める
        delta = x["mz"] - (y["mz"] + shift / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    for relation, loss in INSOURCE_LOSSES:
        delta = x["mz"] - (y["mz"] - loss / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    return found


def _isotope_check(entry: dict, relation: str, x: dict, y: dict) -> None:
    """同位体の行に強度比・期待比・`isotope_consistent` を書く（期待比が求まらなければ consistent は None）。

    期待比は Y のイオンの M+k の和（M+2 の行はどの位置でも同じ M+2 の和）。"""
    k = ISOTOPE_DELTAS[relation][0]
    ratios = expected_isotope_ratios(y.get("formula"), y.get("adduct"))
    ratio = median_ratio(x.get("heights") or [], y.get("heights") or [])
    expected = None if ratios is None else ratios[k - 1]
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
            # 窓が重なる M+2 は一方に決めない。期待比（M+2 の和）に合わないときだけ別物質（情報）とする。
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
