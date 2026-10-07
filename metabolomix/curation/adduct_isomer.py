"""別アダクトの取り違え（純関数）。理由コード `adduct_isomer_of` / `adduct_isomer_minor_of`。

スポット X の m/z が、同時に溶出する（|ΔRT| ≤ `adduct_isomer_drt`）別物質のスポット Y の中性質量の
別アダクトで説明でき（`adduct_isomer_ppm` 以内）、Y の方が証拠が強いとき、X は Y の別アダクトを
別物質として注釈した疑いがある（実例: 20260930_EV neg #173「PI 41:2」[M-H]- は #189 DGDG 35:1
[M+CH3COO]- の [M-H]-）。個々のスポットの証拠（Δppm・MS2）だけでは見えない——#173 は参照一致で
`ok` だった——ので、アラインメント全体の注釈付きスポットを相手にする。

- 仮定するアダクトは X の極性の `DEFAULT_ADDUCTS` と X 自身のアダクト。Y 自身のアダクトは除く
  （同じアダクトで同じ m/z なら重複スポットで、別の話）。
- 組成式もクラスも同じ組は同じ分子種の別アダクト（正しい注釈）なので対象外。
- 強さは段階（`evidence_tier`）で比べる。Y の段階が上なら `strong`、同じ段階なら Y の
  HeightAverage が X の `adduct_isomer_height_ratio` 倍以上のときだけ `minor`。
- 複数の Y が説明するときは、strong → 段階 → 強度の順に強い 1 件を返す。

deps: msdial.adducts、msdial.analysis_params（既定アダクト）、arf2.match_results（名前接頭辞）。
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from functools import lru_cache

from metabolomix.arf2.match_results import name_prefix
from metabolomix.curation.judge import is_chain_level_name
from metabolomix.msdial.adducts import mz_from_neutral, neutral_from_mz, parse_adduct
from metabolomix.msdial.analysis_params import DEFAULT_ADDUCTS

_UNANNOTATED = {"", "unknown"}         # evidence.select_spots と同じ「未注釈」
_POLARITY_KEY = {"+": "Positive", "-": "Negative"}


@lru_cache(maxsize=None)
def _adduct(name):
    return parse_adduct(name)


def evidence_tier(name, match) -> int:
    """注釈の証拠の段階。3 = MS2 で鎖組成まで裏付け（鎖レベルの名前 ∧ is_lipid_chains_match）、
    2 = MS2 が参照と一致、1 = MS2 はある（low score）、0 = MS2 無し・照合結果無し。"""
    if not match or not match.get("has_msms") or name_prefix(name) in ("no MS2", "w/o MS2"):
        return 0
    if match.get("is_lipid_chains_match") and is_chain_level_name(name):
        return 3
    return 2 if match.get("is_reference_matched") else 1


def pool_entry(spot: dict, annotation: dict | None) -> dict | None:
    """`.arf2` のカタログ行（`load_catalog`）と注釈（`load_spot_annotations` の 1 件）を、
    比べる相手の 1 件にする。未注釈・アダクトを解析できないスポットは None。"""
    name = (spot.get("Name") or "").strip()
    if name.lower() in _UNANNOTATED or _adduct(spot.get("AdductType")) is None:
        return None
    if spot.get("MassCenter") is None or spot.get("RT") is None:
        return None
    return {"spot_id": spot["MasterAlignmentID"], "name": name, "adduct": spot["AdductType"],
            "mz": float(spot["MassCenter"]), "rt": float(spot["RT"]),
            "formula": spot.get("Formula") or "", "ontology": spot.get("Ontology") or "",
            "height": float(spot.get("HeightAverage") or 0.0),
            "tier": evidence_tier(name, (annotation or {}).get("representative"))}


def _same_species(x: dict, y: dict) -> bool:
    return bool(x["formula"]) and x["formula"] == y["formula"] and \
        x["ontology"].lower() == y["ontology"].lower()


def _severity(x: dict, y: dict, th: dict) -> str | None:
    if y["tier"] > x["tier"]:
        return "strong"
    if y["tier"] == x["tier"] and x["height"] > 0 and y["height"] >= th["adduct_isomer_height_ratio"] * x["height"]:
        return "minor"
    return None


def _candidate_adducts(x_adduct) -> list:
    names = DEFAULT_ADDUCTS[_POLARITY_KEY[x_adduct.polarity]]
    adducts = [_adduct(n) for n in names]
    return [a for a in adducts if a is not None] + ([x_adduct] if x_adduct.name not in names else [])


def find_adduct_isomer(target: dict, pool: list[dict], th: dict) -> dict | None:
    """`target`（`pool_entry` の形）を別アダクトとして説明する、証拠の強い相手。無ければ None。"""
    x_adduct = _adduct(target["adduct"])
    ordered = sorted(pool, key=lambda p: p["rt"])
    rts = [p["rt"] for p in ordered]
    window = th["adduct_isomer_drt"]
    lo, hi = bisect_left(rts, target["rt"] - window), bisect_right(rts, target["rt"] + window)
    best, best_key = None, None
    for y in ordered[lo:hi]:
        y_adduct = _adduct(y["adduct"])
        if y["spot_id"] == target["spot_id"] or y_adduct.polarity != x_adduct.polarity \
                or _same_species(target, y):
            continue
        severity = _severity(target, y, th)
        if severity is None:
            continue
        neutral = neutral_from_mz(y["mz"], y_adduct)
        for adduct in _candidate_adducts(x_adduct):
            if adduct.name == y_adduct.name:
                continue
            expected = mz_from_neutral(neutral, adduct)
            ppm = (target["mz"] - expected) / expected * 1e6
            if abs(ppm) > th["adduct_isomer_ppm"]:
                continue
            key = (severity == "strong", y["tier"], y["height"], -abs(ppm))
            if best_key is None or key > best_key:
                best_key = key
                best = {"of": y["spot_id"], "of_name": y["name"], "of_adduct": y["adduct"],
                        "as_adduct": adduct.name, "ppm": round(ppm, 2),
                        "drt": round(target["rt"] - y["rt"], 4), "severity": severity,
                        "tier": target["tier"], "of_tier": y["tier"],
                        "height_ratio": (round(y["height"] / target["height"], 2)
                                         if target["height"] > 0 else None)}
    return best
