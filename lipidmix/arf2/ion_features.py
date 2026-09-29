"""`.arf2` Key 10（IonFeatureCharacter）の復号。Key 番号の正準は docs/schema/IonFeatureCharacter.md。

MS-DIAL がアラインメント時に計算したスポット間の関係（同位体・アダクト・クロマトグラムの類似・
上位イオンの MS/MS に断片として出ている・試料間の相関）を、候補付けの ④ の裏付けに使う。
deps: arf2.reader だけ。
"""
from __future__ import annotations

from lipidmix.arf2.reader import load_raw_spots

LINK_KINDS = {0: "same_feature", 1: "isotope", 2: "adduct", 3: "chrom_similar",
              4: "found_in_upper_msms", 5: "correl_similar"}
_EMPTY = {"links": [], "isotope_weight": None, "isotope_parent": None, "peak_group": None}


def _int_at(values: list, key: int):
    value = values[key] if len(values) > key else None
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def decode_ion_features(raw_spot: list) -> dict:
    character = raw_spot[10] if len(raw_spot) > 10 else None
    if not isinstance(character, list):
        return dict(_EMPTY, links=[])
    links = []
    raw_links = character[3] if len(character) > 3 and isinstance(character[3], list) else []
    for link in raw_links:
        if isinstance(link, list) and len(link) >= 2 and isinstance(link[0], int) \
                and isinstance(link[1], int):
            links.append({"spot_id": int(link[0]),
                          "kind": LINK_KINDS.get(int(link[1]), f"unknown_{int(link[1])}")})
    return {"links": links, "isotope_weight": _int_at(character, 4),
            "isotope_parent": _int_at(character, 5), "peak_group": _int_at(character, 6)}


def load_ion_features(arf2_path) -> dict[int, dict]:
    return {int(raw[0]): decode_ion_features(raw) for raw in load_raw_spots(arf2_path)}
