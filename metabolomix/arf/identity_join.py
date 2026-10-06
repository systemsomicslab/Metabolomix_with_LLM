"""差次的結果の spot id に .arf2 の同定情報を結合する(純関数)。

ファイルの読み込みは行わない。呼び出し側が `.arf2` を読んで
`MasterAlignmentID -> スポット辞書` を渡す。純関数にしておくのは、
結合の規則そのものを単独でテストしたいためである。

InChIKey が無い行は落とす。パスウェイ照会の鍵が無い行を残すと、下流で
「照会したがパスウェイが無かった化合物」と区別が付かなくなる。落とした件数は
report で必ず返す(黙って減らさない)。
"""
from __future__ import annotations

import re

_SPOT_ID_RE = re.compile(r"^Spot_(\d+)_")


def spot_id_of(feature_name) -> int | None:
    """feature 名 `Spot_<id>_<prop>` から MasterAlignmentID を取る。"""
    match = _SPOT_ID_RE.match(str(feature_name))
    return int(match.group(1)) if match else None


def join_identity(results: list[dict],
                  catalog: dict[int, dict]) -> tuple[list[dict], dict]:
    """差次的結果に .arf2 の同定情報を結合する。

    引数:
        results: two_group_test + add_fdr の出力。
        catalog: MasterAlignmentID -> .arf2 スポット辞書。

    戻り値: (rows, report)。rows は InChIKey を持つ行だけ。
    """
    rows: list[dict] = []
    n_unannotated = 0
    for result in results:
        spot_id = spot_id_of(result.get("feature"))
        spot = catalog.get(spot_id) if spot_id is not None else None
        inchikey = str((spot or {}).get("InChIKey") or "").strip()
        if not inchikey:
            n_unannotated += 1
            continue
        rows.append({
            "spot_id": spot_id,
            "name": str(spot.get("Name") or "").strip(),
            "ontology": str(spot.get("Ontology") or "").strip(),
            "inchikey": inchikey,
            "smiles": str(spot.get("SMILES") or "").strip(),
            "mz": spot.get("MassCenter"),
            "rt": spot.get("RT"),
            "log2fc": result.get("log2fc"),
            "p_value": result.get("p"),
            "q_value": result.get("q"),
            "mean_a": result.get("mean_a"),
            "mean_b": result.get("mean_b"),
        })
    report = {
        "n_features_total": len(results),
        "n_with_inchikey": len(rows),
        "n_unannotated": n_unannotated,
    }
    return rows, report
