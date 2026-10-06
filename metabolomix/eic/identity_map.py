"""ARF2 の同定情報から EIC 描画対象スポットを選び、rt/mz で対応付けを検証する。

deps: arf2_reader のみ。EIC バイナリ読み出し・matplotlib・MCP には依存しない。
`spot_id = AlignmentID` の対応は MS-DIAL のバッチによって崩れうるため、
`verify_spot_match()` による rt/mz 照合を必須の安全弁として使うこと。
"""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict

RT_TOLERANCE = 0.02   # min
MZ_TOLERANCE = 0.01   # Da
MAX_CANDIDATES = 300


class IdentityCandidate(TypedDict):
    spot_id: int
    name: str
    ontology: str
    adduct: str
    rt: float
    mz: float
    height_average: float


def load_arf2_records(arf2_path: str | Path) -> list[dict]:
    """.arf2 を読み、`extract_arf2_data` 形式の辞書リストを返す。

    実体は `arf2.reader.load_catalog`（同一ファイルなら再パースしない共有キャッシュ）。
    返るリストは共有されるため変更しないこと。
    """
    from metabolomix.arf2.reader import load_catalog

    return load_catalog(Path(arf2_path))


def _matches(record: dict, name_queries: list[str], ontology_set: set[str]) -> bool:
    ontology = str(record.get("Ontology") or "")
    if ontology_set and ontology.casefold() in ontology_set:
        return True
    lowered = str(record.get("Name") or "").casefold()
    return any(query in lowered for query in name_queries)


def select_identity_candidates(
    records: list[dict],
    *,
    names: list[str] | None = None,
    ontologies: list[str] | None = None,
    max_candidates: int = MAX_CANDIDATES,
) -> tuple[list[IdentityCandidate], list[str], int]:
    """クエリに一致する ARF2 スポット候補（spot_id 昇順）と caveat 文字列を返す。

    `names` は Name への大小無視の部分一致、`ontologies` は Ontology への大小無視の
    完全一致。両者の結果は AlignmentID で重複排除した和集合になる。

    3 番目の戻り値 `total_matched` は `max_candidates` による予備選抜より前の
    クエリ一致件数（予備選抜が発生しなければ `len(candidates)` と同じ）。
    """
    name_queries = [str(query).casefold() for query in (names or []) if str(query).strip()]
    ontology_set = {
        str(item).casefold() for item in (ontologies or []) if str(item).strip()
    }
    if not name_queries and not ontology_set:
        raise ValueError("names または ontologies のいずれかを指定してください。")

    candidates: list[IdentityCandidate] = []
    seen: set[int] = set()
    for record in records:
        spot_id = record.get("AlignmentID")
        if isinstance(spot_id, bool) or not isinstance(spot_id, int):
            continue
        if spot_id in seen or not _matches(record, name_queries, ontology_set):
            continue
        seen.add(spot_id)
        candidates.append({
            "spot_id": spot_id,
            "name": str(record.get("Name") or ""),
            "ontology": str(record.get("Ontology") or ""),
            "adduct": str(record.get("AdductType") or ""),
            "rt": float(record.get("RT") or 0.0),
            "mz": float(record.get("MassCenter") or 0.0),
            "height_average": float(record.get("HeightAverage") or 0.0),
        })

    total_matched = len(candidates)
    caveats: list[str] = []
    if len(candidates) > max_candidates:
        candidates.sort(key=lambda item: item["height_average"], reverse=True)
        candidates = candidates[:max_candidates]
        caveats.append(
            f"クエリ一致 {total_matched} 件のうち、ARF2 HeightAverage 上位 {max_candidates} 件だけを "
            "EIC 読み出し対象にしました（クエリを絞ると全件を評価できます）。"
        )

    candidates.sort(key=lambda item: item["spot_id"])
    return candidates, caveats, total_matched


def verify_spot_match(
    candidate: IdentityCandidate,
    spot: dict,
    *,
    rt_tolerance: float = RT_TOLERANCE,
    mz_tolerance: float = MZ_TOLERANCE,
) -> str | None:
    """対応付けが妥当なら None、外れていれば除外理由の語を返す。"""
    if abs(float(spot["rt"]) - float(candidate["rt"])) > rt_tolerance:
        return "rt_mismatch"
    if abs(float(spot["mz"]) - float(candidate["mz"])) > mz_tolerance:
        return "mz_mismatch"
    return None
