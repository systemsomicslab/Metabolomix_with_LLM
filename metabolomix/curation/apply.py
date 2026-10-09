"""差次的エクスポートへのフラグ反映。spec §6。

契約 15 列は変えない。wrong のスポットは同定なしとして扱う＝InChIKey 付きの行だけを出す
エクスポートからは落ちる。suspect は行も値も変えない。redundant は除外、assign は同定を
置き換える（未注釈だったスポットは InChIKey を得てエクスポートに現れる）。
適用状況はメタ行 1 本で宣言する
（massbank-context の `_parse_meta_line` は未知のキーを保持するだけで落ちない）。
フラグ 0 件なら何も足さず、出力は現行と完全に同じ。
"""
from __future__ import annotations

import re
from pathlib import Path

from metabolomix.msdial import tags as msdial_tags
from metabolomix.curation.flags import (
    FlagStore, alignment_key, curation_dir, effective_flags, flags_digest, orphaned_count,
    split_decisions,
)


def flags_for_arf2(arf2_path) -> dict:
    """その `.arf2` に当たる有効フラグ。`orphaned` は同じファイル名で sha256 が違う版に
    付いたまま当たらないフラグの件数（メタ行には出さず、payload の warnings で知らせる）。

    `flags.jsonl` が無ければ `.arf2` を hash しない（数百 MB の全読みを省く）。
    読めない行があれば `flags.FlagFileError` をそのまま投げる（呼び出し側がエラーにする）。
    `assign` は `{spot: 判断の行}`（同定の置き換え）、`redundant` は除外するスポットの集合。
    """
    store = FlagStore(curation_dir(arf2_path))
    if not store.exists():
        return {"wrong": set(), "suspect": set(), "assign": {}, "redundant": set(),
                "digest": flags_digest({}), "n": 0, "orphaned": 0}
    rows = store.rows()
    key = alignment_key(arf2_path)
    effective = effective_flags(rows, key["alignment_sha256"])
    decisions = split_decisions(effective)
    return {"wrong": decisions["wrong"], "suspect": decisions["suspect"],
            "assign": decisions["assign"], "redundant": set(decisions["redundant"]),
            "digest": flags_digest(effective), "n": len(effective),
            "orphaned": orphaned_count(rows, key)}


def confirmed_spots(arf2_path, effective: dict) -> set[int]:
    """注釈が確かめられたスポット: 記録の confirmed と、`_tags.xml` の Confirmed タグ（MS-DIAL の GUI で
    付けたものも含む）の和から、最新の記録が wrong / suspect のスポットを除く。自動判定の likely_wrong を
    当てない（初期選択・exclude_auto_likely_wrong・候補付けの対象）ために使う。タグファイルが無い・
    読めなければ記録だけで決める（判断を止めない）。"""
    tagged: set[int] = set()
    tag_path = msdial_tags.alignment_tag_path(arf2_path)
    if tag_path.is_file():
        try:
            tagged = {int(s) for s, ids in msdial_tags.parse_tag_file(tag_path)["peaks"].items()
                      if msdial_tags.CONFIRMED_TAG_ID in ids}
        except Exception:  # noqa: BLE001 - 壊れたタグファイルでも記録の判断は使う
            tagged = set()
    recorded = {int(s) for s, row in effective.items() if row.get("flag") == "confirmed"}
    overridden = {int(s) for s, row in effective.items() if row.get("flag") in ("wrong", "suspect")}
    return (tagged | recorded) - overridden


def identity_for(spot_id, flag_set) -> dict | None:
    """assign の判断が置き換える同定（`name` / `ontology` / `inchikey`）。無ければ None。"""
    row = (flag_set or {}).get("assign", {}).get(spot_id)
    if row is None:
        return None
    return {"name": row.get("name") or "", "ontology": row.get("ontology") or "",
            "inchikey": row.get("inchikey") or ""}


def override_identity(catalog: dict, flag_set) -> dict:
    """assign のスポットの `Name` / `Ontology` / `InChIKey` を置き換えた新しい catalog。
    置き換えた行には `_curation = "assign"` を付ける。元の dict は変えない。"""
    out = dict(catalog)
    for spot_id in (flag_set or {}).get("assign", {}):
        identity = identity_for(spot_id, flag_set)
        if spot_id in out and identity is not None:
            out[spot_id] = {**out[spot_id], "Name": identity["name"], "Ontology": identity["ontology"],
                            "InChIKey": identity["inchikey"], "_curation": "assign"}
    return out


def arf2_for_mztab(mztab_path) -> Path | None:
    """mzTab-M のファイル名にバッチ語幹（`.arf2` の拡張子を除いた名前）を含む `.arf2`。
    語幹の直後に数字が続くものは別バッチ（`_2_3` と `_2_30`）なので数えない。"""
    mztab_path = Path(mztab_path)
    matches = [p for p in mztab_path.parent.glob("*.arf2")
               if re.search(re.escape(p.name[: -len(".arf2")]) + r"(?!\d)", mztab_path.name)]
    return matches[0] if len(matches) == 1 else None


def payload_summary(state: str | None, flag_set: dict | None, stats: dict | None) -> dict:
    """エクスポートの成功 payload に載せる `curation` の要約（両経路で同じ形）。
    件数は実際に除外を行った `applied` のときだけ数字で、それ以外は None。
    assign / redundant が無ければ wrong / suspect だけの現行と同じ 4 キー。"""
    summary = {"state": state,
               "wrong_excluded": stats["wrong_excluded"] if stats else None,
               "suspect": stats["suspect"] if stats else None,
               "orphaned": (flag_set or {}).get("orphaned", 0)}
    if (flag_set or {}).get("assign") or (flag_set or {}).get("redundant"):
        summary["assigned"] = stats.get("assigned") if stats else None
        summary["redundant_excluded"] = stats.get("redundant_excluded") if stats else None
    return summary


def filter_rows(rows, flag_set: dict, *, key):
    """wrong と redundant の行を落とす。stats は内部の値（メタ行と payload の材料）。"""
    kept, stats = [], {"wrong_excluded": 0, "redundant_excluded": 0, "suspect": 0, "assigned": 0}
    for row in rows:
        spot = key(row)
        if spot in flag_set["wrong"]:
            stats["wrong_excluded"] += 1
            continue
        if spot in flag_set.get("redundant", ()):
            stats["redundant_excluded"] += 1
            continue
        if spot in flag_set["suspect"]:
            stats["suspect"] += 1
        if spot in flag_set.get("assign", {}):
            stats["assigned"] += 1
        kept.append(row)
    return kept, stats


def meta_line(state: str, flag_set: dict | None, stats: dict | None) -> str | None:
    if not flag_set or flag_set["n"] == 0:
        return None
    parts = [f"# curation = {state}", f"curation_flags = {flag_set['n']}"]
    if stats is not None:
        parts += [f"curation_wrong_excluded = {stats['wrong_excluded']}",
                  f"curation_suspect = {stats['suspect']}"]
        if flag_set.get("assign") or flag_set.get("redundant"):
            parts += [f"curation_assigned = {stats.get('assigned', 0)}",
                      f"curation_redundant_excluded = {stats.get('redundant_excluded', 0)}"]
    parts.append(f"curation_flags_sha256 = {flag_set['digest']}")
    return "\t".join(parts)
