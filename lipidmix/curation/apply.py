"""差次的エクスポートへのフラグ反映。spec §6。

契約 15 列は変えない。wrong のスポットは同定なしとして扱う＝InChIKey 付きの行だけを出す
エクスポートからは落ちる。suspect は行も値も変えない。適用状況はメタ行 1 本で宣言する
（massbank-context の `_parse_meta_line` は未知のキーを保持するだけで落ちない）。
フラグ 0 件なら何も足さず、出力は現行と完全に同じ。
"""
from __future__ import annotations

from pathlib import Path

from lipidmix.curation.flags import FlagStore, alignment_key, curation_dir


def flags_for_arf2(arf2_path) -> dict:
    key = alignment_key(arf2_path)
    store = FlagStore(curation_dir(arf2_path))
    effective = store.effective(key["alignment_sha256"])
    return {"wrong": {s for s, r in effective.items() if r["flag"] == "wrong"},
            "suspect": {s for s, r in effective.items() if r["flag"] == "suspect"},
            "digest": store.digest(key["alignment_sha256"]), "n": len(effective)}


def arf2_for_mztab(mztab_path) -> Path | None:
    mztab_path = Path(mztab_path)
    matches = [p for p in mztab_path.parent.glob("*.arf2")
               if p.name[: -len(".arf2")] in mztab_path.name]
    return matches[0] if len(matches) == 1 else None


def filter_rows(rows, flag_set: dict, *, key):
    kept, excluded, suspect = [], 0, 0
    for row in rows:
        spot = key(row)
        if spot in flag_set["wrong"]:
            excluded += 1
            continue
        if spot in flag_set["suspect"]:
            suspect += 1
        kept.append(row)
    return kept, {"wrong_excluded": excluded, "suspect": suspect}


def meta_line(state: str, flag_set: dict | None, stats: dict | None) -> str | None:
    if not flag_set or flag_set["n"] == 0:
        return None
    parts = [f"# curation = {state}", f"curation_flags = {flag_set['n']}"]
    if stats is not None:
        parts += [f"curation_wrong_excluded = {stats['wrong_excluded']}",
                  f"curation_suspect = {stats['suspect']}"]
    parts.append(f"curation_flags_sha256 = {flag_set['digest']}")
    return "\t".join(parts)
