"""差次的エクスポートへのフラグ反映。spec §6。

契約 15 列は変えない。wrong のスポットは同定なしとして扱う＝InChIKey 付きの行だけを出す
エクスポートからは落ちる。suspect は行も値も変えない。適用状況はメタ行 1 本で宣言する
（massbank-context の `_parse_meta_line` は未知のキーを保持するだけで落ちない）。
フラグ 0 件なら何も足さず、出力は現行と完全に同じ。
"""
from __future__ import annotations

import re
from pathlib import Path

from lipidmix.curation.flags import (
    FlagStore, alignment_key, curation_dir, effective_flags, flags_digest, orphaned_count,
)


def flags_for_arf2(arf2_path) -> dict:
    """その `.arf2` に当たる有効フラグ。`orphaned` は同じファイル名で sha256 が違う版に
    付いたまま当たらないフラグの件数（メタ行には出さず、payload の warnings で知らせる）。

    `flags.jsonl` が無ければ `.arf2` を hash しない（数百 MB の全読みを省く）。
    読めない行があれば `flags.FlagFileError` をそのまま投げる（呼び出し側がエラーにする）。
    """
    store = FlagStore(curation_dir(arf2_path))
    if not store.exists():
        return {"wrong": set(), "suspect": set(), "digest": flags_digest({}), "n": 0,
                "orphaned": 0}
    rows = store.rows()
    key = alignment_key(arf2_path)
    effective = effective_flags(rows, key["alignment_sha256"])
    return {"wrong": {s for s, r in effective.items() if r["flag"] == "wrong"},
            "suspect": {s for s, r in effective.items() if r["flag"] == "suspect"},
            "digest": flags_digest(effective), "n": len(effective),
            "orphaned": orphaned_count(rows, key)}


def arf2_for_mztab(mztab_path) -> Path | None:
    """mzTab-M のファイル名にバッチ語幹（`.arf2` の拡張子を除いた名前）を含む `.arf2`。
    語幹の直後に数字が続くものは別バッチ（`_2_3` と `_2_30`）なので数えない。"""
    mztab_path = Path(mztab_path)
    matches = [p for p in mztab_path.parent.glob("*.arf2")
               if re.search(re.escape(p.name[: -len(".arf2")]) + r"(?!\d)", mztab_path.name)]
    return matches[0] if len(matches) == 1 else None


def payload_summary(state: str | None, flag_set: dict | None, stats: dict | None) -> dict:
    """エクスポートの成功 payload に載せる `curation` の要約（両経路で同じ形）。
    件数は実際に除外を行った `applied` のときだけ数字で、それ以外は None。"""
    return {"state": state,
            "wrong_excluded": stats["wrong_excluded"] if stats else None,
            "suspect": stats["suspect"] if stats else None,
            "orphaned": (flag_set or {}).get("orphaned", 0)}


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
