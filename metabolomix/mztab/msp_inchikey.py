"""研究室 MSP の NAME → INCHIKEY 表（spec 2026-10-10 msp-inchikey-design §3.1）。

`library_mode="msp_only"` の Console 実行では、mzTab-M の同定は
`database_identifier = MspDB:<MSP の NAME>` で出て、InChIKey の列は無い。RDKit が
使えない機械（Windows のアプリケーション制御がネイティブ DLL を拒否する）でも
InChIKey に届くよう、Console に渡した MSP から名前で引く。

読むのは各エントリの `NAME:` と `INCHIKEY:` 行だけで、ピーク行は解さない（研究室の
MSP は 1 GB 級・数百万エントリで、全フィールドを解すと遅い）。同じ名前に別の
InChIKey があるものは、名前では決められないので表から外して `conflicting` に入れる。

stdlib と `metabolomix.mztab.identity` だけを import する。`mztab/loading.py`
（pipeline の worker も通る）から呼ぶので、`mcp_core` を引き込む `library.store` の
SQLite キャッシュは使わない。代わりに同じプロセスでは (パス, サイズ, 更新時刻) を
鍵に表を使い回す。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from metabolomix.mztab.identity import _INCHIKEY_RE

__all__ = ["MSP_DB_PREFIX", "MspInchikeyTable", "load_msp_inchikey_table", "lookup"]

#: MS-DIAL が mzTab-M の `database_identifier` に書く MSP 由来の接頭辞。
MSP_DB_PREFIX = "MspDB:"

_BOM = b"\xef\xbb\xbf"


@dataclass(frozen=True)
class MspInchikeyTable:
    """名前 → InChIKey。`file_name` は来歴用のファイル名だけ（置き場所は持たない）。"""
    file_name: str
    by_name: dict[str, str]
    conflicting: frozenset[str]


#: (絶対パス, サイズ, 更新時刻 ns) → 表。陽・陰の 2 本を交互に読んでも追い出さない。
_memo: dict[tuple[str, int, int], MspInchikeyTable] = {}


def load_msp_inchikey_table(path: str | Path) -> MspInchikeyTable:
    """MSP を読み、名前 → InChIKey 表を返す。読めなければ OSError をそのまま投げる。"""
    p = Path(path).resolve()
    st = p.stat()
    key = (str(p), st.st_size, st.st_mtime_ns)
    hit = _memo.get(key)
    if hit is None:
        hit = _scan(p)
        _memo[key] = hit
    return hit


def _scan(path: Path) -> MspInchikeyTable:
    by_name: dict[str, str] = {}
    conflicting: set[str] = set()
    name: str | None = None
    with path.open("rb") as f:
        for i, raw in enumerate(f):
            if i == 0 and raw.startswith(_BOM):
                raw = raw[len(_BOM):]
            line = raw.strip()
            if not line:
                # 空行がエントリの区切り。名前を次のエントリへ持ち越さない。
                name = None
                continue
            head = line[:9].upper()
            if head.startswith(b"NAME:"):
                name = line[5:].decode("utf-8", "replace").strip() or None
            elif head == b"INCHIKEY:" and name is not None:
                inchikey = line[9:].decode("ascii", "replace").strip()
                if not _INCHIKEY_RE.match(inchikey) or name in conflicting:
                    continue
                known = by_name.get(name)
                if known is None:
                    by_name[name] = inchikey
                elif known != inchikey:
                    del by_name[name]
                    conflicting.add(name)
    return MspInchikeyTable(path.name, by_name, frozenset(conflicting))


def lookup(table: MspInchikeyTable, database_identifier: str | None) -> tuple[str | None, str]:
    """`MspDB:<名前>` を表で引く。

    戻り値の第 2 要素: `resolved`（付いた）/ `conflicting`（同名に別キーがあり決めない）/
    `unmatched`（MSP 由来なのに名前が表に無い）/ `not_msp`（MSP 由来の同定ではない）。
    名前は先頭の接頭辞だけを外して引く（`LPC 18:0` のようにコロンを含みうる）。
    """
    text = (database_identifier or "").strip()
    if not text.startswith(MSP_DB_PREFIX):
        return None, "not_msp"
    name = text[len(MSP_DB_PREFIX):].strip()
    inchikey = table.by_name.get(name)
    if inchikey is not None:
        return inchikey, "resolved"
    if name in table.conflicting:
        return None, "conflicting"
    return None, "unmatched"
