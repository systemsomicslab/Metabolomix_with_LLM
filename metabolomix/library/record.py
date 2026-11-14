"""参照ライブラリの正規化レコード（leaf）。

`.dbs` と `.msp` の 2 つの入口が吐く形を 1 つに揃える。これが store のスキーマで
あり、採点エンジンの入力契約でもある。ここは stdlib しか import しない。
"""
from __future__ import annotations

import re

RECORD_FIELDS = (
    "name", "precursor_mz", "ion_mode", "adduct", "rt",
    "formula", "inchikey", "smiles", "compound_class", "ontology",
    "spectrum", "library_id", "record_index",
)

#: 上流 `CommonEnums.cs` の `enum IonMode { Positive, Negative, Both }`。
#: レコードには整数ではなく文字列で持つ（保存形式に上流の列挙を漏らさない）。
ION_MODES = {0: "positive", 1: "negative", 2: "both"}


def make_record(**kwargs) -> dict:
    """欠けたキーを既定値で埋めた正規化レコードを返す。"""
    unknown = set(kwargs) - set(RECORD_FIELDS)
    if unknown:
        raise ValueError(f"未知のフィールド: {sorted(unknown)}")
    record = {field: None for field in RECORD_FIELDS}
    record["spectrum"] = []
    record.update(kwargs)
    return record


#: InChIKey の形（14-10-1）。第 2 ブロックの 9 文字目が標準（`S`）か非標準（`N`）の印。
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{8}([SN])[A-Z]-[A-Z]$")


def inchikey_standard(value) -> bool | None:
    """標準の InChIKey か（第 2 ブロックの 9 文字目が `S`）。非標準（`N`）なら False。

    値が無い・InChIKey の形でない（先頭ブロックだけ・空白入り）なら None（判定できない）。末尾の
    1 文字はプロトン化の印で、標準かどうかとは関係ない。名前から外部（PubChem など）へ照会はしない
    （研究室のライブラリの名前を外へ送らない）。
    """
    if not isinstance(value, str):
        return None
    m = _INCHIKEY_RE.match(value.strip().upper())
    return None if m is None else m.group(1) == "S"
