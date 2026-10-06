"""参照ライブラリの正規化レコード（leaf）。

`.dbs` と `.msp` の 2 つの入口が吐く形を 1 つに揃える。これが store のスキーマで
あり、採点エンジンの入力契約でもある。ここは stdlib しか import しない。
"""
from __future__ import annotations

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
