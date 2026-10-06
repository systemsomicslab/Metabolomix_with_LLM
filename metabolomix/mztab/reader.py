"""mzTab-M 2.0 テキスト形式パーサ（純関数）。

セクション識別規則:
  - MTD 行: メタデータ。tab 区切り [prefix, key, value]。
  - SMH / SFH / SEH 行: ヘッダー（列名）。mzTab-M 2.0 の正規形で、実 MS-DIAL
    出力もこれ。SMH→SML / SFH→SMF / SEH→SME と読み替える。
  - ヘッダー行が無いまま SML/SMF/SME が来た場合は、その最初の行を列名として扱う
    （合成 fixture などヘッダーをデータ接頭辞で書く流儀への後方互換）。**この
    フォールバックは、正規の接頭辞を取りこぼすと先頭データ行を静かに 1 行
    飲み込む**ので、`_HEADER_PREFIX_MAP` に漏れを作らないこと。
  - COM 行: コメント。無視する。
  - "null" / 空文字列 → None に正規化する。
  - SME 行の末尾空欄（既知の MS-DIAL 問題）を除去しwarningを記録する。

spec §9 参照。
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# ヘッダー接頭辞 → データ接頭辞の対応（mzTab-M 2.0 の正規形）
_HEADER_PREFIX_MAP = {"SEH": "SME", "SFH": "SMF", "SMH": "SML"}
# データ行として認識する接頭辞
_DATA_PREFIXES = frozenset({"SML", "SMF", "SME"})
# abundance 列名パターン（命名揺れに対応）
_ABUNDANCE_RE = re.compile(r"abundance_assay\[(\d+)\]", re.IGNORECASE)
# assay MTD キーパターン
_ASSAY_RE = re.compile(r"^assay\[(\d+)\]-ms_run_ref$")


def _normalize_value(v: str | None) -> str | None:
    if v is None:
        return None
    s = v.strip()
    return None if s in ("null", "") else s


def parse_mztab(path: str | Path) -> dict:
    """mzTab-M 2.0 ファイルを解析し構造化辞書を返す。

    戻り値:
        {
          "metadata": {key: value},
          "sections": {
            "SML": {"header": [...], "rows": [...], "warnings": [...]},
            "SMF": {"header": [...], "rows": [...], "warnings": [...]},
            "SME": {"header": [...], "rows": [...], "warnings": [...]},
          },
          "warnings": [...],  # ファイル全体のwarning
        }
    """
    metadata: dict[str, str] = {}
    sections: dict[str, dict] = {}
    global_warnings: list[str] = []

    # セクション別の末尾空列除去タリー: {prefix: {rows, columns, first, last}}
    trailing_tally: dict[str, dict[str, int]] = {}

    with open(path, encoding="utf-8") as f:
        for lineno, raw_line in enumerate(f, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            prefix = parts[0].strip()

            if prefix == "MTD":
                if len(parts) >= 3:
                    metadata[parts[1].strip()] = parts[2].strip()
                elif len(parts) == 2:
                    metadata[parts[1].strip()] = ""
                continue

            if prefix == "COM":
                continue

            # ヘッダー行（mzTab-M 2.0 の正規の接頭辞: SMH / SFH / SEH）
            canonical = _HEADER_PREFIX_MAP.get(prefix)
            if canonical:
                sections.setdefault(canonical, {"header": None, "rows": [], "warnings": []})
                sections[canonical]["header"] = [c.strip() for c in parts[1:]]
                continue

            # データ行
            if prefix in _DATA_PREFIXES:
                sec = sections.setdefault(prefix, {"header": None, "rows": [], "warnings": []})
                # ヘッダー未定義の場合、この行がヘッダー
                if sec["header"] is None:
                    sec["header"] = [c.strip() for c in parts[1:]]
                    continue
                header = sec["header"]
                values = [c.strip() for c in parts[1:]]
                # 末尾空欄の正規化（既知の MS-DIAL SME 問題）。
                # 警告は 1 件ずつ積まずに集約する（下のループ後に 1 件だけ出す）。
                # 良性のこの警告を行×列ぶん積むと実データで数百件になり、
                # 表示上の件数制限に当たって重要な警告が見えなくなる。
                removed = 0
                while values and values[-1] == "":
                    values.pop()
                    removed += 1
                if removed:
                    tally = trailing_tally.setdefault(
                        prefix, {"rows": 0, "columns": 0, "first": lineno, "last": lineno})
                    tally["rows"] += 1
                    tally["columns"] += removed
                    tally["last"] = lineno
                row = {}
                for i, col in enumerate(header):
                    raw = values[i] if i < len(values) else None
                    row[col] = _normalize_value(raw)
                sec["rows"].append(row)

    for prefix, tally in trailing_tally.items():
        sections[prefix]["warnings"].append(
            f"trailing empty column removed: {tally['columns']} columns in "
            f"{tally['rows']} rows (L{tally['first']}..L{tally['last']})"
        )

    return {"metadata": metadata, "sections": sections, "warnings": global_warnings}


def get_assay_count(metadata: dict) -> int:
    """MTD から assay 数を返す。"""
    return sum(1 for k in metadata if _ASSAY_RE.match(k))


def extract_abundance_matrix(
    parse_result: dict,
) -> tuple[np.ndarray, list[str], list[str]]:
    """SMF 行列から abundance 値を numpy 配列に変換する。

    戻り値: (matrix, sample_names, feature_ids)
      matrix shape: (n_features, n_samples)、欠損は np.nan。
      sample_names: abundance 列名をアサイ番号昇順に並べた list。
      feature_ids: SMF_ID を行順に並べた list。
    """
    smf = parse_result["sections"].get("SMF", {})
    header = smf.get("header") or []
    rows = smf.get("rows") or []

    # abundance 列を番号昇順でソート
    abundance_cols = sorted(
        (col for col in header if _ABUNDANCE_RE.search(col)),
        key=lambda c: int(_ABUNDANCE_RE.search(c).group(1)),
    )
    feature_ids = [r.get("SMF_ID") or str(i) for i, r in enumerate(rows)]
    n_feat, n_samp = len(rows), len(abundance_cols)
    matrix = np.full((n_feat, n_samp), np.nan, dtype=float)
    for i, row in enumerate(rows):
        for j, col in enumerate(abundance_cols):
            val = row.get(col)
            if val is not None:
                try:
                    matrix[i, j] = float(val)
                except (ValueError, TypeError):
                    pass

    return matrix, abundance_cols, feature_ids
