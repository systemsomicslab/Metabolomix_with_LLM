"""極性の表記の正規化（stdlib だけの leaf）。

受け口ごとに `pos` を許したり許さなかったりすると、読めない値が黙って別の極性として扱われる
（2026-11-13 の事実確認 A11: library の `ion_mode` と `console_method_candidates` は `pos` を
検証していなかった）。極性を受けるところはここを通し、読めない値は誤りにする。
"""
from __future__ import annotations

__all__ = ["POLARITIES", "canonical_polarity", "normalize_polarity"]

POLARITIES = ("positive", "negative")

_ALIASES = {
    "positive": "positive", "pos": "positive", "+": "positive", "p": "positive",
    "negative": "negative", "neg": "negative", "-": "negative", "n": "negative",
}


def normalize_polarity(value, *, field: str = "polarity") -> str:
    """`positive` / `pos` / `+` / `negative` / `neg` / `-`（大文字小文字・前後の空白を問わない）を
    `"positive"` / `"negative"` にする。それ以外（None・空・別の語）は `ValueError`。"""
    if isinstance(value, str):
        key = value.strip().lower()
        if key in _ALIASES:
            return _ALIASES[key]
    raise ValueError(
        f'{field} は "positive" / "negative"（"pos" / "neg" / "+" / "-" も可）のいずれかです: {value!r}')


def canonical_polarity(value):
    """読める表記（`normalize_polarity` と同じ別名）は `"positive"` / `"negative"` に、それ以外
    （None・読めない値）はそのまま返す。受け口の検証が、読めない値を元の綴りで名指しして断る。"""
    try:
        return normalize_polarity(value)
    except ValueError:
        return value
