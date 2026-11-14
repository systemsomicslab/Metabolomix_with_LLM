"""MS-DIAL のアダクト表記（`[2M+FA-H]-` など）の解析と、m/z ⇄ 中性質量の換算（純関数）。

m/z = (n_mol·M + shift) / charge。shift は付加・脱離する原子団の質量の符号付き和から、
正イオンなら電子 charge 個ぶんを引き、負イオンなら足したもの（既存 `ADDUCT_SHIFTS` と同じ定義）。
略記: FA = ギ酸 CH2O2、Hac = 酢酸 C2H4O2、TFA = トリフルオロ酢酸 C2HF3O2、ACN = C2H3N。

**2 つの表記が同じイオンかは綴りでなく組成で比べる**（`same_adduct`）。MS-DIAL は同じギ酸付加を
`[M+FA-H]-`・`[M+HCOO]-`・`[M+CHO2]-`、酢酸付加を `[M+Hac-H]-`・`[M+CH3COO]-`・`[M+C2H3O2]-` と書き、
mzTab-M は価数 1 を `[M+H]1+` と書く。照合や「自分自身の付加体」の判定で文字列を比べない。
deps: msdial.peak_verification（元素質量と組成式の解析）。
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

from metabolomix.msdial.peak_verification import ELECTRON_MASS, ELEMENT_MASSES, parse_formula

_ALIASES = {"FA": "CH2O2", "Hac": "C2H4O2", "TFA": "C2HF3O2", "ACN": "C2H3N"}
_ADDUCT_RE = re.compile(r"^\[(\d*)M((?:[+-]\d*[A-Z][A-Za-z0-9]*)*)\](\d*)([+-])$")
_TERM_RE = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)")
#: mzTab-M の価数 1 の書き方（`[M+H]1+`）。組成で読めない表記どうしを綴りで比べるときだけ `1` を落とす
#: （`2+` を `+` に潰すと価数が消えるので、落とすのは `1` だけ）。
_UNIT_CHARGE_RE = re.compile(r"\]1([+-])$")


@dataclass(frozen=True)
class Adduct:
    name: str
    n_mol: int
    shift: float
    charge: int
    polarity: str


def _group_mass(token: str) -> float:
    """token の質量を計算する。ALIASES で解決できない場合は parse_formula で試す。
    未知元素は KeyError を発生させる（呼び出し側で例外を掴む）。"""
    counts = parse_formula(_ALIASES.get(token, token))
    return sum(ELEMENT_MASSES[element] * n for element, n in counts.items())


def parse_adduct(name) -> Adduct | None:
    """MS-DIAL のアダクト表記を解析する。

    不正な形式、未知の元素トークン、パース失敗時は None を返す。
    """
    if not name:
        return None
    match = _ADDUCT_RE.match(str(name).strip())
    if not match:
        return None
    n_mol = int(match.group(1) or 1)
    charge = int(match.group(3) or 1)
    polarity = match.group(4)
    shift = 0.0
    try:
        for sign, count, token in _TERM_RE.findall(match.group(2)):
            shift += (1 if sign == "+" else -1) * int(count or 1) * _group_mass(token)
    except (KeyError, ValueError):
        return None
    shift += (-1 if polarity == "+" else 1) * charge * ELECTRON_MASS
    return Adduct(name=str(name).strip(), n_mol=n_mol, shift=shift, charge=charge, polarity=polarity)


def adduct_composition(name) -> tuple[int, int, dict[str, int]] | None:
    """アダクト表記から (n_mol, charge, 付加・脱離する原子の符号付き個数) を返す。

    同位体分布の計算用（電子の出入りは原子数に現れない）。読めなければ None。
    """
    if not name:
        return None
    match = _ADDUCT_RE.match(str(name).strip())
    if not match:
        return None
    delta: dict[str, int] = {}
    try:
        for sign, count, token in _TERM_RE.findall(match.group(2)):
            factor = (1 if sign == "+" else -1) * int(count or 1)
            for element, n in parse_formula(_ALIASES.get(token, token)).items():
                if element not in ELEMENT_MASSES:
                    return None
                delta[element] = delta.get(element, 0) + factor * n
    except ValueError:
        return None
    return int(match.group(1) or 1), int(match.group(3) or 1), {k: v for k, v in delta.items() if v}


@lru_cache(maxsize=1024)
def _identity_key(text: str) -> tuple | None:
    match = _ADDUCT_RE.match(text)
    composition = adduct_composition(text) if match else None
    if composition is None:
        return None
    n_mol, charge, delta = composition
    return n_mol, charge, match.group(4), frozenset(delta.items())


def adduct_key(name) -> tuple | None:
    """アダクトが表すイオンの鍵 `(n_mol, 価数, 極性, 付加・脱離する原子の符号付き個数)`。

    綴りの違い（`[M+FA-H]-` / `[M+HCOO]-` / `[M+CHO2]-`、`[M-H2O-H]-` / `[M-H-H2O]-`、`[M+H]1+` / `[M+H]+`）
    は同じ鍵になる。読めない（形式が違う・未知の元素や略記）なら None。
    """
    if not isinstance(name, str):
        return None
    text = name.strip().replace(" ", "")
    return _identity_key(text) if text else None


def _spelling(name) -> str | None:
    if name is None:
        return None
    text = _UNIT_CHARGE_RE.sub(r"]\1", str(name).strip().replace(" ", "")).casefold()
    return text or None


def same_adduct(a, b) -> bool:
    """2 つのアダクト表記が同じイオンを指すか。

    両方が読めれば組成（`adduct_key`: n_mol・価数・極性・付加・脱離する原子）で比べる。片方でも読めなければ、
    空白を除き価数 1 の `1`（`]1+` / `]1-`）を落とし大文字小文字を無視した綴りで比べる。どちらかが空なら False。
    """
    key_a, key_b = adduct_key(a), adduct_key(b)
    if key_a is not None and key_b is not None:
        return key_a == key_b
    spelling_a = _spelling(a)
    return spelling_a is not None and spelling_a == _spelling(b)


def ion_counts(counts: dict[str, int], adduct) -> dict[str, int] | None:
    """中性分子の元素の数とアダクト表記から、イオンの元素の数（n_mol 倍に付加・脱離する原子を足す）。

    `[2M+H]+` なら中性の 2 倍に H 1 個、`[M+Cl]-` なら Cl 1 個を足す（電子は原子数に現れない）。
    アダクトが読めない、または差し引いて負になる元素があれば None（呼び手は中性の組成式に戻す）。
    """
    composition = adduct_composition(adduct)
    if composition is None:
        return None
    n_mol, _charge, delta = composition
    merged = {element: n * n_mol for element, n in counts.items()}
    for element, n in delta.items():
        merged[element] = merged.get(element, 0) + n
    if any(n < 0 for n in merged.values()):
        return None
    return {element: n for element, n in merged.items() if n}


def mz_from_neutral(neutral: float, adduct: Adduct) -> float:
    """中性質量からアダクト m/z を計算する。

    m/z = (n_mol * neutral + shift) / charge
    """
    return (adduct.n_mol * neutral + adduct.shift) / adduct.charge


def neutral_from_mz(mz: float, adduct: Adduct) -> float:
    """アダクト m/z から中性質量を計算する。

    neutral = (mz * charge - shift) / n_mol
    """
    return (mz * adduct.charge - adduct.shift) / adduct.n_mol
