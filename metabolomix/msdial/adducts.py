"""MS-DIAL のアダクト表記（`[2M+FA-H]-` など）の解析と、m/z ⇄ 中性質量の換算（純関数）。

m/z = (n_mol·M + shift) / charge。shift は付加・脱離する原子団の質量の符号付き和から、
正イオンなら電子 charge 個ぶんを引き、負イオンなら足したもの（既存 `ADDUCT_SHIFTS` と同じ定義）。
略記: FA = ギ酸 CH2O2、Hac = 酢酸 C2H4O2、TFA = トリフルオロ酢酸 C2HF3O2、ACN = C2H3N。
deps: msdial.peak_verification（元素質量と組成式の解析）。
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from metabolomix.msdial.peak_verification import ELECTRON_MASS, ELEMENT_MASSES, parse_formula

_ALIASES = {"FA": "CH2O2", "Hac": "C2H4O2", "TFA": "C2HF3O2", "ACN": "C2H3N"}
_ADDUCT_RE = re.compile(r"^\[(\d*)M((?:[+-]\d*[A-Z][A-Za-z0-9]*)*)\](\d*)([+-])$")
_TERM_RE = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)")


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
