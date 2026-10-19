"""dataset_feature_search の土台: 目標の正規化・理論 m/z・規則 hash（spec 2026-10-19 §4.1）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。判定
（`classify_identity` / `match_target`）と統計の結合は同じモジュールに後から足す。

**既定値を省いた呼び出しと明示した呼び出しを同じ規則として扱う。** `rule_hash` は
正規化後の目標（既定の ppm・RT 許容幅・付加体の層を埋めたもの）から作るので、
`ppm=10` を書いても省いても同じ hash になる。result_id・dataset_hash は入れない
（それらは規則ではなく、規則を当てた対象の出所）。

**丸めない。** 正規化後の目標は理論 m/z を full float のまま持つ。丸めはツール層の仕事。
"""
from __future__ import annotations

import copy
import math
import re

from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

__all__ = [
    "DEFAULT_ADDUCT_TIERS",
    "MAX_TARGETS",
    "normalize_targets",
    "rule_hash",
    "theoretical_mz",
]

# 付加体の層の既定（spec §4.1。手順書 v5 工程 7-1 の第 1 段と同じ集合）。
DEFAULT_ADDUCT_TIERS: dict[str, list[list[str]]] = {
    "positive": [
        ["[M+H]+", "[M]+"],
        ["[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M+H-NH3]+", "[2M+H]+", "[M+2H]2+"],
    ],
    "negative": [
        ["[M-H]-"],
        ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[2M-H]-", "[M+Na-2H]-"],
    ],
}

MAX_TARGETS = 50

# 数値は float に揃える（canonical_hash は 10 と 10.0 を別の JSON に書くため）。
DEFAULT_PPM = 10.0
DEFAULT_RT_TOL_MIN = 0.3

_ALLOWED_KEYS = frozenset({
    "target_id", "inchikey", "formula", "expected_mz", "adduct", "name",
    "ppm", "rt_min", "rt_tol_min", "adduct_tiers", "polarity",
})
_POLARITIES = ("positive", "negative")
_SIGN_OF = {"positive": "+", "negative": "-"}
_INCHIKEY_BLOCK_RE = re.compile(r"[A-Z]{14}")
# parse_formula は読めない部分を黙って捨てるので（"Ca(OH)2" が Ca+O+H2 になる等）、
# 文字列全体が「元素記号+個数」の連なりであることを先に確かめる。
_FORMULA_RE = re.compile(r"(?:[A-Z][a-z]?\d*)+")


def _formula_mass(formula) -> float | None:
    """分子式の中性モノアイソトピック質量。全体が読めない・未知元素なら None。"""
    if not isinstance(formula, str):
        return None
    text = formula.strip()
    if not _FORMULA_RE.fullmatch(text):
        return None
    try:
        return monoisotopic_mass(parse_formula(text))
    except (ValueError, KeyError):
        return None


def theoretical_mz(formula: str, adduct: str) -> float | None:
    """中性分子式と付加体名から理論 m/z を返す。読めない（分子式・付加体・元素）ときは None。"""
    parsed = parse_adduct(adduct)
    if parsed is None:
        return None
    neutral = _formula_mass(formula)
    if neutral is None:
        return None
    return mz_from_neutral(neutral, parsed)


def _invalid(target_id, field: str, message: str) -> DomainError:
    return DomainError(
        "FEATURE_TARGET_INVALID", message, {"target_id": target_id, "field": field})


def _is_number(value) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value))


def _positive_number(target_id, raw: dict, field: str, default):
    if field not in raw or raw[field] is None:
        return default
    value = raw[field]
    if not _is_number(value) or value <= 0:
        raise _invalid(target_id, field, f"{field} は正の数で指定する: {value!r}")
    return float(value)


def _clean_tiers(target_id, tiers) -> list[list[str]]:
    if (not isinstance(tiers, (list, tuple)) or not tiers
            or not all(isinstance(tier, (list, tuple)) and tier for tier in tiers)
            or not all(isinstance(a, str) and a.strip() for tier in tiers for a in tier)):
        raise _invalid(
            target_id, "adduct_tiers",
            "adduct_tiers は付加体名の空でないリストのリスト（例 [[\"[M+H]+\"], [\"[M+Na]+\"]]）で指定する")
    return [[a.strip() for a in tier] for tier in tiers]


def _normalize_one(raw, dataset_polarity: str | None):
    """1 件を正規化する。極性が決まらず層も作れないときは None（呼び出し側が束ねて拒否する）。"""
    if not isinstance(raw, dict):
        raise _invalid(None, "targets", f"目標は辞書で指定する: {raw!r}")
    target_id = raw.get("target_id")
    if not isinstance(target_id, str) or not target_id.strip():
        raise _invalid(None, "target_id", "target_id は空でない文字列で指定する")
    target_id = target_id.strip()
    unknown = sorted(set(raw) - _ALLOWED_KEYS, key=str)
    if unknown:
        raise _invalid(
            target_id, str(unknown[0]),
            f"未知のキー: {unknown}（使えるのは {sorted(_ALLOWED_KEYS)}）")

    inchikey = raw.get("inchikey")
    block = inchikey.strip()[:14].upper() if isinstance(inchikey, str) else ""
    if not _INCHIKEY_BLOCK_RE.fullmatch(block):
        raise _invalid(
            target_id, "inchikey",
            "inchikey は必須（先頭 14 文字の英大文字ブロックを照合に使う。全体を渡してよい）")

    name = raw.get("name")
    if name is not None and not isinstance(name, str):
        raise _invalid(target_id, "name", "name は文字列で指定する")

    ppm = _positive_number(target_id, raw, "ppm", DEFAULT_PPM)
    rt_tol_min = _positive_number(target_id, raw, "rt_tol_min", DEFAULT_RT_TOL_MIN)
    rt_min = raw.get("rt_min")
    if rt_min is not None:
        if not _is_number(rt_min) or rt_min < 0:
            raise _invalid(target_id, "rt_min", f"rt_min は 0 以上の数（分）で指定する: {rt_min!r}")
        rt_min = float(rt_min)

    polarity = raw.get("polarity")
    if polarity is not None and polarity not in _POLARITIES:
        raise _invalid(target_id, "polarity", f"polarity は positive / negative: {polarity!r}")

    formula = raw.get("formula")
    expected_mz = raw.get("expected_mz")
    adduct = raw.get("adduct")
    has_formula = formula is not None
    if has_formula and expected_mz is not None:
        raise _invalid(target_id, "expected_mz",
                       "formula と expected_mz は同時に指定できない（どちらか一方）")
    if has_formula and adduct is not None:
        raise _invalid(target_id, "adduct",
                       "adduct は expected_mz と組で使う（formula には adduct_tiers を使う）")

    explicit_tiers = raw.get("adduct_tiers")
    resolved = polarity or dataset_polarity

    if has_formula:
        if not isinstance(formula, str):
            raise _invalid(target_id, "formula", "formula は文字列で指定する")
        formula = formula.strip()
        if _formula_mass(formula) is None:
            raise _invalid(target_id, "formula", f"分子式を読めない: {formula!r}")
        if explicit_tiers is not None:
            tiers = _clean_tiers(target_id, explicit_tiers)
        elif resolved is not None:
            tiers = copy.deepcopy(DEFAULT_ADDUCT_TIERS[resolved])
        else:
            return None
        mz_by_adduct: dict[str, float] = {}
        for tier in tiers:
            for name_ in tier:
                value = theoretical_mz(formula, name_)
                if value is None:
                    raise _invalid(target_id, "adduct_tiers", f"付加体を読めない: {name_!r}")
                mz_by_adduct[name_] = value
        expected = None
    else:
        if expected_mz is None:
            # adduct だけ渡された目標は、足りない相方（expected_mz）を指す
            raise _invalid(target_id, "expected_mz" if adduct is not None else "formula",
                           "formula か（expected_mz と adduct）のどちらかが必要")
        if not _is_number(expected_mz) or expected_mz <= 0:
            raise _invalid(target_id, "expected_mz", f"expected_mz は正の数: {expected_mz!r}")
        if not isinstance(adduct, str) or not adduct.strip():
            raise _invalid(target_id, "adduct", "expected_mz には adduct が必要")
        if explicit_tiers is not None:
            raise _invalid(target_id, "adduct_tiers",
                           "expected_mz の目標は付加体が 1 つに固定される（adduct_tiers は使えない）")
        expected_mz = float(expected_mz)
        adduct = adduct.strip()
        if parse_adduct(adduct) is None:
            raise _invalid(target_id, "adduct", f"付加体を読めない: {adduct!r}")
        tiers = [[adduct]]
        mz_by_adduct = {adduct: expected_mz}
        formula = None
        expected = expected_mz

    if resolved is not None:
        sign = _SIGN_OF[resolved]
        for tier in tiers:
            for name_ in tier:
                if parse_adduct(name_).polarity != sign:
                    raise _invalid(
                        target_id, "adduct_tiers",
                        f"付加体 {name_} の極性が {resolved} と食い違う")

    return {
        "target_id": target_id,
        "name": name,
        "inchikey_block": block,
        "formula": formula,
        "expected_mz": expected,
        "adduct_tiers": tiers,
        "theoretical_mz": mz_by_adduct,
        "ppm": ppm,
        "rt_min": rt_min,
        "rt_tol_min": rt_tol_min,
        "polarity": resolved,
        "applicable": not (polarity is not None and dataset_polarity is not None
                           and polarity != dataset_polarity),
    }


def normalize_targets(targets: list[dict], *, dataset_polarity: str | None) -> list[dict]:
    """目標のリストを検証して正規化する（既定値を埋め、理論 m/z を全付加体について計算する）。

    検証失敗は ``FEATURE_TARGET_INVALID``（details に ``target_id`` と ``field``）。
    極性が決まらず層も作れない目標があれば ``POLARITY_UNKNOWN``（呼び出し全体のエラー）。
    ``polarity`` が与えられ ``dataset_polarity`` と違う目標は ``applicable=False``。
    """
    if dataset_polarity not in (None, *_POLARITIES):
        raise ValueError(f"dataset_polarity は positive / negative / None: {dataset_polarity!r}")
    if not isinstance(targets, list):
        raise _invalid(None, "targets", "targets はリストで指定する")
    if len(targets) > MAX_TARGETS:
        raise _invalid(None, "targets",
                       f"目標は 1 回 {MAX_TARGETS} 件まで（{len(targets)} 件渡された）")

    normalized: list[dict] = []
    seen: set[str] = set()
    undecided: list[str] = []
    for raw in targets:
        item = _normalize_one(raw, dataset_polarity)
        target_id = raw["target_id"].strip()
        if target_id in seen:
            raise _invalid(target_id, "target_id", f"target_id が重複している: {target_id!r}")
        seen.add(target_id)
        if item is None:
            undecided.append(target_id)
        else:
            normalized.append(item)
    if undecided:
        raise DomainError(
            "POLARITY_UNKNOWN",
            "データセットの極性が不明で、目標にも polarity・adduct_tiers・"
            "（expected_mz と adduct）が無いため付加体の既定の層を決められない",
            {"target_ids": undecided})
    return normalized


def rule_hash(normalized_targets: list[dict]) -> str:
    """正規化後の目標リストから照合規則の hash を作る（result_id・dataset_hash は入れない）。"""
    return canonical_hash(normalized_targets)
