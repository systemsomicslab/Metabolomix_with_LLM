"""dataset_feature_search の純ロジック: 目標の正規化・理論 m/z・規則 hash（spec 2026-10-19 §4.1）と
目標照合の判定（§5: `classify_identity` / `match_target`）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。統計の結合は同じ
モジュールに後から足す。

**判定は強度・ppm・RT の近さ・並び順で 1 件を選ばない**（ユーザー判断 2026-10-19）。
付加体の層と証拠の水準は規則による決定的な順序で、その最上位の箱に qualified な特徴が
ちょうど 1 つあるときだけ ``resolved``。2 つ以上なら ``needs_input``
（``multiple_qualified_candidates``）、無ければ ``needs_input``（``no_qualified_candidate``）。

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

import numpy as np

from metabolomix.analysis.feature_bindings import _normalize_adduct
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

__all__ = [
    "CANDIDATE_COLUMNS",
    "DEFAULT_ADDUCT_TIERS",
    "EVIDENCE_LEVELS",
    "MAX_TARGETS",
    "classify_identity",
    "match_target",
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


# ---------- 目標照合の判定（spec §5） ----------

#: 証拠の水準。判定はこの順に見る（`contradicted` は qualified にならないので含めない）。
EVIDENCE_LEVELS = ("sme_match", "sml_match", "unannotated")

#: 候補行の列（この順で TSV にする）。
CANDIDATE_COLUMNS = (
    "target_id", "feature_id", "adduct", "adduct_tier", "theoretical_mz", "mz", "ppm",
    "rt", "delta_rt", "evidence", "sme_rank", "annotated_name", "annotated_inchikey",
    "qualified", "reasons",
)

# 浮動小数の誤差で境界ちょうどの特徴を落とさないための余裕（ppm・分のどちらにも無視できる大きさ）。
_BOUNDARY_EPS = 1e-9


def _key_block(value) -> str | None:
    """InChIKey（または先頭ブロック）から照合に使う先頭 14 文字を取る。無ければ None。"""
    if not isinstance(value, str):
        return None
    return value.strip()[:14].upper() or None


def _identity(candidates, annotation, *, inchikey_block: str, adduct: str) -> dict:
    """`classify_identity` の本体。表示用の名前・InChIKey も一緒に返す。

    SME 候補は**その候補自身の** ``inchikey`` だけを見る。特徴の rank1 の InChIKey
    （``feature_metadata[fid]["inchikey"]``）を混ぜると、別物質の rank2 候補が rank1 の
    識別子を借りて一致してしまう（``feature_bindings._identifier_values`` の癖。spec §2）。
    """
    wanted_adduct = _normalize_adduct(adduct)
    parsed = parse_adduct(adduct)
    wanted_charge = parsed.charge if parsed is not None else None

    reasons: list[str] = []
    first_informative = None
    for candidate in candidates or []:
        block = _key_block(candidate.get("inchikey"))
        if block is None:
            # InChIKey を解決できない候補は、一致とも矛盾とも言えない（評価できない）。
            continue
        if first_informative is None:
            first_informative = candidate
        if block != inchikey_block:
            reasons.append("identifier_mismatch")
            continue
        if _normalize_adduct(candidate.get("adduct")) != wanted_adduct:
            reasons.append("adduct_mismatch")
            continue
        charge = candidate.get("charge")
        # 符号の書き方（-1 / 1）に依存しないよう価数で比べる。読めない電荷は一致を確かめられない。
        if wanted_charge is None or not _is_number(charge) or abs(charge) != wanted_charge:
            reasons.append("charge_mismatch")
            continue
        return {"evidence": "sme_match", "reasons": [], "sme_rank": candidate.get("rank"),
                "name": candidate.get("chemical_name"), "inchikey": candidate.get("inchikey")}

    if first_informative is not None:
        # 同定（MS/MS を見た SME）が有って矛盾する候補は落とす（ユーザー判断 2026-10-19）。
        return {"evidence": "contradicted", "reasons": list(dict.fromkeys(reasons)),
                "sme_rank": None, "name": first_informative.get("chemical_name"),
                "inchikey": first_informative.get("inchikey")}

    # SME が無い、または識別子を持つ SME 候補が無い → SML へ進む。
    prefix = ["sme_identity_unresolved"] if candidates else []
    annotation = annotation or {}
    if annotation.get("ambiguous"):
        keys = [key for key in annotation.get("inchikeys") or [] if key]
        if any(_key_block(key) == inchikey_block for key in keys):
            return {"evidence": "sml_match", "reasons": prefix + ["sml_ambiguous"],
                    "sme_rank": None, "name": annotation.get("name"),
                    "inchikey": ";".join(keys)}
    elif _key_block(annotation.get("inchikey")) == inchikey_block:
        return {"evidence": "sml_match", "reasons": prefix, "sme_rank": None,
                "name": annotation.get("name"), "inchikey": annotation.get("inchikey")}
    # 注釈なし、または SML が別物質。SML は MS1 照合だけで異性体を区別できないので、
    # 別物質の SML は矛盾にしない（spec §5）。
    return {"evidence": "unannotated", "reasons": prefix, "sme_rank": None,
            "name": None, "inchikey": None}


def classify_identity(candidates: list[dict], annotation: dict | None, *,
                      inchikey_block: str, adduct: str) -> tuple[str, list[str], int | None]:
    """1 つの（特徴, 付加体）候補の同定の証拠を決める（spec §5 の 3）。

    戻り値は ``(evidence, reasons, sme_rank)``。``evidence`` は ``sme_match`` /
    ``sml_match`` / ``unannotated`` / ``contradicted``。``sme_rank`` は ``sme_match`` の
    ときだけ一致した候補の rank、それ以外は None。

    - SME 候補（rank を問わない）のどれかが InChIKey 先頭ブロック・付加体（正規化して比較）・
      価数のすべてで一致 → ``sme_match``。
    - InChIKey を持つ SME 候補がどれも一致しない → ``contradicted``（理由は出た順に一意化）。
    - InChIKey を持つ SME 候補が無ければ SME は無いものとして SML を見る（SME 候補が
      あったなら ``sme_identity_unresolved``）。SML の InChIKey が一致 → ``sml_match``
      （曖昧な注釈はどれか 1 つの一致で足り、``sml_ambiguous`` を付ける）。
    - それ以外 → ``unannotated``。
    """
    result = _identity(candidates, annotation, inchikey_block=inchikey_block, adduct=adduct)
    return result["evidence"], result["reasons"], result["sme_rank"]


def _empty_counts() -> dict:
    return {"mass": 0, "rt": 0, "qualified": 0,
            "by_evidence": {level: 0 for level in EVIDENCE_LEVELS}}


def _feature_mz_array(ds) -> np.ndarray:
    """``ds.feature_ids`` の並びの m/z 配列（読めない値は NaN。NaN はどの窓にも入らない）。"""
    values = []
    for fid in ds.feature_ids:
        mz = (ds.feature_metadata.get(fid) or {}).get("mz")
        values.append(float(mz) if _is_number(mz) else math.nan)
    return np.asarray(values, dtype=float)


def match_target(ds, target: dict, *, max_results: int) -> dict:
    """正規化済みの目標 1 件をデータセットの特徴と照合し、判定と候補の全件を返す（spec §5）。

    候補の単位は（特徴, 付加体）。全層のどれかの付加体の質量窓に入った特徴を、RT で落ちた
    ものも理由つきで行にする。判定は「qualified がある最初の層」→「その層で qualified がある
    最初の証拠水準」の箱の中で qualified を**特徴 ID で一意化して**数え、1 つなら ``resolved``、
    2 つ以上なら ``needs_input``（``multiple_qualified_candidates``）。強度・ppm・RT の近さ・
    並び順では選ばない。``resolved`` の特徴がその箱に複数の付加体で現れたときは、層に並べた
    順で最初の付加体を ``selected_adduct`` にする（近さではなく規則の順）。
    ``needs_input`` で複数のときは ``adduct_tier_used`` だけ埋め、``selected_*`` は None。

    行は（qualified が先、層の昇順、|ppm| の昇順）に並べ、``max_results`` 件で打ち切って
    残りの件数を ``truncated`` に書く。``counts`` は特徴 ID の一意な件数（打ち切り前）。
    ``ppm`` は符号つき（観測 − 理論）、``delta_rt`` は RT − rt_min。数値は丸めない
    （丸めはツール層の仕事）。
    """
    out = {
        "target_id": target["target_id"],
        "status": None,
        "reason": None,
        "selected_feature_id": None,
        "selected_adduct": None,
        "selected_evidence": None,
        "adduct_tier_used": None,
        "counts": _empty_counts(),
        "theoretical_mz": dict(target["theoretical_mz"]),
        "candidates": [],
        "truncated": 0,
    }
    if not target.get("applicable", True):
        out.update(status="not_applicable", reason="polarity_mismatch")
        return out

    block = target["inchikey_block"]
    ppm_tol = target["ppm"]
    rt_min = target["rt_min"]
    rt_tol = target["rt_tol_min"]
    feature_ids = list(ds.feature_ids)
    mz_array = _feature_mz_array(ds)

    rows: list[dict] = []
    rt_ok_of: dict[str, bool] = {}
    for tier_no, tier in enumerate(target["adduct_tiers"], start=1):
        for adduct in tier:
            theo = target["theoretical_mz"][adduct]
            with np.errstate(invalid="ignore"):
                hits = np.nonzero(np.abs((mz_array - theo) / theo * 1e6)
                                  <= ppm_tol + _BOUNDARY_EPS)[0]
            for index in hits:
                fid = feature_ids[index]
                mz = float(mz_array[index])
                rt = (ds.feature_metadata.get(fid) or {}).get("rt")
                if rt_min is None:
                    rt_ok, delta_rt = True, None
                elif _is_number(rt):
                    delta_rt = rt - rt_min
                    rt_ok = abs(delta_rt) <= rt_tol + _BOUNDARY_EPS
                else:
                    rt_ok, delta_rt = False, None
                rt_ok_of[fid] = rt_ok

                identity = _identity(ds.feature_candidates.get(fid) or [],
                                     ds.feature_annotations.get(fid),
                                     inchikey_block=block, adduct=adduct)
                evidence = identity["evidence"]
                reasons = [] if rt_ok else ["rt_outside_tolerance"]
                reasons += identity["reasons"]
                if evidence == "contradicted":
                    qualified = False
                elif evidence == "unannotated" and rt_min is None:
                    # 質量だけでは同定と言わない（feature_bindings の mass_rt が RT を必須にするのと同じ）
                    reasons.append("mass_only_not_identified")
                    qualified = False
                else:
                    qualified = rt_ok
                rows.append({
                    "target_id": target["target_id"],
                    "feature_id": fid,
                    "adduct": adduct,
                    "adduct_tier": tier_no,
                    "theoretical_mz": theo,
                    "mz": mz,
                    "ppm": (mz - theo) / theo * 1e6,
                    "rt": rt,
                    "delta_rt": delta_rt,
                    "evidence": evidence,
                    "sme_rank": identity["sme_rank"],
                    "annotated_name": identity["name"],
                    "annotated_inchikey": identity["inchikey"],
                    "qualified": qualified,
                    "reasons": reasons,
                })

    qualified_rows = [row for row in rows if row["qualified"]]
    counts = out["counts"]
    counts["mass"] = len({row["feature_id"] for row in rows})
    counts["rt"] = sum(1 for ok in rt_ok_of.values() if ok)
    counts["qualified"] = len({row["feature_id"] for row in qualified_rows})
    for level in EVIDENCE_LEVELS:
        counts["by_evidence"][level] = len(
            {row["feature_id"] for row in qualified_rows if row["evidence"] == level})

    if not qualified_rows:
        out.update(status="needs_input", reason="no_qualified_candidate")
    else:
        tier_used = min(row["adduct_tier"] for row in qualified_rows)
        in_tier = [row for row in qualified_rows if row["adduct_tier"] == tier_used]
        level_used = next(level for level in EVIDENCE_LEVELS
                          if any(row["evidence"] == level for row in in_tier))
        # 行は層の中で付加体を並べた順に作ってあるので、bucket の順は層の順を保つ
        bucket = [row for row in in_tier if row["evidence"] == level_used]
        winners = list(dict.fromkeys(row["feature_id"] for row in bucket))
        out["adduct_tier_used"] = tier_used
        if len(winners) == 1:
            out.update(status="resolved", selected_feature_id=winners[0],
                       selected_adduct=bucket[0]["adduct"], selected_evidence=level_used)
        else:
            out.update(status="needs_input", reason="multiple_qualified_candidates")
        for row in qualified_rows:
            if row["adduct_tier"] != tier_used:
                row["reasons"].append("lower_adduct_tier")
            elif row["evidence"] != level_used:
                row["reasons"].append("lower_evidence_tier")

    for row in rows:
        row["reasons"] = ";".join(dict.fromkeys(row["reasons"]))
    rows.sort(key=lambda row: (not row["qualified"], row["adduct_tier"], abs(row["ppm"])))
    limit = max(0, int(max_results))
    out["candidates"] = rows[:limit]
    out["truncated"] = len(rows) - len(out["candidates"])
    return out
