"""dataset_feature_search の純ロジック: 目標の正規化・理論 m/z・規則 hash（spec 2026-10-19 §4.1）、
目標照合の判定（§5: `classify_identity` / `match_target`）、窓検索と差次的統計の結合
（§4.2・§6: `search_window` / `select_results` / `feature_stats` / `result_summaries`）、
割れたスポットをまとめた行列での検定のやり直し（spec 2026-10-28 §4.3 (c): `merge_plan` /
`merged_results` / `merged_feature_stats`）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

**判定は強度・ppm・RT の近さ・並び順で 1 件を選ばない**（ユーザー判断 2026-10-19）。
付加体の層と証拠の水準は規則による決定的な順序で、その最上位の箱に qualified な特徴が
ちょうど 1 つあるときだけ ``resolved``。2 つ以上なら ``needs_input``
（``multiple_qualified_candidates``）、無ければ ``needs_input``（``no_qualified_candidate``）。

**既定値を省いた呼び出しと明示した呼び出しを同じ規則として扱う。** `rule_hash` は
正規化後の目標（既定の ppm・RT 許容幅・付加体の層を埋めたもの）から作るので、
`ppm=10` を書いても省いても同じ hash になる。result_id・dataset_hash は入れない
（それらは規則ではなく、規則を当てた対象の出所）。表示名 `name` も入れず、目標は
`target_id` の順に並べてから作る（改名・並べ替えで規則が変わったように見せない）。

**丸めない。** 正規化後の目標は理論 m/z を full float のまま持つ。丸めはツール層の仕事。
"""
from __future__ import annotations

import copy
import math
import re

import numpy as np

from metabolomix.analysis.dataset_analysis import included_sample_indices
from metabolomix.analysis.result_state import is_current
from metabolomix.analysis.stats_results import as_two_group, contrast_ids, split_result_ref
from metabolomix.analysis.mass_accuracy import offset_for
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

__all__ = [
    "CANDIDATE_COLUMNS",
    "DEFAULT_ADDUCT_TIERS",
    "EVIDENCE_LEVELS",
    "MAX_TARGETS",
    "STATS_COLUMNS",
    "STATS_COLUMNS_LABELED",
    "WINDOW_COLUMNS",
    "classify_identity",
    "feature_stats",
    "match_target",
    "merge_map",
    "merge_plan",
    "merged_dataset",
    "merged_feature_stats",
    "merged_results",
    "normalize_targets",
    "result_summaries",
    "retest_merged",
    "rule_hash",
    "search_window",
    "select_results",
    "theoretical_mz",
]

# 付加体の層の既定（spec §4.1。手順書 v5 工程 7-1 の第 1 段と同じ集合に、陰イオンの酢酸付加
# `[M+CH3COO]-` を足した。MS-DIAL の mzTab は `[M+CH3COO]1-` と書き、比較時に `1` を落とす）。
# ギ酸付加は `[M+FA-H]-`（`[M+HCOO]-` と同じ質量なので別に持たない）。
DEFAULT_ADDUCT_TIERS: dict[str, list[list[str]]] = {
    "positive": [
        ["[M+H]+", "[M]+"],
        ["[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M+H-NH3]+", "[2M+H]+", "[M+2H]2+"],
    ],
    "negative": [
        ["[M-H]-"],
        ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[M+CH3COO]-", "[2M-H]-", "[M+Na-2H]-"],
    ],
}

MAX_TARGETS = 50

# 数値は float に揃える（canonical_hash は 10 と 10.0 を別の JSON に書くため）。
DEFAULT_PPM = 10.0
DEFAULT_RT_TOL_MIN = 0.3

_ALLOWED_KEYS = frozenset({
    "target_id", "inchikey", "lipid_species", "formula", "expected_mz", "adduct", "name",
    "ppm", "rt_min", "rt_tol_min", "adduct_tiers", "polarity",
})
_POLARITIES = ("positive", "negative")
_SIGN_OF = {"positive": "+", "negative": "-"}
_POLARITY_OF_SIGN = {sign: polarity for polarity, sign in _SIGN_OF.items()}
_INCHIKEY_BLOCK_RE = re.compile(r"[A-Z]{14}")

#: `[M+H]1+`（mzTab-M）と `[M+H]+`（利用者の記法）を同じ adduct として扱うための
#: 末尾電荷表記。`1+` / `1-` だけを落とす——`2+` を `+` に潰すと価数が消える。
_ADDUCT_UNIT_CHARGE_RE = re.compile(r"\]1([+-])$")


def _normalize_adduct(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().replace(" ", "")
    return _ADDUCT_UNIT_CHARGE_RE.sub(r"]\1", text).casefold()


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


def _species_target(target_id, raw: dict) -> dict | None:
    """``lipid_species`` を種のレベルに直す。指定が無ければ None（InChIKey の目標）。"""
    from metabolomix.msdial.lipid_identity import lipid_species_of, normalize_lipid_name

    value = raw.get("lipid_species")
    if value is None:
        return None
    if raw.get("inchikey") is not None:
        raise _invalid(target_id, "lipid_species",
                       "inchikey と lipid_species は同時に指定できない（どちらか一方）")
    species = lipid_species_of(value) if isinstance(value, str) else None
    if species is None:
        raise _invalid(
            target_id, "lipid_species",
            f"lipid_species を種のレベルの脂質名として読めない: {value!r}"
            "（例 \"PC 34:4\"・\"PC 18:3_16:1\"。クラスだけの名前は不可）")
    return {"lipid_species": species,
            "lipid_species_input": normalize_lipid_name(value).get("normalized") or species}


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
    species = _species_target(target_id, raw)
    if species is not None:
        block = None
    else:
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

    if polarity is None and (explicit_tiers is not None or not has_formula):
        # 目標が付加体を明示していれば、その符号が目標の極性（polarity の書き忘れを
        # データセットの極性との食い違い＝呼び出し全体の誤りにしない）。符号が
        # 混ざる目標はどちらの極性のデータでも 1 回で照合できないので誤り。
        signs = {parse_adduct(name_).polarity for tier in tiers for name_ in tier}
        if len(signs) > 1:
            field = "adduct_tiers" if explicit_tiers is not None else "adduct"
            raise _invalid(target_id, field,
                           f"付加体の極性が混ざっている: {[a for t in tiers for a in t]}")
        polarity = _POLARITY_OF_SIGN[signs.pop()]
        resolved = polarity

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
        # 種の目標だけが lipid_species* を持つ（InChIKey の目標の辞書と rule_hash を変えない）
        **(species or {}),
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


def normalize_targets(targets: list[dict], *, dataset_polarity: str | None,
                      ppm_offset: dict | None = None) -> list[dict]:
    """目標のリストを検証して正規化する（既定値を埋め、理論 m/z を全付加体について計算する）。

    検証失敗は ``FEATURE_TARGET_INVALID``（details に ``target_id`` と ``field``）。
    極性が決まらず層も作れない目標があれば ``POLARITY_UNKNOWN``（呼び出し全体のエラー）。
    目標の極性は ``polarity``、無ければ明示の付加体（``adduct`` か ``adduct_tiers``）の
    符号（混ざっていれば ``FEATURE_TARGET_INVALID``）、それも無ければ ``dataset_polarity``。
    目標の極性が ``dataset_polarity`` と違う目標は ``applicable=False``（→ ``not_applicable``）。

    ``ppm_offset``（正規化済み。`mass_accuracy.normalize_ppm_offset`）を渡したときだけ、各目標に
    ``ppm_offset``（付加体 → その理論 m/z の帯の中心のずれ）を足す（rule_hash に入る）。渡さなければ
    目標の辞書は今と同じ。
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
            if ppm_offset is not None:
                item["ppm_offset"] = {adduct: offset_for(mz, ppm_offset)
                                      for adduct, mz in item["theoretical_mz"].items()}
            normalized.append(item)
    if undecided:
        raise DomainError(
            "POLARITY_UNKNOWN",
            "データセットの極性が不明で、目標にも polarity・adduct_tiers・"
            "（expected_mz と adduct）が無いため付加体の既定の層を決められない",
            {"target_ids": undecided})
    return normalized


#: rule_hash に入れない欄。``name`` は読み手のためのラベルで照合に使わない
#: （改名しただけで「規則が変わった」ように見せない）。``lipid_species_input`` は
#: 照合に使う ``lipid_species``（種のレベル）の元の綴りで、``PC 18:3_16:1`` と ``PC 34:4`` は同じ規則。
_RULE_HASH_EXCLUDED = frozenset({"name", "lipid_species_input"})


def rule_hash(normalized_targets: list[dict], *, merge_groups_hash: str | None = None) -> str:
    """正規化後の目標リストから照合規則の hash を作る（result_id・dataset_hash は入れない）。

    覆うのは照合を左右する欄だけ: ``target_id``・InChIKey 先頭ブロック（種の目標は
    ``lipid_species``）・分子式か
    expected_mz・付加体の層と理論 m/z・ppm・RT・目標の極性・``applicable``。``name`` は
    入れず、目標は ``target_id`` の順に並べ替えてから作る——表示名の変更や
    ``targets`` の並べ替えでは hash は変わらない。

    ``merge_groups_hash``（`merge_split` で判定に使った割れの組の hash）を渡したときだけ、
    ``merge_split`` と組の hash も入れる。渡さなければ従来（v7）と同じ hash。
    """
    canonical = sorted(
        ({k: v for k, v in target.items() if k not in _RULE_HASH_EXCLUDED}
         for target in normalized_targets),
        key=lambda target: target["target_id"])
    if merge_groups_hash is None:
        return canonical_hash(canonical)
    return canonical_hash({"targets": canonical, "merge_split": True,
                           "merge_groups_hash": merge_groups_hash})


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


def _identity(candidates, annotation, *, inchikey_block: str | None = None,
              lipid_species: str | None = None, adduct: str) -> dict:
    """`classify_identity` の本体。表示用の名前・InChIKey も一緒に返す。

    SME 候補は**その候補自身の** ``inchikey`` だけを見る。特徴の rank1 の InChIKey
    （``feature_metadata[fid]["inchikey"]``）を混ぜると、別物質の rank2 候補が rank1 の
    識別子を借りて一致してしまう（spec §2）。

    ``lipid_species`` を渡すと（種の目標）、識別子は InChIKey の先頭ブロックではなく
    候補自身の ``chemical_name``（SML なら ``name``）を種のレベルに読んだもの。読めない
    名前は InChIKey の無い候補と同じく評価できないものとして飛ばす。
    """
    if (inchikey_block is None) == (lipid_species is None):
        raise ValueError("inchikey_block と lipid_species のどちらか一方を渡す")
    if lipid_species is not None:
        from metabolomix.msdial.lipid_identity import lipid_species_of

        wanted, key_of, name_key = lipid_species, lipid_species_of, "chemical_name"
    else:
        wanted, key_of, name_key = inchikey_block, _key_block, "inchikey"
    wanted_adduct = _normalize_adduct(adduct)
    parsed = parse_adduct(adduct)
    wanted_charge = parsed.charge if parsed is not None else None

    reasons: list[str] = []
    first_informative = None
    for candidate in candidates or []:
        block = key_of(candidate.get(name_key))
        if block is None:
            # InChIKey（種の目標なら種の名前）を解決できない候補は、一致とも矛盾とも言えない（評価できない）。
            continue
        if first_informative is None:
            first_informative = candidate
        if block != wanted:
            reasons.append("identifier_mismatch")
            continue
        if _normalize_adduct(candidate.get("adduct")) != wanted_adduct:
            reasons.append("adduct_mismatch")
            continue
        charge = candidate.get("charge")
        # 符号の書き方（-1 / 1）に依存しないよう価数で比べる。電荷が読めない（None）候補は、
        # 付加体がすでに一致しているので電荷では落とさない（欠測は反証ではない）。
        # 値があるのに数として読めない・価数が違うときだけ不一致。
        if charge is not None and (
                wanted_charge is None or not _is_number(charge) or abs(charge) != wanted_charge):
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
    if lipid_species is not None:
        # 曖昧な SML は名前を持たない（dataset_state が選ばない）ので、種では照らせない
        if not annotation.get("ambiguous") and key_of(annotation.get("name")) == wanted:
            return {"evidence": "sml_match", "reasons": prefix, "sme_rank": None,
                    "name": annotation.get("name"), "inchikey": annotation.get("inchikey")}
    elif annotation.get("ambiguous"):
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
                      inchikey_block: str | None = None, lipid_species: str | None = None,
                      adduct: str) -> tuple[str, list[str], int | None]:
    """1 つの（特徴, 付加体）候補の同定の証拠を決める（spec §5 の 3）。

    戻り値は ``(evidence, reasons, sme_rank)``。``evidence`` は ``sme_match`` /
    ``sml_match`` / ``unannotated`` / ``contradicted``。``sme_rank`` は ``sme_match`` の
    ときだけ一致した候補の rank、それ以外は None。

    - SME 候補（rank を問わない）のどれかが InChIKey 先頭ブロック・付加体（正規化して比較）・
      価数のすべてで一致 → ``sme_match``。電荷が None の候補は、付加体が一致していれば
      電荷の確認を通す（欠測は反証ではない）。
    - InChIKey を持つ SME 候補がどれも一致しない → ``contradicted``（理由は出た順に一意化）。
    - InChIKey を持つ SME 候補が無ければ SME は無いものとして SML を見る（SME 候補が
      あったなら ``sme_identity_unresolved``）。SML の InChIKey が一致 → ``sml_match``
      （曖昧な注釈はどれか 1 つの一致で足り、``sml_ambiguous`` を付ける）。
    - それ以外 → ``unannotated``。

    ``inchikey_block`` の代わりに ``lipid_species``（種のレベルの脂質名）を渡すと、上の
    「InChIKey」を「候補の名前を種のレベルに読んだもの」に置き換えて同じ規則で決める
    （どちらか一方。両方・どちらも無しは ValueError）。曖昧な SML は名前を持たないので
    種では照らせず ``unannotated``。
    """
    result = _identity(candidates, annotation, inchikey_block=inchikey_block,
                       lipid_species=lipid_species, adduct=adduct)
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


def match_target(ds, target: dict, *, max_results: int,
                 mz_array: np.ndarray | None = None,
                 merged_from: dict[str, list[str]] | None = None) -> dict:
    """正規化済みの目標 1 件をデータセットの特徴と照合し、判定と候補の全件を返す（spec §5）。

    候補の単位は（特徴, 付加体）。全層のどれかの付加体の質量窓に入った特徴を、RT で落ちた
    ものも理由つきで行にする。判定は「qualified がある最初の層」→「その層で qualified がある
    最初の証拠水準」の箱の中で qualified を**特徴 ID で一意化して**数え、1 つなら ``resolved``、
    2 つ以上なら ``needs_input``（``multiple_qualified_candidates``）。強度・ppm・RT の近さ・
    並び順では選ばない。``resolved`` の特徴がその箱に複数の付加体で現れたときは、層に並べた
    順で最初の付加体を ``selected_adduct`` にする（近さではなく規則の順）。
    ``needs_input`` で複数のときは ``adduct_tier_used`` と ``selected_evidence`` を埋め、
    ``selected_feature_id`` / ``selected_adduct`` は None。このときの ``selected_evidence`` は
    選んだ水準ではなく**同点で決着が付かなかった水準**（決定を下した水準）で、qualified が
    1 つも無い（``no_qualified_candidate``）ときは None。

    ``mz_array`` は ``_feature_mz_array(ds)`` を事前に作った配列（``ds.feature_ids`` の並び）。
    多数の目標を照合する呼び出し側が 1 回だけ作って使い回すための引数で、None なら内部で作る。

    行は（qualified が先、層の昇順、|ppm| の昇順）に並べ、``max_results`` 件で打ち切って
    残りの件数を ``truncated`` に書く。``counts`` は特徴 ID の一意な件数（打ち切り前）。
    ``ppm`` は符号つき（観測 − 理論）、``delta_rt`` は RT − rt_min。数値は丸めない
    （丸めはツール層の仕事）。``sme_rank`` は ``sme_match`` の行で一致した候補の rank だが、
    その候補が rank を持たなければ ``sme_match`` でも None になりうる。

    ``merged_from``（`merge_map` の ``{残す ID: 群}``）を渡すと、割れたスポットの組を**1 つの候補**
    として数える（spec 2026-10-28 §4.3 (c)・コントローラの決定 2026-10-29）: 行は従来どおり
    スポットごとに作り（``merged_from`` 欄に群）、判定と ``counts`` は残す ID（群の先頭。受け持ちの
    多い側、同数なら小さい ID）に読み替えて一意化する。組の 2 スポットが同じ箱にあっても 1 つと数え、``resolved`` なら
    ``selected_feature_id`` は残す ID、``merged_from`` にその群（まとめないときは None）。
    組の片方だけが質量の窓に入ったときは、相手のスポットも同じ付加体で行に足す（``reasons`` に
    ``merge_partner``。qualified にせず数にも入れない。選んだ残す ID の行が候補に必ずあるように）。
    ``merged_from`` を渡さない（None）ときは戻り値も行も従来と同じ形（``merged_from`` 欄を足さない）。
    """
    collapse = merged_from is not None
    merged_from = merged_from or {}
    kept_of = {fid: kept for kept, group in merged_from.items() for fid in group}
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
    if collapse:
        out["merged_from"] = None
    if not target.get("applicable", True):
        out.update(status="not_applicable", reason="polarity_mismatch")
        return out

    identifier = ({"lipid_species": target["lipid_species"]} if target.get("lipid_species")
                  else {"inchikey_block": target["inchikey_block"]})
    ppm_tol = target["ppm"]
    centers = target.get("ppm_offset") or {}
    rt_min = target["rt_min"]
    rt_tol = target["rt_tol_min"]
    feature_ids = list(ds.feature_ids)
    if mz_array is None:
        mz_array = _feature_mz_array(ds)

    position = {fid: i for i, fid in enumerate(feature_ids)}

    def make_row(index, tier_no, adduct, theo, *, partner=False):
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
        if not partner:
            rt_ok_of[fid] = rt_ok

        identity = _identity(ds.feature_candidates.get(fid) or [],
                             ds.feature_annotations.get(fid),
                             **identifier, adduct=adduct)
        evidence = identity["evidence"]
        # 割れの相手（質量の窓の外）は判定に入れず、行として見せるだけ
        reasons = ["merge_partner"] if partner else []
        reasons += [] if rt_ok else ["rt_outside_tolerance"]
        reasons += identity["reasons"]
        if partner or evidence == "contradicted":
            qualified = False
        elif evidence == "unannotated" and rt_min is None:
            # 質量だけでは同定と言わない
            reasons.append("mass_only_not_identified")
            qualified = False
        else:
            qualified = rt_ok
        row = {
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
        }
        if collapse:
            row["merged_from"] = merged_from.get(kept_of.get(fid))
        return row

    rows: list[dict] = []
    rt_ok_of: dict[str, bool] = {}
    for tier_no, tier in enumerate(target["adduct_tiers"], start=1):
        for adduct in tier:
            theo = target["theoretical_mz"][adduct]
            with np.errstate(invalid="ignore"):
                hits = np.nonzero(np.abs((mz_array - theo) / theo * 1e6 - centers.get(adduct, 0.0))
                                  <= ppm_tol + _BOUNDARY_EPS)[0]
            for index in hits:
                rows.append(make_row(index, tier_no, adduct, theo))

    if collapse:
        # 組の片方だけが質量の窓に入ったとき、相手（組は ≤15 ppm、窓は ppm）も行に足す。
        # 残す ID が窓の外でも、選んだ ID の行が候補に必ずある（`merge_partner`。qualified にしない）。
        present = {row["feature_id"] for row in rows}
        for row in list(rows):
            for member in merged_from.get(kept_of.get(row["feature_id"]), []):
                if member not in present and member in position:
                    present.add(member)
                    rows.append(make_row(position[member], row["adduct_tier"], row["adduct"],
                                         row["theoretical_mz"], partner=True))

    def unit(fid):
        # 判定と数の単位（割れの組は残す ID で 1 つ）
        return kept_of.get(fid, fid)

    qualified_rows = [row for row in rows if row["qualified"]]
    counts = out["counts"]
    counts["mass"] = len({unit(row["feature_id"]) for row in rows})
    counts["rt"] = len({unit(fid) for fid, ok in rt_ok_of.items() if ok})
    counts["qualified"] = len({unit(row["feature_id"]) for row in qualified_rows})
    for level in EVIDENCE_LEVELS:
        counts["by_evidence"][level] = len(
            {unit(row["feature_id"]) for row in qualified_rows if row["evidence"] == level})

    if not qualified_rows:
        out.update(status="needs_input", reason="no_qualified_candidate")
    else:
        tier_used = min(row["adduct_tier"] for row in qualified_rows)
        in_tier = [row for row in qualified_rows if row["adduct_tier"] == tier_used]
        level_used = next(level for level in EVIDENCE_LEVELS
                          if any(row["evidence"] == level for row in in_tier))
        # 行は層の中で付加体を並べた順に作ってあるので、bucket の順は層の順を保つ
        bucket = [row for row in in_tier if row["evidence"] == level_used]
        winners = list(dict.fromkeys(unit(row["feature_id"]) for row in bucket))
        out["adduct_tier_used"] = tier_used
        if len(winners) == 1:
            out.update(status="resolved", selected_feature_id=winners[0],
                       selected_adduct=bucket[0]["adduct"], selected_evidence=level_used)
            if collapse:
                out["merged_from"] = merged_from.get(winners[0])
        else:
            out.update(status="needs_input", reason="multiple_qualified_candidates",
                       selected_evidence=level_used)
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


# ---------- 窓検索と差次的統計の結合（spec §4.2・§6） ----------

#: 窓検索の行の列（この順で TSV にする）。
WINDOW_COLUMNS = ("feature_id", "mz", "rt", "evidence", "name", "inchikey", "adduct", "n_sme",
                  "ontology", "n_detected", "normalized", "msi_level")

#: 特徴ごとの統計の行の列。
STATS_COLUMNS = ("feature_id", "result_id", "status", "log2fc", "p", "q", "mean_a", "mean_b",
                 "detected_a", "detected_b")

#: `with_group_labels=True` のときの列（既存の列の後ろに群の名前と比較の向きを足す）。
STATS_COLUMNS_LABELED = STATS_COLUMNS + ("comparison", "group_a", "group_b",
                                         "log2fc_test_over_ref")


def comparison_name(result: dict) -> str:
    """比較名。pipeline の比較なら `comparison_id`、そうでなければ `<group_b>_vs_<group_a>`
    （v1 の log2fc は正なら group_b が高いので、group_b が「比較」・group_a が「基準」）。"""
    declared = ((result.get("provenance") or {}).get("comparison") or {}).get("comparison_id")
    return declared or f"{result.get('b')}_vs_{result.get('a')}"


def _rank_key(candidate: dict):
    rank = candidate.get("rank")
    return (0, rank) if _is_number(rank) else (1, 0)


def _window_identity(candidates, annotation) -> dict:
    """窓検索の行に載せる証拠（sme / sml / none）と表示用の名前・InChIKey・付加体。

    SME があれば rank が最小の候補（rank を持たない候補は後回し、同順なら出た順）の
    ``chemical_name`` ・**その候補自身の** ``inchikey`` ・``adduct``。SME が無ければ SML 注釈
    （曖昧なら名前は None、InChIKey は ``;`` 連結）。目標の判定（`_identity`）とは別物で、
    ここは「この特徴には何が付いているか」を並べるだけ。
    """
    candidates = list(candidates or [])
    if candidates:
        top = min(candidates, key=_rank_key)
        return {"evidence": "sme", "name": top.get("chemical_name"),
                "inchikey": top.get("inchikey"), "adduct": top.get("adduct"),
                "n_sme": len(candidates)}
    annotation = annotation or {}
    if annotation:
        if annotation.get("ambiguous"):
            keys = [key for key in annotation.get("inchikeys") or [] if key]
            return {"evidence": "sml", "name": None, "inchikey": ";".join(keys) or None,
                    "adduct": annotation.get("adduct"), "n_sme": 0}
        return {"evidence": "sml", "name": annotation.get("name"),
                "inchikey": annotation.get("inchikey"), "adduct": annotation.get("adduct"),
                "n_sme": 0}
    return {"evidence": "none", "name": None, "inchikey": None, "adduct": None, "n_sme": 0}


def _check_range(low, high, label: str) -> None:
    if low is not None and high is not None and low > high:
        raise DomainError(
            "FEATURE_SEARCH_BAD_RANGE", f"{label} の下限が上限より大きい: {low} > {high}",
            {"field": label, "min": low, "max": high})


def _detected_intensity(ds) -> tuple[np.ndarray, np.ndarray | None]:
    """特徴ごとの（検出セルの値の最大, 検出セル数）。解析に残さない試料の列は使わない。

    使わない試料は前処理と同じ定義（`included_sample_indices`: 実験情報シートの include=false と
    `dataset_exclude` の試料の和集合）。

    検出セルだけを数える（gap-fill は含めない）。検出セルが 1 つも無い特徴の最大は NaN。
    検出状態（``detected_mask``）が無いデータセットは、gap-fill を見分けられないので
    除外していない試料の全セルの最大（NaN は無視）を使い、検出数は None とする。
    """
    n_features = len(ds.feature_ids)
    mask = ds.detected_mask
    if ds.feature_matrix is None:
        return np.full(n_features, np.nan), None
    kept = included_sample_indices(ds)
    keep = np.ones(len(ds.sample_names), dtype=bool)
    if kept is not None:
        keep[:] = False
        keep[kept] = True
    values = np.asarray(ds.feature_matrix, dtype=float)[:, keep]
    finite = np.isfinite(values)
    if mask is None:
        # 検出の印が無い: gap-fill を見分けられないので、除外していない試料の全セルの最大
        # （NaN は無視）。検出数は分からないので None。
        counts = None
        detected = finite
    else:
        detected = np.asarray(mask, dtype=bool)[:, keep] & finite
        counts = detected.sum(axis=1)
    if values.shape[1] == 0:
        return np.full(n_features, np.nan), counts
    best = np.where(detected, values, -np.inf).max(axis=1)
    return np.where(np.isfinite(best), best, np.nan), counts


def search_window(ds, *, feature_ids: list[str] | None, min_mz, max_mz, min_rt, max_rt,
                  max_results: int, sort: str = "mz") -> tuple[list[dict], int]:
    """m/z・RT（分）の窓と特徴 ID の指定で特徴を探す（spec §4.2）。

    条件は与えられたものすべての AND で、境界を含む。``feature_ids`` は None のときだけ
    「条件なし」（空リストは「どの ID でもない」＝ 0 件）。条件が 1 つも無ければ
    ``FEATURE_SEARCH_EMPTY_QUERY``（全特徴の列挙を黙って返さない。``sort="intensity"`` は
    上位 ``max_results`` 件で打ち切るので条件なしでもよい）、下限 > 上限は
    ``FEATURE_SEARCH_BAD_RANGE``。m/z や RT を持たない特徴は、その軸の条件があれば外れる。
    ``ds.excluded_features`` の特徴は結果に出さない。

    ``sort="mz"``（既定）は m/z 昇順（m/z 無しは最後。同値は特徴の並び順）。``sort="intensity"``
    は検出セルの値の最大（除外した試料を除く）の降順で、検出セルが無い特徴は最後（同値は特徴の
    並び順）。この強度は**ピークの高さ（データセットの測定値）の最大**で、EIC のクロマトグラムの
    最大ではない。``max_results`` 件に打ち切り、切り捨てた件数を返す。統計はここでは付けない
    （`feature_stats`）。
    """
    if sort not in ("mz", "intensity"):
        raise DomainError("FEATURE_SEARCH_BAD_SORT", f"sort は mz か intensity: {sort!r}",
                          {"field": "sort"})
    if sort == "mz" and all(v is None for v in (feature_ids, min_mz, max_mz, min_rt, max_rt)):
        raise DomainError(
            "FEATURE_SEARCH_EMPTY_QUERY",
            "feature_ids・min_mz/max_mz・min_rt/max_rt のどれか 1 つは指定する"
            "（条件なしで上位を見るなら sort=\"intensity\"）",
            {})
    _check_range(min_mz, max_mz, "mz")
    _check_range(min_rt, max_rt, "rt")
    wanted = set(feature_ids) if feature_ids is not None else None

    def _inside(value, low, high) -> bool:
        if low is None and high is None:
            return True
        if not _is_number(value):
            return False
        return (low is None or value >= low) and (high is None or value <= high)

    hits = []
    for index, fid in enumerate(ds.feature_ids):
        if fid in ds.excluded_features:
            continue
        if wanted is not None and fid not in wanted:
            continue
        meta = ds.feature_metadata.get(fid) or {}
        mz, rt = meta.get("mz"), meta.get("rt")
        if _inside(mz, min_mz, max_mz) and _inside(rt, min_rt, max_rt):
            hits.append((fid, mz, rt, index))
    best, counts = _detected_intensity(ds)
    if sort == "intensity":
        hits.sort(key=lambda hit: (math.isnan(best[hit[3]]),
                                   0.0 if math.isnan(best[hit[3]]) else -best[hit[3]], hit[3]))
    else:
        hits.sort(key=lambda hit: (not _is_number(hit[1]), hit[1] if _is_number(hit[1]) else 0.0))

    limit = max(0, int(max_results))
    rows = []
    for fid, mz, rt, index in hits[:limit]:
        meta = ds.feature_metadata.get(fid) or {}
        annotation = ds.feature_annotations.get(fid) or {}
        # 同定の標準化の列は解析フォルダのデータセットだけが持つ（mzTab-M は None のまま）
        ident = meta.get("identity") or {}
        rows.append({
            "feature_id": fid,
            "mz": float(mz) if _is_number(mz) else None,
            "rt": float(rt) if _is_number(rt) else None,
            **_window_identity(ds.feature_candidates.get(fid), annotation),
            "ontology": meta.get("ontology") or annotation.get("ontology") or None,
            "n_detected": None if counts is None else int(counts[index]),
            "normalized": ident.get("normalized"),
            "msi_level": ident.get("msi_level"),
        })
    return rows, len(hits) - len(rows)


def _result_id(result: dict):
    return (result.get("provenance") or {}).get("result_id")


_ENGINE_SEPARATOR = "#"


def _engine_views(record: dict, rid: str, contrast_id: str | None = None) -> list[tuple[str, dict]]:
    """`compare` の記録を 2 群の結果（`as_two_group`）に展開した ``[(id, 結果)]``。

    ``contrast_id`` が無ければ交互作用を除く全対比。ID は対比が 1 つなら ``result_id``、複数なら
    ``result_id#contrast_id``。``contrast_id`` があればその対比だけ（未知は `CONTRAST_REQUIRED`、
    交互作用は `CONTRAST_NOT_TWO_GROUP`。どちらも `as_two_group` が候補つきで送出）。
    各結果の ``provenance.result_id`` はその ID に置き換える（`feature_stats` や `result_summaries` が
    これで結果を区別する）。元の記録は ``_source`` に ID だけ残す（割れの再検定が記録を引く）。
    """
    ids = contrast_ids(record)
    wanted = [contrast_id] if contrast_id is not None else contrast_ids(record, two_group_only=True)
    views = []
    for cid in wanted:
        view = as_two_group(record, cid)
        view_id = rid if len(ids) == 1 else f"{rid}{_ENGINE_SEPARATOR}{cid}"
        views.append((view_id, {**view, "provenance": {**view["provenance"], "result_id": view_id},
                                "_source": {"result_id": rid, "contrast_id": cid}}))
    return views


def select_results(ds, result_ids: list[str] | None, *,
                   dataset_hash: str | None = None, contrast_id: str | None = None,
                   single: bool = False) -> dict:
    """統計を引く 2 群比較の結果を選ぶ。``{"results", "stale", "missing"}``。

    省略（None）なら ``ds.results`` のうち 2 群比較のすべて。名指しなら、その ID を引く。2 群比較は
    v1 の ``kind == "two_group"`` と、統計エンジンの ``compare`` の記録（対比ごとに `as_two_group` で
    展開する。ID は対比が 1 つなら ``result_id``、複数なら ``result_id#contrast_id``。交互作用は除く）。

    対比の選び方:

    - 名指しの ID は ``result_id`` か ``result_id#contrast_id``。``contrast_id`` 引数があれば ``result_id`` に
      ``#contrast_id`` を付けたのと同じ（食い違えば `DATASET_BAD_REQUEST`）。
    - 名指しした対比が無い・交互作用なら `CONTRAST_REQUIRED` / `CONTRAST_NOT_TWO_GROUP`（候補つき）。
    - 省略（None）は、最初の概観が回した仮の比較（記録の ``origin == "overview"``）を含めない
      （名指しすれば引ける）。
    - 省略（None）で ``contrast_id`` があれば、それを持つ ``compare`` の記録だけがその対比になる。どの記録も
      持たなければ `CONTRAST_REQUIRED`（候補は全記録の対比）。
    - ``single=True``（単一の対比が要る呼び出し）では、対比が複数の記録を ID だけで名指しすると
      `CONTRAST_REQUIRED`。False なら全対比に展開する。

    今も有効でない結果は ``results`` に入れず ``stale`` に名指しした ID を 1 回だけ挙げる。名指しの ID が無い、
    または 2 群比較でない結果（pipeline の v2 統計の記録 ``kind == "welch"`` など）は ``missing``。
    指定を黙って捨てないよう、次の 2 つも ``missing`` に ``"<ID>: <理由>"`` の 1 行で挙げる:
    ``contrast_id`` を渡したが ``compare`` の記録が無く v1 の 2 群の結果だけがある（v1 は対比を持たないので
    効かない。v1 の結果はそのまま使う。行は ``"contrast_id=<ID>: …"``）、``compare`` の記録の対比が交互作用
    だけで 2 群の結果に読み替えられない（統計は付かない）。

    ``dataset_hash`` は計算済みの ``dataset_fingerprint(ds)``（任意）。
    """
    selected, stale, missing = [], [], []
    done: set[str] = set()
    if result_ids is None:
        candidates = {rid: r for rid, r in ds.results.items() if (r or {}).get("origin") != "overview"}
        records = [(rid, r) for rid, r in candidates.items() if (r or {}).get("kind") == "compare"]
        if contrast_id is not None and records and not any(
                contrast_id in contrast_ids(r) for _, r in records):
            all_ids = sorted({c for _, r in records for c in contrast_ids(r)})
            raise DomainError("CONTRAST_REQUIRED",
                              f"どの比較の結果にも対比 {contrast_id!r} がありません。",
                              {"contrast_id": contrast_id, "contrast_ids": all_ids})
        for rid, result in candidates.items():
            kind = (result or {}).get("kind")
            if kind not in ("two_group", "compare"):
                continue
            if kind == "compare" and contrast_id is not None and contrast_id not in contrast_ids(result):
                continue
            if not is_current(ds, result, dataset_hash=dataset_hash):
                stale.append(rid)
            elif kind == "compare":
                views = _engine_views(result, rid, contrast_id)
                if not views:
                    missing.append(_interaction_only(rid, result))
                selected.extend(view for _, view in views)
            else:
                selected.append(result)
        if contrast_id is not None and not records and any(
                (r or {}).get("kind") == "two_group" for r in ds.results.values()):
            missing.append(f"contrast_id={contrast_id}: dataset_compare の結果が無く、v1 の 2 群の結果は"
                           "対比を持たないので効かない（v1 の結果はそのまま使った）")
        return {"results": selected, "stale": stale, "missing": missing}
    for requested in dict.fromkeys(result_ids):
        found = ds.results.get(requested)
        if (found or {}).get("kind") == "two_group":
            base, cid = requested, None
            record = found
        else:
            base, cid = split_result_ref(requested, contrast_id)
            record = ds.results.get(base)
        if (record or {}).get("kind") == "two_group":
            if is_current(ds, record, dataset_hash=dataset_hash):
                selected.append(record)
            else:
                stale.append(requested)
            continue
        if (record or {}).get("kind") != "compare":
            missing.append(requested)
            continue
        if cid is None and single:
            as_two_group(record)                      # 対比が複数なら CONTRAST_REQUIRED（候補つき）
        if not is_current(ds, record, dataset_hash=dataset_hash):
            if requested not in stale:
                stale.append(requested)
            continue
        views = _engine_views(record, base, cid)
        if not views:
            missing.append(_interaction_only(requested, record))
        for view_id, view in views:
            if view_id not in done:
                done.add(view_id)
                selected.append(view)
    return {"results": selected, "stale": stale, "missing": missing}


def _interaction_only(rid: str, record: dict) -> str:
    """対比が交互作用だけの `compare` の記録の `missing` の行（2 群の結果に読み替えられない）。"""
    return (f"{rid}: 対比が交互作用だけ（{', '.join(contrast_ids(record))}）で 2 群の結果に"
            "読み替えられないので統計を付けない")


def _detected_counter(ds, result: dict):
    """特徴 ID → ``(A 群の検出例数, B 群の検出例数)`` を返す関数を作る。作れなければ None。

    マスクの行は ``ds.feature_ids`` の位置、列は結果の ``samples_a`` / ``samples_b`` を
    ``ds.sample_names`` の位置に直したもの。マスクが無い・形が ``feature_ids`` × ``sample_names``
    と違う・結果の試料名が ``ds.sample_names`` に無い（対応を保証できない）ときは None。
    """
    mask = getattr(ds, "detected_mask", None)
    if mask is None:
        return None
    arr = np.asarray(mask)
    if arr.ndim != 2 or arr.shape != (len(ds.feature_ids), len(ds.sample_names)):
        return None
    position = {name: i for i, name in enumerate(ds.sample_names)}
    columns = []
    for key in ("samples_a", "samples_b"):
        names = result.get(key) or []
        if any(name not in position for name in names):
            return None
        columns.append([position[name] for name in names])
    row_of = {fid: i for i, fid in enumerate(ds.feature_ids)}

    def counter(fid):
        if fid not in row_of:
            return None, None
        row = arr[row_of[fid]]
        return tuple(int(np.count_nonzero(row[cols])) for cols in columns)

    return counter


def _num_or_none(value):
    return float(value) if _is_number(value) else None


def feature_stats(ds, feature_ids: list[str], results: list[dict], *,
                  with_labels: bool = False) -> list[dict]:
    """特徴 × 結果の統計行（特徴が外側、結果が内側）。注釈の有無に関わらず全特徴に付ける。

    ``status``: ``tested``（p が有限）/ ``not_testable``（結果の行はあるが p が None・NaN。
    pipeline run から復元した結果の JSON null を含む）/ ``not_in_matrix``（前処理で落ち、
    結果の行列に無い）。後の 2 つは ``log2fc p q mean_a mean_b`` がすべて None。``q`` は結果の
    値をそのまま返し、ここで BH を掛け直さない（母集団は結果が持つ）。``tested`` でも q が
    None のことはある。``detected_a/b`` は検出例数で、マスクから対応を保証して取れないとき
    None（呼び出し側が注記する）。取れるときは status に関わらず付ける。

    ``with_labels=True`` のときだけ ``comparison``・``group_a``・``group_b``・
    ``log2fc_test_over_ref``（v1 の ``log2fc`` と同じ値。group_b / group_a）を足す。
    """
    prepared = []
    for result in results:
        by_feature = {row.get("feature"): row for row in result.get("results") or []}
        prepared.append((_result_id(result), by_feature, _detected_counter(ds, result), result))
    rows = []
    for fid in feature_ids:
        for result_id, by_feature, counter, result in prepared:
            source = by_feature.get(fid)
            if source is None:
                status = "not_in_matrix"
            elif _is_number(source.get("p")):
                status = "tested"
            else:
                status = "not_testable"
            tested = status == "tested"
            det_a, det_b = counter(fid) if counter is not None else (None, None)
            rows.append({
                "feature_id": fid,
                "result_id": result_id,
                "status": status,
                "log2fc": _num_or_none(source.get("log2fc")) if tested else None,
                "p": _num_or_none(source.get("p")) if tested else None,
                "q": _num_or_none(source.get("q")) if tested else None,
                "mean_a": _num_or_none(source.get("mean_a")) if tested else None,
                "mean_b": _num_or_none(source.get("mean_b")) if tested else None,
                "detected_a": det_a,
                "detected_b": det_b,
            })
            if with_labels:
                rows[-1].update({"comparison": comparison_name(result),
                                 "group_a": result.get("a"), "group_b": result.get("b"),
                                 "log2fc_test_over_ref": rows[-1]["log2fc"]})
    return rows


def result_summaries(ds, results: list[dict]) -> list[dict]:
    """結果ごとの要約（群・n・log2FC の向き・BH の母集団・注意書き）。

    ``bh_population``: ``n_features_total`` はデータセットの特徴数、``n_in_matrix`` は結果の
    行数（前処理後に残った特徴）、``n_tested`` は p が有限の行数。q の解釈にはこの母集団が要る
    （窓で絞った特徴の数ではなく、この母集団で BH が掛かっている）。
    """
    summaries = []
    for result in results:
        rows = result.get("results") or []
        summaries.append({
            "result_id": _result_id(result),
            "group_a": result.get("a"),
            "group_b": result.get("b"),
            "n_a": result.get("n_a"),
            "n_b": result.get("n_b"),
            "log2fc_sign": result.get("log2fc_sign"),
            "bh_population": {
                "n_features_total": len(ds.feature_ids),
                "n_in_matrix": len(rows),
                "n_tested": sum(1 for row in rows if _is_number(row.get("p"))),
            },
            "caveats": list(result.get("caveats") or []),
        })
    return summaries


# ---------- 割れたスポットのまとめと検定のやり直し（spec 2026-10-28 §4.3 (c)） ----------

def merge_plan(ds, *, alignment_qc_path: str | None = None, dataset_hash: str | None = None,
               field: str = "alignment_qc_path") -> dict:
    """まとめる割れの群を決める。``{"source", "groups", "skipped", "alignment_qc_sha256"}``。

    群は `alignment_qc` の `merge_groups`（互いに最も近い mergeable な 2 スポット）だけで、
    ここでは足さない・広げない（`analysis_spec` の `features.merge_split` と同じ規則・同じ実装）。
    ``alignment_qc_path`` は `dataset_alignment_qc` が書いた JSON（今のデータセットの指紋で
    作ったものだけを受ける）で、無ければその場で `find_split_pairs(ds)` を既定の閾値で回す。
    ``groups`` は残す ID が先頭の並び、``skipped`` はデータセットに 2 つ揃わない群。
    ``dataset_hash`` は計算済みの `dataset_fingerprint(ds)`（任意）、``field`` はエラーに書く欄の名前。
    """
    from metabolomix.analysis import alignment_qc

    if alignment_qc_path is None:
        qc = alignment_qc.find_split_pairs(ds, dataset_hash=dataset_hash)
        source, sha256 = "find_split_pairs", None
    else:
        qc, sha256 = alignment_qc.read_qc_file(alignment_qc_path, ds, field=field,
                                               dataset_hash=dataset_hash)
        source = "alignment_qc"
    usable, skipped = alignment_qc.usable_merge_groups(qc, ds.feature_ids)
    return {"source": source, "groups": usable, "skipped": skipped,
            "alignment_qc_sha256": sha256}


def merged_dataset(ds, groups) -> tuple[object, dict[str, list[str]]]:
    """割れの群を試料ごとの大きい方の値で 1 列にまとめた ``ds`` の写しと ``merged_from``。

    元の ``ds`` は変えない（浅い写しの定量行列・特徴 ID・検出状態だけを置き換える）。まとめは
    前処理の**前**の生の値に当て、検出状態は 2 スポットの「どちらかで検出」（`analysis_spec` と
    同じ順序。前処理の篩と正規化が 1 つの物質を 1 列として見る）。実装は
    `alignment_qc.merge_matrix_columns`（残す列は群の先頭の特徴 ID）。
    """
    from metabolomix.analysis import alignment_qc

    values = np.asarray(ds.feature_matrix, dtype=float)
    mask = getattr(ds, "detected_mask", None)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != values.shape:
            mask = None
    merged = alignment_qc.merge_matrix_columns(values.T, ds.feature_ids, groups,
                                               detected=None if mask is None else mask.T)
    view = copy.copy(ds)
    view.feature_matrix = merged.matrix.T
    view.feature_ids = list(merged.feature_ids)
    view.detected_mask = None if merged.detected is None else merged.detected.T
    return view, merged.merged_from


def _result_recipe(ds, result: dict) -> dict | None:
    """その 2 群比較を計算した前処理のレシピ（分からなければ None）。

    pipeline の run から引き継いだ結果は来歴の ``handoff.preprocess``
    （`session_handoff._run_preprocess_parameters` が写した run の ``resolved_recipe``）、
    session で計算した結果は親（前処理の結果）の ``effective_parameters``。
    """
    provenance = result.get("provenance") or {}
    handoff = provenance.get("handoff")
    if isinstance(handoff, dict):
        recipe = handoff.get("preprocess")
        return dict(recipe) if isinstance(recipe, dict) else None
    for parent in provenance.get("parent_ids") or []:
        params = (((getattr(ds, "results", None) or {}).get(parent) or {}).get("provenance")
                  or {}).get("effective_parameters")
        if isinstance(params, dict):
            return dict(params)
    return None


def _retest_merged_engine(ds, result: dict, groups) -> dict:
    """統計エンジンの結果（`select_results` が展開したもの）を、割れをまとめた列で検定し直す。

    記録の spec（入力の除外・群の上書き・前処理・設計・対比・検定の設定）をそのまま使い、
    `from_dataset` で作った StudyData の割れた特徴の列を `merge_feature_columns` でまとめてから
    `run_spec_v2` に渡す（StudyData の source には spec の `input` が再適用されないので、ここで当てる）。
    BH の母集団はまとめた後の全特徴。個体・注入の対応は記録の `units` のまま（まとめても変わらない）。
    """
    from metabolomix.analysis.stats.data import from_dataset, merge_feature_columns
    from metabolomix.analysis.stats.spec_v2 import run_spec_v2

    source = result["_source"]
    record = ds.results.get(source["result_id"])
    result_id = _result_id(result)
    if record is None or record.get("kind") != "compare":
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"差次的結果 {result_id} の元の記録が無いので、割れたスポットをまとめて"
                          "検定し直せない。", {"result_id": result_id})
    spec = record["spec"]
    inputs = spec["input"]
    try:
        study = from_dataset(ds, exclude_samples=inputs["exclude_samples"],
                             exclude_features=inputs["exclude_features"],
                             group_overrides=inputs["group_overrides"])
        run = run_spec_v2(spec, merge_feature_columns(study, [list(g) for g in groups]))
        analysis = next((a for a in run.get("analyses") or []
                         if a["analysis_id"] == record["analysis_id"]), None)
    except DomainError as exc:
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"まとめた行列で {result_id} の比較をやり直せない: {exc.message}",
                          {"result_id": result_id, "cause": exc.code}) from exc
    if analysis is None:
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"まとめた行列に {result_id} の前処理を当て直せない"
                          "（利用者への確認が要る前処理）。", {"result_id": result_id})
    redone = as_two_group({**record, "engine": analysis["result"]}, source["contrast_id"])
    return {**result, "results": redone["results"], "volcano": redone["volcano"],
            "statistics": redone["statistics"]}


def retest_merged(ds, view, result: dict, *, groups=()) -> dict:
    """まとめた行列（``view``）で、``result`` と同じ前処理・同じ 2 群検定をやり直す。

    統計エンジンの結果（`select_results` が `compare` の記録から展開したもの）は、``groups`` を
    `merge_feature_columns` で当てて記録の spec のまま `run_spec_v2`（`_retest_merged_engine`）。
    以下は v1 の結果の規則:

    前処理は ``result`` を計算したレシピ（`_result_recipe`）を `run_dataset_preprocess` で
    当て直し、検定は `differential.two_group_test`（``result`` の ``log_transform``。pipeline の
    既定は log 変換の Welch）と `add_fdr`（**まとめた後の全特徴**で BH）。群は ``result`` の
    ``samples_a`` / ``samples_b``（`dataset_service.compare_dataset` が検定に使った試料）。
    検定と補正の種類は変えないので、判定の規則の意味は同じ。戻り値は ``result`` の写しで、
    ``results`` だけがやり直した行。レシピが分からない・前処理を当て直せない・群の試料が
    前処理の後に無いときは ``FEATURE_MERGE_UNAVAILABLE``（まとめる前の数字をまとめた後の顔で
    返さない）。
    """
    from metabolomix.analysis import differential
    from metabolomix.analysis.dataset_analysis import PreconditionError, run_dataset_preprocess

    if result.get("_source"):
        return _retest_merged_engine(ds, result, groups)
    result_id = _result_id(result)
    recipe = _result_recipe(ds, result)
    if recipe is None:
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"差次的結果 {result_id} を計算した前処理のレシピが分からないので、割れた"
                          "スポットをまとめて検定し直せない。merge_split=false で呼ぶか、"
                          "dataset_compare で計算し直す。", {"result_id": result_id})
    try:
        matrix, names, features, *_ = run_dataset_preprocess(view, recipe)
    except PreconditionError as exc:
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"まとめた行列に {result_id} の前処理を当て直せない: {exc.message}",
                          {"result_id": result_id, "recipe": recipe}) from exc
    position = {name: i for i, name in enumerate(names)}
    samples_a = list(result.get("samples_a") or [])
    samples_b = list(result.get("samples_b") or [])
    absent = [name for name in samples_a + samples_b if name not in position]
    if absent or not samples_a or not samples_b:
        raise DomainError("FEATURE_MERGE_UNAVAILABLE",
                          f"差次的結果 {result_id} の群の試料が、前処理を当て直した行列に無い"
                          f"（{absent or 'samples_a / samples_b が空'}）。",
                          {"result_id": result_id, "absent_samples": absent})
    group_a, group_b = result.get("a"), result.get("b")
    labels = [group_a] * len(samples_a) + [group_b] * len(samples_b)
    rows = differential.two_group_test(
        matrix[[position[name] for name in samples_a + samples_b], :], features, labels,
        group_a, group_b, log_transform=bool(result.get("log_transform", True)))
    return {**result, "results": differential.add_fdr(rows)}


def merge_map(groups) -> dict[str, list[str]]:
    """割れの群 → ``{残す ID: 群}``（残す ID は群の先頭。`merge_matrix_columns` と同じ）。"""
    return {group[0]: list(group) for group in groups if len(group) >= 2}


def merged_results(ds, results: list[dict], groups, feature_ids) -> dict:
    """まとめる群が 1 つでもあれば、まとめた行列で ``results`` を全特徴について検定し直す。

    戻り値 ``{"view", "merged_from", "touched", "results"}``。``merged_from`` は全群（残す ID →
    群）、``touched`` はそのうち ``feature_ids`` のどれか（消える側の ID を含む）を含む群（報告用）、
    ``results`` は result_id → やり直した結果。まとめは全群に当て、BH の母集団は「まとめた後の
    全特徴」——1 つの表の q はすべてこの母集団から出す（まとめた行と従来の行を混ぜない。
    コントローラの決定 2026-10-29）。群が無ければ検定はやり直さず、``view`` は ``ds``。
    """
    if not groups:
        return {"view": ds, "merged_from": {}, "touched": {}, "results": {}}
    view, merged_from = merged_dataset(ds, groups)
    wanted = {str(f) for f in feature_ids}
    touched = {kept: group for kept, group in merged_from.items() if wanted.intersection(group)}
    redone = {_result_id(result): retest_merged(ds, view, result, groups=groups) for result in results}
    return {"view": view, "merged_from": merged_from, "touched": touched, "results": redone}


def merged_feature_stats(ds, feature_ids: list[str], results: list[dict], *, groups,
                         with_labels: bool = False) -> dict:
    """`feature_stats` の割れのまとめ版。

    戻り値 ``{"stats", "stats_unmerged", "touched", "renamed", "results", "view"}``。
    ``stats_unmerged`` は従来の ``feature_stats(ds, feature_ids, results)`` そのもの。まとめる群が
    1 つでもあれば、``stats`` は**すべての行**をまとめた行列で検定し直した値（`merged_results`。
    q はまとめた後の全特徴の BH）から出し、まとめて消える特徴 ID は残す ID（群の先頭）に
    置き換える（重複は 1 つ）。群が無ければ ``stats`` は従来の値と同じ。各行の
    ``merged_from`` は群（残す ID が先頭）か None。``renamed`` は置き換えた ID（消える ID → 残す ID）。
    """
    unmerged = feature_stats(ds, feature_ids, results, with_labels=with_labels)
    merged = merged_results(ds, results, groups, feature_ids)
    merged_from = merged["merged_from"]
    rename = {fid: kept for kept, group in merged_from.items() for fid in group if fid != kept}
    folded = list(dict.fromkeys(rename.get(fid, fid) for fid in feature_ids))
    if merged["results"]:
        stats = feature_stats(merged["view"], folded,
                              [merged["results"][_result_id(r)] for r in results],
                              with_labels=with_labels)
    else:
        stats = feature_stats(ds, folded, results, with_labels=with_labels)
    for row in stats:
        row["merged_from"] = merged_from.get(row["feature_id"])
    return {"stats": stats, "stats_unmerged": unmerged, "touched": merged["touched"],
            "renamed": {fid: rename[fid] for fid in feature_ids if fid in rename},
            "results": merged["results"], "view": merged["view"]}
