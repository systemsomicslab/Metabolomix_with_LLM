"""dataset_feature_search の純ロジック: 目標の正規化・理論 m/z・規則 hash（spec 2026-10-19 §4.1）、
目標照合の判定（§5: `classify_identity` / `match_target`）、窓検索と差次的統計の結合
（§4.2・§6: `search_window` / `select_results` / `feature_stats` / `result_summaries`）。

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

from metabolomix.analysis.feature_bindings import _normalize_adduct
from metabolomix.analysis.result_state import is_current
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.msdial.adducts import mz_from_neutral, parse_adduct
from metabolomix.msdial.peak_verification import monoisotopic_mass, parse_formula

__all__ = [
    "CANDIDATE_COLUMNS",
    "DEFAULT_ADDUCT_TIERS",
    "EVIDENCE_LEVELS",
    "MAX_TARGETS",
    "STATS_COLUMNS",
    "WINDOW_COLUMNS",
    "classify_identity",
    "feature_stats",
    "match_target",
    "normalize_targets",
    "result_summaries",
    "rule_hash",
    "search_window",
    "select_results",
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
_POLARITY_OF_SIGN = {sign: polarity for polarity, sign in _SIGN_OF.items()}
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
    目標の極性は ``polarity``、無ければ明示の付加体（``adduct`` か ``adduct_tiers``）の
    符号（混ざっていれば ``FEATURE_TARGET_INVALID``）、それも無ければ ``dataset_polarity``。
    目標の極性が ``dataset_polarity`` と違う目標は ``applicable=False``（→ ``not_applicable``）。
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


#: rule_hash に入れない欄。``name`` は読み手のためのラベルで照合に使わない
#: （改名しただけで「規則が変わった」ように見せない）。
_RULE_HASH_EXCLUDED = frozenset({"name"})


def rule_hash(normalized_targets: list[dict]) -> str:
    """正規化後の目標リストから照合規則の hash を作る（result_id・dataset_hash は入れない）。

    覆うのは照合を左右する欄だけ: ``target_id``・InChIKey 先頭ブロック・分子式か
    expected_mz・付加体の層と理論 m/z・ppm・RT・目標の極性・``applicable``。``name`` は
    入れず、目標は ``target_id`` の順に並べ替えてから作る——表示名の変更や
    ``targets`` の並べ替えでは hash は変わらない。
    """
    canonical = sorted(
        ({k: v for k, v in target.items() if k not in _RULE_HASH_EXCLUDED}
         for target in normalized_targets),
        key=lambda target: target["target_id"])
    return canonical_hash(canonical)


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
      価数のすべてで一致 → ``sme_match``。電荷が None の候補は、付加体が一致していれば
      電荷の確認を通す（欠測は反証ではない）。
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


def match_target(ds, target: dict, *, max_results: int,
                 mz_array: np.ndarray | None = None) -> dict:
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
    if mz_array is None:
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
WINDOW_COLUMNS = ("feature_id", "mz", "rt", "evidence", "name", "inchikey", "adduct", "n_sme")

#: 特徴ごとの統計の行の列。
STATS_COLUMNS = ("feature_id", "result_id", "status", "log2fc", "p", "q", "mean_a", "mean_b",
                 "detected_a", "detected_b")


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


def search_window(ds, *, feature_ids: list[str] | None, min_mz, max_mz, min_rt, max_rt,
                  max_results: int) -> tuple[list[dict], int]:
    """m/z・RT（分）の窓と特徴 ID の指定で特徴を探す（spec §4.2）。

    条件は与えられたものすべての AND で、境界を含む。``feature_ids`` は None のときだけ
    「条件なし」（空リストは「どの ID でもない」＝ 0 件）。条件が 1 つも無ければ
    ``FEATURE_SEARCH_EMPTY_QUERY``（全特徴の列挙を黙って返さない）、下限 > 上限は
    ``FEATURE_SEARCH_BAD_RANGE``。m/z や RT を持たない特徴は、その軸の条件があれば外れる。
    行は m/z 昇順（m/z 無しは最後。同値は特徴の並び順）で、``max_results`` 件に打ち切り、
    切り捨てた件数を返す。統計はここでは付けない（`feature_stats`）。
    """
    if all(v is None for v in (feature_ids, min_mz, max_mz, min_rt, max_rt)):
        raise DomainError(
            "FEATURE_SEARCH_EMPTY_QUERY",
            "feature_ids・min_mz/max_mz・min_rt/max_rt のどれか 1 つは指定する",
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
    for fid in ds.feature_ids:
        if wanted is not None and fid not in wanted:
            continue
        meta = ds.feature_metadata.get(fid) or {}
        mz, rt = meta.get("mz"), meta.get("rt")
        if _inside(mz, min_mz, max_mz) and _inside(rt, min_rt, max_rt):
            hits.append((fid, mz, rt))
    hits.sort(key=lambda hit: (not _is_number(hit[1]), hit[1] if _is_number(hit[1]) else 0.0))

    limit = max(0, int(max_results))
    rows = []
    for fid, mz, rt in hits[:limit]:
        identity = _window_identity(ds.feature_candidates.get(fid), ds.feature_annotations.get(fid))
        rows.append({
            "feature_id": fid,
            "mz": float(mz) if _is_number(mz) else None,
            "rt": float(rt) if _is_number(rt) else None,
            **identity,
        })
    return rows, len(hits) - len(rows)


def _result_id(result: dict):
    return (result.get("provenance") or {}).get("result_id")


def select_results(ds, result_ids: list[str] | None, *,
                   dataset_hash: str | None = None) -> dict:
    """統計を引く 2 群比較の結果を選ぶ。``{"results", "stale", "missing"}``。

    省略（None）なら ``ds.results`` のうち ``kind == "two_group"`` のすべて。名指しなら、その
    ID を引く。どちらも ``is_current`` でない結果は ``results`` に入れず ``stale`` に ID を
    挙げる（古い統計を現在のデータの顔で返さない）。名指しの ID が無い、または 2 群比較で
    ない結果（v2 の ``dataset_statistic`` = ``kind == "welch"`` など）は ``missing``。

    ``dataset_hash`` は計算済みの ``dataset_fingerprint(ds)``（任意）。引き継いだ結果ごとに
    行列全体の hash を計算し直さないよう ``is_current`` へそのまま渡す。
    """
    selected, stale, missing = [], [], []
    if result_ids is None:
        pairs = [(rid, r) for rid, r in ds.results.items() if (r or {}).get("kind") == "two_group"]
    else:
        pairs = []
        for rid in result_ids:
            result = ds.results.get(rid)
            if (result or {}).get("kind") == "two_group":
                pairs.append((rid, result))
            else:
                missing.append(rid)
    for rid, result in pairs:
        if is_current(ds, result, dataset_hash=dataset_hash):
            selected.append(result)
        else:
            stale.append(rid)
    return {"results": selected, "stale": stale, "missing": missing}


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


def feature_stats(ds, feature_ids: list[str], results: list[dict]) -> list[dict]:
    """特徴 × 結果の統計行（特徴が外側、結果が内側）。注釈の有無に関わらず全特徴に付ける。

    ``status``: ``tested``（p が有限）/ ``not_testable``（結果の行はあるが p が None・NaN。
    pipeline run から復元した結果の JSON null を含む）/ ``not_in_matrix``（前処理で落ち、
    結果の行列に無い）。後の 2 つは ``log2fc p q mean_a mean_b`` がすべて None。``q`` は結果の
    値をそのまま返し、ここで BH を掛け直さない（母集団は結果が持つ）。``tested`` でも q が
    None のことはある。``detected_a/b`` は検出例数で、マスクから対応を保証して取れないとき
    None（呼び出し側が注記する）。取れるときは status に関わらず付ける。
    """
    prepared = []
    for result in results:
        by_feature = {row.get("feature"): row for row in result.get("results") or []}
        prepared.append((_result_id(result), by_feature, _detected_counter(ds, result)))
    rows = []
    for fid in feature_ids:
        for result_id, by_feature, counter in prepared:
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
