"""MS/MS の合否（msms-rules.v1）・診断イオン・中性損失（spec 2026-10-22 §4 の 6）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。照合の結果
（`library.matching.score_candidates` の `scored`）と測定スペクトルを受け、spec に書いた
規則を機械的に当てる。**規則に無い判断はしない**（「明らかに別物」かどうかは再解析の
計画の規則で決める。道具は材料と、書かれた規則への合否だけを返す）。

設定していない検査は行を出さない。`verdict` は、すべて pass → `qualified`、
どれか fail → `not_qualified`、fail が無く判定できない検査がある → `not_evaluable`。

既定の規則でも断片の一致を 1 本は求める（`min_fragments_excl_precursor` の既定 1。一致 0 本で
qualified にしない）。断片の一致は、測定の基準ピークに対する強度が `min_fragment_relative_intensity`
（既定 0.01）以上の窓だけを数える（雑音の 1 本で qualified にしない。2026-11-13 の事実確認 A1）。
"""
from __future__ import annotations

import copy
import math
import re
from functools import lru_cache

from metabolomix.analysis.spectral_match import (
    DEFAULT_MIN_FRAGMENT_RELATIVE_INTENSITY, fragment_matches_excluding_precursor)
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = ["CANDIDATE_COLUMNS", "RULES_SCHEMA", "candidate_rows", "find_ions", "isotope_labelled",
           "name_lipid_species", "normalize_ion", "normalize_ions", "normalize_rules", "qualify", "query_window",
           "record_lipid_species", "rules_hash"]

RULES_SCHEMA = "msms-rules.v1"
_INVALID = "MSMS_RULES_INVALID"
_BLOCK_RE = re.compile(r"^[A-Z]{14}$")
DEFAULT_ION_TOL = 0.01

_DEFAULTS = {
    "top_rank": 5,
    "min_fragments_excl_precursor": 1,
    "min_fragment_relative_intensity": DEFAULT_MIN_FRAGMENT_RELATIVE_INTENSITY,
    "precursor_exclusion_tol": None,
    "same_formula_margin": None,
    "lipid_class": None,
    "diagnostic_ions": [],
    "neutral_losses": [],
    "ions_match": "all",
}


def _fail(message: str, **details) -> DomainError:
    return DomainError(_INVALID, message, details)


def _positive(value, label: str, *, allow_zero=False) -> float:
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            or value < 0 or (value == 0 and not allow_zero)):
        raise _fail(f"{label} は正の数で指定する: {value!r}", field=label)
    return float(value)


def inchikey_block(value) -> str:
    block = value.strip()[:14].upper() if isinstance(value, str) else ""
    if not _BLOCK_RE.match(block):
        raise _fail("target.inchikey は InChIKey（先頭 14 文字の英大文字ブロックを使う）",
                    field="target.inchikey")
    return block


def lipid_species(value) -> str:
    """目標の脂質の種（`PC 34:4`・`PC 18:3_16:1` など）を種のレベルに直す。

    同位体の標識のある名前（内部標準 `PC 15:0_18:1(d7)`）は断る: ライブラリの記録も標識付きは種に読まない
    （`name_lipid_species`）ので種として合否を出せず、GOSLIN は標識を落として読む（PC 33:1）ので、受けると
    天然の種の合否を内部標準の合否として返してしまう。"""
    from metabolomix.msdial.lipid_identity import lipid_species_of

    if isotope_labelled(value):
        raise _fail(f"target.lipid_species に同位体の標識（内部標準）がある: {value!r}。ライブラリの記録も標識付きは"
                    "種に読まないので種として合否を出せず、標識を外して読むと天然の種の合否になる"
                    "（天然の種を確かめるなら標識の無い名前を渡す）", field="target.lipid_species")
    species = lipid_species_of(value) if isinstance(value, str) else None
    if species is None:
        raise _fail("target.lipid_species は種のレベル以下の脂質名（例 PC 34:4。クラスだけは不可）",
                    field="target.lipid_species")
    return species


def _normalize_target(raw_target: dict, formula) -> dict:
    """InChIKey の目標は `{inchikey_block, formula}`、種の目標は `{lipid_species, formula}`。"""
    formula = formula.strip() if formula else None
    if raw_target.get("lipid_species") is None:
        return {"inchikey_block": inchikey_block(raw_target.get("inchikey")), "formula": formula}
    if raw_target.get("inchikey") is not None:
        raise _fail("target.inchikey と target.lipid_species は同時に指定できない",
                    field="target.lipid_species")
    return {"lipid_species": lipid_species(raw_target["lipid_species"]), "formula": formula}


def normalize_ion(item, *, kind: str, label: str, loss_key: str = "delta") -> dict:
    """診断イオン（`mz`）か中性損失（`delta`）の 1 件。誤りの文は原因のキーを名指す。

    `min_relative_intensity`（測定の最大ピークに対する比の下限）は両方で受ける。中性損失の比は
    生成イオン（前駆体 − `delta`）の測定強度の比。診断イオンは省くと 0 を埋め、中性損失は
    渡したときだけ載せる（省いた規則の `rule_hash` を変えない）。`loss_key` は呼び手が中性損失の
    値に使うキーの名前（`dcl_diagnostic_ions` は `neutral_loss`）。戻り値はいつも `delta`。
    """
    key = "mz" if kind == "ion" else loss_key
    allowed = {key, "label", "tol", "min_relative_intensity"}
    if not isinstance(item, dict):
        raise _fail(f"{label} は {{{key}, label, tol, min_relative_intensity}} の辞書", field=label)
    unknown = sorted(set(item) - allowed)
    if unknown:
        raise _fail(f"{label} の未知のキー: {unknown}（使えるのは {key}・label・tol・"
                    "min_relative_intensity）", field=label, unknown=unknown)
    if key not in item:
        raise _fail(f"{label} に {key} が無い", field=f"{label}.{key}")
    out = {
        ("mz" if kind == "ion" else "delta"): _positive(item[key], f"{label}.{key}"),
        "label": str(item.get("label") or f"{key}={item[key]}"),
        "tol": _positive(item.get("tol", DEFAULT_ION_TOL), f"{label}.tol"),
    }
    if kind == "ion" or "min_relative_intensity" in item:
        out["min_relative_intensity"] = _positive(item.get("min_relative_intensity", 0.0),
                                                  f"{label}.min_relative_intensity",
                                                  allow_zero=True)
    return out


def normalize_ions(items, *, kind: str, label: str) -> list[dict]:
    """診断イオン（`mz`）か中性損失（`delta`）の一覧。各要素は `normalize_ion`。"""
    if not isinstance(items, list):
        raise _fail(f"{label} はリスト", field=label)
    return [normalize_ion(item, kind=kind, label=f"{label}[{i}]") for i, item in enumerate(items)]


def normalize_rules(rules, *, target: dict | None = None) -> dict:
    """規則を検証して既定を埋める。`target` を渡すと規則の target はそれで決まる（証拠の束）。"""
    if not isinstance(rules, dict):
        raise _fail("rules は辞書で指定する")
    allowed = {"schema", "target", *_DEFAULTS}
    unknown = sorted(set(rules) - allowed)
    if unknown:
        raise _fail(f"rules の未知のキー: {unknown}", unknown=unknown)
    if rules.get("schema", RULES_SCHEMA) != RULES_SCHEMA:
        raise _fail(f"schema は {RULES_SCHEMA}")
    raw_target = target if target is not None else rules.get("target")
    if not isinstance(raw_target, dict) or set(raw_target) - {"inchikey", "lipid_species", "formula"}:
        raise _fail("target は {inchikey | lipid_species, formula} の辞書", field="target")
    formula = raw_target.get("formula")
    if formula is not None and (not isinstance(formula, str) or not formula.strip()):
        raise _fail("target.formula は文字列", field="target.formula")
    out = copy.deepcopy(_DEFAULTS)
    out.update({k: copy.deepcopy(v) for k, v in rules.items() if k in _DEFAULTS})
    top_rank = out["top_rank"]
    if isinstance(top_rank, bool) or not isinstance(top_rank, int) or top_rank < 1:
        raise _fail("top_rank は 1 以上の整数", field="top_rank")
    if out["min_fragments_excl_precursor"] is not None:
        value = out["min_fragments_excl_precursor"]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise _fail("min_fragments_excl_precursor は 0 以上の整数か null",
                        field="min_fragments_excl_precursor")
    floor = out["min_fragment_relative_intensity"]
    if (isinstance(floor, bool) or not isinstance(floor, (int, float)) or not math.isfinite(floor)
            or not 0 <= floor <= 1):
        raise _fail("min_fragment_relative_intensity は 0 以上 1 以下の数（測定の基準ピークに対する比）",
                    field="min_fragment_relative_intensity")
    out["min_fragment_relative_intensity"] = float(floor)
    if out["precursor_exclusion_tol"] is not None:
        out["precursor_exclusion_tol"] = _positive(out["precursor_exclusion_tol"],
                                                   "precursor_exclusion_tol")
    margin = out["same_formula_margin"]
    if margin is not None:
        if not isinstance(margin, dict) or set(margin) - {"score", "min_delta"}:
            raise _fail("same_formula_margin は {score, min_delta}", field="same_formula_margin")
        score = margin.get("score", "total_score")
        if score not in ("total_score", "weighted_dot_product", "reverse_dot_product",
                         "simple_dot_product", "entropy_similarity"):
            raise _fail(f"same_formula_margin.score を読めない: {score!r}",
                        field="same_formula_margin.score")
        out["same_formula_margin"] = {"score": score,
                                      "min_delta": float(margin.get("min_delta", 0.0))}
    if out["lipid_class"] is not None and (not isinstance(out["lipid_class"], str)
                                           or not out["lipid_class"].strip()):
        raise _fail("lipid_class は文字列か null", field="lipid_class")
    if out["ions_match"] not in ("all", "any"):
        raise _fail("ions_match は all / any", field="ions_match")
    out["diagnostic_ions"] = normalize_ions(out["diagnostic_ions"], kind="ion",
                                            label="diagnostic_ions")
    out["neutral_losses"] = normalize_ions(out["neutral_losses"], kind="loss",
                                           label="neutral_losses")
    return {"schema": RULES_SCHEMA, "target": _normalize_target(raw_target, formula), **out}


def rules_hash(normalized: dict) -> str:
    return canonical_hash(normalized)


# ---------- 照会 ----------

def find_ions(measured, ions: list[dict], *, precursor_mz: float | None) -> list[dict]:
    """測定スペクトルで診断イオン（`mz`）・中性損失（`delta`）を探す。

    相対強度は測定スペクトルの最大ピークに対する比。±tol に複数あれば最も強いもの。
    """
    peaks = [(float(mz), float(intensity)) for mz, intensity in (measured or [])]
    top = max((i for _, i in peaks), default=0.0)
    out = []
    for ion in ions:
        if "mz" in ion:
            kind, target = "ion", ion["mz"]
        else:
            kind = "neutral_loss"
            target = None if precursor_mz is None else precursor_mz - ion["delta"]
        hits = [] if target is None else [(mz, i) for mz, i in peaks
                                          if abs(mz - target) <= ion["tol"] + 1e-9]
        best = max(hits, key=lambda p: p[1], default=None)
        relative = best[1] / top if best and top > 0 else None
        minimum = ion.get("min_relative_intensity", 0.0)
        out.append({
            "label": ion["label"], "kind": kind, "target_mz": target,
            "found": bool(best) and (relative or 0.0) >= minimum,
            "observed_mz": best[0] if best else None,
            "relative_intensity": relative,
        })
    return out


def _block(record) -> str | None:
    key = record.get("inchikey")
    return key.strip()[:14].upper() if isinstance(key, str) and key.strip() else None


#: 記録の種を得た経路（`record_lipid_species` の 2 つめの値）。
SPECIES_ROUTES = ("name", "name_cleaned", "ontology_formula")


def _labelled_rest(name) -> tuple[bool, str | None]:
    """限定子を外した名前と、それが同位体の標識（内部標準）を持つか。名前が無ければ (False, None)。"""
    from metabolomix.pathway import msp_lipid

    if not isinstance(name, str) or not name.strip():
        return False, None
    _qualifier, rest = msp_lipid.split_qualifier(name)
    return msp_lipid.has_isotope_label(rest), rest


def isotope_labelled(name) -> bool:
    """名前（MS-DIAL の限定子を外して）に同位体の標識（内部標準。`(d7)`・`-d31`・`13C`）があるか。

    `name_lipid_species` が種に読まない名前。種の目標（`lipid_species`）の側もこれで断る。"""
    return _labelled_rest(name)[0]


@lru_cache(maxsize=65536)
def name_lipid_species(name) -> tuple[str | None, str | None]:
    """名前だけから脂質の種（GOSLIN の種のレベル）と、それを得た経路（`name` / `name_cleaned`）。

    `record_lipid_species` の名前の 2 つの経路そのもので、ライブラリの記録・SME 候補の `chemical_name`・SML と
    MS-DIAL の注釈の名前をすべてこれで読む（目標照合の種の目標・MS/MS の合否・`msdial_agreement` が同じ名前を
    同じ種に読む）。`name`（限定子を外してそのまま）→ `name_cleaned`（注記・付加体の表記を外し古い書式を略記に
    直す。`pathway.msp_lipid.shorthand`。`PE 36:2;O (oxPE)` → `PE 36:2;O`）。同位体の標識のある名前（内部標準）と
    読めない名前は (None, None)。
    """
    from metabolomix.msdial.lipid_identity import lipid_species_of
    from metabolomix.pathway import msp_lipid

    labelled, rest = _labelled_rest(name)
    if rest is None or labelled:
        return None, None          # 同位体の標識のある名前（内部標準）は種に読まない
    species = lipid_species_of(name)
    if species is not None:
        return species, "name"
    shorthand = msp_lipid.shorthand(rest)
    species = lipid_species_of(shorthand) if shorthand else None
    if species is not None:
        return species, "name_cleaned"
    return None, None


@lru_cache(maxsize=65536)
def _species_of(name, ontology, formula) -> tuple[str | None, str | None]:
    from metabolomix.msdial.lipid_identity import lipid_species_of
    from metabolomix.pathway import msp_lipid

    labelled, rest = _labelled_rest(name)
    if rest is None or labelled:
        return None, None          # 同位体の標識のある記録（内部標準）は種に読まない
    species, route = name_lipid_species(name)
    if species is not None:
        return species, route
    built = msp_lipid.from_class_formula(ontology, formula, rest)
    if isinstance(built, tuple):
        species = lipid_species_of(built[0])
        if species is not None:
            return species, "ontology_formula"
    return None, None


def record_lipid_species(record) -> tuple[str | None, str | None]:
    """ライブラリの記録の脂質の種（GOSLIN の種のレベル）と、それを得た経路。読めなければ (None, None)。

    経路は順に `name`（記録の name をそのまま GOSLIN で読む）→ `name_cleaned`（ここまでは `name_lipid_species`。限定子・`|` の後・
    `;` の後の注記（`;` の直後に空白の無い修飾 `;O…` `;2OH` `;(2OH)` `;Ep` … は残す。GOSLIN が読めない修飾なら
    種に読まない）・付加体の表記を外し、古い MS-DIAL の書式と英語のクラス名を略記に直して読む。
    `pathway.msp_lipid.shorthand`）→ `ontology_formula`（ClassyFire の Ontology と分子式から和の組成の
    種を作る。`pathway.msp_lipid.from_class_formula`）。同位体の標識のある名前（内部標準）は読まない。
    """
    if not isinstance(record, dict):
        return None, None
    formula = record.get("formula")
    return _species_of(record.get("name"), record.get("ontology"),
                       formula.strip() if isinstance(formula, str) else None)


def _target_key(rules: dict):
    """記録 → 目標と比べる識別子（InChIKey の先頭ブロックか、記録の種）と、目標の値。"""
    target = rules["target"]
    if "lipid_species" in target:
        return (lambda record: record_lipid_species(record)[0]), target["lipid_species"]
    return _block, target["inchikey_block"]


def _lipid_class(name) -> str | None:
    from metabolomix.msdial.lipid_identity import normalize_lipid_name

    if not isinstance(name, str) or not name.strip():
        return None
    parsed = normalize_lipid_name(name)
    if not parsed.get("parse_ok") or not parsed.get("normalized"):
        return None
    return parsed["normalized"].split(" ", 1)[0]


def _check(check: str, status: str, value=None, detail=None) -> dict:
    return {"check": check, "status": status, "value": value, "detail": detail}


def qualify(scored: list, measured, precursor_mz: float, rules: dict, *, ms2_tol: float) -> dict:
    """照合の結果（`score_candidates` の `scored`。総合スコアの降順）に規則を当てる。

    種の目標（`target.lipid_species`）では、目標か別物かを記録の種（`record_lipid_species`。
    名前 → 注記を外した名前 → Ontology と分子式の順に読む）で決める（同じ種の別の分子種は
    目標として数え、`same_formula_margin` の相手にしない）。目標の候補の種を得た経路は
    `target_in_top_rank` の `detail.species_route`。

    戻り値の `top1_is_target` は 1 位の候補が目標か（目標が候補に無ければ False）。既定の
    `top_rank`（5）では 1 位が別物でも `target_in_top_rank` は pass しうるので、`target_rank` と
    並べて見る。
    """
    key_of, block = _target_key(rules)
    species_target = "lipid_species" in rules["target"]
    checks = []
    target_index = next((i for i, (record, _) in enumerate(scored) if key_of(record) == block),
                        None)
    in_top = target_index is not None and target_index < rules["top_rank"]
    target = scored[target_index] if target_index is not None else None
    detail = {"top_rank": rules["top_rank"]}
    if species_target and target is not None:
        detail["species_route"] = record_lipid_species(target[0])[1]
    checks.append(_check("target_in_top_rank", "pass" if in_top else "fail",
                         None if target_index is None else target_index + 1, detail))

    if rules["min_fragments_excl_precursor"] is not None:
        tol = rules["precursor_exclusion_tol"] or ms2_tol
        floor = rules["min_fragment_relative_intensity"]
        if target is None:
            checks.append(_check("min_fragments", "not_evaluable", None,
                                 {"reason": "target_not_among_candidates"}))
        else:
            count = fragment_matches_excluding_precursor(
                target[1].get("alignment") or [], precursor_mz, tol, min_relative_intensity=floor)
            status = ("not_evaluable" if count is None
                      else "pass" if count >= rules["min_fragments_excl_precursor"] else "fail")
            checks.append(_check("min_fragments", status, count,
                                 {"minimum": rules["min_fragments_excl_precursor"],
                                  "precursor_exclusion_tol": tol,
                                  "min_relative_intensity": floor}))

    margin = rules["same_formula_margin"]
    if margin is not None:
        if target is None:
            checks.append(_check("same_formula_margin", "not_evaluable", None,
                                 {"reason": "target_not_among_candidates"}))
        else:
            formula = rules["target"]["formula"] or target[0].get("formula")
            alternatives = [(r, s) for r, s in scored
                            if formula and (r.get("formula") or "").strip() == formula
                            and key_of(r) not in (None, block)]
            if not formula:
                checks.append(_check("same_formula_margin", "not_evaluable", None,
                                     {"reason": "formula_unknown"}))
            elif not alternatives:
                checks.append(_check("same_formula_margin", "pass", None,
                                     {"reason": "no_same_formula_alternative"}))
            else:
                best = max(alternatives, key=lambda item: item[1][margin["score"]])
                delta = target[1][margin["score"]] - best[1][margin["score"]]
                checks.append(_check("same_formula_margin",
                                     "pass" if delta >= margin["min_delta"] else "fail", delta,
                                     {"score": margin["score"], "min_delta": margin["min_delta"],
                                      "best_alternative": best[0].get("name"),
                                      "best_alternative_inchikey": best[0].get("inchikey")}))

    if rules["lipid_class"] is not None:
        top_name = scored[0][0].get("name") if scored else None
        found = _lipid_class(top_name)
        status = ("not_evaluable" if found is None
                  else "pass" if found == rules["lipid_class"] else "fail")
        checks.append(_check("lipid_class", status, found,
                             {"expected": rules["lipid_class"], "rank1_name": top_name}))

    for check, ions in (("diagnostic_ions", rules["diagnostic_ions"]),
                        ("neutral_losses", rules["neutral_losses"])):
        if not ions:
            continue
        found = find_ions(measured, ions, precursor_mz=precursor_mz)
        hits = [f["found"] for f in found]
        passed = all(hits) if rules["ions_match"] == "all" else any(hits)
        checks.append(_check(check, "pass" if passed else "fail", sum(hits),
                             {"ions": found, "match": rules["ions_match"]}))

    statuses = {c["status"] for c in checks}
    verdict = ("not_qualified" if "fail" in statuses
               else "not_evaluable" if "not_evaluable" in statuses else "qualified")
    return {"verdict": verdict, "checks": checks,
            "target_rank": None if target_index is None else target_index + 1,
            "top1_is_target": target_index == 0}


# ---------- 候補の行（`library_msms_qualify` の表と `dataset_claim_evidence` の `msms.top`） ----------

#: 候補の行の既定の列。`.dbs` の足切りがあれば `passes_dbs_cutoffs`、種の目標なら
#: `lipid_species`・`species_route` を呼び出し側が足す。
CANDIDATE_COLUMNS = ("rank", "name", "formula", "inchikey", "adduct", "delta_ppm", "total_score",
                     "matched_peaks_count", "fragment_matches_excl_precursor")


def candidate_rows(scored, precursor_mz: float, rules: dict, *, ms2_tol: float,
                   cutoffs: dict | None, top_n: int) -> list[dict]:
    """MS/MS の合否の候補の上位 `top_n` 行。

    `delta_ppm` は `(測定 - 参照) / 参照 * 1e6`（小数 2 桁）。`fragment_matches_excl_precursor` は合否と
    同じ規則（前駆体 ± `precursor_exclusion_tol` を除き、測定強度が `min_fragment_relative_intensity`
    以上の窓）で数える。`.dbs` に注釈の足切りがあれば `passes_dbs_cutoffs`（合否には当てない材料）、
    種の目標なら `lipid_species` と `species_route`（記録の種を得た経路）を足す。
    """
    from metabolomix.library import matching

    tol = rules["precursor_exclusion_tol"] or ms2_tol
    species_target = "lipid_species" in rules["target"]
    rows = []
    for rank, (record, result) in enumerate(scored[:top_n], start=1):
        ppm = matching.delta_ppm(precursor_mz, record.get("precursor_mz"))
        row = {"rank": rank, "name": record.get("name"), "formula": record.get("formula"),
               "inchikey": record.get("inchikey"), "adduct": record.get("adduct"),
               "delta_ppm": None if ppm is None else round(ppm, 2),
               "total_score": result["total_score"],
               "matched_peaks_count": result["matched_peaks_count"],
               "fragment_matches_excl_precursor": fragment_matches_excluding_precursor(
                   result["alignment"], precursor_mz, tol,
                   min_relative_intensity=rules["min_fragment_relative_intensity"])}
        if cutoffs:
            row["passes_dbs_cutoffs"] = matching.passes_dbs_cutoffs(result, cutoffs)
        if species_target:
            row["lipid_species"], row["species_route"] = record_lipid_species(record)
        rows.append(row)
    return rows


def query_window(query: dict) -> dict:
    """照合に使った前駆体の窓（`score_candidates` の `query` から）。`mz_tol` は m/z 500 以下の Da の幅、
    `mz_window` は実際に使った Da の幅（m/z 500 を超えたら ppm で広がる）、`mz_tol_ppm` は重ねた ppm の絞り込み。"""
    return {k: query[k] for k in ("mz_tol", "mz_window", "mz_tol_ppm")}
