"""論理target → バッチ内featureの対応付け（spec §8.1 `feature-bindings.v1`）。

数値・ドメイン層の純ロジックで、MCP にもセッションにも依存しない。

profile が持つのは**規則**だけ——「GABA を `[M+H]+` で、104.0706±10ppm、
1.2±0.1 分に、`mass_rt` の証拠で採る」。どの feature_id がそれに当たるかは
MS-DIAL のアライメント結果ごとに変わるので、バッチごとに対応付け直す。
profile 自身は書き換えない（別バッチの結果が規則に混ざると、profile の
同一性——`profile_content_hash`——が意味を失う）。

## 何を AND で見るか

`m/z ppm` ・ `RT 分` ・ 化合物ID ・ adduct/charge ・ 必要証拠。1つでも落ちた候補は
採らない。落ちた理由は候補ごとに全部残す——「なぜ0件だったのか」が分からないと、
規則が厳しすぎるのかバッチが違うのかを人が判断できない。

**同定情報が無い feature をどう扱うか。** 候補（SME）が1件も付いていない feature は、
化合物ID・adduct/charge を*評価できない*（一致も不一致もしていない）。ここを
「不一致」に倒すと `mass_rt`——ライブラリ照合なしで質量と RT だけで採るという
証拠水準——が原理的に成立しなくなる。そこで:

- `mass_rt` … 同定が無くても質量・RT が合えば qualified（`identity_not_evaluable`
  を理由として残す）。ただし同定が**有って矛盾する**なら不一致として落とす。
- `library_match` / `authentic_standard_match` … 同定・証拠そのものを要求するので、
  評価できない時点で qualified にしない。

## 0件と複数件を混ぜない

条件を満たす候補がちょうど1件のときだけ自動確定する。0件と複数件はどちらも
`needs_input` だが理由を分ける——0件は規則かバッチが違う、複数件は規則が緩いか
異性体が居る、で次の一手が正反対になる。最高強度・先頭候補で選ばない
（それは「選んだ」のではなく「並び順に従った」だけ）。

## 手動 override の範囲

override できるのは **qualified 候補の中からの選択**だけ。理由必須、別 dataset の
override は拒否する。0候補のときに無関係な feature を強制採用させない——手動選択は
同定証拠の不足を解消しない（spec §8.1「手動で候補を選択しても同定証拠の不足を
解消したことにはしない。証拠条件を変更する必要があればprofileを改訂する」）。
"""
from __future__ import annotations

import re

from metabolomix.analysis.result_state import dataset_fingerprint
from metabolomix.core.atomic_io import DomainError, canonical_hash

__all__ = [
    "BINDINGS_SCHEMA",
    "INTERNAL_STANDARD_MAP_SCHEMA",
    "bind_features",
    "resolve_candidates",
    "resolved_targets",
]

BINDINGS_SCHEMA = "feature-bindings.v1"
INTERNAL_STANDARD_MAP_SCHEMA = "internal-standard-map.v1"

_OVERRIDE_CODE = "FEATURE_BINDING_OVERRIDE_INVALID"

#: 名前は識別子として数えない（spec §8.1「名前だけは不可」「証拠欠落を名前一致で
#: 補わない」）。同名の異性体・誘導体が同じ規則に当たるため。
_NON_NAME_IDENTIFIERS = ("inchikey", "formula", "cas", "hmdb_id", "kegg_id",
                         "pubchem_cid")

#: `[M+H]1+`（mzTab-M）と `[M+H]+`（profile の記法）を同じ adduct として扱うための
#: 末尾電荷表記。`1+` / `1-` だけを落とす——`2+` を `+` に潰すと価数が消える。
_ADDUCT_UNIT_CHARGE_RE = re.compile(r"\]1([+-])$")


# ---------- 純関数: 候補の絞り込み ----------

def resolve_candidates(candidates: list[dict]) -> dict:
    """qualified な候補がちょうど1件のときだけ確定する。

    0件と複数件はどちらも `needs_input` だが `reason` で区別する。呼び出し側が
    件数から推測しなくて済むように、理由を封筒に載せる。
    """
    valid = [c for c in candidates if c.get("qualified")]
    if len(valid) == 1:
        return {"status": "resolved", "selected": valid[0]["feature_id"],
                "reason": None, "candidates": candidates}
    reason = ("multiple_qualified_candidates" if len(valid) > 1
              else "no_qualified_candidate")
    return {"status": "needs_input", "selected": None, "reason": reason,
            "candidates": candidates}


# ---------- 照合の部品 ----------

def _normalize_adduct(value) -> str | None:
    if value is None:
        return None
    text = str(value).strip().replace(" ", "")
    return _ADDUCT_UNIT_CHARGE_RE.sub(r"]\1", text).casefold()


def _normalize_identifier(value) -> str | None:
    """`HMDB:HMDB0000112` と `HMDB0000112` を同じ値として比較するための正規化。"""
    if value is None:
        return None
    text = str(value).strip()
    if ":" in text:
        text = text.rsplit(":", 1)[-1].strip()
    return text.casefold() or None


def _identifier_values(feature_meta: dict, candidate: dict | None) -> set[str]:
    values = set()
    for raw in (feature_meta.get("inchikey"),
                (candidate or {}).get("database_identifier"),
                (candidate or {}).get("inchi"),
                (candidate or {}).get("smiles")):
        normalized = _normalize_identifier(raw)
        if normalized:
            values.add(normalized)
    return values


def _ppm_delta(observed, expected) -> float | None:
    if observed is None or expected in (None, 0):
        return None
    return abs(float(observed) - float(expected)) / float(expected) * 1e6


def _check_identity(target: dict, feature_meta: dict,
                    candidates: list[dict]) -> tuple[str, dict | None, list[str]]:
    """化合物ID・adduct・charge を評価する。

    戻り値は `(status, matched_candidate, reasons)` で status は
    `matched` / `not_evaluable`（同定が無い） / `mismatch`。
    """
    wanted = {key: _normalize_identifier(value)
              for key, value in (target.get("compound_identifiers") or {}).items()
              if key in _NON_NAME_IDENTIFIERS and value is not None}
    feature_values = _identifier_values(feature_meta, None)

    if not candidates:
        # 同定が1件も無い。名前で埋めないので、ここでは「評価できない」に留める。
        if wanted and feature_values & set(wanted.values()):
            return "matched", None, []
        return "not_evaluable", None, ["identity_not_evaluable"]

    reasons: list[str] = []
    informative = 0
    for candidate in candidates:
        values = _identifier_values(feature_meta, candidate)
        if wanted and not values:
            # 候補は付いているが識別子を1つも持たない（名前だけ、あるいは
            # 全部null）。これは「違う化合物だ」という証拠ではなく「この候補
            # では判定できない」——不一致に倒すと、同定が空の成果物で
            # mass_rt が原理的に成立しなくなる。
            continue
        informative += 1
        if wanted and not (values & set(wanted.values())):
            reasons.append("identifier_mismatch")
            continue
        if (_normalize_adduct(candidate.get("adduct"))
                != _normalize_adduct(target.get("adduct"))):
            reasons.append("adduct_mismatch")
            continue
        if candidate.get("charge") != target.get("charge"):
            reasons.append("charge_mismatch")
            continue
        return "matched", candidate, []
    if not informative:
        # 判定材料を持つ候補が1件も無い。「同定が無い」と同じ扱いにする。
        return "not_evaluable", None, ["identity_not_evaluable"]
    # どの候補も通らなかった。理由は出た順に一意化して返す（件数ではなく種類が要る）。
    return "mismatch", None, list(dict.fromkeys(reasons))


def _check_library_match(spec: dict, candidate: dict | None,
                         observed_hashes: dict) -> list[str]:
    """ライブラリ照合の証拠。source hash と score 条件の**両方**を要求する。"""
    if candidate is None:
        return ["identification_missing"]

    reasons: list[str] = []
    library_id = spec.get("library_id")
    expected_hash = spec.get("library_sha256")
    actual_hash = (observed_hashes or {}).get(library_id)
    if actual_hash is None:
        # 「このバッチがそのライブラリで処理された」ことを確認できない。
        # 未確認を合格扱いにすると、別ライブラリの同定が profile の証拠水準を
        # 満たしたことになる。
        reasons.append("library_hash_unverified")
    elif actual_hash != expected_hash:
        reasons.append("library_hash_mismatch")

    score_field = spec.get("score_field")
    threshold = spec.get("score_threshold")
    if candidate.get("confidence_measure") != score_field:
        reasons.append("library_score_missing")
    else:
        value = candidate.get("confidence_value")
        if value is None:
            reasons.append("library_score_missing")
        elif threshold is not None and float(value) < float(threshold):
            reasons.append("library_score_below_threshold")
    return reasons


def _check_standard_match(spec: dict, target: dict, feature_id: str,
                          evidence: dict, assay_ids: list[str],
                          refs: list[dict]) -> list[str]:
    """標準品注入の実測 RT/m/z を照合する（注入ごとの行だけを使う）。"""
    if not assay_ids:
        return ["standard_assays_missing"]
    if not (evidence or {}).get("availability"):
        return ["assay_evidence_unavailable"]

    rows = {(r.get("feature_id"), r.get("assay_id")): r
            for r in (evidence.get("rows") or [])}
    reasons: list[str] = []
    for assay_id in assay_ids:
        row = rows.get((feature_id, assay_id))
        if row is None:
            reasons.append("standard_evidence_row_missing")
            continue
        if row.get("detection_status") != "detected":
            # gap-fill は補間値であって観測ではない。標準品の一致の根拠にしない。
            reasons.append("standard_injection_not_detected")
            continue
        ppm = _ppm_delta(row.get("observed_mz"), target.get("expected_mz"))
        if ppm is None or ppm > float(spec["mz_tolerance_ppm"]):
            reasons.append("standard_mz_outside_tolerance")
            continue
        rt = row.get("observed_rt_min")
        if rt is None or abs(float(rt) - float(target["expected_rt_min"])) \
                > float(spec["rt_tolerance_min"]):
            reasons.append("standard_rt_outside_tolerance")
            continue
        refs.append({"feature_id": feature_id, "assay_id": assay_id,
                     "observed_rt_min": rt, "observed_mz": row.get("observed_mz")})
    return list(dict.fromkeys(reasons))


def _evaluate_feature(feature_id: str, target: dict, ds, evidence: dict,
                      assay_ids: list[str], observed_hashes: dict) -> dict:
    meta = (getattr(ds, "feature_metadata", {}) or {}).get(feature_id) or {}
    candidates = (getattr(ds, "feature_candidates", {}) or {}).get(feature_id) or []
    reasons: list[str] = []
    refs: list[dict] = []

    ppm = _ppm_delta(meta.get("mz"), target.get("expected_mz"))
    if ppm is None:
        reasons.append("mz_missing")
    elif ppm > float(target["mz_tolerance_ppm"]):
        reasons.append("mz_outside_tolerance")

    rt = meta.get("rt")
    if rt is None:
        reasons.append("rt_missing")
    elif abs(float(rt) - float(target["expected_rt_min"])) \
            > float(target["rt_tolerance_min"]):
        reasons.append("rt_outside_tolerance")

    identity_status, matched, identity_reasons = _check_identity(
        target, meta, candidates)
    reasons.extend(identity_reasons)

    for spec in target.get("required_evidence") or []:
        kind = spec.get("kind")
        if kind == "mass_rt":
            # 質量・RT は上で見た。同定の欠落はこの証拠水準では許容する
            # （矛盾する同定があれば identity 側で既に落ちている）。
            continue
        if identity_status != "matched":
            reasons.append("identification_required")
            continue
        if kind == "library_match":
            reasons.extend(_check_library_match(spec, matched, observed_hashes))
        elif kind == "authentic_standard_match":
            reasons.extend(_check_standard_match(
                spec, target, feature_id, evidence, assay_ids, refs))

    reasons = [r for r in dict.fromkeys(reasons)]
    # `identity_not_evaluable` は理由として残すが、qualified を止めるのは
    # 証拠水準側（`identification_required`）だけ。
    blocking = [r for r in reasons if r != "identity_not_evaluable"]
    return {
        "feature_id": feature_id,
        "qualified": not blocking,
        "reasons": reasons,
        "evidence_refs": refs,
        "observed_mz": meta.get("mz"),
        "observed_rt_min": meta.get("rt"),
        "mz_delta_ppm": None if ppm is None else round(ppm, 3),
        "matched_sme_id": (matched or {}).get("sme_id"),
    }


# ---------- override ----------

def _fail_override(message: str, **details):
    raise DomainError(_OVERRIDE_CODE, message, details)


def _apply_override(target_id: str, override: dict, resolved: dict, ds) -> dict:
    reason = (override or {}).get("reason")
    if not reason or not str(reason).strip():
        _fail_override(
            f"手動選択には理由が必要です（target_id={target_id!r}）。",
            target_id=target_id)

    dataset_id = (override or {}).get("dataset_id")
    if dataset_id != getattr(ds, "dataset_id", None):
        # 別バッチで決めた対応を、同じ feature_id が別の化合物を指す今のバッチへ
        # 持ち込ませない（spec §8.1「別datasetのmapは入力として拒否」と同じ理由）。
        _fail_override(
            f"別データセットのoverrideです（target_id={target_id!r}）。"
            "同じprofile規則でこのバッチに対応付け直してください。",
            target_id=target_id, override_dataset_id=dataset_id,
            dataset_id=getattr(ds, "dataset_id", None))

    feature_id = (override or {}).get("feature_id")
    qualified = {c["feature_id"] for c in resolved["candidates"] if c["qualified"]}
    if feature_id not in qualified:
        _fail_override(
            f"qualifiedでない候補は選べません（target_id={target_id!r}, "
            f"feature_id={feature_id!r}）。手動選択は同定証拠の不足を解消しません"
            "——証拠条件を変える必要があるならprofileを改訂してください。",
            target_id=target_id, feature_id=feature_id,
            qualified_feature_ids=sorted(qualified))

    return {"status": "resolved", "selected_feature_id": feature_id,
            "selection": "manual", "selection_reason": str(reason),
            "reason": None}


# ---------- 公開関数 ----------

def bind_features(ds, profile: dict, evidence: dict,
                  standard_assays: dict, overrides: dict | None) -> dict:
    """profile の `feature_targets` を、このバッチの feature_id へ対応付ける。

    `standard_assays` は `target_id -> [assay_id, ...]`。要求の sample_id から
    assay_id への解決（manifest 経由）は**上流で済ませてある**契約で、ここは
    受け取った assay_id が実在することだけを確認する——存在しない assay を黙って
    読み飛ばすと、標準品の照合が「0注入で全部合格」になる。

    `overrides` は `target_id -> {feature_id, reason, dataset_id}`。
    """
    targets = dict((profile.get("feature_targets") or {}))
    recipe = profile.get("analysis_recipe") or {}

    known_assays = set(getattr(ds, "sample_assay_ids", []) or [])
    for target_id, assay_ids in (standard_assays or {}).items():
        unknown = [a for a in (assay_ids or []) if a not in known_assays]
        if unknown:
            raise DomainError(
                "STANDARD_ASSAY_INVALID",
                f"standard_assaysがこのデータセットに無いassayを指しています"
                f"（target_id={target_id!r}）: {', '.join(unknown)}",
                {"target_id": target_id, "unknown_assays": unknown,
                 "known_assays": sorted(known_assays)},
            )

    for target_id in (overrides or {}):
        if target_id not in targets:
            _fail_override(
                f"profileに無いtargetへのoverrideです: {target_id!r}",
                target_id=target_id, known_target_ids=sorted(targets))

    observed_hashes = dict(getattr(ds, "processing_dependencies", {}) or {})
    feature_ids = list(getattr(ds, "feature_ids", []) or [])

    bindings: dict[str, dict] = {}
    unresolved: list[str] = []
    for target_id, target in targets.items():
        assay_ids = list((standard_assays or {}).get(target_id) or [])
        candidates = [_evaluate_feature(fid, target, ds, evidence, assay_ids,
                                        observed_hashes)
                      for fid in feature_ids]
        resolved = resolve_candidates(candidates)

        entry = {"status": resolved["status"],
                 "reason": resolved["reason"],
                 "selected_feature_id": resolved["selected"],
                 "selection": "automatic" if resolved["selected"] else None,
                 "selection_reason": None,
                 "candidates": candidates}
        if target_id in (overrides or {}):
            entry.update(_apply_override(target_id, overrides[target_id],
                                         resolved, ds))
        selected = entry["selected_feature_id"]
        entry["evidence_refs"] = next(
            (c["evidence_refs"] for c in candidates if c["feature_id"] == selected),
            [])
        if entry["status"] != "resolved":
            unresolved.append(target_id)
        bindings[target_id] = entry

    dataset_hash = dataset_fingerprint(ds)
    # 規則hashは profile 側だけから作る。解決結果を混ぜると、同じ profile が
    # バッチごとに別の規則に見え、profile の同一性検査が使えなくなる。
    rule_hash = canonical_hash({
        "feature_targets": targets,
        "internal_standards": recipe.get("internal_standards") or [],
    })

    pairs = []
    for mapping in (recipe.get("internal_standards") or []):
        target_id = mapping.get("target_id")
        standard_target_id = mapping.get("standard_target_id")
        target_binding = bindings.get(target_id) or {}
        standard_binding = bindings.get(standard_target_id) or {}
        if (target_binding.get("status") != "resolved"
                or standard_binding.get("status") != "resolved"):
            # 片方でも未解決なら対応を作らない。部分的なmapは
            # 「補正済み」と「未補正」を同じ行列へ混ぜる入口になる。
            continue
        pairs.append({
            "target_id": target_id,
            "target_feature_id": target_binding["selected_feature_id"],
            "standard_target_id": standard_target_id,
            "standard_feature_id": standard_binding["selected_feature_id"],
        })

    return {
        "schema": BINDINGS_SCHEMA,
        "dataset_id": getattr(ds, "dataset_id", None),
        "dataset_hash": dataset_hash,
        "rule_hash": rule_hash,
        "status": "resolved" if not unresolved else "needs_input",
        "bindings": bindings,
        "unresolved": unresolved,
        "internal_standard_map": {
            "schema": INTERNAL_STANDARD_MAP_SCHEMA,
            "dataset_id": getattr(ds, "dataset_id", None),
            "dataset_hash": dataset_hash,
            "rule_hash": rule_hash,
            "pairs": pairs,
        },
    }


def resolved_targets(payload: dict | None) -> dict:
    """`feature-bindings.v1` から「解決済み target → feature」だけを取り出す。

    worker（統計の feature_scope 解決）と、run を対話セッションへ引き渡す経路が
    共有する。未解決（候補が複数・不在）の target は**含めない**——解決できて
    いない対象を統計の対象集合へ混ぜない。
    """
    bindings = (payload or {}).get("bindings") or {}
    return {target_id: entry["selected_feature_id"]
            for target_id, entry in bindings.items()
            if entry.get("status") == "resolved"}
