"""論理target → バッチ内featureの対応付け（spec §8.1 `feature-bindings.v1`）。

profileは「GABAを[M+H]+で、104.0706±10ppm、1.2±0.1分に、mass_rtの証拠で採る」と
いう**規則**だけを持つ。どのfeature_idがそれに当たるかはバッチごとに変わるので、
解析のたびに対応付け直す。ここで縛るのは:

1. 許容範囲・化合物ID・adduct/charge・必要証拠を**ANDで**評価する。1つでも
   落ちた候補は採らない。
2. 条件を満たす候補がちょうど1件のときだけ自動確定する。0件と複数件はどちらも
   `needs_input` だが、**理由は区別する**——0件は「規則かバッチが違う」、
   複数件は「規則が緩すぎる/異性体が居る」で、次の一手が正反対になる。
   最高強度や先頭候補で選ばない。
3. 手動overrideは**qualified候補の選択**に限る。理由必須、別datasetのoverrideは
   拒否。0候補のときに無関係なfeatureを強制採用させない——手動選択は同定証拠の
   不足を解消しない（解消したければprofileを改訂する）。

fixtureは全て合成。
"""
from __future__ import annotations

import copy

import pytest

from metabolomix.analysis.feature_bindings import (
    BINDINGS_SCHEMA,
    bind_features,
    resolve_candidates,
)
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab.dataset_state import DatasetState


def _target(**overrides) -> dict:
    target = {
        "compound_identifiers": {"name": "GABA",
                                 "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N"},
        "adduct": "[M+H]+",
        "charge": 1,
        "expected_mz": 104.0706,
        "mz_tolerance_ppm": 10.0,
        "expected_rt_min": 1.2,
        "rt_tolerance_min": 0.1,
        "required_evidence": [{"kind": "mass_rt"}],
    }
    target.update(overrides)
    return target


def _profile(targets: dict | None = None, internal_standards=None) -> dict:
    return {
        "schema": "lcms-profile.v1",
        "profile_id": "p1",
        "revision": 1,
        "feature_targets": targets if targets is not None else {"gaba": _target()},
        "analysis_recipe": {"statistics": [],
                            "internal_standards": internal_standards or []},
    }


def _candidate(**overrides) -> dict:
    candidate = {
        "sme_id": "1", "rank": 1, "charge": 1,
        "exp_mass_to_charge": 104.0706, "theoretical_mass_to_charge": 104.0706,
        "confidence_value": 0.92, "confidence_measure": "MS-DIAL total score",
        "chemical_name": "GABA", "database_identifier": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
        "inchi": None, "smiles": None, "adduct": "[M+H]1+",
        "identification_method": "matched to library", "spectra_ref": None,
    }
    candidate.update(overrides)
    return candidate


def _dataset(features: dict, *, assays=("assay[1]", "assay[2]")) -> DatasetState:
    """`{feature_id: {"mz":, "rt":, "candidates": [...]}}` から最小のdsを作る。"""
    ds = DatasetState()
    ds.feature_ids = list(features)
    ds.sample_assay_ids = list(assays)
    ds.sample_names = [f"S{i + 1}" for i in range(len(assays))]
    for fid, spec in features.items():
        ds.feature_metadata[fid] = {
            "mz": spec.get("mz"), "rt": spec.get("rt"),
            "name": spec.get("name"), "inchikey": spec.get("inchikey"),
        }
        ds.feature_candidates[fid] = list(spec.get("candidates", []))
    return ds


def _one_good_feature() -> DatasetState:
    return _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                           "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                           "candidates": [_candidate()]}})


def _evidence(ds, rows=()) -> dict:
    return {"schema": "assay-feature-evidence.v1", "dataset_id": ds.dataset_id,
            "availability": bool(rows), "reasons": [], "rows": list(rows)}


def _row(feature_id, assay_id, *, rt=1.2, mz=104.0706, status="detected") -> dict:
    return {"feature_id": feature_id, "assay_id": assay_id,
            "observed_rt_min": rt, "observed_mz": mz, "detection_status": status}


# ---------- resolve_candidates ----------

def test_two_qualified_candidates_need_input():
    out = resolve_candidates([{"feature_id": "1", "qualified": True},
                              {"feature_id": "2", "qualified": True}])
    assert out["status"] == "needs_input"
    assert out["selected"] is None


def test_zero_and_multiple_have_different_reasons():
    none_qualified = resolve_candidates([{"feature_id": "1", "qualified": False}])
    two_qualified = resolve_candidates([{"feature_id": "1", "qualified": True},
                                        {"feature_id": "2", "qualified": True}])
    assert none_qualified["status"] == two_qualified["status"] == "needs_input"
    assert none_qualified["reason"] != two_qualified["reason"]
    assert none_qualified["reason"] == "no_qualified_candidate"
    assert two_qualified["reason"] == "multiple_qualified_candidates"


def test_single_qualified_candidate_resolves():
    out = resolve_candidates([{"feature_id": "1", "qualified": True},
                              {"feature_id": "2", "qualified": False}])
    assert out == {"status": "resolved", "selected": "1", "reason": None,
                   "candidates": [{"feature_id": "1", "qualified": True},
                                  {"feature_id": "2", "qualified": False}]}


# ---------- bind_features: 自動確定 ----------

def test_single_match_binds_automatically():
    ds = _one_good_feature()
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)

    assert result["schema"] == BINDINGS_SCHEMA
    assert result["status"] == "resolved"
    binding = result["bindings"]["gaba"]
    assert binding["selected_feature_id"] == "1"
    assert binding["selection"] == "automatic"
    assert binding["selection_reason"] is None


def test_the_same_rules_bind_differently_in_two_batches():
    """規則hashは同じ、dataset hashと採用feature_idはバッチごとに変わる。"""
    profile = _profile()
    batch_a = _one_good_feature()
    batch_b = _dataset({"77": {"mz": 104.0707, "rt": 1.19,
                               "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                               "candidates": [_candidate()]}})

    a = bind_features(batch_a, profile, _evidence(batch_a), {}, None)
    b = bind_features(batch_b, profile, _evidence(batch_b), {}, None)

    assert a["rule_hash"] == b["rule_hash"]
    assert a["dataset_hash"] != b["dataset_hash"]
    assert a["bindings"]["gaba"]["selected_feature_id"] == "1"
    assert b["bindings"]["gaba"]["selected_feature_id"] == "77"


def test_rule_hash_ignores_the_binding_result():
    """解決結果を規則hashへ混ぜない（混ぜると同じprofileが別物に見える）。"""
    profile = _profile()
    ds = _one_good_feature()
    empty = _dataset({"9": {"mz": 999.0, "rt": 9.0}})

    resolved = bind_features(ds, profile, _evidence(ds), {}, None)
    unresolved = bind_features(empty, profile, _evidence(empty), {}, None)
    assert resolved["rule_hash"] == unresolved["rule_hash"]
    assert unresolved["status"] == "needs_input"


# ---------- bind_features: 落ちる条件 ----------

def test_mz_outside_ppm_window_is_not_a_candidate():
    ds = _dataset({"1": {"mz": 104.5, "rt": 1.2,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate()]}})
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "needs_input"
    assert binding["reason"] == "no_qualified_candidate"
    assert "mz_outside_tolerance" in binding["candidates"][0]["reasons"]


def test_rt_outside_window_is_not_a_candidate():
    ds = _dataset({"1": {"mz": 104.0706, "rt": 4.0,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate()]}})
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    assert "rt_outside_tolerance" in result["bindings"]["gaba"]["candidates"][0]["reasons"]


def test_adduct_mismatch_disqualifies():
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate(adduct="[M+Na]1+")]}})
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    assert "adduct_mismatch" in result["bindings"]["gaba"]["candidates"][0]["reasons"]


def test_charge_mismatch_disqualifies():
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate(charge=2)]}})
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    assert "charge_mismatch" in result["bindings"]["gaba"]["candidates"][0]["reasons"]


def test_adduct_charge_notation_difference_is_not_a_mismatch():
    """profileの`[M+H]+`とmzTabの`[M+H]1+`は同じadduct。表記差で落とさない。"""
    ds = _one_good_feature()
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    assert result["bindings"]["gaba"]["selected_feature_id"] == "1"


def test_name_only_match_is_never_identity_evidence():
    """名前だけ一致する候補を同定証拠にしない（spec §8.1）。

    `mass_rt` は質量とRTを証拠水準にするので、同定が空でも対応付け自体は成立
    してよい——ただしそれは「同定できた」ことではないので、理由に
    `identity_not_evaluable` を残す。名前一致で証拠の欠落を埋めない、とは
    「名前が一致したから同定済みと見なさない」という意味であって、
    「名前しか無ければ質量とRTでも採ってはいけない」ではない。
    """
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "candidates": [_candidate(database_identifier=None,
                                                   chemical_name="GABA")]}})
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "resolved"
    assert "identity_not_evaluable" in binding["candidates"][0]["reasons"]
    assert binding["candidates"][0]["matched_sme_id"] is None



def test_same_name_isomers_both_qualify_and_need_input():
    """同名異性体が2件。強度や順序で選ばず、人に返す。"""
    ds = _dataset({
        "1": {"mz": 104.0706, "rt": 1.18,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "candidates": [_candidate()]},
        "2": {"mz": 104.0707, "rt": 1.22,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "candidates": [_candidate(sme_id="2")]},
    })
    result = bind_features(ds, _profile(), _evidence(ds), {}, None)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "needs_input"
    assert binding["reason"] == "multiple_qualified_candidates"
    assert binding["selected_feature_id"] is None
    assert sorted(c["feature_id"] for c in binding["candidates"] if c["qualified"]) \
        == ["1", "2"]


# ---------- 必要証拠 ----------

def _library_target() -> dict:
    return _target(required_evidence=[{
        "kind": "library_match", "library_id": "msp-lib-1",
        "library_sha256": "c" * 64, "score_field": "MS-DIAL total score",
        "score_threshold": 0.8}])


def test_name_only_match_does_not_satisfy_a_library_requirement():
    """同定そのものを要求する証拠水準では、名前だけの候補は通らない。"""
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "candidates": [_candidate(database_identifier=None,
                                                   chemical_name="GABA")]}})
    ds.processing_dependencies = {"msp-lib-1": "c" * 64}
    result = bind_features(ds, _profile({"gaba": _library_target()}),
                           _evidence(ds), {}, None)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "needs_input"
    assert "identification_required" in binding["candidates"][0]["reasons"]


def test_library_match_requires_the_declared_library_hash():
    ds = _one_good_feature()
    ds.processing_dependencies = {"msp-lib-1": "d" * 64}       # 別のライブラリ
    result = bind_features(ds, _profile({"gaba": _library_target()}),
                           _evidence(ds), {}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "library_hash_mismatch" in reasons


def test_library_match_without_any_recorded_hash_is_unverified():
    ds = _one_good_feature()
    result = bind_features(ds, _profile({"gaba": _library_target()}),
                           _evidence(ds), {}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "library_hash_unverified" in reasons


def test_library_match_below_threshold_is_not_evidence():
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate(confidence_value=0.5)]}})
    ds.processing_dependencies = {"msp-lib-1": "c" * 64}
    result = bind_features(ds, _profile({"gaba": _library_target()}),
                           _evidence(ds), {}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "library_score_below_threshold" in reasons


def test_library_match_with_hash_and_score_qualifies():
    ds = _one_good_feature()
    ds.processing_dependencies = {"msp-lib-1": "c" * 64}
    result = bind_features(ds, _profile({"gaba": _library_target()}),
                           _evidence(ds), {}, None)
    assert result["bindings"]["gaba"]["selected_feature_id"] == "1"


def _standard_target() -> dict:
    return _target(required_evidence=[{
        "kind": "authentic_standard_match",
        "mz_tolerance_ppm": 10.0, "rt_tolerance_min": 0.1}])


def test_authentic_standard_match_without_standard_assays_is_unresolved():
    """profileが標準注入を要求するのに、そのroleが無いバッチは未解決。"""
    ds = _one_good_feature()
    result = bind_features(ds, _profile({"gaba": _standard_target()}),
                           _evidence(ds, [_row("1", "assay[1]")]), {}, None)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "needs_input"
    assert "standard_assays_missing" in binding["candidates"][0]["reasons"]


def test_authentic_standard_match_uses_per_injection_rows():
    ds = _one_good_feature()
    rows = [_row("1", "assay[1]", rt=1.21, mz=104.0709),
            _row("1", "assay[2]", rt=1.19, mz=104.0703)]
    result = bind_features(ds, _profile({"gaba": _standard_target()}),
                           _evidence(ds, rows),
                           {"gaba": ["assay[1]", "assay[2]"]}, None)
    binding = result["bindings"]["gaba"]
    assert binding["selected_feature_id"] == "1"
    assert len(binding["evidence_refs"]) == 2


def test_gap_filled_standard_injection_is_not_evidence():
    ds = _one_good_feature()
    rows = [_row("1", "assay[1]", status="gap_filled"),
            _row("1", "assay[2]")]
    result = bind_features(ds, _profile({"gaba": _standard_target()}),
                           _evidence(ds, rows),
                           {"gaba": ["assay[1]", "assay[2]"]}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "standard_injection_not_detected" in reasons


def test_standard_rt_outside_its_own_tolerance_is_not_evidence():
    ds = _one_good_feature()
    rows = [_row("1", "assay[1]", rt=2.5), _row("1", "assay[2]")]
    result = bind_features(ds, _profile({"gaba": _standard_target()}),
                           _evidence(ds, rows),
                           {"gaba": ["assay[1]", "assay[2]"]}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "standard_rt_outside_tolerance" in reasons


def test_unavailable_evidence_blocks_standard_match():
    ds = _one_good_feature()
    result = bind_features(ds, _profile({"gaba": _standard_target()}),
                           _evidence(ds), {"gaba": ["assay[1]"]}, None)
    reasons = result["bindings"]["gaba"]["candidates"][0]["reasons"]
    assert "assay_evidence_unavailable" in reasons


def test_unknown_standard_assay_id_is_rejected():
    ds = _one_good_feature()
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile({"gaba": _standard_target()}),
                      _evidence(ds, [_row("1", "assay[1]")]),
                      {"gaba": ["assay[9]"]}, None)
    assert caught.value.code == "STANDARD_ASSAY_INVALID"


# ---------- override ----------

def test_override_selects_among_qualified_candidates():
    ds = _dataset({
        "1": {"mz": 104.0706, "rt": 1.18,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "candidates": [_candidate()]},
        "2": {"mz": 104.0707, "rt": 1.22,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "candidates": [_candidate(sme_id="2")]},
    })
    overrides = {"gaba": {"feature_id": "2", "reason": "標準品の保持時間と一致",
                          "dataset_id": ds.dataset_id}}
    result = bind_features(ds, _profile(), _evidence(ds), {}, overrides)
    binding = result["bindings"]["gaba"]
    assert binding["status"] == "resolved"
    assert binding["selected_feature_id"] == "2"
    assert binding["selection"] == "manual"
    assert binding["selection_reason"] == "標準品の保持時間と一致"


def test_override_cannot_force_an_unqualified_feature():
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2,
                         "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                         "candidates": [_candidate()]},
                   "9": {"mz": 999.0, "rt": 9.0}})
    overrides = {"gaba": {"feature_id": "9", "reason": "見た目が近い",
                          "dataset_id": ds.dataset_id}}
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile(), _evidence(ds), {}, overrides)
    assert caught.value.code == "FEATURE_BINDING_OVERRIDE_INVALID"


def test_override_without_a_reason_is_rejected():
    ds = _one_good_feature()
    overrides = {"gaba": {"feature_id": "1", "dataset_id": ds.dataset_id}}
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile(), _evidence(ds), {}, overrides)
    assert caught.value.code == "FEATURE_BINDING_OVERRIDE_INVALID"


def test_override_from_another_dataset_is_rejected():
    ds = _one_good_feature()
    overrides = {"gaba": {"feature_id": "1", "reason": "前のバッチで確認済み",
                          "dataset_id": "ds_somewhere_else"}}
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile(), _evidence(ds), {}, overrides)
    assert caught.value.code == "FEATURE_BINDING_OVERRIDE_INVALID"
    assert caught.value.details["target_id"] == "gaba"


def test_override_for_an_unknown_target_is_rejected():
    ds = _one_good_feature()
    overrides = {"not_a_target": {"feature_id": "1", "reason": "x",
                                  "dataset_id": ds.dataset_id}}
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile(), _evidence(ds), {}, overrides)
    assert caught.value.code == "FEATURE_BINDING_OVERRIDE_INVALID"


def test_manual_selection_does_not_add_missing_identification_evidence():
    """手動選択は候補の選択であって、証拠条件の免除ではない。"""
    ds = _one_good_feature()
    ds.processing_dependencies = {"msp-lib-1": "d" * 64}       # hash不一致
    overrides = {"gaba": {"feature_id": "1", "reason": "目視で確認",
                          "dataset_id": ds.dataset_id}}
    with pytest.raises(DomainError) as caught:
        bind_features(ds, _profile({"gaba": _library_target()}),
                      _evidence(ds), {}, overrides)
    assert caught.value.code == "FEATURE_BINDING_OVERRIDE_INVALID"


# ---------- 内部標準map ----------

def test_internal_standard_map_is_built_from_resolved_bindings():
    ds = _dataset({
        "1": {"mz": 104.0706, "rt": 1.2,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "candidates": [_candidate()]},
        "2": {"mz": 110.11, "rt": 1.2,
              "inchikey": "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",
              "candidates": [_candidate(sme_id="2", exp_mass_to_charge=110.11,
                                        database_identifier="ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",
                                        chemical_name="GABA-d6")]},
    })
    d6 = _target(compound_identifiers={"name": "GABA-d6",
                                       "inchikey": "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"},
                 expected_mz=110.11)
    profile = _profile({"gaba": _target(), "gaba_d6": d6},
                       internal_standards=[{"target_id": "gaba",
                                            "standard_target_id": "gaba_d6"}])
    result = bind_features(ds, profile, _evidence(ds), {}, None)

    assert result["status"] == "resolved"
    assert result["internal_standard_map"]["schema"] == "internal-standard-map.v1"
    assert result["internal_standard_map"]["pairs"] == [
        {"target_id": "gaba", "target_feature_id": "1",
         "standard_target_id": "gaba_d6", "standard_feature_id": "2"}]
    assert result["internal_standard_map"]["dataset_hash"] == result["dataset_hash"]


def test_unresolved_standard_leaves_no_internal_standard_pair():
    ds = _one_good_feature()
    d6 = _target(compound_identifiers={"name": "GABA-d6",
                                       "inchikey": "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"},
                 expected_mz=110.11)
    profile = _profile({"gaba": _target(), "gaba_d6": d6},
                       internal_standards=[{"target_id": "gaba",
                                            "standard_target_id": "gaba_d6"}])
    result = bind_features(ds, profile, _evidence(ds), {}, None)
    assert result["status"] == "needs_input"
    assert result["internal_standard_map"]["pairs"] == []
    assert "gaba_d6" in result["unresolved"]


def test_profile_is_not_mutated():
    ds = _one_good_feature()
    profile = _profile()
    before = copy.deepcopy(profile)
    bind_features(ds, profile, _evidence(ds), {}, None)
    assert profile == before


# ---------- SML 由来の注釈は証拠にならない（spec 2026-09-17 §5）----------

def _ms1_annotated_dataset():
    """m/z と RT は合うが、同定は SML 由来のラベルしか無い dataset。"""
    ds = _dataset({"1": {"mz": 104.0706, "rt": 1.2, "candidates": []}})
    ds.feature_annotations = {
        "1": {"sml_id": "1", "ambiguous": False, "name": "GABA",
              "database_identifier": "TextDB:GABA", "chemical_formula": None,
              "smiles": None, "adduct": "[M+H]1+", "reliability": None,
              "confidence_measure": "[,, MS-DIAL algorithm matching score, ]",
              "confidence_value": 0.999982,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "inchikey_source": "smiles_derived"},
    }
    return ds


def test_an_ms1_annotation_does_not_satisfy_a_library_match():
    """SML の score があっても library_match（スペクトル照合）は満たさない。"""
    ds = _ms1_annotated_dataset()
    target = _target(required_evidence=[{
        "kind": "library_match", "library_id": "lib1",
        "library_sha256": "0" * 64,
        "score_field": "[,, MS-DIAL algorithm matching score, ]",
        "score_threshold": 0.8}])

    out = bind_features(ds, _profile({"gaba": target}), _evidence(ds), {}, None)

    entry = out["bindings"]["gaba"]
    assert entry["status"] == "needs_input"
    # 同定が「評価できない」ままであること（matched へ昇格していない）。
    reasons = entry["candidates"][0]["reasons"]
    assert "identity_not_evaluable" in reasons
    assert "identification_required" in reasons


def test_an_ms1_annotation_does_not_unlock_an_authentic_standard_match():
    """ここが緩むと MS1 注釈だけで標準品照合が通る。spec §5 の核心。"""
    ds = _ms1_annotated_dataset()
    target = _target(required_evidence=[{
        "kind": "authentic_standard_match",
        "mz_tolerance_ppm": 10.0, "rt_tolerance_min": 0.1}])
    evidence = _evidence(ds, rows=[_row("1", "assay[1]")])

    out = bind_features(ds, _profile({"gaba": target}),
                        evidence, {"gaba": ["assay[1]"]}, None)

    entry = out["bindings"]["gaba"]
    assert entry["status"] == "needs_input"
    assert "identification_required" in entry["candidates"][0]["reasons"]


def test_an_ms1_annotation_still_allows_a_mass_rt_binding():
    """証拠水準が mass_rt だけなら、同定の有無に関係なく成立する（従来どおり）。"""
    ds = _ms1_annotated_dataset()

    out = bind_features(ds, _profile(), _evidence(ds), {}, None)

    assert out["bindings"]["gaba"]["status"] == "resolved"
