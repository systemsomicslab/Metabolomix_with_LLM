"""統計エンジンの結果の登録と、2 群の結果への読み替え（metabolomix.analysis.stats_results）。

fixture は DatasetState を直に作り、試料の表は `sample_manifest.apply_metadata` で載せる
（`sample_ids` が表示名と違う＝エンジンの試料 ID と表示名がずれる状況を必ず作る）。
"""
from __future__ import annotations

import copy
import json
import math

import numpy as np
import pytest

from metabolomix.analysis import export_contract, result_state, stats_results
from metabolomix.analysis.dataset_exclusions import apply_exclusions
from metabolomix.analysis.result_state import dataset_fingerprint, invalidate_results
from metabolomix.analysis.sample_manifest import apply_metadata
from metabolomix.analysis.stats import data as stats_data
from metabolomix.analysis.stats.pca import pca
from metabolomix.analysis.stats.spec_v2 import normalize_spec_v2, prepare_spec_v2, run_spec_v2
from metabolomix.core.atomic_io import DomainError, canonical_hash
from metabolomix.core.serialization import json_payload, json_safe
from metabolomix.mztab.dataset_state import DatasetState

OFF = {"drift": {"mode": "off"},
       "filters": {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None},
       "normalize": {"method": "none"}}
T_VS_C = {"id": "T_vs_C", "factor": "group", "reference": "C", "test": "T"}


def _row(sample_id, *, role, group, subject, order, **extra):
    row = {"sample_id": sample_id, "source_file": None, "role": role, "group": group,
           "batch": "B1", "injection_order": order, "qc_pool": None, "include": True,
           "subject_id": subject, "provenance": {}, "conflicts": []}
    row.update(extra)
    return row


def _ds(groups=("C", "T"), per_group=4, reps=2, n_qc=4, sex=False, seed=5) -> DatasetState:
    """群ごとに per_group 個体 × reps 注入（技術反復）＋ QC。T 系の群は f0〜f2 が 4 倍。"""
    rng = np.random.default_rng(seed)
    names, rows, order = [], [], 0
    for i in range(n_qc):
        order += 1
        names.append(f"QC{i}")
        rows.append(_row(f"id_QC{i}", role="qc", group=None, subject=None, order=order))
    for g in groups:
        for k in range(per_group):
            for r in range(reps):
                order += 1
                name = f"{g}{k}_r{r}"
                names.append(name)
                extra = {"sex": "F" if k % 2 == 0 else "M"} if sex else {}
                rows.append(_row(f"id_{name}", role="sample", group=g, subject=f"{g}{k}",
                                 order=order, **extra))
    nf = 10
    fm = np.exp2(rng.normal(20, 1, (nf, 1)) + rng.normal(0, 0.2, (nf, len(names))))
    for j, name in enumerate(names):
        if name[0] in ("T", "U"):
            fm[:3, j] *= 4
    ds = DatasetState()
    ds.feature_matrix = fm
    ds.detected_mask = np.ones_like(fm, dtype=bool)
    ds.sample_names = names
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(len(names))]
    ds.feature_ids = [f"f{k}" for k in range(nf)]
    ds.feature_metadata = {f"f{k}": {"name": f"F{k}"} for k in range(nf)}
    apply_metadata(ds, rows)
    return ds


def _run(ds, comparisons, *, factors=("group",), analyses=None, **extra):
    spec = {"schema": "analysis-spec.v2", "preprocess": OFF,
            "design": {"factors": list(factors)}, "comparisons": comparisons,
            "analyses": analyses or [{"analysis_id": "cmp", "kind": "compare"}], **extra}
    n = normalize_spec_v2(spec)
    prepared = prepare_spec_v2(n, ds)
    return n, prepared, run_spec_v2(n, ds, prepared=prepared)


def _register(ds, comparisons=(T_VS_C,), *, threshold=1.0, **kw):
    n, prepared, rr = _run(ds, list(comparisons), **kw)
    rec = stats_results.register_compare(ds, rr, log2fc_threshold=threshold, units=prepared.units)
    return rec, rr, prepared


def _code(fn, *a, **k):
    with pytest.raises(DomainError) as e:
        fn(*a, **k)
    return e.value


def _names(ds, group):
    return [n for n in ds.sample_names if n.startswith(group) and not n.startswith("QC")]


# --- 登録 ---------------------------------------------------------------------

def test_dataset_state_has_stat_matrices_slot():
    assert DatasetState().stat_matrices == {}


def test_register_compare_stores_display_names_and_deterministic_id():
    ds = _ds()
    assert ds.sample_ids[4] == "id_C0_r0" and ds.sample_names[4] == "C0_r0"   # ID と表示名がずれる
    rec, rr, prepared = _register(ds)
    rid = rec["result_id"]
    assert rid == "cmp_" + canonical_hash({"rule_hash": rr["rule_hash"], "matrix_id": rr["matrix_id"],
                                           "log2fc_threshold": 1.0})[:12]
    assert rec["kind"] == "compare" and rec["matrix_id"] == rr["matrix_id"]
    assert rec["rule_hash"] == rr["rule_hash"] and rec["spec"] == rr["spec"]
    assert rec["engine"] is rr["analyses"][0]["result"]
    assert rec["engine"]["schema"] == "stats-compare.v1" and rec["log2fc_threshold"] == 1.0
    assert ds.results[rid] is rec and ds.last_differential is rec
    unit = rec["units"]["C0"]
    assert unit["samples"] == ["C0_r0", "C0_r1"]          # StudyData の ID ではなく表示名（R2）
    assert unit["group"] == "C" and unit["subject_id"] == "C0" and unit["cell"] == "C"
    assert set(rec["units"]) == {u.unit_id for u in prepared.units}
    assert rec["provenance"] == {"result_id": rid, "kind": "compare",
                                 "dataset_hash": dataset_fingerprint(ds),
                                 "metadata_revision": ds.metadata_revision,
                                 "parent_ids": [rr["matrix_id"]]}
    # 同じものの登録は同じ鍵を上書きする。判定に効かない閾値でも ID は分ける
    again, _, _ = _register(ds)
    assert again["result_id"] == rid and list(ds.results) == [rid]
    other, _, _ = _register(ds, threshold=2.0)
    assert other["result_id"] != rid and len(ds.results) == 2 and ds.last_differential is other
    # 記録は JSON にできる（PreparedMatrix の実体を抱えない）
    json.dumps(json_safe(rec), allow_nan=False)


def test_register_compare_requires_exactly_one_compare_analysis():
    ds = _ds()
    n, prepared, rr = _run(ds, [T_VS_C], analyses=[
        {"analysis_id": "a", "kind": "compare"}, {"analysis_id": "b", "kind": "compare"}],
        primary="a")
    with pytest.raises(ValueError):
        stats_results.register_compare(ds, rr, log2fc_threshold=1.0, units=prepared.units)
    with pytest.raises(ValueError):
        stats_results.register_compare(ds, {"status": "needs_input", "field": "x"},
                                       log2fc_threshold=1.0, units=prepared.units)
    assert ds.results == {} and ds.last_differential is None


def test_register_matrix_and_get_matrix_go_stale():
    ds = _ds()
    n = normalize_spec_v2({"schema": "analysis-spec.v2", "preprocess": OFF,
                           "comparisons": [T_VS_C],
                           "analyses": [{"analysis_id": "cmp", "kind": "compare"}]})
    prepared = prepare_spec_v2(n, ds)
    entry = stats_results.register_matrix(ds, prepared, n["preprocess"], {})
    mid = prepared.matrix_id
    assert ds.stat_matrices[mid] is entry and entry["prepared"] is prepared
    assert entry["matrix_id"] == mid and entry["config"] == n["preprocess"]
    assert entry["group_overrides"] == {} and entry["exclusions"] == {"samples": [], "features": []}
    assert entry["dataset_hash"] == dataset_fingerprint(ds)
    assert entry["metadata_revision"] == ds.metadata_revision
    assert stats_results.get_matrix(ds, mid) is entry
    assert stats_results.get_matrix(ds, "mx_nope") is None
    # 除外（metadata_revision が進む）の後は古い。除外は登録時の ds から読む
    apply_exclusions(ds, samples=["C0_r0"], features=["f9"], mode="add", reason="外れ値")
    assert ds.stat_matrices == {}                      # 前処理の入力の変更で空になる
    ds.stat_matrices[mid] = entry
    assert stats_results.get_matrix(ds, mid) is None
    study = stats_data.from_dataset(ds, exclude_samples=sorted(ds.excluded_samples),
                                    exclude_features=sorted(ds.excluded_features))
    from metabolomix.analysis.stats.preprocess import preprocess
    p2 = preprocess(study, n["preprocess"])
    e2 = stats_results.register_matrix(ds, p2, n["preprocess"], {"C1_r0": "T"})
    assert e2["exclusions"] == {"samples": ["C0_r0"], "features": ["f9"]}
    assert e2["group_overrides"] == {"C1_r0": "T"}
    assert mid not in ds.stat_matrices                 # 古い登録は新しい登録のときに捨てる
    assert stats_results.get_matrix(ds, p2.matrix_id) is e2
    e3 = stats_results.register_matrix(ds, p2, n["preprocess"], None,
                                       exclusions={"samples": ["x"], "features": []})
    assert e3["exclusions"] == {"samples": ["x"], "features": []} and e3["group_overrides"] == {}
    # データが変われば（指紋）古い
    ds.feature_matrix = ds.feature_matrix * 2
    assert stats_results.get_matrix(ds, p2.matrix_id) is None
    # 別のデータから作った行列は登録させない
    with pytest.raises(ValueError):
        stats_results.register_matrix(ds, p2, n["preprocess"], {})


def test_register_pca_record_and_staleness():
    ds = _ds()
    n = normalize_spec_v2({"schema": "analysis-spec.v2", "preprocess": OFF,
                           "comparisons": [T_VS_C],
                           "analyses": [{"analysis_id": "cmp", "kind": "compare"}]})
    prepared = prepare_spec_v2(n, ds)
    res = pca(prepared, scaling="autoscale", n_components=3, impute="half_min")
    groups = {u.unit_id: u.group for u in prepared.units}
    rec = stats_results.register_pca(ds, res, groups)
    rid = rec["result_id"]
    assert rid == "pca_" + canonical_hash({"matrix_id": prepared.matrix_id, "scaling": "autoscale",
                                           "n_components": 3, "impute": "half_min",
                                           "groups": groups})[:12]
    assert rec["kind"] == "stat_pca" and rec["engine"] is res and rec["groups"] == groups
    assert ds.last_pca is rec and ds.results[rid] is rec
    assert rec["provenance"] == {"result_id": rid, "kind": "stat_pca",
                                 "dataset_hash": dataset_fingerprint(ds),
                                 "metadata_revision": ds.metadata_revision,
                                 "parent_ids": [prepared.matrix_id]}
    assert stats_results.register_pca(ds, res, groups)["result_id"] == rid
    plain = stats_results.register_pca(ds, res, None)
    assert plain["result_id"] != rid and plain["groups"] is None
    assert result_state.is_current(ds, rec)
    ds.metadata_revision += 1
    assert not result_state.is_current(ds, rec)
    json.dumps(json_safe(rec), allow_nan=False)


# --- 2 群への読み替え -----------------------------------------------------------

def test_as_two_group_single_contrast_matches_rows():
    ds = _ds()
    rec, rr, prepared = _register(ds)
    engine = rr["analyses"][0]["result"]
    view = stats_results.as_two_group(rec)
    assert view["kind"] == "two_group" and view["a"] == "C" and view["b"] == "T"
    assert view["contrast_id"] == "T_vs_C" and view["result_id"] == rec["result_id"]
    assert view["matrix_id"] == rec["matrix_id"] and view["rule_hash"] == rec["rule_hash"]
    assert view["samples_a"] == _names(ds, "C") and view["samples_b"] == _names(ds, "T")
    assert view["n_a"] == 8 and view["n_b"] == 8                    # 注入の数
    assert view["n_units_a"] == 4 == engine["group_n"]["C"]          # 個体の数
    assert view["n_units_b"] == 4 == engine["group_n"]["T"]
    assert view["q_threshold"] == 0.05 and view["log2fc_threshold"] == 1.0
    assert view["log_transform"] is True
    assert view["contract_version"] == export_contract.CONTRACT_VERSION
    assert view["log2fc_sign"] == export_contract.LOG2FC_SIGN
    assert view["provenance"] is rec["provenance"] and view["provenance"]["kind"] == "compare"
    assert view["effect_size_definition"] == engine["effect_size_definition"]
    assert view["statistics"] == {
        "test": "moderated_t", "ci_level": 0.95, "verdict_rule": engine["verdict_rule"],
        "bh_population": {"name": "matrix", "n": engine["bh_population"]["n_by_contrast"]["T_vs_C"],
                          "q_interpretation": "discovery_fdr"},
        "subject_source": "subject_id", "n_units_a": 4, "n_units_b": 4}
    rows = [r for r in engine["rows"] if r["contrast_id"] == "T_vs_C"]
    assert len(view["results"]) == len(rows) == 10
    for got, row in zip(view["results"], rows):
        assert set(got) == {"feature", "mean_a", "mean_b", "log2fc", "t", "p", "q",
                            "ci_low", "ci_high", "verdict"}
        assert got["feature"] == row["feature_id"]
        assert got["mean_a"] == row["mean_reference"] and got["mean_b"] == row["mean_test"]
        assert got["t"] == row["statistic"]
        for key in ("log2fc", "p", "q", "ci_low", "ci_high", "verdict"):
            assert got[key] == row[key]
    up = {r["feature"] for r in view["results"] if r["verdict"] == "changed"}
    assert {"f0", "f1", "f2"} <= up
    finite = [r for r in view["results"] if r["p"] is not None and math.isfinite(r["p"])]
    assert len(view["volcano"]) == len(finite)
    point = view["volcano"][0]
    assert set(point) == {"feature", "log2fc", "neg_log10_p", "sig"}
    sig = {p["feature"]: p["sig"] for p in view["volcano"]}
    assert sig["f0"] == "up" and sig["f9"] == "ns"
    assert view["caveats"] == [w["message"] for w in engine["warnings"] if w.get("message")]
    # 読み替えは記録を書き換えない
    view["results"][0]["log2fc"] = 999
    assert rec["engine"]["rows"][0]["log2fc"] != 999
    json.loads(json_payload(json_safe(view)))


def test_as_two_group_requires_contrast_for_multi():
    ds = _ds(groups=("C", "T", "U"))
    comps = [T_VS_C, {"id": "U_vs_C", "factor": "group", "reference": "C", "test": "U"}]
    rec, rr, _ = _register(ds, comps)
    assert stats_results.contrast_ids(rec) == ["T_vs_C", "U_vs_C"]
    e = _code(stats_results.as_two_group, rec)
    assert e.code == "CONTRAST_REQUIRED" and e.details["contrast_ids"] == ["T_vs_C", "U_vs_C"]
    e = _code(stats_results.as_two_group, rec, "nope")
    assert e.code == "CONTRAST_REQUIRED" and e.details["contrast_ids"] == ["T_vs_C", "U_vs_C"]
    view = stats_results.as_two_group(rec, "U_vs_C")
    assert view["a"] == "C" and view["b"] == "U" and view["contrast_id"] == "U_vs_C"
    assert view["samples_b"] == _names(ds, "U") and view["n_units_b"] == 4
    assert [r["feature"] for r in view["results"]] == [f"f{k}" for k in range(10)]
    engine_u = [r for r in rr["analyses"][0]["result"]["rows"] if r["contrast_id"] == "U_vs_C"]
    assert [r["p"] for r in view["results"]] == [r["p"] for r in engine_u]
    assert view["statistics"]["bh_population"]["n"] == (
        rr["analyses"][0]["result"]["bh_population"]["n_by_contrast"]["U_vs_C"])


def test_as_two_group_rejects_interaction():
    ds = _ds(sex=True)
    at_f = {"id": "T_vs_C_F", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "F"}}
    inter = {"id": "ix", "kind": "interaction", "factors": ["group", "sex"],
             "levels": {"group": ["C", "T"], "sex": ["F", "M"]}}
    rec, rr, _ = _register(ds, [at_f, inter], factors=("group", "sex"))
    assert stats_results.contrast_ids(rec) == ["T_vs_C_F", "ix"]
    assert stats_results.contrast_ids(rec, two_group_only=True) == ["T_vs_C_F"]
    e = _code(stats_results.as_two_group, rec, "ix")
    assert e.code == "CONTRAST_NOT_TWO_GROUP"
    assert e.details["contrast_id"] == "ix" and e.details["two_group_contrast_ids"] == ["T_vs_C_F"]
    # at 付きの対比は、その水準の個体だけを a / b に数える
    view = stats_results.as_two_group(rec, "T_vs_C_F")
    assert view["samples_a"] == ["C0_r0", "C0_r1", "C2_r0", "C2_r1"]
    assert view["samples_b"] == ["T0_r0", "T0_r1", "T2_r0", "T2_r1"]
    assert view["n_units_a"] == 2 == rr["analyses"][0]["result"]["group_n"]["C:F"]
    assert view["n_a"] == 4
    assert rec["units"]["C1"]["cell"] == "C:M"
    # 対比が 1 つだけの交互作用も読み替えない
    rec2, _, _ = _register(ds, [inter], factors=("group", "sex"))
    assert _code(stats_results.as_two_group, rec2).code == "CONTRAST_NOT_TWO_GROUP"


def test_as_two_group_passes_v1_two_group_through():
    v1 = {"kind": "two_group", "a": "C", "b": "T", "results": [], "volcano": [],
          "provenance": {"result_id": "res_x", "kind": "differential"}}
    assert stats_results.as_two_group(v1) is v1
    assert stats_results.contrast_ids(v1) == []
    assert _code(stats_results.as_two_group, {"kind": "stat_pca"}).code == "ANALYSIS_RESULT_NOT_FOUND"


# --- 古さ ---------------------------------------------------------------------

def test_compare_record_stale_after_exclusion_and_manifest_change():
    # 除外
    ds = _ds()
    rec, _, _ = _register(ds)
    assert stats_results.is_stats_current(ds, rec) and result_state.is_current(ds, rec)
    result_state.assert_current(ds, rec)
    apply_exclusions(ds, samples=["C0_r0"], features=None, mode="add", reason="外れ値")
    assert not stats_results.is_stats_current(ds, rec)
    assert not result_state.is_current(ds, rec)
    e = _code(result_state.assert_current, ds, rec)
    assert e.code == "STALE_ANALYSIS_RESULT"
    assert e.details["result_id"] == rec["result_id"] and e.details["kind"] == "compare"
    assert e.details["result_metadata_revision"] == rec["provenance"]["metadata_revision"]
    assert e.details["metadata_revision"] == ds.metadata_revision

    # 試料の表（群だけの訂正でも古い）
    ds = _ds()
    rec, _, _ = _register(ds)
    rows = copy.deepcopy(ds.sample_metadata_rows)
    apply_metadata(ds, copy.deepcopy(rows))            # 何も変わらない再適用は古くしない
    assert result_state.is_current(ds, rec)
    rows[4]["group"] = "T"
    apply_metadata(ds, rows)
    assert ds.last_differential is None                # 群の変更で直近の比較は外れる
    assert rec["result_id"] in ds.results and not result_state.is_current(ds, rec)
    assert _code(result_state.assert_current, ds, rec).code == "STALE_ANALYSIS_RESULT"

    # データ（指紋）の変化
    ds = _ds()
    rec, _, _ = _register(ds)
    h = dataset_fingerprint(ds)
    assert stats_results.is_stats_current(ds, rec, dataset_hash=h)
    ds.feature_matrix = ds.feature_matrix * 2
    assert not stats_results.is_stats_current(ds, rec)
    assert not stats_results.is_stats_current(ds, {"kind": "compare"})   # 来歴が無ければ確かめられない


def test_invalidate_on_exclusion_clears_stat_and_analysis_matrices():
    def fresh():
        ds = DatasetState()
        ds.stat_matrices = {"mx_a": {"matrix_id": "mx_a"}}
        ds.analysis_matrices = {"m": {}}
        return ds

    for changed in ({"exclusion"}, {"role"}, {"batch"}, {"detection"}, {"dataset"},
                    {"matrix_recipe"}):
        ds = fresh()
        invalidate_results(ds, changed)
        assert ds.stat_matrices == {} and ds.analysis_matrices == {}, changed
    # 群・比較・PCA の設定・v1 の前処理のやり直し（recipe）は行列の入力ではない
    for changed in ({"group"}, {"comparison"}, {"pca_settings"}, {"recipe"}, {"note"}):
        ds = fresh()
        invalidate_results(ds, changed)
        assert ds.stat_matrices and ds.analysis_matrices, changed
    # 実際の除外の経路
    ds = _ds()
    n = normalize_spec_v2({"schema": "analysis-spec.v2", "preprocess": OFF,
                           "comparisons": [T_VS_C],
                           "analyses": [{"analysis_id": "cmp", "kind": "compare"}]})
    stats_results.register_matrix(ds, prepare_spec_v2(n, ds), n["preprocess"], {})
    ds.analysis_matrices = {"m": {}}
    assert ds.stat_matrices
    apply_exclusions(ds, samples=None, features=["f1"], mode="add", reason="割れ")
    assert ds.stat_matrices == {} and ds.analysis_matrices == {}


def test_derive_resets_stat_matrices():
    from metabolomix.quant.derive import derive, quant_state

    ds = _ds()
    ds.stat_matrices = {"mx": {"matrix_id": "mx"}}
    new = derive(ds, matrix=ds.feature_matrix * 2, detected=ds.detected_mask.copy(),
                 state=quant_state(normalized=False, method=None, unit=None, record_id="q_1",
                                   steps=[], skipped=[]))
    assert new.stat_matrices == {} and ds.stat_matrices == {"mx": {"matrix_id": "mx"}}


# --- JSON -------------------------------------------------------------------

def test_json_safe_makes_full_engine_result_strict():
    ds = _ds()
    _, _, rr = _run(ds, [T_VS_C], analyses=[{"analysis_id": "cmp", "kind": "compare"},
                                            {"analysis_id": "p", "kind": "pca"}])
    prior = rr["analyses"][0]["result"]["prior"]
    assert prior["df_prior"] == math.inf               # 分散が揃うと事前分布の自由度は無限大
    with pytest.raises(ValueError):
        json.dumps(rr, allow_nan=False)
    safe = json_safe(rr)
    assert safe["analyses"][0]["result"]["prior"]["df_prior"] == "inf"
    text = json.dumps(safe, allow_nan=False)
    json.loads(text, parse_constant=lambda name: (_ for _ in ()).throw(ValueError(name)))
