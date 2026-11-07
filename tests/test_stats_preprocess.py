"""stats.preprocess: 選択・バッチ別ドリフト補正・同時の足切り・正規化・log2・個体でのまとめ。"""
import copy

import numpy as np
import pytest

from metabolomix.analysis.stats import preprocess as P
from metabolomix.analysis.stats.data import QuantState, SampleInfo, StudyData
from metabolomix.core.atomic_io import DomainError

OFF_FILTERS = {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None}


def _study(values, rows, *, detected=None, quant=None, feature_ids=None, dataset_hash="h"):
    """rows: (sample_id, role, group, batch, order, subject) のタプル列。"""
    values = np.asarray(values, dtype=float)
    samples = [SampleInfo(sample_id=r[0], role=r[1], group=r[2], batch=r[3],
                          injection_order=r[4], subject_id=r[5],
                          include=(r[6] if len(r) > 6 else True),
                          attributes=(r[7] if len(r) > 7 else {}))
               for r in rows]
    fids = feature_ids or [f"f{j}" for j in range(values.shape[1])]
    return StudyData(values=values, detected=None if detected is None else np.asarray(detected, bool),
                     sample_ids=[s.sample_id for s in samples], feature_ids=fids,
                     samples=samples, quant=quant, dataset_hash=dataset_hash)


def _simple(n_feat=3, *, with_qc=True, quant=None, detected=None):
    rng = np.random.default_rng(0)
    rows = [(f"a{i}", "sample", "G1", None, None, None) for i in range(3)]
    rows += [(f"b{i}", "sample", "G2", None, None, None) for i in range(3)]
    if with_qc:
        rows += [(f"q{i}", "qc", None, None, None, None) for i in range(3)]
    values = 1000.0 * (1 + 0.05 * rng.standard_normal((len(rows), n_feat)))
    return _study(values, rows, quant=quant, detected=detected)


def _step(prep, name):
    return next(s for s in prep.steps if s["step"] == name)


# --- 設定 -------------------------------------------------------------------

def test_config_defaults_merge_and_reject_unknown():
    cfg = P.normalize_preprocess_config(None)
    assert cfg == P.PREPROCESS_DEFAULTS and cfg is not P.PREPROCESS_DEFAULTS
    cfg = P.normalize_preprocess_config({"drift": {"span": 0.5}, "filters": {"qc_rsd": None}})
    assert cfg["drift"] == {"mode": "auto", "span": 0.5, "min_qc_per_batch": 4}
    assert cfg["filters"]["qc_rsd"] is None and cfg["filters"]["d_ratio"] == {"max": 0.5}
    for bad in ({"nope": 1}, {"drift": {"nope": 1}}, {"filters": {"qc_rsd": {"min": 1}}},
                {"drift": {"mode": "maybe"}}, {"normalize": {"method": "zscore"}},
                {"filters": {"qc_rsd": {"max": "high"}}}, {"drift": {"span": True}},
                {"collapse_replicates": "yes"}, {"filters": {"detection": {"by": "batch"}}}):
        with pytest.raises(DomainError) as e:
            P.normalize_preprocess_config(bad)
        assert e.value.code == "PREPROCESS_CONFIG_INVALID"


# --- 段の並び ---------------------------------------------------------------

def test_step_order_is_select_drift_filter_normalize_log_collapse():
    prep = P.preprocess(_simple(), {})
    assert [s["step"] for s in prep.steps] == ["select", "drift", "filter", "normalize", "log2", "collapse"]
    assert all(s["status"] in ("applied", "skipped") for s in prep.steps)
    assert all("reason" in s for s in prep.steps if s["status"] == "skipped")


def test_select_drops_excluded_and_other_roles():
    rows = [("a", "sample", "G1", None, None, None), ("b", "sample", "G1", None, None, None, False),
            ("c", "sample", "G2", None, None, None), ("d", "sample", "G2", None, None, None),
            ("s", "standard", None, None, None, None), ("e", "sample", "G1", None, None, None)]
    prep = P.preprocess(_study(np.full((6, 2), 10.0), rows), {"filters": OFF_FILTERS})
    sel = _step(prep, "select")
    assert [u.unit_id for u in prep.units] == ["a", "c", "d", "e"]
    assert sel["excluded_by_include"] == ["b"] and sel["excluded_by_role"] == {"standard": ["s"]}


# --- ドリフト ---------------------------------------------------------------

def _drift_study(*, b2_not_interspersed=False, quant=None, qc_step=3):
    """2 バッチ × 16 注入。QC は qc_step ごと。B1 は上向き・B2 は下向きの線形ドリフトで水準も違う。"""
    rows, vals = [], []
    for batch, start, base, slope in (("B1", 1, 1000.0, 0.03), ("B2", 17, 3000.0, -0.02)):
        for o in range(start, start + 16):
            is_qc = (o - start) % qc_step == 0 and not (batch == "B2" and b2_not_interspersed)
            group = None if is_qc else ("G1" if o % 2 else "G2")
            rows.append((f"{batch}_{o}", "qc" if is_qc else "sample", group, batch, float(o), None))
            vals.append([base * (1 + slope * (o - start)), 2 * base * (1 + slope * (o - start))])
    if b2_not_interspersed:
        # QC を挟まない配置: B2 の QC を試料より前の注入順に置く
        for k in range(6):
            rows.append((f"B2_0.{k + 1}", "qc", None, "B2", 0.1 * (k + 1), None))
            vals.append([3000.0, 6000.0])
    return _study(vals, rows, quant=quant)


def test_drift_is_corrected_per_batch_then_levels_aligned():
    prep = P.preprocess(_drift_study(), {"filters": OFF_FILTERS, "normalize": {"method": "none"}})
    drift = _step(prep, "drift")
    assert drift["status"] == "applied" and drift["batches"] == ["B1", "B2"]
    qc_batch = [i.split("_")[0] for i in prep.qc_ids]
    orders = np.array([float(i.split("_")[1]) for i in prep.qc_ids])
    medians = []
    for b in ("B1", "B2"):
        rows = [k for k, x in enumerate(qc_batch) if x == b]
        for j in range(2):
            y = prep.qc_log2[rows, j]
            assert abs(np.polyfit(orders[rows], y, 1)[0]) < 1e-9
        medians.append(np.median(prep.qc_log2[rows], axis=0))
    assert np.allclose(medians[0], medians[1])
    # 生体試料も同じ水準へ（線形のドリフトは局所線形の LOESS で消える）
    assert np.nanstd(prep.log2[:, 0]) < 1e-9


def test_drift_skipped_when_qc_not_interspersed_in_a_batch():
    study = _drift_study(b2_not_interspersed=True)
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    drift = _step(prep, "drift")
    assert drift["status"] == "skipped" and drift["reason"] == "qc_not_interspersed:B2"
    with pytest.raises(DomainError) as e:
        P.preprocess(study, {"filters": OFF_FILTERS, "drift": {"mode": "on"}})
    assert e.value.code == "PREPROCESS_PRECONDITION"


def test_drift_skipped_reasons_no_order_and_too_few_qc():
    prep = P.preprocess(_simple(), {})
    assert _step(prep, "drift")["reason"] == "no_injection_order"
    prep = P.preprocess(_drift_study(), {"filters": OFF_FILTERS, "drift": {"min_qc_per_batch": 7}})
    assert _step(prep, "drift")["reasons"] == ["too_few_qc:B1", "too_few_qc:B2"]


def test_evenly_spaced_four_qc_fit_is_recorded_not_hidden():
    """statsmodels の xvals 経路は等間隔の 4 QC（span 0.75 → 近傍 3）で内側の QC に NaN を返す。
    その特徴はそのバッチで補正せず、数をバッチ別に記録する（補間で埋めない）。"""
    prep = P.preprocess(_drift_study(qc_step=5), {"filters": OFF_FILTERS, "normalize": {"method": "none"}})
    drift = _step(prep, "drift")
    assert drift["status"] == "applied"
    assert drift["uncorrected_fit_invalid"] == {"B1": 2, "B2": 2}


def test_is_normalized_with_qc_and_auto_drift_returns_needs_input():
    quant = QuantState(normalized=True, method="internal_standard", unit="pmol")
    out = P.preprocess(_drift_study(quant=quant), {})
    assert isinstance(out, P.NeedsInput)
    assert out.field == "preprocess.drift.mode"
    assert [o["value"] for o in out.options] == ["on", "off"]
    env = out.to_envelope()
    assert env["status"] == "needs_input" and env["field"] == out.field and env["reason"]
    # 明示すれば進む
    prep = P.preprocess(_drift_study(quant=quant), {"drift": {"mode": "on"}, "filters": OFF_FILTERS})
    assert _step(prep, "drift")["status"] == "applied"


# --- 正規化 -----------------------------------------------------------------

def test_is_normalized_defaults_to_no_normalization():
    quant = QuantState(normalized=True, method="internal_standard", unit="pmol")
    prep = P.preprocess(_simple(quant=quant), {"drift": {"mode": "off"}})
    norm = _step(prep, "normalize")
    assert norm["status"] == "skipped" and norm["method"] == "none" and norm["reason"] == "quant_normalized"


def test_default_normalization_is_pqn_with_qc_reference():
    prep = P.preprocess(_simple(), {})
    norm = _step(prep, "normalize")
    assert norm["status"] == "applied" and norm["requested"] == "auto"
    assert norm["method"] == "pqn" and norm["reference"] == "qc"


def test_normalization_factors_use_only_features_kept_by_filters():
    base = _simple(n_feat=4)
    rows = [(s.sample_id, s.role, s.group, s.batch, s.injection_order, s.subject_id) for s in base.samples]
    outlier = np.r_[1e6 * np.arange(1, 7), 1.0, 1e7, 50.0]
    with_outlier = _study(np.column_stack([base.values, outlier]), rows)
    cfg = {"normalize": {"method": "tic"}}
    a = _step(P.preprocess(base, cfg), "normalize")
    b_prep = P.preprocess(with_outlier, cfg)
    b = _step(b_prep, "normalize")
    assert "f4" in _step(b_prep, "filter")["filters"]["qc_rsd"]["dropped_ids"]
    assert a["factors"] == pytest.approx(b["factors"])
    # 外れ特徴を残すと係数は変わる（この検査が意味を持つことの確認）
    c = _step(P.preprocess(with_outlier, {**cfg, "filters": OFF_FILTERS}), "normalize")
    assert c["factors"] != pytest.approx(a["factors"])


def test_no_qc_skips_qc_steps_with_reasons():
    prep = P.preprocess(_simple(with_qc=False), {})
    assert _step(prep, "drift")["reason"] == "no_qc"
    filt = _step(prep, "filter")["filters"]
    assert filt["qc_rsd"]["status"] == "skipped" and filt["qc_rsd"]["reason"] == "too_few_qc"
    assert filt["d_ratio"]["status"] == "skipped" and filt["d_ratio"]["reason"] == "too_few_qc"
    norm = _step(prep, "normalize")
    assert norm["method"] == "pqn" and norm["reference"] == "sample"
    assert prep.qc_log2.shape == (0, 3) and prep.qc_ids == []


# --- 足切り -----------------------------------------------------------------

def test_filters_masks_computed_after_drift_and_applied_together():
    study = _drift_study()
    rows = [(s.sample_id, s.role, s.group, s.batch, s.injection_order, s.subject_id) for s in study.samples]
    rows.append(("blank1", "blank", None, "B1", 5.5, None))
    order = np.array([r[4] for r in rows])
    # f2: 強いドリフトで生の QC の RSD は大きいが、補正後は小さい → 残る
    start = np.where(order >= 17, 17.0, 1.0)
    strong = 1000.0 * (1 + 0.2 * (order - start))
    strong[-1] = 10.0
    # f3: ブランクで大きく QC がばらつく → ブランク比と QC RSD の両方で落ちる
    is_qc = np.array([r[1] == "qc" for r in rows])
    noisy = np.where(is_qc, np.where(np.cumsum(is_qc) % 2, 100.0, 900.0), 100.0)
    noisy[-1] = 1000.0
    values = np.column_stack([np.vstack([study.values, [10.0, 10.0]]), strong, noisy])
    s2 = _study(values, rows)
    raw_qc = values[is_qc, 2]
    assert np.std(raw_qc, ddof=1) / np.mean(raw_qc) > 0.3
    prep = P.preprocess(s2, {"filters": {"d_ratio": None}})
    filt = _step(prep, "filter")
    assert "f2" in prep.feature_ids
    assert "f3" in filt["filters"]["qc_rsd"]["dropped_ids"]
    assert "f3" in filt["filters"]["blank_ratio"]["dropped_ids"]
    assert filt["dropped"] == 1 and prep.feature_ids == ["f0", "f1", "f2"]


def test_d_ratio_filter_default_half():
    rows = [(f"s{i}", "sample", "G1" if i < 4 else "G2", None, None, None) for i in range(8)]
    rows += [(f"q{i}", "qc", None, None, None, None) for i in range(4)]
    bio = np.array([80.0, 90, 110, 120, 80, 90, 110, 120])            # SD ≈ 17.9
    qc_lo = np.array([96.0, 98, 102, 104])                            # SD ≈ 3.65 → D ≈ 0.20
    qc_hi = np.array([86.0, 92, 108, 114])                            # SD ≈ 12.6 → D ≈ 0.70
    values = np.column_stack([np.r_[bio, qc_lo], np.r_[bio, qc_hi]])
    prep = P.preprocess(_study(values, rows), {"normalize": {"method": "none"}})
    filt = _step(prep, "filter")["filters"]
    assert P.PREPROCESS_DEFAULTS["filters"]["d_ratio"]["max"] == 0.5
    assert filt["d_ratio"]["dropped_ids"] == ["f1"] and filt["qc_rsd"]["dropped_ids"] == []
    assert prep.feature_ids == ["f0"]


def test_detection_filter_by_group_keeps_presence_only_features():
    rows = [(f"a{i}", "sample", "G1", None, None, None) for i in range(4)]
    rows += [(f"b{i}", "sample", "G2", None, None, None) for i in range(4)]
    values = np.full((8, 3), 100.0)
    det = np.ones((8, 3), bool)
    det[4:, 1] = False                       # f1: G1 だけで検出（G2 は gap-fill）→ 残す
    det[:, 2] = [True, False, False, False, True, False, False, False]  # f2: 各群 1/4 → 落とす
    prep = P.preprocess(_study(values, rows, detected=det), {})
    filt = _step(prep, "filter")["filters"]["detection"]
    assert filt["status"] == "applied" and filt["dropped_ids"] == ["f2"]
    assert prep.feature_ids == ["f0", "f1"]
    assert prep.detected_units[:, 1].tolist() == [True] * 4 + [False] * 4
    # detected が無ければ飛ばす
    prep2 = P.preprocess(_study(values, rows), {})
    assert _step(prep2, "filter")["filters"]["detection"]["reason"] == "no_detected_mask"
    assert prep2.detected_units is None


# --- log2・まとめ -----------------------------------------------------------

def test_log2_marks_non_positive_missing_without_pseudocount():
    rows = [(f"s{i}", "sample", "G1" if i < 2 else "G2", None, None, None) for i in range(4)]
    values = np.array([[8.0, 0.0], [8.0, 4.0], [16.0, np.nan], [16.0, -1.0]])
    prep = P.preprocess(_study(values, rows), {"filters": OFF_FILTERS, "normalize": {"method": "none"}})
    assert prep.log2[:, 0].tolist() == [3.0, 3.0, 4.0, 4.0]
    assert np.isnan(prep.log2[[0, 2, 3], 1]).all() and prep.log2[1, 1] == 2.0
    step = _step(prep, "log2")
    assert step["missing_cells"] == {"non_positive": 2, "non_finite": 1}


def test_replicates_averaged_on_log_scale_by_subject():
    rows = [("i1", "sample", "G1", "B", 1.0, "S1", True, {"sex": "F", "day": "1"}),
            ("i2", "sample", "G1", "B", 2.0, "S1", True, {"sex": "F", "day": "2"}),
            ("i3", "sample", "G2", "B", 3.0, "S2"), ("i4", "sample", "G2", "B", 4.0, "S3")]
    values = np.array([[4.0, 5.0], [16.0, 0.0], [2.0, 0.0], [2.0, 3.0]])
    det = np.array([[False, True], [True, False], [True, False], [True, True]])
    prep = P.preprocess(_study(values, rows, detected=det),
                        {"filters": OFF_FILTERS, "normalize": {"method": "none"},
                         "drift": {"mode": "off"}})
    assert prep.subject_source == "subject_id"
    assert [u.unit_id for u in prep.units] == ["S1", "S2", "S3"]
    s1 = prep.units[0]
    assert (s1.n_injections, s1.sample_ids, s1.injection_order, s1.group, s1.batch) == \
        (2, ["i1", "i2"], 1.5, "G1", "B")
    assert s1.attributes == {"sex": "F"}
    assert prep.log2[0, 0] == pytest.approx(3.0)
    assert prep.log2[0, 1] == pytest.approx(np.log2(5.0))     # 有限値だけの平均
    assert np.isnan(prep.log2[1, 1])                          # 全注入が欠測
    assert prep.detected_units[0].tolist() == [True, True]
    assert prep.detected_units[1].tolist() == [True, False]
    assert prep.warnings == []


def test_subject_inconsistent_and_mixed_source():
    rows = [("i1", "sample", "G1", None, None, "S1"), ("i2", "sample", "G2", None, None, "S1"),
            ("i3", "sample", "G2", None, None, None)]
    with pytest.raises(DomainError) as e:
        P.preprocess(_study(np.ones((3, 1)), rows), {"filters": OFF_FILTERS})
    assert e.value.code == "SUBJECT_INCONSISTENT"
    rows[1] = ("i2", "sample", "G1", None, None, "S1")
    prep = P.preprocess(_study(np.ones((3, 1)), rows), {"filters": OFF_FILTERS})
    assert prep.subject_source == "mixed" and [u.unit_id for u in prep.units] == ["S1", "i3"]


def test_injection_source_reports_replicate_name_warnings():
    rows = [(n, "sample", g, None, None, None) for n, g in
            (("P01_1", "G1"), ("P01_2", "G1"), ("P02_1", "G2"), ("P02_2", "G2"))]
    prep = P.preprocess(_study(np.ones((4, 1)), rows), {"filters": OFF_FILTERS})
    assert prep.subject_source == "injection" and len(prep.units) == 4
    assert {w["stem"] for w in prep.warnings} == {"P01", "P02"}


# --- matrix_id・補完 --------------------------------------------------------

def test_matrix_id_is_deterministic_and_config_sensitive():
    a = P.preprocess(_simple(), {}).matrix_id
    assert a == P.preprocess(_simple(), {}).matrix_id
    assert a.startswith("mx_") and len(a) == 19
    assert a == P.preprocess(_simple(), copy.deepcopy(P.PREPROCESS_DEFAULTS)).matrix_id
    assert a != P.preprocess(_simple(), {"normalize": {"method": "median"}}).matrix_id
    s = _simple()
    s.samples[0] = SampleInfo(**{**s.samples[0].__dict__, "group": "G9"})
    assert a != P.preprocess(s, {}).matrix_id
    nohash = _simple()
    nohash.dataset_hash = None
    other = _simple()
    other.dataset_hash = None
    other.values = other.values * 2
    assert P.preprocess(nohash, {}).matrix_id == P.preprocess(_simple_nohash(), {}).matrix_id
    assert P.preprocess(nohash, {}).matrix_id != P.preprocess(other, {}).matrix_id


def _simple_nohash():
    s = _simple()
    s.dataset_hash = None
    return s


def test_impute_for_pca_half_min_and_report():
    log2 = np.array([[1.0, np.nan], [3.0, 2.0], [np.nan, 4.0]])
    out, rep = P.impute_for_pca(log2)
    assert out.tolist() == [[1.0, 1.0], [3.0, 2.0], [0.0, 4.0]]
    assert rep == {"method": "half_min", "imputed_cells": 2}
    assert np.isnan(log2[0, 1])                               # 入力は書き換えない
    out2, _ = P.impute_for_pca(log2, reference_rows=np.array([1, 2]))
    assert out2[2, 0] == 2.0 and out2[0, 1] == 1.0
    out3, rep3 = P.impute_for_pca(np.array([[1.0, 2.0], [np.nan, 3.0], [1.5, 2.5]]), method="knn")
    assert np.isfinite(out3).all() and rep3 == {"method": "knn", "imputed_cells": 1}
