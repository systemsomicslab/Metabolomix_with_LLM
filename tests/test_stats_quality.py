"""品質の純計算（metabolomix.analysis.stats.quality。概観の品質の段）。"""
from __future__ import annotations

import numpy as np

from metabolomix.analysis.stats import data as stats_data
from metabolomix.analysis.stats import quality
from tests.test_stats_results import _ds


def _study(ds):
    return stats_data.from_dataset(ds)


def test_qc_rsd_summary_matches_definition():
    ds = _ds(n_qc=4)
    study = _study(ds)
    qc = [i for i, s in enumerate(study.samples) if s.role == "qc"]
    block = study.values[qc]
    rsd = block.std(axis=0, ddof=1) / block.mean(axis=0)
    out = quality.qc_rsd_summary(study)
    assert out["status"] == "ok"
    assert out["n_qc"] == 4
    assert out["n_computable"] == study.values.shape[1]
    assert np.isclose(out["median"], float(np.median(rsd)))
    assert np.isclose(out["frac_le_030"], float((rsd <= 0.30).mean()))


def test_qc_rsd_not_evaluable_with_one_qc():
    out = quality.qc_rsd_summary(_study(_ds(n_qc=1)))
    assert out["status"] == "not_evaluable"
    assert out["reason"] == "too_few_qc"


def test_failed_qc_flags_low_total():
    ds = _ds(n_qc=4)
    ds.feature_matrix[:, ds.sample_names.index("QC2")] *= 0.1
    out = quality.failed_qc(_study(ds))
    assert out["checked"] == 4
    assert out["failed"] == ["id_QC2"]


def test_failed_qc_needs_three_qc():
    assert quality.failed_qc(_study(_ds(n_qc=2))) == {"checked": 2, "failed": []}


def test_detection_summary_without_mask():
    ds = _ds()
    ds.detected_mask = None
    assert quality.detection_summary(_study(ds))["status"] == "not_evaluable"


def test_detection_summary_counts_gap_fill():
    ds = _ds(n_qc=2)
    bio = [j for j, n in enumerate(ds.sample_names) if not n.startswith("QC")]
    ds.detected_mask[0, bio[0]] = False          # gap-fill（値はあるが未検出）
    ds.feature_matrix[1, bio[1]] = 0.0           # 欠測
    out = quality.detection_summary(_study(ds))
    n_cells = len(bio) * ds.feature_matrix.shape[0]
    assert out["status"] == "ok"
    assert np.isclose(out["gap_filled_frac"], 1 / n_cells)
    assert np.isclose(out["missing_frac"], 1 / n_cells)
    assert np.isclose(out["detected_frac"], (n_cells - 2) / n_cells)


def test_drift_signs_detects_monotone_qc():
    ds = _ds(n_qc=6)
    for k in range(6):
        ds.feature_matrix[:, ds.sample_names.index(f"QC{k}")] *= 2.0 ** k
    signs = quality.drift_signs(_study(ds))
    assert len(signs) == 1 and signs[0]["batch"] == "B1" and signs[0]["n_qc"] == 6
    assert signs[0]["flag"] is True


def test_drift_signs_random_qc_not_flagged_and_skips_small_batches():
    assert quality.drift_signs(_study(_ds(n_qc=6, seed=11)))[0]["flag"] is False
    assert quality.drift_signs(_study(_ds(n_qc=4))) == []


def _conf(status="none", inj="ok"):
    return {"batch_status": status, "injection_order": {"status": inj}, "cramers_v": None}


def test_quality_verdict_rules():
    ok_rsd = {"status": "ok", "median": 0.1, "n_qc": 4}
    none_failed = {"checked": 4, "failed": []}
    verdict, reasons = quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[],
                                               confounding=_conf("complete"), replicate_warnings=[])
    assert verdict == "poor" and "complete_confounding" in reasons
    verdict, reasons = quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[],
                                               confounding=_conf("partial"), replicate_warnings=[])
    assert verdict == "caution" and reasons == ["partial_confounding"]
    assert quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[], confounding=_conf(),
                                   replicate_warnings=[]) == ("ok", [])
    bad_rsd = {"status": "ok", "median": 0.45, "n_qc": 4}
    assert quality.quality_verdict(rsd=bad_rsd, failed=none_failed, drift=[], confounding=_conf(),
                                   replicate_warnings=[])[0] == "poor"
    assert quality.quality_verdict(rsd=ok_rsd, failed={"checked": 4, "failed": ["q"]}, drift=[],
                                   confounding=_conf(), replicate_warnings=[])[0] == "poor"
    assert quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[{"flag": True}],
                                   confounding=_conf(), replicate_warnings=[])[0] == "caution"
    no_qc = {"status": "not_evaluable", "reason": "too_few_qc", "n_qc": 0}
    verdict, _ = quality.quality_verdict(rsd=no_qc, failed={"checked": 0, "failed": []}, drift=[],
                                         confounding=_conf("not_evaluable", "not_evaluable"),
                                         replicate_warnings=[])
    assert verdict == "not_evaluable"
    # QC だけが無いときも ok とは言わない（spec §4.1: 評価できないことをそのまま書く）
    assert quality.quality_verdict(rsd=no_qc, failed={"checked": 0, "failed": []}, drift=[],
                                   confounding=_conf(), replicate_warnings=[]) ==         ("not_evaluable", ["qc_not_evaluable"])
    # 他の項目で caution / poor なら、その判定に QC を評価できなかった理由を添える
    verdict, reasons = quality.quality_verdict(rsd=no_qc, failed={"checked": 0, "failed": []},
                                               drift=[{"flag": True}], confounding=_conf(),
                                               replicate_warnings=[])
    assert verdict == "caution" and reasons == ["drift_sign", "qc_not_evaluable"]
