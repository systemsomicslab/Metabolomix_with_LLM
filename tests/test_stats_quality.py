"""品質の純計算（metabolomix.analysis.stats.quality。概観の品質の段）。"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

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


def test_drift_signs_mark_rows_whose_qc_order_may_be_file_order():
    """ρ に使った QC の注入順が実際の測定順とは限らない（MS-DIAL の既定〔ファイルの順〕かもしれない・試料の表が
    宣言でファイル名の順と書いた）なら、その行に印を付ける。ρ・p・判定は変えない。確かな順の行には鍵ごと付けない。"""
    from scipy import stats as sps

    study = _study(_ds(n_qc=6, seed=11))
    plain = quality.drift_signs(study)
    assert len(plain) == 1 and "injection_order_may_be_file_order" not in plain[0]
    qc = [i for i, s in enumerate(study.samples) if s.role == "qc"]
    rho, p = sps.spearmanr([study.samples[i].injection_order for i in qc],
                           np.log2(study.values[qc].sum(axis=1)))
    for marked in ({"injection_order_may_be_file_order": True}, {"injection_order_source": "filename"}):
        unsure = replace(study, samples=[replace(s, **marked) if s.role == "qc" else s for s in study.samples])
        rows = quality.drift_signs(unsure)
        assert rows[0]["injection_order_may_be_file_order"] is True, marked
        assert np.isclose(rows[0]["rho"], rho) and np.isclose(rows[0]["p"], p)
        assert rows[0]["flag"] == plain[0]["flag"]
    # 生体試料の注入順だけが不確かなら、QC の注入順で測る ρ には関わらないので付けない
    bio = replace(study, samples=[replace(s, injection_order_source="filename") if s.role == "sample" else s
                                  for s in study.samples])
    assert "injection_order_may_be_file_order" not in quality.drift_signs(bio)[0]


# --- 技術反復らしいまとまり（replicate_blocks） ------------------------------------------------------

def _block_study(*, groups=("A", "B"), animals=3, reps=3, replicate=True, layout="consecutive", subject=False,
                 drift=0.0, n_features=300, seed=0):
    """群ごとに animals 匹 × reps 注入（replicate=False なら animals × reps 匹の独立な試料）。試料名に反復の印は無い。

    `layout`: 群の中の注入順。consecutive（同じ個体の注入が続く）/ interleaved（1, 2, 3, 1, 2, 3, …）/
    random（群の中で無作為）。`drift` は注入順に比例して特徴ごとに違う向きに動くドリフト（log2 の傾きの SD）。
    戻り値は StudyData と、期待するまとまり（群 → 試料 ID の組の集合）。
    """
    from metabolomix.analysis.stats.data import SampleInfo, StudyData

    rng = np.random.default_rng(seed)
    mu = rng.normal(18, 2, n_features)
    slope = rng.normal(0, drift, n_features)
    rows, samples, expected, order = [], [], {}, 0
    n_group = animals * reps
    for g in groups:
        if replicate:
            animal = rng.normal(0, 0.4, (animals, n_features))
            subjects = [a for a in range(animals) for _ in range(reps)]
            profiles = [animal[a] + rng.normal(0, 0.12, n_features) for a in subjects]
        else:
            subjects = list(range(n_group))
            profiles = [rng.normal(0, 0.4, n_features) + rng.normal(0, 0.12, n_features) for _ in subjects]
        if layout == "consecutive":
            run = list(range(n_group))
        elif layout == "interleaved":
            run = [r * animals + a for a in range(animals) for r in range(reps)]
        else:
            run = list(rng.permutation(n_group))
        blocks: dict[int, list[str]] = {}
        for j, (subj, prof) in enumerate(zip(subjects, profiles)):
            sid = f"{g}{j + 1:02d}"
            inj = order + run[j]
            rows.append(mu + prof + slope * inj / 50 + rng.normal(0, 0.2))     # 注入ごとの希釈も足す
            samples.append(SampleInfo(sample_id=sid, role="sample", group=g, batch=None,
                                      injection_order=float(inj + 1),
                                      subject_id=f"{g}-{subj}" if subject else None, include=True))
            blocks.setdefault(subj, []).append(sid)
        expected[g] = {frozenset(b) for b in blocks.values()} if replicate else set()
        order += n_group
    values = np.exp2(np.array(rows))
    study = StudyData(values=values, detected=np.ones(values.shape, dtype=bool),
                      sample_ids=[s.sample_id for s in samples],
                      feature_ids=[f"f{k}" for k in range(n_features)], samples=samples)
    return study, expected


def _found_blocks(out) -> dict:
    return {g["group"]: {frozenset(b) for b in g["blocks"]} for g in out["groups"]}


@pytest.mark.parametrize("layout,expected_layout", [("consecutive", "consecutive"),
                                                    ("interleaved", "interleaved"), ("random", "data")])
def test_replicate_blocks_finds_injections_of_the_same_animal(layout, expected_layout):
    """試料名に印が無くても、3 匹 × 3 注入の群は 3 本ずつのまとまりとして見つかる（並びが連続でも巡回でも無作為でも）。"""
    study, expected = _block_study(layout=layout)
    assert stats_data.replicate_name_warnings(study.sample_ids) == []      # 名前からは見つからない
    out = quality.replicate_blocks(study)
    assert out["status"] == "found" and out["k"] == 3
    assert _found_blocks(out) == expected
    for g in out["groups"]:
        assert g["layout"] == expected_layout and g["order"] == "injection_order" and g["cells"] == "detected"
        assert g["margin"] >= 0.3 and g["purity"] >= 0.8 and g["median_icc"] >= 0.4


def test_replicate_blocks_six_groups_of_three_animals():
    """報告のあった形（6 群 × 9 注入 = 3 匹 × 3 注入）。"""
    study, expected = _block_study(groups=tuple("ABCDEF"), seed=3)
    out = quality.replicate_blocks(study)
    assert out["status"] == "found" and out["k"] == 3 and len(out["groups"]) == 6
    assert _found_blocks(out) == expected


@pytest.mark.parametrize("seed", range(5))
def test_replicate_blocks_does_not_fire_on_independent_samples(seed):
    study, _ = _block_study(replicate=False, seed=seed)
    assert quality.replicate_blocks(study)["status"] == "not_found"


def test_replicate_blocks_does_not_mistake_injection_drift_for_replicates():
    """独立な試料でも、注入順のドリフトは並びで連続の注入を似させる。境目をまたぐ近い対も似るので鳴らない
    （群が 1 つだけでも。そのときはまとまりの平均の差が並びの直線で説明できることも見る）。"""
    for groups in (("A", "B"), ("A",)):
        for seed in range(5):
            study, _ = _block_study(groups=groups, replicate=False, drift=6.0, seed=seed)
            assert quality.replicate_blocks(study)["status"] == "not_found", (groups, seed)


def test_replicate_blocks_single_group_still_found_and_reports_the_trend_share():
    """群が 1 つだけでも、同じ個体の注入のまとまりは見つかる（`untrended` が `single_group_untrended` 以上）。"""
    study, expected = _block_study(groups=("A",), seed=4)
    out = quality.replicate_blocks(study)
    assert out["status"] == "found" and _found_blocks(out) == expected
    assert out["groups"][0]["untrended"] >= quality.REPLICATE_BLOCK_THRESHOLDS["single_group_untrended"]


def test_replicate_blocks_skips_when_subject_id_is_set():
    study, _ = _block_study(subject=True)
    assert quality.replicate_blocks(study) == {"status": "subject_id_set"}


def test_replicate_blocks_uses_only_features_detected_in_every_injection():
    """gap-fill のセルを含む特徴は使わない。すべての注入で実検出の特徴が足りなければ見ない（not_evaluable）。"""
    study, _ = _block_study()
    det = np.ones(study.values.shape, dtype=bool)
    det[0, :] = False                          # A の 1 本の注入がすべて gap-fill
    out = quality.replicate_blocks(replace(study, detected=det))
    assert out["status"] == "found" and [g["group"] for g in out["groups"]] == ["B"]
    assert out["skipped"] == [{"group": "A", "n": 9, "reason": "too_few_features", "n_features": 0}]
    det[9, :] = False                          # B も同じなら、見られる群が無い
    assert quality.replicate_blocks(replace(study, detected=det))["status"] == "not_evaluable"


def test_replicate_blocks_needs_three_blocks_and_consistency_across_groups():
    """まとまりが 3 つ以上にならない群（2 匹 × 2 注入）は見ない。片方の群だけがまとまるなら鳴らさない。"""
    out = quality.replicate_blocks(_block_study(animals=2, reps=2)[0])
    assert out["status"] == "not_evaluable"
    assert {s["reason"] for s in out["skipped"]} == {"no_block_size"}
    rep, _ = _block_study(groups=("A",), seed=1)
    ind, _ = _block_study(groups=("B",), replicate=False, seed=2)
    mixed = replace(rep, values=np.vstack([rep.values, ind.values]), samples=rep.samples + ind.samples,
                    detected=np.ones((18, rep.values.shape[1]), dtype=bool),
                    sample_ids=rep.sample_ids + ind.sample_ids)
    assert quality.replicate_blocks(mixed)["status"] == "not_found"


def test_quality_verdict_cautions_on_replicate_blocks():
    ok_rsd = {"status": "ok", "median": 0.1, "n_qc": 4}
    none_failed = {"checked": 4, "failed": []}
    assert quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[], confounding=_conf(),
                                   replicate_warnings=[], replicate_blocks={"status": "found"}) == \
        ("caution", ["possible_replicate_blocks"])
    assert quality.quality_verdict(rsd=ok_rsd, failed=none_failed, drift=[], confounding=_conf(),
                                   replicate_warnings=[], replicate_blocks={"status": "not_found"}) == ("ok", [])
