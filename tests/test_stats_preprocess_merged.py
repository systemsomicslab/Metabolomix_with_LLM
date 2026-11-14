"""陽陰を統合したデータセットの前処理（極性ごとのドリフト補正と正規化）と品質の数え方。

陽イオンと陰イオンは別の注入・別の測定の列なので、統合したデータセット（`StudyData.feature_polarity`）の
ドリフト補正はその極性の注入順・バッチで、正規化の係数はその極性の特徴だけから求める。証明の柱は
「統合したものを前処理した値 = 極性ごとに前処理して並べた値」と、statsmodels・numpy での独立の計算。
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from statsmodels.nonparametric.smoothers_lowess import lowess

from metabolomix.analysis import overview
from metabolomix.analysis.stats import preprocess as P
from metabolomix.analysis.stats import quality
from metabolomix.analysis.stats.data import (SampleInfo, StudyData, from_dataset, polarity_samples,
                                             run_of)
from metabolomix.core.atomic_io import DomainError

OFF_FILTERS = {"blank_ratio": None, "detection": None, "qc_rsd": None, "d_ratio": None}

# 15 注入: QC 6 本（両端と間）・生体 8 本（2 群）・ブランク 1 本。
ROLES = ["qc", "sample", "sample", "qc", "sample", "blank", "qc", "sample", "sample", "qc",
         "sample", "sample", "qc", "sample", "qc"]
GROUPS = [None, "A", "B", None, "A", None, None, "B", "A", None, "B", "A", None, "B", None]
POS_ORDER = list(range(1, 16))
# 陰イオンの測定は別の列: QC は両端に残し、生体の順を入れ替える
NEG_ORDER = [1, 9, 3, 5, 13, 2, 8, 4, 11, 10, 6, 7, 12, 14, 15]
IDS = [f"S{i:02d}" for i in range(15)]


def _samples(order, *, batch=None, source="mddata", flagged=False):
    return [SampleInfo(sample_id=sid, role=role, group=grp, batch=batch, injection_order=float(o),
                       subject_id=None, include=True, injection_order_source=source,
                       injection_order_may_be_file_order=flagged)
            for sid, role, grp, o in zip(IDS, ROLES, GROUPS, order)]


def _values(seed=11):
    """陽イオン 4 特徴（~1e6）と陰イオン 3 特徴（~1e3）。注入ごとの効き（生体試料の量と、極性ごとの注入の
    揺れ）と、極性ごとの注入順に沿ったドリフトを持つ。陽イオンの最後の特徴は QC が揺れる（QC RSD の足切りで落ちる）。"""
    rng = np.random.default_rng(seed)
    bio = np.array([r == "sample" for r in ROLES])
    amount = np.where(bio, rng.lognormal(0, 0.15, 15), 1.0)       # 試料の量（両極性に共通。QC はプール）
    inj_pos, inj_neg = rng.lognormal(0, 0.05, 15), rng.lognormal(0, 0.05, 15)
    drift_pos = 1 + 0.03 * np.array(POS_ORDER)
    drift_neg = 1 - 0.025 * np.array(NEG_ORDER)
    spread = np.where(bio, 0.3, 0.02)[:, None]                       # 生体の個体差・QC の測定の揺れ
    pos = (rng.uniform(5e5, 2e6, 4)[None, :] * (amount * inj_pos * drift_pos)[:, None]
           * rng.lognormal(0, np.broadcast_to(spread, (15, 4))))
    pos[~bio, 3] *= rng.lognormal(0, 0.6, int((~bio).sum()))
    neg = (rng.uniform(5e2, 2e3, 3)[None, :] * (amount * inj_neg * drift_neg)[:, None]
           * rng.lognormal(0, np.broadcast_to(spread, (15, 3))))
    blank = ROLES.index("blank")
    pos[blank] *= 0.01
    neg[blank] *= 0.01
    return pos, neg


POS_IDS = [f"positive:p{j}" for j in range(4)]
NEG_IDS = [f"negative:n{j}" for j in range(3)]


def _merged(*, neg_runs="neg", seed=11, flagged_neg=False):
    """統合した StudyData。`neg_runs`: "neg"（陰イオンの順）/ "same"（陽イオンと同じ）/ "unknown"。"""
    pos, neg = _values(seed)
    main = _samples(POS_ORDER)
    if neg_runs == "unknown":
        runs = [None] * 15
    else:
        table = _samples(NEG_ORDER if neg_runs == "neg" else POS_ORDER, flagged=flagged_neg)
        runs = [run_of(s) for s in table]
    return StudyData(values=np.hstack([pos, neg]), detected=None, sample_ids=list(IDS),
                     feature_ids=POS_IDS + NEG_IDS, samples=main, dataset_hash="merged",
                     feature_polarity={**{f: "positive" for f in POS_IDS},
                                       **{f: "negative" for f in NEG_IDS}},
                     polarity_runs={"negative": runs})


def _alone(study: StudyData, polarity: str) -> StudyData:
    """その極性だけの StudyData（その極性の測定の試料の表で）。"""
    cols = [j for j, f in enumerate(study.feature_ids) if study.feature_polarity[f] == polarity]
    return StudyData(values=study.values[:, cols], detected=None, sample_ids=list(study.sample_ids),
                     feature_ids=[study.feature_ids[j] for j in cols],
                     samples=polarity_samples(study, polarity), dataset_hash=polarity)


def _step(prep, name):
    return next(s for s in prep.steps if s["step"] == name)


def _columns(prep, ids):
    pos = {f: j for j, f in enumerate(prep.feature_ids)}
    return prep.log2[:, [pos[f] for f in ids]]


@pytest.mark.parametrize("method", ["tic", "median", "pqn", "auto"])
@pytest.mark.parametrize("filters", [OFF_FILTERS, None])
def test_merged_equals_each_polarity_alone_then_concatenated(method, filters):
    """統合したものの前処理 = 極性ごとに前処理して並べたもの（ドリフト補正・足切り・正規化・log2）。"""
    study = _merged()
    cfg = {"drift": {"mode": "on"}, "normalize": {"method": method}}
    if filters is not None:
        cfg["filters"] = filters
    merged = P.preprocess(study, cfg)
    alone = {p: P.preprocess(_alone(study, p), cfg) for p in ("positive", "negative")}
    kept = [f for p in ("positive", "negative") for f in alone[p].feature_ids]
    assert merged.feature_ids == kept
    assert [u.unit_id for u in merged.units] == [u.unit_id for u in alone["positive"].units]
    for p in ("positive", "negative"):
        np.testing.assert_allclose(_columns(merged, alone[p].feature_ids), alone[p].log2,
                                   rtol=1e-12, atol=0)
        qc_cols = [merged.feature_ids.index(f) for f in alone[p].feature_ids]
        np.testing.assert_allclose(merged.qc_log2[:, qc_cols], alone[p].qc_log2, rtol=1e-12, atol=0)
        drift = _step(merged, "drift")["by_polarity"][p]
        assert drift == {k: v for k, v in _step(alone[p], "drift").items() if k != "step"}
        norm = _step(merged, "normalize")["by_polarity"][p]
        assert norm["factors"] == pytest.approx(_step(alone[p], "normalize")["factors"], rel=1e-12)
    assert _step(merged, "drift")["status"] == "applied"
    assert _step(merged, "normalize")["method"] == ("pqn" if method == "auto" else method)


def test_merged_negative_block_matches_independent_lowess_and_tic():
    """陰イオンの値を statsmodels の lowess（陰イオンの注入順）と numpy の TIC から計算し直して合わせる。
    陽イオンの強度（~1e6）は陰イオン（~1e3）の係数に効かない。"""
    study = _merged()
    prep = P.preprocess(study, {"drift": {"mode": "on"}, "normalize": {"method": "tic"},
                                "filters": OFF_FILTERS, "collapse_replicates": False})
    rows = [i for i, r in enumerate(ROLES) if r != "blank"]
    qc = [i for i in rows if ROLES[i] == "qc"]
    x_all = np.array([NEG_ORDER[i] for i in rows], dtype=float)
    x_qc = np.array([NEG_ORDER[i] for i in qc], dtype=float)
    neg = study.values[:, 4:]
    corrected = np.empty((len(rows), 3))
    for j in range(3):
        fit = lowess(neg[qc, j], x_qc, frac=0.75, it=0, delta=0.0, xvals=x_all)
        corrected[:, j] = neg[rows, j] / fit * np.median(neg[qc, j])
    # バッチは 1 つなので水準合わせの係数は 1。TIC の係数 = 陰イオンの特徴の和 / その中央値
    totals = corrected.sum(axis=1)
    expected = corrected / (totals / np.median(totals))[:, None]
    bio = [k for k, i in enumerate(rows) if ROLES[i] == "sample"]
    np.testing.assert_allclose(prep.sample_values[:, 4:], expected[bio], rtol=1e-10)
    np.testing.assert_allclose(_step(prep, "normalize")["by_polarity"]["negative"]["factors"],
                               totals / np.median(totals), rtol=1e-10)

    # 陽陰を 1 つの塊として扱う（極性を知らない）前処理では合わない（この試験が直したものを捉えている）
    blind = P.preprocess(replace(study, feature_polarity={}, polarity_runs={}),
                         {"drift": {"mode": "on"}, "normalize": {"method": "tic"},
                          "filters": OFF_FILTERS, "collapse_replicates": False})
    assert not np.allclose(blind.sample_values[:, 4:], expected[bio], rtol=1e-3)


def test_merged_records_per_polarity_drift_and_factors():
    prep = P.preprocess(_merged(), {"drift": {"mode": "on"}, "normalize": {"method": "median"}})
    drift, norm = _step(prep, "drift"), _step(prep, "normalize")
    assert set(drift["by_polarity"]) == {"positive", "negative"}
    assert all(st["status"] == "applied" and st["method"] == "loess"
               for st in drift["by_polarity"].values())
    assert norm["status"] == "applied" and norm["factor_samples"] == [
        s for s, r in zip(IDS, ROLES) if r != "blank"]
    assert set(norm["by_polarity"]) == {"positive", "negative"}
    assert norm["by_polarity"]["positive"]["factors"] != norm["by_polarity"]["negative"]["factors"]
    assert norm["factors_set_to_one"] == 0


def test_merged_warns_when_polarity_runs_differ():
    prep = P.preprocess(_merged(), {"drift": {"mode": "on"}})
    warn = next(w for w in prep.warnings if w["code"] == "polarity_runs_differ")
    differs = [sid for sid, r, a, b in zip(IDS, ROLES, POS_ORDER, NEG_ORDER)
               if r in ("sample", "qc") and a != b]
    assert warn["polarity"] == "negative" and warn["reference_polarity"] == "positive"
    assert warn["n_injection_order_differs"] == len(differs) and warn["n_batch_differs"] == 0
    assert warn["samples"] == differs
    same = P.preprocess(_merged(neg_runs="same"), {"drift": {"mode": "on"}})
    assert not [w for w in same.warnings if w["code"] == "polarity_runs_differ"]


def test_merged_unknown_negative_order_skips_only_negative_drift():
    study = _merged(neg_runs="unknown")
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    drift = _step(prep, "drift")
    assert drift["status"] == "partial" and drift["reasons"] == ["no_injection_order"]
    assert drift["polarities_skipped"] == ["negative"] and "rerun_with" not in drift
    assert drift["by_polarity"]["positive"]["status"] == "applied"
    warn = next(w for w in prep.warnings if w["code"] == "polarity_runs_differ")
    assert warn["n_injection_order_unknown"] == sum(r != "blank" for r in ROLES)
    # 陰イオンの値はドリフト補正なし（その極性だけを前処理したのと同じ）
    alone = P.preprocess(_alone(study, "negative"), {"filters": OFF_FILTERS})
    np.testing.assert_allclose(_columns(prep, NEG_IDS), alone.log2, rtol=1e-12)
    with pytest.raises(DomainError) as exc:
        P.preprocess(study, {"drift": {"mode": "on"}, "filters": OFF_FILTERS})
    assert exc.value.code == "PREPROCESS_PRECONDITION"
    assert exc.value.details["polarity"] == "negative"
    assert exc.value.message.startswith("negative: ")


def test_merged_auto_holds_only_the_polarity_with_file_order():
    prep = P.preprocess(_merged(flagged_neg=True), {"filters": OFF_FILTERS})
    drift = _step(prep, "drift")
    assert drift["status"] == "partial"
    assert drift["reasons"] == ["injection_order_may_be_file_order"]
    assert drift["rerun_with"] == {"preprocess.drift.mode": "on"}
    assert drift["by_polarity"]["negative"]["rerun_with"] == {"preprocess.drift.mode": "on"}
    hold = next(w for w in prep.warnings if w["code"] == "drift_not_applied_by_default")
    assert hold["polarities"] == {"negative": ["injection_order_may_be_file_order"]}
    assert hold["suggestion"] == {"preprocess.drift.mode": "on"}
    order = next(w for w in prep.warnings if w["code"] == "injection_order_may_be_file_order")
    assert order["polarity"] == "negative"
    unflagged = P.preprocess(_merged(), {"filters": OFF_FILTERS})
    assert unflagged.matrix_id != prep.matrix_id


def test_merged_normalization_features_fall_back_per_polarity():
    prep = P.preprocess(_merged(), {"drift": {"mode": "off"}, "filters": OFF_FILTERS,
                                    "normalize": {"method": "tic", "features": ["positive:p0"]}})
    by = _step(prep, "normalize")["by_polarity"]
    assert by["positive"]["factor_feature_basis"] == "features" and by["positive"]["factor_features"] == 1
    assert by["negative"]["factor_feature_basis"] == "kept_features"
    assert by["negative"]["factor_features"] == 3
    note = next(w for w in prep.warnings if w["code"] == "normalization_features_not_in_polarity")
    assert note["polarities"] == ["negative"]
    # 陽イオンの係数は p0 だけから（numpy で）
    rows = [i for i, r in enumerate(ROLES) if r != "blank"]
    p0 = _merged().values[rows, 0]
    np.testing.assert_allclose(by["positive"]["factors"], p0 / np.median(p0), rtol=1e-12)
    with pytest.raises(DomainError) as exc:
        P.preprocess(_merged(), {"normalize": {"method": "tic", "features": ["nope"]}})
    assert exc.value.code == "PREPROCESS_CONFIG_INVALID"


def test_merged_negative_batches_follow_the_negative_run():
    """陰イオンだけが 2 バッチ（陽イオンは 1 バッチ）: 陰イオンのドリフトの可否は陰イオンのバッチで決まり
    （QC が足りず当てない）、陽イオンは補正する。値は極性ごとに前処理して並べたものと同じ。"""
    study = _merged()
    runs = [dict(r, batch="b1" if r["injection_order"] <= 8 else "b2") for r in study.polarity_runs["negative"]]
    study = replace(study, polarity_runs={"negative": runs})
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    drift = _step(prep, "drift")
    assert drift["status"] == "partial" and drift["polarities_skipped"] == ["negative"]
    assert drift["by_polarity"]["negative"]["reasons"] == ["too_few_qc:b1", "too_few_qc:b2"]
    for p in ("positive", "negative"):
        alone = P.preprocess(_alone(study, p), {"filters": OFF_FILTERS})
        np.testing.assert_allclose(_columns(prep, alone.feature_ids), alone.log2, rtol=1e-12)
    warn = next(w for w in prep.warnings if w["code"] == "polarity_runs_differ")
    assert warn["n_batch_differs"] == sum(r != "blank" for r in ROLES)


def test_merged_with_excluded_sample_matches_each_polarity():
    """除外した試料（include=False）があっても、陰イオンの測定は試料に揃ったまま。"""
    study = _merged()
    gone = IDS.index("S07")
    study.samples[gone] = replace(study.samples[gone], include=False)
    cfg = {"drift": {"mode": "on"}, "normalize": {"method": "pqn"}, "filters": OFF_FILTERS}
    prep = P.preprocess(study, cfg)
    assert "S07" not in prep.sample_ids
    for p in ("positive", "negative"):
        alone = P.preprocess(_alone(study, p), cfg)
        np.testing.assert_allclose(_columns(prep, alone.feature_ids), alone.log2, rtol=1e-12)


def test_merged_matrix_id_tracks_negative_runs():
    a = P.preprocess(_merged(), {})
    assert a.matrix_id == P.preprocess(_merged(), {}).matrix_id
    assert a.matrix_id != P.preprocess(_merged(neg_runs="same"), {}).matrix_id
    blind = replace(_merged(), feature_polarity={}, polarity_runs={})
    assert a.matrix_id != P.preprocess(blind, {}).matrix_id


def test_non_merged_matrix_id_is_unchanged(tmp_path):
    """統合でないデータセットは極性の塊を持たず、matrix_id は極性の塊を足す前（48c9c7e）と同じ固定値。"""
    pos, _neg = _values()
    study = StudyData(values=pos, detected=None, sample_ids=list(IDS), feature_ids=["a", "b", "c", "d"],
                      samples=_samples(POS_ORDER), dataset_hash="plain")
    prep = P.preprocess(study, {})
    assert prep.matrix_id == "mx_f3fd5e9470bf0f0c"
    assert [s.get("by_polarity") for s in prep.steps] == [None] * len(prep.steps)
    assert not [w for w in prep.warnings if "polarity" in w]
    # 解析フォルダから読んだ（統合でない）データセットは極性の塊を作らない
    plain = from_dataset(_folder(tmp_path / "plain", POS_ORDER, pos, "PC"))
    assert plain.feature_polarity == {} and plain.polarity_runs == {}


# --- 品質（概観）-------------------------------------------------------------

def test_failed_qc_is_counted_per_polarity():
    """陰イオンだけで失敗した QC 注入は、陽イオンの強度に隠れずに見つかる。"""
    study = _merged()
    bad = IDS.index("S06")
    study.values[bad, 4:] *= 0.05
    out = quality.failed_qc(study)
    assert out["failed"] == ["S06"]
    assert out["by_polarity"] == {"positive": [], "negative": ["S06"]}
    blind = quality.failed_qc(replace(study, feature_polarity={}, polarity_runs={}))
    assert blind["failed"] == []                 # 1 つの塊の総量では陽イオンに隠れる


def test_drift_signs_use_each_polarity_order():
    from scipy import stats as sps
    study = _merged()
    rows = quality.drift_signs(study)
    assert [r["polarity"] for r in rows] == ["positive", "negative"]
    qc = [i for i, r in enumerate(ROLES) if r == "qc"]
    totals = study.values[np.ix_(qc, range(4, 7))].sum(axis=1)
    rho, p = sps.spearmanr([NEG_ORDER[i] for i in qc], np.log2(totals))
    assert rows[1]["rho"] == pytest.approx(rho) and rows[1]["p"] == pytest.approx(p)
    facts = overview.quality_facts(rsd={"status": "not_evaluable"}, failed={"failed": []},
                                   detection={"status": "not_evaluable"},
                                   drift=[{**rows[1], "flag": True}], confounding={},
                                   replicate_warnings=[], verdict="caution", reasons=[], n_qc=6,
                                   n_blank=1, drift_note=None, pca_id=None)
    assert any("negative のバッチ" in item["text"] for item in facts["facts"])


# --- 統合（quant.merge）から統計エンジンまで -----------------------------------

def _point(tmp_path, monkeypatch):
    from tests.quant_fake_rules import write_fake_rules
    folder = write_fake_rules(tmp_path / "rules")
    cfg = tmp_path / "lipidmix.local.toml"
    cfg.write_text(f"[quant]\nrules = '{folder}'\n", encoding="utf-8")
    monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))


def _folder(path, order, values, ion, *, permute=None, names=None):
    from metabolomix.dataset.alignment_files import select_alignment
    from metabolomix.dataset.msdial_adapter import load_msdial_folder
    from tests.msdial_folder_fixtures import write_msdial_folder
    idx = list(range(15)) if permute is None else permute
    names = list(IDS) if names is None else names
    samples = [{"file_id": k, "name": names[i], "class_id": GROUPS[i] or ROLES[i].upper(),
                "file_type": {"sample": 0, "qc": 2, "blank": 3}[ROLES[i]], "order": order[i]}
               for k, i in enumerate(idx)]
    spots = [{"name": f"{ion} {30 + j}:1", "ontology": ion, "mz": 500.0 + j,
              "heights": [float(values[i, j]) for i in idx]} for j in range(values.shape[1])]
    write_msdial_folder(path, spots=spots, samples=samples)
    return load_msdial_folder(select_alignment(path))


def test_merge_keeps_negative_runs_and_preprocess_matches_each_polarity(tmp_path, monkeypatch):
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    pos = _folder(tmp_path / "run" / "pos", POS_ORDER, pos_v, "PC")
    # 陰イオンの解析フォルダは試料の並び（FileID）も違う
    neg = _folder(tmp_path / "run" / "neg", NEG_ORDER, neg_v, "PE",
                  permute=[3, 0, 1, 2, 5, 4, 6, 8, 7, 9, 11, 10, 12, 14, 13])
    merged, info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    assert info["negative_runs"] == {"injection_order": "different", "batch": "same",
                                     "n_samples_differ": sum(a != b for a, b in zip(POS_ORDER, NEG_ORDER))}
    runs = merged.quant_state["polarity_runs"]["negative"]
    assert {sid: runs[sid]["injection_order"] for sid in IDS} == {
        sid: float(o) for sid, o in zip(IDS, NEG_ORDER)}
    assert merged.quant_state["polarity_runs_match"] == {"negative": info["negative_runs"]}

    study = from_dataset(merged)
    assert set(study.feature_polarity.values()) == {"positive", "negative"}
    assert [s.injection_order for s in polarity_samples(study, "negative")] == [
        float(o) for o in NEG_ORDER]
    assert [s.injection_order for s in study.samples] == [float(o) for o in POS_ORDER]

    cfg = {"drift": {"mode": "on"}, "normalize": {"method": "tic"}}
    prep = P.preprocess(study, cfg)
    for prefix, ds in (("positive", pos), ("negative", neg)):
        alone = P.preprocess(from_dataset(ds), cfg)
        by_unit = {u.unit_id: k for k, u in enumerate(alone.units)}
        order = [by_unit[u.unit_id] for u in prep.units]
        ids = [f"{prefix}:{f}" for f in alone.feature_ids]
        np.testing.assert_allclose(_columns(prep, ids), alone.log2[order], rtol=1e-12)
    assert any(w["code"] == "polarity_runs_differ" for w in prep.warnings)


def test_merge_records_same_runs(tmp_path, monkeypatch):
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    pos = _folder(tmp_path / "run" / "pos", POS_ORDER, pos_v, "PC")
    neg = _folder(tmp_path / "run" / "neg", POS_ORDER, neg_v, "PE")
    merged, info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    assert info["negative_runs"] == {"injection_order": "same", "batch": "same", "n_samples_differ": 0}
    prep = P.preprocess(from_dataset(merged), {"drift": {"mode": "on"}})
    assert not [w for w in prep.warnings if w["code"] == "polarity_runs_differ"]


def test_tools_preprocess_and_overview_on_merged_handle(tmp_path, monkeypatch):
    """dataset_merge_polarities → dataset_preprocess → dataset_overview（ツールの戻り値に極性ごとの記録が出る）。"""
    import json

    from metabolomix.core import handles, session_state
    from metabolomix.quant import merge as qm
    from metabolomix.tools.overview_tools import dataset_overview
    from metabolomix.tools.stats_tools import dataset_preprocess
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    neg_v[IDS.index("S06")] *= 0.05                     # 陰イオンだけで失敗した QC 注入
    pos = _folder(tmp_path / "run" / "pos", POS_ORDER, pos_v, "PC")
    neg = _folder(tmp_path / "run" / "neg", NEG_ORDER, neg_v, "PE")
    ids = []
    for ds, polarity in ((pos, "positive"), (neg, "negative")):
        session_state.session.dataset = ds
        holder = handles.registry.get(session_state.session.current_handle_id)
        holder.polarity = polarity
        ids.append(holder.handle_id)
    merged = json.loads(qm.merge_tool(handles.registry.get(ids[0]), handles.registry.get(ids[1]),
                                      None, None))
    assert merged["status"] == "merged" and merged["negative_runs"]["injection_order"] == "different"

    out = json.loads(dataset_preprocess(drift_mode="on", normalize="tic", max_d_ratio=None,
                                        max_qc_rsd=None))
    assert out["status"] == "success"
    steps = {s["step"]: s for s in out["steps"]}
    assert steps["drift"]["status"] == "applied"
    assert set(steps["drift"]["by_polarity"]) == {"positive", "negative"}
    assert steps["normalize"]["by_polarity"]["negative"]["n_factors"] == 14   # ブランクを除く注入
    assert "polarity_runs_differ" in {w["code"] for w in out["warnings"]}

    view = json.loads(dataset_overview())
    assert any("negative S06" in fact for fact in view["quality"]["facts"])


@pytest.mark.parametrize("method", ["tic", "median", "pqn"])
def test_v1_recipe_normalizes_each_polarity(tmp_path, monkeypatch, method):
    """analysis-spec.v1 の recipe の正規化も、統合したデータセットでは係数を極性ごとに求める
    （統合して回した行 = 極性ごとに回した行。係数は numpy の TIC と合う）。"""
    from metabolomix.analysis import analysis_spec as A
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    pos = _folder(tmp_path / "run" / "pos", POS_ORDER, pos_v, "PC")
    neg = _folder(tmp_path / "run" / "neg", NEG_ORDER, neg_v, "PE")
    merged, _info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    groups = tmp_path / "groups.tsv"
    groups.write_text("sample\tgroup\trole\n" + "".join(
        f"{sid}\t{grp or ''}\t{role}\n" for sid, grp, role in zip(IDS, GROUPS, ROLES)), encoding="utf-8")
    recipe = {"steps": ["normalize"], "normalize": {"method": method, "reference": "qc"}}

    def run(ds):
        spec = {"input": {"source": "dataset"},
                "groups": {"source": "tsv", "path": str(groups), "role_column": "role"},
                "comparison": {"reference": "A", "test": "B"}, "features": {"merge_split": False},
                "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}, "matrix": recipe}]}
        return A.run_spec(A.normalize_spec(spec), ds)["analyses"][0]

    both = run(merged)
    rows = {r["feature_id"]: r for r in both["rows"]}
    for prefix, ds in (("positive", pos), ("negative", neg)):
        alone = run(ds)
        for r in alone["rows"]:
            mine = rows[f"{prefix}:{r['feature_id']}"]
            for key in ("log2fc_test_over_ref", "p"):
                assert mine[key] == pytest.approx(r[key], rel=1e-10)
        step = next(s for s in both["matrix"]["steps"] if s["step"] == "normalize")
        own = next(s for s in alone["matrix"]["steps"] if s["step"] == "normalize")
        assert step["by_polarity"][prefix]["factors"] == pytest.approx(own["factors"], rel=1e-12)
    if method == "tic":
        step = next(s for s in both["matrix"]["steps"] if s["step"] == "normalize")
        totals = neg_v.sum(axis=1)
        expected = totals / np.median(totals)
        got = [step["by_polarity"]["negative"]["factors"][sid] for sid in IDS]
        np.testing.assert_allclose(got, expected, rtol=1e-12)


# --- レビュー指摘: 統合は陽イオンの試料の表（シート）を引き継ぎ、陽陰の測定を前処理が読む値で比べる ----------

#: .mddata の注入順（MS-DIAL の既定＝ファイルの順）と、シートで入れた実際の順（陽イオン・陰イオン）。
FILE_ORDER = list(range(1, 16))
POS_TRUE = [1, 9, 3, 5, 13, 2, 8, 4, 11, 10, 6, 7, 12, 14, 15]
NEG_TRUE = [1, 3, 5, 7, 2, 4, 9, 6, 8, 11, 10, 13, 12, 14, 15]
_SHEET_HEADER = "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude"


def _sheet_text(names, order) -> str:
    """raw（fixture の .mddata の file_path は C:/data/<試料名>.wiff）を挙げた sample-manifest.v1。"""
    lines = ["# schema = sample-manifest.v1", _SHEET_HEADER]
    for i, name in enumerate(names):
        lines.append("\t".join([str(name), f"{name}.wiff", ROLES[i], GROUPS[i] or "", "", str(order[i]),
                                "", ""]))
    return "\n".join(lines) + "\n"


def _sheet(ds, order, path):
    from metabolomix.analysis import sample_manifest
    path.write_text(_sheet_text(ds.sample_names, order), encoding="utf-8")
    return sample_manifest.apply_sample_manifest(ds, str(path))


def _sheeted_pair(tmp_path, pos_order, neg_order):
    pos_v, neg_v = _values()
    pos = _folder(tmp_path / "run" / "pos", FILE_ORDER, pos_v, "PC")
    neg = _folder(tmp_path / "run" / "neg", FILE_ORDER, neg_v, "PE")
    _sheet(pos, pos_order, tmp_path / "pos.tsv")
    _sheet(neg, neg_order, tmp_path / "neg.tsv")
    return pos, neg, pos_v, neg_v


def _matches_each_polarity(prep, pos, neg, cfg):
    for prefix, ds in (("positive", pos), ("negative", neg)):
        alone = P.preprocess(from_dataset(ds), cfg)
        by_unit = {u.unit_id: k for k, u in enumerate(alone.units)}
        order = [by_unit[u.unit_id] for u in prep.units]
        ids = [f"{prefix}:{f}" for f in alone.feature_ids]
        np.testing.assert_allclose(_columns(prep, ids), alone.log2[order], rtol=1e-12)
        assert _step(prep, "drift")["by_polarity"][prefix] == {
            k: v for k, v in _step(alone, "drift").items() if k != "step"}


def test_merge_carries_positive_sheet_and_blocks_follow_each_sheet(tmp_path, monkeypatch):
    """陽陰の両方にシート（.mddata のファイルの順と違う実際の順）を当ててから統合: 統合したデータセットの
    試料の表は陽イオンのシート、陰イオンは陰イオンのシートの順で、各極性の値はその極性だけを前処理したものと同じ。
    陽イオンの値は statsmodels の lowess（陽イオンのシートの順）と numpy の TIC から計算し直したものと合う。"""
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos, neg, pos_v, _neg_v = _sheeted_pair(tmp_path, POS_TRUE, NEG_TRUE)
    merged, info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    assert info["negative_runs"] == {"injection_order": "different", "batch": "same",
                                     "n_samples_differ": sum(a != b for a, b in zip(POS_TRUE, NEG_TRUE))}
    assert merged.sample_metadata_rows == pos.sample_metadata_rows      # 陽イオンのシートを引き継ぐ
    assert merged.quant_state["record_id"] is not None
    study = from_dataset(merged)
    assert [s.injection_order for s in study.samples] == [float(o) for o in POS_TRUE]
    assert {s.injection_order_source for s in study.samples} == {"sheet"}
    assert [s.injection_order for s in polarity_samples(study, "negative")] == [float(o) for o in NEG_TRUE]

    for cfg in ({"drift": {"mode": "on"}, "normalize": {"method": "tic"}, "filters": OFF_FILTERS},
                {"drift": {"mode": "on"}, "normalize": {"method": "pqn"}}, {}):
        prep = P.preprocess(study, cfg)
        _matches_each_polarity(prep, pos, neg, cfg)
    # 既定（drift=auto）: シートの順は確かなので両極性とも補正する
    assert _step(P.preprocess(study, {}), "drift")["status"] == "applied"

    # 陽イオンの塊を独立に計算し直す（陽イオンのシートの順で LOESS、陽イオンの特徴の和で TIC）
    prep = P.preprocess(study, {"drift": {"mode": "on"}, "normalize": {"method": "tic"},
                                "filters": OFF_FILTERS, "collapse_replicates": False})
    rows = [i for i, r in enumerate(ROLES) if r != "blank"]
    qc = [i for i in rows if ROLES[i] == "qc"]
    bio = [k for k, i in enumerate(rows) if ROLES[i] == "sample"]

    def expected(order):
        x_all = np.array([order[i] for i in rows], dtype=float)
        x_qc = np.array([order[i] for i in qc], dtype=float)
        out = np.empty((len(rows), 4))
        for j in range(4):
            fit = lowess(pos_v[qc, j], x_qc, frac=0.75, it=0, delta=0.0, xvals=x_all)
            out[:, j] = pos_v[rows, j] / fit * np.median(pos_v[qc, j])
        totals = out.sum(axis=1)
        return (out / (totals / np.median(totals))[:, None])[bio]

    np.testing.assert_allclose(prep.sample_values[:, :4], expected(POS_TRUE), rtol=1e-10)
    # 統合がシートを落とした（陽イオンを .mddata のファイルの順で当てはめた）値とは合わない
    assert not np.allclose(prep.sample_values[:, :4], expected(FILE_ORDER), rtol=1e-3)


def test_merge_compares_the_runs_the_merged_dataset_reads(tmp_path, monkeypatch):
    """両極性に同じ実際の順のシート: 統合の negative_runs は same で、前処理も同じ順を読む
    （陽陰の違い・ファイルの順の注意は出ず、既定の drift=auto は両極性を補正する）。"""
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos, neg, _pv, _nv = _sheeted_pair(tmp_path, POS_TRUE, POS_TRUE)
    merged, info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    assert info["negative_runs"] == {"injection_order": "same", "batch": "same", "n_samples_differ": 0}
    assert merged.quant_state["polarity_runs_match"] == {"negative": info["negative_runs"]}
    prep = P.preprocess(from_dataset(merged), {"filters": OFF_FILTERS})
    codes = {w["code"] for w in prep.warnings}
    assert "polarity_runs_differ" not in codes and "injection_order_may_be_file_order" not in codes
    assert _step(prep, "drift")["status"] == "applied"
    _matches_each_polarity(prep, pos, neg, {"filters": OFF_FILTERS})


def test_merged_file_order_uses_positive_source_names_without_raw_paths(tmp_path, monkeypatch):
    """raw のパスが無く、統合が試料を共通 ID へ改名するとき、陽イオンの注入順がファイルの順かは陽イオンの
    元の試料名で見る（陽イオンだけのデータセットと同じ印・同じ drift=auto の判断）。"""
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    stems = [f"S{i}" for i in range(15)]           # S1 と S10 は _pos を付けると並びが入れ替わる
    pos_names, neg_names = [f"{n}_pos" for n in stems], [f"{n}_neg" for n in stems]
    rank = {n: k + 1 for k, n in enumerate(sorted(pos_names))}
    file_order = [rank[n] for n in pos_names]       # MS-DIAL の既定（陽イオンのファイル名の順）
    pos = _folder(tmp_path / "run" / "pos", file_order, pos_v, "PC", names=pos_names)
    neg = _folder(tmp_path / "run" / "neg", file_order, neg_v, "PE", names=neg_names)
    for ds in (pos, neg):
        for entry in ds.msdial_samples:
            entry["file_path"] = None
    merged, info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    assert info["rule"] == "polarity_token_removed" and merged.sample_names == stems
    assert merged.quant_state["source_sample_names"] == {"positive": pos_names, "negative": neg_names}
    merged_study, pos_study = from_dataset(merged), from_dataset(pos)
    flags = [s.injection_order_may_be_file_order for s in pos_study.samples]
    assert any(flags)
    assert [s.injection_order_may_be_file_order for s in merged_study.samples] == flags
    cfg = {"filters": OFF_FILTERS}
    prep = P.preprocess(merged_study, cfg)
    alone = P.preprocess(pos_study, cfg)
    assert _step(prep, "drift")["by_polarity"]["positive"] == {
        k: v for k, v in _step(alone, "drift").items() if k != "step"}
    assert _step(alone, "drift")["reason"] == "injection_order_may_be_file_order"


# --- レビュー指摘: 統合後に陰イオンの測定を入れ直す手順を、注意とシートの戻り値が示す ----------------------

def test_merged_warnings_name_the_negative_run_remedy():
    prep = P.preprocess(_merged(flagged_neg=True), {"filters": OFF_FILTERS})
    order = next(w for w in prep.warnings if w["code"] == "injection_order_may_be_file_order")
    hold = next(w for w in prep.warnings if w["code"] == "drift_not_applied_by_default")
    for warn in (order, hold):
        assert "negative のハンドルに dataset_set_sample_metadata" in warn["message"]
        assert "dataset_merge_polarities で統合し直す" in warn["message"]
    assert "その注入順のまま" in hold["message"]
    differ = next(w for w in P.preprocess(_merged(), {"drift": {"mode": "on"}}).warnings
                  if w["code"] == "polarity_runs_differ")
    assert "negative のハンドルに dataset_set_sample_metadata" in differ["message"]
    assert "統合したハンドルに dataset_set_sample_metadata" in P.run_remedy("positive")


def test_sheet_on_merged_handle_sets_only_positive_run_and_remerge_sets_negative(tmp_path, monkeypatch):
    """統合したハンドルへのシートは陽イオンの測定だけを決め（戻り値の notes が言う）、陰イオンの測定は
    陰イオンのハンドルにシートを当てて統合し直すと変わる（注意が示す手順が効く）。"""
    import json

    from metabolomix.core import handles, session_state
    from metabolomix.quant import merge as qm
    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata
    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    pos = _folder(tmp_path / "run" / "pos", FILE_ORDER, pos_v, "PC")
    neg = _folder(tmp_path / "run" / "neg", FILE_ORDER, neg_v, "PE")
    ids = []
    for ds, polarity in ((pos, "positive"), (neg, "negative")):
        session_state.session.dataset = ds
        holder = handles.registry.get(session_state.session.current_handle_id)
        holder.polarity = polarity
        ids.append(holder.handle_id)
    first = json.loads(qm.merge_tool(handles.registry.get(ids[0]), handles.registry.get(ids[1]),
                                     None, None))
    sheet = tmp_path / "true.tsv"
    sheet.write_text(_sheet_text(IDS, POS_TRUE), encoding="utf-8")

    out = json.loads(dataset_set_sample_metadata(str(sheet), dataset_handle=first["dataset_handle"]))
    assert out["status"] == "success"
    assert any("陽イオンの測定だけ" in n and "統合し直す" in n for n in out["notes"])
    study = from_dataset(handles.registry.get(first["dataset_handle"]).dataset)
    assert [s.injection_order for s in study.samples] == [float(o) for o in POS_TRUE]
    neg_run = polarity_samples(study, "negative")
    assert [s.injection_order for s in neg_run] == [float(o) for o in FILE_ORDER]
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    hold = next(w for w in prep.warnings if w["code"] == "drift_not_applied_by_default")
    assert hold["polarities"] == {"negative": ["injection_order_may_be_file_order"]}

    # 示された手順: 陰イオンのハンドルにシートを当てて統合し直す
    plain = json.loads(dataset_set_sample_metadata(str(sheet), dataset_handle=ids[1]))
    assert plain["status"] == "success" and "notes" not in plain     # 統合でなければ注記は無い
    assert json.loads(dataset_set_sample_metadata(str(sheet), dataset_handle=ids[0]))["status"] == "success"
    second = json.loads(qm.merge_tool(handles.registry.get(ids[0]), handles.registry.get(ids[1]),
                                      None, None))
    assert second["negative_runs"] == {"injection_order": "same", "batch": "same", "n_samples_differ": 0}
    study = from_dataset(handles.registry.get(second["dataset_handle"]).dataset)
    assert [s.injection_order for s in polarity_samples(study, "negative")] == [float(o) for o in POS_TRUE]
    assert _step(P.preprocess(study, {"filters": OFF_FILTERS}), "drift")["status"] == "applied"


# --- レビュー指摘: analysis-spec v1 の drift_correct も極性ごとの注入順で ------------------------------

def _v1_spec(recipe, order):
    return {"input": {"source": "dataset"}, "groups": {"source": "metadata"},
            "comparison": {"reference": "A", "test": "B"}, "features": {"merge_split": False},
            "injection_order": order,
            "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}, "matrix": recipe}]}


def test_v1_recipe_drift_correct_uses_each_polarity_run(tmp_path, monkeypatch):
    """injection_order.source=metadata: 統合したデータセットの drift_correct は、陽イオンを試料の表（陽イオンの
    シート）の順、陰イオンを統合が残した陰イオンの測定の順で当てはめる（統合して回した行 = 極性ごとに回した行）。
    陰イオンの行は statsmodels の lowess（陰イオンの順）と scipy の Welch の t から計算し直したものと合う。"""
    from scipy import stats as sps

    from metabolomix.analysis import analysis_spec as A
    from metabolomix.quant import merge as qm
    _point(tmp_path, monkeypatch)
    pos, neg, _pos_v, neg_v = _sheeted_pair(tmp_path, POS_TRUE, NEG_TRUE)
    merged, _info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    spec = A.normalize_spec(_v1_spec({"steps": ["drift_correct"]}, {"source": "metadata"}))
    notes: list[str] = []
    both = A.run_spec(spec, merged, notes=notes)["analyses"][0]
    assert notes == []
    rows = {r["feature_id"]: r for r in both["rows"]}
    for prefix, ds in (("positive", pos), ("negative", neg)):
        alone = A.run_spec(spec, ds)["analyses"][0]
        for r in alone["rows"]:
            mine = rows[f"{prefix}:{r['feature_id']}"]
            for key in ("log2fc_test_over_ref", "p", "statistic"):
                assert mine[key] == pytest.approx(r[key], rel=1e-10)
    drift = next(s for s in both["matrix"]["steps"] if s["step"] == "drift_correct")
    assert drift["by_polarity"] == {
        "positive": {"uncorrected_features": 0, "injection_order": "metadata"},
        "negative": {"uncorrected_features": 0, "injection_order": "polarity_runs.negative"}}

    qc = [i for i, r in enumerate(ROLES) if r == "qc"]
    x_all = np.array(NEG_TRUE, dtype=float)
    x_pos = np.array(POS_TRUE, dtype=float)
    a = [i for i, (r, g) in enumerate(zip(ROLES, GROUPS)) if r == "sample" and g == "A"]
    b = [i for i, (r, g) in enumerate(zip(ROLES, GROUPS)) if r == "sample" and g == "B"]
    for j in range(3):
        fit = lowess(neg_v[qc, j], x_all[qc], frac=0.75, it=0, delta=0.0, xvals=x_all)
        corrected = neg_v[:, j] / fit * np.median(neg_v[qc, j])
        t, p = sps.ttest_ind(np.log2(corrected[b]), np.log2(corrected[a]), equal_var=False)
        row = rows[f"negative:{neg.feature_ids[j]}"]
        assert row["p"] == pytest.approx(p, rel=1e-9) and row["statistic"] == pytest.approx(t, rel=1e-9)
        # 陽イオンの順で当てはめた（直す前の）p とは違う
        fit_pos = lowess(neg_v[qc, j], x_pos[qc], frac=0.75, it=0, delta=0.0, xvals=x_pos)
        wrong = neg_v[:, j] / fit_pos * np.median(neg_v[qc, j])
        wrong_p = sps.ttest_ind(np.log2(wrong[b]), np.log2(wrong[a]), equal_var=False).pvalue
        assert row["p"] != pytest.approx(wrong_p, rel=1e-6)

    # 統合が陰イオンの測定を残していない（測定を残す前の版の統合）なら、陰イオンの注入順が無いと極性を添えて断る
    # （陽イオンの順で陰イオンを当てはめない）
    merged.quant_state = {k: v for k, v in merged.quant_state.items() if k != "polarity_runs"}
    with pytest.raises(DomainError) as exc:
        A.run_spec(spec, merged)
    assert exc.value.code == "INJECTION_ORDER_MISSING" and exc.value.details["polarity"] == "negative"
    assert exc.value.message.startswith("negative: ")


def test_v1_recipe_drift_correct_with_tsv_order_notes_the_negative_run(tmp_path, monkeypatch):
    """injection_order が TSV（試料ごとに 1 つ）なら全極性にそれを使い、陰イオンの測定と違えば notes に書く
    （v1 の結果と dataset_analysis_run の戻り値）。同じなら書かない。"""
    import json

    from metabolomix.analysis import analysis_spec as A
    from metabolomix.core import session_state
    from metabolomix.quant import merge as qm
    from metabolomix.tools.reanalysis_tools import dataset_analysis_run
    _point(tmp_path, monkeypatch)
    pos, neg, _pv, _nv = _sheeted_pair(tmp_path, POS_TRUE, NEG_TRUE)
    merged, _info = qm.merge(pos, neg, polarities=("positive", "negative"), mapping_rows=None)
    tsv = tmp_path / "order.tsv"
    tsv.write_text("sample\torder\n" + "".join(f"{sid}\t{o}\n" for sid, o in zip(IDS, POS_TRUE)),
                   encoding="utf-8")
    raw = _v1_spec({"steps": ["drift_correct"]}, {"source": "tsv", "path": str(tsv)})
    notes: list[str] = []
    result = A.run_spec(A.normalize_spec(raw), merged, notes=notes)
    n_differ = sum(a != b for a, b in zip(POS_TRUE, NEG_TRUE))
    assert len(notes) == 1 and f"{n_differ} 試料で違う" in notes[0] and "negative" in notes[0]
    drift = next(s for s in result["analyses"][0]["matrix"]["steps"] if s["step"] == "drift_correct")
    assert drift["by_polarity"]["negative"]["injection_order"] == "tsv"
    session_state.session.dataset = merged
    out = json.loads(dataset_analysis_run(raw))
    assert out["status"] == "success" and any("drift_correct" in n for n in out["notes"])

    same_pos, same_neg, _a, _b = _sheeted_pair(tmp_path / "same", POS_TRUE, POS_TRUE)
    same, _ = qm.merge(same_pos, same_neg, polarities=("positive", "negative"), mapping_rows=None)
    quiet: list[str] = []
    A.run_spec(A.normalize_spec(raw), same, notes=quiet)
    assert quiet == []


# --- レビュー指摘: ある極性の係数の特徴をすべて除く指定は断らず、その極性の残った特徴すべてから --------------

def test_merged_exclusion_covering_a_polarity_falls_back_for_that_polarity():
    prep = P.preprocess(_merged(), {"drift": {"mode": "off"}, "filters": OFF_FILTERS,
                                    "normalize": {"method": "tic", "exclude_features": list(NEG_IDS)}})
    by = _step(prep, "normalize")["by_polarity"]
    assert by["negative"]["factor_features_fallback"] == "all_factor_features_excluded_in_polarity"
    assert by["negative"]["factor_features"] == 3 and by["negative"]["factor_feature_basis"] == "kept_features"
    assert "factor_features_fallback" not in by["positive"]
    note = next(w for w in prep.warnings if w["code"] == "normalization_exclude_features_cover_polarity")
    assert note["polarities"] == ["negative"]
    # 陰イオンの係数は陰イオンの特徴すべての TIC、陽イオンは陽イオンの特徴すべての TIC（numpy で）
    rows = [i for i, r in enumerate(ROLES) if r != "blank"]
    study = _merged()
    neg_tot, pos_tot = study.values[rows, 4:].sum(axis=1), study.values[rows, :4].sum(axis=1)
    np.testing.assert_allclose(by["negative"]["factors"], neg_tot / np.median(neg_tot), rtol=1e-12)
    np.testing.assert_allclose(by["positive"]["factors"], pos_tot / np.median(pos_tot), rtol=1e-12)
    # features と exclude で陰イオンが空になる組も同じ（陽イオンは features のとおり）
    both = P.preprocess(_merged(), {"drift": {"mode": "off"}, "filters": OFF_FILTERS,
                                    "normalize": {"method": "tic", "features": ["positive:p0", NEG_IDS[0]],
                                                  "exclude_features": [NEG_IDS[0]]}})
    by = _step(both, "normalize")["by_polarity"]
    assert by["negative"]["factor_features_fallback"] == "all_factor_features_excluded_in_polarity"
    assert by["positive"]["factor_features"] == 1
    # 全極性で空になる指定は、統合でないデータセットと同じく断る
    with pytest.raises(DomainError) as exc:
        P.preprocess(_merged(), {"normalize": {"method": "tic",
                                               "exclude_features": POS_IDS + NEG_IDS}})
    assert exc.value.code == "PREPROCESS_CONFIG_INVALID"


def test_merged_exclusion_emptying_a_polarity_without_listed_features_reports_the_exclusion():
    """features に陰イオンの特徴が無く、除外が陰イオンの残った特徴をすべて外す: 陰イオンは除外も当てずに残った特徴
    すべてから求め、理由は除外（記録が除外を当てたように読めない）。"""
    prep = P.preprocess(_merged(), {"drift": {"mode": "off"}, "filters": OFF_FILTERS,
                                    "normalize": {"method": "tic", "features": ["positive:p0"],
                                                  "exclude_features": list(NEG_IDS)}})
    by = _step(prep, "normalize")["by_polarity"]
    assert by["negative"]["factor_features_fallback"] == "all_factor_features_excluded_in_polarity"
    assert by["negative"]["factor_features"] == 3 and "excluded_from_factors" not in by["negative"]
    codes = {w["code"] for w in prep.warnings}
    assert "normalization_exclude_features_cover_polarity" in codes
