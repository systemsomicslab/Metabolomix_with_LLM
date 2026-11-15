"""最初の概観のツール `dataset_overview`（metabolomix.tools.overview_tools）。"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import session_state
from metabolomix.core.serialization import json_payload
from metabolomix.tools._holder import resolve_holder
from metabolomix.tools.overview_tools import build_overview_payload, dataset_overview
from metabolomix.tools.sample_tools import dataset_exclude
from tests.test_stats_results import _ds


def _install(state, tmp_path):
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


def _strict(text: str) -> dict:
    return json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))


def _ov(**kw) -> dict:
    return _strict(dataset_overview(**kw))


def _classed(ds):
    for k, fid in enumerate(ds.feature_ids):
        ds.feature_metadata[fid]["ontology"] = "PC" if k < 3 else ("PE" if k < 6 else "TG")
    ds.source_format = "msdial_folder"
    ds.quant_state = {"schema": "lipidmix.quant-state.v1", "normalized": False, "method": None,
                      "unit": None, "steps": ["select"], "skipped": ["normalize"]}
    return ds


def test_overview_two_groups_shape_and_ids(tmp_path):
    ds = _ds()
    ds.source_format = "msdial_folder"
    _install(ds, tmp_path)
    out = _ov()
    assert out["status"] == "success" and out["overview_id"].startswith("ovw_")
    assert out["provisional"] is True
    for key in ("quality", "global", "classes", "species"):
        assert key in out and "next" in out[key]
    assert out["quality"]["facts"][0].startswith("Q1 判定 ")
    assert out["species"]["facts"][0].startswith("S1 T_vs_C（仮の向き）")
    assert "付加体の選別の前" in out["classes"]["omitted"]
    assert "dataset_select_features" in out["classes"]["next"]
    assert out["result_ids"]["compare"] in ds.results and out["result_ids"]["pca"] in ds.results
    assert ds.results[out["overview_id"]]["kind"] == "overview"
    assert ds.last_differential is None and ds.last_pca is None
    assert "counts" in out and out["counts"]["n_features"] == 10
    assert "#S1" in out["note"]


def test_overview_is_deterministic(tmp_path):
    _install(_ds(), tmp_path)
    a, b = _ov(), _ov()
    for key in ("overview_id", "result_ids", "quality", "global", "classes", "species", "state"):
        assert a[key] == b[key]


def test_overview_single_group_skips_comparisons(tmp_path):
    _install(_ds(groups=("C",)), tmp_path)
    out = _ov()
    assert out["status"] == "success"
    assert "群が 1 つ以下" in out["species"]["omitted"]
    assert "compare" not in out["result_ids"] and "pca" in out["result_ids"]
    assert "dataset_set_sample_metadata" in out["global"]["next"]


def test_overview_complete_confounding_does_not_compare(tmp_path):
    ds = _ds()
    for row in ds.sample_metadata_rows:
        if row.get("group") == "T":
            row["batch"] = "B2"
    _install(ds, tmp_path)
    out = _ov()
    assert out["quality"]["verdict"] == "poor"
    assert any("完全交絡" in f for f in out["quality"]["facts"])
    assert "完全交絡" in out["species"]["omitted"]
    assert not any(r.get("kind") == "compare" for r in ds.results.values())


def test_overview_drift_held_by_auto_uses_no_correction_and_says_so(tmp_path):
    """換算済み・QC が試料の間に挟まる（補正の条件を満たす）とき、既定の drift=auto は補正しない。
    概観は補正なしで計算し、その旨と drift_mode="on" での呼び直しを事実に書く。"""
    ds = _ds(n_qc=6)
    qc_slots = [1, 5, 9, 13, 17, 22]
    others = iter(o for o in range(1, 23) if o not in qc_slots)
    slots = iter(qc_slots)
    for row in ds.sample_metadata_rows:
        row["injection_order"] = next(slots) if row["role"] == "qc" else next(others)
    ds.quant_state = {"schema": "lipidmix.quant-state.v1", "normalized": True,
                      "method": "internal_standard", "unit": "pmol", "steps": ["normalize"], "skipped": []}
    _install(ds, tmp_path)
    text = dataset_overview()
    out = _strict(text)
    assert "needs_input" not in out and out["status"] == "success"
    assert "drift skipped" in out["state"]
    assert 'dataset_preprocess(drift_mode="on")' in out["quality"]["next"]
    assert any("既定（drift=auto）では当てなかった（内部標準などで換算済み）" in f
               for f in out["quality"]["facts"])


def test_overview_drift_impossible_has_no_drift_note(tmp_path):
    """補正の条件を満たさない（QC が試料の前に固まっている）なら、drift=on を勧めない。"""
    ds = _ds(n_qc=6)
    ds.quant_state = {"schema": "lipidmix.quant-state.v1", "normalized": True,
                      "method": "internal_standard", "unit": "pmol", "steps": ["normalize"], "skipped": []}
    _install(ds, tmp_path)
    out = _ov()
    assert "drift_mode" not in out["quality"]["next"]
    assert not any("drift=auto" in f for f in out["quality"]["facts"])


def test_overview_multi_group_with_classes_and_size_budget(tmp_path):
    ds = _install(_classed(_ds(groups=("C", "T", "U"), per_group=5)), tmp_path)
    out = _ov()
    assert out["provisional"] is False
    assert "table" in out["classes"] and out["classes"]["table"].startswith("class\t")
    assert "pairs" in out["species"]
    assert out["global"]["facts"][-1].startswith("G3 moderated F")
    text = json_payload(build_overview_payload(session_state.session.current_handle(), embedded=True)) \
        if hasattr(session_state.session, "current_handle") else dataset_overview()
    assert len(text) <= 6000, len(text)
    assert ds.last_differential is None


def test_overview_strict_json(tmp_path):
    _install(_ds(), tmp_path)
    _strict(dataset_overview())


def test_overview_record_written_and_cites_preprocess(tmp_path):
    _install(_ds(), tmp_path)
    out = _ov()
    assert out["record"]["written"] is True
    records = [json.loads(p.read_text(encoding="utf-8")) for p in tmp_path.rglob("stats/*.json")]
    ov = [r for r in records if r["step"] == "overview"]
    pre = [r for r in records if r["step"] == "preprocess"]
    assert len(ov) == 1 and len(pre) == 1
    assert ov[0]["parent_ids"] == [pre[0]["record_id"]]
    assert ov[0]["inputs"]["result_id"] == out["overview_id"]


def test_overview_without_report_dir_still_succeeds():
    session_state.session.dataset = _ds()
    out = _ov()
    assert out["status"] == "success" and out["record"]["written"] is False


def test_overview_with_dataset_handle_uses_that_handle(tmp_path):
    first = _install(_ds(), tmp_path)
    first_id = session_state.session.current_handle_id
    second = _install(_ds(seed=9), tmp_path)
    assert session_state.session.current_handle_id != first_id
    out = _ov(dataset_handle=first_id)
    assert out["overview_id"] in first.results and out["overview_id"] not in second.results


def test_overview_after_exclusion_recomputes(tmp_path):
    ds = _install(_ds(), tmp_path)
    before = _ov()["overview_id"]
    assert "error" not in json.loads(dataset_exclude(samples=["T0_r0"], reason="test"))
    after = _ov()
    assert after["overview_id"] != before and after["overview_id"] in ds.results


def test_overview_bad_top_n_and_missing_dataset():
    """整数でない・負の top_n だけを断る（上限を超える値は下げて続ける）。"""
    for bad in (-1, True, "5", 1.5, None):
        assert _strict(dataset_overview(top_n=bad))["error"]["code"] == "DATASET_BAD_REQUEST", bad
    assert _strict(dataset_overview(top_n=99))["error"]["code"] == "missing_state"
    assert _strict(dataset_overview())["error"]["code"] == "missing_state"


def test_overview_exploratory_dataset_compares_and_carries_the_mark(tmp_path):
    """探索専用のデータセットも dataset_compare・書き出しと同じく比較を回し、比較の数字を使う段の先頭に
    印の事実（理由つき）を置く。探索専用でないデータセットの事実には印が無い。"""
    ds = _ds()
    ds.exploratory_only = True
    ds.exploratory_reason = "analysis job status = partial"   # 理由は 1 行（= は : に）
    _install(ds, tmp_path)
    out = _ov()
    assert "compare" in out["result_ids"] and out["result_ids"]["compare"] in ds.results
    assert "omitted" not in out["species"] and "omitted" not in out["global"]
    for key in ("global", "species"):
        first = out[key]["facts"][0]
        assert first.startswith(("G1 探索専用", "S1 探索専用")), first
        assert "analysis job status : partial" in first
    assert out["species"]["facts"][1].startswith("S2 T_vs_C（仮の向き）")
    record = ds.results[out["overview_id"]]
    assert record["settings"]["exploratory"] == {"exploratory_only": True,
                                                 "reason": "analysis job status : partial"}
    s1 = next(f for f in record["facts"]["species"] if f["id"] == "S1")
    assert s1["data"] == {"exploratory": True, "reason": "analysis job status : partial"}


def test_overview_of_a_complete_dataset_has_no_exploratory_mark(tmp_path):
    ds = _ds()
    _install(ds, tmp_path)
    out = _ov()
    assert not any("探索専用" in f for key in ("global", "species") for f in out[key]["facts"])
    assert "exploratory" not in ds.results[out["overview_id"]]["settings"]


def test_overview_embedded_omits_counts(tmp_path):
    _install(_ds(), tmp_path)
    assert "counts" not in build_overview_payload(resolve_holder(None), embedded=True)


def test_overview_non_folder_source_omits_classes(tmp_path):
    _install(_ds(), tmp_path)                 # source_format が解析フォルダでない（mzTab-M など）
    out = _ov()
    assert "定量の鎖" in out["classes"]["omitted"] and "composition" in out["classes"]["next"]


def test_overview_compare_stays_out_of_default_feature_search(tmp_path):
    """仮の向きの概観の比較は、特徴検索の既定（result_ids 省略）に混ざらない。名指しすれば引ける。"""
    from metabolomix.analysis.feature_search import select_results

    ds = _install(_ds(), tmp_path)
    out = _ov()
    rid = out["result_ids"]["compare"]
    assert ds.results[rid]["origin"] == "overview"
    assert select_results(ds, None)["results"] == []
    assert select_results(ds, [rid])["results"]


def test_overview_does_not_relabel_an_explicit_compare(tmp_path):
    """利用者が同じ比較を先に回していたら、その記録を概観の印で上書きしない（既定の特徴検索から消さない）。"""
    from metabolomix.analysis.feature_search import select_results
    from metabolomix.tools.stats_tools import dataset_compare

    ds = _install(_ds(), tmp_path)
    explicit = json.loads(dataset_compare(reference_group="C", test_group="T"))
    rid = explicit.get("result_id")
    assert rid in ds.results, explicit
    _ov()
    assert ds.results[rid].get("origin") is None
    assert select_results(ds, None)["results"]


def _no_features_left(ds):
    """QC の RSD を大きくして、既定の前処理の足切りで特徴が 1 つも残らないようにする。"""
    import numpy as np

    rng = np.random.default_rng(1)
    for j, name in enumerate(ds.sample_names):
        if name.startswith("QC"):
            ds.feature_matrix[:, j] *= rng.uniform(0.01, 100, ds.feature_matrix.shape[0])
    return ds


def test_overview_keeps_quality_when_preprocessing_leaves_nothing(tmp_path):
    """前処理で特徴が残らなくても、品質の段（前処理前の試料から計算する）は出し、他の段は理由つきで省く。"""
    _install(_no_features_left(_ds()), tmp_path)
    out = _ov()
    assert out["status"] == "success", out
    assert out["quality"]["verdict"] == "poor"
    assert any("RSD" in f for f in out["quality"]["facts"])
    for key in ("global", "classes", "species"):
        assert "PREPROCESS_PRECONDITION" in out[key]["omitted"], out[key]
    assert "qc_rsd" in out["global"]["omitted"]           # 足切りの内訳（details）を落とさない
    assert "pca" not in out["result_ids"] and "compare" not in out["result_ids"]


def test_overview_confounding_fact_carries_the_file_order_caveat(tmp_path):
    """A6: 注入順が MS-DIAL の既定（ファイルの順）かもしれないとき、概観の交絡の事実の「注入順が群で偏っている」に
    「ファイル名の並びを見ているだけかもしれない」を添え、事実のデータと next にも出す。"""
    from scipy import stats as sps

    from tests.test_injection_order_source import _file_order_ds

    ds = _install(_file_order_ds(tmp_path), tmp_path)
    out = _ov()
    assert "injection_order_biased" in out["quality"]["facts"][0]
    conf = next(f for f in out["quality"]["facts"] if "注入順が群で偏っている" in f)
    assert "ファイル名の並びを見ているだけかもしれない" in conf and "injection_order" in conf
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" in out["quality"]["next"]
    data = next(x["data"] for x in ds.results[out["overview_id"]]["facts"]["quality"]
                if "注入順が群で偏っている" in x["text"])
    assert data["injection_order"] == "biased" and data["injection_order_may_be_file_order"] is True
    # 偏りの判定そのものは同じ（個体＝subject の注入順の平均で、Kruskal–Wallis を独立に計算）
    by_subject: dict[str, list[float]] = {}
    for row in ds.sample_metadata_rows:
        if row["role"] == "sample":
            by_subject.setdefault(row["group"] + "|" + row["subject_id"], []).append(row["injection_order"])
    means = {g: [sum(v) / len(v) for k, v in by_subject.items() if k.startswith(g + "|")] for g in ("C", "T")}
    assert sps.kruskal(means["C"], means["T"]).pvalue < 0.05


def test_overview_confounding_fact_without_file_order_has_no_caveat(tmp_path):
    """注入順がファイルの順でない（手で入れた順）なら、偏りの事実に但し書きもデータの印も付けない。"""
    ds = _ds()
    for row in ds.sample_metadata_rows:      # 群ごとに固まった順（偏り）。出所は利用者のシート
        row["provenance"] = {"injection_order": {"value": row["injection_order"], "source": "user_manifest",
                                                 "confidence": "confirmed"}}
    ds = _install(ds, tmp_path)
    out = _ov()
    fact = next(x for x in ds.results[out["overview_id"]]["facts"]["quality"] if "batch_status" in x["data"])
    assert fact["data"]["injection_order"] == "biased" and "注入順が群で偏っている" in fact["text"]
    assert "ファイル名の並び" not in fact["text"]
    assert "injection_order_may_be_file_order" not in fact["data"]
    assert "dataset_set_sample_metadata" not in out["quality"]["next"]


def test_overview_drift_note_names_the_held_polarity_when_another_cannot_be_corrected():
    """陽陰を統合したデータセットで、陽イオンは注入順がファイルの順かもしれず drift=auto が止め、陰イオンは QC が
    足りず補正できない: 段の全体に rerun_with は無いが、概観は止めた極性を名指しし、drift_mode="on" を勧めない。"""
    from dataclasses import replace

    from metabolomix.analysis import overview
    from metabolomix.analysis.stats import preprocess as P
    from metabolomix.tools.overview_tools import _drift_note
    from tests.test_stats_preprocess_merged import OFF_FILTERS, _merged

    study = _merged()
    study.samples = [replace(s, injection_order_may_be_file_order=True) for s in study.samples]
    runs = [dict(r, batch="b1" if r["injection_order"] <= 8 else "b2") for r in study.polarity_runs["negative"]]
    study = replace(study, polarity_runs={"negative": runs})
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    step = next(s for s in prep.steps if s["step"] == "drift")
    assert "rerun_with" not in step and step["by_polarity"]["positive"].get("rerun_with")
    note, rerunnable = _drift_note(prep)
    assert note is not None and rerunnable is False
    assert "positive に当てなかった" in note and "注入順が MS-DIAL の既定" in note
    assert "negative" in note and "PREPROCESS_PRECONDITION" in note
    assert P.run_remedy("positive") in note
    with pytest.raises(Exception):          # 書いたとおり drift=on は陰イオンの条件で止まる
        P.preprocess(study, {"filters": OFF_FILTERS, "drift": {"mode": "on"}})
    facts = overview.quality_facts(rsd={"status": "not_evaluable"}, failed={"failed": []},
                                   detection={"status": "not_evaluable"}, drift=[], confounding={},
                                   replicate_warnings=[], verdict="caution", reasons=[], n_qc=6, n_blank=1,
                                   drift_note=note, drift_rerunnable=rerunnable, pca_id=None)
    assert any(x["text"].endswith(note) for x in facts["facts"])
    assert "drift_mode" not in facts["next"]
    # 全極性が保留（段の全体に rerun_with）なら、これまでどおり drift_mode="on" を勧める
    held_all = P.preprocess(_merged(flagged_neg=True), {"filters": OFF_FILTERS})
    note2, rerun2 = _drift_note(held_all)
    assert rerun2 is True and 'dataset_preprocess(drift_mode="on")' in note2


def _interleaved_file_order_ds(tmp_path, *, source="msdial_mddata"):
    """注入順が raw のファイル名の順位（出所 `source`。既定は .mddata＝MS-DIAL の既定かもしれない順）で、ファイル名の
    順では群が交互に並ぶ（注入順は群で偏らない）DatasetState。QC 6 本（1 バッチ）が試料の間に挟まる。"""
    ds = _ds(n_qc=6)
    qc = [n for n in ds.sample_names if n.startswith("QC")]
    c = [n for n in ds.sample_names if n.startswith("C")]
    t = [n for n in ds.sample_names if n.startswith("T")]
    seq = []
    for k, name in enumerate(x for pair in zip(c, t) for x in pair):
        if k % 3 == 0:
            seq.append(qc[k // 3])
        seq.append(name)
    assert sorted(seq) == sorted(ds.sample_names)
    position = {name: k for k, name in enumerate(seq, start=1)}
    ds.sample_raw_paths = [str(tmp_path / "raw" / f"run{position[n]:02d}.raw") for n in ds.sample_names]
    for name, row in zip(ds.sample_names, ds.sample_metadata_rows):
        row["injection_order"] = position[name]
        row["provenance"] = {"injection_order": {"value": position[name], "source": source,
                                                 "confidence": "declared"}}
    return ds


def _drift_fact(ds, out) -> dict:
    return next(x for x in ds.results[out["overview_id"]]["facts"]["quality"] if "batches" in x["data"])


def test_overview_drift_sign_fact_carries_the_file_order_caveat_when_the_order_is_not_biased(tmp_path):
    """注入順がファイルの順かもしれず、群では偏らない（交絡の事実に但し書きが付かない）ときも、同じ注入順で測る
    ドリフトの兆候の事実に「ρ はファイル名の並びとの相関かもしれない」を添え、データの印と next の入れ方を出す。"""
    import numpy as np
    from scipy import stats as sps

    from metabolomix.analysis.stats import data as stats_data

    ds = _install(_interleaved_file_order_ds(tmp_path), tmp_path)
    assert all(s.injection_order_may_be_file_order for s in stats_data.from_dataset(ds).samples)
    out = _ov()
    assert out["status"] == "success", out
    drift = _drift_fact(ds, out)
    assert "ファイル名の並びとの相関" in drift["text"] and "injection_order に入れる" in drift["text"]
    assert drift["data"]["injection_order_may_be_file_order"] is True
    assert drift["data"]["batches"][0]["injection_order_may_be_file_order"] is True
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" in out["quality"]["next"]
    # ρ は QC の総量（log2）と注入順の Spearman のまま（独立に計算）
    qc = [j for j, n in enumerate(ds.sample_names) if n.startswith("QC")]
    rho, p = sps.spearmanr([ds.sample_metadata_rows[j]["injection_order"] for j in qc],
                           np.log2(ds.feature_matrix[:, qc].sum(axis=0)))
    assert drift["data"]["batches"][0]["rho"] == pytest.approx(rho)
    assert drift["data"]["batches"][0]["p"] == pytest.approx(p)
    # 注入順は群で偏らない: 交絡の事実に但し書きは付かない（個体の注入順の平均で Kruskal–Wallis を独立に計算）
    conf = next(x for x in ds.results[out["overview_id"]]["facts"]["quality"] if "batch_status" in x["data"])
    assert conf["data"]["injection_order"] == "ok" and "ファイル名の並び" not in conf["text"]
    by_subject: dict[tuple, list[float]] = {}
    for row in ds.sample_metadata_rows:
        if row["role"] == "sample":
            by_subject.setdefault((row["group"], row["subject_id"]), []).append(row["injection_order"])
    means = {g: [float(np.mean(v)) for (grp, _s), v in by_subject.items() if grp == g] for g in ("C", "T")}
    assert sps.kruskal(means["C"], means["T"]).pvalue >= 0.05


def test_overview_drift_sign_caveat_reaches_the_quality_only_overview(tmp_path):
    """既定の前処理が失敗し品質の段だけを出す（交絡は評価しない）ときも、ドリフトの兆候の事実に但し書きが付く
    （印は ρ に使った QC の注入順から読み、交絡の判定には頼らない）。"""
    ds = _install(_no_features_left(_interleaved_file_order_ds(tmp_path)), tmp_path)
    out = _ov()
    assert "PREPROCESS_PRECONDITION" in out["global"]["omitted"]
    drift = _drift_fact(ds, out)
    assert "ファイル名の並びとの相関" in drift["text"]
    assert drift["data"]["injection_order_may_be_file_order"] is True
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" in out["quality"]["next"]


def test_overview_drift_sign_fact_with_a_sheet_order_has_no_caveat(tmp_path):
    """同じ並びでも、利用者のシートが入れた注入順（確かな順）なら但し書きもデータの印も付けない。"""
    ds = _install(_interleaved_file_order_ds(tmp_path, source="user_manifest"), tmp_path)
    out = _ov()
    drift = _drift_fact(ds, out)
    assert "ファイル名の並び" not in drift["text"]
    assert "injection_order_may_be_file_order" not in drift["data"]
    assert "injection_order_may_be_file_order" not in drift["data"]["batches"][0]
    assert "dataset_set_sample_metadata" not in out["quality"]["next"]


def test_overview_drift_sign_caveat_names_the_polarity_whose_order_may_be_file_order():
    """陽陰を統合したデータセット: 陽イオンの測定の注入順だけがファイルの順かもしれないなら、その行にだけ印が付き、
    ドリフトの兆候の事実の但し書きはその極性を名指しする。"""
    from dataclasses import replace

    from metabolomix.analysis import overview
    from metabolomix.analysis.stats import quality
    from tests.test_stats_preprocess_merged import _merged

    study = _merged()
    study.samples = [replace(s, injection_order_may_be_file_order=True) for s in study.samples]
    rows = quality.drift_signs(study)
    assert {r["polarity"]: "injection_order_may_be_file_order" in r for r in rows} == {
        "positive": True, "negative": False}
    facts = overview.quality_facts(rsd={"status": "not_evaluable"}, failed={"failed": []},
                                   detection={"status": "not_evaluable"}, drift=rows, confounding={},
                                   replicate_warnings=[], verdict="caution", reasons=[], n_qc=6, n_blank=1,
                                   drift_note=None, pca_id=None)
    drift = next(x for x in facts["facts"] if "batches" in x["data"])
    assert "（ただし positive の QC の注入順は" in drift["text"] and "negative の QC" not in drift["text"]
    assert drift["data"]["injection_order_may_be_file_order"] is True
    # 陽イオン（試料の表の極性）の測定は統合したハンドルへのシートで入る
    from metabolomix.analysis.stats.preprocess import run_remedy
    assert run_remedy("positive").rstrip("。") in drift["text"]
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" in facts["next"]
    assert "dataset_merge_polarities" not in facts["next"]


def _many_changed_ds(n_features=100, n_changed=24):
    """`_ds()` の試料で、特徴 n_features のうち先頭 n_changed が T 系の群で 4 倍（changed が 20 を超える）。変わる特徴は
    少数にし、どの特徴も生体の揺れを QC の揺れより大きくして（既定の D-ratio の足切りで変わらない特徴を落とさない）、
    既定の正規化（PQN）の係数を変わった特徴に引きずらせない。"""
    import numpy as np

    ds = _ds()
    rng = np.random.default_rng(3)
    spread = np.array([0.05 if name.startswith("QC") else 0.4 for name in ds.sample_names])[None, :]
    fm = np.exp2(rng.normal(20, 1, (n_features, 1)) + rng.normal(0, 1, (n_features, len(ds.sample_names))) * spread)
    for j, name in enumerate(ds.sample_names):
        if name.startswith("T"):
            fm[:n_changed, j] *= 4
    ds.feature_matrix = fm
    ds.detected_mask = np.ones_like(fm, dtype=bool)
    ds.feature_ids = [f"f{k}" for k in range(n_features)]
    ds.feature_metadata = {f: {"name": f.upper()} for f in ds.feature_ids}
    return ds


def test_overview_top_n_above_the_cap_is_capped_with_a_note(tmp_path):
    """top_n が上限（20）を超えても断らず 20 に下げて出し、notes に書く。上位の表を 20 行で切ったら
    species に top_truncated と載せなかった行数（top_n_more）。上限ちょうど・上限より小さい指定には印も注記も無い。"""
    ds = _install(_many_changed_ds(), tmp_path)
    out = _ov(top_n=30)
    assert out["status"] == "success", out
    assert any("top_n=30 は上限の 20 に下げた" in n for n in out["notes"])
    assert any("result_ids.compare" in n and "dataset_export_differential" in n for n in out["notes"])
    rows = out["species"]["top"].splitlines()[1:]
    assert len(rows) == 20
    record = ds.results[out["result_ids"]["compare"]]
    n_changed = sum(1 for r in record["engine"]["rows"] if r["verdict"] == "changed" and r["q"] is not None)
    assert n_changed > 20
    assert out["species"]["top_truncated"] is True and out["species"]["top_n_more"] == n_changed - 20
    # 下げた値は上限ちょうどの指定と同じ概観（設定・ID は 20 を記録する）
    at_cap = _ov(top_n=20)
    assert at_cap["overview_id"] == out["overview_id"]
    assert ds.results[out["overview_id"]]["settings"]["top_n"] == 20
    for plain in (at_cap, _ov(top_n=5)):
        assert "notes" not in plain and "top_truncated" not in plain["species"]
        assert "top_n_more" not in plain["species"]


def _merged_quality(study, *, drift_note=None, drift_rerunnable=True):
    from metabolomix.analysis import overview
    from metabolomix.analysis.stats import quality

    rows = quality.drift_signs(study)
    facts = overview.quality_facts(rsd={"status": "not_evaluable"}, failed={"failed": []},
                                   detection={"status": "not_evaluable"}, drift=rows, confounding={},
                                   replicate_warnings=[], verdict="caution", reasons=[], n_qc=6, n_blank=1,
                                   drift_note=drift_note, drift_rerunnable=drift_rerunnable, pca_id=None)
    return rows, facts, next(x for x in facts["facts"] if "batches" in x["data"])


def test_overview_names_the_negative_handle_remedy_when_the_negative_order_may_be_file_order():
    """陽陰を統合したデータセットで陰イオンの測定の注入順だけがファイルの順かもしれない: 統合したハンドルへの
    シートは陽イオンの測定しか変えないので、但し書き・next・ドリフト補正の注記は、陰イオンのハンドルにシートを
    当てて統合し直す手順（run_remedy("negative")）を示し、統合したハンドルへのシートを勧めない。"""
    from metabolomix.analysis.stats import preprocess as P
    from metabolomix.tools.overview_tools import _drift_note
    from tests.test_stats_preprocess_merged import OFF_FILTERS, _merged

    study = _merged(flagged_neg=True)
    prep = P.preprocess(study, {"filters": OFF_FILTERS})
    note, rerunnable = _drift_note(prep, study)
    rows, facts, drift = _merged_quality(study, drift_note=note, drift_rerunnable=rerunnable)
    assert {r["polarity"]: "injection_order_may_be_file_order" in r for r in rows} == {
        "positive": False, "negative": True}
    remedy = P.run_remedy("negative").rstrip("。")
    assert "（ただし negative の QC の注入順は" in drift["text"] and remedy in drift["text"]
    assert "試料の表の injection_order に入れる" not in drift["text"]
    assert "negative のハンドルに dataset_set_sample_metadata" in facts["next"]
    assert "dataset_merge_polarities" in facts["next"]
    # 統合し直すと統合したハンドルに当てたシートは引き継がれない（陽イオンのハンドルにも当てる）
    assert "positive のハンドルにも当てておく" in facts["next"]
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" not in facts["next"]
    # drift_mode="on" で呼び直せる（陽イオンは当たった）が、陰イオンは確かでない注入順のまま当てはめると添える
    assert rerunnable is True and 'dataset_preprocess(drift_mode="on")' in note
    assert "注入順が確かでない極性（negative）はその注入順のまま当てはめる" in note and remedy in note
    assert 'dataset_preprocess(drift_mode="on")' in facts["next"]


def test_overview_merged_remedy_when_both_polarity_orders_may_be_file_order():
    """両極性の注入順がファイルの順かもしれない: next は両極性のハンドルにシートを当てて統合し直す 1 手（統合し直すと
    統合したハンドルへのシートは消えるので、陽イオンも陽イオンのハンドルへ）。但し書きは陰イオンの入れ直し
    （統合し直す）を先に、統合したハンドルへの陽イオンのシートを後に書く。"""
    from dataclasses import replace

    from metabolomix.analysis.stats import preprocess as P
    from metabolomix.tools.overview_tools import _drift_note
    from tests.test_stats_preprocess_merged import OFF_FILTERS, _merged

    study = _merged(flagged_neg=True)
    study.samples = [replace(s, injection_order_may_be_file_order=True) for s in study.samples]
    note, rerunnable = _drift_note(P.preprocess(study, {"filters": OFF_FILTERS}), study)
    _rows, facts, drift = _merged_quality(study, drift_note=note, drift_rerunnable=rerunnable)
    assert "（ただし positive, negative の QC の注入順は" in drift["text"]
    neg, pos = P.run_remedy("negative"), P.run_remedy("positive").rstrip("。")
    assert drift["text"].index(neg) < drift["text"].index(pos)
    # drift=auto の注記も同じ順（両極性とも保留で drift_mode="on" は呼び直せる）
    assert rerunnable is True and "注入順が確かでない極性（positive, negative）" in note
    assert note.index(neg) < note.index(pos)
    assert "positive・negative のハンドルに dataset_set_sample_metadata" in facts["next"]
    assert "当てておく" not in facts["next"]
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" not in facts["next"]


def test_overview_merged_negative_remedy_clears_the_caveat_end_to_end(tmp_path, monkeypatch):
    """統合したハンドルの概観が示す手順が、陰イオンの注入順に実際に届く: 統合したハンドルへのシートの後も陰イオンの
    但し書きは残り（統合したハンドルへのシートは陰イオンを変えない）、概観の next どおりに陰イオン（と陽イオン）の
    ハンドルにシートを当てて統合し直すと、但し書き・印・入れ方の next が消え、ドリフト補正も両極性に当たる。"""
    from metabolomix.analysis.stats import preprocess as P
    from metabolomix.core import handles, session_state
    from metabolomix.quant import merge as qm
    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata
    from tests.test_stats_preprocess_merged import (FILE_ORDER, IDS, POS_TRUE, _folder, _point,
                                                    _sheet_text, _values)

    _point(tmp_path, monkeypatch)
    pos_v, neg_v = _values()
    ids = []
    for sub, polarity, values, ion in (("pos", "positive", pos_v, "PC"), ("neg", "negative", neg_v, "PE")):
        session_state.session.dataset = _folder(tmp_path / "run" / sub, FILE_ORDER, values, ion)
        holder = handles.registry.get(session_state.session.current_handle_id)
        holder.polarity = polarity
        ids.append(holder.handle_id)

    def merged_overview():
        merged = json.loads(qm.merge_tool(handles.registry.get(ids[0]), handles.registry.get(ids[1]),
                                          None, None))["dataset_handle"]
        return merged, _ov(dataset_handle=merged)

    sheet = tmp_path / "true.tsv"
    sheet.write_text(_sheet_text(IDS, POS_TRUE), encoding="utf-8")
    merged, _first = merged_overview()
    assert json.loads(dataset_set_sample_metadata(str(sheet), dataset_handle=merged))["status"] == "success"
    out = _ov(dataset_handle=merged)
    ds = handles.registry.get(merged).dataset
    drift = _drift_fact(ds, out)
    assert {b["polarity"]: "injection_order_may_be_file_order" in b for b in drift["data"]["batches"]} == {
        "positive": False, "negative": True}
    assert P.run_remedy("negative").rstrip("。") in drift["text"]
    assert "negative のハンドルに dataset_set_sample_metadata" in out["quality"]["next"]
    assert "dataset_set_sample_metadata（injection_order に実際の測定順を入れる）" not in out["quality"]["next"]

    # next のとおり: 陰イオンのハンドルにシート（統合したハンドルで入れた陽イオンの順は陽イオンのハンドルにも）→ 統合し直す
    for handle in ids:
        assert json.loads(dataset_set_sample_metadata(str(sheet), dataset_handle=handle))["status"] == "success"
    merged2, out2 = merged_overview()
    drift2 = _drift_fact(handles.registry.get(merged2).dataset, out2)
    assert "injection_order_may_be_file_order" not in drift2["data"]
    assert not any("injection_order_may_be_file_order" in b for b in drift2["data"]["batches"])
    assert "ファイル名の並び" not in drift2["text"]
    assert "dataset_set_sample_metadata" not in out2["quality"]["next"]
    assert "dataset_merge_polarities" not in out2["quality"]["next"]
    assert not any("当てなかった" in f for f in out2["quality"]["facts"])


def _qc_bracketed_file_order_ds(tmp_path):
    """`_interleaved_file_order_ds`（注入順＝raw のファイル名の順位、出所は .mddata）の最後の QC を走りの最後に置く:
    QC が生体を挟むのでドリフト補正の条件を満たし、drift=auto は注入順が確かでないことで補正を止める。"""
    ds = _interleaved_file_order_ds(tmp_path)
    rows = dict(zip(ds.sample_names, ds.sample_metadata_rows))
    last = [n for n in ds.sample_names if n.startswith("QC")][-1]
    moved = rows[last]["injection_order"]
    for name, row in rows.items():
        pos = row["injection_order"]
        new = len(rows) if name == last else (pos - 1 if pos > moved else pos)
        row["injection_order"] = new
        row["provenance"]["injection_order"]["value"] = new
    ds.sample_raw_paths = [str(tmp_path / "raw" / f"run{rows[n]['injection_order']:02d}.raw") for n in ds.sample_names]
    return ds


def _sheet_declared_filename_order(ds):
    """注入順の値はそのまま、試料の表が宣言 `# injection_order_source = filename` で書いた出所の記録にする。"""
    for row in ds.sample_metadata_rows:
        row["provenance"] = {"injection_order": {"value": row["injection_order"], "source": "user_manifest",
                                                 "derived_from": "filename", "confidence": "declared"}}
    return ds


def test_overview_drift_held_note_names_the_order_source(tmp_path):
    """drift=auto が注入順の不確かさで補正を止めた注記は、出所どおりに言う: MS-DIAL が振った順なら「MS-DIAL の既定
    （ファイルの順）かもしれない」、シートの宣言でファイル名の順なら「シートの宣言でファイル名の順」（MS-DIAL の
    既定とは言わない）。"""
    from metabolomix.analysis.stats import data as stats_data

    def held_note():
        out = _ov()
        notes = [f for f in out["quality"]["facts"] if "当てなかった" in f]
        assert len(notes) == 1, out["quality"]["facts"]
        return notes[0]

    ds = _install(_qc_bracketed_file_order_ds(tmp_path), tmp_path)
    assert all(s.injection_order_may_be_file_order for s in stats_data.from_dataset(ds).samples)
    note = held_note()
    assert "注入順が MS-DIAL の既定（ファイルの順）かもしれない" in note and "シートの宣言" not in note

    ds = _install(_sheet_declared_filename_order(_qc_bracketed_file_order_ds(tmp_path / "declared")), tmp_path)
    assert all(s.injection_order_source == "filename" and not s.injection_order_may_be_file_order
               for s in stats_data.from_dataset(ds).samples)
    note = held_note()
    assert "シートの宣言（# injection_order_source = filename）でファイル名の順" in note
    assert "MS-DIAL の既定" not in note


def test_overview_drift_held_note_without_the_sample_table_names_both_sources():
    """試料の表を渡さない（出所を読めない）ときは、両方の出所を挙げる言い方（MS-DIAL の既定とだけは言わない）。"""
    from metabolomix.analysis.stats import preprocess as P
    from metabolomix.tools.overview_tools import _drift_note
    from tests.test_stats_preprocess_merged import OFF_FILTERS, _merged

    note, _ = _drift_note(P.preprocess(_merged(flagged_neg=True), {"filters": OFF_FILTERS}))
    assert "MS-DIAL の既定（ファイルの順）かシートの宣言でファイル名の順かもしれない" in note


def test_overview_cap_note_points_to_a_recovery_path_that_exists(tmp_path):
    """上限に下げた注記の全行の引き方は比べ方に合わせる（最終の確かめ）。比較を回さなかった概観（群が 1 つ）は
    result_ids.compare を案内しない（無い ID を探させない）。多群の上位は F の q の順なので、対比ごとの書き出しでなく
    dataset_compare(output_path=...) の F 検定の表（.ftests）を案内する。"""
    _install(_ds(groups=("C",)), tmp_path)
    single = _ov(top_n=30)
    assert "compare" not in single["result_ids"]
    assert single["notes"] == ["top_n=30 は上限の 20 に下げた（戻り値を肥大させない）。"]

    _install(_ds(groups=("C", "T", "U"), per_group=5), tmp_path)
    multi = _ov(top_n=30)
    assert "compare" in multi["result_ids"]
    (note,) = multi["notes"]
    assert "ftests" in note and "f_q" in note and "dataset_compare(output_path=...)" in note
    assert "dataset_export_differential" not in note
