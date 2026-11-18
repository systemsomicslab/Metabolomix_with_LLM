"""最初の概観の事実の組み立て（metabolomix.analysis.overview）。"""
from __future__ import annotations

from metabolomix.analysis import overview, stats_results
from metabolomix.analysis.stats.composition import composition
from metabolomix.analysis.stats.pca import pca
from tests.test_stats_results import OFF, _ds, _run


def _record(ds, comparisons):
    n, prepared, rr = _run(ds, comparisons)
    rec = stats_results.register_compare(ds, rr, log2fc_threshold=1.0, units=prepared.units,
                                         make_current=False)
    return rec, prepared


def _prepared(ds):
    return _run(ds, [{"id": "x", "factor": "group", "reference": ds_groups(ds)[0],
                      "test": ds_groups(ds)[1]}])[1]


def ds_groups(ds):
    out = []
    for row in ds.sample_metadata_rows:
        if row.get("group") and row["group"] not in out:
            out.append(row["group"])
    return out


def _plan(ds):
    return overview.group_plan(ds, _prepared(ds))


def test_group_plan_two_groups_reference_is_first_seen():
    plan = _plan(_ds(groups=("T", "C")))
    assert plan["mode"] == "two" and plan["levels"] == ["T", "C"]
    assert plan["reference"] == "T" and plan["direction"] == "provisional"
    comps, omitted = overview.group_comparisons(plan)
    assert comps == [{"id": "C_vs_T", "factor": "group", "reference": "T", "test": "C"}] and not omitted


def test_group_plan_single_group_is_none():
    from metabolomix.analysis.stats.data import from_dataset
    from metabolomix.analysis.stats.preprocess import preprocess

    ds = _ds(groups=("C",))
    prepared = preprocess(from_dataset(ds), OFF)
    plan = overview.group_plan(ds, prepared)
    assert plan["mode"] == "none" and overview.group_comparisons(plan) == ([], False)


def test_group_comparisons_caps_pairs_at_six_groups():
    six = {"mode": "multi", "levels": list("ABCDEF"), "reference": "A"}
    comps, omitted = overview.group_comparisons(six)
    assert len(comps) == 15 and not omitted
    seven = {"mode": "multi", "levels": list("ABCDEFG"), "reference": "A"}
    comps, omitted = overview.group_comparisons(seven)
    assert len(comps) == 6 and omitted and all(c["reference"] == "A" for c in comps)


def test_species_facts_counts_verdicts_and_directions():
    ds = _ds()
    plan = _plan(ds)
    comps, _ = overview.group_comparisons(plan)
    rec, _ = _record(ds, comps)
    down = rec["engine"]["rows"][8]
    down["verdict"], down["log2fc"] = "changed", -2.0          # 下降を 1 つ仕込む
    sec = overview.species_facts(ds, record=rec, plan=plan, top_n=5, pairs_omitted=False,
                                 omitted=None, quality_note=None)
    first = sec["facts"][0]
    assert first["id"] == "S1" and first["data"]["up"] == 3 and first["data"]["down"] == 1
    assert "仮の向き" in first["text"]
    header = sec["top"].splitlines()[0].split("\t")
    assert header == ["contrast_id", "feature_id", "name", "id_level", "log2fc", "ci_low", "ci_high", "q"]
    assert len(sec["top"].splitlines()) == 5          # f0〜f2 と仕込んだ f8
    assert f'result_id="{rec["result_id"]}"' in sec["next"]


def test_species_facts_never_turns_not_significant_into_unchanged():
    ds = _ds()
    plan = _plan(ds)
    rec, _ = _record(ds, overview.group_comparisons(plan)[0])
    rows = rec["engine"]["rows"]
    rows[5]["verdict"] = "inconclusive"
    rows[6]["verdict"] = "unchanged"
    sec = overview.species_facts(ds, record=rec, plan=plan, top_n=5, pairs_omitted=False,
                                 omitted=None, quality_note=None)
    data = sec["facts"][0]["data"]
    assert data["unchanged"] == sum(r["verdict"] == "unchanged" for r in rows)
    assert data["inconclusive"] == sum(r["verdict"] == "inconclusive" for r in rows)


def test_species_facts_reports_presence_only_count():
    ds = _ds()
    plan = _plan(ds)
    rec, _ = _record(ds, overview.group_comparisons(plan)[0])
    rows = rec["engine"]["rows"]
    rows[7]["verdict"], rows[7]["presence_q"] = "inconclusive", 0.001
    rows[0]["presence_q"] = 0.001                      # changed の行は「出方だけ」に入らない
    sec = overview.species_facts(ds, record=rec, plan=plan, top_n=5, pairs_omitted=False,
                                 omitted=None, quality_note=None)
    pres = sec["facts"][1]["data"]["feature_ids"]
    assert "f7" in pres and "f0" not in pres
    rec["engine"]["presence"]["detection_available"] = False
    sec = overview.species_facts(ds, record=rec, plan=plan, top_n=5, pairs_omitted=False,
                                 omitted=None, quality_note=None)
    assert sec["facts"][1]["data"]["n"] is None


def test_multi_group_pair_table_counts_only_f_significant():
    ds = _ds(groups=("C", "T", "U"))
    plan = _plan(ds)
    assert plan["mode"] == "multi"
    comps, _ = overview.group_comparisons(plan)
    rec, _ = _record(ds, comps)
    sig = {o["feature_id"] for o in rec["engine"]["omnibus"] if o["f_q"] is not None and o["f_q"] <= 0.05}
    for r in rec["engine"]["rows"]:                     # F で差の無い特徴に changed を仕込む
        if r["feature_id"] not in sig:
            r["verdict"] = "changed"
    sec = overview.species_facts(ds, record=rec, plan=plan, top_n=5, pairs_omitted=False,
                                 omitted=None, quality_note=None)
    pairs = sec["facts"][1]["data"]["pairs"]
    sig_rows = [r for r in rec["engine"]["rows"] if r["feature_id"] in sig]
    assert sum(v[0] for v in pairs.values()) == sum(r["verdict"] == "changed" for r in sig_rows)
    assert sec["top"].splitlines()[0] == "feature_id\tname\tid_level\tf_q"


def test_class_facts_omitted_before_select():
    sec = overview.class_facts(composition_result=None, gate="付加体の選別の前",
                               gate_next="dataset_select_features()", plan={"mode": "two"}, quality_note=None)
    assert sec["omitted"] == "付加体の選別の前" and "dataset_select_features" in sec["next"]
    assert sec["facts"] == []


def _classed(ds):
    for k, fid in enumerate(ds.feature_ids):
        ds.feature_metadata[fid]["ontology"] = "PC" if k < 3 else ("PE" if k < 6 else "TG")
    return ds


def test_class_facts_lists_changed_and_unchanged():
    ds = _classed(_ds(per_group=6))
    plan = _plan(ds)
    comps, _ = overview.group_comparisons(plan)
    _, prepared = _record(ds, comps)
    comp = composition(prepared, assay_kind="lipid", factors=["group"], comparisons=comps,
                       test="moderated_t", effect_size="log2_mean_difference", bh_population="matrix",
                       q_threshold=0.05, equivalence_fc=1.2)
    sec = overview.class_facts(composition_result=comp, gate=None, gate_next="", plan=plan, quality_note=None)
    verdicts = [line.split("\t")[5] for line in sec["table"].splitlines()[1:]]
    assert "changed" in verdicts and "unchanged" in verdicts
    assert "脂質クラス" in sec["facts"][0]["text"] and "group_intensity" in sec["next"]


def _composition(ds, assay_kind="lipid"):
    plan = _plan(ds)
    comps, _ = overview.group_comparisons(plan)
    _, prepared = _record(ds, comps)
    comp = composition(prepared, assay_kind=assay_kind, factors=["group"], comparisons=comps,
                       test="moderated_t", effect_size="log2_mean_difference", bh_population="matrix",
                       q_threshold=0.05, equivalence_fc=1.2)
    return comp, plan


def test_class_facts_report_non_lipid_annotations_separately():
    ds = _classed(_ds(per_group=6))
    for fid, name in (("f6", "Ginsenoside Rb1"), ("f7", "no MS2: Saikosaponin A"), ("f8", "Enniatin B")):
        ds.feature_metadata[fid].update(name=name, ontology="Triterpene saponins" if fid != "f8"
                                        else "Cyclic depsipeptides")
    comp, plan = _composition(ds)
    sec = overview.class_facts(composition_result=comp, gate=None, gate_next="", plan=plan, quality_note=None)
    classes = [line.split("	")[0] for line in sec["table"].splitlines()[1:]]
    assert set(classes) <= {"PC", "PE", "TG"}
    last = sec["facts"][-1]
    assert "非脂質" in last["text"] and "3 個" in last["text"]
    assert last["data"]["n_features"] == 3
    assert last["data"]["by_ontology"] == {"Triterpene saponins": 2, "Cyclic depsipeptides": 1}
    # 一般代謝物では非脂質の事実を出さず、どのクラスも数える
    comp_m, plan_m = _composition(ds, assay_kind="metabolite")
    sec_m = overview.class_facts(composition_result=comp_m, gate=None, gate_next="", plan=plan_m,
                                 quality_note=None)
    assert not any("非脂質" in f["text"] for f in sec_m["facts"])
    assert set(comp_m["classes"]) == {"PC", "PE", "TG", "Triterpene saponins", "Cyclic depsipeptides"}


def test_class_facts_report_lab_msp_lipids_moved_to_the_msdial_class():
    ds = _classed(_ds(per_group=6))
    # TG の 1 つを研究室 MSP の行（ClassyFire のクラス・名前は ID）にする: 組成式から TG に救われる
    ds.feature_metadata["f6"].update(name="MassSpecGymID0000001", ontology="Triacylglycerols", formula="C57H104O6")
    comp, plan = _composition(ds)
    assert "f6" in comp["classes"]["TG"]
    sec = overview.class_facts(composition_result=comp, gate=None, gate_next="", plan=plan, quality_note=None)
    fact = next(f for f in sec["facts"] if "MS-DIAL の脂質クラスへ移して" in f["text"])
    assert fact["data"]["n_features"] == 1 and fact["data"]["n_rescued"] == 1
    assert fact["data"]["moves"] == [{"from": "Triacylglycerols", "to": "TG", "route": "class_formula",
                                      "n_features": 1}]
    assert "Triacylglycerols→TG 1" in fact["text"]
    # 移したものが無ければ事実を出さない
    plain = _classed(_ds(per_group=6))
    comp_p, plan_p = _composition(plain)
    sec_p = overview.class_facts(composition_result=comp_p, gate=None, gate_next="", plan=plan_p, quality_note=None)
    assert not any("移して" in f["text"] for f in sec_p["facts"])


def test_class_facts_omitted_names_the_non_lipid_count_when_no_lipid_class_is_left():
    ds = _ds(per_group=6)
    for fid in ds.feature_ids:
        ds.feature_metadata[fid].update(name="Ginsenoside Rb1", ontology="Triterpene saponins")
    comp, plan = _composition(ds)
    sec = overview.class_facts(composition_result=comp, gate=None, gate_next="n", plan=plan, quality_note=None)
    assert sec["facts"] == [] and "非脂質" in sec["omitted"] and str(len(ds.feature_ids)) in sec["omitted"]


def test_global_facts_eta_squared():
    ds = _ds()
    for j, name in enumerate(ds.sample_names):          # 全特徴が群で大きく分かれる
        if name.startswith("T"):
            ds.feature_matrix[:, j] *= 16
    plan = _plan(ds)
    comps, _ = overview.group_comparisons(plan)
    rec, prepared = _record(ds, comps)
    res = pca(prepared)
    labels = stats_results.unit_group_labels(ds, prepared)
    sec = overview.global_facts(pca_result=res, labels=labels, record=rec, plan=plan,
                                omitted=None, quality_note=None)
    assert sec["facts"][0]["data"]["eta_squared"][0] > 0.9
    assert 'dataset_compare(reference_group="C", test_group="T")' in sec["next"]


def test_fact_ids_are_sequential_per_stage():
    q = overview.quality_facts(
        rsd={"status": "ok", "median": 0.1, "frac_le_030": 0.9, "n_computable": 10, "n_qc": 4},
        failed={"checked": 4, "failed": ["id_QC1"]}, detection={"status": "not_evaluable", "reason": "x"},
        drift=[], confounding={"batch_status": "single_batch", "injection_order": {"status": "ok"}},
        replicate_warnings=[], verdict="poor", reasons=["failed_qc"], n_qc=4, n_blank=0,
        drift_note=None, pca_id="pca_x")
    assert [x["id"] for x in q["facts"]] == [f"Q{k}" for k in range(1, len(q["facts"]) + 1)]
    assert q["facts"][0]["text"].startswith("Q1 判定 poor")
    assert "dataset_exclude" in q["next"] and 'result_id="pca_x"' in q["next"]


def test_chain_state_provisional_until_select():
    ds = _ds()
    assert overview.chain_state(ds, curation_judgments=None)["provisional"] is True
    ds.quant_state = {"normalized": True, "steps": ["normalize", "select"]}
    st = overview.chain_state(ds, curation_judgments=3)
    assert st["provisional"] is False and st["curation"] == 3 and st["pending"] == []


def _found_blocks():
    return {"status": "found", "k": 3, "thresholds": {}, "skipped": [{"group": "U", "n": 8,
                                                                     "reason": "block_size_does_not_divide"}],
            "groups": [{"group": g, "n": 9, "order": "injection_order", "cells": "detected", "n_features": 200,
                        "layout": layout, "margin": 0.8, "purity": 1.0, "median_icc": 0.7, "adjacent_excess": 0.0,
                        "blocks": [[f"{g}{3 * b + j + 1:02d}" for j in range(3)] for b in range(3)]}
                       for g, layout in (("C", "consecutive"), ("T", "data"))]}


def test_quality_facts_lists_replicate_blocks_and_suggests_subject_id_once():
    """データからの技術反復らしいまとまりは、群ごとの候補を挙げた事実になり、next に subject_id の入れ方を 1 つ足す
    （名前からの反復の注意と両方出ても 1 つ）。データにはまとまりの戻り値そのもの。"""
    found = _found_blocks()
    base = dict(rsd={"status": "ok", "median": 0.1, "frac_le_030": 0.9, "n_computable": 10, "n_qc": 4},
                failed={"checked": 4, "failed": []}, detection={"status": "not_evaluable", "reason": "x"},
                drift=[], confounding={"batch_status": "single_batch", "injection_order": {"status": "ok"}},
                verdict="caution", reasons=["possible_replicate_blocks"], n_qc=4, n_blank=0, drift_note=None,
                pca_id=None)
    q = overview.quality_facts(replicate_warnings=[], replicate_blocks=found, **base)
    fact = next(x for x in q["facts"] if "技術反復らしいまとまり" in x["text"])
    assert "3 本ずつ" in fact["text"] and "2 群" in fact["text"]
    assert "C（注入順で連続）[C01, C02, C03] [C04, C05, C06] [C07, C08, C09]" in fact["text"]
    assert "T（相関から束ねた）[T01, T02, T03]" in fact["text"]
    assert "見ていない群" in fact["text"] and "U" in fact["text"]
    assert "自動では入れていない" in fact["text"]
    assert fact["data"]["replicate_blocks"] == found
    assert q["next"].count("dataset_set_sample_metadata") == 1 and "subject_id" in q["next"]
    both = overview.quality_facts(replicate_warnings=[{"stem": "S1", "samples": ["S1_1", "S1_2"]}],
                                  replicate_blocks=found, **base)
    assert both["next"].count("dataset_set_sample_metadata") == 1
    for quiet in (None, {"status": "not_found"}, {"status": "subject_id_set"}):
        q = overview.quality_facts(replicate_warnings=[], replicate_blocks=quiet, **base)
        assert not any("技術反復らしいまとまり" in x["text"] for x in q["facts"])
        assert "dataset_set_sample_metadata" not in q["next"]
