# tests/test_pathway_handoff.py
import math

import numpy as np
import pytest

from metabolomix.pathway import handoff


def _arf_inputs(**over):
    # 4 サンプル（C×2, P×2）＋ QC 1。prop は height と area の 2 種。
    samples = ["C_1", "C_2", "P_1", "P_2", "QC_1"]
    feats = ["Spot_1_height", "Spot_2_height", "Spot_3_height", "Spot_4_height",
             "Spot_5_height", "Spot_1_area"]
    m = np.array([
        # Spot1 Spot2 Spot3 Spot4 Spot5 Spot1area
        [10.0, 1.0, 5.0, 7.0, 3.0, 999.0],
        [12.0, 1.0, 5.0, 7.0, 3.0, 999.0],
        [20.0, 2.0, 6.0, 7.0, np.nan, 999.0],
        [22.0, 2.0, 6.0, 7.0, 4.0, 999.0],
        [15.0, 1.5, 5.5, 7.0, 3.5, 999.0],
    ])
    meta = {s: {"role": "sample", "group": s.split("_")[0]} for s in samples}
    meta["QC_1"] = {"role": "qc", "group": "QC"}
    catalog = {
        1: {"Name": "PC 16:0_18:1", "Ontology": "PC"},
        2: {"Name": "PC 16:0_18:1", "Ontology": "PC"},   # 同名（強度が低い → 落ちる）
        3: {"Name": "Unknown", "Ontology": "Unknown"},
        4: {"Name": "PE 18:0_20:4", "Ontology": "PE"},
        5: {"Name": "Cer 18:1;O2/24:1", "Ontology": "Cer_NS"},
    }
    kw = dict(matrix=m, sample_names=samples, feature_names=feats, sample_meta=meta,
              excluded_samples=set(), excluded_spots=set(), catalog=catalog,
              props=["height", "area"], group_a="C", group_b="P",
              compared={"group_a": ["C_1", "C_2"], "group_b": ["P_1", "P_2"]})
    kw.update(over)
    return kw


def _arf_inputs_p3(**over):
    # _arf_inputs に P_3 を足したもの（群 P を 3 例にして、1 例を除外・比較外にできるようにする）。
    kw = _arf_inputs()
    samples = ["C_1", "C_2", "P_1", "P_2", "P_3", "QC_1"]
    m = np.insert(kw["matrix"], 4, kw["matrix"][3], axis=0)
    meta = dict(kw["sample_meta"])
    meta["P_3"] = {"role": "sample", "group": "P"}
    kw.update(matrix=m, sample_names=samples, sample_meta=meta)
    kw.update(over)
    return kw


def test_arf_keeps_only_biological_samples_with_comparison_labels():
    t = handoff.build_from_arf(**_arf_inputs())
    assert t.sample_names == ["C_1", "C_2", "P_1", "P_2"]
    assert t.groups == ["C", "C", "P", "P"]


def test_arf_pooled_groups_are_labelled_by_the_spec_string():
    kw = _arf_inputs_p3(group_a="ctl", group_b="trt",
                        compared={"group_a": ["C_1", "C_2"], "group_b": ["P_1", "P_2"]})
    t = handoff.build_from_arf(**kw)
    # 比較外の生体試料 P_3 は自分の群 P のまま残る
    assert dict(zip(t.sample_names, t.groups)) == {
        "C_1": "ctl", "C_2": "ctl", "P_1": "trt", "P_2": "trt", "P_3": "P"}


def test_arf_uses_first_prop_only_and_drops_unknown_and_dedups_by_name():
    t = handoff.build_from_arf(**_arf_inputs())
    assert t.names == ["PC 16:0_18:1", "PE 18:0_20:4", "Cer 18:1;O2/24:1"]
    assert t.ontologies == ["PC", "PE", "Cer_NS"]
    assert t.matrix.shape == (3, 4)
    assert t.matrix[0].tolist() == [10.0, 12.0, 20.0, 22.0]   # Spot1(height) が残る
    assert t.report["dropped_unannotated"] == 1
    assert t.report["dropped_duplicate_names"] == 1
    assert t.report["duplicate_name_examples"] == ["PC 16:0_18:1"]
    assert t.report["prop"] == "height"


def test_arf_excluded_spots_and_samples_are_not_written():
    t = handoff.build_from_arf(**_arf_inputs_p3(
        excluded_spots={4}, excluded_samples={"P_3"},
        compared={"group_a": ["C_1", "C_2"], "group_b": ["P_1", "P_2"]}))
    assert "PE 18:0_20:4" not in t.names
    assert "P_3" not in t.sample_names
    assert t.report["dropped_excluded_spots"] == 1


def test_arf_keeps_nan_as_missing():
    t = handoff.build_from_arf(**_arf_inputs())
    cer = t.names.index("Cer 18:1;O2/24:1")
    assert math.isnan(t.matrix[cer, 2])


def test_arf_too_small_group_is_bad_request():
    kw = _arf_inputs(group_a="ctl", group_b="trt",
                     compared={"group_a": ["C_1"], "group_b": ["P_1", "P_2"]})
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_arf(**kw)
    assert exc.value.code == "PATHWAY_BAD_REQUEST"


def test_uncompared_sample_whose_group_collides_with_a_label_is_not_written():
    # C_2 は比較外だが自分の群名 "C" が group_a の指定名と同じ。比較群に紛れ込ませると
    # 群 C が 2 例に見えて検定が通ってしまう → 書かない（結果 C=1 で BAD_REQUEST）。
    kw = _arf_inputs(compared={"group_a": ["C_1"], "group_b": ["P_1", "P_2"]})
    with pytest.raises(handoff.HandoffError):
        handoff.build_from_arf(**kw)


def _ds_inputs(**over):
    kw = dict(
        pp_matrix=np.array([[1.0, 50.0], [2.0, 60.0], [3.0, 70.0], [4.0, 80.0], [9.0, 9.0]]),
        pp_sample_names=["a1", "a2", "b1", "b2", "qc"],
        pp_feature_names=["11", "12"],
        roles={"a1": "sample", "a2": "sample", "b1": "sample", "b2": "sample", "qc": "qc"},
        sample_names=["a1", "a2", "b1", "b2", "qc"],
        sample_ids=[],
        metadata_rows=None,
        annotations={"11": {"name": "LPC 18:1", "ontology": "LPC", "ambiguous": False},
                     "12": {"name": None, "ambiguous": True}},
        group_a=["a1", "a2"], group_b=["b1", "b2"],
        group_a_label="ctl", group_b_label="trt")
    kw.update(over)
    return kw


def test_dataset_uses_sml_ontology_and_labels():
    t = handoff.build_from_dataset(**_ds_inputs())
    assert t.sample_names == ["a1", "a2", "b1", "b2"]
    assert t.groups == ["ctl", "ctl", "trt", "trt"]
    assert t.names == ["LPC 18:1"]
    assert t.matrix.tolist() == [[1.0, 2.0, 3.0, 4.0]]
    assert t.report["dropped_unannotated"] == 1


def test_dataset_metadata_groups_add_uncompared_samples():
    names = ["a1", "a2", "b1", "b2", "z1", "x1", "qc"]
    rows = [{"sample_id": n, "group": g, "role": "sample", "include": inc}
            for n, g, inc in [("a1", "A", True), ("a2", "A", True), ("b1", "B", True),
                              ("b2", "B", True), ("z1", "Z", True), ("x1", "X", False)]]
    kw = _ds_inputs(
        pp_matrix=np.arange(14, dtype=float).reshape(7, 2) + 1,
        pp_sample_names=names, sample_names=names,
        roles={**{n: "sample" for n in names}, "qc": "qc"},
        metadata_rows=rows)
    t = handoff.build_from_dataset(**kw)
    # z1 は実験情報シートの群 Z で載る。x1 は include=False、qc は役割で落ちる。
    assert t.sample_names == ["a1", "a2", "b1", "b2", "z1"]
    assert t.groups == ["ctl", "ctl", "trt", "trt", "Z"]


def test_dataset_too_small_group():
    with pytest.raises(handoff.HandoffError):
        handoff.build_from_dataset(**_ds_inputs(group_b=["b1"]))


def test_write_mztab_matches_upstream_minimal_shape(tmp_path):
    t = handoff.build_from_arf(**_arf_inputs())
    out = tmp_path / "input.mztab"
    handoff.write_mztab(t, out, normalization="pqn")
    lines = out.read_text(encoding="utf-8").splitlines()
    mtd = [l.split("\t") for l in lines if l.startswith("MTD")]
    assert ["MTD", "mzTab-version", "2.0.0-M"] in mtd
    assert ["MTD", "comment[1]", "Normalization: pqn"] in mtd
    assert ["MTD", "assay[1]", "C_1"] in mtd
    assert ["MTD", "study_variable[1]", "C"] in mtd
    assert ["MTD", "study_variable[1]-assay_refs", "assay[1]|assay[2]"] in mtd
    assert ["MTD", "study_variable[2]-assay_refs", "assay[3]|assay[4]"] in mtd
    smh = next(l.split("\t") for l in lines if l.startswith("SMH"))
    assert smh[:9] == ["SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
                       "chemical_formula", "smiles", "inchi", "chemical_name", "uri"]
    assert smh[-1] == "opt_global_Ontology"
    assert smh[14:18] == [f"abundance_assay[{i}]" for i in range(1, 5)]
    sml = [l.split("\t") for l in lines if l.startswith("SML")]
    assert sml[0][7] == "PC 16:0_18:1" and sml[0][-1] == "PC"
    cer = next(r for r in sml if r[7] == "Cer 18:1;O2/24:1")
    assert cer[16] == "null"   # P_1 の NaN
