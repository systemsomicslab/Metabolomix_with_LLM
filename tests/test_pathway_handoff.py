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
        pp_matrix=np.array([[1.0, 50.0, 5.0], [2.0, 60.0, 6.0], [3.0, 70.0, 7.0],
                            [4.0, 80.0, 8.0], [9.0, 9.0, 9.0]]),
        pp_sample_names=["a1", "a2", "b1", "b2", "qc"],
        pp_feature_names=["11", "12", "13"],
        roles={"a1": "sample", "a2": "sample", "b1": "sample", "b2": "sample", "qc": "qc"},
        sample_names=["a1", "a2", "b1", "b2", "qc"],
        sample_ids=[],
        metadata_rows=None,
        annotations={"11": {"name": "LPC 18:1", "ontology": "LPC", "ambiguous": False},
                     "12": {"name": None, "ambiguous": True},
                     "13": {"name": "PC 34:1", "ontology": "PC", "ambiguous": False}},
        group_a=["a1", "a2"], group_b=["b1", "b2"],
        group_a_label="ctl", group_b_label="trt")
    kw.update(over)
    return kw


def test_dataset_uses_sml_ontology_and_labels():
    t = handoff.build_from_dataset(**_ds_inputs())
    assert t.sample_names == ["a1", "a2", "b1", "b2"]
    assert t.groups == ["ctl", "ctl", "trt", "trt"]
    assert t.names == ["LPC 18:1", "PC 34:1"]
    assert t.matrix.tolist() == [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]
    assert t.report["dropped_unannotated"] == 1


def test_dataset_metadata_groups_add_uncompared_samples():
    names = ["a1", "a2", "b1", "b2", "z1", "x1", "qc"]
    rows = [{"sample_id": n, "group": g, "role": "sample", "include": inc}
            for n, g, inc in [("a1", "A", True), ("a2", "A", True), ("b1", "B", True),
                              ("b2", "B", True), ("z1", "Z", True), ("x1", "X", False)]]
    kw = _ds_inputs(
        pp_matrix=np.arange(21, dtype=float).reshape(7, 3) + 1,
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


def test_fewer_than_two_features_is_bad_request():
    # 上流の apply(...) は 1 特徴の表で壊れる → こちらで止める
    ann = {"11": {"name": "LPC 18:1", "ontology": "LPC", "ambiguous": False},
           "12": {"name": None, "ambiguous": True}, "13": {"name": "Unknown", "ontology": "Unknown"}}
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_dataset(**_ds_inputs(annotations=ann))
    assert exc.value.code == "PATHWAY_BAD_REQUEST"
    assert "2" in exc.value.message


def test_same_group_labels_are_bad_request():
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_dataset(**_ds_inputs(group_a_label="X", group_b_label="X"))
    assert exc.value.code == "PATHWAY_BAD_REQUEST"
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_arf(**_arf_inputs(group_a="C", group_b="C"))
    assert exc.value.code == "PATHWAY_BAD_REQUEST"


def test_write_mztab_replaces_hash_in_text_fields(tmp_path):
    # 上流の read.table は # をコメント開始と読む → 書く文字列すべてから # を消す
    t = handoff.HandoffTable(sample_names=["C#1", "C_2", "P_1", "P_2"],
                             groups=["C#", "C#", "P", "P"],
                             names=["PC 34:1#a", "PE 36:2"], ontologies=["PC#", "PE"],
                             matrix=np.ones((2, 4)))
    out = tmp_path / "input.mztab"
    handoff.write_mztab(t, out, normalization="pqn#x")
    text = out.read_text(encoding="utf-8")
    assert "#" not in text
    assert "PC 34:1_a" in text and "assay[1]	C_1" in text and "Normalization: pqn_x" in text


# ---------- 脂質でない行を除く（v8: MS-DIAL の脂質の Ontology ＋ 名前の検査） ----------

def _ds_with(annotations):
    n = len(annotations)
    ids = [str(100 + i) for i in range(n)]
    return _ds_inputs(
        pp_matrix=np.arange(5 * n, dtype=float).reshape(5, n) + 1,
        pp_feature_names=ids,
        annotations={fid: {**a, "ambiguous": False} for fid, a in zip(ids, annotations)})


_KEEP = [{"name": "PC 38:3", "ontology": "PC"}, {"name": "PE 18:0_20:4", "ontology": "PE"}]


def _dropped(report):
    return {k: v for k, v in report.items()
            if k.startswith("dropped_") and k not in ("dropped_duplicate_names", "dropped_excluded_spots")}


def test_handoff_keeps_msdial_lipids():
    t = handoff.build_from_dataset(**_ds_with(_KEEP))
    assert t.names == ["PC 38:3", "PE 18:0_20:4"] and t.ontologies == ["PC", "PE"]
    assert t.report["kept_by_route"] == {"as_is": 2, "name": 0, "class_formula": 0}
    assert not any(_dropped(t.report).values())
    # 限定子接頭辞付きの名前は読める（経路 0 は名前をそのまま通す）
    t = handoff.build_from_dataset(**_ds_with(
        _KEEP + [{"name": "no MS2: PC 34:1", "ontology": "PC"}]))
    assert t.names[-1] == "no MS2: PC 34:1"


def test_handoff_keeps_classyfire_rows_with_a_readable_name():
    # 研究室 MSP の Ontology は ClassyFire 系。名前が MS-DIAL の略記で読めれば経路 1 で残す。
    t = handoff.build_from_dataset(**_ds_with(
        _KEEP + [{"name": "FA 16:0", "ontology": "Long-chain fatty acids"},
                 {"name": "PC 34:1", "ontology": "Phosphatidylcholines"}]))
    assert t.names == ["PC 38:3", "PE 18:0_20:4", "FA 16:0", "PC 34:1"]
    assert t.ontologies == ["PC", "PE", "FA", "PC"]
    assert t.report["kept_by_route"] == {"as_is": 2, "name": 2, "class_formula": 0}
    assert not any(_dropped(t.report).values())


def test_handoff_report_counts_the_kept_features():
    """落とした行の件数と並べて、残った特徴の数（n_kept。重複名の整理の後）を返す。"""
    t = handoff.build_from_dataset(**_ds_with(
        _KEEP + [{"name": "Dibutyl phthalate", "ontology": "Benzoic acid esters"},
                 {"name": "PC 38:3", "ontology": "PC"}]))           # 同名は 1 行にまとめる
    assert t.report["n_kept"] == len(t.names) == 2
    assert t.report["dropped_unmapped_ontology"] == 1
    assert t.report["dropped_duplicate_names"] == 1


def test_handoff_drops_unparsable_names():
    t = handoff.build_from_dataset(**_ds_with(
        _KEEP + [{"name": "3-Sulfamoylbenzoic acid", "ontology": "PC"}]))
    assert t.names == ["PC 38:3", "PE 18:0_20:4"]
    assert t.report["dropped_unparsable_name"] == 1
    assert t.report["dropped_unmapped_ontology"] == 0


def test_arf_handoff_applies_the_same_lipid_checks():
    kw = _arf_inputs()
    kw["catalog"] = {**kw["catalog"],
                     2: {"Name": "Dibutyl phthalate", "Ontology": "Benzoic acid esters", "Formula": "C16H22O4"},
                     3: {"Name": "3-Sulfamoylbenzoic acid", "Ontology": "PC"}}
    t = handoff.build_from_arf(**kw)
    assert t.names == ["PC 16:0_18:1", "PE 18:0_20:4", "Cer 18:1;O2/24:1"]
    assert t.report["dropped_unmapped_ontology"] == 1
    assert t.report["dropped_unparsable_name"] == 1


def test_msp_rows_reach_the_handoff():
    """研究室 MSP の行（慣用名・ID・英語のクラス名 + ClassyFire + 組成式）から脂質が引ける（ユーザーの指示）。"""
    t = handoff.build_from_dataset(**_ds_with(_KEEP + [
        {"name": "Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"},
        {"name": "low score: MassSpecGymID0400493", "ontology": "Phosphatidylcholines", "chemical_formula": "C36H72NO8P"},
        {"name": "Phosphatidylethanolamine alkenyl 16:0-20:2", "ontology": "1-(1Z-alkenyl),2-acylglycerophosphoethanolamines",
         "chemical_formula": "C41H78NO7P"},
    ]))
    assert t.names == ["PC 38:3", "PE 18:0_20:4", "TG 54:3", "low score: PC 28:0", "PE P-16:0_20:2"]
    assert t.ontologies == ["PC", "PE", "TG", "PC", "EtherPE"]
    assert t.report["kept_by_route"] == {"as_is": 2, "name": 1, "class_formula": 2}


def test_name_wins_over_classyfire():
    t = handoff.build_from_dataset(**_ds_with(_KEEP + [
        {"name": "PE 36:7", "ontology": "1-(1Z-alkenyl),2-acyl-glycerophosphocholines", "chemical_formula": "C41H68NO8P"}]))
    assert t.names[-1] == "PE 36:7" and t.ontologies[-1] == "PE"


def test_kept_name_never_starts_with_a_space():
    t = handoff.build_from_dataset(**_ds_with(_KEEP + [
        {"name": "  low score: Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"}]))
    assert t.names[-1] == "low score: TG 54:3"


def test_drop_reasons_partition_the_rows():
    rows = _KEEP + [
        {"name": None, "ontology": "Triacylglycerols"},                                            # unannotated
        {"name": "PE 32:1", "ontology": "null"},                                                  # unannotated（経路 1 に回さない）
        {"name": "3-Sulfamoylbenzoic acid", "ontology": "PC"},                                    # unparsable_name
        {"name": "Cer 34:1;O2", "ontology": "Ceramides", "chemical_formula": "C34H67NO3"},        # ambiguous_class
        {"name": "Dibutyl phthalate", "ontology": "Benzoic acid esters", "chemical_formula": "C16H22O4"},   # unmapped_ontology
        {"name": "N-Palmitoyl-D-sphingosine", "ontology": "Long-chain ceramides", "chemical_formula": "C34H67NO3"},  # ambiguous_ontology
        {"name": "Glyceryl trimyristate", "ontology": "Triacylglycerols"},                       # no_formula
        {"name": "1-Palmitoyl-2-linoleoyl-rac-glycerol", "ontology": "Lineolic acids and derivatives",
         "chemical_formula": "C37H68O5"},                                                         # formula_mismatch
        {"name": "Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"},
        {"name": "TG 54:3", "ontology": "TG"},                                                    # 同じ名前 → 重複
    ]
    t = handoff.build_from_dataset(**_ds_with(rows))
    r = t.report
    dropped = {k: v for k, v in r.items() if k.startswith("dropped_") and k != "dropped_duplicate_names"}
    assert dropped == {"dropped_unannotated": 2, "dropped_unparsable_name": 1, "dropped_ambiguous_class": 1,
                       "dropped_unmapped_ontology": 1, "dropped_ambiguous_ontology": 1, "dropped_no_formula": 1,
                       "dropped_formula_mismatch": 1}
    assert sum(r["kept_by_route"].values()) + sum(dropped.values()) == len(rows)
    assert r["n_kept"] == sum(r["kept_by_route"].values()) - r["dropped_duplicate_names"] == 3
    assert "dropped_non_msdial_ontology" not in r


def test_arf_reads_formula_for_msp_rows():
    kw = _arf_inputs()
    # 鍵 3 は Unknown の行（Spot_3_height の列がある）。組成式は ARF の Formula 列から読む。
    kw["catalog"] = {**kw["catalog"], 3: {"Name": "Glyceryl trioleate", "Ontology": "Triacylglycerols",
                                          "Formula": "C57H104O6"}}
    t = handoff.build_from_arf(**kw)
    assert "TG 54:3" in t.names and t.report["kept_by_route"]["class_formula"] == 1
    assert t.report["kept_by_route"]["as_is"] == 4 and "dropped_non_msdial_ontology" not in t.report


def test_msdial_lipid_class_table_is_the_lbmclass_enum():
    from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES
    assert {"PC", "PE", "TG", "CE", "FA", "SM", "Cer_NS", "EtherPE", "LPC"} <= MSDIAL_LIPID_CLASSES
    assert not {"Undefined", "Others", "Unknown", "Phosphatidylcholines"} & MSDIAL_LIPID_CLASSES
