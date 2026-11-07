# tests/test_pathway_handoff.py
import math

import numpy as np
import pytest

from metabolomix.pathway import handoff




















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
    matrix = np.array([[1.0, 50.0, 5.0], [2.0, 60.0, 6.0], [3.0, 70.0, np.nan],
                       [4.0, 80.0, 8.0], [9.0, 9.0, 9.0]])
    t = handoff.build_from_dataset(**_ds_inputs(pp_matrix=matrix))
    out = tmp_path / "input.mztab"
    handoff.write_mztab(t, out, normalization="pqn")
    lines = out.read_text(encoding="utf-8").splitlines()
    mtd = [l.split("	") for l in lines if l.startswith("MTD")]
    assert ["MTD", "mzTab-version", "2.0.0-M"] in mtd
    assert ["MTD", "comment[1]", "Normalization: pqn"] in mtd
    assert ["MTD", "assay[1]", "a1"] in mtd
    assert ["MTD", "study_variable[1]", "ctl"] in mtd
    assert ["MTD", "study_variable[1]-assay_refs", "assay[1]|assay[2]"] in mtd
    assert ["MTD", "study_variable[2]-assay_refs", "assay[3]|assay[4]"] in mtd
    smh = next(l.split("	") for l in lines if l.startswith("SMH"))
    assert smh[:9] == ["SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
                       "chemical_formula", "smiles", "inchi", "chemical_name", "uri"]
    assert smh[-1] == "opt_global_Ontology"
    assert smh[14:18] == [f"abundance_assay[{i}]" for i in range(1, 5)]
    sml = [l.split("	") for l in lines if l.startswith("SML")]
    lpc = next(r for r in sml if r[7] == "LPC 18:1")
    assert lpc[-1] == "LPC"
    pc = next(r for r in sml if r[7] == "PC 34:1")
    assert pc[16] == "null"   # b1 の NaN


def test_fewer_than_two_features_is_bad_request():
    # 上流の apply(...) は 1 特徴の表で壊れる → こちらで止める
    ann = {"11": {"name": "LPC 18:1", "ontology": "LPC", "ambiguous": False},
           "12": {"name": None, "ambiguous": True}, "13": {"name": "Unknown", "ontology": "Unknown"}}
    with pytest.raises(handoff.HandoffError) as exc:
        handoff.build_from_dataset(**_ds_inputs(annotations=ann))
    assert exc.value.code == "PATHWAY_BAD_REQUEST"
    assert "2" in exc.value.message




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


_DHA_PE_P = "1-(1Z-Octadecenyl)-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-sn-glycero-3-phosphoethanolamine"


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
        {"name": "no MS2: DHA-[21,21,22,22,22-d5]", "ontology": "Very long-chain fatty acids",
         "chemical_formula": "C22H32O2"},                                                         # isotope_label
        {"name": _DHA_PE_P, "ontology": "1-(1Z-alkenyl),2-acyl-glycerophosphocholines",
         "chemical_formula": "C45H78NO7P"},                                                       # class_name_conflict
        {"name": "Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"},
        {"name": "TG 54:3", "ontology": "TG"},                                                    # 同じ名前 → 重複
    ]
    t = handoff.build_from_dataset(**_ds_with(rows))
    r = t.report
    dropped = {k: v for k, v in r.items() if k.startswith("dropped_") and k != "dropped_duplicate_names"}
    assert dropped == {"dropped_unannotated": 2, "dropped_unparsable_name": 1, "dropped_isotope_label": 1,
                       "dropped_ambiguous_class": 1, "dropped_unmapped_ontology": 1, "dropped_ambiguous_ontology": 1,
                       "dropped_class_name_conflict": 1, "dropped_no_formula": 1, "dropped_formula_mismatch": 1}
    assert sum(r["kept_by_route"].values()) + sum(dropped.values()) == len(rows)
    assert r["n_kept"] == sum(r["kept_by_route"].values()) - r["dropped_duplicate_names"] == 3
    assert "dropped_non_msdial_ontology" not in r




def test_msdial_lipid_class_table_is_the_lbmclass_enum():
    from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES
    assert {"PC", "PE", "TG", "CE", "FA", "SM", "Cer_NS", "EtherPE", "LPC"} <= MSDIAL_LIPID_CLASSES
    assert not {"Undefined", "Others", "Unknown", "Phosphatidylcholines"} & MSDIAL_LIPID_CLASSES


def test_msp_rows_with_an_isotope_label_are_dropped():
    # 経路 1・2 の前に落とす（内部標準をクラスの和に入れない）。経路 0（as_is）は変えない。
    for name, onto, formula in (("no MS2: DHA-[21,21,22,22,22-d5]", "Very long-chain fatty acids", "C22H32O2"),
                                ("PC 15:0_18:1(d7)", "Phosphatidylcholines", "C41H73[2H]7NO8P"),
                                ("Palmitic acid-d31", "Long-chain fatty acids", "C16H32O2"),
                                ("13C16-Palmitic acid", "Long-chain fatty acids", "C16H32O2")):
        assert handoff.resolve_row(name, onto, formula).reason == "isotope_label", name
    assert handoff.resolve_row("PS 15:0_18:1(d7)", "PS", "C39H67[2H]7NO10P").route == "as_is"
    assert handoff.resolve_row("SM d18:1/16:0", "Sphingomyelins", "C39H79N2O6P").route == "name"
    assert handoff.resolve_row("Unknown", "Very long-chain fatty acids", "C22H32O2").reason == "unannotated"


def test_class_name_conflict_drops_the_row():
    res = handoff.resolve_row("low score: " + _DHA_PE_P, "1-(1Z-alkenyl),2-acyl-glycerophosphocholines", "C45H78NO7P")
    assert res.reason == "class_name_conflict"
    res = handoff.resolve_row("low score: " + _DHA_PE_P, "1-(1Z-alkenyl),2-acylglycerophosphoethanolamines",
                              "C45H78NO7P")
    assert (res.route, res.name, res.ontology) == ("class_formula", "low score: PE P-40:6", "EtherPE")


def test_every_reason_resolve_row_can_return_is_a_report_key():
    """`resolve_row` が返しうる理由の文字列はすべて report のキーにある（新しい理由で KeyError にしない）。"""
    import ast
    import inspect

    from metabolomix.pathway import msp_lipid
    reasons = set()
    for module in (handoff, msp_lipid):
        for node in ast.walk(ast.parse(inspect.getsource(module))):
            if isinstance(node, ast.FunctionDef) and (node.name == "resolve_row" or node.name.startswith("from_")):
                for sub in ast.walk(node):
                    values = []
                    if isinstance(sub, ast.Return) and sub.value is not None:
                        values.append(sub.value)
                    if isinstance(sub, ast.keyword) and sub.arg == "reason":
                        values.append(sub.value)
                    for v in values:
                        if isinstance(v, ast.Constant) and isinstance(v.value, str):
                            reasons.add(v.value)
                        elif isinstance(v, ast.Attribute) and isinstance(getattr(msp_lipid, v.attr, None), str):
                            reasons.add(getattr(msp_lipid, v.attr))
                        elif isinstance(v, ast.Name) and isinstance(getattr(module, v.id, None), str):
                            reasons.add(getattr(module, v.id))
    assert {"unannotated", "unparsable_name", "isotope_label", "ambiguous_class", "unmapped_ontology",
            "ambiguous_ontology", "class_name_conflict", "no_formula", "formula_mismatch"} <= reasons
    report = handoff._new_report("dataset")
    assert {f"dropped_{r}" for r in reasons} <= set(report), reasons - {k[8:] for k in report}


