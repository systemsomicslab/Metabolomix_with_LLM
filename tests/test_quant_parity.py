"""本物の規則（msdial-workflow の msdial_rules）と同じ値になることの確認。

`--quant-rules=<規則パッケージのフォルダ>` を渡したときだけ走る（渡さなければ skip。本物の規則は
msdial-workflow 側のテストが縛る）。数値（行・係数・期待値）は msdial-workflow の人工データのテスト
（`tests/workflow/test_tables.py` `test_any_others.py` `tests/normalizer/test_core.py`）から写したもの。
MCP の鎖（正規化 → 統合 → 選別）を `write_msdial_folder` の人工フォルダに対して本物の規則で回す。

特徴 ID は MCP では連番（0 = 元の ID 10、1 = 42、2 = 90、3 = 77）。GUI の ID（10 / 42 / 90 / 77）とは
違うので、各テストに対応を書く。意図して違うところ（gap-fill の内部標準・対象の扱い。利用者の判断
2026-11-07）は、そのテストの docstring に書く。
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.dataset.alignment_files import select_alignment
from metabolomix.dataset.msdial_adapter import load_msdial_folder
from metabolomix.quant import inputs as qi
from metabolomix.quant import merge as qm
from metabolomix.quant import normalize as qn
from metabolomix.quant import plugin as qp
from metabolomix.quant import select as qs
from tests.msdial_folder_fixtures import write_msdial_folder

UNIT = "pmol / 10^6 cells"
REFERENCE = (f"TargetClass,StandardName,Coefficient,Unit\nPC,STD-PC,10,{UNIT}\n"
             f"PE,STD-PE,10,{UNIT}\n")
ADDUCT = "Class,Adduct,Ion mode,IsSelected\nPC,[M+H]+,Positive,TRUE\nPE,[M-H]-,Negative,TRUE\n"
POLS = ("positive", "negative")
#: GUI の特徴 ID → MCP の連番 ID（write_project の行の並び: 10, 42, 90）
ID_TARGET, ID_STD, ID_OTHER, ID_CAR = "0", "1", "2", "3"


@pytest.fixture
def point(request, tmp_path, monkeypatch):
    """`--quant-rules` の本物の規則を指す設定を作る。渡されていなければ skip。"""
    folder = request.config.getoption("--quant-rules")
    if not folder:
        pytest.skip("--quant-rules=<規則パッケージのフォルダ> を渡したときだけ走る")

    def _point(adduct: str | None = ADDUCT):
        lines = [f"rules = '{folder}'"]
        if adduct is not None:
            table = tmp_path / "adduct.csv"
            table.write_text(adduct, encoding="utf-8")
            lines.append(f"adducts = '{table}'")
        cfg = tmp_path / "lipidmix.local.toml"
        cfg.write_text("[quant]\n" + "\n".join(lines) + "\n", encoding="utf-8")
        monkeypatch.setenv("LIPIDMIX_CONFIG", str(cfg))

    _point()
    return _point


def _spot(name, ontology, adduct, heights, mz, rt, gap=None):
    return {"name": name, "ontology": ontology, "adduct": adduct, "heights": heights, "mz": mz,
            "rt": rt, "gap_filled": gap}


PAD = "pad"  # 3 つ目の標本（読み込み側が 3 標本未満のスポットを読まない既知の制限を避ける）


def _load(folder, names, spots, group="A"):
    """標本 2 つ（元のテストの A, B）に 3 つ目 PAD を足して読む（高さは最後の値、gap-fill は無し）。"""
    names = [*names, PAD]
    samples = [{"file_id": i, "name": n, "class_id": group} for i, n in enumerate(names)]
    spots = [{**s, "heights": [*s["heights"], s["heights"][-1]],
              "gap_filled": [*s["gap_filled"], False] if s.get("gap_filled") else None}
             for s in spots]
    write_msdial_folder(folder, spots=spots, samples=samples)
    ds = load_msdial_folder(select_alignment(folder))
    for meta in ds.feature_metadata.values():
        meta["msms_assigned"] = True
        meta["comment"] = ""
    return ds


def _mode(root, mode, names, *, extra=()):
    """`tests/workflow/conftest.py write_inputs` の 1 モード（標準 A,B の高さ 100/50 など）。

    positive: 標本 A,B・PC・対象 [200,300]・標準 [100,50]、negative: 標本 B,A・PE・対象 [100,600]・
    標準 [50,100]。3 行目は未割り当てクラス Other（[3,4]）。
    """
    if mode == "positive":
        cls, heights, std, adduct, std_adduct = "PC", [200, 300], [100, 50], "[M+H]+", "[M+Na]+"
    else:
        cls, heights, std, adduct, std_adduct = "PE", [100, 600], [50, 100], "[M-H]-", "[M+Cl]-"
    spots = [_spot(f"{cls} 34:1", cls, adduct, heights, 760.5, 5.0),
             _spot(f"STD-{cls}", cls, std_adduct, std, 790.5, 5.1),
             _spot("unknown", "Other", adduct, [3, 4], 300.1, 1.0), *extra]
    return _load(root / mode, list(names), [{**s, "heights": [float(h) for h in s["heights"]]}
                                            for s in spots])


def _ref(tmp_path, text=REFERENCE, name="normalization.csv", **kw):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return qi.reference_input(str(path), **kw)


def _prepared(tmp_path, pos_names=("A", "B"), neg_names=("B", "A"), *, normalize=True, ref=None):
    """GUI の prepare_mode 相当（正規化。normalize=False なら飛ばす）。"""
    root = tmp_path / "run"
    out = []
    for mode, names in (("positive", pos_names), ("negative", neg_names)):
        ds = _mode(root, mode, names)
        if normalize:
            ds, _ = qn.apply(ds, ref or _ref(tmp_path), None, mode)
        else:
            ds, _ = qn.skip(ds, mode)
        out.append(ds)
    return out


def _row(ds, fid, names):
    return [float(ds.feature_matrix[ds.feature_ids.index(fid), ds.sample_names.index(n)])
            for n in names]


def _decisions(tmp_path, info):
    name = f"{info['record_id']}.audit.json"
    path = next(tmp_path.rglob(name))  # 統合後は共通の親、統合前は解析フォルダ
    return {d["feature_id"]: d for d in json.loads(path.read_text(encoding="utf-8"))["decisions"]}


# ---- tests/workflow/test_tables.py ----

def test_normalize_then_merge_preserves_values_ids_and_sample_identity(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_normalize_then_merge_preserves_values_ids_and_sample_identity
    （positive:10 は MCP では positive:0）。"""
    pos, neg = _prepared(tmp_path)
    merged, _ = qm.merge(pos, neg, polarities=POLS)
    assert _row(merged, f"positive:{ID_TARGET}", ["A", "B"]) == [20, 60]
    assert _row(merged, f"negative:{ID_TARGET}", ["A", "B"]) == [60, 20]
    assert len(merged.feature_ids) == 6
    new, info = qs.select(merged, adduct=qi.adduct_table())
    assert set(new.feature_ids) == {"positive:0", "negative:0"}
    assert _row(new, "negative:0", ["A"]) == [60]
    assert "reference_not_selected" in info["reasons"]
    assert np.isnan(_row(merged, f"positive:{ID_OTHER}", ["A", "B"])).all()


def test_merge_matches_samples_after_removing_polarity_tokens(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_merge_matches_samples_after_removing_polarity_tokens"""
    pos, neg = _prepared(tmp_path, ("20260930_A_POS_1", "20260930_B_POS_1"),
                         ("20260930_B_neg_1", "20260930_A_neg_1"))
    merged, info = qm.merge(pos, neg, polarities=POLS)
    ids = ["20260930_A_1", "20260930_B_1"]
    assert merged.sample_names == [*ids, PAD]
    assert _row(merged, "positive:0", ids) == [20, 60]
    assert _row(merged, "negative:0", ids) == [60, 20]
    assert info["rule"] == "polarity_token_removed"


@pytest.mark.parametrize(("name", "expected"), [
    ("cell_C_POS_1", "cell_C_1"), ("cell-neg-1", "cell-1"), ("QC_pos1", "QC_1"),
    ("Neg_ctrl", "ctrl"), ("positive_ctrl", "positive_ctrl"),
    ("negative_ctrl_NEG", "negative_ctrl")])
def test_polarity_token_removal_keeps_words_containing_pos_or_neg(point, name, expected):
    """写し元: tests/workflow/test_tables.py::test_polarity_token_removal_keeps_words_containing_pos_or_neg
    （strip_polarity の単体を、統合の照合 `match_samples` の結果で確かめる。相手側は期待する名前）。"""
    plugin = qp.load_rules()
    out = qp.call(plugin, "match_samples",
                  names_by_mode={"positive": [name], "negative": [expected]})
    assert out["ids_by_mode"]["positive"] == [expected]
    assert not any(out["mismatched"].values())


def test_identical_names_are_kept_without_token_removal(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_identical_names_are_kept_without_token_removal"""
    pos, neg = _prepared(tmp_path, ("Pos_ctrl_1", "Neg_ctrl_1"), ("Pos_ctrl_1", "Neg_ctrl_1"))
    merged, info = qm.merge(pos, neg, polarities=POLS)
    assert set(merged.sample_names) == {"Pos_ctrl_1", "Neg_ctrl_1", PAD}
    assert info["rule"] == "exact_name"


def test_token_removal_refuses_names_that_collide(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_token_removal_refuses_names_that_collide
    （規則の sample_names_ambiguous は QUANT_RULES_REJECTED の rule_code になる）。"""
    pos, neg = _prepared(tmp_path, ("pos_ctrl_POS", "neg_ctrl_POS"),
                         ("pos_ctrl_NEG", "neg_ctrl_NEG"))
    with pytest.raises(DomainError) as exc:
        qm.merge(pos, neg, polarities=POLS)
    assert exc.value.code == "QUANT_RULES_REJECTED"
    assert exc.value.details["rule_code"] == "sample_names_ambiguous"


def test_filter_excludes_internal_standard_even_when_adduct_selected(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_filter_excludes_internal_standard_even_when_adduct_selected
    （positive:42 は positive:1。クラス合計の部分は msdial-workflow 側の統計）。"""
    point("Class,Adduct,Ion mode,IsSelected\n"
          "PC,[M+H]+,Positive,TRUE\nPC,[M+Na]+,Positive,TRUE\nPE,[M-H]-,Negative,TRUE\n")
    pos, neg = _prepared(tmp_path)
    merged, _ = qm.merge(pos, neg, polarities=POLS)
    new, info = qs.select(merged, adduct=qi.adduct_table())
    assert set(new.feature_ids) == {"positive:0", "negative:0"}
    assert _decisions(tmp_path, info)[f"positive:{ID_STD}"]["reason"] == "internal_standard"


def test_disabled_adduct_filter_still_excludes_internal_standards(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_disabled_adduct_filter_still_excludes_internal_standards"""
    pos, neg = _prepared(tmp_path)
    merged, _ = qm.merge(pos, neg, polarities=POLS)
    new, info = qs.select(merged, adduct=None)
    assert set(new.feature_ids) == {"positive:0", "positive:2", "negative:0", "negative:2"}
    decisions = _decisions(tmp_path, info)
    excluded = {f: d["reason"] for f, d in decisions.items() if d["action"] == "excluded"}
    assert excluded == {"positive:1": "internal_standard", "negative:1": "internal_standard"}


def test_raw_skip_does_not_require_normalization_reference(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_raw_skip_does_not_require_normalization_reference"""
    ds = _mode(tmp_path / "run", "positive", ("A", "B"))
    new, info = qn.skip(ds, "positive")
    assert _row(new, "0", ["A", "B"]) == [200, 300]
    assert new.quant_state["normalized"] is False and info["skipped"] is True


def test_mismatched_names_require_mapping(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_mismatched_names_require_mapping
    （GUI は ValueError、MCP は QUANT_MERGE_SAMPLE_MISMATCH で止まり対応表を求める）。"""
    pos, neg = _prepared(tmp_path, ("A", "B"), ("X", "Y"))
    with pytest.raises(DomainError) as exc:
        qm.merge(pos, neg, polarities=POLS)
    assert exc.value.code == "QUANT_MERGE_SAMPLE_MISMATCH"


def test_filter_zero_results_keep_audit(tmp_path, point):
    """写し元: tests/workflow/test_tables.py::test_filter_zero_results_and_disabled_filter_keep_audit"""
    point("Class,Adduct,Ion mode,IsSelected\nZZ,[M+H]+,Positive,TRUE\n")
    pos, neg = _prepared(tmp_path)
    merged, _ = qm.merge(pos, neg, polarities=POLS)
    empty, info = qs.select(merged, adduct=qi.adduct_table())
    assert empty.feature_ids == [] and len(_decisions(tmp_path, info)) == 6
    assert len(qs.select(merged, adduct=None)[0].feature_ids) == 4  # 内部標準 2 行だけが除かれる


# ---- tests/workflow/test_any_others.py ----

def _any_others(tmp_path, select_car=True):
    """POS に CAR 行（高さ 400/100）を足し、STD-PC（100/50）を Any others の標準にする。"""
    root = tmp_path / "run"
    car = _spot("CAR 16:0", "CAR", "[M+H]+", [400.0, 100.0], 400.3, 2.0)
    pos = _mode(root, "positive", ("A", "B"), extra=[car])
    ref = _ref(tmp_path, REFERENCE + f"Any others,STD-PC,10,{UNIT}\n", "any_others.csv")
    return pos, ref, ADDUCT + ("CAR,[M+H]+,Positive,TRUE\n" if select_car else "")


def test_unlisted_class_is_normalized_by_any_others_and_kept(tmp_path, point):
    """写し元: tests/workflow/test_any_others.py::test_unlisted_class_is_normalized_by_any_others_and_kept
    （CAR 77 は MCP では 3。Comment 欄の文言は GUI の表の部分なので対象外）。"""
    pos, ref, adduct = _any_others(tmp_path)
    point(adduct)
    new, info = qn.apply(pos, ref, None, "positive")
    assert _row(new, ID_CAR, ["A", "B"]) == [40, 20]
    car = next(c for c in info["classes"] if c["target_class"] == "Any others")
    assert car["via_any_others"] is True and car["feature_id"] == ID_STD
    pc = next(c for c in info["classes"] if c["target_class"] == "PC")
    assert pc["via_any_others"] is False
    kept, _ = qs.select(new, adduct=qi.adduct_table(), polarity="positive")
    assert ID_CAR in kept.feature_ids


def test_any_others_row_still_needs_adduct_selection(tmp_path, point):
    """写し元: tests/workflow/test_any_others.py::test_any_others_row_still_needs_adduct_selection"""
    pos, ref, adduct = _any_others(tmp_path, select_car=False)
    point(adduct)
    new, _ = qn.apply(pos, ref, None, "positive")
    kept, info = qs.select(new, adduct=qi.adduct_table(), polarity="positive")
    assert ID_CAR not in kept.feature_ids
    assert _decisions(tmp_path, info)[ID_CAR]["reason"] == "reference_not_selected"


# ---- tests/normalizer/test_core.py ----

def _core(tmp_path, *, std_heights=(100.0, 50.0), std_gap=None, target_gap=None,
          target=(200.0, 300.0), extra=()):
    """tests/normalizer/conftest.py の dataset: 標本 A,B・PC 34:1 [200,300]・標準 [100,50]・Unmapped [5,6]。"""
    spots = [
        _spot("PC 34:1", "PC", "[M+H]+", list(target), 760.5, 5.0, gap=target_gap),
        _spot("PC 33:1(d7)|PC 15:0_18:1(d7)", "PC", "[M+H]+", list(std_heights), 790.5, 5.1,
              gap=std_gap),
        _spot("unassigned", "Unmapped", "[M+H]+", [5.0, 6.0], 300.1, 1.0), *extra]
    return _load(tmp_path / "core", ["A", "B"], spots)


def _core_ref(tmp_path, standard="PC 15:0_18:1(d7)", unit=UNIT, coefficient="10"):
    return _ref(tmp_path, f"TargetClass,StandardName,Coefficient,Unit\nPC,{standard},{coefficient},"
                f"{unit}\n", "reference.csv")


def test_hand_calculated_same_sample_ratio_no_extra_cell_scaling(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_hand_calculated_same_sample_ratio_no_extra_cell_scaling"""
    ds = _core(tmp_path)
    new, info = qn.apply(ds, _core_ref(tmp_path), None, "positive")
    assert _row(new, "0", ["A", "B"]) == [20.0, 60.0]
    # 標準自身の行は値を持たない（GUI の表では標準の自己比 10 が出るが、選別で内部標準として除かれる）
    assert np.isnan(_row(new, "1", ["A", "B"])).all()
    assert np.isnan(_row(new, "2", ["A", "B"])).all()  # class_unmapped
    assert info["status_counts"]["class_unmapped"] == 3  # 2 標本 + PAD
    assert info["classes"][0]["feature_id"] == "1" and info["internal_standard_ids"] == ["1"]


@pytest.mark.parametrize("denominator", [0.0, -2.0, float("nan"), float("inf")])
def test_invalid_standard_remains_missing(tmp_path, point, denominator):
    """写し元: tests/normalizer/test_core.py::test_invalid_standard_remains_missing"""
    ds = _core(tmp_path, std_heights=(denominator, 50.0))
    new, info = qn.apply(ds, _core_ref(tmp_path), None, "positive")
    a, b = _row(new, "0", ["A", "B"])
    assert np.isnan(a) and b == 60.0
    assert info["status_counts"]["standard_invalid_height"] == 1


def test_gap_filled_standard_is_computed_and_marked(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_gap_filled_standard_is_not_used。
    意図して違う（利用者の判断 2026-11-07）: GUI（tables.py）と同じく gap-fill の内部標準でも値を計算し、
    detected_mask で gap-fill 由来の印を引き継ぐ（正規化の単体の「内部標準が gap-fill なら欠測」は使わない）。"""
    ds = _core(tmp_path, std_gap=[True, False])
    new, _ = qn.apply(ds, _core_ref(tmp_path), None, "positive")
    assert _row(new, "0", ["A", "B"]) == [20.0, 60.0]
    row = new.feature_ids.index("0")
    assert [bool(x) for x in new.detected_mask[row]][:2] == [False, True]


def test_ambiguous_standard_has_no_automatic_first_candidate(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_ambiguous_standard_has_no_automatic_first_candidate
    （複製の標準は MCP では 3。対応指定の鍵は内部標準の名前）。"""
    twin = _spot("PC 33:1(d7)|PC 15:0_18:1(d7)", "PC", "[M+H]+", [100.0, 50.0], 790.6, 5.2)
    ds = _core(tmp_path, extra=[twin])
    ref = _core_ref(tmp_path)
    pv = qn.preview(ds, ref, None, "positive")
    assert pv["status"] == "needs_input"
    assert pv["unresolved"][0]["status"] == "standard_ambiguous"
    assert {c["feature_id"] for c in pv["unresolved"][0]["candidates"]} == {"1", "3"}
    assert qn.apply(ds, ref, None, "positive")[0] is None  # 先頭の候補を自動で選ばない
    new, _ = qn.apply(ds, ref, {"PC 15:0_18:1(d7)": 3}, "positive")
    assert _row(new, "0", ["A", "B"]) == [20.0, 60.0]
    with pytest.raises(DomainError) as exc:  # 名前が一致する候補ではない ID
        qn.apply(ds, ref, {"PC 15:0_18:1(d7)": 0}, "positive")
    assert exc.value.code == "QUANT_RULES_REJECTED"


def test_missing_standard_does_not_fall_back_to_another_class(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_missing_standard_does_not_fall_back_to_another_class"""
    ds = _core(tmp_path)
    ref = _core_ref(tmp_path, standard="PC 99:9(d7)")
    pv = qn.preview(ds, ref, None, "positive")
    assert pv["status"] == "needs_input"
    assert pv["unresolved"][0]["status"] == "standard_missing"
    assert qn.apply(ds, ref, None, "positive")[0] is None


def test_target_zero_is_a_real_zero_and_gap_fill_is_flagged(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_target_zero_is_a_real_zero_and_gap_fill_is_flagged
    （値 (0, 60)。gap-fill の印は detected_mask で引き継ぐ）。"""
    ds = _core(tmp_path, target=(0.0, 300.0), target_gap=[False, True])
    new, _ = qn.apply(ds, _core_ref(tmp_path), None, "positive")
    assert _row(new, "0", ["A", "B"]) == [0.0, 60.0]
    row = new.feature_ids.index("0")
    assert [bool(x) for x in new.detected_mask[row]][:2] == [True, False]


@pytest.mark.parametrize("bad", ["-1", "nan", "0"])
def test_reference_rejects_invalid_coefficients(tmp_path, point, bad):
    """写し元: tests/normalizer/test_core.py::test_reference_rejects_invalid_coefficients"""
    with pytest.raises(DomainError) as exc:
        qn.apply(_core(tmp_path), _core_ref(tmp_path, coefficient=bad), None, "positive")
    assert exc.value.code == "QUANT_RULES_REJECTED"


def test_reference_rejects_duplicate_classes_and_wrong_unit(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_reference_rejects_duplicate_classes_and_wrong_unit"""
    dup = REFERENCE + f"PC,STD-PC,10,{UNIT}\n"
    with pytest.raises(DomainError) as exc:
        qn.apply(_core(tmp_path), _ref(tmp_path, dup), None, "positive")
    assert exc.value.code == "QUANT_RULES_REJECTED"
    with pytest.raises(DomainError):
        qn.apply(_core(tmp_path), _core_ref(tmp_path, unit="pmol/uL plasma"), None, "positive")


def test_ambiguous_10e6_unit_spelling_is_not_accepted_as_one_million(tmp_path, point):
    """写し元: tests/normalizer/test_core.py::test_ambiguous_10e6_unit_spelling_is_not_accepted_as_one_million"""
    with pytest.raises(DomainError) as exc:
        qn.apply(_core(tmp_path), _core_ref(tmp_path, unit="pmol_per_10e6_cells"), None, "positive")
    assert exc.value.code == "QUANT_RULES_REJECTED" and "単位" in exc.value.message


@pytest.mark.parametrize(("spelling", "canonical"), [
    ("nmol/10^6 cells", "nmol / 10^6 cells"), ("FMOL / CELL", "fmol / cell"),
    ("pmol / 10^3 cells", "pmol / 10^3 cells"), ("pmol_per_1e6_cells", "pmol / 10^6 cells"),
    ("nmol/10^4cells", "nmol / 10^4 cells"), ("fmol / 10^0 cells", "fmol / cell")])
def test_supported_coefficient_units_are_canonicalised(tmp_path, point, spelling, canonical):
    """写し元: tests/normalizer/test_core.py::test_supported_coefficient_units_are_canonicalised"""
    _new, info = qn.apply(_core(tmp_path), _core_ref(tmp_path, unit=spelling), None, "positive")
    assert info["unit"] == canonical


@pytest.mark.parametrize("unit", ["pmol / 10^13 cells", "pmol / 10^-1 cells", "pmol / 10^1.5 cells"])
def test_out_of_range_cell_basis_is_rejected(tmp_path, point, unit):
    """写し元: tests/normalizer/test_core.py::test_out_of_range_cell_basis_is_rejected"""
    with pytest.raises(DomainError) as exc:
        qn.apply(_core(tmp_path), _core_ref(tmp_path, unit=unit), None, "positive")
    assert exc.value.code == "QUANT_RULES_REJECTED" and "単位" in exc.value.message
