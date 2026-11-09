# tests/test_pathway_compare.py
import pytest
from metabolomix.pathway import compare
from metabolomix.pathway.jsonld import PathwayGraph, Scope


def _graph(values_by_node: dict[str, list], edges, kinds=None):
    """values_by_node: label -> [S1..S6]（S1-3=C, S4-6=P）。None は欠損。"""
    g = PathwayGraph()
    for label in values_by_node:
        nid = f"n/{label}"
        g.node_kinds[nid] = "lipid"
        g.node_labels[nid] = label
        sid = f"{nid}/scope/whole"
        g.scopes[sid] = Scope(sid, nid, "subclass", (), 3)
        g.values[sid] = {f"S{i}": v for i, v in enumerate(values_by_node[label], start=1)}
    for nid, kind in (kinds or {}).items():
        g.node_kinds[nid] = kind
        g.node_labels[nid] = ""
    g.edges = [(f"n/{a}" if not a.startswith("j") else a,
                f"n/{b}" if not b.startswith("j") else b) for a, b in edges]
    g.sample_groups = {f"S{i}": ("C" if i <= 3 else "P") for i in range(1, 7)}
    return g


UP = [10, 11, 10.5, 100, 110, 105]
DOWN = [100, 110, 105, 10, 11, 10.5]
FLAT = [10, 50, 30, 12, 48, 31]


def test_up_down_ns_calls_and_direction():
    g = _graph({"PC": UP, "LPC": DOWN, "PE": FLAT}, [("PC", "LPC")])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert rows["PC"]["call"] == "up" and rows["PC"]["log2fc"] > 0
    assert rows["LPC"]["call"] == "down"
    assert rows["PE"]["call"] == "ns"
    assert rows["PC"]["scope"] == "whole" and rows["PC"]["n_features"] == 3


def test_all_missing_in_one_group_is_na_and_not_in_bh():
    g = _graph({"PC": UP, "TG": [1, 2, 3, None, None, None]}, [])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert rows["TG"]["call"] == "NA"
    assert rows["TG"]["q"] is None
    # TG が BH の分母に入らない＝ PC の q は単独検定の p と一致
    single = compare.compare_scopes(_graph({"PC": UP}, []), "C", "P", 0.05)[0]
    assert rows["PC"]["q"] == single["q"]


def test_bh_is_per_scope_kind():
    g = _graph({"PC": UP}, [])
    g.scopes["n/PC/scope/acyl-18:0"] = Scope("n/PC/scope/acyl-18:0", "n/PC", "acyl_subset", ("18:0",), 1)
    g.values["n/PC/scope/acyl-18:0"] = {f"S{i}": v for i, v in enumerate(FLAT, start=1)}
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    assert [r["scope"] for r in rows] == ["whole", "18:0"]
    whole_alone = compare.compare_scopes(_graph({"PC": UP}, []), "C", "P", 0.05)[0]
    assert rows[0]["q"] == whole_alone["q"]   # 鎖スコープが whole の BH に混ざらない


def test_opposite_pairs_needs_both_significant():
    g = _graph({"PC": UP, "LPC": DOWN, "PE": FLAT, "PS": UP},
               [("PC", "j1"), ("j1", "LPC"), ("PC", "PE"), ("PS", "PC")],
               kinds={"j1": "junction"})
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    pairs, undeterminable = compare.opposite_pairs(g, rows)
    assert [(p["node_1"], p["node_2"], p["via"]) for p in pairs] == [("LPC", "PC", "junction")]
    assert undeterminable == 1          # PC–PE（PE が ns）。PS–PC は同じ向きで不掲載・判定不能でもない


def _empty_scope(g, label, scope_type="subclass", chains=()):
    """上流が出す空スコープ（memberFeature: []、値は全部 null。NotMeasured）を足す。"""
    nid = f"n/{label}"
    g.node_kinds[nid] = "lipid"
    g.node_labels[nid] = label
    sid = f"{nid}/scope/" + ("whole" if scope_type == "subclass" else "acyl-" + "+".join(chains))
    g.scopes[sid] = Scope(sid, nid, scope_type, tuple(chains), 0)
    g.values[sid] = {f"S{i}": None for i in range(1, 7)}
    return nid


def test_empty_scopes_are_not_rows():
    g = _graph({"PC": UP, "LPC": DOWN}, [("PC", "LPC")])
    _empty_scope(g, "PS")                                   # データに無いサブクラスのノード
    _empty_scope(g, "PC", "acyl_subset", ("22:6",))          # データに無い鎖
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    assert [(r["node"], r["scope"]) for r in rows] == [("LPC", "whole"), ("PC", "whole")]
    # 空スコープは BH の分母にも入らない
    alone = compare.compare_scopes(_graph({"PC": UP, "LPC": DOWN}, [("PC", "LPC")]), "C", "P", 0.05)
    assert [r["q"] for r in rows] == [r["q"] for r in alone]


def test_pairs_touching_empty_node_are_excluded_not_undeterminable():
    g = _graph({"PC": UP, "LPC": DOWN}, [("PC", "LPC"), ("PC", "PS"), ("PS", "LPC")])
    _empty_scope(g, "PS")
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    pairs, undeterminable = compare.opposite_pairs(g, rows)
    assert [(p["node_1"], p["node_2"]) for p in pairs] == [("LPC", "PC")]
    assert undeterminable == 0


def test_tsv_rounds_floats():
    tsv = compare.rows_to_tsv([{"a": 1.23456789, "b": None, "c": "x"}], ["a", "b", "c"])
    assert tsv == "a\tb\tc\n1.235\t\tx"


# --- 統計エンジン（moderated t）への移行 ----------------------------------------------

def _raw_moderated(values_by_node, group_a="C", group_b="P"):
    """スコープの値から同じ行列を作って、エンジンの compare を直に呼ぶ（期待値）。"""
    import numpy as np
    from metabolomix.analysis.stats.compare import compare
    from metabolomix.analysis.stats.external import prepared_from_injections

    labels = list(values_by_node)
    log2 = np.array([[np.log2(values_by_node[n][i]) for n in labels] for i in range(6)])
    prepared = prepared_from_injections(log2, sample_ids=[f"S{i}" for i in range(1, 7)],
                                        groups=["C"] * 3 + ["P"] * 3, subject_ids=[None] * 6,
                                        feature_ids=labels)
    return compare(prepared, factors=["group"], comparisons=[
        {"id": "c", "factor": "group", "reference": group_a, "test": group_b}],
        test="moderated_t")["rows"]


def test_compare_scopes_uses_engine_and_bh_per_scope_kind():
    vals = {"PC": UP, "LPC": DOWN, "PE": FLAT}
    g = _graph(vals, [])
    g.scopes["n/PC/scope/acyl-18:0"] = Scope("n/PC/scope/acyl-18:0", "n/PC", "acyl_subset", ("18:0",), 1)
    g.values["n/PC/scope/acyl-18:0"] = {f"S{i}": v for i, v in enumerate(FLAT, start=1)}
    rows = compare.compare_scopes(g, "C", "P", 0.05)
    whole = {r["node"]: r for r in rows if r["scope"] == "whole"}
    expected = {r["feature_id"]: r for r in _raw_moderated(vals)}
    for node, row in whole.items():
        e = expected[node]
        assert row["p"] == pytest.approx(e["p"]) and row["q"] == pytest.approx(e["q"])
        assert row["log2fc"] == pytest.approx(e["log2fc"])
        assert row["ci_low"] == pytest.approx(e["ci_low"]) and row["ci_high"] == pytest.approx(e["ci_high"])
        assert row["verdict"] == e["verdict"]
    # 鎖スコープは別の BH の単位（whole の 3 件に混ざらない）
    chain = [r for r in rows if r["scope"] == "18:0"]
    assert len(chain) == 1 and chain[0]["q"] == pytest.approx(chain[0]["p"])


def test_call_up_down_ns_na_from_verdict():
    g = _graph({"PC": UP, "LPC": DOWN, "PE": FLAT, "TG": [1, 2, 3, None, None, None]}, [])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert (rows["PC"]["verdict"], rows["PC"]["call"]) == ("changed", "up")
    assert (rows["LPC"]["verdict"], rows["LPC"]["call"]) == ("changed", "down")
    assert rows["PE"]["call"] == "ns" and rows["PE"]["verdict"] != "changed"
    assert rows["TG"]["call"] == "NA" and rows["TG"]["p"] is None
    # q 閾値は判定の q ≤ 閾値。厳しくすれば changed が ns に落ちる
    strict = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 1e-12)}
    assert strict["PC"]["call"] == "ns"


def test_nonpositive_scope_values_are_missing():
    g = _graph({"PC": UP, "ZERO": [0, 0, 0, 5, 6, 7]}, [])
    rows = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert rows["ZERO"]["call"] == "NA"      # 0 は欠測（log2 を取れない）。C 群に値が無い


def test_replicates_collapsed_by_subject_of():
    # 6 注入 = 3 個体 × 2 反復にまとめると、各群 1 個体になり比較できない
    g = _graph({"PC": UP, "LPC": DOWN}, [])
    subject_of = {"S1": "x", "S2": "x", "S3": "y", "S4": "z", "S5": "z", "S6": "w"}
    rows = compare.compare_scopes(g, "C", "P", 0.05, subject_of=subject_of)    # C: x, y / P: z, w
    assert len(rows) == 2
    collapsed = {r["node"]: r for r in rows}
    # 手計算: log2 の平均（個体 x = S1,S2、z = S4,S5）で比べた log2fc
    import numpy as np
    pc = np.log2(UP)
    mean_c, mean_p = np.mean([np.mean(pc[0:2]), pc[2]]), np.mean([np.mean(pc[3:5]), pc[5]])
    assert collapsed["PC"]["log2fc"] == pytest.approx(mean_p - mean_c)
    alone = {r["node"]: r for r in compare.compare_scopes(g, "C", "P", 0.05)}
    assert collapsed["PC"]["p"] != pytest.approx(alone["PC"]["p"])    # 個体数が違えば検定も違う
    # 各群 1 個体なら比較できない
    from metabolomix.core.atomic_io import DomainError
    one = {"S1": "x", "S2": "x", "S3": "x", "S4": "z", "S5": "z", "S6": "z"}
    with pytest.raises(DomainError) as exc:
        compare.compare_scopes(g, "C", "P", 0.05, subject_of=one)
    assert exc.value.code == "STATISTIC_GROUP_TOO_SMALL"


def test_pathway_module_does_not_import_differential():
    """`metabolomix/pathway/` のどのモジュールも differential を import せず、属性でも名前でも触らない。"""
    import ast
    from pathlib import Path

    root = Path(compare.__file__).parent
    offenders = []
    for path in sorted(root.glob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                hit = any(a.name.split(".")[-1] == "differential" for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                hit = ((node.module or "").split(".")[-1] == "differential"
                       or any(a.name == "differential" for a in node.names))
            elif isinstance(node, ast.Attribute):
                hit = node.attr == "differential"
            elif isinstance(node, ast.Name):
                hit = node.id == "differential"
            else:
                hit = False
            if hit:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders
