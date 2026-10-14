# tests/test_pathway_compare.py
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


def test_tsv_rounds_floats():
    tsv = compare.rows_to_tsv([{"a": 1.23456789, "b": None, "c": "x"}], ["a", "b", "c"])
    assert tsv == "a\tb\tc\n1.235\t\tx"
