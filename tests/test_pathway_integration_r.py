# tests/test_pathway_integration_r.py
"""実物の R（renv の仮想環境 r/）で MSLipidMapper の measure を 1 回通す（spec §6）。

既定ではスキップする。走らせるには:
    LIPIDMIX_TEST_MSLM=1 LIPIDMIX_TEST_RSCRIPT=<Rscript のパス> pytest tests/test_pathway_integration_r.py
前提: `Rscript r/setup.R` で仮想環境を復元済みであること。ネットワークには出ない。
"""
import os
import sys

import numpy as np
import pytest

from metabolomix.pathway import bridge, compare, handoff, jsonld

pytestmark = [
    pytest.mark.skipif(os.environ.get("LIPIDMIX_TEST_MSLM") != "1",
                       reason="実物の R を使う結合テスト（LIPIDMIX_TEST_MSLM=1 で有効）"),
    pytest.mark.skipif(sys.platform != "win32", reason="start_owned_process は Windows 専用"),
]


@pytest.fixture
def rscript(monkeypatch):
    path = os.environ.get("LIPIDMIX_TEST_RSCRIPT")
    if not path:
        pytest.skip("LIPIDMIX_TEST_RSCRIPT（Rscript のパス）が未設定")
    monkeypatch.setenv("MSLIPIDMAPPER_RSCRIPT", path)
    return path


def _table() -> handoff.HandoffTable:
    names = ["PC 16:0_18:1", "PC 18:0_18:1", "PC 16:0_20:4",
             "PE 16:0_18:1", "PE 18:0_20:4",
             "LPC 18:1", "LPC 16:0"]
    ontologies = ["PC", "PC", "PC", "PE", "PE", "LPC", "LPC"]
    rng = np.random.default_rng(0)
    matrix = rng.uniform(1e4, 5e4, size=(len(names), 6))
    matrix[:, 3:] *= np.array([[2.0], [2.0], [2.0], [0.5], [0.5], [1.0], [1.0]])
    matrix[1, 4] = np.nan                     # 欠損 1 か所
    return handoff.HandoffTable(
        sample_names=["C_1", "C_2", "C_3", "P_1", "P_2", "P_3"],
        groups=["C", "C", "C", "P", "P", "P"],
        names=names, ontologies=ontologies, matrix=matrix)


def test_measure_remodeling_with_real_mslipidmapper(rscript, tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    mztab = run_dir / "input.mztab"
    handoff.write_mztab(_table(), mztab, normalization="none")
    request = {"mode": "measure", "mztab": str(mztab), "out_dir": str(run_dir),
               "pathways": ["remodeling"],
               "networks": {"remodeling": bridge.PATHWAY_NETWORKS["remodeling"]},
               "acyl_chains": ["18:1"], "aggregation": "sum",
               "normalization": "none", "unit": "arbitrary", "dataset_id": "integration"}
    result = bridge.run_bridge(request, run_dir, timeout_s=600)

    assert result["status"] == "ok"
    assert isinstance(result["mapped_nodes"]["remodeling"], int)
    graph = jsonld.load_graph(result["outputs"]["remodeling"])   # schemaVersion を検証する
    rows = compare.compare_scopes(graph, "C", "P", 0.05)

    whole = {r["node"]: r for r in rows if r["scope"] == "whole"}
    for node in ("PC", "PE", "LPC"):
        assert node in whole, node
        assert whole[node]["n_features"] >= 2
        scope_id = next(sid for sid, sc in graph.scopes.items()
                        if sc.node_id == whole[node]["node_id"] and sc.scope_type == "subclass")
        vals = graph.values[scope_id]
        by_group = {g: [vals.get(s) for s, gg in graph.sample_groups.items() if gg == g]
                    for g in ("C", "P")}
        assert all(len(v) == 3 and all(x is not None for x in v) for v in by_group.values()), node
    assert any(r["scope"] == "18:1" and r["n_features"] > 0 for r in rows)
