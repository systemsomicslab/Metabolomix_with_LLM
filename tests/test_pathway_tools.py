# tests/test_pathway_tools.py
import json
import sys

import fitz
import numpy as np
import pytest
from mcp.server.fastmcp import Image

from metabolomix.core import mcp_core, session_state
from metabolomix.pathway import bridge, tools

from tests.test_pathway_jsonld import B, _graph  # 小さな JSON-LD 組み立て器を再利用


@pytest.fixture
def lipid_arf(monkeypatch, tmp_path):
    """前処理済み ARF ＋兄弟 .arf2 カタログを持つセッション（ファイル I/O は差し替え）。"""
    s = session_state.AnalysisSession()
    monkeypatch.setattr(session_state, "session", s)
    monkeypatch.setattr(mcp_core, "DATA_DIR", tmp_path)
    monkeypatch.delenv("LIPIDMIX_REPORTS_DIR", raising=False)
    s.assay_kind = "lipid"
    a = s.arf
    a.features = [{"MasterAlignmentID": i} for i in (1, 2)]
    a.current_file_path = str(tmp_path / "AlignResult-1.arf")
    a.feature_matrix = np.array([[10.0, 1.0], [11.0, 1.1], [30.0, 0.3], [31.0, 0.31]])
    a.pp_sample_names = ["C_1", "C_2", "P_1", "P_2"]
    a.pp_feature_names = ["Spot_1_height", "Spot_2_height"]
    a.sample_meta = {n: {"role": "sample", "group": n.split("_")[0]} for n in a.pp_sample_names}
    a.preprocessing_recipe = {"normalize": "none", "props": ["height"]}
    catalog = [{"MasterAlignmentID": 1, "Name": "PC 34:1", "Ontology": "PC"},
               {"MasterAlignmentID": 2, "Name": "LPC 18:1", "Ontology": "LPC"}]
    monkeypatch.setattr(tools, "_arf2_catalog", lambda: (catalog, None))
    return s


#: fake_measure が JSON-LD に足す要素（テストごとに差し替える）
_EXTRA: list = []


@pytest.fixture
def fake_measure(monkeypatch):
    """run_bridge を差し替え、measure で JSON-LD を書く。"""
    calls = []
    monkeypatch.setattr(sys.modules[__name__], "_EXTRA", [])

    def run(request, run_dir, *, timeout_s=300.0):
        calls.append(request)
        out = {}
        for p in request["pathways"]:
            path = run_dir / f"{p}-measurements.json"
            doc = _graph(extra=_EXTRA)
            # 群名を実際の比較に合わせる
            for x in doc["@graph"]:
                if x["@type"] == "mslm:Sample":
                    x["sampleGroup"] = "C" if x["identifier"] == "S1" else "P"
            path.write_text(json.dumps(doc), encoding="utf-8")
            out[p] = str(path)
        return {"status": "ok", "outputs": out, "mapped_nodes": {},
                "package_version": "0.1", "package_sha": "abc123"}

    monkeypatch.setattr(bridge, "run_bridge", run)
    return calls


def test_requires_lipid_assay(lipid_arf, fake_measure):
    lipid_arf.assay_kind = "metabolite"
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["code"] == "missing_state"
    assert body["error"]["required_tools"] == ["record_objective"]


def test_requires_preprocessed_matrix(lipid_arf, fake_measure):
    lipid_arf.arf.feature_matrix = None
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["required_tools"] == ["arf_preprocess"]


def test_unknown_pathway_is_bad_request(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["nope"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_project_writes_run_and_returns_tables(lipid_arf, fake_measure, tmp_path):
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))
    run_id = body["run_id"]
    run_dir = tmp_path / "reports" / "pathway" / run_id
    assert (run_dir / "input.mztab").is_file()
    assert (run_dir / "nodes.tsv").is_file()
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["upstream_commit"] == "abc123"
    assert manifest["schema_version"] == "2.0.0"
    assert manifest["mapped_nodes"] == {"remodeling": body["pathways"]["remodeling"]["mapped_nodes"]}
    assert manifest["groups"] == {"C": 2, "P": 2}
    assert body["pathways"]["remodeling"]["nodes_tsv"].startswith(
        "node\tscope\tn_features\tmean_a\tmean_b\tlog2fc\tq\tcall")
    assert body["pathways"]["remodeling"]["unmatched_subclasses"] == ["TG"]
    assert "not curated reactions" in body["notes"]["adjacency"] or "反応" in body["notes"]["adjacency"]
    assert lipid_arf.pathway.last_run_id == run_id
    req = fake_measure[0]
    assert req["mode"] == "measure" and req["networks"]["remodeling"] == "remodeling_lipidonly.cyjs"


def _empty_scope_items(node, scope):
    """上流の空スコープ（memberFeature: []、値は全部 null）。"""
    sid = f"{B}n/{node}/scope/{scope}"
    chains = [] if scope == "whole" else [scope.removeprefix("acyl-")]
    return [{"@id": sid, "@type": "mslm:LipidMeasurementScope", "onPathwayNode": {"@id": f"{B}n/{node}"},
             "scopeType": "subclass" if scope == "whole" else "acyl_subset",
             "acylChainSignature": chains, "memberFeature": []},
            {"@id": sid + "/m1", "@type": "mslm:SampleMeasurement", "measuredEntity": {"@id": sid},
             "inSample": {"@id": f"{B}s/1"}, "value": None}]


def test_mapped_nodes_counts_nonempty_whole_over_lipid_nodes(lipid_arf, fake_measure, monkeypatch):
    # PE はノードがあるがデータに無い（空の whole スコープ）。分母は脂質ノード
    # （PC/PE/LPC/DG）だけで、PathwayEntity（PA）と接合ノードは数えない
    monkeypatch.setattr(sys.modules[__name__], "_EXTRA", _empty_scope_items("PE", "whole"))
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))
    res = body["pathways"]["remodeling"]
    assert res["mapped_nodes"] == "1/4"
    assert {l.split("\t")[0] for l in res["nodes_tsv"].splitlines()[1:]} == {"PC"}   # PE は行にならない


def test_requested_chain_absent_everywhere_is_bad_request(lipid_arf, fake_measure, monkeypatch):
    # 22:6 は上流が空スコープとして出すだけ（どのパスウェイにも分子種が無い）
    monkeypatch.setattr(sys.modules[__name__], "_EXTRA", _empty_scope_items("PC", "acyl-22:6"))
    body = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"],
                                            acyl_chains=["18:0", "22:6"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert body["error"]["details"]["missing_chains"] == ["22:6"]
    ok = json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"],
                                          acyl_chains=["18:0"]))
    assert "run_id" in ok


def test_same_groups_are_bad_request(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(group_a="C", group_b="C"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    body = json.loads(tools.pathway_project(source="dataset", group_a=["a"], group_b=["b"],
                                            group_a_label="X", group_b_label="X"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert fake_measure == []


def test_bridge_error_becomes_envelope(lipid_arf, monkeypatch):
    def boom(*a, **k):
        raise bridge.BridgeError("MSLM_NOT_CONFIGURED", "設定して", {"setting": "x"})
    monkeypatch.setattr(bridge, "run_bridge", boom)
    body = json.loads(tools.pathway_project(group_a="C", group_b="P"))
    assert body["error"]["code"] == "MSLM_NOT_CONFIGURED"


# --- pathway_plot ---

def _write_pdf(path):
    doc = fitz.open()
    doc.new_page(width=842, height=595).insert_text((72, 72), "PC")
    doc.save(path)


def test_plot_without_run_is_missing_state(lipid_arf):
    body = json.loads(tools.pathway_plot())
    assert body["error"]["required_tools"] == ["pathway_project"]


def test_plot_renders_pngs(lipid_arf, fake_measure, monkeypatch):
    run_id = json.loads(tools.pathway_project(group_a="C", group_b="P",
                                              pathways=["remodeling", "ceramide"]))["run_id"]
    seen = []

    def render(request, run_dir, *, timeout_s=300.0):
        seen.append(request)
        out = {}
        for p in request["pathways"]:
            pdf = run_dir / f"{p}.pdf"
            _write_pdf(pdf)
            out[p] = str(pdf)
        return {"status": "ok", "outputs": out, "mapped_nodes": {}}

    monkeypatch.setattr(bridge, "run_bridge", render)
    out = tools.pathway_plot(pathways=["ceramide"])
    assert isinstance(out, list) and isinstance(out[0], str) and isinstance(out[1], Image)
    assert len(out) == 2
    assert "ceramide.pdf" in out[0] and tools._FIGURE_NOTE in out[0]
    ceramide_line = next(l for l in out[0].splitlines() if "ceramide.pdf" in l)
    manifest = json.loads((tools._run_dir(run_id) / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["mapped_nodes"]["ceramide"] in ceramide_line
    assert seen[0]["mode"] == "render" and seen[0]["x_order"] == ["C", "P"]
    assert seen[0]["mztab"].endswith("input.mztab")
    assert run_id in seen[0]["out_dir"]


def test_plot_pathway_not_in_run_is_bad_request(lipid_arf, fake_measure):
    json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))
    body = json.loads(tools.pathway_plot(pathways=["global"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_render_failure_envelope(lipid_arf, fake_measure, monkeypatch):
    json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))

    def fail(*a, **k):
        raise bridge.BridgeError("PATHWAY_RENDER_FAILED", "Chromium がありません")
    monkeypatch.setattr(bridge, "run_bridge", fail)
    body = json.loads(tools.pathway_plot())
    assert body["error"]["code"] == "PATHWAY_RENDER_FAILED"


def test_unreadable_pdf_becomes_render_failed(lipid_arf, fake_measure, monkeypatch):
    """PDF が壊れている（0 ページ・読めない）なら例外を漏らさず PATHWAY_RENDER_FAILED。"""
    json.loads(tools.pathway_project(group_a="C", group_b="P", pathways=["remodeling"]))

    def render(request, run_dir, *, timeout_s=300.0):
        pdf = run_dir / "remodeling.pdf"
        pdf.write_bytes(b"not a pdf")
        return {"status": "ok", "outputs": {"remodeling": str(pdf)}}

    monkeypatch.setattr(bridge, "run_bridge", render)
    body = json.loads(tools.pathway_plot())
    assert body["error"]["code"] == "PATHWAY_RENDER_FAILED"
    assert "remodeling" in body["error"]["message"]


@pytest.mark.parametrize("bad", ["../..", r"..\..", "pw_x", "pw_20261014_120000_abcdef/../x"])
def test_plot_rejects_malformed_run_id(lipid_arf, bad):
    body = json.loads(tools.pathway_plot(run_id=bad))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_corrupt_manifest_is_output_invalid(lipid_arf, fake_measure):
    run_id = json.loads(tools.pathway_project(group_a="C", group_b="P",
                                              pathways=["remodeling"]))["run_id"]
    (tools._run_dir(run_id) / "manifest.json").write_text("{not json", encoding="utf-8")
    body = json.loads(tools.pathway_plot(run_id=run_id))
    assert body["error"]["code"] == "MSLM_OUTPUT_INVALID"
