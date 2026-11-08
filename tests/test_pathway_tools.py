# tests/test_pathway_tools.py
import json
import sys
from pathlib import Path

import pytest
from mcp.server.fastmcp import Image

from metabolomix.core import session_state
from metabolomix.pathway import bridge, tools

from tests.pymupdf_guard import require_pymupdf
from tests.test_pathway_jsonld import B, _graph  # 小さな JSON-LD 組み立て器を再利用


@pytest.fixture
def lipid_arf(monkeypatch, tmp_path):
    """統計エンジンの行列を作れるデータセット（脂質アッセイ）を持つセッション。

    名前は旧 ARF 経路の名残（テストの本数が多いので変えない）。中身は `load_dataset` で
    開いたデータセットと同じ `DatasetState`（C0_r0 C1_r0 / P0_r0 P1_r0 の 2 群 × 2 個体。f0・f1 だけが脂質注釈）。"""
    from tests.test_stats_results import _ds

    s = session_state.AnalysisSession()
    monkeypatch.setattr(session_state, "session", s)
    s.assay_kind = "lipid"
    ds = _ds(groups=("C", "P"), per_group=2, reps=1, n_qc=0)
    ds.analysis_dir = str(tmp_path)   # 記録（run）の置き場所の基準: <tmp>/reports/pathway
    ds.feature_annotations = {"f0": {"name": "PC 34:1", "ontology": "PC"},
                              "f1": {"name": "LPC 18:1", "ontology": "LPC"}}
    s.dataset = ds
    return s


#: 2 群の指定（群はサンプル名のリスト。群名はラベルが決める）
_GROUPS = {"group_a": ["C0_r0", "C1_r0"], "group_b": ["P0_r0", "P1_r0"],
           "group_a_label": "C", "group_b_label": "P"}
_RUN_SAMPLES = {"C0_r0": "C", "C1_r0": "C", "P0_r0": "P", "P1_r0": "P"}


def with_run_samples(doc, values=(10.0, 12.0, 30.0, 33.0)):
    """`_graph` の JSON-LD の試料を、pathway_project が書く mzTab-M の試料名（4 注入）に差し替える。
    PC の whole スコープと 18:0 の鎖スコープに値を付ける（鎖は群に差が無い値）。"""
    items = [x for x in doc["@graph"]
             if not (x["@type"] == "mslm:Sample" or x["@type"] == "mslm:SampleMeasurement")]
    for i, (name, group) in enumerate(_RUN_SAMPLES.items(), start=1):
        items.append({"@id": f"{B}rs/{i}", "@type": "mslm:Sample", "identifier": name, "sampleGroup": group})
        for scope, v in (("n/PC/scope/whole", values[i - 1]), ("n/PC/scope/acyl-18%3A0", 5.0 + i % 2)):
            items.append({"@id": f"{B}rm/{scope}/{i}", "@type": "mslm:SampleMeasurement",
                          "measuredEntity": {"@id": B + scope}, "inSample": {"@id": f"{B}rs/{i}"},
                          "value": v})
    return {**doc, "@graph": items}


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
            doc = with_run_samples(_graph(extra=_EXTRA))
            path.write_text(json.dumps(doc), encoding="utf-8")
            out[p] = str(path)
        return {"status": "ok", "outputs": out, "mapped_nodes": {},
                "package_version": "0.1", "package_sha": "abc123"}

    monkeypatch.setattr(bridge, "run_bridge", run)
    return calls


def test_requires_lipid_assay(lipid_arf, fake_measure):
    lipid_arf.assay_kind = "metabolite"
    body = json.loads(tools.pathway_project(**_GROUPS))
    assert body["error"]["code"] == "missing_state"
    assert body["error"]["required_tools"] == ["record_objective"]


def test_unknown_matrix_id_is_missing_state_stat_matrix(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(**_GROUPS, matrix_id="nope"))
    assert body["error"]["code"] == "missing_state"
    assert body["error"]["state"] == "stat_matrix"
    assert body["error"]["required_tools"] == ["dataset_preprocess"]
    assert fake_measure == []


def test_requires_a_dataset(monkeypatch, fake_measure):
    s = session_state.AnalysisSession()
    s.assay_kind = "lipid"
    monkeypatch.setattr(session_state, "session", s)
    body = json.loads(tools.pathway_project(**_GROUPS))
    assert body["error"]["required_tools"] == ["load_dataset"]


def test_pathway_has_no_arf_source():
    """旧 ARF 経路（`source="arf"`）は消えた。`pathway_project` はデータセットだけを読む。"""
    import inspect
    assert "source" not in inspect.signature(tools.pathway_project).parameters
    assert not hasattr(tools, "_build_arf")
    assert not hasattr(tools, "_arf2_catalog")
    assert not hasattr(tools, "_pick_source")


def test_unknown_pathway_is_bad_request(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["nope"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_project_writes_run_and_returns_tables(lipid_arf, fake_measure, tmp_path):
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
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
    assert manifest["notes"]["adjacency"] == body["notes"]["adjacency"]
    assert manifest["notes"].get("handoff") == body["notes"].get("handoff")
    pairs = manifest["opposite_pairs"]["remodeling"]
    assert pairs["undeterminable"] == body["pathways"]["remodeling"]["undeterminable_pairs"]
    assert len(pairs["rows"]) == len(body["pathways"]["remodeling"]["opposite_pairs_tsv"].splitlines()[1:])
    assert lipid_arf.pathway.last_run_id == run_id
    req = fake_measure[0]
    assert req["mode"] == "measure" and req["networks"]["remodeling"] == "remodeling_lipidonly.cyjs"


def test_finite_turns_nonfinite_floats_into_none():
    out = tools._finite({"a": [float("nan"), 1.5, float("inf")], "b": {"c": float("-inf")}})
    assert out == {"a": [None, 1.5, None], "b": {"c": None}}
    json.dumps(out, allow_nan=False)


def _empty_scope_items(node, scope):
    """上流の空スコープ（memberFeature: []、値は全部 null）。"""
    sid = f"{B}n/{node}/scope/{scope}"
    chains = [] if scope == "whole" else [scope.removeprefix("acyl-")]
    return [{"@id": sid, "@type": "mslm:LipidMeasurementScope", "onPathwayNode": {"@id": f"{B}n/{node}"},
             "scopeType": "subclass" if scope == "whole" else "acyl_subset",
             "acylChainSignature": chains, "memberFeature": []},
            {"@id": sid + "/m1", "@type": "mslm:SampleMeasurement", "measuredEntity": {"@id": sid},
             "inSample": {"@id": f"{B}s/1"}, "value": None}]


def test_handoff_note_counts_dropped_rows_by_reason_and_kept_rows_by_route():
    note = tools._handoff_note({"dropped_unannotated": 0, "dropped_unparsable_name": 7, "dropped_no_formula": 120,
                                "kept_by_route": {"as_is": 200, "name": 60, "class_formula": 45}, "n_kept": 300})
    assert "127" in note and "組成式が無い 120" in note and "読めない 7" in note
    assert "200" in note and "60" in note and "45" in note and "300" in note
    assert "unannotated" not in note and "msp_then_lbm" not in note and "除かれ" not in note
    # 落とした行が無ければ注意は出さない
    zero = {f"dropped_{r}": 0 for r in ("unannotated", "unparsable_name", "isotope_label", "ambiguous_class",
                                        "ambiguous_ontology", "unmapped_ontology", "class_name_conflict",
                                        "no_formula", "formula_mismatch")}
    assert tools._handoff_note({**zero, "kept_by_route": {"as_is": 300, "name": 0, "class_formula": 0},
                                "n_kept": 300}) is None
    assert tools._handoff_note({}) is None


def test_handoff_note_names_the_new_reasons_and_excluded_spots():
    from metabolomix.pathway import handoff
    assert {f"dropped_{k}" for k in tools._DROP_LABELS} == {
        k for k in handoff._new_report("arf") if k.startswith("dropped_")}
    note = tools._handoff_note({"dropped_isotope_label": 3, "dropped_class_name_conflict": 2,
                                "kept_by_route": {"as_is": 10, "name": 0, "class_formula": 0}, "n_kept": 10})
    assert "計 5 件" in note and "同位体の標識 3" in note and "頭部" in note and "除外" not in note
    # ARF の手動除外スポットは脂質の理由とは別に述べる（それだけでも注意を出す）
    note = tools._handoff_note({"dropped_excluded_spots": 4,
                                "kept_by_route": {"as_is": 10, "name": 0, "class_formula": 0}, "n_kept": 10})
    assert "手動除外のスポット 4 件" in note and "計" not in note


def test_project_returns_handoff_counts_at_top_level(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
    for key in ("kept_by_route", "dropped_unannotated", "dropped_unparsable_name", "dropped_ambiguous_class",
                "dropped_unmapped_ontology", "dropped_ambiguous_ontology", "dropped_no_formula",
                "dropped_formula_mismatch", "dropped_isotope_label", "dropped_class_name_conflict", "n_kept"):
        assert key in body["handoff"]


def test_mapped_nodes_counts_nonempty_whole_over_lipid_nodes(lipid_arf, fake_measure, monkeypatch):
    # PE はノードがあるがデータに無い（空の whole スコープ）。分母は脂質ノード
    # （PC/PE/LPC/DG）だけで、PathwayEntity（PA）と接合ノードは数えない
    monkeypatch.setattr(sys.modules[__name__], "_EXTRA", _empty_scope_items("PE", "whole"))
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
    res = body["pathways"]["remodeling"]
    assert res["mapped_nodes"] == "1/4"
    assert {l.split("\t")[0] for l in res["nodes_tsv"].splitlines()[1:]} == {"PC"}   # PE は行にならない


def test_requested_chain_absent_everywhere_is_bad_request(lipid_arf, fake_measure, monkeypatch):
    # 22:6 は上流が空スコープとして出すだけ（どのパスウェイにも分子種が無い）
    monkeypatch.setattr(sys.modules[__name__], "_EXTRA", _empty_scope_items("PC", "acyl-22:6"))
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"],
                                            acyl_chains=["18:0", "22:6"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert body["error"]["details"]["missing_chains"] == ["22:6"]
    ok = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"],
                                          acyl_chains=["18:0"]))
    assert "run_id" in ok


def test_same_groups_are_bad_request(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(group_a=["a"], group_b=["b"],
                                            group_a_label="X", group_b_label="X"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert fake_measure == []


def test_bridge_error_becomes_envelope(lipid_arf, monkeypatch):
    def boom(*a, **k):
        raise bridge.BridgeError("MSLM_NOT_CONFIGURED", "設定して", {"setting": "x"})
    monkeypatch.setattr(bridge, "run_bridge", boom)
    body = json.loads(tools.pathway_project(**_GROUPS))
    assert body["error"]["code"] == "MSLM_NOT_CONFIGURED"


# --- pathway_plot ---

def _write_pdf(path):
    fitz = require_pymupdf()
    doc = fitz.open()
    doc.new_page(width=842, height=595).insert_text((72, 72), "PC")
    doc.save(path)


def test_plot_without_run_is_missing_state(lipid_arf):
    body = json.loads(tools.pathway_plot())
    assert body["error"]["required_tools"] == ["pathway_project"]


def _render_stub(monkeypatch, seen):
    def render(request, run_dir, *, timeout_s=300.0):
        seen.append(request)
        out = {}
        for p in request["pathways"]:
            pdf = run_dir / f"{p}.pdf"
            _write_pdf(pdf)
            out[p] = str(pdf)
        return {"status": "ok", "outputs": out, "mapped_nodes": {}}

    monkeypatch.setattr(bridge, "run_bridge", render)


def test_plot_returns_paths_and_summary_by_default(lipid_arf, fake_measure, monkeypatch):
    run_id = json.loads(tools.pathway_project(**_GROUPS,
                                              pathways=["remodeling", "ceramide"]))["run_id"]
    seen = []
    _render_stub(monkeypatch, seen)
    out = tools.pathway_plot(pathways=["ceramide"])
    assert isinstance(out, str)
    body = json.loads(out)
    manifest = json.loads((tools._run_dir(run_id, session_state.session) / "manifest.json").read_text(encoding="utf-8"))
    [figure] = body["figures"]
    assert figure["pathway"] == "ceramide"
    assert figure["pdf"].endswith("ceramide.pdf") and Path(figure["png"]).is_file()
    assert figure["png"].endswith("ceramide.png")
    assert figure["mapped_nodes"] == manifest["mapped_nodes"]["ceramide"]
    assert body["run_id"] == run_id and body["group_a"] == "C" and body["group_b"] == "P"
    assert body["note"] == tools._FIGURE_NOTE and len(out) < 2000
    assert seen[0]["mode"] == "render" and seen[0]["x_order"] == ["C", "P"]
    assert seen[0]["mztab"].endswith("input.mztab")
    assert run_id in seen[0]["out_dir"]


def test_plot_inline_image_keeps_the_old_payload(lipid_arf, fake_measure, monkeypatch):
    json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling", "ceramide"]))
    _render_stub(monkeypatch, [])
    out = tools.pathway_plot(pathways=["ceramide"], inline_image=True)
    assert isinstance(out, list) and isinstance(out[0], str) and isinstance(out[1], Image)
    assert len(out) == 2
    assert json.loads(out[0])["figures"][0]["pathway"] == "ceramide"


def test_plot_inline_image_must_be_bool(lipid_arf, fake_measure):
    json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
    body = json.loads(tools.pathway_plot(inline_image="yes"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_plot_pathway_not_in_run_is_bad_request(lipid_arf, fake_measure):
    json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
    body = json.loads(tools.pathway_plot(pathways=["global"]))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"


def test_render_failure_envelope(lipid_arf, fake_measure, monkeypatch):
    json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))

    def fail(*a, **k):
        raise bridge.BridgeError("PATHWAY_RENDER_FAILED", "Chromium がありません")
    monkeypatch.setattr(bridge, "run_bridge", fail)
    body = json.loads(tools.pathway_plot())
    assert body["error"]["code"] == "PATHWAY_RENDER_FAILED"


def test_unreadable_pdf_becomes_render_failed(lipid_arf, fake_measure, monkeypatch):
    """PDF が壊れている（0 ページ・読めない）なら例外を漏らさず PATHWAY_RENDER_FAILED。"""
    json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))

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
    run_id = json.loads(tools.pathway_project(**_GROUPS,
                                              pathways=["remodeling"]))["run_id"]
    (tools._run_dir(run_id, session_state.session) / "manifest.json").write_text("{not json", encoding="utf-8")
    body = json.loads(tools.pathway_plot(run_id=run_id))
    assert body["error"]["code"] == "MSLM_OUTPUT_INVALID"


# --- 統計エンジンの行列から作る（Task 8） ---------------------------------------------------

def _sml_row(run_dir, name):
    for line in (run_dir / "input.mztab").read_text(encoding="utf-8").splitlines():
        cols = line.split("\t")
        if cols[0] == "SML" and cols[7] == name:
            return cols
    raise AssertionError(f"{name} の SML 行がありません")


def test_pathway_project_default_matrix_and_missing_as_zero_in_mztab(lipid_arf, fake_measure, tmp_path):
    import numpy as np
    ds = lipid_arf.dataset
    ds.feature_matrix[0, ds.sample_names.index("C0_r0")] = np.nan      # f0（PC 34:1）の欠測
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"]))
    assert "error" not in body, body
    run_dir = tmp_path / "reports" / "pathway" / body["run_id"]
    row = _sml_row(run_dir, "PC 34:1")
    assert row[14] == "0.0"                      # 欠測は 0（MSLipidMapper の約束）。null ではない
    assert float(row[15]) > 0 and "null" not in row[14:18]
    # 既定の行列の ID が戻り値と manifest に残る（登録はしない）
    assert body["matrix_id"]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["matrix_id"] == body["matrix_id"] and manifest["test"] == "moderated_t"
    assert getattr(ds, "stat_matrices", None) in (None, {})
    # ノード表に新しい欄
    header = body["pathways"]["remodeling"]["nodes_tsv"].splitlines()[0].split("\t")
    assert header[-4:] == ["p", "ci_low", "ci_high", "verdict"]


def test_pathway_project_uses_registered_matrix_id(lipid_arf, fake_measure):
    from metabolomix.tools.stats_tools import dataset_preprocess
    pre = json.loads(dataset_preprocess(blank_min_ratio=None, detection_min_rate=None,
                                        max_qc_rsd=None, max_d_ratio=None))
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"], matrix_id=pre["matrix_id"]))
    assert "error" not in body, body
    assert body["matrix_id"] == pre["matrix_id"]


def test_samples_outside_the_matrix_are_rejected(lipid_arf, fake_measure):
    body = json.loads(tools.pathway_project(group_a=["C0_r0", "QC0"], group_b=["P0_r0", "P1_r0"],
                                            group_a_label="C", group_b_label="P"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert body["error"]["details"]["absent"] == ["QC0"]
    assert fake_measure == []


def test_qc_injections_never_reach_the_mztab(monkeypatch, tmp_path, fake_measure):
    from tests.test_stats_results import _ds
    s = session_state.AnalysisSession()
    monkeypatch.setattr(session_state, "session", s)
    s.assay_kind = "lipid"
    ds = _ds(groups=("C", "P"), per_group=2, reps=1, n_qc=4)
    ds.analysis_dir = str(tmp_path)
    ds.feature_annotations = {"f0": {"name": "PC 34:1", "ontology": "PC"},
                              "f1": {"name": "LPC 18:1", "ontology": "LPC"}}
    s.dataset = ds
    from metabolomix.tools.stats_tools import dataset_preprocess
    pre = json.loads(dataset_preprocess(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                                        max_qc_rsd=None, max_d_ratio=None, normalize="none"))
    assert pre["n_qc"] == 4
    body = json.loads(tools.pathway_project(**_GROUPS, pathways=["remodeling"], matrix_id=pre["matrix_id"]))
    assert "error" not in body, body
    text = (tmp_path / "reports" / "pathway" / body["run_id"] / "input.mztab").read_text(encoding="utf-8")
    assert "QC" not in text
    assert text.count("assay[") >= 4


def test_group_with_too_few_units_after_collapse_is_bad_request(monkeypatch, tmp_path, fake_measure):
    from tests.test_stats_results import _ds
    s = session_state.AnalysisSession()
    monkeypatch.setattr(session_state, "session", s)
    s.assay_kind = "lipid"
    ds = _ds(groups=("C", "P"), per_group=1, reps=2, n_qc=0)     # 各群 1 個体 × 2 反復
    ds.analysis_dir = str(tmp_path)
    ds.feature_annotations = {"f0": {"name": "PC 34:1", "ontology": "PC"},
                              "f1": {"name": "LPC 18:1", "ontology": "LPC"}}
    s.dataset = ds
    body = json.loads(tools.pathway_project(group_a=["C0_r0", "C0_r1"], group_b=["P0_r0", "P0_r1"],
                                            group_a_label="C", group_b_label="P"))
    assert body["error"]["code"] == "PATHWAY_BAD_REQUEST"
    assert body["error"]["details"]["n_units"] == 1
    assert fake_measure == []
