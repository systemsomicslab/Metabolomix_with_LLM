"""dataset_plot（volcano・pca・loadings・eic・eic_compounds）と図の置き場（B-2b Task 6）。

純粋な組み立て（`plots/*`）は各々のテストが固定している。ここでは `DatasetState` を入力に、
ツール層の仕事——結果の選び方（古い結果は描かない）・別の kind の引数の扱い（断らずに notes）・`.EIC.aef` の出所・
形式名・`save_figure` が `ds.figures` から描くこと——を見る。fixture はテスト自身が作る。
"""
import json

import numpy as np
import pytest
from mcp.server.fastmcp import Image

from metabolomix.analysis import sample_manifest
from metabolomix.core import handles, session_state
from metabolomix.mztab.dataset_state import DatasetState
from metabolomix.tools import plot_tools, reports, sample_tools
from metabolomix.tools.stats_tools import dataset_pca
from tests.config_fixtures import use_config
from tests.curation_fixtures import write_alignment_set
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata

PNG_MAGIC = b"\x89PNG"


@pytest.fixture(autouse=True)
def fresh_state(tmp_path, monkeypatch):
    use_config(monkeypatch, tmp_path, {"paths.reports_dir": tmp_path / "reports"})
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def _synthetic(n_features=12, n_samples=8) -> DatasetState:
    """ctrl_* と treat_* が半々の合成データセット（特徴にクラス `ontology` を持たせる）。"""
    rng = np.random.default_rng(0)
    ds = DatasetState()
    ds.feature_matrix = rng.random((n_features, n_samples)) * 1000.0 + 50.0
    half = n_samples // 2
    ds.sample_names = ([f"ctrl_{i}" for i in range(half)]
                       + [f"treat_{i}" for i in range(n_samples - half)])
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(n_samples)]
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    ds.feature_metadata = {
        f"f{i}": {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                  "ontology": "PC" if i % 2 == 0 else "TG", "adduct": "[M+H]+"}
        for i in range(n_features)}
    ds.quantification_measure = "peak_height"
    ds.validation_result = {"ok": True, "errors": [], "warnings": []}
    session_state.session.dataset = ds
    return ds


def _set_metadata(ds, *, group=None, batch=None):
    rows = sample_manifest.resolve_metadata(ds, None)
    for row in rows:
        if group:
            row["group"] = group(row["sample_id"])
        if batch:
            row["batch"] = batch(row["sample_id"])
    sample_manifest.apply_defaults_and_inputs(ds, rows)


def _analysed() -> DatasetState:
    """統計エンジンの前処理・PCA・2 群比較を済ませたデータセット（`_analysed_engine` の ds）。"""
    return _analysed_engine()[0]


def _analysed_engine():
    """統計エンジンの前処理・PCA・2 群比較（`dataset_compare`）で結果を作る。MCP から作れる状態はこちら。

    戻り値は (ds, pca_id, compare_id)。
    """
    from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess

    ds = _synthetic()
    _set_metadata(ds, group=lambda n: n.split("_")[0],
                  batch=lambda n: "b1" if n.endswith(("0", "1")) else "b2")
    prepared = json.loads(dataset_preprocess())
    assert prepared["status"] == "success", prepared
    pca = json.loads(dataset_pca(n_components=3, matrix_id=prepared["matrix_id"]))
    assert pca["status"] == "success", pca
    compared = json.loads(dataset_compare(reference_group="ctrl", test_group="treat",
                                          matrix_id=prepared["matrix_id"], allow_confounded=True))
    assert compared["status"] == "success", compared
    return ds, pca["result_id"], compared["result_id"]


def _payload(**kwargs) -> dict:
    out = plot_tools.dataset_plot(output="payload", **kwargs)
    assert isinstance(out, str), out
    return json.loads(out)


def _folder_handle(tmp_path, n_files=3):
    """`.arf` `.arf2` `.EIC.aef` を持つ解析フォルダを `load_dataset` で開く。"""
    from metabolomix.tools import handle_tools

    paths = write_alignment_set(tmp_path / "neg", n_files=n_files)
    write_mddata(paths["arf2"].with_suffix(".mddata"), [
        analysis_file_row(i, f"s{i}", "g1" if i < 2 else "g2", file_type=0, order=None, batch=1)
        for i in range(n_files)])
    body = json.loads(handle_tools.load_dataset(directory=str(paths["arf2"].parent)))
    assert body["status"] == "success", body
    return handles.registry.get(body["dataset_handle"])


# ---------- volcano / pca / loadings ----------

def test_volcano_from_a_two_group_result():
    _analysed()
    out = plot_tools.dataset_plot(kind="volcano")
    assert isinstance(out, list) and len(out) == 2
    caption, image = out
    assert isinstance(caption, str) and "ctrl" in caption and "n=4" in caption
    assert isinstance(image, Image) and image.data[:4] == PNG_MAGIC
    payload = _payload(kind="volcano")
    assert payload["plot_schema"] == "lipidmix.volcano.v1"
    assert payload["comparison"]["group_a"] == "ctrl"


def test_pca_color_by_batch():
    _analysed()
    batch = _payload(kind="pca", color_by="batch")
    assert batch["color_by"] == "batch"
    assert {p["group"] for p in batch["points"]} == {"b1", "b2"}
    group = _payload(kind="pca")                       # 既定は group
    assert group["color_by"] == "group"
    assert {p["group"] for p in group["points"]} == {"ctrl", "treat"}
    role = _payload(kind="pca", color_by="role")
    assert {p["group"] for p in role["points"]} == {"sample"}
    image = plot_tools.dataset_plot(kind="pca", color_by="batch")
    assert isinstance(image, list) and image[1].data[:4] == PNG_MAGIC
    bad = json.loads(plot_tools.dataset_plot(kind="pca", color_by="colour"))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST" and "group" in bad["error"]["message"]


def test_loadings_bar_colour_uses_ontology():
    ds = _analysed()
    payload = _payload(kind="loadings", pcs=[1, 2], top_n=3)
    rows = [row for panel in payload["panels"] for row in panel["rows"]]
    assert rows and all(row["ontology"] == ds.feature_metadata[row["feature_id"]]["ontology"]
                        for row in rows)
    assert {row["ontology"] for row in rows} <= {"PC", "TG"}
    image = plot_tools.dataset_plot(kind="loadings", top_n=3)       # 描画まで通る
    assert isinstance(image, list) and image[1].data[:4] == PNG_MAGIC


# ---------- eic / eic_compounds ----------

def test_eic_from_artifact_paths(tmp_path):
    handle = _folder_handle(tmp_path)
    payload = _payload(kind="eic", feature_id="0", dataset_handle=handle.handle_id)
    assert payload["plot_schema"] == "lipidmix.eic.v1"
    assert payload["source"]["file_name"] == "AlignmentResult_x.EIC.aef"
    assert payload["spot"]["spot_id"] == 0 and len(payload["series"]) == 3
    narrowed = _payload(kind="eic", feature_id="1", file_ids=[0, 2], dataset_handle=handle.handle_id)
    assert [s["file_id"] for s in narrowed["series"]] == [0, 2]
    image = plot_tools.dataset_plot(kind="eic", feature_id="0", dataset_handle=handle.handle_id)
    assert isinstance(image, list) and image[1].data[:4] == PNG_MAGIC
    unknown = json.loads(plot_tools.dataset_plot(kind="eic", feature_id="99",
                                                 dataset_handle=handle.handle_id))
    assert unknown["error"]["code"] == "DATASET_BAD_REQUEST"
    missing = json.loads(plot_tools.dataset_plot(kind="eic", dataset_handle=handle.handle_id))
    assert missing["error"]["code"] == "DATASET_BAD_REQUEST" and "feature_id" in missing["error"]["message"]


def test_eic_compounds_by_ontology(tmp_path):
    handle = _folder_handle(tmp_path)
    payload = _payload(kind="eic_compounds", file_id=0, ontologies=["PC"],
                       dataset_handle=handle.handle_id)
    assert payload["plot_schema"] == "lipidmix.eic.multi.v1"
    assert payload["sample"]["file_id"] == 0
    assert {s["spot_id"] for s in payload["series"]} == {0, 1, 2}
    by_name = _payload(kind="eic_compounds", file_id=1, names=["PC 34"],
                       dataset_handle=handle.handle_id)
    assert [s["spot_id"] for s in by_name["series"]] == [0]
    image = plot_tools.dataset_plot(kind="eic_compounds", file_id=0, ontologies=["PC"],
                                    dataset_handle=handle.handle_id)
    assert isinstance(image, list) and image[1].data[:4] == PNG_MAGIC
    nothing = json.loads(plot_tools.dataset_plot(kind="eic_compounds", file_id=0, names=["zzz"],
                                                 dataset_handle=handle.handle_id))
    assert nothing["error"]["code"] == "DATASET_BAD_REQUEST"


def test_eic_default_traces_skip_excluded_samples(tmp_path):
    handle = _folder_handle(tmp_path)
    handle.dataset.excluded_samples = {handle.dataset.sample_names[1]: "failed injection"}
    payload = _payload(kind="eic", feature_id="0", dataset_handle=handle.handle_id)
    assert [s["file_id"] for s in payload["series"]] == [0, 2]


def test_eic_without_artifact_is_explained():
    _synthetic()                                  # mzTab-M 相当: .EIC.aef を持たない
    for kwargs in ({"kind": "eic", "feature_id": "f0"},
                   {"kind": "eic_compounds", "file_id": 0, "names": ["Compound"]}):
        body = json.loads(plot_tools.dataset_plot(**kwargs))
        assert body["error"]["code"] == "COMPANION_ARTIFACT_MISSING"
        assert ".EIC.aef" in body["error"]["message"]
        assert "load_dataset" in body["error"]["message"]


# ---------- 形式名・引数の検査 ----------

def test_payload_format_names_are_kept(tmp_path):
    _analysed()
    assert _payload(kind="volcano")["plot_schema"] == "lipidmix.volcano.v1"
    assert _payload(kind="pca")["plot_schema"] == "lipidmix.pca_scores.v1"
    assert _payload(kind="loadings")["plot_schema"] == "lipidmix.pca_loadings.v1"
    handle = _folder_handle(tmp_path)
    assert _payload(kind="eic", feature_id="0",
                    dataset_handle=handle.handle_id)["plot_schema"] == "lipidmix.eic.v1"
    assert _payload(kind="eic_compounds", file_id=0, ontologies=["PC"],
                    dataset_handle=handle.handle_id)["plot_schema"] == "lipidmix.eic.multi.v1"
    # 項目を選んで群ごとに描く 3 つ（Task 7）。g1 = s0・s1、g2 = s2
    handle_id = handle.handle_id
    assert _payload(kind="group_intensity", items=["PC"], groups=["g1", "g2"],
                    dataset_handle=handle_id)["plot_schema"] == "lipidmix.group_intensity.v1"
    assert _payload(kind="species", items=["PC"], groups=["g1", "g2"],
                    dataset_handle=handle_id)["plot_schema"] == "lipidmix.species_intensity.v1"
    assert _payload(kind="species_pca", items=["PC"], groups=["g1", "g2"],
                    dataset_handle=handle_id)["plot_schema"] == "lipidmix.species_pca.v1"


def test_irrelevant_argument_is_ignored_with_a_note():
    """別の kind の引数は断らずに使わず、図を描いて説明の `notes:` 行（payload は `notes`）にどの kind の引数かを
    書く（save_figure と同じ扱い。レビュー V4）。図はその引数を渡さない呼び出しと同じもの。未知の kind は誤り。"""
    ds = _analysed()
    body = plot_tools.dataset_plot(kind="volcano", top_n=5)
    assert isinstance(body, list) and isinstance(body[1], Image)
    assert "\nnotes: top_n は kind=loadings の引数なので使っていない（kind='volcano'）。" in body[0]
    caption = plot_tools.dataset_plot(kind="pca", feature_id="1")[0]
    assert "notes: feature_id は kind=eic の引数" in caption
    caption = plot_tools.dataset_plot(kind="loadings", color_by="batch")[0]
    assert "notes: color_by は kind=pca の引数" in caption
    # payload は notes 欄（形式名は変えない）。図は引数を渡さない呼び出しと同じ
    plain = _payload(kind="pca")
    noted = _payload(kind="pca", top_n=5, contrast_id="treat_vs_ctrl")
    assert noted.pop("notes") == ["contrast_id は kind=volcano の引数なので使っていない（kind='pca'）。",
                                  "top_n は kind=loadings の引数なので使っていない（kind='pca'）。"]
    assert noted == plain and noted["plot_schema"] == "lipidmix.pca_scores.v1"
    assert "notes" not in plain
    assert ds.figures["pca"]["payload"] == plain       # 置き場の図にも notes を混ぜない
    # 描けない呼び出しの誤りにも、使わなかった引数を添える（ここは .EIC.aef が無い）
    err = json.loads(plot_tools.dataset_plot(kind="eic", result_id="res_x", feature_id="0"))["error"]
    assert err["code"] == "COMPANION_ARTIFACT_MISSING"
    assert err["details"]["notes"] == [
        "result_id は kind=volcano / pca / loadings の引数なので使っていない（kind='eic'）。"]
    err = json.loads(plot_tools.dataset_plot(kind="heatmap"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    assert all(kind in err["message"] for kind in plot_tools._KINDS)


def test_every_error_envelope_carries_the_unused_argument_notes(tmp_path, monkeypatch):
    """kind が分かれば、描けなかったどの誤りの封筒（引数の形の誤り・ハンドルの誤り・データセットや結果の無い
    `missing_state`・書けない `FIGURE_SAVE_FAILED`・図の組み立ての誤り）も `details.notes` に使わなかった引数を書く
    （直して呼び直す前に、その引数が効かないと分かる）。ほかの欄は注記の無い呼び出しと同じ。"""
    top_n = ["top_n は kind=loadings の引数なので使っていない（kind='pca'）。"]
    items = ["items は kind=group_intensity / species / species_pca の引数なので使っていない（kind='eic'）。"]

    def error(**kwargs):
        body = json.loads(plot_tools.dataset_plot(**kwargs))
        assert "error" in body, body
        return body["error"]

    # データセットが無い（missing_state も details.notes を持つ。state / required_tools はそのまま）
    err = error(kind="pca", top_n=5)
    assert err["code"] == "missing_state" and err["required_tools"] == ["load_dataset"]
    assert err["details"] == {"notes": top_n}
    plain = error(kind="pca")
    assert "details" not in plain and {k: v for k, v in err.items() if k != "details"} == plain

    _synthetic()                                     # 結果の無いデータセット
    err = error(kind="volcano", top_n=5)
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_compare"]
    assert err["details"]["notes"] == ["top_n は kind=loadings の引数なので使っていない（kind='volcano'）。"]
    for kwargs, field in (({"output": "svg"}, "output"), ({"output_path": "x.svg"}, "output_path")):
        err = error(kind="pca", top_n=5, **kwargs)
        assert err["code"] == "DATASET_BAD_REQUEST"
        assert err["details"] == {field: kwargs[field], "notes": top_n}
    err = error(kind="eic", feature_id="0", items=["PC"], output_path=str(tmp_path / "fig.svg"))
    assert err["code"] == "DATASET_BAD_REQUEST" and err["details"]["notes"] == items
    err = error(kind="pca", top_n=5, dataset_handle="dh_nope")
    assert err["code"] == "DATASET_HANDLE_NOT_FOUND"
    assert err["details"]["notes"] == top_n and err["details"]["dataset_handle"] == "dh_nope"

    _analysed()
    blocker = tmp_path / "file"
    blocker.write_text("x", encoding="utf-8")         # 親のフォルダがファイルなので書けない
    err = error(kind="pca", top_n=5, output_path=str(blocker / "sub" / "fig.png"))
    assert err["code"] == "FIGURE_SAVE_FAILED" and err["details"]["notes"] == top_n
    # 図の組み立ての誤り（DomainError・ValueError）は今までどおり
    err = error(kind="pca", top_n=5, result_id="res_nope")
    assert err["code"] == "DATASET_BAD_REQUEST" and err["details"]["notes"] == top_n


def test_plot_output_env_is_honoured(monkeypatch):
    _analysed()
    monkeypatch.setenv("LIPIDMIX_PLOT_OUTPUT", "payload")
    out = plot_tools.dataset_plot(kind="volcano")
    assert isinstance(out, str) and json.loads(out)["plot_schema"] == "lipidmix.volcano.v1"
    bad = json.loads(plot_tools.dataset_plot(kind="volcano", output="svg"))
    assert bad["error"]["code"] == "DATASET_BAD_REQUEST"


def test_dataset_plot_is_registered_with_a_dataset_handle():
    import inspect

    import server

    assert "dataset_plot" in {t.name for t in server.mcp._tool_manager.list_tools()}
    params = inspect.signature(plot_tools.dataset_plot).parameters
    for name in ("kind", "dataset_handle", "result_id", "output", "title", "color_by", "pcs",
                 "top_n", "value", "feature_id", "file_ids", "file_id", "names", "ontologies"):
        assert name in params


# ---------- 古い結果と図の置き場 ----------

def test_stale_result_not_drawn_as_current():
    ds, pca_id, diff_id = _analysed_engine()
    assert _payload(kind="pca", result_id=pca_id)["plot_schema"] == "lipidmix.pca_scores.v1"
    assert ds.figures["pca"]["result_id"] == pca_id

    out = json.loads(sample_tools.dataset_exclude(samples=["treat_3"], reason="outlier"))
    assert out["changed"] is True

    for kind, rid in (("pca", pca_id), ("loadings", pca_id), ("volcano", diff_id)):
        err = json.loads(plot_tools.dataset_plot(kind=kind, result_id=rid))["error"]
        assert err["code"] == "STALE_ANALYSIS_RESULT", kind
        assert err["details"]["result_id"] == rid
        assert kind not in ds.figures                 # 古い図は置き場に残らない
    # result_id を省いたときは「直近」を使う。除外で捨てられたので missing_state（作り直しを案内）
    for kind, tool in (("pca", "dataset_pca"), ("loadings", "dataset_pca"),
                       ("volcano", "dataset_compare")):
        err = json.loads(plot_tools.dataset_plot(kind=kind))["error"]
        assert err["code"] == "missing_state" and tool in err["required_tools"], kind
    # 古い図は save_figure からも保存できない。除外が捨てた図は「描いていない」ではなく STALE_FIGURE（理由つき）
    err = json.loads(reports.save_figure("pca", "stale-1"))["error"]
    assert err["code"] == "STALE_FIGURE" and err["details"]["required_tools"] == ["dataset_plot"]
    assert err["details"]["reason"] == ["exclusion"] and "dataset_exclude" in err["message"]


def test_metadata_change_drops_group_dependent_figure_only():
    ds = _analysed()
    _payload(kind="pca")
    _payload(kind="volcano")
    _payload(kind="loadings")
    assert {"pca", "volcano", "loadings"} <= set(ds.figures)
    from metabolomix.analysis.result_state import invalidate_results

    invalidate_results(ds, {"group"})                # 差次的解析だけが古くなる
    assert "volcano" not in ds.figures
    assert {"pca", "loadings"} <= set(ds.figures)
    invalidate_results(ds, {"role"})                 # 前処理から先が古くなる
    assert not ({"pca", "loadings", "volcano"} & set(ds.figures))


def test_group_change_makes_the_pca_figure_refuse_to_save(tmp_path):
    """実験情報シートの適用（`dataset_set_sample_metadata`）の後は、古い色分けの pca 図を保存できない。"""
    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata

    ds = _analysed()
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    ds.assay_sources = {}
    for i, name in enumerate(ds.sample_names):
        (raw_dir / f"{name}.wiff").write_text("x", encoding="utf-8")
        ds.assay_sources[f"abundance_assay[{i + 1}]"] = str(raw_dir / f"{name}.wiff")
    plot_tools.dataset_plot(kind="pca")
    assert "pca" in ds.figures
    rows = ds.sample_metadata_rows
    lines = ["# schema = sample-manifest.v1",
             "sample_id	source_file	role	group	batch	injection_order	qc_pool	include"]
    for row in rows:                                  # 群だけを付け替える（前処理の入力は同じ）
        lines.append("	".join([row["sample_id"], f"{row['sample_id']}.wiff", row["role"], "NEW",
                                 str(row["batch"]), "" if row["injection_order"] is None else str(row["injection_order"]),
                                 "", "true"]))
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert json.loads(dataset_set_sample_metadata(str(manifest)))["status"] == "success"

    err = json.loads(reports.save_figure("pca", "meta-1"))["error"]
    # 試料の表の更新が図ごと捨てた図も、描いたことの無い図（missing_state）と見分けて STALE_FIGURE
    assert err["code"] == "STALE_FIGURE" and err["details"]["required_tools"] == ["dataset_plot"]
    assert "group" in err["details"]["reason"] and "試料の表の更新" in err["message"]
    assert not (tmp_path / "reports" / "figures" / "meta-1_pca.png").exists()
    dataset_pca(n_components=3)
    payload = _payload(kind="pca")                    # 解析と図をやり直せば新しい群で保存できる
    assert {p["group"] for p in payload["points"]} == {"NEW"}
    assert "meta-1_pca.png" in reports.save_figure("pca", "meta-1")


def test_figures_discarded_by_exclusion_are_stale_not_missing(tmp_path):
    """`dataset_exclude` が捨てた図（pca / loadings / volcano）は、save_figure で STALE_FIGURE と理由を返す
    （描いたことの無い図の missing_state と見分ける）。描き直せば印は消え、保存できる。"""
    ds = _analysed()
    for kind, args in (("pca", {}), ("loadings", {"top_n": 3}), ("volcano", {})):
        plot_tools.dataset_plot(kind=kind, **args)
    assert json.loads(sample_tools.dataset_exclude(features=["f0"], reason="noise"))["changed"] is True
    assert not ({"pca", "loadings", "volcano"} & set(ds.figures))
    figures = tmp_path / "reports" / "figures"
    for kind in ("pca", "loadings", "volcano"):
        err = json.loads(reports.save_figure(kind, "ex-1"))["error"]
        assert err["code"] == "STALE_FIGURE", kind
        assert err["details"] == {"kind": kind, "reason": ["exclusion"], "required_tools": ["dataset_plot"]}
        assert "除外" in err["message"] and f'dataset_plot(kind="{kind}")' in err["message"]
        assert not (figures / f"ex-1_{kind}.png").exists()
    # 描いたことの無い図は今までどおり missing_state
    err = json.loads(reports.save_figure("group_intensity", "ex-1"))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_plot"]

    # 解析と図をやり直せば印は消え、保存できる（ほかの捨てた図の印は残る）
    from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess

    prepared = json.loads(dataset_preprocess())
    dataset_pca(n_components=3, matrix_id=prepared["matrix_id"])
    plot_tools.dataset_plot(kind="pca")
    assert "pca" not in ds.discarded_figures and "volcano" in ds.discarded_figures
    assert "ex-1_pca.png" in reports.save_figure("pca", "ex-1")
    err = json.loads(reports.save_figure("volcano", "ex-1"))["error"]
    assert err["code"] == "STALE_FIGURE"
    dataset_compare(reference_group="ctrl", test_group="treat", matrix_id=prepared["matrix_id"],
                    allow_confounded=True)
    plot_tools.dataset_plot(kind="volcano")
    assert "ex-1_volcano.png" in reports.save_figure("volcano", "ex-1")
    assert ds.discarded_figures == {"loadings": ["exclusion"]}


def test_only_held_figures_get_the_discard_mark():
    """置かれていなかった kind には印を残さない（捨てたのではなく、描いていない）。群だけの変更は volcano だけ。"""
    from metabolomix.analysis.result_state import invalidate_results

    ds = _analysed()
    _payload(kind="pca")
    _payload(kind="volcano")
    invalidate_results(ds, {"group"})
    assert ds.discarded_figures == {"volcano": ["group"]} and "pca" in ds.figures
    invalidate_results(ds, {"role", "batch"})
    assert ds.discarded_figures == {"volcano": ["group"], "pca": ["batch", "role"]}
    err = json.loads(reports.save_figure("pca", "m-1"))["error"]
    assert err["code"] == "STALE_FIGURE" and err["details"]["reason"] == ["batch", "role"]
    assert "batch, role" in err["message"]


def test_factor_only_change_names_only_the_changed_column():
    """足した因子の列だけが変わった試料の表の更新は、群と同じく差次的解析を古くする（volcano を捨てる）が、
    捨てた理由に挙げるのは実際に変わった列だけ（群・biological_sample_id は変わっていない。レビュー V2）。"""
    ds = _analysed()
    _payload(kind="pca")
    _payload(kind="volcano")
    rows = [dict(row, dose="hi" if i % 2 else "lo") for i, row in enumerate(ds.sample_metadata_rows)]
    out = sample_manifest.apply_metadata(ds, rows)
    assert out["changed_fields"] == ["dose"]
    assert "volcano" not in ds.figures and "pca" in ds.figures     # 捨てる範囲は群の変更と同じ
    assert ds.discarded_figures == {"volcano": ["dose"]}
    err = json.loads(reports.save_figure("volcano", "f-1"))["error"]
    assert err["code"] == "STALE_FIGURE"
    assert err["details"] == {"kind": "volcano", "reason": ["dose"], "required_tools": ["dataset_plot"]}
    assert "変わった列: dose）" in err["message"] and "group" not in err["message"]


def test_pca_figure_with_baked_metadata_is_stale_when_the_revision_moves(tmp_path):
    """色に焼き込んだ試料の情報が更新されたら、図は古い（統計エンジンの PCA の結果も版が動けば古い）。"""
    from metabolomix.analysis.result_state import is_current

    ds = _analysed()
    plot_tools.dataset_plot(kind="pca")
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    for row in ds.sample_metadata_rows:
        row["group"] = "NEW"
    ds.metadata_revision += 1                         # 情報の更新（apply_metadata が進める版）
    assert not is_current(ds, ds.last_pca) and "pca" in ds.figures   # 統計エンジンの PCA は版が動けば古い
    err = json.loads(reports.save_figure("pca", "rev-1"))["error"]
    assert err["code"] == "STALE_FIGURE" and err["details"]["required_tools"] == ["dataset_plot"]
    assert "作り直" in err["message"] and err["details"]["kind"] == "pca"
    assert not (tmp_path / "reports" / "figures" / "rev-1_pca.png").exists()
    # 試料の情報を焼き込まない図（loadings）も、元の結果が古ければ保存できない（統計エンジンの PCA は版に縛られる）
    err = json.loads(reports.save_figure("loadings", "rev-1"))["error"]
    assert err["code"] == "STALE_FIGURE" and err["details"]["required_tools"] == ["dataset_plot"]
    assert not (tmp_path / "reports" / "figures" / "rev-1_loadings.png").exists()


def test_figures_are_not_part_of_the_fingerprint():
    from metabolomix.analysis.result_state import dataset_fingerprint

    ds = _analysed()
    before = dataset_fingerprint(ds)
    _payload(kind="pca")
    assert ds.figures and dataset_fingerprint(ds) == before


def test_stale_pca_result_after_metadata_change_is_refused():
    ds = _analysed()
    old = ds.last_pca["provenance"]["result_id"]
    ds.metadata_revision += 1                   # 試料の表が変わった（統計エンジンの PCA は群も設計に含む）
    err = json.loads(plot_tools.dataset_plot(kind="pca", result_id=old))["error"]
    assert err["code"] == "STALE_ANALYSIS_RESULT"


def test_result_id_must_match_the_kind():
    ds = _analysed()
    diff_id = ds.last_differential["provenance"]["result_id"]
    err = json.loads(plot_tools.dataset_plot(kind="pca", result_id=diff_id))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    err = json.loads(plot_tools.dataset_plot(kind="pca", result_id="res_nope"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"


def test_no_dataset_is_missing_state():
    err = json.loads(plot_tools.dataset_plot(kind="volcano"))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["load_dataset"]


# ---------- save_figure ----------

def test_save_figure_reads_dataset_figures(tmp_path):
    ds, pca_id, _diff_id = _analysed_engine()
    figures = tmp_path / "reports" / "figures"
    # dataset_plot を呼ぶ前は図が無い
    err = json.loads(reports.save_figure("pca", "a-1"))["error"]
    assert err["code"] == "missing_state" and "dataset_plot" in err["required_tools"]

    for kind in ("pca", "volcano"):
        plot_tools.dataset_plot(kind=kind)
        msg = reports.save_figure(kind, "a-1")
        assert (figures / f"a-1_{kind}.png").is_file() and f"figures/a-1_{kind}.png" in msg
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    msg = reports.save_figure("loadings", "a-1")
    assert (figures / "a-1_loadings.png").is_file() and (figures / "a-1_loadings.svg").is_file()
    assert "figures/a-1_loadings.png" in msg
    assert set(ds.figures) >= {"pca", "volcano", "loadings"}

    # 登録された結果を名指しして描くこともできる（dataset_plot を呼ばなくてよい）
    ds.figures.clear()
    msg = reports.save_figure("pca", "by-id", result_id=pca_id)
    assert (figures / "by-id_pca.png").is_file() and pca_id in msg


def test_save_figure_draws_the_eic_figure(tmp_path):
    handle = _folder_handle(tmp_path)
    plot_tools.dataset_plot(kind="eic_compounds", file_id=0, ontologies=["PC"],
                            dataset_handle=handle.handle_id)
    msg = reports.save_figure("eic_compounds", "e-1", dataset_handle=handle.handle_id)
    assert (tmp_path / "reports" / "figures" / "e-1_eic_compounds.png").is_file()
    assert "e-1_eic_compounds.png" in msg


def test_save_figure_rejects_old_arguments_and_kinds():
    import inspect

    params = inspect.signature(reports.save_figure).parameters
    assert list(params) == ["kind", "analysis_id", "title", "result_id", "dataset_handle", "contrast_id"]
    out = json.loads(reports.save_figure("heatmap", "x"))
    assert out["status"] == "error" and "pca" in out["message"]
    assert "pca_loadings" not in out["message"]
    # 使わない引数は断らない（データセットが無いので、図の前提の missing_state まで進む）
    err = json.loads(reports.save_figure("eic", "x", result_id="r1"))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["load_dataset"]
    err = json.loads(reports.save_figure("pca", "x", contrast_id="T_vs_C"))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["load_dataset"]
    # 説明文に「明示的なユーザー要求時だけ」の縛りは無い（報告書に埋め込む図を保存する）
    assert "明示的" not in reports.save_figure.__doc__


def test_save_figure_ignores_arguments_its_kind_does_not_use(tmp_path):
    """kind が使わない `result_id` / `contrast_id` は断らずに使わず、図を保存して `notes:` 行に書く。"""
    ds, pca_id, diff_id = _analysed_engine()
    figures = tmp_path / "reports" / "figures"
    plot_tools.dataset_plot(kind="pca")
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    msg = reports.save_figure("loadings", "n-1", contrast_id="treat_vs_ctrl")
    assert (figures / "n-1_loadings.png").is_file() and (figures / "n-1_loadings.svg").is_file()
    assert "notes:" in msg and "contrast_id='treat_vs_ctrl'" in msg and "volcano" in msg
    # result_id は pca に効く（名指しの結果から描く）。contrast_id だけを使わない
    msg = reports.save_figure("pca", "n-2", result_id=pca_id, contrast_id="treat_vs_ctrl")
    assert (figures / "n-2_pca.png").is_file() and pca_id in msg and "contrast_id" in msg
    # 効く引数だけなら notes は付かない
    assert "notes:" not in reports.save_figure("pca", "n-4", result_id=pca_id)
    # result_id を使わない kind（eic 系・group_intensity・species 系）では置き場の図を保存する
    handle = _folder_handle(tmp_path / "f")
    plot_tools.dataset_plot(kind="eic", feature_id="0", dataset_handle=handle.handle_id)
    msg = reports.save_figure("eic", "n-3", result_id=diff_id, dataset_handle=handle.handle_id)
    assert "n-3_eic.png" in msg and "notes:" in msg and f"result_id='{diff_id}'" in msg


def test_save_figure_reads_the_chosen_handle(tmp_path):
    first = _folder_handle(tmp_path / "one")
    second = _folder_handle(tmp_path / "two")                 # 直近
    plot_tools.dataset_plot(kind="eic", feature_id="0", dataset_handle=first.handle_id)
    err = json.loads(reports.save_figure("eic", "x"))["error"]          # 直近（second）には図が無い
    assert err["code"] == "missing_state"
    msg = reports.save_figure("eic", "x", dataset_handle=first.handle_id)
    assert "x_eic.png" in msg
    missing = json.loads(reports.save_figure("eic", "x", dataset_handle="dh_nope"))
    assert missing["error"]["code"] == "DATASET_HANDLE_NOT_FOUND"


def test_derived_handle_has_no_figures_or_discard_marks_of_its_parent(tmp_path):
    """定量の派生（`dataset_quant_apply` などの `derive`）は親の図も捨てた図の印も持ち越さない。

    親の図は親の行列から描いたもの（正規化の前の値）で、派生の図ではない。派生でまだ描いていない図は、親が
    持つ図（保存すれば親の値が派生の図として出る）でも、親の除外が捨てた図（STALE_FIGURE・理由 exclusion）でもなく
    missing_state。親の図と印はそのまま。"""
    from metabolomix.quant import tools as quant_tools

    parent = _folder_handle(tmp_path / "f")
    pid = parent.handle_id
    _payload(kind="species_pca", items=["PC"], groups=["g1", "g2"], dataset_handle=pid)
    assert json.loads(sample_tools.dataset_exclude(features=["2"], reason="noise",
                                                   dataset_handle=pid))["changed"] is True
    assert parent.dataset.discarded_figures == {"species_pca": ["exclusion"]}
    _payload(kind="group_intensity", items=["PC"], groups=["g1", "g2"], dataset_handle=pid)
    assert "group_intensity" in parent.dataset.figures

    out = json.loads(quant_tools.dataset_quant_apply(skip=True, dataset_handle=pid))
    assert out["status"] == "skipped", out
    did = out["dataset_handle"]
    derived = handles.registry.get(did).dataset
    assert derived.figures == {} and derived.discarded_figures == {}
    figures = tmp_path / "reports" / "figures"
    for kind in ("group_intensity", "species_pca"):
        err = json.loads(reports.save_figure(kind, "d-1", dataset_handle=did))["error"]
        assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_plot"], (kind, err)
        assert not (figures / f"d-1_{kind}.png").exists()
    # 親の図と印は変わらない
    assert "group_intensity" in parent.dataset.figures
    assert parent.dataset.discarded_figures == {"species_pca": ["exclusion"]}
    assert "p-1_group_intensity.png" in reports.save_figure("group_intensity", "p-1", dataset_handle=pid)
    # 派生で描けば派生の図として保存できる
    _payload(kind="group_intensity", items=["PC"], groups=["g1", "g2"], dataset_handle=did)
    assert "d-1_group_intensity.png" in reports.save_figure("group_intensity", "d-1", dataset_handle=did)


# ---------- save_figure の振る舞い ----------

def test_save_figure_title_override_reaches_the_render(monkeypatch):
    ds = _analysed()
    plot_tools.dataset_plot(kind="loadings", top_n=3, title="drawn title")
    seen = []
    from metabolomix.plots import pca_loadings as loadings_plot

    real = loadings_plot.render_loadings_plot
    monkeypatch.setattr(loadings_plot, "render_loadings_plot",
                        lambda payload, title=None: seen.append(title) or real(payload, title=title))
    reports.save_figure("loadings", "t-1", title="saved title")
    reports.save_figure("loadings", "t-2")                      # 省けば図を作ったときの題
    assert seen == ["saved title", "drawn title"]
    assert ds.figures["loadings"]["title"] == "drawn title"


def test_save_figure_path_traversal_stays_in_figures_dir(tmp_path):
    from metabolomix.tools import notes

    _analysed()
    plot_tools.dataset_plot(kind="pca")
    figures = (tmp_path / "reports" / "figures").resolve()
    for analysis_id in ("A/../x", "..\\x", "../../escape"):
        msg = reports.save_figure("pca", analysis_id)
        written = figures / f"{notes.make_slug(analysis_id)}_pca.png"
        assert written.is_file() and written.parent == figures, analysis_id
        assert str(tmp_path / "reports" / "figures" / written.name) in msg
    assert "/" not in notes.make_slug("A/../x") and "\\" not in notes.make_slug("..\\x")
    assert sorted(p.parent for p in (tmp_path / "reports").rglob("*.png")) == [figures] * 3
    assert not list(tmp_path.glob("*.png"))


def test_save_figure_with_unreadable_config_is_an_envelope(tmp_path, monkeypatch):
    from metabolomix.core import user_config

    _analysed()
    plot_tools.dataset_plot(kind="pca")
    broken = tmp_path / "broken" / "lipidmix.local.toml"
    broken.parent.mkdir()
    broken.write_text("[msdial]\nlbm = 'old.lbm2'\n", encoding="utf-8")        # 旧キー → CONFIG_INVALID
    monkeypatch.setenv(user_config.CONFIG_ENV, str(broken))
    user_config._cache.clear()
    err = json.loads(reports.save_figure("pca", "cfg-1"))["error"]
    assert err["code"] == "CONFIG_INVALID"
    assert not list(tmp_path.rglob("cfg-1_pca.png"))


def test_save_figure_twice_overwrites_and_succeeds(tmp_path):
    _analysed()
    figures = tmp_path / "reports" / "figures"
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    first = reports.save_figure("loadings", "twice")
    png, svg = figures / "twice_loadings.png", figures / "twice_loadings.svg"
    assert str(png) in first and png.is_file() and svg.is_file()
    png.write_bytes(b"old"), svg.write_bytes(b"old")
    second = reports.save_figure("loadings", "twice")
    assert str(png) in second and "保存" in second
    assert png.read_bytes()[:4] == PNG_MAGIC and svg.read_bytes() != b"old"


# ---------- feature_id の型（整数でも文字列でも同じ特徴に当たる） ----------

def _call_through_server(name: str, arguments: dict):
    """FastMCP の層（pydantic の引数検証）を通してツールを呼ぶ。"""
    import asyncio

    import server

    return asyncio.run(server.mcp.call_tool(name, arguments))


def _content_text(result) -> str:
    content = result[0] if isinstance(result, tuple) else result
    return "".join(c.text for c in content if getattr(c, "type", "") == "text")


def test_numeric_feature_id_is_accepted_as_int_or_str_through_fastmcp(tmp_path):
    handle = _folder_handle(tmp_path)
    got = {}
    for label, fid in (("int", 1), ("str", "1")):
        out = _call_through_server("dataset_plot", {
            "kind": "eic", "feature_id": fid, "output": "payload",
            "dataset_handle": handle.handle_id})
        got[label] = json.loads(_content_text(out))
    assert got["int"]["spot"]["spot_id"] == 1
    assert got["int"] == got["str"]


def test_feature_search_accepts_int_ids_through_fastmcp(tmp_path):
    handle = _folder_handle(tmp_path)
    by_type = {}
    for label, ids in (("int", [1]), ("str", ["1"])):
        out = _call_through_server("dataset_feature_search", {
            "feature_ids": ids, "dataset_handle": handle.handle_id})
        by_type[label] = json.loads(_content_text(out))
    assert by_type["int"]["features"] == by_type["str"]["features"]

# ---------- dataset_plot(output_path=...) ----------

def test_output_path_writes_the_drawn_figure(tmp_path):
    ds = _analysed()
    out = tmp_path / "out" / "deep" / "scores.png"
    body = plot_tools.dataset_plot(kind="pca", output_path=str(out))
    assert isinstance(body, list) and str(out) in body[0]
    assert out.is_file() and out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert ds.figures["pca"]["kind"] == "pca"
    # loadings は報告書用と同じく dpi 300 の PNG と同名の SVG
    load = tmp_path / "load.png"
    payload = json.loads(plot_tools.dataset_plot(kind="loadings", top_n=3, output="payload",
                                                 output_path=str(load)))
    assert payload["saved_files"] == [str(load), str(load.with_suffix(".svg"))]
    assert load.is_file() and load.with_suffix(".svg").is_file()
    assert payload["plot_schema"] == "lipidmix.pca_loadings.v1"
    # 同じパスは上書き
    assert plot_tools.dataset_plot(kind="pca", output_path=str(out))


def test_output_path_must_be_a_png_and_writable(tmp_path):
    ds = _analysed()
    err = json.loads(plot_tools.dataset_plot(kind="pca", output_path=str(tmp_path / "x.svg")))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST" and err["details"]["output_path"].endswith("x.svg")
    plot_tools.dataset_plot(kind="pca")
    before = ds.figures["pca"]
    blocker = tmp_path / "afile"
    blocker.write_text("x", encoding="utf-8")                  # 親のフォルダの位置にファイルがある
    err = json.loads(plot_tools.dataset_plot(kind="pca", output_path=str(blocker / "x.png")))["error"]
    assert err["code"] == "FIGURE_SAVE_FAILED"
    assert ds.figures["pca"] is before                         # 書けなければ図の置き場は変えない
    # output_path はどの kind の引数でもある（kind の引数の取り違えにしない）
    assert plot_tools.dataset_plot(kind="volcano", output_path=str(tmp_path / "v.png"))
    assert (tmp_path / "v.png").is_file()


def test_input_read_errors_are_not_reported_as_save_failures(tmp_path, monkeypatch):
    """図を組み立てる間の OSError（`.EIC.aef` が読めない等）は FIGURE_SAVE_FAILED にしない（書くときだけ）。"""
    ds = _analysed()

    def locked(_ds, _request):
        raise PermissionError(13, "The process cannot access the file", "x.EIC.aef")

    before = ds.figures.get("pca")
    monkeypatch.setitem(plot_tools._KINDS, "pca", locked)
    for kwargs in ({}, {"output_path": str(tmp_path / "never.png")}):
        with pytest.raises(PermissionError):
            plot_tools.dataset_plot(kind="pca", **kwargs)
    assert not (tmp_path / "never.png").exists() and ds.figures.get("pca") is before


def test_dataset_plot_png_is_capped_in_pixels(tmp_path):
    import io

    from PIL import Image as PILImage

    from metabolomix.plots import render

    _analysed()
    body = plot_tools.dataset_plot(kind="pca")
    with PILImage.open(io.BytesIO(body[1].data)) as image:
        width, height = image.size
    assert width * height <= render.MAX_IMAGE_PIXELS * 1.05


def test_failed_loadings_call_keeps_the_previous_figure():
    """失敗した呼び出しは図の置き場を変えない（前回の図が残り、save_figure はそれを保存する）。"""
    ds = _analysed()
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    before = ds.figures["loadings"]
    for kwargs in ({"pcs": [1, 9]}, {"value": "nope"}):
        err = json.loads(plot_tools.dataset_plot(kind="loadings", **kwargs))["error"]
        assert err["code"] == "DATASET_BAD_REQUEST", kwargs
        assert ds.figures["loadings"] is before, kwargs
