"""dataset_plot（volcano・pca・loadings・eic・eic_compounds）と図の置き場（B-2b Task 6）。

純粋な組み立て（`plots/*`）は各々のテストが固定している。ここでは `DatasetState` を入力に、
ツール層の仕事——結果の選び方（古い結果は描かない）・引数の取り違えの拒否・`.EIC.aef` の出所・
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
from metabolomix.tools.dataset_analysis_tools import (
    dataset_differential, dataset_pca)
from tests.config_fixtures import use_config
from tests.curation_fixtures import write_alignment_set
from tests.msdial_folder_fixtures import analysis_file_row, write_mddata
from tests.v1_chain import v1_preprocess

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


def _analysed(*, with_differential=True) -> DatasetState:
    ds = _synthetic()
    _set_metadata(ds, group=lambda n: n.split("_")[0],
                  batch=lambda n: "b1" if n.endswith(("0", "1")) else "b2")
    v1_preprocess()
    assert json.loads(dataset_pca(n_components=3))["status"] == "success"
    if with_differential:
        a = [n for n in ds.pp_sample_names if n.startswith("ctrl")]
        b = [n for n in ds.pp_sample_names if n.startswith("treat")]
        assert json.loads(dataset_differential(group_a=a, group_b=b))["status"] == "success"
    return ds


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

def test_volcano_from_dataset_differential():
    _analysed()
    out = plot_tools.dataset_plot(kind="volcano")
    assert isinstance(out, list) and len(out) == 2
    caption, image = out
    assert isinstance(caption, str) and "group_a" in caption and "n=4" in caption
    assert isinstance(image, Image) and image.data[:4] == PNG_MAGIC
    payload = _payload(kind="volcano")
    assert payload["plot_schema"] == "lipidmix.volcano.v1"
    assert payload["comparison"]["group_a"] == "group_a"


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


def test_irrelevant_argument_rejected():
    _analysed()
    err = json.loads(plot_tools.dataset_plot(kind="volcano", top_n=5))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    assert "top_n" in err["message"] and "loadings" in err["message"]
    err = json.loads(plot_tools.dataset_plot(kind="pca", feature_id="1"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST" and "eic" in err["message"]
    err = json.loads(plot_tools.dataset_plot(kind="loadings", color_by="batch"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST" and "pca" in err["message"]
    err = json.loads(plot_tools.dataset_plot(kind="eic", result_id="res_x", feature_id="0"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST" and "result_id" in err["message"]
    err = json.loads(plot_tools.dataset_plot(kind="heatmap"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    assert all(kind in err["message"] for kind in plot_tools._KINDS)


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
    ds = _analysed()
    pca_id = ds.last_pca["provenance"]["result_id"]
    diff_id = ds.last_differential["provenance"]["result_id"]
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
                       ("volcano", "dataset_differential")):
        err = json.loads(plot_tools.dataset_plot(kind=kind))["error"]
        assert err["code"] == "missing_state" and tool in err["required_tools"], kind
    # 古い図は save_figure からも保存できない
    err = json.loads(reports.save_figure("pca", "stale-1"))["error"]
    assert err["code"] == "missing_state" and "dataset_plot" in err["required_tools"]


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
    assert err["code"] == "missing_state" and "dataset_plot" in err["required_tools"]
    assert not (tmp_path / "reports" / "figures" / "meta-1_pca.png").exists()
    v1_preprocess()
    dataset_pca(n_components=3)
    payload = _payload(kind="pca")                    # 解析と図をやり直せば新しい群で保存できる
    assert {p["group"] for p in payload["points"]} == {"NEW"}
    assert "meta-1_pca.png" in reports.save_figure("pca", "meta-1")


def test_pca_figure_with_baked_metadata_is_stale_when_the_revision_moves(tmp_path):
    """結果が現在のままでも（群は前処理の入力ではない）、色に焼き込んだ試料の情報が更新されたら古い。"""
    from metabolomix.analysis.result_state import is_current

    ds = _analysed()
    plot_tools.dataset_plot(kind="pca")
    plot_tools.dataset_plot(kind="loadings", top_n=3)
    for row in ds.sample_metadata_rows:
        row["group"] = "NEW"
    ds.metadata_revision += 1                         # 情報の更新（apply_metadata が進める版）
    assert is_current(ds, ds.last_pca) and "pca" in ds.figures
    err = json.loads(reports.save_figure("pca", "rev-1"))["error"]
    assert err["code"] == "missing_state" and "dataset_plot" in err["required_tools"]
    assert "作り直" in err["message"]
    assert not (tmp_path / "reports" / "figures" / "rev-1_pca.png").exists()
    # 試料の情報を焼き込まない図（loadings）は影響を受けない
    assert "rev-1_loadings.png" in reports.save_figure("loadings", "rev-1")


def test_figures_are_not_part_of_the_fingerprint():
    from metabolomix.analysis.result_state import dataset_fingerprint

    ds = _analysed()
    before = dataset_fingerprint(ds)
    _payload(kind="pca")
    assert ds.figures and dataset_fingerprint(ds) == before


def test_stale_pca_result_after_reprocessing_is_refused():
    ds = _analysed()
    old = ds.last_pca["provenance"]["result_id"]
    v1_preprocess(normalize="median")           # 前処理をやり直す（PCA は捨てられる）
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
    ds = _analysed()
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
    pca_id = ds.last_pca["provenance"]["result_id"]
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
    assert list(params) == ["kind", "analysis_id", "title", "result_id", "dataset_handle"]
    out = json.loads(reports.save_figure("heatmap", "x"))
    assert out["status"] == "error" and "pca" in out["message"]
    assert "pca_loadings" not in out["message"]
    err = json.loads(reports.save_figure("eic", "x", result_id="r1"))
    assert err["status"] == "error" and "result_id" in err["message"]


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
