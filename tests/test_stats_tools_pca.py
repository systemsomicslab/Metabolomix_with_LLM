"""`dataset_pca`（統計エンジンの PCA ツール）と PCA の図。fixture はテスト自身が tmp に作る。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from metabolomix.analysis import stats_results
from metabolomix.analysis.pca import loading_correlations, run_pca_fit_subset
from metabolomix.core import session_state
from metabolomix.core.tool_helpers import _pca_scatter_arrays, dataset_pca_plot
from metabolomix.plots import pca_loadings
from metabolomix.tools import plot_tools, stats_tools
from metabolomix.tools.stats_tools import dataset_pca, dataset_preprocess
from tests.test_stats_results import _ds

FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")
PAYLOAD_KEYS = {"status", "result_id", "matrix_id", "scaling", "log2", "imputation", "n_fit",
                "n_qc_projected", "explained_variance_ratio", "groups", "record", "loadings_note"}


def _install(state, tmp_path):
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


@pytest.fixture
def ds(tmp_path):
    return _install(_ds(), tmp_path)


def _pre(**kw) -> dict:
    return json.loads(dataset_preprocess(**{**FLAT_OFF, **kw}))


def _pca(**kw) -> dict:
    return json.loads(dataset_pca(**kw))


def _strict(text: str) -> dict:
    return json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))


def test_pca_on_matrix_qc_projected_and_groups_counted(ds):
    pre = _pre()
    out = _strict(dataset_pca(matrix_id=pre["matrix_id"], n_components=3))
    assert set(out) == PAYLOAD_KEYS
    assert out["status"] == "success" and out["result_id"].startswith("pca_")
    assert out["matrix_id"] == pre["matrix_id"]
    assert out["scaling"] == "autoscale" and out["log2"] is True
    assert out["n_fit"] == 8 and out["n_qc_projected"] == 4          # 4 個体 × 2 群、QC 4 本
    assert out["groups"] == {"C": 4, "T": 4}
    assert len(out["explained_variance_ratio"]) == 3
    assert out["imputation"]["used"] is False
    assert out["record"]["written"] is True
    record_file = json.loads(Path(out["record"]["path"]).read_text(encoding="utf-8"))
    assert record_file["step"] == "pca"
    assert record_file["parent_ids"] == [pre["record"]["record_id"]]
    assert record_file["inputs"]["result_id"] == out["result_id"]
    assert record_file["settings"]["scaling"] == "autoscale"
    assert record_file["summary"]["groups"] == {"C": 4, "T": 4}
    # 全量は戻り値でなくセッションに
    assert "loadings" not in out and "scores" not in out
    kept = ds.results[out["result_id"]]
    assert kept is ds.last_pca and kept["kind"] == "stat_pca"
    engine = kept["engine"]
    assert np.asarray(engine["scores"]["units"]).shape == (8, 3)
    assert np.asarray(engine["scores"]["qc"]).shape == (4, 3)
    assert len(engine["loadings"]) == 3
    # QC は軸に効かない: QC の値を変えても個体のスコアは同じ
    entry = stats_results.get_matrix(ds, pre["matrix_id"])
    assert engine["unit_ids"] == [u.unit_id for u in entry["prepared"].units]


def test_pca_without_matrix_id_registers_the_default_matrix(ds):
    out = _pca()
    assert out["status"] == "success"
    entry = stats_results.get_matrix(ds, out["matrix_id"])
    assert entry is not None and entry["record_id"]
    assert out["record"]["preprocess_record"]["written"] is True
    assert _pca()["matrix_id"] == out["matrix_id"]                    # 同じ行列を再利用
    # 数字は同じ既定の行列を dataset_preprocess で作って渡したものと同じ
    mx = json.loads(dataset_preprocess())["matrix_id"]
    assert mx == out["matrix_id"]


def test_pca_group_levels_only_affect_labels_not_numbers(ds):
    mx = _pre()["matrix_id"]
    plain = _pca(matrix_id=mx)
    labelled = _pca(matrix_id=mx, group_levels=["C0", "C1"])          # matrix_id との併用は表示だけ
    assert labelled["status"] == "success"
    assert labelled["groups"] == {"C0": 1, "C1": 1, "other": 6}
    assert plain["groups"] == {"C": 4, "T": 4}
    assert labelled["explained_variance_ratio"] == plain["explained_variance_ratio"]
    assert labelled["result_id"] != plain["result_id"]
    a, b = ds.results[plain["result_id"]]["engine"], ds.results[labelled["result_id"]]["engine"]
    assert a["scores"] == b["scores"] and a["loadings"] == b["loadings"]
    # 行列は群の指定で選ばれない: 群を付けた省略呼び出しも、群なしの省略呼び出しと同じ既定の行列
    plain_default = _pca()
    implicit = _pca(group_factors=[["C", "T"], ["r0", "r1"]])
    assert implicit["status"] == "success" and implicit["matrix_id"] == plain_default["matrix_id"]
    assert implicit["groups"] != plain_default["groups"]
    assert (ds.results[implicit["result_id"]]["engine"]["scores"]
            == ds.results[plain_default["result_id"]]["engine"]["scores"])


def test_pca_missing_matrix_is_missing_state(ds):
    err = _pca(matrix_id="mx_nope")["error"]
    assert err["code"] == "missing_state" and err["state"] == "stat_matrix"
    assert err["required_tools"] == ["dataset_preprocess"]
    session_state.session.dataset = None
    assert _pca()["error"]["code"] == "missing_state"


def test_pca_knn_failure_envelope(ds):
    ds.feature_matrix = np.array(ds.feature_matrix, dtype=float)
    cols = [i for i, n in enumerate(ds.sample_names) if n.startswith("C0_")]
    ds.feature_matrix[:8, cols] = 0.0               # 1 個体（2 注入とも）で特徴の 8/10 が欠測（knn の上限を超える）
    out = _pca(impute="knn")
    err = out["error"]
    assert err["code"] == "PCA_IMPUTATION_FAILED"
    assert err["details"]["matrix_id"].startswith("mx_")             # 行列は登録済み。呼び直せる
    assert err["details"]["cause"]["code"] == "RECIPE_PRECONDITION"
    again = _pca(matrix_id=err["details"]["matrix_id"], impute="half_min")
    assert again["status"] == "success" and again["imputation"]["used"] is True


def test_pca_specification_errors_are_envelopes(ds):
    mx = _pre()["matrix_id"]
    assert _pca(matrix_id=mx, scaling="bogus")["error"]["code"] == "PCA_SPECIFICATION_INVALID"
    assert _pca(matrix_id=mx, n_components=0)["error"]["code"] == "PCA_SPECIFICATION_INVALID"


# --- 図 ----------------------------------------------------------------------------------

def _payload(**kwargs) -> dict:
    out = plot_tools.dataset_plot(output="payload", **kwargs)
    assert isinstance(out, str), out
    return json.loads(out)


def test_dataset_plot_pca_and_loadings_from_stat_pca(ds):
    mx = _pre()["matrix_id"]
    out = _pca(matrix_id=mx, n_components=3)
    rid = out["result_id"]
    pay = _payload(kind="pca", result_id=rid)
    assert pay["plot_schema"] == "lipidmix.pca_scores.v1" and pay["result_id"] == rid
    assert len(pay["points"]) == 12                                   # 個体 8 + QC 4
    qc = [p for p in pay["points"] if p["group"] == "QC"]
    assert len(qc) == 4 and all(p["fitted"] is False for p in qc)
    assert {p["group"] for p in pay["points"] if "fitted" not in p} == {"C", "T"}
    assert any("投影" in c for c in pay["caveats"])
    assert {p["group"] for p in _payload(kind="pca", color_by="role")["points"]} == {"sample", "qc"}
    assert {p["group"] for p in _payload(kind="pca", color_by="batch")["points"]} == {"B1", "(unset)"}
    # 直近（result_id 省略）でも描ける。画像も返る
    assert plot_tools.dataset_plot(kind="pca") is not None

    load = _payload(kind="loadings", pcs=[1, 2], top_n=3)
    assert load["plot_schema"] == "lipidmix.pca_loadings.v1" and load["value"] == "r"
    assert [p["pc"] for p in load["panels"]] == ["PC1", "PC2"]
    rows = load["panels"][0]["rows"]
    assert rows and all(r["label"].startswith("F") for r in rows)     # 特徴名は feature_metadata から
    assert all(-1.0 <= r["r"] <= 1.0 for r in rows)
    # 古くなれば描かない
    from metabolomix.tools import sample_tools
    sample_tools.dataset_exclude(samples=["T3_r0"], reason="outlier")
    err = json.loads(plot_tools.dataset_plot(kind="pca", result_id=rid))["error"]
    assert err["code"] == "STALE_ANALYSIS_RESULT"


def test_loadings_r_matches_the_direct_correlation(ds):
    mx = _pre()["matrix_id"]
    rid = _pca(matrix_id=mx, n_components=3)["result_id"]
    record = ds.results[rid]
    feats, r_ok = pca_loadings.dataset_features(ds, record)
    assert r_ok is True
    prepared = stats_results.get_matrix(ds, mx)["prepared"]
    X = np.asarray(prepared.log2, dtype=float)
    ref = run_pca_fit_subset(X, np.ones(len(X), bool), 3, scaling="autoscale", svd_solver="full")
    r = np.asarray(loading_correlations(ref["loadings"], ref["singular_values"], ref["n_fit"]))
    got = np.array([f["r"] for f in feats]).T                         # 成分 × 特徴
    np.testing.assert_allclose(got, r, atol=1e-9)
    pareto = _pca(matrix_id=mx, scaling="pareto", n_components=3)["result_id"]
    assert pca_loadings.dataset_features(ds, ds.results[pareto])[1] is False


def test_full_svd_is_deterministic_over_500_features():
    """500 超の特徴で `svd_solver="full"` は実行ごとに同じ軸を返す（"auto" は randomized を選びうる）。"""
    x = np.random.default_rng(3).normal(size=(40, 600))
    a = run_pca_fit_subset(x, np.ones(40, bool), 5, svd_solver="full")
    b = run_pca_fit_subset(x, np.ones(40, bool), 5, svd_solver="full")
    assert a["components"] == b["components"] and a["loadings"] == b["loadings"]


def test_pipeline_v1_pca_plot_payload_still_builds():
    """pipeline が残す v1 の PCA 結果（scores 行）も散布図に射影できる（D 段まで）。"""
    v1 = {"scores": [{"name": "S1", "role": "sample", "PC1": 1.0, "PC2": 2.0},
                     {"name": "S2", "role": "sample", "PC1": -1.0, "PC2": 0.5},
                     {"name": "S3", "role": "qc", "PC1": 0.0, "PC2": 0.0}],
          "explained_variance_ratio": [0.6, 0.3]}
    xs, ys, labels, x_label, y_label, _t = _pca_scatter_arrays(dataset_pca_plot(v1))
    assert xs == [1.0, -1.0, 0.0] and ys == [2.0, 0.5, 0.0] and labels == ["S1", "S2", "S3"]
    assert x_label == "PC1 (60.00%)" and y_label == "PC2 (30.00%)"
    assert "group" not in dataset_pca_plot(v1)["points"][0]
    assert dataset_pca_plot({"scores": [{"name": "a", "PC1": 1.0}]}) is None


def test_stat_pca_record_projects_to_the_scatter_dict(ds):
    mx = _pre()["matrix_id"]
    record = ds.results[_pca(matrix_id=mx)["result_id"]]
    plot = dataset_pca_plot(record)
    assert len(plot["points"]) == 12
    assert plot["points"][-1]["group"] == "QC" and plot["points"][-1]["fitted"] is False
    assert plot["x_label"].startswith("PC1 (")


def test_dataset_pca_is_registered_once_and_append_typed():
    from metabolomix.core.mcp_core import mcp
    tools = [t for t in mcp._tool_manager.list_tools() if t.name == "dataset_pca"]
    assert len(tools) == 1
    assert "dataset_pca" in stats_tools.__all__
    assert tools[0].annotations.readOnlyHint is False
    assert tools[0].annotations.idempotentHint is False
    assert "log_transform" not in tools[0].parameters["properties"]
