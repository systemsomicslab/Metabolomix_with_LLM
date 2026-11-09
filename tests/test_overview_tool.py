"""最初の概観のツール `dataset_overview`（metabolomix.tools.overview_tools）。"""
from __future__ import annotations

import json

import pytest

from metabolomix.core import session_state
from metabolomix.core.serialization import json_payload
from metabolomix.tools._holder import resolve_holder
from metabolomix.tools.overview_tools import build_overview_payload, dataset_overview
from metabolomix.tools.sample_tools import dataset_exclude
from tests.test_stats_results import _ds


def _install(state, tmp_path):
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


def _strict(text: str) -> dict:
    return json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict constant {c}"))


def _ov(**kw) -> dict:
    return _strict(dataset_overview(**kw))


def _classed(ds):
    for k, fid in enumerate(ds.feature_ids):
        ds.feature_metadata[fid]["ontology"] = "PC" if k < 3 else ("PE" if k < 6 else "TG")
    ds.source_format = "msdial_folder"
    ds.quant_state = {"schema": "lipidmix.quant-state.v1", "normalized": False, "method": None,
                      "unit": None, "steps": ["select"], "skipped": ["normalize"]}
    return ds


def test_overview_two_groups_shape_and_ids(tmp_path):
    ds = _ds()
    ds.source_format = "msdial_folder"
    _install(ds, tmp_path)
    out = _ov()
    assert out["status"] == "success" and out["overview_id"].startswith("ovw_")
    assert out["provisional"] is True
    for key in ("quality", "global", "classes", "species"):
        assert key in out and "next" in out[key]
    assert out["quality"]["facts"][0].startswith("Q1 判定 ")
    assert out["species"]["facts"][0].startswith("S1 T_vs_C（仮の向き）")
    assert "付加体の選別の前" in out["classes"]["omitted"]
    assert "dataset_select_features" in out["classes"]["next"]
    assert out["result_ids"]["compare"] in ds.results and out["result_ids"]["pca"] in ds.results
    assert ds.results[out["overview_id"]]["kind"] == "overview"
    assert ds.last_differential is None and ds.last_pca is None
    assert "counts" in out and out["counts"]["n_features"] == 10
    assert "#S1" in out["note"]


def test_overview_is_deterministic(tmp_path):
    _install(_ds(), tmp_path)
    a, b = _ov(), _ov()
    for key in ("overview_id", "result_ids", "quality", "global", "classes", "species", "state"):
        assert a[key] == b[key]


def test_overview_single_group_skips_comparisons(tmp_path):
    _install(_ds(groups=("C",)), tmp_path)
    out = _ov()
    assert out["status"] == "success"
    assert "群が 1 つ以下" in out["species"]["omitted"]
    assert "compare" not in out["result_ids"] and "pca" in out["result_ids"]
    assert "dataset_set_sample_metadata" in out["global"]["next"]


def test_overview_complete_confounding_does_not_compare(tmp_path):
    ds = _ds()
    for row in ds.sample_metadata_rows:
        if row.get("group") == "T":
            row["batch"] = "B2"
    _install(ds, tmp_path)
    out = _ov()
    assert out["quality"]["verdict"] == "poor"
    assert any("完全交絡" in f for f in out["quality"]["facts"])
    assert "完全交絡" in out["species"]["omitted"]
    assert not any(r.get("kind") == "compare" for r in ds.results.values())


def test_overview_drift_undecided_uses_off_and_says_so(tmp_path):
    ds = _ds(n_qc=6)
    ds.quant_state = {"schema": "lipidmix.quant-state.v1", "normalized": True,
                      "method": "internal_standard", "unit": "pmol", "steps": ["normalize"], "skipped": []}
    _install(ds, tmp_path)
    text = dataset_overview()
    out = _strict(text)
    assert "needs_input" not in out and out["status"] == "success"
    assert "drift skipped" in out["state"]
    assert "dataset_preprocess(drift_mode=" in out["quality"]["next"]
    assert any("ドリフト補正の既定が未決" in f for f in out["quality"]["facts"])


def test_overview_multi_group_with_classes_and_size_budget(tmp_path):
    ds = _install(_classed(_ds(groups=("C", "T", "U"), per_group=5)), tmp_path)
    out = _ov()
    assert out["provisional"] is False
    assert "table" in out["classes"] and out["classes"]["table"].startswith("class\t")
    assert "pairs" in out["species"]
    assert out["global"]["facts"][-1].startswith("G3 moderated F")
    text = json_payload(build_overview_payload(session_state.session.current_handle(), embedded=True)) \
        if hasattr(session_state.session, "current_handle") else dataset_overview()
    assert len(text) <= 6000, len(text)
    assert ds.last_differential is None


def test_overview_strict_json(tmp_path):
    _install(_ds(), tmp_path)
    _strict(dataset_overview())


def test_overview_record_written_and_cites_preprocess(tmp_path):
    _install(_ds(), tmp_path)
    out = _ov()
    assert out["record"]["written"] is True
    records = [json.loads(p.read_text(encoding="utf-8")) for p in tmp_path.rglob("stats/*.json")]
    ov = [r for r in records if r["step"] == "overview"]
    pre = [r for r in records if r["step"] == "preprocess"]
    assert len(ov) == 1 and len(pre) == 1
    assert ov[0]["parent_ids"] == [pre[0]["record_id"]]
    assert ov[0]["inputs"]["result_id"] == out["overview_id"]


def test_overview_without_report_dir_still_succeeds():
    session_state.session.dataset = _ds()
    out = _ov()
    assert out["status"] == "success" and out["record"]["written"] is False


def test_overview_with_dataset_handle_uses_that_handle(tmp_path):
    first = _install(_ds(), tmp_path)
    first_id = session_state.session.current_handle_id
    second = _install(_ds(seed=9), tmp_path)
    assert session_state.session.current_handle_id != first_id
    out = _ov(dataset_handle=first_id)
    assert out["overview_id"] in first.results and out["overview_id"] not in second.results


def test_overview_after_exclusion_recomputes(tmp_path):
    ds = _install(_ds(), tmp_path)
    before = _ov()["overview_id"]
    assert "error" not in json.loads(dataset_exclude(samples=["T0_r0"], reason="test"))
    after = _ov()
    assert after["overview_id"] != before and after["overview_id"] in ds.results


def test_overview_bad_top_n_and_missing_dataset():
    assert _strict(dataset_overview(top_n=99))["error"]["code"] == "DATASET_BAD_REQUEST"
    assert _strict(dataset_overview())["error"]["code"] == "missing_state"


def test_overview_exploratory_dataset_does_not_compare(tmp_path):
    ds = _ds()
    ds.exploratory_only = True
    _install(ds, tmp_path)
    out = _ov()
    assert "探索専用" in out["species"]["omitted"] and "compare" not in out["result_ids"]


def test_overview_embedded_omits_counts(tmp_path):
    _install(_ds(), tmp_path)
    assert "counts" not in build_overview_payload(resolve_holder(None), embedded=True)


def test_overview_non_folder_source_omits_classes(tmp_path):
    _install(_ds(), tmp_path)                 # source_format が解析フォルダでない（mzTab-M など）
    out = _ov()
    assert "定量の鎖" in out["classes"]["omitted"] and "composition" in out["classes"]["next"]


def test_overview_compare_stays_out_of_default_feature_search(tmp_path):
    """仮の向きの概観の比較は、特徴検索の既定（result_ids 省略）に混ざらない。名指しすれば引ける。"""
    from metabolomix.analysis.feature_search import select_results

    ds = _install(_ds(), tmp_path)
    out = _ov()
    rid = out["result_ids"]["compare"]
    assert ds.results[rid]["origin"] == "overview"
    assert select_results(ds, None)["results"] == []
    assert select_results(ds, [rid])["results"]


def test_overview_does_not_relabel_an_explicit_compare(tmp_path):
    """利用者が同じ比較を先に回していたら、その記録を概観の印で上書きしない（既定の特徴検索から消さない）。"""
    from metabolomix.analysis.feature_search import select_results
    from metabolomix.tools.stats_tools import dataset_compare

    ds = _install(_ds(), tmp_path)
    explicit = json.loads(dataset_compare(reference_group="C", test_group="T"))
    rid = explicit.get("result_id")
    assert rid in ds.results, explicit
    _ov()
    assert ds.results[rid].get("origin") is None
    assert select_results(ds, None)["results"]


def _no_features_left(ds):
    """QC の RSD を大きくして、既定の前処理の足切りで特徴が 1 つも残らないようにする。"""
    import numpy as np

    rng = np.random.default_rng(1)
    for j, name in enumerate(ds.sample_names):
        if name.startswith("QC"):
            ds.feature_matrix[:, j] *= rng.uniform(0.01, 100, ds.feature_matrix.shape[0])
    return ds


def test_overview_keeps_quality_when_preprocessing_leaves_nothing(tmp_path):
    """前処理で特徴が残らなくても、品質の段（前処理前の試料から計算する）は出し、他の段は理由つきで省く。"""
    _install(_no_features_left(_ds()), tmp_path)
    out = _ov()
    assert out["status"] == "success", out
    assert out["quality"]["verdict"] == "poor"
    assert any("RSD" in f for f in out["quality"]["facts"])
    for key in ("global", "classes", "species"):
        assert "PREPROCESS_PRECONDITION" in out[key]["omitted"], out[key]
    assert "qc_rsd" in out["global"]["omitted"]           # 足切りの内訳（details）を落とさない
    assert "pca" not in out["result_ids"] and "compare" not in out["result_ids"]
