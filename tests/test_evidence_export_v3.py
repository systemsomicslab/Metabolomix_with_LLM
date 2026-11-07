"""根拠バンドル v3（試料ごとの値と検出・付加体・極性。spec 2026-11-04）。

dataset 経路は常に v3、`samples` と `intensities` を渡さない ARF 経路は v2 のまま。
"""
from __future__ import annotations

import csv
import json

import numpy as np
import pytest

from metabolomix.analysis import dataset_service, feature_search, result_state
from metabolomix.analysis.evidence_export import export_dataset_evidence, write_evidence_bundle
from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools
from tests.reanalysis_fixtures import write_mztab

HC = ["HC1", "HC2", "HC3"]
RA = ["RA1", "RA2", "RA3"]
_COMPARISON = dict(group_a="A", group_b="B", n_a=3, n_b=3,
                   log2fc_sign="positive means group_b is higher", method="Welch",
                   correction="BH", q_threshold=.05, log2fc_threshold=1, preprocess={})
_ONE_FEATURE = [dict(feature_id="f0", spot_id="0", log2fc=1, p_value=.01, q_value=.02,
                     mean_a=10, mean_b=20, significant=True, adduct="[M+H]1+")]


def _tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _bundle_inputs(tmp_path, *, detected: bool):
    """合成の mzTab-M から前処理（正規化なし）と HC 対 RA の比較まで進めた (ds, result)。"""
    ds = loading.load_dataset_state(mztab_path=write_mztab(tmp_path / "in" / "t.mzTab"))
    if detected:
        mask = np.ones((len(ds.feature_ids), len(ds.sample_names)), dtype=bool)
        mask[ds.feature_ids.index("1"), ds.sample_names.index("RA1")] = False
        ds.detected_mask = mask
        ds.feature_matrix[ds.feature_ids.index("2"), ds.sample_names.index("HC1")] = np.nan
    dataset_service.preprocess_dataset(ds, {"normalize": "none"})
    result = dataset_service.compare_dataset(ds, HC, RA, group_a_label="HC", group_b_label="RA")
    return ds, result


def _handle_with_result(tmp_path) -> tuple[str, str]:
    """ハンドルで開き、前処理と HC 対 RA の比較まで進めた (handle, result_id)。"""
    path = write_mztab(tmp_path / "in" / "t.mzTab")
    opened = json.loads(handle_tools.load_dataset(mztab_path=str(path), library="none"))
    handle = opened["dataset_handle"]
    T.dataset_preprocess(dataset_handle=handle)
    result = json.loads(T.dataset_differential(group_a=HC, group_b=RA, dataset_handle=handle))
    return handle, result["result_id"]


def test_v3_manifest_and_roles(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    out = export_dataset_evidence(ds, result, tmp_path / "b", polarity="positive")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert out["contract_version"] == 3 and manifest["contract_version"] == 3
    assert manifest["exporter"] == {"name": "ms-data-parser", "version": "evidence-v3"}
    assert manifest["comparison"]["polarity"] == "positive"
    assert manifest["intensities_stage"] == "pre_imputation"
    assert {f["role"] for f in manifest["files"]} >= {"features", "annotations", "spectrum_links",
                                                      "samples", "intensities"}
    assert manifest["population"]["detection_unavailable_reason"] is None
    assert out["n_samples"] == len(ds.sample_names)
    assert out["n_intensity_rows"] == len(result["results"]) * 6
    assert out["intensities_stage"] == "pre_imputation"


def test_samples_group_is_empty_for_qc(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    # 名前からの推定は `_` 区切りのトークンしか見ず、`QC1` は qc にならない。役割は実験情報シートから取る。
    ds.sample_metadata_rows = [{"sample_id": name, "role": "qc" if name.startswith("QC") else "sample"}
                               for name in ds.sample_names]
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = {r["sample_id"]: r for r in _tsv(tmp_path / "b" / "samples.tsv")}
    assert [rows[s]["group"] for s in ("HC1", "RA1", "QC1")] == ["a", "b", ""]
    assert rows["QC1"]["role"] == "qc" and rows["HC1"]["role"] == "sample"
    assert rows["HC1"]["injection_order"] in ("", "1")      # 分からなければ空、分かれば整数


def test_injection_order_is_integer_or_empty(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    aid = ds.sample_assay_ids[ds.sample_names.index("HC1")]
    other = ds.sample_assay_ids[ds.sample_names.index("HC2")]
    ds.assay_metadata = {aid: {"run_order": 7}, other: {"run_order": "late"}}
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = {r["sample_id"]: r for r in _tsv(tmp_path / "b" / "samples.tsv")}
    assert rows["HC1"]["injection_order"] == "7" and rows["HC2"]["injection_order"] == ""


def test_intensities_long_format_detected_and_pre_imputation_value(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = _tsv(tmp_path / "b" / "intensities.tsv")
    assert len(rows) == len(result["results"]) * 6          # 全特徴 × 比較の試料（QC を含まない）
    cell = {(r["feature_id"], r["sample_id"]): r for r in rows}
    assert cell[("1", "RA1")]["detected"] == "false" and cell[("1", "HC1")]["detected"] == "true"
    assert cell[("1", "HC1")]["value"] == "10"               # normalize none・補完の前
    assert cell[("2", "HC1")]["value"] == ""                 # NaN は補完せず空


def test_reapplying_recipe_does_not_change_dataset_state(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)

    def state():
        return (result_state.array_fingerprint(ds.pp_matrix), ds.preprocess_id,
                result_state.dataset_fingerprint(ds), dict(ds.preprocessing_recipe),
                set(ds.results), list(ds.pp_sample_names), list(ds.pp_feature_names))

    before = state()
    export_dataset_evidence(ds, result, tmp_path / "b")
    assert state() == before


def test_without_detection_state_detected_is_empty_and_reason_recorded(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=False)
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert {r["detected"] for r in _tsv(tmp_path / "b" / "intensities.tsv")} == {""}
    assert manifest["population"]["detection_unavailable_reason"]   # 空でない文字列


def test_mask_of_another_shape_leaves_detected_empty(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    ds.detected_mask = ds.detected_mask[:, :-1]
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["population"]["detection_unavailable_reason"] == "mask_shape_mismatch"
    assert {r["detected"] for r in _tsv(tmp_path / "b" / "intensities.tsv")} == {""}


def test_result_samples_outside_mask_axis_leave_detected_empty(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    result = {**result, "samples_a": ["HC1", "HC2", "ZZ"]}  # 試料軸のずれ
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["population"]["detection_unavailable_reason"] == "sample_axis_mismatch"
    rows = _tsv(tmp_path / "b" / "intensities.tsv")
    assert {r["detected"] for r in rows} == {""}             # false で埋めない
    assert {r["value"] for r in rows if r["sample_id"] == "ZZ"} == {""}


def test_features_have_adduct_column(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    export_dataset_evidence(ds, result, tmp_path / "b")
    features = {r["feature_id"]: r for r in _tsv(tmp_path / "b" / "features.tsv")}
    assert features["1"]["adduct"] == "[M+H]1+"              # SML の adduct_ions


def test_post_imputation_when_recipe_unknown(tmp_path, monkeypatch):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    monkeypatch.setattr(feature_search, "_result_recipe", lambda ds, result: None)
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["intensities_stage"] == "post_imputation"
    cell = {(r["feature_id"], r["sample_id"]): r for r in _tsv(tmp_path / "b" / "intensities.tsv")}
    assert cell[("2", "HC1")]["value"] != ""                 # 親の行列は補完の後


def test_refuses_when_neither_recipe_nor_parent_matrix(tmp_path, monkeypatch):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    monkeypatch.setattr(feature_search, "_result_recipe", lambda ds, result: None)
    ds.preprocess_id = "other"
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, result, tmp_path / "b")
    assert info.value.code == "EVIDENCE_INTENSITIES_UNAVAILABLE"
    assert not (tmp_path / "b").exists()


def test_refuses_handoff_result_without_recipe_on_another_preprocess(tmp_path):
    """`assert_current` を通る形（引き継ぎの結果は親の前処理を照合しない）でも止まる。"""
    ds, result = _bundle_inputs(tmp_path, detected=True)
    provenance = {**result["provenance"], "handoff": {"preprocess": None}, "parent_ids": ["run"],
                  "dataset_hash": result_state.dataset_fingerprint(ds)}
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, {**result, "provenance": provenance}, tmp_path / "b")
    assert info.value.code == "EVIDENCE_INTENSITIES_UNAVAILABLE"
    assert not (tmp_path / "b").exists()


def test_handoff_result_writes_recipe_to_comparison(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    provenance = {**result["provenance"], "handoff": {"preprocess": {"normalize": "none"}},
                  "parent_ids": ["run"], "dataset_hash": result_state.dataset_fingerprint(ds)}
    export_dataset_evidence(ds, {**result, "provenance": provenance}, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["intensities_stage"] == "pre_imputation"
    assert manifest["comparison"]["preprocess"] == {"normalize": "none"}


def test_stale_result_with_known_recipe_stays_stale(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    ds.preprocess_id = "other"
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, result, tmp_path / "b")
    assert info.value.code == "STALE_ANALYSIS_RESULT"


def test_refuses_when_reapplied_features_differ_from_result(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    result = {**result, "results": [*result["results"], {**result["results"][0], "feature": "9"}]}
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, result, tmp_path / "b")
    assert info.value.code == "EVIDENCE_INTENSITIES_MISMATCH"
    assert not (tmp_path / "b").exists()


def test_writer_refuses_samples_without_intensities(tmp_path):
    with pytest.raises(ValueError):
        write_evidence_bundle(comparison=_COMPARISON, features=_ONE_FEATURE, annotations=[],
                              spectrum_links=[], spectra_files=[], output_dir=tmp_path / "o",
                              samples=[])
    assert not (tmp_path / "o").exists()


def test_writer_without_samples_stays_v2(tmp_path):     # ARF 経路
    out = write_evidence_bundle(comparison=_COMPARISON, features=_ONE_FEATURE, annotations=[],
                                spectrum_links=[], spectra_files=[], output_dir=tmp_path / "o")
    assert out["contract_version"] == 2
    assert "adduct" not in _tsv(tmp_path / "o" / "features.tsv")[0]


def test_tool_description_names_v3_and_domain_error_is_envelope(tmp_path):
    assert "v3" in T.dataset_export_evidence_bundle.__doc__
    # 探索専用のデータセットの DomainError が TypeError（exc.details()）にならず封筒で返る
    handle, result_id = _handle_with_result(tmp_path)        # ハンドルで preprocess → differential
    handles.registry.get(handle).dataset.exploratory_only = True
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id,
                                                      dataset_handle=handle))
    assert out["error"]["code"] == "EXPLORATORY_ONLY_DATASET"


def test_tool_passes_handle_polarity(tmp_path):
    handle, result_id = _handle_with_result(tmp_path)
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id,
                                                      dataset_handle=handle))
    assert out["contract_version"] == 3
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["comparison"]["polarity"] == handles.registry.get(handle).polarity == "positive"
