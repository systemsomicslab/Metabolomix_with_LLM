"""根拠バンドル v3（試料ごとの値と検出・付加体・極性。spec 2026-11-04）。

dataset 経路は常に v3、`samples` と `intensities` を渡さない ARF 経路は v2 のまま。
"""
from __future__ import annotations

import copy
import csv
import json

import numpy as np
import pytest

from metabolomix.analysis import result_state
from metabolomix.analysis.stats_results import as_two_group
from metabolomix.analysis.evidence_export import export_dataset_evidence, write_evidence_bundle
from metabolomix.core import handles
from metabolomix.core.atomic_io import DomainError
from metabolomix.mztab import loading
from metabolomix.tools import dataset_analysis_tools as T
from metabolomix.tools import handle_tools
from tests.engine_chain import declare_groups, engine_compare, hc_ra, register_two_group
from tests.reanalysis_fixtures import POS_FEATURES, write_mztab

_COMPARISON = dict(group_a="A", group_b="B", n_a=3, n_b=3,
                   log2fc_sign="positive means group_b is higher", method="Welch",
                   correction="BH", q_threshold=.05, log2fc_threshold=1, preprocess={})
_ONE_FEATURE = [dict(feature_id="f0", spot_id="0", log2fc=1, p_value=.01, q_value=.02,
                     mean_a=10, mean_b=20, significant=True, adduct="[M+H]1+")]


def _tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _bundle_inputs(tmp_path, *, detected: bool):
    """合成の mzTab-M から統計エンジンの前処理（正規化なし）と HC 対 RA の比較まで進めた (ds, 比較の記録)。

    群は MS-DIAL の試料の宣言で与える（試料の表の行は作らない。注入順は mzTab-M の assay から読む）。
    """
    ds = loading.load_dataset_state(mztab_path=write_mztab(tmp_path / "in" / "t.mzTab"))
    if detected:
        mask = np.ones((len(ds.feature_ids), len(ds.sample_names)), dtype=bool)
        mask[ds.feature_ids.index("1"), ds.sample_names.index("RA1")] = False
        ds.detected_mask = mask
        ds.feature_matrix[ds.feature_ids.index("2"), ds.sample_names.index("HC1")] = np.nan
    declare_groups(ds, hc_ra)
    return ds, register_two_group(ds, "HC", "RA")


def _n_rows(record) -> int:
    return len(as_two_group(record)["results"])


def _handle_with_result(tmp_path) -> tuple[str, str]:
    """ハンドルで開き、前処理と HC 対 RA の比較まで進めた (handle, result_id)。"""
    path = write_mztab(tmp_path / "in" / "t.mzTab")
    opened = json.loads(handle_tools.load_dataset(mztab_path=str(path), library="none"))
    handle = opened["dataset_handle"]
    declare_groups(handles.registry.get(handle).dataset, hc_ra)
    return handle, engine_compare("HC", "RA", dataset_handle=handle)["result_id"]


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
    assert out["n_intensity_rows"] == _n_rows(result) * 6
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
    assert len(rows) == _n_rows(result) * 6                 # 全特徴 × 比較の試料（QC を含まない）
    cell = {(r["feature_id"], r["sample_id"]): r for r in rows}
    assert cell[("1", "RA1")]["detected"] == "false" and cell[("1", "HC1")]["detected"] == "true"
    assert cell[("1", "HC1")]["value"] == "10"               # normalize none・補完の前
    assert cell[("2", "HC1")]["value"] == ""                 # NaN は補完せず空


def test_export_does_not_change_dataset_state(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)

    def state():
        return (result_state.dataset_fingerprint(ds), set(ds.results), set(ds.stat_matrices),
                ds.last_differential["result_id"], ds.metadata_revision)

    before = state()
    export_dataset_evidence(ds, result, tmp_path / "b")
    assert state() == before


def test_without_detection_state_detected_is_empty_and_reason_recorded(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=False)
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert {r["detected"] for r in _tsv(tmp_path / "b" / "intensities.tsv")} == {""}
    assert manifest["population"]["detection_unavailable_reason"]   # 空でない文字列


def test_mask_of_another_shape_is_not_used(tmp_path):
    """形の合わない検出マスクは使わない（理由は mask_shape_mismatch）。マスクはデータの指紋に入るので、
    登録の後に差し替えると結果ごと古くなる——ここでは検出の表の選び方だけを見る。"""
    from metabolomix.analysis.evidence_export import _detection_source

    ds, _ = _bundle_inputs(tmp_path, detected=True)
    ds.detected_mask = ds.detected_mask[:, :-1]
    assert _detection_source(ds, ["HC1", "RA1"]) == (None, "mask_shape_mismatch")


def test_result_samples_outside_mask_axis_leave_detected_empty(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    result = copy.deepcopy(result)
    unit = next(u for u in result["units"].values() if "HC1" in u["samples"])
    unit["samples"].append("ZZ")                             # 試料軸のずれ
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


def test_sme_provenance_drops_only_the_computed_theoretical_mz(tmp_path):
    """`theoretical_mz_computed` の特徴（組成式と付加体から計算した理論 m/z）だけ、rank 1 の候補の来歴の
    theoretical_mass_to_charge を null にする（極性の統合で片側だけ計算した `mixed` の形）。ほかの特徴と rank 2 の
    候補は mzTab-M の値のまま。データセットの候補は変えない。"""
    features = (
        {**POS_FEATURES[0], "sme": [*POS_FEATURES[0]["sme"],
                                    {"name": "Betaine", "inchikey": "KWIUHFFTVRNATP-UHFFFAOYSA-N",
                                     "formula": "C5H11NO2", "adduct": "[M+H]1+", "theo": 118.0863}]},
        {**POS_FEATURES[1], "sme": [{"name": "Phosphocholine", "inchikey": "YHHSONZFOIEMCP-UHFFFAOYSA-O",
                                     "formula": "C5H15NO4P", "adduct": "[M]1+", "theo": 184.0733}]},
        POS_FEATURES[2],
    )
    ds = loading.load_dataset_state(mztab_path=write_mztab(tmp_path / "in" / "t.mzTab",
                                                           features=features))
    declare_groups(ds, hc_ra)
    result = register_two_group(ds, "HC", "RA")
    ds.theoretical_mz_computed = {"1"}
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = [r for r in _tsv(tmp_path / "b" / "annotations.tsv") if r["source"] == "msdial_sme"]
    theo = {(r["feature_id"], r["rank"]): json.loads(r["provenance_json"])["theoretical_mass_to_charge"]
            for r in rows}
    assert theo[("1", "1")] is None
    assert theo[("1", "2")] == pytest.approx(118.0863)
    assert theo[("2", "1")] == pytest.approx(184.0733)
    assert ds.feature_candidates["1"][0]["theoretical_mass_to_charge"] == pytest.approx(104.1070)


def test_stale_result_is_refused(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    ds.metadata_revision += 1                                # 試料の表・除外が登録の後に変わった
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, result, tmp_path / "b")
    assert info.value.code == "STALE_ANALYSIS_RESULT"


def test_refuses_when_reapplied_features_differ_from_result(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    result = copy.deepcopy(result)
    rows = result["engine"]["rows"]
    rows.append({**rows[0], "feature_id": "9"})
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
    handle, result_id = _handle_with_result(tmp_path)        # ハンドルで preprocess → differential
    # 結果が見つからない DomainError が TypeError（exc.details()）にならず封筒で返る
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id + "#nope",
                                                      dataset_handle=handle))
    assert out["error"]["code"] == "CONTRAST_REQUIRED"
    assert not (tmp_path / "b").exists()


def test_exploratory_dataset_is_exported_with_the_mark_in_the_manifest(tmp_path):
    """未完了の実行の探索専用のデータセットも断らず、manifest の exploratory に印と理由を書く。"""
    handle, result_id = _handle_with_result(tmp_path)
    ds = handles.registry.get(handle).dataset
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "a"), result_id, dataset_handle=handle))
    manifest = json.loads((tmp_path / "a" / "manifest.json").read_text(encoding="utf-8"))
    assert "exploratory" not in manifest and "exploratory" not in out   # 探索専用でない束は今までどおり
    ds.exploratory_only = True
    ds.exploratory_reason = "解析ジョブが完了していない（status partial）"
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id, dataset_handle=handle))
    assert out["contract_version"] == 3
    assert out["exploratory"] == {"exploratory_only": True, "reason": ds.exploratory_reason}
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["exploratory"] == out["exploratory"]
    assert manifest["contract_version"] == 3 and manifest["comparison"]["correction"] == "BH"


def test_bundle_result_id_defaults_to_the_latest_comparison(tmp_path):
    handle, result_id = _handle_with_result(tmp_path)
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), dataset_handle=handle))
    assert out["contract_version"] == 3
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["source_analysis_id"] == result_id


def test_bundle_without_a_comparison_is_missing_state(tmp_path):
    from metabolomix.core import session_state
    from metabolomix.mztab.dataset_state import DatasetState

    ds = DatasetState()
    ds.analysis_dir = str(tmp_path)
    session_state.session.dataset = ds
    err = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b")))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_compare"]
    assert not (tmp_path / "b").exists()


def test_bundle_correction_follows_the_records_p_adjust_method(tmp_path):
    handle, result_id = _handle_with_result(tmp_path)
    record = handles.registry.get(handle).dataset.results[result_id]
    record["p_adjust_method"] = "bonferroni"                 # 統計エンジンの記録が持つ補正の方法
    T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id, dataset_handle=handle)
    comparison = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))["comparison"]
    assert comparison["correction"] == "Bonferroni" and "p_adjust_method" not in comparison["statistics"]


def test_tool_passes_handle_polarity(tmp_path):
    handle, result_id = _handle_with_result(tmp_path)
    out = json.loads(T.dataset_export_evidence_bundle(str(tmp_path / "b"), result_id,
                                                      dataset_handle=handle))
    assert out["contract_version"] == 3
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["comparison"]["polarity"] == handles.registry.get(handle).polarity == "positive"
