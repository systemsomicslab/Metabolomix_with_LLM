"""統計エンジンの `compare` の記録を読む利用者（書き出し・根拠バンドル・volcano・特徴検索・新しい知見の候補）。

どの利用者も `stats_results.as_two_group` を通して読む。fixture はテスト自身が作る
（`tests/test_stats_results.py` の `_ds`。除外・群の訂正が `metadata_revision` を進める状況も同じ DatasetState で作る）。
"""
from __future__ import annotations

import csv
import json

import numpy as np
import pytest

from metabolomix.analysis import evidence_export, export_contract, feature_search, stats_results
from metabolomix.analysis.dataset_export import export_dataset_result
from metabolomix.analysis.dataset_exclusions import apply_exclusions
from metabolomix.core import session_state
from metabolomix.core.atomic_io import DomainError
from metabolomix.tools import dataset_analysis_tools as DAT
from metabolomix.tools import plot_tools, reanalysis_tools
from metabolomix.tools.stats_tools import dataset_compare, dataset_preprocess
from tests.test_stats_results import _ds
from tests.v1_chain import v1_preprocess

FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")

#: v1 の書き出しがもともと持つメタ行のキー（これ以外を足さない）
V1_META_KEYS = {"contract_version", "exported_at", "source_mztab", "source_job", "id_space", "result_id",
                "preprocess_id", "source_verification", "group_a", "n_a", "group_b", "n_b", "log2fc_sign",
                "q_threshold", "log2fc_threshold", "log_transform", "preprocess", "n_features_total",
                "n_with_inchikey", "n_unannotated"}
NEW_META_KEYS = {"test", "effect_size_definition", "ci_level", "verdict_q_threshold", "equivalence_fc",
                 "bh_population", "bh_n", "q_interpretation", "unit", "n_units_a", "n_units_b",
                 "mean_definition", "matrix_id", "rule_hash"}


def _annotate(ds):
    """全特徴に InChIKey を付ける（書き出しは InChIKey の付いた行だけを出す）。"""
    for k, fid in enumerate(ds.feature_ids):
        ds.feature_metadata[fid].update(inchikey=f"AAAAAAAAAAAAAA-BBBBBBBBBB-{chr(65 + k)}",
                                        mz=100.0 + k, rt=1.0 + k / 10)
    return ds


def _install(state, tmp_path):
    state.analysis_dir = str(tmp_path / "analysis")
    session_state.session.dataset = state
    return state


@pytest.fixture
def ds(tmp_path):
    return _install(_annotate(_ds()), tmp_path)


def _compare(**kw) -> dict:
    out = json.loads(dataset_preprocess(**FLAT_OFF))
    return json.loads(dataset_compare(matrix_id=out["matrix_id"], **kw))


def _read(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    meta = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    return meta, body


def _parse_meta(meta_lines) -> tuple[dict, list[str]]:
    """massbank-context の `parse_differential` と同じ読み方（タブ区切りの `# key = value` を全部拾う）。"""
    out, keys = {}, []
    for line in meta_lines:
        for chunk in line.lstrip("#").split("\t"):
            if "=" in chunk:
                key, _, value = chunk.partition("=")
                out[key.strip()] = value.strip()
                keys.append(key.strip())
    return out, keys


def _tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


# --- 書き出し -------------------------------------------------------------------------

def test_export_compare_result_keeps_15_columns_and_appends_meta_lines(ds, tmp_path):
    out = _compare(reference_group="C", test_group="T")
    path = tmp_path / "out.tsv"
    info = export_dataset_result(ds, ds.last_differential, path)
    meta, body = _read(path)
    assert body[0].split("\t") == export_contract.EXPORT_COLUMNS
    assert len(body) - 1 == info["n_with_inchikey"] == 10
    assert info["contrast_id"] == "T_vs_C" and info["result_id"] == out["result_id"]
    # 新しい行は msi_note の後ろだけ。既存の行の並びは動かない
    msi = next(i for i, line in enumerate(meta) if line.startswith("# ontology / msi_level"))
    assert msi == len(meta) - 1 - 5
    values, _ = _parse_meta(meta[msi + 1:])
    assert values["test"] == "moderated_t" and values["effect_size_definition"] == "log2_mean_difference"
    assert values["ci_level"] == "0.95" and values["verdict_q_threshold"] == "0.05"
    assert values["bh_population"] == "matrix" and values["bh_n"] == "10"
    assert values["q_interpretation"] == "discovery_fdr"
    assert values["unit"] == "subject" and values["n_units_a"] == "4" and values["n_units_b"] == "4"
    assert values["mean_definition"] == "linear arithmetic mean over units after preprocessing"
    assert values["matrix_id"] == out["matrix_id"] and values["rule_hash"] == out["rule_hash"]
    # n_a / n_b は注入の数（4 個体 × 技術反復 2）、単位の数は別に書く
    head, _ = _parse_meta(meta)
    assert head["n_a"] == "8" and head["n_b"] == "8" and head["contract_version"] == "1"
    assert head["log_transform"] == "true" and head["preprocess_id"] == out["matrix_id"]
    # significant は契約の規則（q ≤ しきい値 かつ |log2fc| ≥ log2fc_threshold）
    two = stats_results.as_two_group(ds.last_differential)
    expected = {row["feature"]: export_contract.is_significant(
        q=row["q"], log2fc=row["log2fc"], q_threshold=two["q_threshold"],
        log2fc_threshold=two["log2fc_threshold"]) for row in two["results"]}
    cells = [line.split("\t") for line in body[1:]]
    assert {c[0]: c[14] == "true" for c in cells} == {k: v for k, v in expected.items()}
    assert any(c[14] == "true" for c in cells)


def test_export_meta_parses_with_massbank_rules(ds, tmp_path):
    _compare(reference_group="C", test_group="T")
    path = tmp_path / "out.tsv"
    export_dataset_result(ds, ds.last_differential, path)
    meta, _ = _read(path)
    values, keys = _parse_meta(meta)
    assert len(keys) == len(set(keys)), "メタのキーが重複している（後勝ちで消える）"
    assert NEW_META_KEYS <= set(values)
    assert set(values) <= V1_META_KEYS | NEW_META_KEYS | {"preprocess_id"}
    assert not (NEW_META_KEYS & V1_META_KEYS)
    # 新しい行は「# key = value」の形で、値にタブ・改行を含まない
    assert all(line.startswith("# ") for line in meta)


def test_export_v1_two_group_output_unchanged(tmp_path):
    from metabolomix.analysis.dataset_service import compare_dataset
    ds = _annotate(_ds())
    _install(ds, tmp_path)
    v1_preprocess()
    result = compare_dataset(ds, [n for n in ds.sample_names if n.startswith("C")],
                             [n for n in ds.sample_names if n.startswith("T")])
    assert result["kind"] == "two_group" and "statistics" not in result
    path = tmp_path / "v1.tsv"
    export_dataset_result(ds, result, path)
    meta, body = _read(path)
    values, keys = _parse_meta(meta)
    assert set(values) == V1_META_KEYS and len(keys) == len(set(keys))      # 新しい行は足さない
    assert meta[-1].startswith("# ontology / msi_level")                    # msi_note が最後
    assert body[0].split("\t") == export_contract.EXPORT_COLUMNS
    assert stats_results.as_two_group(result) is result                     # v1 の結果はそのまま


def test_export_multi_contrast_requires_contrast_id(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    assert [c["contrast_id"] for c in out["contrasts"]] == ["T_vs_C", "U_vs_C"]
    path = tmp_path / "out.tsv"
    err = json.loads(DAT.dataset_export_differential(str(path)))["error"]
    assert err["code"] == "CONTRAST_REQUIRED"
    assert err["details"]["two_group_contrast_ids"] == ["T_vs_C", "U_vs_C"]
    assert not path.exists()
    err = json.loads(DAT.dataset_export_differential(str(path), contrast_id="nope"))["error"]
    assert err["code"] == "CONTRAST_REQUIRED"
    ok = json.loads(DAT.dataset_export_differential(str(path), contrast_id="U_vs_C"))
    assert ok["status"] == "success" and ok["group_a"] == "C" and ok["group_b"] == "U"
    assert ok["contrast_id"] == "U_vs_C"
    meta, _ = _read(path)
    assert _parse_meta(meta)[0]["group_b"] == "U"
    # 同じ ds で 2 つの対比は別のファイルになる
    other = tmp_path / "t.tsv"
    json.loads(DAT.dataset_export_differential(str(other), contrast_id="T_vs_C"))
    assert other.read_text(encoding="utf-8") != path.read_text(encoding="utf-8")
    assert ds.last_differential["kind"] == "compare"


def test_export_interaction_contrast_is_not_two_group(tmp_path):
    ds = _install(_annotate(_ds(per_group=6, sex=True)), tmp_path)
    comparisons = [
        {"id": "T_vs_C_F", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "F"}},
        {"id": "inter", "kind": "interaction", "factors": ["group", "sex"],
         "levels": {"group": ["C", "T"], "sex": ["F", "M"]}}]
    json.loads(dataset_preprocess(**FLAT_OFF))
    out = json.loads(dataset_compare(comparisons=comparisons, factors=["group", "sex"]))
    assert out["status"] == "success"
    err = json.loads(DAT.dataset_export_differential(str(tmp_path / "x.tsv"),
                                                     contrast_id="inter"))["error"]
    assert err["code"] == "CONTRAST_NOT_TWO_GROUP" and err["details"]["two_group_contrast_ids"] == ["T_vs_C_F"]
    assert ds.last_differential["kind"] == "compare"


def test_export_after_exclusion_is_stale(ds, tmp_path):
    _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    apply_exclusions(ds, samples=["C0_r0"], features=[], mode="add", reason="outlier")
    with pytest.raises(DomainError) as exc:
        export_dataset_result(ds, record, tmp_path / "x.tsv")
    assert exc.value.code == "STALE_ANALYSIS_RESULT"
    assert not (tmp_path / "x.tsv").exists()
    # ツールは「直近」を捨てるので、古い結果を result_id で名指ししても missing_state（作り直しを案内）
    err = json.loads(DAT.dataset_export_differential(str(tmp_path / "y.tsv"),
                                                     result_id=record["result_id"]))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_compare"]
    assert err["state"] == "dataset_compare_result"


# --- 根拠バンドル v3 --------------------------------------------------------------------

def test_evidence_v3_from_compare_counts_match_samples_tsv(ds, tmp_path):
    out = _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    summary = evidence_export.export_dataset_evidence(ds, record, tmp_path / "b", polarity="positive")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    comparison = manifest["comparison"]
    samples = _tsv(tmp_path / "b" / "samples.tsv")
    assert manifest["contract_version"] == 3 and summary["contract_version"] == 3
    assert comparison["method"] == "moderated_t"
    assert comparison["n_a"] == sum(1 for r in samples if r["group"] == "a") == 8
    assert comparison["n_b"] == sum(1 for r in samples if r["group"] == "b") == 8
    two = stats_results.as_two_group(record)
    assert comparison["statistics"] == json.loads(json.dumps(two["statistics"]))
    assert comparison["statistics"]["n_units_a"] == 4 and comparison["statistics"]["subject_source"] == "subject_id"
    assert comparison["group_a"] == "C" and comparison["group_b"] == "T"
    assert manifest["source_analysis_id"] == out["result_id"]
    # QC は samples.tsv に残り、群は空
    assert {r["group"] for r in samples if r["sample_id"].startswith("QC")} == {""}


def test_evidence_v3_intensities_pre_imputation_from_matrix(ds, tmp_path):
    _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    summary = evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["intensities_stage"] == summary["intensities_stage"] == "pre_imputation"
    rows = _tsv(tmp_path / "b" / "intensities.tsv")
    assert len(rows) == 10 * 16                                   # 全特徴 × 比較の注入
    entry = stats_results.get_matrix(ds, record["matrix_id"])
    prepared = entry["prepared"]
    pos = {sid: i for i, sid in enumerate(prepared.sample_ids)}
    col = {fid: j for j, fid in enumerate(prepared.feature_ids)}
    names = stats_results.sample_ids_by_name(ds)
    cell = {(r["feature_id"], r["sample_id"]): r for r in rows}
    for fid, name in (("f0", "C0_r0"), ("f3", "T2_r1"), ("f9", "T0_r0")):
        expected = prepared.sample_values[pos[names[name]], col[fid]]
        assert float(cell[(fid, name)]["value"]) == pytest.approx(expected, rel=1e-8)
        assert cell[(fid, name)]["detected"] == "true"
    # 線形値（log2 ではない）: 試料の元の値（正規化なし）と同じ尺度
    raw = ds.feature_matrix[ds.feature_ids.index("f0"), ds.sample_names.index("C0_r0")]
    assert float(cell[("f0", "C0_r0")]["value"]) == pytest.approx(raw, rel=1e-8)
    assert float(cell[("f0", "C0_r0")]["value"]) > 1000


def test_evidence_v3_missing_matrix_returns_missing_state(ds, tmp_path):
    out = _compare(reference_group="C", test_group="T")
    ds.stat_matrices.clear()
    err = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), out["result_id"]))["error"]
    assert err["code"] == "missing_state" and err["state"] == "dataset_compare_result"
    assert err["required_tools"] == ["dataset_compare"]
    assert not (tmp_path / "b").exists()


def test_evidence_v3_multi_contrast_requires_contrast_id(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    err = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), out["result_id"]))["error"]
    assert err["code"] == "CONTRAST_REQUIRED"
    ok = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), out["result_id"],
                                                       contrast_id="U_vs_C"))
    assert ok["contract_version"] == 3
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["comparison"]["group_b"] == "U" and manifest["comparison"]["contrast_id"] == "U_vs_C"
    # 比較に入らない群 T の注入は samples.tsv で群が空
    rows = {r["sample_id"]: r for r in _tsv(tmp_path / "b" / "samples.tsv")}
    assert rows["T0_r0"]["group"] == "" and rows["U0_r0"]["group"] == "b" and rows["C0_r0"]["group"] == "a"
    assert ds.last_differential["kind"] == "compare"


# --- volcano ----------------------------------------------------------------------------

def test_volcano_from_compare_result_with_contrast_id(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    err = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))["error"]
    assert err["code"] == "CONTRAST_REQUIRED"
    assert "volcano" not in ds.figures
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload", contrast_id="T_vs_C",
                                                 result_id=out["result_id"]))
    assert payload["plot_schema"] == "lipidmix.volcano.v1"
    assert payload["comparison"]["group_a"] == "C" and payload["comparison"]["group_b"] == "T"
    assert payload["comparison"]["n_a"] == 8
    assert ds.figures["volcano"]["result_id"] == out["result_id"]
    assert payload["points"]
    # 画像も描ける（保存の記録は 2 群の形の結果を持つ）
    body = plot_tools.dataset_plot(kind="volcano", contrast_id="U_vs_C")
    assert isinstance(body, list) and "C" in body[0] and "U" in body[0]
    assert ds.figures["volcano"]["result"]["contrast_id"] == "U_vs_C"
    # contrast_id は volcano 以外の kind の引数ではない
    err = json.loads(plot_tools.dataset_plot(kind="pca", contrast_id="T_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"


def test_volcano_single_contrast_needs_no_contrast_id(ds):
    _compare(reference_group="C", test_group="T")
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    assert payload["comparison"]["group_b"] == "T"
    err = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload", contrast_id="zzz"))["error"]
    assert err["code"] == "CONTRAST_REQUIRED"


# --- 特徴検索 ---------------------------------------------------------------------------

def test_feature_search_stats_from_compare_and_merged_retest_uses_engine(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    rid = out["result_id"]
    chosen = feature_search.select_results(ds, None)
    assert [_r["provenance"]["result_id"] for _r in chosen["results"]] == [f"{rid}#T_vs_C", f"{rid}#U_vs_C"]
    assert chosen["stale"] == [] and chosen["missing"] == []
    # 名指し・contrast_id による絞り込み
    assert [r["contrast_id"] for r in feature_search.select_results(
        ds, [f"{rid}#U_vs_C"])["results"]] == ["U_vs_C"]
    assert [r["contrast_id"] for r in feature_search.select_results(
        ds, None, contrast_id="T_vs_C")["results"]] == ["T_vs_C"]
    assert feature_search.select_results(ds, [rid])["missing"] == [rid]      # 対比が複数なら ID だけでは曖昧
    # 統計の行は as_two_group の値
    rows = feature_search.feature_stats(ds, ["f0", "f5"], chosen["results"], with_labels=True)
    by = {(r["feature_id"], r["result_id"]): r for r in rows}
    two = stats_results.as_two_group(ds.results[rid], "U_vs_C")
    src = next(r for r in two["results"] if r["feature"] == "f0")
    got = by[("f0", f"{rid}#U_vs_C")]
    assert got["status"] == "tested" and got["q"] == pytest.approx(src["q"])
    assert got["log2fc"] == pytest.approx(src["log2fc"]) and got["mean_a"] == pytest.approx(src["mean_a"])
    assert got["group_a"] == "C" and got["group_b"] == "U" and got["comparison"] == "U_vs_C"
    assert got["detected_a"] == 8 and got["detected_b"] == 8
    # 割れた特徴（f1 と f2 を 1 つにまとめる）の再検定はエンジンで、記録の spec のまま
    groups = [["f1", "f2"]]
    merged = feature_search.merged_feature_stats(ds, ["f0", "f1", "f2"], chosen["results"], groups=groups)
    assert merged["renamed"] == {"f2": "f1"}
    redone = merged["results"][f"{rid}#U_vs_C"]
    assert redone["provenance"]["result_id"] == f"{rid}#U_vs_C" and redone["statistics"]["test"] == "moderated_t"
    ids = [r["feature"] for r in redone["results"]]
    assert "f2" not in ids and len(ids) == 9                    # 10 特徴 → まとめて 9
    assert redone["statistics"]["bh_population"]["n"] == 9      # BH はまとめた後の全特徴
    assert redone["samples_a"] == two["samples_a"] and redone["n_b"] == two["n_b"]
    # まとめない特徴の値は、まとめた行列で検定し直した値（記録の値とは一般に違う）。同じ spec なので検定は同じ
    f0_new = next(r for r in redone["results"] if r["feature"] == "f0")
    assert f0_new["log2fc"] == pytest.approx(src["log2fc"], rel=1e-6)      # f0 は割れの組に入らない
    assert {r["status"] for r in merged["stats"] if r["feature_id"] == "f0"} == {"tested"}
    # 元の記録は変えない
    assert len(ds.results[rid]["engine"]["rows"]) == 20
    # 同じ手順は同じ数字（決定論）
    again = feature_search.merged_feature_stats(ds, ["f0", "f1", "f2"], chosen["results"], groups=groups)
    assert again["stats"] == merged["stats"]


def test_feature_search_tool_contrast_and_stale_envelope(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    rid = out["result_id"]
    parsed = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                                   contrast_id="T_vs_C"))
    assert [r["result_id"] for r in parsed["results"]] == [rid + "#T_vs_C"] or \
        [r["result_id"] for r in parsed["results"]] == [rid]
    assert parsed["stats"].splitlines()[1].split("\t")[1] in (rid + "#T_vs_C", rid)
    both = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False))
    assert [r["result_id"] for r in both["results"]] == [f"{rid}#T_vs_C", f"{rid}#U_vs_C"]
    apply_exclusions(ds, samples=["C0_r0"], features=[], mode="add", reason="outlier")
    err = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], result_ids=[f"{rid}#T_vs_C"]))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_compare"]
    assert err["state"] == "dataset_compare_result"
    assert ds.results[rid]["kind"] == "compare"


# --- 新しい知見の候補 ---------------------------------------------------------------------

def _novel_ready(tmp_path, groups=("C", "T")):
    ds = _install(_annotate(_ds(groups=groups)), tmp_path)
    # 質量誤差を評価できる SME 候補（rank 1）を全特徴に付ける
    for k, fid in enumerate(ds.feature_ids):
        ds.feature_candidates[fid] = [{"rank": 1, "chemical_name": f"F{k}",
                                       "inchikey": ds.feature_metadata[fid]["inchikey"],
                                       "exp_mass_to_charge": 100.0 + k,
                                       "theoretical_mass_to_charge": 100.0 + k,
                                       "adduct_ion": "[M+H]1+"}]
    return ds


def test_novel_candidates_reads_compare_result(tmp_path):
    ds = _novel_ready(tmp_path)
    out = _compare(reference_group="C", test_group="T")
    spec = {"statistics": {"result_id": out["result_id"], "merge_split": False},
            "ppm": {"default": 5.0}, "q_max": 0.2}
    parsed = json.loads(reanalysis_tools.dataset_novel_candidates(spec))
    assert parsed.get("status") == "success", parsed
    assert parsed["group_ref"] == "C" and parsed["group_test"] == "T"
    assert parsed["n_rows"] >= 1
    header, *lines = parsed["rows"].splitlines()
    cols = header.split("\t")
    two = stats_results.as_two_group(ds.last_differential)
    q = {r["feature"]: r["q"] for r in two["results"]}
    for line in lines:
        cells = dict(zip(cols, line.split("\t")))
        assert float(cells["q"]) == pytest.approx(q[cells["feature_id"]], rel=1e-3)
        assert float(cells["q"]) < 0.2
    # 対比が複数の結果は `result_id#contrast_id` で名指しする
    ds2 = _novel_ready(tmp_path, groups=("C", "T", "U"))
    multi = _compare(groups=["C", "T", "U"])
    ambiguous = {"statistics": {"result_id": multi["result_id"], "merge_split": False},
                 "ppm": {"default": 5.0}, "q_max": 0.2}
    err = json.loads(reanalysis_tools.dataset_novel_candidates(ambiguous))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_compare"]
    named = {"statistics": {"result_id": f"{multi['result_id']}#U_vs_C", "merge_split": False},
             "ppm": {"default": 5.0}, "q_max": 0.2}
    ok = json.loads(reanalysis_tools.dataset_novel_candidates(named))
    assert ok["status"] == "success" and ok["group_test"] == "U"
    assert ds2.last_differential["kind"] == "compare"
