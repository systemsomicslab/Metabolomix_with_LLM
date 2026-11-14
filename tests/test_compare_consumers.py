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

FLAT_OFF = dict(drift_mode="off", blank_min_ratio=None, detection_min_rate=None,
                max_qc_rsd=None, max_d_ratio=None, normalize="none")

#: 書き出しが契約としてもともと持つメタ行のキー（これ以外を足さない）
V1_META_KEYS = {"contract_version", "exported_at", "source_mztab", "source_job", "id_space", "result_id",
                "preprocess_id", "source_verification", "group_a", "n_a", "group_b", "n_b", "log2fc_sign",
                "q_threshold", "log2fc_threshold", "log_transform", "preprocess", "n_features_total",
                "n_with_inchikey", "n_unannotated"}
NEW_META_KEYS = {"test", "effect_size_definition", "ci_level", "verdict_q_threshold", "equivalence_fc",
                 "bh_population", "bh_n", "q_interpretation", "p_adjust_method", "unit", "n_units_a", "n_units_b",
                 "mean_definition", "matrix_id", "rule_hash", "contrast_id", "significant_rule",
                 "q_value_method", "q_value_source"}
#: 条件つきのメタ行（交互作用・探索専用・検定後の除外・BH の q を計算し直したときだけ）。どれも msi_note の後ろに足す。
OPTIONAL_META_KEYS = {"contrast_kind", "interaction_factors", "interaction_levels", "interaction_definition",
                      "exploratory", "exploratory_reason", "n_excluded_after_test", "excluded_after_test_note",
                      "q_value_note"}


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
    assert msi == len(meta) - 1 - 7
    values, _ = _parse_meta(meta[msi + 1:])
    assert values["test"] == "moderated_t" and values["effect_size_definition"] == "log2_mean_difference"
    assert values["ci_level"] == "0.95" and values["verdict_q_threshold"] == "0.05"
    assert values["bh_population"] == "matrix" and values["bh_n"] == "10"
    assert values["q_interpretation"] == "discovery_fdr" and values["p_adjust_method"] == "bh"
    # q_value は契約の BH の q（比較が BH なのでそのまま）
    assert values["q_value_method"] == "bh" and values["q_value_source"] == "comparison"
    assert "q_value_note" not in values and info["q_value_source"] == "comparison" and "notes" not in info
    assert values["unit"] == "subject" and values["n_units_a"] == "4" and values["n_units_b"] == "4"
    assert values["mean_definition"] == "linear arithmetic mean over units after preprocessing"
    assert values["matrix_id"] == out["matrix_id"] and values["rule_hash"] == out["rule_hash"]
    # n_a / n_b は注入の数（4 個体 × 技術反復 2）、単位の数は別に書く
    head, _ = _parse_meta(meta)
    assert head["n_a"] == "8" and head["n_b"] == "8" and head["contract_version"] == "1"
    assert head["log_transform"] == "true" and head["preprocess_id"] == out["matrix_id"]
    assert values["significant_rule"] == export_contract.SIGNIFICANT_RULE
    # significant は契約の規則（q ≤ しきい値 かつ |log2fc| ≥ log2fc_threshold）に、エンジンの判定が changed（R16）
    two = stats_results.as_two_group(ds.last_differential)
    expected = {row["feature"]: export_contract.is_significant(
        q=row["q"], log2fc=row["log2fc"], q_threshold=two["q_threshold"],
        log2fc_threshold=two["log2fc_threshold"]) and row["verdict"] == "changed"
        for row in two["results"]}
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
    assert not (NEW_META_KEYS & V1_META_KEYS) and not (OPTIONAL_META_KEYS & (V1_META_KEYS | NEW_META_KEYS))
    # 新しい行は「# key = value」の形で、値にタブ・改行を含まない
    assert all(line.startswith("# ") for line in meta)


def _asks_for_contrast(body: str, tool: str, ids: list[str]) -> dict:
    """対比を省いた曖昧さは誤りではなく needs_input（選択肢は対比の一覧、呼び直すのは同じツール）。"""
    parsed = json.loads(body)
    assert "error" not in parsed
    ask = parsed["needs_input"]
    assert ask["field"] == "contrast_id" and ask["argument"] == "contrast_id"
    assert [o["value"] for o in ask["options"]] == ids and all(o["label"] for o in ask["options"])
    assert ask["retry_tools"] == [tool]
    return ask


def test_export_multi_contrast_asks_for_contrast_id(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    assert [c["contrast_id"] for c in out["contrasts"]] == ["T_vs_C", "U_vs_C"]
    path = tmp_path / "out.tsv"
    _asks_for_contrast(DAT.dataset_export_differential(str(path)), "dataset_export_differential",
                       ["T_vs_C", "U_vs_C"])
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


def _interaction_ds(tmp_path):
    ds = _install(_annotate(_ds(per_group=6, sex=True)), tmp_path)
    comparisons = [
        {"id": "T_vs_C_F", "factor": "group", "reference": "C", "test": "T", "at": {"sex": "F"}},
        {"id": "inter", "kind": "interaction", "factors": ["group", "sex"],
         "levels": {"group": ["C", "T"], "sex": ["F", "M"]}}]
    json.loads(dataset_preprocess(**FLAT_OFF))
    out = json.loads(dataset_compare(comparisons=comparisons, factors=["group", "sex"]))
    assert out["status"] == "success"
    return ds, out


def test_export_interaction_contrast_writes_the_difference_of_differences(tmp_path):
    """2×2 の交互作用（差の差）は名指しすれば書き出せる。列は 15 のまま、mean は空欄、定義はメタ行。"""
    ds, out = _interaction_ds(tmp_path)
    path = tmp_path / "x.tsv"
    info = json.loads(DAT.dataset_export_differential(str(path), contrast_id="inter"))
    assert info["status"] == "success" and info["contrast_id"] == "inter"
    assert info["contrast_kind"] == "interaction"
    assert info["group_a"] == "C:M + T:F" and info["group_b"] == "T:M + C:F"
    meta, body = _read(path)
    assert body[0].split("\t") == export_contract.EXPORT_COLUMNS
    values, keys = _parse_meta(meta)
    assert len(keys) == len(set(keys))
    assert values["contrast_kind"] == "interaction" and values["interaction_factors"] == "group x sex"
    # 水準は「基準 → 比較」を言葉で書く（`>` は大小と読めるので使わない）
    assert values["interaction_levels"] == "group reference C test T; sex reference F test M"
    assert ">" not in values["interaction_levels"]
    assert "minus" in values["interaction_definition"]
    assert values["mean_definition"] == export_contract.INTERACTION_MEAN_DEFINITION
    assert set(values) <= V1_META_KEYS | NEW_META_KEYS | OPTIONAL_META_KEYS
    # log2fc はエンジンの交互作用の行そのもの（差の差）、mean_a / mean_b は空欄
    rows = {r["feature_id"]: r for r in ds.last_differential["engine"]["rows"] if r["contrast_id"] == "inter"}
    cells = {line.split("\t")[0]: line.split("\t") for line in body[1:]}
    col = {name: i for i, name in enumerate(export_contract.EXPORT_COLUMNS)}
    assert cells
    for fid, cell in cells.items():
        assert cell[col["mean_a"]] == "" and cell[col["mean_b"]] == ""
        if rows[fid]["log2fc"] is not None:
            assert float(cell[col["log2fc"]]) == pytest.approx(rows[fid]["log2fc"], abs=1e-4)
    # 2 群の形でしか読まない利用者（volcano・根拠バンドル）は交互作用を断る
    err = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload", contrast_id="inter"))["error"]
    assert err["code"] == "CONTRAST_NOT_TWO_GROUP" and err["details"]["two_group_contrast_ids"] == ["T_vs_C_F"]
    err = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), contrast_id="inter"))["error"]
    assert err["code"] == "CONTRAST_NOT_TWO_GROUP"
    assert ds.last_differential["kind"] == "compare"


def test_export_with_an_interaction_and_one_two_group_contrast_asks_which(tmp_path):
    _interaction_ds(tmp_path)
    ask = _asks_for_contrast(DAT.dataset_export_differential(str(tmp_path / "x.tsv")),
                             "dataset_export_differential", ["T_vs_C_F", "inter"])
    assert "interaction" in ask["options"][1]["label"]
    # 2 群の形で読む volcano は、2 群の対比が 1 つだけなので黙ってそれを使う
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    assert payload["comparison"]["group_b"] == "T"


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
    shown = {k: v for k, v in two["statistics"].items() if k != "p_adjust_method"}   # 補正は correction に出す
    assert comparison["statistics"] == json.loads(json.dumps(shown))
    assert comparison["correction"] == "BH"
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


def test_evidence_v3_multi_contrast_asks_for_contrast_id(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    out = _compare(groups=["C", "T", "U"])
    _asks_for_contrast(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), out["result_id"]),
                       "dataset_export_evidence_bundle", ["T_vs_C", "U_vs_C"])
    assert not (tmp_path / "b").exists()
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
    _asks_for_contrast(plot_tools.dataset_plot(kind="volcano", output="payload"), "dataset_plot",
                       ["T_vs_C", "U_vs_C"])
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
    # 対比が複数の記録を ID だけで名指しすれば全対比に展開する（単一の対比が要る呼び出しは single=True で拒む）
    assert [r["contrast_id"] for r in feature_search.select_results(ds, [rid])["results"]] == [
        "T_vs_C", "U_vs_C"]
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
    assert [r["result_id"] for r in parsed["results"]] == [rid + "#T_vs_C"]
    assert parsed["stats"].splitlines()[1].split("\t")[1] == rid + "#T_vs_C"
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
    # ID だけなら誤りでなく問い（spec の statistics.result_id に result_id#contrast_id を入れて呼び直す）
    asked = json.loads(reanalysis_tools.dataset_novel_candidates(ambiguous))
    assert "error" not in asked
    body = asked["needs_input"]
    assert body["field"] == "statistics.result_id" and "argument" not in body
    assert body["retry_tools"] == ["dataset_novel_candidates"]
    assert [o["value"] for o in body["options"]] == [f"{multi['result_id']}#T_vs_C",
                                                    f"{multi['result_id']}#U_vs_C"]
    assert all(o["label"] for o in body["options"])
    # 未知の対比の名指しは誤りのまま（候補つき）
    unknown = {"statistics": {"result_id": f"{multi['result_id']}#nope", "merge_split": False},
               "ppm": {"default": 5.0}, "q_max": 0.2}
    err = json.loads(reanalysis_tools.dataset_novel_candidates(unknown))["error"]
    assert err["code"] == "CONTRAST_REQUIRED" and err["details"]["two_group_contrast_ids"] == ["T_vs_C", "U_vs_C"]
    named = {"statistics": {"result_id": f"{multi['result_id']}#U_vs_C", "merge_split": False},
             "ppm": {"default": 5.0}, "q_max": 0.2}
    ok = json.loads(reanalysis_tools.dataset_novel_candidates(named))
    assert ok["status"] == "success" and ok["group_test"] == "U"
    assert ds2.last_differential["kind"] == "compare"


# --- fix round 1 ---------------------------------------------------------------------------

def _multi(tmp_path):
    ds = _install(_annotate(_ds(groups=("C", "T", "U"))), tmp_path)
    return ds, _compare(groups=["C", "T", "U"])["result_id"]


def test_select_results_contrast_id_rules(tmp_path):
    ds, rid = _multi(tmp_path)
    # (a) ID だけ + contrast_id は rid#contrast_id と同じ
    sel = feature_search.select_results(ds, [rid], contrast_id="T_vs_C")
    assert [r["provenance"]["result_id"] for r in sel["results"]] == [f"{rid}#T_vs_C"]
    assert sel["missing"] == [] and sel["stale"] == []
    # (b) 名指しの対比と contrast_id が食い違えば DATASET_BAD_REQUEST（黙って落とさない）
    with pytest.raises(DomainError) as exc:
        feature_search.select_results(ds, [f"{rid}#U_vs_C"], contrast_id="T_vs_C")
    assert exc.value.code == "DATASET_BAD_REQUEST"
    ok = feature_search.select_results(ds, [f"{rid}#U_vs_C"], contrast_id="U_vs_C")
    assert [r["contrast_id"] for r in ok["results"]] == ["U_vs_C"]
    # (c) 未知の contrast_id は、結果を名指ししても省いても CONTRAST_REQUIRED（候補つき）
    for ids in ([rid], None):
        with pytest.raises(DomainError) as exc:
            feature_search.select_results(ds, ids, contrast_id="nope")
        assert exc.value.code == "CONTRAST_REQUIRED"
        assert set(exc.value.details["contrast_ids"]) >= {"T_vs_C", "U_vs_C"}
    # (d) 単一の対比が要る呼び出しは、複数の記録を ID だけで名指しすると CONTRAST_REQUIRED
    with pytest.raises(DomainError) as exc:
        feature_search.select_results(ds, [rid], single=True)
    assert exc.value.code == "CONTRAST_REQUIRED"
    assert exc.value.details["two_group_contrast_ids"] == ["T_vs_C", "U_vs_C"]
    assert feature_search.select_results(ds, [f"{rid}#T_vs_C"], single=True)["results"]


def test_select_results_stale_lists_named_id_once(tmp_path):
    ds, rid = _multi(tmp_path)
    apply_exclusions(ds, samples=["C0_r0"], features=[], mode="add", reason="x")
    named = f"{rid}#T_vs_C"
    sel = feature_search.select_results(ds, [named, named], contrast_id="T_vs_C")
    assert sel["stale"] == [named] and sel["results"] == []
    assert feature_search.select_results(ds, [rid], contrast_id="T_vs_C")["stale"] == [rid]


def test_feature_search_tool_contrast_errors_are_envelopes(tmp_path):
    ds, rid = _multi(tmp_path)
    ok = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                               result_ids=[rid], contrast_id="T_vs_C"))
    assert [r["result_id"] for r in ok["results"]] == [f"{rid}#T_vs_C"]
    err = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                                result_ids=[f"{rid}#U_vs_C"], contrast_id="T_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    for ids in (None, [rid]):
        err = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                                    result_ids=ids, contrast_id="nope"))["error"]
        assert err["code"] == "CONTRAST_REQUIRED"


def test_feature_search_ambiguous_contrast_is_needs_input(tmp_path, monkeypatch):
    """対比を省いた曖昧さ（`as_two_group` の ambiguous な `CONTRAST_REQUIRED`）は誤りでなく問い。

    ほかの対比を選ぶツール（書き出し・volcano・根拠バンドル）と同じ `contrast_envelope` を通す。
    `select_results` は今は対比ごとに広げるのでこの曖昧さを出さないが、出たときの封筒の形を縛る。
    """
    ds, rid = _multi(tmp_path)

    def _ambiguous(ds_, result_ids, **kwargs):
        stats_results.as_two_group(ds_.results[rid])      # 対比が 2 つ → ambiguous な CONTRAST_REQUIRED

    monkeypatch.setattr(feature_search, "select_results", _ambiguous)
    asked = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                                  result_ids=[rid]))
    assert "error" not in asked
    body = asked["needs_input"]
    assert body["field"] == "contrast_id" and body["argument"] == "contrast_id"
    assert body["retry_tools"] == ["dataset_feature_search"]
    assert [o["value"] for o in body["options"]] == ["T_vs_C", "U_vs_C"]


def test_rid_hash_cid_accepted_by_export_evidence_and_volcano(tmp_path):
    ds, rid = _multi(tmp_path)
    named = f"{rid}#U_vs_C"
    out = json.loads(DAT.dataset_export_differential(str(tmp_path / "e.tsv"), result_id=named))
    assert out["status"] == "success" and out["contrast_id"] == "U_vs_C"
    assert json.loads(DAT.dataset_export_differential(str(tmp_path / "e2.tsv"), result_id=named,
                                                      contrast_id="U_vs_C"))["status"] == "success"
    err = json.loads(DAT.dataset_export_differential(str(tmp_path / "e3.tsv"), result_id=named,
                                                     contrast_id="T_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    ev = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b"), named))
    assert ev["contract_version"] == 3
    err = json.loads(DAT.dataset_export_evidence_bundle(str(tmp_path / "b2"), named,
                                                        contrast_id="T_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload", result_id=named))
    assert payload["comparison"]["group_b"] == "U"
    err = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload", result_id=named,
                                             contrast_id="T_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"


def test_evidence_manifest_carries_engine_definitions(ds, tmp_path):
    out = _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    comparison = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))["comparison"]
    assert comparison["effect_size_definition"]["name"] == "log2_mean_difference"
    assert comparison["mean_definition"] == export_contract.MEAN_DEFINITION
    assert comparison["matrix_id"] == out["matrix_id"] and comparison["rule_hash"] == out["rule_hash"]
    assert comparison["contrast_id"] == "T_vs_C"


def test_evidence_population_exclusions_from_matrix_entry(ds, tmp_path):
    apply_exclusions(ds, samples=["C0_r0", "C0_r1"], features=["f9"], mode="add", reason="x")
    _compare(reference_group="C", test_group="T")
    evidence_export.export_dataset_evidence(ds, ds.last_differential, tmp_path / "b")
    pop = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))["population"]
    by = {e["stage"]: e for e in pop["exclusions"]}
    assert by["excluded_samples"]["removed"] == 2 and by["excluded_samples"]["names"] == ["C0_r0", "C0_r1"]
    assert by["excluded_features"]["names"] == ["f9"]
    assert pop["derivation_status"] == "not_available"


@pytest.mark.parametrize("source,unit", [("subject_id", "subject"), ("injection", "injection"),
                                         ("mixed", "mixed"), (None, "injection"), ("weird", "injection")])
def test_export_unit_meta_values(ds, tmp_path, source, unit):
    _compare(reference_group="C", test_group="T")
    ds.last_differential["engine"]["subject_source"] = source
    path = tmp_path / "u.tsv"
    export_dataset_result(ds, ds.last_differential, path)
    meta, _ = _read(path)
    assert _parse_meta(meta)[0]["unit"] == unit
    assert "unit = None" not in "\n".join(meta)


def test_export_meta_has_contrast_id_on_matrix_line(tmp_path):
    ds, rid = _multi(tmp_path)
    path = tmp_path / "m.tsv"
    export_dataset_result(ds, ds.results[rid], path, contrast_id="U_vs_C")
    meta, _ = _read(path)
    line = next(m for m in meta if m.startswith("# matrix_id"))
    assert [c.split(" = ")[0] for c in line.lstrip("# ").split("\t")] == ["matrix_id", "rule_hash", "contrast_id"]
    values, keys = _parse_meta(meta)
    assert values["contrast_id"] == "U_vs_C" and len(keys) == len(set(keys))


def test_meta_source_has_no_literal_tab_and_dataset_plot_signature():
    import inspect
    from pathlib import Path
    src = Path(export_contract.__file__).read_text(encoding="utf-8")
    assert "\t" not in src
    params = list(inspect.signature(plot_tools.dataset_plot).parameters)
    assert params[:4] == ["kind", "dataset_handle", "result_id", "output"]
    assert params[-2:] == ["contrast_id", "output_path"]


# --- 最終レビューの修正（R16: 判定できない行を有意にしない / 交絡の許可の痕跡） ---------------------

def _not_judgeable(ds, feature="f0", group="T"):
    """`feature` を `group` の注入で実検出なし（gap-fill の値だけ）にする。値は残るので量の q は小さいが、
    検出の割合が足りず判定は inconclusive（quantity_not_judgeable）になる。"""
    row = ds.feature_ids.index(feature)
    for j, name in enumerate(ds.sample_names):
        if name.startswith(group):
            ds.detected_mask[row, j] = False
    return ds


def _body_cells(path) -> dict[str, list[str]]:
    _, body = _read(path)
    return {cells[0]: cells for cells in (line.split("\t") for line in body[1:])}


def test_engine_not_judgeable_row_is_not_significant_anywhere(tmp_path):
    ds = _install(_not_judgeable(_annotate(_ds())), tmp_path)
    out = _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    f0 = next(r for r in record["engine"]["rows"] if r["feature_id"] == "f0")
    assert (f0["verdict"], f0["verdict_reason"]) == ("inconclusive", "quantity_not_judgeable")
    assert f0["q"] < 1e-6 and abs(f0["log2fc"]) > 1.5           # q と |log2fc| だけなら有意に見える
    two = stats_results.as_two_group(record)
    flags = {r["feature"]: r["significant"] for r in two["results"]}
    assert flags["f0"] is False and flags["f1"] is True and flags["f2"] is True
    large = out["contrasts"][0]["large_changes"]
    assert large == 2 == sum(flags.values())
    # 書き出し: significant 列は判定と食い違わない。true の数は要約の large_changes と同じ
    path = tmp_path / "out.tsv"
    export_dataset_result(ds, record, path)
    cells = _body_cells(path)
    assert cells["f0"][14] == "false" and cells["f1"][14] == "true"
    assert sum(c[14] == "true" for c in cells.values()) == large
    meta, _ = _read(path)
    values, keys = _parse_meta(meta)
    assert values["significant_rule"] == "q<=q_threshold & |log2fc|>=log2fc_threshold & verdict=changed"
    assert len(keys) == len(set(keys))
    # volcano: 色の印も同じ（新しい sig の種類は足さない）
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    sig = {p["feature"]: p["sig"] for p in payload["points"]}
    assert sig["f0"] == "ns" and sig["f1"] == "up"
    assert {p["sig"] for p in payload["points"]} <= {"up", "down", "ns"}
    assert payload["selection"]["significant_total"] == large
    # 根拠バンドル: features.tsv の significant も同じ
    evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    features = {r["feature_id"]: r for r in _tsv(tmp_path / "b" / "features.tsv")}
    assert features["f0"]["significant"] == "false" and features["f1"]["significant"] == "true"
    assert sum(r["significant"] == "true" for r in features.values()) == large


def test_engine_significant_flag_is_produced_once_by_as_two_group(ds):
    _compare(reference_group="C", test_group="T")
    two = stats_results.as_two_group(ds.last_differential)
    for row in two["results"]:
        expected = (export_contract.is_significant(q=row["q"], log2fc=row["log2fc"],
                                                   q_threshold=two["q_threshold"],
                                                   log2fc_threshold=two["log2fc_threshold"])
                    and row["verdict"] == "changed")
        assert row["significant"] is expected
        # 消費側は行の印を読む（閾値を変えても行の印が勝つ）
        assert export_contract.row_significant(row, q_threshold=1.0, log2fc_threshold=0.0) is expected


def test_novel_candidates_skip_engine_rows_that_are_not_changed(tmp_path):
    ds = _not_judgeable(_novel_ready(tmp_path))
    out = _compare(reference_group="C", test_group="T")
    spec = {"statistics": {"result_id": out["result_id"], "merge_split": False},
            "ppm": {"default": 5.0}, "q_max": 0.2}
    parsed = json.loads(reanalysis_tools.dataset_novel_candidates(spec))
    assert parsed.get("status") == "success", parsed
    header, *lines = parsed["rows"].splitlines()
    ids = [dict(zip(header.split("\t"), line.split("\t")))["feature_id"] for line in lines]
    assert "f0" not in ids and {"f1", "f2"} <= set(ids)          # f0 は q が小さいが判定は inconclusive
    assert ds.last_differential["kind"] == "compare"


def _confounded(ds, batches=("B1", "B2")):
    """C 群の注入を batches[0]、T 群を batches[1] に置く（QC は両方に分ける）。"""
    from metabolomix.analysis.sample_manifest import apply_metadata
    rows = []
    for k, row in enumerate(ds.sample_metadata_rows):
        r = dict(row)
        name = ds.sample_names[k]
        r["batch"] = (batches[0] if name.startswith("C") else batches[1] if name.startswith("T")
                      else batches[k % 2])
        rows.append(r)
    apply_metadata(ds, rows)
    return ds


def test_allow_confounded_trace_reaches_export_volcano_and_evidence(tmp_path):
    from metabolomix.plots.result_output import _UNADJUSTED_CONFOUNDED_NOTE, figure_annotations
    ds = _install(_confounded(_annotate(_ds())), tmp_path)
    pp = json.loads(dataset_preprocess(**FLAT_OFF))
    err = json.loads(dataset_compare(matrix_id=pp["matrix_id"], reference_group="C", test_group="T",
                                     allow_confounded=False))
    assert err["error"]["code"] == "STATISTIC_CONFOUNDED"
    out = json.loads(dataset_compare(matrix_id=pp["matrix_id"], reference_group="C", test_group="T",
                                     allow_confounded=True))
    assert out["status"] == "success" and out["confounding"]["batch_status"] == "complete"
    record = ds.last_differential
    view = stats_results.as_two_group(record)
    comparison = view["provenance"]["comparison"]
    assert comparison["unadjusted_confounded"] is True and comparison["allow_confounded"] is True
    assert comparison["confounding"]["batch_status"] == "complete"
    assert "comparison" not in record["provenance"]               # 記録そのものは書き換えない
    # 書き出しのメタ行
    path = tmp_path / "out.tsv"
    export_dataset_result(ds, record, path)
    meta, _ = _read(path)
    values, _ = _parse_meta(meta)
    assert values["unadjusted_confounded"] == "true"
    assert "未調整" in values["unadjusted_confounded_note"]
    # volcano の payload と図の但し書き
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    assert any("allow_confounded=true" in c for c in payload["caveats"])
    assert _UNADJUSTED_CONFOUNDED_NOTE in figure_annotations(ds, ds.figures["volcano"]["result"])
    # 根拠バンドルの manifest（comparison と、v1 と同じ来歴の comparison）
    evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["comparison"]["unadjusted_confounded"] is True
    assert manifest["provenance"]["result"]["comparison"]["unadjusted_confounded"] is True


def test_no_confound_trace_without_complete_confounding(ds, tmp_path):
    from metabolomix.plots.result_output import figure_annotations
    out = _compare(reference_group="C", test_group="T", allow_confounded=True)   # 許可しても交絡が無い
    assert out["confounding"]["batch_status"] != "complete"
    record = ds.last_differential
    view = stats_results.as_two_group(record)
    assert "comparison" not in view["provenance"]
    path = tmp_path / "out.tsv"
    export_dataset_result(ds, record, path)
    assert "unadjusted_confounded" not in path.read_text(encoding="utf-8")
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    assert not any("allow_confounded=true" in c for c in payload["caveats"])
    assert figure_annotations(ds, ds.figures["volcano"]["result"]) == []
    evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert "unadjusted_confounded" not in manifest["comparison"]
    assert "comparison" not in manifest["provenance"]["result"]


def test_select_results_interaction_only_record_is_reported_missing(tmp_path):
    ds = _install(_annotate(_ds(per_group=6, sex=True)), tmp_path)
    inter = [{"id": "inter", "kind": "interaction", "factors": ["group", "sex"],
              "levels": {"group": ["C", "T"], "sex": ["F", "M"]}}]
    json.loads(dataset_preprocess(**FLAT_OFF))
    out = json.loads(dataset_compare(comparisons=inter, factors=["group", "sex"]))
    rid = out["result_id"]
    for ids in (None, [rid]):
        sel = feature_search.select_results(ds, ids)
        assert sel["results"] == [] and len(sel["missing"]) == 1
        assert sel["missing"][0].startswith(f"{rid}: ") and "inter" in sel["missing"][0]
    parsed = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False))
    assert any(rid in note for note in parsed["notes"])
    err = json.loads(DAT.dataset_feature_search(feature_ids=["f0"], merge_split=False,
                                                result_ids=[rid]))["error"]
    assert err["code"] == "missing_state" and "inter" in err["message"]

# --- 検定の後の除外・補正の方法・save_figure の対比（担当 7） ---------------------------------

def test_export_exclude_feature_ids_drops_rows_after_testing(ds, tmp_path):
    _compare(reference_group="C", test_group="T")
    full = tmp_path / "full.tsv"
    DAT.dataset_export_differential(str(full))
    meta_full, body_full = _read(full)
    assert not any("n_excluded_after_test" in m for m in meta_full)
    cut = tmp_path / "cut.tsv"
    info = json.loads(DAT.dataset_export_differential(str(cut), exclude_feature_ids=["f0", "f3", "nope"]))
    assert info["status"] == "success" and info["n_excluded_after_test"] == 2
    assert info["exclude_feature_ids_not_found"] == ["nope"]
    meta, body = _read(cut)
    values, keys = _parse_meta(meta)
    assert len(keys) == len(set(keys)) and set(values) <= V1_META_KEYS | NEW_META_KEYS | OPTIONAL_META_KEYS
    assert values["n_excluded_after_test"] == "2" and "not recomputed" in values["excluded_after_test_note"]
    # 注記は補正の方法・母集団によらず正しい書き方（BH と決めつけない）で、件数の関係も書く
    note = values["excluded_after_test_note"]
    assert "BH" not in note and "multiple-testing population is unchanged" in note
    assert "n_features_total still counts them" in note
    assert (int(values["n_features_total"]) == info["n_features_total"]
            == int(values["n_with_inchikey"]) + int(values["n_unannotated"]) + int(values["n_excluded_after_test"]))
    full_rows = {line.split("\t")[0]: line for line in body_full[1:]}
    cut_rows = {line.split("\t")[0]: line for line in body[1:]}
    assert set(full_rows) - set(cut_rows) == {"f0", "f3"}
    # q は計算し直さない（残った行は落とす前と同じ）
    assert all(cut_rows[fid] == full_rows[fid] for fid in cut_rows)
    # 一覧でなければ引数の誤り
    err = json.loads(DAT.dataset_export_differential(str(tmp_path / "x.tsv"), exclude_feature_ids="f0"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"


def _hand_bh(p_by_key: dict) -> dict:
    """Benjamini–Hochberg の q を手で（numpy・エンジンを使わず）: 昇順の i 番目は min_{k>=i} p_(k)·m/k を 1 で頭打ち。"""
    ordered = sorted(p_by_key.items(), key=lambda kv: kv[1])
    m = len(ordered)
    out, running = {}, 1.0
    for rank in range(m, 0, -1):
        key, p = ordered[rank - 1]
        running = min(running, p * m / rank)
        out[key] = running
    return out


def _export_cells(path) -> dict:
    meta, body = _read(path)
    col = {name: i for i, name in enumerate(export_contract.EXPORT_COLUMNS)}
    cells = {line.split("\t")[0]: line.split("\t") for line in body[1:]}
    return _parse_meta(meta)[0], {fid: (float(c[col["q_value"]]) if c[col["q_value"]] else None,
                                        c[col["significant"]] == "true") for fid, c in cells.items()}


@pytest.mark.parametrize("method", ["none", "holm", "bonferroni"])
def test_export_and_bundle_write_the_contract_bh_q_for_a_non_bh_comparison(ds, tmp_path, method):
    """R5/R21: 契約の q_value は BH の q（FDR）。比較が別の補正でも、同じ母集団で p から計算し直し、
    significant もその q で読む（比較の補正の方法はメタ行 p_adjust_method に残す）。"""
    # 中くらいの差を 2 つ作る: f4（p≈0.009。BH では有意、holm / bonferroni では有意でない）と
    # f5（p≈0.033。補正しない p では有意、BH では有意でない）。log2fc の閾値は 0（判定だけで印が決まる）。
    for j, name in enumerate(ds.sample_names):
        if name.startswith("T"):
            ds.feature_matrix[4, j] *= 2 ** 0.08
            ds.feature_matrix[5, j] *= 2 ** 0.27
    bh_out = _compare(reference_group="C", test_group="T", p_adjust="bh", log2fc_threshold=0.0)
    bh_record = ds.results[bh_out["result_id"]]
    other_out = _compare(reference_group="C", test_group="T", p_adjust=method, log2fc_threshold=0.0)
    record = ds.results[other_out["result_id"]]
    assert record["p_adjust_method"] == method
    rows = [r for r in record["engine"]["rows"] if r["contrast_id"] == "T_vs_C"]
    # 手で計算した BH（母集団は比較の q が付いた行。ここでは行列の全特徴）
    population = {r["feature_id"]: r["p"] for r in rows if r["q"] is not None}
    assert len(population) == record["engine"]["bh_population"]["n_by_contrast"]["T_vs_C"]
    expected_q = _hand_bh(population)
    # エンジンが p_adjust="bh" で出した q と一致する（同じ母集団の BH）
    bh_rows = {r["feature_id"]: r for r in bh_record["engine"]["rows"] if r["contrast_id"] == "T_vs_C"}
    assert all(bh_rows[f]["q"] == pytest.approx(q, rel=1e-12) for f, q in expected_q.items())
    bh_significant = {r["feature"]: r["significant"] for r in stats_results.as_two_group(bh_record)["results"]}
    own_significant = {r["feature"]: r["significant"] for r in stats_results.as_two_group(record)["results"]}
    # 比較の補正の印は BH の印と違う（テストが区別できること）: none は f5 を、holm / bonferroni は f4 を変える
    assert own_significant != bh_significant
    assert bh_significant["f4"] is True and bh_significant["f5"] is False

    path = tmp_path / f"{method}.tsv"
    info = export_dataset_result(ds, record, path)
    values, cells = _export_cells(path)
    assert values["p_adjust_method"] == method                    # 比較の補正の方法は残す
    assert values["q_value_method"] == "bh" and values["q_value_source"] == "recomputed_from_p_value"
    assert method in values["q_value_note"] and "Benjamini-Hochberg" in values["q_value_note"]
    assert values["q_interpretation"] == "discovery_fdr"           # q_value の読み方（BH）
    assert values["bh_population"] == "matrix" and values["bh_n"] == str(len(population))
    assert info["p_adjust_method"] == method and info["q_value_source"] == "recomputed_from_p_value"
    assert info["notes"] and method in info["notes"][0]
    for fid, (q, significant) in cells.items():
        assert q == pytest.approx(expected_q[fid], rel=1e-5)      # 書式は有効数字 6 桁
        assert significant == bh_significant[fid]                 # BH の比較と同じ印
    # BH の比較の書き出しと q_value・significant が一致する
    export_dataset_result(ds, bh_record, tmp_path / "bh.tsv")
    _, bh_cells = _export_cells(tmp_path / "bh.tsv")
    assert cells == bh_cells

    # 根拠バンドルも同じ（features.tsv の q_value は BH、correction は BH、比較の方法は comparison に足す）
    summary = evidence_export.export_dataset_evidence(ds, record, tmp_path / f"b-{method}")
    manifest = json.loads((tmp_path / f"b-{method}" / "manifest.json").read_text(encoding="utf-8"))
    comparison = manifest["comparison"]
    assert comparison["correction"] == "BH"
    assert comparison["p_adjust_method"] == method and comparison["q_value_source"] == "recomputed_from_p_value"
    assert comparison["statistics"]["bh_population"]["q_interpretation"] == "discovery_fdr"
    assert "p_adjust_method" not in comparison["statistics"]
    assert summary["q_value_source"] == "recomputed_from_p_value" and summary["notes"]
    features = {r["feature_id"]: r for r in _tsv(tmp_path / f"b-{method}" / "features.tsv")}
    for fid, row in features.items():
        assert float(row["q_value"]) == pytest.approx(expected_q[fid], rel=1e-9)
        assert (row["significant"] == "true") == bh_significant[fid]


def test_bh_comparison_export_and_bundle_are_unchanged_by_the_contract_q(ds, tmp_path):
    """BH の比較は q をそのまま書き、manifest に比較の方法の欄を足さない（共有の見本は変わらない）。"""
    _compare(reference_group="C", test_group="T")
    record = ds.last_differential
    info = export_dataset_result(ds, record, tmp_path / "x.tsv")
    values, cells = _export_cells(tmp_path / "x.tsv")
    assert values["q_value_method"] == "bh" and values["q_value_source"] == "comparison"
    assert "notes" not in info
    own = {r["feature_id"]: r["q"] for r in record["engine"]["rows"]}
    assert all(q == pytest.approx(own[fid], rel=1e-5) for fid, (q, _) in cells.items())
    summary = evidence_export.export_dataset_evidence(ds, record, tmp_path / "b")
    comparison = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))["comparison"]
    assert comparison["correction"] == "BH" and "p_adjust_method" not in comparison
    assert "q_value_source" not in comparison and "notes" not in summary
    # 補正の方法を記録する前の記録（宣言が無い）も BH
    record.pop("p_adjust_method", None)
    export_dataset_result(ds, record, tmp_path / "y.tsv")
    values, _ = _export_cells(tmp_path / "y.tsv")
    assert values["p_adjust_method"] == "bh" and values["q_value_source"] == "comparison"
    assert values["q_interpretation"] == "discovery_fdr"


def test_save_figure_contrast_id_redraws_the_chosen_contrast(tmp_path):
    from metabolomix.tools import reports

    ds, rid = _multi(tmp_path)
    figures = tmp_path / "analysis" / "reports" / "figures"
    msg = reports.save_figure("volcano", "v-1", contrast_id="U_vs_C")       # 直近の差次的結果の対比
    assert (figures / "v-1_volcano.png").is_file() and rid in msg
    msg = reports.save_figure("volcano", "v-2", result_id=f"{rid}#T_vs_C", contrast_id="T_vs_C")
    assert (figures / "v-2_volcano.png").is_file()
    # 名指しの対比と contrast_id が食い違えば DATASET_BAD_REQUEST
    err = json.loads(reports.save_figure("volcano", "v-3", result_id=f"{rid}#T_vs_C",
                                         contrast_id="U_vs_C"))["error"]
    assert err["code"] == "DATASET_BAD_REQUEST"
    assert not (figures / "v-3_volcano.png").exists()
    # 対比を選ばなければ needs_input（呼び直すのは save_figure）
    _asks_for_contrast(reports.save_figure("volcano", "v-4", result_id=rid), "save_figure",
                       ["T_vs_C", "U_vs_C"])
    assert not (figures / "v-4_volcano.png").exists()
    # contrast_id は volcano だけに効く。ほかの kind では断らずに使わない（loadings はまだ描いていないので missing_state）
    err = json.loads(reports.save_figure("loadings", "v-5", contrast_id="T_vs_C"))["error"]
    assert err["code"] == "missing_state" and err["required_tools"] == ["dataset_plot"]


def test_volcano_caption_reads_test_vs_reference(ds):
    _compare(reference_group="C", test_group="T")
    body = plot_tools.dataset_plot(kind="volcano")
    assert body[0].startswith("Volcano: T (n=8) vs C (n=8)")
    payload = json.loads(plot_tools.dataset_plot(kind="volcano", output="payload"))
    assert payload["title"] == "Volcano (T vs C)"
