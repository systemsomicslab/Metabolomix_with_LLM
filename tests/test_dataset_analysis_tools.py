import json

import numpy as np
import pytest

from metabolomix.core import session_state
from metabolomix.core.mcp_errors import MISSING_STATE
from metabolomix.mztab.dataset_state import DatasetState


@pytest.fixture(autouse=True)
def reset_session():
    session_state.session = session_state.AnalysisSession()
    yield
    session_state.session = session_state.AnalysisSession()


def _load_ds(n_features=20, n_samples=8):
    rng = np.random.default_rng(0)
    ds = DatasetState()
    ds.feature_matrix = rng.random((n_features, n_samples)) * 1000.0
    ds.sample_names = ([f"ctrl_{i}" for i in range(n_samples // 2)]
                       + [f"treat_{i}" for i in range(n_samples - n_samples // 2)])
    ds.feature_ids = [f"f{i}" for i in range(n_features)]
    ds.feature_metadata = {
        f"f{i}": {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                  "inchikey": f"AAAAAAAAAAAAAA-BBBBBBBBFB-{i%10}",
                  "inchikey_source": "database_identifier"}
        for i in range(n_features)
    }
    ds.validation_result = {"ok": True, "errors": [], "warnings": []}
    session_state.session.dataset = ds
    return ds


def _groups(ds):
    return ([n for n in ds.pp_sample_names if n.startswith("ctrl")],
            [n for n in ds.pp_sample_names if n.startswith("treat")])


# ---------- dataset_preprocess ----------

def test_dataset_preprocess_without_dataset_returns_missing_state():
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    parsed = json.loads(dataset_preprocess())
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_load"]


def test_dataset_preprocess_success_sets_pp_matrix():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    parsed = json.loads(dataset_preprocess())
    assert parsed["status"] == "success"
    assert parsed["n_samples"] == 8
    assert parsed["n_features"] == 20
    ds = session_state.session.dataset
    assert ds.pp_matrix is not None
    assert ds.preprocessing_recipe["normalize"] == "none"


def test_dataset_preprocess_bad_recipe_is_not_missing_state():
    """引数エラーは missing_state ではない（リプレイしても直らない）。"""
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    parsed = json.loads(dataset_preprocess(normalize="not_a_method"))
    assert parsed["error"]["code"] != MISSING_STATE
    assert "not_a_method" in json.dumps(parsed, ensure_ascii=False)


# ---------- dataset_pca ----------

def test_dataset_pca_requires_preprocess():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_pca
    parsed = json.loads(dataset_pca())
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_preprocess"]


def test_dataset_pca_success_omits_loadings_from_payload():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_pca, dataset_preprocess
    dataset_preprocess()
    parsed = json.loads(dataset_pca(n_components=2))
    assert len(parsed["explained_variance_ratio"]) == 2
    assert len(parsed["scores"]) == 8
    assert "loadings" not in parsed          # 全量は戻り値に載せない
    assert session_state.session.dataset.last_pca["loadings"]  # セッションには保持


# ---------- dataset_differential ----------

def test_dataset_differential_requires_preprocess():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_differential
    parsed = json.loads(dataset_differential(group_a=["ctrl_0"], group_b=["treat_0"]))
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_preprocess"]


def test_dataset_differential_success_returns_summary_only():
    ds = _load_ds()
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_preprocess,
    )
    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    parsed = json.loads(dataset_differential(group_a=a, group_b=b))
    assert parsed["status"] == "success"
    assert "n_tested" in parsed["summary"]
    assert len(parsed["summary"]["top"]) <= 15
    assert "results" not in parsed           # 全量は戻り値に載せない
    assert "volcano" not in parsed
    stored = session_state.session.dataset.last_differential
    assert len(stored["results"]) == 20      # セッションには全量


def test_dataset_differential_small_group_is_bad_request():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_preprocess,
    )
    dataset_preprocess()
    parsed = json.loads(dataset_differential(group_a=["ctrl_0"], group_b=["treat_0"]))
    assert parsed["error"]["code"] != MISSING_STATE
    assert "群サイズ" in parsed["error"]["message"]


def test_dataset_differential_does_not_touch_arf_slot():
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_preprocess,
    )
    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)
    assert session_state.session.arf.feature_matrix is None
    assert getattr(session_state.session.arf, "last_differential", None) is None


# ---------- dataset_export_differential ----------

def test_dataset_export_requires_differential(tmp_path):
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import dataset_export_differential
    parsed = json.loads(dataset_export_differential(str(tmp_path / "out.tsv")))
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_differential"]


def test_dataset_export_writes_contract_format(tmp_path):
    from metabolomix.analysis.export_contract import CONTRACT_VERSION, EXPORT_COLUMNS
    _load_ds()
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )
    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)

    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["status"] == "success"
    assert parsed["contract_version"] == CONTRACT_VERSION

    lines = out.read_text(encoding="utf-8").splitlines()
    meta = [l for l in lines if l.startswith("#")]
    assert f"# contract_version = {CONTRACT_VERSION}" in meta
    assert any("id_space = mztab_smf_id" in l for l in meta)
    header = next(l for l in lines if not l.startswith("#"))
    assert header.split("\t") == EXPORT_COLUMNS
    body = [l for l in lines if not l.startswith("#")][1:]
    assert len(body) == 20
    # p_value / q_value が空欄でないこと（初版の欠陥の回帰テスト）
    cols = body[0].split("\t")
    assert cols[EXPORT_COLUMNS.index("p_value")] != ""
    assert cols[EXPORT_COLUMNS.index("q_value")] != ""
    assert cols[EXPORT_COLUMNS.index("inchikey")] != ""


def test_dataset_export_applies_curation_flags_via_sibling_arf2(tmp_path):
    """Task 9: `dataset_export_differential` がツール層で兄弟 `.arf2` のフラグを
    見つけ出し、`export_dataset_result` へ渡すところまで(結合含む)を確認する。

    feature_qc["source"] == "arf" (evidence.py が SMF_ID == MasterAlignmentID の
    対応を数値で検証できたときだけ立つ) が無いと mztab_smf_id と MasterAlignmentID
    を同じ空間として扱ってよいと言えないため、ここで明示的に立てる。
    """
    from metabolomix.curation import flags as curation_flags
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )

    # feature_ids は mzTab の SMF_ID そのもの(数字の連番文字列)を模す。
    # curation フラグの spot_id は int の MasterAlignmentID なので、両者の
    # 空間が一致するのは `str(int) == SMF_ID 文字列` のときだけ。
    ds = _load_ds()
    ds.feature_ids = [str(i) for i in range(20)]
    ds.feature_metadata = {
        str(i): {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                 "inchikey": f"AAAAAAAAAAAAAA-BBBBBBBBFB-{i % 10}",
                 "inchikey_source": "database_identifier"}
        for i in range(20)
    }
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    ds.source_files = {str(mztab): "sha"}
    ds.feature_qc = {"source": "arf"}

    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
        [{"spot_id": 0, "flag": "wrong"}],
        alignment=curation_flags.alignment_key(arf2), review_id="r", source="user")

    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)

    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["status"] == "success"
    assert parsed["curation"] == {"state": "applied", "wrong_excluded": 1,
                                  "suspect": 0, "orphaned": 0}
    lines = out.read_text(encoding="utf-8").splitlines()
    body = [l for l in lines if l and not l.startswith("#")][1:]
    assert not any(l.split("\t")[0] == "0" for l in body)
    meta = [l for l in lines if l.startswith("# curation")]
    assert len(meta) == 1
    assert meta[0].startswith("# curation = applied\t")

def test_dataset_export_applies_assign_and_redundant(tmp_path):
    from metabolomix.curation import flags as curation_flags
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )

    ds = _load_ds()
    ds.feature_ids = [str(i) for i in range(20)]
    ds.feature_metadata = {
        str(i): {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                 "inchikey": f"AAAAAAAAAAAAAA-BBBBBBBBFB-{i % 10}",
                 "inchikey_source": "database_identifier"}
        for i in range(20)
    }
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    ds.source_files = {str(mztab): "sha"}
    ds.feature_qc = {"source": "arf"}
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
        [{"spot_id": 3, "flag": "assign", "name": "PE 36:2", "level": "sum", "ontology": "PE",
          "inchikey": "PEPEPEPEPEPEPE-XXXXXXXXXX-N"},
         {"spot_id": 4, "flag": "redundant", "of": 3, "relation": "isotope_M+1"}],
        alignment=curation_flags.alignment_key(arf2), review_id="cs-x", source="user")

    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)
    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["status"] == "success"
    assert parsed["curation"]["assigned"] == 1 and parsed["curation"]["redundant_excluded"] == 1
    lines = out.read_text(encoding="utf-8").splitlines()
    body = [l.split("\t") for l in lines if l and not l.startswith("#")]
    header, rows = body[0], [dict(zip(body[0], r)) for r in body[1:]]
    by_id = {r["spot_id"]: r for r in rows}
    assert "4" not in by_id
    assert (by_id["3"]["name"], by_id["3"]["name_source"], by_id["3"]["inchikey"]) ==         ("PE 36:2", "curation", "PEPEPEPEPEPEPE-XXXXXXXXXX-N")


def test_dataset_export_drops_an_assign_without_inchikey(tmp_path):
    """R10: 記録した InChIKey が空の assign は、同定なしの行として落ちる（ARF 経路と同じ）。"""
    from metabolomix.curation import flags as curation_flags
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )

    ds = _load_ds()
    ds.feature_ids = [str(i) for i in range(20)]
    ds.feature_metadata = {
        str(i): {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                 "inchikey": f"AAAAAAAAAAAAAA-BBBBBBBBFB-{i % 10}",
                 "inchikey_source": "database_identifier"}
        for i in range(20)
    }
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    ds.source_files = {str(mztab): "sha"}
    ds.feature_qc = {"source": "arf"}
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
        [{"spot_id": 3, "flag": "assign", "name": "PE 36:2", "level": "sum", "ontology": "PE",
          "inchikey": ""}],
        alignment=curation_flags.alignment_key(arf2), review_id="cs-x", source="user")

    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)
    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["status"] == "success"
    lines = out.read_text(encoding="utf-8").splitlines()
    body = [l.split("	") for l in lines if l and not l.startswith("#")]
    by_id = {r["spot_id"]: r for r in (dict(zip(body[0], r)) for r in body[1:])}
    assert "3" not in by_id and "4" in by_id                 # 3 は落ち、他のスポットは影響を受けない


def test_dataset_export_reports_a_malformed_flags_file(tmp_path):
    """I6: 壊れた flags.jsonl は wrong を黙って落とさず、ファイルと行を名指しで止める。"""
    from metabolomix.curation import flags as curation_flags
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )

    ds = _load_ds()
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    ds.source_files = {str(mztab): "sha"}
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    store = curation_flags.FlagStore(curation_flags.curation_dir(arf2))
    store.path.parent.mkdir(parents=True)
    store.path.write_text('{"spot_id": 0, "flag": "wr\n', encoding="utf-8")

    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)

    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["error"]["code"] == "CURATION_FLAGS_INVALID"
    assert parsed["error"]["details"] == {"flags_file": str(store.path), "line": 1}
    assert not out.exists()


def test_dataset_export_refuses_without_inchikey(tmp_path):
    """InChIKey が 0 件なら書かずに拒否する（arf_export_differential と同じ理由）。"""
    ds = _load_ds()
    ds.feature_metadata = {fid: {"name": None, "mz": None, "rt": None,
                                 "inchikey": None, "inchikey_source": "none"}
                           for fid in ds.feature_ids}
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )
    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)
    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    # 機械可読なコードで返す（"status": "error" だけでは分岐できない）。
    assert parsed["error"]["code"] == "NO_ANNOTATED_FEATURES"
    assert parsed["error"]["details"]["n_with_inchikey"] == 0
    assert not out.exists()


# --- min_detection_rate（evidence sidecar 由来）---

def test_dataset_preprocess_accepts_min_detection_rate():
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    ds = _load_ds(n_features=4, n_samples=4)
    ds.detected_mask = np.array([
        [True, True, True, True],
        [True, True, False, False],
        [True, False, False, False],
        [False, False, False, False],
    ])
    ds.feature_qc = {"source": "arf"}

    payload = json.loads(dataset_preprocess(impute="none", min_detection_rate=0.5))
    assert payload["status"] == "success"
    assert payload["n_features"] == 2
    assert payload["detection_filter"]["features_removed"] == 2
    assert ds.preprocessing_recipe["min_detection_rate"] == 0.5


def test_dataset_preprocess_min_detection_rate_without_state_is_bad_request():
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    _load_ds(n_features=4, n_samples=4)
    payload = json.loads(dataset_preprocess(min_detection_rate=0.5))
    # missing_state ではない——別ツールを呼んでも直らない引数由来のエラー
    assert payload["error"]["code"] != MISSING_STATE
    assert "min_detection_rate" in json.dumps(payload, ensure_ascii=False)


def test_dataset_preprocess_reports_detection_when_available():
    from metabolomix.tools.dataset_analysis_tools import dataset_preprocess
    ds = _load_ds(n_features=2, n_samples=4)
    ds.detected_mask = np.array([[True, True, True, True],
                                 [True, False, False, False]])
    ds.feature_qc = {"source": "arf"}
    payload = json.loads(dataset_preprocess(impute="none"))
    assert payload["detection"]["gap_filled_rate"] == 0.375


# ---------- 結果 ID と無効化（Task 6） ----------
# 前処理をやり直したのに古い差次的結果が残っていると、TSV と図が別々の前処理から
# 出た数字を並べる。ツール層でもその境界が効いていることを確かめる。

def test_tools_report_the_result_id_of_each_computation():
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_pca, dataset_preprocess)
    ds = _load_ds()
    pp = json.loads(dataset_preprocess())
    a, b = _groups(ds)
    pca = json.loads(dataset_pca(n_components=2))
    diff = json.loads(dataset_differential(a, b))
    assert pp["result_id"].startswith("res_")
    assert len({pp["result_id"], pca["result_id"], diff["result_id"]}) == 3


def test_repreprocessing_makes_the_previous_comparison_unavailable(tmp_path):
    """古い結果を「直近の結果」として書き出させない。"""
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess)
    ds = _load_ds()
    dataset_preprocess()
    a, b = _groups(ds)
    dataset_differential(a, b)
    dataset_preprocess(normalize="tic")

    parsed = json.loads(dataset_export_differential(str(tmp_path / "out.tsv")))

    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_differential"]
    assert not (tmp_path / "out.tsv").exists()


# ---------- dataset_set_sample_metadata（Task 9のparse/resolve/applyをMCPへ接続） ----------

def _load_ds_with_raws(tmp_path, n_samples=8):
    """assay_sources（実行時rawの絶対パス）を持つDatasetStateをセッションへ設定する。"""
    ds = _load_ds(n_features=4, n_samples=n_samples)
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    ds.assay_sources = {}
    ds.sample_assay_ids = [f"assay[{i + 1}]" for i in range(n_samples)]
    for i, name in enumerate(ds.sample_names):
        raw_path = raw_dir / f"{name}.wiff"
        raw_path.write_text("x", encoding="utf-8")
        ds.assay_sources[f"abundance_assay[{i + 1}]"] = str(raw_path)
    return ds


def _manifest_text(ds, groups):
    header = ("# schema = sample-manifest.v1\n"
              "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude")
    lines = [header]
    for i, name in enumerate(ds.sample_names):
        lines.append(f"{name}\t{name}.wiff\tsample\t{groups[i]}\tB1\t{i + 1}\t\ttrue")
    return "\n".join(lines) + "\n"


def test_dataset_set_sample_metadata_without_dataset_returns_missing_state():
    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata
    parsed = json.loads(dataset_set_sample_metadata("manifest.tsv"))
    assert parsed["error"]["code"] == MISSING_STATE
    assert parsed["error"]["required_tools"] == ["dataset_load"]


def test_dataset_set_sample_metadata_success_applies_rows(tmp_path):
    ds = _load_ds_with_raws(tmp_path)
    manifest = tmp_path / "sample-manifest.tsv"
    groups = ["ctrl"] * 4 + ["treat"] * 4
    manifest.write_text(_manifest_text(ds, groups), encoding="utf-8")

    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata
    parsed = json.loads(dataset_set_sample_metadata(str(manifest)))

    assert parsed["status"] == "success"
    assert parsed["metadata_revision"] == 1
    assert ds.sample_metadata_rows[0]["group"] == "ctrl"


def test_dataset_set_sample_metadata_invalid_sheet_does_not_touch_session(tmp_path):
    ds = _load_ds_with_raws(tmp_path)
    manifest = tmp_path / "bad-manifest.tsv"
    manifest.write_text(
        "# schema = sample-manifest.v1\n"
        "sample_id\tsource_file\trole\tgroup\tbatch\tinjection_order\tqc_pool\tinclude\n"
        "dup\tctrl_0.wiff\tsample\t\t\t\t\ttrue\n"
        "dup\tctrl_1.wiff\tsample\t\t\t\t\ttrue\n",
        encoding="utf-8",
    )
    before_revision = ds.metadata_revision

    from metabolomix.tools.dataset_analysis_tools import dataset_set_sample_metadata
    parsed = json.loads(dataset_set_sample_metadata(str(manifest)))

    assert parsed["error"]["code"] != MISSING_STATE
    assert ds.metadata_revision == before_revision
    assert ds.sample_metadata_rows is None


def test_export_records_which_result_it_came_from(tmp_path):
    from metabolomix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess)
    ds = _load_ds()
    dataset_preprocess()
    a, b = _groups(ds)
    diff = json.loads(dataset_differential(a, b))
    dataset_export_differential(str(tmp_path / "out.tsv"))

    meta = [l for l in (tmp_path / "out.tsv").read_text(encoding="utf-8").splitlines()
            if l.startswith("#")]
    assert any(f"result_id = {diff['result_id']}" in l for l in meta)
    assert any(l.startswith("# preprocess_id = res_") for l in meta)
