"""図と差次的出力を「どの結果から出たか」に結び付ける（spec §8）。

旧実装は、ARF と mzTab-M の両方に結果があると**黙って ARF を採っていた**。
ARF が古いデータのものでも、新しく読み込んだ mzTab の解析の図として保存される。
戻り値に `source=arf` とは書いてあるが、それは「選んだ結果」の報告であって
「選んでよいか」の確認ではない。どちらを使うかは数字そのものを変えるので、
曖昧なら描かずに止める。

もう一つの穴は、古い結果に新しい前処理条件のラベルが載ることだった。前処理を
やり直した後でも `last_differential` は残り、TSV のメタ行には**現在の**レシピが
書かれていた。出力は結果自身の来歴（provenance）から書く。
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metabolomix.core.atomic_io import DomainError
from metabolomix.plots.result_output import save_result_figure, select_result


def _candidate(source, result_id, dataset_id="d1", valid=True, result=None):
    return {"source": source, "result_id": result_id, "dataset_id": dataset_id,
            "valid": valid, "result": result or {}}


# ---------- select_result ----------

def test_auto_does_not_prefer_old_arf():
    candidates = [
        _candidate("arf", "a", "old"),
        _candidate("mztab", "m", "new"),
    ]
    with pytest.raises(DomainError, match="AMBIGUOUS_RESULT_SOURCE"):
        select_result(candidates)
    assert select_result(candidates, source="mztab", result_id="m")["dataset_id"] == "new"


def test_a_single_candidate_needs_no_choice():
    """候補が 1 件のときの従来どおりの呼び出しは、そのまま通す。"""
    chosen = select_result([_candidate("arf", "a")])
    assert chosen["result_id"] == "a"


def test_invalid_candidates_are_not_selectable():
    """前処理をやり直して古くなった結果は、指名しても採らない。"""
    with pytest.raises(DomainError, match="ANALYSIS_RESULT_NOT_FOUND"):
        select_result([_candidate("mztab", "m", valid=False)])


def test_no_candidates_is_not_found():
    with pytest.raises(DomainError) as exc:
        select_result([])
    assert exc.value.code == "ANALYSIS_RESULT_NOT_FOUND"


def test_source_filter_narrows_to_one():
    candidates = [_candidate("arf", "a"), _candidate("mztab", "m")]
    assert select_result(candidates, source="arf")["result_id"] == "a"


def test_two_results_from_the_same_source_still_need_a_result_id():
    candidates = [_candidate("mztab", "m1"), _candidate("mztab", "m2")]
    with pytest.raises(DomainError, match="AMBIGUOUS_RESULT_SOURCE"):
        select_result(candidates, source="mztab")
    assert select_result(candidates, result_id="m2")["result_id"] == "m2"


def test_an_unknown_result_id_is_not_found():
    """別データセットの ID を指定しても、無いものは無い。"""
    with pytest.raises(DomainError, match="ANALYSIS_RESULT_NOT_FOUND"):
        select_result([_candidate("mztab", "m")], result_id="somewhere_else")


# ---------- 図 ----------

_PCA_PLOT = {"title": "PCA", "x_label": "PC1 (31%)", "y_label": "PC2 (22%)",
             "points": [{"x": 1.0, "y": 2.0, "label": "s1"},
                        {"x": -1.0, "y": 0.5, "label": "s2"}]}

_VOLCANO = {
    "kind": "two_group", "a": "ctrl", "b": "treat", "n_a": 3, "n_b": 3,
    "q_threshold": 0.05, "log2fc_threshold": 1.0, "log_transform": True,
    "volcano": [{"feature": "f1", "log2fc": 1.5, "neg_log10_p": 2.0, "sig": "up"},
                {"feature": "f2", "log2fc": -1.8, "neg_log10_p": 2.2, "sig": "down"}],
}


def test_save_result_figure_writes_a_pca_png(tmp_path):
    out = save_result_figure(None, _PCA_PLOT, tmp_path / "a" / "pca.png", kind="pca")
    assert out.is_file()
    assert out.read_bytes()[:4] == b"\x89PNG"


def test_save_result_figure_accepts_the_dataset_state_pca_shape(tmp_path):
    """DatasetState 経路の PCA 結果（`points` ではなく `scores`）でも点を描く。

    ここを揃えないと、pipeline が必須出力として保存する `results/pca.png` が
    枠だけの空の散布図になる（ファイルは生成され hash も通るので気付けない）。
    """
    import matplotlib.image as mpimg

    dataset_pca = {
        "scores": [{"name": "s1", "role": "sample", "PC1": 1.0, "PC2": 2.0},
                   {"name": "s2", "role": "sample", "PC1": -1.0, "PC2": 0.5}],
        "explained_variance_ratio": [0.31, 0.22],
    }

    out = save_result_figure(None, dataset_pca, tmp_path / "scores.png", kind="pca")

    image = mpimg.imread(str(out))
    # 点の色（既定 C0＝青系）の画素があること。軸枠・文字は黒/灰なので0にはならない。
    assert int(((image[..., 2] - image[..., 0]) > 0.1).sum()) > 100


def test_save_result_figure_refuses_a_pca_result_with_no_points(tmp_path):
    """描ける点が無いなら、空の散布図を保存せずに止める。"""
    out = tmp_path / "empty.png"
    with pytest.raises(DomainError) as exc:
        save_result_figure(None, {"explained_variance_ratio": [0.3]}, out, kind="pca")
    assert exc.value.code == "PCA_FIGURE_EMPTY"
    assert not out.exists()


def test_save_result_figure_writes_a_volcano_png(tmp_path):
    out = save_result_figure(None, _VOLCANO, tmp_path / "v.png", kind="volcano")
    assert out.is_file()
    assert out.read_bytes()[:4] == b"\x89PNG"


def test_save_result_figure_rejects_an_unknown_kind(tmp_path):
    with pytest.raises(DomainError) as exc:
        save_result_figure(None, _PCA_PLOT, tmp_path / "x.png", kind="heatmap")
    assert exc.value.code == "UNSUPPORTED_FIGURE_KIND"


def test_a_failed_render_does_not_leave_a_half_written_png(tmp_path):
    """壊れた PNG が「保存できた図」として残ると、レポートに貼られる。"""
    out = tmp_path / "keep.png"
    out.write_bytes(b"previous")
    broken = {"title": "PCA", "points": [{"x": "not-a-number", "y": 0.0}]}
    with pytest.raises(Exception):
        save_result_figure(None, broken, out, kind="pca")
    assert out.read_bytes() == b"previous"


def test_exploratory_datasets_are_marked_on_the_figure(tmp_path):
    """中断された解析から描いた図は、図の上でそう見えなければ意味がない。

    レポートに貼られた PNG は、貼った後は単体で読まれる。「この図は探索専用の
    データから描いた」という但し書きは、戻り値ではなく図の中に要る。
    """
    from metabolomix.mztab.dataset_state import DatasetState
    from metabolomix.plots.result_output import figure_annotations
    ds = DatasetState()
    ds.exploratory_only = True
    ds.source_verification = "legacy_unverified"

    out = save_result_figure(ds, _PCA_PLOT, tmp_path / "e.png", kind="pca")

    assert out.is_file()
    # 文面は ASCII。matplotlib の既定フォントは日本語字形を持たず、豆腐が並ぶ。
    note = figure_annotations(ds)
    assert any("exploratory only" in a for a in note)
    assert all(a.isascii() for a in note)


def test_a_verified_dataset_gets_no_disclaimer():
    """常に但し書きを付けると、本当に注意が要る図が埋もれる。"""
    from metabolomix.mztab.dataset_state import DatasetState
    from metabolomix.plots.result_output import figure_annotations
    ds = DatasetState()
    ds.source_verification = "verified"
    assert figure_annotations(ds) == []


def test_unadjusted_confounded_comparisons_get_an_ascii_disclaimer_on_the_figure(tmp_path):
    """spec §7.4(R16): allow_confounded=trueで継続した比較は、図の中にもその旨が
    要る。matplotlibの既定フォントは日本語字形を持たないため、他の但し書きと
    同じくASCIIで書く。"""
    from metabolomix.plots.result_output import figure_annotations
    result = dict(_VOLCANO)
    result["provenance"] = {"comparison": {"unadjusted_confounded": True}}

    out = save_result_figure(None, result, tmp_path / "unadjusted.png", kind="volcano")

    assert out.is_file()
    note = figure_annotations(None, result)
    assert any("UNADJUSTED" in a for a in note)
    assert all(a.isascii() for a in note)


def test_adjusted_comparisons_get_no_confounding_disclaimer():
    from metabolomix.plots.result_output import figure_annotations
    result = dict(_VOLCANO)
    result["provenance"] = {"comparison": {"unadjusted_confounded": False}}
    assert figure_annotations(None, result) == []


# ---------- TSV ----------

def _prepared_dataset():
    from metabolomix.analysis.dataset_service import compare_dataset, preprocess_dataset
    from tests.pipeline_fixtures import make_dataset
    ds = make_dataset()
    preprocess_dataset(ds, {"normalize": "none", "impute": "half_min"})
    result = compare_dataset(ds, ds.sample_names[:4], ds.sample_names[4:])
    return ds, result


def test_export_writes_the_contract_columns(tmp_path):
    from metabolomix.analysis import export_contract
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()

    info = export_dataset_result(ds, result, tmp_path / "out.tsv")

    lines = (tmp_path / "out.tsv").read_text(encoding="utf-8").splitlines()
    header = [l for l in lines if not l.startswith("#")][0]
    assert header.split("\t") == export_contract.EXPORT_COLUMNS
    assert info["n_with_inchikey"] == 6
    assert info["contract_version"] == export_contract.CONTRACT_VERSION


def test_export_records_the_preprocessing_that_the_result_actually_used(tmp_path):
    """メタ行の前処理は「今のレシピ」ではなく「その結果が使ったレシピ」。"""
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()

    meta = [l for l in
            _exported_lines(export_dataset_result, ds, result, tmp_path / "o.tsv")
            if l.startswith("#")]

    assert any("'normalize': 'none'" in l for l in meta if "preprocess" in l)
    assert any(f"result_id = {result['provenance']['result_id']}" in l for l in meta)


def test_export_records_unadjusted_confounded_status_in_the_meta(tmp_path):
    """spec §7.4(R16): allow_confounded=trueで継続した比較は、TSV付随メタにもその旨が
    要る。`_meta_lines`は既にprovenanceを読んでいるので確認する（見た目は近いが、
    未調整フラグまで読んでいるかは別に確認しないと分からない）。"""
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    result["provenance"]["comparison"] = {"unadjusted_confounded": True}
    result["provenance"]["warnings"] = [
        "allow_confounded=true が明示されたため、群とバッチの完全交絡を未調整の"
        "まま解析を継続しました（図・TSV・レポートにこの旨を残すこと）。"]

    meta = [l for l in
            _exported_lines(export_dataset_result, ds, result, tmp_path / "u.tsv")
            if l.startswith("#")]

    assert any("unadjusted_confounded = true" in l for l in meta)
    assert any("未調整" in l for l in meta)


def test_export_omits_the_unadjusted_confounded_line_when_not_confounded(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()

    meta = [l for l in
            _exported_lines(export_dataset_result, ds, result, tmp_path / "n.tsv")
            if l.startswith("#")]

    assert not any("unadjusted_confounded" in l for l in meta)


def test_export_refuses_a_result_from_a_superseded_preprocessing(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    from metabolomix.analysis.dataset_service import preprocess_dataset
    ds, result = _prepared_dataset()
    preprocess_dataset(ds, {"normalize": "tic", "impute": "half_min"})

    with pytest.raises(DomainError) as exc:
        export_dataset_result(ds, result, tmp_path / "stale.tsv")

    assert exc.value.code == "STALE_ANALYSIS_RESULT"
    assert not (tmp_path / "stale.tsv").exists()


def test_export_refuses_an_exploratory_dataset(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    ds.exploratory_only = True

    with pytest.raises(DomainError) as exc:
        export_dataset_result(ds, result, tmp_path / "x.tsv")

    assert exc.value.code == "EXPLORATORY_ONLY_DATASET"
    assert not (tmp_path / "x.tsv").exists()


def test_export_refuses_when_no_feature_has_an_inchikey(tmp_path):
    """背景集合が空の TSV は、下流のパスウェイ解析にとって無意味。"""
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    for meta in ds.feature_metadata.values():
        meta["inchikey"] = ""

    with pytest.raises(DomainError) as exc:
        export_dataset_result(ds, result, tmp_path / "empty.tsv")

    assert exc.value.code == "NO_ANNOTATED_FEATURES"
    assert not (tmp_path / "empty.tsv").exists()


def test_a_failed_export_does_not_replace_the_previous_file(tmp_path):
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    out = tmp_path / "prev.tsv"
    export_dataset_result(ds, result, out)
    before = out.read_text(encoding="utf-8")

    ds.exploratory_only = True
    with pytest.raises(DomainError):
        export_dataset_result(ds, result, out)

    assert out.read_text(encoding="utf-8") == before


def _exported_lines(func, ds, result, path: Path) -> list[str]:
    func(ds, result, path)
    return path.read_text(encoding="utf-8").splitlines()


# ---------- SML 由来の同定（spec 2026-09-17-mztab-sml-annotation-design §7.3）----------

def _export_rows(tmp_path, ds, result, name="sml.tsv") -> list[dict]:
    """エクスポート TSV の本体を dict の一覧で返す（`#` メタ行は落とす）。"""
    from metabolomix.analysis.dataset_export import export_dataset_result
    path = tmp_path / name
    export_dataset_result(ds, result, path)
    body = [l for l in path.read_text(encoding="utf-8").splitlines()
            if l and not l.startswith("#")]
    header = body[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in body[1:]]


def test_an_evidence_backed_row_says_sme(tmp_path):
    """既存経路の name_source を "mztab_smf" から "mztab_sme" へ訂正した。

    名前の実体は SMF 行ではなく SME 行なので、従来の値は事実と違っていた。
    """
    ds, result = _prepared_dataset()

    rows = _export_rows(tmp_path, ds, result, "sme.tsv")

    assert {r["name_source"] for r in rows} == {"mztab_sme"}


def test_a_row_sourced_from_an_sml_annotation_says_so(tmp_path):
    """MS1 注釈由来の行は name_source で見分けられる。

    下流（massbank-context）は inchikey を結合キーにするので、証拠水準を
    伝える列が無いと MS/MS 裏付けと区別できなくなる。
    """
    ds, result = _prepared_dataset()
    # InChIKey を持つ特徴を1つ選び、証拠を剥がして注釈だけにする。
    fid = next(f for f in ds.feature_ids
               if (ds.feature_metadata.get(f) or {}).get("inchikey"))
    ds.feature_metadata[fid]["inchikey"] = None
    ds.feature_metadata[fid]["name"] = None
    ds.feature_annotations = {
        fid: {"sml_id": "1", "ambiguous": False, "name": "GABA",
              "database_identifier": "TextDB:GABA", "chemical_formula": None,
              "smiles": "NCCCC(=O)O", "adduct": "[M+H]1+", "reliability": None,
              "confidence_measure": "MS-DIAL algorithm matching score",
              "confidence_value": 0.99,
              "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
              "inchikey_source": "smiles_derived"},
    }

    rows = _export_rows(tmp_path, ds, result)

    row = next(r for r in rows if r["spot_id"] == fid)
    assert row["name"] == "GABA"
    assert row["name_source"] == "mztab_sml"
    assert row["inchikey"] == "BTCSSZJGUNDROE-UHFFFAOYSA-N"
    assert row["inchikey_source"] == "smiles_derived"
    assert row["msi_level"] == ""


def test_a_name_without_an_inchikey_is_still_dropped(tmp_path):
    """InChIKey ゲートは変えない。下流が InChIKey で結合するため。"""
    from metabolomix.analysis.dataset_export import export_dataset_result
    ds, result = _prepared_dataset()
    for fid in ds.feature_ids:
        ds.feature_metadata.setdefault(fid, {})["inchikey"] = None
    # 名前だけの注釈（InChIKey なし）。これでは救われない。
    ds.feature_annotations = {
        fid: {"sml_id": "1", "ambiguous": False, "name": "GABA",
              "database_identifier": "TextDB:GABA", "chemical_formula": None,
              "smiles": None, "adduct": "[M+H]1+", "reliability": None,
              "confidence_measure": "score", "confidence_value": 0.99,
              "inchikey": None, "inchikey_source": "none"}
        for fid in ds.feature_ids}

    with pytest.raises(DomainError, match="NO_ANNOTATED_FEATURES"):
        export_dataset_result(ds, result, tmp_path / "dropped.tsv")


def test_an_ambiguous_annotation_does_not_reach_the_export(tmp_path):
    ds, result = _prepared_dataset()
    fid = next(f for f in ds.feature_ids
               if (ds.feature_metadata.get(f) or {}).get("inchikey"))
    ds.feature_metadata[fid]["inchikey"] = None
    ds.feature_annotations = {
        fid: {"ambiguous": True, "sml_ids": ["1", "2"], "name": None}}

    rows = _export_rows(tmp_path, ds, result, "ambiguous.tsv")

    assert not any(r["spot_id"] == fid for r in rows)
