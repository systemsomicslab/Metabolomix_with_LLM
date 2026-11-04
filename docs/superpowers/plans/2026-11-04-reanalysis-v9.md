# 再解析 v9 の道具（根拠バンドル v3・候補の印・検出の例数・質量の窓の中心ほか）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 手順書 v9 と massbank-context の外来の分子の走査が要る ms-data-parser 側の道具（根拠バンドル v3・`labels`・検出の例数と Fisher・`ppm_offset`・小さな修正 6 件）を入れる。

**Architecture:** 数値は `metabolomix/analysis/` の純ロジック（`evidence_export` `analysis_spec` `spec_statistics` `mass_accuracy` `feature_search` `alignment_qc` `claim_evidence`）に足し、MCP 層（`tools/dataset_analysis_tools.py` `tools/reanalysis_tools.py` `pathway/tools.py`）は引数を通して戻り値を小さく保つだけにする。欄を変えたタスクは、同じタスクの中で出力の文書（`docs/output_format/*.md` ＝ MCP リソース `lipidmix://docs/output-format/{topic}`）と `USAGE.md` を直す。

**Tech Stack:** Python 3.14・numpy・scipy（`scipy.stats.fisher_exact`）・pytest。MCP は FastMCP。

**Spec:**
- `docs/superpowers/specs/2026-11-04-evidence-bundle-v3-samples-design.md`（根拠バンドル v3 と `labels`）
- `C:\Users\yuu18\reanalysis-study\docs\specs\2026-11-04-protocol-v9-design.md`（置き場所が mdp の項目: #4・#6 の道具・#10・#11・#12・#14・#15・#17・#20 の道具）

**依存（ここを先に入れる理由）:** 再解析の手順書 v9（reanalysis-study）と massbank-context の外来の分子の走査は、Task 1〜6 と Task 6b（バンドル v3・`labels`・#4・#6）に依存する。massbank-context は Task 2 の見本のバンドル `tests/fixtures/evidence_bundle_v3/` を**写して**自分のテストで読む（設計 §5）。見本を変えたら massbank-context 側の写しも同時に更新する。手順書 v9 の「始める前に」（protocol #24）は `dataset_export_evidence_bundle` の説明に `v3` があることを目印に読むので、Task 1 で説明文に `v3` を入れる。

## 設計との食い違い（実装の前にコードで確かめたこと）

1. **#14 は実装済み**: `metabolomix/tools/handle_tools.py` の `_pipeline_arf2` と `_alignment_files`（コミット fd70f6d、2026-10-28。テスト `tests/test_dataset_handle_tools.py::test_handle_open_finds_pipeline_arf2`）。issues.md の事象は修正の前の記録。Task 10 で確かめるだけにする。
2. **#17 は実装済み**: `feature_search.DEFAULT_ADDUCT_TIERS["negative"][1]` に `[M+CH3COO]-` がある（コミット 31c4b6f、2026-10-30。テスト `tests/test_feature_search.py::test_negative_tiers_include_acetate`）。Task 10 で確かめるだけにする。
3. **「補完の前の値を持っているなら」**: DatasetState は補完の前の値を持たない（`ds.pp_matrix` は `preprocessing.preprocess` の最後の段 `impute` の後）。`impute` は最後の段なので、結果のレシピを `impute="none"` で当て直せば、正規化などは同じで補完だけが無い値になる。そこで、レシピが分かる結果（session の結果・pipeline から引き継いだ結果のどちらも `feature_search._result_recipe` で引ける）は常に `pre_imputation` を書く。レシピが分からず、`ds.pp_matrix` がその結果の親の前処理のときだけ `post_imputation`。どちらでもなければ書かずにエラー（「必ず書く」を、別の前処理の値で満たさない）。
4. **「書き出しは v3 だけにする」**: 設計が挙げるのは `dataset_export_evidence_bundle` だけで、`arf_export_evidence_bundle`（ARF 経路）は同じ書き出し関数 `write_evidence_bundle` を使うが試料ごとの値の取り出しが無い。ARF 経路は v2 のまま（`samples`・`intensities` を渡さなければ v2。設計者（本体のセッション）の判断 2026-11-04）。
5. **#4 の列名**: 手順書は `detected_test` / `detected_reference`。既存の行の列は `n_ref` `mean_ref` と `_ref` の略。手順書の名前に従う（reanalysis-study が同じ設計書から読むため）。
6. **#4 の `detection_driven`**: 手順書は「規則【手順書】」に置くが、道具が機械的に計算して行に出す（`meets_thresholds` と同じく判定ではなく機械の当てはめ。版のあるテスト済みの実装から出すという再解析の方針）。手順書はこの欄を読んで総合評価に使う。
7. **#11 の「検出の例数（受け持ち）」**: `alignment_qc._pair` の `n_owned`（試料ごとの受け持ちの数。検出状態があれば検出、どちらでも検出されない試料は値で補う）を使う。検出状態だけの数ではない。
8. **#12 の「名前を分ける」**: `counts.msms_available`（中身は規則で qualified の数）を消すと reanalysis-study の triage が壊れるので、同じ値の `rule_qualified` を足して `msms_available` は非推奨の別名として v9 では残す（消すのは次の版。設計者（本体のセッション）の判断 2026-11-04）。新しい `msms_spectrum_available` が代表スペクトルの有無で数える。
9. **既存の不具合**: `dataset_export_evidence_bundle` は `DomainError` を `exc.details()` で受けている（`details` は dict なので TypeError）。Task 1 で `exc.details or None` に直す。
10. **#20 の既定の変更**: `pathway_plot` の戻り値を既定で「場所と要約」にすると、対話（Claude Desktop）で図がチャットに出なくなる。`inline_image=True` で従来の画像を返せるようにする（設計者（本体のセッション）の判断 2026-11-04: 既定の変更はこのまま入れ、出力の文書 `docs/output_format/pathway.md` と `USAGE.md` に書く）。

## Global Constraints

- コード・コメント・docstring・文書・コミットメッセージは日本語。
- テストはリポジトリルートから `.venv/Scripts/python.exe -m pytest tests -q`（worktree では main ツリーの `.venv` の Python）。
- コミットは pre-commit が全テストを回す（数分）。自動化した `git commit` はバックグラウンドで実行し、出力をファイルへ落として `FAILED` を引く。実行中にファイルを編集しない。`git stash` を使わない。
- ツールの戻り値は `json_payload()`、全ツール `structured_output=False`、float は丸めて返す。大きな結果はファイルに書き要約だけを返す。
- 研究室の参照ライブラリのパス・ファイル名を追跡対象（テスト・文書・見本のバンドル）に書かない。見本のバンドルにはローカルの絶対パスも入れない。
- 外部契約の名前は旧名のまま（`lipidmix://`・`LIPIDMIX_*`）。この計画は新しい環境変数を足さない。
- バンドル v3: `contract_version` 3、`exporter.version` `evidence-v3`、`comparison.polarity`（`positive` / `negative`）、`intensities_stage` は `"pre_imputation" | "post_imputation"`、`samples.tsv` の列 `sample_id group role injection_order`、`intensities.tsv` の列 `feature_id sample_id value detected`、`detected` は `true` / `false` / 空、検出状態が無い理由は `population.detection_unavailable_reason`。
- `labels` の既定: `feature_column` `feature_id`、`label_column` `level`、`order` `["A","B","C","D"]`。印は並び・件数・除外を変えない。
- #4: 高い側の群の検出が群の例数の半分に満たない（`<`）か、低い側の群の検出が 0 なら `detection_driven`。
- #10: 外れがすべて 0.02 Da 未満（`<`）なら、外れの割合が 1% を超えても合格（`MZ_HARD_LIMIT` 1 Da の棄却はそのまま）。
- #11: 残す側は受け持ち（`n_owned`）の多い側、同数なら小さい ID。
- 既存の `rule_hash` は、新しい引数を省いた呼び出しでは変えない（`ppm_offset`・`labels` を渡さない spec は v8 と同じ hash）。
- 調査・実装の記録は main ツリーの `docs/HISTRY.md`（末尾へ追記のみ）と `docs/task.md`。

## Review Focus

1. **結果のレシピが分からず、`ds.pp_matrix` が別の前処理のもの** — 別の段の値を `post_imputation` と名乗って書かず、`EVIDENCE_INTENSITIES_UNAVAILABLE` で止まる（Task 1 のテストで固定）。
2. **検出状態はあるが、結果の試料名が `ds.sample_names` に無い**（試料軸のずれ） — `detected` を `false` で埋めず空にし、理由を書く（Task 1・Task 4 のテストで固定）。
3. **`labels` の TSV に、データセットに無い特徴 ID・空の水準・`order` に無い水準がある** — 無い ID は数えて無視（`n_unmatched`）、空の水準は無視、未知の水準は `NOVEL_SPEC_INVALID`（黙って最弱にしない）（Task 3）。
4. **`ppm_offset` の帯の境目** — m/z がちょうど `max_mz` の特徴は次の帯（`_tolerance` の `<` と同じ）（Task 5）。
5. **割れの組を `keep` の無い古い alignment-qc JSON から読む** — 落ちずに小さい ID を残す（従来どおり）（Task 8）。

---

### Task 1: 根拠バンドル v3 の書き出し（samples・intensities・adduct・polarity）

**Files:**
- Modify: `metabolomix/analysis/evidence_export.py`（`write_evidence_bundle` `export_dataset_evidence`、新しい `_samples_table` `_intensity_rows`）
- Modify: `metabolomix/tools/dataset_analysis_tools.py`（`dataset_export_evidence_bundle`）
- Modify: `docs/output_format/reanalysis.md`（新しい §10「根拠バンドル v3（evidence-v3）」）、`docs/workflow/dataset_analysis.md`（`dataset_export_evidence_bundle` の呼び出し順）、`USAGE.md`（`dataset_export_evidence_bundle` の行）
- Test: `tests/test_evidence_export_v3.py`（新規）

**Interfaces:**
- Consumes: `feature_search._result_recipe(ds, result) -> dict | None`、`dataset_analysis.run_dataset_preprocess(ds, recipe)`、`dataset_analysis.build_dataset_pp_inputs(ds)`（`roles` と `sample_meta[name]["run_order"]`）、`alignment_qc._feature_adducts(ds) -> dict[str, str | None]`、`mztab.polarity.dataset_polarity(ds)`
- Produces:
  - 定数 `BUNDLE_CONTRACT_VERSION = 3`、`FEATURE_COLUMNS_V3 = [*FEATURE_COLUMNS, 'adduct']`、`SAMPLE_COLUMNS = ['sample_id','group','role','injection_order']`、`INTENSITY_COLUMNS = ['feature_id','sample_id','value','detected']`
  - `write_evidence_bundle(*, comparison, features, annotations, spectrum_links, spectra_files, output_dir, source_analysis_id=None, exclusions=None, provenance=None, preprocessing_report=None, samples=None, intensities=None, intensities_stage=None, detection_unavailable_reason=None) -> dict` — `samples` と `intensities` を両方渡せば v3（片方だけは `ValueError`）、どちらも無ければ v2（ARF 経路）。戻り値に `contract_version` `n_samples` `n_intensity_rows` `intensities_stage`。
  - `export_dataset_evidence(ds, result, output_dir, *, curation=None, polarity=None) -> dict`（常に v3）
  - `_samples_table(ds, result) -> list[dict]`、`_intensity_rows(ds, result) -> tuple[list[dict], str, str | None]`（行・`intensities_stage`・`detection_unavailable_reason`）
  - 新しいエラー: `DomainError("EVIDENCE_INTENSITIES_UNAVAILABLE", ...)`（レシピが分からず親の前処理の行列も無い）、`DomainError("EVIDENCE_INTENSITIES_MISMATCH", ...)`（当て直した特徴の集合が結果の行と違う）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_evidence_export_v3.py` に合成の入力を作る `_bundle_inputs(tmp_path, *, detected: bool)`（`tests/reanalysis_fixtures.write_mztab` → `mztab.loading.load_dataset_state` → `dataset_service.preprocess_dataset(ds, {"normalize": "none"})` → `dataset_service.compare_dataset(ds, ["HC1","HC2","HC3"], ["RA1","RA2","RA3"], group_a_label="HC", group_b_label="RA")`。`detected=True` なら `ds.detected_mask` を全 True にして (特徴 "1", RA1) だけ False、特徴 "2" の HC1 の値を NaN にしてから前処理）を置く。ツールのテスト用に `_handle_with_result(tmp_path) -> tuple[str, str]`（`write_mztab` → `handle_tools.dataset_handle_open(mztab_path=..., library="none")` → `T.dataset_preprocess(dataset_handle=...)` → `T.dataset_differential(group_a=[...HC], group_b=[...RA], dataset_handle=...)` の `result_id`）も置く。

```python
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

def test_samples_group_is_empty_for_qc(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = {r["sample_id"]: r for r in _tsv(tmp_path / "b" / "samples.tsv")}
    assert [rows[s]["group"] for s in ("HC1", "RA1", "QC1")] == ["a", "b", ""]
    assert rows["QC1"]["role"] == "qc" and rows["HC1"]["role"] == "sample"
    assert rows["HC1"]["injection_order"] in ("", "1")      # 分からなければ空、分かれば整数

def test_intensities_long_format_detected_and_pre_imputation_value(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    export_dataset_evidence(ds, result, tmp_path / "b")
    rows = _tsv(tmp_path / "b" / "intensities.tsv")
    assert len(rows) == len(result["results"]) * 6          # 全特徴 × 比較の試料（QC を含まない）
    cell = {(r["feature_id"], r["sample_id"]): r for r in rows}
    assert cell[("1", "RA1")]["detected"] == "false" and cell[("1", "HC1")]["detected"] == "true"
    assert cell[("1", "HC1")]["value"] == "10"               # normalize none・補完の前
    assert cell[("2", "HC1")]["value"] == ""                 # NaN は補完せず空

def test_without_detection_state_detected_is_empty_and_reason_recorded(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=False)
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert {r["detected"] for r in _tsv(tmp_path / "b" / "intensities.tsv")} == {""}
    assert manifest["population"]["detection_unavailable_reason"]   # 空でない文字列

def test_result_samples_outside_mask_axis_leave_detected_empty(tmp_path):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    result = {**result, "samples_a": ["HC1", "HC2", "ZZ"]}  # 試料軸のずれ
    export_dataset_evidence(ds, result, tmp_path / "b")
    manifest = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["population"]["detection_unavailable_reason"] == "sample_axis_mismatch"

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

def test_refuses_when_neither_recipe_nor_parent_matrix(tmp_path, monkeypatch):
    ds, result = _bundle_inputs(tmp_path, detected=True)
    monkeypatch.setattr(feature_search, "_result_recipe", lambda ds, result: None)
    ds.preprocess_id = "other"
    with pytest.raises(DomainError) as info:
        export_dataset_evidence(ds, result, tmp_path / "b")
    assert info.value.code == "EVIDENCE_INTENSITIES_UNAVAILABLE"
    assert not (tmp_path / "b").exists()

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
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_evidence_export_v3.py -q`
Expected: FAIL（`contract_version` が 2、`samples.tsv` が無い、`polarity` 引数が無い）

- [ ] **Step 3: `evidence_export.py` に v3 を実装する**

- `write_evidence_bundle` は v3 のとき features を `FEATURE_COLUMNS_V3` で、`samples.tsv`・`intensities.tsv` を足して書く（`intensities` は 28 万行になりうるので反復で 1 行ずつ書く）。manifest に `intensities_stage`、`comparison` は呼び出し側が `polarity` を入れたもの、`population.detection_unavailable_reason`。`value` は `export_contract.format_number(v, ".10g")`（NaN・None は空）、`detected` は `true` / `false` / 空。
- `_samples_table`: 試料は `ds.sample_names` の全部。`group` は `result["samples_a"]` なら `a`、`samples_b` なら `b`、ほかは空。`role` と `injection_order` は `build_dataset_pp_inputs(ds)` の `roles` と `sample_meta[name]["run_order"]`（整数でなければ空）。
- `_intensity_rows`: 試料は `samples_a + samples_b`、特徴は `result["results"]` の順。値はレシピが分かれば `run_dataset_preprocess(ds, {**recipe, "impute": "none"})`（`pre_imputation`）、分からず `ds.preprocess_id` が `result["provenance"]["parent_ids"]` にあれば `ds.pp_matrix`（`post_imputation`）、どちらでもなければ `EVIDENCE_INTENSITIES_UNAVAILABLE`。当て直した行列の特徴の集合が結果の行と違えば `EVIDENCE_INTENSITIES_MISMATCH`。検出は `ds.detected_mask`（特徴 × `ds.sample_names`）から引き、マスクが無ければ理由は `(ds.feature_qc or {}).get("reason") or "not_imported"`、形が合わなければ `"mask_shape_mismatch"`、結果の試料が `ds.sample_names` に無ければ `"sample_axis_mismatch"`（どれも `detected` は空）。
- `feature_row` に `adduct`（`alignment_qc._feature_adducts(ds)` の値。無ければ空）を足すのは `export_dataset_evidence` の側で行う（ARF 経路の行は変えない）。
- `export_dataset_evidence(..., polarity=None)`: `comparison["polarity"] = polarity`。

- [ ] **Step 4: ツールを直す**

`dataset_export_evidence_bundle`: 極性はハンドルなら `holder.polarity`、対話のセッションなら `_dataset_polarity(ds)[0]` を渡す。説明文を「指定した 2 群結果を根拠バンドル v3（samples・intensities・adduct・comparison.polarity）へ出力する」に直す。`exc.details()` を `exc.details or None` に直す。

- [ ] **Step 5: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_evidence_export_v3.py tests/test_evidence_export.py tests/test_arf_evidence_export.py tests/test_dataset_handle_tools.py -q`
Expected: PASS（既存の v2 のテスト `test_dataset_export_*` のうち `contract_version == 2` を見るものが dataset 経路なら 3 に直す。ARF 経路は 2 のまま）

- [ ] **Step 6: 文書を直す**

`docs/output_format/reanalysis.md` に §10 を足す: manifest の欄（v2 との差分だけ）・`samples.tsv`・`intensities.tsv` の列と値の意味（`value` は log 変換の前、`pre_imputation` の決め方、`detected` の空は「不明」で gap-fill ではない）・`detection_unavailable_reason` の値・`adduct`・特徴 ID はスポットの ID のまま（まとめた組への読み替えは呼び出し側が `merged_from` で行う）・契約の正は massbank-context の `docs/contracts/evidence-bundle-v3.md`。冒頭の目次にも §10。`USAGE.md` の `dataset_export_evidence_bundle` の行を v3 に。`docs/workflow/dataset_analysis.md` の呼び出し順に `_samples_table` `_intensity_rows` を足す。

Run: `.venv/Scripts/python.exe -m pytest tests/test_workflow_docs.py tests/test_readme_links.py -q`
Expected: PASS

- [ ] **Step 7: コミット**

```bash
git add metabolomix/analysis/evidence_export.py metabolomix/tools/dataset_analysis_tools.py tests/test_evidence_export_v3.py tests/test_evidence_export.py docs/output_format/reanalysis.md docs/workflow/dataset_analysis.md USAGE.md
git commit -m "feat(evidence): 根拠バンドル v3（試料ごとの値と検出・付加体・極性）を書く"
```

---

### Task 2: massbank-context と共有する見本のバンドル

**Files:**
- Create: `tests/evidence_bundle_fixture.py`（見本を作る決定的な関数と、`python -m tests.evidence_bundle_fixture <出力先>` の入口）
- Create: `tests/fixtures/evidence_bundle_v3/`（`manifest.json` `features.tsv` `annotations.tsv` `spectrum_links.tsv` `samples.tsv` `intensities.tsv`）
- Modify: `docs/superpowers/specs/2026-11-04-evidence-bundle-v3-samples-design.md`（§3.3 に規模の実測を追記）
- Test: `tests/test_evidence_bundle_fixture.py`（新規）

**Interfaces:**
- Consumes: Task 1 の `export_dataset_evidence(ds, result, output_dir, *, curation=None, polarity=None)`
- Produces: `build_fixture_bundle(output_dir: Path) -> Path`（manifest のパス）。見本は陽イオン、比較 HC（a）対 RA（b）、QC 2 本（`group` 空）、`detected` に true / false の両方、補完の前の NaN（`value` 空）1 つ、`adduct` のある行と無い行、有意な行と有意でない行、SME 候補 2 件の特徴を含む。`ds.source_files` は `["fixture.mzTab"]` に置き換えてから書く（ローカルの絶対パスを入れない）。

- [ ] **Step 1: 失敗するテストを書く**

```python
FIXTURE = Path(__file__).parent / "fixtures" / "evidence_bundle_v3"
VOLATILE = ("exported_at", "source_analysis_id", "provenance")

def test_fixture_matches_current_exporter(tmp_path):
    build_fixture_bundle(tmp_path / "b")
    for name in ("features.tsv", "annotations.tsv", "spectrum_links.tsv",
                 "samples.tsv", "intensities.tsv"):
        assert (tmp_path / "b" / name).read_bytes() == (FIXTURE / name).read_bytes(), name
    fresh = json.loads((tmp_path / "b" / "manifest.json").read_text(encoding="utf-8"))
    kept = json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8"))
    assert {k: v for k, v in fresh.items() if k not in VOLATILE} == \
           {k: v for k, v in kept.items() if k not in VOLATILE}

def test_fixture_covers_the_cases_massbank_context_reads():
    rows = _tsv(FIXTURE / "intensities.tsv")
    assert {r["detected"] for r in rows} == {"true", "false"}
    assert any(r["value"] == "" for r in rows)
    assert "" in {r["group"] for r in _tsv(FIXTURE / "samples.tsv")}
    manifest = json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["contract_version"] == 3 and manifest["comparison"]["polarity"] == "positive"

def test_fixture_has_no_local_paths():
    text = "".join(p.read_text(encoding="utf-8") for p in FIXTURE.iterdir())
    assert ":\\" not in text and ":/" not in text.replace("file:///", "")
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_evidence_bundle_fixture.py -q`
Expected: FAIL（`tests.evidence_bundle_fixture` が無い）

- [ ] **Step 3: `build_fixture_bundle` を書き、見本を生成する**

入力は Task 1 の `_bundle_inputs` と同じ作り方（`reanalysis_fixtures.write_mztab` を一時ディレクトリに書く）。生成:

Run: `.venv/Scripts/python.exe -m tests.evidence_bundle_fixture tests/fixtures/evidence_bundle_v3`
Expected: `tests/fixtures/evidence_bundle_v3/manifest.json` ほか 5 ファイル

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_evidence_bundle_fixture.py -q`
Expected: PASS

- [ ] **Step 5: 規模を測って設計書に追記する**

`check.py` に、特徴 23,000 × 試料 12（検出状態あり）の合成のデータセットで `export_dataset_evidence` を回し、所要時間と `intensities.tsv` の大きさを出すコードを一時的に書いて実行し、設計書 §3.3 の「規模の目安」の下に「実測（2026-11-xx）: n 秒・n MB」を追記する。終わったら `check.py` の中身を消す。

- [ ] **Step 6: コミット**

```bash
git add tests/evidence_bundle_fixture.py tests/test_evidence_bundle_fixture.py tests/fixtures/evidence_bundle_v3 docs/superpowers/specs/2026-11-04-evidence-bundle-v3-samples-design.md
git commit -m "test(evidence): massbank-context と共有する v3 の見本のバンドルを置く"
```

---

### Task 3: 新しい知見の候補の印（`labels`）

**Files:**
- Modify: `metabolomix/analysis/mass_accuracy.py`（`normalize_novel_spec` `rule_hash`、新しい `read_labels` `attach_labels`）
- Modify: `metabolomix/tools/reanalysis_tools.py`（`dataset_novel_candidates`）
- Modify: `docs/output_format/reanalysis.md` §8、`USAGE.md`（`dataset_novel_candidates` の行）
- Test: `tests/test_novel_labels.py`（新規）、`tests/test_reanalysis_tools.py`

**Interfaces:**
- Consumes: `analysis_spec._read_tsv(path_text, label) -> (columns, rows, sha256)`、`analysis_spec._require_column`
- Produces:
  - spec の `labels`: `{"path": str, "name": str, "feature_column": "feature_id", "label_column": "level", "order": ["A","B","C","D"]}`（`path` と `name` は必須。`name` は `^[A-Za-z0-9_]+$`）。正規化した spec の `labels` に入るが、`rule_hash` には入れない（印は選び方を変えない）。
  - `read_labels(labels: dict) -> dict` → `{"by_feature": {feature_id: 最も強い水準}, "sha256": str, "n_rows": int}`。空の水準の行は無視、`order` に無い水準は `DomainError("NOVEL_SPEC_INVALID", ..., {"field": "labels", "unknown_levels": [...]})`。
  - `attach_labels(rows: list[dict], by_feature: dict[str, str], *, name: str, order: list[str], top_n: int = 20) -> dict` → 各行に `labels.<name>`（水準か None）を足し（行の `merged_from` があれば群のどれかの最も強い水準）、`{"all": {水準: 件数, "none": 件数}, "top20": {...}}` を返す。並びは変えない。
  - ツールの戻り値と結果の JSON に `labels: {"name", "sha256", "order", "n_rows", "n_unmatched", "counts"}`（`n_unmatched` はデータセットの特徴 ID に無い行の ID の数）。行の TSV の末尾に `labels.<name>` 列。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_strongest_level_wins_for_a_feature_with_several_rows(tmp_path):
    path = _write(tmp_path, [("10", "C"), ("10", "A"), ("11", "D")])
    out = read_labels({"path": str(path), "name": "xenobiotic", "feature_column": "feature_id",
                       "label_column": "level", "order": ["A", "B", "C", "D"]})
    assert out["by_feature"] == {"10": "A", "11": "D"} and out["n_rows"] == 3

def test_unknown_level_is_refused_and_empty_level_is_ignored(tmp_path):
    assert read_labels(_spec(_write(tmp_path, [("10", "")])))["by_feature"] == {}
    with pytest.raises(DomainError) as info:
        read_labels(_spec(_write(tmp_path, [("10", "Z")])))
    assert info.value.code == "NOVEL_SPEC_INVALID" and info.value.details["unknown_levels"] == ["Z"]

def test_merged_pair_takes_the_strongest_member_label():
    rows = [{"feature_id": "5", "q": 0.01, "merged_from": ["5", "6"]}]
    attach_labels(rows, {"6": "B"}, name="xenobiotic", order=["A", "B", "C", "D"])
    assert rows[0]["labels.xenobiotic"] == "B"

def test_counts_all_and_top20_and_order_unchanged():
    rows = [{"feature_id": str(i), "q": i / 1000} for i in range(30)]
    before = [r["feature_id"] for r in rows]
    counts = attach_labels(rows, {"0": "A", "25": "B"}, name="x", order=["A", "B", "C", "D"])
    assert [r["feature_id"] for r in rows] == before
    assert counts["all"]["A"] == 1 and counts["all"]["B"] == 1 and counts["all"]["none"] == 28
    assert counts["top20"]["A"] == 1 and counts["top20"]["B"] == 0

def test_labels_do_not_change_rule_hash():
    base = {"statistics": {"result_id": "r"}, "q_max": 0.05}
    with_labels = {**base, "labels": {"path": "x.tsv", "name": "xenobiotic"}}
    assert rule_hash(normalize_novel_spec(base)) == rule_hash(normalize_novel_spec(with_labels))
```

`tests/test_reanalysis_tools.py` に 1 本: `test_novel_candidates_labels_column_and_counts(data, tmp_path)` — `labels` を付けた呼び出しで `n_rows` と行の並びが付けない呼び出しと同じ、`rows` の TSV のヘッダの末尾が `labels.xenobiotic`、`out["labels"]["counts"]["all"]["A"] == 1`、`out["labels"]["n_unmatched"] == 1`（TSV にデータセットに無い ID を 1 つ入れる）。

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_novel_labels.py tests/test_reanalysis_tools.py -k "label" -q`
Expected: FAIL（`read_labels` が無い / spec の `labels` が `NOVEL_SPEC_INVALID`）

- [ ] **Step 3: 実装する**

`normalize_novel_spec` の許すキーに `labels` を足して検証し（`path` は `_PATH_KEYS` の `path` なので `_strip_paths` が落とす）、`rule_hash` は `material.pop("labels", None)` してから hash する。`dataset_novel_candidates` は `mass_accuracy.novel_candidates` の後に `read_labels` → `attach_labels` を呼び、`labels` を結果と戻り値に足す。説明文に `labels` を 1 文足す。

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_novel_labels.py tests/test_reanalysis_tools.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

`docs/output_format/reanalysis.md` §8 の spec に `labels` の形・既定・強さの順・割れの組・`counts`（`all` と `top20`）・`rule_hash` に入らないことを書く。`USAGE.md` の行に 1 文。

```bash
git add metabolomix/analysis/mass_accuracy.py metabolomix/tools/reanalysis_tools.py tests/test_novel_labels.py tests/test_reanalysis_tools.py docs/output_format/reanalysis.md USAGE.md
git commit -m "feat(reanalysis): 新しい知見の候補に外部の印（labels）の列と件数を付ける"
```

---

### Task 4: 統計の行に群ごとの検出の例数・Fisher・detection_driven（#4）

**Files:**
- Modify: `metabolomix/analysis/spec_statistics.py`（新しい `detection_fisher_p` `detection_driven`）
- Modify: `metabolomix/analysis/analysis_spec.py`（`_run_analysis`）
- Modify: `metabolomix/tools/reanalysis_tools.py`（`ROW_COLUMNS` `_SIGNIFICANT` `dataset_analysis_run` の要約）
- Modify: `docs/output_format/reanalysis.md` §2、`USAGE.md`（`dataset_analysis_run` の行）
- Test: `tests/test_analysis_spec_detection.py`（新規）、`tests/test_reanalysis_tools.py`

**Interfaces:**
- Consumes: `_load` の `data["detected"]`（試料 × `data["feature_ids"]`。まとめた後は `_merge_plan` が「どちらかで検出」にした列）
- Produces:
  - `detection_fisher_p(detected_ref: int, n_ref: int, detected_test: int, n_test: int) -> float` — `scipy.stats.fisher_exact([[detected_test, n_test - detected_test], [detected_ref, n_ref - detected_ref]])` の両側 p。
  - `detection_driven(direction: str | None, detected_ref: int | None, n_ref: int, detected_test: int | None, n_test: int) -> bool | None` — `direction` が None か検出の数が None なら None。高い側（`test_higher` なら test）の検出 `< n_high / 2`、または低い側の検出 `== 0` なら True。
  - analysis-result.v1 の各行に `detected_reference` `detected_test`（int | None）・`detection_fisher_p`（float | None）・`detection_driven`（bool | None）。`n` は検定に入れた群の試料の数（`members[label]`）。外部の行列・検出状態の無いデータセット・`aggregate` の行・`not_in_matrix` の行は 4 つとも None（`detection_driven` は `status == "tested"` の行だけ）。
  - `ROW_COLUMNS` の末尾（`merged_from` の後）に `detected_reference` `detected_test` `detection_fisher_p` `detection_driven`。要約の解析ごとに `detection_available`（bool）と `n_detection_driven`。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_fisher_presence_absence():
    assert detection_fisher_p(0, 6, 6, 6) == pytest.approx(2 / 924)

@pytest.mark.parametrize("direction,ref,test,expected", [
    ("test_higher", 0, 2, True),     # 高い側 2/6 < 3
    ("test_higher", 0, 6, True),     # 低い側の検出 0
    ("test_higher", 3, 6, False),
    ("test_higher", 1, 3, False),    # 3 は 3 未満ではない
    ("ref_higher", 2, 0, True),      # 高い側は ref
    (None, 0, 6, None),
    ("test_higher", None, 6, None),
])
def test_detection_driven_rule(direction, ref, test, expected):
    assert detection_driven(direction, ref, 6, test, 6) is expected

def test_rows_carry_detection_counts_from_the_mask():
    ds = _make_ds([("1", 300.0, 2.0, ALL, "X", ADDUCT),
                   ("2", 400.0, 3.0, set(range(21, 42)), "Y", ADDUCT)], detected=True)
    rows = _rows(ds, {"comparison": {"reference": "A", "test": "B"}, "features": {"all": True, "merge_split": False},
                      "analyses": [{"analysis_id": "main", "test": {"kind": "welch"}}]})
    assert (rows["1"]["detected_reference"], rows["1"]["detected_test"]) == (21, 21)
    assert rows["1"]["detection_driven"] is False
    assert (rows["2"]["detected_reference"], rows["2"]["detected_test"]) == (0, 21)
    assert rows["2"]["direction"] == "test_higher" and rows["2"]["detection_driven"] is True
    assert rows["2"]["detection_fisher_p"] < 1e-6

def test_no_detection_state_gives_nulls(files):          # 外部の行列（tests/test_analysis_spec.py の files）
    row = _external_rows(files)[0]
    assert [row[k] for k in ("detected_reference", "detected_test",
                             "detection_fisher_p", "detection_driven")] == [None] * 4

def test_merged_column_counts_either_spot(): ...       # merge_split=True の割れの組: どちらかで検出の数
```

`tests/test_reanalysis_tools.py` に `test_analysis_run_rows_have_detection_columns(data)`: 戻り値の `rows` の TSV のヘッダの末尾 4 列と、要約の `detection_available is False`（合成の mzTab は `.arf` が無い）・`n_detection_driven == 0`。

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_analysis_spec_detection.py -q`
Expected: FAIL（`detection_fisher_p` が無い）

- [ ] **Step 3: 実装する**

`_run_analysis` で `col = {fid: j for j, fid in enumerate(data["feature_ids"])}` を作り、`data["detected"]` があれば `rows_ref` / `rows_test` の行で数える（`aggregate` のときは数えない）。行を組む最後（`direction` を決めた後）に 4 欄を足す。

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_analysis_spec_detection.py tests/test_analysis_spec.py tests/test_reanalysis_tools.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

`docs/output_format/reanalysis.md` §2 の行の欄の表に 4 欄（検出の例数は gap-fill を数えない、Fisher は検出の有無の 2×2 の両側、`detection_driven` の規則そのもの、判定は変えない＝使い方は手順書）。

```bash
git add metabolomix/analysis/spec_statistics.py metabolomix/analysis/analysis_spec.py metabolomix/tools/reanalysis_tools.py tests/test_analysis_spec_detection.py tests/test_reanalysis_tools.py docs/output_format/reanalysis.md USAGE.md
git commit -m "feat(reanalysis): 統計の行に群ごとの検出の例数・Fisher の p・detection_driven を足す"
```

---

### Task 5: 質量の窓の中心のずれ（`ppm_offset`）の正規化と新しい知見の候補（#6 その 1）

**Files:**
- Modify: `metabolomix/analysis/mass_accuracy.py`（新しい `normalize_ppm_offset` `offset_for`、`normalize_novel_spec` `novel_candidates`）
- Modify: `docs/output_format/reanalysis.md` §8
- Test: `tests/test_ppm_offset.py`（新規）

**Interfaces:**
- Produces:
  - `ppm_offset` の形（両ツール共通）: `{"default": 0.0, "bands": [{"max_mz": 200, "offset": -8.3}]}`。帯の選び方は `ppm` の帯と同じ（理論 m/z が `max_mz` 未満の最初の帯、無ければ `default`）。中心は「実測 − 理論」の ppm の符号つきの値。
  - `normalize_ppm_offset(raw, *, code: str, field: str = "ppm_offset") -> dict | None` — None なら None。未知のキー・有限でない数は `DomainError(code, ..., {"field": field})`。`bands` は `max_mz` の昇順、数は float。
  - `offset_for(mz: float, offset_spec: dict | None) -> float` — None なら 0.0。
  - novel spec の任意の `ppm_offset`。渡したときだけ正規化した spec に入る（省けば v8 と同じ `rule_hash`）。判定は `abs(row["ppm"] - offset_for(row["theoretical_mz"], spec)) <= _tolerance(...)`。行の `ppm` は実測 − 理論のまま。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_offset_band_edge_goes_to_next_band():
    spec = normalize_ppm_offset({"default": 1.0, "bands": [{"max_mz": 200, "offset": -8.0}]},
                                code="NOVEL_SPEC_INVALID")
    assert offset_for(199.999, spec) == -8.0 and offset_for(200.0, spec) == 1.0
    assert offset_for(150.0, None) == 0.0

@pytest.mark.parametrize("bad", [{"bands": [{"max_mz": 200}]}, {"default": float("nan")},
                                 {"default": 0, "extra": 1}, "x"])
def test_bad_offset_names_the_field(bad):
    with pytest.raises(DomainError) as info:
        normalize_ppm_offset(bad, code="NOVEL_SPEC_INVALID")
    assert info.value.details["field"] == "ppm_offset"

def test_novel_candidate_inside_shifted_window_only():
    ds, stats = _novel_ds(ppm=-9.0, theo=150.0)                # 同定の誤差 −9 ppm
    plain = normalize_novel_spec({"statistics": {"result_id": "r"}, "ppm": {"default": 5.0}})
    shifted = normalize_novel_spec({**_raw(plain), "ppm_offset": {"bands": [{"max_mz": 200, "offset": -8.3}]}})
    assert novel_candidates(ds, plain, stats, group_ref="A", group_test="B") == []
    [row] = novel_candidates(ds, shifted, stats, group_ref="A", group_test="B")
    assert row["ppm"] == pytest.approx(-9.0)

def test_rule_hash_unchanged_without_offset_and_changed_with_it():
    base = {"statistics": {"result_id": "r"}}
    assert rule_hash(normalize_novel_spec(base)) == _V8_HASH_OF_BASE      # 今の実装で計算した値を貼る
    assert rule_hash(normalize_novel_spec({**base, "ppm_offset": {"default": -3}})) != _V8_HASH_OF_BASE
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_ppm_offset.py -q`
Expected: FAIL（`normalize_ppm_offset` が無い）

- [ ] **Step 3: 実装する**（`_V8_HASH_OF_BASE` は実装の**前**に今のコードで計算してテストに貼っておく）

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_ppm_offset.py tests/test_reanalysis_tools.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

§8 に `ppm_offset` の形・帯・符号（中心 = 実測 − 理論の ppm）・`ppm` の幅はそのまま・`rule_hash` に入る（渡したときだけ）。

```bash
git add metabolomix/analysis/mass_accuracy.py tests/test_ppm_offset.py docs/output_format/reanalysis.md
git commit -m "feat(reanalysis): 新しい知見の候補の質量の窓の中心を帯ごとにずらせるようにする（ppm_offset）"
```

---

### Task 6: `dataset_feature_search` の `ppm_offset`（#6 その 2）

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`（`normalize_targets` `match_target`）
- Modify: `metabolomix/tools/dataset_analysis_tools.py`（`dataset_feature_search`・`_rule_targets`）
- Modify: `docs/output_format/mztab.md` §13.8、`USAGE.md`（`dataset_feature_search` の行）
- Test: `tests/test_ppm_offset.py`、`tests/test_dataset_feature_search_tool.py`

**Interfaces:**
- Consumes: Task 5 の `normalize_ppm_offset(raw, *, code, field)`・`offset_for(mz, spec)`
- Produces:
  - `normalize_targets(targets, *, dataset_polarity, ppm_offset: dict | None = None) -> list[dict]` — `ppm_offset`（正規化済み）を渡したときだけ、各目標に `"ppm_offset": {付加体: offset_for(理論 m/z, spec)}` を足す（`rule_hash` に入る。渡さなければ目標の辞書は今と同じ）。
  - `match_target` の質量の窓は `abs(観測の ppm − ppm_offset.get(adduct, 0.0)) <= ppm + _BOUNDARY_EPS`。候補の行の `ppm` は実測 − 理論のまま。
  - ツールの引数 `ppm_offset: dict | None = None`。`targets` が無い（窓検索）ときに渡せば `DATASET_BAD_REQUEST`（`field: ppm_offset`）。形の誤りは `DATASET_BAD_REQUEST`。`rule` に `ppm_offset`（正規化した spec か None）、`rule.targets[]` に丸めた `ppm_offset`（2 桁）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_match_target_uses_shifted_center():
    ds = _dataset({"20": {"mz": _shift(MH, -9.0), "rt": 2.0, "candidates": [_candidate()]}})
    offset = normalize_ppm_offset({"bands": [{"max_mz": 200, "offset": -8.0}]}, code="X")
    [plain] = normalize_targets([_gaba(ppm=5, rt_min=2.0)], dataset_polarity="positive")
    [shifted] = normalize_targets([_gaba(ppm=5, rt_min=2.0)], dataset_polarity="positive",
                                  ppm_offset=offset)
    assert match_target(ds, plain, max_results=5)["counts"]["mass"] == 0
    out = match_target(ds, shifted, max_results=5)
    assert out["status"] == "resolved" and out["candidates"][0]["ppm"] == pytest.approx(-9.0, abs=1e-6)
    assert rule_hash([plain]) != rule_hash([shifted])
    assert "ppm_offset" not in plain
```

`tests/test_dataset_feature_search_tool.py`: `test_ppm_offset_with_window_search_is_bad_request`（`feature_ids=["1"]` と `ppm_offset` → `DATASET_BAD_REQUEST`、`field == "ppm_offset"`）、`test_ppm_offset_appears_in_rule`（目標照合で `out["rule"]["ppm_offset"]["bands"][0]["offset"] == -8.0`、省いた呼び出しの `rule_hash` は今と同じ）。

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_ppm_offset.py tests/test_dataset_feature_search_tool.py -k offset -q`
Expected: FAIL（`normalize_targets` が `ppm_offset` を受けない）

- [ ] **Step 3: 実装する**

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_ppm_offset.py tests/test_feature_search.py tests/test_dataset_feature_search_tool.py tests/test_reanalysis_tools.py -q`
Expected: PASS（`test_claim_evidence_default_matches_feature_search` は `ppm_offset` を渡さないので変わらない）

- [ ] **Step 5: 文書を直してコミット**

`docs/output_format/mztab.md` §13.8 に `ppm_offset`（形は reanalysis §8 と同じ・窓検索では使えない・`ppm` 列は実測 − 理論のまま）。

```bash
git add metabolomix/analysis/feature_search.py metabolomix/tools/dataset_analysis_tools.py tests/test_ppm_offset.py tests/test_dataset_feature_search_tool.py docs/output_format/mztab.md USAGE.md
git commit -m "feat(feature-search): 目標照合の質量の窓の中心を帯ごとにずらす ppm_offset を足す"
```

---

### Task 6b: `dataset_claim_evidence` の走査にも `ppm_offset`（#6 その 3。設計者（本体のセッション）の判断 2026-11-04）

**Files:**
- Modify: `metabolomix/analysis/claim_evidence.py`（`normalize_claim_spec`・`build` の走査の目標）
- Modify: `metabolomix/tools/reanalysis_tools.py`（`dataset_claim_evidence` の説明文に 1 文）
- Modify: `docs/output_format/reanalysis.md` §7、`USAGE.md`（`dataset_claim_evidence` の行）
- Test: `tests/test_reanalysis_tools.py`

**Interfaces:**
- Consumes: Task 5 の `normalize_ppm_offset(raw, *, code, field)`、Task 6 の `normalize_targets(..., ppm_offset=...)`
- Produces: claim spec の任意の最上位のキー `ppm_offset`（形と帯の規則は Task 5 と同じ）。`normalize_claim_spec` が `normalize_ppm_offset(raw, code="CLAIM_SPEC_INVALID", field="ppm_offset")` で正規化し、**渡したときだけ**正規化した spec の `ppm_offset` に入れ、既定の目標と `scan` の全変種の `normalize_targets(..., ppm_offset=...)` に同じ spec を渡す（中心だけをずらし、幅は `scan.ppm` のまま）。渡さない spec は目標の辞書も正規化した spec も今と同じなので `rule_hash` は変わらない。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_claim_evidence_ppm_offset_matches_feature_search(data):
    handle = _open(data)
    offset = {"bands": [{"max_mz": 200, "offset": -8.0}]}
    claim = {**CLAIM, "ppm_offset": offset}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=handle))
    search = json.loads(T.dataset_feature_search(
        targets=[{"target_id": "C1", **CLAIM["target"]}], ppm_offset=offset, dataset_handle=handle))
    assert out["default_selection"]["status"] == search["targets"][0]["status"]
    assert out["default_selection"]["selected_feature_id"] == search["targets"][0]["selected_feature_id"]
    assert out["rule"]["ppm_offset"]["bands"][0]["offset"] == -8.0

def test_claim_evidence_without_offset_keeps_rule_hash(data):
    handle = _open(data)
    before = _CLAIM_RULE_HASH_V8          # 実装の前に今のコードで計算して貼る
    assert json.loads(R.dataset_claim_evidence(CLAIM, dataset_handle=handle))["rule_hash"] == before
    shifted = json.loads(R.dataset_claim_evidence({**CLAIM, "ppm_offset": {"default": -3}},
                                                  dataset_handle=handle))
    assert shifted["rule_hash"] != before

def test_claim_evidence_bad_offset_is_spec_invalid(data):
    out = json.loads(R.dataset_claim_evidence({**CLAIM, "ppm_offset": {"default": "x"}},
                                              dataset_handle=_open(data)))
    assert out["error"]["code"] == "CLAIM_SPEC_INVALID"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_reanalysis_tools.py -k "claim_evidence and offset" -q`
Expected: FAIL（`spec の未知のキー: ['ppm_offset']`）

- [ ] **Step 3: 実装する**（`allowed` に `ppm_offset`。`build` の走査の `normalize_targets([variant], ...)` にも同じ spec を渡す）

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_reanalysis_tools.py tests/test_ppm_offset.py -q`
Expected: PASS（`test_claim_evidence_default_matches_feature_search` も変わらない）

- [ ] **Step 5: 文書を直してコミット**

§7 の spec に `ppm_offset`（`dataset_feature_search` と同じ形・同じ帯の規則・走査の全変種に効き、幅は変えない・渡したときだけ `rule_hash` に入る）。

```bash
git add metabolomix/analysis/claim_evidence.py metabolomix/tools/reanalysis_tools.py tests/test_reanalysis_tools.py docs/output_format/reanalysis.md USAGE.md
git commit -m "feat(claim-evidence): 証拠の束の走査にも質量の窓の中心のずれ（ppm_offset）を渡せるようにする"
```

---

### Task 7: 検出状態の接合で小さな外れを受け入れる（#10）

**Files:**
- Modify: `metabolomix/mztab/evidence.py`（定数 `MZ_OUTLIER_SMALL_DELTA = 0.02`、`compare_feature_mz`）
- Test: `tests/test_dataset_evidence.py`

**Interfaces:**
- Produces: `compare_feature_mz` の判定 `ok = not over_hard_limit and (over_fraction <= MZ_OUTLIER_MAX_FRACTION or (n_over > 0 and worst < MZ_OUTLIER_SMALL_DELTA))`。`detail` に `small_outlier_limit: 0.02` と `accepted_by`（`"fraction"` / `"small_outliers"` / None）。`assay_evidence` も同じ関数を使うので同時に効く。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_build_evidence_accepts_many_small_outliers():          # 陰イオン 36/2547（1.41%）・最大 0.0188 Da
    spots, feature_mz = _many_spots(100, {i: 0.0188 for i in range(2)})
    result = evidence.build_evidence(evidence.normalize_arf_spots(spots), feature_mz=feature_mz,
                                     sample_names=["s1", "s2"])
    assert result["status"] == "ok" and result["detail"]["accepted_by"] == "small_outliers"

def test_one_outlier_at_or_beyond_002_da_still_rejects():
    spots, feature_mz = _many_spots(100, {0: 0.015, 1: 0.020})
    result = evidence.build_evidence(evidence.normalize_arf_spots(spots), feature_mz=feature_mz,
                                     sample_names=["s1", "s2"])
    assert result["status"] == "rejected" and result["reason"] == "mz_mismatch"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_dataset_evidence.py -k "small_outlier or 002" -q`
Expected: FAIL（1 本目が rejected）

- [ ] **Step 3: 実装する**（定数の上のコメントに MTBLS9365 陰イオンの実測を書く）

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_dataset_evidence.py tests/test_assay_evidence.py -q`
Expected: PASS（`test_build_evidence_rejects_when_many_features_exceed_tolerance` は 0.05 Da なので棄却のまま）

- [ ] **Step 5: コミット**

```bash
git add metabolomix/mztab/evidence.py tests/test_dataset_evidence.py
git commit -m "fix(evidence): m/z の外れがすべて 0.02 Da 未満なら検出状態の接合を受け入れる"
```

---

### Task 8: 割れの組で残す側を受け持ちで決める（#11）

**Files:**
- Modify: `metabolomix/analysis/alignment_qc.py`（`find_split_pairs` の `merge_groups`、新しい `_keep_first`、`usable_merge_groups`、`merge_matrix_columns`）
- Modify: `metabolomix/analysis/feature_search.py`（`merge_map`・`merged_dataset` と `match_target` の説明）
- Modify: `metabolomix/analysis/analysis_spec.py`（`merge_groups_hash`）
- Modify: `metabolomix/tools/dataset_analysis_tools.py`・`metabolomix/tools/reanalysis_tools.py`・`metabolomix/analysis/mass_accuracy.py` の説明文（「残すのは小さい ID」）
- Modify: `docs/output_format/reanalysis.md` §9（と §2・§8 の `merged_from` の記述）
- Test: `tests/test_alignment_qc.py`、`tests/test_analysis_spec.py`

**Interfaces:**
- Produces:
  - 約束: **割れの群の並びは残す ID が先頭**（`merge_groups[].features`・`usable_merge_groups` の戻り・`merged_from` の値）。`merge_groups[]` に `keep`（残す ID）・`keep_basis`（`"n_owned"` / `"smaller_id_tie"`）・`n_owned`（`{特徴 ID: 受け持ちの数}`）を足す。`groups[]`（報告用の連結成分）は ID の小さい順のまま。
  - `_keep_first(pair: dict) -> tuple[list[str], str]` — 組の `features` と `n_owned` から（残す ID が先頭の並び, keep_basis）。
  - `usable_merge_groups(qc, feature_ids)` — `keep` があればそれを先頭、無ければ（古い JSON）`id_order`。
  - `merge_matrix_columns(matrix, feature_ids, groups, *, detected=None)` — 群の中で行列にある最初の ID を残す（渡した並びに従う）。
  - `merge_map(groups) -> dict[str, list[str]]` → `{group[0]: list(group)}`。
  - `merge_groups_hash(groups) -> str` → `canonical_hash(sorted([str(g[0]), *sorted(str(f) for f in g[1:])] for g in groups))`（残す ID が最小なら今と同じ hash）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_keep_side_is_the_one_owning_more_samples():                      # G142（6965-v8）
    ds = _make_ds([("1265", 650.6000, 20.61, FIRST_TWO, None, ADDUCT),
                   ("1266", 650.6050, 20.00, ALL - FIRST_TWO, "CE 18:2", ADDUCT)])
    [merge] = Q.find_split_pairs(ds)["merge_groups"]
    assert merge["features"] == ["1266", "1265"] and merge["keep"] == "1266"
    assert merge["keep_basis"] == "n_owned" and merge["n_owned"] == {"1266": 40, "1265": 2}

def test_tie_keeps_the_smaller_id():
    evens, odds = set(range(0, N, 2)), set(range(1, N, 2))
    ds = _make_ds([("31", 500.0, 5.00, odds, "PC 34:1", ADDUCT),
                   ("30", 500.001, 5.10, evens, "PC 34:1", ADDUCT)], groups=GROUPS_HALVES)
    [merge] = Q.find_split_pairs(ds)["merge_groups"]
    assert merge["features"] == ["30", "31"] and merge["keep_basis"] == "smaller_id_tie"

def test_old_qc_json_without_keep_keeps_the_smaller_id():
    qc = {"merge_groups": [{"features": ["1265", "1266"], "mergeable": True}]}
    assert Q.usable_merge_groups(qc, ["1265", "1266"])[0] == [["1265", "1266"]]

def test_merge_groups_hash_unchanged_when_keep_is_smallest():
    assert merge_groups_hash([["1", "2"]]) == canonical_hash([["1", "2"]])
    assert merge_groups_hash([["2", "1"]]) != merge_groups_hash([["1", "2"]])
```

`tests/test_analysis_spec.py` に `test_merge_split_keeps_the_owning_spot`（上の G142 の ds で行の ID が `"1266"`、`merged_from == ["1266", "1265"]`）。既存の `test_merge_matrix_columns_takes_per_sample_max` は「渡した並びの先頭を残す」に直す（`[["c","a"]]` は `"c"` を残し `merged_from == {"c": ["c","a"]}`、数の ID の例は `[["9","10"]]` を渡す）。

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_alignment_qc.py tests/test_analysis_spec.py -k "keep or hash or merge" -q`
Expected: FAIL（`keep` が無い）

- [ ] **Step 3: 実装する**

`find_split_pairs` で `merges` を作った後、各 merge の組（`pairs` のうち `merged` で features が同じもの）に `_keep_first` を当てて `features` を並べ替え、`keep` `keep_basis` `n_owned` を足す。`merge_matrix_columns` は `present = [f for f in group if str(f) in position]`（`id_order` を外す）。`novel_candidates` の `absorbed = group[1:]` と `dataset_novel_candidates` の `{g[0]: g}` はこの約束のまま正しい（変えない）。

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_alignment_qc.py tests/test_analysis_spec.py tests/test_feature_search.py tests/test_dataset_feature_search_tool.py tests/test_reanalysis_tools.py -q`
Expected: PASS（既存の 1265/1266 の例は 1265 が多く受け持つので結果は変わらない）

- [ ] **Step 5: 文書を直してコミット**

§9 に `keep` `keep_basis` `n_owned`・「群の並びは残す ID が先頭」・古い JSON の扱い。§2・§8 と各ツールの説明の「残すのは小さい ID」を「受け持ちの多い側（同数なら小さい ID）」に。

```bash
git add metabolomix/analysis/alignment_qc.py metabolomix/analysis/feature_search.py metabolomix/analysis/analysis_spec.py metabolomix/analysis/mass_accuracy.py metabolomix/tools/dataset_analysis_tools.py metabolomix/tools/reanalysis_tools.py tests/test_alignment_qc.py tests/test_analysis_spec.py docs/output_format/reanalysis.md
git commit -m "fix(alignment-qc): 割れの組で受け持ちの多いスポットを残す（同数なら小さい ID）"
```

---

### Task 9: 証拠の束の代表スペクトルの数（#12）

**Files:**
- Modify: `metabolomix/analysis/claim_evidence.py`（`_msms_for`・`_selection`・`build`、新しい `_spectrum_available`）
- Modify: `docs/output_format/reanalysis.md` §7、`USAGE.md`（`dataset_claim_evidence` の行）
- Test: `tests/test_reanalysis_tools.py`

**Interfaces:**
- Produces:
  - `_spectrum_available(spot: dict) -> bool | None` — spot の照会の `status` が `ok` なら True、`not_found` なら False、それ以外（`spot_not_found` `dcl_precursor_mismatch` `spot_id_not_integer`、ハンドルに `.arf2` / `.dcl` が無いときの封筒）は None。
  - `_msms_for` は spot の照会を**ライブラリの有無の確認より先に**行い、戻り値に `spectrum_available` を必ず入れる（ライブラリが無くても `{"status": "no_library", "spectrum_available": ...}`）。
  - `_selection(match, verdicts, spectra)`（`spectra` は特徴 ID → `spectrum_available`。MS/MS を頼まなければ None）の `counts` に:
    - `rule_qualified` — 今の `msms_available` と同じ値（規則で qualified の数）。
    - `msms_available` — 非推奨の別名（同じ値。claim-evidence.v1 の互換のため残す）。
    - `msms_spectrum_available` — 質量の窓に入った候補（`counts.mass` と同じ母集団。`merge_partner` の行を除き、割れの組は 1 つ）のうち、組のどれかのスポットの `spectrum_available` が True の数。RT の窓で絞らない（6965-v8 の C4・C13 は RT の窓の外の 1265 に代表 MS/MS があったのに 0 と読めたのが事象）。MS/MS を頼まなければ None。`max_candidates` で打ち切られて照会していない候補は数えない。

- [ ] **Step 1: 失敗するテストを書く**

`test_counts_qualified_is_msms_verdict` に足す:

```python
    assert strict["rule_qualified"] == 1 and strict["msms_spectrum_available"] == 1
    assert bare["default_selection"]["counts"]["msms_spectrum_available"] is None
```

新しい 2 本:

```python
def test_spectrum_available_counts_mass_window_even_when_none_rule_qualified(data, library):
    # 目標の RT を特徴 20（RT 2.0）から外す: 規則で qualified は 0、質量の窓の 20 の代表 MS/MS は数える
    target = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 9.0, "adduct_tiers": [["[M+H]+"]]}
    claim = {"claim_id": "C1", "target": target, "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=_open(data, library="auto")))
    counts = out["default_selection"]["counts"]
    assert counts["rule_qualified"] == 0 and counts["msms_spectrum_available"] == 1

def test_spectrum_available_without_library(data):
    target = {"inchikey": GABA_KEY, "formula": "C4H9NO2", "rt_min": 2.0, "adduct_tiers": [["[M+H]+"]]}
    claim = {"claim_id": "C1", "target": target, "msms": {"top_n": 3, "rules": {"min_fragments_excl_precursor": 2}}}
    out = json.loads(R.dataset_claim_evidence(claim, dataset_handle=_open(data)))
    c = {x["feature_id"]: x for x in out["candidates"]}["20"]
    assert c["msms"]["status"] == "no_library" and c["msms"]["spectrum_available"] is True
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_reanalysis_tools.py -k "counts or spectrum" -q`
Expected: FAIL（`rule_qualified` が無い）

- [ ] **Step 3: 実装する**

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_reanalysis_tools.py tests/test_reanalysis_small_fixes.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

§7 の `counts` の表を 3 つ（`rule_qualified`・`msms_available` は非推奨の別名・`msms_spectrum_available`）と、候補の `msms.spectrum_available` の 3 値に。

```bash
git add metabolomix/analysis/claim_evidence.py tests/test_reanalysis_tools.py docs/output_format/reanalysis.md USAGE.md
git commit -m "feat(claim-evidence): 代表スペクトルの有無で数える msms_spectrum_available と rule_qualified を足す"
```

---

### Task 10: 実装済みの 2 項目を確かめる（#14・#17）

**Files:** なし（確かめるだけ。失敗したら止まって報告する）

- [ ] **Step 1: #14 を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_dataset_handle_tools.py -k "pipeline_arf2 or newest_arf2" -q`
Expected: PASS 2 本（`pipeline_path` で `input/AlignResult-*.arf2` と同じバッチの `.dcl` を見つける）

- [ ] **Step 2: #17 を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_feature_search.py -k "acetate" -q`
Expected: PASS 2 本（`[M+CH3COO]-` が陰イオンの第 2 層にあり、質量が +59.013851）

- [ ] **Step 3: 記録する**

main ツリーの `docs/HISTRY.md` の末尾に、#14 は fd70f6d・#17 は 31c4b6f で済んでいたこと（protocol v9 の設計書の表は修正の前の issues.md から作られた）を追記する。コミットは無い（追跡外）。

---

### Task 11: パスウェイの run の manifest に要約と隣接ペアを書く（#15）

**Files:**
- Modify: `metabolomix/pathway/tools.py`（`pathway_project`）
- Modify: `docs/output_format/pathway.md`
- Test: `tests/test_pathway_tools.py`

**Interfaces:**
- Produces: `manifest.json` に `notes`（`{"adjacency": _ADJACENCY_NOTE, "handoff": _handoff_note(table.report) か None}`）と `opposite_pairs`（`{パスウェイ: {"rows": [ペアの辞書（_PAIR_COLUMNS の欄）], "undeterminable": int}}`）。非有限の float は null。戻り値は変えない。

- [ ] **Step 1: 失敗するテストを書く**

`test_project_writes_run_and_returns_tables` に足す:

```python
    assert manifest["notes"]["adjacency"] == body["notes"]["adjacency"]
    assert manifest["notes"].get("handoff") == body["notes"].get("handoff")
    pairs = manifest["opposite_pairs"]["remodeling"]
    assert pairs["undeterminable"] == body["pathways"]["remodeling"]["undeterminable_pairs"]
    assert len(pairs["rows"]) == len(body["pathways"]["remodeling"]["opposite_pairs_tsv"].splitlines()[1:])
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_tools.py::test_project_writes_run_and_returns_tables -q`
Expected: FAIL（`KeyError: 'notes'`）

- [ ] **Step 3: 実装する**（ペアの行は `compare.opposite_pairs` の戻りをそのまま `per_pathway` を作るループで控える）

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_tools.py tests/test_pathway_remote.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

`docs/output_format/pathway.md` に、manifest の `notes` と `opposite_pairs`（MCP の時間切れで戻り値が失われても run に残る）。

```bash
git add metabolomix/pathway/tools.py tests/test_pathway_tools.py docs/output_format/pathway.md
git commit -m "feat(pathway): run の manifest に handoff の注意と逆向きの隣接ペアを書く"
```

---

### Task 12: `pathway_plot` は既定で場所と要約を返す（#20）

**Files:**
- Modify: `metabolomix/pathway/tools.py`（`pathway_plot`）
- Modify: `docs/output_format/pathway.md`（「図（pathway_plot）」）、`USAGE.md`（`pathway_plot` の行）
- Test: `tests/test_pathway_tools.py`

**Interfaces:**
- Produces: `pathway_plot(run_id: str | None = None, pathways: list[str] | None = None, plot_type: str = "dot", inline_image: bool = False)`。PNG を常に run のディレクトリへ `<pathway>.png` で書く。`inline_image=False`（既定）は `json_payload({"run_id", "group_a", "group_b", "backend", "figures": [{"pathway", "pdf", "png", "width", "height", "mapped_nodes"}], "note": _FIGURE_NOTE})` の文字列だけ。`True` は `[同じ JSON の文字列, *Image]`（従来の画像）。`inline_image` が bool でなければ `PATHWAY_BAD_REQUEST`。

- [ ] **Step 1: 失敗するテストを書く**

`test_plot_renders_pngs` を 2 本に分ける:

```python
def test_plot_returns_paths_and_summary_by_default(lipid_arf, fake_measure, monkeypatch):
    # run の作成と bridge.run_bridge の差し替え（render）は今の test_plot_renders_pngs と同じ
    out = tools.pathway_plot(pathways=["ceramide"])
    assert isinstance(out, str)
    body = json.loads(out)
    [figure] = body["figures"]
    assert figure["pdf"].endswith("ceramide.pdf") and Path(figure["png"]).is_file()
    assert figure["mapped_nodes"] == manifest["mapped_nodes"]["ceramide"]
    assert body["note"] == tools._FIGURE_NOTE and len(out) < 2000

def test_plot_inline_image_keeps_the_old_payload(lipid_arf, fake_measure, monkeypatch):
    out = tools.pathway_plot(pathways=["ceramide"], inline_image=True)
    assert isinstance(out, list) and isinstance(out[1], Image) and len(out) == 2
```

- [ ] **Step 2: 失敗を確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_tools.py -k plot -q`
Expected: FAIL（戻り値が list）

- [ ] **Step 3: 実装する**

- [ ] **Step 4: 通ることを確かめる**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_tools.py tests/test_tool_annotations.py tests/test_server_registration.py -q`
Expected: PASS

- [ ] **Step 5: 文書を直してコミット**

`pathway.md` と `USAGE.md` の `pathway_plot` に、既定が場所と要約に変わったこと（画像はファイルだけ。設計者（本体のセッション）の判断 2026-11-04）、チャットに出すときは `inline_image=true`（1 枚約 2,400 画像トークン）。

```bash
git add metabolomix/pathway/tools.py tests/test_pathway_tools.py docs/output_format/pathway.md USAGE.md
git commit -m "feat(pathway): pathway_plot は既定で図のファイルの場所と要約だけを返す"
```

---

### Task 13: 全体の確認と記録

**Files:**
- Modify: `docs/superpowers/specs/2026-11-04-evidence-bundle-v3-samples-design.md`（状態を「実装済み」に）
- 追跡外: main ツリーの `docs/HISTRY.md`・`docs/task.md`、vault の `ms-data-parser-flow.md`

- [ ] **Step 1: 全テスト**

Run: `.venv/Scripts/python.exe -m pytest tests -q`
Expected: すべて PASS（`FAILED` 0）

- [ ] **Step 2: 実データで v3 を書く（設計 §6）**

ST003841・ST002901 の Console の出力（`C:\Users\yuu18\datasets\reanalysis\<accession>\` 配下。リポジトリに入れない）を `dataset_handle_open(pipeline_path=...)` で開き、2 群比較の結果から `dataset_export_evidence_bundle` で v3 を書く。manifest の `intensities_stage`・`detection_unavailable_reason`・行数を `docs/HISTRY.md` に記録する。massbank-context の `load_evidence_bundle` で読めることの確認は、massbank-context 側の v3 の実装が入ってから（設計書 §8 の陽性の対照）行い、それまでは `docs/task.md` に残す。

- [ ] **Step 3: 記録**

`docs/HISTRY.md` の末尾に日付見出しで、この計画の各タスクの結果・Task 2 の規模の実測・設計との食い違いの決着を追記。`docs/task.md` に「後に回すもの」（下）を足す。`pathway_plot` の既定の戻り値が変わったので、vault の `ms-data-parser-flow.md` の図の段を直し、出典行の日付を書き直す。設計書の状態を「実装済み（2026-11-xx）」に。

- [ ] **Step 4: コミット**

```bash
git add docs/superpowers/specs/2026-11-04-evidence-bundle-v3-samples-design.md
git commit -m "docs: 根拠バンドル v3 の設計書を実装済みにする"
```

## 決定（設計者の判断 2026-11-04。ユーザーの確認待ち）

- ARF 経路（`arf_export_evidence_bundle`）は v2 のまま。
- `dataset_claim_evidence` の走査にも `ppm_offset` を渡せるようにする（Task 6b）。
- `pathway_plot` の既定を「場所と要約」にする変更はそのまま入れ、`inline_image=True` で従来どおり。出力の文書（`pathway.md`）と `USAGE.md` に書く（Task 12）。

## 後に回すもの

- `counts.msms_available`（非推奨の別名）を消す: v9 では残し、次の版（claim-evidence.v2）で消す。reanalysis-study の triage が `rule_qualified` に移ってから。
