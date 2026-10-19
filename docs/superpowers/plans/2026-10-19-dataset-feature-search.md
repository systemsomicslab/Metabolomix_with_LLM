# dataset_feature_search 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** DatasetState の特徴を窓または目標（分子式・InChIKey・規則）で引き、候補・判定（強度で選ばない）・全特徴の差次的統計を返す MCP ツール `dataset_feature_search` を足し、pipeline の差次的結果を session へ引き継げるようにする。

**Architecture:** 判定と統計の結合は新しい純ロジック `metabolomix/analysis/feature_search.py`（session も MCP も import しない）に置き、MCP 公開層は `metabolomix/tools/dataset_analysis_tools.py` に 1 関数足す。pipeline の差次的結果は、来歴に新しく刻む `dataset_hash` と群の指紋の一致を条件に `session_handoff.load_run_dataset` が `ds.results` へ戻し、`result_state.is_current` に引き継ぎ用の分岐を足す。

**Tech Stack:** Python 3.14、FastMCP、numpy、pytest。

**Spec:** `docs/superpowers/specs/2026-10-19-dataset-feature-search-design.md`（以下「spec」。§ は spec の節）

## Global Constraints

- 判定は強度・ppm・RT の近さ・並び順で 1 件を選ばない。qualified が 2 件以上なら `needs_input`（`multiple_qualified_candidates`）、0 件なら `needs_input`（`no_qualified_candidate`）。
- `bind_features`（pipeline v2）の挙動を変えない。`tests/test_feature_bindings.py` は無変更で通る。
- 差次的エクスポートの列・メタ行を変えない（`CONTRACT_VERSION` は 1 のまま）。
- 戻り値は `json_payload()`、一覧は TSV（列名 1 回）、float は丸める（m/z 4 桁・ppm 2 桁・RT 3 桁・log2FC 4 桁・p/q/平均は有効 4 桁）。ツールは `structured_output=False`、`ToolAnnotations(readOnlyHint=True)`。
- 前提状態が無いときは `missing_state` 封筒（`metabolomix/core/mcp_errors.py`）。自分自身を `required_tools` に入れない。
- `metabolomix/pipeline/{engine,worker,service,recovery,store}.py` はグローバル session を import しない（AST テストが固定）。
- 状態の参照は正準モジュール（`session_state.session`）経由。`from ... import DATA_DIR` を書かない。
- fixture はテスト自身が作る。テストはリポジトリルートから `.venv/Scripts/python.exe -m pytest tests -q`。
- コミットは pre-commit が全テストを回す（数分）。自動化された `git commit` はバックグラウンドで実行し、出力をファイルに落とす。実行中は作業ツリーを編集しない。
- 新規の `.txt` は `.gitignore` で無言に追跡漏れする（このタスク群では作らない）。
- 目標は 1 回 50 件まで。目標ごとの候補・窓検索の行は `max_results`（既定 20）件まで。`ppm` 既定 10、`rt_tol_min` 既定 0.3。
- 付加体の層の既定（spec §4.1 そのまま）:
  - positive: `[["[M+H]+", "[M]+"], ["[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M+H-NH3]+", "[2M+H]+", "[M+2H]2+"]]`
  - negative: `[["[M-H]-"], ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[2M-H]-", "[M+Na-2H]-"]]`

## Review Focus

1. **再開・訂正後の run（同じ `differential:<cid>` が複数 ref）**: 引き継ぎは今の ref だけを読み、古い ref の hash 不一致で `RESULT_INTEGRITY_MISMATCH` を出さない → Task 3 のテスト `test_handoff_uses_only_current_refs_after_reexport`。
2. **永続化で NaN が null になった p**: 戻した結果の p が `None` の特徴は `not_testable`、例外にならない → Task 7 のテスト `test_null_p_from_persisted_result_is_not_testable`。
3. **SME の InChIKey が取れない読み込み（MSP の表なし・RDKit なし）**: SME 付きの特徴が全部落ちず、`sme_identity_unresolved` つきで次の水準へ進む → Task 6 のテスト `test_unresolved_sme_identity_is_not_a_contradiction`。
4. **1 つの特徴が 2 つの付加体仮説の質量窓に入る**: 候補は（特徴, 付加体）単位で 2 行になり、同じ特徴が二重に qualified 数えされて偽の `multiple_qualified_candidates` にならない → Task 6 のテスト `test_same_feature_under_two_adducts_counts_once_per_tier`（違う層なら上の層だけが使われる。同じ層なら特徴単位で数える）。
5. **検出マスクの列数が試料数と合わない**: 検出例数は空欄＋注記で、別の試料を数えない → Task 7 のテスト `test_detection_counts_skipped_when_mask_shape_mismatches`。

---

### Task 1: 差次的結果の来歴に `dataset_hash` を刻む

**Files:**
- Modify: `metabolomix/analysis/dataset_service.py`（`compare_dataset`、`group_fingerprint` を刻んでいる箇所の隣）
- Test: `tests/test_result_state.py`

**Interfaces:**
- Produces: すべての差次的結果の `result["provenance"]["dataset_hash"] == result_state.dataset_fingerprint(ds)`（session の `dataset_differential` と pipeline の `run_comparison` の両方が `compare_dataset` を通る）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_differential_provenance_records_the_dataset_fingerprint():
    ds = _preprocessed_ds()            # このファイルの既存ヘルパに倣って前処理済み ds を作る
    a, b = _two_groups(ds)
    result = dataset_service.compare_dataset(ds, a, b)
    assert result["provenance"]["dataset_hash"] == result_state.dataset_fingerprint(ds)
    assert "dataset_hash" not in result["provenance"]["effective_parameters"]
```

- [ ] **Step 2: 落ちることを確認**: `.venv/Scripts/python.exe -m pytest tests/test_result_state.py -k dataset_fingerprint -q` → FAIL（KeyError `dataset_hash`）
- [ ] **Step 3: `compare_dataset` で `result["provenance"]["dataset_hash"] = result_state.dataset_fingerprint(ds)` を足す**（`input_fingerprint` の中身は変えない）
- [ ] **Step 4: 通ることを確認**（同じコマンド → PASS）、続けて `-k "result_state or dataset_analysis"` 相当の既存テストも PASS
- [ ] **Step 5: コミット** `feat(provenance): 差次的結果にデータの指紋を刻む`

### Task 2: 引き継いだ結果の有効性（`is_current` の分岐）

**Files:**
- Modify: `metabolomix/analysis/result_state.py`（`is_current`）
- Test: `tests/test_result_state.py`

**Interfaces:**
- Consumes: Task 1 の `provenance["dataset_hash"]`。
- Produces: `provenance["handoff"]`（dict。中身は Task 3 が書く `{"pipeline_root": str, "comparison_id": str}`）を持つ結果は、`dataset_id` と親の前処理 ID を見ず、`dataset_hash == dataset_fingerprint(ds)` かつ（`group_fingerprint` を持つなら）`group_fingerprint == group_fingerprint(ds)` で有効と判定する。`handoff` の無い結果の判定は今のまま。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_handed_off_result_is_current_without_session_preprocess():
    ds = _ds_with_metadata()
    result = _handed_off_result(ds)    # provenance: dataset_id="ds_other", parent_ids=["pp_other"],
                                       # dataset_hash=dataset_fingerprint(ds), group_fingerprint=group_fingerprint(ds),
                                       # handoff={"pipeline_root": "x", "comparison_id": "c1"}
    assert ds.preprocess_id is None
    assert result_state.is_current(ds, result) is True

def test_handed_off_result_goes_stale_when_groups_change():
    ds = _ds_with_metadata()
    result = _handed_off_result(ds)
    ds.sample_metadata_rows = _rows_with_a_group_changed(ds)
    assert result_state.is_current(ds, result) is False

def test_handed_off_result_with_other_data_is_not_current():
    ds = _ds_with_metadata()
    result = _handed_off_result(ds)
    result["provenance"]["dataset_hash"] = "0" * 64
    assert result_state.is_current(ds, result) is False

def test_session_result_rules_are_unchanged():
    # handoff の無い結果で dataset_id が違えば今まで通り False
    ...
```

- [ ] **Step 2: 落ちることを確認** → 最初のテストが FAIL
- [ ] **Step 3: `is_current` の冒頭（来歴なし判定の直後）に `handoff` 分岐を足す**
- [ ] **Step 4: `tests/test_result_state.py` 全体が PASS**
- [ ] **Step 5: コミット** `feat(result_state): run から引き継いだ差次的結果の有効性をデータと群の指紋で判定する`

### Task 3: run の差次的結果を session へ戻す

**Files:**
- Modify: `metabolomix/pipeline/session_handoff.py`（`load_run_dataset`）
- Modify: `metabolomix/tools/mztab_tools.py`（`dataset_load` の pipeline 経路の要約文）
- Test: Create `tests/test_pipeline_differential_handoff.py`

**Interfaces:**
- Consumes: Task 1・2。`store.current_result_refs(record)`、`store.read_result_data(pipeline_root, refs, output_name)`。
- Produces:
  - `load_run_dataset` が、`current_result_refs(record)` のうち `output_name` が `differential:` で始まる ref を読み、指紋が合うものを `ds.results[result_id]` に登録（result_id は保存時のまま、`provenance["handoff"] = {"pipeline_root": str(pipeline_root), "comparison_id": cid}` を足す、`ds.last_differential` は触らない）。群の指紋の比較は `ds.sample_metadata_rows` を戻した**後**に行う。
  - `ds.handoff_report = {"restored": [{"comparison_id", "result_id", "group_a", "group_b"}], "rejected": [{"comparison_id", "result_id", "reason"}]}`。`reason` は `legacy_result_without_dataset_hash` / `dataset_mismatch` / `group_mismatch`。`DatasetState.__init__` に `handoff_report: dict | None = None` を足す。
  - `dataset_load` の要約に、戻した比較を TSV（`comparison_id result_id group_a group_b`）、戻せなかったものを理由つきで書く。

- [ ] **Step 1: 調査（コードは書かない）**: `pipeline_resume`（`metabolomix/pipeline/recovery.py` の `prepare_resume` と `updates` の扱い）で、完了済みの `differential` 工程を再計算させる手順があるかを確かめ、結果を `docs/HISTRY.md` に 1 行で記録する。あれば `legacy_result_without_dataset_hash` の注記にその手順を、無ければ「session で `dataset_preprocess` → `dataset_differential`」を案内文として使う。
- [ ] **Step 2: 失敗するテストを書く**。run は `tests/test_pipeline_matrix_handoff.py` の `_completed_run` で作り、その run の `console_job_path` と同じ job から `mztab_loading.load_dataset_state(job_path=...)` で ds を作る → run の `sample_manifest` の rows を `ds.sample_metadata_rows` に当て → 前処理と `compare_dataset` で結果を作る → `report.persist_result(pipeline_root, {"output_name": "differential:c1", "kind": "differential", "result_id": ..., "data": result, "parent_ids": [...]})` の ref を `record["results"]` に追記して `store.save_run`。

```python
def test_loading_a_run_restores_its_differential_results(tmp_path): ...
    # ds.results[result_id] があり、provenance.handoff.comparison_id == "c1"、
    # result_state.is_current(ds, that) is True、ds.last_differential is None
def test_legacy_result_without_dataset_hash_is_reported_not_restored(tmp_path): ...
    # provenance から dataset_hash を消して保存 → rejected に reason="legacy_result_without_dataset_hash"
def test_result_from_other_groups_is_reported_not_restored(tmp_path): ...
    # 保存後に run の sample_manifest を群違いで再保存 → reason="group_mismatch"
def test_handoff_uses_only_current_refs_after_reexport(tmp_path): ...
    # 同じ "differential:c1" を 2 回 persist（2 回目は再計算した結果）→ 例外にならず、2 回目の result_id だけが戻る
def test_dataset_load_summary_lists_restored_and_rejected(tmp_path): ...
```

- [ ] **Step 3: 落ちることを確認**: `.venv/Scripts/python.exe -m pytest tests/test_pipeline_differential_handoff.py -q` → FAIL
- [ ] **Step 4: 実装**（比較 ID は `output_name.split(":", 1)[1]`）
- [ ] **Step 5: PASS を確認**。あわせて `tests/test_pipeline_matrix_handoff.py` `tests/test_pipeline_engine.py`（AST 規約）も PASS
- [ ] **Step 6: 文書**: `USAGE.md` の `dataset_load` 行に「pipeline 経路では差次的結果も（データと群の指紋が一致するものだけ）戻す」を足す。`docs/workflow/mztab.md` の `dataset_load` 連鎖に `metabolomix/pipeline/session_handoff.py  load_run_dataset()` 以下の新しい呼び出し（足した関数名）を書く。
- [ ] **Step 7: コミット** `feat(pipeline): run の差次的結果を指紋の一致を条件に session へ戻す`

### Task 4: SME 候補ごとの InChIKey と、曖昧な SML 注釈の InChIKey 一覧

**Files:**
- Modify: `metabolomix/mztab/dataset_state.py`（`build_dataset_state` の SME 候補作成部分、`_build_feature_annotations`）
- Test: `tests/test_dataset_state.py`

**Interfaces:**
- Produces:
  - `ds.feature_candidates[fid][i]["inchikey"]: str | None` と `["inchikey_source"]: str | None`（`_resolve_inchikey` を候補ごとに、特徴単位とは**別の** counts dict で呼ぶ）。
  - 曖昧な SML 注釈 `{"ambiguous": True, "sml_ids": [...], "name": None, "inchikeys": [str, ...]}`（解決できたものだけ。重複なし。`inchikey` キーは今どおり持たない）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_each_sme_candidate_carries_its_own_inchikey(): ...
    # rank1 と rank2 が別の InChIKey 形 database_identifier → 候補ごとに別の inchikey
def test_candidate_inchikeys_do_not_change_feature_coverage_counts(): ...
    # inchikey_coverage（msp_lookup の計数を含む）が変更前と同じ値
def test_ambiguous_sml_annotation_keeps_the_inchikey_list_but_no_inchikey(): ...
    # "inchikeys" に 2 件、"inchikey" キーは無い
```

- [ ] **Step 2: FAIL を確認** → **Step 3: 実装** → **Step 4: `tests/test_dataset_state.py` と `tests/test_dataset_analysis_tools.py`（書き出し）が PASS**
- [ ] **Step 5: コミット** `feat(mztab): SME 候補ごとに InChIKey を解決し、曖昧な SML 注釈の InChIKey を残す`

### Task 5: 目標の正規化・理論 m/z・規則 hash（`feature_search` の土台）

**Files:**
- Create: `metabolomix/analysis/feature_search.py`
- Test: Create `tests/test_feature_search.py`

**Interfaces:**
- Consumes: `msdial/peak_verification.parse_formula` `monoisotopic_mass`、`msdial/adducts.parse_adduct` `mz_from_neutral`、`core/atomic_io.canonical_hash` `DomainError`。
- Produces:
  - `DEFAULT_ADDUCT_TIERS: dict[str, list[list[str]]]`（Global Constraints の値）、`MAX_TARGETS = 50`。
  - `theoretical_mz(formula: str, adduct: str) -> float | None`（読めなければ None）。
  - `normalize_targets(targets: list[dict], *, dataset_polarity: str | None) -> list[dict]`。各要素のキー: `target_id, name, inchikey_block, formula, expected_mz, adduct_tiers, theoretical_mz (dict[adduct, float]), ppm, rt_min, rt_tol_min, polarity, applicable (bool)`。検証失敗は `DomainError("FEATURE_TARGET_INVALID", message, {"target_id", "field"})`。
    - `inchikey` 必須（先頭 14 文字を大文字で `inchikey_block` に）。`formula` か（`expected_mz` と `adduct`）のどちらか必須。`expected_mz` のときは `adduct_tiers = [[adduct]]`。
    - `target_id` の重複、51 件以上、未知のキー、読めない付加体はエラー。
    - `polarity` が与えられ `dataset_polarity` と違えば `applicable=False`。
    - 極性が決まらず（`dataset_polarity is None` かつ目標の `polarity` も無い）、`adduct_tiers` も `expected_mz` も無い目標が 1 つでもあれば `DomainError("POLARITY_UNKNOWN", ..., {"target_ids": [...]})`（呼び出し全体のエラー）。
  - `rule_hash(normalized_targets: list[dict]) -> str`（`canonical_hash` に目標のリストを渡す。result_id・dataset_hash は入れない）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_theoretical_mz_of_gaba_protonated():
    assert theoretical_mz("C4H9NO2", "[M+H]+") == pytest.approx(104.0706, abs=1e-4)

def test_default_tiers_follow_dataset_polarity():
    [t] = normalize_targets([{"target_id": "c1", "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                              "formula": "C4H9NO2"}], dataset_polarity="negative")
    assert t["adduct_tiers"] == DEFAULT_ADDUCT_TIERS["negative"]
    assert t["inchikey_block"] == "BTCSSZJGUNDROE" and t["ppm"] == 10 and t["rt_tol_min"] == 0.3

def test_rule_hash_is_the_same_whether_defaults_are_explicit(): ...
def test_target_of_other_polarity_is_not_applicable(): ...
def test_missing_inchikey_or_mass_is_invalid(): ...           # FEATURE_TARGET_INVALID
def test_more_than_fifty_targets_is_invalid(): ...
def test_unknown_polarity_without_explicit_adducts_is_refused(): ...   # POLARITY_UNKNOWN
```

- [ ] **Step 2: FAIL** → **Step 3: 実装** → **Step 4: PASS**（`.venv/Scripts/python.exe -m pytest tests/test_feature_search.py -q`）
- [ ] **Step 5: コミット** `feat(feature_search): 目標の正規化・理論 m/z・規則 hash`

### Task 6: 目標照合の判定（spec §5）

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`
- Test: `tests/test_feature_search.py`

**Interfaces:**
- Consumes: Task 4 の候補ごとの `inchikey`、曖昧な SML の `inchikeys`。`feature_bindings._normalize_adduct`（import して使う。`bind_features` は変えない）。
- Produces:
  - `classify_identity(candidates: list[dict], annotation: dict | None, *, inchikey_block: str, adduct: str) -> tuple[str, list[str], int | None]` → `(evidence, reasons, sme_rank)`。`evidence` は `sme_match` / `sml_match` / `unannotated` / `contradicted`。SME 候補の InChIKey は**その候補自身の** `inchikey` だけを見る（特徴の rank1 の値を混ぜない）。電荷は `abs(candidate["charge"]) == parse_adduct(adduct).charge` で比べる（符号の書き方の違いに依存しない）。
  - `match_target(ds, target: dict, *, max_results: int) -> dict` → `{"target_id", "status", "reason", "selected_feature_id", "selected_adduct", "selected_evidence", "adduct_tier_used", "counts": {"mass", "rt", "qualified", "by_evidence": {"sme_match", "sml_match", "unannotated"}}, "theoretical_mz", "candidates": [row...], "truncated": int}`。`applicable=False` の目標は `status="not_applicable"`, `reason="polarity_mismatch"`, 候補なし。
  - 候補行のキー: `target_id feature_id adduct adduct_tier theoretical_mz mz ppm rt delta_rt evidence sme_rank annotated_name annotated_inchikey qualified reasons`（`reasons` は `;` 区切り）。並びは（qualified が先、`adduct_tier` 昇順、`abs(ppm)` 昇順）。
  - 判定（spec §5 の 4）: qualified が 1 件以上ある最初の層 → その層で qualified がある最初の証拠水準（`sme_match` → `sml_match` → `unannotated`）→ その水準の qualified を**特徴 ID で一意化して**数え、1 件なら `resolved`、2 件以上なら `needs_input`。使わなかった層・水準の qualified 候補の `reasons` に `lower_adduct_tier` / `lower_evidence_tier`。
  - `reasons` の語彙: `rt_outside_tolerance`、`identifier_mismatch`、`adduct_mismatch`、`charge_mismatch`、`sme_identity_unresolved`、`sml_ambiguous`、`mass_only_not_identified`、`lower_adduct_tier`、`lower_evidence_tier`。

- [ ] **Step 1: 失敗するテストを書く**（ds は `tests/test_feature_bindings.py` の `_dataset` / `_candidate` に倣った最小の合成。候補には Task 4 の `inchikey` を入れる）

```python
def test_main_ion_tier_wins_over_other_adducts(): ...        # [M+H]+ に 1 件、[M+Na]+ に 1 件 → resolved は [M+H]+、Na 側は lower_adduct_tier
def test_sme_match_beats_sml_and_unannotated(): ...
def test_two_qualified_in_best_tier_is_needs_input_even_if_one_is_far_more_intense(): ...
    # feature_matrix で片方を 1000 倍にしても status="needs_input", reason="multiple_qualified_candidates"
def test_no_candidate_is_needs_input_with_counts(): ...      # counts.mass == 0
def test_sme_of_another_compound_is_dropped(): ...           # evidence="contradicted", reasons に identifier_mismatch, qualified False
def test_sml_of_another_isomer_is_not_a_contradiction(): ...  # evidence="unannotated"
def test_unannotated_needs_rt(): ...                         # rt_min なし → mass_only_not_identified
def test_rank2_sme_match_is_found_and_rank1_inchikey_not_borrowed(): ...  # sme_rank == 2 / rank2 が別物質なら contradicted
def test_unresolved_sme_identity_is_not_a_contradiction(): ...  # 候補の inchikey None → sme_identity_unresolved、RT ありで unannotated として qualified
def test_ambiguous_sml_with_a_matching_member_is_sml_match(): ...  # reasons に sml_ambiguous
def test_same_feature_under_two_adducts_counts_once_per_tier(): ...
def test_candidates_are_truncated_with_count(): ...           # max_results=2 で truncated == 総数 - 2
```

- [ ] **Step 2: FAIL** → **Step 3: 実装** → **Step 4: PASS**、`tests/test_feature_bindings.py` も無変更で PASS
- [ ] **Step 5: コミット** `feat(feature_search): 目標照合の判定（付加体の層・証拠の水準・強度で選ばない）`

### Task 7: 窓検索と統計の結合（spec §4.2・§6）

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`
- Test: `tests/test_feature_search.py`

**Interfaces:**
- Consumes: `result_state.is_current`、Task 2 の分岐。`ds.results`、`ds.detected_mask`、`ds.sample_names`、`ds.feature_ids`。
- Produces:
  - `search_window(ds, *, feature_ids: list[str] | None, min_mz, max_mz, min_rt, max_rt, max_results: int) -> tuple[list[dict], int]`（行 `feature_id mz rt evidence name inchikey adduct n_sme`、m/z 昇順、打ち切った件数）。条件が 1 つも無ければ `DomainError("FEATURE_SEARCH_EMPTY_QUERY", ...)`。`min > max` は `DomainError("FEATURE_SEARCH_BAD_RANGE", ...)`。
  - `select_results(ds, result_ids: list[str] | None) -> dict` → `{"results": [result...], "stale": [id...], "missing": [id...]}`。省略時は `kind == "two_group"` で `is_current` のものを全部（古いものは `stale` に）。名指しのときは無い ID を `missing`、古い ID を `stale` に入れ、`results` には有効なものだけ。
  - `feature_stats(ds, feature_ids: list[str], results: list[dict]) -> list[dict]`（行 `feature_id result_id status log2fc p q mean_a mean_b detected_a detected_b`。`status` は `tested` / `not_in_matrix` / `not_testable`。後の 2 つは数値欄 None）。
  - `result_summaries(ds, results: list[dict]) -> list[dict]` → 各 `{"result_id", "group_a", "group_b", "n_a", "n_b", "log2fc_sign", "bh_population": {"n_features_total", "n_in_matrix", "n_tested"}, "caveats"}`。`n_in_matrix` は結果の `results` の件数、`n_tested` は p が有限の件数、`n_features_total` は `len(ds.feature_ids)`。
  - 検出例数: マスクの行は `ds.feature_ids` の位置、列は結果の `samples_a` / `samples_b` を `ds.sample_names` の位置に直す。マスクが無い、または列数が `len(ds.sample_names)` と違えば None（呼び出し側が注記）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_window_and_filters_are_anded_and_inclusive(): ...
def test_empty_query_is_refused(): ...
def test_unannotated_feature_gets_its_statistics(): ...      # inchikey 無しの特徴にも tested 行
def test_feature_dropped_by_preprocess_is_not_in_matrix(): ...
def test_null_p_from_persisted_result_is_not_testable(): ...  # results 行の p=None
def test_q_is_taken_from_the_result_not_recomputed(): ...
def test_bh_population_counts(): ...
def test_default_selection_skips_stale_and_non_two_group(): ...  # kind="welch" の stat_ は含まない
def test_named_stale_result_is_reported(): ...
def test_detection_counts_per_group(): ...
def test_detection_counts_skipped_when_mask_shape_mismatches(): ...
```

- [ ] **Step 2: FAIL** → **Step 3: 実装** → **Step 4: PASS**
- [ ] **Step 5: コミット** `feat(feature_search): 窓検索と、注釈の無い特徴を含む差次的統計の結合`

### Task 8: MCP ツール `dataset_feature_search` と登録・文書

**Files:**
- Modify: `metabolomix/tools/dataset_analysis_tools.py`（関数と `__all__`）
- Modify: `tests/test_server_registration.py`（`EXPECTED_TOOLS` に `dataset_feature_search`）
- Modify: `tests/test_workflow_docs.py`（`IN_SCOPE["dataset_analysis.md"]` に `dataset_feature_search`）
- Modify: `USAGE.md`（ツール表に 1 行、宣言件数を 75）、`CLAUDE.md`（冒頭「ツール 74」→「ツール 75」）、`docs/workflow/dataset_analysis.md`（`## dataset_feature_search` の呼び出し連鎖）、`docs/workflow/index.md`（dataset_analysis.md の件数 5 → 6）、`docs/output_format/mztab.md`（戻り値の欄の意味: `evidence` の値、`status` の 3 値、`bh_population`、`detected_*`、`needs_input` の読み方、`rule_hash` と `stats_source`）
- Test: Create `tests/test_dataset_feature_search_tool.py`

**Interfaces:**
- Consumes: Task 5〜7 のすべて。`console/output_collector.read_adduct_polarity(path)`、`console/job_manager.load_job(Path(ds.job_path)).polarity`。
- Produces: `dataset_feature_search(feature_ids=None, min_mz=None, max_mz=None, min_rt=None, max_rt=None, targets=None, result_ids=None, max_results=20) -> str`（spec §4 の外形と §4.2 の戻り値。`targets` と窓の同時指定は `bad_request` の封筒）。
  - 極性: ジョブの宣言極性があればそれ、無ければ `ds.source_files` の mzTab-M に `read_adduct_polarity` の `adduct_majority`。両方あって食い違えば宣言を使い `notes` に書く。`rule.polarity_source` は `job_declared` / `adduct_majority` / `unknown`。
  - `ds` 無し → `missing_state("dataset", ["dataset_load"])`。名指しの `result_ids` に `missing` / `stale` があれば `missing_state("dataset_differential_result", ["dataset_differential"])`。差次的結果が 1 件も無いときは統計なしで返し `notes` に書く。
  - キュレーションのフラグ（`curation_apply.arf2_for_mztab` → `flags_for_arf2` の `n`）があれば `notes` に「キュレーションを反映していない」。
  - 戻り値のキー: 窓検索 `features stats results stats_source notes`、目標照合 `targets candidates stats results rule stats_source notes`。`rule = {"targets": normalized, "polarity", "polarity_source", "rule_hash"}`、`stats_source = {"result_ids", "dataset_hash"}`。

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_dataset_analysis_tools.py` の `reset_session` と `_load_ds` に倣う）

```python
def test_without_dataset_returns_missing_state(): ...
def test_window_and_targets_together_is_bad_request(): ...
def test_target_mode_payload_shape_and_rule_hash(): ...      # rule.rule_hash が 64 桁、stats_source.result_ids、candidates は TSV
def test_named_stale_result_returns_missing_state(): ...
def test_no_differential_result_returns_without_stats_and_notes_it(): ...
def test_polarity_from_job_declaration_is_preferred(): ...
def test_payload_has_no_structured_duplicate(): ...           # 登録の structured_output=False を確認する既存流儀に倣う
```

- [ ] **Step 2: FAIL** → **Step 3: ツール実装と文書・登録の更新** → **Step 4: PASS**: `.venv/Scripts/python.exe -m pytest tests/test_dataset_feature_search_tool.py tests/test_server_registration.py tests/test_workflow_docs.py tests/test_readme_links.py -q`
- [ ] **Step 5: コミット** `feat(tools): dataset_feature_search（DatasetState の特徴検索と目的の物質の判定）`

### Task 9: 一連の流れの確認と記録

**Files:**
- Test: `tests/test_pipeline_differential_handoff.py`（1 本追加）
- Modify: `docs/HISTRY.md`（追記）、`docs/task.md`（完了項目を消す・残課題を足す）、vault `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`

**Interfaces:**
- Consumes: Task 3 と Task 8。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_feature_search_uses_pipeline_statistics_without_recomputing(tmp_path):
    # Task 3 と同じ手順で differential:c1 を持つ run → dataset_load(pipeline_path=...) →
    # dataset_feature_search(feature_ids=[...]) の stats の q が、保存した結果の q と一致し、
    # stats_source.result_ids が保存時の result_id。compare_dataset は呼ばれない（monkeypatch で検知）
```

- [ ] **Step 2: PASS を確認**（Task 3・8 が済んでいれば実装は要らない。FAIL なら原因の Task に戻る）
- [ ] **Step 3: 全テスト**: `.venv/Scripts/python.exe -m pytest tests -q` → 失敗 0（`tests/test_pathway_tools.py` `tests/test_pdf_raster.py` が PyMuPDF の DLL ブロックで収集エラーになる環境では、その 2 件だけが原因であることを出力で確かめて記録する）
- [ ] **Step 4: 記録**: HISTRY に実装の要点、task.md から 2026-10-19 節の「DatasetState の特徴検索ツール」項目を消す。vault の流れ図に (a) DatasetState 経路の解析枝 `dataset_feature_search`、(b) `dataset_load(pipeline_path=...)` が差次的結果も戻すことを足し、末尾の出典行の日付を直す。
- [ ] **Step 5: コミット** `test: run の統計を再計算せずに特徴検索へ渡す一連の流れ`
