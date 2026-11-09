# D「pipeline を上流に狭める」実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** pipeline を「生データ → MS-DIAL Console → 検証済みの解析フォルダ（と、その受領の記録）」に限り、下流（前処理・比較・PCA・書き出し・レポート）と v2 と v1 の統計の鎖を消す。Console の入口を GUI 解析済みフォルダと対等にそろえる。

**Architecture:** 工程は極性の単位ごとに `prepare_input:<極性>` → `upstream:<極性>` → `validate_outputs:<極性>`、最後に 1 つの `receipt`。要求は `pipeline-request.v3`（1 極性の今の形も受け、内部では `units` に正規化）、run は `pipeline-run.v3`（極性ごとの状態を `record["units"][<極性>]` に持つ）。既存の `pipeline-run.v1`（27 本）は読むだけ（`pipeline_status`・`load_dataset`）。`load_dataset(pipeline_path=…)` は受領（receipt）の指す解析フォルダを解析フォルダのアダプタで読み、試料の対応表とライブラリの記憶も受領から取る。

**Tech Stack:** Python 3.14、FastMCP、pytest（main ツリーの `.venv`）。

**Spec:** `docs/superpowers/specs/2026-11-06-pipeline-upstream-design.md`（全体 `2026-11-06-overhaul-overview-design.md` §2・§4）。

## 範囲

含める:
- spec §4〜§8: 工程の作り直し・両極性の run・既存の run の扱い・v2 の廃止・ツールの整理（旧 Console 6 つの削除）。
- task.md §2「D で消す v1 の関数」の削除（`tests/v1_chain.py` と、それに頼る 12 のテストファイルの書き直しを含む）。v1 の結果（`kind="two_group"`）と pipeline からの引き継ぎ（`handoff`）だけが通る消費側の分岐も消す（旧と新を並べたまま残さない。全体 §7）。
- Console を「補助」と書く箇所を「対等」へ: `MCP_INSTRUCTIONS` の ENTRY POINT、`USAGE.md` §5/§8 と `pipeline_run` の行、`CLAUDE.md` の (d)・(e) と pipeline の決まりごと。
- 文書: `docs/workflow/pipeline.md`（Console の 2 ツールの節を足す）・`metabolomics.md` の削除・`index.md`、`docs/output_format/dataset.md` §10.2.1・§11.2.2 ほか、`docs/schema/` の v2 の 4 文書、end-to-end spec §9 の状態の欄、vault の流れ図（経路 A・D と入口）。
- task.md §9 の振り分け（下の表）。

含めない（task.md に残す）:
- `LIPIDMIX_PIPELINE_INDEX_DIR`（`store` が読むテスト用の環境変数。B-2 の「資産の場所の環境変数」の対象外）。
- `reports/stats/` からの行列の復元（判断表 #10）、定量の `restore_chain`。
- Use-LLLM・massbank-context・reanalysis-study の追従（main の切り替えのとき）。

## Global Constraints

- ツールは **59 → 53**（`console_plan` `console_run` `console_status` `console_prepare_input` `console_cleanup` `job_list` を消す。足すものは無い）。リソース 2・リソーステンプレート 1 は変えない。ツール数を変えるコミットで `tests/test_server_registration.py`・`USAGE.md`（一覧と宣言件数）・`CLAUDE.md` 冒頭・`tests/test_workflow_docs.py` の `OUT_OF_SCOPE` を同時に直す。消した名前は `tests/test_removed_tools_not_mentioned.py` の `REMOVED_TOOLS` に足し、入口の文書と `metabolomix/` の文字列・コメントから消す（`console/` の docstring・エラーの `required_tools` を含む）。
- `export_contract.CONTRACT_VERSION = 1`・15 列・列順・既存のメタ行は変えない（massbank-context との契約）。根拠バンドルは v3 のまま。
- `metabolomix/pipeline/{engine,worker,service,recovery,store}.py`（と新しい `receipt.py`）はグローバル session を import しない（`tests/test_pipeline_engine.py` の AST の検査に新しいモジュールも足す）。
- run の `results` は追記専用（`store._assert_results_append_only`）。永続化する pipeline 内のパスは `pipeline_root` 相対。
- Console の自動再試行はしない。やり直しは `rerun_upstream=true` の明示だけで、新しい attempt ディレクトリを作る。
- `pipeline_status` は読取専用（生存確認の失敗を「死んでいる」と読まない）。
- `STATE_REVISION_CONFLICT` の有限回リトライの規則（`store._patch_request_id_with_retry`・`recovery.prepare_resume`・`engine._save_with_retry`）は保つ。
- 工程の実行順の正準は `engine.build_stages` の計画順。**受領（`receipt`）は必ず最後**（旧「report は必ず最後」を置き換える）。
- 研究室の参照ライブラリ（MSP・`[library] lbm` の LBM）は、受領・戻り値・エラーのどこでもファイル名と sha256 だけを示す（run の記録の `inputs` は絶対パスのまま。今の規則）。
- 全ツールに `structured_output=False`、戻り値は `json_payload`。
- テスト: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。コミットは pre-commit が全テストを回すのでバックグラウンドで実行し出力をファイルへ。フックの実行中はファイルを編集しない。`git stash` を使わない。マージも push もしない。
- 文書・docstring・エラーの文面は日本語。

## ツールの最終形

| ツール | 扱い | 引数・戻り値の変化 |
|---|---|---|
| `pipeline_plan` | 残す | 引数は同じ。`request` は v3（下）。戻り値の `resolved` は極性ごと: `resolved.units.<極性>` に今の欄（method / lbm / polarity / library_mode / msp / alignment_reference / dot_product_cutoffs）。`effective_target` は消す。 |
| `pipeline_run` | 残す | 同上。docstring を「対等な 2 つ目の入口。結果は解析フォルダで、続きは `load_dataset(pipeline_path=…)`」へ。 |
| `pipeline_status` | 残す | `stage_statuses` に加え `units`（極性ごとの `status`: `pending`/`running`/`validated`/`failed`）と、完了後は `receipt`（解析フォルダの場所の要約と `next`）。旧 v1 の run は今のまま読め、`legacy: true` と注意 1 行。 |
| `pipeline_resume` | 残す | `updates` は `sample_manifest`（1 極性）か `units: {<極性>: {"sample_manifest": …}}` だけ。下流のキー（`target` `preprocess` `comparisons`）は `PIPELINE_DOWNSTREAM_REMOVED`（`details.replacements`: `load_dataset` → `dataset_preprocess` / `dataset_compare` / `dataset_pca`）。旧 v1 の run は `PIPELINE_RUN_LEGACY`（読むなら `load_dataset(pipeline_path=…)`、Console をやり直すなら新しい `pipeline_run`）。`rerun_upstream=true` は upstream が成功していない極性だけをやり直す。 |
| `pipeline_cancel` | 残す | 旧 v1 の run は `PIPELINE_RUN_LEGACY`。 |
| `console_method_candidates` | 残す | 変えない。`docs/workflow/pipeline.md` に節を足す（`OUT_OF_SCOPE` から外して `IN_SCOPE["pipeline.md"]` へ）。 |
| `console_method_template` | 残す | 文面の `console_plan` を `pipeline_plan` へ。同上の節。 |
| `console_plan` `console_run` `console_status` `console_prepare_input` `console_cleanup` `job_list` | 消す | 計画・実行・状態は `pipeline_plan` / `pipeline_run` / `pipeline_status` へ吸収。`console_prepare_input` は要求の `keep_extension`（pipeline は `<run>/<極性>/input/` へ配置するので元フォルダの横に作らない）。`console_cleanup` は足さない（判断 7）。`job_list` は無し（既存のジョブは `load_dataset(job_path=…)` で読める）。 |
| `load_dataset` | 残す | 引数 `polarity: str \| None = None` を足す（`pipeline_path` のときだけ有効。2 極性の run で省くと `needs_input`、他の入口で渡すと `PIPELINE_POLARITY_NOT_APPLICABLE`）。`job_path` の経路は既存のジョブを読むために残す。 |

`session.current_job_path` と `console/worker.py`・`console/detached.py`（`console_run` の切り離し監視だけが使う。pipeline は `execution.supervise` を自分の worker から直接呼ぶ）、`store.register_job_owner` / `pipeline_owner_block_reason`（`console_run` / `console_cleanup` から run の job を守るためだけにある）も消す。

## 要求 `pipeline-request.v3`

- `schema` 省略・`"pipeline-request.v1"`・`"pipeline-request.v3"` を受け、保存するのは v3。`"pipeline-request.v2"` は `PIPELINE_REQUEST_INVALID`（`details.removed_schema`）。
- 共有のキー: `library_mode` `measure` `keep_extension` `timeout_s` `save_project` `output_root` `alignment_reference`（意味は今のまま）。
- 1 極性の形（今の形）: `polarity` `method_file` `lbm_file` `msp_file` `sample_manifest` をトップに置く。正規化で `units: [{"raw_dir": ".", …}]` の 1 単位になる。
- 2 極性の形: `units: [{"polarity", "raw_dir", "method_file"?, "lbm_file"?, "msp_file"?, "sample_manifest"?, "keep_extension"?}]`（1〜2 件。`polarity` と `raw_dir` は必須、`raw_dir` は `dataset_root` 相対、極性は重ならない）。`units` があるときにトップへ単位のキーを置くと `PIPELINE_REQUEST_INVALID`。
- 単位の `sample_manifest` の既定は `<raw_dir>/sample-manifest.tsv`（1 極性の形では今と同じ `dataset_root` 直下）。`analysis-request.json` の読み方（明示 > ファイル > 既定）は今のまま。
- 消すキー: `target` `preprocess` `comparisons`（と内部の `effective_target`）。明示した要求にあれば `PIPELINE_DOWNSTREAM_REMOVED`、`analysis-request.json` にあれば読み飛ばして警告 `REQUEST_FILE_DOWNSTREAM_IGNORED`（判断 5）。
- 受付の指紋（`request_fingerprint`）は v3 の正規形の上で計算する。旧 v1 の run とは指紋が別になり、再利用されない（同じ `request_id` を使い回すと `IDEMPOTENCY_CONFLICT`）。

## run `pipeline-run.v3`

```
{schema: "pipeline-run.v3", identity, request: {revision, request_id, content_hash, saved_path},
 status, state_revision,
 stages: {"prepare_input:<p>", "upstream:<p>", "validate_outputs:<p>", …, "receipt"},
 units: {<p>: {"raw_dir", "unit_root": "<p>", "inputs": <入力計画→配置後の snapshot>,
               "upstream": {console_job_path, execution_id, verification}}},
 inputs: {"fingerprint"}, results, needs_input, warnings, worker}
```

- 置き場所: `<pipeline_root>/<p>/input/`（raw の配置。MS-DIAL が `.arf` 等を書く＝解析フォルダ）、`<p>/inputs/`（実効メソッド・注釈器の表）、`<p>/console/attempt-NNNN/`（job と mzTab-M）、`results/receipt.json`。単位の中の相対パスは `unit_root` 基準。
- 単位の順は要求の順。単位の失敗（`failed`）は、その単位の残りの工程を `pending` のまま止め、他の単位へ進む（spec §5「片方が失敗しても、もう片方の結果は残す」）。`receipt` は検証を通った単位が 1 つでもあれば走る。
- 最終の status: 全単位が validated かつ receipt 成功 → `completed`、validated が 1 つ以上で欠けあり → `partial`、0 → `failed`、取消 → `cancelled`。worker の中の `needs_input` は無くなる（受付の事前検査だけが `needs_input` を保存する）。
- `SUPPORTED_SCHEMAS = {"pipeline-run.v1", "pipeline-run.v3"}`。v1 は読むだけ（`save_run`・engine・`prepare_resume`・`request_cancel` は `PIPELINE_RUN_LEGACY`）。`pipeline-run.v2` は読まない（実データに 0 本。読むと `PIPELINE_RUN_INVALID` の今の扱い）。

### 受領 `receipt`（`results/receipt.json`、`lipidmix.pipeline-receipt.v1`）

```
{schema, pipeline_id, created_at, request_revision,
 units: {<p>: {"status": "validated" | "failed" | "not_run",
   "raw_dir", "analysis_dir": "<p>/input", "console_job": "<p>/console/attempt-NNNN/analysis-job.json",
   "mztab": <pipeline_root 相対>, "alignment": {"arf": <名前>, "arf2": <名前>} | null,
   "method": {"file", "sha256", "effective_sha256"},
   "library": {"mode", "msp": {"file", "sha256", "source"} | null, "lbm": {"file", "sha256", "source"} | null},
   "console": {"exe_file", "exe_sha256", "software_version"},
   "alignment_reference": {…今の欄…},
   "samples": {"manifest": {"path", "sha256"} | null, "rows": [parse_manifest の行] | null, "n_raw"},
   "validation": {"job_status", "termination", "exit_code", "gui_project", "warnings"},
   "error": {code, message} | null}},
 next: "load_dataset(pipeline_path=…, polarity=…)"}
```

- `alignment` は `validate_outputs` が、正準 mzTab-M と照合して `.arf` がちょうど 1 本あることを確かめた結果（`mztab.evidence.matching_arfs`、`dataset/opening._matching_arf` と同じ規則）。合わなければ単位は `validated` のまま `alignment: null` と警告 `ANALYSIS_FOLDER_UNMATCHED`、run は `partial`（判断 6）。
- `samples.rows` は受領の時点の実験情報シートを `parse_manifest(…, expected_sources=job の raw の一覧)` で解した行。読めなければ `rows: null` と警告 `SAMPLE_MANIFEST_INVALID`（事前検査を通った後で書き換えられたとき）。

## 工程の新しい列

- `engine.build_stages(request)` が唯一の正準: 単位ごとに `{"stage_id": f"{h}:{p}", "handler": h, "unit": p}`（`h` は `prepare_input` / `upstream` / `validate_outputs`）、最後に `{"stage_id": "receipt", "handler": "receipt", "unit": None}`。`store` が `record["stages"]` を作るときも同じ関数を呼ぶ（`stage_plan.py` は消し、`BASE_V2_STAGE_IDS` も無くなる）。
- `_TRUST_PERSISTED`（成功なら信頼して飛ばす）は handler キーで `prepare_input` / `upstream` / `validate_outputs`。`receipt` は `result_refs` を実ファイルと照合（今の 3 つ目の分類）。`_ALWAYS_RECONSTRUCT` と `_NEEDS_INPUT_CODES` は消す（runtime を組み立てる工程が無くなる）。
- `make_context` は `unit`・`unit_root`・`unit_state`（`record["units"][p]`）を渡す。handler の `record_updates` に加え `unit_updates`（`record["units"][p]` へのキーごとの置き換え）を `commit_stage_outcome` が反映する。
- 再開（`recovery.prepare_resume`）: `rerun_upstream=true` は upstream が `succeeded` でない単位の 3 工程と `receipt` を `pending` に戻す。`sample_manifest` の更新・シートの中身の変化は `receipt` だけを戻す（アライメントの参照ファイルは計画の時点で固定済み。今と同じ）。どれかを戻したら `receipt` も戻す。

## v1 の関数と v2 の削除（一覧）

v2（Task 1）:
- `metabolomix/pipeline/{request_v2,stage_plan,metabolomics_handlers}.py`
- `metabolomix/console/{profile_schema,profiles,profile_adapter}.py`（共有の部品は無い。spec の「`hash_files` は v1 も使う」は誤りだった。調べた）
- `metabolomix/analysis/{assay_evidence,assay_qc,feature_bindings,internal_standards,statistics_v2,matrix_state,multigroup,feature_export}.py`。`feature_bindings._normalize_adduct` と `_ADDUCT_UNIT_CHARGE_RE` は `feature_search.py` へ移す（生きた利用者はこれだけ）。
- 道連れに死ぬもの: `inputs.plan_from_profile`・`inputs.resolve_raw_inventory`（と `test_reanalysis_small_fixes::test_raw_inventory_excludes_cache`）・`method_file.METABOLOMICS_REFERENCE_KEYS`・`sample_manifest.validate_standard_assays` / `validate_independent_samples`・`preprocess_policy.explicit_policy`・`differential.welch_t`・`arf.reader.readable_block_count`・`report.required_outputs_v2` / `analysis_status`・`request._peek_schema_declared_by_request_file`・`DatasetState.analysis_matrices` とその読み手（`result_state` の `standard_assays` / `bindings` / `matrix_recipe` の印、`core/handles.py`・`quant/derive.py`・`tools/mztab_tools.py`）・`analysis-job.v3`（`profile_snapshot`。書くのは v2 だけで実物は 0。`handoff/schema.py`・`mztab/loading._verify_effective_method`・`opening.py`・`library_memory.py` の v3 分岐）。
- テスト（丸ごと）: `metabolomics_fixtures.py`、`test_assay_evidence` `test_assay_qc` `test_feature_bindings` `test_internal_standards` `test_lcms_profile_inputs` `test_lcms_profile_schema` `test_metabolomics_{end_to_end,engine,profile_boundary,public_entry,report,request,stages,tools}` `test_multigroup` `test_statistics_v2` `test_analysis_matrix` `test_pipeline_inputs_v2` `test_pipeline_matrix_handoff`。中の v1 の回帰（`test_v1_request_fingerprint_is_unchanged` の指紋の固定・`test_load_run_still_reads_an_existing_v1_record`・`test_load_run_rejects_an_unknown_schema`・`test_v1_export_meta_is_unchanged_without_an_effect_size_definition`）は `test_pipeline_store.py` / `test_dataset_export*.py` へ移す。`test_pipeline_differential_handoff.py` は `test_pipeline_matrix_handoff._completed_run` に頼るので、Task 1 で v1 の run の 1 件と拒否理由の表のテストだけを残し、Task 5 で丸ごと消す。
- 文書: `docs/workflow/metabolomics.md`（`IN_SCOPE` から外す）・`docs/schema/{pipeline-request-v2,pipeline-run-v2,lcms-profile-v1,analysis-job-v3}.md`・`docs/workflow/index.md` の行・`docs/workflow/dataset.md` の `feature_bindings` / `matrix_state` / `analysis_matrices`・`docs/output_format/{mztab,dataset}.md` の同じ語・CLAUDE.md の (e)。

v1 の統計の鎖（Task 5。pipeline の下流が Task 3 で消えた後）:
- `analysis/dataset_analysis.py` の `run_dataset_preprocess` `run_dataset_pca` `run_dataset_differential` と私的な補助（`included_sample_indices`・`build_dataset_pp_inputs`・`PreconditionError` は残す）。
- `analysis/dataset_service.py` の `preprocess_dataset` `preprocess_auto` `pca_dataset` `compare_dataset` `resolve_comparison` `run_comparison` と補助（`apply_sample_manifest` と `_raw_manifest_layout` は `dataset_set_sample_metadata` が使うので残す）。
- `analysis/differential.py` の `two_group_test` `add_fdr` `summarize_two_group` `check_confounding` `one_way_anova` と内部の分布関数（`volcano_data` と `bh_fdr` は残す。`bh_fdr` の利用者が残るかは実装で確かめ、無ければ消す）。
- `analysis/pca.py` の `run_pca`（CLI の `--pca` も消す。`docs/cli.md` を直す）、`analysis/preprocessing.py` の `preprocess`。
- v1 の結果だけが通る消費側: `stats_results.as_two_group` の `kind == "two_group"` の素通し、`result_state.is_current` の handoff と `parent_ids` の分岐（`is_stats_current` と `species_pca` は残す）・`preprocess_fingerprint` ほか `pp_*` の印、`dataset_export._preprocess_parameters` の handoff 分岐と `HANDOFF_PREPROCESS_UNKNOWN`、`feature_search` の `_result_recipe`・`select_results` の `two_group` 分岐・`retest_merged` の v1 経路、`evidence_export` の v1 の `intensities.tsv` 経路、`export_contract.row_significant` の旗の無い行の後退、`tools/dataset_analysis_tools.py` の `"two_group"`、`tools/mztab_tools._handoff_lines`・`DatasetState.handoff_report`、`core/tool_helpers.dataset_pca_plot` の v1 の分岐、`pipeline/session_handoff.py` の `restore_differentials` 一式と `load_run_dataset` の行列・binding の読み戻し。
- `tests/v1_chain.py` を消し、使っている 12 ファイル（`test_compare_consumers` `test_curation_on_dataset` `test_dataset_analysis_tools` `test_dataset_feature_search_tool` `test_dataset_handle_tools` `test_dataset_plot` `test_evidence_export_links` `test_evidence_export_v3` `test_feature_search_labels` `test_folder_tool_parity` `test_reanalysis_tools` と、Task 5 で消す `test_pipeline_differential_handoff`）と `tests/evidence_bundle_fixture.py` を、統計エンジンの公開の入口（`dataset_preprocess` / `dataset_compare` / `dataset_pca` か `stats_results.register_*`）で状態を作る形へ書き直す。v1 の結果を前提にした断定（v1 の書き出しのゴールデン `tests/fixtures/v1_export_golden.tsv` など）は消す。

## task.md §9 の振り分け

| 項目 | 振り分け |
|---|---|
| `msp_only` の初回実データで Console のログに MSP のパスが出ないか | 残す（pipeline の `MSDIAL_EXECUTION_FAILED` / `CONSOLE_WARNING` の封筒として） |
| `report._section_method` の `lbm_path: (unknown)` | 不要（レポートを消す。受領は `library.lbm: null`） |
| user_config の後続: `LBM_NOT_FOUND` の details / pipeline の `MSDIAL_EXE_NOT_FOUND` の human_action / `_configured_exe` の引き直し / fake Console での受付テスト / `console_run` の `CONFIG_INVALID` | `console_run` の件は不要。残り 4 つは残す |
| `is_console_exe` の起動不能の分類（呼び出し 3 か所） | 残す（`console_tools.py` の 1 か所が消え、`pipeline/inputs.py` の 2 か所） |
| raw 一覧の共通層で実験情報シートの `source_file` に絞る（極性の混在） | 吸収（2 極性の単位は `raw_dir` ごとに既定のシートを持つ。1 つのシートで混在させる使い方は残す: 残す） |
| Console 停止時の理由消失（`JOB_TERMINATE_TIMEOUT` → `worker_lost`） | 残す（`execution.supervise` の中。pipeline でも起きる） |
| `.mcp.json` の per-server timeout と `console_run(detach=True)` | 不要（`console_run` を消す。pipeline は常に切り離し） |
| テキスト DB へ SMILES 列 | 残す |
| Console 4 ツールの `docs/workflow/` 昇格・end-to-end spec の昇格と §9 | 吸収（残る 2 ツールを `pipeline.md` へ。end-to-end spec は昇格せず §9 の状態の欄に「上流は D、下流は B/C/A で置き換えた」） |
| 上流 `lipidome_merge` と「1 解析・1 極性」の案 A/B/C | 吸収（D の両極性の run ＝ 1 run・極性ごとの解析フォルダ。統合は B-3 の `dataset_quant_merge`。案 B 相当に決まった） |
| 解析モードの設計軸を acquisition で | 残す |
| HOLD `analysis-job.json` の `dataset_root` 絶対パス依存 | 残す（`load_dataset(job_path=…)` が残るため） |
| HOLD メソッド候補 spec §9 の残り | 残す（`console_method_candidates` が残るため） |
| Console 層の設計指摘 8 件 | 1・3（`console_plan` / `_resolve_job_path`）は不要、2（`CONSOLE_ERROR_CODES`）・4・5（`load_dataset`）・6（`_job_id` の秒精度。pipeline は attempt で分けるので不要）・7（`output_collector` のサイズのみ）・8（`test_console_runner` の差し替え） → 2・4・5・7・8 を残す、6 は不要 |
| DatasetState 解析層の Minor（`console_tools.py` の bare except・`_OPERATIONAL_FILES` の basename・重複表示名） | 1 つ目は不要、残り 2 つは残す |
| v2 の R3〜R5・独立プロセス試験・Stage B・独立レビュー・`bind_features` の癖・HOLD R2/R5 | 不要（v2 を消す） |
| session_handoff の結果の全履歴の検証 | 不要（引き継ぎを消す） |
| `dataset_feature_search` の実データ確認（MTBLS6039） | 残す（文面を「`dataset_compare` で計算し直してから」へ） |

## タスクの分け方

依存: Task 1 → 2 → 3 → 4 → 5 → 6（3 と 4 は同じ層を触るので直列。どのコミットでも全テスト緑）。

### Task 1: v2 を消す
**Files:** 上の「v2」の一覧。`metabolomix/pipeline/{request,store,engine,recovery,service,report,session_handoff,inputs}.py`、`metabolomix/analysis/feature_search.py`、`metabolomix/handoff/schema.py`、`metabolomix/mztab/{loading,dataset_state}.py`、`metabolomix/dataset/{opening,library_memory}.py`、`metabolomix/core/handles.py`、`metabolomix/analysis/result_state.py`、`metabolomix/quant/derive.py`、`metabolomix/tools/mztab_tools.py`、文書。
- [ ] 移す v1 の回帰テストを先に移し、v2 のテストを消す。
- [ ] `_normalize_adduct` を `feature_search.py` へ移し、v2 のモジュールと分岐を消す（`build_handlers` の `by_schema` をやめ v1 の handler を直に）。
- [ ] `grep -rn "pipeline-request.v2\|pipeline-run.v2\|lcms-profile\|stage_plan\|metabolomics_handlers\|statistics_v2\|matrix_state\|feature_binding\|analysis_matrices\|profile_snapshot" metabolomix tests docs --include=*.py --include=*.md`（`docs/superpowers` を除く）が 0 件。
- [ ] 全テスト緑 → コミット `refactor(pipeline): v2（メタボロミクスの profile・statistics 工程）を消す`。

### Task 2: 旧 Console 6 ツールを消す（59 → 53）
**Files:** `metabolomix/tools/console_tools.py`（6 ツールとその私的な補助）、`metabolomix/console/{worker,detached}.py`（削除）、`metabolomix/pipeline/store.py`（`register_job_owner` / `pipeline_owner_block_reason`）、`metabolomix/core/session_state.py`（`current_job_path`）、`metabolomix/core/mcp_errors.py` と `console/*.py`・`mztab/loading.py`・`tools/{maintenance,mztab_tools}.py` の文面、`server.py`、テスト（`test_console_{cleanup,detached,payload_trim,plan_method,prepare_input_tool,provenance}` を消す。`test_console_runner` `test_console_config` `test_msdial_exe_guidance` `test_mztab_tools` `test_folder_tool_parity` は該当テストだけ消すか `pipeline_plan` / `console_method_template` へ）、`USAGE.md` §5、`CLAUDE.md` 冒頭、`tests/{test_server_registration,test_removed_tools_not_mentioned,test_workflow_docs}.py`。
- [ ] `MIXED_RAW_FORMATS` の `required_tools` を `["pipeline_plan"]` に。`console_method_template` の docstring の `console_plan` を `pipeline_plan` に。
- [ ] `test_removed_tools_not_mentioned` が緑（6 つの名前が `metabolomix/` と入口の文書に単語として残らない）。
- [ ] 全テスト緑 → コミット `feat(console): 旧 Console 6 ツールを消し pipeline へ吸収（53 ツール）`。

### Task 3: pipeline を上流に狭める（v3・単位・受領）
**Files:** `metabolomix/pipeline/{request,store,engine,recovery,service,inputs}.py`、新規 `metabolomix/pipeline/receipt.py`（受領の組み立てと `persist_result` など `report.py` から残す 3 関数）、`metabolomix/pipeline/report.py`（削除）、`metabolomix/tools/pipeline_tools.py`、テスト。
**Interfaces（Produces）:**
- `request.resolve_request(source_root, explicit) -> dict`（v3 正規形。`units` を必ず持つ）、`request.merge_updates(request, updates) -> dict`。
- `engine.build_stages(request) -> list[{"stage_id", "handler", "unit"}]`、`store.stage_ids_for(request) -> list[str]`。
- `store.run_units(record) -> dict[str, dict]`（v3 は `record["units"]`、v1 は `{<極性>: {"raw_dir": ".", "unit_root": ".", "inputs": record["inputs"], "upstream": record["upstream"]}}` の読み替え。読む側はすべてこれを通す）、`store.is_legacy(record) -> bool`。
- `receipt.build_receipt(record, pipeline_root) -> dict`、`receipt.read_receipt(pipeline_root, record) -> dict | None`。
- [ ] テスト（失敗するものから）: `test_pipeline_request.py` — `test_single_polarity_form_normalizes_to_one_unit`、`test_two_units_require_distinct_polarity_and_raw_dir`、`test_downstream_keys_are_rejected_with_replacements`（`PIPELINE_DOWNSTREAM_REMOVED`、`details.replacements` に `dataset_compare`）、`test_request_file_downstream_keys_are_ignored_with_warning`、`test_v2_schema_is_rejected`。`test_pipeline_engine.py` — `test_stage_plan_is_per_unit_then_receipt`、`test_failed_unit_does_not_stop_the_other_unit`（偽 handler で neg の upstream を failed、pos は validated、status `partial`、receipt は走る）、`test_all_units_failed_is_failed_without_receipt`。`test_pipeline_recovery.py` — `test_rerun_upstream_resets_only_failed_units`、`test_manifest_update_resets_only_receipt`、`test_legacy_run_resume_is_refused`（`PIPELINE_RUN_LEGACY`）。`test_pipeline_store.py` — `test_legacy_v1_record_is_read_but_not_saved`、`test_run_units_reads_legacy_flat_record`。`test_pipeline_end_to_end.py`（偽 Console、`tests/msdial_folder_fixtures.write_msdial_folder` で照合できる `.arf` も書かせる）— `test_single_polarity_run_ends_with_receipt`、`test_two_polarity_run_has_two_analysis_folders`、`test_unmatched_arf_is_partial_with_warning`。
- [ ] 実装。下流の handler（`load_dataset`〜`report`）・`_as_needs_input`・`_PRECONDITION_CODES`・`report.py` を消し、`finish_success` は単位の状態から決める。`pipeline_status` に `units` と `receipt` の要約。
- [ ] `tests/test_pipeline_engine.py` の AST の検査に `receipt.py` を足す。
- [ ] 全テスト緑 → コミット `feat(pipeline): 上流に狭める — 極性の単位・受領（pipeline-request.v3 / pipeline-run.v3）`。

### Task 4: `load_dataset(pipeline_path=…, polarity=…)` が受領を読む
**Files:** `metabolomix/dataset/{opening,library_memory}.py`、`metabolomix/pipeline/session_handoff.py`、`metabolomix/tools/handle_tools.py`、`metabolomix/core/mcp_errors.py`（`needs_input` を使う）、テスト `test_dataset_opening.py` `test_load_dataset.py` `test_library_from_record.py`。
**Interfaces（Consumes）:** Task 3 の `store.run_units`・`store.is_legacy`・`receipt.read_receipt`。
- [ ] テスト: `test_v3_run_opens_the_receipt_analysis_folder_with_manifest_rows`、`test_two_polarity_run_without_polarity_needs_input`（`{"needs_input": {"field": "polarity", "argument": "polarity", "retry_tools": ["load_dataset"], "options": [...]}}`）、`test_unit_not_validated_is_refused`（`PIPELINE_UNIT_NOT_VALIDATED`）、`test_polarity_with_directory_is_refused`、`test_legacy_v1_run_still_opens_with_its_manifest`（既存の `_pipeline` の手作りの記録で）、`test_v3_library_record_comes_from_unit_inputs`、`test_v3_run_mztab_source_has_no_restored_results`。
- [ ] 実装。`pipeline_arf2_candidates(pipeline_path, polarity=None)` は単位の `input/` を見る。基のキャッシュの指紋（`_dependencies`）に単位の `input/` と `results/receipt.json`。
- [ ] 全テスト緑 → コミット `feat(load_dataset): pipeline の受領から解析フォルダを開く（polarity）`。

### Task 5: v1 の統計の鎖と引き継ぎを消す
**Files:** 上の「v1 の統計の鎖」の一覧と 12 のテストファイル・`tests/evidence_bundle_fixture.py`・`docs/cli.md`。
- [ ] 12 のテストファイルを統計エンジンで状態を作る形へ書き直す（1 ファイルずつ、書き直した直後にそのファイルを回して緑を確かめる）。
- [ ] v1 の関数と消費側の分岐を消す。`grep -rn "two_group_test\|compare_dataset\|preprocess_dataset\|pca_dataset\|run_dataset_preprocess\|run_dataset_differential\|restore_differentials\|handoff_report\|v1_chain" metabolomix tests` が 0 件。`"two_group"` の残りは `as_two_group` の写しの欄だけ。
- [ ] 全テスト緑 → コミット `refactor(analysis): v1 の統計の鎖と pipeline からの引き継ぎを消す`。

### Task 6: 文書・入口の文言・記録
**Files:** `metabolomix/core/mcp_core.py`（ENTRY POINT）、`metabolomix/tools/pipeline_tools.py` の docstring、`USAGE.md` §5・§8（と §1 の `load_dataset` の `polarity`）、`CLAUDE.md`（冒頭・(d)(e)・構成の `pipeline/` の行・pipeline の決まりごと・ドキュメントの地図の end-to-end の行）、`docs/workflow/{pipeline,index,dataset}.md`、`docs/output_format/{dataset,mztab,reanalysis}.md`（§10.2.1・§11.2.2 の v1 の文面を消す）、`docs/superpowers/specs/2026-09-03-end-to-end-pipeline-design.md` §9、`tests/test_workflow_docs.py`（`console_method_*` を `IN_SCOPE["pipeline.md"]` へ）。リポジトリ外: vault の流れ図、main ツリーの `docs/HISTRY.md`（追記）・`docs/task.md`（末尾に追記だけ）。
- [ ] ENTRY POINT の 2 を「Console の経路は対等な入口。作るのは解析フォルダで、`pipeline_status` が completed / partial になったら `load_dataset(pipeline_path=…, polarity=…)` で 1 の鎖へ合流する」へ。`tests/test_instructions_text.py` が緑。
- [ ] 全テスト緑 → コミット `docs: Console の入口を対等に・pipeline を上流だけの文書へ`。
- [ ] vault の流れ図: 入口（0）、経路 A（A-2 の工程列・A-3 の resume のキー・A-4 の削除）、経路 D（手で回す経路を消し「メソッドを選ぶ・作る」だけに）、早見表、冒頭の規模 53、出典行の日付。
- [ ] HISTRY に節を足し、task.md の末尾に D の残りと §9 の振り分けの結果を足す。

## Review Focus

1. **2 極性の run で片方の Console が失敗する**: もう片方の解析フォルダは受領に載り `load_dataset(…, polarity=<成功した側>)` で開け、run は `partial`。→ Task 3 `test_failed_unit_does_not_stop_the_other_unit`・Task 4 `test_unit_not_validated_is_refused`。
2. **旧 v1 の run（27 本）を開く・再開する**: `load_dataset(pipeline_path=…)` は今までどおり試料の対応表つきで開け、`pipeline_resume` は代わりの道を返して止まり、記録を書き換えない。→ Task 3 `test_legacy_run_resume_is_refused`・Task 4 `test_legacy_v1_run_still_opens_with_its_manifest`。
3. **reanalysis-study の古い要求（`comparisons` や `target` を含む `analysis-request.json`）**: Console は回り、下流のキーは警告で読み飛ばす。→ Task 3 `test_request_file_downstream_keys_are_ignored_with_warning`。
4. **再開を何度も呼ぶ**: Console は `rerun_upstream=true` の明示なしに再起動しない（単位が増えても）。→ Task 3 `test_rerun_upstream_resets_only_failed_units` に「`rerun_upstream=false` では upstream を戻さない」の断定を足す。
5. **受領が研究室のライブラリの置き場所を漏らす**: `receipt.json` と `pipeline_status` の戻り値に MSP / `[library] lbm` の絶対パスが無い。→ Task 3 `test_receipt_names_library_files_only`。

## 判断が要る点（推奨つき）

1. **要求と run の版を上げる（`pipeline-request.v3` / `pipeline-run.v3`）**。推奨: 上げる。下流のキーが消え、極性ごとの状態を持つ形に変わるので、旧 v1 の記録と見分けられる版が要る。v2 は実データ 0 本なので読まず、名前の混同を避けて v3 にする。
2. **2 極性の要求の形**。推奨: `units: [{polarity, raw_dir, …}]`（1 極性は今の形を受けて 1 単位に正規化）。
3. **片方の極性が失敗したとき**。推奨: もう片方を続け、run は `partial`。極性は要求の順に 1 つずつ回す（同時に回すと 1 台の PC で Console が 2 本並ぶ）。
4. **既存の v1 の run は読むだけにする**。推奨: `pipeline_status`・`load_dataset` は今のまま読め、`pipeline_resume` / `pipeline_cancel` は `PIPELINE_RUN_LEGACY` で代わりの道を返す（spec §6 は下流の訂正だけを求めるが、上流の再開まで旧形式で支えると engine が 2 つの形を持ち続ける。27 本は完了済みで、Console をやり直したいなら新しい run で足りる）。
5. **下流のキー**。推奨: 明示した要求は `PIPELINE_DOWNSTREAM_REMOVED` で止め、`analysis-request.json` のものは警告で読み飛ばす（明示した比較が黙って計算されないのは困るが、置いてある古いファイルで Console まで止めるのは過剰）。
6. **`validate_outputs` で解析フォルダを確かめる**（mzTab-M と照合した `.arf` がちょうど 1 本）。推奨: 確かめる。合わなければ失敗にはせず警告と `partial`（`source="mztab"` でなら読めるため）。
7. **`console_cleanup` を pipeline の操作として残すか**。推奨: 残さない。pipeline は raw を `<run>/<極性>/input/` に置いて Console を回すので、元のフォルダに生成物が積もらない（`console_cleanup` が解いていた問題が起きない）。消したい run はフォルダごと消せばよい。
8. **`load_dataset` に `polarity` を足す**。推奨: 足す（2 極性の run はハンドルが 2 つになる。省いたら `needs_input` で選ばせる）。
9. **v1 の統計の鎖に頼るテスト 12 ファイルの扱い**。推奨: 統計エンジンで状態を作る形へ書き直す（消すと書き出し・図・特徴検索・根拠バンドルの消費側のテストが薄くなる）。Task 5 が最も大きくなる。
10. **`analysis-job.v3` と `DatasetState.analysis_matrices` も消す**。推奨: 消す（書くのは v2 だけで実物が無い）。
11. **Console の 2 ツールの文書**。推奨: `docs/workflow/pipeline.md` に節を足し、end-to-end spec は `docs/workflow/Lipidmix/` へ昇格させず §9 に「上流は D、下流は B/C/A で置き換えた」と書く。
