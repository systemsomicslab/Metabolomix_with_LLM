# 再解析 v8 の道具（アライメントの品質・LBM・起動検査・閾値）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** reanalysis-study の手順書 v8 が要る ms-data-parser の道具を足し・直す（割れたスポットの検出とまとめ、検出状態の取り込み、
`msp_then_lbm`、参照ファイルの選び方、起動検査、質量精度の絞り込み、閾値）。

**Architecture:** 純粋な計算は `metabolomix/analysis/` の新しいモジュール（`alignment_qc.py`）と既存モジュールに置き、ツールは
`tools/reanalysis_tools.py` の既存の型（ハンドル・`json_payload`・`mztab_error`）で包む。pipeline の変更は `pipeline/request.py`・
`pipeline/inputs.py`・`console/method_file.py`・`pipeline/service.py` の既存の経路に足す。

**Tech Stack:** Python（main ツリーの `.venv`）、pytest、FastMCP。MS-DIAL Console は現用の Debug ビルド（2026-09-03、ソース 45a531c4d）のまま。

**Spec:** `C:\Users\yuu18\reanalysis-study\.claude\worktrees\protocol-v8\docs\specs\2026-10-28-protocol-v8-design.md`
（4.1 (b)(c)、4.2、4.3、4.4、4.5 のハンドル、4.9）。対になる計画: reanalysis-study の `docs/plans/2026-10-28-protocol-v8.md`。

## Global Constraints

- 作業場所は worktree `.worktrees/feat-reanalysis-v8`（ブランチ `feat/reanalysis-v8`）。テストは
  `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q` を worktree のルートで回す。pre-commit が全数を回すので、
  コミットはバックグラウンドで出力をファイルに書いて待つ。マージ・push はしない。`git stash` は使わない。
- MS-DIAL Console を作り直さない。現用ビルドが読まないキー（`Number of threads`・`LBM annotator priority`）をメソッドに書かない。
- ライブラリ（MSP・LBM）の置き場所をエラー文・receipt・details に出さない（ファイル名と sha256 だけ。既存の `_resolve_msp_pinned` の規約）。
- 既存の要求の hash（`_FINGERPRINT_OMITTED_DEFAULTS`）と、完了済み run の再利用を壊さない。新しい要求キーの既定値は hash から外す。
- 割れのまとめの規則（設計 4.3 (c)）: 型が `injection_block` か `interleaved`、かつ 2 スポットの注釈が両立する組だけをまとめ、試料ごとに
  2 スポットの大きい方の値を使う。`group_aligned` はまとめずに印だけ付ける。
- 新しいツールを足したら、ツール数を数えるテスト・文書（`USAGE.md` 見出し、`docs/workflow/index.md`、`CLAUDE.md` の「ツール 90」、
  `docs/workflow/reanalysis.md`、`docs/output_format/reanalysis.md`）を同じタスクで直す。
- 要求の既定値や段の流れを変えたタスクは、Vault の `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md` も直す。

## Review Focus

1. **完了済みの run を開き直したときの統計の消失**: 検出状態が付くと `dataset_fingerprint` が変わり、v7 で作った差次的結果が
   `dataset_mismatch` で捨てられる。MTBLS6965 のハンドルで `dataset_feature_search` の stats が出続けること（Task 2 のテスト）。
2. **注入順の無い実験情報シート・シートが無い入力**: 参照ファイルはファイル名の順の中央になり、`injection_order_source: filename` が残る（Task 3）。
3. **割れの相手が 3 つ以上の連鎖**（A–B と B–C の両方が組になる）: まとめは連結成分ごとに 1 つ、二重にまとめない（Task 5）。
4. **LBM の設定が無いのに `msp_then_lbm`**: `LBM_NOT_CONFIGURED` で止まり、置き場所を出さない（Task 4）。
5. **陰イオンの付加体の別名**: `[M+FA-H]-` と `[M+HCOO]-` が同じ質量で二重の候補にならない（Task 9）。

---

### Task 1: 起動検査の時間切れを分ける

**Files:**
- Modify: `metabolomix/console/runner.py`（`probe_console_exe` を足し、`is_console_exe` はその包み）
- Modify: `metabolomix/pipeline/inputs.py`（`_resolve_exe`、`plan_from_profile` の 2 か所）、`metabolomix/tools/console_tools.py`（`console_plan`・`console_run`）
- Modify: `metabolomix/core/mcp_errors.py`（`CONSOLE_ERROR_CODES` に `MSDIAL_EXE_PROBE_TIMEOUT`）、`USAGE.md`（`console_plan`・`pipeline_plan` の行）
- Test: `tests/test_console_runner.py`、`tests/test_pipeline_inputs.py`

**Interfaces:**
- Produces: `probe_console_exe(exe_path: str, timeout_s: int = 60, retries: int = 1, *, raise_on_os_error: bool = False) -> Literal["console", "not_console", "timeout"]`。
  `is_console_exe(...) -> bool` は `probe_console_exe(...) == "console"` のまま残す（既存の monkeypatch を壊さない）。
- Produces: エラーコード `MSDIAL_EXE_PROBE_TIMEOUT`（details `{"timeout_s", "attempts"}`）。

- [ ] **Step 1: 失敗するテストを書く** — `test_probe_returns_timeout_after_retry`（`subprocess.run` を毎回 `TimeoutExpired` にする monkeypatch で、
  呼び出し 2 回・戻り値 `"timeout"`）、`test_probe_retry_succeeds`（1 回目時間切れ・2 回目 `lcms` を含む出力 → `"console"`）、
  `test_resolve_exe_reports_probe_timeout`（`probe_console_exe` を `"timeout"` に差し替え → `DomainError.code == "MSDIAL_EXE_PROBE_TIMEOUT"`）、
  `test_console_plan_reports_probe_timeout`。
- [ ] **Step 2: 失敗を確かめる** — Run: `.venv/Scripts/python.exe -m pytest tests/test_console_runner.py tests/test_pipeline_inputs.py -q -k probe` / Expected: FAIL
- [ ] **Step 3: 実装する** — 3 つの呼び出し元は `probe_console_exe` の戻り値で `NOT_CONSOLE` と `PROBE_TIMEOUT` を分ける。
- [ ] **Step 4: 全数を回す** — Expected: PASS
- [ ] **Step 5: コミット** — `fix(console): 起動検査の時間切れを MSDIAL_EXE_PROBE_TIMEOUT として分け、60 秒・1 回やり直す`

### Task 2: 検出状態の取り込み（末尾の空白）と、pipeline の run のハンドルの `.arf2`

**Files:**
- Modify: `metabolomix/mztab/evidence.py`（`build_evidence` の名前の照合を `strip()` した名前どうしで行う）、`metabolomix/analysis/assay_evidence.py:92`（`_normalize_stem` も strip）
- Modify: `metabolomix/pipeline/session_handoff.py`（`_handoff_rejection` / `restore_differentials` の互換）
- Modify: `metabolomix/tools/handle_tools.py`（`_alignment_files` に pipeline の run の `input/AlignResult-*.arf2` を探す段を足す）
- Test: `tests/test_dataset_evidence.py`、`tests/test_session_handoff.py`（無ければ既存の handoff のテストのファイル）、`tests/test_dataset_handle_tools.py`

**Interfaces:**
- Produces: `dataset_handle_open(pipeline_path=…)` で `arf2`・`dcl` が付く（戻り値の `arf2`・`dcl` が null でない）。
- Produces: 差次的結果の互換の規則: 保存された結果の `dataset_fingerprint` が「検出状態を除いて計算した今の fingerprint」と一致すれば受け入れ、
  結果に `evidence_attached_after: true` を付ける。

- [ ] **Step 1: 失敗するテストを書く** — `test_build_evidence_matches_names_with_trailing_space`（`.arf` の FileName が `'X-25 '`、mzTab の試料名が `'X-25'` → status ok）、
  `test_restore_accepts_result_computed_before_evidence`（検出状態なしで計算した結果を、検出状態ありの ds に復元 → 捨てられず `evidence_attached_after` が真）、
  `test_restore_still_rejects_other_mismatch`（行列の値が違えば従来どおり `dataset_mismatch`）、
  `test_handle_open_finds_pipeline_arf2`（`pipeline_path` の `input/AlignResult-1.arf2` と `.dcl` を置いた一時の run → 戻り値に両方）。
- [ ] **Step 2: 失敗を確かめる** — Expected: FAIL
- [ ] **Step 3: 実装する** — `result_state.dataset_fingerprint` に `include_detected: bool = True` を足し、復元の照合で 2 通りを比べる。
- [ ] **Step 4: 実データで確かめる** — MTBLS6965 の陽イオンの run を `dataset_handle_open(pipeline_path=…)` で開くスクリプト（Git の外の scratch）で、
  `feature_qc.source == "arf"`、`arf2`・`dcl` が null でない、`dataset_feature_search(feature_ids=["1266"])` の stats に q があることを確かめる。
- [ ] **Step 5: 全数を回してコミット** — `fix(evidence): 末尾の空白のある試料名で検出状態を取り込み、pipeline の run の .arf2 をハンドルが探す`

### Task 3: アライメントの参照ファイルを注入順の中央にする

**Files:**
- Modify: `metabolomix/pipeline/request.py`（`alignment_reference` を `_TOP_LEVEL_KEYS`・`_NULL_REJECTED_TOP_LEVEL_KEYS` に足し、値 `{"median_injection","method"}`、既定 `"median_injection"`、`_FINGERPRINT_OMITTED_DEFAULTS` に入れる）
- Modify: `metabolomix/pipeline/service.py`（`_prepare_run` で実験情報シートを先に解し、行を `inspect_inputs` に渡す。`_resolved_settings` に `alignment_reference`）
- Modify: `metabolomix/pipeline/inputs.py`（`inspect_inputs(..., manifest_rows: list[dict] | None = None)`、overrides に `Alignment reference file ID`）
- Modify: `metabolomix/console/method_file.py`（`ALIGNMENT_REFERENCE_KEY = "Alignment reference file ID"`）
- Test: `tests/test_pipeline_inputs.py`、`tests/test_pipeline_request.py`、`tests/test_pipeline_service.py`

**Interfaces:**
- Produces: `choose_alignment_reference(primary_names: list[str], manifest_rows: list[dict] | None) -> dict` を `pipeline/inputs.py` に。
  戻り値 `{"file_id": int, "sample_file": str, "injection_order_source": "manifest" | "filename", "method_value": str | None}`。
  `primary_names` は入力フォルダの主ファイル名を `_collect_primaries_and_companions` と同じ順に並べたもの（Console の file ID の順）。
  中央は全試料（include を問わない）の `injection_order` の中央値の試料。偶数なら小さい方。`injection_order` が 1 件でも欠ければファイル名の順で中央。
- Produces: receipt `resolved.alignment_reference` に上の dict。`alignment_reference: "method"` のときは上書きせず `method_value` だけ残す。

- [ ] **Step 1: 失敗するテストを書く** — `test_reference_is_median_injection`（ファイル名順 a,b,c,d,e に注入順 5,1,3,2,4 → 注入順 3 の c、file_id 2）、
  `test_reference_even_count_takes_lower_median`、`test_reference_falls_back_to_filename_order`（注入順の欠け → `injection_order_source == "filename"`）、
  `test_reference_method_mode_keeps_method_value`、`test_default_alignment_reference_does_not_change_the_request_fingerprint`、
  `test_plan_receipt_shows_alignment_reference`。
- [ ] **Step 2: 失敗を確かめる** — Expected: FAIL
- [ ] **Step 3: 実装する**
- [ ] **Step 4: 実データの順を確かめる** — MTBLS6965 の陽イオンの `.arf` の `file_id` の順が `sorted()` の主ファイル名の順と一致することを、
  `evidence.normalize_arf_spots` で読む scratch のスクリプトで確かめ、結果をテストのコメントに 1 行書く。
- [ ] **Step 5: 全数・Vault の流れ図・コミット** — `feat(pipeline): アライメントの参照ファイルを注入順の中央の試料にする（alignment_reference）`

### Task 4: `library_mode = "msp_then_lbm"`

**Files:**
- Modify: `metabolomix/core/user_config.py`（`SETTINGS` に `library.lbm`）、`lipidmix.example.toml`（`[library]` に `# lbm = ...` とコメント）
- Modify: `metabolomix/pipeline/request.py`（`_LIBRARY_MODE_VALUES` に `msp_then_lbm`、`lbm_file` との同時指定を拒否）
- Modify: `metabolomix/console/method_file.py`（`msp_then_lbm_overrides`、注釈器の設定表の書き出し）
- Modify: `metabolomix/pipeline/inputs.py`（`inspect_inputs` の分岐、`_resolve_lbm_from_library_setting`、`stage_inputs` で設定表を `inputs/msp-annotators.tsv` に書く、`_verify_inputs` で再ハッシュ）
- Modify: `metabolomix/pipeline/service.py`（`_resolved_settings` の `lbm.source`・`msp.priority`、`_plan_fingerprint` に lbm の sha256）
- Modify: `metabolomix/core/mcp_errors.py`（`LBM_NOT_CONFIGURED`・`LBM_NOT_FOUND` が無ければ足す）
- Test: `tests/test_pipeline_inputs.py`、`tests/test_console_method_file.py`、`tests/test_pipeline_service.py`、`tests/test_user_config.py`

**Interfaces:**
- Produces: `msp_then_lbm_overrides(method_keys: dict[str, str], annotator_table_path: str, lbm_path: str) -> tuple[dict[str, str], list[str]]`。
  上書き: `Msp file path` = 空、`MSP annotator settings file path` = 設定表、`Lbm file path` = LBM、
  `Use retention information for LBM-based annotation scoring` = `False`、`Use retention information for LBM-based annotation filtering` = `False`。
  `IDENTIFICATION_CLEAR_KEYS` のほかの宣言（Text DB・Text の設定表）は消して、消したキーを返す。
- Produces: `write_msp_annotator_table(path: Path, msp_path: str, *, priority: int = 2, annotator_id: str = "LabMsp") -> None`。
  TSV のヘッダは `annotator_id\tmsp_file_path\tpriority`（ASCII。Console 45a531c4d の `ReadMspAnnotatorSettingsTable` が読む形）。
- Produces: receipt `resolved.library_mode == "msp_then_lbm"`、`resolved.lbm == {"file", "sha256", "source": "config_file"}`（path は出さない）、
  `resolved.msp.priority == 2`。

- [ ] **Step 1: 失敗するテストを書く** — `test_msp_then_lbm_plan_pins_msp_and_lbm`（設定に msp_positive と lbm → plan の `msp.sha256`・`lbm.sha256`、
  上書きに 5 つのキー）、`test_msp_then_lbm_writes_annotator_table_with_priority_2`（stage 後の `inputs/msp-annotators.tsv` の中身）、
  `test_msp_then_lbm_without_lbm_setting_stops`（`LBM_NOT_CONFIGURED`、メッセージと details に置き場所が無い）、
  `test_msp_then_lbm_rejects_lbm_file`、`test_receipt_shows_lbm_without_location`、`test_default_library_mode_does_not_change_the_request_fingerprint`（既存を維持）。
- [ ] **Step 2: 失敗を確かめる** — Expected: FAIL
- [ ] **Step 3: 実装する** — LBM は `[library] lbm` だけから解決する（`[msdial] lbm`・メソッドの宣言・実行体のフォルダは見ない）。
- [ ] **Step 4: 全数・`USAGE.md` の `pipeline_plan` の `library_mode` の説明・Vault・コミット** — `feat(pipeline): 研究室 MSP を優先し LBM で補う library_mode=msp_then_lbm`

### Task 5: 割れたスポットの検出 `dataset_alignment_qc`

**Files:**
- Create: `metabolomix/analysis/alignment_qc.py`
- Modify: `metabolomix/tools/reanalysis_tools.py`（ツール `dataset_alignment_qc`）、ツール数のテストと文書（Global Constraints）
- Test: `tests/test_alignment_qc.py`、`tests/test_reanalysis_tools.py`

**Interfaces:**
- Produces: `find_split_pairs(ds, *, ppm: float = 15.0, min_rt_gap: float = 0.05, max_rt_gap: float = 1.0, max_both_detected: float = 0.10, injection_order: dict[str, int] | None = None) -> dict`。
  戻り値 `{"schema": "alignment-qc.v1", "params", "injection_order_source": "manifest"|"metadata"|"filename", "pairs": [...], "groups": [...], "summary": {...}}`。
  - `pairs[]`: `{"pair_id", "features": [a, b], "mz": [..], "rt": [..], "delta_rt", "adduct", "names": [..], "names_compatible": bool, "type": "injection_block"|"group_aligned"|"interleaved", "owner": {sample: a|b}, "both_detected_fraction", "block": {"spot", "injections": [first, last]} | null, "mergeable": bool}`。
  - `groups[]`: 組の連結成分（3 スポット以上の連鎖を 1 つにする）`{"group_id", "features": [...], "mergeable": bool}`。`mergeable` は成分のすべての組が mergeable のとき真。
  - `summary`: `{"n_pairs", "by_type": {...}, "n_features_affected", "injection_block_median_delta_rt"}`。
  - 受け持ち（owner）は、検出状態があれば検出、無ければ 2 スポットの値の大きい方。`both_detected_fraction` は検出状態が無ければ null で、
    そのとき「補い合い」は「片方の値が他方の 3 倍以上の試料が 90% 以上」で判定する。
  - 型: 片方の受け持ちの試料が注入順で連続した塊（塊の外に受け持ちが無い）→ `injection_block`（群と重なっても優先）、そうでなく受け持ちが群と一致 → `group_aligned`、ほか → `interleaved`。
  - `names_compatible`: 同じ名前、`lipid_identity.normalize_lipid_name` で同じ和組成、または片方が注釈なし。
  - `mergeable` = 型が `injection_block` か `interleaved` で、かつ `names_compatible`。
- Produces: MCP ツール `dataset_alignment_qc(output_path: str | None = None, ppm: float = 15.0, max_rt_gap: float = 1.0, dataset_handle: str | None = None) -> str`
  （要約だけを返し、全体は `output_path` に JSON で書く）。注入順はハンドルの試料のメタデータ（pipeline の実験情報シート）から取り、無ければファイル名の順。
- Produces: `split_partners(qc: dict) -> dict[str, dict]`（特徴 ID → `{"group_id", "partners": [...], "type", "mergeable"}`）。Task 6・7 が使う。

- [ ] **Step 1: 失敗するテストを書く**（合成の ds。`tests/test_alignment_qc.py`）— `test_injection_block_first_two_injections`（42 試料、最初の 2 注入だけが RT +0.6 の別スポット → 型 `injection_block`、block の注入 [1, 2]、mergeable）、
  `test_interleaved_pair`、`test_group_aligned_not_mergeable`（群で受け持ちが分かれ、注入順の塊ではない → mergeable 偽）、
  `test_incompatible_names_not_mergeable`、`test_three_spot_chain_is_one_group`、`test_far_mz_or_rt_not_paired`、
  `test_both_detected_above_threshold_not_split`（検出状態ありで両方検出 40% → 組にしない）、`test_tool_writes_output_and_returns_summary`。
- [ ] **Step 2: 失敗を確かめる** — Expected: FAIL
- [ ] **Step 3: 実装する**
- [ ] **Step 4: 実データで確かめる** — MTBLS6965 の陽イオン・陰イオンのハンドルで回す scratch のスクリプト。期待: 陽イオンで `injection_block` の組が約 230（ΔRT の中央値 +0.5〜+0.6）、
  CE の 1265/1266・1185/1186・1178/1179・1257/1258 がすべて mergeable、陰イオンの 401/403 が組になる。違えば閾値を調整し、調整をテストのコメントに書く。
- [ ] **Step 5: 全数・文書・コミット** — `feat(reanalysis): 割れたスポットを探す dataset_alignment_qc`

### Task 6: spec の統計で割れをまとめる・BH の母集団 `annotated`

**Files:**
- Modify: `metabolomix/analysis/analysis_spec.py`（`_normalize_features` に `merge_split`、`_normalize_analysis` の `bh_population` に `annotated`、`_run_analysis` でまとめる）
- Modify: `docs/output_format/reanalysis.md` §2
- Test: `tests/test_analysis_spec.py`

**Interfaces:**
- Consumes: Task 5 の `find_split_pairs`・`split_partners`。
- Produces: spec `features` は `{"all": true}` か `{"feature_ids": [...]}` に、任意で `"merge_split": true|false`（v8 の既定は **true**。キーが無ければ true）と
  `"alignment_qc": <JSON のパス>`（無ければその場で `find_split_pairs(ds)` を既定の値で回す）。
  まとめは `build_matrix` の後・`_aggregate` と検定の前。mergeable な群ごとに、試料ごとの最大値の 1 列にし、列の ID は群の中で最も小さい特徴 ID。
  出力の行に `merged_from: [...]`。`aggregate` の `select` で選ばれた特徴も同じく 1 つにまとめてから和をとる。
- Produces: `bh_population: "annotated"` = SME のある特徴（`ds.feature_candidates` にある）だけを BH の母集団にする。
- `rule_hash` には `merge_split` と、使った割れの群の一覧の hash を入れる。

- [ ] **Step 1: 失敗するテストを書く** — `test_merge_split_takes_per_sample_max`、`test_merge_split_false_keeps_both`、`test_merge_split_default_is_true`、
  `test_aggregate_counts_split_group_once`（C6 の型: 和に 1185/1186 の両方が選ばれても 1 回）、`test_merged_from_in_rows`、
  `test_bh_population_annotated`、`test_rule_hash_changes_with_merge_split`。
- [ ] **Step 2〜4: 失敗 → 実装 → 全数**
- [ ] **Step 5: コミット** — `feat(analysis_spec): 割れたスポットを 1 つにまとめて検定し、BH の母集団に annotated を足す`

### Task 7: 目標照合・証拠の束・新しい知見の候補の割れ、`counts` の直し

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`（`merge_split` の再検定）、`metabolomix/tools/dataset_analysis_tools.py`（`dataset_feature_search` の引数）
- Modify: `metabolomix/analysis/claim_evidence.py`（候補に `split`、`counts` の直し）、`metabolomix/analysis/mass_accuracy.py`（`novel_candidates` のまとめ）
- Modify: `docs/output_format/reanalysis.md` §7・§8、`docs/workflow/reanalysis.md`
- Test: `tests/test_feature_search.py`、`tests/test_dataset_feature_search_tool.py`、`tests/test_reanalysis_tools.py`

**Interfaces:**
- Produces: `dataset_feature_search(..., merge_split: bool = True, alignment_qc_path: str | None = None, dataset_handle=None)`。
  `merge_split` が真で、照合した特徴に mergeable な割れの群があるとき、その群を試料ごとの最大値でまとめた行列で、run の前処理のレシピ
  （`session_handoff._run_preprocess_parameters` の `resolved_recipe`）を適用し直し、`differential.two_group_test(..., log_transform=True)` と
  `add_fdr` を **まとめた後の全特徴** で回した stats を返す。行に `stats`（まとめた後）と `stats_unmerged`（従来の値）、`merged_from` を並べる。
  割れの無い特徴の `stats` は従来の値と同じになる。
- Produces: 証拠の束の候補に `"split": {"group_id", "partners", "type", "mergeable"} | null`。
  `default_selection.counts`: `qualified` = 候補のうち `msms.qualify.verdict == "qualified"` の数、旧 `qualified` の意味の数は `msms_available` に改名。
- Produces: `dataset_novel_candidates` は mergeable な割れの群を 1 候補に数え、`merged_from` を付ける。

- [ ] **Step 1: 失敗するテストを書く** — `test_feature_search_merge_split_recomputes_q`（合成: 割れた 2 スポットを持つ物質 → まとめた stats の q と `stats_unmerged` が並ぶ）、
  `test_feature_search_no_split_stats_unchanged`、`test_feature_search_merge_split_false`、`test_claim_evidence_candidate_has_split`、
  `test_counts_qualified_is_msms_verdict`（既定の特徴が `not_qualified` → `qualified` 0・`msms_available` 1）、`test_novel_candidates_count_split_once`。
- [ ] **Step 2〜4: 失敗 → 実装 → 全数**
- [ ] **Step 5: 実データで確かめる** — MTBLS6965 の陽イオンで CE 20:3 の目標（`lipid_species: "CE 20:3"`、RT 12.79 ±0.3）の照合が、まとめた後に 1 つの物質として q を出す。
- [ ] **Step 6: コミット** — `feat(reanalysis): 目標照合・証拠の束・新しい知見の候補で割れたスポットをまとめ、counts.qualified を MS/MS の合否にする`

### Task 8: 質量精度の絞り込みと注入順

**Files:**
- Modify: `metabolomix/mztab/dataset_state.py`（SML の `opt_global_spectrum_matched`・`opt_global_ms2_presence` を読む。`_SML_TEXT_FIELDS` と注釈の辞書）
- Modify: `metabolomix/mztab/evidence.py`（`apply_evidence` で試料ごとの m/z を `ds.per_sample_mz: dict[str, dict[str, float]] | None` に残す）
- Modify: `metabolomix/analysis/mass_accuracy.py`、`metabolomix/tools/reanalysis_tools.py`（`dataset_mass_accuracy`）
- Modify: `docs/output_format/mztab.md` §13、`docs/output_format/reanalysis.md` §8
- Test: `tests/test_reanalysis_tools.py`、`tests/test_mztab_*.py`（SML の列を読むテストのあるファイル）

**Interfaces:**
- Produces: `dataset_mass_accuracy(band_edges=None, rank=1, subset: Literal["rank1", "msms_matched"] = "rank1", dataset_handle=None)`。
  `msms_matched` = rank 1 の SME のうち、その特徴の SML の `opt_global_spectrum_matched` が真のもの（LBM・MSP を区別しない）。
- Produces: 戻り値に `injection_trend`: 試料ごとの m/z があれば `{"per_injection": [{"sample", "injection_order", "median_ppm", "n"}], "spearman_rho", "n_samples"}`、無ければ `{"status": "unavailable"}`。
  試料ごとの誤差は、その試料で検出された（gap-fill でない）セルの m/z と理論 m/z の差。

- [ ] **Step 1: 失敗するテストを書く** — `test_mass_accuracy_msms_matched_subset`、`test_sml_spectrum_matched_is_read`、`test_injection_trend_unavailable_without_per_sample_mz`、
  `test_injection_trend_spearman`（合成の注入順に比例する誤差 → ρ ≈ 1）。
- [ ] **Step 2〜4: 失敗 → 実装 → 全数**
- [ ] **Step 5: コミット** — `feat(reanalysis): 質量精度を MS/MS で照合した同定に絞り、注入順との関係を返す`

### Task 9: 陰イオンの付加体の層と、パスウェイの handoff

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`（`DEFAULT_ADDUCT_TIERS["negative"]` の 2 層目に `[M+CH3COO]-`）
- Modify: `metabolomix/pathway/handoff.py`（MS-DIAL の脂質の Ontology の許可表と、名前の検査）
- Create: `metabolomix/pathway/msdial_lipid_classes.py`（MS-DIAL 5 の `LbmClass` の名前の一覧。出所を `MsdialWorkbench` の `LbmClass` の enum として書く）
- Test: `tests/test_feature_search.py::test_default_tiers_values_are_the_spec_ones`（値を直す）、`tests/test_pathway_handoff.py`

**Interfaces:**
- Produces: 陰イオンの既定の層 `[["[M-H]-"], ["[M-H2O-H]-", "[M+Cl]-", "[M+FA-H]-", "[M+CH3COO]-", "[2M-H]-", "[M+Na-2H]-"]]`。
  `[M+HCOO]-` は `[M+FA-H]-` と同じ質量なので足さない（設計 4.9 の「`[M+HCOO]-`」はこれで満たす）。`theoretical_mz` が `[M+CH3COO]-` を読めることをテストで確かめる。
- Produces: handoff の `report` に `dropped_non_msdial_ontology` と `dropped_unparsable_name`。Ontology が許可表に無い行、`normalize_lipid_name(name)["parse_ok"]` が偽の行を落とす。

- [ ] **Step 1: 失敗するテストを書く** — `test_negative_tiers_include_acetate`、`test_acetate_adduct_mass`、`test_handoff_drops_classyfire_ontology`（Ontology `Long-chain fatty acids` の行が落ちる）、
  `test_handoff_drops_unparsable_names`（`3-Sulfamoylbenzoic acid`）、`test_handoff_keeps_msdial_lipids`（Ontology `PC`、名前 `PC 38:3`）。
- [ ] **Step 2〜4: 失敗 → 実装 → 全数**
- [ ] **Step 5: コミット** — `fix: 陰イオンの既定の付加体に酢酸付加を足し、パスウェイの handoff が脂質でない行を除く`

### Task 10: Console の小さな試し（動作確認。コードの変更なし）

**Files:**
- Create（Git の外）: `C:\Users\yuu18\datasets\reanalysis\_v8_trial\`（MTBLS6965 の陽イオン 6 ファイル・陰イオン 6 ファイルのハードリンク、実験情報シート、メソッドの写し）

- [ ] **Step 1**: ユーザーの承認の後に、設定ファイル（追跡外）の `[library] lbm` に conventional の LBM を書く（このタスクの実行役は書かない。オーケストレーターが行う）。
- [ ] **Step 2**: worktree のコードを Python から直接呼ぶ scratch のスクリプトで、陽・陰の `pipeline_plan` → `pipeline_run` を **同時に** 起動する（`library_mode: msp_then_lbm`、既定の `alignment_reference`）。
- [ ] **Step 3**: 確かめる — 両方が完走する。receipt の `alignment_reference.sample_file` が注入順の中央の試料。effective-method に 5 つのキー。
  mzTab に LBM 由来の同定（名前が脂質の略記）と研究室 MSP の同定の両方がある。同じ特徴に両方の候補があるとき、代表（rank 1）が MSP のものになる例を 1 つ以上見つける。
  2 本の同時実行のメモリの最大（タスクマネージャか `Get-Process MSDIALCUI` の WorkingSet の記録）。
- [ ] **Step 4**: 結果を `docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md`（短い記録）に書いてコミット。期待と違えば止めて報告する。
