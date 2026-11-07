# B-2b 公開面の切り替え（新しいツール面・旧ツールとスロットの削除・知識層の削除・data_dir の探索の削除・文書）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** spec B §8 の対応表の B の分を実装し、`session.arf` / `.arf2` / `.pai2` / `.eic` のスロットとそれを読む旧ツール・文献と知識の層・`data_dir` の探索を消し、文書とテストを新しいツール面に合わせる。

**Architecture:** データセットを扱うツールはすべて `DatasetState`（ハンドル）から読む。入口は `load_dataset` 1 本（`opening.open_dataset` を呼び、ハンドルを作って直近にする）、ハンドルの一覧と終了は `dataset_list` / `dataset_close`。パスを省いたツールは直近のデータセットのファイル（`artifact_paths` とハンドルの `.arf2` / `.dcl`）を使い、`data/` を探さない。描画は `dataset_plot(kind=...)` 1 本で、作った図の payload はデータセットの図の置き場（`DatasetState.figures`）に置き、`save_figure` がそこから描く。記録（目的・レポート・図・パスウェイの出力）は直近のデータセットの解析フォルダの `reports/` に書く。

**Tech Stack:** Python 3.14・numpy・matplotlib（Agg）・pytest。MCP は FastMCP。

**Spec:**
- `docs/superpowers/specs/2026-11-06-unified-dataset-and-tools-design.md`（B。この計画は §5 の残り・§7.2 の残り・§8（B の分）・§9・§10）
- `docs/superpowers/specs/2026-11-06-overhaul-overview-design.md`（全体）
- 前提: `docs/superpowers/plans/2026-11-07-b2a-foundation.md`（B-2a）が済んでいること

## 決めたこと（設計が任せた点・利用者の判断）

1. **B-2 の終わりのツール面（56 ツール・リソース 2・テンプレート 1）**。spec §8 の 53 との差: B-3 の 4 つ（`dataset_quant_*` `dataset_merge_polarities` `dataset_select_features`）・C の `dataset_compare`・A の `dataset_overview` はまだ無く、C の旧統計 3 つ（`dataset_differential` `dataset_statistic` `dataset_build_matrix`）と D の旧 Console 6 つ（`console_plan` `console_run` `console_status` `console_prepare_input` `console_cleanup` `job_list`）はそれぞれの段が消すまで残す（C・D の範囲。B が先に消すと、C・D ができるまで差次的解析と Console の単独実行が無くなるため）。
   - 入口: `list_data_files` `load_dataset` `dataset_list` `dataset_close`
   - 試料: `dataset_status` `dataset_set_sample_metadata` `dataset_exclude`
   - キュレーション: `curation_review` `curation_suggest` `curation_submit` `curation_flags` `curation_view_data`
   - 前処理・統計（C で作り直す）: `dataset_preprocess` `dataset_pca` `dataset_analysis_run` `dataset_differential` `dataset_statistic` `dataset_build_matrix`
   - 書き出し: `dataset_export_differential` `dataset_export_evidence_bundle`
   - 特徴と根拠: `dataset_feature_search` `dataset_feature_relations` `dataset_claim_evidence` `dataset_alignment_qc` `dataset_mass_accuracy` `dataset_homolog_check` `dataset_novel_candidates`
   - スペクトルと照合: `dcl_spot_spectrum` `dcl_diagnostic_ions` `dcl_find_msms` `library_load` `library_match_feature` `library_msms_qualify` `library_plot_mirror` `library_search`
   - 描画: `dataset_plot` `save_figure`
   - パスウェイ: `pathway_project` `pathway_plot`
   - 記録: `record_objective` `write_report` `read_report`
   - Console・pipeline: `pipeline_plan` `pipeline_run` `pipeline_status` `pipeline_resume` `pipeline_cancel` `console_method_candidates` `console_method_template` と旧 6 つ
   - 保守: `server_update`
   - リソース: `lipidmix://docs/output-format`、`ui://ms-data-parser/curation-viewer`。テンプレート: `lipidmix://docs/output-format/{topic}`。
2. **`load_dataset` の形**: 引数は今の `dataset_handle_open` と同じ（`directory` `mztab_path` `job_path` `pipeline_path` `allow_incomplete` `source` `measure` `arf2_path` `library="auto"|"none"`）。開いたハンドルは**いつも直近になる**（対話も再解析も同じ。再解析のエージェントはハンドルを渡して呼ぶので直近を使わない。stdio のサーバはクライアントごとに別のプロセス）。戻り値は今の `dataset_handle_open` の JSON に、試料の群・役割の件数（`dataset_status` の要約と同じ形）と `next`（次の一手の例 3 行）を足す。**読み込みの後の概観（`.arf2` の要約と前処理前の PCA）は返さない**（前処理前の PCA は要らない: 利用者の判断。概観は A の `dataset_overview` が来るまで無い）。`ALWAYS_THREADED`。ハンドルの ID の形 `dh_…` と引数名 `dataset_handle` は変えない。
3. **記録の置き場所**: `DatasetState.analysis_dir`（読み込みが決める。解析フォルダのアダプタは選んだ `.arf` のフォルダ、mzTab-M のアダプタは mzTab-M のフォルダ）。書く先の候補は `[paths] reports_dir`（あれば）→ `<analysis_dir>/reports`。どちらも無い（設定が無く、データセットも無い）ときは `missing_state("dataset", ["load_dataset"])`。`record_objective` `write_report` `read_report` `save_figure` `pathway_project` が使う。ハンドルを選べるよう `dataset_handle` を持つ。
4. **`data_dir` の探索を消す（利用者の判断 2026-11-07）**: `mcp_core.DATA_DIR`・`core/data_config.py`・`LIPIDMIX_DATA_DIR` を消す。`list_data_files(directory)` は `directory` 必須。パスを省いたツールは直近のデータセットのファイルを使い、無ければ `missing_state("dataset", ["load_dataset"])`（または必要なファイルが無いときは今のハンドル経路と同じ `no_alignment_files`）。明示のパス引数（`dcl_find_msms(file_path)`・`library_load(file_path)`・`library_match_feature(dcl_file)` など）は残し、探索だけを消す。パーサの CLI（`python -m metabolomix.<形式>.reader`）は `--file` 必須。`console/job_manager._assert_not_in_repo` の許可は `mcp_core.BASE_DIR / "data"` の定数で書く。
5. **`library_load` の極性の経路を消す（spec §7.2）**: `ion_mode` 引数と、設定の極性の MSP・フォルダの探索への後退を消す。`file_path` を省けばデータセットが覚えたライブラリ。覚えていない（`unknown`）なら `LIBRARY_NOT_REMEMBERED`（候補のファイル名と「`file_path` を渡す」の案内）。`mismatch` / `lbm_only` は今どおり理由を返す。極性ごとの MSP の設定は Console にどれを使わせるかにだけ残る。
6. **`record_objective` が `update_objective` を吸収**: 目的が無ければ作る（今の `record_objective` の必須の引数が要る）。あれば渡された欄だけを直す（`status` と `add_subquestions` を受ける）。文献の検索の記録（`searched_labels` など）の欄は消す。**`read_report(analysis_id=None)`**: ID を省けば一覧（今の `list_reports` の行）。
7. **ノートの読み書きの関数を corpus から移す**: `make_slug` `parse_frontmatter` `dump_frontmatter` `write_note` `note_path` と目的の関数（`write_objective` `update_objective_meta` `add_subquestions` `parse_objective` `objective_slug` `build_objective_body`）を新しい `metabolomix/tools/notes.py`（stdlib だけ）へ移し、`metabolomix/corpus/` を消す。
8. **知識の種ノートの移し先（利用者の判断 2026-11-07）**: 追跡されている `knowledge/` の 3 本と `playbook/` の 2 本を `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\archived-knowledge\` へ写してから、リポジトリから消す（`git rm`）。
9. **キュレーションは `DatasetState` の上で動く**: 5 つのツールが `dataset_handle` を持ち、`file_path` を消す。`.arf2` はデータセットの `source_files` の `.arf2`（ハンドルの `arf2_path`）から取る。`curation_submit` が書いた後、同じ `.arf2`（sha256 で照合）を持つ開いているハンドルすべての `curation_flags` を読み直す（書き出しの `apply_curation` が新しい判断を使う）。キュレーションの依存（`session.library.store`）はデータセットのライブラリ（`holder.library`）を読む。
10. **`dataset_status(view="summary"|"samples"|"tags", query=None)`**: `summary` は今の要約に、群（Class ID）ごとの件数・役割ごとの件数・タグの件数を足す（旧 `arf_list_classes` `arf_list_sample_roles` `arf_list_tags`）。`samples` は試料の TSV（name・group・role・batch・injection_order・subject_id・excluded）で、`query` は因子のトークン（旧 `sample_search` の `specs` と同じ文法。`msdial/sample_factors.py`）で絞る。`tags` はタグごとの特徴の ID（上位 50）。
11. **`dataset_exclude(samples=None, features=None, mode="add"|"remove"|"clear", reason=None)`**: `DatasetState.excluded_samples: dict[str, str]` と `excluded_features: dict[str, str]`（値は理由）に持つ。可逆。変えたら `metadata_revision` を進め、前処理済み行列から先の結果を古いとする（`dataset_set_sample_metadata` の前処理入力が変わったときと同じ無効化）。`dataset_preprocess` と `dataset_plot` と `dataset_feature_search` の強度順は除外を当てる（C は同じ欄を読む）。
12. **`dataset_set_sample_metadata` に `subject_id`**: sample-manifest の任意の列 `subject_id` を読み、`sample_meta` に持つ（C §5 が使う。B では持つだけで統計に使わない）。
13. **`dataset_feature_search` が `eic_search` を吸収**: `sort="mz"|"intensity"`（既定 `mz`）。`intensity` は検出セル（`detected_mask`）の値の最大で降順（除外した試料を除く）。条件が無くても `sort="intensity"` なら上位 `max_results` 件を返す（`mz` で条件が無ければ今どおり `FEATURE_SEARCH_EMPTY_QUERY`）。行に `ontology` `n_detected` と、解析フォルダのデータセットでは同定の標準化の列 `normalized` `msi_level`（`feature_metadata["identity"]`。旧 `arf2_annotate_identities` の結果）を足す。EIC のクロマトグラムの最大ではなくピークの高さの最大であることを文書に書く。
14. **`dataset_plot` の `kind`**（形式名は旧ツールのものを保つ。Use-LLLM が読む）:

    | kind | 入力 | 旧ツール | payload の形式名 |
    |---|---|---|---|
    | `volcano` | `result_id` か直近の差次的結果（`ds.last_differential`） | `arf_plot_volcano` | `lipidmix.volcano.v1` |
    | `pca` | `result_id` か直近の `dataset_pca`（`ds.last_pca`）。`color_by="group"|"batch"|"role"` | （新） | `lipidmix.pca_scores.v1`（新設） |
    | `loadings` | 同じ PCA の結果、`pcs` `top_n` `value` | `plot_pca_loadings` | `lipidmix.pca_loadings.v1` |
    | `group_intensity` | `items`（クラス名か特徴 ID）、`groups` | `arf_plot_group_intensity` | `lipidmix.group_intensity.v1` |
    | `species` | `items`、`groups`、`value="share"|"intensity"` | `arf_plot_species` | `lipidmix.species_intensity.v1` |
    | `species_pca` | `items`、`groups`、`normalize` `log_transform` | `arf_pca_species` | `lipidmix.species_pca.v1`（新設） |
    | `eic` | `feature_id`、`file_ids` | `eic_plot_chromatograms` | `lipidmix.eic.v1` |
    | `eic_compounds` | `file_id`、`names` / `ontologies` | `eic_plot_compounds` | `lipidmix.eic.multi.v1` |

    既定は PNG と短い説明、`output="payload"` で点列、`LIPIDMIX_PLOT_OUTPUT`、`structured_output=False`。`group_intensity` `species` `species_pca` のキュレーションの当て方（`apply_curation` `exclude_auto_likely_wrong`）と `low_reliability_samples` `standard_samples` `detection_limit` `require_msms` は旧ツールの意味のまま引数にする。作った payload は `ds.figures[kind]` に置く（`save_figure(kind=...)` が描く）。EIC は `artifact_paths` の `eic` / `chromatogram` から引き、特徴 ID と EIC の AlignmentID の対応は今の `eic` の規則（イオンモビリティで連番でないときの注記）を保つ。
15. **消して代わりを作らないもの**: `verify_peak_annotation` と `pai2_inspect_peak`（spec §8）。失うもの——試料ごとのピークの S/N、付加体とイオンモードと脂質クラスの整合（`peak_verification.adduct_consistency`）、エーテル脂質の注意、MS/MS の `FLAG_ONLY` の帯、1 特徴の質量誤差——を HISTRY に書き、`docs/task.md` に「`dataset_claim_evidence` へ付加体の整合と 1 特徴の質量誤差を移すか」を候補として足す。`docs/output_format/core.md` の MS/MS の証拠の注意（PAI2 の `has_msms` は取得の参照があることだけ）は `.dcl` の注意として残す。
16. **`docs/output_format/` のトピック**: `arf` `arf2` `pai2` を消し、新しい `dataset`（解析フォルダのアダプタ: 行の粒度＝特徴×試料、gap-fill と `detected_mask`、`.arf2` の代表と同定の標準化の列、試料の役割、脂質名の文法への参照）にまとめ直す。`eic` は `dataset_plot` の EIC の説明として書き直して残す。行の粒度・gap-fill・脂質名の文法・必須の注意を失わないこと（各トピックの「必須の注意」の行を移した先を HISTRY に対応表で書く）。
17. **`MCP_INSTRUCTIONS` と注意の要約（`SEMANTICS_CAVEAT`）を書き直す**: ENTRY POINT は `list_data_files` → `load_dataset`、ASSAY KIND は `record_objective(assay_kind=...)`、GATEWAY の索引（knowledge / playbook）・CONFLICTS・LITERATURE DISCOVERY の節を消す、PARALLEL REANALYSIS は `load_dataset` / `dataset_list` / `dataset_close`。MS/MS の証拠の節は `dcl_*` と `library_msms_qualify` と `dataset_claim_evidence` で書く。消したツールの名前を書かない（Task 11 のテストで縛る）。
18. **既存の統計の関数（`analysis/differential.py` `statistics_v2.py` など）は動かさない・消さない**（C の範囲）。`pathway_project` の ARF の分岐（`_build_arf`）だけを消す。
19. **旧 6 ツールの文書**（`console_*` `job_list`）は今のまま（D の範囲）。

## Global Constraints

- コード・コメント・docstring・文書・コミットメッセージは日本語。
- テストは main ツリーの `.venv`: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。
- コミットはバックグラウンドで実行し、出力をファイルへ落とす。実行中にファイルを編集しない。`git stash` を使わない。push・`overhaul/integration` へのマージ・main の操作をしない。
- 旧名は一度に切り替える（移行期間なし、別名を残さない）。外部契約の旧名（`LIPIDMIX_` の接頭辞・`lipidmix://`・`lipidmix.<名前>.vN`・`lipidmix.local.toml`）は変えない。
- 差次的エクスポートの列（`analysis/export_contract.py`）と根拠バンドル v3 の形を変えない。
- 閉じた・存在しないハンドルを、対話のセッションや別のハンドルへ黙って落とさない。`dataset_handle` を持たないのに `session.dataset` / `.library` を読むツールは `handle_dispatch.SESSION_READERS` に載せる。
- ツールの戻り値は `json_payload()`、全ツール `structured_output=False`、float は丸める。図はサーバで描いて画像で返す。
- fixture はテスト自身が tmp に作る。`reports/` `knowledge/` の実ファイルに依存させない。
- 消すツールのテストは、そのツールと一緒に消す。消すツールが使っていた純粋なモジュール（`plots/*` の組み立て、各 `reader.py`）のテストは残す。
- 新しい環境変数を足さない。

## Review Focus

1. **直近のデータセットを閉じた後に、パスを省いて `dcl_spot_spectrum` / `curation_review` / `write_report` を呼ぶ** → `missing_state("dataset", ["load_dataset"])`。`data/` や開いている別のハンドルのファイルを使わない。（Task 9 `test_pathless_tools_do_not_fall_back_after_close`）
2. **`curation_submit` の直後に `dataset_export_differential(apply_curation=True)`** → 出したばかりの `wrong` の判断が書き出しに効く（同じ `.arf2` を持つ別のハンドルにも）。（Task 4 `test_submit_refreshes_flags_on_open_handles`）
3. **`dataset_exclude` の後に、除外前の `dataset_pca` の結果で `dataset_plot(kind="pca")`** → 古い結果として拒むか、古いことを示す（黙って除外前の図を最新として出さない）。（Task 3 `test_exclude_marks_downstream_results_stale`）
4. **mzTab-M だけで読んだデータセット（解析フォルダが無い）で `dataset_plot(kind="eic")` や `curation_review`** → 理由つきで「このデータセットには `.EIC.aef` / `.arf2` が無い」。例外にしない。（Task 6 `test_eic_without_artifact_is_explained`、Task 4 `test_curation_without_arf2_is_explained`）
5. **文書・`MCP_INSTRUCTIONS`・ツールの案内文（`required_tools` を含む）に消したツールの名前が残る** → 腐敗防止のテストで落ちる（`docs/superpowers/` の履歴は対象外）。（Task 11 `test_removed_tools_are_not_mentioned`）

---

### Task 1: 文献と知識の層を消す

**Files:**
- Create: `metabolomix/tools/notes.py`（決めたこと 7）
- Delete: `metabolomix/corpus/`、`tests/test_knowledge_store.py` `tests/test_paper_ingest.py` `tests/test_ingest_tools.py` `tests/test_playbook_tools_consistency.py`、追跡されている `knowledge/*.md` `playbook/*.md`（決めたこと 8 の写しの後）
- Modify: `metabolomix/tools/objective.py`（文献の 7 ツールを消す）、`metabolomix/tools/resources.py`（knowledge / playbook の 5 つを消す）、`metabolomix/tools/reports.py`（`notes` を import）、`metabolomix/core/mcp_core.py`（`KNOWLEDGE_DIR` `PLAYBOOK_DIR` `_state_dir` を消し、`_build_report_meta` の `knowledge_refs` を消す）、`metabolomix/core/tool_helpers.py`（`_build_verification_dossier` の `knowledge_store.coverage` と `candidate_knowledge_slugs`）、`metabolomix/pai2/tools.py`（`load_vocab` の呼び出し。ツールごと Task 8 で消えるが、それまで import が壊れないよう語彙なしにする）、`server.py`、`.gitignore`（knowledge / playbook の行）、`scripts/new-worktree.sh`（`LIPIDMIX_KNOWLEDGE_DIR`）、`tests/test_package_layout.py`（`playbook/` の存在・`KNOWLEDGE_DIR` の検査を外す）、`tests/test_objective_tools.py` `tests/test_report_tools.py`（import 先）、`tests/test_server_registration.py` `tests/test_tool_annotations.py` `tests/test_workflow_docs.py`（7 ツールと 5 リソースを外す）
- Test: `tests/test_notes.py`（`test_report_tools.py` の `make_slug` / `parse_frontmatter` の主張を移す）

- [ ] **Step 1**: vault のフォルダ `...\30_Projects\ms-data-parser\archived-knowledge\` を作り、5 本を写す（中身が同じことを sha256 で確かめる）。
- [ ] **Step 2: 失敗するテストを書く**: `test_notes_module_is_stdlib_only`、`test_no_corpus_imports`（`metabolomix/` に `corpus` の import が無い）、登録の検査（ツールが 7 減り、リソースが 3、テンプレートが 1）。
- [ ] **Step 3: 落ちることを確かめる。**
- [ ] **Step 4: 実装**（`git rm` を含む）。`record_objective` の戻り値の `knowledge_coverage` への案内を消す。
- [ ] **Step 5: 通ることを確かめる**（全体。`USAGE.md` と `docs/workflow/` の行も消すので `tests/test_readme_links.py` `tests/test_workflow_docs.py` が緑。USAGE の件数と `CLAUDE.md` 冒頭の規模表記をこの時点の登録数に直す）。
- [ ] **Step 6: コミット**: `refactor: 文献検索と知識の蓄積を消す（種ノートは vault へ移した）`

### Task 2: 入口を `load_dataset` 1 本にし、`dataset_list` / `dataset_close` へ改名する

**Files:**
- Modify: `metabolomix/tools/dataset.py`（旧 `load_dataset` を消し `list_data_files` を残す）、`metabolomix/tools/handle_tools.py`（`dataset_handle_open` → `load_dataset`、`dataset_handle_list` → `dataset_list`、`dataset_handle_close` → `dataset_close`）、`metabolomix/tools/mztab_tools.py`（`dataset_load` を消す）、`metabolomix/core/handle_dispatch.py`（`ALWAYS_THREADED = {"load_dataset"}`、`_call_pinned` の `required_tools`）、`required_tools` に `dataset_load` / `dataset_handle_open` を書いている全箇所（`grep -rn '"dataset_load"\|"dataset_handle_open"' metabolomix`）、docstring、`server.py`
- Test: `tests/test_load_dataset.py`（書き直し。旧 ARF 経路の 26 テストは消し、`dataset_load` / `dataset_handle_open` のテストの主張を移す）、`tests/test_session_handle_unification.py` `tests/test_dataset_handle_tools.py` ほか `dataset_load` / `dataset_handle_*` を呼ぶテスト（約 30 ファイル。呼び名を変えるだけで主張は変えない）

**Interfaces:**
- Produces: `load_dataset(directory=None, mztab_path=None, job_path=None, pipeline_path=None, allow_incomplete=False, source="auto", measure="peak_height", arf2_path=None, library="auto") -> str`（決めたこと 2）、`dataset_list() -> str`、`dataset_close(dataset_handle: str) -> str`（戻り値は今の `dataset_handle_list` / `_close` と同じ）。`DatasetState.analysis_dir: str | None`（決めたこと 3。`opening.open_dataset` が入れる）。
- [ ] **Step 1: 失敗するテストを書く**: `test_load_dataset_makes_current_handle`（今の `dataset_load` の主張）、`test_load_dataset_opens_remembered_library`（今の `dataset_handle_open` の主張）、`test_load_dataset_returns_group_and_role_counts`、`test_load_dataset_has_no_initial_pca`、`test_load_dataset_is_always_threaded`、`test_analysis_dir_per_adapter`（解析フォルダ → `.arf` のフォルダ、mzTab-M → mzTab-M のフォルダ、job / pipeline の解析フォルダ → 選んだ `.arf` のフォルダ）。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る**（全体。テストの呼び名の置き換えは機械的に行い、主張を変えたテストがあれば HISTRY に書く）。
- [ ] **Step 5: コミット**: `feat(tools): 入口を load_dataset にまとめ、dataset_list / dataset_close へ改名する`

### Task 3: 試料のツール（`dataset_status` の拡張・`dataset_exclude`・`subject_id`）

**Files:**
- Modify: `metabolomix/tools/mztab_tools.py`（`dataset_status`。ファイルに残るツールがこれだけなら `metabolomix/tools/dataset_status_tools.py` へ移して名前を合わせてよい）、`metabolomix/mztab/dataset_state.py`（`excluded_samples` `excluded_features`）、`metabolomix/analysis/dataset_service.py` か前処理の入口（除外を当てる）、`metabolomix/analysis/sample_manifest.py`（`subject_id` の列）、`metabolomix/tools/dataset_analysis_tools.py`（`dataset_set_sample_metadata` の docstring）
- Create: `metabolomix/tools/sample_tools.py`（`dataset_exclude`）
- Test: `tests/test_dataset_status_views.py`、`tests/test_dataset_exclude.py`、`tests/test_sample_manifest.py`（追記）

**Interfaces:**
- Produces: `dataset_status(dataset_handle=None, view="summary", query=None) -> str`（決めたこと 10）、`dataset_exclude(samples=None, features=None, mode="add", reason=None, dataset_handle=None) -> str`（決めたこと 11。存在しない名前・ID は全体を拒み、一部だけ当てない）。`sample_factors` のトークンの照合は旧 `tools/samples.py` が使っていた関数をそのまま使う。
- [ ] **Step 1: 失敗するテストを書く**: `test_summary_counts_groups_roles_tags`、`test_samples_view_filters_by_factor_tokens`（旧 `tests/test_tools_samples.py` の照合の主張を移す）、`test_exclude_is_reversible`、`test_exclude_unknown_name_refuses_whole_request`、`test_exclude_applied_by_preprocess`、`test_exclude_marks_downstream_results_stale`（Review Focus 3）、`test_subject_id_read_from_manifest`。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る。**
- [ ] **Step 5: コミット**: `feat(tools): dataset_status に群・役割・タグ・試料の検索、dataset_exclude、subject_id を足す`

### Task 4: キュレーションを `DatasetState` の上へ

**Files:**
- Modify: `metabolomix/tools/curation_tools.py`（5 ツールに `dataset_handle`、`file_path` を消す、`_candidate_dirs` の既定の探索を消す）、`metabolomix/curation/evidence.py`（`.dcl` / `.EIC.aef` を `artifact_path` から引く）、`metabolomix/dataset/msdial_adapter.py`（`_load_curation` を外から呼べる `reload_curation(ds) -> None` にする）、`metabolomix/core/handles.py`（同じ `.arf2` の sha256 を持つハンドルを探す `handles_for_arf2(sha256) -> list[DatasetHandle]`）
- Test: `tests/test_curation_on_dataset.py`、既存の `tests/test_curation_*.py`（`file_path` を渡していた呼び出しを、tmp の解析フォルダを `load_dataset` で開く形へ）

- [ ] **Step 1: 失敗するテストを書く**: `test_curation_review_reads_dataset_arf2`、`test_submit_refreshes_flags_on_open_handles`（Review Focus 2）、`test_curation_without_arf2_is_explained`（Review Focus 4: mzTab-M だけのデータセット）、`test_curation_uses_dataset_library`。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る。**
- [ ] **Step 5: コミット**: `feat(curation): キュレーションのツールをデータセットの上で動かす`

### Task 5: `dataset_feature_search` が `eic_search` と同定の標準化の列を吸収する

**Files:**
- Modify: `metabolomix/analysis/feature_search.py`、`metabolomix/tools/dataset_analysis_tools.py`（`dataset_feature_search`）
- Test: `tests/test_dataset_feature_search_tool.py`（追記）

- [ ] **Step 1: 失敗するテストを書く**: `test_sort_by_intensity_uses_detected_max`（gap-fill のセルが大きくても順位に効かない、除外した試料も効かない）、`test_intensity_sort_without_conditions_returns_top_n`、`test_mz_sort_without_conditions_still_refused`、`test_rows_carry_ontology_n_detected_identity`（解析フォルダのデータセットで `normalized` `msi_level` が埋まり、mzTab-M のデータセットでは空）。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る。**
- [ ] **Step 5: コミット**: `feat(tools): dataset_feature_search に強度順と同定の標準化の列を足す`

### Task 6: `dataset_plot`（volcano・pca・loadings・eic・eic_compounds）と図の置き場

**Files:**
- Create: `metabolomix/tools/plot_tools.py`（`dataset_plot`。`kind` ごとの関数へ振り分ける）
- Modify: `metabolomix/mztab/dataset_state.py`（`figures: dict[str, dict]`。指紋に含めない）、`metabolomix/tools/reports.py`（`save_figure(kind, analysis_id, title=None, result_id=None, dataset_handle=None)` は `ds.figures` と結果の登録から描く。`source` 引数と ARF の候補を消す）、`metabolomix/plots/pca_scores.py`（`lipidmix.pca_scores.v1` の payload）、`metabolomix/plots/pca_loadings.py`（`_mztab_features` 相当で `feature_metadata` の ontology を使う）、`metabolomix/plots/eic.py`
- Test: `tests/test_dataset_plot.py`

**Interfaces:**
- Produces: `dataset_plot(kind: str, dataset_handle=None, result_id=None, output=None, title=None, **kind ごとの引数) -> list | str`。FastMCP は `**kwargs` を受けないので、表（決めたこと 14）の引数をすべて省略可能な引数として並べ、`kind` に関係ない引数が渡されたら `INVALID_ARGUMENT`（どの `kind` の引数かを書く）。kind ごとの関数は `plot_tools._KINDS: dict[str, Callable]`。
- [ ] **Step 1: 失敗するテストを書く**: `test_volcano_from_dataset_differential`、`test_pca_color_by_batch`、`test_loadings_bar_colour_uses_ontology`、`test_eic_from_artifact_paths`、`test_eic_compounds_by_ontology`、`test_eic_without_artifact_is_explained`（Review Focus 4）、`test_payload_format_names_are_kept`（表の形式名）、`test_irrelevant_argument_rejected`、`test_save_figure_reads_dataset_figures`、`test_stale_result_not_drawn_as_current`（Review Focus 3 の図の側）。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る。**
- [ ] **Step 5: コミット**: `feat(plots): dataset_plot（volcano・pca・loadings・eic）と図の置き場`

### Task 7: `dataset_plot`（group_intensity・species・species_pca）

**Files:**
- Create: `metabolomix/plots/dataset_rows.py`
- Modify: `metabolomix/tools/plot_tools.py`
- Test: `tests/test_dataset_plot_species.py`（旧 `tests/test_group_intensity_tools.py` `tests/test_species_tools.py` の主張のうち、データの意味に関わるものを移す: キュレーションの当て方・検出限界・MS/MS の要求・share の分母）

**Interfaces:**
- Produces: `dataset_rows.select_items(ds, items: list[str], *, apply_curation: bool, exclude_auto_likely_wrong: bool, require_msms: bool) -> list[str]`（クラス名なら `feature_metadata["ontology"]`、特徴 ID ならそのまま。除外した特徴を除く）、`dataset_rows.rows_by_spot(ds, feature_ids, *, groups: list[str], low_reliability_samples, standard_samples) -> dict`（`plots/group_intensity.build_group_intensity_payload` と `plots/species.build_species_payload` が今受けている形。値は `feature_matrix`、gap-fill の印は `detected_mask`、群は `sample_meta` の group、除外した試料を除く）。`species_pca` は `analysis/pca.run_pca_fit_subset` を選んだ特徴の行列に当て、`ds.figures["species_pca"]` と結果の登録に置く。
- [ ] **Step 1: 失敗するテストを書く**（上の主張と、`test_gap_filled_cells_marked_in_group_intensity`）。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る。**
- [ ] **Step 5: コミット**: `feat(plots): dataset_plot に群ごとの強度・分子種・分子種の PCA を足す`

### Task 8: 旧ツールとスロットを消す

**Files:**
- Delete: `metabolomix/arf/tools.py` `metabolomix/arf/species_tools.py` `metabolomix/arf2/tools.py` `metabolomix/pai2/tools.py` `metabolomix/eic/tools.py` `metabolomix/tools/samples.py` `metabolomix/tools/pca_loadings_tools.py`、`dcl_parser`（`metabolomix/dcl/tools.py` の中）、旧ツールだけが使っていたモジュール（`metabolomix/arf/{selection,exclusions,identity_join,evidence_export}.py`、`metabolomix/msdial/peak_verification.py` など。消す前に `grep -rn "import <モジュール>\|from <パッケージ> import" metabolomix server.py` で生きている利用者が 0 であることを確かめ、一覧を HISTRY に書く）、それらのテスト
- Modify: `metabolomix/core/session_state.py`（`ArfState` `Arf2State` `Pai2State` `EicState`、`last_loadings_plot`、`_build_sample_meta` `_new_dataset_id`、それらのための import を消す）、`metabolomix/core/tool_helpers.py`（ARF だけの補助と `_remember_arf_pca_plot` `_pp_*` `_build_verification_dossier` を消す。`_identity_tables` `dataset_pca_plot` `_pca_scatter_arrays` は残す）、`metabolomix/core/handle_dispatch.py`（`SESSION_READERS` から消したツールを外す）、`metabolomix/pathway/tools.py`（`source="arf"` と `_build_arf` `_arf2_catalog`）、`metabolomix/arf/reader.py`（`run_pca` の再エクスポートを消す）、`tests/test_analysis_pca.py`（再エクスポートの主張を消す）、`server.py`（star import と再エクスポート）、`tests/test_missing_state_contract.py` `tests/test_mcp_errors.py` `tests/test_console_payload_trim.py` ほか、消したツール・補助を名指しするテスト
- Test: `tests/test_session_namespaces.py`（スロットが無いことの主張に書き直す）

- [ ] **Step 1: 失敗するテストを書く**: `test_session_has_no_parser_slots`（`AnalysisSession` に `arf` `arf2` `pai2` `eic` `last_loadings_plot` が無い）、`test_registered_tool_set_is_b2_surface`（`tests/test_server_registration.py` の `EXPECTED_TOOLS` を決めたこと 1 の 56 に）、`test_pathway_has_no_arf_source`。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る**（全体）。
- [ ] **Step 5: コミット**: `refactor: 旧 ARF / ARF2 / PAI2 / EIC のツールとセッションのスロットを消す`

### Task 9: `data_dir` の探索を消し、パスを省いたツールは直近のデータセットのファイルを使う

**Files:**
- Delete: `metabolomix/core/data_config.py`
- Modify: `metabolomix/core/mcp_core.py`（`DATA_DIR` を消す。`_report_dir_candidates(holder)` は決めたこと 3）、`metabolomix/core/path_resolvers.py`（`_resolve_data_file` の探索をやめ、パスが無ければ `None`。`list_data_files(directory)` は必須。`resolve_library_path` の探索と `ion_mode` を消す）、`metabolomix/tools/_holder.py`（`alignment_paths`: ハンドルを省いたら直近のハンドルのパス。直近が無ければ `missing_state`）、`metabolomix/dataset/` に `artifact_path(ds, role: str) -> Path | None`（`dcl` / `msms_evidence`、`eic` / `chromatogram`、`arf2` / `spot_catalog` を同じ役割として引く）、`metabolomix/dcl/tools.py`（`dcl_find_msms` の既定はデータセットの `.dcl`）、`metabolomix/library/tools.py`（`library_load` の決めたこと 5、`_spot_query` `_measured_spectrum` の既定）、`metabolomix/tools/reanalysis_tools.py`（`library_msms_qualify`）、`metabolomix/tools/reports.py` `metabolomix/tools/objective.py` `metabolomix/pathway/tools.py`（記録の置き場所。`dataset_handle` を足す）、`metabolomix/library/store.py`（`main` の CLI）、各 `reader.py` の CLI（`find_input_file` `_list_arf_files` `_select_arf_file` `_select_dcl_file` を消し `--file` 必須）、`metabolomix/console/job_manager.py`、`scripts/new-worktree.sh`（`LIPIDMIX_DATA_DIR` を書かない。`LIPIDMIX_CONFIG` だけ）、`CLAUDE.md` `README.md` `DEPLOY.md` `docs/cli.md` `docs/output_format/core.md` `docs/workflow/dataset.md` の `DATA_DIR` / `LIPIDMIX_DATA_DIR`、`mcp_core.DATA_DIR` を差し替えるテスト（約 27 ファイル。解析フォルダを tmp に作って `load_dataset(directory=...)` で開くか、パスを明示する形へ）、`tests/test_package_layout.py`（`DEFAULT_DATA_DIR`）、`tests/test_no_asset_env_vars.py`（`LIPIDMIX_DATA_DIR` を足す）
- Test: `tests/test_pathless_tools.py`

- [ ] **Step 1: 失敗するテストを書く**: `test_pathless_tools_use_current_dataset_files`（解析フォルダを開いた後、`dcl_spot_spectrum(spot_id=...)` `dcl_find_msms(precursor_mz=...)` `library_match_feature(spot_id=...)` がそのフォルダの `.arf2` / `.dcl` を使う）、`test_pathless_tools_do_not_fall_back_after_close`（Review Focus 1。`data/` に同名のファイルを置いても使わない）、`test_list_data_files_requires_directory`、`test_reports_go_to_analysis_dir`（`write_report` → `<analysis_dir>/reports/`、`[paths] reports_dir` があればそちら）、`test_library_load_without_memory_is_refused`（`LIBRARY_NOT_REMEMBERED`。設定の極性の MSP があっても使わない）、`test_library_load_has_no_ion_mode`。
- [ ] **Step 2〜4: 落ちる → 実装 → 通る**（全体）。
- [ ] **Step 5: コミット**: `refactor: data/ の探索を消し、パスを省いたツールは直近のデータセットのファイルを使う`

### Task 10: `MCP_INSTRUCTIONS` と注意の要約を書き直す

**Files:**
- Modify: `metabolomix/core/mcp_core.py`（`MCP_INSTRUCTIONS`、`OUTPUT_FORMAT_SECTIONS`）、`metabolomix/core/session_state.py`（`SEMANTICS_CAVEAT` の行を、`DatasetState` の特徴×試料の行・gap-fill・`.dcl` の MS/MS・EIC の `peak_top` に）、`metabolomix/tools/resources.py`（静的リソースの説明のトピックの一覧に `reanalysis` も）
- Test: `tests/test_output_format_sections.py`、`tests/test_semantics_caveat.py`（書き直し）

- [ ] **Step 1: 失敗するテストを書く**: `test_instructions_name_only_registered_tools`（`MCP_INSTRUCTIONS` の中の `` `名前(` `` と `` `名前` `` のうちツールの名前の形のものが、登録されたツールかワイルドカード `dataset_*` などに限られる）、`test_topics_are_b2_set`（`dataset` があり `arf` `arf2` `pai2` が無い）。
- [ ] **Step 2〜4: 落ちる → 実装（決めたこと 17）→ 通る。** トピックの文書の中身は Task 11。この Task では `docs/output_format/dataset.md` を見出しだけで作り、`arf.md` `arf2.md` `pai2.md` を消す前に Task 11 の写しの表を作っておく（中身を失わないため、Task 10 と 11 は続けて行う）。
- [ ] **Step 5: コミット**: `docs(instructions): MCP_INSTRUCTIONS と注意の要約を新しいツール面で書き直す`

### Task 11: 文書と腐敗防止のテスト

**Files:**
- Modify: `USAGE.md`（新しいツール面で書き直し。冒頭の流れも）、`docs/workflow/`（`arf.md` `arf2.md` `pai2.md` `eic.md` `mztab.md` を消し、`dataset.md` を入口・試料（`list_data_files` `load_dataset` `dataset_list` `dataset_close` `dataset_status` `dataset_exclude`）で、`plots.md` を `dataset_plot` `save_figure` で書き直す。`dataset_set_sample_metadata` は `dataset_analysis.md` のまま。`index.md` の目次と合計）、`tests/test_workflow_docs.py`（`IN_SCOPE` `OUT_OF_SCOPE`、合計の数）、`docs/output_format/`（`dataset.md` の本文、`eic.md` `core.md` `mztab.md` `library.md` `reanalysis.md` `curation.md` の旧ツール名）、`CLAUDE.md`（冒頭の規模表記、「このリポジトリは何か」の主経路の記述を全体 §2 の「2 つの入口を対等に」へ、構成の節の `session.arf` などのスロット分離の鉄則・`run_pca` の再エクスポートの鉄則・`ANALYSES_DIR` / `KNOWLEDGE_DIR` / `DATA_DIR` の記述・知識の蓄積の記述・ドキュメントの地図、worktree の節の蓄積状態）、`README.md` `DEPLOY.md`、`tests/test_tool_annotations.py`
- Create: `tests/test_removed_tools_not_mentioned.py`

- [ ] **Step 1: 失敗するテストを書く**: `test_removed_tools_are_not_mentioned`（Review Focus 5。`REMOVED_TOOLS` = 消した 39 の名前（旧 ARF / ARF2 / PAI2 / EIC / 試料 / 描画 / 文献の各ツール、`dataset_load` `dataset_handle_open` `dataset_handle_list` `dataset_handle_close` `update_objective` `list_reports`）が、`README.md` `USAGE.md` `DEPLOY.md` `CLAUDE.md` `docs/workflow/**` `docs/output_format/**` `docs/cli.md` と `metabolomix/**/*.py` の文字列（`required_tools`・案内・docstring）に単語として現れない。`docs/superpowers/**` と `docs/schema/**` は対象外）。
- [ ] **Step 2: 落ちることを確かめる。**
- [ ] **Step 3: 実装**（文書の書き直し。`docs/output_format/` は決めたこと 16 の写しの表どおりに、必須の注意を 1 行も落とさない）。
- [ ] **Step 4: 通ることを確かめる**（全体。`tests/test_readme_links.py` `tests/test_workflow_docs.py` `tests/test_output_format_sections.py` `tests/test_server_registration.py`）。
- [ ] **Step 5: コミット**: `docs: USAGE・workflow・output_format・CLAUDE.md を新しいツール面に合わせる`

### Task 12: 記録

**Files:**
- Modify（main ツリー側。追跡外）: `C:\Users\yuu18\Metabolomix_with_LLM\docs\HISTRY.md`（末尾に「2026-11-07 B-2b 完了」節）、`C:\Users\yuu18\Metabolomix_with_LLM\docs\task.md`（末尾への追記だけ）

- [ ] **Step 1: HISTRY**: 決めたこと（1〜19）、実装中の裁定、消したモジュールの一覧、`docs/output_format/` の必須の注意の移し先の表、失った機能（決めたこと 15）、**流れが変わる点**（vault の流れ図は親セッションが統合時に直すので列挙だけ）:
  - ENTRY POINT: `load_dataset` が唯一の入口（`dataset_load` `dataset_handle_open` は無い）。読み込み直後の概観と前処理前の PCA が無くなった。
  - 前提状態の連鎖: `missing_state` の `required_tools` が `load_dataset` に。パスを省いたツール（`dcl_*` `library_*` `curation_*` `write_report` `record_objective` `save_figure` `pathway_project`）は直近のデータセットを要求する。
  - 経路の増減: ARF の経路（`arf_*` → `arf_differential` → `arf_export_differential`）が消え、差次的解析は `dataset_*` の経路だけ。文献の発見の経路（`knowledge_coverage` → `paper_search` → `ingest_*`）が消えた。
  - 既定値: 記録は解析フォルダの `reports/`（`data/reports` ではない）。`library_load` は覚えたライブラリだけ（極性の設定に落ちない）。
  - 外部への影響（後から追従）: Use-LLLM のツールカタログ・massbank-context の文書・reanalysis-study の手順書（`dataset_handle_open` → `load_dataset`、`dataset_handle_list` → `dataset_list`、`dataset_handle_close` → `dataset_close`）。
- [ ] **Step 2: task.md の末尾に追記**: B-2 で閉じた項目、決めたこと 15 の候補、C・D へ渡すもの（`dataset_differential` `dataset_statistic` `dataset_build_matrix` の削除は C、旧 Console 6 つは D、`dataset_pca` の群の指定は C）。
- [ ] **Step 3: 全テスト**: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q` の結果を HISTRY に書く。
