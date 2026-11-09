# A 最初の概観（`dataset_overview`）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** データセットを開いた直後に、品質 → 全体 → クラス → 分子種の 4 段で「群の間で何が変わり、何が変わっていないか」を粗く示し、
各段に掘り下げの次の一手を添える。`load_dataset` は読み込みの最後にこれを 1 回呼び、鎖を進めた後は `dataset_overview` で出し直せる。

**Architecture:** 3 層に分ける。(1) 品質の純計算 `metabolomix/analysis/stats/quality.py`（StudyData / PreparedMatrix だけを読む。
セッションにも MCP にも依存しない）。(2) 事実の組み立て `metabolomix/analysis/overview.py`（ds と計算結果から段ごとの事実・判定・次の一手を作る。
`stats_results` と同じ位置づけ）。(3) ツール層 `metabolomix/tools/overview_tools.py`（既定の前処理の登録・PCA・比較・組成を呼び、結果を登録し、
`reports/stats/` に記録を 1 件書き、要約を返す）。統計はすべて C のエンジン（`stats.preprocess` / `stats.pca` / `spec_v2.run_spec_v2` /
`stats.composition`）を通し、概観は新しい統計を発明しない。

**Tech Stack:** Python 3.14、numpy、scipy、FastMCP。依存は足さない。

**Spec:** `docs/superpowers/specs/2026-11-06-dataset-overview-design.md`（A。以下 §n はこの spec）。全体は `2026-11-06-overhaul-overview-design.md`、
C の公開面は `docs/superpowers/plans/2026-11-08-c-wiring.md` と共有ノート `c-wiring-handover.md`（結果の形・`as_two_group`・`needs_input`）。

## Global Constraints

- ツールは 58 → **59**（`dataset_overview` を足す）。リソース 2・リソーステンプレート 1 は変えない。足すコミットで `tests/test_server_registration.py`
  （件数と `ToolAnnotations`）・`USAGE.md`（一覧と宣言件数）・CLAUDE.md 冒頭の規模表記・`docs/workflow/`（`IN_SCOPE`）を同時に直す（C 配線の R5）。
- `export_contract.EXPORT_COLUMNS` と `CONTRACT_VERSION = 1` は変えない。根拠バンドルは v3 のまま。
- `metabolomix/analysis/stats/` は `tools`・`core.session_state`・`core.mcp_core` を import しない。`analysis/overview.py` も `tools.*` を import しない。
- 新しいツールは `dataset_handle` を受ける（`handle_dispatch` が包む。`SESSION_READERS` は空のまま）。`structured_output=False`。
- 戻り値は `json_payload`、float は丸め、行の一覧は TSV（列名 1 回）。非有限は出さない（`json_safe`）。最初の応答の目安は **1,500 トークン**（§10）。
- 判定は C の `verdict`（changed / unchanged / inconclusive）だけを使う。q や log2fc から自分で「変わっていない」を作らない。有意の旗は行の `significant`。
- 群の名前から対照を推測しない。2 群の向きは仮で、仮であることを戻り値に書く（§5）。
- 図は最初の応答に入れない。`dataset_plot` の `kind` と引数を次の一手に書く（§6）。
- テストは main ツリーの `.venv`: `C:/Users/yuu18/Metabolomix_with_LLM/.venv/Scripts/python.exe -m pytest tests -q`。
  コミットは pre-commit が全テストを回すのでバックグラウンドで実行し、出力をファイルへ落とす。`git stash` は使わない。
- 研究室の参照ライブラリのパスを追跡対象に書かない。文書・docstring は日本語。fixture はテスト自身が作る（`tests.test_stats_results._ds` を使う）。

## Review Focus

1. **群が 1 つ・群の宣言が無い（Console 由来で Class ID が空など）**: 概観は失敗せず、品質と全体（PCA）だけを出し、比較の段は
   「群が 1 つなので比較しない」と書いて `dataset_set_sample_metadata` を次の一手にする。→ Task 4 `test_overview_single_group_skips_comparisons`
2. **群とバッチが完全交絡**: 比較を回さない（`allow_confounded` を渡さない）。品質の段の判定が `poor` で、交絡を事実として出す。
   → Task 4 `test_overview_complete_confounding_does_not_compare`
3. **換算済み（`quant_state.normalized`）で QC がある**: `needs_input` を返さず、ドリフト補正なしで計算し、その旨と `dataset_preprocess(drift_mode=...)` を
   次の一手に書く（§7）。→ Task 4 `test_overview_drift_undecided_uses_off_and_says_so`
4. **`allow_incomplete=True` で開いた探索専用のデータセット**: 比較の段を出さない（`dataset_compare` が拒むのと同じ理由）。概観で `load_dataset` が失敗しない。
   → Task 5 `test_load_dataset_incomplete_overview_has_no_comparisons`
5. **概観の途中で例外（行列が作れない・特徴が残らない）**: `load_dataset` は成功し、`overview` に `status:"unavailable"` と理由と次の一手が出る。
   → Task 5 `test_load_dataset_survives_overview_failure`

---

## 判断が要る点（承認のときに決める。推奨を先に書く）

| # | 点 | 推奨 | ほかの案 |
|---|---|---|---|
| 1 | ツールの形 | 新ツール `dataset_overview(top_n=5, dataset_handle=None)`。`load_dataset` は最後にこれと同じ組み立て（`build_overview_payload`）を呼び、戻り値の `overview` 欄に入れる（§2） | `load_dataset` の戻り値だけ（鎖の後に出し直す口が無い） |
| 2 | `load_dataset` で概観を止める口 | `load_dataset(overview=True)` を足す。`False` で今の戻り値のまま（再解析のエージェントが多数のハンドルを開くとき、時間と文脈を食わない） | 止められない（いつも計算） |
| 3 | 概観の比較を「直近」にするか | `ds.results` には登録する（volcano・特徴検索・書き出しが `result_id` で掘れる）が、**`ds.last_differential` / `ds.last_pca` は動かさない**（`register_compare` / `register_pca` に `make_current=False`）。利用者が明示に回した比較を中立の概観で上書きしない。次の一手には `result_id` を必ず書く | 直近にする（素の `dataset_plot(kind="volcano")` が仮の向きの比較を描く） |
| 4 | 2 群の仮の向き | 試料の表（`ds.sample_names` の順）で**先に現れた群を基準**。`direction: "provisional"` と決め方を書き、`dataset_compare(reference_group=…, test_group=…)` を次の一手に | 名前の辞書順 |
| 5 | 3 群以上の対の表 | 群が **6 つまで**は全対（最大 15 対）を 1 回の `run_spec_v2` で比べ（moderated F は同じ当てはめから出る）、F の q ≤ 0.05 の特徴について対ごとの changed の件数を TSV にする。7 群以上は F の件数だけを出し、対の表は省いて理由を書く | 全対をいつも（28 対以上は文脈を食う） |
| 6 | 品質の判定と差の段の関係 | 品質の段に `verdict`（`ok` / `caution` / `poor` / `not_evaluable`）と理由の事実 ID を置く。`poor` でも差の段は計算するが、各段の頭に「品質 Q… のため参考」と書く（ただし完全交絡なら比較そのものを回さない）。判定の規則は下の「品質の判定」 | `poor` なら差の段を出さない（「何が悪いか」以外が何も分からない） |
| 7 | クラスの段の門 | spec §4.3 のとおり、`quant_state.steps` に `"select"` が無ければ出さない。mzTab-M 由来（定量の鎖が使えない）も同じく省き、理由と次の一手（`dataset_select_features` など、使えないなら `dataset_analysis_run` の組成）を書く | 門を外す引数を足す（見せかけの `unchanged` を出しうる。反目標） |
| 8 | 事実の ID と保存 | 段ごとの通し番号 `Q1` `G1` `C1` `S1`。外から引くときは `ovw_<hash>#Q1`（`result_id#contrast_id` と同じ形）。概観は `ds.results[ovw_…]`（`kind:"overview"`）と `reports/stats/<record_id>.json`（`step:"overview"`、親は前処理の記録）に構造化して残す。戻り値は 1 事実 1 文字列 | 事実を dict で返す（キーが文脈を食う） |
| 9 | 範囲の外 | (a) `reports/stats/` からの復元（task.md の判断表 #10。概観は決定論なので出し直せば同じ）、(b) `record_objective` の問いに合わせた並べ替え（§11 は任意）、(c) 2 因子の構造（§5、利用者が因子を指定したときだけ → `dataset_compare`） | (b) を入れる |

### 品質の判定（#6 の規則）

- `poor`: 群とバッチが完全交絡／失敗した QC がある／QC の RSD の中央値 > 0.30。
- `caution`: 部分交絡（Cramér の V ≥ 0.5 か注入順の Kruskal p < 0.05）／ドリフトの兆候（QC が 5 本以上のバッチで、QC の総量の log2 と注入順の
  Spearman の |ρ| ≥ 0.6 かつ p < 0.05）／反復注入らしい名前があるのに `subject_id` が無い（`replicate_name_warnings`）。
- `not_evaluable`: QC が 2 本未満で、かつ交絡も評価できない（バッチ不明）。QC だけが無いときは QC の項目を「評価できない」と書き、他の項目で判定する。
- それ以外は `ok`。閾値は `quality.py` の定数に置き、記録の `settings` に書く。

---

## 戻り値の形（`dataset_overview` と `load_dataset` の `overview` 欄）

```
{"status": "success", "overview_id": "ovw_…", "provisional": true,
 "state": "<1 行: 入口・極性・assay_kind・鎖（キュレーションの判断 N 件／換算・統合・選別の済み）・recipe（matrix_id・drift・足切り・正規化）>",
 "counts": {...},            # 単独で呼んだときだけ（load_dataset の summary と重なるため、埋め込みでは省く）
 "quality": {"verdict": "caution", "facts": ["Q1 …", "Q2 …"], "next": "…"},
 "global":  {"facts": ["G1 …", …], "next": "…"},
 "classes": {"facts": [...], "table": "<TSV>", "next": "…"} | {"omitted": "<理由>", "next": "…"},
 "species": {"facts": [...], "top": "<TSV>", "next": "…"},
 "result_ids": {"compare": "cmp_…", "pca": "pca_…"},
 "record": {...},
 "note": "考察は事実 ID を挙げて事実と別に書く（サーバは考察を書かない）"}
```

- 事実は「条件付きの数字」の 1 文（例 `S1 changed 12（上昇 9・下降 3）／unchanged 140／inconclusive 51（moderated t、q≤0.05、同等 ±1.2 倍、T vs C は仮の向き）`）。
  条件が段で共通なら段の最初の事実にだけ書く。
- 各段の `next` は 1 行（ツール名と引数。§9）。品質: `dataset_plot(kind="pca", result_id="pca_…", color_by="batch")`・`dataset_exclude(samples=[…])`。
  全体: `dataset_compare(reference_group=…, test_group=…)` か `dataset_compare(groups=[…])`。クラス: `dataset_plot(kind="group_intensity", items=[…])`。
  分子種: `dataset_feature_search(result_ids=["cmp_…#<cid>"], feature_ids=[…])`・`dataset_plot(kind="volcano", result_id="cmp_…", contrast_id=…)`。
- `species.top` は changed の上位 `top_n`（既定 5、上限 20）を `contrast_id feature_id name id_level log2fc ci_low ci_high q`。3 群以上は F の q の小さい順に
  `feature_id name id_level f_q` と、対の表 `pair n_changed n_up n_down`。
- `classes.table` は変わったクラスと変わっていないクラスの両方（各最大 5）を `class n_species log2fc_total ci_low ci_high verdict_total clr_log2fc verdict_clr`。
- `id_level` は同定の証拠の水準（`msms`=SME の候補あり / `ms1`=MS1 注釈のみ / `none`。`load_dataset` の summary と同じ 3 分類）。

---

## ファイル構成

```
metabolomix/analysis/stats/quality.py      新規。品質の純計算（QC RSD・失敗 QC・検出率と gap-fill・ドリフトの兆候・判定）
metabolomix/analysis/overview.py           新規。段ごとの事実・判定・次の一手・overview の登録（ds.results）
metabolomix/analysis/stats_results.py      変更。register_compare / register_pca に make_current、compare_spec（stats_tools から移す）
metabolomix/analysis/stats_records.py      変更。STEPS に "overview"
metabolomix/tools/overview_tools.py        新規。dataset_overview と build_overview_payload
metabolomix/tools/stats_tools.py           変更。_spec_for / _contrast_id を stats_results.compare_spec へ（振る舞いは同じ）
metabolomix/tools/handle_tools.py          変更。load_dataset が overview を埋め込む（overview 引数・_NEXT_STEPS・docstring）
metabolomix/core/mcp_core.py               変更。MCP_INSTRUCTIONS の ENTRY POINT と「考察は事実 ID を挙げて別に」
server.py                                  変更。overview_tools の import と再エクスポート
tests/test_stats_quality.py / test_overview_facts.py / test_overview_tool.py / test_load_dataset_overview.py   新規
文書: USAGE.md・CLAUDE.md・docs/workflow/dataset_analysis.md（＋ IN_SCOPE）・docs/output_format/dataset.md（概観の節）・vault の流れ図
```

---

### Task 1: 品質の純計算（`stats/quality.py`）

**Files:** Create `metabolomix/analysis/stats/quality.py`、Test `tests/test_stats_quality.py`

**Interfaces:**
- Produces:
  - `QUALITY_THRESHOLDS: dict`（`qc_rsd_ok=0.30`、`failed_qc_ratio=0.2`、`drift_min_qc=5`、`drift_rho=0.6`、`drift_p=0.05`、`cramers_v=0.5`、`kruskal_p=0.05`）
  - `qc_rsd_summary(study: StudyData) -> dict` — QC の行（`role=="qc"` かつ `include`）の前処理前の値で特徴ごとの RSD（有限正値、ddof=1。`preprocess._sd_mean` と同じ定義）。
    `{"status": "ok"|"not_evaluable", "n_qc", "n_computable", "median", "frac_le_030", "reason"?}`。QC < 2 は `not_evaluable`（`reason:"too_few_qc"`）。
  - `failed_qc(study) -> dict` — QC の総量（有限正値の和）が QC の中央値の 0.2 倍未満の注入。QC < 3 は `{"checked": n, "failed": []}`（`preprocessing.detect_failed_qc` と同じ規則。v1 のモジュールは D で消えるので import しない）。
  - `detection_summary(study) -> dict` — 生体試料のセルの実検出率と gap-fill の割合（`cell_states` の DETECTED / 値ありだが未検出 / 欠測）。`study.detected is None` なら `{"status":"not_evaluable"}`。
  - `drift_signs(study) -> list[dict]` — バッチごと（`preprocess._batch_key` と同じ鍵）に、QC が 5 本以上なら QC の総量の log2 と注入順の Spearman（`scipy.stats.spearmanr`）。`{"batch", "n_qc", "rho", "p", "flag": bool}`。
  - `quality_verdict(*, rsd, failed, drift, confounding: dict, replicate_warnings: list) -> tuple[str, list[str]]` — 上の「品質の判定」。理由はコード（`"complete_confounding"` など）の一覧。

- [ ] **Step 1: 失敗するテストを書く**（`tests.test_stats_results._ds` から `from_dataset` で StudyData を作る）
  - `test_qc_rsd_summary_matches_definition` — QC 4 本の合成で median と `frac_le_030` が手計算と一致。
  - `test_qc_rsd_not_evaluable_with_one_qc` — `n_qc=1` で `status=="not_evaluable"`、`reason=="too_few_qc"`。
  - `test_failed_qc_flags_low_total` — QC の 1 本を全特徴 0.1 倍にすると `failed == ["QC2"]`（試料 ID で返す）。
  - `test_detection_summary_without_mask` — `detected=None` で `not_evaluable`。
  - `test_drift_signs_detects_monotone_qc` — QC 6 本の総量を注入順に単調増加させると `flag is True`、乱数だけなら `False`。
  - `test_quality_verdict_rules` — 完全交絡 → `poor`／部分交絡だけ → `caution`／何も無い → `ok`／QC 無し・バッチ不明 → `not_evaluable`。
- [ ] **Step 2: 走らせて失敗を確かめる** — `…/python.exe -m pytest tests/test_stats_quality.py -q` → ImportError。
- [ ] **Step 3: 実装** — numpy と `scipy.stats.spearmanr` だけ。`stats.data.cell_states` と定数を使う。
- [ ] **Step 4: 通す** — 同じコマンドで PASS。`tests/test_stats_package.py`（import の規則）も通ること。
- [ ] **Step 5: コミット** — `feat(stats): 品質の純計算（QC RSD・失敗 QC・検出率・ドリフトの兆候・判定）`

### Task 2: C の部品の小さな受け口（`make_current`・`compare_spec`・記録の `overview`）

**Files:** Modify `metabolomix/analysis/stats_results.py`・`metabolomix/tools/stats_tools.py`・`metabolomix/analysis/stats_records.py`、Test `tests/test_stats_results.py`・`tests/test_stats_records.py`

**Interfaces:**
- Produces:
  - `register_compare(ds, run_result, *, log2fc_threshold, units, make_current: bool = True) -> dict` — `False` なら `ds.last_differential` を動かさない。
  - `register_pca(ds, pca_result, groups, *, make_current: bool = True) -> dict` — `False` なら `ds.last_pca` を動かさない。
  - `compare_spec(*, factors, comparisons, covariates, block, method, effect_size, bh_population, q_threshold, equivalence_fc, allow_confounded, preprocess, inputs) -> dict`
    （今の `stats_tools._spec_for` をそのまま移す。`stats_tools` はこれを呼ぶ）と `contrast_id_for(test, reference, taken, *, test_pos, ref_pos) -> str`（今の `_contrast_id`）。
  - `stats_records.STEPS = ("preprocess", "compare", "pca", "overview")`。

- [ ] **Step 1: テスト** — `test_register_compare_make_current_false_keeps_last`（先に登録した比較が `last_differential` のまま、新しい記録は `ds.results` にある）、
  同じく PCA、`test_write_stats_record_accepts_overview_step`。
- [ ] **Step 2: 失敗を確かめる。**
- [ ] **Step 3: 実装。** 移した関数の振る舞いは変えない（`tests/test_stats_tools_compare.py` がそのまま通る）。
- [ ] **Step 4: `tests/test_stats_results.py tests/test_stats_records.py tests/test_stats_tools_compare.py` を通す。**
- [ ] **Step 5: コミット** — `refactor(stats): 概観のための受け口（make_current・compare_spec・記録の overview）`

### Task 3: 事実の組み立て（`analysis/overview.py`）

**Files:** Create `metabolomix/analysis/overview.py`、Test `tests/test_overview_facts.py`

**Interfaces:**
- Consumes: Task 1 の関数、`stats_results.as_two_group` は使わない（3 群以上も扱うので compare の記録の `engine.rows` / `engine.omnibus` を直に読む——
  これは C の「消費側は `as_two_group` で読む」の例外で、概観は**件数**を数えるだけで 2 群の形を要らないため。docstring に書く）。
- Produces:
  - `OVERVIEW_SCHEMA = "lipidmix.overview.v1"`
  - `group_plan(ds, prepared) -> dict` — `{"mode": "none"|"two"|"multi", "levels": [...], "reference": str|None, "direction": "provisional"|None, "reason"?}`。
    群は `prepared.units` の `group`（QC・ブランクは入らない）。順序は `ds.sample_names` で先に現れた順（判断 #4）。群が 1 つ以下は `none`。
  - `group_comparisons(plan) -> list[dict]` — `two`: 1 対比。`multi`: 6 群まで全対、7 群以上は基準との対比だけ（omnibus のため）と `pairs_omitted=True`。対比 ID は `contrast_id_for`。
  - `chain_state(ds, *, curation_judgments: int | None) -> dict` — `quant_state.steps` から `{"curation": n|None, "normalize", "merge", "select": bool, "provisional": not select, "pending": [...]}`。
  - `quality_facts(...)`, `global_facts(pca_result, labels, compare_record | None, plan)`, `class_facts(composition_result | None, *, gate: str | None, assay_kind)`,
    `species_facts(ds, compare_record, plan, top_n)` — それぞれ `{"facts": [{"id","text","data"}], "table"?/"top"?, "next": str, "omitted"?}`。
    `text` は戻り値の 1 文、`data` は記録に残す構造（数字・条件）。
  - 全体の要約の数字（§12 で任されたもの）: 主成分 1・2 の寄与率、群ごとの重心（PC1, PC2）、群で説明される分散の割合 η²（PC1・PC2 それぞれ。一元配置の平方和比）。
  - `overview_id(*, dataset_hash, metadata_revision, matrix_id, compare_id, pca_id) -> str` — `"ovw_" + canonical_hash(...)[:12]`。
  - `register_overview(ds, payload_record: dict) -> dict` — `ds.results[overview_id] = {"kind": "overview", "schema", "result_id", "facts": {段: [data…]}, "provenance": …}`。
    古さの判定は他の結果と同じ `provenance`（`dataset_hash` と `metadata_revision`）。

- [ ] **Step 1: テスト**（compare の記録は `tests.test_stats_results._ds` → 既定の前処理 → `run_spec_v2` → `register_compare` で作る）
  - `test_group_plan_two_groups_reference_is_first_seen` — 試料の表で T が先なら基準は T、`direction=="provisional"`。
  - `test_group_plan_single_group_is_none` / `test_group_comparisons_caps_pairs_at_six_groups`（7 群で `pairs_omitted`）。
  - `test_species_facts_counts_verdicts_and_directions` — 合成で f0〜f2 が 4 倍の群について changed の上昇 3、`top` の列名が上の定義どおり。
  - `test_species_facts_never_turns_not_significant_into_unchanged` — `inconclusive` の行が `unchanged` の件数に入らない。
  - `test_species_facts_reports_presence_only_count` — 出方だけが違う特徴（`presence_q` ≤ 閾値で量の判定が `inconclusive`）の件数が事実に出る（§4.4。検出の別が無ければ「評価できない」）。
  - `test_multi_group_pair_table_counts_only_f_significant` — 3 群で、F の q > 0.05 の特徴は対の表に入らない。
  - `test_class_facts_omitted_before_select` — `quant_state=None` で `omitted` に理由、`next` に `dataset_select_features`。
  - `test_class_facts_lists_changed_and_unchanged` — `steps` に `select` を持たせた合成で、両方の行がある。
  - `test_global_facts_eta_squared` — 群で完全に分かれる合成で PC1 の η² > 0.9。
  - `test_fact_ids_are_sequential_per_stage`（`Q1, Q2…`／`S1…`）。
- [ ] **Step 2: 失敗を確かめる。**
- [ ] **Step 3: 実装。** 文は短く、数字は 3 桁に丸める。各段の頭の注意（品質 `poor` のとき「品質 Q… のため参考」、仮の向き）は段の最初の事実に入れる。
- [ ] **Step 4: 通す。**
- [ ] **Step 5: コミット** — `feat(overview): 段ごとの事実・判定・次の一手`

### Task 4: `dataset_overview` ツール

**Files:** Create `metabolomix/tools/overview_tools.py`、Modify `server.py`・`metabolomix/core/mcp_core.py`・`tests/test_server_registration.py`・`USAGE.md`・`CLAUDE.md`（冒頭の規模）・
`docs/workflow/dataset_analysis.md`・`tests/test_workflow_docs.py`（`IN_SCOPE`）、Test `tests/test_overview_tool.py`

**Interfaces:**
- Consumes: Task 1〜3、`stats_results.matrix_inputs` / `prepare_default`、`tools/_matrix_store.register_default`、`stats.pca.pca`、`spec_v2.run_spec_v2`、`stats.composition.composition`、
  `curation.apply.flags_for_arf2`（`holder.arf2_path` があるときだけ。無ければ判断の件数は `None`）。
- Produces:
  - `build_overview_payload(holder, *, top_n: int = 5, embedded: bool = False) -> dict` — 戻り値の形（上）。`DomainError` は送出する（呼び出し側が封筒にする）。
  - `dataset_overview(top_n: int = 5, dataset_handle: str | None = None) -> str` — `ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True)`
    （結果と記録を登録するが、同じ入力なら同じ ID）。

**流れ（`build_overview_payload`）**
1. 既定の前処理: `matrix_inputs(ds, None)` → `from_dataset` → `prepare_default`。`NeedsInput` なら設定の `drift.mode` を `"off"` にして作り直し、
   状態の行と品質の段に「ドリフト補正の既定が未決のため補正なしで計算」と `dataset_preprocess(drift_mode=…)` を書く（§7）。`register_default` で登録（前処理の記録 1 件）。
2. 品質: Task 1 の関数（前処理前の StudyData）＋ `confounding_report(prepared.units, levels)`＋ `prepared.warnings` の反復注入の名前。
3. 全体: `stats.pca.pca(prepared)`（既定の scaling・impute）→ `register_pca(..., make_current=False)`。
4. 比較（`group_plan` が `none` でなく、`ds.exploratory_only` でなく、交絡が complete でないとき）: `compare_spec` で spec を作り `run_spec_v2(spec, study, prepared=prepared)` →
   `register_compare(..., make_current=False)`。`STATISTIC_CONFOUNDED` などの `DomainError` は段の `omitted` に写す（概観全体は止めない）。
5. クラス: 門（`"select" in steps`）を通ったときだけ `composition(prepared, assay_kind=session.assay_kind, factors=["group"], comparisons=…)`。
6. `register_overview` と `stats_records.write_stats_record(reports_dir, "overview", parent_ids=[前処理の記録], settings={thresholds, top_n, plan}, summary=段ごとの data)`。
   書けないときは `record: {written:false, reason}`（C と同じ）。

- [ ] **Step 1: テスト**（`tests/test_stats_tools_compare.py` の `_install` と同じ流儀で tmp の解析フォルダに `reports/` を置く）
  - `test_overview_two_groups_shape_and_ids` — 4 段の鍵、`overview_id` が `ovw_`、`result_ids.compare` が `ds.results` にある、`ds.last_differential` は `None` のまま。
  - `test_overview_is_deterministic` — 2 回呼んで `overview_id`・`result_ids`・事実の文が同じ。
  - `test_overview_single_group_skips_comparisons`（Review Focus 1）。
  - `test_overview_complete_confounding_does_not_compare`（Review Focus 2）— `quality.verdict == "poor"`、`species` が `omitted`、`ds.results` に compare が無い。
  - `test_overview_drift_undecided_uses_off_and_says_so`（Review Focus 3）— `quant_state={"normalized":True,"method":"internal_standard",…}` と QC で、`needs_input` の鍵が無く、
    状態の行に drift off、品質の `next` に `dataset_preprocess(drift_mode=`。
  - `test_overview_payload_size_budget` — 3 群・クラスあり（門を通した合成）で `len(json_payload(...)) <= 6000` 文字（1,500 トークンの代理。実データの実測は HISTRY に書く）。
  - `test_overview_strict_json` — `json.loads(..., parse_constant=<拒む>)` で読める。
  - `test_overview_record_written_and_cites_preprocess` — `reports/stats/` に `step=="overview"` の記録、`parent_ids` に前処理の記録。
  - `test_overview_with_dataset_handle_uses_that_handle` — 2 つ開いて、直近でない方を名指しして直近の ds に結果が載らない。
  - `test_overview_after_exclusion_recomputes` — `dataset_exclude` の後は別の `overview_id`。
- [ ] **Step 2: 失敗を確かめる。**
- [ ] **Step 3: 実装。** `server.py` に import と再エクスポート、`MCP_INSTRUCTIONS` に 1 段落（下の文面の要点: 概観は粗い地図で、判定は C の 3 値、事実は ID で引き、
  考察は事実 ID を挙げて事実と別の欄に書く。サーバは考察を書かない）。登録テストを 59 に、`USAGE.md`・CLAUDE.md 冒頭・`docs/workflow/dataset_analysis.md` に `dataset_overview` の節。
- [ ] **Step 4: 全テストを回す**（`tests/test_readme_links.py` `tests/test_workflow_docs.py` `tests/test_instructions_text.py` `tests/test_server_registration.py` を含む）。
- [ ] **Step 5: コミット** — `feat: dataset_overview（品質→全体→クラス→分子種の最初の概観。59 ツール）`

### Task 5: `load_dataset` に概観を埋め込む

**Files:** Modify `metabolomix/tools/handle_tools.py`・`metabolomix/core/mcp_core.py`（ENTRY POINT の「no overview or PCA is run for you」を直す）、Test `tests/test_load_dataset_overview.py`

**Interfaces:**
- Consumes: `build_overview_payload(holder, embedded=True)`。
- Produces: `load_dataset(..., overview: bool = True)`。成功の戻り値の `summary` の後に `overview` 欄。`overview=False` のときは欄ごと出さない（今の形）。
  概観は新しいハンドルの `handle.lock` の下で作る（`load_dataset` は `handle_dispatch` に包まれないため、直後の別の呼び出しと重ならないように）。
  `DomainError` やその他の例外は `{"status":"unavailable","error":{"code","message"},"next":"dataset_overview() を呼び直す / dataset_preprocess で原因を見る"}` に写し、`load_dataset` は成功のまま。
  `_NEXT_STEPS` は概観を前提に書き直す（`dataset_status` と `curation_review` は残し、3 つ目を「概観の各段の next から掘る」へ）。

- [ ] **Step 1: テスト** — `test_load_dataset_includes_overview`（tmp に作った解析フォルダ。既存の load_dataset のテストの fixture を使う）／
  `test_load_dataset_overview_false_keeps_old_shape`／`test_load_dataset_incomplete_overview_has_no_comparisons`（Review Focus 4）／
  `test_load_dataset_survives_overview_failure`（Review Focus 5。`build_overview_payload` を monkeypatch で例外にする）。
- [ ] **Step 2: 失敗を確かめる。**
- [ ] **Step 3: 実装**（docstring の「`.arf2` の概観や前処理前の PCA は返さない」を「最初の概観（`overview`）を返す」に）。
- [ ] **Step 4: 全テスト。** 既存の `load_dataset` の戻り値を固定したテストが落ちたら、概観の欄を足しただけの差かを確かめて直す（意味を変えたものは HISTRY に書く）。
- [ ] **Step 5: コミット** — `feat: load_dataset が最初の概観を返す（overview=False で止める）`

### Task 6: 文書・記録・流れ図・実データでの大きさの確認

**Files:** `docs/output_format/dataset.md`（概観の節: 段・事実 ID・判定の 3 値・門・仮の向き・品質の判定の規則）、`USAGE.md` の `load_dataset` の行、
vault の流れ図（ENTRY POINT の説明と経路 B: 読み込み → 概観（暫定）→ キュレーション…→ 概観の出し直し）、main ツリーの `docs/HISTRY.md`（追記）と `docs/task.md`（末尾に追記だけ）。

- [ ] **Step 1:** 文書を書く。`tests/test_readme_links.py`（数量表現の禁止を含む）と `tests/test_removed_tools_not_mentioned.py` が通ること。
- [ ] **Step 2:** 利用者の手元の実データがあれば、`load_dataset` の `overview` の文字数・トークンの目安と所要時間を測り HISTRY に書く（無ければ task.md に「実データで測る」を残す）。
- [ ] **Step 3:** 全テストを回し、件数を HISTRY に書く。
- [ ] **Step 4:** コミット — `docs: 最初の概観（dataset_overview）`。統合ブランチへのマージと push はしない（親セッションが行う）。

---

## task.md に残すもの（範囲の外）

- 判断表 #10 の復元（`reports/stats/` の overview の記録も対象に含める）。
- `record_objective` の群・クラスに合わせた次の一手の並べ替え（§11、任意）。
- 実データでの概観の大きさ・判定の閾値（RSD 0.30・ドリフトの |ρ| 0.6・V 0.5）の見直し。
- mzTab-M 由来のデータセットでクラスの段が出ない件（定量の鎖が `QUANT_SOURCE_UNSUPPORTED` のため。B-3 の残りと一緒に）。
