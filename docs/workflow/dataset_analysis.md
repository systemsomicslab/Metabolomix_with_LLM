# DatasetState 解析ツール呼び出し連鎖

解析フォルダ（`.arf` / `.arf2`）も mzTab-M も、`load_dataset` が同じ `DatasetState` にして
このツール群へ渡す。旧 ARF 経路（`arf_*` ツールとそのセッションスロット）は B-2b Task 8 で消えた。

## dataset_preprocess

統計エンジンの前処理（spec 2026-11-06 §4）。行列 ID を返し、`ds.stat_matrices[matrix_id]` に置く。除外（`ds.excluded_*`）と群の上書き（`group_levels` / `group_factors`）は
手順 3 が入力へ写す。手順 6 は問いを返さない: `drift_mode="auto"` は、補正できても内部標準で換算済み・注入順が
ファイルの順かもしれないときは当てずに計算し（`_auto_hold`）、ドリフトの段に理由と `rerun_with`（`drift_mode="on"`）を、
警告に `drift_not_applied_by_default` を載せる。正規化の係数は `normalization_features` の部分集合からも求められる
（`_factor_columns`。除いた特徴も行列に残す）。書いても効かない任意のキー（`normalize="none"` の係数の特徴・空の
`normalization_features`・`drift_mode="off"` の `loess_degree`）は手順 2 が誤りにせず設定から捨て（`_drop_inert`。
`matrix_id` は渡さないときと同じ）、手順 6 が生の設定を正規化し直すときにその注記を `warnings` に載せる。手順 28 の記録は `reports/stats/` への追記で、書き先が決まらない・
書けないときも前処理は成功する。

陽陰を統合したデータセット（`dataset_merge_polarities`。すべての特徴が `origin.polarity` を持つ）は、手順 5 が特徴 → 極性と、
統合が `quant_state.polarity_runs` に残した陰イオンの測定（バッチ・注入順）を `StudyData` に載せ、手順 8 の
`_polarity_blocks` が極性の塊（列と、その極性の測定としての試料の表 `polarity_samples`）に分ける。ドリフト補正
（手順 16）と正規化（手順 22）は塊ごとに同じ `_drift_step` / `_normalize_step` を当て、段に `by_polarity` を残す
（各極性の値は、その極性だけを前処理して並べたものと同じ）。足切り（手順 19）は特徴ごとの計算なので塊に分けない。
統合したデータセットの試料の表は陽イオンのもの（統合が陽イオンに当てたシートごと引き継ぐ）なので、陽イオンの塊は
陽イオンだけのデータセットと同じ試料の表を読む。

1. metabolomix/tools/stats_tools.py  dataset_preprocess()
2. └─ metabolomix/analysis/stats/preprocess.py  normalize_preprocess_config()
3. └─ metabolomix/analysis/stats_results.py  study_inputs()
4. └─ metabolomix/analysis/stats/data.py  from_dataset()
5.    └─ metabolomix/analysis/stats/data.py  merged_feature_polarity()（統合したデータセットだけ。`_polarity_runs()` が陰イオンの測定を引く）
6. └─ metabolomix/analysis/stats/preprocess.py  preprocess()
7.    └─ metabolomix/analysis/stats/preprocess.py  _select()
8.    └─ metabolomix/analysis/stats/preprocess.py  _polarity_blocks()（統合でなければ None）
9.    ├─ [統合でない] metabolomix/analysis/stats/preprocess.py  _drift_step()
10.   │  └─ metabolomix/analysis/stats/preprocess.py  _drift_eligibility()
11.   │  └─ metabolomix/analysis/stats/preprocess.py  _auto_hold()
12.   │  └─ metabolomix/analysis/stats/preprocess.py  _drift()
13.   │     └─ metabolomix/analysis/stats/preprocess.py  _loess()
14.   │        └─ metabolomix/analysis/stats/preprocess.py  _local_poly()
15.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _merged_run_warnings()（陽陰で注入順・バッチが違えば `polarity_runs_differ`。文面に `run_remedy()` の入れ直す手順）
16.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _drift_blocks()（極性ごとに `_drift_step()`）
17.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _drift_hold_warning_blocks()
18.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _polarity_fingerprint()（`matrix_id` に入れる極性の塊）
19.   └─ metabolomix/analysis/stats/preprocess.py  _filters()
20.   ├─ [統合でない] metabolomix/analysis/stats/preprocess.py  _factor_columns()
21.   ├─ [統合でない] metabolomix/analysis/stats/preprocess.py  _normalize_step()
22.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _normalize_blocks()（極性ごとに `_factor_columns()` → `_normalize_step()`）
23.   │  └─ metabolomix/analysis/stats/preprocess.py  _block_factor_config()（係数の特徴の指定をその極性の特徴に絞る。その極性で係数の特徴が残らなければ、その極性は残った特徴すべて）
24.   └─ metabolomix/analysis/stats/preprocess.py  _collapse()
25. └─ metabolomix/tools/_matrix_store.py  store_matrix()
26.    └─ metabolomix/analysis/stats_results.py  register_matrix()
27.    └─ metabolomix/tools/_matrix_store.py  preprocess_summary()
28.    └─ metabolomix/analysis/stats_records.py  write_stats_record()

行列の登録・要約・前処理の記録は `tools/_matrix_store.py`（ツールではない共有部品。`dataset_compare` / `dataset_pca` /
`pathway_project` の既定の前処理も同じ `register_default` → `store_matrix` を通る）。`injection_columns` は手順 2 が設定に写し、
個体のまとめ（手順 24 の `_collapse`）の鍵から外す。LOESS は `loess_degree=1`（既定）なら statsmodels の lowess、
`2` なら同じ近傍と三乗重みの局所 2 次を numpy で解く（`_local_poly`。近傍の最小数は次数 + 3）。

## dataset_compare

統計エンジンの比較（spec 2026-11-06 §5〜§9・§12）。群の指定（手順 2）を spec v2 の対比へ写す。`reported_feature_ids`
（`bh_population="reported"` の母集団。ほかの母集団なら使わずに注記する）と `output_path` の検証は比較の前に済ませる。入力側は手順 3 の `matrix_inputs`
（`prefer_latest=True`）が決める: `matrix_id` を渡せば登録済みの行列（手順 4）の設定・群・除外のまま、省略すれば
このハンドルで `dataset_preprocess` が最後に作った行列（手順 5。登録の `origin` と順番 `seq` で選ぶ。群の引数を渡したときは
群が同じものだけ＝手順 7 の `groups_match`）、それも無ければ今の除外と群の引数から（手順 6）既定の設定で作る（戻り値の
`matrix_source` が `matrix_id` / `latest_preprocess` / `default`）。`matrix_id` と一緒に渡した群の引数は手順 8 が、
`comparisons` の無い `factors` は手順 12 が、効かない（重複した）引数なら注記にし、読み方で数字が変わるなら `needs_input`
（`patch: true`。選択肢の value は替える引数の dict）を返す（断らない）。群が同じかは上書きの dict でなく、前処理に入る試料
ごとの群で比べる（手順 7・9 の `groups_match`。値の行列を読まない `stats/data.py` の `sample_table` で試料の表だけを組む）。
選択肢はその読みで比べられるものだけ: 行列の群と渡した群のどちらに対比の群の名前があるか（手順 10）、`factors` のモデルが
組めるか（手順 13。値の欠け・空の cell）を先に見る。spec は引数から組み立てて（手順 16）正規化する（手順 17）ので、結果の `spec` を
`dataset_analysis_run` に渡せば同じ数字を出し直せる。行列が無いときは既定の前処理を回し（手順 18。尋ねない: drift=auto が補正を
当てなかったときは警告 `drift_not_applied_by_default` が結果に載る）、行列を登録して前処理の記録も 1 件書く（手順 19）。
`run_spec_v2` が登録済みの行列を再利用して比較し（手順 20。交絡は `compare` が対比ごとにも調べ、既定では止めずに印を付ける）、結果は
`ds.results[result_id]` と `ds.last_differential` に置く（手順 21。記録に `p_adjust_method`）。`output_path` があれば全行の
TSV（と F 検定の TSV）を書き（手順 22）、手順 23・24 の記録は親を行列の前処理の記録にして `reports/stats/` へ追記する
（どちらも書けなくても比較は成功）。戻り値の `top` / `f_top` の名前と `id_level` は手順 22 と同じ `identity_of` で引く。

1. metabolomix/tools/stats_tools.py  dataset_compare()
2. └─ metabolomix/tools/stats_tools.py  _design_and_comparisons()
3. └─ metabolomix/analysis/stats_results.py  matrix_inputs()
4.    └─ metabolomix/analysis/stats_results.py  get_matrix()（`matrix_id` あり）
5.    └─ metabolomix/analysis/stats_results.py  latest_preprocess_matrix()（`matrix_id` なし）
6.    └─ metabolomix/analysis/stats_results.py  study_inputs()（`matrix_id` なし）
7.    └─ metabolomix/analysis/stats_results.py  groups_match()（`matrix_id` なし・群の引数あり）
8. └─ metabolomix/tools/stats_tools.py  _matrix_group_conflict()（`matrix_id` あり）
9.    └─ metabolomix/analysis/stats_results.py  groups_match()
10.   └─ metabolomix/tools/stats_tools.py  _group_levels_under()（群が違うとき）
11. └─ metabolomix/analysis/stats/data.py  from_dataset()
12. └─ metabolomix/tools/stats_tools.py  _factors_without_comparisons()（`comparisons` なし）
13.   └─ metabolomix/tools/stats_tools.py  _factor_design_problem()
14.   └─ metabolomix/tools/stats_tools.py  _factor_model_comparisons()
15. └─ metabolomix/tools/stats_tools.py  _covariate_specs()
16. └─ metabolomix/analysis/stats_results.py  compare_spec()
17. └─ metabolomix/analysis/stats/spec_v2.py  normalize_spec_v2()
18. └─ metabolomix/analysis/stats_results.py  prepare_default()（行列なし）
19. └─ metabolomix/tools/_matrix_store.py  register_default()
20. └─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
21. └─ metabolomix/analysis/stats_results.py  register_compare()
22. └─ metabolomix/analysis/compare_table.py  write_compare_tables()（`output_path` あり）
23. └─ metabolomix/tools/stats_tools.py  _write_compare_record()
24.    └─ metabolomix/analysis/stats_records.py  write_stats_record()

## dataset_pca

統計エンジンの PCA（spec 2026-11-06 §10。stats-pca.v1）。軸は生体の個体で決め、QC は同じ軸へ投影する。
行列は群の引数と無関係に決める（`matrix_inputs(prefer_latest=True)`）: `matrix_id` を渡せば登録済みの行列（手順 2）、
省略すればこのハンドルで `dataset_preprocess` が最後に作った行列、それも無ければ既定の前処理で行列を作って登録する
（手順 3。前処理の記録も 1 件書く）。どれを使ったかは戻り値の `matrix_source`。`group_levels` / `group_factors` は手順 5 が個体ごとの
色分けの群に解くだけで、PCA の数字（手順 4）には効かない。結果は `ds.results[result_id]` と `ds.last_pca`
（`kind="stat_pca"`）に置く（手順 6）。drift=auto が補正を当てなかったときは警告 `drift_not_applied_by_default` が載る。
手順 7 の記録は親を行列の前処理の記録にして `reports/stats/` へ追記する（書けなくても PCA は成功）。
図（`dataset_plot` の `pca` / `loadings`）は `tool_helpers.stat_pca_rows`（スコア）と
`plots/pca_loadings.dataset_features`（`kept_features` のローディング）が `stat_pca` の記録を読む。

1. metabolomix/tools/stats_tools.py  dataset_pca()
2. └─ metabolomix/analysis/stats_results.py  matrix_inputs()
3. └─ metabolomix/tools/_matrix_store.py  register_default()
4. └─ metabolomix/analysis/stats/pca.py  pca()
5. └─ metabolomix/analysis/stats_results.py  unit_group_labels()
6. └─ metabolomix/analysis/stats_results.py  register_pca()
7. └─ metabolomix/tools/stats_tools.py  _write_pca_record()

## dataset_overview

最初の概観（spec 2026-11-06-dataset-overview-design）。`load_dataset(overview=True)` のときだけ、読み込みの最後に同じ
`build_overview_payload()` を呼ぶ（`embedded=True`、`dataset.md` の `load_dataset`。既定の `overview=False` では呼ばない）。統計は C のエンジンだけを通す: 既定の前処理（手順 3〜5。尋ねない:
drift=auto がドリフト補正を当てなかったら補正なしで計算し、そう書く）→ 品質の純計算（手順 7〜10。前処理前の値。陽陰を統合したデータセットでは `failed_qc` と `drift_signs` が極性ごとの総量・その極性の注入順で見る）→ PCA（手順 12〜13）→
spec v2 の比較（手順 15〜18。群が 1 つ・全群をまとめた表で完全交絡なら回さない。一部の対だけの完全交絡では
止めず、`allow_confounded=True` で対比に印を付ける。探索専用のデータセットも比較を回し、`analysis/dataset_export.py` の `exploratory_reason()` の理由で
全体・クラス・分子種の段の先頭に印の事実を置く）→ 選別の後だけ組成（手順 20）。比較と PCA は `make_current=False` で
`ds.results` に置き、直近（`last_differential` / `last_pca`）は動かさない。事実の組み立ては `analysis/overview.py`（手順 6・11・19・21〜24）で、
概観は `ds.results[ovw_…]`（`kind="overview"`、手順 25）と `reports/stats/`（step=overview、親は前処理の記録、手順 26）に残す。

1. metabolomix/tools/overview_tools.py  dataset_overview()
2. └─ metabolomix/tools/overview_tools.py  build_overview_payload()
3.    └─ metabolomix/analysis/stats_results.py  matrix_inputs()
4.    └─ metabolomix/analysis/stats_results.py  prepare_default()
5.    └─ metabolomix/tools/_matrix_store.py  register_default()
6.    └─ metabolomix/analysis/overview.py  group_plan()
7.    └─ metabolomix/analysis/stats/confounding.py  confounding_report()
8.    └─ metabolomix/analysis/stats/quality.py  qc_rsd_summary()
9.    └─ metabolomix/analysis/stats/quality.py  failed_qc()
10.   └─ metabolomix/analysis/stats/quality.py  quality_verdict()
11.   └─ metabolomix/analysis/overview.py  group_comparisons()
12.   └─ metabolomix/analysis/stats/pca.py  pca()
13.   └─ metabolomix/analysis/stats_results.py  register_pca()
14.   └─ metabolomix/tools/overview_tools.py  _compare()
15.      └─ metabolomix/analysis/stats_results.py  compare_spec()
16.      └─ metabolomix/analysis/stats/spec_v2.py  normalize_spec_v2()
17.      └─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
18.      └─ metabolomix/analysis/stats_results.py  register_compare()
19.   └─ metabolomix/analysis/overview.py  chain_state()
20.   └─ metabolomix/analysis/stats/composition.py  composition()
21.   └─ metabolomix/analysis/overview.py  quality_facts()
22.   └─ metabolomix/analysis/overview.py  global_facts()
23.   └─ metabolomix/analysis/overview.py  class_facts()
24.   └─ metabolomix/analysis/overview.py  species_facts()
25.   └─ metabolomix/analysis/overview.py  register_overview()
26.   └─ metabolomix/tools/overview_tools.py  _write_record()

## dataset_export_differential

書き出す前に**その結果が今の前処理から出たものか**を確かめる（手順 3。前処理を
やり直した後の古い結果を書くと、TSV の数字と現在の前処理条件が並んだファイルができ、
どちらも正しく見えてずれが分からない）。**探索専用のデータセット**（中断された実行・検証を通っていない pipeline の単位。
欠けた検体が「その群には無い」ようにしか見えない）は断らず、メタ行 `# exploratory = true` と理由（`exploratory_reason`）を足す。
`dataset_compare` の結果（統計エンジン）は、対比を 1 つ選んで（複数あるときは `contrast_id` か `result_id#contrast_id`。`stats_results.split_result_ref`。
省くと `as_two_group` が `CONTRAST_REQUIRED`〔`ambiguous` と選択肢〕を送出し、ツールが `contrast_envelope` で `needs_input` にする）手順 4 の `as_two_group`（`allow_interaction=True`。
2×2 の交互作用の対比も差の差として読む。`contract_q=True`: `q_value` は契約の BH の q で、比較が `bonferroni` / `holm` / `none` なら
手順 6 の `contract_bh_rows` が同じ母集団で p から BH の q を計算し直し、判定と `significant` もその q で読み直す）で 2 群の形に読み替えてから書き、
メタ行に `stats_meta_lines`（検定・判定規則・BH の母集団・比較の補正の方法 `p_adjust_method`・`q_value_method = bh` と `q_value_source`・
単位・`matrix_id`、交互作用ならその定義など）が足される。`exclude_feature_ids` の行は検定の後に落とし
（q は計算し直さない）、件数をメタ行 `n_excluded_after_test` に書く。前処理条件は登録した行列の設定から取る。出力は一時ファイルへ書いてから
置換するので、途中で落ちた出力が「新しい完成品」として残らない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_differential()
2. ├─ metabolomix/analysis/dataset_export.py  export_dataset_result()
3. │  ├─ metabolomix/analysis/result_state.py  assert_current()
4. │  ├─ metabolomix/analysis/dataset_export.py  exploratory_reason()
5. │  ├─ metabolomix/analysis/stats_results.py  as_two_group()（`dataset_compare` の結果の対比 1 つを 2 群の形に読み替える。各行の `significant`（判定が changed のときだけ）と、`allow_confounded` で完全交絡を越えた比較の `provenance.comparison` をここで付ける）
6. │  │  └─ [比較の補正が BH でない] metabolomix/analysis/stats_results.py  contract_bh_rows()（比較の q が付いた行〔補正の母集団〕の p から BH の q、判定はその q で `compare._verdict` を読み直す）
7. │  ├─ metabolomix/analysis/export_contract.py  row_significant()（行の印を読む。印の無い行は `is_significant()` かつ判定が changed）
8. │  ├─ metabolomix/analysis/dataset_export.py  _meta_lines()
9. │  │  ├─ metabolomix/analysis/dataset_export.py  _engine_preprocess()
10.│  │  ├─ metabolomix/analysis/export_contract.py  stats_meta_lines()
11.│  │  │  └─ metabolomix/analysis/export_contract.py  q_value_line()
12.│  │  └─ metabolomix/analysis/export_contract.py  build_meta()
13.│  ├─ metabolomix/analysis/export_contract.py  format_row()
14.│  ├─ metabolomix/analysis/dataset_export.py  _atomic_write_text()
15.│  └─ metabolomix/analysis/dataset_export.py  recomputed_q_note()（BH の q を計算し直したときの戻り値の `notes`）
16.└─ [CONTRAST_REQUIRED ほか] metabolomix/tools/_contrast_choice.py  contrast_envelope()

## dataset_feature_search

読み取り専用。`targets` の有無で窓検索と目標照合に分かれ、どちらも最後に同じ統計の
結合を通る。判定の純ロジックは `metabolomix/analysis/feature_search.py` にあり、session も
MCP も import しない。ツール層の仕事は、極性の出所を決めること（ジョブの宣言 → mzTab-M の
付加体の多数決）、差次的結果を選ぶこと（名指しの ID が無い・古ければ `missing_state`、
省略時は今も有効な 2 群比較だけ）、丸めと TSV 化、注記。キュレーションのフラグは
`dataset_export_differential` と同じ見つけ方で数えるだけで、**反映しない**（読めなくても
検索は止めず、注記にする）。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_feature_search()
2. ├─ [窓検索] metabolomix/analysis/feature_search.py  search_window()
3. │  └─ metabolomix/analysis/feature_search.py  _window_identity()
4. ├─ [目標照合] metabolomix/mztab/polarity.py  dataset_polarity()
5. │  ├─ metabolomix/console/job_manager.py  load_job()
6. │  └─ metabolomix/console/output_collector.py  read_adduct_polarity()
7. ├─ [目標照合] metabolomix/analysis/feature_search.py  normalize_targets()
8. │  └─ metabolomix/analysis/feature_search.py  theoretical_mz()
9. ├─ [目標照合] metabolomix/analysis/feature_search.py  match_target()
10. │  └─ metabolomix/analysis/feature_search.py  _identity()（SME の付加体は綴りでなく組成で比べる: `msdial/adducts.py` `same_adduct()`。種の目標の SME・SML の名前は MS/MS の合否と同じ `analysis/msms_rules.py` `name_lipid_species()` で読む）
11. ├─ [目標照合] metabolomix/analysis/feature_search.py  estimate_rt_offset()
12. ├─ [目標照合] metabolomix/analysis/feature_search.py  rule_hash()
13. ├─ [include_statistics=true] metabolomix/analysis/feature_search.py  select_results()
14. │  └─ metabolomix/analysis/result_state.py  is_current()
15. ├─ [merge_split=false] metabolomix/analysis/feature_search.py  feature_stats()
16. │  └─ metabolomix/analysis/feature_search.py  _detected_counter()
17. ├─ [merge_split=true] metabolomix/analysis/feature_search.py  merge_plan()
18. │  ├─ [パスなし] metabolomix/analysis/alignment_qc.py  remembered()
19. │  ├─ [最後の結果なし] metabolomix/analysis/alignment_qc.py  find_split_pairs()
20. │  ├─ [パスあり] metabolomix/analysis/alignment_qc.py  read_qc_file()
21. │  └─ metabolomix/analysis/alignment_qc.py  usable_merge_groups()
22. ├─ [merge_split=true] metabolomix/analysis/feature_search.py  merged_feature_stats()
23. │  ├─ metabolomix/analysis/feature_search.py  feature_stats()
24. │  └─ metabolomix/analysis/feature_search.py  merged_results()
25. │     ├─ metabolomix/analysis/feature_search.py  merged_dataset()
26. │     │  └─ metabolomix/analysis/alignment_qc.py  merge_matrix_columns()
27. │     └─ metabolomix/analysis/feature_search.py  retest_merged()
28. │        ├─ metabolomix/analysis/stats/data.py  from_dataset()
29. │        ├─ metabolomix/analysis/stats/data.py  merge_feature_columns()
30. │        ├─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
31. │        └─ metabolomix/analysis/stats_results.py  as_two_group()
32. ├─ [DomainError（select_results、統計を付けるときの merge_plan・merged_feature_stats が送出）] metabolomix/tools/_contrast_choice.py  contrast_envelope()
33. ├─ metabolomix/analysis/feature_search.py  result_summaries()
34. ├─ metabolomix/analysis/result_state.py  dataset_fingerprint()
35. └─ metabolomix/tools/dataset_analysis_tools.py  _curation_note()
36.    ├─ metabolomix/curation/apply.py  arf2_for_mztab()
37.    └─ metabolomix/curation/apply.py  flags_for_arf2()

`include_statistics=false` では `select_results` と統計の枝（`feature_stats` 〜 `result_summaries`）を通らない。目標照合の `merge_split=true` は候補を数える前（`match_target` の前）に `merge_plan` を呼ぶので、その組だけを戻り値の `merge_split` に写す（窓検索では統計を付けるときだけ呼ぶ）。`estimate_rt_offset` は `rt_min` を持つ目標と互いに一意に合う特徴の組（目標の付加体ごとに一致がたかだか 1 つで、それらが同時溶出すれば 1 組）から RT のずれを測るだけで、判定には当てない（当てるのは `rt_offset_min` を渡したときの `match_target`）。`feature_stats` の行の `detected_units_a` / `detected_units_b` / `presence_p` / `detection_driven` は `as_two_group` が結果の行から写したもの。窓検索に渡した `ppm_offset` / `rt_offset_min` は効かないので、`notes` に書くだけで照合には通さない（形は検める）。目標照合の `rt_offset_min` は `normalize_targets` が `rt_min` のある目標にだけ足すので、どの目標にも足されなければツールが `notes` に書き、`rule.rt_offset_min` を出さない（`rule_hash` は正規化後の目標から作るので、もともと入らない）。目標照合に渡した `sort="intensity"` と、`merge_split=false` と一緒の `alignment_qc_path` も効かないので、誤りにせず `notes` に書くだけ（`alignment_qc_path` は読まない）。`select_results` と、統計を付けるときの `merge_plan`（窓検索）・`merged_feature_stats` が送出した `DomainError` は、ツール自身がそれぞれの呼び出しの直後に `contrast_envelope` で封筒にする（`select_results` の中から呼ぶのではない。対比を省いた曖昧さなら `needs_input`〔`argument="contrast_id"`〕、未知の対比や `FEATURE_MERGE_UNAVAILABLE` などほかの誤りは `mztab_error`）。目標照合が照合の前に呼ぶ `merge_plan` の誤りは、目標の誤りと同じく直に `mztab_error` にする。

`dataset_compare` の結果（`compare` の記録）は `select_results` が対比ごとに 2 群の形に広げる（対比が複数なら ID は `result_id#contrast_id`、交互作用は除く）。`merge_split` の検定のやり直しは、その結果の spec で `stats/data.py` `from_dataset()` → `merge_feature_columns()` → `spec_v2.run_spec_v2()` → `as_two_group()` を通す（`retest_merged`）。

`merge_split`（既定 true）は、割れの組（`alignment_qc` の `merge_groups`）がデータセットに 1 つでもあれば、
組を試料ごとの大きい方の値で 1 列にまとめた行列（前処理の**前**の生の値・検出状態は OR）に、差次的結果を
計算した spec を当て直し、同じ検定と**まとめた後の全特徴**での BH をやり直して、`stats` の
**すべての行**をそこから出す（1 つの表に 2 つの BH の母集団を混ぜない）。組が無ければ検定はやり直さない。
目標照合では `merge_plan` を照合の前に作り、`match_target(merged_from=merge_map(...))` が組を 1 つの候補
（残す ID）として数える。

`rule_hash` は正規化後の目標（既定値を埋め、付加体の層を解決したもの）から作り、統計の出所は
混ぜない。覆うのは照合を左右する欄（`target_id`・InChIKey 先頭ブロック・分子式か expected_mz・
付加体の層と理論 m/z・ppm・RT・目標の極性・`applicable`）だけで、表示名 `name` は入れず、目標は
`target_id` の順に並べ替えてから作る（改名・`targets` の並べ替えでは変わらない）。hash は表示用の丸め（`_rule_targets()` が理論 m/z を 4 桁に丸める）より前の値から
作るので、返された `rule.targets` から計算し直すことはできない——照合は `rule_hash` の値どうしで行う。
result_id は計算のたびに発行される UUID なので入れず、統計の出所は `stats_source`
（`result_ids` と `dataset_fingerprint` の `dataset_hash`）に分けて返す。

## dataset_set_sample_metadata

実験情報シート(sample-manifest.v1 / v2)を読み、全件検証してから一括反映する。版は
先頭の schema 行で決まり、列見出しからは推測しない。検証で1件でも落ちれば
`apply_metadata` の契約により `session.dataset` は一切変更しない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_set_sample_metadata()
2. └─ metabolomix/analysis/sample_manifest.py  apply_sample_manifest()
3.    └─ metabolomix/analysis/sample_manifest.py  _raw_manifest_layout()
4.    └─ metabolomix/analysis/sample_manifest.py  parse_manifest()（任意の列 `subject_id` を読む。`sample_meta` に持ち、統計エンジンの前処理が反復注入を 1 個体にまとめる鍵に使う。schema の行の直後の宣言 `# injection_order_source = ...` は `_read_declarations()` が読み、注入順の出所に `derived_from` を足す）
5.    └─ metabolomix/analysis/sample_manifest.py  resolve_metadata()
6.    └─ metabolomix/analysis/sample_manifest.py  apply_metadata()
7.       └─ metabolomix/analysis/result_state.py  invalidate_results()
8. └─ metabolomix/analysis/stats/data.py  merged_feature_polarity()（陽陰を統合したデータセットなら戻り値に `notes`）
9.    └─ metabolomix/analysis/stats/preprocess.py  run_remedy()（陰イオンの測定を入れ直す手順の文面）

陽陰を統合したデータセットでは、シートは陽イオンの raw を挙げる（統合したデータセットの raw のパスは陽イオンのもの）
ので、`batch`・`injection_order` は陽イオンの測定だけを決める。陰イオンの測定は統合が残した
`quant_state.polarity_runs.negative` のままで、変えるなら陰イオンのハンドルにシートを当ててから統合し直す
（手順 8・9 が戻り値の `notes` にそう書く）。

対話ツール `dataset_compare`（統計エンジン）は、群名だけで対照/処置の向きを決めない・QC/blank/unknown/standard・`include=false` を混ぜない・完全交絡は既定では比べて印を付け（`confounding` と対比の `confounded`）、`allow_confounded=false` を渡したときだけ止める（`STATISTIC_CONFOUNDED`）。

## dataset_export_evidence_bundle

全特徴量と候補の根拠と、試料ごとの値と検出を新しいv3バンドルへ保存する。旧15列TSVは変更しない。
`result_id` を省けば直近の差次的結果（`ds.last_differential`）。対比は書き出しと同じ選び方（省いて 2 群の対比が複数なら
`needs_input`）。`features.tsv` の `q_value` と `significant` は書き出しと同じく契約の BH の q（`as_two_group(contract_q=True)`）で、
`manifest.json` の `comparison.correction` はいつも `BH`。比較の補正が BH でなければ `comparison.p_adjust_method` と `q_value_source` を足す。
探索専用のデータセットも断らず、manifest の最上位と `comparison` に `exploratory`（印と理由）を足す。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_evidence_bundle()
2.    └─ metabolomix/tools/handle_tools.py  _alignment_files()（対話のセッションのとき。ハンドルはその `.arf2` / `.dcl`）
3.    └─ metabolomix/analysis/evidence_export.py  export_dataset_evidence()
4.       └─ metabolomix/analysis/dataset_export.py  exploratory_reason()
5.       └─ metabolomix/analysis/stats_results.py  as_two_group()（`contract_q=True`。比較の補正が BH でなければ `contract_bh_rows()`）
6.       └─ metabolomix/analysis/evidence_export.py  comparison_of()
7.       └─ metabolomix/analysis/evidence_export.py  _samples_table()
8.       └─ metabolomix/analysis/evidence_export.py  _engine_intensity_rows()
9.       └─ metabolomix/analysis/evidence_export.py  library_candidate()（SME の候補ごと。計算した理論 m/z は写さない）
10.      └─ metabolomix/analysis/evidence_export.py  dataset_spectrum_links()
11.         └─ metabolomix/analysis/evidence_export.py  _id_space_check()
12.         └─ metabolomix/analysis/evidence_export.py  link_alignment_spectra()
13.      └─ metabolomix/analysis/evidence_export.py  write_evidence_bundle()
14.      └─ metabolomix/analysis/dataset_export.py  recomputed_q_note()（BH の q を計算し直したときの戻り値の `notes`）
15.   └─ [CONTRAST_REQUIRED ほか] metabolomix/tools/_contrast_choice.py  contrast_envelope()

spectrum_links は同じアラインメントの `.dcl` に ID（SMF_ID == AlignmentID == `.dcl` の番号）で結び、
m/z・RT は検証にだけ使う（`link_alignment_spectra()` が結ぶ）。
結べないときの状態と理由は manifest の provenance（`docs/output_format/reanalysis.md` §10）。
`annotations.tsv` の `msdial_sme` 行の来歴は SME の候補そのものだが、組成式と付加体から計算した理論 m/z
（`DatasetState.theoretical_mz_computed` の特徴。`directory=` で開いたデータセット）は `library_candidate()` が空にして写す
（`theoretical_mass_to_charge` はライブラリの記録のプリカーサー m/z の意味のまま。データセットの候補は変えない）。
