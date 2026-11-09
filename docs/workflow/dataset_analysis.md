# DatasetState 解析ツール呼び出し連鎖

解析フォルダ（`.arf` / `.arf2`）も mzTab-M も、`load_dataset` が同じ `DatasetState` にして
このツール群へ渡す。旧 ARF 経路（`arf_*` ツールとそのセッションスロット）は B-2b Task 8 で消えた。

## dataset_preprocess

統計エンジンの前処理（spec 2026-11-06 §4）。行列 ID を返し、`ds.stat_matrices[matrix_id]` に置く
（v1 の `pp_*` には書かない）。除外（`ds.excluded_*`）と群の上書き（`group_levels` / `group_factors`）は
手順 3 が入力へ写す。`needs_input`（内部標準で換算済み＋QC のドリフト補正）は手順 5 の戻りから封筒にする。
手順 9 の記録は `reports/stats/` への追記で、書き先が決まらない・書けないときも前処理は成功する。

1. metabolomix/tools/stats_tools.py  dataset_preprocess()
2. └─ metabolomix/analysis/stats/preprocess.py  normalize_preprocess_config()
3. └─ metabolomix/analysis/stats_results.py  study_inputs()
4. └─ metabolomix/analysis/stats/data.py  from_dataset()
5. └─ metabolomix/analysis/stats/preprocess.py  preprocess()
6. └─ metabolomix/tools/_matrix_store.py  store_matrix()
7.    └─ metabolomix/analysis/stats_results.py  register_matrix()
8.    └─ metabolomix/tools/_matrix_store.py  preprocess_summary()
9.    └─ metabolomix/analysis/stats_records.py  write_stats_record()

行列の登録・要約・前処理の記録は `tools/_matrix_store.py`（ツールではない共有部品。`dataset_compare` / `dataset_pca` /
`pathway_project` の既定の前処理も同じ `register_default` → `store_matrix` を通る）。`injection_columns` は手順 2 が設定に写し、
個体のまとめ（手順 5 の `preprocess()` の中の `collapse`）の鍵から外す。

## dataset_compare

統計エンジンの比較（spec 2026-11-06 §5〜§9・§12）。群の指定（手順 2）を spec v2 の対比へ写す。入力側は手順 3 の
`matrix_inputs` が決める: `matrix_id` を渡せば登録済みの行列（手順 4）の設定・群・除外のまま、省略すれば今の除外と
群の引数から（手順 5）既定の設定で作る。spec は引数から組み立てて（手順 8）正規化する（手順 9）ので、結果の `spec` を
`dataset_analysis_run` に渡せば同じ数字を出し直せる。`matrix_id` が無いときは既定の前処理を回し（手順 10・11。`needs_input`
＝換算済み＋QC はここの戻りから封筒にする）、行列を登録して前処理の記録も 1 件書く（手順 12）。`run_spec_v2` が登録済みの
行列を再利用して比較し（手順 13）、結果は `ds.results[result_id]` と `ds.last_differential` に置く（手順 14）。
手順 15・16 の記録は親を行列の前処理の記録にして `reports/stats/` へ追記する（書けなくても比較は成功）。

1. metabolomix/tools/stats_tools.py  dataset_compare()
2. └─ metabolomix/tools/stats_tools.py  _design_and_comparisons()
3. └─ metabolomix/analysis/stats_results.py  matrix_inputs()
4.    └─ metabolomix/analysis/stats_results.py  get_matrix()（`matrix_id` あり）
5.    └─ metabolomix/analysis/stats_results.py  study_inputs()（`matrix_id` なし）
6. └─ metabolomix/analysis/stats/data.py  from_dataset()
7. └─ metabolomix/tools/stats_tools.py  _covariate_specs()
8. └─ metabolomix/analysis/stats_results.py  compare_spec()
9. └─ metabolomix/analysis/stats/spec_v2.py  normalize_spec_v2()
10. └─ metabolomix/tools/stats_tools.py  _prepare_default()（`matrix_id` なし）
11.    └─ metabolomix/analysis/stats_results.py  prepare_default()
12. └─ metabolomix/tools/_matrix_store.py  register_default()
13. └─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
14. └─ metabolomix/analysis/stats_results.py  register_compare()
15. └─ metabolomix/tools/stats_tools.py  _write_compare_record()
16.    └─ metabolomix/analysis/stats_records.py  write_stats_record()

## dataset_pca

統計エンジンの PCA（spec 2026-11-06 §10。stats-pca.v1）。軸は生体の個体で決め、QC は同じ軸へ投影する。
行列は群の引数と無関係に決める: `matrix_id` を渡せば登録済みの行列（手順 2）、省略すれば既定の前処理で
行列を作って登録する（手順 3。前処理の記録も 1 件書く）。`group_levels` / `group_factors` は手順 5 が個体ごとの
色分けの群に解くだけで、PCA の数字（手順 4）には効かない。結果は `ds.results[result_id]` と `ds.last_pca`
（`kind="stat_pca"`）に置く（手順 6）。`needs_input`（換算済み＋QC）は手順 3 の戻りから封筒にする。
手順 7 の記録は親を行列の前処理の記録にして `reports/stats/` へ追記する（書けなくても PCA は成功）。
図（`dataset_plot` の `pca` / `loadings`）は `tool_helpers.stat_pca_rows`（スコア）と
`plots/pca_loadings.dataset_features`（`kept_features` のローディング）が `stat_pca` の記録を読む。

1. metabolomix/tools/stats_tools.py  dataset_pca()
2. └─ metabolomix/analysis/stats_results.py  get_matrix()
3. └─ metabolomix/tools/_matrix_store.py  register_default()
4. └─ metabolomix/analysis/stats/pca.py  pca()
5. └─ metabolomix/analysis/stats_results.py  unit_group_labels()
6. └─ metabolomix/analysis/stats_results.py  register_pca()
7. └─ metabolomix/tools/stats_tools.py  _write_pca_record()

## dataset_overview

最初の概観（spec 2026-11-06-dataset-overview-design）。`load_dataset` も読み込みの最後に同じ `build_overview_payload()` を呼ぶ
（`embedded=True`、`dataset.md` の `load_dataset`）。統計は C のエンジンだけを通す: 既定の前処理（手順 3〜5。ドリフト補正の既定が
未決〔`NeedsInput`〕なら drift を off にして作り直し、そう書く）→ 品質の純計算（手順 7〜10。前処理前の値）→ PCA（手順 12〜13）→
spec v2 の比較（手順 15〜18。群が 1 つ・探索専用・完全交絡なら回さない）→ 選別の後だけ組成（手順 20）。比較と PCA は `make_current=False` で
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

書き出す前に 2 つ確かめる。**その結果が今の前処理から出たものか**（手順 3。前処理を
やり直した後の古い結果を書くと、TSV の数字と現在の前処理条件が並んだファイルができ、
どちらも正しく見えてずれが分からない）と、**探索専用のデータセットでないか**
（中断された実行の出力では、欠けた検体が「その群には無い」ようにしか見えない）。
`dataset_compare` の結果（統計エンジン）は、対比を 1 つ選んで（複数あるときは `contrast_id` か `result_id#contrast_id`。`stats_results.split_result_ref`）手順 4 の `as_two_group` で 2 群の形に読み替えてから書き、メタ行に `stats_meta_lines`（検定・判定規則・BH の母集団・単位・`matrix_id` など）が足される。前処理条件は登録した行列の設定から取る。出力は一時ファイルへ書いてから
置換するので、途中で落ちた出力が「新しい完成品」として残らない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_differential()
2. └─ metabolomix/analysis/dataset_export.py  export_dataset_result()
3. │  └─ metabolomix/analysis/result_state.py  assert_current()
4. │  └─ metabolomix/analysis/stats_results.py  as_two_group()（`dataset_compare` の結果の対比 1 つを 2 群の形に読み替える。各行の `significant`（判定が changed のときだけ）と、`allow_confounded` で完全交絡を越えた比較の `provenance.comparison` をここで付ける）
5. │  └─ metabolomix/analysis/export_contract.py  row_significant()（行の印を読む。印の無い行は `is_significant()` かつ判定が changed）
6. │  └─ metabolomix/analysis/dataset_export.py  _meta_lines()
7. │  │  └─ metabolomix/analysis/dataset_export.py  _engine_preprocess()
8. │  │  └─ metabolomix/analysis/export_contract.py  stats_meta_lines()
9. │  │  └─ metabolomix/analysis/export_contract.py  build_meta()
10. │  └─ metabolomix/analysis/export_contract.py  format_row()
11. │  └─ metabolomix/analysis/dataset_export.py  _atomic_write_text()

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
10. │  └─ metabolomix/analysis/feature_search.py  _identity()
11. ├─ [目標照合] metabolomix/analysis/feature_search.py  rule_hash()
12. ├─ metabolomix/analysis/feature_search.py  select_results()
13. │  └─ metabolomix/analysis/result_state.py  is_current()
14. ├─ [merge_split=false] metabolomix/analysis/feature_search.py  feature_stats()
15. │  └─ metabolomix/analysis/feature_search.py  _detected_counter()
16. ├─ [merge_split=true] metabolomix/analysis/feature_search.py  merge_plan()
17. │  ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
18. │  ├─ metabolomix/analysis/alignment_qc.py  read_qc_file()
19. │  └─ metabolomix/analysis/alignment_qc.py  usable_merge_groups()
20. ├─ [merge_split=true] metabolomix/analysis/feature_search.py  merged_feature_stats()
21. │  ├─ metabolomix/analysis/feature_search.py  feature_stats()
22. │  └─ metabolomix/analysis/feature_search.py  merged_results()
23. │     ├─ metabolomix/analysis/feature_search.py  merged_dataset()
24. │     │  └─ metabolomix/analysis/alignment_qc.py  merge_matrix_columns()
25. │     └─ metabolomix/analysis/feature_search.py  retest_merged()
26. │        ├─ metabolomix/analysis/stats/data.py  from_dataset()
27. │        ├─ metabolomix/analysis/stats/data.py  merge_feature_columns()
28. │        ├─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
29. │        └─ metabolomix/analysis/stats_results.py  as_two_group()
30. ├─ metabolomix/analysis/feature_search.py  result_summaries()
31. ├─ metabolomix/analysis/result_state.py  dataset_fingerprint()
32. └─ metabolomix/tools/dataset_analysis_tools.py  _curation_note()
33.    ├─ metabolomix/curation/apply.py  arf2_for_mztab()
34.    └─ metabolomix/curation/apply.py  flags_for_arf2()

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
4.    └─ metabolomix/analysis/sample_manifest.py  parse_manifest()（任意の列 `subject_id` を読む。`sample_meta` に持ち、統計エンジンの前処理が反復注入を 1 個体にまとめる鍵に使う）
5.    └─ metabolomix/analysis/sample_manifest.py  resolve_metadata()
6.    └─ metabolomix/analysis/sample_manifest.py  apply_metadata()
7.       └─ metabolomix/analysis/result_state.py  invalidate_results()

対話ツール `dataset_compare`（統計エンジン）は、群名だけで対照/処置の向きを決めない・QC/blank/unknown/standard・`include=false` を混ぜない・完全交絡を `allow_confounded=true` の明示なしには通さない（`STATISTIC_CONFOUNDED`）。

## dataset_export_evidence_bundle

全特徴量と候補の根拠と、試料ごとの値と検出を新しいv3バンドルへ保存する。旧15列TSVは変更しない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_evidence_bundle()
2.    └─ metabolomix/tools/handle_tools.py  _alignment_files()（対話のセッションのとき。ハンドルはその `.arf2` / `.dcl`）
3.    └─ metabolomix/analysis/evidence_export.py  export_dataset_evidence()
4.       └─ metabolomix/analysis/evidence_export.py  _samples_table()
5.       └─ metabolomix/analysis/evidence_export.py  _engine_intensity_rows()
6.       └─ metabolomix/analysis/evidence_export.py  dataset_spectrum_links()
7.          └─ metabolomix/analysis/evidence_export.py  _id_space_check()
8.          └─ metabolomix/analysis/evidence_export.py  link_alignment_spectra()
9.       └─ metabolomix/analysis/evidence_export.py  write_evidence_bundle()

spectrum_links は同じアラインメントの `.dcl` に ID（SMF_ID == AlignmentID == `.dcl` の番号）で結び、
m/z・RT は検証にだけ使う（`link_alignment_spectra()` が結ぶ）。
結べないときの状態と理由は manifest の provenance（`docs/output_format/reanalysis.md` §10）。
