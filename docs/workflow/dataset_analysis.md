# DatasetState 解析ツール呼び出し連鎖

解析フォルダ（`.arf` / `.arf2`）も mzTab-M も、`load_dataset` が同じ `DatasetState` にして
このツール群へ渡す。旧 ARF 経路（`arf_*` ツールとそのセッションスロット）は B-2b Task 8 で消えた。

## dataset_preprocess

統計エンジンの前処理（spec 2026-11-06 §4）。行列 ID を返し、`ds.stat_matrices[matrix_id]` に置く
（v1 の `pp_*` には書かない）。除外（`ds.excluded_*`）と群の上書き（`group_levels` / `group_factors`）は
手順 3 が入力へ写す。`needs_input`（内部標準で換算済み＋QC のドリフト補正）は手順 5 の戻りから封筒にする。
手順 7 の記録は `reports/stats/` への追記で、書き先が決まらない・書けないときも前処理は成功する。

1. metabolomix/tools/stats_tools.py  dataset_preprocess()
2. └─ metabolomix/analysis/stats/preprocess.py  normalize_preprocess_config()
3. └─ metabolomix/analysis/stats_results.py  study_inputs()
4. └─ metabolomix/analysis/stats/data.py  from_dataset()
5. └─ metabolomix/analysis/stats/preprocess.py  preprocess()
6. └─ metabolomix/analysis/stats_results.py  register_matrix()
7. └─ metabolomix/analysis/stats_records.py  write_stats_record()

（v1 の前処理 `dataset_analysis.run_dataset_preprocess` は、pipeline の前処理工程と旧
`dataset_differential` / `dataset_pca` がまだ読む `pp_*` を作る。D 段で消える。）

## dataset_compare

統計エンジンの比較（spec 2026-11-06 §5〜§9・§12）。群の指定（手順 2）を spec v2 の対比へ写し、`matrix_id` を
渡せば登録済みの行列（手順 4）の設定・群・除外のまま、省略すれば既定の前処理で行列を作って登録し
（手順 5・6。前処理の記録も 1 件書く）、`run_spec_v2` で比較する（手順 7）。結果は `ds.results[result_id]` と
`ds.last_differential` に置く（手順 8）。`needs_input`（換算済み＋QC）は手順 5 の戻りから封筒にする。
手順 9 の記録は親を行列の前処理の記録にして `reports/stats/` へ追記する（書けなくても比較は成功）。

1. metabolomix/tools/stats_tools.py  dataset_compare()
2. └─ metabolomix/tools/stats_tools.py  _design_and_comparisons()
3. └─ metabolomix/tools/stats_tools.py  _covariate_specs()
4. └─ metabolomix/analysis/stats_results.py  get_matrix()
5. └─ metabolomix/analysis/stats/spec_v2.py  prepare_spec_v2()
6. └─ metabolomix/tools/_matrix_store.py  store_matrix()
7. └─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
8. └─ metabolomix/analysis/stats_results.py  register_compare()
9. └─ metabolomix/tools/stats_tools.py  _write_compare_record()

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
3. └─ metabolomix/tools/_matrix_store.py  store_matrix()
4. └─ metabolomix/analysis/stats/pca.py  pca()
5. └─ metabolomix/analysis/stats_results.py  unit_group_labels()
6. └─ metabolomix/analysis/stats_results.py  register_pca()
7. └─ metabolomix/tools/stats_tools.py  _write_pca_record()

## dataset_differential

1. metabolomix/tools/dataset_analysis_tools.py  dataset_differential()
2. └─ metabolomix/analysis/dataset_analysis.py  run_dataset_differential()
3.    └─ metabolomix/analysis/differential.py  check_confounding()
4.    └─ metabolomix/analysis/differential.py  two_group_test()
5.    └─ metabolomix/analysis/differential.py  add_fdr()
6.    └─ metabolomix/analysis/differential.py  summarize_two_group()
7.    └─ metabolomix/analysis/differential.py  volcano_data()

## dataset_export_differential

書き出す前に 2 つ確かめる。**その結果が今の前処理から出たものか**（手順 3。前処理を
やり直した後の古い結果を書くと、TSV の数字と現在の前処理条件が並んだファイルができ、
どちらも正しく見えてずれが分からない）と、**探索専用のデータセットでないか**
（中断された実行の出力では、欠けた検体が「その群には無い」ようにしか見えない）。
メタ行の前処理条件は `ds.preprocessing_recipe`（今の状態）ではなく、その結果が親に
持つ前処理の `provenance.effective_parameters` から取る。pipeline run から引き継いだ結果
（来歴に `handoff`）は親の前処理が `ds.results` に無いので、引き継ぎ時に写した
`provenance.handoff.preprocess` を使い、無ければ `unknown（…）` と書く——session の
`preprocessing_recipe` へは決してフォールバックしない。出力は一時ファイルへ書いてから
置換するので、途中で落ちた出力が「新しい完成品」として残らない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_differential()
2. └─ metabolomix/analysis/dataset_export.py  export_dataset_result()
3. │  └─ metabolomix/analysis/result_state.py  assert_current()
4. │  └─ metabolomix/analysis/export_contract.py  is_significant()
5. │  └─ metabolomix/analysis/dataset_export.py  _meta_lines()
6. │  │  └─ metabolomix/analysis/dataset_export.py  _preprocess_parameters()
7. │  │  └─ metabolomix/analysis/export_contract.py  build_meta()
8. │  └─ metabolomix/analysis/export_contract.py  format_row()
9. │  └─ metabolomix/analysis/dataset_export.py  _atomic_write_text()

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
26. │        ├─ metabolomix/analysis/feature_search.py  _result_recipe()
27. │        ├─ metabolomix/analysis/dataset_analysis.py  run_dataset_preprocess()
28. │        ├─ metabolomix/analysis/differential.py  two_group_test()
29. │        └─ metabolomix/analysis/differential.py  add_fdr()
30. ├─ metabolomix/analysis/feature_search.py  result_summaries()
31. ├─ metabolomix/analysis/result_state.py  dataset_fingerprint()
32. └─ metabolomix/tools/dataset_analysis_tools.py  _curation_note()
33.    ├─ metabolomix/curation/apply.py  arf2_for_mztab()
34.    └─ metabolomix/curation/apply.py  flags_for_arf2()

`merge_split`（既定 true）は、割れの組（`alignment_qc` の `merge_groups`）がデータセットに 1 つでもあれば、
組を試料ごとの大きい方の値で 1 列にまとめた行列（前処理の**前**の生の値・検出状態は OR）に、差次的結果を
計算したレシピ（pipeline の run から引き継いだ結果は来歴の `handoff.preprocess`、session の結果は親の前処理の
`effective_parameters`）を当て直し、同じ 2 群検定と**まとめた後の全特徴**での BH をやり直して、`stats` の
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
2. └─ metabolomix/analysis/dataset_service.py  apply_sample_manifest()
3.    └─ metabolomix/analysis/dataset_service.py  _raw_manifest_layout()
4.    └─ metabolomix/analysis/sample_manifest.py  parse_manifest()（任意の列 `subject_id` を読む。`sample_meta` に持つだけで統計には使わない）
5.    └─ metabolomix/analysis/sample_manifest.py  resolve_metadata()
6.    └─ metabolomix/analysis/sample_manifest.py  apply_metadata()
7.       └─ metabolomix/analysis/result_state.py  invalidate_results()

pipeline が比較を明示するときに使う `resolve_comparison`/`run_comparison`
（`metabolomix/analysis/dataset_service.py`）は、群名だけで対照/処置の向きを決めない・
QC/blank/unknown/standard・`include=false` を混ぜない（選択規則は
`metabolomix/analysis/sample_manifest.py` `select_statistical_samples()` が唯一の定義）・
完全交絡を `allow_confounded=true` の
明示なしには通さない、という前提検証を `compare_dataset()` の前に挟む。単体ツール
`dataset_differential` はこの前提検証を経ない汎用呼び出し（USAGE.md 参照）。

## dataset_export_evidence_bundle

全特徴量と候補の根拠と、試料ごとの値と検出を新しいv3バンドルへ保存する。旧15列TSVは変更しない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_evidence_bundle()
2.    └─ metabolomix/tools/handle_tools.py  _alignment_files()（対話のセッションのとき。ハンドルはその `.arf2` / `.dcl`）
3.    └─ metabolomix/analysis/evidence_export.py  export_dataset_evidence()
4.       └─ metabolomix/analysis/evidence_export.py  _samples_table()
5.       └─ metabolomix/analysis/evidence_export.py  _intensity_rows()
6.       └─ metabolomix/analysis/evidence_export.py  dataset_spectrum_links()
7.          └─ metabolomix/analysis/evidence_export.py  _id_space_check()
8.          └─ metabolomix/analysis/evidence_export.py  link_alignment_spectra()
9.       └─ metabolomix/analysis/evidence_export.py  write_evidence_bundle()

spectrum_links は同じアラインメントの `.dcl` に ID（SMF_ID == AlignmentID == `.dcl` の番号）で結び、
m/z・RT は検証にだけ使う（`link_alignment_spectra()` が結ぶ）。
結べないときの状態と理由は manifest の provenance（`docs/output_format/reanalysis.md` §10）。
