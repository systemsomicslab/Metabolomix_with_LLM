# ワークフロー: 再解析のデータセットのハンドルと spec のツール

**再解析で使う**経路（spec 2026-10-22）。並行の再解析のエージェントが、対話と同じ MCP の
ツールを**データセットのハンドル**と宣言的な spec で呼ぶ。対話の使い方
（`load_dataset` → `dataset_*`）は変わらない。値の**意味**は
[reanalysis.md（output_format）](../output_format/reanalysis.md) を見ること。

ハンドルの状態は `metabolomix/core/handles.py` の登録簿（`handles.registry`）に住む。
`load_dataset` が開いたデータセットは登録簿のハンドルで、いつも直近になり、`session_state.session` は
その直近のハンドル（`current_handle_id`）を指すだけ（`session.dataset` / `session.library` は
property）。`DatasetHandle` は `dataset` と `library` を `AnalysisSession` と同じ名前で持つので、
既存のツールは `holder_for(dataset_handle)` が返したものを読むだけでハンドルにも対話の
セッションにも効く。

`dataset_handle` を渡した呼び出しは `metabolomix/core/handle_dispatch.py` の
`install()` が包んだ関数がワーカースレッドで回し、そのハンドルの `lock` で直列化する。
渡さない呼び出しも、直近のハンドルがあればそのハンドルの `lock` の下でワーカースレッドで回し、
`session_state.pinned_handle()` でそのハンドルに固定する。ロックは最長 `HANDLE_LOCK_WAIT`（25 秒）だけ待ち
（`_acquire()`）、空かなければツールを呼ばずに `HANDLE_BUSY`（`mcp_errors.handle_busy()`。ロックを持つ長い道具の
仕事があれば、その受付番号と道具。`DatasetHandle.running_job`）を返す。ロックを取った時点でそのハンドルが
閉じられていれば、ツールを呼ばずに `missing_state` を返す。直近のハンドルが無ければ
イベントループのスレッドで同期に呼ぶ。`dataset_handle` を持たないのに `session.dataset` /
`.library` を読むツール（`handle_dispatch.SESSION_READERS`）も同じ経路で包む。

## ハンドルを開く（load_dataset）

`load_dataset` は [dataset.md](dataset.md) の節が連鎖を持つ（`dataset_list` / `dataset_close` も同じ文書）。再解析のエージェントは、返った
`dataset_handle` を各ツールに渡して使い、直近（`current_handle_id`）を当てにしない。
開いたハンドルはいつも直近になるが、`dataset_handle` を渡した呼び出しは直近を読み書きしない（例外は定量の派生で、
`metabolomix/quant/derive.py` の `register()` が `adopt_handle()` で派生のハンドルを直近にする——`dataset_handle` を渡しても。
[quant.md](quant.md)）。
ライブラリはデータセットが覚えた同定のライブラリ（`library_record`）を、sha256 を測り直して
読み取り専用で開く（極性の設定には落ちない）。スポットの照合に使う `.arf2` / `.dcl` は、
解析フォルダで読んだならその選択（`AlignmentSelection`）のもの、mzTab-M で読んだなら
`_alignment_files()` が探す。

## dataset_analysis_run

前提: `input.source=dataset` ならデータセット（ハンドルか対話のセッション。無ければ
`missing_state`）。外部の表（`input.source` が `tsv` / `maf` / `mwtab` / `xlsx`）ならデータセットは要らない。
状態変更: なし（読むだけ。`dataset_alignment_qc` がハンドルに持つ最後の結果も読まない——spec だけで再現する）。
`output_path` を渡すと結果の JSON を書く。

1. metabolomix/tools/reanalysis_tools.py  dataset_analysis_run()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/analysis/analysis_spec.py  normalize_spec()
4. │  ├─ metabolomix/analysis/analysis_spec.py  _normalize_external()
5. │  ├─ metabolomix/analysis/analysis_spec.py  _normalize_features()
6. │  ├─ metabolomix/analysis/analysis_spec.py  _normalize_correlate()
7. │  └─ metabolomix/analysis/recipe_matrix.py  normalize_recipe()
8. ├─ metabolomix/analysis/analysis_spec.py  run_spec()
9. │  ├─ metabolomix/analysis/analysis_spec.py  _load()
10. │  │  └─ metabolomix/analysis/analysis_spec.py  _external_matrix()
11. │  │     └─ [maf / mwtab / xlsx] metabolomix/analysis/external_tables.py  read_table()
12. │  ├─ metabolomix/analysis/analysis_spec.py  _merge_plan()
13. │  │  ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
14. │  │  ├─ metabolomix/analysis/alignment_qc.py  read_qc_file()
15. │  │  ├─ metabolomix/analysis/alignment_qc.py  usable_merge_groups()
16. │  │  └─ metabolomix/analysis/alignment_qc.py  merge_matrix_columns()
17. │  ├─ [統合したデータセット・drift_correct] metabolomix/analysis/analysis_spec.py  _polarity_orders()（極性ごとの注入順。TSV の注入順が陰イオンの測定と違えば `notes`）
18. │  ├─ metabolomix/analysis/analysis_spec.py  _run_analysis()
19. │  │  ├─ metabolomix/analysis/analysis_spec.py  _merged_polarity()（入力が陽陰を統合したデータセットなら 特徴 → 極性）
20. │  │  ├─ metabolomix/analysis/recipe_matrix.py  build_matrix()
21. │  │  │  ├─ metabolomix/analysis/recipe_matrix.py  knn_impute()
22. │  │  │  ├─ [統合したデータセットの normalize] metabolomix/analysis/recipe_matrix.py  _normalize_by_polarity()
23. │  │  │  └─ [統合したデータセットの drift_correct] metabolomix/analysis/recipe_matrix.py  _drift_correct_by_polarity()
24. │  │  ├─ metabolomix/analysis/analysis_spec.py  _aggregate()
25. │  │  ├─ [和] metabolomix/analysis/analysis_spec.py  _contributions()
26. │  │  ├─ metabolomix/analysis/analysis_spec.py  _design()
27. │  │  ├─ metabolomix/analysis/spec_statistics.py  two_group()
28. │  │  ├─ metabolomix/analysis/spec_statistics.py  ols_group()
29. │  │  ├─ metabolomix/analysis/spec_statistics.py  fold_change()
30. │  │  ├─ metabolomix/analysis/analysis_spec.py  _fill_blocks()（読んだ値と、補完の段の直後の値 `imputed_values`）
31. │  │  ├─ metabolomix/analysis/analysis_spec.py  _fill_constant()
32. │  │  │  └─ metabolomix/analysis/spec_statistics.py  constant_groups()
33. │  │  └─ metabolomix/analysis/spec_statistics.py  bh()
34. │  ├─ [correlate] metabolomix/analysis/analysis_spec.py  _correlate()
35. │  │  ├─ metabolomix/analysis/analysis_spec.py  _external_matrix()
36. │  │  ├─ metabolomix/analysis/analysis_spec.py  _sample_map()
37. │  │  └─ metabolomix/analysis/analysis_spec.py  _correlate_rows()
38. │  │     └─ metabolomix/analysis/spec_statistics.py  spearman()
39. │  └─ metabolomix/analysis/analysis_spec.py  rule_hash()
40. ├─ metabolomix/tools/reanalysis_tools.py  write_result()
41. │  └─ metabolomix/tools/reanalysis_tools.py  result_bytes()
42. ├─ [output_path なし] metabolomix/analysis/stats_run_summary.py  first_by_p()
43. ├─ metabolomix/tools/reanalysis_tools.py  _aggregate_summary()
44. ├─ [集約の要約を切った] metabolomix/tools/reanalysis_tools.py  _aggregate_cut_notes()
45. └─ [correlate] metabolomix/tools/reanalysis_tools.py  _correlate_summary()

書いても効かない設定（`steps` に無い段の設定・`merge_split=false` や外部の表での `features.alignment_qc`・
ols でない検定の共変量）は `normalize_spec(spec, notes=...)` が捨てて理由を `notes` に足し、戻り値と JSON の
`notes` に出す（捨てた設定は `rule_hash` に入らない）。`correction=none` の `bh_population`・`thresholds.q` と、
`_uses_injection_order()` が偽（`drift_correct` の解析が無い）の `injection_order`・`_uses_covariates()` が偽（共変量で調整する
ols が無い）の `covariates` も `notes` に書くが、spec には残す（`rule_hash` は前と同じ）。`_load()` はこの 2 つの表を読まない。
入力が陽陰を統合したデータセットなら、`drift_correct` は
`_polarity_orders()` が決めた極性ごとの注入順で当てはめる（`metadata` は陽イオンが試料の表、陰イオンが
`quant_state.polarity_runs.negative`。TSV は 1 つの注入順を全極性に当て、陰イオンの測定と違えば
`run_spec(..., notes=...)` が `notes` に足す）。`output_path` が無く 50 行を超える解析は
`first_by_p()` の先頭 50 行だけを TSV にし、`truncated`・`rows_order` を立てる（断らない）。
集約の要約は `_aggregate_summary()` がリストごとに自分の長さで切り（寄与は和に入ったメンバーが 50 を超えた
ときだけ `share_median` の大きい順に並べ替えて切る）、`_aggregate_cut_notes()` が切ったものだけを `notes` に
名指す。

`spec["schema"] == "analysis-spec.v2"` のときは、`resolve_holder()` のあと v1 の鎖でなく次へ振り分ける
（データセットは常に要る。除外は spec の `input` だけを使い、`dataset_exclude` の除外は読まない。
状態変更・記録は無し）。

1. metabolomix/tools/reanalysis_tools.py  _analysis_run_v2()
2. ├─ metabolomix/analysis/stats/spec_v2.py  normalize_spec_v2()
3. ├─ metabolomix/analysis/stats/spec_v2.py  run_spec_v2()
4. ├─ metabolomix/analysis/stats_run_summary.py  run_summary()
5. │  └─ metabolomix/analysis/stats_run_summary.py  analysis_summary()
6. │     ├─ [compare / composition] metabolomix/analysis/stats_run_summary.py  confounding_summary()
7. │     └─ [composition] metabolomix/analysis/stats_run_summary.py  composition_warnings()
8. ├─ metabolomix/analysis/stats_run_summary.py  inline_rows()
9. └─ metabolomix/tools/reanalysis_tools.py  write_result()

compare / composition の要約は交絡 `confounding`（`batch_status`・`by_contrast`。composition は内側の比較のもの）を持ち、
composition の `warnings` は `composition_warnings()` が内側の比較（絶対量・CLR）の注意をコードで重ねずに足す。

割れたスポットのまとめ（`features.merge_split`、既定 true）は recipe の前（行列を読んだ直後。recipe・集約・検定はまとめた行列に当たる）に、
`alignment_qc` の `merge_groups` だけに行う（`features.alignment_qc` に JSON を渡せばその群、
無ければその場の `find_split_pairs`。JSON は `dataset_fingerprint` が今のデータセットと合うものだけ。
`dataset_alignment_qc` がハンドルに持つ最後の結果は使わない）。外部の表（`input.source` が `dataset` 以外）では何もしない。

## library_msms_qualify

前提: ライブラリ（`library_load` か、ライブラリつきの `load_dataset`）。無ければ
`missing_state`。状態変更: なし（直近の照合 `last_match` も書かない）。

`ion_mode` は手順 4 で正規化し（読めない値は `DATASET_BAD_REQUEST`）、`mz_tol_ppm` は 0 より大きい有限の数だけを
受ける。前駆体の窓は `library_match_feature` と同じ（m/z 500 を超えたら同じ ppm。手順 9 の中）。候補の表の行は
`msms_rules.candidate_rows()`（`dataset_claim_evidence` の `msms.top` と共有）。

1. metabolomix/tools/reanalysis_tools.py  library_msms_qualify()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ [ライブラリ無しで記憶が lbm_only のとき] metabolomix/library/tools.py  lbm_only_not_evaluable()（`not_evaluable` を返す）
4. ├─ metabolomix/library/matching.py  normalize_ion_mode()
5. ├─ metabolomix/analysis/msms_rules.py  normalize_rules()（種の目標に同位体の標識があれば `isotope_labelled()` で見て `MSMS_RULES_INVALID`）
6. │  └─ metabolomix/analysis/msms_rules.py  normalize_ions()
7. ├─ metabolomix/library/tools.py  _spot_query()
8. ├─ metabolomix/library/tools.py  _measured_spectrum()
9. ├─ metabolomix/library/matching.py  score_candidates()
10. ├─ metabolomix/analysis/msms_rules.py  query_window()
11. ├─ metabolomix/analysis/msms_rules.py  qualify()
12. │  ├─ metabolomix/analysis/msms_rules.py  record_lipid_species()  [種の目標のとき]
13. │  ├─ metabolomix/analysis/spectral_match.py  fragment_matches_excluding_precursor()
14. │  ├─ metabolomix/analysis/msms_rules.py  _lipid_class()
15. │  │  └─ metabolomix/msdial/lipid_identity.py  normalize_lipid_name()
16. │  └─ metabolomix/analysis/msms_rules.py  find_ions()
17. ├─ metabolomix/analysis/msms_rules.py  candidate_rows()
18. │  ├─ metabolomix/library/matching.py  delta_ppm()
19. │  ├─ metabolomix/library/matching.py  passes_dbs_cutoffs()  [`.dbs` に足切りがあるとき]
20. │  └─ metabolomix/analysis/msms_rules.py  record_lipid_species()  [種の目標のとき]
21. ├─ metabolomix/library/matching.py  msdial_agreement()  [spot_id のスポットに注釈があるとき]
22. └─ metabolomix/analysis/msms_rules.py  rules_hash()
23. └─ metabolomix/library/tools.py  identified_by_lbm()
24.    └─ metabolomix/dataset/library_memory.py  is_lbm_annotator()  [記憶が lbm_also で spot_id の照合のとき]

## dataset_homolog_check

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_homolog_check()
2. ├─ metabolomix/analysis/homolog.py  normalize_spec()
3. ├─ metabolomix/analysis/homolog.py  members_from_names()（`auto.adduct` は `msdial/adducts.py` `same_adduct()` で組成で比べる）
4. │  └─ metabolomix/analysis/homolog.py  species_of()
5. │     └─ metabolomix/analysis/msms_rules.py  name_lipid_species()（目標照合・MS/MS の合否と同じ読み手。標識付きは読まない）
6. ├─ metabolomix/analysis/homolog.py  check()（比べた隣どうしの対が 0 なら `not_evaluable`。ツールが `notes` に理由を書く）
7. └─ metabolomix/analysis/homolog.py  rule_hash()

## dataset_feature_relations

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_feature_relations()
2. ├─ metabolomix/analysis/claim_evidence.py  feature_relations()
3. │  ├─ metabolomix/analysis/claim_evidence.py  formula_index()
4. │  ├─ metabolomix/curation/relations.py  profile_correlation()
5. │  └─ metabolomix/curation/relations.py  find_relations()（Y 自身の付加体は `msdial/adducts.py` `same_adduct()` で組成で除き、同位体の期待比は Y のイオンの組成の `curation/isotope.py` `isotopologues()` のうち X の窓に入るものの和）
6. └─ metabolomix/tools/reanalysis_tools.py  rows_table()

## dataset_claim_evidence

前提: データセット（無ければ `missing_state`）。`msms` を頼むときはライブラリ（無ければ候補ごとに
`no_library`。同定が LBM だけのデータセット＝記憶が `lbm_only` なら候補ごとに `not_evaluable`
（`reason: lbm_only`）で、理由の文は束の `msms_not_evaluable_reason` に 1 回）。状態変更: なし。`output_path` を渡すと JSON を書く
（省けば束を返し、候補が 10 を超えれば先頭 10 件と `truncated: true`。断らない）。

1. metabolomix/tools/reanalysis_tools.py  dataset_claim_evidence()
2. ├─ metabolomix/analysis/claim_evidence.py  normalize_claim_spec()
3. │  ├─ metabolomix/analysis/feature_search.py  normalize_targets()（種の目標に同位体の標識があれば `FEATURE_TARGET_INVALID`）
4. │  └─ metabolomix/analysis/msms_rules.py  normalize_rules()
5. ├─ metabolomix/analysis/claim_evidence.py  build()
6. │  ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
7. │  ├─ metabolomix/analysis/alignment_qc.py  split_partners()
8. │  ├─ metabolomix/analysis/feature_search.py  match_target()
9. │  ├─ metabolomix/analysis/claim_evidence.py  _msms_for()
10. │  │  ├─ metabolomix/library/tools.py  _spot_query()
11. │  │  ├─ metabolomix/library/matching.py  score_candidates()（`msms.mz_tol_ppm` を渡す）
12. │  │  ├─ metabolomix/analysis/msms_rules.py  candidate_rows()（`library_msms_qualify` と同じ行）
13. │  │  ├─ metabolomix/library/matching.py  msdial_agreement()
14. │  │  └─ metabolomix/analysis/msms_rules.py  qualify()
15. │  ├─ metabolomix/analysis/claim_evidence.py  feature_relations()
16. │  ├─ metabolomix/analysis/claim_evidence.py  _selection()
17. │  └─ metabolomix/analysis/claim_evidence.py  rule_hash()
18. ├─ [output_path あり] metabolomix/tools/reanalysis_tools.py  write_result()
19. └─ [output_path なし・候補 > 10] metabolomix/tools/reanalysis_tools.py  _inline_candidates()

spec の `rt_offset_min`（任意）は `normalize_claim_spec` が形を検め（誤りは `CLAIM_SPEC_INVALID`）、目標に `rt_min` が
あれば既定の目標と走査の全変種の `normalize_targets(rt_offset_min=...)` に渡す（`match_target` が RT の窓の中心を
ずらす）。`rt_min` が無ければ効かないので捨てて `notes` に書く（束の `notes` と戻り値。規則にも `rule_hash` にも入らない）。
`output_path` を省いて候補が 10 を超えるときは断らず、`_inline_candidates` が既定の特徴 → 既定の照合で見えた候補 →
ほかの候補の順（同じ段の中は束の並び）に先頭 10 件を選び、`truncated: true` と全件の数を返す（束そのものは変えない）。

候補ごとの `split` と、既定の特徴・走査の判定でまとめる割れの組は `find_split_pairs`（既定の閾値）の `split_partners` /
`merge_groups`（`dataset_alignment_qc` がハンドルに持つ最後の結果は読まない。束は spec と入力だけで再現する）。`default_selection.counts.qualified`
（と走査ごとの `counts.qualified`）は、規則で qualified になった特徴のうち MS/MS の合否
（`msms.qualify.verdict`）が `qualified` の数（`_selection`。MS/MS を頼んでいない・ライブラリが無ければ null）で、
規則で qualified の数（旧 `qualified` の意味。MS/MS の有無ではない）は `rule_qualified`（`msms_available` は同じ値の非推奨の別名）。
代表 MS/MS の有無で数えるのは `msms_spectrum_available`（質量の窓の候補。割れの組は 1 つ、RT で絞らない）で、候補ごとの
`msms.spectrum_available` は `_msms_for` がライブラリの有無より先に spot を照会して決める（MS/MS を頼まなければ null。`.arf2` / `.dcl` のファイルが無い `not_found`（`reason` が `no_arf2` / `no_dcl`）は False でなく null）。`msms_spectrum_available` は `_spectrum_count` が数え、0 は照会して無かった・null は頼まなかったか分からなかった（True の単位が無く null の単位がある）・数は下限。既定の特徴と走査の判定は
`match_target(merged_from=merge_map(...))` で割れの組を 1 つの候補（残す ID）として数える。

## dataset_mass_accuracy

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_mass_accuracy()
2. ├─ [subset=msms_matched] metabolomix/analysis/mass_accuracy.py  has_match_information()
3. ├─ metabolomix/analysis/mass_accuracy.py  sme_ppm()
4. ├─ metabolomix/analysis/mass_accuracy.py  band_summary()
5. │  └─ metabolomix/analysis/mass_accuracy.py  summarize()
6. └─ metabolomix/analysis/mass_accuracy.py  injection_trend()
7.    └─ metabolomix/analysis/alignment_qc.py  sample_injection_order()

注入順の出所が `filename`（採取時刻が読めずファイル名の順に落ちた）のときは、`injection_trend` の
`spearman_rho` を出さない（`null`。理由は `rho_unavailable_reason: "filename_order"`）。`per_injection` と
`n_samples` は残す。

## dataset_sample_peaks

前提: データセット（無ければ `missing_state`）。状態変更: なし（`output_path` を渡せばそのファイルを書き、既存は置き換える。
注釈は読み取り専用ではない）。解析フォルダ（`.arf`）から読んだ
データセットだけ（`per_sample_rt` が無ければ `DATASET_BAD_REQUEST`）。`output_path` が無く行（特徴 × 試料）が
1,000 を超えれば、先頭の特徴の行を特徴ごとにそろえて載せ `truncated` を立てる（断らない）。

1. metabolomix/tools/reanalysis_tools.py  dataset_sample_peaks()

## dataset_novel_candidates

前提: データセット（無ければ `missing_state`）。`statistics.result_id` なら今も有効な 2 群比較
（無ければ `missing_state` が `dataset_compare` を指す）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_novel_candidates()
2. ├─ metabolomix/tools/reanalysis_tools.py  _drop_ineffective_novel_settings()
3. │  └─ [analysis_spec] metabolomix/analysis/analysis_spec.py  normalize_spec()
4. ├─ metabolomix/analysis/mass_accuracy.py  normalize_novel_spec()
5. │  └─ [statistics.contrast_id] metabolomix/analysis/stats_results.py  split_result_ref()
6. ├─ metabolomix/analysis/feature_search.py  select_results()
7. ├─ [result_id] metabolomix/tools/reanalysis_tools.py  _novel_result_stats()
8. │  ├─ metabolomix/analysis/feature_search.py  merge_plan()
9. │  │  ├─ [パスなし] metabolomix/analysis/alignment_qc.py  remembered()
10. │  │  ├─ [最後の結果なし] metabolomix/analysis/alignment_qc.py  find_split_pairs()
11. │  │  └─ [パスあり] metabolomix/analysis/alignment_qc.py  read_qc_file()
12. │  └─ metabolomix/analysis/feature_search.py  merged_results()
13. │     ├─ metabolomix/analysis/feature_search.py  merged_dataset()
14. │     └─ metabolomix/analysis/feature_search.py  retest_merged()
15. ├─ [analysis_spec] metabolomix/analysis/analysis_spec.py  run_spec()
16. ├─ metabolomix/analysis/mass_accuracy.py  novel_candidates()
17. │  └─ metabolomix/analysis/mass_accuracy.py  sme_ppm()
18. ├─ metabolomix/analysis/mass_accuracy.py  rule_hash()
19. └─ metabolomix/tools/reanalysis_tools.py  write_result()

`statistics.contrast_id`（任意）は `normalize_novel_spec()` が `split_result_ref()` で `result_id` の `#` の後ろと突き合わせ
（食い違えば `DATASET_BAD_REQUEST`）、`result_id#contrast_id` の 1 つの形にしてから規則に載せる（書かない spec の正規化と
`rule_hash` は前と同じ）。対比が複数の結果を ID だけで渡せば `select_results()` の `CONTRAST_REQUIRED` を
`_contrast_choice.contrast_envelope(spec_field="statistics.result_id")` が `needs_input` にする。

`statistics.alignment_qc` を省いた `result_id` の経路は、`dataset_alignment_qc` がこのハンドルに持つ最後の結果
（`alignment_qc.remembered()`。データの指紋と `metadata_revision` が作ったときと同じものだけ）を使い、無ければ
その場で探す（`merge_split.source`）。効かない設定（`merge_split=false` の `statistics.alignment_qc`、`analysis_spec` の隣の `statistics.contrast_id`、
埋め込んだ analysis_spec の効かない設定）は `_drop_ineffective_novel_settings()` が外して `notes` に書く。
analysis_spec で選んだ解析が `correction=none` なら、q が無いので候補が出ないことも `notes` に書く。
`output_path` が無く 50 行を超えれば q の昇順の先頭 50 行だけを TSV にし、`truncated` を立てる（断らない）。

割れたスポットの組は 1 候補（残す ID＝受け持ちの多い側（同数なら小さい ID））に数える。`statistics.result_id` では `statistics.merge_split`
（既定 true）で、組が 1 つでもあればすべての特徴の統計をまとめた行列で検定し直した値にする（`dataset_feature_search` の
`merge_split` と同じ実装）。analysis_spec では spec の `features.merge_split` が決め、`run_spec` の
`merge_split.groups` を使う。まとめた組の hash は `rule_hash` に入る（`merge_split` が `groups_hash` に
出す。`merge_split=false` は v7 と同じ hash）。

`msp_then_lbm` の run をパスウェイ（`pathway_project`）へ渡すと、受け渡しが MSP の行も LBM の行と同じ規則で入れる
（Ontology が MS-DIAL の脂質クラスならそのまま、そうでなければ名前から、名前で読めなければ ClassyFire のクラスと組成式から
MS-DIAL の名前とクラスに直す）。直せない行（クラスが 1 つに決まらない・対応表に無い・名前の頭部の語がクラスと食い違う・
組成式が合わない）と、名前に同位体の標識がある MSP の行（内部標準）は推測せず、
理由ごとの件数（`handoff` の `kept_by_route` と `dropped_*`）つきで落とす。使った行と落とした行の数を `handoff` で確かめる
（試しでの件数は [pathway.md（output_format）](../output_format/pathway.md) の規則と
`docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md`）。

## dataset_alignment_qc

前提: データセット（無ければ `missing_state`）。状態変更: 結果をハンドル（データセット）の
`last_alignment_qc` に最後の 1 つとして持つ（`alignment_qc.remember()`。前の結果は置き換える。解析の結果・行列・
記録は書かない）。`output_path` を渡すと同じ中身（同じバイト列）の JSON を書く。
注入順は `ds.sample_metadata_rows`（pipeline の run では実験情報シート）の `injection_order`、
1 試料でも欠ければファイル名の順。受け持ちは `ds.detected_mask`（`.arf` の検出状態）があれば検出、
無ければ値の大きい方。

1. metabolomix/tools/reanalysis_tools.py  dataset_alignment_qc()
2. ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
3. │  ├─ metabolomix/analysis/alignment_qc.py  _injection_order()
4. │  ├─ metabolomix/analysis/alignment_qc.py  _owners()
5. │  ├─ metabolomix/analysis/alignment_qc.py  _pair()
6. │  │  ├─ metabolomix/analysis/alignment_qc.py  _block()
7. │  │  ├─ metabolomix/analysis/alignment_qc.py  _group_aligned()
8. │  │  └─ metabolomix/analysis/alignment_qc.py  names_compatible()
9. │  │     └─ metabolomix/msdial/lipid_identity.py  lipid_species_of()
10. │  ├─ metabolomix/analysis/alignment_qc.py  _mutual_nearest()
11. │  ├─ metabolomix/analysis/alignment_qc.py  _components()
12. │  │  └─ metabolomix/analysis/alignment_qc.py  id_order()
13. │  └─ metabolomix/analysis/result_state.py  dataset_fingerprint()
14. ├─ metabolomix/tools/reanalysis_tools.py  result_bytes()
15. ├─ metabolomix/analysis/alignment_qc.py  remember()
16. └─ metabolomix/tools/reanalysis_tools.py  rows_table()

ハンドルに持った結果は、`dataset_feature_search` の `alignment_qc_path`・`dataset_novel_candidates` の
`statistics.alignment_qc` を省いたときに `feature_search.merge_plan()` → `alignment_qc.remembered()` が読む
（`merge_split.source: last_alignment_qc`）。データの指紋か `metadata_revision`（`dataset_exclude`・
`dataset_set_sample_metadata` が進める）が変わっていれば、`remembered()` が捨てて None を返し、その場で探す。
`dataset_analysis_run` と `dataset_claim_evidence`（`claim_evidence.build()` はいつも既定の閾値の `find_split_pairs()`）は読まない。

まとめは別のツールの側で行う: `alignment_qc.merge_matrix_columns()`（試料ごとに大きい方の値）と
`alignment_qc.split_partners()`（特徴 → 割れの相手）が唯一の実装。
