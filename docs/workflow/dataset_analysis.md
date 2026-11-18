# DatasetState 解析ツール呼び出し連鎖

解析フォルダ（`.arf` / `.arf2`）も mzTab-M も、`load_dataset` が同じ `DatasetState` にしてこのツール群へ渡す。

## dataset_preprocess

統計エンジンの前処理（引数と戻り値の意味は `docs/output_format/dataset.md` の 10.2）。行列 ID を返し、
`ds.stat_matrices[matrix_id]` に置く。
- 手順 2: 効かない任意のキー（`normalize="none"` の係数の特徴・空の `normalization_features`・`drift_mode="off"` の
  `loess_degree`）を誤りにせず設定から捨てる（`_drop_inert`。`matrix_id` は渡さないときと同じ。注記は手順 6 が生の設定を
  正規化し直すときに `warnings` に載せる）。`injection_columns` を設定に写し、個体のまとめ（手順 24 の `_collapse`）の鍵から外す。
- 手順 3: 除外（`ds.excluded_*`）と群の上書き（`group_levels` / `group_factors`）を入力へ写す。
- 手順 6 は問いを返さない: `drift_mode="auto"` は、補正できても内部標準で換算済み・注入順がファイルの順かもしれないときは当てずに計算し
  （`_auto_hold`）、ドリフトの段に理由と `rerun_with`、警告に `drift_not_applied_by_default` を載せる。正規化の係数は
  `normalization_features` の部分集合からも求める（`_factor_columns`。除いた特徴も行列に残す）。LOESS は `loess_degree=1`
  （既定）なら statsmodels の lowess、`2` なら `_local_poly`（近傍の最小数は次数 + 3）。
- 手順 28 の記録は `reports/stats/` への追記で、書き先が決まらない・書けないときも前処理は成功する。

陽陰を統合したデータセット（`dataset_merge_polarities`。すべての特徴が `origin.polarity` を持つ）: 手順 5 が特徴 → 極性と、統合が
`quant_state.polarity_runs` に残した陰イオンの測定（バッチ・注入順）を `StudyData` に載せ、手順 8 の `_polarity_blocks` が極性の
塊（列と、その極性の測定としての試料の表 `polarity_samples`）に分ける。ドリフト補正（手順 16）と正規化（手順 22）は塊ごとに同じ
`_drift_step` / `_normalize_step` を当て、段に `by_polarity` を残す。足切り（手順 19）は塊に分けない。統合したデータセットの試料の
表は陽イオンのものなので、陽イオンの塊は陽イオンだけのデータセットと同じ試料の表を読む。

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
15.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _merged_run_warnings()（陽陰で注入順・バッチが違えば `polarity_runs_differ`。文面に `run_remedy()` の入れ直す手順。陽イオンの注意を最後に並べる）
16.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _drift_blocks()（極性ごとに `_drift_step()`）
17.   ├─ [統合したデータセット] metabolomix/analysis/stats/preprocess.py  _drift_hold_warning_blocks()（注入順が確かでない極性の入れ直し方は `runs_remedy()`。上から行えば全部が入る順＝陽イオンを最後）
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
`pathway_project` の既定の前処理も同じ `register_default` → `store_matrix` を通る）。

## dataset_compare

統計エンジンの比較（意味は `docs/output_format/dataset.md` の 11.2・11.3）。群名だけで対照/処置の向きを決めず、
QC/blank/unknown/standard・`include=false` を比較群に混ぜない。
- 手順 2: 群の指定を spec v2 の対比へ写す。`reported_feature_ids` と `output_path` の検証は比較の前に済ませる。
- 手順 3: 入力側は `matrix_inputs(prefer_latest=True)` が決める: `matrix_id` を渡せば登録済みの行列（手順 4）の設定・群・除外
  のまま、省けばこのハンドルで `dataset_preprocess` が最後に作った行列（手順 5。登録の `origin` と順番 `seq` で選ぶ。群の引数を
  渡したときは群が同じものだけ＝手順 7）、それも無ければ今の除外と群の引数から既定の設定で（手順 6）。
- `matrix_id` と一緒に渡した群の引数は手順 8 が、`comparisons` の無い `factors` は手順 12 が見て、効かない（重複した）引数なら
  注記、読み方で数字が変わるなら `needs_input`（`patch: true`）を返す（断らない）。群が同じかは前処理に入る試料ごとの群で比べる
  （手順 7・9 の `groups_match`。`stats/data.py` の `sample_table` で試料の表だけを組み、値の行列を読まない）。選択肢はその読みで
  比べられるものだけ: 行列の群と渡した群のどちらに対比の群の名前があるか（手順 10）、`factors` のモデルが組めるか（手順 13。
  値の欠け・空の cell）を先に見る。
- 手順 16・17: spec を引数から組み立てて正規化する（結果の `spec` を `dataset_analysis_run` に渡せば同じ数字を出し直せる）。
- 手順 18・19: 行列が無いときは既定の前処理を回し（尋ねない。drift=auto が補正を当てなかったら警告
  `drift_not_applied_by_default`）、行列を登録して前処理の記録も 1 件書く。
- 手順 20: `run_spec_v2` が登録済みの行列で比較する（交絡は `compare` が対比ごとにも調べ、既定では止めずに印を付ける）。
- 手順 21: 結果を `ds.results[result_id]` と `ds.last_differential` に置く（記録に `p_adjust_method`）。
- 手順 22: `output_path` があれば全行の TSV（と F 検定の TSV）を書く。手順 23・24 の記録は親を行列の前処理の記録にして
  `reports/stats/` へ追記する（どちらも書けなくても比較は成功）。戻り値の `top` / `f_top` の名前と `id_level` は手順 22 と同じ
  `identity_of` で引く。

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

統計エンジンの PCA（stats-pca.v1。意味は `docs/output_format/dataset.md` の 10.4）。行列は群の引数と無関係に
`matrix_inputs(prefer_latest=True)`（手順 2）が決め、無ければ既定の前処理で作って登録する（手順 3。前処理の記録も 1 件書く）。
`group_levels` / `group_factors` は手順 5 が個体ごとの色分けの群に解くだけで、手順 4 の数字には効かない。結果は
`ds.results[result_id]` と `ds.last_pca`（`kind="stat_pca"`）に置く（手順 6）。手順 7 の記録は親を行列の前処理の記録にして
`reports/stats/` へ追記する（書けなくても PCA は成功）。図（`dataset_plot` の `pca` / `loadings`）は `tool_helpers.stat_pca_rows`
（スコア）と `plots/pca_loadings.dataset_features`（`kept_features` のローディング）が `stat_pca` の記録を読む。

1. metabolomix/tools/stats_tools.py  dataset_pca()
2. └─ metabolomix/analysis/stats_results.py  matrix_inputs()
3. └─ metabolomix/tools/_matrix_store.py  register_default()
4. └─ metabolomix/analysis/stats/pca.py  pca()
5. └─ metabolomix/analysis/stats_results.py  unit_group_labels()
6. └─ metabolomix/analysis/stats_results.py  register_pca()
7. └─ metabolomix/tools/stats_tools.py  _write_pca_record()

## dataset_overview

最初の概観（意味は `docs/output_format/dataset.md` の 11.9）。`load_dataset(overview=True)` のときだけ、読み込みの最後に同じ
`build_overview_payload()` を呼ぶ（`embedded=True`。既定の `overview=False` では呼ばない）。統計は C のエンジンだけを、次の順に通す:
- 既定の前処理（手順 3〜5。尋ねない。drift=auto が補正を当てなかったら補正なしで計算し、そう書く。その注記は
  `tools/overview_tools.py` の `_drift_note()` が書き、注入順の出所〔MS-DIAL の既定・シートの宣言〕は前処理に渡した試料の表から読む）
- 品質の純計算（手順 7〜10。前処理前の値。陽陰を統合したデータセットでは `failed_qc` と `drift_signs` が極性ごとの総量・その極性の
  注入順で見る）。ρ に使った QC の注入順がファイルの順かもしれない行には `analysis/stats/data.py` の `order_may_be_file_order()` で
  印を付け、手順 21 の `quality_facts()` がドリフトの兆候の事実の但し書きにする。統合したデータセットではその極性の入れ直し方を
  `analysis/stats/preprocess.py` の `runs_remedy()`（極性ごとの `run_remedy()` を陽イオンを最後に並べる。`dataset_preprocess` の
  注意と同じ文と順）から書く（`next` の手の順は 11.9）。subject_id が無いときの技術反復らしいまとまりも前処理前の値から
  `analysis/stats/quality.py` の `replicate_blocks()` が見つけ（候補を挙げるだけ）、手順 10 の判定と手順 21 の事実に渡す
- PCA（手順 12〜13）
- spec v2 の比較（手順 15〜18。群が 1 つ・全群をまとめた表で完全交絡なら回さない。一部の対だけの完全交絡では止めず、
  `allow_confounded=True` で対比に印を付ける。探索専用のデータセットは `analysis/dataset_export.py` の `exploratory_reason()` の
  理由で段の先頭に印の事実を置く）
- 選別の後だけ組成（手順 20。脂質アッセイでは、注釈が `msdial/lipid_identity.py` の `lipid_annotation_class()`（MS-DIAL の
  脂質クラスか、名前を `lipid_name_parses()` が GOSLIN で読めるか、研究室 MSP の ClassyFire の脂質クラスと組成式から
  `pathway/msp_lipid.py` の `from_class_formula()` が作った略記を読めるか）で脂質と言えない特徴をクラスから外して `non_lipid` に
  数え、MS-DIAL の脂質クラスが決まった特徴はそのクラスへ移して `lipid_class_relabel` に数える。手順 23 の `class_facts()` が
  どちらも段の最後の事実にする）

比較と PCA は `make_current=False` で `ds.results` に置き、直近（`last_differential` / `last_pca`）は動かさない。事実の組み立ては
`analysis/overview.py`（手順 6・11・19・21〜24）。概観は `ds.results[ovw_…]`（`kind="overview"`、手順 25）と `reports/stats/`
（step=overview、親は前処理の記録、手順 26）に残す。`top_n` を上限に下げたときの注記は `tools/overview_tools.py` の `_cap_note()`
（多群なら概観の群〔手順 6〕と行列〔手順 3〜5 の `matrix_id`〕を名指しした `dataset_compare(groups=..., matrix_id=...,
output_path=...)` を挙げる）。

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

TSV の契約は `docs/output_format/dataset.md` の 11.4。
- 手順 3: その結果が今の前処理から出たものかを確かめる。手順 4: 探索専用のデータセット（中断された実行・検証を通っていない
  pipeline の単位）は断らず、メタ行 `# exploratory = true` と理由（`exploratory_reason`）を足す。
- `dataset_compare` の結果は対比を 1 つ選ぶ（複数あるときは `contrast_id` か `result_id#contrast_id`。`stats_results.split_result_ref`。
  省くと `as_two_group` が `CONTRAST_REQUIRED`〔`ambiguous` と選択肢〕を送出し、ツールが `contrast_envelope` で `needs_input` に
  する）。手順 5 の `as_two_group` で 2 群の形に読み替えてから書く（`allow_interaction=True`: 2×2 の交互作用の対比も差の差として
  読む。`contract_q=True`: `q_value` は契約の BH の q で、比較が `bonferroni` / `holm` / `none` なら手順 6 の `contract_bh_rows` が
  同じ母集団で p から BH の q を計算し直し、判定と `significant` もその q で読み直す）。
- メタ行に `stats_meta_lines`（検定・判定規則・BH の母集団・`p_adjust_method`・`q_value_method = bh` と `q_value_source`・単位・
  `matrix_id`、交互作用ならその定義など）が足される。前処理条件は登録した行列の設定から取る。
- `exclude_feature_ids` の行は検定の後に落とし（q は計算し直さない）、件数をメタ行 `n_excluded_after_test` に書く。
- 出力は一時ファイルへ書いてから置換する（途中で落ちた出力は残らない）。

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

読み取り専用（戻り値・`merge_split`・`rule_hash` と `stats_source` の意味は `docs/output_format/mztab.md` の 13.8）。`targets` の
有無で窓検索と目標照合に分かれ、どちらも最後に同じ統計の結合を通る。判定の純ロジックは `metabolomix/analysis/feature_search.py`
にあり、session も MCP も import しない。ツール層の仕事は、極性の出所を決めること（ジョブの宣言 → mzTab-M の付加体の多数決）、
差次的結果を選ぶこと（名指しの ID が無い・古ければ `missing_state`、省略時は今も有効な 2 群比較だけ）、丸めと TSV 化、注記。
キュレーションのフラグは `dataset_export_differential` と同じ見つけ方で数えるだけで**反映しない**（読めなくても検索は止めず、
注記にする）。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_feature_search()
2. ├─ [窓検索] metabolomix/analysis/feature_search.py  search_window()
3. │  └─ metabolomix/analysis/feature_search.py  _window_identity()
4. ├─ [目標照合] metabolomix/mztab/polarity.py  dataset_polarity()
5. │  ├─ metabolomix/console/job_manager.py  load_job()
6. │  └─ metabolomix/console/output_collector.py  read_adduct_polarity()
7. ├─ [目標照合] metabolomix/analysis/feature_search.py  normalize_targets()（種の目標に同位体の標識があれば `analysis/msms_rules.py` `isotope_labelled()` で見て `FEATURE_TARGET_INVALID`）
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

- `include_statistics=false` では `select_results` と統計の枝（`feature_stats` 〜 `result_summaries`）を通らない。
- `select_results` は `dataset_compare` の結果（`compare` の記録）を対比ごとに 2 群の形に広げる（対比が複数なら ID は
  `result_id#contrast_id`、交互作用は除く）。`feature_stats` の行の `detected_units_a` / `detected_units_b` / `presence_p` /
  `detection_driven` は `as_two_group` が結果の行から写したもの。
- `merge_split=true`: 目標照合は候補を数える前（`match_target` の前）に `merge_plan` を呼び、`match_target(merged_from=merge_map(...))`
  が組を 1 つの候補（残す ID）として数え、その組だけを戻り値の `merge_split` に写す。窓検索は統計を付けるときだけ呼ぶ。検定の
  やり直し（`retest_merged`）は、結果の spec で `stats/data.py` `from_dataset()` → `merge_feature_columns()` →
  `spec_v2.run_spec_v2()` → `as_two_group()` を通す。
- `estimate_rt_offset` は RT のずれを測るだけで判定には当てない（当てるのは `rt_offset_min` を渡したときの `match_target`）。
  目標照合の `rt_offset_min` は `normalize_targets` が `rt_min` のある目標にだけ足し、どの目標にも足されなければツールが `notes` に
  書く。
- 効かない引数は誤りにせず `notes` に書くだけ: 窓検索の `ppm_offset` / `rt_offset_min`（照合には通さず、形は検める）、目標照合の
  `sort="intensity"`、`merge_split=false` と一緒の `alignment_qc_path`（読まない）。
- `select_results` と、統計を付けるときの `merge_plan`（窓検索）・`merged_feature_stats` が送出した `DomainError` は、ツールが
  それぞれの呼び出しの直後に `contrast_envelope` で封筒にする（対比を省いた曖昧さなら `needs_input`〔`argument="contrast_id"`〕、
  未知の対比や `FEATURE_MERGE_UNAVAILABLE` などは `mztab_error`）。目標照合が照合の前に呼ぶ `merge_plan` の誤りは、目標の誤りと
  同じく直に `mztab_error`。
- `rule.targets` の表示用の丸め（理論 m/z を 4 桁）は `_rule_targets()` で、`rule_hash` はその前の値から作る。

## dataset_set_sample_metadata

実験情報シート(sample-manifest.v1 / v2)を読み、全件検証してから一括反映する。版は先頭の schema 行で決まり、列見出しからは
推測しない。検証で 1 件でも落ちれば `apply_metadata` は `session.dataset` を変えない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_set_sample_metadata()
2. └─ metabolomix/analysis/sample_manifest.py  apply_sample_manifest()
3.    └─ metabolomix/analysis/sample_manifest.py  _raw_manifest_layout()（シートの置き場所。試料ごとの `_raw_path_for_sample()` から取る＝手順 5 が行を試料へ結ぶのと同じ raw のパス。試料から引けないときだけ `ds.assay_sources`）
4.    └─ metabolomix/analysis/sample_manifest.py  parse_manifest()（任意の列 `subject_id` を読む。`sample_meta` に持ち、統計エンジンの前処理が反復注入を 1 個体にまとめる鍵に使う。schema の行の直後の宣言 `# injection_order_source = ...` は `_read_declarations()` が読み、注入順の出所に `derived_from` を足す）
5.    └─ metabolomix/analysis/sample_manifest.py  resolve_metadata()
6.    └─ metabolomix/analysis/sample_manifest.py  apply_metadata()
7.       └─ metabolomix/analysis/result_state.py  invalidate_results()
8. └─ metabolomix/analysis/stats/data.py  merged_feature_polarity()（陽陰を統合したデータセットなら戻り値に `notes`）
9.    └─ metabolomix/analysis/stats/preprocess.py  run_remedy()（陰イオンの測定を入れ直す手順の文面）

陽陰を統合したデータセット（意味は `docs/output_format/dataset.md` の 10.2）: シートは陽イオンの raw を挙げるので、`batch`・
`injection_order` は陽イオンの測定だけを決める。陰イオンの測定は `quant_state.polarity_runs.negative` のままで、変えるなら陰イオンの
ハンドルにシートを当ててから統合し直す。手順 8・9 が戻り値の `notes` にそう書き、統合し直すと統合したハンドルへ当てたシートと除外
（`dataset_exclude`）は引き継がれないので、同じシートを陽イオンのハンドルにも当ててから統合し直すか、できたハンドルに当て直し、
除外も当て直すとも書く。

## dataset_export_evidence_bundle

全特徴量と候補の根拠と、試料ごとの値と検出を v3 バンドルへ保存する（意味は `docs/output_format/reanalysis.md` の 10）。
15 列の TSV（`dataset_export_differential`）は変えない。`result_id` を省けば直近の差次的結果（`ds.last_differential`）。対比は
書き出しと同じ選び方（省いて 2 群の対比が複数なら `needs_input`）。`features.tsv` の `q_value` と `significant` は書き出しと同じく
契約の BH の q（`as_two_group(contract_q=True)`）で、`manifest.json` の `comparison.correction` はいつも `BH`。比較の補正が BH で
なければ `comparison.p_adjust_method` と `q_value_source` を足す。探索専用のデータセットも断らず、manifest の最上位と `comparison`
に `exploratory`（印と理由）を足す。

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

`link_alignment_spectra()` は spectrum_links を同じアラインメントの `.dcl` に ID（SMF_ID == AlignmentID == `.dcl` の番号）で結び、
m/z・RT は検証にだけ使う（結べないときの状態と理由は manifest の provenance）。`annotations.tsv` の `msdial_sme` 行は SME の候補
そのものだが、組成式と付加体から計算した理論 m/z（`DatasetState.theoretical_mz_computed` の特徴。`directory=` で開いたデータ
セット）は `library_candidate()` が空にして写す（`theoretical_mass_to_charge` はライブラリの記録のプリカーサー m/z のまま。
データセットの候補は変えない）。
