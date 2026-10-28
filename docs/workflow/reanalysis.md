# ワークフロー: 再解析のデータセットのハンドルと spec のツール

**再解析で使う**経路（spec 2026-10-22）。並行の再解析のエージェントが、対話と同じ MCP の
ツールを**データセットのハンドル**と宣言的な spec で呼ぶ。対話の使い方
（`dataset_load` → `dataset_*`）は変わらない。値の**意味**は
[reanalysis.md（output_format）](../output_format/reanalysis.md) を見ること。

ハンドルの状態は `metabolomix/core/handles.py` の登録簿（`handles.registry`）に住み、
`session_state.session` とは共有しない。`DatasetHandle` は `dataset` と `library` を
`AnalysisSession` と同じ名前で持つので、既存のツールは `holder_for(dataset_handle)` が
返したものを読むだけでハンドルにも対話のセッションにも効く。

`dataset_handle` を渡した呼び出しは `metabolomix/core/handle_dispatch.py` の
`install()` が包んだ関数がワーカースレッドで回し、そのハンドルの `lock` で直列化する。
渡さない呼び出しは今までどおりイベントループのスレッドで同期に呼ぶ。

## dataset_handle_open

前提: なし。引数は `dataset_load` と同じ（3 つのパスのうち 1 つ）。
状態変更: `handles.registry` にハンドルを 1 本足す（対話のセッションは変えない）。
読み込み済みの基を `HandleRegistry.put_base` で覚え、同じ入力（size・mtime が同じ）なら
読み直さずに複製を渡す。ライブラリは極性の設定だけから引き、読み取り専用で開く。

1. metabolomix/tools/handle_tools.py  dataset_handle_open()
2. ├─ metabolomix/core/handles.py  HandleRegistry.cached_base()
3. ├─ metabolomix/tools/handle_tools.py  _load()
4. │  ├─ metabolomix/pipeline/session_handoff.py  load_run_dataset()
5. │  └─ metabolomix/mztab/loading.py  load_dataset_state()
6. ├─ metabolomix/core/handles.py  HandleRegistry.put_base()
7. ├─ metabolomix/mztab/polarity.py  dataset_polarity()
8. ├─ metabolomix/tools/handle_tools.py  _alignment_files()
9. ├─ metabolomix/tools/handle_tools.py  _open_library()
10. │  ├─ metabolomix/core/user_config.py  get_setting()
11. │  └─ metabolomix/library/store.py  open_store()
12. └─ metabolomix/core/handles.py  HandleRegistry.add()

## dataset_handle_list

前提: なし。状態変更: なし。

1. metabolomix/tools/handle_tools.py  dataset_handle_list()
2. └─ metabolomix/core/handles.py  HandleRegistry.list()

## dataset_handle_close

前提: そのハンドルが開いていること（無ければ `DATASET_HANDLE_NOT_FOUND`）。
状態変更: 登録簿から外し、ライブラリの接続を閉じる。

1. metabolomix/tools/handle_tools.py  dataset_handle_close()
2. └─ metabolomix/core/handles.py  HandleRegistry.close()
3.    └─ metabolomix/core/handles.py  DatasetHandle.close()

## dataset_analysis_run

前提: `input.source=dataset` ならデータセット（ハンドルか対話のセッション。無ければ
`missing_state`）。外部の行列（`input.source=tsv`）ならデータセットは要らない。
状態変更: なし（読むだけ）。`output_path` を渡すと結果の JSON を書く。

1. metabolomix/tools/reanalysis_tools.py  dataset_analysis_run()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/analysis/analysis_spec.py  normalize_spec()
4. │  └─ metabolomix/analysis/recipe_matrix.py  normalize_recipe()
5. ├─ metabolomix/analysis/analysis_spec.py  run_spec()
6. │  ├─ metabolomix/analysis/analysis_spec.py  _load()
7. │  ├─ metabolomix/analysis/analysis_spec.py  _run_analysis()
8. │  │  ├─ metabolomix/analysis/recipe_matrix.py  build_matrix()
9. │  │  │  └─ metabolomix/analysis/recipe_matrix.py  knn_impute()
10. │  │  ├─ metabolomix/analysis/analysis_spec.py  _aggregate()
11. │  │  ├─ metabolomix/analysis/analysis_spec.py  _design()
12. │  │  ├─ metabolomix/analysis/spec_statistics.py  two_group()
13. │  │  ├─ metabolomix/analysis/spec_statistics.py  ols_group()
14. │  │  ├─ metabolomix/analysis/spec_statistics.py  fold_change()
15. │  │  └─ metabolomix/analysis/spec_statistics.py  bh()
16. │  └─ metabolomix/analysis/analysis_spec.py  rule_hash()
17. └─ metabolomix/tools/reanalysis_tools.py  write_result()

## library_msms_qualify

前提: ライブラリ（`library_load` か、ライブラリつきの `dataset_handle_open`）。無ければ
`missing_state`。状態変更: なし（直近の照合 `last_match` も書かない）。

1. metabolomix/tools/reanalysis_tools.py  library_msms_qualify()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/analysis/msms_rules.py  normalize_rules()
4. ├─ metabolomix/library/tools.py  _spot_query()
5. ├─ metabolomix/library/tools.py  _measured_spectrum()
6. ├─ metabolomix/library/matching.py  score_candidates()
7. ├─ metabolomix/analysis/msms_rules.py  qualify()
8. │  ├─ metabolomix/analysis/spectral_match.py  fragment_matches_excluding_precursor()
9. │  ├─ metabolomix/analysis/msms_rules.py  _lipid_class()
10. │  │  └─ metabolomix/msdial/lipid_identity.py  normalize_lipid_name()
11. │  └─ metabolomix/analysis/msms_rules.py  find_ions()
12. └─ metabolomix/analysis/msms_rules.py  rules_hash()

## dataset_homolog_check

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_homolog_check()
2. ├─ metabolomix/analysis/homolog.py  normalize_spec()
3. ├─ metabolomix/analysis/homolog.py  members_from_names()
4. │  └─ metabolomix/analysis/homolog.py  species_of()
5. │     └─ metabolomix/msdial/lipid_identity.py  normalize_lipid_name()
6. ├─ metabolomix/analysis/homolog.py  check()
7. └─ metabolomix/analysis/homolog.py  rule_hash()

## dataset_feature_relations

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_feature_relations()
2. ├─ metabolomix/analysis/claim_evidence.py  feature_relations()
3. │  ├─ metabolomix/analysis/claim_evidence.py  formula_index()
4. │  ├─ metabolomix/curation/relations.py  profile_correlation()
5. │  └─ metabolomix/curation/relations.py  find_relations()
6. └─ metabolomix/tools/reanalysis_tools.py  rows_table()

## dataset_claim_evidence

前提: データセット（無ければ `missing_state`）。`msms` を頼むときはライブラリ（無ければ候補ごとに
`no_library`）。状態変更: なし。`output_path` を渡すと JSON を書く。

1. metabolomix/tools/reanalysis_tools.py  dataset_claim_evidence()
2. ├─ metabolomix/analysis/claim_evidence.py  normalize_claim_spec()
3. │  ├─ metabolomix/analysis/feature_search.py  normalize_targets()
4. │  └─ metabolomix/analysis/msms_rules.py  normalize_rules()
5. ├─ metabolomix/analysis/claim_evidence.py  build()
6. │  ├─ metabolomix/analysis/feature_search.py  match_target()
7. │  ├─ metabolomix/analysis/claim_evidence.py  _msms_for()
8. │  │  ├─ metabolomix/library/tools.py  _spot_query()
9. │  │  ├─ metabolomix/library/matching.py  score_candidates()
10. │  │  └─ metabolomix/analysis/msms_rules.py  qualify()
11. │  ├─ metabolomix/analysis/claim_evidence.py  feature_relations()
12. │  └─ metabolomix/analysis/claim_evidence.py  rule_hash()
13. └─ metabolomix/tools/reanalysis_tools.py  write_result()

## dataset_mass_accuracy

前提: データセット（無ければ `missing_state`）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_mass_accuracy()
2. ├─ metabolomix/analysis/mass_accuracy.py  sme_ppm()
3. └─ metabolomix/analysis/mass_accuracy.py  band_summary()
4.    └─ metabolomix/analysis/mass_accuracy.py  summarize()

## dataset_novel_candidates

前提: データセット（無ければ `missing_state`）。`statistics.result_id` なら今も有効な 2 群比較
（無ければ `missing_state` が `dataset_differential` を指す）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_novel_candidates()
2. ├─ metabolomix/analysis/mass_accuracy.py  normalize_novel_spec()
3. ├─ metabolomix/analysis/feature_search.py  select_results()
4. ├─ metabolomix/analysis/analysis_spec.py  run_spec()
5. ├─ metabolomix/analysis/mass_accuracy.py  novel_candidates()
6. │  └─ metabolomix/analysis/mass_accuracy.py  sme_ppm()
7. └─ metabolomix/analysis/mass_accuracy.py  rule_hash()

## dataset_alignment_qc

前提: データセット（無ければ `missing_state`）。状態変更: なし。`output_path` を渡すと JSON を書く。
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
10. │  └─ metabolomix/analysis/alignment_qc.py  _components()
11. ├─ metabolomix/tools/reanalysis_tools.py  write_result()
12. └─ metabolomix/tools/reanalysis_tools.py  rows_table()

まとめは別のツールの側で行う: `alignment_qc.merge_matrix_columns()`（試料ごとに大きい方の値）と
`alignment_qc.split_partners()`（特徴 → 割れの相手）が唯一の実装。
