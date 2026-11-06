# ワークフロー: 再解析のデータセットのハンドルと spec のツール

**再解析で使う**経路（spec 2026-10-22）。並行の再解析のエージェントが、対話と同じ MCP の
ツールを**データセットのハンドル**と宣言的な spec で呼ぶ。対話の使い方
（`dataset_load` → `dataset_*`）は変わらない。値の**意味**は
[reanalysis.md（output_format）](../output_format/reanalysis.md) を見ること。

ハンドルの状態は `metabolomix/core/handles.py` の登録簿（`handles.registry`）に住む。
`dataset_load` で読んだデータセットも同じ登録簿のハンドルで、`session_state.session` は
その直近のハンドル（`current_handle_id`）を指すだけ（`session.dataset` / `session.library` は
property）。`DatasetHandle` は `dataset` と `library` を `AnalysisSession` と同じ名前で持つので、
既存のツールは `holder_for(dataset_handle)` が返したものを読むだけでハンドルにも対話の
セッションにも効く。

`dataset_handle` を渡した呼び出しは `metabolomix/core/handle_dispatch.py` の
`install()` が包んだ関数がワーカースレッドで回し、そのハンドルの `lock` で直列化する。
渡さない呼び出しも、直近のハンドルがあればそのハンドルの `lock` の下でワーカースレッドで回し、
`session_state.pinned_handle()` でそのハンドルに固定する。直近のハンドルが無ければ
イベントループのスレッドで同期に呼ぶ。

## dataset_handle_open

前提: なし。引数は `dataset_load` と同じ（`directory` / `mztab_path` / `job_path` /
`pipeline_path` のどれか 1 つ・`source`・`measure`・`allow_incomplete`）。
状態変更: `handles.registry` にハンドルを 1 本足す（対話のセッションは変えない）。
読み込みは `dataset_load` と同じ入口 `open_dataset` が持ち、読み込み済みの基を
`HandleRegistry.put_base` で覚え、同じ入力（読んだ・探したファイルの size・mtime が同じ）なら
読み直さずに複製を渡す。
ライブラリは極性の設定だけから引き、読み取り専用で開く。スポットの照合に使う `.arf2` / `.dcl` は、
解析フォルダで読んだならその選択（`AlignmentSelection`）のもの、mzTab-M で読んだなら
`_alignment_files()` が探す。

1. metabolomix/tools/handle_tools.py  dataset_handle_open()
2. ├─ metabolomix/dataset/opening.py  open_dataset()（連鎖は [mztab.md](mztab.md) の dataset_load）
3. │  ├─ metabolomix/core/handles.py  HandleRegistry.cached_base()
4. │  └─ metabolomix/core/handles.py  HandleRegistry.put_base()
5. ├─ metabolomix/mztab/polarity.py  dataset_polarity()
6. ├─ [mzTab-M で読んだとき] metabolomix/tools/handle_tools.py  _alignment_files()
7. │  └─ [pipeline_path のとき] metabolomix/dataset/opening.py  pipeline_arf2()
8. ├─ metabolomix/tools/handle_tools.py  _open_library()
9. │  ├─ metabolomix/core/user_config.py  get_setting()
10. │  └─ metabolomix/library/store.py  open_store()
11. ├─ metabolomix/dataset/library_memory.py  public_view()
12. └─ metabolomix/core/handles.py  HandleRegistry.add()

## dataset_handle_list

前提: なし。状態変更: なし。`current` は `session.current_handle_id` と同じハンドルか。

1. metabolomix/tools/handle_tools.py  dataset_handle_list()
2. └─ metabolomix/core/handles.py  HandleRegistry.list()

## dataset_handle_close

前提: そのハンドルが開いていること（無ければ `DATASET_HANDLE_NOT_FOUND`）。
状態変更: 登録簿から外し、ライブラリの接続を閉じる。直近のハンドルなら、対話のセッションは
データセットが無い状態になる（次に `session.dataset` を読んだときに `current_handle_id` を
None に戻す。別のハンドルへは落ちない）。

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
7. │  ├─ metabolomix/analysis/analysis_spec.py  _merge_plan()
8. │  │  ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
9. │  │  ├─ metabolomix/analysis/alignment_qc.py  read_qc_file()
10. │  │  ├─ metabolomix/analysis/alignment_qc.py  usable_merge_groups()
11. │  │  └─ metabolomix/analysis/alignment_qc.py  merge_matrix_columns()
12. │  ├─ metabolomix/analysis/analysis_spec.py  _run_analysis()
13. │  │  ├─ metabolomix/analysis/recipe_matrix.py  build_matrix()
14. │  │  │  └─ metabolomix/analysis/recipe_matrix.py  knn_impute()
15. │  │  ├─ metabolomix/analysis/analysis_spec.py  _aggregate()
16. │  │  ├─ metabolomix/analysis/analysis_spec.py  _design()
17. │  │  ├─ metabolomix/analysis/spec_statistics.py  two_group()
18. │  │  ├─ metabolomix/analysis/spec_statistics.py  ols_group()
19. │  │  ├─ metabolomix/analysis/spec_statistics.py  fold_change()
20. │  │  └─ metabolomix/analysis/spec_statistics.py  bh()
21. │  └─ metabolomix/analysis/analysis_spec.py  rule_hash()
22. └─ metabolomix/tools/reanalysis_tools.py  write_result()

割れたスポットのまとめ（`features.merge_split`、既定 true）は recipe の前（行列を読んだ直後。recipe・集約・検定はまとめた行列に当たる）に、
`alignment_qc` の `merge_groups` だけに行う（`features.alignment_qc` に JSON を渡せばその群、
無ければその場の `find_split_pairs`。JSON は `dataset_fingerprint` が今のデータセットと合うものだけ）。外部の行列（`input.source=tsv`）では何もしない。

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
6. │  ├─ metabolomix/analysis/alignment_qc.py  find_split_pairs()
7. │  ├─ metabolomix/analysis/alignment_qc.py  split_partners()
8. │  ├─ metabolomix/analysis/feature_search.py  match_target()
9. │  ├─ metabolomix/analysis/claim_evidence.py  _msms_for()
10. │  │  ├─ metabolomix/library/tools.py  _spot_query()
11. │  │  ├─ metabolomix/library/matching.py  score_candidates()
12. │  │  └─ metabolomix/analysis/msms_rules.py  qualify()
13. │  ├─ metabolomix/analysis/claim_evidence.py  feature_relations()
14. │  ├─ metabolomix/analysis/claim_evidence.py  _selection()
15. │  └─ metabolomix/analysis/claim_evidence.py  rule_hash()
16. └─ metabolomix/tools/reanalysis_tools.py  write_result()

候補ごとの `split` は `find_split_pairs`（既定の閾値）の `split_partners`。`default_selection.counts.qualified`
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

## dataset_novel_candidates

前提: データセット（無ければ `missing_state`）。`statistics.result_id` なら今も有効な 2 群比較
（無ければ `missing_state` が `dataset_differential` を指す）。状態変更: なし。

1. metabolomix/tools/reanalysis_tools.py  dataset_novel_candidates()
2. ├─ metabolomix/analysis/mass_accuracy.py  normalize_novel_spec()
3. ├─ metabolomix/analysis/feature_search.py  select_results()
4. ├─ [result_id] metabolomix/tools/reanalysis_tools.py  _novel_result_stats()
5. │  ├─ metabolomix/analysis/feature_search.py  merge_plan()
6. │  └─ metabolomix/analysis/feature_search.py  merged_results()
7. │     ├─ metabolomix/analysis/feature_search.py  merged_dataset()
8. │     └─ metabolomix/analysis/feature_search.py  retest_merged()
9. ├─ [analysis_spec] metabolomix/analysis/analysis_spec.py  run_spec()
10. ├─ metabolomix/analysis/mass_accuracy.py  novel_candidates()
11. │  └─ metabolomix/analysis/mass_accuracy.py  sme_ppm()
12. └─ metabolomix/analysis/mass_accuracy.py  rule_hash()

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
10. │  ├─ metabolomix/analysis/alignment_qc.py  _mutual_nearest()
11. │  ├─ metabolomix/analysis/alignment_qc.py  _components()
12. │  │  └─ metabolomix/analysis/alignment_qc.py  id_order()
13. │  └─ metabolomix/analysis/result_state.py  dataset_fingerprint()
14. ├─ metabolomix/tools/reanalysis_tools.py  write_result()
15. └─ metabolomix/tools/reanalysis_tools.py  rows_table()

まとめは別のツールの側で行う: `alignment_qc.merge_matrix_columns()`（試料ごとに大きい方の値）と
`alignment_qc.split_partners()`（特徴 → 割れの相手）が唯一の実装。
