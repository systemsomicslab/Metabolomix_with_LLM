# ワークフロー: 参照ライブラリ（MS/MS スペクトル照合）

`.dcl` の測定 MS/MS を参照ライブラリ（データセットが同定に使ったもの、または `file_path` で渡した `.msp` / `.dbs`）と
突き合わせ、MS-DIAL の個別スコア定義で採点し、対向プロットで確かめる経路。前提状態の連鎖は
`library_load` → `library_match_feature` → `library_plot_mirror`。前段が未実行なら `missing_state` 封筒が前段のツール名を指す。
`session.library`（他スロットとは独立）に `store`（開いたままの `LibraryStore`）と `last_match`（直近の照合結果。座標を
含む）を持つ。値の意味は [library.md（output_format）](../output_format/library.md)（以下の §14.N はその節）。

## library_load

前提: `file_path` を渡すか、直近のデータセット（または `dataset_handle`）が同定のライブラリ（`library_record`）を覚えて
いること（`verify_record` が通った `known` だけ）。探索は無い（`ion_mode` 引数・データフォルダの `.dbs` / `.msp`・設定の
極性の MSP〔`[library] msp_positive` / `msp_negative`。Console に使わせるものだけ〕への後退は無い）。記憶が `mismatch` なら
`LIBRARY_RECORD_MISMATCH`、`lbm_only` なら `LIBRARY_RECORD_LBM_ONLY`、`unknown` なら `LIBRARY_NOT_REMEMBERED`、データセットが
無ければ `missing_state("dataset", ["load_dataset"])`、明示の `file_path` が無ければ `LibraryPathError` を `LIBRARY_NOT_FOUND` の
エラー戻り値に（§14.5）。
状態変更: `session.library.store` / `source_path` を格納し、`last_match` を `None` に戻す（別ライブラリの照合結果と取り違えないため）。

store は元ファイルの sha256 をキーにした SQLite キャッシュで、同じファイルは作り直さない（`rebuild=True` で作り直す）。sha256 は
キャッシュ置き場の対応表（`.library-digests.json`）にサイズ・更新時刻と組で覚え、両方が変わらなければ測り直さない
（`_digest_for()`。`rebuild=True` では覚えた値を使わない）。`.msp` は 1 行ずつ読み、文字コードは行ごとに UTF-8 → cp932 → latin-1。
同じ構築を MCP の外で走らせる入口は `store.main()`（[cli.md](../cli.md) の `metabolomix.library.store`）。

1. metabolomix/library/tools.py  library_load()
2. ├─ [file_path を省いたとき] metabolomix/dataset/library_memory.py  verify_record()
3. └─ metabolomix/core/path_resolvers.py  resolve_library_path()
4. └─ metabolomix/library/store.py  open_store()
5. │  └─ metabolomix/library/store.py  _digest_for()
   │     └─ metabolomix/library/store.py  source_sha256()  [対応表に無い／変わった時、または rebuild=True]
6. │  └─ metabolomix/library/store.py  _build()  [未キャッシュ時、または rebuild=True]
7. │     ├─ metabolomix/library/store.py  _iter_records_for()
8. │     │  ├─ [.msp] metabolomix/library/msp.py  iter_records()
9. │     │  └─ [.dbs/.lbm2] metabolomix/library/dbs.py  iter_records()
10. │     └─ metabolomix/library/dbs.py  read_storage_meta()
11. └─ metabolomix/library/store.py  LibraryStore.summary()
12. └─ metabolomix/library/store.py  LibraryStore.compound_class_counts()

## library_match_feature

前提: `library_load` 実行済み（`session.library.store` が `None` なら手順 1 で `MissingState`）。
状態変更: `session.library.last_match` を更新（測定スペクトル・候補ごとの参照スペクトルと alignment。戻り値には座標を含めない）。

- 許容幅（`mz_tol` / `ms2_tol` / `rt_tol`）は `pick_tol()`（`metabolomix/library/defaults.py`）が明示 > store の
  `search_params`（`.dbs` 由来のときだけ）> 既定値の順で決め、`store.candidates()` にだけ効く。`.dcl` の測定 MS/MS を引く窓
  （`_measured_spectrum`）は固定の `mz_tol=0.01` / `rt_tol=0.2`（§14.4）。
- 照合対象は `precursor_mz` か `spot_id` のどちらか一方（両方・どちらも無しは `status="invalid_arguments"`）。`spot_id` では
  手順 2 で precursor m/z・RT・極性を `.arf2`（省けばデータセットの `.arf2`。データフォルダは探さない）から取り、測定スペクトルを
  アラインメント `.dcl` の `dcl[spot_id]` から直接読む（手順 8 の precursor 検索は通らない）。ARF2 と ±0.05 Da・±0.2 分で
  合わなければ `dcl_precursor_mismatch`、ID が無ければ `spot_not_found`。`.arf2` の極性が `Both` / `Unknown` なら絞らない
  （`_arf2_polarity()`）。
- `ion_mode` は手順 4 で正規化（`pos` / `neg` / `+` / `-` も可。読めない値は `invalid_arguments`）。
- `.dcl` に MS/MS が無ければ `status="not_found"`、候補が 0 件なら `status="no_candidates"`（例外にしない）。
- 照合の芯（手順 6〜）は `_run_match()` で、`library_plot_mirror(spot_id=…)` も呼ぶ。照合を試みて失敗したら（ファイルが
  無い・`spot_not_found`・`dcl_precursor_mismatch`・`not_found`・`no_candidates`）`_clear_last_match()` が `last_match` を捨てる。
- 前駆体の窓は手順 18 で MS-DIAL と同じく m/z 500 を超えたら同じ ppm に広げる（`query.mz_window`）。

1. metabolomix/library/tools.py  library_match_feature()
2. └─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ [ライブラリ無しで記憶が lbm_only のとき] metabolomix/library/tools.py  lbm_only_not_evaluable()（`not_evaluable` を返す）
4. └─ metabolomix/library/matching.py  normalize_ion_mode()
5. │  └─ metabolomix/core/polarity_names.py  normalize_polarity()
6. └─ metabolomix/library/tools.py  _run_match()
7. │  └─ metabolomix/tools/_holder.py  alignment_paths()（省いた `.arf2` / `.dcl` はデータセットのもの。無ければ `missing_state` / `no_alignment_files`）
8. │  └─ metabolomix/library/tools.py  _spot_query()
9. │  │  └─ metabolomix/msdial/alignment_spot.py  spot_record()
10. │  │  └─ metabolomix/msdial/alignment_spot.py  alignment_msms()
11. │  │  │  └─ metabolomix/dcl/reader.py  read_dcl_entry()
12. │  │  └─ metabolomix/library/tools.py  _arf2_polarity()
13. │  └─ metabolomix/library/tools.py  _measured_spectrum()
14. │  │  └─ metabolomix/dcl/reader.py  deserialize_dcl()
15. │  │  └─ metabolomix/dcl/reader.py  get_msms_by_precursor()
16. │  └─ metabolomix/library/matching.py  score_candidates()
17. │  │  └─ metabolomix/library/defaults.py  pick_tol()
18. │  │  └─ metabolomix/analysis/spectral_match.py  fix_mass_tolerance()（前駆体の窓）
19. │  │  └─ metabolomix/library/store.py  LibraryStore.candidates()
20. │  │  └─ metabolomix/analysis/spectral_match.py  cutoff_mask()
21. │  │  └─ metabolomix/analysis/spectral_match.py  match_spectrum()
22. │  │  │  └─ metabolomix/analysis/spectral_match.py  normalize_measured()
23. │  │  │  │  └─ metabolomix/analysis/spectral_match.py  cutoff_mask()
24. │  │  │  └─ metabolomix/analysis/spectral_match.py  simple_dot_product()
25. │  │  │  └─ metabolomix/analysis/spectral_match.py  weighted_dot_product()
26. │  │  │  └─ metabolomix/analysis/spectral_match.py  reverse_dot_product()
27. │  │  │  └─ metabolomix/analysis/spectral_match.py  matched_peaks_scores()
28. │  │  │  └─ metabolomix/analysis/spectral_match.py  spectral_entropy_similarity()
29. │  │  └─ metabolomix/analysis/spectral_match.py  total_score()
30. │  │  │  └─ metabolomix/analysis/spectral_match.py  gaussian_similarity()
31. │  │  │  └─ metabolomix/analysis/spectral_match.py  fix_mass_tolerance()
32. │  │  └─ metabolomix/library/matching.py  dbs_cutoffs()
33. │  └─ [失敗したとき] metabolomix/library/tools.py  _clear_last_match()
34. │  └─ [spot_id の照合で注釈があるとき] metabolomix/library/matching.py  msdial_agreement()
35. │     └─ metabolomix/library/matching.py  _same_compound()
36. │        └─ metabolomix/analysis/msms_rules.py  record_lipid_species()（MS-DIAL の注釈の名前は `name_lipid_species()`）
37. └─ metabolomix/library/record.py  inchikey_standard()
38. └─ [detail=true] metabolomix/analysis/spectral_match.py  fragment_matches_excluding_precursor()
39. └─ metabolomix/arf2/reader.py  format_spots_as_table()
40. └─ metabolomix/library/tools.py  identified_by_lbm()
41.    └─ metabolomix/dataset/library_memory.py  is_lbm_annotator()  [記憶が lbm_also で spot_id の照合のとき]

手順 20 の `cutoff_mask()` は候補ループの外で 1 回だけ呼び（測定スペクトルと `.dbs` の足切りだけで決まる）、採点に
入らなかったピークの m/z を `last_match["unscored_mz"]` に持つ（`library_plot_mirror` が区別して描く）。候補は `total_score` の
降順（`_RANK_KEY`）。RT 項を足すかは store の `search_params["use_time_for_annotation_scoring"]`（`.dbs` の `Key(16)`。`.msp` は
`False`）で、戻り値の `scoring` に出る（§14.9）。TSV の列は汎用の整形ヘルパ `arf2/reader.py` の `format_spots_as_table()` で組む。

## library_plot_mirror

前提: `library_match_feature` の直近結果（未実行、または直前の照合が失敗して捨てられたなら手順 1 で `MissingState`）。
`spot_id` を渡したときはライブラリだけが要る（手順 2 で照合し直す）。
状態変更: `spot_id` を渡したときだけ `last_match` を書き換える（全候補。照合に失敗したら捨てて、その封筒を返す）。

- `candidate`（InChIKey か名前）は手順 3 で描く候補を選ぶ（`rank` より優先。選び方は `match.picked_by`）。図の主は
  caption（画像）と `match`（payload）に出る（§14.7）。
- `ms2_tol` と `unscored_mz` は `last_match` から引き継ぐ（足切りで外れたピークは灰色で薄く描き、凡例と caption に件数）。
- 既定は画像（PNG）、`output="payload"` で座標 JSON（`lipidmix.mirror.v2`）。`scale`（`"relative"` / `"sqrt"` / `"log10"`）は
  `render_mirror()` の中で `scale_intensity()` が各点に掛ける（payload の座標は正規化前の生値）。
- `label_policy`（`"auto"` / `"msdial"`）は `_draw_labels()` がスペクトル全点を強度降順に走査して重ならないものだけ描く
  （描画時だけ。`payload["labels"]` は画像には使わない）。

1. metabolomix/library/tools.py  library_plot_mirror()
2. └─ [spot_id] metabolomix/library/tools.py  _run_match()（`library_match_feature` の手順 6 と同じ。`top_n=-1`）
3. └─ [candidate] metabolomix/library/tools.py  _pick_candidate()
4. └─ metabolomix/plots/render.py  resolve_plot_output()
5. └─ metabolomix/plots/mirror.py  build_mirror_payload()
6. └─ [output=image] metabolomix/plots/mirror.py  render_mirror()
7.    └─ metabolomix/plots/mirror.py  scale_intensity()
8.    └─ metabolomix/plots/mirror.py  _draw_labels()
9.    │  └─ metabolomix/plots/mirror.py  _text_extent()
10.    │  └─ metabolomix/plots/mirror.py  _is_overlap()
11.    └─ metabolomix/plots/render.py  figure_to_png()

## library_search

前提: `library_load` 済み（ハンドルならライブラリつきで開いたこと）。無ければ `missing_state`。
状態変更: なし（`SELECT` だけ。ハンドルの store は `mode=ro`）。

1. metabolomix/library/tools.py  library_search()
2. └─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/library/matching.py  normalize_ion_mode()（読めない `ion_mode` は `invalid_arguments`）
4. └─ metabolomix/library/store.py  LibraryStore.search()
5. └─ metabolomix/arf2/reader.py  format_spots_as_table()
