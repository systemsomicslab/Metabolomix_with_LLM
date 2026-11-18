# ワークフロー: 図（dataset_plot）・図の保存と payload 契約

`dataset_plot` が図を描き（既定は画像。`output_path` を渡せば同じ図をそのファイルにも書く）、記録を `ds.figures[kind]`
（kind ごとに最後の 1 つ。指紋に含めない）に置く。`save_figure` はそこから `reports/figures/<analysis_id>_<kind>.png` へ書く。
ファイルへの書き方は両者で `plot_tools.write_figure_files` を共有する。値の意味は
[dataset.md（output_format）](../output_format/dataset.md) の §11.2.1・§11.4.9〜11.8（以下の §11.N はその節）と
[eic.md（output_format）](../output_format/eic.md)。

## payload 契約の比較

生成元は `dataset_plot(kind=<kind>, output="payload")`、保存は `save_figure(kind=<kind>)`。

| payload | kind | 特徴 |
|---|---|---|
| `lipidmix.volcano.v1` | `volcano` | `up`/`down` は全点保持、`ns` のみ間引く。保存は間引き前の全点 |
| `lipidmix.pca_scores.v1` | `pca` | 散布図。点の `group` が `color_by`（group / batch / role）の値 |
| `lipidmix.pca_loadings.v1` | `loadings` | 主成分ごとの横棒。値は r か coefficient。色は `ontology`。PNG（dpi 300）と SVG |
| `lipidmix.eic.v1` | `eic` | 1 特徴 × 試料別の線グラフ |
| `lipidmix.eic.multi.v1` | `eic_compounds` | 1 試料 × 複数物質の重ね描き |
| `lipidmix.group_intensity.v1` | `group_intensity` | 項目（クラス・分子種・特徴 ID）× 群 × 試料の強度合計（log10 の平均 ± SD）。PNG（dpi 300）と SVG。群を焼き込むので試料の情報の更新で古くなる |
| `lipidmix.species_intensity.v1` | `species` | 分子種（特徴）ごとのパネル。`value="share"`（既定。分母は `share_basis`）か `"intensity"`（payload の `value` は `"height"`）。PNG（dpi 300）と SVG |
| `lipidmix.species_pca.v1` | `species_pca` | 選んだ分子種の PCA のスコア散布図（`fitted` が False の点は投影だけ）。ローディングの全行は結果（`ds.results[result_id]`）に持つ。PNG（dpi 300）と SVG |

既定は画像（`output="image"`）。返す PNG は `render_png` が画素数を `MAX_IMAGE_PIXELS`（40 万画素）に抑え、大きな図は dpi を
`MIN_IMAGE_DPI`（72）まで下げる（ファイルに書く図 `render_files` には上限を掛けない）。payload が要るクライアント（Use-LLLM）は
起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、呼び出しごとに `output="payload"` を渡す（判定は `plots/render.py`
`resolve_plot_output()`）。描画は `plots/` 側にあり、画像返しと保存が同じ関数を通る。

保存先は `mcp_core._resolve_report_dir(holder)`（`holder` は `dataset_handle` のハンドル、省けば直近）: `lipidmix.local.toml` の
`[paths] reports_dir` を優先し、書けなければそのデータセットの解析フォルダ（`DatasetState.analysis_dir`）の `reports/`。未指定なら
解析フォルダの `reports/` だけ（リポジトリ側へは退避しない。書けなければエラー）。設定もデータセットも無ければ
`missing_state("dataset", ["load_dataset"])`。

## dataset_plot

前提（`kind` ごと）:
- `volcano` は `dataset_compare`（対比を 2 群の形に読み替えて描く。2 群の対比が複数の結果は `contrast_id` か
  `result_id#contrast_id` が要り、省くと `needs_input`〔`tools/_contrast_choice.py` `contrast_envelope`〕。交互作用の対比は描けない。
  題は `test vs reference`）。
- `pca` / `loadings` は `dataset_pca`。
- `group_intensity` / `species` / `species_pca` は結果が要らない（`items` と `groups` を渡す。キュレーションの判断は
  `load_dataset` が読んだ `ds.curation_flags`）。
- `eic` / `eic_compounds` はデータセットの `.EIC.aef`（`artifact_paths` の `eic`）。
無ければ `missing_state`（`dataset_pca` / `dataset_compare`）か `COMPANION_ARTIFACT_MISSING`（`.EIC.aef` が無い）。

状態変更: `ds.figures[kind]`（メモリ上。描けた呼び出しだけ）。`output_path`（`.png`）を渡せば `write_figure_files` が同じ図を書く
（親のフォルダは作る。書けなければ `FIGURE_SAVE_FAILED`）。`species_pca` は結果を `ds.results` へ登録する（来歴の kind は
`species_pca`。`loadings` の `result_id` には使え、`pca` の `result_id` には使えない。来歴に `metadata_revision` を持ち、試料の情報の
更新・`dataset_exclude` で `is_current` が False〔`STALE_ANALYSIS_RESULT`〕。キュレーションの判断の更新では古くしない）。

- 結果の選び方（`_pick_result`）: `result_id` があれば `ds.results` から引き、現在のものでなければ `STALE_ANALYSIS_RESULT`
  （`result_state.assert_current`）。無ければ直近（`ds.last_pca` / `ds.last_differential`）。除外・前処理・群の変更で古くなった結果は
  `result_state.invalidate_results` が `last_*` と `ds.figures` の該当 kind を同時に捨て、捨てた kind に `ds.discarded_figures` の
  印を残す（`dataset_plot` が同じ kind を描き直すと消える）。
- 別の kind の引数は `_drop_unused_arguments` が外して注記にする（描けなかった誤りの封筒は `_with_notes` が `details.notes` に。
  どの封筒が持つかは §11.4.9）。置き場の図に注記は混ぜない。
- `group_intensity` / `species` / `species_pca` の行は `plots/dataset_rows.py` が `DatasetState` から選ぶ: 値は `feature_matrix`、
  gap-fill の印は `detected_mask`、群は `sample_metadata_rows` の `group` / `role`、残す試料は `included_sample_indices`
  （include=false ∪ `dataset_exclude`）、除外した特徴（`excluded_features`）はクラスにも名前にも入れない。MS/MS の裏付けは
  `feature_metadata["msms_assigned"]`（`.arf` のどれかのセルが MS/MS を持つ）と `.arf2` の Name の `no MS2:` / `w/o MS2:` の接頭辞。
  標準液にだけある分子種（`standard_samples`）は行列から判定。記録は `metadata_revision` を持つ。`species_pca` は前処理の入力・
  除外の変更で `invalidate_results` が `pca` と同じ組で捨てる。kind ごとの引数の持ち主は `_KIND_ARGS`。
- EIC は特徴 ID（= `.arf2` の MasterAlignmentID）を `.EIC.aef` のスポット番号として引く（イオンモビリティで連番にならないときは
  注記）。RT と m/z を特徴表と照合し、`eic` は合わなければ注記、`eic_compounds` は描かず `selection.dropped` に理由を残す。
  `file_ids` を省けば解析に残さない試料（`included_sample_indices`）を外す。

1. metabolomix/tools/plot_tools.py  dataset_plot()
2. ├─ metabolomix/tools/plot_tools.py  _drop_unused_arguments()
3. ├─ metabolomix/tools/_holder.py  resolve_holder()
4. ├─ [kind=volcano] metabolomix/tools/plot_tools.py  _volcano()
5. │  ├─ metabolomix/analysis/stats_results.py  split_result_ref()
6. │  ├─ metabolomix/tools/plot_tools.py  _pick_result()
7. │  ├─ metabolomix/analysis/stats_results.py  as_two_group()
8. │  └─ metabolomix/plots/volcano.py  build_volcano_plot_payload()
9. ├─ [kind=pca] metabolomix/tools/plot_tools.py  _pca()
10.│  └─ metabolomix/plots/pca_scores.py  build_pca_scores_payload()
11.├─ [kind=loadings] metabolomix/tools/plot_tools.py  _loadings()
12.│  ├─ metabolomix/plots/pca_loadings.py  dataset_features()
13.│  └─ metabolomix/plots/pca_loadings.py  build_loadings_payload()
14.├─ [kind=group_intensity] metabolomix/tools/plot_species_kinds.py  group_intensity()
15.│  ├─ metabolomix/tools/plot_species_kinds.py  _common()
16.│  │  ├─ metabolomix/plots/dataset_rows.py  resolve_samples()
17.│  │  ├─ metabolomix/plots/dataset_rows.py  standard_only_spots()
18.│  │  └─ metabolomix/plots/dataset_rows.py  resolve_selection()
19.│  ├─ metabolomix/plots/dataset_rows.py  rows_by_spot()
20.│  └─ metabolomix/plots/group_intensity.py  build_group_intensity_payload()
21.├─ [kind=species] metabolomix/tools/plot_species_kinds.py  species()
22.│  ├─ metabolomix/plots/dataset_rows.py  expand_spots()
23.│  └─ metabolomix/plots/species.py  build_species_payload()
24.├─ [kind=species_pca] metabolomix/tools/plot_species_kinds.py  species_pca()
25.│  ├─ metabolomix/plots/dataset_rows.py  expand_spots()
26.│  ├─ metabolomix/analysis/pca.py  run_pca_fit_subset()
27.│  └─ metabolomix/analysis/result_state.py  register_result()
28.├─ [kind=eic] metabolomix/tools/plot_tools.py  _eic()
29.│  ├─ metabolomix/dataset/artifacts.py  artifact_path()
30.│  ├─ metabolomix/eic/reader.py  read_eic_spot_css1()
31.│  └─ metabolomix/plots/eic.py  build_eic_plot_payload()
32.├─ [kind=eic_compounds] metabolomix/tools/plot_tools.py  _eic_compounds()
33.│  ├─ metabolomix/eic/identity_map.py  select_identity_candidates()
34.│  ├─ metabolomix/eic/reader.py  read_eic_spots_css1()
35.│  └─ metabolomix/plots/eic.py  build_multi_compound_plot_payload()
36.├─ [output_path 指定] metabolomix/tools/plot_tools.py  write_figure_files()
37.├─ [output=payload] metabolomix/core/serialization.py  json_payload()
38.├─ [output=image] metabolomix/tools/plot_tools.py  _png()
39.│  └─ metabolomix/plots/render.py  render_png()
40.│     └─ metabolomix/plots/render.py  image_dpi()
41.├─ [DomainError] metabolomix/tools/_contrast_choice.py  contrast_envelope()
42.└─ [誤りの封筒で、使わなかった引数がある] metabolomix/tools/plot_tools.py  _with_notes()

## save_figure

前提: `dataset_plot` が置いた図（`ds.figures[kind]`）か、`result_id` で登録済みの結果（pca / volcano / loadings。volcano は
`contrast_id` でも対比を選んで描き直せる）。一度も描いていなければ `missing_state`（`required_tools` は `dataset_plot`、
データセットが無ければ `load_dataset`）。置き場の図が古い（`_stale_figure`）か、`invalidate_results` が捨てた印
（`ds.discarded_figures`）があれば保存せず `STALE_FIGURE`（文言は `_discarded_message`。因子の列だけの更新は
`invalidate_results` の `reason` に広げる前の列が渡るので群を挙げない。`details` は §11.4.9）。
状態変更: `reports/figures/` に PNG（`loadings` `group_intensity` `species` `species_pca` は SVG も）を書く。上書きの前の版は
`reports/.history/` へ移す（`_write_keeping_history()`。書けなければ戻し、同じバイト列でも戻す）。

- `kind` は `dataset_plot` の kind（`FIGURE_KINDS`。`source` 引数は無い）。未知の `kind` は `{"status": "error"}`。`result_id` を使わない kind
  （`group_intensity` `species` `eic` 系）への `result_id` と、volcano 以外への `contrast_id` は使わずに置き場の図を保存し、
  `notes:` 行に書く（`_save_record`）。
- `species_pca` の `result_id` は描き直しに使わず、置き場の図（先に `_stale_figure` で確かめる）の結果と照らす: 同じなら保存、
  同じ計算（データセット・入力の行列の指紋・指定・試料の情報の版が同じ。`_same_species_pca`）の別の結果なら保存して `notes:` 行に
  書く。分子種の PCA の結果でなければ `DATASET_BAD_REQUEST`（このハンドルに無い ID は、ほかのハンドルの分子種の PCA の結果なら
  `details.result_handle` にそのハンドル、別の種類の結果ならそのハンドルを文で名指す。呼び直しは案内しない）。ほかは保存せず
  `STALE_FIGURE`（`_other_species_pca`。`details.held_result_id` と `details.mismatch`: 名指しした結果が `is_current` でなければ
  `stale`、指定が置き場の図と同じなら `input`、違えば `settings`。案内と `details.settings` は §11.7 の `result_id`）。
- `contrast_id` と `result_id#contrast_id` が食い違えば `DATASET_BAD_REQUEST`（`split_result_ref`）、2 群の対比が複数で選ばなければ
  `needs_input`（`contrast_envelope`）。
- `pca` / `volcano` は `save_result_figure()`（探索専用データセットの図には図の中に但し書き）、`loadings` `group_intensity`
  `species` `species_pca` は dpi 300 の PNG と SVG、`eic` 系は dpi 120 の PNG。
- 群を焼き込んだ図（`group_intensity` `species` `species_pca`）は `metadata_revision` が変わると拒み、`species_pca` は結果が現在かも
  確かめる（`_stale_figure`）。戻り値に保存した絶対パスを書く（Use-LLLM が画像として表示する）。

1. metabolomix/tools/reports.py  save_figure()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ [result_id（pca / volcano / loadings）/ contrast_id 指定] metabolomix/tools/plot_tools.py  build_figure()
4. ├─ [どちらも省略・置き場に無く、捨てた印がある] metabolomix/tools/reports.py  _discarded_message()
5. ├─ [どちらも省略] metabolomix/tools/reports.py  _stale_figure()
6. ├─ [kind=species_pca で result_id が置き場の図の結果と違う] metabolomix/tools/reports.py  _other_species_pca()
7. │  ├─ metabolomix/tools/reports.py  _same_species_pca()
8. │  └─ metabolomix/analysis/result_state.py  is_current()
9. ├─ [DomainError] metabolomix/tools/_contrast_choice.py  contrast_envelope()
10.└─ metabolomix/tools/reports.py  _save_record()
11.   └─ metabolomix/tools/reports.py  _write_keeping_history()
12.      ├─ metabolomix/tools/notes.py  keep_history()
13.      └─ metabolomix/tools/plot_tools.py  write_figure_files()
14.         ├─ [kind=pca / volcano] metabolomix/plots/result_output.py  save_result_figure()
15.         │  └─ [kind=pca] metabolomix/plots/pca_scores.py  render_pca_scores()
16.         ├─ [kind=loadings / group_intensity / species / species_pca] metabolomix/plots/render.py  render_files()
17.         │  ├─ metabolomix/plots/pca_loadings.py  render_loadings_plot()
18.         │  └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
19.         └─ [kind=eic / eic_compounds] metabolomix/plots/eic.py  render_eic_plot()
