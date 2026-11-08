# ワークフロー: 図（dataset_plot）・図の保存と payload 契約

データセットの図は `dataset_plot` が描く（既定は画像。ファイルは書かない）。作った図の記録は
`ds.figures[kind]`（kind ごとに最後の 1 つ。指紋に含めない）に置かれ、ユーザーが明示的に保存を
求めたときだけ `save_figure` がそこから `reports/figures/<analysis_id>_<kind>.png` に書く。
旧 ARF 系の描画ツール（`arf_plot_*` `dataset_plot(kind="species_pca")` `dataset_plot(kind="loadings")` `eic_plot_*`）は B-2b Task 8 で消えた。

## payload 契約の比較

| payload | 生成元 | 保存 | 特徴 |
|---|---|---|---|
| `lipidmix.volcano.v1` | `dataset_plot(kind="volcano", output="payload")` | `save_figure(kind="volcano")` | `up`/`down` は全点保持、`ns` のみ間引く。保存は間引き前の全点 |
| `lipidmix.pca_scores.v1` | `dataset_plot(kind="pca", output="payload")` | `save_figure(kind="pca")` | 散布図。点の `group` が `color_by`（group / batch / role）の値 |
| `lipidmix.pca_loadings.v1` | `dataset_plot(kind="loadings", output="payload")` | `save_figure(kind="loadings")` | 主成分ごとの横棒。値は r か coefficient。色は `ontology`。PNG（dpi 300）と SVG |
| `lipidmix.eic.v1` | `dataset_plot(kind="eic", output="payload")` | `save_figure(kind="eic")` | 1 特徴 × 試料別の線グラフ |
| `lipidmix.eic.multi.v1` | `dataset_plot(kind="eic_compounds", output="payload")` | `save_figure(kind="eic_compounds")` | 1 試料 × 複数物質の重ね描き |
| `lipidmix.group_intensity.v1` | `dataset_plot(kind="group_intensity", output="payload")` | `save_figure(kind="group_intensity")` | 項目（クラス・分子種・特徴 ID）× 群 × 試料の強度合計（log10 の平均 ± SD）。PNG（dpi 300）と SVG。群を焼き込むので試料の情報の更新で古くなる |
| `lipidmix.species_intensity.v1` | `dataset_plot(kind="species", output="payload")` | `save_figure(kind="species")` | 分子種（特徴）ごとのパネル。`value="share"`（既定。分母は `share_basis`）か `"intensity"`（payload の `value` は `"height"`）。PNG（dpi 300）と SVG |
| `lipidmix.species_pca.v1` | `dataset_plot(kind="species_pca", output="payload")` | `save_figure(kind="species_pca")` | 選んだ分子種の PCA のスコア散布図（`fitted` が False の点は投影だけ）。ローディングの全行は結果（`ds.results[result_id]`）に持つ。PNG（dpi 300）と SVG |

描画ツールの**既定は payload ではなく画像**（`output="image"`）。座標点列を LLM の文脈へ流すと
実測で数万トークンかかるのに対し、サーバ側で描いた PNG は 300〜1,000 画像トークンに収まる。
Plotly で描くクライアント（Use-LLLM）は起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、
呼び出しごとに `output="payload"` を渡す。判定は `plots/render.py resolve_plot_output()`。
描画は `plots/` 側に置き、画像返しと保存が同じ関数を通る（画面の図とレポートの図がずれない）。

保存先は `mcp_core._resolve_report_dir(holder)` が解決する（`holder` は `dataset_handle` で選んだハンドル、
省略時は直近）。`lipidmix.local.toml` の `[paths] reports_dir` の明示先を優先し、書き込めなければ
そのデータセットの解析フォルダ（`DatasetState.analysis_dir`）配下の `reports/` に退避する。未指定時は解析フォルダの
`reports/` だけを使い、リポジトリ側へは退避しない（書けなければエラー）。設定も無く、データセットも無ければ
`missing_state("dataset", ["load_dataset"])`。

## dataset_plot

前提: `kind` ごと——`volcano` は `dataset_compare`（2 群）、`pca` / `loadings` は `dataset_pca`、
`group_intensity` / `species` / `species_pca` は結果を要らない（`items` と `groups` を渡す。キュレーションの
判断は `load_dataset` が読んだ `ds.curation_flags`）、`eic` / `eic_compounds` はデータセットの `.EIC.aef`（`artifact_paths` の `eic`）。無ければ
`missing_state`（`dataset_pca` / `dataset_compare`）か `COMPANION_ARTIFACT_MISSING`（`.EIC.aef` が無い）。
状態変更: `ds.figures[kind]` を更新（メモリ上。ファイルは書かない）。`species_pca` はさらに結果を
`ds.results` へ登録する（来歴の kind は `species_pca`。`loadings` の `result_id` には使え、分子種のローディングを描く。
`pca` の `result_id` には使えない）。この結果は来歴に `metadata_revision` を持ち、試料の情報の更新・`dataset_exclude` で
`is_current` が False になる（`STALE_ANALYSIS_RESULT`）。キュレーションの判断の更新では古くしない（旧 `dataset_plot(kind="species_pca")` と同じ）。

結果の選び方（`_pick_result`）: `result_id` があれば `ds.results` から引き、現在のものでなければ
`STALE_ANALYSIS_RESULT`（`result_state.assert_current`）。無ければ直近（`ds.last_pca` /
`ds.last_differential`）。`dataset_exclude` や前処理・群の変更で古くなった結果は
`result_state.invalidate_results` が `last_*` と `ds.figures` の該当 kind を同時に捨てるので、
古い図は最新として描かれず保存もされない。別の kind の引数が渡されたら `DATASET_BAD_REQUEST`
（`_check_arguments`。どの kind の引数かを書く）。

`group_intensity` / `species` / `species_pca` の項目・行は `plots/dataset_rows.py` が `DatasetState` から選ぶ:
値は `feature_matrix`、gap-fill の印は `detected_mask`、群は試料の対応表（`sample_metadata_rows`）の
`group` / `role`、残す試料は `included_sample_indices`（include=false ∪ `dataset_exclude`）、除外した特徴
（`excluded_features`）はクラスにも名前にも入れない。MS/MS の裏付けは `feature_metadata["msms_assigned"]`
（`.arf` のどれかのセルが MS/MS を持つ）と、`.arf2` の Name の `no MS2:` / `w/o MS2:` の接頭辞。標準液にだけ
ある分子種（`standard_samples`）は行列から判定する。これらの図は群・役割を図に焼き込むので記録は
`metadata_revision` を持ち、試料の情報が更新されると `save_figure` が拒む。`species_pca` は結果から作った
図でもあるので、前処理の入力・除外の変更で `invalidate_results` が `pca` と同じ組で捨てる。
種類の違う 3 つ（kind ごとの引数の持ち主は `_KIND_ARGS`）。

EIC は特徴 ID（= `.arf2` の MasterAlignmentID）をそのまま `.EIC.aef` のスポット番号として引く
（イオンモビリティで連番にならないときは注記）。RT と m/z を特徴表と照合し、`eic` は合わなければ
注記、`eic_compounds` は描かず `selection.dropped` に理由を残す。`file_ids` を省いたときは、解析に残さない
試料（`included_sample_indices`）を外す。

1. metabolomix/tools/plot_tools.py  dataset_plot()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/tools/plot_tools.py  _check_arguments()
4. ├─ [kind=volcano] metabolomix/tools/plot_tools.py  _volcano()
5. │  ├─ metabolomix/tools/plot_tools.py  _pick_result()
6. │  └─ metabolomix/plots/volcano.py  build_volcano_plot_payload()
7. ├─ [kind=pca] metabolomix/tools/plot_tools.py  _pca()
8. │  └─ metabolomix/plots/pca_scores.py  build_pca_scores_payload()
9. ├─ [kind=loadings] metabolomix/tools/plot_tools.py  _loadings()
10.│  ├─ metabolomix/plots/pca_loadings.py  dataset_features()
11.│  └─ metabolomix/plots/pca_loadings.py  build_loadings_payload()
12.├─ [kind=group_intensity] metabolomix/tools/plot_species_kinds.py  group_intensity()
13.│  ├─ metabolomix/tools/plot_species_kinds.py  _common()
14.│  │  ├─ metabolomix/plots/dataset_rows.py  resolve_samples()
15.│  │  ├─ metabolomix/plots/dataset_rows.py  standard_only_spots()
16.│  │  └─ metabolomix/plots/dataset_rows.py  resolve_selection()
17.│  ├─ metabolomix/plots/dataset_rows.py  rows_by_spot()
18.│  └─ metabolomix/plots/group_intensity.py  build_group_intensity_payload()
19.├─ [kind=species] metabolomix/tools/plot_species_kinds.py  species()
20.│  ├─ metabolomix/plots/dataset_rows.py  expand_spots()
21.│  └─ metabolomix/plots/species.py  build_species_payload()
22.├─ [kind=species_pca] metabolomix/tools/plot_species_kinds.py  species_pca()
23.│  ├─ metabolomix/plots/dataset_rows.py  expand_spots()
24.│  ├─ metabolomix/analysis/pca.py  run_pca_fit_subset()
25.│  └─ metabolomix/analysis/result_state.py  register_result()
26.├─ [kind=eic] metabolomix/tools/plot_tools.py  _eic()
27.│  ├─ metabolomix/dataset/artifacts.py  artifact_path()
28.│  ├─ metabolomix/eic/reader.py  read_eic_spot_css1()
29.│  └─ metabolomix/plots/eic.py  build_eic_plot_payload()
30.├─ [kind=eic_compounds] metabolomix/tools/plot_tools.py  _eic_compounds()
31.│  ├─ metabolomix/eic/identity_map.py  select_identity_candidates()
32.│  ├─ metabolomix/eic/reader.py  read_eic_spots_css1()
33.│  └─ metabolomix/plots/eic.py  build_multi_compound_plot_payload()
34.├─ [output=payload] metabolomix/core/serialization.py  json_payload()
35.└─ [output=image] metabolomix/tools/plot_tools.py  _png()
36.   └─ metabolomix/plots/render.py  render_png()

## save_figure

前提: `dataset_plot` で図を作ってある（`ds.figures[kind]`）か、`result_id` で登録済みの結果を名指しする
（pca / volcano / loadings）。無ければ `missing_state`（`required_tools` は `dataset_plot`、データセットが
無ければ `load_dataset`）。結果から作った図は、その結果が現在のものでなければ（`_stale_figure`）保存しない。
状態変更: `reports/figures/` に PNG（`loadings` `group_intensity` `species` `species_pca` は SVG も）を書く。

`kind` は `dataset_plot` の kind（`FIGURE_KINDS`）。`source` 引数と ARF の候補は無い。未知の `kind` と、
`result_id` を使わない `kind`（`group_intensity` `species` `species_pca` `eic` 系）への指定は
`{"status": "error"}`。`pca` / `volcano` は
`save_result_figure()`（探索専用データセットから描いた図には図の中に但し書きが入る）、`loadings` と
`group_intensity` `species` `species_pca` は dpi 300 の PNG と SVG、`eic` 系は dpi 120 の PNG。
`group_intensity` `species` `species_pca` は群を焼き込んだ図なので `metadata_revision` が変わると拒み、
`species_pca` は結果が現在かも確かめる（`_stale_figure`）。戻り値には保存した絶対パスを書く
（Use-LLLM はこれを拾って画像として表示する）。

1. metabolomix/tools/reports.py  save_figure()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ [result_id 指定] metabolomix/tools/plot_tools.py  build_figure()
4. ├─ [result_id 省略] metabolomix/tools/reports.py  _stale_figure()
5. └─ metabolomix/tools/reports.py  _save_record()
6.    ├─ [kind=pca / volcano] metabolomix/plots/result_output.py  save_result_figure()
7.    │  └─ [kind=pca] metabolomix/plots/pca_scores.py  render_pca_scores()
8.    ├─ [kind=loadings] metabolomix/tools/reports.py  _write_png_and_svg()
9.    │  └─ metabolomix/plots/pca_loadings.py  render_loadings_plot()
10.   ├─ [kind=group_intensity / species / species_pca] metabolomix/tools/reports.py  _write_png_and_svg()
11.   │  └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
12.   └─ [kind=eic / eic_compounds] metabolomix/plots/eic.py  render_eic_plot()
