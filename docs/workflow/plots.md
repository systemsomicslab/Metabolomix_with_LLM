# ワークフロー: 図（dataset_plot）・図の保存と payload 契約

データセットの図は `dataset_plot` が描く（既定は画像。ファイルは書かない）。作った図の記録は
`ds.figures[kind]`（kind ごとに最後の 1 つ。指紋に含めない）に置かれ、ユーザーが明示的に保存を
求めたときだけ `save_figure` がそこから `reports/figures/<analysis_id>_<kind>.png` に書く。
旧 ARF 系の描画ツール（`arf_plot_*` `plot_pca_loadings` `eic_plot_*`）は B-2b Task 8 まで残るが、
`save_figure` はもう読まない（画像返しだけ）。

## payload 契約の比較

| payload | 生成元 | 保存 | 特徴 |
|---|---|---|---|
| `lipidmix.volcano.v1` | `dataset_plot(kind="volcano", output="payload")` | `save_figure(kind="volcano")` | `up`/`down` は全点保持、`ns` のみ間引く。保存は間引き前の全点 |
| `lipidmix.pca_scores.v1` | `dataset_plot(kind="pca", output="payload")` | `save_figure(kind="pca")` | 散布図。点の `group` が `color_by`（group / batch / role）の値 |
| `lipidmix.pca_loadings.v1` | `dataset_plot(kind="loadings", output="payload")` | `save_figure(kind="loadings")` | 主成分ごとの横棒。値は r か coefficient。色は `ontology`。PNG（dpi 300）と SVG |
| `lipidmix.eic.v1` | `dataset_plot(kind="eic", output="payload")` | `save_figure(kind="eic")` | 1 特徴 × 試料別の線グラフ |
| `lipidmix.eic.multi.v1` | `dataset_plot(kind="eic_compounds", output="payload")` | `save_figure(kind="eic_compounds")` | 1 試料 × 複数物質の重ね描き |
| `lipidmix.group_intensity.v1` / `lipidmix.species_intensity.v1` / `lipidmix.pca_loadings.v1`（旧） | `arf_plot_group_intensity` / `arf_plot_species` / `plot_pca_loadings` | （Task 7 で `dataset_plot` の kind になるまで保存不可） | 旧 ARF 系 |

描画ツールの**既定は payload ではなく画像**（`output="image"`）。座標点列を LLM の文脈へ流すと
実測で数万トークンかかるのに対し、サーバ側で描いた PNG は 300〜1,000 画像トークンに収まる。
Plotly で描くクライアント（Use-LLLM）は起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、
呼び出しごとに `output="payload"` を渡す。判定は `plots/render.py resolve_plot_output()`。
描画は `plots/` 側に置き、画像返しと保存が同じ関数を通る（画面の図とレポートの図がずれない）。

保存先は `_resolve_report_dir()` が解決する。`lipidmix.local.toml` の `[paths] reports_dir` の明示先を優先し、書き込めなければ
解析フォルダ配下の `reports/` に退避する。未指定時は解析フォルダの `reports/` だけを使い、
リポジトリ側へは退避しない（書けなければエラー）。

## dataset_plot

前提: `kind` ごと——`volcano` は `dataset_differential`（2 群）、`pca` / `loadings` は `dataset_pca`、
`eic` / `eic_compounds` はデータセットの `.EIC.aef`（`artifact_paths` の `eic`）。無ければ
`missing_state`（`dataset_pca` / `dataset_differential`）か `COMPANION_ARTIFACT_MISSING`（`.EIC.aef` が無い）。
状態変更: `ds.figures[kind]` を更新（メモリ上。ファイルは書かない）。

結果の選び方（`_pick_result`）: `result_id` があれば `ds.results` から引き、現在のものでなければ
`STALE_ANALYSIS_RESULT`（`result_state.assert_current`）。無ければ直近（`ds.last_pca` /
`ds.last_differential`）。`dataset_exclude` や前処理・群の変更で古くなった結果は
`result_state.invalidate_results` が `last_*` と `ds.figures` の該当 kind を同時に捨てるので、
古い図は最新として描かれず保存もされない。別の kind の引数が渡されたら `DATASET_BAD_REQUEST`
（`_check_arguments`。どの kind の引数かを書く）。

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
12.├─ [kind=eic] metabolomix/tools/plot_tools.py  _eic()
13.│  ├─ metabolomix/dataset/artifacts.py  artifact_path()
14.│  ├─ metabolomix/eic/reader.py  read_eic_spot_css1()
15.│  └─ metabolomix/plots/eic.py  build_eic_plot_payload()
16.├─ [kind=eic_compounds] metabolomix/tools/plot_tools.py  _eic_compounds()
17.│  ├─ metabolomix/eic/identity_map.py  select_identity_candidates()
18.│  ├─ metabolomix/eic/reader.py  read_eic_spots_css1()
19.│  └─ metabolomix/plots/eic.py  build_multi_compound_plot_payload()
20.├─ [output=payload] metabolomix/core/serialization.py  json_payload()
21.└─ [output=image] metabolomix/tools/plot_tools.py  _png()
22.   └─ metabolomix/plots/render.py  render_png()

## save_figure

前提: `dataset_plot` で図を作ってある（`ds.figures[kind]`）か、`result_id` で登録済みの結果を名指しする
（pca / volcano / loadings）。無ければ `missing_state`（`required_tools` は `dataset_plot`、データセットが
無ければ `load_dataset`）。結果から作った図は、その結果が現在のものでなければ（`_stale_figure`）保存しない。
状態変更: `reports/figures/` に PNG（`loadings` は SVG も）を書く。

`kind` は `dataset_plot` の kind（`FIGURE_KINDS`）。`source` 引数と ARF の候補は無い。未知の `kind` と、
`result_id` を使わない `kind` への指定は `{"status": "error"}`。`pca` / `volcano` は
`save_result_figure()`（探索専用データセットから描いた図には図の中に但し書きが入る）、`loadings` は
dpi 300 の PNG と SVG、`eic` 系は dpi 120 の PNG。戻り値には保存した絶対パスを書く
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
10.   └─ [kind=eic / eic_compounds] metabolomix/plots/eic.py  render_eic_plot()

## plot_pca_loadings

前提: ローディングを持つ PCA 結果がある（`arf_parser` / `arf_pca_preprocessed` / `arf_pca_species` /
`dataset_pca`。無ければ `MissingState`）
状態変更: `session.last_loadings_plot`（payload・title）を更新。ファイルは書かない。

入力元は 3 つ（手順 2）。選び方は `save_figure(kind="pca")` と同じで、どれも優先せず曖昧なら
`AMBIGUOUS_RESULT_SOURCE` で止まる（手順 3）。ARF は上位 N の選び方とスポット情報に
`get_pca_loading_features()` を使い（手順 4〜5）、分子種 PCA は保存済みの行を、mzTab は
`pp_feature_names` と特徴量の注釈を使う（手順 6）。r は autoscale した PCA でだけ
成分 × 特異値 / √n で出し、出せなければ coefficient で描く（手順 7）。

1. metabolomix/tools/pca_loadings_tools.py  plot_pca_loadings()
2. └─ metabolomix/tools/pca_loadings_tools.py  _candidates()
3. └─ metabolomix/plots/result_output.py  select_result()
4. ├─ [source=arf] metabolomix/tools/pca_loadings_tools.py  _arf_features()
5. │  └─ metabolomix/arf/reader.py  get_pca_loading_features()
6. ├─ [source=mztab] metabolomix/tools/pca_loadings_tools.py  _mztab_features()
7. └─ metabolomix/plots/pca_loadings.py  build_loadings_payload()
8. ├─ [output=payload] metabolomix/core/serialization.py  json_payload()
9. └─ [output=image] metabolomix/plots/pca_loadings.py  render_loadings_plot()
10.   └─ metabolomix/plots/render.py  figure_to_png()
