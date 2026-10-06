# ワークフロー: 図の保存と payload 契約

通常の対話では PNG を作らない。ユーザーが明示的に保存を求めたときだけ `save_*_figure`
を呼ぶ。4 ツールは同じ骨格を持つ ——「`session` に記録済みの payload を取り出し、
ファイル名の slug を作り、レポート先ディレクトリを解決し、matplotlib で描いて保存する」。

## payload 契約の比較

| payload | 生成元 | 保存 | 特徴 |
|---|---|---|---|
| PCA スコア（`session.arf.last_pca_plot` / `session.dataset.last_pca`） | `arf_parser` / `arf_pca_preprocessed` / `dataset_pca` | `save_pca_figure` | 散布図。`_pca_scatter_arrays` で座標配列に展開 |
| `lipidmix.volcano.v1` | `arf_plot_volcano(output="payload")` | `save_volcano_figure` | `up`/`down` は全点保持、`ns` のみ間引く |
| `lipidmix.group_intensity.v1` | `arf_plot_group_intensity(output="payload")` | `save_group_intensity_figure` | 項目 × 群 × 試料の PeakHeight 合計。点も平均 ± SD も全量保持（間引かない）。PNG に加えて SVG も書く。描画は `plots/group_intensity.py` に委譲 |
| `lipidmix.eic.v1` / `.multi.v1` | `eic_plot_chromatograms` / `eic_plot_compounds(output="payload")` | `save_eic_figure` | 線グラフ。描画は `plots/eic.py` に委譲 |

`arf_plot_volcano` と `eic_plot_compounds` の**既定は payload ではなく画像**
（`output="image"`）。座標点列を LLM の文脈へ流すと実測で数万トークンかかるのに対し、
サーバ側で描いた PNG は画素課金で 300〜1,000 画像トークンに収まる。Plotly で対話的に
描くクライアント（Use-LLLM）は起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、
呼び出しごとに `output="payload"` を渡す。判定は `plots/render.py resolve_plot_output()`。

描画はいずれも `plots/` 側に置いてある（`render_eic_plot()` / `render_volcano_plot()`）。
画像返しと PNG 保存が同じ関数を通るので、画面の図とレポートに貼る図がずれない。
`save_pca_figure` だけは `tools/reports.py` 内で直接 matplotlib を呼ぶ。

保存先は `_resolve_report_dir()` が解決する。`LIPIDMIX_REPORTS_DIR` の明示先を優先し、
書き込めなければ解析フォルダ配下の `reports/` に退避する。未指定時は解析フォルダの
`reports/` を優先し、不可ならリポジトリ内の `reports/` に退避する
（候補列挙と可否判定が分かれているのは、退避したこと自体を返り値で開示するため）。

## save_pca_figure

前提: `arf_parser` / `arf_pca_preprocessed` / `load_dataset` / `dataset_pca` の
いずれか実行済み（未実行なら手順 3 で `MissingState`）
状態変更: PNG ファイルを書き出す。

入力元は 2 つある。ARF 経路は `session.arf.last_pca_plot`、mzTab-M 経路は
`session.dataset.last_pca`。手順 2 が候補を並べ、手順 3 が 1 件に決める。
**どちらも優先しない**——有効な結果が 2 つ以上あれば `AMBIGUOUS_RESULT_SOURCE` で
止まり、`source`（`arf` / `mztab`）か `result_id` の指定を求める。どちらを描くかは
図の数字そのものを変えるので、優先順位を決め打ちすると「なぜこの図なのか」を
説明できなくなる。前処理をやり直して古くなった結果は候補にならない（`valid=False`）。
生行列 PCA と前処理後 PCA は ARF 側の同じスロットを使うので、ARF 経路で保存されるのは
**直近に実行したほう**の図になる。

1. metabolomix/tools/reports.py  save_pca_figure()
2. └─ metabolomix/tools/reports.py  _pca_candidates()
3. │  └─ metabolomix/core/tool_helpers.py  dataset_pca_plot()
4. │  └─ metabolomix/analysis/result_state.py  is_current()
5. └─ metabolomix/tools/reports.py  _choose_figure_result()
6. │  └─ metabolomix/plots/result_output.py  select_result()
7. └─ metabolomix/core/mcp_errors.py  missing_state()
8. └─ metabolomix/tools/reports.py  _save_figure()
9. │  └─ metabolomix/corpus/knowledge_store.py  make_slug()
10.│  └─ metabolomix/core/mcp_core.py  _resolve_report_dir()
11.│  │  └─ metabolomix/core/mcp_core.py  _report_dir_candidates()
12.│  │  └─ metabolomix/core/mcp_core.py  _first_writable_dir()
13.│  └─ metabolomix/plots/result_output.py  save_result_figure()
14.│     └─ metabolomix/plots/result_output.py  figure_annotations()
15.│     └─ metabolomix/core/tool_helpers.py  _pca_scatter_arrays()

## save_volcano_figure

前提: `arf_differential` または `dataset_differential` 実行済み
（未実行なら手順 3 で `MissingState`）
状態変更: PNG ファイルを書き出す。

`arf_plot_volcano` を経由する必要はない。`session.arf.last_differential` または
`session.dataset.last_differential` に保持された全量 volcano を直接描く
（`arf_plot_volcano` の画像モードと同じ描画関数）。両経路が
`differential.volcano_data()` の点列を共有しているので、描画関数は 1 つで足りる。
入力元の決め方は save_pca_figure と同じ（優先せず、曖昧なら止める）。戻り値には
選ばれた `source=` と `result_id=` を書く。

1. metabolomix/tools/reports.py  save_volcano_figure()
2. └─ metabolomix/tools/reports.py  _differential_candidates()
3. │  └─ metabolomix/analysis/result_state.py  is_current()
4. └─ metabolomix/tools/reports.py  _choose_figure_result()
5. │  └─ metabolomix/plots/result_output.py  select_result()
6. └─ metabolomix/core/mcp_errors.py  missing_state()
7. └─ metabolomix/tools/reports.py  _save_figure()
8. │  └─ metabolomix/corpus/knowledge_store.py  make_slug()
9. │  └─ metabolomix/core/mcp_core.py  _resolve_report_dir()
10.│  │  └─ metabolomix/core/mcp_core.py  _report_dir_candidates()
11.│  │  └─ metabolomix/core/mcp_core.py  _first_writable_dir()
12.│  └─ metabolomix/plots/result_output.py  save_result_figure()
13.│     └─ metabolomix/plots/volcano.py  render_volcano_plot()

## save_eic_figure

前提: `eic_plot_chromatograms` または `eic_plot_compounds` 実行済み
（未実行なら手順 2 で `MissingState`）
状態変更: PNG ファイルを書き出す。

手順 7 は payload の `plot_schema` で描画関数を選ぶ。未知のスキーマは `ValueError`。

1. metabolomix/tools/reports.py  save_eic_figure()
2. └─ metabolomix/core/mcp_errors.py  missing_state()
3. └─ metabolomix/corpus/knowledge_store.py  make_slug()
4. └─ metabolomix/core/mcp_core.py  _resolve_report_dir()
5. │  └─ metabolomix/core/mcp_core.py  _report_dir_candidates()
6. │  └─ metabolomix/core/mcp_core.py  _first_writable_dir()
7. └─ metabolomix/plots/eic.py  render_eic_plot()
8.    ├─ [lipidmix.eic.v1] metabolomix/plots/eic.py  _render_single_spot()
9.    └─ [lipidmix.eic.multi.v1] metabolomix/plots/eic.py  _render_multi_compound()

## save_group_intensity_figure

前提: `arf_plot_group_intensity` 実行済み（未実行なら `MissingState`。`required_tools` は
`arf_plot_group_intensity`）
状態変更: `reports/figures/` に PNG（dpi 300）と同名の SVG を書く。

`session.arf.last_group_intensity` の payload・title・ncols から描き直す。`title` を省くと
`arf_plot_group_intensity` に渡したものを使う。

1. metabolomix/tools/reports.py  save_group_intensity_figure()
2. └─ metabolomix/core/mcp_errors.py  missing_state()
3. └─ metabolomix/corpus/knowledge_store.py  make_slug()
4. └─ metabolomix/core/mcp_core.py  _resolve_report_dir()
5. └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
