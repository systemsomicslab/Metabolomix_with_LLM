# ワークフロー: 図の保存と payload 契約

通常の対話では PNG ファイルを作らない。ユーザーが明示的に保存を求めたときだけ `save_figure`
を呼ぶ。`kind` で「どの直近結果を読むか」を振り分け、`session` に記録済みの結果から描き直して
`reports/figures/<analysis_id>_<kind>.png` に書く。

## payload 契約の比較

| payload | 生成元 | 保存 | 特徴 |
|---|---|---|---|
| PCA スコア（`session.arf.last_pca_plot` / `session.dataset.last_pca`） | `arf_parser` / `arf_pca_preprocessed` / `arf_pca_species`（`session.arf.last_species_pca`）/ `dataset_pca` | `save_figure(kind="pca")` | 散布図 |
| `lipidmix.volcano.v1` | `arf_plot_volcano(output="payload")` | `save_figure(kind="volcano")` | `up`/`down` は全点保持、`ns` のみ間引く。保存は間引き前の全点 |
| `lipidmix.group_intensity.v1` | `arf_plot_group_intensity(output="payload")` | `save_figure(kind="group_intensity")` | 項目 × 群 × 試料の PeakHeight 合計。PNG（dpi 300）に加えて SVG も書く |
| `lipidmix.species_intensity.v1` | `arf_plot_species(output="payload")` | `save_figure(kind="species")` | スポット × 群 × 試料の高さと割合。PNG（dpi 300）と SVG |
| `lipidmix.eic.v1` / `.multi.v1` | `eic_plot_chromatograms` / `eic_plot_compounds(output="payload")` | `save_figure(kind="eic")` | 線グラフ |

描画ツールの**既定は payload ではなく画像**（`output="image"`）。座標点列を LLM の文脈へ流すと
実測で数万トークンかかるのに対し、サーバ側で描いた PNG は 300〜1,000 画像トークンに収まる。
Plotly で描くクライアント（Use-LLLM）は起動 env に `LIPIDMIX_PLOT_OUTPUT=payload` を置くか、
呼び出しごとに `output="payload"` を渡す。判定は `plots/render.py resolve_plot_output()`。
描画は `plots/` 側に置き、画像返しと保存が同じ関数を通る（画面の図とレポートの図がずれない）。

保存先は `_resolve_report_dir()` が解決する。`LIPIDMIX_REPORTS_DIR` の明示先を優先し、書き込めなければ
解析フォルダ配下の `reports/` に退避する。未指定時は解析フォルダの `reports/` を優先し、不可なら
リポジトリ内の `reports/` に退避する。

## save_figure

前提: `kind` に対応する描画・解析を実行済み（無ければ `MissingState`。`required_tools` は kind ごと）
状態変更: `reports/figures/` に PNG（`group_intensity` と `species` は SVG も）を書く。

`kind` が `pca` / `volcano` のときだけ入力元を選ぶ（手順 2〜4）。ARF 経路と mzTab-M 経路の結果を
並べ、**どちらも優先しない**——有効な結果が 2 つ以上あれば `AMBIGUOUS_RESULT_SOURCE` で止まり、
`source` か `result_id` の指定を求める。前処理をやり直して古くなった結果は候補にならない。
未知の `kind` と、`source` / `result_id` を使わない `kind` への指定は `{"status": "error"}`。
戻り値には保存した絶対パスを書く（Use-LLLM はこれを拾って画像として表示する）。

1. metabolomix/tools/reports.py  save_figure()
2. ├─ [kind=pca] metabolomix/tools/reports.py  _pca_candidates()
3. ├─ [kind=volcano] metabolomix/tools/reports.py  _differential_candidates()
4. ├─ [kind=pca / volcano] metabolomix/tools/reports.py  _choose_figure_result()
5. │  └─ metabolomix/plots/result_output.py  select_result()
6. │  └─ metabolomix/tools/reports.py  _save_figure()
7. │     └─ metabolomix/plots/result_output.py  save_result_figure()
   │        └─ [kind=pca] metabolomix/plots/pca_scores.py  render_pca_scores()
8. ├─ [kind=eic] metabolomix/tools/reports.py  _save_eic()
9. │  └─ metabolomix/plots/eic.py  render_eic_plot()
10.├─ [kind=group_intensity] metabolomix/tools/reports.py  _save_group_intensity()
11.│  └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
12.└─ [kind=species] metabolomix/tools/reports.py  _save_species()
13.   └─ metabolomix/plots/species.py  render_species_plot()
