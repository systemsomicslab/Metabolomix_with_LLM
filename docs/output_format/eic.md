# EIC/AEF (`.EIC.aef` / クロマトグラム)

`.EIC.aef` パーサ、EIC 検索・ランキング、`eic_plot_chromatograms` / `eic_plot_compounds` の描画契約。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

## 7. EIC/AEF (`metabolomix/eic/reader.py`)

### 7.1 `parse_eic_aef_css1()` のスポット出力

型は `list[dict]`。**1要素は1アラインメントスポットのEIC集合**で、`samples` に各測定ファイルのクロマトグラム要約を持つ。

| キー | 型 | 意味 |
|---|---|---|
| `spot_id` | int | ファイル中の0始まりスポット順序。今回のデータではARF/ARF2のID順と対応 |
| `rt` | float | スポット代表RT（min） |
| `ri` | float | スポット代表RI。LCデータでは0のことがある |
| `mz` | float | スポット代表 m/z |
| `drift` | float | スポット代表ドリフト時間。未使用時は `-1` |
| `main_type` | int | 横軸種別の列挙値。`0` はRT主軸を表す |
| `num_samples` | int | このスポットに格納されたサンプル数 |
| `samples` | list[dict] | サンプル別EIC要約。リスト長は通常 `num_samples` |

### 7.2 `samples[]` の全フィールド

**1要素は「1スポット × 1測定ファイル」のEIC**である。

| キー | 型 | 意味 |
|---|---|---|
| `file_id` | int | データセット内の測定ファイルID |
| `peak_top` | float | ピーク頂点の**横軸座標**。今回の `main_type=0` データではRT（min）。**強度ではない** |
| `num_peaks` | int | EICに格納されたクロマトグラム点数。名称はピーク数に見えるがデータ点数 |
| `mean_intensity` | float | 全クロマトグラム点の強度算術平均 |
| `max_intensity` | float | 全クロマトグラム点の最大強度 |
| `chromatogram` | list[tuple] | `include_chromatogram=true` のときだけ追加される `(horizontal_coordinate, intensity)` 点列。長さは `num_peaks` |
注意: バイナリにはピーク左端・右端座標もあるが、現パーサーは読み進めるだけで出力辞書には含めない。

### 7.3 `summarize_eic_data()`

| キー | 意味 |
|---|---|
| `total_spots` | EICスポット総数 |
| `rt_range` | スポット代表RTの `(min, max)` |
| `mz_range` | スポット代表 m/z の `(min, max)` |
| `total_samples` | 全スポットの `num_samples` 合計。ユニークサンプル数ではない |
| `total_peaks` | 全 `samples[].num_peaks` の合計。実際には全クロマトグラム点数 |
| `peak_top_mean` | 全サンプルのピーク頂点横軸座標の平均。強度平均ではない |
| `peak_top_max` | 全サンプルのピーク頂点横軸座標の最大。最大強度ではない |
| `unique_file_ids` | 出現した `file_id` の昇順リスト。長さがユニークサンプル数 |

### 7.4 検索・ランキング出力

`search_eic_by_mz_range()` と `search_eic_by_rt_range()` は条件に合う**元のスポット辞書全体**を返す。境界値を含む。

`top_eic_spots_by_max_intensity()` の1行:

| キー | 意味 |
|---|---|
| `spot_id` | スポットID |
| `rt` | スポット代表RT |
| `mz` | スポット代表 m/z |
| `num_samples` | サンプル数 |
| `max_intensity` | サンプル間で最大の「クロマトグラム最大強度」 |

MCPツール `eic_rank_by_max_intensity` はこの `max_intensity` 降順で並べる**強度上位ランキング**である。旧 `top_eic_spots_by_peak_top` / `eicaef_top_peak_tops` は `peak_top`（ピーク頂点の横軸=RT座標。強度ではない）で並べており、実質「遅いRTのスポット一覧」で強度上位ではなかったため置き換えた。

### 8.4 `eic_parser()` と関連ツール

`eic_parser()` は7.3節の要約をJSON文字列として返す。m/z/RT検索の各表示行は `spot_id`, `rt`, `mz`, `num_samples` を持つ。`eic_rank_by_max_intensity()` はさらに `max_intensity`（サンプル間で最大のクロマトグラム最大強度）を表示する強度上位ランキングである（7.4節）。

### 8.5 `eic_plot_chromatograms()` の描画契約

`eic_plot_chromatograms()` は画像ではなく、クライアント中立の構造化JSON
`plot_schema="lipidmix.eic.v1"` を返す。主要フィールドは `axes`、スポットの
`rt`/`mz`、および各試料の `series[].x` / `series[].y` / `file_id` /
`peak_left` / `peak_top` / `peak_right` である。Use-LLLMはこれをPlotlyへ変換し、
Claude Desktop等は各クライアントのUI方式で描画できる。通常の描画ではファイルを
作らない。

PNGが必要だとユーザーが明示した場合に限り、先に得たプロット情報を
`save_eic_figure(analysis_id, title=None)` で `reports/figures/` へ保存する。
この保存は対話描画とは別の書き込み操作である。

### 8.6 `eic_plot_compounds()` の描画契約

`eic_plot_compounds()` は**複数物質 × 1サンプル**のオーバーレイを返す。**既定
（`output="image"`）ではサーバ側で描いた PNG と、何を描き何を落としたかを書いた1行の
キャプション**を返し、ファイルは書かない。以下の構造化JSON
`plot_schema="lipidmix.eic.multi.v1"` は `output="payload"`（または env
`LIPIDMIX_PLOT_OUTPUT=payload`）でクライアントが自分で描くときの契約。

- **1呼び出し1サンプル**。`file_id` は必須。複数サンプルを比べるときはサンプルごとに
  呼び出し、図を並べる。1物質×全サンプルの比較は従来どおり `eic_plot_chromatograms`。
- 物質は ARF2 の同定情報から選ぶ。`names` は `Name` への大小無視の部分一致、
  `ontologies` は `Ontology` への大小無視の完全一致で、両者は和集合。
- 単一の `spot` フィールドは持たない。代わりに `sample`（`file_id` / `sample_name` /
  `class_id`）を持ち、スポット固有の情報は `series[]` の `spot_id` / `name` /
  `ontology` / `adduct` / `mz` / `rt` に入る。`series` はRT昇順。
- 各系列の `annotation` は `{text: "<脂質名> / <peak_top>", x, y}`。`x`/`y` はピーク頂点に
  最も近いデータ点の座標で、`y` は正規化後の値である。
- **`selection.dropped[]` を必ず読むこと。** `spot_id`（= ARF2 `AlignmentID`）で引いた
  EICスポットは `rt`/`mz` の一致（±0.02 min / ±0.01 Da）を検証しており、通らなかった
  物質は描画されない。除外理由は `rt_mismatch` / `mz_mismatch`（ID対応の崩れ）、
  `spot_out_of_range`、`file_id_absent`（その試料にトレースが無い）、
  `below_top_n`（`top_n` 件からあふれた。既定 8 件）の5種。**図に無い＝試料に無い、ではない。**
  `dropped[]` は毎回、除外された物質を1件残らず列挙する。
- `caveats[]` は `dropped[]` を理由ごとに1行へ要約したもの（例:
  `"below_top_n で 276 件を除外しました（詳細は selection.dropped）。"`）。ただし
  `rt_mismatch` / `mz_mismatch` は個別診断に有用なため、最大10件まで物質名を列挙し
  超過分は件数のみ添える。網羅的な一覧が要るときは常に `selection.dropped[]` を見ること。
- `selection.candidates` はクエリに一致した物質の総数（`max_candidates` による予備選抜が
  行われる**前**の件数）。`selection.candidates_evaluated` は実際にEICファイルから読み出し
  rt/mz 検証まで行った件数（予備選抜後）。一致が300件を超えるとARF2 `HeightAverage`
  上位300件に予備選抜され、`candidates > candidates_evaluated` となり `caveats` にも残る。
  `selection.plotted` は実際に描画した数。

PNGが必要だとユーザーが明示した場合に限り `save_eic_figure(analysis_id, title=None)` で
保存する。この保存ツールは `lipidmix.eic.v1` と `lipidmix.eic.multi.v1` の両方に対応する。
