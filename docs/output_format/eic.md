# EIC/AEF (`.EIC.aef` / クロマトグラム)

`.EIC.aef` パーサ、EIC 検索・ランキング、`dataset_plot(kind="eic")` / `dataset_plot(kind="eic_compounds")` の描画契約。

> 先に `lipidmix://docs/output-format`（共通核）を読む。節番号は全トピックの通し番号。

## 7. EIC/AEF (`metabolomix/eic/reader.py`)

### 7.1 `parse_eic_aef_css1()` のスポット出力

型は `list[dict]`。**1要素は1アラインメントスポットのEIC集合**で、`samples` に各測定ファイルのクロマトグラム要約を持つ。

| キー | 型 | 意味 |
|---|---|---|
| `spot_id` | int | ファイル中の0始まりスポット順序。今回のデータではARF/ARF2のID順と対応 |
| `rt` | float | スポット代表RT（min） |
| `ri` | float | スポット代表RI。LCデータでは0のことがある |
| `mz` | float | スポット代表 m/z。**MS-DIAL 2023-11-20 より前の版の出力では全スポット `-1`**（下記） |
| `drift` | float | スポット代表ドリフト時間。未使用時は `-1` |
| `main_type` | int | 横軸種別の列挙値。`0` はRT主軸を表す |
| `num_samples` | int | このスポットに格納されたサンプル数 |
| `samples` | list[dict] | サンプル別EIC要約。リスト長は通常 `num_samples` |

**m/z を欠く古い `.EIC.aef`**: MS-DIAL 2023-11-20 より前の版の LC-MS の出力は全スポットで `mz = -1`（RT は正しい）。EIC 系の処理は読み込み直後に、同じバッチの `<stem>.arf2` の同じ `AlignmentID` の `MassCenter` で補う。RT が ±0.02 min で一致するときだけ補い、補ったスポットには `mz_source = "arf2"` が付き、補った旨が注記（`== 注意 ==` / `caveats`）に出る。ARF2 が無い・RT が合わないスポットは `-1` のままで、m/z 検索に掛からず、`summarize_eic_data()` の `mz_range` からも除かれる。

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

バイナリにはピーク左端・右端座標もあるが、このパーサーの出力辞書には含めない。

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

`search_eic_spots()` は spot ID・m/z 範囲・RT 範囲の AND 検索で、`(一致スポット, 範囲外の spot_ids)` を返す。境界値を含む。m/z を欠く（`-1` のまま補えなかった）スポットは、m/z 条件があるときは一致しない。一致スポットは spot_id 昇順で、各要素は次の形（強度順の `top_eic_spots_by_max_intensity()` も同じ）:

| キー | 意味 |
|---|---|
| `spot_id` | スポットID（= ARF / ARF2 の `MasterAlignmentID`、アラインメント `.dcl` の番号） |
| `rt` | スポット代表RT |
| `mz` | スポット代表 m/z |
| `num_samples` | サンプル数 |
| `max_intensity` | サンプル間で最大の「クロマトグラム最大強度」 |

### 8.4 EIC を検索するツールは無い（特徴の検索は `dataset_feature_search`）

7.4 節の関数は reader のもので、MCP ツールからは呼ばれない。データセットの特徴を窓で探すのは `dataset_feature_search`（`feature_ids`・`min_mz`/`max_mz`・`min_rt`/`max_rt` の AND）で、並びは `sort="mz"`（既定、m/z 昇順）か `sort="intensity"`（強度の降順）だけ。**この強度は特徴行列の検出セル（gap-fill でない実測）の最大で、EIC のクロマトグラムの最大（7.4 節の `max_intensity`）ではない。** 戻り値の欄の意味は `lipidmix://docs/output-format/mztab`。見つけた特徴のクロマトグラムは `dataset_plot(kind="eic")` で描く（特徴 ID は EIC のスポット番号と同じ空間）。

### 8.5 `dataset_plot(kind="eic")` の描画契約

1 特徴（`feature_id` = `dataset_feature_search` の特徴 ID = `.arf2` の MasterAlignmentID。イオンモビリティでは連番にならない）の試料別クロマトグラム。`.EIC.aef` はデータセットが覚えているパス（`artifact_paths` の `eic`）から読む。`.EIC.aef` の無いデータセット（mzTab-M だけから読んだもの）は `COMPANION_ARTIFACT_MISSING`。

既定（`output="image"`）は PNG と 1 行の説明。`output="payload"`（または env `LIPIDMIX_PLOT_OUTPUT=payload`）はクライアント中立の構造化 JSON `plot_schema="lipidmix.eic.v1"` で、主要フィールドは `axes`、スポットの `rt`/`mz`、各試料の `series[].x` / `series[].y` / `file_id` / `peak_left` / `peak_top` / `peak_right`。`output_path` を渡さない描画はファイルを作らない。図は `ds.figures["eic"]` に置かれる。

### 8.6 `dataset_plot(kind="eic_compounds")` の描画契約

**複数物質 × 1サンプル**のオーバーレイ。既定（`output="image"`）はサーバ側で描いた PNG と、何を描き何を落としたかを書いた1行のキャプションで、ファイルは書かない。以下は `output="payload"` の `plot_schema="lipidmix.eic.multi.v1"` の契約。

- **1呼び出し1サンプル**。`file_id` は必須。複数サンプルはサンプルごとに呼ぶ。1物質×全サンプルは `dataset_plot(kind="eic")`。
- 物質は ARF2 の同定情報（データセットの `feature_metadata` の名前・クラス）から選ぶ。`names` は `Name` への大小無視の部分一致、`ontologies` は `Ontology` への大小無視の完全一致で、両者は和集合。
- 単一の `spot` フィールドは持たない。代わりに `sample`（`file_id` / `sample_name` / `class_id`）を持ち、スポット固有の情報は `series[]` の `spot_id` / `name` / `ontology` / `adduct` / `mz` / `rt` に入る。`series` はRT昇順。
- 各系列の `annotation` は `{text: "<脂質名> / <peak_top>", x, y}`。`x`/`y` はピーク頂点に最も近いデータ点の座標で、`y` は正規化後の値。
- **`selection.dropped[]` を必ず読む。** `spot_id`（= ARF2 `AlignmentID`）で引いた EIC スポットは `rt`/`mz` の一致（±0.02 min / ±0.01 Da）を検証し、通らなかった物質は描かない。除外理由は `rt_mismatch` / `mz_mismatch`（ID対応の崩れ）、`spot_out_of_range`、`file_id_absent`（その試料にトレースが無い）、`below_top_n`（`top_n` 件からあふれた。既定 8 件）の5種。**図に無い＝試料に無い、ではない。** m/z を欠く古い EIC（7.1 節）では照合の前に ARF2 の m/z で補うので m/z 照合は自明に通り、取り違えの検出は RT 照合だけが担う（`caveats` にその旨が出る）。`dropped[]` は除外された物質を毎回すべて挙げる。
- `caveats[]` は `dropped[]` を理由ごとに1行へ要約したもの（例: `"below_top_n で 276 件を除外しました（詳細は selection.dropped）。"`）。`rt_mismatch` / `mz_mismatch` は最大10件まで物質名を挙げ、超過分は件数のみ。網羅的な一覧は `selection.dropped[]`。
- `selection.candidates` はクエリに一致した物質の総数（`max_candidates` による予備選抜の**前**）。`selection.candidates_evaluated` は実際にEICファイルから読み rt/mz 検証まで行った件数（予備選抜後）。一致が300件を超えるとARF2 `HeightAverage` 上位300件に予備選抜され、`candidates > candidates_evaluated` となり `caveats` にも残る。`selection.plotted` は実際に描いた数。

**PNG の保存**（`eic` / `eic_compounds` とも）: `save_figure(kind="eic" | "eic_compounds", analysis_id=..., title=None)` で報告書の `reports/figures/` へ、または `dataset_plot(..., output_path=...)` で任意の場所へ書く。どちらも対話描画とは別の書き込み操作で、`lipidmix.eic.v1` と `lipidmix.eic.multi.v1` の両方に対応する。
