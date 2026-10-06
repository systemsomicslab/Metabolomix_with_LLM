# ARF2 (`.arf2` / スポット代表)

`.arf2` パーサと `arf2_parser` の出力定義。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。
>
> 対応する Key 番号表: `docs/schema/AlignmentSpotProperty.md`（MS-DIAL の `[Key(N)]` から抽出した一次資料）。
> リーダーのインデックス定数を変更するときは必ずこちらを先に確認する。

## 4. ARF2 (`metabolomix/arf2/reader.py`)

### 4.1 `deserialize()` / `extract_arf2_data()`

型は `list[dict]`。**1行/1要素は1アラインメントスポットのカタログ情報**で、サンプル別強度は含まない。

| 列/キー | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | MS-DIAL のマスターアラインメントID |
| `AlignmentID` | int | アラインメントID |
| `RT` | float | スポット中心RT（min） |
| `MassCenter` | float | スポット中心 m/z |
| `IonMode` | str | `Positive`、`Negative`、`Both`、`Unknown` |
| `Name` | str | 候補化合物名。欠落時は `Unknown`。`no MS2:` や `low score:` は注釈の確度に関する接頭辞 |
| `HeightAverage` | float | アラインメントスポットのサンプル間平均ピーク高さ |
| `Formula` | str | 分子式。例: `C5H10O2` |
| `Ontology` | str | 化学/脂質クラス。例: `FA`, `PC`, `Cer_NS` |
| `SMILES` | str | 分子構造の SMILES 表現 |
| `InChIKey` | str | 構造同定子 InChIKey |
| `AdductType` | str | 観測イオンの付加体。例: `[M-H]-` |
| `HeightMin` | float | サンプル間ピーク高さの最小値 |
| `HeightMax` | float | サンプル間ピーク高さの最大値 |
| `PeakWidthAverage` | float | サンプル間の平均ピーク幅。LC の場合は通常 min |
| `SignalToNoiseAve` | float | サンプル間 S/N の平均 |
| `SignalToNoiseMax` | float | サンプル間 S/N の最大値 |
| `SignalToNoiseMin` | float | サンプル間 S/N の最小値 |
| `MassMin` | float | サンプル間観測 m/z の最小値 |
| `MassMax` | float | サンプル間観測 m/z の最大値 |
| `FillPercentage` | float | 当該スポットが値を持つサンプルの割合。**0..1 の比率**で、百分率表示には100倍する |
| `MonoIsotopicPercentage` | float | 単同位体ピークとして扱われる割合。**0..1 の比率** |

### 4.2 `summarize_arf2_data()`

| キー | 意味 |
|---|---|
| `total_spots` | スポット総数 |
| `height_average_median` | 0より大きい `HeightAverage` の中央値 |
| `height_average_max` | 0より大きい `HeightAverage` の最大値 |
| `rt_range` | RT の `(min, max)` |
| `mass_range` | m/z の `(min, max)` |
| `ion_modes` | イオンモード別件数 |
| `annotated_count` | `Name` が空、null、`Unknown` ではない件数。注釈確度は考慮しない |
| `annotation_rate` | `annotated_count / total_spots * 100` |
| `ontology_top` | 空でない `Ontology` の件数上位10クラス。文字列 `Unknown` も集計対象 |
| `sn_median` | 0より大きい `SignalToNoiseAve` の中央値 |
| `error` | データが空のときのエラーメッセージ |

### 4.3 `format_spots_as_table()`

1行目が列名、2行目以降が1スポットの TSV/CSV 文字列。既定列順は上記22列の順である。浮動小数は最大4桁程度に丸め、区切り文字・改行・引用符を含むセルは二重引用符で囲む。空値は空セルになる。この丸め済み表はLLM受け渡し用であり、精密な再計算には元の辞書値を使う。

### 8.2 `arf2_parser()`

返り値はMarkdownテキスト1件（`str`）。`total_spots`, `annotated_count`, `annotation_rate`, RT/m/z範囲、強度中央値、イオンモード、S/N中央値、Ontology上位を自然言語で表示する。個々の22列は返さず、セッション内部に保持する。
