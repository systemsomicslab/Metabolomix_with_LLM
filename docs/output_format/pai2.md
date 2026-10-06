# PAI2 (`.pai2` / 単一測定ファイルのピーク)

`.pai2` パーサと `pai2_parser` の出力定義。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。
>
> 対応する Key 番号表: `docs/schema/ChromatogramPeakFeature.md`（MS-DIAL の `[Key(N)]` から抽出した一次資料）。
> リーダーのインデックス定数を変更するときは必ずこちらを先に確認する。

## 5. PAI2 (`metabolomix/pai2/reader.py`)

### 5.1 `deserialize()` のピーク出力

型は `list[dict]`。**1要素は1測定ファイル内の1検出ピーク**である。

| キー | 型 | 意味 |
|---|---|---|
| `time` | dict | ピーク頂点の座標。通常 `{"rt": value}` |
| `time_left` | dict | ピーク左端の座標 |
| `time_right` | dict | ピーク右端の座標 |
| `peak_height` | float | ピーク頂点強度 |
| `peak_height_left` | float | 左端の強度 |
| `peak_height_right` | float | 右端の強度 |
| `peak_area` | float | ピーク面積 |
| `peak_area_above_baseline` | float | ベースラインより上のピーク面積 |
| `m/z` | float | ピーク m/z |
| `S/N` | float | `PeakShape[1]` 由来の S/N |
| `id` | int | 当該 PAI2 ファイル内のピークID |
| `ion_mode` | `IonMode` Enum | `IonMode.Positive`（値0）、`IonMode.Negative`（値1）、`IonMode.Both`（値2）。JSON化時は文字列/数値への変換が必要 |
| `name` | str | 候補化合物名。空文字は未注釈 |
| `formula` | str | 分子式 |
| `ontology` | str | 化学/脂質クラス |
| `smiles` | str | SMILES |
| `inchikey` | str | InChIKey |
| `adduct` | str | 付加体。例: `[M-H]-` |
| `collision_cross_section` | float | CCS。イオンモビリティ未使用時は 0 のことがある |
| `comment` | str | MS-DIAL由来の注釈コメント |
| `has_msms` | bool | MS2 raw ID が非負、または collision-energy map が非空か。MS/MS取得情報の有無 |
| `ms2_raw_id` | int | 元データの MS2 raw spectrum ID。未設定は `-1` |
| `collision_energies` | list[float] | MS/MS取得時の衝突エネルギー一覧。負イオンモードでは `-42.0` など負値もある |
| `msms_peak_count` | int | PAI2 内 `msms_spectrum` の要素数。通常はスペクトル本体がDCL側にあるため0 |
| `msms_spectrum` | list[list] | PAI2自体に格納された `[fragment_mz, intensity]`。通常は空で、本体はDCL側にある |

`time` 系辞書には、入力に存在する非負値だけが入る。

| サブキー | 意味 |
|---|---|
| `rt` | retention time |
| `ri` | retention index |
| `m/z` | 座標表現内の m/z |
| `dt` | drift time |

### 5.2 フィルタ

`filter_features_by_params()` は元と同じピーク辞書を保持したリストを返す。

| パラメータ | 意味 |
|---|---|
| `min_intensity` / `min_height` | `peak_height` の下限 |
| `min_sn` | `S/N` の下限。S/Nが取得できないピークも除外 |

### 5.3 ピーク属性PCA（撤去済み）

**`perform_pca_summary()` は撤去済みで、コード上に存在しない。** 以前この関数は
`[RT, m/z, Height]` の3変数PCAと、その散布図のPNGバイト列を返していた。撤去理由:

- **行がサンプルではなくピーク**だった。PAI2 は単一測定ファイルなのでサンプル間比較
  （オミクスPCA）は原理的にできず、この3変数PCAは生物学的仮説を検定しない
  （RT×m/z散布図の言い換えに近い）。サンプル間の多変量比較は ARF/ARF2 を使う。
- 図をPNGとして返す設計自体も現方針（構造化データを返しクライアントが描画、PNGは
  ユーザーの明示要求時のみ）に反していた。

`pai2_parser()` が返すのは `summarize_pai2_inventory()` による在庫要約（5.5節）のみ。

### 5.4 `inspect_peak_details()`

成功時は `{"status":"success", "matches":[...], "note":...}` を返す。`peak_name` は大文字小文字を無視した部分一致で、複数件返り得る。

| `matches[]` キー | 意味 |
|---|---|
| `id` | PAI2ピークID |
| `name` | 候補化合物名 |
| `m/z` | ピーク m/z |
| `rt` | ピーク頂点RT |
| `height` | ピーク高さ |
| `area` | ピーク面積 |
| `signal_to_noise` | S/N |
| `has_msms_like_fields` | キー名に `msms`、`fragment`、`spectrum` のいずれかを含むフィールドが存在するか。**実スペクトルが非空かは判定しない** |
| `formula` | 分子式 |
| `adduct` | 付加体 |
| `comment` | コメント |

検索条件なしは `status=error`、該当なしは `status=not_found` となる。

### 5.5 `summarize_pai2_inventory()`

`pai2_parser()` が返す在庫要約。PCA を介さず、単一測定ファイルのピーク分布を素直にまとめる。

| キー | 意味 |
|---|---|
| `total_peaks` | フィルタ後のピーク総数 |
| `annotated_peaks` | 確定注釈とみなせるピーク数（空文字・`Unknown`・`no MS2:`・`low score:` は除外） |
| `annotated_fraction` | `annotated_peaks / total_peaks` |
| `ion_mode_counts` | イオンモード別件数 |
| `mz_range` / `rt_range` | m/z・RT の `{min, max}` |
| `height_summary` | ピーク高さの `{min, median, max}` |
| `sn_summary` | S/N の取得割合と `{min, median, max}` |
| `top_by_height` | 強度上位10ピーク（`id`/`name`/`m/z`/`rt`/`height`/`signal_to_noise`） |
| `note` | PAI2は単一サンプルでオミクスPCA不能・MS/MSは`.dcl`参照、の断り書き |

### 8.3 `pai2_parser()` と関連ツール

`pai2_parser()` はテキスト1件を返す（PCA・画像は返さない）。JSONは `summarize_pai2_inventory()` の在庫要約で、`total_peaks` / `annotated_peaks`（`Unknown`・`no MS2:`・`low score:` を除いた確定注釈数）/ `annotated_fraction` / `ion_mode_counts` / `mz_range` / `rt_range` / `height_summary` / `sn_summary` / `top_by_height`（強度上位10ピーク）を含む。PAI2 は単一測定ファイルのピーク一覧なので、サンプル間比較（オミクスPCA）はこの単位では行わない（複数サンプルの多変量比較は ARF/ARF2）。旧 `pai2_get_top_metabolites()` / `pai2_update_analysis_filter()`（いずれもピーク属性PCA依存）は撤去した。`pai2_inspect_peak(peak_id=None, peak_name=None)` は5.4節の辞書をJSON文字列で返す（旧名 `pai2_inspect_metabolite_details` から peak 語彙へ統一）。MS/MS は同名 `.dcl`（`dcl_index` がリスト順に対応）を参照する。
