# データセット（MS-DIAL 解析フォルダから読んだ `DatasetState`）

`load_dataset(directory=...)` で MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.mddata`）を読んだデータセットの出力定義。
`.arf` / `.arf2` / `.pai2` のパーサ出力も含む。

> 先に `lipidmix://docs/output-format`（共通核）を読む。脂質名の文法（§2.1。脂質アッセイのときだけ当てる）と必須注意事項（§9）はそちらで定義する。

## 概要

- **行の粒度は特徴×試料**。特徴＝アラインメントスポット（`MasterAlignmentID`。`.arf2` の `SMF_ID` と同じ軸）。名前・m/z・RT・
  同定は特徴ごとに 1 つ（全試料共通の代表値）で、試料ごとの値は強度の行列（`feature_matrix`。値は `measure` で選んだ `.arf` の欄）だけが持つ。
- **gap-fill と検出の印**。`detected_mask`（特徴 × 試料の bool）が `False` のセルは検出ピークではなく MS-DIAL が埋めた補間値
  （`.arf` の `IsGapFilled` / 負の `MasterPeakID`）。検出数・存在判定・検出率は `detected_mask` で数え、非ゼロで数えない。
  解析フォルダのデータセットでは常に付く。
- **試料ごとのピークの RT**（`per_sample_rt`。特徴 × 試料）。`dataset_sample_peaks(feature_ids)` が `rt` `mz` `height` `detected` の
  行で返す（gap-fill の `rt` は検出の位置ではない）。mzTab-M で読んだデータセットは隣の `.arf` と照合できたときだけ持つ。
  意味は `reanalysis` トピック §6.1。
- **`.arf2` は特徴の代表**（スポット代表カタログ）。`Name` `Ontology` `Formula` `InChIKey` `AdductType` と試料間の統計
  （`HeightAverage` `FillPercentage` `SignalToNoise*` など）を持ち、試料別の強度は持たない。`.arf` の `Name` と食い違えば
  `.arf2` を優先する（`feature_metadata`）。
- **SME の理論 m/z（`directory=` で開いたとき）**。`.arf2` は照合した参照の記録のプリカーサー m/z（mzTab-M の SME の
  `theoretical_mass_to_charge`）を持たないので、SME（rank 1）の `theoretical_mass_to_charge` を `.arf2` の `Formula` と
  `AdductType` のモノアイソトピック m/z として計算する（`dataset_feature_search` の目標の理論 m/z と同じ計算）。
  **注釈が表す m/z で、ライブラリの記録の値ではない**（記録が丸め値・実測値なら同じアラインメントの mzTab-M の値と違う）。
  組成式か付加体が無い・読めない（H C N O P S Na Cl K F Br と D 以外の元素を含む）同定は空のままで、質量精度と新しい知見の
  候補に入らない。件数は `load_notes`、出所は `dataset_mass_accuracy` / `dataset_novel_candidates` の `theoretical_mz_source`
  （`formula_adduct`）。根拠バンドルの `msdial_sme` 行の来歴には写さない（`provenance.theoretical_mass_to_charge` はライブラリの
  記録の値の意味で、この入口では `null`）。`job_path` / `pipeline_path` で解析フォルダを読んだときは、照合した正準 mzTab-M の
  SME の値を同じスポット ID に写す（`mztab`。計算値で上書きしない）。意味は `reanalysis` トピック §8。
- **同定の標準化の列**。`.arf2` の注釈は GOSLIN で正規化し、`feature_metadata["identity"]` に `normalized`（標準化した脂質名）・
  `msi_level`（MSI の同定水準。**`.arf2` 由来のヒューリスティックで、MS/MS の有無ではない**）が入る（`dataset_feature_search` の
  行にも出る。意味は `identity` トピック）。脂質アッセイでないデータセットでは当たらない。
- **試料の役割と群**。試料は役割（`sample` / `qc` / `blank`）・群（Class ID）・バッチ・注入順を持つ（`dataset_status` の `samples`）。
  ブランクは前処理で行列から外れ、QC は残る。QC・ブランクは 2 群比較に入らない。`dataset_exclude` の除外（可逆）は前処理・
  図・強度順の検索が除く。
  - `dataset_status(view="samples")` の TSV の列は `name group role batch batch_source injection_order injection_order_source
    subject_id excluded`、続けて実験情報シートで足した因子の列（`Genotype` など）を名前の順に（値の無い試料は空欄。因子の列が
    無ければ増えない）。
  - `batch_source` / `injection_order_source`（値の出所）: `sheet`（実験情報シート・利用者）/ `mddata`（`.mddata`）/ `mzTab`
    （mzTab-M の MS:4000088・MS:4000089）/ `raw_acquisition_time`（raw の取得時刻の順）/ `filename`（ファイル名の順）。後の 2 つは、
    シートが schema の行の直後に宣言の行 `# injection_order_source = <raw_acquisition_time | filename>` を持つとき
    （`pipeline_sample_sheet` は取得時刻から注入順を作ったシートにこの行を書く）。
  - シートの空欄は `.mddata`（無ければ mzTab-M）の値で埋まり、その出所が出る（統計が使う値と同じ）。ただし宣言のあるシートが
    一部の行にだけ注入順を持つとき（取得時刻の読めない raw を `pipeline_sample_sheet` が空にした）、その空欄は埋めず不明のまま
    （シートの順位と MS-DIAL の番号は付け方が違う）。
  - 統計に入る（`include=true` の）生体試料か QC の注入順が欠ければドリフト補正は行わない（`include=false` の試料は前処理が先に
    外すので、その欠けは補正を止めない）。
  - バッチが空で試料名に 8 桁の数字があれば `batch_source` は `candidate:<8 桁>`（**候補だけで、統計には当てない**。§11.1）。
  - MS-DIAL Console が書く既定の `.mddata`（Class ID = FileID・種類はすべて sample）は宣言とみなさない: `group` は空、役割は
    ファイル名から推定する（役割の出所は `filename_token` / `inferred`）。
  - 解析フォルダが pipeline の run の中にあれば、`directory=` で開いても run の試料の対応表（実験情報シートと受領）から群と役割を
    当てる（注記が出る）。その単位が検証を通っていなければ、`pipeline_path` と同じ探索専用の印（`exploratory_reason` に run と
    単位の名前）も付く。
  - run の単位なら `load_dataset` は `processing` も返す: Console に渡した実効メソッドの `key_params`（最小ピーク高・ピーク幅・
    マススライス幅・RT と m/z の範囲・アライメントの許容幅・gap-fill の強制・`Peak count filter`・`N percent detected in one
    group`・付加体）、`undeclared`（メソッドに宣言が無く Console の既定で走ったキー。既定値そのものは書かない）、`overrides`
    （要求の `method_overrides` が変えたキーの `from` / `to`）。特徴の数・gap-fill の割合・1 化合物の分割はこれらで大きく変わる。
    値を変えて処理し直すのは `pipeline_run` の `request.method_overrides`（新しい run）。
- **MS/MS の実スペクトルは `.dcl`**（`dcl` トピック）。`.arf` のセルの MS/MS の有無は取得の参照があることだけを示す。
- **最初の概観**。`dataset_overview`（`load_dataset(overview=True)` の `overview` も同じ）は品質 → 全体 → クラス → 分子種の
  4 段の事実（ID つき）。2 群の向きは仮、クラスの段は付加体の選別の後だけ。§11.9。
- **EIC**。`.EIC.aef` から `dataset_plot(kind="eic" | "eic_compounds")` が描く（`eic` トピック）。`peak_top` は強度でなく頂点の座標。

以下の旧トピックの節で、撤去済みの関数・ツールは「撤去済み」「旧」と断ってある。

## 旧 `arf` トピック（`.arf` / サンプル別ピーク。§3, §8.1, §10, §11）

`.arf` の読み込み、PCA（`dataset_pca`）、前処理・QC（`dataset_preprocess`）、差次的解析（`dataset_compare`）の出力定義。

> Key 番号表（MS-DIAL の `[Key(N)]` から抽出した一次資料。リーダーのインデックス定数を変える前に確かめる）:
> `docs/schema/AlignmentChromPeakFeature.md`。

### 3. ARF (`metabolomix/arf/reader.py`)

#### 3.1 `deserialize()` のスポット出力

型は `list[dict]`。1 要素は 1 アラインメントスポットで、全サンプルのピークを `AlignedPeakProperties` に持つ。

| キー | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | パーサーがスポット順に付ける 0 始まりの ID |
| `AlignmentID` | int | `MasterAlignmentID` と同じ連番 |
| `RT` | float/null | 代表サンプルから取ったスポット代表 RT（min） |
| `MassCenter` | float/null | 代表サンプルから取ったスポット代表 m/z |
| `IonMode` | str | `Positive` / `Negative` / `Both` / `Unknown` |
| `Name` | str/null | 代表注釈名。空文字は未注釈 |
| `HeightAverage` | float/null | グループ先頭サンプルの `height`（名前に反して全サンプル平均ではない） |
| `AlignedPeakProperties` | list[list] | スポットの全サンプルの MessagePack 配列（`AlignmentChromPeakFeature`）。長さと Key 位置は生のままだが、解析で使わない Key は読み込み時に `null`（Key0〜10・18・20〜22 は生の値、Key15 は RT だけ、Key37 は先頭 2 要素）。通常は次節の表を使う |
| `TagIds` / `Tags` | list | アラインメント結果の `*_tags.xml` から `MasterAlignmentID` で結合したタグ |
| `SamplePeakTags` | dict | サンプル別 `*_tags.xml` から `FileName` と `MasterPeakID` で結合した、タグ付きピークだけの辞書 |
| `SampleClasses` | dict | `.mddata` の `AnalysisFileBean` と `FileID` / `FileName` で結合したサンプルの Class ID メタデータ |

#### 3.2 `extract_peak_properties()` の表/CSV

型は `pandas.DataFrame`。1 行は 1 スポット × 1 サンプルのピーク（同じ `MasterAlignmentID` がサンプル数だけ繰り返す）。

| 列 | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | 行が属するスポットの ID |
| `AlignmentID` | int | `MasterAlignmentID` と同じ値 |
| `SpotRT` | float/null | スポットの代表 RT（min） |
| `SpotMassCenter` | float/null | スポットの代表 m/z |
| `IonMode` | str | スポットのイオンモード |
| `CompoundName` | str/null | スポットの候補化合物名。空文字は未注釈 |
| `AlignmentTags` | list[str] | スポットに付いた MS-DIAL タグ |
| `SampleIndex` | int | `AlignedPeakProperties` 内の 0 始まりの位置（サンプル順） |
| `FileName` | str | 測定ファイル名。取れなければ `Sample_<SampleIndex>` |
| `ClassID` | str/null | MS-DIAL の File property setting の `AnalysisFileClass` |
| `PeakID` | int/null | 測定ファイル内のピーク ID。ギャップフィルでは負値になり得る |
| `FileID` | int/null | データセット内の測定ファイル ID |
| `MasterPeakID` | int/null | 元ピークのマスター ID。負値は未検出を補間したギャップフィル |
| `PeakHeight` | float/null | そのサンプルのピーク頂点強度 |
| `PeakArea` | float/null | そのサンプルのピーク面積 |
| `PeakAreaAboveBaseline` | float/null | ベースラインより上のピーク面積 |
| `PeakMZ` | float/null | そのサンプルで観測したピーク m/z |
| `PeakRT` | float/null | そのサンプルで観測したピーク RT（min） |
| `SignalToNoise` | float/null | `PeakShape[1]` 由来の S/N |
| `IsMsms` | bool | MS2 raw spectrum ID-to-collision-energy map が非空か。MS/MS 取得の情報があるだけで、スペクトル本体は含まない |
| `IsGapFilled` | bool | `MasterPeakID < 0` か。`true` は実検出でなくギャップフィルの値 |
| `PeakTags` | list[str] | そのサンプルピークの `MasterPeakID` に付いた MS-DIAL タグ |

#### 3.3–3.4 `.arf` からの PCA 行列と `run_pca()`（撤去済み）

`build_pca_matrix()` と `run_pca()` は撤去済み（コードに無い）。PCA は `dataset_pca`（§10.4。前処理後の log2 行列で、軸は生体の
個体だけで決め QC は投影する）。

#### 3.5 ARF要約

`summarize_arf_data()` の返す辞書。

| キー | 意味 |
|---|---|
| `total_peaks` | スポット数（サンプル別の行数ではない） |
| `rt_range` | スポット代表 RT の `(min, max)` |
| `mass_range` | スポット代表 m/z の `(min, max)` |
| `height_average_mean` | `HeightAverage` の算術平均 |
| `height_average_max` | `HeightAverage` の最大値 |
| `ion_modes` | イオンモード別のスポット数 |
| `named_compounds` | `Name` が null でない件数。**空文字も数える**ので注釈済みの件数とは限らない |

#### 3.6 MS-DIALタグ

`metabolomix/msdial/tags.py` が ARF と同じディレクトリの `*_tags.xml`（拡張子なしの `*_tags` も可）を読む。サンプル別ファイルは
処理時刻の 12 桁接尾辞を除いた名前で ARF の `FileName` と対応させ、`Peak/@Id` を `MasterPeakID` と結合する。アラインメント結果用
ファイルは `Peak/@Id` を `MasterAlignmentID` と結合する。

タグ条件は `any` / `all` / `none` / `not_all`。`sample_peak` スコープは条件に合わないサンプル別行をスポットから除き、
`alignment_spot` スコープはスポット全体を除く。タグの無いピークは `any` / `all` に一致せず、`none` / `not_all` に一致する。

サンプル名は完全一致を優先し、12 桁接尾辞を除く照合は一意に決まるときだけ使う。重複・曖昧な照合はエラー。タグファイルの無い
サンプルは既定でエラー、`missing_sample_policy=exclude` で除外、`untagged` でタグなし扱い。サンプル数は ARF から検出した数。

#### 3.7 Class ID

`metabolomix/msdial/classes.py` が `.mddata` の `MsdialDataStorageBase.Key0 AnalysisFiles` から各 `AnalysisFileBean` の
`AnalysisFileId`・`AnalysisFileName`・`AnalysisFileClass` を読む。`.mddata` は明示パス・`.mdproject` 内の参照・ARF と同じ
ディレクトリのいずれかから解決する。

ARF サンプルとは `FileID` で、欠ければ正規化した `FileName` で結合する。両方が別のサンプルに解決されればエラー。
`filter_arf_by_class_ids()` は選んだ Class ID 以外のサンプル行を各スポットから除く。

#### 8.1 PCA の返り値と、生スポットの選択（旧 `.arf` パーサ）

旧 `.arf` パーサ（撤去済み）の PCA と生スポットの選択は、次のツールが受け持つ。

- 前処理後の行列の PCA は `dataset_preprocess` → `dataset_pca`（§10.4）。`dataset_exclude` の除外は行列を作る前に反映され、
  変えると前処理済みの行列から先の結果が古くなる（§10.5）。
- 部分集合で計算し直すと、TIC 正規化の分母と BH 補正の検定対象数が変わる。全体で解析した値と混同しない。
- Class ID 別・役割別の試料数とタグの件数は `dataset_status(view="summary")`、タグごとの特徴 ID は `view="tags"`。
- スコア図は `dataset_plot(kind="pca")`（既定は PNG。`output="payload"` で座標の点列を返し、散布図はクライアントが描く。
  `output_path` で任意の場所にも書く）。点列は `ds.figures["pca"]` に残り、`save_figure(kind="pca")` がそこから報告書用の PNG を
  `reports/figures/` に書く。

### 10. 前処理・QC（P2a）

統計エンジン（`metabolomix/analysis/stats/`。MCP 非依存）を `metabolomix/tools/stats_tools.py` の `dataset_preprocess` /
`dataset_compare` / `dataset_pca` が呼ぶ（行列の登録と前処理の記録は `metabolomix/tools/_matrix_store.py`）。読み込んだデータセットには
**何も適用されていない**: `dataset_preprocess` を呼ぶか、`dataset_compare` / `dataset_pca` が既定の設定で自分で回す。

#### 10.1 役割検出（sample/qc/blank）

`preprocessing.detect_sample_roles()` がファイル名と Class ID を `_` で区切り、大小無視で `qc` / `blank` トークンと照合して
`sample` / `qc` / `blank` に分ける（`blank` を先に評価）。`dataset_status` が分類と役割別の件数を返す。

#### 10.2 前処理（`dataset_preprocess`。統計エンジン）

`dataset_preprocess` は前処理済みの行列を `ds.stat_matrices[matrix_id]` に登録し `matrix_id` を返す。
**同じデータ・除外・群・設定なら同じ `matrix_id`**。工程は **選択 → ドリフト補正 → 品質の足切り → 正規化 → log2 → 反復のまとめ**
（個体ごとの行列）。0 以下・非有限は欠測（pseudocount も clip も使わない）。補完は `dataset_pca` が PCA のためだけに行う
（行列は補完しない）。

| 引数 | 既定 | 意味 |
|---|---|---|
| `drift_mode` / `drift_span` | `auto` / 0.75 | QC の LOESS ドリフト補正（`auto` / `on` / `off`）。`auto` は条件（QC の数・挟まり方・注入順）が揃えば当てる。ただし換算済みの値（`quant_normalized`）か、注入順がファイルの順かもしれない（MS-DIAL の既定、またはシートの宣言で出所 `filename`。どちらも `injection_order_may_be_file_order`）ときは**当てずに**計算し、尋ねない（`needs_input` は返さない）: ドリフトの段は `status: skipped`・`reason`・`reasons`・`rerun_with: {"preprocess.drift.mode": "on"}`、`warnings` に `drift_not_applied_by_default`（`reasons`・`message`・`suggestion`）。条件が揃わない `auto` は `reasons` に理由を残して飛ばすだけ（呼び直しの案内は出さない）。`on` は条件が揃えば必ず当て、揃わなければ `PREPROCESS_PRECONDITION` |
| `loess_degree` | 1 | LOESS の局所多項式の次数（`1` = 局所線形。statsmodels の lowess と同じ / `2` = 局所 2 次。同じ近傍・三乗重みを numpy で解き、近傍が 1 点多く要る）。ドリフトの段の記録に `loess_degree`。`drift_mode="off"` では使わず `warnings` に `loess_degree_ignored`（`reason: "drift_off"`。行列と `matrix_id` は渡さないときと同じ） |
| `blank_min_ratio` | 3.0 | 生体試料の平均がブランク平均のこの倍数未満の特徴を除く |
| `detection_min_rate` / `detection_by` | 0.5 / `group` | gap-fill を除いた検出率の下限（`group` は群ごと、`all` は全体）。検出状態が無いデータでは評価できない旨が steps に出る |
| `max_qc_rsd` / `max_d_ratio` | 0.30 / 0.5 | QC の RSD の上限・D-ratio の上限 |
| `normalize` | `auto` | `auto`（換算済みか QC が無ければ正規化しない。それ以外は `pqn`）/ `none` / `tic` / `median` / `pqn`。`auto` が QC の無さで正規化しないときは正規化の段が `status: skipped`・`reason: "no_qc"`、`warnings` に `normalize_auto_no_qc`。PQN の参照スペクトルは QC があれば QC の、QC が無いデータで `pqn` を明示すると**生体試料（全群）**の特徴ごとの中央値で、正規化の段の `reference` は `qc` / `sample`。`sample` のときは `reference_definition: "median_of_samples_all_groups"` と注意 `pqn_reference_is_samples`（群の間で全体の組成が大きく違うと、係数がその差の一部を吸収しうる） |
| `normalization_features` / `normalization_exclude_features` | なし | 正規化の係数を求める特徴を絞る・除く（特徴 ID のリスト）。既定は足切りの後に残った特徴すべて。**係数から外した特徴も行列に残り**、同じ係数で割られる。`normalization_features` の特徴は足切りで落ちていても係数に使う（値はドリフト補正の後）。正規化の段に `factor_features`（係数に使った数）・`factor_feature_basis`（`features` / `kept_features`）・`excluded_from_factors`・`factor_features_filtered_out`（足切りで落ちていた特徴）。無い ID・係数の特徴が 1 つも残らない指定は `PREPROCESS_CONFIG_INVALID`（陽陰を統合したデータセットで一部の極性だけが空になる指定は断らない。下の「陽陰を統合したデータセット」）。空のリストは指定なしと同じ（`normalization_features=[]` は `warnings` に `normalization_features_empty`）。効かないとき（`normalize="none"`、`auto` が換算済み・QC 無しで正規化しない）は使わず、`warnings` に `normalization_features_ignored`（`reason`: `"method_none"` / `"quant_normalized"` / `"no_qc"`。`auto` のときは正規化の段に `factor_features_ignored: true`）。行列と `matrix_id` は指定なしと同じ |
| `collapse_replicates` | true | 同じ個体（`subject_id`）の反復注入を 1 個体にまとめる |
| `injection_columns` | `[]` | 個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う列。技術反復が別の個体に割れるのを防ぐ） |
| `group_levels` / `group_factors` | なし | 試料名（と試料の表の群）の因子トークンで群を解く。解いた群は行列の設定に記録される |

足切り 4 つは `None` で行わない。`dataset_exclude` の除外は自動で反映される（行列にも検出率の分母にも入らない）。

戻り値（要約。全量は `ds.stat_matrices`）: `{status, matrix_id, n_units, n_qc, n_features_before, n_features,
subject_source, groups: {群: 個体数}, quant, steps, warnings, record, next}`。

- `steps[]` は工程ごとに `{step, status: applied|skipped, reason?, ...}`。`filter` は足切りごとの内訳
  （`filters.<blank_ratio|detection|qc_rsd|d_ratio>` の `status` と落とした数。独立のマスクで重なり得るので足さない）と
  全体の `dropped`・`n_features_in`。`collapse` は `subject_source`（`subject_id` / `injection` / `mixed`）・`n_injections`・
  `n_units`・`collapsed_units`。
- `groups` の群が 1 つだけなら警告 `single_group`（`group_levels` / `group_factors` の綴りか試料の表の群を確かめる）。
- 注入順とバッチの注意:
  - `injection_order_may_be_file_order`（`n_samples`）: 注入順が MS-DIAL の振った値（`.mddata` / mzTab-M）で、ファイル名の順に
    並べると狭義に増える（同じ注入順の試料があれば立てない）。MS-DIAL は入力が無ければファイルの順で振るので、実際の測定順とは
    限らない。同じ注意を交絡の判定（`dataset_compare` の `confounding.warnings` に `n_units` 付き。エンジンの結果〔`dataset_analysis_run`
    の比較も〕は `confounding.injection_order.may_be_file_order` も立てる。交絡の判定はシートの宣言
    `# injection_order_source = filename` の注入順〔出所 `filename`〕にも付ける）と、質量精度の注入順の ρ
    （`dataset_mass_accuracy` の `injection_trend.warnings`）も返す。
  - `batch_candidates`（`values`）: バッチが不明な試料の名前の 8 桁が 2 種類以上。ドリフト補正・反復のまとめ・交絡には当てない。
    バッチなら試料の表の `batch` に入れる。
- **陽陰を統合したデータセット**（`dataset_merge_polarities`。全特徴が `feature_metadata[<fid>]["origin"]["polarity"]` を持つ）は、
  **ドリフト補正と正規化を極性ごと**に当てる（各極性の値・係数は、その極性だけを前処理して並べたものと同じ）。足切り
  （ブランク比・検出率・QC RSD・D-ratio）・log2・反復のまとめは特徴ごと・セルごとなので分けなくても同じ。
  - ドリフトはその極性の測定の注入順・バッチで当てはめる。陽イオンは統合したデータセットの試料の表（＝陽イオンの試料の表。
    陽イオンに当てたシートも統合が引き継ぐ）、陰イオンは統合が陰イオンの試料の表（当てたシートを含む）から読んで残した
    `quant_state.polarity_runs.negative`。`drift_mode="on"` は注入順が確かでない極性もその注入順のまま当てはめる。
  - 正規化の係数（`tic` / `median` / `pqn`、`auto` が解いた方法も）はその極性の特徴だけから求め、その極性の列だけを割る。
  - 記録: ドリフトの段は `by_polarity`（極性 → その極性の段）と、`status`（全極性で補正 `applied` / 一部 `partial` / どれも当てない
    `skipped`）・`reasons`（当てなかった極性の理由）・`polarities_skipped`。`rerun_with` は当てなかった極性がすべて drift=auto の
    保留のときだけ。`drift_mode="on"` で条件を満たさない極性があれば `PREPROCESS_PRECONDITION`（message の先頭と `details.polarity`
    に極性）。正規化の段は方法・参照・`factor_samples` を上に、`by_polarity` に極性ごとの `factors`・`factors_set_to_one`・
    `unscaled_samples`・係数の特徴の記録。
  - `normalization_features` / `normalization_exclude_features` は極性ごとにその極性の特徴だけを使う（ほかの極性の特徴では割らない）。
    ある極性に挙げた特徴が 1 つも無ければ、その極性は残った特徴すべてから求め
    （`factor_features_fallback: "no_listed_feature_in_polarity"`）、注意 `normalization_features_not_in_polarity`（`polarities`）。
    除外（と `normalization_features`）でその極性の係数の特徴が残らなければ、同じく残った特徴すべてから求め（その極性に除外は当てない。
    `factor_features_fallback: "all_factor_features_excluded_in_polarity"`）、注意 `normalization_exclude_features_cover_polarity`
    （`polarities`）。全極性で空になる指定は `PREPROCESS_CONFIG_INVALID`。
  - 注意 `injection_order_may_be_file_order` と `drift_not_applied_by_default`（`polarities`: 極性 → 保留の理由）は極性付き。陽陰で
    注入順・バッチが違う（片方だけ分からない）試料があれば `polarity_runs_differ`（`polarity`・`reference_polarity: "positive"`・
    `n_injection_order_differs`・`n_injection_order_unknown`・`n_batch_differs`・`samples`）。統合が測定を残していない（古い統合）なら
    `runs_recorded: false` で、その極性のドリフト補正は当てない（統合し直すと残る）。
  - 統合したハンドルへ当てるシート（`dataset_set_sample_metadata`）は陽イオンの raw を挙げるので、`batch`・`injection_order` は
    **陽イオンの測定だけ**を決める（群・役割・include・subject_id・因子の列は試料ごとなので両極性に効く）。陰イオンの測定を変える
    なら、陰イオンのハンドルに陰イオンの raw を挙げたシートを当ててから `dataset_merge_polarities` で統合し直す（統合したハンドルへの
    シートの戻り値の `notes` も言う）。統合し直すと、前の統合したハンドルへ当てたシートと除外（`dataset_exclude`）は引き継がれない
    （統合が引き継ぐのは陽イオンのハンドルの試料の表と、極性ごとのハンドルの除外）。
  - 上の注意の `message` は、その極性の測定を入れ直す手順を、上から行えば全部が入る順（陽イオンを最後）に添える: 陰イオンは
    統合し直す手順（前の統合したハンドルに当てていれば、統合し直してできたハンドルに当て直す）が先、陽イオンは統合したハンドルへの
    シート（ほかの極性を統合し直すなら、統合し直してできたハンドルに入れ、除外も当て直す）が後。`drift_not_applied_by_default` の
    `message` も極性ごとの `injection_order_may_be_file_order` の並びも同じ（`dataset_overview` の但し書き・注記と同じ文と順）。
  - raw のパスが無いデータセットの「注入順がファイルの順か」は、陽イオンの元の試料名（統合は試料を共通 ID へ改名する。
    `quant_state.source_sample_names`）で見る（陽イオンだけのデータセットと同じ判定）。
  - 交絡の判定・比較の共変量 `batch`・反復のまとめは統合したデータセットの試料の表（陽イオンの測定）を使う。`matrix_id` には
    極性の塊（特徴の分け方と極性ごとの測定）が入る（統合でないデータセットの `matrix_id` は変わらない）。
- `quant` は定量で作ったデータセット（内部標準で換算済み）のときだけ `{normalized, method, unit}`、それ以外は null。
- `record` は `reports/stats/` への記録 1 件（`{written: true, record_id, path}`）。書けなければ `{written: false, reason}` で、
  前処理は成功している。記録は追記専用（`lipidmix.stats-record.v1`: `step` は `preprocess` / `compare` / `pca` / `overview` /
  `analysis_run`〔`dataset_analysis_run`。`reanalysis.md` §2〕、`inputs`・`parent_ids`・`settings`・`summary`）で、比較・PCA の記録は
  前処理の記録を `parent_ids` に持つ。
- **古くなる**: `dataset_exclude`・試料の表の更新（`dataset_set_sample_metadata`。`group` だけの変更も）・データの入れ替えで、
  登録した行列は古くなり（`matrix_id` を渡すと `missing_state`。`dataset_preprocess` からやり直す）、その行列から出した比較・
  PCA の結果も現在のものとして通らない。

> **役割別の行の扱い（重要）**: 前処理後の行列の個体は**生体試料**。QC は行列に残るが軸の計算には使わず、`dataset_pca` が同じ軸へ
> 投影する（スコア図で白抜き、群ラベル `QC`。QC の締まりで品質を見る）。**ブランクは比較にも PCA にも入らない**（`blank_min_ratio`
> の参照にだけ使う。総強度が桁違いに低く、残すと PC1 を支配する）。2 群比較は QC を比較群から外す（§11.1）。

#### 10.3 caveatの扱い

前処理で飛ばした工程と理由は `steps[]`（`status: skipped` と `reason`）に、注意は `warnings[]`（`{code, message, ...}`）に残る
（QC・ブランク・注入順が欠けて行えなかった工程も黙って飛ばさない）。解釈と報告の注意点として引用する。

**失敗した QC 注入**（総量が QC の中央値の 0.2 倍未満）は最初の概観の品質の段（§11.9）が名指しする。残したまま `max_qc_rsd` を
掛けると QC の RSD が全特徴で跳ね上がり、ほぼ全特徴が落ちる（実測 1345 → 51）。足切りが厳しすぎるように見えたら QC 側を疑い、
`dataset_exclude` で除外して前処理し直す。

`load_dataset` の複数バッチの告知は `.arf` / `.arf2` のバッチだけに基づく（`.mddata` / `.mdproject` / `.msp2` / `.pai2` などの
`AlignmentResult` 形式のタイムスタンプは数えない）。

#### 10.4 `dataset_pca`

統計エンジンの PCA（`stats-pca.v1`）。`matrix_id`（`dataset_preprocess` の行列。省くと、そのハンドルで `dataset_preprocess` が最後に
作った行列、無ければ既定の前処理）の log2 行列で、**軸は生体の個体だけで決め、QC は同じ軸へ投影する**。欠測は PCA のためだけに
補う（`impute`）。`group_levels` / `group_factors` は色分けと記録にだけ使い、スコア・寄与率は変わらない。

戻り値は要約 `{status, result_id, matrix_id, matrix_source, scaling, log2: true, imputation, n_fit, n_qc_projected,
explained_variance_ratio, groups: {群: 個体数}, record, loadings_note}`（`matrix_source` は `matrix_id` / `latest_preprocess` /
`default`）で、スコア・ローディングの全量は `result_id`（`ds.last_pca`）。図は `dataset_plot(kind="pca" / "loadings")`（QC は白抜き、
ローディングは軸の計算で値が変わった特徴だけ）。

- 渡した `matrix_id` が古い（除外・試料の表の変更）・無いときは `missing_state`（`dataset_preprocess` を案内）。
- `matrix_id` を省いた行列の前処理は尋ねない: `drift_mode=auto` がドリフト補正を当てなかった行列なら `warnings` に
  `drift_not_applied_by_default`（補正するなら `dataset_preprocess(drift_mode="on")` で作った `matrix_id` を渡す）。既定の前処理で
  行列を作ったときは、その行列を登録して前処理の記録も書く（エラーのときも `details.matrix_id` に残る）。
- 結果（`stat_pca`）は除外・試料の表の**どんな**更新（`group` だけの変更も）でも古くなる（色分けが焼き込まれているため）。

#### 10.5 手動サンプル/ピーク除外（`dataset_exclude`）

外れた試料や特定の特徴（スポット）を名前 / ID で除外する。`metabolomix/analysis/dataset_exclusions.py`（MCP 非依存）と
`metabolomix/tools/sample_tools.py` の `dataset_exclude`。除外は**可逆・非破壊**（読み込んだ定量行列・特徴 ID・試料名は変えない）。

`dataset_exclude(samples=None, features=None, mode="add", reason=None)` は結果の件数を返す。

- **`samples`**: 試料名（`dataset_status(view="samples")` の name、完全一致）のリスト。
- **`features`**: 特徴 ID（`dataset_feature_search` の feature_id。`MasterAlignmentID`）のリスト。
- **`mode`**: `add`（既定。`reason` が要る）/ `remove`（戻す）/ `clear`（全消去）。現状は `dataset_status(view="samples")` の
  `excluded` 列。
- 無い試料名・特徴 ID が 1 つでもあれば全体を拒み、何も変えない。全試料の除外も拒む。
- 除外が実際に変わったときだけ、登録した行列（`stat_matrices`）を空にし、比較・PCA の結果を古くする（`dataset_preprocess` から
  やり直す）。同じ引数の再実行は何も変えない（冪等）。
- 行列・結果の一覧と、それぞれが今も通用するか（`current`）は `dataset_status(view="results")`（列の意味は `mztab` トピック §13.5b）。

除外は `excluded_samples`（試料名 → 理由）と `excluded_features`（特徴 ID → 理由）に持ち、変えると `metadata_revision` が進む。
`dataset_preprocess`（→ `dataset_pca` / `dataset_compare`）・`dataset_plot`・`dataset_feature_search` の強度順が除外を反映する。
前処理の記録に残る（select の工程の `excluded_by_include` の件数、記録の `settings.exclusions`）が、前処理の出力に専用の注記は付かない。

流れ: `dataset_preprocess` → `dataset_pca` → `dataset_exclude(samples=[...], reason=...)` → `dataset_preprocess` ＋ `dataset_pca` →
`dataset_compare`。戻すには `mode="remove"` / `mode="clear"`。

### 11. 差次的解析（P2b）

比較は統計エンジンの `compare`（`metabolomix/analysis/stats/compare.py`。MCP 非依存）を `metabolomix/tools/stats_tools.py` の
`dataset_compare` が呼ぶ。図（`dataset_plot(kind="volcano")`・`metabolomix/tools/reports.py` の `save_figure(kind="volcano")`）と
書き出し（`dataset_export_differential` / `dataset_export_evidence_bundle`）は、`dataset_compare` の結果を対比 1 つの 2 群の形
（`stats_results.as_two_group`）にして作る。`metabolomix/analysis/differential.py` に残るのは volcano の点列 `volcano_data()` だけ。

#### 11.1 群ラベルの由来

群ラベルは試料の情報の `group`（`dataset_status(view="samples")`。ファイル名由来の factor トークン / Class ID、
`metabolomix/msdial/classes.py` の `assign_sample_groups`）。バッチは試料の情報の `batch` で、実験情報シート・利用者の訂正・
MS-DIAL の宣言（`.mddata`、mzTab-M の MS:4000088。全試料で同じ値の宣言は GUI の既定とみなして採らない）の値だけを使う。
**試料名の 8 桁の数字（取得日らしい）はバッチに当てない**（`batch_source` の `candidate:<8 桁>` と注意 `batch_candidates` の候補だけ。
§10.2）。バッチが無ければ交絡は「評価不可」。
比較（`dataset_compare`）は行列（`matrix_id`。省くと、そのハンドルで `dataset_preprocess` が最後に作った行列、無ければ既定の前処理）を使う。

試料の `group` は常に**完全な Class ID**（例 `24M_GF_F`）。`dataset_compare` は群名（`reference_group` / `test_group`）を受ける。
因子トークンで試料名を引くのは `dataset_status(view="samples", query=...)`（`sample_factors.expand_sample_specs` の文法。トークンは
Class ID と試料名の両方から解く）:

- `24M` → `24M_*` のすべて
- `24M_GF` のように `_` で繋ぐと AND（`24M` かつ `GF`）
- 完全な Class ID はその 1 水準だけ
- Class ID が同じで試料名の因子（時点など）だけが違う群も切り出せる（例 `ILG_6h` と `ILG_0h`）

**一致ゼロ・両群が同じ Class ID を掴む指定は `status="error"`**（n=0 を成功で返すと「有意 0 件＝群間差なし」と誤読されるため）。交互作用は
`dataset_compare` で検定する（§11.2: `comparisons` の `kind: "interaction"` の 2×2 の対比〔差の差〕と、因子が 2 つ以上で cell が
そろうときの `factorial` の主効果・交互作用の F）。

**QC・ブランクは比較群に入らない**。`role` が `qc` / `blank` の試料は群ラベルを `None` にしてから展開する（QC の Class ID が群の
トークンを含んでも〔例 `QC_24M` と `24M`〕群平均に混ざらない）。除外の内訳は caveat「比較対象から除外（生体試料でないため）」、
`n_a` / `n_b` は除外後の実 n。

#### 11.1.1 上位ヒットの命名（ARF/ARF2 橋渡し）

`summary.top` の各行に `spot_id` / `name` / `name_source`。`name_source="arf"` は ARF スポットの `Name`、`"arf2"` は ARF が
`Unknown` のときに**同じアラインメントの兄弟 `.arf2`** から補った注釈（`ontology` も併記）。ARF と ARF2 は同じ `MasterAlignmentID`
で代表 `Name` が食い違うことがある（例: ARF `Unknown`、ARF2 `SL 33:0;O|SL 17:0;O/16:0`）。

補完元は GUI の `AlignmentResult_<timestamp>` か Console の `AlignResult-<digits>` の語幹が完全一致する隣の `.arf2` だけ
（`MasterAlignmentID` はアラインメントごとに振り直される）。兄弟が無ければ補わず `name=null`。

#### 11.2 統計（`dataset_compare`。統計エンジン）

`dataset_compare` は 2 群・多群（3 群以上。先頭が基準で各群 − 基準の対比と omnibus F）・2 因子（`comparisons` の `at`・交互作用）を
比べる。群の指定は **reference_group+test_group / groups / comparisons のちょうど 1 つ**。既定の検定は **moderated t**（経験ベイズで
分散を縮める線形モデル。`welch` / `student` / `mannwhitney` は 1 因子・共変量なし・ブロックなしだけ）。`covariates`（試料の表の
列名。値がすべて数値なら numeric、さもなくば categorical）・`block`（`"subject_id"` だけ）を受ける。補正の母集団 `bh_population` は
`matrix`（既定）/ `annotated` / `reported`（`reported_feature_ids` の特徴だけ。ID の一覧は記録の spec の `reported_feature_ids`）。
多重比較の補正 `p_adjust` は `bh`（既定。FDR）/ `bonferroni` / `holm`（FWER）/ `none`（補正なし）で、行の `q` はこの方法の値
（下の表の「BH の `q`」は既定の場合の言い方）。契約の書き出し（`dataset_export_differential`・根拠バンドルの `features.tsv`）の
`q_value` は方法によらず BH の q（§11.4）。

**行列（`matrix_id`）**: 渡せばその行列（その前処理・群・除外）で比べる。省くと、そのハンドルで `dataset_preprocess` が最後に作った
行列（`group_levels` / `group_factors` を渡したときは群が同じ行列だけ。除外・試料の表の変更で古くなった行列は数えない）、無ければ
既定の前処理で行列を作って登録する。使った行列は戻り値の `matrix_id` と `matrix_source`（`matrix_id` / `latest_preprocess` / `default`）。

**判定（`verdict`）は 3 値**（`q_threshold` 既定 0.05、`equivalence_fc` 既定 1.2、CI は 95%）:

| 値 | 条件 |
|---|---|
| `changed` | 検定でき、BH の `q ≤ q_threshold` |
| `unchanged` | `q` は閾値を超え、**95% CI が同等の範囲 `±log2(equivalence_fc)` に完全に収まる** |
| `inconclusive` | それ以外（検定できない・検出の割合が足りず量を判定できない・`q` が出ない・CI が出ない・CI が同等の範囲を超える）。`verdict_reason` に `not_tested` / `quantity_not_judgeable` / `outside_bh_population` / `no_ci` / `neither` |

**「有意でない」は「変わっていない」ではない。** `q` が大きいだけでは `unchanged` にならない。n が小さい・ばらつきが大きいと大半が
`inconclusive`（差の有無でなく、データが足りない）。`log2fc_threshold`（既定 1.0）は判定に使わない（`large_changes` の数え方・
図の色分け・書き出しの `significant` 列だけ）。

戻り値（要約。全行は `result_id` でセッションに保持し、直近の結果＝図・書き出しの既定になる）:
`{status, confounding, result_id, matrix_id, matrix_source, rule_hash, test, p_adjust_method, effect_size_definition: {name},
ci_level, verdict_rule: {q_threshold, equivalence_fc, equivalence_log2}, bh_population: {name, n_by_contrast,
q_interpretation, p_adjust_method}, subject_source, group_n, prior, contrasts, omnibus, factorial, warnings, top, f_top?,
output?, record, next, notes?}`。

- `rule_hash` は結果の spec（analysis-spec v2 に写して正規化したもの）の hash、`result_id` はそれと `matrix_id`・`log2fc_threshold`
  の hash。`p_adjust` が既定の `bh` なら spec に書かない。
- `confounding`（**先頭**）: `{batch_status, batch_table, by_contrast: {対比 ID: {batch_status, confounded}}, warnings}`。
  完全交絡（どのバッチも 1 群だけで、群とバッチを分離できない）でも**既定では止めずに比べ**（`allow_confounded` 既定 true）、全対比を
  まとめた `batch_status = complete` と対比ごとの `by_contrast[...].confounded`（多群では一部の対比だけが完全交絡になりうる）で
  印を付ける。その対比の群差はバッチ差と区別できない（結論に使わない）。`allow_confounded=false` なら、どちらかの完全交絡で
  `STATISTIC_CONFOUNDED`（`details.contrasts` に完全交絡の対比）。`confounded` は `true`（その対比の群だけで見て完全交絡）/
  `false`（部分交絡・交絡なし・単一バッチ）/ `null`（バッチ情報が無く「評価不可」。「交絡あり」でも「なし」でもない。TSV では空欄）。
  部分的な交絡は警告に出る。
- `contrasts[]`: `{contrast_id, reference, test, confounded, batch_status, n_tested, n_not_tested, verdicts: {changed,
  unchanged, inconclusive}, large_changes, presence: {available, n_q_below}}`。`presence` は検出の有無そのものの検定（gap-fill を
  除く検出状態があるときだけ `available`）。`omnibus` は 3 群以上（cell が 3 つ以上）のときの全 cell の F 検定の
  `{n_tested, n_q_below}`（無ければ null）。
- `factorial`: 因子が 2 つ以上で全水準の組（cell）がそろうときの**主効果と交互作用の F**（Type III。多水準を含む moderated F）の
  効果ごとの `{df1, n_tested, n_q_below}`（効果の名前は因子名、交互作用は `group:sex` のように `:` で結ぶ）。無ければ null。F は
  向きと log2fc を持たないので対比の行とは別の欄。
- `top`: 対比ごとに `changed` の上位（q の小さい順。`top_n` 既定 10、50 を超える値は断らず 50 に下げて `notes` に書く。そのとき
  50 行で切った対比があれば `top_truncated: true` と `top_n_more`〔対比 → 載せなかった changed の行数〕、`f_top` で切った効果が
  あれば `f_top_n_more`〔効果 → 行数〕）の TSV（`contrast_id feature_id name id_level log2fc ci_low ci_high q verdict`）。`name` は
  `feature_metadata`、無ければ MS1 注釈（`feature_annotations`）の名前、無ければ空文字。`id_level` は `msms` / `ms1` / `none`。
  `f_top` は omnibus・要因の効果ごとに `f_q ≤ q_threshold` の上位の TSV（`effect feature_id name id_level f df1 df2 f_q`。`effect` は
  `omnibus` か効果の名前）で、F が無い・該当が無ければ欄ごと無い。**全行は `output_path` か `result_id` からの書き出しで取る**。
- `output`（`output_path` を渡したときだけ）: `{written, path, n_rows, ftests_path, n_ftest_rows}`。`output_path` に**全特徴 ×
  全対比の全列**の TSV を書く（列は `contrast_id feature_id` ＋同定の列 `name ontology inchikey id_level annotation_source
  msms_assigned msi_level` ＋比較の行の全欄 `reference test status reason n_reference n_test detected_reference detected_test
  quantity_judgeable constant_in_group confounded mean_reference mean_test log2fc se ci_low ci_high ci_unavailable_reason statistic
  df p q presence_p presence_q verdict verdict_reason`。数は有効 10 桁、真偽は `true` / `false`、無い値は空欄、一覧は `,` で結ぶ。
  メタ行は無い）。omnibus・要因の F があれば `<stem>.ftests<拡張子>` にも書く（`kind effect feature_id name id_level f df1 df2
  f_p f_q`。`kind` は `omnibus` / `factorial`）。`df2`（F の分母の自由度＝`df_prior + df_residual`）が `inf` なのは `df_prior` が ∞
  のとき（F の参照分布は χ²(df1)/df1 の極限。`f_top` も同じ）。空欄は F を計算していない行。`msi_level` は解析フォルダの
  データセットだけ（mzTab-M では空欄）。**差次的エクスポート（`dataset_export_differential`。§11.4 の 15 列の契約）とは別物**。
  書けなくても比較は成功し、`written=false` と `reason`。`output_path` がフォルダ・空文字なら比較の前に `DATASET_BAD_REQUEST`。
- `bh_population.q_interpretation`: `discovery_fdr`（`p_adjust=bh` で `matrix` / `annotated`）・`fwer`（`bonferroni` / `holm`）・
  `unadjusted`（`none`）・`reproduction_only`（`reported`）。`p_adjust_method` は補正の方法（記録にも残る）。
- `effect_size`（`log2_mean_difference` 既定 / `arithmetic_mean_ratio`）の定義は `effect_size_definition`（`name` と、全量の記録の
  `description`）。`log2fc` は**正なら test が高い**。`mean_*` は個体ごとの前処理後の線形値の算術平均。行の `se` は log2fc の標準
  誤差（moderated t・Welch・Student。順位の検定と算術平均の比では空）。
- 行の `constant_in_group`: cell の中で値がすべて同じ（埋め値の疑い。gap-fill を同じ値で埋めた群など）なら、その cell の一覧
  （無ければ null）。調べるのは cell の個体の**読んだ値**＝注入ごとの値（ドリフト補正・正規化・反復のまとめの前）で、有限の正値が
  2 注入以上あり、すべて同じ（相対の丸めの誤差の内）なら挙げる（`dataset_analysis_run` v1 の `constant_in_group` と同じ規則）。
  検定する個体の log2 では探さない（正規化の係数が埋め値を散らし、特徴の少ない PQN は埋め値でない値を一定に見せる）。判定は変えず、
  件数を警告 `constant_in_group`（`{n_rows, n_features, by_contrast, examples}`）に出す。片群が一定だと、群の差は埋め方で決まりうる。
- エラーと問い:
  - `DATASET_BAD_REQUEST`: 群の指定が合わない・`top_n`（整数でない・負）/ `log2fc_threshold` / `p_adjust` / `output_path` が不正・
    `bh_population=reported` なのに `reported_feature_ids` が無い。`reported` 以外で渡した `reported_feature_ids` は誤りにせず、
    使わずに `notes` に書く。
  - `matrix_id` が古い・無い: `missing_state`（`stat_matrix`）。
  - `matrix_id` と一緒に渡した `group_levels` / `group_factors` は誤りにしない: 前処理に入る（除外していない）試料を行列と同じ群に
    分けるなら `notes` に書いて比べる（登録の `group_overrides` の dict ではなく試料ごとの群で比べるので、試料の表の群〔Class ID〕と
    同じ分け方の `group_levels` も同じ。上書きの無い行列と同じ `matrix_id` になる）。違えば `needs_input`（`field` = `matrix_id`、
    `patch: true`。選択肢は「行列の群で比べる」= `{"group_levels": null}`（渡した群の引数を外す）と「渡した群で比べる」=
    `{"matrix_id": null}`（群が同じ直近の `dataset_preprocess` の行列か既定の前処理。行列の前処理の設定を保つなら、`dataset_preprocess`
    にその設定と群を渡して作った `matrix_id` を渡す）。群は行列が決めている）。選択肢は対比の群の名前がその群にある読みだけで、行列の
    群でしか比べられなければ群の引数を使わず `notes` に書いて比べ、どちらでも比べられなければ群の誤り。
  - `comparisons` の無い `factors` も誤りにしない: 空（`[]`）か `["group"]` なら群の指定と同じ設計なので `notes` に書いて比べる。
    ほかの因子を挙げ、群の名前がその因子の水準で `factors` のモデルが組める（比較に入る試料に因子の値が欠けず、2 因子以上なら水準の
    直積のどの cell にも試料がある。組めないモデルは `DESIGN_INVALID`）なら「その因子の水準として比べる（ほかの因子の水準ごとに
    `at` を付ける）」と読め、群の名前が group 列の水準でもあれば「`factors` を外して group 列で比べる」とも読めて数字が変わるので
    `needs_input`（`field` = `comparisons`、`patch: true`。選択肢は因子ごとの
    `{"comparisons": [...], "reference_group": null, "test_group": null}`〔`groups` で指定したなら `"groups": null`〕と、group 列の
    読みが成り立つときだけ `{"factors": null}`。対比が 50 個を超える読みは選択肢にせず `reason` に書く）。群の名前がどの因子の水準でもないか、モデルが組めなければ group 列で比べて
    `notes` に理由を書く。
  - 途中で止まった実行（探索専用）のデータセットも断らない（書き出しと図に探索専用の印。§11.4）。
  - 既定の前処理は尋ねない（`drift_mode=auto` が補正を当てなかったら `warnings` に `drift_not_applied_by_default`）。`matrix_id` を
    省いた呼び出しが行列を決めた後にエラーになったら、`details.matrix_id` と `matrix_source` にその行列と出所が残る（`default` は
    既定の前処理で作って登録した行列で、呼び直しに使える。`latest_preprocess` は直近の `dataset_preprocess` の行列で、群・除外・
    設定はその前処理のもの＝`details.matrix_note`）。
- `record`: `{written, record_id, path}`（比較の記録。親は行列の前処理の記録。設定に `spec`・`log2fc_threshold`・`matrix_source`）。
  既定の前処理で行列を作ったときは `preprocess_record` に前処理の記録も付く。

#### 11.2.1 `dataset_plot(kind="volcano")` の返り値（`lipidmix.volcano.v1`）

`dataset_compare`（2 群）の後に呼ぶ read-only の図。画像には件数入りの 1 行のキャプションが付く（表は payload の契約。§11.4.9）。
PNG ファイルは `output_path`（任意の場所）か `save_figure(kind="volcano")`（`reports/figures/`）。2 群の対比が複数ある結果で `contrast_id`（か `result_id#contrast_id`）を省くと
`needs_input`（選択肢は対比の一覧。`core.md`）。題とキャプションは `test vs reference`（`group_b` vs `group_a`。log2fc は正なら
test が高い）。

画像モードは間引かず全特徴を描く。**有意件数はキャプションの `up=` / `down=` / `ns=` で読む**（図の点を数えない）。

| フィールド | 意味 |
|------------|------|
| `points[].feature` | 特徴量名（`differential` の `feature`） |
| `points[].log2fc` | log2 fold change（x 軸）。**正 = 群 B で高い** |
| `points[].neg_log10_p` | `-log10(p)`（y 軸。`q` ではない）。p が 0.0（アンダーフロー。検定はできている）は最小の正の float の `-log10`（約 307.65。検定不能にしない。結果の p は 0 のまま） |
| `points[].sig` | `up` / `down` / `ns`。書き出しの `significant` と同じ規則（§11.4。判定 `changed` で q と \|log2fc\| が閾値を満たす）の行だけ符号で `up` / `down`、ほかは `ns`（`inconclusive` は q が小さくても `ns`）。q と判定は比較の補正（`p_adjust`）のままなので、`bonferroni` / `holm` / `none` なら書き出しの BH の q と違う |
| `thresholds` | 判定に使った `q` と `log2fc` のしきい値 |
| `comparison` | `group_a` / `group_b` / `n_a` / `n_b`（実 n） |
| `render_hints.guides.x` | 縦破線の位置（`±log2fc_threshold`） |
| `render_hints.guides.y` | 横破線の位置（`-log10(q_threshold)`）。y 軸は `-log10(p)` なので目安で、有意判定とは厳密に一致しない |
| `selection.total` | 元の点数（＝検定対象の特徴量数） |
| `selection.plotted` | 返した点数 |
| `selection.significant_total` / `significant_plotted` | 有意点の総数と返した数。**常に一致する**（有意点は間引かない） |
| `selection.ns_total` / `ns_plotted` | `ns` 点の総数と返した数。間引くと乖離する |
| `selection.max_points` | 上限点数（既定 800。`up` / `down` は上限に関わらず全件残る） |
| `selection.dropped_nonfinite` | `log2fc` か `p` が有限でなく描かなかった件数。**「有意でない」ではない**（分散 0・欠損・片群のみ検出などで検定不能） |

`selection.plotted < selection.total` なら図は全点ではない。有意件数は `selection.significant_total` で読む。間引きは `ns` 点の
等間隔サンプリングで決定的（同じ入力なら同じ点）。

#### 11.3 必須caveat

`dataset_compare` が戻り値に出す（黙って捨てない）注意。解釈と報告の注意点として引用し、交絡・小 n の下での因果的な解釈は保留する:

1. **交絡（群⟂バッチ）**: `confounding` と対比の `confounded`（§11.2）。判定は実際に比べた群（QC・ブランクを除く）で行う。
2. **小 n・検定不能**: `warnings` と判定 `inconclusive`（`verdict_reason`）。`inconclusive` を `unchanged` と読まない。
3. **正規化の状態**: 前処理の `steps[]`（§10.2・§10.3）。正規化しなかった行列の log2FC は測定量の差を含みうる。

#### 11.4 `dataset_export_differential` — 差次的結果のエクスポート契約

`dataset_compare` の後に呼ぶ（`result_id` を省くと直近の結果。対比が複数の結果は `contrast_id` が要り、`result_id#contrast_id` も
受ける。省くと `needs_input`〔`field` / `argument` = `contrast_id`、選択肢は対比の一覧〕で、読める対比が 1 つならそれを使う）。

特徴ごとの同定（`feature_metadata` の名前・InChIKey・m/z・RT）を比較の結果と特徴 ID で結んで 1 ファイルに書く。読み込みの経路
（解析フォルダ・mzTab-M）によらず同じ規則で、行の出所は InChIKey を供給した側: SME 行の InChIKey があれば `name_source = mztab_sme`、
無くて SML の候補が一意なら `mztab_sml`（MS1 照合のみ）、キュレーションの `assign` が当たれば `curation`。**InChIKey の無い特徴は
書き出さない**（下流で「パスウェイが無い化合物」と区別できないため。数はメタ行 `n_unannotated`）。InChIKey 付きの行が 0 件なら
`NO_ANNOTATED_FEATURES` で止まる（本文 0 行のファイルを残さない）。**有意な行だけでなく、InChIKey 付きの全行を書く**（下流の
濃縮解析が背景に使う）。

`#` 始まりのメタ行に来歴（`contract_version` / `group_a` / `group_b` / `log2fc_sign` / 閾値 / `n_features_total` /
`n_with_inchikey` / `n_unannotated`）、続けて TSV 本体。列は `spot_id / name / name_source / ontology / inchikey / inchikey_source /
msi_level / mz / rt / log2fc / p_value / q_value / mean_a / mean_b / significant`（契約の版は `CONTRACT_VERSION`）。

統計エンジンの結果（`dataset_compare`）は `msi_note` の後にメタ行を足す（既存のキーと重ならない。受け手は未知のキーを許す）:

| メタ行 | 意味 |
|---|---|
| `# test = … effect_size_definition = … ci_level = …` | 実際の検定名（`moderated_t` など）・効果量の定義の名前・CI の水準 |
| `# verdict_q_threshold = … equivalence_fc = … bh_population = … bh_n = … q_interpretation = … p_adjust_method = …` | 判定の規則と BH の母集団（比較の補正の母集団と同じ。`bh_n` は選んだ対比の大きさ）。`p_adjust_method` は**比較の**補正（`bh` / `bonferroni` / `holm` / `none`。比較の記録に無ければ `bh`）で、`q_value` 列の方法ではない。`q_interpretation` は `q_value` の読み方: `discovery_fdr`（BH・`matrix` / `annotated`）/ `reproduction_only`（母集団が `reported`） |
| `# q_value_method = bh  q_value_source = …  q_value_note = …` | `q_value` 列はいつも Benjamini–Hochberg の q（FDR）。`q_value_source` は `comparison`（比較が BH で、その q）か `recomputed_from_p_value`（比較が `bonferroni` / `holm` / `none` で、同じ母集団〔`bh_population` / `bh_n`〕の p から計算し直した）。計算し直したときだけ `q_value_note`（比較の方法と、`significant` もこの q で読み直したこと。ASCII） |
| `# unit = … n_units_a = … n_units_b = …` | 比較の単位。`subject`（全試料が被験者 ID を持ち技術反復を 1 個体にまとめた）・`injection`（被験者 ID が無い・まとめを切った）・`mixed`（混在）。`n_units_*` は比較に入った単位の数（`n_a` / `n_b` は注入の数のまま） |
| `# mean_definition = …` | `mean_a` / `mean_b` は個体ごとの前処理後の線形値の算術平均 |
| `# matrix_id = … rule_hash = … contrast_id = …` | どの行列・規則・対比を書いたか（2 群比較の結果の対比は 1 つ） |
| `# significant_rule = q<=q_threshold & \|log2fc\|>=log2fc_threshold & verdict=changed` | 15 列目 `significant` の規則（下記） |
| `# contrast_kind = interaction  interaction_factors = A x B  interaction_levels = A reference a0 test a1; B reference b0 test b1` | 交互作用の対比のときだけ。どの 2 因子のどの水準の差の差か（因子ごとに基準 `reference` と比較 `test` の水準。大小の記号は使わない） |
| `# interaction_definition = …` | 交互作用の `log2fc` の定義（`(a1 vs a0 at B b1) minus (a1 vs a0 at B b0)`）。`mean_a` / `mean_b` は空欄で、`mean_definition` もその旨 |
| `# exploratory = true  exploratory_reason = …` | 探索専用のデータセット（未完了の実行を `allow_incomplete=True` で開いた・検証を通っていない pipeline の単位）のときだけ。理由は 1 行（`=` とタブは含まない） |
| `# n_excluded_after_test = …  excluded_after_test_note = …` | `exclude_feature_ids` を渡したときだけ。検定の後に落とした行の数（q は計算し直さず、多重比較の母集団は検定のときのまま）。`n_features_total` は落とした行も数え、`n_with_inchikey + n_unannotated + n_excluded_after_test` が `n_features_total`（注記にも書く） |

`p_value` はエンジンの p、`q_value` は **BH の q**（上の `q_value_*`）、`log2fc` は正なら `group_b`（比較）が高い。**`significant` は `q ≤ q_threshold`・`|log2fc| ≥ log2fc_threshold`（`dataset_compare` の値）・判定
`verdict` が `changed` のすべてを満たす行だけ `true`**（`significant_rule`）。量を判定できない行（`inconclusive` /
`quantity_not_judgeable`）と BH の母集団の外の行は、q が小さく |log2fc| が大きくても `false`（gap-fill の値で動いた量を主張にしない）。
`verdict` は `q_value` の BH の q で読んだ判定（比較が BH なら比較の判定そのもの、別の補正なら同じ規則〔`compare._verdict`〕で BH の q
から読み直したもの。戻り値の `notes` にも書く）。比較が BH なら `true` の行の数はその対比の `large_changes` と同じ（InChIKey の無い行の
分だけ少ないことがある。`holm` / `bonferroni` / `none` なら違いうる）。同じ印を根拠バンドルの `features.tsv` の `significant` と
volcano の `up` / `down`（§11.2.1）が使う（volcano は比較の補正の q で読むので、比較が BH でなければ件数が違いうる）。**`verdict` の
3 値は TSV の列に出ない**ので、判定は `dataset_compare` の要約・`top` か根拠バンドル（`reanalysis.md` §10）で読む。

**交絡を越えた比較の印**: `allow_confounded=true`（既定）で完全交絡（`confounding.batch_status = complete` か、書き出す対比の
`confounding.by_contrast[対比].confounded = true`）を越えて走った比較は、`source_verification` の後に `# unadjusted_confounded = true`
と `# unadjusted_confounded_note = …`（日本語の説明）を書く。volcano の `caveats` と図の下端（ASCII の但し書き）、根拠バンドルの
`comparison.unadjusted_confounded`（`reanalysis.md` §10）にも同じ印。部分的・評価できない・無い交絡には何も足さない。

**交互作用の対比**（2×2 の 2 因子の設計の `interaction_*`）は名指しすれば書き出せる。`group_a` は基準の 2 セル（`a0:b1 + a1:b0`）、
`group_b` は比較の 2 セル（`a1:b1 + a0:b0`）、`log2fc` は差の差（正なら `a1` の効果が `b1` で大きい）、`mean_a` / `mean_b` は空欄。
根拠バンドルと volcano は 2 群の対比だけ（交互作用は `CONTRAST_NOT_TWO_GROUP`）。`contract_version = 1` と `log2fc_sign` の行は 2 群と
同じなので、交互作用の印はメタ行 `contrast_kind = interaction` と `interaction_definition`（戻り値の `contrast_kind` にも）だけ。受け手は
この印を見て、`log2fc` を 2 群の比、`group_a` / `group_b` を 2 つの群と読まない（2 群の書き出しに `contrast_kind` の行は無い）。

**探索専用のデータセットも断らない**（未完了の実行・検証を通っていない pipeline の単位。`exploratory_only`）: `msi_note` の後に
`# exploratory = true` と `exploratory_reason` を足す（戻り値の `exploratory` にも）。欠けた検体が「その群には無い」ように見えうるので、
受け手は結論の重みを下げる。根拠バンドルでは manifest の最上位の `exploratory` と `comparison.exploratory`（`reanalysis.md` §10）。

**検定の後に行を落とす**（`exclude_feature_ids`）: その特徴 ID の行を書かない。q は計算し直さない（多重比較の母集団は検定のときのまま。
母集団が `annotated` / `reported` なら、落とした行がもともと母集団の外のこともある）。落とした数と注記はメタ行 `n_excluded_after_test`。
`n_features_total` は落とした行も数える（`n_with_inchikey + n_unannotated + n_excluded_after_test`。注記にも書く）。見つからない ID は
戻り値の `exclude_feature_ids_not_found`。

**`ontology` と `msi_level` はいつも空欄**（キュレーションの `assign` が当たった行の `ontology` だけは記録した値）。メタ行 `msi_note` が
宣言し、空欄は「該当なし」でなく「この書き出しでは取っていない」。解析フォルダのデータセットの `msi_level`（`.arf2` 由来の保守的な
クラス上限。MS/MS の有無ではない）は `dataset_feature_search` に出るが、書き出しには載せない。`spot_id` は特徴 ID で、メタ行は
`# id_space = mztab_smf_id`・`# source_mztab`（読んだファイル。キーの名前は経路によらず同じ）。

#### 11.4.9 データセットの上の同じ図（`dataset_plot`。B-2b Task 7）

`dataset_plot(kind="group_intensity" | "species" | "species_pca")` は §11.5〜11.7 の payload 契約（`lipidmix.group_intensity.v1`・
`lipidmix.species_intensity.v1`・`lipidmix.species_pca.v1`）を `load_dataset` のハンドル（`DatasetState`）から作る。`items` は特徴 ID も
受ける。`value` は `"share"` / `"intensity"`（`"intensity"` の payload の `value` は `"height"`）。`ncols` は無い。値は `feature_matrix`、
gap-filled の印は `detected_mask`、群と役割は試料の対応表、MS/MS の裏付けは `.arf` のセルの MS/MS の有無（`.arf2` の `has_msms` では
ない）と `.arf2` の Name の接頭辞、人の判断は `load_dataset` が読んだ `.arf2` のキュレーション、除外は `dataset_exclude`（試料は描かず、
特徴はクラスにも名前にも入れない）。`species_pca` の payload は `points`（`fitted` が False の点は投影だけ）・`explained_variance_ratio`・
`pc1_r_top` / `pc1_r_bottom`・`result_id` などで、ローディングの全行は `ds.results[result_id]`。図は `ds.figures[kind]` に置かれ、
`save_figure(kind=...)` が dpi 300 の PNG と SVG を書く（試料の情報の更新後は作り直しを求められる）。

**共通の形**: 既定（`output="image"`）はサーバ側で描いた PNG と説明。§11.2.1・§11.5〜11.8 の表は `output="payload"`（または env
`LIPIDMIX_PLOT_OUTPUT=payload`）の契約。§11.5〜11.8 の入力の誤りは誤りの封筒
`{"error": {"code": "DATASET_BAD_REQUEST", "message", "details"?}}`（例外ではない）。

**図の保存の誤り**（`dataset_plot` のどの kind でも同じ）:

- `dataset_plot(output_path=...)` に書けなければ `FIGURE_SAVE_FAILED`（`details.output_path`。書く段だけの符号で、図を組み立てる間の
  入力の読み取りの失敗には使わない）。セッションの図（`ds.figures[kind]`）は変えない。
- `save_figure` は置かれた図が古ければ保存せず `STALE_FIGURE`: 試料の情報を焼き込んだ図の後に試料の情報が更新された、または図の元の
  結果がいまの状態から出たものでない（`details` に `kind`・`result_id`・`required_tools: ["dataset_plot"]`。描き直してから保存する）。
  除外（`dataset_exclude`）・試料の表の更新（前処理の入力の列なら pca / loadings / species_pca / volcano、群・足した因子の列だけなら
  volcano）が捨てた図も `STALE_FIGURE`（`details` に `kind`・`reason`〔実際に変わったもの: `exclusion` や試料の表の列名。因子の列
  だけが変わったときは、差次的解析を群と同じく古くしても `group` を挙げずその列だけ〕・`required_tools: ["dataset_plot"]`）。一度も
  描いていない図だけが `missing_state`（`required_tools` は `dataset_plot`）。
- `save_figure` は kind が使わない `result_id` / `contrast_id` を断らずに使わず、戻り値の最後の `notes:` 行に書く（`species_pca` の
  `result_id` は置き場の図の結果と照らし、違えば `STALE_FIGURE`〔`details` に `held_result_id`・`mismatch`〕。§11.7）。
- `dataset_plot` も別の kind の引数は断らずに使わず（図はその引数を渡さないときと同じ）、どの kind の引数かを説明の最後の `notes:` 行に
  書く（payload は最上位の `notes` 欄で、形式名は変えない。描けなかった誤りの封筒は `details.notes`: 引数の形の誤り〔`output` /
  `output_path`〕・ハンドルの誤り・`missing_state`〔データセットや結果が無い。`details` はこのときだけ足す〕・`FIGURE_SAVE_FAILED`・
  図の組み立ての誤りのどれも。kind が未知の誤りと、誤りでない `needs_input`〔対比を選ばせる問い〕は持たない）。

#### 11.5 `dataset_plot(kind="group_intensity")` の返り値（`lipidmix.group_intensity.v1`）

選んだクラス・分子種が、試料群ごとに**見つかるか・どれくらいあるか**を示す。解析フォルダのデータセットが対象（`.arf` と同じ
アラインメントの兄弟 `.arf2` が要る）。画像の 1 行のキャプションは項目ごとのスポット数（MS/MS の件数が分かれば括弧内に。0 も出す）・
群ごとの試料数・除外の 5 種類の件数（手動除外・内部標準・判断・自動判定・標準液）・N.D. の項目名・`caveats`。ファイルは
`save_figure(kind="group_intensity")`（`reports/figures/<analysis_id>_group_intensity.png` を dpi 300 で同名の `.svg` と一緒に書く。
直前の `dataset_plot(kind="group_intensity")` が無ければ `missing_state`）。

**エラーと状態**: 入力の誤り（群が 1 試料にも当たらない・`output` の値が不正・`detection_limit` が 0 以下や非有限など）は
`DATASET_BAD_REQUEST` の封筒（§11.4.9）。失敗した呼び出し（`missing_state` を含む）は**セッションの図
（`ds.figures["group_intensity"]`）を変えない**（`save_figure(kind="group_intensity")` は前回の図を保存する）。図を置くのは、組み立て・
（画像モードなら）描画・`output_path` の書き出しがすべて成功してから。群の指定が `blank` / `qc` の役割で全滅したら、`blank` / `qc` を
書くよう案内を添える。

**検定はしない**。1 パネル = 1 項目、1 点 = 1 試料、群ごとに log10 空間の平均 ± SD を添えるだけ（検定は `dataset_compare`）。

| フィールド | 意味 |
|------------|------|
| `plot_schema` | `lipidmix.group_intensity.v1` |
| `value` | 点の値の定義: `sum of PeakHeight (log10 on the plot)`。1 点 = その試料で項目に当たったスポットの **PeakHeight の合計**（`.arf` の試料別行。gap-fill の値も含める） |
| `detection_limit` | 検出下限（図の破線）。無ければ null |
| `detection_limit_source` | `argument`（`detection_limit` 引数）/ `param_file`（`.arf2` と同じフォルダの param ファイルの `Minimum peak height`）/ null（下限なし）。引数が優先 |
| `groups[]` | `label`（指定した群の文字列）と `samples`（解決した試料名）。指定順。`dataset_exclude` した試料は含まない |
| `low_reliability_samples` | 白抜きで描き、平均 ± SD から外した試料（抽出量が少なかった等）。解決後の試料名 |
| `items[].item` | 指定した項目の文字列（指定順。1〜30 件） |
| `items[].parts[]` | `+` で分けた各部分。`part`・`kind`（`class` = `.arf2` の Ontology に完全一致、`name` = 分子種名に完全一致、`none` = 当たらず。名前が当たっても全スポットが除外で消えたときも `none` で、そのスポットは `excluded` に載る）・`n_spots` |
| `items[].spots[]` | 当たったスポット `spot_id` / `name` / `ontology`（キュレーションの `assign` は付け替え後の値） |
| `items[].n_spots` | 当たったスポット数（付加体違い・重複スポットも合計に入る） |
| `items[].n_spots_msms` | 当たったスポットのうち MS/MS の裏付けがある数。**スポットが 1 つ以上あって 0 のパネルだけ灰色で「MS1-only (unconfirmed)」**（MS1 の一致だけの同定。スポットが 0 の N.D. パネルは灰色にしない） |
| `items[].detected` | どこかの試料で値が 0 より大きいか |
| `items[].groups[].samples[]` | `sample` / `value`（PeakHeight 合計）/ `gap_filled_fraction`（合計のうち gap-fill の値の割合。値 0 なら null）/ `low_reliability` |
| `items[].groups[].log10_mean` / `log10_sd` / `n_in_stats` | log10 空間の平均・SD と、その計算に使った試料数（値 > 0 かつ低信頼でない試料）。n = 1 なら SD は null、0 なら平均も null |
| `excluded` | 除いたスポットの `#<spot_id> <name>`。理由ごとに 5 種類（`manual` / `internal_standard` / `curation` / `auto_likely_wrong` / `standard_only`。下記） |
| `caveats` | 試料が複数の群に当たった（群の組ごとに 1 件。重なった試料数と最大 3 試料名）・同じ群を重複指定した・孤立した判断記録がある・レビューが無く自動判定を除けなかった、等 |

**N.D. は 2 種類**。(a) 項目に当たるスポットが無い（`n_spots` = 0。名前が当たらなければ `parts[].kind` は `none`、クラスに当たっても
除外で全スポットが消えたら `class` のまま 0）。(b) スポットは当たったが、描いた試料のどれでも値が 0（`n_spots` > 0 で `detected` が
false）。どちらも N.D. パネルとして残し、エラーにしない。(b) は「存在しない」ではなく**この群の試料では検出されなかった**。

**点の読み方**: 縦軸は log10(PeakHeight) で**全パネル共通**（パネルごとに拡大するとノイズが信号に見える）。値 0 の点は軸の底に ▽。
`gap_filled_fraction` が **0.5 を超える点は ◆**（合計の過半が gap-fill の値で、検出の証拠として弱い）。

**`excluded` の 5 種類**:

| キー | 除くもの |
|------|----------|
| `manual` | `dataset_exclude` でスポットとして除いたもの（セッションの `excluded_spots`）。クラスにも名前指定の部分にも効く |
| `internal_standard` | 名前に `(d<数字>)`（`PC 33:1(d7)` 等）を持つ標識内部標準。**クラスとして当たった分だけ**除く（名前で指定した部分は描く） |
| `curation` | `apply_curation=True`（既定）のとき、人の判断が `wrong` / `redundant` のスポット（`clear` 済みは無効）。`assign` は記録した名前・クラスで扱う。`suspect` は描く。名前指定の部分にも効く |
| `auto_likely_wrong` | `exclude_auto_likely_wrong=True` のとき、このアラインメント（`.arf2` の sha256 一致）の最新 `curation_review` で `likely_wrong` のスポット。名前指定の部分にも効く。レビューが無ければ除かず `caveats` に出す |
| `standard_only` | `standard_samples` を指定したとき、その試料での最大高さが描く群の試料での最大高さの 10 倍以上のスポット（クラスとして当たった分だけ。標準液の奇数鎖標準物質がクラス合計に混ざるのを防ぐ）。判定には `dataset_exclude` 済みの試料の行も使う |

図を報告に載せるときは `excluded` の件数と理由を併記する。

#### 11.6 `dataset_plot(kind="species")` の返り値（`lipidmix.species_intensity.v1`）

選んだクラスを**分子種（アラインメントスポット）に展開**し、1 分子種 1 パネルで群ごとの試料別の値（割合 % または PeakHeight）を
並べる。解析フォルダのデータセットが対象（`.arf` と兄弟 `.arf2` が要る）。`group_intensity`（クラス合計）と違い**スポットごと**に描く
（同じ名前でも付加イオンが違えば別パネル）。画像の 1 行のキャプションはパネル数・群の試料数・除外の 6 種類の件数・分母が 0 で
描かなかった試料・`caveats`。ファイルは `save_figure(kind="species")`（`reports/figures/<analysis_id>_species.png` を dpi 300 で同名の
`.svg` と一緒に書く）。入力の誤り・パネルが 40 枚を超える・判断の記録ファイルが読めないときは `DATASET_BAD_REQUEST` の封筒。失敗した
呼び出しはセッションの図を変えない（§11.5）。**検定はしない**。

| フィールド | 意味 |
|------------|------|
| `plot_schema` | `lipidmix.species_intensity.v1` |
| `value` | `share`（試料ごとに 高さ ÷ 分母の合計 × 100）または `height`（PeakHeight） |
| `stat` | 平均と SD の定義: `mean ± SD of share (%)` / `mean ± SD of log10(peak height)`。`groups[].mean` / `sd` はこれに従う |
| `share_basis` | `value="share"` のときだけ。`items`（指定した分母の項目。省略時は null で、描く項目が分母）と `spot_ids`（分母にしたスポット）。`height` では null |
| `groups[]` | `label` と `samples`（解決した試料名。指定順。`dataset_exclude` した試料は含まない） |
| `low_reliability_samples` | 白抜きで描き、平均・SD から外した試料 |
| `zero_denominator_samples` | 分母（分母スポットの高さの合計）が 0 で割合を計算できなかった試料。その試料の `share` は null で描かない |
| `spots[]` | パネルごとに `spot_id` / `name` / `label`（表示名）/ `ontology` / `item`（展開元の項目）/ `adduct` / `mz` / `rt` / `msms`（MS/MS の裏付けの有無）と `groups[]` |
| `spots[].groups[].samples[]` | `sample` / `height`（PeakHeight。行が無ければ 0）/ `share`（分母に対する %。`height` 指定時と分母 0 の試料は null）/ `gap_filled` / `low_reliability` |
| `spots[].groups[].mean` / `sd` / `n_in_stats` | `stat` の定義での平均・SD と、その計算に使った試料数（低信頼でなく、`share` なら計算できた試料、`height` なら高さ > 0 の試料）。n = 1 なら SD は null、0 なら平均も null |
| `excluded` | 除いたスポットの `#<spot_id> <name>`。理由ごとに `internal_standard` / `manual` / `curation` / `auto_likely_wrong` / `standard_only`（§11.5 と同じ）と、`require_msms=True` で MS/MS の裏付けがなく除いた `no_msms` |
| `caveats` | `group_intensity` と同じ種類に加え、`share_basis` が描く分子種を含まない（割合が 100 % を超えうる）・`items` に当たらなかった（またはすべて除外された）部品がある・`value="height"` なのに `share_basis` を渡した（無視される） |

**割合の読み方**: `share` は**分母にしたスポットの高さの合計に対する %** で、試料の総量に対する割合ではない。分母は `share_basis`
（省略時は `items`）を同じ除外規則で解いたスポットで、`require_msms=True` は分母にも効く。図を報告に載せるときは `share_basis` を
併記する。`share_basis` の部品が 1 つでも当たらない（綴り違いなど）と図を描かず `{"status": "error"}`（当たらなかった部品名入り）。
`items` の当たらない部品は描いたうえで `caveats` と説明文に書く。

#### 11.7 `dataset_plot(kind="species_pca")` の返り値

選んだ分子種（`species` と同じ項目指定・除外）の PeakHeight 行列（試料 × 分子種）で PCA を回す read-only の図。解析フォルダの
データセットが対象（`.arf` と兄弟 `.arf2` が要る）。画像はスコア図で、説明に寄与率・`pc1_r_top` / `pc1_r_bottom`。ファイルは
`save_figure(kind="species_pca")`。ローディング全量（`loadings_rows`）は戻り値に載せずセッションに残す。入力の誤り（計算に使う試料が
3 未満・合計 0 の試料・未知の `orient_by` / `normalize` など）は `DATASET_BAD_REQUEST` の封筒。失敗した呼び出しは前回の図を変えず
（§11.5）、入力の誤りなら `ds.results` にも何も足さない（結果を登録するのは計算が通ったときだけ）。

| フィールド | 意味 |
|------------|------|
| `points[]` | 試料ごとに `x`（PC1）/ `y`（PC2）/ `label`（試料名）/ `group`（群。2 群に入る試料は最初の群）/ `fitted`。`fitted=false` は低信頼の試料で、**主成分の計算に使わず投影しただけ** |
| `explained_variance_ratio` | 主成分ごとの寄与率（計算に使った試料での値。最大 5 成分） |
| `x_label` / `y_label` | 軸名（寄与率入り） |
| `n_fit` / `n_species` | 計算に使った試料数 / PCA に使った分子種数 |
| `dropped_zero_variance` | 計算に使う試料で値が変わらず外した分子種の数 |
| `excluded` / `caveats` | `species` と同じ（§11.6）。画像の説明文にも除外の件数と caveat を載せる |
| `pc1_r_top` / `pc1_r_bottom` | PC1 との相関 r が正の上位・負の上位の分子種（表示名と r。各最大 5 件。正が無ければ `top` は空）。r は**計算に使った試料**での分子種と主成分スコアの相関で −1〜1 |
| `result_id` | この結果の ID。`save_figure(kind="loadings", result_id=...)` で分子種のローディングを描ける。`save_figure(kind="species_pca", result_id=...)` は描き直さず、置き場の図がこの結果から描いたもの（か同じ指定で描き直しただけの結果。そのときは `notes:` 行に書く）なら保存し、別の結果の図なら保存せず `STALE_FIGURE`（`details.held_result_id` に置き場の図の結果、`details.mismatch` に違いの種類）。`mismatch` は `settings`（項目の指定が違う。`details.settings`〔この結果の項目の指定〕で描き直してから保存する）/ `stale`（この結果が除外・試料の表の更新より前）/ `input`（同じ指定だが、判断の更新などで入力の行列が違う）。`stale` と `input` はこの結果の図をもう描けない（描き直すと今の状態の新しい結果になる）ので、置き場の図を保存するなら `held_result_id` を渡すか省く。`details.settings` と `required_tools: ["dataset_plot"]` は、この結果の指定が置き場の図と違うときだけ付く |

**読み方**: 値は `normalize`（`none` / `total`）→ log10(x + 1)（`log_transform`）→ 計算に使う試料の平均と母標準偏差で autoscale の順。
主成分の**符号は任意**で、`orient_by` を指定したときだけ、その群の平均スコアが正になるようにそろえる。`normalize="total"` は試料量の
差を除いて組成を比べる選択で、合計が 0 の試料があると計算しない。

#### 11.8 `dataset_plot(kind="loadings")` の返り値（`lipidmix.pca_loadings.v1`）

保存済みの PCA（`dataset_pca` と `dataset_plot(kind="species_pca")`）のローディングを、主成分ごとの横棒で描く read-only の図。
ファイルは `save_figure(kind="loadings")`（PNG と同名 `.svg`）。入力元は `result_id`（省けば直近の `dataset_pca`。`species_pca` の結果も
渡せる）。候補（結果の出どころ）が 2 つ以上あってどれも優先できなければ `AMBIGUOUS_RESULT_SOURCE`。ローディングを持たない古い PCA
結果は候補にならない（PCA をやり直す）。入力の誤り（`pcs` が 1〜3 個でない・その PCA にない主成分・未知の `value`・全件表示が
60 特徴量を超える、など）は `DATASET_BAD_REQUEST` の封筒。失敗した呼び出しは前回の図を変えない（§11.5）。

| フィールド | 意味 |
|------------|------|
| `result_id` | どの PCA 結果から描いたか（`save_figure` の `result_id` と同じ名指し） |
| `value` | **実際に描いた値の種類**（`r` / `coefficient`）。`value_requested` と違えば `caveats` に理由 |
| `value_requested` | 求めた値の種類 |
| `top_n` | 主成分ごとに取った正・負それぞれの上位件数。`null` は全件 |
| `aligned` | `true` は全件表示で、全パネルの行が最初の主成分の順にそろう。`false` は主成分ごとに上位を独立に選んだ |
| `panels[]` | 主成分ごとに `pc`（`PC1` …）/ `explained_pct`（寄与率 %）/ `rows[]` |
| `rows[]` | `feature_id` / `label`（表示名）/ `ontology` / `adduct` / `mz` / `rt` / `coefficient`（固有ベクトルの成分）/ `r` / `value`（描いた値） |
| `caveats` | 値を `coefficient` に落とした理由など |

**読み方**: `r` は特徴量とその主成分スコアの相関（−1〜1。成分 × 特異値 / √n）で、autoscale した PCA でだけ出せる（そうでなければ
`coefficient` で描く）。主成分の**符号は任意**で、同じ主成分の中の向きの違いだけを読む。ARF のローディングは、その PCA に使った行列
（`dataset_pca` の前処理後の行列）の全スポットのうち**上位だけ**（`top_n`）。**単一の props の PCA** だけに対応する（スポットあたり
複数の列を並べた PCA は `{"status": "error"}`。props を 1 つにして PCA をやり直す）。棒の色は `.arf2` の Ontology（クラス）。

#### 11.9 最初の概観（`dataset_overview` と `load_dataset` の `overview`。構成見直し A）

読み込んだ直後に、群の間で何が変わり何が変わらなかったかを **品質 → 全体 → クラス → 分子種** の 4 段で粗く示す。
`dataset_overview(top_n=5, dataset_handle=None)` が返し、鎖を進めた後（除外・試料の表の訂正・キュレーション・定量の鎖）に呼び直すと
出し直す。`top_n`（分子種の段の上位の件数）の上限は 20（超える値は断らず 20 に下げ `notes` に書く。整数でない・負は
`DATASET_BAD_REQUEST`）。`load_dataset(overview=True)` は読み込みの最後に同じものを `overview` 欄に入れる（既定 `False` では載せない。
失敗しても読み込みは成功し `{"status": "unavailable", "error": {...}, "next": ...}`）。

**戻り値の形**（`load_dataset` の中では `counts` を省く。`summary` と重なるため）:

| キー | 意味 |
|---|---|
| `overview_id` | `ovw_<12 桁>`。データの指紋・試料の表の版・行列・比較と PCA の結果・設定から決まる（同じ入力なら同じ ID） |
| `assay_kind` | そのハンドルの種別（`lipid` / `metabolite` / `unknown`）。クラスの段と組成のラベルの規則を決める |
| `provisional` | 定量の鎖（付加体の選別 `dataset_select_features`）を通る前なら `true`。分子種の数に同じ分子種の複数の付加体・内部標準の行が混ざりうる |
| `state` | 1 行の状態（キュレーション・正規化・統合・選別の有無、前処理のレシピ、`assay_kind`）。ドリフト補正を飛ばしたら `drift skipped` |
| `counts` | 役割ごとの試料数・特徴数・同定の水準（`msms` / `ms1` / `none`）の内訳（`dataset_overview` だけ） |
| `quality` / `global` / `classes` / `species` | 各段。`facts`（事実の 1 行ずつ。先頭が ID）・`next`（掘る一手。引数は JSON）・省いた段は `omitted`（理由）。品質は `verdict` も |
| `result_ids` | 概観が登録した `compare` / `pca` の結果 ID（回したものだけ） |
| `record` | `reports/stats/` への記録（`written` と `record_id`・`path`、書けなければ `reason`） |
| `notes` | `top_n` を 20 に下げたときだけ。そのとき分子種の上位の表を 20 行で切ったら `species` に `top_truncated: true` と `top_n_more`（載せなかった行数）。全行の引き方を `notes` に書く: 2 群は `result_ids.compare` を `dataset_feature_search`（名指し）・`dataset_export_differential(result_id=...)` で、多群（上位は moderated F の `f_q` の順）は概観の群と行列を名指しした `dataset_compare(groups=[<概観の群>], matrix_id=<概観の行列>, output_path=...)` が書く F 検定の表（`<stem>.ftests.tsv`）で引く（`matrix_id` を省くと直近の `dataset_preprocess` の行列で比べるので、F の母集団と `f_q` が概観と違いうる）。比較を回さなかった概観（群が 1 つ・完全交絡・品質の段だけ）は引き方を書かない |

**事実の ID**: 段ごとに `Q1, Q2, …` / `G1, …` / `C1, …` / `S1, …`。外からは `<overview_id>#Q1`。解釈は事実と分けて書き、根拠の ID を挙げる。

**品質（Q）** — Q1 が判定。

| 判定 | 条件（閾値は記録の `settings.thresholds`） |
|---|---|
| `poor` | 群とバッチの完全交絡、失敗した QC 注入（総量が QC の中央値の 0.2 倍未満。QC 3 本以上のときだけ見る）がある、QC の RSD の中央値が 0.30 を超える |
| `caution` | 部分的な交絡、注入順の偏り、QC の総量（log2）と注入順の Spearman の強い相関（\|ρ\| ≥ 0.6 かつ p < 0.05、QC 5 本以上のバッチ）、技術反復の疑い |
| `not_evaluable` | QC が無い・足りず RSD を評価できない（理由 `qc_not_evaluable`。バッチの交絡も評価できなければ `no_qc_and_batch_unknown`）。**QC が無いときは `ok` と言わない**。`caution` / `poor` の判定にも `qc_not_evaluable` を添える |
| `ok` | 上のどれでもない |

`poor` でも差は出し、各段の事実に「品質 Q1 が poor のため参考」と付く（完全交絡のときだけ比較を回さない）。

- **陽陰を統合したデータセット**は、失敗した QC 注入とドリフトの兆候を**極性ごと**に見る（極性ごとの特徴の総量と、その極性の測定の
  注入順・バッチ。強度の大きい極性がもう一方の失敗やドリフトを隠さない）。失敗した注入はどれかの極性で失敗したもの（データの
  `by_polarity` に内訳）、ドリフトの兆候の行は `polarity` 付き。
- **既定（`drift_mode=auto`）がドリフト補正を当てなかった**（換算済み・注入順がファイルの順かもしれない）ときは、補正なしで計算し、
  そのことを事実と `next`（`dataset_preprocess(drift_mode="on")`）に出す。注入順が確かでない理由は出所どおりに書く（MS-DIAL が
  振った順なら「MS-DIAL の既定（ファイルの順）かもしれない」、シートの宣言 `# injection_order_source = filename` なら「シートの宣言で
  ファイル名の順」）。統合したデータセットでは当てなかった極性を名指しし、注入順が確かでない極性には「drift_mode="on" はその注入順の
  まま当てはめる」と実際の順の入れ方を添える。当てなかった極性のほかに補正の条件を満たさない極性があれば、`drift_mode="on"` はその
  極性のせいで `PREPROCESS_PRECONDITION` になるので `next` に勧めない。
- **注入順がファイルの順かもしれない**（MS-DIAL の既定・シートの宣言 `# injection_order_source = filename`）とき:
  - 交絡の事実「注入順が群で偏っている」に、偏りがファイル名の並びかもしれないと添え、データに
    `injection_order_may_be_file_order: true`、`next` に `dataset_set_sample_metadata`（実際の測定順を入れる）を足す。
  - ドリフトの兆候の事実（兆候あり・なしとも）は同じ注入順で ρ を測るので、ρ に使った QC の注入順がそうなら「ρ はファイル名の
    並びとの相関で、実際の測定順でのドリフトを見ていないかもしれない」と添え、データ（と `batches` のその行）に
    `injection_order_may_be_file_order: true`、`next` に実際の測定順の入れ方を足す。注入順が群で偏っていなくても、品質の段だけの
    概観（交絡を評価しない）でも付く。ρ と判定（`drift_sign`）は変えない。
  - 統合したデータセットでは極性を名指しし、入れ方を極性ごとに書く（統合したハンドルへのシートは陽イオン〔`positive`〕の測定しか
    決めず、統合し直すとそのハンドルのシートと除外は引き継がれない。§10.2 の「陽陰を統合したデータセット」・`quant.md` の
    `polarity_runs`）。陽イオン以外の極性は `next` が `<極性> のハンドルに dataset_set_sample_metadata → dataset_merge_polarities`
    （陽イオンの順も入れるなら陽イオンのハンドルにも当て、統合したハンドルに当てていた除外は統合し直してできたハンドルに当て直す。
    統合し直す手は `next` の中で `dataset_exclude` と `dataset_preprocess(drift_mode="on")` より前に置き、書いた順に行えば除外が新しい
    ハンドルに当たる）。陽イオンだけなら `next` は統合したハンドルへの `dataset_set_sample_metadata`。
- **既定の前処理が失敗する**（足切りで特徴が 1 つも残らない など）ときも概観は失敗しない: 品質の段（前処理の前の試料から計算する）
  だけを出し、全体・クラス・分子種の段は `omitted` に失敗の理由（足切りの内訳 `details` を含む）を書いて省く。交絡はこのとき評価しない。

**全体（G）** — 既定の前処理の行列で PCA（軸は生体の個体、QC は投影）。各軸の寄与率、群が各軸のばらつきをどれだけ説明するか（η²）、
群の重心。比較を回したら判定の数（2 群）か moderated F で差のある特徴の数（多群）。

**クラス（C）** — 組成（クラスの総量と CLR）。**付加体の選別の後だけ**出す（選別の前は同じ分子種の複数の付加体や内部標準の行が
クラスの合計に入り、見せかけの `unchanged` が出るので、`omitted` と鎖を進める `next`）。mzTab-M など解析フォルダ以外から開いた
データは定量の鎖が使えないので、`dataset_analysis_run` の `kind="composition"` を案内する。2 群は `class n_species log2fc_total ci_low
ci_high verdict_total clr_log2fc verdict_clr` の TSV（`changed` と `unchanged` を最大 5 つずつ）、多群は moderated F の q
（`class n_species f_q_total`）。

**分子種（S）** — `changed`（上昇・下降の内訳つき）/ `unchanged` / `inconclusive` の数。**`inconclusive` は「変わっていない」では
ない**（同等性の幅に収まった `unchanged` だけが変わっていない）。出方（検出の有無）だけに差があり強度の判定は `inconclusive` の特徴の
数も別に出す（検出の情報が無いデータでは `null`）。2 群は上位の TSV（`contrast_id feature_id name id_level log2fc ci_low ci_high q`）、
多群は F で差のある特徴の中の対ごとの数（`pair n_changed n_up n_down`）と上位（`feature_id name id_level f_q`）。

**群と向き**:
- 群が 1 つ以下・群とバッチの完全交絡では比較を回さない（`omitted` に理由）。完全交絡は全群をまとめた表で見る: 多群でまとめた表が
  部分交絡なら、一部の対だけが完全交絡でも比較は回し、その対は比較の記録（`result_ids.compare`）の `confounding.by_contrast` と行の
  `confounded` に印が付く。
- 探索専用のデータセット（`allow_incomplete=True` で開いた未完了の実行・検証を通っていない pipeline の単位）も比較を回す。比較の数字を
  使う段（全体・クラス・分子種）の**先頭の事実**が探索専用の印（理由つき。`data` は `{exploratory: true, reason}`）で、欠けた検体が
  「その群には無い」ように見えうるので差は参考。概観の記録の `settings.exploratory`（`{exploratory_only: true, reason}`）にも残る
  （探索専用でなければ欄が無い）。
- 2 群の向きは**仮**: 試料の表で先に出た群を対照にする（事実に「仮の向き」と付く）。向きを決めるのは
  `dataset_compare(reference_group=…, test_group=…)`。
- 3〜6 群は全対、7 群以上は対照との対だけ（と moderated F）。
- 判定は moderated t、q ≤ 0.05、同等性の幅 1.2 倍（`dataset_compare` の既定と同じ）。

**結果の置き場所**: 概観の PCA と比較は `ds.results` に入るが、直近の結果（`dataset_plot` / 書き出しの既定）には**しない**。概観の
比較の記録には `origin: "overview"` が付き、`dataset_feature_search` の既定（`result_ids` 省略）の選択から外れる
（`result_ids=["<result_id>#<contrast_id>"]` と名指しすれば引ける）。同じ比較を先に `dataset_compare` で回していれば、その記録に印は
付けない。概観そのものは `ds.results[overview_id]`（`kind: "overview"`）。記録は解析フォルダの `reports/stats/` に、既定の前処理と
概観の 2 件（概観の記録は前処理の記録を親に持つ）。書き先が無ければ書かずに成功する。

## 旧 `arf2` トピック（`.arf2` / スポット代表。§4, §8.2）

`.arf2`（スポット代表カタログ）のパーサ出力。`load_dataset` が `feature_metadata` に取り込み、`dataset_feature_search` が行で返す。

> Key 番号表（MS-DIAL の `[Key(N)]` から抽出した一次資料。リーダーのインデックス定数を変える前に確かめる）:
> `docs/schema/AlignmentSpotProperty.md`。

### 4. ARF2 (`metabolomix/arf2/reader.py`)

#### 4.1 `deserialize()` / `extract_arf2_data()`

型は `list[dict]`。1 要素は 1 スポットのカタログ情報で、サンプル別の強度は含まない。

| 列/キー | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | MS-DIAL のマスターアラインメント ID |
| `AlignmentID` | int | アラインメント ID |
| `RT` | float | スポット中心 RT（min） |
| `MassCenter` | float | スポット中心 m/z |
| `IonMode` | str | `Positive` / `Negative` / `Both` / `Unknown` |
| `Name` | str | 候補化合物名。欠けていれば `Unknown`。`no MS2:` や `low score:` は注釈の確度の接頭辞 |
| `HeightAverage` | float | サンプル間の平均ピーク高さ |
| `Formula` | str | 分子式 |
| `Ontology` | str | 化学 / 脂質クラス（`FA`, `PC`, `Cer_NS` など） |
| `SMILES` | str | SMILES |
| `InChIKey` | str | InChIKey |
| `AdductType` | str | 観測イオンの付加体（`[M-H]-` など） |
| `HeightMin` | float | サンプル間のピーク高さの最小値 |
| `HeightMax` | float | サンプル間のピーク高さの最大値 |
| `PeakWidthAverage` | float | サンプル間の平均ピーク幅（LC では通常 min） |
| `SignalToNoiseAve` | float | サンプル間の S/N の平均 |
| `SignalToNoiseMax` | float | サンプル間の S/N の最大値 |
| `SignalToNoiseMin` | float | サンプル間の S/N の最小値 |
| `MassMin` | float | サンプル間の観測 m/z の最小値 |
| `MassMax` | float | サンプル間の観測 m/z の最大値 |
| `FillPercentage` | float | スポットが値を持つサンプルの割合。**0..1 の比率**（百分率は 100 倍） |
| `MonoIsotopicPercentage` | float | 単同位体ピークとして扱われる割合。**0..1 の比率** |

#### 4.2 `summarize_arf2_data()`

| キー | 意味 |
|---|---|
| `total_spots` | スポット総数 |
| `height_average_median` | 0 より大きい `HeightAverage` の中央値 |
| `height_average_max` | 0 より大きい `HeightAverage` の最大値 |
| `rt_range` | RT の `(min, max)` |
| `mass_range` | m/z の `(min, max)` |
| `ion_modes` | イオンモード別の件数 |
| `annotated_count` | `Name` が空・null・`Unknown` でない件数（注釈の確度は考えない） |
| `annotation_rate` | `annotated_count / total_spots * 100` |
| `ontology_top` | 空でない `Ontology` の件数の上位 10 クラス（文字列 `Unknown` も数える） |
| `sn_median` | 0 より大きい `SignalToNoiseAve` の中央値 |
| `error` | データが空のときのエラーメッセージ |

#### 4.3 `format_spots_as_table()`

1 行目が列名、2 行目以降が 1 スポットの TSV/CSV 文字列。既定の列順は上の 22 列の順。浮動小数は最大 4 桁程度に丸め、区切り文字・
改行・引用符を含むセルは二重引用符で囲む。空値は空セル。丸め済みなので、精密な再計算には元の辞書値を使う。

#### 8.2 `.arf2` の要約（旧 `.arf2` 要約ツール）

旧 `.arf2` 要約ツールは撤去済み。件数・注釈の内訳は `load_dataset` の `summary`、列は `dataset_feature_search` が返す。

## 旧 `pai2` トピック（`.pai2` / 単一測定ファイルのピーク。§5, §8.3）

`.pai2`（単一測定ファイルのピーク）のパーサ出力。PAI2 を読む MCP ツールは無い（パーサは `metabolomix/pai2/reader.py`）。

> Key 番号表（MS-DIAL の `[Key(N)]` から抽出した一次資料。リーダーのインデックス定数を変える前に確かめる）:
> `docs/schema/ChromatogramPeakFeature.md`。

### 5. PAI2 (`metabolomix/pai2/reader.py`)

#### 5.1 `deserialize()` のピーク出力

型は `list[dict]`。1 要素は 1 測定ファイル内の 1 検出ピーク。

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
| `id` | int | その PAI2 ファイル内のピーク ID |
| `ion_mode` | `IonMode` Enum | `IonMode.Positive`（値 0）/ `IonMode.Negative`（値 1）/ `IonMode.Both`（値 2）。JSON にするには文字列か数値へ変換が要る |
| `name` | str | 候補化合物名。空文字は未注釈 |
| `formula` | str | 分子式 |
| `ontology` | str | 化学 / 脂質クラス |
| `smiles` | str | SMILES |
| `inchikey` | str | InChIKey |
| `adduct` | str | 付加体（`[M-H]-` など） |
| `collision_cross_section` | float | CCS。イオンモビリティを使わなければ 0 のことがある |
| `comment` | str | MS-DIAL 由来の注釈コメント |
| `has_msms` | bool | MS2 raw ID が非負、または collision-energy map が非空か。MS/MS 取得の情報の有無 |
| `ms2_raw_id` | int | 元データの MS2 raw spectrum ID。未設定は `-1` |
| `collision_energies` | list[float] | MS/MS 取得時の衝突エネルギーの一覧。負イオンモードでは `-42.0` など負値もある |
| `msms_peak_count` | int | PAI2 内 `msms_spectrum` の要素数。スペクトル本体は通常 DCL 側にあるので 0 |
| `msms_spectrum` | list[list] | PAI2 自体に格納された `[fragment_mz, intensity]`。通常は空で、本体は DCL 側 |

`time` 系の辞書には、入力にある非負の値だけが入る。

| サブキー | 意味 |
|---|---|
| `rt` | retention time |
| `ri` | retention index |
| `m/z` | 座標表現内の m/z |
| `dt` | drift time |

#### 5.2 フィルタ

`filter_features_by_params()` は元と同じピーク辞書のリストを返す。

| パラメータ | 意味 |
|---|---|
| `min_intensity` / `min_height` | `peak_height` の下限 |
| `min_sn` | `S/N` の下限。S/N の無いピークも除く |

#### 5.3 ピーク属性PCA（撤去済み）

`perform_pca_summary()` は撤去済み（コードに無い）。PAI2 は単一測定ファイルで行がピークなので、サンプル間比較（オミクス PCA）は
できない。サンプル間の多変量比較は ARF/ARF2（`dataset_pca` / `dataset_compare`）で行う。

#### 5.4 `inspect_peak_details()`

成功時は `{"status":"success", "matches":[...], "note":...}`。`peak_name` は大文字小文字を無視した部分一致で、複数件返り得る。
検索条件なしは `status=error`、該当なしは `status=not_found`。

| `matches[]` キー | 意味 |
|---|---|
| `id` | PAI2 のピーク ID |
| `name` | 候補化合物名 |
| `m/z` | ピーク m/z |
| `rt` | ピーク頂点 RT |
| `height` | ピーク高さ |
| `area` | ピーク面積 |
| `signal_to_noise` | S/N |
| `has_msms_like_fields` | キー名に `msms`・`fragment`・`spectrum` のどれかを含むフィールドがあるか。**実スペクトルが非空かは見ない** |
| `formula` | 分子式 |
| `adduct` | 付加体 |
| `comment` | コメント |

#### 5.5 `summarize_pai2_inventory()`

単一測定ファイルのピーク分布の在庫要約（撤去前の PAI2 ツールが返していたもの）。

| キー | 意味 |
|---|---|
| `total_peaks` | フィルタ後のピーク総数 |
| `annotated_peaks` | 確定注釈とみなせるピーク数（空文字・`Unknown`・`no MS2:`・`low score:` を除く） |
| `annotated_fraction` | `annotated_peaks / total_peaks` |
| `ion_mode_counts` | イオンモード別の件数 |
| `mz_range` / `rt_range` | m/z・RT の `{min, max}` |
| `height_summary` | ピーク高さの `{min, median, max}` |
| `sn_summary` | S/N の取得割合と `{min, median, max}` |
| `top_by_height` | 強度上位 10 ピーク（`id` / `name` / `m/z` / `rt` / `height` / `signal_to_noise`） |
| `note` | PAI2 は単一サンプルでオミクス PCA 不能・MS/MS は `.dcl` 参照、の断り書き |

#### 8.3 PAI2 の在庫要約ツール（撤去済み）

PAI2 のツールは撤去済み。在庫要約は `summarize_pai2_inventory()`（§5.5）。複数サンプルの比較はデータセットの `dataset_pca` /
`dataset_compare`。MS/MS は同名の `.dcl`（`dcl_index` がリスト順に対応）を参照する。
