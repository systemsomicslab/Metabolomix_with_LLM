# ARF (`.arf` / サンプル別ピーク)

`.arf` パーサ、`arf_parser`、前処理・QC（`arf_preprocess`）、差次的解析（`arf_differential`）の出力定義。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。
>
> 対応する Key 番号表: `docs/schema/AlignmentChromPeakFeature.md`（MS-DIAL の `[Key(N)]` から抽出した一次資料）。
> リーダーのインデックス定数を変更するときは必ずこちらを先に確認する。

## 3. ARF (`metabolomix/arf/reader.py`)

### 3.1 `deserialize()` のスポット出力

型は `list[dict]`。**1要素は1アラインメントスポット**であり、全サンプルのピークを `AlignedPeakProperties` に保持する。

| キー | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | パーサーがスポット順に付与する 0 始まりのマスターID |
| `AlignmentID` | int | 現実装では `MasterAlignmentID` と同じ連番 |
| `RT` | float/null | 代表サンプルから取得したスポット代表RT（min） |
| `MassCenter` | float/null | 代表サンプルから取得したスポット代表 m/z |
| `IonMode` | str | `Positive`、`Negative`、`Both`、`Unknown` のいずれか |
| `Name` | str/null | 代表アノテーション名。空文字の場合は未注釈 |
| `HeightAverage` | float/null | 現実装ではグループ先頭サンプルの `height` を格納する。名前に反して全サンプル平均を再計算していない |
| `AlignedPeakProperties` | list[list] | そのスポットに属する全サンプルの生 MessagePack 配列。LLM 解析では通常、次節の表形式を使う |
| `TagIds` / `Tags` | list | アラインメント結果の `*_tags.xml` から `MasterAlignmentID` で結合したタグ |
| `SamplePeakTags` | dict | サンプル別 `*_tags.xml` から `FileName` と `MasterPeakID` で結合した、タグ付きピークのみの辞書 |
| `SampleClasses` | dict | `.mddata` の `AnalysisFileBean` と `FileID`/`FileName` で結合したサンプルClass IDメタデータ |

### 3.2 `extract_peak_properties()` の表/CSV

型は `pandas.DataFrame`。**1行は「1アラインメントスポット × 1サンプル」の1ピーク**である。同じ `MasterAlignmentID` がサンプル数だけ繰り返される。

| 列 | 型 | 意味 |
|---|---|---|
| `MasterAlignmentID` | int | 行が属するアラインメントスポットID |
| `AlignmentID` | int | 現実装では `MasterAlignmentID` と同じ値 |
| `SpotRT` | float/null | アラインメントスポット全体の代表RT（min） |
| `SpotMassCenter` | float/null | アラインメントスポット全体の代表 m/z |
| `IonMode` | str | スポットのイオンモード |
| `CompoundName` | str/null | スポットの候補化合物名。空文字は未注釈 |
| `AlignmentTags` | list[str] | アラインメントスポットに付与されたMS-DIALタグ |
| `SampleIndex` | int | `AlignedPeakProperties` 内での 0 始まり位置。サンプル順序を示す |
| `FileName` | str | 測定ファイル名。取得できない場合は `Sample_<SampleIndex>` |
| `ClassID` | str/null | MS-DIALのFile property settingで指定した `AnalysisFileClass` |
| `PeakID` | int/null | 測定ファイル内のピークID。ギャップフィルでは負値になり得る |
| `FileID` | int/null | データセット内の測定ファイルID |
| `MasterPeakID` | int/null | 元ピークのマスターID。負値は未検出を補間したギャップフィルを示す |
| `PeakHeight` | float/null | 当該サンプルのピーク頂点強度 |
| `PeakArea` | float/null | 当該サンプルのピーク面積 |
| `PeakAreaAboveBaseline` | float/null | ベースラインより上のピーク面積 |
| `PeakMZ` | float/null | 当該サンプルで観測されたピーク m/z |
| `PeakRT` | float/null | 当該サンプルで観測されたピークRT（min） |
| `SignalToNoise` | float/null | `PeakShape[1]` 由来の S/N |
| `IsMsms` | bool | MS2 raw spectrum ID-to-collision-energy map が非空か。MS/MS取得情報があることを示すが、スペクトル本体は含まない |
| `IsGapFilled` | bool | `MasterPeakID < 0` か。`true` は実検出ではなくギャップフィルされた値 |
| `PeakTags` | list[str] | 当該サンプルピークの `MasterPeakID` に付与されたMS-DIALタグ |

### 3.3 `build_pca_matrix()` の行列

返り値は `(matrix, sample_names, feature_names)`。

| 出力 | 行・列の意味 |
|---|---|
| `matrix` | 2次元 `numpy.ndarray`。**行=サンプル、列=スポット×選択プロパティ** |
| `sample_names` | 行ラベル。通常は `FileName` |
| `feature_names` | 列ラベル。`Spot_<MasterAlignmentID>_<property>` 形式。例: `Spot_0_height` |

`property` は `_convert_to_alignment_feature()` の `height`, `area`, `area_above_baseline`, `m_z`, `rt`, `signal_to_noise` などを指定できる。欠けた値の扱いは 2 通りある。行はあるが指定プロパティの値が `None`（フィールド欠落・読めない値）のセルは **0.0** になる。ある試料の行がそのスポットに無いセルだけが欠損（NaN）になり、列平均で補完され、全試料で欠損の列は 0 となる。gap-fill（`IsGapFilled=true`）のセルは MS-DIAL の補間値をそのまま値として持つ。分散 0 の列は除外される。`min_detection_rate > 0` の場合は、非ギャップフィルサンプル率が閾値未満の列も除外される。

### 3.4 `run_pca()` と Loading 出力

PCA 前に各列を `StandardScaler` で標準化する。`log_transform=true` の場合は値を 1 以上にクリップして `log10` 変換してから標準化する。

| JSONキー | 形状 | 意味 |
|---|---|---|
| `components` | `[sample][PC]` | 各サンプルの主成分スコア。行順は `sample_names` と同じ |
| `explained_variance_ratio` | `[PC]` | 各主成分が説明する分散の比率。0..1 |
| `singular_values` | `[PC]` | 各主成分に対応する特異値 |
| `loadings` | `[PC][feature]` | 各主成分に対する各入力列の係数。列順は `feature_names` と同じ |

主成分の符号は数学的に反転可能なので、正負そのものよりサンプルと特徴量の相対関係を解釈する。

`get_pca_loading_features()` の各要素:

| キー | 意味 |
|---|---|
| `pc` | `PC1` などの主成分名 |
| `var_ratio` | 説明分散比を百分率にした値 |
| `positive` | Loading 値が大きい側の上位特徴量リスト |
| `negative` | Loading 値が小さい側の上位特徴量リスト |
| `positive/negative[].id` | `MasterAlignmentID` |
| `positive/negative[].value` | Loading 係数 |
| `positive/negative[].annotation` | スポットの `Name` |
| `positive/negative[].m_z` | スポット代表 m/z |
| `positive/negative[].rt` | スポット代表RT（min） |

### 3.5 ARF要約

`summarize_arf_data()` は次の辞書を返す。

| キー | 意味 |
|---|---|
| `total_peaks` | スポット数。名前は peaks だがサンプル別行数ではない |
| `rt_range` | スポット代表RTの `(min, max)` |
| `mass_range` | スポット代表 m/z の `(min, max)` |
| `height_average_mean` | `HeightAverage` の算術平均 |
| `height_average_max` | `HeightAverage` の最大値 |
| `ion_modes` | イオンモード別スポット数 |
| `named_compounds` | `Name` が null でない件数。**空文字も数えるため、真の注釈済み件数とは限らない** |

### 3.6 MS-DIALタグ

`metabolomix/msdial/tags.py` はARFと同じディレクトリの `*_tags.xml`（互換用に拡張子なしの `*_tags` も可）を読む。サンプル別ファイルは処理時刻の12桁接尾辞を除いた名前でARFの `FileName` と対応させ、XMLの `Peak/@Id` をARF行の `MasterPeakID` と結合する。アラインメント結果用ファイルは `Peak/@Id` を `MasterAlignmentID` と結合する。

タグ条件は `any`、`all`、`none`、`not_all` を使用できる。`sample_peak` スコープでは条件に一致しないサンプル別行をスポット内から除外し、`alignment_spot` スコープではスポット全体を除外する。タグ未付与ピークは `any`/`all` には一致せず、`none`/`not_all` には一致する。

サンプル名は完全一致を優先し、MS-DIAL処理時刻の12桁接尾辞を除く補助照合は一意に決まる場合だけ使用する。重複・曖昧照合はエラーとなる。タグファイル未対応サンプルは既定でエラーにし、`missing_sample_policy=exclude` で除外、`untagged` で明示的にタグなし扱いへ変更できる。サンプル数は固定せず、ARFから検出した件数を使用する。

### 3.7 Class ID

`metabolomix/msdial/classes.py` は `.mddata` の `MsdialDataStorageBase.Key0 AnalysisFiles` を読み、各 `AnalysisFileBean` の `AnalysisFileId`、`AnalysisFileName`、`AnalysisFileClass` を抽出する。`.mddata` は明示パス、`.mdproject` 内の参照、またはARFと同じディレクトリから解決する。

ARFサンプルとの結合は `FileID` を優先し、欠損時は正規化した `FileName` を使用する。両方が異なるサンプルへ解決された場合はエラーとする。`filter_arf_by_class_ids()` は選択Class ID以外のサンプル行を各スポットから除外し、その結果を `build_pca_matrix()` に渡すことでPCAの行をClass IDで選別できる。

### 8.1 `arf_parser()`

返り値はMarkdownテキスト1件（`str`）。総スポット数、適用フィルタ/手動除外の注記、Class ID分布、タグファイル対応数、タグ別件数、総サンプル別レコード数、平均サンプル数/スポット、PCA行列形状、PC1/PC2説明分散比、Loading上位、PCAスコアプロット用ブロック（群別サンプル数・描画指示・サンプル別座標の json 点列）を含む。`class_ids` を指定すると、選択したClass IDに属するサンプル行だけを残してPCAを実行する。複数Class IDはOR条件で、照合は大文字小文字を区別しない。`min_intensity`（スポット平均強度の下限）と `annotation_keyword`（脂質クラス/化合物名の部分一致）で生スポットを絞ってPCAをやり直せる（**フィルタ条件を変えたPCAのやり直しは本ツールの再呼び出しで行う**。ファイルはセッションキャッシュされ再パースは走らない）。手動除外（`arf_exclude`）も行列構築前に反映される。正規化・QC・欠損補完を経た「前処理後」行列でのPCAは `arf_preprocess` → `arf_pca_preprocessed`。`arf_list_classes()` は `.mddata` のパスとClass ID別サンプル数をJSONで返す。`arf_list_tags()` は現在のARFセッションについてタグ定義、サンプルファイル対応数、タグ付与数をJSONで返す。

脂質選択では、`annotation_keyword` の名前は同一バッチのARF2を優先し、ARF2がない／未注釈ならARF名へ戻る。`spot_ids=[10, 20]` はMasterAlignmentIDの完全一致、`ontologies=["PC"]` は脂質クラスの完全一致（大小文字を区別しない）で、各引数はAND条件。空リスト・負のIDなど不正な指定はエラーにする。`annotation_keyword="PC"` は従来どおりLPCにも部分一致するため、PCクラスだけを選ぶには `ontologies` を使う。

採用スポット数、注釈の出所、同一バッチARF2のパス、ARF/ARF2の名前不一致数とID先頭10件を表示する。選択用コピーには `arf_name`、`arf2_name`、`annotation_source`、`annotation_conflict` を保持する。元ARFの注釈・ピーク値は変更しない。

これは生スポットの選択であり、選択を指定／変更した場合は古い前処理行列・差次解析結果を無効化する。後続の解析には `arf_preprocess` が必要。部分集合をTIC正規化すれば分母が変わり、部分集合でBH補正すれば検定対象数が変わる。全体の正規化・FDRを維持した候補確認には、全体で解析・エクスポートした結果を参照し、部分集合で再計算した値と混同しない。

脂質選択を指定した場合、分散・検出率フィルタ後にPCAの必要次元数が足りなくても、サンプルが存在すれば生スポット選択は保存する。PCAをスキップした理由と行列形状を返し、PCA図・結果は生成しない。1脂質や定数の脂質もこの経路で後続の前処理へ渡せる。

PCAスコアプロット用ブロック（返り値の末尾。**サンプル別の座標点列を同梱する**）。
散布図は**この点列からクライアント側（チャット）で描く**。`save_figure(kind="pca")` は PNG
ファイルを reports/figures/ へ書き込む操作であり、ユーザーがファイルとしての図を明示的に
求めたときだけ使う——対話の中で図を見せる目的では呼ばない。同じ点列は
`session.arf.last_pca_plot` にも保持され、`save_figure(kind="pca")` はそこから描く。

| キー/行 | 意味 |
|---|---|
| タイトル行 | 図タイトル |
| `PC1 (x%) × PC2 (y%)` | PC1/PC2 説明分散率 |
| 群別サンプル数 | 群ラベルがあれば `群=件数` を列挙、無ければ総サンプル数 |
| 描画指示 | 上記座標から散布図を描く（group があれば群ごとに色分け・凡例付き）。PNG が必要なときのみ `save_figure(kind="pca")` |
| json 点列 | `x_label` / `y_label`（説明分散率つき軸名）と `points`。各点は `pc1` / `pc2` / `sample`、群分け時は `group` |

## 10. 前処理・QC（P2a）

`metabolomix/analysis/preprocessing.py`（MCP非依存の純ロジック層）と `metabolomix/arf/tools.py` の `arf_list_sample_roles()` / `arf_preprocess()` / `arf_pca_preprocessed()` が、ARFロード後のサンプル×特徴量行列に対する前処理・QCを担う。既定では**何も適用されない（opt-in）**。生行列を消費する `arf_parser` の既定挙動は変えない。

### 10.1 役割検出（sample/qc/blank）

`preprocessing.detect_sample_roles()` が、ファイル名と Class ID を `_` 区切りでトークン化し、大小無視で `qc`/`blank` トークンと照合してサンプルを `sample`/`qc`/`blank` に分類する（`blank` を `qc` より優先評価）。`arf_list_sample_roles()` はこの分類結果と役割別件数を、前処理適用前の確認用にJSONで返す。

### 10.2 前処理レシピ（`arf_preprocess`）

`arf_preprocess(normalize, blank_min_fold, drift_correct, max_qc_rsd, impute, props)` が、`preprocessing.preprocess()` に処理を委譲し、以下の順で適用する。

1. **ブランク除去**（`blank_min_fold` 指定時）: 生体試料平均 が `blank_min_fold` × ブランク平均 未満の特徴量を背景として除去。ブランク/生体試料のどちらかが無ければ未実施（caveat）。
2. **正規化**（`normalize="tic"|"median"|"pqn"|"none"`）: 行（サンプル）ごとのスケーリング。`tic`=行総和、`median`=行中央値、`pqn`=Probabilistic Quotient Normalization（参照はQC中央値、QCが無ければ全サンプル中央値）。**正規化係数が0または非有限のサンプル（未検出=0が過半で行中央値=0 になる疎な試料など）は、行全体をNaN化して破棄せず未正規化のまま残置し、`report["unscaled_samples"]` と caveat で明示する**（`median`/`pqn` で起こりやすい。`tic`=行総和は総和>0のため安全）。旧実装は該当行をNaNで全消去し、疎データで多数の試料を無言で失っていた。
3. **QC ドリフト補正**（QC-RLSC 相当の移動中央値版。`drift_correct=True` 指定時）: QCを注入順（`analytical_order`、`.mddata` 由来）に並べ、QC 値の移動中央値（窓 5、QC 4 本以上）で系統ドリフトを推定し、`np.interp` で全注入へ区分線形に補間して `median(QC) / trend` を掛ける。特徴量ごとに全サンプルを 1 系列として補正し、バッチ別には補正しない。LOESS（局所重み付き回帰）は使わない。**注入順が全サンプルで取得できない、またはQCが最小数未満なら未実施**（caveat）。加えて **QC が試料列に挿入されていない設計でも未実施**（`status="skipped"`）: QC-RLSC は QC が試料の前後に散在することを前提とし、QC を全試料の後にまとめて流した設計（実例: 試料 1–48 → Blank 49 → QC 50–56）では `np.interp` が端値で頭打ちになり、補正した外見だけが残る。判定は `report["qc_interspersion"]`（`covered`＝QC注入順区間に入る試料数 / `qc_range` / `sample_range`）。`covered` が試料の半数未満なら適用はするが「外挿補正」caveat を付す。
4. **QC RSDフィルタ**（`max_qc_rsd` 指定時）: QC群での相対標準偏差（SD/mean）が閾値を超える特徴量を除去。QCが無い/不足なら未実施（caveat）。
5. **欠損補完**（`impute="half_min"|"knn"|"column_mean"|"none"`、既定 `half_min`）: 行列生成後に残るNaNを補完。`half_min`=特徴量最小値の半分（既定）、`knn`=sklearn `KNNImputer`、`column_mean`=列平均（旧実装互換）、`none`=補完しない。
6. **ブランク行の除外**: 上記1でブランクを背景除去の**参照**として使い終えたあと、`preprocessing.drop_samples_by_role()` がブランクの**行そのもの**を解析行列から外す。除外内訳は `report["excluded_from_matrix"]`（`{role: [sample_name, ...]}`）と caveat に出る。

> **役割別の行の扱い（重要）**: 前処理後行列 `session.arf.feature_matrix` の行は **生体試料 + QC** であり、**ブランクは含まれない**。ブランクは生体試料と桁違いに総強度が低く、残すと PCA の PC1 を支配して群分離の解釈が壊れるため外す。一方 **QC は残す**——QC クラスタの締まり具合を PCA で見るのは品質確認の定番手段だからである。したがって `arf_pca_preprocessed()` のスコアプロットには QC 点が含まれる（群ラベルは QC の Class ID）。群平均に QC が混ざると困る `arf_differential()` 側は、別途 QC を比較群から外す（§11.1）。

処理結果は `session.arf.feature_matrix`（前処理後行列）・`session.arf.pp_sample_names`・`session.arf.pp_feature_names`・`session.arf.sample_meta`・`session.arf.preprocessing_recipe` に保存され、以降の `arf_pca_preprocessed()` や将来の差次的解析（P2b）はこの前処理後行列を消費する。適用したレシピそのものが `session.arf.preprocessing_recipe` に記録され、`arf_pca_preprocessed()` の出力にも「前処理レシピ」として明示される。

### 10.3 caveatの扱い

QC/ブランク/注入順のいずれかが欠けているためにスキップされたステップは、無言で無視されるのではなく `report["caveats"]`（`arf_preprocess()` のJSON応答）に文言として残る（例:「注入順が欠落、または QC が不足のためドリフト補正は未実施。」）。LLMはこれらのcaveatを解釈結果や報告書の注意点として引用すべきである。追加で前景化される caveat:

- **失敗 QC 注入**: 前処理の最初（正規化でスケールが動く前の生強度）に `preprocessing.detect_failed_qc()` が、総強度が QC 中央値の 20% 未満の QC を名指しする（`steps["qc_health"]`）。失敗注入を残したまま `max_qc_rsd` を掛けると QC の RSD が全特徴で跳ね上がり、ほぼ全特徴が除去される（実測: kidney aging NEG で 1345 → 51）。**「閾値が厳しすぎる」ように見える現象の真因は QC 側にあることが多い**ので、`arf_exclude` で除外してから前処理し直す。
- **プールQC の層別**: QC が複数バッチ（日付）に分かれる場合に加え、**QC 試料名の層別**（部位別 QC 等。`preprocessing.detect_qc_strata` が `20240311_QC_Cerebellum_ICR_NEG_1` → `cerebellum_icr` のように日付/`qc`/極性/数字を除いた残りで判定）も検出し、「全 QC を1系列扱いするドリフト補正/RSD は近似」と警告する。
- **正規化での試料脱落**: `normalize` の `unscaled_samples`（係数0/非有限で未正規化残置した試料数）に対応する caveat（§10.2）。
- **過度な特徴量除去**: フィルタ後に残存0件なら「全特徴が除去（閾値が厳しすぎる可能性、解析不能）」、特徴量の90%超が除去なら残存割合を注記する。除去総数は `report["features_removed_total"]`（= before − after）で参照する。**各 `steps[*]["removed"]` はフィルタごとの独立マスク件数で重複し得るため加算しないこと**（blank と qc_rsd の removed 合計が総数を超えることがある）。

なお `load_dataset` の複数バッチ告知は、解析対象である **`.arf`/`.arf2` のバッチ**にのみ基づく（`.mddata`/`.mdproject`/`.msp2`/`.pai2` 等も `AlignmentResult` 形式のタイムスタンプを持つため、拡張子で限定しないと告知バッチが実際に解析する `.arf` とズレる）。

### 10.4 `arf_pca_preprocessed()`

前処理後行列が無い（`session.arf.feature_matrix is None`）場合はエラーメッセージ1件を返す。あれば `metabolomix/analysis/pca.py` の `run_pca` でPCAを実行し、`arf_parser` と同じ整形ヘルパー（スコアプロット用JSON、Loadings上位）を使って結果を返す。出力テキストの構造・キー意味は8.1節のスコアプロット用JSONと同一（座標点列を同梱し、散布図はそこからチャットで描く。`save_figure(kind="pca")` はPNGファイルが明示的に求められたときだけ）。生スポットの選択とは別の計算経路であり、選択を変更した場合は `arf_preprocess()` を再実行する必要がある。

### 10.5 手動サンプル/ピーク除外（`arf_exclude`）

PCAスコアプロットで明らかに外れた1サンプルや、特定のピーク（スポット）を**名前/IDで手動除外**するためのツール。`metabolomix/arf/exclusions.py`（MCP非依存の純ロジック層、`prune_spots()` / `roster()`）と `metabolomix/arf/tools.py` の `arf_exclude()` が担う。除外は**可逆・非破壊**で、`session.arf.filtered_features` 自体は変更しない。

`arf_exclude(exclude_samples=None, exclude_spots=None, mode="add")` は JSON を返す。

- **`exclude_samples`**: 除外するサンプル名（`file_name`、完全一致）のリスト。
- **`exclude_spots`**: 除外するスポットの `MasterAlignmentID`（int）のリスト。
- **`mode`**: `add`（既定・追加）/ `remove`（再包含）/ `clear`（全消去）/ `list`（現状表示のみ）。
- 現データに存在しない指定は `unmatched_samples` / `unmatched_spots` として警告に載せ、一致分のみ集合へ反映する（タイプミスに寛容）。
- 応答キー: `status` / `mode` / `excluded_samples` / `excluded_spots` / `samples_before` / `samples_after` / `spots_before` / `spots_after` / `unmatched_samples` / `unmatched_spots` / `caveats`。残サンプルまたは残スポットが0件になる指定には caveat が付く。

除外集合は `session.arf.excluded_samples`（`file_name` 集合）と `session.arf.excluded_spots`（`MasterAlignmentID` 集合）に保持され、新ファイルロード（`load_data` のキャッシュミス経路）でリセットされる。行列を組む直前に `exclusions.prune_spots()` が適用され、**`arf_parser`**・**`arf_preprocess`**（→ `arf_pca_preprocessed` / `arf_differential`）がいずれも自動的に除外を反映する。除外が有効なとき、それぞれの出力に「ユーザ手動除外: サンプル N 件 / スポット M 件」の注記が付く。`arf_list_sample_roles()` は各サンプルに `excluded: true/false` を付して現在の除外状態を示す。

運用フロー: `arf_parser`（全体PCAで外れ俯瞰）→ `arf_exclude(exclude_samples=[...])` → `arf_parser` 再呼び出し または `arf_preprocess`＋`arf_pca_preprocessed` で除外後PCAを確認 → `arf_differential`。戻したいときは `mode="remove"` / `mode="clear"`。

## 11. 差次的解析（P2b）

`metabolomix/analysis/differential.py`（MCP非依存の純ロジック層）と、`metabolomix/arf/tools.py` の `arf_differential()` / `metabolomix/tools/reports.py` の `save_figure()`（`kind="volcano"`）が、前処理後のサンプル×特徴量行列に対する群間比較を担う。既定挙動・既存ツールは不変で、明示呼び出し時のみ作用する。

### 11.1 群ラベルの由来

群ラベルは `session.arf.sample_meta[<sample>]["group"]`（ファイル名由来の factor トークン / Class ID 機構、`metabolomix/msdial/classes.py` の `assign_sample_groups`）から取得する。バッチは同 `sample_meta` の `batch`（ファイル名中の8桁日付）。`arf_differential()` は `session.arf.feature_matrix`（前処理後行列）を消費し、無ければエラーを返す（先に `arf_preprocess()` が必要）。

`sample_meta[...]["group"]` は常に**完全な Class ID**（例 `24M_GF_F`）である。一方 `group_a` / `group_b` は**因子トークンによるプール指定**を受け付ける（`sample_factors.expand_sample_specs`。トークンは Class ID とサンプル名の両方から解決される）:

- `group_a="24M", group_b="9w"` → `24M_*` を全てプールし `9w_*` と比較（多因子デザインで主効果を見る正しい経路）
- `group_a="24M_GF"` のように複数トークンを `_` で繋ぐと AND 絞り込み（`24M` かつ `GF`）
- 完全な Class ID を渡せば従来どおりその1水準のみ
- `group_a="ILG_6h", group_b="ILG_0h"` のように、**Class ID が同一でサンプル名の因子（時点等）だけが違う2群**も切り出せる
- 応答の `resolved_samples` が**群の定義そのもの**（実際に比較したサンプル名）、`n_a` / `n_b` が実 n。プールに2件以上のサンプルが入れば「プール群として解決」caveat を付す
- 応答の `resolved_class_ids` は**選択されたサンプルが持つ Class ID**であって群の定義ではない。両群で同じ Class ID になり得る（上の `ILG_6h` vs `ILG_0h`）ため、縮退比較と誤読しないこと。その場合は「Class ID では区別されず」caveat が付く

**一致ゼロ・両群が同じ Class ID を掴む指定は `status="error"` で落とす**（成功扱いで n=0 を返すと「有意0件＝群間差なし」と誤読されるため）。`24M` と `GF` は `24M_GF_*` を共有するので排他ではなくエラーになる。なお交互作用検定は依然として提供しない。

**QC・ブランクは比較群に入らない**。`sample_meta[...]["role"]` が `qc` / `blank` の試料は群ラベルを `None` にしてから展開するため、どちらのプールにも寄らない。QC の Class ID が指定トークンを含む場合（例 `QC_24M` と `group_a="24M"`）に黙って群平均へ混ざるのを防ぐための明示的な封鎖であり、除外内訳は caveat「比較対象から除外（生体試料でないため）」に出る。`n_a` / `n_b` は除外後の実 n。

### 11.1.1 上位ヒットの命名（ARF/ARF2 橋渡し）

`summary.top` の各行には `spot_id` / `name` / `name_source` が付く。`name_source="arf"` は ARF スポットの `Name`、`"arf2"` は ARF が `Unknown` のときに **同一アラインメントの兄弟 `.arf2`** から補った注釈（`ontology` も併記）。

**ARF と ARF2 は同じ `MasterAlignmentID` を指しながら代表 `Name` が食い違うことがある。** 実例（kidney aging, NEG）: Spot 474 は ARF 側 `Unknown`、ARF2 側 `SL 33:0;O|SL 17:0;O/16:0`（Ontology=`SL`）。ARF だけを見ると最大効果量の特徴が無名のまま残り、生物学的解釈に到達できない。

補完元はGUIの `AlignmentResult_<timestamp>` またはConsoleの `AlignResult-<digits>` の語幹が完全一致する隣接 `.arf2` に限定する（`MasterAlignmentID` はアラインメント実行ごとに振り直されるため、別バッチの `.arf2` を引くと ID 対応が黙って崩れる）。兄弟が無ければ補完せず `name=null` のままにする。

### 11.2 統計

- **2群比較**（`group_a` と `group_b` を指定）: 特徴量ごとに Welch t 検定（等分散を仮定しない）と log2 fold change を計算する。**正=群Bで高い（上昇）**。group_a が基準（対照）、group_b が比較対象。2026-08-31 に慣習（log2FC=log2(比較対象/基準)）へ合わせて符号を反転した——それ以前の出力とは符号が逆。擬似カウントは既定 1.0。定義は `log_transform` で 2 通りある:
  - `log_transform=true`（**既定**）: 各値を `log2(max(x, 0) + 擬似カウント)` に変換して Welch 検定し、`log2fc = mean(log2 値_B) − mean(log2 値_A)`（＝ x+1 の幾何平均比の log2）。
  - `log_transform=false`: 生値で Welch 検定し、`log2fc = log2((mean_b + 擬似カウント) / (mean_a + 擬似カウント))`（算術平均比）。
  - `mean_a` / `mean_b` はどちらの場合も生値の算術平均。
  - 結果行の `t` は `(平均_B − 平均_A) / SE` で、**log2fc と同じ向き**（正 = 群Bで高い。2026-09-27 に向きを揃えた——それ以前の `t` は符号が逆）。
  - 小n・分散0・全欠損は `p=NaN`。
- **多群ANOVAは現状非対応**: MS-DIAL メタに「因子（加齢/菌叢等）→水準」の対応が無く、因子を安全に選べない（誤って全 Class ID を水準にした結果を返さないよう封鎖）。3群以上を比べたいときは `group_a`/`group_b` の因子トークン・プール指定で関心のある2群を切り出す。`differential.one_way_anova()` 自体は関数として残るが、MCP からは露出しない。
- **多重検定補正**: いずれも Benjamini-Hochberg で `p → q`（FDR）を付与（NaN は補正から除外し位置は保持）。p値は scipy があればそれで、無ければ自前の t 分布裾確率（正則化不完全ベータ関数）で計算し、どちらも正確（近似ではない）。
- **volcano**: 2群比較のみ。各点は `feature` / `log2fc` / `neg_log10_p` / `sig`（`up`=q≤閾値かつlog2fc≥+閾値（群Bで高い＝上昇） / `down`=q≤閾値かつlog2fc≤−閾値（群Aで高い＝低下） / `ns`）。全特徴分の点列は `session.arf.last_differential["volcano"]` に保持し、`arf_plot_volcano()`（構造化点列）と `save_figure(kind="volcano")`（PNG）の両方がここから読む。**`arf_differential()` の応答 payload には全量 volcano を同梱せず**、`summary`（`n_tested`/`n_significant`/`n_up`/`n_down`＋有意上位 `top`）中心の要約と `volcano_note` のみを返す（先頭の結論が巨大配列＋文脈切り詰めで埋没し「全て ns」と誤読される退行を避けるため）。

### 11.2.1 `arf_plot_volcano` の返り値（`lipidmix.volcano.v1`）

`arf_differential`（2群）の後に呼ぶ read-only ツール。**既定（`output="image"`）では
サーバ側で描いた PNG と件数入りの1行キャプションを返す**ので、以下の表は
`output="payload"`（または env `LIPIDMIX_PLOT_OUTPUT=payload`）でクライアントが自分で
描くときの契約。ファイルとして PNG を残したいときだけ `save_figure(kind="volcano")` を呼ぶ。

画像モードでは間引きをせず全特徴を描く。**有意件数はキャプションの
`up=` / `down=` / `ns=` を読むこと**（図の点を数えない）。

| フィールド | 意味 |
|------------|------|
| `points[].feature` | 特徴量名（`differential` の `feature`） |
| `points[].log2fc` | log2 fold change（x軸）。**正=群Bで高い（上昇）** |
| `points[].neg_log10_p` | `-log10(p)`（y軸）。`q` ではなく **`p`** |
| `points[].sig` | `up` / `down` / `ns` |
| `thresholds` | 判定に使った `q` と `log2fc` のしきい値 |
| `comparison` | `group_a` / `group_b` / `n_a` / `n_b`（プール解決後の実n） |
| `render_hints.guides.x` | 縦破線の位置（`±log2fc_threshold`） |
| `render_hints.guides.y` | 横破線の位置（`-log10(q_threshold)`）。y軸は `-log10(p)` なので**目安**であり有意判定と厳密には一致しない |
| `selection.total` | 元の点数（＝検定対象の特徴量数） |
| `selection.plotted` | 実際に返した点数 |
| `selection.significant_total` / `significant_plotted` | 有意点の総数と返した数。**常に一致する**（有意点は間引かない） |
| `selection.ns_total` / `ns_plotted` | `ns` 点の総数と返した数。間引きが起きるとここが乖離する |
| `selection.max_points` | 要求した上限点数（既定800。`up`/`down` は上限に関わらず全件残る） |
| `selection.dropped_nonfinite` | `log2fc` か `p` が有限でなく描画対象外にした件数。**「有意でない」という意味ではない**（分散0・欠損・片群のみ検出などで検定不能だったもの） |

`selection.plotted < selection.total` のとき、図は全点ではない。有意件数の判断は必ず
`selection.significant_total` を見ること（画面や `points` の点を数えてはいけない）。
間引きは `ns` 点に対する等間隔サンプリングで決定的（同じ入力なら同じ点集合）。

### 11.3 必須caveat

`arf_differential()` の応答 `caveats` には、該当時に以下を前景化する（無言で握りつぶさない）:

1. **交絡（群⟂バッチ）**: 各群が単一バッチに偏る場合、「処理効果と測定バッチを分離できない」旨を警告（`check_confounding`）。例: `2_lipidome_lcms/NEG` は control/LPS=20220901・ILG/G_uralensis=20220902 で交絡。**判定は必ずプール解決後の（＝実際に比較した）2群に対して行う**。プール前の Class ID 単位で見ると細粒度ラベルほど各群が単一バッチになりやすく、プールすれば両群ともバッチ混在という健全な設計を交絡と誤報するため。判定に使うのは群ラベルが `None` でない試料のみ（QC/ブランク除外後）。
2. **正規化状態**: `session.arf.preprocessing_recipe` に正規化が含まれなければ「未正規化データの log2FC は測定量差を含み得る」と警告。
3. **群サイズ不足**: いずれかの群が n<2 なら「各群 n>=2 が必要（群名の誤り／前処理での試料脱落の可能性）」と警告。n>=2 かつ n<4 なら小n（検出力の限界）を注記。
4. **退化（検定不能）**: 検定できた特徴が0件なら「全特徴で p=NaN。群が空・分散0・正規化での試料NaN化の可能性。『有意0件』を『群間差なし』と解釈しない」と警告。0件でなくても特徴数の20%未満しか検定できなければ注記する。これにより「本当に有意差が無い（n_tested 健全）」と「そもそも検定できていない（n_tested≈0）」を区別できる。

LLMはこれらを解釈結果・報告書の注意点として必ず引用すること。統計値は「事実」だが、交絡・小nの下での因果的解釈は保留し、人間の判断に委ねる（既存の分業に整合）。

### 11.4 `arf_export_differential` — 差次的結果のエクスポート契約

`arf_preprocess` → `arf_differential`（2群）の後に呼ぶ。同一アラインメントの
兄弟 `.arf2` から InChIKey・Ontology・m/z・RT を `MasterAlignmentID` で結合し、
1 ファイルに書き出す。**兄弟 `.arf2` が無ければ書き出さない**（InChIKey 空欄の行を
出すと、下流で「パスウェイが無い化合物」と区別が付かなくなるため）。
直近の差次的結果が現行の `contract_version` / `log2fc_sign` と一致しない場合も、
向きを偽装せず `arf_differential` の再実行を要求する。InChIKey 付き行が 0 件なら、
下流が拒否する本文 0 行のファイルを成功扱いで残さない。

`#` 始まりのメタ行に来歴（`contract_version` / `group_a` / `group_b` /
`log2fc_sign` / 閾値 / `n_features_total` / `n_with_inchikey` / `n_unannotated`）を置き、
続けて TSV 本体を置く。列は
`spot_id / name / name_source / ontology / inchikey / inchikey_source / msi_level /
mz / rt / log2fc / p_value / q_value / mean_a / mean_b / significant`。

**有意な行だけでなく、InChIKey が付いた全行を書き出す。** 下流の濃縮解析は
「検出された化合物」を背景に取る必要があり、有意な行だけでは背景が作れない。

`msi_level` は `.arf2` 由来の注釈確度であり、**MS/MS の有無ではない**
（`.arf2` は MS/MS 取得フラグを持たない）。`arf2_annotate_identities` と同じ
保守的なクラス上限を使う。

### 11.5 `arf_plot_group_intensity` の返り値（`lipidmix.group_intensity.v1`）

選んだクラス・分子種が、指定した試料群のそれぞれで**見つかるか・どれくらいあるか**を示す
read-only ツール。ARF 経路のみ（読み込み済み `.arf` と、同じアラインメントの兄弟 `.arf2` が要る）。
**既定（`output="image"`）はサーバ側で描いた PNG と 1 行のキャプション**。キャプションには
項目ごとのスポット数（MS/MS の件数が分かるときは括弧内にその数。0 も出す）、群ごとの試料数、除外の 5 種類の
件数（手動除外・内部標準・判断・自動判定・標準液）、N.D. の項目名、`caveats` が入る。以下の表は `output="payload"`（または env
`LIPIDMIX_PLOT_OUTPUT=payload`）の契約。ファイルとして PNG が要るときは
`save_figure(kind="group_intensity")`（`reports/figures/<analysis_id>_group_intensity.png` を dpi 300 で、
同名の `.svg` と一緒に書く。直前の `arf_plot_group_intensity` が無ければ `missing_state`）。

**エラーと状態**: 入力の誤り（群が 1 試料にも当たらない・`output` の値が不正・`detection_limit` が 0 以下や非有限・`ncols` が 1 未満など）は例外ではなく `{"status": "error", "message": ...}` で返す。呼び出しが（missing_state を含め）失敗したときは、前回の図を `save_figure(kind="group_intensity")` が保存しないようセッションの図を破棄する。画像モードでは描画に成功してから図を保持する。群の指定が `blank` / `qc` の役割で全滅したときは、`blank` / `qc` を書くよう案内を添える。

**検定はしない**。1 パネル = 1 項目、1 点 = 1 試料、群ごとに log10 空間の平均 ± SD を添えるだけで、
群間の有意差は述べない（検定は `arf_differential`。小 n の p 値が図に独り歩きするのを避ける）。

| フィールド | 意味 |
|------------|------|
| `plot_schema` | `lipidmix.group_intensity.v1` |
| `value` | 点の値の定義: `sum of PeakHeight (log10 on the plot)`。1 点 = その試料での、項目に当たったスポットの **PeakHeight の合計**（`.arf` の試料別行。gap-fill の値も含める） |
| `detection_limit` | 検出下限（図の破線）。無ければ null |
| `detection_limit_source` | `argument`（`detection_limit` 引数）/ `param_file`（`.arf2` と同じフォルダの param ファイルの `Minimum peak height`）/ null（下限なし）。引数が優先 |
| `groups[]` | `label`（指定した群の文字列）と `samples`（解決した試料名）。指定順。`arf_exclude` した試料は含まない |
| `low_reliability_samples` | 白抜きで描き、平均 ± SD から外した試料（抽出量が少なかった等）。解決後の試料名 |
| `items[].item` | 指定した項目の文字列（指定順。1〜30 件） |
| `items[].parts[]` | `+` で分けた各部分。`part`・`kind`（`class` = `.arf2` の Ontology に完全一致、`name` = 分子種名に完全一致、`none` = 当たらず。名前が当たっても全スポットが除外で消えたときも `none`。そのスポットは `excluded` に載る）・`n_spots` |
| `items[].spots[]` | 当たったスポット `spot_id` / `name` / `ontology`（キュレーションの `assign` は付け替え後の値） |
| `items[].n_spots` | 当たったスポット数（アダクト違い・重複スポットは合計に入る） |
| `items[].n_spots_msms` | 当たったスポットのうち MS/MS の裏付けがある数。**スポットが 1 つ以上あって 0 のパネルだけ図で灰色になり「MS1-only (unconfirmed)」と書かれる**（MS1 の一致だけで、MS/MS では確かめていない同定。スポットが 0 の N.D. パネルは灰色にしない） |
| `items[].detected` | どこかの試料で値が 0 より大きいか |
| `items[].groups[].samples[]` | `sample` / `value`（PeakHeight 合計） / `gap_filled_fraction`（合計のうち gap-fill の値の割合。値 0 なら null） / `low_reliability` |
| `items[].groups[].log10_mean` / `log10_sd` / `n_in_stats` | log10 空間の平均・SD と、その計算に使った試料数（値 > 0 かつ低信頼でない試料）。n = 1 なら SD は null、0 なら平均も null |
| `excluded` | 除いたスポットの `#<spot_id> <name>`。理由ごとに 5 種類（`manual` / `internal_standard` / `curation` / `auto_likely_wrong` / `standard_only`。下記） |
| `caveats` | 試料が複数の群に当たった（群の組ごとに 1 件。重なった試料数と最大 3 試料名）、同じ群を重複指定した、孤立した判断記録がある、レビューが無く自動判定を除けなかった、等 |

**N.D. には 2 種類ある**。(a) 項目に当たるスポットが無い（`items[].n_spots` = 0。名前が当たらなければ `parts[].kind` は `none`。
クラスに当たっても除外で全スポットが消えた場合は `class` のまま 0 になる）。(b) スポットは当たったが、描いた試料のどれでも値が 0（`n_spots` > 0 で
`detected` が false）。どちらも N.D. パネルとして残し、エラーにしない（「どの群にも無い」ことも結果）。
(b) は「その分子種が存在しない」ではなく、**この群の試料では検出されなかった**という意味。

**点の読み方**: 縦軸は log10(PeakHeight) で**全パネル共通**（パネルごとに拡大すると、強度 10 前後の
ノイズが信号に見えて「見つかるか」の判断を誤らせる）。値 0 の点は軸の底に ▽ で置く。
`gap_filled_fraction` が **0.5 を超える点は ◆**（合計の過半が gap-fill、つまり検出ピークではなく
後から埋めた値）。◆ は検出の証拠として弱いので、群の「見つかった」を読むときに割り引く。

**`excluded` の 5 種類**（キーと条件）:

| キー | 除くもの |
|------|----------|
| `manual` | `arf_exclude` でスポットとして除いたもの（セッションの `excluded_spots`）。クラスにも名前指定の部分にも効く |
| `internal_standard` | 名前に `(d<数字>)`（`PC 33:1(d7)` 等）を持つ標識内部標準。**クラスとして当たった分だけ**除く。名前で直接指定した部分は描く |
| `curation` | `apply_curation=True`（既定）のとき、人の判断が `wrong` / `redundant` のスポット（`clear` 済みは無効）。`assign` は記録した名前・クラスで扱う。`suspect` は描く。名前指定の部分にも効く |
| `auto_likely_wrong` | `exclude_auto_likely_wrong=True` のとき、このアラインメント（`.arf2` の sha256 一致）の最新 `curation_review` で `likely_wrong` のスポット。名前指定の部分にも効く。レビューが無ければ除かず `caveats` に出す |
| `standard_only` | `standard_samples` を指定したとき、その試料での最大高さが、描く群の試料での最大高さの 10 倍以上のスポット（クラスとして当たった分だけ）。システムチェック標準液の奇数鎖標準物質が試料のクラス合計に混ざるのを防ぐ。判定には `arf_exclude` 済みの試料の行も使う |

除外は**解析の前提を変える**ので、図を報告に載せるときは `excluded` の件数と理由を併記すること。

### 11.6 `arf_plot_species` の返り値（`lipidmix.species_intensity.v1`）

選んだクラスを**分子種（アラインメントスポット）に展開**し、1 分子種 1 パネルで、群ごとの試料別の値
（割合 % または PeakHeight）を並べる read-only ツール。ARF 経路のみ（読み込み済み `.arf` と兄弟 `.arf2` が要る）。
`arf_plot_group_intensity`（クラス合計）と違い、**スポットごとに**描く。同じ名前でも付加イオンが違えば別パネル。
**既定（`output="image"`）はサーバ側で描いた PNG と 1 行のキャプション**（パネル数・群の試料数・除外の 6 種類の件数・
分母が 0 で描かなかった試料・`caveats`）。以下は `output="payload"`（または env `LIPIDMIX_PLOT_OUTPUT=payload`）の契約。
ファイルは `save_figure(kind="species")`（`reports/figures/<analysis_id>_species.png` を dpi 300 で、同名の `.svg` と
一緒に書く）。入力の誤り・パネルが 40 枚を超える・判断の記録ファイルが読めないときは
`{"status": "error", "message": ...}`。失敗した呼び出しの後は前回の図を破棄する。**検定はしない**。

| フィールド | 意味 |
|------------|------|
| `plot_schema` | `lipidmix.species_intensity.v1` |
| `value` | `share`（試料ごとに 高さ ÷ 分母の合計 × 100）または `height`（PeakHeight） |
| `stat` | 平均と SD の定義: `mean ± SD of share (%)` / `mean ± SD of log10(peak height)`。`groups[].mean` / `sd` はこれに従う |
| `share_basis` | `value="share"` のときだけ。`items`（指定した分母の項目。省略時は null で、描く項目が分母）と `spot_ids`（分母にしたスポット）。`height` では null |
| `groups[]` | `label` と `samples`（解決した試料名。指定順。`arf_exclude` した試料は含まない） |
| `low_reliability_samples` | 白抜きで描き、平均・SD から外した試料 |
| `zero_denominator_samples` | 分母（分母スポットの高さの合計）が 0 で割合を計算できなかった試料。その試料の `share` は null で描かない |
| `spots[]` | パネルごとに `spot_id` / `name` / `label`（表示名）/ `ontology` / `item`（展開元の項目）/ `adduct` / `mz` / `rt` / `msms`（MS/MS の裏付けの有無）と `groups[]` |
| `spots[].groups[].samples[]` | `sample` / `height`（PeakHeight。行が無ければ 0）/ `share`（分母に対する %。`height` 指定時と分母 0 の試料は null）/ `gap_filled` / `low_reliability` |
| `spots[].groups[].mean` / `sd` / `n_in_stats` | `stat` の定義での平均・SD と、その計算に使った試料数（低信頼でなく、`share` なら計算できた試料、`height` なら高さ > 0 の試料）。n = 1 なら SD は null、0 なら平均も null |
| `excluded` | 除いたスポットの `#<spot_id> <name>`。理由ごとに `internal_standard` / `manual` / `curation` / `auto_likely_wrong` / `standard_only`（`arf_plot_group_intensity` と同じ。§11.5）と、`require_msms=True` で MS/MS の裏付けがないために除いた `no_msms` |
| `caveats` | `arf_plot_group_intensity` と同じ種類に加え、`share_basis` が描く分子種を含まないとき（分母に含まれない分子種があり、割合が 100 % を超えうる） |

**割合の読み方**: `share` は**分母にしたスポットの高さの合計に対する %** で、試料の総量に対する割合ではない。
分母は `share_basis`（省略時は `items`）を同じ除外規則で解決したスポットで、`require_msms=True` なら分母にも効く。
分子種の割合は分母の選び方で変わるので、図を報告に載せるときは `share_basis` を併記すること。

### 11.7 `arf_pca_species` の返り値

選んだ分子種（`arf_plot_species` と同じ項目指定・除外）の PeakHeight 行列（試料 × 分子種）で PCA を回す read-only ツール。
ARF 経路のみ（読み込み済み `.arf` と兄弟 `.arf2` が要る）。既定（`output="image"`）はサーバ側で描いたスコア図の PNG と
説明（寄与率・`pc1_r_top` / `pc1_r_bottom`）。以下は `output="payload"` の契約。ファイルは
`save_figure(kind="pca", source="species")`。ローディング全量（`loadings_rows`）は戻り値に載せずセッションに残す。
入力の誤り（計算に使う試料が 3 未満・合計 0 の試料・未知の `orient_by` / `normalize` など）は
`{"status": "error", "message": ...}`。失敗した呼び出しの後は前回の結果を破棄する。

| フィールド | 意味 |
|------------|------|
| `points[]` | 試料ごとに `x`（PC1）/ `y`（PC2）/ `label`（試料名）/ `group`（群。2 群に入る試料は最初の群）/ `fitted`。`fitted=false` は低信頼の試料で、**主成分の計算に使わず投影しただけ**（スコアは計算に使った試料が決めた軸への射影） |
| `explained_variance_ratio` | 主成分ごとの寄与率（計算に使った試料での値。最大 5 成分） |
| `x_label` / `y_label` | 軸名（寄与率入り） |
| `n_fit` / `n_species` | 計算に使った試料数 / PCA に使った分子種数 |
| `dropped_zero_variance` | 計算に使う試料で値が変わらず外した分子種の数 |
| `excluded` / `caveats` | `arf_plot_species` と同じ（§11.6） |
| `pc1_r_top` / `pc1_r_bottom` | PC1 との相関 r が大きい・小さい分子種（表示名と r、各 5 件）。r は**計算に使った試料**での分子種と主成分スコアの相関で −1〜1 |
| `result_id` | この結果の ID（`save_figure` の `result_id` に使える） |

**読み方**: 値は `normalize`（`none` / `total`）→ log10(x + 1)（`log_transform`）→ 計算に使う試料の平均と母標準偏差で
autoscale の順。主成分の**符号は本来任意**で、`orient_by` を指定したときだけ、その群の平均スコアが正になるようにそろえる。
`normalize="total"` は試料量の差を除いて組成を比べるときの選択で、合計が 0 の試料があると計算しない。
