# 再解析の道具（ハンドル・spec・証拠の束）の出力定義

**再解析で使う**ツールの戻り値と、spec・結果の JSON の欄の意味（spec 2026-10-22）。
対話の使い方の出力は変わらない——既存のツールに足した引数（`dataset_handle` ほか）は、
渡さない限り戻り値を 1 文字も変えない。呼び出し連鎖は `docs/workflow/reanalysis.md`。

共通の約束:

- **spec は正規化して返す。** 既定値を埋め、未知のキーは拒否する（黙って無視しない）。
  結果には正規化した spec と `rule_hash` が入る。`rule_hash` は正規化した spec から
  パスの欄（`path` / `matrix_path`）を除いて作る——置き場所が違っても同じ規則なら同じ hash。
  入力のファイル（行列・群・共変量・注入順の TSV）は `input_fingerprint` に sha256 で入る。
- **同じ spec・同じ入力なら同じ JSON。** 時刻や乱数の ID は結果に入れない。
- **判定（再現したか・明らかに別物か）は返さない。** 道具が返すのは数字と、spec に書いた
  規則を機械的に当てた結果だけ。
- **研究室の参照ライブラリのパスは出さない。** ファイル名と sha256 だけ。

## 1. データセットのハンドル

ハンドル（`dh_` + 12 桁）は、1 つのサーバのプロセスにデータセット（とその極性の
ライブラリ）を複数並べて持つための ID。`dataset_handle_open` が返し、`dataset_*`・
`library_*`・再解析のツールの `dataset_handle` に渡す。

- **対話のセッションとは別の状態。** ハンドルを渡した呼び出しは、ハンドルを渡さない呼び出しが
  使う状態（`dataset_load` で読んだデータセット・`library_load` のライブラリ・直近の照合）を
  読みも書きもしない。閉じた・無いハンドルは `DATASET_HANDLE_NOT_FOUND`（開いている
  ハンドルの一覧つき）で、対話のセッションへ黙って落ちない。
- **ハンドルどうしも別の状態。** 同じ入力を開き直すと読み込み済みのもの（`reused_base`）を
  使い回すが、ハンドルには複製が渡る。片方の前処理はもう片方に見えない。入力のファイルの
  size か更新時刻が変われば読み直す。
- **並行と直列。** ハンドルを渡した呼び出しはサーバの別スレッドで回り、別のハンドルどうしは
  並行に進む。同じハンドルへの呼び出しは 1 本ずつ（到着順は保証しない）。
- **上限。** `LIPIDMIX_MAX_DATASET_HANDLES`（既定 8）。超えると `DATASET_HANDLE_LIMIT`。
  寿命は明示の `dataset_handle_close` か、サーバの終了まで。

`dataset_handle_open` の戻り値:

| 欄 | 意味 |
|---|---|
| `dataset_handle` | ハンドルの ID |
| `polarity` / `polarity_source` | データセットの極性と出所（`job_declared` / `adduct_majority` / `unknown`）。`dataset_feature_search` の目標照合と同じ決め方 |
| `n_features` / `n_samples` | mzTab-M の特徴数・試料数 |
| `detection_available` | 検出状態（gap-fill の区別）を取り込めたか |
| `library` | `{file, sha256, record_count}`。読まなかったときは `null` で、理由が `library_note` |
| `arf2` / `dcl` | スポットの照合（spot_id = SMF_ID）に使うファイル名。決まらなければ `null` |
| `reused_base` | 読み込み済みのものを使い回したか |

ライブラリは**データセットの極性の設定だけ**から引く（`[library] msp_positive` /
`msp_negative`、環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）。極性が決まらない・設定が
無いときはライブラリなしで開き、逆の極性のライブラリを掴まない。読み取り専用で開く。

## 2. 統計の spec（analysis-spec.v1）と結果（analysis-result.v1）

`dataset_analysis_run(spec, dataset_handle=None, output_path=None)`。論文の方法を spec に写し、
名前つきの解析（`analyses`）を並べる——感度解析は同じ spec に並べた別の解析。

```json
{"input": {"source": "dataset"},
 "groups": {"source": "metadata"},
 "comparison": {"reference": "HC", "test": "RA"},
 "injection_order": null,
 "covariates": null,
 "features": {"all": true},
 "aggregate": null,
 "analyses": [{"analysis_id": "main",
               "matrix": {"steps": ["detection_filter", "impute", "normalize", "drift_correct", "qc_rsd_filter"], "...": "§3"},
               "test": {"kind": "student", "transform": "log2", "covariates": []},
               "fold_change": "arithmetic_mean_ratio", "correction": "bh",
               "bh_population": "matrix", "thresholds": {"p": 0.05, "q": null, "fc": 1.5}}]}
```

| 欄 | 既定 | 意味 |
|---|---|---|
| `input` | `{"source": "dataset"}` | データセット（ハンドルか対話のセッション）か、外部の行列 `{"source": "tsv", "matrix_path", "orientation": "features_by_samples" \| "samples_by_features", "id_column"}`（較正用。論文の個体ごとの値など）。欠測は空・`NA`・`NaN`・`null` |
| `groups` | データセットなら `metadata` | `metadata` は実験情報シート（`group`・`role`・`include`）。TSV は `{"source": "tsv", "path", "sample_column": "sample", "group_column": "group", "role_column": null}`。外部の行列には TSV が要る |
| `comparison` | 必須 | `reference`（基準）と `test`（比較）の群の名前。比較名は `<test>_vs_<reference>` |
| `injection_order` | `null` | ドリフト補正の注入順。`{"source": "metadata"}` か TSV `{"source": "tsv", "path", "sample_column", "order_column": "order"}` |
| `covariates` | `null` | OLS の共変量の TSV `{"path", "sample_column": "sample", "missing": "drop" \| "error", "columns": [{"name", "type": "numeric"} \| {"name", "type": "categorical", "reference"}]}`。カテゴリは基準の水準に対する 0/1 の列に符号化する |
| `features` | `{"all": true}` | 報告する特徴（`{"feature_ids": [...]}` で絞る）。BH の母集団は `bh_population` が決める |
| `aggregate` | `null` | 1 本に集約して検定する。`{"kind": "sum", "label", "feature_ids": [...]}`・`{"kind": "sum", "label", "select": {"name_regex", "evidence": ["sme"]}}`（クラスの主張）・`{"kind": "ratio", "label", "numerator", "denominator"}`（比の主張） |

解析（`analyses[]`）:

| 欄 | 既定 | 意味 |
|---|---|---|
| `matrix` | `{"steps": []}`（前処理なし） | 行列の recipe（§3） |
| `test.kind` | 必須 | `student`（等分散の t）・`welch`・`mannwhitney`（両側）・`ols`（`transform(y) ~ 1 + 群 + 共変量`、群の係数の t の両側 p） |
| `test.transform` | `log2` | `log2` は有限の正値だけに当てる（0 以下は欠測。pseudocount も clip も無い）。`none` |
| `test.covariates` | `[]` | `ols` で使う共変量の名前（`covariates.columns` の部分集合） |
| `fold_change` | `arithmetic_mean_ratio` | `arithmetic_mean_ratio`（変換前の算術平均の比 test/ref）・`geometric_mean_ratio`（log2 の平均の差） |
| `correction` | `bh` | `bh`・`none` |
| `bh_population` | `matrix` | `matrix`（その解析の行列で検定できた全特徴）・`reported`（報告する特徴だけ）。検定できない特徴は母集団に入れない |
| `thresholds` | すべて `null` | `p`・`q` は「未満」、`fc` は `|log2fc| ≥ log2(fc)`。null の閾値は見ない |

**試料の照合は完全一致。** 群の TSV は行列の試料と過不足なく一致しなければ
`SAMPLE_MISMATCH`（欠け・余り・重複を名前で挙げる）。注入順・共変量の TSV は行列の試料を
すべて含めばよい（余りは無視）。群に入るのは役割が `sample` で `include` が真の試料だけで、
群の n が 2 未満なら `STATISTIC_GROUP_TOO_SMALL`。データセットの試料名は実験情報シートの
`sample_id`（無ければ assay 名）。

結果の行（解析ごと。TSV・JSON とも同じ欄）:

| 欄 | 意味 |
|---|---|
| `status` / `reason` | `tested`・`not_testable`（`insufficient_finite_values`・`zero_variance`・`rank_deficient`・`no_residual_degrees_of_freedom`・`test_statistic_not_finite`）・`not_in_matrix`（recipe の篩で落ちた） |
| `group_ref` / `group_test` | 群の名前（比較の向き） |
| `n_ref` / `n_test` | 検定に使った値の数（変換後に有限なもの。OLS は共変量の欠測で除いた後） |
| `mean_ref` / `mean_test` | 行列の値（変換前）の算術平均 |
| `fc_test_over_ref` / `log2fc_test_over_ref` | `fold_change` の定義の倍率と、その log2。**正なら test が高い** |
| `statistic` / `p` / `q` | t（student・welch・ols）か U（mannwhitney）、両側 p、BH の q |
| `beta_log2` | OLS の群の係数（変換後の尺度。log2 なら log2 の差） |
| `direction` | `test_higher` / `ref_higher`（OLS は `beta_log2`、ほかは FC の向き） |
| `meets_thresholds` | 与えた閾値をすべて満たすか（閾値が無ければ null）。**再現の判定ではない** |

解析ごとに `matrix`（`n_features_in`・`n_features_out`・段ごとの報告 `steps`）、
`bh_population_n`、集約なら `aggregate`（和: `members`・`summed`・`dropped`）、OLS なら
`covariate_excluded_samples`。結果の全体に `spec`（正規化後）・`rule_hash`・
`input_fingerprint`（`dataset_hash`・入力の TSV ごとの sha256）・`comparison`。

## 3. 行列の recipe（matrix-recipe.v2）

`steps` に挙げた段だけを、挙げた順に走らせる（挙げない段の設定は拒否する）。篩の段
（`detection_filter`・`qc_rsd_filter`）で落とした特徴はその直後に行列から外す——後の KNN の
近傍・PQN の参照に混ざらない。既定（正規化で埋まる）:

| 段 | 既定 | 定義 |
|---|---|---|
| `detection_filter` | `{"basis": "detected", "min_rate": {}}` | 役割（`qc`・`sample`・`blank`）ごとの検出率の最小。挙げた役割のすべてを満たす特徴を残す。`basis`: `detected` は gap-fill を除く実検出（検出状態が無ければ `DETECTION_UNAVAILABLE`）、`positive` は値が有限かつ > 0（gap-fill を含む） |
| `impute` | `{"method": "half_min", "missing": "zero", "k": 10, "scale": "log2", "rowmax": 0.5, "colmax": 0.8}` | 欠測（`missing`: `zero` は ≤0 か非有限、`nan` は非有限だけ）を埋める。`half_min`: 特徴ごとの正の値の最小の半分（正の値が無い特徴は欠測のまま）。`knn`: 下記。`none` |
| `normalize` | `{"method": "pqn", "reference": "qc"}` | `pqn`: 参照（`reference` の役割の試料の、特徴ごとの正の値の中央値）に対する商の中央値で割る。`tic`（有限値の和）・`median`（正の値の中央値）は係数を試料間の中央値で割る。係数が 0・非有限の試料は割らずに報告する |
| `drift_correct` | `{"method": "loess", "span": 0.75, "iterations": 0, "min_qc": 4}` | 特徴ごとに、正の QC 値を注入順に対して LOWESS で当て、`値 / 当てはめ × QC 値の中央値` に直す。正の QC 値が `min_qc` 未満、または当てはめが正の有限でない特徴は直さず数える。注入順はすべての試料に要る（`INJECTION_ORDER_MISSING`） |
| `qc_rsd_filter` | `{"max_rsd": 0.3}` | QC の有限値の `std(ddof=1)/mean` が `max_rsd` 以下の特徴を残す。RSD が出ない特徴も落とす（別に数える） |

**参照する標準の定義:**

- KNN は R の `impute::impute.knn`（Bioconductor）と同じ規則で、`maxp` = 特徴数（2 段の
  分割をしない）に当たる。特徴を行に、`scale="log2"` なら正の値の log2 の上で、欠測のある
  特徴ごとに、その特徴が観測された列のうち相手も観測された列での平方差の平均を距離にし
  （共通の列が無い相手は除く）、近い `k` 個を選ぶ（同距離は特徴の並び順）。欠測の列は近傍の
  うちその列を観測しているものの平均、全員が欠測なら列の平均。欠測の割合が `rowmax` を超える
  特徴は列の平均、`colmax` を超える試料があれば `RECIPE_PRECONDITION`。R の既定
  （`maxp=1500`）は特徴が多いと k-means で分割してから補完するので、特徴が 1500 を超える
  データでは R の既定の結果とは一致しない（`maxp=nrow(data)` を渡した R とは同じ規則）。
- LOWESS は statsmodels 0.15 の
  `statsmodels.nonparametric.smoothers_lowess.lowess(endog=QC 値, exog=QC の注入順,
  frac=span, it=iterations, delta=0, xvals=全試料の注入順)`（1 次の局所重み付き回帰・
  tricube の重み。`iterations=0` は頑健化の反復なし）。R の `loess`（既定は 2 次）とは別物。
  QC が少ないと当てはめが NaN になる点があり、その特徴は補正しない。
