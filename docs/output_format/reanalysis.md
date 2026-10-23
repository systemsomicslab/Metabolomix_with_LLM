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
ハンドルへの `library_load` も同じ規則（`file_path` を省けばハンドルの極性の設定。逆の
`ion_mode` は `LIBRARY_POLARITY_MISMATCH`、設定が無ければ `LIBRARY_NOT_CONFIGURED`）。

**ハンドルは対話のデータフォルダを探さない。** スポットの照合（`library_match_feature` の
spot_id / precursor_mz・`dcl_spot_spectrum`・`dcl_diagnostic_ions`・`library_msms_qualify`・
証拠の束の MS/MS）は、ハンドルの `.arf2` / `.dcl`（か明示した `arf2_path` / `dcl_file`）だけを
使い、無ければ `status: "no_alignment_files"` を返す——別のデータセットの同じ spot_id の
スペクトルを黙って読まないため。

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

- KNN は R の `impute::impute.knn`（Bioconductor）の `knnimp` の規則を写したもので、`maxp` =
  特徴数（2 段の分割をしない）に当たる。特徴を行に、`scale="log2"` なら正の値の log2 の上で:
  欠測の割合が `rowmax` を超える特徴（重い行）は**重い行どうしの**列の平均で埋め、近傍にしない
  （重い行の中にその列の観測が無ければ軽い行の列の平均）。残りの特徴（軽い行）は軽い行だけで、
  欠測のある特徴ごとに、その特徴が観測された列のうち相手も観測された列での平方差の平均を距離に
  し（共通の列が無い相手は除く）、近い `k` 個を選ぶ（同距離は特徴の並び順）。欠測の列は近傍の
  うちその列を観測しているものの平均、全員が欠測なら**軽い行の**列の平均。`colmax` を超える
  試料があれば `RECIPE_PRECONDITION`。R の既定（`maxp=1500`）は特徴が多いと k-means で
  分割してから補完するので、特徴が 1500 を超えるデータでは R の既定の結果とは一致しない。
  **R の実出力との突き合わせはまだしていない**（手元に Bioconductor の `impute` が無い。規則を
  写し、手計算の例でテストしている）。
- LOWESS は statsmodels 0.15 の
  `statsmodels.nonparametric.smoothers_lowess.lowess(endog=QC 値, exog=QC の注入順,
  frac=span, it=iterations, delta=0, xvals=全試料の注入順)`（1 次の局所重み付き回帰・
  tricube の重み。`iterations=0` は頑健化の反復なし）。R の `loess`（既定は 2 次）とは別物。
  QC が少ないと当てはめが NaN になる点があり、その特徴は補正しない。

## 4. MS/MS の合否（msms-rules.v1）と診断イオン

`library_msms_qualify(rules, spot_id | precursor_mz, ...)`。照合は `library_match_feature` と
同じ（採点の芯は `library/matching.py`）。規則（正規化後。既定はこのとおり）:

```json
{"schema": "msms-rules.v1",
 "target": {"inchikey_block": "HVYWMOMLDIMFJA", "formula": null},
 "top_rank": 5, "min_fragments_excl_precursor": null, "precursor_exclusion_tol": null,
 "same_formula_margin": null, "lipid_class": null,
 "diagnostic_ions": [], "neutral_losses": [], "ions_match": "all"}
```

入力の `target` は `{"inchikey": <InChIKey 全体か先頭 14 文字>, "formula": <任意>}`。
設定した検査だけが `checks` に並ぶ（各行 `{check, status: pass|fail|not_evaluable, value, detail}`）:

| 検査 | 合格の条件 | `value` |
|---|---|---|
| `target_in_top_rank`（常に） | 総合スコアの順で `top_rank` 位までに、InChIKey 先頭ブロックが目標と一致する候補がある（最上位のそれを「目標の候補」） | 目標の候補の順位（無ければ null） |
| `min_fragments` | 目標の候補の、前駆体 ± `precursor_exclusion_tol`（省略時は採点の `ms2_tol`）を除いた一致の本数 ≥ 最小 | その本数 |
| `same_formula_margin` | 全候補のうち、目標と同じ分子式（`target.formula`、無ければ目標の候補の分子式）で先頭ブロックが違うものの最高点との差 ≥ `min_delta`（`score` 既定 `total_score`）。別物が無ければ pass | 差（目標 − 別物） |
| `lipid_class` | 1 位の候補の名前を GOSLIN で読んだクラス（正規化名の最初の空白の前）が一致 | 読めたクラス |
| `diagnostic_ions` / `neutral_losses` | 測定スペクトルに ±`tol`（既定 0.01 Da）のピークがあり、最大ピークに対する相対強度 ≥ `min_relative_intensity`。`ions_match` が `all` なら全部、`any` なら 1 つ。中性損失は前駆体 − `delta` を探す | 見つかった数 |

`verdict`: すべて pass → `qualified`、どれか fail → `not_qualified`、fail が無く判定できない
検査がある → `not_evaluable`。目標の候補が無いと、それに依る検査は `not_evaluable`
（`target_in_top_rank` が fail なので `verdict` は `not_qualified`）。MS/MS が無い・候補が無い
ときは検査をせずに `not_evaluable`（「断片が合わなかった」ではない）。**「明らかに別物か」は
判断しない**（再解析の計画の規則で決める）。`rule_hash` は正規化した規則の hash。

`dcl_diagnostic_ions(spot_ids, ions)` は同じ照会（`find_ions`）をスポットの一覧に当てる。
`ions` の各要素は `{"mz"}` か `{"neutral_loss"}`（任意で `label`・`tol`・
`min_relative_intensity`）。行の `status` が `ok` 以外（MS/MS が無い・索引の崩れ）のスポットは
照会していない——`found` が空なのは「イオンが無かった」ではない。

## 5. 同族列の RT（homolog-check.v1）

`dataset_homolog_check(spec)`。`members`（`feature_id`・`carbons`・`double_bonds`）か `auto`
（`lipid_class`・`evidence` 既定 `["sme"]`・`adduct`）のどちらか一方と `rule`
（`carbon_direction` 既定 `increasing`・`double_bond_direction` 既定 `decreasing`・
`max_violations` 既定 0。逆相の LC の向き）。

- `auto` は SME rank 1 の名前（`sml` を足せば SML の名前も）を GOSLIN で読み、分子種
  （`16:0_18:1`）は鎖を足して種のレベル（`34:1`）に直す。エーテル（`O-` / `P-`）は別の
  クラス（`PC O-`）として扱う。
- 同じ（炭素数, 二重結合数）が複数あれば RT の中央値で並べ、`duplicates` に挙げる。
- 同じ二重結合数の中で炭素数の昇順、同じ炭素数の中で二重結合数の昇順に並べた**隣どうし**で、
  規則の向きに反する対を違反と数える（RT が同じ対は違反にしない）。違反が
  `max_violations` 以下なら `pass`。違反に関わった特徴は `members` の `status` が `violation`。

## 6. 特徴どうしの関係

`dataset_feature_relations(feature_id, rt_window=0.05, mz_tol=0.01, min_r=0.7, adducts=None)`。
規則はアラインメントのキュレーションの関係（`curation/relations.py`）と同じ。

- `coeluting`: |ΔRT| ≤ `rt_window` の特徴（`drt` = 相手 − この特徴、`r` = 試料間の log1p 強度の
  Pearson 相関。値が 5 組未満なら空）。
- `as_ion_of`: この特徴が、相手の付加体（`adducts` の組）・同位体（M+1 / M+2）・源内断片
  （`-H2O` ほか）として説明できる関係。`ions_of_this`: 相手がこの特徴のイオンとして説明できる関係。
- 相手の付加体は SME rank 1 の付加体。**同定の無い相手には付加体の関係は出ない**（同位体・
  源内断片だけ）。MS-DIAL のリンク（`.arf2`）は使わないので、裏付け（`strong`）は相関
  `r ≥ min_r` と同位体の強度比だけで決まる。

## 7. 証拠の束（claim-evidence.v1）

`dataset_claim_evidence(spec, dataset_handle, output_path)`。1 主張について、オブザーバーの
機械的な仕事（走らせ直しと照合、候補の網羅の走査）を道具が 1 つの JSON にまとめる。

spec:

```json
{"claim_id": "C1",
 "target": {"inchikey": "...", "formula": "...", "rt_min": 8.4},
 "scan": {"ppm": [10, 20], "rt_tol_min": [0.3, 1.0], "adduct_tiers": null},
 "msms": {"top_n": 10, "rules": {"min_fragments_excl_precursor": 2}},
 "relations": {"rt_window": 0.05, "mz_tol": 0.01, "min_r": 0.7},
 "max_candidates": 50}
```

- `target` は `dataset_feature_search` の目標と同じ形（`target_id` は `claim_id`）。
- `scan` の `ppm` × `rt_tol_min` の全組で目標照合を走らせ直す（`adduct_tiers` を渡すと
  付加体も広げる）。省略時は目標の `ppm` / `rt_tol_min` の 1 組。
- `msms.rules` は msms-rules.v1 から `target` を除いたもの（target は主張の InChIKey・分子式）。
- `msms` / `relations` を省くと、その欄を作らない。

結果:

| 欄 | 意味 |
|---|---|
| `default_selection` | `target` をそのまま当てた目標照合（`dataset_feature_search` と同じ規則）。`status` が `resolved` なら `selected_feature_id` が規則で選んだ既定の特徴、決まらなければ `needs_input`（`reason`） |
| `scans` | 走査ごとの判定と、候補の特徴 ID |
| `candidates` | 全走査の候補の和（特徴の並び順）。`evidence`（`sme_match` など、走査で見えたもの）・`adducts`・`scans`（どの走査で見えたか）・`msms`（spot_id = SMF_ID の照合の上位 `top` と `qualify`。ライブラリが無ければ `{"status": "no_library"}`、MS/MS が無い・索引の崩れはその `status`）・`relations` |
| `truncated` | `max_candidates` で打ち切った候補の数 |
| `rule` / `rule_hash` | 正規化した spec とその hash |
| `input_fingerprint` | `dataset_hash`（データの指紋）・`library_sha256` |

**束は特徴を決めない。** 最終の特徴は再解析の実行役が override の規則の範囲で決め、
その特徴での統計は `dataset_analysis_run`（`features.feature_ids`）で別に出す。同じ spec・
同じ入力なら同じ JSON になる（オブザーバーの「走らせ直しと照合」は JSON の一致で済む）。

## 8. 質量精度と新しい知見の候補

`dataset_mass_accuracy(band_edges=[200], rank=1)`: SME の `rank` 位の候補の
(実測 − 理論)/理論 × 10⁶ を、理論 m/z の帯ごと（`[min_mz, max_mz)`）にまとめる。
`n`・`median`・`q1`・`q3`（`statistics.quantiles(n=4)` の exclusive 法）・`within_5ppm`・
`within_10ppm`（|ppm| ≤ 5・10 の割合）。n<2 の帯は `n` だけ。**許容幅は決めない**。

`dataset_novel_candidates(spec)`（novel-candidates.v1）:

```json
{"statistics": {"result_id": "res_..."},
 "ppm": {"default": 10.0, "bands": [{"max_mz": 200.0, "ppm": 10.0}]},
 "q_max": 0.05,
 "exclude": {"feature_ids": [], "inchikey_blocks": []}}
```

`statistics` は `{"result_id"}`（`dataset_differential` の 2 群比較。`group_a` が基準）か
`{"analysis_spec", "analysis_id"}`（analysis-spec.v1 とその解析）。SME rank 1 の同定があり、
|ppm| が帯の許容幅以内（`max_mz` 未満の最初の帯、無ければ `default`）で、q < `q_max` の特徴から
除外を外した一覧（q の昇順）。行は `feature_id name inchikey adduct mz rt ppm group_ref group_test
log2fc_test_over_ref p q direction`。論文の同定リストとの突き合わせはしない（再解析の側）。
