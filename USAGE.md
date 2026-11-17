# USAGE — ms-data-parser MCP ツール一覧(全58ツール)

MS-DIAL 出力(`.arf` / `.arf2` / `.dcl` / `.EIC.aef`)と mzTab-M を解析し、PCA・差次的解析・
アノテーション検証・レポート記録までを行う MCP サーバーのツール一覧。各ツールの引数・戻り値・断る条件は
ツールの説明(docstring)、出力フィールドの意味は [docs/output_format/](docs/output_format/core.md)(`docs_read`)、
内部の呼び出し順は [docs/workflow/](docs/workflow/index.md) にある。

入口は 2 つで、どちらも MS-DIAL の解析フォルダに合流する。フォルダに MS-DIAL の出力(`.arf` 等)があれば
`load_dataset`。生データだけなら `pipeline_run` が MS-DIAL Console を回す(極性ごとに数十分。run のフォルダを書く)。

```
list_data_files(directory) → load_dataset(directory=...)
   (戻り値の `caveat` は意味論の要約。未読なら 1 回だけ載る。`assay_kind` が未確定のうちは「未確定」の行つき)
   → record_objective(assay_kind=...)(種別は根拠から決めて記録する。読み込み済みのデータセットが要る)
   → dataset_overview()(品質 → 全体 → クラス → 分子種。各段の `next` から掘る)
   → dataset_status / dataset_preprocess
   → dataset_pca / dataset_compare → dataset_plot(`output_path` で保存も) → save_figure(`reports/figures/` へ)
   → write_report

生データから:
pipeline_run → pipeline_status → load_dataset(pipeline_path=..., polarity=...) → dataset_preprocess
   → dataset_pca / dataset_compare → dataset_export_differential
   (訂正・再開は pipeline_resume)
```

---

## 1. データ探索・ロード(入口)

| ツール | 機能 |
|--------|------|
| `list_data_files` | フォルダ(`directory`。必須)のファイルを拡張子別に一覧する。既定は解析できる拡張子(.arf/.arf2/.pai2/.dcl/.EIC.aef/.mddata/.mdproject と、MS-DIAL を回した形跡の `.mzTab`・`analysis-job.json`)だけ。`all_files=True` で生データも並べ、raw ごとの取得時刻・取得順・40 分を超える間隔の数と、`sample-manifest.tsv` があれば群と取得順の重なり(`complete` = 範囲が重ならず測定順と交絡 / `independent` / `partial` = それ以外 / `unknown`)と群ごとの取得順の範囲を返す。`.d` と Waters の `.raw` はフォルダが 1 検体(末尾に `/`。ほかの拡張子の同名フォルダは数えない)。`extension` で絞る。次は `load_dataset(directory=...)`、生データだけなら `pipeline_run`。 |
| `load_dataset` | データセットを読み、新しいハンドル(`dh_…`)を開いて直近にする。入口は `directory`(解析フォルダか中の `.arf` / `.arf2`)・`mztab_path`・`job_path`・`pipeline_path` のどれか 1 つ。Console の実行後は `pipeline_path`(`pipeline_status` の `receipt.next` がそのままの呼び出し。2 極性の run で `polarity` を省くと `needs_input`(`field=polarity`))、既存のジョブは `job_path`(ジョブの宣言で正準 mzTab-M を 1 つに選び、曖昧なら止まる。生成物のパスは記録の `root` から復元する)。run の試料の対応表だけを別に載せる入口は無い。主な引数は `source`(`auto` / `folder` / `mztab`)・`measure`・`allow_incomplete`(読んだデータセットは探索専用で、警告の先頭に中断された実行の生成物である旨が入る)・`library`(`auto` / `none`)・`arf2_path`・`overview`(`dataset_overview` と同じ 4 段から `counts` を省いたもの)・`wait_seconds`(既定 20 秒。終わらなければ `{"status": "loading", "load_id"}` を返し、結果は `dataset_load_status`)。戻り値と注記の意味は [docs/output_format/dataset.md](docs/output_format/dataset.md)。次は `record_objective(assay_kind=...)` と `dataset_overview()`。 |
| `dataset_load_status` | `load_dataset` の受付番号(`ld_…`)の読み込みを引く(`load_id`。省けばまだ結果を渡していない 1 本)。`wait_seconds`(既定 0。最長 30 秒)だけ待つ。読み込み中は `{"status": "loading", ...}`、終われば `load_dataset` の戻り値そのもの。知らない番号は `DATASET_LOAD_NOT_FOUND`。 |
| `job_status` | 長い道具(`pathway_project` / `dataset_homolog_check`)の受付番号(`job_…`)を引く(`job_id`・`wait_seconds` 既定 0、最長 30 秒)。走っている間は `{"status": "running", ...}`、終わればその道具の戻り値。仕事の間、同じハンドルへのほかの呼び出し(`dataset_close` も)は最長 25 秒待って `HANDLE_BUSY` を返すので、このツールで終わりを待ってから呼び直す。知らない番号は `JOB_NOT_FOUND`。 |

## 2. MS/MS(`.dcl`)

| ツール | 機能 |
|--------|------|
| `dcl_find_msms` | precursor m/z(任意で RT)に一致する MS/MS を `.dcl` から引く(アノテーションの裏取り)。`file_path` を省けば直近のデータセット(`dataset_handle`)の `.dcl`。該当ゼロ(`not_found`)は「MS/MS 未取得」で、「フラグメント不在」ではない。 |
| `dcl_spot_spectrum` | スポット ID(mzTab-M の SMF_ID と同じ)で代表 MS/MS の全ピークを引く(`library_match_feature(spot_id=…)` と同じ規則で `.arf2` と照らす)。`arf2_path` / `dcl_file` を省けばデータセットのもの。 |
| `dcl_diagnostic_ions` | スポットの代表 MS/MS(`spot_ids`)に診断イオン・中性損失(`ions`)があるかを照会する。MS/MS が無いスポットは `status` で示し、「イオンが無かった」と区別する。 |

## 3. 実験目的の記録

| ツール | 機能 |
|--------|------|
| `record_objective` | 実験目的を `reports/<analysis_id>.objective.md` に作る・直す(渡した欄だけ。欠けた欄は空のまま作り、戻り値に挙げる)。創発した小問は `add_subquestions` で追記。`assay_kind`(`lipid` / `metabolite` / `unknown`)はハンドルごとで、渡すとそのハンドルの解釈規則が切り替わり、目的ファイルにも残る(同じ解析フォルダを開き直すと読み戻す)。置き場所は `[paths] reports_dir`、無ければデータセットの解析フォルダ配下 `reports/`(どちらも無ければ `missing_state`)。同じ `analysis_id` の目的が別の解析フォルダのものなどのときの置き場所はツールの説明。 |

## 4. 図の保存・レポート

| ツール | 機能 |
|--------|------|
| `dataset_plot` | データセットの図を描く(既定は PNG と 1 行の説明。`output="payload"` で点列 JSON、`output_path` で同じ図をファイルにも)。`kind`: `volcano`(直近の差次的結果。対比が複数なら `contrast_id`)/ `pca` / `loadings` / `group_intensity` / `species` / `species_pca`(分子種の PCA。結果を `result_id` で登録する。試料の情報の更新・`dataset_exclude` で古くなり、キュレーションの更新では古くならない)/ `eic` / `eic_compounds`。kind ごとの引数はツールの説明。古くなった結果は `STALE_ANALYSIS_RESULT`。次は `save_figure`。 |
| `save_figure` | `dataset_plot` が kind ごとに最後に作った図を `reports/figures/<analysis_id>_<kind>.png` に保存する(置き場所は `write_report` と同じ)。`result_id`(`pca` / `volcano` / `loadings`)・`contrast_id`(`volcano`)で結果を名指しして描き直せる。古い図は `STALE_FIGURE`、描いたことの無い図は `missing_state`、未知の `kind` は `{"status":"error"}`。戻り値は保存した絶対パスと、`write_report` の本文に埋め込む相対パス。 |
| `write_report` | 解析・解釈レポート(Markdown 本文。見出しは自由)を `reports/<analysis_id>.md` に上書き保存する(置き場所は `record_objective` と同じ)。 |
| `read_report` | 過去レポートを読み戻す(最新更新のもの)。`analysis_id` を省くと 1 行索引(analysis_id / date / status)。 |

## 5. MS-DIAL Console のメソッドファイル

Console の計画・実行・状態は §8 の pipeline が担う。ここは pipeline の要求の `method_file` に渡すメソッドファイルを用意する 2 ツール。

| ツール | 機能 |
|--------|------|
| `console_method_template` | 既存のパラメータから別極性用のメソッドファイルを作る(`out_path`・`polarity`・`based_on`・`dataset_root`・`omics`・`library_mode`・`msp_file`)。その極性で一度でも GUI 実行があれば、`method_file` を省いた `pipeline_plan` がそのメソッドを見つけるので要らない。差し替えるのは `Ion mode` と `Searched adduct ions` だけで、空の `Lbm file path` も埋める(解決順は `pipeline_plan` と同じで、ビルド生成物が `[library] lbm` より先)。相対で宣言されたパスは元ファイル基準の絶対パスに書き換える(非 ASCII なら `METHOD_ENCODING_UNSUPPORTED`)。 |
| `console_method_candidates` | 使えるメソッドファイルの候補を並べる(`dataset_root`・`polarity`・`omics`・`search_dirs`)。`pipeline_plan` が失敗する前に呼べる(手動実行 UI の「メソッドを選ぶ」画面もこれを呼ぶ)。探すのは `dataset_root` 直下 → 兄弟フォルダ直下 → `runs/*/analysis-job.json` の `software.method_file`(再帰しない。pipeline の run のジョブは 1 段深いので拾わない)。候補ごとに `origin`(`same_dir` / `sibling` / `past_run` / `given`)・`usable`(`direct` / `needs_polarity_conversion`。後者は `console_method_template` を経由する)・`key_params`。 |

## 6. DatasetState 解析(mzTab-M 経路)

読み込んだデータセット(解析フォルダも mzTab-M も同じ `DatasetState`)に前処理・PCA・差次的解析を回す。
流れは `dataset_preprocess`(行列を作り `matrix_id` を返す。省略可)→ `dataset_compare` / `dataset_pca`。
`matrix_id` を省くと、そのハンドルで最後に作った行列、無ければ既定の前処理の行列を使う(戻り値の `matrix_source`)。
判定は `changed` / `unchanged` / `inconclusive` で、「有意でない」は「変わっていない」ではない。
欄の意味は [docs/output_format/dataset.md](docs/output_format/dataset.md)。

| ツール | 機能 |
|--------|------|
| `dataset_status` | データセットの状態を返す。`view`: `summary`(既定。群ごとの件数 `groups` は `dataset_compare` の群名、役割・タグ・除外の件数、検出状態 `detection`、読み込みの注記 `load_notes`、`assay_kind`。解析フォルダのデータセットの `validation` は `{"applicable": false, ...}`)/ `samples`(試料の表。`query` の因子トークンで絞る。§7)/ `tags` / `results`(ハンドルの行列と結果の一覧)。`detection` の `available=false` の間は非ゼロ値に補間値が混ざるので、検出率・欠測率は数えられない。 |
| `dataset_exclude` | 試料・特徴を解析から外す・戻す(可逆。`samples`・`features`・`mode` = `add` / `remove` / `clear`・`reason`)。知らない名前・ID が 1 つでもあれば何も変えずに拒む。除外が変われば登録した行列と結果は古くなる(`dataset_preprocess` からやり直す)。 |
| `dataset_preprocess` | 前処理(選択 → ドリフト補正 → 品質の足切り → 正規化 → log2 → 反復のまとめ)を回し、行列を登録して `matrix_id` を返す。主な引数は `drift_mode`(`auto` / `on` / `off`)・足切り(`blank_min_ratio` / `detection_min_rate` / `max_qc_rsd` / `max_d_ratio`。`None` で行わない)・`normalize`(`auto` / `none` / `tic` / `median` / `pqn`)・`normalization_features`・`collapse_replicates`・`group_levels` / `group_factors`(§7)。次は `matrix_id` を `dataset_compare` / `dataset_pca` / `pathway_project` に渡す。列の意味は `docs/output_format/dataset.md` §10.2。 |
| `dataset_compare` | 統計エンジン(analysis-spec v2)で 2 群・多群・2 因子を比べ、要約を返す。群の指定はちょうど 1 つ: `reference_group`+`test_group` / `groups`(3 群以上)/ `comparisons`(と `factors`)。主な引数は `matrix_id`・`covariates`・`block`・`p_adjust`・`bh_population`・`allow_confounded`(既定 true)・`top_n`・`output_path`(全列の TSV)。判定は q と CI で決まり、`log2fc_threshold` は判定に使わない。全行は `result_id` で保持し、直近の結果になる。次は `dataset_plot(kind="volcano")`・`dataset_feature_search`・`dataset_export_differential`。 |
| `dataset_overview` | 最初の概観。品質 → 全体 → クラス → 分子種の 4 段で、群の間で何が変わり何が変わっていないかを粗く示す(既定の前処理・PCA・moderated t)。2 群の向きは仮で、決めるのは `dataset_compare`。事実には段ごとの ID(`Q1` `G1` `C1` `S1`、外からは `<overview_id>#Q1`)が、各段には掘り下げの `next` が付く。比較と PCA は `result_ids` で掘れるが直近の結果にはしない。`reports/stats/` に step=overview の記録を 1 件追記する(既定の前処理の行列を新しく作ったときは、その記録も)。引数は `top_n`(既定 5、上限 20)と `dataset_handle`。欄の意味は `docs/output_format/dataset.md` の「最初の概観」の節。 |
| `dataset_pca` | 前処理済みの行列で PCA(stats-pca.v1)。軸は生体の個体で決め、QC は同じ軸へ投影する。主な引数は `matrix_id`・`scaling`(`autoscale` / `pareto`)・`n_components`(既定 5)・`impute`。`group_levels` / `group_factors` は色分けと記録だけに使う。次は `dataset_plot(kind="pca")` / `"loadings"`。 |
| `dataset_export_differential` | 差次的結果(`result_id`。省けば直近。対比が複数なら `contrast_id`)を InChIKey 付きの 1 ファイルに書き出す。差次的エクスポートの契約(15 列 + `contract_version` メタ行。`verdict` は列に出ない)なので、下流のパスウェイ解析にそのまま渡せる。`q_value` はいつも BH の q。`apply_curation`(既定 true)でキュレーションの判断を当て、`exclude_feature_ids` で検定の後に行を落とす。古い結果は書き出さない(`STALE_ANALYSIS_RESULT`)。 |
| `dataset_feature_search` | 特徴を探し、目的の物質に当たる特徴を規則で判定する(読み取り専用)。窓検索(`feature_ids`・`min_mz`/`max_mz`・`min_rt`/`max_rt`・`sort`)か、目標照合(`targets`。1 回 50 件まで。各要素は `target_id`・`inchikey` か `lipid_species`・`formula` か `expected_mz`+`adduct`)。目標照合は qualified がちょうど 1 つのときだけ `resolved` で、ほかは `needs_input`(強度・ppm・RT の近さで選ばない)。極性が決まらず付加体の層も無ければ `POLARITY_UNKNOWN`。差次的結果があれば統計を付ける(`result_ids`・`contrast_id`・`include_statistics`・`with_group_labels`・`merge_split`)。照合の窓は `ppm_offset`・`rt_offset_min` でずらせる。欄の意味は `docs_read(topic="mztab")`。次は `library_match_feature(spot_id=…)` / `dataset_plot(kind="eic")`。 |
| `dataset_set_sample_metadata` | 実験情報シート(sample-manifest、TSV)を全件検証してから一括で当て、群・バッチ・注入順・役割などを直す(Console はやり直さない。不備があれば何も変えない)。v2 のシートは `subject_id` と因子の列を足せる。schema の行の直後に宣言の行 `# injection_order_source = <raw_acquisition_time \| filename \| sheet>` を置ける(ほかの `#` の行は `SAMPLE_MANIFEST_INVALID`)。試料の表が変われば(因子の列だけでも。`changed_fields`)行列と結果は古くなる。陽陰を統合したデータセットでは、シートは陽イオンの測定の `batch`・`injection_order` だけを決める(陰イオンの直し方はツールの説明)。次は `dataset_preprocess`。 |

## 7. Class ID に無い因子で絞る・比べる

MS-DIAL の Class ID は入力された 1 文字列で、時点や複製はサンプル名にしか無いことがある(例: Class ID = `ILG`/`control` だけで、時点 `6h` はサンプル名のみ)。
`dataset_status(view="samples", query=...)` の因子トークンは Class ID とサンプル名の両方から解くので、そのまま書ける
(`_` 区切りは順不同の AND、`,` 区切りは OR。一致しなければ使えるトークンの語彙を添えて拒む)。

```python
dataset_status(view="samples", query="ILG_6h,control_6h")   # 6h の 2 群の試料だけを表で見る
dataset_preprocess(normalize="median", group_levels=["control_6h", "ILG_6h"])   # 名前のトークンで群を解く。matrix_id を返す
dataset_compare(reference_group="control_6h", test_group="ILG_6h", matrix_id=<返ってきた matrix_id>)  # 時点を揃えた2群比較
```

行列の群は、既定では試料の `group`(Class ID)。因子トークンで群を作るには `dataset_preprocess` の `group_levels` / `group_factors`。
`dataset_compare` に `matrix_id` を渡すと群は行列が決める(群の引数を一緒に渡すと、行列と同じ分け方なら `notes` に書いて進み、違えば `needs_input`)。

## 8. 生データフォルダからの統括(pipeline。生データ → MS-DIAL Console → 解析フォルダ → 受領)

GUI で解析済みのフォルダを読む入口(§1 の `load_dataset`)と対等な、もう一方の入口。生データフォルダ 1 つから、
raw を run の中(`<run>/<極性>/input/`)へ置き、極性ごとに MS-DIAL Console を回し、できた解析フォルダを mzTab-M と
照合して受領(`results/receipt.json`)を書くまでを 1 本の永続 run(`pipeline-run.json`。既定の置き場所は
`<dataset_root>/runs/pipeline_<UUID>/`)として進める。統計は計算しない(続きは `load_dataset(pipeline_path=..., polarity=...)`)。
同じ入力からは同じ結果になる(fingerprint と `request_id` による冪等性)。メソッドファイルは §5。
内部の呼び出し順は [docs/workflow/pipeline.md](docs/workflow/pipeline.md)。

| ツール | 機能 |
|--------|------|
| `pipeline_run` | 生データフォルダ(`dataset_root`)から解析フォルダを作る。入力検査・計画の保存・worker の起動までを行い、すぐ `pipeline_path` を返す(工程は裏で進む)。引数は `pipeline_plan` と同じ(`request` = `pipeline-request.v3`・`request_id`)。2 極性は要求の順に 1 極性ずつ回し、片方が失敗しても続ける(run は `partial`)。Console は同じ設定ファイルを読むサーバ・worker の間で同時に 1 本(排他のロックは設定ファイルのフォルダの `data/.msdial-console.lock`。別のクローンも `LIPIDMIX_CONFIG` で同じ設定ファイルを指せば 1 本になる)。既存の run を再利用したときは起動しない(`launched=false`。進めるなら `pipeline_resume`)。次は `pipeline_status`。 |
| `pipeline_plan` | `pipeline_run` と同じ入力検査と計画の保存だけを行い、Console を起動しない。戻り値の `resolved.units.<極性>` で実行条件(method・lbm・msp・polarity・library_mode・alignment_reference・`dot_product_cutoffs`。空にした宣言は `removed_declarations`)を確かめる。Console の実行体の起動検査も同じで、未設定は `MSDIAL_EXE_NOT_FOUND`(`human_action_required` / `how_to_set`)、Console でなければ `MSDIAL_EXE_NOT_CONSOLE`、`--help` の時間切れは `MSDIAL_EXE_PROBE_TIMEOUT`。`library_mode`: `msp_only` は研究室 MSP だけ。`msp_then_lbm`(脂質の run)は研究室 MSP を優先し(注釈器設定表 `inputs/msp-annotators.tsv` に priority 2 で 1 行)、MSP に MS/MS で照合した結果が無いものだけを LBM で補う(LBM は `lipidmix.local.toml` の `[library] lbm` だけから。旧 `[msdial] lbm` が残っていれば `CONFIG_INVALID`、無ければ `LBM_NOT_CONFIGURED`、指す先が無ければ `LBM_NOT_FOUND`)。ライブラリの sha256 は計画に固定し、変われば別の run になる(`pipeline_resume` は `INPUT_CHANGED` で止まる)。`alignment_reference`(既定 `"median_injection"`)は MS-DIAL の基準試料を注入順の中央の試料にする(`"method"` は上書きしない)。 `method_overrides`(単位ごと。`{Console のキー: 値}`)はメソッドの原本を変えずに実効メソッドの値を上書きする(例 `{"Minimum peak height": 50000, "N percent detected in one group": 50}`。原本にある別名も同じ値にする)。キーと型は Console の `ConfigParser`(MsdialWorkbench afd5f95)で検査し、知らないキー(近い綴りを `suggestions` に)・現用の Console のビルドが読まないキー(`Number of threads`・`Lbm annotator priority`・`Detailed alignment provenance`・`Annotation candidates` とその別名)・型に合わない値(`Minimum peak height` は整数)・ほかの欄が受け持つキー(`Ion mode`・ライブラリやファイルのパス。`Alignment reference file ID` は `alignment_reference="method"` のときだけ)・計画が自分で書くキーと重なるものは `PIPELINE_REQUEST_INVALID`。上書きが変われば別の run になる(resume では変えられない)。`resolved.units.<極性>.method` は結果を大きく変えるキーの実効値 `key_params`・宣言の無いキー `key_params_undeclared`・`user_overrides`(`from` / `to`)も持つ。 |
| `pipeline_status` | run の状態を読む(読み取り専用。`pipeline_path`・`include_details`。worker は呼ばなくても最後まで進める)。`status`(`completed` / `partial` / `failed` / `needs_input` / `running` ほか)・極性ごとの `units`・`errors`・`warnings`、受領の後は `receipt`(次に呼ぶ `load_dataset(...)` の `next`、極性ごとの実効メソッドの主なキー `method_params`: `key_params` / `undeclared` / `overrides`)。mzTab-M と照合できなかった極性は警告 `ANALYSIS_FOLDER_UNMATCHED` で run は `partial`(`load_dataset(..., source="mztab")` なら読める)。`save_project=true` で GUI プロジェクトが無いのは `completed` のまま警告 `GUI_PROJECT_UNAVAILABLE`。worker が消えていれば `recovery_hint` が `pipeline_resume` を示す。次は `completed` / `partial` なら `load_dataset(pipeline_path=..., polarity=...)`、`needs_input` なら `pipeline_resume`。 |
| `pipeline_resume` | 試料の対応表の訂正(`updates` は `sample_manifest` か `units: {<極性>: {"sample_manifest": ...}}` だけ。受領だけを書き直す)と、止まった工程からの再開。上流の条件を変えるなら新しい run(`NEW_PIPELINE_REQUIRED`)。Console をやり直すのは `rerun_upstream=True` のときだけ(自動再試行しない)。`request_id` で冪等。 |
| `pipeline_sample_sheet` | ISA(`s_*.txt`+`a_*.txt`)か mwTab の試料の表から実験情報シート(`sample-manifest.tsv`。sample-manifest.v2)を raw のフォルダに書く(`raw_dir`・`sample_table`・`factor_columns`・`assay_table`・`qc_selector` / `blank_selector`・`overwrite`)。群は因子の列の値の組で、決まらない生体試料は群 `unassigned`・`include=false`。どれも決まらなければ書かずに `SAMPLE_TABLE_INVALID`(`details.files_without_table_row`)。注入順は raw の取得時刻の順。`unreadable_reasons` の理由にファイル名は入れない。既にあれば `SAMPLE_SHEET_EXISTS`。次は `pipeline_run`(このシートを既定名で読む)。 |
| `pipeline_cancel` | 取消要求を保存する(`pipeline_path`。冪等)。生データは消さない。止まったかは `pipeline_status` で確かめる。 |

## 9. 参照ライブラリ(MS/MS スペクトル照合)

`.dcl` の測定 MS/MS を参照ライブラリ(`*_Loaded.msp2.dbs` か `*.msp`)と照合し、対向プロット(mirror plot)で確かめる。
`library_load` → `library_match_feature` → `library_plot_mirror` の順に呼ぶ(前段の状態が無ければ `missing_state`)。

| ツール | 機能 |
|--------|------|
| `library_load` | 参照ライブラリを開き、照合用の SQLite store を作る(キャッシュがあれば使う)。`file_path`(`.msp` / `.dbs`)を省けば、データセットが同定に使ったライブラリの記憶(`library_record`。sha256 を測り直す)を開く(データフォルダも設定の極性の MSP も探さない)。ライブラリはハンドルごとに 1 つ。キャッシュは `[paths] library_cache`(既定は設定ファイルのフォルダの `data/.library-cache`)。大きな `.msp` の初回構築は `python -m metabolomix.library.store --file <path>`(または `--ion-mode <極性>`)で前もって済ませられる([docs/cli.md](docs/cli.md))。戻り値のライブラリはファイル名と sha256 だけ。 |
| `library_match_feature` | 測定 MS/MS(`.dcl`)を候補と照合し、上位を `total_score`(MS-DIAL の `GetTotalScore`。正規化しない和で 1 を超える)の順に TSV で返す。対象は `precursor_mz` か `spot_id`(`.arf2` から precursor・RT・極性を取り、そのスポットの代表 MS/MS を使う)のどちらか一方。上流の順位付けを全部は再現しないので、上位に化学的にありえない候補が残りうる。`mz_tol_ppm` は窓に重ねる絞り込み、`detail=true` で一致したピークも返す。未注釈スポットの確認は `dataset_feature_search` → `library_match_feature(spot_id=…)` → `library_plot_mirror`。 |
| `library_search` | 読み込んだライブラリを名前(部分一致)・分子式(完全一致)・InChIKey(先頭 14 文字)で引く(読み取りだけ。`limit` 件まで TSV)。 |
| `library_plot_mirror` | 直近の照合から対向プロット(上段 = 測定・下段 = 参照)を描く。`spot_id` で照合し直し、`candidate` で候補を選べる。`scale`(`relative` / `sqrt` / `log10`。precursor がベースピークのスペクトルは `relative` だと診断イオンが潰れる)・`label_policy`(`auto` / `msdial`)・`output`(`image` / `payload`)。採点に入らなかった測定ピークは灰色で描く。 |

## 10. アラインメントのキュレーション(注釈の一覧確認と機械判別)

注釈付きスポットを EIC・対向プロット・Δppm・ΔRT・RT–m/z 傾向と機械判別つきで一覧し、人が付けたフラグを記録する。
順は `load_dataset`(解析フォルダ)→ `library_load`(アラインメントと同じフォルダの `*_Loaded.msp2.dbs` を推奨)→
`curation_review` →(`curation_suggest`)。5 つのツールは `dataset_handle`(省けば直近)を持ち、`file_path` は無い。
`.arf2`・EIC・`.dcl`・ライブラリはそのデータセットのものを使い、`data/` は探さない(`.arf2` の無いデータセットは理由つきのエラー)。

人がビューア(`html_path` をブラウザで開く)でフラグを付けて **Submit** を押すと、ビューアが直接記録し `_tags.xml` に反映する
(LLM を介さない。受け口は `127.0.0.1` で、Finish か 30 分の無通信で閉じる)。受け口が無いとき(時間切れ・サーバ再起動)は
「Copy submission text」のテキストをチャットに貼り、`curation_submit` に渡す。自動判別の `likely_wrong` はビューアで赤の破線枠になり、
最初から「間違い」が選ばれ、メモ欄に判定根拠が入る(suspect もメモ欄に根拠)。フラグの無いスポットは「間違っていない」で、何も記録しない。

| ツール | 機能 |
|--------|------|
| `curation_review` | 対象(既定は注釈付き全部。`ontology`・`name_contains`)の証拠収集・機械判別・HTML ビューアの生成を 1 回で行う。戻り値は suspect 以上とフラグ済みのスポットの TSV(先頭 `max_rows` 行)・判定の件数・クラス別の傾向・`html_path`・`submit`。判定は `likely_wrong` / `suspect` / `ok`(理由のコードは [docs/output_format/curation.md](docs/output_format/curation.md))。`thresholds` で閾値を上書きできる。 |
| `curation_suggest` | 間違いになったスポット(既定は wrong フラグ+元レビューの likely_wrong)と未注釈スポット(`unannotated`)に注釈の候補を並べる。先に `library_load` と `curation_review`(元レビューは `review_id`、省けば最新)。人がビューアで選んで送信する(受け口が無ければ送信用テキスト → `curation_submit`)。 |
| `curation_submit` | フラグを追記する。ビューアの送信用テキストを `submission_text` に渡すか、`review_id` + `flags`(`wrong` / `suspect` / `confirmed` / `clear`。候補付け `cs-…` は `assign` / `redundant` / `clear`)。`source` は判断の出どころ(`user` 既定 / `llm`)。不正な要素が 1 つでもあれば何も書かない。記録の後 `_tags.xml` に反映する(assign / redundant は変えない。MS-DIAL でプロジェクトを開いたままだと GUI の保存で上書きされる)。 |
| `curation_flags` | データセットの `.arf2` で有効な判断を TSV で返す(spot_id・flag・name・of・note・source・ts)。 |
| `curation_view_data` | ビューア(MCP Apps)専用。LLM は呼ばない。 |

フラグは `<アラインメントのフォルダ>\curation\flags.jsonl` に、`.arf2` の sha256 と MasterAlignmentID の組で記録する
(MS-DIAL を再実行すると古いフラグは当たらず、件数を `warnings` で知らせる)。`dataset_export_differential` は既定(`apply_curation=True`)で
`wrong` と `redundant` のスポットを外し、`assign` のスポットの同定を記録した値に置き換える(15 列の契約は変えない)。
`flags.jsonl` に読めない行があれば、どのツールもファイル名と行番号を返して止まる。
MCP Apps で会話内にビューアを出す経路は、Claude Desktop のローカル stdio サーバでは未検証(ブラウザで `html_path` を開く経路が主)。

## 11. パスウェイ投影(MSLipidMapper)

前処理済みの脂質データを MSLipidMapper で remodeling / ceramide / global のパスウェイへ投影する。MSLipidMapper は
研究室サーバの REST API(`[mslipidmapper] mslipidmapper_api_url`)を先に使い、使えなければローカルの Rscript
(`[mslipidmapper] rscript`)に回す。どちらで動いたかは戻り値の `backend`(サーバのジョブの記録と失敗の扱いは
[docs/output_format/pathway.md](docs/output_format/pathway.md))。どちらも使えないときはサーバ側のコード
(`MSLM_REMOTE_UNAVAILABLE` など。`details.local` にローカルが使えない理由)。

**研究室サーバを使う**(R は要らない): `lipidmix.local.toml` の `[mslipidmapper]` に `mslipidmapper_api_url`
(例 `'http://<サーバ>:8000'`)を書く(環境変数では指せない)。認証が要るときだけ `api_token`。
中間 mzTab-M(試料名・群名・脂質名・強度)はサーバへ送られる。

**ローカルの R を使う**ときの初回のセットアップ:

1. R の 4.6 系を入れる。
2. リポジトリ直下で `Rscript r/setup.R` を 1 度実行する(MSLipidMapper と依存を CRAN / Bioconductor / GitHub から
   renv の仮想環境 `r/` にだけ入れる)。未実行なら `MSLM_R_PACKAGES_MISSING`、入っているのに読み込めなければ
   `MSLM_R_PACKAGE_LOAD_FAILED`。`server_update` で `r/renv.lock` が変わったときも再実行する。
3. `lipidmix.example.toml` を `lipidmix.local.toml` に写し、`[mslipidmapper]` の `rscript`(`Rscript.exe` のパス)を書く。
   `browser`(Edge / Chrome)は標準の場所なら自動で見つける。環境変数では指せない。設定は呼ぶたびに読む(サーバの再起動は要らない)。
4. ローカルの `pathway_plot` は Edge か Chrome と、unpkg.com への HTTPS 接続が要る。ローカルの `pathway_project` はネットワークに出ない。

| ツール | 機能 |
|--------|------|
| `pathway_project` | 統計エンジンの行列(`matrix_id`。省略時は `dataset_compare` と同じ)から最小 mzTab-M を書き、MSLipidMapper でサブクラス(と任意の `acyl_chains`)をノードへ集計して 2 群(`group_a` / `group_b`)を `method` と BH で比べる。ノード表と「逆向きに動く隣接ペア」を TSV で返し、JSON-LD・mzTab・manifest は解析フォルダの `reports/pathway/<run_id>/`。ノードの値の和はクラスの総量ではない。`wait_seconds`(既定 20 秒)を超えれば `job_id` を返す(結果は `job_status`)。次は `pathway_plot`。 |
| `pathway_plot` | `pathway_project` の run を PDF に描き、1 ページ目を PNG(長辺 1600px)にする。既定はファイルの場所と要約だけ、`inline_image=true` で画像も返す(1 枚約 2,400 画像トークン)。図中の p 値は MSLipidMapper のもの。 |

## 12. サーバ自身の更新

配布先のクローンを `origin/main` へ追従させる保守ツール(解析の流れには現れない)。遅れていれば `load_dataset` / `dataset_status` に
`update_available` が出る。手順は [DEPLOY.md](DEPLOY.md) の「Receiving Updates」。

| ツール | 機能 |
|--------|------|
| `server_update` | このサーバを `origin/main` へ早送りし、`uv.lock` / `requirements.lock.txt` が変わったときだけ依存も入れ直す。プロセスは触らないので、適用は Claude Desktop(MCP クライアント)を再起動するまで効かない(戻り値の `message`)。未コミット変更・`main` 以外・早送りできない・オフラインでは何もせず `status="refused"` と `reason` を返す。`status` は `up_to_date` / `updated` / `updated_with_warning`(コードだけ新しく、依存の更新に失敗。再起動の前に `setup.ps1`(mac / Linux は `sh setup.sh`)を再実行する)/ `refused`。 |
| `docs_read` | output-format の文書(出力の意味)を返す(リソース `lipidmix://docs/output-format` と `.../{topic}` のツール版)。`topic`(`core` `dataset` `dcl` `eic` `identity` `mztab` `library` `curation` `pathway` `reanalysis` `quant`)・`section`(`## ` 見出し 1 つ)・`outline=True`(見出しと字数だけ)。 |

## 13. 再解析のデータセットのハンドル(並行のエージェント向け)

ハンドルは 1 つのサーバのプロセスにデータセット(とそのライブラリ)を複数並べて持つための ID で、並行のエージェントが
別々のデータセット・極性を互いに干渉せずに扱える。`dataset_handle` を渡した呼び出しは直近を読み書きせず(例外は定量の派生
`dataset_quant_apply` / `dataset_merge_polarities` / `dataset_select_features` で、作った新しいハンドルを直近にする)、
同じハンドルへの呼び出しは 1 本ずつ、別のハンドルどうしは並行に進む。上限は `lipidmix.local.toml` の
`[server] max_dataset_handles`(既定 8)。規則の全体は [docs/output_format/reanalysis.md](docs/output_format/reanalysis.md) の「1. データセットのハンドル」。

| ツール | 機能 |
|--------|------|
| `dataset_list` | 開いているハンドルの一覧(TSV)。`current` は直近のハンドル、`assay_kind` はハンドルごとの種別。 |
| `dataset_close` | ハンドルを閉じる。直近を閉じると、`dataset_handle` を省いた呼び出しは `missing_state` になる(別のハンドルへは落ちない)。長い道具が走っている間は何も閉じずに `HANDLE_BUSY`(`job_status` で待ってから閉じ直す)。戻り値の `source` は開き直すときの手がかり。 |
| `dataset_analysis_run` | 論文の方法を写した spec で行列を作って検定する(数字と `rule_hash` を返し、再現したかの判定は返さない)。`spec.schema` が `analysis-spec.v2` なら統計エンジン(compare / composition / pca / plsda)、`analysis-spec.v1` なら行列の recipe・検定(`student` / `welch` / `mannwhitney` / `ols`)・`aggregate`・`correlate` を spec に書く。spec だけで再現する実行で、`dataset_exclude` の除外とハンドルの `dataset_alignment_qc` の結果は使わない。入力はデータセットか外部の表(TSV・MetaboLights の MAF・mwTab・xlsx)。`output_path` に全量の JSON、省けば TSV(解析ごとに 50 行まで)。spec の形は [docs/output_format/reanalysis.md](docs/output_format/reanalysis.md) の「2. 統計の spec」。 |
| `library_msms_qualify` | MS/MS の照合(`library_match_feature` と同じ)に spec の合否の規則(msms-rules.v1)を当て、`verdict`(`qualified` / `not_qualified` / `not_evaluable`)・検査ごとの結果・`rule_hash` を返す。目標は `inchikey` か `lipid_species`。「明らかに別物か」は判断しない。 |
| `dataset_homolog_check` | 同族列の RT の単調性(同じ二重結合数で炭素数が増えると RT が遅れ、同じ炭素数で二重結合が増えると早まる)を検査する(homolog-check.v1)。特徴を `members` で並べるか、`auto` で脂質クラスの SME 名から集める(プラスマローゲンの列は `lipid_class="PE O-"`)。`verdict`(`pass` / `fail` / `not_evaluable`)・`n_pairs`・違反・重複・`rule_hash`。`wait_seconds` を超えれば `job_id`(結果は `job_status`)。 |
| `dataset_feature_relations` | 特徴 1 つについて、同時溶出の特徴(`coeluting`)と付加体・同位体・源内断片の関係を返す(キュレーションの関係と同じ規則。裏付けは試料間の強度の相関で、MS-DIAL のリンクは使わない)。`as_ion_of`(この特徴が相手のイオン)と `ions_of_this`(相手がこの特徴のイオン)の TSV に `strong` / `ambiguous` / `informational` の印。主な引数は `rt_window`・`mz_tol`・`min_r`(既定 0.7)・`adducts`。 |
| `dataset_claim_evidence` | 1 主張の証拠の束(claim-evidence.v1): 規則で選んだ既定の特徴(`dataset_feature_search` の目標照合と同じ規則)・許容幅や付加体を広げた走査・候補ごとの MS/MS の照合と合否・関係。特徴は決めない。同じ spec・同じ入力なら同じ JSON。`output_path` に束を書いて要約を返す(省けば束を返し、候補は 10 件まで載せる)。欄の意味は [docs/output_format/reanalysis.md](docs/output_format/reanalysis.md) の「7. 証拠の束」。 |
| `dataset_mass_accuracy` | SME の質量誤差(ppm)の分布を m/z 帯(`band_edges`)ごとに返し、勧める許容幅(`recommended_mz_tol_ppm` / `recommended_tolerance`)を添える(当てはしない)。`subset="msms_matched"` で MS/MS が照合された同定に絞る。試料ごとの誤差と注入順の相関は `injection_trend`、理論 m/z の出所は `theoretical_mz_source`。 |
| `dataset_sample_peaks` | 特徴 × 試料ごとのピークの `rt`・`mz`・`height`・`detected` を返す(割れた 2 スポットの判定・RT のドリフトを見る)。戻り値は 1,000 行まで、`output_path` に全量(既存のファイルは置き換える)。試料ごとの RT のあるデータセット(解析フォルダから読んだものか、隣の `.arf` と照合できた mzTab-M)だけ。 |
| `dataset_novel_candidates` | 新しい知見の候補の一覧(novel-candidates.v1): SME・m/z 帯ごとの許容幅・q の閾値で絞り、除く対象を外す。統計は 2 群比較の `result_id`(対比が複数なら `result_id#contrast_id`)か analysis-spec.v1 の解析から。割れたスポットの組は 1 候補に数える。spec の `labels` で外部の印を足し、`ppm_offset` で窓の中心をずらせる。50 行を超えれば先頭 50 行(全量は `output_path`)。 |
| `dataset_alignment_qc` | 1 つの物質が 2 つ以上のアライメントスポットに割れた組を探す(alignment-qc.v1。既定は同じ付加体・\|Δm/z\| ≤ 15 ppm・0.05 < ΔRT ≤ 1.0 分)。型は `injection_block` / `group_aligned`(まとめない)/ `interleaved` で、`mergeable` の組だけを統計の側が試料ごとに大きい方の値でまとめる。全体は `output_path` に JSON。結果はハンドルに最後の 1 つとして持ち、`dataset_feature_search` / `dataset_novel_candidates` が JSON のパスを省いたときに使う(`dataset_analysis_run` と `dataset_claim_evidence` は使わない)。 |
| `dataset_quant_preview` | 定量(1/4)。内部標準の参照表 CSV(`reference_path`。毎回指定)で、データセットのクラスをどの内部標準・係数で割るかを確かめる(値は作らない)。未解決は `needs_input`(`bindings` で選んで呼び直す)。解析フォルダから読んだデータセットだけ。規則プラグイン(`[quant] rules`)が要る([docs/quant-rules-plugin.md](docs/quant-rules-plugin.md))。欄の意味は `lipidmix://docs/output-format/quant`。 |
| `dataset_quant_apply` | 定量(2/4)。preview と同じ引数で内部標準により正規化し、新しいハンドルの派生データセットを作る(元は変えない)。`skip=True` は値を変えず「飛ばした」記録と派生を作る。正規化済み・選別済み・陽陰を統合済みには掛けない(正規化は統合の前に極性ごとに行う)。 |
| `dataset_merge_polarities` | 定量(3/4)。陽・陰イオンのハンドル(`positive_handle` / `negative_handle`。正規化の状態が同じもの)を 1 つに統合し、新しいハンドルを作る。特徴 ID は `positive:<id>` / `negative:<id>`。試料は名前で照合し、食い違えば `needs_input`(`mapping_csv` で対応づける)。試料の表は陽イオンのもの。極性ごとの注入順・バッチは統合の前にそれぞれのハンドルへ `dataset_set_sample_metadata` で入れる。極性ごとのハンドルの除外は引き継ぐが、前の統合したハンドルに当てたシートと除外は統合し直すと引き継がれない。 |
| `dataset_select_features` | 定量(4/4)。キュレーション → 内部標準 → 名前・コメント・MS/MS の除外 → 付加体の採用表(`[quant] adducts`)の順で特徴を選別し、残る行だけの新しいハンドルを作る。採用表が未設定なら `needs_input`(設定するか `skip_adduct_filter=True`)。全判断は `reports/quant/<record_id>.audit.json`。統合の前後どちらにも使える。 |

データセットを読む(または記録の置き場所をデータセットで決める)ツール(`dataset_*`・`curation_*`・`dcl_*`・`library_*`・`pathway_*`・
`record_objective` / `write_report` / `read_report` / `save_figure`)は、省略可能な `dataset_handle` を受ける(省略時は直近のハンドル。
実際の一覧は MCP の `tools/list` の入力スキーマ)。開いていないハンドルは `DATASET_HANDLE_NOT_FOUND` で、直近へは落ちない。

## 根拠バンドル

| ツール | 内容 |
|---|---|
| `dataset_export_evidence_bundle` | 2 群の差次的結果(`result_id`。省けば直近。対比は `contrast_id` か `result_id#contrast_id`)を、全特徴・複数候補付きで根拠バンドル v3(samples・intensities・adduct・comparison.polarity、同じアラインメントの `.dcl` への spectrum_links)に書く。`q_value` はいつも BH の q。探索専用のデータセットも書き、manifest に `exploratory` を足す。旧 v1 は変えない。列の意味は [docs/output_format/reanalysis.md](docs/output_format/reanalysis.md) の「10. 根拠バンドル v3」。 |

出力先は新しいディレクトリを指定する。既存の出力を上書きしない。未同定・手動訂正・前処理の除外の来歴を保持する。
