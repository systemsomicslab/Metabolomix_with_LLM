# USAGE — ms-data-parser MCP ツール一覧(全71ツール)

MS-DIAL 出力(`.arf` / `.arf2` / `.pai2` / `.dcl` / `.EIC.aef`)と mzTab-M を解析し、PCA・差次的解析・
アノテーション検証・文献探索・レポート記録までを行う MCP サーバーのツール群です。
おおまかな標準フロー:

```
list_data_files → load_dataset → arf_list_classes / arf_preprocess
   → arf_pca_preprocessed / arf_differential → save_*_figure
   → record_objective → knowledge_coverage → paper_search → ingest_* → write_report

生データから始める場合(別経路。手動で繋ぐなら):
console_plan → console_run → dataset_load → dataset_preprocess
   → dataset_pca / dataset_differential → dataset_export_differential

生データから始める場合(自動で一括統括するなら):
pipeline_run → pipeline_status(確認) → pipeline_resume(訂正・再開が必要な場合)
```

この文書は各ツールの**外形**（引数と用途）を扱います。内部でどのファイルのどの関数を
どの順に呼ぶかは [docs/workflow/](docs/workflow/index.md) を、出力フィールドの意味は
[docs/output_format/](docs/output_format/core.md) を参照してください。

---

## 1. データ探索・ロード(入口)

| ツール | 機能 |
|--------|------|
| `list_data_files` | 指定ディレクトリ内のファイル一覧を拡張子別に取得。**既定は解析できる拡張子のみ**(.arf/.arf2/.pai2/.dcl/.EIC.aef/.mddata/.mdproject)。測定生データ(MS-DIAL が読む abf/ibf/cdf/mzml/wiff/raw/d/wiff2/qgd/lcd/lrp/imzml)も見たいときは `all_files=True`。Agilent・Bruker の `.d` と Waters の `.raw` は**フォルダ**が 1 検体なので末尾に `/` を付けて示す(`.d`/`.raw` 以外は同名のフォルダがあっても計測データとして数えない — MS-DIAL の `isVendorDirectory` と同じ規則)。`extension`(例 `.arf2`)で絞り込みも可。 |
| `load_dataset` | MS-DIAL 出力フォルダの**入口**。arf2 概観→arf 詳細PCAを一括実行。最新バッチ・PeakProperties.arf を自動選択し、以降のツールの既定探索先を更新。 |
| `sample_search` | 因子トークン(`ILG_6h` 等)でサンプルを検索し、`file_id` と `.pai2`/`.dcl` の実パスを返す。specs 省略で語彙一覧。ARF ロード前でも動く。 |

## 2. ARF2 解析(データセット全体カタログ)

| ツール | 機能 |
|--------|------|
| `arf2_parser` | `.arf2`(全体カタログ)を解析しメタデータ概観を要約。サンプル別強度は含まないため多変量解析は不可。 |
| `arf2_annotate_identities` | ARF2 スポット注釈を GOSLIN 正規化・RefMet/LIPID MAPS ID・MSI レベルで一括標準化し TSV 表で返す(オフライン)。**ファイル先頭から `max_rows` 件**で強度順ではない(総数と未表示件数はヘッダ行に出る)。MSI はクラス上限の保守評価。 |

## 3. ARF 解析(サンプル別強度・PCA・差次的解析)

| ツール | 機能 |
|--------|------|
| `arf_parser` | `.arf` を読み込み PCA を実行。スコアプロット用JSON＋Loading上位を返す。タグ/Class ID/群色分け・log変換・検出率足切りに加え、強度閾値(`min_intensity`)・アノテーションキーワード(`annotation_keyword`)でのフィルタも吸収。フィルタ条件を変えた再PCAは引数を変えて本ツールを再呼び出しする(再パースは走らない)。 |
| `arf_list_classes` | 現在の ARF データセットで利用可能な MS-DIAL Class ID 値を一覧。 |
| `arf_list_tags` | ロード済み ARF データセットの MS-DIAL タグを一覧。 |
| `arf_list_sample_roles` | サンプルを sample/qc/blank に分類し TSV 表で返す(前処理適用前の確認)。列= sample/role/group/batch/run_order/excluded。 |
| `arf_exclude` | PCA 外れサンプルや特定ピークを名前/ID で手動除外・再包含(可逆・非破壊)。 |
| `arf_preprocess` | ロード済み ARF 行列に前処理レシピ(正規化・補完・ブランク/QC RSD 足切り・ドリフト補正)を適用し session を更新。 |
| `arf_pca_preprocessed` | `arf_preprocess` 後の前処理済み行列で PCA を実行(生行列経路とは独立)。 |
| `arf_differential` | 前処理後行列で差次的解析(2群 Welch t 検定＋log2FC、BH 補正)。因子トークンによるプール群指定に対応。多群 ANOVA は MCP から非公開(関心の2群を因子指定で切り出す)。 |
| `arf_export_differential` | 直近の差次的結果を InChIKey 付きの 1 ファイルへ書き出す(`output_path`)。同一アラインメントの兄弟 `.arf2` から同定情報を `MasterAlignmentID` で結合する。濃縮解析の背景を保つため、有意行だけでなく **InChIKey が付いた全行**を出す。列定義は `lipidmix/analysis/export_contract.py` が正準(下流リポジトリとの契約)。`apply_curation`(既定 true)でユーザーのキュレーション判断を反映する。 |

## 4. EIC 解析(`.EIC.aef`)

| ツール | 機能 |
|--------|------|
| `eic_parser` | `.EIC.aef` を解析しテキスト要約を返す。 |
| `eic_search_by_mz_range` | m/z 範囲でスポットを検索。既定上限は **1500**(TG/DGDG/CL など m/z 1000 超の脂質を既定で切り落とさないため)。 |
| `eic_search_by_rt_range` | RT 範囲でスポットを検索。 |
| `eic_rank_by_max_intensity` | 各試料のクロマトグラム最大強度の最大値で降順に並べた強度上位ランキングを返す。 |
| `eic_plot_chromatograms` | 指定 `spot_id` のEIC系列を読み、クライアント中立の構造化プロット情報(`lipidmix.eic.v1`)を返す。`file_ids` で最大12試料を選択でき、画像生成・ファイル保存は行わない。 |
| `eic_plot_compounds` | 脂質名/オントロジーで選んだ複数物質のEICを、指定 `file_id` の**1試料分だけ**同一グラフへ重ねる。**既定はPNG画像＋1行キャプション**(`output="image"`)。点列(`lipidmix.eic.multi.v1`)が要るクライアントは `output="payload"` または env `LIPIDMIX_PLOT_OUTPUT=payload`。重ねる本数の既定は **8本**。ARF2の同定をrt/mzで検証し、除外した物質は `selection.dropped` に理由付きで残す。ファイル保存はしない。 |

### EICの描画フロー

1. `eic_search_by_mz_range` または `eic_search_by_rt_range` で描画対象の `spot_id` を確認する。
2. `eic_plot_chromatograms(spot_id, file_ids=...)` を呼び、`series[].x` / `series[].y`、軸、試料、ピーク範囲を含む構造化プロット情報を取得する(座標は4桁に丸め済み)。
3. 描画方法はクライアントに任せる。Use-LLLMではPlotlyのインタラクティブな線グラフとして表示し、Claude Desktop等は各UIの描画方式を使用する。
4. 通常の対話描画ではPNGを生成しない。ユーザーがPNGの生成・保存を明示的に希望した場合だけ `save_eic_figure` を実行する。

複数物質を1枚に重ねる場合は `eic_plot_compounds(file_id, names=[...], ontologies=[...])` を使う。
1呼び出し1試料なので、試料間で比べるときは試料ごとに呼び出して図を並べる。図に現れない物質は
`selection.dropped` の理由（`rt_mismatch` / `file_id_absent` / `below_top_n` など）を確認する。

`normalize` は既定の `"none"` のほか、系列ごとの最大値を1にする `"per_trace_max"` を指定できる。`file_ids` を省略した場合、対象スポットの試料数が12以下のときだけ全系列を返す。

## 5. PAI2 解析(`.pai2`)

| ツール | 機能 |
|--------|------|
| `pai2_parser` | `.pai2`(単一測定ファイル)を解析し、ピーク在庫要約(注釈状況・m/z・RT・強度・S/N の分布と強度上位ピーク)をテキストで返す。単一サンプルのためオミクスPCAは行わない(多変量比較は ARF/ARF2)。 |
| `pai2_inspect_peak` | `peak_id` または `peak_name` で特定ピークの強度・S/N・MS/MS 相当フィールドの有無を返す。 |

## 6. アノテーション検証

| ツール | 機能 |
|--------|------|
| `verify_peak_annotation` | 1ピークのアノテーション妥当性を検証するドシエ(精密質量誤差ppm・アダクト整合・MS/MS 証拠)を返す。要 `pai2_parser` 先行。`analytical_checks.msms.band` は `PASS`(実スペクトル)/`FLAG_ONLY`(取得フラグのみ)/`ABSENT` の3値。 |
| `dcl_parser` | `.dcl`(MSDecResult＝デコンボリューション済み MS/MS)を解析し在庫要約と先頭数件の主要フラグメントを返す。 |
| `dcl_find_msms` | precursor m/z(任意で RT)に一致する MS/MS を引く。アノテーションの裏取り用。該当ゼロは「MS/MS 未取得」であって「フラグメント不在」ではない。 |

## 7. 実験目的の記録(gap駆動探索の前提)

| ツール | 機能 |
|--------|------|
| `record_objective` | 確定した実験目的を `analyses/<analysis_id>.md` に記録。目的・比較・群・小問(Q1..Qn)を保存。`assay_kind`(`lipid`/`metabolite`/`unknown`)を渡すと以後の解釈規則が切り替わる。 |
| `update_objective` | 目的/文脈/状態/`assay_kind` を更新し、創発的な小問を追記。 |

## 8. 知識カバレッジ・文献探索・取り込み

| ツール | 機能 |
|--------|------|
| `knowledge_coverage` | objective の各小問 Qi を COVERED / WEAK / GAP に分類(未探索 GAP が自動探索対象)。 |
| `paper_search` | Europe PMC を検索し撤回・重複除外した候補を返す(ユーザー確認済みクエリ前提)。 |
| `log_search` | 探索結果を objective の探索ログに記録(既探索 Qi の再探索防止)。 |
| `ingest_stage` | 関連候補を `knowledge/_inbox` に speculative 隔離(出典必須)。 |
| `ingest_review_queue` | `_inbox` の保留ノートを analysis_id×Qi でグルーピングして返す。 |
| `ingest_promote` | `_inbox` の保留ノートを `knowledge/` へ昇格(信頼知識化の唯一の経路)。 |
| `ingest_reject` | `_inbox` の保留ノートを破棄。 |

## 9. 図の保存・レポート

| ツール | 機能 |
|--------|------|
| `arf_plot_volcano` | 直近の2群差次的解析を volcano 図として返す。**既定はPNG画像＋件数入りキャプション**(`output="image"`、全特徴を描画)。点列(`lipidmix.volcano.v1`)が要るクライアントは `output="payload"` または env `LIPIDMIX_PLOT_OUTPUT=payload` —— その場合 `up`/`down` は全件、`ns` は `max_points`(既定800)まで等間隔で間引き、件数は `selection` に出る。ファイル保存はしない。 |
| `save_pca_figure` | ユーザーがPNGを明示的に希望した場合だけ、指定した PCA 結果を `reports/figures/` に保存。通常の描画はクライアントUIに任せる。入力元は ARF 経路と mzTab-M 経路の 2 つあり、**どちらも優先しない** — 有効な結果が 2 つ以上あると `AMBIGUOUS_RESULT_SOURCE` で止まるので `source`(`arf`/`mztab`)か `result_id` を指定する(どちらを描くかは図の数字そのものを変えるため)。前処理をやり直して古くなった結果は選べない。探索専用データセットから描いた図には、図の中に但し書きが入る。 |
| `save_volcano_figure` | ユーザーが**ファイルとしての**PNGを希望した場合だけ、指定した差次的解析(2群)を volcano PNG として保存(レポート埋め込み用)。画面で見るだけなら `arf_plot_volcano` が画像を直接返す。描画関数は共通なので同じ図。入力元の決め方(`source`/`result_id`、曖昧なら停止)は `save_pca_figure` と同じ。 |
| `save_eic_figure` | ユーザーが明示的に保存を希望した場合だけ、直近のEICプロット情報を PNG 化し `reports/figures/<analysis_id>_eic.png` に保存。Use-LLLMではローカル書き込みとして承認が必要。 |
| `write_report` | 解析・解釈レポートを `reports/<analysis_id>.md` に上書き保存。 |
| `read_report` | 過去レポートを読み戻す(最新更新のものを返す。セッション継続用)。 |
| `list_reports` | 既存レポートの1行索引(analysis_id / date / status)を返す。 |

## 10. MS-DIAL Console 実行(生データ → mzTab-M)

MS-DIAL 本体を CLI 実行して解析結果そのものを生成する経路。GUI であらかじめ
メソッドファイルを作っておき、以降をこのサーバから回す。成果物は
`analysis-job.json`(受け渡しスキーマ)と mzTab-M で、`dataset_load` が続きを引き取る。

| ツール | 機能 |
|--------|------|
| `console_plan` | 実行計画を作り `analysis-job.json` を生成(`dataset_root`, `method_file`, `polarity`, `measure`, `omics`, `save_project`, `timeout_s`, `lbm_file`)。`dataset_root` は生データフォルダ(リポジトリ外)。**`method_file` は省略できる** — 省略時は `dataset_root` から MS-DIAL GUI が解析のたびに自動保存する `<project>_param_<終了時刻>.txt` を探し、`Ion mode` が `polarity` と一致する最新のものを使う(極性一致が無くて別極性の候補があれば `METHOD_FILE_CHOICE_REQUIRED` に候補一覧を載せて止まり、候補が 0 件なら `METHOD_FILE_NOT_GIVEN`。どちらも `details.searched` に探した場所が入る。探索は `dataset_root` 直下 → 兄弟フォルダ → 過去 run で、候補の列挙だけなら `console_method_candidates`)。ASCII の `key: value` 形式のみで、`.mdproject`/`.mddata` は ZIP のため `METHOD_FILE_NOT_TEXT` で停止する(渡しても全パラメータが既定値のまま走ってしまうため。`.mdproject` の中身は `.mddata` へのポインタだけで解析パラメータを含まない)。メソッドの `Ion mode` と `polarity` が食い違えば `METHOD_FILE_POLARITY_MISMATCH`。**脂質ライブラリ(`.lbm2`)を自動解決する** — メソッドの宣言 → **MsdialWorkbench のビルド生成物** → `MSDIAL_LBM` → `MSDIAL_EXE` と同じフォルダの順。ビルド生成物は `MSDIAL_EXE` から上へ辿って `src/MSDIAL5/MsdialGuiApp/bin/Debug` を持つ階層を探し、exe 自身の TFM と同じコピー → 素の `Debug/` → mtime 最新 の順に 1 本選ぶ(同名の TFM 別コピーは `LBM_AMBIGUOUS` にしない)。`src/MSDIAL4` は NCDK 無しの別世代なので見ない。Console 実行体の exe フォルダに `.lbm2` は無いため、exe フォルダ探索だけでは当たらない。GUI 同様ちょうど 1 件でなければ `LBM_NOT_FOUND`/`LBM_AMBIGUOUS` で停止する(GUI 由来のパラメータは `Lbm file path:` が構造的に必ず空で、空のまま実行すると**警告なしで同定 0 件**になる)。解決したパスは `run_dir/effective-method.txt` に書き、元ファイルは変更しない。`measure` は `peak_height`(既定)のみ正確で、`peak_area_above_zero` は `UNSUPPORTED_AREA_CONSOLE` で停止する。`save_project`(既定 True)は `-p` を付け GUI で開ける `.mdproject` を出させる。`timeout_s` 既定 21600(6 時間)。実行前に `MSDIAL_EXE` が Console 実行体かを `--help` で確かめ、GUI を指していれば `MSDIAL_EXE_NOT_CONSOLE`、未設定なら `MSDIAL_EXE_NOT_FOUND`(`human_action_required` と候補・設定手順を封筒に載せる)。生データフォルダに MS-DIAL 対象の拡張子が 0 種類または 2 種類以上あると `MIXED_RAW_FORMATS`(`required_tools: [console_prepare_input]` を付す)。既存のアライメント結果があれば `warnings` で知らせる。成功すると `session.current_job_path` が設定される。 |
| `console_run` | MS-DIAL Console を実行(`job_path` 省略時は `session.current_job_path`)。要 `console_plan` 先行。同期・切り離しのどちらも**同じ監視経路**(`execution.supervise`)を通り、成功・非ゼロ終了・タイムアウト・取消・起動不能の**どの終了経路でも終了証跡 `run_dir/execution-result.json` と生成物を残す**。生成物は **`-o`(エクスポート)と生データフォルダ(`.pai2`/`.dcl`/`.arf`/`.arf2`/`.EIC.aef`)の 2 か所**に分かれて出るので両方を収集し、どちらから来たかを `root` に記録する。**終了コード 0 だけでは `completed` にしない** — 主 mzTab-M が一意に選べ、構造が妥当で、定量行列に有限値があり、予定した入力が全て assay に 1 対 1 で対応することまで確認する(満たさなければ生成物があれば `partial`、無ければ `failed`)。`detach=True` にすると**監視ワーカーを親から切り離して起動し、待たずに戻る**(実データ 60 サンプルで約 44 分。同期実行では 1 ツール呼び出しがその間戻らず、呼び出し元が中断されると生成物ごと失う)。返る `pid` は**監視ワーカーのもので Console のものではない**(`pid_of: "worker"`。Console の pid は終了証跡に載る)。ワーカーが最後まで監視して生成物の収集とジョブの確定まで行うため、`console_status` を呼ばなくても結果は確定する。同じジョブを別プロセスが実行中なら `JOB_BUSY` で拒否する(job 単位の OS ロック)。 |
| `console_status` | ジョブの現在のステータスを返す(`job_path`, `include_artifacts`)。生成物は**既定では全文を返さない** — `artifact_count` と 役割別内訳 `artifacts_by_role` だけを返す(1 サンプルにつき 5 件出るため 60 サンプルの実走で 313 件になり、全文 TSV は後ろの `warnings`/`error` を埋没させた)。個々のパスが必要なときだけ `include_artifacts=True` で `path`/`role`/`format`/`root` の TSV(列名 1 回)を得る。`execution` に `save_project`/`timeout_s` が入る。 `server_version` に稼働中のサーバ版数(`<version>+<git SHA>`)が入る——更新後に MCP サーバを再起動し忘れると古いプロセスが黙って動き続けるため。**読み取り専用**で、完了処理はしない(監視ワーカーが確定させる)。終了証跡があれば `execution_receipt` に `termination`/`exit_code`/`execution_id` が入る。旧方式(`.detached-state.json` が残っているジョブ)は終了コードも実行の同一性も記録されていないため、プロセスが消えていても**生成物の有無だけでは完了と判定せず** `EXECUTION_UNRESOLVED` を返す(まだ生きていれば `detached.alive`)。 |
| `console_prepare_input` | 計測フォーマットが 1 種類だけの入力フォルダを作る(`dataset_root`, `keep_extension`(既定 `wiff`), `out_dir`)。`MIXED_RAW_FORMATS` の解消手段。指定拡張子の計測ファイルと**その随伴ファイル**(`.wiff.scan` / `.timeseries.data` 等)だけを**ハードリンク**で集める(同一ボリュームでなければコピー)。Agilent・Bruker の `.d` と Waters の `.raw` のように**フォルダが 1 検体**の形式も扱え、その場合はフォルダ構造を実体として作り直し中のファイルだけをリンクする(フォルダ自体はリンクにしない)。既にある出力先フォルダは中まで見て欠けているファイルだけを補う。元フォルダは一切変更しない。`out_dir` 省略時は `<元フォルダ>_<拡張子>` を兄弟として作る。 |
| `console_method_template` | 既存パラメータから別極性用のメソッドファイルを作る(`out_path`, `polarity`, `based_on`, `dataset_root`, `omics`)。**まず `console_plan` の `method_file` 省略を試すこと** — その極性で一度でも GUI 実行があれば作る必要は無い。`based_on` を省略すると `dataset_root` 直下・兄弟フォルダ・過去 run から候補を探し、1 件なら採用、複数なら `METHOD_FILE_CHOICE_REQUIRED` で選ばせる(土台の選択は解析条件そのものなので最新を黙って採らない)。差し替えるのは `Ion mode` と `Searched adduct ions` だけで、検出・アライメント条件は元のまま引き継ぐ。空の `Lbm file path` も解決して埋める(解決順は `console_plan` と同じで、ビルド生成物が `MSDIAL_LBM` より優先)。 |
| `console_method_candidates` | 使えるメソッドファイルの候補を列挙する(`dataset_root`, `polarity`, `omics`, `search_dirs`)。**`console_plan` が失敗する前に呼べる** — 手動実行 UI の「メソッドを選ぶ」画面はこれを叩く。探すのは `dataset_root` 直下 → 兄弟フォルダ直下 → `runs/*/analysis-job.json` の `software.method_file`(過去 run が実際に使ったメソッド)で、**再帰はしない**。各候補に `origin`(`same_dir`/`sibling`/`past_run`/`given`)、`usable`(`direct` / `needs_polarity_conversion`)、比較用の `key_params` が付く(候補 10 件超では `key_params` を付けない)。別極性の候補は自動採用せず `console_method_template` を経由させる。 |
| `console_cleanup` | あるジョブが生成したファイルだけを一覧・削除する(`job_path`, `dry_run`(既定 True))。再実行のたびに生データフォルダへ別タイムスタンプの一式が積まれる問題への対処。**消すのは `analysis-job.json` が記録した生成物だけ**で、タイムスタンプの推測はしない(記録が無いジョブは `NO_JOB_OUTPUT` で拒否)。生データには触れない。削除後のジョブ status は `cleaned`。**監視ワーカーが実行中のジョブは削除を拒否する**(`JOB_BUSY`。書き込み中のファイルを実行中のプロセスから奪わないため。所有者の生死は pid ではなく process identity で見る)。`dry_run` の一覧はその場合も返すが `warnings` を添える。 |
| `job_list` | `dataset_root/runs/` 以下のジョブ一覧を新しい順に返す。 |

## 11. DatasetState 解析(mzTab-M 経路)

mzTab-M 2.0 を正準状態(`DatasetState`)として読み、ARF 経路と同じ前処理・PCA・
差次的解析を回す。**`session.arf` / `.arf2` / `.pai2` / `.eic` とは独立したスロット**なので、
ARF 経路の状態を壊さない。

| ツール | 機能 |
|--------|------|
| `dataset_load` | mzTab-M 2.0 を読み `session.dataset` を作る。隣接する `.arf` が同一アライメントだと**数値で検証できた場合だけ**、検出状態(gap-fill の区別)を evidence sidecar として取り込む(要約は戻り値の「検出状態」行)。`mztab_path`(絶対パス直指定)・`job_path`(ジョブの宣言 polarity・measure で正準エントリを一意選択。曖昧なら停止)・`pipeline_path`(`pipeline_run` の pipeline-run.json かその親。その run の**解決済み試料対応表・feature binding・解析行列もまとめて**セッションへ載せる。行列だけを別に載せる入口は無い——別 run の行列を混ぜると群の対応が黙ってずれるため)の**いずれか一つ**を渡す(全省略・2つ以上の指定はエラー)。`console_run` 後は `job_path` 推奨(polarity・measure が確定済み)。生成物の絶対パスは記録された `root` から復元する。ジョブが `completed` 以外(`partial` 等)なら**既定で読み込みを拒否する**(`INCOMPLETE_ANALYSIS_JOB`)。探索目的なら `allow_incomplete=True` を明示する — その場合は中断された実行の生成物である旨を警告の先頭に差し、**探索専用**として `dataset_differential` と `dataset_export_differential` を `EXPLORATORY_ONLY_DATASET` で拒否する(欠けた検体が「その群には無い」ように見え、比較がその欠落ごと結論になるため)。出所の裏取り `source_verification` を返す: `verified`(終了証跡と mzTab の hash が一致した完了実行)/`legacy_unverified`(ジョブ経由だが証跡が無い・一致しない)/`direct_unverified`(`mztab_path` 直接指定)。**`completed` という文字列だけでは `verified` にしない**。要約には**同定の出所内訳**の行が出る(MS/MS 証拠あり = SME 由来 / MS1 注釈のみ = SML 由来 / 同定なし)。同じ内訳は `inchikey_coverage["identified_by"]`(`sme` / `sml_only` / `none`)にも入る。**MS1 注釈は m/z 照合だけで付くので MS/MS 裏付けと同一視しないこと** — MS-DIAL は Text DB 由来の同定を SME へ書かないため、Text DB 運用では SME 0 件・SML のみが正常な形になる(`lipidmix://docs/output-format/identity` §13)。 |
| `dataset_status` | 現在の `DatasetState` の概要を返す。`samples` に name/role の TSV が入り、`dataset_differential` の `group_a`/`group_b` はここに出る名前をそのまま使う。`detection` は検出状態の有無(`available`)と gap-fill 率。**`available=false` の間は検出率・欠測率を語ってはいけない**(非ゼロ値に補間値が混在する)。 |
| `dataset_preprocess` | 定量行列に前処理レシピを適用(引数は `arf_preprocess` と同一: `normalize` / `blank_min_fold` / `drift_correct` / `max_qc_rsd` / `impute` / `min_detection_rate`)。`min_detection_rate` は gap-fill を除いた実検出率での足切りで、`dataset_load` が検出状態を取り込めた場合にだけ使える(無い状態で 0 より大きい値を渡すと引数エラー。黙って未検出 0 件として通さない)。フィルタは正規化・補完より**前**に掛ける。注入順とバッチは mzTab-M の `MTD assay[N]-custom[...]`(`MS:4000089` injection sequence label / `MS:4000088` batch label)から読むので **`drift_correct` は実際に適用される**。注入順が無いファイルでは従来どおり未実施の caveat が出る。バッチラベルは 2 値以上あるときだけ採用し(MS-DIAL の既定は全件 `1` で情報を持たないため)、無ければファイル名の日付推定に戻す。採用元は `sample_meta` の `batch_source` / `run_order_source` に入る。 |
| `dataset_pca` | 前処理済み `DatasetState` で PCA(`n_components` 既定 5、`log_transform` 既定 False)。ローディング全量は戻り値に載せず `session.dataset.last_pca` に保持する。 |
| `dataset_differential` | 前処理済み行列で 2 群比較(Welch t 検定＋BH-FDR)。`group_a`/`group_b` は**サンプル名のリスト**(`dataset_status` の `samples` で確認)。**log2FC は正なら `group_b` が高い**(`group_a` が基準)。全特徴量の結果と volcano 点列は `session.dataset.last_differential` に保持する。 |
| `dataset_export_differential` | 指定した差次的結果(`result_id` 省略時は直近)を InChIKey 付きの 1 ファイルへ書き出す。**前処理をやり直した後の古い結果は書き出さない**(古い数字に現在の前処理条件のラベルが付いた TSV は、どちらも正しく見えてずれが分からない)。メタ行に `result_id` / `preprocess_id` / `source_verification` が入り、前処理条件は結果自身の来歴から書く。出力は一時ファイルから置換して確定する。**`arf_export_differential` と同一の契約**(15 列 + `contract_version` メタ行)なので下流のパスウェイ解析にそのまま渡せる。InChIKey は mzTab-M 由来(`.arf2` との結合は不要)。`ontology` と `msi_level` は mzTab-M に対応物が無く空欄で、その旨をメタ行に書く。`apply_curation`(既定 true)でユーザーのキュレーション判断を反映する。 |
| `dataset_build_matrix` | v2 の解析行列(`analysis-matrix.v1`)を、いま session にある DatasetState から1本作る。`recipe` は profile の `matrix_recipes` 1要素と同じ5キー(`base`/`normalize`/`drift_correct`/`filter`/`impute`)で、検証規則は profile と共有する。`matrix_id` は recipe と dataset・metadata から決まるので、同じ入力なら同じ ID になる。揃っていない前提は埋めない — `base="internal_standard_ratio"` は対応付けが要るため `missing_state` で止まり、検出状態が不明なら `filter.min_detection_rate` は `not_evaluable` として feature を落とさない。値は戻り値に載せず session に保持し、`dataset_statistic` へは `matrix_id` を渡す。 |
| `dataset_statistic` | v2(LC–MSメタボロミクス)の統計を、**名指しした解析行列**(`matrix_result_id`、`analysis-matrix.v1`)に対して1件実行する。`specification` は profile の `analysis_recipe.statistics` 1要素と同じ形(`kind` は `welch`/`anova_tukey`/`pca`)。直近の前処理を暗黙に使わない — recipe 違いの行列が2本ある前提の設計で、どの行列の数字かを結果に残す。`dataset_differential`(v1)とは効果量の定義が違う: v2 は**統計変換前の算術平均比の log2**(`effect_size_definition=log2_arithmetic_mean_ratio`)で、log2 変換に pseudocount を足さず、BH の母集団は検定できた feature だけ。行列は pipeline の preprocess/qc_processed か `dataset_build_matrix` が作る。 |
| `dataset_set_sample_metadata` | 実験情報シート(`sample-manifest.v1`、TSV)を読み込み、全件検証してから `session.dataset` へ一括反映する。上流の Console 実行をやり直さずに群・バッチ・注入順・qc_pool 等の誤記入を訂正できる。**一部だけ適用される状態は作らない**(シートに不備があれば `session.dataset` を一切変更せずエラーを返すので、そのまま再送してよい)。前処理入力(role/batch/injection_order/qc_pool/include/sample_id/source_file)が変われば前処理済み行列・PCA・差次的解析まで無効化するが、**`group` だけの変更なら差次的解析だけを無効化し PCA は生かす**(`changed_fields` で確認できる)。**`dataset_differential` を直接呼ぶ経路との違い**: `dataset_differential` は呼び出し側が `group_a`/`group_b` にサンプル名のリストを都度自分で組み立てる汎用ツールで、群名から対照/処置の向きを決めたり QC/blank/unknown・除外行を弾いたりはしない。pipeline 側が使う明示比較(`comparison_id`/`reference_group`/`test_group`)は、群名だけでは向きを推測しない・QC/blank/unknown・`include=false` を混ぜない・完全交絡(群⟂バッチが分離不能)を `allow_confounded=true` の明示なしには実行しない、という前提検証を経てから同じ計算を呼ぶ(内部関数 `resolve_comparison`/`run_comparison`。バッチ情報が無い場合は「評価不可」であって「交絡あり」ではない)。 |

## 12. Class ID に無い因子で絞る・比べる

MS-DIAL の Class ID は入力された1文字列にすぎず、時点や複製がサンプル名にしか
無いことがある(例: Class ID = `ILG`/`control` だけで、時点 `6h` はサンプル名のみ)。
`class_ids` / `group_a` / `group_b` / `group_levels` / `group_factors` のトークンは
Class ID とサンプル名の**両方**から解決されるため、そのまま書ける。

```python
sample_search()                                    # 何で絞れるかの語彙一覧
arf_parser(class_ids=["ILG_6h", "control_6h"])     # 6h だけで PCA
arf_parser(group_factors=[["control", "LPS", "ILG", "G_uralensis"],
                          ["0h", "15min", "1h", "6h", "24h"]])  # 処置x時点の20群で色分け
arf_preprocess(normalize="median")
arf_differential(group_a="ILG_6h", group_b="control_6h")        # 時点を揃えた2群比較
```

QC/blank の扱いはツールごとに異なる。`arf_parser` の `class_ids` フィルタでは既定で除外され
(`include_roles=["sample","qc"]` で戻せる)。`sample_search` は既定で全 role を返す
(role を絞りたいときだけ `include_roles` を渡す)。`arf_differential` は QC/blank を常に
比較対象から除外し、オーバーライドはできない。PCA の色分けは `group_levels`/`group_factors`
で軸を指定した場合に `"qc"`/`"blank"` ラベルとして残る(前処理品質の判断材料になるため。
軸を指定しない場合は生の Class ID でラベルされる)。

## 13. 生データフォルダからの統括(pipeline、生データ → 上流検証 → 前処理 → PCA → 差次的解析 → レポート)

生データフォルダ1つを渡すだけで、上流(MS-DIAL Console)の実行・検証、`dataset_load`
相当の読込、メタデータ解決、前処理、PCA、指定した2群比較の差次的解析、TSV/volcano
書き出し、品質レポート作成までを1本の永続 run(`pipeline-run.json`)として自動で進める
経路。個々のツール(`console_plan`/`console_run`/`dataset_preprocess` 等)を手動で繋ぐ
代わりにこちらを使う。中断・再送・訂正をまたいでも同じ入力からは同じ結果になる
(fingerprint と request_id による冪等性)。詳しい内部の呼び出し順は
[docs/workflow/pipeline.md](docs/workflow/pipeline.md) を参照。

| ツール | 機能 |
|--------|------|
| `pipeline_run` | 生データフォルダ(`dataset_root`)を渡すだけの通常入口。入力検査・計画保存・workerの起動までを一括で行い、短時間で `pipeline_path` を返す(工程自体は非同期に進む)。`request`(省略可)で `target`/`polarity`/`method_file`/`sample_manifest`/`preprocess`/`comparisons` 等を指定できる。省略項目は既定値で埋まる(measure=peak_height、preprocess=conservative-v1、comparisons=空 等)。`dataset_root` 直下の `analysis-request.json` / `sample-manifest.tsv` は既定名として読み、優先順位は「`request` の明示値 > `analysis-request.json` > 既定値」(`preprocess` は子キー単位で重なる)。`request_id` は冪等性キーで、同一内容の再送は同じ結果を返し別内容は `IDEMPOTENCY_CONFLICT`。壊れた実験情報シート等の既知の不正入力は、Console を起動せず `needs_input` を保存する。差次的解析で `comparisons` を指定し忘れた場合だけは、workerを起動したうえで `resolve_comparisons` 工程が `COMPARISON_REQUIRED` の `needs_input` として検出する。 |
| `pipeline_plan` | `pipeline_run` と同じ入力検査・計画保存だけを行い、workerは起動しない。実行前に条件(不足情報と、receiptの `resolved` が返す method/lbm/polarity)を確認したいときに使う。 |
| `pipeline_status` | 実行中の run の状態を読む(`pipeline_path`, `include_details=False`)。**読み取り専用** — 監視や成果物の確定はこの呼び出しに依存しない(呼ばなくても裏の worker が最後まで進める)。`status`/`stage_statuses`/`needs_input`/`warnings` の要約を返し、`include_details=True` で永続レコード全体を追加取得できる。`status="running"` のとき、workerの生存確認 `observed_health`(`ok`/`worker_missing`/`unknown`)を添える。`worker_missing` なら `recovery_hint` で `pipeline_resume` を促す。 |
| `pipeline_resume` | 入力訂正(`updates`: `target`/`sample_manifest`/`preprocess`/`comparisons` に限定)、下流の再計算、または中断した工程からの再開を行う。`rerun_upstream=True` を明示しない限り Console 実行はやり直さない(自動再試行はしない)。`request_id` を指定した再送は、直前と同じ内容なら新しい revision を作らず直前の結果を返す。 |
| `pipeline_cancel` | 取消要求を保存する(`pipeline_path`)。**生データを削除しない** — 保存するのは協調的な取消フラグだけで、実際に停止したかどうかは `pipeline_status` で確認する。二重に呼んでも同じ状態に落ち着く(冪等)。 |

## 14. 参照ライブラリ(MS/MS スペクトル照合)

`.dcl` の測定 MS/MS を参照ライブラリ(`*_Loaded.msp2.dbs` または `*.msp`)と照合し、
対向プロット(mirror plot)で確認する経路。`library_load` → `library_match_feature` →
`library_plot_mirror` の順で呼ぶ(前段の状態が無いと `missing_state` 封筒を返す)。

| ツール | 機能 |
|--------|------|
| `library_load` | 参照ライブラリを解決し、照合用の SQLite store を構築(または既存キャッシュを再利用)して `session.library` へ持つ。要約(`record_count`/`ion_modes`/`compound_classes` 上位10/`search_params`/`records_without_ion_mode`/`non_utf8_lines`)を返す。`.msp` 由来で `search_params` が無いときは既定の許容幅を使う旨を明示する。`ion_mode`(`"positive"`/`"negative"`)を渡すと環境変数 `MSDIAL_MSP_POS`/`MSDIAL_MSP_NEG` が指す研究室ライブラリを読む。解決順は `file_path` → `ion_mode` の環境変数 → `*_Loaded.msp2.dbs` → 設定済みの環境変数 → `*.msp` で、候補が 1 つに決まらなければ `code="MSP_AMBIGUOUS"` のエラーを返す(更新日時で黙って選ばない)。一度に持つライブラリは 1 つ。大きな `.msp` の初回構築は `python -m lipidmix.library.store --ion-mode <極性>` で事前に済ませられる([docs/cli.md](docs/cli.md))。 |
| `library_match_feature` | 測定 MS/MS(`.dcl`、間引かずに全ピーク使用)を候補と照合し、上位をスコア付きの TSV で返す。**並び順は `total_score`**(MS-DIAL の `GetTotalScore` の移植。RT・precursor m/z の一致度を足した**正規化しない和**で 1 を超える)。RT 項を足したかは `scoring.use_rt` に出る(`.dbs` の `IsUseTimeForAnnotationScoring` 次第。`.msp` では常に false)。上流の順位付けタプル全体は再現していないので**上位に化学的にありえない候補が残ることがある** — 1 位の妥当性は対向プロットで確かめること。`.dcl` に MS/MS が無いときは `not_found`(「未取得」であって「合わなかった」ではない)。スペクトル座標は戻り値に含めず `session.library.last_match` へ持つ。 |
| `library_plot_mirror` | 直近の照合結果から対向プロット(上段=測定・下段=参照)を描く。既定は PNG 画像、`output="payload"` で座標 JSON。`scale` で縦軸の写し方を選ぶ(`"relative"` 既定 / `"sqrt"` / `"log10"`)——**precursor がベースピークのスペクトルは `relative` だと診断イオンが潰れて読めない**ので `"sqrt"` を使う。`label_policy` で m/z ラベルの衝突回避を選ぶ(`"auto"` 既定 = 水平・垂直の 2 次元判定 / `"msdial"` = 上流に忠実な水平のみの判定。忠実版のほうがラベルは少ない)。`.dbs` の強度足切りで**採点に入らなかった測定ピーク**は灰色で薄く描き分け、凡例と caption に件数が出る(足切りが 0 の run では現れない)。 |

## 15. アラインメントのキュレーション(注釈の一覧確認と機械判別)

注釈付きスポットを、EIC・対向プロット・Δppm・ΔRT・脂質クラス×不飽和度の RT–m/z 傾向と
機械判別つきで**一覧**確認し、ユーザーが付けたフラグ(「間違い」「疑わしい」)を記録する経路。
`library_load`(アラインメントと同じフォルダの `*_Loaded.msp2.dbs` を推奨) → `curation_review` →
ユーザーがビューアでフラグを付けて「送信用テキストをコピー」→ チャットに貼る → `curation_submit`。
フラグの無いスポットは「間違っていない」の意で、何も記録しない。

| ツール | 機能 |
|--------|------|
| `curation_review` | 対象(既定は注釈付き全部。`ontology` でクラス、`name_contains` で名前の部分一致)を選び、証拠収集・機械判別・HTML ビューア生成を 1 回で行う。戻り値は suspect 以上とフラグ済みのスポットだけの TSV・判定の件数・クラス別の傾向要約・`html_path`。EIC 系列とスペクトルは返さない。判定は `likely_wrong`(polarity_mismatch / precursor_unmatched)/ `suspect`(ppm_out / low_score / drt_out / eic_poor、または弱い兆候の重なり)/ `ok`。RT–m/z 傾向は補強にしか使わない。`thresholds` で既定値を上書きできる。1 回の上限は 3000 スポット。 |
| `curation_submit` | フラグを追記する。ビューアの送信用テキストを `submission_text` にそのまま渡すか、`review_id` + `flags=[{spot_id, flag: wrong\|suspect\|clear, note}]`。不正な要素が 1 つでもあれば何も書かない。`source` は `user`(ユーザー自身の判断)/ `llm`(LLM の提案にユーザーが同意したもの)。 |
| `curation_flags` | 現在のアラインメントで有効なフラグを TSV で返す。 |
| `curation_view_data` | ビューア(MCP Apps)専用。LLM は呼ばない。 |

フラグは `<アラインメントのフォルダ>\curation\flags.jsonl` に、`.arf2` の sha256 と
MasterAlignmentID の組で記録する(MS-DIAL を再実行すると古いフラグは当たらない)。
`arf_export_differential` / `dataset_export_differential` は既定(`apply_curation=True`)で `wrong` のスポットを同定なしとして出力から外し、メタ行 `# curation = ...` で宣言する(15 列の契約は不変)。
MCP Apps で会話内にビューアを出す経路は、Claude Desktop のローカル
stdio サーバでは 2026-09 時点で未検証(ブラウザで `html_path` を開く経路が主)。

## 16. サーバ自身の更新

配布先のクローンを `origin/main` へ追従させるための保守ツール。解析の流れには現れない。
遅れているときは更新通知が出る(`load_dataset` の冒頭行、`dataset_status` /
`console_status` の `update_available`)ので、それを見てから呼ぶ。手順の全体は
[DEPLOY.md](DEPLOY.md) の「Receiving Updates」を参照。

| ツール | 機能 |
|--------|------|
| `server_update` | このMCPサーバ自身を `origin/main` へ早送りし、`requirements.txt` が変わっていたときだけ依存も入れ直す。**プロセスは操作しない** — 適用後に Claude Desktop(MCPクライアント)を再起動するのは利用者。起動中のプロセスは import 済みの古いコードで動き続けるため、戻り値の `message` にある再起動の指示をそのまま利用者へ伝えること。未コミット変更がある・`main` 以外にいる・早送りできない(独自コミットがある)・オフラインのいずれでも**何もせずに** `status="refused"` と機械可読な `reason` を返す(例外は投げない)。`status` は `up_to_date` / `updated` / `updated_with_warning` / `refused`。`updated_with_warning` はコードだけ新しく依存の更新に失敗した状態で、再起動前に手元で `pip install -r requirements.txt` が要る。 |
