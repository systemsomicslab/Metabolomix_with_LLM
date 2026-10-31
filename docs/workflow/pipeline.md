# 生データフォルダpipeline呼び出し連鎖

MCP層は薄い(`metabolomix/tools/pipeline_tools.py`)。実体は `metabolomix/pipeline/service.py`
（受付・工程handler一式）と `metabolomix/pipeline/recovery.py`（読取専用status・取消・
再開判定）にある。工程本体は独立したworkerプロセス（`metabolomix/pipeline/worker.py`）が
`metabolomix/pipeline/engine.py::run_engine()` を通じて進めるため、MCP呼び出し自身の
連鎖と、workerが後で進める連鎖は別物として分けて書く。

## pipeline_run

1. metabolomix/tools/pipeline_tools.py  pipeline_run()
2. └─ metabolomix/pipeline/service.py  start_pipeline()
3.    ├─ metabolomix/pipeline/service.py  _prepare_run()
4.    │  ├─ metabolomix/pipeline/request.py  resolve_request()
5.    │  ├─ metabolomix/pipeline/inputs.py  inspect_inputs()
6.    │  ├─ metabolomix/pipeline/service.py  _precheck_manifest()
7.    │  │  └─ metabolomix/analysis/sample_manifest.py  parse_manifest()
8.    │  ├─ metabolomix/pipeline/store.py  find_or_create_run()
9.    │  └─ [manifest不正を検出した場合] metabolomix/pipeline/service.py  _mark_needs_input_without_launch()
10.   ├─ [manifest不正が無かった場合のみ] metabolomix/pipeline/store.py  load_run()
11.   ├─ [manifest不正あり、または既存run(活動中/completed/失敗/取消/部分完了)を再利用した場合] metabolomix/pipeline/service.py  _dispatch_receipt()
12.   ├─ [新規作成直後の`planned` runのときだけ] metabolomix/pipeline/service.py  _launch_and_await()
13.   │  ├─ metabolomix/pipeline/service.py  launch_pipeline_worker()
14.   │  │  └─ metabolomix/core/process_control.py  launch_detached()
15.   │  └─ metabolomix/pipeline/service.py  _await_launch_handshake()
16.   └─ metabolomix/pipeline/service.py  _dispatch_receipt()

事前検査(手順6-7)で壊れた実験情報シート等の既知の不正入力を検出した場合、
`find_or_create_run` でrunは作るがworkerは起動しない(手順9で`needs_input`のまま
保存し、手順10のstatus読み直し自体を行わず手順11の発送receiptへ直行する)。
差次的解析で `comparisons` を指定し忘れた場合はここでは検出しない——それは唯一、
実際にworkerを起動したうえで下流の `resolve_comparisons` 工程が検出する
(brief「下流の群不足だけはworkerを起動できる」)。

**起動しないもう1つの経路（レビュー指摘1）**: manifestが正常でも、
`find_or_create_run`(手順8)が既存run——活動中／completed／失敗・取消・部分完了
済みのいずれか——を再利用した場合、手順10で読み直したstatusが`planned`で
ないので手順11で起動せず終える。新規に作られたrunだけが`status="planned"`の
まま返るので、それだけを起動条件にする。既存runを本当に進めたい場合は
`pipeline_resume`を使う。

手順12の `_launch_and_await` は「launch直前に読んだ`baseline_record`」を引数に
取る中間関数で、`launch_pipeline_worker`(起動)と`_await_launch_handshake`
(handshake待ち)を1つにまとめている——resumeで前回workerの残骸identityが
既に記録されているケースでも、両呼び出しが同じ基準時点の記録を共有できるように
するため(`pipeline_resume`と共有する実装)。

`launch_pipeline_worker` は `sys.executable` で `python -m metabolomix.pipeline.worker`
を起動し、cwdは呼び出し元のcwdや別checkoutではなくこのcheckout自身(`_REPO_ROOT`)へ
固定する。起動プロセスのstdout/stderrは `launch_detached` がログファイルへ結ぶため、
起動元プロセスがworkerのstdoutパイプを継承して待ち続ける経路(`bInheritHandles=TRUE`
かつハンドル未指定)は無い。`_await_launch_handshake` は短い上限だけ待ち、超えても
起動失敗と決め付けず再起動もしない(`launch.handshake="not_confirmed"`のまま
`pipeline_path`を返す。確認は`pipeline_status`に委ねる)。

## pipeline_plan

1. metabolomix/tools/pipeline_tools.py  pipeline_plan()
2. └─ metabolomix/pipeline/service.py  plan_pipeline()
3.    ├─ metabolomix/pipeline/service.py  _prepare_run()
4.    └─ metabolomix/pipeline/service.py  _dispatch_receipt()
5.       └─ metabolomix/pipeline/service.py  _resolved_settings()

`pipeline_run` と事前検査までは完全に同じ経路(`_prepare_run`)を共有するが、
`launch_pipeline_worker` を一切呼ばない——検査・計画・不足情報の保存だけで終える。

手順5の `_resolved_settings` は `plan_pipeline` の発送receiptにだけ足される
(`_dispatch_receipt(include_resolved=True)`)。載せるのは解決済みの
method(`source_path`/`sha256`)・LBM・polarity(`value`/`source`)・`library_mode`・
MSP(ファイル名・sha256・出どころ。置き場所は載せない)・`alignment_reference` だけで、
`raw_stat`・`entries`・`overrides` のような大きな中間データは `record["inputs"]`
に残したまま出さない。サーバー共通指示(`metabolomix/core/mcp_core.py` の ENTRY POINT)が
このツールを「起動前に解決結果を確認する入口」として案内しているため、receiptが
それを持たないと `pipeline_status(include_details=true)` でrecord全体を引くしかなく、
「receiptはコンパクトに」という拘束と衝突する。起動する側(`pipeline_run`)の
receiptは従来どおり最小のまま。

`library_mode="msp_only"` の要求では LBM を解決せず(`Lbm file path` は空)、研究室の参照ライブラリ(.msp)を
`msp_file` か極性の設定(`[library] msp_positive` / `msp_negative`)から決めて実効メソッドの `Msp file path` に書く。
メソッドが宣言する Text DB と注釈器設定表は空にする。MSP の sha256 は入力計画に固定し、計画の指紋にも入れる
(MSP が変われば `pipeline_run` は別の run を作る)。上流をやり直さない `pipeline_resume` では
`metabolomix/pipeline/inputs.py` の `verify_inputs()` が再検査する(変わっていれば `INPUT_CHANGED`、`which: "msp"`)。
receipt にはファイル名・sha256・出どころだけを載せ、置き場所は載せない。

`library_mode="msp_then_lbm"`(v8。脂質の run 用)は研究室 MSP を優先し、LBM で補う。MSP は `msp_only` と同じ規則で決め、
LBM は `[library] lbm`(環境変数 `MSDIAL_LIBRARY_LBM`)**だけ**から決める(`[msdial] lbm`・メソッドの宣言・実行体のフォルダ・
ビルド生成物は見ない。無ければ `LBM_NOT_CONFIGURED`、指す先が無ければ `LBM_NOT_FOUND`)。実効メソッドは `Msp file path` を空にし、
`MSP annotator settings file path` に `inputs/msp-annotators.tsv`(`stage_inputs` が書く。列 `annotator_id` / `msp_file_path` /
`priority`、研究室 MSP を priority 2 で 1 行)、`Lbm file path` に LBM、`Use retention information for LBM-based annotation
scoring` / `filtering` に `False` を書く。メソッドが宣言する MSP 設定表の別名は同じ設定表へ向け(Console は別名を含めて最初の行を
空でも読むため)、Text DB と Text の設定表は空にする。dot product の足切り(weighted / simple / reverse)は MSP と LBM に同じ値を
`Square root of ... dot product cutoff for MSP-based annotation` / `...for LBM-based annotation` で書き、宣言された
`Weighted dot product cutoff for MSP-based annotation` 等は空にする——Console はこの接頭辞の無いキーを**二乗の**足切りとして読む
(0.6 は実効で 0.775)ので、そのままだと MSP だけが LBM(キー無し=既定の平方根の側 0.6 / 0.8)より厳しくなり、同じ証拠の水準に
並ばない。値は平方根の側の宣言 → 接頭辞の無いキーの宣言(平方根の側の値と読む)→ 既定 0.6 / 0.6 / 0.8 の順で、空にした宣言（原本が値を宣言していたもの）は `removed_declarations` に載り、
実効の足切りは receipt の `resolved.dot_product_cutoffs`(`weighted` / `simple` / `reverse` / `scale: "square_root"` / 種類ごとの `source`: `root_key` / `unprefixed_key_read_as_root` / `default`。msp_then_lbm だけ。計画の指紋の外)に出る。Console(45a531c4d)は設定表から 1 件でも読めると `Msp file path` を使わず、
LBM の注釈器の priority を 1 に固定しているので、代表の同定は同じ証拠の水準なら MSP、MSP に MS/MS の照合が無く LBM にあるときだけ
LBM になる。LBM と MSP の sha256 は計画の指紋に入り、`verify_inputs()` が再検査する(`which: "lbm"` の details もファイル名だけ)。
receipt の `resolved.lbm` は `file` / `sha256` / `source`、`resolved.msp` は `priority`(2)も持つ。
設定表の sha256 は配置時に snapshot(`msp.annotator_table_sha256`)へ固定し、再開時の照合が実効メソッドと同じく確かめる
(変わっていれば `STAGED_INPUT_MISMATCH`)。`pipeline_status(include_details=true)` の record の写しでは、LBM の `path` と
`Lbm file path` の上書き、絶対パスの設定表の上書きもファイル名だけにする。mzTab の読み込みは設定表の MSP で InChIKey を補う
(同定は `MspDB_<n>_<MSP の stem>:<名前>`)。

アライメントの参照ファイル(要求の `alignment_reference`、既定 `"median_injection"`): 実効メソッドの
`Alignment reference file ID` を、全試料(`include` を問わない)の注入順の中央の試料(偶数なら小さい方)の
file ID(入力フォルダの主ファイルを `sorted()` した順の 0 始まりの添字)に書く。`_prepare_run` が実験情報シートを
`inspect_inputs` より先に解して行を渡し、`injection_order` が 1 件でも欠ける・シートが無いときはファイル名順の中央
(receipt の `injection_order_source` が `"filename"`)。`"method"` は上書きせず、メソッドの宣言を `method_value` に残す。
`overrides` は計画の指紋の外なので、既存の要求の hash も完了済み run の再利用も変わらない——ただし
古い要求を再送して新しい run を作ると、実効メソッドの参照ファイルが以前と違う(既定の 0 → 中央)。
`alignment_reference` の既定(`"median_injection"`)は要求の hash の外なので、v7 で完了した run のある元フォルダへ既定のまま
v8 の要求を出すと、その run を再利用する(上書きしたメソッドは使われていない)。このとき receipt の
`resolved.alignment_reference` は `null`（v7 の計画に無い）で、参照ファイルが中央に変わっていないことの印になる——中央の参照で
やり直したいときは `alignment_reference` を明示するか、新しい run を作る。

## pipeline_status

1. metabolomix/tools/pipeline_tools.py  pipeline_status()
2. └─ metabolomix/pipeline/recovery.py  read_status()
3.    └─ metabolomix/pipeline/store.py  load_run()
4.    └─ metabolomix/core/process_control.py  same_process()

`load_run` は読取専用で、このツール自身は `pipeline-run.json` を一切書き換えない。
`status="running"`のときだけ`same_process`でworker identityの生存を確認するが、
確認できなくても`status`フィールドは変えず`observed_health`という別軸で返す。

## pipeline_resume

1. metabolomix/tools/pipeline_tools.py  pipeline_resume()
2. └─ metabolomix/pipeline/service.py  resume_pipeline()
3.    ├─ metabolomix/pipeline/recovery.py  prepare_resume()
4.    │  └─ metabolomix/pipeline/request.py  merge_updates()
5.    ├─ [status=="planned"のときだけ] metabolomix/pipeline/store.py  load_run()
6.    └─ [status=="planned"のときだけ] metabolomix/pipeline/service.py  _launch_and_await()
7.       ├─ metabolomix/pipeline/service.py  launch_pipeline_worker()
8.       │  └─ metabolomix/core/process_control.py  launch_detached()
9.       └─ metabolomix/pipeline/service.py  _await_launch_handshake()

`prepare_resume`は新しいattempt/revisionを用意し、どのstageを`pending`へ戻すかを
決めて保存するだけで、それ自体はworkerを起動しない(no-op resumeで無駄な起動を
避ける)。`resume_pipeline`は`prepare_resume`の戻り値`status`が`"planned"`の
ときだけ、改めて`load_run`でlaunch直前の記録(前回workerの残骸identityを含みうる)を
読み直し、それを基準として`pipeline_run`と共有する`_launch_and_await`を呼ぶ
(`launch_pipeline_worker`と`_await_launch_handshake`を1つにまとめた中間関数。
基準を「launch直前」に揃えないと、前回workerの残骸identityを新規launchの
証拠と誤認しうる)。`rerun_upstream=True`を明示しない限りConsole実行はやり直さない。

## pipeline_cancel

1. metabolomix/tools/pipeline_tools.py  pipeline_cancel()
2. └─ metabolomix/pipeline/recovery.py  request_cancel()
3.    └─ metabolomix/pipeline/store.py  load_run()
4.    └─ metabolomix/pipeline/engine.py  cancel_request_path()

保存するのは協調的な取消フラグ(小さなJSON)だけで、ここでは一切のプロセスを
直接終了させない。**rawデータは削除しない。** 上流Consoleの停止は既存の
`execution.supervise`の`cancel_path`監視に、下流の停止はengineのstage境界
チェックに委ねる。受理(このツールが成功したこと)と実際の停止確定
(`cancelled`)は別状態で、後者は`pipeline_status`で確認する。

## workerが1回分を進める経路(別プロセス。上記5ツールとは別に進行する)

`pipeline_run`/`pipeline_resume`が起動したworkerプロセス自身の経路。MCP接続や
グローバルsessionには一切触れない(`metabolomix.core.session_state` /
`metabolomix.core.mcp_core` / `metabolomix.tools.*` を意図的にimportしない設計——
`tests/test_pipeline_engine.py`のASTテストがこれを固定する)。

1. metabolomix/pipeline/worker.py  run_worker()
2. └─ metabolomix/pipeline/engine.py  run_engine()
3.    └─ metabolomix/pipeline/service.py  build_handlers()
4.    └─ metabolomix/pipeline/engine.py  _run_stage_loop()
5.       ├─ metabolomix/pipeline/service.py  _handle_prepare_input()
6.       ├─ metabolomix/pipeline/service.py  _handle_upstream()
7.       │  └─ metabolomix/console/execution.py  supervise()
8.       ├─ metabolomix/pipeline/service.py  _handle_validate_outputs()
9.       ├─ metabolomix/pipeline/service.py  _handle_load_dataset()
10.      │  └─ metabolomix/mztab/loading.py  load_dataset_state()
11.      ├─ metabolomix/pipeline/service.py  _handle_resolve_metadata()
12.      ├─ metabolomix/pipeline/service.py  _handle_preprocess()
13.      │  └─ metabolomix/analysis/dataset_service.py  preprocess_auto()
14.      ├─ metabolomix/pipeline/service.py  _handle_pca()
15.      │  └─ metabolomix/analysis/dataset_service.py  pca_dataset()
16.      ├─ metabolomix/pipeline/service.py  _handle_resolve_comparisons()
17.      ├─ [comparisonごと] metabolomix/pipeline/service.py  _handle_differential()
18.      │  └─ metabolomix/analysis/dataset_service.py  run_comparison()
19.      ├─ [comparisonごと] metabolomix/pipeline/service.py  _handle_export()
20.      │  └─ metabolomix/analysis/dataset_export.py  export_dataset_result()
21.      └─ metabolomix/pipeline/service.py  _handle_report()
22.         └─ metabolomix/pipeline/report.py  write_pipeline_report()
23.    └─ metabolomix/pipeline/engine.py  finish_success()
24.       └─ metabolomix/pipeline/report.py  evaluate_target()

手順4の `_run_stage_loop` は **`build_stages` の計画順**で回す
(`_stage_order`)。`record["stages"]` の挿入順ではない——`prepare_resume` は
resume時に「recordに無いstage_id」を末尾へ追記するので、比較を後から足したrunでは
`differential:*`/`export:*` が `report` の後ろに並ぶ。その順で回すと手順21のreportが
比較を1件も実行していない時点で書かれる。計画に無いstage_id(exploratory目標の
`resolve_comparisons`)だけが計画順のあとに回り、`_mark_out_of_scope` が畳む。

`_run_stage_loop`は**各stageの直前に最新のrun recordを読み直す**。要求版
(`_request_identity`＝revision/content_hash/saved_path)が変わっているか、
このpassで既に通過したstageが外から`pending`へ戻されている
(`_resume_reset_seen`)場合は、その場でstage計画を組み直してpassをやり直し、
`runtime`(DatasetStateと比較結果のキャッシュ)を捨てる。走っている最中に届いた
`pipeline_resume`は——`prepare_resume`が稼働中workerのstatusを動かさないため
`resume_pipeline`が新しいworkerを起こさない——このworkerが取り込むしかなく、
取り込まなければ旧revisionの条件で計算した結果が新revisionの結果として
completedになる。組み直しは有限回(`_REQUEST_RELOAD_MAX_PASSES`)までで、
超えたら`REQUEST_REVISION_CHURN`で止める。

`_handle_upstream`は`metabolomix.console.worker`（単体Console用の監視ワーカー）を
一切呼ばない——`execution.supervise`を直接呼ぶことで、pipeline用worker自身が
「起動して見張る」役を兼ねる(Console用workerの二重起動を避ける)。
`console_job_path`は`analysis-job.json`そのものへのパス(`run_dir = Path(job_path)
.parent`という、コードベース全体の規約)として記録する。run_dirは
**attemptごとに分ける**(`<pipeline_root>/console/attempt-NNNN/`。NNNNは
`upstream` stageの`attempt`)——固定パスにすると`rerun_upstream=true`の再実行が
前回の終了証跡(`execution-result.json`)を上書きし、`msdial.log`を切り詰め
(`process_control`は`_CREATE_ALWAYS`で開く)、`worker.json`の`recovery.befores`
(再収集に要る実行前スナップショット)を置き換えてしまう。規約そのもの
(`console_job_path`＝job.json、`run_dir`＝その親)は変わらない。

`_handle_resolve_metadata`は解決した行(role/group/batch/injection_order/qc_pool)を
`record_updates`で`record["inputs"]["manifest"]`へ、読んだシート自体の内容hashを
`record["inputs"]["manifest_source"]`(`inputs.manifest_source_record`)へ残す。
後者は`recovery.prepare_resume`が「同じパスのシートを書き直した訂正」を
検出するための基準で、これが無いと要求の`sample_manifest`(パス文字列)は
同じままなので、書き直しは変更として現れず旧群割当のまま`completed`が保たれる。ここが唯一の書き手で、
読み手は品質レポートのサンプル来歴セクション(`report._section_sample_provenance`)・
再開時の絶対化(`recovery._absolutize_inputs_paths`、既に絶対ならno-op)の2つ。
明示シートでも自動生成でも同じ形の行を書く。
**`store._relativize_inputs_paths`はこの`manifest`には効かない**——
`create_run`(このstageより前)でしか呼ばれず、以降の`save_run`は`record`を
deepcopyしてそのまま書くだけなので、ここで書く`source_file`は保存後も
絶対パスのまま残る(既知の制限。
`docs/superpowers/notes/2026-09-05-raw-folder-pipeline-validation.md`参照)。

`_handle_pca`はPCA不成立を検出した`PreconditionError`（`metabolomix.analysis.
dataset_analysis.run_dataset_pca`が送出。`DomainError`の派生ではないため
`engine._invoke_handler`の`_NEEDS_INPUT_CODES`whitelistでは検出できない）を
handler自身が捕らえ、`needs_input`のStageResultへ直接変換する。これにより
`_invoke_handler`の汎用`except Exception`分岐（`failed`扱い）を経由させない。

同じ穴を他のhandlerが開け直さないよう、`build_handlers`は全handlerを
`_as_needs_input`で包む——handlerから漏れた`PreconditionError`は
`ANALYSIS_PRECONDITION_MISSING`(kind=missing_state)/
`ANALYSIS_PRECONDITION_INVALID`(kind=bad_request)の`needs_input`になる。
包まないと`error.code`がPythonのクラス名(`PreconditionError`)のまま`failed`に
なり、入力を直せば再開できる停止が「回復不能な失敗」に見える(例:
検出状態を持たないmzTab-Mへ`min_detection_rate>0`を指定した場合)。

`finish_success`は全stageがsucceeded/skippedで止まらずに終えた時点で
`report.evaluate_target(record)`を呼び、その判定(`completed`/`partial`/
`failed`/`needs_input`)をそのまま`record["status"]`へ採用する。差次的解析で
`comparisons`が空のまま(通常は`resolve_comparisons`が先に検出して
到達しない)は防御的な分岐として残る。
