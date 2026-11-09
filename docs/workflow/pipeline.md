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
5.    │  ├─ [単位（極性）ごと] metabolomix/pipeline/inputs.py  inspect_inputs()
6.    │  ├─ [単位（極性）ごと] metabolomix/pipeline/service.py  _precheck_manifest()
7.    │  │  └─ metabolomix/analysis/sample_manifest.py  parse_manifest()
8.    │  ├─ metabolomix/pipeline/store.py  find_or_create_run()
9.    │  └─ [manifest不正を検出した場合] metabolomix/pipeline/service.py  _mark_needs_input_without_launch()
10.   ├─ [manifest不正が無かった場合のみ] metabolomix/pipeline/store.py  load_run()
11.   ├─ [manifest不正あり、または既存run(活動中/completed/失敗/取消/部分完了)を再利用した場合] metabolomix/pipeline/service.py  _dispatch_summary()
12.   ├─ [新規作成直後の`planned` runのときだけ] metabolomix/pipeline/service.py  _launch_and_await()
13.   │  ├─ metabolomix/pipeline/service.py  launch_pipeline_worker()
14.   │  │  └─ metabolomix/core/process_control.py  launch_detached()
15.   │  └─ metabolomix/pipeline/service.py  _await_launch_handshake()
16.   └─ metabolomix/pipeline/service.py  _dispatch_summary()

事前検査(手順6-7)で壊れた実験情報シート等の既知の不正入力を検出した場合、
`find_or_create_run` でrunは作るがworkerは起動しない(手順9で`needs_input`のまま
保存し、手順10のstatus読み直し自体を行わず手順11の発送receiptへ直行する)。
要求は `pipeline-request.v3`(`request.resolve_request`)で、1 極性の形も 2 極性の
`units` の形も単位の並びに正規化する。入力計画(`inspect_inputs`)と事前検査は単位ごとに
回し、run(`pipeline-run.v3`)は単位ごとの記録 `record["units"][<極性>]` を持つ。
下流のキー(`target` / `preprocess` / `comparisons`)は `PIPELINE_DOWNSTREAM_REMOVED`
(元フォルダの `analysis-request.json` に書かれたものは読み飛ばして警告
`REQUEST_FILE_DOWNSTREAM_IGNORED`)。明示の `units` は単位を丸ごと決めるので、
`analysis-request.json` のトップの単位のキー(`polarity` `method_file` …)は使わず、警告
`REQUEST_FILE_UNIT_KEYS_IGNORED` で伝える。2 単位の `raw_dir` は同じフォルダも、片方が
もう片方の中にある(`.` を含む)のも `PIPELINE_REQUEST_INVALID`。

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
4.    └─ metabolomix/pipeline/service.py  _dispatch_summary()
5.       └─ metabolomix/pipeline/service.py  _resolved_settings()

`pipeline_run` と事前検査までは完全に同じ経路(`_prepare_run`)を共有するが、
`launch_pipeline_worker` を一切呼ばない——検査・計画・不足情報の保存だけで終える。

手順5の `_resolved_settings` は `plan_pipeline` の発送receiptにだけ足される
(`_dispatch_summary(include_resolved=True)`)。載せるのは解決済みの
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
メソッドが宣言する Text DB と注釈器設定表は空にする。dot product の足切りは下の `msp_then_lbm` と同じく平方根の側のキーで書き(2026-11-02 から。それより前の msp_only の run は接頭辞の無いキーを二乗で読ませていた)、receipt の `resolved.dot_product_cutoffs` に出す。MSP の sha256 は入力計画に固定し、計画の指紋にも入れる
(MSP が変われば `pipeline_run` は別の run を作る)。上流をやり直さない `pipeline_resume` では
`metabolomix/pipeline/inputs.py` の `verify_inputs()` が再検査する(変わっていれば `INPUT_CHANGED`、`which: "msp"`)。
receipt にはファイル名・sha256・出どころだけを載せ、置き場所は載せない。

`library_mode="msp_then_lbm"`(v8。脂質の run 用)は研究室 MSP を優先し、LBM で補う。MSP は `msp_only` と同じ規則で決め、
LBM は `[library] lbm` **だけ**から決める(メソッドの宣言・実行体のフォルダ・
ビルド生成物は見ない。無ければ `LBM_NOT_CONFIGURED`、指す先が無ければ `LBM_NOT_FOUND`)。実効メソッドは `Msp file path` を空にし、
`MSP annotator settings file path` に `inputs/msp-annotators.tsv`(`stage_inputs` が書く。列 `annotator_id` / `msp_file_path` /
`priority`、研究室 MSP を priority 2 で 1 行)、`Lbm file path` に LBM、`Use retention information for LBM-based annotation
scoring` / `filtering` に `False` を書く。メソッドが宣言する MSP 設定表の別名は同じ設定表へ向け(Console は別名を含めて最初の行を
空でも読むため)、Text DB と Text の設定表は空にする。dot product の足切り(weighted / simple / reverse)は MSP と LBM に同じ値を
`Square root of ... dot product cutoff for MSP-based annotation` / `...for LBM-based annotation` で書き、宣言された
`Weighted dot product cutoff for MSP-based annotation` 等は空にする——Console はこの接頭辞の無いキーを**二乗の**足切りとして読む
(0.6 は実効で 0.775)ので、そのままだと MSP だけが LBM(キー無し=既定の平方根の側 0.6 / 0.8)より厳しくなり、同じ証拠の水準に
並ばない。値は平方根の側の宣言 → 接頭辞の無いキーの宣言(平方根の側の値と読む)→ 既定 0.6 / 0.6 / 0.8 の順で、空にした宣言（原本が値を宣言していたもの）は `removed_declarations` に載り、
実効の足切りは receipt の `resolved.dot_product_cutoffs`(`weighted` / `simple` / `reverse` / `scale: "square_root"` / 種類ごとの `source`: `root_key` / `unprefixed_key_read_as_root` / `default`。msp_only と msp_then_lbm。計画の指紋の外)に出る。Console(45a531c4d)は設定表から 1 件でも読めると `Msp file path` を使わず、
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
5.    └─ metabolomix/pipeline/receipt.py  read_receipt()

`load_run` は読取専用で、このツール自身は `pipeline-run.json` を一切書き換えない。
`status="running"`のときだけ`same_process`でworker identityの生存を確認するが、
確認できなくても`status`フィールドは変えず`observed_health`という別軸で返す。
単位ごとの状態(`units`)と、受領があればその要約(解析フォルダの場所・照合した `.arf`・
次の呼び出し `load_dataset(pipeline_path=…, polarity=…)`)も返す。旧 v1 の run は
`legacy: true` で読めるだけ(再開・取消は `PIPELINE_RUN_LEGACY`)。

## pipeline_resume

1. metabolomix/tools/pipeline_tools.py  pipeline_resume()
2. └─ metabolomix/pipeline/service.py  resume_pipeline()
3.    ├─ metabolomix/pipeline/recovery.py  prepare_resume()
4.    │  └─ metabolomix/pipeline/request.py  merge_updates()
5.    ├─ [直前が needs_input のときだけ] metabolomix/pipeline/service.py  _recheck_manifests()
6.    │  └─ metabolomix/pipeline/service.py  _precheck_manifest()
7.    ├─ [status=="planned"のときだけ] metabolomix/pipeline/store.py  load_run()
8.    └─ [status=="planned"のときだけ] metabolomix/pipeline/service.py  _launch_and_await()
9.       ├─ metabolomix/pipeline/service.py  launch_pipeline_worker()
10.      │  └─ metabolomix/core/process_control.py  launch_detached()
11.      └─ metabolomix/pipeline/service.py  _await_launch_handshake()

受付の事前検査で `needs_input` に止まった run の resume は、起動する前に保存済みの要求と
入力計画で単位ごとの事前検査(`_precheck_manifest`)をやり直す。直っていなければ
`needs_input` のまま保存し直し、worker は起動しない(`launch: {"launched": false}`)。

`prepare_resume`は新しいattempt/revisionを用意し、どのstageを`pending`へ戻すかを
決めて保存するだけで、それ自体はworkerを起動しない(no-op resumeで無駄な起動を
避ける)。`resume_pipeline`は`prepare_resume`の戻り値`status`が`"planned"`の
ときだけ、改めて`load_run`でlaunch直前の記録(前回workerの残骸identityを含みうる)を
読み直し、それを基準として`pipeline_run`と共有する`_launch_and_await`を呼ぶ
(`launch_pipeline_worker`と`_await_launch_handshake`を1つにまとめた中間関数。
基準を「launch直前」に揃えないと、前回workerの残骸identityを新規launchの
証拠と誤認しうる)。`rerun_upstream=True`を明示しない限りConsole実行はやり直さない。
訂正できるのは試料の対応表(`sample_manifest`、2 極性なら `units.<極性>.sample_manifest`)
だけで、戻すのは受領の工程だけ。`rerun_upstream=True` は失敗した単位の 3 工程と受領を戻す。
ただし worker が生きていて(`_worker_health` が `ok`)上流が `running` の単位は、その
Console を捨てない——やり直しの対象にも、`UPSTREAM_RERUN_REQUIRED` の理由にも数えない
(稼働中の Console の最中に届いたシートの訂正は受け付け、受領だけを戻す)。worker が
居なくなった run の `running` は、今までどおり `rerun_upstream=true` を求める。

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
5.       ├─ [単位（極性）ごと、要求の順] metabolomix/pipeline/service.py  _handle_prepare_input()
6.       │  └─ metabolomix/pipeline/inputs.py  stage_inputs()
7.       ├─ [単位ごと] metabolomix/pipeline/service.py  _handle_upstream()
8.       │  └─ metabolomix/console/execution.py  supervise()
9.       ├─ [単位ごと] metabolomix/pipeline/service.py  _handle_validate_outputs()
10.      │  ├─ metabolomix/mztab/loading.py  job_gate()
11.      │  └─ metabolomix/pipeline/service.py  _matching_alignment()
12.      │     ├─ metabolomix/mztab/evidence.py  first_matching_arfs()
13.      │     │  └─ [候補が合わないとき] metabolomix/dataset/alignment_files.py  batch_arfs_newest_first()
14.      │     └─ metabolomix/mztab/evidence.py  matching_arfs()
15.      └─ metabolomix/pipeline/service.py  _handle_receipt()
16.         └─ metabolomix/pipeline/receipt.py  write_receipt()
17.   └─ metabolomix/pipeline/engine.py  finish_success()
18.      └─ metabolomix/pipeline/engine.py  final_status()

工程の計画(`engine.build_stages`)は単位ごとの `prepare_input:<極性>` →
`upstream:<極性>` → `validate_outputs:<極性>` を要求の順に並べ、最後に `receipt` を
1 つ置く。手順4の `_run_stage_loop` はこの計画順で回す(`record["stages"]` の挿入順では
ない)。ある単位が失敗したら、その単位の残りの工程は飛ばしてもう片方の単位へ進む
(Console を同時に 2 本走らせない)。validated の単位が 1 つも無ければ受領は書かない。

`_run_stage_loop`は**各stageの直前に最新のrun recordを読み直す**。要求版
(`_request_identity`＝revision/content_hash/saved_path)が変わっているか、
このpassで既に通過したstageが外から`pending`へ戻されている
(`_resume_reset_seen`)場合は、その場でpassをやり直す。走っている最中に届いた
`pipeline_resume`は——`prepare_resume`が稼働中workerのstatusを動かさないため
`resume_pipeline`が新しいworkerを起こさない——このworkerが取り込むしかない。
組み直しは有限回(`_REQUEST_RELOAD_MAX_PASSES`)までで、超えたら
`REQUEST_REVISION_CHURN`で止める。

`_handle_upstream`は`execution.supervise`を直接呼ぶ(Console用workerを二重起動しない)。
`console_job_path`は`analysis-job.json`そのものへの、pipeline_root 相対のパス
(`run_dir = Path(job_path).parent`)として単位の `upstream` に記録する。run_dirは
**単位とattemptごとに分ける**(`<pipeline_root>/<極性>/console/attempt-NNNN/`)——固定パス
にすると`rerun_upstream=true`の再実行が前回の終了証跡・ログ・実行前スナップショットを
上書きしてしまう。Console の `-i` は `<pipeline_root>/<極性>/input/`(配置した raw)で、
MS-DIAL はここへ `.arf` などを書く——このフォルダがそのまま解析フォルダになる。

`_handle_validate_outputs` は正準 mzTab-M(`job_gate`。ジョブが完了していなければ失敗)を
確かめ、解析フォルダの `.arf` のうち mzTab-M と同じアラインメントだと確かめられたもの
(`mztab.evidence.matching_arfs`)がちょうど 1 本なら `outputs.alignment` に
`{"arf", "arf2"}` を記録する。候補の選び方は `dataset/opening._matching_arf` と共有の
`mztab.evidence.first_matching_arfs`: まずジョブの `peak_matrix_source`(と mzTab-M の隣)、
どれも通らなければ解析フォルダのバッチを `AlignmentResult_<時刻>` の新しい順に 1 本ずつ試し、
最初に通ったところで止める——`rerun_upstream=true` の後は前の実行のバッチが同じフォルダに
残るので、全バッチを並べると「同じアラインメントが 2 本」になる。合わなければ失敗にせず警告
`ANALYSIS_FOLDER_UNMATCHED`(`source="mztab"` なら読める)で、run は `partial` になる。

`finish_success` の最終判定(`final_status`): 全単位が validated・受領が成功・全単位の
解析フォルダが照合できた → `completed`。validated が 1 つ以上で欠けがある → `partial`。
validated が 0 → `failed`。worker の中の `needs_input` は無い(受付の事前検査だけが保存する)。

受領(`results/receipt.json`、`lipidmix.pipeline-receipt.v1`)は単位ごとの解析フォルダ・
照合した `.arf`・メソッドとライブラリ(研究室の資産はファイル名と sha256 だけ)・試料の
対応表・検証の結果と、次の呼び出し(`load_dataset(pipeline_path=…, polarity=…)`)を載せる。
