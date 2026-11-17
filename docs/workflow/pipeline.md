# 生データフォルダpipeline呼び出し連鎖

生データから MS-DIAL の解析フォルダを作る入口(GUI の入口との関係は [index.md](index.md) の「一気通貫の順序」)。raw を run の
中へ配置し、極性ごとに MS-DIAL Console を回し、解析フォルダを mzTab-M と照合して受領(`results/receipt.json`)を書く。統計は
計算しない。メソッドファイルを選ぶ・作る `console_method_candidates` / `console_method_template` と、受付が読む試料の表
(`sample-manifest.tsv`)を ISA / mwTab から作る `pipeline_sample_sheet` もここに置く。

MCP 層は `metabolomix/tools/pipeline_tools.py` と `metabolomix/tools/console_tools.py`。実体は `metabolomix/pipeline/service.py`
(受付・工程 handler)と `metabolomix/pipeline/recovery.py`(読取専用の status・取消・再開判定)。工程は独立した worker
プロセス(`metabolomix/pipeline/worker.py`)が `metabolomix/pipeline/engine.py::run_engine()` で進めるので、MCP 呼び出しの
連鎖と worker の連鎖を分けて書く。

## pipeline_run

1. metabolomix/tools/pipeline_tools.py  pipeline_run()
2. └─ metabolomix/pipeline/service.py  start_pipeline()
3.    ├─ metabolomix/pipeline/service.py  _prepare_run()
4.    │  ├─ metabolomix/pipeline/request.py  resolve_request()
5.    │  ├─ [単位（極性）ごと] metabolomix/pipeline/inputs.py  inspect_inputs()
6.    │  │  └─ metabolomix/pipeline/inputs.py  _assert_open_format_inputs()
7.    │  │     └─ metabolomix/console/raw_integrity.py  broken_open_format_files()
8.    │  ├─ [単位（極性）ごと] metabolomix/pipeline/service.py  _precheck_manifest()
9.    │  │  └─ metabolomix/analysis/sample_manifest.py  parse_manifest()
10.   │  ├─ metabolomix/pipeline/store.py  find_or_create_run()
11.   │  └─ [manifest不正を検出した場合] metabolomix/pipeline/service.py  _mark_needs_input_without_launch()
12.   ├─ [manifest不正が無かった場合のみ] metabolomix/pipeline/store.py  load_run()
13.   ├─ [manifest不正あり、または既存run(活動中/completed/失敗/取消/部分完了)を再利用した場合] metabolomix/pipeline/service.py  _dispatch_summary()
14.   ├─ [新規作成直後の`planned` runのときだけ] metabolomix/pipeline/service.py  _launch_and_await()
15.   │  ├─ metabolomix/pipeline/service.py  launch_pipeline_worker()
16.   │  │  └─ metabolomix/core/process_control.py  launch_detached()
17.   │  └─ metabolomix/pipeline/service.py  _await_launch_handshake()
18.   └─ metabolomix/pipeline/service.py  _dispatch_summary()
19.      └─ metabolomix/pipeline/service.py  _launch_warnings()

手順 6-7: 開いた形式(mzML・imzML・mzXML)の主ファイルが根の閉じタグ(`</indexedmzML>` / `</mzML>` / `</mzXML>`)で
終わるかと、サイズが下限(`raw_integrity.MIN_OPEN_FORMAT_BYTES`、32 KiB)以上かを確かめる(imzML は閉じタグだけ。`.ibd` は
見ない)。読むのは末尾の数 KB だけで、XML 全体の検証もベンダーの形式の検査もしない。切れた・小さすぎるファイルは名指しして
run を作る前に `PIPELINE_INPUT_INVALID`(`details.files` に `{"name", "size", "reason"}`、`reason` は `truncated` /
`too_small` / `unreadable`)で断る。

事前検査(手順 8-9)が既知の不正入力(壊れた実験情報シートなど)を見つけたら、`find_or_create_run` で run は作るが
worker は起動しない(手順 11 で `needs_input` として保存し、手順 12 を飛ばして手順 13 の発送 receipt を返す)。
manifest が正常でも、`find_or_create_run`(手順 10)が既存 run(活動中・completed・失敗・取消・部分完了)を再利用したときは、
手順 12 で読み直した status が `planned` でないので起動しない(手順 13)。起動するのは新規に作られた `status="planned"` の
run だけ。既存 run を進めるのは `pipeline_resume`。

要求は `pipeline-request.v3`(`request.resolve_request`)で、1 極性の形も 2 極性の `units` の形も単位の並びに正規化する。
入力計画(`inspect_inputs`)と事前検査は単位ごとに回し、run(`pipeline-run.v3`)は単位ごとの記録 `record["units"][<極性>]` を持つ。
- 下流のキー(`target` / `preprocess` / `comparisons`)はほかの知らないキーと同じく `PIPELINE_REQUEST_INVALID`
  (`details.unknown_keys`。下流のキーには文面と `details.replacements` で `load_dataset` → `dataset_preprocess` /
  `dataset_compare` / `dataset_pca` を示す。`request._unknown_keys`)。元フォルダの `analysis-request.json` にある下流のキーは
  読み飛ばして警告 `REQUEST_FILE_DOWNSTREAM_IGNORED`。
- 明示の `units` があると、`analysis-request.json` のトップの単位のキー(`polarity` `method_file` …)は使わず警告
  `REQUEST_FILE_UNIT_KEYS_IGNORED`。どちらの警告も `pipeline_plan` / `pipeline_run` の `warnings` に出る。
- 2 単位の `raw_dir` が同じフォルダ、または片方がもう片方の中(`.` を含む)なら `PIPELINE_REQUEST_INVALID`。

手順 14 の `_launch_and_await`(`pipeline_resume` と共有)は launch 直前に読んだ `baseline_record` を取り、
`launch_pipeline_worker`(起動)と `_await_launch_handshake`(handshake 待ち)を同じ基準の記録で呼ぶ(前回 worker の残骸の
identity を新しい起動の証拠と取り違えない)。`_await_launch_handshake` は短い上限だけ待ち、超えても起動失敗とせず再起動も
しない(`launch.handshake="not_confirmed"` で `pipeline_path` を返す。確認は `pipeline_status`)。

`launch_pipeline_worker` は `sys.executable` で `python -m metabolomix.pipeline.worker` を起こし、cwd はこの checkout
(`_REPO_ROOT`)に固定する。起動プロセスの stdout/stderr は `launch_detached` がログファイルへ結び(起動元が worker の
stdout を継承して待ち続けることは無い)、子へ渡すハンドルはログと NUL の 2 つだけ(`process_control._create_process` の
`PROC_THREAD_ATTRIBUTE_HANDLE_LIST`)。MCP サーバの標準入出力のパイプなどは中継・worker・Console に渡らない。

`launch_detached` は、すぐ終わる中継(`python -m metabolomix.core.detach_relay`。`process_control._launch_via_relay` が起こして
終わりを待つ)に worker を起こさせる。中継は worker の pid と identity(pid と生成時刻)を書いて終わるので、worker の親は残らず、
プロセスの木ごとの終了(`taskkill /T` と同じ)でも worker とその Console(worker が自分の Job で所有する)は残る。2 段(サーバ→
中継・中継→worker)とも同じ規則で切り離す(breakaway のフラグは Job が許すときだけ)。戻り値の pid・identity は worker のもので、
handshake(worker が `pipeline-run.json` に書く identity)と `pipeline_status` / `pipeline_cancel` / `pipeline_resume` の追跡に使う。

Job の制限は直近の 1 段しか読めない(`QueryInformationJobObject(NULL)`)。`DETACH_UNSUPPORTED` で断るのは、直近の Job が
KILL_ON_JOB_CLOSE 付きで breakaway を許さないときだけ。breakaway を許さない Job(`inherited_job`)やその外側の Job に残った
中継・worker は、その Job が閉じる・`TerminateJobObject` されると巻き込まれる。中継は worker を `IsProcessInJob(worker, NULL)`
(何段目の Job でも捕まえる)で確かめ、`launch_detached` の戻り値の `worker_in_job` に返す(`relay_breakaway` は中継の直近の
Job から見た値)。警告は `engine.worker_in_job_warning` が決める:
- 中継が直近の Job を抜けた(`silent_breakaway` / `breakaway_ok`)のに worker が Job にいる(中継から見えない外側の段)→
  `_launch_warnings`(手順 19)が `pipeline_run` / `pipeline_resume` の `warnings` に `WORKER_IN_JOB`(`details.breakaway` /
  `relay_breakaway`)を足す。
- 中継の直近の Job が `inherited_job`(breakaway は許さないが閉じても止めない)→ 警告しない。外側に KILL_ON_JOB_CLOSE の Job が
  あれば、それが閉じたとき worker と Console は止まる(`0x3000 → 0x0 → 0x2000` の入れ子では、中継は 0x3000 を黙って抜け、worker
  は 0x0 と外側の 0x2000 に残る)ので、`_launch_and_await`(手順 14)が起動の結果に注記(`engine.worker_outer_jobs_note`)
  `launch.outer_jobs: "unchecked"` と `launch.outer_jobs_note`(外側の段は確かめていない。止まったら `pipeline_status` →
  `pipeline_resume`)を足す。

worker 自身の確かめは「workerが1回分を進める経路」の節(手順 2・4)。

## pipeline_plan

1. metabolomix/tools/pipeline_tools.py  pipeline_plan()
2. └─ metabolomix/pipeline/service.py  plan_pipeline()
3.    ├─ metabolomix/pipeline/service.py  _prepare_run()
4.    └─ metabolomix/pipeline/service.py  _dispatch_summary()
5.       └─ metabolomix/pipeline/service.py  _resolved_settings()

`pipeline_run` と事前検査まで同じ経路(`_prepare_run`)を通り、`launch_pipeline_worker` を呼ばない(検査・計画・不足情報の
保存だけ)。

手順 5 の `_resolved_settings` は `plan_pipeline` の発送 receipt にだけ足す(`_dispatch_summary(include_resolved=True)`)。
載せるのは解決済みの method(`source_path` / `sha256`)・LBM・polarity(`value` / `source`)・`library_mode`・MSP(ファイル名・
sha256・出どころ。置き場所は載せない)・`alignment_reference` だけで、`raw_stat`・`entries`・`overrides` などの中間データは
`record["inputs"]` に残す。method には結果を大きく変えるキー(`metabolomix/console/method_file.py` の `KEY_PARAM_KEYS`)の実効値
`key_params`・宣言の無いキー `key_params_undeclared`・要求の上書き `user_overrides` も載せる。`pipeline_run` の receipt は最小のまま。

要求の `method_overrides`(単位のキー。`{Console のキー: 値}`): `metabolomix/pipeline/request.py` の `resolve_request()` が
`metabolomix/console/method_keys.py` の `normalize_overrides()` でキーと型を検査し(照合表は MsdialWorkbench afd5f95 の
`ConfigParser.cs` から。Console は知らないキーと型に合わない値を黙って捨てる。現用の Console のビルドが読まない 7 キーは
`NOT_IN_CURRENT_BUILD` で断る)、`{表記の綴り: 文字列}` に正規化して要求に保存する
(要求の hash に入るので、上書きが変われば別の run。空の `{}` は省いたのと同じ)。`metabolomix/pipeline/inputs.py` の
`inspect_inputs()` はメソッドを読んだ直後に上書きを当て(原本に宣言のある別名も同じ値にする。`alias_keys()`)、以後の解決
(LBM・MSP・足切り)は上書き後の値を見る。計画が自分で書くキー(ライブラリ・パスの絶対化・アライメントの参照)と重なれば
`PIPELINE_REQUEST_INVALID`。上書きは計画の `method.overrides` に入り、`stage_inputs()` が実効メソッドへ書く(原本は変えない)。
計画の指紋にも入る。

`library_mode="msp_only"`: LBM を解決せず(`Lbm file path` は空)、研究室の参照ライブラリ(.msp)を `msp_file` か極性の設定
(`[library] msp_positive` / `msp_negative`)から決めて実効メソッドの `Msp file path` に書く。メソッドが宣言する Text DB と
注釈器設定表は空にする。dot product の足切りは下の `msp_then_lbm` と同じく平方根の側のキーで書き、receipt の
`resolved.dot_product_cutoffs` に出す。MSP の sha256 は入力計画と計画の指紋に入る(MSP が変われば `pipeline_run` は別の run を
作る)。上流をやり直さない `pipeline_resume` では `metabolomix/pipeline/inputs.py` の `verify_inputs()` が再検査する(変わって
いれば `INPUT_CHANGED`、`which: "msp"`)。

`library_mode="msp_then_lbm"`(脂質の run 用): 研究室 MSP を優先し、LBM で補う。MSP は `msp_only` と同じ規則で決め、LBM は
`[library] lbm` **だけ**から決める(メソッドの宣言・実行体のフォルダ・ビルド生成物は見ない。無ければ `LBM_NOT_CONFIGURED`、
指す先が無ければ `LBM_NOT_FOUND`)。実効メソッドには次を書く:
- `Msp file path` は空。`MSP annotator settings file path` に `inputs/msp-annotators.tsv`(`stage_inputs` が書く。列
  `annotator_id` / `msp_file_path` / `priority`、研究室 MSP を priority 2 で 1 行)。メソッドが宣言する MSP 設定表の別名も
  同じ設定表へ向ける(Console は別名を含めて最初の行を読む)。Text DB と Text の設定表は空。
- `Lbm file path` に LBM、`Use retention information for LBM-based annotation scoring` / `filtering` に `False`。
- dot product の足切り(weighted / simple / reverse)は MSP と LBM に同じ値を `Square root of ... dot product cutoff for
  MSP-based annotation` / `...for LBM-based annotation` で書き、宣言された接頭辞の無いキー(`Weighted dot product cutoff for
  MSP-based annotation` 等)は空にする(Console は接頭辞の無いキーを**二乗の**足切りとして読む。0.6 は実効で 0.775)。値は
  平方根の側の宣言 → 接頭辞の無いキーの宣言(平方根の側の値と読む)→ 既定 0.6 / 0.6 / 0.8 の順で決め、空にした宣言(原本が値を
  宣言していたもの)は `removed_declarations` に載る。
  実効の足切りは receipt の `resolved.dot_product_cutoffs`(`weighted` / `simple` / `reverse` / `scale: "square_root"` /
  種類ごとの `source`: `root_key` / `unprefixed_key_read_as_root` / `default`。msp_only と msp_then_lbm。計画の指紋の外)。

Console(45a531c4d)は設定表から 1 件でも読めると `Msp file path` を使わず、LBM の注釈器の priority を 1 に固定するので、代表の
同定は同じ証拠の水準なら MSP、MSP に MS/MS の照合が無く LBM にあるときだけ LBM になる。LBM の sha256 も計画の指紋に入り、
`verify_inputs()` が再検査する(`which: "lbm"`。details はファイル名だけ)。receipt の `resolved.lbm` は `file` /
`sha256` / `source`、`resolved.msp` は `priority`(2)も持つ。設定表の sha256 は配置時に snapshot(`msp.annotator_table_sha256`)へ
固定し、再開時に実効メソッドと同じく照合する(変わっていれば `STAGED_INPUT_MISMATCH`)。`pipeline_status(include_details=true)`
の record の写しでは、LBM の `path`・`Lbm file path` の上書き・絶対パスの設定表の上書きもファイル名だけにする。mzTab の読み込みは
設定表の MSP で InChIKey を補う(同定は `MspDB_<n>_<MSP の stem>:<名前>`)。

アライメントの参照ファイル(要求の `alignment_reference`、既定 `"median_injection"`): 実効メソッドの
`Alignment reference file ID` に、全試料(`include` を問わない)の注入順の中央の試料(偶数なら小さい方)の file ID(入力
フォルダの主ファイルを `sorted()` した順の 0 始まりの添字)を書く。`_prepare_run` は実験情報シートを `inspect_inputs` より先に
解して行を渡す。`injection_order` が 1 件でも欠ける・シートが無いときはファイル名順の中央(receipt の `injection_order_source`
が `"filename"`)。`"method"` は上書きせず、メソッドの宣言を `method_value` に残す。`overrides` は計画の指紋の外なので、要求の
hash も完了済み run の再利用も変わらない。ただし古い要求を再送して新しい run を作ると、実効メソッドの参照ファイルが以前
(0)と違う(中央)。`alignment_reference` の既定は要求の hash の外なので、この設定の無い計画で完了した run がある元フォルダへ
既定のまま要求を出すと、その run を再利用する(上書きしたメソッドは使われていない)。このとき receipt の
`resolved.alignment_reference` は `null` で、参照ファイルが中央でないことの印。中央の参照でやり直すには `alignment_reference` を
明示するか新しい run を作る。

## pipeline_status

1. metabolomix/tools/pipeline_tools.py  pipeline_status()
2. └─ metabolomix/pipeline/recovery.py  read_status()
3.    └─ metabolomix/pipeline/store.py  load_run()
4.    └─ metabolomix/core/process_control.py  same_process()
5.    └─ metabolomix/pipeline/receipt.py  read_receipt()
6.    └─ [失敗・取消で終わった工程があるとき] metabolomix/pipeline/recovery.py  _stage_errors()
7.       └─ metabolomix/pipeline/receipt.py  stage_error_summary()

受領の要約の `method_params` は `metabolomix/pipeline/receipt.py` の `unit_method_params()` が単位の実効メソッドを読み直して
作る(受領より前の run も同じに見える)。`load_dataset` の `processing` も同じもの。

`load_run` は読取専用で、このツールは `pipeline-run.json` を書き換えない。`status="running"` のときだけ `same_process` で
worker identity の生存を確かめ、結果は `status` を変えずに `observed_health` で返す。単位ごとの状態(`units`)と、受領があれば
その要約(解析フォルダの場所・照合した `.arf`・次の呼び出し `load_dataset(pipeline_path=…, polarity=…)`)も返す。
- 走っている工程が機械全体の Console の順番を待っていれば `waiting`(`{<stage_id>: {"reason": "console_slot", "since",
  "message", "holder": {"pipeline_id", "unit", "since"}}}`。工程の status は `running` のまま)。待っていた worker が居なければ
  (`worker_missing`)、記録に残った待ちは `waiting` でなく `interrupted_while_waiting` に出し、Console をまだ起こしていない単位
  (`recovery._console_never_started`: 待ちの理由が `console_slot` で、今の attempt の `analysis-job.json`〔順番が来てから
  Console を起こす前に `store.console_attempt_dir` の下へ書く〕が無い。単位の `console_job_path` は見ない)を
  `recovery_hint.console_not_started_units` に挙げる。
- 失敗・取消で終わった工程があれば `errors`(`{<stage_id>: {"code", "message", "cause"?, "suspect_files"?,
  "n_suspect_files"?}}`。手順 6-7。`stage_statuses` は状態の文字列だけ)。要約は `receipt.stage_error_summary` で、受領の単位の
  `error` と同じ形: `MSDIAL_EXECUTION_FAILED` の `details.cause` と `details.diagnosis.suspect_files`、配置の工程で断った
  `PIPELINE_INPUT_INVALID` の `details.files` を `cause` / `suspect_files`(`[{"name", "reason"}]`。10 件まで、超えたら
  `n_suspect_files`)にそろえ、文面は 600 字で切る。ログの行(`details.diagnosis.log_lines`)は載せない(`include_details=true`
  の記録全体にある)。
- 旧 v1 の run は `legacy: true` で読めるだけ(再開・取消は `PIPELINE_RUN_LEGACY`)。

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
12.   └─ [起動したとき] metabolomix/pipeline/service.py  _launch_warnings()

受付の事前検査で `needs_input` に止まった run は、起動の前に保存済みの要求と入力計画で単位ごとの事前検査
(`_precheck_manifest`)をやり直す。直っていなければ `needs_input` のまま保存し直し、worker を起動しない
(`launch: {"launched": false}`)。

`prepare_resume` は新しい attempt / revision を用意し、どの stage を `pending` へ戻すかを決めて保存するだけで、worker を起動
しない。`resume_pipeline` は `prepare_resume` の戻り値の `status` が `"planned"` のときだけ、`load_run` で launch 直前の記録
(前回 worker の残骸の identity を含みうる)を読み直し、それを基準に `_launch_and_await` を呼ぶ(pipeline_run の節)。

- `rerun_upstream=True` を明示しない限り Console はやり直さない。訂正できるのは試料の対応表(`sample_manifest`、2 極性なら
  `units.<極性>.sample_manifest`)だけで、戻すのは受領の工程だけ。`rerun_upstream=True` は失敗した単位の 3 工程と受領を戻す。
- `updates` の知らないキー(消えた下流のキーを含む)は `PIPELINE_REQUEST_INVALID`(下流のキーなら `details.replacements`)、
  上流の条件のキー(`method_file` ほか)は `NEW_PIPELINE_REQUIRED`。
- worker が生きていて(`_worker_health` が `ok`)上流が `running` の単位は、その Console を捨てない: やり直しの対象にも
  `UPSTREAM_RERUN_REQUIRED` の理由にも数えない(`rerun_upstream=True` で外した単位は戻り値の `held_units` と `hint`)。稼働中の
  Console の最中に届いたシートの訂正は受け付け、受領だけを戻す。
- worker が居なくなった run の `running` は `rerun_upstream=true` を求める。ただし Console の順番待ちの間に worker が消え、
  Console をまだ起こしていない単位(`_console_never_started`。判定は pipeline_status の節)は数えず、`upstream` と
  `validate_outputs` だけを戻して `rerun_upstream` なしで最初の Console から続ける(戻り値の `console_not_started_units`)。
- 同じ `request_id` の再送は、`prepare_resume` が `resume_log[request_id]` に残した内訳(`not_retried_units` / `held_units` /
  `console_not_started_units`。空なら持たない)を `_resume_receipt` に渡し、1 回目と同じ内訳を `reused: true` で返す。

## pipeline_cancel

1. metabolomix/tools/pipeline_tools.py  pipeline_cancel()
2. └─ metabolomix/pipeline/recovery.py  request_cancel()
3.    └─ metabolomix/pipeline/store.py  load_run()
4.    └─ metabolomix/pipeline/engine.py  cancel_request_path()

保存するのは協調的な取消フラグ(小さな JSON)だけで、プロセスを直接終了させない。**raw データは削除しない。** Console の停止は
`execution.supervise` の `cancel_path` 監視、工程の停止は engine の stage 境界の確認が行う。旧 v1 の run(`pipeline-run.v1`)は
`PIPELINE_RUN_LEGACY` で断る(再開も同じ)。受理(このツールの成功)と停止の確定(`cancelled`)は別の状態で、後者は
`pipeline_status` で確かめる。

## pipeline_sample_sheet

前提: なし(データセットも run も要らない)。状態変更: なし(セッションに書かない)。raw フォルダに `sample-manifest.tsv`
(`output_name`)を 1 つ書く。raw ごとに取得時刻の見出しを開く(数百のファイルやネットワークの共有では数十秒)ので、
`handle_dispatch.ALWAYS_THREADED` でワーカースレッドに回す(その間もイベントループは `job_status` などに答える)。

1. metabolomix/tools/sample_sheet_tools.py  pipeline_sample_sheet()
2. └─ metabolomix/pipeline/inputs.py  list_primary_names()
3. └─ metabolomix/rawmeta/sample_table.py  build_sample_sheet()
4.    └─ metabolomix/rawmeta/sample_table.py  _sheet_factor_names()
5.    └─ metabolomix/rawmeta/sample_table.py  _load_sample_table()
6.       └─ metabolomix/rawmeta/sample_table.py  read_mwtab_samples()
7.       └─ metabolomix/rawmeta/sample_table.py  read_table()
8.    └─ metabolomix/rawmeta/sample_table.py  resolve_column()
9.    └─ metabolomix/rawmeta/sample_table.py  _table_rows()
10.    └─ metabolomix/rawmeta/sample_table.py  _assay_map()
11.    └─ metabolomix/rawmeta/sample_table.py  _file_index()
12.    └─ metabolomix/rawmeta/sample_table.py  factor_group()
13.    └─ [生体試料のどれも群が決まらないとき] metabolomix/rawmeta/sample_table.py  _unmatched_error()
14.    └─ metabolomix/rawmeta/acquisition.py  acquisition_order()
15.       └─ metabolomix/rawmeta/acquisition.py  read_acquisition_times()
16. └─ metabolomix/analysis/sample_manifest.py  parse_manifest()
17. └─ metabolomix/core/atomic_io.py  atomic_write_bytes()
18. └─ metabolomix/tools/sample_sheet_tools.py  _notes()
19. └─ metabolomix/rawmeta/acquisition.py  gaps_over()

raw の一覧は受付(`pipeline_plan` / `pipeline_run`)と同じ `list_primary_names`。各 raw を試料の表の行へ対応づける(ISA は
assay 表の Raw Spectral Data File と Derived Spectral Data File、mwTab は SUBJECT_SAMPLE_FACTORS の RAW_FILE_NAME、どちらも
無ければ試料名と raw の語幹)。表のファイル名と合わなければ語幹で引く(`_file_index`。名前・語幹が別々の試料に当たるファイルは
対応に入れない)。同じ Sample Name の行どうしで因子(と役割を決める列)の値が食い違う試料は行を決められない(`_table_rows`)。
群は選んだ因子の値の組(`factor_group`)。群の決まらない生体試料は群 `unassigned`・`include=false`。生体試料のどれも群が
決まらなければ書かずに `SAMPLE_TABLE_INVALID`(`_unmatched_error`)。因子の列名がシートの既知の列と同じなら `factor_<名>` に
改める(`_sheet_factor_names`)。

注入順は raw の取得時刻の順(`acquisition_order`。スペクトルは読まず、ファイルの見出しのメタデータだけ)。読めない raw は注入順を
空にしてファイル名と理由を返し、読めた raw の順位は残す。1 つでも読めれば schema の行の直後に宣言の行
`# injection_order_source = raw_acquisition_time` を書き、`parse_manifest` が注入順の出所の記録に `derived_from` を足す(試料の
表と統計の出所が `raw_acquisition_time` になる)。一部だけ空のときは空の行にも `derived_from` を付け
(`sample_manifest.sheet_order_unknown`)、`_merge_row` と統計エンジンはその空欄を `.mddata`・mzTab-M の番号で埋めない。
`.wiff2` の取得時刻は同じ語幹の `.wiff`(SCIEX OS が並べて書く OLE のファイル。`acquisition._read_wiff2`)から読む。

試料の表が UTF-8 で読めなければ `SAMPLE_TABLE_INVALID`。書く前に受付と同じ厳密な読み手(`parse_manifest`)へ一時ファイルで
通し、通ったものだけを原子的に書く。既にあれば `overwrite=True` が無い限り `SAMPLE_SHEET_EXISTS`。`_notes` は、順の欠けた
生体試料か QC が統計に入る(`include=true`)ときだけ「ドリフト補正と注入順の偏りの判定は行われない」と書く(前処理の `_select`
が `include=false` の試料を先に外す)。

## console_method_candidates

前提: なし。状態変更: なし(読取専用)。

1. metabolomix/tools/console_tools.py  console_method_candidates()
2. └─ [polarity を渡したとき] metabolomix/tools/console_tools.py  _polarity_or_error()
3.    └─ metabolomix/core/polarity_names.py  normalize_polarity()
4. └─ metabolomix/console/method_file.py  discover_method_candidates()
5.    └─ metabolomix/console/method_file.py  method_search_dirs()
6.    └─ metabolomix/console/method_file.py  scan_dir_for_method_files()
7.    └─ metabolomix/console/method_file.py  past_run_method_files()
8.       └─ metabolomix/console/method_file.py  read_method_keys()
9.    └─ metabolomix/console/method_file.py  extract_key_params()
10. └─ metabolomix/tools/console_tools.py  _capped_candidates()
11.    └─ metabolomix/tools/console_tools.py  _candidate_payload()
12. └─ metabolomix/tools/console_tools.py  _candidates_next_hint()

`polarity` は `normalize_polarity` で `positive` / `negative` に読み替える(`pos` / `neg` / `+` / `-` も可、大文字小文字を
問わない)。読めない値は `POLARITY_INVALID`。

MS-DIAL GUI が解析のたびに保存する `<project>_param_<終了時刻>.txt` を、`dataset_root` 直下 → 兄弟フォルダ直下 →
`search_dirs`(`method_search_dirs`)の順に探す。**再帰はしない**。`dataset_root/runs/*/analysis-job.json` の
`software.method_file`(`past_run_method_files`。Console 単独実行のジョブが使ったメソッド)も候補にする。pipeline の run の
ジョブ(`runs/pipeline_*/<極性>/console/attempt-*/analysis-job.json`)は 1 段深いので拾わない。極性で候補を落とさない:
`polarity` を渡すと各候補に `usable`(`direct` / `needs_polarity_conversion`)が付き、別極性は `console_method_template` を
通す。候補が `KEY_PARAMS_MAX_CANDIDATES` を超えると `key_params` を付けず、`MAX_REPORTED_CANDIDATES` で切って `truncated: true`
(並びは `direct` 優先・新しい順)。`next` は `pipeline_plan` の `request.method_file` へ渡す道を言う。

## console_method_template

前提: なし。状態変更: `out_path` にメソッドファイルを書く(セッションは変えない)。

1. metabolomix/tools/console_tools.py  console_method_template()
2. └─ metabolomix/tools/console_tools.py  _polarity_or_error()
3. └─ metabolomix/tools/console_tools.py  _library_mode_error()
4. └─ [based_on を省いたとき] metabolomix/console/method_file.py  discover_method_candidates()
5. └─ metabolomix/tools/console_tools.py  _looks_like_method_text()
6. └─ metabolomix/console/method_file.py  read_method_keys()
7. └─ metabolomix/console/runner.py  get_exe_path()
8. └─ metabolomix/tools/console_tools.py  _resolve_libraries()
9.    └─ metabolomix/console/method_file.py  relative_path_overrides()
10.   ├─ [library_mode="msp_only"] metabolomix/tools/console_tools.py  _msp_setting()
11.   ├─ [library_mode="msp_only"] metabolomix/console/method_file.py  resolve_msp()
12.   ├─ [library_mode="msp_only"] metabolomix/console/method_file.py  msp_only_overrides()
13.   ├─ [library_mode="auto"] metabolomix/tools/console_tools.py  _lbm_setting()
14.   └─ [library_mode="auto"] metabolomix/console/method_file.py  resolve_lbm()
15. └─ metabolomix/tools/console_tools.py  _non_ascii_override_error()
16. └─ metabolomix/console/method_file.py  write_effective_method_file()

既存のパラメータから別極性のメソッドファイルを作る。その極性で GUI 実行が一度でもあれば、`method_file` を省いた
`pipeline_plan` で足りる。`polarity` は `console_method_candidates` と同じく読み替え、読めない値は `POLARITY_INVALID`。
- `based_on` を省くと `discover_method_candidates` を**極性で絞らずに**呼び、候補が 1 件ならそれ、複数なら
  `METHOD_FILE_CHOICE_REQUIRED`(最新を黙って採らない)。
- 差し替えるのは `Ion mode` と `Searched adduct ions`(`STANDARD_ADDUCTS`)だけで、検出・アライメント条件は元のまま。
- ASCII テキストでない土台(`.mdproject` / `.mddata` など)は `_looks_like_method_text` が先頭 8192 バイトだけを見て
  `METHOD_FILE_NOT_TEXT` で断る。
- 実行体(`get_exe_path`)はビルド生成物・実行体のフォルダからの LBM の推定にだけ使い、未設定でも止めない(設定ファイルが
  読めない `CONFIG_INVALID` だけは止める)。
- `_resolve_libraries` は相対で宣言されたパスを原本基準の絶対パスへ直し、`auto` では空の `Lbm file path` を pipeline と同じ
  解決順で埋め、`msp_only` では LBM を空にして作りたい極性の研究室 MSP を書き、Text DB と注釈器設定表の宣言を空にする。研究室の
  LBM・MSP は戻り値・エラーにファイル名だけを出す。
- 実効メソッドへ書く値に ASCII で表せない文字があれば `METHOD_ENCODING_UNSUPPORTED`(Console はメソッドを ASCII で読み、黙って
  同定 0 件になる)。
- 戻り値の `next` は `pipeline_plan(dataset_root=..., request={'method_file': ..., 'polarity': ...})`(`msp_only` ならその
  モードも引き継ぐ)。

## workerが1回分を進める経路(別プロセス。pipeline の5ツールとは別に進行する)

`pipeline_run` / `pipeline_resume` が起動した worker プロセスの経路。MCP 接続とグローバル session に触れない
(`metabolomix.core.session_state` / `metabolomix.core.mcp_core` / `metabolomix.tools.*` を import しない。
`tests/test_pipeline_engine.py` の AST テストが固定する)。

1. metabolomix/pipeline/worker.py  run_worker()
2. └─ metabolomix/core/process_control.py  current_job_state()
3. └─ metabolomix/pipeline/engine.py  run_engine()
4.    └─ metabolomix/pipeline/engine.py  _record_worker_identity()
5.    └─ metabolomix/pipeline/service.py  build_handlers()
6.    └─ metabolomix/pipeline/engine.py  _run_stage_loop()
7.       ├─ [単位（極性）ごと、要求の順] metabolomix/pipeline/service.py  _handle_prepare_input()
8.       │  └─ metabolomix/pipeline/inputs.py  stage_inputs()
9.       │     └─ metabolomix/pipeline/inputs.py  _assert_open_format_inputs()
10.      ├─ [単位ごと] metabolomix/pipeline/service.py  _handle_upstream()
11.      │  ├─ metabolomix/console/machine_slot.py  console_slot()
12.      │  │  ├─ metabolomix/core/process_control.py  wait_for_lock()
13.      │  │  └─ [順番を待つとき] metabolomix/pipeline/engine.py  report_stage_waiting()
14.      │  └─ metabolomix/pipeline/service.py  _run_console_unit()
15.      │     ├─ metabolomix/console/machine_slot.py  check_commit_memory()
16.      │     │  └─ metabolomix/core/process_control.py  commit_memory_status()
17.      │     ├─ metabolomix/console/execution.py  supervise()
18.      │     └─ [完了しなかったとき] metabolomix/console/failure_diagnosis.py  diagnose_failure()
19.      ├─ [単位ごと] metabolomix/pipeline/service.py  _handle_validate_outputs()
20.      │  ├─ metabolomix/mztab/loading.py  job_gate()
21.      │  └─ metabolomix/pipeline/service.py  _matching_alignment()
22.      │     ├─ metabolomix/mztab/evidence.py  first_matching_arfs()
23.      │     │  └─ [候補が合わないとき] metabolomix/dataset/alignment_files.py  batch_arfs_newest_first()
24.      │     └─ metabolomix/mztab/evidence.py  matching_arfs()
25.      └─ metabolomix/pipeline/service.py  _handle_receipt()
26.         └─ metabolomix/pipeline/receipt.py  write_receipt()
27.   └─ metabolomix/pipeline/engine.py  finish_success()
28.      └─ metabolomix/pipeline/engine.py  final_status()

手順 2・4: worker は自分が Windows の Job に入っているか(何段目でも)と直近の Job の制限を確かめ
(`process_control.current_job_state`)、run の記録の `worker.in_job` / `worker.nearest_job` に残す。閉じたら止めるか分からない
Job なら警告 `WORKER_IN_JOB`(run に 1 つ。`details.nearest_job`)、直近の Job が `inherited_job` なら警告せず注記
`worker.outer_jobs: "unchecked"` / `worker.outer_jobs_note` を残す。どちらも止めない。試験の harness は `run_engine` を直に
呼ぶので確かめない。
手順 9 は受付と同じ切れた mzML の検査を配置の前にもう一度行う(受付の後に切れたものも Console の前で止める)。

工程の計画(`engine.build_stages`)は単位ごとの `prepare_input:<極性>` → `upstream:<極性>` → `validate_outputs:<極性>` を
要求の順に並べ、最後に `receipt` を 1 つ置く。`_run_stage_loop`(手順 6)はこの計画順で回す(`record["stages"]` の挿入順では
ない)。ある単位が失敗したら、その単位の残りの工程を飛ばして次の単位へ進む(Console を同時に 2 本走らせない)。validated の
単位が 1 つも無ければ受領は書かない。

`_run_stage_loop` は**各 stage の直前に run の記録を読み直し**、要求の版(`_request_identity`＝revision / content_hash /
saved_path)が変わったか、この pass で通過した stage が外から `pending` へ戻された(`_resume_reset_seen`)ら pass をやり直す。
走っている最中の `pipeline_resume` は新しい worker を起こさない(`prepare_resume` が稼働中の worker の status を動かさない)ので、
この worker が取り込む。やり直しは `_REQUEST_RELOAD_MAX_PASSES` 回までで、超えたら `REQUEST_REVISION_CHURN` で止める。

`_handle_upstream` は `execution.supervise` を直接呼ぶ(Console 用の worker を二重に起こさない)。`console_job_path` は
`analysis-job.json` への pipeline_root 相対のパス(`run_dir = Path(job_path).parent`)で、単位の `upstream` に記録する。run_dir は
**単位と attempt ごとに分ける**(`<pipeline_root>/<極性>/console/attempt-NNNN/`。再実行は前回の終了証跡・ログ・実行前
スナップショットを上書きしない)。Console の `-i` は `<pipeline_root>/<極性>/input/`(配置した raw)で、MS-DIAL はここへ `.arf`
などを書く。このフォルダが解析フォルダになる。

**機械全体で Console は 1 本**: run の中では単位を 1 つずつ回すが、極性ごとに別の run を起こすと Console が並ぶ(脂質の run は
1 本で 20 GB 以上の commit を確保し、並ぶと Windows のエラー 1450 で落ちうる)。`_handle_upstream` は Console を起こす前に
`machine_slot.console_slot`(OS のファイルロック
`process_control.wait_for_lock`。持っている worker が死ねば OS が外す)で順番を取る。
- ロックは設定ファイルのフォルダの `data/.msdial-console.lock` で、同じ設定ファイルを読むサーバ・worker の間で取り合う(main
  ツリーのサーバと、`LIPIDMIX_CONFIG` で main ツリーの `lipidmix.local.toml` を指す worktree)。自分の設定ファイルを読む別の
  クローンは別のロックになり、Console が 2 本並びうる(`LIPIDMIX_CONFIG` で同じ設定ファイルを指せば 1 本)。
- 待つ間は engine が context に渡す `report_waiting` から `report_stage_waiting`(手順 13)で工程の `waiting` に書き(status は
  `running` のまま。`pipeline_status` の `waiting`)、順番が来たら消す。待つ間に取消が来たら Console を起こさず
  `CONSOLE_WAIT_CANCELLED` で失敗にする(run は `cancelled`)。
- 順番が来たら commit の空き(`GlobalMemoryStatusEx` の `ullAvailPageFile`)を測り(手順 15-16)、下限(設定
  `[msdial] min_free_commit_gb`、既定 8、0 で測らない)を下回れば警告 `LOW_COMMIT_MEMORY` を残す(止めない)。測った値は単位の
  `upstream.memory_before`、待った秒数は `upstream.console_slot_waited_s`。

Console が完了しなかったとき(終わり方が `exited` / `worker_lost`)は `failure_diagnosis.diagnose_failure`(手順 18)の見立てを
`MSDIAL_EXECUTION_FAILED` の `details.diagnosis` に、原因の分類を `details.cause` に載せる: 終了コード 1450
(`ERROR_NO_SYSTEM_RESOURCES`)かログのメモリ不足は `system_resources`、配置済みの mzML が切れていれば `broken_input`、
`msdial.log` のエラーの行(無ければ最初のエラーより前に最後に名前の出た行)に計測ファイルの名前があれば `input_suspected`
(見立てで確定ではない)。返すログの行は絶対パスをファイル名にする: ドライブ・UNC・`file:` の URL で始まるパスの頭から最後の
区切りまでを消す(フォルダ名の空白で切らない)。原因と疑わしいファイルは `pipeline_status` の `errors` と受領の単位の `error` にも
出る(`receipt.stage_error_summary`)。

`_handle_validate_outputs` は正準 mzTab-M(`job_gate`。ジョブが完了していなければ失敗)を確かめ、解析フォルダの `.arf` のうち
mzTab-M と同じアラインメントと確かめられたもの(`mztab.evidence.matching_arfs`)がちょうど 1 本なら `outputs.alignment` に
`{"arf", "arf2"}` を記録する。候補の選び方は `dataset/opening._matching_arf` と共有の `mztab.evidence.first_matching_arfs`:
まずジョブの `peak_matrix_source`(と mzTab-M の隣)、どれも通らなければ解析フォルダのバッチを `AlignmentResult_<時刻>` の
新しい順に 1 本ずつ試し、最初に通ったところで止める(`rerun_upstream=true` の後も前の実行のバッチが残るため)。
合わなければ失敗にせず警告 `ANALYSIS_FOLDER_UNMATCHED`(`source="mztab"` なら読める)で、run は `partial`。

`finish_success` の最終判定(`final_status`): 全単位が validated・受領が成功・全単位の解析フォルダが照合できた → `completed`。
validated が 1 つ以上で欠けがある → `partial`。validated が 0 → `failed`。worker の中の `needs_input` は無い(受付の事前検査
だけが保存する)。

受領(`results/receipt.json`、`lipidmix.pipeline-receipt.v1`)は単位ごとの解析フォルダ・照合した `.arf`・メソッドとライブラリ
(研究室の資産はファイル名と sha256 だけ)・試料の対応表・検証の結果・次の呼び出し(`load_dataset(pipeline_path=…, polarity=…)`)
を載せる。Console の実行体はファイル名と sha256 で示す(版の欄は無い)。解析フォルダが
照合できなかった単位(`alignment: null`)の次の呼び出しには `source="mztab"` を付ける(既定の解析フォルダでは
`ALIGNMENT_NOT_FOUND` / `ALIGNMENT_AMBIGUOUS` で止まる)。
