# pipeline-run.v2 — 検証済みLC–MSメタボロミクスの実行記録

検証済みLC–MSメタボロミクス実行経路（spec
[2026-09-15-validated-lcms-metabolomics-design.md](../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md)
§6, §6.2）が`pipeline_root/pipeline-run.json`へ書く永続記録。実装は
[`metabolomix/pipeline/store.py`](../../metabolomix/pipeline/store.py)（永続化・schema判定）と
[`metabolomix/pipeline/stage_plan.py`](../../metabolomix/pipeline/stage_plan.py)（stage列・依存無効化表）。
要求契約は[`pipeline-request.v2`](pipeline-request-v2.md)。

既存の`pipeline-run.v1`（リピドミクスv1経路。実装docstringが正準）とは**別の記録**で、
stage IDの綴りもstage集合も異なる。v1のrunはv1のまま新規作成・再開し、
**過去の記録をインプレース変換しない**。

## schema dispatch（readerとwriterの両方）

| 場面 | 規則 |
|---|---|
| 新規作成 | `store.run_schema_for(request)`が要求のschemaで決める。`pipeline-request.v2`なら`pipeline-run.v2`、`pipeline-request.v1`または**schema省略**なら`pipeline-run.v1`。それ以外は`PIPELINE_REQUEST_INVALID`で拒否（推測しない） |
| 読込 | `store.load_run`は`store.SUPPORTED_SCHEMAS`（v1とv2）を受け、**読んだrecordのschemaのまま返す**。未知のschemaは`PIPELINE_RUN_INVALID` |
| 保存 | `store.save_run`は既存recordのschemaと異なるschemaでの上書きを`PIPELINE_RUN_SCHEMA_IMMUTABLE`で拒否する。v1のrunをv2として保存し直す経路は存在しない |

拒否の判定は`create_run`がpipelineディレクトリを掘る**前**に行う——拒否された要求が
空の`pipeline_root`を残さないため。

`state_revision`（楽観的並行制御のカウンタ）・`results`の追記専用性・成果物hashの
照合（`verify_result_refs`）は、v1とまったく同じ仕組みをそのまま使う。

## トップレベル

v1と同じ形。`schema`の値と`stages`のキー集合だけがv2固有になる。

| フィールド | 内容 |
|---|---|
| `schema` | `"pipeline-run.v2"` 固定 |
| `identity` | `pipeline_id`（UUID4）・`source_root`・`pipeline_root`・`created_at`・`updated_at` |
| `request` | `revision` / `request_id` / `content_hash` / `saved_path` / `effective_target` の要約。要求本体は`requests/revision-NNNN.json`（追記のみ） |
| `status` | `planned` / `running` / `needs_input` / `completed` / `partial` / `failed` / `cancelled` |
| `state_revision` | 保存のたびに1増える並行制御カウンタ（`request.revision`とは別物） |
| `stages` | 下記「stage列」のstage IDをキーにしたdict |
| `upstream` | `console_job_path` / `execution_id` / `verification` |
| `inputs` | 固定した入力のスナップショット |
| `results` | 追記専用の結果履歴 |
| `needs_input` | 中断理由（`code` / `stage_id` / `message` / `details`） |
| `warnings` | 受付時から積む警告 |
| `worker` | owner lockを取ったworkerの`identity`と`started_at` |

`request.content_hash`は`metabolomix.pipeline.request.request_fingerprint`が計算する。
**v2要求はv2のキー集合で畳む**——v1のキー集合で畳むと`statistics`等が落ち、統計定義
だけが違う2つの解析が同じhashになって1本のrunへ畳まれてしまう（spec §6.1
「下流結果IDは…変換、群・検定設定に依存する」）。

## stage列

spec §6.2の順序そのまま。正準は`stage_plan.build_v2(request)`で、`store`（`stages`を
作る側）も`engine`（実行順を決める側）も**この1つのbuilderを呼ぶ**。実行順の正準は
`engine.build_stages`の計画順であって`stages`dictの挿入順ではなく、`report`は必ず最後。

| stage ID | handler | statistic_id | 役割 |
|---|---|---|---|
| `prepare_inputs` | 同名 | `null` | 生データ・sample manifest・profile依存の受付と固定 |
| `execute_console` | 同名 | `null` | MS-DIAL Consoleの実行（**自動再試行しない**） |
| `validate_outputs` | 同名 | `null` | 終了証跡・必須生成物・入力不変の検証 |
| `load_dataset` | 同名 | `null` | mzTab-Mの読込とDatasetStateの構築 |
| `resolve_metadata` | 同名 | `null` | sample manifestとassayの対応付け |
| `load_assay_evidence` | 同名 | `null` | 注入ごとの測定証拠の取込（spec §9.1） |
| `resolve_feature_bindings` | 同名 | `null` | profileの論理target_id → バッチ固有feature_id（未解決は`needs_input`） |
| `qc_raw` | 同名 | `null` | 元データQCの評価集合の固定と評価 |
| `preprocess` | 同名 | `null` | 検出filter・内部標準比/正規化・ドリフト補正・補完 |
| `qc_processed` | 同名 | `null` | 処理後QCの評価と最終matrixの確定 |
| `statistics:<statistic_id>` | `statistics` | `<statistic_id>` | 統計1件（`pca` / `welch` / `anova_tukey`）。要求の`statistics`順に展開 |
| `export:<statistic_id>` | `export` | `<statistic_id>` | 同じ統計の図・TSV出力 |
| `report` | 同名 | `null` | 必須出力判定と品質レポート（常に最後） |

`statistic_id`は`pipeline-request.v2`が保証する一意な安全slug。`stage_plan.build_v2`は
空の`statistics`と重複`statistic_id`を`PIPELINE_REQUEST_INVALID`で拒否する
（重複は`stages`のdictキーとして潰れ、stageが1つ黙って消えるため。`request_v2`側の
一意性検査に加えた二重の防御）。

### 上流を自動で再実行しない

`prepare_inputs` / `execute_console` / `validate_outputs`は`engine._TRUST_PERSISTED_STAGE_IDS`
に属する——永続状態が`succeeded`/`skipped`ならそのまま信頼し、engineが自分の判断で
MS-DIALを起動し直すことはない。やり直しには`pipeline_resume(rerun_upstream=true)`の
明示が要り、そのとき新しいattemptディレクトリを作る（前回の終了証跡・ログを
上書きしない）。再開時の「上流をやり直すか」の判定が見るstageは、v1の`upstream`では
なく`execute_console`（`recovery._upstream_stage_id`）。

## 依存無効化

spec §6.2の依存区分。表の正準は`stage_plan.CHANGE_ENTRY_POINTS`と
`stage_plan.invalidated_v2(changed, request)`で、`recovery.prepare_resume`は要求の差分を
下記のtokenへ畳んで渡すだけ（区分そのものを持たない）。

| 変更token | 最初に無効化するstage | 出所 |
|---|---|---|
| `metadata` | `resolve_metadata` | 同じパスのsample manifestを書き直した訂正（要求の値は変わらない） |
| `sample_manifest` | `resolve_metadata` | 要求の`sample_manifest`の変更 |
| `standard_assays` | `resolve_feature_bindings` | 要求の`standard_assays`の変更 |
| `feature_bindings` | `resolve_feature_bindings` | 要求の`feature_bindings`の変更（resumeの候補選択） |
| `preprocess` | `preprocess` | matrix recipeの上書きの変更 |
| `qc_raw_scope` | `qc_raw` | 元データQCの**評価集合**を動かすrecipe変更。要求の差分だけでは判定できないため、集合を組むstage handlerだけがこのtokenを立てる |
| `statistics` | 全`statistics:*` / `export:*` と`report` | 統計定義の全面変更 |
| `target` | 全`statistics:*` / `export:*` と`report` | `exploratory` ↔ `differential`（統計集合の整合そのものが変わる） |
| `statistics:<statistic_id>` | その`statistics:<id>` / `export:<id>` と`report` | 統計1件の追加・変更・**削除** |

入口stageはどれも`load_dataset`より下流にある。したがってbinding・metadata・recipeの
訂正が`prepare_inputs` / `execute_console` / `validate_outputs`を無効化することはない
——**binding訂正で上流Consoleを再実行しない**（spec §6.2）。

削除された統計のstage IDは新しい計画には現れないが`stages`には残るため、
`statistics:<削除されたid>`として無効化し、`report`も作り直す（spec §6.2
「削除済み統計の成果物を必須出力に残さない」）。

未知のtokenは`PIPELINE_STAGE_PLAN_INVALID`で拒否する。黙って無視すると「無効化した
つもりで何もしていない」——古い結果をcurrentのまま残す、最も見つけにくい壊れ方になる。

## 関連

* 要求契約: [pipeline-request.v2](pipeline-request-v2.md)
* プロファイル契約: [lcms-profile.v1](lcms-profile-v1.md)
* Console成果物の受け渡し: [analysis-job.v3](analysis-job-v3.md)
