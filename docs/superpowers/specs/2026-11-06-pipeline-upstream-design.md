# D pipeline を上流に狭める 設計

- 日付: 2026-11-06
- 関連: [全体](2026-11-06-overhaul-overview-design.md)、B [統一データセットとツールの統合](2026-11-06-unified-dataset-and-tools-design.md)、
  C [統計エンジン](2026-11-06-statistics-engine-design.md)。既存の設計
  [2026-09-03-end-to-end-pipeline-design.md](2026-09-03-end-to-end-pipeline-design.md)（このうち下流の部分を置き換える）。
- 状態: 利用者と合意した設計（2026-11-06）。未実装。順番は B → C → A → D の 4 番目。

## 1. 目的

pipeline を「生データから、検証済みの MS-DIAL の解析フォルダを作る」ところまでに狭める。
下流（メタデータの解決・前処理・PCA・差次的解析・書き出し・レポート）は、B・C・A の対話の鎖に任せる。

pipeline にしかない価値は残す: Console の長時間の切り離し実行（実データ 60 試料で約 44 分）、冪等性（`request_id`・指紋）、
中断からの再開、入力（メソッド・MSP・LBM）の sha256 の固定、終了コードだけで完了としない出力の検証。

## 2. 事実（実データの記録）

- reanalysis の pipeline の run は 27 本あり、すべて `pipeline-run.v1`。極性ごとに 1 本ずつ回している
  （例 `datasets/reanalysis/MTBLS3233/pipeline/neg/…`）。
- v2（`pipeline-request.v2`。内部標準の比・feature binding・注入ごとの証拠・固定した QC の母集団・`lcms-profile.v1`）は、
  実データで一度も使われていない。CLAUDE.md の記述でも、合成の入力と偽の Console でしか検証されていない。
- console と pipeline のコードは合わせて約 1.34 万行（パッケージ全体の約 24%）。

## 3. 目標・非目標・反目標

目標:

- pipeline の工程を、入力の準備 → Console の実行 → 出力の検証 → 上流の受領の記録、にする。
- 1 つの run で POS と NEG の両方を回せるようにする（1 極性だけの run も残す）。
- v2 を廃止し、要求の形式を 1 つにする。
- 手で Console を回すツール群を pipeline に吸収する。
- CLAUDE.md の位置づけを「Console は補助」から「対等な 2 つ目の入口」に改める（全体 §2）。

非目標:

- Console の実行・検証・再開の中身を変えること（v8 の `library_mode`・`alignment_reference`・LBM の解決・メソッドの探索は保つ）。
- 下流の工程の新しい作り（B・C・A）。

反目標:

- 既存の 27 本の run が読めなくなる。
- 上流の結果と下流の結果が、pipeline の記録と解析フォルダの `reports/` の 2 か所に分かれて、どちらが正か分からなくなる
  （下流の記録は `reports/` だけに置く）。
- Console の自動再試行を始める（今の規則「やり直しは `rerun_upstream=true` の明示が要る」を保つ）。

## 4. 工程

```
prepare_input → upstream（Console の実行。極性ごと）→ validate_outputs → receipt
```

- `receipt`（上流の受領）: 何をどの入力で回したか（メソッド・ライブラリの sha256・試料の対応・Console の版）と、
  できた解析フォルダの場所・出力の検証の結果を書く。B の `load_dataset(pipeline_path=…)` は、ここが指す解析フォルダを
  解析フォルダのアダプタで読み（B §4.3）、ライブラリの記憶（B §7）と試料の情報もここから取る。
- 出力の検証は今のまま mzTab-M を必須の出力に含める（assay と入力の 1 対 1 の対応の検査に使い、`source="mztab"` での読み込みも
  残すため）。
- 今の `load_dataset`・`resolve_metadata`・`preprocess`・`pca`・`resolve_comparisons`・`differential:*`・`export:*`・
  `report`（解析のレポート）の工程は消す。
- 工程の実行順の正準（`engine.build_stages`）、run の記録の `results` の追記専用、worker がグローバルな session を import しない
  （`tests/test_pipeline_engine.py` の AST の検査）、`pipeline_status` は読取専用、`STATE_REVISION_CONFLICT` の有限回リトライの
  規則は保つ。

## 5. 両極性の run

- 1 つの要求に、極性ごとの単位（生データのフォルダ・メソッド・ライブラリの設定）を並べられるようにする。
  再解析の生データは極性ごとに別のフォルダ（例 `raw/neg`・`raw/pos`）にあることが多い。
- Console は極性ごとに別に実行する（今と同じ）。片方が失敗しても、もう片方の結果は残す。run の状態は極性ごとに読める。
- 要求の形（単位の並べ方、欄の名前）は実装者に任せる。1 極性の要求は今の形のまま受ける。

## 6. 既存の run

- 27 本の v1 の run は、`load_dataset(pipeline_path=…)` で今後もデータセットとして読める（解決済みの試料の対応も読む）。
- run の中の古い差次的結果は引き継がない（今の `session_handoff` の差次的結果を戻す処理は消す）。C のエンジンで計算し直す前提。
- 古い run に `pipeline_resume` で下流の訂正（`comparisons`・`preprocess` の更新）を頼まれたら、下流の工程はもう無いことと、
  代わりのツール（`load_dataset` → `dataset_compare` など）を返して止める。

## 7. v2 の廃止

消すもの（v2 だけが使っているかを実装者が確かめてから消す）:

- `metabolomix/pipeline/request_v2.py`、`stage_plan.py` の v2 の工程計画、`metabolomix/pipeline/metabolomics_handlers.py`
- `lcms-profile.v1`: `metabolomix/console/profile_schema.py`、`profiles.py`、`profile_adapter.py`（`hash_files` のように v1 も使う
  部品は残す場所へ移す）
- v2 の解析の部品: `metabolomix/analysis/` の `assay_evidence.py`・`assay_qc.py`・`feature_bindings.py`・`internal_standards.py`・
  `statistics_v2.py` など（内部標準による正規化は B §6 の定量の工程に一本化する。統計は C）
- 文書: `docs/workflow/metabolomics.md`、CLAUDE.md の「LC–MS メタボロミクス v2」の記述

## 8. ツール

| 残す | 消す（pipeline に吸収） |
|---|---|
| `pipeline_plan`・`pipeline_run`・`pipeline_status`・`pipeline_resume`・`pipeline_cancel`、メソッドを選ぶ `console_method_candidates`・`console_method_template` | `console_plan`、`console_run`、`console_status`、`console_prepare_input`（`keep_extension` で計画の中で行う）、`console_cleanup`（生成物の片付けを pipeline の操作として残すかは実装者の判断）、`job_list` |

`console_*` の内部の部品（`metabolomix/console/` の実行・監視・出力の収集・メソッドの解決）は pipeline が使うので残す。
`analysis-job.json`（受け渡しスキーマ）を読む `load_dataset(job_path=…)` は、既存のジョブを読むために残す。

## 9. 文書とテスト

- `USAGE.md` の §10・§13、`docs/workflow/pipeline.md`、`docs/superpowers/specs/2026-09-03-end-to-end-pipeline-design.md` の
  状態の欄（下流の部分は置き換えた旨）。
- vault の流れ図（pipeline の stage 列が変わるので、同じ作業の中で直す）。
- `MCP_INSTRUCTIONS` の ENTRY POINT の文言（Console の経路が解析フォルダで終わり、その先は同じ鎖を通ること）。
- 両極性の run と、既存の v1 の run の読み込みのテスト（合成の fixture と偽の Console）。

## 10. 実装者に任せること

- 両極性の要求の形と、run の記録の中での極性ごとの状態の持ち方。
- `receipt` の欄の名前。
- 消す順番（ただし旧と新が並んだまま止まらないこと）。
