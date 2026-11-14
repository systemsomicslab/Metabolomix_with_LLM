# 定量の出力（`dataset_quant_preview` / `dataset_quant_apply` / `dataset_merge_polarities` / `dataset_select_features`）

内部標準による正規化・陽陰の統合・特徴の選別。規則は外部の規則プラグイン（`[quant] rules`。
仕様は `docs/quant-rules-plugin.md`）が持ち、未設定なら 4 ツールとも
`{"error": {"code": "quant_rules_unavailable", ...}}`（message に設定の書き方）を返す。
共通核（`lipidmix://docs/output-format`）の必須注意（gap-fill と検出の区別ほか）が先に効く。

**入力は MS-DIAL の解析フォルダ（`.arf` / `.arf2`）から読んだデータセットだけ**
（`source_format == "msdial_folder"`）。規則の入力（クラス・付加体・コメント・MS/MS の有無）を
同じ意味で持たないほかの形式（mzTab-M など）は、空欄のまま黙って定量せず
`QUANT_SOURCE_UNSUPPORTED` で止める。

## 1. 単位と値

- `unit` は規則が返す文字列（例: 細胞あたりの pmol）。**正規化していないデータセット（`skip=True` で
  飛ばしたもの、正規化の工程を通っていないもの）は `normalized=false`・`unit=null`**で、値は MS-DIAL の
  強度のまま。単位が違うデータセット同士・正規化済みと未正規化は統合しない（`QUANT_MERGE_STATE_MISMATCH`）。
  未正規化どうしは、skip を挟まなくてもそのまま統合する。
- **正規化は統合の前に、極性ごとに行う**。内部標準は極性ごとの行で、参照表はクラスで引くので、陽陰を統合した
  データセットを割ると一方の極性の内部標準で他方の極性の同じクラスを割りうる。統合したデータセット（`quant_state.steps`
  に `merge`、または特徴に `origin`）には preview・apply・`skip=True` のどれも掛けず `QUANT_ALREADY_MERGED`
  （`details.parent_handles` に極性ごとの親のハンドル。そのハンドルで正規化してから統合し直す）。極性を尋ねる前に断る。
- 正規化した行列の行は**内部標準以外の特徴**。内部標準の行は値を持たず、後段の選別で除かれる。
- 内部標準の値が 0・負・欠測の試料は、その試料の値を**欠測（NaN）**にして理由を `status_counts`
  に数える。中央値などで埋めない。参照表に無いクラスは `Any others` の行で計算する
  （`classes[].via_any_others`）。

## 2. gap-fill の印（最重要）

gap-fill（検出されず MS-DIAL が補った値）のセルも**値は計算する**が、`detected_mask` は
**対象の特徴か、割った内部標準のどちらかが gap-fill なら False**。つまり
`detected_mask=True` は「対象も内部標準も実際に検出された」セルだけ。値が入っていることを検出の
証拠として読まない。統計側の検出率・検出例数はこの印を使う。

## 3. `quant_state`（派生データセットが持つ状態）

統計（`dataset_*` の前処理・検定）は `getattr(ds, "quant_state", None)` の `normalized` / `method` /
`unit` と、統合したデータセットの `polarity_runs`（陰イオンの測定。ドリフト補正を極性ごとに当てる）を読む。
それ以外のキーは自由。

| キー | 意味 |
|---|---|
| `schema` | `lipidmix.quant-state.v1` |
| `normalized` | 内部標準で正規化済みか。skip した派生は false |
| `method` | 正規化なら `internal_standard`、skip は null |
| `unit` | 値の単位（null = 強度のまま） |
| `record_id` | この工程の記録 ID（`q_<UTC 時刻>_<8 桁>`） |
| `steps` | 通った工程（`normalize` / `merge` / `select`） |
| `skipped` | 飛ばした工程（`normalize` / `adduct_filter`） |
| `internal_standard_ids` | 内部標準の特徴 ID（選別後は残ったものだけ。統合後は `positive:<id>` / `negative:<id>`） |
| `record_folder` | 統合後の記録の置き場（`dataset_merge_polarities` が足す） |
| `arf2_by_polarity` | 統合後に、極性ごとのキュレーションのフラグを引く `.arf2`（`dataset_merge_polarities` が足す。フラグそのものは写さない） |
| `polarity_runs` | 統合後の、陽イオン以外の極性の測定（`{"negative": {<統合後の試料名>: {batch, injection_order, injection_order_source, injection_order_may_be_file_order}}}`。陰イオンのデータセットから統計エンジンと同じ規則で読んだ値。`dataset_merge_polarities` が足す）。統合したデータセットの試料の表は陽イオンの測定 |
| `polarity_runs_match` | `{"negative": <negative_runs と同じもの>}`（陰イオンの測定が陽イオンと同じか） |

## 4. `status` と `needs_input`

| ツール | `status` | 意味 |
|---|---|---|
| preview | `ready` / `needs_input` | `needs_input` は内部標準の候補が複数（`standard_ambiguous`）か無い（`standard_missing`）クラスがある |
| apply | `applied` / `skipped` / `needs_input` | `needs_input` のときは**何も書かず**、preview と同じ封筒 |
| merge | `merged` / `needs_input` | 試料の食い違い（`reason=sample_mismatch`）か記録の置き場（`record_folder_required`） |
| select | `selected` / `needs_input` | 採用表が未設定（`reason=adduct_table_required`） |

`needs_input` の封筒は `how_to_answer` に次の呼び方を書く。答え方:

- **内部標準の曖昧**: `unresolved[].candidates`（`feature_id`・名前・付加体・m/z・RT）から 1 つ選び、
  `bindings={"<standard_name>": <feature_id>}` を足して同じツールを再実行。
- **試料の食い違い**: `only_in` は片側にしか無い試料。`mapping_csv`（列: `mode`, `source_sample`,
  `sample_id`, `group`, `sample_type`。`mode` は positive / negative、pos / neg も可）で対応づけて再実行。
  `group` / `sample_type` 列は規則が検証するだけで適用しない（試料の実験情報は POS のもの）。
- **記録の置き場**: 2 つの解析フォルダの共通の親が決められない（ドライブが違う・親子）とき
  `record_folder`（既存のフォルダ）を渡す。
- **採用表**: `lipidmix.local.toml` の `[quant] adducts` を設定するか `skip_adduct_filter=True`
  （内部標準・名前・MS/MS の除外は行う。`quant_state.skipped` に `adduct_filter` が入る）。

エラー封筒は他の `dataset_*` ツールと同じ形（`{"error": {"code", "message", "details"?}}`）。主な code:
`quant_rules_unavailable`（規則未設定）・`QUANT_RULES_REJECTED`（規則が入力を拒否。`details.rule_code`）・
`QUANT_BAD_REQUEST`（UTF-8 でない表もこれ）・`QUANT_POLARITY_REQUIRED`（`polarity` を渡す）・
`QUANT_SOURCE_UNSUPPORTED`（解析フォルダ以外の入力）・`QUANT_ALREADY_NORMALIZED`・
`QUANT_ALREADY_SELECTED`（選別済みのデータセットは正規化しない）・
`QUANT_ALREADY_MERGED`（陽陰を統合したデータセットは正規化しない。§1。`details` に `record_id` と `parent_handles`）・
`QUANT_CURATION_UNREADABLE`（キュレーションの記録が読めない。`apply_curation=False` で進められる）・
`DATASET_HANDLE_LIMIT`（新しいハンドルを足せない。記録を書く前に止める）・
`QUANT_NEEDS_DETECTION`（検出の印が無い入力）・`QUANT_MERGE_STATE_MISMATCH`（正規化の状態が陽陰で違う。
`details.differs` と `details.normalize`）・
`QUANT_MERGE_SAMPLE_CONFLICT`（同じ試料の群・役割が陽陰で食い違う）・`QUANT_HANDLE_BUSY`・
`QUANT_REFERENCE_NOT_FOUND`。封筒・戻り値には**ローカルのパスを載せない**（参照表は sha256・
単位・ファイル名の末尾だけ）。

## 5. 各ツールの戻り値の欄

- **preview**: `status`・`polarity`・`unit`・`reference`（`sha256` `label` `method`）・
  `classes`（データセットに現れたクラスだけ: `target_class` `standard_name` `feature_id` `adduct`
  `mz` `rt` `coefficient` `via_any_others`）・`unresolved`・`missing_values`（試料の値が欠測になる理由ごとの件数）。
- **apply**: `status`・`dataset_handle`（新しいハンドル）・`record_id`・`status_counts`（セルの状態ごと
  の件数）・`unit`・`polarity`・`reference`・`classes`・`internal_standard_ids`。
- **merge**: `status`・`dataset_handle`・`record_id`・`rule`（試料の照合規則）・`n_features`
  （`positive` `negative` `merged`）・`n_samples`・`unit`・`normalized`・`normalize`（極性ごとに正規化を
  どう通ったか: `applied`＝dataset_quant_apply / `skipped`＝`skip=True` / `not_run`＝どちらも通っていない。
  記録の `settings.normalize` にも同じものを書く）・`negative_runs`（陰イオンの測定が陽イオンと同じか:
  `injection_order` は両方分かる試料のどれかで違えば `different`、違わないが片方でも分からない試料があれば
  `unknown`、それ以外は `same`。`batch` は `same` / `different`。`n_samples_differ` は注入順かバッチの違う試料の数。
  記録の `summary.negative_runs` にも書く。陽陰は別の注入なので、統計の前処理は極性ごとの測定でドリフトを当てる）。
  特徴 ID は `positive:<id>` / `negative:<id>`、行は POS の特徴のあとに NEG、試料の列は共通 ID
  （POS の順）。各特徴の出自は `feature_metadata[<fid>]["origin"]`（極性と元の特徴 ID）。
- **select**: `status`・`dataset_handle`・`record_id`・`audit_file`・`n_retained`・`n_excluded`・
  `reasons`（理由ごとの件数）・`curation`・`curation_by_polarity`・`adducts`（`sha256` `label`）・
  `skipped`・`examples`（除外の例 5 件）。**全判断は audit にあり、戻り値には載せない。**
- **探索専用の印**: 派生は親の `exploratory_only` / `exploratory_reason` を引き継ぐ（未完了の実行・検証を
  通っていない pipeline の単位から開いたデータセットの派生も、書き出しのメタ行・根拠バンドル・図に印が載る）。
  統合はどちらかの親が探索専用なら探索専用で、理由は極性付きで `; ` でつなぐ（`negative: <理由>`。理由の無い
  印は `incomplete analysis run (exploratory only)`）。派生が探索専用のときだけ、apply・merge・select の戻り値に
  `exploratory`（`{exploratory_only: true, reason}`）が載る（でなければ欄が無い）。
- **出所の裏取り**: 統合したデータセットの `source_verification` は、陽陰で同じならその値、違えば弱いほう
  （`verified` > `legacy_unverified` > `direct_unverified`。知らない値はいちばん弱いもの）。書き出しのメタ行
  `# source_verification` と `dataset_status` もこの値。陽陰で違うときだけ、merge の戻り値に
  `source_verification`（`merged` と `positive` `negative`）が載る（同じなら欄が無い）。

## 6. 記録・audit・古さ

- 1 工程 1 記録: 解析フォルダ（統合後は 2 つの解析フォルダの共通の親、または `record_folder`）の
  `reports/quant/<record_id>.json`（schema `lipidmix.quant-record.v1`）。設定の `[paths] reports_dir` は見ない
  （再開がデータと一緒に読み直す記録なので、いつも解析フォルダ側）。入力の指紋（元ファイルの
  sha256・参照表 `reference:<名前>`・採用表 `adducts:<名前>`）・前の工程の記録 ID（`parent_ids`）・
  設定（`settings.rules` に規則の版と digest）・結果の要約・時刻。**上書きも削除もしない**。
  やり直しは新しい記録。リポジトリ側には書かない。
- 選別のキュレーションは、読み込み時の写しではなく**選別の時点の** `flags.jsonl` を読み直して当てる
  （統合後は極性ごとの `.arf2`）。当てた判断（wrong / redundant / assign の name・ontology）の
  digest を極性ごとに `settings.curation_digest` に残す（`apply_curation=False` なら null）。
- 選別の判断の全件は `<record_id>.audit.json`（`lipidmix.quant-select-audit.v1`: `decisions` と
  `curation`）。記録としては読まれない。
- 入力が変わる（元ファイルや参照表の sha256 が違う）か前の記録が無くなる・古くなると、その記録は
  **古い**（`input_changed:` / `parent_stale:` / `parent_missing:`）。選別の記録はキュレーションの
  digest がいまと違っても古い（`curation_changed`）。消さずに示すだけで、規則の digest が
  記録と違えば計算し直した派生は戻さない。途中の工程を戻せなければ、それより後ろの記録は
  `prior_step_not_restored` として挙げる。

## 必須の注意

1. gap-fill のセルは値が入っていても検出ではない（§2）。`detected_mask` を見る。
2. 単位は規則が決める。`normalized=false` の値に単位は無い。単位・normalized が違うものを混ぜない。
   統合したデータセットの `normalize` が `not_run` の側は、正規化の工程を一度も通っていない強度のまま。
3. 内部標準の欠測・0・負の試料の値は欠測。補完していない。
4. 統合後の特徴 ID は `positive:` / `negative:` の接頭辞つき。統合前の ID・キュレーションの
   フラグと同じものではない（フラグは極性ごとの `.arf2` から元の特徴 ID で引く）。
5. 参照表・採用表・内部標準の係数の**中身はこのサーバの戻り値に載らない**。データセットに現れた
   クラスの分だけ `classes` に出る。
6. 選別の除外は audit を見て判断する。付加体の選別を飛ばした派生は `skipped` に `adduct_filter` が入る。
