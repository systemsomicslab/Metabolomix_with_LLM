# パーサー出力フォーマットとオントロジー — 共通核

パーサー出力の行・列・JSON キーの意味の参照資料。**この核だけが全ツール共通**で、個別の定義はトピック別ファイルにある。解釈する出力のトピックだけを追加で引く。

## 0. トピック索引

| トピック | リソース | 収録する節 | 内容 |
|---|---|---|---|
| `core` | `lipidmix://docs/output-format` | §1, §2, §2.1, §9 | 共通オントロジー・脂質名文法（脂質アッセイのみ）・必須注意（**常に先に読む**） |
| `dataset` | `lipidmix://docs/output-format/dataset` | §3, §4, §5, §8, §10, §11 | MS-DIAL 解析フォルダから読んだデータセット。特徴×試料の行・gap-fill と `detected_mask`・`.arf2` の代表と同定の列・試料の役割・前処理と差次的解析 |
| `dcl` | `lipidmix://docs/output-format/dcl` | §6 | `.dcl` パーサ、`dcl_find_msms` ほか `dcl_*` |
| `eic` | `lipidmix://docs/output-format/eic` | §7, §8 | `dataset_plot` の EIC の図（`kind="eic"` / `"eic_compounds"`）。`.EIC.aef` の `peak_top` と描画契約 |
| `identity` | `lipidmix://docs/output-format/identity` | §12 | 同定信頼度・名前正規化・MSI レベル |
| `mztab` | `lipidmix://docs/output-format/mztab` | §13 | mzTab-M から読んだデータセット。SME/SML の同定、特徴表、差次的エクスポートの同定列 |
| `library` | `lipidmix://docs/output-format/library` | §14 | 参照ライブラリ（`.dbs`/`.msp`）照合。スコアの意味、`-1`/`0` の区別、意図的に写した瑕疵 |
| `curation` | `lipidmix://docs/output-format/curation` | — | アラインメントのキュレーション（`curation_review` 等）。判定・理由コードの意味 |
| `pathway` | `lipidmix://docs/output-format/pathway` | — | パスウェイ投影（`pathway_project` / `pathway_plot`、MSLipidMapper 連携）。ノード表・隣接ペア・図の p 値の意味 |
| `reanalysis` | `lipidmix://docs/output-format/reanalysis` | — | 再解析のデータセットのハンドル（`load_dataset` / `dataset_list` / `dataset_close`）と、spec で呼ぶ再解析のツール（行列と統計・MS/MS の合否・証拠の束）。spec と結果の JSON の欄の意味 |
| `quant` | `lipidmix://docs/output-format/quant` | — | 定量（`dataset_quant_preview` / `dataset_quant_apply` / `dataset_merge_polarities` / `dataset_select_features`）。gap-fill の印（`detected_mask`）・単位・`quant_state`・`needs_input` の答え方・記録と audit |

節番号は全トピックの通し番号（`curation`・`pathway`・`reanalysis`・`quant` は節番号を持たない）。`§11.1.1` のような参照はこの表でトピックを引く。

ツール層（`metabolomix/<形式>/tools.py` と `metabolomix/tools/`。`server.py` は登録と再エクスポートだけ）はパーサーの結果を Markdown/JSON 文字列に整形する。表示文ではなく、各トピックが定義する意味で解釈する。

## 1. 目的と対象

| パーサー | 入力 | 主な実装 | トピック |
|---|---|---|---|
| ARF | `*_PeakProperties.arf` | `metabolomix/arf/reader.py` | `dataset` |
| ARF2 | `*.arf2` | `metabolomix/arf2/reader.py` | `dataset` |
| PAI2 | `*.pai2` | `metabolomix/pai2/reader.py` | `dataset`（読むツールは無い。個別測定のピークの定義は `dataset` の PAI2 の節） |
| DCL | `*.dcl` | `metabolomix/dcl/reader.py` | `dcl` |
| EIC/AEF | `*.EIC.aef` | `metabolomix/eic/reader.py` | `eic` |
| mzTab-M | `*.mzTab` | `metabolomix/mztab/reader.py` | `mztab` |

列の定義は `docs/schema/`（`AlignmentSpotProperty.md`・`AlignmentChromPeakFeature.md`・`ChromatogramPeakFeature.md`）。代表テストファイル: `AlignmentResult_2026_05_15_10_13_35{_PeakProperties.arf,.arf2,.EIC.aef}`、`20220901_RAW_control_0h_1_NEG_202605151012{.pai2,.dcl}`。

## 2. 共通オントロジー

> 本書の「オントロジー」はデータ構造の意味体系（行・列の粒度とピークどうしの関係）。フィールド `Ontology` は**脂質クラス**（例 `PC`, `Cer_NS`）で別の概念。名前の記法は §2.1。

同じ LC-MS データを異なる粒度で表す。

```text
データセット
  ├─ アラインメントスポット（全サンプルをまたぐ同一候補ピーク）
  │    ├─ ARF2: スポットの代表値、同定情報、品質統計
  │    ├─ ARF: サンプル別ピークプロパティ
  │    └─ EIC/AEF: サンプル別クロマトグラム点列
  └─ 個別測定ファイル（1サンプル）
       ├─ PAI2: そのサンプルで検出されたピーク
       └─ DCL: PAI2 ピークに対応するデコンボリューション済み MS/MS
```

| エンティティ | 意味 | 主な識別子・対応 |
|---|---|---|
| Alignment spot | 複数サンプル間で同一とみなされたピーク集合 | ARF/ARF2 の `MasterAlignmentID`、`AlignmentID`、EIC の `spot_id`。ID はすべて 0 始まり |
| Aligned sample peak | 1スポットに属する1サンプルのピーク | ARF の1行。`MasterAlignmentID` と `FileID`/`SampleIndex` の組で特定 |
| Raw-file peak feature | 1測定ファイル中の検出ピーク | PAI2 の1要素。`id` は当該 PAI2 内のピークID |
| MSDec result | デコンボリューション済み MS/MS | DCL の1要素。`dcl_index` と PAI2 のリスト順序が対応する設計 |
| EIC sample trace | 1スポット・1サンプルの抽出イオンクロマトグラム | EIC の `spot_id` と `samples[].file_id` の組で特定 |
| Annotation | 候補化合物名 | `Name`/`name`。空文字、`Unknown`、`no MS2:`、`low score:`を含み得るため、存在するだけで確定同定を意味しない |
| Ontology | 化学・脂質クラス | `Ontology`/`ontology`。例: `FA`, `PC`, `PE`, `TG`, `Cer_NS`。化合物名より上位の分類概念 |

| 値 | 意味 |
|---|---|
| RT | retention time、保持時間。通常は分（min） |
| RI | retention index。未使用時は 0 または欠落 |
| m/z | 質量電荷比。単位なし |
| Drift / dt | イオンモビリティのドリフト時間。未使用時は `-1` または欠落 |
| Height | ピーク頂点強度 |
| Area | ピーク面積。`AboveBaseline` はベースラインより上だけを積分した値 |
| S/N | signal-to-noise ratio、信号対雑音比 |
| IonMode | イオン化極性。パーサーにより文字列、Enum、整数のいずれか |

### 2.1 脂質名の記法（ショートハンド文法）

**アッセイ種別が `lipid` のときだけ当てる。** 脂質か一般代謝物かはファイル形式（`.arf` / `.mzTab`）では決まらない。一般代謝物（`assay_kind=metabolite`）では同定を候補集合（アダクト・異性体・候補順位）で読む（`lipidmix://docs/output-format/identity` §12.6）。種別が未確定のあいだはどちらの規則も当てない。

`Name`/`name` と `Ontology`/`ontology` の脂質表記の読み方。構造はこの文法でだけ解釈し、名前から鎖組成やエーテル種別を推定しない（§9-1/§9-2）。

| 表記要素 | 意味 | 例 |
|---|---|---|
| クラストークン | 先頭の脂質クラス＝`Ontology`。headgroup 分類 | `PC`, `PE`, `TG`, `FA`, `Cer_NS`, `SL` |
| `C:D`（sum composition） | 全アシル鎖合計の 炭素数:二重結合数。個々の鎖は未確定＝**species レベル** | `PC 34:1` = 合計C34・二重結合1 |
| `C:D/C:D`（molecular species） | 個々の鎖を確定した表記。合計は species 表記と一致する | `PC 16:0/18:1`（合計 34:1） |
| `/` と `_` | `/`＝sn 位置まで区別、`_`＝鎖は確定だが sn 順不明。MS-DIAL 出力は多くが sn 未確定 | `18:1/16:0` vs `18:1_16:0` |
| `;O`, `;O2`, `;(2OH)` | 追加酸素・水酸基の数。スフィンゴ脂質・酸化脂質で頻出 | `Cer 18:1;O2/16:0`, `SL 33:0;O` |
| `O-` 接頭 | **アルキルエーテル**（plasmanyl） | `PC O-34:1` |
| `P-` 接頭 | **アルケニル（ビニル）エーテル＝プラズマローゲン** | `PC P-34:0` |
| `\|` 区切り | 1スポットに対する**複数候補注釈**の連結。左が代表とは限らない | `SL 33:0;O\|SL 17:0;O/16:0` |
| 確度接頭辞 | `no MS2:` / `low score:` / `w/o MS2:` は注釈の確度限定子（§9-1。§12.1 で除去処理） | `low score: PC 34:1` |

注意点:

- **species と molecular species は別粒度**。`PC 34:1` は鎖組成を主張しない（`16:0/18:1` も `17:0/17:1` も sum 34:1）。
- **`P-` と `O-` は species レベルで曖昧**。pygoslin は `PC P-34:0` を `PC O-34:1` に同一化する。区別は正規化名ではなく、元の `Name` の表記と `.dcl` の MS/MS で読む。この注意を添えるツールは無い（§12.1）。
- **`Ontology` の空文字と文字列 `Unknown` は別扱い**（§9-2）。クラストークンが取れないものをクラス化しない。
- 同じスポットで ARF の `Name` と ARF2 の `Name` が食い違うことがある（§11.1.1）。統合注釈は ARF2 側を優先し、不一致は所見として書く。

## 9. LLM解釈時の必須注意事項

1. `Name` が存在しても確定同定とは限らない。空文字、`Unknown`、`no MS2:`、`low score:`を区別する。
2. `Ontology` は化合物名ではなく分類クラスである。空文字と文字列 `Unknown` も区別する。
3. ARFの1行はサンプル別ピーク、ARF2の1行は全サンプル統合スポットであり、同じ「1行」でも粒度が違う。
4. ARFの `IsGapFilled=true` は実測ピークではなく補間値である。検出率や存在判定では別扱いする。
5. ARFの `HeightAverage` は現実装ではグループ先頭値であり、名前どおりの再計算平均ではない。統合統計にはARF2側を優先する。
6. ARF要約の `named_compounds` は空文字も数えるため、注釈率には使わない。
7. PAI2 は単一サンプルのピーク一覧で、オミクス PCA はできない（PAI2 を読むツールは無い）。複数サンプルの多変量比較は `dataset_pca`。
8. PAI2の `has_msms=true` は取得参照があることだけを示し、PAI2内 `msms_spectrum` が非空とは限らない。実スペクトルはDCL（`dcl_spot_spectrum` / `dcl_diagnostic_ions` / `dcl_find_msms`）。**MSI Level 2 を MS/MS で主張するときは `.dcl` の実スペクトルを見てからにし、取得のフラグ（`has_msms` など）だけでは主張しない。**
9. DCLの `n_msms_peaks` は元本数であり、`top_n_peaks` 適用後の配列長とは異なり得る。
10. EICの `peak_top` は強度ではなく頂点座標である。強度は `max_intensity`、平均強度は `mean_intensity` を使う。
11. `total_samples` はEIC全スポットにわたるサンプルエントリ総数であり、ユニークサンプル数は `len(unique_file_ids)` である。
12. m/zとRTの微小差はファイル形式の浮動小数精度、代表値の定義、アラインメント処理に由来し得る。厳密一致ではなく許容差を用いる。
13. **mzTab-M（`dataset_*` 経路）の `abundance_assay[N]` は、非ゼロでも実測ピークとは限らない。** mzTab-M は gap-fill（未検出セルの補間値）を区別する列を持たず、`load_dataset` は隣接する `.arf` が同一アライメントだと数値で検証できたときだけ検出状態を取り込む。`dataset_status` の `detection.available` が `false` の間は、**検出率・欠測率・「n 件で検出」を主張しない**（実データ 60 サンプル × 714 特徴でセルの 70.0% が gap-fill だった）。検出状態があれば `dataset_preprocess(detection_min_rate=...)` で gap-fill を数えない検出率で足切りできる。無ければこの足切りは行わず、前処理の `steps` の `step: "filter"` に `filters.detection = {status: "skipped", reason: "no_detected_mask"}` と出る。

## 状態不足の伝え方（MCP クライアント向けの契約）

多くのツールは先行ツールが作ったセッション状態に依存する（`dataset_plot(kind="pca")` は `dataset_pca` の結果、`kind="volcano"` は `dataset_compare` の結果）。`dataset_compare`・`dataset_pca`・`pathway_project` は前処理済みの行列（`dataset_preprocess` の `matrix_id`）を受け、省くとそのハンドルで `dataset_preprocess` が最後に作った行列を、それも無ければ既定の前処理で自分で作った行列を使う（戻り値の `matrix_source`）。

状態が無いとき、ツールは本文に次の JSON を返す。

```json
{"error": {"code":           "missing_state",
           "state":          "stat_matrix",
           "required_tools": ["dataset_preprocess"],
           "message":        "行列 'mx_…' がありません（無い・または除外や試料の表の変更で古くなりました）。dataset_preprocess を実行し直して、その matrix_id を渡してください。"}}
```

- `state` は欠けている状態の識別子。**不透明な文字列として扱う**。
- `required_tools` は状態を作れるツールの候補（OR）。サーバは実際の生成元をすべて挙げる。クライアントは候補のうち「リプレイ安全かつ本セッションで成功した呼び出し」の中から**直近に成功したもの**を選んで再実行する（リスト順では選ばない）。
- `message` は LLM と人が読む説明。エンベロープを解釈しないクライアントにはこれだけが見える。
- `details` は任意（無いことが多い）。今は `dataset_plot` だけが、その kind で使わない引数を外したときに `details.notes`（どの kind の引数か）を載せる。復旧の判断には使わない。

エンベロープを解釈するクライアントは、`required_tools` のいずれかを**セッション履歴に記録済みの引数のまま**再実行して状態を復元できる。LLM に再実行させると前処理の引数が変わりうるので、引数の同一性はクライアント側で担保する。

### 利用者に尋ねる封筒（`needs_input`）

`missing_state` と並ぶ機械可読な封筒。**エラーではない**（`error` キーを持たない）。サーバが既定を決めずに `options` からの選択を求める問いで、答えを推測して呼び直さない。

```json
{"needs_input": {"field":       "polarity",
                 "reason":      "この run には極性の単位が 2 つある。どちらを読むか…",
                 "options":     [{"value": "positive", "label": "…"}, {"value": "negative", "label": "…"}],
                 "argument":    "polarity",
                 "retry_tools": ["load_dataset"]}}
```

- `field` は答えが入る設定の欄（spec の路）。`reason` は尋ねる理由、`options` は `{value, label}` の選択肢。
- `argument` は答え（選ばれた `value`）を入れて呼び直すツールの**引数名**。`argument` が無い封筒は、`field` の spec の路に `value` を入れて呼び直す。
- `"patch": true` の封筒は、選択肢ごとに替える引数が違う問い（`argument` を持たない。`field` は問いの主題）。選ばれた `value` は「引数名 → 値」の dict で、**その dict のとおりに引数を置き換え（値が `null` の引数は外し）、ほかの引数はそのまま渡して**呼び直す。例: `dataset_compare` の `{"comparisons": [...], "reference_group": null, "test_group": null}`（群の指定を対比に替える）・`{"factors": null}`・`{"matrix_id": null}`。
- `retry_tools` は答えを入れて呼び直せるツールの候補（OR）。`argument` があればその引数を、`patch` なら `value` の引数をすべて受ける。
- 統計エンジンの前処理（`dataset_preprocess`・`dataset_compare`・`dataset_pca`・`pathway_project`・`dataset_analysis_run` の spec v2）は**尋ねない**: ドリフト補正の既定（`drift_mode=auto`）は、補正できても内部標準で換算済み・注入順がファイルの順かもしれないときは当てずに計算し、戻り値の `warnings`（`pathway_project` は run の manifest にも）の `drift_not_applied_by_default`（`reasons` と `suggestion: {"preprocess.drift.mode": "on"}`）で呼び直し方を示す（`dataset.md` §10.2。`dataset_pca` / `pathway_project` は補正を当てなかった行列を使ったときだけ載せる）。
- 定量の `dataset_quant_preview` / `dataset_quant_apply` / `dataset_merge_polarities` / `dataset_select_features` は、上の封筒ではなく各ツールの戻り値の形で `status="needs_input"` を返す（`quant.md`）。

### ハンドルが塞がっている（`HANDLE_BUSY`）

同じハンドルへの呼び出しは 1 本ずつ動く。前の呼び出し（多くは `pathway_project` / `dataset_homolog_check` が受付番号を返して続けている仕事）がハンドルを持ったまま 25 秒たつと、後の呼び出しは**何もせずに** `{"error": {"code": "HANDLE_BUSY", "details": {"dataset_handle", "waited_s", "job_id", "tool", "running_s", "next"}, "required_tools": ["job_status"], …}}` を返す（`job_id` 以下は長い仕事が持っているときだけ）。`job_status(job_id)` で終わりを待ってから、**同じ呼び出しをし直す**。待ちきれなかった呼び出しが後から効くことはない。`dataset_close` も仕事が走っているハンドルは閉じない（`reanalysis.md` の「並行と直列」）。2 つのハンドルを取る `dataset_merge_polarities` は 2 本**合わせて**最長 25 秒待つ（`details.dataset_handle` は空かなかったほう）。

### 判定 3 値と「有意でない ≠ 変わっていない」

統計エンジンの比較（`dataset_compare`・`dataset_analysis_run` の spec v2）は行ごとに `verdict` を返す: `changed`（検定でき `q ≤ 閾値`）/ `unchanged`（`q` が閾値を超え、95% CI が同等の範囲 `±log2(equivalence_fc)` に収まる）/ `inconclusive`（それ以外すべて）。**「有意でない」は「変わっていない」ではない**: n が小さい・ばらつきが大きいときの大半は `inconclusive`（データが足りない）で、「差なし」と報告しない。規則・`verdict_reason`・列の意味は `reanalysis.md` §2.1 と `dataset.md` §11.2。

`log2fc_threshold` は判定に使わない（図の色分け・`large_changes`・書き出しの `significant` 列だけが使う）。書き出し・根拠バンドルの `significant` と volcano の `up` / `down` は、`q ≤ 閾値` かつ `|log2fc| ≥ log2fc_threshold` かつ `verdict` が `changed` の行だけ（`large_changes` と同じ数え方）。`inconclusive` の行は q が小さくても有意にしない（`dataset.md` §11.4）。書き出しと根拠バンドル（別リポジトリとの契約）の `q_value` は比較の `p_adjust` によらず BH の q で、判定と `significant` もその q で読む（BH でない比較は同じ母集団の p から計算し直し、メタ行 `q_value_source` に書く）。

### annotations の解釈

全ツールが MCP 標準の `ToolAnnotations` を宣言する。

- `readOnlyHint=true` は**サーバの外（ファイルとネットワーク）に副作用が無い**こと。サーバ自身の解析セッション状態の更新は副作用に数えない（その依存は上のエンベロープで伝える）。
- `openWorldHint=true` は外部ネットワークへ出る（`server_update` と `pathway_*`）。
- リプレイ安全（同じ引数での再実行が安全）は `readOnlyHint=true` または `idempotentHint=true`。
