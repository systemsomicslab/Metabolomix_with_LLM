# パーサー出力フォーマットとオントロジー — 共通核

この文書は、LLM が本リポジトリのパーサー出力を解析するときに、各行・列・JSON キーが何を表すかを誤解しないための参照資料である。

**核（この core）だけが全ツール共通**で、パーサー/ツール個別の定義はトピック別ファイルに分かれている。解釈しようとしている出力のトピックだけを追加で引くこと。

## 0. トピック索引

| トピック | リソース | 収録する節 | 内容 |
|---|---|---|---|
| `core` | `lipidmix://docs/output-format` | §1, §2, §2.1, §9 | 共通オントロジー・脂質名文法（脂質アッセイのみ）・必須注意（**常に先に読む**） |
| `arf` | `lipidmix://docs/output-format/arf` | §3, §8.1, §10, §11 | `.arf` パーサ、`arf_parser`、前処理・QC、差次的解析 |
| `arf2` | `lipidmix://docs/output-format/arf2` | §4, §8.2 | `.arf2` パーサ、`arf2_parser` |
| `pai2` | `lipidmix://docs/output-format/pai2` | §5, §8.3 | `.pai2` パーサ、`pai2_parser` |
| `dcl` | `lipidmix://docs/output-format/dcl` | §6 | `.dcl` パーサ、`dcl_parser` / `dcl_find_msms`、PAI2 への MS/MS 付与 |
| `eic` | `lipidmix://docs/output-format/eic` | §7, §8.4, §8.5 | `.EIC.aef` パーサ、EIC 検索、描画契約 |
| `identity` | `lipidmix://docs/output-format/identity` | §12 | 同定信頼度・名前正規化・MSI レベル |
| `mztab` | `lipidmix://docs/output-format/mztab` | §13 | mzTab-M 経路（`dataset_load` 以降）。SME/SML の同定、特徴表、差次的エクスポートの同定列 |
| `library` | `lipidmix://docs/output-format/library` | §14 | 参照ライブラリ（`.dbs`/`.msp`）照合。スコアの意味、`-1`/`0` の区別、意図的に写した瑕疵 |

節番号は分割前の通し番号をそのまま保持している。本文中の `§11.1.1` のような相互参照は、この表からトピックを引いて辿ること。

MCPツール層（形式ごとの `lipidmix/<形式>/tools.py` と `lipidmix/tools/`。`server.py` はそれらを登録・再エクスポートする薄いファサード）が、パーサーの構造化結果を主に Markdown/JSON 文字列へ整形する。LLM は表示文ではなく、各トピックが定義する意味を基準に解釈する。

## 1. 目的と対象

README に記載された主要な MS-DIAL 出力パーサーを対象とする。

| パーサー | 入力 | 主な実装 | トピック |
|---|---|---|---|
| ARF | `*_PeakProperties.arf` | `lipidmix/arf/reader.py` | `arf` |
| ARF2 | `*.arf2` | `lipidmix/arf2/reader.py` | `arf2` |
| PAI2 | `*.pai2` | `lipidmix/pai2/reader.py` | `pai2` |
| DCL | `*.dcl` | `lipidmix/dcl/reader.py` | `dcl` |
| EIC/AEF | `*.EIC.aef` | `lipidmix/eic/reader.py` | `eic` |
| mzTab-M | `*.mzTab` | `lipidmix/mztab/reader.py` | `mztab` |

記載内容は、実装、`docs/schema/AlignmentSpotProperty.md`、`docs/schema/AlignmentChromPeakFeature.md`、`docs/schema/ChromatogramPeakFeature.md`、および NEG 実測データセット（`LIPIDMIX_DATA_DIR` 配下）の実ファイルに対する出力確認に基づく。

代表テストファイル: `AlignmentResult_2026_05_15_10_13_35{_PeakProperties.arf,.arf2,.EIC.aef}` および `20220901_RAW_control_0h_1_NEG_202605151012{.pai2,.dcl}`。

## 2. 共通オントロジー

> **用語注意 — 「オントロジー」の二義性**: 本書の表題・本節でいう「オントロジー」は、データ構造の意味体系（行・列の粒度と、アラインメントスポット／サンプル別ピーク／測定ファイル内ピークの関係）を指す。一方、ARF2/PAI2 等のフィールド `Ontology`（および下のエンティティ表の「Ontology」行）は**脂質クラス**（headgroup 分類。例 `PC`, `Cer_NS`）という別概念である。前者はデータモデル、後者は化学分類であり、混同しない。脂質クラス名・化合物名そのものの記法は §2.1 を参照。

本リポジトリでは、同じ LC-MS データを異なる粒度で表現する。

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

主要エンティティと関係:

| エンティティ | 意味 | 主な識別子・対応 |
|---|---|---|
| Alignment spot | 複数サンプル間で同一とみなされたピーク集合 | ARF/ARF2 の `MasterAlignmentID`、`AlignmentID`、EIC の `spot_id`。ID はすべて 0 始まり |
| Aligned sample peak | 1スポットに属する1サンプルのピーク | ARF の1行。`MasterAlignmentID` と `FileID`/`SampleIndex` の組で特定 |
| Raw-file peak feature | 1測定ファイル中の検出ピーク | PAI2 の1要素。`id` は当該 PAI2 内のピークID |
| MSDec result | デコンボリューション済み MS/MS | DCL の1要素。`dcl_index` と PAI2 のリスト順序が対応する設計 |
| EIC sample trace | 1スポット・1サンプルの抽出イオンクロマトグラム | EIC の `spot_id` と `samples[].file_id` の組で特定 |
| Annotation | 候補化合物名 | `Name`/`name`。空文字、`Unknown`、`no MS2:`、`low score:`を含み得るため、存在するだけで確定同定を意味しない |
| Ontology | 化学・脂質クラス | `Ontology`/`ontology`。例: `FA`, `PC`, `PE`, `TG`, `Cer_NS`。化合物名より上位の分類概念 |

共通の値:

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

**適用はアッセイ種別が `lipid` のときだけ。** MS-DIAL は脂質も一般代謝物も同じ `.arf` / `.mzTab` に書くので、ファイル形式ではどちらか決まらない。一般代謝物（`assay_kind=metabolite`）ではこの文法は成り立たず、同定は候補集合（アダクト・異性体・候補順位）で読む——`lipidmix://docs/output-format/identity` §12.6。種別が未確定のあいだはどちらの規則も当てない。

`Name`/`name` と `Ontology`/`ontology` に現れる脂質表記の読み方。LLM はこの文法に沿ってのみ構造を解釈し、名前から鎖組成やエーテル種別を過剰に推定しない（§9-1/§9-2 と併読）。

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

- **species と molecular species は別粒度**。`PC 34:1` は鎖組成を主張しない。`16:0/18:1` と `17:0/17:1` はいずれも sum 34:1 で、species 名からは判別できない。
- **`P-` と `O-` は species レベルで曖昧**。pygoslin は `PC P-34:0` を `PC O-34:1` に同一化する（プラズマローゲンの二重結合1本 ＝ エーテル鎖の二重結合1本と解釈されるため）。P-/O- を区別したいときは正規化名ではなく `peak_verification.ether_caveats()` の caveat に依拠する（§12.1、[[pe-p-vs-pe-o-annotation]]）。
- **`Ontology` の空文字と文字列 `Unknown` は別扱い**（§9-2）。クラストークンが取れないものを勝手にクラス化しない。
- 同じスポットで ARF の `Name` と ARF2 の `Name` が食い違うことがある（§11.1.1）。統合注釈は ARF2 側を優先しつつ、両者の不一致は所見として明示する。

## 9. LLM解釈時の必須注意事項

1. `Name` が存在しても確定同定とは限らない。空文字、`Unknown`、`no MS2:`、`low score:`を区別する。
2. `Ontology` は化合物名ではなく分類クラスである。空文字と文字列 `Unknown` も区別する。
3. ARFの1行はサンプル別ピーク、ARF2の1行は全サンプル統合スポットであり、同じ「1行」でも粒度が違う。
4. ARFの `IsGapFilled=true` は実測ピークではなく補間値である。検出率や存在判定では別扱いする。
5. ARFの `HeightAverage` は現実装ではグループ先頭値であり、名前どおりの再計算平均ではない。統合統計にはARF2側を優先する。
6. ARF要約の `named_compounds` は空文字も数えるため、注釈率には使わない。
7. PAI2 の `pai2_parser()` は在庫要約のみを返し PCA は行わない（PAI2 は単一サンプルでオミクスPCA不能）。旧ピーク属性PCA系ツール（`pai2_get_top_metabolites` / `pai2_update_analysis_filter`）は撤去済み。複数サンプルの多変量比較は ARF/ARF2 を使う。
8. PAI2の `has_msms=true` は取得参照があることを示すだけで、PAI2内 `msms_spectrum` が非空とは限らない。実スペクトルはDCLを参照する（`pai2_parser` が同名 `.dcl` から自動付与し、`dcl_parser` / `dcl_find_msms` で直接引ける）。**同定確度（MSI Level 2）を MS/MS で主張するときは、`verify_peak_annotation` の `analytical_checks.msms.band` が `PASS`（実スペクトル確認）か `FLAG_ONLY`（フラグのみ）かを必ず見る。**
9. DCLの `n_msms_peaks` は元本数であり、`top_n_peaks` 適用後の配列長とは異なり得る。
10. EICの `peak_top` は強度ではなく頂点座標である。強度は `max_intensity`、平均強度は `mean_intensity` を使う。
11. `total_samples` はEIC全スポットにわたるサンプルエントリ総数であり、ユニークサンプル数は `len(unique_file_ids)` である。
12. m/zとRTの微小差はファイル形式の浮動小数精度、代表値の定義、アラインメント処理に由来し得る。厳密一致ではなく許容差を用いる。
13. **mzTab-M（`dataset_*` 経路）の `abundance_assay[N]` は、非ゼロでも実測ピークとは限らない。** mzTab-M 自体は gap-fill（未検出セルの補間値）を区別する列を持たない。`dataset_load` は隣接する `.arf` が同一アライメントであると数値で検証できた場合にだけ検出状態を取り込む。`dataset_status` の `detection.available` が `false` の間は、**検出率・欠測率・「n 件で検出」といった主張をしてはいけない**。実データでは 60 サンプル × 714 特徴のうち**セルの 70.0% が gap-fill** だったので、非ゼロを検出と数えると検出率を 3 倍以上に過大評価する。検出状態があるときは `dataset_preprocess(min_detection_rate=...)` で実検出率による足切りができる（ARF 経路の同名引数と同義）。

## 状態不足の伝え方（MCP クライアント向けの契約）

多くのツールは、先行するツールが作ったセッション状態に依存する。
`arf_pca_preprocessed` は `arf_preprocess` の行列を、`arf_plot_volcano` は
`arf_differential` の結果を必要とする。

その状態が無いとき、ツールは本文に次の JSON を返す。

```json
{"error": {"code":           "missing_state",
           "state":          "preprocessed_matrix",
           "required_tools": ["arf_preprocess"],
           "message":        "前処理後の行列がありません。先に arf_preprocess を実行してください。"}}
```

- `state` は欠けている状態の識別子。クライアントは**不透明な文字列として扱う**。
- `required_tools` は代替候補（OR）で、クライアントは候補のうち「リプレイ安全かつ本セッションで成功実績のある呼び出し」の中から**直近に成功したもの**を選んで再実行する。リスト順では選ばない。サーバは実際の生成元をすべて列挙する義務がある（取りこぼすと復旧できない、あるいはユーザーが見ていたものと別の状態が復元される）。
- `message` は LLM と人間が読む説明。エンベロープを解釈しないクライアントでは
  これだけが見える。

エンベロープを解釈するクライアントは、`required_tools` のいずれかを
**セッション履歴に記録済みの引数のまま**再実行して状態を復元できる。LLM に
再実行させると前処理の引数が変わって解析条件が黙って変わりうるため、
引数の同一性はクライアント側で担保することが望ましい。

### annotations の解釈

全ツールが MCP 標準の `ToolAnnotations` を宣言する。

- `readOnlyHint=true` は「**サーバの外に副作用が無い**」を意味する。ファイルと
  ネットワークを変更しないこと。**サーバ自身の解析セッション状態の更新は
  副作用に数えない** — その依存関係は上のエンベロープで伝えるため、annotations で
  二重に表現しない。
- `openWorldHint=true` は外部ネットワークへ出ることを示す（`paper_search` のみ）。
- リプレイ安全（同じ引数での再実行が安全）は `readOnlyHint=true` または
  `idempotentHint=true` で判断できる。
