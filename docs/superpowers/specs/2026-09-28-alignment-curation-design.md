# アラインメントのキュレーション（注釈の正しさの一覧確認と機械判別） 設計

- 日付: 2026-09-28
- 発端: 利用者は MS-DIAL GUI で注釈を 1 スポットずつ確かめている——EIC がきれいに描けるか、
  対向プロットが合っているか、脂質クラス×不飽和度でまとめたときに RT–m/z が線形に並ぶか、
  RT・m/z がライブラリから離れすぎていないか。これを本システムで**一覧として**確認でき、
  機械的な判別も付き、判断（フラグ）をシステムへ一斉に送れるようにする。
- 関連: [MS/MS スペクトル照合](2026-09-19-msms-spectral-matching-design.md)、
  [ピーク注釈の検証](2026-07-02-peak-annotation-verification-design.md)、
  [mzTab-M の SML 注釈](2026-09-17-mztab-sml-annotation-design.md)

## 1. 目的と完成条件

キュレーションの入口は 2 つある——(I)「注釈が付いているものは本当に正しいか」、
(II)「注釈が付かないものに妥当な注釈を付ける余地はあるか」。**本 spec は (I) だけを扱う。**
(II) は (I) が作る土台（スポット単位の証拠束・RT–m/z 傾向）の上に後から載せる。

到達点は「判断をこのシステム内に持つ」まで（以下 (b)）。MS-DIAL プロジェクトへの書き戻し（(c)）は
使用感を見てから決める。

完成条件:

1. 対象集合（注釈付き全部 / 脂質クラス指定 / 名前の部分一致）を選ぶと、スポットごとに
   **EIC・対向プロット・Δppm・ΔRT・照合スコア**と、クラス×不飽和度の **RT–m/z 傾向の図**を
   1 つの HTML ビューアに**一覧で**出せる。
2. 同じ数値から機械的な判別（`ok` / `suspect` / `likely_wrong` と理由コード）を出し、
   LLM には数値と理由コードの TSV だけを返す（画像は解析しない）。
3. ビューアで各スポットに「疑わしい」「間違い」のフラグを付け、まとめてシステムへ送れる。
   フラグの無いスポットは「間違っていない」の意で、何も記録しない。
4. 「間違い」を付けたスポットは、差次的エクスポートで同定の無い行として扱われる
   （15 列の契約は変えない）。
5. 実データ（§9）で、GUI が表示する代表の照合結果とスコアを本システムが再現できている。

## 2. 背景と根拠（実測）

上流は `C:\Users\yuu18\source\repos\MsdialWorkbench`（`master`、`afd5f9522`、2026-09-08）。
実データは `C:\Users\yuu18\datasets\a_lipidome_landscape_of_aging_in_mice\rplc\kidney\{neg,pos}`。
以下は実ファイルまたは上流コードで確認した事実で、推測を含まない。

### 2.1 既存資産と足りないもの

| GUI での確認 | 既存 | 足りないもの |
|---|---|---|
| EIC | `eic_plot_chromatograms`（1 スポット・12 試料・点列）、`eic_plot_compounds`（1 試料・複数物質・PNG） | 形状の数値指標。`.EIC.aef` からは頂点・左右 RT・最大/平均強度しか読んでいない |
| 対向プロット | `library_match_feature` → `library_plot_mirror` | スポットとの結び付き（precursor m/z 指定・試料別 `.dcl`）。アラインメント `.dcl` は未使用 |
| RT–m/z 傾向 | なし | 全部。鎖組成は pygoslin（依存に既存）から取れる |
| ライブラリとの差 | 照合の内部で `rt_similarity` / `mass_similarity` を計算 | Δ の出力。MS-DIAL 自身の照合結果（`.arf2` Key 56）は未解析 |

流用できるもの: `lipidmix/msdial/peak_verification.py` の ppm 帯（≤5 PASS / ≤10 BORDERLINE）・
アダクト×イオンモード整合・クラス典型アダクト（いまは `.pai2` の 1 ピーク専用）、
`lipidmix/msdial/lipid_identity.py` の GOSLIN 正規化、`lipidmix/plots/mirror.py` の対向プロット、
`lipidmix/library/store.py` の参照レコード。判断を書き戻す仕組みはどこにも無い。

### 2.2 MatchResults（`AlignmentSpotProperty` Key 56）の構造

`MsScanMatchResultContainer`（`src/MSDIAL5/MsdialCore/DataObj/MsScanMatchResultContainer.cs`）:
Key 0 = `MatchResults`（候補のリスト）、Key 1 = `MSRawID2MspBasedMatchResult`（dict）、
Key 2 = `TextDbBasedMatchResults`（リスト）。

`MsScanMatchResult`（`src/Common/CommonStandard/DataObj/Result/MsScanMatchResult.cs`）の Key:

| Key | 型 | 名前 | Key | 型 | 名前 |
|---|---|---|---|---|---|
| 0 | string | Name | 20 | bool | IsLipidChainsMatch |
| 1 | string | InChIKey | 21 | bool | IsLipidPositionMatch |
| 2 | float | TotalScore | 22 | bool | IsOtherLipidMatch |
| 3 | float | SquaredWeightedDotProduct | 23 | bool | IsRiMatch |
| 4 | float | SquaredSimpleDotProduct | 24 | int | LibraryIDWhenOrdered |
| 5 | float | SquaredReverseDotProduct | 26 | SourceType(byte flags) | Source |
| 6 | float | MatchedPeaksCount | 27 | string | AnnotatorID |
| 7 | float | MatchedPeaksPercentage | 28 | int | SpectrumID |
| 8 | float | EssentialFragmentMatchedScore | 29 | float | AndromedaScore |
| 9 | float | RtSimilarity | 30 | bool | IsDecoy |
| 10 | float | RiSimilarity | 31 | int | Priority |
| 11 | float | CcsSimilarity | 32 | float | PEPScore |
| 12 | float | IsotopeSimilarity | 33 | bool | IsReferenceMatched |
| 13 | float | AcurateMassSimilarity | 34 | bool | IsAnnotationSuggested |
| 14 | int | LibraryID | 35 | bool | IsLipidDoubleBondPositionMatch |
| 15 | bool | IsPrecursorMzMatch | 36 | double | CollisionEnergy |
| 16 | bool | IsSpectrumMatch | 37 | float | EnhancedDotProduct |
| 17 | bool | IsRtMatch | 38 | float | SpectralEntropy |
| 18 | bool | IsCcsMatch | | | |
| 19 | bool | IsLipidClassMatch | | | |

`Source` のビット: Unknown=1, FastaDB=2, MspDB=4, TextDB=16, GeneratedLipid=32, Manual=64。
実データ（neg）の候補は 39 要素で、上表の配置と矛盾しない（MS/MS 無しの候補はドット積が -1）。

**GUI が表示する代表**は `Representative` = 非 decoy の候補の
`argmax(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)`
（`ResultOrder`）。必要な値はすべてシリアライズされているので**忠実に再現できる**。
ここは 09-19 spec が再現しなかった順位付けとは別物で、こちらは既に計算済みの結果から選ぶだけ。

### 2.3 アラインメント `.dcl` の引き方

`AlignmentResult_<ts>.dcl` の結果数はスポット数と等しい（neg: 2207 = 2207）。
`dcl[MasterAlignmentID]` の precursor m/z がスポットの MassCenter と ±0.01 で一致するのは
2207 中 2198。`MSDecResultIdUsed`（Key 59）は -1 が 294、MasterAlignmentID と一致が 2、
それ以外が 1911 で、索引としての意味は未確定。**残り 9 件の不一致と Key 59 の意味は
plan の最初で上流コードから確定させる**（§8）。

### 2.4 注釈名の接頭辞

実データ neg の注釈付き 1468 件のうち `low score:` 380、`no MS2:` 86。既存の
`normalize_lipid_name` はこれを剥がして捨てているが、判別には強い信号なので保持して使う。

## 3. 全体の流れ

```
curation_review(selection, …)
  ├─ 証拠の収集（スポット単位）: arf2 のスポット + MatchResults の代表 / アラインメント .dcl の代表 MS/MS /
  │    ライブラリの参照（LibraryID で引く）/ .EIC.aef の試料別系列 / .arf の試料別行（検出か gap-fill か・S/N）
  ├─ 機械判別（数値に対して）→ 各系統の帯 + 総合判定 + 理由コード
  ├─ RT–m/z 傾向（クラスごと）
  ├─ HTML ビューアを <データフォルダ>\curation\ へ書き出す
  └─ 戻り値: review_id / 件数 / suspect 以上の TSV / クラス別の傾向要約 / HTML のパス

ビューア ──(MCP Apps で表示)──→ 「送信」→ curation_submit を直接呼ぶ
        └─(ブラウザで開く)────→ 「送信用テキストをコピー」→ チャットへ貼る → LLM が curation_submit を呼ぶ
```

## 4. 機械判別

入力はすべて数値。各系統を `PASS` / `BORDERLINE` / `FAIL` / `UNKNOWN` の帯と理由コードにする。
**`UNKNOWN` は `FAIL` に数えない**（「MS/MS 未取得」と「合わなかった」を混ぜない既存規約）。

| 系統 | 入力 | 指標 |
|---|---|---|
| ① MS/MS 照合 | 代表の MatchResult（§2.2）、アラインメント `.dcl` の代表スペクトル、参照スペクトル | Squared{Weighted,Simple,Reverse}DotProduct、MatchedPeaksPercentage、`IsSpectrumMatch`、脂質の `IsLipidClassMatch` / `IsLipidChainsMatch`、名前の接頭辞（`no MS2:` `low score:` `w/o MS2:` `unsettled:`）。本システムの `spectral_match` で再計算した値を併記し、食い違いを理由コードにする |
| ② Δm/z | MassCenter と参照の precursor m/z（参照が無ければ組成式＋アダクトの理論値） | ppm 誤差（既存の帯）。アダクト×イオンモード整合（既存） |
| ③ ΔRT | スポットの RT と参照 RT | 分単位の差。参照 RT が無いライブラリでは `UNKNOWN` |
| ④ EIC 形状 | `.EIC.aef` の試料別系列、`.arf` の試料別行 | 試料ごと: 頂点が積分範囲内か・範囲内の点数・ガウス当てはめの R²・左右対称性・ギザつき（範囲内の極大数）。スポットへ集約: 検出試料のうち形状 PASS の割合、試料間の頂点 RT のばらつき |
| ⑤ RT–m/z 傾向 | 同じ review 内の注釈付きスポット全体（クラス×不飽和度は pygoslin の和組成） | クラスごとに 1 本の加法モデル `RT = a + b·C + c·DB` を頑健回帰し、その残差。点数不足のクラスは `UNKNOWN`。群ごとの R² と点数を併記 |

**⑤ は判断材料の 1 つにとどめる。** 不飽和度の上昇で線形性は下がるため、⑤ 単独で総合判定を
`suspect` 以上にしない。①〜④ のいずれかと重なったときだけ補強として効かせる。当てはまりの悪い
群（R² が低い・点数が少ない）の残差は弱めて表示する。図はクラスごとに不飽和度で線を分けた
RT–m/z 散布図とする（利用者が GUI で見ている形）。判定に加法モデルを使うのは、1 つの不飽和度に
2〜3 点しかない群でも同じクラスの他の不飽和度から傾きを借りられるため。

**総合判定**: `likely_wrong` = 強い `FAIL` が 1 つ以上（例: イオンモード矛盾、Δppm > 10、
MS/MS があるのにスペクトル不一致、`IsLipidClassMatch=false`）。`suspect` = 弱い `FAIL` か
`BORDERLINE` の重なり。それ以外は `ok`。強弱の割り当ては plan で表にして固定する。

**しきい値**: 既定値はその解析自身の設定から取る（`.dbs` の search_params、Console メソッドの
MS1 許容幅・RT 許容幅・同定スコア足切り）。無ければ固定の既定値。すべて引数で上書きでき、
実際に使った値と出所を戻り値に明記する。

## 5. ツールとビューア

### 5.1 ツール（3 つ追加、ほかにアプリ専用 1 つ）

| ツール | 役割 |
|---|---|
| `curation_review(selection, file_ids=None, thresholds=None, …)` | `selection` は `"annotated"` / `{"ontology": [...]}` / `{"name_contains": "..."}`。証拠収集・判別・ビューア生成を 1 回で行う。戻り値は `review_id`・判定の件数・`suspect` 以上の TSV（理由コード・Δppm・ΔRT・スコア）・クラス別の傾向要約（点数・R²・外れ数）・HTML のパス。前提（arf2 未ロード等）が無ければ `missing_state` 封筒 |
| `curation_submit(review_id, flags)` | `flags` の各要素は `{spot_id, flag: "wrong" \| "suspect" \| "clear", note?}`。送り主（`user` / `llm`）を記録する |
| `curation_flags(…)` | 現在有効なフラグ（スポットごとの最新 1 行）を TSV で返す |
| `curation_view_data(review_id, page)` | MCP Apps のビューアが表示後にデータをページ単位で取る。`_meta.ui.visibility = ["app"]` で LLM には見せない |

LLM は機械判別の結果を読んでユーザーに提案する。**LLM が `curation_submit` を呼んでよいのは、
ユーザーがチャットで同意したとき**と、ユーザーが貼った送信用テキストを取り次ぐときだけ。

戻り値の規約（`json_payload`、TSV、float の丸め、`structured_output=False`）は CLAUDE.md に従う。
EIC 系列・スペクトル座標は戻り値に入れず、ビューアと session に持つ。

### 5.2 ビューア

1 つの HTML で、データの渡し方を 2 通り持つ。

- **ブラウザで開く**: データを HTML に埋め込む。送信は「送信用テキストをコピー」→ チャットに貼る。
  フラグを付けたものだけを送るので小さい。
- **MCP Apps（`ui://` リソース）で表示**: HTML は静的リソースとして配信し、データは
  `curation_view_data` から取る。送信は `curation_submit` を直接呼ぶ。
  Claude Desktop はローカル stdio サーバの `_meta.ui` を落とす不具合が 2026-09-22 以降 OPEN
  （anthropics/claude-ai-mcp#1069）なので、当面はブラウザ経路が主。

画面構成: 上部に絞り込み（クラス・判定・フラグの有無）と送信ボタン。その下にクラスごとの傾向
パネル（外れ点をクリックすると該当カードへ）。その下にスポットカードの格子——名前・アダクト・
判定バッジ＋理由コード・Δppm・ΔRT・ドット積、小さな EIC（代表試料は太線、gap-fill は破線、
積分範囲は網掛け）、小さな対向プロット、フラグ（なし / 疑わしい / 間違い）とメモ欄。
カードはクリックで拡大する。

描画は外部 CDN を使わず、埋め込みの小さな canvas 描画にする（MCP Apps の iframe は CSP で
外部読み込みが止まりうる）。数百枚のカードは画面に入った分だけ描く。カードの EIC は既定で
代表試料＋強度上位の検出試料で最大 12 本（`eic_plot_chromatograms` と同じ上限）、`file_ids` で
指定し直せる。

## 6. 判断の記録と下流への反映

- **記録するのはユーザーが付けたフラグだけ**（`wrong` / `suspect`、取り消しは `clear`）。
  フラグの無いスポットは「間違っていない」で、何も書かない。「確認済み」「未確認」は持たない。
- 置き場所は `<データフォルダ>\curation\` の追記専用ファイル（JSON Lines）。キーは
  アラインメントファイルの識別子（ファイル名＋sha256）と `MasterAlignmentID` の組。
  同じスポットは最新の行が勝つ。ローカル単独利用が前提で、共有は考えない。
- **差次的エクスポート**（`arf_export_differential` / `dataset_export_differential`）:
  `wrong` のスポットは同定を外して扱う。エクスポートは InChIKey の付いた行だけを出すので、
  結果としてその行は落ちる。`suspect` は行も値も変えない。
  **15 列の契約（`CONTRACT_VERSION = 1`）は変えない**。代わりにメタ行
  `# curation = applied; wrong_excluded=…; suspect=…; flags_sha256=…` を足す。
  引数 `apply_curation`（既定 true）で切れる。フラグが 1 件も無ければ出力は現行と完全に同じ。
- mzTab-M 経路のスポットとの対応は既存の `.arf`↔SMF 対応（`lipidmix/mztab/evidence.py`）を使う。
  対応が取れないフラグは適用せず、件数をメタ行に出す。
- `arf2_annotate_identities` の出力にもフラグの印を付ける。
- 書き戻し（(c)）は範囲外。キーが `MasterAlignmentID` なので後から経路を足せる。

## 7. 範囲外

- (II) 未注釈スポットへの候補提示（RT–m/z 傾向による RT 予測・MS/MS 照合の組み合わせ）。
- MS-DIAL プロジェクトへの書き戻し。
- Use-LLLM の WebUI でのビューア表示（iframe ホストとツール呼び出しの中継が向こう側に要る）。
- 画像解析による判別。
- `CONTRACT_VERSION` の引き上げを伴うエクスポート列の追加。
- pipeline（`pipeline_run`）への stage 追加。

## 8. 実装前に確定させる事実（plan の先頭タスク）

1. MsScanMatchResult の Key 表（§2.2）を `docs/schema/MsScanMatchResult.md` に起こす
   （上流コミットと欠番 25 を記載）。
2. アラインメント `.dcl` の索引: §2.3 の不一致 9 件と `MSDecResultIdUsed` の意味を上流コード
   （アラインメント `.dcl` を書く側と GUI が読む側）で確定する。
3. `LibraryID` と `.dbs` / `.msp` の参照レコード（`library/store.py` の `library_id` /
   `record_index`）との対応を上流コードと実データで確定する。
4. `.EIC.aef` の試料別系列と `.arf` の試料別行（検出 / gap-fill）の対応付けを実データで確かめる。

## 9. テストと検証

- **fixture はテスト自身が作る**: 合成の `.arf2`（MatchResults 入り）・アラインメント `.dcl`・
  `.EIC.aef`・`.arf`。研究データや研究室ライブラリをテスト・fixture・commit に含めない。
- 純関数の単体テスト: 代表の選択（`ResultOrder`）、形状指標、Δppm / ΔRT、加法モデルの残差と
  R²、帯と総合判定（`UNKNOWN` を `FAIL` に数えない、⑤ 単独で `suspect` にならない）。
- 判断記録: 追記専用、最新が勝つ、`clear` で取り消し。
- エクスポート: フラグ 0 件で現行と完全一致（回帰）、`wrong` の行が落ちる、メタ行。
- ビューア: 埋め込み JSON の形の検証、送信用テキストの往復（生成 → `curation_submit` が受理）。
- 腐敗防止: `USAGE.md`、`docs/workflow/curation.md`（新設）、CLAUDE.md 冒頭の規模表記、
  登録数と `ToolAnnotations`（`tests/test_server_registration.py`）。
- 実データ検証（kidney neg/pos）: 代表の MatchResult の Name・TotalScore が arf2 の Name と
  整合すること、アラインメント `.dcl` の precursor 一致率、判定の分布の目視（内蔵ブラウザで HTML）。
  MCP Apps での表示は Desktop の不具合が直るまで**未検証**と明記する。
- vault 側の流れ図にキュレーションの枝を足す。
