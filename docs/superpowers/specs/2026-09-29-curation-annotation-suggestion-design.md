# キュレーションの候補付け（間違い・未注釈スポットへの注釈候補の提示と採用の記録） 設計

- 日付: 2026-09-29
- 発端: [アラインメントのキュレーション](2026-09-28-alignment-curation-design.md) は入口 (I)
  「注釈が付いているものは本当に正しいか」だけを扱い、(II)「注釈が付かないものに妥当な注釈を
  付ける余地はあるか」を後回しにした。本 spec は (II) と、(I) で**間違い**になったスポットの
  付け直しを、(I) の土台（スポット単位の証拠束・RT–m/z 傾向・判断の記録・ビューア）の上に載せる。
- 関連: [MS/MS スペクトル照合](2026-09-19-msms-spectral-matching-design.md)

## 1. 目的と完成条件

到達点は**判断をこのシステム内に持ち、本システムの下流（差次的エクスポート）で使う**まで（段階 A）。
質を見てから、MS-DIAL への注釈名の書き戻し（B）と、再解析の判断材料としての利用（C）へ広げる。

候補の裏付けは**スペクトル類似度で並べ、制約で削る**方式（P）で始める。MS-DIAL の脂質規則
（診断イオンによるクラス確証・鎖決定）は移植しない。主要クラスの診断イオン確認（Q）へ進むかは、
§11.1 の数字を見て決める。

完成条件:

1. 対象（§4）のスポットごとに、候補源 ①（MS-DIAL 自身の下位候補）・②（閾値を緩めた再検索）・
   ④（注釈付きの別スポットのイオンとしての説明）の候補を、1 つの HTML ビューアに**一覧で**出せる。
2. LLM には候補の数値と理由コードの TSV だけを返す（スペクトル座標・EIC 系列は返さない）。
3. ビューアで各スポットに「候補を採用（和組成／分子種）」「Y の別イオン」「元の注釈に戻す」
   「選ばない」を選び、まとめてシステムへ送れる。
4. 採用した候補は差次的エクスポートで同定を置き換え、「別イオン」は除外される
   （15 列の契約は変えない）。判断が 0 件なら出力は現行と完全に同じ。
5. 実データ（kidney neg/pos）で、§11.1 の正解ありの検証の数字を出す。

## 2. 背景と根拠（実測）

実データは `C:\Users\yuu18\datasets\a_lipidome_landscape_of_aging_in_mice\rplc\kidney\{neg,pos}`、
上流は `C:\Users\yuu18\source\repos\MsdialWorkbench`。以下は 2026-09-29 に実ファイルまたは上流
コードで確認した事実。

### 2.1 候補源の有無

| 事実（kidney neg、2207 スポット） | 値 |
|---|---|
| 注釈付き（代表あり） | 1468 |
| 未注釈 | 739。MatchResults は `Source=Unknown`・`LibraryID=-1`・スコア 0 のダミー 1 件だけ |
| 未注釈のうちアラインメント `.dcl` に MS/MS あり（≥3 ピーク / 1–2 ピーク / なし） | 456 / 75 / 208 |
| 未注釈の m/z（最小 / 四分位 / 中央 / 四分位 / 99%） | 75 / 297 / 499 / 686 / 1234 |
| 注釈付きの候補数（非 decoy・非 Unknown が 1 / 2 / 3 件） | 209 / 61 / 1198 |
| ユーザーの有効な `wrong` フラグ | 31 |

**MS-DIAL は足切りを下回った候補を保存しない**ので、未注釈スポットの候補は自分で再検索しない限り
得られない。注釈付きスポットには MS-DIAL 自身の 2 位・3 位が残っている。

### 2.2 使われた参照ライブラリと解析設定

- 参照はデータフォルダの `Dataset_<ts>_Loaded.msp2.dbs`（MS-DIAL 内蔵の脂質ライブラリを
  検索対象クラスで絞ったもの）。param ファイル（`Dataset_<ts>_param_<ts>.txt`）の
  `Msp file path` / `Lbm file path` は空。
- param ファイルの値: `MS1 tolerance for centroid: 0.01`、`MS2 tolerance for centroid: 0.025`、
  `Retention time tolerance for alignment: 0.1`、`MS1 tolerance for alignment: 0.015`、
  `Searched adduct ions:`（neg は `[M-H]-,[M-H2O-H]-,[M+Na-2H]-,[M+Cl]-,[M+K-2H]-,[M+HCOO]-,
  [M+CH3COO]-,…,[2M-H]-,[2M+FA-H]-,[2M+Hac-H]-,[3M-H]-,[M-2H]2-,[M-3H]3-`）。

### 2.3 既存の MS/MS 検証（再利用するものと、足りないもの）

| 場所 | 何をしているか |
|---|---|
| `analysis/spectral_match.py` `match_spectrum` / `total_score` | MS-DIAL の汎用スペクトル採点の移植（3 種のドット積・一致ピーク率・エントロピー）と、非正規化の総合スコア。**脂質規則（上流 `Lipidomics/`）は移植していない**ため、化学的にありえない候補が上位に残りうる（`docs/output_format/library.md` §14.9） |
| `library/tools.py` `library_match_feature` | 試料別 `.dcl` の 1 スペクトルで store を再検索し、`total_score` 順に上位 N 件を返す |
| `pai2/tools.py` `verify_peak_annotation` → `msdial/peak_verification.py` | `.pai2` の 1 ピークについて ppm・アダクト×極性・MS/MS の有無を判定し、最良 1 件だけ照合する |
| `curation/evidence.py` → `curation/judge.py` | 再検索しない。MS-DIAL の代表の照合結果（`IsReferenceMatched`・脂質規則フラグ）を主にし、代表の参照 1 件との再採点は食い違いの検出にだけ使う |

いまの「間違い」判定の強い根拠は MS-DIAL の脂質規則の**結果を読んでいる**だけなので、こちらの
再検索で出た新しい候補にはその裏付けが無い。§5.3 の制約はこの穴を埋める最小限の代わり。

### 2.4 MS-DIAL のイオン関係（`.arf2` Key 10 `IonFeatureCharacter`）

上流 `src/MSDIAL5/MsdialCore/DataObj/IonFeatureCharacter.cs`: Key 0 `AdductType`、1
`AdductTypeByAmalgamationProgram`、2 `Charge`、3 `PeakLinks`（`LinkedPeakFeature` =
Key 0 `LinkedPeakID` / Key 1 `Character`）、4 `IsotopeWeightNumber`、5 `IsotopeParentPeakID`、
6 `PeakGroupID`、7 `IsLinked`、8 `AdductParent`。`PeakLinkFeatureEnum` =
`SameFeature(0), Isotope(1), Adduct(2), ChromSimilar(3), FoundInUpperMsMs(4), CorrelSimilar(5)`
（`src/Common/CommonStandard/Enum/CommonEnums.cs`）。

| 実データ | neg | pos |
|---|---|---|
| `IsotopeWeightNumber` | 全スポット 0 | 全スポット 0 |
| リンクを持つ未注釈スポット | 531 / 739 | 698 / 1507 |
| 未注釈 → 注釈付きへのリンク（ChromSimilar / FoundInUpperMsMs / CorrelSimilar / Adduct） | 31 / 150 / 522 / 0 | 34 / 202 / 531 / 2 |

アラインメント段階では同位体は別スポットとして残らない。ただし **M+2 が同族体の単同位体ピークに
重なる場合（Type II の重なり）** はこの仕組みで拾えない。

### 2.5 InChIKey の付き方

代表の照合結果の名前と InChIKey はライブラリのレコードのもの（neg: InChIKey あり 1183 / なし 285、
なしはすべて `RIKEN … from …` 名の未知スペクトル）。MS-DIAL は和組成の名前
（`SM 34:2;O3`、`Hex2Cer 34:2;O2`）にも、ライブラリが代表として持つ構造の InChIKey を付けている。
§8.2 の和組成での採用はこの振る舞いに倣う。

## 3. 全体の流れ

```
curation_suggest(review_id=None, wrong=…, unannotated=…, …)
  ├─ 元レビューの読み込み（likely_wrong の判定を読む）＋ 有効フラグ
  ├─ 対象の選定（§4）
  ├─ 証拠の収集（既存 evidence.collect を再利用: 代表試料の m/z・RT、アラインメント .dcl、EIC）
  ├─ ① MS-DIAL の下位候補 ＋ ② 再検索 → 採点 → 制約（§5）
  ├─ ④ イオン関係の探索（§6）
  ├─ HTML ビューアを <データフォルダ>\curation\ へ書き出す
  └─ 戻り値: suggestion_id / 件数 / 1 スポット 1 行の TSV / HTML のパス

ビューア ──(ブラウザ)──→「送信用テキストをコピー」→ チャットへ貼る → LLM が curation_submit を呼ぶ
```

## 4. 対象

- `wrong`（引数）: `"flagged_or_likely"`（既定）= 有効な `wrong` フラグ ∪ 元レビューの判定が
  `likely_wrong` のスポット。`"flagged"` = 有効な `wrong` フラグだけ。
- 元レビュー: `review_id` を省くと、そのアラインメント（sha256 一致）の最新のレビュー。1 つも無ければ
  `missing_state`（`required_tools=["curation_review"]`）。sha256 が現在の `.arf2` と違えばエラー。
  `likely_wrong` は元レビューの対象範囲についてしか分からないことを戻り値に明記する。
- `unannotated`（引数、既定 true）: `.arf2` の代表が無いスポット全部。元レビューの範囲に縛られない。
  MS/MS の無いスポットには ② が働かず ① も無いので、事実上 ④ だけが候補を出す。
- 有効な判断が `assign` / `redundant` のスポットは対象から外す（引数 `include_decided=false`）。

## 5. 候補の生成（① と ②）

### 5.1 測定側

既存 `evidence.collect` と同じく、代表試料の m/z・RT（無ければ MassCenter / RT）と
`dcl[MasterAlignmentID]` を使う。precursor が ±0.01 で一致しなければ MS/MS 無しとして扱う。

### 5.2 候補源

- **① MS-DIAL の下位候補**: MatchResults のうち非 decoy・非 Unknown で、代表を除いたもの。参照は
  `store.record_by_scan_id` で引き、`match_spectrum` / `total_score` で採点し直して ② と同じ物差しに
  揃える。参照が引けない候補は MS-DIAL の値だけで載せ、`reference_unresolved` を付ける。
- **② 閾値を緩めた再検索**: `store.candidates(rep_mz, mz_tol=ms1_tolerance, ion_mode=スポットの極性)`
  を**スコアの足切りなし**で全件採点する。採点の前処理パラメータは `library_match_feature` と同じ集合
  （store の `search_params`）を使う。MS/MS が無いスポットには行わない。
- **ライブラリ**: 既存 `curation_review` と同じく `library_load` 済みの store を使い、無ければ
  `missing_state`（`required_tools=["library_load"]`）。メッセージではデータフォルダの
  `*_Loaded.msp2.dbs`（MS-DIAL が実際に検索した集合）を推奨する。使った store のパスと sha256 を
  戻り値に出す。
- **統合**: 同じレコード（`library_id` と `record_index` が同じ）は 1 つにまとめ、出所を
  `msdial+research` とする。`wrong` フラグのスポットでは、いまの代表と同じレコードを候補から外す。

### 5.3 制約（P）

| 種類 | 理由コード | 条件 |
|---|---|---|
| ハード（削る） | `polarity_mismatch` | 候補のアダクトとスポットの極性が矛盾（既存 `adduct_consistency`） |
| ハード（削る） | `dmz_out` | \|Δm/z\| ≥ 10 mDa（既存 `dmz_fail_mda`） |
| ソフト（順位を下げる） | `adduct_atypical` | そのクラスに典型的でないアダクト（既存 `adduct_consistency` の BORDERLINE） |
| ソフト | `trend_outlier` | RT–m/z 傾向からの残差 \|z\| > `trend_outlier_z` |
| ソフト | `no_matched_peaks` | `matched_peaks_count` = 0 |
| ソフト | `msms_absent` | スポットに MS/MS が無い（① の候補にだけ付く） |
| 情報 | `trend_unknown` | 候補のクラスの傾向が点数不足で当てはまらない。**順位は下げない** |
| 情報 | `reference_unresolved` | ① の候補の参照レコードが store に無い |

**傾向モデル**: そのアラインメントの注釈付きスポットのうち、対象でも `wrong` フラグ付きでもないものから、
既存 `trend.fit_trends`（クラスごとの `RT = a + b·C + c·DB`、Huber）で当てはめる。候補のクラス×和組成
（`trend.composition`）で予測 RT を出し、残差を z にする。元レビューの範囲に縛らない（名前だけで
当てはめられるため）。

**並び順**: ハード制約を通った候補を（ソフト理由の数 昇順, `total_score` 降順）で並べ、スポットごとに
上位 `top_n`（既定 5）件を出す。

### 5.4 記録の粒度

既定は**和組成**（`PC 16:0_18:1` → `PC 34:1`。`trend.composition` と同じ pygoslin の解析）で記録し、
ビューアで分子種に切り替えられる。P は鎖組成を判定しないため。

## 6. ④ イオン関係による説明

### 6.1 判定

対象スポット X と相手 Y の組を、次の 3 条件がそろったときに候補にする。

1. **同時溶出**: \|ΔRT\| ≤ `rt_window`。既定は param ファイルの
   `Retention time tolerance for alignment`（無ければ 0.1 分）。
2. **Δm/z の説明**（許容 `relation_mz_tol`、既定 10 mDa）。関係は次の表から作る。
   - **アダクトの組**: Y の中性質量（Y の m/z とアダクトから逆算）に、その解析の `Searched adduct ions`
     を付けた m/z。二量体・多価はこの一覧から出る。param ファイルが無ければ極性ごとの既定一覧を使う。
   - **同位体**: M+1（+1.003355）・M+2（+2.006710）。強度比 X/Y を、Y の組成式の炭素数から二項分布で
     求めた期待比と比べる。期待比の 1.5 倍以下なら `isotope_consistent`、それを超えれば X は実在の別物質
     とみなし、関係は情報として添えるだけにする（Type II の重なりの判定）。
   - **インソース断片**: 下表の中性損失（表は `lipidmix/curation/relations.py` の定数。出典は plan の
     先頭タスクで確定し、表に併記する）。

     | 関係コード | 質量差 | 想定 |
     |---|---|---|
     | `insource:-H2O` | 18.010565 | 脱水（Cer・DG など） |
     | `insource:-2H2O` | 36.021129 | 二重の脱水（Cer） |
     | `insource:-NH3` | 17.026549 | [M+NH4]+ からの脱アンモニア |
     | `insource:-HCOOCH3` | 60.021129 | PC・SM の [M+HCOO]- → [M-CH3]- |
     | `insource:-CH3COOCH3` | 74.036779 | PC・SM の [M+CH3COO]- → [M-CH3]- |
     | `insource:-C3H5NO2` | 87.032028 | PS の頭部基（セリン）脱離（neg） |
     | `insource:-C2H8NO4P` | 141.019094 | PE の頭部基脱離（pos） |
     | `insource:-C3H8NO6P` | 185.008923 | PS の頭部基脱離（pos） |
     | `insource:-C6H10O5` | 162.052824 | ヘキソース脱離（HexCer・PI など） |

   - **MS-DIAL の `FoundInUpperMsMs` リンク**: 相手が Y なら、Δm/z が表に無くても
     `found_in_upper_msms` の候補にする。
3. **裏付け**: MS-DIAL のリンク（`CorrelSimilar` / `ChromSimilar` / `FoundInUpperMsMs`）があるか、
   `.arf` の試料別の高さの log1p で求めた Pearson の r ≥ `relation_min_r`（既定 0.8、有効試料 5 以上）。
   どちらも無ければ候補には出すが、`no_support` を付けて順位を下げる。

**相手 Y の条件**: 注釈付きで、有効な `wrong` フラグが無く、今回の対象でもないスポット。誤った注釈を
起点にした説明の連鎖を防ぐ。

### 6.2 強い説明

`found_in_upper_msms`、または `isotope_consistent` の M+1/M+2 で裏付けありのものを**強い説明**とする。
④ の候補は ①② と物差しが違うので 1 本の順位に混ぜず、強い説明があるスポットではそれをカードの先頭に
出し、初期選択にする。

## 7. ツールとビューア

### 7.1 ツール

新規は `curation_suggest` の 1 つ。送信と一覧は既存の `curation_submit` / `curation_flags` を広げる。

`curation_suggest(review_id=None, wrong="flagged_or_likely", unannotated=True, include_decided=False,
top_n=5, rt_window=None, relation_mz_tol=None, relation_min_r=None, thresholds=None, file_path=None,
max_rows=100)`（`file_path` と `max_rows` の意味は `curation_review` と同じ）

- 結果は `suggestion_id`（接頭辞 `cs-`）を持ち、`<データフォルダ>\curation\suggest-<id>.json` と
  `.html` に保存する。
- 戻り値: `suggestion_id`、元レビューの ID、ライブラリのパスと sha256、使ったしきい値と出所、件数
  （対象の種類 `flagged` / `likely_wrong` / `unannotated` ごとの対象数、候補のあったスポット数、
  強い説明のあるスポット数）、1 スポット 1 行の TSV（`spot_id`、今の名前、対象の種類、最上位の候補、
  出所、`total_score`、理由コード、④ の関係）、HTML のパス。
- 前提が無ければ `missing_state`（arf2 未ロード、元レビューが無い、`library_load` 未実行）。
- 戻り値の規約（`json_payload`、TSV、丸め、`structured_output=False`）は CLAUDE.md に従う。
  スペクトル座標・EIC 系列は JSON と HTML にだけ置く。
- LLM が `curation_submit` を呼んでよいのは、ユーザーがチャットで同意したときと、ユーザーが貼った
  送信用テキストを取り次ぐときだけ（既存と同じ）。

### 7.2 ビューア

- 既存 `viewer.html` の EIC・対向プロットの canvas 描画を共通 JS に切り出し、`render_html` が両方の
  テンプレートに埋め込む。候補付けは別テンプレート `suggest_viewer.html`。外部 CDN は使わない。
- 上部: 絞り込み（対象の種類、強い説明の有無、候補の出所、選択済みか）と「送信用テキストをコピー」。
- カード（1 スポット 1 枚）:
  - 見出し: `spot_id`、m/z、RT、今の注釈と、対象になった理由（`wrong` フラグ / `likely_wrong` と
    理由コード / 未注釈）。
  - 小さな EIC（既存と同じ描き方）。
  - 候補の一覧（ラジオ）。①② の行: 名前（既定は和組成、「分子種で記録」の切り替え）、アダクト、
    出所バッジ（`msdial` / `research` / `msdial+research`）、`total_score`、ドット積、Δm/z、
    傾向の残差、理由コード。選ぶとその候補の対向プロットを出す。
  - ④ の行: 相手 Y の名前・関係・Δm/z・ΔRT・r。選ぶと X と Y の EIC を重ねて出す。
  - 末尾に「選ばない」（何も記録しない）。`wrong` フラグのスポットには「元の注釈に戻す」（`clear`）。
  - メモ欄（根拠を自動で下書き）。
- 初期選択は強い説明のある ④ の行だけで、それ以外は「選ばない」。
- 送信用テキストは `suggestion_id` と、スポットごとの `{spot_id, flag, candidate, level?, note?}`
  だけを持つ（候補の詳細は持たない）。`curation_submit` が保存済みの `suggest-<id>.json` から候補の
  詳細を引いて記録行に展開する。テキストを小さく保ち、貼り付け時の改変で別の名前が記録されるのを防ぐ。
- MCP Apps での表示は Desktop の不具合（anthropics/claude-ai-mcp#1069）が直るまで未検証とする。

## 8. 判断の記録と下流への反映

### 8.1 記録

`flags.jsonl` に flag の種類を 2 つ足す。キー（アラインメントの sha256、`MasterAlignmentID`）、
「スポットごとに最新の 1 行が勝つ」、`clear` による取り消しは既存どおり。ファイル名は変えず、
「フラグ」は「判断の記録」の意味に広がる。

| flag | 主な項目 |
|---|---|
| `assign` | `name`（記録名）、`level`（`sum` / `species`）、`species_name`、`ontology`、`adduct`、`formula`、`inchikey`、`candidate_source`（`msdial` / `research` / `msdial+research`）、`library`（ファイル名・sha256・`library_id`・`record_index`）、採点のスナップショット、`suggestion_id` |
| `redundant` | `of`（Y の `spot_id`）、`relation`（例 `isotope_M+2`・`adduct:[M+HCOO]-/[M-H]-`・`insource:-H2O`・`found_in_upper_msms`）、根拠（Δm/z・ΔRT・r・MS-DIAL のリンク） |

- 「`wrong` の後に `assign`」は `assign` が勝つ。
- `validate_entries`・送信用テキストの書式・`curation_flags` の TSV を広げる。

### 8.2 差次的エクスポート（両経路、`curation/apply.py` の 1 か所）

- `wrong`: 除外（現行）。`redundant`: 除外。
- `assign`: 行の `name`・`ontology`・`inchikey` を置き換え、`name_source` と `inchikey_source` を
  `curation` にする。`msi_level` は記録名から既存 `build_identity_block` で求める。和組成で採用したときも、
  選んだレコードの InChIKey を持たせる（§2.5 の MS-DIAL の振る舞いと同じ）。
- 未注釈だったスポットに `assign` が付くと、そのスポットはエクスポートに**新しく現れる**。
- メタ行に `curation_assigned` と `curation_redundant_excluded` を足す。15 列と `CONTRACT_VERSION = 1`
  は変えない。判断が 0 件なら出力は現行と完全に同じ。
- mzTab-M 経路は既存の `.arf`↔SMF の対応を使い、対応の取れない判断は適用せず件数をメタ行に出す。
- `arf2_annotate_identities` の出力に `assign` / `redundant` の印を付ける。

### 8.3 `_tags.xml` への書き戻し

`assign` / `redundant` は Misannotation を変えない（MS-DIAL の中の名前はまだ誤ったままのため）。
Misannotation を外すのは B（注釈名そのものの書き戻し）の仕事。

## 9. 範囲外

- B（MS-DIAL への注釈名の書き戻し）、C（再解析の判断材料）。
- Q（主要クラスの診断イオン確認）。
- ③（MS1＋傾向だけによる和組成の提案。MS/MS の無いスポットへの ①② 以外の候補）。
- pipeline（`pipeline_run`）への stage 追加。
- MCP Apps での表示の検証、Use-LLLM の WebUI での表示。
- `CONTRACT_VERSION` の引き上げを伴うエクスポート列の追加。

## 10. 実装前に確定させる事実（plan の先頭タスク）

1. `.arf2` Key 10（`IonFeatureCharacter`・`LinkedPeakFeature`・`AdductIon`）の配置を
   `docs/schema/IonFeatureCharacter.md` に起こす（上流コミット・欠番を記載）。`LinkedPeakID` が
   `MasterAlignmentID` を指すことを実データで確かめる。
2. `Loaded.msp2.dbs` を store で開けること、① の参照が `record_by_scan_id` で引けることを確かめる。
3. param ファイルの `Searched adduct ions` と `Retention time tolerance for alignment` の読み方
   （ファイルの選び方: `Dataset_<ts>_param_*.txt` のうち最新）と、無いときの既定を決める。
   アダクト文字列から中性質量と m/z を相互に計算する関数（既存 `adduct_mz` の対応範囲の確認と拡張）。
4. §6.1 の中性損失表の出典を確定し、表に併記する。
5. massbank-context が `name_source = curation` という新しい値と新しいメタ行のキーをそのまま
   読めることを確かめる。

## 11. テストと検証

### 11.1 候補の質の測り方（実データ、kidney neg/pos）

- **正解ありの検証（自動、スクリプトとして残す）**: MS-DIAL が確信を持って付けた注釈（元レビューの判定
  `ok`、`IsReferenceMatched`、フラグ無し）を正解とみなし、代表を隠して ② を回し、正解のレコードが
  1 位・3 位以内に入る割合を出す。和組成で比べた割合も併記する。
- **④ の再現**: MS-DIAL の `FoundInUpperMsMs` リンクのうち、こちらの Δm/z 表が関係として名指しできた
  割合。
- **利用者の判断**: `wrong` と `likely_wrong` をビューアで見てもらい、候補を採用した割合と採用した候補の
  順位を記録から数える。Q へ進むかはこの数字で決める。

### 11.2 自動テスト（fixture はテスト自身が作る）

研究データと研究室ライブラリをテスト・fixture・commit に含めない。

- 純関数: ①② の統合、ハード／ソフト制約、並び順、傾向の予測 RT と `trend_unknown`、④ の各関係
  （アダクトの組・M+2 と強度比・インソース・リンク起点）、相関、相手 Y の条件（対象や `wrong` を相手に
  しない）、和組成への変換。
- 記録: `assign` が `wrong` に勝つ、`clear` で全部消える、送信用テキストの往復（生成 →
  `curation_submit` が保存済み候補から展開して受理）、存在しない候補 ID の拒否。
- エクスポート: 判断 0 件で現行と完全一致（回帰）、`assign` で同定が置き換わり未注釈の行が現れる、
  `redundant` が除外される、メタ行、両経路で同じ結果。
- ビューア: 埋め込み JSON の形。
- 腐敗防止: `USAGE.md`、`docs/workflow/curation.md`、`docs/output_format/curation.md`、CLAUDE.md
  冒頭の規模表記、登録数と `ToolAnnotations`（`tests/test_server_registration.py`）。
- vault 側の流れ図にキュレーションの候補付けの枝を足す。
