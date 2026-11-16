# アラインメントのキュレーション（`curation_review` 等）の出力フィールド

`curation_review` / `curation_suggest` / `curation_submit` / `curation_flags` が返す値の意味。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項は核に、参照ライブラリ照合のスコアの意味（`-1`/`0` の区別など）は `library` トピックにある。

呼び出しの順序と実装（`metabolomix/curation/`）は `docs/workflow/curation.md`。

## アラインメントのキュレーション

アダクトはこの文書のどこでも、綴りでなく組成で比べる: n_mol・価数・極性・付加と脱離の原子が同じなら同じイオン（`[M+HCOO]-`・`[M+FA-H]-`・`[M+CHO2]-`、`[M+CH3COO]-`・`[M+Hac-H]-`、`[M-H-H2O]-`・`[M-H2O-H]-`、`[M+H-H2O]+`・`[M-H2O+H]+`、mzTab-M の `[M+H]1+` と `[M+H]+`）。

### `curation_review` の戻り値

| キー | 意味 |
|---|---|
| `review_id` | レビューの識別子。`curation_submit` / `curation_view_data` に渡す |
| `n_spots` | 対象にしたスポット数 |
| `counts` | 判定の件数 `{"ok": N, "suspect": N, "likely_wrong": N}` |
| `warnings` | レビュー全体の注意（下記） |
| `trend` | クラス（`ontology`）ごとの RT–m/z 傾向フィットの要約 `{n, r2, n_outliers}`。`r2` は外れ値を除いた当てはまりで、低いクラスでは `trend_outlier` 系の理由コードを割り引いて読む。`n_outliers` は頑健 z が `trend_outlier_z` を超えた点の数（クラスの信頼度によらず数える） |
| `table` | `suspect` 以上とフラグ済みのスポットの TSV（列は「`table`（TSV）の列」。`ok` かつ未フラグは載らない）。`likely_wrong` → `suspect` → フラグ済みの `ok`、同順位は `spot_id` 昇順で、先頭 `max_rows` 行（既定 100）だけ |
| `n_table_rows_total` | 載せる条件に合うスポットの総数（`max_rows` で切る前） |
| `n_table_rows_shown` | `table` に載せた行数。`n_table_rows_total` より小さければ、残りは `html_path` のビューアにだけある |
| `table_note` | `table` が全件か先頭だけかの 1 行の注記 |
| `html_path` | ビューア HTML のパス（英語表示。`max_rows` で切らず全件）。人がブラウザで開いてフラグを付け、**Submit** で直接送る（受け口が無いときは「Copy submission text」の文をチャットに貼る）。画面は次節の「ビューアの画面」 |
| `submit` | ビューアから直接送れるか。`{"via": "viewer", "idle_timeout_min": 30}` ならビューアの Submit（候補付けは「送信」）が記録と `_tags.xml` の反映を行う。`{"via": "copy", "reason"}` なら受け口を立てられず、送信用テキストを貼る経路だけ |
| `thresholds` | 使ったしきい値（`DEFAULT_THRESHOLDS` を `thresholds` 引数で上書きしたもの） |
| `ms2_tol` | 対向照合の MS2 許容幅（`.dbs` の `search_params`、無ければ既定値） |

`warnings` は次の 4 通り:

1. 照合結果を持つスポットのうち参照を引けた割合が半分未満（アラインメントと別のライブラリを読んでいる可能性）。
2. `.dcl` / `.EIC.aef` / `.arf` の兄弟ファイルが無い（その系統の判定は `UNKNOWN`）。
3. **orphaned**: `flags.jsonl` に、同じファイル名で sha256 が違う（以前の版の）`.arf2` に記録された有効フラグがある。件数を出し、スポット番号の対応が保証されないので当てない。
4. 照合結果はあるのに、MS-DIAL の脂質規則フラグ（`is_lipid_class_match` / `is_lipid_chains_match` / `is_other_lipid_match`）がどのスポットにも立っていない。脂質以外の採点器の出力とみなし、規則に基づく 3 コード（`class_rule_rejected` / `class_rules_not_run` / `chains_unsupported`）を使わない（「判定（`verdict`）と理由コード」の末尾）。

エラー（`{"status": "error", "message": ...}`）で止まるのは、`file_ids` に `.arf` に無い試料 ID がある（`missing_file_ids`。重い読み込みの前に止まる）、`flags.jsonl` に読めない行がある（`flags_file` と `line`。「フラグ記録が壊れているとき」）、対象が上限 3000 スポットを超える、など。

`review-<id>.html` / `review-<id>.json` は、読み込んだ参照ライブラリのスペクトル（対向プロットの参照側）をそのまま埋め込む。研究室のライブラリで作ったものは、ライブラリ本体と同じく研究室の外へ出さない。

### ビューアの証拠 payload（`curation_view_data` / HTML ビューア専用）

`curation_review` の戻り値には出ない（座標は LLM に返さない）。`html_path` のビューアと `curation_view_data`（MCP Apps 専用。LLM は呼ばない）がスポットごとに持つ。

EIC 系列（`eic.samples[].points`）と対向スペクトル（`mirror.measured` / `mirror.reference`）は、判定（`eic_shape` / `rescore` / `matched_mz` など）を全点・全スペクトルで計算した後に、返す座標列だけを間引いたもの（判定は変わらない）:

- `eic.samples[].points`: 1 試料あたり最大 `EIC_MAX_POINTS`（既定 40）点。先頭・末尾・頂点（最大強度）・積分範囲（`left` / `right`）に最も近い点を必ず残し、残りはトレース全体に等間隔に散らした点で埋める（アンカーとの和集合が上限に収まるまで等間隔の点数を減らす）。強度は整数に丸める。
- `mirror.measured` / `mirror.reference`: それぞれ最大 `MIRROR_MAX_PEAKS`（既定 150）点。一致した点（測定側は `matched_measured_mz`、参照側は `matched_mz`）を必ず残し、残りは強度の降順で埋める。`scored_peak_count` / `unscored_peak_count` は間引く前の満スペクトルの値。

スポットごとのほかの項目:

- `auto_note`: 判定根拠の文（`likely_wrong` / `suspect` のとき。`ok` は `null`）。`Auto: ` に続けて理由コードごとの根拠（英語）を ` / ` でつなぐ（例 `Auto: m/z — Δm/z 12.4 mDa (≥10 mDa) / MS2 — no match to the reference (low score)`）。ビューアのメモ欄の既定値（記録済みのメモがあればそちら）。
- `flag_cleared`: このアラインメントで最新のフラグ行が `clear`（人が明示的に取り消した）。
- `confirmed`（真偽値）: 記録の `confirmed` か、レビュー作成時の `_tags.xml` の Confirmed（MS-DIAL の GUI で付けたものも）。最新の記録が `wrong` / `suspect` なら `false`。
- `isotopes`: MS1 の同位体パターン（ビューアの MS1 パネル用。判定には使わない）。
  - `measured`: `.arf2` Key 53 `IsotopicPeaks`（代表試料の M, M+1, M+2）の `[m/z, 相対強度 %（M = 100）, 絶対強度]`。絶対強度は要素の Key 4 `AbsoluteAbundance` を整数に丸めたもの、無ければ `null`。m/z は MS-DIAL が書いた単同位体 + 1.00467·k の計算値で、実測の質量ではない。
  - `theoretical`: `{"relative": [M, M+1, M+2 の %], "basis": ...}`。組成式（`Formula`）にアダクトの原子を足し、整数質量の分解能で畳み込んだもの。`basis` は `formula+adduct`、アダクトが読めなければ中性の組成式だけの `formula`。組成式が無い・読めない、または同位体の表に無い元素（Se・Fe など）を含めば `null`（単一同位体の I・As・Co は同位体のピークを作らない元素として数える）。

**ビューアの画面**（`html_path`。人が操作する）:

- 自動判別が `likely_wrong` のカードは赤の破線枠（未確認）で、人が Wrong を付けると実線になる。
- 開いた時点で、記録済みのフラグが無く、`flag_cleared` でも `confirmed` でもない次のスポットを Wrong の未送信の変更にする: `likely_wrong`、理由に `low_score`（MS-DIAL の名前の `low score:`）、情報に `msms_absent`（`no MS2:` など）。メモは `auto_note`（MS/MS なしは判定理由ではないので `Auto: no MS/MS` を書き足す）。そのまま送信すれば記録され、Correct（無印。何も記録しないことの表示名）に戻して送信すれば `clear` が記録される（次のレビューで掛け直さない）。送信は `_tags.xml` の Misannotation にも反映する（「`curation_submit` の戻り値の `tags_xml`」）。
- カード右上の **Confirmed** チェックは判断 `confirmed`（注釈が正しいと確かめた）。Suspect / Wrong と排他（チェックすると Correct に、Suspect / Wrong を選ぶと外れ、外すと Correct＝送信では `clear`）。初期状態はスポットの `confirmed`。メモが自動の判定根拠のままならチェックで空にする。
- 絞り込み: Class（各クラスの下に「<クラス> › likely_wrong」、先頭に全クラス横断の「likely_wrong (all classes, N)」。該当が 1 件以上のときだけ）、Verdict（自動。`ok` / `suspect` / `likely_wrong` / suspect or worse）、Flag（人。Correct（無印か confirmed）/ Confirmed / Suspect / Wrong / Flagged（Suspect か Wrong））。別に Exclude のチェック（low score = 理由 `low_score` / no MS/MS = 情報 `msms_absent`）が当たるスポットを一覧から隠す。隠しても初期選択の Wrong は残り、送信される。
- クラス別傾向のカード: Plotly（cdnjs から読む。届かなければ canvas 描画）で、横軸 RT・縦軸 m/z（範囲はクラスの全点で固定）。回帰直線はサーバで当てたもの。点に重ねると `#spot_id 名前`・DB・RT・m/z が出る。カードを押すとそのクラスに絞り込み（もう一度押すと解除）、点を押すとそのクラスに絞ってスポットのカードへ移る（点として拾うのは hover が出ているときだけ。それ以外はクラスの切り替え）。点は一覧と同じ絞り込み（Verdict・Flag・Exclude、Class の「› likely_wrong」）に従い、Class の選択そのものは他のカードの点を消さない（カードの強調だけ）。
- MS1 パネル: `isotopes` の棒の上に「実測 % (絶対強度) / 理論 %」を出す（絶対強度は 1000 以上なら有効 2 桁の指数表記）。
- ミラープロット: 縦軸は相対強度（上下とも各側の最大値を 100）で目盛り 0/50/100。ピークの m/z（小数 4 桁）のラベルは `library_plot_mirror` と同じ規則（強度降順の貪欲法・縦横とも重なるものを飛ばす・片側 25 本まで・側ごとに独立）。

注釈を確かめたスポット（記録の `confirmed` ∪ `_tags.xml` の Confirmed − 最新の記録が `wrong` / `suspect`）には、自動判定の `likely_wrong` を当てない: ビューアの初期選択、`exclude_auto_likely_wrong`（`dataset_plot` の `kind="species"` / `"species_pca"` / `"group_intensity"`）、`curation_suggest` の `likely_wrong` の対象。

### `table`（TSV）の列

| 列 | 意味 |
|---|---|
| `spot_id` | `MasterAlignmentID` |
| `name` | ARF2 の `Name`（確度接頭辞を含む生の文字列。接頭辞の規則は `core` トピック §9-1） |
| `ontology` | ARF2 の `Ontology`（脂質クラス） |
| `adduct` | ARF2 の `AdductType` |
| `verdict` | `likely_wrong` / `suspect` / `ok`（次節） |
| `reasons` | 立った理由コードをカンマ区切りで、強い → 弱い → 帯のみの順（次節の表の順）。相手のスポットを指すコードは `:<相手の spot_id>` が付く（`adduct_isomer_of:189`） |
| `ppm` | 代表試料の m/z と参照 precursor m/z（引けなければ Formula からの理論値）の相対誤差 [ppm] |
| `dmz_mda` | `ppm` と同じ基準の差 [mDa] |
| `drt` | 代表試料の RT と参照 RT の差 [分]。参照が引けない・参照に RT が無ければ空。情報 `rt_not_used_by_annotation` のときも値は出るが、判定には使わない |
| `wdot` | MS-DIAL 自身が同定時に出した `squared_weighted_dot_product` の平方根（`library` トピック §14.1）。値が無い、または `-1`（比較不能）なら空（必須注意 4） |
| `mpp` | MS-DIAL 自身が出した `matched_peaks_percentage`（必須注意 4） |
| `eic_good` | 検出（非 gap-fill）試料のうち EIC 形状が「良い」割合（`good_fraction`）。良否は頂点がウィンドウ内・点数が十分・理想ガウスとの R² が十分・極大の数が上限以下で決まる。検出試料が無ければ空 |
| `trend_z` | 所属クラスの RT–m/z 傾向モデル（Huber 頑健回帰。RT = a + b·炭素数(+ c·二重結合数)）に対する頑健 z。傾向を当てはめられないクラス（点数不足・炭素数のばらつき不足）では空 |
| `flag` | 既存のフラグ（`wrong` / `suspect`）。未フラグは空。有効フラグ（スポットごとの最新 1 行）が `clear` のスポットは未フラグとして扱うので、`clear` はここに出ない |

### 判定（`verdict`）と理由コード

`verdict` は 6 系統の帯（`checks.{msms,mz,rt,eic,ion,trend}`。各 `PASS` / `BORDERLINE` / `FAIL` / `UNKNOWN` と理由コード）から決める。`ion` は別アダクトの取り違え（下記）で、`FAIL` か `PASS`（`adduct_isomer` を評価していないときだけ `UNKNOWN`）。

- **`likely_wrong`**: 強い理由コードが 1 つでも立った。
- **`suspect`**: 強い理由が無く、(a) msms/mz/rt/eic/ion のどれかが `FAIL`、(b) msms/mz/rt/eic/ion のうち `BORDERLINE` が 2 つ以上、(c) `BORDERLINE` が 1 つ以上あり RT–m/z 傾向も `BORDERLINE`（`trend_outlier`）、のどれか。
- **`ok`**: それ以外。

理由コードの区分: **強い**（`STRONG_REASONS`）は単独で `likely_wrong` にする。**弱い**（`WEAK_REASONS`）は当該系統を `FAIL` にするが、単独では `suspect` 止まり。**帯のみ**は `BORDERLINE` の内訳で、`suspect` の積み上げにだけ数える。**情報**は `verdict` を動かさない注記（`info`。`checks` とは別枠で常に出る）。

| 区分 | コード | 意味 |
|---|---|---|
| 強い | `polarity_mismatch` | アダクトの電荷符号と実測のイオン化極性（`IonMode`）が不一致（`adduct_consistency` の `band` が `FAIL`）。符号はアダクトの表記を組成で読んで決める（MS-DIAL の `[M+K]+` なども読む。読めない表記だけ `UNKNOWN`）。クラスに典型的なアダクトか（`class_typical`）は advisory で `band` に効かず、非典型のアダクトだけでは立たない |
| 強い | `precursor_unmatched` | MS-DIAL 自身の `is_precursor_mz_match` が `False` |
| 強い | `dmz_out` | \|Δm/z\|（`dmz_mda`）が `dmz_fail_mda`（既定 10 mDa）以上。ppm は m/z に比例して緩むので絶対差で見る（Δppm の `ppm_out` は弱い） |
| 強い | `class_rule_rejected` | MS/MS ありで、MS-DIAL の脂質クラス規則（診断イオン）が棄却した（`is_lipid_class_match=False` かつ `is_other_lipid_match=False`）。実測では全件 `low score:` なので `low_score` も立つ。脂質規則が走ったデータだけ（この節の末尾） |
| 強い | `adduct_isomer_of:<spot>` | m/z が、同時に溶出する別物質のスポット `<spot>` の中性質量の別のアダクトで説明でき、`<spot>` の証拠の段階が上（下記「別アダクトの取り違え」）。単独の証拠では `ok` のスポットにも立つ |
| 弱い | `ppm_out` | Δppm が `ppm_borderline`（既定 10）を超えた。アダクトによらない。mz 系統は `FAIL` になるが、単独では `suspect` 止まり |
| 弱い | `low_score` | MS/MS はあるが MS-DIAL 自身の `is_reference_matched` が `False` |
| 弱い | `drt_out` | ΔRT が `drt_borderline`（既定 1.0 分）を超えた |
| 弱い | `eic_poor` | EIC 形状帯が `FAIL`（検出試料が 1 件以上で、`good_fraction` が `eic_borderline_frac` 未満）。検出試料が 0 件なら帯は `UNKNOWN` で、このコードは立たない |
| 弱い | `adduct_isomer_minor_of:<spot>` | `adduct_isomer_of` と同じ条件で説明できるが、`<spot>` の証拠の段階は同じで、強度（`HeightAverage`）が `adduct_isomer_height_ratio`（既定 3）倍以上 |
| 帯のみ | `ppm_borderline` | Δppm が `ppm_pass`〜`ppm_borderline`（既定 5〜10） |
| 帯のみ | `drt_borderline` | ΔRT が `drt_pass`〜`drt_borderline`（既定 0.5〜1.0 分） |
| 帯のみ | `eic_borderline` | EIC 形状帯が `BORDERLINE`（検出試料が 1 件以上で、`good_fraction` が `eic_borderline_frac`〜`eic_pass_frac`）。検出 0 件では立たない |
| 帯のみ | `rt_scatter` | 検出試料間の EIC 頂点 RT の母標準偏差（検出 2 件以上のときだけ計算）が `eic_rt_scatter_sd`（既定 0.1 分）を超えた。EIC 帯が `PASS` なら、これだけで `BORDERLINE` に下げる |
| 帯のみ | `trend_outlier` | 所属クラスの RT–m/z 傾向で頑健 z が `trend_outlier_z`（既定 3.0）を超え、そのクラスの傾向フィットが信頼できる（`r2 >= trend_min_r2`、既定 0.7） |
| 情報 | `msms_absent` | MS/MS 未取得、または名前接頭辞が `no MS2` / `w/o MS2`（`msms` 系統は `UNKNOWN`） |
| 情報 | `reference_not_found` | ライブラリから参照レコードを引けなかった。`rt` 系統は `UNKNOWN`。`mz` 系統は Formula / AdductType からの理論値（`ppm_basis="formula"`）で計算を続け、組成式かアダクトが読めないときだけ `UNKNOWN` |
| 情報 | `reference_rt_absent` | 参照は引けたが RT を持たない（`rt` 系統は `UNKNOWN`） |
| 情報 | `rt_not_used_by_annotation` | MS-DIAL がこの照合で RT を絞り込みにも採点にも使っていない（`.dbs` の注釈器の `use_time_for_annotation_filtering` / `use_time_for_annotation_scoring` がどちらも False）。参照 RT は別のクロマトグラフィーの予測値でありうる（LBM は約 18 分系の予測 RT）ので、`rt` 系統を `UNKNOWN` にして `drt_out` / `drt_borderline` を立てない（`drt` の値は残る）。スイッチが分からない store（`.msp` / `.lbm2` を直接読んだもの）では RT で判定する |
| 情報 | `no_match_result` | ARF2 に MS-DIAL の照合結果（`representative`）が無い（`msms` 系統は `UNKNOWN`） |
| 情報 | `class_rules_not_run` | 脂質クラス規則が評価されていない（`is_lipid_class_match=False` かつ `is_other_lipid_match=True`。CompoundClass が Unknown / Others・SPLASH・名前の解析失敗など）。参照一致でも規則の裏付けが無い（未検証）という意味で、誤りの意味ではない。脂質規則が走ったデータだけ |
| 情報 | `chains_unsupported` | 名前が鎖レベル（`16:0_18:1`・`18:1;O2/16:0` のように鎖を `_` / `/` で区切る。`PC 34:1\|PC 16:0_18:1` の `\|` 以降も含む。単鎖は数えない）なのに `is_lipid_chains_match=False`。MS-DIAL は照合に失敗しても参照名から `\|` 付きの名前を作るので、鎖組成の裏付けは `is_lipid_chains_match=True` だけ。MS/MS ありのときだけ。脂質規則が走ったデータだけ |
| 情報 | `rescore_discrepancy` | `curation_review` が対向照合で出した `weighted_dot_product` が、MS-DIAL 自身の値（平方根換算）と `rescore_tolerance`（既定 0.1）を超えてずれた |
| 情報 | `adduct_differs_from_reference` | 参照レコードのアダクトと注釈のアダクトが別のイオン（組成で比べるので、同じイオンの書き方の違いでは立たない） |
| 情報 | `manually_modified` | MS-DIAL 側で `is_manually_modified` が立っている（GUI で人手修正済み） |
| 情報 | `manually_unsettled` | 名前接頭辞が `unsettled`（MS-DIAL の Alignment Viewer で未確定のまま） |
| 情報 | `trend_outlier_unreliable` | 頑健 z はしきい値を超えたが、そのクラスの傾向フィットが信頼できない（`trend_outlier` に数えない） |
| 情報 | `dcl_precursor_mismatch` | `.dcl` の該当インデックスの precursor m/z が代表試料の m/z と 0.01 以上ずれる（順序対応が崩れている可能性があるので、MS/MS 系統の実測スペクトルは使わない） |

**`likely_wrong` はまれ。** MS-DIAL は同定時に precursor の許容幅を課している（`precursor_unmatched` は MS-DIAL 自身の `IsPrecursorMzMatch`）ので、許容幅の外の候補は代表に残りにくい。`likely_wrong` が 0 件でも「全部正しい」ではなく、`suspect` の理由コードを読む。この絞り込みの後にも残る強い理由は `class_rule_rejected`。

**別アダクトの取り違え（`adduct_isomer_of` / `adduct_isomer_minor_of`）。** 別の物質の別アダクトを近い質量の別物質として注釈した誤り（個々のスポットの Δppm・MS2 では見えない）を、同時に溶出するスポットどうしで見る。対象 X と相手 Y について:

- **相手の範囲**: レビューの対象（`ontology` / `name_contains`）に絞らず、アラインメントの注釈付きスポット全部。有効な `wrong` フラグのあるスポットは相手にしない。
- **条件**: Y が X と同じ極性で |ΔRT| ≤ `adduct_isomer_drt`（既定 0.05 分）。X と Y は別物質（組成式とクラス（`Ontology`）が両方同じなら、同じ分子種の別アダクト＝正しい注釈なので対象外）。Y の m/z と注釈アダクトから中性質量を出し、仮定アダクト A での m/z が X の m/z と `adduct_isomer_ppm`（既定 5 ppm）以内。A は X の極性の既定アダクト一覧（`msdial.analysis_params.DEFAULT_ADDUCTS`）と X 自身のアダクトから、Y 自身のアダクト（組成で見分ける）と、Y と同じ m/z になるアダクトを除いたもの（どちらも X が Y と同じ m/z にいるだけで、重複スポットの話）。m/z・RT はどちらもスポットの `MassCenter` / `RT`。
- **証拠の段階**: 3 = MS2 で鎖組成まで裏付け（鎖レベルの名前 ∧ `is_lipid_chains_match`）、2 = MS2 が参照と一致（`is_reference_matched`）、1 = MS2 はある（low score）、0 = MS2 無し（名前接頭辞 `no MS2` / `w/o MS2` を含む）・照合結果無し。Y の段階が上なら `adduct_isomer_of`（強い）、同じ段階で強度比（Y / X）が `adduct_isomer_height_ratio` 以上なら `adduct_isomer_minor_of`（弱い）。Y の段階が下なら立たない（逆向きの組は Y の側に立つ）。
- 複数の Y が説明するときは、強い → 段階 → 強度の順で 1 件だけ。

レビュー JSON（とビューア）のスポットの `adduct_isomer`: `of`・`of_name`・`of_adduct`（相手とその注釈アダクト）・`as_adduct`（X を説明する仮定アダクト A）・`ppm`（X の m/z と A の理論 m/z の差）・`drt`（X − Y の RT [分]）・`severity`（`strong` / `minor`）・`tier` / `of_tier`・`height_ratio`（Y / X）。説明する相手が無ければ `null`。

**脂質規則フラグは、脂質規則が走ったデータでだけ読む。** 規則フラグは MS-DIAL の Lipidomics 採点器でしか立たず、ほかの採点器では全部 `False`（`class_rule_rejected` の条件と区別できない）。そこで、1 回のレビューの対象スポットに規則フラグ（class / chains / other）が `True` のものが 1 件でもあるときだけ上の 3 コードを使い、無ければ `warnings` に出す。

### `curation_suggest` の戻り値

間違いになったスポットと未注釈スポットに、注釈の候補（①②）と、別スポットのイオンとしての説明（④）を並べる。先に `library_load` と `curation_review`（無ければ `missing_state`）。対象のデータセットは `dataset_handle`（省けば直近）。

| キー | 意味 |
|---|---|
| `suggestion_id` | 候補付けの識別子（`cs-YYYYMMDD-HHMMSS-xxxx`）。`curation_submit` の送信用テキストが運ぶ |
| `base_review_id` | `likely_wrong` の判定を読んだ元のレビュー。`review_id` を省くと、このアラインメント（sha256 一致）の最新のレビュー |
| `n_spots` | 対象のスポット数 |
| `counts` | `targets`（対象の種類別 `{flagged, likely_wrong, unannotated}`）・`with_candidates`（候補か、情報でないイオン関係があるスポット数）・`with_strong_relation`（強い説明のあるスポット数）・`hard_removed`（ハード制約で削った候補の延べ数） |
| `warnings` | param ファイルが無く既定の検索アダクトと RT 窓を使った、`PeakProperties.arf` が無く相関を計算できない、兄弟ファイルが無い、など。「`likely_wrong` は元レビューの対象範囲についてしか分からない」は常に出る（必須注意 9） |
| `table` | 1 スポット 1 行の TSV（列は下表）。強い説明のあるスポット → 対象の種類（`flagged` → `likely_wrong` → `unannotated`）→ `spot_id` 昇順で、先頭 `max_rows` 行（既定 100）だけ。載せた行数は `n_table_rows_shown`。全件は `html_path` |
| `html_path` | 候補付けビューア HTML のパス。人がブラウザで開いて選び、「送信」で直接送る（受け口が無いときは「送信用テキストをコピー」の文をチャットに貼る） |
| `submit` | `curation_review` の `submit` と同じ |
| `library` | 使った参照ライブラリの `file`（ファイル名だけ）と `sha256` |
| `analysis_params` | 検索アダクトと RT 窓の出所。`source`（param ファイルか既定か）・`path`・`rt_window`（分） |
| `options` / `thresholds` | 使った引数（`rt_window` は解決後の値）としきい値 |

`table`（TSV）の列:

| 列 | 意味 |
|---|---|
| `spot_id` | `MasterAlignmentID` |
| `name` | 今の名前（未注釈は `Unknown` 系） |
| `target` | 対象の種類。`flagged`＝有効な `wrong` フラグのあるスポット、`likely_wrong`＝元レビューの判定が `likely_wrong`（フラグ無し）、`unannotated`＝`.arf2` に代表が無いスポット。`include_decided=True` なら、判断済み（assign / redundant）の注釈付きスポットも `flagged` で戻る |
| `top_candidate` | 候補（①②）の先頭の名前（和組成があれば和組成）。候補が無ければ空 |
| `source` | 先頭の候補の出所。`msdial`（MS-DIAL 自身の下位候補）・`research`（閾値を緩めた再検索だけで出たもの。前駆体の窓は `library_match_feature` と同じ `ms1_tolerance`（既定 0.01 Da。m/z 500 を超えたら同じ ppm に広げる）で、窓に入っても `dmz_out` はハードで落とす）・`msdial+research`（両方） |
| `total_score` | 先頭の候補の総合スコア（`library` トピックの `total_score` と同じ物差し。RT・precursor の一致度を足した非正規化の和で 1 を超える）。比較できなければ `-1` |
| `reasons` | 先頭の候補のソフト理由・情報の理由コード（カンマ区切り） |
| `relation` | 先頭のイオン関係（④）を `<関係コード>(#<相手スポット>)` で。無ければ空 |
| `strong` | 強い説明があるか（`True` / `False`） |

### 候補付けの payload（HTML ビューア専用）

`curation_suggest` の戻り値には出ない（座標は LLM に返さない）。`suggest-<id>.json`（正本）と `suggest-<id>.html` のビューアだけが持つ。スポットごとに `eic`（対象の EIC）・`measured`（測定 MS/MS）・`candidates`（①②）・`relations`（④）・`strong`・`preset`（強い説明があるときの初期選択の `R<n>`）。`curation_submit` はここから候補の詳細を引いて記録行に展開する。

**候補（`candidates[]`、①②）**: `candidate_id`（`L1`…、並び順）・`source`・`name`（分子種）・`sum_name`（和組成。解析できなければ `null`）・`ontology`・`adduct`・`formula`・`inchikey`・`library_id` / `record_index`（参照レコード）・`dmz_mda` / `ppm`（代表試料の m/z と参照 precursor m/z の差）・`scores`（`total_score` `weighted_dot_product` `simple_dot_product` `reverse_dot_product` `matched_peaks_percentage` `matched_peaks_count`）・`trend`（RT–m/z 傾向の予測 RT・残差・z。当てはまらなければ `null`）・理由コード（`soft` / `info`。ハードは削られるので無い）・`mirror`（対向プロットの参照側）。`total_score` などの `-1` は比較していない（スポットに MS/MS が無い、または参照が引けずスペクトルが無い。`matched_peaks_count` は 0）で、ビューアは `–` と表示する（必須注意 8）。

**イオン関係（`relations[]`、④）**: 次々節。ビューアは関係の行を選ぶ（または指す）と、相手 Y の代表試料の EIC（`partner_eic`）を X の代表試料の EIC に重ね、X の代表試料の最大値に合わせて縮めた点線で描く（必須注意 6）。

### 候補の理由コード

候補（①②）ごとに次の制約を当てる。**ハード**は候補から削る（`hard_removed` に数だけ残る）、**ソフト**は順位を下げる、**情報**は何もしない。並び順は（ソフト理由の数 昇順、`total_score` 降順）で、スポットごとに上位 `top_n`（既定 5）件。

| 区分 | コード | 意味 |
|---|---|---|
| ハード | `polarity_mismatch` | 候補のアダクトの電荷符号とスポットのイオン化極性が矛盾 |
| ハード | `dmz_out` | \|Δm/z\|（代表試料の m/z と参照 precursor m/z の差）が `dmz_fail_mda`（既定 10 mDa）以上 |
| ソフト | `adduct_atypical` | そのクラスに典型的でないアダクト（組成で比べる: PC の `[M+FA-H]-` は典型の `[M+HCOO]-` と同じ） |
| ソフト | `trend_outlier` | 候補のクラス×和組成から予測した RT との残差の \|z\| が `trend_outlier_z`（既定 3.0）を超えた |
| ソフト | `no_matched_peaks` | MS/MS はあるが、比べた結果の一致ピークが 0（参照が引けず比べていない候補は `reference_unresolved`） |
| ソフト | `msms_absent` | スポットに MS/MS が無い（①の候補にだけ付く。②は MS/MS が無いと働かない） |
| ソフト | `no_support`（イオン関係） | ④の関係に裏付けが無い（次節） |
| 情報 | `trend_unknown` | 候補のクラスの傾向が点数不足などで当てはまらない。順位は下げない（「傾向から外れた」ではなく「傾向を持たない」） |
| 情報 | `reference_unresolved` | ①の候補の参照レコードがライブラリに無い（スペクトルを採点し直せず、`scores` は `-1`） |

傾向モデルは、そのアラインメントの注釈付きスポットから、対象・`wrong` フラグ付き・`redundant` 判断済みを除き（`assign` 済みは記録した名前・クラスで）、クラスごとに `RT = a + b·炭素数 + c·二重結合数`（Huber）で当てはめる。元レビューの範囲には縛られない。

### イオン関係（`relation`）

対象 X を、注釈付きの別スポット Y のイオンとして説明できるとき、`relations[]` の行（`candidate_id` は `R1`…）になる。①② とは物差しが違うので、候補と 1 本の順位には混ぜない。

**相手 Y**: 注釈付きで、有効な `wrong` フラグが無く、今回の対象でも `redundant` 判断済みでもないスポット。Y が `assign` 済みなら、記録した名前・組成式・アダクトで扱う（同位体の期待比も記録した組成式で計算する）。

**関係になる条件**（両方）: (1) \|ΔRT\| ≤ `rt_window`（既定は param ファイルの `Retention time tolerance for alignment`、無ければ 0.1 分）、(2) Δm/z が下の関係のどれかで `relation_mz_tol`（既定 10 mDa）以内、または MS-DIAL の `found_in_upper_msms` リンクがある。

| 関係コード | 意味 |
|---|---|
| `adduct:<X のアダクト>/<Y のアダクト>` | 同じ中性分子の別アダクト。Y の中性質量に、解析で検索したアダクト（param ファイルの `Searched adduct ions`。二量体・多価を含む。無ければ極性ごとの既定）を付けた m/z が X に一致。Y 自身のアダクト（組成で見分ける）と、Y と同じ m/z になるアダクト（`[M+H]+` に対する `[2M+2H]2+`。Δm/z が区別しない）は除く |
| `isotope_M+1` / `isotope_M+2` | Y の M+1（+1.003355）/ M+2 の 13C2（+2.006710）。どちらも電荷で割る。18O・41K を含むことがある（表の下） |
| `isotope_M+2_18O` | Y の 18O 同位体（+2.004245 / 電荷）。Y のイオンが O を持たなければ出さない |
| `isotope_M+2_41K` | Y の 41K 同位体（+1.998119 / 電荷）。Y のイオンが K を持たなければ出さない（`[M+K]+` は K を持つ） |
| `isotope_M+2_34S` / `isotope_M+2_37Cl` / `isotope_M+2_81Br` | Y の 34S / 37Cl / 81Br 同位体（+1.995796 / +1.997050 / +1.997952 を電荷で割る）。Y のイオンの組成（組成式の n_mol 倍にアダクトの原子を足したもの。`[M+Cl]-` は Cl を持つ）が分かっていてその元素を持たなければ出さない（組成式が無い・読めないなら出す。同位体の表に無い元素（Se など）を含んでも、組成式が読めればこの門には使う） |
| `insource:-H2O` `-2H2O` `-NH3` `-HCOOCH3` `-CH3COOCH3` `-C3H5NO2` `-C2H8NO4P` `-C3H8NO6P` `-C6H10O5` | Y からの中性損失によるインソース断片（質量はそれぞれの組成式から計算） |
| `found_in_upper_msms` | 質量差が上のどれでも説明できないが、MS-DIAL が「Y の MS/MS に X が見える」とリンクしている |

**18O・41K と 13C2**: `relation_mz_tol` が 13C2 との差（18O は 2.46 mDa / 電荷、41K は 8.59 mDa / 電荷）より狭ければ、18O・41K の行を独立に見る。差以上の窓（既定の 10 mDa はどちらも）では、X が 13C2 の窓に入るかぎり `isotope_M+2` に含める（`dmz_mda` は 13C2 の位置からの差なので、41K のピークなら約 −8.6 mDa）。X が 13C2 の窓の外なら、X が当たった M+2 の位置のうち 18O（41K）がいちばん近いときだけ `isotope_M+2_18O`（`_41K`）を出す（分離された 34S のように、より近い位置で説明される X には出さない）。41K で「窓の外」は、41K の位置から 13C2 と逆の側へ「窓 − 8.59 mDa」（既定の窓では 1.41 mDa）より離れた X。

**`found_in_upper_msms` リンクの向き**: MS-DIAL は向きの無い対で保存するので、前駆体（upper）側は m/z の大小で決める。リンクだけで関係にするのも、リンクで強い説明にするのも、X の m/z が Y の m/z より小さい（X が断片）ときだけ。裏付け（下）には向きによらず数える。

**裏付け**: MS-DIAL のリンク（`correl_similar` / `chrom_similar` / `found_in_upper_msms`）があるか、`.arf` の試料別の高さの log1p で求めた Pearson の `r` ≥ `relation_min_r`（既定 0.8。両方の高さが正の試料が 5 以上のときだけ計算し、足りなければ `r` は空）。どちらも無い関係は載せるが `no_support`（ソフト）を付けて順位を下げる。`msdial_links` は X–Y 間の MS-DIAL のリンク種別。

**強い説明**（`strong = true`）: 情報でなく、次のどれか。(a) 裏付けがあり、`found_in_upper_msms` の関係か、Y が前駆体側（上の m/z の規則）の `found_in_upper_msms` リンクを持つ関係。(b) 同位体で、裏付けがあり、`isotope_consistent = true` で、`ambiguous` でない。(c) アダクトの組・インソース断片で、Δm/z と同時溶出（`rt_window`）に加えて試料間の相関 `r` ≥ `relation_min_r`（MS-DIAL の `correl_similar` / `chrom_similar` リンクだけでは足りない。r を測れないときは (a) だけ）。強い説明のあるスポットでは、それがカードの先頭で初期選択（`preset`）になる。

**同位体の強度比**: `isotope_ratio` は試料ごとの X の高さ / Y の高さの中央値（両方が正の試料だけ）。`expected_ratio` は、Y のイオンの組成（中性の組成式の n_mol 倍にアダクトの原子を足したもの。`[2M+H]+` は炭素も 2 倍、`[M+Cl]-` は Cl を 1 個、`[M+HCOO]-` はギ酸の原子を持つ。アダクトが読めない・差し引いて負になるなら中性の組成式のまま）と元素ごとの天然存在比（IUPAC の代表値。H・C・N・O・P・S・Na・K・Cl・F・Br と単一同位体の I・As・Co）から求めた、M+k の同位体の組み合わせ（M+1 は 13C・15N・2H・17O・33S・40K、M+2 は 13C2・18O・34S・37Cl・81Br・41K と 13C+15N などの交差項。質量差は NIST の原子質量から）のうち、m/z（Y + 質量差 / 電荷）が X から `relation_mz_tol` 以内に来るものの、M に対する比の和。窓の外の同位体は足さない: 分離された 1 つの位置（Cl を持つイオンの 13C2 だけ、スルファチドの 34S だけ）はその近くの同位体とだけ比べ、窓が M+k の位置を全部覆えば（既定の 10 mDa で X が 13C の位置にいる M+1 など）M+k の和になる。同じ X と Y の M+2 の行は窓が同じ（X の Δm/z で決まる）なので、同じ期待比を持つ。

- `isotope_ratio` ≤ 1.5 × `expected_ratio` なら `isotope_consistent = true`（Y の同位体として矛盾しない）。
- 1.5 倍を超えれば、X は実在の別物質（同位体の重なり）とみなして `informational = true`（添えるだけで、`redundant` に選べない）。重なった M+2 の行はそろって `informational` になるかならないか。
- 比が求められない（組成式が無い・読めない、同位体の表に無い元素（Se・Fe など、最も多い同位体の下に軽い同位体がある元素）を含む、共通の試料が無い）ときは `isotope_consistent = null` で、強い説明にならない。表に無い元素を含んでも、組成式が読めれば元素の門（Cl・Br・S・O を持たなければその行を出さない）には使う。

**M+2 の取り違え（`ambiguous`）**: M+2 の位置（13C2・18O・34S・37Cl・81Br・41K）は互いに 11 mDa 以内（13C2–37Cl 9.66、34S–37Cl 1.25、13C2–34S 10.91、41K–37Cl 1.07、41K–81Br 0.17 mDa）なので、`relation_mz_tol` の窓に同じ Y の M+2 が 2 つ以上入ることがある。

- 分離できない M+2 の束（13C2・34S・18O が重なった重心。タウロコール酸 [M-H]- で約 +2.0015）は、窓が重心から 13C2 と 34S の両方までの差（重心が真ん中なら 10.91 mDa の半分の約 5.5 mDa / 電荷。既定の 10 mDa は覆う）以上なら両方の窓に入る。そのときは一方に決めず、どの行にも `ambiguous = true` を付けて並べ、強い説明にしない。期待比が束全体の和になるのも、X から束のどの同位体までも `relation_mz_tol` 以内のときだけ。
- それより狭い窓では、近い方の位置の行（`ambiguous` でない 1 行）だけに当たり、期待比は遠い方の同位体を欠くので、Y 自身の束が `informational` になりうる（タウロコール酸 [M-H]- の束で窓 5.5 mDa: 期待比 0.056、束の強度比 0.103）。重心から近い方の位置までの差（タウロコール酸で約 5.2 mDa）より狭ければ 13C2 と 34S の行は出ず、窓が重心から 18O の位置までの差（約 2.7 mDa）以上なら 18O の行だけが立ち、期待比が 13C2 と 34S を欠くので `informational` になる。束を分離しない装置のデータでは `relation_mz_tol` を約 6 mDa / 電荷以上にする。
- Y のイオンが S・Cl・Br・K を持たなければその行は出ないので 13C2 だけになり、`ambiguous` にならない（18O・41K は「18O・41K と 13C2」のとおり、13C2 の行に含むか、いちばん近いときだけ出る）。
- `curation_suggest` のビューアは `ambiguous` の行の見出しに「（M+2 の候補が重なる）」と添える（どれを `redundant` に選んでも、どの同位体かは決まっていない）。

関係の並び順は 強い説明 → 情報でないもの → 裏付けのあるもの → `|Δm/z|` 昇順。上位 3 件の関係は相手 Y の代表試料の EIC（`partner_eic`）を持つ。

### 判断の記録（assign / redundant）

`flags.jsonl`（追記専用。「フラグ」は判断の記録の意味）の `flag` は `wrong` / `suspect` / `confirmed` / `clear` / `assign` / `redundant`。キーはアラインメントの sha256 と `MasterAlignmentID` で、スポットごとに最新の 1 行が勝ち（`wrong` の後の `assign` は `assign`）、`clear` は取り消し。

| flag | 主な項目 |
|---|---|
| `assign` | `name`（記録名）・`level`（`sum` / `species`）・`species_name`（レコード本来の分子種名）・`ontology`・`adduct`・`formula`・`inchikey`・`candidate_source`・`library`（`sha256` `library_id` `record_index`）・`scores`（採点のスナップショット）・`suggestion_id` |
| `redundant` | `of`（Y の `spot_id`）・`relation`（関係コード）・`evidence`（`dmz_mda` `drt` `r` `msdial_links`）・`suggestion_id` |

- 候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign|redundant|clear, candidate, level, note}]`。`candidate` は `assign` なら `L<n>`、`redundant` なら `R<n>`（情報のイオン関係は `redundant` に選べない）。名前・InChIKey は送信側で書かず、保存済みの `suggest-<id>.json` から展開する。`level` の既定は `sum`（`name` = `sum_name`、解析できなければ分子種名）で、`species` なら `name` はレコードの分子種名。和組成で記録しても選んだレコードの InChIKey を持たせる（MS-DIAL が和組成の名前にもライブラリの構造の InChIKey を付けるのと同じ）。この InChIKey は MS/MS で鎖組成まで確かめたことを意味しない。
- `curation_submit` の戻り値は `n_assign` / `n_redundant` / `n_confirmed`（有効な判断の件数）を持つ。`_tags.xml` への反映は「`curation_submit` の戻り値の `tags_xml`」。
- 候補付けビューアの「元の注釈に戻す」は `clear` を送る。`clear` はそのスポットの有効な判断を丸ごと（`wrong` でも `assign` でも）消す。
- `curation_flags` の TSV は 7 列: `spot_id`・`flag`・`name`（`assign` の記録名）・`of`（`redundant` の相手）・`note`・`source`・`ts`。`source` は判断の出どころで、`user`（人の判断。ビューアの送信はいつもこれ）/ `llm`（LLM の判断）。`curation_submit` の呼び出し側が申告し、1 回の送信の全行に同じ値が付く。
- `dataset_feature_search` は判断を反映しない（判断の列も無く、SME / SML の生の同定だけを見る）。判断の記録があれば、件数と「反映していない」を `notes` に 1 行書く（`dataset_export_differential` の同定とは食い違いうる）。
- `curation_review` のスポットの `flag` / `flag_note` は `wrong` / `suspect` / `confirmed` だけ。`assign` / `redundant` は別項目 `decision`（`{flag, name}` または `{flag, of}`）に出る。
- エクスポートでは `assign` は同定を置き換え、`redundant` は除外する（「エクスポートでの `assign` / `redundant` の扱い」）。

### `curation_submit` の `review_id` を省いたとき

送信用テキストも `review_id` も無い送信（`flags` を直接渡す）は、データセット（`dataset_handle`、省けば直近）の `.arf2` に対する最新のレビュー（`cr-…`。`curation_suggest` の既定と同じ規則で、同じ秒の ID は `created_at` の新しいほう）に記録し、戻り値の `review_id` と `review_id_source: "latest"` で示す。レビューが無ければ `missing_state("curation_review")`。候補付け（`assign` / `redundant`）の送信は `review_id`（`cs-…`）を省けない。最新のレビューは対象のスポットを問わず選ぶので、`name_contains` や `ontology` で絞ったレビューが最新だと、その外のスポットは「このレビューの対象外」で断られる（その `status: "error"` にも `review_id` と `review_id_source: "latest"` が付く）。別のレビューに送るなら `review_id` を渡す。

### `curation_submit` の戻り値の `open_handles`（開いているデータセットへの反映）

記録の後、同じ `.arf2`（sha256 で照合）を開いているハンドルすべての `curation_flags`（データセットが持つ有効な判断の写し）を読み直す。書き出しの `apply_curation` は出したばかりの `wrong` / `redundant` をすぐ使う。`refreshed` は読み直したハンドル ID、`not_refreshed` は別の呼び出しがロックを持っていて数秒待っても読み直せなかったもの（その呼び出しが終わってから `dataset_close` → `load_dataset` で開き直すと新しい判断を読む）。同じ `.arf2` を開いているハンドルが無ければこのキーは出ない。

### `curation_submit` の戻り値の `tags_xml`（MS-DIAL への反映）

記録（`flags.jsonl`、正本）の後、アラインメントの `<.arf2 の stem>_tags.xml` の **Misannotation**（タグ Id 3）と **Confirmed**（タグ Id 1）に反映する:

| 判断 | Confirmed | Misannotation |
|---|---|---|
| `confirmed` | 付ける | 外す |
| `wrong` | 外す | 付ける |
| `suspect` | 外す | 触らない |
| `clear` | 外す | 外す |

`assign` / `redundant` は触らない（MS-DIAL の中の名前はまだ誤ったまま）。候補付けビューアの「元の注釈に戻す」（`clear`）は、有効な判断が `wrong` / `assign` / `redundant` のどれでも両方のタグを外す。2 つのタグは 1 回の書き込みで変える。ほかのタグと定義は残し、ファイルが無ければ MS-DIAL と同じ 5 定義で作る。形式は MS-DIAL と同じ（`<Peak Id="<MasterAlignmentID>"><Tag>3</Tag></Peak>`、UTF-8 BOM・CRLF）。

| フィールド | 意味 |
|---|---|
| `path` | 書いた（書こうとした）`_tags.xml` |
| `added` / `removed` | 実際に Misannotation が付いた／外れた `spot_id`（元から同じ状態のものは含まない） |
| `confirmed` | `{"added": [...], "removed": [...]}`。実際に Confirmed が付いた／外れた `spot_id` |
| `created` | ファイルを新しく作ったか |
| `backup` | 書く前の控え（`curation/tags-backup/<名前>.<UTC 時刻>`）。元のファイルが無ければ `null` |
| `note` | MS-DIAL でプロジェクトを開いたままだと GUI の保存で上書きされる旨 |
| `error` / `message` | 反映に失敗したとき（壊れた XML・書けないなど）。フラグは記録済みで、`_tags.xml` は触っていない |

MS-DIAL はアラインメントを保存するたびにメモリ上のタグで `_tags.xml` を丸ごと書き直し、読むのはプロジェクトを開くときだけ。反映を見るにはプロジェクトを閉じて開き直す。

### ビューアからの直接送信（2026-10-09）

`curation_review` / `curation_suggest` は `127.0.0.1` の空きポートに受け口を立て（プロセスで 1 つを共有）、レビューごとのトークンを HTML にだけ埋め込む（レビュー・候補付けの JSON とツールの戻り値には入れない）。ビューアの Submit は `curation_submit` と同じ検証・記録・`_tags.xml` 反映を行い、戻り値も同じ形。`source` は常に `user`。

- 受け口は Finish か、最後の通信（ping は 5 分おき）から 30 分で登録を外し、登録が 0 件になれば止まる。MCP サーバが終われば消える。ハンドラのソケットは 30 秒で時間切れ。
- 応答の状態コード: 401（トークン不明・期限切れ → ビューアは送信用テキストの欄を開く）/ 403（`review_id` 不一致・`Host` 不正）/ 400（検証エラー・JSON でない・`Content-Length` が無い／不正（負を含む））/ 404（不明なパス）/ 405（GET）/ 409（レビューの後でアラインメントが変わった）/ 410（保存済みレビューが無い）/ 413（本文が 8 MB 超）/ 500（`flags.jsonl` が読めない・想定外の例外）。本文は `{status: "error", message, ...}`。
- 送信中は Submit を無効にする。送信中に加えた変更は未送信のまま残り、送った時点の写しだけが反映・消去される。受け口を失った理由（時間切れ・不達など）は結果の行に出て、送信用テキストの欄が開く。
- 送った結果は `curation_flags` で確かめられる。
- 受け口が記録した後で応答が届かなかった場合（ビューアは Copy に切り替わる）、その Copy の文を貼ると同じ行がもう一度記録されるが、スポットごとに最新の 1 行が効き値も同じなので結果は変わらない。

### エクスポートのメタ行（`# curation = ...`）

`dataset_export_differential` は既定（`apply_curation=True`）で有効な判断（`curation_flags` と同じ: スポットごとの最新 1 行、`clear` 済みを除く）を読み、契約の 15 列（`export_contract.EXPORT_COLUMNS`）は変えずに、メタ行の `source_lines` の末尾（`# source_arf` / `# source_mztab` など経路ごとのメタ行の直後）に 1 行足す。有効な判断が 0 件ならこの行を出さず、出力は判断が無いときと同じ。

タブ区切りの 1 行で、外側の角括弧内は `state` が `applied` のときだけ、内側の角括弧内はさらに有効な `assign` / `redundant` が 1 件以上あるときだけ出る:

```
# curation = <state>\tcuration_flags = N[\tcuration_wrong_excluded = N\tcuration_suspect = N[\tcuration_assigned = N\tcuration_redundant_excluded = N]]\tcuration_flags_sha256 = <hex>
```

| フィールド | 意味 |
|---|---|
| `<state>` | 下表 |
| `curation_flags = N` | 有効な判断の総数（`clear` 済みを除くスポットごとの最新 1 行の数。`flags_for_arf2()["n"]`） |
| `curation_wrong_excluded = N` | `state == "applied"` のときだけ。出力の行から除外した `wrong` のスポット数 |
| `curation_suspect = N` | 同上。除外せず残した行のうち `suspect` のものの数（値は変えない） |
| `curation_assigned = N` | `applied` かつ有効な `assign` / `redundant` があるときだけ。`assign` で同定を置き換えて出力に残った行数 |
| `curation_redundant_excluded = N` | 同上。`redundant` で出力の行から除外したスポット数 |
| `curation_flags_sha256 = <hex>` | 判断集合のダイジェスト（`flags_for_arf2()["digest"]`）。再エクスポートの間で判断の内容が変わっていないかを照合できる。`assign` は記録名・`level`・`inchikey`（同じ和組成で別レコードへの付け替えも検出する）、`redundant` は相手と関係コードまで含む |

`state` は経路と引数で決まる:

| `state` | 経路 | 条件 |
|---|---|---|
| `applied` | ARF 経路（解析フォルダのデータセット）は `apply_curation=True`（既定）なら常にこれ（`.arf2` の `MasterAlignmentID` をそのままキーに使うので対応の検証が要らない）。mzTab 経路は `apply_curation=True` かつ、隣接する `.arf` との対応が検証済み（`feature_qc["source"] == "arf"`。SMF_ID が `MasterAlignmentID` と同じ空間） | `wrong` と `redundant` のスポットを出力から除外し、`assign` のスポットの同定を置き換え、件数のフィールドを出す |
| `not_applied` | 両経路 | `apply_curation=False`。判断を無視し、行を落とさず同定も置き換えない |
| `unmapped` | mzTab 経路だけ | `apply_curation=True` だが SMF_ID と `MasterAlignmentID` の対応が未検証（`feature_qc["source"] != "arf"`）。判断があっても行を落とさず `assign` も当てない（mzTab の別の特徴を `.arf2` のスポットと取り違えないため） |

`not_applied` と `unmapped` が件数のフィールドを省くのは、「0 件除外した」と「除外していない」を数字の 0 で混同させないため。

### エクスポートでの `assign` / `redundant` の扱い

両経路で、`state == "applied"` のとき:

- `wrong`: 除外。`redundant`: 除外（`n_unannotated` に数える。「調べていない」行）。
- `assign`: 行の `name`・`ontology`・`inchikey` を記録した値で置き換え、`name_source` と `inchikey_source` を `curation` にする。ARF 経路の `msi_level` は記録名から求め直す（mzTab 経路は空欄）。未注釈だったスポットも InChIKey を得て新しく出力に現れる。
- 記録した InChIKey が空の `assign`（参照レコードに InChIKey が無かった候補。`RIKEN … from …` 名の未知スペクトルなど）は、その行を出力から落とす（MS-DIAL の InChIKey は残さない。`n_unannotated` に数える）。
- 15 列と `CONTRACT_VERSION` は変わらない。

### エクスポートの戻り値の `curation`

成功 payload は両経路で同じ形の `curation` `{"state": ..., "wrong_excluded": N, "suspect": N, "orphaned": N}` を持ち、有効な `assign` / `redundant` があるときは `assigned` と `redundant_excluded` が足される。

| キー | 意味 |
|---|---|
| `state` | 上表の `applied` / `not_applied` / `unmapped`。有効な判断が 0 件（メタ行を出さない）なら `null` |
| `wrong_excluded` / `suspect` | `state == "applied"` のときだけ数字（メタ行の同名フィールドと同じ値）。それ以外は `null` |
| `assigned` / `redundant_excluded` | 有効な `assign` / `redundant` があるときだけ出るキー。`applied` のときだけ数字（メタ行の `curation_assigned` / `curation_redundant_excluded` と同じ値）、それ以外は `null` |
| `orphaned` | 同じファイル名で sha256 が違う（以前の版の）`.arf2` に付いた有効フラグの件数。当てないので行は落とさない。1 件以上なら payload の `warnings` にも出る。メタ行には出さない |

パイプライン（`pipeline_run`）は差次的エクスポートを書かない。フラグを反映したファイルは、`load_dataset(pipeline_path=…)` で開いてから `dataset_compare` → `dataset_export_differential` で書く。

### フラグ記録が壊れているとき

`flags.jsonl` に読めない行（途中で切れた JSON、`spot_id` が整数でない、`flag` が `wrong` / `suspect` / `confirmed` / `clear` / `assign` / `redundant` のどれでもない）があると、読み飛ばさずに止める（読み飛ばすと、その行の `wrong` がエクスポートから外れないまま気づかれない）。

| ツール | 振る舞い |
|---|---|
| `curation_review` / `curation_submit` / `curation_flags` | `{"status": "error", "message", "flags_file", "line"}`。`curation_submit` は壊れた記録に追記しない |
| `dataset_export_differential` | `{"error": {"code": "CURATION_FLAGS_INVALID", "message", "details": {"flags_file", "line"}}}`。ファイルは書かない |
| `dataset_feature_search` | 一覧は返す（判断はもともと反映しない）。読めなかった記録とその理由を名指しし、「キュレーションを反映していない」と書いた注記を `notes` に 1 行足す |

`line` は 1 始まりの行番号。その行を直すか消してから、もう一度実行する。

## 必須注意事項

1. **`UNKNOWN` は「不一致」ではない。** MS/MS 未取得・参照が引けない・参照に RT が無いなどは「確かめられなかった」で、「食い違った」ではない。`UNKNOWN` の系統は `verdict` の `FAIL` / `BORDERLINE` の集計に数えない（`core` トピック §9 のキュレーション版）。
2. **RT–m/z 傾向（`trend_outlier`）は単独では `verdict` を動かさない。** msms/mz/rt/eic の `BORDERLINE` と重なったときの補強に留まる。傾向フィットの `r2` が低いクラスでは `trend_outlier_unreliable`（情報）になる。
3. **Δ（`ppm` / `drt`）は代表試料（`RepresentativeFileID`）の行が基準。** ほかの試料の m/z・RT は使わない。EIC は複数試料を見るが、Δppm / ΔRT は 1 行の値。
4. **`wdot` / `mpp` は MS-DIAL 自身の値で、`curation_review` が照合し直した値ではない。** 再照合の結果は `rescore_discrepancy`（ずれたときだけ）にしか出ない。`wdot` の `-1`（比較不能。`library` トピック §14.2 の番兵）と欠測は TSV ではどちらも空セルなので、「値が無い」と「一致度ゼロ」を分けるには `library_match_feature` を使う。
5. **候補の並び順はスペクトル類似度と制約だけで決まり、MS-DIAL の脂質規則（診断イオン・鎖決定）を評価していない。** 上位が化学的にありえない鎖組成でもありうるので、記録の既定は和組成（`level = sum`）で、`assign` は鎖組成まで確かめたことを意味しない。`likely_wrong` の判定に使う脂質規則フラグは MS-DIAL の照合結果のもので、再検索で出た新しい候補にはその裏付けが無い。
6. **④ の EIC の重ね描きは形の比較用。** 相手 Y のトレースは X の代表試料の最大値に合わせて縮めてあるので、重なって見えても強度が同程度ではない。強度の比は `isotope_ratio`（と `expected_ratio`）で見る。
7. **`redundant` は「X は Y のイオンだ」という判断で、X が実在しないという判断ではない。** `informational`（同位体の強度比が期待の 1.5 倍超）は実在の別物質でありうるので `redundant` に選べない。
8. **`total_score` の `-1` は「比較していない」で、「一致度ゼロ」ではない。** ビューアは `–` と表示する。
9. **`curation_suggest` の `likely_wrong` は元レビューの対象範囲についてしか分からない。** 元レビューを `ontology` などで絞っていれば、範囲外の `likely_wrong` は拾えない（範囲外でも `wrong` フラグと未注釈は対象になる）。
