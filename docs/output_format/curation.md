# アラインメントのキュレーション（`curation_review` 等）の出力フィールド

`curation_review` / `curation_suggest` / `curation_submit` / `curation_flags` が返す値の意味。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。参照ライブラリ照合のスコアの意味（`-1`/`0` の区別など）は `library` トピックが定義し、ここでは繰り返さない。

ツールが**どのファイルのどの関数をどの順に呼ぶか**は `docs/workflow/curation.md`。
ここは**値の意味**だけを定義する。実装は `metabolomix/curation/`
（`judge.py` 機械判別・`evidence.py` 証拠収集・`eic_shape.py` EIC 形状・`trend.py` RT–m/z 傾向・
`adduct_isomer.py` 別アダクトの取り違え・`review.py` レビュー生成/要約・`flags.py` 判断の記録の永続化・`apply.py` エクスポート反映・
`suggest.py` 候補付けの組み立て/保存/送信内容の展開・`candidates.py` 注釈候補 ①②・`relations.py` イオン関係 ④）。

## アラインメントのキュレーション

### `curation_review` の戻り値

| キー | 意味 |
|---|---|
| `review_id` | このレビューの識別子。`curation_submit` / `curation_view_data` に渡す |
| `n_spots` | 対象にしたスポット数 |
| `counts` | 判定の件数内訳 `{"ok": N, "suspect": N, "likely_wrong": N}` |
| `warnings` | レビュー全体にかかる注意（下記） |
| `trend` | クラス（`ontology`）ごとの RT–m/z 傾向フィットの要約 `{n, r2, n_outliers}`。`r2` は外れ値を除いた当てはまりで、低いクラスでは `trend_outlier` 系の理由コードを割り引いて読む。`n_outliers` はそのクラスで頑健 z が `trend_outlier_z` を超えた点の数（クラスの信頼度とは無関係に数える） |
| `table` | `suspect` 以上とフラグ済みのスポットだけの TSV（列は次表）。`ok` かつ未フラグのスポットは載らない。判定の重い順（`likely_wrong` → `suspect` → フラグ済みの `ok`、同順位は `spot_id` 昇順）の**先頭 `max_rows` 行だけ**（引数、既定 100） |
| `n_table_rows_total` | 載せる条件に合うスポットの総数（`max_rows` で切る前） |
| `n_table_rows_shown` | `table` に実際に載せた行数。`n_table_rows_total` より小さければ残りは `html_path` のビューアでしか見られない |
| `table_note` | `table` が全件か先頭だけかの 1 行の注記 |
| `html_path` | ビューア HTML のパス（表示は英語）。ユーザーがブラウザで開き、フラグを付けて「Copy submission text」する。**全件**（`max_rows` で切らない）が入る。自動判別が `likely_wrong` のカードは赤の破線枠（未確認）で、ユーザーが Wrong を付けると実線になる。Class の選択肢には、各クラスの下に「<クラス> › likely_wrong」、先頭に全クラス横断の「likely_wrong (all classes, N)」が並ぶ（該当が 1 件以上あるときだけ）。Verdict（自動）は `ok` / `suspect` / `likely_wrong` / suspect or worse、Flag（人）は Correct / Suspect / Wrong / Flagged で個別に絞り込める。クラス別傾向のカードは Plotly（cdnjs から読む。届かなければ canvas 描画に落ちる）で、**横軸 RT・縦軸 m/z**（範囲はクラスの全点で固定）。点に重ねると `#spot_id 名前`・DB・RT・m/z が出る。カードを押すとそのクラスに絞り込み（もう一度押すと解除）、点を押すとそのクラスに絞ってスポットのカードへ移る（点として拾うのはカーソルが点の上＝hover が出ているときだけ。それ以外の場所はクラスの切り替え）。絞り込みとは別に Exclude のチェック（low score = 理由 `low_score` / no MS/MS = 情報 `msms_absent`）があり、当たるスポットを一覧から隠す。**傾向カードの点は一覧と同じ絞り込み（Verdict・Flag・Exclude、Class の「› likely_wrong」下位項目）に従う**。Class の選択そのものは他のカードの点を消さない（カードの強調だけ）。回帰直線はサーバで当てたまま。隠しても初期選択の Wrong は残り送信される |
| `thresholds` | 実際に使ったしきい値（既定 `DEFAULT_THRESHOLDS` に `thresholds` 引数を上書きしたもの） |
| `ms2_tol` | 対向照合に使った MS2 許容幅（`.dbs` の `search_params` があればそこから、無ければ既定値） |

`warnings` に出る 4 通り: (1) 照合結果を持つスポットのうち参照を引けた割合が半分未満
（アラインメントと別のライブラリを読んでいる可能性）、(2) `.dcl`/`.EIC.aef`/`.arf` の
兄弟ファイルが見つからない（その系統の判定は `UNKNOWN` になる）、(3) **orphaned**——
`flags.jsonl` に、同じファイル名で sha256 が違う（＝以前の版の）`.arf2` に対して記録された
有効フラグがある（件数を出す。スポット番号が今の版と対応する保証が無いので当てない）、
(4) 照合結果はあるのに MS-DIAL の脂質規則フラグ（`is_lipid_class_match` /
`is_lipid_chains_match` / `is_other_lipid_match`）がどのスポットにも立っていない——脂質以外の
採点器の出力とみなし、規則に基づく 3 コード（`class_rule_rejected` / `class_rules_not_run` /
`chains_unsupported`）を使わない（次節）。

エラーで止まる場合（`{"status": "error", "message": ...}`）:
`file_ids` に `.arf` の行に無い試料 ID がある（`missing_file_ids` にその ID を出す。重い
読み込みの前に止まる）、`flags.jsonl` に読めない行がある（`flags_file` と `line` を出す。
下記「フラグ記録が壊れているとき」）、対象が上限（3000 スポット）を超える、など。

**共有の注意**: `review-<id>.html` / `review-<id>.json` は、読み込んだ参照ライブラリの
スペクトル（対向プロットの参照側）をそのまま埋め込む。研究室のライブラリを読んで作った
レビューは、ライブラリ本体と同じく**研究室の外へ出さない**。

### ビューアの証拠 payload（`curation_view_data` / HTML ビューア専用）

`curation_review` 自体の戻り値には出ない（座標を LLM に返さない規約。上記 `table` と
`html_path` のみ）が、`html_path` のビューアと `curation_view_data`（MCP Apps 専用、
LLM は呼ばない）はスポットごとに EIC 系列（`eic.samples[].points`）と対向スペクトル
（`mirror.measured` / `mirror.reference`）を持つ。判定（`eic_shape` / `rescore` /
`matched_mz` 等）は必ず全点・全スペクトルで計算した**後**に、返す座標列だけを
間引く（payload only。判定は変えない。定数・実装は `metabolomix/curation/evidence.py`）:

- `eic.samples[].points`: 1 サンプルあたり最大 `EIC_MAX_POINTS`（既定 40）点。先頭・
  末尾・頂点（最大強度）・積分範囲（`left`/`right`）に最も近い点は必ず残し、残りは
  トレース全体に等間隔に散らした位置で埋める（アンカーとの和集合が上限に収まるまで
  等間隔の点数を減らす。先頭から詰めないので、ピークの片側だけが落ちることはない）。
  強度は整数に丸める。
- `mirror.measured` / `mirror.reference`: それぞれ最大 `MIRROR_MAX_PEAKS`（既定 150）点。
  一致した点（測定側は `matched_measured_mz`、参照側は `matched_mz`）は必ず残し、
  残りは強度降順で埋める。`scored_peak_count` / `unscored_peak_count` は**間引き前の
  満スペクトルの値のまま**返す（間引き後の表示点数には合わせない）。

実データ check（kidney pos, 2196 spots）で JSON が 41.7 MB（内訳 EIC 31.0 MB・
mirror 5.5 MB）になったため導入した上限。

ビューア用にスポットごとに次も持つ:

- `auto_note`: 判定根拠の文（`likely_wrong` / `suspect` のとき。`ok` は `null`）。
  `Auto: ` に続けて理由コードごとの根拠（英語）を ` / ` でつなぐ（例
  `Auto: m/z — Δm/z 12.4 mDa (≥10 mDa) / MS2 — no match to the reference (low score)`。
  文面は `judge.REASON_TEXT`。2026-10-09 に日本語から英語へ変えた。それより前に記録された
  `flags.jsonl` のメモは日本語のまま）。ビューアのメモ欄の既定値（記録済みのメモがあればそちら）。
- `flag_cleared`: このアラインメントで最新のフラグ行が `clear`（人が明示的に取り消した）。

- `isotopes`: MS1 の同位体パターン（ビューアの MS1 パネル用。判定には使わない）。
  `measured` は `.arf2` Key 53 `IsotopicPeaks`（代表試料の M, M+1, M+2）の `[m/z, 相対強度 %]`
  （M = 100）。MS-DIAL が書く m/z は単同位体 + 1.00467·k の計算値で、実測の質量ではない。
  `theoretical` は `{"relative": [M, M+1, M+2 の %], "basis": ...}` で、組成式（`Formula`）に
  アダクトの原子を足して整数質量の分解能で畳み込んだもの（`curation/isotope.py`）。
  `basis` は `formula+adduct`、アダクトが読めなければ中性の組成式だけの `formula`。
  組成式が無い・読めなければ `null`。

ビューアは、記録済みのフラグが無く `flag_cleared` でもない次のスポットを、開いた時点で
Wrong の未送信の変更にする: `likely_wrong`、理由に `low_score`（MS-DIAL の名前の
`low score:`）、情報に `msms_absent`（`no MS2:` など。MS/MS なし）。メモは `auto_note`
（MS/MS なしは判定理由ではないので `Auto: no MS/MS` を書き足す）。そのまま送信すれば記録
され、Correct（無印。何も記録しないことの表示名）に戻して送信すれば `clear` が記録
される（次のレビューで掛け直さない）。**送信すると `_tags.xml` の Misannotation にも
反映される**ので、残したいスポットは送信前に Correct へ戻す。
各スポットのカード右上の **Confirmed** チェック（2026-10-09）は判断 `confirmed`（注釈が正しいと確かめた）で、
Suspect / Wrong と排他（チェックすると Correct に、Suspect / Wrong を選ぶと外れる。外すと Correct＝送信では `clear`）。
初期状態は記録の `confirmed` か、レビュー作成時の `_tags.xml` の Confirmed（MS-DIAL の GUI で付けたものも）で、
スポットの `confirmed`（真偽値）に出る（最新の記録が `wrong` / `suspect` なら `false`）。`confirmed` のスポットは
Wrong を初期選択しない。メモが自動の判定根拠のままならチェックしたときに空にする。Flag の絞り込みは
Correct（無印か confirmed）/ Confirmed / Suspect / Wrong / Flagged（Suspect か Wrong）。
注釈を確かめたスポット（`curation/apply.py` `confirmed_spots`: 記録の `confirmed` ∪ `_tags.xml` の Confirmed −
最新の記録が `wrong` / `suspect`）には自動判定の likely_wrong を当てない: 初期選択・`exclude_auto_likely_wrong`
（`arf_plot_species` / `arf_pca_species` / `arf_plot_group_intensity`）・`curation_suggest` の likely_wrong 対象。
ミラープロットは縦軸に相対強度（上下とも各側の最大値を 100）の目盛り 0/50/100 を付け、
ピークに m/z（小数 4 桁）のラベルを `library_plot_mirror` と同じ規則（強度降順の貪欲法・
縦横とも重なるものを飛ばす・片側 25 本まで・側ごとに独立）で付ける。

### `table`（TSV）の列

| 列 | 意味 |
|---|---|
| `spot_id` | `MasterAlignmentID` |
| `name` | ARF2 の `Name`（確度接頭辞を含む生の文字列。§9-1 の接頭辞規則は `core` トピック） |
| `ontology` | ARF2 の `Ontology`（脂質クラス） |
| `adduct` | ARF2 の `AdductType` |
| `verdict` | `likely_wrong` / `suspect` / `ok`（次節） |
| `reasons` | 立った理由コードをカンマ区切りで、強い→弱い→帯のみの順に並べたもの（次節の表と同じ順）。相手のスポットを指すコードは `adduct_isomer_of:189` のように `:<相手の spot_id>` が付く |
| `ppm` | 代表試料の m/z と参照 precursor m/z（引けなければ Formula からの理論値）の相対誤差 [ppm] |
| `dmz_mda` | 代表試料の m/z と参照 precursor m/z（引けなければ Formula からの理論値）の差 [mDa]。`ppm` と同じ基準 |
| `drt` | 代表試料の RT と参照 RT の差 [分]。参照が引けない、または参照に RT が無ければ空。MS-DIAL が RT を使っていない照合（情報 `rt_not_used_by_annotation`）でも値は出るが、判定には使わない |
| `wdot` | **MS-DIAL 自身が出した** `squared_weighted_dot_product` の平方根（`library` トピック §14.1 と同じ規約: 平方根側の値）。この照合結果は同定時点でアラインメントに保存済みのもので、`curation_review` がここで再照合した値ではない（後述の必須注意）。値が無い、または MS-DIAL 側が `-1`（比較不能）なら空 |
| `mpp` | **MS-DIAL 自身が出した** `matched_peaks_percentage`。`wdot` と同じく再照合値ではない |
| `eic_good` | 検出（非 gap-fill）サンプルのうち EIC 形状が「良い」と判定された割合（`eic_shape.spot_shape` の `good_fraction`）。良否は頂点がウィンドウ内・点数十分・理想ガウスとの R² 十分・極大数が上限以下で決まる。検出サンプルが無ければ空 |
| `trend_z` | 所属クラスの RT–m/z 傾向モデル（Huber 頑健回帰、RT = a + b·炭素数(+ c·二重結合数)）に対する頑健 z 値。傾向がフィットできないクラス（点数不足・炭素数のばらつき不足）では空 |
| `flag` | 既存のキュレーションフラグ（`wrong` / `suspect`）。未フラグは空。**`clear` はここに出ない**——`clear` は取り消しの記録で、有効フラグ（スポットごとの最新 1 行）が `clear` のスポットは未フラグとして扱うため |

### 判定（`verdict`）と理由コード

`verdict` は 6 系統の帯（`checks.{msms,mz,rt,eic,ion,trend}`。各 `PASS`/`BORDERLINE`/`FAIL`/`UNKNOWN`
と理由コード）から決める。`ion` は別アダクトの取り違え（下記）で、`FAIL` か `PASS` のどちらか
（`adduct_isomer` が無い＝評価していないときだけ `UNKNOWN`）:

- **`likely_wrong`**: 強い理由コードが 1 つでも立った場合。
- **`suspect`**: 強い理由が無く、かつ (a) msms/mz/rt/eic/ion のいずれかが `FAIL`、または
  (b) msms/mz/rt/eic/ion のうち `BORDERLINE` が 2 つ以上、または (c) `BORDERLINE` が 1 つ以上
  あって RT–m/z 傾向も `BORDERLINE`（`trend_outlier`）の場合。
- **`ok`**: それ以外。

理由コードは 3 つの区分に分かれる。**強い**（`STRONG_REASONS`）はそれ単独で `likely_wrong`
にする。**弱い**（`WEAK_REASONS`）は当該系統を `FAIL` にするが単独では `suspect` 止まり。
**帯のみ**は `BORDERLINE` の内訳で、`suspect` の積み上げにのみ数える。**情報**は
`verdict` を動かさない注記（`info`。`checks` とは別枠で常に出る）。

| 区分 | コード | 意味 |
|---|---|---|
| 強い | `polarity_mismatch` | アダクトの電荷符号と実測イオン化極性（`IonMode`）が不一致（`adduct_consistency` の `band` が `FAIL`）。脂質クラスとの典型性（`class_typical`）は advisory のみで `band` には効かない——非典型アダクトだけでは立たない |
| 強い | `precursor_unmatched` | MS-DIAL 自身の `is_precursor_mz_match` が `False` |
| 強い | `dmz_out` | \|Δm/z\|（`dmz_mda`）が `dmz_fail_mda`（既定 10 mDa）**以上**。ppm は m/z に比例して緩むので絶対差で見る（ユーザー決定 2026-09-29）。Δppm の `ppm_out` は弱いまま |
| 強い | `class_rule_rejected` | MS/MS ありで、MS-DIAL の脂質クラス規則（診断イオン）を評価して棄却した（`is_lipid_class_match=False` かつ `is_other_lipid_match=False`）。実測では全件 `low score:` なので `low_score` も同時に立つ。**脂質規則が走ったデータに限る**（下記） |
| 強い | `adduct_isomer_of:<spot>` | m/z が、同時に溶出する別物質のスポット `<spot>` の中性質量の**別のアダクト**で説明でき、かつ `<spot>` の証拠の段階が上（下記「別アダクトの取り違え」）。例: neg #173「PI 41:2」[M-H]- は #189 DGDG 35:1 [M+CH3COO]- の [M-H]-（−1.2 ppm）。単独の証拠では `ok` だったスポットにも立つ |
| 弱い | `ppm_out` | Δppm が `ppm_borderline` しきい値（既定 10）を超えた。**adduct 非依存**（実測: kidney neg/pos で全 adduct の中央値が約 −0.8 ppm）で、単独では強い理由に数えない（ユーザー決定 2026-09-29）——mz 系統自体は `FAIL` になるが、単独では `suspect` 止まり。`polarity_mismatch`/`precursor_unmatched` が別途立てば、そちらの強さで `likely_wrong` になる |
| 弱い | `low_score` | MS/MS はあるが MS-DIAL 自身の `is_reference_matched` が `False` |
| 弱い | `drt_out` | ΔRT が `drt_borderline` しきい値（既定 1.0 分）を超えた |
| 弱い | `eic_poor` | EIC 形状帯が `FAIL`（検出サンプルが 1 件以上あり、かつ `good_fraction` が `eic_borderline_frac` 未満）。**検出サンプルが 0 件のときは `good_fraction` が計算できず帯は `UNKNOWN` になり、このコードは立たない**（`spot_shape` は `n_detected==0` なら `fraction=None`→`band="UNKNOWN"`） |
| 弱い | `adduct_isomer_minor_of:<spot>` | `adduct_isomer_of` と同じ条件で説明できるが、`<spot>` の証拠の段階は同じで、強度（`HeightAverage`）が `adduct_isomer_height_ratio`（既定 3）倍以上。例: neg の RIKEN ID 注釈は FA 18:0 の [2M-H]-（強度 41 倍） |
| 帯のみ | `ppm_borderline` | Δppm が `ppm_pass`〜`ppm_borderline`（既定 5〜10）の帯 |
| 帯のみ | `drt_borderline` | ΔRT が `drt_pass`〜`drt_borderline`（既定 0.5〜1.0 分）の帯 |
| 帯のみ | `eic_borderline` | EIC 形状帯が `BORDERLINE`（検出サンプルが 1 件以上あり、かつ `good_fraction` が `eic_borderline_frac`〜`eic_pass_frac` の帯）。これも検出 0 件では立たず `UNKNOWN` になる |
| 帯のみ | `rt_scatter` | 検出サンプル間で EIC 頂点 RT のばらつき（母標準偏差 `statistics.pstdev`。検出 2 件以上のときだけ計算）が `eic_rt_scatter_sd`（既定 0.1 分）を超えた。EIC 帯が `PASS` ならこの 1 件だけで `BORDERLINE` に格下げする |
| 帯のみ | `trend_outlier` | 所属クラスの RT–m/z 傾向で頑健 z が `trend_outlier_z`（既定 3.0）を超え、かつそのクラスの傾向フィット自体が信頼できる（`r2 >= trend_min_r2`、既定 0.7） |
| 情報 | `msms_absent` | MS/MS 未取得、または名前接頭辞が `no MS2`/`w/o MS2`（`msms` 系統は `UNKNOWN`） |
| 情報 | `reference_not_found` | ライブラリから参照レコードを引けなかった。**`rt` 系統だけが `UNKNOWN`** になる。`mz` 系統は Formula/AdductType からの理論値（`mass_error_ppm`、`ppm_basis="formula"`）にフォールバックして計算を続け、それも失敗したときだけ `UNKNOWN` になる（`rt` と違って自動的に `UNKNOWN` にはならない） |
| 情報 | `reference_rt_absent` | 参照は引けたが RT を持たない（`rt` 系統は `UNKNOWN`） |
| 情報 | `rt_not_used_by_annotation` | MS-DIAL がこの照合で RT を絞り込みにも採点にも使っていない（`.dbs` に保存された注釈器の `use_time_for_annotation_filtering` / `use_time_for_annotation_scoring` がどちらも False）。参照 RT は別のクロマトグラフィーを前提にした予測値でありうる（例: LBM は約 18 分系の予測 RT を持つ）ので、`rt` 系統は `UNKNOWN` にして `drt_out` / `drt_borderline` を立てない。`drt` の値は表示用に残る。`.msp`/`.lbm2` を直接読んだ store のようにスイッチが分からないときは、従来どおり RT で判定する |
| 情報 | `no_match_result` | ARF2 に MS-DIAL の照合結果（`representative`）自体が無い（`msms` 系統は `UNKNOWN`） |
| 情報 | `class_rules_not_run` | 脂質クラス規則が評価されていない（`is_lipid_class_match=False` かつ `is_other_lipid_match=True`。CompoundClass が Unknown/Others・SPLASH・名前解析失敗など）。参照一致でも規則の裏付けは無い＝**誤りの意味ではない**（未検証）。脂質規則が走ったデータに限る |
| 情報 | `chains_unsupported` | 名前が鎖レベル（`16:0_18:1` や `18:1;O2/16:0` のように鎖を `_`/`/` で区切る。`PC 34:1\|PC 16:0_18:1` の `\|` 以降も含む。単鎖は数えない）なのに `is_lipid_chains_match=False`。MS-DIAL は照合に失敗しても参照名から `\|` 付きの名前を作るので、鎖組成の裏付けは `is_lipid_chains_match=True` だけ。MS/MS ありのときだけ。脂質規則が走ったデータに限る |
| 情報 | `rescore_discrepancy` | `curation_review` が対向照合で出した `weighted_dot_product` が MS-DIAL 自身の値（平方根換算）と `rescore_tolerance`（既定 0.1）を超えてずれた |
| 情報 | `adduct_differs_from_reference` | 参照レコードのアダクトと注釈のアダクトが食い違う |
| 情報 | `manually_modified` | MS-DIAL 側で `is_manually_modified` が立っている（GUI で人手修正済み） |
| 情報 | `manually_unsettled` | 名前接頭辞が `unsettled`（MS-DIAL Alignment Viewer で未確定のまま） |
| 情報 | `trend_outlier_unreliable` | 頑健 z はしきい値を超えたが、そのクラスの傾向フィット自体が信頼できない（`trend_outlier` には数えない） |
| 情報 | `dcl_precursor_mismatch` | `.dcl` の該当インデックスの precursor m/z が代表試料の m/z と 0.01 以上ずれる（順序対応の前提が崩れている可能性。その場合 MS/MS 系統の実測スペクトルは使わない） |

**MS-DIAL で絞り込まれたデータでは `likely_wrong` はまれ。** 強い理由の
`precursor_unmatched` は MS-DIAL 自身の `IsPrecursorMzMatch` で、MS-DIAL は同定時に
precursor の許容幅をすでに課している——許容幅の外の候補はそもそも代表に残りにくい。
`likely_wrong` が 0 件でも「全部正しい」ではなく、`suspect` の中身（理由コード）を読むこと。
`class_rule_rejected`（MS-DIAL の脂質クラス規則による棄却）は、そうした絞り込みの後でも
残る強い理由で、kidney neg/pos の実測では 370 / 404 件あった。

**別アダクトの取り違え（`adduct_isomer_of` / `adduct_isomer_minor_of`）。** 個々のスポットの
証拠（Δppm・MS2）だけでは見えない誤り——別の物質の別アダクトを、たまたま近い質量の別物質として
注釈したもの（PI 41:2 と DGDG 35:1 は [M-H]- で 8 ppm しか違わない）——を、同時に溶出する
スポット同士で見る（`metabolomix/curation/adduct_isomer.py`）。対象 X と相手 Y について:

- **相手の範囲**: レビューの対象（`ontology` / `name_contains`）に絞らず、アラインメントの
  注釈付きスポット全部。有効な `wrong` フラグのあるスポットは相手にしない。
- **条件**: Y が X と同じ極性で |ΔRT| ≤ `adduct_isomer_drt`（既定 0.05 分）。X と Y は別物質
  （組成式とクラス（`Ontology`）が両方同じなら同じ分子種の別アダクト＝正しい注釈なので対象外）。
  Y の m/z と注釈アダクトから中性質量を出し、仮定アダクト A（X の極性の既定アダクト一覧
  `msdial.analysis_params.DEFAULT_ADDUCTS` と X 自身のアダクト。Y 自身のアダクトは除く）での
  m/z が X の m/z と `adduct_isomer_ppm`（既定 5 ppm）以内。m/z・RT はどちらもスポットの
  `MassCenter` / `RT`。
- **証拠の段階**: 3 = MS2 で鎖組成まで裏付け（鎖レベルの名前 ∧ `is_lipid_chains_match`）、
  2 = MS2 が参照と一致（`is_reference_matched`）、1 = MS2 はある（low score）、0 = MS2 無し
  （名前接頭辞 `no MS2`/`w/o MS2` を含む）・照合結果無し。Y の段階が上なら `adduct_isomer_of`
  （強い）、同じ段階で強度比が `adduct_isomer_height_ratio` 以上なら `adduct_isomer_minor_of`
  （弱い）。Y の段階が下なら立たない（逆向きの組は Y の側に立つ）。
- 複数の Y が説明するときは、強い → 段階 → 強度の順に 1 件だけ。

レビュー JSON（とビューア）のスポットは `adduct_isomer` に根拠を持つ: `of`・`of_name`・`of_adduct`
（相手と、その注釈アダクト）・`as_adduct`（X を説明する仮定アダクト A）・`ppm`（X の m/z と A の
理論 m/z の差）・`drt`（X − Y の RT [分]）・`severity`（`strong` / `minor`）・`tier` / `of_tier`・
`height_ratio`（Y / X）。説明する相手が無ければ `null`。**同じ組成式・同じクラスの組（例: DG の
[M+Na]+ と [M+NH4]+ の二重計上）はこの理由の対象外**（重複の整理は別の話）。20260930_EV の実測で
neg 154 件中 9 件・pos 469 件中 27 件に立った。

**脂質規則フラグは脂質規則が走ったデータに限って読む。** 規則フラグは MS-DIAL の Lipidomics
採点器でしか立たず、それ以外の採点器では全部 `False` になる（`class_rule_rejected` の条件と
区別できない）。そこで 1 回のレビューの対象スポットに規則フラグ（class / chains / other）が
1 件でも `True` のものがあるときだけ、上の 3 コードを使う（`judge.lipid_rules_active`）。
1 件も無ければ `warnings` にその旨を出す。フラグの意味は上流 MsdialWorkbench `afd5f9522` の
`MsReferenceScorer` / `LipidMsmsCharacterization` と kidney の実測で確かめた。

### `curation_suggest` の戻り値

間違いになったスポットと未注釈スポットに、注釈の候補（①②）と別スポットのイオンとしての説明（④）を並べる。
**先に `library_load` と `curation_review`**（どちらも無ければ `missing_state`）。

| キー | 意味 |
|---|---|
| `suggestion_id` | この候補付けの識別子（`cs-YYYYMMDD-HHMMSS-xxxx`）。`curation_submit` の送信用テキストが運ぶ |
| `base_review_id` | `likely_wrong` の判定を読んだ元のレビュー。`review_id` を省くと、このアラインメント（sha256 一致）の最新のレビュー |
| `n_spots` | 対象になったスポット数 |
| `counts` | `targets`（対象の種類別 `{flagged, likely_wrong, unannotated}`）・`with_candidates`（候補または情報でないイオン関係があるスポット数）・`with_strong_relation`（強い説明（後述）のあるスポット数）・`hard_removed`（ハード制約で削った候補の延べ数） |
| `warnings` | param ファイルが無く既定の検索アダクトと RT 窓を使った、`PeakProperties.arf` が無く相関を計算できない、兄弟ファイルが無い、など。**`likely_wrong` は元レビューの対象範囲についてしか分からない**旨は常に出る |
| `table` | 1 スポット 1 行の TSV（列は下表）。強い説明のあるスポットが先、次に対象の種類（`flagged` → `likely_wrong` → `unannotated`）、同順位は `spot_id` 昇順の**先頭 `max_rows` 行だけ**（引数、既定 100）。載せた行数は `n_table_rows_shown`。全件は `html_path` |
| `html_path` | 候補付けビューア HTML のパス。ユーザーがブラウザで開いて選び、「送信用テキストをコピー」でチャットへ貼る |
| `library` | 使った参照ライブラリの `path` と `sha256` |
| `analysis_params` | 検索アダクトと RT 窓の出所。`source`（param ファイルか既定か）・`path`・`rt_window`（分） |
| `options` / `thresholds` | 実際に使った引数（`rt_window` は解決後の値）としきい値 |

`table`（TSV）の列:

| 列 | 意味 |
|---|---|
| `spot_id` | `MasterAlignmentID` |
| `name` | 今の名前（未注釈は `Unknown` 系） |
| `target` | 対象の種類。`flagged`＝有効な `wrong` フラグのあるスポット、`likely_wrong`＝元レビューの判定が `likely_wrong`（フラグ無し）、`unannotated`＝`.arf2` に代表が無いスポット。`include_decided=True` のときは、判断済み（assign / redundant）の注釈付きスポットも `flagged` で戻る |
| `top_candidate` | 候補（①②）の先頭の名前（和組成があれば和組成）。候補が無ければ空 |
| `source` | 先頭の候補の出所。`msdial`（MS-DIAL 自身の下位候補）・`research`（閾値を緩めた再検索だけで出たもの）・`msdial+research`（両方） |
| `total_score` | 先頭の候補の総合スコア（`library` トピックの `total_score` と同じ物差し。RT・precursor の一致度を足した**非正規化の和**で 1 を超える）。比較できなければ `-1` |
| `reasons` | 先頭の候補のソフト理由・情報の理由コードをカンマ区切り |
| `relation` | 先頭のイオン関係（④）を `<関係コード>(#<相手スポット>)` で。無ければ空 |
| `strong` | 強い説明があるか（`True`/`False`） |

### 候補付けの payload（HTML ビューア専用）

`curation_suggest` の戻り値には出ない（座標を LLM に返さない規約）。`suggest-<id>.json`（正準）と
`suggest-<id>.html` のビューアだけが持つ。スポットごとに `eic`（対象の EIC）・`measured`（測定 MS/MS）・
`candidates`（①② の候補）・`relations`（④ の関係）・`strong`・`preset`（強い説明があるときの初期選択の `R<n>`）を持つ。
`curation_submit` はここから候補の詳細を引いて記録行に展開する。

**候補（`candidates[]`、①②）**: `candidate_id`（`L1`…、並び順）・`source`・`name`（分子種）・`sum_name`
（和組成。解析できなければ `null`）・`ontology`・`adduct`・`formula`・`inchikey`・`library_id`/`record_index`
（参照レコード）・`dmz_mda`/`ppm`（代表試料の m/z と参照 precursor m/z の差）・`scores`
（`total_score` `weighted_dot_product` `simple_dot_product` `reverse_dot_product`
`matched_peaks_percentage` `matched_peaks_count`）・`trend`（RT–m/z 傾向の予測 RT・残差・z。当てはまらなければ `null`）・
理由コード（`soft` / `info`。ハードは削られるので残らない）・`mirror`（対向プロットの参照側）。
`total_score` 等の `-1` は**比較していない**（スポットに MS/MS が無い、または参照が引けずスペクトルが無い。`matched_peaks_count` は 0 のまま）で、
ビューアは `–` と表示する（「一致度ゼロ」ではない。`library` トピックの `-1`/`0` の区別と同じ）。

**イオン関係（`relations[]`、④）**: 次々節。

ビューアは、イオン関係の行を選ぶ（または指す）と、相手 Y の代表試料の EIC（`partner_eic`）を X の代表試料の EIC に
重ねて描く。Y のトレースは X の代表試料の最大値に合わせて縮め、点線で描く（形の比較用。強度の比は `isotope_ratio` を見る）。

### 候補の理由コード

候補（①②）ごとに、次の制約を当てる。**ハード**は候補から削る（戻り値の `hard_removed` に数だけ残る）、
**ソフト**は順位を下げる、**情報**は何もしない。並び順は（ソフト理由の数 昇順、`total_score` 降順）で、
スポットごとに上位 `top_n`（既定 5）件を出す。

| 区分 | コード | 意味 |
|---|---|---|
| ハード | `polarity_mismatch` | 候補のアダクトの電荷符号とスポットのイオン化極性が矛盾 |
| ハード | `dmz_out` | \|Δm/z\|（代表試料の m/z と参照 precursor m/z の差）が `dmz_fail_mda`（既定 10 mDa）**以上** |
| ソフト | `adduct_atypical` | そのクラスに典型的でないアダクト |
| ソフト | `trend_outlier` | 候補のクラス×和組成から予測した RT との残差の \|z\| が `trend_outlier_z`（既定 3.0）を超えた |
| ソフト | `no_matched_peaks` | MS/MS はあるが、比べた結果の一致ピークが 0（参照が引けず比べていない候補には付けない。それは `reference_unresolved`） |
| ソフト | `msms_absent` | スポットに MS/MS が無い（①の候補にだけ付く。②は MS/MS が無いと働かない） |
| ソフト | `no_support`（イオン関係） | ④の関係に裏付けが無い（次節） |
| 情報 | `trend_unknown` | 候補のクラスの傾向が点数不足などで当てはまらない。**順位は下げない**（「傾向から外れた」ではなく「傾向を持たない」） |
| 情報 | `reference_unresolved` | ①の候補の参照レコードがライブラリに無い（スペクトルの再採点ができず、`scores` は `-1`） |

傾向モデルは、そのアラインメントの注釈付きスポットのうち、対象・`wrong` フラグ付き・判断済みの `redundant` のスポットを除いたものから、
（`assign` を決めたスポットは記録した名前・クラスで）クラスごとに
`RT = a + b·炭素数 + c·二重結合数`（Huber）で当てはめる。元レビューの範囲には縛られない。

### イオン関係（`relation`）

対象スポット X を、注釈付きの別スポット Y のイオンとして説明できるとき、`relations[]` の行（`candidate_id` は
`R1`…）になる。①② とは物差しが違うので、候補と 1 本の順位には混ぜない。

**相手 Y の条件**: 注釈付きで、有効な `wrong` フラグが無く、今回の対象でも `redundant` 判断済みでもないスポット。
誤った注釈や消した行を起点にした説明の連鎖を防ぐ。**Y が `assign` 済みなら、記録した名前・組成式・アダクトで扱う**
（置き換えた MS-DIAL の名前や組成式で同位体の期待比を計算しない）。

**関係になる条件**（すべて）: (1) \|ΔRT\| ≤ `rt_window`（既定は param ファイルの
`Retention time tolerance for alignment`、無ければ 0.1 分）、(2) Δm/z が下の関係のどれかで `relation_mz_tol`
（既定 10 mDa）以内、または MS-DIAL の `found_in_upper_msms` リンクがある。

| 関係コード | 意味 |
|---|---|
| `adduct:<X のアダクト>/<Y のアダクト>` | 同じ中性分子の別アダクト。Y の中性質量に、解析で検索したアダクト（param ファイルの `Searched adduct ions`。二量体・多価を含む。無ければ極性ごとの既定）を付けた m/z が X に一致 |
| `isotope_M+1` / `isotope_M+2` | Y の 13C 同位体（+1.003355 / +2.006710 を電荷で割る） |
| `insource:-H2O` `-2H2O` `-NH3` `-HCOOCH3` `-CH3COOCH3` `-C3H5NO2` `-C2H8NO4P` `-C3H8NO6P` `-C6H10O5` | Y からの中性損失によるインソース断片（`metabolomix/curation/relations.py` の表。質量はそれぞれの組成式から計算） |
| `found_in_upper_msms` | 質量差が上のどれでも説明できないが、MS-DIAL が「Y の MS/MS に X が見える」とリンクしている |

**MS-DIAL の `found_in_upper_msms` リンクは向きの無い対で保存される**（実データで相互に 100%、相手が高 m/z 側になるのは半々）。
そこで前駆体（upper）側は m/z の大小で決める: **リンクだけで関係にするのも、リンクで強い説明にするのも、X の m/z が Y の m/z
より小さい（X が断片＝軽い側）ときに限る**。リンクは裏付け（下）には向きによらず数える。

**裏付け**: MS-DIAL のリンク（`correl_similar` / `chrom_similar` / `found_in_upper_msms`）があるか、`.arf` の試料別の高さの
log1p で求めた Pearson の `r` ≥ `relation_min_r`（既定 0.8。両方の高さが正の試料が 5 以上のときだけ計算、足りなければ `r` は空）。
どちらも無い関係は載せるが `no_support`（ソフト）を付けて順位を下げる。`msdial_links` は X–Y 間の MS-DIAL のリンク種別。

**強い説明**（`strong = true`）: 裏付けがあり、情報でなく、かつ (a) `found_in_upper_msms` の関係、(b) Y が前駆体側（上の m/z の規則）の
`found_in_upper_msms` リンクを持つ関係、または (c) 同位体で `isotope_consistent` のもの。強い説明のあるスポットでは、それがカードの先頭で
初期選択（`preset`）になる。**アダクトの組・インソース断片は、Y が前駆体側のリンクを伴わない限り強い説明にならない**（Δm/z が合うだけでは
断定しない）。

**同位体の強度比の規則**: `isotope_ratio` は試料ごとの X の高さ / Y の高さの中央値（両方が正の試料だけ）。`expected_ratio` は Y の組成式の
炭素数 n から二項分布で求めた期待比 `C(n,k)·q^k`（k = 1 または 2、q = 0.0107 / (1 − 0.0107)）。`isotope_ratio` が `expected_ratio` の
**1.5 倍以下**なら `isotope_consistent = true`（Y の同位体として矛盾しない）。**1.5 倍を超えれば X は実在の別物質**（同位体の重なり）と
みなし、`informational = true`（情報として添えるだけで、`redundant` には選べない）。比が求められない（炭素数不明・共通試料なし）ときは
`isotope_consistent` は `null` で、強い説明にはならない。

関係の並び順は（強い説明 → 情報でないもの → 裏付けありのもの、`|Δm/z|` 昇順）。ビューアのために、上位 3 件の関係は相手 Y の代表試料の
EIC（`partner_eic`）を持つ。

### 判断の記録（assign / redundant）

`flags.jsonl`（追記専用。ファイル名は変わらず、「フラグ」は「判断の記録」の意味に広がった）の `flag` は
`wrong` / `suspect` / `clear` / `assign` / `redundant`。キー（アラインメントの sha256 と `MasterAlignmentID`）、
**スポットごとに最新の 1 行が勝つ**、`clear` は取り消し、は既存どおり。`wrong` の後の `assign` は `assign` が勝つ。

| flag | 主な項目 |
|---|---|
| `assign` | `name`（記録名）・`level`（`sum` / `species`）・`species_name`（レコード本来の分子種名）・`ontology`・`adduct`・`formula`・`inchikey`・`candidate_source`・`library`（`sha256` `library_id` `record_index`）・`scores`（採点のスナップショット）・`suggestion_id` |
| `redundant` | `of`（Y の `spot_id`）・`relation`（関係コード）・`evidence`（`dmz_mda` `drt` `r` `msdial_links`）・`suggestion_id` |

- 候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign|redundant|clear, candidate, level, note}]`。`candidate` は `assign` なら
  `L<n>`、`redundant` なら `R<n>`（**情報のイオン関係は `redundant` に選べない**）。名前・InChIKey は送信側で書かず、
  保存済みの `suggest-<id>.json` から展開する。`level` の既定は `sum`（和組成。`name` = `sum_name`、解析できなければ分子種名）で、
  `species` なら `name` はレコードの分子種名。**和組成で記録しても、選んだレコードの InChIKey を持たせる**（MS-DIAL が和組成の名前にも
  ライブラリの構造の InChIKey を付けるのと同じ）。この InChIKey は MS/MS で鎖組成まで確かめたことを意味しない。
- `curation_submit` の戻り値に `n_assign` / `n_redundant` / `n_confirmed`（有効な判断の件数）が加わる。**`assign` / `redundant` は `_tags.xml` の
  Misannotation を変えない**（MS-DIAL の中の名前はまだ誤ったままのため）。
- 候補付けビューアの「元の注釈に戻す」は `clear` を送る。`clear` はそのスポットの**有効な判断を丸ごと**（`wrong` でも `assign` でも）消し、
  `_tags.xml` の Misannotation も外す。
- `curation_flags` の TSV は 7 列: `spot_id`・`flag`・`name`（`assign` の記録名）・`of`（`redundant` の相手）・`note`・`source`・`ts`。
- `arf2_annotate_identities` の `curation_flag` 列は有効な判断: `wrong` / `suspect` / `assign:<記録名>` / `redundant`（無ければ空）。
- `curation_review` のスポットの `flag` / `flag_note` は `wrong` / `suspect` / `confirmed` だけ。`assign` / `redundant` は別項目 `decision`
  （`{flag, name}` または `{flag, of}`）に出る。
- エクスポートでの扱い（下の「エクスポートのメタ行」以降）: `assign` は同定を置き換え、`redundant` は除外する。

### `curation_submit` の戻り値の `tags_xml`（MS-DIAL への反映）

記録（`flags.jsonl`、正本）の後、アラインメントの `<.arf2 の stem>_tags.xml` の
**Misannotation**（タグ Id 3）と **Confirmed**（タグ Id 1）に反映する（2026-10-09 に Confirmed を追加）:

| 判断 | Confirmed | Misannotation |
|---|---|---|
| `confirmed` | 付ける | 外す |
| `wrong` | 外す | 付ける |
| `suspect` | 外す | 触らない |
| `clear` | 外す | 外す |

**`assign` / `redundant` は触らない**（MS-DIAL の中の名前はまだ誤ったままで、Misannotation を外す条件は満たさない）。
候補付けビューアの「元の注釈に戻す」（`clear`）は、そのスポットの有効な判断が `wrong` でも `assign` でも `redundant` でも両方のタグを外す。
2 つのタグは 1 回の書き込みで変える（`msdial/tags.py` `update_alignment_tags`）。
他のタグと定義は残し、ファイルが無ければ MS-DIAL と同じ 5 定義で作る。形式は上流
`AlignmentResultContainer.Save` と同じ（`<Peak Id="<MasterAlignmentID>"><Tag>3</Tag></Peak>`、
UTF-8 BOM・CRLF）。

| フィールド | 意味 |
|---|---|
| `path` | 書いた（書こうとした）`_tags.xml` |
| `added` / `removed` | 実際に Misannotation が付いた／外れた `spot_id`（元から同じ状態のものは含まない） |
| `confirmed` | `{"added": [...], "removed": [...]}`。実際に Confirmed が付いた／外れた `spot_id` |
| `created` | ファイルを新しく作ったか |
| `backup` | 書く前の控え（`curation/tags-backup/<名前>.<UTC 時刻>`）。元のファイルが無ければ `null` |
| `note` | **MS-DIAL でプロジェクトを開いたままだと GUI の保存で上書きされる**旨。ユーザーに伝える |
| `error` / `message` | 反映に失敗したとき（壊れた XML・書けない等）。**フラグは記録済み**で、`_tags.xml` は触っていない |

MS-DIAL はアラインメントを保存するたびにメモリ上のタグで `_tags.xml` を丸ごと書き直し、
読むのはプロジェクトを開くときだけ。反映を見るにはプロジェクトを閉じてから開き直す。

### エクスポートのメタ行（`# curation = ...`）

`arf_export_differential` / `dataset_export_differential` は既定（`apply_curation=True`）で
有効な判断（`curation_flags` と同じ、スポットごとの最新 1 行・`clear` 済みは除く。`wrong` / `suspect` /
`assign` / `redundant`）を読み、差次的エクスポートの契約 15 列（`export_contract.EXPORT_COLUMNS`）自体は変えずに、
メタ行ブロックの `source_lines` スロットの**末尾**（`# source_arf`/`# source_mztab` 等、
経路固有のメタ行のすぐ後ろ。`export_contract.build_meta` の順序契約でスロット 3）に
1 行だけ足す（`metabolomix/curation/apply.py` の `meta_line()`）。**有効な判断が 0 件**なら
`meta_line()` は `None` を返し、この行自体を出さない——出力は現行と完全に同じになる。

タブ区切りの 1 行で、形は次のとおり（外側の角括弧内は `state` が `applied` のときだけ、内側の角括弧内はさらに
有効な `assign` / `redundant` が 1 件以上あるときだけ出る）:

```
# curation = <state>\tcuration_flags = N[\tcuration_wrong_excluded = N\tcuration_suspect = N[\tcuration_assigned = N\tcuration_redundant_excluded = N]]\tcuration_flags_sha256 = <hex>
```

| フィールド | 意味 |
|---|---|
| `<state>` | 下表 |
| `curation_flags = N` | 有効な判断の総数（`wrong` + `suspect` + `assign` + `redundant`。`flags_for_arf2()["n"]`） |
| `curation_wrong_excluded = N` | **`state == "applied"` のときだけ**出る。実際に出力の行から除外したスポット数 |
| `curation_suspect = N` | 同上。除外せず残した行のうち `suspect` フラグが付いているものの数（値そのものは変えていない） |
| `curation_assigned = N` | **`applied` かつ有効な `assign` / `redundant` があるときだけ**出る。`assign` で同定を置き換えて**出力に残った**行数 |
| `curation_redundant_excluded = N` | 同上。`redundant` で実際に出力の行から除外したスポット数 |
| `curation_flags_sha256 = <hex>` | 判断集合のダイジェスト（`flags_for_arf2()["digest"]`）。再エクスポートを跨いで判断の内容が変わっていないかを機械的に照合できる。`assign` は記録名・`level`・`inchikey`（同じ和組成でも別レコードへの付け替えを検出する）、`redundant` は相手と関係コードまで含む（`wrong` / `suspect` だけの集合のダイジェストは従来と同じ） |

`state` は経路とその引数で決まる:

| `state` | 経路 | 条件 |
|---|---|---|
| `applied` | ARF 経路（`arf_export_differential`）は `apply_curation=True`（既定）なら常にこれ。mzTab 経路（`dataset_export_differential`）は `apply_curation=True` かつ、隣接する `.arf` との対応が検証済み（`DatasetState.feature_qc["source"] == "arf"`、SMF_ID が `MasterAlignmentID` と同じ空間だと確認できている） | `wrong` と `redundant` のスポットを実際に出力から除外し、`assign` のスポットの同定を置き換え、`curation_wrong_excluded`/`curation_suspect`（と `assign` / `redundant` があれば `curation_assigned`/`curation_redundant_excluded`）を出す |
| `not_applied` | 両経路 | `apply_curation=False` を明示した。判断はあっても意図的に無視——行は 1 つも落とさず、同定も置き換えない |
| `unmapped` | mzTab 経路のみ | `apply_curation=True` だが SMF_ID と `MasterAlignmentID` の対応が未検証（`feature_qc["source"] != "arf"`）。**判断が 1 件でもあっても行は 1 つも落とさず、`assign` も当てない**——対応が確認できないまま除外・置き換えすると、mzTab の別 feature を `.arf2` のスポット番号と取り違えて黙って消す（または別の同定を付ける）しかねないため |

### エクスポートでの `assign` / `redundant` の扱い

両経路（`curation/apply.py` の 1 か所を共有）で、`state == "applied"` のとき:

- `wrong`: 除外（従来どおり）。`redundant`: 除外（`n_unannotated` に数える。「調べていない」行）。
- `assign`: 行の `name`・`ontology`・`inchikey` を記録した値で置き換え、`name_source` と `inchikey_source` を `curation` にする。
  ARF 経路の `msi_level` は記録名から求め直す（mzTab 経路は従来どおり空欄）。**未注釈だったスポットも InChIKey を得てエクスポートに新しく現れる**。
- **記録した InChIKey が空の `assign` は、その行を出力から落とす**（置き換えた MS-DIAL の InChIKey は残さない。InChIKey の無い行は
  そもそも書き出さないので、`n_unannotated` に数えられる）。ライブラリの参照レコードに InChIKey が無かった候補（`RIKEN … from …` 名の
  未知スペクトルなど）で起こる。
- 15 列と `CONTRACT_VERSION` は変わらない。判断が 0 件なら出力は従来と完全に同じ。

### エクスポートの戻り値の `curation`

両経路の成功 payload は同じ形の `curation` を持つ:
`{"state": ..., "wrong_excluded": N, "suspect": N, "orphaned": N}`。有効な `assign` / `redundant` があるときは
`assigned` と `redundant_excluded` が足される。

| キー | 意味 |
|---|---|
| `state` | 上表の `applied` / `not_applied` / `unmapped`。有効な判断が 0 件（メタ行を出さない）なら `null` |
| `wrong_excluded` / `suspect` | `state == "applied"` のときだけ数字（メタ行の同名フィールドと同じ値）。それ以外は `null`——メタ行と同じく「0 件除外した」と「除外していない」を混同させない |
| `assigned` / `redundant_excluded` | 有効な `assign` / `redundant` があるときだけ出るキー。`state == "applied"` のときだけ数字（メタ行の `curation_assigned` / `curation_redundant_excluded` と同じ値）、それ以外は `null` |
| `orphaned` | 同じファイル名で sha256 が違う（以前の版の）`.arf2` に付いた有効フラグの件数。当てないので行は落とさない。1 件以上なら payload の `warnings` にも同じ内容が出る。**メタ行には出さない** |

**パイプライン（`pipeline_run`）が書く差次的エクスポートはフラグを当てない。** キュレーションは
pipeline に入れていない（worker は session もフラグ記録も読まない）ので、フラグを反映した
ファイルが要るときは `dataset_export_differential` / `arf_export_differential` で書き出し直す。

### フラグ記録が壊れているとき

`flags.jsonl` に読めない行（途中で切れた JSON、`spot_id` が整数でない、`flag` が
`wrong`/`suspect`/`clear`/`assign`/`redundant` のどれでもない）があると、黙って読み飛ばさずに止める——
読み飛ばすとその行の `wrong` がエクスポートから外れないまま、誰も気づかないため。

| ツール | 振る舞い |
|---|---|
| `curation_review` / `curation_submit` / `curation_flags` | `{"status": "error", "message", "flags_file", "line"}`。`curation_submit` は壊れた記録に追記しない |
| `arf_export_differential` | 同上の形でエラー。ファイルは書かない |
| `dataset_export_differential` | `{"error": {"code": "CURATION_FLAGS_INVALID", "message", "details": {"flags_file", "line"}}}`。ファイルは書かない |
| `arf2_annotate_identities` | 一覧は返す。`curation_flag` 列を空にし、ヘッダに壊れた記録を名指しする注記を 1 行足す |

`line` は 1 始まりの行番号。その行を直すか消してから、もう一度実行する。

ARF 経路は `.arf2` の `MasterAlignmentID` をそのままキーに使うので対応の検証が要らず、
`unmapped` にはならない（`applied` か `not_applied` の 2 通りのみ）。`not_applied` と
`unmapped` はどちらも行を落とさないため、件数フィールド（`curation_wrong_excluded`/
`curation_suspect`）自体を省く——「0 件除外した」と「検証できないので除外していない」を
数字の 0 で混同させないため。

## 必須注意事項

1. **`UNKNOWN` は「不一致」ではない。** MS/MS 未取得・参照が引けない・参照に RT が無い等は
   「確かめられなかった」のであって「食い違った」のではない。`UNKNOWN` の系統は
   `verdict` の `FAIL`/`BORDERLINE` 集計に数えない（`core` トピック §9 の一般原則の
   キュレーション版）。
2. **RT–m/z 傾向（`trend_outlier`）は単独では `verdict` を動かさない。** ①〜④
   （msms/mz/rt/eic）の `BORDERLINE` と重なったときの補強材料に留める（judge.py の設計
   コメントのとおり）。傾向フィット自体の信頼度（`r2`）が低いクラスでは `trend_outlier`
   より `trend_outlier_unreliable`（情報止まり）に落ちる。
3. **Δ（`ppm`/`drt`）は代表試料（`RepresentativeFileID`）の行が基準。** 他サンプルの
   行の m/z・RT は使わない。EIC 系列は複数サンプルを見るが、Δppm/ΔRT は 1 行だけの値。
4. **`wdot` は MS-DIAL 自身の値の平方根であり、`-1` は「比較不能」（`library` トピック
   §14.2 と同じ番兵規約）。** TSV では `-1`（または欠測）はどちらも空セルにまとめて
   出すため、「値が無い」ことと「本当に一致度ゼロだった」ことを区別したいときは
   `curation_review` ではなく `library_match_feature` を直接使う。また `wdot`/`mpp` は
   `curation_review` が対向照合し直した値ではない——再照合の結果は理由コード
   `rescore_discrepancy`（両者が乖離したときだけ）にしか現れない。
5. **候補の並び順はスペクトル類似度と制約だけで決まり、MS-DIAL の脂質規則（診断イオン・鎖決定）は評価していない。**
   候補の上位は化学的にありえない鎖組成でもありうる。そのため既定の記録は和組成（`level = sum`）で、`assign` は
   「鎖組成まで確かめた」という意味を持たない。`likely_wrong` の判定に使う脂質規則フラグは MS-DIAL の照合結果を読んだものであって、
   こちらの再検索で出た新しい候補にはその裏付けが無い。
6. **④ の行の EIC の重ね描きで、相手 Y の代表試料のトレースは X の代表試料の最大値に合わせて縮めて描く**（形の比較用）。
   重なって見えても強度が同程度という意味ではない。強度の比は `isotope_ratio`（と `expected_ratio`）を見る。
7. **`redundant` は「X は Y のイオンだ」という判断であって、X が実在しないという判断ではない。** 同位体で強度比が期待の 1.5 倍を超えたもの
   （`informational`）は実在の別物質の可能性があるので、`redundant` には選べない。`found_in_upper_msms` リンクは向きの無い対なので、
   前駆体側は m/z の大小で決めている（軽い側が断片）。
8. **`total_score` の `-1` は「比較していない」であって「一致度ゼロ」ではない。** ビューアは `–` と表示する。
9. **`curation_suggest` の `likely_wrong` は元レビューの対象範囲についてしか分からない。** 元レビューを `ontology` などで絞っていた場合、
   範囲外の `likely_wrong` は拾えない（範囲外でも `wrong` フラグと未注釈は対象になる）。
