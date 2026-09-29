# アラインメントのキュレーション（`curation_review` 等）の出力フィールド

`curation_review` / `curation_submit` / `curation_flags` が返す値の意味。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。参照ライブラリ照合のスコアの意味（`-1`/`0` の区別など）は `library` トピックが定義し、ここでは繰り返さない。

ツールが**どのファイルのどの関数をどの順に呼ぶか**は `docs/workflow/curation.md`。
ここは**値の意味**だけを定義する。実装は `lipidmix/curation/`
（`judge.py` 機械判別・`evidence.py` 証拠収集・`eic_shape.py` EIC 形状・`trend.py` RT–m/z 傾向・
`review.py` レビュー生成/要約・`flags.py` フラグ永続化・`apply.py` エクスポート反映）。

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
| `html_path` | ビューア HTML のパス。ユーザーがブラウザで開き、フラグを付けて「送信用テキストをコピー」する。**全件**（`max_rows` で切らない）が入る。自動判別が `likely_wrong` のカードは赤の破線枠（未確認）で、ユーザーが「間違い」を付けると実線になる。クラスの選択肢には、各クラスの下に「<クラス> › 自動判別: 間違い」、先頭に全クラス横断の「自動判別: 間違い」が並ぶ（該当が 1 件以上あるときだけ） |
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
間引く（payload only。判定は変えない。定数・実装は `lipidmix/curation/evidence.py`）:

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

### `table`（TSV）の列

| 列 | 意味 |
|---|---|
| `spot_id` | `MasterAlignmentID` |
| `name` | ARF2 の `Name`（確度接頭辞を含む生の文字列。§9-1 の接頭辞規則は `core` トピック） |
| `ontology` | ARF2 の `Ontology`（脂質クラス） |
| `adduct` | ARF2 の `AdductType` |
| `verdict` | `likely_wrong` / `suspect` / `ok`（次節） |
| `reasons` | 立った理由コードをカンマ区切りで、強い→弱い→帯のみの順に並べたもの（次節の表と同じ順） |
| `ppm` | 代表試料の m/z と参照 precursor m/z（引けなければ Formula からの理論値）の相対誤差 [ppm] |
| `drt` | 代表試料の RT と参照 RT の差 [分]。参照が引けない、または参照に RT が無ければ空 |
| `wdot` | **MS-DIAL 自身が出した** `squared_weighted_dot_product` の平方根（`library` トピック §14.1 と同じ規約: 平方根側の値）。この照合結果は同定時点でアラインメントに保存済みのもので、`curation_review` がここで再照合した値ではない（後述の必須注意）。値が無い、または MS-DIAL 側が `-1`（比較不能）なら空 |
| `mpp` | **MS-DIAL 自身が出した** `matched_peaks_percentage`。`wdot` と同じく再照合値ではない |
| `eic_good` | 検出（非 gap-fill）サンプルのうち EIC 形状が「良い」と判定された割合（`eic_shape.spot_shape` の `good_fraction`）。良否は頂点がウィンドウ内・点数十分・理想ガウスとの R² 十分・極大数が上限以下で決まる。検出サンプルが無ければ空 |
| `trend_z` | 所属クラスの RT–m/z 傾向モデル（Huber 頑健回帰、RT = a + b·炭素数(+ c·二重結合数)）に対する頑健 z 値。傾向がフィットできないクラス（点数不足・炭素数のばらつき不足）では空 |
| `flag` | 既存のキュレーションフラグ（`wrong` / `suspect`）。未フラグは空。**`clear` はここに出ない**——`clear` は取り消しの記録で、有効フラグ（スポットごとの最新 1 行）が `clear` のスポットは未フラグとして扱うため |

### 判定（`verdict`）と理由コード

`verdict` は 5 系統の帯（`checks.{msms,mz,rt,eic,trend}`。各 `PASS`/`BORDERLINE`/`FAIL`/`UNKNOWN`
と理由コード）から決める:

- **`likely_wrong`**: 強い理由コードが 1 つでも立った場合。
- **`suspect`**: 強い理由が無く、かつ (a) msms/mz/rt/eic のいずれかが `FAIL`、または
  (b) msms/mz/rt/eic のうち `BORDERLINE` が 2 つ以上、または (c) `BORDERLINE` が 1 つ以上
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
| 強い | `class_rule_rejected` | MS/MS ありで、MS-DIAL の脂質クラス規則（診断イオン）を評価して棄却した（`is_lipid_class_match=False` かつ `is_other_lipid_match=False`）。実測では全件 `low score:` なので `low_score` も同時に立つ。**脂質規則が走ったデータに限る**（下記） |
| 弱い | `ppm_out` | Δppm が `ppm_borderline` しきい値（既定 10）を超えた。**adduct 非依存**（実測: kidney neg/pos で全 adduct の中央値が約 −0.8 ppm）で、単独では強い理由に数えない（ユーザー決定 2026-09-29）——mz 系統自体は `FAIL` になるが、単独では `suspect` 止まり。`polarity_mismatch`/`precursor_unmatched` が別途立てば、そちらの強さで `likely_wrong` になる |
| 弱い | `low_score` | MS/MS はあるが MS-DIAL 自身の `is_reference_matched` が `False` |
| 弱い | `drt_out` | ΔRT が `drt_borderline` しきい値（既定 1.0 分）を超えた |
| 弱い | `eic_poor` | EIC 形状帯が `FAIL`（検出サンプルが 1 件以上あり、かつ `good_fraction` が `eic_borderline_frac` 未満）。**検出サンプルが 0 件のときは `good_fraction` が計算できず帯は `UNKNOWN` になり、このコードは立たない**（`spot_shape` は `n_detected==0` なら `fraction=None`→`band="UNKNOWN"`） |
| 帯のみ | `ppm_borderline` | Δppm が `ppm_pass`〜`ppm_borderline`（既定 5〜10）の帯 |
| 帯のみ | `drt_borderline` | ΔRT が `drt_pass`〜`drt_borderline`（既定 0.5〜1.0 分）の帯 |
| 帯のみ | `eic_borderline` | EIC 形状帯が `BORDERLINE`（検出サンプルが 1 件以上あり、かつ `good_fraction` が `eic_borderline_frac`〜`eic_pass_frac` の帯）。これも検出 0 件では立たず `UNKNOWN` になる |
| 帯のみ | `rt_scatter` | 検出サンプル間で EIC 頂点 RT のばらつき（母標準偏差 `statistics.pstdev`。検出 2 件以上のときだけ計算）が `eic_rt_scatter_sd`（既定 0.1 分）を超えた。EIC 帯が `PASS` ならこの 1 件だけで `BORDERLINE` に格下げする |
| 帯のみ | `trend_outlier` | 所属クラスの RT–m/z 傾向で頑健 z が `trend_outlier_z`（既定 3.0）を超え、かつそのクラスの傾向フィット自体が信頼できる（`r2 >= trend_min_r2`、既定 0.7） |
| 情報 | `msms_absent` | MS/MS 未取得、または名前接頭辞が `no MS2`/`w/o MS2`（`msms` 系統は `UNKNOWN`） |
| 情報 | `reference_not_found` | ライブラリから参照レコードを引けなかった。**`rt` 系統だけが `UNKNOWN`** になる。`mz` 系統は Formula/AdductType からの理論値（`mass_error_ppm`、`ppm_basis="formula"`）にフォールバックして計算を続け、それも失敗したときだけ `UNKNOWN` になる（`rt` と違って自動的に `UNKNOWN` にはならない） |
| 情報 | `reference_rt_absent` | 参照は引けたが RT を持たない（`rt` 系統は `UNKNOWN`） |
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

**脂質規則フラグは脂質規則が走ったデータに限って読む。** 規則フラグは MS-DIAL の Lipidomics
採点器でしか立たず、それ以外の採点器では全部 `False` になる（`class_rule_rejected` の条件と
区別できない）。そこで 1 回のレビューの対象スポットに規則フラグ（class / chains / other）が
1 件でも `True` のものがあるときだけ、上の 3 コードを使う（`judge.lipid_rules_active`）。
1 件も無ければ `warnings` にその旨を出す。フラグの意味は上流 MsdialWorkbench `afd5f9522` の
`MsReferenceScorer` / `LipidMsmsCharacterization` と kidney の実測で確かめた。

### エクスポートのメタ行（`# curation = ...`）

`arf_export_differential` / `dataset_export_differential` は既定（`apply_curation=True`）で
有効フラグ（`curation_flags` と同じ、スポットごとの最新 1 行・`clear` 済みは除く）を読み、
差次的エクスポートの契約 15 列（`export_contract.EXPORT_COLUMNS`）自体は変えずに、
メタ行ブロックの `source_lines` スロットの**末尾**（`# source_arf`/`# source_mztab` 等、
経路固有のメタ行のすぐ後ろ。`export_contract.build_meta` の順序契約でスロット 3）に
1 行だけ足す（`lipidmix/curation/apply.py` の `meta_line()`）。**有効フラグが 0 件**なら
`meta_line()` は `None` を返し、この行自体を出さない——出力は現行と完全に同じになる。

タブ区切りの 1 行で、形は次のとおり（角括弧内は `state` が `applied` のときだけ出る）:

```
# curation = <state>\tcuration_flags = N[\tcuration_wrong_excluded = N\tcuration_suspect = N]\tcuration_flags_sha256 = <hex>
```

| フィールド | 意味 |
|---|---|
| `<state>` | 下表 |
| `curation_flags = N` | 有効フラグの総数（`wrong` + `suspect`。`flags_for_arf2()["n"]`） |
| `curation_wrong_excluded = N` | **`state == "applied"` のときだけ**出る。実際に出力の行から除外したスポット数 |
| `curation_suspect = N` | 同上。除外せず残した行のうち `suspect` フラグが付いているものの数（値そのものは変えていない） |
| `curation_flags_sha256 = <hex>` | フラグ集合のダイジェスト（`flags_for_arf2()["digest"]`）。再エクスポートを跨いでフラグ内容が変わっていないかを機械的に照合できる |

`state` は経路とその引数で決まる:

| `state` | 経路 | 条件 |
|---|---|---|
| `applied` | ARF 経路（`arf_export_differential`）は `apply_curation=True`（既定）なら常にこれ。mzTab 経路（`dataset_export_differential`）は `apply_curation=True` かつ、隣接する `.arf` との対応が検証済み（`DatasetState.feature_qc["source"] == "arf"`、SMF_ID が `MasterAlignmentID` と同じ空間だと確認できている） | `wrong` のスポットを実際に出力から除外し、`curation_wrong_excluded`/`curation_suspect` を出す |
| `not_applied` | 両経路 | `apply_curation=False` を明示した。フラグはあっても意図的に無視——行は 1 つも落とさない |
| `unmapped` | mzTab 経路のみ | `apply_curation=True` だが SMF_ID と `MasterAlignmentID` の対応が未検証（`feature_qc["source"] != "arf"`）。**フラグが 1 件でもあっても行は 1 つも落とさない**——対応が確認できないまま除外すると、mzTab の別 feature を `.arf2` のスポット番号と取り違えて黙って消しかねないため |

### エクスポートの戻り値の `curation`

両経路の成功 payload は同じ形の `curation` を持つ:
`{"state": ..., "wrong_excluded": N, "suspect": N, "orphaned": N}`。

| キー | 意味 |
|---|---|
| `state` | 上表の `applied` / `not_applied` / `unmapped`。有効フラグが 0 件（メタ行を出さない）なら `null` |
| `wrong_excluded` / `suspect` | `state == "applied"` のときだけ数字（メタ行の同名フィールドと同じ値）。それ以外は `null`——メタ行と同じく「0 件除外した」と「除外していない」を混同させない |
| `orphaned` | 同じファイル名で sha256 が違う（以前の版の）`.arf2` に付いた有効フラグの件数。当てないので行は落とさない。1 件以上なら payload の `warnings` にも同じ内容が出る。**メタ行には出さない** |

**パイプライン（`pipeline_run`）が書く差次的エクスポートはフラグを当てない。** キュレーションは
pipeline に入れていない（worker は session もフラグ記録も読まない）ので、フラグを反映した
ファイルが要るときは `dataset_export_differential` / `arf_export_differential` で書き出し直す。

### フラグ記録が壊れているとき

`flags.jsonl` に読めない行（途中で切れた JSON、`spot_id` が整数でない、`flag` が
`wrong`/`suspect`/`clear` のどれでもない）があると、黙って読み飛ばさずに止める——
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
