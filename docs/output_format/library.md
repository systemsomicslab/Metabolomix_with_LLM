# 参照ライブラリ照合（MS/MS スペクトル照合）の出力フィールド

`library_load` / `library_match_feature` / `library_plot_mirror` / `library_search` が返す値の意味。呼ぶ順は `docs/workflow/library.md`。

> 先に `lipidmix://docs/output-format`（共通核）を読む。節番号は全トピックの通し番号。

## 14. 参照ライブラリ照合（`metabolomix/library/`, `metabolomix/analysis/spectral_match.py`）

`.dcl` の測定 MS/MS を参照ライブラリ（`*_Loaded.msp2.dbs` 優先、無ければ `*.msp`）と突き合わせ、MS-DIAL の個別スコアの定義を移植して数値を出す。候補の並び順も MS-DIAL の `TotalScore` に合わせるが、上流の順位付けタプル全体は再現しない（§14.9）。

### 14.1 各スコアの意味

`library_match_feature` の候補一覧・`library_plot_mirror` の材料は同じ指標を使う。

| キー | MS-DIAL 対応 | 意味 |
|---|---|---|
| `simple_dot_product` | `GetSimpleDotProduct` | m/z 重みなし・ペナルティなしの単純ドット積 |
| `weighted_dot_product` | `GetWeightedDotProduct` | m/z 重みとピーク数ペナルティを掛けたドット積。mzTab の `id_confidence_measure[5]`。**単独では並び順に使わない**（§14.9） |
| `reverse_dot_product` | `GetReverseDotProduct` | 参照側のピークだけを基準に測定側で説明できているかを問う非対称スコア |
| `matched_peaks_percentage` | `GetMatchedPeaksScores`（比） | 参照ピークのうち測定側でも観測された割合 |
| `matched_peaks_count` | 同上（count） | 一致した参照ピークの本数。**前駆体（precursor）の一致も 1 本に数える**（前駆体しか合っていない候補でも 1）。断片だけの本数は `detail=true` の `fragment_matches_excl_precursor` |
| `entropy_similarity` | `GetSpectralEntropySimilarity` | スペクトルエントロピーに基づく類似度。**mzTab には出ない** |
| `total_score` | `GetTotalScore` | 上の指標に RT・precursor m/z の一致度を足した**候補の並び順（`ranked_by`）**。**0〜1 ではない**（§14.9） |

**`simple` / `weighted` / `reverse` の 3 つは平方根側の値**（`match_spectrum` が二乗値から `sqrt` を取る。モジュール内の単体関数 `weighted_dot_product()` 等は二乗値を返す）。mzTab-M の `id_confidence_measure[4..6]` と同じ尺度で、`0.0`〜`1.0` の類似度として読める。

### 14.2 `-1` と `0` の区別（最重要）

**`-1` は「比較していない」、`0` は「合わなかった」。**

- `-1`（3 つの dot product・`matched_peaks_percentage`・`matched_peaks_count`・`entropy_similarity`）: 測定か参照のスペクトルが空で比較できなかった。`entropy_similarity` の `-1` は `sqrt` を経ずに `candidates_table` に `-1.0` で出る。
- `0`: 比較できたが一致度がゼロ。dot product 3 種は窓合算後に信号が無かったときも `0.0`。`entropy_similarity` は、両方非空でもどちらかの総強度が 0（例: `[[100.0, 0.0]]`）なら `0.0`。

mzTab-M の値とは番兵が違う:

- `id_confidence_measure[4..6]`（simple/weighted/reverse）は、MS-DIAL が `-1` を `0` にクランプしてから `sqrt` を取る（`Math.Sqrt(Math.Max(Squared*, 0f))`）ので、比較不能は `0` で出る。比較不能に限り、こちらの `-1` と mzTab の `0` は同じ。
- `id_confidence_measure[7..8]`（matched peaks count / percentage）はクランプせず `-1` のまま出る。

### 14.3 意図的に写した瑕疵（`metabolomix/analysis/spectral_match.py`）

MS-DIAL との数値の一致のため、上流の次の 5 点をそのまま残している（直すと数値が変わる）。

1. weighted/reverse の `wM` / `wR` は上流でも使われない。移植していない。
2. weighted の `if (sumM <= 0 && sumR > 0) {...} else {...}` は両枝が同一処理。条件ごと落として共通処理だけ残す。
3. simple の `× 999` は比で打ち消えるので数値に影響しない。
4. 窓の走査は固定幅グリッドではなくカーソル走査で、同じピークが隣り合う 2 つの窓に二重計上されうる。
5. entropy に Li et al. 2021 の低エントロピー重み変換は入っていない（上流にも無い）。

### 14.4 `.msp` 由来には許容幅が同梱されない

`.dbs` の `Storage` エントリには、その run が使った許容幅（`mz_tol` / `ms2_tol` / `rt_tol` 等。`docs/schema/molecule_ms_reference.md` の `MsRefSearchParameterBase`）が入っている。`library_load` はこれを store のメタの `search_params` に保存し、`library_match_feature` は 明示指定 > `search_params` > 既定値 の順で使う。

**`.msp` にはこの情報が無い。** `.msp` 由来の store は `search_params=None` で、`library_load` の戻り値に既定値を使う旨の `note` が付く。既定値は `metabolomix/library/defaults.py` の `DEFAULT_MZ_TOL` (0.01) / `DEFAULT_MS2_TOL` (0.025) / `DEFAULT_RT_TOL` (0.2)。`DEFAULT_RT_TOL` は上流 `RtTolerance` の既定（100.0）ではなく、`.dcl` 検索（`dcl/reader.py` / `dcl/tools.py`）の RT 窓の既定と同じ値。

**ライブラリ照合の窓と `.dcl` 検索の窓は別。** `library_match_feature` の `mz_tol` / `rt_tol` と `search_params` 由来の値は、ライブラリ候補の検索（`store.candidates()`）にだけ効く。`.dcl` から測定 MS/MS を引く窓（`_measured_spectrum`）は常に `dcl/reader.py` と同じ `mz_tol=0.01` / `rt_tol=0.2` で、`search_params` を見ない（`.dbs` の `RtTolerance` は RT フィルタを実質無効にする 100.0 が典型で、流用すると別ピークのスペクトルを拾う）。RT を渡したとき `.dcl` で複数ヒットすれば、RT が最も近いものを選ぶ。

### 14.5 `library_load`

`.dbs`/`.msp` を解決して SQLite store を構築（または既存キャッシュを再利用）した要約を返す。

| キー | 意味 |
|---|---|
| `record_count` | store に格納した参照レコード件数 |
| `ion_modes` | `{"positive": N, "negative": N}` のような件数内訳 |
| `compound_classes` | 化合物クラス別件数の上位 10（`compound_class IS NULL` は除外） |
| `search_params` | `.dbs` 由来なら実測の許容幅と注釈スコアリングのフラグ、`.msp` 由来なら `null`（§14.4） |
| `source_sha256` | 元ファイルの sha256（store のキャッシュキーと同じ） |
| `skipped_no_precursor_mz` | precursor m/z が無い（または解釈できない）ため読み飛ばしたレコード件数（`record_count` に含まれない。公開 `.msp` に普通に混ざる）。0 より多ければ `note` にも件数が出る |
| `records_without_ion_mode` | 極性（`IONMODE`）を持たないレコード件数（`record_count` に含まれる）。`library_match_feature` の `ion_mode` 絞り込みはこれらを**極性不明として候補に残す**ので、極性違いのファイルを読んでも候補が残る。測定の極性に合ったファイルかは呼び出し側が確かめる（0 より多ければ `note` にも出る） |
| `non_utf8_lines` | UTF-8 で読めず cp932（無理なら latin-1）で読んだ `.msp` の行数。0 より多ければ化合物名が化けている可能性がある |

**`file_path` を省いたとき（記憶から開く）**: 直近のデータセット（`dataset_handle` を渡せばそのハンドル）が同定に使ったライブラリの記憶（下の「データセットが覚えているライブラリ（LibraryRecord）」）が `known` で、sha256 を測り直して通れば、その実体を開く（戻り値に `library_source: "dataset_record"`）。**記憶以外には落ちない**（データフォルダの探索・設定の極性の MSP・`ion_mode` 引数は無い）。記憶が `unknown` か無ければ `LIBRARY_NOT_REMEMBERED`、`mismatch` / `lbm_only` は `LIBRARY_RECORD_MISMATCH` / `LIBRARY_RECORD_LBM_ONLY`（理由つき。`file_path` で明示すれば開ける）、データセットが無ければ `missing_state`（`load_dataset`）。

`file_path` で開いたものが記憶と違うライブラリ（中身の sha256 で比べる。同じ内容のコピーは違わない）なら、そのまま開いて `differs_from_record: {"remembered": <記憶のファイル名>, "opened": <開いたファイル名>}` を足す。記憶が `lbm_only` のまま別のライブラリを明示して開いたときは、同定に使ったものではない旨が `note` に付く。

**照合側の LBM の扱い**: `library_match_feature` / `library_msms_qualify` は、ライブラリが開かれておらず記憶が `lbm_only` のとき、`missing_state` ではなく `{"status": "not_evaluable", "reason": ...}`（qualify は `verdict` も `not_evaluable`）を返す。記憶が `lbm_also` で、`spot_id` で照合した特徴の代表の注釈器が LBM 由来（`is_lbm_annotator` が True）なら、結果に `identified_by_lbm: true` を足す（そうでなければキーが無い。`false` は書かない）。この印は解析フォルダ（`.arf2`）から読んだデータセットにだけ付く（`annotator_id` は `.arf2` の値で、mzTab-M のデータセットは持たない）。印が無いことは LBM 由来でない証拠にならない。

**ライブラリの示し方**: 戻り値はライブラリをいつもファイル名（`file`）と中身の sha256（`source_sha256`）で示し、置き場所（パス）は出さない（伏せる切り替えの引数は無い）。`load_dataset` の戻り値の `library` / `library_record` も同じ（[reanalysis.md](reanalysis.md) §1 の `load_dataset` の表）。

**エラー戻り値**: ライブラリを 1 つに決められないときは `{"status": "error", "code": ..., "message": ..., "details": ...}` を返し、store は差し替えない。`code`: `LIBRARY_NOT_REMEMBERED`（`file_path` を省いたが、データセットが同定に使ったライブラリを覚えていない。候補のファイル名と「`file_path` を渡す」の案内）/ `LIBRARY_RECORD_MISMATCH` / `LIBRARY_RECORD_LBM_ONLY` / `LIBRARY_NOT_FOUND`（明示した `file_path` が無い）/ `CONFIG_INVALID`（設定ファイルが読めない。`details` に `config_file` `line` `column`）。`message` と `details` にはファイル名・設定キー・設定ファイルのパスだけが入り、ライブラリの置き場所は入らない。

### 14.6 `library_match_feature`

照合対象は `precursor_mz` か `spot_id` のどちらか一方で指定する。

**候補ごとの質量誤差**: `candidates_table` の各行に `delta_mz`（= 測定 − 参照、Da、符号付き）と `delta_ppm`（= `(測定 − 参照) / 参照 * 1e6`、測定が高いと正）が付く。`precursor_mz` 列は**参照**レコードの前駆体で、6 桁で出す。

**前駆体の窓（MS-DIAL と同じ）**: `mz_tol`（明示 → `.dbs` の `Ms1Tolerance` → 既定 0.01 Da）は m/z 500 以下での Da の幅で、m/z 500 を超えたら同じ ppm に広げる（上流 `FixMassTolerance`。既定なら `max(0.01 Da, m/z × 20 ppm)`、m/z 760 で約 0.0152 Da）。実際に使った Da の幅は `query.mz_window`（`query.mz_tol` は広げる前の値）。低い m/z では 0.01 Da が数十 ppm に当たり +44 ppm の候補も残るので、受け入れる前に `delta_ppm` を見る。`mz_tol_ppm`（任意）は窓に重ねる追加の絞り込みで、`|delta_ppm|` がこれを超える候補を落とす（窓は広げない。`mz_tol_ppm=10` なら +44 ppm の候補は落ちる）。`query.mz_tol_ppm` に指定値（省略時 `null`）が残る。

**InChIKey の列**: `formula` の後に `inchikey`（ライブラリの値のまま。無ければ空）と `inchikey_standard`（第 2 ブロックの 9 文字目が `S` なら `true`、`N` なら `false`、InChIKey の形でなければ空）。同じ名前の立体異性体を見分けるのに使う。外部（PubChem など）へは照会しない。

**`ion_mode`**: `positive` / `negative` のほか `pos` / `neg` / `+` / `-`（大小を問わない。`core.polarity_names.normalize_polarity`）を読む。読めない値は `invalid_arguments`。`library_search` も同じ。store の比較（`store.candidates()`）も大小を問わない（`COLLATE NOCASE`。`"positive"` と `"Positive"` は一致する）。

| `status` | 意味 |
|---|---|
| `not_found` | `.dcl` に該当 precursor（`spot_id` のときはそのスポット）の MS/MS が無い（**未取得**であって「合わなかった」ではない） |
| `no_candidates` | 測定 MS/MS はあるが、m/z 窓・極性・RT に該当する参照レコードが無い |
| `success` | 候補を採点して `candidates_table`（TSV、`rank`/`name`/... 列。§14.1 の指標を含む）を返した |
| `spot_not_found` | `spot_id` が `.arf2` に無い |
| `dcl_precursor_mismatch` | アラインメント `.dcl` の `spot_id` 番の precursor・RT が ARF2 の `MassCenter`・`RT` と ±0.05 Da・±0.2 分で合わない（索引対応の崩れ。`.dcl` の値は代表試料の実測値なので完全一致はしない）。照合しない |
| `invalid_arguments` | `precursor_mz` と `spot_id` を両方渡した、またはどちらも渡さなかった。`ion_mode` が読めない・`mz_tol_ppm` が正の有限の数でない |

照合を試みて `success` 以外（`not_found` / `no_candidates` / `spot_not_found` / `dcl_precursor_mismatch` / `.arf2`・`.dcl` が無い）に終わったら、直前の照合（`last_match`）を捨てる（`library_plot_mirror` が前の照合を描かない）。`invalid_arguments` では照合を試みないので残る。

**`spot_id` 経由の照合**: precursor m/z・RT・極性（`IonMode` を正規化。`Both` / `Unknown` なら極性で絞らない）は同じバッチの `.arf2` から取る（`rt` / `ion_mode` を明示すればそちらが優先）。測定スペクトルはアラインメント `.dcl` の `dcl[spot_id]`＝代表試料の MS/MS（`curation` と同じ対応）で、precursor m/z による `.dcl` 検索は通らない（同じ m/z の別スポットを拾わない）。`query.spot_id` と `session.library.last_match["spot_id"]` に照合したスポットが残る（precursor 経由では `null`）。

**MS-DIAL の注釈との食い違い**: `spot_id` で照合したスポットに MS-DIAL の注釈（`.arf2` の代表の `Name`。`Unknown` や空は注釈なし）があれば、戻り値に 3 つの欄を足す（注釈が無ければキーごと無い）。

| 欄 | 意味 |
|---|---|
| `msdial_name` | MS-DIAL の注釈の名前（`\|` があればその後の区間） |
| `msdial_rank` | 注釈と同じ物質の候補の最上位の順位（全候補の中。`top_n` で切る前）。候補に無ければ `null` |
| `top_differs_from_msdial` | 1 位の候補が注釈と別の物質か（`msdial_rank` が 1 でなければ `true`） |

同じ物質かは、InChIKey の先頭ブロック（両方にあれば）→ 脂質の種（両方が読めれば。`msms_rules.record_lipid_species` と同じ読み方）→ 名前（限定子を外し、大小を問わない）の順に比べる。`true` は「MS-DIAL が誤った」ではない（こちらの順位は `total_score` 単独で上流と違う。§14.9）。対向プロットで両方を見る。

候補は `total_score` の降順（`ranked_by` が基準を示す。限界は §14.9）。スペクトル座標・alignment は戻り値に含めず `session.library.last_match` に持つ（`library_plot_mirror` がそこから読む）。

`scoring` は総合スコアの組み立て（`rule` / `use_rt` / `ms1_tol` / `rt_tol`）を 1 回だけ載せる（RT 項が入ったかは候補の数値からは分からない）。総合スコアの内訳（`rt_similarity` / `mass_similarity` / `spectrum_score`）は候補表に出さず、`session.library.last_match` の各候補の `scores` にだけ持つ。

### 14.7 `library_plot_mirror`

**何を描くか**: 既定は直近の `library_match_feature` の照合（`last_match`）。`spot_id` を渡すと、そのスポットを `library_match_feature(spot_id=…)` と同じ規則（許容幅は `.dbs` の `search_params` か既定）で照合し直し、全候補を `last_match` に持ってから描く。照合に失敗したら直前の照合を捨て、その封筒（`spot_not_found` / `not_found` / `no_candidates` など）を返す。`candidate`（任意）は描く候補を InChIKey（全体か先頭 14 文字）か名前（大小を問わない完全一致）で選び、`rank` より優先する（照合の候補に無ければ `status: "error"`）。**全体の InChIKey は全体で一致する候補を選ぶ**（先頭ブロックが同じ別の立体異性体を選ばない）。全体で合う候補が無いときだけ、先頭ブロックが同じ最上位の候補を選び、caption に 1 文を足す。

どの図かは、画像は caption（`Mirror: spot <id> (precursor m/z=…) vs <名前> (rank <n>, InChIKey=…, precursor m/z=<参照>)`。precursor で照合したときは `measured`）、`output="payload"` は `match`（`spot_id`・`precursor_mz`・`rank`・`name`・`inchikey`・`reference_precursor_mz`・`picked_by`）に出る。`picked_by`: `rank`（`candidate` を省いた）/ `inchikey`（全体の一致）/ `inchikey_block`（先頭ブロックの一致。全体の InChIKey を渡して全体で合わなかったときも）/ `name`。

**座標契約** `lipidmix.mirror.v2`（`build_mirror_payload` / `render_mirror`、`metabolomix/plots/mirror.py`）: 上段が測定（上向き）、下段が参照（下向き）、横軸 m/z 共通。`matched_mz` は**常に参照側の m/z**（参照グリッドの窓中心）。測定側の一致 m/z は `matched_measured_mz` で、`ms2_tol` を渡したときだけ計算する（渡さなければ空で、測定側は一致の色分けをしない）。

**採点に入らなかった測定ピーク（`unscored_mz`）**: 採点は `.dbs` の `relative_amp_cutoff` / `absolute_amp_cutoff` で足切りした後のスペクトル（`normalize_measured`）に対して行い、図は足切り前の生ピークを描く。足切りで落ちたピーク（判定は `spectral_match.cutoff_mask`）は消さずに灰色で薄く描き、ラベルを付けない。凡例の `Below cutoff (not scored)` と画像 caption の件数は該当ピークがあるときだけ出る。`payload` の `scored_peak_count` / `unscored_peak_count` は常に載る。足切りが 0（`MsRefSearchParameterBase` の Key 7/8 の既定 0.0）の run ではこの層は出ない。

ピークは素のステムだけで、先端にマーカーは打たない（上流 `LineSpectrumControlSlim` と同じ）。

`label_policy` は m/z ラベルの衝突回避の方式。どちらも上流 `Annotator.OnRender` と同じ貪欲法（強度降順に走査し、既に置いたラベルと重なるものを飛ばす）で、重なりの見方が違う:

| | `"auto"`（既定） | `"msdial"` |
|---|---|---|
| 判定 | 水平と垂直の両方が近いときだけ飛ばす | **水平距離だけ**見る |
| 箱の大きさ | ラベル自身の文字列を実測 | 代表文字列 `"1000.00000"` 1 つで全ラベル共通 |
| 結果 | ラベルが多く載る | 少なくなる |

`"msdial"` が水平だけを見るのは上流の縮退を写したもの（`Overlap="Horizontal, Direct"` を全ラベル共通の代表箱で判定すると `Horizontal` 単独と同じになる）。縦に離れていても m/z が近ければ飛ばす。どちらも**上下は別枠**で数える。片側 25 本の上限は安全弁で、`"msdial"` では水平棄却が先に効いて届かない。

**画像のラベルはスペクトル全点から選ぶ**（衝突判定は描画時）。`payload["labels"]` は自前で描くクライアント向けの要約で、`top_labels` は**側ごと**の上限（`labels` の長さは最大 `2 × top_labels`）。

`scale`（`"relative"` 既定 / `"sqrt"` / `"log10"`）は縦軸の写し方で、上下それぞれ自分の最大値で正規化したあとに掛ける。**precursor がベースピークのスペクトル**（脂質の `[M-H]-` など）は `relative` だと診断イオンが相対数 % に潰れるので `"sqrt"` を使う（大小関係は保つ）。`log10` は相対 0.1%（3 桁）を下限にする（こちら独自の値。上流の Log10 軸の下限 `lowBound=1` は生強度 1 カウントで、同じ値ではない）。軸ラベルに選んだスケールが出る。**`output="payload"` の座標は正規化前の生値**で `scale` の影響を受けない。`Absolute` 軸は無い（単位の違う 2 つのスペクトルを生強度で並べても比較にならない）。

### 14.8 （撤去）ピーク裏取りツールへの統合

旧ピーク裏取りツールの `analytical_checks.msms.spectral_match` ブロックは無い。ピークごとの裏取りは `dcl_spot_spectrum`（スペクトルを読む）と `library_match_feature`（候補ごとの照合。「最良」は `total_score` 最大）を組み合わせる。

### 14.9 候補の並び順（`total_score`）と、その限界

`total_score` は上流 `MsScanMatching.GetTotalScore` の移植で、**正規化しない和**:

```
RtSimilarity + AcurateMassSimilarity + (Weighted + Simple + Reverse)/3
+ MatchedPeaksPercentage
```

**各項は `> 0` のときだけ加算する**（`-1` の番兵を足さない。§14.2）。**上限は 1 ではなく項数（この経路では最大 4）**で、0〜1 の類似度として読まない。RT と precursor m/z の一致度は `exp(-0.5*((actual-reference)/tolerance)^2)`（`gaussian_similarity`）で、m/z の許容幅は 500 を超えると ppm 換算で伸びる（`fix_mass_tolerance`）。

**RT 項が入る条件**: store の `search_params["use_time_for_annotation_scoring"]`（`.dbs` の `MsRefSearchParameterBase` `Key(16)` = `IsUseTimeForAnnotationScoring`。上流の既定は `False`）が真で、測定と参照の RT が両方あるときだけ。`.msp` 由来のライブラリでは入らない。結果は `library_match_feature` の `scoring.use_rt`。

**`total_score` の絶対値は MS-DIAL の `TotalScore` と比べられない。** 上流の脂質経路（`GetLipidomicsMatchedPeaksScores`）は、脂質クラスの注釈レベルが立つと matched peaks の項を固定のボーナス値（`2.0` など。1 を超える）で上書きする。こちらは汎用の割合（≤ 1.0）を足すので、確度の高い脂質では上流より約 1.0 低い。順位付けには支障ないが、mzTab の `id_confidence_measure[1]`（MS-DIAL algorithm matching score）と数値を突き合わせない。`[8]` の一致率が 68% に留まるのもこの仕様差による。

**上流の順位付けとは同じでない。** 上流 `MsScanMatchResultContainer.ResultOrder` は `(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)` の辞書順で、`TotalScore` は最後のタイブレーク。`IsReferenceMatched` は脂質クラス固有の判定（`MsReferenceScorer.ValidateOnLipidomics` → `GetRefinedLipidAnnotationLevel`）に依存するので移植していない（この判定が無いと実データで全候補が `IsReferenceMatched=True` になり、ゲートとして働かない）。

そのため**上位に化学的にありえない候補が残ることがある**（例: PI 18:0_20:4 の照合で `SMGDG 18:0_20:4` が 2 位、PC 18:0_22:6 の照合で `SHexCer 41:1;O3` が 2 位）。1 位の妥当性は対向プロット（`library_plot_mirror`）と `matched_peaks_count` で確かめる。

実データ（aging mice kidney neg、120 特徴、正解は mzTab の `chemical_name`）の top-1 一致は `weighted_dot_product` 単独 106/120（88.3%）、`total_score` 113/120（94.2%）。変わった 7 件はすべて `total_score` 側が正解だった。

### 14.10 MS-DIAL との一致点・相違点（まとめ）

MS-DIAL の採点・順位付けは GUI 固有ではなく共有ライブラリ（`MsScanMatching`・`MsReferenceScorer`・`MsScanMatchResultContainer`・`MztabFormatExport`）にあり、GUI と Console で同じ。描画は GUI だけにある（Console は図を描かない）。

#### ロジック

同じもの: 個別スコア 5 種（瑕疵ごと。実データで `[4][5][6][7]` が 194/194 一致。§14.3）・採点前処理（`relative_amp_cutoff` / `absolute_amp_cutoff` / 質量範囲を `.dbs` の実値から）・総合スコアの式と RT / m/z の一致度・RT 項を足す条件（§14.9）・候補検索の窓（precursor m/z ± `Ms1Tolerance`〔m/z 500 超は同じ ppm〕、RT ± `RtTolerance`）。

| | こちら | MS-DIAL | 備考 |
|---|---|---|---|
| **並び順** | `total_score` 単独 | `(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)` の辞書順 | §14.9 |
| **候補の足切り** | しない（採点した全件を返す。`library_msms_qualify` は `.dbs` の足切りを通るかを `passes_dbs_cutoffs` の列で示すが、候補は消さない） | `GetRefinedLipidAnnotationLevel` が空文字を返す候補は消える。各種 cutoff・`IsSpectrumMatch` でも絞る | 上位に化学的にありえない候補が残る直接の原因 |
| **名前** | ライブラリのレコード名のまま | 注釈レベルに応じて付け直す（クラス止まり / 鎖組成まで） | |
| **`matched_peaks` の 2 値** | 汎用 `GetMatchedPeaksScores` | 脂質は `GetLipidomicsMatchedPeaksScores`（ボーナス値で上書き） | §14.9 |
| **足さない項** | RT・m/z・スペクトル・matched% の 4 項のみ | CCS / isotope / Andromeda も足しうる | この経路（LC-MS、MS2、脂質）には無い項 |
| **採点対象のファイル** | 呼び出し側が指定した 1 つの `.dcl` | アラインメント spot の「代表ファイル」1 つ | mzTab の値と 0.3〜2% ずれる主因 |

#### 描画（対向プロット）

同じもの: 上=測定 / 下=参照で m/z 軸共有・上下それぞれ自分の最大値で正規化・素のステム（先端マーカー無し）・m/z ラベルは小数 4 桁で強度降順（上流 `Format="F4"`）・ラベルの衝突回避の貪欲法。

| | こちら | MS-DIAL |
|---|---|---|
| 縦軸スケール | `relative` / `sqrt` / `log10`（`Absolute` 無し） | Relative / Absolute / Log10 / Sqrt（上下独立） |
| **既定の衝突判定** | 水平＋垂直（`"auto"`）。忠実版は `label_policy="msdial"` | 水平のみ |
| ラベルの本数上限 | 片側 25 本（安全弁） | 指定なし（MS2 ビューは `TopN` なし） |
| **一致ピークの色** | 緑でハイライト＋ガイド線 | 概念が無い（色は `SpectrumPeak.SpectrumComment` 由来。既定は側ごとに 1 色） |
| **採点対象外ピーク** | 灰色で薄く描き、凡例と caption に件数（足切りが非ゼロのときだけ。§14.7） | 概念が無い（足切り後のスペクトルを描く） |
| 色 | `#2471a3` / `#c0392b`（volcano と同系統） | 純 Blue / 純 Red |
| 差分・積の重ね描き | 無い | ある（`UpperDifferenceSpectrumModel` / `UpperProductSpectrumModel`） |
| 対話操作 | 無い（静止 PNG、`output="payload"` で座標） | ツールチップ（m/z・強度・`SpectrumComment`）・ズーム |

### 照合の詳細（`detail=true`）とライブラリの検索

`library_match_feature(detail=true)`（opt-in）は、既存の列の並びを変えずに TSV の末尾へ 1 列を足し、payload に 4 つの欄を足す。`detail` を渡さない戻り値は変わらない。

| 欄 | 意味 |
|---|---|
| `fragment_matches_excl_precursor`（列） | 一致した参照ピーク（`matched_peaks_count` と同じ規則）のうち、前駆体 m/z ± `precursor_exclusion_tol` を除き、測定強度が基準ピークの `min_fragment_relative_intensity` 以上の窓の本数。下限の既定は MS/MS の合否の規則と同じ 0.01 なので、既定のままなら `library_msms_qualify` / `dataset_claim_evidence` の同じ名前の列と同じ本数になる（reanalysis.md §4）。比較できない（測定か参照が空）なら空 |
| `precursor_exclusion_tol` | 前駆体を除く幅（Da）。省略時は採点に使った `ms2_tol` |
| `min_fragment_relative_intensity` | 上の列に使った測定強度の下限（基準ピークに対する比）。引数 `min_fragment_relative_intensity`（0 以上 1 以下、`0` で下限なし）で変えられ、省略時は 0.01 |
| `measured_peaks` | 照合に使った測定スペクトル（全ピーク `[m/z, 強度]`） |
| `alignments` | 候補ごと（`rank` 順）の一致したピーク `[m/z, 測定の強度, 参照の強度]`。強度は採点の前処理（足切り・正規化）の後の値 |

`library_search` は読み込んだライブラリのレコードを名前（大小を問わない部分一致）・分子式（完全一致）・InChIKey（先頭 14 文字）で引く。スペクトルは返さず、ピーク数 `n_peaks` だけ。`dcl_spot_spectrum`（`dcl` トピック）は spot_id で代表 MS/MS の全ピークを返す。

### データセットが覚えているライブラリ（LibraryRecord）

データセットは同定に使われた参照ライブラリを 1 件覚える（`metabolomix.dataset.library_memory`）。後段の MS/MS 照合はライブラリを名指しせずにこれを引く。戻り値にはファイル名と sha256 だけを出し、置き場所（`path`）は出さない。

解析フォルダで読んだとき、記憶はそのフォルダの**今の** `*_Loaded.*.dbs`（MS-DIAL が保存した写し）に従い、読み込みのたびに測り直す（dbs が書き換わっていれば、読み直した記憶は新しい中身に付く）。`mismatch` になるのは、読み込み済みのデータセットの記憶と今のファイルが食い違うとき、pipeline の run が固定した sha256 と今のファイルが違うとき、ジョブの実効メソッドが宣言するファイルが無い・読めないとき。

| 欄 | 意味 |
|---|---|
| `status` | `known`（1 本に決まった）/ `unknown`（決められない。`candidates` に候補のファイル名）/ `mismatch`（記録した sha256 と今のファイルが違う、またはファイルが無い）/ `lbm_only`（LBM だけで MSP は使っていない） |
| `kind` | `dbs`（解析フォルダの `*_Loaded.*.dbs`）/ `msp` / `lbm` / 無し |
| `origin` | `loaded_dbs`（フォルダ）/ `effective_method`（ジョブの実効メソッド。sha256 は読んだ時点で測る）/ `run_inputs`（pipeline の run が固定した入力。固定時の sha256 と今を突き合わせる）/ 無し |
| `lbm_also` | MSP に加えて LBM も照合に使った（`library_mode="msp_then_lbm"`、フォルダに LBM の dbs もある、またはジョブの実効メソッドが MSP と LBM の両方を宣言している） |

- 解析フォルダに非 LBM の `*_Loaded.*.dbs` が 2 本以上あると `unknown`。LBM の dbs（`.lbm2.dbs` など）は数えず、`lbm_also` にする。
- LBM の dbs だけが 2 本以上あるときは `lbm_only` で `file` / `path` は無し、`candidates` に LBM のファイル名。
- ジョブの実効メソッドが宣言するファイルが無い・読めないときは `mismatch`。
- mzTab-M を直接読んだときは常に `unknown`（mzTab-M は使ったライブラリを持たない）。
- 入口ごとの出どころ（`load_dataset`）: `directory` はフォルダの `*_Loaded.*.dbs`、`job_path` はジョブの実効メソッド、`pipeline_path` は run が固定した入力（どちらのアダプタで読んでも同じ）、`mztab_path` は `unknown`。`load_dataset` の `library_record` がこの記録を（`path` を外して）返す。

#### LBM が代表になった特徴の印（`is_lbm_annotator`）

`.arf2` の代表の `annotator_id` から LBM 由来かを判定する。Console の注釈器キーは、LBM なら常にライブラリのパス（末尾が `.lbm` + 数字）、MSP は `library_mode="msp_only"` で `Msp file path` のパス、`msp_then_lbm` では設定表の名前（pipeline の既定は `LabMsp`）。

- 末尾（ファイル名の部分）が `.lbm` + 数字なら **True**（LBM 由来）、`.msp` + 数字（`msp_only` の MSP のパス）なら **False**（MSP 由来）、それ以外は **None**（判定できない）。GUI の注釈器 ID は自由名で、`msp_then_lbm` の設定表の名前も None。None は MSP 由来の証拠にならない。
- **注釈器 ID はライブラリのフルパスを含みうる**（`msp_only` の MSP、LBM は常に）。ツールの戻り値に `annotator_id` の生の値を出さない（出すのは導出した真偽かファイル名だけ）。データセットの `feature_metadata[].annotator_id` はファイル名の部分だけを持つ。
