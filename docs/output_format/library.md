# 参照ライブラリ照合（MS/MS スペクトル照合）の出力フィールド

`library_load` / `library_match_feature` / `library_plot_mirror` が返す値の意味と、
`verify_peak_annotation`（`analytical_checks.msms.spectral_match`）に載る同じ形の値。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

ツールが**どのファイルのどの関数をどの順に呼ぶか**は `docs/workflow/library.md`。
ここは**値の意味**だけを定義する。設計の根拠・実測は
[docs/superpowers/specs/2026-09-19-msms-spectral-matching-design.md](../superpowers/specs/2026-09-19-msms-spectral-matching-design.md)
にある。

## 14. 参照ライブラリ照合（`metabolomix/library/`, `metabolomix/analysis/spectral_match.py`）

`.dcl` の測定 MS/MS を参照ライブラリ（`*_Loaded.msp2.dbs` 優先、無ければ `*.msp`）と
突き合わせ、MS-DIAL の個別スコア定義をそのまま移植して数値を出す。目的は「付いている
名前は本当か」を人が対向プロットと数値の両方で判断できるようにすること。個別スコアの
定義一致に加え、**候補の並び順も MS-DIAL の `TotalScore` に合わせてある**（§14.9）。
ただし上流の順位付けタプル全体を再現しているわけではない——同じ §14.9 に限界を書く。

### 14.1 各スコアの意味

`library_match_feature` の候補一覧・`library_plot_mirror` の材料・
`verify_peak_annotation` の `spectral_match.best_match` はすべて同じ 5 指標を使う。

| キー | MS-DIAL 対応 | 意味 |
|---|---|---|
| `simple_dot_product` | `GetSimpleDotProduct` | m/z 重みなし・ペナルティなしの単純ドット積 |
| `weighted_dot_product` | `GetWeightedDotProduct` | m/z 重みとピーク数ペナルティを掛けたドット積。mzTab の `id_confidence_measure[5]`。**単独では並び順に使わない**（§14.9） |
| `reverse_dot_product` | `GetReverseDotProduct` | 参照側のピークだけを基準に測定側で説明できているかを問う非対称スコア |
| `matched_peaks_percentage` | `GetMatchedPeaksScores`（比） | 参照ピークのうち測定側でも観測された割合 |
| `matched_peaks_count` | 同上（count） | 一致した参照ピークの本数 |
| `entropy_similarity` | `GetSpectralEntropySimilarity` | スペクトルエントロピーに基づく類似度。**mzTab には出ないため実データとの突き合わせ対象外**（設計 §2.7） |
| `total_score` | `GetTotalScore` | 上の指標に RT・precursor m/z の一致度を足した**候補の並び順（`ranked_by`）**。**0〜1 ではない**（§14.9） |

**`simple` / `weighted` / `reverse` の 3 つの dot product は平方根側の値**である
（`match_spectrum` が二乗値から `sqrt` を取って返す。単体関数
`weighted_dot_product()` 等はモジュール内部の二乗値を返す——混同しないこと）。
これは mzTab-M の `id_confidence_measure[4..6]` と同じ土俵で比較するための変換で、
`0.0`〜`1.0` の類似度として読める。

### 14.2 `-1` と `0` の区別（最重要）

**`-1` は「比較していない」、`0` は「合わなかった」。混同すると「照合していない」が
「合わなかった」に化ける。**

- `-1`（3 つの dot product・`matched_peaks_percentage`・`matched_peaks_count`・
  `entropy_similarity` いずれも）: 測定・参照どちらかのスペクトルが空で、
  そもそも比較できなかった。`entropy_similarity` も他の 5 指標と同じ
  `_is_compared_available` ガードを使い（`spectral_match.py` の
  `spectral_entropy_similarity`）、`match_spectrum` では `sqrt` を経ずに素通し
  するため `candidates_table` にも `-1.0` がそのまま出る。
- `0`: 比較はできたが、一致度がゼロだった（本当に合わなかった）。**`entropy_similarity`
  はこれにもう1パターン足す**: 測定・参照どちらも非空だが、どちらかの総強度が 0
  （例: `[[100.0, 0.0]]`）という縮退入力も `-1`（比較不能）ではなく `0.0` を返す
  （最終レビュー Important 4）。空ではないので `-1` を流用すると意味が変わるため、
  dot product 3 種が「窓合算後に信号が無かった」ときに `0.0` を返すのと同じ規約に
  揃えた。

**mzTab の値とは番兵の扱いが違う**ので、実データ突き合わせの際に注意する:

- `id_confidence_measure[4..6]`（simple/weighted/reverse dot product）は
  MS-DIAL 側で非二乗 getter が `Math.Sqrt(Math.Max(Squared*, 0f))` として `-1` を
  `0` にクランプしてから `sqrt` を取る。**したがって mzTab 上のこの 3 列では、
  比較不能は `-1` ではなく `0` として出ている。** この関数の `-1` と mzTab の `0` は、
  比較不能ケースに限り同一視してよい。
- `id_confidence_measure[7..8]`（matched peaks count / percentage）は素のフィールドで、
  クランプを経ないため `-1` のまま出る。こちらは `-1` 同士で素直に比較できる。

### 14.3 意図的に写した瑕疵（`metabolomix/analysis/spectral_match.py`）

移植の目的は MS-DIAL との**数値の一致**であって、実装の改善ではない。次の 5 点は
上流の「明らかに変な点」を**意図的にそのまま残している**。素直な実装に書き直すと
比較の土俵が消えるため、直さないこと（docstring にも明記済み）。

1. **`wM` / `wR` の計算が未使用のまま残っている（weighted/reverse）。** 上流にも
   これらの値を使う後続コードが無い。移植していない。
2. **weighted の `if (sumM <= 0 && sumR > 0) {...} else {...}` は両枝が同一処理。**
   条件ごと落として共通処理だけ残している。
3. **simple の `× 999` は比を取る時点で打ち消えるので数値に影響しない。** 忠実に
   写しているが意味は無い。
4. **窓の走査は固定幅グリッドではなく、同じピークが隣り合う 2 つの窓に二重計上され
   得るカーソル走査。** 素直な 1 対 1 アラインメントに書き直してはいけない
   （数値が変わる）。
5. **entropy に Li et al. 2021 の低エントロピー重み変換は入っていない。** 上流にも
   無いので足さない。

### 14.4 `.msp` 由来には許容幅が同梱されない

`.dbs` の `Storage` エントリには、その run が実際に使った許容幅
（`mz_tol` / `ms2_tol` / `rt_tol` 等、`docs/schema/molecule_ms_reference.md` の
`MsRefSearchParameterBase`）が入っている。`library_load` はこれを
`search_params` として store のメタに保存し、`library_match_feature` が
明示指定 > `search_params` > 既定値の順で使う。

**`.msp` にはこの情報が無い。** `.msp` 由来の store は `search_params=None` を持ち、
`library_load` の戻り値に既定値を使う旨の `note` が付く。既定値は
`metabolomix/library/defaults.py` の `DEFAULT_MZ_TOL` (0.01) / `DEFAULT_MS2_TOL` (0.025) /
`DEFAULT_RT_TOL` (0.2)。

**`DEFAULT_RT_TOL` の出所訂正（最終レビュー Important 6）**: これは
`MsRefSearchParameterBase` の `RtTolerance` の既定値ではない
（`docs/schema/molecule_ms_reference.md` によれば上流既定は **100.0**）。
`DEFAULT_RT_TOL` (0.2) は `metabolomix/dcl/reader.py` / `dcl/tools.py` の `.dcl` 検索
（precursor m/z から測定 MS/MS を引くときの RT 窓）の既定値であって、ライブラリ
候補検索用の値ではない。値そのもの（0.2 分）は `.dcl` の窓として妥当。

**ライブラリ照合用の窓と `.dcl` 検索用の窓は別物（最終レビュー Important 7）**:
`library_match_feature` の `mz_tol` / `rt_tol` 引数と、それらが解決した
`search_params` 由来の値は**ライブラリ候補検索（`store.candidates()`）にだけ**
効く。`.dcl` から測定 MS/MS を引く窓（`_measured_spectrum`）は常に
`dcl/reader.py` と同じ固定既定（`mz_tol=0.01` / `rt_tol=0.2`）を使い、
`search_params` の影響を受けない。`.dbs` の `RtTolerance`（実質 RT フィルタを
無効化する 100.0 が典型）をこの窓に流用すると、precursor m/z が近い別ピークの
測定スペクトルを黙って拾ってしまうため。RT が渡されたときは、`.dcl` 側で複数
ヒットしても **RT 距離が最も近いもの**を選ぶ（以前はファイル内の出現順の先頭を
無条件に採っていた）。

### 14.5 `library_load`

`.dbs`/`.msp` を解決して SQLite store を構築（または既存キャッシュを再利用）した
要約を返す。

| キー | 意味 |
|---|---|
| `record_count` | store に格納した参照レコード件数 |
| `ion_modes` | `{"positive": N, "negative": N}` のような件数内訳 |
| `compound_classes` | 化合物クラス別件数の上位 10（`compound_class IS NULL` は除外） |
| `search_params` | `.dbs` 由来なら実測の許容幅と注釈スコアリングのフラグ、`.msp` 由来なら `null`（§14.4） |
| `source_sha256` | 元ファイルの sha256（store のキャッシュキーと同じ） |
| `skipped_no_precursor_mz` | precursor m/z が無い（または解釈できなかった）ため読み飛ばしたレコード件数（最終レビュー Important 3。`record_count` には含まれない） |
| `records_without_ion_mode` | 極性（`IONMODE`）を持たないレコード件数（`record_count` に含まれる）。`library_match_feature` の `ion_mode` 絞り込みはこれらを**極性不明として候補に残す** |
| `non_utf8_lines` | UTF-8 で読めず cp932（無理なら latin-1）で読んだ `.msp` の行数。0 より多ければ化合物名が化けている可能性がある |

**エラー戻り値**: ライブラリを 1 つに決められないときは `{"status": "error", "code": ..., "message": ..., "details": ...}`
を返し、store は差し替えない。`code` は `MSP_AMBIGUOUS`（候補が複数。`ion_mode` か
`file_path` を指定する）/ `MSP_ENV_NOT_FOUND`（極性の設定——環境変数 `MSDIAL_MSP_POS` /
`MSDIAL_MSP_NEG` か `lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`——の
指す先が無い。設定ファイル由来でも同じコード）/ `LIBRARY_NOT_FOUND`（明示した `file_path` が無い）/
`INVALID_ION_MODE` / `CONFIG_INVALID`（設定ファイルが読めない。`details` に `config_file`
`line` `column`）。`details` は `MSP_ENV_NOT_FOUND` と、何も見つからなかったときに付き、
`setting`（設定キー）`env_var` `source`（`env` / `config_file` / null）`config_file`
`config_file_exists` `example`、あれば `unknown_keys`（打ち間違いの疑い）と `config_invalid`（設定ファイルが読めないときの理由）を持つ。
`message` と `details` にはファイル名・設定キー・環境変数名・設定ファイルのパスだけが入り、
ライブラリの置き場所は入らない。

**`records_without_ion_mode`**: 研究室ライブラリのように極性ごとにファイルが分かれて
いると、ファイル側に `IONMODE` 欄が無いことがある。NULL を絞り込みで弾くと候補が
黙って 0 件になるので残す。その代わり、極性違いのファイルを読んでも候補が消えずに
残るので、測定の極性に合ったファイルを読んでいるかは呼び出し側が確かめる
（0 件より多いときは `note` にもそう出る）。

**`skipped_no_precursor_mz`（Important 3）**: reader（`.msp`/`.dbs`）は
PRECURSORMZ 欠損・パース失敗を契約どおり `None` に潰すだけで例外にしないが、
store のスキーマは `precursor_mz REAL NOT NULL`。公開 `.msp`（MassBank 由来など）
には precursor を持たないレコードが普通に混ざるため、そうしたレコードは
1 件ずつ読み飛ばして件数だけ残す（1 件の欠損でライブラリ全体を使用不能に
しない）。0 件より多いときは `library_load` の `note` にも件数が出る。

### 14.6 `library_match_feature`

`status` は 3 通り。

| `status` | 意味 |
|---|---|
| `not_found` | `.dcl` に該当 precursor の MS/MS が無い（**未取得**であって「合わなかった」ではない。`dcl_find_msms` の `not_found` と同じ文言方針） |
| `no_candidates` | 測定 MS/MS はあるが、m/z 窓・極性・RT に該当する参照レコードが無い |
| `success` | 候補を採点して `candidates_table`（TSV、`rank`/`name`/... 列。§14.1 の指標を含む）を返した |

候補は `total_score` の降順（`ranked_by` フィールドが基準を明示）。詳細と限界は §14.9。
スペクトル座標・alignment は戻り値に含めず `session.library.last_match` に持つ
（`library_plot_mirror` がそこから読む）。

`scoring` フィールドは総合スコアの組み立て方を 1 回だけ載せる
（`rule` / `use_rt` / `ms1_tol` / `rt_tol`）。**RT 項が入ったかどうかは候補の数値だけ
見ても分からない**ので、呼び出し側が順位を再現できるようここに出す。
総合スコアの内訳（`rt_similarity` / `mass_similarity` / `spectrum_score`）は
**候補表には出さない**——候補数ぶん掛け算で効いて戻り値が膨らむため、
`session.library.last_match` の各候補の `scores` にだけ持つ。

**`ion_mode` の大小は問わない（最終レビュー Critical 1）**: `store.candidates()`
の比較は `COLLATE NOCASE`。`"positive"`（store の正規化表記）と `"Positive"`
（PAI2 由来の `IonMode.Positive.name`）はどちらも一致する。以前は SQLite の既定
BINARY 照合のままで、`verify_peak_annotation` の統合経路（§14.8）が
`ion_mode.name`（大文字始まり）をそのまま渡すため常に 0 件になっていた
——`no_candidates` はここでは「該当レコードが無い」という実質的な主張なので、
この不一致は「ライブラリに無い化合物だ」という誤った結論に直結していた。

### 14.7 `library_plot_mirror`

`build_mirror_payload` / `render_mirror`（`metabolomix/plots/mirror.py`）の座標契約は
`lipidmix.mirror.v2`。上段が測定（上向き）、下段が参照（下向き）、横軸 m/z 共通。
`matched_mz` は**常に参照側の m/z**（参照グリッドの窓中心）で、測定側の一致 m/z は
別フィールド `matched_measured_mz`（`ms2_tol` を渡したときだけ計算、渡さなければ空）。
`ms2_tol` を渡さない呼び出しでは、無根拠な厳密一致で色を付けないという設計判断により
測定側は一致色分けされない。

**採点に入らなかった測定ピークは灰色で薄く描く（`unscored_mz`、v2 で追加）。**
採点は `normalize_measured` が `.dbs` の `relative_amp_cutoff` / `absolute_amp_cutoff`
で足切りした後のスペクトルに対して行うのに対し、図は**足切り前の生ピーク**を描く。
この差は今まで図の上で見分けられず、「描かれているのに採点されていない」ピークが
黙って混ざっていた。判定の正準は `spectral_match.cutoff_mask`（`normalize_measured`
も同じ mask の上に載る）。落ちたピークも**消さずに描く**——消すと「MS/MS が
取れていない」と読めてしまうため。これらはラベル枠を取らず（一致候補として見た、と
誤読されるため）、凡例の `Below cutoff (not scored)` と画像 caption の件数は
**該当ピークがあるときだけ**出る。`payload` の `scored_peak_count` /
`unscored_peak_count` は常に載る（層が空でも数字で気づけるように）。

**足切りが 0 の run では何も変わらない。** `MsRefSearchParameterBase` の Key 7/8 の
既定は 0.0 で、実データ検証に使った `.dbs`（aging mice kidney neg）も 0.0/0.0 だった
——この層が現れるのは足切りを非ゼロに設定したプロジェクトだけ。

**ピークは素のステムだけで、先端にマーカーは打たない**（上流の
`LineSpectrumControlSlim` も `DrawLine` だけで描く）。マーカーは m/z 軸上の見かけの
太さを増やして近接ピークを潰すだけで、情報を足さないため。

`label_policy` は m/z ラベルの**衝突回避**の方式。どちらも上流 `Annotator.OnRender`
と同じ貪欲法（強度降順に走査し、既に置いたラベルと重なるものを飛ばす）だが、
重なりの見方が違う:

| | `"auto"`（既定） | `"msdial"` |
|---|---|---|
| 判定 | 水平と垂直の両方が近いときだけ飛ばす | **水平距離だけ**見る |
| 箱の大きさ | ラベル自身の文字列を実測 | 代表文字列 `"1000.00000"` 1 つで全ラベル共通 |
| 結果 | ラベルが多く載る | 少なくなる |

**`"msdial"` が水平のみなのは上流の縮退を写したもの。** MS2 ビューの設定は
`Overlap="Horizontal, Direct"` だが、この合成は OR で、しかも全ラベルが同じ
代表箱を使うため `(w1+w2)/2 = repWidth` となり `Direct`（水平 AND 垂直）は
`Horizontal` の部分集合になる——実質 `Horizontal` 単独に縮退している。
縦にどれだけ離れていても m/z が近ければ飛ばす。

どちらの方式も**上下は別枠**で数える（上流は上下で別々の `Annotator` を持つので
側をまたぐ衝突は起きない）。片側 25 本の上限は暴走止めの安全弁で、
`"msdial"` では水平棄却が先に効いて到達しない。

**画像のラベルは `payload["labels"]` ではなくスペクトル全点から選ぶ。**
衝突判定はピクセル座標が要るので描画時にしかできず、`payload["labels"]` は
自前で描くクライアント向けの要約（`top_labels` 件で切った版）という役割分担。
その `top_labels` は**側ごと**の上限なので、`labels` の長さは最大 `2 × top_labels`。

`scale`（`"relative"` 既定 / `"sqrt"` / `"log10"`）は**縦軸の写し方**で、上下それぞれ
自分の最大値で正規化したあとに掛ける。**precursor がベースピークのスペクトル**
（脂質の `[M-H]-` など）は `relative` だと診断イオンが相対数 % に潰れて読めないので、
そのとき `"sqrt"` を使う。大小関係は保つので「どちらが高いか」の読みは変わらない。
`log10` は相対 0.1%（3 桁）を下限に潰す——対数は 0 へ向けて発散するため。
**この下限はこちら独自の値**で、上流の Log10 軸の下限（`ToReactiveLogScaleAxisManager`
の `lowBound=1`）は**生強度 1 カウント**を指す。上流は生強度を軸に載せるのに対し
こちらは先に相対化しているので、同じ土俵の値にはならない。
軸ラベルに選んだスケールが出る。**`output="payload"` の座標は正規化前の生値**
なので `scale` の影響を受けない（自前で描くクライアントが自分で決める）。

上流 MS-DIAL も同じ問題を軸の切り替えで解いている
（`ObservableMsSpectrum.CreateAxisPropertySelectors2` の Relative / Absolute /
Log10 / Sqrt を**上下独立**に選ばせる）。**`Absolute` は用意しない**——対向プロットは
単位の違う 2 つのスペクトルを上下に並べるので、生の強度で並べても比較にならない
（上流は片側ずつ別の図として見られるので成立している）。

### 14.8 `verify_peak_annotation` への統合（`analytical_checks.msms.spectral_match`）

`peak_verification.msms_evidence()` の 3 状態契約（`PASS`/`FLAG_ONLY`/`ABSENT`、
`identity` トピック §12.4）は変えない。`band == "PASS"` かつ参照ライブラリが
読み込み済み（`session.library.store` あり）のときだけ、`spectral_match` ブロックが
追加で載る。ライブラリ未読み込みなら `spectral_match` キー自体を持たない。

**`status` の語彙は暫定である。** `library_match_feature` とは独立に実装されており
（`_spectral_match_for_feature`、`metabolomix/msdial/peak_verification.py`）、
実際の `verify_peak_annotation` の使われ方を見る前に固まらないよう、
現時点では次の 4 値を使う（将来 `library_match_feature` 側の語彙
`not_found`/`no_candidates`/`success` と統一する可能性がある）:

| `status` | 意味 |
|---|---|
| `unavailable` | precursor m/z が feature に無く照合できない |
| `error` | ライブラリ照会自体が例外を投げた（`store.candidates()` の失敗など） |
| `no_candidates` | 候補が 0 件 |
| `matched` | 最良候補（`total_score` 最大。`library_match_feature` と**同じ基準**）のスコアを `best_match` に載せた |

どの段階で失敗しても例外は投げない——`spectral_match` は msms_evidence 本来の目的
（`band` 判定）を止めてはいけないため、失敗は `status` の値として開示するだけに留める。

`best_match` は候補一覧を持たない（候補比較が要るときは `library_match_feature` を
直接使うこと）。`weighted_dot_product` 等は §14.1 と同じ丸め済みの平方根側の値で、
`total_score` も載る。**「最良」の選び方は `library_match_feature` と必ず一致させる**
——食い違わせると、同じ feature について 2 つのツールが別の候補を名指しすることになる。

### 14.9 候補の並び順（`total_score`）と、その限界

`total_score` は上流 `MsScanMatching.GetTotalScore` の移植で、**正規化しない和**:

```
RtSimilarity + AcurateMassSimilarity + (Weighted + Simple + Reverse)/3
+ MatchedPeaksPercentage
```

**各項は `> 0` のときだけ加算する**（`-1` の番兵＝比較不能・欠測を足して総合スコアを
下げないため。§14.2）。したがって**上限は 1 ではなく、項数ぶん（この経路では最大 4）**。
0〜1 の類似度として読んではいけない。RT と precursor m/z の一致度は
`exp(-0.5*((actual-reference)/tolerance)^2)`（`gaussian_similarity`）で、m/z の
許容幅は 500 を超えると ppm 換算で伸びる（`fix_mass_tolerance`）。

**RT 項が入る条件**: store の `search_params["use_time_for_annotation_scoring"]`
（`.dbs` の `MsRefSearchParameterBase` `Key(16)` = `IsUseTimeForAnnotationScoring`）が
真で、かつ測定側の RT と参照レコードの RT が両方ある場合だけ。**上流の既定は
`False`** なので、`.msp` 由来のライブラリでは RT 項は入らない。判断結果は
`library_match_feature` の `scoring.use_rt` に出る。

**`total_score` の絶対値は MS-DIAL の `TotalScore` と比較できない。** 式は同じだが
`matched_peaks_percentage` の**項の意味が違う**——上流の脂質経路
（`GetLipidomicsMatchedPeaksScores`、`MsScanMatching.cs:1169`）は、脂質クラスの
注釈レベルが立つと `resultArray[0] = 2.0` のように**固定のボーナス値で上書きする**
（`:1191`）。つまり上流のこの項は割合ではなく**クラス確証のボーナス点**で 1 を超える。
こちらは汎用の割合（≤ 1.0）を足すので、確度の高い脂質では上流より約 1.0 低く出る。
順位付けの用途では内部で一貫しているので支障はないが、**mzTab の
`id_confidence_measure[1]`（MS-DIAL algorithm matching score）と数値を突き合わせては
いけない**。これが `[8]` の一致率が 68% に留まる理由でもある（仕様差であって不具合ではない）。

**上流の順位付けと完全に同じではない。** `MsScanMatchResultContainer.ResultOrder` は

```
(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)
```

の**辞書順タプル**で、`TotalScore` は最後のタイブレークにすぎない。支配項の
`IsReferenceMatched` は `MsReferenceScorer.ValidateOnLipidomics` の
`IsSpectrumMatch &= isLipidChainsMatch | isLipidClassMatch | ...` を経由して
**脂質クラス固有の判定**（`GetRefinedLipidAnnotationLevel` → `Lipidomics/`、
上流 66,932 行。`MsmsCharacterization.cs` だけで 20,177 行）に依存するため
移植していない。実データで確認したところ、この判定が無いと候補が全件
`IsReferenceMatched=True` になりゲートとして働かない（脂質分岐は 3 つの
dot product の **OR** 判定で、実 run の cutoff が緩いため）。

**帰結**: `total_score` 単独比較は現実的な最良の近似であって、**上位に化学的に
ありえない候補が残ることがある**（例: PI 18:0_20:4 の照合で `SMGDG 18:0_20:4` が
2 位、PC 18:0_22:6 の照合で `SHexCer 41:1;O3` が 2 位）。1 位の妥当性は
対向プロット（`library_plot_mirror`）と `matched_peaks_count` で人が確かめること。

**実データでの効き**（aging mice kidney neg、120 feature、正解は mzTab の
`chemical_name`）: top-1 一致は `weighted_dot_product` 単独の 106/120 (88.3%) に対し
`total_score` は 113/120 (94.2%)。判定が変わった 7 件は全て総合スコア側が正解で、
悪化はゼロだった（HISTRY 2026-09-22(1)）。

### 14.10 MS-DIAL との一致点・相違点（まとめ）

**採点・順位付けは GUI 固有ではない。** 実体は共有ライブラリにある——
`MsScanMatching`（`src/Common/CommonStandard/Algorithm/Scoring/`）、
`MsReferenceScorer` / `MsScanMatchResultContainer` / `MztabFormatExport`
（`src/MSDIAL5/MsdialCore/`）。`MsdialGuiApp` には無い。したがって GUI でも
Console でも同じ判定が走る。一方**描画は GUI だけ**にある
（`src/Common/ChartDrawing/` の `LineSpectrumControlSlim` / `Annotator` と
`MsdialGuiApp/Model/Chart/`）——Console は図を描かない。

#### ロジック

| | こちら | MS-DIAL | 備考 |
|---|---|---|---|
| 個別スコア 5 種 | 同じ | 同じ | 瑕疵ごと移植。実データで `[4][5][6][7]` **194/194 一致**（§14.3） |
| 採点前処理 | 同じ | 同じ | `relative_amp_cutoff` / `absolute_amp_cutoff` / 質量範囲を `.dbs` の実値から |
| 総合スコアの式 | 同じ | 同じ | `GetTotalScore`。正規化しない和、各項 `> 0` のときだけ加算（§14.9） |
| RT / m/z の一致度 | 同じ | 同じ | `gaussian_similarity`、500 超は ppm 換算（`fix_mass_tolerance`） |
| RT 項を足す条件 | 同じ | 同じ | `.dbs` の `Key(16)` を読む。既定 `False` |
| 候補検索の窓 | 同じ | 同じ | precursor m/z ± `Ms1Tolerance`、RT ± `RtTolerance` |
| **並び順** | `total_score` 単独 | `(IsManuallyModified, IsReferenceMatched, IsAnnotationSuggested, Priority, TotalScore)` の辞書順 | ブール項が脂質クラス判定に依存するため未移植（§14.9） |
| **候補の足切り** | しない（採点した全件を返す） | `GetRefinedLipidAnnotationLevel` が空文字を返す候補は `null` で消える。各種 cutoff・`IsSpectrumMatch` でも絞る | **上位に化学的にありえない候補が残る**直接の原因 |
| **名前** | ライブラリのレコード名をそのまま | 注釈レベルに応じて付け直す（クラス止まり / 鎖組成まで） | `GetRefinedLipidAnnotationLevel` の戻り値 |
| **`matched_peaks` の 2 値** | 汎用 `GetMatchedPeaksScores` | 脂質は `GetLipidomicsMatchedPeaksScores`。注釈レベルが立つと `2.0` 等の**ボーナス値で上書き**する（割合ではない） | `[8]` が 68% に留まる理由。`total_score` の絶対値が比較できない理由でもある（§14.9） |
| **足さない項** | RT・m/z・スペクトル・matched% の 4 項のみ | CCS / isotope / Andromeda も足しうる | この経路（LC-MS、MS2、脂質）には存在しない項 |
| **採点対象のファイル** | 呼び出し側が指定した 1 つの `.dcl` | アラインメント spot の「代表ファイル」1 つ | mzTab の値と 0.3〜2% ずれる主因 |

#### 描画（対向プロット）

| | こちら | MS-DIAL | 備考 |
|---|---|---|---|
| レイアウト | 上=測定 / 下=参照、m/z 軸共有 | 同じ | |
| 正規化 | 上下それぞれ自分の最大値 | 同じ（上下で別の軸） | |
| ピークの描き方 | 素のステムのみ | 同じ（`DrawLine` だけ） | 先端マーカーは両者とも描かない |
| m/z ラベルの書式・順序 | 小数 4 桁・強度降順 | 同じ（`Format="F4"`、`OrderingPropertyName`=Intensity） | |
| 縦軸スケール | `relative` / `sqrt` / `log10` | Relative / Absolute / Log10 / Sqrt | |
| **`Absolute` 軸** | 無い | ある（上下独立に選べる） | 対向プロットで単位の違う 2 つを生強度で並べても比較にならないため |
| ラベルの衝突回避 | する（`label_policy="auto"` / `"msdial"`） | する（外接矩形が重なるものを飛ばす） | 貪欲法は同じ。`"msdial"` は判定も箱の測り方も忠実 |
| **既定の判定** | 水平＋垂直の 2 次元（`"auto"`） | 水平のみ（`Overlap="Horizontal, Direct"` の縮退） | `"auto"` のほうがラベルが多く載る。忠実版は `label_policy="msdial"` |
| ラベルの本数上限 | 片側 25 本（安全弁） | 指定なし（MS2 ビューは `TopN` なし） | `"msdial"` では水平棄却が先に効いて上限に届かない |
| **一致ピークの色** | 緑でハイライト＋ガイド線 | 概念が無い（色は `SpectrumPeak.SpectrumComment` 由来。既定は側ごとに 1 色） | こちらの追加。静止画では情報量が多い |
| **採点対象外ピーク** | 灰色で薄く描き、凡例と caption に件数（足切りが非ゼロのときだけ） | 概念が無い（足切り後のスペクトルをそのまま描く） | こちらの追加。§14.7。上流は採点対象だけを見せるので「描かれているのに採点外」が起きない |
| 色 | `#2471a3` / `#c0392b`（volcano と同系統） | 純 Blue / 純 Red | |
| **差分・積の重ね描き** | 無い | ある（`UpperDifferenceSpectrumModel` / `UpperProductSpectrumModel`） | |
| **対話操作** | 無い（静止 PNG、`output="payload"` で座標） | ツールチップ（m/z・強度・`SpectrumComment`）・ズーム | |
