# 再解析の道具（ハンドル・spec・証拠の束）の出力定義

**再解析で使う**ツールの戻り値と、spec・結果の JSON の欄の意味（spec 2026-10-22）。
対話の使い方の出力は変わらない——既存のツールに足した引数（`dataset_handle` ほか）は、
渡さない限り戻り値を 1 文字も変えない。呼び出し連鎖は `docs/workflow/reanalysis.md`。

共通の約束:

- **spec は正規化して返す。** 既定値を埋め、未知のキーは拒否する（黙って無視しない）。
  結果には正規化した spec と `rule_hash` が入る。`rule_hash` は正規化した spec から
  パスの欄（`path` / `matrix_path`）を除いて作る——置き場所が違っても同じ規則なら同じ hash。
  入力のファイル（行列・群・共変量・注入順の TSV）は `input_fingerprint` に sha256 で入る。
- **同じ spec・同じ入力なら同じ JSON。** 時刻や乱数の ID は結果に入れない。
- **大きな結果は `output_path` の JSON に書く。** キーは昇順で、1 行に収まらない辞書・リストは要素ごとに
  改行する（行は辞書のキーと区切りの `,` を含めて 1,000 字以下。数や文字列だけのリストは 1,000 字までずつ
  詰める。字下げはしない。分けられない 1 つの値——1,000 字を超える文字列など——の行だけは超える）——1 行の
  JSON は行の範囲で区切って読めないため。`output_path` を省いて行が上限を超えても**断らない**: 先頭の行
  だけを載せ、`truncated: true` と全件の数を返す（上限と並びはツールごと。下の各節）。
- **効かない設定は誤りにしない。** 書いても使われない設定（recipe の `steps` に無い段の設定、
  `merge_split=false` の `alignment_qc`、ols でない検定の共変量、`rt_min` の無い主張の `rt_offset_min` など）は
  捨てて `notes` に理由を書く。捨てた
  設定は正規化した spec にも `rule_hash` にも入らない。形の誤りと、併用すると意味がぶつかる組み合わせ
  （`aggregate` と `features.feature_ids`）は従来どおり誤り。
- **判定（再現したか・明らかに別物か）は返さない。** 道具が返すのは数字と、spec に書いた
  規則を機械的に当てた結果だけ。
- **研究室の参照ライブラリのパスは出さない。** ファイル名と sha256 だけ。

## 1. データセットのハンドル

ハンドル（`dh_` + 12 桁）は、1 つのサーバのプロセスにデータセット（とその極性の
ライブラリ）を複数並べて持つための ID。`load_dataset` が返し、`dataset_*`・
`library_*`・再解析のツールの `dataset_handle` に渡す。

- **対話のセッションは直近のハンドルを指すだけ。** `load_dataset` が開いたハンドルは**いつも
  直近になる**（`load_dataset` の戻り値の `dataset_handle`、`dataset_list` の `current`）。
  ハンドルを渡さない呼び出しはこの直近のハンドル（のデータセット・`library_load` の
  ライブラリ・直近の照合）を使う。再解析のエージェントは、返ったハンドルを渡して呼び、直近を
  当てにしない——ハンドルを渡した呼び出しは対話の状態を読みも書きもしない。閉じた・無い
  ハンドルは `DATASET_HANDLE_NOT_FOUND`（開いているハンドルの一覧つき）で、対話のセッションへ
  黙って落ちない。直近のハンドルを閉じると、渡さない呼び出しはデータセットもライブラリも無い状態
  （`missing_state`）になり、開いている別のハンドルへも、前に読んだライブラリへも落ちない。
  順番を待つ間に直近のハンドルが閉じられた呼び出しも、何もせずに `missing_state` を返す。
  `load_dataset` より前に `library_load` で読んだライブラリは、このデータセットのライブラリが
  決まらなかった（記憶が無い・`library="auto"`）とき、そのハンドルへ移る（写しではない）。
- **ハンドルどうしも別の状態。** 同じ入力を開き直すと読み込み済みのもの（`reused_base`）を
  使い回すが、ハンドルには複製が渡る。片方の前処理はもう片方に見えない。読み込みが読んだ・
  探したファイル（入力・`.dcl` などの付随ファイル・ライブラリ・ジョブの記録と実効メソッド・
  後から置かれうる `.arf` の候補）の size か更新時刻が変われば読み直す。
- **並行と直列。** ハンドルを渡した呼び出しはサーバの別スレッドで回り、別のハンドルどうしは
  並行に進む。同じハンドルへの呼び出しは 1 本ずつ（到着順は保証しない）。ハンドルを渡さない
  呼び出しも、直近のハンドルがあればそのハンドルの呼び出しとして 1 本ずつに並ぶ（途中で
  `load_dataset` が直近を切り替えても、始めたときの直近のハンドルで最後まで動く）。
- **上限。** `lipidmix.local.toml` の `[server] max_dataset_handles`（既定 8。1 未満は 1）。`load_dataset` がハンドルを 1 本使う
  （前に読んだものは閉じずに残る）。超えると `DATASET_HANDLE_LIMIT` で止め、古いハンドルを
  黙って閉じない。寿命は明示の `dataset_close` か、サーバの終了まで。
- **受付番号。** `load_dataset` は最長 `wait_seconds`（既定 20。0〜30 秒。範囲の外の数は丸めて `wait_seconds_note`
  を付ける）待ち、終わらなければ読み込みを続けたまま `{"status": "loading", "load_id": "ld_…"}` を返す。
  `dataset_load_status(load_id)` は読み込み中なら `status: "loading"`、終われば `load_dataset` の戻り値そのもの
  （1 時間覚える）。`load_id` を省くと、まだ結果を渡していない読み込みが 1 本ならそれ、2 本以上なら `needs_input`
  （受付の一覧つき）、無ければ `missing_state`。`wait_seconds`（既定 0）を渡すと最長 30 秒まで終わりを待ってから
  返す。受付番号を返した後に別のハンドルが直近になっていたら直近にせず、`current: false` が付く。上限は読み込み中の
  受付も 1 本と数える。
- **長い道具の受付番号。** `pathway_project` と `dataset_homolog_check` も `wait_seconds`（既定 20。0〜30 秒）を持ち、
  その中で終わらなければ仕事を続けたまま `{"status": "running", "job_id": "job_…", "tool", "next"}` を返す。ハンドルの
  ロックは仕事が終わるまでその仕事が持つ（同じハンドルへの次の呼び出しは仕事の後に並ぶ）。結果は
  `job_status(job_id, wait_seconds)` で引く（走っている間は `status: "running"` と `elapsed_s`、終われば道具の戻り値
  そのもの。`job_id` の省き方は `dataset_load_status` と同じ。知らない番号は `JOB_NOT_FOUND`（終わって 1 時間を
  過ぎた・件数の上限で古いものから捨てた・サーバが再起動した。道具を呼び直す）。長い道具を呼んだとき走っている受付が
  上限なら、その道具が `JOB_LIMIT`（`details.running` と `max_entries`。`job_status` でどれかが終わるのを待って呼び直す）。
  道具が `DomainError` 以外で止まれば、その結果（道具の戻り値か `job_status`）が `JOB_FAILED`（`details.exception` は
  例外の種類だけ。本文はサーバのログ）。
- **ハンドルを見失う条件。** ハンドルは次のときに無くなり、以後そのハンドルを渡した呼び出しは
  `DATASET_HANDLE_NOT_FOUND`（`details.open_handles` に開いているハンドル）を返す（直近や別のハンドルへは落ちない）:
  (1) `dataset_close` で閉じた（戻り値の `source` が閉じたものの入口の種類・入力の名前・極性）、(2) サーバのプロセスが
  終わった・再起動した（`server_update` の後の再起動、クライアントの再接続でサーバが起動し直したときも。ハンドルは
  プロセスの中にだけある）、(3) 読み込みが受付の後に失敗した（ハンドルは登録されない）。上限に達しても古いハンドルは
  閉じない。見失ったら `dataset_list` で開いているものを確かめ、`load_dataset` で開き直す（同じ入力なら読み込み済みの
  基を使い回すが、前処理・比較の結果はハンドルごとなので作り直す）。

`load_dataset` の戻り値:

| 欄 | 意味 |
|---|---|
| `dataset_handle` | ハンドルの ID |
| `assay_kind` | このハンドルの種別（`lipid` / `metabolite` / `unknown`）。ハンドルごとに持ち、解析フォルダ配下 `reports/` の目的ファイルに記録があれば読み戻す（`notes` に注記）。データセットを開く前に（`[paths] reports_dir` の目的へ）記録した種別は、次に直近になるハンドルへ引き継ぐ（目的ファイルの種別と食い違えば `unknown`。`notes` に注記）。`record_objective(assay_kind=...)` で確定・変更する |
| `adapter` | 読んだアダプタ。`msdial_folder`（解析フォルダの `.arf` / `.arf2` を直接）/ `mztab`（mzTab-M）。`job_path` / `pipeline_path` の既定は `msdial_folder`（`source="mztab"` で `mztab`） |
| `polarity` / `polarity_source` | データセットの極性と出所（`job_declared` / `adduct_majority` / `arf2_ion_mode` / `unknown`）。`dataset_feature_search` の目標照合と同じ決め方 |
| `n_features` / `n_samples` | 特徴数・試料数 |
| `detection_available` | 検出状態（gap-fill の区別）を取り込めたか |
| `arf` / `not_chosen` | 解析フォルダで読んだとき、選んだ `.arf` と選ばなかった `.arf` のファイル名。mzTab-M で読んだときは `null` |
| `library` | `{file, sha256, record_count}`。読まなかったときは `null` で、理由が `library_note` |
| `library_record` | データセットが同定に使ったライブラリの記録（[library.md](library.md) の LibraryRecord。置き場所は出さない） |
| `arf2` / `dcl` | スポットの照合（spot_id = SMF_ID）に使うファイル名。解析フォルダで読んだならその `.arf2` / `.dcl`。決まらなければ `null` |
| `reused_base` | 読み込み済みのものを使い回したか |
| `summary` | 特徴・試料の数、定量種別、検出状態、同定の内訳、出所の裏取り、読み込みの注記（行の配列） |
| `sample_counts` | `{groups, roles, n_ungrouped}`: 群（Class ID）ごと・役割ごとの試料数と、群の宣言が無い試料の数 |
| `overview` | `overview=True` のときだけ。最初の概観（[dataset.md](dataset.md) §11.9。既定は載せない） |
| `wait_seconds_note` | `wait_seconds` を 0〜30 秒に丸めたときだけ |
| `next` | 次の一手の例（3 行。最初が `dataset_overview()`） |

ライブラリはいつもファイル名と sha256 だけを返す（置き場所は出さない）。

ライブラリは**データセットが同定に使ったものの記憶**（`library_record`。[library.md](library.md)）から
引く。実体の sha256 を測り直して `known` のままなら読み取り専用で開き、`mismatch`（内容が変わった・
ファイルが無い）・`unknown`（決められない）・`lbm_only`（LBM だけ）は開かずに `library_note` に
理由を書く。**極性の設定（`[library] msp_positive` / `msp_negative`）には落ちない**ので、逆の
極性のライブラリを掴まない。
ハンドルへの `library_load`（`file_path` を省く）も、記憶が `known` ならそれを開く。`mismatch` /
`lbm_only` は理由（`LIBRARY_RECORD_MISMATCH` / `LIBRARY_RECORD_LBM_ONLY`）を返し、記憶が `unknown` なら
`LIBRARY_NOT_REMEMBERED`（極性の設定には落ちない）。`file_path` を明示すれば
記憶と違っても開き、戻り値の `differs_from_record` に記憶したファイル名と開いたファイル名を書く。

**ハンドルは対話のデータフォルダを探さない。** スポットの照合（`library_match_feature` の
spot_id / precursor_mz・`dcl_spot_spectrum`・`dcl_diagnostic_ions`・`library_msms_qualify`・
証拠の束の MS/MS）は、ハンドルの `.arf2` / `.dcl`（か明示した `arf2_path` / `dcl_file`）だけを
使い、無ければ `status: "no_alignment_files"` を返す——別のデータセットの同じ spot_id の
スペクトルを黙って読まないため。

## 2. 統計の spec（analysis-spec.v1 / v2）と結果（analysis-result.v1 / v2）

**この節の前半（表と欄）は v1。** v2（統計エンジン）は §2.1。`dataset_analysis_run` は `spec.schema` で振り分ける。

`dataset_analysis_run(spec, dataset_handle=None, output_path=None)`。論文の方法を spec に写し、
名前つきの解析（`analyses`）を並べる——感度解析は同じ spec に並べた別の解析。

```json
{"input": {"source": "dataset"},
 "groups": {"source": "metadata"},
 "comparison": {"reference": "HC", "test": "RA"},
 "injection_order": null,
 "covariates": null,
 "features": {"all": true, "merge_split": true},
 "aggregate": null,
 "analyses": [{"analysis_id": "main",
               "matrix": {"steps": ["detection_filter", "impute", "normalize", "drift_correct", "qc_rsd_filter"], "...": "§3"},
               "test": {"kind": "student", "transform": "log2", "covariates": []},
               "fold_change": "arithmetic_mean_ratio", "correction": "bh",
               "bh_population": "matrix", "thresholds": {"p": 0.05, "q": null, "fc": 1.5}}]}
```

| 欄 | 既定 | 意味 |
|---|---|---|
| `input` | `{"source": "dataset"}` | データセット（ハンドルか対話のセッション）か、外部の行列。TSV は `{"source": "tsv", "matrix_path", "orientation": "features_by_samples" \| "samples_by_features", "id_column"}`（較正用。論文の個体ごとの値など。欠測は空・`NA`・`NaN`・`null`、数として読めない値は誤り）。論文側の表は MetaboLights の MAF `{"source": "maf", "path"}`・Metabolomics Workbench の mwTab `{"source": "mwtab", "path"}`・補足の xlsx `{"source": "xlsx", "path", "sheet": null, "feature_column": "A", "first_sample_column": "B"}`（下の「外部の表」） |
| `groups` | データセットなら `metadata` | `metadata` は実験情報シート（`group`・`role`・`include`）。TSV は `{"source": "tsv", "path", "sample_column": "sample", "group_column": "group", "role_column": null}`。外部の行列には TSV が要る |
| `comparison` | 必須 | `reference`（基準）と `test`（比較）の群の名前。比較名は `<test>_vs_<reference>` |
| `injection_order` | `null` | ドリフト補正の注入順。`{"source": "metadata"}` か TSV `{"source": "tsv", "path", "sample_column", "order_column": "order"}` |
| `covariates` | `null` | OLS の共変量の TSV `{"path", "sample_column": "sample", "missing": "drop" \| "error", "columns": [{"name", "type": "numeric"} \| {"name", "type": "categorical", "reference"}]}`。カテゴリは基準の水準に対する 0/1 の列に符号化する |
| `features` | `{"all": true, "merge_split": true}` | 報告する特徴（`{"feature_ids": [...]}` で絞る）と、割れたスポットのまとめ（下の「割れたスポットのまとめ」）。BH の母集団は `bh_population` が決める。`merge_split` は無ければ true（v8 の既定）、`alignment_qc` は `dataset_alignment_qc` の `output_path` の JSON のパス（任意。`dataset_alignment_qc` がハンドルに持つ最後の結果は**使わない**——spec だけで再現するため） |
| `aggregate` | `null` | 1 本に集約して検定する。`{"kind": "sum", "label", "feature_ids": [...]}`・`{"kind": "sum", "label", "select": {"name_regex", "evidence": ["sme"]}}`（クラスの主張）・`{"kind": "ratio", "label", "numerator", "denominator"}`（比の主張）。和はメンバーごとの寄与（`contributions`）も返す。**`features.feature_ids` とは併用できない**（報告するのは集約した 1 本なので、報告する特徴の定義がぶつかる。`ANALYSIS_SPEC_INVALID`・`field: features`。種ごとの行は `aggregate` の無い別の spec で回す） |
| `correlate` | `null` | 外部の表の値との Spearman の ρ（下の「外部の表との相関」）。報告の数字（行・q）は変えない |

解析（`analyses[]`）:

| 欄 | 既定 | 意味 |
|---|---|---|
| `matrix` | `{"steps": []}`（前処理なし） | 行列の recipe（§3） |
| `test.kind` | 必須 | `student`（等分散の t）・`welch`・`mannwhitney`（両側）・`ols`（`transform(y) ~ 1 + 群 + 共変量`、群の係数の t の両側 p） |
| `test.transform` | `log2` | `log2` は有限の正値だけに当てる（0 以下は欠測。pseudocount も clip も無い）。`none` |
| `test.covariates` | `[]` | `ols` で使う共変量の名前（`covariates.columns` の部分集合） |
| `fold_change` | `arithmetic_mean_ratio` | `arithmetic_mean_ratio`（変換前の算術平均の比 test/ref）・`geometric_mean_ratio`（log2 の平均の差） |
| `correction` | `bh` | `bh`・`none` |
| `bh_population` | `matrix` | `matrix`（その解析の行列で検定できた全特徴）・`reported`（報告する特徴だけ）・`annotated`（SME のある特徴 = `ds.feature_candidates` に候補のある特徴だけ。論文が同定した種だけで補正したと読めるときの二次解析用）。検定できない特徴は母集団に入れない。`annotated` の母集団の外の特徴は、報告しても `q` は null。まとめた特徴は元のどちらかに SME があれば母集団に入る。外部の行列（SME を持たない）には使えない |
| `thresholds` | すべて `null` | `p`・`q` は「未満」、`fc` は `|log2fc| ≥ log2(fc)`。null の閾値は見ない |

**割れたスポットのまとめ（`features.merge_split`、既定 true）。** MS-DIAL のアライメントが 1 つの
物質を 2 スポットに割ると、片方の「無い」試料に gap-fill の背景が入り、群の差でもない差が出る。
`merge_split` が真なら、**recipe の前**（定量行列と検出状態を読んだ直後）に、割れた 2 スポットを
試料ごとの大きい方の値（NaN は無視）の 1 列にまとめ、列の ID は 2 スポットのうち受け持ちの多い側（同数なら小さい ID。§9）に
する。まとめた列の検出状態は「どちらかで検出」。recipe はそのまとめた行列に当たる——割れの少数側は
検出された試料が少なく、まとめる前に検出率の篩を通すと落ちてまとめが打ち消されるため、また
正規化（TIC・PQN・KNN）が 1 つの物質を 1 列として見るため。`features.alignment_qc` の JSON は
`dataset_fingerprint` が今のデータセットと合わなければ `ANALYSIS_INPUT_INVALID`（§9）。
まとめるのは `alignment_qc`（§9）の `merge_groups`（互いに最も近い mergeable な組。常に 2 スポット）
だけで、`features.alignment_qc` に `dataset_alignment_qc` の JSON を渡せばその群、無ければその場で
`find_split_pairs`（既定の値）を回す。`group_aligned` などの mergeable でない組はまとめない。
まとめるのは行列全体なので、`bh_population=matrix` の母集団も、まとめた後の特徴の数になる。
報告する特徴（`features.feature_ids`）に消えた側の ID を書けば、残した 1 行（`merged_from` 付き）に当たる
（両方書いても 1 行）。`aggregate` の `feature_ids`・`select`・比の分子と分母も同じ置き換えをしてから和・
比をとるので、割れた 2 スポットが両方選ばれても 1 回だけ数える（`aggregate.members` は選んだままの
並び、`summed` はまとめた後の並び、`aggregate.merged_from` は選んだ特徴のうちまとめたもの。比の分子と
分母が同じ物質にまとまるなら `ANALYSIS_SPEC_INVALID`）。行列に 2 つ揃わない群（別の run の JSON など）は
まとめず、結果の `merge_split.skipped` に挙げる。外部の行列（`input.source` が `dataset` 以外）にはアライメントが無く、
`merge_split` は何もしない（`features.alignment_qc` は効かないので捨てて `notes` に書く。`bh_population=annotated` は
`ANALYSIS_SPEC_INVALID`）。`merge_split=false` にすると 2 スポットのまま検定する（`alignment_qc` は効かないので捨てて
`notes` に書く）。`dataset_alignment_qc` がハンドルに持つ最後の結果（§9）は使わない——同じ spec なら誰が回しても
同じ組でまとめるため。既定でない閾値の組を使うなら、その JSON を `features.alignment_qc` に書く。

**外部の表（MAF・mwTab・xlsx）。** 論文側の値の表を、特徴 × 試料の行列として読む（標準ライブラリだけ。
reanalysis-study の読み取りを移したもの）。数として読めないセル（`ND`・`<LOQ` など）と非有限の値は欠測にし、
**空でない読めないセルの数と例**（`input_read.n_unreadable_cells`・`unreadable_examples`。例は最大 5 種）を返す
（黙って捨てない。空・`NA`・`NaN`・`null`・`none` は欠測で数えない）。試料・特徴の重複は `ANALYSIS_INPUT_INVALID`、
ファイルが無ければ `ANALYSIS_INPUT_NOT_FOUND`。群は TSV で渡す（`groups.source=metadata` は使えない）。

- **MAF**（`m_*.tsv`）: 見出しの `smallmolecule_abundance_std_error_sub` より右の列が試料。特徴 ID は
  `row_<n>:<metabolite_identification>`（`n` は見出しの次の行を 1 とするデータの行の番号。空の行も数えるので、
  同じ名前の行も別の特徴）。見出しにこの列が無ければ `ANALYSIS_INPUT_INVALID`。
- **mwTab**（テキスト形式）: `MS_METABOLITE_DATA_START` 〜 `MS_METABOLITE_DATA_END` の区画。1 行目が `Samples`、
  `Factors` の行は読まない。特徴 ID は 1 列目の名前。
- **xlsx**: 1 枚のシート（`sheet` を省けば先頭。無い名前はあるシートの名前を添えて `ANALYSIS_INPUT_INVALID`）の
  横長の表。列は Excel の文字（`A`・`B`・`AA`。`first_sample_column` は `feature_column` より右）。見出しの行は
  2 つの列の両方に値がある最初の行で、それより上の表題・注記は読まない。試料は見出しの行の
  `first_sample_column` から右の空でないセル。試料のセルに数が 1 つも無く、欠測の印でない文字がある行
  （群の見出しなど）は特徴にせず、その数と名前の例（最大 5）を `input_read.n_rows_skipped_as_text`・
  `rows_skipped_as_text` に返す（xlsx だけの欄）。欠測の印（`ND`・`n.d.`・`NQ`・`BDL`・`LOD`・`LOQ`・`tr`・
  `-`・`n/a`・`not detected`、`<`・`>`・`≤`・`≥` で始まるもの（`<LOQ`・`<0.1`）、`#` で始まる Excel のエラー値）
  だけの行は、どの試料でも検出されなかった物質として特徴に読み、そのセルを読めないセルに数える。spec の
  特徴 ID（`features.feature_ids`・`aggregate`・`correlate.pairs`）が飛ばした行の名前なら、`FEATURE_NOT_FOUND`
  がそう書き添える（`details.skipped_as_text`）。式は計算せず、保存された値を読む。

`rule_hash` にパス（`path`）は入らない。ファイルの sha256 は `input_fingerprint.files.matrix`。

**外部の表との相関（`correlate`、任意）。** 論文の表の値とこちらの値の一致を、特徴の対応ごとの Spearman の ρ で
見る（較正用。判定はしない）。

```json
"correlate": {"external": {"source": "mwtab", "path": "ST000000.txt"},
              "pairs": [{"feature_id": "1265", "external_id": "Glucose"}],
              "sample_map": {"path": "map.tsv", "sample_column": "sample", "external_column": "external_sample"},
              "stage": "input"}
```

- `external` は外部の表（`input` と同じ形: TSV・MAF・mwTab・xlsx）。`pairs` はこちらの `feature_id` と外部の表の
  特徴 ID（`external_id`）の組で、**対応は spec が渡す**（名前の照合はしない）。どちらかの側に無い ID は
  `FEATURE_NOT_FOUND`（割れの消える側の ID は残す ID に当たる。行の `merged_into`）。
- `sample_map`（任意）は試料名の対応の TSV（こちらの試料 → 外部の表の試料）。無ければ同じ名前で照合する。
  当たった試料だけを使い、1 つも当たらなければ `SAMPLE_MISMATCH`。対応は 1 対 1 で、こちらの試料の重複も、
  外部の同じ試料を複数の試料が指す対応（外部の 1 つの値を 2 回数えて n と p が膨らむ）も `SAMPLE_MISMATCH`
  （`details.duplicated` / `duplicated_external`。反復の注入は対応表に 1 つだけ書く）。
- `stage`: `input`（既定。読んだ行列＝割れのまとめの後・recipe の前。群・役割を問わず当たった全試料）か
  `matrix`（解析ごとの recipe の後の行列・集約の前。行に `analysis_id`。recipe が落とした特徴は
  `reason: not_in_matrix`）。
- ρ は両方が有限の試料だけで、同順位は平均の順位、p は t 近似（`scipy.stats.spearmanr` と同じ）。組が 5 未満なら
  ρ と p は null（`reason: too_few_pairs`）、どちらかが一定なら null（`constant`）。
- 結果の `correlate` は `{stage, external（読み取りの要約）, n_samples_matched, unmatched_samples,
  unmatched_external_samples, rows: [{analysis_id, feature_id, external_id, merged_into, n, rho, p, reason}]}`。
  戻り値の TSV の行は ρ を小数 4 桁、p を有効 4 桁にする（小さい p も 0 にしない。JSON は丸めない）。
  外部の表と対応表の sha256 は `input_fingerprint.files` の `correlate.external`・`correlate.sample_map`。
  書いたときだけ正規化した spec と `rule_hash` に入る（書かない spec の hash は変わらない）。

**試料の照合は完全一致。** 群の TSV は行列の試料と過不足なく一致しなければ
`SAMPLE_MISMATCH`（欠け・余り・重複を名前で挙げる）。注入順・共変量の TSV は行列の試料を
すべて含めばよい（余りは無視）。群に入るのは役割が `sample` で `include` が真の試料だけで、
群の n が 2 未満なら `STATISTIC_GROUP_TOO_SMALL`。データセットの試料名は実験情報シートの
`sample_id`（無ければ assay 名）。

結果の行（解析ごと。TSV・JSON とも同じ欄）:

| 欄 | 意味 |
|---|---|
| `status` / `reason` | `tested`・`not_testable`（`insufficient_finite_values`・`zero_variance`・`rank_deficient`・`no_residual_degrees_of_freedom`・`test_statistic_not_finite`）・`not_in_matrix`（recipe の篩で落ちた） |
| `group_ref` / `group_test` | 群の名前（比較の向き） |
| `n_ref` / `n_test` | 検定に使った値の数（変換後に有限なもの。OLS は共変量の欠測で除いた後） |
| `mean_ref` / `mean_test` | 行列の値（変換前）の算術平均 |
| `fc_test_over_ref` / `log2fc_test_over_ref` | `fold_change` の定義の倍率と、その log2。**正なら test が高い** |
| `statistic` / `p` / `q` | t（student・welch・ols）か U（mannwhitney）、両側 p、BH の q |
| `beta_log2` | OLS の群の係数（変換後の尺度。log2 なら log2 の差） |
| `direction` | `test_higher` / `ref_higher`（OLS は `beta_log2`、ほかは FC の向き） |
| `meets_thresholds` | 与えた閾値をすべて満たすか（閾値が無ければ null）。**再現の判定ではない** |
| `merged_from` | 割れた 2 スポットをまとめた特徴だけ、まとめた特徴 ID（先頭が `feature_id` = 残した ID。受け持ちの多い側、同数なら小さい ID）。ほかは null（TSV では `;` 区切り） |
| `detected_reference` / `detected_test` | 検定に入れた群の試料のうち、検出された（gap-fill でない）数。gap-fill は数えない。まとめた特徴は「どちらかのスポットで検出」の試料を数える。外部の行列・検出状態の無いデータセット・`aggregate` の行・`not_in_matrix` の行は null。群の試料数は `n_ref` / `n_test` ではなく群の `members`（値が欠けた試料も入る） |
| `detection_fisher_p` | 検出の有無（検出 / 未検出）の 2×2 に対する Fisher の正確検定の両側 p（上の 2 欄が null なら null） |
| `constant_in_group` | 一方の群の値（検定に入れる行列の値。変換の前）が、有限の値 2 つ以上ですべて同じなら、その群の名前（TSV では `;` 区切り。両群なら 2 つ）。ほかは null。未検出を同じ値で埋めた群（補完の `half_min` など）の印で、検定の値は変えない（**判定には使わない**。差が埋め値から来ていないかを見る）。行の最後の欄 |
| `detection_driven` | 差が値の大小より検出の有無から来ているか。高い側の群（`direction`）の検出が群の試料数の半分に満たない（`<`）か、低い側の群の検出が 0 なら true。`status` が `tested` で向きがある行だけ（ほかは null）。`direction` が null の行は `detection_driven` も null（一方の群の平均がちょうど 0 で、変換も補完も無い行は向きが決まらないので、検出の差が大きくても null）。**判定（`meets_thresholds` など）は変えない。使い方は手順書** |

解析ごとに `matrix`（`n_features_in`（まとめた後、recipe の前）・`n_features_out`（recipe の後）・段ごとの報告
`steps`）、`bh_population_n`、`detection_available`（検出の例数を出せたか）・`n_detection_driven`・`n_constant_in_group`（要約のみ）、集約なら `aggregate`（和: `members`・`summed`・`dropped`、まとめがあれば
`merged_from`、寄与 `contributions` と `n_samples_for_shares`。比: `numerator`・`denominator`・`dropped`）、OLS なら `covariate_excluded_samples`。

**和の寄与（`aggregate.contributions`）。** 和の各メンバー（`summed` の順）が和に占める割合（試料ごとの
メンバー ÷ 和。値は行列の値＝変換の前）の要約 `{feature_id, share_median, share_median_ref, share_median_test,
share_max}`。使うのは比べる 2 群の試料のうち和が有限で正のもの（その数が `n_samples_for_shares`）で、
`share_median` は 2 群を合わせた中央値、`_ref` / `_test` は群ごと、`share_max` は最大。少数の種が和を支配して
いないか・群で寄与が入れ替わっていないかを見る材料で、判定には使わない（試料ごとの割合の和は 1 だが、
中央値の和は 1 とは限らない）。
結果の全体に `spec`（正規化後）・`rule_hash`・`merge_split`（`enabled`・`source`: `find_split_pairs` か
`alignment_qc`・`groups`: まとめた群（残す ID が先頭）・`groups_hash`・`n_groups_given`・`skipped`）・
`input_fingerprint`（`dataset_hash`・入力の TSV と `alignment_qc` の JSON ごとの sha256）・`comparison`、
外部の表（MAF・mwTab・xlsx）の入力なら `input_read`（`{source, n_features, n_samples, n_unreadable_cells,
unreadable_examples}`。xlsx は `n_rows_skipped_as_text`・`rows_skipped_as_text` も）、`correlate` を書いたなら `correlate`、効かない設定を捨てたなら `notes`。
`rule_hash` は正規化した spec（パスを除く）に、データセットの入力で `merge_split` が真なら使った群の
hash も入れる——同じ spec でも割れの組が変われば別の規則で、`merge_split` を変えても変わる
（外部の行列では `merge_split` は hash に入らない）。

**戻り値（v1）。** `{status, schema, rule_hash, merge_split: {enabled, source, n_groups_used, n_groups_skipped,
groups_hash}, input_fingerprint, input_read?, comparison, group_ref, group_test, output_path?, analyses, correlate?,
notes?, note}`。`analyses[]` は解析ごとの要約（`n_rows`・`n_tested`・`n_meets_thresholds`・`n_constant_in_group` ほか。
和の `aggregate.contributions` は TSV）。`output_path` を省けば解析ごとに行を TSV（`rows`。上の行の欄の順）で返し、
**50 行を超える解析は p の小さい順の先頭 50 行だけ**を載せる（`truncated: true`・`rows_order: "p_ascending"`。
`n_rows` は全件。p が null の行は後ろ。断らない）。`correlate` の行も TSV で、50 を超えれば先頭だけ（`truncated`）。
全量は `output_path` の JSON。

### 2.1 統計エンジンの spec（analysis-spec.v2）と結果（analysis-result.v2）

`spec["schema"] == "analysis-spec.v2"` の `dataset_analysis_run` は統計エンジン（`metabolomix/analysis/stats/`）へ回る
（v1 の上の表・行の欄とは別の系。`analysis-spec.v1` の spec は今までどおり v1 が実行する）。前処理を **1 回**
回し、`analyses` の順に compare / composition / pca をすべて実行して**全部**を返す（`primary` は報告の主役を
指すだけ）。どれかの解析が失敗したら全体を失敗にする（一部だけ載せた結果は返さない）。**状態も記録も書かない**
（対話の `dataset_compare` と違い、`result_id` も `matrix_id` の登録も作らない）。

```json
{"schema": "analysis-spec.v2",
 "input": {"source": "dataset", "exclude_samples": [], "exclude_features": [], "group_overrides": {}},
 "preprocess": {"drift": {"mode": "auto", "span": 0.75, "min_qc_per_batch": 4},
                "filters": {"blank_ratio": {"min": 3.0}, "detection": {"by": "group", "min_rate": 0.5},
                            "qc_rsd": {"max": 0.30}, "d_ratio": {"max": 0.5}},
                "normalize": {"method": "auto"}, "collapse_replicates": true, "injection_columns": []},
 "design": {"factors": ["group"], "covariates": [], "block": null},
 "comparisons": [{"id": "RA_vs_HC", "factor": "group", "reference": "HC", "test": "RA"}],
 "analyses": [{"analysis_id": "main", "kind": "compare", "test": "moderated_t",
               "effect_size": "log2_mean_difference", "bh_population": "matrix", "p_adjust": "bh",
               "verdict": {"q": 0.05, "equivalence_fc": 1.2}, "min_detected_fraction": 0.5,
               "allow_confounded": true}],
 "primary": "main"}
```

| 欄 | 既定 | 意味 |
|---|---|---|
| `input` | `{"source": "dataset"}` | `source` はデータセットだけ。`exclude_samples` / `exclude_features`（空でない文字列のリスト。重複は順序を保って落とす）と `group_overrides`（試料 → 群）は**この spec に書いたものだけ**が効く（下の「自己完結」） |
| `preprocess` | 上の値 | 前処理の設定（`dataset_preprocess` の引数と同じ。`drift.mode` は `auto` / `on` / `off`、足切りは `null` でその足切りを行わない、`normalize.method` は `auto` / `none` / `tic` / `median` / `pqn`。任意の欄 `drift.loess_degree`（`1` 既定 / `2`）と `normalize.features` / `normalize.exclude_features`（正規化の係数を求める特徴 ID。既定の値のときは spec に現れない）は `dataset_preprocess` の `loess_degree` / `normalization_features` / `normalization_exclude_features` と同じ）。`injection_columns` は個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う属性）。不正なら `ANALYSIS_SPEC_INVALID`（`details.field` が `preprocess.` から始まる） |
| `design` | `{"factors": ["group"], "covariates": [], "block": null}` | `factors`（重複のない文字列。`group`・`batch` か、試料の表の列名。sample-manifest v2 の後ろに足した因子の列も使える。交互作用の対比は 2 因子のモデルだけ）・`covariates`（`{name, type: numeric \| categorical, reference?}`）・`block`（`"subject_id"` など） |
| `comparisons` | 必須 | 対比のリスト。`{id, factor, reference, test, at?}`（`at` は他の因子の水準で層別）か、交互作用 `{id, kind: "interaction", factors: [2 因子], levels: {因子: [基準, 比較]}}`。`id` は英数字と `_ . -` で重複不可 |
| `analyses[].kind` | 必須 | `compare` / `composition` / `pca` / `plsda`。`analysis_id` は重複不可 |
| compare / composition | | `test`（既定 `moderated_t`。`welch` / `student` / `mannwhitney` は 1 因子・共変量なし・ブロックなしだけ）・`effect_size`（`log2_mean_difference` 既定 / `arithmetic_mean_ratio`）・`bh_population`（`matrix` 既定 / `annotated` / `reported`。`reported` は `reported_feature_ids` が要り、そのときだけ書ける）・`p_adjust`（多重比較の補正。`bh` 既定 / `bonferroni` / `holm` / `none`。`q` はこの方法の値）・`verdict.q`（既定 0.05）・`verdict.equivalence_fc`（既定 1.2）・`min_detected_fraction`（既定 0.5）・`allow_confounded`（既定 true。完全交絡でも比べて `confounding` と行の `confounded` に印を付ける。false なら完全交絡で `STATISTIC_CONFOUNDED` で止める） |
| pca | | `scaling`（`autoscale` 既定 / `pareto`）・`n_components`（既定 5）・`impute`（`half_min` 既定 / `knn` / `none`） |
| plsda | | PLS 判別分析。クラスは `design.factors` の水準の組（cell。`classes` に cell 名を並べて絞れる。2 つ以上）。`n_components`（既定 2）・`scaling`（`autoscale` 既定 / `pareto`）・`impute`（`half_min` 既定 / `knn` / `none`）・`cv_folds`（層化 k 分割の数。既定 7、2 以上）・`n_permutations`（並べ替え検定の回数。既定 200、1〜10000）・`seed`（分割と並べ替えの乱数。既定 0） |
| `primary` | compare / composition が 1 つならそれ | 2 つ以上あるときは必須（`ANALYSIS_SPEC_INVALID`）。pca / plsda は選べない |

**自己完結（再現のための規則）。** v2 は spec だけで同じ数字を出し直せる実行で、`dataset_exclude` の除外は読まない。
spec の `input` に無い除外がデータセットにあれば、結果の `warnings` に
`{"code": "dataset_exclusions_not_in_spec", "n_samples", "n_features"}` が載る（効いていない除外の件数）。対話の
`dataset_compare` が結果の spec に除外を写して持つのとは逆の向き。組成のクラスの名前は、そのハンドルの
`assay_kind`（`record_objective(assay_kind=...)`。spec には載らず、戻り値の `assay_kind` に出る）で決まる（`class_label_kind`: `lipid_class` /
`msdial_ontology`）。同じ spec でも `assay_kind` が違えばクラスの切り方が違うので、組成を再現するときは種別も揃える。

**前処理は尋ねない。** `drift.mode="auto"` は、補正できても内部標準で換算済みの値か、注入順が MS-DIAL の既定
（ファイルの順）かもしれないときは補正を当てずに回す。前処理の `steps` のドリフトの段に `reason`・`reasons`・
`rerun_with: {"preprocess.drift.mode": "on"}` が、`warnings` に `drift_not_applied_by_default` が載る。補正するなら
spec の `preprocess.drift.mode` を `"on"` にして呼び直す（`rule_hash` が変わる）。

**戻り値は要約だけ。** `{status, schema: "analysis-result.v2", rule_hash, matrix_id, input_fingerprint: {dataset_hash,
sample_table_hash}, primary, warnings, analyses}`。`analyses[]` は解析ごとの要約で、compare は `n_rows`・`n_tested`・
`verdicts: {changed, unchanged, inconclusive}`・`n_q_below`（`q ≤ 閾値`）・`test`・`warnings`、composition はそれに
`tables`（`clr` / `absolute` ごとの同じ件数）・`class_label_kind`・`n_classes`、compare はさらに `p_adjust_method` と、
要因の F があれば `factorial`（効果ごとの `{df1, n_tested, n_q_below}`）、pca は `scaling`・`imputation`・`n_fit`・
`n_units`・`n_qc_projected`・`n_features`・`explained_variance_ratio`、plsda は `classes`（cell → 個体数）・`n_units`・
`n_features`・`scaling`・`imputation`・`n_components`・`r2x_cum`・`r2y_cum`・`q2_cum`・`q2`・`cv`・`permutation`・`vip_top`
（VIP の上位 10）・`warnings`（pca と plsda は行を持たない）。長い一覧は切る。全量（行・CI・事前分布・交絡の
表・PCA のスコア）は `output_path` の JSON に書く。`output_path` を省けば行を TSV で返し（列は `contrast_id feature_id
status reason log2fc ci_low ci_high p q verdict`、composition は先頭に `table` 列）、**50 行を超える解析は p の小さい順の
先頭 50 行だけ**を載せる（v1 と同じ。`truncated: true`・`rows_order: "p_ascending"`・`notes`。断らない）。`rule_hash` は正規化した spec の hash
（実行のたびに変わる値を含まない）、`matrix_id` は前処理の入力と設定の hash（同じなら同じ ID）。

**`analysis-result.v2` の全量（`output_path` の JSON）。** `{schema, spec（正規化後）, rule_hash, input_fingerprint,
matrix_id, preprocess: {steps, warnings, subject_source}, primary, analyses: [{analysis_id, kind, result}]}`。`result`
は kind ごとに次のスキーマ。

- **`stats-compare.v1`**（compare）: `{matrix_id, test（実際に使った検定名）, effect_size_definition: {name, description},
  ci_level, p_adjust_method, verdict_rule: {q_threshold, equivalence_fc, equivalence_log2}, bh_population: {name,
  n_by_contrast, p_adjust_method, q_interpretation}, presence, prior, design, group_n, subject_source, confounding,
  warnings, rows, omnibus, factorial}`。
  `q_interpretation` は `discovery_fdr`（`p_adjust=bh` で `matrix` / `annotated`）・`fwer`（`bonferroni` / `holm`）・
  `unadjusted`（`none`）か `reproduction_only`（`reported`。報告した特徴だけの補正で、探索の FDR ではない）。
  `omnibus` は cell が 3 つ以上のときの全 cell の F 検定の行 `{feature_id, f, f_p, f_q, df1, df2}`（無ければ null）。
  `df2` は `df_prior + df_residual` で、事前分布の自由度が ∞ なら `"inf"`（JSON。F の参照分布は χ²(df1)/df1 の極限）、
  F を計算していない行では null。
  `factorial` は因子が 2 つ以上で全 cell がそろうときの主効果と交互作用の F（Type III。多水準を含む moderated F）で、
  `{effects: {効果の名前: {factors, df1}}, hypothesis, rows: [{effect, feature_id, f, f_p, f_q, df1, df2}]}`（無ければ null。
  効果の名前は因子名、交互作用は `:` で結ぶ。`f_q` は効果ごとに `p_adjust` で補正）。
  `confounding` は `{batch_table, batch_status, cramers_v, injection_order, warnings, by_contrast: {対比 ID: {batch_status,
  confounded}}}`（全対比をまとめた表と、対比ごとの表。`confounded` は `true`＝完全交絡 / `false` / `null`＝バッチ不明で
  評価できない。交絡なしではない）。
  行（対比 × 特徴）: `feature_id`・`contrast_id`・`reference`・`test`・`status`（`tested` / `not_tested`）・`reason`・
  `n_reference` / `n_test`（有限値の個体数）・`detected_reference` / `detected_test`（検出された個体数。検出状態が
  無ければ null）・`quantity_judgeable`（群のどの cell も検出の割合が `min_detected_fraction` 以上か。null は検出状態なし）・
  `constant_in_group`（群の cell の中で値がすべて同じ＝埋め値の疑いの cell の一覧。無ければ null。判定は変えず、件数を警告
  `constant_in_group` に出す）・`confounded`（その対比が群とバッチの完全交絡か。評価できなければ null）・
  `log2fc`（正なら test が高い）・`se`（log2fc の標準誤差。順位の検定・算術平均の比では null）・`ci_low` / `ci_high`・
  `ci_unavailable_reason`・`statistic` / `df` / `p`・`mean_reference` /
  `mean_test`（個体ごとの前処理後の線形値の算術平均。単純な 2 群の対比だけ、交互作用は null）・`q`・`presence_p` /
  `presence_q`（検出の有無そのものの検定）・`verdict`・`verdict_reason`。
- **判定（`verdict`）は 3 値**で、規則は機械的:
  `changed`＝検定でき、`q ≤ verdict.q`。`unchanged`＝`q` が閾値を超え、95% CI が同等の範囲
  （`|log2fc| < log2(equivalence_fc)`）に完全に収まる。それ以外はすべて `inconclusive`（`verdict_reason`:
  `not_tested` / `quantity_not_judgeable` / `outside_bh_population` / `no_ci` / `neither`）。
  **「有意でない」は「変わっていない」ではない**——`q` が大きいだけでは `unchanged` にならず、差が小さいと言える
  CI があるときだけ `unchanged`。`inconclusive` は「差があるともないとも言えない」で、n が小さい・ばらつきが大きい
  ときの大半はここに入る。この判定は spec に書いた規則を当てた結果で、論文を再現したかの判定ではない。
- **`stats-pca.v1`**（pca）: `{matrix_id, scaling, log2: true, imputation, n_fit, unit_ids, qc_ids, scores: {units, qc},
  explained_variance_ratio, loadings, kept_features}`。軸は生体の個体だけで決め、QC は同じ軸へ**投影**する（軸の計算に
  使わない）。欠測は PCA のためだけに補い（`imputation`）、行列そのものは変えない。
- **`stats-plsda.v1`**（plsda）: `{matrix_id, factors, classes: {cell: 個体数}, y_coding, n_units, unit_ids, scaling,
  imputation, n_features, n_features_dropped_constant, n_components_requested, n_components, r2x_cum, r2y_cum, q2_cum, q2,
  cv: {method: "stratified_kfold", folds, seed, scaling_within_fold: true, q2_definition: "1 - PRESS/TSS"},
  permutation: {n, n_finite, p, statistic: "q2", q2_null_mean, q2_null_q95}, vip: [{feature_id, vip}]（降順）,
  scores: {unit_ids, values}, warnings}`。**成分数・交差検証の Q2・並べ替え検定の p を必ず返す**（PLS-DA は群の札で
  軸を選ぶので、当てはめの R2Y やスコアの分離だけでは過学習と区別できない）。Q2 は層化 k 分割で、分割ごとに学習の個体
  だけで中心と尺度を決めて予測した `1 − PRESS/TSS`（`q2_cum` は成分 1〜A それぞれ、`q2` は A 成分）。並べ替え検定は
  クラスの札を入れ替えて同じ分割で Q2 を出し直した帰無分布で `p = (1 + #{Q2_perm ≥ Q2}) / (1 + n)`。VIP は全個体の
  当てはめの変数重要度（二乗の平均が 1。目安の 1 は経験則で、検定ではない）。Y は 2 クラスなら後のクラスを 1 とした
  1 列、3 クラス以上はクラスごとの指示変数。全個体で値の変わらない特徴は外す。成分数・分割数を個体数と特徴数で
  頭打ちにしたら警告（`n_components_capped` / `cv_folds_capped`）、Q2 が正でなければ `plsda_q2_not_positive`（予測
  できない模型の VIP は判別の根拠にならない）。値の欠けた個体は外す（`units_missing_model_values`）。エラーは
  `PLSDA_SPECIFICATION_INVALID`（引数）・`STATISTIC_GROUP_TOO_SMALL`（2 個体以上のクラスが 2 つ未満）・
  `PLSDA_MISSING_VALUES`（`impute=none` で欠測あり）・`PLSDA_TOO_FEW_FEATURES`（値の変わる特徴が 2 未満）。
- **`stats-composition.v1`**（composition）: `{class_label_kind, classes, unclassified_features, zero_handling, clr, absolute,
  unit, warnings}`。`clr` / `absolute` はそれぞれ stats-compare.v1。クラスの合計が 0 以下の個体があるクラスは `clr` から
  除く（`zero_handling.excluded_classes`）。欠測のメンバーは合計に 0 として入り、その件数と警告
  `class_total_includes_missing_members` を出す。
- 中の `warnings[]` は `{code, message, ...}`。長い一覧（`samples` など）は要約で切る。

対話の `dataset_compare` / `dataset_pca` が保持する結果は、同じ stats-compare.v1 / stats-pca.v1 を
`ds.results[result_id]` の記録（`kind` が `compare` / `stat_pca`）の `engine` に持つ。これらの対話の要約と、記録（`reports/stats/`）の
形は `dataset.md` の §10・§11。

## 3. 行列の recipe（matrix-recipe.v2）

`steps` に挙げた段だけを、挙げた順に走らせる（挙げない段の設定は効かないので捨て、`notes` に書く。
正規化した recipe にも `rule_hash` にも入らない）。篩の段
（`detection_filter`・`qc_rsd_filter`）で落とした特徴はその直後に行列から外す——後の KNN の
近傍・PQN の参照に混ざらない。既定（正規化で埋まる）:

| 段 | 既定 | 定義 |
|---|---|---|
| `detection_filter` | `{"basis": "detected", "min_rate": {}}` | 役割（`qc`・`sample`・`blank`）ごとの検出率の最小。挙げた役割のすべてを満たす特徴を残す。`basis`: `detected` は gap-fill を除く実検出（検出状態が無ければ `DETECTION_UNAVAILABLE`）、`positive` は値が有限かつ > 0（gap-fill を含む） |
| `impute` | `{"method": "half_min", "missing": "zero", "k": 10, "scale": "log2", "rowmax": 0.5, "colmax": 0.8}` | 欠測（`missing`: `zero` は ≤0 か非有限、`nan` は非有限だけ）を埋める。`half_min`: 特徴ごとの正の値の最小の半分（正の値が無い特徴は欠測のまま）。`knn`: 下記。`none` |
| `normalize` | `{"method": "pqn", "reference": "qc"}` | `pqn`: 参照（`reference` の役割の試料の、特徴ごとの正の値の中央値）に対する商の中央値で割る。`tic`（有限値の和）・`median`（正の値の中央値）は係数を試料間の中央値で割る。係数が 0・非有限の試料は割らずに報告する |
| `drift_correct` | `{"method": "loess", "span": 0.75, "iterations": 0, "min_qc": 4}` | 特徴ごとに、正の QC 値を注入順に対して LOWESS で当て、`値 / 当てはめ × QC 値の中央値` に直す。正の QC 値が `min_qc` 未満、または当てはめが正の有限でない特徴は直さず数える。注入順はすべての試料に要る（`INJECTION_ORDER_MISSING`） |
| `qc_rsd_filter` | `{"max_rsd": 0.3}` | QC の有限値の `std(ddof=1)/mean` が `max_rsd` 以下の特徴を残す。RSD が出ない特徴も落とす（別に数える） |

**参照する標準の定義:**

- KNN は R の `impute::impute.knn`（Bioconductor）の `knnimp` の規則を写したもので、`maxp` =
  特徴数（2 段の分割をしない）に当たる。特徴を行に、`scale="log2"` なら正の値の log2 の上で:
  欠測の割合が `rowmax` を超える特徴（重い行）は**重い行どうしの**列の平均で埋め、近傍にしない
  （重い行の中にその列の観測が無ければ軽い行の列の平均）。残りの特徴（軽い行）は軽い行だけで、
  欠測のある特徴ごとに、その特徴が観測された列のうち相手も観測された列での平方差の平均を距離に
  し（共通の列が無い相手は除く）、近い `k` 個を選ぶ（同距離は特徴の並び順）。欠測の列は近傍の
  うちその列を観測しているものの平均、全員が欠測なら**軽い行の**列の平均。`colmax` を超える
  試料があれば `RECIPE_PRECONDITION`。R の既定（`maxp=1500`）は特徴が多いと k-means で
  分割してから補完するので、特徴が 1500 を超えるデータでは R の既定の結果とは一致しない。
  **R の実出力との突き合わせはまだしていない**（手元に Bioconductor の `impute` が無い。規則を
  写し、手計算の例でテストしている）。
- LOWESS は statsmodels 0.15 の
  `statsmodels.nonparametric.smoothers_lowess.lowess(endog=QC 値, exog=QC の注入順,
  frac=span, it=iterations, delta=0, xvals=全試料の注入順)`（1 次の局所重み付き回帰・
  tricube の重み。`iterations=0` は頑健化の反復なし）。R の `loess`（既定は 2 次）とは別物。
  QC が少ないと当てはめが NaN になる点があり、その特徴は補正しない。

## 4. MS/MS の合否（msms-rules.v1）と診断イオン

`library_msms_qualify(rules, spot_id | precursor_mz, ..., mz_tol_ppm)`。照合は `library_match_feature` と
同じ（採点の芯は `library/matching.py`。前駆体の窓は MS-DIAL と同じく m/z 500 を超えたら同じ ppm に広がる。
`mz_tol_ppm`（任意）はその窓に重ねる ppm の絞り込み。使った窓は `query` の `mz_tol`・`mz_window`・
`mz_tol_ppm`。`ion_mode` は `pos` / `neg` / `+` / `-` も読み、読めない値は `DATASET_BAD_REQUEST`）。
規則（正規化後。既定はこのとおり）:

```json
{"schema": "msms-rules.v1",
 "target": {"inchikey_block": "HVYWMOMLDIMFJA", "formula": null},
 "top_rank": 5, "min_fragments_excl_precursor": 1, "min_fragment_relative_intensity": 0.01,
 "precursor_exclusion_tol": null, "same_formula_margin": null, "lipid_class": null,
 "diagnostic_ions": [], "neutral_losses": [], "ions_match": "all"}
```

**2026-11-14 に既定を変えた（事実確認 A1）**: `min_fragments_excl_precursor` の既定を `null` → `1`
（前駆体だけが合った候補を既定で qualified にしない）、`min_fragment_relative_intensity`（0 以上 1 以下。
測定の基準ピークに対する比）を足して既定 0.01 にした。断片の一致は、測定強度が基準ピークの
`min_fragment_relative_intensity` 以上の窓だけを数える（雑音の窓を「一致」に数えない）。既定が
変わったので、同じ規則でも `rule_hash` は前と違う（意図した変化）。`min_fragments_excl_precursor: null`
を明示すれば断片の検査を外せる。`matched_peaks_count`（MS-DIAL の写し）の意味は変えていない。

入力の `target` は `{"inchikey": <InChIKey 全体か先頭 14 文字>, "formula": <任意>}`。
InChIKey が一意に決まらない脂質の種は `{"lipid_species": "PC 18:3_16:1", "formula": <任意>}`
（`inchikey` と排他。opt-in）で、正規化後は `{"lipid_species": "PC 34:4", "formula": …}`（GOSLIN の
種のレベル。読み方は `lipidmix://docs/output-format/mztab` §13.8）。種の目標では、下の
`target_in_top_rank` と `same_formula_margin` の「目標か別物か」をライブラリの記録の種
（`msms_rules.record_lipid_species`）で決める——同じ種の別の分子種（`PC 14:0_20:4`）は目標の候補で、
別物には数えない。記録の種は次の順に読み、最初に読めた経路を `species_route` とする:

| `species_route` | 読み方 |
|---|---|
| `name` | 記録の `name`（MS-DIAL の限定子を外して）をそのまま GOSLIN で読む |
| `name_cleaned` | `\|` の後・`; ` の後の注記（`PC 34:1; ExampleDB ID-0000`）・末尾の付加体の表記（`[M+H]+`）を外し、古い MS-DIAL の書式（`Cer[NS] d18:1/16:0`）・英語のクラス名（`Phosphatidylcholine 16:0-18:1`）を略記に直して読む（`pathway.msp_lipid.shorthand`） |
| `ontology_formula` | 名前で読めなければ、ClassyFire の `Ontology` と分子式から和の組成の種を作る（`pathway.msp_lipid.from_class_formula`。LPS・LPG も） |

同位体の標識のある名前（`(d7)` など、内部標準）と、どの経路でも読めない記録は目標とも別物とも
数えない。目標の候補の経路は `target_in_top_rank` の `detail.species_route`、候補の表では
`lipid_species`・`species_route` の列に出る。
設定した検査が `checks` に並ぶ（各行 `{check, status: pass|fail|not_evaluable, value, detail}`。
既定の規則では `target_in_top_rank` と `min_fragments`）:

| 検査 | 合格の条件 | `value` |
|---|---|---|
| `target_in_top_rank`（常に） | 総合スコアの順で `top_rank` 位までに、InChIKey 先頭ブロックが目標と一致する候補がある（最上位のそれを「目標の候補」） | 目標の候補の順位（無ければ null） |
| `min_fragments`（既定で入る） | 目標の候補の、前駆体 ± `precursor_exclusion_tol`（省略時は採点の `ms2_tol`）を除き、測定強度が基準ピークの `min_fragment_relative_intensity` 以上の一致の本数 ≥ 最小（`detail.min_relative_intensity` に下限） | その本数 |
| `same_formula_margin` | 全候補のうち、目標と同じ分子式（`target.formula`、無ければ目標の候補の分子式）で先頭ブロックが違うものの最高点との差 ≥ `min_delta`（`score` 既定 `total_score`）。別物が無ければ pass | 差（目標 − 別物） |
| `lipid_class` | 1 位の候補の名前を GOSLIN で読んだクラス（正規化名の最初の空白の前）が一致 | 読めたクラス |
| `diagnostic_ions` / `neutral_losses` | 測定スペクトルに ±`tol`（既定 0.01 Da）のピークがあり、最大ピークに対する相対強度 ≥ `min_relative_intensity`。`ions_match` が `all` なら全部、`any` なら 1 つ。中性損失は前駆体 − `delta`（生成イオン）を探し、`min_relative_intensity` はその生成イオンの比 | 見つかった数 |

要素は診断イオンが `{mz, label, tol, min_relative_intensity}`（`min_relative_intensity` 既定 0）、中性損失が
`{delta, label, tol, min_relative_intensity}`（中性損失の `min_relative_intensity` は渡したときだけ規則に載る——
省いた規則の `rule_hash` は変わらない）。知らないキー・欠けたキーは `MSMS_RULES_INVALID` で、誤りの文が
そのキーと要素（`neutral_losses[0]` など）を名指す。

戻り値は `verdict`・`checks` のほかに、`target_rank`（目標の候補の順位。無ければ null）と
`top1_is_target`（1 位が目標か）を返す。既定の `top_rank` 5 では、目標が 2 位でも `target_in_top_rank` は
pass するので、1 位かどうかはこの欄で読む。`spot_id` のスポットに MS-DIAL の注釈があれば、
`library_match_feature` と同じ `msdial_name`・`msdial_rank`・`top_differs_from_msdial`（library.md §14.6）。
候補の表（`candidates`、TSV）の列は `rank name formula inchikey adduct delta_ppm total_score
matched_peaks_count fragment_matches_excl_precursor`（`delta_ppm` は `(測定 − 参照) / 参照 × 1e6`、
`fragment_matches_excl_precursor` は合否と同じ規則＝強度の下限つきで数える。`library_match_feature(detail=true)` の
同じ名前の列も既定は同じ下限 0.01 で、使った下限を `min_fragment_relative_intensity` に返す）。`.dbs` に注釈の足切り
（dot・一致率・総合スコア・`minimum_spectrum_match`）があれば、値を `dbs_cutoffs` に 1 回載せ、
候補ごとに `passes_dbs_cutoffs`（MS-DIAL の `IsSpectrumMatch` と総合スコアの足切りを全部通るか。
dot の足切りは二乗の値と比べる）の列を足す——**合否には当てない材料**（`.msp` には足切りが無いので列も無い）。
種の目標なら `lipid_species`・`species_route` の列も足す。

`verdict`: すべて pass → `qualified`、どれか fail → `not_qualified`、fail が無く判定できない
検査がある → `not_evaluable`。目標の候補が無いと、それに依る検査は `not_evaluable`
（`target_in_top_rank` が fail なので `verdict` は `not_qualified`）。MS/MS が無い・候補が無い
ときは検査をせずに `not_evaluable`（「断片が合わなかった」ではない）。**「明らかに別物か」は
判断しない**（再解析の計画の規則で決める）。`rule_hash` は正規化した規則の hash。

`dcl_diagnostic_ions(spot_ids, ions)` は同じ照会（`find_ions`）をスポットの一覧に当てる。
`ions` の各要素は `{"mz"}` か `{"neutral_loss"}`（任意で `label`・`tol`・
`min_relative_intensity`。中性損失の比は生成イオン＝前駆体 − `neutral_loss` の測定強度の比。知らない
キーは `MSMS_RULES_INVALID` で、誤りの文が `ions[i]` とそのキーを名指す）。行の `status` が `ok` 以外（MS/MS が無い・索引の崩れ）のスポットは
照会していない——`found` が空なのは「イオンが無かった」ではない。

## 5. 同族列の RT（homolog-check.v1）

`dataset_homolog_check(spec)`。`members`（`feature_id`・`carbons`・`double_bonds`）か `auto`
（`lipid_class`・`evidence` 既定 `["sme"]`・`adduct`）のどちらか一方と `rule`
（`carbon_direction` 既定 `increasing`・`double_bond_direction` 既定 `decreasing`・
`max_violations` 既定 0。逆相の LC の向き）。

- `auto` は SME rank 1 の名前（`sml` を足せば SML の名前も）を GOSLIN で読み、分子種
  （`16:0_18:1`）は鎖を足して種のレベル（`34:1`）に直す。エーテル（`O-` / `P-`）は別の
  クラス（`PC O-`）として扱う。
- 同じ（炭素数, 二重結合数）が複数あれば RT の中央値で並べ、`duplicates` に挙げる。
- 同じ二重結合数の中で炭素数の昇順、同じ炭素数の中で二重結合数の昇順に並べた**隣どうし**で、
  規則の向きに反する対を違反と数える（RT が同じ対は違反にしない）。違反が
  `max_violations` 以下なら `pass`。違反に関わった特徴は `members` の `status` が `violation`。

## 6. 特徴どうしの関係

`dataset_feature_relations(feature_id, rt_window=0.05, mz_tol=0.01, min_r=0.7, adducts=None)`。
規則はアラインメントのキュレーションの関係（`curation/relations.py`）と同じ。

- `coeluting`: |ΔRT| ≤ `rt_window` の特徴（`drt` = 相手 − この特徴、`r` = 試料間の log1p 強度の
  Pearson 相関。値が 5 組未満なら空）。
- `as_ion_of`: この特徴が、相手の付加体（`adducts` の組）・同位体（M+1、M+2 は 13C2・37Cl・81Br）・
  源内断片（`-H2O` ほか）として説明できる関係。列は `of of_name relation dmz_mda drt r isotope_ratio
  expected_ratio isotope_consistent strong ambiguous informational`。
- `ions_of_this`: 相手がこの特徴のイオンとして説明できる関係。列は `ion relation dmz_mda drt r
  isotope_consistent strong ambiguous informational`。付加体の関係はこの特徴の付加体から計算するので、
  **注釈の無い相手も付加体・同位体・源内断片として出る**。
- 付加体の関係は説明される側の付加体から計算する。その付加体は SME rank 1、無ければ SML 注釈（MS1 だけの
  同定）のもの。どちらも無い相手の `as_ion_of` には付加体の関係は出ない（同位体・源内断片だけ）。
  同位体の期待比に使う組成式も、SME の行 → SML 注釈 → 同定のある特徴の `formula`（解析フォルダの `.arf2`）の順に採る。
- **`strong` は Δm/z が合うことだけでは立たない**。MS-DIAL のリンク（`.arf2`）は使わないので、裏付けは
  試料間の相関 `r ≥ min_r` だけ。付加体・源内断片は Δm/z と同時溶出（`rt_window`）に加えて `r ≥ min_r`
  のとき、同位体はさらに強度比が組成式の期待比の 1.5 倍以下（`isotope_consistent=true`）で、M+2 の
  取り違えの恐れ（`ambiguous`）が無いときだけ `strong`。
- `ambiguous`: 同じ相手に 2 つ以上の M+2 の説明（13C2 +2.006710・37Cl +1.997050・81Br +1.997952）が
  同じ窓に入った行。一方に決めず並べ、どれも `strong` にしない。相手の組成式が分かっていて Cl・Br を
  持たなければそのハロゲンの行は出さない（脂質の M+2 は 13C2 だけ）。
- `informational`: 同位体の強度比が期待比を超えた行（その強度は同位体だけでは説明できず、別物質が
  重なっている目印）。`ambiguous` の行は、すべての説明が期待比を超えたときだけ立つ——13C2 の期待比に
  合わないことだけで「別物質」と言わない。

### 6.1 試料ごとのピーク（`dataset_sample_peaks`）

`dataset_sample_peaks(feature_ids, dataset_handle=None, output_path=None)`: 割れた 2 スポットが 1 つの物質か、
RT が注入ごとにずれていくかを、試料ごとの値で見るための読み取り。行は列 `feature_id, sample, rt, mz, height,
detected` の **TSV**（列名 1 回、`rows` の文字列。依頼の特徴の順 × データセットの試料の順）。戻り値は
`{status, n_features, n_samples, n_rows, rows}`。数値は小数 6 桁まで、空欄が null（gap-fill の `mz`）、
`detected` は `true` / `false`。

- `rt`: その試料のピークの RT（分）。gap-fill のセルも `.arf` が持つ値（補間の位置）を返すので、
  **`detected: false` の `rt` は検出された位置ではない**。ドリフトは `detected: true` の行だけで見る。
- `mz`: 検出されたセルの実測 m/z。gap-fill は空欄（null）。
- `height`: データセットの強度の行列（`feature_matrix`）の値。`measure` で選んだ `.arf` の欄で、定量
  （`dataset_quant_apply` ほか）で値を変えたデータセットではその値。
- `detected`: `detected_mask`（gap-fill なら false）。
- 特徴の数に上限は無い。知らない特徴 ID・空のリストは `DATASET_BAD_REQUEST`。ID は int でも str でもよい。
- 戻り値に載せる行（特徴 × 試料）は **1,000 行まで**。超えるとき `output_path` が無ければ**断らず**、依頼の順の
  先頭の特徴の行を特徴ごとにそろえて（1 つの特徴の試料を途中で切らない）floor(1000 / 試料の数) 特徴分だけ載せる
  （試料が 1,000 を超えれば最初の特徴の先頭 1,000 行）。戻り値は `truncated: true`・`n_rows`（全件）・
  `n_rows_shown`・`n_features_shown`・`notes`（`truncated` はいつも付く。超えなければ false）。
  `output_path` があれば全量の TSV（同じ列・列名 1 回）をそのファイルへ原子的に書き、戻り値は要約
  `{status, n_features, n_samples, n_rows, output_path, sha256}`（`rows` は付かない。ほかの再解析のツールの
  `output_path` と同じ扱い。1,000 行以下でも `output_path` を渡せば同じく要約）。
- 試料ごとの RT（`per_sample_rt`）のあるデータセットだけ: 解析フォルダ（`.arf`）から読んだもの、または mzTab-M で
  読み、隣の `.arf` と同じアラインメントだと照合できたもの（検出状態を取り込めたとき。`.arf` の列を mzTab-M の試料の順に
  並べ替えて持つ）。照合できなかった mzTab-M と、統合などで作り直したデータセットは `DATASET_BAD_REQUEST`
  （`field: per_sample_rt`）。

## 7. 証拠の束（claim-evidence.v1）

`dataset_claim_evidence(spec, dataset_handle, output_path)`。1 主張について、オブザーバーの
機械的な仕事（走らせ直しと照合、候補の網羅の走査）を道具が 1 つの JSON にまとめる。

spec:

```json
{"claim_id": "C1",
 "target": {"inchikey": "...", "formula": "...", "rt_min": 8.4},
 "scan": {"ppm": [10, 20], "rt_tol_min": [0.3, 1.0], "adduct_tiers": null},
 "msms": {"top_n": 10, "rules": {"min_fragments_excl_precursor": 2}, "mz_tol_ppm": 10},
 "relations": {"rt_window": 0.05, "mz_tol": 0.01, "min_r": 0.7},
 "max_candidates": 50}
```

- `ppm_offset`（任意）は質量の窓の中心のずれ。形と帯の規則は `dataset_feature_search` と同じ
  （`{"default": 0.0, "bands": [{"max_mz": 200, "offset": -8.0}]}`。中心は実測 − 理論の ppm、帯は理論 m/z が
  `max_mz` 未満の最初の帯）。既定の特徴と `scan` の全変種に同じ spec が効き、**窓の幅（`scan.ppm`）は変えない**。
  渡したときだけ `rule.ppm_offset` と `rule_hash` に入る（省けば v8 と同じ hash）。形の誤りは `CLAIM_SPEC_INVALID`。
- `rt_offset_min`（任意、分）は RT の窓の中心のずれ。意味は `dataset_feature_search` の `rt_offset_min` と同じで、
  中心を `rt_min + rt_offset_min` に置き `|RT − 中心| ≤ rt_tol_min` で見る（`mztab` §13.8）。既定の特徴と `scan` の全変種に
  同じ値が効き、**窓の幅（`scan.rt_tol_min`）は変えない**。渡したときだけ `rule.rt_offset_min` と `rule.target.rt_offset_min`、
  `rule_hash` に入る（0 でも入る。省けば従来と同じ hash）。`target` に `rt_min` が無ければ RT の窓が無く効かないので、
  誤りにせず捨てて `notes`（束と戻り値）に書く（規則にも hash にも入らない）。有限の数でなければ `CLAIM_SPEC_INVALID`
  （`field: rt_offset_min`）。ずれの見積もりは `dataset_feature_search` の目標照合が返す `rt_offset`（当てはしない）。
- `target` は `dataset_feature_search` の目標と同じ形（`target_id` は `claim_id`）。脂質の種の主張
  （InChIKey が無い C1〜C4 型）は `"target": {"lipid_species": "PC 18:3_16:1", "formula": …}`。
  `rule.target` に `lipid_species`（種のレベル）と `lipid_species_input` が載り、`rule_hash` は
  `lipid_species_input` を入れない（綴りの違う同じ種は同じ hash）。**束の既定の特徴は種のレベルの
  一致で選ぶので、分子種（どの鎖の組か）は確かめていない**。
- `scan` の `ppm` × `rt_tol_min` の全組で目標照合を走らせ直す（`adduct_tiers` を渡すと
  付加体も広げる）。省略時は目標の `ppm` / `rt_tol_min` の 1 組。
- `msms.rules` は msms-rules.v1 から `target` を除いたもの（target は主張の InChIKey か脂質の種と分子式）。
  `msms.mz_tol_ppm`（任意、0 より大きい数）は前駆体の窓に重ねる ppm の絞り込み（`library_msms_qualify` と同じ）。
  渡したときだけ `rule.msms` と `rule_hash` に入る。形の誤りは `CLAIM_SPEC_INVALID`。
- `msms` / `relations` を省くと、その欄を作らない。

結果:

| 欄 | 意味 |
|---|---|
| `default_selection` | `target` をそのまま当てた目標照合（`dataset_feature_search` と同じ規則）。`status` が `resolved` なら `selected_feature_id` が規則で選んだ既定の特徴、決まらなければ `needs_input`（`reason`）。割れたスポットの組（§9 の `merge_groups`）は 1 つの候補として数え、`selected_feature_id` は残す ID（受け持ちの多い側、同数なら小さい ID）、`merged_from` にその組（それ以外は null。走査も同じ）。`counts` は下の表（組は 1 つと数える。MS/MS の合否はどちらかのスポットが `qualified` なら 1） |
| `scans` | 走査ごとの判定（`counts` は `default_selection` と同じ意味）と、候補の特徴 ID |
| `candidates` | 全走査の候補の和（特徴の並び順）。`evidence`（`sme_match` など、走査で見えたもの）・`adducts`・`scans`（どの走査で見えたか。組の片方だけが質量の窓に入ったときに足した相手は `merge_partner`）・`split`（下）・`msms`（spot_id = SMF_ID の照合の上位 `top`（`library_msms_qualify` の候補の表と同じ欄。`passes_dbs_cutoffs`・種の欄も同じ条件で）と `qualify`（`target_rank`・`top1_is_target` を含む）、使った前駆体の窓 `window`（`mz_tol`・`mz_window`・`mz_tol_ppm`）、スポットに MS-DIAL の注釈があれば `msdial_name`・`msdial_rank`・`top_differs_from_msdial`。ライブラリが無ければ `{"status": "no_library", "spectrum_available": ...}`（データセットの同定が LBM だけ＝記憶が `lbm_only` のときは、MSP のスペクトルが無いので評価できない: 候補ごとは `{"status": "not_evaluable", "reason": "lbm_only", "spectrum_available": ...}` で、理由の文は束の `msms_not_evaluable_reason` に 1 回だけ。`no_library` は `library_load` を促すが、このデータセットには照合に使える MSP が無い）、MS/MS が無い・索引の崩れはその `status`。`spectrum_available` は代表スペクトルの有無で、ライブラリの有無によらず spot を先に照会して決める: `true`=ある、`false`=照会して無かった（`not_found` で `reason` が `no_spectrum`）、`null`=分からない（`.arf2` / `.dcl` のファイル自体が無い `not_found` の `reason` `no_arf2` / `no_dcl`・`spot_not_found`・`dcl_precursor_mismatch`・`spot_id_not_integer`、ハンドルに `.arf2` / `.dcl` が無いとき））・`relations` |
| `truncated` | `max_candidates` で打ち切った候補の数（束の欄。`output_path` を省いて候補を切った戻り値だけは意味が違う——下の「戻り値」） |
| `notes` | 効かずに捨てた spec の設定とその理由（`target` に `rt_min` が無いときの `rt_offset_min`）。無ければ欄ごと無い |
| `rule` / `rule_hash` | 正規化した spec とその hash。spec の `merge_split`（既定 true）は真のときだけ `rule` に載り、`rule_hash` には判定でまとめた組の hash も入る（false は v7 と同じ規則。hash が v7 と同じなのは spec に `msms` が無いときだけ: `msms` があると、2026-11-14 に MS/MS の規則の既定（`min_fragments_excl_precursor` 1・`min_fragment_relative_intensity` 0.01）を変えたので hash も変わる） |
| `merge_split` / `merge_groups_hash` | 割れの組を判定で 1 つの候補にしたか（spec の `merge_split`）と、その組（§9 の `merge_groups`）の hash（false なら null） |
| `input_fingerprint` | `dataset_hash`（データの指紋）・`library_sha256` |

`counts`（`default_selection` と `scans[]`。どれも特徴の一意な件数）:

| 欄 | 意味 |
|---|---|
| `mass` / `rt` / `by_evidence` | `dataset_feature_search` の `counts` と同じ（質量の窓・RT の窓に入った数、水準ごとの規則の qualified 数） |
| `qualified` | 規則（質量・RT・注釈）で qualified になった特徴のうち、**MS/MS の合否 `msms.qualify.verdict` が `qualified`** の数。`msms` を頼んでいない・ライブラリが無いときは `null`（0 ではない）。`max_candidates` で打ち切った候補は数えない |
| `rule_qualified` | 規則（質量の窓・RT の窓・注釈）で qualified になった特徴の数（`dataset_feature_search` の `counts.qualified`、v7 までの束の `counts.qualified` と同じ数）。**MS/MS が取れているかどうかの数ではない** |
| `msms_available` | **非推奨の別名**（`rule_qualified` と同じ値。claim-evidence.v1 の互換のため残し、次の版で消す）。名前に反して MS/MS の有無の数ではない |
| `msms_spectrum_available` | 質量の窓に入った候補（`mass` と同じ母集団。`merge_partner` の行を除き、割れの組は 1 つ）のうち、組のどれかのスポットに**代表 MS/MS がある**（候補の `msms.spectrum_available` が `true`）数。RT の窓では絞らない（RT の窓の外のスポットの代表 MS/MS も数える）。**0 は照会して無かった**（全部の単位が `false`）、**`null` は頼まなかった（`msms` なし）か分からなかった**（`true` の単位が無く、`null` の単位が 1 つでもある）、**数は下限**（`null` の単位は数に入らない）。単位の中は、どれか `true` なら `true`、`true` が無く `null` があれば `null`、全部 `false` なら `false`。`max_candidates` で打ち切って照会していない候補は数えない |

`split`（候補が割れたアライメントスポットの組に入るときだけ。入らなければ `null`。§9 の
`find_split_pairs` を既定の閾値で回し、`split_partners` を写す）:

| 欄 | 意味 |
|---|---|
| `group_id` | 全体の連結成分（§9 の `groups[].group_id`） |
| `partners` | まとめる相手（§9 の `merge_groups`。互いに最も近い mergeable な組の相手 1 つ）。まとめない特徴（`group_aligned`・最も近くない組など）は空 |
| `type` | まとめる特徴はその組の型、まとめない特徴は成分の型（`group_aligned` > `interleaved` > `injection_block` の順で強い方） |
| `mergeable` | 統計の側（`dataset_feature_search` / `dataset_analysis_run` の `merge_split`）が 1 つにまとめるか |
| `component` | 成分のほかの特徴（報告用。まとめない割れの相手もここで分かる） |

**束は特徴を決めない。** 最終の特徴は再解析の実行役が override の規則の範囲で決め、
その特徴での統計は `dataset_analysis_run`（`features.feature_ids`）で別に出す。同じ spec・
同じ入力なら同じ JSON になる（オブザーバーの「走らせ直しと照合」は JSON の一致で済む）。

**戻り値。** `output_path` を省けば束をそのまま返す。候補が 10 を超えるときも**断らない**: 束の候補を
`candidates_order: "default_selection_first"` の順——既定の特徴（`default_selection` の `selected_feature_id` と
`merged_from`）→ 既定の照合で見えた候補（`scans` に `default`）→ ほかの候補、同じ段の中は束の並び（特徴の並び順）——に
並べた先頭 10 件だけを `candidates` に載せ、`truncated: true`・`n_candidates`（束の候補の全件＝`output_path` が書く数）・
`n_candidates_shown`・`candidates_order`・`notes`（`output_path` で全件を書ける旨）を返す。このとき束の `truncated`
（`max_candidates` で打ち切った数）は `max_candidates_truncated` に移る（0 でなければ `notes` にも書く。打ち切った候補は
`output_path` の JSON にも入らないので、要るなら `max_candidates` を増やす）。候補の一覧のほか（`default_selection`・
`scans` の `candidate_feature_ids` を含む判定・`rule`・`rule_hash`・`input_fingerprint`）は束と同じ。`output_path` が
あれば JSON を書き、要約（`claim_id`・`rule_hash`・`input_fingerprint`・`default_selection`・`n_candidates`・`truncated`
（束と同じ `max_candidates` の打ち切りの数）・走査ごとの判定、捨てた設定があれば `notes`）だけを返す。

## 8. 質量精度と新しい知見の候補

`dataset_mass_accuracy(band_edges=[200], rank=1)`: SME の `rank` 位の候補の
(実測 − 理論)/理論 × 10⁶ を、理論 m/z の帯ごと（`[min_mz, max_mz)`）にまとめる。
`n`・`median`・`q1`・`q3`（`statistics.quantiles(n=4)` の exclusive 法）・`within_5ppm`・
`within_10ppm`（|ppm| ≤ 5・10 の割合）。n<2 の帯は `n` だけ。**許容幅は当てない**が、材料として勧める幅を
添える: 帯ごとの `recommended_mz_tol_ppm` と、全体の `recommended_tolerance`
（`{n, basis, mz_tol_ppm, median_ppm, mad_ppm}`）。幅は `ceil(|中央値| + 3 × MAD)`（MAD は中央値からの絶対偏差の
中央値で、尺度の補正 1.4826 は掛けない。1 ppm 未満にしない）。中心を 0 に置いた ±窓で系統的なずれと散らばりの
両方を覆う幅なので、中心をずらすなら中央値を `ppm_offset` に、`3 × MAD` を幅に分けて渡してよい。同定が 10 件未満
なら `mz_tol_ppm`（帯は `recommended_mz_tol_ppm`）は null で、`recommended_tolerance` に `reason: too_few` と `min_n`。
どれを当てるかは再解析の規則で決める（`dataset_feature_search` の目標の `ppm`、`library_match_feature` の
`mz_tol_ppm`）。

`subset`（既定 `rank1`）: `rank1` は `rank` 位の SME すべて。`msms_matched` はそのうち、その特徴の
MS/MS が参照スペクトルに照合されたもの（SML の `opt_global_spectrum_matched` が真。mzTab-M の
SML の列は `docs/output_format/mztab.md` §13.3。**研究室 MSP と LBM の同定を区別せず両方を数える**。
偽・不明の特徴は入らない）。照合の旗は rank 1 の同定を表すので、**`msms_matched` は `rank=1` のときだけ**
（それ以外は `DATASET_BAD_REQUEST`、`field: rank`。`subset` の値が不正でも `DATASET_BAD_REQUEST`）。
SML にその列が無い mzTab では**エラーにせず**、`subset_available: false` と `note` を付けて
`overall` `{"n": 0}`・帯の `n` 0 を返す（手順の後退——帯の n 不足 → 極性全体 → `rank1`——を止めない。
`subset_available` は列が無いときだけ出る）。戻り値に `subset` を書く。

`injection_trend`（選んだ `subset` の SME について）: 試料ごとの質量誤差の、注入順との関係。
**全試料を数える**（ブランク・QC があれば含む）。試料の誤差は、**その試料で検出された（gap-fill でない）セル**の実測 m/z と理論 m/z の差
（ppm。`ds.per_sample_mz`。隣接 `.arf` を取り込めたときだけある）。

| キー | 意味 |
|---|---|
| `per_injection` | `[{sample, injection_order, median_ppm, n}]`（注入順の昇順。`n` はその試料で使えた特徴の数。1 つも無い試料は載せない） |
| `spearman_rho` | `injection_order` と `median_ppm` の Spearman の ρ。**`injection_order_source` が `filename`（採取時刻・注入順が読めず、ファイル名の順に落ちた）なら出さず `null`**（ファイル名の順は注入順とは限らない）。データのある試料が 5 未満、または片方が定数でも `null` |
| `rho_unavailable_reason` | `spearman_rho` を出さない理由（`filename_order` / `too_few_samples`）。ρ が出ているとき、片方が定数で `null` のときは無い |
| `n_samples` | `per_injection` の件数 |
| `injection_order_source` | 注入順の出所（`manifest` / `metadata` / `filename`。§9 の `alignment_qc` と同じ規則）。`filename` のとき `per_injection` はファイル名の順で、ρ は無い |
| `warnings` | 注入順の注意（無ければ鍵ごと無い）。`injection_order_may_be_file_order`（`n_samples`・`message`）: 注入順が MS-DIAL の振った値（`.mddata` / mzTab-M の MS:4000089）で、ファイル名の順と一致する。MS-DIAL は利用者が入れない限りファイルの順で振るので、実際の測定順とは限らない。**ρ はそのまま出す**（ρ はファイル名の並びとの相関かもしれない）。`injection_order_source` が `filename` のときは出さない |

試料ごとの m/z が無い（`.arf` を取り込めていない）ときは `{"status": "unavailable"}` だけ。
ρ が 0 に近いことは「ずれが無い」ではない——ずれが注入順に単調でないだけのこともある
（`per_injection` の `median_ppm` を見る）。

`dataset_novel_candidates(spec)`（novel-candidates.v1）:

```json
{"statistics": {"result_id": "res_..."},
 "ppm": {"default": 10.0, "bands": [{"max_mz": 200.0, "ppm": 10.0}]},
 "q_max": 0.05,
 "exclude": {"feature_ids": [], "inchikey_blocks": []}}
```

`statistics` は `{"result_id"}`（`dataset_compare` の 2 群比較。`reference_group` が基準。対比が複数の結果は
`result_id#contrast_id` で 1 つ選ぶ。ID だけなら誤りでなく `needs_input`——`field` は `statistics.result_id`、選択肢の
`value` は `result_id#contrast_id`、`argument` は無い）か
`{"analysis_spec", "analysis_id"}`（**analysis-spec.v1 だけ**とその解析。v2 の spec はここでは通らないので、統計エンジンの
結果は `dataset_compare` の `result_id` で渡す）。SME rank 1 の同定があり、
|ppm| が帯の許容幅以内（`max_mz` 未満の最初の帯、無ければ `default`）で、q < `q_max` の特徴から
除外を外した一覧（q の昇順）。`dataset_compare` の結果は、加えて判定 `verdict` が `changed` の行だけ
（検出が足りず量を判定できない `inconclusive` の行は q が小さくても候補にしない。`q_max` が判定の q の閾値より
大きくても、判定が `changed` でない行は通らない）。analysis-spec.v1 の結果（判定を持たない）は q だけで選ぶ。行は `feature_id name inchikey adduct mz rt ppm group_ref group_test
log2fc_test_over_ref p q direction merged_from`（JSON には `sme_feature_id` も）。論文の同定リストとの
突き合わせはしない（再解析の側）。

**割れたスポットは 1 候補に数える**（spec 2026-10-28 §4.3 (c)）。§9 の `merge_groups` の 2 スポットは
1 つの物質で、`feature_id` は残す ID（受け持ちの多い側、同数なら小さい ID）、`merged_from` は 2 つの ID（残す ID が先頭）。統計はまとめた行列で検定し直した
値（`statistics.result_id` では `dataset_feature_search` の `merge_split` と同じ実装で、組が 1 つでもあれば
すべての特徴の統計をまとめた行列の検定・BH から出す。analysis_spec では spec の `features.merge_split`）。**同定と質量誤差**（`name` `inchikey`
`adduct` `mz` `ppm` `rt`）は SME（rank 1）を持つ方のスポットから取り、その ID を `sme_feature_id` に書く
（両方にあれば残す ID）。`exclude.feature_ids` は組のどちらかの ID に当たれば組ごと外す。
`merge_split=false`（または analysis_spec の `features.merge_split=false`）では、行・戻り値・`rule`・`rule_hash` が
v7 と同じ（`merged_from`・`sme_feature_id`・`merge_split` を出さず、正規化した spec にも `merge_split` を載せない）。
`statistics.result_id` の spec は `{"result_id", "merge_split": true（既定）, "alignment_qc": <JSON のパス、任意>}`
（`alignment_qc` は `dataset_alignment_qc` の `output_path`。省略時はこのハンドルで `dataset_alignment_qc` が最後に
出した結果（§9。データ・試料の表・除外が変わっていないもの）、それも無ければその場で探す。`rule_hash` にパスは
入らない）。`merge_split=false` の `alignment_qc` は効かないので捨てて `notes` に書く（`rule_hash` は書かないときと同じ）。
`statistics.analysis_spec` の効かない設定も、欄の名前の前に `statistics.analysis_spec: ` を付けて `notes` に書く
（analysis_spec では最後の結果は使わない。spec だけで再現する）。
戻り値と JSON の `merge_split` は `{enabled, source, n_groups}`（`n_groups` はまとめた組の数。`source` は
`alignment_qc`（ファイル）・`last_alignment_qc`（ハンドルの最後の結果）・`find_split_pairs`（その場））。

**行の上限。** `output_path` を省けば行を TSV で返し、50 を超えれば q の昇順の先頭 50 行だけを載せる
（`truncated: true`・`n_rows` が全件・`notes`。断らない）。`truncated` は行を返すときいつも付く。全量は `output_path` の JSON。

**外部の印（`labels`、任意）**。spec に `"labels": {"path": "<TSV>", "name": "xenobiotic", "feature_column": "feature_id",
"label_column": "level", "order": ["A","B","C","D"]}` を足すと、行に `labels.<name>` の列（TSV の末尾。印が無い行は空）を足す。
`path` と `name`（英数字と `_`）は必須、ほかの既定は上のとおり。`order` は強い順で、同じ特徴に複数の行があれば最も強い水準を取る。
空の水準の行は無視、`order` に無い水準は `NOVEL_SPEC_INVALID`（`details.unknown_levels`。黙って最弱にしない）。
**割れの組**（`merged_from`）は組のどれかの最も強い水準。戻り値と JSON の
`labels` は `{name, sha256, order, n_rows, n_unmatched, counts}`（`n_unmatched` はデータセットの特徴 ID に無い、印のある ID の数。無視する）で、
`counts` は `{"all": {水準: 件数, "none": 件数}, "top20": {...}}`（`top20` は q の昇順の先頭 20 行）。
**印は列と件数を足すだけ**で、並び・行数・除外は変えず、`rule_hash` にも入れない（`path` を含む `labels` 全体が hash の材料から外れる）。

**質量の窓の中心のずれ（`ppm_offset`、任意）**。spec に `"ppm_offset": {"default": 0.0, "bands": [{"max_mz": 200, "offset": -8.3}]}`
を足すと、判定が `|ppm − 中心| ≤ 許容幅` になる。中心は「実測 − 理論」の ppm の符号つきの値で、帯は `ppm` の帯と同じ選び方
（理論 m/z が `max_mz` **未満**の最初の帯、無ければ `default`。ちょうど `max_mz` の特徴は次の帯）。**窓の幅（`ppm`）は変えず、
行の `ppm` 列も実測 − 理論のまま**。渡したときだけ正規化した spec と `rule_hash` に入る（省けば v8 と同じ）。形の誤りは
`NOVEL_SPEC_INVALID`（`field: ppm_offset`）。

## 9. 割れたアライメントスポット（alignment-qc.v1）

`dataset_alignment_qc(output_path=None, ppm=15.0, max_rt_gap=1.0, max_interleaved_rt_gap=0.15,
min_block_size=2, dataset_handle=None)`: 1 つの物質が 2 つ以上の
アライメントスポットに割れた組を探す（spec 2026-10-28 §4.3）。**組を見つけて型を付けるだけで、
行列は変えない。** まとめは統計の側が `mergeable` の成分だけに行う（試料ごとに大きい方の値。
gap-fill の値は背景で、本物のピークより大きくならないため）。

- **候補**: 同じ付加体（SML の `adduct_ions`。付加体が分からないスポットは組にしない）・|Δm/z| ≤ `ppm`（m/z の小さい方に対して）・
  0.05 < |ΔRT| ≤ `max_rt_gap` 分の 2 スポット。
- **受け持ち（`owner`）**: 試料ごとに、その試料の信号を持つスポット。検出状態（`.arf` の
  `detected_mask`）があれば検出された方（両方で検出なら値の大きい方）で、どちらでも検出されて
  いない試料は、片方の値が他方の `min_dominance`（3）倍以上ならその方（gap-fill の値で補う。
  検出の疎な run で本物のピークの多くが gap-fill に入るため）、そうでなければ `null`。
  検出状態が無ければ値の大きい方。どちらにも信号（値 > 0）が無い試料は `null`。
- **割れの判定**: 2 スポットがそれぞれ 1 試料以上を受け持ち、かつ補い合う——検出状態があれば
  両方で検出された試料（`both_detected_fraction`。分母はどちらかで検出された試料）が
  `max_both_detected`（0.10）以下で、かつ信号か検出のある試料のうち受け持ちの決まった試料が
  `min_complementary`（0.9）以上。検出状態が無ければ（`both_detected_fraction` は `null`）、片方の値が
  他方の `min_dominance` 倍以上の試料が、信号のある試料の `min_complementary` 以上。
  検出だけで判定しないのは、どちらも数試料でしか検出されない無関係な 2 スポットが「両方で検出
  0%」で割れに見えるから。3 スポットの連鎖は、3 つ目が受け持つ試料で 2 スポットの背景が同程度だと
  その組が成り立たず、隣り合う組だけで成分になる。
- **型**（`type`）:
  - `injection_block` — 片方の受け持ちが注入順で連続した塊（信号の無い試料は飛ばす）。
    両スポットとも `min_block_size`（2）試料以上を受け持つときだけ（1 試料だけの受け持ちは
    どこにあっても「連続」に見えるので塊にせず、`interleaved` か `group_aligned` に回す）。
    群と重なっても優先する（測定順と群が重なる run では、測定順のずれは必ず群とも重なる）。
    両方が塊なら受け持ちの少ない方を塊とする。
  - `group_aligned` — 塊でなく、受け持ちが群と一致する（どの群でも多い方のスポットの受け持ちが
    `group_purity`（0.9）以上で、2 スポットとも少なくとも 1 つの群の多い方）。群ごとに別の
    異性体が出ている可能性があるので**まとめない**。印を付けてチームの監査とレポートに回す。
  - `interleaved` — ほか（試料ごとに入れ替わる）。
- **注釈の両立**（`names_compatible`）: 同じ名前（大文字小文字を無視）、GOSLIN の種のレベルで同じ
  和組成（`PC 16:0_18:1` と `PC 34:1`）、または片方が注釈なし。名前は SME の rank 1、無ければ SML。
  加えて、組の出力の `inchikey_match`（`feature_metadata[fid]["inchikey"]` どうしの照合）が `full`
  （InChIKey が全体で一致）か `first_block`（先頭の 14 文字の塊 = 骨格だけが一致。立体異性体の可能性）
  なら、名前が違っても両立とみなす（組の `names_compatible` はこの合算の真偽）。どちらかが無い・空・
  一致しないときは `none` で、名前だけで判定する。
- **受け持ちの数**: `n_owned` は `[a の数, b の数]` で、**gap-fill で値を補った試料の受け持ちを含む**
  （意味は変えない）。`n_owned_detected` は、そのうちどちらかのスポットで**検出された試料**
  だけで数えた同じ並びの数（検出状態が無い run では null）。`n_owned` が釣り合って見えても
  `n_owned_detected` が片側 0 なら、その側は gap-fill だけで支えられている。
- **`mergeable`** = `names_compatible` で、型が `injection_block`、または ΔRT ≤ `max_interleaved_rt_gap`
  （既定 0.15 分 = MS-DIAL のアライメントの RT 許容 0.1 分の 1.5 倍。アライメントの窓で割れたもので、
  別の異性体ではない。MTBLS6965 の陰イオンの FA 20:3 の異性体 401/402 は 0.20 分離れており、
  0.2 では別の異性体をまとめてしまう）の `interleaved`。ΔRT の大きい `interleaved` は載せるが
  まとめない。まとめると判定に使うデータが変わるので、迷う組はまとめない側に倒す。
- **まとめる組**（`merged`・`merge_groups`）: mergeable の組のうち**互いに最も近い相手**どうしだけ——
  A の mergeable な相手の中で |ΔRT| が最小なのが B で、かつ B のそれが A（最小が同点なら相手なし）。
  まとめは常に 2 スポットで、1 つのスポットは高々 1 組にしか入らない。ほかの mergeable な組は
  `merge_reason: "not_mutual_nearest"` で載せるがまとめない。例: MTBLS6965 の陰イオンの FA 20:3 は
  2 つの異性体が同じ後半の注入で一緒にずれ、401（後半の n-9）が 403（前半の n-9、ΔRT 0.088）とも
  402（前半の n-3/n-6、ΔRT 0.199）とも組になる——互いに最も近い 401↔403 だけをまとめる。
- **注入順**（`injection_order_source`）: `manifest`（pipeline の実験情報シートの `injection_order`）・
  `metadata`（シート以外から推した値、または呼び出し側が渡した順。渡した順に欠ける試料があれば
  エラー）・`filename`（1 試料でも
  欠ければ試料名の順。`block.injections` はその順位）。

`output_path` の JSON（キーは昇順。長いリストは要素ごとに改行する——`pairs` は 1 組 1 行で、行の範囲で読める。
共通の約束）:

```json
{"schema": "alignment-qc.v1",
 "dataset_fingerprint": "<解析の入力としてのデータセットの指紋。result_state.dataset_fingerprint>",
 "params": {"ppm": 15.0, "min_rt_gap": 0.05, "max_rt_gap": 1.0, "max_both_detected": 0.1,
            "max_interleaved_rt_gap": 0.15, "min_block_size": 2,
            "min_dominance": 3.0, "min_complementary": 0.9, "group_purity": 0.9,
            "signal": "detected"},
 "injection_order_source": "manifest",
 "pairs": [{"pair_id": "P1", "features": ["1265", "1266"], "mz": [..], "rt": [..],
            "delta_rt": 0.61, "delta_ppm": 7.7, "adduct": "[M+NH4]1+",
            "names": ["CE 18:2", null], "names_compatible": true, "type": "injection_block",
            "owner": {"<試料名>": "1265"}, "n_owned": [40, 2], "n_owned_detected": [40, 1], "inchikey_match": "none", "both_detected_fraction": 0.0,
            "block": {"spot": "1266", "injections": [1, 2]}, "mergeable": true,
            "merged": true, "merge_reason": null}],
 "groups": [{"group_id": "G1", "features": ["1265", "1266"], "pair_ids": ["P1"],
             "type": "injection_block", "types": ["injection_block"], "mergeable": true}],
 "merge_groups": [{"merge_group_id": "M1", "group_id": "G1", "features": ["1265", "1266"],
                   "pair_ids": ["P1"], "type": "injection_block", "types": ["injection_block"],
                   "mergeable": true}],
 "summary": {"n_pairs": 1, "by_type": {"injection_block": 1, "group_aligned": 0, "interleaved": 0},
             "n_features_affected": 2, "injection_block_median_delta_rt": 0.61,
             "n_mergeable_pairs": 1, "n_merged_pairs": 1, "n_groups": 1}}
```

- `dataset_fingerprint` は組を探したデータセットの指紋（`dataset_analysis_run` の `input_fingerprint.dataset_hash` と
  同じ関数。定量行列・特徴 ID・試料・検出状態・出所から作る。同じ run・同じハンドルなら同じ値）。
  `dataset_analysis_run` の `features.alignment_qc` はこれが今のデータセットと違う JSON（別の極性・別の run。
  特徴 ID は小さい整数でどちらにもあるので、黙ってまとめてしまう）、または無い JSON（この欄を足す前の版）を
  `ANALYSIS_INPUT_INVALID`（`field: features.alignment_qc`）で拒む。
- `pairs[].features` は `injection_block` なら [塊でない側, 塊の側]、ほかは RT の昇順。
  `delta_rt = rt[1] − rt[0]`（`injection_block` では塊の RT のずれ。符号つき）、
  `delta_ppm = (mz[1] − mz[0]) / mz[0] × 10⁶`。`block.injections` は塊の最初と最後の注入順。
- `groups[]` は**すべての組**の連結成分（報告用。3 スポット以上の連鎖を 1 つにする）。`features` は
  特徴 ID の小さい順（すべて整数なら数の順）。`type` は成分の型で、組の
  型が混ざれば `group_aligned` > `interleaved` > `injection_block` の順で強い方。`mergeable` は成分
  全体が 1 つのまとめる組のときだけ真。
- `merge_groups[]` はまとめる組ごとの 2 スポット（`merge_group_id`・`group_id`＝含まれる `groups[]` の
  成分・`features`・`pair_ids`・`type`・`types`・`mergeable: true`・`keep`・`keep_basis`・`n_owned`）。
  残すのは**受け持ちの多い側**、同数なら小さい ID（すべて整数なら数の順）で、`keep` がその ID、
  `keep_basis` が決め手（`n_owned` / `smaller_id_tie`）、`n_owned` が `{特徴 ID: 受け持ちの試料の数}`。
  **群の並びは残す ID が先頭**（`features`・`merged_from` の値。`groups[]` は ID の小さい順のまま）で、
  `merge_matrix_columns` は渡した並びで行列にある最初の ID を残す。`keep` の無い古い JSON の群は
  従来どおり小さい ID を残す。統計の側（`dataset_analysis_run` など）がまとめるのは
  `merge_groups` だけ。`summary.n_merged_pairs` はその数。`rule_hash` に入る群の hash
  （`merge_groups_hash`）は、残す ID が群の最小なら前の版と同じで、ほかの ID を残す群では変わる。
- `split_partners(qc)`（統計の側が使う）: 特徴 → `{group_id, merge_group_id, partners, component,
  type, mergeable}`。まとめる特徴は `partners` がその相手 1 つで `mergeable` が真、`type` はその組の型。
  まとめない特徴は `partners` が空・`mergeable` が偽・`merge_group_id` が null。`component` は全体の
  成分のほかの特徴（報告用）。
- 戻り値は `status`・`schema`・`injection_order_source`・`signal`（`detected` / `values`）・`summary`・
  `kept_on_handle: true`・`output_path`（渡したとき）と、`group_aligned` の組があれば `n_group_aligned`（全件）・
  `group_aligned_truncated`・先頭 50 組の TSV `group_aligned`（`pair_id feature_a feature_b mz rt_a rt_b delta_rt
  adduct name_a name_b`。超えれば `notes` に、全組は `output_path` の JSON の `pairs` にあると書く）。
- **最後の結果はハンドルに持つ**（`kept_on_handle`）。`output_path` と同じ中身（同じバイト列の sha256）を、
  そのハンドル（データセット）に最後の 1 つとして持ち、前の結果は置き換える。`dataset_feature_search` の
  `alignment_qc_path`・`dataset_novel_candidates` の `statistics.alignment_qc` を省いたときは、まずこれを使う
  （`merge_split.source: last_alignment_qc`。`alignment_qc_sha256` はファイルで渡したときと同じ値）。既定でない
  閾値（`ppm` ほか）で回した組を、パスを渡し直さずに使うため。**使わないのは**、データの指紋
  （`dataset_fingerprint`）が作ったときと違うとき、試料の表・除外が変わったとき（`metadata_revision`。
  `dataset_exclude`・`dataset_set_sample_metadata`）——群と注入順で型が変わるので、そのときは捨ててその場で探す。
  `dataset_analysis_run` は使わない（spec だけで再現する決まり。§2）。

## 10. 根拠バンドル v3（evidence-v3）

`dataset_export_evidence_bundle` が書く。v2 に試料表・試料ごとの値と検出・付加体・極性を足したもの。
`result_id` を省けば直近の差次的結果。対比は `contrast_id` か `result_id#contrast_id` で選び、2 群の対比が複数で選ばなければ
`needs_input`（選択肢は対比の一覧）。交互作用の対比は 2 群の形にならないので書けない（`CONTRAST_NOT_TWO_GROUP`）。
`dataset_export_evidence_bundle`（ARF 経路）は試料ごとの値を取り出さないので v2 のまま。契約の正は
massbank-context の `docs/contracts/evidence-bundle-v3.md`。TSV は全ファイルとも同じ書き方（タブ区切り・
UTF-8・改行は CRLF）。

**特徴 ID はスポットの ID のまま。** 割れの組（§9）へまとめた値は書かない。まとめた組への読み替えは
呼び出し側が `merged_from` で行う。

`manifest.json` の v2 との差分:

| 欄 | 意味 |
|---|---|
| `contract_version` / `exporter.version` | `3` / `evidence-v3` |
| `comparison.polarity` | データセットの極性（`positive` / `negative`）。ハンドルなら `load_dataset` の `polarity`、対話のセッションなら同じ規則で決める。決まらなければ `null` |
| `comparison.preprocess` | `value` に当たった前処理の設定（結果が読んだ行列の前処理の設定） |
| `comparison.method` | 実際に使った検定名（`moderated_t` / `welch` / `student` / `mannwhitney`） |
| `comparison.statistics` | `{test, ci_level, verdict_rule: {q_threshold, equivalence_fc, equivalence_log2}, bh_population: {name, n, q_interpretation}, subject_source, n_units_a, n_units_b}`。`n` は選んだ対比の多重比較の補正の母集団の大きさ（キーの名前は BH のときのまま）。`q_interpretation` は `correction` の方法での `q` の読み方（`discovery_fdr` / `fwer` / `unadjusted` / `reproduction_only`。`dataset.md` の書き出しのメタ行と同じ） |
| `comparison.effect_size_definition` / `mean_definition` | 効果量の定義（`{name, description}`）と、`mean_a` / `mean_b` が「個体ごとの前処理後の線形値の算術平均」であること |
| `comparison.matrix_id` / `rule_hash` / `contrast_id` | どの行列・どの規則・どの対比を書いたか（`contrast_id` は 2 群の結果なら 1 つの対比の ID） |
| `comparison.correction` | `q` を作った多重比較の補正。比較の記録の `p_adjust_method`（`bh` / `bonferroni` / `holm` / `none`）を `BH` / `Bonferroni` / `Holm` / `none` で書く。記録に無ければ `BH` |
| `exploratory` | 最上位の欄。探索専用のデータセット（未完了の実行を `allow_incomplete=True` で開いた・検証を通っていない pipeline の単位）から書いたときだけ `{exploratory_only: true, reason}`（書き出しは断らない）。探索専用でなければ欄が無い |
| `comparison.unadjusted_confounded` | `allow_confounded=true`（既定）で完全交絡（全対比をまとめた表か、その対比の表）を越えて走った比較だけ `true`（それ以外は欄が無い）。同じ印は来歴 `provenance.result.comparison.unadjusted_confounded`（と `provenance.result.warnings` の説明）にも載る |
| `intensities_stage` | `pre_imputation`（補完の前）/ `post_imputation`（補完の後） |
| `population.detection_unavailable_reason` | `detected` が全行空のときの理由。検出状態を使えれば `null` |
| `population.exclusions` | 除外の来歴。エンジンの結果は、その行列を作ったときの除外（`dataset_exclude` の試料・特徴。名前は 20 件で切り、残りは `n_more`）を、結果が読んだ行列の登録から引く（今の `ds.excluded_*` ではない） |
| `files[].role` | `samples` と `intensities` が増える |

`samples.tsv`（データセットの全試料）:

| 列 | 意味 |
|---|---|
| `sample_id` | 行列の試料名 |
| `group` | 比較の群 A の試料は `a`、群 B は `b`、比較に入らない試料（QC・ブランク・群なし・`include=false` で除いた試料）は空。`samples.tsv` には全試料が載るが、`intensities.tsv` には比較の試料だけ |
| `role` | 実験情報シートの `role`（除外された試料も `sample` のことがある）。シートが無ければ名前の `_` 区切りのトークンからの推定で、`QC1` のようにトークンにならない名前は `sample` になる（`group` が空でも QC とは限らない） |
| `injection_order` | 整数。分からなければ空 |

`intensities.tsv`（縦持ち。全特徴 × 比較に入る試料。QC・ブランクの行は無い）:

| 列 | 意味 |
|---|---|
| `feature_id` | `features.tsv` の行と同じ |
| `sample_id` | `samples.tsv` の試料名 |
| `value` | log 変換の前の値（有効数字 10 桁）。欠測は空。結果の試料が値の行列に無いときも空。実データでは `0` が多い（実測で全セルの 17〜28%、すべて `detected=false`）。補完は NaN だけを触るので `0` は実際の値で、欠測の空とは別。この場合 `pre_imputation` と `post_imputation` の値は同じになりうる。検出の有無は `value` でなく `detected` で読む |
| `detected` | `true` / `false` / 空。`false` は gap-fill か補完で、実測のピークではない。**空は「不明」で、gap-fill の意味ではない** |

`intensities_stage` の決め方:

- エンジンの結果（`dataset_compare`）: 常に `pre_imputation`。値は結果が読んだ行列（`ds.stat_matrices[matrix_id]`）の `sample_values`——**正規化の後・log2 の前の線形値**で、欠測は空（エンジンは検定に補完しない）。前処理のレシピを当て直さない。行列が登録に無い（除外・試料の表の変更で古くなった）ときは `missing_state`（`dataset_preprocess` と `dataset_compare` をやり直す）、結果の行と行列の特徴の集合が違えば `EVIDENCE_INTENSITIES_MISMATCH`。

`detection_unavailable_reason` の値:

| 値 | 意味 |
|---|---|
| `dataset_status` の `detection` と同じ理由（`sample_axis_mismatch`・`mz_mismatch` など） | 検出状態を取り込めなかった。取り込みのときの理由 |
| `not_imported` | 検出状態が無く、理由も記録されていない |
| `mask_shape_mismatch` | 検出状態の形が特徴 × 試料と合わない |
| `sample_axis_mismatch` | 結果の試料が `ds.sample_names` に無い（検出状態の試料軸とずれている） |

`features.tsv` に任意の列 `adduct` が増える。mzTab-M の SML の `adduct_ions`（SML が無い特徴は rank 1 の
SME の付加体）。無ければ空。`features.tsv` の `significant` は差次的の書き出しの `significant` 列と同じ印で、
エンジンの結果は判定 `verdict` が `changed` の行だけ `true` になりうる（`dataset.md` §11.4。判定できない行は q が
小さくても `false`）。

戻り値には v2 の欄に `n_samples`・`n_intensity_rows`・`intensities_stage` が足される（探索専用のデータセットなら
`exploratory` も）。

### `spectrum_links.tsv` と `.dcl` の出どころ

dataset 経路も ARF 経路（`dataset_export_evidence_bundle`）と同じ規則で、**同じアラインメントの `.dcl`** の
代表スペクトルに結ぶ。結ぶのは ID だけ（**SMF_ID == `.arf2` の AlignmentID == `.dcl` の番号**）で、
m/z・RT は結びの検証にだけ使う（m/z の近さでは結ばない）。結べた特徴は 1 行（`spectrum_file` は `.dcl`
のファイル名、`spectrum_id` は `<.dcl の語幹>:<番号>`、`mapping_evidence` は
`same-alignment-stem; SMF_ID == AlignmentID == DCL index; precursor/RT verified`）。結べた行が 1 件でも
あれば `.dcl` をバンドルへ写し、`files[]` に `role: spectra` で載る。

使う `.arf2` / `.dcl`: ハンドルなら `load_dataset` が決めたもの（§1）。対話のセッションなら
同じ探し方（検出状態を取り込んだ `.arf` の同じバッチ → mzTab-M の隣）。

検証は 2 段:

1. **番号の空間**: SMF_ID の m/z・RT（mzTab-M）と、同じ番号の `.arf2` のスポットの m/z・RT を比べる。
   1 件でも合わなければ全体を結ばない（番号の付け方が違う別のアラインメントとみなす）。
2. **特徴ごと**: `.arf2` のスポットと `.dcl` の同じ番号の precursor・RT を比べる（ARF 経路と同じ）。

`manifest.json` の `provenance` の欄:

| 欄 | 意味 |
|---|---|
| `spectrum_links_status` | `verified_alignment`（結んだ）/ `not_provided`（`.arf2` か `.dcl` が見つからない）/ `id_space_mismatch`（1 の検証で合わない特徴があった）/ `id_space_unverified`（m/z を比べられる特徴が無く、1 を確かめられない）。`verified_alignment` 以外は links が空で、**空は「MS/MS が無い」の意味ではない** |
| `spectrum_links_reason` | `verified_alignment` 以外のときの理由（文）。結べたときは `null` |
| `alignment_files` | 使った `{arf2, dcl}` のファイル名。見つからなければ `null` |
| `id_space_check` | 1 の検証の `{n_checked, n_mismatch, examples}`（`examples` は合わない特徴の先頭 20 件。SMF と `.arf2` の m/z・RT）。`.arf2` / `.dcl` が無ければ `null` |
| `msms_link_statuses` | 特徴ごとの `{status, dcl_index}`（ARF 経路と同じ形）。`status` は `ok`（結んだ）/ `no_msms`（エントリはあるがピークが無い＝MS/MS 未取得）/ `precursor_mismatch` / `metadata_missing`（`.arf2` にそのスポットが無いか m/z・RT が無い）/ `out_of_range`（`.dcl` にその番号が無い）/ `spot_id_not_integer`。`verified_alignment` 以外では `{}` |
| `msms_link_status_counts` | `msms_link_statuses` の `status` ごとの件数（0 件の状態は載らない） |
