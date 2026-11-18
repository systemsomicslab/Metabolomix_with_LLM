# 再解析の道具（ハンドル・spec・証拠の束）の出力定義

再解析で使うツールの戻り値と、spec・結果の JSON の欄の意味。既存のツールに足した引数（`dataset_handle` ほか）は、
渡さない限り戻り値を変えない。呼び出し連鎖は `docs/workflow/reanalysis.md`。

共通の約束:

- **spec は正規化して返す。** 既定値を埋め、未知のキーは拒否する。結果には正規化した spec と `rule_hash` が入る。
  `rule_hash` は正規化した spec からパスの欄（`path` / `matrix_path`）を除いて作る。入力のファイル（行列・群・共変量・
  注入順の TSV）は `input_fingerprint` に sha256 で入る。
- **同じ spec・同じ入力なら同じ JSON。** 時刻や乱数の ID は結果に入れない（例外は `dataset_analysis_run` の戻り値の
  `record` だけ。`output_path` の JSON には入らない）。
- **大きな結果は `output_path` の JSON に書く。** キーは昇順で字下げなし。1 行に収まらない辞書・リストは要素ごとに改行し、
  1 行は辞書のキーと区切りの `,` を含めて 1,000 字以下（数や文字列だけのリストは 1,000 字までずつ詰める。分けられない
  1 つの値の行だけは超える）。`output_path` を省いて上限を超えても断らず、先頭の行だけを載せて `truncated: true` と
  全件の数を返す（上限と並びは各節）。
- **効かない設定は誤りにしない。** 書いても使われない設定（recipe の `steps` に無い段の設定、`merge_split=false` の
  `alignment_qc`、ols でない検定の共変量、`rt_min` の無い主張の `rt_offset_min` など）は捨てて `notes` に理由を書き、
  正規化した spec にも `rule_hash` にも入れない。例外: analysis-spec.v1 の次の設定は `notes` に書くが spec と `rule_hash` に
  残す——`correction="none"` の解析の `bh_population` と `thresholds.q`（q が null なので `meets_thresholds` はどの行も
  false）、`drift_correct` を `steps` に持つ解析の無い spec の `injection_order`、`test.kind=ols` で `test.covariates` の
  ある解析の無い spec の `covariates`（この 2 つの TSV は読まないので、ファイルが無い・試料が欠けても断らず、
  `input_fingerprint.files` にも載らない）。`dataset_novel_candidates` の `statistics.analysis_spec` で選んだ解析が
  `correction="none"` なら、q が無いので候補が出ないことも `notes` に書く。形の誤りと、意味がぶつかる組み合わせ
  （`aggregate` と `features.feature_ids`）は誤り。
- **判定（再現したか・明らかに別物か）は返さない。** 返すのは数字と、spec の規則を機械的に当てた結果だけ。
- **研究室の参照ライブラリのパスは出さない。** ファイル名と sha256 だけ。
- **記録は解析フォルダの `reports/` に残る。** `dataset_analysis_run` は呼ぶたびに `reports/stats/` へ 1 件（§2 の「記録」）。
  `write_report`・`record_objective`・`save_figure` が上書きするファイルの前の版は `reports/.history/<名前>.<旧版の mtime の
  UTC。例 20261118T143012Z>[-n]<拡張子>` へ移す（mtime は旧版のまま。中身が同じなら書かない・戻す）。`read_report`・その索引・
  目的の読み戻しは `.history/` を見ない。

## 1. データセットのハンドル

ハンドル（`dh_` + 12 桁）は、1 つのサーバのプロセスにデータセット（とその極性のライブラリ）を複数並べて持つための ID。
`load_dataset` が返し、`dataset_*`・`library_*`・再解析のツールの `dataset_handle` に渡す。

- **対話のセッションは直近のハンドルを指すだけ。** `load_dataset` が開いたハンドルはいつも直近になる（戻り値の
  `dataset_handle`、`dataset_list` の `current`）。ハンドルを渡さない呼び出しは直近のハンドル（のデータセット・
  `library_load` のライブラリ・直近の照合）を使う。ハンドルを渡した呼び出しはそのハンドルだけを読み書きし、直近を
  読みも切り替えもしない。例外:
  - 定量の派生（`dataset_quant_apply`・`dataset_merge_polarities`・`dataset_select_features`）は、`dataset_handle` を
    渡しても作った派生の新しいハンドルを直近にする（上限で止まれば直近は動かない）。以後ハンドルを渡さない呼び出しは
    派生に効く（直近は `dataset_list` の `current` で確かめる）。
  - `pathway_project` は、ハンドルを渡しても対話のセッションの「直近のパスウェイの run」（`pathway_plot` が `run_id` を
    省いたときに描くもの）を書き換える。

  直近のハンドルを閉じると、渡さない呼び出しは `missing_state`（開いている別のハンドルにも、前に読んだライブラリにも
  落ちない）。順番を待つ間に直近のハンドルが閉じられた呼び出しも、何もせず `missing_state`。`load_dataset` より前に
  `library_load` で読んだライブラリは、このデータセットのライブラリが決まらなかった（記憶が無い・`library="auto"`）
  とき、そのハンドルへ移る（写しではない）。
- **ハンドルどうしも別の状態。** 同じ入力を開き直すと読み込み済みのもの（`reused_base`）を使い回すが、ハンドルには
  複製が渡り、片方の前処理はもう片方に見えない。読み込みが読んだ・探したファイル（入力・`.dcl` などの付随ファイル・
  ライブラリ・ジョブの記録と実効メソッド・後から置かれうる `.arf` の候補）の size か更新時刻が変われば読み直す。
- **並行と直列。** ハンドルを渡した呼び出しはサーバの別スレッドで回り、別のハンドルどうしは並行。同じハンドルへの
  呼び出しは 1 本ずつ（到着順は保証しない）。ハンドルを渡さない呼び出しも、直近のハンドルがあればそのハンドルの
  呼び出しとして並ぶ（途中で `load_dataset` が直近を切り替えても、始めたときのハンドルで最後まで動く）。
  - 順番は最長 25 秒待つ。それでもハンドルが空かなければ**何もせず** `HANDLE_BUSY`（`details.dataset_handle`・`waited_s`。
    持っているのが長い道具の仕事なら `details.job_id` / `tool` / `running_s` / `next` と `required_tools: ["job_status"]`）。
    `job_status(job_id)` で仕事の終わりを待ってから、同じ呼び出しをし直す。`HANDLE_BUSY` を返した呼び出しが後で効く
    ことはない。
  - 上限があるのはロックの待ちだけで、道具の走る時間には上限が無い。25 秒の内にハンドルを取れた呼び出しは、
    クライアントの中継が 36〜60 秒ほどで打ち切っても、サーバでは最後まで動いて効く（行列の登録や結果の記録は残り、
    `matrix_id` を省いた後の `dataset_compare` / `dataset_pca` はその行列を使う）。状態を変える呼び出しが打ち切られたら、
    し直す前に `dataset_status(view="results")` などで効いたかを確かめる。
  - `dataset_close` も、仕事が走っているハンドルは閉じずに `HANDLE_BUSY`（仕事は止めない）。2 つのハンドルを取る
    `dataset_merge_polarities` も同じ封筒で、待つのは 2 本合わせて最長 25 秒（`details.dataset_handle` は空かなかったほう）。
- **上限。** `lipidmix.local.toml` の `[server] max_dataset_handles`（既定 8。1 未満は 1）。`load_dataset` がハンドルを
  1 本使う（前のものは閉じずに残る）。超えると `DATASET_HANDLE_LIMIT` で止め、古いハンドルを閉じない。寿命は
  `dataset_close` か、サーバの終了まで。
- **受付番号。** `load_dataset` は最長 `wait_seconds`（既定 20。0〜30 秒。範囲の外は丸めて `wait_seconds_note`）待ち、
  終わらなければ読み込みを続けたまま `{"status": "loading", "load_id": "ld_…"}` を返す。`dataset_load_status(load_id)` は
  読み込み中なら `status: "loading"`、終われば `load_dataset` の戻り値そのもの（1 時間覚える）。`load_id` を省くと、
  受付番号を渡し済みで結果を渡していない読み込み（自分の `wait_seconds` の中で待っている呼び出しのものは除く）が
  1 本ならそれ、2 本以上なら `needs_input`（受付の一覧つき）、無ければ `missing_state`。`wait_seconds`（既定 0）を
  渡すと最長 30 秒まで終わりを待つ。受付番号を返した後に別のハンドルが直近になっていたら直近にせず、`current: false` を
  付ける。上限は読み込み中の受付も 1 本と数える。
- **長い道具の受付番号。** `pathway_project` と `dataset_homolog_check` も `wait_seconds`（既定 20。0〜30 秒）を持ち、
  その中で終わらなければ仕事を続けたまま `{"status": "running", "job_id": "job_…", "tool", "dataset_handle", "next"}` を
  返す（`dataset_handle` はその仕事のハンドルで、ハンドルなしの仕事には付かない。同じハンドルで先の長い仕事を待って
  いれば `queued_behind: {job_id, tool}`）。ロックは仕事が終わるまで仕事が持つ（ほかの呼び出しは「並行と直列」のとおり
  `HANDLE_BUSY`。後から呼んだ長い道具の仕事だけは、受付番号を返してからロックを待つ）。
  - 結果は `job_status(job_id, wait_seconds)`: 走っている間は `status: "running"`・`dataset_handle`・`elapsed_s`、終われば
    道具の戻り値（JSON のオブジェクトなら先頭に `job_id` と `dataset_handle` を足す）。
  - `job_id` を省くと、受付番号を渡し済みで結果を渡していない仕事が 1 本で、それが直近のハンドルかハンドルなしの
    仕事ならそれ。直近でないハンドルの仕事 1 本や 2 本以上なら `needs_input`（選択肢の label にハンドル）、無ければ
    `missing_state`。
  - 知らない番号は `JOB_NOT_FOUND`（終わって 1 時間を過ぎた・件数の上限で古いものから捨てた・サーバが再起動した。
    道具を呼び直す）。走っている受付が上限なら、長い道具が `JOB_LIMIT`（`details.running` と `max_entries`。
    `job_status` でどれかが終わるのを待って呼び直す）。道具が `DomainError` 以外で止まれば、その結果（道具の戻り値か
    `job_status`）が `JOB_FAILED`（`details.exception` は例外の種類だけ。本文はサーバのログ）。
- **ハンドルを見失う条件。** 閉じた・無いハンドルを渡した呼び出しは `DATASET_HANDLE_NOT_FOUND`（`details.open_handles` に
  開いているハンドル）で、直近や別のハンドルへは落ちない。ハンドルが無くなるのは:
  (1) `dataset_close` で閉じた（戻り値の `source` が閉じたものの入口の種類・入力の名前・極性）、
  (2) サーバのプロセスが終わった・再起動した（`server_update` の後、クライアントの再接続で起動し直したときも）、
  (3) 読み込みが受付の後に失敗した（ハンドルは登録されない）、
  (4) 別のサーバのプロセスで開いた（既定の stdio では MCP クライアント（Claude Desktop・Claude Code の各セッション・
  WebUI など）ごとにプロセスが 1 本起き、登録簿はプロセスごとに別。ハンドルを受け渡せるのは同じプロセスにつないだ
  者どうしだけ。`LIPIDMIX_TRANSPORT=streamable-http` の HTTP のサーバは、すべての接続が 1 つの登録簿と 1 つの直近を
  共有し、別の接続の `dataset_close`・直近の切り替えも効く）。

  ハンドルを渡さない呼び出しの行き先は、`load_dataset`（受付番号を返した後に終わった読み込みも、その間に別のハンドルが
  直近になっていなければ、呼んだ側が待つのをやめていても終わったときに直近になる。直近のハンドルを閉じても世代は
  進まないので、閉じた後に終わった読み込みは直近になる）と定量の派生で変わる。並行のエージェントはいつもハンドルを渡す。
  受付番号（`load_id`・`job_id`）は終わって 1 時間・件数の上限・再起動で忘れられ、その後は自分の読み込みのハンドルを
  `dataset_list` の `source`・`polarity` で見分ける。見失ったら `dataset_list` で確かめ、`load_dataset` で開き直す（同じ
  入力なら読み込み済みの基を使い回すが、前処理・比較の結果はハンドルごとなので作り直す）。設定ファイル
  （`lipidmix.local.toml`）が読めないとき、`dataset_list` は `CONFIG_INVALID` の封筒だけを返す（開いているハンドルは
  残る。設定を直して引き直す。ハンドルの ID は設定を読まない道具にはそのまま渡せる）。

`load_dataset` の戻り値:

| 欄 | 意味 |
|---|---|
| `dataset_handle` | ハンドルの ID |
| `assay_kind` | このハンドルの種別（`lipid` / `metabolite` / `unknown`）。ハンドルごとに持ち、`record_objective(assay_kind=...)` で確定・変更する。解析フォルダ配下 `reports/` の目的ファイルに記録があれば読み戻す（`notes` に注記）。データセットを開く前に（`[paths] reports_dir` の目的へ）記録した種別は、次に直近になるハンドルへ引き継ぎ、その目的にこのハンドルの解析フォルダ（`analysis_dir`）を書く（書けなかった・どの目的の種別か分からないときは `notes` が記録し直すよう言う）。目的ファイルの種別と食い違えば `unknown`（`notes` に注記）。共有の `[paths] reports_dir` の目的は、frontmatter の `analysis_dir` がこの解析フォルダのものだけを読む（`record_objective` は別の解析フォルダに結び付いた目的を書き換えず、自分の `reports/` に、そこに書けなければ共有の置き場所に `<analysis_id>-2.objective.md` … の名前で作る。開く前の記録もまだ結び付いていない別の名前の目的に書く） |
| `adapter` | 読んだアダプタ。`msdial_folder`（解析フォルダの `.arf` / `.arf2` を直接）/ `mztab`（mzTab-M）。`job_path` / `pipeline_path` の既定は `msdial_folder`（`source="mztab"` で `mztab`） |
| `polarity` / `polarity_source` | データセットの極性と出所（`job_declared` / `adduct_majority` / `arf2_ion_mode` / `unknown`）。`dataset_feature_search` の目標照合と同じ決め方 |
| `n_features` / `n_samples` | 特徴数・試料数 |
| `detection_available` | 検出状態（gap-fill の区別）を取り込めたか |
| `arf` / `not_chosen` | 解析フォルダで読んだとき、選んだ `.arf` と選ばなかった `.arf` のファイル名。mzTab-M で読んだときは `null` |
| `library` | `{file, sha256, record_count}`（置き場所は出さない）。読まなかったときは `null` で、理由が `library_note` |
| `library_record` | データセットが同定に使ったライブラリの記録（[library.md](library.md) の LibraryRecord。置き場所は出さない） |
| `arf2` / `dcl` | スポットの照合（spot_id = SMF_ID）に使うファイル名。解析フォルダで読んだならその `.arf2` / `.dcl`。決まらなければ `null` |
| `reused_base` | 読み込み済みのものを使い回したか |
| `summary` | 特徴・試料の数、定量種別、検出状態、同定の内訳、出所の裏取り、読み込みの注記（行の配列） |
| `sample_counts` | `{groups, roles, n_ungrouped}`: 群（Class ID）ごと・役割ごとの試料数と、群の宣言が無い試料の数。`groups` と `n_ungrouped` は統計に入る試料だけを数え、シートの `include=false` の試料は `excluded_by_sheet`（シートの群 → 件数。無ければ鍵ごと無い）に分ける（`dataset_status` の要約と同じ。`excluded_by_dataset_exclude` は `dataset_status` にだけ出る） |
| `overview` | `overview=True` のときだけ。最初の概観（[dataset.md](dataset.md) §11.9） |
| `wait_seconds_note` | `wait_seconds` を 0〜30 秒に丸めたときだけ |
| `next` | 次の一手の例（3 行。最初が `dataset_overview()`） |

ライブラリは、データセットが同定に使ったものの記憶（`library_record`。[library.md](library.md)）から引く。実体の sha256 を
測り直して `known` のままなら読み取り専用で開き、`mismatch`（内容が変わった・ファイルが無い）・`unknown`（決められない）・
`lbm_only`（LBM だけ）は開かずに `library_note` に理由を書く。**極性の設定（`[library] msp_positive` / `msp_negative`）には
落ちない。** ハンドルへの `library_load`（`file_path` を省く）も記憶が `known` ならそれを開き、`mismatch` / `lbm_only` は
`LIBRARY_RECORD_MISMATCH` / `LIBRARY_RECORD_LBM_ONLY`、`unknown` は `LIBRARY_NOT_REMEMBERED`。`file_path` を明示すれば記憶と
違っても開き、戻り値の `differs_from_record` に記憶したファイル名と開いたファイル名を書く。

**ハンドルは対話のデータフォルダを探さない。** スポットの照合（`library_match_feature` の spot_id / precursor_mz・
`dcl_spot_spectrum`・`dcl_diagnostic_ions`・`library_msms_qualify`・証拠の束の MS/MS）は、ハンドルの `.arf2` / `.dcl`
（か明示した `arf2_path` / `dcl_file`）だけを使い、無ければ `status: "no_alignment_files"`。

## 2. 統計の spec（analysis-spec.v1 / v2）と結果（analysis-result.v1 / v2）

**この節の前半（表と欄）は v1。** v2（統計エンジン）は §2.1。`dataset_analysis_run` は `spec.schema` で振り分ける。

`dataset_analysis_run(spec, dataset_handle=None, output_path=None)`。名前つきの解析（`analyses`）を並べる（感度解析は
同じ spec の別の解析）。

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
| `input` | `{"source": "dataset"}` | データセット（ハンドルか対話のセッション）か、外部の行列。TSV は `{"source": "tsv", "matrix_path", "orientation": "features_by_samples" \| "samples_by_features", "id_column"}`（欠測は空・`NA`・`NaN`・`null`、数として読めない値は誤り）。MetaboLights の MAF `{"source": "maf", "path"}`・Metabolomics Workbench の mwTab `{"source": "mwtab", "path"}`・補足の xlsx `{"source": "xlsx", "path", "sheet": null, "feature_column": "A", "first_sample_column": "B"}`（下の「外部の表」） |
| `groups` | データセットなら `metadata` | `metadata` は実験情報シート（`group`・`role`・`include`）。TSV は `{"source": "tsv", "path", "sample_column": "sample", "group_column": "group", "role_column": null}`。外部の行列には TSV が要る |
| `comparison` | 必須 | `reference`（基準）と `test`（比較）の群の名前。比較名は `<test>_vs_<reference>` |
| `injection_order` | `null` | ドリフト補正の注入順。`{"source": "metadata"}` か TSV `{"source": "tsv", "path", "sample_column", "order_column": "order"}`。`metadata` はデータセットの入力でだけ使える（外部の行列で `drift_correct` の解析があれば `ANALYSIS_SPEC_INVALID`。無ければ `notes` に書いて進む） |
| `covariates` | `null` | OLS の共変量の TSV `{"path", "sample_column": "sample", "missing": "drop" \| "error", "columns": [{"name", "type": "numeric"} \| {"name", "type": "categorical", "reference"}]}`。カテゴリは基準の水準に対する 0/1 の列に符号化する |
| `features` | `{"all": true, "merge_split": true}` | 報告する特徴（`{"feature_ids": [...]}` で絞る）と割れたスポットのまとめ（下の「割れたスポットのまとめ」）。BH の母集団は `bh_population` が決める。`merge_split` は無ければ true。`alignment_qc` は `dataset_alignment_qc` の `output_path` の JSON のパス（任意） |
| `aggregate` | `null` | 1 本に集約して検定する。`{"kind": "sum", "label", "feature_ids": [...]}`・`{"kind": "sum", "label", "select": {"name_regex", "evidence": ["sme"]}}`（クラスの主張）・`{"kind": "ratio", "label", "numerator", "denominator"}`（比の主張）。和はメンバーごとの寄与（`contributions`）も返す。**`features.feature_ids` とは併用できない**（`ANALYSIS_SPEC_INVALID`・`field: features`。種ごとの行は `aggregate` の無い別の spec で回す） |
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
| `bh_population` | `matrix` | `matrix`（その解析の行列で検定できた全特徴）・`reported`（報告する特徴だけ）・`annotated`（SME のある特徴 = `ds.feature_candidates` に候補のある特徴だけ）。検定できない特徴は母集団に入れない。`annotated` の母集団の外の特徴は報告しても `q` は null。まとめた特徴は元のどちらかに SME があれば母集団に入る。外部の行列では `correction=bh` の解析に `annotated` は使えない（`ANALYSIS_SPEC_INVALID`。`correction=none` なら `notes` に書いて進む） |
| `thresholds` | すべて `null` | `p`・`q` は「未満」、`fc` は `|log2fc| ≥ log2(fc)`。null の閾値は見ない |

**割れたスポットのまとめ（`features.merge_split`、既定 true）。** 真なら recipe の前（定量行列と検出状態を読んだ直後）に、
§9 の `merge_groups`（互いに最も近い mergeable な組。常に 2 スポット）を試料ごとの大きい方の値（NaN は無視）の 1 列に
まとめる。列の ID は受け持ちの多い側（同数なら小さい ID。§9）、検出状態は「どちらかで検出」。recipe はまとめた行列に
当たる（`bh_population=matrix` の母集団もまとめた後の数）。

- 組は `features.alignment_qc` に渡した `dataset_alignment_qc` の JSON から、無ければその場で `find_split_pairs`（既定の値）を
  回して取る。JSON の `dataset_fingerprint` が今のデータセットと合わなければ `ANALYSIS_INPUT_INVALID`（§9）。
  `dataset_alignment_qc` がハンドルに持つ最後の結果（§9）は使わない（spec だけで再現する。既定でない閾値の組を使うなら
  その JSON を書く）。`group_aligned` など mergeable でない組はまとめない。
- 報告する特徴（`features.feature_ids`）に消えた側の ID を書けば、残した 1 行（`merged_from` 付き）に当たる（両方書いても
  1 行）。`aggregate` の `feature_ids`・`select`・比の分子と分母も同じ置き換えの後で和・比をとる（両方選ばれても 1 回だけ
  数える）。`aggregate.members` は選んだままの並び、`summed` はまとめた後の並び、`aggregate.merged_from` は選んだ特徴の
  うちまとめたもの。比の分子と分母が同じ物質にまとまるなら `ANALYSIS_SPEC_INVALID`。
- 行列に 2 つ揃わない組（別の run の JSON など）はまとめず、結果の `merge_split.skipped` に挙げる。
- 外部の行列（`input.source` が `dataset` 以外）では何もしない（`features.alignment_qc` は捨てて `notes`）。
  `merge_split=false` なら 2 スポットのまま検定する（`alignment_qc` は捨てて `notes`）。

**外部の表（MAF・mwTab・xlsx）。** 論文側の値の表を特徴 × 試料の行列として読む。数として読めないセル（`ND`・`<LOQ`
など）と非有限の値は欠測にし、空でない読めないセルの数と例（`input_read.n_unreadable_cells`・`unreadable_examples`。
例は最大 5 種）を返す（空・`NA`・`NaN`・`null`・`none` は欠測で、数えない）。試料・特徴の重複は `ANALYSIS_INPUT_INVALID`、
ファイルが無ければ `ANALYSIS_INPUT_NOT_FOUND`。群は TSV で渡す（`groups.source=metadata` は使えない）。

- **MAF**（`m_*.tsv`）: 見出しの `smallmolecule_abundance_std_error_sub` より右の列が試料。特徴 ID は
  `row_<n>:<metabolite_identification>`（`n` は見出しの次の行を 1 とするデータの行の番号。空の行も数えるので、同じ名前の
  行も別の特徴）。見出しにこの列が無ければ `ANALYSIS_INPUT_INVALID`。
- **mwTab**（テキスト形式）: `MS_METABOLITE_DATA_START` 〜 `MS_METABOLITE_DATA_END` の区画。1 行目が `Samples`、`Factors` の
  行は読まない。特徴 ID は 1 列目の名前。
- **xlsx**: 1 枚のシート（`sheet` を省けば先頭。無い名前はあるシートの名前を添えて `ANALYSIS_INPUT_INVALID`）の横長の表。
  列は Excel の文字（`A`・`B`・`AA`。`first_sample_column` は `feature_column` より右）。見出しの行は 2 つの列の両方に値が
  ある最初の行で、それより上は読まない。試料は見出しの行の `first_sample_column` から右の空でないセル。試料のセルに数が
  1 つも無く、欠測の印でない文字がある行（群の見出しなど）は特徴にせず、その数と名前の例（最大 5）を
  `input_read.n_rows_skipped_as_text`・`rows_skipped_as_text` に返す（xlsx だけの欄）。欠測の印（`ND`・`n.d.`・`NQ`・`BDL`・
  `LOD`・`LOQ`・`tr`・`-`・`n/a`・`not detected`、`<`・`>`・`≤`・`≥` で始まるもの、`#` で始まる Excel のエラー値）だけの行は、
  どの試料でも検出されなかった物質として特徴に読み、そのセルを読めないセルに数える。spec の特徴 ID
  （`features.feature_ids`・`aggregate`・`correlate.pairs`）が飛ばした行の名前なら、`FEATURE_NOT_FOUND` がそう書き添える
  （`details.skipped_as_text`）。式は計算せず、保存された値を読む。

ファイルの sha256 は `input_fingerprint.files.matrix`。

**外部の表との相関（`correlate`、任意）。** 論文の表の値とこちらの値の一致を、特徴の対応ごとの Spearman の ρ で見る
（判定はしない）。

```json
"correlate": {"external": {"source": "mwtab", "path": "ST000000.txt"},
              "pairs": [{"feature_id": "1265", "external_id": "Glucose"}],
              "sample_map": {"path": "map.tsv", "sample_column": "sample", "external_column": "external_sample"},
              "stage": "input"}
```

- `external` は外部の表（`input` と同じ形: TSV・MAF・mwTab・xlsx）。`pairs` はこちらの `feature_id` と外部の表の特徴 ID
  （`external_id`）の組で、**対応は spec が渡す**（名前の照合はしない）。どちらかの側に無い ID は `FEATURE_NOT_FOUND`
  （割れの消える側の ID は残す ID に当たる。行の `merged_into`）。
- `sample_map`（任意）は試料名の対応の TSV（こちらの試料 → 外部の表の試料）。無ければ同じ名前で照合する。当たった試料
  だけを使い、1 つも当たらなければ `SAMPLE_MISMATCH`。対応は 1 対 1 で、こちらの試料の重複も、外部の同じ試料を複数の
  試料が指す対応も `SAMPLE_MISMATCH`（`details.duplicated` / `duplicated_external`。反復の注入は 1 つだけ書く）。
- `stage`: `input`（既定。読んだ行列＝割れのまとめの後・recipe の前。群・役割を問わず当たった全試料）か `matrix`（解析
  ごとの recipe の後・集約の前の行列。行に `analysis_id`。recipe が落とした特徴は `reason: not_in_matrix`）。
- ρ は両方が有限の試料だけで、同順位は平均の順位。p は両側で、出し方は行の `p_method`:
  - 9 組以下、または同順位なしで 14 組以下: `exact_permutation`（全順列を数えた正確な p。同順位が無ければ S = Σd² の正確な
    分布で、9 組以下は R の `cor.test(method = "spearman")` と同じ値。同順位があれば平均の順位のまま並べ替える条件付きの
    正確な p）。
  - 15〜1290 組で同順位なし: `edgeworth`（AS 89 の Edgeworth 展開。R と同じ）。
  - 1290 組を超える: `t_approx`（t 近似、自由度 n−2）。
  - 10〜1290 組で同順位あり（未検出の 0 など）、または近似（`edgeworth` / `t_approx`）の両側 p が 1e-3 未満:
    `monte_carlo_permutation`（100,000 回の並べ替え。乱数の種は固定で、同じ入力なら同じ p）。
  - **並べ替えの p の最小は 2/(回数+1) ≈ 2.0e-5**。この値は観測より極端な並べ替えが 1 つも無かった意味で、本当の p は
    それ以下（上限として読む。論文の p がこれより小さくても食い違いではない）。1e-4 を下回る辺りは桁の目安。

  p は 0 にならない（正確な p の最小は 2/n!、近似は 1e-3 以上だけを使う）。組が 5 未満なら ρ と p は null
  （`reason: too_few_pairs`）、どちらかが一定なら null（`constant`）。どちらも `p_method` は null。
- 結果の `correlate` は `{stage, external（読み取りの要約）, n_samples_matched, unmatched_samples, unmatched_external_samples,
  rows: [{analysis_id, feature_id, external_id, merged_into, n, rho, p, reason, p_method}]}`（戻り値の TSV の列は
  `analysis_id feature_id external_id merged_into n rho p p_method reason`。ρ は小数 4 桁、p は有効 4 桁で 0 にしない。
  JSON は丸めない）。外部の表と対応表の sha256 は `input_fingerprint.files` の `correlate.external`・`correlate.sample_map`。
  書いたときだけ正規化した spec と `rule_hash` に入る。

**試料の照合は完全一致。** 群の TSV は行列の試料と過不足なく一致しなければ `SAMPLE_MISMATCH`（欠け・余り・重複を名前で
挙げる）。注入順・共変量の TSV は行列の試料をすべて含めばよい（余りは無視）。群に入るのは役割が `sample` で `include` が
真の試料だけで、群の n が 2 未満なら `STATISTIC_GROUP_TOO_SMALL`。データセットの試料名は実験情報シートの `sample_id`
（無ければ assay 名）。

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
| `merged_from` | 割れた 2 スポットをまとめた特徴だけ、まとめた特徴 ID（先頭が `feature_id` = 残した ID）。ほかは null（TSV では `;` 区切り） |
| `detected_reference` / `detected_test` | 検定に入れた群の試料のうち、検出された（gap-fill でない）数。まとめた特徴は「どちらかのスポットで検出」の試料を数える。外部の行列・検出状態の無いデータセット・`aggregate` の行・`not_in_matrix` の行は null。群の試料数は `n_ref` / `n_test` ではなく群の `members`（値が欠けた試料も入る） |
| `detection_fisher_p` | 検出の有無の 2×2 に対する Fisher の正確検定の両側 p（上の 2 欄が null なら null） |
| `constant_in_group` | 一方の群の**読んだ値**（割れのまとめの後・recipe の前。recipe に補完の段があれば補完の直後＝後の正規化・ドリフト補正の前の値）で、観測の値が 2 つ以上あり、すべて同じ（相対の丸めの誤差の内）なら、その群の名前（TSV では `;` 区切り。両群なら 2 つ）。ほかは null。観測の値は `transform="log2"` なら有限の正値（0 だけの群は印を付けない。その行は検定できない）、`none` なら有限の値（0 も値）。recipe の後の値では探さない。統計エンジンの `constant_in_group` と同じ規則（`spec_statistics.constant_columns`）。埋め値（`half_min`・LOD など）の印で、検定の値も判定も変えない。行の最後の欄 |
| `detection_driven` | 差が値の大小より検出の有無から来ているか。高い側の群（`direction`）の検出が群の試料数の半分未満か、低い側の群の検出が 0 なら true。`status` が `tested` で `direction` のある行だけ（ほかは null。一方の群の平均がちょうど 0 で変換も補完も無い行は向きが決まらず null）。判定（`meets_thresholds` など）は変えない |

解析ごとに `matrix`（`n_features_in`（まとめた後、recipe の前）・`n_features_out`（recipe の後）・段ごとの報告 `steps`）、
`bh_population_n`、`detection_available`（検出の例数を出せたか）・`n_detection_driven`・`n_constant_in_group`（要約のみ）、
集約なら `aggregate`（和: `members`・`summed`・`dropped`、まとめがあれば `merged_from`、寄与 `contributions` と
`n_samples_for_shares`。比: `numerator`・`denominator`・`dropped`）、OLS なら `covariate_excluded_samples`。

**和の寄与（`aggregate.contributions`）。** 和の各メンバー（`summed` の順）が和に占める割合（試料ごとのメンバー ÷ 和。
値は変換の前）の要約 `{feature_id, share_median, share_median_ref, share_median_test, share_max}`。使うのは比べる 2 群の試料の
うち和が有限で正のもの（その数が `n_samples_for_shares`）。`share_median` は 2 群を合わせた中央値、`_ref` / `_test` は群
ごと、`share_max` は最大。判定には使わない（中央値の和は 1 とは限らない）。

結果の全体に `spec`（正規化後）・`rule_hash`・`merge_split`（`enabled`・`source`: `find_split_pairs` か `alignment_qc`・
`groups`: まとめた組（残す ID が先頭）・`groups_hash`・`n_groups_given`・`skipped`）・`input_fingerprint`（`dataset_hash`・
入力の TSV と `alignment_qc` の JSON ごとの sha256）・`comparison`、外部の表（MAF・mwTab・xlsx）の入力なら `input_read`
（`{source, n_features, n_samples, n_unreadable_cells, unreadable_examples}`。xlsx は `n_rows_skipped_as_text`・
`rows_skipped_as_text` も）、`correlate` を書いたなら `correlate`、`notes`（捨てた設定と、§3 の `drift_correct` の注入順の
食い違い）。
`rule_hash` は正規化した spec（パスを除く）に、データセットの入力で `merge_split` が真なら使った組の hash も入れる
（外部の行列では `merge_split` は hash に入らない）。

**戻り値（v1）。** `{status, schema, rule_hash, merge_split: {enabled, source, n_groups_used, n_groups_skipped, groups_hash},
input_fingerprint, input_read?, comparison, group_ref, group_test, output_path?, analyses, correlate?, record, notes?, note}`。
`analyses[]` は解析ごとの要約（`n_rows`・`n_tested`・`n_meets_thresholds`・`n_constant_in_group` ほか）。`output_path` を
省けば解析ごとに行を TSV（`rows`。上の行の欄の順）で返し、**50 行を超える解析は p の小さい順の先頭 50 行だけ**を載せる
（`truncated: true`・`rows_order: "p_ascending"`。`n_rows` は全件。p が null の行は後ろ。断らない）。`correlate` の行も TSV で、
50 を超えれば先頭だけ（`truncated`）。全量は `output_path` の JSON。

集約の要約（`analyses[].aggregate`）は `output_path` の有無によらず、載せたリストごとに件数 `n_<名前>`（`n_members`・
`n_summed`・`n_dropped`・`n_merged_from` のうち載るもの。和なら `n_contributions` も）を持ち、リストをそれぞれ自分の長さで
断らずに切る:

- `contributions`（TSV）は、和に入ったメンバー（`n_contributions` ＝ `n_summed`）が 50 を超えたときだけ `share_median` の大きい順
  （null は後ろ、同じなら `summed` の順）の先頭 50 行（`contributions_order: "share_median_desc"`）。50 以下なら `summed` の順の
  全行（`contributions_order` は無い）。
- `members`・`summed`・`dropped`・`merged_from` は 50 件を超えたものだけ先頭 50 件。
- 切ったものに `<名前>_truncated: true`（切らなければ false）。`notes` は切ったものだけを名指す（寄与の並べ替えは寄与を
  切ったときだけ書く）。

切るのは戻り値の要約だけで、`output_path` の JSON・`rule_hash` は全メンバーのまま。

**記録（v1・v2 とも）。** 呼ぶたびに `reports/stats/<record_id>.json` を 1 件書く（`lipidmix.stats-record.v1`、`step:
"analysis_run"`。`record_id` は時刻と乱数で、呼び直しても上書きしない）。`inputs` はハンドルの身元（`dataset_handle`・
`source_kind`・`source`〔入力の名前だけ〕・`polarity`・`assay_kind`）と `input_fingerprint`、`settings` は `spec`（渡したまま）・
`normalized_spec`・`rule_hash`、`summary` は戻り値（`record` を除く）、`created_at` が時刻。戻り値の `record` は
`{written: true, record_id, path}`、書き先が決まらない・書けなければ `{written: false, reason}`（解析は成功）。

### 2.1 統計エンジンの spec（analysis-spec.v2）と結果（analysis-result.v2）

`spec["schema"] == "analysis-spec.v2"` の `dataset_analysis_run` は統計エンジン（`metabolomix/analysis/stats/`）で動く
（v1 の表・行の欄とは別の系。`analysis-spec.v1` の spec は v1 が実行する）。前処理を **1 回**回し、`analyses` の順にすべての
解析を実行して**全部**を返す（`primary` は報告の主役を指すだけ）。どれかの解析が失敗したら全体を失敗にする。**状態も
記録も書かない**（`result_id` も `matrix_id` の登録も作らない）。

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
               "effect_size": "log2_mean_difference", "bh_population": "matrix",
               "verdict": {"q": 0.05, "equivalence_fc": 1.2}, "min_detected_fraction": 0.5,
               "allow_confounded": true}],
 "primary": "main"}
```

| 欄 | 既定 | 意味 |
|---|---|---|
| `input` | `{"source": "dataset"}` | `source` はデータセットだけ。`exclude_samples` / `exclude_features`（空でない文字列のリスト。重複は順序を保って落とす）と `group_overrides`（試料 → 群）は**この spec に書いたものだけ**が効く（下の「自己完結」） |
| `preprocess` | 上の値 | `dataset_preprocess` の引数と同じ設定。`drift.mode` は `auto` / `on` / `off`、足切りは `null` でその足切りを行わない、`normalize.method` は `auto` / `none` / `tic` / `median` / `pqn`。任意の `drift.loess_degree`（`1` 既定 / `2`）と `normalize.features` / `normalize.exclude_features`（正規化の係数を求める特徴 ID）は `dataset_preprocess` の `loess_degree` / `normalization_features` / `normalization_exclude_features` と同じで、既定の値なら spec に現れない。`normalize.method="none"` の `features` / `exclude_features` と `drift.mode="off"` の `loess_degree` は捨てて `notes`、空の `features` は書かないのと同じ。`injection_columns` は個体のまとめの鍵から外す試料の表の列（取得日・注入番号など注入ごとに違う属性）。不正なら `ANALYSIS_SPEC_INVALID`（`details.field` が `preprocess.` で始まる） |
| `design` | `{"factors": ["group"], "covariates": [], "block": null}` | `factors`（重複のない文字列。`group`・`batch` か、試料の表の列名。sample-manifest v2 の後ろに足した因子の列も使える。交互作用の対比は 2 因子のモデルだけ）・`covariates`（`{name, type: numeric \| categorical, reference?}`）・`block`（`"subject_id"` など） |
| `comparisons` | 必須 | 対比のリスト。`{id, factor, reference, test, at?}`（`at` は他の因子の水準で層別）か、交互作用 `{id, kind: "interaction", factors: [2 因子], levels: {因子: [基準, 比較]}}`。`id` は英数字と `_ . -` で重複不可 |
| `analyses[].kind` | 必須 | `compare` / `composition` / `pca` / `plsda`。`analysis_id` は重複不可 |
| compare / composition | | `test`（既定 `moderated_t`。`welch` / `student` / `mannwhitney` は 1 因子・共変量なし・ブロックなしだけ）・`effect_size`（`log2_mean_difference` 既定 / `arithmetic_mean_ratio`）・`bh_population`（`matrix` 既定 / `annotated` / `reported`。`reported` は `reported_feature_ids` が要る。ほかの母集団で書いた `reported_feature_ids` は捨てて `notes`）・`p_adjust`（多重比較の補正。`bh` 既定 / `bonferroni` / `holm` / `none`。`q` はこの方法の値。既定の `bh` は正規化した spec に書かない＝明示しても同じ hash）・`verdict.q`（既定 0.05）・`verdict.equivalence_fc`（既定 1.2）・`min_detected_fraction`（既定 0.5）・`allow_confounded`（既定 true。完全交絡でも比べて `confounding` と行の `confounded` に印を付ける。false なら完全交絡で `STATISTIC_CONFOUNDED`。正規化した spec にいつも書く） |
| pca | | `scaling`（`autoscale` 既定 / `pareto`）・`n_components`（既定 5）・`impute`（`half_min` 既定 / `knn` / `none`） |
| plsda | | PLS 判別分析。クラスは `design.factors` の水準の組（cell。`classes` に cell 名を並べて絞れる。2 つ以上）。`n_components`（既定 2）・`scaling`（`autoscale` 既定 / `pareto`）・`impute`（`half_min` 既定 / `knn` / `none`）・`cv_folds`（層化グループ k 分割の数。既定 7、2 以上。上限は独立な単位＝グループの数）・`n_permutations`（並べ替え検定の回数。既定 200、1 以上。10000 を超えれば 10000 に下げて `notes`。正規化した spec と `rule_hash` も 10000）・`seed`（分割と並べ替えの乱数。既定 0） |
| `primary` | compare / composition が 1 つならそれ | 2 つ以上あるときは必須（`ANALYSIS_SPEC_INVALID`）。pca / plsda は選べない |

**自己完結（再現のための規則）。** v2 は `dataset_exclude` の除外を読まない。spec の `input` に無い除外がデータセットに
あれば、結果の `warnings` に `{"code": "dataset_exclusions_not_in_spec", "n_samples", "n_features"}`（効いていない除外の件数）。
組成のクラスは種別によらず MS-DIAL の ontology（特徴の `ontology`。無い特徴は `unclassified_features` に数える）で切る。
ハンドルの `assay_kind`（`record_objective(assay_kind=...)`。spec には載らず、戻り値の `assay_kind` に出る）が変えるのは
ラベル `class_label_kind`（`lipid` なら `lipid_class`、ほかは `msdial_ontology`）だけで、クラス・CLR・絶対量の数字と
`rule_hash` は種別によらず同じ。

**前処理は尋ねない。** `drift.mode="auto"` は、補正できても内部標準で換算済みの値か、注入順が MS-DIAL の既定（ファイルの
順）かもしれないときは補正を当てない。そのとき前処理の `steps` のドリフトの段に `reason`・`reasons`・
`rerun_with: {"preprocess.drift.mode": "on"}`、`warnings` に `drift_not_applied_by_default` が載る（補正するなら
`preprocess.drift.mode` を `"on"` にして呼び直す。`rule_hash` が変わる）。
陽陰を統合したデータセット（`dataset_merge_polarities`）では、ドリフト補正と正規化を極性ごとに当てる（`dataset.md`
§10.2。段の `by_polarity`・注意 `polarity_runs_differ`）。spec の欄は増えず `rule_hash` は変わらないが、
`input_fingerprint.sample_table_hash` に陰イオンの測定（`quant_state.polarity_runs`）が、`matrix_id` に極性の塊が入る。

**戻り値は要約だけ。** `{status, schema: "analysis-result.v2", rule_hash, matrix_id, input_fingerprint: {dataset_hash,
sample_table_hash}, primary, warnings, analyses, record}`（`record` は §2 の「記録」）。`analyses[]` は解析ごとの要約:

- compare: `n_rows`・`n_tested`・`verdicts: {changed, unchanged, inconclusive}`・`n_q_below`（`q ≤ 閾値`）・`test`・`warnings`・
  `confounding`（群とバッチの交絡 `{batch_status, by_contrast}`）・`p_adjust_method`、要因の F があれば `factorial`（効果ごとの
  `{df1, n_tested, n_q_below}`）。
- composition: compare の欄に `tables`（`clr` / `absolute` ごとの同じ件数）・`class_label_kind`・`n_classes`。`confounding` は
  内側の比較のもの（`absolute`、無ければ `clr`。同じ個体なので同じ値）。`warnings` は組成の注意に内側の比較の注意（完全交絡・
  注入順の偏り・`constant_in_group`・`units_missing_model_values` ほか）を足し、コードで重ねない。内側から来た注意は
  `tables`（`absolute` / `clr`）にどの表のものかを持つ。
- pca: `scaling`・`imputation`・`n_fit`・`n_units`・`n_qc_projected`・`n_features`・`explained_variance_ratio`。
- plsda: `classes`（cell → 個体数）・`n_units`・`n_features`・`scaling`・`imputation`・`n_components`・`r2x_cum`・`r2y_cum`・
  `q2_cum`・`q2`・`cv`・`permutation`・`vip_top`（VIP の上位 10）・`warnings`。

pca と plsda は行を持たない。長い一覧は切る。全量（行・CI・事前分布・交絡の表・PCA のスコア）は `output_path` の JSON。
`output_path` を省けば行を TSV で返し（列は `contrast_id feature_id status reason log2fc ci_low ci_high p q verdict`、
composition は先頭に `table` 列）、**50 行を超える解析は p の小さい順の先頭 50 行だけ**を載せる（`truncated: true`・
`rows_order: "p_ascending"`・`notes`。断らない）。`rule_hash` は正規化した spec の hash、`matrix_id` は前処理の入力と設定の hash（同じなら同じ ID）。

**`analysis-result.v2` の全量（`output_path` の JSON）。** `{schema, spec（正規化後）, rule_hash, input_fingerprint, matrix_id,
preprocess: {steps, warnings, subject_source}, primary, analyses: [{analysis_id, kind, result}]}`。`result` は kind ごとに次の
スキーマ。

- **`stats-compare.v1`**（compare）: `{matrix_id, test（実際に使った検定名）, effect_size_definition: {name, description},
  ci_level, p_adjust_method, verdict_rule: {q_threshold, equivalence_fc, equivalence_log2}, bh_population: {name,
  n_by_contrast, p_adjust_method, q_interpretation}, presence, prior, design, group_n, subject_source, confounding,
  warnings, rows, omnibus, factorial}`。
  - `q_interpretation`: `discovery_fdr`（`p_adjust=bh` で `matrix` / `annotated`）・`fwer`（`bonferroni` / `holm`）・
    `unadjusted`（`none`）・`reproduction_only`（`reported`。報告した特徴だけの補正で、探索の FDR ではない）。
  - `omnibus`: cell が 3 つ以上のときの全 cell の F 検定の行 `{feature_id, f, f_p, f_q, df1, df2}`（無ければ null）。`df2` は
    `df_prior + df_residual` で、事前分布の自由度が ∞ なら `"inf"`（F の参照分布は χ²(df1)/df1 の極限）、F を計算していない
    行では null。
  - `factorial`: 因子が 2 つ以上で全 cell がそろうときの主効果と交互作用の F（Type III。多水準を含む moderated F）。
    `{effects: {効果の名前: {factors, df1}}, hypothesis, rows: [{effect, feature_id, f, f_p, f_q, df1, df2}]}`（無ければ null。
    効果の名前は因子名、交互作用は `:` で結ぶ。`f_q` は効果ごとに `p_adjust` で補正）。
  - `confounding`: `{batch_table, batch_status, cramers_v, injection_order, warnings, by_contrast: {対比 ID: {batch_status,
    confounded}}}`（全対比をまとめた表と対比ごとの表。`confounded` は `true`＝完全交絡 / `false` / `null`＝バッチ不明で評価
    できない。交絡なしではない）。
  - 行（対比 × 特徴）: `feature_id`・`contrast_id`・`reference`・`test`・`status`（`tested` / `not_tested`）・`reason`・
    `n_reference` / `n_test`（有限値の個体数）・`detected_reference` / `detected_test`（検出された個体数。検出状態が無ければ
    null）・`quantity_judgeable`（群のどの cell も検出の割合が `min_detected_fraction` 以上か。null は検出状態なし）・
    `constant_in_group`（埋め値の疑いの cell の一覧。cell の個体の**読んだ値**＝注入ごとの値（ドリフト補正・正規化・反復の
    まとめの前）の有限の正値が 2 注入以上あり、すべて同じ（相対の丸めの誤差の内）なら挙げる。v1 と同じ規則。無ければ null。
    判定は変えず、件数を警告 `constant_in_group` に出す）・`confounded`（その対比が群とバッチの完全交絡か。評価できなければ
    null）・`log2fc`（正なら test が高い）・`se`（log2fc の標準誤差。順位の検定・算術平均の比では null）・`ci_low` / `ci_high`・
    `ci_unavailable_reason`・`statistic` / `df` / `p`・`mean_reference` / `mean_test`（個体ごとの前処理後の線形値の算術平均。
    単純な 2 群の対比だけ、交互作用は null）・`q`・`presence_p` / `presence_q`（検出の有無そのものの検定）・`verdict`・
    `verdict_reason`。
- **判定（`verdict`）は 3 値**で、規則は機械的: `changed`＝検定でき、`q ≤ verdict.q`。`unchanged`＝`q` が閾値を超え、95% CI が
  同等の範囲（`|log2fc| < log2(equivalence_fc)`）に完全に収まる。それ以外はすべて `inconclusive`（差があるともないとも言えない。`verdict_reason`:
  `not_tested` / `quantity_not_judgeable` / `outside_bh_population` / `no_ci` / `neither`）。**「有意でない」は「変わって
  いない」ではない**: `q` が大きいだけでは `unchanged` にならない。n が小さい・ばらつきが大きいときの大半は `inconclusive`。
  論文を再現したかの判定ではない。
- **`stats-pca.v1`**（pca）: `{matrix_id, scaling, log2: true, imputation, n_fit, unit_ids, qc_ids, scores: {units, qc},
  explained_variance_ratio, loadings, kept_features}`。軸は生体の個体だけで決め、QC は同じ軸へ**投影**する。欠測は PCA の
  ためだけに補い（`imputation`）、行列は変えない。
- **`stats-plsda.v1`**（plsda）: `{matrix_id, factors, classes: {cell: 個体数}, y_coding, n_units, unit_ids, scaling,
  imputation, n_features, n_features_dropped_constant, n_components_requested, n_components, r2x_cum, r2y_cum, q2_cum, q2,
  cv: {method: "stratified_group_kfold", group_by, n_groups, class_groups, folds, seed, scaling_within_fold: true,
  q2_definition: "1 - PRESS/TSS"}, permutation: {n, n_finite, p, statistic: "q2", scheme, min_attainable_p,
  q2_null_mean, q2_null_q95}, vip: [{feature_id, vip}]（降順）, scores: {unit_ids, values}, warnings}`。成分数・交差検証の
  Q2・並べ替え検定の p を必ず返す（R2Y やスコアの分離だけでは過学習と区別できない）。
  - Q2 は層化グループ k 分割で、分割ごとに学習の個体だけで中心と尺度を決めて予測した `1 − PRESS/TSS`（`q2_cum` は成分
    1〜A それぞれ、`q2` は A 成分）。並べ替え検定はクラスの札を入れ替えて同じ分割で Q2 を出し直した帰無分布で
    `p = (1 + #{Q2_perm ≥ Q2}) / (1 + n)`。
  - **分割と並べ替えは独立な単位（グループ）で行う**: `design.block` があればその値、無ければ `subject_id`、それも無い単位は
    1 つずつ（`cv.group_by` = block の名前 / `subject_id` / `unit`、`n_groups`、クラスごとのグループ数 `class_groups`）。分割は
    同じグループの単位を同じ分割に入れ、グループをクラスの組で層化して配る。分割数と成分数の上限はグループの数で決まる。
    並べ替えは、1 つのクラスだけのグループは札をグループごと、複数のクラスにまたがるグループ（対応のある計画・交差した計画）
    はグループの中で入れ替える（`permutation.scheme`: `between_groups` / `within_groups` / `mixed`）。各単位が 1 グループなら
    単位ごとの層化 k 分割と同じ。
  - `min_attainable_p` は並べ替えの p がおよそ下回れない値（観測と同じ Q2 になる札の並びの数 / 入れ替えで作れる並びの数。
    2 クラスは名前の入れ替えも同じ Q2）。
  - VIP は全個体の当てはめの変数重要度（二乗の平均が 1。目安の 1 は経験則で、検定ではない）。Y は 2 クラスなら後のクラスを
    1 とした 1 列、3 クラス以上はクラスごとの指示変数。全個体で値の変わらない特徴は外す。
  - 警告: 成分数・分割数を個体数と特徴数で頭打ちにした（`n_components_capped` / `cv_folds_capped`）、Q2 が正でない
    （`plsda_q2_not_positive`。VIP は判別の根拠にならない）、`min_attainable_p` が 0.05 を超える（`permutation_few_labelings`
    `{n_labelings, min_attainable_p}`。p が大きくても判別が無いとは言えない）、単位に `subject_id` が無く前処理が反復注入らしい
    名前（`possible_technical_replicates`）を見つけた（`plsda_units_may_not_be_independent`。Q2 と p が甘く出うる）。因子・
    block の値の欠けた個体は外す（`units_missing_model_values`。`variables` に欠けた変数）。
  - エラー: `PLSDA_SPECIFICATION_INVALID`（引数）・`STATISTIC_GROUP_TOO_SMALL`（2 個体以上のクラスが 2 つ未満、またはグループが
    1 つだけのクラスがある。`details.class_groups`・`group_by`）・`PLSDA_MISSING_VALUES`（`impute=none` で欠測あり）・
    `PLSDA_TOO_FEW_FEATURES`（値の変わる特徴が 2 未満）。
- **`stats-composition.v1`**（composition）: `{class_label_kind, classes, unclassified_features, zero_handling, clr, absolute,
  unit, warnings}`。`clr` / `absolute` はそれぞれ stats-compare.v1。クラスの合計が 0 以下の個体があるクラスは `clr` から除く
  （`zero_handling.excluded_classes`）。欠測のメンバーは合計に 0 として入り、その件数と警告
  `class_total_includes_missing_members` を出す。
- 中の `warnings[]` は `{code, message, ...}`。長い一覧（`samples` など）は要約で切る。

対話の `dataset_compare` / `dataset_pca` が保持する結果は、同じ stats-compare.v1 / stats-pca.v1 を `ds.results[result_id]` の
記録（`kind` が `compare` / `stat_pca`）の `engine` に持つ。対話の要約と記録（`reports/stats/`）の形は `dataset.md` の §10・§11。

## 3. 行列の recipe（matrix-recipe.v2）

`steps` に挙げた段だけを、挙げた順に走らせる（挙げない段の設定は捨てて `notes`）。篩の段（`detection_filter`・
`qc_rsd_filter`）で落とした特徴はその直後に行列から外す（後の KNN の近傍・PQN の参照に混ざらない）。既定（正規化で埋まる）:

| 段 | 既定 | 定義 |
|---|---|---|
| `detection_filter` | `{"basis": "detected", "min_rate": {}}` | 役割（`qc`・`sample`・`blank`）ごとの検出率の最小。挙げた役割のすべてを満たす特徴を残す。`basis`: `detected` は gap-fill を除く実検出（検出状態が無ければ `DETECTION_UNAVAILABLE`）、`positive` は値が有限かつ > 0（gap-fill を含む） |
| `impute` | `{"method": "half_min", "missing": "zero", "k": 10, "scale": "log2", "rowmax": 0.5, "colmax": 0.8}` | 欠測（`missing`: `zero` は ≤0 か非有限、`nan` は非有限だけ）を埋める。`half_min`: 特徴ごとの正の値の最小の半分（正の値が無い特徴は欠測のまま）。`knn`: 下記。`none` |
| `normalize` | `{"method": "pqn", "reference": "qc"}` | `pqn`: 参照（`reference` の役割の試料の、特徴ごとの正の値の中央値）に対する商の中央値で割る。`tic`（有限値の和）・`median`（正の値の中央値）は係数を試料間の中央値で割る。係数が 0・非有限の試料は割らずに報告する。陽陰を統合したデータセットでは、係数を極性ごとにその極性の特徴だけから求め、その極性の列だけを割る（段の報告は `by_polarity`: 極性 → `factors`・`unscaled_samples`） |
| `drift_correct` | `{"method": "loess", "span": 0.75, "iterations": 0, "min_qc": 4}` | 特徴ごとに、正の QC 値を注入順に対して LOWESS で当て、`値 / 当てはめ × QC 値の中央値` に直す。正の QC 値が `min_qc` 未満、または当てはめが正の有限でない特徴は直さず数える。注入順はすべての試料に要る（`INJECTION_ORDER_MISSING`）。陽陰を統合したデータセットでは極性ごとにその極性の注入順で当てはめる: `injection_order.source=metadata` は陽イオンを試料の表（当てたシートを含む）、陰イオンを統合が残した陰イオンの測定（`quant_state.polarity_runs.negative`）の注入順で（その極性の注入順の無い試料があれば、message の先頭と `details.polarity` に極性を添えて `INJECTION_ORDER_MISSING`）。TSV の注入順は全極性に当て、陰イオンの測定の注入順と違う（分からない）試料があれば、結果と戻り値の `notes` にその数と試料を書く（陰イオンをその測定の順で補正するなら `metadata` か v2・`dataset_preprocess`）。段の報告の `by_polarity` は極性 → `uncorrected_features`・`injection_order`（使った注入順: `metadata` / `tsv` / `polarity_runs.negative`） |
| `qc_rsd_filter` | `{"max_rsd": 0.3}` | QC の有限値の `std(ddof=1)/mean` が `max_rsd` 以下の特徴を残す。RSD が出ない特徴も落とす（別に数える） |

**参照する標準の定義:**

- KNN は R の `impute::impute.knn`（Bioconductor）の `knnimp` の規則を写したもので、`maxp` = 特徴数（2 段の分割をしない）に
  当たる。特徴を行に、`scale="log2"` なら正の値の log2 の上で: 欠測の割合が `rowmax` を超える特徴（重い行）は重い行どうしの
  列の平均で埋め、近傍にしない（重い行の中にその列の観測が無ければ軽い行の列の平均）。残りの特徴（軽い行）は軽い行だけで、
  欠測のある特徴ごとに、その特徴が観測された列のうち相手も観測された列での平方差の平均を距離にし（共通の列が無い相手は
  除く）、近い `k` 個を選ぶ（同距離は特徴の並び順）。欠測の列は近傍のうちその列を観測しているものの平均、全員が欠測なら
  軽い行の列の平均。`colmax` を超える試料があれば `RECIPE_PRECONDITION`。特徴が 1500 を超えるデータでは R の既定
  （`maxp=1500`。k-means で分割してから補完）の結果とは一致しない。R の実出力とは突き合わせていない（手計算の例でテスト）。
- LOWESS は statsmodels 0.15 の `statsmodels.nonparametric.smoothers_lowess.lowess(endog=QC 値, exog=QC の注入順,
  frac=span, it=iterations, delta=0, xvals=全試料の注入順)`（1 次の局所重み付き回帰・tricube の重み。`iterations=0` は
  頑健化の反復なし）。R の `loess`（既定は 2 次）とは別物。QC が少なく当てはめが NaN になる特徴は補正しない。

## 4. MS/MS の合否（msms-rules.v1）と診断イオン

`library_msms_qualify(rules, spot_id | precursor_mz, ..., mz_tol_ppm)`。照合は `library_match_feature` と同じ（採点の芯は
`library/matching.py`。前駆体の窓は MS-DIAL と同じく m/z 500 を超えたら同じ ppm に広がる。`mz_tol_ppm`（任意）はその窓に
重ねる ppm の絞り込みで、使った窓は `query` の `mz_tol`・`mz_window`・`mz_tol_ppm`。`ion_mode` は `pos` / `neg` / `+` / `-`
も読み、読めない値は `DATASET_BAD_REQUEST`）。規則（正規化後。既定はこのとおり）:

```json
{"schema": "msms-rules.v1",
 "target": {"inchikey_block": "HVYWMOMLDIMFJA", "formula": null},
 "top_rank": 5, "min_fragments_excl_precursor": 1, "min_fragment_relative_intensity": 0.01,
 "precursor_exclusion_tol": null, "same_formula_margin": null, "lipid_class": null,
 "diagnostic_ions": [], "neutral_losses": [], "ions_match": "all"}
```

`min_fragments_excl_precursor` は 0 以上の整数か `null`（`min_fragments` の検査をしない）。`min_fragment_relative_intensity`
（0 以上 1 以下。測定の基準ピークに対する比。正規化した規則にいつも載る）: 断片の一致は、測定強度が基準ピークのこの比以上の
窓だけを数える。`matched_peaks_count`（MS-DIAL の写し）にはこの下限を当てない。

入力の `target` は `{"inchikey": <InChIKey 全体か先頭 14 文字>, "formula": <任意>}`。InChIKey が一意に決まらない脂質の種は
`{"lipid_species": "PC 18:3_16:1", "formula": <任意>}`（`inchikey` と排他。opt-in）で、正規化後は
`{"lipid_species": "PC 34:4", "formula": …}`（GOSLIN の種のレベル。読み方は `lipidmix://docs/output-format/mztab` §13.8）。
種の目標では、`target_in_top_rank` と `same_formula_margin` の「目標か別物か」をライブラリの記録の種
（`msms_rules.record_lipid_species`）で決める（同じ種の別の分子種 `PC 14:0_20:4` は目標の候補で、別物に数えない）。記録の
種は次の順に読み、最初に読めた経路を `species_route` とする:

| `species_route` | 読み方 |
|---|---|
| `name` | 記録の `name`（MS-DIAL の限定子を外して）をそのまま GOSLIN で読む |
| `name_cleaned` | `\|` の後・`; ` の後の注記（`PC 34:1; ExampleDB ID-0000`）・末尾の付加体の表記（`[M+H]+`）を外し、古い MS-DIAL の書式（`Cer[NS] d18:1/16:0`）・英語のクラス名（`Phosphatidylcholine 16:0-18:1`）を略記に直して読む（`pathway.msp_lipid.shorthand`）。`;` の直後に空白の無い修飾（`;O` `;O2` `;O3` `;OH` `;OOH` `;2OH` `;(2OH)` `;(3OH)` `;Ep` `;O(FA 16:0)` など。知らない語も）は注記が続いても残し、括弧の外の次の空白から後ろだけを外す（`PE 36:2;O (oxPE)` → `PE 36:2;O`）。GOSLIN が読めない修飾（`;(2OH)` `;(OH)`）を残した名前は種に読まず、次の経路へ回る |
| `ontology_formula` | 名前で読めなければ、ClassyFire の `Ontology` と分子式から和の組成の種を作る（`pathway.msp_lipid.from_class_formula`。LPS・LPG も） |

同位体の標識のある名前（`(d7)` など、内部標準）と、どの経路でも読めない記録は目標とも別物とも数えない。目標の
`target.lipid_species` に標識があれば `MSMS_RULES_INVALID`（`field: target.lipid_species`。
天然の種なら標識の無い名前を渡す）。名前の 2 つの経路（`name` → `name_cleaned`）は、目標照合（`dataset_feature_search` と束の
既定の特徴）が SME 候補・SML の名前を読むのとも、MS-DIAL の注釈（`msdial_rank`）を読むのとも同じ
（`msms_rules.name_lipid_species`）。目標の候補の経路は `target_in_top_rank` の `detail.species_route`、候補の表では
`lipid_species`・`species_route` の列に出る。

設定した検査が `checks` に並ぶ（各行 `{check, status: pass|fail|not_evaluable, value, detail}`。既定の規則では
`target_in_top_rank` と `min_fragments`）:

| 検査 | 合格の条件 | `value` |
|---|---|---|
| `target_in_top_rank`（常に） | 総合スコアの順で `top_rank` 位までに、InChIKey 先頭ブロックが目標と一致する候補がある（最上位のそれを「目標の候補」） | 目標の候補の順位（無ければ null） |
| `min_fragments`（既定で入る） | 目標の候補の、前駆体 ± `precursor_exclusion_tol`（省略時は採点の `ms2_tol`）を除き、測定強度が基準ピークの `min_fragment_relative_intensity` 以上の一致の本数 ≥ 最小（`detail.min_relative_intensity` に下限） | その本数 |
| `same_formula_margin` | 全候補のうち、目標と同じ分子式（`target.formula`、無ければ目標の候補の分子式）で先頭ブロックが違うものの最高点との差 ≥ `min_delta`（`score` 既定 `total_score`）。別物が無ければ pass | 差（目標 − 別物） |
| `lipid_class` | 1 位の候補の名前を GOSLIN で読んだクラス（正規化名の最初の空白の前）が一致 | 読めたクラス |
| `diagnostic_ions` / `neutral_losses` | 測定スペクトルに ±`tol`（既定 0.01 Da）のピークがあり、最大ピークに対する相対強度 ≥ `min_relative_intensity`。`ions_match` が `all` なら全部、`any` なら 1 つ。中性損失は前駆体 − `delta`（生成イオン）を探し、`min_relative_intensity` はその生成イオンの比 | 見つかった数 |

要素は診断イオンが `{mz, label, tol, min_relative_intensity}`（`min_relative_intensity` 既定 0）、中性損失が
`{delta, label, tol, min_relative_intensity}`（中性損失の `min_relative_intensity` は渡したときだけ規則に載る）。知らない
キー・欠けたキーは `MSMS_RULES_INVALID` で、誤りの文がそのキーと要素（`neutral_losses[0]` など）を名指す。

戻り値は `verdict`・`checks` のほかに、`target_rank`（目標の候補の順位。無ければ null）と `top1_is_target`（1 位が目標か。
`target_in_top_rank` は `top_rank` 位までなら pass するので、1 位かはこの欄で読む）。`spot_id` のスポットに MS-DIAL の注釈が
あれば、`library_match_feature` と同じ `msdial_name`・`msdial_rank`・`top_differs_from_msdial`（library.md §14.6）。候補の表
（`candidates`、TSV）の列は `rank name formula inchikey adduct delta_ppm total_score matched_peaks_count
fragment_matches_excl_precursor`（`delta_ppm` は `(測定 − 参照) / 参照 × 1e6`、`fragment_matches_excl_precursor` は合否と
同じ規則＝強度の下限つきで数える。`library_match_feature(detail=true)` の同じ名前の列も既定は同じ下限 0.01 で、使った下限を
`min_fragment_relative_intensity` に返す）。`.dbs` に注釈の足切り（dot・一致率・総合スコア・`minimum_spectrum_match`）が
あれば、値を `dbs_cutoffs` に 1 回載せ、候補ごとに `passes_dbs_cutoffs`（MS-DIAL の `IsSpectrumMatch` と総合スコアの足切りを
全部通るか。dot の足切りは二乗の値と比べる）の列を足す。**合否には当てない**（`.msp` では列が無い）。種の
目標なら `lipid_species`・`species_route` の列も足す。

`verdict`: すべて pass → `qualified`、どれか fail → `not_qualified`、fail が無く判定できない検査がある → `not_evaluable`。
目標の候補が無いと、それに依る検査は `not_evaluable`（`target_in_top_rank` が fail なので `verdict` は `not_qualified`）。
MS/MS が無い・候補が無いときは検査をせずに `not_evaluable`（「断片が合わなかった」ではない）。**「明らかに別物か」は判断
しない。** `rule_hash` は正規化した規則の hash。

`dcl_diagnostic_ions(spot_ids, ions)` は同じ照会（`find_ions`）をスポットの一覧に当てる。`ions` の各要素は `{"mz"}` か
`{"neutral_loss"}`（任意で `label`・`tol`・`min_relative_intensity`。中性損失の比は生成イオン＝前駆体 − `neutral_loss` の
測定強度の比。知らないキーは `MSMS_RULES_INVALID` で、誤りの文が `ions[i]` とそのキーを名指す）。行の `status` が `ok` 以外
（MS/MS が無い・索引の崩れ）のスポットは照会していない（`found` が空なのは「イオンが無かった」ではない）。

## 5. 同族列の RT（homolog-check.v1）

`dataset_homolog_check(spec)`。`members`（`feature_id`・`carbons`・`double_bonds`）か `auto`（`lipid_class`・`evidence` 既定
`["sme"]`・`adduct`）のどちらか一方と `rule`（`carbon_direction` 既定 `increasing`・`double_bond_direction` 既定
`decreasing`・`max_violations` 既定 0。逆相の LC の向き）。

- `auto` は SME rank 1 の名前（`sml` を足せば SML の名前も）を GOSLIN の種のレベルに読む。読み手は §4 と同じ
  `msms_rules.name_lipid_species`（`name` → `name_cleaned`。Console / mzTab-M の `PC 34:1; ExampleDB ID-0001` は `PC 34:1`）。
  同位体の標識のある名前（`PC 15:0_18:1(d7)` などの内部標準）は集めない。分子種（`16:0_18:1`）は鎖を足して種のレベル
  （`34:1`）に直る。エーテルは別のクラス（`PC O-`）で、`P-` は GOSLIN の規約でビニルエーテルの二重結合を数えて `O-` に直る
  （`PE P-18:0_20:4` も `PE P-38:4` も `PE O-38:5`）。プラスマローゲンの列は `lipid_class: "PE O-"` で集める（`"PE P-"` では
  集まらず、`notes` に `O-` の指定を書く。修飾つきの `"PE P-;O"` なら `"PE O-;O"`）。`adduct` は SME rank 1 の付加体と**綴りで
  なく組成で**比べる（`msdial.adducts.same_adduct`。`[M+H]1+` と `[M+H]+`、`[M+HCOO]-` と `[M+FA-H]-` は同じ）。
- **種のレベルの修飾（`;` の後ろ）は列を分ける。** 種の鍵は（修飾, 炭素数, 二重結合数）で、修飾の違う種どうしは重複にも
  隣どうしにもしない。例: `PE 36:2;O (oxPE)` は `PE 36:2;O`、`PC 16:0_18:2;OOH` は `PC 34:2;O2`、`Cer 18:1;O2/16:0`（NS）は
  `Cer 34:1;O2`、`Cer 18:0;O3/16:0`（NP）は `Cer 34:0;O3`。`lipid_class: "PE"` はそのクラスの修飾をすべて集めて修飾ごとの列で
  比べ、2 つ以上の修飾を集めたら内訳を `notes` に書く。`"PE;O"`・`"Cer;O2"` のように修飾を書けば、その修飾の種だけを集める。
  修飾は種のレベルの合計なので、同じ数の酸素を持つ別の亜クラス（Cer NP と AS はどちらも `;O3`）は同じ列になる。修飾のある
  種を集めたときは `members` に `modification` の列（無修飾は空）が付き、`duplicates` のキーと違反の `low` / `high` は
  `36:2;O`・`34:0;O3` のように修飾を添える（無修飾は `36:2`）。`members`（手で並べる）には修飾が無い。
- 同じ（修飾, 炭素数, 二重結合数）が複数あれば RT の中央値で並べ、`duplicates` に挙げる。
- 同じ修飾・同じ二重結合数の中で炭素数の昇順、同じ修飾・同じ炭素数の中で二重結合数の昇順に並べた**隣どうし**で、規則の
  向きに反する対を違反と数える（RT が同じ対は違反にしない）。比べた対の数が `n_pairs`。違反が `max_violations` 以下なら
  `pass`。違反に関わった特徴は `members` の `status` が `violation`。
- **比べた対が 0（`n_pairs: 0`）なら `verdict` は `not_evaluable`**（`pass` ではない。違反が無いという意味でもない）。`auto` が
  集めた種が 2 つ未満、RT のある種が 2 つ未満、または隣どうしになる種が無い（`34:1` と `36:2` だけ）とき。集めた件数と理由は
  `notes` に書く。

## 6. 特徴どうしの関係

`dataset_feature_relations(feature_id, rt_window=0.05, mz_tol=0.01, min_r=0.7, adducts=None)`。規則はアラインメントの
キュレーションの関係（`curation/relations.py`）と同じ。

- `coeluting`: |ΔRT| ≤ `rt_window` の特徴（`drt` = 相手 − この特徴、`r` = 試料間の log1p 強度の Pearson 相関。値が 5 組未満
  なら空）。
- `as_ion_of`: この特徴が、相手の付加体（`adducts` の組）・同位体・源内断片（`-H2O` ほか）として説明できる関係。同位体は
  M+1 と M+2（13C2・34S・37Cl・81Br。`mz_tol` が 13C2 との差 2.46 mDa より狭いときは 18O、8.59 mDa より狭いときは 41K も別の
  行。それより広い窓では、この特徴が 13C2 の窓に入ればどちらも 13C2 の行 `isotope_M+2` に含み、窓の外なら当たった M+2 の
  位置のうち最も近いときだけ 18O・41K の行を出す）。列は `of of_name relation dmz_mda drt r isotope_ratio expected_ratio
  isotope_consistent strong ambiguous informational`。
- `ions_of_this`: 相手がこの特徴のイオンとして説明できる関係。列は `ion relation dmz_mda drt r isotope_consistent strong
  ambiguous informational`。付加体の関係はこの特徴の付加体から計算するので、**注釈の無い相手も付加体・同位体・源内断片と
  して出る**。
- 付加体の関係は説明される側の付加体から計算する。その付加体は SME rank 1、無ければ SML 注釈（MS1 だけの同定）のもの。
  どちらも無い相手の `as_ion_of` には付加体の関係は出ない（同位体・源内断片だけ）。同位体の期待比に使う組成式は、SME の行 →
  SML 注釈 → 同定のある特徴の `formula`（解析フォルダの `.arf2`）の順に採る。
- **相手自身の付加体は綴りでなく組成で見分けて除く**（n_mol・価数・極性・付加と脱離の原子が同じなら同じ付加体。`[M+H]1+` と
  `[M+H]+`、`[M+HCOO]1-`・`[M+CHO2]-` と `[M+FA-H]-`、`[M+Hac-H]-` と `[M+CH3COO]-`）。相手と同じ m/z になる付加体
  （`[M+H]+` に対する `[2M+2H]2+`）も除く。
- **`expected_ratio` は相手のイオンの組成から求める**: 中性の組成式の n_mol 倍に付加体の原子を足した組成（`[2M+H]+` は炭素も
  2 倍、`[M+Cl]-` は Cl を 1 個持つ）の M+k の同位体の組み合わせ（M+1 は 13C・15N・2H・17O・33S・40K、M+2 は 13C2・18O・34S・
  37Cl・81Br・41K と交差項。IUPAC の天然存在比と NIST の質量差）のうち、**m/z がこの特徴から `mz_tol` 以内に来るものの、M に
  対する比の和**。
  - 分離された 1 つの M+2（Cl を持つイオンの 13C2 だけ、硫黄を持つイオンの 34S だけ）は、その位置の近くの同位体だけと比べる。
  - 分離できない束（重心が 13C2 と 34S の間）が束の全体と比べられるのは、この特徴から束のどの同位体までの差も `mz_tol` 以内の
    ときだけ（重心が真ん中なら約 5.5 mDa / 電荷以上。既定の 0.01 Da は覆う）。それより狭い窓では、重心は近い方の位置の 1 行
    だけに当たり、期待比が遠い方の同位体を欠くので、相手自身の束が `informational` になりうる（別物質の印ではない）。重心から
    近い方の位置までの差より狭ければ 13C2 と 34S の行は出ず、18O の位置までの差（約 2.7 mDa）以上なら 18O の行
    `isotope_M+2_18O` だけが立って、同じく `informational` になりうる。それより狭ければ M+2 の行は出ない。束を分離しない
    装置では `mz_tol` を約 0.006 Da / 電荷以上にする。
  - 同じ組の M+2 の行は同じ期待比。付加体が読めなければ中性の組成式のまま。組成式が無い・読めない・同位体の表に無い元素
    （Se・Fe など）を含むなら空（`isotope_consistent` も空）。I・As・Co は単一同位体として数える。表に無い元素を含んでも
    組成式が読めれば、O・S・Cl・Br・K を持たない相手にその M+2 の行は出さない。
- **`strong` は Δm/z が合うことだけでは立たない**。MS-DIAL のリンク（`.arf2`）は使わず、裏付けは試料間の相関 `r ≥ min_r`
  だけ。付加体・源内断片は Δm/z と同時溶出（`rt_window`）に加えて `r ≥ min_r` のとき、同位体はさらに強度比が期待比の 1.5 倍
  以下（`isotope_consistent=true`）で `ambiguous` でないときだけ `strong`。
- `ambiguous`: 同じ相手に 2 つ以上の M+2 の説明（13C2 +2.006710・18O +2.004245・34S +1.995796・37Cl +1.997050・
  81Br +1.997952・41K +1.998119）が同じ窓に入った行。一方に決めず並べ、どれも `strong` にしない。分離できない M+2 の束は、
  窓が重心から両方までの差以上なら 13C2 と 34S の両方に入る。相手のイオンの組成が分かっていて O・S・Cl・Br・K を持たなければ
  その行は出さない（`[M+K]+` は K を持つ）。18O・41K は `mz_tol` が 13C2 との差（2.46 / 8.59 mDa）以上なら、この特徴が 13C2 の
  窓に入るかぎり 13C2 の行が含み（既定の 0.01 Da で 41K の位置のピークは `isotope_M+2` の行で `dmz_mda` が約 −8.6）、並べて
  曖昧にはしない。13C2 の窓の外（既定の窓で 41K の位置から 1.41 mDa より低い側）なら、当たった M+2 の位置のうち最も近い
  ときだけ 18O・41K の行を出す（より近い 34S などで説明されるなら出さない）。
- `informational`: 同位体の強度比が期待比を超えた行（同位体だけでは説明できず、別物質が重なっている目印）。重なった M+2 の
  行は期待比（この特徴の窓に入る M+2 の和）が同じなので、そろって立つか立たないか。

### 6.1 試料ごとのピーク（`dataset_sample_peaks`）

`dataset_sample_peaks(feature_ids, dataset_handle=None, output_path=None)`: 割れた 2 スポットが 1 つの物質か、RT が注入ごとに
ずれていくかを試料ごとの値で見る読み取り。行は列 `feature_id, sample, rt, mz, height, detected` の **TSV**（列名 1 回、`rows` の
文字列。依頼の特徴の順 × データセットの試料の順）。戻り値は `{status, n_features, n_samples, n_rows, rows}`。数値は小数 6 桁
まで、空欄が null（gap-fill の `mz`）、`detected` は `true` / `false`。

- `rt`: その試料のピークの RT（分）。gap-fill のセルも `.arf` が持つ値（補間の位置）を返すので、**`detected: false` の `rt` は
  検出された位置ではない**。ドリフトは `detected: true` の行だけで見る。
- `mz`: 検出されたセルの実測 m/z。gap-fill は空欄（null）。
- `height`: データセットの強度の行列（`feature_matrix`）の値。`measure` で選んだ `.arf` の欄で、定量（`dataset_quant_apply`
  ほか）で値を変えたデータセットではその値。
- `detected`: `detected_mask`（gap-fill なら false）。
- 特徴の数に上限は無い。知らない特徴 ID・空のリストは `DATASET_BAD_REQUEST`。ID は int でも str でもよい。
- 戻り値に載せる行は **1,000 行まで**。超えて `output_path` が無ければ断らず、依頼の順の先頭から特徴ごとにそろえて（1 つの
  特徴の試料を途中で切らない）floor(1000 / 試料の数) 特徴分だけ載せる（試料が 1,000 を超えれば最初の特徴の先頭 1,000 行）。
  戻り値は `truncated`（いつも付く。超えなければ false）・`n_rows`（全件）・`n_rows_shown`・`n_features_shown`・`notes`。
- `output_path` があれば全量の TSV（同じ列・列名 1 回）をそのファイルへ原子的に書き（既にあるファイルは置き換える）、戻り値は
  要約 `{status, n_features, n_samples, n_rows, output_path, sha256}`（`rows` は付かない。1,000 行以下でも同じ）。ツールの注釈は
  `readOnlyHint=false`・冪等。
- 使えるのは試料ごとの RT（`per_sample_rt`）のあるデータセットだけ: 解析フォルダ（`.arf`）から読んだもの、または mzTab-M で
  読み、隣の `.arf` と同じアラインメントだと照合できたもの（検出状態を取り込めたとき。`.arf` の列を mzTab-M の試料の順に並べ
  替えて持つ）。照合できなかった mzTab-M と、統合などで作り直したデータセットは `DATASET_BAD_REQUEST`（`field: per_sample_rt`）。

## 7. 証拠の束（claim-evidence.v1）

`dataset_claim_evidence(spec, dataset_handle, output_path)`。1 主張について、走らせ直しと照合・候補の網羅の走査を 1 つの
JSON にまとめる。

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
  （`{"default": 0.0, "bands": [{"max_mz": 200, "offset": -8.0}]}`。中心は実測 − 理論の ppm、帯は理論 m/z が `max_mz` 未満の
  最初の帯）。既定の特徴と `scan` の全変種に効き、**窓の幅（`scan.ppm`）は変えない**。渡したときだけ `rule.ppm_offset` と
  `rule_hash` に入る。形の誤りは `CLAIM_SPEC_INVALID`。
- `rt_offset_min`（任意、分）は RT の窓の中心のずれ。`dataset_feature_search` の `rt_offset_min` と同じく中心を
  `rt_min + rt_offset_min` に置き `|RT − 中心| ≤ rt_tol_min` で見る（`mztab` §13.8）。既定の特徴と `scan` の全変種に効き、
  **窓の幅（`scan.rt_tol_min`）は変えない**。渡したときだけ（0 でも）`rule.rt_offset_min`・`rule.target.rt_offset_min`・
  `rule_hash` に入る。`target` に `rt_min` が無ければ捨てて `notes`（束と戻り値）に書く。有限の数でなければ
  `CLAIM_SPEC_INVALID`（`field: rt_offset_min`）。ずれの見積もりは `dataset_feature_search` の目標照合が返す `rt_offset`。
- `target` は `dataset_feature_search` の目標と同じ形（`target_id` は `claim_id`）。脂質の種の主張（InChIKey が無いもの）は
  `"target": {"lipid_species": "PC 18:3_16:1", "formula": …}`。`rule.target` に `lipid_species`（種のレベル）と
  `lipid_species_input` が載り、`rule_hash` は `lipid_species_input` を入れない（綴りの違う同じ種は同じ hash）。**束の既定の
  特徴は種のレベルの一致で選ぶので、分子種（どの鎖の組か）は確かめていない**。
- `scan` の `ppm` × `rt_tol_min` の全組で目標照合を走らせ直す（`adduct_tiers` を渡すと付加体も広げる）。省略時は目標の
  `ppm` / `rt_tol_min` の 1 組。
- `msms.rules` は msms-rules.v1 から `target` を除いたもの（target は主張の InChIKey か脂質の種と分子式）。`msms.mz_tol_ppm`
  （任意、0 より大きい数）は前駆体の窓に重ねる ppm の絞り込み（`library_msms_qualify` と同じ）で、渡したときだけ `rule.msms`
  と `rule_hash` に入る。形の誤りは `CLAIM_SPEC_INVALID`。
- `msms` / `relations` を省くと、その欄を作らない。

結果:

| 欄 | 意味 |
|---|---|
| `default_selection` | `target` をそのまま当てた目標照合（`dataset_feature_search` と同じ規則。ただし割れの組はいつも `find_split_pairs` を既定の閾値で回したもの。§9 の「最後の結果はハンドルに持つ」）。`status` が `resolved` なら `selected_feature_id` が規則で選んだ既定の特徴、決まらなければ `needs_input`（`reason`）。割れたスポットの組（§9 の `merge_groups`）は 1 つの候補として数え、`selected_feature_id` は残す ID、`merged_from` にその組（それ以外は null。走査も同じ）。`counts` は下の表（組は 1 つと数え、MS/MS の合否はどちらかのスポットが `qualified` なら 1） |
| `scans` | 走査ごとの判定（`counts` は `default_selection` と同じ意味）と、候補の特徴 ID |
| `candidates` | 全走査の候補の和（特徴の並び順）。`evidence`（`sme_match` など、走査で見えたもの）・`adducts`・`scans`（どの走査で見えたか。組の片方だけが質量の窓に入ったときに足した相手は `merge_partner`）・`split`（下）・`msms`・`relations`。`msms` は spot_id = SMF_ID の照合の上位 `top`（`library_msms_qualify` の候補の表と同じ欄。`passes_dbs_cutoffs`・種の欄も同じ条件で）と `qualify`（`target_rank`・`top1_is_target` を含む）、使った前駆体の窓 `window`（`mz_tol`・`mz_window`・`mz_tol_ppm`）、スポットに MS-DIAL の注釈があれば `msdial_name`・`msdial_rank`・`top_differs_from_msdial`。ライブラリが無ければ `{"status": "no_library", "spectrum_available": ...}`。データセットの同定が LBM だけ（記憶が `lbm_only`）なら照合に使える MSP が無く、候補ごとに `{"status": "not_evaluable", "reason": "lbm_only", "spectrum_available": ...}`（理由の文は束の `msms_not_evaluable_reason` に 1 回だけ）。MS/MS が無い・索引の崩れはその `status`。`spectrum_available` は代表スペクトルの有無で、ライブラリの有無によらず spot を照会して決める: `true`=ある、`false`=照会して無かった（`not_found` で `reason` が `no_spectrum`）、`null`=分からない（`not_found` の `reason` `no_arf2` / `no_dcl`・`spot_not_found`・`dcl_precursor_mismatch`・`spot_id_not_integer`、ハンドルに `.arf2` / `.dcl` が無いとき） |
| `truncated` | `max_candidates` で打ち切った候補の数（`output_path` を省いて候補を切った戻り値では意味が違う。下の「戻り値」） |
| `notes` | 効かずに捨てた spec の設定とその理由（`rt_min` の無い `target` の `rt_offset_min`）。無ければ欄ごと無い |
| `rule` / `rule_hash` | 正規化した spec とその hash。spec の `merge_split`（既定 true）は真のときだけ `rule` に載り、`rule_hash` には判定でまとめた組の hash も入る。目標は付加体の層を解決した正規化後のもの（層と理論 m/z が hash に入る。陽イオンで分子式を渡し `adduct_tiers` を省いた目標の既定の層は `dataset_feature_search` と同じで、`[M+CH3OH+H]+`・`[M+H-2H2O]+` を含む）、`msms` があれば既定を埋めた MS/MS の規則（§4）が hash に入る |
| `merge_split` / `merge_groups_hash` | 割れの組を判定で 1 つの候補にしたか（spec の `merge_split`）と、その組（§9 の `merge_groups`）の hash（false なら null） |
| `input_fingerprint` | `dataset_hash`（データの指紋）・`library_sha256` |

`counts`（`default_selection` と `scans[]`。どれも特徴の一意な件数）:

| 欄 | 意味 |
|---|---|
| `mass` / `rt` / `by_evidence` | `dataset_feature_search` の `counts` と同じ（質量の窓・RT の窓に入った数、水準ごとの規則の qualified 数） |
| `qualified` | 規則（質量・RT・注釈）で qualified になった特徴のうち、**MS/MS の合否 `msms.qualify.verdict` が `qualified`** の数。`msms` を頼んでいない・ライブラリが無いときは `null`（0 ではない）。`max_candidates` で打ち切った候補は数えない |
| `rule_qualified` | 規則（質量の窓・RT の窓・注釈）で qualified になった特徴の数（`dataset_feature_search` の `counts.qualified` と同じ数）。**MS/MS が取れているかどうかの数ではない** |
| `msms_available` | **非推奨の別名**（`rule_qualified` と同じ値。次の版で消す）。名前に反して MS/MS の有無の数ではない |
| `msms_spectrum_available` | 質量の窓に入った候補（`mass` と同じ母集団。`merge_partner` の行を除き、割れの組は 1 つ）のうち、組のどれかのスポットに**代表 MS/MS がある**（`msms.spectrum_available` が `true`）数。RT の窓では絞らない。**0 は照会して無かった**（全部の単位が `false`）、**`null` は頼まなかった（`msms` なし）か分からなかった**（`true` の単位が無く、`null` の単位がある）、**数は下限**（`null` の単位は数えない）。単位（組）の中は、どれか `true` なら `true`、`true` が無く `null` があれば `null`、全部 `false` なら `false`。`max_candidates` で打ち切って照会していない候補は数えない |

`split`（候補が割れたアライメントスポットの組に入るときだけ。入らなければ `null`。既定の閾値の `find_split_pairs` の
`split_partners`（§9）を写す）:

| 欄 | 意味 |
|---|---|
| `group_id` | 全体の連結成分（§9 の `groups[].group_id`） |
| `partners` | まとめる相手（§9 の `merge_groups`。互いに最も近い mergeable な組の相手 1 つ）。まとめない特徴（`group_aligned`・最も近くない組など）は空 |
| `type` | まとめる特徴はその組の型、まとめない特徴は成分の型（`group_aligned` > `interleaved` > `injection_block` の順で強い方） |
| `mergeable` | 統計の側（`dataset_feature_search` / `dataset_analysis_run` の `merge_split`）が 1 つにまとめるか |
| `component` | 成分のほかの特徴（報告用。まとめない割れの相手もここで分かる） |

**束は特徴を決めない。** 最終の特徴を決めた後の統計は `dataset_analysis_run`（`features.feature_ids`）で別に出す。同じ spec・
同じ入力なら同じ JSON になる。

**戻り値。** `output_path` を省けば束をそのまま返す。候補が 10 を超えても断らず、束の候補を
`candidates_order: "default_selection_first"` の順（既定の特徴（`default_selection` の `selected_feature_id` と `merged_from`）→
既定の照合で見えた候補（`scans` に `default`）→ ほかの候補。同じ段の中は束の並び）に並べた先頭 10 件だけを `candidates` に
載せ、`truncated: true`・`n_candidates`（束の候補の全件＝`output_path` が書く数）・`n_candidates_shown`・`candidates_order`・
`notes`（`output_path` で全件を書ける旨）を返す。このとき束の `truncated`（`max_candidates` で打ち切った数）は
`max_candidates_truncated` に移る（0 でなければ `notes` にも書く。打ち切った候補は `output_path` の JSON にも入らないので、
要るなら `max_candidates` を増やす）。候補の一覧のほか（`default_selection`・`scans` の `candidate_feature_ids` を含む判定・
`rule`・`rule_hash`・`input_fingerprint`）は束と同じ。`output_path` があれば JSON を書き、要約（`claim_id`・`rule_hash`・
`input_fingerprint`・`default_selection`・`n_candidates`・`truncated`（束と同じ `max_candidates` の打ち切りの数）・走査ごとの
判定、捨てた設定があれば `notes`）だけを返す。

## 8. 質量精度と新しい知見の候補

`dataset_mass_accuracy(band_edges=[200], rank=1)`: SME の `rank` 位の候補の (実測 − 理論)/理論 × 10⁶ を、理論 m/z の帯ごと
（`[min_mz, max_mz)`）にまとめる。`n`・`median`・`q1`・`q3`（`statistics.quantiles(n=4)` の exclusive 法）・`within_5ppm`・
`within_10ppm`（|ppm| ≤ 5・10 の割合）。n<2 の帯は `n` だけ。

**理論 m/z の出所（`theoretical_mz_source`）**:

- `mztab`: mzTab-M の SME の `theoretical_mass_to_charge`（上流の Reference m/z＝ライブラリの記録のプリカーサー m/z。job /
  pipeline を解析フォルダで読んだときも、照合した mzTab-M から写した値）。
- `formula_adduct`: mzTab-M と照合しない解析フォルダ（`load_dataset(directory=...)`。`.arf2` は Reference m/z を持たない）。注釈の組成式（`Formula`）と付加体（`AdductType`）から計算したモノアイソトピックの m/z。**注釈を表す m/z で、ライブラリの記録の
  プリカーサー m/z ではない**（記録の値が丸めてある・実測値であるライブラリでは `mztab` と少し違う）。組成式か付加体が無い・
  読めない同定は数えない（件数は読み込みの注記）。
- `mixed`: 出所の違う 2 極性を `dataset_merge_polarities` で統合したもの。

`dataset_novel_candidates` の行の `ppm` も同じ理論 m/z から出し、同じ `theoretical_mz_source` を書く。根拠バンドル（§10）の
`msdial_sme` 行の来歴には計算した値を写さない（`null`）。

**許容幅は当てない**が、材料として帯ごとの `recommended_mz_tol_ppm` と全体の `recommended_tolerance`
（`{n, basis, mz_tol_ppm, median_ppm, mad_ppm}`）を添える。幅は `ceil(|中央値| + 3 × MAD)`（MAD は中央値からの絶対偏差の中央値で、
1.4826 は掛けない。1 ppm 未満にしない）。中心 0 の ±窓の幅で、中心をずらすなら中央値を `ppm_offset` に、`3 × MAD` を
幅に分けて渡してよい。同定が 10 件未満なら `mz_tol_ppm`（帯は `recommended_mz_tol_ppm`）は null で、`recommended_tolerance` に
`reason: too_few` と `min_n`。当てる先は `dataset_feature_search` の目標の `ppm`、`library_match_feature` の `mz_tol_ppm`。

`subset`（既定 `rank1`）: `rank1` は `rank` 位の SME すべて。`msms_matched` はそのうち、MS/MS が参照スペクトルに照合された特徴
（SML の `opt_global_spectrum_matched` が真。列は `mztab.md` §13.3。**研究室 MSP と LBM の同定を区別せず両方を数える**。偽・不明は
入らない）。`msms_matched` は `rank=1` のときだけ（それ以外は `DATASET_BAD_REQUEST`、`field: rank`。`subset` の値が不正でも
`DATASET_BAD_REQUEST`）。SML にその列が無い mzTab では誤りにせず、`subset_available: false` と `note` を付けて `overall`
`{"n": 0}`・帯の `n` 0 を返す（`subset_available` は列が無いときだけ出る）。戻り値に `subset` を書く。

`injection_trend`（選んだ `subset` の SME について）: 試料ごとの質量誤差の、注入順との関係。**全試料を数える**（ブランク・QC が
あれば含む）。試料の誤差は、その試料で検出された（gap-fill でない）セルの実測 m/z と理論 m/z の差（ppm。`ds.per_sample_mz`。
隣接 `.arf` を取り込めたときだけある）。

| キー | 意味 |
|---|---|
| `per_injection` | `[{sample, injection_order, median_ppm, n}]`（注入順の昇順。`n` はその試料で使えた特徴の数。1 つも無い試料は載せない） |
| `spearman_rho` | `injection_order` と `median_ppm` の Spearman の ρ。**`injection_order_source` が `filename`（採取時刻・注入順が読めずファイル名の順に落ちた、または試料の表が宣言 `# injection_order_source = filename` でファイル名の順と書いた）なら `null`**。データのある試料が 5 未満、または片方が定数でも `null` |
| `rho_unavailable_reason` | `spearman_rho` を出さない理由（`filename_order` / `too_few_samples`）。ρ が出ているとき、片方が定数で `null` のときは無い |
| `n_samples` | `per_injection` の件数 |
| `injection_order_source` | 注入順の出所（`manifest` / `metadata` / `filename`。§9 と同じ規則）。`filename` のとき `per_injection` はファイル名の順（シートの宣言なら、シートの値） |
| `warnings` | 注入順の注意（無ければ鍵ごと無い）。`injection_order_may_be_file_order`（`n_samples`・`message`）: 注入順が MS-DIAL の振った値（`.mddata` / mzTab-M の MS:4000089）で、ファイル名の順と一致する（MS-DIAL は注入順を与えなければファイルの順で振るので、実際の測定順とは限らない）。ρ はそのまま出す。`injection_order_source` が `filename` のときは出さない |

試料ごとの m/z が無い（`.arf` を取り込めていない）ときは `{"status": "unavailable"}` だけ。ρ が 0 に近いことは「ずれが無い」では
ない（単調でないだけのこともある。`per_injection` の `median_ppm` を見る）。

`dataset_novel_candidates(spec)`（novel-candidates.v1）:

```json
{"statistics": {"result_id": "res_..."},
 "ppm": {"default": 10.0, "bands": [{"max_mz": 200.0, "ppm": 10.0}]},
 "q_max": 0.05,
 "exclude": {"feature_ids": [], "inchikey_blocks": []}}
```

`statistics` は次のどちらか:

- `{"result_id", "contrast_id": <任意>, "merge_split": true（既定）, "alignment_qc": <JSON のパス、任意>}`: `dataset_compare` の
  2 群比較（`reference_group` が基準）。対比が複数の結果は `result_id#contrast_id` か `statistics.contrast_id` で 1 つ選ぶ
  （`{"result_id": "res_x", "contrast_id": "T_vs_C"}` は `"res_x#T_vs_C"` に正規化するので同じ `rule_hash`）。`#` の後ろと
  `contrast_id` が食い違えば `DATASET_BAD_REQUEST`、`contrast_id` が空でない文字列でなければ `NOVEL_SPEC_INVALID`。ID だけなら
  `needs_input`（`field` は `statistics.result_id`、選択肢の `value` は `result_id#contrast_id`、`argument` は無い）。
  `alignment_qc` は `dataset_alignment_qc` の `output_path`。省略時はこのハンドルで `dataset_alignment_qc` が最後に出した結果
  （§9。データ・試料の表・除外が変わっていないもの）、それも無ければその場で探す（パスは `rule_hash` に入らない）。
- `{"analysis_spec", "analysis_id"}`: **analysis-spec.v1 だけ**とその解析（v2 の spec は通らない。統計エンジンの結果は
  `dataset_compare` の `result_id` で渡す）。割れのまとめは spec の `features.merge_split` / `features.alignment_qc` が決め、
  最後の結果は使わない。

SME rank 1 の同定があり、|ppm| が帯の許容幅以内（`max_mz` 未満の最初の帯、無ければ `default`）で、q < `q_max` の特徴から除外を
外した一覧（q の昇順。`ppm` の理論 m/z と `theoretical_mz_source` は `dataset_mass_accuracy` と同じ）。`dataset_compare` の
結果は、加えて判定 `verdict` が `changed` の行だけ（`q_max` が判定の q の閾値より大きくても、`inconclusive` の行は通らない）。
analysis-spec.v1 の結果（判定を持たない）は q だけで選ぶ。行は `feature_id name inchikey adduct mz rt ppm group_ref group_test
log2fc_test_over_ref p q direction merged_from`（JSON には `sme_feature_id` も）。論文の同定リストとは突き合わせない。

**割れたスポットは 1 候補に数える。** §9 の `merge_groups` の 2 スポットは 1 つの物質で、`feature_id` は残す ID（§9）、
`merged_from` は 2 つの ID（残す ID が先頭）。統計はまとめた行列で検定し直した値（`statistics.result_id` では
`dataset_feature_search` の `merge_split` と同じ実装で、組が 1 つでもあればすべての特徴の統計をまとめた行列の検定・BH から出す。
analysis_spec では spec の `features.merge_split`）。**同定と質量誤差**（`name` `inchikey` `adduct` `mz` `ppm` `rt`）は SME（rank 1）を
持つ方のスポットから取り、その ID を `sme_feature_id` に書く（両方にあれば残す ID）。`exclude.feature_ids` は組のどちらかの ID に
当たれば組ごと外す。`merge_split=false`（または analysis_spec の `features.merge_split=false`）では `merged_from`・
`sme_feature_id`・`merge_split` を出さず、正規化した spec にも `merge_split` を載せない。

効かない設定（`merge_split=false` の `alignment_qc`、`analysis_spec` の隣に書いた `result_id` の経路の設定 `contrast_id`・
`merge_split`・`alignment_qc`）は捨てて欄ごとに `notes` に書く（`rule_hash` は書かないときと同じ）。`statistics.analysis_spec` の
効かない設定と実行時の注記（統合したデータセットの `drift_correct` の注入順。§3）も、欄の名前の前に `statistics.analysis_spec: ` を
付けて `notes` に書く。戻り値と JSON の `merge_split` は `{enabled, source, n_groups}`（`n_groups` はまとめた組の数。`source` は
`alignment_qc`（ファイル）・`last_alignment_qc`（ハンドルの最後の結果）・`find_split_pairs`（その場））。

**行の上限。** `output_path` を省けば行を TSV で返し、50 を超えれば q の昇順の先頭 50 行だけを載せる（`truncated: true`・
`n_rows` が全件・`notes`。断らない）。`truncated` は行を返すときいつも付く。全量は `output_path` の JSON。

**外部の印（`labels`、任意）。** spec に `"labels": {"path": "<TSV>", "name": "xenobiotic", "feature_column": "feature_id",
"label_column": "level", "order": ["A","B","C","D"]}` を足すと、行に `labels.<name>` の列（TSV の末尾。印が無い行は空）を足す。
`path` と `name`（英数字と `_`）は必須、ほかの既定は上のとおり。`order` は強い順で、同じ特徴に複数の行があれば最も強い水準を
取る。空の水準の行は無視、`order` に無い水準は `NOVEL_SPEC_INVALID`（`details.unknown_levels`）。割れの組（`merged_from`）は
組のどれかの最も強い水準。戻り値と JSON の `labels` は `{name, sha256, order, n_rows, n_unmatched, counts}`（`n_unmatched` は
データセットの特徴 ID に無い、印のある ID の数。無視する）で、`counts` は `{"all": {水準: 件数, "none": 件数}, "top20": {...}}`
（`top20` は q の昇順の先頭 20 行）。印は列と件数を足すだけで、並び・行数・除外は変えず、`labels` 全体（`path` を含む）を
`rule_hash` に入れない。

**質量の窓の中心のずれ（`ppm_offset`、任意）。** spec に `"ppm_offset": {"default": 0.0, "bands": [{"max_mz": 200, "offset": -8.3}]}`
を足すと、判定が `|ppm − 中心| ≤ 許容幅` になる。中心は「実測 − 理論」の ppm の符号つきの値で、帯は `ppm` の帯と同じ選び方
（理論 m/z が `max_mz` **未満**の最初の帯、無ければ `default`。ちょうど `max_mz` の特徴は次の帯）。**窓の幅（`ppm`）は変えず、
行の `ppm` 列も実測 − 理論のまま**。渡したときだけ正規化した spec と `rule_hash` に入る。形の誤りは `NOVEL_SPEC_INVALID`
（`field: ppm_offset`）。

## 9. 割れたアライメントスポット（alignment-qc.v1）

`dataset_alignment_qc(output_path=None, ppm=15.0, max_rt_gap=1.0, max_interleaved_rt_gap=0.15, min_block_size=2,
dataset_handle=None)`: 1 つの物質が 2 つ以上のアライメントスポットに割れた組を探す。**組を見つけて型を付けるだけで、行列は
変えない。** まとめは統計の側が下の `merge_groups` だけに行う（試料ごとに大きい方の値）。

- **候補**: 同じ付加体（SML の `adduct_ions`。綴りでなく組成で比べ、`[M+HCOO]-` と `[M+CHO2]-`、`[M+H]1+` と `[M+H]+` は同じ。
  付加体が分からないスポットは組にしない）・|Δm/z| ≤ `ppm`（m/z の小さい方に対して）・0.05 < |ΔRT| ≤ `max_rt_gap` 分の 2 スポット。
- **受け持ち（`owner`）**: 試料ごとに、その試料の信号を持つスポット。検出状態（`.arf` の `detected_mask`）があれば検出された方
  （両方で検出なら値の大きい方）。どちらでも検出されていない試料は、片方の値が他方の `min_dominance`（3）倍以上ならその方
  （gap-fill の値で補う）、そうでなければ `null`。検出状態が無ければ値の大きい方。どちらにも信号（値 > 0）が無い試料は `null`。
- **割れの判定**: 2 スポットがそれぞれ 1 試料以上を受け持ち、かつ補い合う。検出状態があれば、両方で検出された試料
  （`both_detected_fraction`。分母はどちらかで検出された試料）が `max_both_detected`（0.10）以下で、かつ信号か検出のある試料の
  うち受け持ちの決まった試料が `min_complementary`（0.9）以上。検出状態が無ければ（`both_detected_fraction` は `null`）、片方の値が
  他方の `min_dominance` 倍以上の試料が、信号のある試料の `min_complementary` 以上。
- **型**（`type`）:
  - `injection_block` — 片方の受け持ちが注入順で連続した塊（信号の無い試料は飛ばす）。両スポットとも `min_block_size`（2）試料
    以上を受け持つときだけ（そうでなければ `interleaved` か `group_aligned`）。群と重なっても優先する。両方が塊なら受け持ちの
    少ない方を塊とする。
  - `group_aligned` — 塊でなく、受け持ちが群と一致する（どの群でも多い方のスポットの受け持ちが `group_purity`（0.9）以上で、
    2 スポットとも少なくとも 1 つの群の多い方）。群ごとに別の異性体の可能性があるので**まとめない**。
  - `interleaved` — ほか（試料ごとに入れ替わる）。
- **注釈の両立**（`names_compatible`）: 同じ名前（大文字小文字を無視）、GOSLIN の種のレベルで同じ和組成（`PC 16:0_18:1` と
  `PC 34:1`）、または片方が注釈なし。名前は SME の rank 1、無ければ SML。加えて組の `inchikey_match`
  （`feature_metadata[fid]["inchikey"]` どうしの照合）が `full`（全体で一致）か `first_block`（先頭の 14 文字の塊＝骨格だけが一致。
  立体異性体の可能性）なら、名前が違っても両立（`names_compatible` はこの合算の真偽）。どちらかが無い・空・一致しないときは
  `none` で、名前だけで判定する。
- **受け持ちの数**: `n_owned` は `[a の数, b の数]` で、**gap-fill で値を補った試料の受け持ちを含む**。`n_owned_detected` は、
  そのうちどちらかのスポットで検出された試料だけの数（検出状態が無い run では null）。`n_owned_detected` が片側 0 なら、その側は
  gap-fill だけで支えられている。
- **`mergeable`** = `names_compatible` で、型が `injection_block`、または ΔRT ≤ `max_interleaved_rt_gap`（既定 0.15 分）の
  `interleaved`。ΔRT の大きい `interleaved` は載せるがまとめない。
- **まとめる組**（`merged`・`merge_groups`）: mergeable の組のうち**互いに最も近い相手**どうしだけ（A の mergeable な相手の中で
  |ΔRT| が最小なのが B で、かつ B のそれが A。最小が同点なら相手なし）。まとめは常に 2 スポットで、1 つのスポットは高々 1 組に
  しか入らない。ほかの mergeable な組は `merge_reason: "not_mutual_nearest"` で載せるがまとめない。
- **注入順**（`injection_order_source`）: `manifest`（pipeline の実験情報シートの `injection_order`）・`metadata`（シート以外から
  推した値、または呼び出し側が渡した順。渡した順に欠ける試料があればエラー）・`filename`（1 試料でも欠ければ試料名の順。
  `block.injections` はその順位。全試料の値をシートが宣言 `# injection_order_source = filename` でファイル名の順と書いたときも
  `filename` で、値はシートのまま）。

`output_path` の JSON（共通の約束の書き方。`pairs` は 1 組 1 行）:

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

- `dataset_fingerprint` は組を探したデータセットの指紋（`dataset_analysis_run` の `input_fingerprint.dataset_hash` と同じ関数。
  定量行列・特徴 ID・試料・検出状態・出所から作る。同じ run・同じハンドルなら同じ値）。`dataset_analysis_run` の
  `features.alignment_qc` は、これが今のデータセットと違う JSON（別の極性・別の run）か、この欄の無い JSON を
  `ANALYSIS_INPUT_INVALID`（`field: features.alignment_qc`）で拒む。
- `pairs[].features` は `injection_block` なら [塊でない側, 塊の側]、ほかは RT の昇順。`delta_rt = rt[1] − rt[0]`
  （`injection_block` では塊の RT のずれ。符号つき）、`delta_ppm = (mz[1] − mz[0]) / mz[0] × 10⁶`。`block.injections` は塊の
  最初と最後の注入順。
- `groups[]` は**すべての組**の連結成分（報告用。3 スポット以上の連鎖を 1 つにする）。`features` は特徴 ID の小さい順（すべて
  整数なら数の順）。`type` は成分の型で、組の型が混ざれば `group_aligned` > `interleaved` > `injection_block` の順で強い方。
  `mergeable` は成分全体が 1 つのまとめる組のときだけ真。
- `merge_groups[]` はまとめる組ごとの 2 スポット（`merge_group_id`・`group_id`＝含まれる `groups[]` の成分・`features`・
  `pair_ids`・`type`・`types`・`mergeable: true`・`keep`・`keep_basis`・`n_owned`）。残すのは**受け持ちの多い側**、同数なら
  小さい ID（すべて整数なら数の順）で、`keep` がその ID、`keep_basis` が決め手（`n_owned` / `smaller_id_tie`）、`n_owned` が
  `{特徴 ID: 受け持ちの試料の数}`。**組の並びは残す ID が先頭**（`features`・`merged_from` の値。`groups[]` は ID の小さい順の
  まま）で、`merge_matrix_columns` は渡した並びで行列にある最初の ID を残す。`keep` の無い JSON の組は小さい ID を残す。統計の
  側（`dataset_analysis_run` など）がまとめるのは `merge_groups` だけ。`summary.n_merged_pairs` はその数。
- `split_partners(qc)`（統計の側が使う）: 特徴 → `{group_id, merge_group_id, partners, component, type, mergeable}`。まとめる
  特徴は `partners` がその相手 1 つで `mergeable` が真、`type` はその組の型。まとめない特徴は `partners` が空・`mergeable` が偽・
  `merge_group_id` が null。`component` は全体の成分のほかの特徴（報告用）。
- 戻り値は `status`・`schema`・`injection_order_source`・`signal`（`detected` / `values`）・`summary`・`kept_on_handle: true`・
  `output_path`（渡したとき）と、`group_aligned` の組があれば `n_group_aligned`（全件）・`group_aligned_truncated`・先頭 50 組の
  TSV `group_aligned`（`pair_id feature_a feature_b mz rt_a rt_b delta_rt adduct name_a name_b`。超えれば `notes` に、全組は
  `output_path` の JSON の `pairs` にあると書く）。
- **最後の結果はハンドルに持つ**（`kept_on_handle`）。`output_path` と同じ中身（同じバイト列の sha256）を、そのハンドル
  （データセット）に最後の 1 つとして持ち、前の結果は置き換える。`dataset_feature_search` の `alignment_qc_path`・
  `dataset_novel_candidates` の `statistics.alignment_qc` を省いたときは、まずこれを使う
  （`merge_split.source: last_alignment_qc`。`alignment_qc_sha256` はファイルで渡したときと同じ値）。データの指紋（`dataset_fingerprint`）が作ったとき
  と違うとき、試料の表・除外が変わったとき（`metadata_revision`。`dataset_exclude`・`dataset_set_sample_metadata`）は使わず、
  その場で探す。`dataset_analysis_run`（§2）と `dataset_claim_evidence`（§7）はいつも使わない（spec だけで再現する）。束の割れの
  組は既定の閾値の `find_split_pairs` なので、既定でない閾値で `dataset_alignment_qc` を回した後は、同じ目標でも
  `dataset_feature_search` と組・`counts`・`status`・`rule_hash` が違いうる。

## 10. 根拠バンドル v3（evidence-v3）

`dataset_export_evidence_bundle` が書く。v2 に試料表・試料ごとの値と検出・付加体・極性を足したもの。`result_id` を省けば
直近の差次的結果。対比は `contrast_id` か `result_id#contrast_id` で選び、2 群の対比が複数で選ばなければ `needs_input`
（選択肢は対比の一覧）。交互作用の対比は書けない（`CONTRAST_NOT_TWO_GROUP`）。ARF 経路（試料ごとの値を取り出さない）は v2 の
まま。契約の正は massbank-context の `docs/contracts/evidence-bundle-v3.md`。TSV は全ファイルとも同じ書き方（タブ区切り・
UTF-8・改行は CRLF）。

**特徴 ID はスポットの ID のまま。** 割れの組（§9）へまとめた値は書かない。まとめた組への読み替えは呼び出し側が
`merged_from` で行う。

`manifest.json` の v2 との差分:

| 欄 | 意味 |
|---|---|
| `contract_version` / `exporter.version` | `3` / `evidence-v3` |
| `comparison.polarity` | データセットの極性（`positive` / `negative`）。ハンドルなら `load_dataset` の `polarity`、対話のセッションなら同じ規則で決める。決まらなければ `null` |
| `comparison.preprocess` | `value` に当たった前処理の設定（結果が読んだ行列の前処理の設定） |
| `comparison.method` | 実際に使った検定名（`moderated_t` / `welch` / `student` / `mannwhitney`） |
| `comparison.statistics` | `{test, ci_level, verdict_rule: {q_threshold, equivalence_fc, equivalence_log2}, bh_population: {name, n, q_interpretation}, subject_source, n_units_a, n_units_b}`。`n` は選んだ対比の BH の母集団の大きさ（比較の補正の母集団と同じ）。`q_interpretation` は `features.tsv` の `q_value`（BH）の読み方（`discovery_fdr`、母集団が `reported` なら `reproduction_only`。`dataset.md` の書き出しのメタ行と同じ） |
| `comparison.effect_size_definition` / `mean_definition` | 効果量の定義（`{name, description}`）と、`mean_a` / `mean_b` が「個体ごとの前処理後の線形値の算術平均」であること |
| `comparison.matrix_id` / `rule_hash` / `contrast_id` | どの行列・どの規則・どの対比を書いたか（`contrast_id` は 2 群の結果なら 1 つの対比の ID） |
| `comparison.correction` | `features.tsv` の `q_value` を作った多重比較の補正で、いつも `BH`（`significant` もその q で読む） |
| `comparison.p_adjust_method` / `q_value_source` | 比較の補正が BH でない（`dataset_compare(p_adjust=bonferroni / holm / none)`）ときだけ。比較の方法（`bonferroni` / `holm` / `none`）と `recomputed_from_p_value`（同じ母集団の p から BH の q を計算し直し、判定と `significant` もその q で読み直した）。戻り値にも同じ 2 つと `notes` |
| `exploratory` | 最上位の欄。探索専用のデータセット（未完了の実行を `allow_incomplete=True` で開いた・検証を通っていない pipeline の単位）から書いたときだけ `{exploratory_only: true, reason}`（書き出しは断らない）。同じものを `comparison.exploratory` にも写す。探索専用でなければどちらの欄も無い |
| `comparison.unadjusted_confounded` | `allow_confounded=true`（既定）で完全交絡（全対比をまとめた表か、その対比の表）を越えて走った比較だけ `true`（それ以外は欄が無い）。同じ印は来歴 `provenance.result.comparison.unadjusted_confounded`（と `provenance.result.warnings` の説明）にも載る |
| `intensities_stage` | `pre_imputation`（補完の前）/ `post_imputation`（補完の後） |
| `population.detection_unavailable_reason` | `detected` が全行空のときの理由。検出状態を使えれば `null` |
| `population.exclusions` | 除外の来歴。エンジンの結果は、その行列を作ったときの除外（`dataset_exclude` の試料・特徴。名前は 20 件で切り、残りは `n_more`）を、結果が読んだ行列の登録から引く（今の `ds.excluded_*` ではない） |
| `files[].role` | `samples` と `intensities` が増える |

`samples.tsv`（データセットの全試料）:

| 列 | 意味 |
|---|---|
| `sample_id` | 行列の試料名 |
| `group` | 比較の群 A の試料は `a`、群 B は `b`、比較に入らない試料（QC・ブランク・群なし・`include=false` で除いた試料）は空。`intensities.tsv` には比較の試料だけ |
| `role` | 実験情報シートの `role`（除外された試料も `sample` のことがある）。シートが無ければ名前の `_` 区切りのトークンからの推定で、`QC1` のようにトークンにならない名前は `sample` になる（`group` が空でも QC とは限らない） |
| `injection_order` | 整数。分からなければ空 |

`intensities.tsv`（縦持ち。全特徴 × 比較に入る試料。QC・ブランクの行は無い）:

| 列 | 意味 |
|---|---|
| `feature_id` | `features.tsv` の行と同じ |
| `sample_id` | `samples.tsv` の試料名 |
| `value` | log 変換の前の値（有効数字 10 桁）。欠測は空。結果の試料が値の行列に無いときも空。`0` は実際の値（補完は NaN だけを触る。実データでは多く、`detected=false`）で、欠測の空とは別。そのため `pre_imputation` と `post_imputation` の値は同じになりうる。検出の有無は `value` でなく `detected` で読む |
| `detected` | `true` / `false` / 空。`false` は gap-fill か補完で、実測のピークではない。**空は「不明」で、gap-fill の意味ではない** |

`intensities_stage` の決め方:

- エンジンの結果（`dataset_compare`）: 常に `pre_imputation`。値は結果が読んだ行列（`ds.stat_matrices[matrix_id]`）の
  `sample_values`（**正規化の後・log2 の前の線形値**。欠測は空。エンジンは検定に補完しない）。前処理のレシピを当て直さない。
  行列が登録に無い（除外・試料の表の変更で古くなった）ときは `missing_state`（`dataset_preprocess` と `dataset_compare` を
  やり直す）、結果の行と行列の特徴の集合が違えば `EVIDENCE_INTENSITIES_MISMATCH`。

`detection_unavailable_reason` の値:

| 値 | 意味 |
|---|---|
| `dataset_status` の `detection` と同じ理由（`sample_axis_mismatch`・`mz_mismatch` など） | 検出状態を取り込めなかった。取り込みのときの理由 |
| `not_imported` | 検出状態が無く、理由も記録されていない |
| `mask_shape_mismatch` | 検出状態の形が特徴 × 試料と合わない |
| `sample_axis_mismatch` | 結果の試料が `ds.sample_names` に無い（検出状態の試料軸とずれている） |

`features.tsv` に任意の列 `adduct` が増える（mzTab-M の SML の `adduct_ions`。SML が無い特徴は rank 1 の SME の付加体。無ければ
空）。`features.tsv` の `significant` は差次的の書き出しの `significant` 列と同じ印で、エンジンの結果は判定 `verdict` が
`changed` の行だけ `true` になりうる（`dataset.md` §11.4。判定できない行は q が小さくても `false`）。`q_value` と `significant` の
補正は上の `comparison.correction`・`p_adjust_method` / `q_value_source`。

`annotations.tsv` の `source` が `msdial_sme` の行の `provenance_json` は SME の候補そのもの（mzTab-M の SME の欄）。その
`theoretical_mass_to_charge` は**ライブラリの記録のプリカーサー m/z**（mzTab-M の SME の値。job / pipeline を解析フォルダで
読んだときは照合した mzTab-M から写した値）だけを表す。`load_dataset(directory=...)` で開いたデータセットが組成式と付加体から
計算した理論 m/z（§8 の `formula_adduct`）は**写さず `null`**（欄も足さない）。極性を統合したデータセットでも、計算した側の
特徴だけが `null`。計算した値と出所は `dataset_mass_accuracy` / `dataset_novel_candidates` で読む。

戻り値には v2 の欄に `n_samples`・`n_intensity_rows`・`intensities_stage` が足される（探索専用のデータセットなら `exploratory`
も。比較の補正が BH でなければ `p_adjust_method`・`q_value_source`・`notes` も）。

### `spectrum_links.tsv` と `.dcl` の出どころ

dataset 経路も ARF 経路と同じ規則で、**同じアラインメントの `.dcl`** の代表スペクトルに結ぶ。結ぶのは ID だけ
（**SMF_ID == `.arf2` の AlignmentID == `.dcl` の番号**）で、m/z・RT は結びの検証にだけ使う（m/z の近さでは結ばない）。結べた
特徴は 1 行（`spectrum_file` は `.dcl` のファイル名、`spectrum_id` は `<.dcl の語幹>:<番号>`、`mapping_evidence` は
`same-alignment-stem; SMF_ID == AlignmentID == DCL index; precursor/RT verified`）。結べた行が 1 件でもあれば `.dcl` をバンドルへ
写し、`files[]` に `role: spectra` で載る。

使う `.arf2` / `.dcl`: ハンドルなら `load_dataset` が決めたもの（§1）。対話のセッションなら同じ探し方（検出状態を取り込んだ
`.arf` の同じバッチ → mzTab-M の隣）。

検証は 2 段:

1. **番号の空間**: SMF_ID の m/z・RT（mzTab-M）と、同じ番号の `.arf2` のスポットの m/z・RT を比べる。1 件でも合わなければ全体を
   結ばない（番号の付け方が違う別のアラインメントとみなす）。
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
