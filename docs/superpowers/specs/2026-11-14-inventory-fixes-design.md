# 仕分けの一覧（2026-11-13）の指摘を直す設計

- 日付: 2026-11-14
- 出どころ: reanalysis-study の仕分けの一覧（2026-11-13）のうち ms-data-parser に向いた 74 件と、その事実確認
  （HISTRY 2026-11-13、全文は追跡外の `data/factcheck/2026-11-13-inventory-factcheck.md`）で新しく見つかった欠陥
- 利用者の指示（2026-11-14）: 指摘をすべて直す。`MCP_INSTRUCTIONS` は無駄を削ぎ落として最適化する。`save_figure` は
  推奨の直し（`contrast_id` の引数を足す）。Claude Desktop のサーバは `.venv` で動かす
- 方針（利用者、2026-11-13）: ツールは LLM の行動を**コードで縛らない**。手順・作法は自然言語のルール（研究側の文書）に置き、
  ツールは事実と、選択肢の材料と、記録を返す。解析的に誤る点と足りない機能は直す
- 例外として残す縛り: Console を回す前に利用者の許可を取る（ENTRY POINT。2026-10-07 の利用者判断）。研究室の参照ライブラリを
  外へ出さない（パスを出さない）。massbank-context との契約（差次的エクスポートの 15 列 `CONTRACT_VERSION=1`、
  根拠バンドル `BUNDLE_CONTRACT_VERSION=3`）は**列を変えない**。足すのはメタ行と、契約外の出力（`output_path` の全列 TSV など）だけ

中心の主張が誤りだった 4 件（A9 組成・B4 クラスの和と比・D5 割れの組・B11-T46 handoff の件数）は直す対象から外す。
ただし事実確認で見つかった、それぞれの小さな残り（文書の食い違い・1 行の注記など）は直す。

## 担当（worktree 1 本 = 担当 1 つ）と持ち分

担当は**ファイルの持ち分**で分ける。持ち分の外のファイルは、ここに書いた最小の変更だけをする。
`docs/` と `USAGE.md` は各担当が自分のツールの節だけを直す（統合のとき main 側で衝突を解く）。
ツール数が変わる担当（3・8）は、自分のツールだけを足した数で `tests/test_server_registration.py`・`USAGE.md`・
CLAUDE.md 冒頭を直す（統合で合算する）。`MCP_INSTRUCTIONS` は統合の後に main で書き直す（担当は触らない）。

### 担当 1 `fix/inv-msms-library` — MS/MS の合否と参照ライブラリ

持ち分: `analysis/msms_rules.py` `analysis/spectral_match.py` `library/*` `pathway/msp_lipid.py`
`analysis/claim_evidence.py`（MS/MS の部分）`tools/reanalysis_tools.py` の `library_msms_qualify`。

- **A1-T14**: 断片の一致に測定側の強度の下限を足す。規則のキー `min_fragment_relative_intensity`（測定の基準ピークに対する比、
  既定 0.01）。断片の一致数（前駆体を除く）はこの下限以上の測定強度を持つ窓だけを数える。`matched_peaks_count`（MS-DIAL の写し）の
  意味は変えない。既定の `min_fragments_excl_precursor` を `null` → `1` にする（一致 0 本で qualified にならない）。既定が変わるので
  `rule_hash` は変わる（意図どおり。固定値のテストと reanalysis.md §4 の例を直す）。`.dbs` の `search_params` の足切り
  （dot・一致率・総合スコア・minimum_spectrum_match）は合否には当てず、候補ごとに「`.dbs` の足切りを通るか」の欄で返す（材料）。
- **A2-T13**: 脂質の種の目標を、pygoslin で読めない名前でも引けるようにする。記録の name の注記（付加体・`|` の後など）を外して
  読み直し、それでも読めなければ Ontology と分子式（`msp_lipid.from_class_formula`）から種を作る。`CLASS_FORMS` に LPS・LPG を足す。
  どの経路で種を得たかを行に出す。reanalysis.md の「記録の name の種で決める」規則も書き直す。
- **A3-T12**: 前駆体の既定の窓を MS-DIAL と同じ `max(0.01 Da, m/z × 20 ppm)` にする（m/z 500 以下は今と同じ）。候補の表に
  `delta_ppm`（= (測定 − 参照) / 参照 × 1e6。`library_match_feature` の同名の列と同じ名前・符号。列名は `ppm_error` では
  ない）を出す。`library_msms_qualify` と `dataset_claim_evidence` も `mz_tol_ppm` を受け、使った窓を返り値に出す。
  library.md:380 の「窓は MS-DIAL と同じ」を正しくする。
- **E6-T25**: 照合の上位が MS-DIAL の注釈と食い違うとき、`library_match_feature`・`library_msms_qualify`・
  `dataset_claim_evidence` が印（`top_differs_from_msdial` と、MS-DIAL の注釈の候補順位）を返す。claim の合否で目標が
  top1 でないとき、目標の順位を返す。
- **E7-T26**: `min_relative_intensity` を中性損失の要素でも受ける（意味: 生成イオン＝前駆体−損失の測定強度の比）。既定値は
  足さない（省いた規則の `rule_hash` を変えない）。エラー文は原因のキーを名指す。reanalysis.md:346,355-356 を実装に合わせる。
- **B6-G04**: `library_match_feature` の既定の候補表に `inchikey` を出し、標準の InChIKey か（第 2 ブロックの 9 文字目が `S`）を
  `inchikey_standard` で示す。外部への照会（PubChem）は足さない（研究室のライブラリの名前を外へ送らない）。
- **C16-M18**: `library_plot_mirror` が `spot_id`（と任意の候補の指定）を受け、自分で照合し直せる。照合に失敗したら直前の照合を
  捨てる（古い照合を黙って描かない）。図の返り値にどの spot・候補かを出す。
- **A11-V68（library 側）**: `ion_mode` を `core/polarity_names.normalize_polarity` で正規化し、読めない値は誤りにする（黙って誤らない）。
- **C5-M08（library 側）**: `library_load` の `library_names` を外す。返り値のライブラリはいつもファイル名と sha256（パスは出さない）。

### 担当 2 `fix/inv-relations-search` — 特徴の関係と特徴検索

持ち分: `curation/relations.py` `analysis/feature_search.py` `analysis/mass_accuracy.py`
`tools/dataset_analysis_tools.py` の `dataset_feature_search`、`tools/reanalysis_tools.py` の `dataset_feature_relations` と
`dataset_mass_accuracy`。

- **A4-T15**: 付加体・源内断片の関係でも、Δm/z に加えて同時溶出（RT）と試料間の強度の相関 r が閾値以上なら `strong` を立てる。
  既定の解析フォルダの読み込みでも同位体の `strong` が立つようにする（MS-DIAL のリンクが無くても、Δm/z・RT・r と期待比で）。
  curation.md:324-327 の規則は「Δm/z だけでは断定しない。r と RT がそろえば strong」と書き直す。`curation_suggest` への影響を
  テストで確かめる。reanalysis.md:384-385 の「裏付け＝strong」を直す。
- **A5-T16**: 37Cl（+1.99705）と 81Br（+1.99795）の M+2 を同位体の種類として足す。13C2（+2.00671）と窓が重なるときは
  一方に決めず `ambiguous` の印と両方の候補を返す（13C2 の期待比に合わないことを「別物質」と言わない）。curation.md:312 を直す。
- **A12-T18**: 陽イオンの既定の付加体の層に `[M+CH3OH+H]+` と `[M+H-2H2O]+` を足す（`rule_hash` は変わる。テストと spec の
  既定の表を直す）。
- **B5-V67**: 特徴検索の統計に Fisher の p（`presence_p`）と `detection_driven` の印を足す。`as_two_group` が記録の
  `presence_p` と個体単位の検出数を 2 群の形に引き継ぎ `detection_driven` の印を足す部分も**この担当で入れた**（当初は
  担当 7 の持ち分だったが、列にする側と同じ変更にまとめた）。列の並びを固定したテストを直す。
- **B1-T17**: `dataset_feature_search` に `include_statistics`（既定 True）を足し、False なら比較の統計を付けない。
- **D3-rt-offset**: 特徴検索が目標の RT と特徴の一意の一致から系統的な RT のずれを推定して返す（一意の一致 n ≥ 10 かつ
  |中央値| > 0.2 分のときだけ `shift_min` を出す。適用はしない）。目標の RT の中心をずらす `rt_offset_min` を受ける
  （省けば `rule_hash` は変わらない）。
- **A3（材料）**: `dataset_mass_accuracy` が、測った質量誤差から勧める `mz_tol_ppm`（|中央値| + 3×MAD を切り上げ）を返す。
- **C7（特徴検索の部分）**: 比較の結果が無いと効かない引数を黙って受けない。効かなかった引数を `notes` に書く（誤りにはしない）。

### 担当 3 `fix/inv-stats-preprocess` — 前処理・試料の表・取得時刻

持ち分: `analysis/stats/preprocess.py` `analysis/stats/data.py` `analysis/sample_manifest.py`
`analysis/stats/confounding.py`（注入順の出所の警告だけ）、新しい `metabolomix/rawmeta/`、`tools/dataset.py`（`list_data_files`）、
`tools/stats_tools.py` の `dataset_preprocess`、新しい `tools/sample_sheet_tools.py`。

- **A6-T47**: 注入順の出所を試料ごとに持つ（`sheet` / `mzTab`（MS:4000089）/ `mddata` / `raw_acquisition_time` / `filename`）。
  MS-DIAL 由来の値がファイル名の順と一致するときは「MS-DIAL の既定（ファイル順）かもしれない」の警告
  `injection_order_may_be_file_order` を、交絡の判定・ドリフト補正・ρ の返り値に出す。`drift_mode="auto"` は、注入順の出所が
  確かでない（ファイル順と一致する MS-DIAL 由来・ファイル名）とき補正しない（理由を返す。`on` を明示すれば補正する）。
  シートの空欄を `.mddata` で埋めた値も、試料の表に値と出所で出す。
- **A7-V65**: 試料名の 8 桁からのバッチの推定は、黙って当てない。推定は候補（`batch_candidates`、出所付き）として返すだけにし、
  ドリフト補正・反復のまとめ・交絡の判定には、シート・宣言・明示のバッチだけを使う。PQN の挙動は変えない（QC が無いときの
  参照を返り値で明示する）。テストと dataset.md:184・248 を直す。
- **B2-T21（正規化の部分）**: 正規化の係数を特徴の部分集合から計算できるようにする（`normalization_features` /
  `normalization_exclude_features`。外した特徴も解析には残る）。
- **C3-M25**: 換算済みで QC があるときの `drift_mode="auto"` は `needs_input` で尋ねず、`off` で計算して、選んだ理由と
  選択肢（`on` で呼び直せること）を返す。補正できないデータでは尋ねない。テストを直す。
- **E10-T21b**: ドリフト補正の LOESS に局所 2 次（`loess_degree`: 1 か 2、既定 1）を足す（2 次は numpy で書く）。
- **D1-acquisition**: reanalysis-study の `acquisition.py` を移す（.wiff の OLE の時刻・Thermo .raw の先頭のファイルヘッダ・
  Bruker .d の tdf/sqlite・mzML の startTimeStamp。**取得時刻のメタデータだけを読み、スペクトルは解さない**）。
  `list_data_files(all_files=True)` が raw ごとの取得時刻と、そこからの順を返す。群が分かれば（シートがあれば）群と取得順の
  重なりの判定と 40 分を超える間隔の数も返す。CLAUDE.md の「生の測定ファイルは自分では解さない」を「スペクトルは解さない。
  取得時刻のメタデータだけ読む」に直す。
- **D2-sample-table**: 新しいツール `pipeline_sample_sheet`。ISA の試料の表（`s_*.txt` と `a_*.txt`）か mwTab
  （`SUBJECT_SAMPLE_FACTORS`）から、因子の列付きの実験情報シート（sample-manifest v2）を raw のフォルダに書く。群は指定した因子の
  列の値の組。注入順は raw の取得時刻の順（1 つでも読めなければ空にし、出所を返す）。

### 担当 4 `fix/inv-stats-compare` — 比較・検定・多変量

持ち分: `analysis/stats/compare.py` `limma.py` `design.py` `spec_v2.py` `confounding.py`（関門の部分）、新しい
`analysis/stats/plsda.py`、`tools/stats_tools.py` の `dataset_compare` と `dataset_pca` と行列の解決（`matrix_inputs`）。

- **B2-T21（補正と多変量）**: 多重比較の補正に `bonferroni`・`holm` を足す（`p_adjust`: `bh` 既定 / `bonferroni` / `holm` /
  `none`。compare と spec v2）。記録に `p_adjust_method` を書く（担当 7 が書き出しのメタ行と根拠バンドルの `correction` に使う）。
  PLS-DA と VIP を spec v2 の解析の種類 `plsda` として足す（成分数・交差検証の Q2・並べ替え検定の p を必ず返す。ツールは増やさない）。
- **B3-interaction / E1-T07**: 要因モデルの主効果と交互作用の F（多水準を含む moderated F）を返す。向きと log2fc を持たないので、
  比較の行とは別の欄（omnibus と同じ扱い）で返す。
- **C4-M23**: 完全交絡でも比較を止めない（`allow_confounded` の既定を True）。交絡は返り値の先頭の `confounding` 欄と、
  対比ごとの `confounded` の印で示す。テストを直す。
- **C12-M13**: `matrix_id` を省いたときは、そのハンドルで `dataset_preprocess` が最後に作った行列を使う（無ければ既定の前処理）。
  返り値に使った `matrix_id` と出所（`latest_preprocess` / `default`）を出す。
- **B6-G01 / G02 / G03**: `dataset_compare` に `output_path` を足し、全特徴 × 全対比の全列（統計量・df・CI・SE・群ごとの n・
  検出数・同定の列（名前・ontology・注釈の出所・`msms_assigned`・道具の msi_level））を TSV で書く。行に `se` を足す（定義できない
  検定では空）。top の行に名前を足す。差次的エクスポートの 15 列は変えない。
- **B6-G07**: `dataset_compare` に `reported_feature_ids` を足し、`bh_population="reported"` を受ける。記録の spec にも書く。
- **D4（片群一定の検査）**: 一方の群の値がすべて同じ（埋め値）の特徴に `constant_in_group` の印を付け、件数を警告で返す。

### 担当 5 `fix/inv-pipeline` — pipeline と Console

持ち分: `pipeline/*` `console/*` `core/process_control.py` `core/detach_relay.py` `tools/pipeline_tools.py`
`tools/console_tools.py`。

- **B7-T63**: 受付と prepare_input で mzML/mzXML の末尾（閉じタグ）とサイズの下限を確かめ、壊れたファイルを名指しして
  断る（`PIPELINE_INPUT_INVALID`）。Console が落ちたら `msdial.log` から原因のファイルを探して返り値に出す。
- **B8-T66**: Console を起こす前に commit の空き（`GlobalMemoryStatusEx` の `ullAvailPageFile`）を測り、足りなければ警告を
  記録に残す（止めない。閾値は設定 `[msdial] min_free_commit_gb`、既定 8）。run を跨いで Console を同時に 1 本にする
  機械全体のロック（OS のファイルロック。プロセスが死ねば外れる）を入れ、待っている間は工程の状態にそう出す。
  終了コード 1450（`ERROR_NO_SYSTEM_RESOURCES`）を資源不足と分類して返す。stdlib だけの leaf に置く。
- **E2-T02**: worker を起こした後、`IsProcessInJob(worker, NULL)` で worker がまだ Job に入っているかを確かめる。入っていて
  kill-on-close かどうか分からない Job なら、run の記録と返り値に警告を出す（何段目の Job でも捕まえる）。
- **C21-M27**: 明示の要求・単位・resume の updates の下流キーは、専用のコード `PIPELINE_DOWNSTREAM_REMOVED` をやめ、未知キーの
  一般の誤りに寄せる（文面で代わりのツールを案内する）。`analysis-request.json` の読み飛ばし（警告）は残す。
- **A11-V68（Console 側）**: `console_method_candidates` の極性を `normalize_polarity` で正規化し、読めない値は誤りにする。

### 担当 6 `fix/inv-reanalysis-tools` — spec 駆動の道具と QC

持ち分: `tools/reanalysis_tools.py`（担当 1・2・8 が触る関数を除く）`analysis/analysis_spec.py` `analysis/alignment_qc.py`
`analysis/homolog.py` と、新しい外部の行列の読み取り。

- **C7-M30/M29**: 行数の上限で断らない。先頭 N 行と `truncated`・全件数を返し、`output_path` があれば全量をそこへ書く
  （`dataset_analysis_run` の `ANALYSIS_OUTPUT_REQUIRED`、`dataset_novel_candidates`、`dataset_sample_peaks` の 1,000 行）。
  効かない引数の組は誤りにせず、`notes` に書く。`dataset_alignment_qc` が 50 組を超えて外したものは件数を返す。
- **C18-M20**: `dataset_alignment_qc` の最後の結果をハンドルに持ち（試料の表・除外が変われば捨てる）、`dataset_feature_search`・
  `dataset_novel_candidates` が `alignment_qc_path` を省いたときに使う（どちらを使ったかを返す）。`dataset_analysis_run` は
  spec だけで再現する決まりなので使わない（文書に書く）。
- **B12-V77/T72**: `dataset_alignment_qc` の `output_path` の JSON を改行入りで書く（Read で区切れるように）。
- **E9-T33**: v1 の `aggregate` と `features.feature_ids` の併用は誤りのまま、説明文・reanalysis.md・USAGE.md に書き、
  テストで固定する。
- **B4（残り）**: v1 の `aggregate` の返り値にメンバーごとの寄与（和に占める割合）を足す。
- **D4-external-matrix**: v1 の `input.source` に `maf`（MetaboLights）・`mwtab`・`xlsx`（標準ライブラリで読む）を足す
  （reanalysis-study の `paper_values.py` と `xlsx.py` を移す）。外部の値との Spearman 相関を v1 spec の `correlate` で返す
  （特徴の対応は spec で渡す）。

### 担当 7 `fix/inv-results-consumers` — 結果の消費側（書き出し・図・記録）

持ち分: `analysis/stats_results.py` `analysis/dataset_export.py` `analysis/evidence_export.py` `analysis/export_contract.py`
（メタ行だけ）`tools/dataset_analysis_tools.py` の書き出し・根拠バンドル・試料情報の関数、`tools/plot_tools.py` `plots/*`
`tools/reports.py`。

- **B6-G08（利用者が推奨を選んだ）**: `save_figure` に `contrast_id` の引数を足す。`result_id#contrast_id` と食い違えば
  `DATASET_BAD_REQUEST`。テスト・reports.py の説明・USAGE.md を直す。
- **C17-M19**: `dataset_plot` に `output_path` を足して、描いた図をそのまま保存できるようにする。説明文の「明示的なユーザー要求時
  だけ」の縛りを外す。図が古いときは `missing_state` でなく理由の分かるコード（`STALE_FIGURE`）を返す。
- **C13-M14**: `dataset_export_evidence_bundle` の `result_id` を省けるようにする（直近の結果）。
- **C14-M15**: すべての対比を選ぶツールが `contrast_id` の引数と `result_id#contrast_id` の両方を受ける。省いて対比が複数なら、
  誤りでなく `needs_input`（対比の一覧付き）。1 つなら黙ってそれを使う。
- **C19-M24**: 未完了の run のデータセットでも書き出しを断らず、メタ行（と根拠バンドルの manifest の欄）に `exploratory` と理由を
  書く。印はジョブの status だけでなく pipeline の単位の状態も見る。docstring の「dataset_compare も拒む」を直す。
- **D6-drop-features**: `dataset_export_differential` に `exclude_feature_ids` を足す（検定の後に行を落とす。q は計算し直さない
  ことと落とした件数をメタ行に書く）。
- **E1-T07（書き出し）**: 2×2 の交互作用の対比（差の差。向きと log2fc を持つ）を書き出せるようにする。
- **B2（消費側）**: 書き出しのメタ行と根拠バンドルの `correction` に、記録の `p_adjust_method`（無ければ BH）を出す。
- ~~**B5（消費側）**~~: `as_two_group` の `presence_p`・個体単位の検出数・`detection_driven` は担当 2 へ移した（担当 2 の
  B5-V67）。
- **A11-V68（図）**: volcano の題の「vs」の並びを `test vs reference` にそろえる。
- **B12（図の大きさ）**: `dataset_plot` と（担当 1 の）mirror の既定の PNG の大きさを見直す（読める範囲で小さく）。
- **C11-M06**: `write_report` の説明文から推奨の見出しの一覧を外す（見出しは自由）。

### 担当 8 `fix/inv-dataset-loading` — 読み込み・ハンドル・状態・長い呼び出し

持ち分: `dataset/*` `mztab/*` `tools/handle_tools.py` `tools/mztab_tools.py`（`dataset_status`）`tools/overview_tools.py`
`tools/objective.py` `core/pending_loads.py` `core/handle_dispatch.py` `core/handles.py` `core/session_state.py`、新しいツール
`job_status`。

- **C1-M02/M03・A9 の残り**: `assay_kind` をハンドル（データセット）に持たせる。`record_objective(assay_kind=...)` は直近の
  （または渡した）ハンドルに書き、目的ファイルにも残し、`load_dataset` は解析フォルダの目的ファイルにあれば読み戻す。
  spec v2 と概観はハンドルの `assay_kind` を使い、返り値に出す。グローバルの `session.assay_kind` はダイジェストの切り替えのため
  直近のハンドルの値を読む property にする。`pathway_project` の lipid の関門は担当 9 が外す。
- **C2-M07**: `load_dataset` の `overview` の既定を False にする（概観は `dataset_overview` で出す）。テスト・docstring・USAGE・
  dataset.md・CLAUDE.md の地図を直す。
- **C5-M08（読み込み側）**: `load_dataset` の `library_names` を外す。ライブラリはいつもファイル名と sha256。
- **C6-M09/M10**: `wait_seconds` は 0〜30 の外を断らず丸めて、丸めたことを返す。`dataset_load_status` の `load_id` を省いたら、
  受付中の読み込みが 1 つならそれ、複数なら `needs_input`（一覧付き）、無ければ `missing_state`。
- **B10-T27/T43**: 長い呼び出しの仕組みを `handle_dispatch` に一般化する（`load_dataset` の受付と同じ形）。`pathway_project` と
  `dataset_homolog_check` に `wait_seconds`（既定 20）を足し、超えたら `status: "running"` と `job_id` を返して続ける（ハンドルの
  ロックは仕事が終わるまでその仕事が持つ）。新しいツール `job_status(job_id, wait_seconds)` で結果を受け取る。
- **A10-T23**: `dataset_status` の summary の `samples` TSV を試料の表の役割から作る（見出し `name\trole` は変えない）。
- **B6-G05**: `dataset_status(view="results")` で、ハンドルの行列（`matrix_id`・作り方の要約）と結果（`result_id`・kind・対比・
  古いか）を一覧で返す。
- **B6-G09**: mzTab-M のアダプタが隣の `.arf` と照合できたとき、試料ごとの RT（`per_sample_rt`）を捨てずに持つ。
  `dataset_sample_peaks` が断らなくなる。文書の「mzTab-M は持たない」を直す。
- **B9-T06**: Console の `.mddata`（class=FileID・type=0 の既定）を群・役割として信じない。解析フォルダが pipeline の run の中に
  あれば、その run の実験情報シートと受領から群と役割を当てる（`directory=` で開いても）。ライブラリも run の記録から引く。
- **E5-T24**: 名前から InChIKey を引いたとき同名衝突があった特徴に、特徴単位の曖昧の印を付ける（書き出しの列は変えない。
  `dataset_feature_search` などの同定の欄と mztab.md に出す）。
- **E8-T31**: ハンドルを見失う条件を reanalysis.md に書く。`dataset_close` の返り値に閉じたハンドルの出所を出す。
- **C10-M05**: `record_objective` が目的ファイルを作るとき、渡した欄だけで作れる（欠けた欄は空のまま。後で足せる）。
  Q の採番は残す。

### 担当 9 `fix/inv-pathway-curation-quant` — パスウェイ・キュレーション・定量

持ち分: `pathway/*`（`msp_lipid.py` は担当 1）`tools/curation_tools.py` `curation/*`（`relations.py` は担当 2）`quant/*`。

- **A8-V72/V60**: `pathway_project` が `dataset_compare` と同じ検定の選択を `method`（既定 `moderated_t`）で受ける
  （引数名は `test` ではなく `method`）。「常に反復をまとめる」と書いた docstring・
  pathway.md:71 を、行列の units に従うと直す。
- **C1（パスウェイ側）**: `pathway_project` の lipid の関門（`missing_state`）を外す。脂質のアッセイでないときは、名前の文法が
  前提であることを `notes` で返す。
- **B11-T46（残り）**: 返り値に「ノードの和はクラスの総量ではない」の 1 文を足す。
- **B11-T44**: リモートの失敗で `job_id`・ジョブの記録の場所・`run_dir` を返り値に残す（上流の worker.log は上流の改修が要る）。
  記録はジョブの種類ごとに `remote-job-measure.json`（`pathway_project` の計測）と `remote-job-render.json`
  （`pathway_plot` の描画）に分けた（1 つの `remote-job.json` では描画のジョブが計測の記録を上書きする）。
- **B11-T45**: `mslm_bridge.R` の `requireNamespace` の失敗を、未導入と読み込みの失敗に分け、失敗の本文を返す
  （新しいコード `MSLM_R_PACKAGE_LOAD_FAILED`）。
- **B11-T42**: handoff の as_is の経路でも名前の頭（`(` で始まる FA 名など）を確かめて落とし、理由別の件数に入れる。pathway.md:17 を合わせる。
- **C8-M31**: `dataset_merge_polarities` は `dataset_quant_apply` / skip を通っていないデータセットも受ける（両方が未正規化なら
  そのまま統合し、記録に書く）。正規化の状態が食い違うときの `STATE_MISMATCH` は残す（混ぜない）。
- **C15-M16**: `curation_submit` の `review_id` を省いたら、そのアラインメントの最新のレビューを使う（返り値に出す）。
- **C22-M34**: 「同意なしに呼ばない」を説明文から外し、`source` の定義（`llm`）を「LLM が出した」に直す。
- **C5-M08（キュレーション側）**: `curation_suggest` の `library.path` をファイル名と sha256 にする。

### main（統合の後、自分で）

- **新・MCP_INSTRUCTIONS**: Claude Code は 4,096 字で切る。意味の定義と入口だけを 4,000 字以内に収める。手順・作法
  （C20 の比較の前の `dataset_status`・C23 の「判定の数字は道具から」・FIRST OVERVIEW の提示の順・PARALLEL の作法）は
  output-format の文書と研究側の文書へ移す。残す縛りは「Console の前に許可」だけ。「MUST」「never verdicts」をやめる（C9・C23）。
- **C9-M01**: ダイジェスト（`caveat`）は情報なので残す（種別ごとに 1 回）。止める手段は `caveat_mode=off`。
- **E4-T32**: `docs_read` に任意の `section`（見出しで 1 節だけ返す。省けば今どおり全文）を足す。
- **E3-T22**: `tests/test_sample_manifest_v2.py` の古い説明（REPEATED_MEASURES_UNSUPPORTED）を直す。
- **文書とコードの食い違い**: reanalysis.md:232-233（組成は assay_kind でラベルだけ変わる）ほか。
- **B13-T68**: Claude Desktop の設定のサーバを `.venv` の Python にする（追跡外・利用者の設定ファイル）。
- 統合の後、CLAUDE.md・USAGE.md・vault の流れ図・HISTRY・task.md を直す。

## 共有の部品（main に先に置く）

- `metabolomix/core/polarity_names.py`（stdlib だけの leaf）: `normalize_polarity(value, *, field="polarity")` は
  `positive`/`pos`/`+`/`negative`/`neg`/`-`（大文字小文字を問わない）を `"positive"`/`"negative"` にし、それ以外は
  `ValueError`。担当 1・5（とほかの受け口）が使う。

## 統合

1. 担当ごとに worktree で実装し、関係するテストと全テストを通してコミットする（pre-commit）。
2. 担当ごとに別の目でレビューし、指摘を確かめてから直す。
3. main で直列に `--no-ff` で統合し（衝突は main で解く）、全テストを通す。
4. MCP_INSTRUCTIONS・文書・記録を main で直す。

## 統合で決めたこと

9 本を `fix/inventory` へ統合するときに、担当をまたいで決めたこと（実装はこのとおり）。

- **`drift_mode="auto"` は尋ねない**: 既定の前処理は `needs_input` を出さない。補正の条件を満たしていても既定で
  当てないとき（換算済みの値・注入順の出所が確かでない など）は補正せずに作り、理由と `drift_mode="on"` で
  呼び直せることを返す（担当 3 の A6・C3）。
- **交絡は印であって関門ではない**: `dataset_compare` と spec v2 の比較は `allow_confounded` の既定を True にし、完全交絡でも
  比べて先頭の `confounding` と対比ごとの `confounded` に印を付ける。止めるのは `allow_confounded=false` を明示したとき
  だけ（`STATISTIC_CONFOUNDED`。担当 4 の C4）。最初の概観（`dataset_overview`）だけは、全群をまとめた表で完全交絡
  （と探索専用のデータセット）なら比較の段を回さず `omitted` に理由を書く。
- **行列は直近の前処理のものを使う**: `matrix_id` を省いた `dataset_compare` / `dataset_pca` は、そのハンドルで
  `dataset_preprocess` が最後に作った行列を使い（無ければ既定の前処理）、`matrix_source`（`matrix_id` /
  `latest_preprocess` / `default`）で示す（担当 4 の C12）。
- **探索専用のデータセットも書き出せる**: 未完了の実行・検証を通っていない pipeline の単位から開いたデータセットの
  書き出しと根拠バンドルを断らず、メタ行 `# exploratory = true` と manifest の `exploratory` に理由を書く（担当 7 の C19。
  `EXPLORATORY_ONLY_DATASET` は消した）。定量の派生（正規化・統合・選別）も親の印を引き継ぐ（統合はどちらかの親が
  探索専用なら探索専用。統合後の追補）。
- **`assay_kind` はハンドルごと**: アッセイ種別はデータセット（ハンドル）が持ち、`record_objective` は直近か渡した
  ハンドルへ書き、`load_dataset` は解析フォルダの目的ファイルから読み戻す。`session.assay_kind` は直近のハンドルの値を
  読むだけ（担当 8 の C1）。
- **ツールは 58**: 担当 3 の `pipeline_sample_sheet` と担当 8 の `job_status` を足した合計（登録数・`USAGE.md`・
  CLAUDE.md の冒頭は統合で合算した）。
