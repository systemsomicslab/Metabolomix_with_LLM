# mzTab-M ツール呼び出し連鎖

ここは**どのファイルのどの関数をどの順に呼ぶか**。出てくる値の**意味**
（SME と SML の区別、`identification_status`、`name_source`、`inchikey_coverage`、
gap-fill）は `docs/output_format/mztab.md`（`lipidmix://docs/output-format/mztab`）。

## dataset_load

読み込みの実体は `metabolomix/mztab/loading.py`（session に依存しない層）。MCP ツールは
戻ってきた DatasetState を成功時だけ `session.dataset` へ置く。pipeline のワーカーは
同じ関数を呼び、自分が持つ DatasetState として使う。

mztab_path 経路（直接指定）:
1. metabolomix/tools/mztab_tools.py dataset_load()
2. └─ metabolomix/mztab/loading.py load_dataset_state()
3. └─ metabolomix/mztab/reader.py parse_mztab()
3. └─ metabolomix/mztab/validator.py validate_mztab()
4. └─ metabolomix/mztab/validator.py detect_quantification_measure()
5. └─ metabolomix/mztab/dataset_state.py build_dataset_state()
6.    └─ metabolomix/mztab/reader.py extract_abundance_matrix()
7.    └─ metabolomix/mztab/dataset_state.py _index_sme_rows()
8.    └─ metabolomix/mztab/dataset_state.py _best_evidence()
9.    └─ metabolomix/mztab/dataset_state.py _resolve_inchikey()（MSP 表なし。中で derive_inchikey）
10.   └─ metabolomix/mztab/dataset_state.py _build_feature_annotations()
11.      └─ metabolomix/mztab/dataset_state.py _resolve_inchikey()
12.   └─ metabolomix/mztab/dataset_state.py _resolve_sample_names()
13. └─ metabolomix/mztab/evidence.py attach_to_dataset()
14.    └─ metabolomix/mztab/evidence.py arf_candidates()
15.    └─ metabolomix/mztab/evidence.py load_arf_evidence()
16.       └─ metabolomix/arf/reader.py deserialize()
17.       └─ metabolomix/mztab/evidence.py normalize_arf_spots()
18.       └─ metabolomix/mztab/evidence.py build_evidence()
19.    └─ metabolomix/mztab/evidence.py apply_evidence()

job_path 経路（analysis-job.json の宣言 polarity + measure で正準を選ぶ）:
1. metabolomix/tools/mztab_tools.py dataset_load()
2. └─ metabolomix/mztab/loading.py load_dataset_state()
3.    └─ metabolomix/handoff/schema.py AnalysisJob.load()
4.    └─ metabolomix/mztab/loading.py select_primary_entry()
5.    └─ metabolomix/mztab/reader.py parse_mztab()
6.    └─ metabolomix/mztab/validator.py validate_mztab()
7.    └─ metabolomix/mztab/loading.py _msp_table_for_job()
8.       └─ metabolomix/console/method_file.py read_method_keys()
9.       └─ metabolomix/mztab/msp_inchikey.py load_msp_inchikey_table()
10.   └─ metabolomix/mztab/validator.py detect_quantification_measure()
11.   └─ metabolomix/mztab/dataset_state.py build_dataset_state()
12.      └─ metabolomix/mztab/reader.py extract_abundance_matrix()
13.      └─ metabolomix/mztab/dataset_state.py _index_sme_rows()
14.      └─ metabolomix/mztab/dataset_state.py _best_evidence()
15.      └─ metabolomix/mztab/dataset_state.py _resolve_inchikey()
16.         └─ metabolomix/mztab/msp_inchikey.py lookup()
17.         └─ metabolomix/mztab/identity.py derive_inchikey()
18.      └─ metabolomix/mztab/dataset_state.py _build_feature_annotations()
19.         └─ metabolomix/mztab/dataset_state.py _resolve_inchikey()
20.      └─ metabolomix/mztab/dataset_state.py _resolve_sample_names()
21. └─ metabolomix/mztab/evidence.py attach_to_dataset()

pipeline_path 経路（完了した run の成果物をまとめて載せる）:
1. metabolomix/tools/mztab_tools.py dataset_load()
2. └─ metabolomix/pipeline/session_handoff.py load_run_dataset()
3.    └─ metabolomix/pipeline/recovery.py normalize_pipeline_root()
4.    └─ metabolomix/pipeline/store.py load_run()
5.    └─ metabolomix/mztab/loading.py load_dataset_state()
6.    └─ metabolomix/pipeline/store.py read_result_data()
7.       └─ metabolomix/pipeline/store.py verify_result_refs()
8.    └─ metabolomix/analysis/feature_bindings.py resolved_targets()
9.    └─ metabolomix/analysis/matrix_state.py load_matrix()
10.   └─ metabolomix/pipeline/session_handoff.py restore_differentials()
11.      └─ metabolomix/pipeline/store.py current_result_refs()
12.      └─ metabolomix/pipeline/store.py read_result_data()
13.      └─ metabolomix/pipeline/session_handoff.py _handoff_rejection()
14.         └─ metabolomix/analysis/result_state.py dataset_fingerprint()
15.         └─ metabolomix/analysis/result_state.py group_fingerprint()
16.      └─ metabolomix/pipeline/session_handoff.py _run_preprocess_parameters()
17.         └─ metabolomix/pipeline/store.py read_result_data()
18. └─ metabolomix/tools/mztab_tools.py _handoff_lines()

dataset・解決済み試料対応表・binding・解析行列を**同じ run から**載せる。
行列だけを別ツールで載せる入口は置いていない——run A の dataset に run B の
行列を混ぜられると、群の対応が黙ってずれた統計が出る。行列の値は
`save_matrix` が残した npz から読み、meta・配列のどちらが書き換わっても
`load_matrix` が拒否する。

run の差次的結果（`differential:<comparison_id>`）は、試料対応表を戻した**後**に
データの指紋（`provenance.dataset_hash`）と群の指紋（`group_fingerprint`）を今の
DatasetState と照合し、一致したものだけを `ds.results` へ戻す。result_id は付け替えない
（pipeline の TSV のメタ行と同じ ID）。来歴に `handoff`（pipeline_root・comparison_id）を
足すので、`result_state.is_current` は ID と前処理ではなく指紋で有効性を見る。
`ds.last_differential` は設定しない（戻した結果は result_id で名指しして使う）。
親の前処理（run の `preprocess` 成果物）は `ds.results` へ戻さず、そのレシピ
（`resolved_recipe`）だけを `provenance.handoff.preprocess` へ写す。差次的結果の親は
最後の前処理とは限らない（resume で前処理は作り直されるが、差次的工程は skip されうる）
ので、ここだけは履歴を result_id で引いて hash を照合する。読めなければ None で、
書き出しのメタ行 `# preprocess = ...` は session の前処理へフォールバックせず
`unknown（…）` と書く（`dataset_export._preprocess_parameters`）。
読むのは `current_result_refs()` が選ぶ今の ref だけ（export が同じ比較を再永続化すると
履歴に古い ref が残るため）で、今の ref の成果物が改変されていれば
`RESULT_INTEGRITY_MISMATCH` で読み込みごと止める。戻さなかったものは
`ds.handoff_report["rejected"]` に理由（`legacy_result_without_dataset_hash` /
`dataset_mismatch` / `group_mismatch`）つきで残り、要約に並ぶ。指紋を持たない旧 run の
結果は `pipeline_resume` で差次的工程だけを再計算させる手順が無いため、session で
`dataset_preprocess` → `dataset_differential` を案内する。

読み込みは「読めたか」だけでなく「どれだけ信用してよいか」も返す。
`source_verification` は終了証跡（exit code・identity）とファイル hash が揃って初めて
`verified` になり、ジョブ経由でも証跡が無ければ `legacy_unverified`、直接読みは
`direct_unverified`。**`status == "completed"` という文字列は根拠にしない**——旧経路の
completed は「実行後にファイルが増えた」以上の意味を持たなかった。

完了していない実行（partial / failed / running）の出力は既定で読まない
（`INCOMPLETE_ANALYSIS_JOB`）。`allow_incomplete=True` で読んだ場合は
`exploratory_only=True` になり、2 群比較（dataset_differential）と差次的エクスポートは
`EXPLORATORY_ONLY_DATASET` で拒否する。欠けた検体は mzTab 上「その群に無い」ように
しか見えず、比較はその欠落ごと結論にしてしまう。

`select_primary_entry()` は候補が一意に決まらなければ読まずに停止する
（`QUANTIFICATION_CONFLICT` / `POLARITY_MISMATCH` / `AMBIGUOUS_PRIMARY_MZTAB`）。
どのファイルを読むかは解析結果そのものを変えるため、辞書順にも LLM にも決めさせない。

job_path 経路では手順 16（`attach_to_dataset`）が `artifact_paths` を入れ終えた
**後**に来る。handoff が
記録した `peak_matrix_source` を `.arf` 候補の先頭に使えるのはその時点以降だけ。

### evidence sidecar（検出状態）の取り込み

mzTab-M の `abundance_assay[N]` は非ゼロでも実測ピークか gap-fill 補間かを区別しない。
実データ（60 サンプル × 714 特徴）では**セルの 70.0% が gap-fill** だったので、
非ゼロを検出と数えると検出率を 3 倍以上に過大評価する。隣接する `.arf` は
(スポット × サンプル) の粒度で `MasterPeakID < 0` = gap-fill を持つので、そこから補う。

接合は**名前ではなく数値で検証する**。`.arf` のスポット順が mzTab の `SMF_ID` と
一致することを、スポット数の一致と全特徴の m/z 差（既定 0.01 Da 以内）で確かめる。
実データでは位置一致で最大 6.7 mDa、1 つずらすと 87 Da に爆発するので偶然は起きない。
確認できなければ取り込まず、`feature_qc` に理由を残して warning を出す
——誤接合は検出/未検出を特徴間で入れ替えたまま静かに嘘をつくため。

結果は `ds.detected_mask`（(特徴 × サンプル) の bool 行列）と `ds.feature_qc`（要約）。
TSV は書かない（消費者のいない `feature-qc.tsv` を 2026-09-03 に廃止した経緯）。

## dataset_status

1. metabolomix/tools/mztab_tools.py dataset_status()
2. └─ metabolomix/tools/mztab_tools.py _samples_tsv()
3.    └─ metabolomix/analysis/preprocessing.py detect_sample_roles()
4. └─ metabolomix/tools/mztab_tools.py _detection_summary()
