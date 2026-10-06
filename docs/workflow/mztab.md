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
9.    └─ metabolomix/mztab/identity.py derive_inchikey()
10.   └─ metabolomix/mztab/dataset_state.py _build_feature_annotations()
11.      └─ metabolomix/mztab/identity.py derive_inchikey()
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
7.    └─ metabolomix/mztab/validator.py detect_quantification_measure()
8.    └─ metabolomix/mztab/dataset_state.py build_dataset_state()
9.       └─ metabolomix/mztab/reader.py extract_abundance_matrix()
10.      └─ metabolomix/mztab/dataset_state.py _index_sme_rows()
11.      └─ metabolomix/mztab/dataset_state.py _best_evidence()
12.      └─ metabolomix/mztab/identity.py derive_inchikey()
13.      └─ metabolomix/mztab/dataset_state.py _build_feature_annotations()
14.         └─ metabolomix/mztab/identity.py derive_inchikey()
15.      └─ metabolomix/mztab/dataset_state.py _resolve_sample_names()
16. └─ metabolomix/mztab/evidence.py attach_to_dataset()

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

dataset・解決済み試料対応表・binding・解析行列を**同じ run から**載せる。
行列だけを別ツールで載せる入口は置いていない——run A の dataset に run B の
行列を混ぜられると、群の対応が黙ってずれた統計が出る。行列の値は
`save_matrix` が残した npz から読み、meta・配列のどちらが書き換わっても
`load_matrix` が拒否する。

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
