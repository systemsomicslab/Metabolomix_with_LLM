# ワークフロー: LC–MS メタボロミクス v2

profile（`lcms-profile.v1`）を前提に、注入ごとの測定証拠・feature binding・
内部標準比・固定母集団 QC・v2 統計までを回す経路。

**公開 `pipeline_plan` / `pipeline_run` に v2 要求を渡す経路は接続済み**で、
profile が method・依存・実行体・極性の唯一の情報源になる（v1 の推定経路へは落ちない）。
検証は合成入力と fake Console まで——**実 Console 接続と profile の科学的検証は未了**で、
合成入力での完走をそれらの合格と読まないこと。

**routine は引き続き `PROFILE_VALIDATION_INVALID` で停止する**（証明書の独立した
実測hashの受付が未接続）。draft profile + `execution_purpose="validation"` を使う。

統計までの対話経路は pipeline 抜きでも通る:
`dataset_load` → `dataset_set_sample_metadata` → `dataset_build_matrix` → `dataset_statistic`。
この経路は binding と注入証拠を持たないぶん、内部標準比を作れず検出率 filter も評価不能になる。
完了した v2 run がある場合は `dataset_load(pipeline_path=...)` が、その run の dataset・
試料対応表・binding・解析行列をまとめてセッションへ載せる（[mztab.md](mztab.md)）。
検証範囲と残タスクは[実装監査](../superpowers/plans/2026-09-16-validated-lcms-metabolomics-audit.md)を参照。

v1（`dataset_differential` ほか、[dataset_analysis.md](dataset_analysis.md)）とは
**数値契約が違う**。同じ「log2FC」という名前でも、v1 は `log2(x + pseudo_count)` 空間の
平均差、v2 は統計変換前の算術平均比の log2（`effect_size_definition =
log2_arithmetic_mean_ratio`）。v1 の既定値・契約は一切変更していない。

向きは v1・v2 とも揃えてある: `log2fc`、Welch の `t_statistic`（v1 は結果行の `t`）、
Tukey の `mean_difference` は、いずれも**正なら test（v1 は group_b）が高い**。
t は `(平均_test − 平均_reference) / SE`（2026-09-27 に向きを揃えた——それ以前の
t は符号が逆）。

## プロジェクト保存（`save_project`）と Application Control

要求の `save_project` は**既定 true**（v1/v2 とも）。MS-DIAL はこのとき `-p` を
受け取り、全検体の解析とアライメントを終えた**後**にプロジェクトを書く。その段で
`MsdialLcImMsApi` を読むため、**Windows の Smart App Control が有効かつ MS-DIAL が
ローカルビルド（未署名）**だと、ここだけがポリシーに弾かれる。

症状は最悪の形になる——長い解析を全部終えてから落ちる。実測では 41 検体で
13 分半を費やした後に `MSDIAL_EXECUTION_FAILED exit_code=1` になった。

**計画時に止めるようにしてある。** `inputs._assert_project_save_possible` が
「ポリシーが Enforce」かつ「当該アセンブリが未署名」の両方を確認したときだけ
`PROJECT_SAVE_BLOCKED` で停止する（`lipidmix/core/app_control.py`）。ポリシー単独では
止めない——署名済みの公式配布版を使う正当な構成まで塞いでしまうため。

**対処は要求に `"save_project": false` を書く。** `.mdproject` は GUI で開くための
便宜で、v2 の必須成果物 10 種には含まれない（`report._V2_REQUIRED_OUTPUTS`）。

自分の環境を確かめるには:

```powershell
# 0=Off / 1=Enforce / 2=評価
(Get-ItemProperty "HKLM:\SYSTEM\CurrentControlSet\Control\CI\Policy" `
  -Name VerifiedAndReputablePolicyState).VerifiedAndReputablePolicyState
```

## dataset_build_matrix

セッションの `DatasetState` から解析行列（`analysis-matrix.v1`）を1本作り、
`ds.analysis_matrices` へ登録する。pipeline を回さずに v2 統計へ進むための入口。

- **前提** — `session.dataset` のみ（`dataset_load`）。群・batch・注入順・include は
  `dataset_set_sample_metadata` で入れておく（正規化の参照試料とドリフト補正の
  前提判定がこれを読む）。
- **recipe の検証** — profile の `matrix_recipes` と**同じ規則**を共有する
  （`profile_schema.validate_matrix_recipe`）。経路によって通る recipe が変わると、
  同じ `recipe_id` の行列が別物になる。
- **揃っていない前提を埋めない** — `base="internal_standard_ratio"` は対象feature↔
  内部標準の対応付けを要求するので `missing_state`（`feature_bindings`）で止める。
  検出状態が不明な dataset では `filter.min_detection_rate` は `not_evaluable` として
  記録され、feature を落とさない。
- **状態変更** — `ds.analysis_matrices[matrix_id]`。値は戻り値に載せない。
- **呼び出し連鎖**

```
1. lipidmix/tools/dataset_analysis_tools.py  dataset_build_matrix()
2. ├─ lipidmix/console/profile_schema.py  validate_matrix_recipe()
3. └─ lipidmix/analysis/dataset_service.py  build_analysis_matrix()
4.    └─ lipidmix/analysis/matrix_state.py  make_matrix()
5.       └─ lipidmix/analysis/matrix_state.py  _detected_mask()
6.       └─ lipidmix/analysis/matrix_state.py  _apply_recipe()
7.       └─ lipidmix/analysis/matrix_state.py  _detection_eligibility()
8.       └─ lipidmix/analysis/internal_standards.py  apply_internal_standards()
```

`matrix_id` は値ではなく**入力の identity**（dataset fingerprint・recipe・binding・
証拠・metadata・eligibility）の hash。同じ入力なら同じ ID になるので、pipeline が
作った行列と突き合わせられる。

## dataset_statistic

名指しした解析行列（`analysis-matrix.v1`）に対して統計を1件実行する。

- **前提** — `session.dataset`（`dataset_load`）と、`ds.analysis_matrices` に登録済みの
  解析行列。行列の出どころは3つ: pipeline の `preprocess` / `qc_processed` 工程が作った
  ものを `dataset_load(pipeline_path=...)` で載せる、`dataset_build_matrix` でこの場で作る、
  のいずれか。どちらが欠けても `missing_state` 封筒を返す。
- **状態変更** — `session.dataset.results["stat_<statistic_id>"]` に結果全量を置く。
  戻り値には要約だけを載せる。
- **呼び出し連鎖**

```
1. lipidmix/tools/dataset_analysis_tools.py  dataset_statistic()
2. └─ lipidmix/analysis/dataset_service.py  statistic_dataset()
3.    └─ lipidmix/analysis/statistics_v2.py  run_statistic()
4.       └─ lipidmix/analysis/statistics_v2.py  transform_values()
5.       └─ lipidmix/analysis/sample_manifest.py  select_statistical_samples()
6.       └─ lipidmix/analysis/sample_manifest.py  validate_independent_samples()
7.       └─ lipidmix/analysis/differential.py  welch_t()
8.       └─ lipidmix/analysis/differential.py  bh_fdr()
9.       └─ lipidmix/analysis/multigroup.py  test_feature()
10.      └─ lipidmix/analysis/statistics_v2.py  arithmetic_log2fc()
11.      └─ lipidmix/analysis/pca.py  run_pca()
```

`matrix_result_id` が `ds.analysis_matrices` に無ければ `ANALYSIS_RESULT_NOT_FOUND` で
止まる。**近い行列で代用しない**——recipe 違いの行列が2本ある前提の設計で、
どちらの数字かを言えない結果を返すくらいなら止めるほうが安全。

## pipeline の v2 工程（参考）

公開受付接続後にv2 workerが回す予定の順序（handler の実体は
`lipidmix/pipeline/metabolomics_handlers.py`）:

```
prepare_inputs → execute_console → validate_outputs → load_dataset
  → resolve_metadata → load_assay_evidence → resolve_feature_bindings
  → qc_raw → preprocess → qc_processed → statistics:<id> → export:<id> → report
```

- `load_assay_evidence` — `.arf` の `AlignedPeakProperties` から注入ごとの RT/m/z を
  読む（`lipidmix/analysis/assay_evidence.py` `build_assay_evidence()`）。取得できない
  場合も理由 JSON を残して先へ進む。
- `resolve_feature_bindings` — profile の `feature_targets` をこのバッチの feature へ
  対応付ける（`lipidmix/analysis/feature_bindings.py` `bind_features()`）。0件・複数件は
  `FEATURE_BINDING_UNRESOLVED` で停止し、自動では選ばない。
- `qc_raw` — filter 前に QC 評価集合を固定する（`lipidmix/analysis/assay_qc.py`
  `evaluate_qc()`）。
- `preprocess` — recipe ごとに行列を作る（`lipidmix/analysis/matrix_state.py`
  `make_matrix()`）。補完はしない。
- `qc_processed` — 処理後 QC を評価してから補完する（`finalize_matrix()`）。QC で落ちた
  feature は eligibility にだけ反映し、QC は再集計しない。
- `report` — 実行・QC・解析の3軸を分けて書く（`lipidmix/pipeline/report.py`
  `analysis_status()`）。

工程ごとの停止理由・再開の規則は
[../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md](../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md)
（spec §6.2・§8〜§11）が正準。
