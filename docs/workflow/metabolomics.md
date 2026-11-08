# ワークフロー: LC–MS メタボロミクス v2

profile（`lcms-profile.v1`）を前提に、注入ごとの測定証拠・feature binding・
内部標準比・固定母集団 QC・v2 統計までを回す経路。

**公開 `pipeline_plan` / `pipeline_run` に v2 要求を渡す経路は接続済み**で、
profile が method・依存・実行体・極性の唯一の情報源になる（v1 の推定経路へは落ちない）。
検証は合成入力と fake Console まで——**実 Console 接続と profile の科学的検証は未了**で、
合成入力での完走をそれらの合格と読まないこと。

**routine は引き続き `PROFILE_VALIDATION_INVALID` で停止する**（証明書の独立した
実測hashの受付が未接続）。draft profile + `execution_purpose="validation"` を使う。

**この文書にツールの節は無い**（pipeline 専用の文書）。v2 の解析行列（`analysis-matrix.v1`）と統計
（welch / anova_tukey / pca）は `pipeline_run` が回す工程の中だけで作られ、公開ツールからこの
行列を作る・統計を回す入口はない（旧 `analysis-matrix` / v2 統計の単体ツールは 2026-11 に消えた）。
対話で統計を回すなら統計エンジンの `dataset_preprocess` → `dataset_compare` / `dataset_pca`、
spec で回すなら `dataset_analysis_run`（[dataset_analysis.md](dataset_analysis.md)・
[reanalysis.md](reanalysis.md)）を使う。完了した v2 run がある場合は
`load_dataset(pipeline_path=...)` が、その run の dataset・試料対応表・binding・解析行列を
まとめてセッションへ載せる（[dataset.md](dataset.md)）。
検証範囲と残タスクは[実装監査](../superpowers/plans/2026-09-16-validated-lcms-metabolomics-audit.md)を参照。

pipeline の v2 統計は、v1（pipeline の 2 群比較）や統計エンジン（`dataset_compare`）とは
**数値契約が違う**。同じ「log2FC」という名前でも、v1 は `log2(x + pseudo_count)` 空間の
平均差、v2 は統計変換前の算術平均比の log2（`effect_size_definition =
log2_arithmetic_mean_ratio`）。

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
`PROJECT_SAVE_BLOCKED` で停止する（`metabolomix/core/app_control.py`）。ポリシー単独では
止めない——署名済みの公式配布版を使う正当な構成まで塞いでしまうため。

**対処は要求に `"save_project": false` を書く。** `.mdproject` は GUI で開くための
便宜で、v2 の必須成果物 10 種には含まれない（`report._V2_REQUIRED_OUTPUTS`）。

自分の環境を確かめるには:

```powershell
# 0=Off / 1=Enforce / 2=評価
(Get-ItemProperty "HKLM:\SYSTEM\CurrentControlSet\Control\CI\Policy" `
  -Name VerifiedAndReputablePolicyState).VerifiedAndReputablePolicyState
```

## pipeline の v2 工程（参考）

v2 worker が回す順序（handler の実体は
`metabolomix/pipeline/metabolomics_handlers.py`）:

```
prepare_inputs → execute_console → validate_outputs → load_dataset
  → resolve_metadata → load_assay_evidence → resolve_feature_bindings
  → qc_raw → preprocess → qc_processed → statistics:<id> → export:<id> → report
```

- `load_assay_evidence` — `.arf` の `AlignedPeakProperties` から注入ごとの RT/m/z を
  読む（`metabolomix/analysis/assay_evidence.py` `build_assay_evidence()`）。取得できない
  場合も理由 JSON を残して先へ進む。
- `resolve_feature_bindings` — profile の `feature_targets` をこのバッチの feature へ
  対応付ける（`metabolomix/analysis/feature_bindings.py` `bind_features()`）。0件・複数件は
  `FEATURE_BINDING_UNRESOLVED` で停止し、自動では選ばない。
- `qc_raw` — filter 前に QC 評価集合を固定する（`metabolomix/analysis/assay_qc.py`
  `evaluate_qc()`）。
- `preprocess` — recipe ごとに行列を作る（`metabolomix/analysis/matrix_state.py`
  `make_matrix()`）。補完はしない。
- `qc_processed` — 処理後 QC を評価してから補完する（`finalize_matrix()`）。QC で落ちた
  feature は eligibility にだけ反映し、QC は再集計しない。
- `statistics:<id>` — 指定統計を**最終行列**（`preprocess` / `qc_processed` が作ったもの）へ
  1 件ずつ実行する。評価不能でも結果 JSON は残す（近い行列で代用しない）。

```
1. metabolomix/pipeline/metabolomics_handlers.py  _handle_statistics()
2. └─ metabolomix/analysis/statistics_v2.py  run_statistic()
3.    └─ metabolomix/analysis/statistics_v2.py  transform_values()
4.    └─ metabolomix/analysis/sample_manifest.py  select_statistical_samples()
5.    └─ metabolomix/analysis/sample_manifest.py  validate_independent_samples()
6.    └─ metabolomix/analysis/differential.py  welch_t()
7.    └─ metabolomix/analysis/differential.py  bh_fdr()
8.    └─ metabolomix/analysis/multigroup.py  test_feature()
9.    └─ metabolomix/analysis/statistics_v2.py  arithmetic_log2fc()
10.   └─ metabolomix/analysis/pca.py  run_pca()
```

  この工程のテストは `tests/test_metabolomics_end_to_end.py`（手計算の ANOVA・内部標準比行列への
  Welch）。

- `report` — 実行・QC・解析の3軸を分けて書く（`metabolomix/pipeline/report.py`
  `analysis_status()`）。

工程ごとの停止理由・再開の規則は
[../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md](../superpowers/specs/2026-09-15-validated-lcms-metabolomics-design.md)
（spec §6.2・§8〜§11）が正準。
