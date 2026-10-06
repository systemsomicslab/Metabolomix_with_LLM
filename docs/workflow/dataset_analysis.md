# DatasetState 解析ツール呼び出し連鎖

ARF 経路（`docs/workflow/arf.md`）と同じ純関数を共有している。同じ入力からは
同じ数字が出る。違いは入口（mzTab-M か .arf か）とセッションスロットだけ。

## dataset_preprocess

1. metabolomix/tools/dataset_analysis_tools.py  dataset_preprocess()
2. └─ metabolomix/analysis/dataset_analysis.py  run_dataset_preprocess()
3.    └─ metabolomix/analysis/dataset_analysis.py  build_dataset_pp_inputs()
4.       └─ metabolomix/analysis/preprocessing.py  detect_sample_roles()
5.    └─ metabolomix/analysis/dataset_analysis.py  _apply_detection_filter()
6.       └─ metabolomix/analysis/preprocessing.py  detection_rates()
7.    └─ metabolomix/analysis/preprocessing.py  preprocess()
8.    └─ metabolomix/analysis/preprocessing.py  drop_samples_by_role()
9.    └─ metabolomix/analysis/preprocessing.py  detect_qc_strata()

手順 5 は前処理（手順 7）より**前**に来る。正規化・補完のあとでは gap-fill セルが
実測値と区別できなくなり、何を根拠に特徴を残したかが言えなくなるため。
`min_detection_rate > 0` なのに `ds.detected_mask` が無い場合は `bad_request` で
拒否する（検出状態が「無い」のと「全部未検出」は解釈が正反対で、0 扱いで通すと
gap-fill だけの特徴を実測として数えた行列が黙って下流に流れる）。

検出率の分母は**実際に前処理へ渡す試料だけ**。手順 3 が `include=false` の行を
行列から落とすのと同じ位置（`_included_sample_indices`）でマスクの列も絞る
（`_restrict_mask_to_samples`）——絞らないと、解析から除外した試料の未検出が
分母に残り、`min_detection_rate=1.0` が「残す試料では全件検出されている」正当な
特徴量を削り落とす。何試料分を数えたかは `report["detection"]["n_samples"]`。

## dataset_pca

1. metabolomix/tools/dataset_analysis_tools.py  dataset_pca()
2. └─ metabolomix/analysis/dataset_analysis.py  run_dataset_pca()
3.    └─ metabolomix/analysis/pca.py  run_pca()

## dataset_differential

1. metabolomix/tools/dataset_analysis_tools.py  dataset_differential()
2. └─ metabolomix/analysis/dataset_analysis.py  run_dataset_differential()
3.    └─ metabolomix/analysis/differential.py  check_confounding()
4.    └─ metabolomix/analysis/differential.py  two_group_test()
5.    └─ metabolomix/analysis/differential.py  add_fdr()
6.    └─ metabolomix/analysis/differential.py  summarize_two_group()
7.    └─ metabolomix/analysis/differential.py  volcano_data()

## dataset_export_differential

書き出す前に 2 つ確かめる。**その結果が今の前処理から出たものか**（手順 3。前処理を
やり直した後の古い結果を書くと、TSV の数字と現在の前処理条件が並んだファイルができ、
どちらも正しく見えてずれが分からない）と、**探索専用のデータセットでないか**
（中断された実行の出力では、欠けた検体が「その群には無い」ようにしか見えない）。
メタ行の前処理条件は `ds.preprocessing_recipe`（今の状態）ではなく、その結果が親に
持つ前処理の `provenance.effective_parameters` から取る。出力は一時ファイルへ書いてから
置換するので、途中で落ちた出力が「新しい完成品」として残らない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_export_differential()
2. └─ metabolomix/analysis/dataset_export.py  export_dataset_result()
3. │  └─ metabolomix/analysis/result_state.py  assert_current()
4. │  └─ metabolomix/analysis/export_contract.py  is_significant()
5. │  └─ metabolomix/analysis/dataset_export.py  _meta_lines()
6. │  │  └─ metabolomix/analysis/dataset_export.py  _preprocess_parameters()
7. │  │  └─ metabolomix/analysis/export_contract.py  build_meta()
8. │  └─ metabolomix/analysis/export_contract.py  format_row()
9. │  └─ metabolomix/analysis/dataset_export.py  _atomic_write_text()

## dataset_set_sample_metadata

実験情報シート(sample-manifest.v1 / v2)を読み、全件検証してから一括反映する。版は
先頭の schema 行で決まり、列見出しからは推測しない。検証で1件でも落ちれば
`apply_metadata` の契約により `session.dataset` は一切変更しない。

1. metabolomix/tools/dataset_analysis_tools.py  dataset_set_sample_metadata()
2. └─ metabolomix/analysis/dataset_service.py  apply_sample_manifest()
3.    └─ metabolomix/analysis/dataset_service.py  _raw_manifest_layout()
4.    └─ metabolomix/analysis/sample_manifest.py  parse_manifest()
5.    └─ metabolomix/analysis/sample_manifest.py  resolve_metadata()
6.    └─ metabolomix/analysis/sample_manifest.py  apply_metadata()
7.       └─ metabolomix/analysis/result_state.py  invalidate_results()

pipeline が比較を明示するときに使う `resolve_comparison`/`run_comparison`
（`metabolomix/analysis/dataset_service.py`）は、群名だけで対照/処置の向きを決めない・
QC/blank/unknown/standard・`include=false` を混ぜない（選択規則は
`metabolomix/analysis/sample_manifest.py` `select_statistical_samples()` が唯一の定義）・
完全交絡を `allow_confounded=true` の
明示なしには通さない、という前提検証を `compare_dataset()` の前に挟む。単体ツール
`dataset_differential` はこの前提検証を経ない汎用呼び出し（USAGE.md 参照）。
