# ワークフロー: パスウェイ投影（MSLipidMapper）

前処理済みの脂質データを、ローカルの Rscript で MSLipidMapper（R パッケージ）に渡し、
同梱パスウェイ（remodeling / ceramide / global）へ投影する経路。前提状態の連鎖は
`record_objective(assay_kind="lipid")` ＋ `arf_preprocess`（または `dataset_preprocess`）→
`pathway_project` → `pathway_plot`。値の**意味**は
[pathway.md（output_format）](../output_format/pathway.md) を見ること。

`session.pathway`（他スロットとは独立）は直近の `run_id` だけを持つ。run の正準は
`reports/pathway/<run_id>/manifest.json`。

R 側は `metabolomix/pathway/mslm_bridge.R` で、起動直後に renv の仮想環境 `r/` を有効化し、
上流の公開関数（`MSLipidMapper::` で始まるもの）だけを呼ぶ。Python との契約は
run ディレクトリの `request.json` → `bridge-result.json`。

## pathway_project

前提: `session.assay_kind == "lipid"`（違えば `missing_state` → `record_objective`）。
arf は `session.arf.feature_matrix` と兄弟 `.arf2`、dataset は `DatasetState.pp_matrix`。
`source="auto"` は前処理済み ARF があれば arf、無ければ `dataset_load` 済みなら dataset、
どちらも無ければ arf（`arf_preprocess` を求める `missing_state`）を選ぶ。
状態変更: `session.pathway.last_run_id` / `last_run_dir`。
R の仮想環境は起動前に確かめる（`_check_environment()`）: `r/renv/library` に MSLipidMapper が無ければ
`MSLM_R_PACKAGES_MISSING`（`reason=not_restored`）、`r/setup.R` が残す目印（復元に使った `renv.lock` の md5、
`renv/library/restored-lock.md5`）が今の `renv.lock` と違えば同じコードで `reason=stale`。どちらも `Rscript r/setup.R` で直る。
描画のブラウザは設定が無ければ Edge / Chrome の標準の置き場所から探す（`browser_path()`）。

1. metabolomix/pathway/tools.py  pathway_project()
2. ├─ [source=arf] metabolomix/pathway/tools.py  _build_arf()
3. │  ├─ metabolomix/arf/tools.py  _pool_group_labels()
4. │  ├─ metabolomix/pathway/tools.py  _arf2_catalog()
5. │  └─ metabolomix/pathway/handoff.py  build_from_arf()
6. ├─ [source=dataset] metabolomix/pathway/tools.py  _build_dataset()
7. │  └─ metabolomix/pathway/handoff.py  build_from_dataset()
8. ├─ metabolomix/pathway/handoff.py  write_mztab()
9. ├─ metabolomix/pathway/bridge.py  run_bridge()
10. │  ├─ metabolomix/pathway/bridge.py  _check_environment()
11. │  └─ （R 側。mode=measure で `mslm_bridge.R` が `MSLipidMapper::mztab_pathway_measurements` を呼ぶ）
12. ├─ metabolomix/pathway/jsonld.py  load_graph()
13. ├─ metabolomix/pathway/compare.py  compare_scopes()
14. │  └─ metabolomix/analysis/differential.py  two_group_test()
15. └─ metabolomix/pathway/compare.py  opposite_pairs()
16.    └─ metabolomix/pathway/jsonld.py  adjacent_lipid_pairs()

## pathway_plot

前提: `pathway_project` の run（`run_id` 省略時は `session.pathway.last_run_id`）。
無ければ `missing_state` → `pathway_project`。`run_id` は `pathway_project` が採番する形
（`pw_YYYYmmdd_HHMMSS_xxxxxx`）以外を `PATHWAY_BAD_REQUEST` で拒む（`reports/pathway` の外を
指させない）。状態変更: なし（PDF を run ディレクトリに書く）。

1. metabolomix/pathway/tools.py  pathway_plot()
2. ├─ metabolomix/pathway/tools.py  _load_manifest()
3. ├─ metabolomix/pathway/bridge.py  run_bridge()
4. │  └─ （R 側。mode=render で `mslm_bridge.R` が `MSLipidMapper::render_pathway_from_config` を呼ぶ）
5. └─ metabolomix/plots/pdf_raster.py  pdf_first_page_png()
