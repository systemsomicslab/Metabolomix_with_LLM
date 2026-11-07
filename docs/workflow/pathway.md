# ワークフロー: パスウェイ投影（MSLipidMapper）

前処理済みの脂質データを MSLipidMapper（R パッケージ）に渡し、
同梱パスウェイ（remodeling / ceramide / global）へ投影する経路。MSLipidMapper は
`[mslipidmapper] mslipidmapper_api_url` の研究室サーバ（REST API）を先に使い、使えなければローカルの Rscript へ回す。前提状態の連鎖は
`record_objective(assay_kind="lipid")` ＋ `load_dataset` → `dataset_preprocess` →
`pathway_project` → `pathway_plot`。値の**意味**は
[pathway.md（output_format）](../output_format/pathway.md) を見ること。

`session.pathway`（他スロットとは独立）は直近の `run_id` だけを持つ。run の正準は
`reports/pathway/<run_id>/manifest.json`（置き場所は `[paths] reports_dir` → データセットの `analysis_dir` 配下。
`pathway_project` / `pathway_plot` は `dataset_handle` を持つ）。

入口は `bridge.run_bridge()` だけで、ここがサーバかローカルかを決める（`api_url()`＝`[mslipidmapper] mslipidmapper_api_url` が未設定ならローカルだけ）。
サーバ側は `metabolomix/pathway/remote.py`（`GET /healthz` → `POST /v1/jobs` → 状態を succeeded / failed まで
→ 成果物の取得。config は JSON で送る）。request はローカルと同じ dict で、`_remote_config()` が API の
config に写す（`analysis.normalization` は常に none。サーバは 3 パスウェイを全部計算し、頼んだ分だけ落とす）。
サーバのどんな失敗でもローカルへ回し、理由を戻り値の `remote_fallback` に残す。ローカルも使えない
（`MSLM_NOT_CONFIGURED` / `MSLM_R_PACKAGES_MISSING`）ときはサーバ側のエラーを上げる
（`details.local` にローカルの理由）。戻り値の `backend` がどちらで動いたかを示す。

ローカルの R 側は `metabolomix/pathway/mslm_bridge.R` で、起動直後に renv の仮想環境 `r/` を有効化し、
上流の公開関数（`MSLipidMapper::` で始まるもの）だけを呼ぶ。Python との契約は
run ディレクトリの `request.json` → `bridge-result.json`。

## pathway_project

前提: `session.assay_kind == "lipid"`（違えば `missing_state` → `record_objective`）。
入力は直近のデータセットの `DatasetState.pp_matrix`（無ければ `missing_state` → `load_dataset`、
前処理が無ければ `dataset_preprocess`）。旧 ARF 経路（`source="arf"`）は B-2b Task 8 で消えた。
状態変更: `session.pathway.last_run_id` / `last_run_dir`。
R の仮想環境は起動前に確かめる（`_check_environment()`）: `r/renv/library` に MSLipidMapper が無ければ
`MSLM_R_PACKAGES_MISSING`（`reason=not_restored`）、`r/setup.R` が残す目印（復元に使った `renv.lock` の md5、
`renv/library/restored-lock.md5`）が今の `renv.lock` と違えば同じコードで `reason=stale`。どちらも `Rscript r/setup.R` で直る。
描画のブラウザは設定が無ければ Edge / Chrome の標準の置き場所から探す（`browser_path()`）。

1. metabolomix/pathway/tools.py  pathway_project()
2. ├─ metabolomix/pathway/tools.py  _build_dataset()
3. │  └─ metabolomix/pathway/handoff.py  build_from_dataset()
4. ├─ metabolomix/pathway/handoff.py  write_mztab()
5. ├─ metabolomix/pathway/bridge.py  run_bridge()
6. │  ├─ [mslipidmapper_api_url あり] metabolomix/pathway/bridge.py  _run_remote()
7. │  │  ├─ metabolomix/pathway/bridge.py  _remote_config()
8. │  │  └─ metabolomix/pathway/remote.py  run_job()（output.format=jsonld）
9. │  └─ [mslipidmapper_api_url なし・サーバ失敗] metabolomix/pathway/bridge.py  _run_local()
10. │     ├─ metabolomix/pathway/bridge.py  _check_environment()
11. │     └─ （R 側。mode=measure で `mslm_bridge.R` が `MSLipidMapper::mztab_pathway_measurements` を呼ぶ）
12. ├─ metabolomix/pathway/jsonld.py  load_graph()
13. ├─ metabolomix/pathway/compare.py  compare_scopes()
14. │  └─ metabolomix/analysis/differential.py  two_group_test()
15. └─ metabolomix/pathway/compare.py  opposite_pairs()
16.    └─ metabolomix/pathway/jsonld.py  adjacent_lipid_pairs()

## pathway_plot

前提: `pathway_project` の run（`run_id` 省略時は `session.pathway.last_run_id`）。
無ければ `missing_state` → `pathway_project`。`run_id` は `pathway_project` が採番する形
（`pw_YYYYmmdd_HHMMSS_xxxxxx`）以外を `PATHWAY_BAD_REQUEST` で拒む（`reports/pathway` の外を
指させない）。`inline_image` が bool でなければ `PATHWAY_BAD_REQUEST`。状態変更: なし（PDF と PNG を run ディレクトリに書く）。

1. metabolomix/pathway/tools.py  pathway_plot()
2. ├─ metabolomix/pathway/tools.py  _load_manifest()
3. ├─ metabolomix/pathway/bridge.py  run_bridge()
4. │  ├─ [mslipidmapper_api_url あり] metabolomix/pathway/remote.py  run_job()（output.format=pdf。ジョブ失敗は `PATHWAY_RENDER_FAILED`）
5. │  └─ [mslipidmapper_api_url なし・サーバ失敗] （R 側。mode=render で `mslm_bridge.R` が `MSLipidMapper::render_pathway_from_config` を呼ぶ）
6. └─ metabolomix/plots/pdf_raster.py  pdf_first_page_png()（PNG を `<pathway>.png` に書く。戻り値は既定で場所と要約の JSON、`inline_image=true` は [JSON, 画像…]）
