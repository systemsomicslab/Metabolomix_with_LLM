# ワークフロー: パスウェイ投影（MSLipidMapper）

統計エンジンの行列（前処理済みの脂質データ）を MSLipidMapper（R パッケージ）に渡し、同梱パスウェイ（remodeling / ceramide /
global）へ投影する経路。連鎖は `load_dataset` →（`dataset_preprocess`。省略可）→ `pathway_project` → `pathway_plot`
（`assay_kind` は問わない）。値の意味は [pathway.md（output_format）](../output_format/pathway.md)。

`session.pathway`（他スロットとは独立）は直近の `run_id` だけを持つ。run の正準は `reports/pathway/<run_id>/manifest.json`
（置き場所は `[paths] reports_dir` → データセットの `analysis_dir` 配下。`pathway_project` / `pathway_plot` は `dataset_handle` を持つ）。

入口は `bridge.run_bridge()` だけで、研究室サーバ（REST API）かローカルの Rscript かを決める（`api_url()`＝
`[mslipidmapper] mslipidmapper_api_url` が未設定ならローカルだけ）。
- サーバ側は `metabolomix/pathway/remote.py`（`GET /healthz` → `POST /v1/jobs` → 状態を succeeded / failed まで → 成果物の取得。
  config は JSON）。request はローカルと同じ dict で、`_remote_config()` が API の config に写す（`analysis.normalization` は常に
  none。サーバは 3 パスウェイを全部計算し、頼んだ分だけ落とす）。
- サーバのどんな失敗でもローカルへ回し、理由を `run_bridge()` の戻り値の `remote_fallback`（`job_id` は投入まで進んだときだけ、
  `record` はそのジョブの送受信の記録の場所）に残す。記録の名前は `job_record_name()` が mode ごとに決める
  （`remote-job-measure.json` / `remote-job-render.json`）。ローカルも使えない（`MSLM_NOT_CONFIGURED` / `MSLM_R_PACKAGES_MISSING`）
  ときはサーバ側のエラーを上げ（`details.remote` にサーバの理由、`details.local` にローカルの理由）、ローカルの R が動いて失敗した
  とき（`MSLM_R_PACKAGE_LOAD_FAILED` を含む）はローカルのエラーを上げる（サーバの理由は `details.remote_fallback`）。戻り値の
  `backend`・`details.run_dir` は output_format の「失敗したとき」。
- ローカルの R 側は `metabolomix/pathway/mslm_bridge.R` で、起動直後に renv の仮想環境 `r/` を有効化し、上流の公開関数
  （`MSLipidMapper::`）だけを呼ぶ。Python との契約は run ディレクトリの `request.json` → `bridge-result.json`。MSLipidMapper が
  仮想環境に無ければ（`find.package`）`MSLM_R_PACKAGES_MISSING`、読み込めなければ（`loadNamespace`）`MSLM_R_PACKAGE_LOAD_FAILED`
  （message に R の本文）。

## pathway_project

前提: なし（`assay_kind` が lipid でなければ戻り値と manifest の `notes.assay_kind` に注意）。データセットが無ければ
`missing_state` → `load_dataset`、名指しの行列が無い・古いなら `missing_state("stat_matrix")` → `dataset_preprocess`。
状態変更: `session.pathway.last_run_id` / `last_run_dir`。

- 行列は `dataset_compare` / `dataset_pca` と同じ決め方（手順 3 の `matrix_inputs(prefer_latest=True)`）: `matrix_id` → このハンドルで
  `dataset_preprocess` が最後に作った行列（手順 4）→ 既定の前処理（手順 5）。既定の行列は要求の検査（試料・単位の数・
  aggregation・q_threshold・method）がすべて済んでから登録し、前処理の記録を書く（手順 8 の `register_default`）。
  `matrix_source`・`drift_not_applied_by_default`・エラーの `details.matrix_id` / `matrix_note`（手順 6 の受け渡しの検査を含む）は
  output_format の「ノード表（nodes_tsv）」。
- MSLipidMapper へは行列の注入ごとの値（生体のみ・欠測は 0）を渡す。分子種の入った鎖スコープがどのパスウェイにも無い
  `acyl_chains` は、手順 17 の後で数えて `missing_chains`・`notes.acyl_chains` に書くだけ（空スコープは手順 18 の比較から外れる）。
- R の仮想環境は起動前に `_check_environment()` が確かめる: `r/renv/library` に MSLipidMapper が無ければ
  `MSLM_R_PACKAGES_MISSING`（`reason=not_restored`）、`r/setup.R` が残す目印（`renv/library/restored-lock.md5`。復元に使った
  `renv.lock` の md5）が今の `renv.lock` と違えば同じコードで `reason=stale`。どちらも `Rscript r/setup.R` で直る。
- 描画のブラウザは設定が無ければ Edge / Chrome の標準の置き場所から探す（`browser_path()`）。

1. metabolomix/pathway/tools.py  pathway_project()
2. ├─ metabolomix/pathway/tools.py  _build_dataset()
3. │  ├─ metabolomix/analysis/stats_results.py  matrix_inputs()
4. │  │  └─ [matrix_id なし] metabolomix/analysis/stats_results.py  latest_preprocess_matrix()
5. │  ├─ [行列なし] metabolomix/analysis/stats_results.py  prepare_default()
6. │  └─ metabolomix/pathway/handoff.py  build_from_dataset()
7. │     └─ metabolomix/pathway/handoff.py  resolve_row()（経路か落とす理由。残す名前の頭は upstream_head_ok() で確かめる）
8. ├─ [行列なし] metabolomix/tools/_matrix_store.py  register_default()
9. ├─ metabolomix/pathway/handoff.py  write_mztab()
10. ├─ metabolomix/pathway/bridge.py  run_bridge()
11. │  ├─ [mslipidmapper_api_url あり] metabolomix/pathway/bridge.py  _run_remote()
12. │  │  ├─ metabolomix/pathway/bridge.py  _remote_config()
13. │  │  └─ metabolomix/pathway/remote.py  run_job()（output.format=jsonld。記録は `remote-job-measure.json`）
14. │  └─ [mslipidmapper_api_url なし・サーバ失敗] metabolomix/pathway/bridge.py  _run_local()
15. │     ├─ metabolomix/pathway/bridge.py  _check_environment()
16. │     └─ （R 側。mode=measure で `mslm_bridge.R` が `MSLipidMapper::mztab_pathway_measurements` を呼ぶ）
17. ├─ metabolomix/pathway/jsonld.py  load_graph()
18. ├─ metabolomix/pathway/compare.py  compare_scopes()
19. │  ├─ metabolomix/analysis/stats/external.py  prepared_from_injections()（注入を行列の単位にまとめる。まとめない行列なら注入ごと）
20. │  └─ metabolomix/analysis/stats/compare.py  compare()（`method` の検定。スコープ種別ごとに 1 回）
21. └─ metabolomix/pathway/compare.py  opposite_pairs()
22.    └─ metabolomix/pathway/jsonld.py  adjacent_lipid_pairs()

## pathway_plot

前提: `pathway_project` の run（`run_id` 省略時は `session.pathway.last_run_id`。無ければ `missing_state` → `pathway_project`）。
`run_id` は `pw_YYYYmmdd_HHMMSS_xxxxxx` の形以外を、`inline_image` は bool 以外を `PATHWAY_BAD_REQUEST` で拒む。
状態変更: なし（PDF と PNG を run ディレクトリに書く）。

1. metabolomix/pathway/tools.py  pathway_plot()
2. ├─ metabolomix/pathway/tools.py  _load_manifest()
3. ├─ metabolomix/pathway/bridge.py  run_bridge()
4. │  ├─ [mslipidmapper_api_url あり] metabolomix/pathway/remote.py  run_job()（output.format=pdf。ジョブ失敗は `PATHWAY_RENDER_FAILED`。記録は `remote-job-render.json`）
5. │  └─ [mslipidmapper_api_url なし・サーバ失敗] （R 側。mode=render で `mslm_bridge.R` が `MSLipidMapper::render_pathway_from_config` を呼ぶ）
6. └─ metabolomix/plots/pdf_raster.py  pdf_first_page_png()（PNG を `<pathway>.png` に書く。戻り値は既定で場所と要約の JSON、`inline_image=true` は [JSON, 画像…]）
