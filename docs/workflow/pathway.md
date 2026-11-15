# ワークフロー: パスウェイ投影（MSLipidMapper）

統計エンジンの行列（前処理済みの脂質データ）を MSLipidMapper（R パッケージ）に渡し、
同梱パスウェイ（remodeling / ceramide / global）へ投影する経路。MSLipidMapper は
`[mslipidmapper] mslipidmapper_api_url` の研究室サーバ（REST API）を先に使い、使えなければローカルの Rscript へ回す。前提状態の連鎖は
`load_dataset` →（`dataset_preprocess`。省略時は直近の `dataset_preprocess` の行列、それも無ければ既定の前処理を自分で回す）→
`pathway_project` → `pathway_plot`（`assay_kind` は問わない。lipid でなければ `notes.assay_kind` に名前の文法が前提である旨）。値の**意味**は
[pathway.md（output_format）](../output_format/pathway.md) を見ること。

`session.pathway`（他スロットとは独立）は直近の `run_id` だけを持つ。run の正準は
`reports/pathway/<run_id>/manifest.json`（置き場所は `[paths] reports_dir` → データセットの `analysis_dir` 配下。
`pathway_project` / `pathway_plot` は `dataset_handle` を持つ）。

入口は `bridge.run_bridge()` だけで、ここがサーバかローカルかを決める（`api_url()`＝`[mslipidmapper] mslipidmapper_api_url` が未設定ならローカルだけ）。
サーバ側は `metabolomix/pathway/remote.py`（`GET /healthz` → `POST /v1/jobs` → 状態を succeeded / failed まで
→ 成果物の取得。config は JSON で送る）。request はローカルと同じ dict で、`_remote_config()` が API の
config に写す（`analysis.normalization` は常に none。サーバは 3 パスウェイを全部計算し、頼んだ分だけ落とす）。
サーバのどんな失敗でもローカルへ回し、理由を戻り値の `remote_fallback` に残す（`job_id` は投入まで進んだときだけ、
`record` はそのジョブの送受信の記録の場所。サーバは失敗の本文を返さず、本文は上流の worker.log にだけある）。
記録の名前は `job_record_name()` が request の mode ごとに決める（`remote-job-measure.json` / `remote-job-render.json`。
同じ run で `pathway_plot` の描画のジョブが、manifest の `remote_fallback.record` が指す計測の記録を上書きしない）。ローカルも使えない
（`MSLM_NOT_CONFIGURED` / `MSLM_R_PACKAGES_MISSING`）ときはサーバ側のエラーを上げる
（`details.remote` にサーバの理由、`details.local` にローカルの理由）。ローカルの R が動いて失敗したとき
（`MSLM_R_PACKAGE_LOAD_FAILED` を含む）はローカルのエラーを上げ、サーバの理由は `details.remote_fallback`。
戻り値の `backend` がどちらで動いたかを示す。`pathway_project` / `pathway_plot` はブリッジを呼んでからの失敗の封筒
（ブリッジ・JSON-LD の読み取り・`pathway_plot` の PDF の欠落と画像化の失敗）に `details.run_dir` を足す。

ローカルの R 側は `metabolomix/pathway/mslm_bridge.R` で、起動直後に renv の仮想環境 `r/` を有効化し、
上流の公開関数（`MSLipidMapper::` で始まるもの）だけを呼ぶ。Python との契約は
run ディレクトリの `request.json` → `bridge-result.json`。MSLipidMapper が仮想環境に無ければ（`find.package`）
`MSLM_R_PACKAGES_MISSING`、あるのに読み込めなければ（`loadNamespace` の失敗）`MSLM_R_PACKAGE_LOAD_FAILED` で、
message に R の本文を入れる。

## pathway_project

前提: なし（データセット（ハンドル）の `assay_kind` は問わない。lipid でなければ戻り値と manifest の `notes.assay_kind` に注意を書く）。
入力は直近のデータセットの統計エンジンの行列。行列は `dataset_compare` / `dataset_pca` と同じ決め方（手順 3 の
`matrix_inputs(prefer_latest=True)`）: `matrix_id` を渡せば登録済みのその行列、省けばこのハンドルで `dataset_preprocess` が
最後に作った行列（手順 4）、それも無ければ既定の前処理（手順 5。drift=auto は換算済み・注入順がファイルの順かもしれないとき
ドリフト補正を当てずに計算する。使った行列がそうなら戻り値と manifest の `warnings` に `drift_not_applied_by_default`）。どれを使ったかは戻り値と manifest の `matrix_source`（`matrix_id` / `latest_preprocess` / `default`）。
データセットが無ければ `missing_state` → `load_dataset`、名指しの行列が無い・古いなら `missing_state("stat_matrix")` →
`dataset_preprocess`。既定の行列は、要求の検査（試料・単位の数・aggregation・q_threshold・method）がすべて済んでから登録し、前処理の記録を書く（手順 8 の `register_default`）。`matrix_id` を省いた呼び出しのエラーには、使った行列の `details.matrix_id` と `matrix_source` を添える（手順 6 の受け渡しの検査——サブクラス注釈の付いた特徴の数・群の大きさ——も含む。既定の行列は登録の後のエラーだけ。`latest_preprocess` なら `details.matrix_note` も）。
旧 ARF 経路（`source="arf"`）は B-2b Task 8 で消えた。MSLipidMapper へは行列の注入ごとの値（生体のみ・欠測は 0）を渡す。
指定した `acyl_chains` のうち分子種の入った鎖スコープがどのパスウェイにも無い鎖は、手順 17 の後で数えて戻り値と manifest の
`missing_chains`・`notes.acyl_chains` に書くだけで、run は捨てない（空スコープは手順 18 の比較から外れる）。
状態変更: `session.pathway.last_run_id` / `last_run_dir`。
R の仮想環境は起動前に確かめる（`_check_environment()`）: `r/renv/library` に MSLipidMapper が無ければ
`MSLM_R_PACKAGES_MISSING`（`reason=not_restored`）、`r/setup.R` が残す目印（復元に使った `renv.lock` の md5、
`renv/library/restored-lock.md5`）が今の `renv.lock` と違えば同じコードで `reason=stale`。どちらも `Rscript r/setup.R` で直る。
描画のブラウザは設定が無ければ Edge / Chrome の標準の置き場所から探す（`browser_path()`）。

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

前提: `pathway_project` の run（`run_id` 省略時は `session.pathway.last_run_id`）。
無ければ `missing_state` → `pathway_project`。`run_id` は `pathway_project` が採番する形
（`pw_YYYYmmdd_HHMMSS_xxxxxx`）以外を `PATHWAY_BAD_REQUEST` で拒む（`reports/pathway` の外を
指させない）。`inline_image` が bool でなければ `PATHWAY_BAD_REQUEST`。状態変更: なし（PDF と PNG を run ディレクトリに書く）。

1. metabolomix/pathway/tools.py  pathway_plot()
2. ├─ metabolomix/pathway/tools.py  _load_manifest()
3. ├─ metabolomix/pathway/bridge.py  run_bridge()
4. │  ├─ [mslipidmapper_api_url あり] metabolomix/pathway/remote.py  run_job()（output.format=pdf。ジョブ失敗は `PATHWAY_RENDER_FAILED`。記録は `remote-job-render.json`）
5. │  └─ [mslipidmapper_api_url なし・サーバ失敗] （R 側。mode=render で `mslm_bridge.R` が `MSLipidMapper::render_pathway_from_config` を呼ぶ）
6. └─ metabolomix/plots/pdf_raster.py  pdf_first_page_png()（PNG を `<pathway>.png` に書く。戻り値は既定で場所と要約の JSON、`inline_image=true` は [JSON, 画像…]）
