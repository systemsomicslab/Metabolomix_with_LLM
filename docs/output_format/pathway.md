# 出力形式: パスウェイ投影（pathway_project / pathway_plot）

## 何の値か

- MSLipidMapper は**正規化しない**。値はこちらの前処理後行列（`normalization` 欄にレシピ）を、パスウェイノードに対応するサブクラスごとに `aggregation`（既定 sum）で集計したもの。
  研究室サーバ（`backend` が研究室サーバ）で動いた run の JSON-LD は `normalization` が `none`（API が前処理の記述を受け付けない）。前処理の正準は manifest.json の `normalization`。
  サーバは版を返さないので、その run の `upstream_commit` / `package_version` は空（ローカルの `r/renv.lock` の版と同じとは限らない）。
- **ノードの値（とその和）はクラスの総量ではない**。入るのは前処理の足切りを通った特徴のうち、下の受け渡しで脂質として残した行だけ（落とした行は `handoff` の `dropped_*`、同名は 1 行）。戻り値の `notes.values` もそう書く。
- 中間 mzTab-M に載るのは生体試料（role=sample）だけ。サブクラスが空・Unknown・Others の行と、**脂質として届かない行**は載らない（件数は `handoff`）。
  行はちょうど 1 つの区分に入る（`handoff.kept_by_route` の経路か `dropped_*` の理由。全行 = 経路の合計 + `dropped_duplicate_names` 以外の `dropped_*` の合計）。
  判定の順は `dropped_unannotated` → 経路 `as_is`（と `dropped_unparsable_name`）→ `dropped_isotope_label` → 経路 `name` → 経路 `class_formula`
  で、経路が決まった行は最後に名前の頭を確かめる（`dropped_name_head`）:
  - 経路 `as_is`: Ontology が MS-DIAL の脂質クラス（`LbmClass`。`PC` `EtherPE` `Cer_NS` …。表は
    `metabolomix/pathway/msdial_lipid_classes.py`）で、名前が GOSLIN で読める行（`normalize_lipid_name` の
    `parse_ok`。GOSLIN は略記のほか `(9Z)-octadecenoic acid` のような慣用名の一部も読む）。名前も Ontology もそのまま。
    名前が読めない行は `dropped_unparsable_name`（限定子 `no MS2:` 等は除いて読む）。
    同位体標識の内部標準（名前に `(d7)` など）もここに入れば、そのサブクラスに合算される（合算させたくなければ事前に `dataset_exclude` などで除く）。
  - `dropped_isotope_label`: Ontology が MS-DIAL の脂質クラスでない行のうち、名前に同位体の標識（重水素
    `-d31`・`(d7)`・`-[…-d5]`、`13C`）がある行（内部標準）。経路 `name` / `class_formula` より先に見る。
    長鎖塩基の `d18:1` は標識ではない（`:` が続く）。経路 `as_is` の行には当てない。
  - 経路 `name`: Ontology が MS-DIAL の脂質クラスでない行（研究室 MSP の ClassyFire 系の Ontology。
    `Phosphatidylcholines` 等）のうち、名前から脂質として読める行。古い MS-DIAL の書式・英語のクラス名と鎖も
    MS-DIAL 5 の略記と `LbmClass` に直して入れる。ClassyFire の Ontology は使わず、名前と食い違えば名前を採る。
    クラスが 1 つに決まらない行（`Cer 34:1;O2` など）は `dropped_ambiguous_class`。Cer / HexCer / SM 以外の
    頭部で酸素を足す修飾（`;O…` `;OH` `;2OH` `;(2OH)` `;Ep` `;oxo` 等）がある名前（`FA 18:1;O`、`PC 16:0_18:2;O`、
    `PC 16:0_20:4;Ep`、`FA 18:1;O (9-HODE)` 等）はこの経路で読まず、経路 `class_formula` に回す。
    名前の後ろの空白を含む注記は切るが、`;` の直後に空白の無い修飾は残す（`SM 18:1;O2/16:0 (x)` → `SM 18:1;O2/16:0`）。
    注記があっても行き先は変わらない（`Cer 18:1;O2/16:0;(2OH) (Cer AS)` は注記の無い名前と同じく `dropped_ambiguous_class`）。
  - 経路 `class_formula`: 名前が読めない行（慣用名・IUPAC 名・ID）のうち、ClassyFire の Ontology が 1 つの
    `LbmClass` に決まり、組成式がそのクラスの基準の形に合う行。組成式から**和の組成**の名前（`TG 54:3` 等。
    鎖ごとの内訳は無い）を作る。作れない行は、この順で最初に当たった理由で落ちる:
    `dropped_ambiguous_ontology`（MS-DIAL で複数のクラスに分かれる。`Ceramides` 等）・`dropped_unmapped_ontology`
    （Ontology が対応表に無い）・`dropped_class_name_conflict`（名前に頭部の語 `cholin` / `ethanolamin` / `serin` /
    `inositol` があり、Ontology のクラスのものと食い違う。PC と PE は組成式では C3H6 しか違わない。
    PC 系と SM は `cholin`、PE 系と NAE は `ethanolamin`、PS は `serin`、PI は `inositol`、ほかのクラスは語を持たない）・
    `dropped_no_formula`・`dropped_formula_mismatch`。
  - `dropped_unannotated`: 名前か Ontology が空・`null`・Unknown・Others の行と、dataset で 1 つの特徴を
    複数の SML が指す（注釈が `ambiguous`。名前を選ばない）行。
  - `dropped_name_head`: どの経路でも、残す名前（mzTab-M に書く形。前後の空白を除く）の頭の 1 文字が上流の
    クラスの読み取り（`R/lipid_chain_parser_se.R` の `^[A-Za-z0-9\-+;_]+`）に合わない行。`(` や `[` で始まる名前
    （経路 `as_is` の `(9Z)-octadecenoic acid` / FA など）は上流でクラスが NA になり、計測全体が落ちる。
  クラスや組成は推測で埋めない（決まらない行は落として数える）。経路 `name` / `class_formula` の名前は、元の
  `low score: ` / `no MS2: ` を先頭に残す。どの経路も MSP 由来か LBM 由来かは問わず、出力にも書かない。
  **`handoff`（`kept_by_route`・`n_kept`・`dropped_*`）と、落とした行か同名の重複があるときの
  `notes.handoff` で、使った行と落とした行の数を確かめる。** 件数は理由ごとで、どのクラス（ノード）から
  何行落ちたかの内訳は返さない。
  同名の行は平均強度最大の 1 行だけ残す（sum がアダクト違いを二重計上しないため。`dropped_duplicate_names` は落とした行の数。
  `kept_by_route` は絞る前の数で、`n_kept` = その合計 − `dropped_duplicate_names`）。
- 中間 mzTab-M に書く文字列（試料名・群名・脂質名・Ontology）の `#` は `_` に置き換わる
  （上流の `read.table` が `#` 以降を捨てるため）。群名も置換後の形で返る。
  注釈付きの特徴が 2 件未満なら `PATHWAY_BAD_REQUEST`（上流の読み込みが 1 特徴の表で壊れる）。
- **sum は有限値だけを足す**（上流の仕様）。ある試料で一部の分子種が欠損していると、その試料の
  サブクラス合計はそのぶん小さくなる。JSON-LD には試料ごとの有効分子種数（`validFeatureCount`）が
  載るが、ノード表の比較には使っていない。欠損の多いサブクラスの sum の差は割り引いて読む。

## ノード表（nodes_tsv）

| 列 | 意味 |
|---|---|
| node | パスウェイのノード名（＝MS-DIAL Ontology と完全一致したサブクラス） |
| scope | `whole`（サブクラス全体）か、アシル鎖（`18:0` 等。そのノードの鎖サブセット） |
| n_features | そのスコープに入った分子種の数（常に 1 以上。下記） |
| mean_a / mean_b | 各群の単位（`matrix_id` の行列の units。反復をまとめた行列なら個体、まとめない行列なら注入）の `2**log2`（スコープの集計値を log2 にしたもの。個体にまとめた単位は注入の log2 の平均）の算術平均（エンジンの定義） |
| log2fc | log2 空間の群平均差（統計エンジンの `log2_mean_difference`）。**正なら group_b が高い** |
| q | BH-FDR。補正単位は**パスウェイ × スコープ種別**（whole と各鎖は別々） |
| call | `verdict` が `changed`（q ≤ 閾値）なら符号で up / down、検定できなければ NA、それ以外は ns |
| p / ci_low / ci_high | 選んだ検定（`method`）の p と、log2fc の 95% 区間（`mannwhitney` は区間なし） |
| verdict | changed / unchanged / inconclusive（エンジンの判定。区間が同等の範囲に入れば unchanged） |

検定は統計エンジンの `method`（`dataset_compare` と同じ値: `moderated_t` 既定 / `welch` / `student` /
`mannwhitney`。使った検定は戻り値の `method` と manifest の `test`）。スコープの値（正の有限値）を log2 にして比べ、
0 以下・欠測は欠測。**比べる単位は `matrix_id` の行列の units に従う**: `collapse_replicates=True`（既定）の行列で
試料の表に `subject_id` がある注入だけが log2 の平均で 1 個体にまとまり、`collapse_replicates=False` の行列なら注入ごと。
群ごとの単位の数は戻り値・manifest の `n_units`（`groups` は注入の数）で、各群 2 単位以上が要る。
moderated t は、あるスコープ種別のスコープ（ノード）が 3 未満だと事前分布を推定できず縮小が効かない（プールした分散の t と同じ）。

MSLipidMapper へ渡す mzTab-M の値は行列の注入ごとの線形値で、欠測・0 以下・非有限は 0 で書く（MSLipidMapper の約束。エンジンも 0 以下を欠測と読む）。
**注意**: 上流の集計は 0 を値として扱うので、`aggregation=mean|median` では欠測の 0 がサブクラスの値を引き下げる（`sum` なら足さないのと同じ）。欠測の多いサブクラスの mean / median の差は割り引いて読む。比較側（ノード表）は集計後の 0 以下を欠測として除く。

`matrix_id` を省いたときの行列の選び方は `dataset_compare` / `dataset_pca` と同じ（共通核の「状態不足の伝え方」。除外・足切り・設定はその前処理のもの）。どれを使ったかは戻り値と manifest の
`matrix_id` と `matrix_source`（`matrix_id`＝名指し / `latest_preprocess`＝直近の `dataset_preprocess` / `default`＝既定の前処理）。
行列を決めたあとにエラーになったときは `details.matrix_id` と `matrix_source` に残る
（`latest_preprocess` なら `details.matrix_note` も。受け渡しの検査——サブクラス注釈の付いた特徴が 2 未満・群の大きさ——も
使った行列の足切りで決まるので同じく添える。既定の行列は登録の後のエラーだけ）。
使った行列が drift=auto でドリフト補正を当てなかったものなら、戻り値と manifest の `warnings` に `drift_not_applied_by_default`（共通核の `needs_input` の節。当てなかった行列でなければ `warnings` は無い）。補正するなら `dataset_preprocess(drift_mode="on")` の `matrix_id` を渡す。
既定の行列は登録され、前処理の記録が `reports/stats/` に 1 件書かれる（manifest の `preprocess_record_id`。同じ行列の 2 回目以降は記録を足さない。
直近の `dataset_preprocess` の行列ならその前処理の記録）。返った `matrix_id` は dataset_compare / dataset_pca にそのまま渡せる。

上流はデータに無いサブクラスのノード（対応表経由の PC/PE/PS など）と、指定した各鎖について
全ノードに**空のスコープ**（分子種 0・値は全部 null）を出す。これは行にせず、BH の分母にも入れない。
指定した `acyl_chains` のうち、分子種の入った鎖スコープがどのパスウェイにも無い鎖があっても run は返す（その鎖のスコープが行にならないだけで、ほかのスコープの
比較は鎖を指定しなかった run と同じ値）。無かった鎖は戻り値と manifest の `missing_chains` と `notes.acyl_chains` に出る
（綴りの誤りなら正しい綴りで呼び直す）。

## 逆向きに動く隣接ペア（opposite_pairs_tsv）

- **隣接**＝間に接合ノード（PathwayJunction）だけを挟んでつながる 2 つの脂質ノード（`via` が
  direct / junction）。値の付かない化合物ノード（CDP-DAG 等）を挟むものは隣接としない。
- 掲載は**両端とも q ≤ 閾値（verdict が changed）かつ向きが逆**の組だけ。片方が ns/NA の組は `undeterminable_pairs` に件数のみ。
  片方に値（分子種の入った whole スコープ）が無い組は数えない（判定不能にも入れない）。
- MSLipidMapper のつながりは描画用テンプレートで、反応の主張ではない（上流も自動で精査した反応ではないと書く）。観察の提示として扱う。

## run の manifest.json に残る要約

`pathway_project` は戻り値と同じ注意と隣接ペアを `reports/pathway/<run_id>/manifest.json` にも書く
（MCP の時間切れで戻り値が失われても run に残る）。

- `notes`: `adjacency`（隣接の注意）・`handoff`（行を落としたか同名を絞ったときの注意。無ければ null）・
  `values`（値の意味。クラスの総量ではない旨）・`assay_kind`（下記。lipid なら null）・`acyl_chains`（データに無かった
  鎖があるときだけ）。戻り値の `notes` も同じ（null のものは載せず、`figure` を足す）。
- `test`（使った検定）と `n_units`（群ごとの単位の数）。`matrix_id` と `matrix_source`。
  ドリフト補正を当てなかった行列なら `warnings`。
  `acyl_chains` は指定した鎖、`missing_chains` はそのうちデータに無かった鎖（あるときだけ。戻り値にも同じ欄）。
- `opposite_pairs`: パスウェイ名 → `{"rows": [...], "undeterminable": 件数}`。`rows` は
  `opposite_pairs_tsv` と同じ欄（`node_1` `call_1` `log2fc_1` `node_2` `call_2` `log2fc_2` `via`）の辞書。
- 非有限の数値（NaN・inf）は null で書く。

## 対応付けの見出し

- `mapped_nodes` = 値の付いたノード（分子種の入った whole スコープを持つ脂質ノード）の数 /
  脂質ノード（`PathwayLipidNode`）の数。化合物ノード（`PathwayEntity`。CDP-DAG 等）と接合ノードは
  分母に入らない。
- `unmatched_subclasses` = データにあるがどのノードにも乗らなかったサブクラス
  （ノード名との完全一致のみ。`EtherPC` などの表記揺れはここに出る）。

## 図（pathway_plot）

- 各パスウェイの PDF（`reports/pathway/<run_id>/<pathway>.pdf`）と、その 1 ページ目の PNG
  （`<pathway>.png`）を run のディレクトリに書く。
- **既定の戻り値は図のファイルの場所と要約だけ**（`run_id` `group_a` `group_b` `backend` と、
  `figures`＝`pathway` `pdf` `png` `width` `height` `mapped_nodes`、`note`）。画像はファイルにしか出ない。
  `inline_image=true` なら要約の JSON に続けて PNG を画像で返す（1 枚約 2,400 画像トークン）。
- ノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、ノード表の q とは**別物**。
- ノード内グラフの横軸は group_a, group_b, その他の群の順（上流が config の `plot.x_order` を軸順と配色に使う）。

## assay_kind が lipid でないとき

`pathway_project` は `assay_kind`（`record_objective(assay_kind=...)`）を問わず動く。MSLipidMapper は
クラス（Ontology）と脂質の略記の名前を前提にするので、`assay_kind` が `lipid` でないときは
`notes.assay_kind` にその前提を書く。脂質として読めない行は上の受け渡しで落ち、`handoff` の `dropped_*` に数える。

## 失敗したとき

- MSLipidMapper を呼んでからの失敗（ブリッジの `MSLM_*`、JSON-LD を読めない `MSLM_OUTPUT_INVALID`、
  `pathway_plot` の `PATHWAY_RENDER_FAILED`——PDF が出ない・画像にできないときも）の封筒は `details.run_dir`
  （その run のディレクトリ）を持つ。呼ぶ前に断る `pathway_plot` の manifest.json の読み取りの失敗は `details.manifest`。
- run のディレクトリには、サーバの送受信の記録がジョブの種類ごとに残る: `remote-job-measure.json`（`pathway_project` の
  計測のジョブ。run に 1 つ）と `remote-job-render.json`（`pathway_plot` の描画のジョブ。直近の 1 回だけ）。ローカルの
  `bridge.log` は直近のローカルの R の呼び出しのもの（`pathway_plot` が R で描けば計測の分は上書きされる）。
- 研究室サーバのジョブが失敗したとき、サーバは理由を固定文でしか返さない（本文は上流の worker.log にだけあり、API では取れない）。
  - **ローカルの R へ回って成功したとき**: `pathway_project` / `pathway_plot` の戻り値の `backend` は 1 行の文で、
    ローカルの R で動いたことと、使えなかったサーバの api_url・`job <job_id>`（投入まで進んだときだけ）・コードと
    message を書く。戻り値に `remote_fallback` の鍵は無く、記録の場所（`record`）も出ない。`pathway_project` は
    manifest の `remote_fallback`（`api_url`・`code`・`message` と、あれば `job_id`・`record`。サーバで動いた・
    サーバを設定していないなら null）に残し、その `record` は計測の記録 `<run_dir>/remote-job-measure.json` を指す
    （描画では上書きされない）。`pathway_plot` は manifest を書き直さないので、描画のジョブの記録はいつも
    `<run_dir>/remote-job-render.json`（`run_dir` は `pathway_project` の戻り値の `run_dir`。`pathway_plot` の `figures[].pdf` の置かれたフォルダでもある）。
  - **失敗の封筒**: `details.remote`（ローカルの R も使えなかったとき）か `details.remote_fallback`（ローカルの R が
    動いて失敗したとき）に `job_id`（投入まで進んだときだけ）と `record`（そのジョブの `remote-job-<mode>.json` の
    場所）が残り、`details.run_dir` も付く。
- `MSLM_R_PACKAGES_MISSING` は R の仮想環境に MSLipidMapper が無い（未復元・古い・公開関数が無い）。
  `MSLM_R_PACKAGE_LOAD_FAILED` は入っているのに R が読み込めない（message に R の `loadNamespace` の本文）。
  後者はローカルの R を「動かして失敗した」ものとして扱い、サーバが先に失敗していればその理由は `details.remote_fallback`。
