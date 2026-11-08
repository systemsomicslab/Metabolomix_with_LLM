# 出力形式: パスウェイ投影（pathway_project / pathway_plot）

## 何の値か

- MSLipidMapper は**正規化しない**。値はこちらの前処理後行列（`normalization` 欄にレシピ）を、
  パスウェイノードに対応するサブクラスごとに `aggregation`（既定 sum）で集計したもの。
  研究室サーバ（`backend` が研究室サーバ）で動いた run の JSON-LD は `normalization` が `none` になる
  （API が前処理の記述を受け付けないため）。前処理の正準は manifest.json の `normalization`。
  サーバは版を返さないので、その run の `upstream_commit` / `package_version` は空（ローカルの
  `r/renv.lock` の版と同じとは限らない）。
- 中間 mzTab-M に載るのは生体試料（role=sample）だけ。サブクラスが空・Unknown・Others の行と、
  手動除外スポット（arf のみ。`dropped_excluded_spots`）は載らない（件数は `handoff`）。さらに**脂質として届かない行**も載らない。
  行はちょうど 1 つの区分に入る（`handoff.kept_by_route` の経路か `dropped_*` の理由。
  全行 = 経路の合計 + `dropped_duplicate_names` 以外の `dropped_*` の合計。arf は `dropped_excluded_spots` も含む）。
  判定の順は `dropped_unannotated` → 経路 `as_is`（と `dropped_unparsable_name`）→ `dropped_isotope_label` → 経路 `name` → 経路 `class_formula`:
  - 経路 `as_is`: Ontology が MS-DIAL の脂質クラス（`LbmClass`。`PC` `EtherPE` `Cer_NS` …。表は
    `metabolomix/pathway/msdial_lipid_classes.py`）で、名前が脂質ショートハンドとして読める行。名前も Ontology もそのまま。
    名前が読めない行は `dropped_unparsable_name`（限定子 `no MS2:` 等は除いて読む）。
  - `dropped_isotope_label`: Ontology が MS-DIAL の脂質クラスでない行のうち、名前に同位体の標識（重水素
    `-d31`・`(d7)`・`-[…-d5]`、`13C`）がある行（内部標準）。経路 `name` / `class_formula` より先に見る。
    長鎖塩基の `d18:1` は標識ではない（`:` が続く）。経路 `as_is` の行には当てない（下の内部標準の注意）。
  - 経路 `name`: Ontology が MS-DIAL の脂質クラスでない行（研究室 MSP の ClassyFire 系の Ontology。
    `Phosphatidylcholines` 等）のうち、名前から脂質として読める行。古い MS-DIAL の書式・英語のクラス名と鎖も
    MS-DIAL 5 の略記と `LbmClass` に直して入れる。ClassyFire の Ontology は使わず、名前と食い違えば名前を採る。
    クラスが 1 つに決まらない行（`Cer 34:1;O2` など）は `dropped_ambiguous_class`。Cer / HexCer / SM 以外の
    頭部で `;O…` の修飾がある名前（`FA 18:1;O`、`PC 16:0_18:2;O`、`FA 18:1;O (9-HODE)` 等）はこの経路で読まず、経路 `class_formula` に回す。
  - 経路 `class_formula`: 名前が読めない行（慣用名・IUPAC 名・ID）のうち、ClassyFire の Ontology が 1 つの
    `LbmClass` に決まり、組成式がそのクラスの基準の形に合う行。組成式から**和の組成**の名前（`TG 54:3` 等。
    鎖ごとの内訳は無い）を作る。作れない行は、この順で最初に当たった理由で落ちる:
    `dropped_ambiguous_ontology`（MS-DIAL で複数のクラスに分かれる。`Ceramides` 等）・`dropped_unmapped_ontology`
    （Ontology が対応表に無い）・`dropped_class_name_conflict`（名前に頭部の語 `cholin` / `ethanolamin` / `serin` /
    `inositol` があり、Ontology のクラスのものと食い違う。PC と PE は C3H6 しか違わず組成式では分からないため。
    PC 系と SM は `cholin`、PE 系と NAE は `ethanolamin`、PS は `serin`、PI は `inositol`、ほかのクラスは語を持たない）・
    `dropped_no_formula`・`dropped_formula_mismatch`。
  - `dropped_unannotated`: 名前か Ontology が空・`null`・Unknown・Others の行と、dataset で 1 つの特徴を
    複数の SML が指す（注釈が `ambiguous`。名前を選ばない）行。
  クラスや組成は推測で埋めない（決まらない行は落として数える）。経路 `name` / `class_formula` の名前は、元の
  `low score: ` / `no MS2: ` を先頭に残す。どの経路も MSP 由来か LBM 由来かは問わず、出力にも書かない。
  **`handoff`（`kept_by_route`・`n_kept`・`dropped_*`）と、落とした行（arf の除外スポットを含む）があるときの
  `notes.handoff` で、使った行と落とした行の数を必ず確かめる。**
  試し（MTBLS6965 の 2 群比較、`msp_then_lbm`）の件数は
  `docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md` の節「MSP の脂質の handoff（2026-11-02）」。
  同名の行は平均強度最大の 1 行だけ残す（`dropped_duplicate_names`。`kept_by_route` は絞る前の数で、
  `n_kept` = その合計 − `dropped_duplicate_names`）——sum がアダクト違いを二重計上しないため。
- arf はレシピ先頭の prop（`handoff.prop`）だけを使う。
- 中間 mzTab-M に書く文字列（試料名・群名・脂質名・Ontology）の `#` は `_` に置き換わる
  （上流の `read.table` は `#` 以降をコメントとして捨てるため）。群名も置換後の形で返る。
  注釈付きの特徴が 2 件未満なら `PATHWAY_BAD_REQUEST`（上流の読み込みが 1 特徴の表で壊れる）。
- **sum は有限値だけを足す**（上流の仕様）。ある試料で一部の分子種が欠損していると、その試料の
  サブクラス合計はそのぶん小さくなる。JSON-LD には試料ごとの有効分子種数（`validFeatureCount`）が
  載るが、ノード表の比較には使っていない。欠損の多いサブクラスの sum の差は割り引いて読む。
- **同位体標識の内部標準**（名前に `(d7)` などを含むもの）も、Ontology が MS-DIAL の脂質クラスなら（経路 `as_is`）
  そのサブクラスに合算される。合算させたくないときは事前に除外する（例: `dataset_exclude`）。Ontology が MS-DIAL の
  脂質クラスでない行（研究室 MSP の行）は `dropped_isotope_label` で落ちる。

## ノード表（nodes_tsv）

| 列 | 意味 |
|---|---|
| node | パスウェイのノード名（＝MS-DIAL Ontology と完全一致したサブクラス） |
| scope | `whole`（サブクラス全体）か、アシル鎖（`18:0` 等。そのノードの鎖サブセット） |
| n_features | そのスコープに入った分子種の数（常に 1 以上。下記） |
| mean_a / mean_b | 各群の集計値の平均（生スケール） |
| log2fc | log2(x+1) 空間の群平均差。**正なら group_b が高い** |
| q | BH-FDR。補正単位は**パスウェイ × スコープ種別**（whole と各鎖は別々） |
| call | up / down（q < 閾値）、ns、NA（片群の値が足りない等で検定不能） |

上流はデータに無いサブクラスのノード（対応表経由の PC/PE/PS など）と、指定した各鎖について
全ノードに**空のスコープ**（分子種 0・値は全部 null）を出す。これは「値の付いたノード」ではないので
行にせず、BH の分母にも入れない。指定した `acyl_chains` のうち、分子種の入った鎖スコープが
どのパスウェイにも無い鎖があれば `PATHWAY_BAD_REQUEST`（`details.missing_chains`）。

## 逆向きに動く隣接ペア（opposite_pairs_tsv）

- **隣接**＝間に接合ノード（PathwayJunction）だけを挟んでつながる 2 つの脂質ノード（`via` が
  direct / junction）。値の付かない化合物ノード（CDP-DAG 等）を挟むものは隣接としない。
- 掲載は**両端とも q < 閾値かつ向きが逆**の組だけ。片方が ns/NA の組は `undeterminable_pairs` に件数のみ。
  片方に値（分子種の入った whole スコープ）が無い組は数えない（判定不能にも入れない）。
- MSLipidMapper のつながりは描画用テンプレートで、反応の主張ではない（上流文書:
  "not automatically curated reactions"）。**観察の提示**として扱い、基質・生成物の解釈は文献で裏を取る。

## run の manifest.json に残る要約

`pathway_project` は戻り値と同じ注意と隣接ペアを `reports/pathway/<run_id>/manifest.json` にも書く
（MCP の時間切れで戻り値が失われても run に残る）。

- `notes`: `adjacency`（隣接の注意）と `handoff`（行を落としたときの注意。落とした行が無ければ null）。
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
  `figures`＝`pathway` `pdf` `png` `width` `height` `mapped_nodes`、`note`）。画像はファイルにしか
  出ない（設計者の判断 2026-11-04）。チャットに図を出すときは `inline_image=true`
  （要約の JSON に続けて PNG を画像で返す。1 枚約 2,400 画像トークン）。
- ノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、ノード表の q とは**別物**。
- ノード内グラフの横軸は group_a, group_b, その他の群の順（上流が config の `plot.x_order` を軸順と配色に使う）。
