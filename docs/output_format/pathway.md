# 出力形式: パスウェイ投影（pathway_project / pathway_plot）

## 何の値か

- MSLipidMapper は**正規化しない**。値はこちらの前処理後行列（`normalization` 欄にレシピ）を、
  パスウェイノードに対応するサブクラスごとに `aggregation`（既定 sum）で集計したもの。
- 中間 mzTab-M に載るのは生体試料（role=sample）だけ。サブクラスが空・Unknown・Others の行と、
  手動除外スポットは載らない（件数は `handoff`）。同名の行は平均強度最大の 1 行だけ残す
  （`dropped_duplicate_names`）——sum がアダクト違いを二重計上しないため。
- arf はレシピ先頭の prop（`handoff.prop`）だけを使う。
- 中間 mzTab-M に書く文字列（試料名・群名・脂質名・Ontology）の `#` は `_` に置き換わる
  （上流の `read.table` は `#` 以降をコメントとして捨てるため）。群名も置換後の形で返る。
  注釈付きの特徴が 2 件未満なら `PATHWAY_BAD_REQUEST`（上流の読み込みが 1 特徴の表で壊れる）。
- **sum は有限値だけを足す**（上流の仕様）。ある試料で一部の分子種が欠損していると、その試料の
  サブクラス合計はそのぶん小さくなる。JSON-LD には試料ごとの有効分子種数（`validFeatureCount`）が
  載るが、ノード表の比較には使っていない。欠損の多いサブクラスの sum の差は割り引いて読む。
- **同位体標識の内部標準**（名前に `(d7)` などを含むもの）も、実在の Ontology が付いていれば
  そのサブクラスに合算される。合算させたくないときは事前に除外する（例: `arf_exclude`）。

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

## 対応付けの見出し

- `mapped_nodes` = 値の付いたノード（分子種の入った whole スコープを持つ脂質ノード）の数 /
  脂質ノード（`PathwayLipidNode`）の数。化合物ノード（`PathwayEntity`。CDP-DAG 等）と接合ノードは
  分母に入らない。
- `unmatched_subclasses` = データにあるがどのノードにも乗らなかったサブクラス
  （ノード名との完全一致のみ。`EtherPC` などの表記揺れはここに出る）。

## 図（pathway_plot）

- 各パスウェイの PDF（`reports/pathway/<run_id>/<pathway>.pdf`）の 1 ページ目を PNG にしたもの。
- ノード内の p 値ラベルは MSLipidMapper（Shiny と同じ計算）のもので、ノード表の q とは**別物**。
- ノード内グラフの横軸は group_a, group_b, その他の群の順（上流が config の `plot.x_order` を軸順と配色に使う）。
