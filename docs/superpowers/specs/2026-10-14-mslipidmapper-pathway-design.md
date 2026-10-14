# MSLipidMapper を通したパスウェイ投影ツール —— 設計

日付: 2026-10-14 / 状態: 設計承認済み（実装前）
前提調査: [2026-09-04-mslipidmapper-transferable-assets.md](2026-09-04-mslipidmapper-transferable-assets.md)
（候補 9 = 相互運用、候補 10 = 隣接不一致検出）。本設計は両者を、上流の JSON-LD 出力
（MSLipidMapper `097a88b`、2026-10-13 "Add RDF JSON-LD pathway measurements"）の上に組み直したもの。

## 1. 目的と合意事項

- MS-DIAL 出力（ARF）と mzTab-M（DatasetState）を入力に、MSLipidMapper のパスウェイ
  （remodeling / ceramide / global）へ脂質サブクラス（＋任意でアシル鎖サブセット）を投影する。
- **主な受け手は LLM**: ノード表と「逆向きに動く隣接ペア」を小さな TSV で返す。
  **副産物として** JSON-LD と PDF をファイルに残し、**PDF はチャットに画像（PNG）で表示する**。
- 実行方式は**ローカル `Rscript` のサブプロセス**（データをマシン外へ出さない）。
  上流 REST API（Docker / 研究室サーバ）は将来の別バックエンド候補で、本設計の対象外。
- `pipeline_run` には組み込まない（補助経路のため）。

### 上流について確認済みの事実

- JSON-LD は**出力**。入力は MS-DIAL Alignment CSV または mzTab-M。
- `mztab_pathway_measurements(x, network, ...)` は mzTab-M パスか `SummarizedExperiment` を取り、
  サンプル別の値をパスウェイノードへ写した JSON-LD（schema `2.0.0`）を返す。**FC・p 値・正規化はしない**。
- `render_pathway_from_config(config, ...)` が PDF を描く。YAML＋入力ファイル前提で SE は受けない。
  ヘッドレス Chromium（`MSLIPIDMAPPER_BROWSER` で指定可。このマシンは Edge）と
  `unpkg.com` の Cytoscape.js が要る。
- `write_mztab_m(se)` が上流自身の読み込みと往復する最小 mzTab-M を書く
  （MTD `assay[n]` / `study_variable[n]-assay_refs`、SML `chemical_name` `abundance_assay[n]` `opt_global_Ontology`）。
- 上流の mzTab 取り込みは群を MTD `study_variable` から、サブクラスを `opt_global_Ontology`
  （無ければ rgoslin）から取る。手元の Console 出力は `opt_global_Ontology` を持つが、
  `study_variable` がサンプル 1 個ずつになっている例がある。
- ノードとサブクラスの対応は空白除去後の完全一致。
- ライセンス MIT。

## 2. 全体構成とデータの流れ

```
session.arf（前処理後行列・sample_meta・除外・Ontology）
  または DatasetState（前処理後行列・試料メタ・SML の opt_global_Ontology ← 新たに読む）
        │  Python: 上流 write_mztab_m と同形の最小 mzTab-M を書く（群 = こちらで設定した群）
        ▼
<解析フォルダ>/reports/pathway/<run_id>/input.mztab
        │  subprocess: Rscript metabolomix/pathway/mslm_bridge.R <mode> <args.json>
        │    ├ measure: mztab_pathway_measurements() → <pathway>-measurements.json
        │    └ render : render_pathway_from_config()  → <pathway>.pdf
        ▼
Python: JSON-LD を読む → ノード別・スコープ別のサンプル値 → Welch＋BH（既存）→ 隣接の向き判定
        ▼
LLM へ: 見出し＋ノード表＋逆向き隣接ペア（TSV）／チャットへ: PDF 1 ページ目 → PNG（PyMuPDF）
```

中間形式を TSV でなく mzTab-M にする理由: 測定（`mztab_pathway_measurements`）と描画
（`render_pathway_from_config`、`input.format: mztab-m`）の**両方を上流の公開入口のまま**呼べ、
上流の内部関数（`.render_pathway_cytoscape_pdf` 等）に依存しないため。群・除外・前処理の
決定権はこちらに残る。

### 配置

| 場所 | 中身 |
|---|---|
| `metabolomix/pathway/handoff.py` | ARF / DatasetState → 最小 mzTab-M の書き出し（純ロジック） |
| `metabolomix/pathway/bridge.py` | 設定解決・Rscript 起動・タイムアウト・終了判定 |
| `metabolomix/pathway/mslm_bridge.R` | 起動直後に renv の仮想環境 `r/` を有効化し、公開関数（`MSLipidMapper::` で始まるもの）だけを呼んで `check` / `measure` / `render` を実行する薄いラッパ（上流の R ファイルを `source()` しない） |
| `metabolomix/pathway/jsonld.py` | JSON-LD の読み取り（ノード・接続・スコープ・サンプル値） |
| `metabolomix/pathway/compare.py` | 2 群比較（既存 Welch＋BH を呼ぶ）と隣接の向き判定 |
| `metabolomix/pathway/tools.py` | MCP ツール `pathway_project` / `pathway_plot` |

- 上流は **R パッケージとして renv の仮想環境 `r/` にインストールする**。
  上流の R ファイルの多くがトップレベルで `library()` を呼ぶため、PDF 描画に要る関数群
  （`aggregate_to_class_se` `get_plot_fun` `.render_pathway_cytoscape_pdf` ほか）を
  部分的に `source()` することはできない（上流 REST ワーカーも `library(MSLipidMapper)` で動く）。
  入れ元は GitHub `systemsomicslab/MSLipidMapper@<commit>` で、`r/renv.lock` に固定する。
  利用者は `Rscript r/setup.R` の 1 コマンドで同じ環境を復元する。R 本体のユーザー
  ライブラリには入れない。ブリッジは起動直後に `r/renv/activate.R` を読み、
  `MSLipidMapper::` の公開関数だけを呼ぶ。上流のコミットはインストール済みパッケージの
  `RemoteSha` から読む（ローカルクローンは不要）。
  設定は `[mslipidmapper] rscript`（必須）・`browser`（描画時）。環境変数
  `MSLIPIDMAPPER_RSCRIPT` `MSLIPIDMAPPER_BROWSER`。ブラウザは描画 YAML の `pathway.browser` で渡す。
- 成果物は解析フォルダ配下 `reports/pathway/<run_id>/`。リポジトリ側には書かない。

## 3. ツール

### `pathway_project`（LLM 向け）

```
pathway_project(
  source="auto",            # auto | arf | dataset（auto: 直近に読んだ方、両方あれば arf）
  group_a=..., group_b=..., # source の流儀（arf: Class ID / 因子トークン、dataset: サンプル名リスト）
  group_a_label="group_a", group_b_label="group_b",   # dataset のときだけ使う（arf は指定文字列）
  pathways=["remodeling","ceramide","global"],
  acyl_chains=[],           # 指定時は whole に加えて鎖サブセットのスコープを作る
  aggregation="sum",        # sum | mean | median（上流の集計）
  q_threshold=0.05,
)
```

- mzTab には除外されていない生体試料を**全群**書く。比較はこちら側で group_a 対 group_b。
- 戻り値（`json_payload` の見出し＋TSV。`structured_output=False`）:
  1. 見出し: `run_id`、source、群ごとのサンプル数、前処理レシピ（MSLipidMapper は正規化しないので
     実施内容を申告）、パスウェイごとの対応付いたノード数 / 全ノード数、どのノードにも乗らなかった
     サブクラス、§4 で落とした行の件数。
  2. ノード表（パスウェイごと、値の付いたノードのみ）:
     `node scope n_features mean_a mean_b log2fc q call`。`scope` は `whole` または鎖、
     `call` は `up` / `down` / `ns` / `NA`。検定は log2(x+1) 空間の Welch、BH は**パスウェイ×スコープごと**。
  3. 逆向きに動く隣接ペア: **隣接**＝間に `PathwayJunction` だけを挟んで `PathwayConnection` で
     つながる 2 つの `PathwayLipidNode`。値の付かない `PathwayEntity`（CDP-DAG 等）を挟めば隣接としない。
     掲載条件は**両端とも q < 閾値かつ向きが逆**。片方が `ns` / `NA` の組は判定不能として件数のみ。
     エッジ方向は上流の `sourceEntity`→`targetEntity` をそのまま示すが、上流は反応とは主張しない
     （"not automatically curated reactions"）ので、生化学的解釈は LLM に委ね、出力にもその旨を書く。
- 保存: `input.mztab`、`<pathway>-measurements.json`、`nodes.tsv`、`manifest.json`
  （引数・群・source・`dataset_id`・上流クローンのコミット・schemaVersion）。

### `pathway_plot`（人間向け、チャットに表示）

```
pathway_plot(run_id=None, pathways=None, plot_type="dot")   # run_id 省略時は直近
```

- その run の `input.mztab` で上流 PDF 描画を呼び、各 PDF の 1 ページ目を PyMuPDF で PNG 化。
- 戻り値: 短い説明文（PDF パス・対応付いたノード数）＋ `Image` をパスウェイごと。長辺約 1600px
  （1 枚約 2,400 画像トークン）。`pathways` で絞れる。
- 図中のノード内 p 値ラベルは上流 Shiny と同じ計算で、`pathway_project` の q とは別物だと説明文に明記。
- run が無ければ `missing_state`（`required_tools: ["pathway_project"]`）。

### セッション

`session.pathway` スロット（新設）は直近の `run_id` と run ディレクトリだけを持つ。正準は
`manifest.json` で、再起動後も `run_id` 指定で `pathway_plot` できる。他パーサのスロットに触らない。

## 4. 中間 mzTab-M の中身

上流 `write_mztab_m` と同形。正規化レシピは `MTD comment[1] Normalization: <レシピ>`。

**前提（無ければ `missing_state`）**
- `session.assay_kind == "lipid"`。`metabolite` / `unknown` は `required_tools: ["record_objective"]`。
- 前処理後の行列: arf は `session.arf.feature_matrix`（`arf_preprocess` 後）、dataset は DatasetState の
  前処理後行列。既存の差次的解析と同じ前提。行列は線形スケール（`preprocess` は対数化しない）。

**サンプル（列）**
- 役割 `sample` のみ。QC・ブランク・手動除外は書かない。
- 群: 比較に入るサンプルは**比較の指定名**（arf は `group_a`/`group_b` の文字列、dataset は
  `group_a_label`/`group_b_label`）を群ラベルにする（因子トークンでプールした群が図でも 1 群に並ぶ）。
  比較に入らない生体試料は、arf は `sample_meta[...]["group"]`、dataset は実験情報シート
  （`sample_metadata_rows` の `group`、`include` が偽のものは除く）があればその群で書き、
  無ければ書かない。

**特徴（行）**
- サブクラス: arf は**同じアラインメントの兄弟 `.arf2`** のスポット（`MasterAlignmentID` で結合）の
  `Ontology` と `Name`（ARF のスポットは Ontology を持たない）。兄弟 `.arf2` が無ければ
  `missing_state("sibling_arf2", ["load_dataset"])`（既存 `_require_arf_with_arf2` と同じ）。
  前処理行列が複数の prop（`height` と `area` 等）を持つときはレシピ `props` の先頭だけを使う。
  dataset は SML の `opt_global_Ontology` を新たに読み、
  SMF↔SML を結んで特徴メタデータに `ontology` を足す（差次的エクスポートの列は不変、
  `CONTRACT_VERSION` は据え置き）。
- サブクラスが空・`Unknown`・`Others` の行と手動除外スポットは書かない（件数は見出しで報告）。
- 同じ脂質名が複数行（アダクト違い等）にあるときは、sum の二重計上を避けるため
  **名前ごとに平均強度最大の 1 行だけ残す**。落とした行数と名前の例を見出しに出す。
- 補完済みの値はそのまま、未補完の欠損は `null`（上流で MissingValue）。

## 5. エラー処理

`mcp_errors.py` に `pathway_error()` と `PATHWAY_ERROR_CODES` を足す（`console_error` と同型）。

| code | いつ | 対処 |
|---|---|---|
| `MSLM_NOT_CONFIGURED` | repo / rscript 未設定・不在 | 設定ファイルに書く（ユーザー） |
| `MSLM_R_PACKAGES_MISSING` | `check` で R パッケージ不足（測定用・描画用を分けて判定） | 不足名とインストールコマンドを `details` に。自動インストールしない |
| `MSLM_RUN_FAILED` | Rscript が非ゼロ終了 | stderr 末尾 20 行、全文は run ディレクトリの `bridge.log` |
| `MSLM_TIMEOUT` | 既定 300 秒超過 | プロセスツリーごと終了（`core/process_control.py`） |
| `MSLM_OUTPUT_INVALID` | JSON-LD 不在・破損・`schemaVersion != "2.0.0"` | 上流の仕様変更の検知 |
| `PATHWAY_RENDER_FAILED` | ブラウザ不在・unpkg 不達 | `pathway_plot` 固有。`pathway_project` の結果は有効のまま |
| `PATHWAY_BAD_REQUEST` | 未知のパスウェイ名、群のサンプル < 2、指定鎖がデータに無い | 引数を直す |

子プロセスは引数リストで起動（シェル不使用）、cwd は run ディレクトリ。再試行しない。

## 6. テスト

R が無くても `pytest tests` は全部緑。fixture はテスト自身が作る（上流の例ファイルに依存しない）。

- `handoff`: ARF / DatasetState の fixture から、役割絞り込み・除外・`Unknown` 落とし・同名重複の
  間引き・`null`・MTD の群割り当てを検証。
- `jsonld`: テスト内で組んだ小さな `@graph` で、接合ノード経由の隣接、`PathwayEntity` で切れる隣接、
  MissingValue を検証。
- `compare`: 逆向きペア抽出、判定不能件数、BH の単位（パスウェイ×スコープ）。
- PDF→PNG: fitz で作った PDF で寸法上限を検証。
- `bridge`: 偽 Rscript（Python スクリプト）で正常・非ゼロ・タイムアウト・出力不在。
- 実 R＋クローンの結合テストは `LIPIDMIX_TEST_MSLM=1` のときだけ走る。
- 登録数（`test_server_registration`、CLAUDE.md 冒頭）、`USAGE.md`、`docs/workflow/`、
  新規 `docs/output_format/pathway.md` とリソーストピック `pathway` を更新。
- vault 側の流れ図に新しい解析枝を追記。

## 7. 実データでの検証（実装の最後）

データ: `C:\Users\yuu18\datasets\20260930_EV\neg_wiff` / `pos_wiff`（GUI 解析済み ARF。
`cell_C` 対 `cell_PlnA` 各 3 例。ほかに EV・blank・MeOH）。極性ごとに別 run（両極性の統合は対象外）。

1. R 依存を入れ、`mslm_bridge.R check` が通ること。
2. 各極性で `pathway_project` が走り、ノード対応付けの件数を記録する
   （2026-09-04 spec §4-2「`.cyjs` のクラス語彙と MS-DIAL Ontology の一致」をここで実測して閉じる）。
3. `pathway_plot` で PNG がチャットに出ること。

## 8. 対象外

- 上流 REST API バックエンド、`pipeline_run` への組み込み、両極性の統合、独自 `.cyjs` の指定、
  遺伝子発現の重ね合わせ、metabolite アッセイへの適用。
