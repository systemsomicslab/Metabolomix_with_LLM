# 分子種ごとの図・分子種 PCA・ローディング図・図の保存の一本化 設計

日付: 2026-10-08。設計はユーザー承認済み（会話内、節ごとに承認）。
試作は `C:\Users\yuu18\datasets\20260930_EV\v3_20261008\scripts\`（`species_plots.py` / `pca_loadings.py`）。
リポジトリの外にあるので、本実装はそれを参照するだけで取り込まない。

## 1. 目的

20260930_EV で上司から求められた資料（分子種ごとの割合・強度の図、PCA のローディング）を、
別のデータセットでも LLM が MCP ツールを呼ぶだけで作れるようにする。

- 生物系の議論では「その分子種がどれだけ変わったか」を群ごとの点と棒で見たい。
- PCA はスコアだけでなく、どの分子種が主成分を作っているか（ローディング）を見たい。
- 図を PDF に束ねる機能は作らない（ユーザー判断）。

## 2. 範囲と、公開面の変化

| 区分 | ツール |
|---|---|
| 足す | `arf_plot_species`、`arf_pca_species`、`plot_pca_loadings`、`save_figure` |
| 消す | `save_pca_figure`、`save_volcano_figure`、`save_eic_figure`、`save_group_intensity_figure` |

- ツールの総数は変わらない（4 足して 4 消す）。
- 分子種ごとの図と分子種 PCA は **ARF 経路のみ**（`load_dataset` / `arf_parser` で読んだ `.arf` と兄弟の `.arf2`）。
- ローディング図は ARF（`arf_parser` / `arf_pca_preprocessed`）・分子種 PCA・mzTab（`dataset_pca`）の 3 つの PCA を描く。
- 既存の 4 つの保存ツールは**移行期間を置かずに消す**（ユーザー判断）。下流の Use-LLLM は
  本体がツール名を見ず `annotations`（`may_write_files`）で PNG を拾っているので、直すのはテストの
  ツール名一覧と README だけ。同じ作業の中で直す（§9）。
- 描画 payload の形式名は外部契約の旧名規則に従う: `lipidmix.species_intensity.v1`、`lipidmix.pca_loadings.v1`。

## 3. 共通部品: 項目から分子種（スポット）の一覧を作る

`arf_plot_group_intensity` の中にある解決処理を、MCP 非依存の関数に切り出す（新モジュール
`metabolomix/plots/item_selection.py`）。入力と除外の規則は群別強度と同じ:

- 項目: `.arf2` の Ontology に完全一致すればクラス、しなければ分子種名（`|` 区切りの候補名、`low score:` 等の接頭辞は無視）。`+` で合算。
- 除外（payload の `excluded` に理由別で残す）: 標識内部標準（`(d7)` 等）、`standard_samples` の試料にだけある分子種、
  キュレーションの wrong / redundant（`apply_curation`、既定 True）、最新レビューの likely_wrong（`exclude_auto_likely_wrong`）、
  `arf_exclude` で除いたスポット（manual）。assign は名前・クラスの付け替え。
- MS/MS の裏付け: 照合結果に MS/MS があり、かつ `.arf2` の Name が `no MS2:` / `w/o MS2:` でない。

群別強度との違い:

- 群別強度は項目を 1 点（当たったスポットの合計）にまとめる。新しい関数は**アラインメントスポットごとに展開**して返す。
  同じ名前で付加イオンが違うスポットは別の要素になる。
- `require_msms`（既定 False）を足す。True のとき MS/MS の裏付けがないスポットを `excluded["no_msms"]` に入れて外す。
- 付加イオンの絞り込み（例: DG の [M+Na]+ を外す）は引数にしない。`arf_exclude` でスポットを除く。

`arf_plot_group_intensity` は中身をこの関数に置き換えるだけで、振る舞いと payload（`lipidmix.group_intensity.v1`）を変えない。
既存のテストがそのまま通ることで確かめる。`arf_plot_species` と `arf_pca_species` は同じ関数で解決するので、
同じ指定なら図と PCA が同じ分子種の集合を使う。

## 4. `arf_plot_species`

### 4.1 引数

| 引数 | 既定 | 意味 |
|---|---|---|
| `items` | 必須 | クラス名・分子種名。クラスは配下のスポットに展開する |
| `groups` | 必須 | `arf_plot_group_intensity` と同じトークン規則。QC・ブランクは明示したときだけ入る |
| `value` | `"share"` | `"share"`（割合 %）または `"height"`（PeakHeight） |
| `share_basis` | `None` | 割合の分母にする項目。省略時は `items`。分母にも同じ除外を適用する |
| `low_reliability_samples` | `None` | 白抜きで描き、平均と SD から外す |
| `apply_curation` / `exclude_auto_likely_wrong` / `standard_samples` / `require_msms` | §3 | §3 のとおり |
| `title` / `ncols` / `output` | | `output` は `"image"`（既定）/ `"payload"`。`LIPIDMIX_PLOT_OUTPUT` に従う |

### 4.2 値と描き方

- 1 パネル = 1 スポット。見出しは「分子種名、付加イオン（ID）」。パネルの順はクラス（`items` の順）→ 分子種名。
- `share`: 試料ごとに「スポットの高さ ÷ 分母（`share_basis` を解決したスポットの高さの合計）× 100」。
  棒 = 群平均、ひげ = SD（n ≥ 2 の群だけ）、点 = 各試料。分母が 0 の試料は描かず、説明文に試料名を書く。
- `height`: 対数軸に点。群の平均 ± SD は群別強度と同じ規則（log10 の平均）。0 は軸の下端に置く。
- 点の形: 低信頼の試料は白抜き（平均・SD に入れない）、gap-filled の値は菱形。群の色と並びは `groups` の順。
- パネルは 1 回 40 枚まで。超えたら `{"status": "error"}` で、項目を分けて呼ぶよう促す（画像トークンの抑制）。

### 4.3 戻り値と状態

- 既定: `[説明文, Image(png)]`。説明文は群の n、パネル数、除外の内訳、分母が 0 だった試料。
- `payload`: `lipidmix.species_intensity.v1`。スポットごとに ID・名前・クラス・付加イオン・m/z・RT・MS/MS の有無、
  試料ごとの高さ・割合・gap-filled・低信頼、群ごとの n・平均・SD。分母の定義（`share_basis` を解決したスポット ID）も持つ。
- `session.arf.last_species_plot` に payload と描画引数を残す。失敗した呼び出しの後は `None` にする。

## 5. `arf_pca_species`

### 5.1 引数

`items` / `groups` / `low_reliability_samples` / `apply_curation` / `exclude_auto_likely_wrong` /
`standard_samples` / `require_msms` は §4 と同じ（同じ関数で解決）。加えて:

| 引数 | 既定 | 意味 |
|---|---|---|
| `normalize` | `"none"` | `"total"` は試料ごとに選んだスポットの合計で割り、全試料の合計の中央値を掛ける |
| `log_transform` | `True` | log10(x + 1) |
| `orient_by` | `None` | 群名。その群（計算に使った試料）の平均スコアが正になるよう各主成分の符号をそろえる |
| `title` / `output` | | |

### 5.2 計算

- 対象の試料は `groups` に入る試料。低信頼の試料は主成分の計算に使わず、投影だけする。
- 前処理: 正規化 → log → 計算に使う試料で autoscale（平均 0、SD 1。母標準偏差）。投影する試料にも同じ変換を当てる。
- 計算に使う試料で分散が 0 のスポットは外し、件数を説明文に書く。計算に使う試料は 3 以上（満たさなければ `{"status": "error"}`）。
- `analysis/pca.py` に「計算に使う試料を選べる」関数 `run_pca_fit_subset` を新しく足す。既存の `run_pca` は変えない
  （`arf/reader.py` の再エクスポートとテストのモックを守る）。
- ローディング: 主成分ごとに固有ベクトルの成分（coefficient）と相関 r（= 成分 × 特異値 / √n、n は計算に使った試料数。
  母標準偏差で autoscale したときに成り立つ）。

### 5.3 戻り値と状態

- 既定: `[説明文, Image(png)]`。スコア図（PC1 × PC2、群で色分け、低信頼は白抜き、投影点が枠外なら端に矢印）。
  説明文は寄与率、分子種の数、外したスポット数、PC1 の r の上位・下位 5 件。
- `payload`: 既存の PCA の点列と同じ形（`x` / `y` / `label` / `group`）に、寄与率とローディングの要約を添える。
- `session.arf.last_species_pca` にスコア・ローディング全量・特異値・特徴量（スポット）情報・計算に使った試料・
  scaling・provenance を残す。PCA の保存候補に `source="species"` として加わる（§7）。

## 6. `plot_pca_loadings`

### 6.1 引数

| 引数 | 既定 | 意味 |
|---|---|---|
| `source` | `"auto"` | `"auto"` / `"arf"` / `"mztab"` / `"species"`。`auto` は候補が 2 つ以上なら `AMBIGUOUS_RESULT_SOURCE` で止まる（既存の規則） |
| `result_id` | `None` | 特定の結果を名指しする |
| `pcs` | `[1, 2]` | 描く主成分（1〜3 個） |
| `top_n` | `15` | 主成分ごとに正の上位 N と負の上位 N。`None` は全件 |
| `value` | `"r"` | `"r"` または `"coefficient"` |
| `title` / `output` | | |

### 6.2 規則

- `top_n=None` で特徴量が 60 を超えるときは `{"status": "error"}` で絞るよう促す。
- 全件のときは全パネルの行を最初の主成分の値の順にそろえる。上位 N のときは主成分ごとに選ぶ。
  特徴量が 2N 未満なら正と負の重複を除く。
- r は autoscale した PCA でだけ出せる。scaling が autoscale でない結果で `value="r"` のときは coefficient で描き、説明文に書く。
- 棒の色は Ontology（クラス）。ARF と分子種 PCA は `.arf2`、mzTab は特徴量の注釈から取る。分からなければ灰色。
- 描画は新設の `plots/pca_loadings.py`（MCP 非依存）。3 つの PCA を同じ行の形（特徴量 ID・表示名・クラス・付加イオン・
  主成分ごとの coefficient と r）に変換してから渡す。

### 6.3 各 PCA が残すものの変更

- ARF（`arf_parser` / `arf_pca_preprocessed`）: `_remember_arf_pca_plot` がスコアに加えてローディング全量・特異値・
  特徴量名・計算に使った試料数・scaling を `last_pca_plot` に残す。LLM への戻り値は変えない。上位 N の選び方と
  スポット情報の付与は既存の `get_pca_loading_features` を使う。
- mzTab（`dataset_pca`）: `last_pca` はすでにローディングを持つ。r のために特異値を足す（要約の戻り値には載せない）。
  特徴量 ID は `pp_feature_names`、表示名とクラスは特徴量の注釈から。
- 分子種 PCA: §5.3。

### 6.4 戻り値と状態

既定 `[説明文, Image(png)]`、`payload` は `lipidmix.pca_loadings.v1`（どの PCA か・寄与率・値の種類・主成分ごとの行）。
`session.last_loadings_plot` に残す（どの経路の PCA でも 1 か所）。

## 7. `save_figure`

`save_figure(kind, analysis_id, title=None, source="auto", result_id=None)`。`annotations` は既存の保存ツールと同じ
（`readOnlyHint=False, destructiveHint=False, idempotentHint=True`）。

| `kind` | 読む結果 | 書くもの | `source` / `result_id` |
|---|---|---|---|
| `pca` | ARF / mzTab / 分子種 PCA のスコア | `<id>_pca.png` | 使う（`arf` / `mztab` / `species`） |
| `volcano` | ARF / mzTab の差次的解析 | `<id>_volcano.png` | 使う |
| `eic` | 直近の EIC 描画 | `<id>_eic.png` | 使わない |
| `group_intensity` | `last_group_intensity` | `<id>_group_intensity.png`（300 dpi）＋ `.svg` | 使わない |
| `species` | `last_species_plot` | `<id>_species.png`（300 dpi）＋ `.svg` | 使わない |
| `pca_loadings` | `last_loadings_plot` | `<id>_pca_loadings.png`（300 dpi）＋ `.svg` | 使わない（描いたときに決まっている） |

- 書き出し先は今と同じ `reports/figures/`。`pca` / `volcano` / `eic` の中身（解像度・描画関数）は今の保存ツールと同じ。
- 未知の `kind`、`source` / `result_id` を使わない `kind` への指定は `{"status": "error"}`。
- 対象の結果がなければ `missing_state`。`required_tools` は kind ごと（例 `species` → `["arf_plot_species"]`、
  `pca` → 今の `save_pca_figure` と同じ）。
- 戻り値の文面は今と同じ形（保存した**絶対パス**と、本文に埋め込む `![...](figures/...)`）。Use-LLLM は絶対パスから PNG を拾う。

## 8. エラーの返し方

- 前提の状態がない: `missing_state` の封筒（ARF 未読込・兄弟 `.arf2` なし・PCA なし等）。
- 入力の誤り: 群別強度と同じ `{"status": "error", "message": ...}`（群が当たらない、パネル超過、`top_n` < 1、
  `pcs` の範囲外、計算に使う試料が 3 未満など）。
- 失敗した描画呼び出しの後は、その種類の直近結果を捨てる（古い図を保存させない）。

## 9. テストと文書

### 9.1 テスト（テストを先に書く）

- 共通部品: 群別強度の既存テストがそのまま通る。展開・`require_msms`・除外理由の振り分け。
- `arf_plot_species`: 割合の計算と `share_basis`、低信頼を平均から外す、分母 0 の試料、パネル上限、payload の形、失敗時の状態破棄。
- `arf_pca_species`: 低信頼の試料が計算に入らず投影される、r が計算に使った試料での相関に一致、`orient_by`、分散 0 の除外、試料 3 未満。
- `plot_pca_loadings`: 上位 N の選び方と重複除去、全件時の行の並び、60 超のエラー、r が出せない場合、ARF / mzTab / 分子種の変換、`AMBIGUOUS_RESULT_SOURCE`。
- `save_figure`: 6 種類の保存、未知の kind、使わない `source` の指定、kind ごとの `missing_state`。
- 既存テストの付け替え: 消す 4 ツールを参照するテスト（`test_report_tools` `test_result_output` `test_differential_tools`
  `test_eic_plot` `test_eic_multi_plot` `test_group_intensity_tools` `test_missing_state_contract`）を `save_figure` に移す。
- 腐敗防止テスト: `test_server_registration`（登録数・`annotations`）、`test_tool_annotations`、`test_workflow_docs`、`test_readme_links`。

### 9.2 文書

- `USAGE.md`（ツール一覧と件数）、`docs/output_format/arf.md`（新しい payload 2 つ、保存の記述）、`docs/output_format/eic.md`、
  `docs/workflow/{arf,eic,plots}.md`（新ツールの呼び出し順。`docs/workflow/index.md` の対象範囲も確認）。
- ツールの説明文・エラー文面で消す 4 ツールを案内している箇所（`arf/tools.py` `eic/tools.py` `core/tool_helpers.py`
  `core/session_state.py` `mztab/dataset_state.py` `plots/volcano.py`）。
- vault の流れ図（新しい描画の枝と保存ツールの一本化）、`docs/HISTRY.md`、`docs/task.md`。
- Use-LLLM: `tests/test_tool_catalog_completeness.py`・`tests/test_approval_snapshot.py`・`tests/test_general_agent.py` の
  ツール名、`README.md` の 1 行。Use-LLLM 側のテストも通す。

## 10. やらないこと

- 図を PDF に束ねる機能。
- mzTab 経路の分子種ごとの図・分子種 PCA。
- 付加イオンでの絞り込み引数（`arf_exclude` で足りる）。
- 古い保存ツールの移行期間（ユーザー判断で同時に消す）。
