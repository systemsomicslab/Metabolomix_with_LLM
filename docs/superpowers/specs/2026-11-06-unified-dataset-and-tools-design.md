# B 統一データセットとツールの統合 設計

- 日付: 2026-11-06
- 関連: [全体](2026-11-06-overhaul-overview-design.md)（目的・流れ・横断の決定はそちら）、
  C [統計エンジン](2026-11-06-statistics-engine-design.md)、A [最初の概観](2026-11-06-dataset-overview-design.md)、
  D [pipeline](2026-11-06-pipeline-upstream-design.md)。msdial-workflow（`C:\Users\yuu18\msdial-workflow`）。
- 状態: 利用者と合意した設計（2026-11-06）。未実装。順番は B → C → A → D の 1 番目。

## 1. 目的

データの持ち方を `DatasetState` の 1 つにそろえ、その上でツールを作業の単位に並べ直す。
あわせて、msdial-workflow が持つ定量の工程（内部標準による正規化・POS と NEG の統合・付加体の選別）を
MCP の鎖に入れ、データセットが自分を同定したライブラリを覚えるようにする。

## 2. 目標・非目標・反目標

目標:

- MS-DIAL の解析フォルダ（`.arf` / `.arf2` / `.dcl` / `.EIC.aef` / `_tags.xml` / `.mddata`）から `DatasetState` を作るアダプタ。
  mzTab-M のアダプタ（今の `dataset_load`）と並ぶ 2 本目で、解析フォルダの側を主とする。
- 対話のセッションとハンドルを 1 つの仕組みにまとめる。
- 旧経路（`session.arf` / `.arf2` / `.pai2` / `.eic` のスロットと、それを読むツール）を消す。
- 定量の工程（§6）とライブラリの記憶（§7）。
- ツールの対応表（§8）どおりに公開面を作り直す。

非目標:

- 統計の計算そのもの（C）。B は C のツールが乗る土台と、旧統計ツールを消すところまで。
- 最初の概観の中身（A）。
- pipeline の工程の整理（D）。B は `load_dataset` が pipeline の run を読めることだけを保つ。
- LBM（`.lbm2`）を下流の照合で読むこと（§7.3）。

反目標:

- 解析フォルダから作った `DatasetState` と、同じアラインメントの mzTab-M から作った `DatasetState` が、
  特徴 ID・試料・値で食い違う（§4.4 の照合で縛る）。
- gap-fill の値を、実測と区別できない形で次の工程に渡す。
- 閉じたデータセットや別のデータセットの状態を、黙って読む（今の `_holder.resolve_holder` の性質を保つ）。
- 研究室の参照データ（内部標準の係数・付加体の採用表・参照ライブラリ）の場所や中身を、追跡対象の文書・テスト・戻り値に書く
  （既存の外部流出禁止の規則。係数の値や採用表の行も同じ扱いにする）。

## 3. 今の状態（実装者が読む場所）

- `DatasetState`: `metabolomix/mztab/dataset_state.py`。`source_format` の欄があり、mzTab-M 以外の入力も受けられる作り。
  検出状態 `detected_mask`、試料ごとの m/z `per_sample_mz`、SME 候補 `feature_candidates`、SML のラベル `feature_annotations`、
  結果の登録 `results`、指紋 `result_state.dataset_fingerprint` を持つ。
- mzTab-M の読み込み: `metabolomix/mztab/loading.py`、隣の `.arf` から検出状態を取り込む `metabolomix/mztab/evidence.py`。
- ハンドル: `metabolomix/core/handles.py`（登録簿）、`metabolomix/tools/_holder.py`（`resolve_holder`）、
  `metabolomix/core/handle_dispatch.py`（ハンドルを渡した呼び出しだけをワーカースレッドで回す）。
- 対話のセッション: `metabolomix/core/session_state.py` の `session`（`arf` / `arf2` / `pai2` / `eic` / `dataset` / `library`）。
- ARF のパーサ: `metabolomix/arf/reader.py`、`metabolomix/arf2/reader.py`、`metabolomix/dcl/reader.py`、`metabolomix/eic/reader.py`。
  MessagePack の Key 番号の正解表は `docs/schema/*.md`（推測で直さない）。
- 同定の標準化（GOSLIN・RefMet/LIPID MAPS・MSI レベル）: `arf2_annotate_identities` の実体（`metabolomix/msdial/` 配下）。
- キュレーション: `metabolomix/curation/`。判断は `<アラインメントのフォルダ>\curation\flags.jsonl` に、`.arf2` の sha256 と
  MasterAlignmentID の組で記録する。
- 特徴 ID: ARF の MasterAlignmentID と mzTab-M の SMF_ID は同じ番号の空間（`analysis/evidence_export.py` の
  `link_alignment_spectra` と `_id_space_check` が m/z・RT で確かめている）。

## 4. 解析フォルダのアダプタ

### 4.1 入力の選び方

今の `load_dataset` の規則を引き継ぐ: 同じフォルダに `DriftSpots.arf` と `PeakProperties.arf` があれば `PeakProperties.arf`、
複数のバッチがあればファイル名の `AlignmentResult_<timestamp>` で最新を選ぶ。選んだものと、選ばなかった候補を戻り値に書く。

### 4.2 `DatasetState` の欄への写し方

| 欄 | 出どころ |
|---|---|
| `feature_matrix` | `.arf` の試料ごとの PeakHeight（既定）。面積を選べるようにするかは今の `arf_parser` の `props` に合わせる |
| `detected_mask` | `.arf` の各セルの gap-fill の印（検出 = gap-fill でない）。mzTab 経路の evidence sidecar と同じ意味 |
| `per_sample_mz` | `.arf` の検出セルの m/z |
| `feature_ids` | MasterAlignmentID（mzTab-M の SMF_ID と同じ番号） |
| `feature_metadata` | `.arf2` の名前・オントロジー・付加体・m/z・RT・InChIKey・Comment・MS/MS 照合の有無とスコア。同定の標準化（GOSLIN・MSI レベル）の結果もここに足す |
| 同定の証拠の水準 | mzTab 経路の `sme` / `sml` / `none` と同じ 3 つに写す。MS/MS で照合された代表 → `sme` 相当、MS1 だけの注釈 → `sml` 相当、無し → `none`。`.arf2` のどの欄で決めるかは `docs/schema/` で確かめて決める |
| 試料 | `.arf2` / `.mddata` の試料名・Class ID（群）・ファイルの種別（Sample / QC / Blank / Standard）・注入順・バッチ。種別があるときは、名前から役割を推す今の `detect_sample_roles` より優先する |
| 極性 | `.arf2` / `.mddata` のイオンモード |
| ライブラリ | §7 |
| タグとキュレーション | `_tags.xml`（MS-DIAL GUI の Confirmed / Misannotation など）と `flags.jsonl` |
| `source_format` | `"msdial_folder"` |

`.dcl`（代表 MS/MS）と `.EIC.aef` は値として取り込まず、場所と sha256 を `artifact_paths` に持つ（スペクトル・EIC のツールが引く）。

### 4.3 mzTab-M のアダプタ

今の `dataset_load` の中身（mzTab-M・`job_path`・`pipeline_path`、出所の裏取り `source_verification`、未完了の run を探索専用にする
`exploratory_only`、検出状態の取り込み）をそのまま `load_dataset` の 1 経路として残す。
pipeline の run を読んだとき、run の中の古い差次的結果は引き継がない（D §6）。

`load_dataset(pipeline_path=…)` と `load_dataset(job_path=…)` の既定は、run（またはジョブ）が作った**解析フォルダのアダプタ**
（§4.2）にする（どちらの入口からも同じ鎖を通すため。私の判断）。`source="mztab"` を渡すと、今と同じ mzTab-M のアダプタで読む。
再解析研究が既存の spec（C §12）で同じ数字を出し直すときは、こちらを明示して使う。mzTab-M を直接渡したときは mzTab-M のアダプタ。

### 4.4 2 つのアダプタの照合

Console の run は、同じアラインメントについて `.arf` 一式と mzTab-M の両方を出す。この組で、2 つのアダプタが作った
`DatasetState` を突き合わせる検査を作る: 特徴 ID の集合、試料の集合と並び、値（許容幅つき）、検出状態、同定の名前。
リポジトリのテストは合成の fixture で行い（テストの規約: fixture はテスト自身が作る）、実データでの照合は
reanalysis の run（`C:\Users\yuu18\datasets\reanalysis\`）で 1 度回して `docs/HISTRY.md` に結果を書く。

## 5. データセットの登録の一本化

- 読み込んだデータセットはすべてハンドルの登録簿に入る。対話のセッションは「直近のデータセット」を指すだけになる。
- データセットを扱うツールは省略可能な引数 `dataset_handle` を持つ。省略すると直近のデータセットを使う。
  閉じた・存在しないハンドルは `DATASET_HANDLE_NOT_FOUND` で止め、直近のデータセットへ黙って落とさない。
- 引数の名前 `dataset_handle` とハンドルの ID の形（`dh_…`）は今のまま（再解析の手順書の書き換えを減らす。実装者の判断で可）。
- 並行の規則は今の `handle_dispatch` を全ツールへ広げる: 同じデータセットへの呼び出しは 1 本ずつ、別のデータセットどうしは並行。
  対話のセッションも同じ規則で動く（イベントループ上で同期に動かす特例をやめてよい。実装者の判断）。
- `session.arf` / `.arf2` / `.pai2` / `.eic` のスロットは消す。図を保存するツールなどが読んでいた中間データは、
  データセットの結果の登録（`results`）か、データセットに付く図の置き場へ移す。
- ライブラリは今のとおりデータセットごとに持つ（ハンドルの `library`）。

## 6. 定量の工程（msdial-workflow と共有）

### 6.1 共有の形

msdial-workflow（Tk の GUI と CLI。LLM を使わない人に配る）は残す。定量の規則は 1 か所に置き、GUI と MCP の両方がそれを使う。

- 規則の置き場所は msdial-workflow の側。MCP がその版を固定して依存する（逆向きにすると、オフラインで配る GUI の zip に
  FastMCP や scikit-learn などの MCP の依存が入るため）。
- msdial-workflow の側で、規則を入出力から切り離した純粋な関数として切り出す。入力は「特徴の表（名前・オントロジー・付加体・
  Comment・MS/MS 照合の有無・試料ごとの値）と参照表」で、MS-DIAL のテキスト出力（GUI の今の入力）からも
  `DatasetState`（MCP）からも作れる形にする。今の実体は次のあたり:
  - 内部標準の正規化: `src/msdial_normalizer/core.py` の `normalize`、`service.py` の `preview` / `reference_context`
  - 試料の照合と統合: `src/msdial_workflow/tables.py` の `strip_polarity` / `match_samples` / `merge_tables`
  - 付加体の選別と除外規則: `src/msdial_workflow/tables.py` の `filter_data` と `exclusion_reason`
- **依存の取り方（利用者の判断）**: msdial-workflow を GitHub に置き、ms-data-parser は commit を固定して依存する
  （`pyproject.toml` の git の依存と `uv lock`。`requirements*.lock.txt` の書き出しも合わせる。CLAUDE.md の依存の手順）。
  msdial-workflow には今リモートが無く、公開ライセンスも未指定。置き場所（systemsomicslab の下か）・公開か非公開か・
  ライセンスは、置くときに利用者が決める（push は都度確認）。非公開なら、研究室の各機の `server_update` が依存を取れるか
  （GitHub の認証）を確かめる。研究用の参照 CSV と provenance は今のとおり Git に入れない。
- 研究室の参照データはリポジトリの外に置き、中身を戻り値に全量で出さない。
  - **付加体の採用表**は実験をまたいで共通なので、`lipidmix.local.toml` で指す（例 `[quant] adducts = ...`）。
  - **内部標準の参照表（係数）は実験ごとに変わる**（係数 = 濃度 × 溶媒量 ÷ (バイアル中の細胞数 ÷ 10^6)。通常版と ST 版、
    濃度方式もある）。そのため設定ファイルに既定を置かず、`dataset_quant_preview` / `dataset_quant_apply` の呼び出しごとに
    指定させる（私の判断。黙って違う実験の係数を使い、単位だけ正しく見える事故を避ける）。指定の形は 2 つ:
    参照表のファイルの場所か、内部標準ごとの濃度・溶媒量・細胞数（msdial-workflow の GUI の「内部標準濃度を設定」と同じ計算を、
    共有する規則の関数で行う）。後者で作った表は解析フォルダの `reports/` の記録の横に保存する。どちらも、使った表の
    sha256 と単位を記録に書く。

### 6.2 工程

鎖の順番は全体の設計書 §2 のとおり（キュレーション → 極性ごとの正規化 → 統合 → 選別）。各工程は記録を書く（全体 §5）。

1. **内部標準の確認（止まり場所）**: 極性ごとに、参照表のクラスごとの内部標準の候補（スポット ID・付加体・m/z・RT）と係数、
   欠けそうな値を示す。候補が 1 つに決まらないクラスがあれば、`needs_input` の封筒（候補つき）で止め、利用者の選択
   （msdial-workflow の `bindings` と同じ形）を受けて進む。
2. **内部標準による正規化**: 極性ごとに `対象の値 / 同じ試料の内部標準の値 × 係数`。単位は参照表のもの（既定 pmol / 10^6 cells）。
   内部標準の値が 0・負・欠測なら、その試料の値は欠測にして理由を残す（中央値などで埋めない）。参照表に無いクラスは
   `Any others` の行で計算する（msdial-workflow と同じ）。
   - gap-fill の扱い: 対象か内部標準のどちらかが gap-fill のセルは、正規化後の値に「gap-fill 由来」の印を引き継ぐ
     （msdial-workflow はテキスト出力を読むためこれができず、gap-fill を実測として使っている。MCP では `detected_mask` で区別できる）。
   - 正規化した後のデータセットには「正規化済み」「方法 = 内部標準」「単位」を持たせる。C の前処理はこれを読んで、
     TIC や PQN を既定で重ねない（C §4）。
3. **POS と NEG の統合**: 試料の照合は msdial-workflow と同じ規則（名前の完全一致 → `pos` / `neg` の語を除いて一致 →
   対応表）。食い違えば止まる（止まり場所）。統合したデータセットは新しいハンドルになり、特徴 ID に極性を付ける
   （例 `pos:123`）。スペクトル・EIC を引くときに元の極性のフォルダへ戻れるよう、特徴ごとに元のデータセットを持つ。
4. **付加体の選別と除外**: 参照表（クラス × 付加体 × 極性 × 採用）で選び、クラスの別名（`Cer_NS` / `Cer-NS` など）を解く。
   除外規則（`Unknown`・`w/o`・`RIKEN`・重水素標識の `(d`・`low score`・Comment が `x` で始まる・MS/MS 照合なし）と
   内部標準そのものの行を除く。MCP のキュレーションの判断（`wrong` / `redundant` は除く、`assign` は同定を置き換える）も
   ここで当てる（今の書き出しの `apply_curation` と同じ意味）。除いた行と理由はすべて記録に残す。

工程を飛ばす（正規化しない・選別しない）ことも選べる。飛ばしたことは記録とデータセットの状態に残る。

### 6.3 ツール

工程ごとにツールを分ける（名前は案。実装者の判断で変えてよい）: `dataset_quant_preview`（1）、`dataset_quant_apply`（2）、
`dataset_merge_polarities`（3）、`dataset_select_features`（4）。

## 7. ライブラリの記憶と設定の一本化

### 7.1 データセットが自分を同定したライブラリを覚える

`load_dataset` が、そのデータセットの同定に使ったライブラリを `DatasetState` に記録する:
ファイル名・sha256・種類（`dbs` / `msp` / `lbm`）・出どころ。

| 入口 | 出どころ |
|---|---|
| GUI で解析したフォルダ | MS-DIAL がフォルダに残す `*_Loaded.msp2.dbs`（読み込んだライブラリの写しと照合条件）。複数あって決まらなければ「不明」とし、候補を戻り値に書く |
| Console の run | 実効メソッド（`effective-method.txt`）と run の入力計画（MSP・LBM の sha256 を固定済み） |
| mzTab-M を直接読んだ | 「不明」 |

記録した場所のファイルが消えている・sha256 が違うときは、別のファイルを黙って使わず「不一致」と記録する。

### 7.2 下流は覚えたライブラリを使う

`library_load` は引数なしで、データセットが覚えたライブラリを読む。`library_match_feature`・`library_msms_qualify`・
根拠の束（`dataset_claim_evidence`）も既定でそれを使う。`file_path` を明示すれば別のライブラリも使える（記録に残す）。
今の `library="auto"`（極性の設定から引く）と `ion_mode` の設定から引く経路は、下流からは消す
（極性ごとの MSP の設定は、Console にどれを使わせるかを決めるためだけに残る）。

### 7.3 LBM は後回し

LBM（`.lbm2`）は今の照合の道具で読めない。

- LBM だけで同定したデータセット（`library_mode` が `auto` の脂質の run など）は、下流の照合を `not_evaluable` と理由で返す。
- `msp_then_lbm` のデータセットが覚えるライブラリは MSP（読める方）にし、**すべての特徴をその MSP と照合する**（今の
  `library="auto"` と同じ結果になる。再解析の既存の spec の数字を変えないため。私の判断）。LBM 側が代表になった特徴には、
  照合の結果に「同定は LBM による」の印を付ける。特徴ごとにどちらが代表かを `.arf2` から読めるかは、実装者が `docs/schema/` で確かめる。
  読めなければ印は付けず、そのことを文書に書く。
- 記録した sha256 とファイルが合わないときの「不一致」は、今は黙って使っていた場合だけ挙動が変わる（そのとき今の数字はもともと
  再現できていない）。

### 7.4 設定と環境変数

- 外部資産の場所は `lipidmix.local.toml` だけで指す。資産の場所を上書きする環境変数を消す:
  `MSDIAL_EXE`・`MSDIAL_LBM`・`MSDIAL_MSP_POS`・`MSDIAL_MSP_NEG`・`MSDIAL_LIBRARY_LBM`・`MSLIPIDMAPPER_RSCRIPT`・
  `MSLIPIDMAPPER_BROWSER`・`MSLIPIDMAPPER_API_URL`・`MSLIPIDMAPPER_API_TOKEN`（`metabolomix/core/user_config.py` の対応表）。
- LBM の設定 `[msdial] lbm` と `[library] lbm` を 1 つ（`[library] lbm`）にまとめる。Console の LBM の自動解決
  （メソッドの宣言 → MsdialWorkbench のビルド生成物 → 設定 → 実行体のフォルダ）の順は変えず、設定のキーだけを 1 つにする。
  `msp_then_lbm` が `[library] lbm` だけを読む規則は残す。
- 知識の蓄積を消すので `LIPIDMIX_KNOWLEDGE_DIR`・`LIPIDMIX_PLAYBOOK_DIR` も消える。
- 置き場所の環境変数（`LIPIDMIX_DATA_DIR`・`LIPIDMIX_REPORTS_DIR`・`LIPIDMIX_LIBRARY_CACHE_DIR`・
  `LIPIDMIX_MAX_DATASET_HANDLES`・`LIPIDMIX_CAVEAT_MODE`）は設定ファイルへ移す（私の判断。利用者が確認する）。
- 環境変数として残すのは、起動ごとに変わるものだけ: `LIPIDMIX_CONFIG`（テストと worktree）、`LIPIDMIX_TRANSPORT` /
  `LIPIDMIX_HOST` / `LIPIDMIX_PORT`、`LIPIDMIX_PLOT_OUTPUT`（Use-LLLM がクライアントごとに置く）。
- テストの conftest（環境変数を消し `LIPIDMIX_CONFIG` を存在しないパスへ向ける）と `scripts/new-worktree.sh`
  （worktree の `.mcp.json` に置き場所の環境変数を書く）を合わせて直す。worktree は main ツリーの蓄積状態を向く、という性質は保つ。

## 8. ツールの対応表

旧名は一度に切り替える（移行期間なし）。C・A・D のツールは名前だけ載せる（中身はそれぞれの設計書）。

| 作業 | 新しい公開面 | 吸収・削除する旧ツール |
|---|---|---|
| 入口 | `list_data_files`、`load_dataset`（解析フォルダ・mzTab-M・`job_path`・`pipeline_path` のどれでも受け、ハンドルを返し、最後に A の概観を 1 回呼ぶ）、`dataset_list`、`dataset_close` | `dataset_load`、`dataset_handle_open`・`_list`・`_close`（改名）、`arf_parser`、`arf2_parser`、`pai2_parser`、`eic_parser`、`dcl_parser` |
| 試料 | `dataset_status`（群・タグ・試料の役割・因子トークンでの試料の検索を含む）、`dataset_set_sample_metadata`（個体 ID `subject_id` を足す。C §5）、`dataset_exclude`（試料と特徴の除外。可逆） | `arf_list_classes`、`arf_list_tags`、`arf_list_sample_roles`、`sample_search`、`arf_exclude` |
| 定量 | `dataset_quant_preview`、`dataset_quant_apply`、`dataset_merge_polarities`、`dataset_select_features`（§6） | — |
| キュレーション | `curation_review`、`curation_suggest`、`curation_submit`、`curation_flags`、`curation_view_data`（`DatasetState` の上で動かす） | — |
| 前処理・統計（C） | `dataset_preprocess`、`dataset_compare`、`dataset_pca`、`dataset_analysis_run` | `arf_preprocess`、`arf_pca_preprocessed`、`arf_differential`、`dataset_differential`、`dataset_statistic`、`dataset_build_matrix`（旧 `dataset_preprocess` と `dataset_pca` は中身を C で作り直す） |
| 概観（A） | `dataset_overview` | — |
| 書き出し | `dataset_export_differential`（契約の列は `analysis/export_contract.py`）、`dataset_export_evidence_bundle` | `arf_export_differential`、`arf_export_evidence_bundle` |
| 特徴と根拠 | `dataset_feature_search`（`eic_search` の範囲・ID・強度順の検索を吸収）、`dataset_feature_relations`、`dataset_claim_evidence`、`dataset_alignment_qc`、`dataset_mass_accuracy`、`dataset_homolog_check`、`dataset_novel_candidates` | `eic_search`、`arf2_annotate_identities`（標準化は読み込みに入れ、結果は特徴の表の列で見せる） |
| スペクトルと照合 | `dcl_spot_spectrum`、`dcl_diagnostic_ions`、`dcl_find_msms`、`library_load`、`library_match_feature`、`library_msms_qualify`、`library_plot_mirror`、`library_search` | `pai2_inspect_peak`、`verify_peak_annotation` |
| 描画 | `dataset_plot`（`kind` で volcano・PCA のスコア・ローディング・群ごとの強度・分子種・分子種の PCA・EIC・複数物質の EIC を切り替える）、`save_figure` | `arf_plot_volcano`、`arf_plot_group_intensity`、`arf_plot_species`、`arf_pca_species`、`plot_pca_loadings`、`eic_plot_chromatograms`、`eic_plot_compounds` |
| パスウェイ | `pathway_project`（2 群の比較は C のエンジンを通す）、`pathway_plot` | — |
| 記録 | `record_objective`（`update_objective` を吸収。`assay_kind` の設定に要る）、`write_report`、`read_report`（ID を省けば一覧） | `update_objective`、`list_reports` |
| 文献と知識 | — | `paper_search`、`ingest_stage`、`ingest_review_queue`、`ingest_promote`、`ingest_reject`、`knowledge_coverage`、`log_search` |
| Console・pipeline（D） | `pipeline_plan`・`run`・`status`・`resume`・`cancel`、`console_method_candidates`、`console_method_template` | `console_plan`、`console_run`、`console_status`、`console_prepare_input`、`console_cleanup`、`job_list` |
| 保守 | `server_update` | — |

新しい公開面は 53 ツールになる。MCP リソースは `lipidmix://docs/output-format` とトピック別のテンプレートだけを残し、
`lipidmix://knowledge/*`・`lipidmix://playbook/*` を消す。

描画ツールを 1 本にまとめるときも、今の決まりは保つ: 既定は PNG と短い説明、`output="payload"` で点列、`LIPIDMIX_PLOT_OUTPUT`、
形式名 `lipidmix.<名前>.vN`（Use-LLLM が読む）、`structured_output=False`、図は座標を返さずサーバで描く。

## 9. 消すもの

- `session.arf` / `.arf2` / `.pai2` / `.eic` と、それだけを読むコード（`metabolomix/arf/tools.py` のツール群、
  `metabolomix/pai2/tools.py`、`metabolomix/eic/tools.py` の要約、`metabolomix/arf2/tools.py`）。パーサ（`reader.py`）は
  アダプタが使うので残す。`metabolomix/arf/reader.py` の `run_pca` の再エクスポート（テストのモックのための束縛）は、
  それを差し替えていたテストと一緒に整理する。
- 文献と知識: `metabolomix/corpus/`、`metabolomix/tools/` の該当ツール、`knowledge/` `playbook/` の配信と
  `LIPIDMIX_KNOWLEDGE_DIR` / `LIPIDMIX_PLAYBOOK_DIR`。追跡されている `knowledge/` の種ノートと `playbook/` のファイルは、
  リポジトリの外（Obsidian の vault など。移す先は移すときに利用者に確かめる）へ移してから消す（利用者の判断）。
- 旧ツールの文書: `docs/output_format/` の `arf` `arf2` `pai2` `eic` のトピックは、`DatasetState` の側の説明へまとめ直す
  （行の粒度・gap-fill・脂質名の文法・必須の注意は失わない）。

## 10. 文書とテスト

- `USAGE.md` を新しい公開面で書き直す。`docs/workflow/` を書き直す（`tests/test_workflow_docs.py` が実登録と突き合わせる）。
  CLAUDE.md 冒頭の規模表記・構成の節を直す（`tests/test_readme_links.py`）。`tests/test_server_registration.py` の登録数。
- 腐敗防止のテストは、消したツールの名前が文書に残っていないことも縛る。
- アダプタの照合（§4.4）、定量の工程（§6）、ライブラリの記憶（§7）、ハンドルの一本化（§5）のテストを、合成の fixture で作る。
- 定量の規則は、msdial-workflow の既存のテスト（人工データ）と同じ入力で MCP 側が同じ値を出すことを確かめる。

## 11. 実装者に任せること

- ツールと引数の正確な名前（`dataset_handle` を保つかを含む）、戻り値の欄の名前。
- 工程の記録のスキーマ。
- 統合後の特徴 ID の書式。
- 描画ツールの `kind` の名前と、各 `kind` の引数の切り分け。
- 旧スロットを消す順番（ただし、旧と新が並んだまま止まらないこと）。
