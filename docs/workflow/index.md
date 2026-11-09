# ワークフロー: MCP ツールの呼び出し連鎖

各 MCP ツールが、どのファイルのどの関数を、どの順に呼ぶかを記録する。
出力フィールドの**意味**は `docs/output_format/` に、ツールの**外形**（引数・用途）は
`USAGE.md` にある。ここに書くのは**内部の呼び出し順**だけ。

## 読み方

各節は 3 要素からなる。

- **前提** — 呼ぶ前に成立していなければならないセッション状態。満たされない場合に
  どの手順で `MissingState` を返すかを書く。
- **状態変更** — そのツールが `session` に書き込む内容。次のツールの前提になる。
- **呼び出し連鎖** — 番号付きリスト。`└─` の深さがコールスタックの深さに対応する。

呼び出し連鎖には**行番号を書かない**。リファクタで静かに嘘になるため。パスと関数名は
`tests/test_workflow_docs.py` が AST で実在を検証しており、実装から乖離すると
テストが落ちる。

メソッドは `ClassName.method()` と書く。`session`（`AnalysisSession`）の状態は、データセットと
ライブラリがハンドルの登録簿（`core/handles.py`）に、キュレーション・パスウェイが
専用スロットにある（Console ジョブのスロット `current_job_path` は D 段で消えた）。旧パーサ別のスロット（`session.arf` / `arf2` / `pai2` / `eic`）は B-2b Task 8 で消えた。

## 層の構成

- `metabolomix/core/` — FastMCP インスタンス、セッション状態、パス解決、共通ヘルパ
- `metabolomix/msdial/` — MS-DIAL 固有のサイドカー（`*_tags.xml` / `.mddata`）と同定・検証
- `metabolomix/arf/` `arf2/` `pai2/` `dcl/` `eic/` — 形式ごとのパーサ（`reader.py`）。MCP ツールを持つのは `dcl/` だけ
- `metabolomix/analysis/` — 入力形式に依存しない数値処理
- `metabolomix/plots/` — レンダラ中立の payload 組み立て
- `metabolomix/pathway/` — MSLipidMapper 連携（最小 mzTab-M の受け渡し・研究室サーバの REST API とローカル Rscript のブリッジ・JSON-LD 読み取り・隣接比較）
- `metabolomix/tools/` — 形式に紐づかない MCP 公開層（`notes.py` はレポートと objective の読み書き）

## 目次

| 文書 | 対象 | ツール数 |
|---|---|---|
| [dataset.md](dataset.md) | データセット投入の入口 — ファイルの一覧・読み込み（`load_dataset`）・一覧と閉じる・試料の状態と除外 | 6 |
| [dcl.md](dcl.md) | `.dcl` — デコンボリューション済み MS/MS | 3 |
| [plots.md](plots.md) | 図（dataset_plot）と図の保存（save_figure）と payload 契約の比較 | 2 |
| [dataset_analysis.md](dataset_analysis.md) | DatasetState — 前処理・PCA・差次的解析・エクスポート・実験情報訂正・特徴検索 | 7 |
| [pipeline.md](pipeline.md) | 生データフォルダ起点の統括(pipeline) — Console の実行・解析フォルダの検証・受領までと、メソッドファイルを選ぶ・作る 2 ツール | 7 |
| [library.md](library.md) | 参照ライブラリ（`.dbs`/`.msp`）— MS/MS スペクトル照合と対向プロット・ライブラリの検索 | 4 |
| [curation.md](curation.md) | アラインメントのキュレーション — 注釈の一覧確認・機械判別・フラグ記録・候補付け | 5 |
| [pathway.md](pathway.md) | パスウェイ投影（MSLipidMapper）— サブクラスのノード集計・2 群比較・図 | 2 |
| [reanalysis.md](reanalysis.md) | 再解析 — spec で呼ぶツール（ハンドルの一覧・クローズは dataset.md） | 8 |
| [quant.md](quant.md) | 定量 — 内部標準の正規化・極性の統合・特徴の選別（規則プラグインが要る） | 4 |

合計 49 ツール。対象外は 4 ツール——レポート記録系 3
（`record_objective` `write_report` `read_report`）＋
サーバ自身の保守 1（`server_update`。解析の流れに現れず、配布先クローンの更新だけを行う）。
登録ツール総数は 53（49 + 4）。

## 一気通貫の順序

入口は 2 つで、合流点は MS-DIAL の解析フォルダ（2 つは対等）。GUI で解析済みの出力は `list_data_files` →
`load_dataset` で読む（[dataset.md](dataset.md)）。生データフォルダから解析フォルダを作るのは `pipeline_run`
（[pipeline.md](pipeline.md)）で、フォルダに MS-DIAL を動かした形跡が無く、ユーザーが Console を回すことに
同意したときだけ使う。pipeline が回すのは上流だけ——worker が極性ごとに `prepare_input → upstream →
validate_outputs` を要求の順に進め、最後に `receipt`（受領）を書く。統計は計算しない。続きは
`load_dataset(pipeline_path=…, polarity=…)` で GUI の入口と同じ鎖（[dataset_analysis.md](dataset_analysis.md)
の前処理・比較・PCA）へ合流する。

各関数の引数・戻り値のシグネチャ表（個々のツールを手動で順に呼ぶ経路の設計根拠）を書いた
[../superpowers/specs/2026-09-03-end-to-end-pipeline-design.md](../superpowers/specs/2026-09-03-end-to-end-pipeline-design.md)
は、`docs/workflow/` へ昇格させずに過去の設計記録として残す——上流は構成見直しの D 段（2026-11）の pipeline、
下流は B/C/A 段で置き換えた（同文書 §9 の状態の欄）。現行の順序の正準は
`pipeline.md` と `dataset.md` / `dataset_analysis.md` にある。
