# ワークフロー: MCP ツールの呼び出し連鎖

各 MCP ツールが、どのファイルのどの関数を、どの順に呼ぶか。出力フィールドの**意味**は `docs/output_format/`、
ツールの**外形**（引数・用途）は `USAGE.md`。

## 読み方

各節は 3 要素からなる。

- **前提** — 呼ぶ前に成り立っているべきセッション状態と、満たされないときにどの手順で `MissingState` を返すか。
- **状態変更** — そのツールが `session` に書く内容（次のツールの前提）。
- **呼び出し連鎖** — 番号付きリスト。`└─` の深さがコールスタックの深さ。

呼び出し連鎖に行番号は書かない。パスと関数名は `tests/test_workflow_docs.py` が AST で実在を検証する。

メソッドは `ClassName.method()` と書く。`session`（`AnalysisSession`）の状態は、データセットとライブラリが
ハンドルの登録簿（`core/handles.py`）に、キュレーション・パスウェイが専用スロットにある。

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
| [dataset.md](dataset.md) | データセット投入の入口 — ファイルの一覧・読み込み（`load_dataset`）と受付番号（`dataset_load_status`・長い道具の `job_status`）・一覧と閉じる・試料の状態と除外 | 8 |
| [dcl.md](dcl.md) | `.dcl` — デコンボリューション済み MS/MS | 3 |
| [plots.md](plots.md) | 図（dataset_plot）と図の保存（save_figure）と payload 契約の比較 | 2 |
| [dataset_analysis.md](dataset_analysis.md) | DatasetState — 前処理・PCA・差次的解析・エクスポート・実験情報訂正・特徴検索 | 8 |
| [pipeline.md](pipeline.md) | 生データフォルダ起点の統括(pipeline) — Console の実行・解析フォルダの検証・受領までと、メソッドファイルを選ぶ・作る 2 ツール・試料の表を ISA / mwTab から作る 1 ツール | 8 |
| [library.md](library.md) | 参照ライブラリ（`.dbs`/`.msp`）— MS/MS スペクトル照合と対向プロット・ライブラリの検索 | 4 |
| [curation.md](curation.md) | アラインメントのキュレーション — 注釈の一覧確認・機械判別・フラグ記録・候補付け | 5 |
| [pathway.md](pathway.md) | パスウェイ投影（MSLipidMapper）— サブクラスのノード集計・2 群比較・図 | 2 |
| [reanalysis.md](reanalysis.md) | 再解析 — spec で呼ぶツール（ハンドルの一覧・クローズは dataset.md） | 9 |
| [quant.md](quant.md) | 定量 — 内部標準の正規化・極性の統合・特徴の選別（規則プラグインが要る） | 4 |

合計 53 ツール。対象外は 5 ツール——レポート記録系 3（`record_objective` `write_report` `read_report`）＋
サーバ自身の保守 1（`server_update`。配布先クローンの更新だけ）＋出力の意味の文書を読む 1（`docs_read`。リソースと同じ本文を返す）。
登録ツール総数は 58（53 + 5）。

## 一気通貫の順序

入口は 2 つで対等、合流点は MS-DIAL の解析フォルダ。GUI で解析済みの出力は `list_data_files` → `load_dataset` で読む
（[dataset.md](dataset.md)）。MS-DIAL を動かした形跡の無い生データフォルダから解析フォルダを作るのは `pipeline_run`（[pipeline.md](pipeline.md)）で、
MS-DIAL Console を極性ごとに回し（極性ごとに数十分）、run のフォルダを書く。pipeline が回すのは上流だけ: worker が極性ごとに
`prepare_input → upstream → validate_outputs` を要求の順に進め、最後に `receipt`（受領）を書く。統計は計算しない。続きは
`load_dataset(pipeline_path=…, polarity=…)` で GUI の入口と同じ鎖（[dataset_analysis.md](dataset_analysis.md) の前処理・比較・PCA）へ合流する。
