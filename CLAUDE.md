# CLAUDE.md

このファイルは、このリポジトリで作業する Claude 向けの地図。**全体スキャンの代わりに読む**。

## このリポジトリは何か

MS-DIAL（LC-MS リピドミクス／メタボロミクス解析ソフト）の出力を読み、LLM から使える形で
公開する **MCP サーバ `ms-data-parser`**。ツール 72・リソース 5・リソーステンプレート 3。
GUI の外では読めない圧縮 MessagePack と独自バイナリを解し、解析そのものをサーバ側で回して
要約だけを返す（生の行列を LLM の文脈に載せない）。

**読む形式** — `.arf`（アライン後・サンプル別ピーク。サンプル間比較の主データ源）/
`.arf2`（スポット代表カタログと注釈）/ `.pai2`（個別測定のピーク）/ `.dcl`（デコンボリュー
ション済み MS/MS。MessagePack ではない独自バイナリ）/ `.EIC.aef`（クロマトグラム）/
**mzTab-M 2.0**（標準交換形式。Console 経路が産み、`DatasetState` の正準入力になる）。
サイドカー（`*_tags.xml` `.mddata` `.mdproject`）はタグと試料クラスの取得に読む。
生の測定ファイル（`.wiff`/`.raw` 等）は自分では解さず、MS-DIAL Console への入力として扱う。

**やること** — (a) 前処理/QC（正規化・ブランク減算・QC RSD 足切り・ドリフト補正）・PCA・
2 群差次的解析（Welch + BH-FDR）、(b) 同定の標準化と裏取り（GOSLIN 正規化・RefMet/LIPID MAPS
ID・MSI レベル・`.dcl` の実スペクトル照合）、(c) EIC の検索と描画、(d) **生データ起点の統括**
——MS-DIAL Console を計画・実行し、mzTab-M を経て前処理・PCA・比較・品質レポートまでを
永続 run として自動で進める（`pipeline_run`。中断・訂正をまたいで冪等）、(e) LC–MS メタボロ
ミクス v2（内部標準比・feature binding・注入ごとの測定証拠・固定母集団 QC）。**v2 は合成入力と
fake Console までしか検証されていない**——`docs/workflow/metabolomics.md` を先に読む。

差次的エクスポートの列定義は別リポジトリ（massbank-context）との契約を兼ねる。加えて、
文献知識と再利用手順の蓄積層を持つ——`knowledge/` `playbook/` は **MCP リソース**として
索引と本文を配信し（`lipidmix://{knowledge,playbook}/{index,expand/<slug>}`）、
目的の `analyses/` とレポートの `reports/` はツール経由で読み書きする。

クライアントは Claude Desktop / Claude Code と、別リポの WebUI（下記「関連リポジトリ」）。

## コマンド

```bash
C:/Python314/python.exe -m pip install -r requirements-dev.txt   # 初回のみ。pytest はここにだけ入る
C:/Python314/python.exe -m pytest tests -q
```

- **リポジトリルートから `pytest` で実行する**。`-m unittest discover -s tests -t .` は
  これより少なく拾う（`test_dataset_state.py` `test_console_runner.py` `test_mztab_tools.py` 等、
  pytest 関数形式で書かれたファイルは `unittest.TestCase` を継承しないため拾えない）。
  **テストの数はここに書かない**。pytest の出力が正準（写しだけが腐るため）。
- MCP サーバ起動: `C:/Python314/python.exe server.py`（既定 stdio）。
- パーサ単体の CLI: `python -m metabolomix.arf.reader --file <path> --pca` など（README「Command-line examples」）。

**Python は `C:/Python314/python.exe` を使う**（PATH 上の `python` も同じ実体で、依存はここに入っている）。
Claude Code 用の `.mcp.json` はこれを絶対パスで指すが、`.git/info/exclude` でローカル除外されており
追跡されない。追跡されている `.vscode/mcp.json` は PATH 上の `python` を使う。`.venv/` は存在しない。
`.venv-1/` は 2026-07 のローカル LLM 検証用に残っている別環境で、通常の開発・テストでは使わない。

環境変数（すべて任意・上書き用）: `LIPIDMIX_DATA_DIR`（データ探索先。既定 `<project>/data`）/
`LIPIDMIX_KNOWLEDGE_DIR` `LIPIDMIX_PLAYBOOK_DIR` `LIPIDMIX_ANALYSES_DIR` `LIPIDMIX_REPORTS_DIR`（蓄積先の上書き用）/
`LIPIDMIX_TRANSPORT` `LIPIDMIX_HOST` `LIPIDMIX_PORT`（HTTP 待受）/ `LIPIDMIX_CAVEAT_MODE` /
`LIPIDMIX_PLOT_OUTPUT`（描画系ツールの戻り値。既定 `image`。Plotly で自分で描く
クライアント＝Use-LLLM は `payload` を置く）/
`LIPIDMIX_LIBRARY_CACHE_DIR`（照合用 SQLite キャッシュの置き場所。既定 `data/.library-cache`）/
`LIPIDMIX_CONFIG`（下記の設定ファイルの場所の上書き。worktree とテストが使う）。

**外部資産の場所は設定ファイルか環境変数で指す**（spec 2026-10-06）。リポジトリ直下の
`lipidmix.local.toml`（追跡外。雛形 `lipidmix.example.toml`）の `[msdial] exe` `lbm` と
`[library] msp_positive` `msp_negative`、または環境変数 `MSDIAL_EXE` `MSDIAL_LBM`
`MSDIAL_MSP_POS` `MSDIAL_MSP_NEG`（環境変数が優先）。読むのは `metabolomix/core/user_config.py`
の `get_setting` だけで、各リゾルバは `os.environ` を直に読まない。設定ファイルは呼ばれる
たびに読むので再起動は要らない。テストは conftest がこれらの環境変数を消し
`LIPIDMIX_CONFIG` を存在しないパスへ向けている（手元の設定ファイルを掴ませない）。

**研究室の参照ライブラリは外部流出禁止**。本体はリポジトリの外に置いて設定ファイルか環境変数で指し、
パスを追跡対象の文書・テストに書かない。本体とキャッシュの拡張子は `.gitignore` が
置き場所を問わず無視し、`tests/test_gitignore_library.py` がその網を縛っている。
照合の上位候補（化合物名・スコア）がツール戻り値に出るのは許容されている。

## 構成と、触るときの鉄則

ルート直下の `.py` は `server.py` と `check.py` の 2 つだけ。実装は `metabolomix/` にある。

**パッケージ名は `metabolomix`（2026-10-06 に `lipidmix` から改名）だが、外部契約は旧名のまま**:
環境変数 `LIPIDMIX_*`、MCP リソース URI `lipidmix://`、描画データの形式名
`lipidmix.<名前>.vN`（Use-LLLM が読む）、設定ファイル `lipidmix.local.toml`。
これらを `metabolomix` に揃えるのは将来の検討事項で、変えるなら移行期間と下流の同時更新が要る。
新しく足すものもこの規則に合わせる（新しい環境変数も `LIPIDMIX_` で始める）。

```
metabolomix/core/      FastMCP インスタンス・設定・セッション状態・パス解決・共通ヘルパ（依存グラフの leaf）
metabolomix/msdial/    MS-DIAL 固有サイドカー（*_tags.xml / .mddata）・同定・ピーク検証・サンプル因子
metabolomix/analysis/  入力形式に依存しない数値処理（前処理/QC・PCA・差次的解析・エクスポート契約）
metabolomix/plots/     描画 payload の組み立てと matplotlib 描画（volcano / eic / render）
metabolomix/{arf,arf2,pai2,dcl,eic}/   形式ごとの reader.py（パーサ）と tools.py（MCP ツール）
metabolomix/mztab/     mzTab-M リーダ・DatasetState 構築
metabolomix/library/   参照ライブラリ（.dbs / .msp）の読み取りと永続 store。
                    形式ごとに分けないのは、2 つが同じレコード形と同じ store を
                    共有する 2 つの入口にすぎないため（spec 2026-09-19）
metabolomix/console/   MS-DIAL Console 実行層（job_manager / runner / output_collector）
metabolomix/pipeline/  生データフォルダ起点の統括（受付 request/inputs/store・独立 worker が回す engine・
                    工程 handler の service・再開/取消の recovery・必須出力判定と品質レポートの report）
metabolomix/handoff/   Console 成果物の受け渡しスキーマ（analysis-job.json）
metabolomix/curation/  アラインメントのキュレーション（証拠収集・機械判別・フラグ記録・ビューア HTML）
metabolomix/corpus/    蓄積ノートの純ロジック（knowledge_store / paper_ingest）
metabolomix/tools/     形式に紐づかない MCP 公開層（入口・サンプル検索・目的・レポート・リソース）
```

- **`server.py` は薄いファサード**。import 副作用で `@mcp.tool` を登録し、`from ... import *`
  （各モジュールの `__all__` ＝そのモジュールのツール名）で公開面を再エクスポートするだけ。
  新しいツールを足すときは、実体を該当モジュールに書き `__all__` に載せる。
- **`server.py` をルートから動かさない**。`.mcp.json` / `.vscode/mcp.json` が絶対パスで指している。
- **可変状態の正準は `metabolomix.core.mcp_core`（`DATA_DIR` `KNOWLEDGE_DIR` `ANALYSES_DIR`）と
  `metabolomix.core.session_state`（`session`）**。`server.<name>` はスナップショット束縛にすぎないので、
  差し替え・モンキーパッチは必ず正準モジュール側に当てる。参照も `mcp_core.DATA_DIR` の
  module 修飾で行い、`from ... import DATA_DIR` を書かない。
- **`metabolomix.core.mcp_core` は leaf**。ここから `metabolomix.<形式>.tools` や `metabolomix.tools.*` を
  import してはならない（循環）。
- **セッション状態はパーサ別に分離済み**: `session.arf` / `.arf2` / `.pai2` / `.eic`。あるパーサが
  別スロットを触ってはいけない（`pai2_parser` が ARF の前処理行列を無言破棄した過去のバグの再発防止）。
- **解釈規則の選択軸は `session.assay_kind`**（`lipid` / `metabolite` / `unknown`、既定
  `unknown`）。脂質名文法は脂質アッセイでしか成り立たないので、ダイジェスト
  （`assay_digest()`）も `section_hint` もここで分岐する。**`caveat_emitted` の bool だけに
  戻さない**——`load_dataset` は GATEWAY より先に走るのでダイジェストは必ず `unknown` で
  1 回出る。`caveat_emitted_kind` があるから、種別が確定したときにもう 1 回だけ届く。
- **前提状態が無いときは例外でなく機械可読な封筒を返す**:
  `{"error": {"code": "missing_state", "state": ..., "required_tools": [...], "message": ...}}`
  （`metabolomix/core/mcp_errors.py` の `missing_state`）。クライアントはこれを読んでリプレイする契約。
- **reader の MessagePack Key インデックスの正解表は `docs/schema/*.md`**（MS-DIAL の C# クラス
  `[Key(N)]` から抽出した Key 番号表。上流のコミットと欠番も記載してある）。
  インデックス定数を変える前に必ず参照する。推測で直さない。
- **`run_pca` の正準は `metabolomix/analysis/pca.py`**。`metabolomix/arf/reader.py` の同名は後方互換の
  再エクスポートで、ARF テストが `patch.object(server.arf_reader, "run_pca", ...)` で module 属性
  としてこの束縛を差し替えてモックしている。**消すとモックが効かなくなり、テストは緑のまま
  実物の scikit-learn PCA が走り出す。**
- **差次的エクスポートの列定義は `metabolomix/analysis/export_contract.py` が唯一の正準**。
  ARF（`arf_export_differential`）と mzTab-M（`dataset_export_differential`）の両経路がここを
  共有しており、**別リポジトリ（massbank-context）との契約**でもある。列の追加・改名・並べ替えは
  `CONTRACT_VERSION` の引き上げと下流の同時更新なしにやってはいけない。
- **ツールの戻り値を肥大させない**。戻り値はそのまま LLM の文脈を占め、結論が埋没する。
  守るべき決まりごと:
  - JSON は `metabolomix.core.serialization.json_payload()` で返す（`json.dumps(..., indent=2)`
    を書かない。実測で戻り値の 15〜57% が空白だった）。
  - 行が並ぶ一覧は `arf2/reader.py format_spots_as_table()` の TSV（列名 1 回）。
  - float は丸めてから返す（既定 repr は 17 桁出る）。座標点列は `round_floats()`。
  - 全ツールに `structured_output=False` を付ける。付けないと FastMCP が outputSchema を
    導出し、MCP が同じ内容を content と structuredContent の**両方**で送る（＝2 倍）。
  - 座標配列など巨大な中間データは payload から外してセッションに保持する
    （図保存ツールがそこから読む）。
  - 図は座標を LLM に渡すより**サーバで描いて画像で返す**ほうが 2 桁安い
    （volcano 実測: 点列 183,578 字 ≒数万トークン → PNG 327 画像トークン）。
    画像トークンは `幅×高さ/750` なので dpi は上げない（`plots/render.py`）。
- **`metabolomix/pipeline/` は独立 worker プロセスが回す層**。守るべき決まりごと:
  - **`pipeline/{engine,worker,service,recovery,store}.py` はグローバル session を import しない**
    （`metabolomix.core.session_state` / `metabolomix.core.mcp_core` / `metabolomix.tools.*`。
    `tests/test_pipeline_engine.py` の AST テストが固定している）。進行状況は
    `pipeline-run.json` と、そのプロセスだけが持つ `runtime` dict に住む。
  - **`core/atomic_io.py` と `core/process_control.py` と `core/user_config.py` は stdlib だけの leaf**。
    別プロセスの worker が最初に import するので、重い依存を持ち込まない。
  - **run record の `results` は追記専用**（`store._assert_results_append_only` が
    既存要素の書換えを拒否する）。永続化する pipeline 内のパスは `pipeline_root` 相対。
  - **`pipeline_status` は読取専用**。生存確認 probe の失敗を「死んでいる」と読まない。
  - **Console の自動再試行は行わない**。やり直しは `rerun_upstream=true` の明示が要り、
    そのときは新しい attempt ディレクトリを作る（前回の終了証跡・ログを上書きしない）。
  - 工程の実行順の正準は `engine.build_stages` の計画順（`record["stages"]` の挿入順ではない）。
    report は必ず最後。
  - `store`（`_patch_request_id_with_retry`）と `recovery`（`prepare_resume`）は
    それぞれ自分の保存経路を `STATE_REVISION_CONFLICT` の有限回リトライで覆う。
    `engine` は長い区間を跨ぐ保存（`mark_stage_running` /
    `commit_stage_outcome`、`_save_with_retry`）だけをリトライで覆い、
    それ以外の `save_run` 直呼び出しは読み直しを挟まずコンフリクトを
    そのまま伝播する——ただしいずれも load してすぐ save するだけの
    短い区間なので、伝播した先は worker が死に、次の resume が状態を
    作り直すことで自己修復する。

## ドキュメントの地図（用途別に読み分ける）

| 知りたいこと | 見る場所 |
|---|---|
| ツールの引数・用途（`tests/test_readme_links.py` が一覧と件数を実登録と突き合わせている） | `USAGE.md` |
| 出力フィールドの**意味**（行の粒度・脂質名文法・必須注意） | `docs/output_format/core.md` ＋ トピック別（`arf` `arf2` `pai2` `dcl` `eic` `identity` `mztab` `library` `curation`）。MCP リソース `lipidmix://docs/output-format[/{topic}]` としても配信 |
| ツールが**どのファイルのどの関数をどの順に呼ぶか** | `docs/workflow/`（対象範囲の線引きと内訳は `index.md` が正準。`tests/test_workflow_docs.py` が実登録と突き合わせている）。行番号は書かない規約 |
| **生データ → Console → mzTab-M → 差次的解析 → パスウェイ**の一気通貫の順序と、内部関数の引数・戻り値 | [docs/superpowers/specs/2026-09-03-end-to-end-pipeline-design.md](docs/superpowers/specs/2026-09-03-end-to-end-pipeline-design.md)（**目標状態**の記述。実装状況は同文書 §9。完成後 `docs/workflow/Lipidmix/` へ昇格） |
| MessagePack の Key 番号 | `docs/schema/*.md` |
| パーサ単体の CLI（フラグ一覧と実行例） | `docs/cli.md` |
| 文献由来の設計候補（採択済み） | [docs/research-backlog.md](docs/research-backlog.md)（正準は vault 側。仕組みは spec 2026-09-14） |
| 設計判断の経緯・調査で判明した事実 | `docs/HISTRY.md`（綴りはこのまま。**追跡外＝ローカル専用ログ**） |
| 進行中/完了タスク | `docs/task.md`（**追跡外**。ステータス = TODO/DOING/DONE/HOLD） |
| 過去の設計書・計画書 | `docs/superpowers/{specs,plans,notes}/` |
| **利用者から見た**入口の分岐と、経路ごとの強制順序／オプション分岐の全体像 | `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`（**リポジトリ外**・Obsidian vault。維持義務は「作業の記録と Git」） |

## テストの規約

- 腐敗防止テストが効いている。壊すと落ちる:
  - `tests/test_workflow_docs.py` — `docs/workflow/` が挙げるパス・関数名を AST で実在検証し、
    対象ツール数を実登録数と突き合わせる。
  - `tests/test_package_layout.py` — 移動で静かに壊れる `BASE_DIR` 起点の解決と、
    ルート直下の `.py` が `server.py` `check.py` だけであることを縛る。
  - `tests/test_readme_links.py` — 入口文書（`README.md` `USAGE.md` `DEPLOY.md` `CLAUDE.md`）の
    相対リンクが実在するか、`USAGE.md` のツール集合・宣言件数が実登録と一致するか、
    CLAUDE.md 冒頭の規模表記が実登録と一致するかを検証する。
    加えて **README.md / CLAUDE.md に数量表現を書かせない**（正準は別文書にあり、
    写しだけが腐るため）。README は詳細を他文書へ委譲した要約なので、
    ポインタが切れると案内そのものが壊れる。
- `tests/test_server_registration.py` がツール/リソースの登録数と `ToolAnnotations` を検証する。
- **fixture はテスト自身が作る**。`analyses/` `knowledge/` の実ファイルに依存させない
  （追跡外なのでユーザ環境依存の不安定テストになる。tmp に作って `ANALYSES_DIR` を差し替える流儀）。

## 作業の記録と Git

- 調査・実装をしたら `docs/HISTRY.md` に追記し、`docs/task.md` のステータスを更新する。
- **この 2 つは追記専用**。日付見出しで区切って末尾に足し、既存の節は書き換えない。
  どちらも追跡外で git が競合を検出しないため、複数のエージェントが同時に走ると
  書き換えは後勝ちで静かに消える。
- **フローが変わったら vault 側の流れ図も同じ作業の中で直す**。反映先は
  `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`。
  対象は「利用者から見て何がどの順に起きるか」が変わる改修——入口の判定規則
  （`metabolomix/core/mcp_core.py` の `MCP_INSTRUCTIONS` の ENTRY POINT / GATEWAY）、
  pipeline の stage 列（`pipeline/engine.py` `build_stages` /
  `pipeline/stage_plan.py` `BASE_V2_STAGE_IDS`）、前提状態の連鎖（どのツールが
  どのツールを `missing_state` で要求するか）、経路の増減（新しい入口・新しい解析枝）、
  既定値の変更。**ツールの引数だけの変更は対象外**（正準は `USAGE.md`）。
  直したら末尾の出典行の日付も書き直す。
- この流れ図は**追跡外かつリポジトリ外**で、テストも git も腐敗を検出しない。
  後回しにすると、次のエージェントが古い流れ図を正しいものとして読む。
- `check.py` はスクラッチ（疑似ワークスペース）。一時検証コードをここに書き、通ったら適切な
  モジュールへ移して中身を消す。何もここに依存させない。
- **コミット時に全テストが自動で走る**（`.githooks/pre-commit`）。**所要は全テストの
  実行時間そのもので、数分かかる**。エージェントの既定コマンドタイムアウトを超えるので、
  自動化された `git commit` は最初からバックグラウンド実行にし、出力をファイルへ
  落として `FAILED` を後から引けるようにする。フックは index ではなく**作業ツリー**を
  検証するので、実行中にファイルを編集しない（その編集がそのテスト実行に混ざる）。
  クローン直後は `git config core.hooksPath .githooks` を 1 度実行して有効化する。
  迂回は `git commit --no-verify`（緊急時のみ）。
- `main` へのマージは `--no-ff`、`Merge <branch>: <日本語の要約>` 形式のマージコミット。
- マージ済みブランチは**ローカルのみ削除**し `origin` 側は残す。`push` は指示があっても都度確認する。

## 並列でエージェントを走らせるとき

作業ツリーを共有すると、片方の未コミット変更がもう片方のテスト結果に混ざり、HEAD も
断りなく動く。**エージェント 1 体につき worktree 1 本**に分ける。

```bash
sh scripts/new-worktree.sh <branch>    # .worktrees/ 配下に作る（/ は - に潰す）
git worktree remove .worktrees/<slug>  # 片付け
```

- **追跡外のものは worktree に来ない**: `data/` `analyses/` `.mcp.json`
  `docs/HISTRY.md` `docs/task.md`。テストはこれらに依存しない規約なので、
  新品の worktree でもテストは全数緑になる（実測で確認済み）。
- スクリプトが worktree 専用の `.mcp.json` を生成する。`server.py` は**その worktree の
  もの**を指す（main を指すと、worktree のコードを編集しながら main の実装を試すことに
  なり、最も気づきにくい形で嘘をつく）。蓄積状態（`LIPIDMIX_DATA_DIR`
  `LIPIDMIX_ANALYSES_DIR` `LIPIDMIX_KNOWLEDGE_DIR` `LIPIDMIX_REPORTS_DIR`）は
  **main ツリー**を向けて分裂させない。版管理対象の `playbook/` だけは worktree ローカル
  （変更対象そのものなので）。
- **記録は main ツリー側の `docs/HISTRY.md` / `docs/task.md` へ追記する**。
  worktree には存在しない。
- `core.hooksPath` は worktree 間で共有され `.githooks/` は追跡対象なので、
  pre-commit は worktree でもそのまま効く。
- **`git stash` を使わない**。stash スタックは main チェックアウトと全 worktree で
  共有されており、別のエージェントの退避を pop しうる。退避が要るなら WIP コミットにする。
- 統合は main ツリーで直列に行う。

## 踏みやすい罠

- `.gitignore` に `*.txt` があるため、**新規の .txt は無言で追跡漏れする**（`git add -f` が要る）。
- `docs/HISTRY.md` `docs/task.md` `data/` `analyses/` `archives/` は追跡外。`knowledge/` は
  種ノートのみ追跡、`playbook/` は版管理対象。クリーンチェックアウトに無い前提で書く。
- `archives/` は 2026-07 に退避したローカル LLM エージェント・解釈精度評価のコード置き場
  （`local_llm_agent/` `interp_eval/` ほか）。**ライブ側からの参照はゼロ**。動かすには
  ルート＋当該フォルダを `PYTHONPATH` に載せる必要がある参照用の退避物。現行開発では触らない。
- `.arf` は同一フォルダに `DriftSpots.arf` と `PeakProperties.arf` が併存しうる。リゾルバは
  `PeakProperties.arf` を自動選択し、複数バッチが混在する場合はファイル名の
  `AlignmentResult_<timestamp>` で最新バッチを選ぶ。

## 関連リポジトリ

`Use-LLLM`（本リポジトリと同じ親ディレクトリに置く別クローン） — このサーバをローカル LLM ＋ WebUI から操作する**別リポの独立システム**。
本体の変更では基本的に触らない。往復（`missing_state` 封筒の解釈やツールカタログ）を疑うときだけ、
そちらのメモリ（`use-lllm-*`）を参照する。
