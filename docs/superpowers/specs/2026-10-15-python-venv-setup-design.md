# Python の仮想環境化と 1 コマンドセットアップ（設計）

日付: 2026-10-15 / ブランチ: `feat/python-venv`

## 1. 目的と合意事項

**目的**: 配布先ユーザー（研究室の他メンバーに限らない）が、クローン後に迷わず
MCP サーバを動かせるようにする。現状は「Python を各自で入れ、グローバルへ
`pip install -r requirements.txt`、MCP クライアント設定を手で書く」で、
版の固定も仮想環境も無い（ユーザーから「セットアップがわかりにくい」との不満。task.md 2026-10-14）。

ユーザー判断（2026-10-15）:

- 対象は**配布先ユーザー**。
- 道具は **uv 優先、無ければ標準 venv + pip** へ落ちる。
- 1 コマンドに含めるもの: **Python `.venv` の構築 / MCP クライアント設定の生成 /
  `lipidmix.local.toml` 雛形の配置**。R（`Rscript r/setup.R`）は含めず別手順のまま。
- 入口は **`setup.ps1`（Windows）と `setup.sh`（mac/Linux・Git Bash）の 2 つ**。
  共通の後半は Python スクリプトに寄せて二重実装しない。
- 版は **`uv.lock` を正準**とし、pip 経路用に書き出したロックを追跡する。

前提（設計者の仮定、ユーザー承認済み）: 開発側（pre-commit・worktree・`.mcp.json`・
CLAUDE.md）も同じ `.venv` へ移す。`.venv-1/` は削除する。

## 2. ファイル構成

| ファイル | 役割 |
|---|---|
| `pyproject.toml`（新規） | `[project]` に直接依存（現 `requirements.txt` の内容。`mcp<2` を保持）、`requires-python = ">=3.13"`。`[dependency-groups] dev = ["pytest"]`。**`[tool.uv] package = false`**。 |
| `.python-version`（新規） | `3.14`。uv だけ持つ利用者には uv がこの版を取得する。 |
| `uv.lock`（新規・追跡） | 全依存の版の正準。 |
| `requirements.lock.txt` / `requirements-dev.lock.txt`（新規・追跡） | `uv export --no-hashes --no-emit-project`（dev は `--group dev` 付き、本体は `--no-dev`）で `uv.lock` から生成。venv+pip 経路が同じ版を入れるため。 |
| `requirements.txt` / `requirements-dev.txt` | それぞれ `-r requirements.lock.txt` / `-r requirements-dev.lock.txt` の 1 行に縮める。既存利用者の**旧版** `server_update` がこのファイル名を読むので移行期間は残す。 |
| `.gitignore` | `*.txt` 無視に対し `!requirements*.txt` の例外を足す。 |
| `setup.ps1` / `setup.sh`（新規・ルート） | 下記 §3。 |
| `scripts/setup_finish.py`（新規） | 下記 §3 の後半。 |

**`package = false` の理由**: プロジェクトを editable install すると `.pth` が main ツリーの
`metabolomix` を指し、worktree から起動した `server.py` / pytest が main のコードを
拾いうる。インストールしなければ import は常に `server.py` のあるツリー（＝cwd）から
解決されるので、全 worktree が main ツリーの `.venv` を安全に共有できる。

## 3. セットアップの流れ

`setup.ps1 [-Dev] [-NoUv]` / `sh setup.sh [--dev] [--no-uv]`（リポジトリ直下で実行）:

1. **uv があり `-NoUv` でない** → `uv sync`（`-Dev` なら `--group dev`、無ければ `--no-dev`）。
   Python が無ければ uv が `.python-version` の版を取得する。
2. **uv が無い** → Python 3.13 以上を探す（Windows: `py -3.14` → `py -3` → `python`、
   sh: `python3` → `python`。版を `sys.version_info` で確認）。見つからなければ
   「uv か Python 3.13+ を入れてください」と導入先 URL を示して終了。
   見つかれば `python -m venv .venv` → `.venv の python -m pip install -r requirements(-dev).lock.txt`。
3. 共通: `.venv の python scripts/setup_finish.py` を実行。

`-NoUv` は uv があっても 2 へ進ませる検証用フラグ（利用者向けにも残してよい）。
どの段で失敗しても非 0 で終わり、何が失敗したかを 1 行で述べる。

`scripts/setup_finish.py`（`.venv` の Python で走る。標準ライブラリ＋依存の import のみ）:

1. **import の自己確認**: `pyproject.toml` の直接依存を import し、失敗したものを列挙して非 0 終了。
   （パッケージ名と import 名の対応表をスクリプトに持つ: `PyMuPDF→fitz`、`scikit-learn→sklearn`、`pygoslin→pygoslin` 等。）
2. **`lipidmix.local.toml`**: 無ければ `lipidmix.example.toml` を複写し、場所と書く項目を案内。あれば触らない。
3. **`.mcp.json`**: 無ければ `ms-data-parser` を `.venv` の Python の絶対パス＋`server.py` の
   絶対パスで書く（`env` に `PYTHONIOENCODING=utf-8`）。**あれば上書きしない**。
   既存の `ms-data-parser.command` が生成値と違えば、その差を表示して手で直すよう案内する。
4. **Claude Desktop 用の貼り付け断片**（`claude_desktop_config.json` の `mcpServers` 項）を表示する。
   Desktop の設定ファイル自体は書き換えない。同名サーバの重複に注意する旨も添える（task.md の既知問題）。
5. 次の手順を表示: パスウェイを使うなら `Rscript r/setup.R`、更新は `server_update`。

## 4. 開発側と `server_update` の切替

**解釈系の選択**（どこも同じ規則: `.venv` があれば使い、無ければ従来どおり）:

- `.githooks/pre-commit`: 自ツリーの `.venv` → main ツリー（`git rev-parse --git-common-dir`
  の親）の `.venv` → `C:/Python314/python.exe` → PATH の `python`。
  `Scripts/python.exe` と `bin/python` の両方を見る。
- `scripts/new-worktree.sh`: 生成する `.mcp.json` の `command` を **main ツリーの `.venv` の Python**
  にする（無ければ従来の `C:/Python314/python.exe`）。
- `.vscode/mcp.json`: `${workspaceFolder}\\.venv\\Scripts\\python.exe`。
- 開発者自身の `.mcp.json`（追跡外）: 既存の env（MSP 等）を保つため、setup_finish の差分表示に
  従い `command` だけを手で直す。

**`server_update`（`metabolomix/core/version.py`）**:

- 依存の変化の判定対象を `uv.lock` / `requirements.lock.txt` にする（`requirements.txt` は
  ポインタになり、それ単独では変化を表さない）。
- 入れ直し: 走行中の Python で `pip` が import できれば `sys.executable -m pip install -r requirements.lock.txt`。
  できなければ（uv 製の venv は pip を持たない）`uv pip install --python <sys.executable> -r requirements.lock.txt`。
  どちらも使えなければ `updated_with_warning` で手順を案内する。
  `uv sync` は使わない（dev group の pytest など余分なものを消すため）。
- 失敗時の案内文と、`C:/Python314` 固定・`.venv` 無しと書いた古いコメントを直す。

**pipeline / Console worker**: 既に `sys.executable` で起動し、`core/process_control.py` の
`resolve_python_launch` が venv のリダイレクタを処理しているので**変更しない**（§6 の手動確認で裏を取る）。

**整理**:

- `.venv-1/` を削除（追跡外。参照は `archives/` と `.superpowers/sdd/` の過去記録だけ。過去記録は書き換えない）。
- 文書: CLAUDE.md「コマンド」（`.venv` の Python を使う、`.venv/` は存在しない の記述を削除、
  `pip install -r requirements-dev.txt` を `setup.ps1 -Dev` に）、README Quick start、
  DEPLOY Member Setup / Receiving Updates（既存利用者の移行: 1 度 `setup.ps1` を実行し、
  クライアント設定の `command` を `.venv` に変える）、USAGE §16 の該当箇所。
- vault の流れ図（`ms-data-parser-flow.md`）のセットアップ入口を更新し、出典行の日付を直す。

## 5. 自動テスト

- **ロックの整合**（`tests/test_dependency_lock.py`）: `uv.lock` を `tomllib` で読み、
  `requirements.lock.txt` / `requirements-dev.lock.txt` の `name==version` 集合が `uv.lock` の
  該当パッケージの版と一致すること。uv 不要で走る。書き出し忘れを捕まえる。
- **ポインタと設定**: `requirements.txt` が `-r requirements.lock.txt` だけであること、
  `pyproject.toml` に `[tool.uv] package = false` が残っていること。
- **`setup_finish.py`**: tmp で、無い→生成（`command` が `.venv` の Python の絶対パス）、
  ある→バイト不変かつ差分を表示。import 自己確認が欠けた依存を列挙して非 0 になること（import を差し替えて）。
- **`server_update` の依存段**: lock の変化で入れ直しが走る／走らない、pip の有無で
  `-m pip` と `uv pip` が選ばれる、両方無いと `updated_with_warning`。既存のモック流儀で、実際には入れない。
- **gitignore**: `tests/test_gitignore_library.py` と同じ流儀で `requirements*.lock.txt` が無視されないこと。

## 6. 手動検証（新品の worktree で）

1. **uv 経路**: `setup.ps1 -Dev` → `.venv` で全 pytest が緑（pipeline / Console worker の
   テストが venv 下で通る＝リダイレクタ処理の裏取り）。
2. **venv フォールバック**: 別の新品 worktree で `setup.ps1 -NoUv -Dev` → 標準 venv+pip で構築され
   import 自己確認が通る。
3. **`setup.sh`**: Git Bash で uv 経路を 1 回。mac/Linux 実機は無いので**未検証と記録**する。
4. **MCP 起動**: 生成された `.mcp.json` で Claude Code から接続し、`list_data_files` が応答する。

## 7. 範囲外

- R 環境の構築を 1 コマンドに含めること（`Rscript r/setup.R` のまま）。
- Claude Desktop の設定ファイルの自動書き換え。
- 外部契約名（`LIPIDMIX_*` 等）の変更、Python パッケージとしての配布（PyPI）。
- `archives/` を新しい環境で動かすこと。
