# 外部資産のパスを設定ファイルから読む（環境変数に頼らない導入） 設計

- 日付: 2026-10-06
- 発端: 外部の利用者から、クリーンインストールで動かない・テストが落ちるという issue
  （GitHub #1〜#3）が届いた。導入の障害を見直すなかで、MS-DIAL Console と参照ライブラリの
  場所を**環境変数でしか渡せない**ことが、クライアントごとの設定の重複と、外部の利用者の
  導入しづらさの両方の原因になっていると判断した。
- 範囲外の関連課題: `data_dir` の名残の撤去（`docs/task.md` 2026-10-06 節）、
  依存の宣言・Windows の Job Object・ミラー図（同節の A〜C。本 spec とは独立）。

## 1. 目的と完成条件

目的は 2 つ。

1. **クライアント間の重複をなくす**。Claude Code の `.mcp.json`、Claude Desktop の
   `claude_desktop_config.json`、VS Code の `.vscode/mcp.json`、Use-LLLM がそれぞれ同じパスを
   `env` ブロックに書いていて、値が食い違う（実際に食い違っている。`docs/task.md`
   2026-10-05 節）。サーバ自身が 1 か所から読めば、どのクライアントから起動しても同じになる。
2. **外部の利用者が導入しやすい**。clone した利用者が、Console やライブラリの場所を
   どこに書けばよいか迷わない。

方針: **設定ファイルを土台にし、未設定のときのエラーでファイルの場所と書き方を示す**。
対話式のセットアップ CLI や会話から設定するツールは後から足せる形にして、今回は作らない。

完成条件:

1. リポジトリ直下の `lipidmix.local.toml` に `[msdial] exe` を書けば、`env` ブロックの無い
   クライアント（`.vscode/mcp.json` など）からでも `console_plan` / `console_run` /
   `pipeline_run` が Console を見つける。
2. `[library] msp_positive` / `msp_negative` を書けば、`library_load(ion_mode=...)` が
   それを読む。
3. 環境変数を設定済みの既存環境（この機械の `.mcp.json`）は、何も変えずに今まで通り動く。
4. Console やライブラリが要る場面で何も設定されていなければ、エラーが「どのファイルの、
   どのキーに、何を書くか」を示す。
5. 設定ファイルを書き換えたら、MCP サーバを再起動しなくても次の呼び出しから効く。
6. テストは開発者の手元の `lipidmix.local.toml` を読まない。

## 2. 背景（2026-10-06 時点の実装）

| 設定 | いまの読み方 | 読んでいる箇所 |
|---|---|---|
| Console 本体 | 環境変数 `MSDIAL_EXE` | `console/runner.get_exe_path()`（`console/execution.py` `pipeline/service.py` `tools/console_tools.py` から呼ばれる）。**`tools/console_tools.py` の 2 か所は `get_exe_path` を迂回して `os.environ` を直に読む**（`console_method_template` の LBM 解決と、起動不能時の `details`） |
| 脂質ライブラリ | 環境変数 `MSDIAL_LBM` | `console/method_file.resolve_lbm(env=...)`。呼び出し側 3 か所（`tools/console_tools.py` 2 か所、`pipeline/inputs.py`）が `env=os.environ` を渡す |
| 研究室の参照ライブラリ | 環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG` | `core/path_resolvers.py` の `LIBRARY_ENV_VARS`・`_library_from_env`・`resolve_library_path` |

- `requirements.txt` に `python-dotenv` があるが、ライブコードは使っていない（`archives/` だけ）。
  `.env.example` は誰にも読まれないうえ、`LIPIDMIX_DATA_DIR` を最初の設定項目として勧めている。
- `.vscode/mcp.json`（追跡対象）には `env` が無く、VS Code から起動すると Console も
  研究室ライブラリも使えない。
- v2 の profile（`console/profiles.py`）は exe のパスと sha256 を profile ファイルに宣言して
  固定する。これは再現性のための仕組みで、本 spec の対象外。

## 3. 設定ファイルに載せるもの

**線引き**: 設定ファイルには「**リポジトリの外にあって、リポジトリに入れられない資産の
場所**」だけを載せる。クライアントや起動のしかたで変わる振る舞いは環境変数に残す。

| キー | 意味 | 対応する環境変数 |
|---|---|---|
| `[msdial] exe` | MS-DIAL Console の実行体（`MsdialConsoleApp.exe`） | `MSDIAL_EXE` |
| `[msdial] lbm` | 脂質ライブラリ `.lbm2`（任意） | `MSDIAL_LBM` |
| `[library] msp_positive` | 研究室の参照ライブラリ（正イオン） | `MSDIAL_MSP_POS` |
| `[library] msp_negative` | 研究室の参照ライブラリ（負イオン） | `MSDIAL_MSP_NEG` |

載せないもの（環境変数のまま。コードからも消さない）:

- `LIPIDMIX_DATA_DIR`: データフォルダは対話で渡すのが想定用途。この環境変数は CLI で
  ツールを関数として試していた頃の名残で、撤去は別課題。
- `LIPIDMIX_{KNOWLEDGE,PLAYBOOK,ANALYSES,REPORTS}_DIR`: 既定（リポジトリ直下）で外部の
  利用者は困らない。上書きが要るのは worktree 運用だけで、`scripts/new-worktree.sh` が
  環境変数で済ませている。
- `LIPIDMIX_LIBRARY_CACHE_DIR`: 設定する利用者がいない。
- `LIPIDMIX_TRANSPORT` / `HOST` / `PORT` / `PLOT_OUTPUT` / `CAVEAT_MODE`: クライアントや
  起動のしかたで変わる振る舞いそのもの（例: `PLOT_OUTPUT=payload` は Use-LLLM だけが置く）。
- pipeline の受付索引の置き場所（内部・テスト用）。

`lbm` は任意。書かなければ現行の推定（明示引数 → メソッドの宣言 → ビルド生成物 →
exe と同じフォルダ）がそのまま働く。設定ファイルの値は、現行で `MSDIAL_LBM` が占めている
位置に入る。

## 4. ファイルの形と規則

雛形 `lipidmix.example.toml`（追跡対象）:

```toml
# 複製して lipidmix.local.toml という名前でリポジトリ直下に置く（追跡されない）。
# Windows のパスは単一引用符で囲む（二重引用符だと \ がエスケープとして読まれる）。
# 相対パスは、このファイルのあるフォルダ基準で解く。

[msdial]
# MS-DIAL Console を使うときだけ必要。
exe = 'C:\MS-DIAL\MsdialConsoleApp.exe'
# lbm = 'C:\MS-DIAL\Lipids.lbm2'

[library]
# library_load(ion_mode=...) が読む参照ライブラリ。
# msp_positive = 'D:\lab-library\pos.msp'
# msp_negative = 'D:\lab-library\neg.msp'
```

規則:

- **探す順**: `LIPIDMIX_CONFIG` が指すファイル → `<repo>/lipidmix.local.toml` → どちらも
  無ければ「設定なし」（現行の既定のまま動く）。読むのは 1 ファイルだけで、複数をマージしない。
  `LIPIDMIX_CONFIG` が指すファイルが無いときは「設定なし」として扱い、未設定エラーの
  `details.config_file` にそのパスを載せる（worktree やテストが意図して空を指すため）。
- **値を引く順**: 環境変数 → 設定ファイル → 各リゾルバの既存の推定。環境変数が空文字や
  空白だけなら未設定とみなす（現行の `.strip()` と同じ）。
- 設定ファイルの値の**相対パス**は設定ファイルのあるフォルダ基準で解き、絶対パスにして
  返す。`~` は展開する。環境変数の値は加工しない（§5 `Setting`）。
- 設定ファイルは UTF-8 で読み、先頭の BOM は許す（メモ帳が付けることがあり、`tomllib` は
  BOM を構文エラーにする）。
- **環境変数が設定されているキーについては設定ファイルを読まない**。設定ファイルが壊れて
  いても、環境変数で渡している値は今まで通り使える。
- **`<repo>`** は `core/data_config.py` と同じく `Path(__file__).resolve().parents[2]`。
  worktree ではその worktree のルートになる（worktree へ main の設定を届けるのは §8）。
- `lipidmix.local.toml` は `.gitignore` に入れる（研究室ライブラリのパスが入るため）。

## 5. モジュール構成

新設は `lipidmix/core/user_config.py` の 1 つだけ。

- **stdlib だけの leaf**（`tomllib` `os` `pathlib` `dataclasses`）。pipeline の worker も
  import するので、`core/atomic_io.py` / `core/process_control.py` と同じ扱いにする。
  `mcp_core` / `session_state` / `lipidmix.tools.*` を import しない。
- 公開面:
  - `SETTINGS`: キー → 環境変数名の対応表（`"msdial.exe": "MSDIAL_EXE"` など 4 件）。
    環境変数名の正準はここに置き、`path_resolvers.LIBRARY_ENV_VARS` はここから導く。
  - `config_file_path() -> Path`: §4 の探す順で決めた設定ファイルのパス（存在しなくても返す）。
  - `get_setting(key) -> Setting | None`: §4 の引く順で値を返す。未設定なら `None`。
  - `Setting`: `key`、`value`（`str`。環境変数の値は前後の空白を除いてそのまま——
    既存の利用者とテストが相対名や `/` 区切りを置いているので形を変えない。設定ファイルの
    値は §4 の規則で絶対パスにする）、`source`（`"env"` / `"config_file"`）、`env_var`、
    `config_file`。
  - `setting_label(setting) -> str`: エラー文で値の出どころを言う句
    （「環境変数 MSDIAL_EXE」/「lipidmix.local.toml の [msdial] exe」）。
  - `missing_hint(key, what) -> str`: 未設定のとき利用者に伝える文（§7.1）。
  - `describe_missing(key, setting=None) -> dict`: 未設定・指す先が無いエラーの
    `details` に足す内容（§7.1）。値そのものは載せない。
  - `ConfigInvalidError(Exception)`: `code="CONFIG_INVALID"`、`message`、`line` / `column`
    （分かれば）。設定ファイルが読めないとき `get_setting` が送出する。
- **読むのは呼ばれるたび**。設定ファイルのパスと `st_mtime_ns` をキーに、解析結果を
  モジュール内にキャッシュする。書き換えれば次の呼び出しから効く（再起動不要）。
- 環境変数の優先順は `get_setting` の中だけに書く。各リゾルバは自分で `os.environ` を
  読まなくなる。

## 6. つなぎ方

| 箇所 | 変更後 |
|---|---|
| `console/runner.get_exe_path()` | `user_config.get_setting("msdial.exe")` を読む。未設定なら `MsdialExeNotFoundError`（メッセージは §7.1）。戻り値は `str` のまま |
| `tools/console_tools.py` の `os.environ["MSDIAL_EXE"]` 直読み 2 か所 | `get_exe_path` に寄せる。`console_method_template` は exe が無くても動く現行の振る舞い（`exe=None`）を保つため、未設定の例外を捕まえて `None` にする |
| `method_file.resolve_lbm(..., env=...)` | 引数 `env: dict` を `lbm_setting: user_config.Setting \| None` に置き換える（`MSDIAL_LBM` の位置に入る値と、その出どころ）。呼び出し側 3 か所は `user_config.get_setting("msdial.lbm")` の戻り値をそのまま渡す。`LBM_NOT_FOUND` の文面は `Setting.source` で言い分ける。`user_config` は leaf なので `method_file` から import してよい |
| `path_resolvers._library_from_env` と、両極性の判定（`configured`） | `user_config.get_setting("library.msp_positive" / "library.msp_negative")` を読む。関数名は `_library_from_setting` に改める。`resolve_library_path` の解決順（明示 → `ion_mode` の極性 → `*_Loaded.msp2.dbs` → 極性の設定 → `*.msp`）は変えない |

`ConfigInvalidError` は、各ツールの既存のエラー封筒（`console_error` / `LibraryPathError`
を包む既存の経路）で `CONFIG_INVALID` として返す。pipeline の受付（`pipeline/service.py`
`_resolve_exe_path`）では `DomainError("CONFIG_INVALID", ...)` に包む。

reanalysis-study plan Task 9（Console に渡す MSP を、明示 → メソッド → 極性の設定の順で
決める）は、本 spec の後に `user_config.get_setting("library.msp_<極性>")` を引くだけで
載る。`msp_positive` / `msp_negative` は `library_load` と Console の両方の正になる。

## 7. エラーの扱い

**原則**: 設定の不備で**サーバの起動は止めない**。止まるのはその設定を必要とするツール
だけ。設定は呼ばれるたびに読むので、壊れていてもパーサや解析のツールには影響しない。

### 7.1 未設定・指す先が無い

- エラーコードは**変えない**（`MSDIAL_EXE_NOT_FOUND` / `LBM_NOT_FOUND` /
  `MSP_ENV_NOT_FOUND`）。Use-LLLM がコードで分岐しているため、コードは契約として据え置く。
  `MSP_ENV_NOT_FOUND` は値が設定ファイルから来た場合に名前と実態が少しずれるが、
  互換性を優先する。
- `details` に足す（`user_config.describe_missing`）:

  ```json
  {"setting": "msdial.exe", "env_var": "MSDIAL_EXE", "source": null,
   "config_file": "C:\\...\\lipidmix.local.toml", "config_file_exists": false,
   "example": "lipidmix.example.toml"}
  ```

  `source` は値があったときの出どころ（指す先が無い場合）、未設定なら `null`。
- メッセージは利用者にそのまま伝えられる文にする。
  - 未設定: 「`lipidmix.example.toml` を `lipidmix.local.toml` として複製し、`[msdial] exe`
    に MsdialConsoleApp.exe のパスを書いてください（環境変数 `MSDIAL_EXE` でも指定できます）」。
    設定ファイルが既にあれば「`lipidmix.local.toml` の `[msdial] exe` に…」とする。
  - 指す先が無い: 出どころに合わせて「環境変数 `MSDIAL_EXE` が指すファイルがありません」か
    「`lipidmix.local.toml` の `[msdial] exe` が指すファイルがありません」。
- 研究室ライブラリは**置き場所を出さない規則を守る**。`msp_*` の値については、メッセージにも
  `details` にもファイル名とキー名だけを載せ、ディレクトリを出さない。設定ファイル自身の
  パスは出してよい。

### 7.2 設定ファイルが読めない

TOML の構文エラー、または値の型違い（文字列でない）なら `CONFIG_INVALID`。メッセージに
設定ファイルの名前と行・列を載せる（`details.config_file` にパス）。構文エラーで
ファイルに `\` が含まれるときは、「Windows のパスは単一引用符で囲んでください」と添える
（`"C:\MS-DIAL"` は `Unescaped '\'`、`"C:\Users"` は `Invalid hex value` になる。
2026-10-06 に Python 3.14 の `tomllib` で確認）。

構文としては通ってしまう誤りも拾う: `"C:\new\tool.exe"` の `\n` `\t` はエスケープとして
読まれ、値に制御文字が入る。値に改行・タブ・`\b` `\f` `\r` が含まれていたら
`CONFIG_INVALID` にして同じ案内を添える。

### 7.3 打ち間違い

未知の節・キー（`[msdail]`、`msp_postive` など）は、§7.1 の未設定エラーの中で
「認識できないキーがあります: `[msdail]`」と添え、`details.unknown_keys` に列挙する。値が
見つからない原因はたいてい打ち間違いなので、エラーが出たその場で気づけるようにする。
読むたびの警告は出さない（戻り値を肥大させない）。

### 7.4 `required_tools`

付けない。直すのは利用者がファイルを編集することで、呼べば解決する MCP ツールは今は無い。
会話から設定するツールを足したときにここへ載せる。

## 8. 移行

- この機械の `.mcp.json` の `env` ブロックはそのままでも動く（環境変数が優先）。設定
  ファイルへ移すかは実装後にユーザーが判断し、`docs/task.md` 2026-10-05 節の
  「`.mcp.json`・User 環境変数・`claude_desktop_config.json` の値の食い違い」の整理と
  一緒に片づける。
- `scripts/new-worktree.sh`: main ツリーに `lipidmix.local.toml` があれば、生成する
  `.mcp.json` の `env` に `LIPIDMIX_CONFIG`（main の設定ファイル）を入れる。
- `.vscode/mcp.json` と Use-LLLM: 変更不要。設定ファイルがあれば自動で効く。
- `.env.example` は削除し、`lipidmix.example.toml` に置き換える。

## 9. テスト

- **隔離**（`tests/conftest.py` の autouse fixture を広げる）: `MSDIAL_MSP_*` に加えて
  `MSDIAL_EXE` / `MSDIAL_LBM` を消し、`LIPIDMIX_CONFIG` を存在しない tmp のパスに向ける。
  開発者の手元の `lipidmix.local.toml`（研究室ライブラリのパス入り）をテストが掴まない
  ようにする。設定ファイルを使うテストは tmp に toml を書いて `LIPIDMIX_CONFIG` で指す。
  subprocess の worker は環境変数を受け継ぐので同じく隔離される。
- **新規** `tests/test_user_config.py`:
  - 探す順（`LIPIDMIX_CONFIG` → `<repo>` → 設定なし）、`LIPIDMIX_CONFIG` が無いファイルを
    指すときは設定なし
  - 引く順（環境変数が設定ファイルに勝つ、空文字の環境変数は未設定）
  - 相対パスを設定ファイルのフォルダ基準で解く、`~` の展開
  - 書き換えると次の呼び出しから効く（`st_mtime_ns` による再読込）
  - `CONFIG_INVALID`（構文エラー・型違い・二重引用符の `\` の案内、行と列）
  - 未知の節・キーが `describe_missing` に載る
  - `core/user_config.py` が stdlib 以外を import しない（AST で縛る）
- **既存テストの調整**: `resolve_lbm(env=...)` の呼び出し（`tests/test_console_runner.py`
  ほか）を `lbm_setting` に直す。各リゾルバに「設定ファイルから値が来る」「出どころで
  文面が変わる」ケースを足す（exe・lbm・msp_* それぞれ）。
- **網**: `tests/test_gitignore_library.py` に、`lipidmix.local.toml` が `.gitignore` で
  無視されることを加える。
- **pipeline**: `pipeline_run` の受付が設定ファイルの exe で通ること（fake Console）。

## 10. 文書

- `lipidmix.example.toml`（新規・追跡対象）: §4 の雛形。
- `README.md`: 導入手順に「Console や研究室ライブラリを使うなら、雛形を複製してパスを書く」
  を 1 段落足す。数量表現を書かない規約に従う。
- `USAGE.md`、`docs/output_format/library.md`、`docs/workflow/`（Console・library の節）:
  環境変数の説明を「環境変数か設定ファイル」に改め、`CONFIG_INVALID` と `details` の新しい
  項目を載せる。
- `CLAUDE.md`: 環境変数の段落を書き換え、設定ファイルの節を足し、`core/user_config.py` を
  stdlib だけの leaf の一覧に加える。
- vault の流れ図（`ms-data-parser-flow.md`）: Console・ライブラリ入口の前提条件の書き方が
  変わるので更新し、出典行の日付を直す。

## 11. 範囲外

- 会話から設定を書き込む MCP ツール（`config_set` など）。
- 対話式のセットアップ CLI と、Console の自動探索（MS-DIAL は zip を任意の場所に展開して
  使うことが多く、当たり外れが大きい）。
- `data_dir` の名残の撤去（`docs/task.md` 2026-10-06 節）。
- v2 profile の exe 宣言の扱い（profile の宣言を正とし、設定ファイルは関与しない）。
- `python-dotenv` の依存からの削除（依存の宣言を直す issue A で扱う）。
