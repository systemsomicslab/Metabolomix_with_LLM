# パーサ単体の CLI

MCP サーバを起動せずに reader を直接呼んで `.arf` / `.arf2` / `.dcl` / `.EIC.aef` を確かめる入口。MCP ツールは
[../USAGE.md](../USAGE.md)。

## 前提

- Python は `.venv/Scripts/python.exe`（`.mcp.json` / `.vscode/mcp.json` と同じ環境。setup が作る。mac / Linux は `.venv/bin/python`）。
- 読むファイルは `--file` で渡す（必須）。データフォルダの探索も環境変数 `LIPIDMIX_DATA_DIR` も無い。
- 出力（CSV / JSON / PNG）は、パスを渡さなければカレントディレクトリに書く。生成物はリポジトリで追跡しない。

## `metabolomix.arf.reader` — `.arf`

`--file` は必須。PCA のフラグは無い（PCA は MCP の `dataset_pca`）。

| フラグ | 用途 |
|---|---|
| `--file` / `-f` | 解析する `.arf` のパス（必須） |
| `--export` / `-e` | ピークプロパティを CSV に書き出すパス |
| `--group-replicates` | ファイル名の複製番号（`_1`, `_2` …）を畳んで群にする |
| `--group-regex` | 畳む文字列を正規表現で直接指定する |
| `--plot-distribution` | 脂質クラス別の PeakHeight 分布プロットを作る |
| `--filter-name` | 分布プロットで抽出するキーワード（例 `TG`, `EtherPE_P`） |
| `--output-dist-plot` | 分布プロットの PNG 保存先（未指定なら作らない） |

```bash
.venv/Scripts/python.exe -m metabolomix.arf.reader --file "data/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf" --export output_peaks.csv
```

## `metabolomix.arf2.reader` / `metabolomix.dcl.reader` / `metabolomix.eic.reader`

`--file`（必須）のファイルの要約を標準出力に出す。

```bash
.venv/Scripts/python.exe -m metabolomix.arf2.reader --file "data/AlignmentResult_2026_05_15_10_13_35.arf2"
```

## `metabolomix.library.store` — 参照ライブラリの store 事前構築

`--file`（または `--ion-mode` の設定）のライブラリから照合用の SQLite store を作る（既にあれば開くだけ）。大きな `.msp` の初回構築が
MCP クライアントの時間切れに当たるときに手元で 1 度走らせると、以後の `library_load` はキャッシュを開くだけになる。標準出力に
ファイル名・件数・所要秒数を 1 行の JSON で出す（置き場所は出さない）。キャッシュの置き場所は `library_load` と同じ
（`lipidmix.local.toml` の `[paths] library_cache`、既定は設定ファイルのフォルダの `data/.library-cache`）。

| フラグ | 意味 |
|---|---|
| `--ion-mode {positive,negative}` | どちらの極性の設定（`lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`）を使うか（`--file` が無いとき） |
| `--file <path>` | ライブラリのパス（設定より優先。どちらか一方が要る。データフォルダは探さない） |
| `--rebuild` | キャッシュと記憶した sha256 を無視して作り直す |

決まらないときは `LIBRARY_NOT_FOUND` などのコードを標準エラーに出し、終了コード 2。

```bash
.venv/Scripts/python.exe -m metabolomix.library.store --ion-mode negative
```

## テスト

```bash
.venv/Scripts/python.exe -m pytest tests -q
```

リポジトリルートから `pytest` で走らせる（`unittest discover` は pytest 関数形式のテストを拾わない）。
