# パーサ単体の CLI

MCP サーバを起動せずに、各 reader を直接叩いて `.arf` / `.arf2` / `.dcl` / `.EIC.aef` を
確認するための入口。MCP ツール経由の使い方は [../USAGE.md](../USAGE.md) を参照。

## 前提

- Python は `.venv/Scripts/python.exe`（`.mcp.json` / `.vscode/mcp.json` が指しているのと同じ環境。setup が作る。mac / Linux は `.venv/bin/python`）。
- 読むファイルは必ず `--file` で渡す（必須）。**データフォルダの探索は無い**（B-2b で、環境変数
  `LIPIDMIX_DATA_DIR` と自動選択を消した）。
- 出力（CSV / JSON / PNG）は、パスを渡さない限りカレントディレクトリに書かれる。
  生成物はリポジトリに追跡させない。

## `metabolomix.arf.reader` — `.arf`

引数を取る本格的な CLI。`--file` は必須。PCA のフラグは無い（D 段で消した。PCA は MCP の `dataset_pca`）。

| フラグ | 用途 |
|---|---|
| `--file` / `-f` | 解析する `.arf` のパス（必須） |
| `--export` / `-e` | ピークプロパティを CSV に書き出すパス |
| `--group-replicates` | ファイル名の複製番号（`_1`, `_2` …）を畳んで群にする |
| `--group-regex` | 畳む文字列を正規表現で直接指定する |
| `--plot-distribution` | 脂質クラス別の PeakHeight 分布プロットを作る |
| `--filter-name` | 分布プロットで抽出するキーワード（例 `TG`, `EtherPE_P`） |
| `--output-dist-plot` | 分布プロットの PNG 保存先（未指定なら作らない） |

ピークプロパティを CSV に書き出す:

```bash
.venv/Scripts/python.exe -m metabolomix.arf.reader --file "data/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf" --export output_peaks.csv
```

## `metabolomix.arf2.reader` / `metabolomix.dcl.reader` / `metabolomix.eic.reader`

いずれも `--file`（必須）で渡したファイルの要約を標準出力に出す簡易確認用。

```bash
.venv/Scripts/python.exe -m metabolomix.arf2.reader --file "data/AlignmentResult_2026_05_15_10_13_35.arf2"
```

## `metabolomix.library.store` — 参照ライブラリの store 事前構築

`--file` で渡した（または `--ion-mode` の設定が指す）ライブラリから、照合用の SQLite store を構築する
（既にあれば開くだけ）。大きな `.msp` の初回構築は MCP クライアントのタイムアウトに
当たりうるので、手元で 1 度走らせておくと以後の `library_load` はキャッシュを開くだけで
済む。標準出力にはファイル名・件数・所要秒数を 1 行の JSON で出す（置き場所は出さない）。
キャッシュの置き場所は `library_load` と同じ（`lipidmix.local.toml` の `[paths] library_cache`、
既定は設定ファイルのフォルダの `data/.library-cache`）。

| フラグ | 意味 |
|---|---|
| `--ion-mode {positive,negative}` | どちらの極性の設定（`lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`）を使うか（`--file` が無いとき） |
| `--file <path>` | ライブラリのパス（設定より優先。どちらか一方が要る。データフォルダは探さない） |
| `--rebuild` | キャッシュと記憶した sha256 を無視して作り直す |

決まらないときは `LIBRARY_NOT_FOUND` などのコードを標準エラーに出して終了コード 2。

```bash
.venv/Scripts/python.exe -m metabolomix.library.store --ion-mode negative
```

## テスト

```bash
.venv/Scripts/python.exe -m pytest tests -q
```

リポジトリルートから `pytest` で実行する。`unittest discover` は pytest 関数形式で
書かれたテストを拾わないため件数が合わない。
