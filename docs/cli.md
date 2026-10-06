# パーサ単体の CLI

MCP サーバを起動せずに、各 reader を直接叩いて `.arf` / `.arf2` / `.dcl` / `.EIC.aef` を
確認するための入口。MCP ツール経由の使い方は [../USAGE.md](../USAGE.md) を参照。

## 前提

- Python は `C:/Python314/python.exe`（`.mcp.json` / `.vscode/mcp.json` が指しているのと同じ環境）。
- 探索先はすべて `metabolomix.core.data_config.get_data_dir()` 経由。環境変数
  `LIPIDMIX_DATA_DIR` を設定するとその MS-DIAL 出力フォルダが対象になる（未設定時は `<project>/data`）。
- 出力（CSV / JSON / PNG）は、パスを渡さない限りカレントディレクトリに書かれる。
  生成物はリポジトリに追跡させない。

## `metabolomix.arf.reader` — `.arf`

唯一、引数を取る本格的な CLI。`--file` を省略するとデータディレクトリから自動選択する
（`--index` で候補の何番目かを指定）。

| フラグ | 用途 |
|---|---|
| `--file` / `-f` | 解析する `.arf` のパス |
| `--index` / `-i` | `--file` 省略時、自動検出した候補の何番目を使うか |
| `--export` / `-e` | ピークプロパティを CSV に書き出すパス |
| `--pca` | PCA を実行する |
| `--components` | 主成分数 |
| `--log-transform` | 標準化の前に log 変換する |
| `--min-detection-rate` | 検出率がこの値未満の特徴量を落とす |
| `--output-pca` | PCA 結果を JSON 保存するパス |
| `--output-plot` | PCA プロットの PNG 保存先（**未指定なら PNG を作らない**） |
| `--output-sample-scores` | サンプル別スコアプロットの保存先 |
| `--top-features` / `-t` | PCA Loading の上位・下位を何件表示するか |
| `--props` | 対象プロパティ（既定 `height`） |
| `--group-replicates` | ファイル名の複製番号（`_1`, `_2` …）を畳んで群にする |
| `--group-regex` | 畳む文字列を正規表現で直接指定する |
| `--plot-distribution` | 脂質クラス別の PeakHeight 分布プロットを作る |
| `--filter-name` | 分布プロットで抽出するキーワード（例 `TG`, `EtherPE_P`） |
| `--output-dist-plot` | 分布プロットの PNG 保存先（未指定なら作らない） |

ピークプロパティを CSV に書き出す:

```bash
C:/Python314/python.exe -m metabolomix.arf.reader --file "data/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf" --export output_peaks.csv
```

PCA を実行して結果と図を保存する:

```bash
C:/Python314/python.exe -m metabolomix.arf.reader --file "data/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf" --pca --output-pca pca_result.json --output-plot pca_plot.png --output-sample-scores sample_scores.png
```

PCA Loading の上位特徴量を見る:

```bash
C:/Python314/python.exe -m metabolomix.arf.reader --file "data/AlignmentResult_2026_05_15_10_13_35_PeakProperties.arf" --pca --top-features 10 --props height
```

## `metabolomix.arf2.reader` / `metabolomix.dcl.reader` / `metabolomix.eic.reader`

いずれも引数を取らず、データディレクトリから対象ファイルを自動選択して要約を標準出力に出す
簡易確認用。対象を切り替えるときは `LIPIDMIX_DATA_DIR` を変える。

```bash
C:/Python314/python.exe -m metabolomix.arf2.reader
```

## `metabolomix.library.store` — 参照ライブラリの store 事前構築

`library_load` と同じ解決規則でライブラリを選び、照合用の SQLite store を構築する
（既にあれば開くだけ）。大きな `.msp` の初回構築は MCP クライアントのタイムアウトに
当たりうるので、手元で 1 度走らせておくと以後の `library_load` はキャッシュを開くだけで
済む。標準出力にはファイル名・件数・所要秒数を 1 行の JSON で出す（置き場所は出さない）。

| フラグ | 意味 |
|---|---|
| `--ion-mode {positive,negative}` | どちらの極性の設定（環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`、無ければ `lipidmix.local.toml` の `[library] msp_positive` / `msp_negative`）を使うか |
| `--file <path>` | ライブラリのパス（環境変数より優先） |
| `--rebuild` | キャッシュと記憶した sha256 を無視して作り直す |

解決できないときは `MSP_AMBIGUOUS` などのコードを標準エラーに出して終了コード 2。

```bash
C:/Python314/python.exe -m metabolomix.library.store --ion-mode negative
```

## テスト

```bash
C:/Python314/python.exe -m pytest tests -q
```

リポジトリルートから `pytest` で実行する。`unittest discover` は pytest 関数形式で
書かれたテストを拾わないため件数が合わない。
