# 再解析の道具（ハンドル・spec・証拠の束）の出力定義

**再解析で使う**ツールの戻り値と、spec・結果の JSON の欄の意味（spec 2026-10-22）。
対話の使い方の出力は変わらない——既存のツールに足した引数（`dataset_handle` ほか）は、
渡さない限り戻り値を 1 文字も変えない。呼び出し連鎖は `docs/workflow/reanalysis.md`。

共通の約束:

- **spec は正規化して返す。** 既定値を埋め、未知のキーは拒否する（黙って無視しない）。
  結果には正規化した spec と `rule_hash` が入る。`rule_hash` は正規化した spec から
  パスの欄（`path` / `matrix_path`）を除いて作る——置き場所が違っても同じ規則なら同じ hash。
  入力のファイル（行列・群・共変量・注入順の TSV）は `input_fingerprint` に sha256 で入る。
- **同じ spec・同じ入力なら同じ JSON。** 時刻や乱数の ID は結果に入れない。
- **判定（再現したか・明らかに別物か）は返さない。** 道具が返すのは数字と、spec に書いた
  規則を機械的に当てた結果だけ。
- **研究室の参照ライブラリのパスは出さない。** ファイル名と sha256 だけ。

## 1. データセットのハンドル

ハンドル（`dh_` + 12 桁）は、1 つのサーバのプロセスにデータセット（とその極性の
ライブラリ）を複数並べて持つための ID。`dataset_handle_open` が返し、`dataset_*`・
`library_*`・再解析のツールの `dataset_handle` に渡す。

- **対話のセッションとは別の状態。** ハンドルを渡した呼び出しは、ハンドルを渡さない呼び出しが
  使う状態（`dataset_load` で読んだデータセット・`library_load` のライブラリ・直近の照合）を
  読みも書きもしない。閉じた・無いハンドルは `DATASET_HANDLE_NOT_FOUND`（開いている
  ハンドルの一覧つき）で、対話のセッションへ黙って落ちない。
- **ハンドルどうしも別の状態。** 同じ入力を開き直すと読み込み済みのもの（`reused_base`）を
  使い回すが、ハンドルには複製が渡る。片方の前処理はもう片方に見えない。入力のファイルの
  size か更新時刻が変われば読み直す。
- **並行と直列。** ハンドルを渡した呼び出しはサーバの別スレッドで回り、別のハンドルどうしは
  並行に進む。同じハンドルへの呼び出しは 1 本ずつ（到着順は保証しない）。
- **上限。** `LIPIDMIX_MAX_DATASET_HANDLES`（既定 8）。超えると `DATASET_HANDLE_LIMIT`。
  寿命は明示の `dataset_handle_close` か、サーバの終了まで。

`dataset_handle_open` の戻り値:

| 欄 | 意味 |
|---|---|
| `dataset_handle` | ハンドルの ID |
| `polarity` / `polarity_source` | データセットの極性と出所（`job_declared` / `adduct_majority` / `unknown`）。`dataset_feature_search` の目標照合と同じ決め方 |
| `n_features` / `n_samples` | mzTab-M の特徴数・試料数 |
| `detection_available` | 検出状態（gap-fill の区別）を取り込めたか |
| `library` | `{file, sha256, record_count}`。読まなかったときは `null` で、理由が `library_note` |
| `arf2` / `dcl` | スポットの照合（spot_id = SMF_ID）に使うファイル名。決まらなければ `null` |
| `reused_base` | 読み込み済みのものを使い回したか |

ライブラリは**データセットの極性の設定だけ**から引く（`[library] msp_positive` /
`msp_negative`、環境変数 `MSDIAL_MSP_POS` / `MSDIAL_MSP_NEG`）。極性が決まらない・設定が
無いときはライブラリなしで開き、逆の極性のライブラリを掴まない。読み取り専用で開く。
