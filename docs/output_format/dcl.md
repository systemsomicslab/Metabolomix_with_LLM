# DCL (`.dcl` / デコンボリューション済み MS/MS)

`.dcl` パーサ、`dcl_parser` / `dcl_find_msms`、および PAI2 ピークへの MS/MS 付与の出力定義。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

> **MS/MS の所在（重要）**: 実スペクトルを持つのは `.dcl` であって `.pai2` ではない。PAI2 の `has_msms` は**取得参照の有無**を示すだけで、`msms_spectrum` が非空である保証はない（共通核 §9-8）。同定確度（MSI Level 2）は MS/MS を根拠にするため、フラグだけで確度を主張してはならない。`pai2_parser` は同名 `.dcl` を自動で読み、索引対応で PAI2 ピークへ MS/MS を付与する。

## 6. DCL (`lipidmix/dcl/reader.py`)

### 6.1 `deserialize_dcl()`

型は `list[dict]`。**1要素は1 MSDecResult、すなわち1プリカーサーピークに対応するデコンボリューション済みMS/MS結果**である。件数は同名PAI2のピーク数と一致する設計である。

| キー | 型 | 意味 |
|---|---|---|
| `dcl_index` | int | DCL内の0始まり順序。PAI2のリスト順序/MasterPeakIDに対応する設計 |
| `scan_id` | int | MSDec結果のスキャンID |
| `raw_spec_id` | int | 元のraw spectrum ID。PAI2の `ms2_raw_id` と対応し得る |
| `precursor_mz` | float | プリカーサーイオン m/z |
| `ion_mode` | int | 生のイオンモード列挙値。`0=Positive`, `1=Negative`, `2=Both` |
| `rt` | float | プリカーサーピークRT（min） |
| `model_peak_height` | float | デコンボリューションモデルピーク高さ |
| `signal_to_noise` | float | DCL scoring blockのS/N。PAI2のピーク検出S/Nとは別フィールド |
| `estimated_noise` | float | DCL scoring blockの推定ノイズ |
| `n_msms_peaks` | int | 元スペクトルに格納されたフラグメントピーク数 |
| `msms_spectrum` | list[list] | `[fragment_mz, intensity]` のリスト |

`include_spectrum=false` では `msms_spectrum=[]` だが、`n_msms_peaks` は元の本数を保持する。`top_n_peaks=N` では `msms_spectrum` だけを強度上位N本に縮め、m/z昇順に戻す。したがって **`len(msms_spectrum)` と `n_msms_peaks` は一致しないことがある**。

### 6.2 `summarize_dcl()`

| キー | 意味 |
|---|---|
| `total_results` | MSDecResult総数 |
| `with_msms` | `n_msms_peaks > 0` の件数 |
| `msms_rate_pct` | `with_msms / total_results * 100` |
| `msms_peak_count_median` | MS/MS保有結果におけるフラグメント数中央値。偶数件でも中央2値平均ではなく上側の中央要素 |
| `msms_peak_count_max` | フラグメント数最大値 |
| `precursor_mz_range` | 正の precursor m/z の `(min, max)`。小数4桁丸め |
| `rt_range` | RT の `(min, max)`。小数2桁丸め |
| `error` | 空入力時のメッセージ |

### 6.3 検索とPAI2への付与

`get_msms_by_precursor()` は `abs(result.precursor_mz - query) <= tol`、任意で `abs(result.rt - query_rt) <= rt_tol`、かつ `n_msms_peaks > 0` のDCL辞書をそのまま返す。

`attach_msms_to_features()` は同じリスト位置の PAI2 と DCL を対応付け、m/z差が `mz_tol` 以下ならPAI2辞書へ次を追加/上書きする。

| 追加キー | 意味 |
|---|---|
| `msms_spectrum` | DCL由来フラグメント配列 |
| `n_msms_peaks` | DCLに記録された元フラグメント数 |
| `msms_peak_count` | `n_msms_peaks` と同じ値 |

関数の返り値は、付与に成功したうち `n_msms_peaks > 0` だったピーク数である。

### 6.4 MCPツール

#### `dcl_parser(file_path=None, top_n_peaks=None, preview=5)`

1つの `.dcl` を解析し、`summary`（§6.2 と同一キー）と `preview`（先頭 `preview` 件）を JSON で返す。`preview[]` の各要素は `dcl_index` / `precursor_mz` / `rt` / `n_msms_peaks` / `signal_to_noise` / `top_fragments`（`[[mz, intensity], ...]`）。

- `top_n_peaks` の既定は **10**（全ピークは容易に数百本になり文脈を食うため）。全件が必要なときだけ `0` を渡す。
- `n_msms_peaks` は**間引き前の元本数**であり `top_fragments` の長さとは異なり得る（共通核 §9-9）。
- `file_path` 省略時は最新バッチの `.dcl` を自動選択する（他パーサと同じ `_select_latest_batch`）。`.dcl` は測定ファイル1つにつき1個なので、**特定サンプルの MS/MS が欲しいときは `file_path` を明示すること**。

#### `dcl_find_msms(precursor_mz, file_path=None, rt=None, mz_tol=0.01, rt_tol=0.2, top_n_peaks=None)`

precursor m/z（任意で RT）に一致する MS/MS を引く。アノテーションの裏取り（期待されるフラグメントが実際に出ているか）に使う。RT を渡すと同一 m/z の別溶出ピークを分離できる。該当なしは `status="not_found"` を返し、「MS/MS 未取得のピークである可能性を同定確度に反映せよ」と促す——**該当ゼロを「フラグメントが無い＝別化合物」と解釈してはならない**（そもそも MS/MS が取られていない場合と区別できない）。

### 6.5 PAI2 ピークへの MS/MS 付与

`pai2_parser` は解析時に `lipidmix/pai2/tools.py` の `_attach_sibling_msms()` を通じて同名 `.dcl` を読み、`lipidmix/dcl/reader.py` の `attach_msms_to_features()` で索引対応（`dcl_index` ↔ PAI2 のピーク順）により `msms_spectrum` / `n_msms_peaks` を書き込む。安全弁として **precursor m/z が feature の `m/z` と 0.01 を超えて食い違う場合は付与しない**（索引対応が崩れている疑いがあるため、誤った MS/MS を貼るより欠測にする）。

結果は `pai2_parser` の応答 `summary.msms_attachment` に出る。

| キー | 意味 |
|---|---|
| `attached` | MS/MS を付与できたピーク数 |
| `dcl_file` | 読んだ `.dcl` のファイル名。隣接ファイルが無ければ `null` |
| `caveat` | `.dcl` 不在／読み込み失敗／付与0件のときの説明。`null` なら正常 |

`.dcl` が無い・壊れている場合でも PAI2 の在庫要約は成立するため解析は続行する。ただし caveat が出ている状態を「MS/MS 無し」と読み替えないこと（**未確認**であって**不在**ではない）。
