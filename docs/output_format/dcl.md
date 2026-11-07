# DCL (`.dcl` / デコンボリューション済み MS/MS)

`.dcl` パーサと `dcl_*` ツール（`dcl_spot_spectrum` / `dcl_diagnostic_ions` / `dcl_find_msms`）の出力定義。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

> **MS/MS の所在（重要）**: 実スペクトルを持つのは `.dcl` であって `.pai2` ではない。PAI2 の `has_msms` は**取得参照の有無**を示すだけで、`msms_spectrum` が非空である保証はない（共通核 §9-8）。同定確度（MSI Level 2）は MS/MS を根拠にするため、フラグだけで確度を主張してはならない。`.dcl` を読むツールは `dcl_*` で、データセットのアラインメントの `.dcl`（スポットの代表スペクトル）を既定に読む。

## 6. DCL (`metabolomix/dcl/reader.py`)

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

#### `dcl_find_msms(precursor_mz, file_path=None, rt=None, mz_tol=0.01, rt_tol=0.2, top_n_peaks=None, dataset_handle=None)`

precursor m/z（任意で RT）に一致する MS/MS を引く。アノテーションの裏取り（期待されるフラグメントが実際に出ているか）に使う。RT を渡すと同一 m/z の別溶出ピークを分離できる。該当なしは `status="not_found"` を返し、「MS/MS 未取得のピークである可能性を同定確度に反映せよ」と促す——**該当ゼロを「フラグメントが無い＝別化合物」と解釈してはならない**（そもそも MS/MS が取られていない場合と区別できない）。

- `top_n_peaks` の既定は **10**（全ピークは容易に数百本になり文脈を食うため）。全件が必要なときだけ `0` を渡す。`n_msms_peaks` は**間引き前の元本数**であり、返るフラグメントの数とは異なり得る（共通核 §9-9）。
- **`dataset_handle` を渡せばそのデータセットの、省けば直近のデータセットの `.dcl`（アラインメント結果の `.dcl`＝スポットの代表スペクトル）を読む**。データフォルダは探さない。データセットが無ければ `missing_state`。
- `.dcl` は測定ファイル 1 つにつき 1 個なので、**特定サンプルの MS/MS が欲しいときは `file_path` を明示すること**（省くとアラインメントの `.dcl` を読む）。

`dcl_spot_spectrum` と `dcl_diagnostic_ions` も同じ `.dcl` を `dataset_handle` / 直近のデータセットから読む（引数は `USAGE.md`）。

### 6.5 `.dcl` が読めないとき

`.dcl` が見つからない・壊れている場合に「MS/MS 無し」と読み替えないこと（**未確認**であって**不在**ではない）。`attach_msms_to_features()`（`metabolomix/dcl/reader.py`）は同じリスト位置の PAI2 ピークと DCL を対応付ける関数として残るが、PAI2 を読むツールは無い。その安全弁——**precursor m/z がピークの `m/z` と 0.01 を超えて食い違う場合は付与しない**（索引対応が崩れている疑いがあるため、誤った MS/MS を貼るより欠測にする）——は、対応づけを自前で行うときにも守る。
