# ワークフロー: `.dcl`（デコンボリューション済み MS/MS）

MS/MS スペクトルの実体。msgpack/lz4 ではない独自バイナリ（MS-DIAL の
`MsdecResultsReader` 相当）で、`dcl_index` が同名 `.pai2` のピーク順に一致する。

`dcl_find_msms` の `not_found` は「そのプリカーサで MS/MS が取得されていない」という
意味であって、「期待したフラグメントが無い」ではない。同定確度を否定する根拠には
できない。

個々のスポットの同定確度に効かせたいときは、スポットの代表 MS/MS を `dcl_spot_spectrum` で
引き、`library_msms_qualify` と `dataset_claim_evidence`（[reanalysis.md](reanalysis.md)）で
照合・合否を見る。

返すスペクトルは既定で強度上位 10 本に間引く（全ピークは容易に数百本になるため）。
`top_n_peaks=0` で全件。

## dcl_find_msms

前提: `file_path` を渡すか、直近のデータセット（または `dataset_handle`）が `.dcl` を持つこと。
`file_path` を省くとそのデータセットの `.dcl`（`handles` の `dcl_path`、無ければ `artifact_paths` の
`msms_evidence`）を使う。**データフォルダは探さない**。データセットが無い・閉じたなら
`missing_state("dataset", ["load_dataset"])`、あるが `.dcl` を持たなければ `no_alignment_files`。
状態変更: なし

RT を渡すと、同一 m/z の別溶出ピークを分離できる。

1. metabolomix/dcl/tools.py  dcl_find_msms()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/tools/_holder.py  alignment_paths()
4. │  └─ metabolomix/dataset/artifacts.py  artifact_path()
5. ├─ metabolomix/dcl/reader.py  deserialize_dcl()
6. └─ metabolomix/dcl/reader.py  get_msms_by_precursor()

## dcl_spot_spectrum

前提: `arf2_path` / `dcl_file` を渡すか、直近のデータセット（または `dataset_handle`）が `.arf2` / `.dcl` を
持つこと（省いたものはそのデータセットのファイル。データフォルダは探さない。データセットが無ければ
`missing_state("dataset", ["load_dataset"])`、ファイルが無ければ `no_alignment_files`）。
状態変更: なし。`library_match_feature(spot_id=...)` と同じ規則でスポットの代表 MS/MS を引く。

1. metabolomix/dcl/tools.py  dcl_spot_spectrum()
2. └─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/tools/_holder.py  alignment_paths()
4. └─ metabolomix/library/tools.py  _spot_query()
5.    └─ metabolomix/msdial/alignment_spot.py  spot_record()
6.    └─ metabolomix/msdial/alignment_spot.py  alignment_msms()

## dcl_diagnostic_ions

前提: `dcl_spot_spectrum` と同じ（省いたパスは直近のデータセットのファイル）。
状態変更: なし。

1. metabolomix/dcl/tools.py  dcl_diagnostic_ions()
2. ├─ metabolomix/analysis/msms_rules.py  normalize_ions()
3. ├─ metabolomix/library/tools.py  _spot_query()
4. └─ metabolomix/analysis/msms_rules.py  find_ions()
