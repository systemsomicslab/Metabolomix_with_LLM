# ワークフロー: `.dcl`（デコンボリューション済み MS/MS）

MS/MS スペクトルの実体（MS-DIAL の `MsdecResultsReader` 相当の独自バイナリ。msgpack/lz4 ではない。`dcl_index` が同名 `.pai2` の
ピーク順に一致する）。`dcl_find_msms` の `not_found` は MS/MS が取得されていないという意味で、フラグメントが無い根拠にはならない。
返すスペクトルは既定で強度上位 10 本（`top_n_peaks=0` で全件）。値の意味は [dcl.md（output_format）](../output_format/dcl.md)。
スポットの同定確度は、代表 MS/MS を `dcl_spot_spectrum` で引き、`library_msms_qualify` と `dataset_claim_evidence`
（[reanalysis.md](reanalysis.md)）で照合・合否を見る。

## dcl_find_msms

前提: `file_path` を渡すか、直近のデータセット（または `dataset_handle`）が `.dcl` を持つこと（省けば `handles` の `dcl_path`、
無ければ `artifact_paths` の `msms_evidence`。データフォルダは探さない）。データセットが無い・閉じたなら
`missing_state("dataset", ["load_dataset"])`、`.dcl` を持たなければ `no_alignment_files`。状態変更: なし。
RT を渡すと同一 m/z の別溶出ピークを分離できる。

1. metabolomix/dcl/tools.py  dcl_find_msms()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. ├─ metabolomix/tools/_holder.py  alignment_paths()
4. │  └─ metabolomix/dataset/artifacts.py  artifact_path()
5. ├─ metabolomix/dcl/reader.py  deserialize_dcl()
6. └─ metabolomix/dcl/reader.py  get_msms_by_precursor()

## dcl_spot_spectrum

前提: `arf2_path` / `dcl_file` を渡すか、直近のデータセット（または `dataset_handle`）が `.arf2` / `.dcl` を持つこと（省いたものは
そのデータセットのファイル。データフォルダは探さない。データセットが無ければ `missing_state("dataset", ["load_dataset"])`、
ファイルが無ければ `no_alignment_files`）。状態変更: なし。`library_match_feature(spot_id=...)` と同じ規則で代表 MS/MS を引く。

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
