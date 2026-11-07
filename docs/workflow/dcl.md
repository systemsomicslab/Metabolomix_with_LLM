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

前提: なし
状態変更: なし

RT を渡すと、同一 m/z の別溶出ピークを分離できる。

1. metabolomix/dcl/tools.py  dcl_find_msms()
2. └─ metabolomix/core/path_resolvers.py  resolve_dcl_file_path()
3. └─ metabolomix/dcl/reader.py  deserialize_dcl()
4. └─ metabolomix/dcl/reader.py  get_msms_by_precursor()

## dcl_spot_spectrum

前提: なし（`.arf2` は `arf2_path`、ハンドルならハンドルの `.arf2`、無ければ現在のデータフォルダ）。
状態変更: なし。`library_match_feature(spot_id=...)` と同じ規則でスポットの代表 MS/MS を引く。

1. metabolomix/dcl/tools.py  dcl_spot_spectrum()
2. └─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/library/tools.py  _spot_query()
4.    └─ metabolomix/msdial/alignment_spot.py  spot_record()
5.    └─ metabolomix/msdial/alignment_spot.py  alignment_msms()

## dcl_diagnostic_ions

前提: なし（`.arf2` は `arf2_path`、ハンドルならハンドルの `.arf2`、無ければ現在のデータフォルダ）。
状態変更: なし。

1. metabolomix/dcl/tools.py  dcl_diagnostic_ions()
2. ├─ metabolomix/analysis/msms_rules.py  normalize_ions()
3. ├─ metabolomix/library/tools.py  _spot_query()
4. └─ metabolomix/analysis/msms_rules.py  find_ions()
