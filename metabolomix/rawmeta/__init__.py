"""生の測定ファイルのメタデータと、実験情報シートの前段（stdlib だけ）。

- `acquisition`: raw の**取得時刻のメタデータだけ**を読む（.wiff の OLE・Thermo .raw のファイルヘッダ・
  Bruker .d の tdf/sqlite・mzML の startTimeStamp）。スペクトルは解さない。取得の順・40 分を超える
  間隔の数・群と取得順の重なりの判定もここ。
- `sample_table`: ISA の試料の表（`s_*.txt` / `a_*.txt`）と mwTab（`SUBJECT_SAMPLE_FACTORS`）を読み、
  因子の列付きの実験情報シート（sample-manifest.v2）を組む。

MCP にもセッションにも依存しない（`list_data_files` と `pipeline_sample_sheet` が呼ぶ）。
"""
