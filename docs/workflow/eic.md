# ワークフロー: `.EIC.aef`（抽出イオンクロマトグラム）

描画系 2 ツールは**画像を作らない**。レンダラ中立の構造化 payload を返すだけで、
描画方法はクライアント（Use-LLLM は Plotly、Claude Desktop は各 UI の方式）に委ねる。
PNG が要るときだけ `save_figure(kind="eic")`（[plots.md](plots.md)）を明示的に呼ぶ。

`peak_top` は座標であって強度ではない。混同しないこと。

MS-DIAL 2023-11-20 より前の版が書いた `.EIC.aef` は、全スポットで m/z が `-1` になる
（上流 `DataObjConverter` の既知の不具合。RT は正しい）。全ツールは読み込み直後に
`mz_backfill.backfill_missing_mz()` を通し、同じバッチの ARF2 の `MassCenter` で補う
（同じ AlignmentID の RT 一致が条件）。補ったことは注記に出る。

```mermaid
flowchart TD
    ES[eic_search<br/>spot ID・m/z・RT の AND / 強度順] --> SPOT[spot_id を特定]
    SPOT --> PC[eic_plot_chromatograms<br/>lipidmix.eic.v1]
    SPOT --> PM[eic_plot_compounds<br/>lipidmix.eic.multi.v1]
    SPOT --> LM[library_match_feature spot_id=…<br/>→ library_plot_mirror（library.md）]
    PC --> SAVE[save_figure kind=eic → plots.md]
    PM --> SAVE
```

`eic_parser` と `eic_search` は `EicState.load_data()` を通るため、どれから呼んでも
未ロードなら読み込まれる（`MissingState` にはならない）。描画系 2 ツールは
`session.eic` を経由せず、必要なスポットだけをファイルから直接読む。

## eic_parser

前提: なし
状態変更: `session.eic` に EIC スポットを格納。

1. metabolomix/eic/tools.py  eic_parser()
2. └─ metabolomix/core/path_resolvers.py  resolve_eicaef_file_path()
3. └─ metabolomix/eic/tools.py  _load_session_spots()
4. │  ├─ metabolomix/core/session_state.py  EicState.load_data()
5. │  └─ metabolomix/eic/tools.py  _backfill_mz()
6. │     └─ metabolomix/eic/mz_backfill.py  backfill_missing_mz()
7. └─ metabolomix/eic/reader.py  summarize_eic_data()
8. └─ metabolomix/core/session_state.py  AnalysisSession.maybe_prepend_caveat()

## eic_search

前提: なし
状態変更: `session.eic` に EIC スポットを格納（未ロードなら読み込む）。

`spot_ids`・m/z 範囲・RT 範囲を AND で掛ける（何も指定しなければ全件）。
`sort="max_intensity"` は各試料のクロマトグラム最大強度の、そのまた最大値で降順
（RT 座標順ではない）。表示する行にだけ、同じバッチの `.arf2` の同定（手順 11）と
アラインメント `.dcl` の MS/MS ピーク数（手順 13–14）を付ける。旧
`eic_search_by_mz_range` / `eic_search_by_rt_range` / `eic_rank_by_max_intensity` を統合した。

1. metabolomix/eic/tools.py  eic_search()
2. └─ metabolomix/core/path_resolvers.py  resolve_eicaef_file_path()
3. └─ metabolomix/eic/tools.py  _load_session_spots()
4. │  ├─ metabolomix/core/session_state.py  EicState.load_data()
5. │  └─ metabolomix/eic/tools.py  _backfill_mz()
6. │     └─ metabolomix/eic/mz_backfill.py  backfill_missing_mz()
7. └─ metabolomix/eic/reader.py  search_eic_spots()
8. └─ metabolomix/eic/tools.py  _annotate_rows()
9. │  ├─ metabolomix/eic/mz_backfill.py  sibling_arf2_path()
10. │  ├─ metabolomix/eic/identity_map.py  load_arf2_records()
11. │  │  └─ metabolomix/arf2/reader.py  load_catalog()
12. │  ├─ metabolomix/msdial/alignment_spot.py  alignment_sibling()
13. │  ├─ metabolomix/msdial/alignment_spot.py  msms_peak_counts()
14. │  └─ metabolomix/msdial/alignment_spot.py  precursor_agrees()
15. └─ metabolomix/arf2/reader.py  format_spots_as_table()

## eic_plot_chromatograms

前提: 描画対象の `spot_id`（検索系ツールで特定する）
状態変更: `session` に `lipidmix.eic.v1` payload を記録。画像は作らない。

`file_ids` で最大 12 試料まで選べる。1 スポット分だけを CSS1 から直接読むので、
`session.eic` の全件ロードは経由しない。座標は 4 桁に丸めて返す（RT 0.0001 分・
m/z 4 桁で図にも解釈にも足り、float の既定 repr の 17 桁は空白同然のため）。

1. metabolomix/eic/tools.py  eic_plot_chromatograms()
2. └─ metabolomix/core/path_resolvers.py  resolve_eicaef_file_path()
3. └─ metabolomix/eic/reader.py  read_eic_spot_css1()
4. └─ metabolomix/eic/tools.py  _backfill_mz()
5. │  └─ metabolomix/eic/mz_backfill.py  backfill_missing_mz()
6. └─ metabolomix/plots/eic.py  build_eic_plot_payload()

## eic_plot_compounds

前提: 描画対象の脂質名/オントロジーと `file_id`（**1 試料分だけ**を重ねる）
状態変更: `session` に `lipidmix.eic.multi.v1` payload を記録。ファイルは書かない。

**既定は画像**（`output="image"`）。サーバ側で描いた PNG を、何を描き何を落としたかを
書いた 1 行のキャプションと一緒に返す。点列が要るクライアントは `output="payload"` か
env `LIPIDMIX_PLOT_OUTPUT=payload`（[plots.md](plots.md) 参照）。重ねる本数の既定は
8 本 —— それ以上は 1 枚の図として判読できず、payload も比例して膨らむ。

ARF2 の同定候補を rt/mz で照合し、外れた物質は `selection.dropped` に理由付きで残る。
照合（手順 12）は候補選抜の中ではなく payload 組み立ての中で走る —— 実際に読み出した
EIC スポットと突き合わせないと検証にならないため。

1. metabolomix/eic/tools.py  eic_plot_compounds()
2. └─ metabolomix/plots/render.py  resolve_plot_output()
3. └─ metabolomix/core/path_resolvers.py  resolve_eicaef_file_path()
4. └─ metabolomix/core/path_resolvers.py  resolve_arf2_file_path()
5. └─ metabolomix/eic/identity_map.py  load_arf2_records()
6. │  └─ metabolomix/arf2/reader.py  load_catalog()
7. └─ metabolomix/eic/identity_map.py  select_identity_candidates()
8. └─ metabolomix/eic/reader.py  read_eic_spots_css1()
9. └─ metabolomix/eic/tools.py  _backfill_mz()
10. │  └─ metabolomix/eic/mz_backfill.py  backfill_missing_mz()
11. └─ metabolomix/plots/eic.py  build_multi_compound_plot_payload()
12. │  └─ metabolomix/eic/identity_map.py  verify_spot_match()
