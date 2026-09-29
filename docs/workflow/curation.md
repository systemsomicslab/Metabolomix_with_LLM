# ワークフロー: アラインメントのキュレーション

注釈付きスポットの証拠を集めて機械判別し、ビューアを書き出し、ユーザーのフラグを記録する経路。
値の**意味**は [curation.md（output_format）](../output_format/curation.md)。

```mermaid
flowchart TD
    CR[curation_review] --> LOAD[arf2.reader.load_catalog]
    LOAD --> SEL[curation.evidence.select_spots]
    CR --> RUN[curation.review.run_review]
    RUN --> COL[curation.evidence.collect]
    COL --> ANN[arf2.match_results.load_spot_annotations]
    COL --> REF[library.store.LibraryStore.record_by_scan_id]
    COL --> EIC[eic.reader.read_eic_spot_css1]
    COL --> SHP[curation.eic_shape.spot_shape]
    COL --> MS[analysis.spectral_match.match_spectrum]
    COL --> MIR[plots.mirror.build_mirror_payload]
    RUN --> TR[curation.trend.fit_trends]
    RUN --> JG[curation.judge.judge_spot]
    CR --> SAV[curation.review.save_review]
    SAV --> HTML[curation.viewer.render_html]
    CS[curation_submit] --> PAR[curation.flags.parse_submission_text]
    CS --> VAL[curation.flags.validate_entries]
    CS --> APP[curation.flags.FlagStore.append]
    CS --> WB[curation.msdial_writeback.sync_misannotation]
    WB --> TAG[msdial.tags.update_alignment_tag]
    CG[curation_suggest] --> LRV[curation.suggest.latest_review]
    CG --> RS[curation.suggest.run_suggestion]
    RS --> TG[curation.suggest.select_targets]
    RS --> IF[arf2.ion_features.load_ion_features]
    RS --> AP[msdial.analysis_params.resolve_analysis_params]
    RS --> COL
    RS --> LC[curation.candidates.build_library_candidates]
    RS --> FR[curation.relations.find_relations]
    CG --> SS[curation.suggest.save_suggestion]
    SS --> SH[curation.viewer.render_suggest_html]
    CS --> EX[curation.suggest.expand_entries]
```

## curation_review

前提: `.arf2` が解決できること（無ければ `missing_state("arf2_file", ["load_dataset", "arf2_parser"])`）、
`session.library.store` があること（無ければ `missing_state("library", ["library_load"])`）。
状態変更: `session.curation.last_review_id` と `review_dirs[review_id]`。ディスクに
`<arf2 のフォルダ>/curation/review-<id>.json` と `.html` を書く。

1. lipidmix/core/path_resolvers.py  resolve_arf2_file_path()
2. lipidmix/curation/judge.py  resolve_thresholds()
3. lipidmix/arf2/reader.py  load_catalog()
4. lipidmix/curation/evidence.py  select_spots()（load_catalog() が返したカタログを渡す）
5. lipidmix/curation/review.py  run_review()
6. │  └─ lipidmix/curation/flags.py  alignment_key()
7. │  └─ lipidmix/curation/flags.py  FlagStore.rows()（重い収集の**前**に読む。読めない行は `FlagFileError` → エラー payload）
8. │  └─ lipidmix/curation/flags.py  effective_flags()
9. │  └─ lipidmix/curation/flags.py  orphaned_count()（以前の版の `.arf2` に付いたフラグ → `warnings`）
10. │  └─ lipidmix/curation/flags.py  cleared_spots()（人が取り消したスポット → `flag_cleared`）
11. │  └─ lipidmix/curation/evidence.py  collect()
12. │     ├─ lipidmix/curation/evidence.py  _arf_rows()
13. │     │  └─ lipidmix/arf/reader.py  deserialize()
14. │     ├─ lipidmix/curation/evidence.py  _check_file_ids()（`.arf` に無い `file_ids` は `UnknownFileIdsError` → エラー payload）
15. │     ├─ lipidmix/arf2/match_results.py  load_spot_annotations()
16. │     ├─ lipidmix/dcl/reader.py  deserialize_dcl()
17. │     ├─ lipidmix/curation/evidence.py  _reference()
18. │     │  └─ lipidmix/library/store.py  LibraryStore.record_by_scan_id()
19. │     ├─ lipidmix/analysis/spectral_match.py  match_spectrum()
20. │     ├─ lipidmix/plots/mirror.py  build_mirror_payload()
21. │     ├─ lipidmix/eic/reader.py  read_eic_spot_css1()
22. │     ├─ lipidmix/curation/eic_shape.py  spot_shape()
23. │     └─ lipidmix/curation/evidence.py  _downsample_points()（形状計算の後に payload だけ間引く）
24. │  └─ lipidmix/curation/trend.py  composition()（`LipidParser` はモジュールで 1 つだけ作る）
25. │  └─ lipidmix/curation/trend.py  fit_trends()
26. │  └─ lipidmix/curation/judge.py  lipid_rules_active()（脂質規則フラグが 1 件でも True か）
27. │  └─ lipidmix/curation/judge.py  judge_spot()
28. │  └─ lipidmix/curation/judge.py  auto_note()（判定根拠の文 → `auto_note`）
29. lipidmix/curation/review.py  save_review()
30. └─ lipidmix/curation/viewer.py  render_html()
31. lipidmix/curation/review.py  n_summary_rows()
32. lipidmix/curation/review.py  summary_tsv()（`max_rows` で先頭だけ）
33. lipidmix/curation/review.py  trend_summary()

## curation_suggest

前提: `.arf2` が解決できること、`session.library.store` があること（無ければ
`missing_state("library", ["library_load"])`）、このアラインメントのレビューがあること（無ければ
`missing_state("curation_review", ["curation_review"])`）。状態変更: `session.curation.review_dirs[suggestion_id]`。
ディスクに `<arf2 のフォルダ>/curation/suggest-<id>.json` と `.html` を書く。

1. lipidmix/core/path_resolvers.py  resolve_arf2_file_path()
2. lipidmix/curation/judge.py  resolve_thresholds()
3. lipidmix/curation/flags.py  alignment_key()
4. lipidmix/curation/suggest.py  latest_review()（同じ sha256 のアラインメントに対する最新のレビュー。`review_id` を渡せば `_find_review()` でそれ）
5. lipidmix/curation/suggest.py  run_suggestion()
6. │  └─ lipidmix/curation/flags.py  FlagStore.effective()（有効な判断。読めない行は `FlagFileError` → エラー payload）
7. │  └─ lipidmix/arf2/ion_features.py  load_ion_features()（同位体・アダクト・インソースのリンク）
8. │  └─ lipidmix/msdial/analysis_params.py  resolve_analysis_params()（検索アダクトと RT 窓）
9. │  └─ lipidmix/curation/suggest.py  select_targets()（wrong フラグ・likely_wrong・未注釈から対象を選ぶ）
10. │  └─ lipidmix/curation/trend.py  fit_trends()（相手 Y = 注釈付きで対象でも wrong / redundant でもないスポットの傾向）
11. │  └─ lipidmix/curation/evidence.py  collect()（`keep_measured=True` で対象の証拠を集める）
12. │  └─ lipidmix/curation/candidates.py  build_library_candidates()（① MS-DIAL の下位候補 ② 閾値を緩めた再検索）
13. │  └─ lipidmix/curation/relations.py  find_relations()（④ 別スポットの同位体・アダクト・インソース断片としての説明）
14. lipidmix/curation/suggest.py  save_suggestion()
15. └─ lipidmix/curation/viewer.py  render_suggest_html()
16. lipidmix/curation/suggest.py  summary_tsv()（`max_rows` で先頭だけ）

`curation_submit` は `review_id` が `cs-…` のとき、`curation.flags.validate_entries` の代わりに
`curation.suggest.expand_entries` で候補 ID を記録行へ展開する。

## curation_submit

前提: `review_id` のレビューがディスクにあること。状態変更: `flags.jsonl` へ追記し、
アラインメントの `_tags.xml` の Misannotation を書き換える（書く前に `curation/tags-backup/` へ控え）。
レビューの探し先は、セッションの `review_dirs` → 送信用テキストの `arf2_path` →
引数 `file_path` → 既定の `.arf2` の順（サーバ再起動の後でも貼った文で送れる）。

1. lipidmix/curation/flags.py  parse_submission_text()（`submission_text` のとき）
2. lipidmix/curation/review.py  is_valid_review_id()（形が違えばパスに使う前にエラー）
3. lipidmix/tools/curation_tools.py  _find_review()
4. └─ lipidmix/curation/review.py  load_review()（候補フォルダを順に）
5. lipidmix/curation/flags.py  validate_entries()（`cs-…` のときは lipidmix/curation/suggest.py  expand_entries()）
6. lipidmix/curation/flags.py  alignment_key()（レビュー時の sha256 と一致しなければ拒否）
7. lipidmix/curation/flags.py  FlagStore.rows()（読めない行があれば追記せずにエラー）
8. lipidmix/curation/flags.py  FlagStore.append()
9. lipidmix/curation/msdial_writeback.py  sync_misannotation()（失敗しても記録は残し `tags_xml.error`）
10. └─ lipidmix/msdial/tags.py  update_alignment_tag()
11.    └─ lipidmix/core/atomic_io.py  atomic_write_bytes()

## curation_flags

1. lipidmix/curation/flags.py  alignment_key()
2. lipidmix/curation/flags.py  FlagStore.effective()

## curation_view_data

1. lipidmix/curation/review.py  is_valid_review_id()
2. lipidmix/tools/curation_tools.py  _find_review()
3. └─ lipidmix/curation/review.py  load_review()
4. lipidmix/curation/review.py  page()
