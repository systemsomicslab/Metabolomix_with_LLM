# ワークフロー: 定量（内部標準の正規化・極性の統合・特徴の選別）

出力の意味は [quant.md（output_format）](../output_format/quant.md)（`lipidmix://docs/output-format/quant`）。規則は外部の
規則プラグイン（`[quant] rules`。面の仕様は [../quant-rules-plugin.md](../quant-rules-plugin.md)）が持ち、未設定なら 4 ツールとも
`quant_rules_unavailable` の封筒を返す（ほかのツールは影響を受けない）。

順序: `dataset_quant_preview`（確認）→ `dataset_quant_apply`（正規化。`skip=True` で飛ばせる）→（陽・陰の 2 データセットなら）
`dataset_merge_polarities` → `dataset_select_features`。正規化しないなら apply を通さずに merge へ進んでよい（未正規化どうしとして
統合する）。どの工程も元のハンドルは変えず、新しいハンドルの派生データセットを作って直近にする（`metabolomix/quant/derive.py` の
`register()`）。

`dataset_handle` を持つ 3 ツールは `handle_dispatch.install` が包み、そのハンドルの `lock` の下でワーカースレッドに回す。
`dataset_merge_polarities` は `dataset_handle` も `session.dataset` も持たず、`handle_dispatch.ALWAYS_THREADED` に載せてセッションの
ロックを取らずにワーカースレッドで回る（`SESSION_READERS` には載せない）。2 つのハンドルのロックは `handle_dispatch.acquire_all()` が
`handle_id` の順に取り、待つのは全部で最長 `handle_dispatch.HANDLE_LOCK_WAIT`（25 秒）。待ちきれなければ取ったロックを放し、何もせずに
`HANDLE_BUSY`（`mcp_errors.handle_busy()`。ロックを持つ長い道具の仕事があればその `job_id` と道具）。

4 ツールとも入力は MS-DIAL の解析フォルダから読んだデータセットだけ（`metabolomix/quant/table.py` の `require_msdial_source()`。
ほかは `QUANT_SOURCE_UNSUPPORTED`）。記録を書く工程は、書く前に `derive.py` の `check_capacity()` で新しいハンドルを足せるかを確かめる。

## dataset_quant_preview

前提: データセット（無ければ `missing_state`）。陽陰を統合したデータセットは極性を解く前に `check_not_merged()` が
`QUANT_ALREADY_MERGED` で断る（`_check_not_normalized()` も同じ検査を通す）。極性が決まらなければ `polarity` が要る
（`QUANT_POLARITY_REQUIRED`）。状態変更: なし（記録も書かない）。

1. metabolomix/quant/tools.py  dataset_quant_preview()
2. ├─ metabolomix/quant/tools.py  _context()
3. │  ├─ metabolomix/tools/_holder.py  resolve_holder()
4. │  ├─ metabolomix/quant/normalize.py  check_not_merged()
5. │  └─ metabolomix/quant/table.py  resolve_quant_polarity()
6. ├─ metabolomix/quant/tools.py  _reference()
7. │  └─ metabolomix/quant/inputs.py  reference_input()
8. └─ metabolomix/quant/normalize.py  preview_tool()
9.    ├─ metabolomix/quant/normalize.py  preview()
10.    │  ├─ metabolomix/quant/normalize.py  _compute()
11.    │  │  ├─ metabolomix/quant/plugin.py  load_rules()
12.    │  │  ├─ metabolomix/quant/normalize.py  _check_not_normalized()（解析フォルダ由来・未統合・未選別・未正規化）
13.    │  │  ├─ metabolomix/quant/table.py  feature_rows()
14.    │  │  └─ metabolomix/quant/plugin.py  call()
15.    │  └─ metabolomix/quant/normalize.py  _preview_dict()
16.    └─ metabolomix/quant/inputs.py  quant_error()

## dataset_quant_apply

前提: preview と同じ。`reference_path` は `skip=True` のときだけ省ける。正規化済み・選別済み・統合済みには掛けない
（`QUANT_ALREADY_NORMALIZED` / `QUANT_ALREADY_SELECTED` / `QUANT_ALREADY_MERGED`。`skip()` も `_check_not_normalized()` を通す）。
状態変更: 解析フォルダの `reports/quant/` に記録を 1 件書き、新しいハンドルを登録して直近にする。内部標準が未解決
（`needs_input`）なら何も書かない。`derive()` は `carry_exploratory()` で親の探索専用の印を引き継ぎ、親の結果・行列・図
（`figures`）・捨てた図の印（`discarded_figures`）を空にする（統合・選別も同じ。output_format §5）。

1. metabolomix/quant/tools.py  dataset_quant_apply()
2. ├─ metabolomix/quant/tools.py  _context()
3. ├─ metabolomix/quant/tools.py  _reference()
4. └─ metabolomix/quant/normalize.py  apply_tool()
5.    ├─ [skip=False] metabolomix/quant/normalize.py  apply()
6.    │  ├─ metabolomix/quant/normalize.py  _compute()（連鎖は dataset_quant_preview）
7.    │  ├─ metabolomix/quant/normalize.py  _preview_dict()
8.    │  ├─ metabolomix/quant/normalize.py  analysis_folder()
9.    │  ├─ metabolomix/quant/records.py  write_record()
10.    │  └─ metabolomix/quant/derive.py  derive()
11.    ├─ [skip=True] metabolomix/quant/normalize.py  skip()
12.    │  ├─ metabolomix/quant/normalize.py  _check_not_normalized()
13.    │  ├─ metabolomix/quant/records.py  write_record()
14.    │  └─ metabolomix/quant/derive.py  derive()
15.    └─ metabolomix/quant/derive.py  register()

## dataset_merge_polarities

前提: 2 つのハンドルのどちらにもデータセットがあり、正規化の状態（normalized・method・unit）が同じ（`_check_states()`。
apply を通していない側は未正規化〔`skip=True` と同じ値〕。違えば `QUANT_MERGE_STATE_MISMATCH`）。各側の通り方
（`applied` / `skipped` / `not_run`）は `normalize_status()` が決め、記録の `settings.normalize` と戻り値の `normalize` に書く。
試料名が食い違えば `needs_input`（`only_in`）、記録の置き場が決まらなければ `needs_input`（`record_folder`）。
状態変更: 2 つの解析フォルダの共通の親（または `record_folder`）の `reports/quant/` に記録を書き、新しいハンドルを登録する。

- 試料の表とシート（`sample_metadata_rows`）は `derive()` が陽イオンのものを引き継ぐ。陰イオンの試料ごとのバッチ・注入順は
  手順 16 で陰イオンの試料の表（当てたシートを含む）から読んで `quant_state.polarity_runs.negative` に、統合前の試料名は
  `quant_state.source_sample_names` に残す。統合したデータセットが陽イオンの測定として読む値（手順 22）と同じかを手順 23 で
  `negative_runs` にするので、記録は統合したデータセットを作ってから書く（手順 24。記録 ID はそのあと `quant_state` に入れる）。
  統計エンジンの前処理はこれでドリフト補正を極性ごとに当てる（[dataset_analysis.md](dataset_analysis.md) の `dataset_preprocess`）。
- アッセイ種別は `merged_assay_kind()`（手順 21）が両側から決める（`derive()` が写す陽イオンの値のままにしない。`derive()` が写す
  `analysis_dir` は陽イオンのものなので、統合したハンドルへの `record_objective` は陽イオンの解析フォルダの目的に書く）。
- 除外は `merged_exclusions()`（手順 18）が両側から統合後の ID に写す。`derive()` は持ち越す ID の除外だけを残すので、特徴の除外は
  すべて落ち、試料の除外は陽イオンの側のものだけが共通 ID へ改名されて残るため、作り直す。
- 決め方・写し方と、統合し直しで引き継がれないもの（統合したハンドルに当てたシートと除外）は output_format の §3 の
  `polarity_runs` と §5 の merge。

1. metabolomix/quant/tools.py  dataset_merge_polarities()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/quant/merge.py  merge_tool()
4.    ├─ metabolomix/core/handle_dispatch.py  acquire_all()（2 つのハンドルのロック。全部で最長 `HANDLE_LOCK_WAIT`。待ちきれなければ `HANDLE_BUSY`）
5.    ├─ metabolomix/quant/merge.py  read_mapping()
6.    ├─ metabolomix/quant/table.py  resolve_quant_polarity()
7.    ├─ metabolomix/quant/merge.py  merge()
8.    │  ├─ metabolomix/quant/plugin.py  load_rules()
9.    │  ├─ metabolomix/quant/merge.py  _check_polarities()
10.    │  ├─ metabolomix/quant/merge.py  _check_states()
11.    │  │  └─ metabolomix/quant/merge.py  normalize_status()（食い違いの details と記録・戻り値）
12.    │  ├─ metabolomix/quant/merge.py  _match()
13.    │  │  └─ metabolomix/quant/plugin.py  call()
14.    │  ├─ metabolomix/quant/merge.py  _check_groups()
15.    │  ├─ metabolomix/quant/records.py  merged_record_folder()
16.    │  ├─ metabolomix/analysis/stats/data.py  sample_runs()（陰イオンの試料の測定。統計エンジンと同じ読み方）
17.    │  ├─ metabolomix/quant/derive.py  derive()（陽イオンの試料の表・シートを引き継ぐ）
18.    │  ├─ metabolomix/quant/merge.py  merged_exclusions()（両側の除外を統合後の ID に写す）
19.    │  ├─ metabolomix/quant/derive.py  carry_exploratory()（両側の探索専用の印を極性付きで合わせる）
20.    │  ├─ metabolomix/quant/merge.py  merged_source_verification()（出所の裏取りは陽陰の弱いほう）
21.    │  ├─ metabolomix/quant/merge.py  merged_assay_kind()（アッセイ種別を両側から）
22.    │  ├─ metabolomix/analysis/stats/data.py  sample_runs()（統合したデータセットが陽イオンの測定として読む値）
23.    │  ├─ metabolomix/quant/merge.py  _runs_match()（陰イオンの注入順・バッチがそれと同じか → `negative_runs`）
24.    │  └─ metabolomix/quant/records.py  write_record()
25.    └─ metabolomix/quant/derive.py  register()

## dataset_select_features

前提: データセット（統合前でも後でもよい）。付加体の採用表（`[quant] adducts`）が無ければ
`needs_input`（`skip_adduct_filter=True` で付加体の選別だけ飛ばせる）。状態変更: `reports/quant/`
に記録と `<record_id>.audit.json` を書き、新しいハンドルを登録する。

1. metabolomix/quant/tools.py  dataset_select_features()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/quant/select.py  select_tool()
4.    ├─ metabolomix/quant/inputs.py  adduct_table()
5.    ├─ metabolomix/quant/table.py  resolve_quant_polarity()
6.    ├─ metabolomix/quant/select.py  select()
7.    │  ├─ metabolomix/quant/plugin.py  load_rules()
8.    │  ├─ metabolomix/quant/table.py  feature_rows()
9.    │  ├─ metabolomix/quant/select.py  _curation_groups()
10.    │  │  ├─ metabolomix/quant/select.py  _flags_for()
11.    │  │  │  └─ metabolomix/curation/apply.py  flags_for_arf2()
12.    │  │  ├─ metabolomix/quant/table.py  curation_digest()
13.    │  │  └─ metabolomix/quant/table.py  apply_curation()
14.    │  ├─ metabolomix/quant/select.py  _decide()
15.    │  │  └─ metabolomix/quant/plugin.py  call()
16.    │  ├─ metabolomix/quant/records.py  write_record()
17.    │  ├─ metabolomix/quant/records.py  write_audit()
18.    │  └─ metabolomix/quant/derive.py  derive()
19.    └─ metabolomix/quant/derive.py  register()
