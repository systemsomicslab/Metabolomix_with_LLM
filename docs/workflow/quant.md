# ワークフロー: 定量（内部標準の正規化・極性の統合・特徴の選別）

出力の**意味**（`quant_state`・gap-fill の印・単位・`needs_input`・記録と audit）は
[quant.md（output_format）](../output_format/quant.md)
（`lipidmix://docs/output-format/quant`）。規則は ms-data-parser の外の任意の規則プラグイン
（`[quant] rules`、面の仕様は [../quant-rules-plugin.md](../quant-rules-plugin.md)）が持つ。
未設定なら 4 ツールとも `quant_rules_unavailable` の封筒を返す（ほかのツールは影響を受けない）。

標準の順序は `dataset_quant_preview`（確認）→ `dataset_quant_apply`（正規化。`skip=True` で飛ばせる）
→（陽・陰の 2 データセットなら）`dataset_merge_polarities` → `dataset_select_features`。
どの工程も元のハンドルは変えず、**新しいハンドル**の派生データセットを作って直近にする
（`metabolomix/quant/derive.py` の `register()`）。

`dataset_handle` を持つ 3 ツールは `handle_dispatch.install` が包み、そのハンドルの `lock` の下で
ワーカースレッドに回す。`dataset_merge_polarities` は `dataset_handle` も `session.dataset` も持たないが、
`handle_dispatch.ALWAYS_THREADED` に載せてあり、セッションのロックを取らずにワーカースレッドで
回る（ロックの待ちでイベントループを塞がない）。2 つのハンドルのロックは `handle_id` の順に
自分で取る（`SESSION_READERS` には載せない）。

4 ツールとも入力は MS-DIAL の解析フォルダから読んだデータセットだけ
（`metabolomix/quant/table.py` の `require_msdial_source()`。ほかは `QUANT_SOURCE_UNSUPPORTED`）。
記録を書く工程は、書く前に `metabolomix/quant/derive.py` の `check_capacity()` で新しいハンドルを
足せるかを確かめる。

## dataset_quant_preview

前提: データセット（無ければ `missing_state`）。極性が決まらなければ `polarity` が要る
（`QUANT_POLARITY_REQUIRED`）。状態変更: なし（読み取り専用。記録も書かない）。

1. metabolomix/quant/tools.py  dataset_quant_preview()
2. ├─ metabolomix/quant/tools.py  _context()
3. │  ├─ metabolomix/tools/_holder.py  resolve_holder()
4. │  └─ metabolomix/quant/table.py  resolve_quant_polarity()
5. ├─ metabolomix/quant/tools.py  _reference()
6. │  └─ metabolomix/quant/inputs.py  reference_input()
7. └─ metabolomix/quant/normalize.py  preview_tool()
8.    ├─ metabolomix/quant/normalize.py  preview()
9.    │  ├─ metabolomix/quant/normalize.py  _compute()
10.    │  │  ├─ metabolomix/quant/plugin.py  load_rules()
11.    │  │  ├─ metabolomix/quant/table.py  feature_rows()
12.    │  │  └─ metabolomix/quant/plugin.py  call()
13.    │  └─ metabolomix/quant/normalize.py  _preview_dict()
14.    └─ metabolomix/quant/inputs.py  quant_error()

## dataset_quant_apply

前提: preview と同じ。`reference_path` は `skip=True` のときだけ省ける。正規化済みには掛けない
（`QUANT_ALREADY_NORMALIZED`）。状態変更: 解析フォルダの `reports/quant/` に記録を 1 件書き、
新しいハンドルを登録して直近にする。内部標準が未解決（`needs_input`）なら何も書かない。

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
12.    │  ├─ metabolomix/quant/records.py  write_record()
13.    │  └─ metabolomix/quant/derive.py  derive()
14.    └─ metabolomix/quant/derive.py  register()

## dataset_merge_polarities

前提: 2 つのハンドルのどちらにもデータセットがあり、`dataset_quant_apply` を通している
（normalized・method・unit が同じ）。試料名が食い違えば `needs_input`（`only_in`）で止まり、
記録の置き場が決まらなければ `needs_input`（`record_folder`）。状態変更: 2 つの解析フォルダの
共通の親（または `record_folder`）の `reports/quant/` に記録を書き、新しいハンドルを登録する。

1. metabolomix/quant/tools.py  dataset_merge_polarities()
2. ├─ metabolomix/tools/_holder.py  resolve_holder()
3. └─ metabolomix/quant/merge.py  merge_tool()
4.    ├─ metabolomix/quant/merge.py  read_mapping()
5.    ├─ metabolomix/quant/table.py  resolve_quant_polarity()
6.    ├─ metabolomix/quant/merge.py  merge()
7.    │  ├─ metabolomix/quant/plugin.py  load_rules()
8.    │  ├─ metabolomix/quant/merge.py  _check_polarities()
9.    │  ├─ metabolomix/quant/merge.py  _check_states()
10.    │  ├─ metabolomix/quant/merge.py  _match()
11.    │  │  └─ metabolomix/quant/plugin.py  call()
12.    │  ├─ metabolomix/quant/merge.py  _check_groups()
13.    │  ├─ metabolomix/quant/records.py  merged_record_folder()
14.    │  ├─ metabolomix/quant/records.py  write_record()
15.    │  └─ metabolomix/quant/derive.py  derive()
16.    └─ metabolomix/quant/derive.py  register()

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
