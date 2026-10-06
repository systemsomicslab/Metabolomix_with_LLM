# ワークフロー: `.arf2`（スポット代表カタログ）

サンプル別の生データを持たない軽量なメタ層。1 スポット = 1 行の代表値だけを持つため
**PCA はできない**。サンプル間比較が要るときは `.arf`（[arf.md](arf.md)）へ。

`session.arf2` は `session.arf` とは別スロット。`arf2_parser` を呼んでも進行中の
ARF 解析（前処理行列・差次的結果・手動除外）は壊れない。

## arf2_parser

前提: なし（`file_path` 省略時は最新バッチを自動選択）
状態変更: `session.arf2` にカタログを格納。ARF 側の解析基盤には触れない。

要約テキストの生成（手順 4）はセッションへの格納（手順 5）より先に走る。格納は
「将来の検索用に保持する」だけで、返り値の組み立てには使われないため。

1. metabolomix/arf2/tools.py  arf2_parser()
2. └─ metabolomix/core/path_resolvers.py  resolve_arf2_file_path()
3. └─ metabolomix/arf2/reader.py  deserialize()
4. └─ metabolomix/arf2/reader.py  generate_text_summary()
5. └─ metabolomix/core/session_state.py  Arf2State.load()
6. └─ metabolomix/core/session_state.py  AnalysisSession.maybe_prepend_caveat()

## arf2_annotate_identities

前提: なし（`.arf2` を直接読む。`arf2_parser` の実行結果には依存しない）
状態変更: なし

GOSLIN 正規化・RefMet / LIPID MAPS ID・MSI レベルをオフラインで付与する。
ARF2 には MS/MS 取得フラグも精密質量誤差も無いので、`has_msms=False`・バンド
`UNKNOWN` の保守評価になる。MSI レベルはクラス上限の見積もりであって、MS/MS 実測の
裏付けとは別物。個別ピークの確度は `verify_peak_annotation`（[pai2.md](pai2.md)）を見る。

参照表（手順 4）はモジュールレベルで 1 回だけ読み込んでキャッシュする。`.arf2` 本体も
`load_catalog()`（パス＋mtime＋サイズをキーにした共有キャッシュ）経由なので、差次的解析の
名前補完・EIC の同定照合と同じファイルを読み直さない。

返り値は列名を 1 回だけ出す TSV 表で、ヘッダ行に**総スポット数と未表示件数**を明記する
（返すのはファイル先頭から `max_rows` 件で、強度順でも MSI 順でもない）。

1. metabolomix/arf2/tools.py  arf2_annotate_identities()
2. └─ metabolomix/core/path_resolvers.py  resolve_arf2_file_path()
3. └─ metabolomix/arf2/reader.py  load_catalog()
4. └─ metabolomix/core/tool_helpers.py  _identity_tables()
5. │  └─ metabolomix/msdial/lipid_identity.py  load_reference_tables()
6. └─ metabolomix/msdial/lipid_identity.py  build_identity_block()
