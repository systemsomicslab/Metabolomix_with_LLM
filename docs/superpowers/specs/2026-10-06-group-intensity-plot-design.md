# 群別強度プロット（`arf_plot_group_intensity`）設計

日付: 2026-10-06。設計はユーザー承認済み（会話内）。試作は `check.py`（`build_class_intensity_payload` /
`render_class_intensity_plot`）にあり、本実装で置き換えて `check.py` から消す。

## 1. 目的

選んだ脂質（代謝物）クラスや分子種が、指定した試料群のそれぞれで**見つかるか・どれくらいあるか**を、
報告書に載せられる図で示す。1 パネル = 1 項目、1 点 = 1 試料、群ごとに平均 ± SD。
検定はしない（群間の検定は `arf_differential` の役割。小 n の p 値が図に独り歩きするのを避ける）。

実例（20260930_EV）: 膜脂質 15 クラス × 群（Blank / EV_C / EV_P / Cell_C / Cell_PlnA）。
PE のように「どの群にも無い」ことも結果なので、該当なしの項目も N.D. パネルとして残す。

## 2. 範囲

- **ARF 経路のみ**（`load_dataset` / `arf_parser` で読んだ `.arf` と、同じアラインメントの兄弟 `.arf2`）。
  mzTab-M 経路（`dataset_*`）は対象外。
- 新しいツールは 2 つ: `arf_plot_group_intensity`（図を返す）と `save_group_intensity_figure`（ファイルに保存）。
  ツール数 72 → 74。

## 3. ユーザーが選ぶもの

### 3.1 `items`（必須、1〜30 件。1 項目 = 1 パネル、指定順）

項目の文字列を `+` で区切った各部分を次の規則で解決し、当たったスポットの和集合を 1 パネルにする。

- **クラス**: 部分が `.arf2` の `Ontology` のどれかと（大文字小文字を無視して）完全一致すれば、そのクラスの
  全スポット。例 `"PG"`、`"Cer_NS"`。
- **名前**: クラスに当たらなければ分子種名として扱う。`.arf2` の `Name` を `|` で分けた各候補名
  （先頭の確度接頭辞 `low score:` / `no MS2:` / `w/o MS2:` / `unsettled:` を除く）のどれかと、大文字小文字を
  無視して完全一致するスポット。例 `"PG 34:1"`（和組成）、`"PG 16:0_18:1"`（分子種）。
  アダクト違い・重複スポットで複数当たれば合計する。部分一致はしない（`"PG 34:1"` が `LPG …` 等に当たるのを防ぐ）。
- **合算**: `"PE+EtherPE"`、`"Cer_NS+Cer_AP"`、`"PG 34:1+PG 35:1"`。クラスと名前を混ぜてもよい。
- どの部分にも当たらない項目は **N.D.（該当なし）** のパネルになる。全項目が N.D. でもエラーにしない。

### 3.2 `groups`（必須、1 群以上。1 群 = 横軸の 1 列、指定順）

`arf_differential` の `group_a` / `group_b` と同じトークン規則（`metabolomix/msdial/sample_factors.py` の
`expand_sample_specs`）。`"KO"`、`"9w"`、`"KO_9w"`（両方を含む試料）。例 `["control", "KO"]`、
`["9w", "12M", "24M"]`。

- 既定では QC とブランクを群に入れない（差次的解析と同じ）。**群の指定自体にトークン `blank` / `qc` を含むときだけ**
  その役割を入れる（`["blank", "control", "KO"]` で基準として並べる用途）。
- 1 つの試料が複数の群に当たったら両方に描き、`caveats` に出す。
- 一致ゼロの群はエラー（`expand_sample_specs` の ValueError をそのまま伝える。利用可能トークンが載る）。
- `arf_exclude` で除外した試料は描かない（他の ARF ツールと同じ）。

## 4. 値の定義

- 1 点 = その試料での、項目に当たったスポットの **PeakHeight の合計**（`.arf` の試料別行。gap-fill の値も含める）。
- `gap_filled_fraction` = 合計のうち gap-fill の値が占める割合。0.5 を超える点は形を変える（◆）。
- 値が 0 の点は軸の底（log10 = 0）に ▽ で置く。
- 群の平均 ± SD は **log10 空間**で、値 > 0 かつ低信頼でない試料から計算（n = 1 なら平均だけ）。
- 縦軸は **全パネル共通**（log10(PeakHeight)、0〜最大値の切り上げ）。パネルごとに拡大すると、強度 10 前後の
  ノイズが信号に見えて「見つかるか」の判断を誤らせる（試作で実際に起きた）。
- 検出下限の破線: 引数 `detection_limit` > param ファイルの `Minimum peak height:`（`analysis_params` で読む。
  `.arf2` と同じフォルダの `*_param_*.txt`）> 無し。

## 5. 自動で除くもの

除いたスポットは理由ごとに戻り値 `excluded` に残す（名前と spot_id）。

1. **標識内部標準**: 名前に `(d<数字>)`（`PC 33:1(d7)`、`FA 16:0(d3)`）。**クラスとして当たった分だけ除く**。
   名前で直接指定した部分は描く（内部標準そのものを見たい用途）。
2. **キュレーション（`apply_curation=True` 既定）**: 有効な判断（`curation_flags` と同じ、スポットごとの最新 1 行、
   `clear` 済みを除く）。`wrong` / `redundant` は除外。`assign` は記録した名前・クラスで扱う（クラスの帰属も変わる）。
   `suspect` は描く。規則はエクスポート（`curation/apply.py`）と同じ。
3. **自動判定（`exclude_auto_likely_wrong=False` 既定）**: True なら、このアラインメント（`.arf2` の sha256 一致）の
   最新の `curation_review` の `likely_wrong` も除く。レビューが無ければ除かず `caveats` に出す。
4. **標準液にだけある分子種（`standard_samples` 任意）**: 指定トークンの試料（例 `["systemlot8"]`）での最大高さが、
   描く群の試料での最大高さの 10 倍以上のスポットを除く（クラスとして当たった分だけ）。システムチェック標準液の
   奇数鎖標準物質が試料のクラス合計に混ざるのを防ぐ（実例: 17:0_14:1 リン脂質・SM d18:1/12:0 等）。

## 6. その他の引数

- `low_reliability_samples`（任意、試料名またはトークン）: 白抜きで描き、平均 ± SD から外す（抽出量が少なかった等）。
- `title`、`ncols`（既定 `min(5, 項目数)`）、`output`（`image` 既定 / `payload`。`LIPIDMIX_PLOT_OUTPUT` と同じ規則）。

## 7. 戻り値

- `image`（既定）: PNG（dpi 100、`plots/render.figure_to_png`）と 1 行の説明（項目ごとのスポット数、群ごとの試料数、
  除外の件数）。
- `payload`: `plot_schema = "lipidmix.group_intensity.v1"`（形式名の接頭辞は外部契約として `lipidmix.` のまま）。

```text
{plot_schema, value: "sum of PeakHeight (log10 on the plot)", detection_limit, detection_limit_source,
 groups: [{label, samples: [name...]}], low_reliability_samples: [...],
 items: [{item, parts: [{part, kind: "class"|"name"|"none", n_spots}], n_spots, n_spots_msms,
          spots: [{spot_id, name, ontology}], detected,
          groups: [{label, samples: [{sample, value, gap_filled_fraction, low_reliability}],
                    log10_mean, log10_sd, n_in_stats}]}],
 excluded: {internal_standard: [...], curation: [...], auto_likely_wrong: [...], standard_only: [...]},
 caveats: [...]}
```

- `n_spots_msms`: 当たったスポットのうち MS/MS の裏付けがあるもの（`.arf2` の照合結果に MS/MS があり、名前が
  `no MS2:` でない）。0 のパネルは灰色にし「MS1-only (unconfirmed)」と書く（試作と同じ）。
- payload は `session.arf.last_group_intensity` に保持し、保存ツールが読む。

## 8. 保存ツール `save_group_intensity_figure(analysis_id, title=None)`

直前の payload から描き直し、`reports/figures/<slug>_group_intensity.png`（dpi 300）と同名の `.svg` を書く。
置き場所の解決は既存の `save_*_figure` と同じ（`metabolomix/tools/reports.py`）。payload が無ければ
`missing_state`（`required_tools: ["arf_plot_group_intensity"]`）。

## 9. コードの置き場所

- `metabolomix/plots/group_intensity.py`（MCP 非依存の純関数）:
  `resolve_items`（項目 → スポット集合と除外）、`build_group_intensity_payload`、`render_group_intensity_plot`。
- `metabolomix/msdial/analysis_params.py`: `Minimum peak height` を読むキーを足す。
- `metabolomix/arf/tools.py`: `arf_plot_group_intensity`（セッションの ARF・兄弟 `.arf2`・判断・レビューを集めて渡す）。
- `metabolomix/tools/reports.py`: `save_group_intensity_figure`。
- 前提状態が無いとき（ARF 未読み込み、兄弟 `.arf2` 無し）は `missing_state` の封筒。
- `check.py` の試作は消す。

## 10. テスト（fixture はテスト自身が作る）

- 項目の解決: クラス（大小無視）、名前（`|` の各候補・接頭辞除去・完全一致で部分一致しない）、`+` 合算、N.D.。
- 群の解決: 複数群・指定順、`blank` を書いたときだけ blank が入る、重なりの caveat、一致ゼロのエラー、`arf_exclude` の反映。
- 値: 合計、gap-fill 割合、0 の扱い、log10 平均 ± SD（低信頼の試料を外す、n = 1）。
- 除外: 内部標準（クラスでは除き名前指定では描く）、判断（wrong / redundant / assign の付け替え）、
  自動判定（レビュー有り・無し）、標準液だけの分子種。
- 描画: 縦軸の共通化、N.D. パネル、MS1-only パネル、検出下限の線（引数・param ファイル・無し）。
- ツール層: 画像と payload、`missing_state`、保存ツールのファイル出力。
- 登録数（72 → 74）、`USAGE.md`、`docs/workflow/`、出力形式の文書の各検査テスト。

## 11. 文書

- `USAGE.md`: 2 ツールの行。
- `docs/output_format/arf.md`（arf トピック）: payload の各フィールドの意味と、除外・N.D.・MS1-only の読み方。
- `docs/workflow/`（arf か plots）: 呼び出し順。
- vault の流れ図: ARF 経路の任意の描画ツールとして追加。

## 12. 対象外・将来

- mzTab-M 経路、群間検定、クラスの別名表（`PE` に `EtherPE` を自動で含める等。今は `+` で明示する）。
- 外部契約の名前（`LIPIDMIX_*`、`lipidmix://`、`lipidmix.*.vN`）を `metabolomix` に揃えるかは別途検討。
  本ツールの形式名も揃えるときに一緒に移す。
