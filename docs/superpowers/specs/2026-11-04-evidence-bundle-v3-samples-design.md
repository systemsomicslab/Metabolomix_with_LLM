# 根拠バンドル v3（試料ごとの値と検出）と新しい知見の候補の印 設計

- 日付: 2026-11-04
- 関連: massbank-context `docs/superpowers/specs/2026-11-04-xenobiotic-screen-design.md`（外来の分子の走査。このバンドルを読む）、
  reanalysis-study `docs/specs/2026-11-04-xenobiotic-signals-design.md`（手順書での使い方）
- 状態: 設計（未実装）

## 1. 目的

massbank-context に、薬剤・麻酔薬・溶媒など外から入った分子に関係する特徴を探す走査を足す。その 3 つの筋のうち
「出方」（処理群だけで検出される、親の薬剤と連動する）を見るには、**試料ごとの値と、その値が実測か補完か**が要る。
今の根拠バンドル v2（`dataset_export_evidence_bundle`）は群の平均（`mean_a` / `mean_b`）しか持たない。

あわせて、再解析の新しい知見の候補（`dataset_novel_candidates`）に、走査の結果の印を付けられるようにする。
v8 までの run では、薬剤由来とみられる特徴が候補の上位に並び、印付けが人の目（名前の推測）に頼っていた。

## 2. 目標・非目標・反目標

目標:

- `dataset_export_evidence_bundle` が v3 のバンドルを書く。v3 は v2 に `samples` と `intensities` の 2 役割と、
  features の任意の列 `adduct` を足したもの。
- `dataset_novel_candidates` が、特徴ごとの印（外部のファイル）を受け取り、候補の行に印の列を付け、印ごとの件数を出す。

非目標:

- 外来の分子の判定そのもの（massbank-context の仕事）。
- v2 の書き出しを残すこと。書き出しは v3 だけにする（読む側の massbank-context は 2 と 3 の両方を読む）。

反目標:

- 補完した値を実測と区別できない形で出す（「処理群だけで検出」を、補完の値の差で偽って作ってしまう）。
- 印の付いた候補を黙って除く（印は列で示すだけで、並びと件数は変えない）。

## 3. バンドル v3

### 3.1 manifest

`contract_version` を 3 にする。`exporter.version` は `evidence-v3`。`comparison` に `polarity`（`positive` / `negative`。
データセットの極性）を足す。それ以外の欄は v2 と同じ。

書き出す結果は呼び出しの `result_id` の 2 群比較そのもの（今と同じ）。再解析では工程 5 の主の比較（割れの組をまとめる前の
スポット単位）を渡す。特徴 ID はスポットの ID のままで、まとめた組への読み替えは呼び出し側（reanalysis-study）が `merged_from` で行う。

### 3.2 `samples.tsv`（役割 `samples`、必ず書く）

| 列 | 内容 |
|---|---|
| `sample_id` | 行列の試料名（`intensities` と同じ） |
| `group` | 比較の `group_a` の試料は `a`、`group_b` の試料は `b`、比較に入らない試料（QC・ブランク・群なし）は空 |
| `role` | 実験情報シートの `role`（`sample` / `qc` / `blank` など）。分からなければ空 |
| `injection_order` | 分かれば整数、分からなければ空 |

### 3.3 `intensities.tsv`（役割 `intensities`、必ず書く）

縦持ちで `feature_id, sample_id, value, detected`。

- `feature_id` は features の行と同じ。全特徴 × 比較に入る試料（`group` が `a` か `b`）を書く。
- `value` は、比較の統計に入った値と同じ前処理の段（正規化・変換の扱いを manifest の `comparison.preprocess` に書く）。
  **補完の前の値を持っているなら補完の前の値**、持っていなければ補完の後の値。どちらかを manifest に
  `intensities_stage: "pre_imputation" | "post_imputation"` で書く。
- `detected` は、検出状態（`.arf` から取り込んだ検出の印）があれば `true` / `false`（gap-fill と補完は `false`）、
  無ければ空（不明）。検出状態を取り込めなかった理由（`sample_axis_mismatch`・`mz_mismatch` など）は manifest の
  `population` に `detection_unavailable_reason` で書く。

規模の目安: 特徴 2.3 万 × 試料 12 で約 28 万行。問題にならない大きさとみるが、書き出しの時間とファイルの大きさを
テストで測って設計書に追記する。

### 3.4 features の `adduct` 列

MS-DIAL の付加体の表記（mzTab の SML/SMF の adduct）を任意の列として足す。無ければ空。

## 4. 新しい知見の候補の印

`dataset_novel_candidates` の spec に任意の `labels` を足す:

```json
"labels": {"path": "<TSV のパス>", "feature_column": "feature_id", "label_column": "level", "name": "xenobiotic"}
```

- TSV の 1 つの特徴に複数の行があれば、最も強い印を取る（強さの順は呼び出し側が `order` で渡す。既定は
  `["A","B","C","D"]`）。
- 候補の各行に `labels.<name>` の列を足し、要約に印ごとの件数（全候補と、上位 20 件の中）を出す。
- 並び・件数・除外（`exclude`）は変えない。
- 割れの組（`merged_from`）は、組のどれかに付いた最も強い印を取る。

## 5. 作りの決まり

- 書き出しは `metabolomix/analysis/evidence_export.py`。v2 の列定義はそのまま使い、新しい役割を足す。
- massbank-context との契約の文書は massbank-context 側（`docs/contracts/evidence-bundle-v3.md`）が正。ここの
  `export_contract` の版と食い違わないよう、両リポのテストで同じ見本のバンドルを読む（見本はこのリポの `tests/fixtures/` に置き、
  massbank-context 側は写しを持つ）。
- 出力の文書（`lipidmix://docs/output-format/reanalysis`）に v3 と `labels` を書く。

## 6. 検証

- 単体テスト: v3 の manifest・2 役割・`intensities_stage`・検出状態がある場合と無い場合・QC がある場合の `group` の空。
- `labels` のテスト: 複数行の特徴、割れの組、上位 20 件の件数、並びが変わらないこと。
- 実データ: ST003841・ST002901 の Console の出力から v3 を書き、massbank-context の走査が読めることを確かめる
  （massbank-context の設計書 §8 の陽性の対照）。
