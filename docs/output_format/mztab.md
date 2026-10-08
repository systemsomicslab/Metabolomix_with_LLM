# mzTab-M 経路の出力フィールド

`dataset_load` が作る `DatasetState` と、その下流（特徴表・差次的エクスポート）が
出す値の意味。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

ツールが**どのファイルのどの関数をどの順に呼ぶか**は `docs/workflow/mztab.md`。
ここは**値の意味**だけを定義する。

## 13. mzTab-M 経路（`dataset_load` 以降）

mzTab-M 2.0.0-M は MS-DIAL Console が産む標準交換形式で、`DatasetState` の正準入力。
ARF 経路（`.arf` を直接読む）とは**別の ID 空間・別の粒度**なので、両者の値を
混ぜて読まないこと。

| | mzTab-M 経路 | ARF 経路 |
|---|---|---|
| 行の識別子 | `SMF_ID`（= `feature_id` / エクスポートの `spot_id`） | `MasterAlignmentID` |
| RT の単位 | ファイル上は**秒**（`retention_time_in_seconds`）。`DatasetState` は**分**に揃える | 分 |
| 同定の出所 | SME 行 / SML 行（§13.1） | `.arf2` の注釈 |
| MSI レベル | **常に空欄**（§13.5） | ヒューリスティック推定あり（`identity` トピック §12.3） |

### 13.1 同定の出所: SME と SML

mzTab-M では同定が 2 か所に出る。**意味が違うので混ぜてはいけない。**

| 出所 | 何か | `DatasetState` の置き場所 |
|---|---|---|
| SME 行（Small Molecule Evidence） | **スペクトル照合の証拠**。MS/MS が割り当てられた同定だけが載る | `feature_metadata` / `feature_candidates` |
| SML 行（Small Molecule） | 代表同定。MS1 だけの照合（Text DB）もここに載る | `feature_annotations` |

MS-DIAL は Text DB 由来の同定を **SME へ書かない**（`MztabFormatExport.cs` の
`ShouldWriteSmeLine` が `IsTextDbBasedRepresentative` を除外する）。したがって
Text DB 運用のメタボロミクスでは **SME が 0 行**になり、同定は SML にしか無い。
これは異常ではない。

**構造・名称は SMF 行には無い。** SMF が持つのは `SMF_ID` / `SME_ID_REFS` /
`exp_mass_to_charge` / `retention_time_in_seconds` / `abundance_assay[N]` だけで、
`database_identifier` / `smiles` / `inchi` / `chemical_name` は SME 専用の列。
SMF からこれらを読もうとすると常に `None` になり、「この測定には同定が無い」と
誤読される。

### 13.2 `feature_metadata` / `feature_candidates`（SME 由来＝証拠）

`feature_metadata[feature_id]` は rank 最上位の SME 1 件を投影したもので、
`name` / `mz` / `rt`（**分**）/ `inchikey` / `inchikey_source` / `smiles` / `inchi` を持つ。

`feature_candidates[feature_id]` は `SME_ID_REFS` が指す**全候補**を rank 昇順で
並べた一覧（`sme_id` / `rank` / `charge` / `exp_mass_to_charge` /
`theoretical_mass_to_charge` / `confidence_value` ＋ SME のテキスト列）。
実データでは 8 割の特徴が複数候補を持つので、**`feature_metadata` の 1 件を
「確定同定」と読まないこと**。rank 1 でも候補にすぎない。

### 13.3 `feature_annotations`（SML 由来＝ラベル。証拠ではない）

`feature_id`（= `SMF_ID`）→ 注釈の辞書。**証拠スロットではない**
（`feature_bindings` はこれを読まない。読ませると MS1 注釈だけで
`authentic_standard_match` が通ってしまう）。曖昧でない注釈は次のキーを持つ。

| キー | 由来（SML 列） | 備考 |
|---|---|---|
| `sml_id` | `SML_ID` | どの SML 行から来たか |
| `ambiguous` | — | 常に `False`（曖昧な場合は下記の別形） |
| `name` | `chemical_name` | |
| `database_identifier` | `database_identifier` | MS-DIAL は必ず `<db>:<name>` 形式で書く |
| `chemical_formula` | `chemical_formula` | |
| `smiles` | `smiles` | |
| `adduct` | `adduct_ions` | |
| `reliability` | `reliability` | 例: `annotated by user-defined text library` |
| `confidence_measure` | `best_id_confidence_measure` | |
| `confidence_value` | `best_id_confidence_value`（無ければ SME 由来の total score でフォールバック。下記） | float |
| `confidence_measures` | SME の `id_confidence_measure[2..8]` | dict。**どのツール戻り値にも出ない**（下記） |
| `inchikey` | 導出 | `smiles` 経由でのみ到達しうる（下記） |
| `inchikey_source` | 導出 | `smiles_derived` / `none` |

**`confidence_value` のフォールバック（Minor 10）**: SML 行が `best_id_confidence_value`
を持たない（実 mzTab で稀ではない: Text DB 経由でない同定はスコアが SME 側にしか
無い）ときは、対応する SME の個別スコア（`id_confidence_measure[2..8]` の合成
total score）を代わりに入れる。SML 側の値がある場合はそちらを優先し上書きしない。
このフォールバックが入ったことで、`dataset_export_differential` /
`arf_export_differential`（`metabolomix/analysis/feature_export.py`）の同じ列が、
以前は `null` だった行に SME 由来の値を出すようになった（列定義
`export_contract.py` 自体は無変更）。

**`confidence_measures` はどこからも読まれていない（Minor 9）**: `feature_annotations[fid]`
に SME の個別スコア（`id_confidence_measure[2..8]`、`_extract_confidence_measures`）が
入るが、`dataset_export_differential` / `arf_export_differential`
（`metabolomix/analysis/{dataset_export,feature_export}.py`）は列定義で明示された
キーだけを拾うため、この辞書は現状どのツール戻り値にも出ない。参照するには
`DatasetState.feature_annotations[fid]["confidence_measures"]` を直接読むこと。

**`inchi` キーは存在しない。** MS-DIAL は SML の `inchi` を常に `null` で書くため、
キーを置くと「取得していない」と「無い」の区別を偽ることになる。
`database_identifier` も必ず `<db>:<name>` 形式なので InChIKey にはならない。
したがって mzTab-M の SML から InChIKey に到達できるのは `smiles` 経由（RDKit）だけで、
テキスト DB が SMILES 列を持たなければ `inchikey` は `None` のままになる。

**対応付けの規則**:

- 1 つの SML が複数の SMF を指す（`SMF_ID_REFS` が `1|2`）場合は、同一分子が
  複数 adduct で検出された形。**全 feature に同じ注釈が付く**のが正しい。
- 複数の SML が 1 つの SMF を指す場合は、どれが正しいか決められないので名前を
  選ばない。その形は `{"ambiguous": True, "sml_ids": [...], "name": None}`。
- `SMF_ID_REFS` が存在しない特徴を指す行は捨て、警告を積む。

警告は**種類**で 1 件に集約する（行ごとに積むと実データで数百件になり、表示制限で
重要な警告が埋もれるため）。

### 13.4 `.annotations.tsv` の `identification_status`

特徴表（`export_features` が書く `<name>.annotations.tsv`）は §13.1 の区別を
`identification_status` で伝える。列（`_ANNOTATION_COLUMNS`）は 3 値のどれでも同じ。

| 値 | 意味 |
|---|---|
| `candidate` | SME 由来。スペクトル照合の候補（rank 1 でも確定同定は名乗らない）。候補の数だけ行が出る |
| `ms1_annotation` | **SML 由来。m/z 照合のみ。MS/MS の裏付けは無い**。1 特徴 1 行 |
| `unidentified` | 同定情報が無い。複数の SML が 1 特徴に当たって決められない場合もここ |

**`ms1_annotation` を `candidate` と同一視しないこと。** MSI レベルで言えば
MS/MS 照合とは別水準で、v2 の `required_evidence` では `library_match` も
`authentic_standard_match` も満たさない。`identity` トピック §12.4 の
「`FLAG_ONLY` を `PASS` と同等に扱わない」と同じ趣旨。

未同定の特徴も 1 行残る。表から消すと「同定できたものだけの世界」になるため。

### 13.5 差次的エクスポートの同定列

差次的エクスポート（`dataset_export_differential`）は**行単位で出所を 1 つに決める**。
決め手は InChIKey を供給した側で、列ごとに選ぶと name が SML 由来・inchikey が
SME 由来という食い違った行ができる。

| 列 | 値 | 意味 |
|---|---|---|
| `name_source` | `mztab_sme` | SME 由来（スペクトル照合の証拠） |
| | `mztab_sml` | **SML 由来。MS1 照合のみ** |
| `inchikey_source` | `database_identifier` / `inchi_derived` / `smiles_derived` | どこから導出したか |
| `msi_level` | **常に空欄** | 下記 |
| `spot_id` | `SMF_ID` | メタ行の `# id_space = mztab_smf_id` が ID 空間を宣言する |

`msi_level` が空欄なのは、`.arf2` 由来の MSI ヒューリスティック（`identity`
トピック §12.3）が Level 2 に MS/MS の取得を要件としており、別ルールの値を同じ列へ
入れると比較不能な 2 つの意味が同居するため。**空欄は「該当なし」ではなく
「この経路では取得していない」。**

**InChIKey が無い行は捨てられる**（`NO_ANNOTATED_FEATURES`）。下流
（別リポ massbank-context）が InChIKey を結合キーにするため。名前だけの SML 注釈は
このゲートを越えられない。

### 13.6 `inchikey_coverage`

`dataset_load` の要約と `dataset_status` が読む集計。

| キー | 意味 |
|---|---|
| `total_features` | SMF 行数 |
| `with_inchikey` | InChIKey が付いた特徴数（SME・SML 合算） |
| `by_source` | `database_identifier` / `inchi_derived` / `smiles_derived` / `none` の内訳。全特徴を 1 回ずつ数える |
| `rdkit_available` | RDKit で構造から導出できる環境か。`False` なら `smiles`/`inchi` 経路は静かに落ちる |
| `identified_by` | **同定の出所内訳**。`sme`（MS/MS 証拠あり）/ `sml_only`（MS1 注釈のみ）/ `none` |

`identified_by` は `dataset_load` の要約にも 1 行で出る。
**`sml_only` を `sme` と足して「同定できた件数」にしないこと** — 証拠水準が違う。

### 13.7 `abundance` は実測か gap-fill かを区別しない

mzTab-M の `abundance_assay[N]` は非ゼロでも実測ピークか gap-fill 補間かを
区別しない。実データ（60 サンプル × 714 特徴）では**セルの 70.0% が gap-fill**
だったので、非ゼロを検出と数えると検出率を 3 倍以上に過大評価する。

隣接する `.arf` が同一アライメントだと**数値で検証できた場合だけ**、
`ds.detected_mask`（(特徴 × サンプル) の bool 行列）に検出状態が入る。
要約は `ds.feature_qc`:

| キー | 意味 |
|---|---|
| `source` | `"arf"`（取り込めた）/ `None`（取り込めていない） |
| `reason` | `source is None` のときだけ。`no_candidate` ほか |
| `n_cells` / `n_detected` / `gap_filled_rate` | 取り込めたときの内訳 |
| `join` | 接合の検証内訳（スポット数一致・m/z 差） |

**`source is None` を「gap-fill が無い」と読まないこと。** 取り込めていないだけで、
その mzTab-M の非ゼロ値には補間値が混じり得る。この状態では検出率・欠測率を語れず、
`dataset_preprocess` の `min_detection_rate` も使えない（0 より大きい値を渡すと
引数エラーになる。黙って未検出 0 件として通さない）。

## `plot_pca_loadings(source="mztab")`

`dataset_pca` の結果は `plot_pca_loadings(source="mztab")` で主成分ごとのローディング図に描ける（`source="auto"` でも、候補がこの 1 つだけなら選ばれる）。
特徴量は前処理後の `pp_feature_names`、表示名は `feature_annotations`（SML 由来）か `feature_metadata`（SME 由来）の名前で、
無ければ `feature_id`。mzTab-M にはクラス（Ontology）が無いので棒は 1 色になる。r は `dataset_pca` が autoscale で回したときの
成分 × 特異値 / √n（`n_samples`）で出す。返り値の意味は `arf` トピックの §11.8 と同じ。
