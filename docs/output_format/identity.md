# 同定信頼度・標準化

名前正規化（GOSLIN）、参照表マッピング、MSI レベル推定。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

## 12. 同定信頼度・標準化（P2c）

`metabolomix/msdial/lipid_identity.py`（MCP非依存の純ロジック層、**完全オフライン**）と `metabolomix/msdial/peak_verification.py` の拡張が、脂質同定名の標準化と信頼度レベルの推定を担う。外部識別子の取得はネットワークを一切使わず、`pygoslin`（同梱・純Python）と同梱 TSV 表のみで行う。既存ツール・既定挙動・`verify_peak_annotation` の既存キーは不変で、新データは新ブロックに追加する。

### 12.1 GOSLIN 名正規化（`normalize_lipid_name`）

`pygoslin` で脂質ショートハンド名を正規化し、`normalized`（正規化名）/ `level`（構造レベル: SPECIES/MOLECULAR_SPECIES など）/ `lipid_maps_category` / `parse_ok`（失敗時 False＋`error`）/ `stripped`（下記の前処理をしたか）を返す。`pygoslin` 未導入や解析不能でも例外を投げず `parse_ok=False` を返す（グレースフルデグレード）。

**MS-DIAL 限定子接頭辞の除去**: MS-DIAL は Name に信頼度の限定子（`"no MS2: "` / `"low score: "` / `"w/o MS2:"` 等）を前置し、複数候補を `|` で連結する。解析前にこれら接頭辞と `|` 以降を除去して単一の species 表記に整える（除去した場合 `stripped=True`）。この修正で正規化名の付与が増える。`"RIKEN N-VS1 ID-…"` 等の真の未同定名は接頭辞除去後も解析不能のまま（`parse_ok=False`）で正しい。

**プラズマローゲン注意**: pygoslin は species レベルで `PC P-34:0` を `PC O-34:1` に**同一化**する（P-/O- エーテルの曖昧性）。P- と O- を区別したい場合は正規化名に頼らず、元の `Name` の表記と MS/MS の裏付け（`.dcl` のスペクトル。[[pe-p-vs-pe-o-annotation]] / [[plasmalogen-oxidation]] に直結）で読むこと。この曖昧さを注意として添えるツールは無いので、読み手が覚えておく。

### 12.2 同梱マッピング表（`load_reference_tables` / `map_to_reference`）

`reference/lipidmaps_classes.tsv`（クラス→LIPID MAPS カテゴリ/メインクラス）と `reference/refmet_map.tsv`（クラス→RefMet 名）は**キュレート済みの部分集合**（一般的な脂質クラスを網羅）。置き場所はリポジトリルート直下の `reference/` 固定（`lipid_identity.REFERENCE_DIR`）で、サーバの起動ディレクトリに依存しない。クラストークンで写像し、`matched`（bool）/ `lipid_maps_category` / `lipid_maps_main_class` / `refmet_name` / `caveat` を返す。表に無いクラスは `matched=False`＋caveat「同梱マッピング表に無いため ID 未付与」を返し、**推測はしない**。

### 12.3 MSI レベル推定（`msi_level`、ヒューリスティック）

決定論的シグナル（名称の有無・MS/MS取得・精密質量誤差バンド・アダクト整合バンド）を組み合わせて MSI 同定信頼度を推定する。`level`（2/3/4）/ `label` / `rationale` / `heuristic=True` を返す。

- **Level 2**（putative annotated compound）: 名称あり＋MS/MS取得＋精密質量整合（PASS）＋アダクト非FAIL。
- **Level 3**（putative class-level）: クラス（ontology）は判別できるが上記を満たさない。
- **Level 4**（unknown）: 名称・クラスとも無し。
- **Level 1（標準品照合）は決して主張しない**。返り値の `heuristic=True` が示すとおり、これは決定論的推定であって同定の確定ではない。

### 12.4 MS/MS の証拠の読み方と、`.arf2` の同定の標準化

旧 `verify_peak_annotation`（ピークごとの裏取りのドシエ）と `arf2_annotate_identities` は撤去した
（B-2b）。その `analytical_checks.msms` の 3 状態は、**ツールが返す値としては残っていない**が、
MS/MS を同定確度の根拠にするときの読み方の規則としては変わらない:

| 状態 | 意味 | MSI への効き方 |
|---|---|---|
| 実スペクトルあり（`.dcl` 由来） | `dcl_spot_spectrum` / `dcl_find_msms` でフラグメントの強度上位を確かめられる | Level 2 の根拠として有効 |
| フラグのみ（取得の参照はあるがスペクトル本体を読んでいない） | `.arf` のセルの MS/MS の有無（`has_msms` 系）が立つだけ | Level 2 にはなり得るが**根拠が弱い**。`.dcl` を読んで確かめる |
| どちらも無い | MS/MS を取得していない（または未確認） | Level 2 に到達しない |

**フラグのみを実スペクトルと同等に扱わないこと**——`.dcl` を読んでいないだけかもしれず、「MS/MS で
裏が取れている」とは言えない。参照ライブラリとの照合スコアは `library_match_feature` /
`library_msms_qualify`、主張の合否の数値は `dataset_claim_evidence`（`library` / `reanalysis` トピック）。

`.arf2` スポットカタログの注釈は、データセットを読むときに GOSLIN で正規化され、ID と MSI レベルが付く
（`feature_metadata["identity"]`。`dataset_feature_search` の行の `normalized` / `msi_level` 列に出る）。
**`.arf2` には MS/MS 取得フラグ・精密質量誤差が無いため MSI は保守的にクラス上限で評価される**
（`has_msms=False`、バンド UNKNOWN）。より確度の高い MSI の主張には `.dcl` のスペクトルと質量誤差を
別に確かめる。キュレーションの判断（`curation_submit` で付けた `wrong` / `suspect` / `assign:<記録名>` /
`redundant`）は `curation_flags` で読む。`.arf2` の sha256 が記録時と違う（MS-DIAL を再実行した）
フラグは当てない。

### 12.5 アダクト/元素表の拡張（`metabolomix/msdial/peak_verification.py`）

`ADDUCT_SHIFTS` を `(sign, shift, charge, n_mol)` の4タプル化し、多量体 `[2M-H]-`・多価 `[M-2H]2-`・`[M+FA-H]-`（`[M+HCOO]-` の別名）を追加。`adduct_mz` は `m/z = (n_mol×neutral + shift) / charge` で多量体・多価に対応する（既存1価アダクトの数値挙動は不変）。元素表に D(²H)/F/Br/¹³C を追加（標識・ハロゲン対応）。CCS/RT 参照照合・同位体パターン照合は参照表未同梱のため v1 対象外。

### 12.6 一般代謝物（`assay_kind=metabolite`）では §12 の大半が適用外

§12.1〜§12.4 は脂質名文法を前提にした層で、親水性・一般代謝物のアッセイでは成り立たない。

| 層 | `metabolite` での扱い |
|---|---|
| GOSLIN 正規化（§12.1） | **不適用**。`parse_ok=False` は「未同定」ではなく「脂質名ではない」 |
| RefMet/LIPID MAPS 表（§12.2） | 脂質のみ収録。非脂質は未収載で当然（欠落を同定失敗と読まない） |
| MSI レベル（§12.3） | 脂質クラス整合を前提にしたヒューリスティックなので**参照しない** |
| 質量・アダクト整合（§12.5） | 再利用可。ただし元素・アダクト表は固定表で全化学種を網羅しない |

代わりに同定の実体になるのは**候補集合**である。同一化合物が複数アダクト・複数極性・in-source fragment として別 feature に現れうるため、名前1件を確定同定として扱わない。mzTab-M 経路では SME が1特徴に複数行を持ち、`rank` 1 が代表にすぎない（`lipidmix://docs/output-format/mztab` §13）。異性体は分離条件で決まり、名前では決まらない。

## 13. mzTab-M 経路の同定は別トピック

この §12 は `.arf2` / `.pai2` 由来の注釈（脂質名・MSI ヒューリスティック）を扱う。
**mzTab-M 経路の同定は規則が違う**——同定の出所が SME 行と SML 行の 2 つに割れ、
`msi_level` は常に空欄になる（Level 2 の要件が噛み合わないため）。
定義は `lipidmix://docs/output-format/mztab` §13 にある。
