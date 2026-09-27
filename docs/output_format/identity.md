# 同定信頼度・標準化

名前正規化（GOSLIN）、参照表マッピング、MSI レベル推定、`verify_peak_annotation` の材料。

> 先に `lipidmix://docs/output-format`（共通核）を読むこと。行・列の粒度、脂質名文法、必須注意事項はそちらで定義され、ここでは繰り返さない。節番号は分割前の通し番号。

## 12. 同定信頼度・標準化（P2c）

`lipidmix/msdial/lipid_identity.py`（MCP非依存の純ロジック層、**完全オフライン**）と `lipidmix/msdial/peak_verification.py` の拡張が、脂質同定名の標準化と信頼度レベルの推定を担う。外部識別子の取得はネットワークを一切使わず、`pygoslin`（同梱・純Python）と同梱 TSV 表のみで行う。既存ツール・既定挙動・`verify_peak_annotation` の既存キーは不変で、新データは新ブロックに追加する。

### 12.1 GOSLIN 名正規化（`normalize_lipid_name`）

`pygoslin` で脂質ショートハンド名を正規化し、`normalized`（正規化名）/ `level`（構造レベル: SPECIES/MOLECULAR_SPECIES など）/ `lipid_maps_category` / `parse_ok`（失敗時 False＋`error`）/ `stripped`（下記の前処理をしたか）を返す。`pygoslin` 未導入や解析不能でも例外を投げず `parse_ok=False` を返す（グレースフルデグレード）。

**MS-DIAL 限定子接頭辞の除去**: MS-DIAL は Name に信頼度の限定子（`"no MS2: "` / `"low score: "` / `"w/o MS2:"` 等）を前置し、複数候補を `|` で連結する。解析前にこれら接頭辞と `|` 以降を除去して単一の species 表記に整える（除去した場合 `stripped=True`）。この修正で正規化名の付与が増える。`"RIKEN N-VS1 ID-…"` 等の真の未同定名は接頭辞除去後も解析不能のまま（`parse_ok=False`）で正しい。

**プラズマローゲン注意**: pygoslin は species レベルで `PC P-34:0` を `PC O-34:1` に**同一化**する（P-/O- エーテルの曖昧性）。P- と O- を区別したい場合は正規化名ではなく `peak_verification.ether_caveats()` の caveat（[[pe-p-vs-pe-o-annotation]] / [[plasmalogen-oxidation]] に直結）に依拠すること。

### 12.2 同梱マッピング表（`load_reference_tables` / `map_to_reference`）

`reference/lipidmaps_classes.tsv`（クラス→LIPID MAPS カテゴリ/メインクラス）と `reference/refmet_map.tsv`（クラス→RefMet 名）は**キュレート済みの部分集合**（一般的な脂質クラスを網羅）。置き場所はリポジトリルート直下の `reference/` 固定（`lipid_identity.REFERENCE_DIR`）で、サーバの起動ディレクトリに依存しない。クラストークンで写像し、`matched`（bool）/ `lipid_maps_category` / `lipid_maps_main_class` / `refmet_name` / `caveat` を返す。表に無いクラスは `matched=False`＋caveat「同梱マッピング表に無いため ID 未付与」を返し、**推測はしない**。

### 12.3 MSI レベル推定（`msi_level`、ヒューリスティック）

決定論的シグナル（名称の有無・MS/MS取得・精密質量誤差バンド・アダクト整合バンド）を組み合わせて MSI 同定信頼度を推定する。`level`（2/3/4）/ `label` / `rationale` / `heuristic=True` を返す。

- **Level 2**（putative annotated compound）: 名称あり＋MS/MS取得＋精密質量整合（PASS）＋アダクト非FAIL。
- **Level 3**（putative class-level）: クラス（ontology）は判別できるが上記を満たさない。
- **Level 4**（unknown）: 名称・クラスとも無し。
- **Level 1（標準品照合）は決して主張しない**。返り値の `heuristic=True` が示すとおり、これは決定論的推定であって同定の確定ではない。

### 12.4 統合ツール

- `verify_peak_annotation` のドシエに `identity_normalization` ブロック（`goslin` / `reference` / `msi` / `class_token`）を追加。既存の `analytical_checks`（精密質量誤差・アダクト整合・**MS/MS 証拠**）から算出したバンドを MSI 推定に流用する（質量・アダクトロジックの二重化を回避）。
- **MS/MS 証拠 `analytical_checks.msms`**（`peak_verification.msms_evidence`）。`band` は3値で、実スペクトルと取得フラグを峻別する:

  | `band` | `source` | 意味 | MSI への効き方 |
  |---|---|---|---|
  | `PASS` | `spectrum` | 実スペクトルあり（`.dcl` 由来）。`top_fragments` に強度上位が入る | Level 2 の根拠として有効。rationale に「実 MS/MS スペクトル確認」と明記 |
  | `FLAG_ONLY` | `flag` | `has_msms` は立つがスペクトル本体が無い | Level 2 にはなるが **rationale で「根拠が弱い」と開示**し、`dcl_parser` / `dcl_find_msms` での確認を促す。caveat にも出る |
  | `ABSENT` | `null` | フラグもスペクトルも無い | Level 2 に到達しない |

  `n_peaks` は間引き前の元本数、`top_fragments` は強度降順の上位5本（`[[mz, intensity], ...]`）。`llm_decision.deterministic_summary` にも `msms=<band>` が出る。**`FLAG_ONLY` を `PASS` と同等に扱わないこと**——`.dcl` を読んでいないだけかもしれず、「MS/MS で裏が取れている」とは言えない。

  `band == "PASS"` かつ参照ライブラリが読み込み済み（`library_load` 実行済み）のときは、この内側に `spectral_match` ブロックが追加で載る（参照ライブラリとの照合スコア）。3状態契約自体は変えない追加情報であり、`status` 語彙は暫定。定義は `lipidmix://docs/output-format/library` §14.8 を見ること。
- `arf2_annotate_identities(file_path=None, max_rows=50)`: ARF2 スポットカタログの注釈を一括で正規化・ID/レベル付与し、上位 `max_rows` 件を返す。**ARF2 には MS/MS 取得フラグ・精密質量誤差が無いため MSI は保守的にクラス上限で評価**（`has_msms=False`、バンド UNKNOWN）。より確度の高い MSI 評価は個別ピークの `verify_peak_annotation` を用いること。

### 12.5 アダクト/元素表の拡張（`lipidmix/msdial/peak_verification.py`）

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
