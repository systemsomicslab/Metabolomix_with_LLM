# 脂質の種の目標（lipid_species）— 計画（2026-10-23）

再解析（reanalysis-study, MTBLS3233 C1〜C4）の主張は「PC 18:3_16:1」「PC 34:4」のような脂質の種で、
InChIKey が一意に決まらない（sn 位・二重結合の位置）。目標照合（`dataset_feature_search(targets=…)`）と
証拠の束（`dataset_claim_evidence`）はどちらも `inchikey` 必須なので、この型の主張を扱えない。

制約（ユーザー判断 2026-10-22）: 既存のツールの挙動・戻り値・既定値を変えない。opt-in の引数か新しいツールだけ。

## 決めたこと

1. **目標に新しいキー `lipid_species`**（`inchikey` と排他。どちらか一方が必須）。
   - GOSLIN で読み、**種のレベル**（`LipidLevel.SPECIES`）の文字列に直して照合に使う。
     `PC 18:3_16:1` → `PC 34:4`、`PC 16:0/18:1` → `PC 34:1`、`PE 16:0_18:1;O` → `PE 34:1;O`。
   - 修飾（`;O` など）は種の文字列に残るので `PE 34:1;O` ≠ `PE 34:1`。エーテルは別クラス
     （`PC O-34:1` ≠ `PC 34:1`）。GOSLIN は `P-` をビニルエーテルの二重結合を数えて `O-` に直す
     （`PC P-34:1` = `PC O-34:2`）——これは GOSLIN の種のレベルの規約で、そのまま使う。
   - 種より粗い（クラスだけ `PC`）・読めない名前は `FEATURE_TARGET_INVALID`（field `lipid_species`）。
   - `homolog.species_of` は `;O` を落とす（同族列の検査には十分だが同定には甘い）ので使わず、
     `metabolomix/msdial/lipid_identity.py` に `lipid_species_of(name)` を足す（MS-DIAL の限定子・
     `|` の除去は `normalize_lipid_name` と同じ `_clean_msdial_name`）。
2. **正規化後の目標**: 種の目標だけ `lipid_species`（種のレベル）と `lipid_species_input`（入力を GOSLIN で
   読んだ文字列。分子種のまま）を持ち、`inchikey_block` は None。InChIKey の目標には**キーを足さない**
   （正規化後の辞書・`rule_hash`・戻り値がバイト一致）。`rule_hash` は `lipid_species` を含み、
   `lipid_species_input` は `name` と同じく外す（`PC 18:3_16:1` と `PC 34:4` は同じ規則）。
3. **同定の照合**（`_identity`）: SME 候補の `chemical_name` を種のレベルに読み、目標の種と比べる。
   読めない候補は InChIKey の無い候補と同じく「評価できない」。不一致は `identifier_mismatch`。
   その後の付加体・価数・SML・層と水準の規則は InChIKey の目標と同じ。SML は単一注釈の `name` を読む。
   曖昧な SML（名前を持たない）は種の目標では照合できず `unannotated`。
4. **msms-rules.v1**: `target` に `lipid_species` を足す（`inchikey` と排他）。正規化後は
   `{"lipid_species", "formula"}`（InChIKey の規則は今のまま `{"inchikey_block", "formula"}`）。
   `target_in_top_rank` と `same_formula_margin` の「目標か別物か」をライブラリの記録の `name` の種で判定する。
   `library_msms_qualify` も同じ正規化を通るので、そのまま種の目標を受ける。
5. **証拠の束**: 主張の `target` の `lipid_species` を msms の規則の target へ渡す。

## テスト（先に書く）

- 種の一致（分子種の候補 ↔ 種の目標）・不一致（別クラス・同じ m/z）・エーテル・修飾・曖昧（2 特徴が
  qualified → `needs_input`）・SML の名前・読めない候補名・入力の誤り。
- InChIKey の目標の `rule_hash`・`match_target` の出力・msms の `rules_hash`・証拠の束の `rule_hash` を
  変更前の値で固定する。
- 証拠の束と `library_msms_qualify` を種の目標で通す（MCP ツール層）。
