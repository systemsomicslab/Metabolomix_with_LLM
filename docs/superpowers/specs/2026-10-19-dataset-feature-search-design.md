# DatasetState の特徴検索と、目的の物質の判定スポット（`dataset_feature_search`）

- 日付: 2026-10-19
- 状態: 設計（ユーザーは質問を省略して設計を依頼。判定方針（§3）と pipeline 結果の引き継ぎ（§11、案 A）はユーザー判断済み。書面レビュー（builder check）反映済み）
- 由来: reanalysis-study の run MTBLS6039 の教訓 `analysis/20261017_reanalysis_MTBLS6039_lessons.md`
  §3.1（目的の特徴を引く処理が毎回手書き）・§3.2（注釈の無い特徴の統計が引けない）。検討は `docs/HISTRY.md` 2026-10-19。

## 1. 目的

MTBLS6039 の再解析では、「論文が主張する代謝物 X は、この解析結果のどの特徴か。その特徴の差次的統計はいくつか」を
引く処理に既存の道具が無かった。オーケストレーターが 3 本のスクリプトを書き、実行役 35 人が同じ処理を
それぞれ書き直した（計 10,691 行）。質量窓・候補の絞り方・付加体の条件の解釈も実行役ごとに割れた。
さらに、注釈（InChIKey）の無い特徴は差次的エクスポートに出ないので、統計を再実装した（約 120 行と自己検査）。

このツールで次を満たす:

1. 同じ物質・同じ規則なら、誰が呼んでも同じ候補と同じ判定が返る（規則は引数で明示し、戻り値に写す）。
2. 注釈の有無に関係なく、候補の特徴の差次的統計（log2FC・p・q・群平均）が引ける。BH の母集団も明示する。
3. 規則で 1 件に決まらないときは、勝手に選ばずに「人（実行役）が決める」と返す。
   チームのトリアージ（教訓 §5.2）は、この戻りをそのまま使える。

## 2. 確認した事実（2026-10-19 時点のコード）

- **判定の部品は v2 にある**: `metabolomix/analysis/feature_bindings.py` の `bind_features` は、profile の `feature_targets` を
  特徴へ対応付ける。条件（m/z の ppm・RT・化合物 ID・付加体・電荷・証拠の水準）を AND で見て、候補ごとに落ちた理由を残し、
  qualified がちょうど 1 件のときだけ確定する（0 件／複数件はどちらも `needs_input` で、`reason` で区別）。
  強度や並び順では選ばない（docstring に明記）。
- **ただしそのままは使えない**:
  - pipeline v2 の `resolve_feature_bindings` 工程からしか呼ばれず、MCP から単独で呼べない。
  - 目標の形（`console/profile_schema.py` の `_validate_feature_target`）は RT が必須で、付加体は 1 つだけ。
    分子式から m/z を計算しない（`expected_mz` を直に要求する）。
  - `bind_features` は **全特徴**を評価して全候補を返す（陰で 11 万件）。対話ツールの戻り値には載せられない。
  - InChIKey の照合は正規化した全体一致だけ（先頭ブロック一致が無い）。
  - `_identifier_values(feature_meta, candidate)` は、特徴の InChIKey（rank1 の SME から解決した値）を**どの候補の値にも混ぜる**。
    rank2 以下の候補が別物質でも、rank1 の InChIKey で一致したことになる。
  - `compound_identifiers` に `formula` だけを書いた目標は、`_identifier_values` が分子式を値に入れないため、
    注釈付きの候補がすべて `identifier_mismatch` になる。
- **理論 m/z の部品はある**: `msdial/peak_verification.parse_formula` → `monoisotopic_mass` → `msdial/adducts.parse_adduct` →
  `mz_from_neutral`。付加体の表記は `[2M+H]+` `[M+2H]2+` `[M+FA-H]-` `[M+Na-2H]-` `[M+H-H2O]+` を読める。
- **SME 候補は候補ごとの InChIKey を持たない**: `ds.feature_candidates[fid]` の各候補は `database_identifier`
  （`MspDB:<名前>` 形）・`inchi`・`smiles`・`adduct`・`charge`・`exp_mass_to_charge`・`theoretical_mass_to_charge`・
  `confidence_measure`・`confidence_value`・`spectra_ref` を持つが、InChIKey は rank1 から解決した特徴単位の値
  （`feature_metadata[fid]["inchikey"]`）しか無い。解決（`dataset_state._resolve_inchikey`）は読み込み時にだけ走り、
  研究室 MSP の表（`msp_table`）は保持されない。
- **SML 由来の注釈は別スロット**: `ds.feature_annotations[fid]`（`inchikey` を持つ）。docstring のとおり「証拠ではない」
  （MS1 照合のみ。Text DB 運用では SML が唯一の同定情報源）。`ambiguous` のものがある。
- **差次的統計は、session で計算したものに限り全特徴ぶん session にある**: `dataset_differential` の結果は
  `ds.results[result_id]` に登録され、`results` に前処理済み行列の全特徴（`pp_feature_names`）の
  `{feature, mean_a, mean_b, log2fc, t, p, q}` を持つ。BH（`differential.add_fdr` → `bh_fdr`）は p が有限の全特徴に掛かる。
  前処理（`min_detection_rate` など）で落ちた特徴は `results` に居ない。
  書き出し（`dataset_export.py`）は InChIKey の無い行を捨て、加えてキュレーション（`apply_curation` 既定 true）の
  `wrong` を同定なし扱い・`assign` で同定を置き換え・`redundant` を除外する。
- **pipeline の差次的結果は session に戻らない**: `pipeline/service.py` の `_handle_differential` は結果を run 配下に
  `differential:<cid>` として永続化するが、`dataset_load(pipeline_path=...)`（`pipeline/session_handoff.load_run_dataset`）が
  `ds` に戻すのは試料対応表・binding・解析行列だけ。戻したとしても来歴の `dataset_id` が新しい DatasetState と違うので
  `is_current` は False になる。MTBLS6039 はこの経路（pipeline の `differential_*.tsv`）を使った。§11 で引き継ぐ。
- **v2 の `dataset_statistic` の結果**も `ds.results["stat_<id>"]` に入る（`kind="welch"`、`features` 形、来歴なし）。
- **SML の曖昧な注釈**: `_build_feature_annotations` は 1 特徴に複数の SML 注釈があると個々の InChIKey を捨て、
  `{ambiguous, sml_ids, name: None}` だけを残す。
- `_resolve_inchikey` は引数の `msp_counts` を加算する（候補ごとに呼ぶなら別の dict を渡さないと特徴単位の計数が変わる）。
- **古い結果の判定はある**: `result_state.is_current(ds, result)`（データセット・前処理・群の指紋）。
- **検出マスクがある**: `ds.detected_mask`（特徴 × 試料の bool。隣接 `.arf` と数値で対応が取れたときだけ。無ければ None）。
- **検索の軸の前例**: ARF 経路の `eic_search`（`spot_ids`・`min_mz`/`max_mz`・`min_rt`/`max_rt` の AND、TSV）。
- `dataset_fingerprint` は SME 候補を含まない。候補に項目を足しても、既存の結果は古くならない。

## 3. 目標・非目標・反目標

**ユーザー判断（2026-10-19）**: 判定の方針は `feature_bindings` に揃える。**強度最大では選ばず、複数件なら `needs_input` を返す**。
同定が有って矛盾する候補は落とす。reanalysis-study の手順書 v5（工程 7-1 第 4 段「両群平均強度の和が最大」）の側を v6 で改訂する
（別リポ・別作業）。

目標:
- 窓検索: 特徴 ID・m/z 範囲・RT 範囲の AND で特徴を並べ、注釈と差次的統計を付ける（注釈の無い特徴も）。
- 目標照合: 分子式（または期待 m/z）と InChIKey と規則から、候補の全件・落ちた理由・判定（`resolved` / `needs_input`）を返す。
- 1 回の呼び出しで複数の目標を処理できる（MTBLS6039 は 63 主張 × 6 比較）。
- 規則を正規化して戻り値に写し、その hash を返す（判断ログに「この規則で呼んだ」とだけ書けるように）。

非目標（この版ではやらない）:
- ARF 経路（GUI 解析済みの `.arf`）への対応。純ロジックは入力形式に依存しない形にしておき、ARF の受け口は後で足す。
  GUI の利用者も、MS-DIAL GUI から mzTab-M を書き出して `dataset_load` すれば使える。
- 選択の永続化・手動 override の受け付け（このツールは読み取り専用）。実行役が選んだ結果は reanalysis-study 側が記録する。
- 同定の裏取り（MS/MS・共溶出・同位体）。それは `dcl_find_msms` `library_match_feature` `verify_peak_annotation` の仕事で、
  このツールは「誰に何を確かめさせるか」を返すところまで。
- 差次的エクスポートの列・メタ行の変更（massbank-context との契約。`CONTRACT_VERSION` は上げない）。
- キュレーション（`curation_submit` のフラグ）の反映。この版は SME / SML の生の同定だけを見る。フラグがある
  データセットでは、書き出しと検索で同定が食い違いうるので、`notes` に「キュレーションを反映していない」と書く。
  反映するなら §5 の証拠の水準に `curation` を足す形になる（次の版）。
- pipeline v2 の `bind_features` の挙動の変更。§2 の 2 つの癖（rank1 InChIKey の混入、`formula` だけの目標）は
  `docs/task.md` に残課題として足し、ここでは直さない。

反目標（技術的に動いても失敗とみなすもの）:
- 候補が複数あるのに、強度・並び順・m/z の近さ・RT の近さで 1 件を選んで `resolved` と返す。
- 群の差や論文の値との相関で候補を選び直す（相関は選んだ後の確認であって、選ぶ根拠ではない）。
- 前処理で落ちた特徴・検定できなかった特徴の統計欄を空欄だけで返し、「差が無かった」と読める形にする。
- 古い差次的結果（前処理・群の訂正前）の数字を、今の結果として返す。
- 戻り値が全特徴ぶん膨らむ（`bind_features` の全候補をそのまま返す）。

## 4. ツールの外形

名前: `dataset_feature_search`（`metabolomix/tools/` に置く。読み取り専用 `readOnlyHint=True`、`structured_output=False`）。

```
dataset_feature_search(
    feature_ids: list[str] | None = None,
    min_mz: float | None = None, max_mz: float | None = None,
    min_rt: float | None = None, max_rt: float | None = None,
    targets: list[dict] | None = None,
    result_ids: list[str] | None = None,
    max_results: int = 20,
) -> str
```

- **窓検索**: `targets` を渡さないとき。`feature_ids`（SMF_ID の完全一致）・m/z・RT（分）の AND。境界を含む。
  何も指定しなければエラー（全件 11 万行を返さない）。並びは m/z 昇順。`max_results` で打ち切り、打ち切った件数を書く。
- **目標照合**: `targets` を渡したとき。窓の引数と同時に渡したらエラー（どちらの意味か曖昧になる）。目標は 1 回 50 件まで。
- `result_ids`: 統計を付ける差次的結果。省略時は `ds.results` のうち今も有効な（`is_current`）2 群比較
  （`kind == "two_group"`）の**全部**。v2 の `dataset_statistic` の結果（`kind="welch"`）は対象外。
  名指しした ID が無い・古いときは `missing_state`（`dataset_differential_result`、`required_tools=["dataset_differential"]`）。
  省略時に古いものがあれば、使わずに ID を注記する。差次的結果が 1 件も無いときは統計なしで返し、そう注記する（エラーにしない）。
- `dataset` が無ければ `missing_state`（`dataset`、`["dataset_load"]`）。

### 4.1 目標（`targets` の 1 要素）

| キー | 必須 | 意味 |
|---|---|---|
| `target_id` | ○ | 呼び出し側の識別子（主張 ID など）。戻り値の各行に付く |
| `inchikey` | ○ | 照合は**先頭ブロック**（14 文字）で行う。全体を渡してよい（先頭だけ使う） |
| `formula` | どちらか | 中性分子の分子式。付加体ごとの理論 m/z を計算する |
| `expected_mz` + `adduct` | どちらか | 分子式が無いとき。付加体は 1 つに固定される（層は使わない） |
| `name` | | 表示用。照合には使わない |
| `ppm` | | 質量の許容幅（既定 10） |
| `rt_min` | | 期待 RT（分）。無ければ RT の条件を掛けない |
| `rt_tol_min` | | RT の許容幅（既定 0.3） |
| `adduct_tiers` | | 付加体の層（下記）。省略時は極性の既定 |
| `polarity` | | `positive` / `negative`。読み込んだデータセットの極性と違えば、その目標は `status: not_applicable` で照合しない |

`polarity` を持たせるのは、session が持てる DatasetState が 1 つだけで、陽陰が別データセットになるため。呼び出し側は
陽陰の目標をまとめて渡し、データセットを切り替えながら同じ呼び出しを繰り返せばよい。

InChIKey を必須にするのは、名前や分子式だけでは同定の一致・矛盾を判定できないため（`feature_bindings` の「名前だけは不可」と同じ考え）。

**付加体の層の既定**（手順書 v5 工程 7-1 の第 1 段と同じ集合。極性は読み込んだ mzTab-M から決める——§7）:
- 陽: 第 1 層 `[M+H]+` `[M]+` ／ 第 2 層 `[M+Na]+` `[M+NH4]+` `[M+H-H2O]+` `[M+H-NH3]+` `[2M+H]+` `[M+2H]2+`
- 陰: 第 1 層 `[M-H]-` ／ 第 2 層 `[M-H2O-H]-` `[M+Cl]-` `[M+FA-H]-` `[2M-H]-` `[M+Na-2H]-`

極性が決まらないときに `adduct_tiers` も `expected_mz`+`adduct` も無ければ、その目標はエラー（`POLARITY_UNKNOWN`）。

### 4.2 戻り値（`json_payload`）

窓検索:
- `features`: TSV（列名 1 回）。列 `feature_id mz rt evidence name inchikey adduct n_sme`。
  `evidence` は `sme`（SME 候補あり）/ `sml`（SML 注釈のみ）/ `none`。`name` `inchikey` `adduct` は SME rank1、無ければ SML。
- `stats`・`results`・`notes`（下記と同じ）。

目標照合:
- `targets`: 目標ごとに `{target_id, status, reason, selected_feature_id, selected_adduct, selected_evidence,
  adduct_tier_used, counts, theoretical_mz}`。
  - `status` は `resolved` / `needs_input` / `not_applicable`（極性違い）。`reason` は `no_qualified_candidate` /
    `multiple_qualified_candidates` / `polarity_mismatch`（`resolved` なら null）。
  - `counts` は段ごとの件数 `{mass, rt, qualified, by_evidence: {sme_match, sml_match, unannotated}}`（0 件の理由を人が読めるように）。
  - `theoretical_mz` は付加体ごとの理論 m/z（丸めて）。
- `candidates`: TSV。列 `target_id feature_id adduct adduct_tier theoretical_mz mz ppm rt delta_rt evidence sme_rank
  annotated_name annotated_inchikey qualified reasons`。
  質量窓（全層のいずれかの付加体）に入った特徴だけを載せる（RT で落ちたものも理由つきで載せる）。
  目標ごとに `max_results` 件まで、超えた件数を `targets[*].truncated` に書く。並びは（qualified 先、付加体の層、|ppm|）。
- `stats`: TSV。`candidates`（目標照合）または `features`（窓検索）に出た特徴だけ。列
  `feature_id result_id status log2fc p q mean_a mean_b detected_a detected_b`。
  - `status`: `tested` / `not_in_matrix`（前処理で落ちた）/ `not_testable`（p が NaN）。後の 2 つは数値欄を空にする。
  - `detected_a` / `detected_b`: 群ごとの検出例数（`detected_mask` があるときだけ。無ければ列を空にし、`notes` に理由）。
    gap-fill の値と実測を区別するための列（教訓 §3.7: 論文の AUC=1 の中身が埋め値だった）。
- `results`: 統計を付けた結果ごとに `{result_id, group_a, group_b, n_a, n_b, log2fc_sign, bh_population, caveats}`。
  `bh_population = {n_features_total, n_in_matrix, n_tested}`（BH は `n_tested` 件に掛かっている）。
  `caveats` は差次的結果の caveat をそのまま（交絡の警告をここで落とさない）。
- `rule`: 正規化した照合規則（既定値を埋めた目標・**解決後の**付加体の層・許容幅・極性）と `rule_hash`（`canonical_hash`）。
  result_id は計算のたびに発行される UUID なので **`rule_hash` に入れない**。統計の出所は別に
  `stats_source = {result_ids, dataset_hash}`（`dataset_hash` は `dataset_fingerprint`）として返す。
  同じ目標・同じ極性なら、別の session・別の実行役でも `rule_hash` は一致する。
- `notes`: 打ち切り・古い結果の除外・検出マスクが無い・極性の出所・キュレーションのフラグがあるのに反映していないこと など。

数値は丸める（m/z 4 桁、ppm 2 桁、RT 3 桁、log2FC 4 桁、p・q・平均は有効 4 桁）。

## 5. 判定の規則（目標照合）

各目標について、付加体の層を上から順に見る。

1. **質量**: 層の各付加体について理論 m/z（分子式から、または `expected_mz`）を出し、特徴の m/z（`feature_metadata[fid]["mz"]`）が
   ±`ppm` 以内の特徴を候補にする。候補の単位は（特徴, 付加体）。
2. **RT**: `rt_min` があれば |RT − rt_min| ≤ `rt_tol_min`。外れたら `rt_outside_tolerance`。
3. **同定**（候補ごとに 1 つの `evidence` を決める）:
   - SME 候補のどれか（rank を問わない）が、InChIKey 先頭ブロックが一致し、付加体（正規化して比較）と電荷も一致 → `sme_match`
     （一致した rank を `sme_rank` に書く）。
   - SME 候補があり、識別子を持つ候補がどれも一致しない → **落とす**（`identifier_mismatch` / `adduct_mismatch` / `charge_mismatch`）。
     ユーザー判断「同定が有って矛盾する候補は落とす」。
   - 「識別子を持つ候補」とは、**InChIKey を解決できた** SME 候補のこと（§7 で読み込み時に候補ごとに解決する）。
     `MspDB:<名前>` だけで MSP の表も RDKit も無く InChIKey が取れない候補は、一致とも矛盾とも言えない
     （`feature_bindings` の「判定材料を持たない候補は評価できない」と同じ）。こうした候補しか無い特徴は SME が無いものとして
     次へ進み、`reasons` に `sme_identity_unresolved` を書く（名前上は別物質の可能性があることを人が見られるように）。
   - SME 候補が無い（または識別子を持たない）で、SML 注釈の InChIKey 先頭ブロックが一致 → `sml_match`。
     SML が曖昧（1 特徴に複数の注釈）なときは、そのどれかが一致すれば `sml_match` とし、`reasons` に `sml_ambiguous` を書く
     （§7 で、曖昧な注釈の InChIKey の一覧を読み込み時に残す）。
   - それ以外（注釈なし、または SML 注釈が別の物質）→ `unannotated`。SML が別物質でも落とさない（下の判断を参照）。
     `unannotated` は **`rt_min` があるときだけ** qualified にする（質量だけで同定とは言わない。`feature_bindings` の
     `mass_rt` が RT を必須にしているのと同じ）。`rt_min` が無ければ `mass_only_not_identified` で落とす。
4. **判定**: 「qualified の候補が 1 件以上ある最初の付加体の層」を使い、その層の中で証拠の水準
   `sme_match` → `sml_match` → `unannotated` の順に、qualified が 1 件以上ある最初の水準を見る。
   - その水準の qualified がちょうど 1 件 → `resolved`。
   - 2 件以上 → `needs_input`（`multiple_qualified_candidates`）。強度・ppm・RT の近さでは選ばない。
   - どの層にも qualified が無い → `needs_input`（`no_qualified_candidate`）。
   - 使わなかった層・水準の候補も `candidates` に残す（qualified のまま。`reasons` に `lower_adduct_tier` / `lower_evidence_tier` を情報として書く）。

**SML の別物質注釈を矛盾として扱わない理由**: SML は MS1 の照合だけなので、同じ分子式の異性体はどれも同じ m/z に当たる。
SML が別の異性体を付けていても、それは「この物質ではない」証拠にならない。SME（MS/MS を見た同定）だけを矛盾の根拠にする。

## 6. 統計の付け方

- `result_ids` の各結果について、`results` を `feature` で引く。居なければ `not_in_matrix`、p が有限でなければ `not_testable`。
- `q` はその結果に付いた値をそのまま使う（再計算しない）。BH の母集団は `n_tested`。
- 検出例数: `detected_mask` の行（`ds.feature_ids` の位置）と列（結果の `samples_a` / `samples_b` を `ds.sample_names` で位置に直す）で数える。
  列の対応が取れないとき（マスクの列数が `ds.sample_names` と合わない等）は数えずに注記する
  （既存 `dataset_analysis._restrict_mask_to_samples` と同じ考え方）。

## 7. 既存コードへの変更

- **新規 `metabolomix/analysis/feature_search.py`**（純ロジック。session も MCP も import しない）:
  窓検索、目標の正規化と検証、理論 m/z、§5 の判定、§6 の統計の結合。入力は DatasetState 相当の属性
  （`feature_ids` `feature_metadata` `feature_candidates` `feature_annotations` `results` `detected_mask` `sample_names`）だけを読む。
- **`metabolomix/analysis/feature_bindings.py`**: 共有したい部品（`_normalize_adduct` `_ppm_delta` `resolve_candidates` の考え方）を
  import して使う。公開名にするか、`feature_search` 側で import するだけにするかは実装者に任せる。**`bind_features` の挙動は変えない。**
- **`metabolomix/mztab/dataset_state.py`**:
  - 読み込み時に、SME 候補ごとに `inchikey` と `inchikey_source` を解決して候補 dict に足す（`_resolve_inchikey` を候補ごとに呼ぶ）。
    `inchikey_coverage` の計数（特徴単位）は変えない——`_resolve_inchikey` には特徴単位とは**別の** `msp_counts` の dict を渡す。
    RDKit 経路の費用は SME 行数に比例する（実測例 ST003841 で SME 1,165 行）。
  - 曖昧な SML 注釈（`ambiguous`）に、各注釈の InChIKey の一覧（例 `inchikeys`）を足す。既存の `inchikey` を入れないことは変えない
    （書き出しは今どおり曖昧な注釈を同定なしとして扱う）。
- **極性**: 読み込んだ mzTab-M について `console/output_collector.read_adduct_polarity` の多数決を使う。ジョブ由来なら
  ジョブの宣言極性があればそれを優先し、両者が食い違えば `notes` に書く。どこから取ったかを `rule.polarity_source` に書く。
- **新規ツール**: `metabolomix/tools/` の既存モジュール（`dataset_analysis_tools.py`）に足すか、新モジュールにするかは実装者に任せる。
  `__all__` に載せ、`server.py` の再エクスポートに乗せる。

## 8. 文書・テスト・登録

- `USAGE.md` のツール表と件数、CLAUDE.md 冒頭の規模表記（74 → 75）、`tests/test_server_registration.py` の件数。
- `docs/workflow/`（`tests/test_workflow_docs.py` が実登録と突き合わせる）に呼び出しの流れを足す。
- `docs/output_format/mztab.md` に戻り値の欄の意味（`evidence` の 3 値、`status` の 3 値、`bh_population`、`detected_*`、
  `needs_input` の読み方）を足す。
- vault の流れ図（`ms-data-parser-flow.md`）に、DatasetState 経路の解析枝として足す（CLAUDE.md「作業の記録と Git」）。
- テスト（fixture はテスト自身が作る）:
  - 窓検索: AND・境界・打ち切り・空指定のエラー・注釈の無い特徴に統計が付く。
  - 目標照合: 層の順（第 1 層に qualified があれば第 2 層を使わない）、証拠の水準の順、複数件で `needs_input`、0 件で `needs_input`、
    SME の矛盾で落ちる、SML の別物質では落ちない、`rt_min` 無しで `unannotated` が qualified にならない、
    rank2 の SME 一致を `sme_rank=2` で拾い、rank1 の InChIKey を rank2 に混ぜない。
  - 統計: `not_in_matrix` / `not_testable`、古い結果の除外と名指し時の `missing_state`、`bh_population`、検出マスクの有無。
  - 規則の hash が同じ入力で同じ、既定値を明示しても同じ。
  - 強度が最大の候補が 1 件あっても、qualified が 2 件なら `needs_input` のまま（反目標の回帰防止）。
  - `bind_features` の既存テストが変わらず通る。

## 9. 実装者に任せること

- 内部の関数分割・モジュール名・ツールを置くファイル。
- `feature_bindings` の部品を公開名にするかどうか。
- TSV の組み立て（`format_spots_as_table` を流用するか）。
- 目標の検証エラーの封筒の形（既存の `mztab_error` / `DomainError` の流儀に合わせる）。

## 10. 設計者（Claude）の判断として置いたもの

ユーザーは質問を省略したので、次は設計者の判断で決めた。読みながら確認してほしい。

1. 窓検索と目標照合を 1 ツールにまとめ、`targets` の有無で切り替える（両方の同時指定はエラー）。
2. 付加体の層・証拠の水準の順による段階的な絞り込みは残す（v5 の順序を、強度による決着だけ外して保つ）。
   どちらも規則による決定的な順序で、「並び順に従った選択」ではないと判断した。
3. SME 候補は rank を問わず一致を認める（`feature_bindings._check_identity` と同じ）。
4. SML の別物質注釈は矛盾にしない（§5 の理由）。
5. `unannotated` は RT 条件があるときだけ qualified。
6. InChIKey を目標の必須にする（分子式だけの目標は受け付けない）。
7. `result_ids` の省略時は、有効な 2 群比較を全部付ける。
8. 群ごとの検出例数を統計に含める（教訓の埋め値の件から）。
9. 目標は 1 回 50 件まで、候補は目標ごとに `max_results`（既定 20）件まで。
10. InChIKey を解決できない SME 候補は「評価できない」とし、矛盾として落とさない（`sme_identity_unresolved` を書く）。
    落とすと、MSP の表が無い読み込み（mzTab-M の直読み）で SME 付きの特徴がほぼ全部落ちる。
11. 曖昧な SML 注釈は、どれかが一致すれば `sml_match`（`sml_ambiguous` を書く）。読み込み時に InChIKey の一覧を残す。
12. 目標に任意の `polarity` を持たせ、違えば `not_applicable`。
13. `rule_hash` は照合規則だけから作り、統計の出所（result_ids・dataset_hash）は `stats_source` に分ける。
14. 対象は `kind == "two_group"` の結果だけ（v2 の `dataset_statistic` は除く）。
15. キュレーションはこの版では反映せず、フラグがあれば注記する。

## 11. pipeline の差次的結果を session へ戻す（ユーザー判断 2026-10-19: 案 A）

session で `dataset_differential` を回した結果しか `ds.results` に無い（§2）ので、pipeline の run を読み込んだだけでは
このツールは「統計なし」を返す。比べた案:

- **A（採用）**: run に永続化された差次的結果を、`dataset_load(pipeline_path=...)` が `ds.results` へ戻す。
  データの指紋（`dataset_fingerprint`）と群の指紋が一致するときだけ有効とみなす。
- B（不採用）: 対象を session に限り、利用者が比較の数だけ session で回し直す。レシピの写し間違いで黙ってずれ、
  陽陰を切り替えるたびにやり直しになる（MTBLS6039 なら 6 比較 × 2 極性）。

A を採る理由: 再解析の記録で「pipeline が書き出した TSV の数字」と「検索で見た数字」が同じ計算から出ることが保証される。

### 11.1 確認した事実

- pipeline の `differential` 工程（`pipeline/service.py` の `_handle_differential` → `dataset_service.run_comparison` →
  `compare_dataset`）は、結果の全量（`results` を含む）を `report.persist_result` で JSON に書き、`record["results"]` に
  `output_name = "differential:<cid>"` の ref を追記する（追記専用）。`export` 工程は、同じ pass に結果が無いとき再計算して
  **もう一度** `differential:<cid>` として永続化し、TSV のメタ行 `# result_id` はその新しい ID を指す。
- 来歴（`result_state.new_provenance`）は `dataset_id`・`input_fingerprint`・`parent_ids`（前処理の ID）・`effective_parameters`・
  `group_fingerprint` を持つが、**データの指紋（`dataset_fingerprint`）は持たない**。`input_fingerprint` は前処理済み行列と設定の
  指紋（`dataset_service._derived_fingerprint`）。
- `session_handoff.load_run_dataset` は `mztab_loading.load_dataset_state(job_path=...)` で DatasetState を作り直す
  （`dataset_id` は新しく振られ、`preprocess_id` は None）。試料対応表は `ds.sample_metadata_rows` に戻すので、
  `group_fingerprint(ds)` は run と同じ規則で計算できる。
- `result_state.is_current` は `dataset_id` の一致と、親の前処理 ID が `ds.preprocess_id` であることを要求する。
  戻した結果はこのままでは必ず「古い」になる。
- 成果物の読み出しは `store.read_result_data(pipeline_root, results, output_name)`（hash を照合してから読む。完全一致の
  `output_name` で引く）。

### 11.2 変更

1. **来歴にデータの指紋を刻む**: `dataset_service.compare_dataset` が、`group_fingerprint` と同じ場所で
   `provenance["dataset_hash"] = result_state.dataset_fingerprint(ds)` を足す。session・pipeline の両経路に効く。
   `input_fingerprint` の中身は変えない（既存の照合を動かさない）。
2. **引き継ぎ**: `session_handoff.load_run_dataset` が、`record["results"]` の `output_name` が `differential:` で始まる ref を
   比較 ID（`<cid>`）ごとに集め、**各比較の最後に永続化されたもの**を読む（`export` が再計算した場合、TSV が名指しするのは
   最後のもの）。次の 2 つが一致したものだけを `ds.results[result_id]` に登録する:
   - `provenance.dataset_hash == dataset_fingerprint(ds)`
   - `provenance.group_fingerprint == group_fingerprint(ds)`（試料対応表を戻した後に計算する）

   result_id は**付け替えない**（pipeline の TSV のメタ行と同じ ID で引けるように）。来歴に
   `provenance["handoff"] = {"pipeline_root", "comparison_id"}` を足す。`ds.last_differential` は設定しない
   （「直近の結果」は session で計算したものに限る。戻した結果は result_id で名指しして使う）。
3. **戻せなかったもの**: 理由つきで `dataset_load` の戻り値に並べる。理由は
   `legacy_result_without_dataset_hash`（この変更より前の run）/ `dataset_mismatch` / `group_mismatch`。
   旧 run は差次的解析を計算し直せば指紋つきの結果ができる。`pipeline_resume` で完了済みの差次的工程を再計算させる手順
   （比較の再送など）が使えるかは**未確認**——実装の最初に確かめ、使えるならその手順を、使えなければ session での
   `dataset_differential` を戻り値の注記で案内する。
4. **有効性の判定**: `result_state.is_current` に分岐を足す。来歴に `handoff` がある結果は、`dataset_id` と親の前処理 ID を
   見ずに、`dataset_hash == dataset_fingerprint(ds)` と群の指紋の一致で判定する（前処理は run が行ったもので、そのレシピは
   `effective_parameters` と run の記録にある。session で後から前処理をやり直しても、戻した結果の数字は変わらない）。
   群を訂正すれば（`dataset_set_sample_metadata`）群の指紋が変わり、戻した結果は古いと判定される。
   `handoff` の無い結果の判定は今のまま。
5. **`dataset_load` の戻り値**: 戻した結果を `比較 ID・result_id・group_a・group_b` の 1 行ずつ（TSV）、戻せなかったものを理由つきで。

### 11.3 副作用（意図したもの）

- 戻した結果は `dataset_export_differential(result_id=...)` や図のツールからも「有効な結果」として使える。
  書き出せば pipeline の TSV と同じ数字になる。
- メモリ: 比較ごとに全特徴ぶんの結果を持つ（陰 11 万特徴 × 6 比較）。pipeline の worker が持つのと同じ量。

### 11.4 テスト

- 来歴に `dataset_hash` が刻まれる（session・pipeline の両経路）。
- 引き継ぎ: 指紋が一致する結果が result_id そのままで `ds.results` に入る。比較ごとに最後の永続化を選ぶ。
  旧 run（`dataset_hash` 無し）・データ違い・群違いは戻さず理由を返す。hash 不一致の成果物は既存どおり `RESULT_INTEGRITY_MISMATCH`。
- `is_current`: 戻した結果は session の前処理の有無に関係なく有効。群を訂正すると古くなる。`handoff` の無い結果の判定は変わらない。
- `pipeline/{engine,worker,service,recovery,store}.py` がグローバル session を import しない規約（AST テスト）を破らない。
- 一連の流れ: run を読み込む → `dataset_feature_search` が pipeline の比較の統計を付ける（再計算しない）。

### 11.5 文書

- `USAGE.md` の `dataset_load` の説明と、`docs/workflow/` の `dataset_load(pipeline_path=...)` の流れ。
- vault の流れ図: run を読み込んだあと、差次的結果が session で使えるようになる（前提状態の連鎖が変わる）。
