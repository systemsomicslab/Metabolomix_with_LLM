# 再解析 v8: Console の試し（Task 10）の記録

2026-10-30〜31。MTBLS6965 の R17_1-3・P17_1-3（陽陰 6 ファイルずつ）を `library_mode=msp_then_lbm` で回した。
MS-DIAL Console は MsdialCoreTestApp（Debug net48、45a531c4d）。作業場所は `datasets/reanalysis/_v8_trial/`。

## 1 回目（修正前、31c4b6f）

- 未署名の `MsdialGcMsApi.dll` が Smart App Control に拒否され、2 本とも failed になった（CodeIntegrity 3033/3077）。
  PC の再起動の後、neg は通った。pos は再起動の直後にもう一度拒否され、次の再開で通った。
- 2 本とも completed になった。2 本を同時に回し、約 33 分かかった。
  MSDIALCUI の WorkingSet は、1 本の最大が 25.8 GB、2 本の合計の最大が 27.1 GB だった。
- 結果は期待と違った。研究室 MSP の候補が足切りを越えるのに LBM が代表になる特徴が、26 件あった（pos 12、neg 14）。
  例は TG 54:3 と Glyceryl trioleate。
- 原因: Console は `Weighted/Simple/Reverse dot product cutoff for MSP-based annotation` の値を**二乗値**として読む。
  そのため、メソッドの 0.6/0.6/0.8 は実効で 0.775/0.775/0.894 になる。
  LBM には足切りのキーが無いので既定の 0.6/0.6/0.8 が使われ、MSP だけが厳しくなっていた。
  MSP だけの 1 本（R17_3 pos）でこの規則を確かめた。.mdpeak の `MS/MS matched` は、13 ファイルで食い違いが 0 だった。
- 修正 de70012: `msp_then_lbm` は、MSP と LBM の dot product の足切りを `Square root of ... cutoff for {MSP,LBM}-based annotation` のキーで同じ値に書き、二乗のキーは空にする。

## 回し直し（de70012）

| 確かめること | pos | neg |
|---|---|---|
| 完走 | completed | completed |
| アライメントの基準（`alignment_reference`） | file_id 2 = R17_3（P17-R-33） | file_id 2 = R17_3（P17-R-33） |
| effective-method | 5 キーと平方根の足切り 6 キー（0.6/0.6/0.8）。二乗のキー 3 つは空 | 同じ |
| mzTab の SME（MSP / LBM） | 835 / 528（1 回目 804 / 559） | 602 / 571（1 回目 563 / 610） |
| 代表が MSP で matched | 553（1 回目 448） | 356（1 回目 251） |
| 代表が LBM で matched | 261（1 回目 292） | 223（1 回目 262） |
| LBM が代表なのに、MSP の候補が足切りを越える特徴 | 0（1 回目 12） | 1（1 回目 14） |

effective-method の 5 キーは、`Msp file path` が空、`MSP annotator settings file path`、`Lbm file path`、LBM の RT のキー 2 つが False。

- **両方の候補があるとき、rank 1 が MSP になる例**: pos の Alignment ID 1284（m/z 902.822、RT 11.99、[M+NH4]+）。
  1 回目の代表は LBM の `TG 54:3|TG 18:0_18:1_18:2` だった。
  回し直しでは研究室 MSP の `Glyceryl trioleate`（matched、wdp 0.666・sdp 0.712・rdp 0.863）になった。
- 残る 1 件は neg の 1108（`PC 35:1|PC 17:0_18:1`）。こちらで採点し直すと、MSP の PC 35:1 は wdp 0.748 で足切りを越える。
  ただし、この採点は .mdalign の代表スペクトルに対するもので、MS-DIAL の判定（代表ピークのスペクトル）とは少しずれ得る。
  1 回目に採点が一致したのは 593 件中 592 件だった。この 1 件はそのずれの範囲とみなし、追っていない。
- メモリ: 回し直しの最大は記録できていない（見張りのシェルがメモリ不足で止められた）。設定はキーの値しか変わらないので、1 回目の値を使う。

## 残り

- 過去の accession（MTBLS3233、MTBLS6039、MTBLS6965、ST002901、ST003841）の effective-method にも、同じ二乗のキー（0.6/0.6/0.8）がある。
  そのため、研究室 MSP の同定は意図より厳しかった可能性がある。どう扱うかは別に決める。

## パスウェイの受け渡しで落ちる行（最終レビューの確認）

`msp_then_lbm` の試し（`_v8_trial/pipeline_r2` の pos / neg。MTBLS6965 の R17 と P17 の 2 群比較）の mzTab を
`dataset_handle_open` で開き、`dataset_preprocess` の後に `pathway_project` が使う `handoff.build_from_dataset` を回した
（MSLipidMapper・R は呼ばない）。MSP 由来かは `database_identifier` が `MspDB_` で始まるかで分けた。

| 数（前処理後の特徴） | pos | neg |
|---|---|---|
| 特徴 | 1666 | 1619 |
| 残った行（同名の整理の後。`n_kept`） | 251 | 388 |
| `dropped_unannotated` | 461 | 350 |
| `dropped_non_msdial_ontology` | 864 | 758 |
| うち MSP 由来（`MspDB_`） | 864 | 758 |
| うち MSP 由来で MS/MS 照合の旗が真 | 512 | 349 |
| `dropped_unparsable_name` | 44 | 109 |
| うち MSP 由来 | 0 | 0 |
| `dropped_duplicate_names` | 46 | 14 |
| 残った行のうち MSP 由来 | 0 | 0 |

- 落ちた MSP 由来の Ontology は ClassyFire 系（Phosphatidylcholines、Phosphatidylethanolamines、Triacylglycerols、Acyl carnitines、
  Long-chain ceramides など）で、MS-DIAL の脂質クラスと一致しない。
- 残った行はすべて LBM 側の同定。MSP で代表になった同定は 1 行も残らない。クラスごとの合計は LBM で同定した分子種が中心になる偏りがある。
- 扱い: 名前→クラスの対応表は作らない（判断）。`docs/output_format/pathway.md` と `docs/workflow/reanalysis.md` に偏りを書き、
  `pathway_project` の `handoff`（`n_kept` を追加）と、落とした行があるときの `notes.handoff` で件数を読めるようにした。
- 追記（2026-11-02）: この「対応表は作らない」は、その後の変更で置き換わった。MSP の行は handoff が名前か ClassyFire のクラスと組成式から直して入れる。変更後の数は次の節。

## MSP の脂質の handoff（2026-11-02）

ブランチ `feat/msp-pathway-handoff`（handoff が研究室 MSP の行を名前か ClassyFire のクラスと組成式から MS-DIAL の名前とクラスに直して入れる。
設計は `2026-11-02-msp-lipid-pathway-handoff-design.md`）の後で、上の節と同じ手順を回した。同じ mzTab（`_v8_trial/pipeline_r2` の pos / neg）を
`dataset_handle_open(pipeline_path=…)` で開き、`dataset_preprocess`（既定のレシピ。pipeline の `conservative-v1` が解決した
レシピも `impute=half_min` だけで同じ）の後、`pathway_project` と同じ引数で `handoff.build_from_dataset` を回した（群は R17 と P17 の各 3 試料。
MSLipidMapper・R は呼ばない）。MSP 由来かは `database_identifier` が `MspDB_` で始まるかで分けた（残った行は、同名の整理で勝った行の注釈で数えた）。

| 数（前処理後の特徴） | 前 pos | 前 neg | 後 pos | 後 neg |
|---|---|---|---|---|
| 特徴 | 1666 | 1619 | 1666 | 1619 |
| 残った行（同名の整理の後。`n_kept`） | 251 | 388 | 448 | 672 |
| 残った行の絞る前の数（`kept_by_route` の合計。前は `n_kept` + `dropped_duplicate_names` で出した） | 297 | 402 | 668 | 814 |
| `kept_by_route.as_is` | – | – | 297 | 402 |
| `kept_by_route.name` | – | – | 1 | 227 |
| `kept_by_route.class_formula` | – | – | 370 | 185 |
| `dropped_unannotated` | 461 | 350 | 461 | 350 |
| `dropped_non_msdial_ontology` | 864 | 758 | （キーなし） | （キーなし） |
| `dropped_unparsable_name` | 44 | 109 | 44 | 109 |
| `dropped_isotope_label` | – | – | 0 | 1 |
| `dropped_ambiguous_class` | – | – | 0 | 0 |
| `dropped_ambiguous_ontology` | – | – | 62 | 27 |
| `dropped_unmapped_ontology` | – | – | 414 | 287 |
| `dropped_class_name_conflict` | – | – | 2 | 3 |
| `dropped_no_formula` | – | – | 0 | 0 |
| `dropped_formula_mismatch` | – | – | 15 | 28 |
| `dropped_duplicate_names` | 46 | 14 | 220 | 142 |
| 残った行（`n_kept`）のうち MSP 由来 | 0 | 0 | 199 | 284 |
| うち MS/MS 照合の旗（`spectrum_matched`）が真 | – | – | 140 | 175 |

- 「前」は上の節の数（`dropped_non_msdial_ontology` は 864 / 758 で、すべて MSP 由来）。「後」の `dropped_non_msdial_ontology` は新しい report に無い
  （`isotope_label`・`ambiguous_class`・`ambiguous_ontology`・`unmapped_ontology`・`class_name_conflict`・`formula_mismatch` に分かれた）。
  「–」は前の report に無かった項目。「後」は最終レビューの修正（`isotope_label`・`class_name_conflict`・経路 `name` の `;O…`）の後に測り直した数。
- 区分は重ならず、全行 = `kept_by_route` の合計 + 落ちた理由の合計（`dropped_duplicate_names` を除く）:
  pos 668 + 461 + 44 + 0 + 0 + 62 + 414 + 2 + 0 + 15 = 1666、neg 814 + 350 + 109 + 1 + 0 + 27 + 287 + 3 + 0 + 28 = 1619。
  `n_kept` = 合計 − `dropped_duplicate_names`（pos 668 − 220 = 448、neg 814 − 142 = 672）。
- 前の 864 / 758 の MSP 行（ClassyFire 系の Ontology）の行き先: 経路 `name` と `class_formula` で pos 371（1 + 370）/ neg 412（227 + 185）が入り、
  落ちたのは pos 493（62 + 414 + 2 + 15）/ neg 346（1 + 27 + 287 + 3 + 28）。経路 `as_is` の 297 / 402 は、前の絞る前の数（251 + 46、388 + 14）と同じ。
- MSP の SML 行は pos 980・neg 811。うち `dropped_unannotated` が pos 116・neg 53（Ontology が空・Unknown・Others）。
- 修正の前に同じ手順で測った数（`kept_by_route.class_formula` 372 / 189、`n_kept` 449 / 674、`dropped_duplicate_names` 221 / 144、
  MSP 由来の `n_kept` 200 / 286）から、行の行き先が変わったのは 6 行（特徴ごとの解決を突き合わせた）。
  `class_name_conflict` の 5 行（pos 2・neg 3）はすべて、名前 `1-(1Z-Octadecenyl)-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-sn-glycero-3-phosphoethanolamine`・
  Ontology `1-(1Z-alkenyl),2-acyl-glycerophosphocholines`・組成式 C45H78NO7P で、前は経路 `class_formula` の `PC P-37:6`。
  `isotope_label` の 1 行は neg の `no MS2: DHA-[21,21,22,22,22-d5]`（`Very long-chain fatty acids`、C22H32O2）で、前は `no MS2: FA 22:6`。
  経路 `name` の `;O…` の修正で行き先が変わった行は無かった。
- `dropped_duplicate_names` を、同名になった行の群の由来（`database_identifier` が `MspDB_` で始まるか）で分けた数
  （群の数 / 落ちた行）: pos は MSP だけ 98 / 169、LBM だけ 38 / 44、混在 3 / 7。neg は MSP だけ 66 / 128、LBM だけ 14 / 14、混在 0 / 0。
  MSP だけの群のうち、元の MSP の化合物名が 2 つ以上ある群は pos 55 群（154 行、落ちた行 99）、neg 21 群（81 行、落ちた行 60）。
  同名になった名前の例（report の `duplicate_name_examples`。名前順の先頭 5 件）は `CL …`・`Cer 18:1;O2/…` など。

`unmapped_ontology` の多い ClassyFire の Ontology 名（上位 5 つ。件数は行数）:

| 順位 | pos | neg |
|---|---|---|
| 1 | Benzoic acid esters（18） | Dihydroxy bile acids, alcohols and derivatives（25） |
| 2 | Polyethylene glycols（14） | Triterpenoids（12） |
| 3 | N-acyl amines（13） | Dipeptides（11） |
| 4 | N-acyl-alpha amino acids（12） | Prostaglandins and related compounds（10） |
| 5 | Triterpenoids（12） | Nitrophenols（9） |

`ambiguous_ontology` の ClassyFire の Ontology 名（pos 62・neg 27 行は、次の 5 種類ですべて。上位 5 つ。件数は行数）:

| 順位 | pos | neg |
|---|---|---|
| 1 | Long-chain ceramides（22） | Phosphosphingolipids（8） |
| 2 | Phosphosphingolipids（14） | Ceramides（6） |
| 3 | Ceramides（12） | Long-chain ceramides（6） |
| 4 | Cholesterols and derivatives（7） | Glycosyl-N-acylsphingosines（6） |
| 5 | Glycosyl-N-acylsphingosines（7） | Cholesterols and derivatives（1） |

- `unmapped_ontology` と `ambiguous_ontology` の行はすべて MSP 由来（pos 414 / 62、neg 287 / 27）。件数が同じ Ontology の並びは、出現順（`Counter.most_common`）。
- 測定のスクリプトは追跡外（`.superpowers/sdd/2026-11-02-msp-lipid-pathway-handoff/measure_handoff.py`）。結果は観測のまま。
