# 研究室 MSP で同定した脂質をパスウェイの handoff に入れる — 設計

日付: 2026-11-02。ブランチ `feat/msp-pathway-handoff`（main 4e7cc04 から）。

## 1. 目的と理由

再解析 v8 の脂質の run は、研究室 MSP を優先し LBM で補って同定する（`library_mode="msp_then_lbm"`）。
同じ証拠の水準なら MSP が代表になるので、脂質の多くが MSP の行として mzTab に出る。ところがパスウェイの
handoff（`metabolomix/pathway/handoff.py`）は、Ontology が MS-DIAL の脂質クラス（`LbmClass`）でない行と、
名前が脂質の略記として読めない行を落とす。MSP の行の Ontology は ClassyFire 系の名前（`Phosphatidylcholines` 等）で、
名前は慣用名・IUPAC 名・`MassSpecGymID…` などなので、ほぼすべて落ちる。

2026-10-30 の試し（MTBLS6965 の 6 試料、`_v8_trial/pipeline_r2`）では、handoff に残ったのは陽 251・陰 388 行で、
すべて LBM の行だった。MSP の行は陽 864・陰 758 行が落ち、そのうち MS/MS で照合したもの（`spectrum_matched`）は
陽 512・陰 349 行だった。パスウェイのノードとクラスの和は、LBM で同定した種に偏っている。

ユーザーの決定（2026-10-28）は「LBM と MSP の同定を区別しない」。handoff だけがこれに反して MSP の行を
落としているので、MSP の行を LBM の行と同じ形（MS-DIAL のクラスと略記の名前）に直して handoff に入れる。

## 2. ゴール・ゴールでないもの

ゴール:
- MSP の脂質の行が、LBM の行と同じ規則で handoff に入る（同じクラス名・同じ書式の名前）。
- 直せない行は、理由ごとの件数とともに落とす。何件がどの経路で入り、何件がどの理由で落ちたかを出力に出す。
- テストで、MSP の行（慣用名・ClassyFire のクラス・組成式）から脂質が引けることを確かめる。

ゴールでないもの:
- MSP 由来かどうかを出力に書くこと（ユーザーの指示: テストで MSP から引けることを確かめれば、明示は要らない）。
- パスウェイの結果を「良くする」こと。結果は観測どおりに報告し、変化が無いことも結果の一つとする（ユーザーの指示）。
- 脂質でない代謝物（ペプチド・テルペン等）をパスウェイに入れること。
- 外部のデータベース（LIPID MAPS の構造の表など）を引くこと。オフラインで、リポジトリの中の表だけで決める。
- SMILES / InChIKey からの構造の分類（RDKit の部分構造で頭部を決める等）。今回はしない。

してはいけないこと:
- クラスや組成を推測で埋めること。ClassyFire のクラスから MS-DIAL のクラスが一意に決まらない行（例 `Ceramides`）、
  組成式がそのクラスの形に合わない行は、名前で決まらない限り落とす。間違ったノードに数えるより、落として数える方がよい。
- LBM の行の扱いを変えること。いま残っている LBM の行は、そのまま同じ名前・クラスで残る。

## 3. 試しの mzTab で見えたこと（設計の前提）

陽の MSP の SML は 980 行（MS/MS 照合 553）、Ontology は 211 種類。陰は 811 行（照合 356）、169 種類。
脂質で件数の多い Ontology（陽 / 陰、括弧は照合）:

- `Phosphatidylcholines` 142（123）/ 89（59）、`Phosphatidylethanolamines` 39（37）/ 61（47）、
  `Triacylglycerols` 34（23）、`Phosphatidylserines` 18 / 25、`Phosphatidylglycerols` 5 / 36、`Phosphatidylinositols` – / 23
- `1-(1Z-alkenyl),2-acylglycerophosphoethanolamines` 23 / 17、`1-alkyl,2-acylglycero-3-phosphocholines` 8 / 5 ほか、エーテル・プラズマローゲンの系
- `Long-chain ceramides` 22 / 25、`Ceramides` 12 / 16、`Phosphosphingolipids` 14 / 19、`Sphingomyelins` 10、`Glycosyl-N-acylsphingosines` 7 / 6
- `1,2-diacylglycerols` 15、`1,3-diacylglycerols` 12、`1-monoacylglycerols` 5
- `Long-chain fatty acids` – / 61、`Very long-chain fatty acids` – / 24、`Lineolic acids and derivatives` 11 / 4
- `Acyl carnitines` 23、`N-acylethanolamines` 7、`1-acyl-sn-glycero-3-phosphocholines` 12 / 10 ほかリゾ体、`Cholesterols and derivatives` 7

名前の形は混ざっている:
- すでに略記: `PE 32:1`、`PI 32:1`、`SM d34:2`、`Cer[NS] d30:1`、`LPC 16:0; PlaSMA ID-1486`
- 英語のクラス名と鎖: `Phosphatidylcholine 16:0-22:6`、`Phosphatidylcholine lyso 18:2`、`Phosphatidylethanolamine alkyl 16:0-22:6`
- 慣用名・IUPAC 名: `Glyceryl trimyristate`、`1-Palmitoyl-2-oleoyl-sn-glycerol`、`N-Palmitoyl-D-sphingosine`、`Docosahexaenoic acid`
- 名前の無い ID: `MassSpecGymID0400493`
- 先頭に `low score: ` / `no MS2: ` が付くもの（LBM の行と同じ）

ClassyFire のクラスと実際の物質が食い違う例もある: `1-(1Z-alkenyl),2-acyl-glycerophosphocholines` の行の名前が
`…-sn-glycero-3-phosphoethanolamine`。`Lineolic acids and derivatives` に DG（`1-Palmitoyl-2-linoleoyl-rac-glycerol`、C37H68O5）が入る。

LBM の行の名前（MS-DIAL 5 の書式）: `PC 33:6`、`PE P-16:0_16:1`、`PC O-31:3`、`TG 8:0_8:0_15:2`、`DG 16:0_16:0`、
`FA 25:0`、`Cer 18:2;O2/16:0`、`SM 32:7;O3`、`LPE 17:1`、`CAR 28:6`、`CE 22:6`、`HexCer 16:2;O2/20:5`、`NAE 19:4`。

## 4. 設計

### 4.1 どこに入れるか

handoff の行ごとの判定（`_drop_reason`）の前に、行を MS-DIAL の形に直す段を 1 つ置く。入力は行の名前・Ontology・組成式、
出力は `(MS-DIAL の名前, LbmClass)` か、落とす理由。dataset の経路（`build_from_dataset`。`ds.feature_annotations`）と
ARF の経路（`build_from_arf`。catalog の行）の両方に通す。組成式はどちらの経路にもすでにある（dataset は
`feature_annotations` の `chemical_formula`、`mztab/dataset_state.py:633`。ARF は catalog の `Formula`、`arf2/reader.py:110`）。

経路 0: 行の Ontology が LbmClass なら、いまの規則だけで決める（名前が読めれば残し、読めなければ `unparsable_name` で落とす）。
経路 1・2 には回さない。LBM の行の件数と名前を変えないため（設計者の判断。2026-10-30 の試しでは、この `unparsable_name` は
陽 44・陰 109 行で、すべて LBM の行）。経路 1・2 は、Ontology が LbmClass でない行（いまの `non_msdial_ontology`）だけに使う。

### 4.2 経路 1: 名前から読む

名前から先頭の `low score: ` / `no MS2: ` と `;` 以降（`; PlaSMA ID-1486` 等）を外し、次の順に試す:

1. そのまま、または小さな書き換え（`Cer[NS] d30:1` → `Cer 30:1;O2` のような、MS-DIAL の古い書式から 5 の書式へ）で
   `normalize_lipid_name`（pygoslin）が読める。
2. 英語のクラス名と鎖の形（`Phosphatidylcholine 16:0-22:6` → `PC 16:0_22:6`、`Phosphatidylcholine lyso 18:2` → `LPC 18:2`、
   `Phosphatidylethanolamine alkyl 16:0-22:6` → `PE O-16:0_22:6`、`… alkenyl …` → `PE P-…`）を、小さな対応表で略記に直してから読む。

読めたら、クラスは読めた名前から LbmClass を決める。このとき行の Ontology（ClassyFire）は使わない。名前と ClassyFire の
クラスが食い違うとき（3 節の PE / PC の例）、名前を採る。GOSLIN のクラスが LbmClass の複数に当たるとき（Cer → Cer_NS / NDS / AS …、
HexCer、SM の `;O3` など）は、名前が示す構造（`Cer 18:1;O2/16:0` の長鎖塩基と脂肪酸の水酸基、古い書式の `Cer[NS]`）で一意に
決まるときだけ決め、決まらなければ推測せず落とす（`ambiguous_class`。和の組成だけの `Cer 34:1;O2` は落ちる）。

名前の書式: 出す名前は MS-DIAL 5 の書式（3 節の LBM の行と同じ。`PE P-…` を GOSLIN の `PE O-…` に変えない）で、元の行の
`low score: ` / `no MS2: ` は外さずに前に付け直す（LBM の行もこれを付けたまま handoff を通っている。`handoff.py:156, 206`）。
`normalize_lipid_name` は読めるかどうかと構造を取るのに使い、その出力の文字列（GOSLIN の書式）をそのまま名前にはしない。

### 4.3 経路 2: クラスと組成式から作る

名前で読めない行（慣用名・IUPAC 名・ID）に使う。

1. ClassyFire の Ontology → LbmClass の対応表（リポジトリに置く小さな TSV）で、LbmClass を 1 つに決める。
   一意に決まらない Ontology（`Ceramides`、`Long-chain ceramides`、`Phosphosphingolipids` など、MS-DIAL では複数のクラスに分かれるもの）と、
   表に無い Ontology は、この経路では決めない（落とす理由 `unmapped_ontology` / `ambiguous_ontology`）。
2. LbmClass ごとの基準の種（例: PC は `PC 34:1` = C42H82NO8P）から、組成式で炭素数の合計と二重結合の数を出す:
   鎖の炭素 n = C − (C₀ − n₀)、二重結合 = DB₀ + ((H₀ + 2(n − n₀)) − H) / 2。C・H 以外の元素の数が基準と違う、
   二重結合が整数でない・負になる行は、そのクラスの形に合わないとして落とす（`formula_mismatch`）。
   3 節の `Lineolic acids and derivatives` に入った DG（O5）は、FA（O2）の形に合わないので落ちる。
3. 名前は和の組成（`TG 54:3`、`PE 38:4`、`PE P-36:4`、`SM 34:1;O2`、`CAR 18:1`、`FA 22:6`）。鎖ごとの内訳は出さない。
   接頭辞は経路 1 と同じく元の行のものを付け直す。組成式の無い行は `no_formula` で落とす。

例: `Glyceryl trioleate`（C57H104O6、`Triacylglycerols`）→ `TG 54:3`。`1-Palmitoyl-2-oleoyl-sn-glycerol`（C37H70O5、
`1,2-diacylglycerols`）→ `DG 34:1`。`MassSpecGymID0400493`（C36H72NO8P、`Phosphatidylcholines`）→ `PC 28:0`。

### 4.4 そのあと

直した行は、いまの `_drop_reason` と重複の扱い（同じ名前は平均強度が最大の 1 行）にそのまま通す。LBM の行と MSP の行が
同じ名前（例 `TG 54:3`）になったら、いまの重複の規則で 1 行に絞る（どちらの由来かは問わない）。

### 4.5 出力の記録

行はどれか 1 つの区分にだけ入る（設計者の判断）: 残る行は `kept_by_route` の `as_is`（経路 0）・`name`（経路 1）・
`class_formula`（経路 2）のどれか。落ちる行は次の理由のどれか 1 つ:
- `unannotated`（いまと同じ）、`unparsable_name`（経路 0 の行で名前が読めない。いまと同じ定義）、
- 経路 1 で読めて経路 2 に回らなかった行の `ambiguous_class`、
- 経路 1 で読めず経路 2 で決まらなかった行の `unmapped_ontology`・`ambiguous_ontology`・`no_formula`・`formula_mismatch`
  （経路 2 の最初に当たった理由）。
`dropped_non_msdial_ontology` はこれらに分かれるので出さない（受け取る側の文書を直す）。重複で絞った行は、いまと同じ
`dropped_duplicate_names` に数える。`kept_by_route` は重複で絞る前の数で、`n_kept` = その合計 − `dropped_duplicate_names`。
全行 = 残る行（絞る前）+ 落ちる理由の合計、が成り立つことをテストで確かめる。

重複は、いまと同じく名前の文字列で判定する（設計者の判断。`low score: TG 54:3` と `TG 54:3` は別の名前のまま。LBM の行どうしでも
いまそうなっていて、変えると LBM の行の結果が変わる）。

## 5. テスト

- 単体: 3 節の名前の形ごとに、期待する `(名前, LbmClass)` か落とす理由を固定する。少なくとも:
  略記（`PE 32:1`）、古い書式（`Cer[NS] d30:1`）、英語のクラス名（`Phosphatidylcholine 16:0-22:6`、`… lyso 18:2`、`… alkyl …`）、
  慣用名 + ClassyFire + 組成式（`Glyceryl trioleate`）、ID + ClassyFire + 組成式（`MassSpecGymID…`）、
  一意に決まらない Ontology（`Ceramides`）、組成式が合わない行（`Lineolic acids` の DG）、名前と ClassyFire が食い違う行（名前を採る）。
- handoff の経路: MSP の行（慣用名・ClassyFire・組成式）を含む `feature_annotations` から `build_from_dataset` を通し、
  その行が期待する名前とクラスで残ることを確かめる（ユーザーの指示: MSP から脂質が引けることをテストで確かめる）。
  ARF の経路も同じ形で 1 本。
- 対応表の整合: 対応表の LbmClass がすべて `MSDIAL_LIPID_CLASSES` にある。基準の種の組成式が、その名前の組成と合う。
- 試しの mzTab（`_v8_trial/pipeline_r2`）で handoff を回し、経路ごとの件数と落とした理由の件数を、MTBLS6965 の試しの記録
  （`docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md`）に足す。数字は観測どおりに書く。

## 6. 文書

- ms-data-parser: `docs/` の pathway の説明（handoff の落とす規則と report のキー）、`msdial_lipid_classes.py` の docstring
  （「MSP の行は落ちる」の記述）、Vault の流れの文書。
- reanalysis-study: 手順書の工程 6 の脂質の項、工程 7 の 1 のパスウェイの主張、レポートの雛形の工程 6 の注意、
  引き継ぎのメモ（v8）にある「MSP の行は落ちる・クラスの和は LBM に偏る」を、新しい記録（`kept_by_route` と `dropped_*`）を
  書く形に直す。ライブラリ契約の「LBM と MSP を区別しない」への例外の記述は消す。手順書の版は v8 のまま直す
（ユーザーの決定 2026-11-02。v8 で回した run はまだ無い）。

## 7. 決まっていること・作る人に任せること

ユーザーが決めたこと:
- MSP と LBM を区別しない。MSP 由来の明示は出力に要らず、テストで確かめる。
- 結果は観測どおりに報告する（変化なしも結果）。
- reanalysis-study の手順書は v8 のまま直す（版を上げない）。

私（設計者）の判断:
- 経路 1（名前）を経路 2（クラスと組成式）より先にする。名前の方が鎖まで分かり、ClassyFire の食い違いに強い。
- 一意に決まらない ClassyFire のクラスは推測しないで落とす。
- 外部データベースと構造からの分類は使わない。

設計者の判断（作る人の確認を受けて足したもの）:
- 経路 1・2 は Ontology が LbmClass でない行だけに使う（4.1）。
- 経路 1 でも、LbmClass が一意に決まらない名前は推測しないで落とす（4.2）。
- 名前は MS-DIAL 5 の書式で、元の接頭辞を付け直す。重複は文字列で判定する（4.2、4.5）。
- 落とす理由は互いに重ならない区分にする（4.5）。

作る人に任せること:
- ClassyFire → LbmClass の対応表と、英語のクラス名の対応表の中身（3 節の脂質の Ontology を覆い、一意に決まらないものは載せない）。
  エーテル・プラズマローゲンの ClassyFire 名から `PE O-` / `PE P-` をどう決めるか。
- LbmClass ごとの基準の種と、MS-DIAL 5 の和の組成の書式（`;O2` 等の付け方）。LBM の行の名前（3 節）に合わせる。
- 組成式の読み方（同位体の標識 `[2H]7` などを含む式は経路 2 で落としてよい）。
