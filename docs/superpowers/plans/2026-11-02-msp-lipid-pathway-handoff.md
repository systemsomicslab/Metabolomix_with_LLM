# 研究室 MSP の脂質をパスウェイの handoff に入れる — 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ontology が MS-DIAL の脂質クラスでない行（研究室 MSP の行）を、名前か「ClassyFire のクラス＋組成式」から MS-DIAL のクラスと略記の名前に直して handoff に入れ、直せない行は重ならない理由で数える。

**Architecture:** 新しい純ロジックのモジュール `metabolomix/pathway/msp_lipid.py` が 2 つの経路（名前から読む・クラスと組成式から作る）と対応表を持つ。`handoff.py` は行ごとの判定を `_drop_reason` から `resolve_row` に替え、経路 0（いまの LBM の行）・1・2 のどれで残したか、どの理由で落としたかを report に数える。最後に試しの mzTab で件数を測り、両リポの文書を直す。

**Tech Stack:** Python 3、pygoslin 2.2.5（`normalize_lipid_name`）、pytest。外部データベースは使わない。

**Spec:** `docs/superpowers/specs/2026-11-02-msp-lipid-pathway-handoff-design.md`

## Global Constraints

- 経路 0（Ontology が `MSDIAL_LIPID_CLASSES` にある行）の判定・名前・件数は変えない。既存の `tests/test_pathway_handoff.py` の LBM の行の期待値はそのまま通る。
- 一意に決まらないクラスは推測しないで落とす（ClassyFire の `Ceramides` 等、名前の `Cer 34:1;O2` 等）。
- 出す名前は MS-DIAL 5 の書式（`PC 34:1`、`PE P-32:1`、`PC O-34:1`、`SM 34:1;O2`、`Cer 18:1;O2/16:0`）。GOSLIN の出力の文字列を名前にしない。元の名前の限定子（`low score: ` / `no MS2: ` 等、`lipid_identity._MSDIAL_QUALIFIER_RE` が外すもの）は元の綴りのまま前に付け直す。
- 重複は名前の文字列で判定する（いまの `_dedup_by_name` のまま）。
- 落とす理由は互いに重ならない: `unannotated`・`unparsable_name`（経路 0）・`ambiguous_class`（経路 1）・`unmapped_ontology`・`ambiguous_ontology`・`no_formula`・`formula_mismatch`（経路 2）。`dropped_non_msdial_ontology` は出さない。
- MSP 由来かどうかは出力に書かない。オフライン（ネットワーク・外部 DB なし）。
- 研究室のライブラリのパス・ファイル名を追跡対象に書かない。
- reanalysis-study の手順書は v8 のまま直す（版を上げない）。

## Review Focus

1. 名前に候補の区切り `|` や `; PlaSMA ID-1486` のような後ろの注記がある — 先頭の候補・`;` の前だけで読む（Task 2 のテスト）。
2. 組成式に同位体の標識（`C39H67[2H]7NO10P`）や読めない文字がある — 例外を投げず `formula_mismatch` で落とす（Task 1 のテスト）。
3. ClassyFire の名前の前後の空白 — `strip()` してから引く。大文字小文字は区別する（表の綴りどおり）（Task 1 のテスト）。
4. MSP の行と LBM の行が同じ名前になる — 1 行に絞られ、`kept_by_route` の合計 − `dropped_duplicate_names` = `n_kept`（Task 3 のテスト）。
5. Ontology が `null` / 空で名前が略記の MSP の行 — いまと同じく `unannotated`（経路 1 に回さない）（Task 3 のテスト）。

---

### Task 1: クラスと組成式から作る（経路 2）

**Files:**
- Create: `metabolomix/pathway/msp_lipid.py`
- Test: `tests/test_pathway_msp_lipid.py`

**Interfaces:**
- Produces:
  - `parse_formula(formula: str | None) -> dict[str, int] | None` — `C57H104O6` → `{"C": 57, "H": 104, "O": 6}`。元素記号＋数だけの式でなければ None（`[2H]`、電荷、空）。
  - `CLASSYFIRE_TO_CLASS: dict[str, str]` — ClassyFire の名前 → クラスの鍵（下の `CLASS_FORMS` の鍵）。
  - `AMBIGUOUS_CLASSYFIRE: frozenset[str]` — MS-DIAL で複数のクラスに分かれる ClassyFire の名前。
  - `CLASS_FORMS: dict[str, ClassForm]`、`@dataclass(frozen=True) class ClassForm: lbm_class: str; prefix: str; suffix: str; ref_n: int; ref_db: int; ref_formula: str` — 名前は `f"{prefix}{n}:{db}{suffix}"`。
  - `from_class_formula_key(key: str, formula: str | None) -> tuple[str, str] | str` — 鍵（`CLASS_FORMS` の鍵）から直接作る。理由は `"no_formula" | "formula_mismatch"`。
  - `from_class_formula(ontology: str | None, formula: str | None) -> tuple[str, str] | str` — 成功は `(名前（限定子なし）, LbmClass)`、失敗は理由 `"unmapped_ontology" | "ambiguous_ontology" | "no_formula" | "formula_mismatch"`（この順に見る）。

最低限の表（`CLASSYFIRE_TO_CLASS` → `CLASS_FORMS`。基準の種と組成式。設計書 3 節の脂質の Ontology を覆う。足すのは自由、ここに無い ClassyFire 名で迷うものは載せない）:

| ClassyFire | 鍵 | LbmClass | 名前の形 | 基準 |
|---|---|---|---|---|
| Phosphatidylcholines | PC | PC | `PC n:db` | 34:1 = C42H82NO8P |
| Phosphatidylethanolamines | PE | PE | `PE n:db` | 34:1 = C39H76NO8P |
| Phosphatidylserines | PS | PS | `PS n:db` | 34:1 = C40H76NO10P |
| Phosphatidylglycerols | PG | PG | `PG n:db` | 34:1 = C40H77O10P |
| Phosphatidylinositols, 1-phosphatidyl-1D-myo-inositols | PI | PI | `PI n:db` | 34:1 = C43H81O13P |
| 1,2-diacylglycerol-3-phosphates | PA | PA | `PA n:db` | 34:1 = C37H71O8P |
| 1-acyl-sn-glycero-3-phosphocholines, 2-acyl-sn-glycero-3-phosphocholines | LPC | LPC | `LPC n:db` | 16:0 = C24H50NO7P |
| 1-acyl-sn-glycero-3-phosphoethanolamines | LPE | LPE | `LPE n:db` | 16:0 = C21H44NO7P |
| 1-(1Z-alkenyl),2-acylglycerophosphoethanolamines | PE P- | EtherPE | `PE P-n:db` | 32:1 = C37H72NO7P |
| 1-alkyl,2-acylglycerophosphoethanolamines | PE O- | EtherPE | `PE O-n:db` | 32:1 = C37H74NO7P |
| 1-(1Z-alkenyl),2-acyl-glycerophosphocholines | PC P- | EtherPC | `PC P-n:db` | 34:1 = C42H82NO7P |
| 1-alkyl,2-acylglycero-3-phosphocholines | PC O- | EtherPC | `PC O-n:db` | 34:1 = C42H84NO7P |
| Triacylglycerols | TG | TG | `TG n:db` | 54:3 = C57H104O6 |
| 1,2-diacylglycerols, 1,3-diacylglycerols | DG | DG | `DG n:db` | 34:1 = C37H70O5 |
| 1-monoacylglycerols, 2-monoacylglycerols | MG | MG | `MG n:db` | 18:1 = C21H40O4 |
| Long-chain fatty acids, Very long-chain fatty acids, Medium-chain fatty acids, Lineolic acids and derivatives | FA | FA | `FA n:db` | 18:1 = C18H34O2 |
| Acyl carnitines | CAR | CAR | `CAR n:db` | 16:0 = C23H45NO4 |
| N-acylethanolamines | NAE | NAE | `NAE n:db` | 16:0 = C18H37NO2 |
| Sphingomyelins | SM | SM | `SM n:db;O2` | 34:1 = C39H79N2O6P |

`AMBIGUOUS_CLASSYFIRE` は少なくとも `Ceramides`、`Long-chain ceramides`、`Very long-chain ceramides`、`Phosphosphingolipids`、`Glycosyl-N-acylsphingosines`、`Cholesterols and derivatives`。

- [ ] **Step 1: Write the failing tests** in `tests/test_pathway_msp_lipid.py`

```python
from metabolomix.pathway import msp_lipid as ml
from metabolomix.pathway.msdial_lipid_classes import MSDIAL_LIPID_CLASSES

def test_parse_formula():
    assert ml.parse_formula("C57H104O6") == {"C": 57, "H": 104, "O": 6}
    assert ml.parse_formula("C39H79N2O6P") == {"C": 39, "H": 79, "N": 2, "O": 6, "P": 1}
    for bad in (None, "", "C39H67[2H]7NO10P", "C6H12O6+", "abc"):
        assert ml.parse_formula(bad) is None

def test_from_class_formula_examples():
    assert ml.from_class_formula("Triacylglycerols", "C57H104O6") == ("TG 54:3", "TG")
    assert ml.from_class_formula("1,2-diacylglycerols", "C37H70O5") == ("DG 34:1", "DG")
    assert ml.from_class_formula("Phosphatidylcholines", "C36H72NO8P") == ("PC 28:0", "PC")
    assert ml.from_class_formula("Very long-chain fatty acids", "C22H32O2") == ("FA 22:6", "FA")
    assert ml.from_class_formula("1-(1Z-alkenyl),2-acylglycerophosphoethanolamines", "C37H72NO7P") == ("PE P-32:1", "EtherPE")
    assert ml.from_class_formula("1-alkyl,2-acylglycero-3-phosphocholines", "C42H84NO7P") == ("PC O-34:1", "EtherPC")
    assert ml.from_class_formula("Sphingomyelins", "C39H77N2O6P") == ("SM 34:2;O2", "SM")
    assert ml.from_class_formula("  Triacylglycerols ", "C57H104O6") == ("TG 54:3", "TG")   # 前後の空白

def test_from_class_formula_reasons():
    assert ml.from_class_formula("Benzoic acid esters", "C16H22O4") == "unmapped_ontology"
    assert ml.from_class_formula("Ceramides", "C34H67NO3") == "ambiguous_ontology"
    assert ml.from_class_formula("Triacylglycerols", None) == "no_formula"
    assert ml.from_class_formula("Lineolic acids and derivatives", "C37H68O5") == "formula_mismatch"   # DG は FA の形でない
    assert ml.from_class_formula("Phosphatidylserines", "C39H67[2H]7NO10P") == "formula_mismatch"
    assert ml.from_class_formula("Triacylglycerols", "C57H103O6") == "formula_mismatch"               # 二重結合が整数でない

def test_tables_are_consistent():
    for key, form in ml.CLASS_FORMS.items():
        assert form.lbm_class in MSDIAL_LIPID_CLASSES, key
        assert ml.from_class_formula_key(key, form.ref_formula) == (
            f"{form.prefix}{form.ref_n}:{form.ref_db}{form.suffix}", form.lbm_class)
    assert set(ml.CLASSYFIRE_TO_CLASS.values()) <= set(ml.CLASS_FORMS)
    assert not set(ml.CLASSYFIRE_TO_CLASS) & ml.AMBIGUOUS_CLASSYFIRE
```

（`from_class_formula` は ClassyFire の名前を鍵に直してから `from_class_formula_key` を呼ぶ。）

- [ ] **Step 2: Run to verify they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_msp_lipid.py -q`
Expected: FAIL（`ModuleNotFoundError: metabolomix.pathway.msp_lipid`）

- [ ] **Step 3: Implement** `metabolomix/pathway/msp_lipid.py`（モジュールの docstring に設計書の 4.3 の式と「推測しない」を書く）

組成: 鎖の炭素 n = C − (C₀ − n₀)、二重結合 db = db₀ + ((H₀ + 2(n − n₀)) − H) / 2。C・H 以外の元素の数が基準と違う、db が整数でない・負、n < 1 なら `formula_mismatch`。

- [ ] **Step 4: Run to verify they pass**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_msp_lipid.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add metabolomix/pathway/msp_lipid.py tests/test_pathway_msp_lipid.py
git commit -m "feat(pathway): ClassyFire のクラスと組成式から MS-DIAL の脂質クラスと和の組成を作る"
```

---

### Task 2: 名前から読む（経路 1）

**Files:**
- Modify: `metabolomix/pathway/msp_lipid.py`
- Test: `tests/test_pathway_msp_lipid.py`

**Interfaces:**
- Consumes: `normalize_lipid_name(name) -> dict`（`metabolomix/msdial/lipid_identity.py`。`parse_ok` だけを使う）。
- Produces:
  - `split_qualifier(name: str) -> tuple[str, str]` — `("low score: ", "PC 34:1")`。限定子が無ければ `("", name.strip())`。限定子は元の綴り（`lipid_identity._MSDIAL_QUALIFIER_RE` と同じ範囲）。
  - `from_name(name: str | None) -> tuple[str, str] | str | None` — 成功は `(名前（限定子なし）, LbmClass)`、クラスが一意でなければ `"ambiguous_class"`、読めなければ None（経路 2 に回す）。入力は限定子を外した名前。

決め方:
1. `|` の前、`;` の前（ただし `;O2` のような MS-DIAL 5 の修飾は残す: `;` の後ろが空白を含む注記のときだけ切る。`LPC 16:0; PlaSMA ID-1486` → `LPC 16:0`）。
2. 古い書式: `Cer[NS] d30:1` → `Cer 30:1;O2`（`[NS]` は Cer_NS、`[NDS]` は Cer_NDS として覚える）、`SM d34:2` → `SM 34:2;O2`、`t` は `;O3`。
3. 英語のクラス名: `Phosphatidylcholine` → `PC`、`Phosphatidylethanolamine` → `PE`、`Phosphatidylserine` → `PS`、`Phosphatidylglyceride` / `Phosphatidylglycerol` → `PG`、`Phosphatidylinositol` → `PI`、`Phosphatidic acid` → `PA`。続く `lyso` は `L` を前に（`LPC`）、`alkyl` は `O-`、`alkenyl` は `P-`（`PE O-16:0_22:6`）。鎖の区切り `-` は `_`。
4. 直した文字列を `normalize_lipid_name` で読めなければ None。
5. LbmClass は先頭の頭部（と `O-` / `P-`）から: `PC`/`PE`/`PS`/`PG`/`PI`/`PA`/`LPC`/`LPE`/`LPG`/`LPI`/`LPS`/`LPA`/`TG`/`DG`/`MG`/`FA`/`CAR`/`CE`/`NAE`/`SM` は同名、`PC O-`/`PC P-` は EtherPC、`PE O-`/`PE P-` は EtherPE。`Cer` と `HexCer` は、古い書式の `[NS]`/`[NDS]` か、`Cer 18:1;O2/16:0` 形（長鎖塩基が `;O2`、脂肪酸に `(2OH)` 等の水酸基の注記が無い）で、塩基の二重結合 ≥ 1 なら `Cer_NS`/`HexCer_NS`、0 なら `Cer_NDS`/`HexCer_NDS`。それ以外の `Cer`/`HexCer`（和の組成だけ・`;O3`・水酸基付き）は `"ambiguous_class"`。表に無い頭部は None。

- [ ] **Step 1: Write the failing tests**

```python
def test_split_qualifier():
    assert ml.split_qualifier("low score: PC 34:1") == ("low score: ", "PC 34:1")
    assert ml.split_qualifier("no MS2: Myristic acid") == ("no MS2: ", "Myristic acid")
    assert ml.split_qualifier(" TG 54:3 ") == ("", "TG 54:3")

def test_from_name_examples():
    assert ml.from_name("PE 32:1") == ("PE 32:1", "PE")
    assert ml.from_name("LPC 16:0; PlaSMA ID-1486") == ("LPC 16:0", "LPC")
    assert ml.from_name("Cer 24:1;O2|Cer 12:0;O2/12:1") == "ambiguous_class"
    assert ml.from_name("Cer[NS] d30:1") == ("Cer 30:1;O2", "Cer_NS")
    assert ml.from_name("SM d34:2") == ("SM 34:2;O2", "SM")
    assert ml.from_name("Phosphatidylcholine 16:0-22:6") == ("PC 16:0_22:6", "PC")
    assert ml.from_name("Phosphatidylcholine lyso 18:2") == ("LPC 18:2", "LPC")
    assert ml.from_name("Phosphatidylethanolamine alkyl 16:0-22:6") == ("PE O-16:0_22:6", "EtherPE")
    assert ml.from_name("Phosphatidylethanolamine alkenyl 16:0-20:2") == ("PE P-16:0_20:2", "EtherPE")
    assert ml.from_name("Cer 18:1;O2/16:0") == ("Cer 18:1;O2/16:0", "Cer_NS")
    assert ml.from_name("Cer 18:0;O2/16:0") == ("Cer 18:0;O2/16:0", "Cer_NDS")

def test_from_name_does_not_guess():
    assert ml.from_name("Cer 34:1;O2") == "ambiguous_class"
    assert ml.from_name("Cer 18:1;O2/18:4(2OH)") == "ambiguous_class"
    for unreadable in ("Glyceryl trioleate", "MassSpecGymID0400493", "Docosahexaenoic acid", "", None):
        assert ml.from_name(unreadable) is None
```

- [ ] **Step 2: Run to verify they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_msp_lipid.py -q`
Expected: FAIL（`AttributeError: ... split_qualifier`）

- [ ] **Step 3: Implement** `split_qualifier` と `from_name`（上の決め方。書き換えの表はモジュールの定数）

- [ ] **Step 4: Run to verify they pass**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_msp_lipid.py tests/test_lipid_identity.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add metabolomix/pathway/msp_lipid.py tests/test_pathway_msp_lipid.py
git commit -m "feat(pathway): 研究室 MSP の脂質名を MS-DIAL 5 の略記と LbmClass に読む"
```

---

### Task 3: handoff に組み込む

**Files:**
- Modify: `metabolomix/pathway/handoff.py`（`_drop_reason` → `resolve_row`、`build_from_arf`・`build_from_dataset` の report）
- Modify: `metabolomix/pathway/tools.py:35-50`（`_handoff_note`）
- Test: `tests/test_pathway_handoff.py`、`tests/test_pathway_tools.py`

**Interfaces:**
- Consumes: Task 1・2 の `split_qualifier`・`from_name`・`from_class_formula`。
- Produces:
  - `@dataclass(frozen=True) class RowResolution: route: str | None; name: str | None; ontology: str | None; reason: str | None` — 残す行は `route` が `"as_is" | "name" | "class_formula"`、落とす行は `reason`。
  - `resolve_row(name, ontology, formula) -> RowResolution` — 順に: `_annotated` でなければ `unannotated` → Ontology が `MSDIAL_LIPID_CLASSES` にあれば経路 0（読めれば `as_is` で名前・Ontology はそのまま、読めなければ `unparsable_name`）→ `from_name(限定子を外した名前)`（タプルなら `name`、`"ambiguous_class"` ならその理由）→ `from_class_formula(ontology, formula)`（タプルなら `class_formula`、文字列ならその理由）。残す名前は `限定子 + 直した名前`。
  - report のキー（両経路）: `kept_by_route: {"as_is", "name", "class_formula"}`、`dropped_unannotated`・`dropped_unparsable_name`・`dropped_ambiguous_class`・`dropped_unmapped_ontology`・`dropped_ambiguous_ontology`・`dropped_no_formula`・`dropped_formula_mismatch`（いずれも 0 から）、既存の `dropped_duplicate_names`・`duplicate_name_examples`・`n_kept`、ARF は `dropped_excluded_spots`。`dropped_non_msdial_ontology` は無い。
  - 組成式: dataset は `ann.get("chemical_formula")`、ARF は `spot.get("Formula")`。

- [ ] **Step 1: Write / update the failing tests** in `tests/test_pathway_handoff.py`

既存のテストで `dropped_non_msdial_ontology` を見るもの（`test_handoff_keeps_msdial_lipids`・`test_handoff_drops_classyfire_ontology`・`test_handoff_report_counts_the_kept_features`・`test_handoff_drops_unparsable_names`・`test_arf_handoff_applies_the_same_lipid_checks`）は、新しい理由に書き直す（`FA 16:0` + `Long-chain fatty acids` は経路 1 で `FA 16:0` として残る点に注意）。新しいテスト:

```python
def test_msp_rows_reach_the_handoff():
    """研究室 MSP の行（慣用名・ID・英語のクラス名 + ClassyFire + 組成式）から脂質が引ける（ユーザーの指示）。"""
    t = handoff.build_from_dataset(**_ds_with(_KEEP + [
        {"name": "Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"},
        {"name": "low score: MassSpecGymID0400493", "ontology": "Phosphatidylcholines", "chemical_formula": "C36H72NO8P"},
        {"name": "Phosphatidylethanolamine alkenyl 16:0-20:2", "ontology": "1-(1Z-alkenyl),2-acylglycerophosphoethanolamines",
         "chemical_formula": "C41H78NO7P"},
    ]))
    assert t.names == ["PC 38:3", "PE 18:0_20:4", "TG 54:3", "low score: PC 28:0", "PE P-16:0_20:2"]
    assert t.ontologies == ["PC", "PE", "TG", "PC", "EtherPE"]
    assert t.report["kept_by_route"] == {"as_is": 2, "name": 1, "class_formula": 2}

def test_name_wins_over_classyfire():
    t = handoff.build_from_dataset(**_ds_with(_KEEP + [
        {"name": "PE 36:7", "ontology": "1-(1Z-alkenyl),2-acyl-glycerophosphocholines", "chemical_formula": "C41H68NO8P"}]))
    assert t.names[-1] == "PE 36:7" and t.ontologies[-1] == "PE"

def test_drop_reasons_partition_the_rows():
    rows = _KEEP + [
        {"name": None, "ontology": "Triacylglycerols"},                                            # unannotated
        {"name": "PE 32:1", "ontology": "null"},                                                  # unannotated（経路 1 に回さない）
        {"name": "3-Sulfamoylbenzoic acid", "ontology": "PC"},                                    # unparsable_name
        {"name": "Cer 34:1;O2", "ontology": "Ceramides", "chemical_formula": "C34H67NO3"},        # ambiguous_class
        {"name": "Dibutyl phthalate", "ontology": "Benzoic acid esters", "chemical_formula": "C16H22O4"},   # unmapped_ontology
        {"name": "N-Palmitoyl-D-sphingosine", "ontology": "Long-chain ceramides", "chemical_formula": "C34H67NO3"},  # ambiguous_ontology
        {"name": "Glyceryl trimyristate", "ontology": "Triacylglycerols"},                       # no_formula
        {"name": "1-Palmitoyl-2-linoleoyl-rac-glycerol", "ontology": "Lineolic acids and derivatives",
         "chemical_formula": "C37H68O5"},                                                         # formula_mismatch
        {"name": "Glyceryl trioleate", "ontology": "Triacylglycerols", "chemical_formula": "C57H104O6"},
        {"name": "TG 54:3", "ontology": "TG"},                                                    # 同じ名前 → 重複
    ]
    t = handoff.build_from_dataset(**_ds_with(rows))
    r = t.report
    dropped = {k: v for k, v in r.items() if k.startswith("dropped_") and k != "dropped_duplicate_names"}
    assert dropped == {"dropped_unannotated": 2, "dropped_unparsable_name": 1, "dropped_ambiguous_class": 1,
                       "dropped_unmapped_ontology": 1, "dropped_ambiguous_ontology": 1, "dropped_no_formula": 1,
                       "dropped_formula_mismatch": 1}
    assert sum(r["kept_by_route"].values()) + sum(dropped.values()) == len(rows)
    assert r["n_kept"] == sum(r["kept_by_route"].values()) - r["dropped_duplicate_names"] == 3
    assert "dropped_non_msdial_ontology" not in r

def test_arf_reads_formula_for_msp_rows():
    kw = _arf_inputs()
    kw["catalog"] = {**kw["catalog"], 2: {"Name": "Glyceryl trioleate", "Ontology": "Triacylglycerols", "Formula": "C57H104O6"}}
    t = handoff.build_from_arf(**kw)
    assert "TG 54:3" in t.names and t.report["kept_by_route"]["class_formula"] == 1
```

（`_ds_with` は `{**a, "ambiguous": False}` で注釈を作るので、`chemical_formula` はそのまま通る。`_arf_inputs` の catalog の鍵 2 が既存のどの行かは実装時に確かめ、空いている鍵を使う。）

`tests/test_pathway_tools.py` の `_handoff_note` のテストは、新しいキーで書き直す: 注意の文は落とした件数の合計と `kept_by_route` を示し、すべて 0 のときは None。

- [ ] **Step 2: Run to verify they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_handoff.py tests/test_pathway_tools.py -q`
Expected: FAIL（`kept_by_route` が無い、MSP の行が残らない）

- [ ] **Step 3: Implement** `resolve_row` と両経路の組み込み、`_handoff_note` の文（「MSP の同定の多くが除かれる」の文は消し、経路ごとの残りと理由ごとの落ちた件数を述べる）。`_drop_reason` は消すか `resolve_row` の経路 0 に吸収する。docstring（モジュール・`msdial_lipid_classes.py`）の「MSP の行は落ちる」を直す。

- [ ] **Step 4: Run to verify they pass**

Run: `.venv/Scripts/python.exe -m pytest tests/test_pathway_handoff.py tests/test_pathway_tools.py tests/test_pathway_msp_lipid.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add metabolomix/pathway/handoff.py metabolomix/pathway/tools.py metabolomix/pathway/msdial_lipid_classes.py tests/test_pathway_handoff.py tests/test_pathway_tools.py
git commit -m "feat(pathway): handoff が研究室 MSP の脂質を名前かクラスと組成式から入れ、落とす理由を重ならない区分で数える"
```

---

### Task 4: 試しのデータで測り、ms-data-parser の文書を直す

**Files:**
- Modify: `docs/output_format/pathway.md:12-27`（落とす規則・report のキー）
- Modify: `docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md`（測った件数を足す）
- Modify: Vault `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`（handoff の記述。リポジトリの外。コミットしない）

**Interfaces:**
- Consumes: Task 3 の report のキー。

- [ ] **Step 1: 試しの mzTab で handoff の件数を測る**

`_v8_trial\pipeline_r2\{pos,neg}\pipeline_*` を `dataset_handle_open(pipeline_path=...)` で開き（`LIPIDMIX_CONFIG=C:\Users\yuu18\Metabolomix_with_LLM\lipidmix.local.toml`）、前処理（`conservative-v1`）の後、`metabolomix.pathway.handoff.build_from_dataset` を `pathway/tools.py` の `_build_dataset` と同じ引数で呼ぶ（群は R17 と P17）。上流（MSLipidMapper）は呼ばない。スクリプトはセッションの scratchpad に置き、追跡対象にしない。
Expected: pos・neg それぞれの `kept_by_route`・`dropped_*`・`n_kept` が出る。

- [ ] **Step 2: 件数を試しの記録に足す**

`2026-10-28-reanalysis-v8-trial.md` に節「MSP の脂質の handoff（2026-11-02）」を足し、前（陽 251・陰 388、MSP の行は陽 864・陰 758 が落ちる）と後の数を、観測どおりに並べる（増えても増えなくても、そのまま書く）。理由ごとの件数と、`unmapped_ontology` の多い ClassyFire 名の上位 5 つ（件数つき）も書く。

- [ ] **Step 3: `pathway.md` と Vault を新しい規則・キーに直す**

- [ ] **Step 4: 全テスト**

Run: `.venv/Scripts/python.exe -m pytest tests -q`
Expected: PASS（3914 + 新しいテスト）

- [ ] **Step 5: Commit**

```bash
git add docs/output_format/pathway.md docs/superpowers/specs/2026-10-28-reanalysis-v8-trial.md
git commit -m "docs: MSP の脂質の handoff の規則と、試しのデータでの件数を記す"
```

---

### Task 5: reanalysis-study の手順書 v8 を直す

**Files**（リポジトリ `C:\Users\yuu18\reanalysis-study`。main から新しいブランチ `feat/protocol-v8-msp-pathway` の worktree で作業する）:
- Modify: `protocol/steps/step-6.md:18-25`
- Modify: `protocol/steps/step-7.md:58-62`
- Modify: `protocol/templates/report.md:69`
- Modify: `protocol/reanalysis-protocol.md:118-126`（ライブラリ契約の「区別しない」への例外の記述を消す）
- Modify: `docs/notes/2026-10-28-v8-handoff.md:55-62`
- Test: `tests/test_protocol_docs.py`（変える必要があれば）

**Interfaces:**
- Consumes: Task 3 の report のキー（`kept_by_route`・7 つの `dropped_*`・`dropped_duplicate_names`・`n_kept`、`notes.handoff`）。

- [ ] **Step 1: 文書を直す** — 「MSP の行は ClassyFire の Ontology で落ちる・クラスの和は LBM に偏る」を、「handoff の `kept_by_route` と理由ごとの `dropped_*` を判断ログとレポートに書く。クラスの和は handoff に残った行の和で、落ちた行（とくに `ambiguous_*`・`unmapped_ontology`）の分だけクラスの総量より小さい」に替える。`dropped_non_msdial_ontology` の名前を消す。版の表記は v8 のまま。
- [ ] **Step 2: テスト**

Run: `C:/Python314/python.exe -m pytest -q`
Expected: PASS（383）

- [ ] **Step 3: Commit**

```bash
git add protocol docs/notes/2026-10-28-v8-handoff.md
git commit -m "docs: 手順書 v8 の handoff の記述を、MSP の脂質も入る新しい件数の記録に合わせる"
```
