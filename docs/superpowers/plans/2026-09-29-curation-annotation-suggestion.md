# キュレーションの候補付け Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** キュレーションで「間違い」になったスポットと未注釈スポットに、① MS-DIAL の下位候補・② 閾値を緩めた再検索・④ 別スポットのイオンとしての説明、の候補を一覧で出し、ユーザーが選んだ判断（`assign` / `redundant`）を記録して差次的エクスポートへ反映する。

**Architecture:** 既存の `lipidmix/curation/` に「候補付けレビュー」（`suggestion_id` = `cs-…`）を足す。証拠収集（`evidence.collect`）・傾向（`trend`）・フラグ記録（`flags.jsonl`、スポットごとに最新 1 行が勝つ）・ビューアの描画を再利用し、純ロジックを `candidates.py`（①②と制約）と `relations.py`（④）に、組み立てを `suggest.py` に置く。MCP 面は `curation_suggest` 1 本を足し、`curation_submit` / `curation_flags` を広げる。エクスポートへの反映は `curation/apply.py` の 1 か所から両経路に効かせる。

**Tech Stack:** Python 3.14（`C:/Python314/python.exe`）、pytest、msgpack / lz4（fixture）、pygoslin、numpy、埋め込み canvas の HTML（node があれば純関数ブロックを node で実行するテスト）。

**Spec:** `docs/superpowers/specs/2026-09-29-curation-annotation-suggestion-design.md`（先に全文を読む）。前段の設計は `docs/superpowers/specs/2026-09-28-alignment-curation-design.md`。

## Global Constraints

- テストはリポジトリルート（worktree のルート）から `C:/Python314/python.exe -m pytest tests -q` で回す。
- ツールの戻り値は `lipidmix.core.serialization.json_payload()`、行の一覧は TSV（列名 1 回）、float は丸める、`@mcp.tool` には必ず `structured_output=False`。スペクトル座標・EIC 系列は戻り値に入れない。
- 前提状態が無いときは例外でなく `mcp_errors.missing_state(state, required_tools, message)` を返す。
- 差次的エクスポートの 15 列と `CONTRACT_VERSION = 1` は変えない。判断が 0 件なら出力は現行と完全に同じ。wrong / suspect だけのときのメタ行も現行と完全に同じ。
- `flags.jsonl` は追記専用。キーは（`alignment_sha256`、`MasterAlignmentID`）、スポットごとに最新の 1 行が勝つ、`clear` で取り消す。
- fixture はテスト自身が作る。研究データ・研究室ライブラリのパスをテスト・文書に書かない。
- 既定値（spec の値そのまま）: `top_n=5`、`relation_mz_tol=0.010` Da、`rt_window` は param ファイルの `Retention time tolerance for alignment`（無ければ 0.1 分）、`relation_min_r=0.8`（有効試料 5 以上）、同位体の整合は「実測比 ≤ 期待比 × 1.5」、|Δm/z| ≥ `dmz_fail_mda`（10 mDa）はハード。
- 記録の粒度は既定で和組成（`trend.sum_composition`）。和組成で採用しても InChIKey は選んだレコードのもの。
- `assign` / `redundant` は `_tags.xml` の Misannotation を変えない（`sync_misannotation` は wrong / clear だけを見る現行のまま）。
- LLM が `curation_submit` を呼んでよいのは、ユーザーがチャットで同意したときと、ユーザーが貼った送信用テキストを取り次ぐときだけ（ツールの docstring に書く）。
- `lipidmix.core.mcp_core` から `lipidmix.tools.*` を import しない。`lipidmix/curation/*` は `lipidmix.tools.*` を import しない。
- コミットは pre-commit で全テストが走り数分かかる。自動化されたコミットはバックグラウンドで実行し、出力をファイルに落として `FAILED` を後から引く。実行中は作業ツリーを編集しない。`git stash` は使わない。
- `.gitignore` に `*.txt` がある。fixture で `.txt` を書くのは tmp だけにする（リポジトリに `.txt` を足さない）。

## Review Focus

1. **元レビューが別版の `.arf2` に対して作られている**（MS-DIAL を再実行した後）: `curation_suggest` は黙って古い `likely_wrong` を使わず、sha256 不一致のエラーを返す。→ Task 6 のテスト `test_run_suggestion_rejects_a_base_review_of_another_alignment`。
2. **貼られた送信用テキストが、保存済みの候補付けに無い `candidate` ID や対象外の `spot_id` を名指しする**（手で書き換えた・別の候補付けの文を貼った）: 何も書かずに全体をエラーにする。→ Task 6 `test_expand_entries_rejects_unknown_candidate_and_spot`。
3. **① の候補の参照レコードが store に無い**（`reference_unresolved`）のに採用される: MS-DIAL の照合結果の名前と InChIKey で記録が展開でき、エクスポートに届く。→ Task 6 `test_expand_entries_uses_msdial_identity_when_reference_unresolved`。
4. **param ファイルや `.arf` が無いフォルダ**: 既定のアダクト一覧・RT 窓で動き、相関は `None` になり MS-DIAL のリンクだけで裏付けを判定する。警告に出る。→ Task 6 `test_run_suggestion_without_param_file_and_arf_uses_defaults`。
5. **wrong / suspect だけの既存ユーザーの出力が変わる**: メタ行・エクスポート行・`arf2_annotate_identities` の `curation_flag` 列は現行と同じ。→ Task 9 `test_meta_line_unchanged_without_decisions` と `test_export_unchanged_with_wrong_only`。

---

## File Structure

| ファイル | 種別 | 責務 |
|---|---|---|
| `lipidmix/arf2/ion_features.py` | 新規 | `.arf2` Key 10（`IonFeatureCharacter`）の復号。リンクの種類を名前にする |
| `docs/schema/IonFeatureCharacter.md` | 新規 | Key 番号の正準表（上流コミットつき） |
| `lipidmix/arf2/match_results.py` | 変更 | 代表以外の候補を公開する `usable_candidates` / `load_spot_candidates` |
| `lipidmix/msdial/adducts.py` | 新規 | アダクト文字列（`[2M+FA-H]-` など）の解析と m/z ⇄ 中性質量 |
| `lipidmix/msdial/analysis_params.py` | 新規 | param ファイル（`*_param_*.txt`）から検索アダクト・RT 窓・MS1 許容幅を読む |
| `lipidmix/curation/flags.py` | 変更 | flag の種類に `assign` / `redundant` を足す。digest・判断の要約 |
| `lipidmix/curation/review.py` | 変更 | 既存レビューで `assign` / `redundant` を `decision` として持つ（flag 欄に混ぜない） |
| `lipidmix/curation/trend.py` | 変更 | `sum_composition`・`predict_rt` を足す |
| `lipidmix/curation/evidence.py` | 変更 | `collect(keep_measured=True)` と公開の `arf_rows_by_spot` |
| `lipidmix/curation/candidates.py` | 新規 | ①② の候補の採点・統合・制約・並び順（純ロジック＋store 呼び出し） |
| `lipidmix/curation/relations.py` | 新規 | ④ の関係表・同位体比・相関・ペア判定（純関数） |
| `lipidmix/curation/suggest.py` | 新規 | 対象の選定・組み立て・保存・読み込み・要約 TSV・送信内容の展開 |
| `lipidmix/curation/viewer_common.js` | 新規 | 両ビューア共通の canvas 描画（既存 `viewer.html` から移す） |
| `lipidmix/curation/viewer.html` | 変更 | 共通 JS を差し込む印に置き換える。`decision` の表示 |
| `lipidmix/curation/suggest_viewer.html` | 新規 | 候補付けのビューア |
| `lipidmix/curation/viewer.py` | 変更 | 共通 JS の差し込み、`render_suggest_html` |
| `lipidmix/curation/apply.py` | 変更 | `assign` による同定の置き換え、`redundant` の除外、メタ行 |
| `lipidmix/tools/curation_tools.py` | 変更 | `curation_suggest`、`curation_submit` / `curation_flags` の拡張 |
| `lipidmix/arf/tools.py` | 変更 | `arf_export_differential` に置き換えを適用 |
| `lipidmix/analysis/dataset_export.py` | 変更 | mzTab-M 経路に置き換えと除外を適用 |
| `lipidmix/arf2/tools.py` | 変更 | `arf2_annotate_identities` の `curation_flag` に `assign` / `redundant` |
| `scripts/verify_curation_suggest.py` | 新規 | 実データでの正解ありの検証（② の順位）と ④ の再現率 |
| `tests/curation_fixtures.py` | 変更 | `arf2_spot_raw(peak_links=…)`、`write_suggest_set` |
| `tests/test_arf2_ion_features.py` ほか | 新規 | 各タスクのテスト（タスク内に明記） |
| `USAGE.md` `CLAUDE.md` `docs/workflow/{index,curation}.md` `docs/output_format/curation.md` | 変更 | 腐敗防止テストが縛る文書 |

---

### Task 1: 実装前の事実確認と `.arf2` Key 10 の復号

**Files:**
- Create: `docs/schema/IonFeatureCharacter.md`
- Create: `lipidmix/arf2/ion_features.py`
- Modify: `tests/curation_fixtures.py`（`arf2_spot_raw` に `peak_links`）
- Test: `tests/test_arf2_ion_features.py`

**Interfaces:**
- Produces: `LINK_KINDS: dict[int, str]`、`decode_ion_features(raw_spot: list) -> dict`（`{"links": [{"spot_id": int, "kind": str}], "isotope_weight": int|None, "isotope_parent": int|None, "peak_group": int|None}`）、`load_ion_features(arf2_path) -> dict[int, dict]`（キーは MasterAlignmentID）。fixture `arf2_spot_raw(..., peak_links: list[tuple[int, int]] | None = None)`（`(linked_spot_id, character_int)`）。

- [ ] **Step 1: 上流の定義を確かめて schema 表を書く**

`C:\Users\yuu18\source\repos\MsdialWorkbench` で次を読み、Key 番号を確認する（2026-09-29 時点の HEAD は `afd5f9522`）。

```bash
grep -n "Key(" /c/Users/yuu18/source/repos/MsdialWorkbench/src/MSDIAL5/MsdialCore/DataObj/IonFeatureCharacter.cs
grep -n "Key(" -A 1 /c/Users/yuu18/source/repos/MsdialWorkbench/src/MSDIAL5/MsdialCore/DataObj/ChromatogramPeakFeature.cs | sed -n '/LinkedPeakFeature/,$p' | head
grep -n "enum PeakLinkFeatureEnum" -A 2 /c/Users/yuu18/source/repos/MsdialWorkbench/src/Common/CommonStandard/Enum/CommonEnums.cs
```

`docs/schema/IonFeatureCharacter.md` を書く（既存 `docs/schema/MsScanMatchResult.md` の体裁に合わせる）:

```markdown
# IonFeatureCharacter（`AlignmentSpotProperty` Key 10）

上流: `MsdialWorkbench` `afd5f9522`（2026-09-08）
`src/MSDIAL5/MsdialCore/DataObj/IonFeatureCharacter.cs`、
`LinkedPeakFeature` は `src/MSDIAL5/MsdialCore/DataObj/ChromatogramPeakFeature.cs`、
`PeakLinkFeatureEnum` は `src/Common/CommonStandard/Enum/CommonEnums.cs`。

| Key | 型 | 名前 | 備考 |
|---|---|---|---|
| 0 | AdductIon | AdductType | 本システムは `.arf2` Key 54 の AdductType を使う |
| 1 | AdductIon | AdductTypeByAmalgamationProgram | |
| 2 | int | Charge | 既定 1 |
| 3 | List<LinkedPeakFeature> | PeakLinks | 下表 |
| 4 | int | IsotopeWeightNumber | 0 = 単同位体。-1 既定 |
| 5 | int | IsotopeParentPeakID | |
| 6 | int | PeakGroupID | |
| 7 | bool | IsLinked | |
| 8 | int | AdductParent | |

LinkedPeakFeature: Key 0 `LinkedPeakID`（int）、Key 1 `Character`（PeakLinkFeatureEnum）。

PeakLinkFeatureEnum: 0 SameFeature, 1 Isotope, 2 Adduct, 3 ChromSimilar, 4 FoundInUpperMsMs, 5 CorrelSimilar。

## 実データでの確認（2026-09-29、kidney neg/pos）

- IsotopeWeightNumber はアラインメントの全スポットで 0（同位体は別スポットとして残らない）。
- LinkedPeakID は <Step 2 の結果を書く> を指す。
- FoundInUpperMsMs の向き: <Step 2 の結果を書く>。
```

- [ ] **Step 2: `LinkedPeakID` の意味と FoundInUpperMsMs の向きを実データで確かめる**

スクラッチに次を書いて実行する（scratchpad か `check.py`。コミットしない）。

```python
# check.py
import glob, statistics
from lipidmix.arf2.reader import load_raw_spots
D = r"C:\Users\yuu18\datasets\a_lipidome_landscape_of_aging_in_mice\rplc\kidney\neg"
raw = load_raw_spots(sorted(glob.glob(D + r"\*.arf2"))[-1])
by_id = {int(r[0]): r for r in raw}
rt = lambda r: r[4][0][1][0]
drts, mz_higher = [], {"link_target_higher": 0, "link_target_lower": 0}
for r in raw:
    for link in (r[10][3] or []):
        other = by_id.get(int(link[0]))
        if other is None:
            continue
        drts.append(abs(rt(r) - rt(other)))
        if link[1] == 4:
            mz_higher["link_target_higher" if other[5] > r[5] else "link_target_lower"] += 1
print("n", len(drts), "median |dRT|", statistics.median(drts), "p95", sorted(drts)[int(len(drts) * .95)])
print(mz_higher)
```

Run: `PYTHONPATH=. C:/Python314/python.exe check.py`
期待: |ΔRT| の中央値が 0.1 分未満なら `LinkedPeakID` は MasterAlignmentID を指す。`link_target_higher` が多ければ「断片側のスポットが上位の precursor へリンクを張る」。結果を Step 1 の表の 2 行に書き、`check.py` の中身を消す。中央値が 0.1 分以上なら止めてユーザーに報告する（前提が崩れる）。

- [ ] **Step 3: 残りの事実を確かめて記録する**（spec §10 の 2・3・5）

1. `library_load` で `Dataset_*_Loaded.msp2.dbs` を開き、kidney neg の注釈付きスポットの ① の下位候補（代表以外）のうち `store.record_by_scan_id` で引ける割合を数える（`check.py` で、`lipidmix.arf2.match_results.load_spot_annotations` と `library.store.open_store` を使う）。
2. param ファイルの行 `Searched adduct ions:` と `Retention time tolerance for alignment:` と `MS1 tolerance for centroid:` の書式を目で確認する（Task 2 のパーサの入力）。
3. massbank-context（`C:\Users\yuu18\` 直下の兄弟リポジトリ。場所は `ls /c/Users/yuu18 | grep -i massbank` で探す）の差次的エクスポートの読み取り（`load_differential` と `_parse_meta_line`）を読み、`name_source` の値を列挙して検査しているか、未知のメタ行キーを落とさず保持するかを確かめる。値を検査していたら止めてユーザーに報告する。

結果は main ツリーの `docs/HISTRY.md` に日付見出しで追記する。

- [ ] **Step 4: 失敗するテストを書く**

`tests/curation_fixtures.py` の `arf2_spot_raw` に引数 `peak_links` を足す前提でテストを書く。

```python
# tests/test_arf2_ion_features.py
from lipidmix.arf2.ion_features import LINK_KINDS, decode_ion_features, load_ion_features
from tests.curation_fixtures import arf2_spot_raw, write_arf2


def test_decode_links_and_isotope_fields():
    raw = arf2_spot_raw(spot_id=3, peak_links=[(0, 4), (1, 5), (2, 9)])
    decoded = decode_ion_features(raw)
    assert decoded["links"] == [{"spot_id": 0, "kind": "found_in_upper_msms"},
                                {"spot_id": 1, "kind": "correl_similar"},
                                {"spot_id": 2, "kind": "unknown_9"}]
    assert decoded["isotope_weight"] == 0
    assert decoded["peak_group"] == -1


def test_decode_tolerates_missing_key_10():
    raw = arf2_spot_raw(spot_id=0)
    raw[10] = None
    assert decode_ion_features(raw) == {"links": [], "isotope_weight": None,
                                        "isotope_parent": None, "peak_group": None}


def test_load_ion_features_is_keyed_by_master_alignment_id(tmp_path):
    path = write_arf2(tmp_path / "a.arf2", [arf2_spot_raw(spot_id=0),
                                            arf2_spot_raw(spot_id=1, peak_links=[(0, 3)])])
    features = load_ion_features(path)
    assert features[1]["links"] == [{"spot_id": 0, "kind": "chrom_similar"}]
    assert features[0]["links"] == []


def test_link_kinds_match_the_upstream_enum():
    assert LINK_KINDS == {0: "same_feature", 1: "isotope", 2: "adduct", 3: "chrom_similar",
                          4: "found_in_upper_msms", 5: "correl_similar"}
```

- [ ] **Step 5: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_arf2_ion_features.py -q`
Expected: FAIL（`ModuleNotFoundError: lipidmix.arf2.ion_features`）

- [ ] **Step 6: fixture と実装を書く**

`tests/curation_fixtures.py` の `arf2_spot_raw` のシグネチャに `peak_links: list | None = None` を足し、`values` に次を足す（既存の呼び出しは引数なしなので空リストになる）:

```python
        10: [[-1.00782503207, 1, adduct, 1, 1 if ion_mode == 1 else 0, True, 0.0, 0.0, False, True],
             [0.0, 0, "", 0, 0, False, 0.0, 0.0, False, False], 1,
             [[int(i), int(c)] for i, c in (peak_links or [])], 0, -1, -1, bool(peak_links), -1],
```

```python
# lipidmix/arf2/ion_features.py
"""`.arf2` Key 10（IonFeatureCharacter）の復号。Key 番号の正準は docs/schema/IonFeatureCharacter.md。

MS-DIAL がアラインメント時に計算したスポット間の関係（同位体・アダクト・クロマトグラムの類似・
上位イオンの MS/MS に断片として出ている・試料間の相関）を、候補付けの ④ の裏付けに使う。
deps: arf2.reader だけ。
"""
from __future__ import annotations

from lipidmix.arf2.reader import load_raw_spots

LINK_KINDS = {0: "same_feature", 1: "isotope", 2: "adduct", 3: "chrom_similar",
              4: "found_in_upper_msms", 5: "correl_similar"}
_EMPTY = {"links": [], "isotope_weight": None, "isotope_parent": None, "peak_group": None}


def _int_at(values: list, key: int):
    value = values[key] if len(values) > key else None
    return int(value) if isinstance(value, int) and not isinstance(value, bool) else None


def decode_ion_features(raw_spot: list) -> dict:
    character = raw_spot[10] if len(raw_spot) > 10 else None
    if not isinstance(character, list):
        return dict(_EMPTY, links=[])
    links = []
    raw_links = character[3] if len(character) > 3 and isinstance(character[3], list) else []
    for link in raw_links:
        if isinstance(link, list) and len(link) >= 2 and isinstance(link[0], int) \
                and isinstance(link[1], int):
            links.append({"spot_id": int(link[0]),
                          "kind": LINK_KINDS.get(int(link[1]), f"unknown_{int(link[1])}")})
    return {"links": links, "isotope_weight": _int_at(character, 4),
            "isotope_parent": _int_at(character, 5), "peak_group": _int_at(character, 6)}


def load_ion_features(arf2_path) -> dict[int, dict]:
    return {int(raw[0]): decode_ion_features(raw) for raw in load_raw_spots(arf2_path)}
```

- [ ] **Step 7: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_arf2_ion_features.py tests/test_curation_evidence.py tests/test_curation_review.py -q`
Expected: PASS（fixture の変更で既存テストが壊れていないこと）

- [ ] **Step 8: コミット**

```bash
git add docs/schema/IonFeatureCharacter.md lipidmix/arf2/ion_features.py tests/curation_fixtures.py tests/test_arf2_ion_features.py
git commit -m "feat(arf2): IonFeatureCharacter（Key 10）のリンクを復号する"
```

---

### Task 2: アダクトの解析と param ファイルの読み取り

**Files:**
- Create: `lipidmix/msdial/adducts.py`
- Create: `lipidmix/msdial/analysis_params.py`
- Test: `tests/test_msdial_adducts.py`、`tests/test_msdial_analysis_params.py`

**Interfaces:**
- Consumes: `lipidmix.msdial.peak_verification.parse_formula`、`ELEMENT_MASSES`、`ELECTRON_MASS`、`ADDUCT_SHIFTS`（テストの突き合わせ用）。
- Produces:
  - `Adduct`（frozen dataclass: `name: str, n_mol: int, shift: float, charge: int, polarity: str`）
  - `parse_adduct(name: str | None) -> Adduct | None`
  - `mz_from_neutral(neutral: float, adduct: Adduct) -> float`、`neutral_from_mz(mz: float, adduct: Adduct) -> float`
  - `DEFAULT_ADDUCTS: dict[str, list[str]]`（キー `"Positive"` / `"Negative"`）
  - `find_param_file(arf2_path) -> Path | None`、`read_analysis_params(path) -> dict`（`{"ion_mode": str|None, "searched_adducts": list[str], "rt_tolerance_alignment": float|None, "ms1_tolerance": float|None}`）、`resolve_analysis_params(arf2_path, ion_mode: str | None) -> dict`（`{"searched_adducts": [...], "rt_window": float, "source": "param_file"|"default", "path": str|None}`）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_msdial_adducts.py
import pytest

from lipidmix.msdial.adducts import mz_from_neutral, neutral_from_mz, parse_adduct
from lipidmix.msdial.peak_verification import ADDUCT_SHIFTS


@pytest.mark.parametrize("name", sorted(ADDUCT_SHIFTS))
def test_parse_agrees_with_the_existing_shift_table(name):
    sign, shift, charge, n_mol = ADDUCT_SHIFTS[name]
    adduct = parse_adduct(name)
    assert adduct is not None
    assert adduct.polarity == sign and adduct.charge == charge and adduct.n_mol == n_mol
    assert adduct.shift == pytest.approx(shift, abs=2e-5)


@pytest.mark.parametrize("name, n_mol, charge, polarity", [
    ("[M+C2H3N+Na-2H]-", 1, 1, "-"), ("[M+CH3COONa-H]-", 1, 1, "-"), ("[2M+FA-H]-", 2, 1, "-"),
    ("[2M+Hac-H]-", 2, 1, "-"), ("[3M-H]-", 3, 1, "-"), ("[M-3H]3-", 1, 3, "-"),
    ("[M+TFA-H]-", 1, 1, "-"), ("[M-C6H10O5-H]-", 1, 1, "-"), ("[M+2H]2+", 1, 2, "+"),
    ("[M+H-H2O]+", 1, 1, "+"), ("[M+K]+", 1, 1, "+"), ("[M+2NH4]2+", 1, 2, "+"),
])
def test_parse_msdial_searched_adducts(name, n_mol, charge, polarity):
    adduct = parse_adduct(name)
    assert (adduct.n_mol, adduct.charge, adduct.polarity) == (n_mol, charge, polarity)


def test_formic_acid_alias_equals_formate():
    # [M+FA-H]- と [M+HCOO]- は同じイオン
    assert parse_adduct("[M+FA-H]-").shift == pytest.approx(parse_adduct("[M+HCOO]-").shift, abs=1e-6)


@pytest.mark.parametrize("name", [None, "", "Unknown", "[M+Xx]+", "M+H", "[M+H]"])
def test_unparseable_is_none(name):
    assert parse_adduct(name) is None


def test_round_trip_between_mz_and_neutral():
    adduct = parse_adduct("[2M+HCOO]-")
    assert neutral_from_mz(mz_from_neutral(759.5778, adduct), adduct) == pytest.approx(759.5778)
```

```python
# tests/test_msdial_analysis_params.py
from lipidmix.msdial.analysis_params import (
    DEFAULT_ADDUCTS, find_param_file, read_analysis_params, resolve_analysis_params)

PARAM = """# Project information
Ion mode: Negative
Searched adduct ions: [M-H]-,[M+HCOO]-,[M+CH3COO]-,[2M-H]-
MS1 tolerance for centroid: 0.01
Retention time tolerance for alignment: 0.1
"""


def test_read_params(tmp_path):
    path = tmp_path / "Dataset_2026_09_09_17_28_59_param_202609091800.txt"
    path.write_text(PARAM, encoding="utf-8")
    assert read_analysis_params(path) == {
        "ion_mode": "Negative", "searched_adducts": ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[2M-H]-"],
        "rt_tolerance_alignment": 0.1, "ms1_tolerance": 0.01}


def test_find_param_file_picks_the_newest_name(tmp_path):
    (tmp_path / "Dataset_a_param_202601010000.txt").write_text(PARAM, encoding="utf-8")
    (tmp_path / "Dataset_a_param_202609091800.txt").write_text(PARAM, encoding="utf-8")
    arf2 = tmp_path / "AlignmentResult_x.arf2"
    assert find_param_file(arf2).name == "Dataset_a_param_202609091800.txt"


def test_resolve_falls_back_to_defaults_without_a_param_file(tmp_path):
    resolved = resolve_analysis_params(tmp_path / "AlignmentResult_x.arf2", "Negative")
    assert resolved == {"searched_adducts": DEFAULT_ADDUCTS["Negative"], "rt_window": 0.1,
                        "source": "default", "path": None}


def test_resolve_uses_the_param_file(tmp_path):
    (tmp_path / "Dataset_a_param_1.txt").write_text(PARAM.replace("0.1\n", "0.08\n"), encoding="utf-8")
    resolved = resolve_analysis_params(tmp_path / "AlignmentResult_x.arf2", "Negative")
    assert resolved["source"] == "param_file" and resolved["rt_window"] == 0.08
    assert resolved["searched_adducts"][0] == "[M-H]-"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_msdial_adducts.py tests/test_msdial_analysis_params.py -q`
Expected: FAIL（モジュールが無い）

- [ ] **Step 3: 実装を書く**

```python
# lipidmix/msdial/adducts.py
"""MS-DIAL のアダクト表記（`[2M+FA-H]-` など）の解析と、m/z ⇄ 中性質量の換算（純関数）。

m/z = (n_mol·M + shift) / charge。shift は付加・脱離する原子団の質量の符号付き和から、
正イオンなら電子 charge 個ぶんを引き、負イオンなら足したもの（既存 `ADDUCT_SHIFTS` と同じ定義）。
略記: FA = ギ酸 CH2O2、Hac = 酢酸 C2H4O2、TFA = トリフルオロ酢酸 C2HF3O2、ACN = C2H3N。
deps: msdial.peak_verification（元素質量と組成式の解析）。
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from lipidmix.msdial.peak_verification import ELECTRON_MASS, ELEMENT_MASSES, parse_formula

_ALIASES = {"FA": "CH2O2", "Hac": "C2H4O2", "TFA": "C2HF3O2", "ACN": "C2H3N"}
_ADDUCT_RE = re.compile(r"^\[(\d*)M((?:[+-]\d*[A-Z][A-Za-z0-9]*)*)\](\d*)([+-])$")
_TERM_RE = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)")


@dataclass(frozen=True)
class Adduct:
    name: str
    n_mol: int
    shift: float
    charge: int
    polarity: str


def _group_mass(token: str) -> float:
    counts = parse_formula(_ALIASES.get(token, token))
    return sum(ELEMENT_MASSES[element] * n for element, n in counts.items())


def parse_adduct(name) -> Adduct | None:
    if not name:
        return None
    match = _ADDUCT_RE.match(str(name).strip())
    if not match:
        return None
    n_mol = int(match.group(1) or 1)
    charge = int(match.group(3) or 1)
    polarity = match.group(4)
    shift = 0.0
    try:
        for sign, count, token in _TERM_RE.findall(match.group(2)):
            shift += (1 if sign == "+" else -1) * int(count or 1) * _group_mass(token)
    except (KeyError, ValueError):
        return None
    shift += (-1 if polarity == "+" else 1) * charge * ELECTRON_MASS
    return Adduct(name=str(name).strip(), n_mol=n_mol, shift=shift, charge=charge, polarity=polarity)


def mz_from_neutral(neutral: float, adduct: Adduct) -> float:
    return (adduct.n_mol * neutral + adduct.shift) / adduct.charge


def neutral_from_mz(mz: float, adduct: Adduct) -> float:
    return (mz * adduct.charge - adduct.shift) / adduct.n_mol
```

`ELEMENT_MASSES` に無い元素（例: `Li`）を含むアダクトは `KeyError` で `None` になる。これは意図どおりで、関係表に入らないだけ。

```python
# lipidmix/msdial/analysis_params.py
"""MS-DIAL の解析 param ファイル（`<Dataset>_param_<ts>.txt`）から、候補付けに要る値だけを読む。

読む行: `Ion mode:`、`Searched adduct ions:`（その解析の極性の一覧）、`MS1 tolerance for centroid:`、
`Retention time tolerance for alignment:`。ファイルが無ければ極性ごとの既定値に落とし、どちらを
使ったかを `source` で返す。deps: stdlib だけ。
"""
from __future__ import annotations

from pathlib import Path

DEFAULT_ADDUCTS = {
    "Negative": ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Cl]-", "[M-H2O-H]-", "[2M-H]-", "[M-2H]2-"],
    "Positive": ["[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+", "[2M+H]+", "[M+2H]2+"],
}
DEFAULT_RT_WINDOW = 0.1
_KEYS = {"Ion mode": "ion_mode", "Searched adduct ions": "searched_adducts",
         "MS1 tolerance for centroid": "ms1_tolerance",
         "Retention time tolerance for alignment": "rt_tolerance_alignment"}


def find_param_file(arf2_path) -> Path | None:
    folder = Path(arf2_path).parent
    files = sorted(folder.glob("*_param_*.txt"), key=lambda p: p.name)
    return files[-1] if files else None


def _float(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def read_analysis_params(path) -> dict:
    found = {"ion_mode": None, "searched_adducts": [], "rt_tolerance_alignment": None,
             "ms1_tolerance": None}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        key, sep, value = line.partition(":")
        name = _KEYS.get(key.strip())
        if not sep or name is None:
            continue
        value = value.strip()
        if name == "searched_adducts":
            found[name] = [v.strip() for v in value.split(",") if v.strip()]
        elif name == "ion_mode":
            found[name] = value or None
        else:
            found[name] = _float(value)
    return found


def resolve_analysis_params(arf2_path, ion_mode: str | None) -> dict:
    path = find_param_file(arf2_path)
    params = read_analysis_params(path) if path else None
    mode = (params or {}).get("ion_mode") or ion_mode or "Negative"
    adducts = (params or {}).get("searched_adducts") or DEFAULT_ADDUCTS.get(mode, DEFAULT_ADDUCTS["Negative"])
    rt_window = (params or {}).get("rt_tolerance_alignment") or DEFAULT_RT_WINDOW
    return {"searched_adducts": list(adducts), "rt_window": rt_window,
            "source": "param_file" if params else "default", "path": str(path) if path else None}
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_msdial_adducts.py tests/test_msdial_analysis_params.py -q`
Expected: PASS。`ADDUCT_SHIFTS` との突き合わせが 2e-5 を超えて外れたら、`ADDUCT_SHIFTS` 側の定数ではなく `_group_mass` を疑う（既存表は実データで検証済み）。

- [ ] **Step 5: コミット**

```bash
git add lipidmix/msdial/adducts.py lipidmix/msdial/analysis_params.py tests/test_msdial_adducts.py tests/test_msdial_analysis_params.py
git commit -m "feat(msdial): アダクト表記の解析と param ファイルの読み取りを足す"
```

---

### Task 3: 判断の記録（`assign` / `redundant`）

**Files:**
- Modify: `lipidmix/curation/flags.py`
- Modify: `lipidmix/curation/review.py`（`run_review` で `decision` を分ける）
- Modify: `lipidmix/curation/viewer.html`（カードに記録済みの判断を 1 行出す）
- Test: `tests/test_curation_flags.py`、`tests/test_curation_review.py`

**Interfaces:**
- Produces:
  - `FLAG_VALUES = ("wrong", "suspect", "clear", "assign", "redundant")`、`REVIEW_FLAG_VALUES = ("wrong", "suspect", "clear")`、`SUGGEST_FLAG_VALUES = ("assign", "redundant", "clear")`
  - `validate_entries`（既存。`REVIEW_FLAG_VALUES` だけを受ける。メッセージは現行と同じ文面で値の列挙だけ変わらない）
  - `flags_digest(effective)`: wrong / suspect は現行と同じ `(spot, flag)`、assign は `(spot, flag, name, level)`、redundant は `(spot, flag, of, relation)`
  - `split_decisions(effective: dict[int, dict]) -> dict`（`{"wrong": set, "suspect": set, "assign": {spot: row}, "redundant": {spot: row}}`）
  - レビューの各スポットの `ev["flag"]` は wrong / suspect のときだけ値を持ち、`ev["decision"]` は `{"flag": "assign", "name": ...}` か `{"flag": "redundant", "of": ...}` か `None`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_curation_flags.py` の末尾に足す:

```python
from lipidmix.curation.flags import (
    FlagStore, effective_flags, flags_digest, split_decisions, validate_entries)

_ALIGN = {"alignment_file": "a.arf2", "alignment_sha256": "s1"}


def test_assign_after_wrong_wins_and_clear_removes_it(tmp_path):
    store = FlagStore(tmp_path)
    store.append([{"spot_id": 4, "flag": "wrong", "note": ""}], alignment=_ALIGN, review_id="cr-1", source="user")
    store.append([{"spot_id": 4, "flag": "assign", "name": "PC 34:1", "level": "sum", "note": ""}],
                 alignment=_ALIGN, review_id="cs-1", source="user")
    effective = store.effective("s1")
    assert effective[4]["flag"] == "assign"
    assert split_decisions(effective)["assign"][4]["name"] == "PC 34:1"
    assert split_decisions(effective)["wrong"] == set()
    store.append([{"spot_id": 4, "flag": "clear", "note": ""}], alignment=_ALIGN, review_id="cs-1", source="user")
    assert store.effective("s1") == {}


def test_rows_accept_the_new_flag_values(tmp_path):
    store = FlagStore(tmp_path)
    store.append([{"spot_id": 1, "flag": "redundant", "of": 0, "relation": "isotope_M+2", "note": ""}],
                 alignment=_ALIGN, review_id="cs-1", source="user")
    assert store.rows()[0]["flag"] == "redundant"


def test_digest_for_wrong_and_suspect_is_unchanged():
    # 既存ユーザーの curation_flags_sha256 を変えない（メタ行が現行と同じであること）
    import hashlib, json
    effective = {1: {"flag": "wrong"}, 2: {"flag": "suspect"}}
    legacy = hashlib.sha256(json.dumps(sorted([(1, "wrong"), (2, "suspect")])).encode()).hexdigest()
    assert flags_digest(effective) == legacy


def test_digest_changes_when_the_assigned_name_changes():
    a = {1: {"flag": "assign", "name": "PC 34:1", "level": "sum"}}
    b = {1: {"flag": "assign", "name": "PC 34:2", "level": "sum"}}
    assert flags_digest(a) != flags_digest(b)


def test_review_submissions_still_reject_assign():
    with pytest.raises(ValueError, match="flag"):
        validate_entries([{"spot_id": 1, "flag": "assign"}], allowed_spot_ids={1})
```

`tests/test_curation_review.py` に足す（既存のヘルパ `_review(paths)` と fixture を使う）:

```python
def test_assign_and_redundant_are_decisions_not_flags(dataset_paths):
    paths = dataset_paths
    first = _review(paths)
    store = flags.FlagStore(flags.curation_dir(paths["arf2"]))
    store.append([{"spot_id": 0, "flag": "assign", "name": "PC 34:1", "level": "sum", "note": ""},
                  {"spot_id": 1, "flag": "redundant", "of": 0, "relation": "isotope_M+1", "note": ""}],
                 alignment=first["alignment"], review_id="cs-x", source="user")
    again = {s["spot_id"]: s for s in _review(paths)["spots"]}
    assert again[0]["flag"] is None and again[0]["decision"] == {"flag": "assign", "name": "PC 34:1"}
    assert again[1]["flag"] is None and again[1]["decision"] == {"flag": "redundant", "of": 0}
```

（`dataset_paths` はこのファイルの既存 fixture 名に合わせる。既存は `write_alignment_set` を呼ぶ fixture なので、その名前を使う。）

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_flags.py tests/test_curation_review.py -q`
Expected: FAIL（`split_decisions` が無い、`redundant` の行で `FlagFileError`）

- [ ] **Step 3: 実装を書く**

`lipidmix/curation/flags.py`:

```python
FLAG_VALUES = ("wrong", "suspect", "clear", "assign", "redundant")
REVIEW_FLAG_VALUES = ("wrong", "suspect", "clear")
SUGGEST_FLAG_VALUES = ("assign", "redundant", "clear")
```

`validate_entries` の `if flag not in FLAG_VALUES:` を `if flag not in REVIEW_FLAG_VALUES:` にし、メッセージの `{FLAG_VALUES}` も `{REVIEW_FLAG_VALUES}` にする。モジュール docstring の 2 行目を「記録するのはユーザーの判断: レビューの wrong / suspect、候補付けの assign / redundant、取り消しの clear。」に直す。

`flags_digest` を置き換える:

```python
def _digest_item(spot: int, row: dict) -> list:
    flag = row["flag"]
    if flag == "assign":
        return [spot, flag, row.get("name"), row.get("level")]
    if flag == "redundant":
        return [spot, flag, row.get("of"), row.get("relation")]
    return [spot, flag]


def flags_digest(effective: dict[int, dict]) -> str:
    # wrong / suspect は以前と同じ [spot, flag]（既存の curation_flags_sha256 を変えない）
    canonical = json.dumps(sorted((_digest_item(spot, row) for spot, row in effective.items()),
                                  key=lambda item: item[0]))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def split_decisions(effective: dict[int, dict]) -> dict:
    return {"wrong": {s for s, r in effective.items() if r["flag"] == "wrong"},
            "suspect": {s for s, r in effective.items() if r["flag"] == "suspect"},
            "assign": {s: r for s, r in effective.items() if r["flag"] == "assign"},
            "redundant": {s: r for s, r in effective.items() if r["flag"] == "redundant"}}
```

注意: 以前は `sorted((spot, flag) ...)` のタプルを `json.dumps` しており、JSON では `[[1, "wrong"], ...]` になる。上の実装も `[[1, "wrong"], ...]` を出すので、テスト `test_digest_for_wrong_and_suspect_is_unchanged` で一致を確かめる。

`lipidmix/curation/review.py` の `run_review` のループで、`ev["flag"]` と `ev["flag_note"]` を設定している 2 行を置き換える:

```python
        row = existing.get(ev["spot_id"]) or {}
        is_flag = row.get("flag") in ("wrong", "suspect")
        ev["flag"] = row.get("flag") if is_flag else None
        ev["flag_note"] = row.get("note") if is_flag else None
        ev["decision"] = ({"flag": "assign", "name": row.get("name")} if row.get("flag") == "assign"
                          else {"flag": "redundant", "of": row.get("of")} if row.get("flag") == "redundant"
                          else None)
```

`lipidmix/curation/viewer.html` の `spotCard` で、理由コードの行の直後に足す:

```javascript
  if (spot.decision) card.appendChild(el("div", {class: "meta"}, spot.decision.flag === "assign"
    ? `記録済み: 候補を採用（${spot.decision.name}）` : `記録済み: #${spot.decision.of} の別イオン`));
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_flags.py tests/test_curation_review.py tests/test_curation_tools.py tests/test_curation_apply.py -q`
Expected: PASS

- [ ] **Step 5: コミット**

```bash
git add lipidmix/curation/flags.py lipidmix/curation/review.py lipidmix/curation/viewer.html tests/test_curation_flags.py tests/test_curation_review.py
git commit -m "feat(curation): 判断の記録に assign / redundant を足す"
```

---

### Task 4: ①② の候補（下位候補・再検索・制約・並び順）

**Files:**
- Modify: `lipidmix/arf2/match_results.py`（`usable_candidates`、`load_spot_candidates`）
- Modify: `lipidmix/curation/trend.py`（`sum_composition`、`predict_rt`）
- Create: `lipidmix/curation/candidates.py`
- Test: `tests/test_curation_candidates.py`、`tests/test_curation_trend.py`、`tests/test_arf2_match_results.py`（既存が無ければ新規）

**Interfaces:**
- Consumes: `analysis.spectral_match.match_spectrum` / `total_score`、`library.store.LibraryStore.candidates` / `record_by_scan_id`、`library.store.library_id_from_annotator`、`msdial.adducts.parse_adduct`、`msdial.peak_verification.adduct_consistency`、`plots.mirror.build_mirror_payload`、`curation.evidence._keep_most_intense_peaks` / `MIRROR_MAX_PEAKS`、`trend.composition`。
- Produces:
  - `match_results.usable_candidates(container) -> list[dict]`（非 decoy・Source≠Unknown）、`load_spot_candidates(path) -> dict[int, list[dict]]`
  - `trend.sum_composition(name) -> str | None`（`PC 16:0_18:1` → `PC 34:1`）
  - `trend.predict_rt(classes: dict, ontology: str, carbon: int, db: int) -> dict | None`（`{"predicted_rt": float, "scale": float}`）
  - `candidates.record_key(record) -> tuple`
  - `candidates.build_library_candidates(ev: dict, *, msdial_matches: list[dict], representative: dict | None, store, scoring: dict, trends: dict, th: dict, exclude_current: bool, top_n: int) -> dict`。戻り値 `{"candidates": [...], "n_hard_removed": int, "measured": [[mz, int], ...]}`。`ev` は `evidence.collect(keep_measured=True)` の 1 要素（`_measured` を持つ）。`scoring` は `{"mz_tol", "ms2_tol", "mass_begin", "mass_end", "relative_amp_cutoff", "absolute_amp_cutoff", "use_rt"}`。
  - 候補 dict のキー: `candidate_id`（`"L1"`, `"L2"`, …、並べた後に振る）、`kind="library"`、`source`、`name`、`sum_name`、`ontology`、`adduct`、`formula`、`inchikey`、`precursor_mz`、`ref_rt`、`library_id`、`record_index`、`dmz_mda`、`ppm`、`scores`（`total_score` `weighted_dot_product` `simple_dot_product` `reverse_dot_product` `matched_peaks_percentage` `matched_peaks_count`）、`msdial_total_score`、`trend`（`{"predicted_rt","residual","z"}` か `None`）、`hard`、`soft`、`info`、`mirror`（`{"reference", "matched_mz", "matched_measured_mz"}` か `None`）。

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_curation_trend.py に足す
from lipidmix.curation.trend import predict_rt, sum_composition


def test_sum_composition_of_species_and_msdial_names():
    assert sum_composition("PC 16:0_18:1") == "PC 34:1"
    assert sum_composition("PC 34:1|PC 16:0_18:1") == "PC 34:1"
    assert sum_composition("Cer 18:1;O2/24:0") == "Cer 42:1;O2"
    assert sum_composition("RIKEN N-VS1 ID-45 from x") is None


def test_predict_rt_uses_the_additive_model():
    classes = {"PC": {"coef": [1.0, 0.5, -0.3], "scale": 0.05},
               "PE": {"coef": [2.0, 0.4], "scale": 0.1}}          # DB が 1 種類しかなかったクラス
    assert predict_rt(classes, "PC", 34, 1) == {"predicted_rt": 1.0 + 17.0 - 0.3, "scale": 0.05}
    assert predict_rt(classes, "PE", 36, 2) == {"predicted_rt": 2.0 + 14.4, "scale": 0.1}
    assert predict_rt(classes, "TG", 50, 1) is None
```

```python
# tests/test_arf2_match_results.py（既存ファイルがあれば末尾に足す）
from lipidmix.arf2.match_results import usable_candidates
from tests.curation_fixtures import match_result


def test_usable_candidates_drop_decoys_and_unknown_source():
    container = [[match_result({0: "A", 31: 2}), match_result({0: "B", 30: True}),
                  match_result({0: "C", 26: 1})], {}, []]
    assert [c["name"] for c in usable_candidates(container)] == ["A"]
```

```python
# tests/test_curation_candidates.py
import pytest

from lipidmix.curation import candidates, judge
from lipidmix.library import store as library_store

TH = judge.resolve_thresholds(None)
SCORING = {"mz_tol": 0.01, "ms2_tol": 0.025, "mass_begin": 0.0, "mass_end": 2000.0,
           "relative_amp_cutoff": 0.0, "absolute_amp_cutoff": 0.0, "use_rt": False}

MSP = """NAME: PC 16:0_18:1
PRECURSORMZ: 804.5760
PRECURSORTYPE: [M+HCOO]-
IONMODE: Negative
INCHIKEY: KEY-PC341
FORMULA: C42H82NO8P
Num Peaks: 2
255.23 999
281.25 800

NAME: PE 18:0_18:2
PRECURSORMZ: 804.5750
PRECURSORTYPE: [M+H]+
INCHIKEY: KEY-PE-POS
Num Peaks: 1
184.07 999

NAME: PS 16:0_18:1
PRECURSORMZ: 804.5700
PRECURSORTYPE: [M-H]-
IONMODE: Negative
INCHIKEY: KEY-PS
Num Peaks: 1
100.00 999
"""


@pytest.fixture()
def store(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    path = tmp_path / "lib.msp"
    path.write_text(MSP, encoding="utf-8")
    s = library_store.open_store(path)
    yield s
    s.close()


def _ev(measured, *, mz=804.5762, ion_mode="Negative"):
    return {"spot_id": 7, "rep_mz": mz, "mz": mz, "rep_rt": 12.0, "rt": 12.0, "ion_mode": ion_mode,
            "_measured": measured}


def test_research_scores_all_in_window_and_filters_polarity(store):
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0], [281.25, 700.0]]), msdial_matches=[], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    names = [c["name"] for c in out["candidates"]]
    # PE [M+H]+ は IONMODE 欄が無い（極性不明）ので store の絞り込みを抜けて検索に掛かり、
    # アダクトの極性が負イオンのスポットと矛盾する（ハード）ので消える。窓は ±0.01 で 3 件とも掛かる
    assert names == ["PC 16:0_18:1", "PS 16:0_18:1"]
    assert out["n_hard_removed"] == 1
    top = out["candidates"][0]
    assert top["candidate_id"] == "L1" and top["source"] == "research"
    assert top["sum_name"] == "PC 34:1" and top["scores"]["matched_peaks_count"] == 2
    assert "no_matched_peaks" in out["candidates"][1]["soft"]
    assert "trend_unknown" in top["info"]


def test_msdial_candidate_merges_with_the_same_research_record(store):
    rec = store.candidates(804.5762, mz_tol=0.01, ion_mode="Negative")
    pc = next(r for r in rec if r["name"] == "PC 16:0_18:1")
    match = {"name": "PC 16:0_18:1", "inchikey": "KEY-PC341", "library_id": pc["record_index"],
             "annotator_id": "lib_1", "total_score": 3.1, "has_msms": True}
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0]]), msdial_matches=[match], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    merged = [c for c in out["candidates"] if c["name"] == "PC 16:0_18:1"]
    assert len(merged) == 1 and merged[0]["source"] == "msdial+research"
    assert merged[0]["msdial_total_score"] == 3.1


def test_current_representative_is_excluded_when_requested(store):
    rec = store.candidates(804.5762, mz_tol=0.01, ion_mode="Negative")
    pc = next(r for r in rec if r["name"] == "PC 16:0_18:1")
    representative = {"name": "PC 16:0_18:1", "library_id": pc["record_index"], "annotator_id": "lib_1",
                      "inchikey": "KEY-PC341"}
    out = candidates.build_library_candidates(
        _ev([[255.23, 900.0]]), msdial_matches=[], representative=representative, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=True, top_n=5)
    assert [c["name"] for c in out["candidates"]] == ["PS 16:0_18:1"]


def test_unresolved_msdial_candidate_is_kept_with_its_own_identity(store):
    match = {"name": "PG 34:1", "inchikey": "KEY-PG", "library_id": 999, "annotator_id": "lib_1",
             "total_score": 2.0, "has_msms": True}
    out = candidates.build_library_candidates(
        _ev([]), msdial_matches=[match], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    pg = next(c for c in out["candidates"] if c["name"] == "PG 34:1")
    assert pg["source"] == "msdial" and "reference_unresolved" in pg["info"]
    assert "msms_absent" in pg["soft"] and pg["inchikey"] == "KEY-PG"


def test_no_research_without_msms(store):
    out = candidates.build_library_candidates(
        _ev([]), msdial_matches=[], representative=None, store=store,
        scoring=SCORING, trends={"classes": {}}, th=TH, exclude_current=False, top_n=5)
    assert out["candidates"] == []


def test_trend_outlier_is_soft_and_ranks_below_clean():
    clean = {"soft": [], "scores": {"total_score": 1.0}}
    outlier = {"soft": ["trend_outlier"], "scores": {"total_score": 3.0}}
    assert candidates.rank([outlier, clean]) == [clean, outlier]


def test_dmz_out_is_hard():
    reasons = candidates.constraint_reasons(
        {"adduct": "[M-H]-", "ontology": "PS", "precursor_mz": 804.555}, rep_mz=804.5762,
        ion_mode="Negative", measured=[[1.0, 1.0]], scores={"matched_peaks_count": 1}, trend_entry=None,
        th=TH)
    assert "dmz_out" in reasons["hard"]
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_candidates.py tests/test_curation_trend.py tests/test_arf2_match_results.py -q`
Expected: FAIL（`candidates` / `sum_composition` / `usable_candidates` が無い）

- [ ] **Step 3: `match_results` と `trend` を書く**

`lipidmix/arf2/match_results.py` で、`representative` の中のフィルタを関数に出す:

```python
def usable_candidates(container) -> list[dict]:
    """非 decoy・Source != Unknown の候補（GUI の代表を選ぶ母集団）。並びはファイル順。"""
    return [c for c in _candidates(container)
            if not c.get("is_decoy") and (c.get("source") or 0) != SOURCE_UNKNOWN]


def representative(container) -> dict | None:
    """GUI が表示する代表（非 decoy・Source != Unknown の ResultOrder 最大）。無ければ None。"""
    usable = usable_candidates(container)
    if not usable:
        return None
    return max(usable, key=_order_key)


def load_spot_candidates(file_path) -> dict[int, list[dict]]:
    """MasterAlignmentID → 使える候補の一覧（代表を含む）。"""
    return {int(raw[0]): usable_candidates(raw[56] if len(raw) > 56 else None)
            for raw in load_raw_spots(file_path)}
```

`lipidmix/curation/trend.py` の末尾に足す:

```python
def sum_composition(name) -> str | None:
    """和組成の名前（pygoslin の SPECIES レベル。`PC 16:0_18:1` → `PC 34:1`）。読めなければ None。"""
    if not name or not str(name).strip():
        return None
    clean, _ = _clean_msdial_name(name)
    parser = _lipid_parser()
    if not clean or parser is None:
        return None
    try:
        from pygoslin.domain.LipidLevel import LipidLevel
        with _PARSER_LOCK:
            lipid = parser.parse(clean)
        return lipid.get_lipid_string(LipidLevel.SPECIES)
    except Exception:  # noqa: BLE001 - 脂質名として読めない名前は和組成を持たない
        return None


def predict_rt(classes: dict, ontology: str, carbon: int, db: int) -> dict | None:
    """`fit_trends` の `classes` からクラスの加法モデルで RT を予測する。当てはめの無いクラスは None。
    DB が 1 種類しかなかったクラスは係数が 2 つ（DB の項なし）。"""
    entry = classes.get(ontology or "")
    if not entry:
        return None
    coef = entry["coef"]
    predicted = coef[0] + coef[1] * carbon + (coef[2] * db if len(coef) > 2 else 0.0)
    return {"predicted_rt": predicted, "scale": entry["scale"]}
```

`predict_rt` のテストは浮動小数の和なので、`test_predict_rt_uses_the_additive_model` の比較は `pytest.approx` に直す（`assert predict_rt(...)["predicted_rt"] == pytest.approx(17.7)` と `["scale"] == 0.05` の 2 行ずつ）。

- [ ] **Step 4: `candidates.py` を書く**

```python
# lipidmix/curation/candidates.py
"""候補付けの ① MS-DIAL の下位候補と ② 閾値を緩めた再検索。spec §5。

両方を同じ物差し（`match_spectrum` + `total_score`、採点の前処理は store の search_params）で
採点し、同じレコードは 1 つにまとめる。制約（P）: ハードは削る、ソフトは順位を下げる、情報は
何もしない。並び順は（ソフト理由の数 昇順, total_score 降順）。
deps: analysis.spectral_match、library.store、msdial.adducts、msdial.peak_verification、
plots.mirror、curation.evidence（payload の間引き）、curation.trend。tools_* は import しない。
"""
from __future__ import annotations

from lipidmix.analysis.spectral_match import match_spectrum, total_score
from lipidmix.curation import trend
from lipidmix.curation.evidence import MIRROR_MAX_PEAKS, REFERENCE_MZ_WINDOW, _keep_most_intense_peaks
from lipidmix.library.store import library_id_from_annotator
from lipidmix.msdial.adducts import parse_adduct
from lipidmix.msdial.peak_verification import adduct_consistency
from lipidmix.plots.mirror import build_mirror_payload

SCORE_KEYS = ("total_score", "weighted_dot_product", "simple_dot_product", "reverse_dot_product",
              "matched_peaks_percentage", "matched_peaks_count")


def record_key(record: dict) -> tuple:
    return (record.get("library_id"), int(record["record_index"]))


def _polarity(ion_mode) -> str | None:
    mode = str(ion_mode or "").strip().lower()
    return "+" if mode.startswith("pos") else "-" if mode.startswith("neg") else None


def constraint_reasons(identity: dict, *, rep_mz, ion_mode, measured, scores, trend_entry, th) -> dict:
    hard, soft, info = [], [], []
    adduct = parse_adduct(identity.get("adduct"))
    polarity = _polarity(ion_mode)
    if adduct is not None and polarity is not None and adduct.polarity != polarity:
        hard.append("polarity_mismatch")
    ref_mz = identity.get("precursor_mz")
    if ref_mz is not None and rep_mz is not None and abs(float(rep_mz) - float(ref_mz)) * 1000 >= th["dmz_fail_mda"]:
        hard.append("dmz_out")
    if adduct_consistency(identity.get("adduct"), ion_mode, identity.get("ontology")).get("class_typical") is False:
        soft.append("adduct_atypical")
    if not measured:
        soft.append("msms_absent")
    elif (scores or {}).get("matched_peaks_count", 0) <= 0:
        soft.append("no_matched_peaks")
    if trend_entry is None:
        info.append("trend_unknown")
    elif abs(trend_entry["z"]) > th["trend_outlier_z"]:
        soft.append("trend_outlier")
    return {"hard": hard, "soft": soft, "info": info}


def rank(items: list[dict]) -> list[dict]:
    return sorted(items, key=lambda c: (len(c["soft"]), -(c["scores"].get("total_score") or 0.0)))


def _trend_entry(identity: dict, rt, trends: dict) -> dict | None:
    comp = trend.composition(identity.get("name"))
    if comp is None or rt is None:
        return None
    predicted = trend.predict_rt(trends.get("classes") or {}, identity.get("ontology") or "", *comp)
    if predicted is None:
        return None
    residual = float(rt) - predicted["predicted_rt"]
    return {"predicted_rt": round(predicted["predicted_rt"], 3), "residual": round(residual, 3),
            "z": round(residual / predicted["scale"], 2)}


def _score(measured, reference_spectrum, *, rep_mz, rt, ref_mz, ref_rt, scoring) -> dict:
    if not measured or not reference_spectrum:
        return {**{k: -1.0 for k in SCORE_KEYS}, "matched_peaks_count": 0, "alignment": []}
    result = match_spectrum(measured, reference_spectrum, ms2_tol=scoring["ms2_tol"],
                            mass_begin=scoring["mass_begin"], mass_end=scoring["mass_end"],
                            relative_amp_cutoff=scoring["relative_amp_cutoff"],
                            absolute_amp_cutoff=scoring["absolute_amp_cutoff"])
    result.update(total_score(result, precursor_mz=rep_mz, reference_precursor_mz=ref_mz,
                              ms1_tol=scoring["mz_tol"], rt=rt, reference_rt=ref_rt,
                              rt_tol=scoring.get("rt_tol"), use_rt=scoring["use_rt"]))
    return result


def _from_record(record: dict, source: str) -> dict:
    return {"source": source, "name": record["name"], "ontology": record.get("ontology") or record.get("compound_class"),
            "adduct": record.get("adduct"), "formula": record.get("formula"), "inchikey": record.get("inchikey"),
            "precursor_mz": record.get("precursor_mz"), "ref_rt": record.get("rt"),
            "library_id": record.get("library_id"), "record_index": record.get("record_index"),
            "_spectrum": record.get("spectrum") or [], "msdial_total_score": None}


def _from_match(match: dict) -> dict:
    return {"source": "msdial", "name": match.get("name"), "ontology": None, "adduct": None,
            "formula": None, "inchikey": match.get("inchikey"), "precursor_mz": None, "ref_rt": None,
            "library_id": library_id_from_annotator(match.get("annotator_id")),
            "record_index": match.get("library_id"), "_spectrum": [],
            "msdial_total_score": match.get("total_score"), "_unresolved": True}


def _resolve_match(store, match: dict, rep_mz) -> dict | None:
    if match.get("library_id") is None or rep_mz is None:
        return None
    record = store.record_by_scan_id(int(match["library_id"]),
                                     library_id=library_id_from_annotator(match.get("annotator_id")),
                                     precursor_mz=float(rep_mz), mz_tol=REFERENCE_MZ_WINDOW)
    if record is None:
        return None
    if match.get("inchikey") and record.get("inchikey") and match["inchikey"] != record["inchikey"]:
        return None
    return record


def build_library_candidates(ev: dict, *, msdial_matches, representative, store, scoring, trends, th,
                             exclude_current, top_n) -> dict:
    measured = ev.get("_measured") or []
    rep_mz, rt = ev.get("rep_mz"), ev.get("rep_rt")
    pool: dict[tuple, dict] = {}
    unresolved: list[dict] = []
    current = _resolve_match(store, representative, rep_mz) if representative else None
    current_key = record_key(current) if current else None

    for match in msdial_matches:
        record = _resolve_match(store, match, rep_mz)
        if record is None:
            unresolved.append(_from_match(match))
            continue
        item = pool.setdefault(record_key(record), _from_record(record, "msdial"))
        item["msdial_total_score"] = match.get("total_score")
    if measured and rep_mz is not None:
        for record in store.candidates(float(rep_mz), mz_tol=scoring["mz_tol"], ion_mode=ev.get("ion_mode")):
            key = record_key(record)
            if key in pool:
                if pool[key]["source"] == "msdial":
                    pool[key]["source"] = "msdial+research"
            else:
                pool[key] = _from_record(record, "research")
    if exclude_current and current_key is not None:
        pool.pop(current_key, None)
    if exclude_current and representative is not None:
        unresolved = [u for u in unresolved if u["name"] != representative.get("name")]

    kept, n_hard, matched_measured = [], 0, set()
    for item in [*pool.values(), *unresolved]:
        result = _score(measured, item["_spectrum"], rep_mz=rep_mz, rt=rt, ref_mz=item["precursor_mz"],
                        ref_rt=item["ref_rt"], scoring=scoring)
        trend_entry = _trend_entry(item, rt, trends)
        reasons = constraint_reasons(item, rep_mz=rep_mz, ion_mode=ev.get("ion_mode"), measured=measured,
                                     scores=result, trend_entry=trend_entry, th=th)
        if item.pop("_unresolved", False):
            reasons["info"].append("reference_unresolved")
        if reasons["hard"]:
            n_hard += 1
            continue
        mirror = None
        if measured and item["_spectrum"]:
            payload = build_mirror_payload(measured, item["_spectrum"], result["alignment"],
                                           title=item["name"] or "", ms2_tol=scoring["ms2_tol"])
            matched_measured.update(payload.get("matched_measured_mz") or ())
            mirror = {"reference": _keep_most_intense_peaks(
                          payload["reference"], lambda p, m=set(payload.get("matched_mz") or ()): p[0] in m,
                          MIRROR_MAX_PEAKS),
                      "matched_mz": payload.get("matched_mz") or [],
                      "matched_measured_mz": payload.get("matched_measured_mz") or []}
        item.pop("_spectrum")
        ppm = dmz = None
        if item["precursor_mz"] and rep_mz is not None:
            dmz = round((float(rep_mz) - item["precursor_mz"]) * 1000, 2)
            ppm = round((float(rep_mz) - item["precursor_mz"]) / item["precursor_mz"] * 1e6, 2)
        kept.append({**item, "kind": "library", "sum_name": trend.sum_composition(item["name"]),
                     "dmz_mda": dmz, "ppm": ppm,
                     "scores": {k: (round(result[k], 4) if isinstance(result[k], float) else result[k])
                                for k in SCORE_KEYS},
                     "trend": trend_entry, **reasons, "mirror": mirror})
    ranked = rank(kept)[:top_n]
    for index, candidate in enumerate(ranked, 1):
        candidate["candidate_id"] = f"L{index}"
    cut_measured = _keep_most_intense_peaks(measured, lambda p: p[0] in matched_measured, MIRROR_MAX_PEAKS)
    return {"candidates": ranked, "n_hard_removed": n_hard, "measured": cut_measured}
```

`build_mirror_payload` の戻り値のキー（`measured` / `reference` / `matched_mz` / `matched_measured_mz`）は `lipidmix/plots/mirror.py` を開いて確かめ、名前が違えば上のコードをそちらに合わせる（evidence.py の `_cut_mirror_for_payload` が同じキーを読んでいるので一致しているはず）。`total_score` の引数 `rt_tol` は `scoring.get("rt_tol")` で、`use_rt=False` のときは使われない。

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_candidates.py tests/test_curation_trend.py tests/test_arf2_match_results.py tests/test_curation_evidence.py -q`
Expected: PASS。`test_research_scores_all_in_window_and_filters_polarity` で `n_hard_removed` が 0 なら、`.msp` の PE レコードが IONMODE 欄なしで読まれて `ion_mode IS NULL` になっているか（`store.candidates` は NULL を残す）を `library.msp.iter_records` で確かめる。

- [ ] **Step 6: コミット**

```bash
git add lipidmix/arf2/match_results.py lipidmix/curation/trend.py lipidmix/curation/candidates.py tests/test_curation_candidates.py tests/test_curation_trend.py tests/test_arf2_match_results.py
git commit -m "feat(curation): 下位候補と再検索の候補を同じ物差しで採点し制約で並べる"
```

---

### Task 5: ④ イオン関係による説明

**Files:**
- Create: `lipidmix/curation/relations.py`
- Test: `tests/test_curation_relations.py`

**Interfaces:**
- Consumes: `msdial.adducts.parse_adduct` / `mz_from_neutral` / `neutral_from_mz`、`msdial.peak_verification.parse_formula`。
- Produces:
  - `INSOURCE_LOSSES: list[tuple[str, float]]`、`ISOTOPE_DELTAS: dict[str, tuple[int, float]]`（`{"isotope_M+1": (1, 1.003355), "isotope_M+2": (2, 2.006710)}`）
  - `expected_isotope_ratio(n_carbon: int, k: int) -> float`
  - `profile_correlation(x: list[float|None], y: list[float|None], min_n: int = 5) -> float | None`
  - `median_ratio(x, y) -> float | None`
  - `find_relations(target: dict, partners: list[dict], *, adducts: list[str], rt_window: float, mz_tol: float, min_r: float, links: dict[int, list[dict]]) -> list[dict]`
  - スポット dict（target / partners）: `{"spot_id", "name", "mz", "rt", "adduct", "formula", "heights": [float|None, ...]}`（`heights` は `.arf` の試料順。無ければ `[]`）
  - 関係 dict: `candidate_id`（`"R1"`…）、`kind="ion"`、`of`、`of_name`、`relation`、`dmz_mda`、`drt`、`r`、`msdial_links`、`isotope_ratio`、`expected_ratio`、`isotope_consistent`、`informational`、`strong`、`soft`（`["no_support"]` か `[]`）

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_curation_relations.py
import pytest

from lipidmix.curation import relations
from lipidmix.msdial.adducts import mz_from_neutral, parse_adduct
from lipidmix.msdial.peak_verification import monoisotopic_mass, parse_formula

PC342 = monoisotopic_mass(parse_formula("C42H80NO8P"))
Y_MZ = mz_from_neutral(PC342, parse_adduct("[M+HCOO]-"))
Y = {"spot_id": 0, "name": "PC 34:2", "mz": Y_MZ, "rt": 10.00, "adduct": "[M+HCOO]-",
     "formula": "C42H80NO8P", "heights": [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0]}
ADDUCTS = ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Cl]-"]
KW = {"adducts": ADDUCTS, "rt_window": 0.1, "mz_tol": 0.010, "min_r": 0.8}


def _x(mz, *, rt=10.01, heights=None):
    return {"spot_id": 9, "name": None, "mz": mz, "rt": rt, "adduct": None, "formula": None,
            "heights": heights if heights is not None else [h * 0.1 for h in Y["heights"]]}


def test_expected_isotope_ratio_is_binomial():
    assert relations.expected_isotope_ratio(43, 1) == pytest.approx(43 * 0.0107 / 0.9893, rel=1e-6)
    assert relations.expected_isotope_ratio(43, 2) == pytest.approx(903 * (0.0107 / 0.9893) ** 2, rel=1e-6)


def test_m_plus_2_with_a_consistent_ratio_is_a_strong_explanation():
    found = relations.find_relations(_x(Y_MZ + 2.006710), [Y], links={}, **KW)
    iso = next(r for r in found if r["relation"] == "isotope_M+2")
    assert iso["isotope_consistent"] is True and iso["strong"] is True and iso["of"] == 0
    assert iso["r"] == pytest.approx(1.0, abs=1e-2) and iso["candidate_id"] == "R1"   # log1p なので厳密に 1 ではない


def test_m_plus_2_that_is_too_intense_is_informational_only():
    found = relations.find_relations(_x(Y_MZ + 2.006710, heights=list(Y["heights"])), [Y], links={}, **KW)
    iso = next(r for r in found if r["relation"] == "isotope_M+2")
    assert iso["isotope_consistent"] is False and iso["informational"] is True and iso["strong"] is False


def test_adduct_pair_is_explained():
    x_mz = mz_from_neutral(PC342, parse_adduct("[M+CH3COO]-"))
    found = relations.find_relations(_x(x_mz), [Y], links={}, **KW)
    assert [r["relation"] for r in found] == ["adduct:[M+CH3COO]-/[M+HCOO]-"]


def test_insource_demethylation_from_formate():
    found = relations.find_relations(_x(Y_MZ - 60.021129), [Y], links={}, **KW)
    assert [r["relation"] for r in found] == ["insource:-HCOOCH3"]


def test_found_in_upper_msms_link_alone_makes_a_strong_candidate():
    found = relations.find_relations(_x(500.0), [Y], links={9: [{"spot_id": 0, "kind": "found_in_upper_msms"}]}, **KW)
    assert [(r["relation"], r["strong"]) for r in found] == [("found_in_upper_msms", True)]


def test_not_coeluting_is_not_related():
    assert relations.find_relations(_x(Y_MZ + 2.006710, rt=10.5), [Y], links={}, **KW) == []


def test_unsupported_relation_is_soft():
    found = relations.find_relations(_x(Y_MZ - 60.021129, heights=[5, 1, 4, 2, 6, 3]), [Y], links={}, **KW)
    assert found[0]["soft"] == ["no_support"] and found[0]["strong"] is False


def test_correlation_needs_five_samples_with_both_values():
    assert relations.profile_correlation([1, 2, 3, 4], [1, 2, 3, 4]) is None
    assert relations.profile_correlation([1, 2, None, 4, 5, 6], [2, 4, 6, 8, 10, 12]) == pytest.approx(1.0, abs=2e-2)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_relations.py -q`
Expected: FAIL（モジュールが無い）

- [ ] **Step 3: 実装を書く**

```python
# lipidmix/curation/relations.py
"""④ 注釈付きの別スポット Y のイオンとしての説明（純関数）。spec §6。

X と Y が同時に溶出し（|ΔRT| ≤ rt_window）、Δm/z が既知の関係（アダクトの組・同位体・インソース断片）
で説明できるか、MS-DIAL の FoundInUpperMsMs リンクがあるときに候補にする。裏付けは MS-DIAL の
リンク（correl_similar / chrom_similar / found_in_upper_msms）か、試料間の log 強度の相関 r ≥ min_r。
M+1/M+2 は強度比が炭素数から期待される比の 1.5 倍以下のときだけ「同位体で説明できる」とし、
それを超えれば X は実在の別物質とみなして情報として添えるだけにする（Type II の重なり）。
deps: msdial.adducts、msdial.peak_verification（組成式）。
"""
from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from statistics import median

from lipidmix.msdial.adducts import mz_from_neutral, neutral_from_mz, parse_adduct
from lipidmix.msdial.peak_verification import parse_formula

C13_ABUNDANCE = 0.0107
ISOTOPE_RATIO_FACTOR = 1.5
SUPPORT_LINKS = frozenset({"correl_similar", "chrom_similar", "found_in_upper_msms"})
ISOTOPE_DELTAS = {"isotope_M+1": (1, 1.003355), "isotope_M+2": (2, 2.006710)}
#: 出典は Task 1 Step 3 で確かめたうえで、各行の末尾コメントに書く。
INSOURCE_LOSSES = [
    ("insource:-H2O", 18.010565),
    ("insource:-2H2O", 36.021129),
    ("insource:-NH3", 17.026549),
    ("insource:-HCOOCH3", 60.021129),
    ("insource:-CH3COOCH3", 74.036779),
    ("insource:-C3H5NO2", 87.032028),
    ("insource:-C2H8NO4P", 141.019094),
    ("insource:-C3H8NO6P", 185.008923),
    ("insource:-C6H10O5", 162.052824),
]


def expected_isotope_ratio(n_carbon: int, k: int) -> float:
    q = C13_ABUNDANCE / (1 - C13_ABUNDANCE)
    return math.comb(n_carbon, k) * q ** k


def _carbon(formula) -> int | None:
    try:
        return parse_formula(formula).get("C") if formula else None
    except ValueError:
        return None


def _pairs(x, y):
    return [(a, b) for a, b in zip(x, y) if a is not None and b is not None and a > 0 and b > 0]


def profile_correlation(x, y, min_n: int = 5) -> float | None:
    pairs = _pairs(x, y)
    if len(pairs) < min_n:
        return None
    lx = [math.log1p(a) for a, _ in pairs]
    ly = [math.log1p(b) for _, b in pairs]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    sx = math.sqrt(sum((v - mx) ** 2 for v in lx))
    sy = math.sqrt(sum((v - my) ** 2 for v in ly))
    if sx == 0 or sy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / (sx * sy)


def median_ratio(x, y) -> float | None:
    pairs = _pairs(x, y)
    return median(a / b for a, b in pairs) if pairs else None


def _link_kinds(links: dict, a: int, b: int) -> list[str]:
    kinds = {l["kind"] for l in links.get(a, []) if l["spot_id"] == b}
    kinds |= {l["kind"] for l in links.get(b, []) if l["spot_id"] == a}
    return sorted(kinds)


def _explanations(x: dict, y: dict, adducts: list[str], mz_tol: float) -> list[tuple[str, float]]:
    """(relation, Δm/z [Da]) の一覧。"""
    found = []
    y_adduct = parse_adduct(y.get("adduct"))
    if y_adduct is not None:
        neutral = neutral_from_mz(y["mz"], y_adduct)
        for name in adducts:
            other = parse_adduct(name)
            if other is None or other.name == y_adduct.name:
                continue
            delta = x["mz"] - mz_from_neutral(neutral, other)
            if abs(delta) <= mz_tol:
                found.append((f"adduct:{other.name}/{y_adduct.name}", delta))
    charge = y_adduct.charge if y_adduct else 1
    for relation, (_, shift) in ISOTOPE_DELTAS.items():
        delta = x["mz"] - (y["mz"] + shift / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    for relation, loss in INSOURCE_LOSSES:
        delta = x["mz"] - (y["mz"] - loss / charge)
        if abs(delta) <= mz_tol:
            found.append((relation, delta))
    return found


def find_relations(target: dict, partners: list[dict], *, adducts, rt_window, mz_tol, min_r, links) -> list[dict]:
    ordered = sorted(partners, key=lambda p: p["rt"])
    rts = [p["rt"] for p in ordered]
    lo, hi = bisect_left(rts, target["rt"] - rt_window), bisect_right(rts, target["rt"] + rt_window)
    out = []
    for y in ordered[lo:hi]:
        if y["spot_id"] == target["spot_id"]:
            continue
        kinds = _link_kinds(links, target["spot_id"], y["spot_id"])
        r = profile_correlation(target.get("heights") or [], y.get("heights") or [])
        supported = bool(SUPPORT_LINKS & set(kinds)) or (r is not None and r >= min_r)
        explained = _explanations(target, y, adducts, mz_tol)
        if not explained and "found_in_upper_msms" in kinds:
            explained = [("found_in_upper_msms", None)]
        for relation, delta in explained:
            entry = {"kind": "ion", "of": y["spot_id"], "of_name": y.get("name"), "relation": relation,
                     "dmz_mda": None if delta is None else round(delta * 1000, 2),
                     "drt": round(target["rt"] - y["rt"], 3), "r": None if r is None else round(r, 3),
                     "msdial_links": kinds, "isotope_ratio": None, "expected_ratio": None,
                     "isotope_consistent": None, "informational": False}
            if relation in ISOTOPE_DELTAS:
                k = ISOTOPE_DELTAS[relation][0]
                carbon = _carbon(y.get("formula"))
                ratio = median_ratio(target.get("heights") or [], y.get("heights") or [])
                expected = expected_isotope_ratio(carbon, k) if carbon else None
                entry.update(isotope_ratio=None if ratio is None else round(ratio, 4),
                             expected_ratio=None if expected is None else round(expected, 4))
                if ratio is not None and expected is not None:
                    entry["isotope_consistent"] = ratio <= expected * ISOTOPE_RATIO_FACTOR
                    entry["informational"] = not entry["isotope_consistent"]
            entry["strong"] = bool(supported and not entry["informational"] and (
                relation == "found_in_upper_msms" or "found_in_upper_msms" in kinds
                or entry["isotope_consistent"] is True))
            entry["soft"] = [] if supported else ["no_support"]
            out.append(entry)
    out.sort(key=lambda e: (not e["strong"], e["informational"], bool(e["soft"]),
                            abs(e["dmz_mda"]) if e["dmz_mda"] is not None else 0.0))
    for index, entry in enumerate(out, 1):
        entry["candidate_id"] = f"R{index}"
    return out
```

`test_found_in_upper_msms_link_alone_makes_a_strong_candidate` では X の m/z が Y と関係の無い 500.0 で、リンクは裏付け（`SUPPORT_LINKS`）にも数えられるので `strong=True`。

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_relations.py -q`
Expected: PASS。`test_adduct_pair_is_explained` で同位体やインソースも同時に拾われて一覧が 1 件でなくなったら、Δm/z の偶然の一致なので、テストの期待を「`adduct:[M+CH3COO]-/[M+HCOO]-` を含む」に緩めず、まず 3 つの質量差（アダクト差 14.0157、同位体、損失）が 10 mDa 以内で重なっていないかを手計算で確かめる。

- [ ] **Step 5: コミット**

```bash
git add lipidmix/curation/relations.py tests/test_curation_relations.py
git commit -m "feat(curation): 別スポットの同位体・アダクト・インソース断片として説明する"
```

---

### Task 6: 候補付けの組み立て・保存・送信内容の展開

**Files:**
- Modify: `lipidmix/curation/evidence.py`（`collect(..., keep_measured=False)`、公開の `arf_rows_by_spot`）
- Create: `lipidmix/curation/suggest.py`
- Modify: `tests/curation_fixtures.py`（`write_suggest_set`）
- Test: `tests/test_curation_suggest.py`

**Interfaces:**
- Consumes: Task 1〜5 のすべて、`review.load_review` / `review.is_valid_review_id`、`flags.*`、`arf2.reader.load_catalog`、`arf2.match_results.load_spot_annotations` / `load_spot_candidates`。
- Produces:
  - `evidence.collect(..., keep_measured: bool = False)`: True のとき各 ev に `_measured`（全ピーク）を入れる。既存の呼び出しは変わらない。
  - `evidence.arf_rows_by_spot(arf_path) -> dict[int, list[dict]]`（既存 `_arf_rows` の公開名。`_arf_rows = arf_rows_by_spot` を残す）
  - `suggest.SUGGESTION_ID_RE`、`new_suggestion_id()`、`is_valid_suggestion_id(x)`
  - `suggest.latest_review(arf2_path, alignment_sha256) -> dict | None`
  - `suggest.select_targets(catalog, annotations, base_review, effective, *, wrong, unannotated, include_decided) -> list[dict]`（各要素 `{"spot": catalog spot, "target_kind": "flagged"|"likely_wrong"|"unannotated", "target_reasons": list[str]}`）
  - `suggest.run_suggestion(arf2_path, *, base_review, store, th, options) -> dict`（`options` = `{"wrong", "unannotated", "include_decided", "top_n", "rt_window", "relation_mz_tol", "relation_min_r"}`、`None` の値は既定に置き換える）
  - `suggest.save_suggestion(s) -> {"json": Path, "html": Path}`、`load_suggestion(dir_or_arf2, suggestion_id) -> dict`
  - `suggest.summary_tsv(s, max_rows) -> str`、`suggest.TSV_COLUMNS`
  - `suggest.expand_entries(entries, suggestion) -> list[dict]`（`flags.jsonl` に書く行の本体。不正なら `ValueError`）
  - 保存形: `{"suggestion_id", "created_at", "arf2_path", "alignment", "base_review_id", "library": {"sha256", "path"}, "thresholds", "options", "analysis_params", "counts", "warnings", "spots": [...]}`。各スポット: `spot_id, name, ontology, adduct, ion_mode, mz, rt, rep_mz, rep_rt, target_kind, target_reasons, current (代表の {name, inchikey} か None), eic, measured, candidates (library), relations (ion。上位 `PARTNER_EIC_RELATIONS`=3 件は `partner_eic` = Y の代表試料のトレース `{file_id, left, top, right, points}` か None を持つ), strong, preset (candidate_id か None)`

- [ ] **Step 1: fixture を書く**

`tests/curation_fixtures.py` の末尾に足す。負イオンの 5 スポット（spec の 5 つの状況をそれぞれ 1 つずつ）:

```python
from lipidmix.msdial.adducts import mz_from_neutral, parse_adduct
from lipidmix.msdial.peak_verification import monoisotopic_mass, parse_formula

SUGGEST_MSP = textwrap.dedent("""\
    NAME: PC 16:0_18:2
    PRECURSORMZ: {pc342:.4f}
    PRECURSORTYPE: [M+HCOO]-
    IONMODE: Negative
    INCHIKEY: KEY-PC342
    FORMULA: C42H80NO8P
    Num Peaks: 2
    255.23 999
    279.23 800

    NAME: PC 16:0_18:1
    PRECURSORMZ: {pc341:.4f}
    PRECURSORTYPE: [M+HCOO]-
    IONMODE: Negative
    INCHIKEY: KEY-PC341
    FORMULA: C42H82NO8P
    Num Peaks: 2
    255.23 999
    281.25 800

    NAME: PE 18:0_18:2
    PRECURSORMZ: {pe362:.4f}
    PRECURSORTYPE: [M-H]-
    IONMODE: Negative
    INCHIKEY: KEY-PE362
    FORMULA: C41H78NO8P
    Num Peaks: 2
    283.26 999
    279.23 700
""")


def _mz(formula: str, adduct: str) -> float:
    return mz_from_neutral(monoisotopic_mass(parse_formula(formula)), parse_adduct(adduct))


def write_suggest_set(folder: Path, *, with_param: bool = True, with_arf: bool = True) -> dict:
    """候補付け用の一式（負イオン、6 試料）。

    spot 0: PC 34:2 [M+HCOO]-、注釈付き（相手 Y）
    spot 1: PC 34:1 と注釈されているが、実は spot 0 の M+2（強度 0.1 倍・同時溶出）。wrong フラグの対象。
            MatchResults は代表（PC 16:0_18:1、ライブラリ 1 番）と、store に無い 2 位（PG 34:1、LibraryID 99）
    spot 2: 未注釈、MS/MS が PE 18:0_18:2 [M-H]- と合う（② で候補が出る）
    spot 3: 未注釈、spot 0 の [M+HCOO]- → [M-CH3]-（-60.0211）、spot 0 への FoundInUpperMsMs リンク
    spot 4: 未注釈、MS/MS なし、関係なし（候補は出ない）
    """
    folder.mkdir(parents=True, exist_ok=True)
    pc342 = _mz("C42H80NO8P", "[M+HCOO]-")
    pc341 = _mz("C42H82NO8P", "[M+HCOO]-")
    pe362 = _mz("C41H78NO8P", "[M-H]-")
    unknown = match_result({0: None, 1: "", 2: 0.0, 3: 0.0, 14: -1, 26: 1, 33: False})
    specs = [
        {"spot_id": 0, "name": "PC 34:2", "mz": pc342, "rt": 10.00, "adduct": "[M+HCOO]-",
         "formula": "C42H80NO8P", "ontology": "PC", "scale": 1.0, "spectrum": [(255.23, 999.0), (279.23, 800.0)],
         "matches": [match_result({0: "PC 16:0_18:2", 1: "KEY-PC342", 14: 0, 27: "lib_1"})], "links": []},
        {"spot_id": 1, "name": "PC 34:1", "mz": pc342 + 2.006710, "rt": 10.01, "adduct": "[M+HCOO]-",
         "formula": "C42H82NO8P", "ontology": "PC", "scale": 0.1, "spectrum": [(255.23, 999.0), (279.23, 700.0)],
         "matches": [match_result({0: "PC 16:0_18:1", 1: "KEY-PC341", 14: 1, 27: "lib_1", 31: 2}),
                     match_result({0: "PG 34:1", 1: "KEY-PG341", 14: 99, 27: "lib_1", 31: 1, 2: 2.0})],
         "links": [(0, 5)]},
        {"spot_id": 2, "name": "Unknown", "mz": pe362, "rt": 11.50, "adduct": "[M-H]-", "formula": "",
         "ontology": "", "scale": 0.5, "spectrum": [(283.26, 999.0), (279.23, 650.0)], "matches": [unknown],
         "links": []},
        {"spot_id": 3, "name": "Unknown", "mz": pc342 - 60.021129, "rt": 10.00, "adduct": "[M-H]-",
         "formula": "", "ontology": "", "scale": 0.3, "spectrum": [(255.23, 999.0)], "matches": [unknown],
         "links": [(0, 4)]},
        {"spot_id": 4, "name": "Unknown", "mz": 432.1234, "rt": 3.00, "adduct": "[M-H]-", "formula": "",
         "ontology": "", "scale": 0.7, "spectrum": [], "matches": [unknown], "links": []},
    ]
    write_arf2(folder / "AlignmentResult_x.arf2", [
        arf2_spot_raw(spot_id=s["spot_id"], name=s["name"], mz=s["mz"], rt=s["rt"], ontology=s["ontology"],
                      adduct=s["adduct"], formula=s["formula"], ion_mode=1, representative_file_id=0,
                      matches=s["matches"], peak_links=s["links"])
        for s in specs])
    (folder / "AlignmentResult_x.dcl").write_bytes(build_dcl_bytes([
        {"precursor_mz": s["mz"], "rt": s["rt"], "spectrum": s["spectrum"]} for s in specs]))
    eic_spots, groups = [], []
    for s in specs:
        samples, rows = [], []
        for file_id in range(6):
            height = 1000.0 * (file_id + 1) * s["scale"]
            samples.append({"file_id": file_id, "top": s["rt"], "left": s["rt"] - 0.1, "right": s["rt"] + 0.1,
                            "points": gaussian_points(s["rt"], height=height)})
            rows.append(arf_row(file_id=file_id, mz=s["mz"], rt=s["rt"], height=height, gap_filled=False))
        eic_spots.append({"rt": s["rt"], "mz": s["mz"], "samples": samples})
        groups.append(rows)
    (folder / "AlignmentResult_x.EIC.aef").write_bytes(css1_bytes(eic_spots))
    if with_arf:
        write_arf(folder / "AlignmentResult_x_PeakProperties.arf", groups)
    if with_param:
        (folder / "Dataset_x_param_1.txt").write_text(
            "Ion mode: Negative\nSearched adduct ions: [M-H]-,[M+HCOO]-,[M+CH3COO]-\n"
            "MS1 tolerance for centroid: 0.01\nRetention time tolerance for alignment: 0.1\n", encoding="utf-8")
    (folder / "lib.msp").write_text(SUGGEST_MSP.format(pc342=pc342, pc341=pc341, pe362=pe362), encoding="utf-8")
    return {"arf2": folder / "AlignmentResult_x.arf2", "msp": folder / "lib.msp"}
```

`arf2_spot_raw` の `formula=""` のとき Key 13 が `["", 0.0]` になり、`load_catalog` の `Formula` が空になることを確かめる（既存の `Unknown` スポットと同じ扱い）。`.msp` の record_index が 0 始まりで、PC 16:0_18:2 = 0、PC 16:0_18:1 = 1 であることを `library.msp.iter_records` で確かめ、違えば `match_result` の Key 14 を合わせる。

- [ ] **Step 2: 失敗するテストを書く**

```python
# tests/test_curation_suggest.py
import json

import pytest

from lipidmix.arf2.reader import load_catalog
from lipidmix.curation import evidence, flags, judge, review, suggest
from lipidmix.library import store as library_store
from tests.curation_fixtures import write_suggest_set

TH = judge.resolve_thresholds(None)
OPTIONS = {"wrong": "flagged_or_likely", "unannotated": True, "include_decided": False, "top_n": 5,
           "rt_window": None, "relation_mz_tol": None, "relation_min_r": None}


def _setup(tmp_path, monkeypatch, **kw):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_suggest_set(tmp_path / "neg", **kw)
    s = library_store.open_store(paths["msp"])
    spots = evidence.select_spots(load_catalog(paths["arf2"]), ontology=None, name_contains=None)
    base = review.run_review(paths["arf2"], spots, store=s, ms2_tol=0.025, th=TH, file_ids=None,
                             max_traces=12, selection={})
    review.save_review(base)
    flags.FlagStore(flags.curation_dir(paths["arf2"])).append(
        [{"spot_id": 1, "flag": "wrong", "note": ""}], alignment=base["alignment"],
        review_id=base["review_id"], source="user")
    return paths, s, base


@pytest.fixture()
def built(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)
    yield paths, s, base, result
    s.close()


def test_targets_are_the_flagged_and_the_unannotated(built):
    _, _, _, result = built
    kinds = {sp["spot_id"]: sp["target_kind"] for sp in result["spots"]}
    assert kinds == {1: "flagged", 2: "unannotated", 3: "unannotated", 4: "unannotated"}
    assert result["counts"]["targets"] == {"flagged": 1, "likely_wrong": 0, "unannotated": 3}


def test_isotope_explanation_is_preset_for_the_wrong_spot(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    assert spot1["strong"] is True
    preset = next(r for r in spot1["relations"] if r["candidate_id"] == spot1["preset"])
    assert preset["relation"] == "isotope_M+2" and preset["of"] == 0
    # 代表（PC 16:0_18:1）は除かれ、store に無い 2 位（PG 34:1）は reference_unresolved で残る
    names = [c["name"] for c in spot1["candidates"]]
    assert "PC 16:0_18:1" not in names and "PG 34:1" in names


def test_relation_carries_the_partner_eic(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    trace = spot1["relations"][0]["partner_eic"]
    assert trace["file_id"] == 0 and trace["left"] < trace["top"] < trace["right"]
    assert 0 < len(trace["points"]) <= evidence.EIC_MAX_POINTS


def test_research_candidate_for_an_unannotated_spot(built):
    _, _, _, result = built
    spot2 = next(sp for sp in result["spots"] if sp["spot_id"] == 2)
    assert spot2["candidates"][0]["name"] == "PE 18:0_18:2"
    assert spot2["candidates"][0]["source"] == "research"
    assert spot2["candidates"][0]["sum_name"] == "PE 36:2"
    assert spot2["preset"] is None


def test_insource_link_is_strong(built):
    _, _, _, result = built
    spot3 = next(sp for sp in result["spots"] if sp["spot_id"] == 3)
    assert spot3["relations"][0]["relation"] == "insource:-HCOOCH3" and spot3["relations"][0]["strong"]


def test_spot_without_msms_and_relation_has_nothing(built):
    _, _, _, result = built
    spot4 = next(sp for sp in result["spots"] if sp["spot_id"] == 4)
    assert spot4["candidates"] == [] and spot4["relations"] == []


def test_flagged_only_mode(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH,
                                    options={**OPTIONS, "wrong": "flagged", "unannotated": False})
    assert [sp["spot_id"] for sp in result["spots"]] == [1]
    s.close()


def test_run_suggestion_rejects_a_base_review_of_another_alignment(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    stale = {**base, "alignment": {**base["alignment"], "alignment_sha256": "other"}}
    with pytest.raises(ValueError, match="アラインメント"):
        suggest.run_suggestion(paths["arf2"], base_review=stale, store=s, th=TH, options=OPTIONS)
    s.close()


def test_run_suggestion_without_param_file_and_arf_uses_defaults(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch, with_param=False, with_arf=False)
    result = suggest.run_suggestion(paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)
    assert result["analysis_params"]["source"] == "default"
    assert any("param" in w for w in result["warnings"])
    spot3 = next(sp for sp in result["spots"] if sp["spot_id"] == 3)
    assert spot3["relations"][0]["r"] is None and spot3["relations"][0]["strong"]   # リンクだけで裏付け
    s.close()


def test_latest_review_finds_the_base(built):
    paths, _, base, _ = built
    found = suggest.latest_review(paths["arf2"], base["alignment"]["alignment_sha256"])
    assert found["review_id"] == base["review_id"]


def test_save_load_and_summary(built):
    paths, _, _, result = built
    saved = suggest.save_suggestion(result)
    assert saved["json"].name == f"suggest-{result['suggestion_id']}.json"
    loaded = suggest.load_suggestion(paths["arf2"], result["suggestion_id"])
    assert loaded["suggestion_id"] == result["suggestion_id"]
    lines = suggest.summary_tsv(loaded, max_rows=10).splitlines()
    assert lines[0].split("\t") == suggest.TSV_COLUMNS and len(lines) == 5
    with pytest.raises(ValueError):
        suggest.load_suggestion(paths["arf2"], "../evil")


def test_expand_entries_builds_assign_and_redundant(built):
    _, _, _, result = built
    rows = suggest.expand_entries([
        {"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "sum", "note": "ok"},
        {"spot_id": 1, "flag": "redundant", "candidate": "R1"},
        {"spot_id": 4, "flag": "clear"}], result)
    assign, redundant, clear = rows
    assert assign["name"] == "PE 36:2" and assign["species_name"] == "PE 18:0_18:2"
    assert assign["inchikey"] == "KEY-PE362" and assign["candidate_source"] == "research"
    assert assign["suggestion_id"] == result["suggestion_id"]
    assert redundant["of"] == 0 and redundant["relation"] == "isotope_M+2"
    assert clear == {"spot_id": 4, "flag": "clear", "note": ""}


def test_expand_entries_species_level_keeps_the_library_name(built):
    _, _, _, result = built
    row = suggest.expand_entries([{"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "species"}],
                                 result)[0]
    assert row["name"] == "PE 18:0_18:2" and row["level"] == "species"


def test_expand_entries_uses_msdial_identity_when_reference_unresolved(built):
    _, _, _, result = built
    spot1 = next(sp for sp in result["spots"] if sp["spot_id"] == 1)
    pg = next(c for c in spot1["candidates"] if c["name"] == "PG 34:1")
    row = suggest.expand_entries([{"spot_id": 1, "flag": "assign", "candidate": pg["candidate_id"],
                                   "level": "sum"}], result)[0]
    assert row["name"] == "PG 34:1" and row["inchikey"] == "KEY-PG341"


@pytest.mark.parametrize("entry, message", [
    ({"spot_id": 2, "flag": "assign", "candidate": "L9", "level": "sum"}, "candidate"),
    ({"spot_id": 0, "flag": "assign", "candidate": "L1", "level": "sum"}, "対象外"),
    ({"spot_id": 2, "flag": "wrong"}, "flag"),
    ({"spot_id": 2, "flag": "redundant", "candidate": "L1"}, "candidate"),
    ({"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "chain"}, "level"),
])
def test_expand_entries_rejects_unknown_candidate_and_spot(built, entry, message):
    _, _, _, result = built
    with pytest.raises(ValueError, match=message):
        suggest.expand_entries([entry], result)


def test_decided_spots_are_skipped_by_default(tmp_path, monkeypatch):
    paths, s, base = _setup(tmp_path, monkeypatch)
    flags.FlagStore(flags.curation_dir(paths["arf2"])).append(
        [{"spot_id": 2, "flag": "assign", "name": "PE 36:2", "level": "sum", "note": ""}],
        alignment=base["alignment"], review_id="cs-x", source="user")
    ids = [sp["spot_id"] for sp in suggest.run_suggestion(
        paths["arf2"], base_review=base, store=s, th=TH, options=OPTIONS)["spots"]]
    assert 2 not in ids
    ids = [sp["spot_id"] for sp in suggest.run_suggestion(
        paths["arf2"], base_review=base, store=s, th=TH, options={**OPTIONS, "include_decided": True})["spots"]]
    assert 2 in ids
    s.close()
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest.py -q`
Expected: FAIL（`suggest` が無い）

- [ ] **Step 4: `evidence.collect` を広げる**

`lipidmix/curation/evidence.py`: `_arf_rows` を `arf_rows_by_spot` に改名し、直後に `_arf_rows = arf_rows_by_spot` を置く。`collect` のシグネチャに `keep_measured: bool = False` を足し、`results.append({...})` の dict に 1 行足す:

```python
            **({"_measured": measured} if keep_measured else {}),
```

- [ ] **Step 5: `suggest.py` を書く**

```python
# lipidmix/curation/suggest.py
"""候補付けレビューの組み立て・保存・要約・送信内容の展開。spec §3〜§8。

1 回の候補付け = 対象（wrong フラグ・likely_wrong・未注釈）のスポットごとの ①②（candidates）と
④（relations）。`<arf2 のフォルダ>/curation/` に `suggest-<id>.json`（正準）と `.html` を書く。
送信用テキストは候補 ID だけを運び、記録行の中身はここで保存済みの候補から展開する。
deps: curation.{evidence,candidates,relations,flags,review,trend,viewer}、arf2.{reader,match_results,
ion_features}、msdial.analysis_params。tools_* は import しない。
"""
from __future__ import annotations

import json
import re
import secrets
from datetime import datetime, timezone
from pathlib import Path

from lipidmix.arf2.ion_features import load_ion_features
from lipidmix.arf2.match_results import load_spot_annotations, load_spot_candidates
from lipidmix.arf2.reader import load_catalog
from lipidmix.core.atomic_io import atomic_write_json
from lipidmix.curation import candidates, evidence, flags, relations, review, trend, viewer
from lipidmix.eic.reader import read_eic_spot_css1
from lipidmix.library.defaults import pick_tol
from lipidmix.msdial.analysis_params import resolve_analysis_params

SUGGESTION_ID_RE = re.compile(r"^cs-\d{8}-\d{6}-[0-9a-f]{4}$")
#: 相手 Y の EIC を持たせる関係の数（スポットごと、並び順の上から）。JSON を膨らませないため。
PARTNER_EIC_RELATIONS = 3
WRONG_MODES = ("flagged_or_likely", "flagged")
LEVELS = ("sum", "species")
DEFAULTS = {"wrong": "flagged_or_likely", "unannotated": True, "include_decided": False, "top_n": 5,
            "relation_mz_tol": 0.010, "relation_min_r": 0.8}
TSV_COLUMNS = ["spot_id", "name", "target", "top_candidate", "source", "total_score", "reasons",
               "relation", "strong"]


def new_suggestion_id() -> str:
    return f"cs-{datetime.now():%Y%m%d-%H%M%S}-{secrets.token_hex(2)}"


def is_valid_suggestion_id(value) -> bool:
    return isinstance(value, str) and SUGGESTION_ID_RE.fullmatch(value) is not None


def latest_review(arf2_path, alignment_sha256: str) -> dict | None:
    """そのアラインメント（sha256 一致）の最新のレビュー。ID の時刻順に新しいものから読む。"""
    directory = flags.curation_dir(arf2_path)
    for path in sorted(directory.glob("review-cr-*.json"), key=lambda p: p.name, reverse=True):
        review_id = path.stem[len("review-"):]
        if not review.is_valid_review_id(review_id):
            continue
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved.get("alignment", {}).get("alignment_sha256") == alignment_sha256:
            return saved
    return None


def _is_unknown_name(name) -> bool:
    return (name or "").strip().lower() in ("", "unknown")


def select_targets(catalog, annotations, base_review, effective, *, wrong, unannotated, include_decided):
    decisions = flags.split_decisions(effective)
    decided = set(decisions["assign"]) | set(decisions["redundant"])
    verdicts = {s["spot_id"]: s for s in base_review["spots"]}
    targets = []
    for spot in catalog:
        spot_id = spot["MasterAlignmentID"]
        if spot_id in decided and not include_decided:
            continue
        annotated = (annotations.get(spot_id) or {}).get("representative") is not None
        judged = verdicts.get(spot_id)
        if annotated and spot_id in decisions["wrong"]:
            targets.append({"spot": spot, "target_kind": "flagged",
                            "target_reasons": list((judged or {}).get("reasons") or [])})
        elif annotated and wrong == "flagged_or_likely" and judged and judged["verdict"] == "likely_wrong":
            targets.append({"spot": spot, "target_kind": "likely_wrong",
                            "target_reasons": list(judged.get("reasons") or [])})
        elif not annotated and unannotated:
            targets.append({"spot": spot, "target_kind": "unannotated", "target_reasons": []})
    return targets


def _heights(rows: list[dict]) -> list:
    return [row.get("height") for row in sorted(rows, key=lambda r: r.get("file_id") or 0)]


def _partner(spot: dict, annotation: dict, rows) -> dict:
    rep = annotation.get("representative") or {}
    return {"spot_id": spot["MasterAlignmentID"], "name": spot.get("Name") or rep.get("name"),
            "mz": spot.get("MassCenter"), "rt": spot.get("RT"), "adduct": spot.get("AdductType"),
            "formula": spot.get("Formula"), "heights": _heights(rows)}


def _partner_trace(eic_path, spot_id: int, annotations: dict) -> dict | None:
    """相手 Y の代表試料の EIC 1 本（ビューアの重ね描き用。spec §7.2）。読めなければ None。"""
    rep_file = (annotations.get(spot_id) or {}).get("representative_file_id")
    if eic_path is None or rep_file is None:
        return None
    try:
        eic = read_eic_spot_css1(str(eic_path), spot_id, [rep_file], max_traces=1, max_total_points=5000)
    except (ValueError, OSError):
        return None
    if not eic["samples"]:
        return None
    sample = eic["samples"][0]
    left, right = sample["peak_left"], sample["peak_right"]
    points = evidence._trim(sample["chromatogram"], left, right)
    return {"file_id": sample["file_id"], "left": round(left, 4), "top": round(sample["peak_top"], 4),
            "right": round(right, 4), "points": evidence._downsample_points(points, left, right)}


def _scoring(store, rt_tol_default=None) -> dict:
    sp = store.summary().get("search_params") or {}
    return {"mz_tol": pick_tol(None, sp, "ms1_tolerance", 0.01),
            "ms2_tol": pick_tol(None, sp, "ms2_tolerance", 0.025),
            "rt_tol": pick_tol(None, sp, "rt_tolerance", 0.2),
            "mass_begin": pick_tol(None, sp, "mass_range_begin", 0.0),
            "mass_end": pick_tol(None, sp, "mass_range_end", 2000.0),
            "relative_amp_cutoff": pick_tol(None, sp, "relative_amp_cutoff", 0.0),
            "absolute_amp_cutoff": pick_tol(None, sp, "absolute_amp_cutoff", 0.0),
            "use_rt": bool(sp.get("use_time_for_annotation_scoring", False))}


def run_suggestion(arf2_path, *, base_review, store, th, options) -> dict:
    arf2_path = Path(arf2_path).resolve()
    opts = {**DEFAULTS, **{k: v for k, v in (options or {}).items() if v is not None}}
    if opts["wrong"] not in WRONG_MODES:
        raise ValueError(f"wrong は {WRONG_MODES} のいずれかにしてください。")
    alignment = flags.alignment_key(arf2_path)
    if base_review["alignment"]["alignment_sha256"] != alignment["alignment_sha256"]:
        raise ValueError("元のレビューは別の版のアラインメント（.arf2）に対して作られています。"
                         "curation_review をやり直してください。")
    effective = flags.FlagStore(flags.curation_dir(arf2_path)).effective(alignment["alignment_sha256"])
    catalog = load_catalog(arf2_path)
    annotations = load_spot_annotations(arf2_path)
    all_candidates = load_spot_candidates(arf2_path)
    ion_features = load_ion_features(arf2_path)
    links = {spot_id: f["links"] for spot_id, f in ion_features.items()}
    ion_mode = next((s.get("IonMode") for s in catalog if s.get("IonMode")), None)
    params = resolve_analysis_params(arf2_path, ion_mode)
    rt_window = opts.get("rt_window") or params["rt_window"]
    warnings = []
    if params["source"] == "default":
        warnings.append("解析の param ファイル（*_param_*.txt）が見つからないので、既定の検索アダクトと "
                        f"RT 窓 {rt_window} 分を使いました。")

    targets = select_targets(catalog, annotations, base_review, effective, wrong=opts["wrong"],
                             unannotated=opts["unannotated"], include_decided=opts["include_decided"])
    if len(targets) > evidence.MAX_SPOTS:
        raise ValueError(f"対象が {len(targets)} 件あり、上限 {evidence.MAX_SPOTS} を超えています。")
    target_ids = {t["spot"]["MasterAlignmentID"] for t in targets}
    decisions = flags.split_decisions(effective)
    excluded_partners = target_ids | decisions["wrong"] | set(decisions["redundant"])

    files = evidence.sibling_files(arf2_path)
    rows_by_spot = evidence.arf_rows_by_spot(files["arf"])
    if files["arf"] is None:
        warnings.append("PeakProperties.arf が無いので試料間の相関を計算できません（MS-DIAL のリンクだけで裏付けを判定）。")
    partners, points = [], []
    for spot in catalog:
        spot_id = spot["MasterAlignmentID"]
        annotation = annotations.get(spot_id) or {}
        if annotation.get("representative") is None or spot_id in excluded_partners:
            continue
        partners.append(_partner(spot, annotation, rows_by_spot.get(spot_id, [])))
        comp = trend.composition(spot.get("Name"))
        if comp is not None and spot.get("RT") is not None:
            points.append({"spot_id": spot_id, "ontology": spot.get("Ontology") or "", "rt": spot["RT"],
                           "mz": spot.get("MassCenter"), "carbon": comp[0], "db": comp[1]})
    trends = trend.fit_trends(points, th)

    scoring = _scoring(store)
    evs, stats = evidence.collect(arf2_path, [t["spot"] for t in targets], store=store,
                                  ms2_tol=scoring["ms2_tol"], th=th, keep_measured=True)
    by_id = {t["spot"]["MasterAlignmentID"]: t for t in targets}
    spots_out, n_hard = [], 0
    for ev in evs:
        target = by_id[ev["spot_id"]]
        annotation = annotations.get(ev["spot_id"]) or {}
        rep = annotation.get("representative")
        msdial = [c for c in all_candidates.get(ev["spot_id"], []) if c is not rep and c != rep]
        lib = candidates.build_library_candidates(
            ev, msdial_matches=msdial, representative=rep, store=store, scoring=scoring, trends=trends,
            th=th, exclude_current=target["target_kind"] != "unannotated", top_n=opts["top_n"])
        n_hard += lib["n_hard_removed"]
        me = _partner(target["spot"], annotation, rows_by_spot.get(ev["spot_id"], []))
        ion = relations.find_relations(me, partners, adducts=params["searched_adducts"], rt_window=rt_window,
                                       mz_tol=opts["relation_mz_tol"], min_r=opts["relation_min_r"],
                                       links=links)
        for relation in ion[:PARTNER_EIC_RELATIONS]:
            relation["partner_eic"] = _partner_trace(files["eic"], relation["of"], annotations)
        strong = next((r for r in ion if r["strong"]), None)
        ev.pop("_measured", None)
        spots_out.append({
            "spot_id": ev["spot_id"], "name": ev["name"], "ontology": ev["ontology"], "adduct": ev["adduct"],
            "ion_mode": ev["ion_mode"], "mz": ev["mz"], "rt": ev["rt"], "rep_mz": ev["rep_mz"],
            "rep_rt": ev["rep_rt"], "target_kind": target["target_kind"],
            "target_reasons": target["target_reasons"],
            "current": None if rep is None else {"name": rep.get("name"), "inchikey": rep.get("inchikey")},
            "eic": ev["eic"], "measured": lib["measured"], "candidates": lib["candidates"],
            "relations": ion, "strong": strong is not None,
            "preset": strong["candidate_id"] if strong else None})

    counts = {"targets": {k: sum(1 for s in spots_out if s["target_kind"] == k)
                          for k in ("flagged", "likely_wrong", "unannotated")},
              "with_candidates": sum(1 for s in spots_out if s["candidates"] or
                                     any(not r["informational"] for r in s["relations"])),
              "with_strong_relation": sum(1 for s in spots_out if s["strong"]),
              "hard_removed": n_hard}
    if stats["missing_files"]:
        warnings.append(f"兄弟ファイルがありません: {', '.join(stats['missing_files'])}。")
    warnings.append("likely_wrong は元のレビュー（" + base_review["review_id"] + "）の対象範囲についてだけ分かります。")
    return {"suggestion_id": new_suggestion_id(), "created_at": datetime.now(timezone.utc).isoformat(),
            "arf2_path": str(arf2_path), "alignment": alignment, "base_review_id": base_review["review_id"],
            "library": {"sha256": store.summary().get("source_sha256")},
            "thresholds": th, "options": {**opts, "rt_window": rt_window}, "analysis_params": params,
            "scoring": scoring, "counts": counts, "warnings": warnings, "spots": spots_out}


def save_suggestion(s: dict) -> dict:
    directory = flags.curation_dir(s["arf2_path"])
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / f"suggest-{s['suggestion_id']}.json"
    html_path = directory / f"suggest-{s['suggestion_id']}.html"
    atomic_write_json(json_path, s)
    html_path.write_text(viewer.render_suggest_html(s), encoding="utf-8")
    return {"json": json_path, "html": html_path}


def load_suggestion(arf2_path_or_dir, suggestion_id: str) -> dict:
    if not is_valid_suggestion_id(suggestion_id):
        raise ValueError(f"suggestion_id={suggestion_id!r} の形が不正です（cs-YYYYMMDD-HHMMSS-xxxx）。")
    base = Path(arf2_path_or_dir)
    directory = base if base.is_dir() else flags.curation_dir(base)
    path = directory / f"suggest-{suggestion_id}.json"
    if not path.is_file():
        raise FileNotFoundError(f"suggestion_id={suggestion_id} の候補付けがありません: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def summary_tsv(s: dict, max_rows: int | None = None) -> str:
    order = {"flagged": 0, "likely_wrong": 1, "unannotated": 2}
    spots = sorted(s["spots"], key=lambda sp: (not sp["strong"], order[sp["target_kind"]], sp["spot_id"]))
    if max_rows is not None:
        spots = spots[:max(0, max_rows)]
    lines = ["\t".join(TSV_COLUMNS)]
    for sp in spots:
        top = sp["candidates"][0] if sp["candidates"] else None
        rel = sp["relations"][0] if sp["relations"] else None
        cells = [sp["spot_id"], sp["name"], sp["target_kind"],
                 top and (top["sum_name"] or top["name"]), top and top["source"],
                 top and top["scores"]["total_score"], top and ",".join(top["soft"] + top["info"]),
                 rel and f"{rel['relation']}(#{rel['of']})", sp["strong"]]
        lines.append("\t".join(review._cell(c) for c in cells))
    return "\n".join(lines)


def _note(entry) -> str:
    note = str(entry.get("note") or "")
    for control in ("\r", "\n", "\t"):
        note = note.replace(control, " ")
    return note[:500]


def expand_entries(entries, suggestion: dict) -> list[dict]:
    if not isinstance(entries, list) or not entries:
        raise ValueError("flags は 1 件以上のリストで渡してください。")
    spots = {sp["spot_id"]: sp for sp in suggestion["spots"]}
    out = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"flags[{index}] が dict ではありません。")
        spot_id, flag = entry.get("spot_id"), entry.get("flag")
        if isinstance(spot_id, bool) or not isinstance(spot_id, int):
            raise ValueError(f"flags[{index}].spot_id が整数ではありません: {spot_id!r}")
        if spot_id not in spots:
            raise ValueError(f"flags[{index}].spot_id={spot_id} はこの候補付けの対象外です。")
        if flag not in flags.SUGGEST_FLAG_VALUES:
            raise ValueError(f"flags[{index}].flag={flag!r} は {flags.SUGGEST_FLAG_VALUES} のいずれかにしてください。")
        base = {"spot_id": spot_id, "flag": flag, "note": _note(entry)}
        if flag == "clear":
            out.append(base)
            continue
        pool = spots[spot_id]["candidates"] if flag == "assign" else spots[spot_id]["relations"]
        chosen = next((c for c in pool if c["candidate_id"] == entry.get("candidate")), None)
        if chosen is None or (flag == "redundant" and chosen.get("informational")):
            raise ValueError(f"flags[{index}].candidate={entry.get('candidate')!r} はスポット {spot_id} の"
                             f"{'候補' if flag == 'assign' else 'イオン関係'}にありません。")
        if flag == "assign":
            level = entry.get("level") or "sum"
            if level not in LEVELS:
                raise ValueError(f"flags[{index}].level={level!r} は {LEVELS} のいずれかにしてください。")
            name = (chosen["sum_name"] or chosen["name"]) if level == "sum" else chosen["name"]
            out.append({**base, "name": name, "level": level, "species_name": chosen["name"],
                        "ontology": chosen.get("ontology"), "adduct": chosen.get("adduct"),
                        "formula": chosen.get("formula"), "inchikey": chosen.get("inchikey"),
                        "candidate_source": chosen["source"],
                        "library": {"sha256": suggestion["library"].get("sha256"),
                                    "library_id": chosen.get("library_id"),
                                    "record_index": chosen.get("record_index")},
                        "scores": chosen["scores"], "suggestion_id": suggestion["suggestion_id"]})
        else:
            out.append({**base, "of": chosen["of"], "relation": chosen["relation"],
                        "evidence": {k: chosen.get(k) for k in ("dmz_mda", "drt", "r", "msdial_links")},
                        "suggestion_id": suggestion["suggestion_id"]})
    return out
```

`viewer.render_suggest_html` は Task 7 で作る。Task 6 のテストを先に通すため、この時点で `lipidmix/curation/viewer.py` に仮の実装を置く:

```python
def render_suggest_html(suggestion: dict | None) -> str:
    return "<!doctype html><title>suggest</title>"   # Task 7 で本実装に置き換える
```

`msdial = [c for c in ... if c is not rep and c != rep]`: `load_spot_annotations` と `load_spot_candidates` は別々に復号するので、同一性（`is`）では外れない。`!=`（dict の等値）で代表を除く。

- [ ] **Step 6: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest.py tests/test_curation_evidence.py tests/test_curation_review.py -q`
Expected: PASS。`test_isotope_explanation_is_preset_for_the_wrong_spot` が落ちたら、fixture の spot 1 と spot 0 の高さの比（0.1）が C43 の期待比 0.1056 × 1.5 以下か、spot 0 の `Formula` が `load_catalog` で `C42H80NO8P` として読めているか（`_carbon` は Y の組成式から炭素数を取る）を確かめる。ここで使う炭素数はアダクト（HCOO）の炭素を含まない 42 なので、期待比は `expected_isotope_ratio(42, 2)` ≈ 0.1006 で、0.1 ≤ 0.151 を満たす。

- [ ] **Step 7: コミット**

```bash
git add lipidmix/curation/evidence.py lipidmix/curation/suggest.py lipidmix/curation/viewer.py tests/curation_fixtures.py tests/test_curation_suggest.py
git commit -m "feat(curation): 候補付けレビューを組み立てて保存し、送信内容を候補から展開する"
```

---

### Task 7: ビューア（共通描画の切り出しと候補付けビューア）

**Files:**
- Create: `lipidmix/curation/viewer_common.js`
- Modify: `lipidmix/curation/viewer.html`、`lipidmix/curation/viewer.py`
- Create: `lipidmix/curation/suggest_viewer.html`
- Test: `tests/test_curation_review.py`（既存のビューアのテストがそのまま通ること）、`tests/test_curation_suggest_viewer.py`

**Interfaces:**
- Consumes: Task 6 の保存形（`spots[].candidates` / `relations` / `measured` / `eic` / `preset`）。
- Produces: `viewer.render_html(review)`（挙動は現行と同じ）、`viewer.render_suggest_html(suggestion | None) -> str`。JS の純関数ブロック（`// --- suggest (pure) ---` 〜 `// --- end suggest ---`）: `initialChoice(spot)`、`choiceEntry(spot, choice)`、`suggestSubmission(suggestion, choices)`、`overlayFor(spot, relation)`（X のトレース＋Y の `partner_eic` を X の代表の最大値に合わせて縮めたもの）、`mirrorFor(spot, candidate)`。共通の `drawEic` は `partner: true` のトレースを参照色の点線で描く。

- [ ] **Step 1: 共通 JS を切り出す**

`viewer.html` の `<script>` の中から、`const css = ...` / `const el = ...`（2 つの定数）と、`function preparerFor` から `function drawMirror` の閉じ括弧まで（`const fixed`、`drawEic`、`// --- mirror labels (pure) ---` ブロック、`drawMirror` を含む）をそのまま `lipidmix/curation/viewer_common.js` に移し、元の場所に 1 行 `/*__VIEWER_COMMON__*/` を置く。`drawMirror` の中の「参照スペクトルなし」の分岐は `spot.match` を読むので、共通側では次の形に直す（`spot.match` が無いときも動く）:

```javascript
  if (!m) { ctx.fillStyle = css("--muted"); ctx.font = "12px system-ui";
    ctx.fillText(spot.mirror_empty_text || (spot.match && !spot.match.has_msms ? "MS/MS なし" : "参照スペクトルなし"), 8, h / 2); return; }
```

同じく共通側の `drawEic` のトレースを描くループで、相手 Y のトレース（`s.partner`、Step 4 の `overlayFor` が作る）を描き分ける。既存レビューのデータには `partner` が無いので、既存ビューアの見た目は変わらない:

```javascript
  for (const s of samples) {
    ctx.beginPath();
    ctx.strokeStyle = s.partner ? css("--ref") : s.representative ? css("--meas") : css("--muted");
    ctx.lineWidth = s.representative || s.partner ? 2 : 1;
    ctx.setLineDash(s.partner ? [2, 2] : s.detected ? [] : [4, 3]);
    s.points.forEach((p, i) => i ? ctx.lineTo(px(p[0]), py(p[1])) : ctx.moveTo(px(p[0]), py(p[1])));
    ctx.stroke();
  }
```

`viewer.py`:

```python
_TEMPLATE = Path(__file__).with_name("viewer.html")
_SUGGEST_TEMPLATE = Path(__file__).with_name("suggest_viewer.html")
_COMMON_JS = Path(__file__).with_name("viewer_common.js")
_PLACEHOLDER = "/*__CURATION_DATA__*/null"
_COMMON_PLACEHOLDER = "/*__VIEWER_COMMON__*/"


def _embed(data) -> str:
    if data is None:
        return "null"
    return json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False,
                      default=str).replace("<", "\\u003c")


def _render(template_path: Path, data) -> str:
    template = template_path.read_text(encoding="utf-8")
    template = template.replace(_COMMON_PLACEHOLDER, _COMMON_JS.read_text(encoding="utf-8"), 1)
    template = template.replace("__SUBMISSION_PREFIX__", SUBMISSION_PREFIX.strip())
    return template.replace(_PLACEHOLDER, _embed(data), 1)


def render_html(review: dict | None) -> str:
    return _render(_TEMPLATE, review)


def render_suggest_html(suggestion: dict | None) -> str:
    return _render(_SUGGEST_TEMPLATE, suggestion)
```

置換の順序は docstring のとおり（テンプレート側の置換を先に、データの埋め込みを最後に）。共通 JS に `__SUBMISSION_PREFIX__` や `/*__CURATION_DATA__*/` を書かないこと。

Run: `C:/Python314/python.exe -m pytest tests/test_curation_review.py tests/test_curation_tools.py -q`
Expected: PASS（切り出しで既存ビューアの挙動が変わっていない。`function drawMirror` 〜 `function drawTrend` を切り出す既存テストは、差し込み後の HTML で `drawMirror` の直後に `drawTrend` が来るので通る——来なければ差し込みの位置を `drawTrend` の直前に直す）。

- [ ] **Step 2: 失敗するテストを書く**

```python
# tests/test_curation_suggest_viewer.py
import json
import shutil
import subprocess

import pytest

from lipidmix.curation import viewer

SPOT = {"spot_id": 1, "target_kind": "flagged", "preset": "R1", "measured": [[255.23, 999.0]],
        "candidates": [{"candidate_id": "L1", "name": "PG 16:0_18:1", "sum_name": "PG 34:1",
                        "mirror": {"reference": [[255.23, 999.0]], "matched_mz": [255.23],
                                   "matched_measured_mz": [255.23]}}],
        "relations": [{"candidate_id": "R1", "informational": False},
                      {"candidate_id": "R2", "informational": True}]}


def _run(tmp_path, expression, prelude=""):
    node = shutil.which("node")
    if not node:
        pytest.skip("node が無いのでビューアの純関数を実行できない")
    html = viewer.render_suggest_html(None)
    begin, end = "// --- suggest (pure) ---", "// --- end suggest ---"
    source = html[html.index(begin):html.index(end)]
    script = tmp_path / "s.js"
    script.write_text(source + "\n" + prelude + f"\nconsole.log(JSON.stringify({expression}));\n", encoding="utf-8")
    result = subprocess.run([node, str(script)], capture_output=True, text=True, encoding="utf-8")
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_initial_choice_is_the_strong_relation_only(tmp_path):
    spot = json.dumps(SPOT)
    assert _run(tmp_path, f"initialChoice({spot})") == {"candidate": "R1", "flag": "redundant", "level": "sum"}
    no_preset = json.dumps({**SPOT, "preset": None})
    assert _run(tmp_path, f"initialChoice({no_preset})") is None


def test_choice_entries_carry_only_ids(tmp_path):
    spot = json.dumps(SPOT)
    entry = _run(tmp_path, f"choiceEntry({spot}, {{candidate: 'L1', flag: 'assign', level: 'species', note: 'x'}})")
    assert entry == {"spot_id": 1, "flag": "assign", "candidate": "L1", "level": "species", "note": "x"}
    restore = _run(tmp_path, f"choiceEntry({spot}, {{flag: 'clear', note: ''}})")
    assert restore == {"spot_id": 1, "flag": "clear", "note": ""}


def test_submission_uses_the_suggestion_id(tmp_path):
    out = _run(tmp_path, "suggestSubmission({suggestion_id: 'cs-1', arf2_path: 'C:/a.arf2', spots: []}, new Map([[1, {candidate: 'R1', flag: 'redundant'}]]))")
    assert out == {"review_id": "cs-1", "arf2_path": "C:/a.arf2",
                   "flags": [{"spot_id": 1, "flag": "redundant", "candidate": "R1", "note": ""}]}


def test_mirror_for_combines_spot_measured_with_candidate_reference(tmp_path):
    spot = json.dumps(SPOT)
    out = _run(tmp_path, f"mirrorFor({spot}, {spot}.candidates[0])")
    assert out == {"mirror": {"measured": [[255.23, 999.0]], "reference": [[255.23, 999.0]],
                              "matched_mz": [255.23]}}


def test_overlay_adds_the_partner_trace_scaled_to_the_representative(tmp_path):
    spot = json.dumps({"eic": {"samples": [
        {"file_id": 0, "representative": True, "detected": True, "points": [[10.0, 0.0], [10.1, 100.0]]},
        {"file_id": 1, "representative": False, "detected": True, "points": [[10.0, 0.0], [10.1, 50.0]]}]}})
    relation = json.dumps({"partner_eic": {"file_id": 0, "left": 9.9, "top": 10.1, "right": 10.3,
                                           "points": [[10.0, 0.0], [10.1, 1000.0]]}})
    out = _run(tmp_path, f"overlayFor({spot}, {relation}).eic.samples")
    assert len(out) == 3
    assert out[2]["partner"] is True and out[2]["representative"] is False
    assert out[2]["points"][1] == [10.1, 100.0]                     # Y の最大 1000 → X の代表の最大 100
    assert _run(tmp_path, f"overlayFor({spot}, {{partner_eic: null}}).eic.samples.length") == 2


def test_eic_draws_the_partner_in_its_own_style():
    html = viewer.render_suggest_html(None)
    source = html[html.index("function drawEic"):html.index("// --- mirror labels (pure) ---")]
    assert "s.partner" in source


def test_embedded_data_is_escaped():
    html = viewer.render_suggest_html({"suggestion_id": "cs-1", "spots": [{"name": "</script>"}]})
    assert "</script>\"" not in html and "\\u003c/script>" in html
```

- [ ] **Step 3: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest_viewer.py -q`
Expected: FAIL（仮の `render_suggest_html` に純関数ブロックが無い）

- [ ] **Step 4: `suggest_viewer.html` を書く**

`viewer.html` の `<head>`（`<style>` の `:root` 変数・ダークモード・`header` / `main` / `.card` / `.badge` / `canvas` / `#zoom` の規則）をそのまま写し、次の規則を足す:

```css
.cand { display:grid; grid-template-columns:auto 1fr; gap:2px 8px; align-items:start; padding:4px 0;
        border-top:1px solid var(--line); font-size:12px; }
.cand .src { display:inline-block; border-radius:4px; padding:0 5px; background:var(--line); }
.cand.info { opacity:.6; }
.card.strong { border:2px solid var(--match); }
.card.chosen { outline:2px solid var(--meas); }
```

本体:

```html
<body>
<header>
  <h1 id="title">Annotation suggestions</h1>
  <label>対象 <select id="f-target"><option value="">すべて</option><option value="flagged">wrong フラグ</option>
    <option value="likely_wrong">likely_wrong</option><option value="unannotated">未注釈</option></select></label>
  <label>説明 <select id="f-strong"><option value="">すべて</option><option value="strong">強い説明あり</option></select></label>
  <label>出所 <select id="f-source"><option value="">すべて</option><option value="msdial">msdial を含む</option>
    <option value="research">research を含む</option><option value="ion">イオン関係あり</option></select></label>
  <label>選択 <select id="f-chosen"><option value="">すべて</option><option value="chosen">選んだもの</option>
    <option value="open">未選択</option></select></label>
  <button id="copy" class="primary">送信用テキストをコピー</button>
  <span id="status"></span>
</header>
<main>
  <div id="warnings"></div>
  <section id="grid"></section>
  <textarea id="fallback" hidden readonly></textarea>
</main>
<dialog id="zoom" aria-labelledby="zoom-title">
  <div class="zoom-head"><h2 id="zoom-title"></h2><button id="zoom-close" type="button">閉じる（Esc）</button></div>
  <canvas id="zoom-eic"></canvas>
  <canvas id="zoom-mirror"></canvas>
</dialog>
<script>
const EMBEDDED = /*__CURATION_DATA__*/null;
const PREFIX = "__SUBMISSION_PREFIX__";
let DATA = EMBEDDED;
const choices = new Map();   // spot_id -> {candidate?, flag, level?, note}

/*__VIEWER_COMMON__*/

// --- suggest (pure) ---
// 初期選択は強い説明のある ④ の行だけ（spec §6.2）。それ以外は「選ばない」。
function initialChoice(spot) { return spot.preset ? {candidate: spot.preset, flag: "redundant", level: "sum"} : null; }
// 送信用テキストには候補 ID だけを入れる（中身は curation_submit が保存済みの候補から展開する）。
function choiceEntry(spot, choice) {
  const entry = {spot_id: spot.spot_id, flag: choice.flag};
  if (choice.flag === "assign") { entry.candidate = choice.candidate; entry.level = choice.level || "sum"; }
  if (choice.flag === "redundant") entry.candidate = choice.candidate;
  entry.note = choice.note || "";
  return entry;
}
function suggestSubmission(suggestion, choiceMap) {
  const byId = new Map((suggestion.spots || []).map(s => [s.spot_id, s]));
  const flags = [...choiceMap.entries()].map(([id, c]) => choiceEntry(byId.get(id) || {spot_id: id}, c));
  return {review_id: suggestion.suggestion_id, arf2_path: suggestion.arf2_path, flags};
}
// ④ の行: X のトレースに相手 Y の代表試料のトレースを足す。Y は X の代表試料の最大値に合わせて
// 縮める（M+2 は Y の 1 割ほどしかなく、そのままでは X が潰れて形を比べられないため）。
function overlayFor(spot, relation) {
  const samples = (spot.eic && spot.eic.samples) || [];
  const trace = relation && relation.partner_eic;
  if (!trace || !trace.points.length) return {eic: {samples}};
  const rep = samples.find(s => s.representative) || samples[0];
  const xMax = rep ? Math.max(...rep.points.map(p => p[1])) : 0;
  const yMax = Math.max(...trace.points.map(p => p[1])) || 1;
  const k = xMax > 0 ? xMax / yMax : 1;
  const partner = {...trace, detected: true, representative: false, partner: true,
                   points: trace.points.map(p => [p[0], p[1] * k])};
  return {eic: {samples: [...samples, partner]}};
}
function mirrorFor(spot, candidate) {
  if (!candidate || !candidate.mirror) return {mirror: null, mirror_empty_text: spot.measured && spot.measured.length ? "参照スペクトルなし" : "MS/MS なし"};
  return {mirror: {measured: spot.measured, reference: candidate.mirror.reference, matched_mz: candidate.mirror.matched_mz}};
}
// --- end suggest ---

function candidateLabel(c) {
  return `${c.sum_name || c.name} ${c.adduct || ""}`;
}

function spotCard(spot) {
  const card = el("article", {class: "card" + (spot.strong ? " strong" : ""), id: "spot-" + spot.spot_id});
  card.appendChild(el("h2", {}, spot.name && spot.target_kind !== "unannotated" ? spot.name : `#${spot.spot_id}（未注釈）`));
  const kind = {flagged: "wrong フラグ", likely_wrong: "likely_wrong", unannotated: "未注釈"}[spot.target_kind];
  card.appendChild(el("div", {class: "meta"},
    `#${spot.spot_id} ${spot.adduct || ""} m/z ${fixed(spot.mz, 4)} RT ${fixed(spot.rt, 2)} / ${kind}` +
    (spot.target_reasons.length ? ` / ${spot.target_reasons.join(", ")}` : "")));
  const eic = canvasFor(card), mirror = canvasFor(card);
  let focused = spot.candidates[0] || null;
  const draw = () => {
    const chosen = spot.relations.find(r => r.candidate_id === (choices.get(spot.spot_id) || {}).candidate);
    drawEic(eic, chosen ? {...spot, ...overlayFor(spot, chosen)} : spot);   // 初期選択の ④ は重ねて描く
    drawMirror(mirror, {...spot, ...mirrorFor(spot, focused)});
  };
  for (const prepare of [eic, mirror]) prepare.canvas.addEventListener("click", () => openZoom(spot, focused));
  const list = el("div");
  const radioName = "c-" + spot.spot_id;
  const current = () => choices.get(spot.spot_id) || null;
  const note = el("input", {type: "text", placeholder: "メモ"});
  const setChoice = choice => {
    if (choice) choices.set(spot.spot_id, {...choice, note: note.value}); else choices.delete(spot.spot_id);
    card.classList.toggle("chosen", !!choice);
    document.getElementById("status").textContent = statusText();
  };
  const addRow = (value, labelText, detail, onPick, extraClass) => {
    const row = el("label", {class: "cand" + (extraClass ? " " + extraClass : "")});
    const radio = el("input", {type: "radio", name: radioName, value});
    const c = current();
    radio.checked = c ? (c.candidate ? c.candidate === value : c.flag === value) : value === "none";
    radio.addEventListener("change", onPick);
    const text = el("div"); text.append(el("div", {}, labelText));
    if (detail) text.append(el("div", {class: "meta"}, detail));
    row.append(radio, text); list.appendChild(row); return row;
  };
  for (const r of spot.relations) {
    const detail = `Δm/z ${r.dmz_mda ?? "–"} mDa / ΔRT ${r.drt} / r ${r.r ?? "–"}` +
      (r.msdial_links.length ? ` / MS-DIAL: ${r.msdial_links.join(",")}` : "") +
      (r.isotope_ratio != null ? ` / 強度比 ${r.isotope_ratio}（期待 ${r.expected_ratio}）` : "") +
      (r.soft.length ? ` / ${r.soft.join(",")}` : "") + (r.strong ? " / 強い説明" : "");
    const label = `#${r.of}（${r.of_name || "?"}）の別イオン: ${r.relation}` + (r.informational ? "（情報のみ）" : "");
    const row = addRow(r.candidate_id, label, detail, () => {
      drawEic(eic, {...spot, ...overlayFor(spot, r)});
      setChoice({candidate: r.candidate_id, flag: "redundant"});
    }, r.informational ? "info" : "");
    row.addEventListener("mouseenter", () => drawEic(eic, {...spot, ...overlayFor(spot, r)}));
    if (r.informational) row.querySelector("input").disabled = true;
  }
  for (const c of spot.candidates) {
    const src = el("span", {class: "src"}, c.source);
    const detail = `score ${fixed(c.scores.total_score, 2)} / wdot ${fixed(c.scores.weighted_dot_product, 2)} / ` +
      `一致 ${c.scores.matched_peaks_count} / Δm/z ${c.dmz_mda ?? "–"} mDa` +
      (c.trend ? ` / 傾向 z ${c.trend.z}` : "") + ([...c.soft, ...c.info].length ? ` / ${[...c.soft, ...c.info].join(",")}` : "");
    const row = addRow(c.candidate_id, candidateLabel(c), detail, () => {
      focused = c; drawMirror(mirror, {...spot, ...mirrorFor(spot, c)});
      setChoice({candidate: c.candidate_id, flag: "assign", level: species.checked ? "species" : "sum"});
    });
    row.firstChild.nextSibling.firstChild.prepend(src, " ");
    row.addEventListener("mouseenter", () => { focused = c; drawMirror(mirror, {...spot, ...mirrorFor(spot, c)}); });
  }
  if (spot.target_kind === "flagged")
    addRow("clear", "元の注釈に戻す（wrong を取り消す）", null, () => setChoice({flag: "clear"}));
  addRow("none", "選ばない", null, () => setChoice(null));
  card.appendChild(list);
  const species = el("input", {type: "checkbox"});
  const speciesLabel = el("label", {class: "meta"}); speciesLabel.append(species, " 分子種で記録（既定は和組成）");
  species.addEventListener("change", () => { const c = current();
    if (c && c.flag === "assign") setChoice({...c, level: species.checked ? "species" : "sum"}); });
  note.addEventListener("change", () => { const c = current(); if (c) setChoice(c); });
  const foot = el("div", {class: "flags"}); foot.append(speciesLabel, note); card.appendChild(foot);
  if (current()) card.classList.add("chosen");
  card._draw = draw;
  return card;
}

function statusText() {
  const c = DATA.counts;
  return `対象 ${DATA.spots.length}（wrong ${c.targets.flagged} / likely_wrong ${c.targets.likely_wrong} / 未注釈 ${c.targets.unannotated}）` +
    ` / 強い説明 ${c.with_strong_relation} / 未送信の選択 ${choices.size}`;
}

function passes(spot) {
  const t = document.getElementById("f-target").value, s = document.getElementById("f-strong").value;
  const src = document.getElementById("f-source").value, ch = document.getElementById("f-chosen").value;
  if (t && spot.target_kind !== t) return false;
  if (s && !spot.strong) return false;
  if (src === "ion" && !spot.relations.length) return false;
  if ((src === "msdial" || src === "research") && !spot.candidates.some(c => c.source.includes(src))) return false;
  if (ch === "chosen" && !choices.has(spot.spot_id)) return false;
  if (ch === "open" && choices.has(spot.spot_id)) return false;
  return true;
}

const observer = new IntersectionObserver(entries => { for (const e of entries)
  if (e.isIntersecting && e.target._draw) { e.target._draw(); observer.unobserve(e.target); } }, {rootMargin: "200px"});

function render() {
  const grid = document.getElementById("grid"); grid.textContent = "";
  for (const spot of DATA.spots.filter(passes)) { const card = spotCard(spot); grid.appendChild(card); observer.observe(card); }
  document.getElementById("status").textContent = statusText();
}

const zoomEic = preparerFor(document.getElementById("zoom-eic"));
const zoomMirror = preparerFor(document.getElementById("zoom-mirror"));
function closeZoom() { const d = document.getElementById("zoom"); if (typeof d.close === "function") d.close(); else d.removeAttribute("open"); }
function openZoom(spot, candidate) {
  const d = document.getElementById("zoom");
  document.getElementById("zoom-title").textContent = `#${spot.spot_id} ${candidate ? candidateLabel(candidate) : ""}`;
  if (typeof d.showModal === "function") { if (!d.open) d.showModal(); } else d.setAttribute("open", "");
  drawEic(zoomEic, spot); drawMirror(zoomMirror, {...spot, ...mirrorFor(spot, candidate)});
}
document.getElementById("zoom-close").addEventListener("click", closeZoom);
document.addEventListener("keydown", e => { if (e.key === "Escape") closeZoom(); });

function start() {
  document.getElementById("title").textContent = `Annotation suggestions — ${DATA.alignment.alignment_file}`;
  const warnings = document.getElementById("warnings");
  for (const w of DATA.warnings || []) warnings.appendChild(el("div", {class: "warn"}, w));
  for (const s of DATA.spots) { const c = initialChoice(s); if (c && !choices.has(s.spot_id)) choices.set(s.spot_id, {...c, note: ""}); }
  for (const id of ["f-target", "f-strong", "f-source", "f-chosen"]) document.getElementById(id).addEventListener("change", render);
  render();
}

document.getElementById("copy").addEventListener("click", async () => {
  const text = PREFIX + " " + JSON.stringify(suggestSubmission(DATA, choices));
  const area = document.getElementById("fallback");
  try { await navigator.clipboard.writeText(text); document.getElementById("status").textContent = "コピーしました。チャットに貼って送ってください。"; }
  catch { area.hidden = false; area.value = text; area.select(); document.getElementById("status").textContent = "下の欄を選択してコピーしてください。"; }
});

if (DATA) start(); else document.getElementById("status").textContent = "データがありません。";
</script>
</body>
</html>
```

先頭は `viewer.html` と同じ `<!doctype html>` 〜 `<title>Annotation suggestions</title>` と `<style>`。spec §7.2 の「④ の行を選ぶと X と Y の EIC を重ねる」は、Task 6 が `relations[].partner_eic` に入れた Y の代表試料のトレースを `overlayFor` で X のトレースに足して描く（Step 4）。

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_suggest_viewer.py tests/test_curation_review.py tests/test_curation_suggest.py -q`
Expected: PASS

- [ ] **Step 6: 内蔵ブラウザで目視する**

`check.py` で `write_suggest_set` を tmp に書き、`run_suggestion` → `save_suggestion` した HTML のパスを出力し、Browser pane（`mcp__Claude_Browser__navigate` に `file:///` のパス）で開く。確かめること: spot 1 のカードが緑枠で「#0 の別イオン: isotope_M+2」が初期選択、その EIC に spot 0 のトレースが参照色の点線で重なっている、候補の行にマウスを載せるとミラーが替わる、「送信用テキストをコピー」の文に `review_id":"cs-` と `candidate` だけが入り名前が入っていない、携帯幅（`resize_window` の mobile）で横スクロールが出ない。終わったら `check.py` を空にする。

- [ ] **Step 7: コミット**

```bash
git add lipidmix/curation/viewer_common.js lipidmix/curation/viewer.html lipidmix/curation/viewer.py lipidmix/curation/suggest_viewer.html tests/test_curation_suggest_viewer.py
git commit -m "feat(curation): 候補付けのビューアを足し、描画を両ビューアで共有する"
```

---

### Task 8: MCP ツール（`curation_suggest` と送信・一覧の拡張）

**Files:**
- Modify: `lipidmix/tools/curation_tools.py`
- Modify: `tests/test_server_registration.py`、`tests/test_tool_annotations.py`、`tests/test_workflow_docs.py`
- Modify: `CLAUDE.md`、`USAGE.md`、`docs/workflow/index.md`、`docs/workflow/curation.md`（腐敗防止テストが縛る部分だけ。Step 4）
- Test: `tests/test_curation_tools.py`

**Interfaces:**
- Consumes: `suggest.*`（Task 6）、`viewer.render_suggest_html`（Task 7）。
- Produces: `curation_suggest(review_id=None, wrong="flagged_or_likely", unannotated=True, include_decided=False, top_n=5, rt_window=None, relation_mz_tol=None, relation_min_r=None, thresholds=None, file_path=None, max_rows=100) -> str`。`curation_submit` は `review_id` が `cs-` のとき `suggest.expand_entries` を使う。`curation_flags` の TSV 列は `spot_id, flag, name, of, note, source, ts`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_curation_tools.py` の既存の流儀（`curation_tools` の関数を直接呼び、`json.loads` で読む。`session_state.session.library.store` を差し替える）に合わせて足す:

```python
from lipidmix.curation import suggest
from tests.curation_fixtures import write_suggest_set


@pytest.fixture()
def suggest_env(tmp_path, monkeypatch):
    monkeypatch.setenv(library_store.LIBRARY_CACHE_ENV, str(tmp_path / "cache"))
    paths = write_suggest_set(tmp_path / "neg")
    s = library_store.open_store(paths["msp"])
    monkeypatch.setattr(session_state.session.library, "store", s)
    yield paths
    s.close()


def test_suggest_requires_a_review(suggest_env):
    out = json.loads(curation_tools.curation_suggest(file_path=str(suggest_env["arf2"])))
    assert out["error"]["code"] == "missing_state" and out["error"]["required_tools"] == ["curation_review"]


def test_suggest_requires_the_library(suggest_env, monkeypatch):
    monkeypatch.setattr(session_state.session.library, "store", None)
    out = json.loads(curation_tools.curation_suggest(file_path=str(suggest_env["arf2"])))
    assert out["error"]["required_tools"] == ["library_load"]


def test_suggest_then_submit_assign_and_redundant(suggest_env):
    arf2 = str(suggest_env["arf2"])
    json.loads(curation_tools.curation_review(file_path=arf2))
    out = json.loads(curation_tools.curation_suggest(file_path=arf2))
    sid = out["suggestion_id"]
    assert out["counts"]["targets"]["unannotated"] == 3
    assert out["table"].splitlines()[0].split("\t") == suggest.TSV_COLUMNS
    assert out["html_path"].endswith(f"suggest-{sid}.html")
    assert "spots" not in out                                     # 座標・EIC は戻り値に入れない
    text = 'CURATION_SUBMIT ' + json.dumps({"review_id": sid, "arf2_path": arf2, "flags": [
        {"spot_id": 2, "flag": "assign", "candidate": "L1", "level": "sum", "note": ""},
        {"spot_id": 3, "flag": "redundant", "candidate": "R1", "note": ""}]})
    submitted = json.loads(curation_tools.curation_submit(submission_text=text))
    assert submitted["status"] == "ok" and submitted["recorded"] == 2
    assert submitted["n_assign"] == 1 and submitted["n_redundant"] == 1
    assert submitted["tags_xml"]["added"] == [] and submitted["tags_xml"]["removed"] == []   # _tags.xml は不変
    table = json.loads(curation_tools.curation_flags(file_path=arf2))["table"].splitlines()
    assert table[0] == "spot_id\tflag\tname\tof\tnote\tsource\tts"
    assert table[1].split("\t")[:4] == ["2", "assign", "PE 36:2", ""]
    assert table[2].split("\t")[:4] == ["3", "redundant", "", "0"]


def test_submit_rejects_a_forged_candidate(suggest_env):
    arf2 = str(suggest_env["arf2"])
    json.loads(curation_tools.curation_review(file_path=arf2))
    sid = json.loads(curation_tools.curation_suggest(file_path=arf2))["suggestion_id"]
    out = json.loads(curation_tools.curation_submit(review_id=sid, file_path=arf2, flags=[
        {"spot_id": 2, "flag": "assign", "candidate": "L7", "level": "sum"}]))
    assert out["status"] == "error"
    flags_path = suggest_env["arf2"].parent / "curation" / "flags.jsonl"
    assert not flags_path.exists() or '"assign"' not in flags_path.read_text(encoding="utf-8")
```

（ファイル冒頭の既存 import に `library_store` / `session_state` / `curation_tools` があれば再 import しない。）

`tests/test_server_registration.py`: `EXPECTED_TOOLS` に `"curation_suggest"` を足し、コメントを 1 行足す: `# キュレーションの候補付け: curation_suggest を追加(71→72)。`。登録数を数える箇所（`len(EXPECTED_TOOLS)` と比べる箇所）は数を直書きしていなければ触らない。

`tests/test_tool_annotations.py` の表に `"curation_suggest": LOCAL_WRITE_APPEND,` を足す。

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_tools.py tests/test_server_registration.py tests/test_tool_annotations.py -q`
Expected: FAIL（`curation_suggest` が無い）

- [ ] **Step 3: 実装を書く**

`lipidmix/tools/curation_tools.py`:

`__all__` に `"curation_suggest"` を足す。import に `suggest` を足す（`from lipidmix.curation import evidence, judge, msdial_writeback, review, suggest, viewer`）。

`_bad_review_id` と `_find_review` を、レビューと候補付けの両方を扱う形に直す:

```python
def _bad_review_id(review_id) -> str | None:
    if review.is_valid_review_id(review_id) or suggest.is_valid_suggestion_id(review_id):
        return None
    return _error(f"review_id={review_id!r} の形が不正です（cr-… か cs-…）。"
                  "ビューアの送信用テキストをそのまま貼ってください。")


def _load_any(directory, review_id: str) -> dict:
    if suggest.is_valid_suggestion_id(review_id):
        return suggest.load_suggestion(directory, review_id)
    return review.load_review(directory, review_id)


def _find_review(review_id: str, *arf2_paths) -> tuple[dict | None, list[Path]]:
    searched = _candidate_dirs(review_id, *arf2_paths)
    for directory in searched:
        try:
            return _load_any(directory, review_id), searched
        except FileNotFoundError:
            continue
    return None, searched
```

`curation_view_data` は `review.page(saved, page)` を呼ぶので、候補付けの ID が来たら `_error("curation_view_data はレビュー（cr-…）専用です。")` を返す（先頭で `if suggest.is_valid_suggestion_id(review_id): return _error(...)`）。

`curation_submit` の検証部分（`try: cleaned = flag_log.validate_entries(...)`）を置き換える:

```python
    try:
        if suggest.is_valid_suggestion_id(review_id):
            cleaned = suggest.expand_entries(entries, saved)
        else:
            cleaned = flag_log.validate_entries(
                entries, allowed_spot_ids={s["spot_id"] for s in saved["spots"]})
    except ValueError as exc:
        return _error(str(exc))
```

sha256 不一致のメッセージは「curation_review（候補付けなら curation_suggest）をやり直してください。」にする。戻り値に `"n_assign"` と `"n_redundant"` を足す（`effective` から数える）。docstring に次を足す: 「候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign|redundant|clear, candidate: "L1"|"R1", level: sum|species, note}]`。候補の中身（名前・InChIKey）は保存済みの候補付けから展開するので、送信側で名前を書かない。assign / redundant は `_tags.xml` を変えない。」

`curation_flags` の表を 7 列にする:

```python
    lines = ["spot_id\tflag\tname\tof\tnote\tsource\tts"]
    for spot_id in sorted(effective):
        row = effective[spot_id]
        lines.append("\t".join([str(spot_id), row["flag"], str(row.get("name") or ""),
                                "" if row.get("of") is None else str(row["of"]),
                                (row.get("note") or "").replace("\t", " "),
                                row.get("source", ""), row.get("ts", "")]))
```

新しいツール:

```python
@mcp.tool(annotations=_LOCAL_WRITE_APPEND, structured_output=False)
def curation_suggest(review_id: str | None = None, wrong: str = "flagged_or_likely",
                     unannotated: bool = True, include_decided: bool = False, top_n: int = 5,
                     rt_window: float | None = None, relation_mz_tol: float | None = None,
                     relation_min_r: float | None = None, thresholds: dict | None = None,
                     file_path: str | None = None, max_rows: int = 100) -> str:
    """キュレーションで「間違い」になったスポットと未注釈スポットに、注釈の候補を並べる
    （候補付けレビュー）。**先に library_load と curation_review**。

    対象: `wrong="flagged_or_likely"`（既定）= wrong フラグ ＋ 元レビューの likely_wrong、
    `"flagged"` = wrong フラグだけ。`unannotated=True` で未注釈スポットも。判断済み
    （assign / redundant）のスポットは `include_decided=True` のときだけ含める。
    元レビューは `review_id`（省略時はこのアラインメントの最新のレビュー）。
    候補: ① MS-DIAL の下位候補、② 閾値を緩めた再検索（スコアの足切りなし）、
    ④ 注釈付きの別スポットの同位体・アダクト・インソース断片としての説明。
    スペクトル類似度で並べ、極性矛盾・|Δm/z| ≥ 10 mDa は削り、非典型アダクト・RT–m/z 傾向の外れ・
    一致ピーク 0 は順位を下げる。MS-DIAL の脂質規則は評価していない（鎖組成は保証しない）ので、
    既定では和組成で記録する。
    戻り値は 1 スポット 1 行の TSV（強い説明のあるものが先、先頭 `max_rows` 行）・件数・`html_path`。
    ユーザーには `html_path` をブラウザで開いてもらい、選んだ内容を「送信用テキストをコピー」で
    チャットに貼ってもらう。貼られたら curation_submit(submission_text=...) に渡す。
    **ユーザーの同意なしに curation_submit を呼ばない。**
    """
    arf2_path = resolve_arf2_file_path(file_path)
    if not arf2_path:
        return mcp_errors.missing_state(
            "arf2_file", ["load_dataset", "arf2_parser"],
            ".arf2 が見つかりません。先に load_dataset で MS-DIAL の出力フォルダを指定してください。")
    store = session_state.session.library.store
    if store is None:
        return mcp_errors.missing_state(
            "library", ["library_load"],
            "参照ライブラリが読み込まれていません。先に library_load を実行してください"
            "（アラインメントと同じフォルダの *_Loaded.msp2.dbs を推奨）。")
    try:
        th = judge.resolve_thresholds(thresholds)
    except ValueError as exc:
        return _error(str(exc))
    if review_id is not None:
        bad = None if review.is_valid_review_id(review_id) else _error(
            f"review_id={review_id!r} の形が不正です（cr-YYYYMMDD-HHMMSS-xxxx）。")
        if bad:
            return bad
        base, searched = _find_review(review_id, arf2_path)
        if base is None:
            return _not_found(review_id, searched)
    else:
        base = suggest.latest_review(arf2_path, flag_log.alignment_key(arf2_path)["alignment_sha256"])
        if base is None:
            return mcp_errors.missing_state(
                "curation_review", ["curation_review"],
                "このアラインメントのレビューがありません。先に curation_review を実行してください"
                "（likely_wrong の判定をそこから読みます）。")
    options = {"wrong": wrong, "unannotated": unannotated, "include_decided": include_decided,
               "top_n": top_n, "rt_window": rt_window, "relation_mz_tol": relation_mz_tol,
               "relation_min_r": relation_min_r}
    try:
        result = suggest.run_suggestion(arf2_path, base_review=base, store=store, th=th, options=options)
    except flag_log.FlagFileError as exc:
        return _flag_file_error(exc)
    except ValueError as exc:
        return _error(str(exc))
    saved = suggest.save_suggestion(result)
    session_state.session.curation.review_dirs[result["suggestion_id"]] = str(saved["json"].parent)
    table = suggest.summary_tsv(result, max_rows=max_rows)
    return json_payload(round_floats({
        "suggestion_id": result["suggestion_id"], "base_review_id": result["base_review_id"],
        "n_spots": len(result["spots"]), "counts": result["counts"], "warnings": result["warnings"],
        "table": table, "n_table_rows_shown": len(table.splitlines()) - 1,
        "html_path": str(saved["html"]), "library": result["library"],
        "analysis_params": {k: result["analysis_params"][k] for k in ("source", "path", "rt_window")},
        "options": result["options"], "thresholds": th}, 4))
```

`_find_review(review_id, arf2_path)` の第 2 引数は `.arf2` のパス（`_candidate_dirs` が `curation_dir` にする）。

- [ ] **Step 4: 腐敗防止テストが縛る文書を同じタスクで直す**

ツールを 1 本足すと、pre-commit の全テストのうち `test_readme_links.py`（USAGE.md のツール集合と件数、CLAUDE.md 冒頭の規模表記）と `test_workflow_docs.py`（対象ツール数）が落ちる。コミットを通すため、公開面に縛られた次の文書をこのタスクで直す（値の意味の詳細は Task 10）。

- `CLAUDE.md` 冒頭: 「ツール 71」を「ツール 72」に（ほかの数字は変えない）。
- `USAGE.md` のキュレーションの表に `curation_suggest` の行を足し、ツール件数の宣言（`test_readme_links.py` が読む行）を 1 増やす。行の文面:
  「`curation_suggest` | キュレーションで間違いになったスポット（既定は wrong フラグ＋元レビューの likely_wrong、`wrong="flagged"` でフラグだけ）と未注釈スポット（`unannotated`）に注釈の候補を並べる。**先に `library_load` と `curation_review`**（元レビューは `review_id`、省略時は最新）。候補は MS-DIAL の下位候補・閾値を緩めた再検索・別スポットの同位体／アダクト／インソース断片としての説明。スペクトル類似度で並べ、極性矛盾と \|Δm/z\| ≥ 10 mDa は削る。脂質規則は評価しないので既定は和組成で記録。戻り値は 1 スポット 1 行の TSV・件数・`html_path`。ビューアで選んで送信用テキストを貼る → `curation_submit`。|」
  `curation_submit` の行に「候補付け（`cs-…`）の送信は `flags=[{spot_id, flag: assign\|redundant\|clear, candidate, level, note}]`。中身は保存済みの候補から展開する。assign / redundant は `_tags.xml` を変えない。」を足す。`curation_flags` の行を「現在のアラインメントで有効な判断を TSV で返す（列は spot_id・flag・name（assign の記録名）・of（redundant の相手）・note・source・ts）。」にする。
- `docs/workflow/index.md`: curation.md の行の説明を「アラインメントのキュレーション — 注釈の一覧確認・機械判別・フラグ記録・候補付け」、対象ツール数を 5 に、合計行（「合計 50 ツール。対象外は 21 ツール」）を「合計 51 ツール。対象外は 21 ツール」に。
- `tests/test_workflow_docs.py` の `"curation.md": (...)` に `"curation_suggest"` を足す。
- `docs/workflow/curation.md`: mermaid 図に次を足し、`## curation_suggest` の節を `## curation_submit` の前に足す（行番号は書かない規約。挙げる関数はすべて Task 1〜7 で実在する）:

```
    CG[curation_suggest] --> LR[curation.suggest.latest_review]
    CG --> RS[curation.suggest.run_suggestion]
    RS --> TG[curation.suggest.select_targets]
    RS --> IF[arf2.ion_features.load_ion_features]
    RS --> AP[msdial.analysis_params.resolve_analysis_params]
    RS --> COL
    RS --> LC[curation.candidates.build_library_candidates]
    RS --> FR[curation.relations.find_relations]
    CG --> SS[curation.suggest.save_suggestion]
    SS --> SH[curation.viewer.render_suggest_html]
    CS --> EX[curation.suggest.expand_entries]
```

  節の本文（この形で書く）:

```markdown
## curation_suggest

1. `curation.suggest.latest_review` が、同じ sha256 のアラインメントに対する最新のレビューを探す（`review_id` を渡せばそれ）。無ければ `missing_state`。
2. `curation.suggest.run_suggestion` が有効な判断（`curation.flags.FlagStore.effective`）を読み、`curation.suggest.select_targets` で対象（wrong フラグ・likely_wrong・未注釈）を選ぶ。
3. 相手 Y の候補（注釈付きで、対象でも wrong / redundant でもないスポット）から `curation.trend.fit_trends` で傾向を当てはめる。
4. `curation.evidence.collect`（`keep_measured=True`）で対象の証拠を集め、`curation.candidates.build_library_candidates`（①②）と `curation.relations.find_relations`（④、`arf2.ion_features.load_ion_features` のリンクと `msdial.analysis_params.resolve_analysis_params` の検索アダクト）を回す。
5. `curation.suggest.save_suggestion` が `suggest-<id>.json` と `curation.viewer.render_suggest_html` の HTML を書く。

`curation_submit` は `review_id` が `cs-…` のとき、`curation.flags.validate_entries` の代わりに `curation.suggest.expand_entries` で候補 ID を記録行へ展開する。
```

- [ ] **Step 5: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS

- [ ] **Step 6: コミット（バックグラウンド、出力をファイルへ）**

```bash
git add lipidmix/tools/curation_tools.py tests/test_curation_tools.py tests/test_server_registration.py tests/test_tool_annotations.py tests/test_workflow_docs.py CLAUDE.md USAGE.md docs/workflow/index.md docs/workflow/curation.md
git commit -m "feat(curation): curation_suggest を足し、送信と一覧を候補付けに広げる"
```

---

### Task 9: 差次的エクスポートと同定一覧への反映

**Files:**
- Modify: `lipidmix/curation/apply.py`
- Modify: `lipidmix/arf/tools.py`（`arf_export_differential`）
- Modify: `lipidmix/analysis/dataset_export.py`、`lipidmix/tools/dataset_analysis_tools.py`（必要なら）
- Modify: `lipidmix/arf2/tools.py`（`arf2_annotate_identities`）
- Test: `tests/test_curation_apply.py`、既存のエクスポートのテスト（`tests/test_export_contract.py` など、`apply_curation` を試している既存テストを `grep -rln apply_curation tests` で探す）

**Interfaces:**
- Consumes: `flags.split_decisions`。
- Produces:
  - `flags_for_arf2(path)` の戻り値に `"assign": {spot: row}` と `"redundant": set` を足す（既存キーは同じ）
  - `override_identity(catalog: dict[int, dict], flag_set) -> dict[int, dict]`（assign のスポットの `Name` / `Ontology` / `InChIKey` を置き換え、`_curation="assign"` を付けた新しい dict）
  - `identity_for(spot_id, flag_set) -> dict | None`（`{"name", "ontology", "inchikey"}`）
  - `filter_rows(rows, flag_set, *, key)` の stats に `"redundant_excluded"` と `"assigned"` を足す（内部の値。既存テスト `tests/test_curation_apply.py::test_wrong_rows_are_dropped_and_suspect_kept` の `stats ==` の期待だけを 4 キーに直す）
  - `meta_line(...)`: assign / redundant が 1 件以上のときだけ `curation_assigned = N` と `curation_redundant_excluded = N` を `curation_suspect` の後に足す
  - `payload_summary(...)`: assign / redundant が 1 件以上のときだけ `"assigned"` と `"redundant_excluded"` を足す（wrong / suspect だけの payload は現行と同じ形のまま。`tests/test_export_differential.py` と `tests/test_dataset_analysis_tools.py` の `payload["curation"] == {...}` の既存の期待は無修正で通ること）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_curation_apply.py` に足す:

```python
from lipidmix.curation import apply as curation_apply


def _flag_set(**kw):
    base = {"wrong": set(), "suspect": set(), "assign": {}, "redundant": set(), "digest": "d", "n": 0,
            "orphaned": 0}
    base.update(kw)
    base["n"] = len(base["wrong"]) + len(base["suspect"]) + len(base["assign"]) + len(base["redundant"])
    return base


def test_meta_line_unchanged_without_decisions():
    fs = _flag_set(wrong={1}, suspect={2})
    line = curation_apply.meta_line("applied", fs, {"wrong_excluded": 1, "suspect": 1,
                                                    "redundant_excluded": 0, "assigned": 0})
    assert line == ("# curation = applied\tcuration_flags = 2\tcuration_wrong_excluded = 1\t"
                    "curation_suspect = 1\tcuration_flags_sha256 = d")


def test_meta_line_with_decisions():
    fs = _flag_set(assign={3: {"name": "PC 34:1", "ontology": "PC", "inchikey": "K"}}, redundant={4})
    line = curation_apply.meta_line("applied", fs, {"wrong_excluded": 0, "suspect": 0,
                                                    "redundant_excluded": 1, "assigned": 1})
    assert "curation_assigned = 1\tcuration_redundant_excluded = 1\tcuration_flags_sha256 = d" in line


def test_override_identity_makes_an_unannotated_spot_exportable():
    catalog = {3: {"MasterAlignmentID": 3, "Name": "Unknown", "Ontology": "", "InChIKey": ""}}
    fs = _flag_set(assign={3: {"name": "PE 36:2", "ontology": "PE", "inchikey": "KEY-PE362"}})
    out = curation_apply.override_identity(catalog, fs)
    assert out[3]["Name"] == "PE 36:2" and out[3]["InChIKey"] == "KEY-PE362" and out[3]["_curation"] == "assign"
    assert catalog[3]["Name"] == "Unknown"                       # 元の dict は変えない


def test_filter_rows_excludes_wrong_and_redundant_and_counts_assigned():
    rows = [{"spot_id": i} for i in range(5)]
    fs = _flag_set(wrong={0}, redundant={1}, suspect={2}, assign={3: {"name": "x", "ontology": "", "inchikey": "k"}})
    kept, stats = curation_apply.filter_rows(rows, fs, key=lambda r: r["spot_id"])
    assert [r["spot_id"] for r in kept] == [2, 3, 4]
    assert stats == {"wrong_excluded": 1, "redundant_excluded": 1, "suspect": 1, "assigned": 1}
```

`tests/test_curation_apply.py` の既存テスト `test_wrong_rows_are_dropped_and_suspect_kept` の期待を
`assert stats == {"wrong_excluded": 1, "redundant_excluded": 0, "suspect": 1, "assigned": 0}` に直す（stats は内部の値）。

**「wrong / suspect だけなら現行と同じ」の証拠（Review Focus 5）**: 新しいテストは書かない。既存の
`tests/test_export_differential.py` の `test_no_curation_flags_means_no_curation_meta_line`・
`test_wrong_flag_drops_the_row_and_declares_curation_applied`・`test_orphaned_flags_are_reported_and_their_rows_are_kept`・
`test_apply_curation_false_declares_not_applied_and_keeps_the_rows`・`test_applied_curation_is_summarised_in_the_payload`、
`tests/test_dataset_analysis_tools.py::test_dataset_export_applies_curation_flags_via_sibling_arf2`、
`tests/test_parser_decode.py` の `test_arf2_annotate_identities_shows_the_curation_flag` が**無修正で**通ることを証拠にする。
これらの期待値を書き換えたくなったら、それは回帰なので実装を直す。

ARF 経路（`tests/test_export_differential.py` の `TestExportDifferential` に足す。setUp の catalog は spot 1 = PC 34:1（InChIKey あり）、spot 2 = Unknown（InChIKey なし））:

```python
    def test_assign_adds_an_unannotated_spot_and_redundant_drops_one(self):
        from lipidmix.curation import flags as curation_flags

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "differential.tsv"
            arf2 = Path(tmp) / "AlignmentResult_2026_01_01_00_00_00.arf2"
            arf2.write_bytes(b"")
            curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
                [{"spot_id": 2, "flag": "assign", "name": "PE 36:2", "level": "sum", "ontology": "PE",
                  "inchikey": "PEPEPEPEPEPEPE-XXXXXXXXXX-N"},
                 {"spot_id": 1, "flag": "redundant", "of": 9, "relation": "isotope_M+1"}],
                alignment=curation_flags.alignment_key(arf2), review_id="cs-x", source="user")
            with patch("lipidmix.arf.tools._sibling_arf2_path", return_value=arf2), \
                 patch("lipidmix.arf2.reader.load_catalog", return_value=self.catalog):
                payload = json.loads(server.arf_export_differential(str(out)))
            text = out.read_text(encoding="utf-8")
        self.assertEqual(payload["status"], "success")
        body = [l for l in text.splitlines() if not l.startswith("#")]
        header = body[0].split("\t")
        rows = [dict(zip(header, l.split("\t"))) for l in body[1:]]
        self.assertEqual([r["name"] for r in rows], ["PE 36:2"])
        self.assertEqual((rows[0]["name_source"], rows[0]["inchikey_source"], rows[0]["inchikey"]),
                         ("curation", "curation", "PEPEPEPEPEPEPE-XXXXXXXXXX-N"))
        meta = next(l for l in text.splitlines() if l.startswith("# curation"))
        self.assertIn("curation_assigned = 1\tcuration_redundant_excluded = 1", meta)
        self.assertEqual(payload["curation"]["assigned"], 1)
        self.assertEqual(payload["curation"]["redundant_excluded"], 1)
```

mzTab-M 経路（`tests/test_dataset_analysis_tools.py` に足す。既存の `test_dataset_export_applies_curation_flags_via_sibling_arf2` と同じ準備で、記録する判断だけが違う）:

```python
def test_dataset_export_applies_assign_and_redundant(tmp_path):
    from lipidmix.curation import flags as curation_flags
    from lipidmix.tools.dataset_analysis_tools import (
        dataset_differential, dataset_export_differential, dataset_preprocess,
    )

    ds = _load_ds()
    ds.feature_ids = [str(i) for i in range(20)]
    ds.feature_metadata = {
        str(i): {"name": f"Compound {i}", "mz": 100.0 + i, "rt": 1.0 + i * 0.1,
                 "inchikey": f"AAAAAAAAAAAAAA-BBBBBBBBFB-{i % 10}",
                 "inchikey_source": "database_identifier"}
        for i in range(20)
    }
    mztab = tmp_path / "Height_AlignmentResult_2026_01_01_2026_01_01_09.mzTab"
    mztab.write_text("", encoding="utf-8")
    ds.source_files = {str(mztab): "sha"}
    ds.feature_qc = {"source": "arf"}
    arf2 = tmp_path / "AlignmentResult_2026_01_01.arf2"
    arf2.write_bytes(b"x")
    curation_flags.FlagStore(curation_flags.curation_dir(arf2)).append(
        [{"spot_id": 3, "flag": "assign", "name": "PE 36:2", "level": "sum", "ontology": "PE",
          "inchikey": "PEPEPEPEPEPEPE-XXXXXXXXXX-N"},
         {"spot_id": 4, "flag": "redundant", "of": 3, "relation": "isotope_M+1"}],
        alignment=curation_flags.alignment_key(arf2), review_id="cs-x", source="user")

    dataset_preprocess()
    a, b = _groups(session_state.session.dataset)
    dataset_differential(group_a=a, group_b=b)
    out = tmp_path / "diff.tsv"
    parsed = json.loads(dataset_export_differential(str(out)))
    assert parsed["status"] == "success"
    assert parsed["curation"]["assigned"] == 1 and parsed["curation"]["redundant_excluded"] == 1
    lines = out.read_text(encoding="utf-8").splitlines()
    body = [l.split("\t") for l in lines if l and not l.startswith("#")]
    header, rows = body[0], [dict(zip(body[0], r)) for r in body[1:]]
    by_id = {r["spot_id"]: r for r in rows}
    assert "4" not in by_id
    assert (by_id["3"]["name"], by_id["3"]["name_source"], by_id["3"]["inchikey"]) == \
        ("PE 36:2", "curation", "PEPEPEPEPEPEPE-XXXXXXXXXX-N")
```

（`_load_ds` がセッションに ds を置いているかは既存テストの冒頭で確かめ、置いていなければ既存テストと同じ 1 行 `session_state.session.dataset = ds` を足す。）

`arf2_annotate_identities`（`tests/test_parser_decode.py` の既存 `test_arf2_annotate_identities_shows_the_curation_flag` の隣に足す。`arf2_spot()` は MasterAlignmentID 7 のスポット）:

```python
    def test_arf2_annotate_identities_shows_assign(self):
        from lipidmix.arf2.tools import arf2_annotate_identities
        from lipidmix.curation import flags as curation_flags

        path = Path(self.write(
            "synthetic.arf2", lz4_container(msgpack_stream([arf2_spot()]))))
        curation_flags.FlagStore(curation_flags.curation_dir(path)).append(
            [{"spot_id": 7, "flag": "assign", "name": "PC 34:2", "level": "sum"}],
            alignment=curation_flags.alignment_key(path), review_id="cs-x", source="user")

        result = arf2_annotate_identities(file_path=str(path), max_rows=10)

        header, *body = result.splitlines()[-2:]
        column = header.split("\t").index("curation_flag")
        self.assertEqual(body[0].split("\t")[column], "assign:PC 34:2")
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_curation_apply.py tests/test_export_differential.py tests/test_dataset_analysis_tools.py tests/test_parser_decode.py -q`
Expected: FAIL（`override_identity` が無い、stats のキーが足りない、新しい 3 つのテスト）

- [ ] **Step 3: `apply.py` を書く**

```python
def flags_for_arf2(arf2_path) -> dict:
    store = FlagStore(curation_dir(arf2_path))
    if not store.exists():
        return {"wrong": set(), "suspect": set(), "assign": {}, "redundant": set(),
                "digest": flags_digest({}), "n": 0, "orphaned": 0}
    rows = store.rows()
    key = alignment_key(arf2_path)
    effective = effective_flags(rows, key["alignment_sha256"])
    decisions = split_decisions(effective)
    return {"wrong": decisions["wrong"], "suspect": decisions["suspect"],
            "assign": decisions["assign"], "redundant": set(decisions["redundant"]),
            "digest": flags_digest(effective), "n": len(effective),
            "orphaned": orphaned_count(rows, key)}


def identity_for(spot_id, flag_set) -> dict | None:
    row = (flag_set or {}).get("assign", {}).get(spot_id)
    if row is None:
        return None
    return {"name": row.get("name") or "", "ontology": row.get("ontology") or "",
            "inchikey": row.get("inchikey") or ""}


def override_identity(catalog: dict, flag_set) -> dict:
    out = dict(catalog)
    for spot_id in (flag_set or {}).get("assign", {}):
        identity = identity_for(spot_id, flag_set)
        if spot_id in out and identity is not None:
            out[spot_id] = {**out[spot_id], "Name": identity["name"], "Ontology": identity["ontology"],
                            "InChIKey": identity["inchikey"], "_curation": "assign"}
    return out


def filter_rows(rows, flag_set: dict, *, key):
    kept, stats = [], {"wrong_excluded": 0, "redundant_excluded": 0, "suspect": 0, "assigned": 0}
    for row in rows:
        spot = key(row)
        if spot in flag_set["wrong"]:
            stats["wrong_excluded"] += 1
            continue
        if spot in flag_set.get("redundant", ()):
            stats["redundant_excluded"] += 1
            continue
        if spot in flag_set["suspect"]:
            stats["suspect"] += 1
        if spot in flag_set.get("assign", {}):
            stats["assigned"] += 1
        kept.append(row)
    return kept, stats


def meta_line(state: str, flag_set: dict | None, stats: dict | None) -> str | None:
    if not flag_set or flag_set["n"] == 0:
        return None
    parts = [f"# curation = {state}", f"curation_flags = {flag_set['n']}"]
    if stats is not None:
        parts += [f"curation_wrong_excluded = {stats['wrong_excluded']}",
                  f"curation_suspect = {stats['suspect']}"]
        if flag_set.get("assign") or flag_set.get("redundant"):
            parts += [f"curation_assigned = {stats.get('assigned', 0)}",
                      f"curation_redundant_excluded = {stats.get('redundant_excluded', 0)}"]
    parts.append(f"curation_flags_sha256 = {flag_set['digest']}")
    return "\t".join(parts)


def payload_summary(state, flag_set, stats) -> dict:
    """エクスポートの成功 payload の `curation`。assign / redundant が無ければ現行と同じ 4 キー。"""
    summary = {"state": state,
               "wrong_excluded": stats["wrong_excluded"] if stats else None,
               "suspect": stats["suspect"] if stats else None,
               "orphaned": (flag_set or {}).get("orphaned", 0)}
    if (flag_set or {}).get("assign") or (flag_set or {}).get("redundant"):
        summary["assigned"] = stats.get("assigned") if stats else None
        summary["redundant_excluded"] = stats.get("redundant_excluded") if stats else None
    return summary
```

import に `split_decisions` を足す。モジュール docstring の「wrong のスポットは同定なしとして扱う」の後に「redundant は除外、assign は同定を置き換える（未注釈だったスポットは InChIKey を得てエクスポートに現れる）」を足す。

- [ ] **Step 4: ARF 経路（`lipidmix/arf/tools.py`）**

`arf_export_differential` の順序を「フラグを読む → 置き換えた catalog で結合 → 除外」にする。現行の

```python
    catalog = {spot.get("MasterAlignmentID"): spot for spot in load_catalog(arf2_path)}
    rows, report = identity_join.join_identity(last.get("results") or [], catalog)
    ... flag_set = curation_apply.flags_for_arf2(arf2_path) ...
    if apply_curation and flag_set["n"]:
        rows, curation_stats = curation_apply.filter_rows(...)
        report = {**report, "n_with_inchikey": len(rows),
                  "n_unannotated": report["n_unannotated"] + curation_stats["wrong_excluded"]}
```

を

```python
    catalog = {spot.get("MasterAlignmentID"): spot for spot in load_catalog(arf2_path)}
    from lipidmix.curation import apply as curation_apply
    from lipidmix.curation.flags import FlagFileError, orphaned_warning
    try:
        flag_set = curation_apply.flags_for_arf2(arf2_path)
    except FlagFileError as exc:
        return json_payload({"status": "error", "message": str(exc), **exc.details()})
    if apply_curation and flag_set["n"]:
        catalog = curation_apply.override_identity(catalog, flag_set)
    rows, report = identity_join.join_identity(last.get("results") or [], catalog)
    curation_stats = None
    if apply_curation and flag_set["n"]:
        rows, curation_stats = curation_apply.filter_rows(rows, flag_set, key=lambda r: r["spot_id"])
        report = {**report, "n_with_inchikey": len(rows),
                  "n_unannotated": report["n_unannotated"] + curation_stats["wrong_excluded"]
                  + curation_stats["redundant_excluded"]}
```

にする。行を書く箇所の `"name_source": "arf2"` と `"inchikey_source": "arf2"` を次にする:

```python
        by_curation = apply_curation and row["spot_id"] in flag_set["assign"]
        ...
            "name_source": "curation" if by_curation else "arf2",
            "inchikey_source": "curation" if by_curation else "arf2",
```

docstring の `apply_curation` の段落に「assign は同定を置き換え（name_source / inchikey_source = curation）、redundant は除外する」を足す。

- [ ] **Step 5: mzTab-M 経路（`lipidmix/analysis/dataset_export.py`）**

行を作るループの先頭（`inchikey = str(meta.get("inchikey") or "").strip()` の前）に足す:

```python
        override = None
        if curation is not None and curation["state"] == "applied":
            from lipidmix.curation.apply import identity_for
            override = identity_for(_as_spot_id(fid), curation["flag_set"])
```

`_as_spot_id` はこのファイルで `str(r["spot_id"])` と比べている流儀（SMF_ID は文字列）に合わせ、`flag_set["assign"]` のキー（int）と突き合わせるための小関数としてファイル内に置く:

```python
def _as_spot_id(fid):
    try:
        return int(fid)
    except (TypeError, ValueError):
        return fid
```

`override` があるときは、その行の `inchikey` / `name` / `name_source="curation"` / `inchikey_source="curation"` / `ontology` を `override` から取り、SME / SML の分岐を飛ばす。除外は現行の `wrong` の比較に `redundant` を足す:

```python
            wrong = {str(s) for s in flag_set["wrong"]}
            redundant = {str(s) for s in flag_set.get("redundant", ())}
            suspect = {str(s) for s in flag_set["suspect"]}
            before = len(rows)
            rows_kept = [r for r in rows if str(r["spot_id"]) not in wrong]
            n_wrong = before - len(rows_kept)
            rows = [r for r in rows_kept if str(r["spot_id"]) not in redundant]
            n_redundant = len(rows_kept) - len(rows)
            n_unannotated += n_wrong + n_redundant
            stats = {"wrong_excluded": n_wrong, "redundant_excluded": n_redundant,
                     "suspect": sum(1 for r in rows if str(r["spot_id"]) in suspect),
                     "assigned": sum(1 for r in rows if r.get("name_source") == "curation")}
```

- [ ] **Step 6: `arf2_annotate_identities`（`lipidmix/arf2/tools.py`）**

壊れた記録のときの既定を `{"wrong": set(), "suspect": set(), "assign": {}, "redundant": set()}` にし、`curation_flag` の決め方を:

```python
        if spot_id in flag_set["wrong"]:
            curation_flag = "wrong"
        elif spot_id in flag_set["suspect"]:
            curation_flag = "suspect"
        elif spot_id in flag_set.get("assign", {}):
            curation_flag = f"assign:{flag_set['assign'][spot_id].get('name') or ''}"
        elif spot_id in flag_set.get("redundant", ()):
            curation_flag = "redundant"
        else:
            curation_flag = ""
```

にする（`redundant` の相手の番号は `flags_for_arf2` が set しか返さないので出さない。テストの期待も `"redundant"` にする）。

- [ ] **Step 7: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: 全件 PASS。既存の wrong 除外・メタ行・payload のテストが**無修正で**通ること（Review Focus 5）。

- [ ] **Step 8: コミット（バックグラウンド、出力をファイルへ）**

```bash
git add lipidmix/curation/apply.py lipidmix/arf/tools.py lipidmix/analysis/dataset_export.py lipidmix/arf2/tools.py tests/test_curation_apply.py tests/test_export_differential.py tests/test_dataset_analysis_tools.py tests/test_parser_decode.py
git commit -m "feat(curation): assign を同定の置き換え、redundant を除外としてエクスポートへ反映する"
```

---

### Task 10: 値の意味の文書と vault の流れ図

**Files:**
- Modify: `USAGE.md`（差次的エクスポート 2 行と `arf2_annotate_identities` の行）、`docs/output_format/curation.md`
- 外部: `C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md`

**Interfaces:**
- Consumes: Task 6〜9 の保存形・戻り値・メタ行。

- [ ] **Step 1: 文書を直す**

- `USAGE.md`: 差次的エクスポート 2 行（`arf_export_differential` / `dataset_export_differential`）の `apply_curation` の説明に「assign は同定を置き換え（`name_source` / `inchikey_source` = `curation`）、redundant は除外する」を足す。`arf2_annotate_identities` の行の `curation_flag` の説明を「有効な判断（`wrong` / `suspect` / `assign:<記録名>` / `redundant`、無ければ空）」にする。
- `docs/output_format/curation.md`: 「### `curation_suggest` の戻り値」「### 候補付けの payload（HTML ビューア専用）」「### 候補の理由コード」（`polarity_mismatch` `dmz_out` がハード、`adduct_atypical` `trend_outlier` `no_matched_peaks` `msms_absent` `no_support` がソフト、`trend_unknown` `reference_unresolved` が情報）、「### イオン関係（`relation`）」（関係コードの一覧と、強い説明の定義、同位体の強度比の規則）、「### 判断の記録（assign / redundant）」の節を足し、メタ行の書式行を
  `# curation = <state>\tcuration_flags = N[\tcuration_wrong_excluded = N\tcuration_suspect = N[\tcuration_assigned = N\tcuration_redundant_excluded = N]]\tcuration_flags_sha256 = <hex>` に直す。必須注意事項に「候補の並び順はスペクトル類似度と制約だけで決まり、MS-DIAL の脂質規則（診断イオン・鎖決定）は評価していない」「④ の行の EIC の重ね描きで、相手 Y の代表試料のトレースは X の代表試料の最大値に合わせて縮めて描く（形の比較用。強度の比は `isotope_ratio` を見る）」を足す。

- [ ] **Step 2: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_readme_links.py tests/test_workflow_docs.py -q`
Expected: PASS（`docs/output_format/curation.md` は MCP リソースとして配信されるので、リソースのテスト `grep -rln "output-format" tests` も併せて回す）

- [ ] **Step 3: vault の流れ図を直す**

`C:\Users\yuu18\Documents\KnowledgeVault\30_Projects\ms-data-parser\ms-data-parser-flow.md` のキュレーションの枝に、`curation_review` の後の分岐として `curation_suggest`（前提: `library_load`・`curation_review`、どちらも `missing_state` で要求）→ ビューア → `curation_submit`（assign / redundant。`_tags.xml` は不変）→ エクスポートの置き換え（assign）・除外（redundant）、を足す。末尾の出典行の日付を今日に直す。

- [ ] **Step 4: コミット（バックグラウンド、出力をファイルへ）**

```bash
git add USAGE.md docs/output_format/curation.md
git commit -m "docs(curation): 候補付けの値の意味と判断の反映を文書にする"
```

出力の `FAILED` を `grep` で確かめる。

---

### Task 11: 実データでの検証

**Files:**
- Create: `scripts/verify_curation_suggest.py`
- 記録: main ツリーの `docs/HISTRY.md` / `docs/task.md`

**Interfaces:**
- Consumes: `suggest.run_suggestion`、`candidates.build_library_candidates`、`relations.find_relations`、`review.run_review`。

- [ ] **Step 1: 検証スクリプトを書く**

`scripts/verify_spectral_match.py` の体裁（argparse・標準出力に表）に合わせる。

```python
# scripts/verify_curation_suggest.py
"""候補付けの質を実データで測る（spec §11.1）。研究データのパスは引数で渡す（ここに書かない）。

1. 正解ありの検証: 元レビューの判定が ok・IsReferenceMatched・フラグ無しのスポットを正解とみなし、
   代表を隠して ② を回したとき、正解のレコードが 1 位・3 位以内に入る割合（分子種 / 和組成）。
2. ④ の再現: MS-DIAL の FoundInUpperMsMs リンク（注釈付き同士）のうち、Δm/z の関係表で名指しできた割合。
使い方: python scripts/verify_curation_suggest.py --arf2 <AlignmentResult_*.arf2> --library <*_Loaded.msp2.dbs>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lipidmix.arf2.ion_features import load_ion_features            # noqa: E402
from lipidmix.arf2.match_results import load_spot_annotations       # noqa: E402
from lipidmix.arf2.reader import load_catalog                       # noqa: E402
from lipidmix.curation import candidates, evidence, judge, relations, review, suggest, trend  # noqa: E402
from lipidmix.library.store import open_store                       # noqa: E402
from lipidmix.msdial.analysis_params import resolve_analysis_params # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arf2", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--limit", type=int, default=400)
    args = parser.parse_args(argv)
    th = judge.resolve_thresholds(None)
    store = open_store(args.library)
    catalog = load_catalog(args.arf2)
    annotations = load_spot_annotations(args.arf2)
    spots = evidence.select_spots(catalog, ontology=None, name_contains=None)
    base = review.run_review(args.arf2, spots, store=store, ms2_tol=0.025, th=th, file_ids=None,
                             max_traces=1, selection={})
    truth = [s for s in base["spots"] if s["verdict"] == "ok" and (s.get("match") or {}).get("is_reference_matched")
             and not s.get("flag")][: args.limit]
    scoring = suggest._scoring(store)
    evs, _ = evidence.collect(args.arf2, [next(c for c in catalog if c["MasterAlignmentID"] == t["spot_id"])
                                          for t in truth], store=store, ms2_tol=scoring["ms2_tol"], th=th,
                              max_traces=1, keep_measured=True)
    hits = {"species_top1": 0, "species_top3": 0, "sum_top1": 0, "sum_top3": 0, "n": 0}
    for ev in evs:
        rep = annotations[ev["spot_id"]]["representative"]
        out = candidates.build_library_candidates(ev, msdial_matches=[], representative=None, store=store,
                                                  scoring=scoring, trends={"classes": {}}, th=th,
                                                  exclude_current=False, top_n=3)
        names = [c["name"] for c in out["candidates"]]
        sums = [c["sum_name"] for c in out["candidates"]]
        truth_sum = trend.sum_composition(rep["name"])
        hits["n"] += 1
        hits["species_top1"] += names[:1] == [rep["name"]]
        hits["species_top3"] += rep["name"] in names
        hits["sum_top1"] += bool(truth_sum) and sums[:1] == [truth_sum]
        hits["sum_top3"] += bool(truth_sum) and truth_sum in sums
    n = max(hits["n"], 1)
    print("② 正解ありの検証（代表を隠した再検索）")
    for key in ("species_top1", "species_top3", "sum_top1", "sum_top3"):
        print(f"  {key}\t{hits[key]}/{hits['n']}\t{hits[key] / n:.1%}")

    features = load_ion_features(args.arf2)
    by_id = {s["MasterAlignmentID"]: s for s in catalog}
    params = resolve_analysis_params(args.arf2, next((s.get("IonMode") for s in catalog), None))
    named = total = 0
    for spot_id, feature in features.items():
        for link in feature["links"]:
            if link["kind"] != "found_in_upper_msms" or link["spot_id"] not in by_id:
                continue
            x, y = by_id[spot_id], by_id[link["spot_id"]]
            if (annotations.get(y["MasterAlignmentID"]) or {}).get("representative") is None:
                continue
            total += 1
            found = relations._explanations(
                {"mz": x["MassCenter"]}, {"mz": y["MassCenter"], "adduct": y.get("AdductType")},
                params["searched_adducts"], 0.010)
            named += bool(found)
    print(f"④ FoundInUpperMsMs（相手が注釈付き）のうち関係表で名指しできた割合\t{named}/{total}\t{named / max(total, 1):.1%}")
    store.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: kidney neg / pos で回す**

```bash
C:/Python314/python.exe scripts/verify_curation_suggest.py --arf2 "C:/Users/yuu18/datasets/a_lipidome_landscape_of_aging_in_mice/rplc/kidney/neg/AlignmentResult_2026_09_09_17_31_52.arf2" --library "C:/Users/yuu18/datasets/a_lipidome_landscape_of_aging_in_mice/rplc/kidney/neg/Dataset_2026_09_09_17_28_59_Loaded.msp2.dbs"
```

pos も同じフォルダ構成で回す（ファイル名は `ls` で確かめる）。数字を main ツリーの `docs/HISTRY.md` に日付見出しで追記する。**研究データの `curation/` には書かない**（スクリプトは保存しない）。

- [ ] **Step 3: MCP の実経路で候補付けを 1 回回し、ビューアを目視する**

worktree の `.mcp.json` のサーバで（または `check.py` から同じ関数を呼んで）kidney neg に対して `library_load(<Loaded.msp2.dbs>)` → `curation_review()`（既存のレビューがあればそれでよい）→ `curation_suggest()` を回す。これは研究データの `curation/` に `suggest-cs-*.json/.html` を書く（`curation/` はこの機能が書くフォルダで、ユーザーの判断記録 `flags.jsonl` には何も書かない）。書く前にユーザーに「kidney neg の curation フォルダに候補付けの結果を書いてよいか」を確かめる。
内蔵ブラウザで HTML を開き、wrong 31 件のカードで ④ の初期選択が付いたもの・① の候補・② の候補が見えること、未注釈の MS/MS ありのカードに候補が出ること、表示が重すぎないこと（JSON のサイズを `ls -l` で記録）を確かめる。`curation_submit` は**呼ばない**（ユーザーの判断なので）。

- [ ] **Step 4: 報告と記録**

`docs/task.md` に DONE / TODO を追記する。ユーザーへの報告には次を含める: ② の top1 / top3（分子種・和組成）、④ の再現率、kidney neg / pos の対象数と強い説明の件数、JSON の大きさ、未検証のもの（MCP Apps 表示、実際の採用率はユーザーの判断待ち）。
スクリプトのコミット:

```bash
git add scripts/verify_curation_suggest.py
git commit -m "chore(curation): 候補付けの質を実データで測る検証スクリプトを足す"
```

---

## Self-Review（書いた側の点検結果）

- spec §4 対象: Task 6 `select_targets`（flagged / likely_wrong / unannotated、include_decided）。§5 ①②と制約: Task 4。§5.4 和組成: Task 4 `sum_composition`、Task 6 `expand_entries`。§6 ④: Task 5（関係表・同位体比・相関・リンク）、Task 1（リンクの復号）、Task 2（アダクト・param）。§7 ツールとビューア: Task 7・8。§8 記録と下流: Task 3・9。§8.3 `_tags.xml` 不変: Task 8（`sync_misannotation` は wrong / clear だけを見る現行のまま。`test_suggest_then_submit_assign_and_redundant` が `tags_xml.added == []` と `removed == []` を確かめる）。§10 事実: Task 1 Step 1〜3。§11: Task 10・11。
- 腐敗防止テスト（ツール数・USAGE・workflow）が縛る文書は、ツールを足す Task 8 の中で直す（pre-commit の全テストを `--no-verify` なしで通すため）。値の意味の文書（output_format）と vault は Task 10。
- spec §7.2 の Y の EIC 重ね描き: Task 6（`relations[].partner_eic`、関係の上位 3 件だけ）と Task 7（`overlayFor`、`drawEic` の partner 描き分け）。ユーザー決定 2026-09-29 で縮小せず実装する。
- 型の一貫性: 候補 ID は library が `L<n>`、ion が `R<n>`。`preset` は ion の ID。`expand_entries` は assign に library の ID、redundant に ion の ID だけを受ける。`flags_for_arf2` の `redundant` は set、`assign` は dict。
