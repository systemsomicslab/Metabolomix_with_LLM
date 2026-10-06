# 群別強度プロット Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 選んだ脂質（代謝物）クラス・分子種ごとに、指定した試料群の試料別強度を 1 パネル 1 項目で描く MCP ツール
`arf_plot_group_intensity` と、その図を保存する `save_group_intensity_figure` を足す。

**Architecture:** MCP 非依存の純関数（項目の解決・群の解決・payload 組み立て・描画）を
`metabolomix/plots/group_intensity.py` に置き、ツール層（`metabolomix/arf/tools.py`）がセッションの ARF・兄弟 `.arf2`・
キュレーションの判断・最新レビューを集めて渡す。保存は既存の `save_*_figure` と同じ `metabolomix/tools/reports.py`。

**Tech Stack:** Python 3.14、FastMCP、matplotlib、pytest。

**Spec:** `docs/superpowers/specs/2026-10-06-group-intensity-plot-design.md`

## Global Constraints

- ARF 経路のみ。mzTab-M（`dataset_*`）は対象外。
- 描画データの形式名は `lipidmix.group_intensity.v1`（外部契約の接頭辞は `lipidmix.` のまま）。
- 全ツールに `structured_output=False`。JSON は `metabolomix.core.serialization.json_payload()` で返す。
- 前提状態が無いときは例外でなく `mcp_errors.missing_state(state, required_tools, message)` を返す。
- 縦軸は log10(PeakHeight) で全パネル共通。群の平均 ± SD は log10 空間、値 > 0 かつ低信頼でない試料だけ。
- 検定はしない。
- `items` は 1〜30 件、`groups` は 1 件以上。
- 標準液だけの分子種の判定倍率は 10（`STANDARD_ONLY_FOLD = 10.0`）。gap-fill 過半の判定は 0.5（`GAPFILL_MAJORITY = 0.5`）。
- 画像返しは dpi 100（`plots/render.figure_to_png` の既定）。保存は PNG dpi 300 と SVG。
- fixture はテスト自身が作る（`data/` `analyses/` の実ファイルに依存しない）。
- コミットは pre-commit で全テスト（約 4 分）が走る。自動化では背景実行にして出力をファイルへ落とす。

## Review Focus

- **同じ名前が複数スポットに当たる**（アダクト違い・重複スポット）→ 合計して 1 パネル。`n_spots` に数が出る。Task 2 のテストで固定。
- **群の指定が別の群の部分集合**（`"KO"` と `"KO_9w"`）→ 重なる試料は両方に描き、`caveats` に重なりを出す。Task 2 で固定。
- **全試料で値 0 の項目**（スポットはあるが強度が無い）→ N.D. パネルに「intensity 0 in all samples」と書き、「no annotated species」と区別。Task 4 で固定。
- **`assign` でクラスが変わったスポット**（`PC` と注釈されていたものを `PE` に付け替え）→ 付け替え後のクラスに数える。Task 2 で固定。
- **`arf_exclude` で除いた試料が群の指定に当たる** → 描かず、群の試料数にも数えない。全部除かれた群はエラー。Task 5 で固定。

---

## File Structure

| ファイル | 役割 |
|---|---|
| Create `metabolomix/plots/group_intensity.py` | 項目の解決（`resolve_items`）、標準液だけのスポット（`standard_only_spots`）、群の解決（`resolve_groups` / `resolve_sample_specs`）、payload（`build_group_intensity_payload`）、描画（`render_group_intensity_plot`） |
| Modify `metabolomix/msdial/analysis_params.py` | `Minimum peak height:` を `min_peak_height` として読む |
| Modify `metabolomix/core/session_state.py` | `ArfState.last_group_intensity` を足し、`reset_analysis` で消す |
| Modify `metabolomix/arf/tools.py` | ツール `arf_plot_group_intensity` |
| Modify `metabolomix/tools/reports.py` | ツール `save_group_intensity_figure` |
| Create `tests/test_group_intensity.py` | 純関数のテスト |
| Create `tests/test_group_intensity_tools.py` | ツール層のテスト |
| Modify テスト定数・文書 | 登録数・注釈・USAGE・CLAUDE の規模・workflow 文書・出力形式・`check.py` の片付け |

---

### Task 1: param ファイルから Minimum peak height を読む

**Files:**
- Modify: `metabolomix/msdial/analysis_params.py`
- Test: `tests/test_group_intensity.py`（新規）

**Interfaces:**
- Produces: `read_analysis_params(path) -> dict` に `"min_peak_height": float | None` が加わる。

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_group_intensity.py
"""群別強度プロット（metabolomix/plots/group_intensity.py ほか）の純関数テスト。"""
from __future__ import annotations

from metabolomix.msdial import analysis_params


def test_param_file_min_peak_height_is_read(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\nMinimum peak height: 1000\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] == 1000.0


def test_param_file_without_min_peak_height_gives_none(tmp_path):
    p = tmp_path / "Dataset_x_param_202610061200.txt"
    p.write_text("Ion mode: Negative\n", encoding="utf-8")
    assert analysis_params.read_analysis_params(p)["min_peak_height"] is None
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: FAIL（`KeyError: 'min_peak_height'`）

- [ ] **Step 3: 実装**

`metabolomix/msdial/analysis_params.py` の `_KEYS` に 1 行足し、`read_analysis_params` の初期値に足す。docstring の「読む行」にも `Minimum peak height:` を足す。

```python
_KEYS = {"Ion mode": "ion_mode", "Searched adduct ions": "searched_adducts",
         "MS1 tolerance for centroid": "ms1_tolerance",
         "Retention time tolerance for alignment": "rt_tolerance_alignment",
         "Minimum peak height": "min_peak_height"}
```

```python
    found = {"ion_mode": None, "searched_adducts": [], "rt_tolerance_alignment": None,
             "ms1_tolerance": None, "min_peak_height": None}
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py tests -k "analysis_params or group_intensity or suggest" -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/msdial/analysis_params.py tests/test_group_intensity.py
git commit -m "feat(msdial): param ファイルの Minimum peak height を読む"
```

---

### Task 2: 項目と群の解決（純関数）

**Files:**
- Create: `metabolomix/plots/group_intensity.py`
- Test: `tests/test_group_intensity.py`

**Interfaces:**
- Consumes: `metabolomix.msdial.sample_factors.build_sample_facets`, `expand_sample_specs`, `split_tokens`。
- Produces:
  - `INTERNAL_STANDARD: re.Pattern`（`\(d\d+\)`）、`STANDARD_ONLY_FOLD = 10.0`、`MAX_ITEMS = 30`
  - `resolve_items(items: list[str], catalog: list[dict], *, curation: dict | None = None, standard_only: frozenset[int] = frozenset()) -> tuple[list[dict], dict[str, list[str]]]`
    - `catalog` の要素は `{"MasterAlignmentID": int, "Name": str, "Ontology": str}`。
    - `curation` は `{"wrong": set[int], "redundant": set[int], "assign": {sid: {"name": str, "ontology": str}}, "auto_likely_wrong": set[int]}`（無いキーは空扱い）。
    - 戻り値の各項目: `{"item": str, "parts": [{"part": str, "kind": "class"|"name"|"none", "n_spots": int}], "spots": [{"spot_id": int, "name": str, "ontology": str}]}`（spots は spot_id 昇順・重複なし）。
    - 除外: `{"internal_standard": [...], "curation": [...], "auto_likely_wrong": [...], "standard_only": [...]}`（`"#<sid> <name>"`、重複なし・昇順）。
  - `standard_only_spots(rows_by_spot: dict[int, list[dict]], standard_samples: set[str], plotted_samples: set[str], fold: float = STANDARD_ONLY_FOLD) -> frozenset[int]`
  - `resolve_groups(group_specs: list[str], facets: dict) -> tuple[list[dict], list[str]]` → `([{"label": spec, "samples": [name, ...]}], caveats)`
  - `resolve_sample_specs(specs: list[str], facets: dict) -> set[str]`（低信頼の試料の解決。完全な試料名またはトークン、役割で絞らない）

- [ ] **Step 1: 失敗するテストを書く**（`tests/test_group_intensity.py` に追記）

```python
import pytest

from metabolomix.msdial import sample_factors
from metabolomix.plots import group_intensity as gi

CATALOG = [
    {"MasterAlignmentID": 0, "Name": "PG 34:1|PG 16:0_18:1", "Ontology": "PG"},
    {"MasterAlignmentID": 1, "Name": "low score: PG 35:1|PG 16:0_19:1", "Ontology": "PG"},
    {"MasterAlignmentID": 2, "Name": "PC 33:1(d7)|PC 15:0_18:1(d7)", "Ontology": "PC"},
    {"MasterAlignmentID": 3, "Name": "PC 34:1", "Ontology": "PC"},
    {"MasterAlignmentID": 4, "Name": "LPG 34:1", "Ontology": "LPG"},
    {"MasterAlignmentID": 5, "Name": "PG 34:1", "Ontology": "PG"},          # 同じ名前の別スポット
    {"MasterAlignmentID": 6, "Name": "PE O-34:1", "Ontology": "EtherPE"},
    {"MasterAlignmentID": 7, "Name": "Unknown", "Ontology": ""},
]


def _ids(item):
    return [s["spot_id"] for s in item["spots"]]


def test_class_item_matches_ontology_case_insensitively():
    items, _ = gi.resolve_items(["pg"], CATALOG)
    assert _ids(items[0]) == [0, 1, 5]
    assert items[0]["parts"] == [{"part": "pg", "kind": "class", "n_spots": 3}]


def test_name_item_matches_any_candidate_exactly_and_sums_duplicates():
    items, _ = gi.resolve_items(["PG 34:1", "pg 16:0_19:1"], CATALOG)
    assert _ids(items[0]) == [0, 5]               # LPG 34:1 には当たらない（部分一致しない）
    assert _ids(items[1]) == [1]                  # 接頭辞 low score: を除いた候補名
    assert items[0]["parts"][0]["kind"] == "name"


def test_plus_joins_parts_and_unknown_part_is_none():
    items, _ = gi.resolve_items(["PE+EtherPE", "SQDG"], CATALOG)
    assert _ids(items[0]) == [6]
    assert [p["kind"] for p in items[0]["parts"]] == ["none", "class"]
    assert items[1]["spots"] == [] and items[1]["parts"][0]["kind"] == "none"


def test_internal_standard_is_dropped_from_class_but_kept_by_name():
    items, excluded = gi.resolve_items(["PC", "PC 33:1(d7)"], CATALOG)
    assert _ids(items[0]) == [3]
    assert _ids(items[1]) == [2]
    assert excluded["internal_standard"] == ["#2 PC 33:1(d7)|PC 15:0_18:1(d7)"]


def test_curation_wrong_redundant_and_auto_likely_wrong_are_excluded():
    curation = {"wrong": {0}, "redundant": {5}, "assign": {}, "auto_likely_wrong": {1}}
    items, excluded = gi.resolve_items(["PG"], CATALOG, curation=curation)
    assert items[0]["spots"] == []
    assert excluded["curation"] == ["#0 PG 34:1|PG 16:0_18:1", "#5 PG 34:1"]
    assert excluded["auto_likely_wrong"] == ["#1 low score: PG 35:1|PG 16:0_19:1"]


def test_assign_moves_spot_to_the_recorded_class():
    curation = {"assign": {3: {"name": "PE 34:1", "ontology": "PE"}}}
    items, _ = gi.resolve_items(["PE", "PC"], CATALOG, curation=curation)
    assert items[0]["spots"] == [{"spot_id": 3, "name": "PE 34:1", "ontology": "PE"}]
    assert _ids(items[1]) == []                   # 内部標準 #2 は除かれ、#3 は PE へ移った


def test_standard_only_is_dropped_from_class_only():
    items, excluded = gi.resolve_items(["PG", "PG 34:1"], CATALOG, standard_only=frozenset({0}))
    assert _ids(items[0]) == [1, 5]
    assert _ids(items[1]) == [0, 5]
    assert excluded["standard_only"] == ["#0 PG 34:1|PG 16:0_18:1"]


def test_item_count_limits():
    with pytest.raises(ValueError):
        gi.resolve_items([], CATALOG)
    with pytest.raises(ValueError):
        gi.resolve_items(["PG"] * 31, CATALOG)
    with pytest.raises(ValueError):
        gi.resolve_items(["+"], CATALOG)


def test_standard_only_spots_compares_max_heights():
    rows = {0: [{"file_name": "std_1", "height": 5000.0}, {"file_name": "ctrl_1", "height": 100.0}],
            1: [{"file_name": "std_1", "height": 5000.0}, {"file_name": "ctrl_1", "height": 900.0}]}
    assert gi.standard_only_spots(rows, {"std_1"}, {"ctrl_1"}) == frozenset({0})


def _facets(names):
    return sample_factors.build_sample_facets(names)


def test_groups_keep_order_and_skip_blank_unless_named():
    facets = _facets(["x_blank_1", "x_ctrl_1", "x_ctrl_2", "x_ko_1", "x_ko_2"])
    groups, caveats = gi.resolve_groups(["ko", "ctrl"], facets)
    assert [g["label"] for g in groups] == ["ko", "ctrl"]
    assert groups[1]["samples"] == ["x_ctrl_1", "x_ctrl_2"]
    groups, _ = gi.resolve_groups(["blank", "ctrl"], facets)
    assert groups[0]["samples"] == ["x_blank_1"]
    assert caveats == []


def test_overlapping_groups_are_drawn_in_both_with_a_caveat():
    facets = _facets(["x_ko_9w_1", "x_ko_24m_1", "x_ctrl_9w_1"])
    groups, caveats = gi.resolve_groups(["ko", "ko_9w"], facets)
    assert groups[1]["samples"] == ["x_ko_9w_1"]
    assert any("x_ko_9w_1" in c for c in caveats)


def test_unmatched_or_empty_group_raises():
    facets = _facets(["x_ctrl_1"])
    with pytest.raises(ValueError):
        gi.resolve_groups(["ko"], facets)
    with pytest.raises(ValueError):
        gi.resolve_groups([], facets)
    with pytest.raises(ValueError):
        gi.resolve_groups(["_"], facets)


def test_low_reliability_specs_accept_names_and_tokens():
    facets = _facets(["x_ctrl_1", "x_ko_1", "x_ko_2"])
    assert gi.resolve_sample_specs(["x_ko_1"], facets) == {"x_ko_1"}
    assert gi.resolve_sample_specs(["ko"], facets) == {"x_ko_1", "x_ko_2"}
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: FAIL（`ModuleNotFoundError: metabolomix.plots.group_intensity`）

- [ ] **Step 3: 実装**（`metabolomix/plots/group_intensity.py` を新規作成）

```python
"""群別強度プロット: 選んだクラス・分子種ごとに、試料群の試料別強度を並べる（MCP 非依存）。

1 パネル = 1 項目、1 点 = 1 試料（項目に当たったスポットの PeakHeight の合計）、群ごとに
log10 の平均 ± SD。検定はしない（`arf_differential` の役割）。spec:
docs/superpowers/specs/2026-10-06-group-intensity-plot-design.md。
"""
from __future__ import annotations

import math
import re
import statistics

from metabolomix.msdial import sample_factors

GROUP_INTENSITY_SCHEMA = "lipidmix.group_intensity.v1"
MAX_ITEMS = 30
#: 重水素などで標識した内部標準（`PC 33:1(d7)`、`FA 16:0(d3)`）。クラスの合計に入れない。
INTERNAL_STANDARD = re.compile(r"\(d\d+\)")
#: 標準液での最大高さが、描く試料での最大高さのこの倍数以上なら「標準液にだけある」。
STANDARD_ONLY_FOLD = 10.0
#: 合計のうち gap-fill の値がこの割合を超える点は形を変える。
GAPFILL_MAJORITY = 0.5
_PREFIX = re.compile(r"^\s*(?:low score|no MS2|w/o MS2|unsettled)\s*:\s*", re.IGNORECASE)
_EXCLUDED_KEYS = ("internal_standard", "curation", "auto_likely_wrong", "standard_only")


def _candidate_names(name: str) -> set[str]:
    return {_PREFIX.sub("", part).strip().casefold() for part in (name or "").split("|") if part.strip()}


def resolve_items(items, catalog, *, curation=None, standard_only=frozenset()):
    """項目（クラス・名前、`+` で合算）をスポット集合に解決する。spec §3.1 / §5。"""
    items = [str(item) for item in (items or [])]
    if not 1 <= len(items) <= MAX_ITEMS:
        raise ValueError(f"items は 1〜{MAX_ITEMS} 件で指定してください（受け取った数: {len(items)}）。")
    curation = curation or {}
    wrong = set(curation.get("wrong") or ()) | set(curation.get("redundant") or ())
    auto = set(curation.get("auto_likely_wrong") or ())
    assign = curation.get("assign") or {}

    spots = []
    for row in catalog:
        sid = int(row["MasterAlignmentID"])
        name, ontology = row.get("Name") or "", (row.get("Ontology") or "").strip()
        if sid in assign:
            name = assign[sid].get("name") or name
            ontology = (assign[sid].get("ontology") or ontology).strip()
        spots.append({"spot_id": sid, "name": name, "ontology": ontology})
    ontologies = {s["ontology"].casefold() for s in spots if s["ontology"]}

    excluded = {key: set() for key in _EXCLUDED_KEYS}

    def keep(spot, *, as_class):
        label = f"#{spot['spot_id']} {spot['name']}"
        if spot["spot_id"] in wrong:
            excluded["curation"].add(label)
            return False
        if spot["spot_id"] in auto:
            excluded["auto_likely_wrong"].add(label)
            return False
        if as_class and INTERNAL_STANDARD.search(spot["name"]):
            excluded["internal_standard"].add(label)
            return False
        if as_class and spot["spot_id"] in standard_only:
            excluded["standard_only"].add(label)
            return False
        return True

    resolved = []
    for item in items:
        parts = [p.strip() for p in item.split("+") if p.strip()]
        if not parts:
            raise ValueError(f"空の項目は指定できません: {item!r}")
        chosen: dict[int, dict] = {}
        out_parts = []
        for part in parts:
            key = part.casefold()
            if key in ontologies:
                hits = [s for s in spots if s["ontology"].casefold() == key and keep(s, as_class=True)]
                kind = "class"
            else:
                hits = [s for s in spots if key in _candidate_names(s["name"]) and keep(s, as_class=False)]
                kind = "name" if hits else "none"
            for s in hits:
                chosen[s["spot_id"]] = s
            out_parts.append({"part": part, "kind": kind, "n_spots": len(hits)})
        resolved.append({"item": item, "parts": out_parts,
                         "spots": [chosen[k] for k in sorted(chosen)]})
    return resolved, {key: sorted(values) for key, values in excluded.items()}


def standard_only_spots(rows_by_spot, standard_samples, plotted_samples, fold=STANDARD_ONLY_FOLD):
    """標準液の試料にだけ実質的な強度があるスポット。spec §5-4。"""
    found = set()
    for sid, rows in rows_by_spot.items():
        std = max((float(r.get("height") or 0.0) for r in rows if r.get("file_name") in standard_samples),
                  default=0.0)
        grp = max((float(r.get("height") or 0.0) for r in rows if r.get("file_name") in plotted_samples),
                  default=0.0)
        if std > 0 and std >= fold * max(grp, 1.0):
            found.add(sid)
    return frozenset(found)


def resolve_groups(group_specs, facets):
    """群の指定を試料名に解決する（指定順）。spec §3.2。"""
    specs = [str(s) for s in (group_specs or [])]
    if not specs:
        raise ValueError("groups を 1 群以上指定してください。")
    groups = []
    for spec in specs:
        tokens = sample_factors.split_tokens(spec)
        if not tokens:
            raise ValueError(f"空の群指定は使えません: {spec!r}")
        roles = ["sample"] + [role for role in ("blank", "qc") if role in tokens]
        matches, _ = sample_factors.expand_sample_specs([spec], facets, include_roles=tuple(roles))
        groups.append({"label": spec, "samples": list(matches[spec])})
    caveats = []
    seen: dict[str, str] = {}
    for group in groups:
        for name in group["samples"]:
            if name in seen and seen[name] != group["label"]:
                caveats.append(f"試料 {name} は群 {seen[name]} と {group['label']} の両方に当たり、両方に描きました。")
            seen.setdefault(name, group["label"])
    return groups, caveats


def resolve_sample_specs(specs, facets):
    """試料名（完全一致）またはトークン指定を試料名の集合にする（役割で絞らない）。"""
    found = set()
    for spec in specs or []:
        if spec in facets:
            found.add(spec)
            continue
        matches, _ = sample_factors.expand_sample_specs([spec], facets, include_roles=None)
        found.update(matches.get(spec, []))
    return found
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/plots/group_intensity.py tests/test_group_intensity.py
git commit -m "feat(plots): 群別強度プロットの項目と群の解決"
```

---

### Task 3: payload の組み立て

**Files:**
- Modify: `metabolomix/plots/group_intensity.py`
- Test: `tests/test_group_intensity.py`

**Interfaces:**
- Consumes: Task 2 の `resolve_items` / `resolve_groups` の戻り値。
- Produces: `build_group_intensity_payload(resolved_items, groups, rows_by_spot, *, msms=None, low_reliability=frozenset(), excluded=None, detection_limit=None, detection_limit_source=None, caveats=None) -> dict`（spec §7 の形。`rows_by_spot` の行は `file_name` `height` `is_gap_filled` を持つ）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def _rows():
    return {
        0: [{"file_name": "c1", "height": 100.0, "is_gap_filled": False},
            {"file_name": "c2", "height": 1000.0, "is_gap_filled": False},
            {"file_name": "k1", "height": 10.0, "is_gap_filled": True}],
        1: [{"file_name": "c1", "height": 900.0, "is_gap_filled": True},
            {"file_name": "c2", "height": 0.0, "is_gap_filled": True},
            {"file_name": "k1", "height": 0.0, "is_gap_filled": True}],
    }


def _payload(**kw):
    items = [{"item": "PG", "parts": [{"part": "PG", "kind": "class", "n_spots": 2}],
              "spots": [{"spot_id": 0, "name": "PG 34:1", "ontology": "PG"},
                        {"spot_id": 1, "name": "no MS2: PG 35:1", "ontology": "PG"}]},
             {"item": "PE", "parts": [{"part": "PE", "kind": "none", "n_spots": 0}], "spots": []}]
    groups = [{"label": "ctrl", "samples": ["c1", "c2"]}, {"label": "ko", "samples": ["k1"]}]
    return gi.build_group_intensity_payload(items, groups, _rows(), msms={0: True, 1: False}, **kw)


def test_values_are_sums_with_gap_fill_fraction():
    p = _payload()
    ctrl = p["items"][0]["groups"][0]["samples"]
    assert ctrl[0] == {"sample": "c1", "value": 1000.0, "gap_filled_fraction": 0.9, "low_reliability": False}
    assert ctrl[1]["value"] == 1000.0 and ctrl[1]["gap_filled_fraction"] == 0.0
    assert p["plot_schema"] == "lipidmix.group_intensity.v1"


def test_log10_mean_sd_and_low_reliability_exclusion():
    p = _payload()
    g = p["items"][0]["groups"][0]
    assert g["log10_mean"] == 3.0 and g["log10_sd"] == 0.0 and g["n_in_stats"] == 2
    p = _payload(low_reliability=frozenset({"c2"}))
    g = p["items"][0]["groups"][0]
    assert g["n_in_stats"] == 1 and g["log10_sd"] is None
    assert g["samples"][1]["low_reliability"] is True


def test_zero_values_are_kept_and_not_in_stats():
    p = _payload()
    k = p["items"][0]["groups"][1]
    assert k["samples"][0]["value"] == 10.0
    assert k["n_in_stats"] == 1
    pe = p["items"][1]
    assert pe["detected"] is False and pe["n_spots"] == 0
    assert pe["groups"][0]["samples"][0]["value"] == 0.0 and pe["groups"][0]["log10_mean"] is None


def test_msms_count_and_metadata_are_carried():
    p = _payload(detection_limit=1000.0, detection_limit_source="argument",
                 excluded={"internal_standard": ["#2 x"]}, caveats=["c"])
    assert p["items"][0]["n_spots_msms"] == 1
    assert p["detection_limit"] == 1000.0 and p["detection_limit_source"] == "argument"
    assert p["excluded"]["internal_standard"] == ["#2 x"] and p["caveats"] == ["c"]
    assert p["groups"] == [{"label": "ctrl", "samples": ["c1", "c2"]}, {"label": "ko", "samples": ["k1"]}]
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: FAIL（`AttributeError: ... build_group_intensity_payload`）

- [ ] **Step 3: 実装**（`group_intensity.py` に追記）

```python
def build_group_intensity_payload(resolved_items, groups, rows_by_spot, *, msms=None,
                                  low_reliability=frozenset(), excluded=None,
                                  detection_limit=None, detection_limit_source=None, caveats=None):
    """項目 × 群 × 試料のクラス合計強度を組み立てる（描画しない）。spec §4 / §7。"""
    index = {sid: {r.get("file_name"): r for r in rows} for sid, rows in rows_by_spot.items()}
    out_items = []
    for item in resolved_items:
        sids = [s["spot_id"] for s in item["spots"]]
        out_groups = []
        for group in groups:
            samples = []
            for name in group["samples"]:
                total = gap = 0.0
                for sid in sids:
                    row = index.get(sid, {}).get(name)
                    if row is None:
                        continue
                    h = float(row.get("height") or 0.0)
                    total += h
                    if row.get("is_gap_filled"):
                        gap += h
                samples.append({"sample": name, "value": round(total, 1),
                                "gap_filled_fraction": round(gap / total, 3) if total > 0 else None,
                                "low_reliability": name in low_reliability})
            logs = [math.log10(s["value"]) for s in samples if s["value"] > 0 and not s["low_reliability"]]
            out_groups.append({
                "label": group["label"], "samples": samples,
                "log10_mean": round(statistics.mean(logs), 4) if logs else None,
                "log10_sd": round(statistics.stdev(logs), 4) if len(logs) > 1 else None,
                "n_in_stats": len(logs),
            })
        out_items.append({
            "item": item["item"], "parts": item["parts"], "spots": item["spots"],
            "n_spots": len(sids),
            "n_spots_msms": sum(1 for sid in sids if (msms or {}).get(sid)) if msms is not None else None,
            "detected": any(s["value"] > 0 for g in out_groups for s in g["samples"]),
            "groups": out_groups,
        })
    return {
        "plot_schema": GROUP_INTENSITY_SCHEMA,
        "value": "sum of PeakHeight (log10 on the plot)",
        "detection_limit": detection_limit, "detection_limit_source": detection_limit_source,
        "groups": [{"label": g["label"], "samples": list(g["samples"])} for g in groups],
        "low_reliability_samples": sorted(low_reliability),
        "items": out_items,
        "excluded": {key: list((excluded or {}).get(key, [])) for key in _EXCLUDED_KEYS},
        "caveats": list(caveats or []),
    }
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/plots/group_intensity.py tests/test_group_intensity.py
git commit -m "feat(plots): 群別強度プロットの payload を組み立てる"
```

---

### Task 4: 描画

**Files:**
- Modify: `metabolomix/plots/group_intensity.py`
- Test: `tests/test_group_intensity.py`

**Interfaces:**
- Consumes: Task 3 の payload。
- Produces: `render_group_intensity_plot(payload: dict, *, title: str | None = None, ncols: int | None = None)` → matplotlib Figure（呼び出し側が閉じる。`figure_to_png` は閉じる）。

- [ ] **Step 1: 失敗するテストを書く**

```python
def _texts(ax):
    return [t.get_text() for t in ax.texts]


def test_render_shares_the_y_axis_and_marks_nd_and_ms1_only():
    import matplotlib.pyplot as plt
    p = _payload(detection_limit=1000.0)
    fig = gi.render_group_intensity_plot(p)
    try:
        axes = [ax for ax in fig.axes if ax.get_visible() and ax.axison]
        assert len(axes) == 2
        assert axes[0].get_ylim() == axes[1].get_ylim()
        assert any("no annotated species" in t for t in _texts(axes[1]))
        assert any(l.get_linestyle() == "--" for l in axes[0].get_lines())   # 検出下限の破線
    finally:
        plt.close(fig)


def test_render_distinguishes_all_zero_from_no_species():
    import matplotlib.pyplot as plt
    rows = {9: [{"file_name": "c1", "height": 0.0, "is_gap_filled": True}]}
    items = [{"item": "PS", "parts": [{"part": "PS", "kind": "class", "n_spots": 1}],
              "spots": [{"spot_id": 9, "name": "PS 34:1", "ontology": "PS"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "ctrl", "samples": ["c1"]}], rows,
                                         msms={9: False})
    fig = gi.render_group_intensity_plot(p)
    try:
        assert any("intensity 0 in all samples" in t for t in _texts(fig.axes[0]))
    finally:
        plt.close(fig)


def test_render_marks_ms1_only_panels():
    import matplotlib.pyplot as plt
    rows = {1: [{"file_name": "c1", "height": 500.0, "is_gap_filled": False}]}
    items = [{"item": "PG 35:1", "parts": [{"part": "PG 35:1", "kind": "name", "n_spots": 1}],
              "spots": [{"spot_id": 1, "name": "no MS2: PG 35:1", "ontology": "PG"}]}]
    p = gi.build_group_intensity_payload(items, [{"label": "ctrl", "samples": ["c1"]}], rows,
                                         msms={1: False})
    fig = gi.render_group_intensity_plot(p)
    try:
        assert any("MS1-only" in t for t in _texts(fig.axes[0]))
    finally:
        plt.close(fig)
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: FAIL（`AttributeError: ... render_group_intensity_plot`）

- [ ] **Step 3: 実装**（`group_intensity.py` に追記。`check.py` の試作 `render_class_intensity_plot` の移植。matplotlib は関数内 import）

```python
GROUP_COLORS = ["#bdbdbd", "#9ecae1", "#3182bd", "#de2d26", "#fd8d3c", "#31a354", "#756bb1", "#8c564b"]


def render_group_intensity_plot(payload, *, title=None, ncols=None):
    """payload を 1 項目 1 パネルの点 + 平均 ± SD 図にする。spec §4。

    縦軸は全パネル共通（log10）。パネルごとに拡大すると強度 10 前後のノイズが信号に見え、
    「見つかるか」の判断を誤らせるため。0 の試料は軸の底に ▽。
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    items = payload["items"]
    labels = [g["label"] for g in payload["groups"]]
    ncols = ncols or min(5, len(items))
    nrows = math.ceil(len(items) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 2.7 * nrows), squeeze=False)
    positive = [s["value"] for it in items for g in it["groups"] for s in g["samples"] if s["value"] > 0]
    floor = 0.0
    top = math.ceil(math.log10(max(positive)) + 0.2) if positive else 1.0
    limit = payload.get("detection_limit")
    for ax, item in zip(axes.flat, items):
        ms1_only = item["n_spots"] > 0 and item.get("n_spots_msms") == 0
        msms_note = "" if item.get("n_spots_msms") is None else f", MS/MS {item['n_spots_msms']}"
        ax.set_title(f"{item['item']}  (spots {item['n_spots']}{msms_note})", fontsize=9,
                     color="#777777" if ms1_only else "black", style="italic" if ms1_only else "normal")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
        ax.set_xlim(-0.6, len(labels) - 0.4)
        ax.set_ylim(floor - 0.3, top)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylabel("log10(Intensity)", fontsize=7)
        if limit:
            ax.axhline(math.log10(limit), color="#888888", linestyle="--", linewidth=0.8, zorder=1)
        if ms1_only:
            ax.set_facecolor("#f4f4f4")
            ax.text(0.03, 0.95, "MS1-only (unconfirmed)", transform=ax.transAxes, ha="left", va="top",
                    fontsize=7, color="#777777", style="italic")
        if not item["detected"]:
            reason = "no annotated species" if item["n_spots"] == 0 else "intensity 0 in all samples"
            ax.text(0.5, 0.5, f"N.D.\n({reason})", transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="#555555",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 2}, zorder=5)
            continue
        for i, group in enumerate(item["groups"]):
            color = GROUP_COLORS[i % len(GROUP_COLORS)]
            n = len(group["samples"])
            for j, s in enumerate(group["samples"]):
                x = i + (j - (n - 1) / 2) * 0.12
                if s["value"] <= 0:
                    ax.scatter(x, floor, marker="v", s=22, facecolor="white", edgecolor="#555555", zorder=3)
                    continue
                gapfilled = (s["gap_filled_fraction"] or 0) > GAPFILL_MAJORITY
                ax.scatter(x, math.log10(s["value"]), s=30, marker="D" if gapfilled else "o",
                           facecolor="white" if s["low_reliability"] else color,
                           edgecolor="black", linewidth=0.6, zorder=3)
            if group["log10_mean"] is not None:
                m, sd = group["log10_mean"], group["log10_sd"]
                ax.hlines(m, i - 0.28, i + 0.28, color="black", linewidth=1.6, zorder=4)
                if sd is not None:
                    ax.errorbar(i, m, yerr=sd, color="black", capsize=4, linewidth=1.0, zorder=4)
        ax.axhline(floor, color="#dddddd", linewidth=0.6, zorder=1)
    for ax in list(axes.flat)[len(items):]:
        ax.axis("off")
    legend = [
        Line2D([], [], marker="o", ls="", mfc="#888888", mec="black", label="sample (mostly detected peaks)"),
        Line2D([], [], marker="D", ls="", mfc="#888888", mec="black", label="sample (mostly gap-filled values)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="black", label="low-reliability sample (not in mean/SD)"),
        Line2D([], [], marker="v", ls="", mfc="white", mec="#555555", label="zero (at the axis floor)"),
        Line2D([], [], color="black", lw=1.6, label="mean ± SD of log10"),
    ]
    if limit:
        legend.append(Line2D([], [], color="#888888", ls="--", lw=0.8,
                             label=f"detection limit ({limit:,.0f})"))
    fig.legend(handles=legend, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.suptitle(title or "Intensity by sample group", fontsize=11)
    fig.tight_layout(rect=(0, 0.08, 1, 0.96))
    return fig
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/plots/group_intensity.py tests/test_group_intensity.py
git commit -m "feat(plots): 群別強度プロットを描く"
```

---

### Task 5: ツール `arf_plot_group_intensity`

**Files:**
- Modify: `metabolomix/core/session_state.py`（`ArfState.__init__` と `reset_analysis` に `self.last_group_intensity = None`）
- Modify: `metabolomix/arf/tools.py`（ツール本体、`__all__` に追加）
- Test: `tests/test_group_intensity_tools.py`（新規）

**Interfaces:**
- Consumes: Task 1〜4 の関数、`metabolomix.arf.tools._sibling_arf2_path()`、`metabolomix.arf2.reader.load_catalog`、
  `metabolomix.arf2.match_results.load_spot_annotations`、`metabolomix.curation.apply.flags_for_arf2`、
  `metabolomix.curation.flags.alignment_key` / `FlagFileError`、`metabolomix.curation.suggest.latest_review`、
  `metabolomix.arf.reader.alignment_feature_row`、`metabolomix.msdial.sample_factors.arf_sample_names` / `build_sample_facets`、
  `metabolomix.msdial.analysis_params.find_param_file` / `read_analysis_params`、`metabolomix.plots.render`。
- Produces: `arf_plot_group_intensity(items: list[str], groups: list[str], low_reliability_samples: list[str] | None = None, apply_curation: bool = True, exclude_auto_likely_wrong: bool = False, standard_samples: list[str] | None = None, detection_limit: float | None = None, title: str | None = None, ncols: int | None = None, output: str | None = None) -> list | str`。
  `session.arf.last_group_intensity` に payload（と `title`）を保持: `{"payload": dict, "title": str | None, "ncols": int | None}`。

- [ ] **Step 1: 失敗するテストを書く**

```python
# tests/test_group_intensity_tools.py
"""arf_plot_group_intensity / save_group_intensity_figure（ツール層）。fixture はテストが作る。"""
from __future__ import annotations

import json

import pytest

from metabolomix.arf import reader as arf_reader
from metabolomix.core import session_state
from metabolomix.curation import flags as curation_flags
from tests.curation_fixtures import arf2_spot_raw, arf_row, match_result, write_arf, write_arf2

STEM = "AlignmentResult_2026_01_01_00_00_00"
SAMPLES = ["blank_1", "ctrl_1", "ctrl_2", "ko_1", "ko_2", "std_1"]


def _group(heights, gap=()):
    rows = []
    for fid, (name, h) in enumerate(zip(SAMPLES, heights)):
        row = arf_row(file_id=fid, mz=700.0, rt=3.0, height=h, gap_filled=name in gap)
        row[1] = f"20261006_{name}"
        rows.append(row)
    return rows


@pytest.fixture()
def loaded(tmp_path):
    arf = write_arf(tmp_path / f"{STEM}_PeakProperties.arf", [
        _group([10, 5000, 6000, 2000, 2500, 10]),            # 0 PG 34:1
        _group([10, 3000, 3500, 4000, 4200, 10]),            # 1 PG 35:1
        _group([1500, 1400, 1300, 1450, 1500, 10]),          # 2 PC 33:1(d7) 内部標準
        _group([5, 800, 900, 700, 600, 5]),                  # 3 no MS2: PC 34:1
        _group([0, 20, 30, 25, 20, 9000]),                   # 4 PG 31:1 標準液だけ
    ])
    arf2 = write_arf2(tmp_path / f"{STEM}.arf2", [
        arf2_spot_raw(spot_id=0, name="PG 34:1|PG 16:0_18:1", ontology="PG"),
        arf2_spot_raw(spot_id=1, name="PG 35:1|PG 16:0_19:1", ontology="PG"),
        arf2_spot_raw(spot_id=2, name="PC 33:1(d7)|PC 15:0_18:1(d7)", ontology="PC"),
        arf2_spot_raw(spot_id=3, name="no MS2: PC 34:1", ontology="PC",
                      matches=[match_result({3: -1.0})]),
        arf2_spot_raw(spot_id=4, name="PG 31:1|PG 17:0_14:1", ontology="PG"),
    ])
    session_state.session = session_state.AnalysisSession()
    with open(arf, "rb") as handle:
        session_state.session.arf.features = arf_reader.deserialize(handle)
    session_state.session.arf.current_file_path = str(arf)
    session_state.session.arf.class_index = None
    return {"arf": arf, "arf2": arf2}


def _payload(**kw):
    from metabolomix.arf.tools import arf_plot_group_intensity
    return json.loads(arf_plot_group_intensity(output="payload", **kw))


def test_missing_arf_returns_missing_state():
    from metabolomix.arf.tools import arf_plot_group_intensity
    session_state.session = session_state.AnalysisSession()
    out = json.loads(arf_plot_group_intensity(items=["PG"], groups=["ctrl"], output="payload"))
    assert out["error"]["code"] == "missing_state"


def test_payload_classes_names_and_nd(loaded):
    p = _payload(items=["PG", "PC", "PE", "PG 34:1"], groups=["ctrl", "ko"])
    by = {it["item"]: it for it in p["items"]}
    assert [s["spot_id"] for s in by["PG"]["spots"]] == [0, 1, 4]
    assert [s["spot_id"] for s in by["PC"]["spots"]] == [3]          # 内部標準は除く
    assert by["PC"]["n_spots_msms"] == 0
    assert by["PE"]["detected"] is False
    assert [s["spot_id"] for s in by["PG 34:1"]["spots"]] == [0]
    assert [g["label"] for g in p["groups"]] == ["ctrl", "ko"]
    assert p["groups"][0]["samples"] == ["20261006_ctrl_1", "20261006_ctrl_2"]


def test_standard_samples_drop_standard_only_spots(loaded):
    p = _payload(items=["PG"], groups=["ctrl", "ko"], standard_samples=["std"])
    assert [s["spot_id"] for s in p["items"][0]["spots"]] == [0, 1]
    assert p["excluded"]["standard_only"] == ["#4 PG 31:1|PG 17:0_14:1"]


def test_curation_flags_are_applied_by_default(loaded):
    store = curation_flags.FlagStore(curation_flags.curation_dir(loaded["arf2"]))
    key = curation_flags.alignment_key(loaded["arf2"])
    store.append([{"spot_id": 1, "flag": "wrong", "note": ""}], alignment=key,
                 review_id=None, source="test")
    p = _payload(items=["PG"], groups=["ctrl"])
    assert 1 not in [s["spot_id"] for s in p["items"][0]["spots"]]
    p = _payload(items=["PG"], groups=["ctrl"], apply_curation=False)
    assert 1 in [s["spot_id"] for s in p["items"][0]["spots"]]


def test_auto_likely_wrong_without_review_adds_caveat(loaded):
    p = _payload(items=["PG"], groups=["ctrl"], exclude_auto_likely_wrong=True)
    assert any("レビュー" in c for c in p["caveats"])


def test_excluded_samples_and_low_reliability(loaded):
    session_state.session.arf.excluded_samples = {"20261006_ko_2"}
    p = _payload(items=["PG"], groups=["ctrl", "ko"], low_reliability_samples=["ctrl_2"])
    assert p["groups"][1]["samples"] == ["20261006_ko_1"]
    ctrl = p["items"][0]["groups"][0]
    assert ctrl["n_in_stats"] == 1 and ctrl["samples"][1]["low_reliability"] is True


def test_group_fully_excluded_raises(loaded):
    session_state.session.arf.excluded_samples = {"20261006_ko_1", "20261006_ko_2"}
    with pytest.raises(ValueError):
        _payload(items=["PG"], groups=["ctrl", "ko"])


def test_detection_limit_from_param_file_and_argument(loaded, tmp_path):
    (tmp_path / "Dataset_x_param_202610061200.txt").write_text(
        "Minimum peak height: 1000\n", encoding="utf-8")
    p = _payload(items=["PG"], groups=["ctrl"])
    assert p["detection_limit"] == 1000.0 and p["detection_limit_source"] == "param_file"
    p = _payload(items=["PG"], groups=["ctrl"], detection_limit=500)
    assert p["detection_limit"] == 500.0 and p["detection_limit_source"] == "argument"


def test_image_output_and_session_keeps_payload(loaded):
    from mcp.server.fastmcp import Image
    from metabolomix.arf.tools import arf_plot_group_intensity
    out = arf_plot_group_intensity(items=["PG", "PE"], groups=["ctrl", "ko"])
    assert isinstance(out, list) and isinstance(out[1], Image)
    assert "2 項目" in out[0] and "N.D.: PE" in out[0]
    assert session_state.session.arf.last_group_intensity["payload"]["items"][0]["item"] == "PG"
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity_tools.py -q`
Expected: FAIL（`ImportError: cannot import name 'arf_plot_group_intensity'`）

- [ ] **Step 3: 実装**

`metabolomix/core/session_state.py` の `ArfState.__init__`（`last_differential` の次）と `reset_analysis`（`last_differential = None` の次）に:

```python
        self.last_group_intensity = None  # arf_plot_group_intensity / save_group_intensity_figure が参照
```

`metabolomix/arf/tools.py` の `__all__` に `"arf_plot_group_intensity"` を足し、`arf_plot_volcano` の後ろに:

```python
@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True), structured_output=False)
def arf_plot_group_intensity(
    items: list[str],
    groups: list[str],
    low_reliability_samples: list[str] | None = None,
    apply_curation: bool = True,
    exclude_auto_likely_wrong: bool = False,
    standard_samples: list[str] | None = None,
    detection_limit: float | None = None,
    title: str | None = None,
    ncols: int | None = None,
    output: str | None = None,
) -> list | str:
    """選んだクラス・分子種ごとに、試料群の試料別強度を並べる（1 項目 1 パネル。既定は PNG 画像）。

    - items: 1 項目 1 パネル（1〜30 件）。`.arf2` の Ontology に完全一致すればクラス（全分子種の合計）、
      しなければ分子種名（`|` で分けた候補名に完全一致。`low score:` 等の接頭辞は無視）。
      `"PE+EtherPE"` のように `+` で合算。当たらない項目は N.D. のパネルになる。
    - groups: 1 群以上。`arf_differential` と同じトークン規則（`"KO_9w"` は両方を含む試料）。
      QC・ブランクは群の指定に `qc` / `blank` を書いたときだけ入る。
    - 1 点 = 1 試料の PeakHeight 合計（gap-fill 含む）。縦軸は log10 で全パネル共通、群は log10 の平均 ± SD。
      検定はしない（`arf_differential` を使う）。
    - クラスの合計からは標識内部標準（`(d7)` 等）と、standard_samples の試料にだけある分子種を除く。
      apply_curation（既定 True）で人の判断（wrong / redundant 除外、assign 付け替え）、
      exclude_auto_likely_wrong で最新レビューの likely_wrong も除く。除いたものは payload の excluded。
    - low_reliability_samples: 白抜きで描き、平均 ± SD から外す。
    - detection_limit: 検出下限の破線。省略時は param ファイルの Minimum peak height。
    - output: "image"（既定）/ "payload"（`lipidmix.group_intensity.v1`）。PNG ファイルが要るときは
      save_group_intensity_figure。
    """
    from metabolomix.arf2.match_results import load_spot_annotations, name_prefix
    from metabolomix.curation import apply as curation_apply
    from metabolomix.curation import flags as curation_flags
    from metabolomix.curation import suggest as curation_suggest
    from metabolomix.msdial import analysis_params
    from metabolomix.plots import group_intensity as gi

    mode = plot_render.resolve_plot_output(output)
    arf_state = session_state.session.arf
    if arf_state.features is None or not str(arf_state.current_file_path or "").lower().endswith(".arf"):
        return mcp_errors.missing_state(
            "arf_dataset", ["arf_parser", "load_dataset"], "先に load_dataset で ARF データを読み込んでください。")
    arf2_path = _sibling_arf2_path()
    if arf2_path is None:
        return mcp_errors.missing_state(
            "sibling_arf2", ["arf_parser", "load_dataset"],
            "同じアラインメントの .arf2 が見つかりません（名前とクラスの解決に要ります）。")

    excluded_samples = set(arf_state.excluded_samples or ())
    rows_by_spot = {}
    for feature in arf_state.features:
        rows = [arf_reader.alignment_feature_row(raw) for raw in feature["AlignedPeakProperties"]]
        rows_by_spot[feature["MasterAlignmentID"]] = [
            r for r in rows if r and r.get("file_name") not in excluded_samples]
    names = [n for n in sample_factors.arf_sample_names(arf_state.features) if n not in excluded_samples]
    facets = sample_factors.build_sample_facets(names, arf_state.class_index)
    resolved_groups, caveats = gi.resolve_groups(groups, facets)
    plotted = {n for g in resolved_groups for n in g["samples"]}
    low = gi.resolve_sample_specs(low_reliability_samples or [], facets)

    catalog = arf2_reader.load_catalog(str(arf2_path))
    curation = {}
    if apply_curation:
        try:
            flag_set = curation_apply.flags_for_arf2(arf2_path)
        except curation_flags.FlagFileError as exc:
            return json_payload({"status": "error", **exc.details()})
        curation = {"wrong": flag_set["wrong"], "redundant": flag_set["redundant"],
                    "assign": flag_set["assign"]}
        if flag_set["orphaned"]:
            caveats.append(curation_flags.orphaned_warning(flag_set["orphaned"]))
    if exclude_auto_likely_wrong:
        review = curation_suggest.latest_review(
            arf2_path, curation_flags.alignment_key(arf2_path)["alignment_sha256"])
        if review is None:
            caveats.append("このアラインメントの curation_review のレビューが無いため、自動判定 likely_wrong は除いていません。")
        else:
            curation["auto_likely_wrong"] = {s["spot_id"] for s in review["spots"] if s.get("verdict") == "likely_wrong"}

    standard = set()
    if standard_samples:
        standard = gi.resolve_sample_specs(standard_samples, sample_factors.build_sample_facets(
            sample_factors.arf_sample_names(arf_state.features), arf_state.class_index))
    standard_only = gi.standard_only_spots(rows_by_spot, standard, plotted) if standard else frozenset()
    resolved_items, excluded = gi.resolve_items(items, catalog, curation=curation, standard_only=standard_only)

    # MS/MS の裏付け: 照合結果に MS/MS があり、かつ .arf2 の Name が `no MS2:` / `w/o MS2:` でない
    # （接頭辞は照合結果ではなく .arf2 の Name に付く）。assign 済みでも元のスペクトルの有無で判断する。
    annotations = load_spot_annotations(str(arf2_path))
    name_by_id = {int(r["MasterAlignmentID"]): r.get("Name") or "" for r in catalog}
    msms = {}
    for sid, ann in annotations.items():
        rep = ann.get("representative") or {}
        msms[sid] = (bool(rep.get("has_msms"))
                     and name_prefix(name_by_id.get(sid, "")) not in ("no MS2", "w/o MS2"))

    source = None
    if detection_limit is not None:
        detection_limit, source = float(detection_limit), "argument"
    else:
        param = analysis_params.find_param_file(arf2_path)
        value = analysis_params.read_analysis_params(param)["min_peak_height"] if param else None
        if value:
            detection_limit, source = value, "param_file"

    payload = gi.build_group_intensity_payload(
        resolved_items, resolved_groups, rows_by_spot, msms=msms, low_reliability=frozenset(low),
        excluded=excluded, detection_limit=detection_limit, detection_limit_source=source, caveats=caveats)
    arf_state.last_group_intensity = {"payload": payload, "title": title, "ncols": ncols}
    if mode == plot_render.PAYLOAD:
        return json_payload(payload)
    png = plot_render.figure_to_png(gi.render_group_intensity_plot(payload, title=title, ncols=ncols))
    return [_group_intensity_caption(payload), Image(data=png, format="png")]


def _group_intensity_caption(payload: dict) -> str:
    """画像に添える 1 行。図から読めない内訳（試料数・除外数・N.D.）を言葉で残す。"""
    groups = "、".join(f"{g['label']} n={len(g['samples'])}" for g in payload["groups"])
    nd = [it["item"] for it in payload["items"] if not it["detected"]]
    ex = payload["excluded"]
    caption = (f"群別強度: {len(payload['items'])} 項目（{groups}）。除外 — 内部標準 {len(ex['internal_standard'])}・"
               f"判断 {len(ex['curation'])}・自動判定 {len(ex['auto_likely_wrong'])}・標準液 {len(ex['standard_only'])}。")
    if nd:
        caption += f" N.D.: {'、'.join(nd)}。"
    if payload["caveats"]:
        caption += " 注意: " + " / ".join(payload["caveats"])
    return caption
```

（`arf2_reader`・`sample_factors`・`json_payload`・`Image`・`plot_render`・`mcp_errors` が `arf/tools.py` の既存 import に無ければ、モジュール先頭の import に足す。`name_prefix` は `metabolomix/arf2/match_results.py` の既存関数。）

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity_tools.py tests/test_group_intensity.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/core/session_state.py metabolomix/arf/tools.py tests/test_group_intensity_tools.py
git commit -m "feat(arf): arf_plot_group_intensity を足す"
```

---

### Task 6: 保存ツール `save_group_intensity_figure`

**Files:**
- Modify: `metabolomix/tools/reports.py`（`__all__` に追加）
- Test: `tests/test_group_intensity_tools.py`

**Interfaces:**
- Consumes: `session.arf.last_group_intensity`（Task 5）、`group_intensity.render_group_intensity_plot`。
- Produces: `save_group_intensity_figure(analysis_id: str, title: str | None = None) -> str`。
  `reports/figures/<slug>_group_intensity.png`（dpi 300）と `.svg` を書く。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_save_without_plot_is_missing_state():
    from metabolomix.tools.reports import save_group_intensity_figure
    session_state.session = session_state.AnalysisSession()
    out = json.loads(save_group_intensity_figure("x"))
    assert out["error"]["required_tools"] == ["arf_plot_group_intensity"]


def test_save_writes_png_and_svg(loaded, tmp_path, monkeypatch):
    from metabolomix.arf.tools import arf_plot_group_intensity
    from metabolomix.tools.reports import save_group_intensity_figure
    monkeypatch.setenv("LIPIDMIX_REPORTS_DIR", str(tmp_path / "reports"))
    arf_plot_group_intensity(items=["PG"], groups=["ctrl", "ko"], output="payload")
    msg = save_group_intensity_figure("EV membrane")
    figures = tmp_path / "reports" / "figures"
    pngs = list(figures.glob("*_group_intensity.png"))
    assert len(pngs) == 1 and pngs[0].with_suffix(".svg").is_file()
    assert "group_intensity.png" in msg
```

- [ ] **Step 2: 失敗を確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity_tools.py -q -k save`
Expected: FAIL（`ImportError`）

- [ ] **Step 3: 実装**（`reports.py`、`save_eic_figure` の後ろ。`__all__` に `"save_group_intensity_figure"`）

```python
@mcp.tool(annotations=ToolAnnotations(
    readOnlyHint=False, destructiveHint=False, idempotentHint=True), structured_output=False)
def save_group_intensity_figure(analysis_id: str, title: str | None = None) -> str:
    """直前の arf_plot_group_intensity の図を reports/figures/<analysis_id>_group_intensity.png（dpi 300）と
    同名 .svg に保存し、相対パスを返す。レポート本文に `![...](figures/<...>.png)` で埋め込める。"""
    from metabolomix.plots.group_intensity import render_group_intensity_plot

    last = getattr(session_state.session.arf, "last_group_intensity", None)
    if not last or not last.get("payload"):
        return mcp_errors.missing_state(
            "group_intensity_plot", ["arf_plot_group_intensity"],
            "先に arf_plot_group_intensity を実行してください（群別強度の図がありません）。")
    slug = knowledge_store.make_slug(analysis_id)
    figures_dir = _resolve_report_dir() / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_path = figures_dir / f"{slug}_group_intensity.png"
    fig = render_group_intensity_plot(last["payload"], title=title or last.get("title"), ncols=last.get("ncols"))
    try:
        fig.savefig(out_path, dpi=300, format="png", bbox_inches="tight")
        fig.savefig(out_path.with_suffix(".svg"), format="svg", bbox_inches="tight")
    finally:
        plt.close(fig)
    rel = f"figures/{out_path.name}"
    return (f"群別強度の図を保存: {out_path}（同名の .svg も保存）\n"
            f"本文に ![group intensity]({rel}) で埋め込めます。")
```

- [ ] **Step 4: 通ることを確かめる**

Run: `C:/Python314/python.exe -m pytest tests/test_group_intensity_tools.py -q`
Expected: PASS

- [ ] **Step 5: コミット**（背景実行）

```bash
git add metabolomix/tools/reports.py tests/test_group_intensity_tools.py
git commit -m "feat(reports): save_group_intensity_figure を足す"
```

---

### Task 7: 登録・文書・試作の片付け

**Files:**
- Modify: `tests/test_server_registration.py`（`EXPECTED_TOOLS` に 2 件、変更履歴コメント「(72→74)」）
- Modify: `tests/test_tool_annotations.py`（`EXPECTED_ANNOTATIONS`: `arf_plot_group_intensity` は `READ_ONLY`、`save_group_intensity_figure` は `LOCAL_WRITE`）
- Modify: `tests/test_missing_state_contract.py`（空セッションで 2 ツールが `missing_state` を返すケース。docstring の件数も更新）
- Modify: `tests/test_workflow_docs.py`（`IN_SCOPE` の `"arf.md"` に `arf_plot_group_intensity`、`"plots.md"` に `save_group_intensity_figure`。合計の定数 51 → 53）
- Modify: `USAGE.md`（見出し「全72ツール」→「全74ツール」。「## 9. 図の保存・レポート」の表に 2 行）
- Modify: `CLAUDE.md`（冒頭の規模「ツール 72」→「ツール 74」）
- Modify: `docs/workflow/index.md`（arf.md 10 → 11、plots.md 3 → 4、合計 51 → 53、登録総数 72 → 74（53 + 21））
- Modify: `docs/workflow/arf.md`（`## arf_plot_group_intensity` の節。mermaid に `ARF -->|session.arf.last_group_intensity| GI[arf_plot_group_intensity]`）
- Modify: `docs/workflow/plots.md`（`## save_group_intensity_figure` の節。payload 契約の比較表に `lipidmix.group_intensity.v1` の行）
- Modify: `docs/output_format/arf.md`（`### 11.5 arf_plot_group_intensity の返り値`: spec §4・§5・§7 のフィールドの意味、N.D. の 2 種類、MS1-only、gap-fill 過半の ◆、除外の 4 種類）
- Modify: `check.py`（試作を消し、冒頭 3 行のコメントだけに戻す）

- [ ] **Step 1: 検査テストの定数を先に変えて落ちることを確かめる**

`test_server_registration.py` の `EXPECTED_TOOLS` に 2 件、`test_tool_annotations.py` に 2 件、`test_workflow_docs.py` の `IN_SCOPE` と合計 53 を入れる。

Run: `C:/Python314/python.exe -m pytest tests/test_server_registration.py tests/test_tool_annotations.py tests/test_readme_links.py tests/test_workflow_docs.py -q`
Expected: FAIL（USAGE.md の行が無い、件数 72 の宣言、workflow の節が無い）

- [ ] **Step 2: USAGE.md の 2 行**（「## 9. 図の保存・レポート」の表）

```markdown
| `arf_plot_group_intensity` | 選んだクラス・分子種ごとに、試料群の試料別強度を並べる(`items`, `groups`, `low_reliability_samples`, `apply_curation`, `exclude_auto_likely_wrong`, `standard_samples`, `detection_limit`, `title`, `ncols`, `output`)。1 項目 1 パネル(1〜30)。`items` は Ontology に完全一致ならクラス合計、しなければ分子種名(`\|` で分けた候補名に完全一致)、`+` で合算、当たらなければ N.D. パネル。`groups` は `arf_differential` と同じトークン規則で 1 群以上(QC・ブランクは `qc`/`blank` を書いたときだけ)。縦軸は log10(PeakHeight 合計) で全パネル共通、群は log10 の平均 ± SD、検定はしない。クラス合計から標識内部標準と `standard_samples` の試料にだけある分子種を除き、`apply_curation` で人の判断、`exclude_auto_likely_wrong` で最新レビューの likely_wrong を除く。既定は PNG 画像、`output="payload"` で `lipidmix.group_intensity.v1`。前提: ARF の読み込みと兄弟 `.arf2`(無ければ `missing_state`)。 |
| `save_group_intensity_figure` | 直前の `arf_plot_group_intensity` の図を `reports/figures/<analysis_id>_group_intensity.png`(dpi 300)と同名 `.svg` に保存(`analysis_id`, `title`)。図が無ければ `missing_state`。 |
```

- [ ] **Step 3: 件数の宣言**: `USAGE.md` 1 行目「全72ツール」→「全74ツール」、`CLAUDE.md` の「ツール 72・リソース 5・リソーステンプレート 3」→「ツール 74・…」、`docs/workflow/index.md` の表と合計。

- [ ] **Step 4: workflow 文書**

`docs/workflow/arf.md`（`## arf_plot_volcano` の後ろ）:

```markdown
## arf_plot_group_intensity

前提: ARF の読み込み（`load_dataset` / `arf_parser`）と兄弟 `.arf2`。
状態変更: `session.arf.last_group_intensity`（`save_group_intensity_figure` が読む）。

選んだクラス・分子種ごとに試料群の試料別強度を並べる。項目と群の解決、payload、描画は
`metabolomix/plots/group_intensity.py` の純関数。キュレーションの判断と最新レビューを読む。

1. metabolomix/arf/tools.py  arf_plot_group_intensity()
2. └─ metabolomix/plots/render.py  resolve_plot_output()
3. └─ metabolomix/core/mcp_errors.py  missing_state()
4. └─ metabolomix/arf/tools.py  _sibling_arf2_path()
5. └─ metabolomix/arf/reader.py  alignment_feature_row()
6. └─ metabolomix/msdial/sample_factors.py  arf_sample_names()
7. └─ metabolomix/msdial/sample_factors.py  build_sample_facets()
8. └─ metabolomix/plots/group_intensity.py  resolve_groups()
9. │  └─ metabolomix/msdial/sample_factors.py  expand_sample_specs()
10.└─ metabolomix/plots/group_intensity.py  resolve_sample_specs()
11.└─ metabolomix/arf2/reader.py  load_catalog()
12.└─ [apply_curation] metabolomix/curation/apply.py  flags_for_arf2()
13.└─ [exclude_auto_likely_wrong] metabolomix/curation/suggest.py  latest_review()
14.└─ [standard_samples] metabolomix/plots/group_intensity.py  standard_only_spots()
15.└─ metabolomix/plots/group_intensity.py  resolve_items()
16.└─ metabolomix/arf2/match_results.py  load_spot_annotations()
17.└─ metabolomix/msdial/analysis_params.py  find_param_file()
18.└─ metabolomix/msdial/analysis_params.py  read_analysis_params()
19.└─ metabolomix/plots/group_intensity.py  build_group_intensity_payload()
20.└─ [output=image] metabolomix/plots/group_intensity.py  render_group_intensity_plot()
21.   └─ metabolomix/plots/render.py  figure_to_png()
```

（2 桁の番号の行は `CHAIN_RE` に合わせ、既存の 2 桁の行と同じ書き方にそろえる。）

`docs/workflow/plots.md`（`## save_eic_figure` の後ろ）:

```markdown
## save_group_intensity_figure

前提: `session.arf.last_group_intensity`（`arf_plot_group_intensity`）。
状態変更: なし（ファイルを書く）。

1. metabolomix/tools/reports.py  save_group_intensity_figure()
2. └─ metabolomix/core/mcp_errors.py  missing_state()
3. └─ metabolomix/corpus/knowledge_store.py  make_slug()
4. └─ metabolomix/core/mcp_core.py  _resolve_report_dir()
5. └─ metabolomix/plots/group_intensity.py  render_group_intensity_plot()
```

- [ ] **Step 5: 出力形式の文書** `docs/output_format/arf.md` に `### 11.5 arf_plot_group_intensity の返り値` を足す（spec §7 のフィールド表、`detected` / N.D. の 2 種類、`n_spots_msms` = 0 の MS1-only、`gap_filled_fraction` > 0.5 の ◆、`excluded` の 4 種類、`detection_limit_source` の `argument` / `param_file` / null、検定をしないこと）。

- [ ] **Step 6: `test_missing_state_contract.py` に 2 ケース**

```python
    def test_group_intensity_without_arf(self):
        from metabolomix.arf.tools import arf_plot_group_intensity
        out = json.loads(arf_plot_group_intensity(items=["PG"], groups=["x"], output="payload"))
        self.assertEqual(out["error"]["code"], "missing_state")

    def test_save_group_intensity_without_plot(self):
        from metabolomix.tools.reports import save_group_intensity_figure
        out = json.loads(save_group_intensity_figure("x"))
        self.assertEqual(out["error"]["required_tools"], ["arf_plot_group_intensity"])
```

（既存ケースと同じクラスに足し、空セッションの準備は既存の setUp に従う。）

- [ ] **Step 7: `check.py` を片付ける**（冒頭 3 行のコメントだけ残す）

```python
# スクラッチ用ファイル（pseudo-workspace）。
# 一時的な検証・確認コードをここに書き、うまくいったら適切なモジュールへ移し、
# このファイルの中身は消す。詳細は README / docs/HISTRY.md を参照。
```

- [ ] **Step 8: 全テスト**

Run: `C:/Python314/python.exe -m pytest tests -q`
Expected: PASS（件数は増える）

- [ ] **Step 9: コミット**（背景実行）

```bash
git add -A
git commit -m "docs: 群別強度プロットの文書・登録テストを足し、check.py の試作を片付ける"
```

- [ ] **Step 10: 実データで確かめる**（研究フォルダには書かない。`LIPIDMIX_REPORTS_DIR` を scratchpad に向ける）

20260930_EV の neg/pos で `items=["PC","PE","PS","PI","PG","PA","CL","LPC","LPE","LPG","SM","Cer_NS+Cer_AP+Cer_NDS+Cer_HS","MGDG","DGDG"]`、
`groups=["blank","EV_cell","EV_PlnA","cell_C","cell_PlnA"]`、`low_reliability_samples=["cell_PlnA_1"]`、`standard_samples=["systemlot8"]`、
`exclude_auto_likely_wrong=True` を呼び、試作の図（`figures/membrane_class_intensity_{neg,pos}.png`）と同じ読み取りになることを目で確かめる。

- [ ] **Step 11: vault の流れ図**: ARF 経路に任意の描画ツールとして `arf_plot_group_intensity` → `save_group_intensity_figure` を足し、出典行の日付を更新する（リポジトリ外・追跡外）。
